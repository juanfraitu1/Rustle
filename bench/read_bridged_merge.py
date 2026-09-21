#!/usr/bin/env python3
"""Read-bridged node merge (RBM), per `docs/PREREG_read_bridged_node_merge_2026-09-21.md`.

The de novo split-locus pathology is not a distance problem: 53.7% of adjacent same-gene locus pairs
OVERLAP in span (median gap -1,575 bp), and 94.1% of the intra-gene spliced pairs overlap only
PARTIALLY in exon -- the break falls mid-exon. Junction-disjointness between loci is BY CONSTRUCTION
(0 of 13,685 junctions are used by two loci), so a junction rule cannot see the split at all.

    H(r) = loci where primary read r (-F 2308, MAPQ >= Q) covers >= 25 bp of exonic sequence UNIQUE to
           that locus relative to the partner
    merge A,B  iff  same strand AND |{r : {A,B} subset of H(r)}| >= N

No distance term. Register 331 refuted the SUPPLEMENTARY-alignment form of this ("every one at a single
read"); these are primary alignments with a count requirement.

Scoring is one-to-one with collisions as misses, over a universe FIXED on the unmerged output -- see the
prereg for why each of those is load-bearing.
"""
import argparse
import bisect
import collections
import re
import sys

MIN_UNIQ_BP = 25


def load_gtf(path, chrom):
    tx = collections.defaultdict(lambda: collections.defaultdict(list))
    strand = {}
    for line in open(path):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] != 'exon':
            continue
        g = re.search(r'gene_id "([^"]+)"', f[8])
        t = re.search(r'transcript_id "([^"]+)"', f[8])
        if not g or not t:
            continue
        tx[g.group(1)][t.group(1)].append((int(f[3]), int(f[4])))
        strand[g.group(1)] = f[6]
    return tx, strand


def merge_iv(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def ov_bp(A, B):
    """overlap bp between two sorted disjoint interval lists (1-based inclusive)."""
    i = j = t = 0
    while i < len(A) and j < len(B):
        s = max(A[i][0], B[j][0]); e = min(A[i][1], B[j][1])
        if e >= s:
            t += e - s + 1
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return t


def subtract(A, B):
    """A minus B, both sorted disjoint."""
    out = []
    for s, e in A:
        cur = [(s, e)]
        for bs, be in B:
            if be < s or bs > e:
                continue
            nxt = []
            for cs, ce in cur:
                if be < cs or bs > ce:
                    nxt.append((cs, ce)); continue
                if cs < bs:
                    nxt.append((cs, bs - 1))
                if ce > be:
                    nxt.append((be + 1, ce))
            cur = nxt
        out += cur
    return [(s, e) for s, e in out if e >= s]


def bridge_votes(bam_path, chrom, exm, spans, tiers):
    """ONE pass over the BAM; returns {mapq_floor: Counter over unordered locus pairs}.

    Three passes would triple the only expensive step, so every MAPQ tier is accumulated at once."""
    import pysam
    bam = pysam.AlignmentFile(bam_path, 'rb')
    order = sorted(spans.items(), key=lambda kv: kv[1][0])
    ids = [k for k, _ in order]
    starts = [v[0] for _, v in order]
    ends = [v[1] for _, v in order]
    maxend = []
    m = 0
    for e in ends:
        m = max(m, e); maxend.append(m)
    votes = {t: collections.Counter() for t in tiers}
    uniq_cache = {}
    n_read = 0
    for r in bam.fetch(chrom):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        mq = r.mapping_quality
        if mq < tiers[0]:
            continue
        n_read += 1
        rs, re_ = r.reference_start + 1, r.reference_end
        # candidate loci whose SPAN overlaps the read
        hi = bisect.bisect_right(starts, re_)
        cand = []
        for j in range(hi - 1, -1, -1):
            if maxend[j] < rs:
                break
            if ends[j] >= rs:
                cand.append(ids[j])
        if len(cand) < 2:
            continue
        blocks = [(b + 1, e) for b, e in r.get_blocks()]
        blocks = merge_iv(blocks)
        hit = [g for g in cand if ov_bp(blocks, exm[g]) >= MIN_UNIQ_BP]
        if len(hit) < 2:
            continue
        hit.sort()
        for a_i in range(len(hit)):
            for b_i in range(a_i + 1, len(hit)):
                a, b = hit[a_i], hit[b_i]
                key = (a, b)
                if key not in uniq_cache:
                    uniq_cache[key] = (subtract(exm[a], exm[b]), subtract(exm[b], exm[a]))
                ua, ub = uniq_cache[key]
                if not ua or not ub:
                    continue
                if ov_bp(blocks, ua) >= MIN_UNIQ_BP and ov_bp(blocks, ub) >= MIN_UNIQ_BP:
                    for t in tiers:
                        if mq >= t:
                            votes[t][key] += 1
    return votes, n_read


def apply_merge(votes, strand, n_min, exm=None, require_overlap=False):
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    applied = 0
    for (a, b), v in sorted(votes.items()):
        if v < n_min or strand[a] != strand[b]:
            continue
        if require_overlap and ov_bp(exm[a], exm[b]) <= 0:
            continue
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb; applied += 1
    groups = collections.defaultdict(list)
    for g in strand:
        groups[find(g)].append(g)
    return groups, applied


def load_genes(gff, chrom):
    g = {}
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] not in ('gene', 'pseudogene'):
            continue
        n = re.search(r'Name=([^;]+)', f[8])
        if n:
            g[n.group(1)] = (int(f[3]), int(f[4]))
    return g


def score(groups, exm, genes, universe):
    """one-to-one, collisions are misses, universe FIXED by the caller."""
    node_ex = {}
    for root, mem in groups.items():
        node_ex[root] = merge_iv([iv for g in mem for iv in exm[g]])
    order = sorted(node_ex.items(), key=lambda kv: kv[1][0][0])
    nids = [k for k, _ in order]
    nstart = [v[0][0] for _, v in order]
    nend = [v[-1][1] for _, v in order]
    maxend = []
    m = 0
    for e in nend:
        m = max(m, e); maxend.append(m)

    def overlapping(s, e):
        hi = bisect.bisect_right(nstart, e)
        out = []
        for j in range(hi - 1, -1, -1):
            if maxend[j] < s:
                break
            if nend[j] >= s:
                out.append(nids[j])
        return out

    claim = {}
    for gname in universe:
        gs, ge = genes[gname]
        best, bo = None, 0
        # ties broken by node id, deterministically: dict iteration order otherwise moved 5 genes
        # between runs and made the headline irreproducible.
        for nid in sorted(overlapping(gs, ge)):
            o = ov_bp(node_ex[nid], [(gs, ge)])
            if o > bo:
                bo, best = o, nid
        claim[gname] = best
    owner = collections.defaultdict(list)
    for gname, nid in claim.items():
        if nid:
            owner[nid].append(gname)
    correct = 0
    for gname, nid in claim.items():
        if not nid or len(owner[nid]) > 1:
            continue
        gs, ge = genes[gname]
        tot = sum(e - s + 1 for s, e in node_ex[nid])
        if tot and ov_bp(node_ex[nid], [(gs, ge)]) / tot >= 0.50:
            correct += 1
    # false merge / junk over ALL nodes
    fm = junk = 0
    for nid, iv in node_ex.items():
        hits = 0
        for gname, (gs, ge) in genes.items():
            if ge < iv[0][0] or gs > iv[-1][1]:
                continue
            if ov_bp(iv, [(gs, ge)]) >= 100:
                hits += 1
        if hits >= 2:
            fm += 1
        elif hits == 0:
            junk += 1
    return dict(correct=correct, n=len(universe), nodes=len(node_ex),
                false_merge=fm, junk=junk, collisions=sum(1 for v in owner.values() if len(v) > 1))


def main():
    ap = argparse.ArgumentParser()
    for x in ('--gtf', '--bam', '--gff', '--chrom'):
        ap.add_argument(x, required=True)
    ap.add_argument('--mapq-tiers', default='0,1,60', help='MAPQ floors to report, one table each')
    ap.add_argument('--sweep', default='1,3,5,10')
    ap.add_argument('--emit', help='write the merged GTF at --emit for the given --n')
    ap.add_argument('--n', type=int)
    ap.add_argument('--emit-q', type=int, default=60)
    ap.add_argument('--require-exon-overlap', action='store_true',
                    help='POST-HOC (not in the prereg): merge only pairs whose exons overlap. The '
                         'pre-registered unrestricted rule was refuted because 50.7%% of bridged pairs '
                         'join DIFFERENT genes at a 30 kb median distance -- readthrough transcription.')
    a = ap.parse_args()

    tx, strand = load_gtf(a.gtf, a.chrom)
    exm = {g: merge_iv([iv for e in t.values() for iv in e]) for g, t in tx.items()}
    spans = {g: (v[0][0], v[-1][1]) for g, v in exm.items()}
    genes = load_genes(a.gff, a.chrom)
    print(f'{a.chrom}: {len(exm)} de novo loci | {len(genes)} annotated genes', file=sys.stderr)

    base_groups = {g: [g] for g in exm}
    # UNIVERSE fixed on the UNMERGED output (register 770)
    universe = []
    order = sorted(spans.items(), key=lambda kv: kv[1][0])
    for gname, (gs, ge) in genes.items():
        for lid, (s, e) in spans.items():
            if not (e < gs or s > ge) and ov_bp(exm[lid], [(gs, ge)]) > 0:
                universe.append(gname); break
    universe = sorted(set(universe))
    print(f'  fixed universe: {len(universe)} genes', file=sys.stderr)

    tiers = [int(x) for x in a.mapq_tiers.split(',')]
    votes, nread = bridge_votes(a.bam, a.chrom, exm, spans, tiers)
    print(f'  {nread} primary reads | bridged pairs per MAPQ floor: '
          + ', '.join(f'Q>={t}: {len(votes[t])}' for t in tiers), file=sys.stderr)

    base = score(base_groups, exm, genes, universe)
    b = base['correct'] / base['n']
    print(f'\n{a.chrom}   universe {base["n"]} genes (FIXED on the unmerged output)')
    print(f'  {"arm":14s} {"merges":>7} {"nodes":>7} {"correct":>9} {"rate":>8} {"delta":>9} '
          f'{"collis":>7} {"falsemrg%":>10} {"junk":>6}')
    s0 = base
    print(f'  {"baseline":14s} {0:>7} {s0["nodes"]:>7} {s0["correct"]:>9} {b:>8.4f} '
          f'{0.0:>+8.2f}pp {s0["collisions"]:>7} {100*s0["false_merge"]/s0["nodes"]:>9.1f}% {s0["junk"]:>6}')
    for t in tiers:
        for n in [int(x) for x in a.sweep.split(',')]:
            groups, applied = apply_merge(votes[t], strand, n, exm, a.require_exon_overlap)
            s = score(groups, exm, genes, universe)
            r = s['correct'] / s['n']
            print(f'  {"Q>=%d N>=%d" % (t, n):14s} {applied:>7} {s["nodes"]:>7} {s["correct"]:>9} '
                  f'{r:>8.4f} {(r-b)*100:>+8.2f}pp {s["collisions"]:>7} '
                  f'{100*s["false_merge"]/s["nodes"]:>9.1f}% {s["junk"]:>6}')
            if a.emit and a.n == n and a.emit_q == t:
                write_gtf(a.emit, a.gtf, a.chrom, groups)


def write_gtf(out, src, chrom, groups):
    root = {}
    for r, mem in groups.items():
        for g in mem:
            root[g] = r
    with open(out, 'w') as fh:
        for line in open(src):
            if line.startswith('#'):
                fh.write(line); continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[0] != chrom:
                fh.write(line); continue
            m = re.search(r'gene_id "([^"]+)"', f[8])
            if m and m.group(1) in root:
                f[8] = f[8].replace(f'gene_id "{m.group(1)}"', f'gene_id "{root[m.group(1)]}"')
            fh.write('\t'.join(f) + '\n')


if __name__ == '__main__':
    main()
