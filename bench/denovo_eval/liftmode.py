"""LIFTMODE prototype (PENDING_2026-10-27 item: Liftoff-style de novo family mode; template rule from ASMTPL-3, docs/LEDGER_O1.md).

  liftmode.py templates GTF FASTA OUTDIR      every unique intron chain of the assembled GTF (mono-exonic: unique span) -> templates.fa + meta
  liftmode.py families OUTDIR PAF...          Liftoff chaining of each template (coverage >= 0.5, identity >= 0.5), hit judged on the target
                                              LOCUS's exon bases (>= 0.5 of the chained bases on exons of one gene_id); family = connected
                                              component of loci joined by a template placement; intronic/intergenic landings = putative copies.
                                              Writes <OUTDIR>/copies.tsv (famtruth catalog format) and putative.bed.
No MCL, no thresholds beyond Liftoff's -a/-s defaults.
"""
import bisect
import collections
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'guided_eval'))
import asmtpl as A  # noqa: E402
import cdsseed as C  # noqa: E402
import guided3 as g3  # noqa: E402

COV, IDENT = 0.5, 0.5
RECIPROCAL = os.environ.get('LIFT_RECIPROCAL') == '1'  # LIFTMODE-GW2 arm R
RETRO = os.environ.get('LIFT_RETRO') == '1'  # LIFTMODE-GW3
MASK = os.environ.get('LIFT_MASK', '')  # LIFTMODE-GW4: 'rep' (one chain per locus) or 'shared' (locus-shared exons)


def templates(gtf, fasta, outdir):
    tx = A.read_gtf(gtf)
    gene_of_tx = {}
    for line in open(gtf):
        f = line.split('\t', 8)
        if len(f) > 8 and f[2] == 'transcript':
            a = dict(__import__('re').findall(r'(\S+) "([^"]*)"', f[8]))
            gene_of_tx[a['transcript_id']] = a['gene_id']
    chains = {}
    for tid, t in tx.items():
        ex = t['exons']
        key = (t['chrom'], t['strand'], tuple((a, b) for (_, a), (b, _) in zip(ex, ex[1:])) or tuple(ex))
        cur = chains.get(key)
        ln = sum(b - a for a, b in ex)
        if cur is None or (t['reads'], ln) > (cur['reads'], cur['len']):
            chains[key] = dict(locus=gene_of_tx[tid], tid=tid, chrom=t['chrom'], strand=t['strand'], exons=ex, reads=t['reads'], len=ln)
    C.FA = fasta
    regs = sorted({f"{t['chrom']}:{a + 1}-{b}" for t in chains.values() for a, b in t['exons']})
    seqs = C.faidx(regs, f'{outdir}/ex')
    meta = {}
    with open(f'{outdir}/templates.fa', 'w') as w:
        for i, t in enumerate(chains.values()):
            s = ''.join(seqs[f"{t['chrom']}:{a + 1}-{b}"] for a, b in t['exons'])
            if t['strand'] == '-':
                s = s.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]
            name = f'T{i}'
            w.write(f'>{name}\n{s}\n')
            meta[name] = dict(locus=t['locus'], tid=t['tid'], chrom=t['chrom'], strand=t['strand'], exons=t['exons'], reads=t['reads'], len=len(s))
    loci = collections.defaultdict(lambda: dict(chrom=None, strand=None, exons=[], reads=0, tid=None, best=0))
    for t in chains.values():
        L = loci[t['locus']]
        L['chrom'], L['strand'] = t['chrom'], t['strand']
        L['exons'] += t['exons']
        L['reads'] += t['reads']
        if t['len'] > L['best']:
            L['best'], L['tid'], L['rep'] = t['len'], t['tid'], t['exons']
    json.dump(dict(templates=meta, loci=loci), open(f'{outdir}/meta.json', 'w'))
    print(f'transcripts {len(tx)}, unique chains {len(chains)}, loci {len(loci)}, template bp {sum(m["len"] for m in meta.values())}')


def locus_index(loci):
    idx = collections.defaultdict(list)
    for g, L in loci.items():
        for a, b in L['exons']:
            idx[(L['chrom'], L['strand'])].append((a, b, g))
    out = {}
    for k, v in idx.items():
        v.sort()
        out[k] = (v, [x[0] for x in v], max(b - a for a, b, _ in v))
    return out


def locus_of(idx, chrom, strand, bl):
    """Locus whose exons share most bases with the blocks, plus the exon fraction of the blocks."""
    if (chrom, strand) not in idx:
        return None, 0.0
    v, st, m = idx[(chrom, strand)]
    share = collections.Counter()
    for s, e in bl:
        i = bisect.bisect_left(st, s - m)
        while i < len(v) and v[i][0] < e:
            x, y, g = v[i]
            if y > s:
                share[g] += min(y, e) - max(x, s)
            i += 1
    if not share:
        return None, 0.0
    g = max(share, key=lambda k: (share[k], k))
    return g, share[g] / max(1, sum(b - a for a, b in bl))


def masks(tm):
    """GW4: template -> kept intervals in template coordinates (None = dropped); 'rep' keeps whole templates of the most-reads chain."""
    by_locus = collections.defaultdict(list)
    for n, t in tm.items():
        by_locus[t['locus']].append(n)
    out = {}
    if MASK == 'rep':
        for L, ns in by_locus.items():
            best = max(ns, key=lambda n: (tm[n]['reads'], tm[n]['len']))
            out[best] = [(0, tm[best]['len'])]
        return out
    for L, ns in by_locus.items():
        if len(ns) == 1:
            n = ns[0]
            out[n] = [(0, tm[n]['len'])]
            continue
        exact, inner = collections.Counter(), collections.Counter()
        for n in ns:
            ex = tm[n]['exons']
            for e in ex:
                exact[tuple(e)] += 1
            if len(ex) > 1:
                inner[('first', ex[0][1])] += 1
                inner[('last', ex[-1][0])] += 1
        for n in ns:
            ex = tm[n]['exons']
            keep = [i for i, e in enumerate(ex) if len(ex) == 1 or exact[tuple(e)] >= 2
                    or (i == 0 and inner[('first', e[1])] >= 2) or (i == len(ex) - 1 and inner[('last', e[0])] >= 2)]
            if not keep:
                continue
            # template coordinates run 5'->3' in transcript orientation
            order = ex if tm[n]['strand'] == '+' else ex[::-1]
            pos, iv = 0, []
            for e in order:
                ln = e[1] - e[0]
                if (ex.index(e) in keep):
                    iv.append((pos, pos + ln))
                pos += ln
            out[n] = iv
    return out


def families(outdir, pafs):
    M = json.load(open(f'{outdir}/meta.json'))
    tm, loci = M['templates'], M['loci']
    mk = masks(tm) if MASK else None
    g3.MASK_ALL = mk is not None
    if mk is not None:
        tm = {n: dict(t, len=sum(b - a for a, b in mk[n])) for n, t in tm.items() if n in mk}
        print(f'mask {MASK}: templates kept {len(tm)}, bp {sum(t["len"] for t in tm.values())}')
    idx = locus_index(loci)
    groups = collections.defaultdict(list)
    land = collections.defaultdict(list)
    nblk = {}
    spliced_locus = {g: len(L['exons']) > 1 for g, L in loci.items()}
    n_retro = 0
    for q, t, st, bl, qi, ti, cdsb, ident in g3.records(pafs, mk or {}):
        if q not in tm or not bl or cdsb <= 0:
            continue
        g, frac = locus_of(idx, t, st, bl)
        if g == tm[q]['locus']:
            continue
        key = (q, g if g and frac >= 0.5 else f'{t}:{bl[0][0]}-{bl[-1][1]}', st)
        groups[key].append((qi, ti, cdsb, ident))
        nblk[key] = max(nblk.get(key, 0), len(bl))
        land[key].append((t, bl[0][0], bl[-1][1]))
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    putative = []
    retro = []
    n_place = 0
    placed = set()  # (source locus, target locus)
    edges = {}
    for (q, g, st), rs in groups.items():
        bases, ident = g3.chain(rs, st)
        if bases / tm[q]['len'] < COV or ident < IDENT:
            continue
        n_place += 1
        if g in loci:
            a, b = tm[q]['locus'], g
            if RETRO and len(tm[q]['exons']) > 1 and nblk[(q, g, st)] == 1 and spliced_locus[b]:
                n_retro += 1  # GW3: a spliced template landing intronless inside a spliced locus = processed copy, not a family edge
                retro.append((a, b))
                continue
            placed.add((a, b))
            w = ident * min(1.0, bases / tm[q]['len'])
            edges[(min(a, b), max(a, b))] = max(edges.get((min(a, b), max(a, b)), 0.0), w)
            if not RECIPROCAL or (b, a) in placed:
                parent[find(a)] = find(b)
        else:
            c, s, e = min(land[(q, g, st)]), min(x[1] for x in land[(q, g, st)]), max(x[2] for x in land[(q, g, st)])
            putative.append((land[(q, g, st)][0][0], s, e, tm[q]['locus'], round(ident, 3)))
    with open(f'{outdir}/edges.tsv', 'w') as w:  # locus pairs with a placement, for --extra-edges (arm X)
        for (a, b), wt in sorted(edges.items()):
            if not RECIPROCAL or ((a, b) in placed and (b, a) in placed):
                w.write(f'{a}\t{b}\t{wt:.6f}\n')
    with open(f'{outdir}/retro.tsv', 'w') as w:
        for a, b in sorted(set(retro)):
            w.write(f'{a}\t{b}\n')
    comp = collections.defaultdict(set)
    for g in loci:
        comp[find(g)].add(g)
    fams = [sorted(v) for v in comp.values() if len(v) >= 2]
    fams.sort(key=lambda v: (-len(v), v[0]))
    with open(f'{outdir}/copies.tsv', 'w') as w:
        w.write('family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\tmax_family_identity\n')
        for i, v in enumerate(fams):
            for j, g in enumerate(v):
                L = loci[g]
                ex = L['rep']
                w.write(f"LIFT{i}\t{j}\t{L['tid']}\t{L['chrom']}\t{ex[0][0]}\t{ex[-1][1]}\t{len(ex)}\t{L['strand']}\t{L['reads']}\t"
                        f"{','.join(f'{a}-{b}' for a, b in ex)}\t0\n")
    merged = []
    for c, s, e, src, ident in sorted(putative):
        if merged and merged[-1][0] == c and s <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], e)
            merged[-1][3].add(src)
        else:
            merged.append([c, s, e, {src}])
    with open(f'{outdir}/putative.bed', 'w') as w:
        for c, s, e, srcs in merged:
            w.write(f'{c}\t{s}\t{e}\t{";".join(sorted(srcs))}\n')
    sz = [len(v) for v in fams]
    print(f'templates {len(tm)}, placements {n_place}, families {len(fams)} (loci in families {sum(sz)}, largest {sz[:5]}), '
          f'putative copies {len(merged)} (from {len(putative)} landings), retro placements {n_retro}')


if __name__ == '__main__':
    if sys.argv[1] == 'templates':
        os.makedirs(sys.argv[4], exist_ok=True)
        templates(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        families(sys.argv[2], sys.argv[3:])
