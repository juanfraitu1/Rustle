#!/usr/bin/env python3
"""Advisor's proposal (2026-09-25): replace the exon-sum edge rule with (a) read-set JACCARD between loci and
(b) MULTIMAPPING READS AS CHAINING ANCHORS. Guided nodes (annotated genes, exon unions) so every truth gene is a
node in every arm; only the EDGE rule changes. Report: docs/archive/2026-09/ADVISOR_JACCARD_ANCHORS_2026-09-25.md

  S      shipped: asm20 all-vs-all of gene spans -> mcl_families --min-exonic-bp 1 --min-shared-exon-frac 0.60
  J(t)   Jaccard of molecule sets (primary OR secondary block on the exons), counting a shared molecule only when
         two DIFFERENT records reach g and h (multimapping, not readthrough); Jany(t) counts any sharing
  C      shared reads' alignments give anchors (pos in g, pos in h) at the same read base; LIS-chain them;
         edge iff the chain covers >= 0.60 of the SMALLER gene's exonic bases (the exon-sum rule, read-built)
  J+C    both
Every arm is clustered twice from its edge list: mcl_port (I=2.8) and connected components (.cc, the operator-free control)
and S is re-clustered from its dumped graph so the operator is held fixed. Exon-overlapping pairs are never read edges.

usage: jaccard_anchor_test.py TAG BAM GFF FASTA CHROM TRUTH OUTDIR
"""
import sys, os, re, bisect, subprocess, collections, itertools

tag, bam, gff, fasta, chrom, truth, out = sys.argv[1:8]
os.makedirs(out, exist_ok=True)
P = lambda s: os.path.join(out, f'{tag}.{s}')
R = '/mnt/linuxdisk/home/juanfraitu/rustle_target/release'
BIN = 50           # anchor / coverage bin (bp)
MIN_EX_OV = 20     # a block must put >= this many bases on a gene's exons to count as support
COV = 0.60         # same C as the shipped --min-shared-exon-frac


import time, pickle
T0 = time.time()


def log(*a):
    print(f'[{time.time() - T0:6.0f}s]', *a, file=sys.stderr, flush=True)


# ---------- 1. gene nodes: gene/pseudogene records, exon union (exonless -> span), as in mcl_families ----------
genes, rna2gene, exons = {}, {}, collections.defaultdict(list)
for ln in open(gff):
    if ln[0] == '#':
        continue
    r = ln.rstrip('\n').split('\t')
    if len(r) < 9 or r[0] != chrom:
        continue
    at = dict(kv.split('=', 1) for kv in r[8].split(';') if '=' in kv)
    if r[2] in ('gene', 'pseudogene'):
        genes[at['ID']] = dict(name=at.get('Name', at['ID']), s=int(r[3]), e=int(r[4]), strand=r[6])
    elif 'Parent' in at and r[2] not in ('exon', 'CDS'):
        if 'ID' in at:
            rna2gene[at['ID']] = at['Parent']
    if r[2] == 'exon' and 'Parent' in at:
        exons[at['Parent']].append((int(r[3]), int(r[4])))
gex = collections.defaultdict(list)
for par, xs in exons.items():
    g = par if par in genes else rna2gene.get(par)
    while g is not None and g not in genes and g in rna2gene:
        g = rna2gene[g]
    if g in genes:
        gex[g].extend(xs)
for g, d in genes.items():
    xs = sorted(gex.get(g) or [(d['s'], d['e'])])
    m = []
    for a, b in xs:
        if m and a <= m[-1][1] + 1:
            m[-1][1] = max(m[-1][1], b)
        else:
            m.append([a, b])
    d['ex'] = m
    d['exlen'] = sum(b - a + 1 for a, b in m)
    d['key'] = f"{chrom}:{d['s']}-{d['e']}"
G = sorted(genes, key=lambda g: genes[g]['s'])
log(f'[{tag}] {len(G)} gene nodes')
with open(P('nodes.gff3'), 'w') as f:
    f.write('##gff-version 3\n')
    for g in G:
        d = genes[g]
        f.write(f"{chrom}\t.\tgene\t{d['s']}\t{d['e']}\t.\t{d['strand']}\t.\tID={g};Name={d['name']}\n")
        for a, b in d['ex']:
            f.write(f"{chrom}\t.\texon\t{a}\t{b}\t.\t{d['strand']}\t.\tParent={g};gene={d['name']}\n")
with open(P('regions'), 'w') as f:
    f.write('\n'.join(genes[g]['key'] for g in G) + '\n')

# exon interval index: sorted (start, end, gene)
EXI = sorted((a, b, g) for g in G for a, b in genes[g]['ex'])
EXS = [x[0] for x in EXI]
BINW = 10000
EXB = collections.defaultdict(list)                # 10-kb bin -> exon intervals touching it
for iv in EXI:
    for k in range(iv[0] // BINW, iv[1] // BINW + 1):
        EXB[k].append(iv)


def ex_hits(a, b):
    """genes with exonic bases in [a,b] -> {gene: overlap}"""
    out_ = collections.Counter()
    seen = set()
    for k in range(a // BINW, b // BINW + 1):
        for iv in EXB.get(k, ()):
            if iv in seen:
                continue
            seen.add(iv)
            s, e, g = iv
            ov = min(b, e) - max(a, s) + 1
            if ov > 0:
                out_[g] += ov
    return out_


def ex_ivs(a, b):
    """exon intervals overlapping [a,b]"""
    out_ = set()
    for k in range(a // BINW, b // BINW + 1):
        for iv in EXB.get(k, ()):
            if iv[0] <= b and iv[1] >= a:
                out_.add(iv)
    return out_


# exon-overlapping gene pairs = one locus, never a read edge
same_locus = set()
for i, (a, b, g) in enumerate(EXI):
    j = i + 1
    while j < len(EXI) and EXI[j][0] <= b:
        if EXI[j][2] != g:
            same_locus.add(frozenset((g, EXI[j][2])))
        j += 1

# ---------- 2. reads: primary + secondary (drop unmapped/supplementary), pysam; cached ----------
import pysam
cache = P('reads.pkl')
if os.path.exists(cache):
    read_genes, read_recs, prim_reads = pickle.load(open(cache, 'rb'))
else:
    read_genes = collections.defaultdict(set)      # molecule -> genes
    read_recs = collections.defaultdict(list)      # molecule -> [genes of each record]
    prim_reads = collections.Counter()             # gene -> primary molecules (expression)
    n = 0
    for rd in pysam.AlignmentFile(bam).fetch(chrom):
        if rd.flag & 2052:
            continue
        hits = collections.Counter()
        for a, b in rd.get_blocks():
            for g, ov in ex_hits(a + 1, b).items():
                hits[g] += ov
        gs = {g for g, ov in hits.items() if ov >= MIN_EX_OV}
        if gs:
            read_genes[rd.query_name] |= gs
            read_recs[rd.query_name].append(frozenset(gs))
            if not rd.is_secondary:
                for g in gs:
                    prim_reads[g] += 1
        n += 1
    read_genes = dict(read_genes); read_recs = dict(read_recs)
    pickle.dump((read_genes, read_recs, prim_reads), open(cache, 'wb'))
    log(f'[{tag}] {n} records')
log(f'[{tag}] {len(read_genes)} molecules on exons')

gene_reads = collections.defaultdict(set)
for q, gs in read_genes.items():
    for g in gs:
        gene_reads[g].add(q)
shared = collections.Counter()      # any sharing: incl. ONE record spanning both genes (readthrough / adjacency)
shared_mm = collections.Counter()   # MULTIMAPPING only: g and h reached by DIFFERENT records, neither spanning both
multi = set()
for q, gs in read_genes.items():
    if len(gs) < 2:
        continue
    rs = read_recs[q]
    for a, b in itertools.combinations(sorted(gs), 2):
        if frozenset((a, b)) in same_locus:
            continue
        shared[(a, b)] += 1
        if any(a in r and b not in r for r in rs) and any(b in r and a not in r for r in rs):
            shared_mm[(a, b)] += 1
            multi.add(q)
log(f'[{tag}] {len(shared)} read-linked pairs (any); {len(shared_mm)} multimapper-linked pairs from {len(multi)} molecules')

# pass 2: anchors for linking molecules only: per record, per gene, read base (original orientation, 10-bp bins)
# -> reference position, sampled every BIN bp of the gene's exons
cache2 = P('anchors2.pkl')
if os.path.exists(cache2):
    recs = pickle.load(open(cache2, 'rb'))
else:
    recs = collections.defaultdict(list)
    rid = 0
    for rd in pysam.AlignmentFile(bam).fetch(chrom):
        if rd.flag & 2052 or rd.query_name not in multi:
            continue
        q = rd.query_name; rev = rd.is_reverse
        ct = rd.cigartuples
        hl = ct[0][1] if ct[0][0] == 5 else 0
        qlen = rd.infer_read_length()
        per = collections.defaultdict(dict)
        want = read_genes[q]
        for qp, rp in rd.get_aligned_pairs(matches_only=True):
            r1 = rp + 1
            if r1 % BIN:
                continue
            for s, e, g in EXB.get(r1 // BINW, ()):
                if s <= r1 <= e and g in want:
                    qq = qp + hl
                    per[g][((qlen - 1 - qq) if rev else qq) // 10] = r1
        rid += 1
        for g, m in per.items():
            recs[q].append((rid, g, rev, m))
    recs = dict(recs)
    pickle.dump(recs, open(cache2, 'wb'))
log(f'[{tag}] anchors extracted for {len(recs)} molecules')


def lis_len(pairs, decreasing):
    """size of the longest collinear subset; returns chained pairs"""
    pairs = sorted(pairs, key=lambda x: (x[0], -x[1] if not decreasing else x[1]))
    ys = [(-y if decreasing else y) for _, y in pairs]
    tails, tidx, prev = [], [], [-1] * len(ys)
    for i, y in enumerate(ys):
        k = bisect.bisect_left(tails, y)
        if k == len(tails):
            tails.append(y); tidx.append(i)
        else:
            tails[k] = y; tidx[k] = i
        prev[i] = tidx[k - 1] if k else -1
    out_, i = [], tidx[-1] if tidx else -1
    while i >= 0:
        out_.append(pairs[i]); i = prev[i]
    return out_


anchors = collections.defaultdict(lambda: ([], []))   # (a,b) -> (same-orient anchors, opposite)
for q, rl in recs.items():
    for (i1, g1, r1, m1), (i2, g2, r2, m2) in itertools.combinations(rl, 2):
        if i1 == i2 or g1 == g2 or frozenset((g1, g2)) in same_locus:
            continue
        if g1 > g2:
            g1, r1, m1, g2, r2, m2 = g2, r2, m2, g1, r1, m1
        tgt = anchors[(g1, g2)][0 if r1 == r2 else 1]
        for qb in m1.keys() & m2.keys():
            tgt.append((m1[qb], m2[qb]))

cov = {}
for (a, b), (same, opp) in anchors.items():
    best = []
    for arr, dec in ((same, genes[a]['strand'] != genes[b]['strand']), (opp, genes[a]['strand'] == genes[b]['strand'])):
        if arr:
            ch = lis_len(list(set(arr)), dec)
            if len(ch) > len(best):
                best = ch
    ca = len({x // BIN for x, _ in best}) * BIN / genes[a]['exlen']
    cb = len({y // BIN for _, y in best}) * BIN / genes[b]['exlen']
    small = ca if genes[a]['exlen'] <= genes[b]['exlen'] else cb
    cov[(a, b)] = (min(1.0, small), min(1.0, ca), min(1.0, cb))

log(f'[{tag}] chained {len(cov)} pairs')
# ---------- 3. shipped arm ----------
if not os.path.exists(P('paf')):
    subprocess.run(f"samtools faidx {fasta} -r {P('regions')} -o {P('loci.fa')}", shell=True, check=True)
    subprocess.run(f"minimap2 -x asm20 -c -X -N 50 -p 0.1 --secondary=yes -t 2 {P('loci.fa')} {P('loci.fa')} > {P('paf')} 2> {P('mm2.log')}",
                   shell=True, check=True)
subprocess.run(f"{R}/mcl_families --paf {P('paf')} --gff {P('nodes.gff3')} --min-exonic-bp 1 --min-shared-exon-frac 0.60 "
               f"--dump-graph {P('S.graph')} --out {P('S')} > {P('S.mcl.log')} 2>&1", shell=True, check=True)
key2g = {genes[g]['key']: g for g in G}
S_edges = {}
for ln in open(P('S.graph')):
    f = ln.rstrip('\n').split('\t')
    if len(f) == 3 and f[0] in key2g and f[1] in key2g:
        a, b = sorted((key2g[f[0]], key2g[f[1]]))
        S_edges[(a, b)] = float(f[2])

# ---------- 4. arms -> MCL -> clusters -> family_score ----------
tfam = collections.defaultdict(set)
for ln in list(open(truth))[1:]:
    nm, fam = ln.rstrip('\n').split('\t')[:2]
    tfam[nm].add(fam)
name2g = collections.defaultdict(list)
for g in G:
    name2g[genes[g]['name']].append(g)


def jac(a, b, mm=True):
    ra, rb = gene_reads[a], gene_reads[b]
    return (shared_mm if mm else shared)[(a, b)] / max(1, len(ra | rb))


arms = {'S': dict(S_edges)}
for t in (0.01, 0.05, 0.10, 0.20, 0.30):
    arms[f'J{t:.2f}'] = {k: jac(*k) for k in shared_mm if jac(*k) >= t}
    arms[f'Jany{t:.2f}'] = {k: jac(*k, mm=False) for k in shared if jac(*k, mm=False) >= t}
arms['C'] = {k: v[0] for k, v in cov.items() if v[0] >= COV}
arms['J0.05+C'] = {k: jac(*k) for k in arms['C'] if jac(*k) >= 0.05}
# DESCRIPTIVE (added after the dev run, not in the prereg): the anchor chain as a CERTIFICATE rather than a coverage
# floor -- keep a multimapper link only if the two alignments share read bases at all (cov > 0), or cover >= 0.30
arms['Cany'] = {k: v[0] for k, v in cov.items() if v[0] > 0}
arms['C0.30'] = {k: v[0] for k, v in cov.items() if v[0] >= 0.30}
arms['J0.01+Cany'] = {k: jac(*k) for k in arms['Cany'] if jac(*k) >= 0.01}
arms['S_or_C'] = dict(S_edges); arms['S_or_C'].update({k: v for k, v in arms['C'].items() if k not in S_edges})


def score(arm, clusters_path):
    r = subprocess.run(f"{R}/family_score --clusters {clusters_path} --gff {P('nodes.gff3')} --soto {truth} --chrom {chrom} "
                       f"--label {arm} --pairwise", shell=True, capture_output=True, text=True)
    return ' | '.join(l.strip() for l in r.stdout.strip().splitlines()[-2:])


def edge_prf(E):
    """edge-level, over truth genes only: an edge is TP if its two genes share a truth family"""
    tp = fp = 0
    for a, b in E:
        na, nb = genes[a]['name'], genes[b]['name']
        if na in tfam and nb in tfam:
            if tfam[na] & tfam[nb]:
                tp += 1
            else:
                fp += 1
    return tp, fp


truth_pairs = set()
fam2g = collections.defaultdict(set)
for nm, fs in tfam.items():
    for f in fs:
        for g in name2g.get(nm, []):
            fam2g[f].add(g)
for f, gs in fam2g.items():
    for a, b in itertools.combinations(sorted(gs), 2):
        truth_pairs.add((a, b))

def emit(fo, arm, E, cl):
    cp = P(f'{arm}.clusters.tsv')
    with open(cp, 'w') as f:
        f.write('cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n')
        for i, line in enumerate(l for l in cl.splitlines() if l.strip()):
            mem = line.split('\t')
            if len(mem) < 2:
                continue
            for k in mem:
                c, se = k.rsplit(':', 1); s, e = se.split('-')
                f.write(f'RC{i}\t{len(mem)}\tNA\tNA\tNA\t{c}\t{s}\t{e}\n')
    tp, fp = edge_prf(E)
    hit = len(truth_pairs & set(E))
    sc = score(arm, cp)
    fo.write(f'{arm}\t{len(E)}\t{tp}\t{fp}\t{hit}\t{sc}\n')
    log(f'{arm:12s} edges {len(E):6d} TP {tp:5d} FP {fp:5d} truth-pairs {hit}/{len(truth_pairs)} || {sc}')


def components(E):
    par = {}

    def fd(x):
        while par.setdefault(x, x) != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for a, b in E:
        par[fd(a)] = fd(b)
    comp = collections.defaultdict(list)
    for x in list(par):
        comp[fd(x)].append(genes[x]['key'])
    return '\n'.join('\t'.join(v) for v in comp.values())


with open(P('summary.tsv'), 'w') as fo:
    fo.write('arm\tedges\tedge_TP\tedge_FP\ttruth_pairs_hit\tfamily_score\n')
    for arm, E in arms.items():
        gp = P(f'{arm}.graph')
        with open(gp, 'w') as f:
            for (a, b), w in E.items():
                f.write(f"{genes[a]['key']}\t{genes[b]['key']}\t{max(w, 1e-6):.6f}\n")
        emit(fo, arm + '.mcl', E, subprocess.run(f"{R}/mcl_port --graph {gp}", shell=True, capture_output=True,
                                                 text=True, check=True).stdout)
        emit(fo, arm + '.cc', E, components(E))       # operator-free control
    sc = score('S_shipped', P('S.clusters.tsv'))
    tp, fp = edge_prf(S_edges)
    fo.write(f'S_shipped\t{len(S_edges)}\t{tp}\t{fp}\t{len(truth_pairs & set(S_edges))}\t{sc}\n')
    log(f'S_shipped  {sc}')

# ---------- 5. ceilings: which truth pairs could a read rule ever see? ----------
tp_all = len(truth_pairs)
both_expr = [k for k in truth_pairs if prim_reads[k[0]] >= 3 and prim_reads[k[1]] >= 3]
any_shared = [k for k in truth_pairs if shared_mm.get(k, 0) >= 1]
ge3_shared = [k for k in truth_pairs if shared_mm.get(k, 0) >= 3]
rt_only = [k for k in truth_pairs if shared.get(k, 0) >= 1 and not shared_mm.get(k, 0)]
s_hit = [k for k in truth_pairs if k in S_edges]
with open(P('ceiling.txt'), 'w') as f:
    msg = (f'{tag}: truth same-family gene pairs {tp_all}\n'
           f'  both genes expressed (>=3 primary molecules)  {len(both_expr)} ({len(both_expr)/tp_all:.3f})\n'
           f'  >=1 shared molecule, ANY (incl. readthrough)  {len([k for k in truth_pairs if shared.get(k, 0)])}\n'
           f'    of which ONLY via a record spanning both    {len(rt_only)}\n'
           f'  >=1 MULTIMAPPING molecule (distinct records)  {len(any_shared)} ({len(any_shared)/tp_all:.3f})\n'
           f'  >=3 multimapping molecules                    {len(ge3_shared)} ({len(ge3_shared)/tp_all:.3f})\n'
           f'  shipped DNA edge (pre-MCL)                    {len(s_hit)} ({len(s_hit)/tp_all:.3f})\n'
           f'  shipped edge AND >=1 shared molecule          {len(set(s_hit) & set(any_shared))}\n'
           f'  shared molecule but NO shipped edge           {len(set(any_shared) - set(s_hit))}\n'
           f'  read-linked NON-family pairs (both in truth)  '
           f'{sum(1 for (a, b) in shared_mm if genes[a]["name"] in tfam and genes[b]["name"] in tfam and not tfam[genes[a]["name"]] & tfam[genes[b]["name"]])}\n')
    f.write(msg)
log(msg)
# per-pair dump for the doc
with open(P('pairs.tsv'), 'w') as f:
    f.write('a\tb\tshared_any\tshared_mm\tjaccard_mm\tcov_small\tcov_a\tcov_b\tS_edge\ttruth_same\n')
    for k in set(shared) | set(S_edges):
        a, b = k
        c = cov.get(k, (0, 0, 0))
        na, nb = genes[a]['name'], genes[b]['name']
        same = 'NA' if not (na in tfam and nb in tfam) else int(bool(tfam[na] & tfam[nb]))
        f.write(f"{na}\t{nb}\t{shared.get(k, 0)}\t{shared_mm.get(k, 0)}\t{jac(a, b):.4f}\t{c[0]:.3f}\t{c[1]:.3f}\t{c[2]:.3f}\t"
                f"{int(k in S_edges)}\t{same}\n")
