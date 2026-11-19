#!/usr/bin/env python3
"""PREREG r1777 (docs/PREREG_allele_graph_2026-11-17.md): score near-tied reads on locus graphs = exons of the produced GTF
+ the alleles the sample's own confidently placed reads show there; judged with the read table of neutral_truth.py.

  allele_graph.py regions READS.tsv OUT.bed
  allele_graph.py score READS.tsv SEL.bam REGION.bam PRODUCED.gtf
        SEL.bam = all records of the table's reads (neutral_truth's MM.bam); REGION.bam = every record of the sample on the regions"""
import sys, os, bisect, collections
import numpy as np, pysam

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan"); GUARD = int(os.environ.get("AG_GUARD", 5000)); GAIN = 5  # AG_GUARD: r1829, default unchanged

def table(path):
    """scored reads of the table: all S2 + the first GUARD S1; read -> dict(st, status, primary, origin, loci [(chrom, s, e, strand, AS)])"""
    T = {}; n1 = 0
    for l in open(path):
        f = l.rstrip('\n').split('\t')
        if f[0] == 'read' or f[1] == 'S0': continue
        if f[1] == 'S1':
            n1 += 1
            if n1 > GUARD: continue
        loci = []
        for x in f[6].split(';'):
            i, c, se, st, a = x.rsplit(':', 4); s, e = se.split('-'); loci.append((c, int(s), int(e), st, int(a[2:])))
        T[f[0]] = dict(st=f[1], status=f[3], primary=int(f[4].split(';')[0].split('=')[1]), origin={int(i) for i in f[5].split(',') if i}, loci=loci)
    return T

def regions(tsv, out):
    R = collections.defaultdict(list)
    for d in table(tsv).values():
        for c, s, e, st, a in d['loci']:
            if a >= 0.8 * d['loci'][0][4]: R[c].append((s, e))
    n = 0
    with open(out, 'w') as w:
        for c in sorted(R):
            cur = None
            for s, e in sorted(R[c]):
                if cur and s <= cur[1]: cur[1] = max(cur[1], e)
                else:
                    if cur: w.write(f"{c}\t{cur[0]}\t{cur[1]}\n"); n += 1
                    cur = [s, e]
            w.write(f"{c}\t{cur[0]}\t{cur[1]}\n"); n += 1
    print(f"{n} regions")

def mismatches(start, cigar, seq):
    """(reference position, read base) of every X column; seq = the whole read in the alignment's orientation, clips included."""
    out = []; p = start; q = 0
    for op, n in cigar:
        if op == 8: out += [(p + i, seq[q + i]) for i in range(n)]
        if op in (0, 7, 8): p += n; q += n
        elif op in (1, 4, 5): q += n
        elif op in (2, 3): p += n
    return out
assert mismatches(10, [(4, 2), (7, 3), (8, 1), (3, 100), (1, 1), (7, 2), (8, 2)], 'NNAAAGTCCTA') == [(13, 'G'), (116, 'T'), (117, 'A')]

def score(tsv, sel, region, gtf):
    T = table(tsv)
    nodes = collections.defaultdict(list)  # chrom -> merged exon intervals of the produced GTF
    for l in open(gtf):
        f = l.split('\t')
        if len(f) > 4 and f[2] == 'exon': nodes[f[0]].append((int(f[3]) - 1, int(f[4])))
    for c in nodes:
        m = []
        for s, e in sorted(nodes[c]):
            if m and s <= m[-1][1]: m[-1][1] = max(m[-1][1], e)
            else: m.append([s, e])
        nodes[c] = (np.array([x[0] for x in m]), np.array([x[1] for x in m]))
    def in_node(c, p):
        if c not in nodes: return False
        i = bisect.bisect_right(nodes[c][0], p) - 1
        return i >= 0 and p < nodes[c][1][i]
    cnt = collections.Counter(); bs = collections.defaultdict(list); be = collections.defaultdict(list); na = 0
    for r in pysam.AlignmentFile(region):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 60 or r.query_sequence is None: continue
        na += 1; c = r.reference_name
        for p, b in mismatches(r.reference_start, r.cigartuples, r.query_sequence): cnt[(c, p, b)] += 1
        for a, b in r.get_blocks(): bs[c].append(a); be[c].append(b)
    alleles = set()
    for c in bs:  # anchor coverage at the candidate positions only
        P = np.array(sorted({p for (cc, p, b) in cnt if cc == c}))
        if not len(P): continue
        d = np.zeros(len(P) + 1, dtype=np.int64)
        np.add.at(d, np.searchsorted(P, np.array(bs[c])), 1); np.add.at(d, np.searchsorted(P, np.array(be[c])), -1)
        cov = dict(zip(P.tolist(), np.cumsum(d)[:-1].tolist()))
        for (cc, p, b), n in cnt.items():
            if cc == c and n >= 3 and n >= 0.2 * cov[p] and in_node(c, p): alleles.add((c, p, b))
    print(f"anchor reads {na}; exon nodes {sum(len(v[0]) for v in nodes.values())}; substitution columns seen {len(cnt)}; alleles on the graphs {len(alleles)}")
    recs = collections.defaultdict(list); seqs = {}
    for r in pysam.AlignmentFile(sel, check_sq=False):
        if r.is_unmapped or r.is_supplementary or r.query_name not in T: continue
        recs[r.query_name].append((r.reference_name, '-' if r.is_reverse else '+', r.reference_start, r.reference_end, r.get_tag('AS'), r.cigartuples))
        if not r.is_secondary: seqs[r.query_name] = r.query_sequence.translate(COMP)[::-1] if r.is_reverse else r.query_sequence
    C = collections.Counter(); bad = 0
    for q, d in T.items():
        fw = seqs[q]; rv = fw.translate(COMP)[::-1]; sc = []
        for c, s, e, st, a in d['loci']:  # the locus's best-AS record, as neutral_truth merged them
            rr = [x for x in recs[q] if x[0] == c and x[1] == st and x[2] < e and s < x[3]]
            if not rr: bad += 1; sc.append(a); continue
            x = max(rr, key=lambda x: x[4])
            sc.append(a + GAIN * sum((c, p, b) in alleles for p, b in mismatches(x[2], x[5], fw if st == '+' else rv)))
        o = sorted(range(len(sc)), key=lambda i: -sc[i]); lead = sc[o[0]] - sc[o[1]] if len(o) > 1 else 10 ** 6
        for arm, call in (('primary', d['primary']), ('graph', o[0] if lead > 10 else None), ('graph5', o[0] if lead >= GAIN else None)):
            k = (d['st'], d['status'], arm)
            if call is None: C[k + ('abstain',)] += 1; call = d['primary'] if arm != 'primary' else None
            else:
                C[k + ('right' if call in d['origin'] else 'wrong',)] += 1
                if call != d['primary']: C[k + ('moved', call in d['origin'])] += 1
            if arm != 'primary': C[(d['st'], d['status'], arm + '-or-primary', 'right' if call in d['origin'] else 'wrong')] += 1
    print(f"loci without a record in SEL.bam (scored at their AS): {bad}")
    for st in ('S2', 'S1'):
        n = sum(1 for d in T.values() if d['st'] == st and d['status'] == 'decidable')
        print(f"\n{st}: {n} decidable reads\narm\tright\twrong\tabstain\twrong_share_of_calls\tmoved_off_primary_right\tmoved_off_primary_wrong")
        for arm in ('primary', 'graph', 'graph-or-primary', 'graph5', 'graph5-or-primary'):
            r, x, ab = (C[(st, 'decidable', arm, k)] for k in ('right', 'wrong', 'abstain'))
            print(f"{arm}\t{r}\t{x}\t{ab}\t{x / (r + x) if r + x else 0:.4f}\t{C[(st, 'decidable', arm, 'moved', True)]}\t{C[(st, 'decidable', arm, 'moved', False)]}")
        print("class\treads\tgraph_calls\tgraph_moves_off_primary\tgraph5_calls")
        for cl in ('outside', 'tie', 'no near-perfect hit'):
            n = sum(1 for d in T.values() if d['st'] == st and d['status'] == cl)
            g = n - C[(st, cl, 'graph', 'abstain')]; g5 = n - C[(st, cl, 'graph5', 'abstain')]
            print(f"{cl}\t{n}\t{g}\t{C[(st, cl, 'graph', 'moved', True)] + C[(st, cl, 'graph', 'moved', False)]}\t{g5}")

if __name__ == '__main__':
    regions(*sys.argv[2:4]) if sys.argv[1] == 'regions' else score(*sys.argv[2:6])
