#!/usr/bin/env python3
"""Tied reads decided for a CLUSTER of near-identical loci (register 1974 / 1975; stage o2, RUSTLE_O2_CLUSTERS=1).

    cluster_verdict.py ASSEMBLY.gtf GENOME.fa TIES.tsv.gz OUT_PREFIX [--identity 0.99] [--lead 10] [--threads N]

The Eichler lab's 10-point filter (Guitart 2024) compares paralog CLUSTERS and ignores the reads that fall between them.
Here, on one reference: the loci that tied reads touch are compared with each other (the spliced sequence of each
locus's longest transcript, `minimap2 -cx asm20 -X`) and joined when they align over >= half of the shorter at
--identity or more (single linkage). A tied read is DECIDED for a cluster when every alignment within --lead `ms` of its
best lies in that one cluster, AMBIGUOUS otherwise (also when one of them lies on no locus).
TIES.tsv.gz = tools/tie_table.py's table (read, locus, gap, best_ms) with its .placed.tsv beside it.
Writes OUT_PREFIX.clusters.tsv (per cluster: loci, reads placed alone, tied reads decided, ambiguous reads touching it)
and OUT_PREFIX.cluster_loci.tsv (locus, cluster, loci in the cluster, and the same three counts for the locus).
A verdict is "one of N loci": N is in both tables and belongs to every statement made from them.
On HG002, KB3781 and planted reads this decides 74-85% of the reads tied between loci at 96.6-100% right; per locus the
same lead decides 8-16% (docs/PREREG_cluster_filter_denovo_2026-11-26.md). Needs pysam and minimap2 on PATH."""
import sys, os, re, gzip, argparse, subprocess, collections, tempfile
OUTSIDE = 'outside'
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

def read_ties(path):
    """-> {read: {locus: smallest gap}}"""
    t = collections.defaultdict(dict)
    with gzip.open(path, 'rt') as f:
        for l in f:
            if l[0] == '#' or l.startswith('read\t'): continue
            q, g, gap, _ = l.rstrip('\n').split('\t'); gap = int(gap)
            if gap < t[q].get(g, 10 ** 9): t[q][g] = gap
    return t

def longest(gtf, want):
    """-> {locus: (chrom, strand, exons of its transcript with the most exon bases)} for the loci in `want`"""
    ex = collections.defaultdict(list); meta = {}
    for l in open(gtf):
        f = l.split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        g = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        if g not in want: continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); ex[(g, t)].append((int(f[3]) - 1, int(f[4]))); meta[g] = (f[0], f[6])
    best = {}
    for (g, t), e in ex.items():
        n = sum(b - a for a, b in e)
        if n > best.get(g, (0,))[0]: best[g] = (n, sorted(e))
    return {g: (meta[g][0], meta[g][1], e) for g, (n, e) in best.items()}

def join(loci, paf_lines, identity):
    """single linkage over PAF lines -> {locus: cluster representative}"""
    p = {}
    def find(x):
        while p.setdefault(x, x) != x: p[x] = p[p[x]]; x = p[x]
        return x
    for l in paf_lines:
        f = l.split('\t'); a, b = f[0], f[5]
        if a != b and int(f[10]) >= 0.5 * min(int(f[1]), int(f[6])) and int(f[9]) / max(1, int(f[10])) >= identity:
            ra, rb = sorted((find(a), find(b))); p[rb] = ra
    return {g: find(g) for g in loci}

def verdict(rows, cluster, lead):
    """rows {locus: gap} of one read -> the cluster it is decided for, or None (ambiguous)"""
    cs = {cluster.get(g) if g != OUTSIDE else None for g, gap in rows.items() if gap < lead}
    return next(iter(cs)) if len(cs) == 1 and None not in cs else None

def run(gtf, fasta, ties, out, identity=0.99, lead=10, threads=8):
    import pysam
    T = read_ties(ties); placed = collections.Counter()
    for l in list(open(ties + '.placed.tsv'))[1:]:
        g, n = l.rstrip('\n').split('\t'); placed[g] = int(n)
    want = {g for r in T.values() for g in r if g != OUTSIDE}; L = longest(gtf, want); fa = pysam.FastaFile(fasta)
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(out)) or '.') as d:
        with open(d + '/loci.fa', 'w') as w:
            for g, (c, s, e) in L.items():
                seq = ''.join(fa.fetch(c, a, b) for a, b in e); w.write(f'>{g}\n{seq[::-1].translate(COMP) if s == "-" else seq}\n')
        r = subprocess.run(['minimap2', '-cx', 'asm20', '-X', '--eqx', '-N', '100', '-p', '0.3', '-t', str(threads), d + '/loci.fa', d + '/loci.fa'], capture_output=True, text=True)
        if r.returncode != 0: sys.exit(f'cluster_verdict: minimap2 failed: {r.stderr[-300:]}')
        cl = join(want, r.stdout.splitlines(), identity)
    members = collections.defaultdict(list)
    for g, c in cl.items(): members[c].append(g)
    dec = collections.Counter(); amb = collections.Counter(); ldec = collections.Counter(); lamb = collections.Counter(); nd = na = 0; one = 0
    for q, rows in T.items():
        v = verdict(rows, cl, lead); near = [g for g, gap in rows.items() if gap < lead and g != OUTSIDE]
        if v is None:
            na += 1
            for g in near: lamb[g] += 1
            for c in {cl[g] for g in near}: amb[c] += 1
        else:
            nd += 1; dec[v] += 1; one += len(near) == 1
            for g in near: ldec[g] += 1
    ids = {c: f'C{i + 1}' for i, c in enumerate(sorted(members, key=lambda c: (-len(members[c]), c)))}
    with open(out + '.clusters.tsv.tmp', 'w') as w:
        w.write('cluster\tloci\treads_placed_alone\ttied_reads_decided\tambiguous_reads_touching\tmembers\n')
        for c in sorted(members, key=lambda c: (-len(members[c]), c)):
            m = sorted(members[c]); w.write(f"{ids[c]}\t{len(m)}\t{sum(placed[g] for g in m)}\t{dec[c]}\t{amb[c]}\t{','.join(m)}\n")
    with open(out + '.cluster_loci.tsv.tmp', 'w') as w:
        w.write('locus\tcluster\tloci_in_cluster\treads_placed_alone\ttied_reads_decided_for_the_cluster\tambiguous_reads\n')
        for g in sorted(cl): w.write(f'{g}\t{ids[cl[g]]}\t{len(members[cl[g]])}\t{placed[g]}\t{ldec[g]}\t{lamb[g]}\n')
    os.replace(out + '.clusters.tsv.tmp', out + '.clusters.tsv'); os.replace(out + '.cluster_loci.tsv.tmp', out + '.cluster_loci.tsv')
    sz = sorted(len(members[v]) for v in dec.elements()); med = sz[len(sz) // 2] if sz else 0
    print(f'cluster_verdict: {len(T)} tied reads: {nd} decided for a cluster ({nd / max(1, len(T)):.3f}; {one} of them for a single locus), {na} ambiguous; '
          f'{len(members)} clusters over {len(cl)} loci, {sum(1 for m in members.values() if len(m) > 1)} with several loci, largest {max((len(m) for m in members.values()), default=0)}; '
          f'a decided read names {med} loci in the median -> {out}.clusters.tsv, .cluster_loci.tsv')

def selftest():
    cl = join(['A', 'B', 'C', 'D'], ['A\t1000\t0\t1000\t+\tB\t1000\t0\t1000\t995\t1000\t60', 'A\t1000\t0\t400\t+\tC\t1000\t0\t400\t400\t400\t60', 'C\t1000\t0\t1000\t+\tD\t1000\t0\t1000\t980\t1000\t60'], 0.99)
    assert cl['A'] == cl['B'] == 'A' and cl['C'] == 'C' and cl['D'] == 'D', cl          # short hit and 98% hit join nothing
    assert verdict({'A': 0, 'B': 0}, cl, 10) == 'A' and verdict({'A': 0, 'C': 4}, cl, 10) is None and verdict({'A': 0, 'C': 10}, cl, 10) == 'A'
    assert verdict({'A': 0, OUTSIDE: 3}, cl, 10) is None and verdict({'A': 0, OUTSIDE: 10}, cl, 10) == 'A' and verdict({OUTSIDE: 0}, cl, 10) is None
    print('cluster_verdict selftest ok')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: selftest()
    else:
        ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
        ap.add_argument('gtf'); ap.add_argument('fasta'); ap.add_argument('ties'); ap.add_argument('out')
        ap.add_argument('--identity', type=float, default=0.99); ap.add_argument('--lead', type=int, default=10); ap.add_argument('--threads', type=int, default=8)
        a = ap.parse_args(); run(a.gtf, a.fasta, a.ties, a.out, a.identity, a.lead, a.threads)
