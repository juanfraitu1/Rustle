#!/usr/bin/env python3
"""More expressed sequences than two alleles per locus: a sign of a copy the reference lacks (register 1983).

    locus_sequences.py ASSEMBLY.gtf BAM AS_TABLE GENOME.fa OUT.tsv [--min-seq-reads 3] [--alpha 0.01] [--names] [--threads 8]

recurrence_test.py (register 1926, IsoCon's test) reads only the tied reads of a group. Here EVERY read of a locus is
read: its primaries plus the secondaries exactly tied on ms with their read's best. Per locus and isoform (intron chain)
the same decision gives the distinct sequences; a sequence counts with >= --min-seq-reads reads. Loci joined by exactly
tied reads form a group (tie_table.groups; a locus with no tie is its own group). Sequences of different loci or
isoforms that share at least half of the smaller one's reads are one sequence (a tied read is read at each of its loci).
Cells (locus, isoform) joined that way form a class, with slots = 2 x the loci its cells are on (two alleles per
reference copy; the loci of the whole group would be too generous: groups chain through a few tied reads).
--count conflicts (register 1985, opt-in): per cell the count is the largest set of sequences that pairwise disagree at a
site both cover (a shorter isoform or shorter reads of an allele agree with it and add nothing).
flag = some class of the group has more sequences than slots; `sequences` and `slots` are those of the class with the
largest excess. OUT.tsv: group, loci, n_loci, slots, sequences, flag, primaries, per_locus (locus:isoform:sequences), and with
--names the names of the group's primary reads (planted sets: the scorer reads the truth from them).
It cannot tell a third allele of a hidden copy from a copy with three real alleles, and says nothing about which copy.
# shortcut: one serial pass over the BAM and names in memory, fine for planted sets; parallel per contig before a whole sample."""
import sys, os, argparse, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
import eq_classes as E, tie_table as TT, recurrence_test as R

def count_classes(cands):
    """cands [(cell, set of reads)], cell = (locus, isoform) -> {cell: (distinct sequences in the cell's class, loci of the
    class)}; two candidates are one sequence when they share >= half of the smaller one's reads; cells holding one
    sequence are one class"""
    p = list(range(len(cands)))
    def find(x):
        while p[x] != x: p[x] = p[p[x]]; x = p[x]
        return x
    by = collections.defaultdict(list)
    for i, (_, rs) in enumerate(cands):
        for q in rs: by[q].append(i)
    share = collections.Counter()
    for v in by.values():
        for i in v:
            for j in v:
                if i < j: share[(i, j)] += 1
    for (i, j), n in share.items():
        if 2 * n >= min(len(cands[i][1]), len(cands[j][1])): p[find(j)] = find(i)
    cp = {}
    def cfind(c):
        while cp.setdefault(c, c) != c: cp[c] = cp[cp[c]]; c = cp[c]
        return c
    comp = collections.defaultdict(list)
    for i in range(len(cands)): comp[find(i)].append(i)
    for v in comp.values():
        for i in v[1:]: cp[cfind(cands[i][0])] = cfind(cands[v[0]][0])
    n = collections.Counter(cfind(cands[v[0]][0]) for v in comp.values()); loci = collections.defaultdict(set)
    for c, _ in cands: loci[cfind(c)].add(c[0])
    return {c: (n[cfind(c)], len(loci[cfind(c)])) for c, _ in cands}

def max_conflict(pats):
    """pats [{site: allele}] -> size of the largest set of sequences that pairwise disagree at a site both cover (register
    1985). Two sequences that agree wherever both are covered can be one haplotype seen through two isoforms or two read
    lengths; only sequences that contradict each other must be different haplotypes."""
    n = len(pats); adj = [{j for j in range(n) if j != i and any(pats[j].get(x, v) != v for x, v in pats[i].items())} for i in range(n)]
    best = 0
    def grow(size, cand):
        nonlocal best
        best = max(best, size)
        if size + len(cand) <= best: return
        for v in sorted(cand): grow(size + 1, {u for u in cand if u > v and u in adj[v]})
    grow(0, set(range(n))); return best

def main():
    if sys.argv[1:] == ['--selftest']:
        h1, h2 = {1: True, 2: False}, {1: False, 2: True}
        assert max_conflict([h1, h2, {2: False, 3: True}, {2: True, 3: False}]) == 2          # each allele seen through two isoforms
        assert max_conflict([h1, h2, {1: True, 2: True}]) == 3 and max_conflict([h1, h2, {3: True}]) == 2 and max_conflict([]) == 0
        a, b, h = {'a1', 'a2', 'a3'}, {'b1', 'b2', 'b3'}, {'h1', 'h2', 'h3', 'h4'}
        c = count_classes([(('A', 1), a), (('A', 1), a | {'x'} - {'a1'}), (('A', 1), h), (('B', 1), h - {'h4'}), (('B', 1), b), (('C', 2), {'c1', 'c2', 'c3'})])
        assert c[('A', 1)] == (3, 2) and c[('B', 1)] == (3, 2) and c[('C', 2)] == (1, 1), c          # A and B share the hidden sequence: one class of 3 (a, h, b)
        print('locus_sequences selftest ok'); return
    ap = argparse.ArgumentParser(); [ap.add_argument(x) for x in ('gtf', 'bam', 'tab', 'fasta', 'out')]
    ap.add_argument('--min-seq-reads', type=int, default=3); ap.add_argument('--alpha', type=float, default=0.01); ap.add_argument('--names', action='store_true'); ap.add_argument('--count', choices=('sequences', 'conflicts'), default='sequences', help='conflicts (register 1985): per cell, the largest set of sequences that pairwise disagree, instead of the number of sequences'); ap.add_argument('--dump', metavar='FILE', help='every counted sequence: locus, isoform, reads, alleles, read names (register 1985; changes no decision)'); ap.add_argument('--threads', type=int, default=8)
    a = ap.parse_args()
    import pysam
    meta, sets, gl, _ = TT.tie_groups(a.gtf, a.bam, a.tab, a.threads)
    grp = {g: k for k, v in gl.items() for g in v}
    reads = collections.defaultdict(dict); ln = {}; prim = collections.defaultdict(list); B = pysam.AlignmentFile(a.bam); idx = {}
    for r in B:
        if r.is_unmapped or r.is_supplementary: continue
        q = r.query_name
        if r.is_secondary and (E.G_TIES.get(q) is None or r.get_tag('ms') != E.G_TIES[q]): continue
        c, s = r.reference_name, E.tx_strand(r)
        if (c, s) not in E.G_GEX: continue
        if (c, s) not in idx: idx[(c, s)] = ([x for x, _, _ in E.G_GEX[(c, s)]], max(e - x for x, e, _ in E.G_GEX[(c, s)]))
        ov = E.overlap_genes([(x + 1, e) for x, e in r.get_blocks()], E.G_GEX[(c, s)], *idx[(c, s)])
        if not ov: continue
        g = max(ov, key=lambda k: (ov[k], k))
        if not r.is_secondary: prim[g].append(q)
        if q not in reads[g] or r.query_alignment_length > ln[(g, q)]: reads[g][q] = (*R.variants(r), R.introns(r)); ln[(g, q)] = r.query_alignment_length
    F = pysam.FastaFile(a.fasta); cands = []; cells = collections.defaultdict(list); dump = []; pats = collections.defaultdict(list)
    for g in sorted(reads):
        rd = reads[g]; c = meta[g][0]; mult = {}
        if len(rd) < a.min_seq_reads: continue
        for x in {x for v in rd.values() for x in v[0] if x[1] in 'ID'}:          # homopolymer length at each indel, as recurrence_test._collect_group
            seq = F.fetch(c, max(0, x[0] - 20), x[0] + 21).upper(); i = min(20, x[0]); b = seq[i] if i < len(seq) else 'N'; u = 1; j = i + 1
            while j < len(seq) and seq[j] == b: u += 1; j += 1
            j = i - 1
            while j >= 0 and seq[j] == b: u += 1; j -= 1
            mult[x] = u
        iso, _ = R.by_isoform({q: v[2] for q, v in rd.items()})
        for ci, ch in enumerate(sorted(iso), 1):
            cs, _, _ = R.decide({q: rd[q][:2] for q in iso[ch]}, mult, a.alpha, 5, 'fold', None, None, 1000)
            for x in cs or ():
                if len(x['reads']) >= a.min_seq_reads:
                    cands.append(((g, ci), set(x['reads']))); cells[g].append(ci); pats[(g, ci)].append(x['pat'])
                    dump.append(f"{g}\t{ci}\t{len(iso[ch])}\t{len(x['reads'])}\t{','.join(f'{y[0] + 1}{y[1]}{y[2]}' for y in sorted(x['pat']) if x['pat'][y]) or 'reference'}\t{','.join(f'{y[0] + 1}{y[1]}{y[2]}' for y in sorted(x['pat']))or '.'}\t{'NA' if x['p'] is None else format(x['p'], '.3g')}\t{','.join(sorted(x['reads']))}\n")
    if a.dump:
        with open(a.dump, 'w') as w: w.write('locus\tisoform\tisoform_reads\treads\talleles\tsites_covered\tp\tread_names\n'); w.writelines(dump)
    n = count_classes(cands); out = collections.defaultdict(list)
    for g in sorted(set(reads) | set(prim)): out[grp.get(g, g)].append(g)
    flagged = 0
    with open(a.out + '.tmp', 'w') as w:
        w.write('group\tloci\tn_loci\tslots\tsequences\tflag\tprimaries\tper_locus' + ('\tprimary_names' if a.names else '') + '\n')
        for k in sorted(out):
            ls = out[k]; nl = len(gl.get(k, ls)); per = [(g, ci, *n[(g, ci)]) for g in ls for ci in sorted(set(cells[g]))]
            if a.count == 'conflicts': per = [(g, ci, max_conflict(pats[(g, ci)]), m) for g, ci, _, m in per]
            s, sl = max(((x[2], 2 * x[3]) for x in per), key=lambda v: (v[0] - v[1], v[0]), default=(0, 2)); f = int(s > sl); flagged += f
            w.write(f"{k}\t{','.join(ls)}\t{nl}\t{sl}\t{s}\t{f}\t{sum(len(prim[g]) for g in ls)}\t{','.join(f'{g}:{ci}:{m}' for g, ci, m, _ in per) or '.'}" + ('\t' + ','.join(q for g in ls for q in prim[g]) if a.names else '') + '\n')
    os.replace(a.out + '.tmp', a.out)
    print(f'locus_sequences: {len(out)} groups, {flagged} with more sequences than two per locus -> {a.out}')

if __name__ == '__main__': main()
