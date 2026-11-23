#!/usr/bin/env python3
"""How many distinct sequences are expressed among copies the reference cannot tell apart (register 1924).

    sample_psv.py ASSEMBLY.gtf BAM AS_TABLE OUT.tsv [--min-reads 3] [--min-share 0.1] [--threads N]

Loci and exact `ms` ties as tools/eq_classes.py. Loci joined by tied reads form a group. For each group the tied reads
are read at ONE locus of it (the one most of them reach; minimap2 -Y keeps the sequence in secondary records), and the
positions where they differ from the reference are collected. A difference seen in >= --min-reads reads and >= --min-share
of the reads covering it is a site of the SAMPLE (an error does not recur). Reads are then grouped greedily into
haplotypes that agree at every site both cover; haplotypes of >= --min-reads reads are counted.
OUT.tsv: group, loci, locus read at, tied reads, sites, haplotypes, haplotype sizes, reads per haplotype (names).
No read is assigned to a locus: the count says how many different copies are expressed, not which."""
import sys, os, argparse, collections, multiprocessing
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
import eq_classes as E

def variants(r):
    """-> set of (pos, kind, alt) on the reference, list of aligned blocks (0-based half open)"""
    out = set(); rp, qp = r.reference_start, 0; seq = r.query_sequence
    for op, n in r.cigartuples:
        if op in (0, 7): rp += n; qp += n
        elif op == 8:
            for i in range(n): out.add((rp + i, 'S', seq[qp + i] if seq else 'N'))
            rp += n; qp += n
        elif op == 1: out.add((rp, 'I', seq[qp:qp + n] if seq else str(n))); qp += n
        elif op == 2: out.add((rp, 'D', n)); rp += n
        elif op == 3: rp += n
        elif op == 4: qp += n
    return out, r.get_blocks()

def haplotypes(reads, min_reads, min_share):
    """reads {name: (variants, blocks)} -> sites, [haplotype read lists]"""
    cov = collections.Counter(); sup = collections.Counter()
    def covered(blocks, pos): return any(a <= pos < b for a, b in blocks)
    for v, bl in reads.values():
        for x in v: sup[x] += 1
    cand = [x for x, n in sup.items() if n >= min_reads]
    sites = []
    for x in cand:
        c = sum(1 for v, bl in reads.values() if covered(bl, x[0]))
        if c and sup[x] / c >= min_share: sites.append(x)
    sites.sort()
    vec = {}
    for q, (v, bl) in reads.items():
        vec[q] = {x: (x in v) for x in sites if covered(bl, x[0])}
    haps = []          # [consensus dict, [reads]]
    for q in sorted(vec, key=lambda q: (-len(vec[q]), q)):
        for h in haps:
            if all(h[0].get(x, a) == a for x, a in vec[q].items()): h[0].update(vec[q]); h[1].append(q); break
        else: haps.append([dict(vec[q]), [q]])
    return sites, [h[1] for h in haps]

def selftest():
    bl = [(0, 100)]; R = {}
    for i in range(10): R[f'a{i}'] = ({(10, 'S', 'T')}, bl)
    for i in range(8): R[f'b{i}'] = ({(50, 'D', 1)}, bl)
    R['e'] = ({(10, 'S', 'T'), (70, 'S', 'G')}, bl)          # an error: not a site
    R['p'] = (set(), [(60, 100)])                             # a partial read covering no site: joins the first haplotype
    s, h = haplotypes(R, 3, 0.1)
    assert s == [(10, 'S', 'T'), (50, 'D', 1)], s
    assert sorted(len(x) for x in h if len(x) >= 3) == [8, 12], [len(x) for x in h]
    s, h = haplotypes({f'r{i}': ({(5, 'S', 'A')} if i == 0 else set(), bl) for i in range(20)}, 3, 0.1)
    assert s == [] and len(h) == 1
    print('sample_psv selftest ok')

def main():
    if sys.argv[1:] == ['--selftest']: return selftest()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('gtf'); ap.add_argument('bam'); ap.add_argument('tab'); ap.add_argument('out')
    ap.add_argument('--min-reads', type=int, default=3); ap.add_argument('--min-share', type=float, default=0.1); ap.add_argument('--threads', type=int, default=8)
    a = ap.parse_args()
    import pysam
    meta, E.G_GEX = E.load_loci(a.gtf)
    E.G_FLOOR = E.floor_fn(10); E.G_BAM, E.G_TAG, E.G_TIES = a.bam, 'ms', E.load_ties(a.tab, E.G_FLOOR, 'ms')
    with multiprocessing.get_context('fork').Pool(a.threads) as pool:
        res = pool.map(E.scan, list(pysam.AlignmentFile(a.bam).references), chunksize=1)
    sets = collections.defaultdict(set)
    for _, pairs, _ in res:
        for q, g, exact, *x in pairs:
            if exact and g != E.OUTSIDE: sets[q].add(g)
    parent = {}
    def find(x):
        while parent.setdefault(x, x) != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for s in sets.values():
        if len(s) > 1:
            s = sorted(s)
            for y in s[1:]: parent[find(y)] = find(s[0])
    groups = collections.defaultdict(set); gr = collections.defaultdict(set)
    for q, s in sets.items():
        if len(s) > 1: k = find(next(iter(s))); groups[k] |= s; gr[k].add(q)
    B = pysam.AlignmentFile(a.bam)
    with open(a.out + '.tmp', 'w') as w:
        w.write('group\tloci\tread_at\ttied_reads\tread_there\tsites\thaplotypes\tsizes\treads\n')
        for k in sorted(groups, key=lambda k: (meta[k][0], meta[k][1])):
            reach = collections.Counter(g for q in gr[k] for g in sets[q]); rep = max(sorted(reach), key=lambda g: reach[g])
            c, s, e, st = meta[rep][:4]; reads = {}
            for r in B.fetch(c, s - 1, e):
                q = r.query_name
                if r.is_unmapped or r.is_supplementary or q not in gr[k] or rep not in sets[q]: continue
                if r.is_secondary and r.get_tag('ms') != E.G_TIES.get(q): continue
                if E.tx_strand(r) != st and st in '+-': continue
                if q not in reads or r.query_alignment_length > reads[q][2]: reads[q] = (*variants(r), r.query_alignment_length)
            sites, haps = haplotypes({q: v[:2] for q, v in reads.items()}, a.min_reads, a.min_share)
            big = [h for h in haps if len(h) >= a.min_reads]
            w.write(f"{k}\t{','.join(sorted(groups[k]))}\t{rep}\t{len(gr[k])}\t{len(reads)}\t{len(sites)}\t{len(big)}\t{','.join(str(len(h)) for h in sorted(haps, key=len, reverse=True))}\t{';'.join(','.join(h) for h in big)}\n")
    os.replace(a.out + '.tmp', a.out)
    print(f'sample_psv: {len(groups)} groups of loci joined by tied reads -> {a.out}')

if __name__ == '__main__':
    main()
