#!/usr/bin/env python3
"""One pass over the BAM for every step that works on tied reads: which loci each multimapping read can sit at.

    tie_table.py ASSEMBLY.gtf BAM AS_TABLE OUT.tsv.gz [--threads N]

The provenance labels (bench/o2_families/em_filter.py) and the sequence test (bench/o2_families/recurrence_test.py) both
need, for every read with a second alignment within 10 `ms` of its best, the loci of ASSEMBLY.gtf its alignments fall on.
Each used to make that pass itself; this writes it once, readable, and both read it (their `--ties FILE`).

OUT.tsv.gz, one row per (read, locus):
  read      the read's name
  locus     gene_id of the locus the alignment overlaps most (as tools/eq_classes.py), or `outside` (no locus there)
  gap       how far that alignment's `ms` is below the read's best: 0 = it ties the best (or is the primary), 1..10 = near
  best_ms   the read's best `ms` over all its alignments (from the score table)
OUT.tsv.gz.placed.tsv: per locus, the reads with NO second alignment in the band, i.e. placed there without question.
The first line of OUT.tsv.gz is a comment naming the three inputs with size and time; a table whose inputs changed is
not used. Only reads with a rival in the band are listed: the others are in the `placed` counts."""
import sys, os, re, gzip, json, argparse, collections, multiprocessing
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import eq_classes as E

def stamp(*paths):
    """path, size and time of each input: what a cached product is checked against"""
    return [[os.path.abspath(p), os.path.getsize(p), int(os.path.getmtime(p))] for p in paths]

def signature(gtf, bam, tab): return stamp(gtf, bam, tab) + ['ms', 10, 'v1']

def check_table(tab, bam):
    """stop when the score table was made from another BAM (its header records the BAM's size): the scores of another BAM give
    wrong ties without any error (2026-11-24: a table overwritten for a reduced BAM). Size, not path, so a moved BAM passes."""
    with open(tab) as f: m = re.search(r'\tbam_bytes=(\d+)', f.readline())
    if m and int(m.group(1)) != os.path.getsize(bam): sys.exit(f'tie_table: {tab} was made from a BAM of {m.group(1)} bytes; {bam} has {os.path.getsize(bam)}')

def scan(contig):
    """as eq_classes.scan, with every candidate in the band kept together with its ms gap: -> placed, [(read, locus, gap)]"""
    import pysam
    placed, pairs = collections.Counter(), []
    idx = {s: ([x for x, _, _ in E.G_GEX[(contig, s)]], max(b - x for x, b, _ in E.G_GEX[(contig, s)])) for s in "+-" if (contig, s) in E.G_GEX}
    for r in pysam.AlignmentFile(E.G_BAM).fetch(contig):
        if r.is_unmapped or r.is_supplementary: continue
        best = E.G_TIES.get(r.query_name)
        if best is None:
            if r.is_secondary: continue
        elif r.is_secondary and r.get_tag(E.G_TAG) < E.G_FLOOR(best): continue
        s = E.tx_strand(r)
        ov = E.overlap_genes([(x + 1, b) for x, b in r.get_blocks()], E.G_GEX[(contig, s)], *idx[s]) if s in idx else None
        g = max(ov, key=lambda k: (ov[k], k)) if ov else E.OUTSIDE
        if best is None: placed[g] += 1
        else: pairs.append((r.query_name, g, 0 if not r.is_secondary else max(0, best - r.get_tag(E.G_TAG))))
    return placed, pairs

def compute(gtf, bam, tab, threads=8):
    """-> meta (loci), placed {locus: n}, pairs [(read, locus, gap)] in BAM order, best {read: best ms}"""
    import pysam
    meta, E.G_GEX = E.load_loci(gtf)
    E.G_FLOOR = E.floor_fn(10); E.G_BAM, E.G_TAG, E.G_TIES = bam, 'ms', E.load_ties(tab, E.G_FLOOR, 'ms')
    with multiprocessing.get_context('fork').Pool(threads) as pool:
        res = pool.map(scan, list(pysam.AlignmentFile(bam).references), chunksize=1)
    placed = collections.Counter(); pairs = []
    for pl, pr in res: placed.update(pl); pairs += pr
    return meta, placed, pairs, E.G_TIES

def write(path, sig, placed, pairs, best):
    with gzip.open(path + '.tmp', 'wt', compresslevel=3) as w:
        w.write('#tie_table\t' + json.dumps(sig) + '\n'); w.write('read\tlocus\tgap\tbest_ms\n')
        for q, g, gap in pairs: w.write(f'{q}\t{g}\t{gap}\t{best[q]}\n')
    with open(path + '.placed.tsv.tmp', 'w') as w:
        w.write('locus\treads_placed_alone\n'); w.writelines(f'{g}\t{n}\n' for g, n in sorted(placed.items()))
    os.replace(path + '.placed.tsv.tmp', path + '.placed.tsv'); os.replace(path + '.tmp', path)

def read(path, sig):
    """-> placed, pairs, best when the table was made from the same inputs, else None"""
    if not (path and os.path.exists(path) and os.path.exists(path + '.placed.tsv')): return None
    placed = collections.Counter(); pairs = []; best = {}
    with gzip.open(path, 'rt') as f:
        head = f.readline().rstrip('\n').split('\t')
        if head[0] != '#tie_table' or json.loads(head[1]) != sig: return None
        f.readline()
        for l in f:
            q, g, gap, b = l.rstrip('\n').split('\t'); pairs.append((q, g, int(gap))); best[q] = int(b)
    with open(path + '.placed.tsv') as f:
        f.readline()
        for l in f:
            g, n = l.rstrip('\n').split('\t'); placed[g] = int(n)
    return placed, pairs, best

def get(gtf, bam, tab, threads=8, ties=None):
    """The shared pass: read `ties` when it matches the inputs, otherwise scan (and write it when a path is given).
    -> meta, placed, pairs, best; also leaves eq_classes' globals set as the scan would."""
    check_table(tab, bam); sig = signature(gtf, bam, tab)
    got = read(ties, sig) if ties else None
    if got is not None:
        meta, E.G_GEX = E.load_loci(gtf); E.G_FLOOR = E.floor_fn(10); E.G_BAM, E.G_TAG, E.G_TIES = bam, 'ms', got[2]
        print(f'tie_table: read from {ties}', file=sys.stderr)
        return meta, got[0], got[1], got[2]
    meta, placed, pairs, best = compute(gtf, bam, tab, threads)
    if ties:
        write(ties, sig, placed, pairs, best); print(f'tie_table: written to {ties}', file=sys.stderr)
    return meta, placed, pairs, best

def groups(pairs):
    """-> {read: loci it is exactly tied between}, {group: loci}, {group: tied reads}; a group = loci joined by tied reads"""
    sets = collections.defaultdict(set)
    for q, g, gap in pairs:
        if gap == 0 and g != E.OUTSIDE: sets[q].add(g)
    parent = {}
    def find(x):
        while parent.setdefault(x, x) != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for s in sets.values():
        if len(s) > 1:
            s = sorted(s)
            for y in s[1:]: parent[find(y)] = find(s[0])
    loci = collections.defaultdict(set); reads = collections.defaultdict(set)
    for q, s in sets.items():
        if len(s) > 1: k = find(next(iter(s))); loci[k] |= s; reads[k].add(q)
    return sets, loci, reads

def tie_groups(gtf, bam, tab, threads=8, ties=None):
    """the shared pass, then `groups`: -> meta, sets, {group: loci}, {group: tied reads}"""
    meta, _, pairs, _ = get(gtf, bam, tab, threads, ties)
    return (meta, *groups(pairs))

def selftest():
    import tempfile
    s_, lo, rd = groups([('r1', 'A', 0), ('r1', 'B', 0), ('r2', 'B', 0), ('r2', 'C', 0), ('r3', 'D', 0), ('r3', 'E', 4), ('r4', 'F', 0), ('r4', 'outside', 0)])
    assert [sorted(v) for v in lo.values()] == [['A', 'B', 'C']] and [sorted(v) for v in rd.values()] == [['r1', 'r2']], (lo, rd)   # a near tie (gap 4) and `outside` join nothing
    d = tempfile.mkdtemp(); p = os.path.join(d, 't.tsv.gz'); sig = [['a', 1, 2], 'ms', 10, 'v1']
    placed = collections.Counter({'G1': 7, 'outside': 2}); pairs = [('r1', 'G1', 0), ('r1', 'G2', 0), ('r2', 'G1', 0), ('r2', 'outside', 5)]; best = {'r1': 900, 'r2': 450}
    write(p, sig, placed, pairs, best)
    assert read(p, sig) == (placed, pairs, best)
    assert read(p, [['a', 1, 3], 'ms', 10, 'v1']) is None and read(os.path.join(d, 'none.gz'), sig) is None
    bam = os.path.join(d, 'x.bam'); open(bam, 'wb').write(b'12345'); tab = os.path.join(d, 't.tsv')
    open(tab, 'w').write('#as_table\tbam=/x.bam\tbam_bytes=5\trecords=1\n'); check_table(tab, bam)
    open(tab, 'w').write('#as_table\tbam=/x.bam\tbam_bytes=6\trecords=1\n')
    try: check_table(tab, bam); raise AssertionError('a table of another BAM was accepted')
    except SystemExit: pass
    print('tie_table selftest ok')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: selftest()
    else:
        ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
        ap.add_argument('gtf'); ap.add_argument('bam'); ap.add_argument('tab'); ap.add_argument('out'); ap.add_argument('--threads', type=int, default=8)
        a = ap.parse_args()
        meta, placed, pairs, best = get(a.gtf, a.bam, a.tab, a.threads, a.out)
        print(f'tie_table: {len(best)} reads with a rival within 10 ms, {len(pairs)} (read, locus) rows, {sum(placed.values())} reads placed alone -> {a.out}')
