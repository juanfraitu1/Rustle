# P19 exploration (dev only): which tag / CIGAR feature separates the true copy from a wrong band candidate?
import sys, math, collections, pysam, json
SPLIT = sys.argv[1]; BAND = 10
FA = pysam.FastaFile('/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta')
exec(open(__file__.replace('features', 'score')).read().split("truth = ")[0].split("FA = ")[0])  # R, LX, LM
exec("def run" + open(__file__.replace('features', 'score')).read().split("def run")[1].split("truth = ")[0])  # run(), ll()
# annotated introns
intr = set(); ex = collections.defaultdict(list); par = {}
for l in open('/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff'):
    if l[0] == '#': continue
    f = l.split('\t')
    if len(f) > 8 and f[2] == 'exon':
        p = next((x[7:] for x in f[8].split(';') if x.startswith('Parent=')), None)
        ex[(f[0], p)].append((int(f[3]) - 1, int(f[4])))
for (c, p), e in ex.items():
    e.sort()
    for a, b in zip(e, e[1:]): intr.add((c, a[1], b[0]))
CANON = {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')}
comp = str.maketrans('ACGTN', 'TGCAN'); rc = lambda s: s.translate(comp)[::-1]
def feats(r):
    f = {k: r.get_tag(k) for k in ('ms', 'AS', 'NM', 'cm', 's1', 'de')}
    f['ll'] = ll(r)[0]
    js = []; rp = r.reference_start; clip = 0
    for op, L in r.cigartuples:
        if op == 3: js.append((rp, rp + L)); rp += L
        elif op in (0, 2, 7, 8): rp += L
        elif op == 4: clip += L
    f['introns'] = len(js); f['clip'] = -clip
    f['annot_junctions'] = sum((r.reference_name, a, b) in intr for a, b in js)
    nc = 0
    for a, b in js:
        d, acc = FA.fetch(r.reference_name, a, a + 2).upper(), FA.fetch(r.reference_name, b - 2, b).upper()
        if r.is_reverse: d, acc = rc(acc), rc(d)
        nc += (d, acc) not in CANON
    f['noncanon_junctions'] = -nc
    # genomic A-run just past the transcript 3' end (internal priming / processed copy)
    if r.is_reverse: g = rc(FA.fetch(r.reference_name, max(0, r.reference_start - 20), r.reference_start).upper())
    else: g = FA.fetch(r.reference_name, r.reference_end, r.reference_end + 20).upper()
    f['downstream_A'] = -g.count('A')
    f['spliced_len'] = r.query_alignment_length
    return f
truth = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open(f'{SPLIT}.truth.tsv') if not l.startswith('read\t')}
recs = collections.defaultdict(list)
for r in pysam.AlignmentFile(f'{SPLIT}.bam'):
    if not (r.is_unmapped or r.is_supplementary): recs[r.query_name].append(r)
def right(r, t): return r.reference_name == t[3] and r.reference_start < int(t[5]) and r.reference_end > int(t[4])
W = collections.defaultdict(lambda: collections.Counter())
for name, rs in recs.items():
    t = truth[name]; p = next(r for r in rs if not r.is_secondary)
    band = [r for r in rs if r.get_tag('ms') >= p.get_tag('ms') - BAND]
    if len(band) < 2: continue
    for r in band:
        if r.query_sequence is None:
            r.query_sequence = p.query_sequence if r.is_reverse == p.is_reverse else rc(p.query_sequence)
    good = [r for r in band if right(r, t)]; bad = [r for r in band if not right(r, t)]
    if not good or not bad: continue
    F = {id(r): feats(r) for r in band}
    g = max(good, key=lambda r: F[id(r)]['ms'])
    tie = any(abs(F[id(b)]['ll'] - F[id(g)]['ll']) < 1e-9 for b in bad)
    grp = 'LL-tie reads' if tie else 'LL-decidable reads'
    W[grp]['reads'] += 1; W[grp]['primary right'] += right(p, t)
    for b in bad:
        for k, v in F[id(g)].items():
            u = F[id(b)][k]; W[grp][(k, 'true better' if v > u else 'wrong better' if v < u else 'equal')] += 1
for grp, c in W.items():
    print(f'== {grp}: {c["reads"]} reads, primary right {c["primary right"]}  (pairs true copy vs each wrong band candidate)')
    for k in ('ms', 'AS', 'NM', 'll', 'cm', 's1', 'de', 'introns', 'annot_junctions', 'noncanon_junctions', 'clip', 'downstream_A', 'spliced_len'):
        tb, wb, eq = c[(k, 'true better')], c[(k, 'wrong better')], c[(k, 'equal')]
        print(f'  {k:20s} true better {tb:5d}  wrong better {wb:5d}  equal {eq:5d}')
