#!/usr/bin/env python3
"""The 10-AS filter against the truth, cell by cell, and what a secondary alignment is worth (register 1954).

    as10_truth.py real READS.tsv AS_TABLE S1_TOTAL S2_TOTAL [LABEL]    READS.tsv = bench/o2_truth/neutral_truth.py's table (row 1785)
    as10_truth.py planted PLANT_DIR SCENARIO BAM [THREADS]             reads named SCENARIO|copy|transcript|i

Cells as as10_census.py: the rule fires when the primary's AS is >= 10 above every other alignment; near tie on ms = another
alignment within 10 ms of the read's best, or the primary below the best ms.
real: loci and their AS come from READS.tsv (locus 0 = best AS; `primary=k` names the primary's locus; origin_loci = the true
ones), the ms cell from AS_TABLE; counts are weighted to the sample (strata S1 clear lead / S2 near tie on AS were capped).
planted: every alignment is given the planted copy it shares >= 50 exon bases with; a copy with the same sequence as the
true one (same unit, PLANT_DIR/units.tsv) is counted apart."""
import sys, os, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'o2_families'))
DB = ('above the primary', '0', '1-9', '10-49', '>= 50')
def dbin(d): return DB[0] if d < 0 else DB[1] if d == 0 else DB[2] if d < 10 else DB[3] if d < 50 else DB[4]
def cell(lead, near): return ('fires' if lead >= 10 else 'abstains', 'near' if near else 'clear')
CELLS = [(a, b) for a in ('fires', 'abstains') for b in ('clear', 'near')]

def report(C, S, lab, kinds):
    print(f'## {lab}')
    print('cell\treads\t' + '\t'.join(kinds) + '\tprimary at the origin, of the judged')
    for c in CELLS:
        n = sum(C[(c, k)] for k in kinds); j = C[(c, kinds[0])] + C[(c, kinds[1])] + (C[(c, kinds[2])] if kinds[2].startswith('primary') else 0)
        good = C[(c, kinds[0])] + (C[(c, kinds[1])] if kinds[1].startswith('primary at an identical') else 0)
        print(f'rule {c[0]}, ms {c[1]}\t{n:.0f}\t' + '\t'.join(f'{C[(c, k)]:.0f}' for k in kinds) + f'\t{good / j:.4f}' if j else f'rule {c[0]}, ms {c[1]}\t{n:.0f}')
    print('where the origin is when the primary is NOT at it\t' + '\t'.join(f'{k}: {v:.0f}' for k, v in sorted(C.items(), key=lambda kv: str(kv[0])) if isinstance(k, tuple) and k[0] == 'origin_at'))
    print('secondary alignments by their AS distance below the primary\tread near-tied on ms?\tsecondaries\tat the origin\tshare')
    for b in DB:
        for t in ('near', 'clear'):
            n, y = S[(b, t, 'n')], S[(b, t, 'y')]
            if n: print(f'{b}\t{t}\t{n:.0f}\t{y:.0f}\t{y / n:.4f}')

def real(tsv, table, n1, n2, lab):
    rows = []; took = collections.Counter()
    for l in open(tsv):
        f = l.rstrip('\n').split('\t')
        if f[0] == 'read' or f[1] == 'S0': continue
        took[f[1]] += 1; rows.append(f)
    w = {'S1': int(n1) / max(1, took['S1']), 'S2': int(n2) / max(1, took['S2'])}; want = {f[0] for f in rows}; ms = {}
    for l in open(table):
        if l[0] == '#': continue
        q = l.split('\t', 1)[0]
        if q in want:
            f = l.rstrip('\n').split('\t'); ms[q] = (int(f[6]), int(f[7]), int(f[8]))
    C = collections.Counter(); S = collections.Counter(); miss = 0
    kinds = ('primary at the origin', 'primary NOT at the origin', 'several copies fit equally (truth: tie)', 'truth cannot judge')
    for f in rows:
        if f[0] not in ms: miss += 1; continue
        AS = [int(x.rsplit(':AS', 1)[1]) for x in f[6].split(';')]; p = int(f[4].split(';')[0].split('=')[1]); org = {int(i) for i in f[5].split(',') if i}
        if len(AS) < 2 or p >= len(AS): continue
        lead = AS[p] - max(a for i, a in enumerate(AS) if i != p); b, s, pm = ms[f[0]]; near = s >= b - 10 or pm < b; c = cell(lead, near); x = w[f[1]]
        k = kinds[0] if (f[3] == 'decidable' and p in org) else kinds[1] if f[3] == 'decidable' else kinds[2] if f[3] == 'tie' else kinds[3]
        C[(c, k)] += x
        if f[3] == 'decidable':
            if p not in org:
                d = min(AS[p] - AS[i] for i in org) if org else None
                C[('origin_at', 'a locus ' + dbin(d) + ' AS' if d is not None else 'no listed locus')] += x
            for i, a in enumerate(AS):
                if i != p: S[(dbin(AS[p] - a), c[1], 'n')] += x; S[(dbin(AS[p] - a), c[1], 'y')] += x * (i in org)
    report(C, S, f'{lab}: {len(rows)} sampled reads with several loci (weights S1 {w["S1"]:.2f}, S2 {w["S2"]:.2f}; {miss} not in the score table)', kinds)

G = {}
def scan(contig):
    import pysam, score as SC
    out = []
    for r in pysam.AlignmentFile(G['bam']).fetch(contig):
        if r.is_unmapped or r.is_supplementary: continue
        bl = r.get_blocks(); best, bg = 0, None; rev = r.is_reverse
        if r.has_tag('ts') and r.get_tag('ts') == '-': rev = not rev
        st = '-' if rev else '+'
        for g in G['by'].get(contig, ()):
            if G['cloc'][g][1] != st: continue
            ce = G['cex'][g]
            if ce[-1][1] <= bl[0][0] or ce[0][0] >= bl[-1][1]: continue
            n = SC.ov(sorted(bl), ce)
            if n >= 50 and n > best: best, bg = n, g
        out.append((r.query_name, r.is_secondary, r.get_tag('AS'), r.get_tag('ms'), bg))
    return out

def planted(d, sc, bam, thr=8):
    import pysam, multiprocessing, score as SC
    truth, pl, fam = SC.load_plan(d, sc); cex = collections.defaultdict(list); cloc = {}
    for t, (c, st, ex, g) in truth.items(): cex[g] += ex; cloc[g] = (c, st)
    by = collections.defaultdict(list)
    for g in cex: cex[g] = SC.merge(cex[g]); by[cloc[g][0]].append(g)
    U = {l.split('\t')[1]: l.split('\t')[3] for l in list(open(f'{d}/units.tsv'))[1:]}
    G.update(bam=bam, by=by, cex=cex, cloc=cloc)
    with multiprocessing.get_context('fork').Pool(int(thr)) as pool: res = pool.map(scan, list(pysam.AlignmentFile(bam).references), chunksize=1)
    R = collections.defaultdict(list)
    for o in res:
        for x in o: R[x[0]].append(x)
    C = collections.Counter(); S = collections.Counter(); single = 0
    kinds = ('primary at the true copy', 'primary at an identical copy (same unit)', 'primary at another copy', 'primary on no planted copy')
    for q, v in R.items():
        pr = [x for x in v if not x[1]]; se = [x for x in v if x[1]]
        if not pr: continue
        if not se: single += 1; continue
        p = pr[0]; true = q.split('|')[1]; lead = p[2] - max(x[2] for x in se); bm = max(x[3] for x in v); near = max(x[3] for x in se) >= bm - 10 or p[3] < bm; c = cell(lead, near)
        k = kinds[0] if p[4] == true else kinds[1] if p[4] is not None and U.get(p[4]) == U.get(true) else kinds[2] if p[4] is not None else kinds[3]
        C[(c, k)] += 1
        if p[4] != true:
            ds = [p[2] - x[2] for x in se if x[4] == true]
            C[('origin_at', ('a secondary ' + dbin(min(ds)) + ' AS') if ds else 'no alignment at the true copy')] += 1
        for x in se: S[(dbin(p[2] - x[2]), 'near' if x[3] >= bm - 10 else 'clear', 'n')] += 1; S[(dbin(p[2] - x[2]), 'near' if x[3] >= bm - 10 else 'clear', 'y')] += x[4] == true
    report(C, S, f'planted {d} {sc}: {len(R)} reads, {single} with one alignment (here the ms column of the secondaries table is per ALIGNMENT)', kinds)

def selftest():
    assert dbin(-1) == DB[0] and dbin(9) == '1-9' and dbin(10) == '10-49' and cell(10, True) == ('fires', 'near') and cell(9, False) == ('abstains', 'clear')
    print('as10_truth selftest ok')

if __name__ == '__main__':
    a = sys.argv[1:]
    if a == ['--selftest']: selftest()
    elif a and a[0] == 'real': real(a[1], a[2], a[3], a[4], a[5] if len(a) > 5 else a[1])
    elif a and a[0] == 'planted': planted(*a[1:])
    else: sys.exit(__doc__)
