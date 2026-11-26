#!/usr/bin/env python3
"""Where the secondary alignments of each score band land, and whether they could support a true transcript (register 1961).

    band_census.py PLANT_DIR SCENARIO BAM AS_TABLE [K=3]

Planted reads are named SCENARIO|copy|transcript|i. For every secondary alignment that is not at its read's true copy:
band on AS and on ms (ratio to the read's best from AS_TABLE: >= 0.98, 0.9-0.98, under 0.9); landing = an expressed copy,
a silent copy (0 planted reads in the scenario), or no planted copy (>= 50 exon bases decide); anchored = >= 50 aligned
bases where >= K primary alignments with mapping quality >= 1 are stacked (pool_bam.anchored_bases); fits = its introns
(at least one) are all introns of one isoform that is expressed at the landing copy."""
import sys, os, re, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'o2_families'))
import pool_bam as PB
BANDS = ('>= 0.98', '0.9-0.98', '< 0.9')
def band(x, top): return BANDS[0] if x >= 0.98 * top else BANDS[1] if x >= 0.9 * top else BANDS[2]

def introns(blocks, slack=20):
    m = []
    for a, b in blocks:
        if m and a - m[-1][1] < slack: m[-1][1] = b
        else: m.append([a, b])
    return {(m[i][1], m[i + 1][0]) for i in range(len(m) - 1)}

def main(plant, sc, bam, tab, k=3):
    import pysam
    plan = [l.rstrip('\n').split('\t') for l in open(plant + '/plan.tsv')]; h = plan[0].index(sc); expr_tx = {r[3] for r in plan[1:] if int(r[h]) > 0}; expr = {r[1] for r in plan[1:] if int(r[h]) > 0}
    ex = collections.defaultdict(list); tx = collections.defaultdict(list); loc = {}; tg = {}
    for l in open(plant + '/loci.gtf'):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon':
            g = re.search(r'gene_id "([^"]+)"', f[8]).group(1); t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            ex[g].append((int(f[3]) - 1, int(f[4]))); tx[t].append((int(f[3]) - 1, int(f[4]))); loc[g] = f[0]; tg[t] = g
    by = collections.defaultdict(list); iso = collections.defaultdict(list)
    for g in ex: ex[g] = sorted(set(ex[g])); by[loc[g]].append(g)
    for t, e in tx.items():
        if t in expr_tx: iso[tg[t]].append(introns(sorted(e), 1))
    best = PB.load_best(tab, bam); B = pysam.AlignmentFile(bam); C = collections.Counter(); N = collections.Counter()
    for c in B.references:
        if not B.count(c): continue
        deep = PB.anchored_bases(B, c, int(k))
        for r in B.fetch(c):
            if not r.is_secondary: continue
            src = r.query_name.split('|')[1]; b = best.get(r.query_name)
            if b is None or b[0] is None or b[1] is None: continue
            bl = r.get_blocks(); land = None
            for g in by.get(c, ()):
                if sum(max(0, min(y, q) - max(x, p)) for x, y in bl for p, q in ex[g]) >= 50: land = g; break
            if land == src: continue
            where = 'no planted copy' if land is None else 'expressed copy' if land in expr else 'silent copy'
            anch = sum(int(deep[y] - deep[x]) for x, y in bl) >= 50; I = introns(bl); fits = bool(I) and any(I <= s for s in iso.get(land, ()))
            for score, key in (('AS', band(r.get_tag('AS'), b[0])), ('ms', band(r.get_tag('ms'), b[1]))):
                N[(score, key)] += 1; C[(score, key, where, anch, fits, bool(I))] += 1
    print(f'## {sc} {bam}: secondary alignments away from the true copy')
    print('score\tband\talignments\ton an expressed copy\ton a silent copy\tno planted copy\tspliced\tanchored\tfits a true isoform of the landing copy\tanchored and fits\tanchored, on a silent copy')
    for score in ('AS', 'ms'):
        for key in BANDS:
            n = N[(score, key)]; s = lambda f: sum(v for kk, v in C.items() if kk[0] == score and kk[1] == key and f(kk))
            print(f"{score}\t{key}\t{n}\t{s(lambda q: q[2] == 'expressed copy')}\t{s(lambda q: q[2] == 'silent copy')}\t{s(lambda q: q[2] == 'no planted copy')}\t{s(lambda q: q[5])}\t{s(lambda q: q[3])}\t{s(lambda q: q[4])}\t{s(lambda q: q[3] and q[4])}\t{s(lambda q: q[3] and q[2] == 'silent copy')}")

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        assert band(98, 100) == '>= 0.98' and band(90, 100) == '0.9-0.98' and band(89, 100) == '< 0.9' and introns([(0, 10), (15, 30), (100, 120)]) == {(30, 100)}; print('band_census selftest ok')
    else: main(*sys.argv[1:])
