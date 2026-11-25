#!/usr/bin/env python3
"""Are the family copies that only a wider seeding pool models expressed? (register 1959)

    gained_copies.py human|gorilla BAM TRUTH_READS.tsv|- TODAY.gtf ARM=GTF ...

A copy is modelled when a spliced transcript maps to it (as bench/denovo_eval/toolcompare.py). For the copies modelled
today, and for those each arm adds: primary alignments with mapping quality >= 1 on >= 50 exon bases of the copy (reads
placed there, not by a tie-break), and from neutral_truth.py's table (a capped sample of family reads) the reads whose
ORIGIN is the copy against the reads that only have it as another candidate. Env as famtruth."""
import sys, os, re, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))

def ov(bl, ex):
    return sum(max(0, min(b, y) - max(a, x)) for a, b in bl for x, y in ex)

def bucket(n): return '0' if n == 0 else '1-4' if n < 5 else '5 or more'

def main(species, bam, truth, today, *arms):
    import pysam, famtruth as F
    G = F.load_ggo(F.load_cat()) if species == 'gorilla' else F.load_cat(); idx = F.exon_index(G)
    pats = dict(F.GGO_MAINSTAYS) if species == 'gorilla' else {**F.MAINSTAYS, **F.YAGS}; members = set()
    for pat in pats.values():
        rx = re.compile(pat)
        for g, d in G.items():
            if (rx.search(d['name']) or rx.search(d.get('desc', ''))) if species == 'gorilla' else (rx.match(d['name']) and d['biotype'] == 'protein_coding'): members.add(g)
    def modelled(path):
        tx = collections.defaultdict(list); meta = {}
        for l in open(path):
            f = l.split('\t')
            if len(f) > 8 and f[2] == 'exon':
                t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); tx[t].append((int(f[3]) - 1, int(f[4]))); meta[t] = (f[0], f[6] if f[6] in '+-' else '+')
        return {g for t, ex in tx.items() if len(ex) > 1 for g in [F.map_copy(idx, meta[t][0], meta[t][1], sorted(ex))] if g in members}
    base = modelled(today); B = pysam.AlignmentFile(bam); orig = collections.Counter(); cand = collections.Counter()
    span = {g: (G[g]['chrom'], min(a for a, _ in G[g]['exons']), max(b for _, b in G[g]['exons'])) for g in members}
    if truth != '-':
        for l in list(open(truth))[1:]:
            f = l.rstrip('\n').split('\t'); o = set(f[5].split(',')) if f[5] else set()
            for x in f[6].split(';'):
                k, c, se = x.split(':')[:3]; s, e = map(int, se.split('-'))
                for g, (gc, gs, ge) in span.items():
                    if gc == c and s < ge and e > gs: (orig if k in o else cand)[g] += 1
    def placed(g):
        c, s, e = span[g]; ex = G[g]['exons']
        return sum(1 for r in B.fetch(c, s, e) if not (r.is_secondary or r.is_supplementary) and r.mapping_quality >= 1 and ov(r.get_blocks(), ex) >= 50)
    def line(lab, gs):
        gs = sorted(gs); p = [placed(g) for g in gs]; pc = collections.Counter(bucket(x) for x in p)
        t = f"\twith a truth-origin read {sum(orig[g] > 0 for g in gs)}, only as another candidate {sum(orig[g] == 0 and cand[g] > 0 for g in gs)}, no sampled read {sum(orig[g] == 0 and cand[g] == 0 for g in gs)}" if truth != '-' else ''
        print(f"{lab}\t{len(gs)} copies\treads placed there without a tie-break: 0 -> {pc['0']}, 1-4 -> {pc['1-4']}, 5 or more -> {pc['5 or more']}{t}")
        return [(G[g]['name'], x, orig[g], cand[g]) for g, x in zip(gs, p)]
    line('modelled today', base)
    for a in arms:
        lab, path = a.split('=', 1); m = modelled(path); d = line(f'added by {lab}', m - base)
        print('   ' + '; '.join(f'{n}: placed {x}, origin {o}, candidate {c}' for n, x, o, c in d)[:900])
        if base - m: print(f'   lost by {lab}: {len(base - m)}')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        assert ov([(0, 10), (20, 30)], [(5, 25)]) == 10 and bucket(0) == '0' and bucket(4) == '1-4' and bucket(5) == '5 or more'; print('gained_copies selftest ok')
    else: main(*sys.argv[1:])
