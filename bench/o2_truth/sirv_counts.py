#!/usr/bin/env python3
"""Counts against an equimolar spike-in (SIRV mix E0): exact-chain reads, with and without the unspliced reads inside the
isoform's exons (register 1989, second addendum).     sirv_counts.py BAM TRUTH.gtf [CONTIG_PREFIX=SIRV]"""
import sys, os, collections, statistics, math, importlib.util
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../tools'))
import oneread, one_read_loci as O
_s = importlib.util.spec_from_file_location('us', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'unspliced_support.py')); US = importlib.util.module_from_spec(_s); _s.loader.exec_module(US)

def main(bam, truth, prefix='SIRV'):
    import pysam
    T = {t: v for t, v in oneread.models(truth).items() if v[0].startswith(prefix) and len(v[2]) > 1}
    by = collections.defaultdict(list)
    for t, (c, s, ex, _) in T.items(): by[(c, oneread.chain(ex))].append(t)
    exact = collections.Counter(); B = pysam.AlignmentFile(bam); contigs = sorted({v[0] for v in T.values()} & set(B.references)); n = collections.Counter(); uns = []
    for c in contigs:
        for r in B.fetch(c):
            if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
            n['primary reads on the spliced SIRV genes'] += 1
            if 'N' not in r.cigarstring: uns.append((c, r.reference_start, r.reference_end)); n['unspliced'] += 1; continue
            ts = by.get((c, O.read_chain(r)))
            if ts:
                for t in ts: exact[t] += 1 / len(ts)
                n['spliced, chain of an isoform exactly'] += 1
            else: n['spliced, other chain'] += 1
    exons = collections.defaultdict(lambda: collections.defaultdict(list))
    for t, (c, s, ex, _) in T.items():
        for e in ex: exons[c][e].append(t)
    US._G['maxlen'] = max(b - a for v in exons.values() for a, b in v); plus = collections.Counter(exact)
    idx = {c: sorted(v) for c, v in exons.items()}
    for c, a, b in uns:
        ex = idx[c]; ts = {t for i in US.inside([e[0] for e in ex], ex, a, b) for t in exons[c][ex[i]]}
        if not ts: continue
        n['unspliced, inside an exon of an isoform'] += 1; tot = sum(exact[t] for t in ts)
        for t in ts: plus[t] += exact[t] / tot if tot else 1 / len(ts)
    for k, v in n.items(): print(f'{k}\t{v}')
    if not n: print('no read on these contigs'); return
    def disp(C):
        v = [C[t] for t in T]; m = statistics.median(v) or 1
        return statistics.median(abs(math.log2((x + 0.5) / (m + 0.5))) for x in v), sum(x == 0 for x in v)
    genes = collections.defaultdict(list)
    for t, v in T.items(): genes[v[0]].append(t)
    def cv(C):
        g = [sum(C[t] for t in ts) / len(ts) for ts in genes.values()]; return statistics.pstdev(g) / (statistics.mean(g) or 1)
    print(f'spliced isoforms {len(T)} in {len(genes)} genes (truth: all equally abundant)')
    for lab, C in (('exact chain', exact), ('exact chain + unspliced', plus)):
        d, z = disp(C); print(f'{lab}\tmedian |log2(count / median)| {d:.3f}\tisoforms with no read {z}\tCV over genes of count per isoform {cv(C):.3f}\ttotal {sum(C.values()):.0f}')

if __name__ == '__main__': main(*sys.argv[1:])
