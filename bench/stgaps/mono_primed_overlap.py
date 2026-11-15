#!/usr/bin/env python3
"""Row 1666 step 0: recount the single-exon priming filter's evidence under three read-overlap rules.

  mono_primed_overlap.py BAM FASTA FUNNEL.tsv CONTIG DEFAULT.gtf FILTER_OFF.gtf     (env as pairscore.py; --selftest)
"""
import sys, os, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))

def rich(w, b): return bool(w) and (w.count(b) >= 0.6 * len(w) or b * 6 in w)

def verdict(model, reads):
    """reads [(a, b, primed)] -> {rule: (counted, primed)}"""
    s, e = model; out = {k: [0, 0] for k in ('model', 'shorter', 'any')}
    for a, b, p in reads:
        o = min(e, b) - max(s, a)
        for k, need in (('model', 0.5 * (e - s)), ('shorter', 0.5 * min(e - s, b - a)), ('any', 1)):
            if o >= need: out[k][0] += 1; out[k][1] += p
    return out

dropped = lambda c: c[0] > 0 and c[1] / c[0] >= 0.8

def monos(gtf, contig):
    n = collections.Counter(); sp = {}
    for l in open(gtf):
        f = l.split('\t')
        if len(f) > 8 and f[0] == contig and f[2] == 'exon':
            t = f[8].split('transcript_id "')[1].split('"')[0]; n[t] += 1; sp[t] = (int(f[3]) - 1, int(f[4]), f[6])
    return {sp[t] for t in n if n[t] == 1}

def main(bam, fa, funnel, contig, dflt, off):
    import pysam, famtruth as F
    G = F.genes_from_argv([], 'human'); B = pysam.AlignmentFile(bam); S = pysam.FastaFile(fa).fetch(contig).upper()
    nomodel = {l.split('\t')[0] for l in open(funnel) if l.split('\t')[5] == 'no_model'}
    genes = sorted((g['exons'][0][0], g['exons'][-1][1], k) for k, g in G.items() if g['chrom'] == contig)
    def cls(s, e):
        best = 'none'
        for a, b, k in genes:
            if a >= e: break
            if b <= s: continue
            ex = G[k]['exons']
            if 2 * sum(max(0, min(y, e) - max(x, s)) for x, y in ex) >= sum(y - x for x, y in ex):
                if k in nomodel: return 'member'
                best = 'annotated'
        return best
    drops = sorted(monos(off, contig) - monos(dflt, contig)); C = collections.Counter()
    for s, e, st in drops:
        reads = [(r.reference_start, r.reference_end, rich(S[r.reference_end:r.reference_end + 20], 'A') or rich(S[max(0, r.reference_start - 20):r.reference_start], 'T'))
                 for r in B.fetch(contig, s, e) if not r.flag & 2820 and 'N' not in r.cigarstring]
        v = verdict((s, e), reads); c = cls(s, e); C['drops'] += 1; C['drops:' + c] += 1
        C['reproduced'] += dropped(v['model']); C['decided_by<=3'] += v['model'][0] <= 3
        for k in ('shorter', 'any'):
            if not dropped(v[k]): C[f'flip:{k}'] += 1; C[f'flip:{k}:{c}'] += 1
        if c == 'member': print('#member', contig, s, e, v)
    for k in sorted(C): print(contig, k, C[k])

def selftest():
    long_primed = [(0, 5000, True)]; short_clean = [(0, 1500, False)] * 50
    v = verdict((0, 5000), long_primed + short_clean)
    assert v['model'] == [1, 1] and v['shorter'] == [51, 1] and v['any'] == [51, 1]
    assert dropped(v['model']) and not dropped(v['shorter']) and not dropped([0, 0])
    assert rich('AAAAAAGCGCGCGCGCGCGC', 'A') and not rich('ACGT' * 5, 'A')
    print('selftest ok')

if __name__ == '__main__':
    selftest() if '--selftest' in sys.argv else main(*sys.argv[1:7])
