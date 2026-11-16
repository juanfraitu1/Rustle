#!/usr/bin/env python3
"""Register 1761: the frozen r1677 chain score as a filter (driver stage `chainscore`, opt-in).

  chainscore.py fit   FEAT_A.tsv FEAT_B.tsv MODEL.json         pooled fit on the default chains of two samples + the cut TAU
  chainscore.py eval  MODEL.json NAME=FEAT.tsv ...             kept correct share and precision vs the read-count curve
  chainscore.py apply MODEL.json GTF BAM FASTA OUT_GTF [THREADS [GTF2 OUT_GTF2]]

apply: OUT_GTF = GTF without the multi-exon transcripts whose intron chain scores below TAU (features of
chain_rules_feat.py, standardised within the sample). Single-exon transcripts are kept. GTF2 (the layered PREFIX.final.gtf: same
chains, other ends) is filtered by the same chain list into OUT_GTF2. Last stdout line is the summary."""
import sys, os, re, json, collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chain_rules as R

def prob(d, keep, m):
    X, names = R.design(d, 'default_reads', keep); assert names == m['names'], 'feature set differs from the model'
    return 1 / (1 + np.exp(-np.clip(X @ np.array(m['w']), -30, 30)))

cmd = sys.argv[1]
if cmd == 'fit':
    S = [R.load(p) for p in sys.argv[2:4]]; K = [d['in_default'] == 1 for d in S]
    Xs = [R.design(d, 'default_reads', k) for d, k in zip(S, K)]; names = Xs[0][1]
    w = R.fit(np.vstack([x for x, _ in Xs]), np.concatenate([d['in_ref'][k] for d, k in zip(S, K)]))
    m = {'names': names, 'w': [float(x) for x in w]}; P = [prob(d, k, m) for d, k in zip(S, K)]
    # largest cut (grid 0.001) at which both samples keep >= 0.95 of their correct default chains
    tau = max(t for t in np.arange(0, 1, 0.001) if all(d['in_ref'][k][p >= t].sum() >= 0.95 * d['in_ref'][k].sum() for d, k, p in zip(S, K, P)))
    m['tau'] = round(float(tau), 3); m['fit'] = [os.path.basename(p) for p in sys.argv[2:4]]
    json.dump(m, open(sys.argv[4], 'w'), indent=1); print(f"TAU {m['tau']}")
elif cmd == 'eval':
    m = json.load(open(sys.argv[2])); print('sample\tdefault chains\tcorrect\tprecision\tkept\tcorrect kept (share)\tprecision kept\tread curve there\tdifference')
    for a in sys.argv[3:]:
        n, p = a.split('='); d = R.load(p); k = d['in_default'] == 1; y = d['in_ref'][k]; s = prob(d, k, m) >= m['tau']; c = y[s].sum()
        ref = R.prec_at(d['default_reads'][k], y, c)
        print(f'{n}\t{k.sum()}\t{y.sum():.0f}\t{y.mean():.3f}\t{s.sum()}\t{c:.0f} ({c / y.sum():.3f})\t{c / s.sum():.3f}\t{ref:.3f}\t{c / s.sum() - ref:+.3f}')
else:
    import multiprocessing, chain_rules_feat as F
    m = json.load(open(sys.argv[2])); gtf, bam, fa, out = sys.argv[3:7]; thr = int(sys.argv[7]) if len(sys.argv) > 7 else 16
    D = F.load(gtf); per = collections.defaultdict(dict)
    for (ctg, c), (r, s, g, e) in D.items(): per[ctg][c] = (r, s, g, e, 0, r)
    with multiprocessing.Pool(thr) as pool:
        rows = [r for rs in pool.imap_unordered(F.work, [(bam, fa, ctg, cand) for ctg, cand in sorted(per.items(), key=lambda x: -len(x[1]))]) for r in rs]
    d = {k: np.array([r[i] for r in rows]) if k in ('chrom', 'strand', 'key') else np.array([float(r[i]) for r in rows]) for i, k in enumerate(F.COLS)}
    p = prob(d, np.ones(len(rows), bool), m); bad = {(c, k) for c, k, x in zip(d['chrom'], d['key'], p) if x < m['tau']}
    def strip(src, dst):
        ex = collections.defaultdict(list)
        for l in open(src):
            f = l.split('\t')
            if len(f) > 8 and f[2] == 'exon': ex[(f[0], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]) - 1, int(f[4])))
        drop = set()
        for (ctg, t), e in ex.items():
            e.sort()
            if len(e) > 1 and (ctg, ','.join(f'{e[i][1]}-{e[i + 1][0]}' for i in range(len(e) - 1))) in bad: drop.add((ctg, t))
        with open(dst, 'w') as w:
            for l in open(src):
                f = l.split('\t'); t = re.search(r'transcript_id "([^"]+)"', f[8]) if len(f) > 8 else None
                if not (t and (f[0], t.group(1)) in drop): w.write(l)
        return drop
    drop = strip(gtf, out)
    if len(sys.argv) > 9: print(f'{sys.argv[8]}: dropped {len(strip(sys.argv[8], sys.argv[9]))} multi-exon transcripts by the same chain list')
    print(f'dropped {len(drop)} multi-exon transcripts ({len(bad)} of {len(rows)} chains) below score {m["tau"]}')
