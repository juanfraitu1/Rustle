#!/usr/bin/env python3
"""Do secondary alignments find silent relatives, and can those places be told from the others? (register 1955)

    silent_reach.py COUNTS.tsv LOCI.gtf BAM SCENARIO OUT_PREFIX [THREADS]

COUNTS.tsv (row 1909's planted set): set, unit, gene, name, role, then one column of planted reads per scenario; reads are
named SCENARIO|gene|... A PLACE = secondary alignments of the reads of one source gene whose aligned blocks overlap, off the
source locus, >= 2 reads. Class: `silent relative` (a silent member of the source's unit, >= 50 exon bases), `other planted
locus`, `elsewhere` (split by the CAT annotation when DENOVO_CAT_GFF / DENOVO_CAT_TSV are set).
Writes OUT_PREFIX.places.tsv; prints the near-tie reach, each feature's AUC and a cross-validated logistic model."""
import sys, os, re, collections, multiprocessing, math
G = {}
FEATS = ('reads', 'best_ratio', 'median_ratio', 'covered', 'identity', 'spliced', 'bases')

def ovl(bl, ex):
    n = i = j = 0
    while i < len(bl) and j < len(ex):
        lo, hi = max(bl[i][0], ex[j][0]), min(bl[i][1], ex[j][1])
        if hi > lo: n += hi - lo
        if bl[i][1] < ex[j][1]: i += 1
        else: j += 1
    return n

def scan(c):
    import pysam
    out = []
    for r in pysam.AlignmentFile(G['bam']).fetch(c):
        if r.is_unmapped or r.is_supplementary: continue
        bl = r.get_blocks(); al = sum(b - a for a, b in bl); mis = sum(n for o, n in r.cigartuples if o in (1, 2, 8))
        hit = None
        for g in G['by'].get(c, ()):
            e = G['ex'][g]
            if e[-1][1] <= bl[0][0] or e[0][0] >= bl[-1][1]: continue
            if ovl(bl, e) >= 50: hit = g; break
        out.append((r.query_name, c, r.is_secondary, r.get_tag('ms'), r.query_alignment_length / max(1, r.infer_read_length() or 1), 1 - mis / max(1, al + mis), 'N' in r.cigarstring, al, bl, hit))
    return out

def auc(pos, neg):
    """probability that a positive scores above a negative (ties count half)"""
    v = sorted([(x, 1) for x in pos] + [(x, 0) for x in neg]); r = 0.0; i = 0
    while i < len(v):
        j = i
        while j < len(v) and v[j][0] == v[i][0]: j += 1
        r += sum(k for _, k in v[i:j]) * ((i + j + 1) / 2.0); i = j
    return (r - len(pos) * (len(pos) + 1) / 2.0) / max(1, len(pos) * len(neg))

def logit_cv(X, y, groups, folds=5):
    import numpy as np
    X = np.array(X, float); y = np.array(y, float); X[:, 0] = np.log1p(X[:, 0]); X[:, 6] = np.log1p(X[:, 6]); p = np.zeros(len(y))
    keys = sorted(set(groups)); fold = {k: i % folds for i, k in enumerate(keys)}; f = np.array([fold[g] for g in groups])
    for k in range(folds):
        tr = f != k; mu = X[tr].mean(0); sd = X[tr].std(0) + 1e-9; A = np.c_[(X - mu) / sd, np.ones(len(y))]; w = np.zeros(A.shape[1])
        for _ in range(400):
            z = 1 / (1 + np.exp(-np.clip(A[tr] @ w, -30, 30))); w -= 0.5 * (A[tr].T @ (z - y[tr]) / tr.sum() + 1e-3 * w)
        p[~tr] = 1 / (1 + np.exp(-np.clip(A[~tr] @ w, -30, 30)))
    return p

def prec_at(p, y, recall):
    o = sorted(zip(p, y), reverse=True); need = recall * sum(y); tp = 0
    for i, (_, t) in enumerate(o, 1):
        tp += t
        if tp >= need: return tp / i, i
    return float('nan'), 0

def main(cnt, gtf, bam, sc, out, thr=16):
    import pysam
    rows = [l.rstrip('\n').split('\t') for l in open(cnt)]; h = rows[0]; rows = rows[1:]; ci = h.index(sc)
    info = {r[2]: r for r in rows}; ex = collections.defaultdict(list); loc = {}
    for l in open(gtf):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon':
            g = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
            if g in info: ex[g].append((int(f[3]) - 1, int(f[4]))); loc[g] = f[0]
    by = collections.defaultdict(list)
    for g in ex: ex[g] = sorted(set(ex[g])); by[loc[g]].append(g)
    G.update(bam=bam, ex=ex, by=by)
    with multiprocessing.get_context('fork').Pool(int(thr)) as pool: res = pool.map(scan, list(pysam.AlignmentFile(bam).references), chunksize=1)
    recs = [x for o in res for x in o]; best = collections.defaultdict(int)
    for x in recs: best[x[0]] = max(best[x[0]], x[3])
    # places: per (source gene, contig) merge the overlapping blocks of secondary alignments off the source locus
    per = collections.defaultdict(list)
    for x in recs:
        src = x[0].split('|')[1]
        if x[2] and x[9] != src and src in info: per[(src, x[1])].append(x)
    P = []
    for (src, c), v in per.items():
        v.sort(key=lambda x: x[8][0][0]); cur = None
        for x in v:
            s, e = x[8][0][0], x[8][-1][1]
            if cur and s < cur[1]: cur[1] = max(cur[1], e); cur[2].append(x)
            else: cur = [s, e, [x]]; P.append((src, c, cur))
    idx = None
    if os.environ.get('DENOVO_CAT_GFF') and os.environ.get('DENOVO_CAT_TSV'):
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
        import famtruth as F
        CAT = F.load_cat(); idx = F.exon_index(CAT)
    T = []; med = lambda a: sorted(a)[len(a) // 2]
    for src, c, (s, e, xs) in P:
        reads = {x[0] for x in xs}
        if len(reads) < 2: continue
        hits = collections.Counter(x[9] for x in xs if x[9]); g = hits.most_common(1)[0][0] if hits else None
        if g and info[g][1] == info[src][1] and info[g][ci] == '0': k = 'silent relative'
        elif g: k = 'other planted locus'
        else:
            k = 'elsewhere'
            if idx is not None:
                a = None
                for st in '+-':
                    a = a or F.map_copy(idx, c, st, [b for x in xs[:20] for b in x[8]])
                bt = CAT[a]['biotype'] if a else None
                k += ': ' + ('nothing annotated' if not a else 'protein-coding gene' if bt == 'protein_coding' else 'pseudogene' if 'pseudo' in bt else 'other gene')
        ratio = [x[3] / max(1, best[x[0]]) for x in xs]
        T.append(dict(src=src, unit=info[src][1], set=info[src][0], role=(info[g][4] if g else '.'), locus=g or '.', contig=c, start=s, end=e, klass=k, reads=len(reads), best_ratio=max(ratio), median_ratio=med(ratio),
                      covered=med([x[4] for x in xs]), identity=med([x[5] for x in xs]), spliced=sum(x[6] for x in xs) / len(xs), bases=med([x[7] for x in xs]), near=int(any(best[x[0]] - x[3] <= 10 for x in xs))))
    cols = ['src', 'unit', 'set', 'locus', 'role', 'contig', 'start', 'end', 'klass', 'near'] + list(FEATS)
    with open(out + '.places.tsv', 'w') as w:
        w.write('\t'.join(cols) + '\n')
        for t in T: w.write('\t'.join(f'{t[c]:.4f}' if isinstance(t[c], float) else str(t[c]) for c in cols) + '\n')
    K = collections.Counter(t['klass'] for t in T); print(f'places with >= 2 reads: {len(T)}; ' + '; '.join(f'{k} {v}' for k, v in K.most_common()))
    print('silent relatives reached (distinct loci), and by a read within 10 ms of its best:')
    by_role = collections.defaultdict(dict)
    for t in T:
        if t['klass'] == 'silent relative': d = by_role[(t['set'], t['role'])]; d[t['locus']] = max(d.get(t['locus'], 0), t['near'])
    for k, d in sorted(by_role.items()): print(f'  {k[0]} {k[1]}: reached {len(d)}, in the near-tie band {sum(d.values())} ({sum(d.values()) / max(1, len(d)):.3f})')
    pos = [t for t in T if t['klass'] == 'silent relative']; neg = [t for t in T if t['klass'].startswith('elsewhere')]
    print(f'pattern: silent relatives {len(pos)} places against {len(neg)} places elsewhere')
    print('feature\tmedian, silent relatives\tmedian, elsewhere\tAUC')
    for f in FEATS: print(f'{f}\t{med([t[f] for t in pos]):.3f}\t{med([t[f] for t in neg]):.3f}\t{auc([t[f] for t in pos], [t[f] for t in neg]):.3f}')
    both = pos + neg; y = [1] * len(pos) + [0] * len(neg); p = logit_cv([[t[f] for f in FEATS] for t in both], y, [t['unit'] for t in both])
    print(f'logistic model, 5-fold by source unit: AUC {auc([a for a, b in zip(p, y) if b], [a for a, b in zip(p, y) if not b]):.3f}')
    for r in (0.5, 0.8, 0.9):
        pr, n = prec_at(p, y, r); print(f'  at recall {r}: precision {pr:.3f} ({n} places called)')
    top = sorted(zip(p, both), key=lambda z: -z[0]); fp = collections.Counter(t['klass'] for a, t in top[:len(pos)] if not t['klass'].startswith('silent'))
    print('the highest-scoring places that are NOT listed relatives, by annotation: ' + '; '.join(f'{k} {v}' for k, v in fp.most_common()))
    for role in sorted({t['role'] for t in pos}):
        pp = [a for a, t in zip(p, both) if t['klass'] == 'silent relative' and t['role'] == role]; nn = [a for a, b in zip(p, y) if not b]
        print(f'  AUC for silent {role}s alone ({len(pp)} places): {auc(pp, nn):.3f}')

def selftest():
    assert abs(auc([3, 4], [1, 2]) - 1) < 1e-9 and abs(auc([1, 2], [1, 2]) - 0.5) < 1e-9 and ovl([(0, 10), (20, 30)], [(5, 25)]) == 10
    assert prec_at([0.9, 0.8, 0.1], [1, 0, 1], 0.5)[0] == 1.0
    print('silent_reach selftest ok')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: selftest()
    else: main(*sys.argv[1:])
