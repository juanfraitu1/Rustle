#!/usr/bin/env python3
"""Register 1806: every one-read chain the default lacks, ranked by a score from the reads (driver stage `oneread`).

  oneread_tier.py fit   MODEL.json FEAT.tsv ...                                  pooled fit on the one-read chains of the row 1806 tables
  oneread_tier.py apply MODEL.json DEFAULT.gtf FLOOR1.gtf BAM FASTA OUT.gtf [THREADS]
  oneread_tier.py promote RANKED.gtf STRICT.gtf OUT.gtf PARTNER.bam ...          register 1808 (opt-in)

promote: OUT.gtf = the models of RANKED.gtf (apply's output) in the best 20% of the ranking whose exact intron chain is in at
least one primary alignment of a PARTNER library of the same species (chains compared without strand), tagged oneread_promoted "1".
Chains of STRICT.gtf (PREFIX.oneread.gtf) are left out, so both files can be appended to one GTF. Rows 1807 / 1808: these chains
match the annotation more often than the default's two-read chains; they do not reach the default's overall precision.

apply: OUT.gtf = the multi-exon models of FLOOR1 (the same assembly at a 1-read floor) with reads "1" whose intron chain DEFAULT
lacks, best score first, transcript lines tagged oneread_score and oneread_rank. Features = chain_rules_feat.py over all chains of
FLOOR1, design as chain score r1769, standardised over the sample's one-read chains. A RANKING, NOT A FILTER: row 1806 found no
cut at which these chains are as precise as the chains with two reads (AUC 0.59-0.81 across datasets). Last stdout line = summary."""
import sys, os, re, json, collections
import numpy as np
B = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, B)
import chain_rules as R
src = open(f'{B}/chainscore.py').read(); exec(src[:src.index('\ncmd = ')])  # intron_extras (chainscore.py is a script)

def design(d, one, noncanon, min_intron):
    X0, names = R.design(d, 'default_reads', one)
    return R.design_decile(X0, noncanon, min_intron), names

TOP = 0.2  # share of the ranking that may be promoted (rows 1807 / 1808)
def promote(ranked, strict, out, bams):
    import pysam, chain_rules_feat as F
    skip = set(F.load(strict)); T = collections.OrderedDict()  # tid -> [contig, exons, lines, rank]
    for l in open(ranked):
        f = l.split('\t'); t = T.setdefault(re.search(r'transcript_id "([^"]+)"', f[8]).group(1), [f[0], [], [], 0])
        t[2].append(l)
        if f[2] == 'exon': t[1].append((int(f[3]) - 1, int(f[4])))
        else: t[3] = int(re.search(r'oneread_rank "(\d+)"', f[8]).group(1))
    cand = {}
    for tid, (c, e, lines, rank) in T.items():
        e.sort(); k = (c, tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1)))
        if rank <= TOP * len(T) and k not in skip: cand[k] = tid
    first = {(c, ch[0]) for c, ch in cand}; seen = set()
    for b in bams:
        for r in pysam.AlignmentFile(b, threads=8):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or 'N' not in r.cigarstring: continue
            p = r.reference_start; ch = []
            for op, n in r.cigartuples:
                if op == 3: ch.append((p, p + n))
                if op in (0, 2, 3, 7, 8): p += n
            if (r.reference_name, ch[0]) in first and (r.reference_name, tuple(ch)) in cand: seen.add((r.reference_name, tuple(ch)))
    with open(out, 'w') as w:
        for k, tid in cand.items():
            if k in seen:
                for l in T[tid][2]: w.write(l.rstrip() + ' oneread_promoted "1";\n' if l.split('\t')[2] == 'transcript' else l)
    print(f'{len(seen)} one-read chains promoted of {len(cand)} in the best {TOP:.0%} of the ranking ({len(bams)} partner libraries)')

cmd = sys.argv[1]
if cmd == 'promote': promote(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:])
elif cmd == 'fit':
    Xs, ys = [], []
    for p in sys.argv[3:]:
        d = R.load(p); one = d['reads'] == 1
        X, names = design(d, one, d['noncanon'][one], d['min_intron'][one]); Xs.append(X); ys.append(d['y'][one])
    w = R.fit(np.vstack(Xs), np.concatenate(ys))
    json.dump({'names': names, 'w': [float(x) for x in w], 'fit': [os.path.basename(p) for p in sys.argv[3:]]}, open(sys.argv[2], 'w'), indent=1)
    print(f'{sum(map(len, ys))} one-read chains, {int(sum(y.sum() for y in ys))} correct, {len(w)} weights')
else:
    import multiprocessing, chain_rules_feat as F
    m = json.load(open(sys.argv[2])); default, floor1, bam, fa, out = sys.argv[3:8]; thr = int(sys.argv[8]) if len(sys.argv) > 8 else 16
    have = set(F.load(default)); D = F.load(floor1); per = collections.defaultdict(dict)
    for (ctg, c), (r, s, g, e) in D.items(): per[ctg][c] = (r, s, g, e, 0, r)
    with multiprocessing.Pool(thr) as pool:
        rows = [r for rs in pool.imap(F.work, [(bam, fa, ctg, cand) for ctg, cand in sorted(per.items(), key=lambda x: -len(x[1]))]) for r in rs]
    d = {k: np.array([r[i] for r in rows]) if k in ('chrom', 'strand', 'key') else np.array([float(r[i]) for r in rows]) for i, k in enumerate(F.COLS)}
    one = d['reads'] == 1
    if not one.any(): open(out, 'w').close(); print('0 one-read chains'); sys.exit(0)
    X, names = design(d, one, *intron_extras(d, one, fa)); assert names == m['names'] and X.shape[1] == len(m['w']), 'feature set differs from the model'
    p = 1 / (1 + np.exp(-np.clip(X @ np.array(m['w']), -30, 30)))
    key = lambda c, k: (c, tuple(tuple(map(int, iv.split('-'))) for iv in k.split(',')))
    score = {key(c, k): s for c, k, s in zip(d['chrom'][one], d['key'][one], p) if key(c, k) not in have}
    ex = collections.defaultdict(list); lines = collections.defaultdict(list)
    for l in open(floor1):
        f = l.split('\t')
        if len(f) < 9 or f[2] not in ('transcript', 'exon'): continue
        t = (f[0], re.search(r'transcript_id "([^"]+)"', f[8]).group(1)); lines[t].append(l)
        if f[2] == 'exon': ex[t].append((int(f[3]) - 1, int(f[4])))
    best = {}  # one model per chain: the first in file order with reads "1"
    for t, e in ex.items():
        e.sort(); k = (t[0], tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1)))
        if k in score and k not in best and ' reads "1";' in lines[t][0]: best[k] = t
    with open(out, 'w') as w:
        for rank, (k, t) in enumerate(sorted(best.items(), key=lambda x: -score[x[0]]), 1):
            for l in lines[t]:
                l = re.sub(r'(transcript_id "[^"]+)', r'\1.oneread', l, 1)  # the 1-read run reuses ids of default transcripts
                if l.split('\t')[2] == 'transcript': l = l.rstrip().rstrip(';') + f'; oneread_score "{score[k]:.4f}"; oneread_rank "{rank}";\n'
                w.write(l)
    print(f'{len(best)} one-read chains the default lacks, ranked (model {os.path.basename(sys.argv[2])})')
