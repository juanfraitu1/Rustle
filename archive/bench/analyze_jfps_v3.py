#!/usr/bin/env python3
"""
Junction uniqueness analysis: for each j-FP junction, check how many assembled
transcripts share it. De novo proxy for "novel": junction unique to a single transcript.
"""
import sys, re, json, collections

TMAP   = "jfp_study_cmp.jfp_study_denovo.gtf.tmap"
RGTF   = "/mnt/c/Users/jfris/Desktop/Rustle/bench/jfp_study_denovo.gtf"
REFGTF = "/mnt/c/Users/jfris/Desktop/GGO_19.gtf"
PARITY = "/mnt/c/Users/jfris/Desktop/Rustle/bench/jfp_study_parity.jsonl"

def parse_attr(s, k):
    m = re.search(rf'{k}\s+"([^"]+)"', s)
    return m.group(1) if m else None

def parse_float(s, k):
    v = parse_attr(s, k); return float(v) if v else None

# ---- load reference junctions ----
ref_juncs = set()
cur_exons = []
with open(REFGTF) as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip().split('\t')
        if len(f) < 9: continue
        feat, s, e = f[2], int(f[3])-1, int(f[4])
        if feat == 'transcript':
            if len(cur_exons) > 1:
                cur_exons.sort()
                for i in range(len(cur_exons)-1):
                    ref_juncs.add((cur_exons[i][1], cur_exons[i+1][0]))
            cur_exons = []
        elif feat == 'exon': cur_exons.append((s, e))
if len(cur_exons) > 1:
    cur_exons.sort()
    for i in range(len(cur_exons)-1): ref_juncs.add((cur_exons[i][1], cur_exons[i+1][0]))

# ---- load Rustle GTF ----
tx_meta = {}; cur_tid = None; cur_exons = []
with open(RGTF) as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip().split('\t')
        if len(f) < 9: continue
        feat, s, e, attrs = f[2], int(f[3])-1, int(f[4]), f[8]
        if feat == 'transcript':
            if cur_tid:
                cur_exons.sort(); tx_meta[cur_tid]['exons'] = cur_exons
            cur_tid = parse_attr(attrs, 'transcript_id')
            tx_meta[cur_tid] = {
                'gene_id': parse_attr(attrs, 'gene_id'),
                'cov': parse_float(attrs, 'cov') or 0.0,
                'longcov': parse_float(attrs, 'longcov') or 0.0,
                'source': parse_attr(attrs, 'source') or '?',
                'strand': f[6], 'chrom': f[0], 'exons': []
            }; cur_exons = []
        elif feat == 'exon' and cur_tid: cur_exons.append((s, e))
if cur_tid:
    cur_exons.sort(); tx_meta[cur_tid]['exons'] = cur_exons

# ---- load tmap ----
tx_class = {}; tx_best_ref = {}
with open(TMAP) as fh:
    for i, line in enumerate(fh):
        if i == 0: continue
        f = line.rstrip().split('\t')
        if len(f) < 5: continue
        tx_class[f[4]] = f[2]; tx_best_ref[f[4]] = f[0]

# ---- load parity log ----
junc_mm = {}
with open(PARITY) as fh:
    for line in fh:
        try: d = json.loads(line.strip())
        except: continue
        if d.get('step') == 'junction_accept':
            junc_mm[(d['start'], d['end'])] = d['payload']['mm']

def junc_mm_lookup(d, a): return junc_mm.get((d, a+1), -1.0)

# ---- build junction usage count across ALL assembled transcripts ----
junc_tx_count = collections.Counter()  # (donor, acceptor) -> n transcripts using it
junc_gene_count = collections.Counter()  # (gene_id, donor, acceptor) -> n transcripts at gene

for tid, meta in tx_meta.items():
    exons = meta['exons']
    gene = meta['gene_id']
    for i in range(len(exons)-1):
        d = exons[i][1]; a = exons[i+1][0]
        junc_tx_count[(d, a)] += 1
        junc_gene_count[(gene, d, a)] += 1

# ---- build j-FP rows ----
def all_juncs(meta):
    ex = meta['exons']
    return [(ex[i][1], ex[i+1][0]) for i in range(len(ex)-1)]

def novel_juncs(meta):
    return [(d, a) for d, a in all_juncs(meta) if (d, a) not in ref_juncs]

gene_eq_count = collections.Counter(
    meta['gene_id'] for tid, meta in tx_meta.items() if tx_class.get(tid) == '='
)

rows = []
for tid, meta in tx_meta.items():
    if tx_class.get(tid) != 'j': continue
    exons = meta['exons']
    if len(exons) < 2: continue
    juncs = all_juncs(meta)
    novel = novel_juncs(meta)
    gene = meta['gene_id']
    # per-junction stats
    jstats = []
    for d, a in juncs:
        mm = junc_mm_lookup(d, a)
        in_ref = (d, a) in ref_juncs
        tx_count = junc_tx_count[(d, a)]  # how many assembled transcripts use this junction
        gene_tx_count = junc_gene_count[(gene, d, a)]
        jstats.append({'d': d, 'a': a, 'mm': mm, 'in_ref': in_ref, 'tx_count': tx_count, 'gene_count': gene_tx_count})
    # "singleton" junctions: used by only 1 assembled transcript at this gene
    singleton_juncs = [j for j in jstats if j['gene_count'] == 1]
    # unique-and-novel: not in ref AND singleton at this gene
    novel_singleton_juncs = [j for j in jstats if not j['in_ref'] and j['gene_count'] == 1]
    min_singleton_mm = min((j['mm'] for j in singleton_juncs if j['mm'] >= 0), default=-1.0)
    min_novel_singleton_mm = min((j['mm'] for j in novel_singleton_juncs if j['mm'] >= 0), default=-1.0)
    min_novel_mm_all = min((j['mm'] for j in jstats if not j['in_ref'] and j['mm'] >= 0), default=-1.0)
    min_mm_all = min((j['mm'] for j in jstats if j['mm'] >= 0), default=-1.0)
    rows.append({
        'tid': tid, 'gene': gene,
        'source': meta['source'],
        'cov': meta['cov'], 'longcov': meta['longcov'],
        'n_juncs': len(juncs), 'n_novel': len(novel),
        'n_singleton': len(singleton_juncs),
        'n_novel_singleton': len(novel_singleton_juncs),
        'min_mm_all': min_mm_all,
        'min_novel_mm_all': min_novel_mm_all,
        'min_singleton_mm': min_singleton_mm,
        'min_novel_singleton_mm': min_novel_singleton_mm,
        'eq_at_locus': gene_eq_count[gene],
        'jstats': jstats,
    })

print(f"Total j-FPs: {len(rows)}")
print(f"By n_novel: " + ", ".join(f"{k}:{v}" for k,v in sorted(collections.Counter(r['n_novel'] for r in rows).items())))
print(f"By n_novel_singleton: " + ", ".join(f"{k}:{v}" for k,v in sorted(collections.Counter(r['n_novel_singleton'] for r in rows).items())))

# ---- TP analogous stats ----
tp_rows = []
for tid, meta in tx_meta.items():
    if tx_class.get(tid) != '=' or len(meta['exons']) < 2: continue
    gene = meta['gene_id']
    juncs = all_juncs(meta)
    jstats = []
    for d, a in juncs:
        mm = junc_mm_lookup(d, a)
        in_ref = (d, a) in ref_juncs
        gene_tx_count = junc_gene_count[(gene, d, a)]
        jstats.append({'mm': mm, 'in_ref': in_ref, 'gene_count': gene_tx_count})
    singleton_juncs = [j for j in jstats if j['gene_count'] == 1]
    novel_singleton_juncs = [j for j in jstats if not j['in_ref'] and j['gene_count'] == 1]
    min_singleton_mm = min((j['mm'] for j in singleton_juncs if j['mm'] >= 0), default=-1.0)
    min_novel_singleton_mm = min((j['mm'] for j in novel_singleton_juncs if j['mm'] >= 0), default=-1.0)
    tp_rows.append({
        'n_novel_singleton': len(novel_singleton_juncs),
        'min_singleton_mm': min_singleton_mm,
        'min_novel_singleton_mm': min_novel_singleton_mm,
        'longcov': meta['longcov'],
    })

print(f"\nTPs: {len(tp_rows)}")
print(f"TPs with n_novel_singleton > 0: {sum(1 for r in tp_rows if r['n_novel_singleton'] > 0)}")
print(f"j-FPs with n_novel_singleton > 0: {sum(1 for r in rows if r['n_novel_singleton'] > 0)}")
print(f"j-FPs with n_novel_singleton == 0 (chimeric, all shared juncs): {sum(1 for r in rows if r['n_novel_singleton'] == 0)}")

# ---- FILTER SIMULATION using novel singleton junctions ----
print(f"\n{'='*80}")
print(f"FILTER: require novel singleton junctions (if any) to have mm >= threshold")
print(f"Apply only to transcripts that HAVE novel singleton junctions")
print(f"{'='*80}")

def pct(n, d): return f"{100*n/d:.1f}%" if d else "0%"
def p(vals, ptile):
    vals = sorted(vals)
    return vals[int(len(vals)*ptile/100)] if vals else 0

for thr in [2, 3, 5, 10]:
    # j-FPs removed: have at least 1 novel singleton junction with mm < thr
    fp_rm = sum(1 for r in rows if r['n_novel_singleton'] > 0 and 0 <= r['min_novel_singleton_mm'] < thr)
    # TPs removed: same condition
    tp_rm = sum(1 for r in tp_rows if r['n_novel_singleton'] > 0 and 0 <= r['min_novel_singleton_mm'] < thr)
    print(f"  min_novel_singleton_mm < {thr}: removes {fp_rm} j-FPs, {tp_rm} TPs  "
          f"(ratio {fp_rm/(tp_rm+1e-9):.2f}x)")

print(f"\n  Note: novel_singleton = junction not in reference AND used by only 1 transcript at this gene")

# ---- also try: any singleton junction (in ref or not) ----
print(f"\n{'='*80}")
print(f"FILTER: singleton junctions (any, not just novel) with low mm")
print(f"{'='*80}")
for thr in [2, 3, 5]:
    fp_rm = sum(1 for r in rows if r['n_singleton'] > 0 and 0 <= r['min_singleton_mm'] < thr)
    tp_rm = sum(1 for r in tp_rows if r['min_singleton_mm'] >= 0 and r['min_singleton_mm'] < thr)
    # For TPs: singleton junction = used by only 1 tx at the gene; could be real but unique
    print(f"  min_singleton_mm < {thr}: removes {fp_rm} j-FPs, {tp_rm} TPs  "
          f"(ratio {fp_rm/(tp_rm+1e-9):.2f}x)")

# ---- distribution of singleton/novel-singleton junction mm ----
print(f"\n{'='*80}")
print(f"MIN_NOVEL_SINGLETON_MM distribution: j-FPs vs TPs")
print(f"{'='*80}")
fp_nsmm = [r['min_novel_singleton_mm'] for r in rows if r['min_novel_singleton_mm'] >= 0]
tp_nsmm = [r['min_novel_singleton_mm'] for r in tp_rows if r['min_novel_singleton_mm'] >= 0]
fp_nsmm.sort(); tp_nsmm.sort()
for label, vals in [("j-FPs", fp_nsmm), ("TPs", tp_nsmm)]:
    if vals:
        print(f"  {label}: n={len(vals)}  p10={p(vals,10)}  p25={p(vals,25)}  median={p(vals,50)}  "
              f"p75={p(vals,75)}  p90={p(vals,90)}  max={max(vals)}")
    else:
        print(f"  {label}: n=0")

# ---- F1 impact ----
print(f"\n{'='*80}")
print(f"F1 IMPACT of novel_singleton_mm filter")
print(f"{'='*80}")
total_refs = 1814  # from gffcompare
for thr in [2, 3, 5]:
    fp_rm = sum(1 for r in rows if r['n_novel_singleton'] > 0 and 0 <= r['min_novel_singleton_mm'] < thr)
    tp_rm = sum(1 for r in tp_rows if r['n_novel_singleton'] > 0 and 0 <= r['min_novel_singleton_mm'] < thr)
    new_tp = 1738 - tp_rm
    new_query = 1929 - fp_rm - tp_rm
    new_pr = new_tp / new_query if new_query > 0 else 0
    new_sn = new_tp / total_refs
    new_f1 = 2*new_pr*new_sn/(new_pr+new_sn) if new_pr+new_sn > 0 else 0
    old_pr, old_sn = 0.9012, 0.9579
    old_f1 = 2*old_pr*old_sn/(old_pr+old_sn)
    print(f"  thr={thr}: Pr={100*new_pr:.2f}%  Sn={100*new_sn:.2f}%  F1={100*new_f1:.2f}%  "
          f"(delta_F1={100*(new_f1-old_f1):+.2f}pp)  removes {fp_rm} j-FPs, {tp_rm} TPs")

# ---- detailed j-FP table showing singleton junction info ----
print(f"\n{'='*80}")
print(f"j-FPs with novel singleton junctions (n_novel_singleton > 0), sorted by mm")
print(f"{'='*80}")
novel_single = [r for r in rows if r['n_novel_singleton'] > 0]
novel_single.sort(key=lambda r: r['min_novel_singleton_mm'])
print(f"  n={len(novel_single)}")
print(f"  {'tid':<16} {'src':<18} {'lc':>5} {'cov':>6} {'nv_sng':>7} {'min_ns_mm':>10} {'min_all_mm':>11}")
for r in novel_single[:40]:
    print(f"  {r['tid']:<16} {r['source']:<18} {r['longcov']:>5.1f} {r['cov']:>6.2f} "
          f"{r['n_novel_singleton']:>7} {r['min_novel_singleton_mm']:>10.0f} {r['min_mm_all']:>11.0f}")

# ---- analysis of why chimeric j-FPs can't be filtered with novelty approach ----
print(f"\n{'='*80}")
print(f"CHIMERIC j-FPs (n_novel=0): why they survive novelty filter")
print(f"{'='*80}")
chim = [r for r in rows if r['n_novel'] == 0]
print(f"  n_chimeric={len(chim)}")
print(f"  All junctions in ref: by definition for chimeric")
print(f"  Junction uniqueness at gene level:")
has_singleton = sum(1 for r in chim if r['n_singleton'] > 0)
print(f"  chimeric j-FPs with singleton junctions (within assembled set): {has_singleton}")
print(f"  chimeric j-FPs with NO singleton junctions: {len(chim) - has_singleton}")

print(f"\n  Distribution of min_singleton_mm for chimeric j-FPs with singletons:")
chim_sng = [r['min_singleton_mm'] for r in chim if r['min_singleton_mm'] >= 0]
chim_sng.sort()
if chim_sng:
    for thr in [2, 3, 5, 10]:
        cnt = sum(1 for v in chim_sng if v < thr)
        print(f"    min_singleton_mm < {thr}: {cnt} / {len(chim_sng)} chimeric j-FPs")

# TP chimeric comparison (TPs that have singleton junctions)
tp_sng = [r['min_singleton_mm'] for r in tp_rows if r['min_singleton_mm'] >= 0]
tp_sng.sort()
print(f"\n  Distribution of min_singleton_mm for TPs with singleton junctions: n={len(tp_sng)}")
for thr in [2, 3, 5, 10]:
    cnt = sum(1 for v in tp_sng if v < thr)
    print(f"    TP min_singleton_mm < {thr}: {cnt} / {len(tp_sng)}")

print(f"\n  Ratio analysis (chimeric j-FPs per TP removed):")
for thr in [2, 3, 5]:
    fp_rm = sum(1 for r in chim if 0 <= r['min_singleton_mm'] < thr)
    tp_rm = sum(1 for r in tp_rows if 0 <= r['min_singleton_mm'] < thr)
    print(f"    singleton_mm < {thr}: removes {fp_rm} chimeric j-FPs, {tp_rm} TPs  "
          f"ratio={fp_rm/(tp_rm+1e-9):.2f}x")
