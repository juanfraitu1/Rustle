#!/usr/bin/env python3
"""
Comprehensive j-FP analysis.
Extracts all 119 class-'j' Rustle transcripts and annotates with:
  - source (flow / checktrf_rescue)
  - cov, longcov, n_exons
  - n_novel junctions (not in any reference transcript)
  - novel junction coordinates + approximate read support from GTF exon cov
  - locus (gene) and how many TPs are at the same locus
Prints a sorted table and summary statistics.
"""
import sys, re, collections

TMAP  = "jfp_study_cmp.jfp_study_denovo.gtf.tmap"
RGTF  = "/mnt/c/Users/jfris/Desktop/Rustle/bench/jfp_study_denovo.gtf"
REFGTF = "/mnt/c/Users/jfris/Desktop/GGO_19.gtf"

# ---------- helpers ----------
def parse_attr(attr_str, key):
    m = re.search(rf'{key}\s+"([^"]+)"', attr_str)
    return m.group(1) if m else None

def parse_attr_float(attr_str, key):
    v = parse_attr(attr_str, key)
    return float(v) if v else None

# ---------- load reference junctions ----------
ref_juncs = set()           # (chrom, donor, acceptor)  0-based half-open [donor, acceptor)
ref_tx_juncs = {}           # ref_tid -> list of (donor, acceptor)
cur_ref_tid = None
cur_ref_exons = []
print("Loading reference GTF ...", file=sys.stderr)
with open(REFGTF) as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip().split('\t')
        if len(f) < 9: continue
        chrom, feat, start, end, attrs = f[0], f[2], int(f[3])-1, int(f[4]), f[8]
        if feat == 'transcript':
            if cur_ref_tid and len(cur_ref_exons) > 1:
                cur_ref_exons.sort()
                for i in range(len(cur_ref_exons)-1):
                    d = cur_ref_exons[i][1]   # end of exon i = donor (0-based exclusive)
                    a = cur_ref_exons[i+1][0] # start of exon i+1 = acceptor (0-based inclusive)
                    ref_juncs.add((cur_ref_chrom, d, a))
                    ref_tx_juncs.setdefault(cur_ref_tid, []).append((d, a))
            cur_ref_tid = parse_attr(attrs, 'transcript_id')
            cur_ref_chrom = chrom
            cur_ref_exons = []
        elif feat == 'exon':
            cur_ref_exons.append((start, end))
# flush last
if cur_ref_tid and len(cur_ref_exons) > 1:
    cur_ref_exons.sort()
    for i in range(len(cur_ref_exons)-1):
        d = cur_ref_exons[i][1]
        a = cur_ref_exons[i+1][0]
        ref_juncs.add((cur_ref_chrom, d, a))
        ref_tx_juncs.setdefault(cur_ref_tid, []).append((d, a))
print(f"  {len(ref_juncs)} reference junctions across {len(ref_tx_juncs)} transcripts", file=sys.stderr)

# ---------- load Rustle GTF ----------
print("Loading Rustle GTF ...", file=sys.stderr)
tx_meta = {}   # tid -> {gene_id, cov, longcov, source, chrom, strand, exons: [(s,e)], exon_covs}
cur_tid = None
cur_exons = []
cur_exon_covs = []
with open(RGTF) as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip().split('\t')
        if len(f) < 9: continue
        chrom, feat, start, end, attrs = f[0], f[2], int(f[3])-1, int(f[4]), f[8]
        if feat == 'transcript':
            if cur_tid:
                cur_exons.sort()
                tx_meta[cur_tid]['exons'] = cur_exons
                tx_meta[cur_tid]['exon_covs'] = cur_exon_covs
            cur_tid = parse_attr(attrs, 'transcript_id')
            tx_meta[cur_tid] = {
                'gene_id': parse_attr(attrs, 'gene_id'),
                'cov': parse_attr_float(attrs, 'cov') or 0.0,
                'longcov': parse_attr_float(attrs, 'longcov') or 0.0,
                'source': parse_attr(attrs, 'source') or '?',
                'chrom': chrom,
                'strand': f[6],
                'exons': [],
                'exon_covs': [],
            }
            cur_exons = []
            cur_exon_covs = []
        elif feat == 'exon' and cur_tid:
            cur_exons.append((start, end))
            ec = parse_attr_float(attrs, 'cov')
            cur_exon_covs.append(ec or 0.0)
if cur_tid:
    cur_exons.sort()
    tx_meta[cur_tid]['exons'] = cur_exons
    tx_meta[cur_tid]['exon_covs'] = cur_exon_covs
print(f"  {len(tx_meta)} Rustle transcripts", file=sys.stderr)

# ---------- load tmap — get class codes ----------
print("Loading tmap ...", file=sys.stderr)
tx_class = {}   # tid -> class_code
tx_best_ref = {}  # tid -> ref_tid
with open(TMAP) as fh:
    for i, line in enumerate(fh):
        if i == 0: continue
        f = line.rstrip().split('\t')
        if len(f) < 5: continue
        ref_tid, ref_gene, cls, qgene, qtid = f[0], f[1][:20], f[2], f[3], f[4]
        tx_class[qtid] = cls
        tx_best_ref[qtid] = ref_tid

# ---------- compute per-tx junction stats ----------
def tx_juncs(meta):
    exons = meta['exons']
    juncs = []
    for i in range(len(exons)-1):
        d = exons[i][1]
        a = exons[i+1][0]
        juncs.append((meta['chrom'], d, a))
    return juncs

# gene -> list of tids
gene_tids = collections.defaultdict(list)
gene_eq_count = collections.defaultdict(int)
for tid, meta in tx_meta.items():
    gene_tids[meta['gene_id']].append(tid)
    if tx_class.get(tid) == '=':
        gene_eq_count[meta['gene_id']] += 1

# ---------- build j-FP rows ----------
rows = []
for tid, meta in tx_meta.items():
    if tx_class.get(tid) != 'j':
        continue
    juncs = tx_juncs(meta)
    novel = [(c, d, a) for c, d, a in juncs if (c, d, a) not in ref_juncs]
    n_total = len(juncs)
    n_novel = len(novel)
    gene = meta['gene_id']
    n_eq = gene_eq_count[gene]
    n_gene = len(gene_tids[gene])
    rows.append({
        'tid': tid,
        'gene': gene,
        'source': meta['source'],
        'cov': meta['cov'],
        'longcov': meta['longcov'],
        'n_exons': len(meta['exons']),
        'n_juncs': n_total,
        'n_novel': n_novel,
        'n_eq_at_locus': n_eq,
        'n_total_at_locus': n_gene,
        'chrom': meta['chrom'],
        'strand': meta['strand'],
        'novel_juncs': novel,
        'best_ref': tx_best_ref.get(tid, '?'),
    })

print(f"\n{'='*80}")
print(f"TOTAL j-FPs: {len(rows)}")
print(f"{'='*80}\n")

# ---------- summary breakdowns ----------
from collections import Counter
src_counts = Counter(r['source'] for r in rows)
novel_counts = Counter(r['n_novel'] for r in rows)
lc1_count = sum(1 for r in rows if r['longcov'] <= 1.0)
lc2_count = sum(1 for r in rows if r['longcov'] <= 2.0)

print("By source:")
for k, v in sorted(src_counts.items(), key=lambda x: -x[1]):
    print(f"  {k:<20} {v}")

print("\nBy n_novel junctions:")
for k, v in sorted(novel_counts.items()):
    print(f"  n_novel={k}  {v}")

print(f"\nLongcov distribution:")
print(f"  longcov <= 1.0: {lc1_count}")
print(f"  longcov <= 2.0: {lc2_count}")
print(f"  longcov >  2.0: {sum(1 for r in rows if r['longcov'] > 2.0)}")

lc_vals = [r['longcov'] for r in rows]
lc_vals.sort()
print(f"  median: {lc_vals[len(lc_vals)//2]:.2f}  max: {lc_vals[-1]:.2f}")

print("\nBy source x n_novel:")
counts = Counter((r['source'], r['n_novel']) for r in rows)
print(f"  {'source':<22} {'n_novel':>7} {'count':>6}")
for (src, nn), cnt in sorted(counts.items(), key=lambda x: (-x[1], x[0][0])):
    print(f"  {src:<22} {nn:>7} {cnt:>6}")

# ---------- locus analysis ----------
locus_counts = Counter(r['gene'] for r in rows)
print(f"\nTop loci by j-FP count:")
print(f"  {'gene':<14} {'jFPs':>5} {'TPs_at_locus':>13} {'total_at_locus':>15}")
for gene, cnt in locus_counts.most_common(15):
    n_eq = gene_eq_count[gene]
    n_tot = len(gene_tids[gene])
    print(f"  {gene:<14} {cnt:>5} {n_eq:>13} {n_tot:>15}")

# ---------- novel junction analysis ----------
print(f"\nAll novel junctions in j-FPs (n_novel >= 1):")
novel_junc_count = Counter()
for r in rows:
    for j in r['novel_juncs']:
        novel_junc_count[j] += 1

print(f"  Total unique novel junctions: {len(novel_junc_count)}")
print(f"  Shared by 2+ j-FPs: {sum(1 for v in novel_junc_count.values() if v >= 2)}")
print(f"\n  Most shared novel junctions:")
for junc, cnt in novel_junc_count.most_common(10):
    chrom, d, a = junc
    print(f"    {chrom}:{d}-{a}  shared by {cnt} j-FPs")

# ---------- detailed table for chimeric j-FPs (n_novel=0) ----------
chimeric = [r for r in rows if r['n_novel'] == 0]
print(f"\n{'='*80}")
print(f"CHIMERIC j-FPs (n_novel=0): {len(chimeric)}")
print(f"{'='*80}")
chimeric.sort(key=lambda r: -r['longcov'])
print(f"  {'tid':<14} {'gene':<14} {'src':<22} {'cov':>6} {'lc':>5} {'ni':>3} {'eq@locus':>9} {'best_ref'}")
for r in chimeric:
    print(f"  {r['tid']:<14} {r['gene']:<14} {r['source']:<22} {r['cov']:>6.2f} {r['longcov']:>5.1f} {r['n_juncs']:>3} {r['n_eq_at_locus']:>9} {r['best_ref']}")

# ---------- single-novel j-FPs ----------
single = [r for r in rows if r['n_novel'] == 1]
print(f"\n{'='*80}")
print(f"SINGLE-NOVEL j-FPs (n_novel=1): {len(single)}")
print(f"{'='*80}")
single.sort(key=lambda r: -r['longcov'])
print(f"  {'tid':<14} {'gene':<14} {'src':<22} {'cov':>6} {'lc':>5} {'ni':>3} {'novel_junc'}")
for r in single[:30]:  # top 30 by longcov
    nj = r['novel_juncs'][0]
    print(f"  {r['tid']:<14} {r['gene']:<14} {r['source']:<22} {r['cov']:>6.2f} {r['longcov']:>5.1f} {r['n_juncs']:>3}  {r['chrom']}:{nj[1]}-{nj[2]}")

# ---------- multi-novel j-FPs ----------
multi = [r for r in rows if r['n_novel'] >= 2]
print(f"\n{'='*80}")
print(f"MULTI-NOVEL j-FPs (n_novel>=2): {len(multi)}")
print(f"{'='*80}")
multi.sort(key=lambda r: -r['n_novel'])
print(f"  {'tid':<14} {'gene':<14} {'src':<22} {'cov':>6} {'lc':>5} {'ni':>3} {'n_novel':>7}")
for r in multi:
    print(f"  {r['tid']:<14} {r['gene']:<14} {r['source']:<22} {r['cov']:>6.2f} {r['longcov']:>5.1f} {r['n_juncs']:>3} {r['n_novel']:>7}")

# ---------- coverage distribution comparison j-FP vs TP ----------
tp_longcovs = [meta['longcov'] for tid, meta in tx_meta.items()
               if tx_class.get(tid) == '=' and len(meta['exons']) > 1]
fp_longcovs = [r['longcov'] for r in rows]

def percentile(lst, p):
    lst = sorted(lst)
    i = int(len(lst) * p / 100)
    return lst[min(i, len(lst)-1)]

print(f"\n{'='*80}")
print(f"LONGCOV DISTRIBUTION: j-FPs vs TPs")
print(f"{'='*80}")
for label, vals in [("j-FPs", fp_longcovs), ("TPs", tp_longcovs)]:
    print(f"  {label}: n={len(vals)}  p10={percentile(vals,10):.1f}  p25={percentile(vals,25):.1f}  "
          f"median={percentile(vals,50):.1f}  p75={percentile(vals,75):.1f}  "
          f"p90={percentile(vals,90):.1f}  max={max(vals):.1f}")

print(f"\n  j-FP longcov <= 1.0: {sum(1 for v in fp_longcovs if v <= 1.0)} / {len(fp_longcovs)}")
print(f"  TP  longcov <= 1.0: {sum(1 for v in tp_longcovs if v <= 1.0)} / {len(tp_longcovs)}")
print(f"\n  j-FP longcov <= 2.0: {sum(1 for v in fp_longcovs if v <= 2.0)} / {len(fp_longcovs)}")
print(f"  TP  longcov <= 2.0: {sum(1 for v in tp_longcovs if v <= 2.0)} / {len(tp_longcovs)}")

# ---------- cov distribution ----------
fp_covs = [r['cov'] for r in rows]
tp_covs = [meta['cov'] for tid, meta in tx_meta.items()
           if tx_class.get(tid) == '=' and len(meta['exons']) > 1]
print(f"\n{'='*80}")
print(f"COV DISTRIBUTION: j-FPs vs TPs")
print(f"{'='*80}")
for label, vals in [("j-FPs", fp_covs), ("TPs", tp_covs)]:
    print(f"  {label}: n={len(vals)}  p10={percentile(vals,10):.2f}  p25={percentile(vals,25):.2f}  "
          f"median={percentile(vals,50):.2f}  p75={percentile(vals,75):.2f}  "
          f"p90={percentile(vals,90):.2f}  max={max(vals):.2f}")

# ---------- locus complexity ----------
print(f"\n{'='*80}")
print(f"LOCUS COMPLEXITY: how many TPs exist at j-FP loci")
print(f"{'='*80}")
eq_at_locus_dist = Counter(r['n_eq_at_locus'] for r in rows)
for k in sorted(eq_at_locus_dist):
    print(f"  n_TPs_at_locus={k}: {eq_at_locus_dist[k]} j-FPs")

# ---------- checktrf vs flow by longcov threshold ----------
print(f"\n{'='*80}")
print(f"FILTER SIMULATION: what if we apply longcov > X threshold?")
print(f"{'='*80}")
for thresh in [1.0, 1.5, 2.0, 3.0, 5.0]:
    fp_removed = sum(1 for r in rows if r['longcov'] <= thresh)
    tp_removed = sum(1 for tid, meta in tx_meta.items()
                     if tx_class.get(tid) == '=' and len(meta['exons']) > 1
                     and meta['longcov'] <= thresh)
    print(f"  longcov > {thresh:.1f}: removes {fp_removed} j-FPs, {tp_removed} TPs  "
          f"(FP:TP ratio = {fp_removed/(tp_removed+1e-9):.1f}x)")

print(f"\n{'='*80}")
print(f"FILTER SIMULATION: by source + longcov threshold")
print(f"{'='*80}")
for src in ['checktrf_rescue', 'flow']:
    for thresh in [1.0, 1.5, 2.0]:
        fp_src = sum(1 for r in rows if r['source'] == src and r['longcov'] <= thresh)
        tp_src = sum(1 for tid, meta in tx_meta.items()
                     if tx_class.get(tid) == '=' and len(meta['exons']) > 1
                     and meta['source'] == src and meta['longcov'] <= thresh)
        print(f"  {src} longcov>{thresh:.1f}: removes {fp_src} j-FPs, {tp_src} TPs  "
              f"(ratio {fp_src/(tp_src+1e-9):.2f}x)")
