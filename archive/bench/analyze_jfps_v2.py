#!/usr/bin/env python3
"""
Enhanced j-FP analysis with junction mm values from parity log.
For each j-FP: source, cov, longcov, n_novel, entry_abund, min_mm, novel_junc_mm.
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
    v = parse_attr(s, k)
    return float(v) if v else None

# -------- load reference junctions --------
ref_juncs = set()   # (donor, acceptor) — using GTF exon-boundary convention
print("Loading reference GTF...", file=sys.stderr)
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
        elif feat == 'exon':
            cur_exons.append((s, e))
if len(cur_exons) > 1:
    cur_exons.sort()
    for i in range(len(cur_exons)-1):
        ref_juncs.add((cur_exons[i][1], cur_exons[i+1][0]))
print(f"  {len(ref_juncs)} reference junctions", file=sys.stderr)

# -------- load Rustle GTF --------
print("Loading Rustle GTF...", file=sys.stderr)
tx_meta = {}
cur_tid = None; cur_exons = []
with open(RGTF) as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip().split('\t')
        if len(f) < 9: continue
        feat, s, e, attrs = f[2], int(f[3])-1, int(f[4]), f[8]
        if feat == 'transcript':
            if cur_tid:
                cur_exons.sort()
                tx_meta[cur_tid]['exons'] = cur_exons
            cur_tid = parse_attr(attrs, 'transcript_id')
            tx_meta[cur_tid] = {
                'gene_id': parse_attr(attrs, 'gene_id'),
                'cov': parse_float(attrs, 'cov') or 0.0,
                'longcov': parse_float(attrs, 'longcov') or 0.0,
                'source': parse_attr(attrs, 'source') or '?',
                'strand': f[6], 'chrom': f[0],
                'exons': [],
            }
            cur_exons = []
        elif feat == 'exon' and cur_tid:
            cur_exons.append((s, e))
if cur_tid:
    cur_exons.sort(); tx_meta[cur_tid]['exons'] = cur_exons
print(f"  {len(tx_meta)} transcripts", file=sys.stderr)

# -------- load tmap --------
tx_class = {}; tx_best_ref = {}
with open(TMAP) as fh:
    for i, line in enumerate(fh):
        if i == 0: continue
        f = line.rstrip().split('\t')
        if len(f) < 5: continue
        tx_class[f[4]] = f[2]; tx_best_ref[f[4]] = f[0]

# -------- load parity log: junction mm and path_extracted --------
print("Loading parity log...", file=sys.stderr)
# junction_accept: key = (donor, acceptor) — from parity start/end
# The parity log uses start = j.donor (0-based first intron base)
# and end = j.acceptor+1 ... let's determine from context
# path_extracted introns format: "d1-a1,d2-a2,..." where d and a are 0-based?
# We'll figure out the offset empirically.
junc_mm = {}       # (start, end) → mm  [raw parity coords]
path_events = []   # list of path_extracted dicts

with open(PARITY) as fh:
    for line in fh:
        try:
            d = json.loads(line.strip())
        except: continue
        step = d.get('step','')
        if step == 'junction_accept':
            key = (d['start'], d['end'])
            junc_mm[key] = d['payload']['mm']
        elif step == 'path_extracted':
            path_events.append(d)

print(f"  {len(junc_mm)} junctions, {len(path_events)} path_extracted events", file=sys.stderr)

# Determine parity coordinate convention by checking a known junction
# from ref_juncs against junc_mm keys.
# Try a sample: find first junction that's in both ref_juncs and junc_mm.
# In ref_juncs: (donor, acceptor) where donor=exon_end(0-based-excl), acceptor=next_exon_start(0-based)
# parity log: start=donor, end=?
# Try (start, end) = (donor, acceptor) first, then (donor, acceptor+1)
sample_junc = next(iter(ref_juncs))
donor, acceptor = sample_junc
match_direct = (donor, acceptor) in junc_mm
match_plus1 = (donor, acceptor+1) in junc_mm
print(f"  Coord check - direct match: {match_direct}, +1 match: {match_plus1}", file=sys.stderr)
PARITY_OFFSET = 1 if match_plus1 else 0   # add this to acceptor when looking up junc_mm

def junc_mm_lookup(donor, acceptor):
    return junc_mm.get((donor, acceptor + PARITY_OFFSET), -1.0)

# -------- build intron-chain index from path_extracted --------
# path introns string: "d1-a1,d2-a2,..." where d=donor, a=acceptor in parity coords
def intron_chain_key(exons):
    """Compute intron chain key matching path_extracted 'introns' field format.
    path_extracted uses {d+1}-{a}: left = 1-based first intronic base, right = 0-based acceptor."""
    if len(exons) < 2: return ""
    parts = []
    for i in range(len(exons)-1):
        donor = exons[i][1]            # 0-based exclusive end of exon = first intronic base
        acceptor = exons[i+1][0]       # 0-based inclusive start of next exon
        parts.append(f"{donor+1}-{acceptor}")
    return ",".join(parts)

# Build lookup from intron chain to entry_abund
chain_to_entry_abund = {}
for ev in path_events:
    chain = ev['payload'].get('introns', '')
    if chain:
        chain_to_entry_abund[chain] = ev['payload'].get('entry_abund', -1.0)

# Verify with a known 2-exon transcript
sample_tid = next((tid for tid, m in tx_meta.items() if len(m['exons'])==2), None)
if sample_tid:
    chain = intron_chain_key(tx_meta[sample_tid]['exons'])
    hit = chain in chain_to_entry_abund
    print(f"  Chain lookup check (2-exon tx {sample_tid}): chain={chain} found={hit}", file=sys.stderr)

# -------- build j-FP rows --------
def tx_novel_juncs(meta):
    exons = meta['exons']
    novel = []
    for i in range(len(exons)-1):
        d = exons[i][1]; a = exons[i+1][0]
        if (d, a) not in ref_juncs:
            novel.append((d, a))
    return novel

def tx_all_juncs(meta):
    exons = meta['exons']
    return [(exons[i][1], exons[i+1][0]) for i in range(len(exons)-1)]

gene_eq_count = collections.Counter(
    meta['gene_id'] for tid, meta in tx_meta.items() if tx_class.get(tid) == '='
)

rows = []
for tid, meta in tx_meta.items():
    if tx_class.get(tid) != 'j': continue
    exons = meta['exons']
    if len(exons) < 2: continue
    novel = tx_novel_juncs(meta)
    all_juncs = tx_all_juncs(meta)
    chain = intron_chain_key(exons)
    entry_abund = chain_to_entry_abund.get(chain, -1.0)
    # mm for all junctions
    mm_vals = [junc_mm_lookup(d, a) for d, a in all_juncs]
    novel_mm = [junc_mm_lookup(d, a) for d, a in novel]
    min_mm = min(mm_vals) if mm_vals else -1.0
    rows.append({
        'tid': tid,
        'gene': meta['gene_id'],
        'source': meta['source'],
        'cov': meta['cov'],
        'longcov': meta['longcov'],
        'n_exons': len(exons),
        'n_juncs': len(all_juncs),
        'n_novel': len(novel),
        'novel_juncs': novel,
        'novel_mm': novel_mm,
        'min_mm_all': min_mm,
        'min_novel_mm': min(novel_mm) if novel_mm else -1.0,
        'entry_abund': entry_abund,
        'eq_at_locus': gene_eq_count[meta['gene_id']],
        'chrom': meta['chrom'],
    })

print(f"\n{'='*80}")
print(f"TOTAL j-FPs: {len(rows)}")
print(f"{'='*80}\n")

# -------- summary --------
from collections import Counter

def pct(n, d): return f"{100*n/d:.0f}%" if d else "0%"
def p(vals, ptile): return vals[int(len(vals)*ptile/100)] if vals else 0

rows_flow = [r for r in rows if r['source'] == 'flow']
rows_ctrf = [r for r in rows if r['source'] == 'checktrf_rescue']

print(f"By source: flow={len(rows_flow)}, checktrf={len(rows_ctrf)}")
print(f"By n_novel: " + ", ".join(f"{k}:{v}" for k,v in sorted(Counter(r['n_novel'] for r in rows).items())))
print()

# -------- entry_abund analysis --------
print("ENTRY_ABUND distribution (j-FPs):")
ea_vals = [r['entry_abund'] for r in rows if r['entry_abund'] >= 0]
ea_missing = sum(1 for r in rows if r['entry_abund'] < 0)
if ea_vals:
    ea_vals.sort()
    print(f"  n={len(ea_vals)} matched, {ea_missing} unmatched")
    print(f"  p10={p(ea_vals,10):.2f} p25={p(ea_vals,25):.2f} median={p(ea_vals,50):.2f} "
          f"p75={p(ea_vals,75):.2f} p90={p(ea_vals,90):.2f} max={max(ea_vals):.2f}")
    print(f"  entry_abund < 2: {sum(1 for v in ea_vals if v < 2)} ({pct(sum(1 for v in ea_vals if v < 2), len(ea_vals))})")
    print(f"  entry_abund < 5: {sum(1 for v in ea_vals if v < 5)} ({pct(sum(1 for v in ea_vals if v < 5), len(ea_vals))})")

# TP entry_abund comparison
tp_chains = set()
for tid, meta in tx_meta.items():
    if tx_class.get(tid) == '=' and len(meta['exons']) >= 2:
        chain = intron_chain_key(meta['exons'])
        tp_chains.add(chain)
tp_ea = [chain_to_entry_abund[c] for c in tp_chains if c in chain_to_entry_abund]
tp_ea.sort()
print(f"\nENTRY_ABUND distribution (TPs):")
if tp_ea:
    print(f"  n={len(tp_ea)} matched")
    print(f"  p10={p(tp_ea,10):.2f} p25={p(tp_ea,25):.2f} median={p(tp_ea,50):.2f} "
          f"p75={p(tp_ea,75):.2f} p90={p(tp_ea,90):.2f} max={max(tp_ea):.2f}")

# -------- min_mm analysis --------
print(f"\nMIN_MM (across ALL junctions in j-FPs):")
mm_vals = [r['min_mm_all'] for r in rows if r['min_mm_all'] >= 0]
mm_vals.sort()
if mm_vals:
    print(f"  n={len(mm_vals)}  min={mm_vals[0]:.0f}  "
          f"p25={p(mm_vals,25):.0f}  median={p(mm_vals,50):.0f}  "
          f"p75={p(mm_vals,75):.0f}  max={max(mm_vals):.0f}")
    for thr in [1, 2, 3, 5, 10]:
        cnt = sum(1 for v in mm_vals if v < thr)
        print(f"  min_mm < {thr}: {cnt} j-FPs ({pct(cnt, len(mm_vals))})")

print(f"\nMIN_NOVEL_MM (across novel junctions only, n_novel>0 cases):")
novel_rows = [r for r in rows if r['n_novel'] > 0 and r['min_novel_mm'] >= 0]
novel_mm_vals = [r['min_novel_mm'] for r in novel_rows]
novel_mm_vals.sort()
if novel_mm_vals:
    print(f"  n={len(novel_mm_vals)}")
    print(f"  min={min(novel_mm_vals):.0f}  p25={p(novel_mm_vals,25):.0f}  "
          f"median={p(novel_mm_vals,50):.0f}  p75={p(novel_mm_vals,75):.0f}  max={max(novel_mm_vals):.0f}")
    for thr in [1, 2, 3, 5]:
        cnt = sum(1 for v in novel_mm_vals if v < thr)
        print(f"  min_novel_mm < {thr}: {cnt} j-FPs ({pct(cnt, len(novel_mm_vals))})")

# TP min_mm for comparison
tp_min_mms = []
for tid, meta in tx_meta.items():
    if tx_class.get(tid) != '=': continue
    juncs = tx_all_juncs(meta)
    mms = [junc_mm_lookup(d, a) for d, a in juncs]
    if mms and min(mms) >= 0:
        tp_min_mms.append(min(mms))
tp_min_mms.sort()
print(f"\nMIN_MM distribution (TPs):")
if tp_min_mms:
    print(f"  n={len(tp_min_mms)}  min={tp_min_mms[0]:.0f}  "
          f"p25={p(tp_min_mms,25):.0f}  median={p(tp_min_mms,50):.0f}  max={max(tp_min_mms):.0f}")
    for thr in [1, 2, 3, 5]:
        cnt = sum(1 for v in tp_min_mms if v < thr)
        print(f"  TP min_mm < {thr}: {cnt} ({pct(cnt, len(tp_min_mms))})")

# -------- filter simulation with mm --------
print(f"\n{'='*80}")
print(f"FILTER SIMULATION: novel junction quality (n_novel>0 j-FPs only)")
print(f"{'='*80}")
novel_count = sum(1 for r in rows if r['n_novel'] > 0)
print(f"Total novel j-FPs: {novel_count}")
for mm_thr in [1, 2, 3, 5]:
    fp_removed = sum(1 for r in rows if r['n_novel'] > 0 and 0 <= r['min_novel_mm'] < mm_thr)
    # TPs that would be lost: TPs with any novel junction below threshold
    # (TPs by definition have all junctions in reference, so min_novel_mm is irrelevant)
    # BUT: if we apply this filter to ALL transcripts with novel junctions, we'd also remove
    # TPs that happen to have junctions "novel relative to reference" — but TPs are exact matches,
    # so they have no novel junctions. This filter only applies to j-class transcripts.
    print(f"  novel_mm < {mm_thr}: removes {fp_removed} j-FPs, 0 TPs (TPs have no novel junctions)")

print(f"\n{'='*80}")
print(f"FILTER SIMULATION: min_mm across all junctions (affects chimeric + novel)")
print(f"{'='*80}")
for mm_thr in [2, 3, 5]:
    fp_removed = sum(1 for r in rows if 0 <= r['min_mm_all'] < mm_thr)
    tp_removed = sum(1 for v in tp_min_mms if v < mm_thr)
    total_tp = len(tp_min_mms)
    print(f"  min_mm < {mm_thr}: removes {fp_removed} j-FPs ({pct(fp_removed,len(rows))}), "
          f"{tp_removed} TPs ({pct(tp_removed, total_tp)})  ratio={fp_removed/(tp_removed+1e-9):.2f}x")

# -------- entry_abund + min_mm combined --------
print(f"\n{'='*80}")
print(f"COMBINED FILTER: entry_abund < X AND min_mm < Y")
print(f"{'='*80}")
for ea_thr in [2.0, 3.0, 5.0]:
    for mm_thr in [2, 3, 5]:
        fp_rm = sum(1 for r in rows if r['entry_abund'] >= 0 and r['entry_abund'] < ea_thr
                    and r['min_mm_all'] >= 0 and r['min_mm_all'] < mm_thr)
        # TPs removed: TPs with entry_abund < ea_thr AND min_mm < mm_thr
        tp_rm = 0
        for tid, meta in tx_meta.items():
            if tx_class.get(tid) != '=': continue
            chain = intron_chain_key(meta['exons'])
            ea = chain_to_entry_abund.get(chain, -1.0)
            juncs = tx_all_juncs(meta)
            mms = [junc_mm_lookup(d, a) for d, a in juncs]
            min_mm = min(mms) if mms else -1.0
            if 0 <= ea < ea_thr and 0 <= min_mm < mm_thr:
                tp_rm += 1
        print(f"  ea<{ea_thr:.0f} AND min_mm<{mm_thr}: removes {fp_rm} j-FPs, {tp_rm} TPs  ratio={fp_rm/(tp_rm+1e-9):.2f}x")

# -------- detailed table for top j-FPs by longcov --------
print(f"\n{'='*80}")
print(f"TOP 30 j-FPs BY LONGCOV (with mm and entry_abund)")
print(f"{'='*80}")
top = sorted(rows, key=lambda r: -r['longcov'])[:30]
print(f"  {'tid':<16} {'src':<18} {'cov':>6} {'lc':>5} {'ea':>5} {'ni':>3} {'nv':>3} {'min_mm':>7} {'novel_mm'}")
for r in top:
    nv_mm_str = ','.join(f"{v:.0f}" for v in r['novel_mm']) if r['novel_mm'] else '-'
    print(f"  {r['tid']:<16} {r['source']:<18} {r['cov']:>6.2f} {r['longcov']:>5.1f} "
          f"{r['entry_abund']:>5.1f} {r['n_juncs']:>3} {r['n_novel']:>3} {r['min_mm_all']:>7.0f} {nv_mm_str}")

# -------- multi-novel j-FPs detail --------
print(f"\n{'='*80}")
print(f"MULTI-NOVEL j-FPs (n_novel >= 2) — full detail")
print(f"{'='*80}")
multi = sorted([r for r in rows if r['n_novel'] >= 2], key=lambda r: -r['longcov'])
print(f"  {'tid':<16} {'src':<18} {'lc':>5} {'nv':>3} {'min_novel_mm':>13}")
for r in multi:
    nv_mm = [f"{v:.0f}" for v in r['novel_mm']]
    print(f"  {r['tid']:<16} {r['source']:<18} {r['longcov']:>5.1f} {r['n_novel']:>3} "
          f"  novel_mm=[{','.join(nv_mm)}]  {r['gene']}")

# -------- chimeric j-FPs: how many locus TPs? --------
print(f"\n{'='*80}")
print(f"CHIMERIC j-FPs: locus complexity vs coverage")
print(f"{'='*80}")
chim = sorted([r for r in rows if r['n_novel'] == 0], key=lambda r: -r['longcov'])
print(f"  {'tid':<16} {'src':<18} {'lc':>5} {'cov':>6} {'ea':>5} {'ni':>3} {'eq_at_locus':>12} {'min_mm':>7}")
for r in chim:
    print(f"  {r['tid']:<16} {r['source']:<18} {r['longcov']:>5.1f} {r['cov']:>6.2f} "
          f"{r['entry_abund']:>5.1f} {r['n_juncs']:>3} {r['eq_at_locus']:>12} {r['min_mm_all']:>7.0f}")

print(f"\nChimeric j-FPs by locus-TPs vs longcov<=1:")
for eq_thresh in [0, 1, 3, 5, 10]:
    n_lc1 = sum(1 for r in chim if r['eq_at_locus'] >= eq_thresh and r['longcov'] <= 1.0)
    n_all = sum(1 for r in chim if r['eq_at_locus'] >= eq_thresh)
    print(f"  eq_at_locus>={eq_thresh}: {n_all} chimeric j-FPs, {n_lc1} have longcov<=1.0")
