#!/usr/bin/env python3
"""
analyze_jfps_v4.py
For each j-FP, find the pipeline stage where StringTie diverges.

  JUNCTION_DEMOTE  - ST did NOT accept one or more of the j-FP's junctions
                     (the junction is in Rustle's accepted set but absent from ST's)
  PATH_ENUM        - ST accepted all junctions but never built this intron chain
                     (chain in Rustle path_extracted but not ST path_extracted)
  PRED_FILTER      - ST built the same chain but pred-killed it; Rustle did not

For JUNCTION_DEMOTE we distinguish novel (not-in-reference) junctions from
reference junctions (parity gap if reference junctions are demoted by ST).
For novel-demoted junctions we compare the Rustle mm vs ST nreads_good.
"""

import json, sys, re, collections

TMAP      = "se_last_v3_cmp.se_last_v3.gtf.tmap"
RGTF      = "se_last_v3.gtf"
REFGTF    = "/mnt/c/Users/jfris/Desktop/GGO_19.gtf"
RU_PARITY = "strg343_depl/rustle_selast3.jsonl"
ST_PARITY = "jfp_study_st_parity.jsonl"


def parse_attr(s, k):
    m = re.search(rf'{k}\s+"([^"]+)"', s)
    return m.group(1) if m else None

def parse_float(s, k):
    v = parse_attr(s, k); return float(v) if v else None


# ---------- load reference junctions ----------
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
        elif feat == 'exon':
            cur_exons.append((s, e))
if len(cur_exons) > 1:
    cur_exons.sort()
    for i in range(len(cur_exons)-1):
        ref_juncs.add((cur_exons[i][1], cur_exons[i+1][0]))
print(f"Reference junctions: {len(ref_juncs)}", file=sys.stderr)


# ---------- load Rustle GTF ----------
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
                'cov':     parse_float(attrs, 'cov') or 0.0,
                'longcov': parse_float(attrs, 'longcov') or 0.0,
                'source':  parse_attr(attrs, 'source') or '?',
                'strand':  f[6], 'chrom': f[0], 'exons': []
            }; cur_exons = []
        elif feat == 'exon' and cur_tid:
            cur_exons.append((s, e))
if cur_tid:
    cur_exons.sort(); tx_meta[cur_tid]['exons'] = cur_exons
print(f"Rustle transcripts: {len(tx_meta)}", file=sys.stderr)


def intron_chain_key(exons):
    """path_extracted format: '{d+1}-{a}' per junction, comma-joined."""
    if len(exons) < 2: return ""
    parts = []
    for i in range(len(exons)-1):
        d = exons[i][1]       # 0-based exclusive end of exon = first intronic base
        a = exons[i+1][0]     # 0-based inclusive start of next exon
        parts.append(f"{d+1}-{a}")
    return ",".join(parts)

def get_junctions(exons):
    """Returns (donor, acceptor) pairs: donor=exon[i].end (0-based), acceptor=exon[i+1].start (0-based)."""
    if len(exons) < 2: return []
    return [(exons[i][1], exons[i+1][0]) for i in range(len(exons)-1)]

def junc_key(d, a):
    """Parity-log junction key: (donor_0based, acceptor_0based+1). PARITY_OFFSET=1."""
    return (d, a + 1)


# ---------- load parity logs ----------
# Rustle: junc_mm[(key)] = mm   (mm = anchor-filtered read count)
# ST:     st_junc[(key)] = (mm, nreads_good)
ru_junc    = {}  # key -> mm
ru_chain_ea = {}  # chain -> entry_abund
st_junc    = {}  # key -> (mm, nreads_good)
st_chain_ea = {}  # chain -> entry_abund

print("Loading Rustle parity log...", file=sys.stderr)
with open(RU_PARITY) as fh:
    for line in fh:
        ev = json.loads(line)
        step = ev.get("step", "")
        if step == "junction_accept":
            p = ev.get("payload", {})
            if p.get("accepted", False):
                ru_junc[(ev["start"], ev["end"])] = p.get("mm", 0)
        elif step == "path_extracted":
            chain = ev.get("payload", {}).get("introns", "")
            ea    = ev.get("payload", {}).get("entry_abund",
                          ev.get("payload", {}).get("cov", 0))
            if chain:
                ru_chain_ea[chain] = ea

print("Loading ST parity log...", file=sys.stderr)
with open(ST_PARITY) as fh:
    for line in fh:
        ev = json.loads(line)
        step = ev.get("step", "")
        if step == "junction_accept":
            s, e = ev.get("start", 0), ev.get("end", 0)
            if s == 0 and e == 0: continue  # padding entry
            p = ev.get("payload", {})
            if p.get("accepted", False):
                st_junc[(s, e)] = (p.get("mm", 0), p.get("nreads_good", 0))
        elif step == "path_extracted":
            chain = ev.get("payload", {}).get("introns", "")
            ea    = ev.get("payload", {}).get("entry_abund",
                          ev.get("payload", {}).get("cov", 0))
            if chain:
                st_chain_ea[chain] = ea

print(f"  Rustle: {len(ru_junc)} accepted junctions, {len(ru_chain_ea)} chains", file=sys.stderr)
print(f"  ST:     {len(st_junc)} accepted junctions, {len(st_chain_ea)} chains", file=sys.stderr)


# ---------- load tmap ----------
jfp_tids = set()
tp_tids  = set()
with open(TMAP) as fh:
    next(fh)
    for line in fh:
        p = line.rstrip().split('\t')
        if len(p) < 6: continue
        cls  = p[2]   # class_code
        qtid = p[4]   # qry_id
        try: nexons = int(p[5])   # num_exons is column 5 (0-indexed)
        except: continue
        if nexons <= 1: continue
        if cls == 'j':   jfp_tids.add(qtid)
        elif cls == '=': tp_tids.add(qtid)
print(f"j-FPs: {len(jfp_tids)}, TPs: {len(tp_tids)}", file=sys.stderr)


# ---------- classify each j-FP ----------
categories = collections.defaultdict(list)

for tid in sorted(jfp_tids):
    if tid not in tx_meta: continue
    tx    = tx_meta[tid]
    exns  = tx['exons']
    juncs = get_junctions(exns)
    chain = intron_chain_key(exns)
    novel_juncs = [(d, a) for (d, a) in juncs if (d, a) not in ref_juncs]

    # Check each junction: was it accepted by ST?
    novel_demoted = []
    ref_demoted   = []
    for d, a in juncs:
        k = junc_key(d, a)
        if k not in st_junc:
            is_novel = (d, a) not in ref_juncs
            ru_mm = ru_junc.get(k, -1)
            rec = {"d": d, "a": a, "ru_mm": ru_mm}
            if is_novel:
                novel_demoted.append(rec)
            else:
                ref_demoted.append(rec)

    if novel_demoted or ref_demoted:
        categories["JUNCTION_DEMOTE"].append({
            "tid": tid, "chain": chain,
            "src": tx['source'], "cov": tx['cov'], "lc": tx['longcov'],
            "n_juncs": len(juncs), "n_novel": len(novel_juncs),
            "novel_demoted": novel_demoted,
            "ref_demoted":   ref_demoted,
            "ru_ea": ru_chain_ea.get(chain, -1),
        })
    elif chain in st_chain_ea:
        categories["PRED_FILTER"].append({
            "tid": tid, "chain": chain,
            "src": tx['source'], "cov": tx['cov'], "lc": tx['longcov'],
            "n_novel": len(novel_juncs),
            "st_ea": st_chain_ea[chain],
            "ru_ea": ru_chain_ea.get(chain, -1),
        })
    else:
        # ST accepted all junctions but never built this path
        n_novel_st_acc = sum(1 for (d,a) in novel_juncs if junc_key(d,a) in st_junc)
        novel_mm = [ru_junc.get(junc_key(d,a), -1) for (d,a) in novel_juncs]
        categories["PATH_ENUM"].append({
            "tid": tid, "chain": chain,
            "src": tx['source'], "cov": tx['cov'], "lc": tx['longcov'],
            "n_novel": len(novel_juncs),
            "n_novel_st_accepted": n_novel_st_acc,
            "novel_mm": novel_mm,
            "ru_ea": ru_chain_ea.get(chain, -1),
        })


# =====================================================================
print("\n" + "="*80)
print("j-FP DIVERGENCE POINT CLASSIFICATION")
print("="*80)

total = sum(len(v) for v in categories.values())
print(f"\nTotal j-FPs classified: {total}")
for cat in ["JUNCTION_DEMOTE", "PATH_ENUM", "PRED_FILTER"]:
    n = len(categories.get(cat, []))
    pct = 100*n/total if total else 0
    print(f"  {cat:20s}: {n:4d}  ({pct:.0f}%)")


# ---- JUNCTION_DEMOTE detail ----
jd = categories.get("JUNCTION_DEMOTE", [])
if jd:
    print(f"\n{'='*80}")
    print(f"JUNCTION_DEMOTE ({len(jd)} j-FPs)")
    print(f"{'='*80}")
    print("ST did NOT accept one or more of the j-FP's junctions.\n")

    only_novel  = [x for x in jd if x["novel_demoted"] and not x["ref_demoted"]]
    only_ref    = [x for x in jd if x["ref_demoted"]   and not x["novel_demoted"]]
    mixed       = [x for x in jd if x["novel_demoted"] and x["ref_demoted"]]
    print(f"  Only novel junctions demoted: {len(only_novel)}")
    print(f"  Only ref junctions demoted:   {len(only_ref)}")
    print(f"  Both novel+ref demoted:       {len(mixed)}")

    # mm distribution of demoted novel junctions
    all_novel_mm = []
    for item in jd:
        for nd in item["novel_demoted"]:
            if nd["ru_mm"] >= 0:
                all_novel_mm.append(nd["ru_mm"])
    if all_novel_mm:
        all_novel_mm.sort()
        def pct(vals, p):
            i = max(0, int(len(vals)*p/100)-1)
            return vals[min(i, len(vals)-1)]
        print(f"\n  Rustle-mm of novel junctions demoted by ST (n={len(all_novel_mm)}):")
        print(f"    min={int(all_novel_mm[0])} p10={int(pct(all_novel_mm,10))} "
              f"p25={int(pct(all_novel_mm,25))} median={int(pct(all_novel_mm,50))} "
              f"p75={int(pct(all_novel_mm,75))} max={int(all_novel_mm[-1])}")

    # Check if demoted junctions appear in ST's log as accepted=false
    # (ST might have logged them with accepted=false — we only stored accepted=true above)
    # Re-scan to find any ST events for these junction keys
    demoted_keys = set()
    for item in jd:
        for nd in item["novel_demoted"] + item["ref_demoted"]:
            demoted_keys.add(junc_key(nd["d"], nd["a"]))

    st_demoted_detail = {}  # key -> payload
    with open(ST_PARITY) as fh:
        for line in fh:
            ev = json.loads(line)
            if ev.get("step") == "junction_accept":
                s, e = ev.get("start", 0), ev.get("end", 0)
                if (s, e) in demoted_keys:
                    st_demoted_detail[(s, e)] = ev.get("payload", {})

    if st_demoted_detail:
        print(f"\n  ST log entries found for {len(st_demoted_detail)} of {len(demoted_keys)} demoted junction keys:")
        # Show what ST said about them
        found_accepted = sum(1 for p in st_demoted_detail.values() if p.get("accepted"))
        found_rejected = sum(1 for p in st_demoted_detail.values() if not p.get("accepted"))
        print(f"    ST accepted: {found_accepted}   ST rejected: {found_rejected}")
        if found_rejected:
            print(f"\n    Sample ST-rejected entries:")
            shown = 0
            for (s,e), p in st_demoted_detail.items():
                if not p.get("accepted") and shown < 5:
                    print(f"      ({s},{e}) mm={p.get('mm')} nreads={p.get('nreads')} nreads_good={p.get('nreads_good')} strand={p.get('strand')}")
                    shown += 1
    else:
        print(f"\n  None of the {len(demoted_keys)} demoted junction keys appear in ST's log at all")
        print(f"  => These junctions were never seen by ST (different BAM parsing or below ST's scan threshold)")

    # Per-j-FP table
    print(f"\n  {'tid':<18} {'src':<18} {'lc':>5} {'cov':>7} nj nd_nov nd_ref  ru_ea  demoted_ru_mm")
    print(f"  {'-'*100}")
    for item in sorted(jd, key=lambda x: min([nd['ru_mm'] for nd in x['novel_demoted']] or [9999])):
        nd_str = ";".join(f"mm={nd['ru_mm']}" for nd in item["novel_demoted"]) or "-"
        rd_n   = len(item["ref_demoted"])
        print(f"  {item['tid']:<18} {item['src']:<18} {item['lc']:>5.1f} {item['cov']:>7.2f}"
              f"  {item['n_juncs']:2d} {len(item['novel_demoted']):>6} {rd_n:>6}"
              f"  {item['ru_ea']:>6.2f}  {nd_str}")


# ---- PATH_ENUM detail ----
pe = categories.get("PATH_ENUM", [])
if pe:
    print(f"\n{'='*80}")
    print(f"PATH_ENUM ({len(pe)} j-FPs)")
    print(f"{'='*80}")
    print("ST accepted all junctions but never built this intron chain.\n")

    chimeric = [x for x in pe if x["n_novel"] == 0]
    w_novel  = [x for x in pe if x["n_novel"] > 0]
    print(f"  Chimeric (all-reference junctions, novel combination): {len(chimeric)}")
    print(f"  With novel junctions ST accepted:                       {len(w_novel)}")
    if w_novel:
        n_all_acc = sum(1 for x in w_novel if x["n_novel_st_accepted"] == x["n_novel"])
        print(f"    ...of which all novel junctions accepted by ST:        {n_all_acc}")
        n_partial  = len(w_novel) - n_all_acc
        print(f"    ...partial acceptance (some novel juncs accepted by ST): {n_partial}")

    lc_vals = sorted(x["lc"] for x in pe)
    ea_vals = sorted(x["ru_ea"] for x in pe if x["ru_ea"] >= 0)
    p50 = lambda v: v[len(v)//2] if v else -1
    print(f"\n  longcov distribution: median={p50(lc_vals):.1f}  "
          f"values={lc_vals[:30]}{'...' if len(lc_vals)>30 else ''}")
    print(f"  ru_ea distribution:  median={p50(ea_vals):.2f}  "
          f"min={ea_vals[0]:.2f}  max={ea_vals[-1]:.2f}" if ea_vals else "  no ea data")

    if chimeric:
        print(f"\n  Chimeric j-FPs (ST built different path combinations at same locus):")
        print(f"  {'tid':<18} {'src':<18} {'lc':>5} {'cov':>7}  ru_ea")
        print(f"  {'-'*65}")
        for item in sorted(chimeric, key=lambda x: x["cov"], reverse=True):
            print(f"  {item['tid']:<18} {item['src']:<18} {item['lc']:>5.1f} {item['cov']:>7.2f} {item['ru_ea']:>7.2f}")

    if w_novel:
        print(f"\n  Novel-junction j-FPs where ST accepted junctions but not the path:")
        print(f"  {'tid':<18} {'src':<18} {'lc':>5} {'cov':>7} n_novel n_st_acc  novel_mm")
        print(f"  {'-'*85}")
        for item in sorted(w_novel, key=lambda x: min(x["novel_mm"]) if x["novel_mm"] else 999):
            mm_str = ",".join(str(int(m)) for m in item["novel_mm"])
            print(f"  {item['tid']:<18} {item['src']:<18} {item['lc']:>5.1f} {item['cov']:>7.2f}"
                  f"  {item['n_novel']:>4}    {item['n_novel_st_accepted']:>4}    {mm_str}")


# ---- PRED_FILTER detail ----
pf = categories.get("PRED_FILTER", [])
if pf:
    print(f"\n{'='*80}")
    print(f"PRED_FILTER ({len(pf)} j-FPs)")
    print(f"{'='*80}")
    print("ST built the same chain but pred-killed it; Rustle kept it.\n")
    print(f"  {'tid':<18} {'src':<18} {'lc':>5} {'cov':>7} n_novel  st_ea  ru_ea")
    print(f"  {'-'*80}")
    for item in sorted(pf, key=lambda x: x["cov"], reverse=True):
        print(f"  {item['tid']:<18} {item['src']:<18} {item['lc']:>5.1f} {item['cov']:>7.2f}"
              f"  {item['n_novel']:>4}  {item['st_ea']:>7.2f} {item['ru_ea']:>7.2f}")
