#!/usr/bin/env python3
"""retained_intron_geometry_oracle.py

CHEAP analysis-only oracle that bounds the OVER-KILL risk of a proposed
`retainedintron_like` PREDICATE-GEOMETRY change in Rustle, BEFORE any code is
written. No production code is touched and no behavior changes.

THE PROPOSED PREDICATE (the per-locus-traced miss)
--------------------------------------------------
StringTie's pairwise `retained_intron` filter kills a victim V when V's FIRST
exon lies fully INSIDE a higher-coverage local killer K's LOW-COVERAGE intron
(and ends before K's next exon) -- a "spurious-donor-within-intron" /
contained-first-exon case. Rustle's current predicate misses this. We want to
know whether adding that kill would realize ST's measured net +25 (27 RU-only
FP chains + 2 RU TPs killed) at LOW TP cost, or whether it would OVER-KILL real
alternative isoforms (TP loss). A cruder pairing attempt already cost -4 TPs;
that is the warning we are bounding against.

WHAT THIS TOOL DOES
-------------------
1. Loads RU final predictions from the RU GTF: per transcript ->
   (strand, sorted exons, transcript cov, FPKM/longcov, per-exon cov).
2. Determines per-intron LOW-NESS for each potential killer K from RU's own
   `pred_intron_low` events (the PREFERRED source). pred_intron_low carries no
   chain, so each event is matched to a GTF transcript by
   (exact span [start+1,end], strand, intron-count, cov within 0.02). This is
   unique for 1879/1903 multi-exon RU transcripts (see report). The positional
   `intron_low` mask is then aligned onto the transcript's introns in genomic
   order. (Fallback, only if a K has no matched mask: derive low-ness from the
   intron being absent-as-exon in K but present-as-exon in a higher-cov
   overlapping prediction. The tool reports how many Ks used each source.)
3. For every ORDERED same-strand overlapping pair (K, V) with cov(K) > cov(V):
   for each LOW intron [d,a] of K, fires the kill if V's relevant exon is fully
   contained: vs >= d and ve <= a. Two variants are computed:
     - CONTAINED-FIRST-EXON: only V's FIRST exon tested (the specific predicate).
     - ANY-EXON-CONTAINED: any V exon tested (broader; shows over-kill blow-up).
4. KILL SET = union of victims. Each V classified FP vs TP by whether its
   (strand, intron-chain) appears in the reference GTF. Intron-token convention:
   (exon[i-1].end+1, exon[i].start-1).  Report |kill|, FP (prize), TP (cost),
   NET = FP - TP, vs ST's +25.
5. SANITY: ST's retained_intron pred_kill `chain`s use the SAME intron-token
   convention (verified: 3205/3422 ST kill tokens match RU introns under
   (end+1,start-1), 0 under donor-acceptor). Intersect ST's killed-RU-FP set
   with the geometry kill set; report the overlap and the EXTRA (over-kill)
   chains the geometry kills beyond ST, split FP vs TP.

Stdlib only.
"""

import json
import sys
import collections

RU_GTF   = "/tmp/ri_ru.gtf"
RU_JSONL = "/tmp/ri_ru.jsonl"
ST_JSONL = "/tmp/ri_st.jsonl"
REF_GTF  = "/mnt/c/Users/jfris/Desktop/GGO_19.gtf"


# ---------------------------------------------------------------------------
# GTF parsing
# ---------------------------------------------------------------------------
def _attr(attrs, key):
    needle = key + ' "'
    if needle not in attrs:
        return None
    return attrs.split(needle, 1)[1].split('"', 1)[0]


def load_gtf(path):
    """Return {tid: {'chrom','strand','exons':[(s,e)...] sorted,
                     'cov':float|None,'fpkm','longcov','excov':[..]}}."""
    tx = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("transcript", "exon"):
            continue
        tid = _attr(f[8], "transcript_id")
        if tid is None:
            continue
        rec = tx.setdefault(tid, {"chrom": f[0], "strand": f[6],
                                  "exons": [], "cov": None, "fpkm": None,
                                  "longcov": None, "excov": []})
        if f[2] == "transcript":
            c = _attr(f[8], "cov")
            rec["cov"] = float(c) if c is not None else None
            fp = _attr(f[8], "FPKM")
            rec["fpkm"] = float(fp) if fp is not None else None
            lc = _attr(f[8], "longcov")
            rec["longcov"] = float(lc) if lc is not None else None
        else:  # exon
            rec["exons"].append((int(f[3]), int(f[4])))
            ec = _attr(f[8], "cov")
            rec["excov"].append(float(ec) if ec is not None else 0.0)
    # sort exons + parallel excov by genomic start
    for rec in tx.values():
        if rec["exons"]:
            order = sorted(range(len(rec["exons"])), key=lambda i: rec["exons"][i])
            rec["exons"] = [rec["exons"][i] for i in order]
            if len(rec["excov"]) == len(order):
                rec["excov"] = [rec["excov"][i] for i in order]
    return tx


def intron_tokens(exons):
    """(exon[i-1].end+1, exon[i].start-1) per intron, genomic order."""
    return [(exons[i - 1][1] + 1, exons[i][0] - 1) for i in range(1, len(exons))]


def chain_key(strand, exons):
    """Canonical (strand, intron-chain) reference-matching key."""
    return (strand, tuple(intron_tokens(exons)))


# ---------------------------------------------------------------------------
# Low-intron source: match RU pred_intron_low events -> GTF transcripts
# ---------------------------------------------------------------------------
def load_pred_intron_low(path):
    evs = []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        try:
            d = json.loads(line)
        except Exception:
            continue
        if d.get("step") != "pred_intron_low":
            continue
        p = d.get("payload", {})
        mask = [x for x in p.get("intron_low", "").split(",") if x != ""]
        icov = [x for x in p.get("intron_covs", "").split(",") if x != ""]
        ecov = [x for x in p.get("exon_covs", "").split(",") if x != ""]
        evs.append({
            "start": d["start"] + 1,   # jsonl start is 0-based; GTF is 1-based
            "end":   d["end"],
            "strand": d["strand"],
            "nintron": len(mask),
            "mask": [m == "1" for m in mask],
            "icov": [float(x) for x in icov] if icov else [],
            "ecov": [float(x) for x in ecov] if ecov else [],
            "cov":  p.get("cov"),
        })
    return evs


def build_low_intron_map(tx, evs):
    """For each multi-exon GTF tid, attach rec['low'] = list[bool] aligned to
    its introns (genomic order), and rec['low_src'] in {'mask','none'}.

    Match key: (exact span start+1/end, strand, intron-count, cov within 0.02).
    Reports counts of unique/multi/zero matches.
    """
    # index events for quick lookup
    by_key = collections.defaultdict(list)
    for ev in evs:
        by_key[(ev["start"], ev["end"], ev["strand"], ev["nintron"])].append(ev)

    stats = {"unique": 0, "multi": 0, "zero": 0, "mask_used": 0, "no_mask": 0}
    for tid, rec in tx.items():
        ni = len(rec["exons"]) - 1
        rec["low"] = None
        rec["low_src"] = "none"
        if ni < 1:
            continue
        cands = by_key.get((rec["exons"][0][0], rec["exons"][-1][1],
                            rec["strand"], ni), [])
        # disambiguate by cov when several share the span+strand+nintron
        if len(cands) > 1 and rec["cov"] is not None:
            tight = [e for e in cands
                     if e["cov"] is not None and abs(e["cov"] - rec["cov"]) < 0.02]
            if tight:
                cands = tight
        if len(cands) == 1:
            stats["unique"] += 1
            rec["low"] = list(cands[0]["mask"])
            rec["low_src"] = "mask"
            stats["mask_used"] += 1
        elif len(cands) > 1:
            stats["multi"] += 1
            # ambiguous: take the closest-cov candidate's mask (still mask-sourced)
            best = min(cands, key=lambda e: abs((e["cov"] or 0) - (rec["cov"] or 0)))
            rec["low"] = list(best["mask"])
            rec["low_src"] = "mask"
            stats["mask_used"] += 1
        else:
            stats["zero"] += 1
            stats["no_mask"] += 1
    return stats


def add_fallback_low(tx, fb_overlap_index):
    """Fallback for Ks with no matched mask: an intron of K is 'low' if some
    higher-cov overlapping prediction has an EXON covering that intron region
    (intron absent-as-exon in K but present-as-exon in a higher-cov pred).

    Returns count of Ks that received a fallback mask. Only applied where
    rec['low'] is None and the tx is multi-exon.
    """
    applied = 0
    for tid, rec in tx.items():
        if rec.get("low") is not None or len(rec["exons"]) < 2:
            continue
        introns = intron_tokens(rec["exons"])
        low = []
        kcov = rec["cov"] or 0.0
        for (d, a) in introns:
            covered = False
            for oid in fb_overlap_index(rec["chrom"], rec["strand"], d, a):
                o = tx[oid]
                if oid == tid:
                    continue
                if (o["cov"] or 0.0) <= kcov:
                    continue
                # any exon of o covering >50% of K's intron span?
                ilen = a - d + 1
                for (es, ee) in o["exons"]:
                    ov = min(a, ee) - max(d, es) + 1
                    if ov > 0 and ov >= 0.5 * ilen:
                        covered = True
                        break
                if covered:
                    break
            low.append(covered)
        rec["low"] = low
        rec["low_src"] = "fallback"
        applied += 1
    return applied


# ---------------------------------------------------------------------------
# Geometry kill computation
# ---------------------------------------------------------------------------
def compute_kill(tx, ref_chains, first_only):
    """Return (kill_set_tids, fp_chains, tp_chains, detail).
    first_only=True -> contained-FIRST-exon variant; False -> any-exon variant.
    """
    # bucket transcripts by (chrom, strand) for overlap pairing
    buckets = collections.defaultdict(list)
    for tid, rec in tx.items():
        if not rec["exons"]:
            continue
        buckets[(rec["chrom"], rec["strand"])].append(tid)

    kill = set()
    killed_by = {}  # vid -> (kid, intron)
    for key, tids in buckets.items():
        # sort by genomic start for windowed overlap
        tids.sort(key=lambda t: tx[t]["exons"][0][0])
        n = len(tids)
        for i in range(n):
            kid = tids[i]
            K = tx[kid]
            if len(K["exons"]) < 2 or not K.get("low"):
                continue
            kstart, kend = K["exons"][0][0], K["exons"][-1][1]
            kcov = K["cov"] if K["cov"] is not None else 0.0
            kintrons = intron_tokens(K["exons"])
            low = K["low"]
            low_introns = [kintrons[j] for j in range(len(kintrons))
                           if j < len(low) and low[j]]
            if not low_introns:
                continue
            for j in range(n):
                if j == i:
                    continue
                vid = tids[j]
                V = tx[vid]
                vstart, vend = V["exons"][0][0], V["exons"][-1][1]
                if vend < kstart or vstart > kend:
                    continue  # no overlap
                vcov = V["cov"] if V["cov"] is not None else 0.0
                if not (kcov > vcov):
                    continue  # K must be the higher-cov local killer
                # test exon containment in a low K intron
                test_exons = [V["exons"][0]] if first_only else V["exons"]
                fired = None
                for (vs, ve) in test_exons:
                    for (d, a) in low_introns:
                        if vs >= d and ve <= a:
                            fired = (d, a)
                            break
                    if fired:
                        break
                if fired:
                    kill.add(vid)
                    killed_by.setdefault(vid, (kid, fired))

    fp_chains = []
    tp_chains = []
    kill_chains = set()
    for vid in kill:
        ck = chain_key(tx[vid]["strand"], tx[vid]["exons"])
        kill_chains.add(ck)
        if ck in ref_chains:
            tp_chains.append((vid, ck))
        else:
            fp_chains.append((vid, ck))
    detail = {"killed_by": killed_by, "kill_chains": kill_chains}
    return kill, fp_chains, tp_chains, detail


# ---------------------------------------------------------------------------
# ST sanity set: ST's retained_intron-killed RU-FP chains
# ---------------------------------------------------------------------------
def load_st_retained_intron_chains(path):
    """Return set of (strand, intron-chain-tuple) ST killed via retained_intron,
    plus the raw count of pred_kill events."""
    chains = set()
    nev = 0
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        try:
            d = json.loads(line)
        except Exception:
            continue
        if d.get("step") != "pred_kill":
            continue
        p = d.get("payload", {})
        if p.get("reason") != "retained_intron":
            continue
        nev += 1
        ch = p.get("chain", "")
        if not ch:
            continue  # single-exon kill; no intron chain
        toks = []
        ok = True
        for t in ch.split(","):
            if "-" not in t:
                ok = False
                break
            a, b = t.split("-")
            toks.append((int(a), int(b)))
        if ok and toks:
            chains.add((d["strand"], tuple(toks)))
    return chains, nev


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 78)
    print("retained_intron geometry oracle  --  over-kill bound for the predicate fix")
    print("=" * 78)

    ru = load_gtf(RU_GTF)
    ref = load_gtf(REF_GTF)
    ru_multi = sum(1 for r in ru.values() if len(r["exons"]) > 1)
    print(f"\nRU final predictions: {len(ru)} transcripts ({ru_multi} multi-exon)")
    print(f"Reference (GGO_19) transcripts: {len(ref)}")

    # RU cov source
    n_cov = sum(1 for r in ru.values() if r["cov"] is not None)
    print(f"\n[RU coverage source] transcript-level cov attr present for "
          f"{n_cov}/{len(ru)} RU tx (from GTF `cov`).")
    if n_cov < len(ru):
        print("  WARNING: some RU tx lack a cov attr; those fall back to cov=0 "
              "(cannot be a killer; weak victim). Reported above.")

    # reference chain set
    ref_chains = set()
    for r in ref.values():
        if len(r["exons"]) >= 2:
            ref_chains.add(chain_key(r["strand"], r["exons"]))
    print(f"Reference multi-exon intron-chains: {len(ref_chains)}")

    # low-intron source
    evs = load_pred_intron_low(RU_JSONL)
    print(f"\n[Low-intron source] RU pred_intron_low events: {len(evs)} "
          "(NO chain field -- matched to GTF by span+strand+nintron+cov).")
    stats = build_low_intron_map(ru, evs)
    print(f"  unique mask match: {stats['unique']}  ambiguous(multi, cov-broken): "
          f"{stats['multi']}  no match: {stats['zero']}  (of {ru_multi} multi-exon)")

    # overlap index for fallback
    by_cs = collections.defaultdict(list)
    for tid, rec in ru.items():
        by_cs[(rec["chrom"], rec["strand"])].append(tid)

    def fb_overlap_index(chrom, strand, d, a):
        out = []
        for tid in by_cs.get((chrom, strand), []):
            r = ru[tid]
            if r["exons"][0][0] <= a and r["exons"][-1][1] >= d:
                out.append(tid)
        return out

    fb = add_fallback_low(ru, fb_overlap_index)
    print(f"  fallback (absent-as-exon vs higher-cov pred) applied to {fb} K(s) "
          "with no matched mask.")
    src_counts = collections.Counter(r.get("low_src", "none")
                                     for r in ru.values() if len(r["exons"]) > 1)
    print(f"  low-ness source per multi-exon tx: {dict(src_counts)}")

    # --- variant A: contained-FIRST-exon (the proposed predicate) ---
    print("\n" + "-" * 78)
    print("VARIANT A: contained-FIRST-exon  (the proposed predicate)")
    print("-" * 78)
    killA, fpA, tpA, detA = compute_kill(ru, ref_chains, first_only=True)
    netA = len(fpA) - len(tpA)
    print(f"  kill-set size : {len(killA)}")
    print(f"  FP (prize)    : {len(fpA)}")
    print(f"  TP (cost)     : {len(tpA)}")
    print(f"  NET (FP-TP)   : {netA:+d}   (ST's measured net = +25; 27 FP + 2 TP)")
    if tpA:
        print("  TP victims (real isoforms this predicate would lose):")
        for vid, ck in sorted(tpA):
            kid = detA["killed_by"][vid][0]
            print(f"    {vid}  killed_by={kid}  cov={ru[vid]['cov']}")

    # --- variant B: ANY-exon-contained (broader predicate) ---
    print("\n" + "-" * 78)
    print("VARIANT B: ANY-exon-contained  (broader predicate; over-kill demo)")
    print("-" * 78)
    killB, fpB, tpB, detB = compute_kill(ru, ref_chains, first_only=False)
    netB = len(fpB) - len(tpB)
    print(f"  kill-set size : {len(killB)}")
    print(f"  FP (prize)    : {len(fpB)}")
    print(f"  TP (cost)     : {len(tpB)}")
    print(f"  NET (FP-TP)   : {netB:+d}")

    # --- sanity vs ST ---
    print("\n" + "-" * 78)
    print("SANITY: geometry kill-set vs ST retained_intron kills")
    print("-" * 78)
    st_chains, st_nev = load_st_retained_intron_chains(ST_JSONL)
    print(f"  ST pred_kill(retained_intron) events: {st_nev}  "
          f"(multi-exon intron-chains: {len(st_chains)})")

    # ST's killed RU-FP set: ST chains that are (a) in RU final, (b) NOT in ref
    ru_chain_to_tids = collections.defaultdict(list)
    for tid, rec in ru.items():
        if len(rec["exons"]) >= 2:
            ru_chain_to_tids[chain_key(rec["strand"], rec["exons"])].append(tid)
    ru_final_chains = set(ru_chain_to_tids.keys())

    st_in_ru = st_chains & ru_final_chains
    st_killed_ru_fp = {c for c in st_in_ru if c not in ref_chains}
    st_killed_ru_tp = {c for c in st_in_ru if c in ref_chains}
    print(f"  ST chains present in RU final: {len(st_in_ru)} "
          f"(RU-FP: {len(st_killed_ru_fp)}, RU-TP: {len(st_killed_ru_tp)})")
    print("    (Baseline narrative: ST kills 27 RU-only FP chains + 2 RU TPs, "
          "net +25.)")

    for name, kill_chains in (("A (first-exon)", detA["kill_chains"]),
                              ("B (any-exon)", detB["kill_chains"])):
        inter = kill_chains & st_killed_ru_fp
        extra = kill_chains - st_in_ru          # geometry kills ST did NOT
        extra_fp = {c for c in extra if c not in ref_chains}
        extra_tp = {c for c in extra if c in ref_chains}
        recall = (len(inter) / len(st_killed_ru_fp) * 100.0
                  if st_killed_ru_fp else 0.0)
        print(f"\n  [Variant {name}]")
        print(f"    geometry kills inter ST's killed-RU-FP : "
              f"{len(inter)}/{len(st_killed_ru_fp)}  ({recall:.0f}% of ST's FP kills)")
        print(f"    EXTRA chains geometry kills beyond ST   : {len(extra)}")
        print(f"      of which FP (bonus prize) : {len(extra_fp)}")
        print(f"      of which TP (OVER-KILL)   : {len(extra_tp)}")
        if extra_tp and name.startswith("A"):
            print("      TP over-kills (real isoforms ST keeps but geometry "
                  "would drop):")
            for c in sorted(extra_tp):
                for vid in ru_chain_to_tids.get(c, []):
                    print(f"        {vid}  cov={ru[vid]['cov']}")

    # --- verdict ---
    print("\n" + "=" * 78)
    print("VERDICT")
    print("=" * 78)
    extraA = detA["kill_chains"] - st_in_ru
    extraA_tp = {c for c in extraA if c in ref_chains}
    n_extra_tp = len(extraA_tp)
    tp_cost = len(tpA)
    print(f"  contained-FIRST-exon: kill={len(killA)} FP={len(fpA)} "
          f"TP={len(tpA)} NET={netA:+d}")
    print(f"  TP over-kills beyond ST: {n_extra_tp}")
    safe = (netA >= 20) and (tp_cost <= 2) and (n_extra_tp <= 2)
    if safe:
        print(f"  --> SAFE to build: realizes net {netA:+d} at TP cost "
              f"{tp_cost} (<=~2), only {n_extra_tp} TP over-kill beyond ST.")
    else:
        reasons = []
        if netA < 20:
            reasons.append(f"net only {netA:+d} (< ~+25 target)")
        if tp_cost > 2:
            reasons.append(f"TP cost {tp_cost} (> ~2)")
        if n_extra_tp > 2:
            reasons.append(f"{n_extra_tp} TP over-kills ST does NOT make")
        print(f"  --> RISKY / reconsider: " + "; ".join(reasons))
    print("=" * 78)


if __name__ == "__main__":
    main()
