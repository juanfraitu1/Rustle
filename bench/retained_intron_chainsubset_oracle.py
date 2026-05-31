#!/usr/bin/env python3
"""retained_intron_chainsubset_oracle.py

CHEAP analysis-only oracle that bounds the OVER-KILL risk of the FAITHFUL
StringTie `retained_intron` predicate in Rustle, BEFORE any code is written. No
production code is touched and no behavior changes.

WHY THIS, NOT THE GEOMETRY ORACLE
---------------------------------
The prior oracle (bench/retained_intron_geometry_oracle.py, ce8a0c0) proved the
"exon-contained-in-low-intron" GEOMETRY over-kills catastrophically
(contained-first-exon: net -114, 133 TP killed; any-exon: net -556) and catches
only 5/27 of ST's actual retained_intron kills. It concluded ST's true mechanism
is a CHAIN-SUBSET / refinement relationship plus high killer-coverage dominance:

  Victim V is killed by killer K (same strand, overlapping) when V splices the
  SAME WAY as K everywhere except V MERGES ACROSS exactly one of K's introns --
  i.e. V's intron chain == K's intron chain with exactly ONE (low-coverage) K
  intron REMOVED ("retained" as exon in V) -- AND K dominates V in coverage
  (ST median killer/victim cov ratio ~= 114.8).

THE FAITHFUL PREDICATE (what this tool replicates on RU predictions)
--------------------------------------------------------------------
For each ORDERED pair (K killer, V victim), same strand, overlapping, with
cov(K)/cov(V) >= DOMINANCE (sweep DOMINANCE in {2,5,20,50,100}):

  STRICT (V chain = K chain minus one LOW intron):
    There is exactly one index r such that
       V.chain == K.chain[:r] + K.chain[r+1:]
    AND K's intron r is flagged LOW in K's mask.
    (V is K with intron r retained as exon. All other introns identical.)

  LOOSE (V retains >=1 low K intron; all V introns are K introns):
    Every V intron is also a K intron, and the set of K introns absent from V
    is non-empty and ALL of them are flagged LOW in K. (Allows >1 retained
    intron, but every retention must be a low K intron.)

A fired pair => V in KILL SET (killed by K). KILL SET classified FP (chain not
in ref) vs TP (in ref). NET = FP - TP, per DOMINANCE level. Sanity: overlap with
ST's 27-FP set and extra-beyond-ST over-kill (FP vs TP).

REUSED MACHINERY (identical to the geometry oracle)
---------------------------------------------------
GTF parse, intron-token convention (exon[i-1].end+1, exon[i].start-1), the
RU pred_intron_low span+strand+nintron+cov matching, and the ST
retained_intron chain loader. Stdlib only.
"""

import json
import collections

RU_GTF   = "/tmp/ri_ru.gtf"
RU_JSONL = "/tmp/ri_ru.jsonl"
ST_JSONL = "/tmp/ri_st.jsonl"
REF_GTF  = "/mnt/c/Users/jfris/Desktop/GGO_19.gtf"

DOMINANCE_LEVELS = [2, 5, 20, 50, 100]


# ---------------------------------------------------------------------------
# GTF parsing  (identical to geometry oracle)
# ---------------------------------------------------------------------------
def _attr(attrs, key):
    needle = key + ' "'
    if needle not in attrs:
        return None
    return attrs.split(needle, 1)[1].split('"', 1)[0]


def load_gtf(path):
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
        else:
            rec["exons"].append((int(f[3]), int(f[4])))
            ec = _attr(f[8], "cov")
            rec["excov"].append(float(ec) if ec is not None else 0.0)
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
    return (strand, tuple(intron_tokens(exons)))


# ---------------------------------------------------------------------------
# Low-intron mask: match RU pred_intron_low -> GTF tx  (identical to geometry)
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
            "start": d["start"] + 1,
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
            best = min(cands, key=lambda e: abs((e["cov"] or 0) - (rec["cov"] or 0)))
            rec["low"] = list(best["mask"])
            rec["low_src"] = "mask"
            stats["mask_used"] += 1
        else:
            stats["zero"] += 1
            stats["no_mask"] += 1
    return stats


# ---------------------------------------------------------------------------
# ST sanity set  (identical to geometry oracle)
# ---------------------------------------------------------------------------
def load_st_retained_intron_chains(path):
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
            continue
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
# FAITHFUL chain-subset kill computation
# ---------------------------------------------------------------------------
def _strict_retained(kchain, vchain, klow):
    """Return the removed-intron index r if V.chain == K.chain minus exactly one
    LOW K intron, else None.

    Strict definition: there exists exactly one r in [0,len(kchain)) such that
    kchain[:r]+kchain[r+1:] == vchain, and klow[r] is True.
    Requires len(vchain) == len(kchain) - 1.
    """
    if len(vchain) != len(kchain) - 1:
        return None
    # find the first index where they diverge
    r = None
    nk = len(kchain)
    i = 0  # index into kchain
    j = 0  # index into vchain
    skipped = None
    while i < nk:
        if j < len(vchain) and kchain[i] == vchain[j]:
            i += 1
            j += 1
        else:
            if skipped is not None:
                return None  # more than one mismatch
            skipped = i
            i += 1  # skip this K intron, do NOT advance j
    if skipped is None:
        return None
    if j != len(vchain):
        return None
    if skipped < len(klow) and klow[skipped]:
        return skipped
    return None


def _loose_retained(kset, vchain, kchain, klow_set):
    """LOOSE: every V intron is a K intron, the K introns absent from V are
    non-empty, and ALL of them are flagged LOW. klow_set = set of LOW K introns.
    Returns the list of retained (removed) introns if it fires, else None.
    """
    vset = set(vchain)
    if not vset.issubset(kset):
        return None
    removed = kset - vset
    if not removed:
        return None
    # all removed must be low
    if not removed.issubset(klow_set):
        return None
    return removed


def compute_kill(tx, ref_chains, dominance, mode):
    """Return (kill_set_tids, fp_chains, tp_chains, detail).
    mode in {'strict','loose'}.
    Fires K->V when cov(K)/cov(V) >= dominance, same strand, overlapping, and the
    chain-subset/retained-low-intron test passes.
    """
    buckets = collections.defaultdict(list)
    for tid, rec in tx.items():
        if len(rec["exons"]) >= 2:
            buckets[(rec["chrom"], rec["strand"])].append(tid)

    kill = set()
    killed_by = {}
    ratios = []
    for key, tids in buckets.items():
        tids.sort(key=lambda t: tx[t]["exons"][0][0])
        # precompute chains + low for each
        chain = {}
        low = {}
        cset = {}
        clowset = {}
        for t in tids:
            ch = tuple(intron_tokens(tx[t]["exons"]))
            chain[t] = ch
            lw = tx[t].get("low") or []
            low[t] = lw
            cset[t] = set(ch)
            clowset[t] = {ch[i] for i in range(len(ch)) if i < len(lw) and lw[i]}
        n = len(tids)
        for i in range(n):
            kid = tids[i]
            K = tx[kid]
            kcov = K["cov"] if K["cov"] is not None else 0.0
            if kcov <= 0:
                continue
            kstart, kend = K["exons"][0][0], K["exons"][-1][1]
            kchain = chain[kid]
            klow = low[kid]
            # killer must have at least one LOW intron to retain
            if not clowset[kid]:
                continue
            for j in range(n):
                if j == i:
                    continue
                vid = tids[j]
                V = tx[vid]
                vstart, vend = V["exons"][0][0], V["exons"][-1][1]
                if vend < kstart or vstart > kend:
                    continue
                vcov = V["cov"] if V["cov"] is not None else 0.0
                if vcov <= 0:
                    # treat as infinitely dominated (cov 0 victim); allow kill
                    ratio = float("inf")
                else:
                    ratio = kcov / vcov
                if ratio < dominance:
                    continue
                vchain = chain[vid]
                fired = False
                if mode == "strict":
                    r = _strict_retained(kchain, vchain, klow)
                    if r is not None:
                        fired = True
                        retained = (kchain[r],)
                else:  # loose
                    removed = _loose_retained(cset[kid], vchain, kchain,
                                              clowset[kid])
                    if removed is not None:
                        fired = True
                        retained = tuple(sorted(removed))
                if fired:
                    if vid not in kill:
                        kill.add(vid)
                        killed_by[vid] = (kid, retained, ratio)
                        ratios.append(ratio if ratio != float("inf") else None)

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
    detail = {"killed_by": killed_by, "kill_chains": kill_chains, "ratios": ratios}
    return kill, fp_chains, tp_chains, detail


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 78)
    print("retained_intron CHAIN-SUBSET + DOMINANCE oracle  --  faithful-predicate bound")
    print("=" * 78)

    ru = load_gtf(RU_GTF)
    ref = load_gtf(REF_GTF)
    ru_multi = sum(1 for r in ru.values() if len(r["exons"]) > 1)
    print(f"\nRU final predictions: {len(ru)} transcripts ({ru_multi} multi-exon)")
    print(f"Reference (GGO_19) transcripts: {len(ref)}")

    n_cov = sum(1 for r in ru.values() if r["cov"] is not None)
    print(f"[RU coverage source] transcript cov attr present for {n_cov}/{len(ru)} tx.")

    ref_chains = set()
    for r in ref.values():
        if len(r["exons"]) >= 2:
            ref_chains.add(chain_key(r["strand"], r["exons"]))
    print(f"Reference multi-exon intron-chains: {len(ref_chains)}")

    evs = load_pred_intron_low(RU_JSONL)
    print(f"\n[Low-intron mask] RU pred_intron_low events: {len(evs)} "
          "(matched to GTF by span+strand+nintron+cov).")
    stats = build_low_intron_map(ru, evs)
    print(f"  unique mask match: {stats['unique']}  ambiguous(cov-broken): "
          f"{stats['multi']}  no match: {stats['zero']}  (of {ru_multi} multi-exon)")
    src_counts = collections.Counter(r.get("low_src", "none")
                                     for r in ru.values() if len(r["exons"]) > 1)
    print(f"  low-ness source per multi-exon tx: {dict(src_counts)}  "
          "(NOTE: this oracle does NOT use the geometry fallback; mask-only).")

    # ST sanity prize set
    st_chains, st_nev = load_st_retained_intron_chains(ST_JSONL)
    ru_chain_to_tids = collections.defaultdict(list)
    for tid, rec in ru.items():
        if len(rec["exons"]) >= 2:
            ru_chain_to_tids[chain_key(rec["strand"], rec["exons"])].append(tid)
    ru_final_chains = set(ru_chain_to_tids.keys())
    st_in_ru = st_chains & ru_final_chains
    st_killed_ru_fp = {c for c in st_in_ru if c not in ref_chains}
    st_killed_ru_tp = {c for c in st_in_ru if c in ref_chains}
    print(f"\n[ST sanity prize] ST pred_kill(retained_intron) chains: {len(st_chains)}; "
          f"present in RU final: {len(st_in_ru)} "
          f"(RU-FP={len(st_killed_ru_fp)}, RU-TP={len(st_killed_ru_tp)}); "
          f"narrative net = +{len(st_killed_ru_fp)-len(st_killed_ru_tp)}.")

    # ---------------------------------------------------------------
    # DOMINANCE sweep, both modes
    # ---------------------------------------------------------------
    for mode in ("strict", "loose"):
        title = ("STRICT (V chain = K chain minus exactly one LOW intron)"
                 if mode == "strict"
                 else "LOOSE (V retains >=1 low K intron; all V introns in K)")
        print("\n" + "=" * 78)
        print(f"MODE: {title}")
        print("=" * 78)
        hdr = (f"{'DOM':>5} | {'kill':>5} {'FP':>4} {'TP':>4} {'NET':>5} | "
               f"{'capST27':>8} {'%ST':>4} | {'xkFP':>5} {'xkTP':>5}")
        print(hdr)
        print("-" * len(hdr))
        best = None
        for dom in DOMINANCE_LEVELS:
            kill, fp, tp, det = compute_kill(ru, ref_chains, dom, mode)
            net = len(fp) - len(tp)
            kc = det["kill_chains"]
            inter = kc & st_killed_ru_fp
            recall = (100.0 * len(inter) / len(st_killed_ru_fp)
                      if st_killed_ru_fp else 0.0)
            extra = kc - st_in_ru
            xk_fp = sum(1 for c in extra if c not in ref_chains)
            xk_tp = sum(1 for c in extra if c in ref_chains)
            print(f"{dom:>5} | {len(kill):>5} {len(fp):>4} {len(tp):>4} {net:>+5} | "
                  f"{len(inter):>4}/{len(st_killed_ru_fp):<3} {recall:>3.0f}% | "
                  f"{xk_fp:>5} {xk_tp:>5}")
            cand = {"dom": dom, "kill": len(kill), "fp": len(fp), "tp": len(tp),
                    "net": net, "inter": len(inter), "xk_tp": xk_tp,
                    "xk_fp": xk_fp, "det": det, "tp_list": tp}
            # best = max NET subject to TP cost <= 2
            if len(tp) <= 2:
                if best is None or net > best["net"]:
                    best = cand
        # report best for this mode
        print("-" * len(hdr))
        if best is not None:
            print(f"  best @TP<=2: DOM={best['dom']}  NET={best['net']:+d}  "
                  f"(FP={best['fp']}, TP={best['tp']}, capST27={best['inter']}/"
                  f"{len(st_killed_ru_fp)}, over-kill TP={best['xk_tp']})")
            if best["tp_list"]:
                print("    TP victims lost at best DOM:")
                for vid, ck in sorted(best["tp_list"]):
                    kid, retained, ratio = best["det"]["killed_by"][vid]
                    rs = "inf" if ratio == float("inf") else f"{ratio:.1f}"
                    print(f"      {vid} cov={ru[vid]['cov']} killed_by={kid} "
                          f"ratio={rs}")
        else:
            print("  no DOM level achieves TP cost <= 2.")
        # stash strict-best for the verdict
        if mode == "strict":
            strict_best = best

    # ---------------------------------------------------------------
    # VERDICT
    # ---------------------------------------------------------------
    print("\n" + "=" * 78)
    print("VERDICT")
    print("=" * 78)
    b = strict_best
    if b is None:
        print("  STRICT predicate: no DOM level reaches TP cost <= 2.")
        print("  --> SHELVE: cannot realize the prize without TP loss.")
    else:
        net = b["net"]
        tpc = b["tp"]
        xk = b["xk_tp"]
        cap = b["inter"]
        print(f"  STRICT best: DOM={b['dom']} NET={net:+d} FP={b['fp']} TP={tpc} "
              f"capST27={cap}/{len(st_killed_ru_fp)} over-kill-TP={xk}")
        # SAFE if clearly net-positive (>=+10), TP<=2, and reproduces most of ST 27
        safe = (net >= 10) and (tpc <= 2) and \
               (cap >= 0.6 * len(st_killed_ru_fp))
        if safe:
            print(f"  --> SAFE to build. The faithful chain-subset+dominance "
                  f"predicate realizes net {net:+d} at TP cost {tpc}, "
                  f"reproducing {cap}/{len(st_killed_ru_fp)} of ST's FP kills.")
            print(f"      retainedintron_like rule: kill V when same-strand "
                  f"overlapping K has chain == V.chain + one LOW intron AND "
                  f"cov(K)/cov(V) >= {b['dom']}.")
        else:
            reasons = []
            if net < 10:
                reasons.append(f"net only {net:+d} (< +10)")
            if tpc > 2:
                reasons.append(f"TP cost {tpc} (> 2)")
            if cap < 0.6 * len(st_killed_ru_fp):
                reasons.append(f"reproduces only {cap}/{len(st_killed_ru_fp)} "
                               "of ST's 27 FP kills (< 60%)")
            print("  --> SHELVE / not cleanly realizable: " + "; ".join(reasons))
    print("=" * 78)


if __name__ == "__main__":
    main()
