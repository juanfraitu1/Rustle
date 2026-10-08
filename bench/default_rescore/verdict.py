#!/usr/bin/env python3
"""The decision of docs/PREREG_default_rescore_npip_2026-10-06.md, computed from the runner's products (nothing is judged by hand).

    verdict.py --work DIR

Reads DIR/{g0.json, g1.json, g3.json, support.json, support.copies.tsv, fs_DEF.json, fs_PRE.json} (written by run.sh) and the stored
reference values below; prints the gates, the four rules, the verdict (INVALID if a gate fails, else PASS iff R1-R4 hold, else FAIL) and the
reported-beside tables; writes DIR/verdict.json.
"""
import argparse
import csv
import json

# stored reference values (docs/archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md Amendment E table; docs/archive/2026-09/CONTAINER_HEADROOM_2026-09-30.md, out/scores.human.chr16.json)
REG_E = {"P": dict(own=23, found=10, locus=19, found_any=10), "GOOD": dict(own=21, found=6, locus=18, found_any=6),
         "ALL": dict(own=24, found=2, locus=16, found_any=2)}
D_STORED = {"u2": (0.588, 0.714, 0.645), "compara": (0.5, 1.0, 0.667)}
BARS = dict(R1_found=6, R2_locus=18, R3_u2_f=0.645)
KEYS = (("old_overlap_in_npip_nodes", "own"), ("tc_found_in_npip_nodes", "found"), ("locus_tc_found", "locus"), ("tc_found", "found_any"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True)
    w = ap.parse_args().work
    ld = lambda n: json.load(open(f"{w}/{n}"))
    g0, g1, g3, sup = ld("g0.json"), ld("g1.json"), ld("g3.json"), ld("support.json")
    fsd, fsp = ld("fs_DEF.json"), ld("fs_PRE.json")
    out = {"gates": {}, "rules": {}, "beside": {}}

    # G0: the HEAD scorer on the frozen read-pool arms reproduces the registered Amendment E numbers
    bad = [(arm, k, g0["arms"][arm].get(k), want[nm]) for arm, want in REG_E.items() for k, nm in KEYS if g0["arms"][arm].get(k) != want[nm]]
    out["gates"]["G0"] = not bad
    # G1: the own-node rule reproduces the stored pagedata flags on the frozen arms
    out["gates"]["G1"] = g1["mismatches"] == 0
    # G3: the HEAD family_score reproduces the stored default-arm scores on the stored clusters
    t = g3["truths"]
    g3ok = all((round(t[k]["sens"], 3), round(t[k]["prec"], 3), round(t[k]["f"], 3)) == D_STORED[k] for k in D_STORED)
    out["gates"]["G3"] = g3ok
    print("GATES  G0 scorer on frozen arms:", "PASS" if out["gates"]["G0"] else f"FAIL {bad}")
    print("       G1 own-node rule vs pagedata:", "PASS" if out["gates"]["G1"] else f"FAIL ({g1['mismatches']} mismatches)")
    print("       G3 family_score vs stored D:", "PASS" if g3ok else "FAIL " + str({k: (t[k]['sens'], t[k]['prec'], t[k]['f']) for k in D_STORED}))

    arms = sup["arms"]
    rows = list(csv.DictReader(open(f"{w}/support.copies.tsv"), delimiter="\t"))
    byname = {r["name"]: r for r in rows}
    D, P = arms["DEF"], arms["PRE"]
    u2 = fsd["truths"]["u2"]
    r1 = D["tc_found_in_npip_nodes"] >= BARS["R1_found"]
    r2 = D["locus_tc_found"] >= BARS["R2_locus"]
    r3 = round(u2["f"], 3) >= BARS["R3_u2_f"]
    own = {n: byname[n]["DEF_page_own_node"] == "1" for n in ("NPIPB2", "NPIPB6")}
    own_pre = {n: byname[n]["PRE_page_own_node"] == "1" for n in ("NPIPB2", "NPIPB6")}
    r4 = all(own.values())
    out["rules"] = dict(R1=dict(value=D["tc_found_in_npip_nodes"], bar=BARS["R1_found"], ok=r1), R2=dict(value=D["locus_tc_found"], bar=BARS["R2_locus"], ok=r2),
                        R3=dict(value=u2["f"], bar=BARS["R3_u2_f"], ok=r3), R4=dict(own_node_DEF=own, own_node_PRE=own_pre, ok=r4))
    print("RULES  R1 strict E-found within own nodes (DEF) %d of 25  >= %d: %s" % (D["tc_found_in_npip_nodes"], BARS["R1_found"], "ok" if r1 else "FAIL"))
    print("       R2 locus-level E-found (DEF) %d of 25  >= %d: %s" % (D["locus_tc_found"], BARS["R2_locus"], "ok" if r2 else "FAIL"))
    print("       R3 U2 bipartite F (DEF) %.3f  >= %.3f: %s   (sens %.3f prec %.3f; pairs tp %s of %s truth, %s predicted)" % (
        u2["f"], BARS["R3_u2_f"], "ok" if r3 else "FAIL", u2["sens"], u2["prec"], u2.get("pair_tp"), u2.get("pair_truth"), u2.get("pair_pred")))
    print("       R4 own node of NPIPB2, NPIPB6 in DEF: %s  (PRE: %s): %s" % (own, own_pre, "ok" if r4 else "FAIL"))
    if not all(out["gates"].values()):
        verdict = "INVALID (a gate failed: the instrument is not the registered one)"
    else:
        verdict = "PASS" if (r1 and r2 and r3 and r4) else "FAIL (" + ", ".join(n for n, ok in (("R1", r1), ("R2", r2), ("R3", r3), ("R4", r4)) if not ok) + ")"
    out["verdict"] = verdict
    print("VERDICT", verdict)

    print("\nBESIDE (no bar)")
    reg = REG_E["GOOD"]
    g2 = dict(own=P["old_overlap_in_npip_nodes"], found=P["tc_found_in_npip_nodes"], locus=P["locus_tc_found"], found_any=P["tc_found"])
    out["beside"]["G2_PRE_vs_registered_GOOD"] = dict(PRE=g2, registered=reg, equal=g2 == reg)
    print("  G2 diagnostic: PRE (HEAD, pre-flip) own/found/locus/found_any", g2, "registered GOOD", reg, "-> " + ("reproduced" if g2 == reg else "DRIFT"))
    for arm in ("PRE", "DEF"):
        a = arms[arm]
        print("  %s: own node %d | E-found in own nodes %d | locus-level %d | E-found chr16-wide %d | spliced-expressed (E) %s" % (
            arm, a["old_overlap_in_npip_nodes"], a["tc_found_in_npip_nodes"], a["locus_tc_found"], a["tc_found"], sup.get("tc_expressed")))
    flips = [(r["name"], r["PRE_tc_found"], r["DEF_tc_found"], r["PRE_page_own_node"], r["DEF_page_own_node"]) for r in rows
             if (r["PRE_tc_found"], r["PRE_page_own_node"]) != (r["DEF_tc_found"], r["DEF_page_own_node"])]
    out["beside"]["copy_changes_PRE_to_DEF"] = flips
    print("  copies that change between PRE and DEF (name, found PRE->DEF, own node PRE->DEF):")
    for n, fp, fd, op, od in flips:
        print(f"     {n:14s} found {fp}->{fd}  own node {op}->{od}")
    for arm, fs in (("PRE", fsp), ("DEF", fsd)):
        for k in ("u2", "compara", "soto"):
            r = fs["truths"][k]
            print("  %s %-7s sens %.3f prec %.3f F %.3f | pairs tp %s of %s truth, %s predicted | clusters %s" % (
                arm, k, r["sens"], r["prec"], r["f"], r.get("pair_tp"), r.get("pair_truth"), r.get("pair_pred"), r["clusters"]))
        cf = [x for x in fs["truths"]["compara"]["per_family"] if x["family_id"] == "CF153"]
        for x in cf:
            print("  %s CF153: %s of %s genes hit, cluster %s, sens %s prec %s F %s" % (arm, x["hit"], x["n_truth"], x["cluster"], x["sens"], x["prec"], x["f"]))
            out["beside"][f"CF153_{arm}"] = x
    json.dump(out, open(f"{w}/verdict.json", "w"), indent=1)


if __name__ == "__main__":
    main()
