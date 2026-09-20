#!/usr/bin/env python3
"""psv_graph_combine.py — final corrected genome-wide result.

The genomic-span scan (psv_graph_genomewide.py + psv_graph_verify.py) is monotone-safe
for the 121 RESOLVABLE families: exon-union alignment only ever finds MORE distinguishing
columns, and adding columns cannot lower K, so a resolvable family stays resolvable (the
controls fam18/0/4 confirm this directly). Only the 24 FRONTIER candidates (18 genomic-K0
+ 6 align_fail) can move, so only those were re-run on exon-union copy models
(psv_graph_exonunion_custom.json). This script combines:
  - 121 resolvable families  -> kept from genomic-span (verify.json)
  - 24 frontier candidates   -> reclassified from exon-union
into the authoritative 145-family table, written in the figure schema.
Run: /home/juanfra/miniforge3/bin/python bench/psv_graph_combine.py
"""
import collections
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
verify = json.load(open(os.path.join(HERE, "psv_graph_verify.json")))
gw = {r["family"]: r for r in json.load(open(os.path.join(HERE, "psv_graph_genomewide.json")))["families"]}
exon = {r["family"]: r for r in json.load(open(os.path.join(HERE, "psv_graph_exonunion_custom.json")))["results"]}

final = []
moved = []
for f in verify["families_refined"]:
    fid, cls = f["family"], f["cls_refined"]
    if cls in ("fully_resolvable", "partial"):
        g = gw[fid]
        final.append(dict(family=fid, n_copies=f["n_copies"], K=f["K"], psvs=f["psvs"],
                          reads=f["reads"], cls_refined=cls, method="genomic_span",
                          assign=g["assign"], agree=g["agree"], utot=g["utot"]))
    else:  # frontier candidate -> exon-union reclassification
        e = exon[fid]
        ecls = e.get("cls", "still_failed")
        if ecls == "no_psv":
            ecls = "genuine_K0" if e["copies_aligned"] == e["n_copies"] else "align_fail"
        if ecls != cls:
            moved.append((fid, cls, ecls))
        final.append(dict(family=fid, n_copies=e["n_copies"], K=e["K"], psvs=e["psvs"],
                          reads=e["reads"], cls_refined=ecls, method="exon_union",
                          assign=e.get("assign", {}), agree=e.get("agree", 0), utot=e.get("utot", 0)))

cnt = collections.Counter(f["cls_refined"] for f in final)
N = len(final)
resolvable = cnt["fully_resolvable"] + cnt["partial"]
frontier = cnt["genuine_K0"]
tot_reads = sum(f["reads"] for f in final)
tot_one = sum(f["assign"].get("one_copy", 0) for f in final)
tot_collapsed = sum(f["assign"].get("collapsed_group", 0) for f in final)
tot_unassign = sum(f["assign"].get("unassignable_no_psv", 0) for f in final)
agree = sum(f["agree"] for f in final)
utot = sum(f["utot"] for f in final)

out = dict(
    method="exon_union_corrected", unique_families=N, by_class=dict(cnt),
    resolvable=resolvable, frontier=frontier,
    pct_resolvable=round(100 * resolvable / N, 1), pct_frontier=round(100 * frontier / N, 1),
    no_psv_triage={"genuine_K0": cnt["genuine_K0"], "align_fail": cnt["align_fail"]},
    total_reads=tot_reads, reads_to_one_copy=tot_one, reads_to_collapsed_group=tot_collapsed,
    reads_unassignable=tot_unassign, single_copy_validation=f"{agree}/{utot}",
    single_copy_validation_pct=round(100 * agree / max(utot, 1), 1),
    frontier_reclassified=moved,
    families_refined=[{k: f[k] for k in ("family", "n_copies", "K", "psvs", "reads", "cls_refined")}
                      for f in final],
)
json.dump(out, open(os.path.join(HERE, "psv_graph_exonunion_all.json"), "w"), indent=2)

print("=== FINAL corrected genome-wide result (exon-union) ===")
print(f"  unique families : {N}")
print(f"  by class        : {dict(cnt)}")
print(f"  RESOLVABLE : {resolvable} ({out['pct_resolvable']}%)   "
      f"genuine K0 : {frontier} ({out['pct_frontier']}%)   align_fail : {cnt['align_fail']}")
print(f"  reads {tot_reads}: {tot_one} one-copy, {tot_collapsed} collapsed-group, "
      f"{tot_unassign} no-PSV; single-copy val {agree}/{utot} ({out['single_copy_validation_pct']}%)")
print(f"\n  frontier candidates reclassified by exon-union ({len(moved)}):")
for fid, old, new in sorted(moved):
    print(f"    fam{fid:>3}: {old} -> {new}")
print(f"\n[+] wrote psv_graph_exonunion_all.json")
