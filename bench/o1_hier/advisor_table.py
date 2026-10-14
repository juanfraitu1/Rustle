#!/usr/bin/env python3
"""The advisor-facing O1 table: flat vs hierarchy-aware matching, where each truth family lands, and coverage.

    advisor_table.py OUT.tsv NAME=hier_match.json [NAME=hier_match.json ...]

One row per sample; also printed as a markdown table for the results section.
"""
import json
import sys

COLS = ["sample", "truth_families", "evaluable", "one_member", "no_member", "records_on_members",
        "records_on_unclustered_locus", "records_without_locus", "flat_exact", "flat_exact_rate", "hier_exact",
        "hier_exact_rate", "hier_exact_of_all", "flat_mean_jaccard", "hier_mean_jaccard", "flat_macro_recall",
        "hier_macro_recall", "flat_macro_precision", "hier_macro_precision", "hier_gene_recall", "lands_family",
        "lands_expansion", "lands_superfamily", "lands_sd_complex", "exact_any_expansion", "p_global", "p_within_family"]


def row(name, d):
    c, h, f = d["coverage"], d["hierarchy"], d["flat"]
    lv, nl = h["by_level"], d.get("null", {})
    r3 = lambda x: "NA" if x is None else round(x, 3)  # noqa: E731
    return {"sample": name, "truth_families": c["labels"], "evaluable": c.get("labels_evaluable", 0),
            "one_member": c.get("labels_one_member", 0), "no_member": c.get("labels_no_member", 0),
            "records_on_members": c.get("records_on_members", 0),
            "records_on_unclustered_locus": c.get("records_on_unclustered_locus", "NA"),
            "records_without_locus": c.get("records_without_locus", c.get("records_without_member", 0)),
            "flat_exact": f["exact"], "flat_exact_rate": r3(f["exact_rate"]), "hier_exact": h["exact"],
            "hier_exact_rate": r3(h["exact_rate"]), "hier_exact_of_all": r3(h["exact"] / c["labels"]),
            "flat_mean_jaccard": r3(f["mean_jaccard"]), "hier_mean_jaccard": r3(h["mean_jaccard"]),
            "flat_macro_recall": r3(f["macro_recall"]), "hier_macro_recall": r3(h["macro_recall"]),
            "flat_macro_precision": r3(f["macro_precision"]), "hier_macro_precision": r3(h["macro_precision"]),
            "hier_gene_recall": r3(h.get("macro_gene_recall")),
            "lands_family": lv.get("family", 0), "lands_expansion": lv.get("expansion", 0),
            "lands_superfamily": lv.get("superfamily", 0), "lands_sd_complex": lv.get("sd_complex", 0),
            "exact_any_expansion": d["exact_any"].get("expansion", 0),
            "p_global": r3(nl.get("global", {}).get("p_exact")),
            "p_within_family": r3(nl.get("within_family", {}).get("p_exact"))}


def main(argv=None):
    a = argv if argv is not None else sys.argv[1:]
    rows = [row(n, json.load(open(p))) for n, p in (x.split("=", 1) for x in a[1:])]
    with open(a[0], "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        fh.writelines("\t".join(str(r[k]) for k in COLS) + "\n" for r in rows)
    print("| " + " | ".join(COLS) + " |")
    print("|" + "---|" * len(COLS))
    for r in rows:
        print("| " + " | ".join(str(r[k]) for k in COLS) + " |")


if __name__ == "__main__":
    sys.exit(main())
