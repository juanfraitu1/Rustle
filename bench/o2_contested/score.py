#!/usr/bin/env python3
"""r1509 per-class readout. Joins haplotype_truth's per_read.tsv (one row per (read, family) with a near-perfect
haplotype hit) to the arm's assignments table and the extended as_table (splice shapes of the genome-wide primary and
its best-AS rival). Classes (margins = the assignment row's region-local as_margin / ms_margin):
  wall      AS margin 0 and ms margin 0, same splice shape
  both      |AS margin| <= 10 and ms margin <= 10, same shape (the near tie on both scores)
  shape     primary and rival differ in splice shape (spliced vs unspliced)
  other     the rest (AS margin beyond 10 with ms within 10, or no ms column = BASE arm)
Per class and arm: rows, decidable (truth = one locus), assigned-right, assigned-wrong, abstained, outside the
catalog; the aligner's primary right on the decidable rows; the 10-AS rule (as_margin > 10 picks the primary) right.
Abstained / outside rows come from haplotype_truth's 3,000-read sample plus every assigned read, so they are counts on
that sample, not on the arm. Usage: score.py ARM_PREFIX TRUTH_DIR MOLECULES.tsv"""
import collections, csv, sys

def main(prefix, truth_dir, molecules):
    shape = {}
    with open(molecules) as fh:
        for line in fh:
            if line[0] == "#": continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 15 and f[11] != "NA": shape[f[0]] = (f[10], f[13])
    rows = {}
    with open(prefix + ".assignments.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows[(r["read_name"], r["family_id"])] = r
    def cls(r):
        a = int(r["as_margin"]) if r["as_margin"] not in ("", "NA") else None
        m = int(r["ms_margin"]) if r.get("ms_margin", "NA") not in ("", "NA", None) else None
        sh = shape.get(r["read_name"]); same = sh is not None and sh[0] == sh[1]
        if sh is not None and not same: return "shape"
        if a == 0 and m == 0 and same: return "wall"
        if a is not None and m is not None and abs(a) <= 10 and m <= 10 and same: return "both"
        return "other"
    c = collections.defaultdict(collections.Counter)
    with open(f"{truth_dir}/per_read.tsv") as fh:
        for t in csv.DictReader(fh, delimiter="\t"):
            r = rows.get((t["read"], t["family"]))
            if r is None: continue
            k = cls(r); c[k]["rows"] += 1
            if t["truth"] == "": c[k]["outside"] += 1; continue
            c[k]["decidable"] += 1
            if t["status"] == "assigned": c[k]["assigned_right" if t["o2_right"] == "1" else "assigned_wrong"] += 1
            else: c[k]["abstained"] += 1
            c[k]["aligner_right"] += int(t["aligner_right"])
            if int(r["as_margin"]) > 10: c[k]["rule10_fires"] += 1; c[k]["rule10_right"] += int(t["aligner_right"])
    print("class\trows\tdecidable\tassigned_right\tassigned_wrong\tabstained\toutside\to2_wrong_rate\taligner_wrong_rate\trule10_fires\trule10_wrong_rate")
    for k in ("both", "shape", "wall", "other"):
        x = c[k]; asg = x["assigned_right"] + x["assigned_wrong"]; d = x["decidable"]
        print("\t".join(map(str, [k, x["rows"], d, x["assigned_right"], x["assigned_wrong"], x["abstained"], x["outside"],
              f"{x['assigned_wrong']/asg:.3f}" if asg else "NA", f"{1-x['aligner_right']/d:.3f}" if d else "NA",
              x["rule10_fires"], f"{1-x['rule10_right']/x['rule10_fires']:.3f}" if x["rule10_fires"] else "NA"])))

if __name__ == "__main__": main(*sys.argv[1:4])
