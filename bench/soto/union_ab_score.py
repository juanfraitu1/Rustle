#!/usr/bin/env python3
"""Score two catalog runs against Soto truth: does the exon-union substrate produce PURER families?

This is the end-to-end version of the offline experiment. The offline test built exon-sums inside Soto TRUTH
WINDOWS, so it knew where the loci were; the pipeline discovers them from reads. If the union's benefit came
from truth-window scoping rather than from the substrate, it disappears here.

Purity is scored by SEVERITY, not by counting impure families. A predicted family fusing 40 truth families and
one fusing 2 both count as "1 impure" — that metric is severity-blind, and using it earlier led me to call
this result falsified when it was not. The headline numbers are therefore:

  WORST     the largest number of truth families absorbed by any single predicted family
  FUSED     total truth families sitting inside a predicted family that spans more than one
  LOCI_IMP  truth members trapped inside such a family

Recall is reported beside them because the union trades loci for purity (measured ~12%): a purity gain paid
for by dropping every hard locus is not a gain, and the two must be read together.
"""
import argparse, csv
from collections import defaultdict, Counter


def truth(bed):
    out = []
    for ln in open(bed):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 4 or "|" not in f[3]:
            continue
        name, fam = f[3].split("|")[0], f[3].split("|")[1]
        out.append((fam, name, f[0], int(f[1]), int(f[2])))
    return out


def copies(path):
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r.get("chrom"):
                continue
            try:
                out.append((r["family_id"], r["chrom"], int(r["start"]), int(r["end"])))
            except (ValueError, KeyError):
                continue
    return out


def score(label, cps, tr):
    # map each predicted copy to every truth member it overlaps
    fam_to_truth = defaultdict(set)
    covered = set()
    for fid, c, s, e in cps:
        for tfam, tname, tc, ts, te in tr:
            if tc == c and s < te and e > ts:
                fam_to_truth[fid].add(tfam)
                covered.add((tfam, tname))
    truth_fams = {t[0] for t in tr}
    impure = {k: v for k, v in fam_to_truth.items() if len(v) > 1}
    worst = max((len(v) for v in fam_to_truth.values()), default=0)
    fused = len({f for v in impure.values() for f in v})
    loci_imp = sum(1 for tfam, tname in covered if any(tfam in v for v in impure.values()))
    pred_fams = {f for f, _, _, _ in cps}
    return {
        "label": label,
        "pred_families": len(pred_fams),
        "copies": len(cps),
        "members_covered": len(covered),
        "member_recall": 100.0 * len(covered) / max(len(tr), 1),
        "truth_families_hit": len({t for v in fam_to_truth.values() for t in v}),
        "impure_families": len(impure),
        "worst": worst,
        "fused": fused,
        "loci_imp": loci_imp,
        "pure_families": len(fam_to_truth) - len(impure),
        "_impure": impure,
        "_n_truth_fams": len(truth_fams),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True, help="copies.tsv, arm A (single-chain)")
    ap.add_argument("--b", required=True, help="copies.tsv, arm B (exon-union)")
    ap.add_argument("--bed", default="bench/soto/80_fams.chr.bed")
    ap.add_argument("--label-a", default="arm A")
    ap.add_argument("--label-b", default="arm B")
    a = ap.parse_args()

    tr = truth(a.bed)
    A = score(a.label_a, copies(a.a), tr)
    B = score(a.label_b, copies(a.b), tr)
    print(f"Soto truth: {len({t[0] for t in tr})} families, {len(tr)} members\n")

    rows = [
        ("predicted families", "pred_families", ""),
        ("copies emitted", "copies", ""),
        ("truth members covered", "members_covered", ""),
        ("member recall %", "member_recall", "recall — read WITH purity, not instead of it"),
        ("", "", ""),
        ("WORST: truth families in one predicted family", "worst", "lower is better"),
        ("FUSED: truth families inside impure predictions", "fused", "lower is better"),
        ("LOCI_IMP: members trapped in impure families", "loci_imp", "lower is better"),
        ("impure predicted families", "impure_families", "count only — severity-blind"),
        ("pure predicted families", "pure_families", "higher is better"),
    ]
    print(f"{'metric':<46}{a.label_a:>12}{a.label_b:>12}{'delta':>10}   note")
    print("-" * 100)
    for lab, k, note in rows:
        if not k:
            print()
            continue
        va, vb = A[k], B[k]
        if isinstance(va, float):
            d = f"{vb - va:+.1f}"
            va_s, vb_s = f"{va:.1f}", f"{vb:.1f}"
        else:
            d = f"{vb - va:+d}"
            va_s, vb_s = str(va), str(vb)
        print(f"{lab:<46}{va_s:>12}{vb_s:>12}{d:>10}   {note}")

    for arm in (A, B):
        if arm["_impure"]:
            big = sorted(arm["_impure"].items(), key=lambda kv: -len(kv[1]))[:3]
            print(f"\n{arm['label']}: largest impure predictions")
            for fid, fams in big:
                shown = sorted(fams)[:8]
                print(f"  {fid}: {len(fams)} truth families  {shown}{' ...' if len(fams) > 8 else ''}")

    print("\nVERDICT")
    if B["worst"] < A["worst"] and B["fused"] < A["fused"]:
        print(f"  {a.label_b} reduces both the worst fusion and the total families fused.")
    elif B["worst"] > A["worst"]:
        print(f"  {a.label_b} has a WORSE largest fusion than {a.label_a}.")
    d = B["members_covered"] - A["members_covered"]
    if d:
        verb = "gains" if d > 0 else "loses"
        print(f"  {a.label_b} {verb} {abs(d)} truth members vs {a.label_a} "
              f"({100*abs(d)/max(A['members_covered'],1):.0f}%).")


if __name__ == "__main__":
    main()
