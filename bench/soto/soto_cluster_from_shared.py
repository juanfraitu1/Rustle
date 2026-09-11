#!/usr/bin/env python3
"""Steps 5-6 of Soto's clustering (connected components -> famCN MAD split -> family/singleton call),
factored out so both soto_replicate_clustering.py (PAF-based, minimap2 map-back) and
soto_replicate_from_sedef.py (SEDEF-CIGAR-based, this replication round) share ONE implementation of the
famCN-split rule instead of two copies that could drift -- this file's step5_step6() is byte-identical
in behaviour to soto_replicate_clustering.py's inline steps 5-6 (same MAD formula, same biotype
eligibility set, same singleton bookkeeping), verified by diffing this script's output against the
original on the SAME shared-exon input (see the "verify" run in the ledger for that check).

Usage: soto_cluster_from_shared.py --shared shared_exons_from_sedef.tsv --geneset sd98_geneset_v1.tsv
                                    --famcn soto_famCN_S1C.tsv --out replicated_families_sedef.tsv
"""
import argparse, csv, sys
from collections import defaultdict

ELIGIBLE = {"protein_coding", "unprocessed_pseudogene",
            "transcribed_unprocessed_pseudogene", "translated_unprocessed_pseudogene"}


def mad(vals):
    if len(vals) < 2:
        return 0.0
    m = sum(vals) / len(vals)
    return sum(abs(v - m) for v in vals) / len(vals)


def step5_step6(shared, genes, biotype, famcn, mad_threshold, out_path):
    """shared: gene -> set(partner genes). genes: set of all SD98 gene ids (for singleton bookkeeping).
    biotype: gene -> biotype string. famcn: gene -> float famCN (genes absent are treated as un-splittable).
    Writes `out_path` in the same TSV shape soto_replicate_clustering.py emits.
    """
    seen, comps = set(), []
    for g in shared:
        if g in seen:
            continue
        stack, comp = [g], []
        seen.add(g)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in shared[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        comps.append(comp)

    final = []
    for comp in comps:
        vals = [(famcn[g], g) for g in comp if g in famcn]
        if len(vals) < 2 or mad([v for v, _ in vals]) < mad_threshold:
            final.append(comp)
            continue
        vals.sort()
        cur, groups = [], []
        for v, g in vals:
            if cur and mad([x for x, _ in cur] + [v]) >= mad_threshold:
                groups.append([g2 for _, g2 in cur])
                cur = []
            cur.append((v, g))
        if cur:
            groups.append([g2 for _, g2 in cur])
        nocn = [g for g in comp if g not in famcn]
        if groups and nocn:
            groups[0].extend(nocn)
        final.extend(groups)

    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "biotype", "family_id", "n_members", "famCN", "status"])
        fam_i = n_fam = n_single = 0
        for comp in sorted(final, key=lambda c: -len(c)):
            eligible = any(biotype.get(g) in ELIGIBLE for g in comp)
            if len(comp) >= 2 and eligible:
                fid = f"SEDEFFAM{fam_i}"
                fam_i += 1
                n_fam += 1
                status = "family"
            else:
                fid = ""
                status = "singleton" if len(comp) < 2 else "no_coding_member"
                n_single += len(comp)
            for g in sorted(comp):
                w.writerow([g, biotype.get(g, ""), fid, len(comp),
                            f"{famcn[g]:.2f}" if g in famcn else "", status])
        for g in sorted(genes - seen):
            w.writerow([g, biotype.get(g, ""), "", 1, f"{famcn[g]:.2f}" if g in famcn else "", "singleton"])
            n_single += 1
    return n_fam, n_single, len(comps)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shared", required=True, help="gene_a<TAB>gene_b TSV (soto_replicate_from_sedef.py output)")
    ap.add_argument("--geneset", required=True)
    ap.add_argument("--famcn", required=True)
    ap.add_argument("--mad", type=float, default=1.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    genes, biotype = set(), {}
    with open(a.geneset) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genes.add(r["gene_id"])
            biotype[r["gene_id"]] = r.get("biotype", "")

    shared = defaultdict(set)
    with open(a.shared) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            shared[r["gene_a"]].add(r["gene_b"])
            shared[r["gene_b"]].add(r["gene_a"])

    famcn = {}
    with open(a.famcn) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gid = r.get("gene_id") or r.get("Gene ID")
            v = r.get("famCN") or r.get("famCN_median") or r.get("Median famCN")
            try:
                famcn[gid] = float(v)
            except (TypeError, ValueError):
                pass

    n_fam, n_single, n_comps = step5_step6(shared, genes, biotype, famcn, a.mad, a.out)
    print(f"[step4-input] {len(genes)} SD98 genes, {len(shared)} genes with >=1 shared exon, "
          f"{n_comps} raw components", file=sys.stderr)
    print(f"[done] {n_fam} families, {n_single} singleton/ineligible genes -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
