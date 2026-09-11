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

Add --full-geneset (+ build --shared over that same full universe, see soto_replicate_from_sedef.py's own
--geneset) to ALSO admit non-family-eligible genes (lncRNA, processed_pseudogene, etc.) as MEMBERS of an
already-formed family -- Soto's own definition includes them ("SD98 genes associated with other gene
features... were also assigned a gene family ID [when they join one]"), and the advisor's own framing is
specifically "replicate Soto, whose definition includes pseudogenes/lncRNAs". Opt-in, off by default: the
--full-geneset omitted case is BYTE-IDENTICAL to before this flag existed (docs/o1_ledger.md §6ih).

This does NOT add those genes as full graph nodes able to bridge two families together -- that was tried
and rejected (§6ih: precision 0.925->0.730, the same promiscuous-bridge-gene failure as §6ie's own
mega-component bug). It attaches each extra gene to whichever already-formed family it shares an exon
with (ties broken by edge count), via soto_attach_noncoding_members.py's attach() -- one implementation,
imported here, not duplicated.
"""
import argparse, csv, os, sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from soto_attach_noncoding_members import attach  # noqa: E402

ELIGIBLE = {"protein_coding", "unprocessed_pseudogene",
            "transcribed_unprocessed_pseudogene", "translated_unprocessed_pseudogene"}


def mad_mean(vals):
    """Mean absolute deviation about the mean -- the paper's own METHODS-text wording ("mean absolute
    deviation")."""
    if len(vals) < 2:
        return 0.0
    m = sum(vals) / len(vals)
    return sum(abs(v - m) for v in vals) / len(vals)


def mad_median(vals):
    """Median absolute deviation about the median, UNSCALED (scale=1.0) -- what their actual released
    code computes (B_SD98_families.ipynb: `stats.median_abs_deviation(wssd_clust_median)`, no `scale=`
    override, so scipy's own default of 1.0 applies -- NOT the "mean absolute deviation" the paper's own
    prose describes). Confirmed via two independent reads of the notebook's raw source. Robust to the
    extreme-CN outliers this project has already documented in famCN distributions (e.g. BET1L=702),
    which a mean-based statistic is not.
    """
    if len(vals) < 2:
        return 0.0
    s = sorted(vals)
    n = len(s)
    med = s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2
    dev = sorted(abs(v - med) for v in vals)
    return dev[n // 2] if n % 2 else (dev[n // 2 - 1] + dev[n // 2]) / 2


def step5_step6(shared, genes, biotype, famcn, mad_threshold, out_path, mad_fn=mad_mean):
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
        if len(vals) < 2 or mad_fn([v for v, _ in vals]) < mad_threshold:
            final.append(comp)
            continue
        vals.sort()
        cur, groups = [], []
        for v, g in vals:
            if cur and mad_fn([x for x, _ in cur] + [v]) >= mad_threshold:
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
    ap.add_argument("--mad-statistic", choices=["mean", "median"], default="mean",
                     help="mean = the paper's own METHODS-text wording; median = what their released "
                          "code (B_SD98_families.ipynb) actually computes (scipy median_abs_deviation, "
                          "unscaled) -- default stays 'mean' so this flag is opt-in, not a silent change")
    ap.add_argument("--full-geneset",
                     help="OPT-IN: also admit genes in this (larger) geneset that are NOT in --geneset "
                          "as MEMBERS of an already-formed family, via a shared exon -- never as a way to "
                          "found a family or merge two together. --shared must have been built over this "
                          "SAME full geneset (soto_replicate_from_sedef.py --geneset <this file>), or the "
                          "extra genes will have no edges to attach through. Omit for the original, "
                          "eligible-only behaviour (byte-identical to before this flag existed).")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    mad_fn = mad_median if a.mad_statistic == "median" else mad_mean

    genes, biotype = set(), {}
    with open(a.geneset) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genes.add(r["gene_id"])
            biotype[r["gene_id"]] = r.get("biotype", "")

    # the BACKBONE clustering step must only ever see eligible-eligible edges -- an edge touching a
    # --full-geneset-only gene must NOT let step5_step6 treat that gene as a graph node (that is exactly
    # the naive, rejected approach: docs/o1_ledger.md §6ih measured it collapsing precision 0.925->0.730
    # by letting non-eligible genes bridge two components together). The unfiltered file is still used
    # as-is for the attach() call below, which needs the extra genes' own edges.
    shared = defaultdict(set)
    with open(a.shared) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ga, gb = r["gene_a"], r["gene_b"]
            if ga in genes and gb in genes:
                shared[ga].add(gb)
                shared[gb].add(ga)

    famcn = {}
    with open(a.famcn) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gid = r.get("gene_id") or r.get("Gene ID")
            v = r.get("famCN") or r.get("famCN_median") or r.get("Median famCN")
            try:
                famcn[gid] = float(v)
            except (TypeError, ValueError):
                pass

    n_fam, n_single, n_comps = step5_step6(shared, genes, biotype, famcn, a.mad, a.out, mad_fn=mad_fn)
    print(f"[step4-input] {len(genes)} SD98 genes, {len(shared)} genes with >=1 shared exon, "
          f"{n_comps} raw components", file=sys.stderr)
    print(f"[done] {n_fam} families, {n_single} singleton/ineligible genes -> {a.out}", file=sys.stderr)

    if not a.full_geneset:
        return

    # --full-geneset given: attach non-eligible members onto the backbone families just written, then
    # rewrite --out with the combined result (same 6-column shape, extra genes marked in `status`).
    gene_family = {}
    rows = []
    with open(a.out) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            rows.append(r)
            if r["family_id"]:
                gene_family[r["gene_id"]] = r["family_id"]

    with open(a.full_geneset) as fh:
        full_genes, full_biotype = set(), {}
        for r in csv.DictReader(fh, delimiter="\t"):
            full_genes.add(r["gene_id"])
            full_biotype[r["gene_id"]] = r.get("biotype", "")
    extra_genes = full_genes - genes
    attached = attach(gene_family, extra_genes, a.shared)

    fam_size = defaultdict(int)
    for r in rows:
        if r["family_id"]:
            fam_size[r["family_id"]] += 1
    for f in attached.values():
        fam_size[f] += 1

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "biotype", "family_id", "n_members", "famCN", "status"])
        for r in rows:
            n = fam_size[r["family_id"]] if r["family_id"] else int(r["n_members"])
            w.writerow([r["gene_id"], r["biotype"], r["family_id"], n, r["famCN"], r["status"]])
        for g in sorted(attached):
            f = attached[g]
            w.writerow([g, full_biotype.get(g, ""), f, fam_size[f],
                        f"{famcn[g]:.2f}" if g in famcn else "", "attached_noncoding_member"])
        for g in sorted(extra_genes - set(attached)):
            w.writerow([g, full_biotype.get(g, ""), "", 1,
                        f"{famcn[g]:.2f}" if g in famcn else "", "extra_gene_no_attachment"])

    print(f"[attach] {len(extra_genes)} extra genes from --full-geneset considered, "
          f"{len(attached)} attached to an existing family -> {a.out} (rewritten)", file=sys.stderr)


if __name__ == "__main__":
    main()
