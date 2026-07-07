#!/usr/bin/env python3
"""Evaluate the existing `complete_poa_core` baseline.

Compares the cross-chromosome read-conflict catalog WITHOUT completion
(gw_xcbase) against WITH completion (gw_comp2). Both outputs already exist in
`/home/juanfra/winloci_scratch/`.  The script maps every copy to a gene via
max-overlap against `annot_intervals.tsv`, then reports:

  * overall family/copy counts and how many copies completion added
  * per-family added-copy annotations
  * coverage of known multi-copy gene-family windows
  * coverage of diploid-CN oracle multi-copy genes

This is a lightweight coordinate-based evaluation that does NOT require the
shipped family_rna_refine DN-id namespace.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/eval_complete_core_baseline.py
"""
import csv
import os
import sys
from bisect import bisect_right
from collections import defaultdict

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

SCRATCH = "/home/juanfra/winloci_scratch"
BASE_FAM = os.path.join(SCRATCH, "gw_xcbase.families.tsv")
BASE_COPY = os.path.join(SCRATCH, "gw_xcbase.copies.tsv")
COMP_FAM = os.path.join(SCRATCH, "gw_comp2.families.tsv")
COMP_COPY = os.path.join(SCRATCH, "gw_comp2.copies.tsv")
ANNOT = os.path.join(SCRATCH, "annot_intervals.tsv")
ORACLE = os.path.join("bench", "diploid_cn_oracle.tsv")

KNOWN_WINDOWS = {
    # literature-known families: seed gene symbols -> expanded by annotation lookup
    "RABL2": {"seeds": {"RABL2A", "RABL2B"}, "prefix": "RABL2"},
    "APOBEC3": {"seeds": {"APOBEC3C", "APOBEC3D", "APOBEC3F"}, "prefix": "APOBEC3"},
    "MAGEA": {"seeds": {"MAGEA1", "MAGEA4", "MAGEA9", "MAGEA12"}, "prefix": "MAGEA"},
    "ANKRD18": {"seeds": {"ANKRD18A", "ANKRD18B"}, "prefix": "ANKRD18"},
    "RGPD8": {"seeds": {"RGPD8", "RANBP2"}, "prefix": "RGPD"},
    "ZNF92": {"seeds": {"ZNF92"}, "prefix": "ZNF92"},
    "GSTM": {"seeds": {"GSTM1", "GSTM2", "GSTM4", "GSTM5"}, "prefix": "GSTM"},
    "HERC2": {"seeds": {"HERC2"}, "prefix": "HERC2"},
}


def load_annot():
    by_c = defaultdict(list)
    with open(ANNOT) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            by_c[f[0]].append((int(f[1]), int(f[2]), f[3], f[4]))
    for c, ivals in by_c.items():
        ivals.sort()
    return by_c


def gene_of_factory(by_c):
    # build index: chrom -> (starts, sorted_ivals, max_len)
    idx = {}
    for c, ivals in by_c.items():
        starts = [x[0] for x in ivals]
        maxlen = max((e - s) for s, e, _, _ in ivals) if ivals else 0
        idx[c] = (starts, ivals, maxlen)

    def gene_of(chrom, qs, qe):
        rec = idx.get(chrom)
        if rec is None:
            return None, None, 0.0
        starts, ivals, maxlen = rec
        hi = bisect_right(starts, qe)
        best_ov = 0
        best = (None, None)
        best_span = 0
        i = hi - 1
        lo_bound = qs - maxlen
        while i >= 0:
            s, e, bt, g = ivals[i]
            if s < lo_bound:
                break
            if e > qs:
                ov = min(e, qe) - max(s, qs)
                if ov > best_ov:
                    best_ov = ov
                    best = (g, bt)
                    best_span = e - s
            i -= 1
        if best[0] is None:
            return None, None, 0.0
        dn_span = qe - qs
        frac_dn = (best_ov / dn_span) if dn_span > 0 else 0.0
        frac_gene = (best_ov / best_span) if best_span > 0 else 0.0
        return best[0], best[1], max(frac_dn, frac_gene)
    return gene_of


def load_catalog(copies_tsv):
    copies = []
    with open(copies_tsv) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            copies.append(dict(
                family_id=row["family_id"],
                tid=row["tid"],
                chrom=row["chrom"],
                start=int(row["start"]),
                end=int(row["end"]),
                n_exon=int(row["n_exon"]),
                strand=row["strand"],
                n_reads=int(row["n_reads"]),
            ))
    return copies


def families_from_copies(copies):
    fam = defaultdict(list)
    for c in copies:
        fam[c["family_id"]].append(c)
    return fam


def annotate(copies, gene_of):
    for c in copies:
        g, bt, ov = gene_of(c["chrom"], c["start"], c["end"])
        c["gene"] = g
        c["biotype"] = bt
        c["gene_ov"] = ov


def known_family_windows(annot, gene_of):
    """Build per-known-family (chrom, start, end, genes) windows from annotation."""
    windows = {}
    for name, spec in KNOWN_WINDOWS.items():
        seeds = spec["seeds"]
        prefix = spec["prefix"]
        genes = set()
        # collect all annotation intervals whose gene matches seeds or prefix
        for chrom, ivals in annot.items():
            for s, e, bt, g in ivals:
                if g in seeds or (prefix and g.startswith(prefix)):
                    genes.add((chrom, s, e, g))
        if not genes:
            windows[name] = None
            continue
        by_chrom = defaultdict(lambda: [1 << 60, -1, set()])
        for chrom, s, e, g in genes:
            by_chrom[chrom][0] = min(by_chrom[chrom][0], s)
            by_chrom[chrom][1] = max(by_chrom[chrom][1], e)
            by_chrom[chrom][2].add(g)
        windows[name] = [
            dict(chrom=c, start=lo, end=hi, genes=gs)
            for c, (lo, hi, gs) in by_chrom.items()
        ]
    return windows


def copies_in_window(copies, windows, pad=0):
    """Return copies whose span overlaps any known-family window (with padding)."""
    out = []
    for c in copies:
        for w in windows or []:
            if c["chrom"] == w["chrom"]:
                if c["end"] >= w["start"] - pad and c["start"] <= w["end"] + pad:
                    out.append(c)
                    break
    return out


def load_oracle_genes():
    genes = []
    with open(ORACLE) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            fam, gene, cls = f[0], f[1], f[12] if len(f) > 12 else ""
            if gene and gene != "NA" and cls == "multi_copy":
                genes.append(gene)
    return set(genes)


def main():
    print("[load] annotation ...")
    annot = load_annot()
    gene_of = gene_of_factory(annot)

    print("[load] baseline and complete-core catalogs ...")
    base_copies = load_catalog(BASE_COPY)
    comp_copies = load_catalog(COMP_COPY)
    annotate(base_copies, gene_of)
    annotate(comp_copies, gene_of)

    base_fam = families_from_copies(base_copies)
    comp_fam = families_from_copies(comp_copies)

    # ---- overall counts ----
    print("\n================ OVERALL COUNTS ================")
    print(f"Baseline families: {len(base_fam)}  copies: {len(base_copies)}")
    print(f"Complete-core families: {len(comp_fam)}  copies: {len(comp_copies)}")

    # Family IDs are re-sorted when copies are added, so matching by family_id string
    # is unreliable. Identify added copies by global tid set difference instead.
    base_tid_set = {c["tid"] for c in base_copies}
    added = [c for c in comp_copies if c["tid"] not in base_tid_set]
    print(f"Added copies: {len(added)}")

    # group added copies by their complete-catalog family
    gains = defaultdict(list)
    for c in added:
        gains[c["family_id"]].append(c)
    gains = sorted(gains.items(), key=lambda x: -len(x[1]))
    print(f"Families with added copies: {len(gains)}")

    # ---- per-added-copy report ----
    print("\n================ ADDED COPY ANNOTATIONS ================")
    print("family_id\tchrom\tstart\tend\tstrand\tn_reads\tn_exon\tgene\tbiotype\tgene_ov")
    all_added = []
    for fid, added in sorted(gains, key=lambda x: -len(x[1])):
        for c in added:
            all_added.append(c)
            print(f"{fid}\t{c['chrom']}\t{c['start']}\t{c['end']}\t{c['strand']}\t"
                  f"{c['n_reads']}\t{c['n_exon']}\t{c['gene'] or 'NA'}\t{c['biotype'] or 'NA'}\t"
                  f"{c['gene_ov']:.3f}")

    # biotype breakdown
    bt = defaultdict(int)
    for c in all_added:
        bt[c["biotype"] or "unannotated"] += 1
    print("\nAdded-copy biotype breakdown:")
    for k, v in sorted(bt.items(), key=lambda x: -x[1]):
        print(f"  {k}: {v}")

    # ---- known-family window coverage ----
    print("\n================ KNOWN FAMILY WINDOW COVERAGE ================")
    windows = known_family_windows(annot, gene_of)
    for name, wins in windows.items():
        if wins is None:
            print(f"{name}: no annotated genes found")
            continue
        base_in = copies_in_window(base_copies, wins, pad=1_500_000)
        comp_in = copies_in_window(comp_copies, wins, pad=1_500_000)
        base_genes = {c["gene"] for c in base_in if c["gene"]}
        comp_genes = {c["gene"] for c in comp_in if c["gene"]}
        added_genes = comp_genes - base_genes
        print(f"{name}: windows={len(wins)}  baseline copies={len(base_in)}  "
              f"complete copies={len(comp_in)}  (+{len(comp_in)-len(base_in)})  "
              f"new genes covered: {sorted(added_genes) or 'none'}")

    # ---- diploid-CN oracle multi-copy gene coverage ----
    print("\n================ DIPLOID-CN ORACLE MULTI-COPY GENES ================")
    oracle_genes = load_oracle_genes()
    # map gene symbol -> coords via annotation
    gene_coords = {}
    for chrom, ivals in annot.items():
        for s, e, bt, g in ivals:
            if g in oracle_genes:
                gene_coords[g] = (chrom, s, e)

    def covers(copies, gene):
        chrom, s, e = gene_coords[gene]
        for c in copies:
            if c["chrom"] == chrom and c["end"] >= s and c["start"] <= e:
                return True
        return False

    covered_base = [g for g in sorted(gene_coords) if covers(base_copies, g)]
    covered_comp = [g for g in sorted(gene_coords) if covers(comp_copies, g)]
    newly_covered = sorted(set(covered_comp) - set(covered_base))
    print(f"Oracle multi-copy genes with coordinates: {len(gene_coords)}")
    print(f"Covered baseline: {len(covered_base)}  complete: {len(covered_comp)}  "
          f"newly covered: {len(newly_covered)}")
    if newly_covered:
        print("Newly covered genes:", newly_covered)

    # ---- gene-annotation summary of added copies ----
    print("\n================ ADDED-COPY GENE SET ================")
    added_genes = sorted({c["gene"] for c in all_added if c["gene"]})
    print(f"Distinct gene symbols among added copies: {len(added_genes)}")
    print(", ".join(added_genes[:50]) + (" ..." if len(added_genes) > 50 else ""))


if __name__ == "__main__":
    main()
