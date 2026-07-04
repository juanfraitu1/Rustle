#!/usr/bin/env python3
"""rna_copy_number_depth.py -- the THIRD O1 copy-number estimator: read-DEPTH.

Completes the RNA copy-number estimator stack for the TRULY UNRESOLVABLE (K=0,
exonically-identical, collapsed) copies -- the O2 assignment floor. Those copies share
one hap-vector, so they raise NO read-conflict edge; the chromatic number chi(H) cannot
count them (it collapses the group to ONE color) and family_copy_number.py flags the
shortfall as `collapsed_excess = sum(size-1)` per collapse group, noting "need depth/DNA".
THIS script builds the depth signal that fills that gap, from RNA ALONE.

  Estimator stack (each counts a STRICTLY DIFFERENT set of copies):
    (1) n_loci   -- reference-locus count. Counts every co-located reference copy fed to
                    the graph (Tier1+Tier2+Tier3). MISSES copies the *reference itself*
                    collapses into one locus (reference-collapsed / SDA-Vollger regime).
    (2) chi(H)   -- conflict chromatic number. Counts every DISTINGUISHABLE copy (a private
                    hap-vector => its own color), INCLUDING unassignable Tier-2 copies.
                    MISSES truly-identical copies (Tier-3 collapse) = `collapsed_excess`.
    (3) depth-CN -- THIS estimator. Counts EXPRESSED copies, identical ones INCLUDED, by
                    aggregate expression / single-copy expression. MISSES silent copies
                    (expressed ~0) and is confounded by per-copy expression variation.
                    => depth-CN is a NOISY LOWER BOUND on genomic CN, not an unbiased count.

  depth-CN(family) = E_fam / lambda,  E_fam = total reads over the family's reference copies.
    lambda = single-copy expression unit. Two anchors, for the two regimes:
      lambda_within  = median per-reference-copy read count WITHIN the family. Valid on the
                       REFERENCE-RESOLVED substrate (GGO): even Tier-3 hap-identical copies
                       sit at DISTINCT reference loci, each with its own read pile, so a
                       per-copy expression sample exists. depth-CN_within then RECOVERS
                       ~n_loci when expression is uniform, and EXCEEDS it exactly when some
                       reference locus carries anomalous depth (a *within-reference* collapse
                       hiding several genomic copies at one position). It CANNOT reveal a copy
                       the reference resolved into its own zero-extra-depth locus -- so on GGO
                       it is mostly a cross-check on n_loci (honest: depth adds little here).
      lambda_global  = genome-wide single-copy RNA expression floor = median n_reads over
                       single-copy (non-family) de-novo transcripts. This is the anchor for
                       the REFERENCE-COLLAPSED regime, where the whole family piles onto ONE
                       reference locus (no sibling loci to form lambda_within) and depth is
                       the ONLY multiplicity signal (parCN's genome-wide single-copy-depth
                       normalization, Sudmant 2010 / QuicK-mer2, ported to RNA).

  HONESTY (load-bearing). RNA depth counts EXPRESSED copies, NOT genomic copies:
    (C1) per-copy expression VARIES -- one highly-expressed copy mimics several (over-count
         per locus) while co-regulated copies mimic one; we report the within-family
         expression CV and the dominant-copy share as the confound magnitude.
    (C2) SILENT copies (expressed ~0) are INVISIBLE to RNA -> systematic UNDER-count.
    (C3) => depth-CN is a noisy LOWER bound on genomic CN. The DNA per-haplotype CN
         (diploid_cn_oracle.tsv asm_hapCN = max(mat,pat) contig copies; Verkko/hifiasm
         trio) is the genomic-CN TRUTH and is used for VALIDATION ONLY (correlation +
         under-count quantification). DNA NEVER gates; RNA is primary.

  Bound structure this estimator lives in (family_copy_number.py / THEORY.md, unchanged):
      n_resolvable <= chi(H) <= n_loci = true_copy_lower_bound <= true_copy_number
      general form (both regimes):  max(n_loci, chi(H)) <= true_copy_number
    depth-CN is an INDEPENDENT lower-bound probe of the same true_copy_number, agreeing with
    n_loci where the reference is resolved and reaching PAST it (toward asm_hapCN) where the
    reference is collapsed -- but capped below asm_hapCN by silent copies (C2).

Inputs (read-only; does NOT modify family_copy_number.* shipped output):
  bench/family_copy_number.tsv          n_loci, chi_H, collapsed_excess, n_resolvable, gene
  bench/sun_identifiability.tsv         per-copy tier (1/2/3), chrom/start/end
  /home/juanfra/winloci_scratch/validated_families.tsv   per-locus nreads (deduped via pg)
  /home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv   single-copy expression anchor
  bench/diploid_cn_oracle.tsv           asm_hapCN  (DNA truth; VALIDATION ONLY)
Outputs: bench/rna_copy_number_depth.tsv (per family) + .json (summary + confounds + valid.)
Run from repo root: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/rna_copy_number_depth.py
"""
import collections
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import psv_graph_genomewide as pg  # noqa: E402  (dedup_copies, FAM_TSV)

FCN_TSV = os.path.join(HERE, "family_copy_number.tsv")
SUN_TSV = os.path.join(HERE, "sun_identifiability.tsv")
ORACLE = os.path.join(HERE, "diploid_cn_oracle.tsv")
META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
OUT_TSV = os.path.join(HERE, "rna_copy_number_depth.tsv")
OUT_JSON = os.path.join(HERE, "rna_copy_number_depth.json")


# --------------------------------------------------------------------------- small stats
def median(xs):
    s = sorted(xs)
    n = len(s)
    if n == 0:
        return 0.0
    m = n // 2
    return float(s[m]) if n % 2 else (s[m - 1] + s[m]) / 2.0


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def _ranks(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def spearman(xs, ys):
    return pearson(_ranks(list(xs)), _ranks(list(ys)))


# --------------------------------------------------------------------------- load
def load_fcn():
    rows = {}
    for r in csv.DictReader(open(FCN_TSV), delimiter="\t"):
        rows[int(r["family"])] = {k: r[k] for k in r}
    return rows


def load_tiers():
    """family -> {(chrom,start,end): tier}."""
    t = collections.defaultdict(dict)
    for r in csv.DictReader(open(SUN_TSV), delimiter="\t"):
        t[int(r["family"])][(r["chrom"], int(r["start"]), int(r["end"]))] = int(r["tier"])
    return t


def load_oracle():
    o = {}
    for r in csv.DictReader(open(ORACLE), delimiter="\t"):
        try:
            o[int(r["family"])] = float(r["asm_hapCN"])
        except (ValueError, KeyError):
            pass
    return o


def load_family_loci():
    fams = collections.defaultdict(list)
    with open(pg.FAM_TSV) as f:
        next(f)
        for line in f:
            fi, lid, c, s, e, nr = line.rstrip("\n").split("\t")
            fams[int(fi)].append((lid, c, int(s), int(e), int(nr)))
    return fams


def global_single_copy_anchor():
    """median n_reads over single-copy (non-family) multi-exon de-novo transcripts."""
    famids = set()
    with open(pg.FAM_TSV) as f:
        next(f)
        for line in f:
            famids.add(line.split("\t")[1])
    vals = []
    with open(META) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] in famids:
                continue
            if int(p[5]) >= 2:                 # multi-exon, matches family transcripts
                vals.append(int(p[6]))
    return median(vals), len(vals)


# --------------------------------------------------------------------------- per family
def family_depth(fam, fcn, tiers, loci, lam_global):
    ded = pg.dedup_copies(loci)                 # [(chrom,start,end,reads_summed)] == n_loci copies
    depths = [r[3] for r in ded]
    E_fam = sum(depths)
    tmap = tiers.get(fam, {})
    per_copy = []
    for (c, s, e, nr) in ded:
        per_copy.append({"chrom": c, "start": s, "end": e, "reads": nr,
                         "tier": tmap.get((c, s, e), 0)})
    lam_within = median(depths) if depths else 0.0

    n_loci = int(fcn["n_loci"])
    chi_H = int(fcn["chi_H"])
    collapsed_excess = int(fcn["collapsed_excess"])

    depth_cn_within = (E_fam / lam_within) if lam_within > 0 else float("nan")
    depth_cn_global = (E_fam / lam_global) if lam_global > 0 else float("nan")

    # confound magnitudes (C1)
    mean_d = E_fam / len(depths) if depths else 0.0
    var = sum((d - mean_d) ** 2 for d in depths) / len(depths) if depths else 0.0
    cv = (math.sqrt(var) / mean_d) if mean_d > 0 else 0.0
    dominant_share = (max(depths) / E_fam) if E_fam else 0.0

    # within-reference hidden multiplicity: reference copies whose depth is a MULTIPLE of the
    # family single-copy unit -> that one reference locus likely hides >1 genomic copy.
    hidden = sum(max(0, round(d / lam_within) - 1) for d in depths) if lam_within > 0 else 0

    return {
        "family": fam, "gene": fcn["gene"],
        "n_loci": n_loci, "chi_H": chi_H, "n_resolvable": int(fcn["n_resolvable"]),
        "collapsed_excess": collapsed_excess,
        "n_dedup_copies": len(depths), "E_fam_reads": E_fam,
        "lambda_within": round(lam_within, 3),
        "depth_cn_within": round(depth_cn_within, 2) if depth_cn_within == depth_cn_within else None,
        "depth_cn_global": round(depth_cn_global, 2) if depth_cn_global == depth_cn_global else None,
        "within_ref_hidden_copies": hidden,
        "depth_cn_within_plus_hidden": None,   # filled below
        "expr_cv": round(cv, 3), "dominant_copy_share": round(dominant_share, 3),
        "all_collapsed_no_tier1": int(sum(1 for pc in per_copy if pc["tier"] == 1) == 0),
        "per_copy_depths": sorted(depths, reverse=True),
        "_n_dedup_matches_n_loci": int(len(depths) == n_loci),
    }


def main():
    fcn = load_fcn()
    tiers = load_tiers()
    loci = load_family_loci()
    oracle = load_oracle()
    lam_global, n_singlecopy = global_single_copy_anchor()

    rows = []
    for fam in sorted(fcn):
        if fam not in loci:
            continue
        r = family_depth(fam, fcn[fam], tiers, loci[fam], lam_global)
        r["asm_hapCN_dna"] = oracle.get(fam)             # VALIDATION ONLY
        rows.append(r)

    # cross-check: dedup copy count must equal the shipped n_loci (same pipeline path)
    mismatch = [r["family"] for r in rows if not r["_n_dedup_matches_n_loci"]]

    collapsed = [r for r in rows if r["collapsed_excess"] > 0]     # the 21 "need depth/DNA"

    # ---------------- validation vs DNA (asm_hapCN), RNA is primary; DNA never gates
    val = [r for r in rows if r["asm_hapCN_dna"] and r["asm_hapCN_dna"] > 0]

    def corr(getter, label, subset):
        xs = [getter(r) for r in subset]
        ys = [r["asm_hapCN_dna"] for r in subset]
        keep = [(x, y) for x, y in zip(xs, ys) if x is not None and not (isinstance(x, float) and x != x)]
        if len(keep) < 3:
            return {"label": label, "n": len(keep), "pearson": None, "spearman": None}
        xs2 = [k[0] for k in keep]
        ys2 = [k[1] for k in keep]
        return {"label": label, "n": len(keep),
                "pearson": round(pearson(xs2, ys2), 3),
                "spearman": round(spearman(xs2, ys2), 3)}

    corrs = [
        corr(lambda r: r["n_loci"], "n_loci vs asm_hapCN", val),
        corr(lambda r: r["chi_H"], "chi_H vs asm_hapCN", val),
        corr(lambda r: r["depth_cn_within"], "depth_cn_within vs asm_hapCN", val),
        corr(lambda r: r["depth_cn_global"], "depth_cn_global vs asm_hapCN", val),
        corr(lambda r: r["depth_cn_within"], "depth_cn_within vs asm_hapCN (COLLAPSED fams)",
             [r for r in val if r["collapsed_excess"] > 0]),
    ]

    # ---------------- confound quantification (C1/C2)
    # C2 silent-copy undercount: fraction of families where depth_cn < asm_hapCN
    under = [r for r in val if r["depth_cn_within"] and r["depth_cn_within"] < r["asm_hapCN_dna"]]
    over = [r for r in val if r["depth_cn_within"] and r["depth_cn_within"] > r["asm_hapCN_dna"] + 0.5]
    # median ratio depth/DNA (how far below truth on average)
    ratios = sorted(r["depth_cn_within"] / r["asm_hapCN_dna"] for r in val if r["depth_cn_within"])
    med_ratio = median(ratios)
    # C1 dominant-copy families
    dom = [r for r in rows if r["dominant_copy_share"] >= 0.5 and r["n_loci"] >= 3]

    # ---------------- fill combined estimate + honest lower-bound field
    for r in rows:
        dw = r["depth_cn_within"]
        r["depth_cn_within_plus_hidden"] = (
            round(dw + r["within_ref_hidden_copies"], 2) if dw is not None else None)
        # honest RNA lower bound on genomic CN = max of the reference count and the depth count
        base = max(r["n_loci"], r["chi_H"])
        r["rna_cn_lower_bound"] = max(base, int(math.floor(dw)) if dw else base)

    # write TSV
    cols = ["family", "gene", "n_loci", "chi_H", "n_resolvable", "collapsed_excess",
            "E_fam_reads", "lambda_within", "depth_cn_within", "within_ref_hidden_copies",
            "depth_cn_within_plus_hidden", "depth_cn_global", "rna_cn_lower_bound",
            "expr_cv", "dominant_copy_share", "all_collapsed_no_tier1", "asm_hapCN_dna"]
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # named flagship collapsed families
    def _named(fam):
        r = next((x for x in rows if x["family"] == fam), None)
        if not r:
            return None
        return {k: r[k] for k in ("family", "gene", "n_loci", "chi_H", "collapsed_excess",
                                  "E_fam_reads", "lambda_within", "depth_cn_within",
                                  "depth_cn_within_plus_hidden", "depth_cn_global",
                                  "per_copy_depths", "asm_hapCN_dna")}

    summary = {
        "n_families": len(rows),
        "global_single_copy_anchor_lambda": lam_global,
        "n_singlecopy_transcripts_for_anchor": n_singlecopy,
        "dedup_count_matches_shipped_n_loci": len(mismatch) == 0,
        "n_dedup_mismatch": len(mismatch),
        "estimator_stack": {
            "n_loci": "reference-locus count (Tier1+2+3); misses reference-collapsed copies",
            "chi_H": "conflict chromatic; counts distinguishable copies; misses Tier-3 collapsed_excess",
            "depth_cn": "expressed-copy count incl identical; confounded by expression; blind to silent",
        },
        "totals": {
            "sum_n_loci": sum(r["n_loci"] for r in rows),
            "sum_chi_H": sum(r["chi_H"] for r in rows),
            "sum_collapsed_excess": sum(r["collapsed_excess"] for r in rows),
            "sum_within_ref_hidden_copies": sum(r["within_ref_hidden_copies"] for r in rows),
        },
        "collapsed_families_need_depth": {
            "n": len(collapsed),
            "note": "families with collapsed_excess>0 -- chi(H) undercounts; depth is the only RNA signal",
            "families": [{"family": r["family"], "gene": r["gene"], "n_loci": r["n_loci"],
                          "chi_H": r["chi_H"], "collapsed_excess": r["collapsed_excess"],
                          "depth_cn_within": r["depth_cn_within"],
                          "depth_cn_global": r["depth_cn_global"],
                          "asm_hapCN_dna": r["asm_hapCN_dna"]} for r in
                         sorted(collapsed, key=lambda x: -x["collapsed_excess"])],
        },
        "validation_vs_DNA": {
            "note": "asm_hapCN = DNA per-haplotype genomic CN (max(mat,pat)); VALIDATION ONLY, "
                    "never gates. RNA depth is a LOWER bound (silent copies invisible).",
            "n_families_with_dna": len(val),
            "correlations": corrs,
            "silent_copy_undercount": {
                "n_depth_below_dna": len(under), "n_depth_above_dna": len(over),
                "median_depth_over_dna_ratio": round(med_ratio, 3),
                "interpretation": "median depth-CN is this fraction of the DNA copy number; "
                                  "<1 => systematic under-count = confound C2 (silent copies)",
            },
        },
        "confounds": {
            "C1_expression_variation": {
                "n_families_dominant_copy_ge50pct_and_nloci_ge3": len(dom),
                "examples": [{"family": r["family"], "gene": r["gene"],
                              "dominant_copy_share": r["dominant_copy_share"],
                              "expr_cv": r["expr_cv"], "per_copy_depths": r["per_copy_depths"]}
                             for r in sorted(dom, key=lambda x: -x["dominant_copy_share"])[:6]],
            },
            "C2_silent_copies": "see validation_vs_DNA.silent_copy_undercount (median ratio < 1)",
            "C3_conclusion": "depth-CN is a NOISY LOWER BOUND on genomic CN, not an unbiased estimate",
        },
        "flagship_examples": {
            "GSTM2_fam0": _named(0),
            "LOC109025447_fam1_all_identical": _named(1),
            "LOC115932259_fam69_K0_copy_vs_allele": _named(69),
            "MAGEA9_fam94": _named(94),
        },
        "deterministic": True,
    }

    json.dump({"summary": summary, "families": rows}, open(OUT_JSON, "w"), indent=1)

    print(f"wrote {OUT_TSV} and {OUT_JSON}  ({len(rows)} families)")
    print(f"global single-copy anchor lambda = {lam_global} reads "
          f"(median over {n_singlecopy} single-copy multi-exon transcripts)")
    print(f"dedup count == shipped n_loci for all families: {len(mismatch) == 0}")
    print(f"collapsed families (collapsed_excess>0, need depth): {len(collapsed)}")
    print("VALIDATION vs DNA asm_hapCN (RNA primary, DNA never gates):")
    for c in corrs:
        print(f"  {c['label']:52s} n={c['n']:3d}  pearson={c['pearson']}  spearman={c['spearman']}")
    sv = summary["validation_vs_DNA"]["silent_copy_undercount"]
    print(f"silent-copy undercount: {sv['n_depth_below_dna']} families depth<DNA, "
          f"{sv['n_depth_above_dna']} depth>DNA, median depth/DNA ratio={sv['median_depth_over_dna_ratio']}")
    print("FLAGSHIP collapsed families:")
    for k, v in summary["flagship_examples"].items():
        if v:
            print(f"  {k}: n_loci={v['n_loci']} chi_H={v['chi_H']} collapsed_excess={v['collapsed_excess']} "
                  f"depth_cn_within={v['depth_cn_within']} depth_cn_global={v['depth_cn_global']} "
                  f"DNA_asm_hapCN={v['asm_hapCN_dna']} per_copy={v['per_copy_depths']}")


if __name__ == "__main__":
    main()
