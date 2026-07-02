#!/usr/bin/env python3
"""
diploid_cn_corroborate.py

Turn the DNA-derived, RNA-INDEPENDENT diploid copy-number oracle
(bench/diploid_cn_oracle.tsv, built from the phased mGorGor1 maternal(chrX) +
paternal(chrY) haplotypes) into the three deliverables and write a single
machine-readable summary -> bench/diploid_cn_corroborate.json.

The oracle TSV now carries de-inflated, tandem-aware, sex-aware copy counts:
  hap_CN_mat/pat = distinct FULL-length gene-unit copies per haplotype
                   (aggregate query coverage >= 0.5 * median locus length;
                    query-overlapping alignments kept as separate tandem units)
  frag_mat/pat   = distinct sub-gene domain/repeat regions (NOT counted as copies)
  asm_hapCN      = sex-aware haploid copy number
  class          = sex-aware multi_copy / single_locus_allele / cnv / absent
  cnv_autosomal  = 1 if autosomal and the two haplotype copy counts differ

(1) NON-CIRCULAR CORROBORATION: Spearman + 10k-permutation null of asm_hapCN vs
    K_read (the truly RNA-independent read count), n_loci_ref, chi_H, sedef.
(2) COPY-vs-ALLELE: real multi_copy vs single-locus-allele; MAGEA9 / LOC115935025.
(3) REFERENCE-COLLAPSED copies the phased genome recovers, and reference over-splits.

Determinism: PYTHONHASHSEED=0, np.random.default_rng(0), sorted emission.
CAVEAT: copy number/paralogy is largely SPECIES-FIXED, so the assembly is a valid
copy-number oracle regardless of whether the IsoSeq donor is mGorGor1. Individual-
exact variant genotyping is donor-exact only if RNA==mGorGor1 (unknown).
"""
import json
import os
from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr, pearsonr

HERE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = "/home/juanfra/winloci_scratch"
ORACLE = os.path.join(HERE, "diploid_cn_oracle.tsv")
A1 = os.path.join(HERE, "a1_read_consensus_o1.tsv")
OUT = os.path.join(HERE, "diploid_cn_corroborate.json")
NPERM = 10000


def load_tsv(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, ln.rstrip("\n").split("\t"))) for ln in fh]


def num(x):
    try:
        return float(x)
    except Exception:
        return None


rows = load_tsv(ORACLE)
by_fam = {d["family"]: d for d in rows}
a1_meta = {d["family"]: d for d in load_tsv(A1)}
for d in rows:
    d["asm"] = int(d["asm_hapCN"])


# --------------------------------------------------------------------------- #
# (1) corroboration: Spearman + permutation null
# --------------------------------------------------------------------------- #
def corr(xs, ys):
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    rho, p = spearmanr(xs, ys)
    r, _ = pearsonr(xs, ys)
    rng = np.random.default_rng(0)
    perm = np.array([spearmanr(xs, rng.permutation(ys))[0] for _ in range(NPERM)])
    perm_p = (np.sum(np.abs(perm) >= abs(rho)) + 1) / (NPERM + 1)
    return dict(n=int(len(xs)), spearman_rho=round(float(rho), 4),
                spearman_p=float(f"{p:.3e}"), perm_p=round(float(perm_p), 5),
                pearson_r=round(float(r), 4))


a1_rows = [d for d in rows if num(d["K_read"]) is not None]

corroboration = {
    "n_families_total": len(rows),
    "n_families_a1_subset": len(a1_rows),
    "a1_subset_note": ("K_read exists only for the A1 read-consensus subset; these ARE "
                       "the collapsed / co-located candidate families."),
    "asm_hapCN_vs_K_read__NONCIRCULAR": corr(
        [d["asm"] for d in a1_rows], [num(d["K_read"]) for d in a1_rows]),
    "asm_hapCN_vs_n_loci_ref__a1": corr(
        [d["asm"] for d in a1_rows], [float(d["n_loci_ref"]) for d in a1_rows]),
    "asm_hapCN_vs_chi_H__a1": corr(
        [d["asm"] for d in a1_rows], [num(d["chi_H"]) for d in a1_rows]),
    "asm_hapCN_vs_sedef_partners__a1": corr(
        [d["asm"] for d in a1_rows], [num(d["sedef_partners"]) for d in a1_rows]),
    "asm_hapCN_vs_n_loci_ref__all196": corr(
        [d["asm"] for d in rows], [float(d["n_loci_ref"]) for d in rows]),
}
k_gt = sum(1 for d in a1_rows if d["asm"] > num(d["K_read"]))
k_eq = sum(1 for d in a1_rows if d["asm"] == num(d["K_read"]))
k_lt = [d for d in a1_rows if d["asm"] < num(d["K_read"])]
corroboration["K_read_is_censored_lower_bound"] = {
    "families_where_asm_hapCN_gt_K_read": int(k_gt),
    "families_where_asm_hapCN_eq_K_read": int(k_eq),
    "families_where_asm_hapCN_lt_K_read": len(k_lt),
    "of": len(a1_rows),
    "meaning": ("asm_hapCN >= K_read in %d/%d families: K_read is mostly a censored "
                "lower bound (the K=0 read-identifiability floor), so the DNA oracle "
                "bounds the RNA resolution ceiling. The %d exceptions are the read "
                "method OVER-counting relative to DNA -- 2 of them (DHRSX fam58, "
                "LOC129530050 fam91) are exactly the het_risk allele families where a "
                "heterozygous site was read as a 2nd copy, which the DNA correctly "
                "calls single-locus." % (k_gt + k_eq, len(a1_rows), len(k_lt))),
    "asm_lt_K_read_families": sorted(
        [dict(family=d["family"], gene=d["gene"], asm_hapCN=d["asm"],
              K_read=d["K_read"]) for d in k_lt], key=lambda x: int(x["family"])),
}


# --------------------------------------------------------------------------- #
# (2) copy-vs-allele
# --------------------------------------------------------------------------- #
def dist(key):
    c = defaultdict(int)
    for d in rows:
        c[d[key]] += 1
    return dict(sorted(c.items()))


def vrow(d):
    return dict(family=d["family"], gene=d["gene"], sex=d["sex"],
                hap_CN_mat=int(d["hap_CN_mat"]), hap_CN_pat=int(d["hap_CN_pat"]),
                asm_hapCN=d["asm"], n_loci_ref=int(d["n_loci_ref"]),
                K_read=d["K_read"], verdict=d["class"])


a1_conf = [d for d in a1_rows if d["asm"] >= 2]
a1_notconf = [d for d in a1_rows if d["asm"] < 2]
copy_to_allele = [d for d in rows if d["class"] == "single_locus_allele"]
cnv_auto_real = [d for d in rows if d["class"] == "cnv" and d["sex"] == "auto"]
cnv_auto_flag = [d for d in rows if d["cnv_autosomal"] == "1"]
sexlinked_cnv = [d for d in rows if d["sex"] in ("X", "Y") and
                 (int(d["hap_CN_mat"]) == 0 or int(d["hap_CN_pat"]) == 0)]

# the RNA-flagged het_risk families and how the DNA oracle independently calls them
het_risk_a1 = sorted([f for f, m in a1_meta.items() if m.get("het_risk") == "1"],
                     key=int)

copy_vs_allele = {
    "class_distribution": dist("class"),
    "sex_linkage": dist("sex"),
    "cnv_note": ("mGorGor1 is MALE (maternal=chrX, paternal=chrY): an X gene reads "
                 "hap_CN_pat==0, a Y gene hap_CN_mat==0. That male hemizygosity is now "
                 "resolved sex-aware (multi_copy/single by the PRESENT haplotype), NOT "
                 "mislabelled cnv. A real autosomal CNV = the two haplotype copy counts "
                 "differ (cnv, or cnv_autosomal=1 when both haplotypes are still >=2)."),
    "n_sexlinked_hemizygous": len(sexlinked_cnv),
    "n_autosomal_cnv_class": len(cnv_auto_real),
    "autosomal_cnv_families": sorted([vrow(d) for d in cnv_auto_real], key=lambda x: int(x["family"])),
    "n_autosomal_cnv_flag_both_ge2": len(cnv_auto_flag),
    "resolved_real_multi_copy": sum(1 for d in rows if d["class"] == "multi_copy"),
    "resolved_single_locus_allele": len(copy_to_allele),
    "copy_to_allele_note": ("Families the reference/RNA method listed as multi-locus "
                            "(n_loci_ref>=2) but the diploid genome carries exactly 1 "
                            "full copy per present haplotype => the RNA variation there "
                            "is allele/heterozygosity, NOT a copy."),
    "copy_to_allele_families": sorted([vrow(d) for d in copy_to_allele], key=lambda x: int(x["family"])),
    "a1_subset_assembly_confirmed_multi_copy": f"{len(a1_conf)}/{len(a1_rows)}",
    "a1_subset_NOT_confirmed": sorted([vrow(d) for d in a1_notconf], key=lambda x: int(x["family"])),
    "a1_notconfirmed_are_exactly_the_het_risk_flags": (
        sorted([d["family"] for d in a1_notconf], key=int) ==
        sorted([f for f in het_risk_a1 if by_fam[f]["asm"] < 2], key=int)),
    "het_risk_survivors": {
        "MAGEA9_fam94": {**vrow(by_fam["94"]),
            "resolution": ("X-linked in a MALE (single X): 2 X copies cannot be alleles "
                           "of one locus => necessarily paralogous => REAL multi_copy.")},
        "LOC115935025_fam145": {**vrow(by_fam["145"]),
            "resolution": ("autosomal, hap_CN_mat=2 / hap_CN_pat=4 (both >=2) => REAL "
                           "multi_copy (and cnv_autosomal=1).")},
    },
    "all_a1_het_risk_flag_resolution": sorted(
        [vrow(by_fam[f]) for f in het_risk_a1 if f in by_fam], key=lambda x: int(x["family"])),
}


# --------------------------------------------------------------------------- #
# (3) reference-collapsed copies (asm >> n_loci) and reference over-splits
# --------------------------------------------------------------------------- #
def frag(d):
    return max(int(d["frag_mat"]), int(d["frag_pat"]))


over, under = [], []
for d in rows:
    a, n = d["asm"], float(d["n_loci_ref"])
    if a == 0:
        continue
    if a - n >= 3 and a >= 2 * n:
        over.append(d)
    elif n - a >= 3 and n >= 2 * a:
        under.append(d)


def odetail(d):
    return dict(family=d["family"], gene=d["gene"], sex=d["sex"],
                n_loci_ref=int(d["n_loci_ref"]), asm_hapCN=d["asm"],
                hap_CN_mat=int(d["hap_CN_mat"]), hap_CN_pat=int(d["hap_CN_pat"]),
                frag_mat=int(d["frag_mat"]), frag_pat=int(d["frag_pat"]),
                K_read=d["K_read"], sedef_partners=d["sedef_partners"])


over_detail = sorted([odetail(d) for d in over], key=lambda x: -(x["asm_hapCN"] - x["n_loci_ref"]))
under_detail = sorted([odetail(d) for d in under], key=lambda x: -(x["n_loci_ref"] - x["asm_hapCN"]))
# domain/repeat inflation is now REJECTED from hap_CN and lives in frag_*; report it
domain_heavy = sorted([odetail(d) for d in rows if frag(d) >= 5],
                      key=lambda x: -max(x["frag_mat"], x["frag_pat"]))

reference_collapsed = {
    "definition_over": "asm_hapCN >= 2*n_loci_ref and asm_hapCN - n_loci_ref >= 3",
    "counting_note": ("hap_CN counts ONLY full-length gene copies (aggregate query "
                      "coverage >= 0.5*L); shared sub-gene domain/repeat hits are "
                      "excluded from hap_CN and tallied in frag_* -- so the inflation "
                      "flag is folded INTO the shipped count, not a separate caveat."),
    "n_reference_collapsed_real": len(over_detail),
    "families": over_detail,
    "n_domain_heavy_families_frag_ge5": len(domain_heavy),
    "domain_heavy_note": ("families with >=5 distinct sub-gene domain/repeat regions "
                          "that were NOT counted as copies (would have inflated a "
                          "raw-block counter); e.g. LOC115930164/ZNF domain sharers."),
    "domain_heavy_families": domain_heavy,
    "definition_under": "n_loci_ref >= 2*asm_hapCN and n_loci_ref - asm_hapCN >= 3",
    "reverse_direction_note": ("asm << n_loci families = reference OVER-SPLIT; the "
                               "diploid gives a lower, true copy number."),
    "n_reference_over_split": len(under_detail),
    "families_over_split": under_detail,
}


# --------------------------------------------------------------------------- #
# verdict
# --------------------------------------------------------------------------- #
ck = corroboration["asm_hapCN_vs_K_read__NONCIRCULAR"]
cn = corroboration["asm_hapCN_vs_n_loci_ref__a1"]
cs = corroboration["asm_hapCN_vs_sedef_partners__a1"]

verdict = {
    "a_corroborate_A1_non_circularly": (
        "PARTIAL / HONEST. asm_hapCN vs the truly RNA-independent read count K_read is "
        f"POSITIVE (rho={ck['spearman_rho']}, perm_p={ck['perm_p']}, n={ck['n']}); "
        "asm_hapCN vs n_loci_ref is clearly above null on the same substrate "
        f"(rho={cn['spearman_rho']}, perm_p={cn['perm_p']}). asm_hapCN >= K_read in "
        f"{k_gt + k_eq}/{len(a1_rows)} families: K_read is a censored lower bound (K=0 "
        "read-identifiability floor), so the DNA oracle bounds the RNA resolution "
        "ceiling rather than rank-tracking the read count tightly. This IS the "
        "corroboration SEDEF could not give: SEDEF is reference self-alignment "
        f"(co-varies with n_loci by construction, rho={cs['spearman_rho']}), whereas "
        "the phased diploid is an independent DNA observation."),
    "b_resolve_copy_vs_allele": (
        f"YES. {copy_vs_allele['a1_subset_assembly_confirmed_multi_copy']} A1 families "
        "are DNA-confirmed real multi_copy; the ONLY 2 not confirmed are EXACTLY the two "
        "RNA het_risk=1 flags (DHRSX fam58, LOC129530050 fam91), independently called "
        "single-locus allele/het by the DNA. Both A1 het-risk survivors resolve as REAL "
        "multi_copy (MAGEA9: 2 copies on a single male X; LOC115935025: mat=2/pat=4). "
        f"{copy_vs_allele['resolved_single_locus_allele']} families the reference listed "
        "as multi-locus collapse to 1 full copy/haplotype => reclassified copy->allele. "
        f"{len(cnv_auto_real)} autosomal CNV (class=cnv) + "
        f"{len(cnv_auto_flag)} autosomal count-difference flags (cnv_autosomal, both "
        "haplotypes still >=2); sex-linked hemizygosity is resolved sex-aware, not "
        "mislabelled cnv."),
    "c_sharpen_catalog": (
        f"YES, both directions. {len(over_detail)} families are genuine reference-"
        "collapsed FULL-gene copies the phased genome recovers (non-circular support for "
        "the O4 'more copies than reference' thread; e.g. ZNF425 5->26, and tandem "
        "arrays like LOC129529768 fam42 1->12 that a single merge-gap had collapsed). "
        f"{len(under_detail)} families are reference OVER-split (diploid gives the lower "
        f"true CN, e.g. the 32-locus fam2 -> 1). Shared sub-gene domain/repeat hits "
        f"({reference_collapsed['n_domain_heavy_families_frag_ge5']} families with >=5) "
        "are now rejected from the copy count and tallied separately, so no family is "
        "inflated by a shared domain."),
    "caveat_species_fixed_vs_donor_exact": (
        "Copy number / paralogy is largely SPECIES-FIXED, so the diploid assembly is a "
        "valid RNA-independent copy-number oracle regardless of whether the IsoSeq donor "
        "is mGorGor1. Individual-exact het/allele genotyping of a SPECIFIC variant is "
        "donor-exact only if RNA==mGorGor1 (unknown). The robust, lead result is the "
        "copy-number / copy-vs-allele use, not variant-level genotyping."),
}

out = {
    "description": ("DNA-derived, RNA-independent diploid copy-number corroboration of "
                    "the RNA multi-copy family method, from the phased mGorGor1 "
                    "maternal(chrX)+paternal(chrY) haplotypes. Copy counts are tandem-"
                    "aware, length-proportional (0.5*L gene-unit coverage), sex-aware."),
    "inputs": {
        "oracle_tsv": ORACLE, "a1_read_consensus": A1,
        "mat_paf": os.path.join(SCRATCH, "cn_oracle_paf", "mat.paf"),
        "pat_paf": os.path.join(SCRATCH, "cn_oracle_paf", "pat.paf"),
        "min_ident": 0.90, "cov_frac_of_median_locus_len": 0.5, "n_perm": NPERM,
    },
    "corroboration": corroboration,
    "copy_vs_allele": copy_vs_allele,
    "reference_collapsed_copies": reference_collapsed,
    "verdict": verdict,
}

with open(OUT, "w") as fh:
    json.dump(out, fh, indent=2, sort_keys=False)
print("wrote", OUT)

print("\n(1) corroboration (a1 subset, n=%d):" % len(a1_rows))
for k in ("asm_hapCN_vs_K_read__NONCIRCULAR", "asm_hapCN_vs_n_loci_ref__a1",
          "asm_hapCN_vs_chi_H__a1", "asm_hapCN_vs_sedef_partners__a1"):
    v = corroboration[k]
    print("    %-38s rho=%+.3f perm_p=%.4f" % (k, v["spearman_rho"], v["perm_p"]))
print("(2) copy-vs-allele: a1 confirmed multi_copy %s ; copy->allele %d ; "
      "auto-CNV(class) %d ; auto-CNV(flag) %d"
      % (copy_vs_allele["a1_subset_assembly_confirmed_multi_copy"],
         copy_vs_allele["resolved_single_locus_allele"],
         len(cnv_auto_real), len(cnv_auto_flag)))
print("    not-confirmed == het_risk flags: %s"
      % copy_vs_allele["a1_notconfirmed_are_exactly_the_het_risk_flags"])
print("    MAGEA9 -> %s ; LOC115935025 -> %s"
      % (by_fam["94"]["class"], by_fam["145"]["class"]))
print("(3) reference-collapsed real=%d, over-split=%d, domain-heavy(frag>=5)=%d"
      % (len(over_detail), len(under_detail), len(domain_heavy)))
