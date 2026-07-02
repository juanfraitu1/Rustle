#!/usr/bin/env python
"""Self-check for the GENE-PROJECTION RELABEL fix in bench/family_level_pr_current.py.

Context (bench/CANDIDATE_GENERATION_RECALL.md, diagnosis d292e65): the diploid-oracle recall matches
oracle genes by NAME (the max-overlap gene_of projection). In the ARL17A / LRRC37A segdup cluster a
de-novo copy overlaps several paralogous genes, so max-overlap labels a recovered copy after a
NEIGHBOUR -- and two oracle multi-copy genes whose OWN copy IS grouped into a recovered multi-copy
family (LOC115934629, LOC129534585) are not credited by name -> counted as missed. The catalog
GROUPING is correct; only the gene-name projection/measurement misses them.

The fix (relabel_recovered_genes) credits an oracle gene by GENOMIC / copy linkage instead of name:
its own de-novo copy sits in a POA-passing cross-copy edge with an endpoint in a multi-copy catalog
family. It augments the recall NUMERATOR ONLY.

Asserts:
  (a) LOC115934629 and LOC129534585 are now credited as recovered multi_copy oracle genes; the
      diploid-oracle recall goes 48/57 -> 50/57 (recall_without_relabel 0.8421 -> 0.8772).
  (b) The family GROUPING is byte-identical: family_rna_refine.tsv md5 == the committed golden
      (the fix never touches the catalog / gene_of / block projection).
  (c) Precision / FP are UNCHANGED by the relabel: diploid-oracle distinct FP blocks == 6,
      P_oracle(dedup) == 0.875, P_oracle(task) == 0.8542, oracle-mapped blocks == 48, and the E_p
      block purity == 0.868 (526/606, 80 impure) -- the relabel raises recall only.
  (d) Determinism: two runs (PYTHONHASHSEED=0) produce a byte-identical JSON.
  (e) denovo_families.py and genome_family_def.py are untouched by the run.

Run: /home/juanfra/miniforge3/bin/python bench/test_family_level_pr_relabel.py
"""
import hashlib
import json
import os
import subprocess
import sys

BENCH = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(BENCH, "family_level_pr_current.py")
OUT_JSON = os.path.join(BENCH, "family_level_pr_current.json")
CATALOG_TSV = os.path.join(BENCH, "family_rna_refine.tsv")
DENOVO = os.path.join(BENCH, "denovo_families.py")
GENOME = os.path.join(BENCH, "genome_family_def.py")
PY = sys.executable

# md5 of the committed DEFAULT grouping catalog (== GOLDEN_DEFAULT_TSV_MD5 in test_family_rna_refine.py).
GROUPING_GOLDEN_MD5 = "f94f387e4f3a53a69e9d8a2d0b1f497a"
TARGETS = ["LOC115934629", "LOC129534585"]


def _md5(path):
    with open(path, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


def _run():
    env = dict(os.environ)
    env["PYTHONHASHSEED"] = "0"
    r = subprocess.run([PY, SCRIPT], cwd=BENCH, env=env, capture_output=True, text=True)
    assert r.returncode == 0, f"family_level_pr_current.py exit={r.returncode}\n{r.stderr[-2000:]}"
    return json.load(open(OUT_JSON)), _md5(OUT_JSON)


def main():
    # capture the two 'untouched' scripts before running
    denovo_before, genome_before = _md5(DENOVO), _md5(GENOME)

    d1, json_md5_1 = _run()
    cur = d1["truth3_diploid_oracle"]["current"]
    rel = d1["gene_projection_relabel"]

    # (a) the two targets credited; recall 48/57 -> 50/57
    assert rel["credited_genes"] == TARGETS, f"relabel credited_genes != {TARGETS}: {rel['credited_genes']}"
    assert rel["newly_credited"] == TARGETS, f"newly_credited != {TARGETS}: {rel['newly_credited']}"
    assert cur["relabel_credited_genes"] == TARGETS, cur["relabel_credited_genes"]
    assert cur["oracle_genes_recovered_multicopy"] == 50, \
        f"recall numerator != 50: {cur['oracle_genes_recovered_multicopy']}"
    assert cur["oracle_multicopy_genes"] == 57, cur["oracle_multicopy_genes"]
    assert cur["recall_oracle"] == 0.8772, f"R_oracle != 0.8772 (50/57): {cur['recall_oracle']}"
    assert rel["recall_without_relabel"] == 0.8421, \
        f"pre-relabel recall != 0.8421 (48/57): {rel['recall_without_relabel']}"
    assert rel["recall_with_relabel"] == 0.8772, rel["recall_with_relabel"]
    print(f"(a) LOC115934629 + LOC129534585 credited; R_oracle 48/57 (0.8421) -> 50/57 (0.8772) : OK")

    # (b) family GROUPING byte-identical (catalog untouched by the relabel)
    got = _md5(CATALOG_TSV)
    assert got == GROUPING_GOLDEN_MD5, \
        f"family_rna_refine.tsv grouping drifted from golden {GROUPING_GOLDEN_MD5}: {got}"
    print(f"(b) family_rna_refine.tsv grouping byte-identical (md5 {GROUPING_GOLDEN_MD5}) : OK")

    # (c) precision / FP unchanged by the relabel (recall-numerator-only fix)
    assert cur["distinct_fp_blocks"] == 6, f"diploid-oracle distinct FP blocks != 6: {cur['distinct_fp_blocks']}"
    assert cur["precision_oracle_dedup"] == 0.875, f"P_oracle(dedup) != 0.875: {cur['precision_oracle_dedup']}"
    assert cur["precision_oracle_taskformula"] == 0.8542, \
        f"P_oracle(task) != 0.8542: {cur['precision_oracle_taskformula']}"
    assert cur["oracle_mapped_families"] == 48, f"oracle-mapped blocks != 48: {cur['oracle_mapped_families']}"
    assert (cur["n_allele"], cur["n_oversize"], cur["n_multifam"]) == (0, 3, 4), \
        f"FP flag composition changed: {(cur['n_allele'], cur['n_oversize'], cur['n_multifam'])}"
    ep = d1["truth1_protein_Ep"]["current"]
    assert ep["precision_Ep"] == 0.868 and ep["total_blocks"] == 606 and ep["impure_blocks"] == 80, \
        f"E_p purity changed: {ep}"
    print("(c) precision/FP unchanged: diploid distinct-FP 6, P_oracle(dedup) 0.875, oracle-mapped 48, "
          "E_p purity 526/606=0.868 : OK")

    # (d) determinism: a second run yields a byte-identical JSON
    _, json_md5_2 = _run()
    assert json_md5_1 == json_md5_2, f"JSON not byte-identical across runs: {json_md5_1} vs {json_md5_2}"
    print(f"(d) deterministic (JSON byte-identical across runs; md5 {json_md5_1}) : OK")

    # (e) denovo_families.py / genome_family_def.py untouched
    assert _md5(DENOVO) == denovo_before, "denovo_families.py was modified by the run"
    assert _md5(GENOME) == genome_before, "genome_family_def.py was modified by the run"
    print("(e) denovo_families.py + genome_family_def.py untouched by the run : OK")


if __name__ == "__main__":
    main()
    print("\nALL TESTS PASSED")
