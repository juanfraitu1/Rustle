# Gene-family-level (block-level) precision / recall of the CURRENT DEFAULT family definition

**Current default = the FLOORED catalog** `bench/family_rna_refine.tsv` — the RNA-only refined
oracle SCOPED by the ABSOLUTE DIVERGENCE FLOOR:
`core_recip >= 0.19 AND aln_frac >= 0.24 AND recip_id_best >= 0.80` (divergence floor, default-ON)
`-> gamma-quasi-clique(0.20) -> recombinant-split -> multi-repeat-bridge -> allele-demote`.
**307 multi-copy families** (md5 `e84dc2bc`). The catalog is **loaded, not re-derived**; all numbers
come from the shipped eval machinery (`graph_def_refine_sweep.eval_partition`,
`rna_only_edge_oracle.oracle_residuals`) so they match the committed catalogs.

The divergence floor SCOPES the definition to gene families whose members are **>= 80% reciprocal
whole-transcript identity** (Soto SD98-style; see `bench/DIVERGENCE_FLOOR.md`). It is **RNA-only**
(the reciprocal identity is from the RNA spliced-exon alignment, NOT DNA). It DROPS recall vs the
divergent-inclusive catalog **BY DESIGN**; that catalog is recovered exactly via the opt-out.

- **Opt-out** (`family_rna_refine.py --no-divergence-floor` / `--min-family-identity=0` /
  `RUSTLE_NO_DIVERGENCE_FLOOR=1`) = the divergent-inclusive catalog **573 families, md5 `dca64cbd`,
  byte-identical** (R_oracle 50/57 = 0.8772, E_p 0.8918, distinct-FP 4).
- **Legacy** (`--legacy`) = the shipped `core_recip >= 0.13` gamma-quasi-clique, **858 families**
  (reproduces the committed `graph_def_refine_sweep.json` baseline exactly).

Reproduce: `PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/family_level_pr_current.py`
Outputs (byte-identical across runs): `bench/family_level_pr_current.{json,tsv}`. The JSON is
self-documenting about the scoping regime (`divergence_floor.min_family_identity = 0.80`). All three
primary metrics are computed at the family/block level (they loop over blocks), NOT at the edge level.

---

## The three truths, clearly labelled (FLOORED default vs the floor-off opt-out)

### Truth 1 — PROTEIN E_p purity (block precision; conservative lower bound)
A family block is **PURE** iff its member genes span at most one non-mega protein family E_p;
**IMPURE = over-merge**. `Precision_Ep = pure_blocks / total_blocks`.

| | pure / total | impure | **P_Ep** |
|---|---|---|---|
| **DEFAULT (floor 0.80)** | 297 / 307 | 10 | **0.9674** |
| opt-out (floor off) | 511 / 573 | 62 | 0.8918 |
| legacy | 743 / 858 | 115 | 0.866 |

The floor lifts E_p purity **0.892 -> 0.967** (impure blocks 62 -> 10) by cutting divergent /
domain-share over-merges. **Caveat:** E_p is still a conservative lower bound (it splits genuinely-real
divergent-paralog families into separate PRFAMs, counted as "impurity").

### Truth 2 — DNA-loose cDNA (`in_dna_loose`)
Gene-pair truth from `family_er_pr.tsv`, cDNA id >= 0.90 & cov >= 0.30. The floor removes many
divergent/partial homolog pairs, so the block-level component-recovery recall of real cDNA families
DROPS (this is the scoping, expected): **0.923 (opt-out) -> ~0.62 (floored)**. Reported for continuity;
the load-bearing precision signal is the diploid oracle (Truth 3) + E_p (Truth 1).

### Truth 3 — INDEPENDENT DIPLOID DNA ORACLE (the gold, assembly-independent)
The only assembly-independent scalar (`asm_hapCN` = maternal + paternal diploid copy number).
`P_oracle = 1 - (multifam + oversize + allele) / oracle_mapped`;
`R_oracle = oracle_multicopy_genes_recovered / oracle_multicopy_genes`.

| metric | **DEFAULT (floor 0.80)** | opt-out (floor off) | legacy |
|---|---|---|---|
| oracle-mapped families | 31 | 48 | 50 |
| FP (allele / oversize / multifam) | **0 / 1 / 2** | 0 / 1 / 4 | 2 / 4 / 6 |
| flag-instances / **distinct FP blocks** | 3 / **2** | 5 / 4 | 12 / 11 |
| **P_oracle (task formula 1 - flags/mapped)** | **0.9032** | 0.896 | 0.76 |
| **P_oracle (dedup 1 - distinct/mapped)** | **0.9355** | 0.917 | 0.78 |
| **R_oracle (multi-copy genes)** | **32 / 57 = 0.5614** | 50 / 57 = 0.8772 | 50 / 57 = 0.877 |

**The floor is materially cleaner on the gold oracle: distinct FP 4 -> 2, P_dedup 0.917 -> 0.936, E_p
0.892 -> 0.967 — at a recall cost 0.877 -> 0.561 (18 divergent oracle genes lost BY DESIGN).** The
recall is recoverable via the opt-out (`--no-divergence-floor` = dca64cbd, R_oracle 0.877).

---

## Residual false-positive set the FLOORED default STILL gets (named)

### (a) DNA-CONFIRMED over-merges — 2 distinct blocks (the whole distinct-FP residual)
1. **`LOC129529978 + LOC129529986`** — the irreducible **MAGE-X** DNA-only array over-size floor
   (multifam + oversize; no RNA metric removes it — it is a genuine collapsed/inflated-copy DNA fact).
2. **`GSTM2 + LOC101129940`** — one **GST protein-domain hub** (multifam).

### (b) Residual allele-as-copy: **0** (demote removes DHRSX + LOC129530050; neither recurs).

### (c) E_p-impure blocks NOT reachable by the DNA oracle: **10** (the oracle covers 31/307 blocks).
Salient (from the roster JSON): FAM228B/MTHFD2L/ZNF626, SERINC4/TTPAL, GSPT2/SNX29, MORC2/OSBP2,
ALG1L1P/LOC101127338, FAM236B/LOC129530010, LOC115935268/LOC129534750, LOC115931432/LOC129525065,
LOC129526391/LOC129526471/LOC129531453 — small residual protein-family mixes (some may be real
divergent paralogs mis-split by E_p; do NOT all count against biological precision).

The named domain-shares the floor EXCLUDES (co-membered floor-off, gone floored): **RHD+SDHD**
(recip 0.785), **MOV10+RHOC** (0.11), **RBBP4+SYNC** (0.35), **DEDD+NIT1** (0.53).

---

## One-line takeaway
On the gold diploid oracle the FLOORED default (307 families, >= 80% reciprocal identity) reaches
**P_oracle(dedup) 0.936 / E_p purity 0.967 / distinct-FP 2** — the two residuals being the
irreducible MAGE-X DNA-only over-size floor and one GST domain hub — at **R_oracle 0.561 (32/57)**.
The recall is materially reduced vs the divergent-inclusive opt-out (0.877, dca64cbd, one flag away)
**by design** — this is the SD98-style scoping (`bench/DIVERGENCE_FLOOR.md`).

---

**Determinism:** `PYTHONHASHSEED=0`; `family_level_pr_current.{json,tsv}` byte-identical across
repeated runs; the floored catalog is md5 `e84dc2bc`, the opt-out is `dca64cbd`.

**Files (absolute):**
- `/mnt/c/Users/jfris/Desktop/Rustle/bench/family_level_pr_current.py`
- `/mnt/c/Users/jfris/Desktop/Rustle/bench/family_level_pr_current.tsv` (per-family flags)
- `/mnt/c/Users/jfris/Desktop/Rustle/bench/family_level_pr_current.json` (full P/R table +
  `divergence_floor` regime stamp + named FP roster)
- `/mnt/c/Users/jfris/Desktop/Rustle/bench/DIVERGENCE_FLOOR.md` (the P/R-vs-floor sweep + scoped claim)
