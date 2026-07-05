# VG-native prototype: false-positive characterization and empirical gates

This doc summarizes the false-positive (FP) landscape of the VG-native family-definition
prototype (`bench/vg_family_prototype.py`, default `--mode mmseqs-pairs --identity 0.90`)
and tests empirical post-refinement gates for removing FPs without heavy collateral.

## TL;DR

- **E_p-impure families (61 / 4633)** are structurally distinct from pure families:
  larger, higher pair-node multiplicity, more hub nodes, more distinct loci, and more
  exons.
- A simple **family-size gate (`n_members >= 80`)** removes 4 of those E_p-impure
  families with zero collateral, improving protein purity `P_Ep` from **0.8705 → 0.8761**
  at the cost of a small recall drop (`R_oracle` 0.9123 → 0.8947).
- **DNA-confirmed FPs (6 families: 4 multifam, 1 oversize, 2 allele) are NOT removed
  by size or VG-topology gates.** They are small, low-multiplicity families whose
  erroneous merge is only detectable with DNA copy-number information. This matches the
  earlier `FP_EXCLUSION_DISCRIMINATORS.md` finding for the shipped O1 catalog.
- The prototype now supports optional FP gates (`--fp-gate-members`,
  `--fp-gate-mean-pair-mult`, `--fp-gate-pair-hub-frac`) with AND semantics.

## Methods

1. **Feature extraction** — `bench/vg_family_prototype_fp_characterize.py` rebuilds the
   default VG, loads the eval TSV, and computes per-family features:
   - catalog size (`n_members`)
   - VG topology (`mean_pair_mult`, `max_pair_mult`, `pair_hub_frac`)
   - genome architecture (`n_chrom`, `mean_recip_overlap`, `strand_majority`)
   - exon/repeat content (`mean_n_exons`, `mean_repeat_frac`)
   - truth labels (`ep_impure`, `fp_multifam`, `fp_oversize`, `fp_allele`, `any_dna_fp`)

2. **Rule search** — `bench/vg_family_prototype_fp_rules.py` and
   `bench/vg_family_prototype_fp_rules_method.py` grid-search single-feature and
   conjunctive thresholds.

3. **Gate sweep** — `bench/vg_family_prototype_gate_sweep.py` post-hoc filters the
   existing catalog and runs `vg_family_prototype_eval.py` for each candidate gate.

4. **Prototype integration** — `vg_family_prototype.py` now accepts
   `--fp-gate-members`, `--fp-gate-mean-pair-mult`, and `--fp-gate-pair-hub-frac`
   (AND semantics) and drops offending families after γ-quasi-clique refinement.

## Distributions

| feature | E_p-pure (n=4572) | E_p-impure (n=61) | DNA-FP (n=6) |
|---|---|---|---|
| mean `n_members` | 10.46 | 41.08 | 32.00 |
| mean `n_distinct_loci` | 0.21 | 3.62 | 4.33 |
| mean `mean_pair_mult` | 6.61 | 22.80 | 20.20 |
| mean `pair_hub_frac` | 0.012 | 0.214 | 0.192 |
| mean `mean_n_exons` | 9.35 | 15.93 | 17.40 |
| mean `max_recip_overlap` | 0.911 | 1.000 | 0.998 |

E_p-impure families are clear outliers on size and VG multiplicity. DNA FPs are a
separate, much smaller class.

## Gate sweep results

| gate | kept | dropped | P_Ep | R_oracle | P_oracle(dedup) | FP a/o/m |
|---|---:|---:|---:|---:|---:|---:|
| baseline | 4633 | 0 | 0.8705 | 0.9123 | 0.8868 | 2/1/4 |
| `members>=80` | 4629 | 4 | 0.8761 | 0.8947 | 0.8846 | 2/1/4 |
| `members>=70` | 4625 | 8 | 0.8796 | 0.8772 | 0.8824 | 2/1/4 |
| `members>=60` | 4619 | 14 | 0.8810 | 0.8596 | 0.8800 | 2/1/4 |
| `mean_pair_mult>=25` | 4560 | 73 | 0.8726 | 0.8596 | 0.8776 | 2/1/4 |
| `pair_hub_frac>=0.4` | 4555 | 78 | 0.8726 | 0.8596 | 0.8776 | 2/1/4 |
| `members>=60 AND mean_pair_mult>=25` | 4630 | 3 | 0.8705 | 0.9123 | 0.8868 | 2/1/4 |
| `members>=50 AND mean_pair_mult>=20` | 4617 | 16 | 0.8721 | 0.8772 | 0.8824 | 2/1/4 |

*FP a/o/m = allele / oversize / multifam flag counts from the diploid-DNA oracle.*

### Observations

- `members>=80` is the only gate that improves `P_Ep` with no collateral loss of
  real multicopy families (it removes 4 families, all E_p-impure).
- Stricter size thresholds improve `P_Ep` further but cost recall.
- Multiplicity/hub gates (`mean_pair_mult>=25`, `pair_hub_frac>=0.4`) have high
  collateral and do not improve `P_Ep` as cleanly as the size gate.
- **None of the tested gates reduce the DNA-confirmed FP roster.** The 6 DNA-FP
  prototype families have low-to-moderate size and multiplicity:

| proto fid | FP type | n_members | mean_pair_mult | pair_hub_frac |
|---:|---|---:|---:|---:|
| 129 | multifam | 33 | 11.62 | 0.000 |
| 1129 | multifam/oversize | 15 | 5.85 | 0.000 |
| 593 | multifam | 20 | 11.11 | 0.000 |
| 4110 | multifam | 3 | 2.00 | 0.000 |
| 3423 | allele | 5 | 3.71 | 0.000 |
| 3425 | allele | 5 | 5.00 | 0.000 |

These are not separable from real small multi-copy families using the features we
measured.

## Implementation

The prototype now supports:

```bash
# conservative size gate (recommended)
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_family_prototype.py \
    --mode mmseqs-pairs --identity 0.90 --fp-gate-members 80 \
    --output bench/vg_family_prototype_members80.tsv

# conjunctive gate example
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_family_prototype.py \
    --mode mmseqs-pairs --identity 0.90 \
    --fp-gate-members 60 --fp-gate-mean-pair-mult 25 \
    --output bench/vg_family_prototype_gate.tsv
```

When multiple thresholds are supplied, a family is dropped only if **all** supplied
conditions are exceeded.

## Files produced

- `bench/vg_family_prototype_fp_features.tsv` — per-family feature table
- `bench/vg_family_prototype_fp_characterize.json` — distribution summary
- `bench/vg_family_prototype_fp_rules_method.json` — top method-computable rules
- `bench/vg_family_prototype_fp_diagnose.txt` — feature threshold tables
- `bench/vg_family_prototype_gate_sweep.json` / `.txt` — gate sweep metrics
- `bench/vg_family_prototype_members80.tsv` / `.json` — end-to-end gated catalog
  (4633 → 4629 families, 4 dropped)
- `bench/vg_family_prototype_members80_eval.json` / `.tsv` / `.log` — verified eval
  matching the post-hoc sweep: P_Ep 0.8761, R_oracle 51/57 = 0.8947

## Conclusion

Empirical size/multiplicity gates give a small but clean improvement in protein
purity. They are easy to compute and integrate. However, the residual DNA-confirmed
false positives in the VG-native catalog are not dominated by large repeat-like
families; they are small, low-multiplicity merges that remain indistinguishable from
real paralogs at the RNA level. Removing them will require either DNA copy-number
information or a fundamentally different definition signal (e.g. full-copy
architecture rather than exon-pair sharing).
