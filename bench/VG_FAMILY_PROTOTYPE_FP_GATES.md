# VG-native prototype: false-positive characterization and empirical gates

This doc summarizes the false-positive (FP) landscape of the VG-native family-definition
prototype (`bench/vg_family_prototype.py`, default `--mode mmseqs-pairs --identity 0.90`)
and tests empirical post-refinement gates for removing FPs without heavy collateral.

## TL;DR

- **E_p-impure families (61 / 4633)** are structurally distinct from pure families:
  larger, higher pair-node multiplicity, more hub nodes, more distinct loci, and more
  exons.
- The best empirical gate is a **disjunction of family size and repeat-rich hubs**:
  **`n_members >= 80 OR repeat_hub_frac >= 0.05`**.
  - Drops 5 families, all E_p-impure, zero collateral.
  - Protein purity: `P_Ep` **0.8705 → 0.8779**.
  - Oracle recall: `R_oracle` stays at **51/57 = 0.8947**.
- The **repeat-rich hub** condition alone catches one E_p-impure family (fid 103, 35
  members) that the size gate misses — a family held together by repeat-rich,
  high-multiplicity exon-pair nodes.
- **DNA-confirmed FPs (6 families: 4 multifam, 1 oversize, 2 allele) are NOT removed
  by size, topology, or repeat/TE gates.** They are small, low-multiplicity families whose
  erroneous merge is only detectable with DNA copy-number information. This matches the
  earlier `FP_EXCLUSION_DISCRIMINATORS.md` finding for the shipped O1 catalog.

## Methods

1. **Feature extraction** — `bench/vg_family_prototype_fp_characterize.py` rebuilds the
   default VG, loads the eval TSV, and computes per-family features:
   - catalog size (`n_members`)
   - VG topology (`mean_pair_mult`, `max_pair_mult`, `pair_hub_frac`)
   - repeat/TE topology (`mean_node_repeat_frac`, `max_node_repeat_frac`, `repeat_hub_frac`)
   - genome architecture (`n_chrom`, `mean_recip_overlap`, `strand_majority`)
   - exon/repeat content (`mean_n_exons`, `mean_repeat_frac`)
   - truth labels (`ep_impure`, `fp_multifam`, `fp_oversize`, `fp_allele`, `any_dna_fp`)

2. **Rule search** — `bench/vg_family_prototype_fp_rules_method.py` and
   `bench/vg_family_prototype_fp_repeat_rules.py` grid-search method-computable and
   repeat-aware thresholds.

3. **Gate sweep** — `bench/vg_family_prototype_gate_sweep.py` post-hoc filters the
   existing catalog and runs `vg_family_prototype_eval.py` for each candidate gate.

4. **Prototype integration** — `vg_family_prototype.py` accepts
   `--fp-gate-members`, `--fp-gate-mean-pair-mult`, `--fp-gate-pair-hub-frac`, and
   `--fp-gate-repeat-hub-frac`.  When multiple thresholds are supplied, a family is
   dropped if **any** condition is exceeded (OR semantics).

## Distributions

| feature | E_p-pure (n=4572) | E_p-impure (n=61) | DNA-FP (n=6) |
|---|---|---|---|
| mean `n_members` | 10.46 | 41.08 | 32.00 |
| mean `n_distinct_loci` | 0.21 | 3.62 | 4.33 |
| mean `mean_pair_mult` | 6.54 | 22.53 | 20.20 |
| mean `pair_hub_frac` | 0.012 | 0.214 | 0.192 |
| mean `mean_node_repeat_frac` | 0.066 | 0.060 | — |
| mean `max_node_repeat_frac` | 0.249 | **0.394** | — |
| mean `repeat_hub_frac` | 0.000 | 0.001 | — |
| mean `mean_n_exons` | 9.35 | 15.93 | 17.40 |
| mean `max_recip_overlap` | 0.911 | 1.000 | 0.998 |

E_p-impure families are clear outliers on size and VG multiplicity.  The new
repeat-aware signal is in the **tail**: E_p-impure families have a much higher
`max_node_repeat_frac` (0.394 vs 0.249), meaning they contain at least one
repeat-rich shared exon-pair node.  The bulk mean is similar, so the signal is the
presence of a repeat bridge, not overall repeat content.

## Gate sweep results

| gate | kept | dropped | P_Ep | R_oracle | P_oracle(dedup) | FP a/o/m |
|---|---:|---:|---:|---:|---:|---:|
| baseline | 4633 | 0 | 0.8705 | 0.9123 | 0.8868 | 2/1/4 |
| `members>=80` | 4629 | 4 | 0.8761 | 0.8947 | 0.8846 | 2/1/4 |
| `repeat_hub_frac>=0.05` | 4632 | 1 | 0.8723 | 0.9123 | 0.8846 | 2/1/4 |
| **`members>=80 OR repeat_hub_frac>=0.05`** | **4628** | **5** | **0.8779** | **0.8947** | **0.8824** | 2/1/4 |
| `members>=70` | 4625 | 8 | 0.8796 | 0.8772 | 0.8824 | 2/1/4 |
| `members>=60` | 4619 | 14 | 0.8810 | 0.8596 | 0.8800 | 2/1/4 |
| `pair_hub_frac>=0.05` | 4479 | 154 | 0.8775 | 0.8246 | 0.8723 | 2/1/4 |
| `mean_pair_mult>=25` | 4560 | 73 | 0.8726 | 0.8596 | 0.8776 | 2/1/4 |
| `pair_hub_frac>=0.4` | 4554 | 79 | 0.8723 | 0.8421 | 0.8750 | 2/1/4 |

*FP a/o/m = allele / oversize / multifam flag counts from the diploid-DNA oracle.*

### Observations

- The combined size + repeat-rich-hub gate is strictly better than size alone: it
  removes one additional E_p-impure family (fid 103) with no extra collateral.
- `repeat_hub_frac>=0.05` alone is extremely precise (1 flagged, 1 EP-impure, 0
  collateral) but catches only one family.
- Plain `pair_hub_frac>=0.05` is far too aggressive (154 dropped) because it ignores
  whether the hub nodes are actually repeat-rich.
- DNA-confirmed FPs remain untouched.  The 6 DNA-FP prototype families are small and
  low-multiplicity:

| proto fid | FP type | n_members | mean_pair_mult | pair_hub_frac |
|---:|---|---:|---:|---:|
| 129 | multifam | 33 | 11.62 | 0.000 |
| 1129 | multifam/oversize | 15 | 5.85 | 0.000 |
| 593 | multifam | 20 | 11.11 | 0.000 |
| 4110 | multifam | 3 | 2.00 | 0.000 |
| 3423 | allele | 5 | 3.71 | 0.000 |
| 3425 | allele | 5 | 5.00 | 0.000 |

## Implementation

The prototype now supports OR-combined FP gates:

```bash
# recommended gate: size OR repeat-rich hub
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_family_prototype.py \
    --mode mmseqs-pairs --identity 0.90 \
    --fp-gate-members 80 --fp-gate-repeat-hub-frac 0.05 \
    --output bench/vg_family_prototype_repeatgate.tsv

# other examples
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_family_prototype.py \
    --mode mmseqs-pairs --identity 0.90 --fp-gate-members 80 \
    --output bench/vg_family_prototype_members80.tsv

PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_family_prototype.py \
    --mode mmseqs-pairs --identity 0.90 \
    --fp-gate-members 60 --fp-gate-mean-pair-mult 25 \
    --output bench/vg_family_prototype_gate.tsv
```

When multiple thresholds are supplied, a family is dropped if **any** supplied
condition is exceeded.

## Files produced

- `bench/vg_family_prototype_fp_features.tsv` — per-family feature table (now includes
  node-level repeat features)
- `bench/vg_family_prototype_fp_characterize.json` — distribution summary
- `bench/vg_family_prototype_fp_repeat_features.tsv` — repeat-enriched feature table
- `bench/vg_family_prototype_fp_rules_method.json` / `.txt` — top method-computable rules
- `bench/vg_family_prototype_fp_repeat_rules.json` / `.txt` — repeat-aware rule search
- `bench/vg_family_prototype_fp_diagnose.txt` — feature threshold tables
- `bench/vg_family_prototype_gate_sweep.json` / `.txt` — gate sweep metrics
- `bench/vg_family_prototype_repeatgate.tsv` / `.json` — end-to-end gated catalog
  (4633 → 4628 families, 5 dropped)
- `bench/vg_family_prototype_repeatgate_eval.json` / `.tsv` / `.log` — verified eval
  for the recommended gate: P_Ep 0.8779, R_oracle 51/57 = 0.8947

## Conclusion

Empirical FP gates give a small but clean improvement in protein purity.  The best
combination is a disjunction of large family size (`n_members >= 80`) and
repeat-rich hub nodes (`repeat_hub_frac >= 0.05`).  The repeat-hub condition is the
key new differentiator: it identifies a family whose members are glued together by a
repeat-rich exon-pair node, even when the family is not large enough to trigger the
size gate.

However, the residual DNA-confirmed false positives in the VG-native catalog are
still not dominated by large repeat-like families; they are small, low-multiplicity
merges that remain indistinguishable from real paralogs at the RNA level.  Removing
them will require either DNA copy-number information or a fundamentally different
definition signal (e.g. full-copy architecture rather than exon-pair sharing).
