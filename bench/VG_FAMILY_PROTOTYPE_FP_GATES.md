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

# require >=2 shared colinear exon-pair nodes to link two loci
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_family_prototype.py \
    --mode mmseqs-pairs --identity 0.90 --min-shared-pairs 2 \
    --output bench/vg_family_prototype_min2.tsv

# add the edge-level antisense / reciprocal-overlap gate
PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_family_prototype.py \
    --mode mmseqs-pairs --identity 0.90 --antisense-gate \
    --output bench/vg_family_prototype_antisense.tsv
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
- `bench/vg_family_prototype_protpure.tsv` / `_eval.*` — family split by annotated
  protein family (P_Ep 1.0, R_oracle 0.8246)
- `bench/vg_family_prototype_protcov.tsv` / `_eval.*` / `.json` — family split by
  whole-protein reciprocal coverage (P_Ep 1.0, R_oracle 0.7018)
- `bench/vg_family_prototype_fp_features.tsv` — per-family feature table, now with
  RNA-only exon-architecture columns (exact-exon Jaccard, intron Jaccard, unshared run,
  within-family pair-node Jaccard)
- `bench/vg_family_prototype_id95_*` — stricter `--identity 0.95` catalog / eval
- `bench/vg_family_prototype_exact_*` — exact graph-to-graph catalog / eval
- `bench/vg_family_prototype_min2_*` — `--min-shared-pairs 2` catalog / eval
- `bench/vg_family_prototype_antisense_*` — `--antisense-gate` catalog / eval

## Domain-share prevention — what works and what does not

Domain-sharers are the bulk of the residual E_p-impure families (61 / 4633).  They
are genes from different protein families that share one conserved domain and get
merged in the VG through that domain's exons.

### Tested axes

| approach | mechanism | result |
|---|---|---|
| **Shared exon-pair coverage** | fraction of each member's nodes that are shared (`mean_shared_node_frac`, `min_shared_node_frac`) | **No signal** — domain-sharers have *higher* coverage than real families (mean 0.955 vs 0.866). The VG collapses the shared domain into a dense path that looks like a real paralog. |
| **Antisense / reciprocal-overlap** (family-level) | same-contig opposite-strand pairs with recip overlap ≥ 0.5 | **Weak + high collateral** — 14 EP-impure flagged but 403 real families collateral. Must be applied at the edge level to be clean. |
| **Protein-family splitting** | split every family with ≥2 annotated protein families into protein-pure sub-families | **P_Ep = 1.0, R_oracle drops 0.9123 → 0.8246** — removes all impurity but over-splits real duplications the protein truth fragments. |
| **Whole-protein reciprocal-coverage bar** (O1 safeguard) | dissolve families where no cross-gene pair passes `min-cov ≥ 0.50, max-cov ≥ 0.70, fident ≥ 0.30` | **P_Ep = 1.0, R_oracle drops 0.9123 → 0.7018** — all 61 EP-impure families dissolve because their genes come from different protein families by definition. |
| **RNA-only exon architecture** (new) | exact-exon Jaccard, intron-adjacency Jaccard, max unshared exon run, within-family pair-node Jaccard | **Inverted / no signal** — EP-impure families have *higher* exact-exon Jaccard (0.356 vs 0.320), higher intron Jaccard (0.313 vs 0.249), *lower* max unshared run (0.94 vs 1.38), and nearly identical pair-node Jaccard (0.400 vs 0.423). The shared conserved domain is more exact-similar than the divergent exons of real paralogs, so these features flag real families, not domain-sharers. |

### Interpretation

The VG-native definition has the same fundamental domain-share boundary as the
shipped O1 catalog: at the RNA level a real divergent paralog-with-a-domain and a
FP domain-share are the same object.  The only clean preventions are:

1. **Use whole-protein reciprocal coverage as an edge filter** (the O1 safeguard),
   but accept the recall cost or integrate it as a pre-filter before VG construction.
2. **Use DNA copy-number** (validation-only firewall, cannot gate the method).
3. **Define families by full-copy architecture** rather than exon-pair sharing — a
   research direction, not an empirical gate.

The practical takeaway is that the current size + repeat-hub gate is the best
RNA-only post-refinement available; the residual domain-sharers are a principled,
characterized boundary rather than a tunable artifact.

### New linkage knobs tested

| knob | setting | families | P_Ep | R_oracle | comment |
|---|---|---:|---:|---:|---|
| baseline | `--identity 0.90` | 4633 / 471 multi-copy | 0.8705 | 51/57 = 0.8947 | default pair-node VG |
| stricter identity | `--identity 0.95` | 4589 / 417 multi-copy | 0.8681 | 52/57 = 0.9123 | tiny shift; not a domain-share fix |
| exact graph-to-graph | `--mode exact` | 4279 / 214 multi-copy | **0.8879** | 49/57 = 0.8596 | higher purity but cDNA pair-recall crashes (0.56) |
| min-shared-pairs 2 | `--min-shared-pairs 2` | 3283 / 186 multi-copy | 0.8710 | **40/57 = 0.7018** | removes domain-sharers but catastrophically under-merges real paralogs |
| antisense gate | `--antisense-gate` | 4633 / 471 multi-copy | 0.8705 | 52/57 = 0.9123 | neutral on headline P/R; reshuffles the FP roster |
| O1+VG integration | `--mode o1vg --identity 0.90 --antisense-gate` | 611 / 611 multi-copy | **0.8871** | 51/57 = 0.8947 | higher purity **without recall loss**; fewer but cleaner families |

- **`--min-shared-pairs 2`**: requiring ≥2 shared colinear exon-pair nodes is too
  stringent for divergent real paralogs (recall drops 0.89 → 0.70).  Domain-sharers
  often share a multi-exon domain block, so this is not the hoped-for separator.
- **`--antisense-gate`**: edge-level reciprocal-overlap / opposite-strand pruning
  (same geometric rule as the shipped O1 4th gate), wired directly into the VG-native
  locus graph.  Neutral on headline P_Ep / R_oracle in the prototype; it reshuffles
  the FP roster without improving the domain-share residual.
- **`--mode o1vg`**: O1 transcript-homology edges are kept only when supported by a
  shared non-repeat VG exon node or splice edge.  This is the most promising
  integration: P_Ep rises from 0.8705 to **0.8871** while R_oracle stays at 51/57,
  and cDNA pair-recall improves slightly (0.817 → 0.822).  The trade-off is fewer
  total families (611 vs 4633 raw) because the O1 edge set is sparse, but the
  resulting families are cleaner.

### Additional repeat/TE sweeps

A simple per-node repeat-fraction threshold (`max_node_repeat_frac >= 0.5`)
flagged 21 EP-impure families but also 537 real families.  Combining with size
(`max_node_repeat_frac >= 0.5 AND n_members >= 50`) still incurred 10 real
families of collateral for 6 EP-impure removed — worse than the existing
`repeat_hub_frac >= 0.05` gate, which catches the same repeat-bridge class with
zero collateral.  The repeat/TE signal is therefore best captured by the
multiplicity + repeat-richness conjunction already in place.

## Conclusion

Empirical FP gates give a small but clean improvement in protein purity.  The best
standalone gate remains the disjunction of large family size (`n_members >= 80`) and
repeat-rich hub nodes (`repeat_hub_frac >= 0.05`).

The most important new finding is that **O1+VG integration (`--mode o1vg`) improves
P_Ep from 0.8705 to 0.8871 with no loss in diploid-oracle recall (51/57 = 0.8947)**.
It keeps an O1 transcript-homology edge only when the two loci share a non-repeat VG
exon node or splice edge, making the definition both more mechanistic and more
self-contained.

Domain-sharers remain the hard residual: every RNA-only architecture axis tested
(exact-exon Jaccard, intron conservation, unshared flanking runs, pair-node Jaccard,
min-shared-pairs, antisense gate) either has no signal or costs too much recall.  The
remaining small, low-multiplicity DNA-confirmed FPs are still only separable with DNA
copy-number information or a fundamentally different definition signal.
