# PREREG — O3 (advisor's definition): PSV clustering of the reads O2 places on a family, held-out simulation (2026-10-15)

**Written after development on seeds 31, 32, 33 and before any held-out world is built.** Design:
`docs/superpowers/specs/2026-10-15-o3-unexplained-multimappers-design.md`. Frozen code (md5): `bench/o3_unexplained/psvclust.py`
e99d49b6bd3d9d1ad9f61ae19fa70c98, `score.py` 852aa244a00967c153abb6f2b9bc2cc1, `world.py` dfe160287b5fe50e07a81487d48eb237;
`copy_assign` built `--release` from 5022350c (sha1 87eed3f8); minimap2 2.30, samtools 1.22.1.

## Claim under test

Given the reads of a multi-copy family mapped to a reference that LACKS one of its expressed sequences (a copy, or an allele), O2
places those reads on the remaining copies; clustering each copy's reads by read-derived variant sites (psvclust.py) recovers the
missing sequence's reads as one cluster that no reference copy explains, does not produce such a cluster when nothing is missing,
separates two missing sequences absorbed onto one copy, and labels the cluster `copy` when its estimated divergence exceeds
delta = 0.00958 (else `allele_or_copy`: from RNA alone an allele and a copy are the same signal below delta).

## Method (frozen; psvclust.py docstring)

Per host class (copies with an identical spliced sequence are one class), the reads whose O2 best copy (`catalog_copy_idx`) is in
it and every family copy's spliced sequence are aligned to the class's padded locus; sites = positions where >= 3 reads carry the
same non-reference base; reads join the compatible cluster with the most co-covered sites; clusters < 3 reads dissolve; clusters
merge unless a site separates them with >= 3 reads on each side; a cluster is UNEXPLAINED when >= 2 of its alt sites are carried by
no family copy (and no copy explains all of them); est_div = unexplained alt sites / bases covered.

## Held-out worlds and arms

`world.py build OUT --seed S --strands` for **S = 41, 42, 43, 44, 45**; read variants **E0** and **E2**; `world.py run` with arm
**P** (`copy_assign --families --flag-missing-copies --no-as-tied-only`: every molecule placed on the family) — the primary arm — and
arm **T** (the shipped exact-tie gate) reported beside; then `psvclust.py` and `score.py rows`. 51 families per world.

## Bars (arm P, pooled over 5 seeds × 2 variants)

| id | measure | bar | dev (s32-33, arm P) |
|---|---|---|---|
| **D1** | absent copies with d >= 0.6% recovered (pure unexplained cluster: >= 80% of the reads O2 placed, purity >= 0.9) | >= 0.90 | 112/112 |
| **D2** | control families (NULL_TWIN, NULL_MULTI, NULL_CLOSE, ADDBACK, DECOY) with any unexplained cluster | <= 2 of 120 | 0/48 |
| **D3** | COMBO_EH and DOUBLE families with BOTH missing sequences recovered as separate clusters | >= 0.80 | 8/8 |
| **D4** | recovered alleles labelled `allele_or_copy`; recovered copies with d >= 2% labelled `copy` | >= 0.90 each | 28/28; 60/60 |

## Reported beside (no bar)

Copies at 0.3% (below the origin certificate's absorption floor) and at 1% (at delta: either label); every scenario; arm T on every
measure; the existing detector (`o3_flag == missing_copy` at the host copy) on D1 and D2 as the baseline; est_div against d.

## Reading, fixed now

- **D1-D4 pass** → PSV clustering of the reads O2 places on a family is the O3 method in its advisor's definition, at simulation
  level: it finds missing copies down to the allele range, does not flag intact families, separates two missing sequences, and
  states which ones only DNA can tell from an allele. Next: the same with real families on the main machine — remove real copies from
  the GENOME, re-map the real reads, add them back — and a single-copy-gene control for the real allele background (pre-registered
  separately). The method is then ported to Rust inside `copy_assign` under its own spec.
- **D1 or D2 fails** → the method is not adopted; the failing scenarios are named; nothing is re-tuned on the held-out worlds.
- **D3 or D4 fails alone** → detection is claimed, separation / labelling is not.
- Arm T is not a bar: on dev it sees only the exact-tie (twin) scenarios (16/16 there, 0 elsewhere), which is itself the finding to
  report — the absorbed reads of a missing copy are mostly UNIQUE mappers, not ties.

## Known limits, stated before the run

Simulated reads carry uniform HiFi errors, equal expression and one or two absent sequences per family; real data add gene
conversion, many-copy families, real error structure, and a heterozygous allele background at every expressed gene of a diploid
individual, which this method reports as `allele_or_copy` clusters — the real-data specificity of the `copy` label is therefore the
open question for the next level.

# OUTCOME (2026-10-15) — D1-D4 all PASS on the held-out worlds

`pixi run bash bench/o3_unexplained/heldout.sh ~/rustle_data/o3_unexplained/heldout` (seeds 41-45 `--strands`, E0 / E2, arms P and
T); frozen code unchanged (md5s as above). Arm P:

| id | measure | value | bar | verdict |
|---|---|---|---|---|
| D1 | absent copies >= 0.6% recovered as a pure unexplained cluster | **310 / 310 = 1.000** (existing detector 292/310) | >= 0.90 | ✅ PASS |
| D2 | control families with an unexplained cluster | **0 / 120** (existing detector 0/120) | <= 2 | ✅ PASS |
| D3 | two missing sequences on one host recovered separately | **40 / 40** | >= 0.80 | ✅ PASS |
| D4 | alleles labelled `allele_or_copy`; copies >= 2% labelled `copy` | **74 / 78 = 0.949; 170 / 170** | >= 0.90 each | ✅ PASS |

**Reported beside.**
- Below the origin certificate's absorption floor: copies at 0.3% recovered 40/40 (existing detector 4/40); alleles at 0.25% 18/20
  (the two misses, both seed 43, are the only unrecovered sources of the run).
- Labels at delta: copies at 1% split between `copy` and `allele_or_copy` (TWIN 10/10, UNIQ 14/6, MULTI 12/8, DOUBLE 6/14); the 4
  alleles labelled `copy` are all at 0.8%, next to delta. The label is a divergence estimate and is honest only away from delta.
- Low expression: copies at 8, 12 and 16 reads 30/30.
- **Arm T** (the shipped exact-tie gate, the advisor's literal "multi-mappers") recovers only the twin scenario (100/100 there, 0
  elsewhere; D1 80/310): the reads of a missing copy are exact ties only when the reference holds identical twins of it; otherwise they
  map UNIQUELY to the nearest copy with extra mismatches. O3 must therefore cluster every read O2 places on a family, not only ties.

**Reading (as fixed).** At simulation level, PSV clustering of the reads O2 places on a family is the O3 method in its advisor's
definition: it finds every missing copy from 0.3% upward, flags no intact family, separates two missing sequences on one copy, and
labels by divergence which ones RNA cannot tell from an allele. Next: real families on the main machine (remove real copies from the
GENOME, re-map the real reads, add them back) with a single-copy-gene control for the real allele background, pre-registered before
it runs; then the Rust port inside `copy_assign`.
