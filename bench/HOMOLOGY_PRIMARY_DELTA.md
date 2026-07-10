# `--homology-primary`: the E_r membership delta

**Date:** 2026-07-09. **Implements:** `docs/superpowers/specs/2026-07-09-o2-homology-primary-families-design.md`.
**Question:** does defining family *membership* by transcript homology (E_r) instead of by read ambiguity
(E_c) recover the copies that the conflict graph drops — and is it safe to make the default?

**Answer:** on the planted sim, decisively yes. On real GGO data, **no — not as a pure replacement.** Pure E_r
recovers nothing the conflict graph missed on four probe regions, and on one of them it *destroys* a family the
conflict graph found. The mechanism is proven; the default must not flip. Details below.

## 1. The sim — the mechanism is real and complete

`bench/sim_k0_flank.py`: four co-located copies of one gene. A, B have diverged exons (assignable). C, D are
**exonically identical** (K = 0) with divergent introns and 3′ flanks. Copy **B** maps uniquely (MAPQ 60), so it
forms no de-tie conflict edge — the exact shape of the drop bug.

```
copy_assign --bam k0.bam --fasta k0.fasta --region k0chr:2500-31284 --min-copies 2 --skip-poa-diagnostic [--homology-primary]
```

| | E_c (conflict) | E_r (`--homology-primary`) |
|---|---|---|
| copies detected | 3 | **4** (B admitted) |
| PSV columns | 56 | **102** |
| unique reads seen | 458 | **480** (22 of B's reads never even reached assignment) |
| assigned | 240 (52%) | **360 (75%)** |
| tied | 218 | **120** |

Per true copy (unique reads, `assigned` wins over `tied`):

| true copy / class | E_c | E_r |
|---|---|---|
| A_exon / A_flank | 60/60, 60/60 → copy 0 | 60/60, 60/60 → copy 0 |
| **B_exon** | **0/59 — 100% tied** | **60/60 → copy 1** |
| **B_flank** | **0/39 — 100% tied** | **60/60 → copy 1** |
| C_exon (K=0) | 0/60 tied | 0/60 tied |
| C_flank | 60/60 → copy 1 | 60/60 → copy 2 |
| D_exon (K=0) | 0/60 tied | 0/60 tied |
| D_flank | 60/60 → copy 2 | 60/60 → copy 3 |

Two things to note. **Zero cross-assignment** in either mode — enlarging the copy set did not corrupt any
existing call, even though the Bonferroni certificate tightened from α/2 to α/3. And under E_r the residual
tied mass is **exactly** C_exon + D_exon = 120 reads: the abstained set is now precisely the provable K = 0
wall, with the "missing copy" cause eliminated. That is the whole claim of the spec, demonstrated.

## 2. Real GGO data — the mechanism does not fire, and once it backfires

Four probe regions, `GGO_mm.bam` against `GGO.fasta`, `--min-copies 2`, both modes, run serially.

| region | E_c edges → families | E_r edges → families | assigned (E_c → E_r) | tied % (E_c → E_r) |
|---|---|---|---|---|
| `NC_073224.2:101578582-101607889` | 1 → 1 | 1 → 2 | 620 → 620 | 6.8% → 6.8% |
| `NC_073224.2:5565365-5721372` | 0 → 0 | 0 → 3 | 0 → 0 | — |
| `NC_073228.2:144930369-145083603` | 0 → 0 | **1 → 1** | 0 → 0 | — |
| `NC_073228.2:182473722-182663103` | 3 → 1 | **1 → 3** | **412 → 0** | 65.9% → — |

(`families` here is the count before `colocated_families`, which is unchanged by design; the "assigned" column
is what actually came out of assignment.)

Region 1 is identical in both modes: its two copies were already conflict-linked, so the oracles agree.
Region 2 is a true negative in both — E_r finds no homology edge either, so the three de-novo transcripts there
are genuinely unrelated. **E_r does not invent families**, which is the precision result we want.

Region 3 is the one place E_r beats E_c on real data: E_c finds **0** conflict edges, E_r finds **1** homology
edge and unions the pair. But the family is still dropped downstream by `colocated_families`. Both reps are on
`NC_073228.2` and the 153 kb region sits well inside the 5 Mb window with `min_copies = 2` satisfied, so by
elimination the splitter is **strand** — an antisense homologous pair, which the 4th FP gate deliberately
excludes. (This is an inference from the surviving filters, not a direct strand read-out; it should be
confirmed before being cited.) So the recovered edge is real and the drop is, most likely, correct behavior.

### Region 4 is a regression, and it is the finding that matters

E_c finds 3 conflict edges among 4 reps → 1 component → 2 co-located families (`CAFAM0`: 2 copies, both
*collapsed*, 816 reads; `CAFAM1`: 2 copies, 393 reads) → **1209 read assignments, 412 of them `assigned`.**

Pure E_r finds only **1** homology edge among the same 4 reps → 3 components → **0 families, 0 assignments.**

Two of E_c's three edges are rejected by the E_r homology gate (asm20 ≥ 0.80 ∪ sensitive ≥ 0.60, coverage ≥ 0.50).
So **`E_c ⊆ E_r` does not hold operationally** — the spec's "pure E_r, no union" decision assumed it did. There
are two live readings and this run does not distinguish them:

- **E_r is losing a real family.** The two reps are homologous below the gate (ancient paralogs, or exon-sum
  coverage < 0.50 because the de-novo transcripts are fragmentary), and 412 correct assignments are destroyed.
- **E_c's extra edges are repeat-bridge false positives.** This is precisely the failure the `refine` gate
  (id ≥ 0.80, cov ≥ 0.50) was built to kill — raw conflict-graph precision is 0.64 pre-refine. Under this
  reading, `CAFAM0` is a spurious family and its 412 assignments are assignments *into a fake copy set*. Its
  profile is suspicious: 816 reads, both copies *collapsed* rather than distinct reps, and 90% tied — this is
  the same outlier family that skewed the earlier tied-fraction estimate from 25% to 44.9%.

**The settling test:** dump the four reps' exon-sum sequences and compute pairwise identity and coverage for
the two E_c-only edges. If identity is high but coverage < 0.50, the gate's coverage term is the culprit and
E_r is under-calling. If identity is low and the shared span is a known repeat, E_c is over-calling and the
region-4 "regression" is actually E_r removing a false positive. Until that is run, region 4 must not be cited
in either direction.

## 3. What this changes

- **Do not flip the default.** The spec anticipated an assignment shift from a larger copy set; it did not
  anticipate E_r *shrinking* the copy set below E_c on real data. `--homology-primary` stays opt-in.
- **The "missing copy" cause of abstention is real and fixable** — proven end to end on the sim, where it
  accounts for 98 of 218 tied reads. Its **prevalence on Iso-Seq is not yet established**: on four real probe
  regions it fired once (region 3), and that family was correctly filtered on strand. Four regions is far too
  few to call the rate. This is the honest position, and it is weaker than the sim alone suggests.
- **"Pure E_r, no union with E_c" is now questionable.** The spec ruled a union out as a non-goal on the
  premise `E_c ⊆ E_r`. Region 4 falsifies that premise operationally. Whether the right object is `E_r ∪ E_c`,
  or `E_r` with a relaxed coverage term, depends on the settling test above — not on taste.
- The residual tied mass under E_r on the sim is **exactly** the K = 0 set. The three causes of abstention
  (K = 0 / coverage / missing copy) are now separable in a controlled setting, which is what makes the K = 0
  wall a measurement rather than an assumption.

## Reproduce

```
python bench/sim_k0_flank.py    # note: -N 50 is required; secondaries build the conflict graph
copy_assign --bam .../k0.bam --fasta .../k0.fasta --region k0chr:2500-31284 \
    --min-copies 2 --skip-poa-diagnostic --out delta_ec
copy_assign ... --homology-primary --out delta_er
```

Related: `bench/K0_FLANK_EXPERIMENT.md` (the three causes of abstention), `bench/RNA_FAMILY_HOMOLOGY_REFRAME.md`
(why E_c is an ambiguity oracle, not a homology oracle), `bench/family_definition_formal.md`.
