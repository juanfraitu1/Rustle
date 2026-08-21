# ⚠ The shipped catalog was built by a path that is no longer the default

> ⚠ **CATALOG PROVENANCE.** Figures quoted per **/494** (families) or **/1415** (copies) describe the
> **superseded** 2026-07-17 catalog, which **no invocation of the current binary reproduces** — it was
> built with `refine` on, a default removed on 2026-08-20. The current default emits **627 families /
> 2,019 copies**. Re-measure before quoting: see [`NUMBERS.md`](NUMBERS.md) and
> [`o1_catalog_provenance.md`](o1_catalog_provenance.md).

**Status 2026-08-19. Discovered while trying to confirm the orientation-guard default on real data.**

## 1. The finding

`GGO_gwcat.log`, the log of the catalog every O1 number in this project is measured against:

```text
[gw-catalog-homology] 79569 skeletons -> 12415 reps over 26 contigs
[gw-catalog-homology] 3042 E_r edges -> 526 families (>= 2 distinct loci)
[gw-catalog] exon-sum +sensitive refine: 526 raw families -> 494 refined
[gw-catalog] wrote 494 families (309 cross-chromosome)
```

**The shipped 494-family catalog went through `--refine`.** It was built on 2026-07-17, when refine was
still the default on the homology path. `refine_enabled(o1_homology = true, false, false)` now returns
**false** — refine has been off by default since 2026-08-09.

⟹ **Re-running the current binary at current defaults on the same BAM would emit 526 families, not
494** — a **32-family, 6.5%** difference, before any of today's changes are counted.

## 2. What this touches

Every O1 headline is measured on `GGO_gwcat.*`:

| number | measured on |
|---|---|
| false-merge **1.33%**, false-omission **5.6%** | the 494-family catalog |
| definitional-hole exposure ceiling **41/494 = 8.30%** | ” |
| **70.45%** of families are 2-copy; hierarchy reach **19.23%**; λ=1 shape **15.99%** | ” |
| the orientation guard's offline whole-genome analysis (2.22% of edges, 31 dissolutions) | ” |

None of these are *wrong* — they describe a real catalog. But they describe an artifact of a
**superseded default**, and a reader who reruns the tool will not reproduce the denominator.

## 3. ⚠ Two corrections to statements made earlier today

1. **"`--refine` is not the shipped O1 path — the catalog never calls it."** That is true of the
   *current* default and false of the *shipped catalog*. The λ-recomputation work therefore matters
   more than I said: had the shipped catalog been built after the certificate columns existed, all
   four would have printed `NA` for all 494 families.
2. **The offline orientation analysis modelled the wrong object.** It computed `E_r` over the 1,415
   *emitted copies*. The pipeline computes `E_r` over **12,415 reps → 526 families → refine → 494**.
   The guard acts at the 12,415-rep stage, and refine then re-clusters. So "31 families dissolve" is a
   prediction about a simplified model, not about the pipeline. **The T8 gap is wider than stated.**

## 4. Why it had not been caught

The genome-wide run has been unrunnable: it OOM'd at **23.7 GB RSS + 10.2 GB swap** and was killed at
1h34m. Nobody could rebuild the catalog, so nobody compared it to its own log. The streaming fix in
`37a4e01` removes that blocker — the buffered all-vs-all PAF was the cause.

## 5. ⭐ RESOLVED 2026-08-20 — the O1 catalog now has exactly ONE path

`--refine` and `--refine-introns` are **rejected** on the O1 homology catalog. They are not ignored and
they are not honoured: passing either is an error that explains why. The rationale is not that refine is
useless — it is that **a flag which changes what the catalog *is* is a provenance hazard**, and this one
proved it: it let the shipped 494-family catalog be built as `refine(γ-QC(E_r))` while the default
emitted `γ-QC(E_r)`, and the discrepancy survived six weeks.

The evidence says refine does not belong in the O1 definition anyway: it re-clusters by **connected
components** over its own substrates (no core-coverage denominator, no stub guard), which is not the
γ-quasi-clique(E_r) object `seeded_family_definition.md` §1 names, and under `--homology-genomic-span`
it split two recovered MAGEA copies into a spurious family and dropped GSTM4.

Refine keeps its real home — the legacy conflict catalogs (`--window-catalog` / `--cross-chrom`), where
it was written to run, is on by default, and was measured to help (APOBEC3 1f/2c → 0f/0c, SHARP
2f/4c → 1f/2c). `--no-refine` is still the escape hatch there.

⚠ **Consequence: the 494-family catalog is no longer reproducible by any invocation of the current
binary.** That is deliberate. It is not the defined object, and pretending it could be regenerated was
the hazard. Numbers measured on it must be re-measured, not re-derived.

## 5b. ⭐ THE GENOME-WIDE RUN COMPLETED (2026-08-20) — and two of my claims about it were wrong

First successful genome-wide catalog at current defaults: **2h18m, 627 families (445 cross-chromosome)**.

```text
94257 skeletons -> 17924 reps over 26 contigs      (Jul-17: 79569 -> 12415)
E_r edges by tier: sensitive=4778 (sole 4778)      (Jul-17: 3042)
16483 γ-quasi-clique blocks -> 627 families        (Jul-17: 526 raw -> 494 refined)
```

| | shipped Jul-17 | current defaults |
|---|---:|---:|
| families | 494 | **627** |
| 2-copy share | 348/494 = **0.7045** | 440/627 = **0.7018** |
| λ populated | (column did not exist) | **627/627 = 1.0000**, 75 `cut_certified` |

⭐ **The two-copy dominance REPLICATES** — 0.7018 vs 0.7045 on a catalog with 27% more families built
by a different path. That property is robust, not an artifact of the old catalog.

### ⚠⚠ Correction 1: I killed a run that was probably going to finish

The 2026-08-19 run was killed at 1h34m at **23.7 GB RSS + 10.2 GB swap, state `D`, ~11% CPU**, on the
judgement that it was swap-thrashing and would not complete. This run passed through **the same state**
— 24.46 GB RSS + 7.89 GB swap, state `D` — at ~2h07m and **finished 11 minutes later**. That phase is
the normal heavy end of the E_r stage, not a death spiral. The kill was wrong.

### ⚠⚠ Correction 2: the streamed PAF was not the unblocker

`37a4e01` was presented as the fix for the OOM. The completed run shows **our own process still reaches
24.46 GB**, and minimap2 separately reports a **24.479 GB** peak, on a 25 GB machine with 16 GB swap.
The buffered PAF was one allocation among much larger ones. Streaming is still correct — bounded beats
unbounded, and it removes one contributor at the exact moment memory is tightest — but it did not cause
the completion, and the E_r all-vs-all over 17,924 reps still runs **at the edge of the machine**.

**What is actually true:** the genome-wide catalog was always runnable; it takes ~2h20m and peaks near
25 GB. Anyone running it should expect swap and not interpret state `D` as failure.

## 6. What to do

1. ✅ **DONE 2026-08-20 — catalog rebuilt at current defaults: 627 families.** `o1_gw/ggo_gw.*`.
2. ⚠ **PARTLY DONE.** *Structural* properties are re-measured and agree with the offline estimates to
   ±0.019 (see §5b and `ONE_METHOD.md`). The *rates* — false-merge 1.33%, false-omission 5.6%, the
   8.30% exposure ceiling — are **NOT** re-measured: each needs truth labels (curation, excision), not
   just a recount. Until then every one must be quoted **with the 494-family catalog named**.
3. **State the provenance wherever a number is quoted** — which binary, which defaults, which date.
   The `refine` default changed under these numbers once already without anyone noticing.
4. Do **not** silently substitute the 526-family catalog for the 494 one. They are different objects;
   the difference (32 families that refine's "homology component AND ≥2 distinct loci" rule removed)
   should be characterised before either is called the catalog.
