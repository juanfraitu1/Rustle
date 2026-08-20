# ⚠ The shipped catalog was built by a path that is no longer the default

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

## 6. What to do

1. **Rebuild the catalog at current defaults** (running now) and record both numbers.
2. **Re-measure the headline rates on it.** They are the thesis's numbers and they are currently
   attached to a catalog the tool no longer produces.
3. **State the provenance wherever a number is quoted** — which binary, which defaults, which date.
   The `refine` default changed under these numbers once already without anyone noticing.
4. Do **not** silently substitute the 526-family catalog for the 494 one. They are different objects;
   the difference (32 families that refine's "homology component AND ≥2 distinct loci" rule removed)
   should be characterised before either is called the catalog.
