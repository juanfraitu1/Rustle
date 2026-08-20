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

## 5. What to do

1. **Rebuild the catalog at current defaults** (running now) and record both numbers.
2. **Re-measure the headline rates on it.** They are the thesis's numbers and they are currently
   attached to a catalog the tool no longer produces.
3. **State the provenance wherever a number is quoted** — which binary, which defaults, which date.
   The `refine` default changed under these numbers once already without anyone noticing.
4. Do **not** silently substitute the 526-family catalog for the 494 one. They are different objects;
   the difference (32 families that refine's "homology component AND ≥2 distinct loci" rule removed)
   should be characterised before either is called the catalog.
