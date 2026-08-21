# HANDOFF — state as of 2026-08-19

> ⚠ **CATALOG PROVENANCE.** Figures quoted per **/494** (families) or **/1415** (copies) describe the
> **superseded** 2026-07-17 catalog, which **no invocation of the current binary reproduces** — it was
> built with `refine` on, a default removed on 2026-08-20. The current default emits **627 families /
> 2,019 copies**. Re-measure before quoting: see [`NUMBERS.md`](NUMBERS.md) and
> [`o1_catalog_provenance.md`](o1_catalog_provenance.md).

Supersedes the 2026-08-15 handoff (which recorded HEAD `4285b8f`, 0 unpushed — both now stale).

Branch `dna-from-genome`, HEAD **`1247d67`**, **7 unpushed commits**. Only `tools/stringtie`
(third-party submodule) is dirty — leave it. Source changes in `gw_family_catalog.rs` /
`denovo_pipeline.rs` (the λ certificate work) remain **uncommitted** — see §4.

## 1. What changed on 2026-08-19

Seven commits, all documentation and all negative-or-bounded. **The O1 definition is UNCHANGED and
no default moved.**

| commit | finding |
|---|---|
| `e8998a6` | Block-aware provenance model **assessed**: keep the typed representation, **defer the rooting**. `repeat_multiplicity` renamed `n_occurrences`, made descriptive-only (it was R5 by another name). Rooting pilot certified **3/18** family probes; its substrate is inverted (human-study/gorilla-outgroup). Hierarchy ceiling **29.55%**, real reach **19.23%**; all 79 λ=1 families are core-plus-one-satellite = **15.99%** of the catalog. |
| `3ca4713` | ⭐ **A P1-safe repeat gate.** R5's defect was the **universe, not the statistic**. `min_shared_gmult` counts a shared 21-mer's occurrences in **GGO.fasta**, so it depends only on `(a,b,genome)`. Seed-invariance: **0/147** changes vs the catalog-counted analogue's **94/147**. FP median 182 vs TP median 2, AUC 0.9429; at M=50 **10/12 FPs at 0/135 TP cost**. |
| `3e413f8` | It **cannot replace** coverage — no operating point. Parity with γ's 253 cross-family edges needs M≈3 and discards **48%** of the shipped edge set. **It is a veto, never an admission criterion.** |
| `4c9fc70` | Full-length-read **edge** guard: no separation. Blocks shorter than a median FLNC read (2,772 bp) are FP **14/14** but TP **130/144 = 0.9028**. |
| `2153561` | Full-length-read **node** guard: **inert**, 1,412/1,415 pass. Cause: the catalog has **0/1415 single-exon nodes**, so census pathology (b) cannot occur — a min-2-exon rule already does the job. ⟹ the 47 node-construction failures are dominated by pathology (a), one-locus-cut-in-two, which is a coordinate signature, not a read signature. |
| `ed86742` | Pre-registration for the haplotype copy-number run. |
| `1247d67` | ⚠ Run is **UNINFORMATIVE** by its own criterion (control floor 0.1512 vs signal 0.0278). **And it puts a prior number in doubt — see §2.** |

## 2. ⚠⚠ READ THIS BEFORE QUOTING ROW 3.10

`OBJECTIVES_AND_VERIFICATION.md` row 3.10 reports MAPKBP1/PLA2G4B/SPTBN5 at **8/9/9** (`_pri`/`_pat`)
vs **5/6/8** (`_mat`), and that "several-copy difference between one animal's two haplotypes" has been
the standing refutation of *"copy number is stable; only SNPs and indels differ."*

**2026-08-19 cannot reproduce it.** At minimap2's default `-p 0.8` MAPKBP1 returns **1/1** (its
paralogues score under 0.8× a 58 kb perfect self-hit and are discarded); at `-p 0.1` it returns
**9/8**. Neither is 8/5. SPTBN5 reproduces exactly; PLA2G4B is off by one. With secondaries retained
the mat deficits shrink from **3, 3, 1** to about **1, 2, 0**.

⟹ **Do not quote "1–3 whole gene copies differ between KB3781's haplotypes"** until row 3.10's
instrument is re-derived with its `-p`/`-N` settings recorded. Direction survives; magnitude does not.
**Copy counting is far more sensitive to secondary-alignment filtering than any prior doc states.**

## 3. Where O1 stands

The definition is unchanged and **P1 remains a theorem**. Measured: false-merge **1.33%**, false-omission
**5.6%**, identity-clause failures **0/728**, reach ~0.55. One named hole — the **scale-free min-length
coverage denominator**, exposure ceiling **41/494 = 8.30%**, and 30 is a floor. **Five** repair routes
are now closed (see `NEGATIVE_RESULTS_REGISTER.md`, 569 entries). **Two guards exist, both pair-local,
neither through the shipped binary.** Canonical statement: `ONE_METHOD.md`.

## 4. Open items, in priority order

1. **Whole-genome GGO/HSA comparison for the transcript-orientation guard** — the single thing blocking
   the strongest available FP reduction (29/74 for 4/9,032). `o1_false_positive_rules.md`.
2. **Held-out FP set** — without it, the genome-anchored gate's 10/12 is a description of a known set,
   not a rate. The 0/135 TP cost is the load-bearing half and does not need it.
3. **Commit or revert the λ work** — `gw_family_catalog.rs` and `denovo_pipeline.rs` are modified and
   uncommitted; `families.tsv` gains `n_edges/density/lambda/cut_certified`. Note the offline
   re-derivation in `/mnt/linuxdisk/home/juanfraitu/o1_lambda/offline/` already answered the reach
   question (control: 145/146 connected), so a genome-wide re-run is **not** required for that.
   ⚠ A previous attempt at that re-run ballooned to **23.7 GB RSS with 10.2 GB in swap** and had to be
   killed; it is single-threaded despite `--threads 4`.
4. **Re-derive row 3.10's instrument** (§2).
5. **Composition-matched controls** if the haplotype CNV run is retried — span-matching alone put the
   floor 5.4× above the signal because random intervals are far more repeat-rich than gene bodies.
6. Older, still open: fibroblast O1 replication (`/mnt/linuxdisk/home/juanfraitu/o1_replicate/`,
   ⚠ fix `build.sh`'s stderr/stdout redirect first); the Illumina/HiFi S2 depth check; sweep `c` at
   0.45/0.55/0.75; audit the 30/194 overlapping-copy cases.

## 5. Scratch directories from this session

```text
/mnt/linuxdisk/home/juanfraitu/o1_lambda/offline/   reach.py, cutshape.py, reach_{pre,post}_M1.tsv
/mnt/linuxdisk/home/juanfraitu/o1_gmult/            blocks.py, gmult.py, eval.py, covfree.py, fullread.py
/mnt/linuxdisk/home/juanfraitu/o3_hapcnv/           probes.py, align.sh, count.py, pri_provenance.tsv
```

⚠ **Crash rule still applies**: heavy runs one at a time (5 cores / 25 GB), never `pkill -f`, big
outputs to `/mnt/linuxdisk` (`winloci_scratch` free space is a sparse-vhdx fiction).
