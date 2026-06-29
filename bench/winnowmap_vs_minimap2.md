# Would winnowmap (or a minimap2∪winnowmap merge) improve copy-assignability? — No (measured)

**Test.** Same 2,783 candidate reads (every read minimap2 placed over four hard tandem families on
NC_073247.2), re-aligned to the same NC_073247.2 sub-reference by **minimap2** vs **winnowmap** with
*identical* params (`-ax splice:hq -uf --eqx -Y -N50 -p0.1 --secondary=yes`; winnowmap adds its
repeat-weighted-minimizer DB `-W repetitive_k15.txt`, meryl k=15 distinct=0.9998, 15,985 repetitive
k-mers). Only the aligner differs. Assignability scored with the production engine
(`copy_assign.py::assign_family`). Families: DSFAM10 MAGEA (12 copies), DSFAM12 MAGEB (7), DSFAM42
(5 copies / **3 PSV columns** — near-identical), DSFAM43 (5 / 12).

**Result — winnowmap matches minimap2, never beats it:**

| family | aligner | reads | MQ0% | clip% | PSVc | resolv% | assigned% | silver agree |
|---|---|---|---|---|---|---|---|---|
| DSFAM10 (MAGEA) | minimap2 | 1311 | 37% | 0.1% | 2013 | 99.5 | 93.5 | 0.9939 |
| | winnowmap | 1310 | 37% | 0.2% | 2013 | 99.6 | 93.6 | 0.9939 |
| DSFAM12 (MAGEB) | minimap2 | 1199 | 0% | 0.1% | 1683 | 99.3 | 79.2 | 0.9926 |
| | winnowmap | 1200 | 0% | 0.2% | 1683 | 99.3 | 79.5 | 0.9895 |
| DSFAM42 (5cp/3PSV) | minimap2 | 253 | **95%** | 0.2% | 3 | 98.9 | 21.5 | 1.000 |
| | winnowmap | 247 | **95%** | 0.4% | 3 | 98.4 | 18.7 | 1.000 |
| DSFAM43 (5cp/12PSV) | minimap2 | 276 | 42% | 0.6% | 12 | 96.1 | 96.1 | 1.000 |
| | winnowmap | 274 | 42% | 0.4% | 12 | 96.1 | 93.4 | 1.000 |

**Read placement (recovery check):** of the 2,783 reads, minimap2 mapped 2,783, winnowmap 2,782 —
**winnowmap-only recoveries = 0, minimap2-only = 1.** Soft-clipping is negligible for *both* (0.1–0.6%),
so there is no "clipped around the PSV" problem for winnowmap to fix.

## Why (and it confirms the identifiability thesis)

1. **minimap2 `-N50 -p0.1` already saturates placement for long HiFi IsoSeq reads.** Every read is
   placed, with full PSV-column coverage and ~no clipping. A long read anchors fine even in a
   near-identical array, so winnowmap's repeat-aware placement has nothing to recover. (Soto's
   short-read 0.85% sensitivity in SD98 is a *short-read* phenomenon; long reads don't have it.)
2. **Our assignment uses read BASES, not MAPQ/AS.** minimap2's arbitrary MAPQ-0 primary among
   near-identical copies is irrelevant — the read's bases at PSV columns decide the copy regardless of
   which paralog it aligned to. So a "better primary" buys nothing.
3. **The residual hardness is the identifiability floor, which is aligner-invariant.** DSFAM42 sits at
   **95% MAPQ-0 under *both* aligners** because 5 copies share all but 3 PSV columns — that is
   information-theoretic non-separability, exactly the bound proven elsewhere. No aligner crosses it.

## Verdict

**No** — neither winnowmap nor a minimap2∪winnowmap merge improves copy-assignability here. A merge would
add only its **double-counting hazard** (our vote aggregation is global-per-read; the same PSV base
placed by both aligners would be counted twice, inflating the dominance margin and risking *false*
confidence) and an **AS-scale mismatch**, for zero placement upside. Keep the current minimap2
`-ax splice:hq -N50 -p0.1` BAM. The lever is the assignment gate (n_decisive≥1 + τ), not the aligner.

## Honest caveats

- One chromosome (NC_073247.2 / X), 4 families, 2,783 reads — not genome-wide.
- Per-chromosome meryl DB (the array repeats are local, so they are captured; a whole-genome DB would
  mostly add cross-chromosome repeats irrelevant to these same-chromosome tandem arrays).
- Tested reads minimap2 already placed near the loci; the explicit dropped-read recovery check
  (0 winnowmap-only) covers the "winnowmap rescues SD reads" rationale and finds it absent for long
  reads here. An extreme array not in this set could in principle differ, but DSFAM42 (3 PSV / 5 copies)
  is already near the hardest case and showed no benefit.

Artifacts: `winnowmap_vs_minimap2_assign.py` · `winnowmap_vs_minimap2_summary.json` ·
`/home/juanfra/winloci_scratch/win_test/{mm2,wm}.bam`.
