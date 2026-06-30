# Performance And Io (consolidated)

> Merged from 3 source docs (verbatim, git keeps the originals' history). Each section below was a separate `bench/*.md`.

**Contents:** [PERFORMANCE_OPTIMIZATION](#performance-optimization) · [input_formats_and_ties](#input-formats-and-ties) · [winnowmap_vs_minimap2](#winnowmap-vs-minimap2)


---

## PERFORMANCE_OPTIMIZATION

# Copy-assignment pipeline speedup (2026-06-27)

Goal: make the genome-wide copy-assignment sweep faster **without changing a single assignment** (the O2
headline is byte-frozen). Every change below was verified output-identical by re-running real regions and
diffing the assignments + abundance (`diff` line count 0), plus the 659-test suite stays green.

## Where the time actually went (measured, not guessed)

Profiled the heaviest real family (GAGE, `NC_073228.2:144928569-145085403`, 20,067 mapped reads → 8,461
assignments, 5 copies, 3,238 PSV columns) with env-gated stage timers (`RUSTLE_TIMING=1`):

| stage | time | note |
|-------|------|------|
| `reads_in_region` (BAM I/O, 20k reads) | 0.6 s | cheap |
| pass1 + gate + collapse | 0.0 s | cheap |
| **`build_read_placements` (20k × 10 reps)** | 0.0 s | the audit feared O(B·S); it's negligible |
| **POA homology *diagnostic*** | **~458 s** | dominant — and *purely diagnostic* |
| **`discover_psvs` (poasta, 4 pairwise DP)** | **69.7 s** | the real compute after the diagnostic is gone |
| per-read assignment (8,461 reads) | 1.1 s | already parallel; was never the bottleneck |

The original audit's asymptotic guess (per-read assignment is the O(R·C·N) hot path) was **wrong in
practice**: assignment is ~1 s. The wall-clock is the two alignment stages.

## The changes (all output-identical)

1. **Skip the POA homology *diagnostic* on sweeps — `--skip-poa-diagnostic`.** Families come from the de-tie
   conflict graph; the POA pairwise pass only writes `.fallback.tsv` and the diagnostic edge counts. It was
   left ON in the definitive run and cost ~85% of wall-clock. Skipping it: **537 s → 79 s (6.8×)**,
   assignments + abundance byte-identical. (Wraps the existing `RUSTLE_SKIP_POA_DIAGNOSTIC` env as a real,
   documented flag.)
2. **Parallelize `discover_psvs`** (`copy_assign_pipeline.rs`). The (n−1) ref-vs-copy poasta alignments are
   independent; run them concurrently and merge (per-index amaps + set-union of diffs = identical result).
   poasta is a pure function (thread-safe); the opt-in minimap2 path stays serial (its temp-file nonce can
   collide between equal-length copies under threads). **79 s → 43 s**, peak RSS 1.66 GB (safe on the 19 GB
   box), diff 0.
3. **Parallelize the per-read assignment loop** + **iterate only spanned PSV columns** in the significance
   gate and the editing pre-pass (`copy_assign.rs`, `copy_assign_pipeline.rs`). Reads are independent;
   rayon's indexed collect preserves order. The gate previously re-scanned all 3,238 columns per competitor
   per read — now it walks only the ~tens a read spans. Marginal on the heavy region (assignment is ~1 s
   there) but real on small/dense families and at cluster scale.
4. **`BamIndexCache`** (`denovo_assemble.rs`): parse the multi-MB `.bai` + header ONCE and reuse across all
   regions instead of per-region. ~free on 79 regions; compounds over thousands at genome scale.

**Heaviest region: 537 s → 43 s = 12.5×, byte-identical.** End-to-end 79-region sweep (`o2_verify`):
**3,854 s vs the old 8,130 s = 2.1×**, and the full result is **byte-identical** — all 74 families / 104,147
assignments diff = 0, abundance diff = 0, silver 24,660/24,682 = 99.9% unchanged. (The 2.1× is a *conservative
lower bound*: that sweep ran contended with the A\*/old-binary benchmark runs competing for the 5 cores; the
clean per-region 12.5× is the reliable figure. A clean re-run would land between the two.)

## Also tested, gated OFF (changes output → revalidate-grade, not a free speedup)

**Exact-A\* (`RUSTLE_POA_ASTAR=1`)** swaps poasta's `AffineDijkstra` for `AffineMinGapCost` (same optimal
*cost*, admissible heuristic, fewer explored states). Faster, but A\* returns a *co-optimal-but-different*
traceback: on the heavy family it changed **1,254 / 8,461 assignments** and found 3,239 vs 3,238 columns. So
it is **not** byte-identical and is gated **off** by default (the O2 result stays exactly poasta-Dijkstra). The
same co-optimal-traceback hazard applies to banded/WFA aligners — so the only *truly* output-identical
`discover_psvs` speedup is the parallelization above; going faster means accepting a small, revalidate-grade
output change (benchmark vs the poasta 3,238-col / silver-99.9 headline before adopting).

## Rejected (changes output)

**minimap2 PSV discovery** (`RUSTLE_PSV_MINIMAP2=1`) is ~270× faster (43 s → 2 s) but is a heuristic,
coarser aligner: on this family it found **109 PSV columns vs poasta's 3,238**, and **all 8,461 assignments
changed**. It trades accuracy for speed and would move the O2 headline, so it is NOT used for the frozen
result. (It remains a valid opt-in for a fast/approximate pass.)

## Recommendation

Genome-wide / cluster sweeps: run `copy_assign … --skip-poa-diagnostic`. The parallel `discover_psvs` and
per-read assignment are automatic. Keep poasta (default) for the PSV columns — it is the accurate engine the
O2 numbers are computed with. Use `RUSTLE_TIMING=1` to re-profile if a region is unexpectedly slow.

## Region-parallel sweep — `--region-threads N` (2026-06-29, byte-identical, opt-in)

After `--skip-poa-diagnostic`, the dominant remaining cost is the **exact poasta DP in `discover_psvs`** (~37 s
for the 5-copy/3,238-column GAGE family; ~64 s for a 2-copy/3,874-column chr19 family). That DP is irreducible
for byte-identical output (the heuristic minimap2 PSV engine changes the columns) and its `(n−1)` per-copy
alignments are *already* internally parallel — so a *single* heavy family cannot be sped up further.

The untapped axis is **across families**: the sweep processed regions serially, so on a many-core box a 2-copy
family ran one alignment on one core while the rest idled. `--region-threads N` (default 1 = the exact serial
path) processes a contig's regions in **chunks of N concurrently** — each family is independent
(`detect_and_assign` is pure; `BamIndexCache` opens a fresh reader per call and is `Sync`; the genome is
read-only). Output is collected and `CAFAM` ids assigned **in region order afterward**, so the result is
**byte-identical** (verified: serial vs `--region-threads 3/6` → all of families/assignments/quant/posterior/
the 3 PSV dumps diff = 0). Chunking bounds peak memory to ~N regions' reads (the documented genome-wide OOM is
why it is opt-in). The temp files of the opt-in minimap2 PSV path were made process-unique (atomic nonce) so
all alignment helpers are thread-safe.

**CROSS-CONTIG (the shipped design).** Regions are FLATTENED across all contigs and processed out-of-order in
one pool, so the globally-heaviest families — which live on *different* contigs — overlap (a single-family
contig no longer wastes the pool). A **bounded LRU cache of loaded contig genomes** (`lru`, capacity ≈
`region_threads`) lets any worker reuse an already-loaded chromosome and caps resident genome sequences (the
memory bound; `Arc` keeps a genome alive while an evicting worker still holds it). The heavy read SEQUENCES are
dropped inside each worker — the output stage needs only read NAMES (everything else is in the computed
assignments) — so collecting every region's result out-of-order stays lightweight. Output is still drained in
the original flat order, so `CAFAM` ids and all rows are byte-identical.

**The ceiling is the single globally-heaviest family** (≈ one 80 s poasta DP here): no region-parallelism can
split it (it is already internally parallel). So `parallel_wall ≈ max(heaviest_family, total_work / N)`:
- small/dominated sets (one giant family) → modest ratio (a 3-family chr19 set and an 8-contig set were
  bounded near the lone heavy family);
- the **full genome-wide sweep** (~3854 s total, heaviest family ≈ 80 s) → `total/N` dominates, so it
  approaches **~N×** up to the heaviest-family floor — the real win.

**Verified byte-identical** at every step: serial vs `--region-threads 3/6` on chr19, an 11-family contig, and
an **8-DISTINCT-contig** set → families / assignments / quant / posterior / mosaic / all 3 PSV dumps diff = 0.
Use the largest `N` your memory allows on a heavy sweep (peak ≈ `N` regions' reads + ≈`N` contig genomes);
keep `N=1` (the exact serial path) when memory-bound. `--dump-psv` at genome scale accumulates all regions'
PSV matrices in memory — prefer smaller `N` or per-region runs for that opt-in mode.


---

## input_formats_and_ties

# Input formats (SAM/BAM/CRAM) + aligner ties

## Native SAM / BAM / CRAM input (SHIPPED)

rustle now accepts **SAM, BAM, and CRAM** input natively (pure noodles — no samtools/htslib needed).

* **Detection** is by magic bytes (`bam.rs::detect_format`): `CRAM` container magic → CRAM; BGZF/gzip
  (`0x1f 0x8b`) with first decompressed bytes `BAM\1` → BAM, else bgzipped SAM; plain text → SAM.
* **BAM** uses the existing fast multithreaded BGZF reader, **completely unchanged** (the hot path —
  BAM decode is ~40 % of wall-clock — is untouched).
* **SAM / CRAM** are transcoded **once** to a temporary BAM at startup (`bam.rs::transcode_to_temp_bam`)
  and the unchanged multi-pass pipeline runs on it. This is far less invasive than re-plumbing every
  read loop, and the pipeline reads the alignments several times, so one up-front transcode is also
  cheaper than repeatedly streaming a slow SAM/CRAM reader. The temp BAM is deleted on exit
  (`TempBam` drop guard).
* **CRAM** needs an indexed reference FASTA: `--cram-ref <ref.fa>` (defaults to `--genome-fasta`); the
  FASTA must have a sibling `.fai`. Missing reference → a clear, actionable error.

Why this is safe: BAM is a lossless container of the same alignment records, so the assembly is
**identical** to running on an equivalent BAM.

**Verified (chr19 / NC_073243.2):** running rustle on `GGO_19.bam`, a SAM of it, and a CRAM of it
(`--cram-ref GGO.fasta`) produced **byte-identical GTFs** (BAM-vs-SAM and BAM-vs-CRAM = 0 differing
lines; 2006 transcripts / 23565 exons each). Temp files auto-cleaned; 618-test suite green.

**Operational caveats (SAM/CRAM only — BAM is untouched):**
- The transcode writes a temporary BAM **the same size class as the input** (a multi-GB SAM → a
  multi-GB temp BAM); its size is printed to stderr. It goes to `$TMPDIR` (falls back to `/tmp`) — point
  `$TMPDIR` at a roomy disk for whole-genome SAM/CRAM.
- The transcode is currently **single-threaded** (the multi-pass assembly that follows is parallel; only
  the one-time up-front transcode is serial). For a large whole-genome CRAM this adds a few minutes.
- The `TempBam` drop guard deletes the temp on normal exit and on errors, but **not on `SIGKILL`/OOM** —
  a hard kill can leave a `rustle_input_<pid>_*.bam` behind; clean `$TMPDIR` if a run is force-killed.
- Byte-identity is verified on chr19 only; it holds by construction (BAM losslessly carries the same
  records), but a whole-genome equality check has not been run.

```bash
rustle -L reads.bam  -o out.gtf                       # BAM   (fast path)
rustle -L reads.sam  -o out.gtf                       # SAM   (native)
rustle -L reads.cram --cram-ref ref.fa -o out.gtf     # CRAM  (native; ref must be .fai-indexed)
```

## Can we stop minimap2 / winnowmap from making ties?

**Short answer: no — and that's a feature, not a gap.** A genuine tie means the read is
*information-theoretically identical* across the copies over its mapped span (no PSV in the read → no
distinguishing column). That is the identifiability floor (the `n_decisive = 0` / K-frontier case): no
aligner score, `--eqx`, scoring matrix, or winnowmap weighting can break it, because there is no
sequence difference to score. Forcing the aligner to "pick one" doesn't remove the tie — it **hides**
it behind an arbitrary primary flag, which is exactly the failure mode to avoid.

What you **can** (and should) do — expose ties and resolve them, abstain only at the floor:

1. **Expose all tied/near-tied placements** instead of hiding them:
   `minimap2 -ax map-hifi --secondary=yes -N <K> -p <low> --eqx` (Eichler used `-p0.5 -N20`; for full
   tie exposure `-N50 -p0.1`). Then every co-/near-optimal hit is in the BAM with its `AS`, and rustle's
   PSV layer adjudicates.
2. **Don't trust the arbitrary primary among MAPQ-0 paralogs** — rustle already doesn't. Copy
   assignment uses the **AS / decisive-margin gate** (`psv_linkage::assign_read_to_copy`), which
   *abstains on a tie* ("no copy beats the runner-up by the margin → `None`"), the PSV-space restatement
   of Eichler's AS ≥ 10, with **no 1/k** splitting.
3. **Reduce the *breakable* (spurious) ties**: `--eqx` for base-level PSVs; soft-mask repeats to cut
   repeat-driven multimaps; a merged minimap2(`-N50 -p0.1`)+winnowmap BAM maximizes exposed candidates
   in segdups (winnowmap changes *which* copy is primary / reduces mismapping, but identical-copy ties
   remain). Longer reads / less 5′ degradation span more PSVs → fewer ties (an input property, not an
   aligner knob).

So the lever is never "make the aligner not tie"; it is "expose the ties, resolve with PSVs, and
**abstain** exactly when the molecule carries no distinguishing information" — which is precisely the
thesis's identifiability boundary (and Canzar-aligned: resolve ambiguity, no arbitrary 1/k).


---

## winnowmap_vs_minimap2

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

