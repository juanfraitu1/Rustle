# The K=0 wall is a property of the READ, not the COPY

**Date:** 2026-07-09. **Sim:** `bench/sim_k0_flank.py` (planted, deterministic). **Question (advisor):** for
exonically-identical copies (K=0), is there *absolutely nothing* in the BAM that determines a read's provenance?

## Design — perfect coverage isolates the wall

One chromosome, four co-located copies of one gene, **all reads full-length (perfect coverage)** so no read can
fail for lack of coverage:

| copies | exons | introns / 3′ flank | expectation |
|---|---|---|---|
| A, B | diverged (d=0.01) → **exonic PSVs** | diverged | assignable |
| **C, D** | **IDENTICAL (d=0)** → **no exonic PSV** | diverged (introns d=0.02, flank d=0.05) | **K=0** |

Two read classes per copy (60 each): **`_exon`** = the spliced transcript only; **`_flank`** = the transcript
plus that copy's own 3′ flank (readthrough). Aligned `minimap2 -ax splice -N 50` (secondaries, as in GGO_mm.bam).

## Result

Per **unique read** (placements deduplicated):

| true copy / class | reads | assigned | tied |
|---|---|---|---|
| A_exon | 60 | 60 (100%) | 0 |
| **C_exon** (identical exons, no flank) | 60 | **0** | **60 (100% tied)** |
| **D_exon** | 60 | **0** | **60 (100% tied)** |
| **C_flank** | 60 | **60 (100%)** | 0 |
| **D_flank** | 60 | **60 (100%)** | 0 |

**Separation is exact:** all 60 assigned `C_flank` → copy 1; all 60 `D_flank` → copy 2. Zero cross-assignment.

**Direct information-content test** (pipeline-independent: align only the read's flank portion to flank_C vs
flank_D): **40/40 = 100%** recover the true copy.

## Conclusions

1. **With perfect coverage, K=0 is the ONLY wall.** No read abstained for coverage reasons; the tied mass is
   exactly the exon-only reads of the exonically-identical copies. This isolates the two causes of abstention
   observed on real data — (i) K=0 (no distinguishing column exists) vs (ii) coverage (columns exist, reads
   don't reach them) — and shows (ii) vanishes under full-length reads.
2. **The wall belongs to the READ, not the COPY.** The *same* K=0 copies are perfectly separable from reads that
   carry flanking sequence. For a read confined to the identical exonic span, provenance is
   information-theoretically absent (100% tied, as it must be). For a read that reads through into the divergent
   flank, provenance is fully recoverable (100% assigned, 100% correct).
3. **The pipeline already exploits it.** C and D exist as *separate copies at all* only because the flank-bearing
   reads extended the de-novo transcripts into the divergent flank, yielding 56 PSV columns there. Nothing new
   had to be built — the flank became exonic in the de-novo model.

## Bearing on real GGO data — MEASURED, and the rescue is RARE

⚠ **Correction of an earlier error.** An initial audit reported "47.7% of reads carry ≥20 bp soft-clips." That
number counted **secondary alignments**, which are clipped *by construction* (a secondary is a partial
re-alignment of the same read, not extra sequence). Filtering to **primary** records:

| record type | records | ≥20 bp soft-clip |
|---|---|---|
| all (unfiltered — the wrong number) | 747 | 47.7% |
| **primary only (the read)** | 373 | **1.1%** |
| secondary | 373 | 94.1% (artifact) |
| supplementary | 1 | 100% (artifact) |

The BAM was built with **`minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`** — the **`-Y`**
flag makes secondary/supplementary alignments **soft-clipped and carry the full `SEQ`**, which is exactly why
they show a ~94% "clip" rate. They are alternative placements of the *same* read, not extra sequence. Primary
records always soft-clip their unaligned ends, so the **1.1% primary figure is unaffected by `-Y`.**
**Methodological rule: filter `-F 2308` before any per-read CIGAR/clip statistic.**

**Splitting the real abstained mass** (p6 run, per unique read, primary records):

| | flank-bearing (≥20 bp clip) → rescuable | exon-confined → provably lost |
|---|---|---|
| **TIED (abstained)** — 140 reads | **7 (5.0%)** | **133 (95.0%)** |
| assigned — 636 reads | 27 (4.2%) | 609 (95.8%) |

When a clip *is* present it is substantial (median 325 bp, max 1030 bp) — but it is rare.

**Conclusion: on Iso-Seq the K=0 wall holds in practice.** ~95% of tied reads are exon-confined and genuinely
carry no provenance. Only ~5% could be rescued by the flank mechanism the sim proves. This is biologically
expected: mature mRNA ends at the polyA site and does not read through into genomic flanks; only intron-retained
or readthrough transcripts would. The sim establishes the *mechanism* (the wall is per-read, and flank breaks
it); the real data establishes the *magnitude* (flank is rarely present, so the wall stands).

## ⚠ A THIRD cause of abstention, exposed by this sim: a MISSING COPY

**B was dropped from the family** — its reads map at MAPQ 60 (uniquely placed), so they form no de-tie conflict
edge. This is **not acceptable behavior**: a family is defined by **homology**, not by read ambiguity. B is a
2%-diverged copy of the same gene; it is a family member, and its reads being *easily assignable* is a virtue,
not grounds for exclusion. This is exactly the `E_c`-no-edge flaw the `E_r` homology-primary definition was
built to fix.

**The damage is concrete:** B's 98 reads were still pulled into the family and came back **100% "tied"** — not
because they are unassignable, but because **their true copy was never admitted to the copy set.**

So abstention has **three** causes, not two:

| cause | mechanism | fixable? |
|---|---|---|
| **K = 0** | no distinguishing column exists (copies exonically identical) | **No** — proven wall; needs DNA/longer reads |
| **Coverage** | columns exist, the read doesn't span them | Yes — longer reads / deeper sampling |
| **Missing copy** | the read's true copy was never admitted to the family (`E_c` conflict-graph definition drops uniquely-mappable members) | **Yes — a definition defect** |

**Verified defect (not sim-only):** `copy_assign` (O2) forms families via `read_conflict::conflict_edges` /
`conflict_families`. The `E_r` homology-primary definition (`detect_homology_catalog_genome_wide`,
`gamma_quasi_clique_partition`) is reachable **only from `gw_family_catalog`**; `copy_assign` never references
it. **O2 therefore still runs on the family definition O1 already replaced** — so uniquely-mappable copies
(globin problem, GSTM1/2/4/5) are systematically dropped from assignment, and their reads become spurious
"tied" mass. Wiring `E_r` into `copy_assign`'s family detection is the outstanding O1↔O2 harmony fix.

## Other honest caveats
- Planted flank divergence is 5% over 300 bp — generous. Recent duplicates may have near-identical flanks, in
  which case the flank carries less (or no) signal. The real-data soft-clip audit shows the sequence is *there*;
  whether it *discriminates* on GGO is the open measurement.
- Introns diverge too (d=0.02) but spliced reads carry no intron sequence (`N` skips it); only intron-**retained**
  reads would expose intron PSVs. Junction *positions* are already used (O3).

## Reproduce

```
python bench/sim_k0_flank.py          # builds genome + reads + BAM, prints the direct flank test
copy_assign --bam .../k0.bam --fasta .../k0.fasta --region <printed> --min-copies 2 --skip-poa-diagnostic --out .../k0
```
Related: `bench/READONLY_COPY_NUMBER.md`, `bench/em_consistency.md` (the K=0/SoftZone floor), `project_k0_frontier_unresolvable`.
