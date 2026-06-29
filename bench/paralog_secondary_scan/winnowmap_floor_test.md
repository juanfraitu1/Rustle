# Can a better aligner break paralog ties? Empirical test: minimap2 vs Winnowmap on DAZ

**Question (advisor, 2026-06-08).** Could we improve copy identifiability by tuning
minimap2's levers or switching aligner (e.g. Winnowmap)?

**Answer: No aligner/lever can break a tie among reads already in hand** — a read whose
span is byte-identical between two copies matches both with NM=0; no scoring scheme
(identical alignments → identical score), no minimizer setting, no aligner can break it.
Levers/Winnowmap help a *different* axis — mappability into repeats — not discrimination.

**BUT for DAZ specifically the ties are not an intrinsic floor — they are 5′-truncated
reads (see refinement at end).** Full-length reads resolve the copies; the fix is better
library prep (full-length capture), not a different aligner. A true intrinsic floor applies
only to copies identical over the *entire* molecule.

## Two supporting results

### 1. The tie-zone is byte-identical and invariant to minimizer params (`/tmp/floor_test.py`)
| | longest EXACT-identical substring (tool-independent) | minimizer tie-zone k19w19 → **k11w1** (max sensitivity) |
|---|---|---|
| RABL2 | 210 bp | 196 → 201 bp |
| DAZ | **2875 bp** | 2891 → **2876 bp** |

The DAZ tie-zone stays ~2,876 bp even at k=11/w=1 (every k-mer a seed) and equals the
exact-identical substring — so it is sequence, not seeding.

### 2. Winnowmap resolves 0 of minimap2's tied DAZ reads
217 primary reads over the DAZ span aligned **separately** to `DAZ1.fa` and `DAZ3.fa`
(so AS to each copy is unambiguous); a read is *tied* if AS is identical on both.

| aligner | preset | tied / 217 | decisive DAZ1 / DAZ3 |
|---|---|---|---|
| minimap2 | `-ax map-hifi` | **31 (14%)** | 98 / 88 |
| Winnowmap | `-ax map-pb`, `-W` meryl k15 distinct=0.9998 | **32 (15%)** | 103 / 82 |

- Reads Winnowmap **resolved** that minimap2 tied: **0**
- Reads Winnowmap newly tied that minimap2 resolved: 1 (i.e. marginally *worse*)
- The same 31 reads tie under both tools.

### 3. Lowering `-p` (report secondaries down to 60% of primary) breaks 0 ties
Same DAZ reads vs combined `DAZ1+DAZ3` reference; `-p` is the secondary-to-primary
score-ratio threshold for *reporting* a secondary (default 0.8).

| `-p` | secondaries reported | tied (sec AS == prim AS) | new spillover surfaced (ratio < 0.80) |
|---|---|---|---|
| 0.8 | 207 | **31** | 43 |
| 0.6 | 220 | **31** | 53 |
| 0.3 | 224 | **31** | 53 |

Ties are invariant to `-p`: the tied reads already sit at score-ratio ≈ 1.0, so they are
reported at any threshold. Lowering `-p` only adds placements at 60–80% of the primary
score — reads that **decisively belong to their primary copy** (spillover, the
LOC115934290 failure mode). There is no band of genuine-but-suppressed copy reads,
because a read's true copy is its highest-scoring (primary) alignment, never a suppressed
low-ratio secondary. Net effect of `-p 0.6`: more phantom risk, zero identifiability gain.

### 4. Raising `-N` (max secondaries reported) — matters only for high-copy arrays, still no ties broken
`-N` caps the number of secondary alignments reported (default 5).

| family | `-N 5` | `-N 500` |
|---|---|---|
| DAZ (2 copies) | 31 tied, ≤4 sec/read | identical (31 tied) |
| n=24 array (NC_073244.2) | copies/read median 4, max 6 | copies/read median 13, max 24 |
| n=24 **best-AS tie fan-out** | median 1, **max 6** | median 1, **max 6** (unchanged) |

`-N 500` surfaces many more placements in high-copy arrays (real, unlike `-p`), but every
added placement is lower-scoring spillover — the **best-AS tie set is unchanged**. Useful
only as *completeness of the tie set for EM apportionment* when a tie set exceeds ~6
copies; it is input quality, not resolving power. Breaks zero ties.

## Summary — four levers, none breaks the floor
| lever | moves | breaks ties? |
|---|---|---|
| minimizer k/w | seeding sensitivity | no (invariant to k11/w1) |
| Winnowmap | mappability into repeats | no (0 of 31 resolved) |
| `-p` ↓ | reports lower-ratio secondaries | no (spillover) |
| `-N` ↑ | reports more secondaries (high-copy) | no (completes tie set / spillover) |

## Where levers/Winnowmap DO help (mappability, not discrimination)
- `-N` / `-p`: emit secondaries so the cross-copy comparison happens at all.
- `-f` / Winnowmap weighted minimizers: map reads **into** highly repetitive segmental
  duplications that minimap2 hard-filters to MAPQ 0 (a recall gain for big n≈24–32
  paralog clusters) — but still cannot distinguish identical copies.
- scoring `-B`: sharpens the **near-tied** margin (~84 DAZ reads at AS-ratio 0.95–0.999),
  never the ~31–39 byte-identical ones.

## Does the T2T reference help? Yes — representation, not discrimination (and it strengthens the floor)
The assembly is T2T-grade: DAZ region and all of chrY have **0 Ns** (gapless), and DAZ1 /
DAZ3 are represented as **two separate complete loci** (un-collapsed). T2T's real
contribution is making the problem *well-posed*: a draft assembly could have collapsed the
copies into one (no resolution possible) or left gaps (reads unmappable). But T2T
faithfully represents the *true* sequence — the 2,875 bp identical core is measured **on
the T2T assembly itself**, so the floor is confirmed REAL biology, not a collapse/gap/error
artifact. (A *less* accurate assembly could even look more distinguishable, via random
errors acting as false PSVs; T2T removes those, so discrimination uses only real PSVs.)
Better reference ⇒ more honest, not more identifiable.

## What actually breaks the floor
Longer reads that anchor in a **discriminative flank** and span *through* the identical
core (this is why DAZ's ~174 decisive reads resolve while the core-trapped ~31–39 never
can). The aligner is already optimal; the limit is in the data.

Reproduce: `/tmp/wmtest/` (daz_reads.fa, DAZ1.fa, DAZ3.fa, mm_*.sam, wm_*.sam).

---

## IMPORTANT REFINEMENT: DAZ ties are 5′-truncated reads, not an intrinsic floor (2026-06-08)
The reads are IsoSeq FLNC (meant to be full-length), so a read covering the whole
transcript should hit any divergent exon. Measured: they don't, because the tied ones are
**truncated**.

| read class | median read len | median exonic-aligned | reaches 5′ end of transcript |
|---|---|---|---|
| TIED (n=41) | 2,410 bp | 2,337 bp | 3 / 41 |
| DECISIVE (n=176) | 4,224 bp | 4,153 bp | 133 / 174 |

DAZ1/DAZ3 **diverge at the 5′ end**; the tied reads are systematically **5′-truncated**
(have 3′/polyA, RT didn't reach the 5′ cap), so they fall in the identical 3′ core. The
full-length reads reach the 5′ divergent exons and **are decisive** — 176 of them resolve
the copy. So DAZ's copies ARE distinguishable; the ~19% ties are a **library-prep
(5′-truncation) artifact**, not intrinsic copy identity.

Correction to the earlier framing: full-length HiFi (~4.2 kb) DOES exceed DAZ's 2.9 kb
identical core; only truncated reads don't. The genuine information floor is reserved for
copies identical over the *entire* transcript (no distinguishing position anywhere) — DAZ
is not that.

**Can we complete the reads?** Computationally on existing reads: no (a 5′-truncated read
lacks the divergent 5′ sequence; recovering it is circular). Prospectively: yes —
**5′-cap-selected full-length cDNA** (TeloPrime / Cap-trapping / SMARTer) yields complete
reads that reach the 5′ end → decisive; this would also settle whether DAZ3 is expressed
or a phantom. The four aligner levers still can't help the truncated reads already in hand
(can't align sequence that wasn't read) — the lever here is read COMPLETENESS (data).
