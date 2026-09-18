# Pre-registration — NPIP algorithmic-ceiling simulation (canonical transcripts, distinguishable reads)

**Written 2026-09-18 BEFORE generating a single read.** Follows §6m9: the existing `npip_ideal` substrate
cannot exercise the assembler at all, because its 30 byte-identical reads per transcript collapse to ONE
witness per distinct (position, CIGAR) and every chain then falls below `pass1_min_reads`.

## What this substrate is, and the ONE thing it can answer

**Question:** given reads that cover every CANONICAL junction and are distinguishable, does the assembler
reach 26/26 complete NPIP chains?

⚠ **It is an ALGORITHMIC CEILING and nothing else.** Transcripts are derived from the same RefSeq
annotation that defines the truth, so — per `npip_ideal/DECLARATIONS.txt` C0, which binds here too — it
**cannot** measure how faithfully node construction recovers annotation. A high score is close to the
definition of the substrate. The informative outcome is a NEGATIVE: a failure here is an algorithmic defect
with no data excuse.

## Construction (fixed now)

1. **Transcripts = the CANONICAL chain of each of the 26 spliced NPIP copies.** Each copy's annotated exon
   blocks are merged across every NON-canonical junction (a non-canonical "intron" becomes retained
   sequence), so the transcript's junction set is exactly that copy's canonical junctions. Rationale
   (§6m9): 10 of the 40 non-canonical junctions are 1-29 bp — adjacent exon blocks split by a base or two,
   not introns. Simulating them would encode the artifact into the reads and then "recover" it.
   Target truth = **209 canonical junctions across 26 copies**, the same truth §6m7 onward uses.
2. **Reads are DISTINGUISHABLE.** 30 reads per transcript, each with independently jittered ends, so no two
   records share a (position, CIGAR) and the dedup cannot collapse them.
   - **Arm FL (full length):** start jittered 0-30 bp into the 5' end, end jittered 0-30 bp before the 3'
     end. Every read still spans every junction.
   - **Arm TR (truncated):** 3' end anchored (Iso-Seq is polyA-anchored), 5' start drawn uniformly over the
     transcript — the realistic failure mode.
3. **Error-free.** Sequence is taken verbatim from `chm13v2.0.fa` at the exon coordinates, reverse
   complemented for `-` strand copies.
4. **Alignment:** `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes` against the prebuilt
   whole-genome index `npip_ladder/idx/target.splice.mmi` — identical settings to both existing BAMs.
5. **Assembler arm:** `copy_assign --regions npip.regions --gtf --read-isoform-k 3` with
   `RUSTLE_JUNCTION_MAJORITY=1` (the §6m8 recommended arm). Scored exactly as §6m7 onward.

## Pre-registered decision rules

**SIM-0 (substrate sanity, checked FIRST).** Every one of the 26 transcripts must align back to its own
locus, and the aligned reads must between them carry **209/209 canonical junctions**. If the substrate
itself does not present every junction, nothing below is interpretable and the run is void.

**SIM-1 (the ceiling, arm FL).** With full-length distinguishable reads the assembler should reach
**>= 24/26 complete chains**. **A result below that is an ALGORITHMIC DEFECT with no data excuse** and is
the informative outcome. Real data currently gives 10/26 (11/26 at `GATE_MIN_READS=1`).

**SIM-2 (the cost of 5' truncation, arm TR).** Reported as the gap FL - TR in complete chains and in
canonical junctions. This isolates how much of the real-data deficit is truncation alone.

**SIM-3 (no circular claim).** Neither arm may be cited as evidence that node construction, the family
definition, or O1 works. Only as an upper bound on the assembler, and only for NPIP.

**SIM-4 (dedup guard).** Report the number of distinct (position, CIGAR) records vs total records per arm.
If the ratio is not ~1.0, the jitter failed and the run repeats SIM-0.
