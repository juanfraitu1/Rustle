# Pre-registration: a known family, SMN1/SMN2: does O2's pick of the real primary give better loci than minimap2's, StringTie's and FLAIR's? (2026-10-20)

**User (2026-10-20):** "before moving to O3, can we ensure it works on a known gene family? One that has members that are
similar to each other by >90% but less than 100%, can we outperform minimap2 in picking which is the real primary or can we
reduce the pool from several multimappers to just a handful of more likely candidate reads that we can use to make new loci
and thus outperform stringtie and flair?"

## The family and why it was chosen

**SMN1 / SMN2** (pool family S19SMN2: `.0` = SMN2, chr5:70,808,816-70,837,687 −; `.1` = SMN1, chr5:71,381,817-71,409,805 +).
The two copies are 99.8% identical.

The biology gives a truth that does not depend on any aligner. SMN2 skips exon 7 in most transcripts (Δ7, its c.840C>T), and
SMN1 includes it. A read minimap2 gives to the wrong copy therefore mixes the two copies' exon-7 patterns.

**Not held out by family.** On HG002, which was development data for the corrector, its version-1 and version-2 calls on
SMN's tied reads were already seen: 53% → 99% against the donor assembly.

**What is new:**
- what that does to **loci**: StringTie and FLAIR run on each read's chosen alignment;
- the **exon-7 readout**;
- **A119b's native SMN reads**, never used: the A119b fragment sets cut other reads.

## Methods: which alignment is each read's primary

**Reads:** every read whose minimap2 primary lies on SMN1 or SMN2 (`summary.attribute`).

**Arms:**
- **mm2:** minimap2's primary.
- **O2:** for a read in the corrector's domain (tied within 10 AS on both copies, nothing off the family), the corrector's
  copy; otherwise minimap2's primary. The corrector is `o2correct.py donor`, version 2, commit bfa15f12.
- **truth:** HG002 only. The donor-assembly copy (`hg002proof` truth classes unique / unique_elsewhere), else minimap2's
  primary.

**Building each arm's BAM.** One record per read: its best-AS alignment on the chosen copy, given the read's sequence and
qualities and made primary. A read whose chosen copy has no alignment keeps minimap2's primary (counted).

**Assemblers**, on each arm's BAM (the lab's recipes, `bench/entangled/run.sh`):
- **StringTie** 3.0.3: `-L`, no `-G`;
- **FLAIR** 3.0.0: `bam2bed` + `collapse --trust_ends --generate_map` (the genome `chm13v2.0.fa`).

## Readout: exon-7 inclusion per copy

Exons are taken from the copy table (0-based, half-open):

| copy | exon 6 | exon 7 | exon 8 |
|---|---|---|---|
| SMN2 | 70,815,947-70,816,056 | 70,810,763-70,810,817 | 70,808,816-70,810,319 |
| SMN1 | 71,403,494-71,403,603 | 71,408,728-71,408,782 | 71,409,227-71,409,805 |

**Classifying a read or transcript on a copy.** It **spans** when its blocks overlap exon 6 and exon 8. Then it is:
- **FL** when a block overlaps exon 7 by at least 20 bp;
- **Δ7** otherwise.

**Per copy,** inclusion = FL / (FL + Δ7):
- reads: counted one each;
- StringTie: transcripts weighted by `cov`;
- FLAIR: isoforms weighted by their mapped reads.

## Bars

- **L1 (HG002, reads):** for SMN1 and for SMN2, O2's inclusion is closer to the truth arm's than mm2's is.
- **L2 (HG002, loci):** each assembler is compared with itself on the truth arm. On each copy, the assembler on O2's BAM must
  be no farther from that than the assembler on mm2's BAM, and strictly closer on at least 2 of the 4 (tool, copy) pairs.
- **L3 (A119b, held out, biology):**
  - SMN1's inclusion minus SMN2's is at least as large under O2 as under mm2, for the reads and for both assemblers;
  - under O2, SMN2's inclusion is below SMN1's.
  - The minimap2 coin toss is expected to shrink that gap.
  - **UNDERPOWERED** if fewer than 20 of the corrector's domain reads span exons 6 to 8.

**Verdict:**
- **PASS:** L1, L2 and L3 hold.
- **FAIL:** any of them fails.
- **UNDERPOWERED:** L3 lacks power; L1 and L2 are then reported alone.

**A PASS would claim:** on a known 99.8%-identical family, O2's pick of the real primary turns minimap2's coin-tossed SMN reads
into each copy's own isoforms, and StringTie and FLAIR built on O2's picks get closer to the donor's truth than on
minimap2's.

**It would not claim:** other families, or new copies.

**Rules:**
- The tool is `bench/o3_graph/o2loci.py`, written and tested before any arm is scored.
- Run script: `/mnt/linuxdisk/tmp/o3_graph/o2loci_smn/run.sh`.
- A bug found during the run is fixed and recorded as an amendment before any rerun.
