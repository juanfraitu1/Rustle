# PREREG row 1711 (2026-11-14): a blind exon-level detector for exon gains

User go 2026-11-14 ("Pre-register a new detector"), after row 1694: `copy_phasing` fires on 2.3% of the gaining exon
units reads can see. Written before any detector code ran on the units. Seen so far: row 1694's aggregate tables only.

## Setting (stated because it limits the reading)
Units, reads and truth are row 1694's rebuild (workspace `r1694_exgain`): 15,469 exon units, 469 gaining (a KB3781
haplotype holds more copies of the exon than the reference). The reference `GGO.fa` is the same individual's primary
assembly (mGorGor1-v2.0_pri), so non-reference exon sequences in the reads are mostly this individual's other alleles.
Row 1694's oracle count is the ceiling for plain novelty: reads nearer a haplotype-only sequence in 45.4% of gaining
units and 21.0% of controls (about 2:1; at 3% prevalence that is about 6% precision). A detector has to tell a copy
from an allele, not find novelty. 58 gaining units have no haplotype-only sequence: no sequence rule can reach them.

## Detector (`bench/exgain/exseq.py`; reads only the BAM, `copies.tsv`, `net.tsv` and the reference)
Per unit: for each net read (primary record) over a unit locus, the read's exon sequence = query bases between the
aligned positions of the locus's first and last base (both must be aligned; else the read is skipped), on the locus
strand. Reads with identical exon sequence form one sequence. SUPPORTED = at least 3 reads and at least 10% of the
unit's exon-sequence reads. Two arms, thresholds fixed here, no tuning, no dev set:
- **Arm H (haplotype count):** fires when SUPPORTED sequences > 2 x (number of reference loci of the unit).
- **Arm N (divergent novelty):** fires when a SUPPORTED sequence is at edit distance >= 2 from every reference locus
  sequence of the unit.

## Scoring (`bench/exgain/exseq_score.py`, the only script that reads truth), once, on every unit with net reads
Per arm: fires in gaining / non-gaining units, precision, recall of the 469, share of the 213 VISIBLE gaining units
(row 1694's level B) fired. Baseline: `copy_phasing` putative tier, 7 gaining / 22 non-gaining fires.
- **PASS (per arm):** precision >= 0.50 over all units AND fires on >= 25% of the VISIBLE gaining units.
- Right reason, reported: of an arm's fires in gaining units, the share where a non-reference SUPPORTED sequence equals
  a haplotype-only copy sequence. Below 0.5, a PASS is recorded as PASS (WRONG REASON) and nothing is proposed.
- Stability, reported: both numbers on the two halves of the families (parity of sha1 of the family id).
- **Stop rule:** if neither arm passes, the item is recorded as a stated limit of the thesis; no third arm unless the
  user asks. A PASS proposes an opt-in report only; wiring is the user's decision.

## Outcome (2026-11-14, job 4471248, 49 s): FAIL on both arms; exon gains recorded as a stated limit
Scored once on the 14,815 units with net reads (457 gaining; 14,258 with at least one exon-sequence read).
- **Arm H:** 0 gaining / 18 non-gaining fires. Precision 0.000, VISIBLE 0 of 213. FAIL.
- **Arm N:** 47 gaining / 380 non-gaining fires. Precision 0.110 (bar 0.50), VISIBLE 44 of 213 = 0.207 (bar 0.25).
  FAIL on both bars. Right reason 41 of 47. Halves: 15 / 186 and 32 / 194.
- Baseline (`copy_phasing` putative tier): 7 / 22.
Reading: a divergent read-supported exon sequence is enriched in gaining units (10% of them against 2.6% of the
others, about 4:1) and is usually the gained copy's own sequence, but at 3% prevalence that is one right call in nine;
more than two haplotypes at one exon never appears where a gain is. By the stop rule no third arm is run: a detector
for exon gains is a stated limit of the thesis. Nothing is wired.
Post hoc look at arm N's 380 non-gaining fires (descriptive, after the verdict): in 245 the divergent supported
sequence equals a copy sequence on the maternal or paternal haplotype with no count difference, i.e. an allele of an
existing copy that differs from the primary assembly by >= 2 edits (the allele-versus-copy confusion the arm was meant
to avoid); in 135 it equals no copy sequence of the exon on any of the three assemblies (not resolved: reads from
loci below the exon mapping's 0.8 coverage / 0.90 identity, or errors). 26 of the 380 units have fewer copies on one
haplotype than the reference. Median edit distance 3; 94 at >= 5.
