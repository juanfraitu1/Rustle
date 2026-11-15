# PREREG r1645: does FLAIR's realignment place tied multimapping reads better than the primary? (2026-11-12)

Question (user): FLAIR realigns every read to its first-pass isoforms; does that second look rescue tied reads?
Sample and truth: gorilla KB3781 (the assembled individual), r1509's set: the `o2mstie/ms10` arm's contested rows, truth =
`bench/o2_contested/strict_truth.py` at g = 0 (the read's best hit on its own haplotypes is unique and names one catalog
copy locus). On that set the primary is right on 280 of 286 assigned rows and O2 on 284.
Arm: the FLAIR emulation E3 of r1640 (`bench/famcmp/flair_like.py`): first pass from the whole KB3781 BAM, then only the
truth reads are realigned to all first-pass sequences with FLAIR's command and assigned by FLAIR's rule. A read's FLAIR
placement = the genomic position of the isoform it is assigned to, judged exactly as the primary's position is.
Readout, on decidable rows (assigned and abstained) where the read has a FLAIR placement: FLAIR right / wrong, the primary
right / wrong on the same rows, and how many reads FLAIR moved off the primary's copy.
- "FLAIR places tied reads better": FLAIR wrong <= half the primary's wrong AND at least 5 fewer wrong.
- "Worse": the mirror. Anything else: no measurable difference on this set.
Limits fixed now: emulation, not FLAIR's own binary (gorilla Jaccard 0.732 with FLAIR's output); only reads the
haplotypes can decide are scored; one fibroblast sample.

## Result (2026-11-12, `r1645_flairtied/`, job 4466720): WORSE than the primary on this set
Truth dir `o2mstie/truth2_ms10` (`truth_ms10` has no hits table). strict_truth today gives 315 decidable assigned rows
(not r1509's 286; its copy-locus rule was tightened after r1509): O2 299 right / 16 wrong, primary 296 / 19.
First pass 483,687 chains; 3,641 of 3,699 truth reads assigned.
| rows (g = 0, decidable) | FLAIR-like right / wrong | primary right / wrong, same rows |
|---|---|---|
| all, as registered (assigned + abstained) | 315 / 155 | 446 / 31 |
| spliced reads only (352 rows with a placement) | 298 / 53 | 337 / 14 |
| unspliced reads (119 rows; not a fair test, see limit) | 17 / 102 | 102 / 17 |
Spliced reads: the realignment moved 39 reads off the right copy (25 more than 100 kb away, 12 at 5-100 kb, 2 to another
chromosome; none is a position artefact) and corrected none of the primary's 14 wrong ones. Verdict by the registered
rule: worse (primary wrong 14 <= half of 53, 39 fewer). Unspliced reads are confounded: the emulation's first pass has no
single-exon isoforms, so they can only land on a spliced isoform; on them it did correct 8 of the primary's 17 wrong.
Limits as registered: emulation, not FLAIR's binary; one sample; only haplotype-decidable reads.
