# PREREG row 1745 (2026-11-15): is the retained-intron filter what loses the annotated retained-intron isoforms, and what does relaxing it cost?

Row 1744: annotated pairs that differ by one retained internal intron, both expressed (minor >= 5 reads), are both
emitted 0.744 (OR6737) / 0.682 (A119b) / 0.584 (KB3781) of the time against a cassette control of 0.975 / 0.912 / 0.942;
FLAIR 0.953 / 0.906. Seen before this registration (from row 1744's pair tables): the lost form is the retained one when
it is the minor (lost in 10 of 19, 49 of 96, 33 of 58 such pairs); the spliced form as minor is almost never lost.
Written and committed before any arm was assembled.

## The rule under test

`copy_assign --polish-retained-ratio R` (default 10, `src/copy_assign_bin.rs` "§6za retained-intron filter"): a transcript
is dropped when one of its exons contains an intron of its gene whose summed read support is >= R x the transcript's
reads. New driver switch `RUSTLE_POLISH_RETAINED_RATIO` (unset = 10, the same command).

## Arms

`assemble` through the driver, today's release binaries, `RUSTLE_ASSEMBLY_LAYERS=off` (only `PREFIX.gtf` is read):
R10 (default, re-run so every arm shares one binary), R30, R0 (filter off). Samples OR6737 and A119b decide; KB3781
reported. Workspace `r1745_retained`.

## Measures

- `alt_events.py` (row 1744 classes) on each arm: RI and CASSETTE both-emitted share at minor >= 5 reads.
- `chain_prec.py`: multi-exon chains, RefSeq-correct chains, precision; new and lost chains of R30 / R0 against R10 with
  their correct counts.

## Reading rules

1. **Attribution (R0):** RI share within 0.05 of the same arm's CASSETTE share on both deciding samples = the filter is
   the cause. Less than half of the R10 gap recovered on either sample = another rule drops them; report what the
   remaining lost pairs have in common and stop.
2. **Cost (R30 and R0, each against R10, both deciding samples):** an arm is eligible as a default change (user
   decision, no flip in this row) when it closes at least half of the RI gap, chain precision falls by no more than
   0.005, and correct chains are not fewer. Otherwise the filter stays at 10 and the gap is a stated limit of the
   default, with the switch as the opt-in.

## Limits

RefSeq scores an unannotated retained-intron chain as false although it can be real or precursor RNA; the precision
cost is therefore an upper bound on the harm and says nothing on which. RI pairs were seen on all three samples in row
1744, so no sample is held out for the RI measure; the precision cost of R30 / R0 is unseen.
