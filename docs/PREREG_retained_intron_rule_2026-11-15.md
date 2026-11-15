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

## OUTCOME (2026-11-15, job 4472243, `r1745_retained/`; binaries of 2026-11-14 22:01)

R10 reproduces the base runs (OR6737 25,794 correct at 0.352; A119b 38,715 at 0.176).

| sample | arm | RI both emitted | CASSETTE | multi-exon chains | correct | precision | new chains | new correct |
|---|---|---|---|---|---|---|---|---|
| OR6737 | R10 | 0.744 (32 of 43) | 0.975 | 73,348 | 25,794 | 0.3517 | | |
| OR6737 | R30 | 0.953 | 0.977 | 74,806 | 25,865 | 0.3458 | 1,458 | 71 (0.049) |
| OR6737 | R0 | 1.000 | 0.977 | 76,381 | 25,881 | 0.3388 | 3,033 | 87 (0.029) |
| A119b | R10 | 0.682 (131 of 192) | 0.912 | 219,517 | 38,715 | 0.1764 | | |
| A119b | R30 | 0.812 | 0.917 | 226,156 | 38,951 | 0.1722 | 6,641 | 236 (0.036) |
| A119b | R0 | 0.932 | 0.919 | 235,432 | 39,070 | 0.1660 | 15,916 | 355 (0.022) |
| KB3781 | R10 | 0.584 (52 of 89) | 0.942 | 72,102 | 25,844 | 0.3584 | | |
| KB3781 | R30 | 0.764 | 0.943 | 74,444 | 25,911 | 0.3481 | 2,342 | 67 (0.029) |
| KB3781 | R0 | 0.933 | 0.945 | 78,130 | 25,951 | 0.3322 | 6,028 | 107 (0.018) |

1. **Attribution: the filter is the cause.** With it off the RI share is at or above the cassette control on both
   deciding samples (1.000 vs 0.977; 0.932 vs 0.919) and on KB3781 (0.933 vs 0.945). No arm loses a chain (0-2 lost).
2. **Cost: neither arm is eligible.** R0 costs 0.013 and 0.010 precision. R30 closes more than half of the gap on both
   samples and adds correct chains, but costs 0.0059 on OR6737 (bar 0.005; A119b 0.0042). FAIL as registered, by 0.0009
   on one sample. The filter stays at 10; the retained-intron gap is a stated limit of the default and
   `RUSTLE_POLISH_RETAINED_RATIO` is the opt-in.

The chains the filter removes are 2-5% RefSeq-correct, so it is doing its job on the bulk; the annotated ones it also
removes are 11 / 48 / 31 pairs here. Not known: how many of the unannotated retained-intron chains are real isoforms
and how many precursor RNA (RefSeq cannot say).

**User ruling (2026-11-15):** the filter stays at 10. `RUSTLE_POLISH_RETAINED_RATIO` remains the opt-in.
