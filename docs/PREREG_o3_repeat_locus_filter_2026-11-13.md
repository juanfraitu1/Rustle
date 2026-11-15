# PREREG row 1698 (2026-11-13): a repeat-locus filter, step 0

User go 2026-11-13 ("O3: repeat-locus filter"). The O3 ledger lists "the repeat-locus filter" as open without a
definition; the definition below is this session's, not an earlier registered one.

## Definition (frozen from rows 1690 and 1691)
A locus is REPEAT when the exons of its longest transcript are >= 0.9 soft-masked, or >= 50% covered by minimap2 asm20
alignments to `tools/data/rdna_human.fa`. No new parameter.

## Labelled set
HG002 default assembly on the paternal assembly (`r1669_monofeed_hg/base/hg.gtf`), labels of row 1672 / 1673
(`r1672_shadowflag/hg/loci.tsv`): the 1,518 secondary-only loci with haplotype truth, among them the 106 "real" loci
that broke the shadow flag, 103 of which sit on unplaced scaffolds and are rDNA-array sequence (row 1673 follow-up).
Row 1673 already showed the soft mask alone removes none of them; row 1691 reached 92 of 92 CHM13-placed ones with the
rDNA unit, from another sequence set. Neither measured what the combined rule removes among ordinary loci.

## Bars (written before the count)
- Reach: REPEAT covers at least 80% of the 103 scaffold loci labelled L.
- Cost: REPEAT covers at most 2% of the chromosome loci that are not flagged (about 17,700), and at most 2% of the
  chromosome loci labelled P (shadows of real genes: removing those is harmless, but a filter that hits them is not
  specific to repeats).
All three hold = PASS: the rule is a usable repeat-locus filter for locus lists on this sample, to be tested on a second
sample before any wiring. Otherwise STOP. Reported without a bar: every where x label cell, and row 1672's L share
among judged flagged loci with REPEAT loci set aside.
Limits: one sample; the paternal assembly's soft mask; nothing is wired by this row.

## Outcome (2026-11-13, job 4470933): PASS as registered
| where | label | loci | soft mask >= 0.9 | rDNA >= 0.5 | REPEAT |
|---|---|---|---|---|---|
| scaffold | L (the "real" loci of row 1672) | 103 | 0 | 102 | 102 (0.990) |
| chromosome | not flagged | 17,676 | 268 | 41 | 297 (0.017) |
| chromosome | P (shadows of genes) | 135 | 0 | 1 | 1 (0.007) |
| chromosome | other | 120 | 78 | 5 | 83 (0.692) |
| chromosome | L / both / mixed / thin | 3 / 12 / 3 / 8 | | | 1 / 0 / 2 / 5 |
| scaffold | P / both / mixed / other / thin | 98 / 76 / 160 / 260 / 540 | | | 97 / 76 / 160 / 243 / 527 |
| scaffold | not flagged | 699 | 5 | 605 | 605 (0.866) |

- Reach 102 of 103 (bar 80%); cost 1.7% of unflagged and 0.7% of P chromosome loci (bars 2%). All three hold.
- The two halves do different work: the rDNA unit takes the unplaced scaffolds almost whole (1,810 of 1,936 loci,
  flagged or not), the soft mask takes the chromosome `other` class (78 of 120) and 268 unflagged chromosome loci.
- Row 1672 re-read with REPEAT loci set aside (descriptive): judged flagged loci 970 -> 206; L 106 -> 3 (0.015; row
  1672's bar <= 5 and <= 0.03 would hold); P 233 -> 135 (0.66; its bar 0.75 still fails). The shadow flag stays closed.
- Not done: a second sample (the bar for wiring), and what the 297 unflagged chromosome loci are (268 by the mask:
  expressed repeats; row 1690 found 2.6% of assigned rows in LINE families). Nothing is wired.
Table: workspace `r1698_repeatloci/hg.tsv`.

## Row 1702 (2026-11-13, registered before the run): the second sample
Row 1698's rule, unchanged (exons >= 0.9 soft-masked, or >= 50% covered by `tools/data/rdna_human.fa`), on gorilla
KB3781: default assembly `r1665_monofeed/k_base/k.gtf` on `GGO.fa`, labels of row 1671 (`r1671_seconly/kb0/loci.tsv`:
P 155, both 8, other 10, thin 21). This sample has no labelled real repeat loci (no L class, and the reference has no
unplaced scaffolds), so reach cannot be tested here; it is a test of the cost on another genome's mask and another
species with the human rDNA unit.
Bars: REPEAT covers at most 2% of the unflagged loci and at most 2% of the P loci. Both hold = PASS: the filter is
safe on a second sample. Reported without a bar: every label cell, and how many loci each half takes.
Wiring stays the user's decision either way.
Outcome: **PASS as registered.** KB3781: REPEAT covers 38 of 14,666 unflagged loci (0.003; 12 by the mask, 26 by the rDNA unit) and 1 of 155 P loci (0.006); bars 2%. Also 6 of 10 `other`, 0 of 8 `both`, 1 of 21 `thin`. The cost is lower than on HG002 (1.7%), where the paternal assembly's mask takes 268 expressed-repeat loci. Reach is still measured on one sample only (HG002, 102 of 103). Nothing is wired.

Wiring (user ruling 2026-11-13, "Wire as opt-in"): driver stage `repeatloci` writes `PREFIX.repeat_loci.tsv` (needs `PREFIX.gtf` and `--fasta`; rDNA unit `RUSTLE_REPEAT_FASTA`, default `tools/data/rdna_human.fa`). It marks loci and removes nothing; no default stage runs it. KB3781 through the driver: 46 of 14,860 loci marked.
