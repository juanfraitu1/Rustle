# PREREG: IsoCon's base-quality model in the sequence test, on the two samples with haplotype truth (row 1933, 2026-11-23)

The port of IsoCon's test (row 1926) uses its EMPIRICAL error model (a read's error rate from its own alignment). IsoCon's
main mode uses the base qualities of the CCS reads. Asked what of IsoCon is not ported, this was the first candidate;
user: try the quality-value model. Checked first: the aligned BAMs of KB3781 and HG002 carry base qualities in primary
and secondary records (KB3781 up to 93; HG002 Kinnex binned, top bin 50, mostly 40).

## What is ported (`recurrence_test.py --qual`, opt-in; from `functions.get_read_ccs_probabilities_c / _t`, `ccs_info.py`)
Per read and per separating variant the error probability is 10^(-Q/10) of the read base at that position, times the
share of that error type among the reads' errors (/ 3 substitution, / 4 insertion); in a homopolymer the full
10^(-Q/10); qualities inside a homopolymer run of the read are sorted first; with qualities IsoCon applies NO placement
factor. Two fixed forms, no tuning:
- `isocon`: Q rescaled as IsoCon does, (Q - 3) x 40 / 90 + 3 (written for qualities that run to 93).
- `phred`: Q used as it is, capped at 43 (for binned qualities that top out near 40-50).
Planted reads have one constant quality, so this cannot be developed or tested on them: real samples only.

## Test
KB3781 and HG002, the runs of rows 1928 and 1931 (same assembly, BAM, groups); the sequence test rerun with each form;
the reads of the kept sequences aligned to both haplotypes (`-N 200 -p 0.2`) and the site-by-site lookup of rows 1930 /
1931 (`explain.py`) unchanged. Baseline (empirical model), known: KB3781 107 non-reference sequences, 52 in the DNA
(47 best hit + 5 one lower place), 47 nowhere (0.44); HG002 689, 517 in the DNA, 98 nowhere (0.14).

## Claim (per form)
BETTER when, on BOTH samples, the share of non-reference sequences found nowhere falls by >= 0.05 and the form keeps
>= 0.90 of the baseline's sequences confirmed in the DNA (KB3781 >= 47, HG002 >= 466). Otherwise not adopted.

## Predictions (mine)
Neither form is better. `phred` gives smaller error probabilities than the empirical model and drops the placement
factor, so it should keep MORE sequences and not a lower unexplained share. `isocon` rescales Q 40 to about Q 19, which
is far above the real error rate of these reads (one per kb), so it should keep fewer sequences and lose more than a
tenth of the confirmed ones. If the unexplained class on KB3781 were low-quality bases, `isocon` would remove it: that
is the one outcome that would teach something new.

## What happened on the way
The first job (4498880) finished KB3781 in 7 minutes and was cancelled at its 8-hour limit inside HG002's sequence
test. Cause, found with the cache and a profile: ONE isoform at HG002's immunoglobulin lambda locus (chr22 23.18-23.33
Mb; 26,776 tied reads, 4,922 distinct variant patterns; HG002's RNA is from a B-cell line) took 43 of the tool's 47
minutes. A hypervariable locus is not a multi-copy family; the user agreed to a guard: `--max-patterns` (default 1000),
such an isoform is not tested and is listed in `OUT.skipped.tsv`. HG002 from the cache then takes 29 s. The guard
changes nothing on the planted sets, KB3781 or OR6737 (tables identical). HG002's arms below are with the guard, the
lookup on the haplotype alignments of row 1931.

## Result (non-reference sequences; in the DNA = at the read's best hit or one lower-scoring place)

| | KB3781 empirical | isocon | phred | HG002 empirical | isocon | phred |
|---|---|---|---|---|---|---|
| sequences scored | 107 | 91 | 73 | 648 | 530 | 489 |
| in the DNA | 52 | 50 | 46 | 476 | 419 | 393 |
| partly | 8 | 2 | 2 | 74 | 33 | 27 |
| nowhere | 47 (0.44) | 39 (0.43) | 25 (0.34) | 98 (0.151) | 78 (0.147) | 69 (0.141) |

- **Neither form is BETTER.** `isocon`: the share found nowhere moves by 0.01 or less on both samples. `phred`: it falls
  by 0.10 on KB3781 but keeps 46 of the 47 confirmed sequences required; on HG002 it falls by 0.01 and keeps 393 of the
  428 required (0.90 of 476). Not adopted; both stay as options, off.
- What the qualities do: they remove sequences (HG002 648 -> 530 / 489), confirmed and unconfirmed alike. On KB3781 the
  plain form removes proportionally more unconfirmed ones, on HG002 it does not. Base quality is not what separates the
  sequences the DNA confirms from those it does not.
- Predictions: "neither is better" right; the direction for `phred` (I expected more sequences kept) wrong on both samples.
- **Row 1931's HG002 figures with the guard:** 648 non-reference sequences, 476 in the DNA (0.735), 74 partly, 98 nowhere
  (0.151). The 41 sequences the guard removes were all confirmed, so the unexplained share moves from 0.142 to 0.151:
  against row 1931's fixed reading (<= 0.15 specific to KB3781, >= 0.30 general) it now sits ON the line, i.e. "low, and
  far below the gorilla's 0.44", not a clean "<= 0.15". Count check with the guard: purity 789 of 861, count 574 of 677.
