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
