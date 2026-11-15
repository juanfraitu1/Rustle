# PREREG row 1716 (2026-11-14): are the nested chains real? SQANTI3 with CAGE peaks, A119b

User 2026-11-14: "what if we try those chains and use SQANTI3 to see if they are real". Reading taken of "those
chains": row 1715's nested candidates on A119b (read chains with >= 5 reads, terminal sub-run of an emitted chain),
in four groups: emitted unannotated, not emitted unannotated, and annotated kept / annotated lost as positive
controls. Human only (the CAGE data is human). Written before any CAGE count. Seen: rows 1714 / 1715.

## What SQANTI3 can and cannot say here
Against RefSeq its structural category is fixed by construction (an unannotated nested chain is an incomplete splice
match, an annotated one a full match): the category is not a verdict and is not read. The evidence is CAGE support of
the 5' end (`within_CAGE_peak`, `dist_to_CAGE_peak`) for nested starts. For nested early ends SQANTI3's polyA motif and
downstream A content repeat row 1715's test with its own definitions; a polyA-site atlas could not be fetched (the
source's certificate has expired; not bypassed), so early ends get no new evidence and are reported only.

## Inputs
- SQANTI3 6.1 (`sqanti3_qc.py`, RefSeq CHM13 GTF, `chm13v2.0.fa`, `--CAGE_peak`, `--polyA_motif_list`); the peak file
  and motif list are the ones SQANTI3 ships (refTSS v3.1 hg38, taken from the v5.2.2 tag), peaks lifted to CHM13 with
  UCSC liftOver `hg38ToHs1` (223,146 of 224,694 lift, 99.3%).
- One model per candidate, same construction in every group (`nested_sub.py --models`): the candidate's intron chain,
  variable end = the commonest end of its own reads, other end = the container's. Containers that are annotated get a
  model too, with the commonest start of their own reads (group FULL).

## Gates (the row is UNDERPOWERED, not negative, when one fails)
- G1: CAGE support (`within_CAGE_peak`) of FULL models >= 0.40.
- G2: annotated nested starts (kept + lost) show CAGE support >= 3 x the background.
Background = the same test (5' end inside a lifted peak on the model's strand; and within 50 bp) at each candidate's
start shifted 100 bp downstream in the transcript, computed by `sq_nested.py`; the script's own flag on the unshifted
start must agree with SQANTI3's on >= 0.99 of models or the background is not used.

## Reported and read
Per group, all and `peaked` (>= 0.5 of reads within 10 bp of the commonest start) separately: models, share within a
peak (SQANTI3), share within 50 bp, background shares. Reading (fixed here), on the unannotated groups:
- "mostly real starts" when the share within a peak is >= half the annotated nested share;
- "mostly not" when it is <= 2 x background;
- between: mixed, and the count with CAGE support is the number RefSeq understates by.
Limits stated up front: A119b is testis and the peaks are a multi-tissue aggregate, so a peak is evidence and its
absence is weak; refTSS includes annotation-derived promoter sets, so it is not independent of annotation at known
genes. CAGE is external human data: no reference-free default can follow from this row. It can only correct how the
RefSeq-based rows 1714 / 1715 are read.

## Outcome (2026-11-14, job 4471366; `r1716_sqanti/score.txt`, SQANTI3 run 1,197 s on 148,724 models)
Gates: G1 PASS (FULL models within a peak 0.807); G2 PASS (annotated nested starts 0.557 = 1,165 of 2,093, background
0.066, 8.4 x). The script's flag equals SQANTI3's `within_CAGE_peak` on 148,724 of 148,724 models.

Nested starts (NEST5), share of 5' ends within a CAGE peak / shifted background:

| group | models | within peak | background | peaked stratum |
|---|---|---|---|---|
| annotated, emitted | 1,666 | 0.632 | 0.072 | 0.803 (441) |
| annotated, not emitted | 427 | 0.262 | 0.042 | 0.279 (86) |
| unannotated, emitted | 19,625 | 0.066 | 0.016 | 0.147 (2,917) |
| unannotated, not emitted | 77,756 | 0.028 | 0.017 | 0.030 (16,754) |

Reading by the fixed rule (half the annotated nested share = 0.278):
- **Unannotated, not emitted: mostly not real starts** (0.028 = 1.6 x background, under the 2 x line). The default
  ISM rule is keeping out chains that CAGE does not support either.
- **Unannotated, emitted: mixed** (0.066 = 4.1 x background, far below 0.278): about 1,300 of 19,625 have CAGE
  support (429 of the 2,917 peaked ones); the other 93% have none. RefSeq understates this class by roughly that
  1,300; most of what the default emits here is not a CAGE-supported start.
- The annotated chains the default loses are themselves weakly supported (0.262 against 0.632 for the kept ones).
Early ends (NEST3, SQANTI3's own columns, reported only): polyA motif 0.451 annotated-emitted, 0.376 annotated-lost,
0.272 unannotated-emitted, 0.230 unannotated-not-emitted; A-rich downstream (>= 60%, internal-priming signature)
0.173 / 0.261 / 0.563 / 0.620. Unannotated early ends are mostly internal priming by SQANTI3's test.
Limits as registered (testis sample vs aggregate peaks: absence of a peak is weak; human only). Disclosed: data files
from the SQANTI3 v5.2.2 tag, tool 6.1; no polyA-site atlas.
