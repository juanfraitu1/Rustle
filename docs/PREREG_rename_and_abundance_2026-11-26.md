# PREREG: two shortcomings of primary-only assemblers at tied reads (rows 1969 and 1970, 2026-11-26)

The user: are there other shortcomings of other assemblers that multimapping reads let us overcome; try both of the two
proposed. Both follow from one fact (row 1954 and the worked read of 2026-11-25): at an exact tie minimap2 orders the
alignments by a hash of the read's name, length, a seed and the first anchor, so the primary is a reproducible draw.

## Data (both rows)
New planted sets, seeds 1971-1973 (CHM13 families, 60 reads per isoform, ragged ends, as rows 1917 / 1952): scenarios
E (every copy expressed), S (half silent), D (one copy 10 x, the others a third). Tools as in `run_tools.sh`: Rustle on
primaries, Rustle mode with the labels, StringTie -L, FLAIR, IsoQuant, bambu (NDR 1).

## Row 1969: does the assembly change when the reads are renamed or the seed is changed?
Arms per set (E and S): `orig`; `again` (the same BAM, every tool run a second time: the tool's own run-to-run noise);
`renamed` (every read name prefixed with x, mapped again); `seed7` (`minimap2 --seed 7`).
`bench/o2_truth/stability.py`: per tool, against `orig`: spliced intron chains only in one of the two, planted copies
with a transcript only in one of the two (>= 50 exon bases), and in S the silent copies called in each. For Rustle's mode
also the copies labelled `expressed`.
- **R1:** for StringTie, FLAIR and bambu the set of copies with a transcript differs between `orig` and `renamed` by more
  copies than between `orig` and `again`, in at least 4 of 6 sets each.
- **R2:** the copies Rustle labels `expressed` differ between `orig` and `renamed` by fewer copies than the called copies of
  each of StringTie, FLAIR and bambu, in at least 5 of 6 sets.
- **R3:** in every S set, for each of StringTie, FLAIR and bambu, the silent copies called under `orig` and under
  `renamed` are not the same set.
IsoQuant is reported and not in the claims: it drops part of the mapping-quality-0 reads, so less is expected to move.
Predictions (mine): R1, R2, R3 hold. Rustle on primaries behaves like the others.

## Row 1970: per-copy abundance against the planted truth
`bench/o2_truth/abundance.py`: reads per planted copy from each tool's own quantification (Rustle: the EM of the labels,
`em` of `PREFIX.o2.em.tsv`; StringTie `longcov`; FLAIR `isoform.counts`; IsoQuant discovered transcript counts; bambu
`readCount`; and the plain count of primary alignments as the baseline every primary-only tool starts from), each
transcript given to the copy it shares the most exon bases with. Within every family of two or more copies: a copy's
SHARE of the family's reads against its planted share. Measures: mean absolute error of the share; in S the share of a
family's reads put on silent copies; both also by divergence class of the copy's unit (identical to another copy, 1-5
edits from the nearest, more).
- **Q1:** in at least 2 of 3 D sets Rustle's mean absolute share error is below that of every other tool.
- **Q2:** in at least 2 of 3 S sets Rustle puts a smaller share of the reads on silent copies than every other tool.
- **Q3:** for copies identical to another copy no tool is better than the plain primary count by more than 0.02 in D
  (nothing in the reads can split identical copies).
Predictions (mine): Q2 holds clearly; Q1 holds by a small margin; Q3 holds.
Limits fixed beforehand: sample = reference (no private variants), no read is a chimera, expression per copy is planted,
not learned from a real split; a result here is not a result on real tissue.

## Result, row 1969 (jobs 4510224-4510226)
Primaries that move: renamed 2,286 / 2,275 / 2,311 of 7,980 (E) and 1,271 / 1,160 / 1,262 of 4,200 (S); `--seed 7` the
same size. Run twice on the same BAM, every tool gives the same assembly (0 chains, 0 copies differ, all seven arms).
Original against renamed reads, seeds 1971 / 1972 / 1973 (seed 7 gives the same picture):

| tool | E: intron chains in one run only (of about) | E: copies in one run only | S: intron chains in one run only | S: copies in one run only | S: silent copies in one run only |
|---|---|---|---|---|---|
| Rustle, mode | 8 / 3 / 4 (365) | 0 / 0 / 0 | 0 / 3 / 0 (280) | 0 / 0 / 0 | 0 / 0 / 0 |
| Rustle, mode, kept by the EM | 8 / 3 / 8 (305) | 0 / 0 / 2 | 0 / 3 / 0 (170-200) | 0 / 0 / 0 | 0 / 0 / 0 |
| Rustle on primaries | 18 / 22 / 22 (265) | 0 / 0 / 0 | 23 / 28 / 24 (180) | 3 / 2 / 4 | 3 / 2 / 4 |
| IsoQuant | 8 / 13 / 11 (185) | 0 / 0 / 0 | 7 / 16 / 12 (125) | 2 / 2 / 3 | 2 / 2 / 3 |
| StringTie | 16 / 13 / 15 (220) | 0 / 0 / 0 | 19 / 34 / 32 (140) | 2 / 3 / 3 | 2 / 3 / 3 |
| FLAIR | 32 / 36 / 45 (265) | 1 / 1 / 0 | 27 / 42 / 42 (175) | 4 / 5 / 6 | 2 / 5 / 5 |
| bambu | 80 / 69 / 63 (305) | 0 / 0 / 0 | 54 / 61 / 52 (190) | 3 / 2 / 5 | 3 / 2 / 5 |
- **R1 FAILS as registered:** the set of COPIES changes in the half-silent sets only; where every copy is expressed
  every tool finds the same copies either way (StringTie 3 of 6 sets, bambu 3 of 6, FLAIR 5 of 6; bar 4 each).
- **R2 FAILS as registered** for the same reason: in the all-expressed sets nobody's copies change (and ours by 2 in one).
- **R3 HOLDS:** in all three half-silent sets StringTie, FLAIR and bambu call different silent copies after renaming.
- What the table shows beside the claims: the TRANSCRIPTS change in every primary-only tool in every set (6-25% of the
  intron chains; bambu most, IsoQuant least), and in the half-silent sets 2-6 copies appear or vanish, every one of them
  a silent copy. Rustle's mode changes 0-8 chains and no copy in any of the six sets; Rustle on primaries behaves like
  the other primary-only tools. My predictions for R1 and R2 were wrong in their unit (copies in E sets), right in kind.

## Result, row 1970. Mean absolute error of a copy's share of its family's reads (D = one copy 10 x; S = half silent)
| tool | D: all copies, seeds 1971 / 1972 / 1973 | D: copies identical to another | D: 1-5 edits | S: share of reads on silent copies |
|---|---|---|---|---|
| primary alignments, counted | 0.039 / 0.071 / 0.052 | 0.065 / 0.139 / 0.120 | 0.019 / 0.042 / 0.011 | 0.155 / 0.185 / 0.161 |
| StringTie | 0.040 / 0.074 / 0.054 | 0.064 / 0.143 / 0.117 | 0.020 / 0.041 / 0.015 | 0.156 / 0.178 / 0.153 |
| Rustle, EM of the labels | 0.050 / 0.071 / 0.064 | 0.095 / 0.141 / 0.141 | 0.004 / 0.036 / 0.004 | 0.160 / 0.190 / 0.135 |
| FLAIR | 0.053 / 0.081 / 0.069 | 0.089 / 0.151 / 0.150 | 0.029 / 0.051 / 0.031 | 0.188 / 0.246 / 0.172 |
| IsoQuant | 0.057 / 0.062 / 0.077 | 0.092 / 0.145 / 0.128 | 0.012 / 0.006 / 0.008 | 0.174 / 0.177 / 0.188 |
| bambu | 0.069 / 0.074 / 0.068 | 0.088 / 0.154 / 0.142 | 0.024 / 0.020 / 0.009 | 0.160 / 0.189 / 0.188 |
- **Q1 FAILS** (0 of 3: StringTie, IsoQuant or the plain count is lower in each). **Q2 FAILS** (lowest in 1 of 3).
  **Q3 HOLDS** (no tool is more than 0.02 better than the plain count on identical copies).
- Both of my predictions for a gain were wrong. A plain count of primaries is as accurate as anything, ours included:
  the tie-break is an unbiased draw, so over 60 reads per isoform the shares average out. Our EM is the best only for
  copies 1-5 edits from their nearest relative (0.004 / 0.036 / 0.004) and is worse than the plain count on identical
  copies and overall. About 16-19% of a family's reads sit on silent copies in every tool, ours too: the labels mark
  those copies, the abundance column does not move the reads away.

## Reading of both rows
- Supported: a primary-only assembly is not a function of the data alone. Renaming the reads changes 6-25% of its
  transcripts in family regions and, where some copies are silent, WHICH silent copies it reports. Rustle's mode does not
  have this (0-8 chains, 0 copies), because it seeds every tied place instead of the drawn one.
- Not supported: that primary-only tools quantify copies worse. At this depth they do not; the draw is fair. Not tested:
  low depth, where a fair draw has large variance (3 reads per isoform), and real tissue.
