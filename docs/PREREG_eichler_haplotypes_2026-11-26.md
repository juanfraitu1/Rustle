# PREREG: the Eichler lab's 10-point filter run as the paper describes it, on an individual's own haplotypes (row 1972, 2026-11-26)

The user: run that situation, to show the advisor that the lab's rule was followed. Until now only our ADAPTATION to a
single reference was run (rows 1775-1786, 1954). The paper (Guitart et al. 2024, Methods, "Iso-Seq and transcript
analyses"): FLNC reads were mapped to all validated human haplotypes with `minimap2 -ax splice --sam-hit-only
--secondary=yes -p 0.5 --eqx -K 2G -G 8k -N 20`; "we compared these primary alignments to one another and considered the
cluster paralog from which they were derived. Any Iso-Seq read with primary minimap2 alignment scores of 10 or greater
for a given paralog cluster relative to all other cluster mappings were retained, while other mappings were marked as
ambiguous and ignored."

## What is run (`bench/o2_truth/eichler_haplotypes.py`)
HG002: its own maternal and paternal assemblies (hg002v1.1, workspace `r1776_hg002_room`), the 136,074 family reads of
row 1776's selection (a 25% hash sample of the reads on the truth families).
1. The reads mapped to each haplotype separately with the paper's minimap2 parameters. Per read, its PRIMARY alignment on
   each haplotype and that alignment's AS.
2. The paralog a haplotype locus belongs to: the CHM13 family genes' genomic sequences aligned to each haplotype
   (`minimap2 -x asm5`); a haplotype locus is labelled with the reference copy that covers it at the highest identity.
3. The rule: both primaries on the same paralog (or only one haplotype has a mapping) = RETAINED for it; on different
   paralogs = retained for the higher one when its AS leads by 10 or more, otherwise AMBIGUOUS.
4. Beside it, the same reads on CHM13 as our pipeline maps them (`-N 50 -p 0.1`): the paralog under the primary; near tie
   (another alignment within 10 `ms`); our single-reference adaptation (primary AS >= 10 above every other alignment).

## Where this still differs from the paper (said in every statement of the result)
- Two haplotypes of ONE individual, the one the reads come from; the paper maps to all validated human haplotypes and
  its reads come from other individuals.
- The paper's paralog CLUSTERS come from its own phylogeny across haplotypes; here a cluster is the nearest reference
  copy. Finer than theirs, so more reads fall between clusters here than would in the paper.
- The paper keeps reads of >= 1,000 bp and >= 99.9% quality; these reads carry no quality value, the length rule is applied.
- All truth families, with TBC1D3 (the paper's family) and NPIP reported on their own.

## Claims
- **E1:** of the reads the procedure retains for a paralog, our CHM13 primary is on that same paralog for at least 95%.
- **E2:** of the reads it marks ambiguous, at least 70% are near ties on `ms` on the single reference.
- **E3:** the single-reference adaptation abstains on at least 90% of the reads the paper's procedure marks ambiguous.
## Predictions (mine)
E1 holds (97% or more). E2 holds. E3 holds. So: followed as written, the rule agrees with the primary where it speaks
and is silent where the reads are tied; the adaptation was a fair stand-in.

## Result (job 4510253; 1,459 family genes; maternal 133,021 and paternal 131,845 reads with a primary)
One change AFTER the first scoring, said plainly: the first paralog label took the reference copy with the highest
identity over its whole alignment, and a short local alignment of a neighbouring paralog could beat the full one. The
label now counts the matching bases of each reference copy INSIDE the read's span on the haplotype. First version:
retained 97.9%, same paralog 90.0%, TBC1D3 retained 97 of 97. Figures below are the second version
(`r1972_eichler_haplotypes/score_window.txt`); the verdict on the three claims is the same under both.

| | all truth families | TBC1D3 | NPIP |
|---|---|---|---|
| reads with a haplotype primary on a labelled family locus | 116,099 | 97 | 1,463 |
| RETAINED by the procedure | 114,850 (98.9%) | 52 (53.6%) | 1,346 (92.0%) |
| AMBIGUOUS | 1,249 (1.1%) | 45 (46.4%) | 117 (8.0%) |
| retained: our CHM13 primary on the same paralog | 91.5% | 26.9% | 50.9% |
| retained: on another paralog of the family | 8.3% | 73.1% | 44.7% |
| retained: on the single reference clear / near tie / exact tie (`ms`) | 94.7% / 2.4% / 2.9% | 3.8% / 53.8% / 42.3% | 55.1% / 33.5% / 11.4% |
| retained: the single-reference adaptation abstains | 6.2% | 94.2% | 40.4% |
| ambiguous: on the single reference clear / near tie / exact tie | 42.0% / 37.3% / 20.7% | 2.2% / 28.9% / 68.9% | 23.1% / 56.4% / 20.5% |
| ambiguous: the single-reference adaptation abstains | 50.2% | 88.9% | 67.5% |
Retained reads by their state on the single reference (all families): clear 108,729, same paralog 93.5%; near tie
2,783, 62.0%; exact tie 3,338, 51.8%.

### Claims
- **E1 FAILS:** same paralog 91.5% (bar 95%).
- **E2 FAILS:** of the ambiguous reads 58.0% are near or exact ties on the single reference (bar 70%).
- **E3 FAILS:** the single-reference adaptation abstains on 50.2% of the ambiguous reads (bar 90%).
All three of my predictions were wrong.

### Reading
- **Run as the paper describes it, the rule decides almost every read** (98.9% retained): on the individual's own
  haplotypes a read has a home that fits it better than any other paralog.
- **On a single reference the same reads are often ties.** TBC1D3, the paper's family: 96% of its reads are near or
  exact ties on CHM13 and our adaptation abstains on 92%; the paper's procedure retains 54% of them. Of those it retains,
  our CHM13 primary is on the same paralog for 27%: among ten or so near-identical reference copies the primary is a draw.
- **So the single-reference adaptation was not a fair stand-in for the paper's rule** in the families the rule was made
  for. What it tests is a different thing (does the primary lead on one reference), and its abstentions say "tied on
  this reference", not "ambiguous in this genome". The rule's power comes from the haplotype assemblies, not from the
  number 10.
- **Where the reads are clear on the single reference** (94.7% of the retained), the primary is on the paralog the
  haplotypes name for 93.5%. The other 6.5% are not all errors of the primary: a label is "the nearest reference copy",
  and in arrays whose copies differ between people the individual's copy can be nearest to another reference copy (row
  1971 found the same thing from the other side). This run cannot split that 6.5% into the two.
- Limits as registered: one individual's two haplotypes; clusters = nearest reference copy, not the paper's phylogeny
  (so MORE reads fall between clusters here than would there); no quality filter; TBC1D3 has only 97 reads in this sample.
