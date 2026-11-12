# PREREG EXON-SPACE PLACEMENT, UNSPLICED-PRIMARY CLASS ON HG002 (register 1644, 2026-11-12)

Why: r1643 (multi-isoform reference, KB3781) equalled the primary overall but split by the primary's shape: unspliced
primary 12 vs 3 for the arm, spliced primary 19 vs 33 against. That split was post hoc. This row registers it on a
second sample.

Sample: HG002 Kinnex on its paternal haplotype (`hg002.pat.bam`, base run `~/runs/hg002_kinnex/hg.*`, 4,066 copies).
Truth: HG002's own two haplotypes (pat `~/data/hg002/hg002.pat.fa`, mat GCA_021951015.1 downloaded for this row),
`haplotype_truth.py` + `strict_truth.py` at slack 0, as r1514. Same caveat as KB3781: reads, reference and truth are one
individual.

Set (fast-feedback panel, no genome-wide run): extended `as_table` of the BAM; `panel.py` top 40 families by ms-near-tie
reads (W 10) among the top 80, skipping families whose regions sum to > 2 Mb (cost guard, hotspot rule); `copy_assign --tie-ms 10` on the panel (primary-on-family guard on, binary `target_r17`); contested rows
with truth. CLASS = molecules whose primary record is unspliced (as_table column 11 = 0). Per molecule.

Arm: r1643 unchanged (`multi_iso_ref.py` on `hg.fam.copies.*` + `hg.gtf`; `minimap2 -cx map-hifi --eqx -N 100 -p 0.1`;
best AS per copy, margin >= 1, cover >= 0.90), `EXON_PLACE_SHAPE=0`.

Bars, on the class:
0. Power: >= 100 decidable class molecules, else UNDERPOWERED (no verdict).
1. Reach: arm assigns >= 30%.
2. Arm right / primary wrong >= 2 x (arm wrong / primary right), first count >= 10.
3. Bar 2's ratio with the largest family removed (first count >= 5).
4. Negative class molecules assigned <= 10% (r1642 0.456, r1643 0.587: expected to FAIL; O2's rate beside it).
5. Wall <= 3 h for the whole job (limit = FAIL).
Readout: the spliced-primary class and all molecules, same tables; vs O2.

Verdicts: 1-4 pass = PASS. 1-3 pass, 4 fails = ACCURACY ONLY: the class effect holds on a second sample, usable only
behind an outside-catalog guard, which would be its own row. 2 or 3 fails = FAIL, class closed with the rest.

Outputs: `/anvme/workspace/v138da13-rustle/r1644_hg002_unspliced/`, job script `job.sbatch` there.

## RESULT r1644 (2026-11-12, job 4466707, 56 min; `r1644_hg002_unspliced/result.shape0.txt`, `split.txt`)
Panel: 37 families (43 of the top 80 exceed the 2 Mb cost guard), 232 regions, 13.2 Mb; `copy_assign --tie-ms 10`
44,651 contested rows, 2,694 assigned; truth 257,962 hits. Class (unspliced primary): 1,056 decidable molecules, 3,137
negative. Reference 26,424 sequences (2,404 transcripts over two loci left out).
- Bar 0 PASS (1,056). Bar 1 PASS: assigns 430 (0.407).
- Bar 2 FAIL: arm right / primary wrong 4, the reverse 214. Arm 194 / 236, primary 404 / 26.
- Bar 3 FAIL: minus MCL131 4 vs 204.
- Bar 4 FAIL: 2,062 of 3,137 negative molecules assigned (0.657); O2 assigns 344 of the same (0.110).
- Bar 5 PASS.
**FAIL as registered.**
POST HOC (flagged): 171 of the 214 are a scoring asymmetry, not placement errors. The pick is a copy of ANOTHER family
whose span overlaps the truth copy (AS gap 0, 188 of 214 picks are `pre` spans); the registered scorer wants the pick
itself in the truth set, while the primary is right when its position lies in a truth copy's span. Scored the same way
(`EXON_PLACE_RIGHT=locus`, `result.locus.shape*.txt`): class 11 vs 43 (arm 372 / 58, primary 404 / 26; O2 208 vs arm 210
on 210 shared), spliced primaries 14 vs 91. **Still FAIL**, by four to one. On KB3781 the same correction changed
nothing (r1642 check).
Reading: the r1643 unspliced-primary split (12 vs 3) does not hold on a second sample. Exon-space re-alignment is closed
for placement, both references, both classes. Side finding for the catalog: HG002 families MCL32/MCL27 vs MCL373,
MCL65 vs MCL191, MCL12 vs MCL338 hold copies on overlapping spans in different families.
