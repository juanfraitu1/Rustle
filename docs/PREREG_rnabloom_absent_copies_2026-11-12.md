# PREREG REFERENCE-FREE ASSEMBLY FOR ABSENT COPIES, RNA-Bloom (register 1661, 2026-11-12, before any output)

Question (user 2026-11-12, after the LRGASP discussion: "lets try the step 0 you suggest"): does a reference-free long-read
assembler give one clean contig per expressed gene copy that a draft genome lacks, so that "a contig that does not fit the
draft" is a better absent-copy call than the read divergence pile of r1539 (recall 13 of 17, precision 11 of 24, 1 of 7
present copies called)? O3 measurement, descriptive, opt-in only, not for `advisor_share`.

Step 0 = a CEILING test on the gorilla truth panel only. Data: `/anvme/workspace/v138da13-rustle/nont2t/gorilla/`
(`reads.fq` 13,041 KB3781 reads, `read_origin.tsv`, `truth2.tsv`, `panel.genes.tsv`; draft index
`annot/nont2t/kamilah.splice.mmi`). The reads were chosen with the T2T (panel regions), so off-panel reads are absent and any
precision here is an upper bound. The reads are Jim's, Jim's T2T (`GGO.fa`, plus the `KB3781.mat` / `.pat` haplotypes) is truth
for sequence.

Method, fixed now:
- RNA-Bloom 2.0.1, `rnabloom -long reads.fq -lrpb` (PacBio preset), every other option default, not stranded. No labels.
- Contigs aligned with `minimap2 -ax splice:hq --eqx -Y -N 50 -p 0.1` to GGO.fa, the two haplotypes, and the draft.
- ARTIFACT contig: best primary over the three truth genomes has de > 0.002 or covers < 0.95 of the contig. Scored by
  sequence, never by re-placing reads on contigs.
- A contig's gene: the panel gene whose exons hold >= 50% of the contig's aligned bases on GGO.fa. TIE RULE for identical
  copies: every other record of the contig with AS >= the primary's AS gives its gene too; a contig is STRICT for a gene when
  it has that gene and no tied gene, TIED otherwise.
- Blind CALL (draft only): the contig is unmapped on the draft, or its primary has de > 0.005 (r1539's TAU, not retuned) or
  covers < 0.90 of the contig.
- Scorable genes: >= 10 origin reads: 17 missing, 7 present (the r1539 sets). A missing gene is FOUND when a non-artifact
  called contig has it (strict or tied). A present gene is a FALSE ALARM when a called contig has it strictly.

Bars (all four, else STOP): found >= 13 of 17; present false alarms <= 1 of 7; artifact rate <= 0.10 of contigs; strict
found >= 9 of 17 (the part piles cannot give: one contig per copy).
Stop rule: PASS -> register the same recipe on the 12 strong OR6737 loci (no truth there; reach and agreement only) and ask
the user. FAIL -> close; the OR6737 half is NOT run. Either way report the present contigs' draft identities (the draft is
another individual, so present copies diverge too: the allele confound) and the calls outside the panel genes.
Limits: null of 7 genes; panel reads only; one assembler setting.
Script: `bench/nont2t/bloom_contigs.py` (selftest), workspace `r1661_rnabloom/`.

## Result (2026-11-12, job 4468530, workspace `r1661_rnabloom/`): FAIL on all four bars, STOP
RNA-Bloom ran in 79 s and wrote 614 contigs from the 13,041 reads.
| bar | result |
|---|---|
| artifact rate <= 0.10 | 429 of 614 = 0.699 |
| missing found >= 13 of 17 | 8 |
| strict found >= 9 of 17 | 6 |
| present false alarms <= 1 of 7 | 4 |
- The contigs are full length on Jim's own genome (coverage 1.0 at the median, 28 of 614 under 0.95) but not accurate: de
  against the best of GGO.fa / mat / pat is 0.004 at the median and above 0.002 for 421 contigs, where the HiFi reads
  themselves sit near 0.001. Post hoc, the artifact rate at a looser de bar is 0.417 (0.005) and 0.21 (0.01): the verdict
  does not depend on the bar.
- So the blind call (de > 0.005 on the draft) fires on the assembler's own error: 448 of 614 contigs called, 371 of them
  artifacts, 322 on no scorable panel gene, 4 of 7 present genes called.
- Not found: LRPAP1 and its LOC copies, the JMJD7-PLA2G4B LOC copies (the near-identical families).
Reading: reference-free assembly of HiFi Iso-Seq reads loses accuracy relative to the reads, in the families where copies
differ by less than the error it adds. It is worse than the read divergence pile on every count (pile 13 of 17, 1 of 7).
Same outcome class as IsoCon (r1478, r1515). The OR6737 half is not run. Reference-free assembly is closed for O3.
Limits: one assembler, one setting (PacBio preset); whether the contig error is consensus error or copy blending was not
separated.

## Follow-up, descriptive (2026-11-12, user asked: consensus error or copy blending?)
`r1661_rnabloom/errsrc.py`: every difference between a contig and its best alignment on the sample's own genomes
(13,858 sites), tested by whether the contig's 31-mer around the site occurs in the panel reads. In >= 2 reads = real
sequence that belongs somewhere else (another copy: a blend); in exactly 1 read = one read's error carried into the
contig; in 0 reads = made by the assembler. Control, matching positions: 0.931 in >= 2 reads, none in 0.
| difference | sites | 0 reads | 1 read | >= 2 reads |
|---|---|---|---|---|
| substitution | 8,131 | 255 | 2,171 | 5,705 (0.70) |
| insertion | 4,902 | 264 | 4,063 (0.83) | 575 |
| deletion | 825 | 14 | 559 (0.68) | 252 |
| all | 13,858 | 533 (0.04) | 6,793 (0.49) | 6,532 (0.47) |
Per artifact contig (427): mostly read-supported differences 214, mostly unsupported 86, mixed 119, none scorable 8.
Reading: both, split by type. Substitutions are mostly blending of near-identical copies. Indels are mostly a single
read's error left unpolished (thin contigs: `-lrrd 3`). The assembler invents little (4%). Raising the read depth would
remove the single-read errors but not the blends, which are half the artifact contigs on their own.
Limit: a site in >= 2 reads could also be a recurrent HiFi error (homopolymers); not separated.

## Follow-up 2, descriptive (2026-11-12, user asked: blend or recurring HiFi error?)
`r1661_rnabloom/kscan.py` (exact, every overlapping window; selftest): is the contig's 31-mer at each differing site
present anywhere in the sample's genomes (GGO.fa, KB3781 mat, pat)? Control (matching positions in >= 2 reads): 566 of 566.
| sites | n | in a genome (another copy) | in no genome | of those, in a run of >= 4 identical bases |
|---|---|---|---|---|
| in >= 2 reads, substitution | 5,705 | 4,606 (0.81) | 1,099 | 123 |
| in >= 2 reads, insertion | 575 | 301 (0.52) | 274 | 128 |
| in >= 2 reads, deletion | 252 | 179 (0.71) | 73 | 10 |
| in >= 2 reads, all | 6,532 | 5,086 (0.78) | 1,446 | 261 |
| in exactly 1 read | 6,793 | 342 (0.05) | 6,451 | 1,464 |
| in 0 reads | 533 | 35 | 498 | 61 |
Of all 13,858 differences: 5,463 (0.39) are another copy's real sequence (blend), 6,451 (0.47) a single read's error,
1,446 (0.10) recur in reads but are in no genome, 498 (0.04) made by the assembler. Per artifact contig (427): mostly
another copy 122, mostly not in any genome 104, mixed 193, no site 8.
The "recurring, in no genome" class is not mostly homopolymer error (261 of 1,446). Limit: a 31-mer within 15 bp of an
exon junction is never in the genome, so part of that class is junction-adjacent sites, not errors (control windows did
not span junctions); RNA editing would land there too. So 0.78 is a lower bound for blending among read-supported sites.
