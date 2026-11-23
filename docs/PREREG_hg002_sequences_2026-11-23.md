# PREREG: the non-reference sequences on HG002, a second sample with both haplotypes (row 1931, 2026-11-23)

Rows 1928-1930 (KB3781, the gorilla T2T individual): of 107 recurrent non-reference sequences among tied copies, 52
are in the animal's assembled DNA, 8 partly, 47 in neither haplotype anywhere the read aligns; errors, editing,
mosaics, duplicates and a lower-ranked origin were tested and removed. No genomic reads exist to settle the rest.
User: run it on HG002. A different species, library (Kinnex), reference (CHM13, not the sample) and pair of
assemblies (HG002 T2T v1.1 maternal and paternal, `r1776_hg002_room/{mat,pat}.mmi`).

## What is run (nothing looked at yet; HG002's families were used by earlier placement rows, not by this tool)
`HG002_kinnex.t2t.bam` on CHM13: assemble (good pool), families, `o2` with `RUSTLE_O2_SEQUENCES=1`, at today's
defaults; then the reads of the kept sequences aligned to both haplotypes with `-N 200 -p 0.2` and the lookup of row
1930 unchanged (`explain.py`): a sequence is IN THE DNA when >= 0.8 of its reads match >= 0.8 of their variant positions
at the best hit, or at one lower-scoring place; partly; or nowhere.

## The question and how the answer will be read (fixed now)
Share of non-reference sequences explained NOWHERE on either haplotype (KB3781: 47 of 107 = 0.44):
- >= 0.30: the unexplained class is general: it follows the method or long-read libraries, not the gorilla assembly.
- <= 0.15: it was specific to KB3781 (its assembly or library); on HG002 the sequences the tool reports are real.
- between: undecided.
Also reported: sequences per locus, the split of the confirmed ones by haplotype (HG002 differs from CHM13 in both
haplotypes, so confirmed non-reference sequences should be common), and the count check of row 1928 (purity, count,
pairs) with the same script.

## Prediction (mine)
Between 0.15 and 0.30, leaning low: on HG002 most non-reference patterns should be HG002's own variants against
CHM13, which the assemblies hold, so the unexplained share should fall well below the gorilla's; I do not expect it
to vanish. Four of my last five predictions in this line were wrong.

## What happened first (job 4497070)
Assemble 6,796 s, families, `o2` to 18,182 s (5 h in all; the sequence test alone 57 min): 40,627 loci labelled, 490
families, 1,953 sequences kept. Then the stage FAILED: `o2_family_report.py --sequences` read the table with the csv
module and a read-name field of a deep locus exceeded its limit. A defect of the opt-in wired this morning (it must
never fail the stage). Fixed: the column is read by hand; self-test with a 300,000-character field. The scoring scripts
had the same limit. Report and scoring rerun on the files of the run (job 4497529); no haplotype number existed before.

## Result (HG002 on CHM13 against its own maternal and paternal assemblies)
Step 0: 1,953 sequences at 1,591 loci x isoforms, 185 with >= 2 (one locus with 42); 989 with >= 3 reads; 1,455 with a
non-reference pattern.

| non-reference sequences (689 with a substitution and reads on a haplotype) | HG002 | KB3781 (rows 1929, 1930) |
|---|---|---|
| in the DNA at the read's best hit | 459 (0.67) | 47 (0.44) |
| at ONE lower-scoring place | 58 (0.08) | 5 (0.05) |
| partly | 74 (0.11) | 8 (0.07) |
| **nowhere on either haplotype** | **98 (0.14)** | **47 (0.44)** |

Per variant position: 0.989 match at the best hit (KB3781 0.71).
- **By the reading fixed beforehand (<= 0.15): the large unexplained class was specific to KB3781.** On HG002 three
  quarters of the non-reference sequences are in the individual's assembled DNA and one in seven is nowhere. 0.142 is
  just under the line; my prediction (0.15-0.30, leaning low) missed by that margin.
- Count check with row 1928's script: purity 830 of 902 (0.92); count equal to the haplotypes' at 574 of 678 loci x
  isoforms (0.85). Its pair measure again says "not separated" for 1,232 of 1,318 pairs, but 861 of those pairs come
  from the one locus with 42 sequences, and the site-by-site lookup above contradicts it: that pair measure (copies by
  overlapping hit regions, allele by cost) is too coarse and should not be quoted, here or in row 1928.
- What this says about KB3781: its 44% is not a property of the method. Something about that sample differs: the
  gorilla assemblies at multi-copy loci, or the fibroblast library. Not resolved.
- What can be said about the opt-in output now: on HG002 a recurrent non-reference sequence among tied copies is in
  the individual's DNA three times out of four (517 of 689), partly in one of ten, absent in one of seven.
- Cost: the sequence test took 57 minutes on this library (python, one group at a time); the three stages 5 hours.
