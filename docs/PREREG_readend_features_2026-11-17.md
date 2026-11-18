# Row 1771: do read-end features (cap G, clip content, end sites shared across the locus) add precision? (pre-registered 2026-11-17)

User (2026-11-17): "can we try the read-end direction?", after row 1770 (BAM tags as chain aggregates, FAIL). Registered before
any feature table with truth labels was made. Nothing here changes a default.

**Seen before registering (no truth labels):** the soft-clipped 5' end of the reads is mostly empty or a single `G`
(PTR testis, HG002, PPY AG05252: first 13-40 thousand spliced primary reads of one contig), i.e. the non-templated G that
reverse transcription adds at a capped 5' end. Row 1770 saw, with labels, that chains whose reads carry no 5' clip are rarely
correct in the human samples (0.04 against 0.22-0.28); that observation is the reason for this row and is not independent of it.

**Already in the score:** share of full-chain reads within 10 bp of the commonest 5' / 3' end, polyA signal and A-richness at
the 3' end, container and partial reads. Not in it: what the clipped bases are, and whether other chains' reads end at the same site.

**Features, per multi-exon chain, over its full-chain primary records, oriented by the model strand** (`r1771_readends/endfeat.py`):
1. `cap_g` share of reads whose 5' clip is 1-4 nt, all G; `clip5_none` share without a 5' clip; `clip5_other` the rest
2. `clip3_none` share without a 3' clip; `clip3_a` share with a 3' clip of >= 3 nt that is >= 80% A
3. `tss_other`, `tes_other` log1p of spliced primary reads of OTHER chains whose end on that side lies within 10 bp of the chain's
   commonest 5' / 3' end; `tss_cap_other` share of all spliced reads starting at that 5' site (the chain's own included) that carry the cap G
4. `end5_iqr`, `end3_iqr` log1p of the interquartile range of the chain's read ends
5. `first_exon`, `last_exon` log length from the commonest end to the first / last splice site

**Test** (same as row 1770). 12 samples, six species. Base = register 1769 C1 design refitted without the sample's species;
arm = the same plus the 12 features, each as within-sample deciles plus the value. Both cut to the number of correct default
chains the r1677 filter keeps on the sample; precision of the kept chains, arm minus base.
- PASS: gain >= +0.005 on at least 10 of 12 samples and no sample below -0.002. A PASS makes a candidate weight file.
- Otherwise FAIL.
Descriptive, no bar: lowest | highest decile precision per feature; gain of each of the five groups alone (post hoc, labelled).

**Limits.** All 12 samples are spent; transfer across species only. The cap G depends on the library protocol (template
switching) and on soft clipping being kept in the BAM (`-Y`); a rule built on it holds for such data only. RefSeq 5' ends on the
apes are predicted, but scoring is by intron chain, so truth does not depend on them.

## Result (2026-11-17): FAIL, decided on 11 samples (`r1771_readends/loso11.txt`, job 4478528)
Deviation: HG002's BAM exists only on the atuin file server, which read at ~6 MB/s that morning, so the comparison was first run
without HG002 (and without HG002 in any training set). Gain >= +0.005 on 5 of 11 (KB3781 +0.008, PTR AG18354 +0.010, PPA testis
+0.006, PPY AG05252 +0.011, PAB AG06213 +0.011; PPA PR00251 +0.0050 counted), minimum -0.0005, median +0.0050. With at most one
more sample the bar of 10 of 12 cannot be reached: FAIL whatever HG002 gives. The 12-sample run (job 4478498) is appended below when done.
Post hoc, one group alone (median; samples >= +0.005): 5' clip content -0.0002 (0); 3' clip content -0.0001 (0); shared end sites
+0.0016 (4); end spread +0.0008 (0); terminal exons +0.0033 (3). The cap G adds nothing on top of the score.
Descriptive: the cap G separates strongly only in A119b (0.037 | 0.135); the consistent raw contrasts are terminal exon lengths
(long first exon: 0.05-0.43 against 0.23-0.48; short last exon low) and 3' end spread. As in row 1770 the gain sits in the cell lines.

**12-sample run (job 4478498, `loso.txt`), as registered: FAIL.** Gain >= +0.005 on 7 of 12, minimum -0.0011, median +0.0056.
HG002 +0.0139; A119b -0.0002, PTR testis +0.0016, PPY testis -0.0001, PAB testis -0.0011, OR6737 +0.0041. Post hoc groups alone
(median; samples >= +0.005): 5' clip content -0.0003 (0), 3' clip content -0.0001 (0), shared end sites +0.0024 (2), end spread
+0.0004 (0), terminal exons +0.0034 (4).

**Post hoc arm, registered 2026-11-17 before running (`loso_both.py`, `both.txt`): rows 1770 and 1771 together.** Both feature
sets fail alone with similar gains on the same samples (cell lines, HG002). Question: do the gains add? Same protocol and bar
(+0.005 on >= 10 of 12, none below -0.002), base = 1769 C1 design, arm = base + 11 tag features + 12 read-end features. Post hoc
on spent samples: a PASS would be weak evidence and makes at most a candidate.

**Result of the post hoc arm (job 4478626, `both.txt`): FAIL by the bar, gains roughly add.** Tags + read ends: gain >= +0.005 on
8 of 12 (bar 10), minimum +0.0002, median +0.0103. HG002 +0.026, PPY AG05252 +0.017, PTR AG18354 +0.015, PPA PR00251 +0.015,
PAB AG06213 +0.014, KB3781 +0.012, PPA testis +0.008, OR6737 +0.006; under the bar: A119b +0.0047, PTR testis +0.0039, PAB testis
+0.0011, PPY testis +0.0002. No sample loses. The four failures are testis samples; every cell line and HG002 pass.
