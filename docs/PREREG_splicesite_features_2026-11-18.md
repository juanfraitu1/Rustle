# Row 1789: does splice-site sequence strength add precision? (pre-registered 2026-11-18)

Goal (user, standing): more empirical precision rules on the 12 samples; no annotation at run time. Rows 1770 (BAM tags) and
1771 (read ends) each FAIL; together 8 of 12. Registered before any splice-site feature was computed. Nothing here changes a default.

**Already in the score:** whether an intron is non-canonical (`noncanon_capped2`), junction wobble, minimum intron support.
Not in it: how well the sequence around each junction matches a splice site, and whether a stronger site lies a few bases away.

**Features, per multi-exon chain, from the genome sequence only** (`r1789_splicesites/ssfeat.py`):
- Position weight matrices built per sample from the distinct introns of that sample's own default assembly (no annotation):
  donor = 3 exonic + 6 intronic nt, acceptor = 20 intronic + 3 exonic nt, log2 odds against 0.25, pseudocount 1.
1. `don_min`, `acc_min` lowest donor / acceptor score over the chain's introns
2. `pair_min`, `pair_mean` lowest and mean of donor + acceptor per intron
3. `ppt_min` lowest pyrimidine share of the 20 intronic nt before the acceptor
4. `don_alt`, `acc_alt` largest excess, over the chain's introns, of the best-scoring other GT (AG) within 10 nt over the used site
   (0 when none is better)

**Test** (same as rows 1770 and 1771). Base = register 1769 C1 design refitted without the sample's species; arm = the same
plus the 7 features, each as within-sample deciles plus the value. Both cut to the number of correct default chains the r1677
filter keeps; precision of the kept chains, arm minus base.
- PASS: gain >= +0.005 on at least 10 of 12 samples and no sample below -0.002. A PASS makes a candidate weight file.
- Otherwise FAIL.
- Declared second arm, same bar: base + BAM tags (row 1770) + read ends (row 1771) + these 7. Its samples have been used for
  two failed rows already, so a PASS there is weak evidence and at most an opt-in candidate; the user decides.
Descriptive, no bar: lowest | highest decile precision per feature.

**Limits.** All 12 samples are spent; transfer across species only. RefSeq as truth: a real unannotated transcript with good
splice sites counts as wrong, so the feature can only help where wrong chains have weak sites.

## Result (2026-11-18, job 4480449, `r1789_splicesites/loso.txt`): FAIL, both arms
| arm | samples with gain >= +0.005 | minimum | median |
|---|---|---|---|
| splice sites alone | 4 of 12 | -0.0017 | +0.0022 |
| splice sites + tags + read ends | 9 of 12 | +0.0012 | +0.0121 |
| tags + read ends (rows 1770, 1771, for reference) | 8 of 12 | +0.0002 | +0.0103 |

- Gains of splice sites alone: HG002 +0.0123, PPA PR00251 +0.0103, PPY AG05252 +0.0091, PTR AG18354 +0.0059; every testis sample
  and both gorillas under +0.003; PAB testis -0.0017.
- Three-set arm misses the bar by one sample: under +0.005 are PTR testis +0.0036, PAB testis +0.0028, PPY testis +0.0012.
  A119b moves from +0.0047 to +0.0072.
- Descriptive: chains in the lowest decile of mean donor + acceptor score are almost never correct in the human samples
  (HG002 0.001, A119b 0.041 against baselines 0.274 and 0.203) and only mildly worse in the apes (for example OR6737 0.396
  against 0.416). A stronger GT / AG within 10 nt lowers precision by 0.02 to 0.11 in every sample but PAB testis.
- Nothing ported. Third feature family in a row whose gain is in cell lines and HG002 and near zero in ape testis.
