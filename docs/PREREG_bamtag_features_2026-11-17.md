# Row 1770: do alignment tags and CIGAR details add precision to the chain score? (pre-registered 2026-11-17)

User (2026-11-17): keep looking for empirical precision rules on the 12 samples we have ("I doubt my advisor would allow me to
use more than these"), and look elsewhere, "like in the bam tags". Registered before any tag feature table was made.
Nothing here changes a default.

**What the BAMs carry.** All 12 are minimap2 2.30 `-ax splice:hq -uf --eqx -N 50 -p 0.1` output: `NM ms AS nn ts tp cm s1 s2
de rl SA` (no PacBio read-quality or pass-count tags survive). The chain score already uses `de`, the share of MAPQ 0 reads
and, for seeding, `AS`. Not used so far: `s1` / `s2` (chaining score and that of the best other chain), `rl` (query length
in repetitive seeds), `cm` (minimizers on the chain), `SA` (supplementary alignment), soft clips, and the `=` / `X` / indel
detail of the CIGAR around each splice junction.

**Features, per multi-exon chain, over its full-chain primary records** (`r1770_bamtags/tagfeat.py`; a chain without such a
record gets the sample median, the score already has an indicator for that case):
1. `s2_share` share of records with an `s2` tag; `s2_ratio` mean `s2 / s1` (0 without `s2`)
2. `rl_frac` mean `rl` / read length
3. `cm_dens` mean `cm` / aligned read bases
4. `as_base` mean `AS` / aligned read bases
5. `sa_share` share of records with `SA`
6. `clip5`, `clip3` mean soft-clipped share of the read at the transcript's 5' and 3' end; `clip_big` share of records with a
   clip > 50 nt at either end
7. `jflank_min` over the chain's junctions, the lowest share of records with no mismatch or indel within 8 aligned bases on
   both sides of the junction; `jflank_mean` the mean over junctions

**Test.** 12 samples, six species (human = A119b + HG002). Base = the register 1769 C1 design (decile model of the 23
features + non-GT-AG count + shortest intron), refitted leaving the sample's species out. Arm = the same plus the 11 features
above, each as within-sample deciles plus the value, same fit routine. Both are cut to the number of correct default chains
the shipped r1677 filter keeps on the sample; compare the precision of the kept chains, arm minus base.
- PASS: gain >= +0.005 on at least 10 of 12 samples and no sample below -0.002. A PASS makes a candidate weight file.
- Otherwise FAIL: closed for these tags as chain-level aggregates.
Descriptive, no bar: per feature, precision of the base-kept chains in its lowest and highest decile, per sample; and the gain of
each feature group added alone (post hoc ablation, labelled as such).

**Limits.** All 12 samples are spent for the base score; the test is transfer to a species absent from the fit, and no unseen
sample remains. RefSeq is predicted and incomplete on the apes. `s2` depends on `-N` / `-p` of the alignment: the rule would
hold for BAMs made with the documented minimap2 settings only.

## Result (2026-11-17, job 4478473, `r1770_bamtags/loso.txt`): FAIL by the registered bar
Precision at the r1677 filter's correct count, model fitted without the sample's species.
| sample | base (1769 C1 design) | + all tag features | gain |
|---|---|---|---|
| A119b | 0.2033 | 0.2068 | +0.0034 |
| HG002 | 0.2743 | 0.2925 | +0.0182 |
| OR6737 | 0.4165 | 0.4179 | +0.0015 |
| KB3781 | 0.4385 | 0.4451 | +0.0066 |
| PTR testis | 0.4538 | 0.4554 | +0.0016 |
| PTR AG18354 | 0.4205 | 0.4279 | +0.0073 |
| PPA testis | 0.3726 | 0.3777 | +0.0051 |
| PPA PR00251 | 0.2708 | 0.2824 | +0.0116 |
| PPY testis | 0.2347 | 0.2347 | +0.0001 |
| PPY AG05252 | 0.3791 | 0.3887 | +0.0096 |
| PAB testis | 0.2804 | 0.2791 | -0.0013 |
| PAB AG06213 | 0.4340 | 0.4405 | +0.0066 |
Gain >= +0.005 on 7 of 12 (bar 10), minimum -0.0013 (bar -0.002 met), median +0.0058. Positive on 11 of 12, but under the bar on
five of the six testis samples; the gain sits in the cell lines and HG002.
Post hoc, one group alone (median gain; samples >= +0.005): s2 +0.0004 (0); rl / cm / AS +0.0012 (0); SA + clips +0.0018 (1);
junction flanks +0.0022 (3; PAB testis -0.0035). No single group carries it.
Descriptive (base-kept chains, lowest | highest decile): strong raw contrasts for `rl_frac` (HG002 0.264 | 0.080), `cm_dens`
(0.127 | 0.299), `clip5` (A119b 0.039 | 0.207, HG002 0.036 | 0.242: chains whose reads have no 5' clip at all are rarely correct in
the human samples) and `jflank_mean` (HG002 0.022 | 0.155), i.e. the tags do separate, but most of that is already in the score.
`s2_share`, `sa_share`, `clip_big` are near constant. Closed for these tags as chain-level aggregates; nothing ported.
