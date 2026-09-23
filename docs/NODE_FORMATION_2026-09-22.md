# Node formation — empirical rules for FP/FN (§6z4, 2026-09-22)

Session goal: *look for empirical rules for reducing false positives and false negatives by performing
better node formation.* Every claim below is a measurement on the chr20 de novo assembly
(`a119b_polished.gtf`, 1,779 loci / 5,844 transcripts; chr16 where stated), register rows **1049–1054**
plus the r1041 chain. Follows `docs/NODE_DEFINITION_2026-09-22.md` (§6x3–§6x7), which had already closed
splitting, boundary pull-in and the local denominator.

## TL;DR — what the rules turned out to be

| candidate node rule | verdict | the number |
|---|---|---|
| drop / trim single-exon loci (intronic pileups) | ⛔ not a lever | 25.6% of loci, **7.1% of nodes, 0 cross-family edges** (r1049) |
| pick the transcript that maximises corroborated exon coverage | ⚠ candidate, unscored | better in 60.4% of loci, median gain **15%**, at **0.23×** the rep's reads (r1050) |
| complete the rep (union / corroborated union) to fix FN | ⛔ wrong sequence | rep incompleteness predicts eviction (**AUC .605**, survives depth) **but** the missing exons are non-homologous: under the partner in **48%** of evicted loci, median **0.00** (r1051/r1052) |
| fuller exonic denominator | ⛔ evicts more | union **−47**, corroborated **−23** net admissions (r1053) |
| smaller (ORF) exonic denominator | ⛔ hubs | edges **4,284 → 8,378**, largest component **330 → 480** (r1054) |
| corroborated union rep, end to end (r1041 chain) | ⛔ U loses vs B on every axis; C2/C3 stopped (r1055) | §5 |

⭐**The rule that survives is negative and precise: at the DNA level the node's exon model is only a
coverage DENOMINATOR, and the shipped single-chain rep already sits at the usable point between
eviction (fuller) and hubs (smaller).** The population r1002 identified as false negatives is real and is
concentrated in high-isoform loci (r1051), but its remedy is the asymmetric containment escape on the EDGE
(§6x4/§6x5, shipped as `--min-cov-shorter`), not a node change.

## 1. Single-exon loci are not a node lever (r1049)

They are 8 kb genomic spans with no splicing evidence (r1043) — 46.7% inside a multi-exon gene at exonic
fraction 0.00, 34.3% intergenic. On chr16 de novo they are 25.6% of loci but **7.1% of graph nodes**, and
among referee-labelled edges they carry **0 cross-family edges** (all 13 are spliced↔spliced). Only 12.4%
of chr20's sit inside a spliced locus's intron. With r356 (dropping them costs precision *and* recall),
there is nothing to gain on either side.

## 2. Rep incompleteness is a real FN mechanism (r1051) — but completion is the wrong remedy (r1052)

Admitted loci have median rep completeness **0.709**, evicted **0.580**; admission rate rises
**0.29 → 0.35 → 0.48** across completeness bands; **AUC 0.605**, and within read-depth tertiles
**0.564 / 0.542 / 0.677** (strongest where isoform diversity is highest). Depth (**0.441**) and isoform
count (**0.404**) are *anti*-predictive: deeper, more-isoform loci are evicted more.

But the exons the rep dropped are, for evicted loci, mostly **not** where the partner aligns — under the
partner's best alignment in **48%** of cases, median **0.00**, ≥0.5 in **10%** (admitted: 87% / 0.19 / 35%).
A union adds locus-specific sequence that enlarges the denominator without adding aligned bases —
register 303's failure, exon by exon — and corroboration does not select for homology.

## 3. Every denominator direction, measured (r1053, r1054)

Offline ceiling, chr20, numeric gates only (absolute counts exceed the shipped graph; the deltas are the
point):

| exonic denominator | admitted / 1,067 | edges | largest component |
|---|---|---|---|
| union | 882 (**−47** vs rep) | | |
| corroborated ≥2 tx | 906 (**−23**) | | |
| **single-chain rep (shipped)** | **929** | 4,284 | 330 |
| longest ORF (CDS proxy, 0.16× the rep) | 1,033 | **8,378** | **480** |

Fuller → eviction. Smaller → hubs (r913). The shipped rep is the usable point.

## 4. Rep selection (r1050) — the one candidate left, and why it is parked

Choosing the *existing* transcript that maximises corroborated coverage cannot inflate a denominator, but
it helps in 60.4% of multi-transcript loci by a median **15%** while taking a transcript with **0.23×** the
rep's read support (208/305 under 50%). Its family-level effect needs the full path; sequence after §5.

## 5. The r1041 chain — corroborated union rep, end to end

**Stopped by the user after arms B and U** (C2 killed mid-run, C3 never started; register 1055). On
`gw_family_catalog` chr20, scored on `copies.tsv` against the chr20 protein referee (Soto has no scoreable
chr20 family):

| arm | families | copies | mean copy span | sens | prec | F | collapsed |
|---|---|---|---|---|---|---|---|
| **B** shipped single-chain rep | **127** | **800** | 5.0 kb | **0.158** | **0.667** | **0.255** | **0** |
| U plain union rep | 84 | 579 | 18.5 kb | 0.132 | 0.488 | 0.207 | 5 |

⛔ U loses on every axis, as the refined prediction said (the union is 3.7× longer and the length is
denominator, not homology — r1052). The corroborated arms C2/C3 are **unmeasured on this driver**; their
prediction ("lose less, do not beat B") stands unverified. Rule 5 therefore rests on B vs U here plus the
offline r1053.

## Reproduce

```sh
# r1051-r1054 all read /mnt/linuxdisk/tmp/se/{loci.json,all.paf,all.graph.tsv} + a119b_polished.gtf;
# the probes are inline in the §6z4 session and the register rows carry every number.
python3 bench/edge_probes.py poset ...      # comparability/containment machinery reused by r1051
target/release/family_score --clusters X.clusters.tsv --gff chr20.genes.gff --soto TRUTH --chrom chr20
```
