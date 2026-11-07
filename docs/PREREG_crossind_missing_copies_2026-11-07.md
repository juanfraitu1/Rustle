# PREREG rows 1578-1579: does a second gorilla's Iso-Seq on Jim's T2T show copies the reference lacks?

## Why
Advisor (2026-11-07, via the user): any sample aligned to a T2T genome from another individual misses many copies
(reads of copies the reference lacks have nowhere to go); only Jim's own fibroblast (KB3781) is exact on Jim's T2T, the
OR6737 testis sample is not. The user's view: SNPs and small changes yes, whole missing paralogs few. Direction tested:
copies the INDIVIDUAL has and the REFERENCE lacks (RNA cannot tell a reference copy absent in the individual from an
unexpressed one).

## Provenance (read first)
Both decision rules were stated to the user in the session before the job that decided them was run; this file was
written after the results (2026-11-08). The commit does NOT prove precedence. Post hoc parts are marked.

## Data
OR6737 testis `GGO_mm.bam` (4.40 M primaries) and KB3781 fibroblast `GCA_029281585.2_flnc_mm.bam` on GGO
(GCF_029281585.2), both `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1`. KB3781 subsampled by read-name hash to
keep 0.35 (4.73 M primaries) = depth match. Tissues differ (testis vs fibroblast): NOT the matched-tissue design the O3
ledger asks for. No third HiFi Iso-Seq gorilla is public: DIAN / OMOYE (PRJNA797478, LCL) are Sequel subreads only
(no BAM, no CCS possible); PRJNA1443378 is ONT PCB109 iPSC cDNA.

## Row 1578: divergence piles (`bench/crossind/piles.py`)
Primary reads clustered by overlapping landing; clusters >= 10 reads. Read multi-copy = a secondary with AS >= 0.9 x
primary AS; cluster multi-copy when >= 30% of its reads are. Pile = >= 10 reads at de > T and >= 25% of the cluster.
**Rule (stated before results):** the claim "a lot of missing copies" FAILS if OR6737's multi-copy pile rate minus
Jim's (depth-matched) is < 5 percentage points at T = 0.01.

## Row 1579: three-haplotype test (`bench/crossind/threehap.py`)
A diploid shows <= 2 haplotypes at a locus. Reads de <= 0.02; site = alt base on >= 3 reads and >= 15% of covering
reads (indels ignored); a site pair within 2 kb has >= 3 haplotypes when each has >= 3 reads and >= 10% of the
spanning reads; locus flagged at >= 2 site-disjoint such pairs. Jim on his own genome = the noise floor.
**Rule (stated before the run):** claim FAILS if OR6737's flagged rate minus Jim's is < 5 points on loci expressed in
both samples. (The single-copy stratum was dropped as a null before the run: reads of a reference-absent copy have no
secondary, so its landing locus looks single-copy.)

## Results (2026-11-07/08, `/anvme/workspace/v138da13-rustle/crossind/`)
Row 1578, T = 0.01, as registered: OR6737 45/601 multi-copy clusters (7.5%) vs Jim 8/438 (1.8%) = +5.7 points:
**claim NOT refuted as registered (borderline).** Post hoc, loci expressed in both samples: 8/316 (2.5%) vs 1/316 (0.3%)
= +2.2; the excess sits in testis-only loci (individual and tissue confounded there). OR6737 reads are noisier
(read de p95 0.0285 vs 0.0085); at T = 0.005 OR6737 single-copy clusters pile 3x Jim's too (SNPs + error).

Row 1579, 7,056 shared loci: OR6737 499 (7.1%) vs Jim 324 (4.6%) = +2.5 points: **claim FAILS as registered.**
OR-only flags 361, Jim-only 186.

Post hoc: strong candidates = OR-only three-haplotype flag AND a 1% pile: 12 loci (`strong.tsv`), Jim 0 independent
pairs at each (thin Jim depth at some). 7 of 12 lie in OR6737 de novo families (`or.fam.copies.tsv`, 446 families):

| locus | indep. pairs | family | genes |
|---|---|---|---|
| NC_073224.2:117.29-117.53 Mb | 60 | MCL5 (15 copies) | NBPF / NOTCH2NL-like |
| NC_073242.2:18.35-18.60 Mb | 22 | MCL67 (5) + MCL3 (16) | IGHV3-11-like |
| NC_073242.2:22.13 Mb | 18 | MCL7 (10) | NACA-like |
| NC_073230.2:74.82-74.95 Mb | 14 | MCL35 (7) + MCL271 | INTS4-like / CCT6A |
| NC_073235.2:31.68 Mb | 3 | MCL310 (2) | SMPD4-like |
| NC_073224.2:223.16-223.94 Mb | 2 | MCL5 + MCL224 | NBPF / FAM231 |
| NC_086017.1:85.70 Mb | 2 | MCL197 (3) | ZNF717 |

No family (weaker): TMT1A (22 pairs), LOC129527636 (17), APOOL/SATL1, GLIPR1/KRR1, PKNOX1/WDR4 region.

## Reading
A second individual brings expressed copies Jim's T2T lacks, concentrated in classic copy-number-variable families
(NBPF, IGHV): about 7 of 446 families, 12 loci, not "a lot". Not shown: testis-only loci (no Jim comparison), copies
present in Jim and absent in OR6737 (untestable from RNA), copy-level confirmation (no IGV / DNA check of the 12).
Do not quote as "multi-copy families do not differ between individuals".
