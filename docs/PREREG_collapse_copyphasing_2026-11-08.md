# PREREG: does copy_phasing's allele-slot certificate detect collapsed copies on GRCh38? (register 1583, 2026-11-08, before any number)

## Question
On a non-T2T reference, a collapsed tandem copy's reads land on its surviving relative at that relative's own identity
(r1537), so there is no per-read sign. r1538 counted variant sites per kb (S) and the share of reads on a second haplotype (H).
Those blind statistics failed: 4 of 13 autosomal loci, against null loci at heterozygosity levels.

No test so far has counted haplotypes beyond two with linked reads. In a diploid, a locus holding two collapsed copies carries up
to 4 haplotypes in 2 allele slots. That is what copy_phasing's certificate flags, and what it flagged in OR6737 (row 1582).

## Data (r1537/r1538 products, human A119b testis on GRCh38)
- BAM `/anvme/workspace/v138da13-rustle/nont2t/human/nt.bam`, aligned with `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1`.
- Genome `~/annot/nont2t/grch38.fa`.
- Loci and labels from `psv.loci.tsv`: 16 scorable missing copies (absorbing loci) and 73 null loci. Origin labels come from
  `read_origin.tsv` and are used for scoring only.

## Run
- Families are the pipeline's own de novo families on GRCh38: `x.fam.copies.tsv`, 60 families, 319 copies. This is blind; the
  catalog is not built from truth.
- `bench/o3_graph/real_net.py nt.bam x.fam.copies.tsv cp` builds the >= 2-copy families and the net, using the r1483 convention.
- The binary is the same frozen `~/o2ties_real/copy_phasing.port` that row 1582 read:
  `copy_phasing.port run --bam nt.bam --copies cp.copies.tsv --genome grch38.fa --net cp.net.tsv --out cp.A --threads 24`.
  The run uses defaults (cc clusterer); no `--semi` or `--guard`, which affect only O2.

## Definitions (as row 1582)
- **Testable.** At least one `cp.copies.tsv` copy overlaps the locus window (`chrom:start-end` of `psv.loci.tsv`).
- **Hit.** A `flagged=1` row in `cp.A.candidates.tsv` or any row of `cp.A.putative.tsv` whose home copy overlaps the window.
- **Single-origin null locus.**
  - It is a `kind=null` row.
  - Count the primary reads in its window that carry an origin label.
  - The row's own gene must be the top origin, and no other gene may hold >= 10 % of those reads.
  - This is the null r1538 "after the numbers" item 2 named as the next row.
- **Strata.** Autosomal and chrX/chrY are reported separately. Every copy gets 2 slots, so a haploid X/Y locus needs 2 extra
  haplotypes to certify; X/Y is reported but not barred.

## Bars (autosomal)
1. Hits at >= 50 % of the testable missing copies, counted per copy as in r1538. Each collapsed / absent split is also reported.
2. Hits at <= 25 % of the testable single-origin null loci. The full null (all `kind=null`) is reported, not barred.
3. Neither bar is decided on fewer than 5 testable units. With fewer, the verdict is UNDERPOWERED and nothing is closed.

**PASS:** the certificate is a blind collapse detector on non-T2T references, and an HG002 / second-genome row comes next.
**FAIL on bar 1:** the linked-haplotype signal is too weak at RNA depth for these loci. This closes the collapsed regime on reads
for the certificate as well as for S/H; collapsed copies then need a better assembly or DNA.

## Expected (written now)
r1537 bounds the missing-to-absorber divergence at <= 0.002 for 19 of 21 copies. copy_phasing needs >= 2 significant private
columns and >= 6 reads per haplotype. The loci with labelled purity 1.00 (GOLGA6L9, NPIPB1P, GOLGA8M, GOLGA6L17P) are the
likeliest hits. The loci below 10 % of the pile (GOLGA6L1, NPIPA3) are the likeliest misses.

Scorer: `bench/nont2t/collapse_cp.py`.

## RESULT r1583 (2026-11-08): PASS as registered; the labelled check halves it
- **Run.** Slurm 4452735, 158 s. 60 families, 225,601 net rows, 319 copies; 94 certified and 23 putative copies. Outputs are in
  `/anvme/workspace/v138da13-rustle/nont2t/human/cp/`: `collapse_cp.txt` for the verdict and `purity.txt` for the labelled check.
- **Results by stratum.**

| stratum | testable | hit |
|---|---|---|
| autosomal missing | 13 | 10 (0.77) |
| … of which collapsed | 7 | 5 |
| … of which absent | 6 | 5 |
| autosomal null, single-origin | 30 | 7 (0.23) |
| autosomal null, all | 51 | 15 (0.29) |
| chrX/chrY missing | 3 | 0 |
| chrX/chrY null | 16 | 0 |

- **Bars.** Bar 1 passes (0.77 >= 0.50) and bar 2 passes (0.23 <= 0.25), so the verdict is PASS.
- **Per copy.** Copies inside absorbing windows are hit at 24/47 = 0.51; single-origin null windows at 7/53 = 0.13, which matches
  the catalog background of 43/319 = 0.135. Window sizes are comparable: median 108 kb for absorbing windows, 94 kb for null.

### After the numbers (labelled; the verdict stands as registered)
1. **Many hits are someone else's haplotype.** For each certified copy homed in an absorbing window, I looked up the origin
   labels of its reads (`purity.txt`).
   - A copy whose reads are mostly (>= 0.5) the missing gene's own exists at 5 of the 13 loci:
     - NPIPB1P: 1.00, 128 reads.
     - GOLGA8DP: 0.99, 394 reads.
     - GOLGA6L4: 0.97, 32 reads.
     - NPIPB9: 0.58, 158 reads.
     - NPIPB4: 1.00, but only 9 reads; its large candidates are 0.19.
   - The other hits are certified copies of neighbouring genes, which the large window overlaps:
     - GOLGA6L1's hit is GOLGA8DP's copy.
     - NPIPB15's hit is mostly NPIPB14P.
     - NPIPA3's and NPIPB11's hits are unlabelled reads or other NPIP genes.
   - **Labelled-correct: 5 of 13** (collapsed 3 of 7, absent 2 of 6). That is no better than r1538's 4-6 of 13.
2. **The misses are where r1538 succeeded.** GOLGA6L9, GOLGA8M and GOLGA6L17P have small windows and a clean second haplotype
   (r1538 labelled purity 1.00). They draw no certificate. The likely reason is Theorem A: if the surviving copy is homozygous, the
   missing copy's haplotype fills its free second slot. This is not verified per locus.
3. **The two methods catch different copies.** The certificate catches collapses where the survivor is heterozygous (GOLGA8DP,
   GOLGA6L4). The S/H pileup catches pure second haplotypes on a homozygous survivor. A combined rule is a new row, not this one.
4. **Null hits are not necessarily false.** 7 single-origin null loci carry a certificate. A119b could carry a copy that GRCh38
   collapsed and that CHM13 also lacks; the panel cannot tell.
5. **X/Y.** 0 of 3 missing copies are hit. Two slots on a haploid locus need two extra haplotypes; this was predicted.

**Reading.** The registered blind rule passes, but the windows are large and sit in NPIP/GOLGA clusters dense with certificates,
and that is what lifts it. Only about 5 of 13 hits carry the missing copy's own reads, which is the same reach as the pileup row.
The collapsed regime is partly visible from reads: 3 of 7 by certificate and about 4 of 7 by labelled purity. No single blind
rule reaches half yet. A next row would combine the two rules, with per-copy home matching instead of window overlap; it is not
run.
