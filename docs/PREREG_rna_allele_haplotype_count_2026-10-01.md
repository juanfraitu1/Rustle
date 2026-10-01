# Pre-registration: counting haplotype copies from RNA alone ("transcripts that are alleles"), gorilla KB3781, 2026-10-01

Written before any truth table was built and before any RNA call was made. Un-parks `project_rna_haplotype_count_parked` (30 Sep) with
the advisor's framing (1 Oct): the reads carry enough information to tell how many haplotype copies of each gene copy an individual has,
because some distinct transcripts are alleles of one copy rather than separate copies.

## Question

From the RNA of one individual, with no DNA of that individual: for each copy of a gene family, can we tell that the copy is carried
on both haplotypes (two distinguishable alleles in the reads), and can we count the haplotype copies of a family as a lower bound?

## What RNA can and cannot say (stated before any data)

- A copy whose reads split into two phased versions at positions where all of its paralogs agree is carried on both haplotypes.
- A copy whose reads show one version is either on both haplotypes with identical alleles (homozygous) or on one haplotype only
  (hemizygous). Expression level is not dosage, so RNA cannot separate these. Every RNA count is therefore a lower bound.
- Silent copies are invisible; copies with identical exons merge. Both make the bound lower, never higher.
- What can push a count ABOVE the truth, and is therefore what this test hunts for: a paralog read as an allele (PSV or mis-assignment),
  sequencing error, RNA editing.

## Substrate

- **Individual:** KB3781 (Jim, mGorGor1), male. The fibroblast Iso-Seq is the same animal (shown 2026-08-13 from SRA run accessions and
  a homozygous-alt drop with a heterozygous internal control).
- **Reference:** `GGO.fasta` = GCF_029281585.2 (`_pri`), a mosaic of whole chromosomes, 16 paternal and 9 maternal.
- **Haplotypes (truth only, never read by the RNA side):** `gorilla_haps/mat.fa` (GCA_028885495.2), `gorilla_haps/pat.fa` (GCA_028885475.2).
- **RNA:** `fibroblasts/GCA_029281585.2_flnc_mm.bam` (4 runs; `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1` to `_pri`).
- **Annotation:** NCBI RefSeq `GGO_genomic.gff` (gorilla keeps RefSeq; CAT applies to human only).

## Gene sets

| set | role | definition |
|---|---|---|
| S_fam | test (named) | the 25 NPIP and 14 TBC1D3 gorilla copies of the 29 Sep copy-recovery truth (`copy_recovery_tools/ann/copies.ggo.tsv`) |
| S_multi | test (genome-wide) | autosomal RefSeq genes whose exon union has another `_pri` locus at >= 0.90 identity over >= 0.50 of it (`minimap2 -x splice -N 50` of the exon-sum sequence against `_pri`) |
| S_single | calibration | autosomal RefSeq genes with no such second locus |
| S_X | negative control | chrX genes outside the pseudo-autosomal region (exon sum maps to chrY at < 0.90 identity), single-copy by the same test. KB3781 is male, so each has one haplotype |

A gene or copy is **expressed** when at least one exonic position has >= 10 of its assigned reads.

## Truth (built from the assemblies before any RNA call; frozen by sha1 in the result)

1. Chromosome correspondence `_pri` <-> mat / pat by exact sequence (each `_pri` chromosome is identical to one haplotype's; the other
   haplotype is its **B** chromosome).
2. Each `_pri` chromosome aligned to its B chromosome: `minimap2 -x asm5 -c --cs`, query sharded (`tools/mm2_shard.sh`), primary
   alignments only (`tp:A:P`).
3. For each gene or copy G on `_pri`, lift its exon union through those alignments to B:
   - **T2d** (both haplotypes, distinguishable): >= 95% of G's exonic bases lift, and B differs from `_pri` at >= 1 lifted exonic base.
   - **T2i** (both haplotypes, identical exons): >= 95% lift, 0 differences.
   - **T1** (one haplotype): < 50% lifts, and no B locus elsewhere is closer to G's exon sequence than G's closest `_pri` paralog is
     (`minimap2 -x splice -N 50` of G's exon sum against B).
   - **T?** (excluded and counted): 50-95% lifts, or a closer B locus exists away from the syntenic position.
   - chrX genes outside the PAR are T1 by construction (male).
4. **B-only copies** of the S_fam families: loci of the family's copy exon sums on B (`-x splice -N 50`, >= 0.90 identity, >= 0.80
   coverage) that are not the lift of any `_pri` copy.
5. **Haplotype copies of a family** T = sum over its `_pri` copies of (2 if T2d or T2i, 1 if T1) + its B-only copies.

## The RNA-only caller (fixed now; reads `_pri`, the BAM and the annotation only)

1. **Reads of G.** Single-copy sets: primary records (`-F 2308`) with MAPQ >= 20 overlapping G's exons. Family sets: primaries that are
   not AS-tied genome-wide (second AS < 0.98 x best, from an `as_table` scan of this BAM) overlapping G's exons, plus AS-tied reads that
   the shipped O2 (`copy_assign`, default settings) assigns to G. Reads O2 abstains on are dropped.
2. **Allele sites** in G's exons (`_pri` coordinates), all required:
   - >= 10 of G's reads cover it; the minor base is a substitution seen in >= 3 reads and >= 0.20 of the covering reads (indels ignored);
   - not a PSV: G's paralog loci (other `_pri` hits of G's exon sum, as above) all carry G's reference base at that column;
   - not editing-like: not A>G on the transcript strand;
   - not within 3 bp of a homopolymer >= 5 bp.
3. **Phasing.** With >= 2 sites, the reads spanning >= 2 sites must split into two haplotypes with >= 80% of them consistent; otherwise
   the gene is "inconsistent" and gets no "2" call.
4. **Call per G:** **2** (>= 1 allele site, phasing consistent), **1+** (expressed, no site), **NA** (not expressed).
5. **Novel haplotypes** (S_fam): family reads whose bases at the family's PSV columns are at Hamming distance >= 2 from every `_pri` copy,
   grouped by identical PSV pattern, >= 2 reads per group, the pattern consistent along each read.
6. **RNA lower bound for a family:** L = sum over expressed copies of (2 if called 2, else 1) + novel groups.

## Hypotheses and decision rules

- **H1 detector noise (S_X).** f_X = expressed S_X genes called 2 / expressed S_X genes. **PASS f_X <= 0.02; FAIL f_X > 0.05** (between:
  marginal). On FAIL, H2 and H3 are reported but not interpreted.
- **H1b calibration (S_single).** Precision of "2" against T2d, and recall among expressed T2d. Descriptive; sets the expected recall.
- **H2 the advisor's claim (S_fam and S_multi, expressed copies with truth T2d, T2i or T1).** precision = T2d among 2-calls;
  false-2 on T1 = T1 copies called 2 / expressed T1 copies.
  - **HOLDS:** precision >= 0.90 and false-2 on T1 <= 0.10.
  - **PARTIAL:** precision >= 0.75 and false-2 on T1 <= 0.25.
  - **FAILS:** otherwise.
  Recall on T2d is reported beside the verdict and is not part of it. S_fam (NPIP, TBC1D3) is reported on its own; with n < 40 it is
  descriptive.
- **H3 lower bound.** Per family (S_fam) and per S_multi cluster with >= 1 expressed copy: **VALID** if L <= T in >= 95% of them.
- **H4 novel haplotypes (S_fam), descriptive.** Each novel group matched (identity at the PSV columns) to a B-only copy, to the B allele
  of a T2d copy, or to nothing.
- **The blind spot, reported as a number:** expressed T1 copies called 1+ against expressed T2i copies called 1+ (RNA cannot tell these
  apart; their sizes say how much of the truth RNA can never reach).

**Overall:** RNA-only allele counting is usable as a lower bound if H1 PASS, H2 HOLDS and H3 VALID; usable with caution if H2 is PARTIAL;
not usable if H2 FAILS or H1 FAILS.

## Execution (in order; heavy steps foreground under `tools/rlock.sh heavy`, < 9 min per call)

1. Truth: chromosome correspondence; sharded `_pri` -> B alignments; lift tables; B-only copies; freeze with sha1.
2. Gene sets: exon sums from `GGO_genomic.gff`; second-locus test against `_pri`; PAR test against chrY.
3. RNA: `as_table` scan of the fibroblast BAM; `copy_assign` (O2 default) on the S_fam and S_multi regions; per-gene pileups; caller.
4. Score H1-H4 against the frozen truth; result in `docs/RNA_ALLELE_HAPLOTYPE_COUNT_2026-10-01.md`. No threshold above is changed after
   step 1; any deviation is written down with its reason.
