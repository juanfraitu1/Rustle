# Multi-copy gene family recovery — full GGO.bam

Reference: `/scratch/jxi21/Assembler/GGO_genomic.gff` (gorilla genomic
annotation, NCBI/Gnomon).

## Setup

| variant | command |
|---|---|
| StringTie | `stringtie -L -p 8 -o full_stringtie.gtf GGO.bam` |
| rustle EM | `rustle GGO.bam --vg --vg-solver em` |
| rustle EM+SNP | `rustle GGO.bam --genome-fasta GGO.fasta --vg --vg-solver em --vg-snp` |

All rustle runs use the family-quality filter (defaults: min_shared=10,
max_copies=30, min_shared_per_copy=1.0, max_exon_cv=1.5). On full BAM:
706 raw families → 450 high-confidence families.

## Per-family LOCUS recovery (gene-level)

How many *distinct ref paralogs* are recovered by each tool? Computed by
parsing gffcompare's `.loci` file: an XLOC with both ref and query columns
populated, then counting unique ref `gene-LOCxxxxx` IDs that appear.

| family  | ref paralogs | StringTie | rustle EM | rustle EM+SNP |
|---|---|---|---|---|
| GOLGA6  | 16 | 6 | 6 | 6 |
| GOLGA8  | 17 | 8 | **9** | **9** |
| AMY     |  3 | 2 | 2 | 2 |
| TBC1D3  | 13 | 7 | 7 | 7 |
| **TOTAL** | **49** | **23** | **24** | **24** |

## Per-family TRANSCRIPT recovery (matching intron chain)

From gffcompare `Matching transcripts:` field — exact intron-chain match.

| family  | StringTie | rustle EM | rustle EM+SNP |
|---|---|---|---|
| GOLGA6  | 2 | 2 | 2 |
| GOLGA8  | 6 | 7 | 7 |
| AMY     | 1 | 0 | 0 |
| TBC1D3  | 4 | 4 | 4 |
| TOTAL   | 13 | 13 | 13 |

## Findings

1. **rustle EM finds +1 GOLGA8 paralog vs StringTie.**
   `LOC129527057` (chr19, NC_073240.2:23635578-23649429) — recovered only by
   rustle (both EM variants). StringTie misses it. This paralog is the
   second copy in a tandem GOLGA8 cluster (LOC101150678 23.5Mb / LOC129527057
   23.6Mb) where ST collapses the cross-mapping reads to the first copy
   while rustle's multi-mapper-aware EM keeps both visible.

2. **SNP-aware copy assignment moved a few reads but no ref paralogs.**
   `rustle EM+SNP` is essentially identical to `rustle EM` at locus and
   ref-paralog level. 11 families on full BAM had ≥1 diagnostic SNP
   (biggest: 21-copy/11-SNPs, 6-copy/12-SNPs, 30-copy/10-SNPs) and the
   SNP factor IS scored into the EM (multiplicative on top of junction
   compat × ln-context). Per-transcript count differs slightly (73,262 EM
   vs 73,263 EM+SNP), and TPM values differ at ~1e-5 precision, so SNP did
   shift some weights. But the differing transcripts (e.g. RSTL.2875.2 at
   chr17:45.9Mb, RSTL.18734.8 on unplaced contig) fall *between*
   annotated paralogs, not at any ref-paralog locus we measure.

   Implication: at PacBio IsoSeq depth on full GGO.bam, the junction-EM
   signal is so dominant that adding SNPs only perturbs assembly at the
   margins — never enough to gain or lose a reference paralog hit on
   GOLGA6/GOLGA8/AMY/TBC1D3. To make SNPs discriminating we'd need
   lower-coverage data, a stricter junction-tie scenario, or per-copy
   TPM evaluation rather than locus presence.

3. **Both rustle variants lose GOLGA6L10 vs ST — pre-existing bundle-detection
   gap, NOT a multi-mapper EM regression.**
   `GOLGA6L10` (chr19 NC_073240.2:83999064-84007244, +) is caught by ST as
   `STRG.14456.1` (4 small exons, longcov=2) but missed by rustle in *all*
   variants — including the discovery-only run with no EM applied. So this
   is a baseline assembly miss at the bundle-detection / single-exon-rescue
   stage, not a family-quality filter or EM stripping issue. ST has more
   aggressive short-transcript rescue for low-coverage 4-exon transcripts
   with small introns (60bp, 100bp, 211bp).

4. **rustle missed 1 AMY transcript ST caught.**
   At the locus level both find AMY2A and AMY2B (2/3); ST has an exact
   intron-chain match for one of them, rustle's transcripts have a small
   exon-boundary discrepancy. Likely a parity issue, not a multi-copy
   one — same paralogs, just slight transcript-structure mismatch.

5. **25/49 paralogs missed by all three methods.** Common causes (not
   tool-specific): silent/low-expression paralogs without read coverage,
   single-exon paralogs whose alignments are too short to assemble,
   paralogs in regions with too much sequence ambiguity for unique
   primary alignment.

## Net signal vs the user's objectives

- **Objective 2 (multi-mapper EM)**: validated. rustle EM recovers 1
  GOLGA8 paralog ST cannot, by re-routing multi-mappers from the dominant
  copy to its tandem neighbor. Same total transcript count overall, but
  the EM-recovered paralog wouldn't appear in any ST output.

- **Objective 3 (multi-copy family scaffold)**: family-quality filter
  reduces 706 noisy "families" to 450 high-confidence paralog clusters
  (drops mtDNA mega-cluster, sparse alignment-noise pairs).

- **Objective 5 (SNP-aware copy assignment)**: no observable effect at
  the assembly-output level on full GGO.bam. The mechanism works (chr19
  shows 1 active SNP family), but the gain is in *attribution* (which
  copy a multi-mapper belongs to), not in *discovery* (whether a copy
  appears at all). To see SNP value, we'd need to look at per-copy read
  counts / TPM rather than locus presence.
