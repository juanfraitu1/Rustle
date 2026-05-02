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

Pre-fix:

| family  | StringTie | rustle EM | rustle EM+SNP |
|---|---|---|---|
| GOLGA6  | 2 | 2 | 2 |
| GOLGA8  | 6 | 7 | 7 |
| AMY     | 1 | 0 | 0 |
| TBC1D3  | 4 | 4 | 4 |
| TOTAL   | 13 | 13 | 13 |

Post-fix (commit `c69e7fd`, dropped-family-exempt):

| family  | StringTie | rustle EM-v3 | rustle EM+SNP-v3 |
|---|---|---|---|
| GOLGA6  | 2 | **3**  | **3**  |
| GOLGA8  | 6 | **12** | **12** |
| AMY     | 1 | **2**  | **2**  |
| TBC1D3  | 4 | **8**  | **8**  |
| **TOTAL** | **13** | **25** | **25** |

**rustle gains +12 exact-intron-chain transcripts over StringTie** —
even better than the +10 paralog gain. GOLGA8 nearly doubles (6 → 12).
Note GOLGA6 picks up +1 transcript at intron-chain level even though
paralog count tied at 6 (different alternative isoforms found).

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

3. **`--vg` cross-mapping noise: complete fix landed in `c69e7fd`.**

   Investigation chain:
   - Original behavior: stripped non-family secondaries, didn't rebuild
     `junction_stats` → stale junctions polluted splice graphs.
   - Conservative fix (`5c8e401`): strip + rebuild for non-family
     bundles. Helped chr19-isolated runs (+3 GOLGA8 paralogs) but
     *didn't generalize* to full GGO.bam.
   - Root-cause investigation: traced full-BAM behavior on the
     LOC134757307 bundle. On chr19-only, the bundle was in a 27-copy
     family that passed the filter and kept its secondaries. On
     chr19+chr17, family discovery merged chr19 GOLGA8 with chr17
     TBC1D3 into one 56-copy family (cross-mapping via repetitive
     elements / TEs), `>` max_copies=30 → dropped as megafamily →
     bundle became "non-family" → the conservative fix stripped
     secondaries → low-coverage primaries (1 read at LOC134757307)
     could no longer support path emission.
   - **Final fix (`c69e7fd`)**: track the union of bundle indices
     across ALL *raw* (pre-filter) families. Any bundle that was in a
     raw family gets its secondaries preserved at strip time, even if
     its family was filter-dropped. Stripping is reserved for bundles
     that never participated in family discovery at all.

   **Full GGO.bam per-paralog recovery (the headline result):**

   | family   | ref | StringTie | rustle EM-v3 | rustle EM+SNP-v3 |
   |----------|-----|-----------|--------------|------------------|
   | GOLGA6   | 16  | 6         | 6            | 6                |
   | GOLGA8   | 17  | 8         | **13**       | **13**           |
   | AMY      |  3  | 2         | **3**        | **3**            |
   | TBC1D3   | 13  | 7         | **11**       | **11**           |
   | **TOTAL**| 49  | **23**    | **33**       | **33**           |

   **rustle gains +10 paralogs over StringTie** across the four
   multi-copy proof families. New unique recoveries (vs ST):
   - GOLGA8 +5: LOC129527057, LOC101151050, LOC134757307,
     LOC129527021, LOC115930779
   - AMY +1: LOC101133335 (alpha-amylase 1B)
   - TBC1D3 +4: LOC101151653, LOC129533808, LOC129533806, LOC129533813

   GOLGA6L10 (the original "missing paralog" mystery) is still missed
   by default — its bundle was already in a kept family (no megafamily
   drop), so the secondary-noise issue there is different in nature.
   Available as opt-in via `RUSTLE_VG_STRIP_FAMILY_SECONDARY=1`
   (recovers GOLGA6L10 but trades LOC129523543; net 0).

4. **rustle missed 1 AMY transcript ST caught.**
   At the locus level both find AMY2A and AMY2B (2/3); ST has an exact
   intron-chain match for one of them, rustle's transcripts have a small
   exon-boundary discrepancy. Likely a parity issue, not a multi-copy
   one — same paralogs, just slight transcript-structure mismatch.

5. **15/49 paralogs still missed by all methods (post-fix)** — broken down
   by failure mode (BAM coverage probed via samtools at the ref ranges):

   - **Truly silent (≤3 total alignments, 6 paralogs)**: cannot be
     recovered without more depth or external annotation guides:
     LOC134757233 (0/0), LOC134759231 (0/0), LOC134757232 (0/1),
     LOC101153918 (1/3), LOC101153208 (1/2), LOC109026840 (1/1).
   - **Low-cov assembly gap (4-7 primaries, 7 paralogs)**: reads
     present but bundle thresholds reject them:
     LOC101138059 (4/7), LOC115933515 (4/8), LOC101133389 (4/8),
     LOC115930772 (5/7), LOC115930818 (2/5), LOC101140620 (7/7),
     LOC134757625 (6/67, where 67 includes secondaries),
     LOC101137218 (5/66). These are the depth-aware-bundling targets.
   - **Silent paralog with antisense expression (1 paralog)**: GOLGA6L7
     (chr19 29848-29854, `-` strand) has 52 primary reads at its range
     but ALL on `+` strand (alignment AND `ts:A:+` tag). They belong to
     `LOC115930831` — an overlapping `+`-strand lncRNA at 29838-29861.
     Both rustle and ST correctly assemble the lncRNA on `+` and
     correctly emit nothing for GOLGA6L7 (which has 0 `-`-strand reads
     in this BAM). gffcompare's locus-overlap matching can't tell the
     two genes apart, so it appears as a "miss". Not a tool bug.

   The +10 paralog gain from rustle's multi-mapper EM (vs ST) directly
   addresses the low-cov-assembly-gap class — by preserving cross-
   mapping evidence even when family discovery misclassifies a cluster.

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
