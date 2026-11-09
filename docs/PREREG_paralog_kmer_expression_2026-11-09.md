# PREREG row 1596: paralog-specific k-mer expression, Soto's parCN rule moved from DNA reads to RNA reads

Written 2026-11-09, before any k-mer is extracted or counted. Advisor context: Soto 2025 measures paralog copy number
(parCN, QuicK-mer2) by counting, in DNA reads, k-mers that occur once in the reference (paralog-specific). Row 1594 showed
that counting distinct expressed haplotypes does not recover Soto's DNA copy number (FAIL). This row asks a different
question: can the same k-mer rule, applied to RNA reads, give each copy's expression without relying on where the aligner
placed a read? It does **not** claim to measure DNA copy number.

## Inputs (A119b testis Iso-Seq on CHM13 v2.0)

- Families: `/anvme/workspace/v138da13-rustle/port1593/Ab.fam.copies.tsv` (the shipped default, row 1595's families):
  1,395 families, 6,495 copies, each with its exons.
- Genome: `$HOME/data/chm13/chm13v2.0.fa` (uncompressed, as row 1594), **chrY kept**: A119b is male testis.
- Reads: `~/work/jfim/rustle_data/reads/A119b.t2t.bam` (`A119b.flnc.bam` is not on this host).

## Definitions (fixed now; `bench/soto/rnakmer.py`)

- **Copy sequence**: the copy's exons joined in order (spliced), forward strand.
- **k-mers**: k = 30, canonical, case-insensitive, non-ACGT resets, the all-A k-mer skipped: `parcn_assembly.py`'s
  `canonical_kmers` (Soto's QuicK-mer2 settings).
- **SPEC(copy)**: the copy's distinct k-mers that occur in no other copy of any family and occur at most once in the
  whole genome (`parcn_assembly.stream_kmer_counts`, every contig incl. chrY). Exon-internal k-mers then occur exactly
  once in the genome. Splice-junction k-mers occur 0 or 1 times. Deviation from Soto: no edit-depth filter (it drops
  0.44% of once-only k-mers in `parcn_assembly`).
- **Resolved copy**: |SPEC| >= 100 (Soto's RESOLVED threshold).
- **Reads in scope**: every read with a non-supplementary BAM record (primary or secondary) overlapping a copy span
  +- 1 kb. **Copy span** = `copies.tsv` `locus_start`..`locus_end` (the node's genomic extent), here and below. One sequence per read (its records carry SEQ; the primary is taken when present). Out of scope: reads whose
  every record lies off the families. That is r1510's off-catalog class.
- **k-mer assignment**: per read and per family, the number of distinct SPEC k-mers of each copy it contains. The read is
  assigned to copy c iff c has >= 3 hits and every other copy of that family has 0. A read can be assigned in more than
  one family; each family is scored on its own. Read-level assignment is used, not Soto's median over k-mers: RNA depth
  along a copy varies with isoform and 3' bias, so a per-k-mer median does not estimate expression.
- **Unique-placement reads**: primary record MAPQ 60 overlapping exactly one copy's span of that family. Their copy is
  that copy.

## Bars

- **STOP, no verdict**, if fewer than 200 copies are resolved, or fewer than 1,000 unique-placement reads on resolved
  copies are k-mer assigned.
- **PASS** iff all three hold on resolved copies:
  1. **Agreement**: of unique-placement reads that are k-mer assigned, >= 0.98 go to their primary copy.
  2. **Reach**: of unique-placement reads on resolved copies, >= 0.70 are k-mer assigned.
  3. **Counts**: Spearman >= 0.90 between each copy's k-mer assigned reads (all reads in scope) and its unique-placement
     reads, over resolved copies with >= 10 unique-placement reads.

**This is necessary, not sufficient.** SPEC k-mers come from the reference, so they will agree with the aligner where
the aligner already resolves the read. The bars show the counter is correct and reaches most reads. They do not show it
is right where the aligner ties. A119b has no truth for tied reads.

## Reported, no bar

- **Read buckets**: every in-scope read falls in exactly one bucket: unique-placement; tied (primary on a family copy,
  MAPQ < 60, or MAPQ 60 overlapping two copies' spans); primary off the families (fetched through a secondary record);
  unmapped primary. Each bucket is reported with its count and its k-mer assigned share.

- **Tied reads** (primary on a family copy with MAPQ < 60): share that are k-mer assigned, and share of those whose
  k-mer copy is the primary's copy. This is a comparison, not a truth: a disagreement can be the aligner's error or the
  counter's. O2 `--tie-ms` is not run on these families for this row.
- **Effective expressed copies** per family: N_eff = 1 / sum(s_c^2), where s_c = copy c's share of the family's k-mer
  assigned reads. Reported with the family's copy count and the share of its reads on unresolved copies. Families with
  < 20 assigned reads are listed but get no N_eff.
- For families holding an S1E gene (row 1594's mapping by Gene ID through CAT): N_eff beside the CHM13 copy count and
  the S1E median parCN, descriptive only. Row 1594 already showed RNA does not give DNA copy number.
- Per copy (`r1596/per_copy.tsv`): |SPEC|, k-mer assigned reads, unique-placement reads, tied reads assigned.

## What follows

- PASS: the counter is a valid per-copy expression measure where SPEC exists. Next row is the real test, on tied reads
  with truth: KB3781 against Jim's haplotype truth (`haplotype_truth`, the O2 ms-tie truth), k-mer counter vs O2
  `--tie-ms`. N_eff becomes the RNA-side family statistic (expressed copies), reported beside Soto's DNA copy number,
  never as it.
- FAIL: register which bar failed. Agreement < 0.98 points at SPEC k-mers shared with off-family sequence (the genome
  count misses read errors) and is the first thing to fix. Reach < 0.70 points at too few SPEC k-mers per read.
- STOP: register the resolved-copy count; the family copies are too similar for k = 30 at this depth.

Outputs: `/anvme/workspace/v138da13-rustle/r1596/`.

## Outcome (2026-11-09, job 4453985, 1,240 s): FAIL

SPEC k-mers 2,671,098 (of 4,328,059 single-copy k-mers, genome count <= 1); resolved copies **2,571** of 6,495 (STOP bar
200 cleared). 1,944,969 reads in scope.

| bar | value | registered | |
|---|---|---|---|
| 1 agreement | 373,869 / 379,623 = **0.9848** | >= 0.98 | pass |
| 2 reach | 379,623 / 892,507 = **0.4253** | >= 0.70 | **fail** |
| 3 Spearman (1,995 copies) | **0.6715** | >= 0.90 | **fail** |

**FAIL** as registered (bars 2 and 3).

Read buckets, counted over (read, family) pairs; a read is counted once per family it was fetched for or hits:

| bucket | pairs | k-mer assigned |
|---|---|---|
| unique | 970,763 | 390,907 |
| tied | 1,497,433 | 129,402 |
| off | 15,514,193 | 177,001 |

Tied reads: 129,402 / 1,497,433 (8.6%) are k-mer assigned. Of those, 98,222 (76%) go to a copy the primary overlaps.

Reading:
- Where the counter calls a copy, it agrees with the aligner (98.5%).
- It calls fewer than half of the uniquely placed reads on resolved copies.
- Its per-copy counts track the aligner's only loosely (0.67). Counts from all buckets enter bar 3, so off-family reads
  with hits shift the counts.
- Not diagnosed here: whether unassigned reads carry < 3 SPEC hits (SPEC k-mers clustered in a few exons) or hits on
  another copy (A119b alleles matching another copy's reference k-mer).

Per the registered follow-up, agreement held, so there is no sign of off-family SPEC contamination. Low reach means too
few SPEC k-mers per read. The KB3781 tied-read test is not run on a FAIL. N_eff is in `r1596/rk.per_family.tsv`
(descriptive; it rests on a counter that failed bars 2 and 3). Per copy: `r1596/rk.per_copy.tsv`.
