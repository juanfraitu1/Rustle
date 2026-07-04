# Do non-expressed genes still carry lots of multimapping reads? — Only if they are multi-copy (or repeat-overlapping)

**The disagreement.** *Advisor:* genes that are not expressed still have lots of multimapping reads
covering them. *Student:* not so, **if** the non-expressed gene is **unique** (single-copy, not part of a
multi-copy family) — a unique non-expressed gene has ~nothing mapping to it.

**Verdict: the student is right for the clean case; the advisor is right *specifically* for multi-copy
genes; the repeat confound is real but secondary.**

## The logic

A non-expressed gene transcribes ~nothing, so any read *covering* it is **spillover**, and spillover has
exactly two sources: **paralog** spillover (reads from an expressed copy of a multi-copy family multimap
onto a silent copy) and **repeat** spillover (reads from other TE-containing loci multimap onto a gene
overlapping a transposable element). A gene that is **unique AND repeat-free has neither source → ~0**.

## The decisive 2×2 (gorilla `GGO_mm.bam`, 34,114 RefSeq genes)

non-expressed = ≤1 unique primary read (MAPQ>0); repeat-free = soft-mask fraction <0.10; multimap =
secondary (0x100) + primary-MAPQ0 + supplementary (0x800), splice-aware block overlap.

| cell | n | multimap median | mean | **frac ≥1** | frac ≥5 |
|---|---|---|---|---|---|
| **UNIQUE · repeat-free** (the student's clean case) | 2,082 | **0** | 10.5 | **0.097** | 0.057 |
| UNIQUE · repeat-overlap | 5,421 | 0 | 28.1 | 0.286 | 0.152 |
| **MULTI-COPY · repeat-free** | 1,407 | **4** | 139 | **0.652** | 0.474 |
| MULTI-COPY · repeat-overlap | 2,037 | 11 | 137 | 0.836 | 0.649 |

**Each axis isolated (clean dissociation):**
- **Copy-status, holding repeat-free: 0.097 → 0.652 = 6.7×.** Paralog spillover dominates, and it's
  genuine — it holds *even when repeat-free*, so it is not a repeat artifact. This is the advisor's case.
- **Repeat, holding unique: 0.097 → 0.286 = 2.9×.** Real but secondary — and the **median stays 0**. A
  repeat-overlapping unique gene gets *some* spillover but is never *broadly* multimapped. This is almost
  certainly the source of the advisor's intuition (the genome is TE-dense), but it's ~half the copy-status
  effect and never covers a gene at the median.
- Robust across non-expressed thresholds 0/1/2 (unique-repeat-free frac≥1 = 0.093/0.097/0.100).

So a non-expressed gene that is genuinely **unique and repeat-free attracts ~nothing** (median 0, 90% get
zero). **The student is right.**

## The residual in the clean cell is uncatalogued paralogs — not counterexamples

In UNIQUE·repeat-free the top 1% carry 78% of the cell's multimap mass, and those outliers are
**uncatalogued paralogs**: `LOC109027830` (ferritin-heavy-like), ferritin-light-like, SEC22b-like,
UFM1-like — retrocopy / pseudo-paralog "-like" genes that ARE real family members the catalog missed
(single-copy by name, <50% segdup). The multimapping signal **correctly flagged family members the
definition didn't catch.** Removing them makes the student's claim cleaner *and* proves the deeper point:

⭐ **Multimapping coverage of a silent locus is itself a copy-family signal — a paralog or a repeat, never
noise.** This is why the definition can treat "unique non-expressed gene = ~0 reads" as clean, and why a
silent locus lit up by spillover is evidence of hidden paralogy (which is how the ferritin-like genes
announced themselves).

## Corroborating detail (for the advisor conversation)

**49% of multi-copy genes (3,444 / 7,021) score "non-expressed" by unique-primary reads while being
saturated by multimapping** — reads are ambiguous across paralogs and land as secondary. e.g. **ACTA1**
(skeletal actin, in segdup+family): ~1 unique primary read but **19,992 secondary spillover reads**. So
"looks silent, heavily multimapped" is real and common — *for multi-copy genes* — which is likely the
pattern the advisor generalized to all non-expressed genes.

## Definitions / caveats

- minimap2 `-N50 -p0.1` expresses ambiguity as **secondary** alignments, not MAPQ0 (BAM totals: secondary
  4.04 M vs primary-MAPQ0 only 6 k) — so multimap is captured via the 0x100 flag.
- copy-status multi = gene in `genome_families_refined.tsv` ∪ `family_rna_refine.tsv` ∪ SEDEF `final.bed`
  ≥50%-covered → 7,021 multi / 27,093 unique.
- repeat_frac = soft-mask fraction of the **gene span** (genomic, intron-inclusive); RNA spillover is
  exonic, so this overstates the mechanistic repeat load for large intron-heavy genes (e.g. GLRA1).

Scripts: `nonexpr_build_gene_annot.py` → `nonexpr_measure_bam.py` → `nonexpr_build_final_table.py` →
`nonexpr_write_summary.py`. Per-gene table (34,114 rows) + summary in `winloci_scratch/`
(`nonexpr_multimap_pergene.tsv`). Deterministic (`PYTHONHASHSEED=0`); adversarially verified (spillover
measure genuine, repeat confound isolated, hidden-paralog pollution of the unique cell audited).
