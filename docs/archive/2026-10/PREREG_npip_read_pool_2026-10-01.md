# Pre-registration: which alignments build the loci, shown on NPIP (human chr16, A119b), 2026-10-01

Written before any of the three arms below was run on chr16. Follows `docs/PREREG_locus_read_pool_2026-09-22.md` (§6z7, register
1059-1061: chr20 + gorilla NC_073244.2), which compared the same three read pools genome-scale but never on NPIP.

## Question (advisor)

Primaries-only loci are unsafe because the primary of a tied multi-mapper is a coin toss. The advisor proposes building loci from ALL
alignments (primaries + every secondary) and expects the later all-vs-all minimap2 step to remove the false positives. On NPIP, with
real data: what does each read pool do to locus definition, and does the all-vs-all + clustering step remove what ALL adds?

## Arms (one frozen binary, one command, only the read pool differs)

`/mnt/linuxdisk/tmp/rustle_figures/cc_bin_frozen/copy_assign --assemble-only --regions chr16 --assembly-junctions strict` + the
shipped polish (`--assembly-polish full --polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82
--polish-ism-ratio 0.7 --polish-retained-ratio 10`), BAM `A119b.t2t.bam` (`-N 50 -p 0.1` secondaries present).

| arm | read pool | env |
|---|---|---|
| P | primaries only | none |
| GOOD | primaries + secondaries with AS >= 0.98 x the read's genome-wide best AS (the pipeline default since 2026-09-24) | `RUSTLE_GTF_SECONDARY=1 RUSTLE_GTF_SECONDARY_AS_RATIO=0.98 RUSTLE_GTF_SECONDARY_AS_TABLE=human_A119b.molecules.tsv` |
| ALL | primaries + every secondary | `RUSTLE_GTF_SECONDARY=1` |

Downstream, identical for all three: one locus per `gene_id` (`gw22/sec/loci_from_gtf.py`: span = hull of its transcripts, rep = most
reads), `minimap2 -x asm20 -c -X -N 50 -p 0.1 --secondary=yes` all-vs-all of the locus spans, `mcl_families --min-exonic-bp 1
--min-shared-exon-frac 0.60` with the binary's current defaults (`--min-cov-shorter 0.70`, `--bridge-regroup f1v2`).

## Definitions (fixed now)

- **NPIP copies:** the 25 chr16 copies of the CAT/Liftoff v2.0 truth (`copy_recovery_tools_cat/sens_a4_npipp1/ann/copies.hsa.tsv`,
  protocol Amendment 2), each as its CAT gene's exon union. RefSeq NPIPB3 (no CAT gene) is reported separately by its RefSeq span.
- **Locus on a copy:** the locus's rep exons overlap the copy's exon union on the same strand by >= 1 bp.
- **Fused locus:** on >= 2 copies. **Copy covered:** >= 1 locus on it.
- **Echo locus** (the §6z7 definition): zero primary records (`-F 2308`) over the locus span.
- **Annotation class of a locus:** on an NPIP copy / on another CAT gene (rep exons overlap its exons, same strand) / no CAT gene.
- **NPIP family (per arm):** the `mcl_families` cluster holding the most loci that are on NPIP copies.
- **Added by ALL:** an ALL locus whose rep exons overlap no GOOD locus's rep exons on the same strand.

## The advisor's claim, and the decision rule

Claim: the all-vs-all + clustering step removes the false positives ALL adds. The candidate false positives are ALL-added loci that are
echoes or are not on an NPIP copy. Let f = (such loci that end up in the NPIP family) / (such loci).

- f <= 0.10: **the claim holds** (the family step removes them).
- 0.10 < f < 0.50: **partial**.
- f >= 0.50: **the claim fails** (the family step keeps them; homology cannot tell an echo from a copy).

Reported beside f, for all three arms: loci on chr16; loci on NPIP copies; copies covered (of 25); fused loci; echo loci; NPIP-family
size and its composition by annotation class and echo status; all-vs-all records and run time. Descriptive only beyond f.
