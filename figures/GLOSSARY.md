# Figure glossary: exact definitions and legend wording

This file is the reference for every technical term used in the figure captions (`captions/*.md`), the `META`
claims and on-figure text of `fig_*.py`, and the table notes of `data/*.tsv`. Each entry gives:

- **Code.** The operational definition as the code computes it, with every threshold and flag, and the
  `file:function` it was checked against (repository commit 628894f9, `figures/` as of 2026-09-25).
- **Legend.** Wording for a figure legend or caption. It uses no project jargon, no internal labels and no
  register numbers, and it says "truth" only for the simulation, where the source of each read is known.
- **Occurs.** Where the term appears, as `file:line`. The lists were produced by grep over `captions/`,
  `fig_*.py`, `_o1.py`, `_o1_recovery.py`, `_o2.py`, `_sqanti.py`, `figlib.py`, `README.md` and the `# note`
  lines of `data/*.tsv`. A trailing `…` means the list was cut at 25 lines per file.

Part A lists the words that mean different things in different figures, and the factual mismatches found while
checking the code. Part B is the glossary. Part C is the style guide.

> **2026-09-25 (user decision).** Rustle has ONE default de novo family definition: reads → seeded assembly loci →
> one representative per locus (its *positional exon sum*) → families (the driver's `families` stage). Copy
> assignment consumes the same families (their copy table). The `gw_family_catalog` copy catalog is **legacy**, and
> protein is an optional, manually invoked extra-sensitive step shown only in a supplementary figure. Figures 6–8 and
> their captions use the default families and external references (Ensembl Compara, Soto 2025, Liftoff copies); fig.
> 6d moved to the supplementary fig. 6s-seeding and the translated protein tier to fig. 6s-protein. Entries below
> that describe the catalog, the protein-homology families or old panel letters say so. See Part B §4, *Default de
> novo families*.

---

## Part A. Problems to fix before the figures go to the advisor

### A1. One word, several meanings

| word | meanings in the current text | proposal |
|---|---|---|
| **tied / tie** | (1) Figs 1–3 and 6d ("tied secondaries", "near-tied", "tie fraction"): a secondary alignment, or a read's second-best alignment, scores **≥ 0.98 ×** the read's genome-wide best alignment score (AS). (2) Figs 4–5a,b ("tied reads", "aligner tie"): the simulated read's **primary alignment has MAPQ 0**. It is not always an AS tie: 3 of the 1,263 human MAPQ-0 reads have the best AS on the true copy alone (`fig4_assignability_upset.tsv:31`). (3) Fig. 5c ("AS-tied", hard set): `copy_assign`'s gate, **runner-up AS = best AS exactly** (ratio 1.0), counted only inside the swept chr16 window. (4) Figs 4–5, "stay tied": the `copy_assign` verdict `tied`, meaning the read cannot reach significance against some candidate copy. (5) Statistical equality: "the two arms are tied in human" (fig1.md:10), "in this range they are a tie" (fig6.md:7), "a tie on 3 families" (fig7.md:118, 144, 257). | Keep "tie" for alignments only, and always say which kind: "MAPQ 0", "equal best alignment score" or "second alignment ≥ 98% of the best score". For (4) write "left unassigned". For (5) write "equal" or "not separable". |
| **expressed** | (1) Fig. 3: a multi-exon RefSeq transcript with **≥ 1 read** that has a candidate placement (primary, or secondary ≥ 0.98 × best AS) whose aligned block overlaps one of its exons (fig3.md:50–51, 174). (2) Fig. 6d recall universe: a gene with **≥ 2 reads whose primary alignment** has an aligned block on its annotated exons. Files and flags call it `expressed` (`--expressed`, `NC_073244.2.expressed.exon.tsv`, fig6.md:195). (3) Recorded span rule, used by fig. 7's caveat "the referee's expressed same-family pairs" (fig7.md:221): **≥ 2 primary records overlapping the gene span**, spliced-over reads included. (4) Fig. 6a–c universe: a gene has **a locus in Rustle's chr16 assembly** with a representative of ≥ 200 bp. `score.py spectrum` prints it as "both expressed", and fig6.md:234 ("no expressed Compara pair") uses that sense. (5) Informal: fig1.md:96 "unexpressed transcripts" (no rule given), fig7.md:214 "expressed or not", fig3.md:207 "which copy is expressed" (biological). | Never use "expressed" without its rule. Use "with ≥ 1 read (candidate placement on an exon)", "with ≥ 2 primary reads on its exons", "with a locus in the assembly". For fig. 1, "transcripts that no read carries exactly". |
| **held out** | (1) Fig. 7 status: "held out, untouched" (no rule, threshold or verdict ever scored on it) and "held out, reused verdict set" (verdict set of about 30 earlier guided-mode tests). (2) Figs 4–5: "gorilla (held out)" for **NC_073244.2** (fig4.md:18, 120; fig_assignability.py:42). The same contig is **development** in figs 6 and 7, because the seeding default was decided on it. Both statements are true, but for different decisions: the copy-assignment rules were developed on human chr16, and the seeding rule on NC_073244.2. (3) Fig. 6: "held-out substrate of the seeding pre-registration" (fig6.md:32). | Always name the decision: "not used to develop the copy-assignment rules" or "used to choose the seeding rule". Do not call NC_073244.2 "held out" anywhere without that qualifier. |
| **development** | (1) Substrate status: rules were chosen on it (fig. 7; fig. 6 "development chromosome"). (2) An earlier recorded run: "hash-seeded development simulation" (fig4.md:17; fig5.md:17), "recorded development runs" (README.md:14; the `PROVISIONAL` stamp, figlib.py:295). | Use "development" only for (1). For (2), write "earlier run (random seed not fixed)". |
| **seed / seeded** | (1) Admitting secondary alignments into the read pool the assembler builds loci from (figs 1–3, 6d, 7). (2) A random-number seed: "hash-seeded", "stable-seed", "seed 20260925" (fig4.md:17–18; fig5.md:38). | Keep "seeding" for (1). For (2), write "random seed 20260925" and "earlier run without a fixed random seed". |
| **family** | (1) **Default family** (since 2026-09-25): a family of the driver's `families` stage (`mcl_families --from-gtf`), the one default de novo definition (figs 6–8; its copy table feeds figs 4–5). (2) **Catalog family** (LEGACY): a group of copies in `gw_family_catalog`'s output (the figs 4–6 tables built before 2026-09-25). (3) **Reference families**: Compara families (fig. 7), Soto 2025 families, the NPIP union set, Liftoff copy pairs (a relation, not a partition), protein-homology families (secondary, supplementary only), and "paralogue groups" (fig. 6), which are connected components of Compara pairs. (4) Fig. 5c's "NPIP family, 26 copies" (`copies16.tsv`). | Say "Rustle's default family" (or "predicted family") for (1); "legacy catalog family" for (2), only where it is compared. Say "reference family" for (3), always with its source. |
| **copy** | (1) **Copy** of a default family: one member locus, represented by its most-read transcript's exons (the copy table `<id>.fam.copies.tsv`; figs 4–8); formerly a **catalog copy** (legacy). (2) **Source copy**: the copy a simulated read was drawn from. (3) The 26 NPIP copies of the hard-locus benchmark (fig. 5c). (4) The biological paralogue at a genomic position (fig3.md:190, 207). | Define "copy" once per figure. In figs 4–5 write "source copy" for (2). |
| **catalog / catalogue** | (1) `gw_family_catalog` output: LEGACY since 2026-09-25 (the default families replace it in figs 6–8). (2) Fig. 7 "Catalogue sizes" and "share a catalogue" mean the **predicted clusters of `mcl_families`**, and use British spelling. | "legacy copy catalog" only for (1). For (2), "predicted families". Spell it "catalog" throughout. |
| **locus** | (1) Fig. 6a–c: one `gene_id` group of Rustle's assembled GTF, whose representative is its transcript with the most reads, kept if ≥ 200 bp. (2) Fig. 7 de novo: the same `gene_id` group, taken over its genomic span. (3) Fig. 7 guided: an annotated gene or pseudogene body; bodies whose exon unions overlap are folded into one locus. (4) "Same locus" for scoring (figs 4–5): two catalog spans on one chromosome that overlap by ≥ 50% of the shorter one (`score.py:same_locus`). (5) A genomic position in general. | Define per figure. |
| **cluster** | (1) Bootstrap unit: "cluster interval", "cluster bootstrap" (figs 3, 5), where a cluster is a read-sharing group or gene (fig. 3) or a source copy (fig. 5). (2) A predicted family of `mcl_families` (figs 6d, 7). | For (1), write "95% interval, resampling read-sharing groups" (or "…source copies"). Keep "cluster" for (2). |
| **group** | (1) Fig. 3 "read group": a **read-sharing group**. This collides with the SAM/BAM read group (`@RG`). (2) Fig. 6 "Compara group": a connected component of the universe's Compara pairs, built by us and not an Ensembl object. (3) "gene group" (fig6.md:9). | Always write "read-sharing group" (the panels print "read groups": fig_secondary.py:1029–1030, 1109, 1133). Write "paralogue group (connected component of the Compara pairs)". |
| **identity** | (1) Figs 4–5 bands: 1 − minimap2 `de` (gap-compressed divergence) between a copy and its most similar **directly aligned** copy in its catalog family (`max_family_identity`). (2) Fig. 6a–c x-axis: Ensembl Compara **protein** percent identity, the larger of the two directional values. (3) Fig. 6c edge rows and fig. 6a direct edge: PAF matching bases ÷ alignment block length (nucleotide tiers), or mmseqs `pident` (protein tier). (4) Fig. 6d: best minimap2 (`asm20 -k11 -w5`) matches ÷ block length between the two genes' **annotated mRNAs**. (5) Protein-homology families: **no identity used**. | State the kind in every axis label: "nucleotide identity (…)" or "protein identity (Compara)". |
| **coverage** | (1) Fig. 5 "coverage": assigned tied reads ÷ tied reads, which is an assignment rate. (2) Alignment coverage: ≥ 50% of the shorter copy (catalog edge, "aligned core"), ≥ 30% of the longer protein (protein-homology families), ≥ 30% of the longer locus plus ≥ 60% shared-exon fraction (family rule). (3) Short-read coverage (fig2.md:118). | For (1), write "fraction assigned". Keep "coverage" for alignments. |
| **accuracy** | (1) Fig. 5: correct ÷ assigned, which is a precision among assigned reads. (2) Fig. 1 title "Intron-chain accuracy", meaning sensitivity and precision together. | Fig. 5: "fraction correct among assigned (precision)". Fig. 1: "Intron-chain sensitivity and precision". |
| **sensitivity / recall** | Fig. 6 says "recall"; figs 1, 3 and 7 say "sensitivity"; the fig. 7 table notes and module docstring say "recall = matched / truth members" (fig7_summary.tsv:54, fig_family_recovery.py:13). | Use "sensitivity" everywhere (the user's reporting rule is sensitivity, precision and bipartite matching). |
| **universe** | (1) Fig. 6a–c: 290 Compara pairs whose two genes both have an assembled locus. (2) Fig. 6d: 576 genes (recall universe). (3) Fig. 7: the genes of the reference families with ≥ 2 members on the chromosome; predicted clusters are intersected with it. (4) Fig. 5c: molecules whose primary alignment lies inside a copy. | Replace with "the pairs scored", "the genes scored" and so on, each with its rule. |
| **mode** | (1) Transcript modes (figs 1–3): annotation-free (de novo) vs annotation-guided tool runs. (2) Family modes (figs 6–7): de novo (loci assembled from reads) vs guided (loci = annotated gene and pseudogene bodies). (3) Locus modes (fig. 8): Rustle's guided candidate search vs its de novo loci, in the Liftoff framework. (4) Count mode `reads` / `distinct_ends` (fig. 1). (5) "seeding modes" (README.md; now "seeding configurations" there). | Reserve "mode" for (1)–(3) and always qualify it: "transcript mode", "family mode", "locus mode" (Part B, *Modes*). (1) and (2) are different axes: every transcript figure is annotation-free on every side, while fig. 7's guided family mode reads no RNA at all. Write "count rule" for (4) and "configuration" for (5). |
| **certificate** | (1) PSV certificate, the per-family assignment test. (2) Union certificate. (3) Origin certificate (inside `copy_assign`; it shows up as "origin-rejected"). (4) "within-family distance certificate" = the catalog's `pairs.tsv` (fig6.md:188). | Use it only with the definition below, or replace with "test". |
| **arm** | Every method, tool, Rustle configuration and mode (the whole of figs 1–7). | "method" for tools; "configuration" for Rustle's two seedings; "mode" for de novo / guided. |
| **row** | (1) A row of `copy_assign`'s per-family table (figs 4–5, "true family's row"). (2) A comparison row in fig. 7. (3) A plot row in fig. 5. | For (1), write "the read's result for its source family". |
| **F** | Fig. 7b: **pooled** F, the harmonic mean of pooled sensitivity and precision. Fig. 7c–f: **per-family** F. | Label them "F (pooled over families)" and "per-family F". |

### A2. Factual mismatches between text and code

1. **Fig. 6 says "1,418 copies in 258 families". The catalog has 290 families** (fig6.md:104,
   `fig6_chr16_precision.tsv:19`, `fig6_chr16_groups.tsv:23`). `${work}/fig6/chr16.cat.copies.tsv` and the
   fig. 4 catalog `chr16_arm/on.copies.tsv` have identical rows: 1,418 copies and 290 family ids. 258 is the
   number of families with at least one copy that overlaps an annotated RefSeq gene on chr16. That is what
   `score.py:load_members` counts: it keeps only gene-mapped copies. I re-counted it and got 258. Write "290
   families, 258 of which have a copy on an annotated gene".
2. **"Decisive PSV" also counts splice junctions.** `n_decisive` counts the PSV columns the read spans where the
   candidate copies differ, **plus** the read's splice-junction boundaries that some candidate copies have and
   others lack (`copy_assign.rs:read_copy_evidence`). The UpSet set names, the fig. 5 legend and the claims say
   "decisive PSV"; parts of the fig. 4 caption say "decisive site". Use "decisive site" everywhere, with the
   definition in Part B.
3. **The fig. 6a nucleotide-edge coverage is not "of the shorter representative".** The spectrum scorer's default
   is query span ÷ shorter length (`lib.paf_coverage 'query_over_min'`), which exceeds 1 when the query is the
   longer sequence. The protein tier uses **max(query coverage, target coverage)** ≥ 0.50 (`score.py:cmd_spectrum`,
   `cov = max(qcov, tcov)`); the pre-registration said query coverage. fig6.md:52–53 and
   `fig6_chr16_recall.tsv:17` ("qcov>=0.50") state the rule wrongly. The caveat at fig6.md:246–248 describes the
   nucleotide rule correctly.
4. **"identity nm/bl"** (fig6.md:247) is PAF column 10 ÷ column 11 (**matching bases ÷ alignment block length**).
   It is not NM, the edit distance. Rename it.
5. **"equally good placements"** (fig3.md:41) is wrong: tied means a second alignment scores **within 2%**
   (≥ 0.98 × best).
6. **"Placed uniquely (MAPQ > 0)"** (`_o2.SET_LABEL`, the fig. 4 UpSet) does not follow from MAPQ > 0. MAPQ 1–59
   still allows other placements. Write "MAPQ > 0".
7. **"Compara group"** is our connected component of the 290 universe pairs (`_o1.chr16_attribution`,
   `fig6_chr16_groups.tsv:17`), not an Ensembl Compara gene tree or family. Say so.
8. **The Ensembl Compara release is not recorded.** The table was exported from BioMart (`useast.ensembl.org`,
   `hsapiens_paralog_*`) on 2026-09-24 (docs/PREREG_identity_spectrum_2026-09-24.md). A publication needs the
   release number. Record it, or re-export from a numbered release.
9. **(Resolved 2026-09-25: IsoSeq collapse 26.2.0, `isoseq_upload/isoseq_A119b/logs/collapse_A119b_55435595.out:7`
   and the gorilla log.) The IsoSeq collapse version is stated nowhere.** SQANTI3 is "5.5" in README.md:21 and "5.5.4" in fig2.md:3.
10. **Judgeable pairs make precision an upper bound in fig. 6 too.** For Compara, a pair is judged only when both
    genes have at least one Compara paralogue row (any chromosome), so pairs that involve a gene Compara lists
    with no paralogue are never judged. For the protein-homology families, a pair is judged only when both genes
    are members of a family of ≥ 2 genes on that contig. Pairs that involve a non-coding gene, a pseudogene or a
    protein with no homologue are never judged. Fig. 7 states this ("precision is an upper bound"); fig. 6 does
    not.
11. **Closest sibling** is the most similar copy **among the copy's direct alignment neighbours** in the catalog
    family (`max_family_identity` = the highest E_r identity on an incident edge). It is not the most similar
    copy over all family members.
12. **Soto 2025 genes are matched by exact gene name.** Table S1C uses GENCODE-style names (e.g. `AC134879.1`),
    and they are matched to RefSeq `Name=` values in the CHM13 RefSeq GFF (`_o1_recovery.truth_table`,
    `gene_spans`). Genes named differently in the two annotations drop out silently. The caption does not say
    this.
13. Fig. 7 keeps the **first** Soto family ID of each gene (`_o1_recovery.truth_table`, as `family_score` does),
    but `bench/lib.py:soto_gene_family` **drops** genes with more than one ID. Fig. 7's text matches its own code.
    Other bench scorers treat Soto differently, so numbers from them are not interchangeable.
14. **The protein-homology families exclude more than V(D)J segments.** They also exclude constant regions
    (`C_region`) and every `IG_*` / `TR_*` biotype (`truth.py:excluded`, rule 2). The legend should say
    "immunoglobulin and T-cell-receptor gene segments".

### A3. Scope words (relevant to the genome-wide caching decision)

"Genome-wide" currently has two meanings. **(a) Every annotated contig was assembled and scored:** figs 1 and 3,
gorilla only. **(b) Reads were mapped to the whole genome, but only one chromosome's catalog is scored:** figs 4–5
("Genome-wide assignment of the tied reads is abstained", fig4.md:18; human chr16 catalog, gorilla NC_073244.2
catalog). Scopes today:

| figure | human (A119b) | gorilla (OR6737) |
|---|---|---|
| 1 | chr20–22 (every method) in the current tables; genome-wide (24 annotated contigs) after the rebuild | genome-wide |
| 2 | chr20–22; genome-wide after the rebuild | chr20, chr22, chrY; genome-wide after the rebuild |
| 3 | chr20–22; genome-wide after the rebuild | genome-wide |
| 4, 5a–b | chr16 catalog | chr20 (NC_073244.2) catalog |
| 5c | chr16 NPIP locus | none |
| 6 | chr16 | chr20 (NC_073244.2) |
| 7 | chr2, chr6, chr8, chr10, chr16 | chr10 (NC_073234.2), chr20 (NC_073244.2) |

Always print the scope next to "genome-wide" when it means (b). Gorilla contig accessions, checked against
`GGO_genomic.gff` region lines: NC_073224.2 = chr1, NC_073234.2 = chr10, NC_073244.2 = chr20, NC_073246.2 =
chr22, NC_073248.2 = chrY.

---

## Part B. Glossary

Entries are grouped by figure family: 1. data and methods, 2. transcripts (figs 1–3), 3. copy assignment
(figs 4–5), 4. families (figs 6–7).

### 1. Data, references and methods

#### Rustle (default configuration)
- **Code.** `tools/rustle_pipeline.sh assemble` (stage_assemble, lines 76–87). It runs
  `copy_assign --assemble-only --genome-wide --assembly-junctions strict --assembly-polish full
  --polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7
  --polish-retained-ratio 10 --gtf-tpm`. By default it also sets `RUSTLE_GTF_SECONDARY=1`,
  `RUSTLE_GTF_SECONDARY_AS_RATIO=0.98` and `RUSTLE_GTF_SECONDARY_AS_TABLE=<prefix>.molecules.tsv` (see
  *Seeding loci with secondary alignments*). A locus-building read chain needs ≥ 2 reads
  (`denovo_assemble.rs`, `PASS1_MIN_READS = 2`; the O1 node floor is 2).
- **Legend.** "Rustle: loci built from primary alignments plus secondary alignments that score within 2% of the
  read's best alignment score."
- **Occurs.** Every figure; labels in `figlib.TOOL_LABEL` (figlib.py:28–35).

#### Rustle (primaries only), also "Rustle (prim.)"
- **Code.** The same driver with `--no-seed-secondaries` (= `--seed-pool primary`, `SEED_POOL=primary` in rustle_pipeline.sh). The seeding
  variables are unset, so the assembler's read pool holds primary alignments only.
- **Legend.** "Rustle, primary alignments only (no secondary alignments used)". In dense panels, "Rustle (prim.)"
  is defined once in the legend.
- **Occurs.** fig1.md:29,51,64,68,79; fig2.md:29,49,62; fig3.md:9,16,22,73,80,91,99,112,153,184,192;
  fig6.md:40,126,128,134,156,159; fig_intron_chain.py:481–482 (label), 553, 559, 650; fig_secondary.py:1146;
  fig_family_spectrum.py:589, 612 (label "primaries only"); `data/fig1_*.tsv:33`; fig6_gorilla_recall.tsv:31.

#### Seeding loci with secondary alignments (the 0.98 × genome-wide best AS rule); "tied secondaries", "near-tied secondaries"
- **Code.** `as_table` (src/bin/as_table.rs:1–17) reads the **whole** BAM once: every mapped non-supplementary
  record, primary and secondary. For each read it writes `best_as` (the maximum AS over those records) and
  `second_as`; a missing AS counts as 0. The streaming assembler (denovo_assemble.rs:1452–1515) **drops a
  secondary record only when** it has an AS, the read's best AS is known and > 0, and
  `AS < 0.98 × best_as`. Secondaries of reads missing from the table are kept.
  `fig_secondary.is_candidate` (lines 236–241) mirrors this rule exactly. The rule is only genome-wide when the
  table is built from the full BAM (rustle_pipeline.sh:31–32).
- **Legend.** "Secondary alignments whose alignment score is at least 98% of the read's best alignment score
  anywhere in the genome are added to the reads used to build loci."
- **Occurs.** fig1.md:16,27,28; fig2.md:27,119; fig3.md:3,7,35,40,128,138,165,214; fig6.md:17,127,170;
  fig7.md:47,175; fig_family_spectrum.py:589, 612 (on-figure "+ tied secondaries"); fig_secondary.py:1178
  (on-figure "Secondary ≥ 0.98 × best"); `data/fig1_*.tsv:33`; `data/fig3_*.tsv:23`; fig6_gorilla_recall.tsv:31.
- **Flag.** "Tied secondaries" (fig. 6d panel, fig. 7) and "near-tied secondaries" (fig2.md:119) name the same
  rule. Pick one: "secondary alignments within 2% of the best score".

#### StringTie, FLAIR, IsoSeq collapse (the other methods)
- **Code.** GTFs that the lab produced on the same BAMs (`samples.tsv` columns `stringtie_gtf`, `flair_gtf`,
  `isoseq_gff`; human A119b and gorilla OR6737 only). All three are **annotation-free** runs, checked in the
  scripts on 2026-09-25:
  - StringTie 3.0.1: `stringtie -L -p N -o OUT.gtf -A OUT.abund BAM` with **no `-G`**, on the full BAM with its
    secondary alignments (`benchmark_collapse/run_stringtie.sbatch:44`).
  - FLAIR 3.0.1: FLAIR's own `dofiltering()` turns the existing BAM into BED12 (`bin/flair_bam2bed.py`; primary
    alignments only, so the alignments are the same as everyone's), then `flair collapse -q BED -g GENOME -r READS
    --trust_ends --generate_map` with **no `-f`/`--gtf`**; `flair correct` is skipped because it needs an annotation
    or short-read junctions (`run_flair.sbatch:9–21, 73–86`).
  - IsoSeq collapse 26.2.0 (`isoseq collapse`, default settings), which reads no annotation.
  Every method's GTF and the annotation are cut to the same contigs (`assembly.py`).
- **Legend.** "StringTie 3.0.1 (`-L`, no annotation), FLAIR 3.0.1 (collapse without annotation; `flair correct`
  skipped) and PacBio IsoSeq collapse 26.2.0, run by the lab on the same alignments."
- **Occurs.** Figs 1–3 and 5c; labels in figlib.py:28–35.

#### Modes: transcript modes (figs 1–3), family modes (figs 6–7), locus modes (fig. 8)
Three different things are called a "mode". They are never mixed, and every panel names its own.
- **Transcript modes (figs 1–3; `assembly.MODE_DENOVO_METHODS`, `figlib.mode_lines`).**
  - *Annotation-free (de novo)*: the method sees the reads and the genome, never the annotation. Every tool
    comparison of figs 1–3 is in this mode on every side: Rustle `assemble` (`tools/rustle_pipeline.sh assemble`,
    reads + genome), StringTie `-L` without `-G`, FLAIR collapse without annotation (`flair correct` skipped),
    IsoSeq collapse. Printed on each figure: "Annotation-free (de novo) comparison: Rustle assemble (reads +
    genome), StringTie -L without -G, FLAIR collapse without annotation (flair correct skipped), IsoSeq collapse."
  - *Annotation-guided*: the tool is given the annotation (StringTie `-G`, FLAIR with the annotation). These runs
    are supplied by the user (`samples.tsv` columns `stringtie_guided_gtf`, `flair_guided_gtf`; `-` for every
    sample today) and are scored by the same code into **separate** tables and figures (`fig1_guided`,
    `fig2_guided_*`, `fig3_guided_bins`; `fig1g_guided`, `fig2g_guided`, `fig3g_guided`), never in a panel with the
    annotation-free methods (docs/PREREG_guided_transcript_comparison_2026-09-25.md). Until a guided GTF exists
    every figure prints "Guided comparison: not available (guided StringTie/FLAIR GTFs not supplied)."
  - **Rustle has no annotation-guided transcript assembly.** `copy_assign` reads no annotation while assembling
    (its `--gff` only tags catalog copies as annotated or not in the `--phase` copy graph), no `RUSTLE_*` switch
    feeds one to the assembler, and the driver's `assemble` takes `--bam --fasta` only. The guided transcript
    comparison therefore has no Rustle row, and it describes the guided tools against the annotation they were
    given (not comparable with any annotation-free number).
- **Family modes (figs 6–7; O1): de novo vs guided.** See *De novo mode; guided mode* in section 4. Here "guided"
  means the loci are the annotated gene and pseudogene bodies (annotation + genome, **no reads**); it is
  Rustle-internal (two Rustle modes compared with each other against reference families), not a tool comparison,
  and it has nothing to do with StringTie's `-G`.
- **Locus modes (fig. 8).** Rustle's guided candidate search (annotation + genome) is compared like for like with
  Liftoff's self-lift (`-copies`) of the same annotation; Rustle's de novo loci (reads + genome) are scored against
  that baseline as a reference, never as a competitor (captions/fig8.md).
- **Legend.** Name the mode in every panel title or caption sentence that compares methods: "annotation-free (de
  novo)" or "annotation-guided" for transcripts; "de novo family mode" / "guided family mode" for families.
- **Occurs.** captions/fig1.md, fig2.md, fig3.md (Mode paragraph); on-figure first line of figs 1–3
  (`figlib.mode_lines`); table notes `mode:` of `fig1_*`, `fig2_*`, `fig3_*` (after the next build); fig7.md;
  fig8.md.

#### Rustle (copy_assign) (fig. 5c)
- **Code.** `copy_assign --bam hsa16.bam --fasta chm13v2.0.fa --region chr16:11963320-80438591 --families
  copies16.tsv --copies-fa copies16.fa --gtf --origin-drop-indels --threads 4` (fig5_hard_locus.tsv:15). It uses
  the transcripts of the copy-assignment program, not the assembler of figs 1–3.
- **Legend.** "Rustle, copy-assignment step (its own transcripts; not the assembler of Figs 1–3)".
- **Occurs.** fig5.md:21,76,148; fig_assign_accuracy.py:12, 41, 88 (tick label), 252, 543;
  fig5_hard_locus.tsv:24.

#### Reference annotation (what figs 1–3 call "truth")
- **Code.** Gorilla: `GGO.GCF_029281585.2_RefSeq.gtf.gz`. Human: `A119b.chm13v2.0_RefSeq.gtf.gz`. Chimpanzee and
  orangutan: their RefSeq GFF3 converted with `gff_to_gtf` (`assembly.annotation_gtf`). Every figure 1–3 build is
  genome-wide: every contig the sample's annotation covers (human: all but chrM, which it does not annotate); the
  current human tables are still chr20–22 until the rebuild (captions' Status boxes).
- **Legend.** "RefSeq annotation (gorilla GCF_029281585.2; human T2T-CHM13 v2.0)". Do not call it truth: an
  annotation is a reference, not ground truth.
- **Occurs.** fig3.md:5,146,148; fig_secondary.py META (line 118, "truth: GGO RefSeq").

#### Read, molecule, alignment, placement, primary, secondary, supplementary
- **Code.** A read is one IsoSeq full-length read, i.e. one cDNA molecule. A placement is one BAM record.
  "Primary" means flag 0x100 and 0x800 unset; per-read statistics filter with `-F 2308` (fig_intron_chain.py
  `BAM_EXCLUDE = 2308`). Supplementary (0x800) and unmapped (0x4) records never count (fig_secondary.py). Fig. 4
  placements are every `-F 2052` record (`_o2.placements`).
- **Legend.** Use "read" for the molecule and "alignment" for the record. Figs 1 and 3 use "reads" and
  "molecules" for the same thing: pick "reads".

### 2. Transcripts (figs 1–3)

#### Intron chain; exact intron-chain match; gffcompare `=`
- **Code.** A read's intron chain is the ordered list of CIGAR `N` spans, 0-based half-open
  (`fig_intron_chain.read_intron_chain`). A reference chain is matched by a method when gffcompare 0.12.10
  (`-r ref.gtf`) gives class code `=` (every intron identical, ends free) to any reference transcript that
  carries that chain. In fig. 3, a match also counts for reference transcripts with an identical chain,
  keyed on chrom + strand + introns (`fig_secondary.matched_by_tool`). In fig. 5c a mono-exonic read matches by
  span containment.
- **Legend.** "Exact intron-chain match: the same introns, with identical coordinates, in the same order;
  transcript ends may differ."
- **Occurs.** fig1.md:4,40,43,48,81; fig2.md:14,107; fig3.md:5,27,50,112,132,158,216; fig5.md:74,97,99;
  fig_secondary.py:1063 (on-figure "gffcompare '='"); fig_assign_accuracy.py:499, 507 (on-figure);
  `data/fig3_*.tsv:24`; fig5_hard_locus.tsv:22.

#### Intron-chain sensitivity / precision (fig. 1a,d,b,e)
- **Code.** gffcompare's intron-chain level. Sensitivity = matching reference chains ÷ multi-exon reference
  transcripts. Precision = matching query transcripts ÷ multi-exon query transcripts
  (`fig_intron_chain.py` docstring lines 5–10).
- **Legend.** As above. The denominator is multi-exon transcripts only; state it.
- **Occurs.** fig1.md:4–17, 35–41; fig_intron_chain.py:571–572, 645–646 (axis labels).

#### Read-support strata and "Rustle's admission floor" (fig. 1c,f,g,h); on the figure: "Rustle's minimum read support"
- **Code.** For each distinct multi-exon reference chain (contig + chain), the number of **primary** alignments
  (`-F 2308`) whose intron chain equals it; strand is not used (`count_exact_chains`, `support_rows`). The
  strata are all, ≥ 1, ≥ 2 and ≥ 5. In count mode `reads` each alignment counts once. In count mode
  `distinct_ends` each distinct (start, end) span counts once, which is the key the assembler uses to remove
  coordinate duplicates. "Admission floor" means Rustle reports a chain only when ≥ 2 reads support it (its own
  count after junction correction, which is not this count).
- **Legend.** "Reference intron chains carried exactly by ≥ k primary alignments. The ≥ 2 stratum equals
  Rustle's minimum read support, so this stratum is matched to Rustle's design."
- **Occurs.** fig1.md:11,43–49,86,108; fig_intron_chain.py:657, 690 (on-figure "Rustle's floor"), 700 (axis);
  fig1_support.tsv:35.

#### Paired comparison: Tango interval, exact McNemar p
- **Code.** For each ≥ k stratum: 2×2 table of chains matched by both / Rustle only / the other method only /
  neither. The difference is Rustle's sensitivity minus the other method's, with Tango's asymptotic score 95%
  interval. p is the exact two-sided McNemar test, 2·P(X ≤ min(b, c)) with X ~ Bin(b + c, ½)
  (`fig_intron_chain.tango_ci`, `mcnemar_log10p`). Chains are treated as independent and p values are not
  adjusted.
- **Legend.** "Difference in sensitivity on the same reference chains (percentage points; Tango 95% interval);
  exact McNemar test on the discordant chains."
- **Occurs.** fig1.md:56–71; fig_intron_chain.py:757 (axis), 827 (panel title).

#### SQANTI3 structural categories: FSM / ISM / NIC / NNC / Other (incl. genic intron)
- **Code.** SQANTI3 5.5.4 `sqanti3_qc.py` against the per-contig RefSeq GTF (`_sqanti.py`). The categories are
  folded by `figlib.SQANTI_CATEGORY_MAP` (figlib.py:71–76):
  - FSM, full-splice_match: every splice junction matches a reference transcript with the same number of
    junctions.
  - ISM, incomplete-splice_match: consecutive reference junctions, but fewer 5′ or 3′ exons.
  - NIC, novel_in_catalog: a new combination of annotated splice sites or junctions.
  - NNC, novel_not_in_catalog: at least one unannotated donor or acceptor.
  - Other: antisense, intergenic, genic (overlaps exons and introns), genic_intron (entirely inside an intron),
    fusion.
  
  Subcategories: an FSM `reference_match` has both ends within 50 bp of the reference transcript's ends. An ISM
  `3prime_fragment` shares the reference's last junction but not its first, so it lacks the 5′ exons.
  Multi-exon means the SQANTI3 `exons` column is > 1 (fig2_sqanti_subcategories.tsv:95).
- **Legend.** Spell out the four names at first use (figlib.SQANTI_LABEL is already on the figure). Write
  "genic intron (entirely within an annotated intron)".
- **Occurs.** fig2.md throughout (5–13, 34–70, 150–160, 178–243); fig_sqanti.py:262, 286, 321 (on-figure);
  figlib.py:63–76.

#### Rules-filter PASS (fig. 2b,e)
- **Code.** `sqanti3_filter.py rules` with SQANTI3's default `filter_default.json`, read from
  `SQANTI3/src/utilities/filter/filter_default.json`. PASS means `filter_result == "Isoform"`. An FSM passes if
  `perc_A_downstream_TTS` ∈ [0, 59]. Any other transcript passes if `perc_A_downstream_TTS` ∈ [0, 59],
  `RTS_stage` is FALSE and all junctions are canonical, **or** if it meets the same two conditions and has
  short-read junction coverage ≥ 3 (`min_cov`). No short reads were supplied, so only the first branch can pass.
- **Legend.** "Passes SQANTI3's default rules filter: ≤ 59% A in the 20 bp downstream of the 3′ end; for
  non-FSM transcripts also all-canonical junctions and no reverse-transcriptase template-switching signature.
  Not a correctness measure."
- **Occurs.** fig2.md:1,5,44,53,57,76,116,120; fig_sqanti.py:276 (on-figure "Rules-filter PASS").

#### Candidate placement (fig. 3)
- **Code.** A read's primary record, or a secondary with AS ≥ 0.98 × the read's genome-wide best AS
  (`fig_secondary.is_candidate`). This is exactly the pool Rustle seeds from. Supplementary and unmapped records
  never count.
- **Legend.** "Candidate alignment: the primary alignment, or a secondary alignment scoring ≥ 98% of the read's
  best."
- **Occurs.** fig3.md:35,37,137,214; fig_secondary.py:1178 (on-figure); `data/fig3_*.tsv:23`.

#### Counting a read at a transcript (fig. 3; the same rule at gene level in fig. 6d)
- **Code.** A read counts at a transcript when one of its candidate placements has an **aligned block** (a
  gapless M/=/X run) that overlaps one of the transcript's exons. Intron `N` and deletion `D` spans do not count
  (`fig_secondary._count_contig`; `lib.aligned_blocks`).
- **Legend.** "A read counts at a transcript when its aligned bases overlap one of the transcript's exons. A read
  that splices over the transcript does not count."

#### Tie fraction; tied read (fig. 3); tie bins; high-tie
- **Code.** A read is tied when `best_as > 0` and `second_as ≥ 0.98 × best_as`, both genome-wide from
  `as_table` (`fig_secondary.molecule_tied`, lines 135–136). Reads with a single record are never tied. The tie
  fraction of a transcript is n_tied ÷ n_mol over the reads counted at it. The bins are 0, (0, 0.1], (0.1, 0.5],
  (0.5, 0.9] and > 0.9 (`TIE_BINS`, lines 99–105), fixed before any match rate was examined (fig3.md:173–179).
  High-tie means tie fraction > 0.5.
- **Legend.** "Tie fraction: the share of a transcript's reads that have a second alignment scoring ≥ 98% of
  their best alignment score (genome-wide)." Do not write "equally good".
- **Occurs.** fig3.md:14,40,44,54,56,64,66,67,97,108,119,127,183; fig_secondary.py:983–984, 1055, 1133
  (on-figure); fig3_bins.tsv:28; fig3_example.tsv:23.

#### Expressed (fig. 3)
- **Code.** A multi-exon RefSeq transcript with n_mol ≥ 1, i.e. ≥ 1 read counted at it
  (`fig3_ref_tie` keeps only these). Gorilla has 88,387 of 95,833; human 9,711 of 9,846
  (`fig3_bins.tsv:26–27`).
- **Legend.** "Transcripts with ≥ 1 read (primary or candidate secondary alignment on an exon)". See A1 for the
  other meanings of "expressed".
- **Occurs.** fig3.md:6,50,51,166,174; fig_secondary.py:1060 (on-figure "expressed multi-exon RefSeq
  transcripts"); `data/fig3_*.tsv:26–27`.

#### Reached by a primary alignment / "reached by no primary alignment" (no-primary)
- **Code.** `n_mol_primary` is the number of reads counted at the transcript through their **primary** record.
  "No primary" means `n_mol_primary == 0`. Those transcripts form one facet, whatever their tie fraction
  (`PRIMARY_STRATA`, `PANEL_A`).
- **Legend.** "No primary alignment overlaps the transcript's exons; only secondary alignments do."
- **Occurs.** fig3.md:3,8,20,30,42,54,64,94,103,106,108,118,132; fig_secondary.py:1056, 1106 (on-figure).

#### Read-sharing group ("read group"); independent unit
- **Code.** High-tie transcripts (tie fraction > 0.5) that share any tied read are joined by union-find over the
  whole streamed genome (`fig_secondary.read_sharing_groups`, lines 394–414). A transcript's **unit** is its
  read-sharing group if it has one, otherwise its gene (`cluster_of`, lines 518–522).
- **Legend.** "Read-sharing group: transcripts with > 50% tied reads that share at least one tied read. Counts are
  given as transcripts (read-sharing groups or genes)." Do not print "read group" (the SAM term).
- **Occurs.** fig3.md:9,15,17,22,25,29,46,61,62,65,67,81,94,97,108,128,183,198,201,203,204,213; fig6.md:25,161;
  fig_secondary.py:1029–1030, 1109, 1133 (on-figure "read group(s)"); fig3_gain.tsv:29; fig3_bins.tsv:29.

#### Cluster interval / cluster bootstrap (figs 3, 5)
- **Code.** Fig. 3: a 95% percentile bootstrap that resamples units (read-sharing groups or genes) with
  replacement, 2,000 resamples shared by all methods. When no transcript or every transcript is matched, a
  Wilson interval with n = number of units (`cluster_intervals`; fig3_bins.tsv:29). Fig. 5: the unit is the
  **source copy**; see *Copy-level 95% interval*.
- **Legend.** "95% interval from resampling read-sharing groups (or genes)."
- **Occurs.** fig3.md:8,17,76,184,212; fig5.md:38; fig_secondary.py:1063 (on-figure "95% interval over read
  groups or genes").

#### "Any baseline" (fig. 3c,d)
- **Code.** The union of StringTie, FLAIR and IsoSeq collapse: a transcript is matched by "any baseline" if at
  least one of the three matches it (fig3_gain.tsv:28).
- **Legend.** "Matched by at least one of StringTie, FLAIR or IsoSeq collapse." The panel already prints this at
  fig_secondary.py:1125.
- **Occurs.** fig3.md:18,90,92,99,184,209; fig_secondary.py:1079, 1125.

#### Example-locus rule (fig. 3e)
- **Code.** `fig_secondary.pick_example`: among no-primary transcripts that Rustle matches with `=`, keep those
  with tie fraction > 0.5 and all annotated introns ≥ 20 bp, and take the one with the most reads. Ties go to the
  shortest span, then the id.
- **Legend.** "Chosen by a fixed rule: the no-primary transcript Rustle rebuilds exactly with the most reads
  (tie fraction > 0.5, all introns ≥ 20 bp)."
- **Occurs.** fig3.md:118,185; fig3_example.tsv:23.

### 3. Copy assignment (figs 4–5)

#### Default de novo families; family copy table; positional exon sum (figs 4–8)
- **Code.** The ONE default de novo family definition (user decision 2026-09-25; `docs/seeded_family_definition.md`
  §0★; `docs/THESIS_OBJECTIVES.md` O1):
  1. `tools/rustle_pipeline.sh assemble`: loci from primary alignments plus secondary alignments scoring ≥ 98% of the
     read's genome-wide best alignment score; a locus is one `gene_id`.
  2. The **representative** of a locus is its transcript with the most reads (ties: the longer span). Its **positional
     exon sum** = its exon coordinates (from the reads) with the genome's bases at them (`docs/PSEUDOCODE_2026-09-08.md`:
     "the reads supply the coordinates, the genome supplies the bases").
  3. `rustle_pipeline.sh families` = `mcl_families --from-gtf --min-exonic-bp 1 --min-shared-exon-frac 0.60`: the
     family rule of *De novo mode; guided mode* below, on the loci's genomic spans; MCL (inflation 2.8); families of
     ≥ 2 loci.
  4. With `--emit-units` the stage writes the **copy table** `<id>.fam.copies.tsv` / `.fa`: one copy per member locus
     (its representative's exons, 0-based half-open, and its spliced exon sum), in the `copy_assign --families`
     contract (columns 1–11 as the legacy catalog's: `family_id copy_idx tid chrom start end n_exon strand n_reads
     exons max_family_identity`; `max_family_identity` = the best families-stage edge identity, genomic spans, not the
     catalog's exon-sum identity). Copy assignment and every figure that needs copies read this table
     (`samples.product(cfg, id, 'families', 'copies')`).
- **Legend.** "Rustle's default de novo families: IsoSeq reads are assembled into loci (secondary alignments within 2%
  of the read's best score also seed loci); each locus is represented by its most-read transcript's exons; loci are
  grouped when their genomic sequences align (≥ 70% identity) and share ≥ 60% of the smaller locus's exonic sequence
  (Markov clustering). No protein is used."
- **Occurs.** fig6.md, fig7.md, fig8.md; fig_family_spectrum.py (`FAMILIES_WHAT`, `FIG_NOTE`); fig_loci.py (`families`
  rows); `_o1.families_copies`; README.md.

#### Copy catalog (LEGACY); catalog copy; catalog family
- **Status.** Legacy since 2026-09-25: kept runnable (`catalog` stage) and scored only as a labelled comparison (fig. 8
  table rows `catalog`); not the default, not a headline. The fig. 4–5 development tables built before that date used
  it.
- **Code.** `gw_family_catalog` (`tools/rustle_pipeline.sh catalog`; src/bin/gw_family_catalog.rs:1–12). Its
  nodes are collapsed read-supported locus representatives (spliced sequences). A family is a γ-quasi-clique of
  the homology graph E_r spanning ≥ 2 distinct loci. There is an E_r edge when minimap2 aligns two
  representatives with coverage ≥ 0.50 of the shorter sequence and identity (1 − `de`) ≥ 0.80 (`asm20`) or
  ≥ 0.60 (sensitive tier, `-k11 -w5`) (`denovo_pipeline.rs` `RefineParams` defaults, lines 5148–5162). Output
  `copies.tsv` has one row per copy: `family_id copy_idx tid chrom start end n_exon strand n_reads exons
  max_family_identity`. The fig. 4 human catalog has 1,418 copies in 290 families; gorilla NC_073244.2 has 357
  copies in 54 families.
- **Legend.** "Copy catalog: loci built from the reads and grouped into families of copies whose spliced
  sequences align (≥ 60% identity over ≥ 50% of the shorter copy)."
- **Occurs.** fig4.md:24,29,48,56,70,78,87,88,89,95,98,119,120; fig5.md:48,49,53,64–81,145;
  fig6.md:6–14,31,35,51–70,76,83,104–105,171,186–187,221,236,238; fig_assignability.py:98–100 (notes);
  fig_assign_accuracy.py:15; `data/fig4*`, `data/fig5*`, `data/fig6_chr16_*`.

#### Simulated read; source copy; truth (figs 4–5a,b; the only literal ground truth in the figures)
- **Code.** `bench/sim.py copies` (random seed 20260925) simulates reads from the spliced sequence of every
  catalog copy ≥ 300 bp that is in a family of ≥ 2 copies. Each copy gets min(100, max(10, n_reads)) reads with
  the HiFi error model (0.1% substitutions, 0.03% deletions, 0.03% insertions, up to 10% trimmed from one end,
  then 0–30 bp from each end). Reads are named `family|copy|i` and mapped genome-wide with
  `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes` (fig4.md:94–102).
- **Legend.** "Truth = the copy each simulated read was drawn from (its source copy)."
- **Occurs.** fig4.md:25,96; fig5.md:3,32; fig_assignability.py:218 (axis "Simulated reads").

#### MAPQ 0 ("tied reads", "Aligner tie (MAPQ 0)")
- **Code.** The simulated read's primary record has MAPQ 0 (`score.py:cmd_reads` keeps MAPQ-0 reads only;
  `_o2.SET_LABEL['mapq0']`). MAPQ 0 is the aligner's own mapping-quality estimate. It does **not** always mean
  equal alignment scores: 3 of the 1,263 human MAPQ-0 reads have the best AS on the true copy alone
  (`fig4_assignability_upset.tsv:31`, `no_equal_as_elsewhere`).
- **Legend.** "Reads whose primary alignment has mapping quality 0 (the aligner could not choose a placement)."
- **Occurs.** fig4.md:5,39,97; fig5.md:31; `_o2.py:61` (UpSet label); fig_assign_accuracy.py:532 (panel titles
  "simulated reads tied at MAPQ 0").

#### AS-tied gate (copy_assign; fig. 5c hard set, and why MAPQ > 0 reads get no row)
- **Code.** `copy_assign` admits a read to copy assignment only when it has ≥ 2 placements in the swept region
  (primary and secondary; supplementaries excluded) and `second_AS ≥ as_tie_ratio × best_AS`, with
  `as_tie_ratio` = **1.0 by default** (an exact tie) (src/bin/copy_assign.rs:537–541 flag, 2027–2032 `as_tied`,
  3260 gate). Other reads are dropped before the test and get no row. Placements on other contigs are invisible
  to the gate.
- **Legend.** "Reads with ≥ 2 alignments in the region whose best alignment scores are exactly equal."
- **Occurs.** fig4.md:37; fig5.md:79; fig_assign_accuracy.py:541 (on-figure "AS-tied"); fig5_hard_locus.tsv:21.

#### PSV (paralogous sequence variant)
- **Code.** For each catalog family, `copy_assign` aligns every copy's spliced sequence to the family's first
  copy (`copy_assign_pipeline.rs:discover_psvs`, line 1067). A **PSV column** is an alignment position where at
  least one copy's base differs. Every copy carries its own base there (a gap counts as missing). The family's
  PSV columns form a bubble graph (`copy_assign.rs:BubbleGraph`).
- **Legend.** "PSV (paralogous sequence variant): a position at which the copies of a family carry different
  bases."
- **Occurs.** fig4.md:3,4,6,7,8,9,42,44,60,61,62,64,72; fig5.md:7,9,10,11,36,43,55–59,68,108;
  fig_assign_accuracy.py:98 (legend); `_o2.py:63–64` (UpSet labels).

#### Decisive site ("decisive PSV", "decisive family PSV")
- **Code.** `n_decisive` (`copy_assign.rs:read_copy_evidence`, lines 492–567) counts two things. (a) The PSV
  columns the read spans at which the candidate copies carry ≥ 2 different non-missing bases (`Bubble.decisive`,
  line 39). (b) The read's splice-junction boundaries that are present (within 4 bp) in some candidate copies
  and absent in others. **Junctions are included.** "≥ 1 decisive PSV within its family" means `own_n_decisive
  ≥ 1` in the per-family run's result for the read's source family. "No decisive PSV over all placements" means
  the union run's result for that family has `n_decisive == 0` (`_o2.per_read`, lines 402–414).
- **Legend.** "Decisive site: a PSV or a splice junction, covered by the read, at which the candidate copies
  differ." Write "decisive site", not "decisive PSV", because junctions count too.
- **Occurs.** fig4.md:4,6,8,9,42,43,44,61–64,72; fig5.md:7,9,10,11,36,55–59,68,108; fig_assignability.py:269–270
  (on-figure); fig_assign_accuracy.py:98 (legend); `_o2.py:63–64` (UpSet); fig4_assignability_upset.tsv:30,41,44;
  fig5_assign_accuracy_bands.tsv:41.

#### PSV certificate (the per-family copy-assignment test); statuses assigned / ambiguous / tied
- **Code.** `copy_assign.rs:assign_read_editing` (lines 581–695) with the defaults in `AssignParams`
  (lines 336–341): error rate 0.003, junction weight ±5, junction tolerance 4 bp, α = 1e-3, junction error 1e-4.
  1. **Best copy.** Log-likelihood over the decisive sites (a base match scores ln(1 − e), a mismatch
     ln(e/3), using per-base quality when present), plus ±5 per junction. The best copy is the argmax.
  2. **Test.** Against each competitor c, the p-value is the Poisson-binomial probability that sequencing
     errors alone would make a read from c agree with the best copy at at least k of their distinguishing sites
     (`copy_pair_significance`, line 392). p_read is the maximum over competitors.
  3. **Decision.** Let thr = α/(n − 1) (Bonferroni over n candidate copies). The read is **assigned** if
     p_read < thr and the log-likelihood margin > 0. It is **tied** (unresolvable) if even perfect agreement
     could not reach thr (`min_p ≥ thr`); this includes reads with no decisive site. Otherwise it is
     **ambiguous**.
  4. **Origin check.** A read whose substitutions (plus indels, unless `--origin-drop-indels`, plus unaligned
     read bases) against the best copy's unit exceed Binomial(aligned bases, error rate) at one-sided α is
     **origin-rejected** and reported as ambiguous (`copy_assign_pipeline.rs:2804–2830`).
- **Legend.** "Copy-assignment test: a read is assigned to its most likely copy only if, against every other
  candidate copy, sequencing error alone would explain its agreement with that copy at the distinguishing sites
  with probability < 0.001/(n − 1). Otherwise the read is left unassigned." Say "test", or define "certificate"
  once with this sentence.
- **Occurs.** fig4.md:3,12,41,42,47,52,67,92,111,112,116; fig5.md:5,10,14,36,43,47,53,63,82,106,134,146;
  fig_assign_accuracy.py:98 (legend); `_o2.py:75` (legend "PSV certificate, row of the read's true family").

#### Sole-candidate assignment
- **Code.** When exactly one copy of the family has an alignment of the read (`n_candidates == 1`) and the origin
  check passes, the read is **assigned** with `n_decisive = 0` (`copy_assign.rs` `sole_candidate` doc,
  lines 285–289; `copy_assign_pipeline.rs:2898`). With ≥ 2 candidates and no decisive site the test cannot
  succeed (`min_p = 1`), so an assigned read with 0 decisive sites is always a sole-candidate read. There were
  10 human reads, all from one copy, `GWFAM221|4` (fig4_assignability_upset.tsv:30).
- **Legend.** "Only one copy of the family aligns the read, so no distinguishing site is needed."
- **Occurs.** fig4.md:6,60; fig5.md:7; META claims fig_assign_accuracy.py:36, fig_assignability.py:37;
  fig4_assignability_upset.tsv:30,41.

#### Per-family table; default output; "no cross-family arbitration"; claimed for a wrong locus; conflict
- **Code.** Default `copy_assign --families` writes `<out>.assignments.tsv` with **one row per (read, catalog
  family)** for every family whose copies the read's tied placements touch. Each family is tested
  independently, so one read can be `assigned` in two families' rows. `score.py reads` ANY (lines 1112–1200)
  reads **any** assigned, not-origin-rejected row. Correct: every claimed copy is the source copy or at the same
  locus (`same_locus`: same chromosome, overlap ≥ 50% of the shorter span). Wrong: none is. **Conflict**: two
  or more loci are claimed and some but not all are correct. "Claimed for a wrong locus by another family's
  row" is `wrong_locus_rows_other_family > 0` (fig4_assignability_upset.tsv:45).
- **Legend.** "Default output: each read is tested separately within each family whose copies it aligns to, so
  one read can be assigned in more than one family."
- **Occurs.** fig4.md:12,40,52,65,66; fig5.md:12,45,46,61,62; `_o2.py:76` (legend "Per-family table, any
  assigned row (default output)"); fig_assignability.py:271, 276 (on-figure);
  fig4_assignability_upset.tsv:44–46; fig5_assign_accuracy_bands.tsv:38.

#### Truth-selected reading; "true family's row" (OWN)
- **Code.** `score.py reads` OWN: only the rows whose `family_id` is the read's **source family** are read. That
  choice uses the simulation's knowledge of the source (`cmd_reads`, `own = [r for r in rows if
  r['family_id'] == t[0]]`). "Assigned in true family's row" means OWN verdict ∈ {correct, wrong, conflict},
  i.e. status `assigned` and not origin-rejected. "… and to the true copy" means OWN = correct.
- **Legend.** "Scored in the read's source family (possible only in simulation; a user cannot choose this
  row)."
- **Occurs.** fig4.md:11,46; fig5.md:5,12,36,43,142; `_o2.py:65,75` (UpSet and legend labels);
  fig_assignability.py:266 (on-figure "true-family row"); fig4_assignability_upset.tsv:30,41;
  fig5_assign_accuracy_bands.tsv:41.

#### Union certificate (`--union-certificate`)
- **Code.** `copy_assign --families --union-certificate` (src/bin/copy_assign.rs:572–590). A read whose tied
  placements touch copies of ≥ 2 families, or ≥ 1 locus outside every catalog copy, gets **one** test over the
  union of candidates: all copies of every family touched, plus one pseudo-copy per outside locus built from the
  genome sequence over the placement's aligned blocks. The winner family's row is `assigned`; every other row is
  `tied`. An outside winner, or no winner, leaves every row unassigned. The option requires the AS-tied gate.
- **Legend.** "Union test: one test per read over all its candidate placements, across families and outside the
  catalog, so a read is assigned at most once."
- **Occurs.** fig4.md:12,42,52,67,92,116; fig5.md:14,47,106; `_o2.py:77` (legend "Union certificate
  (--union-certificate)"); fig_assignability.py:272 (on-figure); fig4_assignability_upset.tsv:44–46;
  fig5_assign_accuracy_bands.tsv:38.
- **Flag.** Drop the flag name `--union-certificate` from the legend (fig5 legend line) and keep it in Methods.

#### Aligner's primary placement (baseline)
- **Code.** The catalog copy with the largest raw overlap of the primary alignment's reference span. Correct if
  it is the source copy or at the same locus. A primary outside every catalog copy is not placed
  (`_o2.per_read`; fig5_assign_accuracy_bands.tsv:38).
- **Legend.** "Aligner: the copy its primary alignment overlaps most."

#### Coverage and accuracy (fig. 5a,b)
- **Code.** coverage = (correct + wrong + conflict) ÷ tied reads. accuracy = correct ÷ (correct + wrong +
  conflict) (fig5_assign_accuracy_bands.tsv:39).
- **Legend.** "Fraction assigned" and "fraction correct among assigned". The current axis labels,
  "Coverage (assigned / tied reads)" and "Accuracy (correct / assigned)" at fig_assign_accuracy.py:443–444, are
  acceptable if the words are kept.

#### Copy-level 95% interval (fig. 5)
- **Code.** `_o2.copy_level_ci` (line 486): a percentile bootstrap over **source copies** (2,000 resamples, seed
  20260925). If every copy has the same ratio, a Wilson interval at the pooled ratio with n = copies assigned.
  No interval when < 5 source copies are assigned.
- **Legend.** "95% interval with the source copy as the unit." The axis already says this
  (fig_assign_accuracy.py:444).

#### Identity band; closest sibling; "100%* over the aligned core"
- **Code.** A read's band is `figlib.identity_band(max_family_identity)` of its source copy (figlib.py:79–93).
  The bands are "identical" (≥ 1.0), 99.5–100%, 99–99.5%, 98–99% and < 98%.
  `max_family_identity` is the highest E_r identity, 1 − minimap2 `de` (gap-compressed per-base divergence;
  fallback matches ÷ block length), on any catalog edge incident to the copy (gw_family_catalog.rs:351–357). An
  E_r edge needs coverage ≥ 0.50 of the shorter spliced sequence. So "100%" means no substitution and no gap over
  the aligned block, which covers ≥ 50% of the shorter copy. It does not mean full-length identity, and it
  compares the copy only with its **directly aligned** neighbours. `score.py reads` bins the same column as
  divergence: < 0.5%, 0.5–1%, 1–2%, 2–5% and ≥ 5%.
- **Legend.** "Identity of the source copy to its most similar directly aligned copy of the same family
  (minimap2, over their aligned segment, which covers ≥ 50% of the shorter copy). 100%* = identical over that
  segment, not necessarily over the full length." The fig. 4 legend already says most of this
  (fig_assignability.py:328–329).
- **Occurs.** fig4.md:8,28–30,59,63,71,104; fig5.md:8,11,32,33; fig_assign_accuracy.py:99 (tick "100% (aligned
  core)"), 441 (axis "Source copy vs closest sibling"); fig_assignability.py:49, 185, 194, 328–329 (on-figure);
  fig4_assignability_upset.tsv:49.

#### NM-identical twin
- **Code.** `_o2.twin_state` (lines 331–354), over every `-F 2052` alignment of the read. A true-copy placement
  overlaps the source copy's catalog span; every other placement is "elsewhere". The read has an NM-identical
  twin when some placement elsewhere has AS equal to the read's best AS **and** the same NM as the true-copy
  placement (the highest-AS true-copy alignment, then lowest NM). The other states are:
  `as_tie_nm_differs` (placements elsewhere reach the best AS, none with that NM), `no_equal_as_elsewhere` (the
  true copy holds the best AS alone), `true_below_best` (the best AS is elsewhere) and `no_true_placement`. The
  measure needs the truth, so a user cannot compute it.
- **Legend.** "NM-identical twin: an alignment at another locus with the same best alignment score and the same
  number of mismatches and indels (NM) as the alignment at the source copy. No base in the read can tell the two
  places apart."
- **Occurs.** fig4.md:13,16,52,73,104,105,116; fig5.md:15,63,68,109; fig_assignability.py:273 (on-figure);
  fig4_assignability_upset.tsv:31,42,47; fig5_assign_accuracy_bands.tsv:41.

#### Hard-locus benchmark; hard molecules / hard set; contested stratum; support-matched control (fig. 5c)
- **Code.** Real A119b reads at the chr16 NPIP family (`copies16.tsv`, 26 copies). `hsa16.bam` =
  `samtools view -b -M -L copyregions.bed A119b.t2t.bam`.
  - **Hard molecules**: the rows of Rustle's own `copy_assign` run, i.e. reads admitted by the AS-tied gate (≥ 2
    placements in the swept window with runner-up AS = best, one of them in a catalog copy), restricted to reads
    whose primary alignment (`-F 2308`) lies in the merged copy intervals (fig5_hard_locus.tsv:21). This set is
    defined by Rustle, not by the data alone.
  - **Contested**: hard rows with `origin_rejected == 0` and `n_candidates ≥ 2`
    (`score.py:cmd_bakeoff_compare`, lines 1482–1490).
  - **Support-matched control (post hoc)**: `--min-mult 2`. Keep only reads whose exact primary-alignment chain
    (chrom + introns) is carried by ≥ 2 reads in `hsa16.bam`. This was added after the pre-registered
    prediction P5 failed.
  - **Fraction carried**: a read is carried when the method's GTF has a transcript with the read's exact intron
    chain, and that transcript overlaps a catalog copy (copy = the one with the largest raw overlap).
    `derived_one` or `derived_multi` in `score.py:cmd_bakeoff_calls`. The junction tolerance is 0 bp for
    Rustle, StringTie and FLAIR, and 5 bp for IsoSeq collapse.
  - **Transcripts overlapping the copies / fraction with ≥ 1 read's exact chain**: the method's transcripts that
    overlap a catalog copy, and the share whose exact intron chain (ends free) is carried by ≥ 1 read of the
    scored set. **Conflation** is a transcript that overlaps ≥ 2 copy intervals.
- **Legend.** "Hard reads: reads at the NPIP locus with ≥ 2 alignments of equal best score, one inside a
  catalogued copy (defined by Rustle's copy-assignment step). Contested: hard reads with ≥ 2 candidate copies
  whose best candidate passes Rustle's origin check. Reads whose exact intron chain occurs in ≥ 2 reads (post hoc
  control). Carried: the method reports a transcript with the read's exact intron chain that overlaps a copy."
  Replace "Support-matched" with "chain carried by ≥ 2 reads".
- **Occurs.** fig5.md:1,20–24,29,72–101,110–141; fig_assign_accuracy.py:84–86 (stratum labels), 499, 506–507
  (axes), 541 (title); fig5_hard_locus.tsv:20–23; fig5_hard_locus_transcripts.tsv:15–16.

### 4. Families (figs 6–7)

#### Protein-homology families (replaces "protein referee"); SECONDARY reference, supplementary figures only
- **Status (2026-09-25).** Not a main-figure reference any more: they appear only as a labelled secondary reference in
  the supplementary figures (fig. 6s-seeding development panel; fig. 7s-protein-homology) and are built genome-wide
  only on request (`fig6s_protein_homology 1`, `fig7_protein_homology 1`). Reasons: most of their same-family pairs
  have no nucleotide alignment (register 1101), they mix fold sharers with paralogues (registers 1030/1031), they
  exclude pseudogenes, and they are circular for the protein step.
- **Code.** `bench/truth.py:protein_referee` (lines 161–192) with `protein_edges` (122–147), `write_proteins`
  (150–158), `excluded` (48–55) and `blastp_all_vs_all` (72–83); `bench/lib.py:longest_cds` (148–175),
  `gene_biotypes` (116–129), `translate_refseq` (74–81) and `mcl` → `bench/mcl_port.py` (Rust `mcl_port`). One
  chromosome at a time:
  1. **One protein per gene.** For every `gene=` symbol with ≥ 1 CDS feature on the chromosome, take the
     transcript with the most CDS bases. Exclude genes whose `gene_biotype` contains `pseudogene`, is
     `V_segment` / `D_segment` / `J_segment` / `C_region`, or starts with `IG_` / `TR_`.
  2. **Translate.** Concatenate the CDS segments, ignoring phase, and reverse-complement on the minus strand. A
     protein starting with M is cut at its first stop; otherwise every stop becomes X. Keep proteins ≥ 10 aa.
  3. **Compare.** All-vs-all BLASTP: `-evalue 1e-5 -max_target_seqs 100000`.
  4. **Link.** For each unordered gene pair, take the HSPs whose query is the **longer** protein, greedily by
     bitscore, keeping only HSPs that do not overlap on the longer protein. The pair is linked when their summed
     length is ≥ **0.30 × the longer protein's length**. There is **no identity floor** and no bitscore floor
     beyond E ≤ 1e-5.
  5. **Cluster.** MCL on the unweighted link graph (self-loops 1, inflation 2.8, prune 1e-9, ≤ 100 iterations).
     This is the bench comparator MCL, not the shipped one.
  6. **Report.** Families are MCL clusters with **≥ 2 genes**, named `PF<i>` (MCL cluster index).

  Families are within one chromosome, so inter-chromosomal paralogues are out of scope. MCL can put two linked
  genes in different families. Fig. 6 reads a recorded gorilla table (`NC_073244.2.tsv`, 775 genes in 113
  families, built by the predecessor script); fig. 7 rebuilds each contig with `truth.py protein-referee` and
  gets the same 113 / 775 on NC_073244.2.
- **Legend.** "Protein-homology families: genes whose longest annotated proteins are similar by BLASTP (E ≤ 10⁻⁵;
  non-overlapping aligned segments covering ≥ 30% of the longer protein; no identity threshold), grouped by
  Markov clustering (inflation 2.8); pseudogenes and immunoglobulin / T-cell-receptor gene segments excluded;
  families of ≥ 2 genes on one chromosome. Built from the annotation's protein sequences, independently of both
  modes' family rule." Replace "referee", "PF82" and so on with "protein-homology family" and a gene name
  (e.g. "the RFPL4A-like array").
- **Occurs.** fig3.md:197; fig6.md:16–19,109–163,174,194,199–218; fig7.md:1,6–14,48,64,82,91,99–132,173,190–192,
  218–245; fig_family_recovery.py:66 (`TRUTH_SHORT`, printed in the row labels and panel titles 543–544);
  `_o1_recovery.py:63` (`TRUTH_LABEL`); fig_family_spectrum.py:598 (axis "referee pairs"), 676 (title
  "protein-referee pairs"), 536 (table row "referee families"); README.md:21,64,65,82;
  fig6_gorilla_clusters.tsv:27–30; fig6_gorilla_recall.tsv:27–33; `data/fig7_*.tsv:45,50`.

#### Recall universe of the seeding supplement (fig. 6s-seeding, development panel; formerly fig. 6d) (exon rule, plotted; span rule, recorded)
- **Code.** `_o1.primary_counts` (lines 141–188). **Exon rule**: a gene counts when ≥ 2 reads have their
  **primary** record (not 0x100, 0x800 or 0x4) with an aligned block on one of the gene's annotated exons. 576
  of 775 genes qualify. **Span rule**: ≥ 2 primary records overlap the gene span (`samtools view -c -F 0x904`),
  spliced-over reads included. 613 genes qualify, which reproduces the recorded list. Truth pairs are
  same-family pairs with both genes in the universe.
- **Legend.** "Genes with ≥ 2 reads whose primary alignment lies on the gene's exons."
- **Occurs.** fig6.md:16,19,113–121,129,141,193,203,225,226,254; fig_family_spectrum.py:676–677 (title);
  fig6_gorilla_recall.tsv:27–28; fig6_gorilla_clusters.tsv:28.

#### Annotated-mRNA nucleotide identity band (fig. 6s-seeding development panel, x-axis; formerly fig. 6d); "none"
- **Code.** For each pair of reference genes, the best `matches ÷ block length` (PAF columns 10/11) of any
  record of `minimap2 -x asm20 -k11 -w5 -c -X -N 100 -p 0.1 --secondary=yes` between their spliced annotated
  mRNAs (`score.py:pairs_referee_bands`, `BANDS_MRNA` ≥ 90, 80–90, 70–80, 60–70, < 60). "none" means no
  alignment record at all.
- **Legend.** "Best nucleotide identity between the two genes' annotated mRNAs (minimap2); none = the mRNAs do
  not align."
- **Occurs.** fig6.md:122–125; fig_family_spectrum.py:599 (axis).

#### Judgeable pairs
- **Code.** Compara (fig. 6c): a predicted gene pair is judged when **both genes appear in the Compara export
  with ≥ 1 paralogue row** on any chromosome (`lib.load_compara` → `genes_with_data`). Protein-homology families
  (fig. 6d): both genes are members of a family of ≥ 2 genes on the contig (`lib.read_referee`,
  `score.py:pairs_referee_bands`). Precision = same-family pairs ÷ judgeable predicted pairs. Pairs that
  involve an unlisted gene are not counted, so precision is an **upper bound**.
- **Legend.** "Precision over pairs whose two genes both belong to a reference family (pairs with a gene outside
  every reference family cannot be judged; the value is an upper bound)."
- **Occurs.** fig6.md:20,93,140,141,144,163,253; fig_family_spectrum.py:533, 617, 674 (on-figure);
  fig6_gorilla_clusters.tsv:30; fig6_gorilla_recall.tsv:33; fig6_chr16_precision.tsv:16.

#### Ensembl Compara paralogue pairs; Compara protein identity; universe of 290 pairs
- **Code.** BioMart `hsapiens_paralog_*` export, columns gene, paralog, perc_id, perc_id_r1, subtype,
  paralog_chromosome (`/mnt/linuxdisk/tmp/gw22/spectrum/compara_chr16.tsv`, 88,449 rows, no header, exported 2026-09-24
  from the `useast` mirror; release not recorded). `lib.load_compara` (lines 212–235) keeps same-chromosome pairs
  with distinct genes; a pair's identity is **max(perc_id, perc_id_r1)**. Compara genes are matched to our loci
  by gene symbol through the RefSeq CHM13 annotation. The **universe** is the 290 of 1,324 chr16 pairs whose two
  genes both have a locus in Rustle's chr16 assembly (`score.py:cmd_spectrum`). Bands are ≥ 90, 80–90, 70–80,
  60–70, 50–60, 30–50 and < 30%.
- **Legend.** "Ensembl Compara (release N) human paralogue pairs with both genes on chr16. Protein identity = the
  larger of Compara's two percent identities."
- **Occurs.** fig6.md:3,5,10,12,43,44,69,73–76,93,94,105,163,172,222,234,241,249; fig_family_spectrum.py:432–433,
  496–497, 533, 669, 672 (on-figure); fig6_chr16_recall.tsv:16; fig6_chr16_groups.tsv:16; README.md:64.

#### "Compara group" (fig. 6a,b)
- **Code.** A connected component of the universe's pair graph over all bands (fig6_chr16_groups.tsv:17). This is
  **our** grouping of Compara's pairs, not a Compara gene tree.
- **Legend.** "Paralogue group: genes connected by Compara pairs in this set."
- **Occurs.** fig6.md:10,44,69,73,76,222; fig_family_spectrum.py:434, 496–497, 672 (on-figure).

#### Locus and locus representative (fig. 6a–c); "direct edge" (spectrum tiers)
- **Code.** `score.py:cmd_spectrum` (lines 276–515). A locus is one `gene_id` group of Rustle's assembled GTF
  on chr16. Its representative is the transcript with the most `reads` (ties: longer span), kept when its spliced
  sequence is ≥ 200 bp, and mapped to the RefSeq gene with the largest exonic overlap. All-vs-all tiers run on
  the representatives with `-c -X --no-long-join -N 50 -p 0.1 --secondary=yes`:
  - T1 `minimap2 -x asm20`: identity ≥ 0.80;
  - T2 `-x asm20 -k11 -w5`: identity ≥ 0.60;
  - T3 `mmseqs easy-search --search-type 2 -e 1e-5` (translated): protein identity ≥ 0.30.

  Each tier also needs coverage ≥ 0.50. For T1/T2, identity = matches ÷ block length and coverage = query span
  ÷ shorter length (can exceed 1). For T3, coverage = max(qcov, tcov). A pair of genes has a **direct edge**
  when any locus of one and any locus of the other have a tier hit. The on-figure "nucleotide (k11 sensitive
  pass)" is T2. **The main figure uses T1 or T2 only** (`score.py spectrum --skip-t3`; the per-pair `recovered` flag is
  then nucleotide only, checked against spectrum.tsv). T3 ("translated protein search") is a **supplementary
  comparator** (fig. 6s-protein, the chr16 development spectrum); it is not part of Rustle's rule and is never run
  genome-wide.
- **Legend.** "Direct alignment: the two genes' assembled transcripts align directly (minimap2, ≥ 60% nucleotide
  identity over ≥ 50% of the shorter)." Supplement only: "+ translated protein search (≥ 30% protein identity; not
  part of Rustle's rule)."
- **Occurs.** fig6.md:10,11,12,39,49–56,66,76,81–95,238,239; fig_family_spectrum.py:499, 511–512 (on-figure);
  fig6_chr16_recall.tsv:17; fig6_chr16_groups.tsv:19,21.

#### Same default family; with / without a direct alignment; families' ceiling; "fam." (fig. 6a–c, since 2026-09-25)
- **Code.** "Same family": both genes have a copy in one default family (`score.py pairs --members
  <id>.fam.copies.tsv --universe <spectrum truth_pairs>`; copies mapped to genes by the largest span overlap). The
  split of the family bar uses the SAME pair's direct alignment (the spectrum's per-pair T1-or-T2 flag, the grey bar):
  **with a direct alignment** (blue) or **without one** (navy: the two genes are joined through other loci of the
  family, or by the family rule's genomic-span alignment where the representatives' direct alignment falls below its
  floor). The families stage writes no within-family path table, so the catalog's path-length split below no longer
  applies. **Families' ceiling** (dashed outline): both genes have a copy in some default family (the copy table holds
  family members only). "k/m" = pairs in one family ÷ pairs under the ceiling; "fam." = families holding them.
  Both bars now use the same assembly's loci (the direct alignment on spliced representatives ≥ 200 bp, the families
  on genomic spans).
- **Legend.** "Same default family: both genes have loci in one of Rustle's default families. Blue: their loci also
  align directly; navy: in one family without a direct alignment. Dashed outline: both genes are family members (the
  most the families can recover)."
- **Occurs.** fig6.md (Definitions, panel a); fig_family_spectrum.py (`_band_bars`, panel b/c labels);
  `fig6_*_recall` / `fig6_*_groups` views `families`, `families_edge`, `families_no_edge`, `families_both_members`.

#### (Legacy) Same shipped family; direct catalog edge; "only via other copies"; families' ceiling; "fam."
- **Status.** The catalog-based fig. 6a–c of 2026-09-24/25 (superseded by the entry above).
- **Code.** "Same family": both genes have a copy in one catalog family (`score.py pairs --members
  chr16.cat.copies.tsv`; copies mapped to genes by the largest span overlap). The split uses the catalog's
  `pairs.tsv` (gw_family_catalog.rs:358–368). **Direct catalog edge**: the shortest within-family E_r path
  between any copy of gene A and any copy of gene B has length 1. **Only via other copies**: the path has
  length ≥ 2. The catalog's copies are a different node set from the spectrum's loci, so a direct catalog edge
  can exist where the spectrum has none. **Families' ceiling** (dashed outline): both genes have ≥ 1 catalog
  copy. "k/m" means pairs in one family ÷ pairs under the ceiling. "fam." is the number of catalog families
  holding those pairs.
- **Legend.** "Same Rustle family: both genes have copies in one family of Rustle's copy catalog. Blue: two of
  their copies align directly. Navy: they are linked only through other copies of the family. Dashed outline:
  pairs whose two genes both have catalog copies (the most the families can recover)." Replace "shipped" with
  "Rustle's default".
- **Occurs.** fig6.md:6,39,54–71,76,81–89,99–102; fig_family_spectrum.py:436, 499–503, 514–515 (on-figure
  "Same shipped family", "only via other copies", "fam."); fig6_chr16_recall.tsv:18–19;
  fig6_chr16_groups.tsv:20–22.

#### Leave-largest-cluster-out (fig. 6c open circles; fig. 6s-seeding black ticks; formerly fig. 6d)
- **Code.** Rescore with `score.py pairs` after deleting the predicted cluster that holds the most recovered
  pairs (≥ 90% band) or judged pairs (precision) (fig6_gorilla_clusters.tsv:31).
- **Legend.** "Value without each configuration's largest predicted family."
- **Occurs.** fig6.md:131,203–208,226; fig_family_spectrum.py:537–538, 577–579 (on-figure).

#### De novo mode; guided mode; "the same family rule"
(Family modes. Not the transcript modes of figs 1–3: see *Modes* in section 1.) **The de novo mode IS the default de
novo family definition** (*Default de novo families*, section 3); the guided mode applies the same rule to the
annotation's gene bodies.
- **Code.**
  - **De novo** (`_o1_recovery.denovo_families`): the genome-wide Rustle default assembly, cut to the contig,
    then `rustle_pipeline.sh families` = `mcl_families --from-gtf --min-exonic-bp 1
    --min-shared-exon-frac 0.60`. A locus is a `gene_id` group, its span the min/max over its transcripts, its
    exons those of the transcript with the most reads. The **genomic spans** are aligned all-vs-all with
    `minimap2 -x asm20 -c -X -N 50 -p 0.1 --secondary=yes` (mcl_families.rs:34–41).
  - **Guided** (`_o1_recovery.guided_families`): the annotated gene and pseudogene bodies of the contig →
    `samtools faidx -r` → `minimap2 -x asm20 -c --eqx -P` all-vs-all →
    `mcl_families --paf --gff FULL_GFF --min-exonic-bp 1 --min-shared-exon-frac 0.60`.
  - **Family rule** (both modes; mcl_families.rs defaults 60–123, 262). An edge needs identity ≥ 0.70,
    an alignment block (PAF column 11) ≥ 300 bp, coverage ≥ 0.30 of the **longer** locus's exon-union length, ≥ 1 exonic base on both
    sides, and a best record covering ≥ 0.60 of the smaller gene's exonic length with exon-to-exon evidence.
    Clustering is MCL with inflation 2.8 and prune 1e-9. Clusters of ≥ 2 loci are reported.

  The two modes differ in node set **and** in minimap2 flags; the figure does not isolate the flags.
- **Legend.** "De novo: loci assembled from the reads. Guided: loci are the annotated gene and pseudogene bodies.
  Both are grouped by the same rule: minimap2 all-vs-all; loci linked at ≥ 70% identity over ≥ 30% of the longer
  locus, sharing ≥ 60% of the smaller locus's exonic sequence; Markov clustering (inflation 2.8)."
- **Occurs.** fig6.md:35; fig7.md:1–24,43,48,51,64–95,113,138,143,146,174–180,212–216,248–251;
  fig_family_recovery.py:530 (on-figure "same family rule in both modes"), legend "de novo (upper dot and bar)"
  and "guided (lower dot and bar)"; `data/fig7_*.tsv:49`; README.md:65,82.

#### Predicted cluster; catalog sizes (fig. 7)
- **Code.** `mcl_families` `clusters.tsv`. Each member locus is labelled with **one** gene, the one with the
  largest span overlap (`_o1_recovery.gene_at`; gene, pseudogene and ncRNA_gene records by `Name=`). Clusters are
  intersected with the reference universe before scoring.
- **Legend.** "Predicted family". Write "Numbers of predicted families and of loci in them", not "catalogue
  sizes".

#### Soto 2025 families
- **Code.** Soto et al. 2025, *Cell* 188:5363–5383, Table S1C (`bench/soto/soto_famCN_S1C.tsv`). These are
  human-specific segmental-duplication (SD98) gene families from T2T-CHM13, grouped by shared exons
  (minimap2 map-back) and family copy number. A gene can carry several Family IDs: 149 of 2,334 gene IDs (6.4%)
  do, so the table is a **cover, not a partition**. Fig. 7 keeps each gene's **first** Family ID in file order,
  skips blank and `N/A` IDs (`_o1_recovery.truth_table`), matches genes to RefSeq `Name=` values by exact name,
  and scores families with ≥ 2 genes on the chromosome. Not independent of the family rule: the 0.60
  shared-exon threshold was chosen against Soto families on chr5, chr7 and chr21.
- **Legend.** "Soto et al. 2025 human segmental-duplication gene families (Table S1C; first family per gene). Not
  independent: the family rule's exon threshold was chosen against these families on other chromosomes."
- **Occurs.** fig6.md:258–260; fig7.md:1,4,10,14,33,87,100–133,142–154,181,193–194,238–245;
  fig_family_recovery.py:66 (on-figure "Soto 2025"), 543; `_o1_recovery.py:22–27,63`; `data/fig7_*.tsv:42,50`.

#### NPIP union set ("NPIP union truth U2")
- **Code.** The table recorded at `/mnt/linuxdisk/tmp/union2/U2_truth.tsv` (37 genes, Soto family IDs). It is
  the union of Soto's NPIP families, RefSeq NPIP-named genes, and genes that align over ≥ 95% of their own length
  (length floor 10,633 bp) to a Soto NPIP member. It was built by `bench/union_truth_npip.py` (retired in
  8db314c7). 34 genes in 3 families are scored on chr16. It inherits Soto's link to the family rule and covers
  one gene family only.
- **Legend.** "NPIP reference set: Soto's NPIP families extended with RefSeq NPIP genes and genes aligning over
  ≥ 95% of their length to a Soto NPIP gene (3 families, chr16 only)." Drop "U2".
- **Occurs.** fig7.md:1,5,26,27,109,134,150,166–167,182,195–198,210,233; fig_family_recovery.py:66 (on-figure
  "NPIP union U2"); `_o1_recovery.py:26,63`; `data/fig7_*.tsv:50`.

#### Pairwise sensitivity / precision (fig. 7a)
- **Code.** `_o1_recovery.score_arm` (lines 290–357). Truth pairs are within-family gene pairs of reference
  families with ≥ 2 members on the chromosome. Predicted pairs are within-cluster gene pairs, with clusters
  intersected with the reference universe. tp = the intersection. Sensitivity = tp ÷ truth pairs; precision =
  tp ÷ predicted pairs. Large families dominate because pairs grow quadratically.
- **Legend.** "Pairwise: fraction of same-family gene pairs placed in one predicted family (sensitivity), and
  fraction of predicted same-family pairs that are same-family in the reference (precision), over genes in the
  reference families."
- **Occurs.** fig7.md:56–58,120–136,185,201,226–236; fig_family_recovery.py:520 (header "Pairwise (gene
  pairs)"); fig7_summary.tsv:53.

#### One-to-one bipartite sensitivity / precision / F (fig. 7b); `family_score` semantics
- **Code.** `target/release/family_score`, re-derived by `score_arm` (scipy `linear_sum_assignment`, asserted
  equal for all 26 rows). Families are matched one-to-one to clusters, maximising total shared members.
  Sensitivity = matched members ÷ reference members. Precision = matched members ÷ members of the matched
  clusters (clusters with overlap > 0); unmatched clusters are ignored. F is the harmonic mean of the **pooled**
  sensitivity and precision, not a mean of per-family F. Cluster members that are not in the reference universe
  are removed first, so precision is an **upper bound**.
- **Legend.** "One-to-one matching of reference families to predicted families (maximising shared genes).
  Sensitivity = matched genes ÷ reference genes; precision = matched genes ÷ genes in the matched predicted
  families (an upper bound: genes without a reference family are not scored); F = harmonic mean."
- **Occurs.** fig7.md:3,12,60–65,95–118,183,201,209–212,223,236; fig_family_recovery.py:521 (header "One-to-one
  bipartite (family members)"); fig7_summary.tsv:54; fig7_per_family.tsv:53.

#### exact / partial / missed; per-family F; flagship families (fig. 7c–f)
- **Code.** For each reference family, `hit` = members shared with its matched cluster. Per-family
  sensitivity = hit ÷ family size; precision = hit ÷ matched cluster size; F = harmonic mean. **exact**: F = 1.
  **partial**: 0 < F < 1. **missed**: hit = 0, meaning no match or a match with no shared member.
  Points in d–f are families with F > 0 in at least one mode. "Flagship" is the first of NPIP, TBC1D3, FAM90A or
  AGAP that prefixes a member name (`FLAGSHIP`). Labels are the longest common prefix of the named members.
- **Legend.** "Exact: the matched predicted family equals the reference family. Partial: they overlap. Missed:
  no member recovered."
- **Occurs.** fig7.md:14,67–72,74–93,97,146–157; fig_family_recovery.py:522 (header "Truth families, per mode"),
  532–533 (legend), 549 (on-figure "Per-family F").
- **Flag.** The panel c header "Truth families, per mode" (fig_family_recovery.py:522) should read "Reference
  families, per mode".

#### Substrate status: development / held out, untouched / held out, reused verdict set
- **Code.** `_o1_recovery.DEFAULT_STATUS` (line 69) and `STATUS_NOTE`:
  - **development** (human chr16; gorilla chr20 = NC_073244.2): rules were chosen on it. NC_073244.2 was the
    pre-registered verdict set of the seeding decision.
  - **held out, untouched** (human chr6; gorilla chr10 = NC_073234.2): no rule, threshold or verdict was ever
    scored on it.
  - **held out, reused verdict set** (human chr2, chr8, chr10): not used to develop the early family rules, but
    used since 2026-09-20 as the verdict set of about 30 pre-registered guided-mode tests. Exposure may favour
    guided.
- **Legend.** "Development: family rules were tuned here. Held out, never used: no rule, threshold or test result
  was ever computed here. Held out, reused: not used for tuning, but used to evaluate about 30 earlier
  guided-mode tests."
- **Occurs.** fig7.md:5–13,30–49,53–54,76–93,99–113,117,124–125,135,142–154,164–166,252–258;
  fig_family_recovery.py (row group labels and legend "held out, untouched", "held out, reused verdict set",
  "development"); `data/fig7_*.tsv:40–46`. See A1 for the conflicting use in figs 4–6.

#### Seeding-decision contig (fig. 6s-seeding development title; formerly fig. 6d)
- **Code.** NC_073244.2 (gorilla chr20), the contig on which the 0.98 seeding default was chosen
  (fig_family_spectrum.py:676).
- **Legend.** "Gorilla chr20 (NC_073244.2), the chromosome used to choose the seeding rule (not independent
  evidence for it)".

#### Compara families (duplications within primates) (fig. 7; the headline reference)
- **Code.** `_o1_recovery.compara_families`: nodes = genes whose Compara symbol is a RefSeq CHM13 `gene` record with
  `gene_biotype=protein_coding` on the chromosome Compara gives; a family = a connected component of the Compara
  release 116 paralogue pairs whose duplication node (`subtype`) is at or below **Primates** (Homo sapiens …
  Primates). Within one reconciled gene tree this relation is ultrametric, so the components are the gene-tree clades
  below duplication nodes of that age; no identity threshold. 426 families, 1,313 genes. Genome-wide, cross-chromosome;
  the development tables restrict it to each chromosome (a family keeps its members there, ≥ 2 to be scored).
- **Legend.** "Ensembl Compara families: protein-coding genes joined by Compara paralogue pairs whose duplication lies
  within the primates (gene-tree clades; no identity threshold)."

#### Liftoff copy pairs (family reference: figs 6s-seeding, 7; fig. 8 claim F1)
- **Code.** `_liftoff.copy_pairs` / `pair_families`. From Liftoff's self-lift of the species' annotation (`-copies -sc
  0.95`; fig. 8): pairs (a source record's annotated placement, one of its extra copies), extra copy `sequence_ID` ≥
  0.95 (rows also at 0.98 / 0.99 / 1.00), both exon unions ≥ 200 bp. Figs 6s and 7 keep the pairs whose two loci are
  both **read-supported** in the sample (≥ 2 reads whose primary alignment, `-F 2308`, has an aligned block on the exon
  union: the fig. 8 C2 rule), a universe fixed before scoring. A pair is **in one family** when a copy of the family
  covers ≥ 50% of each locus's exon union (Liftoff's `-a` on exon bases) and the two covering copies share a family.
  Fig. 8's F1 conditions on both loci being covered (`covered`). Pairwise sensitivity only: the relation certifies
  pairs and is not a partition (no precision, no bipartite matching). The guided family mode is not scored on it: its
  loci are the annotation, and an extra copy is unannotated by construction.
- **Legend.** "Liftoff copies: each annotated gene or pseudogene and the unannotated copies (≥ 95% identical) that
  Liftoff finds for it in the same genome; a pair counts when both loci carry reads and one Rustle family covers both."

#### Translated protein search (T3; supplementary comparator, fig. 6s-protein)
- **Code.** `score.py spectrum` tier T3: `mmseqs easy-search --search-type 2 -e 1e-5`, protein identity ≥ 0.30,
  max(query, target coverage) ≥ 0.50, on the chr16 development loci; never run genome-wide (hours, ≥ 13 GB).
- **Legend.** "Translated protein search (a comparator; not part of Rustle's family rule)."

#### Extra-sensitive (protein) step; family calibration I_F; missing member (supplementary figure S-P)
- **Code.** `tools/protein_attach.py` (`docs/PREREG_protein_attach_2026-09-25.md`), run by hand after the default
  families exist; never part of the default. Each locus's protein = the longest stop-to-stop frame of its
  representative's spliced exons (≥ 100 aa, < 50% soft-masked); BLASTP (E ≤ 1e-5) against the members' proteins; a
  hit counts when non-overlapping alignments cover ≥ 30% of the longer protein. **Family calibration I_F** = the
  lowest, over a family's members, of each member's best identity to another member (a family needs ≥ 2 members with
  such a hit to be calibrated). A locus the default rule left out joins the family of its best hit only if its
  identity is ≥ max(I_F, 0.60) and no other family passes (never merges families; a second passing family is
  reported as a merge proposal). A **missing member**: an unattached locus whose gene the external reference relates
  to a member's gene (sensitivity = attached ÷ missing members).
- **Legend.** "Optional protein step (run by hand; not part of the default): a locus the RNA rule left out joins a
  family only if its protein is at least as close to the family as the family's own loosest member (and ≥ 60%
  identical)."

---

## Part C. Style guide

**Words to keep out of figure text and legends.** The internal names are "arm", "referee", "truth" (outside the
simulation), "certificate" (without its definition), "shipped", "driver default", O1/O2/O3, "register row NNNN",
§-numbers, PREREG file names, E_r, U2, GWFAM ids, PF ids, T1/T2/T3, "k11 sensitive pass", "read group" and
"catalogue". Methods may keep file names and flags, but must give them a key once.

**Modes.** Every tool comparison states its mode on the figure and in the caption: figs 1–3 are annotation-free
(de novo) on every side; an annotation-guided tool run is only ever shown in its own figure, next to no
annotation-free method. Fig. 7 is Rustle-internal (two family modes), not a tool comparison: say so. Never write
"guided" without the level (transcript, family, locus).

**Method names.**
- Rustle's configurations: "Rustle" (default) and "Rustle, primary alignments only". The short form
  "Rustle (prim.)" is allowed only where it is defined in the same panel.
- Fig. 5c: "Rustle (copy-assignment step)".
- Other methods, with versions at first mention in each caption: StringTie 3.0.1, FLAIR 3.0.1, IsoSeq collapse
  26.2.0 (PacBio), minimap2 2.31, gffcompare 0.12.10, SQANTI3 5.5.4, BLASTP (BLAST+ version to
  be recorded), MMseqs2, MCL. After that, the bare name.
- Treat the other methods as baselines run by the lab, not as competitors (AGENTS.md §6).

**Samples and genomes.** Write "Human A119b (T2T-CHM13 v2.0)" and "Gorilla OR6737, testis (mGorGor1,
GCF_029281585.2)". Name gorilla chromosomes as "chr20 (NC_073244.2)" at first mention in a figure, then "chr20".
Never pool species, and say so once per figure.

**Counts.**
- In text, write "k of n" ("163 of 1,263"). Use "k/n" only inside bar labels.
- Use a thousands comma (27,447) and "n = 1,263" with spaces.
- Ranges take an en dash: "chr20–22", "0.85–1".
- Differences in proportions are in percentage points ("+17.4 points").
- Report a proportion as a percentage with one decimal when it is a share of transcripts or reads (figs 1–3), or
  as a decimal with two or three places when it is a sensitivity, precision, F, coverage or accuracy of families
  or assignments (figs 4–7). Do not mix the two in one panel.
- Give intervals as "95% CI a–b" and name the unit resampled ("resampling source copies").

**Metrics.** Say "sensitivity", never "recall". Say "precision" and name its denominator (all predicted pairs,
judgeable pairs, or assigned reads), and write "upper bound" whenever unlabelled predictions are dropped. Say
"fraction assigned" instead of fig. 5's "coverage". Write "F (pooled)" and "per-family F".

**Ties.** Use "tie" only for alignments, and always name the rule. Figs 1–3 and 6s-seeding: "a second alignment
≥ 98% of the best score". Figs 4–5a,b: "MAPQ 0". Fig. 5c: "equal best alignment score". Say "left unassigned"
for the copy-assignment status `tied`. Say "equal" or "not separable" for statistical comparisons.

**Expression and scope.** Never write "expressed" without the rule (see A1). Write "genome-wide" only when every
contig was scored; otherwise name the contigs.

**Families and references (since 2026-09-25).** Call Rustle's families "Rustle's default (de novo) families" and
state the definition once per caption (reads → seeded loci → one representative per locus → families); never call
the legacy copy catalog "shipped" or "default". Main figures use external references only (Ensembl Compara, Soto
2025 with "not independent", Liftoff copies); protein-homology families and the translated protein search appear only
in supplementary figures, labelled "secondary reference" or "comparator, not part of Rustle's rule".

**Independence and exposure.** Every figure names its unit of independence (chain, read-sharing group, source
copy, reference family) and the exposure of each substrate: development, held out and never used, or held out
and reused. Name the decision the exposure refers to.

**Cross-references.** Write "Fig. 3" (capital F, period, space) and "Figs 1–3" in text and panel notes.
