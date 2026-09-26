# Figure 2 — SQANTI3 structural categories and rules filter per method (the Figure 1 GTFs, all annotation-free)

> **Status (2026-09-25).** The build is now genome-wide for both samples and all five methods (every contig the
> sample's annotation covers, one SQANTI3 call per contig and part), plus Rustle's two configurations on the four
> other samples (supplementary figure fig2s_samples). The current tables predate that change: gorilla covers chr20,
> chr22 and chrY (NC_073244.2, NC_073246.2, NC_073248.2) and human chr20–22. **Every number below is from those
> contigs** and is replaced by the genome-wide rebuild (`make.py data fig2`, about 13–25 h of SQANTI3 in bounded
> calls); `python3 figures/fig_sqanti.py summary` then prints every number this caption quotes.

> **Mode (apples to apples).** Annotation-free (de novo) on every side: Rustle assemble (reads + genome), StringTie
> `-L` without `-G`, FLAIR collapse without annotation (`flair correct` skipped), IsoSeq collapse (checked in
> `benchmark_collapse/run_stringtie.sbatch` and `run_flair.sbatch`; the figure's first line says the same). **Guided
> comparison: not available (guided StringTie/FLAIR GTFs not supplied).** When the user registers annotation-guided
> runs (`samples.tsv` columns `stringtie_guided_gtf`, `flair_guided_gtf`), the same build scores them into a separate
> table and figure (`fig2_guided_categories` and `fig2_guided_filter`, drawn as fig2g_guided), never in a panel with the annotation-free methods. Rustle has no annotation-guided
> transcript assembly, so that figure has no Rustle row, and the guided tools are there scored against the annotation
> they were given (docs/PREREG_guided_transcript_comparison_2026-09-25.md; GLOSSARY *Modes*).

**Claim.** SQANTI3 5.5.4 classifies every method's Figure 1 transcripts against the RefSeq annotation. In both
species, Rustle's two configurations have the highest full-splice-match (FSM) share and the highest share passing
SQANTI3's rules filter of the five methods. These are point estimates from one run each: gorilla FSM 39.8% against
36.6% for StringTie, the next method; human 17.0% against 13.5% (all values in the table below). The FSM order is the
same among multi-exon transcripts only, so it is not an effect of how many single-exon transcripts a method reports.
IsoSeq collapse has the most FSM transcripts in both species (3,414 and 5,376, against Rustle's 1,931 and 2,531), and
in human FLAIR has slightly more than Rustle (2,572). **Completeness is mixed.** More of Rustle's FSM transcripts than
StringTie's match both reference ends within 50 bp (gorilla 42.1% vs 35.0%, human 26.8% vs 22.5%; in human FLAIR's
share is the highest, 37.9%). But Rustle has more 3′-fragment ISM transcripts, which lack the reference's 5′ exons:
7.2% of its transcripts against 3.7% for StringTie in gorilla, and 5.8% against 1.5% in human. Human FSM shares are
low for every method, and so is Figure 1's human intron-chain precision. On multi-exon transcripts the two count the
same thing (an exact intron-chain match), so that agreement is partly by construction. The independent evidence is
the intergenic share: 9.6–20.9% of every human method against 1.2–5.7% in gorilla. This suggests the gap lies in the
library–annotation pair, not in one method.

## Panels

There are two sample rows, never pooled. The legend at the top names the five folded categories.
- **a–c**: gorilla OR6737, testis (mGorGor1, GCF_029281585.2; RefSeq GCF_029281585.2 annotation). Current table:
  chr20 (NC_073244.2), chr22 (NC_073246.2) and chrY (NC_073248.2), named in the row header; 7,630 reference
  transcripts. After the rebuild: genome-wide.
- **d–f**: human A119b (T2T-CHM13 v2.0; CHM13 v2.0 RefSeq annotation). Current table: chr20–22, 10,851 reference
  transcripts. After the rebuild: the 24 contigs the annotation covers (chrM has no annotation and is left out for
  every method, as in Figure 1).

**Methods.** These are the Figure 1 GTFs, cut to the same contigs:
- Rustle: loci built from the primary alignments plus the secondary alignments whose alignment score is at least 98%
  of the read's best;
- Rustle, primary alignments only (`--no-seed-secondaries`), drawn hatched in b, c, e and f;
- StringTie 3.0.1, FLAIR 3.0.1 and IsoSeq collapse 26.2.0: baselines the lab ran on the same alignments.

**a, d** SQANTI3 structural categories, as a 100% stacked bar per method. Every segment is named with its category
and share, inside the segment or, when it is too narrow, in a call-out above the bar.
- FSM (full splice match): every splice junction matches a reference transcript with the same number of junctions.
  ISM (incomplete splice match): consecutive reference junctions, but fewer 5′ or 3′ exons. NIC (novel in catalog):
  a new combination of annotated splice sites or junctions. NNC (novel not in catalog): at least one unannotated
  splice site.
- **Other** holds SQANTI3's antisense, intergenic, genic (overlapping exons and introns), genic intron (entirely
  inside an annotated intron) and fusion transcripts. No method has a moreJunctions transcript.
  - Other is 4.1–17.5% of each gorilla method and 22.6–50.3% of each human method.
  - Intergenic transcripts are its largest part in every human method except FLAIR, where genic intron (18.3%) is
    slightly above intergenic (17.5%). Rustle: 14.6% intergenic, 5.8% genic intron, 3.7% antisense, 1.7% genic,
    1.6% fusion.
- The n under each method is the number of transcripts SQANTI3 classified on these contigs.

**b, e** The share of the method's transcripts that pass SQANTI3's default rules filter (`filter_result` =
"Isoform"): at most 59% A in the 20 bp downstream of the 3′ end; for non-FSM transcripts also all-canonical junctions
and no reverse-transcriptase template-switching signature. It is not a correctness measure.

**c, f** The number of FSM transcripts.

| current tables | Rustle | Rustle, primary alignments only | StringTie | FLAIR | IsoSeq collapse |
|---|---|---|---|---|---|
| gorilla n | 4,850 | 4,727 | 4,647 | 7,998 | 28,216 |
| gorilla FSM share (a) | 39.8% | 39.9% | 36.6% | 21.8% | 12.1% |
| gorilla PASS (b) | 79.3% | 79.5% | 75.0% | 59.8% | 49.4% |
| gorilla FSM transcripts (c) | 1,931 | 1,885 | 1,703 | 1,747 | 3,414 |
| human n | 14,906 | 13,386 | 14,164 | 54,765 | 169,422 |
| human FSM share (d) | 17.0% | 17.5% | 13.5% | 4.7% | 3.2% |
| human PASS (e) | 59.9% | 63.5% | 48.3% | 29.0% | 27.8% |
| human FSM transcripts (f) | 2,531 | 2,347 | 1,910 | 2,572 | 5,376 |

Not drawn; from `fig2_sqanti_filter` (multi-exon columns) and `fig2_sqanti_subcategories`:

| current tables | Rustle | Rustle, primary alignments only | StringTie | FLAIR | IsoSeq collapse |
|---|---|---|---|---|---|
| gorilla FSM share, multi-exon transcripts only | 39.9% | 39.9% | 36.8% | 22.6% | 13.5% |
| gorilla FSM with both ends within 50 bp (share of FSM) | 42.1% | 41.9% | 35.0% | 30.1% | 16.0% |
| gorilla ISM, all (share of transcripts) | 9.9% | 10.2% | 7.1% | 3.5% | 12.4% |
| gorilla ISM 3′ fragment (share of transcripts) | 7.2% | 7.4% | 3.7% | 1.7% | 5.5% |
| human FSM share, multi-exon transcripts only | 17.3% | 18.7% | 14.8% | 6.4% | 3.8% |
| human FSM with both ends within 50 bp (share of FSM) | 26.8% | 28.9% | 22.5% | 37.9% | 18.8% |
| human ISM, all (share of transcripts) | 9.9% | 10.8% | 5.7% | 2.2% | 7.9% |
| human ISM 3′ fragment (share of transcripts) | 5.8% | 6.5% | 1.5% | 0.6% | 2.7% |

## Supplementary figure fig2s_samples: Rustle on all six samples

The same three panels (structural categories, rules-filter PASS, FSM transcripts) for Rustle's two configurations on
every sample of the registry, each against its own RefSeq annotation, genome-wide: human A119b, human testis
(ENA ERR13885926), gorilla OR6737 (testis), gorilla KB3781 (fibroblast cell line), chimpanzee PTR (mPanTro3,
GCF_028858775.2) and orangutan PPY (mPonPyg2, GCF_028885625.2). Two bars per sample (Rustle, then hatched: primary
alignments only). Chimpanzee and orangutan use their RefSeq GFF3 converted to GTF by `gff_to_gtf`, the converter of
the human and gorilla GTFs (it reproduces the gorilla GTF byte for byte). This is Rustle only; there are no lab
baselines for these samples. Tables: `fig2_samples_categories`, `fig2_samples_filter`.

## Methods

`python3 figures/make.py data fig2` runs `fig_sqanti.build`. It writes `fig2_sqanti_categories` (the counts per
structural category), `fig2_sqanti_filter` (n, PASS and FSM per method, all transcripts and multi-exon ones),
`fig2_sqanti_subcategories` (SQANTI3's `subcategory` counts per structural category) and the two supplementary
tables.
- **Reference.** The Figure 1 annotation GTFs, split per contig: gorilla `GGO.GCF_029281585.2_RefSeq.gtf.gz`; human
  `A119b.chm13v2.0_RefSeq.gtf.gz`; chimpanzee and orangutan converted from `PTR_genomic.gff` and `PPY_genomic.gff`.
  Gene-level RefSeq records with an empty `gene_id` (480 gorilla, 287 human on the current contigs) take their
  transcript id as `gene_id`. Without this, SQANTI3 would merge all of them into one gene.
- **Per sample × method × contig:**
  - `sqanti3_qc.py --isoforms <method>.gtf --refGTF ref.gtf --refFasta genome.fa --report skip -t N`
  - `sqanti3_filter.py rules --sqanti_class <method>_classification.txt --filter_gtf <method>_corrected.gtf
    --skip_report`
  - Each output directory is emptied before a run. A per-contig cut is replaced only when its content changes, so an
    unchanged contig is never re-classified.
- **Large methods in parts.** A method with more than 25,000 transcripts on a contig is run in parts of at most
  25,000 whole transcripts (`fig2_chunk_tx`).
- **Parts and contigs sum exactly.** A structural category is decided per transcript against its own contig's
  annotation. Every default filter rule is also per transcript (`perc_A_downstream_TTS`, `all_canonical`,
  `RTS_stage`, `min_cov`). So the per-contig and per-part tallies sum to one run over all the contigs.
- **Subcategories.** An FSM is `reference_match` when both its ends are within 50 bp of the matched reference
  transcript's; otherwise it has alternative 5′ and/or 3′ ends. An ISM `3prime_fragment` shares the reference's last
  junction but not its first, so it lacks the reference's 5′ exons. Multi-exon means more than one exon (SQANTI3's
  `exons` column).
- **Cost.** 36 recorded calls (170–21,462 transcripts) took 6–274 s each. Genome-wide, about 6.2 M transcripts remain
  in 573 calls, estimated 13–25 h in total (Rustle on all six samples 3–5 h; the gorilla baselines 1.5–3 h; the human
  baselines 8.5–17 h, IsoSeq collapse alone 6–12 h); `figs_plan=1` prints the list, `fig2_part=samples` builds the
  supplementary tables first, and `figs_budget_s` bounds each call.

## Caveats

- **The gorilla rows of the current table cover three chromosomes, not Figure 1's genome-wide gorilla scope.** They
  are chr20, chr22 and chrY, a subset kept small when SQANTI3 was run once. chrY carries testis-expressed multi-copy
  gene families. The rebuild removes this caveat.
- **FSM shares and counts depend on the annotation.** NIC and NNC transcripts are novel junction combinations, either
  unannotated transcripts or assembly errors. SQANTI3 cannot tell them apart.
- **FSM share and Figure 1's precision are not independent.** On multi-exon transcripts both count exact
  intron-chain matches to the annotation. Human multi-exon FSM shares are 17.3, 18.7, 14.8, 6.4 and 3.8% (in the
  method order of the table), against Figure 1 precisions of 17.3, 18.7, 14.8, 6.3 and 2.4%.
- **Genic-intron transcripts are ambiguous.** They are 5.8–18.3% of every human method. They can be unannotated
  transcripts, but also pre-mRNA or genomic priming, so only the intergenic share argues for missing annotation.
- **A 3′ fragment is not always an incomplete assembly.** It can also be an alternative first exon downstream or a
  5′-truncated cDNA. SQANTI3 does not tell these apart.
- **FSM counts are per transcript, not per reference transcript.** Several transcripts of one method can match the
  same reference, which favours the methods that report more transcripts. IsoSeq collapse reports 28,216 (gorilla)
  and 169,422 (human) transcripts on the current contigs, against Rustle's 4,850 and 14,906.
- **PASS is not a correctness measure.** The default rules keep an FSM with at most 59% A in the 20 bp downstream of
  its 3′ end. They keep any other transcript that passes the same test, is not a reverse-transcriptase
  template-switching candidate and has all-canonical junctions (the alternative, short-read junction coverage ≥ 3,
  is not available here).
- **Adding the secondary alignments within 2% of the best score lowers both shares in human.** It adds 1,520
  transcripts: 899 intergenic and 184 FSM. FSM share falls from 17.5% to 17.0% and PASS from 63.5% to 59.9%. In
  gorilla it adds 123 transcripts, 46 of them FSM.
- **SQANTI3 does not classify transcripts without a strand.** So StringTie's n here is 14,164 of the 14,169 human
  transcripts in Figure 1, and 4,647 of 4,648 in gorilla. No other method has transcripts without a strand.
- **The human reference count differs from Figure 1.** gffcompare's parse in Figure 1 loads 10,836 of the 10,851
  reference transcripts from the same restricted file.
