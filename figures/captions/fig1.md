# Figure 1 — Intron-chain sensitivity and precision against the annotation: Rustle, StringTie, FLAIR and IsoSeq collapse on the same alignments, all annotation-free (de novo)

> **Status (2026-09-25).** The build is now genome-wide for every sample: every contig the sample's annotation
> covers (`assembly.evaluation_contigs`). The current tables predate that change: gorilla is already genome-wide,
> human A119b is still chr20–22. The human numbers below are therefore **chr20–22 numbers** and are replaced by the
> genome-wide rebuild (`make.py data fig1`), after which `python3 figures/fig_intron_chain.py summary` prints every
> number this caption quotes. The supplementary figure (fig1s_samples) appears with that rebuild.

> **Mode (apples to apples).** Annotation-free (de novo) on every side: Rustle assemble (reads + genome), StringTie
> `-L` without `-G`, FLAIR collapse without annotation (`flair correct` skipped), IsoSeq collapse (checked in
> `benchmark_collapse/run_stringtie.sbatch` and `run_flair.sbatch`; the figure's first line says the same). **Guided
> comparison: not available (guided StringTie/FLAIR GTFs not supplied).** When the user registers annotation-guided
> runs (`samples.tsv` columns `stringtie_guided_gtf`, `flair_guided_gtf`), the same build scores them into a separate
> table and figure (`fig1_guided`, drawn as fig1g_guided), never in a panel with the annotation-free methods. Rustle has no annotation-guided
> transcript assembly, so that figure has no Rustle row, and the guided tools are there scored against the annotation
> they were given (docs/PREREG_guided_transcript_comparison_2026-09-25.md; GLOSSARY *Modes*).

**Claim.** Every method's transcripts are built from the same alignments, restricted to the same contigs and scored
by gffcompare 0.12.10 against the RefSeq annotation of those contigs. Rustle's intron-chain precision is slightly
above StringTie's in gorilla (35.3% vs 34.6%) and above it in human (17.3% vs 14.8%), and 2.0–7.2 times that of
FLAIR (17.6%, 6.3%) and IsoSeq collapse (6.7%, 2.4%). Its sensitivity is above StringTie's in both species, above
FLAIR's in gorilla (27.1% vs 25.5%) and equal to it in human (23.0% vs 23.1%). IsoSeq collapse has the highest
sensitivity over the whole annotation (30.8%, 27.7%), with 7.4 (gorilla) and 11.4 (human) times as many transcripts as
Rustle. On the reference intron chains that at least 2 primary reads carry exactly (gorilla 27,447 chains, human
2,677), Rustle's two configurations have the highest sensitivity point estimates (93.4%, 84.4%; the two
configurations are equal in human and 3 chains apart in gorilla). This stratum is Rustle's minimum read support: it
does not depend on any method's output, but it is matched to Rustle's design. Paired on those chains, Rustle's
sensitivity is higher than that of StringTie, FLAIR and IsoSeq collapse in gorilla (+17.4, +11.9 and +10.9 points)
and than that of StringTie and FLAIR in human (+20.5, +8.6). In human it is **not separable** from IsoSeq collapse's
(+1.5 points, Tango 95% CI −0.4 to +3.5, exact McNemar p = 0.13). On chains carried by at least 1 read, IsoSeq
collapse has the higher sensitivity in both species (−10.1 and −13.6 points). Adding the secondary alignments that
score within 2% of the read's best alignment score adds 162 (gorilla) and 4 (human) matched chains, at a small
precision cost (35.6% → 35.3%, 18.7% → 17.3%).

## Samples and scope

- **Gorilla OR6737, testis** (mGorGor1, GCF_029281585.2): IsoSeq reads aligned with minimap2 2.31 (`GGO_mm.bam`);
  RefSeq annotation GCF_029281585.2. Scope: genome-wide (the annotation covers all 26 contigs).
- **Human A119b** (T2T-CHM13 v2.0; the tissue is not recorded): IsoSeq reads aligned with minimap2 2.31
  (`A119b.t2t.bam`); CHM13 v2.0 RefSeq annotation. Genome-wide scope = the 24 contigs the annotation covers. The
  annotation has no record on chrM, so chrM is left out for every method alike (chrM transcripts: Rustle 39 in each
  configuration, StringTie 9, FLAIR 47, IsoSeq collapse 1,509; also in the table notes). Current tables: chr20, chr21 and chr22.
- Species and samples are never pooled.

## Panels

Rows 1 and 2 are the two samples; row 3 is the paired comparison, one panel per sample.
- **a–c, g**: gorilla OR6737.
- **d–f, h**: human A119b.

**Methods.**
- **Rustle** (filled circle) is the pipeline's default. It builds loci from the primary alignments plus the
  secondary alignments whose alignment score is at least 98% of the read's best alignment score anywhere in the
  genome.
- **Rustle, primary alignments only** ("Rustle (prim.)" in the panels) is the same run with
  `--no-seed-secondaries`. It has the same hue, drawn as an open circle, a dashed line or a hatched bar; the open
  marker and the hatch mean this configuration only.
- **StringTie** 3.0.1 (square), **FLAIR** 3.0.1 (triangle) and **IsoSeq collapse** 26.2.0 (`isoseq collapse`,
  default settings; diamond) are baselines the lab ran on the same alignments, genome-wide.

**a, d** Intron-chain precision against sensitivity, one marker per method. Sensitivity = reference intron chains
matched ÷ multi-exon reference transcripts; precision = matching transcripts ÷ multi-exon transcripts of the method.
- The reference denominator is 95,575 multi-exon transcripts (gorilla) or 9,834 (human), given in the panel title.
- In a, the two Rustle markers are 0.2 points apart in sensitivity and 0.3 in precision ("Δ sens 0.2 / prec 0.3
  pts"). A framed key names both.

**b, e** Reference intron chains matched exactly (gffcompare "Matching intron chains"). The n under each method is
its number of transcripts (all of them, not only multi-exon ones).

**c, f** Intron-chain sensitivity stratified by the number of primary reads (`samtools -F 2308`) whose intron chain
equals the reference chain exactly.
- The unit is a distinct multi-exon reference chain (contig + chain). The strata are all chains, ≥ 1, ≥ 2 and ≥ 5
  reads.
- The chains per stratum are under the ticks: gorilla 95,810, 36,408, 27,447 and 18,927 (shown as 95.8k and so on);
  human 9,833, 3,384, 2,677 and 1,854.
- A chain counts as matched when gffcompare gives class code `=` to any reference transcript that carries it
  (exact intron-chain match: the same introns with identical coordinates in the same order; ends may differ).
- The shaded ≥ 2 column is Rustle's minimum read support ("Rustle's minimum").
- The two Rustle configurations coincide at almost every stratum (at most 0.2 points apart; identical at ≥ 2 and ≥ 5
  in human and at ≥ 5 in gorilla). So the primary-alignments-only configuration is drawn as a larger open ring around
  the filled Rustle marker, and where their last points are less than 0.5 points apart one label, "Rustle (both)",
  names both.
- The methods are shifted slightly left and right of each stratum so that coinciding markers stay visible (for
  example StringTie 83.6% and IsoSeq collapse 83.3% at ≥ 5 in gorilla).

**g, h** Paired comparison on the ≥ 2-read stratum, with the same reference chains for every method.
- Left: Rustle's sensitivity minus the other method's, in percentage points, with Tango's asymptotic score 95%
  interval. In g the intervals are narrower than the markers (each at most 1.1 points wide; see the table).
- Right: the chains matched by both, by Rustle only and by the other method only, and the exact two-sided McNemar p.
- p values below 1e-300 are printed as "<1e-300"; the exact values are in the table below.

| ≥ 2 reads (count mode `reads`) | vs | both | Rustle only | other only | Δ (points) | 95% CI | McNemar p |
|---|---|---|---|---|---|---|---|
| gorilla, genome-wide, 27,447 chains | Rustle (prim.) | 25,630 | 6 | 3 | +0.01 | −0.01 to +0.04 | 0.51 |
| | StringTie | 20,440 | 5,196 | 414 | +17.4 | +16.9 to +17.9 | 3.1e-1049 |
| | FLAIR | 21,519 | 4,117 | 861 | +11.9 | +11.4 to +12.3 | 5.3e-505 |
| | IsoSeq collapse | 21,162 | 4,474 | 1,485 | +10.9 | +10.4 to +11.4 | 5.6e-343 |
| human, chr20–22 (current table), 2,677 chains | Rustle (prim.) | 2,259 | 0 | 0 | 0.0 | −0.1 to +0.1 | 1 |
| | StringTie | 1,617 | 642 | 92 | +20.5 | +18.7 to +22.4 | 2.4e-102 |
| | FLAIR | 1,824 | 435 | 206 | +8.6 | +6.7 to +10.4 | 8.7e-20 |
| | IsoSeq collapse | 1,888 | 371 | 330 | +1.5 | −0.4 to +3.5 | 0.13 |

## Supplementary figure fig1s_samples: Rustle on all six samples

Rustle's two configurations on every sample of the registry (`figures/samples.tsv`), each scored genome-wide
against its own RefSeq annotation, with the same rules as panels a, d and the ≥ 2-read stratum of c, f. There are no
lab baselines for these samples, so this is Rustle only and is not a comparison between methods.
- Samples: human A119b and human testis (public library, ENA ERR13885926; T2T-CHM13 v2.0); gorilla OR6737 (testis)
  and gorilla KB3781 (fibroblast cell line of the individual the mGorGor1 assembly was built from; mGorGor1,
  GCF_029281585.2); chimpanzee PTR (mPanTro3, GCF_028858775.2) and orangutan PPY (mPonPyg2, GCF_028885625.2). The
  tissue is not recorded for human A119b, chimpanzee and orangutan.
- Annotations: chimpanzee and orangutan have a RefSeq GFF3 only. It is converted to GTF by the converter that built
  the human and gorilla GTFs (`gff_to_gtf`, one contig at a time; every transcript type kept). Run on the gorilla
  GFF3, the same procedure reproduces the gorilla GTF byte for byte (md5 17052eec).
- **a** Intron-chain sensitivity over every multi-exon reference transcript. **b** Intron-chain precision. **c**
  Sensitivity on the reference chains carried exactly by ≥ 2 primary reads. One row per sample, grouped by species;
  filled = Rustle, open ring = primary alignments only; printed = Rustle.
- The human testis library was aligned without minimap2's `-uf` (strand guessed, not forced) and with minimap2 2.30.
- Table: `fig1_samples`. Numbers: `python3 figures/fig_intron_chain.py summary`.

## Methods

`python3 figures/make.py data fig1` runs `fig_intron_chain.build`, which writes `fig1_gffcompare`, `fig1_support`,
`fig1_paired` and `fig1_samples` to `figures/data/`.
1. **Assemblies** (`figures/assembly.py`, run cache `figures/samples.py`). Rustle is run genome-wide through the
   pipeline driver: `bash tools/rustle_pipeline.sh assemble --bam B --fasta G --out <sample> --bin BIN --threads N
   [--no-seed-secondaries]` (`make.py runs --stage assemble|assemble_primary`). The figure build never re-assembles:
   it stops and names the `make.py runs` command when an assembly is missing or out of date. Every method's GTF and
   the annotation are then restricted to the evaluation contigs.
2. **Scoring.** Each method is scored with `gffcompare -r ref.gtf` (0.12.10). A reference chain is matched by a
   method when the method's `.tmap` gives class code `=`.
3. **Exact-chain support.** One pysam pass is made over each sample's BAM on the evaluation contigs, one contig at a
   time (cached per contig; the parts add up to the same counts as one pass, checked on two gorilla contigs).
   - Only primary alignments are used (`-F 2308`); every CIGAR `N` is an intron.
   - A read is counted only when its chain equals a reference chain. The read's strand is not used.
   - Two count modes are kept. `reads` counts alignments. `distinct_ends` counts distinct (start, end) spans, which
     is the assembler's coordinate-duplicate key. The panels plot `reads`.
4. **Paired statistics.** The Tango interval and the exact McNemar p are computed without scipy. On every row of
   `fig1_paired`, the McNemar p agrees with `scipy.stats.binomtest` and the Tango interval is within 0.07 points of
   Newcombe's hybrid score interval.
5. **Cost** (measured 2026-09-25): genome-wide gffcompare of the 3.2 M-transcript human IsoSeq collapse GTF, 107 s /
   5.1 GB. `figs_budget_s` splits a build into calls of bounded length; `figs_plan=1` lists what a build would run.

## Caveats

- **The annotation-wide sensitivity counts every annotated transcript**, including those that no read carries. This
  applies to a, d and the "all" points of c and f.
- **The ≥ 2 stratum is not neutral.** Rustle reports a chain only when at least 2 reads carry it (its own count, after
  junction correction, not the exact-chain primary reads counted here). So Rustle also matches 243 gorilla reference
  chains that fewer than 2 exact-chain primary reads carry (192 with none, 51 with one).
  - On chains carried by ≥ 1 read, IsoSeq collapse has the higher sensitivity: gorilla −10.1 points (−10.7 to −9.5),
    human −13.6 (−15.5 to −11.6).
  - On chains carried by ≥ 1 read, human Rustle and FLAIR are not separable: −0.2 points (−1.9 to +1.5), p = 0.84.
- **IsoSeq collapse's higher sensitivity at ≥ 1 read comes from single-read chains.** Take the ≥ 1-read chains that
  IsoSeq collapse matches and Rustle does not. Exactly one read carries 6,663 of 8,148 of them in gorilla and 501 of
  831 in human. At ≥ 5 reads, Rustle's sensitivity is higher than IsoSeq collapse's in both species (+14.4 and +8.6
  points).
- **The count mode matters for human vs IsoSeq collapse.** The `distinct_ends` mode is Rustle's minimum after its
  coordinate-duplicate filter. In that mode the human ≥ 2 stratum has 2,596 chains, and Rustle's sensitivity is
  higher: 87.0% vs 82.8%, +4.2 points (+2.3 to +6.1), p = 1.7e-5.
- **Chains are treated as independent, and p values are unadjusted** for the four comparisons per panel. Chains of
  one gene share reads, so the true uncertainty is somewhat larger. With discordant counts such as 5,196 vs 414, this
  changes no conclusion.
- **The two Rustle configurations barely differ on these strata.** At ≥ 2 reads, gorilla has 6 vs 3 discordant
  chains (p = 0.51) and human has none. At ≥ 1 read, gorilla has 57 vs 3 (+0.15 points, p = 6.3e-14).
- **Two parses of the annotation give slightly different counts.** Our parse gives the chain denominators of c and f
  (95,810 and 9,833). gffcompare's parse gives 95,575 and 9,834 multi-exon reference transcripts. The "all"
  sensitivities agree with a and d to within 0.1 point.
