# Pre-registration: a read-end readthrough junction filter before locus formation (better locus representatives)

**Written 2026-09-25, before the filter exists in the assembler and before any number of the kind decided here
exists.** User request: *"improve the locus representatives to try to avoid some readthroughs and account for more
TES/TSS; the info is in `docs/READTHROUGH_G50K_AND_LAST_EXON_2026-09-25.md`; let's just ensure it is beneficial."*
This file fixes the intervention, the metrics, the substrates and a multi-part bar for "beneficial". A default flip
is the user's call whatever the outcome.

## 0. What existed and what was seen before this file

- **Seen (junction level only):** every number in `docs/READTHROUGH_G50K_AND_LAST_EXON_2026-09-25.md` (the "doc"),
  i.e. the R, R_peak, R & Q, R & Q1 and BAM-tag results on human A119b chr16 (dev) and chr20, and gorilla OR6737
  NC_073244.2. Code read: `bench/mechanism/readthrough_rules.py`, `copy_assign --assemble-only` (pass-1 skeletons,
  `collapse_loci_groups`), `family_detect::{collapse_parent, pick_locus_rep}`, `mcl_families::{gtf_loci,
  loci_from_gtf}`, `denovo_pipeline::{retain_non_readthrough, is_unspliced_readthrough, is_giant_intron_mischain}`,
  `src/bin/readthrough_filter.rs`, `tools/rustle_pipeline.sh` (`assemble`, `families`).
- **Not seen / does not exist:** no implementation of R or Q1 in `src/` (only the Python instrument); no assembly,
  family table, gffcompare, TES/TSS, fused-locus or Liftoff number for any arm below; no families-stage output for any
  sample (`${work}/runs/*/*.fam.*`: 0 files); no Liftoff self-lift merged for any species (checked by listing).
  The genome-wide BASE assemblies (`${work}/runs/<id>/<id>.gtf`, 2026-09-25 ~13:10) exist; I read their file lists,
  run times, and one unrelated count (`low_confidence "true"` tags in the gorilla_OR6737 GTF: 844 of 74,045
  transcripts, the §6hs boundary tag). Nothing here was computed on them.
- **Where the representative sits today.** A de novo locus = one `gene_id` of the assembly GTF (transcripts joined by
  a shared junction, then span-aware merging). The families stage (`mcl_families --from-gtf`) aligns the locus's
  **genomic span** (union of all its transcripts) all-vs-all and uses the **representative** (most-read transcript,
  ties to the longer span) for exonic bases and the copy table. A readthrough chain in the `gene_id` therefore (i)
  fuses gene A and gene B into one locus, (ii) stretches the span over B, and (iii) leaves B without a representative
  of its own (register 1015: "a locus carries ONE label"). No default filter removes SPLICED readthroughs today.

## 1. Audit: every prior result on representatives, readthrough loci and TSS/TES

Implication column = what it means for this intervention and how this one differs.

| source | what it measured | implication here |
|---|---|---|
| doc §1-§3, r749 (09-07) | `-G 50k` keeps 100% of ≤50 kb readthrough read-junctions, breaks 99.6% of real ≥50 kb introns (r749: 1 of 46 RT units > 50 kb; 8.2% of gorilla transcripts damaged) | no aligner lever; the filter works on read ENDS, never on `-G`/`-r` (r752) |
| doc §3/§8/§14 | R (θ 0.9524, dev chr16) held-out AUC .960/.964; canonical R recall .520/.420, FPR .0040/.0012; R & Q1 (θ 0.625) recall .358/.250, FPR .0011/.0001, precision .86/.96 | thresholds frozen from here. ⚠ Scored ONLY on junctions labelled RT or ANN: novel junctions (neither) were never scored, but an assembler filter will see them → flag census by class (§5) and guards G2/G4 |
| doc §6, §10, §12 | 3′-end density, 9 BAM tags, tag combinations: refuted; read count alone AUC .85 but a floor catching half the RT costs 12-14% of real junctions | no quality/depth gate is added; R is a RATIO at the donor, not a count floor |
| r808 (§6hr) | per-transcript relative depth flags 75% of family transcripts; "the real signal is JUNCTION-level: crossing reads vs reads at the same donor that terminate before it" | R is exactly that signal, never built into the assembler |
| §6hs `--min-boundary-fraction 0.10` (default) | isolated far-boundary outliers tagged `low_confidence`, 0.5-1% of transcripts | TAG only; `gtf_loci` reads every transcript, so it never reaches the families. The filter here changes the read set, hence the loci |
| R4 `is_unspliced_readthrough` / `retain_non_readthrough` (07-09, default ON on family paths) | single-exon transcripts engulfing ≥ 5 distinct supported junctions = unspliced pre-mRNA; 4 rival rules failed | different class (unspliced). Spliced readthroughs pass R4 by construction (`introns` non-empty) |
| `is_giant_intron_mischain` | giant intron carried by < gate reads | different class (sub-gate mis-chain); RT introns are ordinary (median 18 kb) |
| `readthrough_filter.rs` / §6n8-§6o0, r852 | flags ANNOTATED readthrough records whose fusion reads are mostly SECONDARY (default no-op); fusion reads are primary (supplementary 0.1%), 97.2% MAPQ 60, `de` .0020 | readthroughs are real resident molecules: nothing to fix by alignment quality (doc §10 agrees). This filter is de novo, primary-read, junction-level |
| r750-r754, r761, doc §8 | readthroughs canonical (97-100%), recur across molecules, 32/42 annotated, `flair correct` has nothing to fix | R keeps dominant fusions (PKD1P6-NPIPP1, SLX1B-SULT1A4) by design; it only removes minor run-on past a strong polyA site |
| r845 (§6m2) chimera policy | truth-side lever; every gain was denominator shrinkage (CP-5 guard) | fixed universes and denominator guards below (G2, G4) |
| r846 (§6m3) node cut; r967/r968 (§6w0) read splits; r1001 (§6x3) | cut DOUBLES rather than separates, short pieces become hubs; chimeric-bridge split +0.00 pp, turnover +0.78 pp (bar 2.0); a binary split gives ≤ 2× shrink where 4.83× is needed | this filter cuts nothing at a coordinate and splits no read: reads carrying a flagged junction are REMOVED before pass-1. Same residual risk (pieces fail the coverage gate, r979/r1018) → G5 |
| r815 (§6ji) `RUSTLE_LOCUS_BRIDGE_CUT` | cutting spliced bridge transcripts (topology: junction touching ≥ 2 better-supported groups) won dev, failed gorilla hold-out (R_G .0734 vs .0736) | a local win did not reach family metrics → family clause is NON-INFERIORITY, not efficacy. Predicate differs (read-end geometry, validated per junction against annotation) |
| r828 (§6ke) | splitting chained read groups: human hold-out R_G lower | same lesson: hold-out at the family level decides |
| H:453 | read min-cut between fused loci median 12, max 273: fusion is REDUNDANT | a locus separates only when EVERY bridging chain goes ⇒ locus-level effect < junction-level recall (sets X, §4) |
| H:274 | read-through filter on blind DNA nodes removed 0 of 6,824 links | effect must be shown on the loci actually built, not assumed |
| r937-r940 (§6u7) | defect is OVER-merge (172/369 collisions = disjoint genes fused by readthrough); fusions carry median 38 MAPQ-60 reads (≥3-read floor leaves 84.9%); RefSeq label reaches 4.7%; junction-disjointness between loci is a tautology | fused loci are defined by ANNOTATION overlap (§3a), never by junction sharing; no support floor, no RefSeq label |
| r951/r953/r954 (§6v1/§6v2) | readthrough alone explains the over-merge (sim); a readthrough-free ideal gives family F +0.0395 Soto / +0.0139 referee on chr16; guided still ahead through node COVERAGE | ceiling for removing ALL readthrough reads; the filter removes 24-40% of RT reads (doc §3) ⇒ expected family ΔF ≤ ~+0.01 |
| r971, r1012, r1013, r1015-r1018 (§6w2, §6x6, §6x7) | fused loci cost 28/306 referee genes (NPIPB4/B5/B12 lost to LOC/SMG1P loci); FP and FN of nodes independent; bimodal (49.6% trim-type, 34.2% comparable genes); **oracle perfect split of 233 chr16 fused loci: referee F 0.214 → 0.235 (+0.021), NPIP Soto F 0.727 → 0.687**; pieces fail the coverage gate | family F gain is not expected and not required; NPIP is reported and caps the verdict (§6) |
| r303, r1039-r1041, r1050-r1055 (`RUSTLE_LOCUS_EXON_UNION`, union/corroborated/ORF reps, selection) | rep completeness median 0.627; union widening EVICTS (−47), ORF denominator makes HUBS; "the shipped single-chain rep sits at the usable point"; selection never scored at family level | the representative RULE is unchanged (most-read transcript). Only the read set feeding the loci changes; no denominator is widened or shrunk by fiat |
| §6w3 r973-r975, r978, r979 (§6w6) | 5′ error is DISPERSION (median \|d5\| 229 bp human, 136 gorilla; 3′ 22/5 bp); bias per library (+1 vs +67 bp); k = 2 terminal support is a local optimum; quantile boundary worse at every q (shrinking a locus costs its own alignments) | no end is moved or extended. TSS judged at ±250 bp, TES at ±25 bp; TES is primary |
| H:396, H:397, H:543, M:492, M:912, r806 | trimming to a read core refuted; extending loci to the TES: 0 verdict changes; TES differs in 33% of sibling pairs; `tss_tts.rs` dead; `RUSTLE_TSS_SNAP` opt-in (n.s.); exon-sum 3′ exact < 30 bp on all tools | TES carries real signal and is precise; no snapping/extension arm is added |
| r1067-r1073, r1072, r1078/r1079 (§6z9/§6za) | single-read floor stays 2; 5′ fold-in and terminal trims manufacture chains no read carries | reads with a flagged junction are dropped whole, never trimmed/split at the junction (a piece would carry a 3′ end at the donor that no molecule has) |
| r1075 (adopted) | retained-intron drop by a junction-support RATIO (10) passes held-out | ratio rules on junction support can work; different target (retained introns) |
| r760 / feedback "hold a substrate back" | 3 dev substrates agreed, the held-out lost NPIPA1/A6 → reverted | held-out samples decide; dev never enters the verdict |
| §6m1 (no-readthrough counterfactual) | deleting all 209 readthrough records does not certify NPIP; the boundary is internal (PKD1P6-NPIPP1) | no family rescue is claimed for NPIP |
| metric traps | "never judge a node change on node-level metrics"; denominators conditioned on the prediction; deleting evaluation lever → denominator guard; shrunk nodes → hubs; bipartite precision not tie-invariant (r1045) | a family clause is required for "adopt"; all efficacy counts have fixed-universe companions; largest-family size reported |

## 2. The intervention

**Switch:** environment variable `RUSTLE_READTHROUGH_JUNCTIONS` read by the assembler (`copy_assign --assemble-only`
and every path that builds pass-1 skeletons). Unset or `off` = BASE, byte-identical. Values `r` and `rq1`.

**Per contig, one pass, before pass-1 skeletons (so before locus formation):**
1. Reads: primary alignments only (`-F 2308`) with ≥ 1 intron (CIGAR `N`). Transcript strand = `+` iff (`ts` is `+`)
   XOR (read reverse), `ts` absent = `+` (the instrument's rule). 5′/3′ ends strand-aware.
2. Junctions (donor, acceptor, strand) with S ≥ 2 reads, **canonical** motif only (GT-AG, GC-AG, AT-AC on the
   junction's strand, from the genome FASTA). Junctions with S < 2 or a non-canonical motif are never flagged.
3. `U` = spliced reads whose 3′ end is inside the intron and whose 5′ end is upstream of the donor; `R = U/(U+S)`.
4. Start clusters: spliced reads' 5′ ends per strand, sorted; a new cluster when the gap to the previous 5′ end
   exceeds 100 bp; a real start = a cluster of ≥ 3 reads. `V1` = reads whose 5′ end is inside the intron and in a real
   start, whose 3′ end is beyond the acceptor, and whose FIRST intron (transcript orientation) has its donor inside
   the intron. `Q1 = V1/(V1+S)`.
5. Flag J: arm **R** iff `R ≥ 0.9524`; arm **RQ1** iff `R ≥ 0.9524` AND `Q1 ≥ 0.625`. Thresholds, gaps (100 bp),
   cluster size (3), scope (S ≥ 2, canonical) are FROZEN from the doc; no re-tuning at any stage.
6. **Action:** every alignment the assembler admits on that contig (primary or seeded good secondary) that carries a
   flagged (donor, acceptor, strand) is removed from the assembler's read set: pass-1, read-isoform widening,
   tied-secondary seeding and polish support. Flags are computed once from primaries and not recomputed.
   *Why removal, not a split at J:* a split would give the upstream piece a 3′ end at J's donor that no molecule
   has (r1078/r1079: trims manufacture chains no read carries), and r967's read splits moved nothing.

**Arms:** BASE (current default: `tools/rustle_pipeline.sh assemble` then `families`), R, RQ1, and a descriptive
**NULL**: the same number of alignments removed per contig, taken from reads carrying randomly chosen canonical
S ≥ 2 junctions that neither arm flags, matched to arm R per contig by junction read-count bin (log2), seed 20260925.
NULL is run for the assembly-level metrics (a)-(d) only.

**Implementation gates (before any arm is scored; a failure is fixed in the code, never in a threshold):**
- IV1 unset = BASE: `<id>.gtf` `cmp`-identical to a build of the previous commit on the gorilla OR6737 slice and
  human chr16; test suite passes.
- IV2 port = instrument: on the three dev BAMs of the doc (`advisor_jaccard/chr16.bam`, `ppar/chr20.bam`,
  `gw22/sec/ggo44.bam`), the Rust per-junction S, U, R, V1, Q1 equal `readthrough_rules.py`'s `junctions.tsv` on
  every row it writes (R and Q1 to 4 decimals). Where the script has a quirk, the port matches the script (the
  thresholds were frozen on its values) and the quirk is written into Amendment 1.
- IV3 action: 0 transcripts in an arm's GTF use a flagged junction; the log reports flagged junctions and removed
  alignments per contig; each arm runs under its own output prefix under `/mnt/linuxdisk/tmp/rustle_figures_dev/
  readthrough_ends/<arm>/<sample>/` with the run cache off (or keyed on the switch), and the log line naming the arm
  is checked, so no arm is served BASE's cached products. The molecules table (`<id>.molecules.tsv`) is shared.
- The binary's sha1 is recorded in Amendment 1 after IV1-IV3 and dev, BEFORE any held-out arm runs.

## 3. What "better locus representatives" means, and how each part is measured

Annotation = the species' RefSeq GFF (human `chm13v2.0_RefSeq_full`, never `HSA_genomic.gff`); gene set = `gene` +
`pseudogene` records minus records whose `description` contains "readthrough" (as the doc; gorilla has none).
Read support of an annotated record = ≥ 2 primary reads (`-F 2308`) with an aligned block (M/=/X) on its exon union
(the C2 / `_o1.primary_counts` rule); it depends on the BAM only, so every universe below is **identical across
arms**. Locus = one `gene_id` of the arm's GTF; its representative = `gtf_loci`'s choice (most reads, then longer
span, then last id); its span = union of its transcripts. Owner of a locus = the same-strand annotated gene with the
largest exonic overlap with the representative (ties: larger overlap with the locus's exon union, then gene ID).

**(a) Readthrough-fused loci — the efficacy metric. Fewer is better.**
- **FUSED (primary):** a locus holding ≥ 1 spliced transcript whose exon union (all its transcripts) overlaps by
  ≥ 1 bp the exon unions of ≥ 2 annotated genes on the representative's strand whose spans do not overlap each other
  (the locus-level form of the doc's RT junction). Count per substrate.
- Secondary: REP-FUSED (same, representative's exons only); SPAN-COVER (the span contains ≥ 50% of the exon union
  of ≥ 2 such genes; the families stage aligns the span); the junction-fused subclass (some transcript of the locus
  has a junction with donor in one gene's exons and acceptor in another's; the class the filter can reach); strata
  both-protein-coding vs any lncRNA/LOC; **absorbed genes** (fixed universe: read-supported genes that own no locus
  and whose exons are overlapped by a locus owned by another gene; r1015's lost label).

**(b) Representative ends vs TSS/TES.** Tolerances from §6w3's dispersion: **TES ±25 bp, TSS ±250 bp**.
- **Primary: TES-recovered genes** = read-supported annotated genes g (fixed universe) owning a locus whose SPAN 3′
  end lies within 25 bp of the 3′ end of any annotated transcript of g. Reported as a count and a fraction of the
  universe. *Why primary:* the filter acts at gene A's polyA site; 3′ ends are precise (§6w3, r806); the annotation
  is external; a fixed universe cannot be inflated by an arm changing its loci. The span is the object the families
  align.
- Secondary: TSS-recovered genes (span 5′ end, ±250 bp; 5′ is dispersed and library-biased, r973/r974, so never
  primary); the same two with the representative transcript's own ends; **read-derived** ends (annotation-free,
  partly circular because the assembler used the same reads; needed where the annotation is Gnomon-only): TES
  clusters = spliced primaries' 3′ ends, new cluster at a gap > 25 bp, ≥ 3 reads, not internal priming (< 12 A in the
  20 genomic bases downstream, strand-aware); TSS clusters = 5′ ends, gap > 100 bp, ≥ 3 reads. Reported as the number
  of read-derived clusters (fixed per sample) that are a locus span end within 25 / 250 bp, and (descriptive,
  prediction-conditioned) the fraction of loci whose span ends hit a cluster.

**(c) Transcript accuracy, as StringTie/FLAIR are judged.** `gffcompare -r <annotation GTF> <arm GTF>` through the
figures' restriction machinery (`figures/assembly.py`: same contigs for every arm; chimp/orangutan GTF from
`gff_to_gtf`). Decisions on raw counts: intron-chain precision = query multi-exon transcripts with class `=` / query
multi-exon transcripts (`.tmap`); intron-chain sensitivity = matching reference intron chains (`.stats`). Also
reported: transcript and locus level, total transcripts, and multi-gene transcripts (the transcript-level form of
FUSED).

**(d) Loci in the Liftoff framework** (`docs/PREREG_liftoff_loci_2026-09-25.md` §3, C2). Reference loci = the
species' Liftoff table (in place + moved + extra copies, exon union ≥ 200 bp), restricted to read-supported loci of
the sample. Found iff cov(G | R) ≥ 0.5 for some de novo locus exon union R. **Guard metric:** found fraction,
annotated and extra copies separately. Descriptive: fraction of de novo loci lying at a Liftoff locus (cov(R | G)
≥ 0.5) and reciprocal one-to-one matches (both ≥ 0.5), where splitting fused loci should show. **Until a species'
Liftoff table is merged,** the annotated part is computed on the annotation's own records (Liftoff places ≥ 99% in
place by claim L1) and labelled "annotation-only"; extra copies wait for the table.

**(e) Families vs external references** (the families stage unchanged on each arm's GTF: `mcl_families --from-gtf
--min-exonic-bp 1 --min-shared-exon-frac 0.60`, MCL 2.8, through `tools/mm2_shard.sh`). Restriction to a substrate
by the rule of `PREREG_genome_wide_families_2026-09-25.md` §1 (families built genome-wide, then restricted).
- **Human:** Ensembl Compara families at Primates (that prereg §3.2), `family_score --chrom ALL --per-family
  --pairwise`: one-to-one bipartite sensitivity, precision (an upper bound; not tie-invariant, r1045, same scorer
  build for every arm), F. Also Liftoff copy pairs (pair recall), Soto 2025 (labelled not independent; descriptive).
- **Apes:** Liftoff copy pairs (`sequence_ID` ≥ 0.95, both loci read-supported; recovered iff covering copies share a
  family): pair sensitivity only (the relation is not a partition: no precision, no bipartite matching).
- Also reported: number of families, largest family (hub check, r846/r913), loci entering the graph.
- ⚠ r1017's oracle (perfect split of every chr16 fused locus) bounds the gain at **+0.021** referee F, and r953's
  readthrough-free simulation at +0.0139. The filter removes a minority of readthrough reads. **A family gain is
  therefore not expected and not required**; the family clause is non-inferiority.

## 4. Substrates

**DEV (reported, never in the verdict):** human A119b chr16 and chr20, gorilla OR6737 NC_073244.2, i.e. the doc's
BAMs, as restrictions of the genome-wide arms. ⚠ These are **no longer untouched**: chr16 chose θ and θ_Q1, and chr20
and NC_073244.2 were read as held-out by five successive refinements (R, R_peak, R & Q, R & Q1, tags). They are
development contigs now. A bug found on dev may be fixed (Amendment 1); no threshold may move.

**HELD-OUT verdict substrates (each judged on its own; never pooled, species never pooled):**

| id | substrate | why |
|---|---|---|
| V1 | human_A119b genome minus chr16, chr20 (annotated contigs) | same library as dev, chromosomes this rule never saw |
| V2 | human_testis whole genome | a second human library, never used for any decision (mapped without `-uf`; strand from `ts` handles it) |
| V3 | gorilla_OR6737 genome minus NC_073244.2 | same library as the gorilla dev contig, unseen contigs |
| V4 | gorilla_KB3781 whole genome | fibroblast: many downstream genes silent, the RQ1 recall-loss class (doc §8) |
| V5 | chimp_PTR whole genome | never used |
| V6 | orangutan_PPY whole genome | never used |

Secondary rows (no verdict): human A119b chr6 alone (scored once, no rule chosen, r1113) and chr19 alone (reserved
unexamined by §6u7 for node construction); human A119b genome minus every contig used for a family decision (S2 of
the genome-wide families prereg); whole genomes (S0).

**Power floors (fixed now; a clause below its floor on a substrate is reported, not judged there):** (a) BASE FUSED
≥ 50; (b) universe ≥ 1,000 genes; (c) BASE matching reference chains ≥ 500; (d) ≥ 1,000 read-supported reference
loci; (e) human ≥ 30 scored Compara families, apes ≥ 30 read-supported Liftoff pairs. A substrate is a verdict
substrate when (a) qualifies. Any verdict other than "undecided" needs ≥ 1 human and ≥ 1 ape verdict substrate.

## 5. Bar for "beneficial" (per arm, each clause on every qualifying verdict substrate)

| clause | measure | passes iff |
|---|---|---|
| **A1** efficacy | FUSED loci, arm vs BASE | reduction ≥ **X = 10%** relative, and NULL's reduction < half the arm's |
| **G1** precision | intron-chain precision (c) | arm ≥ BASE (raw counts) |
| **G2** sensitivity | matching reference intron chains (c) | arm ≥ BASE × (1 − **Y**), **Y = 1%** |
| **G3** ends | TES-recovered genes (b) | arm ≥ BASE |
| **G4** loci kept | found read-supported Liftoff loci (d), annotated and extra copies separately | arm ≥ BASE × (1 − Y) |
| **G5** families | human: bipartite F, sensitivity, precision vs Compara; apes: Liftoff pair recall (e) | each ≥ BASE − **0.005** |

**Why X = 10%.** The doc's held-out junction recall is .42-.52 (R) and .25-.36 (RQ1), and R removes 24-40% of
readthrough READS. A locus separates only when every bridging chain goes (H:453: median read min-cut 12), and R
keeps the dominant, well-read fusions by design, so the locus-level effect must be smaller than the junction-level
one. 10% sits below every junction-level rate (lowest: RQ1 gorilla .25) and near what the rates give for a fused
locus with two bridging junctions removed independently (.25 × .42 ≈ .10): a miss means the junction effect does not
reach the loci, not that the bar asked for the junction rate. Below 10%, r1017's +0.021 ceiling for fixing ALL fused loci implies < ~0.002 F, not worth a
default. **Why Y = 1%.** The doc's worst real-junction removal is .0044 (R, dev) and .0017 (RQ1): a matched chain or
locus is lost only through a removed junction, and the removed real junctions are alternative last exons used by one
long isoform; Y = 2 × .0044 rounded up allows each false removal to cost ~2 matched objects and no more. **Why
0.005.** About ¼ of r1017's oracle gain (+0.021) and ⅓ of r953's no-readthrough gain (+0.0139): an arm costing more
than that has spent more than the best case of any readthrough repair. Family efficacy is not a clause (§3e).

**Verdicts (per arm):**
- **Adopt as default — recommendation only (the flip is the user's call):** A1 and G1-G5 pass on every verdict
  substrate where they qualify; G5 qualifies on ≥ 1 human and ≥ 1 ape verdict substrate; no NPIP cap (§6). A flip
  changes the O1 default loci, hence every downstream family and copy-assignment figure (re-run required).
- **Refute (as a representative lever):** A1 fails on more than half of the verdict substrates (the junction effect
  does not reach the loci), OR on any verdict substrate G2 or G4 loses more than 2Y (2%) or a G5 measure drops by more
  than 0.010, OR G1 or G3 is lower on more than half of them. Register row; the switch stays opt-in until a cleanup
  wave removes it with the user's agreement.
- **Keep opt-in:** everything else (including: G5 not measurable because the Liftoff tables or the Compara families
  are not built; the NPIP cap). The numbers are recorded with the switch.
- **If both arms reach "adopt"**, the recommendation is the arm with the larger median A1 reduction over the verdict
  substrates (RQ1 on a tie, the more conservative). Arms are compared with BASE, never tuned against each other.

## 6. Reported beside the bar

- NPIP (thesis family): human A119b and testis, whole genome — the Compara Primates family holding NPIP-named genes
  and the U2 NPIP set (chr16, a dev contig): members found, bipartite F. **Cap:** if NPIP's F or members found drop in
  both human samples, the verdict cannot exceed "keep opt-in" (r1017/r1018 predict this risk). PKD1P6-NPIPP1 is
  expected to be kept by the rule (dominant fusion).
- Flag census per substrate: flagged junctions and removed alignments, by annotation class (RT / ANN / novel, labels
  only); flagged junctions by intron length; the 20 most-read flagged junctions with gene pairs.
- Loci and transcripts per arm; loci lost (BASE loci with no arm locus overlapping their exons) and gained.
- Cost: wall time and peak RSS per arm and sample vs BASE.
- Copy assignment (O2) consumes the families but is **not** measured here.

## 7. Predictions (before any number)

1. A1: R reduces FUSED by 10-25% (median over verdict substrates), RQ1 by 5-15%; RQ1 fails A1 on ≥ 1 gorilla
   substrate. NULL ≈ 0. Most removed fusions are the lncRNA/LOC and minor run-on class, not the comparable-genes class.
2. G1: precision up 0.0-0.3 pt (removed chains are unannotated). G2: R loses 0.1-0.5%, RQ1 < 0.1%.
3. G3: TES-recovered genes up 0.2-1%; TSS-recovered up in RQ1 (B's own start cluster is a condition of the flag).
4. G4: R loses silent downstream genes most on KB3781 (fibroblast) and may fail there; RQ1 passes.
5. G5: |ΔF| ≤ 0.005 everywhere; NPIP unchanged.
6. Overall: RQ1 "keep opt-in" (A1 short somewhere, or G5 unmeasurable on apes), R "keep opt-in" or refuted on G4.

## 8. Order and stop rules

1. Implement behind the switch; IV1-IV3. 2. Dev arms and every metric on dev (reported). 3. Record the binary sha1
(Amendment 1). 4. Held-out arms, then scoring; no code or threshold change after step 3 except a bug fix that
re-runs every arm, recorded as an amendment. Machine rules: one heavy process at a time in the foreground; outputs and
`TMPDIR` under `/mnt/linuxdisk`; human A119b assembly ~7.5 min per arm. A step that cannot fit is split or stops with
a note; a clause not computed is "not measured" and caps the verdict at "keep opt-in". Every number is a new register
row; the doc's rows are not overwritten.

## Amendments

### Amendment 1 — implementation (2026-09-25, written before any arm number exists)

Written by the implementation step before the port was run on any BAM with the switch set. The only runs so far are
IV1's unset / `off` byte-identity runs; no flag census, removal count, assembly or metric number of any arm exists.
Nothing here changes an arm, a threshold, a metric or the bar. Each item states how a line of §2 is realised where
the text is ambiguous or cannot hold literally.

1. **Values and scope.** `RUSTLE_READTHROUGH_JUNCTIONS` unset, empty or `off` = BASE; `r` or `rq1` (any case); any
   other value is fatal at startup. It is accepted only with `copy_assign --assemble-only` (fatal otherwise). §2's
   "every path that builds pass-1 skeletons" is realised as the assembler's pass 1, both the streaming reader (the
   pipeline's default) and the buffered one (`--materialize-reads` and every other non-streaming assemble-only
   setting). That is the only pass 1 the arms run (`rustle_pipeline.sh assemble`, then `families` on its GTF). The
   family-detection pass 1 inside `detect_and_assign` / `gw_family_catalog` is not touched.
2. **Thresholds as integers.** "R ≥ 0.9524" is implemented as `U ≥ 20·S`, and "Q1 ≥ 0.625" as `3·V1 ≥ 5·S`.
   θ = 0.9524 is 20/21 printed to 4 decimals: the doc read θ off the TSV's 4-dp `R`. A literal float test
   `U/(U+S) ≥ 0.9524` would drop every `U = 20·S` tie, and the doc's flagged sets contain such ties (chr16
   SMG1P5→NPIPB13: S = 8, U = 160, printed R 0.9524). `3·V1 ≥ 5·S` is exactly `Q1 ≥ 5/8`. The integer forms differ
   from the script's rounded comparison only where S ≥ 74 (`U = 20·S − 1` prints as .9524) and, for Q1, S ≥ 938.
   The doc's tables hold no canonical R-pass junction with S ≥ 74 on any dev contig. IV2 reports any junction in
   those gaps.
3. **The script's first-donor lookup is reproduced (V1).** The script stores each read's first donor under (strand,
   5′ end, 3′ end, read name), then counts a read when ANY stored donor for its (strand, 5′ end, 3′ end) lies inside
   the intron. So a read counts toward V1 when any read with the same strand and both ends has its first exon ending
   inside J's intron. IV2 says the port matches the script, so the port does the same. It differs from the per-read
   wording of §2.4 only for reads that share both ends with another read.
4. **Removal matches the junction's coordinates, not the read's own `ts` strand.** Pass 1 groups reads by (contig,
   intron chain) with no strand, and a spliced model takes its strand from its junction motifs. A flagged junction
   is canonical on exactly one strand, because the two strands' canonical motif sets are disjoint. A chain holding
   its (donor, acceptor) can therefore only be emitted on the flagged strand. §2.6's "carries a flagged (donor,
   acceptor, strand)" is realised as: the alignment's intron chain holds a flagged (donor, acceptor), whatever the
   alignment's own `ts`/FLAG strand. This is what IV3 (0 transcripts use a flagged junction) needs.
5. **Population of the statistics.** Statistics are computed per REGION from the records overlapping it, which is
   per contig under `--genome-wide` (the pipeline). A record overlapping an earlier window of the same region is
   counted once. Primary means mapped, not secondary, not supplementary (the script's `-F 2308`; both keep
   QC-fail and duplicate flags). Spliced means ≥ 1 intron from the pool's own CIGAR parser
   (`bam::exons_from_cigar`). Its edge cases can differ from pysam's one-intron-per-`N`: an `N I N` run becomes one
   intron, and a leading or trailing `N` is dropped. IV2 lists every row where this matters. Removed alignments
   are counted in the pool, after the pool's own rules (seeded good secondaries included, coordinate duplicates
   once). The statistics are counted before those rules, as §2.1 says.
6. **Outputs, only when set.**
   - `<out>.readthrough_junctions.tsv` has one row per flagged junction:
     - `contig`;
     - `donor`/`acceptor`: the 1-based inclusive intron, i.e. the script's `start`/`end`, genomic left/right on
       both strands;
     - `strand`, `S`, `U`;
     - `V1`: Q1's read count, as §2 names it;
     - `rule`: `rq1` passes R and Q1, `r` passes R only. `r` rows appear only in arm R.
   - `params.tsv` gets the rows `readthrough_junctions`, `readthrough_junctions_flagged`,
     `readthrough_junctions_chains_removed` and `readthrough_junctions_alignments_removed`.
   - The log prints one `[readthrough]` line per region: junctions with S ≥ 2, canonical R-pass, flagged, pool
     chains and alignments removed (IV3).
   - `RUSTLE_READTHROUGH_JUNCTIONS_ALL=1` adds `<out>.readthrough_junctions.all.tsv`, the IV2 instrument: every
     junction with S ≥ 2, with S, T, U, V1 and canonical.
7. **Not in this step.** The NULL arm (§2) has no switch value yet. It must be implemented, under its own value,
   before any arm is scored. The binary sha1 (§2, after IV1-IV3 and dev) is appended to this amendment when dev is
   done.

### Amendment 2 — the NULL list (2026-09-25, written before any R, RQ1 or NULL number exists)

Written by the scoring step (`bench/mechanism/readthrough_eval.py`). No assembly and no metric of arm R, RQ1 or NULL
exists. The only numbers so far are the scorer's end-to-end test on the development contigs chr20 and chr16: BASE's
metrics there (development, never in the verdict, §4) beside synthetic TEST arms (transcripts dropped from BASE after
assembly; not an arm of §2). Nothing below depends on them. No arm, threshold, metric or bar changes. This fixes how
§2's NULL is drawn, because "the same number of alignments removed per contig" cannot hold literally:
- the NULL removes every alignment carrying a drawn junction, so its count moves in steps of one junction's reads;
- arm R's removal count is taken in the assembler's pool (seeded good secondaries, duplicates once; Amendment 1
  item 5). The same count for the NULL's junctions is not known before the NULL runs.

1. **Population.** Per contig, from the sample's BAM: primary spliced alignments (`-F 2308`, ≥ 1 `N`, pysam's
   CIGAR parse). Strand = `ts` (absent = `+`) XOR FLAG 0x10, the instrument's rule. A junction is (1-based
   inclusive intron start, end, strand). S = the alignments carrying it.
2. **Sets.** F = every row of arm R's `<out>.readthrough_junctions.tsv`. RQ1 ⊆ R, so "neither arm flags" means
   "not in F". The candidates are the junctions with S ≥ 2, canonical on their strand (GT-AG, GC-AG, AT-AC), and not
   in F.
3. **Target.** T = the distinct primary alignments carrying ≥ 1 junction of F. Both sides are counted in the same
   population. The pool counts of R and NULL (params.tsv) are reported beside T; they are not matched.
4. **Draw.** The RNG is `random.Random("20260925:<sample>:<contig>")`.
   - F is sorted, shuffled once, then visited in that order, cycling.
   - For each visited f, one candidate is drawn uniformly, without replacement, from the candidates whose
     floor(log2 S) equals f's. An empty bin falls back to the nearest non-empty bin (the lower one on a tie).
   - The draw stops as soon as the distinct primaries carrying the drawn junctions reach T (the last junction may
     overshoot by its own reads), or when the candidates run out.
5. **Use.** `readthrough_eval.py null` writes the list and a per-contig summary (T, achieved, junctions).
   - The list's columns are chrom, intron_start_1b, intron_end_1b, strand, S, bin, and the R junction it was
     matched to.
   - The NULL arm removes the alignments whose intron chain holds a listed (donor, acceptor), exactly as Amendment 1
     item 4 removes flagged junctions. The switch value that reads the list is the implementer's (Amendment 1
     item 7).
   - One NULL serves both arms (§2 matches it to R). A1's null part for RQ1 is judged against that same NULL.
6. **Procedure check (dev chr20; F = the design's scratch R set, since no arm dump exists).**
   - 540 flagged junctions and 15,050 candidates; T = 2,405, achieved 2,407.
   - Junctions per log2 bin 1-5: R 380/97/34/23/6, NULL 378/96/34/23/6.
   - The list is identical on a rerun (md5 87b4316e).

### Amendment 3 — the NULL arm's switch value and the frozen binary (2026-09-25, 20:10; before any arm metric)

Written by the implementation step that closes Amendment 1 item 7. The file's md5 before this edit was ab86850a. Only
dev slices were run (gorilla `gw22/sec/ggo44.bam` = NC_073244.2, human `ppar/chr20.bam`). The numbers that exist are
listed in item 4: byte-identity checks, the NULL draws on dev and pool removal counts. No metric (a)-(e) exists for
R, RQ1 or NULL on any substrate. Nothing here changes an arm, a threshold, a metric or the bar.

1. **Switch value.** `RUSTLE_READTHROUGH_JUNCTIONS=list:<path>`. The `list:` prefix may be in any case. As in
   Amendment 1 item 1, it is accepted only with `copy_assign --assemble-only`, streaming and buffered.
   - **File:** tab-separated. `#` lines and blank lines are skipped. The first other line is the header. It must name:
     - a contig column (`contig` or `chrom`);
     - a donor column (`donor` or `intron_start_1b`);
     - an acceptor column (`acceptor` or `intron_end_1b`);
     - `strand`.

     Donor and acceptor are the 1-based inclusive intron. The header therefore reads both
     `readthrough_eval.py null`'s list and an arm's own `<out>.readthrough_junctions.tsv`. Other columns are ignored
     and duplicate rows count once.
   - **Fatal at startup, before any read:** any of these ends the run.
     - an unreadable file;
     - no header, or a missing column;
     - a non-integer coordinate;
     - `donor < 1` or `acceptor < donor`;
     - a strand other than `+`/`-`;
     - no junction row;
     - a contig that is not a reference sequence of the BAM header;
     - an intron beyond the contig's length.

     A list for another assembly or naming scheme would remove nothing and run BASE under the name NULL.
2. **Semantics.**
   - **Flagged set.** A region's flagged set is every listed junction whose intron overlaps one of the region's
     windows. Under `--genome-wide` that is every listed junction of the contig, whatever its statistics.
   - **Removal.** Removal is Amendment 1 item 4 unchanged: the same `ReadthroughFlags` set drives `drop_chains_with`
     (streaming) and `filter_pool` (buffered).
   - **Statistics and outputs.** The statistics are still collected. Each listed junction is scored in full, so the
     outputs are Amendment 1 item 6's, with these differences:
     - `<out>.readthrough_junctions.tsv` has one row per listed junction of an assembled region, with its S, U and V1
       and `rule` = `list`;
     - `params.tsv` has `readthrough_junctions` = `list`;
     - `params.tsv` also has two rows written only by this arm: `readthrough_junctions_list` (the path) and
       `readthrough_junctions_listed` (distinct junctions in the file).
   - `[readthrough]` log lines keep their fields. "canonical R-pass" is still the R count, so a NULL junction that
     passes R would show; none does on dev.
3. **A consequence of Amendment 1 item 5, reported and not fixed.**
   - The pool parses an `N I N` CIGAR run as ONE intron.
   - So a listed junction seen only inside such a run is absent from every pool chain and cannot be removed.
   - Arm R cannot flag such a junction either, because it uses the same parse.
   - On dev chr20 this is 1 of 537 NULL junctions: 31,363,270-31,491,038 +. The scorer counts S = 2 for it and the
     assembler S = 0; its 2 reads carry `1195N 280I 127769N`. Gorilla has none.
4. **Checks (dev only).** The previous build is `copy_assign` 3abab3df, Amendment 1's binary, rebuilt bit-identical
   from the tree before this change.
   - **IV1, streaming.** Settings unset, empty, `off`, `r` and `rq1`; ggo44 and chr20; seeded and
     `--no-seed-secondaries`. Driver `tools/rustle_pipeline.sh assemble --no-cache`. All 128 products (gtf, params,
     quant, families, assignments, famcn_readonly, readthrough TSV) are `cmp`-identical to the previous build, and the
     file sets are identical.
   - **IV1, buffered.** `--materialize-reads`, unset: 12 of 12 identical.
   - **Same drop path.** `list:` of an arm's own dump reproduces that arm in all 8 cases (r and rq1 × 2 slices × 2
     modes):
     - gtf, quant, families, assignments and famcn are identical;
     - the TSV's columns 1-7 are identical;
     - params differ only in the rule value and the two list rows.
   - **NULL lists.** Drawn by `readthrough_eval.py null` (sha1 0a50ab11) from arm R's dump.
     - chr20: T = 2,405, achieved 2,407 with 537 junctions. The rows are identical to Amendment 2's procedure-check
       list.
     - NC_073244.2: T = 506, achieved 512 with 137 junctions.
     - The assembler's S equals the list's S on 536/537 junctions (chr20; the exception is item 3) and on 137/137
       (gorilla). No listed junction has `U >= 20·S`.
   - **NULL pool removals, alignments / chains** (the pool counts Amendment 2 item 3 reports beside T):

     | | NULL | R |
     |---|---|---|
     | gorilla, seeded | 538 / 386 | 499 / 319 |
     | gorilla, not seeded | 504 / 356 | |
     | chr20, seeded | 2,096 / 1,488 | 2,062 / 1,553 |
     | chr20, not seeded | 2,069 / 1,478 | |

   - **IV3 for NULL.** 0 transcripts use a listed junction. BASE had 72 / 70 (gorilla, seeded / not) and 210 / 211
     (chr20).
   - **Streaming = buffered for NULL** (seeded, ggo44 and chr20): params, TSV, quant, families, assignments and famcn
     are identical. The GTF is identical except the O2-only `matched_reads` attribute, which differs in the same way
     when the switch is unset.
   - **Tests.** `cargo test --release`: 892 passed, 0 failed, 13 ignored. That is 5 new tests: list parsing (both
     formats, round trip), fatal values and files, listed junctions dropped and unlisted kept, the window restriction,
     and streaming = buffered = the rq1 arm on the BAM fixture.
5. **Frozen binary for the held-out arms (§2, §8 step 3).**
   - `copy_assign` sha1 **452454c68b14e7541e6a2037bd6fdec0e0a2ea47**, from `cargo build --release --bin copy_assign
     --bin as_table` with target `rustle_target_dev`, a reproducible build.
   - `as_table` is f9556f36, unchanged.
   - Sources: the uncommitted tree over 2e64a5a3, with `src/bin/copy_assign.rs` sha1 0471cc85 and
     `src/rustle/vg_family/denovo_assemble.rs` f8fd8ca2.
   - Copies are kept in `/mnt/linuxdisk/tmp/rustle_figures_dev/rtB_null/bin_new/`.
   - ⚠ The `copy_assign` left in `target/release` by `cargo test` has a different sha1 (cd789e9b), probably because
     the test build unifies dev-dependency features. Use the kept copy; do not use whatever a test run leaves behind.
   - This binary serves every arm: BASE (unset), R, RQ1 and NULL (`list:<NULL list>`). A held-out NULL list is drawn
     from that substrate's own arm-R dump.

## Outcome (2026-09-25 21:55; frozen binary copy_assign sha1 452454c6; held-out substrates only)

`bench/mechanism/readthrough_eval.py verdict` over the six held-out tables (`/mnt/linuxdisk/tmp/rustle_figures/rt_arms/tables/`):
**R = keep opt-in** (median fused-locus reduction 0.174; passes every measured clause on 5/6 substrates, A1 misses on
chimp_PTR at −9.5% vs the 10% bar; G4 extra copies and G5 families NOT MEASURED ⇒ capped at opt-in by §5);
**RQ1 = refute** (A1 fails on 5/6 substrates; median reduction 0.093).

| substrate | fused loci BASE → R (RQ1, NULL) | intron-chain precision | TES-recovered genes | matched chains | loci kept |
|---|---|---|---|---|---|
| human_testis | 200 → 165 (187, 198) | .4596 → .4616 | 5,069 → 5,112 | 11,399 → 11,357 | 5,709 → 5,684 |
| human_A119b −chr16/chr20 | 2,050 → 1,626 (1,827, 1,913) | .1765 → .1795 | 9,792 → 10,216 | 35,956 → 35,796 | 21,547 → 21,461 |
| gorilla_OR6737 −NC_073244.2 | 549 → 454 (495, 538) | .3496 → .3511 | 6,985 → 7,077 | 24,277 → 24,180 | 13,178 → 13,134 |
| gorilla_KB3781 | 662 → 466 (560, 646) | .3594 → .3632 | 6,177 → 6,295 | 25,911 → 25,709 | 12,398 → 12,324 |
| chimp_PTR | 420 → 380 (405, 411) | .3915 → .3927 | 8,378 → 8,419 | 21,031 → 20,952 | 12,469 → 12,441 |
| orangutan_PPY | 687 → 579 (627, 674) | .1931 → .1942 | 9,053 → 9,221 | 20,243 → 20,175 | 13,750 → 13,709 |

Reading: at the transcript and locus level R is beneficial everywhere measured (fewer readthrough-fused loci, higher
precision, more genes whose annotated 3′ end is recovered, ≤ 0.8% of matched chains and ≤ 0.6% of loci lost; the
matched random removal changes fused loci by 1-7%). The family-level clause is unmeasured, so no default flip is
recommended by this prereg; measuring G5 needs genome-wide families of BASE and R on ≥ 1 human and ≥ 1 ape sample.

### Outcome addendum (2026-09-26 01:46): the family clause G5 on human_testis and chimp_PTR

Genome-wide families (`tools/rustle_pipeline.sh families`, frozen binaries, shard wrapper) of BASE and R:
- **human_testis vs Ensembl Compara families at Primates (426 families):** bipartite F .2740 → .2751, sensitivity
  .1599 → .1607, precision .9545 → .9548 (Soto 2025, descriptive: unchanged .2385). Non-inferior (tolerance −.005): pass.
- **chimp_PTR vs Liftoff (record, extra copy) pairs:** pair recall .1172 → .1172 (equal): pass. Loci found in the Liftoff
  framework .5464 → .5451 annotated (−0.2%), extra copies equal; one-to-one reciprocal loci 10,156 → 10,192.
- Families 338 → 337 (human testis), 379 → 378 (chimp); largest family unchanged (no hub).
**Verdict (unchanged): R = keep opt-in.** The family clause passes where measured, but A1 misses on chimp_PTR
(−9.5% vs the 10% bar), so "passing everything" cannot be reached; G4 extra copies / G5 remain unmeasured on the other
four substrates and would not change that. The default flip is the user's call.
