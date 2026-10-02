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

## Amendment 1 (2026-10-01): IsoCon, the reference-free arm (written after the chr20 + chrX pilot of the caller, before any IsoCon run)

The user asked to add IsoCon (Sahlin et al. 2018), the reference-free clustering of Iso-Seq reads into distinct transcripts, as a
head-to-head with the PSV allele caller, on the same frozen truth, and to keep it to what this machine can run.

**Pilot facts that set the scope (seen before writing this):** fibroblast reads touching the S_fam copies (any alignment record): NPIP
2,951 distinct reads (1,040 primaries), TBC1D3 1,264 (5 primaries: TBC1D3 is essentially not expressed in fibroblasts). The caller pilot
(chr20 + chrX) gave f_X = 2/439 and S_single precision 588/590; it does not touch these rulings.

- **Arm I, real reads (NPIP; TBC1D3 reported for completeness).** Reads = every distinct read with any alignment record (primary or
  secondary) overlapping an S_fam copy span of the family: a broad net, so a copy the reference lacks still contributes its reads through
  whichever copy they best resemble. Sequence = the read's primary record, restored to the read's own orientation. Above 3,000 reads a
  family is downsampled at random (seed 1) to 3,000. IsoCon `pipeline` with its default parameters (version 0.3.3, PyPI), no reference.
- **Arm S, simulation (NPIP and TBC1D3).** Truth transcripts = each S_fam copy's exon union on each haplotype: haplotype A = the `_pri`
  exon sequence; haplotype B = the same exon blocks lifted through the frozen `_pri` -> B alignments; plus the B-only locus (lifted
  from the family copy whose hit it is). 20 reads per haplotype-copy transcript, 1% errors (0.5% substitutions, 0.25% insertions, 0.25%
  deletions), 5' end shortened by an exponential with mean 100 bp, 3' end by one with mean 10 bp; seed 1. IsoCon `pipeline`, defaults.
- **Scoring, both arms** (descriptive; no pass/fail):
  - every IsoCon transcript is mapped to the maternal and paternal assemblies (`minimap2 -c -x splice` to the frozen splice indexes) and
    assigned to the copy whose locus its best hit overlaps and to the haplotype with the higher identity (equal = both);
  - **recovered** = haplotype copies with a transcript at identity >= 0.999; **alleles separated** = T2d copies with one transcript on
    each haplotype; **alleles merged** = T2d copies with transcripts on one haplotype only; **paralogs merged** = a transcript whose best
    hits on two copies tie; **unmatched** = transcripts below 0.99 to every copy;
  - **lower bound check**: IsoCon's count of distinct coding groups (its own copy estimate) against T (NPIP 41, TBC1D3 28; real arm:
    against the haplotype copies with reads);
  - in Arm S, recall and precision against the known truth transcripts (exact match, then >= 0.999).

**Amendment 2 (2026-10-01, after building the Arm S truth and before any IsoCon run).** Two defects in Amendment 1's truth transcripts:
(a) an exon union merges every annotated isoform of a copy, reaching 36 kb for NPIP, far beyond Iso-Seq read lengths; (b) lifting a block
by the min/max target position over every record that touches it produced a 1.6 Mb "transcript" when one block touched two records far
apart. Changes: the simulated transcript of a copy is its primary annotated transcript in the 2026-09-29 truth GTF (most junctions, ties:
longest); each exon block is lifted through the single primary record that covers most of its bases, and the lift is rejected when that
record covers < 95% of the block or the lifted span exceeds 1.2 x the block length + 100 bp. The caller's truth (T2d/T2i/T1 per copy) is
unchanged. Scoring of both arms is unchanged.

**Amendment 3 (2026-10-01, after the TBC1D3 simulation ran and the NPIP simulation timed out; no NPIP IsoCon output exists).** The gorilla
RefSeq NPIP models are 6-27 kb, while the real fibroblast reads at NPIP have median 3,820 bp (p90 4,367); IsoCon on 880 reads of up to
27 kb did not finish in the 10-minute window. NPIP simulated transcripts are therefore the 3'-most 4,000 bp of each truth transcript (A, B
and B-only alike); TBC1D3 is unchanged. Real arm: a family whose IsoCon run does not finish in one 10-minute call is downsampled at random
(seed 1) by halves until it does; the size used is reported. The TBC1D3 simulation result stands as run.

## Amendment 4 (2026-10-01): deletion test on NPIP — does a copy missing from the reference pass as an allele? (written before any masking)

User request: remove a copy that is well covered by its own reads and see what the RNA-only methods make of it. Prior art in this project:
the 2026-08-14 whole-genome excision of 162 two-copy families (a deleted copy's reads are ABSORBED by one paralog 64% / ORPHANED 33%;
the S2 divergence detector TPR 0.27). New here: NPIP, and the allele framing.

- **Copy:** the S_fam NPIP copy with the most primary reads over its exons in the fibroblast BAM = **NPIPA2** (gN15, 183 primaries;
  NC_073242.2:32,426,793-32,456,482, minus strand; truth class T2i). Its span is hard-masked to N in a copy of `_pri` (no other NPIP copy
  overlaps it; it also covers RefSeq LOC101141855 = the NPIPA2 record, LOC115932701, and the first 12.7 kb of LOC129527628, disclosed).
- **Reads:** the 2,951 NPIP-net reads, realigned with the fibroblast BAM's `@PG` command (`-ax splice:hq -uf --eqx -Y -N 50 -p 0.1
  --secondary=yes`) to the masked genome AND to the unmasked genome (local baseline, same minimap2, so the comparison is paired).
- **R1 fate:** the reads whose baseline primary lies on NPIPA2's exons: unmapped / landing copy (the copy holding most of them) /
  concentration on that copy / their mismatch rate there.
- **R2 fake allele:** at the landing copy, the registered allele-site test (>= 10 reads, minor >= 3 and >= 0.20, not a PSV among the
  copies left in the reference, not A>G, not near a homopolymer) on its primary reads, baseline vs masked. A site that appears only
  in the masked arm and whose minor base is NPIPA2's base is a fake allele.
- **R3 no-reference-match detector:** H4's read groups at the landing copy, with NPIPA2 removed from the paralog list, baseline vs masked.
- **R4 IsoCon:** net reads re-selected from the masked alignment (any record on a remaining NPIP copy, plus the net reads it leaves
  unmapped), IsoCon `pipeline` defaults; an output is flagged "not in the reference" when it has no hit at identity x coverage >= 0.999
  in the masked genome. True flag = the output matches NPIPA2 at >= 0.999 in the unmasked genome. Report true and false flags.

## Amendment 5 (2026-10-01): IsoCon's reference-free transcripts as extra copies for O2 (written before any of these runs)

User request: use IsoCon's machinery for O2. The part tested: the copy set. O2 can only choose among reference copies, so a missing copy's
reads are forced onto a paralog (Amendment 4: 97% of NPIPA2's reads on LOC124907808).

- **Extra copies:** the 35 IsoCon outputs flagged "not in the reference" in Amendment 4's R4, appended to the NPIPA2-masked genome as
  contigs `iso_<k>` (transcript orientation). No other change to the genome.
- **Reads:** the 2,951 NPIP-net reads, aligned with the BAM's `@PG` command to (R) the masked genome [`masked.bam`, already made] and
  (R+I) the masked genome + the 35 contigs.
- **O2:** shipped `copy_assign --families` (defaults; O2 scope = AS-tied reads). Catalog R = the 24 NPIP copies left in the reference
  (exon unions, as in the main catalog). Catalog R+I = the same 24 + the 35 contigs, each one exon spanning the contig.
- **Final call per read:** O2's copy where O2 rows the read (`assigned`); abstain where O2 rows it as `ambiguous`/`tied`; the copy holding
  its primary alignment where O2 does not touch it (not AS-tied).
- **Truth:** NPIPA2's 176 reads (baseline primary on its exons) belong to NPIPA2; an `iso_<k>` copy "is NPIPA2" when its sequence matches
  NPIPA2 at >= 0.999 in the unmasked genome (16 of the 35, from R4). Every other read belongs to the copy of its baseline primary.
- **Readouts:** for NPIPA2 reads in R and R+I: right (an NPIPA2 `iso` copy) / wrong (LOC124907808, another copy, or a non-NPIPA2 `iso`) /
  abstain. For the other NPIP reads: how many move onto an `iso` copy that is not theirs (false moves) or onto any `iso` copy at all.
- **Rule (fixed now):** the extra copies **help** if, for NPIPA2 reads, wrong falls by >= 50% from R to R+I, and false moves of the other
  reads stay <= 5% of them; **hurt** if false moves exceed 10%; otherwise **mixed**.

## Amendment 6 (2026-10-01): held-out test of Amendment 5 on the 2026-08-14 excision panel (written before any of these runs)

Amendment 5 was developed on NPIP (one deletion). Held-out substrate, never used for it: the 162 two-copy gorilla families of the
2026-08-14 whole-genome excision (`winloci_scratch/o3_excise/PREREG.md`): in each family the copy with the larger start was hard-masked
(all 162 in one genome, `o3_excise/GGO.masked.fasta`), the other kept; 348,046 reads realigned (`o3_excise/masked.bam`, minimap2
2.30 = the version used here). Known fates (2026-08-14): absorbed 104, orphaned 54, scattered 4.

- **Reads per family f (scored set S_f):** baseline primaries (`o3_excise/panel_primary.bam`) on the masked copy D (D reads) and on the
  kept copy K (K reads); up to 500 of each, random, seed 1. Sequences from `o3_excise/panel_reads.fq`.
- **IsoCon input (net_f):** the reads of S_f with any masked-arm record overlapping K, plus the reads of S_f the masked arm leaves
  unmapped. D reads absorbed by another locus are not given to IsoCon (the method is family-scoped; disclosed). IsoCon `pipeline` defaults.
- **Flags and contigs:** IsoCon outputs with no hit at identity x coverage >= 0.999 in the masked genome (`minimap2 -c -x splice:hq -uf`)
  become contigs `iso_<fam>_<k>` appended to the masked genome. An output "is D" / "is K" when its best hit in the unmasked genome is
  >= 0.999 and overlaps D's / K's span.
- **Arms:** R = the existing masked-arm alignments; R+I = S_f realigned (same `@PG` command) to the masked genome + all contigs.
- **Final call per read (deviation from Amendment 5, declared):** the locus of its primary alignment; a read whose best and second-best AS
  over its records tie (second >= 0.98 x best) is "abstain". O2 is not run: in R every family has one copy left (O2 needs >= 2), and in
  Amendment 5 O2 assigned none of the tied reads in either arm.
- **Readouts, pooled over families:** D reads: right (an "is D" contig of its family) / wrong (K, another locus, or a contig that is not D)
  / unplaced (unmapped or abstain). K reads: false move = primary on any contig that is not "is K".
- **Rule (Amendment 5's, unchanged):** HELP if wrong D calls fall by >= 50% from R to R+I and false moves <= 5% of K reads; HURT if false
  moves > 10%; otherwise MIXED. Also reported: right D calls R -> R+I; per family; by fate (absorbed / orphaned) and by D-K divergence
  (median `de` of D reads on K in the masked arm, >= 0.01 vs < 0.01).

## Amendment 7 (2026-10-01): linking IsoCon transcripts to their source locus, held-out on multi-copy families (written before any run)

Amendment 6 showed that most added transcripts come from a family's own copies (alleles, ends), so reads tie between a locus and its own
transcript and abstain. The linking rule is fixed here from an independent source, then tested on deletions never used before.

- **Linking rule.** For each flagged IsoCon output, d = 1 - (matching bases of its best hit in the masked genome) / (output length).
  If d <= delta the output is an allele/variant of that hit's locus and is NOT added as a copy (its reads stay with the locus); otherwise
  it is a new copy and is added as a contig. **delta = 0.00958** = the 99th percentile of per-gene exonic divergence between KB3781's two
  haplotypes over 28,541 single-copy genes (>= 500 exonic bp, on both haplotypes; frozen truth `lift.tsv`), computed before this test.
- **Held-out substrate (new deletions).** Families of `o3_collapse/method/intervals/data/intervals.tsv` with >= 3 copies, every copy listed,
  >= 20 clean reads on every copy and every clean span <= 200 kb: **55 families, 207 copies**. In each family the copy last by (chrom,
  clean_start) is hard-masked to N over its clean interval; a family whose masked interval overlaps any other interval of the table is
  dropped (G3). Fibroblast Iso-Seq, KB3781.
- **Reads:** baseline primaries (the fibroblast BAM) on the masked copy (D reads) and on each surviving copy (S reads), up to 500 per copy,
  random, seed 1. Arms: R (masked genome), R+I (masked + every flagged output, as Amendment 6), R+I+L (masked + only the outputs the
  linking rule keeps as new copies). IsoCon input per family: the reads with an R-arm record overlapping a surviving copy, plus the
  reads R leaves unmapped; at most 1,000 (random, seed 1).
- **Final call:** the locus of the primary alignment (a contig is its own locus); abstain when best and second-best AS over the read's
  records tie (second >= 0.98 x best) and lie on different loci.
- **Labels:** a contig is D-derived when its best hit in the unmasked genome lies on D's interval (any identity). D read: right = a
  D-derived contig; wrong = anything else that is placed; unplaced = unmapped / abstain. S read: false move = placed on a contig that is
  not derived from its own copy (best unmasked hit on its own interval); stay = its own copy or a contig derived from it.
- **Rules (fixed now):**
  - **Linking works** if, from R+I to R+I+L, S reads unplaced fall by >= 50% AND D reads right fall by <= 20%.
  - **Overall (Amendment 5's rule, R vs R+I+L):** HELP if wrong D falls by >= 50% and false moves <= 5% of S reads; HURT if > 10%.
  - Reported beside: right D, unplaced, by D-to-nearest-survivor divergence (< 0.01 vs >= 0.01).

## Amendment 8 (2026-10-01): merging a missing copy's new-copy transcripts into one candidate copy, without truth (written before any run)

Amendment 7 leaves the several new-copy transcripts of one missing copy as separate contigs, so a deleted copy's reads tie among its own
transcripts and abstain (D right 1,616; counted as one locus WITH truth, 12,879). The merge rule is fixed here and scored on the same 53
families and the same R+I+L alignments: no realignment, a read's call changes only through which contigs count as one locus.

- **Merge rule.** Within a family, the new-copy contigs (`contigs_L.fa`, 565 in 46 families) are aligned all-vs-all
  (`minimap2 -c -x asm20 --cs -N 200 -p 0.1`, self hits dropped, the best alignment per unordered pair by matching bases). Two contigs are
  joined when that alignment covers >= 50% of the SHORTER contig's length AND its gap-compressed divergence (`de`) <= delta = 0.00958
  (Amendment 7's delta, unchanged). Components of the joined pairs (transitive) are the candidate copies; a contig joined to nothing is a
  candidate copy on its own. Gap-compressed divergence is used so that isoform differences (a skipped exon = one gap event) do not count as
  sequence divergence; the overlap floor stops two transcripts from being joined on a shared fragment.
- **Final call (Amendment 7's, with components as loci):** the locus of the primary alignment; a contig's locus is its component; abstain
  when records within 0.98 x best AS lie on different loci. Reads tied between two contigs of one component are placed in that component.
- **Labels (truth, scoring only).** A component is D-derived when it holds >= 1 D-derived contig (best unmasked hit on D's interval),
  S:g-derived when it holds >= 1 contig derived from surviving copy g, mixed when both. D read: right = placed on a D-derived component;
  wrong = any other placement; unplaced = unmapped / abstain. S read: stay = its own copy or a component holding a contig derived from its
  own copy; false move = any other contig / component; elsewhere = another reference locus.
- **Rules (fixed now):**
  - **M1 (merging works, read level):** from R+I+L to R+I+L+M, D right rises to >= 50% of Amendment 7's truth-grouped value (12,879,
    i.e. >= 6,440) AND false moves stay <= 5% of S reads. Reported beside it, Amendment 5's overall rule for R vs R+I+L+M (HELP if wrong D
    falls >= 50% and false moves <= 5%; HURT if > 10%).
  - **M2 (one missing copy, one candidate — the O3 count):** among the families with >= 2 D-derived new-copy contigs (41), the D-derived
    contigs fall in ONE component in >= 2/3 of them AND mixed components are <= 10% of all components.
  - Reported beside (not the verdict): components per family (total, D-derived, S-only = spurious new copies, mixed); the same counts at
    delta/2 and 2 x delta (sensitivity); the truth-grouped ceiling; the cause of each over-split (no alignment covering >= 50% vs
    divergence > delta).
- **Not tested here:** families with no missing copy (the false-flag rate without a deletion) — a separate control.
