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

## Amendment 9 (2026-10-01): the no-deletion control — flags raised when nothing is missing (written before any run)

Amendments 7-8 measured the chain where one copy per family is missing. Its specificity is measured here on the same families with
nothing masked: the chain runs against the full `_pri`, and every candidate copy it names is a flag raised without a deletion.

- **Substrate:** the 53 families of Amendment 7 and the same 59,013 scored reads (labels: family, copy), the unmasked `_pri`
  (`GGO.fasta`; splice index `winloci_data/GGO.splice.mmi`, the one the deletion run used for its unmasked arm).
- **Procedure, identical to Amendments 7-8 with `_pri` as the reference:** arm R0 = the scored reads aligned to `_pri` (the same minimap2
  command); IsoCon input per family = the reads with an R0 record on ANY copy of the family plus the reads R0 leaves unmapped, <= 1,000
  (random, seed 1); IsoCon; flag = output with identity x coverage < 0.999 against `_pri`; link = d <= delta (0.00958) to its best `_pri`
  hit -> allele of that locus, not a copy; merge = Amendment 8's components at delta -> candidate copies; arm C = `_pri` + the new-copy
  contigs, components as loci. "Derived from copy g" = best `_pri` hit overlapping g's clean interval (Amendment 7's "source").
- **Classification of each candidate against the diploid truth** (KB3781's own mat / pat assemblies; the component's contigs aligned with
  the same splice command; the component takes the class of its best contig): (a) **haplotype-only locus** — best hit at identity x
  coverage >= 0.999 on the haplotype `_pri` did NOT take that chromosome from, OUTSIDE the lifted B interval of every copy of the family
  (copy intervals lifted through the frozen truth's asm5 `_pri` -> B alignments, `out/chr*.paf`, primary records; a copy with < 50% of its
  interval lifted has no B interval) — a genuine reference-absent locus, a TRUE flag; (b) **allele** — best hit >= 0.999 inside a lifted
  copy interval of the family — a false flag (an allele beyond the 99th percentile); (c) **unmatched** — no hit >= 0.999 on either
  haplotype — a false flag (error / chimera).
- **Rules (fixed now):**
  - **C1 (specificity):** the fraction of the 53 families with >= 1 FALSE candidate (classes b + c) is <= 0.28 — one third of the
    deletion run's detection rate (44/53 = 0.83), i.e. a flag carries a likelihood ratio >= 3. Reported beside: the same fraction counting
    every candidate (a + b + c), the specificity against the haploid reference alone.
  - **C2 (cost without a deletion):** from R0 to C, reads placed on a candidate not derived from their own copy (false moves) <= 5% of
    all reads. Reported beside: reads that become unplaced, reads moving onto candidates derived from their own copy (harmless).
  - Reported: candidates per family by class; overlap with the deletion run's 21 survivor-derived candidates (same family and a contig at
    identity >= 0.999 to one of them); the counts at delta/2 and 2 x delta.
- **Not tested:** other individuals or tissues; the 1,000-read cap per family.

## Amendment 10 (2026-10-01): real reference-absent copies of KB3781 — does the chain flag them? (truth built first, rules written before the chain ran on these families)

**Truth (built independently of the chain; `bench/rna_allele/refabsent_truth.py`, work dir `/mnt/linuxdisk/tmp/rna_allele/refabsent`).**
Every copy of the 2026-08-14 interval table (378 families, 915 copies) was lifted to its B haplotype (asm5 alignments of the frozen truth;
913/915 lifted >= 50%) and its clean interval aligned to both haplotype assemblies (`minimap2 -c -x asm20 -p 0.1 -N 100`). A **B-only
locus** = a hit at identity >= 0.90 and coverage >= 0.80 on the haplotype `_pri` did NOT take that chromosome from, overlapping the lifted
B interval of no copy of the family (overlapping hits merged): 127 loci in 34 families. Each locus sequence was aligned back to `_pri`
(asm20): **absent beyond delta** (best `_pri` identity < 1 - 0.00958 over >= 50% of the locus; detectable in principle): 13 loci in 6
families; **absent within delta** (an allele of an unlisted `_pri` locus or a near-identical duplicate; indistinguishable from an allele
by construction): 114 loci in 33 families. **Expressed** = >= 3 fibroblast reads of the family whose best record over both haplotype
assemblies (AS, untied at 0.98) lies on the locus: **11 loci in 8 families — beyond delta: GWFAM175_B0 (281 reads), GWFAM26_B2 (13),
GWFAM26_B3 (9), GWFAM175_B1 (6); within delta: 7 loci (4-77 reads).** Two further beyond-delta loci have 2 reads (GWFAM175_B2,
GWFAM26_B6), below the floor. One locus has >= 20 reads: this is a demonstration on the biology of one individual and one tissue, not a
rate.

- **Chain:** Amendments 7-9's chain on the 34 B-only families (reads: every baseline record on a copy of the family, <= 2,000 per family,
  seed 1; IsoCon input <= 1,000); **flag = a candidate copy holding >= 2 IsoCon transcripts** (the support floor Amendment 9 found post
  hoc, pre-registered here); candidates classified against the haplotypes as in Amendment 9 (a / b / c).
- **Rules (fixed now):**
  - **D1 (detection):** every expressed beyond-delta locus with >= 20 reads (n = 1, GWFAM175_B0) receives a flag whose best haplotype hit
    at identity x coverage >= 0.999 overlaps the locus. Reported beside: the three beyond-delta loci with 6-13 reads, and the two with 2.
  - **D2 (delta does what it says):** none of the 7 expressed within-delta loci receives a new-copy flag (their transcripts link to a
    `_pri` locus); a flag there is reported as a false flag of the design.
  - **D3 (false flags in these families):** families among the 34 with a >= 2-transcript candidate of class b or c <= 20% (twice the
    control's 5/49 = 10.2% on families without any B-only locus).
  - Reported: for each expressed beyond-delta locus, where its reads sit in arm R0 (`_pri` copy, median `de`) and in arm C (the
    candidate), i.e. the O2 side; the single-transcript candidates (below the flag floor) matching a truth locus.

## Amendment 11 (2026-10-02): the in-house consensus in place of IsoCon, on the same held-out (written before any run)

The user wants no IsoCon dependency. The chain used IsoCon for one thing: reference-free consensus transcripts of a family's read net.
The in-house counterpart is `missing_copy_flag` (`src/rustle/vg_family/missing_copy.rs`): per locus, the per-read `de` mixture splits
the pile into host and divergent sub-pile, the sub-pile's mismatches are tested for PSV consistency, and a spliced consensus of the hidden
copy is built by patching the reference at those sites. This amendment swaps it in and scores it against IsoCon's own numbers.

- **Substrate:** Amendment 7's 53 families, the same masked genome, the same scored reads and the same R-arm alignments (`R.bam`).
- **In-house consensus:** `missing_copy_flag --bam R.bam --fasta masked.fa --loci surviving_copies.bed --index masked.splice.mmi` with
  defaults (`--delta-min 0.01 --min-reads 10 --min-sub 3 --max-reads 2000`), `--scan-only` then `--from-scan`. Loci = the 201 surviving
  copies' clean intervals (BED, name = family:copy). Its `consensus.fa` (one sequence per locus with a consistent divergent sub-pile)
  replaces IsoCon's outputs; everything downstream is Amendments 7-8 unchanged: flag (identity x coverage < 0.999 in the masked genome),
  link (d <= 0.00958), merge (components at delta), labels (D-derived = best unmasked hit on the masked interval), arm M = masked genome +
  the new-copy consensus contigs with components as loci.
- **Rules (fixed now), IsoCon's Amendment 8 result as the comparator (D right 12,787 = 74.0%, false moves 25 = 0.06%, one candidate per
  deleted copy in 33/41 families with >= 2 D contigs):**
  - **IH1 (adopt):** D right >= 80% of IsoCon's (>= 10,230) AND false moves <= 5% of S reads -> the in-house consensus replaces IsoCon in
    the `candidates` stage. Otherwise the gap is reported by cause (deleted copies with no consensus at all vs consensus present but reads
    not placed).
  - **IH2 (one copy, one candidate):** among families with >= 1 D-derived candidate, exactly one D-derived component in >= 2/3.
  - Reported beside: consensus count and `missing_copy_flag`'s own verdict classes; per deleted copy, whether a consensus exists and its
    d to the nearest survivor (copies closer than `--delta-min` 0.01 are below the S2 statistic's floor by construction, as they are
    below delta for the link rule); the 2-member floor's effect; structure (`n_bigins`, `n_rearr`) in the scan table.
- **Not tested here:** the structural extension (keeping clipped / inserted segments) — IsoCon's transcripts carry structure, the patched
  consensus does not; if IH1 fails, that is the first suspect, and the human DAZ3 panel (`yag_hsa`) is the structural check to run next.

## Amendment 12 (2026-10-02): `o3_candidates` (the in-house, IsoCon-free stage) on the same held-out (written before the stage ran on it)

Spec `docs/superpowers/specs/2026-10-02-o3-candidates-design.md`. Substrate: Amendment 7's 53 families, masked genome, the same
59,013 scored reads (`linktest/scored.fa`) given to the stage as the BAM `linktest/R.bam` (reads on the masked genome) with a copies
table made from `panel.json`'s surviving copies (clean intervals, `n_reads` from the BAM) in `P.fam.copies.tsv` format, and the
masked splice index. Arm M = masked genome + `P.cand.contigs.fa` (one union per flagged candidate), components as loci, scored by
`merge_test.py score` semantics (D right / wrong / unplaced; S false moves) with the contigs' D/S labels from their best unmasked hit.
- Flag = a component with >= 6 supporting reads over its clusters (`--min-support 6`, the >= 2-transcript floor translated: 2 x IsoCon's
  3-read minimum; the stage's clusters merge a copy's isoforms, so cluster counts cannot play that role). The >= 2-cluster count is reported beside.
- **A12-1 (adopt):** D right >= 80% of IsoCon's 12,787 (>= 10,230) AND false moves <= 5% of S reads.
- **A12-2 (representative):** >= 95% of the reads of each flagged component's clusters keep an AS on the union >= 0.98 x their best
  AS over the component's cluster consensuses (the `rep_choice.py` measure, run on the stage's own alignments).
- **A12-3 (cost):** wall time of the stage on the 53 families <= 40 min (2 x IsoCon's ~20 min) on this machine, 4 threads.
- Reported: candidates per family, clusters per candidate, the deleted copies with no candidate by cause (no reads in the net / clusters
  below the floor / linked to a survivor), the same numbers at delta/2 and 2 x delta, and the flag counts under the alternative >= 2-cluster floor.

## Amendment 13 (2026-10-03): the stage's net and template — the two causes Amendment 12 measured (written before any change runs)

Amendment 12 failed on two measured causes (`docs/O3_CANDIDATES_ACCEPTANCE_2026-10-02.md`): (1) the k-mer attribution placed 0 of 5,312
unmapped reads, so 25 of 53 deleted copies had no read in their family's net; (2) a longest-read template that retains an intron drags the
union (two-cluster unions kept 68.7% of reads). This amendment fixes how the stage builds those two things and re-runs Amendment 12's
rules unchanged. The chain's rules (delta link, component merge at delta, `--min-support 6`, 0.98 tie ratio) do not move.

- **Net attribution by alignment (replaces §5.2's k-mer rule for unmapped reads).** The unmapped records >= 300 bp of BAM pass B are
  written to a FASTA and aligned once against `--copies-fa` (every family's copy sequences) with `MM2_ATTRIB = [-x splice:hq -uf -c -N 5
  -p 0.5]`; a read joins the family of its best hit (most matches) when that hit covers >= 50% of the read and its `de` <= 0.15; no
  other read enters a net this way (reads whose primary lies at a locus outside every family and that carry no secondary on a family
  copy stay out of reach — stated, not fixed here). The k-mer index is retired from the net; `ATTRIB_*` constants go with it.
- **Template = the structurally central member (replaces "longest").** For each cluster, from the net's all-vs-all PAF (`--cs`), score
  every member by the total bases of indels >= 20 bp in its alignments to the other members (insertions and deletions alike: a retained
  intron shows as an insertion in every pair, a skipped exon as a deletion); the template is the member with the LOWEST score, ties
  broken by length (longest) then name. Clusters of one member keep that member.
- **Consensus details carried from the reviews (named so they are not silent):** the vote and refinement alignments use the splice preset
  (`MM2_MEMBERS` becomes `[-x splice:hq -uf -c --cs -N 5 -p 0.5]`; asm20 cut alignments at exon skips, as R6 found for the union); the
  insertion vote considers insertions >= 20 bp with >= 3 carriers before the < 20 bp majority rule at the same position; when refinement
  splits the template off, the kept set is re-templated by the same structural rule; an absorbing cluster whose re-polished consensus is
  empty keeps its absorbed clusters separate instead of dropping them.
- **Substrate and scoring: Amendment 12's, unchanged** — the 53 families, `linktest/R.bam`, the masked genome and index, `A12.copies.*`,
  arm M realigned with the pipeline's flags, `merge_test.py score` semantics, the A12-2 union measure with the arm-M preset, the summed
  batch wall time. **A13-1 = A12-1 (D right >= 10,230 and false moves <= 5%), A13-2 = A12-2 (>= 95% kept), A13-3 = A12-3 (<= 40 min).**
  Adopt (flip the `candidates` stage to default-on) iff all three hold. Reported beside: the same cause table as Amendment 12 (how many
  deleted copies reach their net now), the attribution counts (unmapped reads aligned / attributed / to the right family by label),
  and the delta/2 and 2 x delta reruns.
- **Not changed, by design:** `--min-support 6`, delta, the 0.98 tie ratio, the 1,000-read cap, R13 (reads placed uniquely on a candidate
  are the aligner's result).

### Amendment 13b (2026-10-03, before the A13 run): the attribution rule measured at the attribution step, and the comparator re-registered

Measured on the attribution step only (never on the outcome metric), with the 5,312 unmapped and the un-netted D reads of Amendment 12's
run aligned (`minimap2 -x map-ont -c`) against the families' own mapped net reads: the best hits of the deleted copies' reads to their
family's transcripts sit at `de` 0.15-0.30 (902 reads at coverage >= 0.5) or cover < 50% of the read (1,824); at the family definition's
own identity floor (clause 2: identity >= 0.80, coverage >= 0.50, i.e. `de <= 0.20`, read coverage >= 0.5) **135 of 5,312 unmapped
reads** (97.8% to the right family, 11 families) and **302 of 3,462 poorly placed un-netted reads** (92.1%, 14 families) are attributable;
Amendment 13's first draft (`splice:hq` against the copies, `de <= 0.15`) attributes 24. So: (1) no truth-free rule reaches the ~5,200
reads IsoCon received by truth label in Amendment 8 — the Amendment 12 bar (80% of IsoCon's 12,787) was partly unattainable by
construction, and it stands as recorded; (2) the rule is fixed as the family's own edge rule in read space.

- **Attribution set:** unmapped records >= 300 bp, PLUS mapped reads with no record on any family copy whose primary record has `de >
  0.02` or MAPQ 0 ("poorly placed"), both collected in BAM pass B. **Targets:** the families' mapped net reads (every read pass A put in a
  net, tagged `<family>|<read>`) together with `--copies-fa`. **Preset:** `MM2_ATTRIB = [-x map-ont -c -N 5 -p 0.5]`. **Rule:** a read
  joins the family of its best hit (most matches) iff the hit covers >= 50% of the READ and `de <= 0.20`. The chain's genome check
  (`InReference` at 0.999) remains the guard against reads of foreign genes pulled in this way.
- **Comparator re-registered for A13-1:** C = IsoCon's right D reads of Amendment 8 (`linktest/merge/score.out` semantics, per read)
  counted over the TRUTH-FREE ATTAINABLE D reads — D reads with any record on a surviving copy of their family in `R.bam`, or attributable
  by the rule above. **A13-1: stage D right (all reads) >= 0.80 x C AND false moves <= 5% of S reads.** A12-1's original bar (10,230) is
  reported beside, not decided on. A13-2 and A13-3 unchanged.
- Nothing else moves: delta, `--min-support 6`, 0.98, the 1,000-read cap, R13; the template rule of Amendment 13 stands.

#### Amendment 13c (2026-10-03, before the A13 run): two clarifications of 13b's attribution set. The 300-bp floor applies to BOTH classes
(unmapped and poorly placed reads). "No record on any family copy" is read as "in no net of this run" (pass A's scope under `--families`,
ruling R18), so a read whose only copy record is supplementary is eligible when its primary is poorly placed.

#### Amendment 13d (2026-10-03, before the A13 run): the template rule corrected at the step itself. Amendment 13's "lowest total bases of
indels >= 20 bp" favours fragments (a short read's alignments hold no indels) and, under the all-vs-all's `-N 100`, members with few aligned
partners (smoke on five A12 families: an 877-bp read, 48th longest of 50, became a cluster's template; a consensus never extends past its
template). The rule becomes the **medoid under a structural distance**: d(m, p) = bases of indels >= 20 bp (insertions, deletions and `~`
introns alike) in m's best alignment to p, PLUS p's terminal bases (>= 20 bp at either end) that m's alignment leaves uncovered; a member is
**eligible** when it is aligned to >= 50% of the cluster's other members (every member is eligible in clusters of < 4); template = the
eligible member with the lowest MEAN d over its aligned partners, ties -> longest -> smallest name; a cluster with no eligible member takes
its longest member. Unchanged: everything else in Amendments 13-13c.

#### Amendment 13e (2026-10-03, before the A13 run): eligibility made cap-aware. The all-vs-all keeps at most 100 hits per query (`-N 100`,
`--dual=no`), so in clusters of several hundred reads no member can be aligned to 50% of the others and eligibility followed read-name order
(smoke: 0% of the first name quartile vs ~50% of the last in 474- and 482-read clusters). Eligible = aligned to >= min(0.5 x (n - 1), 50)
other members. A member with no aligned partner has no mean and is never chosen while another member has one; if no member has a mean, the
longest is taken. Everything else as 13d.

## Amendment 14 (2026-10-03): the no-deletion control for the Rust stage, before its default-on ships (written before the run)

A13 passed its three rules (`docs/O3_CANDIDATES_ACCEPTANCE_A13_2026-10-03.md`) and the driver's `candidates` stage was flipped to
default-on in commit 1f49d0f0 — but the same run raised survivor-derived flags from 12 (A12) to 46 and left 3,527 survivor reads unplaced,
and the stage has never been run on Amendment 9's no-deletion control. The flip is NOT pushed until this control is measured.

- **Substrate and procedure: Amendment 9's, with the stage in place of the IsoCon chain** — the 53 families with nothing masked, the
  unmasked `_pri` (`control/R0.bam`, `winloci_data/GGO.splice.mmi`), a copies table of all 201 copies (`panel_to_copies.py` over `mask` +
  `keep`), the stage at its defaults (`--min-support 6`, batched as A13), every flagged candidate classified against the diploid truth as
  in Amendment 9 (a: haplotype-only locus outside every lifted copy interval = true flag; b: allele; c: unmatched = false flags), arm C =
  `_pri` + the flagged unions, components as loci.
- **Rules (fixed now):**
  - **C1' (specificity at the stage's operating point):** the fraction of the 53 families with >= 1 FALSE flag (b + c) <= 1/3 of A13's
    family-level detection rate 25/53 = 0.472, i.e. **<= 0.157 (<= 8 families)** — a flag carries a likelihood ratio >= 3, the bar
    Amendment 9 set. Reported beside: the rate counting every flag (a + b + c), and A9's 16/53 for the IsoCon chain at any support.
  - **C2' (cost without a deletion):** reads placed on a candidate not derived from their own copy <= 5% of all reads.
  - **Decision:** the default-on flip (1f49d0f0) ships iff C1' and C2' hold; otherwise it is reverted (the stage stays opt-in) and the
    control's numbers are the next prereg's starting point.
- Reported: candidates per family by class; the attribution counts without a deletion (how many unmapped / poorly placed reads join, and
  what they become); the per-family counters (refine splits, templates by kind).

## Amendment 15 (2026-10-03, written before any re-run): the consensus defect behind Amendment 14's false flags — the correction, the re-runs, a held-out

**What was found, post hoc (2026-10-03, after Amendment 14's verdict; `docs/O3_CANDIDATES_CONTROL_A14_2026-10-03.md`; Task 5 report in the
A13 ledger directory, reproduction byte-identical for 33 clusters of 5 families).** The 54 class-c flags of the control diverge from the
primary assembly by insertions only (no mismatches, no deletions; 99% of the inserted bases in runs >= 20 bp), and 41 of the 54 carry a
duplication signature: 94% of their genome-inserted bases (22,144 of 23,553) are copies of the union's own sequence (A13: 35 of the 46
survivor-derived flagged unions, 4 of the 30 D-derived; 3 of 159 sound linked consensus sequences). Mechanism (GWFAM37:c1, 368 reads):
the template is the medoid, a 3,007-bp read of the majority structure; 85 reads carry a 603-bp segment it lacks; `minimap2 -x splice:hq`
places each carrier's insertion at a column that depends on where the read starts, paired with a `~` over the template's own bases, so the
one segment appears as >= 3 identical carriers at SIX columns (618-747); `consensus_from_template` (Amendment 13's long-insertion rule:
the most frequent identical insertion >= 20 bp with >= 3 carriers, whatever its share of the covering members) inserts it at each, each
copy also re-inserting 71-200 bp of template that the `~` skipped and that a >= 20-bp deletion never removes (R2). This happens at the
first polish; the refinement rebuilds the same bytes (its fit test reads `de`, which is gap-compressed, and the read's own coverage, so
extra consensus sequence is invisible to it). The other 13 class-c unions have other causes (10 mismatch-dominated with few reads, median
7.5; 2 exact unions whose genome hit is split across ~204 kb).

**The correction (one rule made uniform; no new constant).** A long insertion (>= 20 bp) enters the consensus only when its carriers are a
majority of the members covering that column (`2 x count >= covering`, with `STRUCT_MIN_SUPPORT` = 3 kept as the floor) — the rule the
< 20 bp insertions, the substitutions and the deletions already obey. Measured offline before this amendment (Task 5): it clears all six
stage flags of the three worst control families (GWFAM37 8,201 bp / 0.485 -> 2,988 / 0.9997; GWFAM425 11,627 / 0.494 -> 3,529 / 0.9994;
GWFAM99 4,240 / 0.591 -> 1,894 / 0.9952) and adds none; two sound A13 consensus sequences stay sound and lose 128-216 bp (minority exons,
which the cluster consensus no longer carries; the candidate-level union across clusters is unchanged). The alternatives measured and NOT
taken: the longest eligible member as template (adds a false flag in GWFAM99; breaks A13 GWFAM37 to 0.856) and both together (adds the same
flag). Known residual, conservative by construction: an exon carried by a majority whose carriers minimap2 splits across columns is
under-counted and left out (fewer insertions, never duplications). Rulings R2/R5 and Amendment 13's vote are amended to this rule; delta,
the merge rule, the flag floor, the tie ratio, the net (13b/13c), the template (13d/13e) and the attribution rule are unchanged.

**Re-runs, same bars, nothing re-tuned.**
- A13 again (the 53-family deletion held-out, the same five batches, `ACC=a13`, the corrected binary): A13-1 (D right >= 0.80 x C, C
  re-measured from the run's own nets as Amendment 13b says; false moves <= 5%), A13-2 (union keeps >= 95%, pooled reading per R24),
  A13-3 (<= 40 min).
- A14 again (the same 53 families, nothing deleted, `ACC=a14`): C1' <= 8 of 53 families with a class b/c/pri flag; C2' false moves <= 5%.
- **Held-out, never run by the stage: Amendment 10's read set with nothing deleted** — inputs that exist: `refabsent/R0.bam` (Amendment
  10's 32,219 scored reads of 34 families, cap 2,000 per family, aligned to `_pri` with the pipeline flags) and the 915-copy table
  `a14/wholebam/W.copies.{tsv,fa}` restricted to those families; truth = Amendment 10's haplotype-only loci (`refabsent/bonly.tsv`, lift +
  asm20 to mat/pat), classification as Amendment 9. Four of the 34 (GWFAM4, GWFAM169, GWFAM175, GWFAM402) are among the 53 development
  families, so **the registered held-out is the 30 disjoint families**. H1: families with >= 1 class b/c/pri flag <= 15.7% of 30 = **at
  most 4** (the same LR >= 3 construction as C1' at A13's detection 25/53). H3: false moves <= 5% of all reads of the 30 families (arm C
  of Amendment 9 on this read set). Reported, not decided on (dev-overlapping, n = 1): H2, the expressed beyond-delta locus GWFAM175_B0 is
  flagged (Amendment 10's D1).
- **Decision:** the default-on flip is re-made iff A13-1/2/3, C1'/C2' and H1/H3 all hold; otherwise the stage stays opt-in (R14) and the
  outcome is recorded. Reported beside: the identity x coverage of every flagged union to its best genome hit (the number that exposed
  the defect) and the duplication signature count; the whole-BAM cost of the corrected stage on one 50-family batch (R23 repeated).
- Not yet scheduled (user's decision on cost, 2026-10-03 23:50): this amendment binds whenever the re-runs happen.
