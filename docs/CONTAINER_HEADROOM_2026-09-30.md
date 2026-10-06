# Container headroom under the current defaults: what a fused-locus container can still find (2026-09-30)

Diagnosis only. No rule was built, nothing in `src/`, `tools/` or `bench/` was edited, nothing was committed or pushed, and
no register row was appended (drafts with suffix H are in §9). Repository `main` at `e163d955`, clean at the start.

**Question.** Under the current default pipeline, how much headroom does a fused-locus "container" have for FINDING FAMILY
MEMBERS, and where do fused loci still cost members? "Container" is used in the two senses of the record: container v1
(`docs/PREREG_fusion_container_sim_2026-09-28.md`: a post-hoc core / accessory / relation view, it never changes a family),
and a UNIT-AWARE container (families found on units, the extra pieces of a locus kept as relations). The headroom of the second
is bounded here by an oracle (r1017's "perfect split", rebuilt for the current pipeline); the oracle reads the annotation and
can never ship.

**Conventions.**
- Sources are keys in brackets, listed in §10 with file paths; "computed here" marks a number derived from a listed product.
- **Species are never pooled.** Every count in §0, §2, §3, §6 and §7 is given per species (human: A119b chr16, A119b chr17, testis;
  gorilla: OR6737, KB3781) and the two species are never added or turned into one rate. A species sum adds *copy-by-cell
  observations*: the same annotated copy seen in two libraries counts twice (the 146 observations are 81 distinct annotated copies:
  26 human NPIP copies in two human libraries plus 16 human TBC1D3 records, 25 gorilla NPIP copies and 14 gorilla TBC1D3 records,
  each in two gorilla libraries). Libraries differ in depth, so a per-cell count is a statement about that library.
- **Current defaults** = the `e163d955` build in `rustle_target/release`: `assemble` runs `--bridge-regroup f1v2`, `families`
  runs `mcl_families --min-cov-shorter 0.70` (binaries copied to scratch, sha1 in §10; full command lines in §1). The pre-flip
  pipeline is `RUSTLE_BRIDGE_REGROUP=off RUSTLE_MIN_COV_SHORTER=0`.
- **Two instruments.** The HOLDER RULE (§1 Definitions: per truth copy, the locus with the most same-strand exon overlap) produces the
  fate table and every "copies in the family" count. The FAMILY SCORER (`family_score`: per truth gene, by gene name) produces every
  Compara / U2 / Soto / Liftoff number. They can disagree on the same copy (§3.7); every claim below names its instrument.
- **Nothing here is held out.** The two flips were decided on products that include these samples (F1v2 and COVER verdicts,
  09-29), and human chr16 is the development block of every fusion rule. These are measurements, not tests.
- **Exposure.** No chimp (PTR) or orangutan (PPY) family product was opened. One early directory listing of
  `rt_arms/chimp_PTR` showed file names only; no file in it was read.

## 0. Answer

1. **Per species, the membership headroom of a unit-aware container is 2 of 68 observations in human and 0 of 78 in gorilla**
   (holder-rule instrument; the oracle splits every fused locus into one node per constituent annotated gene; A and C agree
   where both were run). *Human* (chr16 NPIP 26, chr17 TBC1D3 16, testis NPIP 26 observations): under the current defaults 40 are
   in their family, 6 have a locus in another cluster, 22 have no locus. The oracle recovers **NPIPB5 on chr16** (24 → 25 of 26,
   sensitivity .923 → .962, F .857 → .909) and **LOC100420311 on chr17** (11 → 12 of 16), and nothing in testis (5 of 26
   before and after). *Gorilla* (OR6737 NPIP 25 and TBC1D3 14, KB3781 NPIP 25 and TBC1D3 14): 25 in their family, 1 in another
   cluster, 52 with no locus; the oracle recovers none (OR6737 10 of 25 and 7 of 14; KB3781 8 of 25 and 0 of 14). "In their
   family" is presence (the holder shares >= 1 same-strand exon base with the copy): the holder covers >= 50% of the copy
   in 34 of the 40 human and 9 of the 25 gorilla observations (§2.1). [fate, §2]
2. **The oracle does not recover the other five wrong-cluster observations (four human, one gorilla).** Human PKD1P6-NPIPP1 (chr16; a
   RefSeq readthrough record, not fused under any reading) and TBC1D3P5 (chr17; not fused) have every alignment to the family
   rejected by the edge gates (`no_exonic`, `low_shared`); human testis NPIPB4 and NPIPB5 have 529-bp single-exon holders with no
   alignment to any family member; gorilla OR6737 NPIPB7's holder is fused (LOC129527696 / 695 dominant) but overlaps 1.4% of the
   copy (touch only). [mech, §2.2]
3. **The 09-29 flip already moved 2 of the 4 copies that change fate in any arm.** Pre-flip, human chr16 had 22 of 26 copies in
   NPIP; under the current defaults 24: NPIPB2 and NPIPB6, whose holders are now `.rg2` pieces separated from GSPT1 and EIF3CL by
   the regroup of `--bridge-regroup f1v2` (no bridge was cut at those two loci, §2.4). The other two copies that change fate in any
   arm are the oracle's (NPIPB5, LOC100420311). Containment added no NPIP or TBC1D3 copy on any cell that has a containment-only
   arm (chr16 `Bc` = `B0`; OR6737 has no `Bc` arm, and D and Gc give the same counts there). Compara F (scorer instrument) on the human substrates: chr16
   .615 → .667, chr17 .412 → .475, the 19 A119b chromosomes .445 → .459 (pooled 2,746-pair truth, §3.4), testis .264 → .265.
   [sc16, sc17, sca19, sctes; §3]
4. **Genome-wide, fused loci cost few multi-copy genes on the holder rule; the family scorer loses more on the same substrate.** Of
   1,062 Compara multi-copy genes on the 19 A119b chromosomes (per-chromosome truth: a family with >= 2 genes on the chromosome), 400
   are in their matched cluster, 412 have no locus, 216 are unclustered and 34 are in another cluster; only **41 (3.9%)** are
   present, unplaced and in a fused locus (the same 41 under both fused readings). Human testis: 1 of 1,313 (family with >= 2
   genes anywhere). On chr16 the scorer finds more losses than the holder rule (CF153 misses 6 of the 19 genes the holder rule
   places, §3.7), so 3.9% is a holder-rule number. The oracle raises per-chromosome Compara F by +.047 on chr1 and +.054 (A) /
   +.042 (C) on chr19, the two chromosomes that hold 24 of the 41, and pooled-truth F over the 19 chromosomes from .459 to .472;
   part of the oracle's Compara gain is the scorer's own collision count (collapsed genes 3 → 1 → 0 on chr16, §4.2). [gf, sca19o; §4]
5. **The links that join two multi-copy families are in loci, not F1v2 bridges.** Loci joining two multi-copy Compara families: 3 on
   chr16, 15 on the 19 chromosomes (5 dominant, 5 mixed, 5 minority; all present in `families.gtf`), 3 in testis; Liftoff: 7 on the 19
   chromosomes, 1 in testis, 1 in gorilla OR6737, 0 in KB3781. F1v2 relation records joining two multi-copy families: 0 on the 19
   chromosomes and in testis under Compara or Liftoff, and 1 on chr16 under U2 (NPIPA6 | NPIPA7 | PKD1P1: one relation of two
   transcripts, 6 reads). [jn, §4.4]
6. **Container v1 puts the partner's bases in core blocks where they align to other family members, and leaves most of them in no
   block.** Twelve fused observations are already in their family (human 11, gorilla 1). The exon bases of their partner records:
   human chr16 113,007, of which 28,293 fall in core blocks of the holder unit, 8,327 in accessory blocks and 76,387 in no block (77% of
   the bases that are in a block are core); chr17 2,981 / 1,951 / 515 / 515; testis 4,194 / 846 / 61 / 3,287; gorilla OR6737 TBC1D3
   3,692 / 0 / 885 / 2,807. The container is empty (no accessory block) at NPIPA1, A6 and A9 (their PKD1P-side pieces align to other
   members, as the 09-28 Outcome found), at NPIPA7 (after F1v2's split its own locus, 11 transcripts and 184 reads, touches PKD1P2 with 4 link reads and
   every block of its unit is core: 53 of PKD1P2's 9,905 exon bases are in any block) and at TBC1D3G. [fate, §2.4]
7. **Simulation baseline under the current defaults** (prereg §3 measures, gorilla, S = 20 copies, 10 fused). Extended pool
   (realistic read ends; the pool of F1v2's dev numbers): fused copies back in NPIP at f = 0 / .1 / .5 / .9 / 1.0 = **9 / 9 / 5 /
   0 / 1**; pre-flip 9 / 7 / 1 / 0 / 1, so F1v2's dev numbers (9 / 5 / 0 at f = .1 / .5 / .9) are reproduced by the Rust default.
   Container-sim pool (error-free reads, 0-30 bp end jitter): 9 / 2 / 0 / 0 / 3; pre-flip 9 / 1 / 0 / 0 / 4. The current defaults
   fail the prereg bars (§5.3). Extended pool: M1 at f = .5 (5 against the f = 0 control 9), M2 at f = 1.0 (1 against 0), M3 at f = .5, .9, 1.0
   (precision / recall .19 / .12, .31 / .42, .63 / .44; bar .90 / .90). Container-sim pool: M1 at f = .1 and .5 (2 and 0 against 9), M2 at f = 1.0 (3), M3 at the same three f. M4 is 0 in every arm, by
   construction wherever the partner's holder is the fused member's own locus. [sim, §5]
8. **Ranked mechanisms, per species** (§6). Human, 28 observations not in their family: (e) no locus 22; (e) fragment 2; (d) admission gate
   2; (a) unit not separated 2; plus (c) 11 fused observations that are in their family, for which membership is unchanged.
   Gorilla, 53: (e) no locus 52; (e) touch-only holder 1; plus (c) 1. A unit-aware container changes (a) by at most +2 human
   observations (+1 of 26 on chr16, +1 of 16 on chr17) and nothing in gorilla; (b) placement by representative is suspected at NPIPB5 and
   was not tested; (d) and (e) are unchanged; (c) changes only as relations (on chr16 the oracle puts pieces of 12 to 14 of 83 fused loci in
   two or more families, 5 to 6 of them involving NPIP; an oracle count, scored against no truth). If every observation that has a locus
   were placed in its family (no container does that) membership would rise by 6 of 68 in human and 1 of 78 in gorilla. [mech, rel, §6, §7]

Independent check [ver]: a separate verifier with its own code recomputed from the raw products the fate counts of all seven cells (item 1 and the counts
in item 2), the Compara scores of chr16, chr17, the 19 chromosomes and testis under B0 and D (item 3), the loci joining two multi-copy Compara families
(item 5), the simulation's M1 and copies in NPIP (item 7), the fused-locus counts and the oracle GTF's integrity; every one reproduces. A second, hostile
read of the finished document [rev] re-derived the tables, the register rows and ten spot checks, and its findings are applied in this version. A third check [ver2]
recomputed from the raw products the claims that changed in this version (fate coverage, NPIPB5's bridge zone, testis MCL2, the container sums, the pooled truth, the per-gene Compara fates of §4.3, the simulation bars): 13 of 15 groups reproduced at once, and the
two numeric slips (the OR6737 median coverage, the exon count of the PKD1P6-NPIPP1 record) and the TBC1D3G zone label it found are corrected here. Not recomputed independently by anyone: the graph and alignment evidence of §2.2 beyond the degree of NPIPB5's unit, and the scorer runs themselves
(the `family_score` outputs were read, not re-run).

## 1. What was run

Every product below was produced by the current-default binaries unless the table says it was replayed or reused; nothing
was re-run that an earlier product reproduced byte for byte.

| cell | truth copies | products |
|---|---|---|
| human A119b chr16 (dev), NPIP | 26 Dishuck-checked copies (`truth.hsa.json`; NPIPB1P on chr18 excluded) | regional `copy_assign --assemble-only` on `chr16:0-96330374` with the stored genome-wide best-AS table, then `mcl_families --from-gtf` (commands below), direct minimap2 2.30 |
| human A119b chr17, TBC1D3 | 16 RefSeq records (9 protein-coding, 2 transcribed pseudogenes, 5 pseudogene spans without exon features) | same, chr17 |
| human testis, NPIP | the same 26 copies | the 09-29 `defaults_flip` genome-wide run (current defaults, all-vs-all through `tools/mm2_shard.sh`); chr16 alone re-run for the container and graph; the other arms by replay |
| gorilla OR6737 and KB3781, NPIP | 25 T_member copies (2026-09-17 proxy; NC_073241.2 1, NC_073242.2 23, NC_073244.2 1) | regional assemble of five contigs, then families on those five contigs through `tools/mm2_shard.sh` |
| gorilla OR6737 and KB3781, TBC1D3 | 14 RefSeq records (NC_073228.2 12, NC_073224.2 2) | same five contigs |
| 19 A119b chromosomes (chr1-12, 14, 15, 17, 19, X, Y, M) | Compara Primates, Liftoff pairs | the o1_cover per-chromosome GTFs (`F1v2.gtf` minus bridges = the families input) through the current binary with the stored all-vs-all replayed (below) |

**Commands** (copied from `tools/rustle_pipeline.sh`; scratch `code/asm.sh`, `code/fam.sh`).
Assemble: `RUSTLE_GTF_SECONDARY=1 RUSTLE_GTF_SECONDARY_AS_RATIO=0.98 RUSTLE_GTF_SECONDARY_AS_TABLE=<molecules table> copy_assign --assemble-only
--region CONTIG:0-LEN --assembly-junctions strict --assembly-polish full --polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82
--polish-ism-ratio 0.7 --polish-retained-ratio 10 --bridge-regroup f1v2 --gtf-tpm --bam BAM --fasta FASTA` (`--bridge-regroup off` in the pre-flip arms).
Families: `mcl_families --from-gtf GTF --fasta FASTA --threads 4 --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units --emit-container --dump-graph G
--min-cov-shorter 0.70` (0 in the pre-flip arms); `--emit-container` and `--dump-graph` only add outputs. **Annotations:** human
`/mnt/linuxdisk/tmp/regress/chm13.gff` is byte-identical (1,689,421,265 bytes, md5 691d676f) to the decompressed
`winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz` (RefSeq RS_2025_08; never `HSA_genomic.gff`); gorilla `winloci_data/GGO_genomic.gff`.

**Arms** (one families-input GTF and one `--min-cov-shorter` each):

| arm | families-input GTF | `--min-cov-shorter` | role |
|---|---|---|---|
| B0 | assembly with `--bridge-regroup off` | 0 | the pre-flip pipeline (09-25 products) |
| Gc | `f1v2` `families.gtf` (bridges removed) | 0 | F1v2 alone |
| Bc | `off` | 0.70 | containment alone (not run for OR6737) |
| **D** | `f1v2` `families.gtf` | 0.70 | **the current defaults** |
| OA | D's GTF, each fused locus replaced by one node per constituent record, each transcript clipped to the record's span (r1017 `oracle_a`) | 0.70 (and 0) | oracle |
| OC | as OA, each exon whole to the record it overlaps most; the piece spans its exons (r1017 `oracle_c`) | 0.70 (and 0) | oracle (not run for OR6737) |
| OS | r1017's selection (locus span contains >= 2 records, strand-blind, readthrough records included), clipped as OA | 0.70 (and 0) | oracle (chr16, chr17 only) |

**Aligner invocation (r1018).** The aligner call is inside `mcl_families --from-gtf` and is identical in every arm, the oracle
arms included: `minimap2 -x asm20 -c -X -N 50 -p 0.1 --secondary=yes -t 4 loci.fa loci.fa`. Only the GTF differs between an
arm and its baseline. The sharded wrapper `tools/mm2_shard.sh` (gorilla, testis) gave a byte-identical PAF to direct minimap2 on
testis chr16 (checked here: 1,333,010 bytes, `cmp` equal).

**Gates** (each passed before the product it covers was used; outputs in `out/gates.txt`) [gates]:
- The regional current-default GTF equals the stored genome-wide GTF (TPM and cov attributes stripped): A119b chr16 (F1v2.gtf, 91,283
  lines); gorilla OR6737 and KB3781 on all five contigs in both modes (`off` = stored BASE, `f1v2` = stored F1v2.gtf; 20 of 20).
  So the Rust default reproduces the Python F1v2 on these inputs.
- chr16 `B0` equals the stored dev BASE and `Gc` the stored dev CORE: `clusters.tsv`, `loci.paf`, `loci.tsv`, `copies.tsv` byte-identical.
- A119b per-chromosome replay: the stored ALL-run PAF with every record that has a bridge locus at either end removed is fed to
  the current binary through a stand-in `RUSTLE_MINIMAP2`. `Gc` equals the o1_cover COVER core and `B0` equals the stored BASE,
  `clusters.tsv` byte-identical on 19 of 19 chromosomes; on chr17 the replayed `D` and `Gc` equal the real minimap2 runs. Testis: replayed
  `D` equals the real genome-wide current-default run, `B0` the stored BASE, `Gc` the COVER core (byte-identical).
- Gorilla contig restriction: KB3781 `B0` on the five contigs vs the stored genome-wide BASE: 84 of 86 cluster sets identical on those contigs
  (the two that differ are on NC_073244.2 at 17.7-18.5 Mb and 27.7-29.0 Mb, not at NPIPB15, 21.1 Mb); the stored BASE NPIP and TBC1D3 clusters lie wholly inside the five contigs
  (OR MCL4, MCL34; KB MCL8). OR6737 `B0` is the stored genome-wide BASE restricted to the five contigs, not a re-run (flagged † in §2.1 and §3.1). Testis chr16 alone equals the genome-wide run on chr16 (20 of 20 cluster sets).
- Simulation: the pre-flip pipeline on the container-sim arms reproduces the stored Outcome (M1 9 / 1 / 0 / 0 / 4, M3 .3173 / .5328 at f = .1);
  on the extended arms it reproduces BASE (9 / 7 / 1 / 0 / 1; copies 23 / 22 / 16 / 15 / 16) and, with containment off, F1v2 (9 / 9 / 5 / 0 / 1; 23 / 23 / 19 / 15 / 16).
- The fused-locus instrument reproduces the F1v2 dev numbers on chr16: 83 (locus union) and 82 (per transcript).

**Definitions.**
- *Holder* of a copy: the locus (gene_id) with the most same-strand exon overlap with the copy's territory (ties: reads); `npf_audit`'s rule. A copy is
  *present* when its holder shares >= 1 same-strand exon base with it (an unstated 1-bp threshold: Appendix A gives the share of the copy the holder covers).
  *Fate*: right family (holder in the cluster that holds the most copies of the family), wrong family (clustered elsewhere), unclustered, no locus. A holder folded by
  `loci.tsv` counts in its absorber's cluster. `out/fate_table.tsv` also gives how much of the copy's territory the family-cluster loci cover and the clusters the copy touches, because
  the holder alone can hide a copy that sits in several clusters (NPIPB5: .056 in family-cluster loci).
- *Fused*, three columns P / L / R, each flagging a second annotated record at >= 1 bp of exon overlap: **P** per-transcript reading (one spliced transcript whose exons overlap two annotated gene / pseudogene
  records with disjoint spans, readthrough-described records excluded), **L** locus-union reading (the same on the union of the locus's exons; readthrough_eval `a.fused`, register 1147),
  **R** the audit reading (the holder carries a non-copy annotated record outside the copy's span, readthrough-described records included). R adds the four RefSeq
  readthrough copies that the instrument's gene set cannot see (NPIPA9 is R only).
- *Link* between a copy and a partner: L = reads of the holder's transcripts that touch both, C = reads of transcripts touching the copy and nothing else outside
  its span, P = reads of transcripts touching the partner and not the copy. **dominant** L >= C and L >= P; **minority** L < C and L < P; **mixed** otherwise.
  F1 status is read from the assembler's own tables. *F1v2 bridge (cut)*: a `bridges.tsv` row with share < 1/2 (its `keep` column is True: F1v2 KEEPS the junction as a bridge, that is, cuts the
  locus there and turns the transition transcripts into a `fusion_of` relation record). *F1 bridge candidate, not cut*: share >= 1/2 (`keep` False). *structural* with the proofs that
  hold (up / down; both shown when different junctions carry them, so no single junction is a bridge). This document says "cut" and "not cut", never "kept".
  A removed bridge is a `fusion_of` relation record and is not in the families input.
- *Container v1*: `--emit-container` (the Rust port, byte-identical to the frozen script): a block of the holder unit is core when one aligned CIGAR column joins it to an
  exon base of another member of the family, else accessory.

## 2. Fate table (deliverable 1)

### 2.1 Per cell and per species (current defaults D; B0 = pre-flip; oracle A / C)

| species | cell | copy-by-cell observations | right family | of which the holder covers >= 50% of the copy | wrong family | unclustered | no locus | right under B0 (pre-flip) | right under oracle A / C |
|---|---|---|---|---|---|---|---|---|---|
| human | human A119b chr16, NPIP (Dishuck 26) | 26 | **24** | 23 | 2 | 0 | 0 | 22 | 25 / 25 |
| human | human A119b chr17, TBC1D3 (16 records) | 16 | **11** | 10 | 2 | 0 | 3 | 11 | 12 / 12 |
| human | human testis, NPIP (Dishuck 26) | 26 | **5** | 1 | 2 | 0 | 19 | 5 | 5 / 5 |
| **human** | **sum of its 3 cells (not a rate)** | **68** | **40** | **34** | **6** | **0** | **22** | **38** | **A: 42** |
| gorilla | gorilla OR6737, NPIP (T_member 25) | 25 | **10** | 0 | 1 | 0 | 14 | 10 (†) | 10 / - |
| gorilla | gorilla OR6737, TBC1D3 (14 records) | 14 | **7** | 7 | 0 | 0 | 7 | 7 (†) | 7 / - |
| gorilla | gorilla KB3781, NPIP (T_member 25) | 25 | **8** | 2 | 0 | 0 | 17 | 8 | 8 / 8 |
| gorilla | gorilla KB3781, TBC1D3 (14 records) | 14 | **0** | 0 | 0 | 0 | 14 | 0 | 0 / 0 |
| **gorilla** | **sum of its 4 cells (not a rate)** | **78** | **25** | **9** | **1** | **0** | **52** | **25** | **A: 25** |

Source: [fate]. Only four copies change fate among B0, D and the oracles: NPIPB2 and NPIPB6 (B0 → D, both chr16) and NPIPB5 (chr16) and LOC100420311 (chr17) (D → oracle). (†) OR6737's B0 is the stored genome-wide BASE
restricted to the five contigs, not a five-contig re-run; KB3781's B0 is a re-run (84 of 86 cluster sets equal to the stored BASE). The full 146-row table is Appendix A.

"Right family" is presence, not recovery: the fifth column counts the right-family observations whose holder covers >= 50% of the copy. Human 34 of 40 (chr16 23 of 24, chr17 10 of 11, testis 1 of 5);
gorilla 9 of 25 (OR6737 NPIP 0 of 10, median holder coverage .23; OR6737 TBC1D3 7 of 7; KB3781 NPIP 2 of 8, minimum .115). The oracle comparison is unaffected (it asks whether the holder's cluster changes), but the absolute
"in family" counts are presence counts.

### 2.2 The seven observations with a locus outside their family (six human and one gorilla; the two testis observations share a row)

| cell | copy | holder locus (gene_id; node; transcripts / reads; share of territory) | fused P / L / R, link | representative | cluster it lands in, and what that cluster is | graph and alignment evidence | oracle A / C | container v1 for its unit | mechanism |
|---|---|---|---|---|---|---|---|---|---|
| human chr16 | NPIPB5 | `DN_chr16_22714381_13`; chr16:22,696,405-22,854,661 (158 kb); 66 tx / 503 reads; holder covers .658 of the copy, family-cluster (MCL1) loci cover .056; the territory is also touched by five MCL3 loci (1.4-2.5 kb of overlap each) and by MCL12 (385 bp) | P / L / R; SMG1P1 link L4 / C32 / P456 (minority, 2 transcripts). The locus's one F1 bridge candidate (chr16:22,753,211-22,756,690; 24 transition transcripts, 188 reads, UP 166, DOWN 149, share .5579 >= 1/2, **not cut**) lies inside SMG1P1 (22,714,642-22,770,227), 25 kb upstream of NPIPB5 (22,781,607-22,814,310): it is not the NPIPB5 link, and NPIPB5's own link has no bridge candidate (its shared junction row, 22,801,954-22,805,581, has 6 transcripts and neither proof). OTOAP1 minority L2 / C32 / P11, no proof | SMG1P1 transcript, 73 reads (NPIPB5's own transcripts carry 32 reads in the locus) | MCL26 (4 units: SMG1, SMG1P6, the SMG1P4 locus, this): the SMG1P family | degree 5 over four clusters, a hub: 2 admitted edges into MCL1 (w .991, .990), 1 into MCL26, 1 into MCL84, 1 folded into MCL10; the locus aligns end to end (117,861 bp at 99.2%) to the NPIPB3~SMG1P3 locus and (139,925 bp) to the LOC128966608~SMG1-like locus, both in MCL1 | right / right | acc 17,886 bp in 23 blocks, relations to 12 families including MCL1 (NPIP); partner bases core 4,389 / acc 3,825 | (a); hub over four clusters; (b) suspected, untested |
| human chr17 | LOC100420311 | `DN_chr17_31496546_13`; chr17:31,496,547-31,513,862; 21 tx / 72 reads; holder covers .266 | P / L / R; TBC1D29P dominant L59 / C0 / P11, structural without proofs | the fusion transcript (15 reads) | MCL62 (2 units: TBC1D3P5 and this): a second, small cluster of two of the 16 TBC1D3 records; the TBC1D3 family cluster is MCL3 (11 units) | degree 1 (to TBC1D3P5, w .862); 24 alignments to the TBC1D3 family, all 11 evaluated pairs fail `low_shared` (the shared-exon fraction is not recorded, so the distance from .60 is unknown) | right / right | acc 2,854 bp, partner core 666 / acc 0, relation to MCL28 only, none to the TBC1D3 family | (a); `low_shared` is the proximate gate |
| human chr16 | PKD1P6-NPIPP1 | `DN_chr16_15127017_8`; chr16:15,127,011-15,137,722; 9 tx / 27 reads; holder covers .541 | not fused | the copy (5 reads) | MCL25 (3 units: this, a PKD1P2 piece, an NPIPA8 / PKD1P4-NPIPA8 piece): a PKD1P-homologous cluster | degree 3, all inside MCL25; 22 alignments to 3 NPIP members (the PKD1P-homologous part), all `no_exonic`; 27 of the family's 30 units have no alignment record | wrong / wrong | acc 2,409 bp, relations to MCL1 and MCL26 | (d); open: the territory chr16:15,126,651-15,141,806 lies inside the RefSeq PKD1P6 record (15,126,099-15,159,720), so a truth-labelling question is not excluded (§8) |
| human chr17 | TBC1D3P5 | `DN_chr17_28369204_2`; chr17:28,359,452-28,373,096; 12 tx / 34 reads; holder covers 1.0 | not fused | the copy (6 reads) | MCL62 (with LOC100420311): the small TBC1D3 cluster above | degree 1; 26 alignments to the TBC1D3 family, all 11 evaluated pairs fail `low_shared` (fraction not recorded) | wrong / wrong | acc 820 bp, no relation | (d) |
| human testis | NPIPB4 and NPIPB5 | `DN_chr16_22381847_1` and `DN_chr16_22813781_1`; 529-bp single-exon loci, 1 tx / 18 reads each; holders cover .142 and .057 | not fused | the copy | MCL2: 7 single-transcript `+`-strand loci of 18-19 reads (529 or 1,274 bp), each overlapping an NPIP copy: two on the copy's own strand (NPIPB4 515 bp, NPIPB5 529 bp), five antisense to it (NPIPB3, LOC128966608, NPIPB12, LOC124907834, NPIPB13); not shown to be NPIP transcripts | degree 6, all inside MCL2; **no alignment record** to any of the 5 units of the NPIP family (MCL3) | wrong / wrong | none | (e) fragment |
| gorilla OR6737 | NPIPB7 | `DN_NC_073242.2_99251348_15`; NC_073242.2:99,251,349-99,267,592; 2 tx / 4 reads; holder covers .014 | P / L / R; LOC129527696 / 695 dominant L4 / C0 / P0 | the fusion transcript (2 reads) | MCL6 (8 units, the SMG1-like family) | degree 7, all inside MCL6; no alignment record to the NPIP family's 14 units | A: wrong (MCL3, 13 units) at cov .70, unclustered at cov 0; C not run | partner core 1,729 | (e) touch-only holder (fused; the copy is touched only by the two fusion transcripts, 1.4% of its territory) |

Source: [fate, mech, reads]. U2 and Soto file PKD1P6-NPIPP1 in ID_149 (the PKD1 block). Under D the scorer's matched cluster for ID_149 is MCL23 (hits PKD1, PKD1P2, PKD1P4-NPIPA8), not MCL25, and PKD1P6-NPIPP1 is not among the
hits under U2 either; Dishuck calls it NPIP (the truths disagree on the fused members, as `npf_critique` §5.2 found). Hub: NPIPB5's locus is the only holder of the seven with admitted edges into more than one cluster; the other six have all
their edges inside their own cluster (degree 1 to 7), so hubness is not the mechanism for them. The link transcripts of NPIPB5 are two (by start coordinate, DN_chr16_22752169_21 and DN_chr16_22761474_23, both starting inside SMG1P1); only the first is among the 24 transition transcripts of the bridge candidate.

### 2.3 The 74 observations with no locus (human 22, gorilla 52)

| cell | copies | no locus | of which < 2 same-strand primary reads on the territory | of which >= 2 reads (none assembled) |
|---|---|---|---|---|
| human A119b chr16, NPIP (Dishuck 26) | 26 | 0 | 0 | 0 |
| human A119b chr17, TBC1D3 (16 records) | 16 | 3 | 0 | 3 |
| human testis, NPIP (Dishuck 26) | 26 | 19 | 10 | 9 |
| gorilla OR6737, NPIP (T_member 25) | 25 | 14 | 2 | 12 |
| gorilla OR6737, TBC1D3 (14 records) | 14 | 7 | 4 | 3 |
| gorilla KB3781, NPIP (T_member 25) | 25 | 17 | 2 | 15 |
| gorilla KB3781, TBC1D3 (14 records) | 14 | 14 | 14 | 0 |

Source: [reads] (primary alignments, flag & 2308 = 0, with >= 1 aligned base on the territory exons; the same-strand count is the gate because the assembler's pass-1 floor is 2
reads; same strand = the read's transcript strand, the minimap2 `ts` tag combined with the alignment orientation, or the alignment strand when there is no tag, equals the copy's strand; on the plain alignment strand instead, the human split would be 14 with >= 2 reads and 8 with < 2, not 12 and 10). 12 of the 74 have a locus on the opposite strand only (human testis 5; gorilla OR6737 NPIP 2, KB3781 NPIP 5); they count as no locus because a copy is scored on its own strand. A first losing step is cited, not recomputed, for 24 of the 74
[CRT]: human chr17 3 (TBC1D3P4 mono shadow, TBC1D3P3 gate, TBC1D3P7 mono floor); gorilla OR6737 NPIP 14 = pass-1 floor 3 + gate 8 + mono floor 3; gorilla OR6737 TBC1D3 7 = pass-1 floor 4 + seeding 3; none of these steps is a fusion step. The other 50
(human testis 19, gorilla KB3781 31) have only read counts, and the 42 with >= 2 same-strand primary reads are labelled floor / gate without a step; the [CRT] pass-1-floor count for OR6737 NPIP (3) differs by one from the read-count split here (2 with < 2 reads), not reconciled. The
assembler seeds GOOD secondaries (register 1192), so primary counts do not separate "no reads" from "gate". Fusion involvement was not examined for any of the 74; there is no locus for a container to split, so it cannot reach them.

### 2.4 Fused copies, the link, and container v1

| cell | copy | fused P / L / R | partners and link (partner: class, link reads L / copy-side reads C / partner-side reads P; F1 status) | F1v2 bridge relations at this locus | fate |
|---|---|---|---|---|---|
| human chr16 | NPIPA1 | P / L / R | PKD1P3-NPIPA1:dominant L193/C2/P17 [struct up/-]; PKD1P3:dominant L8/C2/P2 [struct up/-] | - | right family |
| human chr16 | NPIPA6 | P / L / R | PKD1P1:dominant L401/C128/P0 [no row]; LOC131696449:dominant L401/C128/P0 [no row]; MIR6511A2:mixed L3/C128/P0 [no row] | 2 | right family |
| human chr16 | NPIPA7 | P / L / R | PKD1P2:mixed L4/C180/P0 [struct -/down] | 2 | right family |
| human chr16 | NPIPA9 | - / - / R | PKD1P5-LOC105376752:dominant L403/C126/P0 [no row] | - | right family |
| human chr16 | NPIPB3 | P / L / R | LOC100190986:mixed L10/C357/P2 [struct -/-]; SMG1P3:minority L2/C357/P466 [struct up/down] | - | right family |
| human chr16 | LOC128966608 | P / L / R | LOC124900576:minority L2/C134/P8 [struct up/down]; LOC128966632:minority L67/C134/P256 [struct up/down] | - | right family |
| human chr16 | NPIPB4 | P / L / R | RRN3P1:minority L2/C351/P131 [struct -/down] | - | right family |
| human chr16 | NPIPB5 | P / L / R | SMG1P1:minority L4/C32/P456 [F1 bridge zone chr16:22,753,211-22,756,690 lies INSIDE the partner record (not at the copy boundary), not cut (share 0.5579 >= 1/2)]; OTOAP1:minority L2/C32/P11 [struct -/-] | - | wrong family |
| human chr16 | NPIPB14P | P / L / R | PDXDC2P-NPIPB14P:dominant L94/C0/P8 [struct up/down]; PDXDC2P:dominant L57/C0/P8 [struct up/down] | - | right family |
| human chr17 | LOC100420311 | P / L / R | TBC1D29P:dominant L59/C0/P11 [struct -/-] | - | wrong family |
| human chr17 | TBC1D3G | P / L / R | LOC101060212:mixed L53/C237/P0 [F1 bridge zone chr17:37,268,487-37,273,349 lies INSIDE the copy's span (an intron of the copy, not at the copy-partner boundary), not cut (share 0.6226 >= 1/2)] | - | right family |
| human chr17 | TBC1D3D | P / L / - | none | - | right family |
| human chr17 | TBC1D3 | P / L / R | NPEPPSP1:mixed L54/C231/P44 [F1 bridge zone chr17:39,054,087-39,058,950 spans the gap between the copy and the partner, not cut (share 0.5510 >= 1/2)] | - | right family |
| human testis | NPIPB14P | - / - / R | PDXDC2P-NPIPB14P:dominant L2/C0/P0 [no row] | - | right family |
| gorilla OR NPIP | NPIPB7 | P / L / R | LOC129527696:dominant L4/C0/P0 [no row]; LOC129527695:dominant L4/C0/P0 [no row] | - | wrong family |
| gorilla OR TBC1D3 | LOC129533797 | P / L / R | LOC115933338:mixed L2/C7/P0 [no row]; LOC115933339:dominant L15/C7/P0 [no row] | - | right family |

Source: [fate]. "struct up / down" = a structural junction row exists, the UP proof holds at one junction and the DOWN proof at another, none at both, so the link is
not an F1 bridge. "F1 bridge zone ... not cut" = the link transcripts are transition transcripts of a bridge zone with read share >= 1/2, so F1v2 leaves the locus whole; the label says where the zone lies. Only TBC1D3's zone (chr17:39,054,087-39,058,950) spans the gap between the copy and
the partner (NPEPPSP1); TBC1D3G's (chr17:37,268,487-37,273,349) lies inside the copy's own span, an intron of the copy 526 bp inside it, and NPIPB5's lies inside the partner record (SMG1P1). Across the 146 observations' loci, F1v2's cut bridges (`fusion_of` records) touch two loci only, NPIPA6 and NPIPA7 (one relation of two transcripts, 6 reads),
which is why those two are now separate loci. NPIPB2 and NPIPB6 were separated by the regroup into connected pieces (`.rg2`), not by a bridge: under B0 their holders were unions of disjoint transcript sets (link reads L = 0 at both; GSPT1's 925 reads against the copy's 23 at NPIPB2, EIF3CL's 697 against
207 at NPIPB6; register 1191 names the same NPIPB2 / GSPT1 fusion) [fate]. TBC1D3D's locus carries transcripts that touch TBC1D3E (another truth copy, L7 / C241; TBC1D3E has its own unfused holder): a copy-copy fusion, counted in the P and L tallies and not a partner.

Container v1 on the fused observations that are already in their family (R reading; a fused copy whose partner bases fall in no block at all would be omitted, none is) [fate]:

| cell | fused copies already in the family | partner records: exon bases | in core blocks | in accessory blocks (the container) | in no block of the holder unit | core share of the bases that are in a block | copies with no accessory block |
|---|---|---|---|---|---|---|---|
| human A119b chr16, NPIP (Dishuck 26) | 8 (NPIPA1, NPIPA6, NPIPA7, NPIPA9, NPIPB3, LOC128966608, NPIPB4, NPIPB14P) | 113,007 | 28,293 | 8,327 | 76,387 | 77.3% | 4 (NPIPA1, NPIPA6, NPIPA7, NPIPA9) |
| human A119b chr17, TBC1D3 (16 records) | 2 (TBC1D3G, TBC1D3) | 2,981 | 1,951 | 515 | 515 | 79.1% | 1 (TBC1D3G) |
| human testis, NPIP (Dishuck 26) | 1 (NPIPB14P) | 4,194 | 846 | 61 | 3,287 | 93.3% | 0 |
| gorilla OR6737, TBC1D3 (14 records) | 1 (LOC129533797) | 3,692 | 0 | 885 | 2,807 | 0.0% | 0 |

Per copy on human chr16:

| copy (human chr16) | partner record(s) | partner exon bases | in core | in accessory | in no block of the unit | accessory bp of the unit (all bases) |
|---|---|---|---|---|---|---|
| NPIPA1 | PKD1P3-NPIPA1, PKD1P3 | 7,508 | 3,715 | 0 | 3,793 | 0 |
| NPIPA6 | PKD1P1, LOC131696449, MIR6511A2 | 8,311 | 4,473 | 0 | 3,838 | 0 |
| NPIPA7 | PKD1P2 | 9,905 | 53 | 0 | 9,852 | 0 |
| NPIPA9 | PKD1P5-LOC105376752 | 8,111 | 3,841 | 0 | 4,270 | 0 |
| NPIPB3 | LOC100190986, SMG1P3, SLC7A5P2 | 10,047 | 7,316 | 617 | 2,114 | 2,040 |
| LOC128966608 | LOC124900576, LOC128966632, LOC124905420, LOC128966680 | 9,433 | 6,269 | 1,401 | 1,763 | 2,830 |
| NPIPB4 | RRN3P1 | 1,609 | 0 | 1,609 | 0 | 10,044 |
| NPIPB14P | PDXDC2P-NPIPB14P, PDXDC2P | 58,083 | 2,626 | 4,700 | 50,757 | 5,178 |

Reading. (i) The partner bases that container v1 puts in accessory blocks are few: 8,327 of 113,007 on chr16 (7.4%), 515 of 2,981 on chr17, 61 of 4,194 in testis; most partner bases (76,387 on chr16) are in no block of the holder unit (the unit's
blocks are built from the holder's own transcripts; NPIPB14P's partner record alone has 58,083 exon bases). (ii) Of the partner bases that do fall in a block, 77.3% (chr16) are core: the partner piece aligns to another family member, as the PKD1P pieces of NPIPA1, A6 and A9 (3,715,
4,473 and 3,841 bases core, none accessory) and the SMG1P pieces of NPIPB3 and LOC128966608 (7,316 and 6,269 core) do. Core is defined by alignment to another member, so this is by construction; whether the partner is co-duplicated inside the family was not measured independently.
(iii) The container is informative where the partner does not align to a family member (NPIPB4~RRN3P1: all 1,609 bases accessory; gorilla LOC129533797: 885 bases accessory), but the accessory set is mostly not partner bases: NPIPB4's container is 10,044 bp of which 1,609 (16%) are
partner bases; NPIPB3 617 of 2,040; LOC128966608 1,401 of 2,830; NPIPB14P 4,700 of 5,178 (91%); LOC129533797 885 of 1,603. (iv) NPIPA7's container is empty for a different reason than NPIPA1, A6 and A9: under B0 NPIPA6 and NPIPA7 were one locus (87 transcripts, 719 reads; container empty, all core); F1v2 split it into NPIPA6's piece (`.rg2`, 74 transcripts, 529 reads, which keeps the PKD1P1 link, L401) and
NPIPA7's own locus (11 transcripts, 184 reads), which touches PKD1P2 with 4 link reads, so only 53 of PKD1P2's 9,905 exon bases lie in any block of its unit and all of its blocks are core. The 09-28 Outcome (`PREREG_fusion_container_sim`, Outcome, human descriptive) found NPIPA1, A6 and A9 empty and all core; that holds on the current defaults.

## 3. Family scores (deliverable 2)

Instrument: `family_score` (copied, sha1 in §10) exactly as `o1_cover_frozen/score.py` calls it (`--chrom ALL --pairwise`, clusters and truth restricted to the substrate's
contigs first). Truths: Compara Primates families (human); U2, the NPIP union truth (chr16); Soto 2025, **labelled not independent**; the Compara chr16 paralogue table
`compara_chr16.tsv` (pair level; its same-chromosome pairs, 1,324 on chr16, are scored as pairs whose two genes are both in some cluster, and it reproduces the family
scorer's true-pair counts exactly); Liftoff copy pairs (both species). **RefSeq gene families:** there is no genome-wide RefSeq family truth in this stack and none was scored (no description-stem family was used because no earlier test did); the only
RefSeq-name family scored is the TBC1D3 record set (16 human records, 14 gorilla records), scored per copy in §3.1 where the truth is the record list itself. For the gorilla samples no family truth exists except the
T_member and TBC1D3 copy sets, so they are scored per copy (§3.1). The Dishuck NPIP set is the human NPIP truth of §3.1. For a one-family truth the matched-pair F of the bipartite score reduces to the harmonic mean of sensitivity and unit precision printed in §3.1; the family scorer's own bipartite numbers
are the Compara / U2 / Soto columns of §3.2-§3.5.

Three traps of this scorer apply to every Compara / U2 / Soto row below and were not removed: (1) predictions are intersected with the truth universe before precision is taken (an unlabelled real member is invisible, so
precision is inflated, most for the arm with the most unlabelled members; register 991 / 992); (2) one-to-one bipartite precision depends on scipy's assignment tie policy (register 1045), the pairwise true-positive /
predicted counts are tie-free and are printed beside every bipartite number; (3) the scorer resolves each truth gene to ONE locus by name, so in a locus that spans two genes the other gene is invisible (the 09-21 one-name-per-locus trap; `collapsed` counts these). Differences between
arms of the same substrate share traps (1) and (2) for the non-oracle arms; they do NOT for the oracle arms, because the oracle changes the node set, hence the predicted universe (chr16 Soto truth genes with no locus: D 17, oracle A 18, C 17, S 19; collapsed Compara genes 3 → 1 / 1 / 0) and names its nodes from the
same annotation the truths derive from. Oracle numbers are upper bounds, not measurements of family finding.

### 3.1 The NPIP and TBC1D3 families, per copy (sensitivity = copies in the family cluster / copies; unit precision = family-cluster units that hold a copy / units)

| cell | arm | copies in the family cluster | family cluster units (holding a copy) | sensitivity | unit precision | F (harmonic mean of the two) |
|---|---|---|---|---|---|---|
| human chr16 NPIP | B0 pre-flip | 22 / 26 | 27 (21) | 0.846 | 0.778 | 0.811 |
| human chr16 NPIP | Gc F1v2 only | 24 / 26 | 30 (24) | 0.923 | 0.800 | 0.857 |
| human chr16 NPIP | Bc containment only | 22 / 26 | 27 (21) | 0.846 | 0.778 | 0.811 |
| human chr16 NPIP | **D current** | 24 / 26 | 30 (24) | 0.923 | 0.800 | 0.857 |
| human chr16 NPIP | oracle A | 25 / 26 | 29 (25) | 0.962 | 0.862 | 0.909 |
| human chr16 NPIP | oracle C | 25 / 26 | 29 (25) | 0.962 | 0.862 | 0.909 |
| human chr17 TBC1D3 | B0 pre-flip | 11 / 16 | 11 (11) | 0.688 | 1.000 | 0.815 |
| human chr17 TBC1D3 | Gc F1v2 only | 11 / 16 | 11 (11) | 0.688 | 1.000 | 0.815 |
| human chr17 TBC1D3 | Bc containment only | 11 / 16 | 11 (11) | 0.688 | 1.000 | 0.815 |
| human chr17 TBC1D3 | **D current** | 11 / 16 | 11 (11) | 0.688 | 1.000 | 0.815 |
| human chr17 TBC1D3 | oracle A | 12 / 16 | 12 (12) | 0.750 | 1.000 | 0.857 |
| human chr17 TBC1D3 | oracle C | 12 / 16 | 13 (12) | 0.750 | 0.923 | 0.828 |
| human testis NPIP | B0 pre-flip | 5 / 26 | 5 (5) | 0.192 | 1.000 | 0.323 |
| human testis NPIP | Gc F1v2 only | 5 / 26 | 5 (5) | 0.192 | 1.000 | 0.323 |
| human testis NPIP | Bc containment only | 5 / 26 | 5 (5) | 0.192 | 1.000 | 0.323 |
| human testis NPIP | **D current** | 5 / 26 | 5 (5) | 0.192 | 1.000 | 0.323 |
| human testis NPIP | oracle A | 5 / 26 | 5 (5) | 0.192 | 1.000 | 0.323 |
| human testis NPIP | oracle C | 5 / 26 | 5 (5) | 0.192 | 1.000 | 0.323 |
| gorilla OR6737 NPIP | B0 pre-flip (†) | 10 / 25 | 14 (10) | 0.400 | 0.714 | 0.513 |
| gorilla OR6737 NPIP | Gc F1v2 only | 10 / 25 | 14 (10) | 0.400 | 0.714 | 0.513 |
| gorilla OR6737 NPIP | **D current** | 10 / 25 | 14 (10) | 0.400 | 0.714 | 0.513 |
| gorilla OR6737 NPIP | oracle A | 10 / 25 | 14 (10) | 0.400 | 0.714 | 0.513 |
| gorilla OR6737 TBC1D3 | B0 pre-flip (†) | 7 / 14 | 7 (7) | 0.500 | 1.000 | 0.667 |
| gorilla OR6737 TBC1D3 | Gc F1v2 only | 7 / 14 | 7 (7) | 0.500 | 1.000 | 0.667 |
| gorilla OR6737 TBC1D3 | **D current** | 7 / 14 | 7 (7) | 0.500 | 1.000 | 0.667 |
| gorilla OR6737 TBC1D3 | oracle A | 7 / 14 | 7 (7) | 0.500 | 1.000 | 0.667 |
| gorilla KB3781 NPIP | B0 pre-flip | 8 / 25 | 10 (8) | 0.320 | 0.800 | 0.457 |
| gorilla KB3781 NPIP | Gc F1v2 only | 8 / 25 | 10 (8) | 0.320 | 0.800 | 0.457 |
| gorilla KB3781 NPIP | Bc containment only | 8 / 25 | 10 (8) | 0.320 | 0.800 | 0.457 |
| gorilla KB3781 NPIP | **D current** | 8 / 25 | 10 (8) | 0.320 | 0.800 | 0.457 |
| gorilla KB3781 NPIP | oracle A | 8 / 25 | 10 (8) | 0.320 | 0.800 | 0.457 |
| gorilla KB3781 NPIP | oracle C | 8 / 25 | 10 (8) | 0.320 | 0.800 | 0.457 |
| gorilla KB3781 TBC1D3 | B0 pre-flip | 0 / 14 | 0 (0) | 0.000 | n/a | n/a |
| gorilla KB3781 TBC1D3 | Gc F1v2 only | 0 / 14 | 0 (0) | 0.000 | n/a | n/a |
| gorilla KB3781 TBC1D3 | Bc containment only | 0 / 14 | 0 (0) | 0.000 | n/a | n/a |
| gorilla KB3781 TBC1D3 | **D current** | 0 / 14 | 0 (0) | 0.000 | n/a | n/a |
| gorilla KB3781 TBC1D3 | oracle A | 0 / 14 | 0 (0) | 0.000 | n/a | n/a |
| gorilla KB3781 TBC1D3 | oracle C | 0 / 14 | 0 (0) | 0.000 | n/a | n/a |

Source: [cell]. (†) OR6737's B0 row is the stored genome-wide BASE restricted to the five contigs, not a five-contig re-run. "n/a" = the family cluster does not exist (no copy has a locus). On chr16 the six non-copy units of MCL1 under D are 1-2 transcript loci of 2-10 reads; by exon overlap with the copies'
territories, three overlap NPIPB5, NPIPB12 or NPIPB13 on the copy's own strand (72, 443 and 443 bp), two overlap NPIPB5 antisense and one overlaps no copy, so unit precision .800 (24 of 30) understates NPIP membership (27 of 30 by same-strand overlap, 29 of 30 by overlap on either strand).

### 3.2 Human chr16 (development block), NPIP; truths Compara, U2, Soto, Liftoff

| arm | Compara bipartite sens / prec / F | Compara pairs TP / predicted (truth 193) | Compara paralogue table (1,324 same-chromosome pairs): pairs recovered / pairs with both genes clustered (recall of all 1,324) | Compara CF153 (NPIP) sens / prec / F | U2 pooled F | U2 ID_154 sens / prec / F | Soto pooled F (not independent) | Liftoff recall (36 pairs) |
|---|---|---|---|---|---|---|---|---|
| B0 pre-flip | 0.444 / 1.000 / 0.615 | 66 / 66 | 66 / 89 (5.0%) | 0.526 / 1.000 / 0.690 | 0.576 | 0.632 / 0.706 / 0.667 | 0.617 | 15 / 36 |
| Bc containment only | 0.444 / 1.000 / 0.615 | 66 / 66 | 66 / 89 (5.0%) | 0.526 / 1.000 / 0.690 | 0.576 | 0.632 / 0.706 / 0.667 | 0.617 | 15 / 36 |
| Gc F1v2 only | 0.500 / 1.000 / 0.667 | 99 / 99 | 99 / 131 (7.5%) | 0.684 / 1.000 / 0.812 | 0.645 | 0.789 / 0.750 / 0.769 | 0.655 | 17 / 36 |
| **D current default** | 0.500 / 1.000 / 0.667 | 99 / 99 | 99 / 131 (7.5%) | 0.684 / 1.000 / 0.812 | 0.645 | 0.789 / 0.750 / 0.769 | 0.655 | 17 / 36 |
| oracle A (cov .70) | 0.537 / 1.000 / 0.699 | 109 / 109 | 109 / 132 (8.2%) | 0.737 / 1.000 / 0.848 | 0.635 | 0.789 / 0.714 / 0.750 | 0.667 | 15 / 36 |
| oracle C (cov .70) | 0.537 / 1.000 / 0.699 | 109 / 109 | 109 / 133 (8.2%) | 0.737 / 1.000 / 0.848 | 0.656 | 0.789 / 0.714 / 0.750 | 0.679 | 17 / 36 |
| oracle S (r1017 set, A) | 0.556 / 1.000 / 0.714 | 122 / 122 | 122 / 132 (9.2%) | 0.789 / 1.000 / 0.882 | 0.667 | 0.842 / 0.762 / 0.800 | 0.667 | 16 / 36 |
| oracle A (cov 0) | 0.537 / 1.000 / 0.699 | 111 / 111 | n/a | 0.737 / 1.000 / 0.848 | 0.635 | 0.789 / 0.714 / 0.750 | 0.667 | 15 / 36 |
| oracle C (cov 0) | 0.537 / 1.000 / 0.699 | 111 / 111 | n/a | 0.737 / 1.000 / 0.848 | 0.656 | 0.789 / 0.714 / 0.750 | 0.679 | 17 / 36 |
| oracle S (cov 0) | 0.556 / 1.000 / 0.714 | 122 / 122 | n/a | 0.789 / 1.000 / 0.882 | 0.667 | 0.842 / 0.762 / 0.800 | 0.667 | 16 / 36 |

Source: [sc16]; [pairs16]. Decomposition: F1v2 alone (Gc) gives the whole chr16 gain (Compara F .615 → .667, true pairs 66 → 99, CF153 F .690 → .812); containment alone (Bc)
gives none; D = Gc. The U2 NPIP family ID_154 goes .667 → .769. The paralogue-table column gives both the conditional count (pairs with both genes clustered, a denominator set by the prediction) and the recall of all 1,324 pairs. §3.1's .923 and this section's CF153 sensitivity
.684 are the same family and arm measured by two instruments (26 Dishuck copies by holder, 19 Compara genes by name); §3.7 reconciles them.

### 3.3 Human chr17, TBC1D3

| arm | Compara bipartite sens / prec / F | Compara pairs TP / predicted | Compara CF185 (TBC1D3, 9 genes) sens / prec / F | Soto ID_468 (TBC1D3, 10 genes; not independent) sens / prec / F | Liftoff recall (14 pairs) |
|---|---|---|---|---|---|
| B0 pre-flip | 0.267 / 0.909 / 0.412 | 41 / 44 | 1.000 / 1.000 / 1.000 | 0.800 / 1.000 / 0.889 | 1 / 14 |
| Bc containment only | 0.280 / 0.913 / 0.429 | 41 / 44 | 1.000 / 1.000 / 1.000 | 0.800 / 1.000 / 0.889 | 1 / 14 |
| Gc F1v2 only | 0.307 / 0.920 / 0.460 | 44 / 47 | 1.000 / 1.000 / 1.000 | 0.800 / 1.000 / 0.889 | 1 / 14 |
| **D current default** | 0.320 / 0.923 / 0.475 | 44 / 47 | 1.000 / 1.000 / 1.000 | 0.800 / 1.000 / 0.889 | 1 / 14 |
| oracle A (cov .70) | 0.320 / 0.923 / 0.475 | 44 / 47 | 1.000 / 1.000 / 1.000 | 0.900 / 0.900 / 0.900 | 1 / 14 |
| oracle C (cov .70) | 0.320 / 0.923 / 0.475 | 44 / 47 | 1.000 / 1.000 / 1.000 | 0.900 / 0.900 / 0.900 | 1 / 14 |
| oracle S (r1017 set, A) | 0.320 / 0.923 / 0.475 | 42 / 45 | 1.000 / 1.000 / 1.000 | 0.900 / 0.900 / 0.900 | 1 / 14 |
| oracle A (cov 0) | 0.307 / 0.920 / 0.460 | 44 / 47 | 1.000 / 1.000 / 1.000 | 0.900 / 0.900 / 0.900 | 1 / 14 |
| oracle C (cov 0) | 0.307 / 0.920 / 0.460 | 44 / 47 | 1.000 / 1.000 / 1.000 | 0.900 / 0.900 / 0.900 | 1 / 14 |
| oracle S (cov 0) | 0.293 / 0.917 / 0.444 | 42 / 45 | 1.000 / 1.000 / 1.000 | 0.900 / 0.900 / 0.900 | 1 / 14 |

Source: [sc17]. The Compara TBC1D3 family (CF185, 9 genes) is recovered perfectly in every arm (36 of 36 pairs); the chr17 pooled gain comes from other families
(F1v2 +.048, containment +.015).

### 3.4 The 19 A119b chromosomes already run (per-chromosome families; cross-chromosome edges absent, as in o1_cover)

| arm (19 A119b chromosomes, per-chromosome families) | Compara bipartite sens / prec / F | pairs TP / predicted (truth 2,746) | Liftoff recall (322 pairs) |
|---|---|---|---|
| B0 pre-flip | 0.291 / 0.936 / 0.445 | 551 / 643 | 24 / 322 |
| Bc containment only | 0.302 / 0.921 / 0.455 | 558 / 680 | 23 / 322 |
| Gc F1v2 only | 0.297 / 0.938 / 0.451 | 554 / 646 | 24 / 322 |
| **D current default** | 0.305 / 0.924 / 0.459 | 561 / 681 | 23 / 322 |
| D with oracle A on chr1 and chr19 only | 0.316 / 0.936 / 0.472 | 585 / 648 | n/a |
| D with oracle A on chr1, oracle C on chr19 | 0.314 / 0.943 / 0.471 | 583 / 635 | n/a |

| chromosome | arm | Compara sens / prec / F | pairs TP / predicted |
|---|---|---|---|
| chr1 | B0 pre-flip | 0.246 / 1.000 / 0.395 | 49 / 52 |
| chr1 | Bc containment only | 0.262 / 1.000 / 0.415 | 50 / 53 |
| chr1 | Gc F1v2 only | 0.246 / 1.000 / 0.395 | 49 / 52 |
| chr1 | **D current default** | 0.262 / 1.000 / 0.415 | 50 / 53 |
| chr1 | oracle A (cov .70) | 0.300 / 1.000 / 0.462 | 77 / 80 |
| chr1 | oracle A (cov 0) | 0.285 / 1.000 / 0.443 | 76 / 79 |
| chr19 | B0 pre-flip | 0.170 / 0.737 / 0.276 | 23 / 60 |
| chr19 | Bc containment only | 0.194 / 0.727 / 0.306 | 28 / 79 |
| chr19 | Gc F1v2 only | 0.176 / 0.744 / 0.284 | 23 / 60 |
| chr19 | **D current default** | 0.188 / 0.721 / 0.298 | 28 / 79 |
| chr19 | oracle A (cov .70) | 0.224 / 0.822 / 0.352 | 25 / 45 |
| chr19 | oracle C (cov .70) | 0.212 / 0.854 / 0.340 | 23 / 38 |
| chr19 | oracle A (cov 0) | 0.200 / 0.868 / 0.325 | 22 / 33 |
| chr19 | oracle C (cov 0) | 0.188 / 0.838 / 0.307 | 20 / 31 |

Source: [sca19]. Scored minus chr16, chr18, chr13, chr20, chr21, chr22 (the o1_cover substrate); 372 truth families, 1,163 genes. The B0 and Gc rows equal the o1_cover BASE and
COVER-core rows (.2915 / .9365 / .4446 and .2966 / .9375 / .4507 in the rounding of that report). The current defaults add +.014 F over the pre-flip pipeline (F1v2 +.006,
containment +.008) and 10 true pairs; Liftoff recall is 24 → 23 of 322.

**Pooled truth versus per-chromosome truth.** The 19-chromosome rows score the concatenated per-chromosome products against the pooled Compara truth: 372 families with >= 2 genes, 1,163 genes, 2,746 pairs, of which 2,394 lie on one chromosome (the per-chromosome
universe of §4.3: 1,062 genes) and 352 join genes on different chromosomes (63 families span more than one chromosome). Per-chromosome families can contain none of those 352, so pooled pair sensitivity is capped at 2,394 / 2,746 = 87.2%; predicted pairs are taken over the pooled
1,163-gene universe and are therefore not the sum of the per-chromosome rows (chr1 and chr19 rows of §4.2 against the pooled oracle row). Levels are not comparable between the pooled rows and the per-chromosome rows; differences between arms inside one table are.

### 3.5 Human testis

| arm | Compara genome-wide minus chr16/chr18: sens / prec / F | pairs TP / predicted | Liftoff recall genome-wide (33 pairs) | chr16 block: Compara F (pairs) | chr16 U2 F | chr16 Liftoff (5 pairs) |
|---|---|---|---|---|---|---|
| B0 pre-flip | 0.153 / 0.950 / 0.264 | 210 / 227 | 10 / 33 | 0.500 (12/12) | 0.256 | 1 / 5 |
| Bc containment only | 0.153 / 0.950 / 0.264 | 210 / 227 | 10 / 33 | 0.500 (12/12) | 0.256 | 1 / 5 |
| Gc F1v2 only | 0.154 / 0.950 / 0.265 | 212 / 229 | 11 / 33 | 0.500 (12/12) | 0.256 | 1 / 5 |
| **D current default** | 0.154 / 0.950 / 0.265 | 212 / 229 | 11 / 33 | 0.500 (12/12) | 0.256 | 1 / 5 |

Source: [sctes]. Testis is shallow (7 of 26 NPIP copies have a locus, 5 in the family cluster); the flip moves genome-wide Compara F .264 → .265 and Liftoff recall 10 → 11 of 33.

### 3.6 Gorilla

| sample | arm | Liftoff pairs recovered / pairs (five contigs) |
|---|---|---|
| gorilla_OR6737 | B0 pre-flip | 1 / 27 |
| gorilla_OR6737 | Gc F1v2 only | 1 / 27 |
| gorilla_OR6737 | **D current default** | 1 / 27 |
| gorilla_OR6737 | oracle A (cov .70) | 1 / 27 |
| gorilla_OR6737 | oracle A (cov 0) | 1 / 27 |
| gorilla_KB3781 | B0 pre-flip | 1 / 32 |
| gorilla_KB3781 | Bc containment only | 1 / 32 |
| gorilla_KB3781 | Gc F1v2 only | 1 / 32 |
| gorilla_KB3781 | **D current default** | 1 / 32 |
| gorilla_KB3781 | oracle A (cov .70) | 1 / 32 |
| gorilla_KB3781 | oracle C (cov .70) | 1 / 32 |
| gorilla_KB3781 | oracle A (cov 0) | 1 / 32 |
| gorilla_KB3781 | oracle C (cov 0) | 1 / 32 |

Source: [scgor]. Liftoff recovers 1 of 27 (OR6737) and 1 of 32 (KB3781) pairs on the five contigs in every arm: the gorilla family references have no power at
these substrates, which is why §3.1 scores the gorilla families per copy.

### 3.7 Two instruments on human chr16 NPIP: why §3.1 and §3.2 differ

§3.1 counts Dishuck copies whose holder locus is in the NPIP cluster; §3.2 counts truth genes of Compara CF153 (19 genes), U2 ID_154 (19) and Soto ID_154 (14) that the family scorer places in the matched cluster by gene name. The truth sets differ and so does the resolution of a
copy to a locus.

| arm | holder rule: Dishuck copies whose holder is in the NPIP cluster (of 26) | scorer: Compara CF153 genes hit (of 19) | scorer: U2 ID_154 genes hit (of 19) | scorer: Soto ID_154 genes hit (of 14) | scorer: genes collapsed into another gene's locus (Compara / U2 / Soto) |
|---|---|---|---|---|---|
| B0 pre-flip | 22 | 10 | 12 | 10 | 5 / 9 / 8 |
| Gc F1v2 only | 24 | 13 | 15 | 13 | 3 / 7 / 6 |
| **D current default** | 24 | 13 | 15 | 13 | 3 / 7 / 6 |
| oracle A | 25 | 14 | 15 | 13 | 1 / 6 / 5 |
| oracle C | 25 | 14 | 15 | 13 | 1 / 6 / 5 |
| oracle S | 25 | 15 | 16 | 13 | 0 / 7 / 5 |

Per-copy disagreements (D; oracle gains against D) [fate, sc16]:

| truth family (scorer) | genes hit under D that the holder rule calls wrong / unclustered | genes missed under D that the holder rule calls right | oracle A gain: scorer hits gained | oracle A gain: holder-rule copies gained | oracle S gain: scorer hits gained |
|---|---|---|---|---|---|
| Compara CF153 | NPIPB5 (wrong family) | NPIPA1, NPIPA6, NPIPA9, NPIPB13, NPIPB3, NPIPB4 | NPIPB3 | NPIPB5 | NPIPA9, NPIPB3 |
| U2 ID_154 | none | LOC124907834, NPIPA6, NPIPA9, NPIPB14P | none | NPIPB5 | NPIPA9 |
| Soto ID_154 | none | NPIPB14P | none | NPIPB5 | none |

Reading. (i) The instruments agree on the size of each step: B0 → D is +2 copies by the holder rule and +3 hits under each of CF153, U2 and Soto; D → oracle A is +1 copy by the holder rule, +1 hit under CF153, 0 under U2 and Soto. (ii) They disagree on WHICH copy: the holder rule
recovers NPIPB5 under the oracle, while the scorer already counts NPIPB5 as a CF153 hit under D and gains NPIPB3 (A, C) or NPIPB3 and NPIPA9 (S). (iii) Under D the scorer misses six CF153 genes the holder rule counts as placed (NPIPA1, A6, A9, B13, B3, B4) and four U2 genes
(LOC124907834, NPIPA6, NPIPA9, NPIPB14P): the 09-21 one-name-per-locus trap, in which a locus that spans two genes is named for one of them; which locus each missed name resolves to was not examined. (iv) Every "copies recovered" and "lost to fused loci" number in §0-§2, §4.1, §4.3 and §6-§7 is a holder-rule
number, and every F in §3.2-§3.5 and §4.2 is a scorer number; a statement about fused-locus losses has to carry its instrument, and the holder rule's 41 of 1,062 (§4.3) is not the scorer's count.

## 4. Headroom (deliverable 3)

### 4.1 The oracle, rebuilt for the current pipeline

r1017 (`docs/PREREG_overmerge_ceiling_2026-09-22.md`) replaced every fused chr16 locus by one node per constituent annotated gene on a pre-flip, pre-containment baseline and
read +.021 referee F and a worse NPIP. Here the nodes are gene_id groups of the families-input GTF, so the shipped families stage runs on the oracle GTF unchanged
(§1, same aligner call). A fused locus is the locus-union reading above (readthrough_eval gene set); constituents are the records its exon union overlaps on its strand.
Variants A (clip to the annotated span) and C (exons whole to the record they overlap most) are r1017's `oracle_a` and `oracle_c`; S is r1017's own selection.

| substrate | loci | fused loci (locus union) | nodes after A / C | r1017 selection (S) loci → nodes |
|---|---|---|---|---|
| human A119b chr16 | 2,832 | 83 | 189 / 175 | 240 → 351 (r1017 counted 233 on its pre-flip baseline) |
| human A119b chr17 | 2,880 | 97 | 220 / 199 | 319 → 414 |
| human A119b chr19 | 2,415 | 118 | 259 / 236 | not run |
| human A119b chr1 | 7,153 | 171 | 383 (A) | not run |
| human testis chr16 | 586 | 11 | 23 / 20 | not run |
| gorilla OR6737 (five contigs) | 5,601 | 183 | 410 (A) | not run |
| gorilla KB3781 (five contigs) | 5,222 | 167 | 357 / 318 | not run |

What the oracle does, checked independently on chr16 OC [ver]: the 83 fused loci become 175 pieces (2,832 → 2,924 loci), the other 2,749 loci are byte-identical, and every original transcript's exons are conserved across
its pieces (1,125 transcripts, 10,687 exons, none lost or duplicated). "Overlaps most" is overlap with the record's span. Six of the 83 loci stay one piece because the second constituent wins no exon (NDE1, two
VKORC1 loci, KAT8, SETD6, LOC124903710) and eight more lose one constituent (NPIPB5's locus loses LOC124907830); 672 of the 10,687 exons (6.3%) overlap no record and go to the nearest constituent. So C is
"one node per constituent that wins an exon", A "one node per constituent with an exon inside its span".

**Coverage of the oracle.** Of the 19 A119b chromosomes, chr1, chr17 and chr19 were split (the pooled 19-chromosome oracle row splits chr1 and chr19 only, so 17 of the 19 are unsplit in that row and 16 were never split); outside the 19, chr16 (the development block), testis chr16 and both gorillas were also split. chr1 and chr19 were chosen because they hold 24 of the 41 fused-and-unplaced
Compara genes and 27 of the 34 genes in another cluster, a choice conditioned on the outcome, and the oracle can lower scores (chr19 true pairs 28 → 25 / 23; U2 ID_154 -.019 under A and C; Liftoff 17 → 15 on chr16 under A). The sign of the effect on the unsplit chromosomes is therefore unknown, and the genome-wide gain is an
estimate from two chromosomes, neither a lower nor an upper bound. Gorilla OR6737 was run with A only (the one copy with a fused locus outside its family, NPIPB7, is a touch-only holder).

### 4.2 Gain per family and truth (oracle minus D)

| measure | D current default | oracle A | oracle C | oracle S (r1017 selection, A) |
|---|---|---|---|---|
| chr16 NPIP: Dishuck copies in the family cluster (of 26) | 24 | 25 | 25 | 25 |
| chr16 NPIP: unit precision (units holding a copy / units) | 0.800 (24/30) | 0.862 (25/29) | 0.862 (25/29) | 0.833 (25/30) |
| chr16 Compara pooled F (true pairs) | 0.667 (99) | 0.699 (109) | 0.699 (109) | 0.714 (122) |
| chr16 Compara CF153 (NPIP) F | 0.812 | 0.848 | 0.848 | 0.882 |
| chr16 U2 pooled F | 0.645 | 0.635 | 0.656 | 0.667 |
| chr16 U2 ID_154 (NPIP main) F | 0.769 | 0.750 | 0.750 | 0.800 |
| chr16 Soto pooled F (not independent) | 0.655 | 0.667 | 0.679 | 0.667 |
| chr16 Soto ID_154 (NPIP, 14 genes) F (not independent; r1017's NPIP-family metric) | 0.867 | 0.929 | 0.929 | 0.929 |
| chr16 Compara: genes collapsed into another gene's locus (scorer) | 3 | 1 | 1 | 0 |
| chr16 Soto: truth genes with no locus (scorer) | 17 | 18 | 17 | 19 |
| chr16 Liftoff pairs recovered (of 36) | 17 | 15 | 17 | 16 |
| chr17 TBC1D3: copies in the family cluster (of 16) | 11 | 12 | 12 | 12 |
| chr17 Compara pooled F (true pairs) | 0.475 (44) | 0.475 (44) | 0.475 (44) | 0.475 (42) |
| chr17 Compara CF185 (TBC1D3) F | 1.000 | 1.000 | 1.000 | 1.000 |
| chr17 Soto pooled F (not independent) | 0.487 | 0.557 | 0.531 | 0.559 |
| chr17 Soto ID_468 (TBC1D3) F (not independent) | 0.889 | 0.900 | 0.900 | 0.900 |
| chr1 Compara pooled F (true pairs) | 0.415 (50) | 0.462 (77) | not run | not run |
| chr19 Compara pooled F (true pairs) | 0.298 (28) | 0.352 (25) | 0.340 (23) | not run |
| 19 A119b chromosomes pooled Compara F, oracle on chr1 + chr19 only | 0.459 | 0.472 | 0.471 (A on chr1, C on chr19) | not run |
| gorilla OR NPIP: copies in the family cluster (of 25) | 10 | 10 | not run | not run |
| gorilla OR TBC1D3 (of 14) | 7 | 7 | not run | not run |
| gorilla KB NPIP (of 25) | 8 | 8 | 8 | not run |
| gorilla KB TBC1D3 (of 14) | 0 | 0 | 0 | not run |
| human testis NPIP (of 26) | 5 | 5 | 5 | not run |

Source: [sc16, sc17, sc1, sc19, sca19o, cell]. Reading. (i) Copies (holder rule): +1 on chr16 and +1 on chr17, the two fusion-related copies of §2.2, under all three oracles; nothing elsewhere.
(ii) Compara F (scorer) rises where fused loci hold family members (chr16 +.032 to +.047, chr1 +.047, chr19 +.042 to +.054) and is flat on chr17 (.475 = .475; Soto pooled F on chr17 does rise, .487 → .557 / .531 / .559, because the oracle changes Soto's node universe); the Compara gains are mostly
fewer false pairs on chr19 (predicted pairs 79 → 45 / 38 at a cost of 3 to 5 true pairs) and more true pairs on chr1 (50 → 77, no new false pair). Part of this is the scorer's own collision count (collapsed genes 3 → 1 → 0 on chr16, row above): names made exact, not families found.
(iii) r1017's "NPIP worse" was Soto's NPIP-family F (.727 → .687) on a pre-flip, pre-containment baseline. Here Soto ID_154 rises .867 → .929 (A, C, S) and Compara CF153 .812 → .848 / .882, so it does not reproduce under those two truths; U2 ID_154 falls .769 → .750 under A and C (one extra false unit, precision .750 → .714) and rises to .800
under S, which is r1017's direction, small. On the copy metric every oracle adds a copy. r1018's mechanism ("pieces fail the coverage gate") was not re-examined here, so "does not reproduce" says nothing about why. (iv) Liftoff recall does not rise (chr16 17 → 15 under A, 17 under C, 16 under S, of 36; chr17 1 of 14; gorilla 1 of 27 and 1 of 32).
(v) Pooled truth over the 19 A119b chromosomes with only chr1 and chr19 split: F .459 → .472, sensitivity .305 → .316, precision .924 → .936 (§3.4 on the pooled truth).

### 4.3 Where fused loci still cost members genome-wide (Compara multi-copy genes, per chromosome, current defaults D; holder rule)

| chromosome | fused loci (locus-union / per-transcript) of all loci | Compara multi-copy genes scored | in matched cluster | no locus | unclustered | other cluster | present but not placed: fused locus-union (per-transcript) |
|---|---|---|---|---|---|---|---|
| chr16 (dev) | 83 / 82 of 2832 | 54 | 36 | 9 | 7 | 2 | 3 (3) |
| chr1 | 171 / 167 of 7153 | 130 | 28 | 61 | 28 | 13 | 4 (4) |
| chr2 | 113 / 113 of 7215 | 48 | 23 | 13 | 12 | 0 | 1 (1) |
| chr3 | 96 / 93 of 6060 | 12 | 3 | 8 | 1 | 0 | 0 (0) |
| chr4 | 65 / 64 of 5181 | 42 | 4 | 11 | 27 | 0 | 2 (2) |
| chr5 | 90 / 89 of 5291 | 31 | 7 | 14 | 9 | 1 | 2 (2) |
| chr6 | 104 / 103 of 5042 | 35 | 10 | 13 | 12 | 0 | 2 (2) |
| chr7 | 123 / 122 of 5173 | 60 | 43 | 8 | 9 | 0 | 1 (1) |
| chr8 | 54 / 54 of 4110 | 53 | 5 | 43 | 5 | 0 | 0 (0) |
| chr9 | 79 / 76 of 3830 | 42 | 14 | 22 | 5 | 1 | 1 (1) |
| chr10 | 74 / 71 of 3906 | 33 | 16 | 7 | 10 | 0 | 2 (2) |
| chr11 | 90 / 89 of 3844 | 77 | 13 | 51 | 13 | 0 | 1 (1) |
| chr12 | 88 / 88 of 4114 | 27 | 3 | 15 | 9 | 0 | 2 (2) |
| chr14 | 74 / 72 of 3140 | 14 | 7 | 6 | 1 | 0 | 0 (0) |
| chr15 | 81 / 80 of 3438 | 35 | 28 | 4 | 0 | 3 | 0 (0) |
| chr17 | 97 / 94 of 2880 | 75 | 28 | 41 | 6 | 0 | 3 (3) |
| chr19 | 118 / 117 of 2415 | 165 | 35 | 66 | 50 | 14 | 20 (20) |
| chrX | 20 / 20 of 1819 | 153 | 106 | 27 | 18 | 2 | 0 (0) |
| chrY | 13 / 13 of 664 | 30 | 27 | 2 | 1 | 0 | 0 (0) |
| chrM | 0 / 0 of 17 | 0 | 0 | 0 | 0 | 0 | 0 (0) |
| **19 A119b chromosomes (sum)** | **1550 / 1525 of 75292** | **1062** | **400** | **412** | **216** | **34** | **41 (41)** |
| human testis, genome-wide | 184 / 183 of 13034 | 1313 | 247 | 885 | 166 | 15 | 1 (0) |

Source: [gf, fused]. A gene is scored when its Compara family has >= 2 genes on the chromosome (family_score's rule; testis, run genome-wide: >= 2 genes anywhere in the annotation, 1,313 genes, a third rule); "in matched cluster" = its holder is in the cluster that holds most of its
family's genes; fused columns count the present-but-unplaced genes whose holder is fused. Among the 250 present-but-unplaced genes of the 19 chromosomes, 41 sit in a fused locus
(35 unclustered, 6 in another cluster); 10 genes have a partner's transcript as the locus representative (8 of the 250, all unclustered, and 2 placed right) and all 10 are in fused loci; 40 of the 400 correctly placed genes also sit in fused loci,
so a fused locus does not by itself prevent placement. The 41 are 3.9% of the 1,062 scored genes and 16.4% of the 250 present-but-unplaced (a share conditioned on being unplaced). The 1,550 fused loci of the 19 chromosomes are 2.1% of their 75,292 loci.

The fused loci that remain, by the class of their strongest link (locus-union reading; [fl]): 19 chromosomes 1,550 = 664 dominant (link reads >= both sides) + 678 mixed + 183 minority (link reads < both sides; F1v2 would
remove them if the END proofs held) + 25 with no link transcript; chr16 83 = 40 / 36 / 6 / 1; testis 184 = 102 / 71 / 10 / 1. F1v2 cut 336 bridge loci on the 19 chromosomes (16 of them on chr17) and, separately, 12 on chr16 and 14 in the genome-wide testis run (`bridges.tsv` rows with `keep` True).
So the links F1v2 leaves in loci are mostly dominant or mixed, which is the premise of the record's container question: they cannot be cut by a read-share rule.

### 4.4 Real loci that join two multi-copy families, and how they are linked

| substrate | truth | loci joining two multi-copy families | dominant | mixed | minority (in a locus: the END proofs are absent or fail, inferred from presence in `families.gtf`) | no link transcript | F1v2 bridge relations joining two | examples |
|---|---|---|---|---|---|---|---|---|
| A119b chr16 (dev) | Compara | 3 | 3 | 0 | 0 | 0 | 0 | EIF3C\|NPIPB9, SLX1B\|SULT1A4, SLX1A\|SULT1A3 |
| A119b chr16 (dev) | U2 (NPIP union) | 3 | 2 | 1 | 0 | 0 | 1 | PKD1P1\|NPIPA6, PKD1P2\|NPIPA7, NPIPB14P\|PDXDC2P |
| A119b chr16 (dev) | Liftoff multi-copy | 0 | 0 | 0 | 0 | 0 | 0 |  |
| A119b 19 chromosomes | Compara | 15 | 5 | 5 | 5 | 0 | 0 | NOTCH2NLB\|NBPF14, NOTCH2NLC\|NBPF19, NOTCH2NLR\|NBPF26, UGT1A8\|UGT1A5, UGT1A8\|UGT1A5 ... |
| A119b 19 chromosomes | Liftoff multi-copy | 7 | 4 | 2 | 0 | 1 | 0 | SEPTIN14P2\|RPL23AP88, SEPTIN14P3\|LOC124906332, SEPTIN14P3\|LOC124906332, LOC124906410\|USP17L15, LOC124906515\|LOC124906510 ... |
| human testis (genome-wide) | Compara | 3 | 2 | 1 | 0 | 0 | 0 | NOTCH2NLR\|NBPF26, PRH1\|TAS2R19, EIF3C\|NPIPB9 |
| human testis (genome-wide) | U2 (NPIP union) | 0 | 0 | 0 | 0 | 0 | 0 |  |
| human testis (genome-wide) | Liftoff multi-copy | 1 | 0 | 1 | 0 | 0 | 0 | LOC284412\|LOC124908084 |
| gorilla_OR6737 (genome-wide) | Liftoff multi-copy | 1 | 0 | 0 | 0 | 1 | 0 | LOC115931102\|LOC129527589 |
| gorilla_KB3781 (genome-wide) | Liftoff multi-copy | 0 | 0 | 0 | 0 | 0 | 0 |  |

Source: [jn]. A locus joins two multi-copy families under a truth when its genes (exon-overlap, locus strand, readthrough records excluded) include two genes of two different multi-copy
families of that truth (Compara: >= 2 genes present; U2: its NPIP and PKD1 groups; Liftoff: a record with >= 1 extra copy at sequence_ID >= 0.95 and both exon unions >= 200 bp, genes
joined by those copies form a family). Soto is not in this table: it is labelled not independent and its families are a cover, which makes "joining two families" ambiguous. The class is that of the pair with the most link reads in the locus; "no link transcript" = the genes are joined through overlapping
transcripts of different isoforms only. "Minority" loci are in `families.gtf`, which is the only evidence that F1v2 did not cut them (the joins product holds no proof field). F1v2 relation records are the `fusion_of` transcripts of the assembly GTF tested the same way (the U2 record's join uses its first transcript, 2 reads; the relation has two transcripts, 6 reads). No F1v2 bridge joins two multi-copy families
on any substrate here except the one U2 record; the joins that exist are loci (dominant, mixed, or minority). On chr16 the Compara joins are EIF3C | NPIPB9
(697 reads, dominant), SLX1B | SULT1A4 and SLX1A | SULT1A3; the U2 joins are the PKD1P1 | NPIPA6, PKD1P2 | NPIPA7 and NPIPB14P | PDXDC2P loci. The count is specific to the exon-union, same-strand definition [ver]: the
15 loci on the 19 chromosomes are 13 distinct gene-set joins (UGT1A on chr2 and ZNF587 / ZNF587B on chr19 occur in two loci each); with the records' spans instead of their exon unions the same loci give 38
(chr16 3, testis 5), and ignoring strand 45 (chr16 5, testis 12); those two sensitivity readings were run on loci, not on the F1v2 relation records, for which the support is register 1185 (0 of 348 A119b bridges and 0 of 12 testis bridges have parents in two Compara families).

### 4.5 What a partition of units plus a derived relation could express (oracle)

| cell / chromosome | oracle | fused loci split | loci whose pieces land in >= 2 families (a relation) | of these, involving the NPIP / TBC1D3 family | distinct family pairs |
|---|---|---|---|---|---|
| human chr16 NPIP | OA | 83 | 12 | 5 | 22 |
| human chr16 NPIP | OC | 83 | 14 | 6 | 24 |
| human chr17 TBC1D3 | OA | 97 | 9 | 3 | 5 |
| human chr17 TBC1D3 | OC | 97 | 5 | 3 | 4 |
| gorilla KB3781 NPIP | OA | 167 | 14 | 0 | 9 |
| gorilla OR6737 NPIP | OA | 183 | 19 | 0 | 16 |
| human chr19 (A119b) | OA | 118 | 5 | n/a (no family cluster given) | 5 |
| human chr19 (A119b) | OC | 118 | 4 | n/a (no family cluster given) | 4 |
| human chr1 (A119b) | OA | 171 | 7 | n/a (no family cluster given) | 5 |

Source: [rel]. A fused locus "has a relation" when the nodes of its constituent genes fall in two or more clusters. On human chr16, 12 to 14 of 83 fused loci do (22 to 24 distinct
family pairs) and 5 to 6 involve the NPIP family; on chr17, 5 to 9 of 97 and 3 involve TBC1D3. In gorilla 14 to 19 of 167 to 183 do and none involves the NPIP family. These are counts of oracle cluster assignments; no truth scores a relation, so they are not a validated gain.

## 5. Fusion simulation under the current defaults (deliverable 4)

Gorilla simulation of `PREREG_fusion_container_sim_2026-09-28.md`: 25 T_member copies, 20 with a fusable neighbour, 10 fused with their partner, fusion share f of the copy's 30 reads
(f = 0, .1, .5, .9, 1.0), 30 reads per partner transcript. Scorer: the frozen `score.py` (8f512fb8) with paths and the families-input GTF changed only (`fsim_score.py`, e2b1f07b;
diff = path constants, `families.gtf` for the f1v2 arms, the Rust container file). M1 = fused copies whose holder is in the NPIP cluster (/10); M2 = S partners whose holder is in the NPIP
cluster (/20); M3 = container precision / recall against the partner's bases; M4 = whether a fused member's container relates to its partner's family; M5 = accessory bp on unfused holders.
Driver: `tools/rustle_pipeline.sh assemble` then `families` with `RUSTLE_FAMILY_CONTAINER=1`, `--bin` = the copied binaries.

### 5.1 Extended pool (`locus_fix_design` arms: realistic read ends; this is the pool of F1v2's dev numbers 9 / 5 / 0 at f = .1 / .5 / .9 and BASE 7 / 1 / 0)

| f | pipeline | NPIP cluster units | copies in NPIP /25 | **M1** fused copies in NPIP /10 | unfused controls in NPIP /15 | **M2** S partners in NPIP /20 (fused partners /10) | **M3** prec / rec | M3 members | **M4** (no-reason counts) | **M5** unfused accessory bp (units) | M5 fused@f0 bp |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.0 | pre-flip (off, cov 0) | 26 | 23 | 9 | 14 | 0 (0) | 0.0 / 0.0 | 10 | partner elsewhere, unrelated 9; partner unclustered 1 | 11472 (4) | 132 |
| 0.0 | F1v2 only (cov 0) | 26 | 23 | 9 | 14 | 0 (0) | 0.0 / 0.0 | 10 | partner elsewhere, unrelated 9; partner unclustered 1 | 11472 (4) | 132 |
| 0.0 | containment only (off, cov .70) | 23 | 22 | 9 | 13 | 0 (0) | 0.0 / 0.0 | 10 | partner elsewhere, unrelated 9; partner unclustered 1 | 11370 (3) | 8956 |
| 0.0 | **current default** | 23 | 22 | 9 | 13 | 0 (0) | 0.0 / 0.0 | 10 | partner elsewhere, unrelated 9; partner unclustered 1 | 11370 (3) | 8956 |
| 0.1 | pre-flip (off, cov 0) | 26 | 22 | 7 | 15 | 0 (0) | 0.013 / 0.0141 | 10 | partner elsewhere, unrelated 6; partner IS the fused member 3; partner unclustered 1 | 11370 (3) | - |
| 0.1 | F1v2 only (cov 0) | 26 | 23 | 9 | 14 | 0 (0) | 0.0 / 0.0 | 10 | partner elsewhere, unrelated 9; partner unclustered 1 | 11472 (4) | - |
| 0.1 | containment only (off, cov .70) | 22 | 21 | 7 | 14 | 0 (0) | 0.013 / 0.0141 | 10 | partner elsewhere, unrelated 6; partner IS the fused member 3; partner unclustered 1 | 20364 (4) | - |
| 0.1 | **current default** | 23 | 22 | 9 | 13 | 0 (0) | 0.0 / 0.0 | 10 | partner elsewhere, unrelated 9; partner unclustered 1 | 11370 (3) | - |
| 0.5 | pre-flip (off, cov 0) | 19 | 16 | 1 | 15 | 0 (0) | 0.195 / 0.4174 | 9 | partner elsewhere, unrelated 1; partner IS the fused member 8; member not clustered 1 | 11370 (3) | - |
| 0.5 | F1v2 only (cov 0) | 21 | 19 | 5 | 14 | 0 (0) | 0.4335 / 0.1152 | 9 | partner elsewhere, unrelated 6; partner IS the fused member 3; member not clustered 1 | 11472 (4) | - |
| 0.5 | containment only (off, cov .70) | 14 | 14 | 1 | 13 | 0 (0) | 0.195 / 0.4174 | 9 | partner elsewhere, unrelated 1; partner IS the fused member 8; member not clustered 1 | 11370 (3) | - |
| 0.5 | **current default** | 18 | 18 | 5 | 13 | 0 (0) | 0.1938 / 0.1152 | 9 | partner elsewhere, unrelated 6; partner IS the fused member 3; member not clustered 1 | 11370 (3) | - |
| 0.9 | pre-flip (off, cov 0) | 18 | 15 | 0 | 15 | 0 (0) | 0.3145 / 0.4174 | 9 | partner IS the fused member 9; member not clustered 1 | 11370 (3) | - |
| 0.9 | F1v2 only (cov 0) | 18 | 15 | 0 | 15 | 0 (0) | 0.3145 / 0.4174 | 9 | partner IS the fused member 9; member not clustered 1 | 11370 (3) | - |
| 0.9 | containment only (off, cov .70) | 12 | 12 | 0 | 12 | 0 (0) | 0.3145 / 0.4174 | 9 | partner IS the fused member 9; member not clustered 1 | 11646 (4) | - |
| 0.9 | **current default** | 12 | 12 | 0 | 12 | 0 (0) | 0.3145 / 0.4174 | 9 | partner IS the fused member 9; member not clustered 1 | 11646 (4) | - |
| 1.0 | pre-flip (off, cov 0) | 19 | 16 | 1 | 15 | 1 (1) | 0.628 / 0.4405 | 9 | partner IS the fused member 9; member not clustered 1 | 11370 (3) | - |
| 1.0 | F1v2 only (cov 0) | 19 | 16 | 1 | 15 | 1 (1) | 0.628 / 0.4405 | 9 | partner IS the fused member 9; member not clustered 1 | 11370 (3) | - |
| 1.0 | containment only (off, cov .70) | 15 | 14 | 1 | 13 | 1 (1) | 0.628 / 0.4405 | 9 | partner IS the fused member 9; member not clustered 1 | 11370 (3) | - |
| 1.0 | **current default** | 15 | 14 | 1 | 13 | 1 (1) | 0.628 / 0.4405 | 9 | partner IS the fused member 9; member not clustered 1 | 11370 (3) | - |

### 5.2 Container-sim pool (the pool of the record's container Outcome: error-free reads with 0-30 bp end jitter, which moves the 3' mode out of the PAS window of F1's proofs)

| f | pipeline | NPIP cluster units | copies in NPIP /25 | **M1** fused copies in NPIP /10 | unfused controls in NPIP /15 | **M2** S partners in NPIP /20 (fused partners /10) | **M3** prec / rec | M3 members | **M4** (no-reason counts) | **M5** unfused accessory bp (units) | M5 fused@f0 bp |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.0 | pre-flip (off, cov 0) | 27 | 23 | 9 | 14 | 0 (0) | 0.0 / 0.0 | 10 | partner elsewhere, unrelated 9; partner unclustered 1 | 10712 (4) | 132 |
| 0.0 | **current default** | 25 | 24 | 9 | 15 | 0 (0) | 0.0 / 0.0 | 10 | partner elsewhere, unrelated 9; partner unclustered 1 | 10610 (3) | 8950 |
| 0.1 | pre-flip (off, cov 0) | 20 | 15 | 1 | 14 | 1 (1) | 0.3173 / 0.5328 | 9 | partner IS the fused member 9; member not clustered 1 | 11908 (4) | - |
| 0.1 | **current default** | 15 | 15 | 2 | 13 | 0 (0) | 0.2439 / 0.3984 | 9 | partner IS the fused member 5; partner elsewhere, unrelated 4; member not clustered 1 | 10798 (5) | - |
| 0.5 | pre-flip (off, cov 0) | 20 | 14 | 0 | 14 | 0 (0) | 0.2136 / 0.4174 | 9 | partner IS the fused member 9; member not clustered 1 | 11908 (4) | - |
| 0.5 | **current default** | 13 | 13 | 0 | 13 | 0 (0) | 0.215 / 0.4174 | 9 | partner IS the fused member 8; partner elsewhere, unrelated 1; member not clustered 1 | 10798 (5) | - |
| 0.9 | pre-flip (off, cov 0) | 20 | 14 | 0 | 14 | 0 (0) | 0.2137 / 0.4174 | 9 | partner IS the fused member 9; member not clustered 1 | 11908 (4) | - |
| 0.9 | **current default** | 13 | 13 | 0 | 13 | 0 (0) | 0.2295 / 0.4174 | 9 | partner IS the fused member 8; partner elsewhere, unrelated 1; member not clustered 1 | 10798 (5) | - |
| 1.0 | pre-flip (off, cov 0) | 24 | 18 | 4 | 14 | 4 (4) | 0.3891 / 0.6212 | 10 | partner IS the fused member 10 | 11386 (4) | - |
| 1.0 | **current default** | 18 | 18 | 3 | 15 | 3 (3) | 0.4092 / 0.6212 | 10 | partner IS the fused member 10 | 10610 (3) | - |

### 5.3 The prereg bars (prereg §3), per pipeline

| pool | pipeline | M1 bar (at f <= .5 the count equals the f = 0 control) | M2 bar (at every f no more than the f = 0 control) | M3 bar (pooled precision and recall >= .90 at f >= .5) |
|---|---|---|---|---|
| extended | pre-flip | f = 0.1: FAIL (7 vs 9), f = 0.5: FAIL (1 vs 9) | f = 0.1: pass (0), 0.5: pass (0), 0.9: pass (0), 1.0: FAIL (1) | f = 0.5: FAIL (0.195 / 0.4174), 0.9: FAIL (0.3145 / 0.4174), 1.0: FAIL (0.628 / 0.4405) |
| extended | F1v2 only | f = 0.1: pass (9 vs 9), f = 0.5: FAIL (5 vs 9) | f = 0.1: pass (0), 0.5: pass (0), 0.9: pass (0), 1.0: FAIL (1) | f = 0.5: FAIL (0.4335 / 0.1152), 0.9: FAIL (0.3145 / 0.4174), 1.0: FAIL (0.628 / 0.4405) |
| extended | containment only | f = 0.1: FAIL (7 vs 9), f = 0.5: FAIL (1 vs 9) | f = 0.1: pass (0), 0.5: pass (0), 0.9: pass (0), 1.0: FAIL (1) | f = 0.5: FAIL (0.195 / 0.4174), 0.9: FAIL (0.3145 / 0.4174), 1.0: FAIL (0.628 / 0.4405) |
| extended | **current default** | f = 0.1: pass (9 vs 9), f = 0.5: FAIL (5 vs 9) | f = 0.1: pass (0), 0.5: pass (0), 0.9: pass (0), 1.0: FAIL (1) | f = 0.5: FAIL (0.1938 / 0.1152), 0.9: FAIL (0.3145 / 0.4174), 1.0: FAIL (0.628 / 0.4405) |
| container-sim | pre-flip | f = 0.1: FAIL (1 vs 9), f = 0.5: FAIL (0 vs 9) | f = 0.1: FAIL (1), 0.5: pass (0), 0.9: pass (0), 1.0: FAIL (4) | f = 0.5: FAIL (0.2136 / 0.4174), 0.9: FAIL (0.2137 / 0.4174), 1.0: FAIL (0.3891 / 0.6212) |
| container-sim | **current default** | f = 0.1: FAIL (2 vs 9), f = 0.5: FAIL (0 vs 9) | f = 0.1: pass (0), 0.5: pass (0), 0.9: pass (0), 1.0: FAIL (3) | f = 0.5: FAIL (0.215 / 0.4174), 0.9: FAIL (0.2295 / 0.4174), 1.0: FAIL (0.4092 / 0.6212) |

Source: [sim]. M1's ceiling is 9, not 10: NPIPA7 is already held by a locus outside the NPIP cluster at f = 0 in every current-default arm of both pools (the extended pools also have the unfused NPIPB3
and NPIPB6 outside at f = 0) [ver]. M3 is undefined (printed 0.0 / 0.0 by the scorer) where a pool has no accessory partner base. Reading.
- F1v2 reproduces its dev numbers on the extended pool (9 / 9 / 5 / 0 / 1 and, with containment off, copies 23 / 23 / 19 / 15 / 16). It returns fusions at f <= .5 because it cuts
  the minority links into relation records: of the 4 / 7 / 4 / 1 F1 bridge junctions at f = .1 / .5 / .9 / 1.0 it cuts 4 / 6 / 1 / 1 and leaves the others (share >= 1/2) in the locus, 3 of 4 at f = .9; M1 is 9 / 5 / 0 / 1 at f = .1 / .5 / .9 / 1.0 (`bridges.tsv` `keep` True counts).
- Current defaults: M2 is 0 at every f < 1 (no partner enters NPIP) and 1 at f = 1.0 (extended; 3 on the container-sim pool: the dominant fusion takes the partner into the family). M4 is 0 in every arm: where the fused locus survives the partner's holder is the fused member's own locus (extended pool 3 / 9 / 9 of 10
  at f = .5 / .9 / 1.0; container-sim pool 8 / 8 / 10), for which M4 is 0 by construction (prereg Outcome); where F1v2 separated the pair (extended f = .1, 9 of 10; most of f = .5) the partner's unit is clustered elsewhere and unrelated to NPIP.
- M3 / M4 are not a container signal under the current defaults: where F1v2 removed the bridge (f = .1, extended) no fused member has any accessory partner base (M3 0 / 0), and where the fusion stays in the locus (f >= .5) the members that are still fused are the ones container v1 already failed on
  (extended pool M3 precision .19 / .31 / .63 and recall .12 / .42 / .44 at f = .5 / .9 / 1.0; the record's Outcome on the container-sim pool had .21 / .42 at f = .5 and .9 and .39 / .62 at f = 1.0).
- **The containment default costs unfused NPIP copies in this simulation.** Unfused controls in NPIP /15: pre-flip 14 / 15 / 15 / 15 / 15; current defaults 13 / 13 / 13 / 12 / 13. The isolating pairs: containment alone (pre-flip vs `bc`) 14 / 15 / 15 / 15 / 15 → 13 / 14 / 13 / 12 / 13; containment on top of F1v2 (`gc` vs `cur`) 14 / 14 / 14 / 15 / 15 → 13 / 13 / 13 / 12 / 13; F1v2 alone
  (pre-flip vs `gc`) moves NPIPB10P (outside NPIP at f = .1 and .5 in `gc`, inside in `pre`). NPIPB3 and NPIPB6 are outside NPIP in every `cur` arm and in every `bc` arm except NPIPB6 at f = .1 (inside), NPIPB5 leaves at f = .9 in both, and NPIPB10P, outside in `pre` at f = 0, is inside in `cur`.
  On real human chr16 the containment-only arm equals B0 (Bc = B0 = 22, D = 24; NPIPB3 stays in NPIP in both, NPIPB6 is outside under B0 and Bc and inside under D and Gc), and on the real
  substrates the containment default added no NPIP or TBC1D3 copy and lost none; the cost is therefore a simulation fact, not measured on a real library, and its mechanism was not examined. It is reported, not judged.
- M5: at f = 0 the fused copies' units carry 8,950-8,956 accessory bp under containment .70 (132 pre-flip): NPIPA7's unit, already outside NPIP at f = 0, joins a cluster to which its
  NPIP exons are accessory.

## 6. Ranked failure mechanisms (deliverable 5)

Counted per species over the observations of §2 (current defaults D; holder-rule instrument). (a) unit not separated at node formation; (b) placement by representative; (c) container / relation definition; (d) admission /
coverage gate; (e) not fusion-related. "Would a unit-aware container change it": the oracle arms of §4 say whether the copy changes family, which is the most a container that finds
families on units can do for membership. Ranked by the number of observations not in their family; (c) concerns observations that are in it.

| rank | class | human (68 observations; 28 not in family) | gorilla (78 observations; 53 not in family) | examples | unit-aware container changes the outcome? | at most |
|---|---|---|---|---|---|---|
| 1 | (e) no locus | 22 | 52 | human: 12 with >= 2 same-strand primary reads and no assembled locus (testis 9, chr17 3), 10 with < 2 (testis 10); gorilla: 30 with >= 2 (KB NPIP 15, OR NPIP 12, OR TBC1D3 3), 22 with < 2 (KB TBC1D3 14, OR TBC1D3 4, OR NPIP 2, KB NPIP 2) | no: there is no locus to split (fusion involvement not examined) | 0 |
| 2 | (e) fragment or touch-only holder | 2 | 1 | testis NPIPB4, NPIPB5 (529-bp fragments, no alignment to the family); OR NPIPB7 (fused, holder overlaps 1.4%) | no | 0 |
| 3 | (d) admission gate | 2 | 0 | PKD1P6-NPIPP1 (`no_exonic` against the three members it aligns to; a truth-labelling question is open, §8); TBC1D3P5 (`low_shared` against all 11 evaluated members) | no (oracle A, C, S leave both unchanged) | 0 |
| 4 | (a) unit not separated | 2 | 0 | NPIPB5 (158-kb NPIP + SMG1P block, one locus, a hub over four clusters); LOC100420311 (fused with TBC1D29P; its fused exon model fails `low_shared`) | yes: both become right under A, C and S | +2 observations (+3.8 and +6.3 points of sensitivity on chr16 and chr17); 0 in gorilla |
| - | (b) placement by representative | suspected at NPIPB5 (same observation as (a)) | 0 | NPIPB5: the locus's representative is the SMG1P1 transcript (73 reads) and the locus goes to the SMG1P cluster; the locus has admitted edges into both clusters | same observation; not isolated from (a), not tested | (inside the +2: it is the NPIPB5 observation) |
| - | (c) container / relation definition | 11 fused observations already in their family (chr16 8, chr17 2, testis 1) | 1 (LOC129533797) | NPIPA1, NPIPA6, NPIPA7, NPIPA9, NPIPB3, LOC128966608, NPIPB4, NPIPB14P; TBC1D3G, TBC1D3; testis NPIPB14P; LOC129533797; partner bases in core blocks 77% of those in a block (chr16) | membership: none. Relations: yes; on chr16 pieces of 12-14 of 83 fused loci land in 2+ families, 5-6 involve NPIP (oracle count, unscored) | 0 observations; relations only |

Notes.
- (b) is not isolated by these data. The representative sets the exon model on which every edge of the locus is tested (`npf_audit`'s finding), but NPIPB5's locus also aligns end to end to the two
  fused NPIP loci (99.2% over 117,861 bp) and keeps admitted edges into MCL1, so MCL is choosing between two families that the duplicated NPIP + SMG1P block belongs to. What the oracle shows is that once
  the pieces are separate nodes the NPIPB5 piece is placed in NPIP; which of the two choices (representative or block) is the cause was not tested by changing the representative.
- Why F1v2 leaves NPIPB5's locus whole: the locus has one F1 bridge candidate, inside SMG1P1 (188 transition reads, UP 166, DOWN 149, share .5579 >= 1/2, not cut), and NPIPB5's own link to SMG1P1 (L4, minority, 2 transcripts) has no bridge candidate. Whether a read-share rule at the NPIPB5 boundary would separate the copy
  is untested.
- Genome-wide (Compara, A119b, 19 chromosomes, 1,062 multi-copy genes; holder rule): no locus 412 (38.8%) > unclustered 216 (20.3%) > other cluster 34 (3.2%); the fused share of the present-but-unplaced 250
  is 41 (16.4%, conditioned on being unplaced) and the partner-representative cases are 10 (8 unplaced, 2 placed).

## 7. What a unit-aware container would change, at most

- **Membership (holder rule), per species.** Human: +2 of 68 observations (+1 of 26 NPIP on chr16, +1 of 16 TBC1D3 on chr17, 0 in testis). Gorilla: 0 of 78 (all four cells). No cross-species total. Compara multi-copy genes (human):
  at most 41 of 1,062 (3.9%) on the 19 A119b chromosomes, 1 of 1,313 in testis; the oracle realises +24 true pairs (+27 on chr1, -3 on chr19) and raises the pooled-truth F by .013 with two chromosomes split.
- **Family scores (scorer).** Compara F: chr16 +.032 (A, C) to +.047 (S), chr1 +.047, chr19 +.042 to +.054, chr17 0, pooled 19 chromosomes +.013 (two split); U2 ID_154 -.019 (A, C) to +.031 (S); Liftoff
  recall 0 to -2 pairs on chr16. Part of the Compara gain is the scorer's collision count (§4.2).
- **Relations.** The one thing a unit-aware container adds that no membership metric scores: 12 to 14 of 83 fused loci on chr16 (5 to 6 involving NPIP) and 5 to 9 of 97 on chr17 (3 involving TBC1D3)
  have constituents in two or more families, where container v1 leaves partner pieces as core or in no block (4 of the 8 fused in-family chr16 observations have an empty container); the relation count is an oracle assignment, not scored against a truth.
- **What it cannot change.** The no-locus observations (human 22, gorilla 52), the gate failures and the copies whose holder is a fragment or touch-only. These are 26 of the 28 human and 53 of the 53 gorilla observations that are not in their family.
- **Envelope beyond units.** If every observation that has a locus were in its family (a perfect family assigner, which no container is), membership would rise by 6 of 68 in human (the 2 above, the 2 gate failures and the 2
  testis fragments) and by 1 of 78 in gorilla (the touch-only holder); the no-locus observations would stay out. The 4 beyond the units in human are edge problems (exonic evidence for PKD1P6-NPIPP1, a shared-exon fraction for TBC1D3P5, alignments for
  two fragments), not container problems.

## 8. Limits and traps

- **The oracle reads the annotation** (r1017's own caveat). A clips to annotated spans (it cut about 2 kb of read-supported NPIP exon outside RefSeq's span at NPIPB7 in r1017; here NPIPB7 stays in NPIP in
  every oracle and on D); C keeps read extents. Both recover the same two copies; the family-level numbers differ by .012 to .021 between them. The oracle's piece names also make the scorer's gene-to-locus mapping exact (collapsed genes 3 → 1 → 0), which is a scoring gain and not family finding.
- **Coverage of the oracle.** chr16, chr17, chr19, chr1, testis chr16, both gorillas (OR: A only; KB: A and C); the statement "A and C agree" holds where both were run (chr16, chr17, chr19, testis, KB3781). Chromosomes were split where fused loci concentrate, so the genome-wide gain is an estimate of unknown sign (§4.1).
- **"Right family" is presence.** The holder shares >= 1 same-strand exon base with the copy; it covers < 50% of the copy in 22 of the 65 right-family observations (chr16 1, chr17 1, testis 4; OR6737 NPIP 10, KB3781 NPIP 6). A 50% floor on the right-family column moves 6 of the 40 human and 16 of the 25 gorilla right-family observations into a "present, holder covers < 50%" bucket, leaving human 34 right / 6 wrong / 6 low-coverage / 22 no locus and gorilla 9 / 1 / 16 / 52; three of the seven wrong-cluster observations are holders covering 1.4 to 14%.
- **The no-locus count is strand-specific.** A copy is scored on its own strand. Ignoring strand would change the holder of 13 copies (testis 5, OR6737 NPIP 3, KB3781 NPIP 5; 12 of the 74 have an
  opposite-strand locus only) [ver], and 4 of those 12 (OR6737 NPIPB8, NPIPB2; KB3781 NPIPB4, NPIPB2) would land in the NPIP family through an opposite-strand locus of 3.1-8.6 kb on their territory (a neighbouring
  inverted member or antisense transcription; not examined). In human testis a strand-blind reading would put 7 copies in MCL2 (each of its seven `+`-strand loci overlaps a copy, five of them antisense) against 5 in MCL3, so the argmax family cluster would flip to MCL2; the strand-specific reading is the one used.
- **Holder-based fate.** One holder per copy, `npf_audit`'s rule; no holder is shared by two copies and ties never fire [ver]. NPIPB5's territory is 5.6% covered by family-cluster loci although its holder is the big SMG1P1 locus; the touch columns are in `out/fate_table.tsv`.
  Copies that sit in several loci or clusters are not double counted. The instrument disagrees with the family scorer on single copies (§3.7).
- **Truth proxies.** Gorilla NPIP is the 2026-09-17 T_member proxy, human NPIP the Dishuck-checked set (NPIPB14P as the readthrough exons inside its span), TBC1D3 the
  RefSeq records (five human pseudogene spans without exon features use the span as one exon). PKD1P6-NPIPP1's territory is chr16:15,126,651-15,141,806 (1-based; 20 of the readthrough record's 30 merged exon segments; `truth.hsa.json` h03, whose territory start is the 0-based 15,126,650), which lies entirely inside the RefSeq PKD1P6 record (15,126,099-15,159,720); the GFF has no separate NPIPP1 record, and whether the territory is the NPIPP1 half was not re-derived here.
  Dishuck calls PKD1P6-NPIPP1 and the readthrough copies NPIP, while Soto and U2 file them in ID_149 (the PKD1 block) (`npf_critique` §5.2).
- **Instrument blind spot and thresholds.** `a.fused` drops readthrough-described records, so NPIPA9 (403 fusion reads to PKD1P5-LOC105376752) is fused under R only; "12 fused observations already in their family" is the R reading (P or L gives 10). Every fused reading flags a second record at >= 1 bp of exon overlap (NPIPA7 is fused by 53 bp of PKD1P2). The P and L columns can differ by up to 25 loci of 1,550
  on the 19 chromosomes; the lost-member counts are identical under both. Register 1193 counts 13 fused human NPIP copies for our models (11 excluding the two readthrough-defined copies); this document counts 9 fused NPIP holders on chr16 (R reading: 8 under P and L, plus NPIPA9) after F1v2's removal of bridge transcripts; model chains and locus holders are different instruments and were not reconciled copy by copy.
- **Per-chromosome families** (A119b) have no cross-chromosome edges (declared in o1_cover); the pooled 19-chromosome scores use a truth that contains 352 cross-chromosome pairs (§3.4), so they are comparable across arms, not to a genome-wide run or to the per-chromosome rows. chr13 is not in the 19 chromosomes.
- **Gorilla five-contig runs** are validated on KB BASE (84 of 86 cluster sets) and by the stored NPIP / TBC1D3 clusters lying inside the contigs; OR's B0 is the stored genome-wide BASE, not a re-run (†).
- **Libraries of different depth.** The no-locus share is a depth statement (gorilla KB3781 TBC1D3 14 of 14 with < 2 primary reads; testis 19 of 26 NPIP copies without a locus; human A119b 3 of 42). Per-species counts add cells of different depth.
- **Simulation.** n = 10 fusions, annotation-built reads (circular by construction), two pools that differ in the read-end model; M1 counts move by whole fusions.
- **Everything is in-sample** for the two flips (Conventions). Both recovered copies (NPIPB5, LOC100420311) are in the development block (human chr16, chr17); the substrates the fusion rules were not developed on (human testis, both gorillas) show 0 recoveries, which is a descriptive contrast, not a held-out test.

## 9. Draft register rows (suffix H; not appended)

| # | date | area | claim | verdict |
|---|---|---|---|---|
| 1194H | 2026-09-30 | container headroom | Under the current defaults (`assemble` f1v2, `families` `--min-cov-shorter 0.70`) a unit-aware container (families found on units, extra pieces kept as relations) has material headroom for FINDING NPIP / TBC1D3 members | **NO, per species: human 2 of 68 observations, gorilla 0 of 78 (holder-rule instrument).** Fates of the 146 copy-by-cell observations (81 distinct copies; seven cells, species never added): human 40 in family / 6 in another cluster / 22 no locus (68); gorilla 25 / 1 / 52 (78). The oracle (every fused locus split into one node per constituent annotated gene; A and C agree where both were run) recovers NPIPB5 (human chr16, 24 → 25 of 26) and LOC100420311 (human chr17, 11 → 12 of 16) and nothing in human testis or either gorilla. The other five wrong-cluster observations fail admission gates (PKD1P6-NPIPP1 `no_exonic`, TBC1D3P5 `low_shared`) or are fragments / touch-only holders (testis NPIPB4, NPIPB5; gorilla OR6737 NPIPB7). "Right family" is presence: the holder covers >= 50% of the copy in 34 of the 40 human and 9 of the 25 gorilla right-family observations. `docs/CONTAINER_HEADROOM_2026-09-30.md` |
| 1195H | 2026-09-30 | container headroom | Fused loci are the main reason NPIP / TBC1D3 members are missing | **NO: 74 observations have no locus at all** (human 22 of 68, gorilla 52 of 78; the species are not added): 42 have >= 2 same-strand primary reads and no assembled locus (floor / gate; human 12, gorilla 30), 32 have < 2 same-strand primary reads (human 10, gorilla 22). A container cannot reach them; fusion involvement was not examined for any of the 74. Another 3 present observations are fragments or touch-only holders (no alignment to any family member). |
| 1196H | 2026-09-30 | container headroom (flip decomposition) | The 09-29 defaults flip (F1v2 + `--min-cov-shorter 0.70`) left the fusion-related NPIP headroom untouched | **PARTLY: it moved 2 of the 4 copies that change fate in any arm**: chr16 22 → 24 of 26 (NPIPB2 and NPIPB6, through the regroup of F1v2; link reads L = 0, no bridge cut), and containment added 0 NPIP / TBC1D3 copies on every cell that has a containment-only arm (none for OR6737, where D and Gc agree). Compara F (scorer) decomposition (pre-flip → F1v2 alone → current): chr16 .615 → .667 → .667, chr17 .412 → .460 → .475, 19 A119b chromosomes .445 → .451 → .459 (pooled truth), testis .264 → .265 → .265. The fusion-related copies still lost are NPIPB5 and LOC100420311. |
| 1197H | 2026-09-30 | container headroom (oracle) | A perfect split of fused loci raises family-level scores materially, and makes NPIP worse (r1017) | **PARTLY: only where fused loci hold members, and the NPIP result depends on the truth and the instrument.** Compara F +.032 (chr16, A and C; +.047 with r1017's selection), +.047 (chr1), +.042 / +.054 (chr19), 0 (chr17); pooled truth over the 19 A119b chromosomes .459 → .472 with only chr1 and chr19 split (chromosomes chosen where fused loci concentrate, so the sign elsewhere is unknown). Copies in the NPIP cluster 24 → 25 of 26 under A, C and S (holder rule). Soto ID_154 F .867 → .929 and Compara CF153 .812 → .848 (no 'NPIP worse' there); U2 ID_154 F .769 → .750 (A, C) / .800 (S), r1017's direction, small. Liftoff recall 17 → 15 (A) / 17 (C) / 16 (S) of 36. Part of the Compara gain is the scorer's collision count (collapsed genes 3 → 1 → 0). r1017's mechanism was not re-examined. Oracle rebuilt for the `--from-gtf` GTF interface; same aligner call in every arm (r1018). |
| 1198H | 2026-09-30 | container headroom (genome-wide) | Fused loci are where multi-copy (Compara) members are lost genome-wide | **NO (minor), holder rule: 41 of 1,062 (3.9%)** Compara multi-copy genes of the 19 A119b chromosomes are present, unplaced and in a fused locus (the same 41 under the per-transcript and locus-union readings), against 412 with no locus and 216 unclustered; 34 are in another cluster (6 fused). Testis 1 of 1,313. All 10 genes whose locus representative is a partner's transcript (8 of the 250 present-but-unplaced, 2 placed right) are in fused loci. 40 of the 400 correctly placed genes also sit in fused loci. Fused loci are 1,550 (locus union) / 1,525 (per transcript) of 75,292 loci (2.1%). The scorer instrument finds more losses on chr16 (CF153 misses 6 of 19 genes that the holder rule places). |
| 1199H | 2026-09-30 | container headroom (real data) | Real loci that join two multi-copy families are F1v2 bridges (extends 1185) | **NO.** Loci joining two multi-copy Compara families: 3 (chr16), 15 (19 chromosomes: 5 dominant, 5 mixed, 5 minority, all present in `families.gtf`), 3 (testis); Liftoff: 7 + 1 + 1 (OR6737) + 0 (KB3781). F1v2 relation records joining two multi-copy families: 0 on the 19 chromosomes and in testis under Compara or Liftoff, 1 on chr16 under U2 (NPIPA6 \| NPIPA7 \| PKD1P1: one relation of two transcripts, 6 reads). The links that join multi-copy families are dominant, mixed or minority and stay in their loci. |
| 1200H | 2026-09-30 | container v1 (current defaults) | Container v1 expresses the partner of a fused member that is in its family | **NO: mostly core or in no block.** 12 fused observations (human 11, gorilla 1) are already in their family; of the exon bases of their partner records 28,293 (chr16) fall in core blocks, 8,327 in accessory blocks and 76,387 in no block (chr17 1,951 / 515 / 515; testis 846 / 61 / 3,287); 77% of the bases that are in a block are core (by construction: core = aligned to another member). The container is empty at NPIPA1, A6, A9 (PKD1P pieces align to other members), NPIPA7 (after F1v2's split its locus holds 4 link reads to PKD1P2; 53 of 9,905 PKD1P2 bases in any block) and TBC1D3G. It is informative where the partner does not align to a family member (NPIPB4~RRN3P1 1,609 of 1,609 bases accessory; gorilla LOC129533797 885 of 1,603 accessory bp), and its accessory set is mostly not partner bases (NPIPB4 16%). The container Outcome of 09-28 holds on the current defaults. |
| 1201H | 2026-09-30 | container headroom (simulation) | Dev baseline of the fusion simulation under the current defaults (prereg §3 measures) | **MEASURED.** Extended pool: fused copies in NPIP at f = 0 / .1 / .5 / .9 / 1.0 = 9 / 9 / 5 / 0 / 1 (pre-flip 9 / 7 / 1 / 0 / 1; F1v2 alone 9 / 9 / 5 / 0 / 1), so F1v2's dev numbers 9 / 5 / 0 are reproduced by the Rust default; copies in NPIP /25 22 / 22 / 18 / 12 / 14. Container-sim pool: 9 / 2 / 0 / 0 / 3 (pre-flip 9 / 1 / 0 / 0 / 4). M2 (partners in NPIP) 0 for f < 1 and 1 (extended) / 3 (container-sim) at f = 1.0; M3 precision / recall .19 / .12, .31 / .42, .63 / .44 at f = .5 / .9 / 1.0 (extended); M4 0 by construction where the partner's holder is the fused member. Prereg bars failed by the current defaults: M1 at f = .5 (extended; container-sim also f = .1), M2 at f = 1.0, M3 at f >= .5. |
| 1202H | 2026-09-30 | container headroom (simulation, side result) | The `--min-cov-shorter 0.70` default is free for NPIP | **PARTLY: on real substrates yes, on the extended simulation no.** Real: 0 NPIP / TBC1D3 copies changed between Gc and D on any of the seven cells, nor between B0 and Bc on the five cells that have a Bc arm. Simulation: unfused NPIP controls in the family 14 / 15 / 15 / 15 / 15 (pre-flip) → 13 / 13 / 13 / 12 / 13 (current defaults); containment alone (pre-flip vs `bc`) 13 / 14 / 13 / 12 / 13, containment on top of F1v2 (`gc` vs `cur`) 14 / 14 / 14 / 15 / 15 → 13 / 13 / 13 / 12 / 13; F1v2 alone moves NPIPB10P. NPIPB3 and NPIPB6 are outside NPIP in every `cur` arm (and in `bc`, except NPIPB6 at f = .1) and NPIPB5 at f = .9. Reported, not judged; the mechanism was not examined (the NPIP-guided regression of r1007 / r1009 is the known relative). |
| 1203H | 2026-09-30 | instrument | `a.fused` counts every fusion at NPIP | **NO: readthrough-described RefSeq records are outside its gene set**, so NPIPA9 (403 fusion reads to PKD1P5-LOC105376752) is fused under the audit reading only. Human chr16 holders: 9 fused NPIP copies (8 under the per-transcript and locus-union readings, 9 with readthrough records). chr16 fused loci 83 (locus union) / 82 (per transcript) of 2,832, equal to the F1v2 dev numbers; register 1147's warning to name the reading applies (the lost-member counts are the same under both readings). Every reading flags a second record at >= 1 bp of exon overlap. |
| 1204H | 2026-09-30 | instrument (replay) | An all-vs-all PAF computed with bridge loci present can stand in for the current defaults' all-vs-all | **YES, byte for byte on every check made**: the ALL-run PAF minus records with a bridge locus at either end, fed to the current binary through a stand-in `RUSTLE_MINIMAP2`, gives `clusters.tsv` identical to the real run on chr17 (D and Gc), identical to the stored BASE / COVER core on 19 of 19 A119b chromosomes (B0, Gc), and identical to the real genome-wide current-default testis run (D). The gorilla five-contig restriction reproduces the stored genome-wide BASE on 84 of 86 cluster sets (KB) and leaves the stored NPIP and TBC1D3 clusters whole; `tools/mm2_shard.sh` equals direct minimap2 on testis chr16. |
| 1205H | 2026-09-30 | instrument (holder rule vs family scorer) | Holder-rule copy recovery and family-scorer gene hits are interchangeable on human chr16 NPIP | **NO.** The steps agree in size (B0 → D: +2 copies by the holder rule, +3 hits under each of Compara CF153, U2 ID_154 and Soto ID_154; D → oracle A: +1 copy, +1 CF153 hit, 0 under U2 and Soto) but not in identity (holder rule: NPIPB5; CF153: NPIPB3). Under D the scorer counts NPIPB5 as a CF153 hit (holder rule: wrong cluster) and misses six CF153 genes (NPIPA1, A6, A9, B13, B3, B4) and four U2 genes (LOC124907834, NPIPA6, NPIPA9, NPIPB14P) that the holder rule places in the family: the 09-21 one-name-per-locus trap. Every lost-member count names its instrument. |
| 1206H | 2026-09-30 | container headroom (F1 attribution) | F1v2 leaves NPIPB5 fused because its read-share test refuses to cut the NPIPB5 link | **NO (attribution).** The locus's one F1 bridge candidate (chr16:22,753,211-22,756,690; 24 transition transcripts, 188 reads, UP 166, DOWN 149, share .5579, not cut) lies inside SMG1P1, 25 kb upstream of NPIPB5; NPIPB5's own link to SMG1P1 (L4, 2 transcripts, minority) has no bridge candidate. Whether a read-share rule at the NPIPB5 boundary would separate the copy is untested. In `bridges.tsv` `keep` True means the bridge IS cut. |

## 10. Provenance and source keys

Scratch `/mnt/linuxdisk/tmp/rustle_figures_dev/container_headroom/` (`code/`, `data/`, `out/`, `fsim/`, `fsim_ext/`, `NOTES.md`; 24 GB before and 7.5 GB after the clean-up of regenerable loci FASTAs, PAFs (four `D` PAFs kept) and PAF caches; this run's nine entries of the shared `mm2_shard_cache` were removed, so a re-run re-aligns). Repository `main` at `e163d955`.

**Binaries** (copied to `bin/`, `SHA1SUMS`): `copy_assign` 1325e9d1, `mcl_families` 67b2d40f, `as_table` 0452b1b9, `family_score` 27aa9445 (built 2026-09-29 23:00-23:03; no `src/` file newer).
minimap2 2.30-r1287, samtools 1.22.1, python3 3.14.4 (scipy, pysam). The shared heavy lock was free throughout; every heavy call went through `tools/rlock.sh heavy`, the table rebuilds through `tools/rlock.sh light`.
While this document was being revised other agents left uncommitted edits in the working tree (`src/bin/copy_assign.rs`, `src/bin/mcl_families.rs`, `src/rustle/vg_family/{bridge_regroup,mod}.rs`, `bench/soto/soto_replication.py`, a new `src/rustle/vg_family/family_relations.rs`); none of it was used. Every product here comes from the binaries copied to scratch `bin/` (sha1 above), and the one source statement quoted (the `keep` semantics of `bridges.tsv`, `bridge_regroup.rs` lines 32 and 673) was checked against `git show HEAD`.

**Source keys.**
- [fate] `out/fate.<cell>.<arm>.json`, `out/fate_table.tsv`, `code/hl_core.py`, `code/fate.py`, `code/assemble_tables.py`.
- [cell] `out/cell_summary.json`. [mech] `out/mech.<cell>.D.json`, `out/mech_summary.json`, `code/mech.py`. [reads] `out/reads_absent.json`.
- [sc16] `out/scores.human.chr16.json`; [sc17] `.chr17.json`; [sc1] `.chr1.json`; [sc19] `.chr19.json`; [sca19] `.merged.json`; [sca19o] `.a19_oracle.json`; [sctes] `.testis.json`;
  [scgor] `out/scores.gorilla.json`; [pairs16] `out/compara_pairs_chr16.json`; per-family hit sets `out/fs/*.pf.tsv`. Code `code/score_human.py`, `score_fam.py`, `score_gorilla.py`.
- [gf] `out/gene_fates.a19.json`, `.testis.json`, `code/gene_fates.py`; [fused] `out/fused_counts.a19.json`; [fl] `out/fused_links.json`, `code/fused_links.py`; [jn] `out/joins.human.json`, `.gorilla.json`, `code/joins.py`; [rel] `out/oracle_relations.json`.
- [sim] `fsim/score_{pre,cur}/score.out|json`, `fsim_ext/score_{pre,gc,bc,cur}/score.out|json`, `code/fsim_score.py`, `code/fsim_run_arm.sh`.
- [gates] `out/gates.txt` (the recorded comparisons), `code/cmp_gtf.py`, `cmp_part.py`, `replay_all.sh`, `prep_replay.py`, `mm2_replay.sh`.
- Oracle: `code/oracle.py`. Reused unchanged: `o1_cover_frozen/{score.py, lo_score.py, emu.py}`, `f1_frozen/rg3_lib/{npf.py, audit.py}`, `fusion_container_sim/score.py` 8f512fb8 (copied with paths changed), truths
  `copy_recovery_tools/ann/truth.{hsa,ggo}.json`, `ggo_npip_sim/ann/truth.json`, `families_gw/species/human/{compara.Primates,npip_u2,soto}.families.tsv`, `gw22/spectrum/compara_chr16.tsv`, `liftoff/<species>/liftoff_loci.tsv`.
  Annotations: human `/mnt/linuxdisk/tmp/regress/chm13.gff` (= `winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz` decompressed, md5 691d676f); gorilla `winloci_data/GGO_genomic.gff`.
- [ver] `figs/container_headroom_verify.md` and `/mnt/linuxdisk/tmp/rustle_figures_dev/container_headroom_verify/` (own parser and scoring code; claims 1-6 of the brief reproduced exactly, 0 FAIL; it also ran the frozen
  `run_family_score` with the copied binary and with `fj_bin_frozen/family_score` 7723029b, identical output on the eight Compara checks). [rev] `figs/container_headroom_review.md` (hostile read of the first full version; its findings are applied here). [ver2] `figs/container_headroom_verify2.md` and `/mnt/linuxdisk/tmp/rustle_figures_dev/container_headroom_verify2/` (own code; recomputation of the claims that changed in the revision).
- [CRT] `docs/` and scratchpad `figs/copy_recovery_tools.md` (2026-09-30). Earlier records: `docs/PREREG_fusion_container_sim_2026-09-28.md`, `PREREG_f1v2_readshare_2026-09-29.md`, `PREREG_o1_cover_growth_2026-09-29.md`,
  `PREREG_overmerge_ceiling_2026-09-22.md`, register rows 845, 846, 967, 971, 1013-1018, 1134, 1145-1151, 1184-1193.

## Appendix A. Per-copy fate table (146 observations, current defaults D)

Columns: cell (h16 human chr16 NPIP, h17 human chr17 TBC1D3, tes human testis NPIP, OR-N / OR-T gorilla OR6737 NPIP / TBC1D3, KB-N / KB-T gorilla KB3781); the copy's other truth families (Compara multi-copy family, U2, Soto; every copy is also in the Dishuck NPIP / RefSeq TBC1D3 set by construction; gorilla has none); locus gene_id and node (`contig:start-end`); fused P / L / R;
strongest link; cluster and size, with the labels of the cluster's members when it is not the family cluster; the share of the copy's territory the holder covers; fate under B0 > D > oracle A / C (right, wrong, uncl = unclustered, none = no locus; `-` = arm not run); mechanism; container v1 of the holder unit.

| cell | copy | other truth families of the copy | locus gene_id | node (contig:start-end) | fused P/L/R | strongest link (partner: class L/C/P [F1 status]) | cluster (size); what the cluster is, when it is not the family | holder covers (share of the copy) | fate B0 > D > oracle A/C | mechanism | container v1 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| h16 | NPIPB2 | Compara CF153; U2 ID_154; Soto ID_154 | DN_chr16_11906177_15.rg2 | chr16:11963295-11977765 | --- | - | MCL1 (30) | 0.698 | uncl > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h16 | NPIPA2 | Compara CF153; U2 ID_154; Soto ID_154 | DN_chr16_14749526_8 | chr16:14740775-14763868 | --- | - | MCL1 (30) | 0.847 | right > right > right/right | - | acc 707 bp; partner core 0 / acc 0 |
| h16 | NPIPA1 | Compara CF153; U2 ID_149; Soto ID_149 | DN_chr16_14938768_7 | chr16:14924775-14953128 | PLR | PKD1P3-NPIPA1:dominant L193/C2/P17 [struct up/-] | MCL1 (30) | 1.0 | right > right > right/right | (c) in family; partner bases carried: core 3715 / container 0 | acc 0 bp; partner core 3715 / acc 0 |
| h16 | PKD1P6-NPIPP1 | U2 ID_149; Soto ID_149 | DN_chr16_15127017_8 | chr16:15127011-15137722 | --- | - | MCL25 (3): PKD1P6-NPIPP1,PKD1P6,PKD1P2 | 0.541 | wrong > wrong > wrong/wrong | (d) alignments to the family exist, every one rejected: no_exonic | acc 2409 bp; partner core 0 / acc 0; rel MCL1,MCL26 |
| h16 | NPIPA5 | Compara CF153; U2 ID_154; Soto ID_154 | DN_chr16_15368421_8 | chr16:15368407-15401113 | --- | - | MCL1 (30) | 0.851 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h16 | NPIPA6 | Compara CF153; U2 ID_154 | DN_chr16_16402040_2.rg2 | chr16:16329977-16359025 | PLR | PKD1P1:dominant L401/C128/P0 [no-row] | MCL1 (30) | 1.0 | right > right > right/right | (c) in family; partner bases carried: core 4473 / container 0 | acc 0 bp; partner core 4473 / acc 0 |
| h16 | NPIPA7 | Compara CF153; U2 ID_154; Soto ID_154 | DN_chr16_16402040_2 | chr16:16389357-16406194 | PLR | PKD1P2:mixed L4/C180/P0 [struct -/down] | MCL1 (30) | 1.0 | right > right > right/right | (c) in family; partner bases carried: core 53 / container 0 | acc 0 bp; partner core 53 / acc 0 |
| h16 | NPIPA8 | Compara CF153; U2 ID_154; Soto ID_154 | DN_chr16_18325162_2 | chr16:18325158-18341989 | --- | - | MCL1 (30) | 0.795 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h16 | NPIPA9 | Compara CF153; U2 ID_154 | DN_chr16_18372354_2 | chr16:18372351-18401447 | --R | PKD1P5-LOC105376752:dominant L403/C126/P0 [no-row] | MCL1 (30) | 1.0 | right > right > right/right | (c) in family; partner bases carried: core 3841 / container 0 | acc 0 bp; partner core 3841 / acc 0 |
| h16 | NPIPB3 | Compara CF153; U2 ID_151; Soto ID_151 | DN_chr16_21337422_3 | chr16:21337406-21455266 | PLR | LOC100190986:mixed L10/C357/P2 [struct -/-] | MCL1 (30) | 0.998 | right > right > right/right | (c) in family; partner bases carried: core 7316 / container 617 | acc 2040 bp; partner core 7316 / acc 617; rel MCL26,MCL77,MCL84 |
| h16 | LOC128966608 | U2 ID_151 | DN_chr16_21680693_3 | chr16:21640357-21780522 | PLR | LOC128966632:minority L67/C134/P256 [struct up/down] | MCL1 (30) | 0.998 | right > right > right/right | (c) in family; partner bases carried: core 6269 / container 1401 | acc 2830 bp; partner core 6269 / acc 1401; rel MCL26 |
| h16 | NPIPB4 | Compara CF153; U2 ID_152; Soto ID_152 | DN_chr16_22374851_3 | chr16:22350725-22422849 | PLR | RRN3P1:minority L2/C351/P131 [struct -/down] | MCL1 (30) | 1.0 | right > right > right/right | (c) in family; partner bases carried: core 0 / container 1609 | acc 10044 bp; partner core 0 / acc 1609 |
| h16 | NPIPB5 | Compara CF153; U2 ID_153; Soto ID_153 | DN_chr16_22714381_13 | chr16:22696405-22854661 | PLR | SMG1P1:minority L4/C32/P456 [F1 zone inside partner, not cut (share 0.5579)] | MCL26 (4): SMG1,SMG1P6,LOC124903796 | 0.658 | wrong > wrong > right/right | (a) fused locus (the oracle recovers it); (b) suspected, untested: the representative is the partner transcript | acc 17886 bp; partner core 4389 / acc 3825; rel MCL1,MCL12,MCL14,MCL25,MCL3,MCL5,MCL6,MCL64,MCL75,MCL76,MCL77,MCL84 |
| h16 | NPIPB6 | Compara CF153; U2 ID_154; Soto ID_154 | DN_chr16_28659942_21.rg2 | chr16:28623010-28637560 | --- | - | MCL1 (30) | 0.601 | wrong > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h16 | NPIPB7 | Compara CF153; U2 ID_154; Soto ID_154 | DN_chr16_28736995_2 | chr16:28736992-28772243 | --- | - | MCL1 (30) | 0.786 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h16 | NPIPB8 | Compara CF153; U2 ID_154; Soto ID_154 | DN_chr16_28934209_2 | chr16:28924759-28939456 | --- | - | MCL1 (30) | 0.999 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h16 | NPIPB9 | Compara CF153; U2 ID_154; Soto ID_154 | DN_chr16_29048209_2 | chr16:29038768-29053456 | --- | - | MCL1 (30) | 0.526 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h16 | NPIPB10P | U2 ID_154; Soto ID_154 | DN_chr16_29319348_7 | chr16:29319291-29333927 | --- | - | MCL1 (30) | 0.806 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h16 | NPIPB11 | Compara CF153; U2 ID_154; Soto ID_154 | DN_chr16_29663336_7 | chr16:29663337-29679826 | --- | - | MCL1 (30) | 0.531 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h16 | NPIPB12 | Compara CF153; U2 ID_154; Soto ID_154 | DN_chr16_29776171_2 | chr16:29768895-29789163 | --- | - | MCL1 (30) | 0.347 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h16 | LOC124907834 | U2 ID_154 | DN_chr16_30507437_3 | chr16:30507421-30531871 | --- | - | MCL1 (30) | 1.0 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h16 | NPIPB13 | Compara CF153; U2 ID_155; Soto ID_155 | DN_chr16_30620297_2 | chr16:30609375-30626632 | --- | - | MCL1 (30) | 0.507 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h16 | NPIPB14P | U2 ID_154; Soto ID_154 | DN_chr16_75785724_2 | chr16:75785704-75876271 | PLR | PDXDC2P-NPIPB14P:dominant L94/C0/P8 [struct up/down] | MCL1 (30) | 1.0 | right > right > right/right | (c) in family; partner bases carried: core 2626 / container 4700 | acc 5178 bp; partner core 2626 / acc 4700 |
| h16 | NPIPB15 | Compara CF153; U2 ID_154; Soto ID_154 | DN_chr16_80204820_2 | chr16:80194084-80209901 | --- | - | MCL1 (30) | 1.0 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h16 | LOC124907808 | U2 ID_154 | DN_chr16_80319195_2 | chr16:80308466-80324276 | --- | - | MCL1 (30) | 0.962 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h16 | LOC124907807 | U2 ID_154 | DN_chr16_80433526_2 | chr16:80422832-80438607 | --- | - | MCL1 (30) | 0.721 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h17 | TBC1D3P4 | Soto ID_469 | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (2 same-strand reads) | - |
| h17 | TBC1D3P3 | Soto ID_469 | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (4 same-strand reads) | - |
| h17 | TBC1D3P5 | - | DN_chr17_28369204_2 | chr17:28359452-28373096 | --- | - | MCL62 (2): TBC1D3P5,LOC100420311,TBC1D29P | 1.0 | wrong > wrong > wrong/wrong | (d) alignments to the family exist, every one rejected: low_shared | acc 820 bp; partner core 0 / acc 0 |
| h17 | LOC100420311 | - | DN_chr17_31496546_13 | chr17:31496547-31513862 | PLR | TBC1D29P:dominant L59/C0/P11 [struct -/-] | MCL62 (2): TBC1D3P5,LOC100420311,TBC1D29P | 0.266 | wrong > wrong > right/right | (a) fused locus; fused exon model fails the shared-exon gate | acc 2854 bp; partner core 666 / acc 0; rel MCL28 |
| h17 | TBC1D3B | Compara CF185; Soto ID_468 | DN_chr17_37113568_14 | chr17:37113569-37124531 | --- | - | MCL3 (11) | 0.799 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h17 | TBC1D3I | Compara CF185; Soto ID_468 | DN_chr17_37201555_14 | chr17:37201556-37213294 | --- | - | MCL3 (11) | 1.0 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h17 | TBC1D3G | Compara CF185; Soto ID_468 | DN_chr17_37271741_14 | chr17:37257065-37283189 | PLR | LOC101060212:mixed L53/C237/P0 [F1 zone inside the copy, not cut (share 0.6226)] | MCL3 (11) | 1.0 | right > right > right/right | (c) in family; partner bases carried: core 814 / container 0 | acc 0 bp; partner core 814 / acc 0 |
| h17 | TBC1D3H | Compara CF185; Soto ID_468 | DN_chr17_37325487_14 | chr17:37325488-37336470 | --- | - | MCL3 (11) | 0.811 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h17 | TBC1D3F | Compara CF185 | DN_chr17_37427051_14 | chr17:37427052-37438030 | --- | - | MCL3 (11) | 1.0 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h17 | TBC1D3E | Compara CF185; Soto ID_468 | DN_chr17_38911370_14 | chr17:38910592-38922356 | --- | - | MCL3 (11) | 0.942 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h17 | TBC1D3K | Compara CF185; Soto ID_468 | DN_chr17_38965100_14 | chr17:38964636-38977095 | --- | - | MCL3 (11) | 0.829 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h17 | TBC1D3D | Compara CF185; Soto ID_468 | DN_chr17_38990950_14 | chr17:38918049-39002399 | PL- | - | MCL3 (11) | 0.862 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| h17 | TBC1D3 | Compara CF185; Soto ID_468 | DN_chr17_39044717_14 | chr17:39044253-39120416 | PLR | NPEPPSP1:mixed L54/C231/P44 [F1 zone across the copy-partner gap, not cut (share 0.5510)] | MCL3 (11) | 1.0 | right > right > right/right | (c) in family; partner bases carried: core 1137 / container 515 | acc 784 bp; partner core 1137 / acc 515; rel MCL67 |
| h17 | TBC1D3P7 | - | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (2 same-strand reads) | - |
| h17 | TBC1D3P1 | Soto ID_468 | DN_chr17_60876851_14 | chr17:60876848-60891779 | --- | - | MCL3 (11) | 0.371 | right > right > right/right | - | acc 272 bp; partner core 0 / acc 0; rel MCL26 |
| h17 | TBC1D3P2 | Soto ID_468 | DN_chr17_63134600_14 | chr17:63134601-63145611 | --- | - | MCL3 (11) | 1.0 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| tes | NPIPB2 | Compara CF153; U2 ID_154; Soto ID_154 | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (2 same-strand reads) | - |
| tes | NPIPA2 | Compara CF153; U2 ID_154; Soto ID_154 | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (3 same-strand reads) | - |
| tes | NPIPA1 | Compara CF153; U2 ID_149; Soto ID_149 | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (9 same-strand reads) | - |
| tes | PKD1P6-NPIPP1 | U2 ID_149; Soto ID_149 | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (1 same-strand reads) | - |
| tes | NPIPA5 | Compara CF153; U2 ID_154; Soto ID_154 | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (2 same-strand reads) | - |
| tes | NPIPA6 | Compara CF153; U2 ID_154 | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| tes | NPIPA7 | Compara CF153; U2 ID_154; Soto ID_154 | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (3 same-strand reads) | - |
| tes | NPIPA8 | Compara CF153; U2 ID_154; Soto ID_154 | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| tes | NPIPA9 | Compara CF153; U2 ID_154 | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (3 same-strand reads) | - |
| tes | NPIPB3 | Compara CF153; U2 ID_151; Soto ID_151 | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (2 same-strand reads) | - |
| tes | LOC128966608 | U2 ID_151 | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (2 same-strand reads) | - |
| tes | NPIPB4 | Compara CF153; U2 ID_152; Soto ID_152 | DN_chr16_22381847_1 | chr16:22381848-22382376 | --- | - | MCL2 (7): NPIPB4,NPIPB5 | 0.142 | wrong > wrong > wrong/wrong | (e) no alignment record to any family member (fragment) | acc 0 bp; partner core 0 / acc 0 |
| tes | NPIPB5 | Compara CF153; U2 ID_153; Soto ID_153 | DN_chr16_22813781_1 | chr16:22813782-22814310 | --- | - | MCL2 (7): NPIPB4,NPIPB5 | 0.057 | wrong > wrong > wrong/wrong | (e) no alignment record to any family member (fragment) | acc 0 bp; partner core 0 / acc 0 |
| tes | NPIPB6 | Compara CF153; U2 ID_154; Soto ID_154 | DN_chr16_28623027_8 | chr16:28623028-28637416 | --- | - | MCL3 (5) | 0.198 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| tes | NPIPB7 | Compara CF153; U2 ID_154; Soto ID_154 | DN_chr16_28737008_7 | chr16:28737009-28751525 | --- | - | MCL3 (5) | 0.472 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| tes | NPIPB8 | Compara CF153; U2 ID_154; Soto ID_154 | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (1 same-strand reads) | - |
| tes | NPIPB9 | Compara CF153; U2 ID_154; Soto ID_154 | DN_chr16_29038827_8 | chr16:29038828-29053438 | --- | - | MCL3 (5) | 0.375 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| tes | NPIPB10P | U2 ID_154; Soto ID_154 | DN_chr16_29319686_7 | chr16:29319687-29333905 | --- | - | MCL3 (5) | 0.806 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| tes | NPIPB11 | Compara CF153; U2 ID_154; Soto ID_154 | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (5 same-strand reads) | - |
| tes | NPIPB12 | Compara CF153; U2 ID_154; Soto ID_154 | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (1 same-strand reads) | - |
| tes | LOC124907834 | U2 ID_154 | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (1 same-strand reads) | - |
| tes | NPIPB13 | Compara CF153; U2 ID_155; Soto ID_155 | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (1 same-strand reads) | - |
| tes | NPIPB14P | U2 ID_154; Soto ID_154 | DN_chr16_75785714_6 | chr16:75785715-75800090 | --R | PDXDC2P-NPIPB14P:dominant L2/C0/P0 [no-row] | MCL3 (5) | 0.337 | right > right > right/right | (c) in family; partner bases carried: core 846 / container 61 | acc 61 bp; partner core 846 / acc 61 |
| tes | NPIPB15 | Compara CF153; U2 ID_154; Soto ID_154 | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| tes | LOC124907808 | U2 ID_154 | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| tes | LOC124907807 | U2 ID_154 | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| OR-N | NPIPB10P | - | - | - | - | - | - | - | none > none > none/- | (e) >= 2 reads, none assembled (floor / gate) (3 same-strand reads) | - |
| OR-N | NPIPA1 | - | DN_NC_073242.2_15567462_1 | NC_073242.2:15567463-15571513 | --- | - | MCL1 (14) | 0.447 | right > right > right/- | - | acc 0 bp; partner core 0 / acc 0 |
| OR-N | NPIPB1P | - | - | - | - | - | - | - | none > none > none/- | (e) >= 2 reads, none assembled (floor / gate) (9 same-strand reads) | - |
| OR-N | NPIPB4 | - | DN_NC_073242.2_21074244_2 | NC_073242.2:21074245-21077017 | --- | - | MCL1 (14) | 0.212 | right > right > right/- | - | acc 0 bp; partner core 0 / acc 0 |
| OR-N | NPIPB11 | - | DN_NC_073242.2_21764017_2 | NC_073242.2:21764018-21769412 | --- | - | MCL1 (14) | 0.44 | right > right > right/- | - | acc 0 bp; partner core 0 / acc 0 |
| OR-N | NPIPB13 | - | DN_NC_073242.2_22132879_5 | NC_073242.2:22132880-22136157 | --- | - | MCL1 (14) | 0.241 | right > right > right/- | - | acc 0 bp; partner core 0 / acc 0 |
| OR-N | NPIPA8 | - | - | - | - | - | - | - | none > none > none/- | (e) < 2 primary reads on the territory (1 same-strand reads) | - |
| OR-N | NPIPA7 | - | - | - | - | - | - | - | none > none > none/- | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| OR-N | NPIPB8 | - | - | - | - | - | - | - | none > none > none/- | (e) >= 2 reads, none assembled (floor / gate) (66 same-strand reads) | - |
| OR-N | LOC124907807 | - | DN_NC_073242.2_28995202_2 | NC_073242.2:28995203-28997273 | --- | - | MCL1 (14) | 0.19 | right > right > right/- | - | acc 0 bp; partner core 0 / acc 0 |
| OR-N | NPIPA5 | - | - | - | - | - | - | - | none > none > none/- | (e) >= 2 reads, none assembled (floor / gate) (7 same-strand reads) | - |
| OR-N | NPIPA9 | - | - | - | - | - | - | - | none > none > none/- | (e) >= 2 reads, none assembled (floor / gate) (5 same-strand reads) | - |
| OR-N | NPIPA6 | - | DN_NC_073242.2_30595916_2 | NC_073242.2:30595917-30598035 | --- | - | MCL1 (14) | 0.211 | right > right > right/- | - | acc 0 bp; partner core 0 / acc 0 |
| OR-N | LOC124907834 | - | DN_NC_073242.2_31456826_2 | NC_073242.2:31456827-31462078 | --- | - | MCL1 (14) | 0.492 | right > right > right/- | - | acc 0 bp; partner core 0 / acc 0 |
| OR-N | LOC124907808 | - | - | - | - | - | - | - | none > none > none/- | (e) >= 2 reads, none assembled (floor / gate) (6 same-strand reads) | - |
| OR-N | NPIPA2 | - | DN_NC_073242.2_32426876_3 | NC_073242.2:32426877-32432044 | --- | - | MCL1 (14) | 0.456 | right > right > right/- | - | acc 0 bp; partner core 0 / acc 0 |
| OR-N | NPIPB2 | - | - | - | - | - | - | - | none > none > none/- | (e) >= 2 reads, none assembled (floor / gate) (40 same-strand reads) | - |
| OR-N | LOC128966608 | - | DN_NC_073242.2_35555594_2 | NC_073242.2:35555595-35558181 | --- | - | MCL1 (14) | 0.186 | right > right > right/- | - | acc 0 bp; partner core 0 / acc 0 |
| OR-N | NPIPB3 | - | DN_NC_073242.2_35962645_2 | NC_073242.2:35962646-35965232 | --- | - | MCL1 (14) | 0.198 | right > right > right/- | - | acc 0 bp; partner core 0 / acc 0 |
| OR-N | NPIPB7 | - | DN_NC_073242.2_99251348_15 | NC_073242.2:99251349-99267592 | PLR | LOC129527696:dominant L4/C0/P0 [no-row] | MCL6 (8): SMG1,LOC115931102,LOC129527585 | 0.014 | wrong > wrong > wrong/- | (e) no alignment record to any family member (touch-only holder) | acc 0 bp; partner core 1729 / acc 0 |
| OR-N | NPIPB6 | - | - | - | - | - | - | - | none > none > none/- | (e) >= 2 reads, none assembled (floor / gate) (4 same-strand reads) | - |
| OR-N | NPIPB14P | - | - | - | - | - | - | - | none > none > none/- | (e) >= 2 reads, none assembled (floor / gate) (27 same-strand reads) | - |
| OR-N | NPIPB12 | - | - | - | - | - | - | - | none > none > none/- | (e) >= 2 reads, none assembled (floor / gate) (17 same-strand reads) | - |
| OR-N | NPIPB5 | - | - | - | - | - | - | - | none > none > none/- | (e) >= 2 reads, none assembled (floor / gate) (19 same-strand reads) | - |
| OR-N | NPIPB15 | - | - | - | - | - | - | - | none > none > none/- | (e) >= 2 reads, none assembled (floor / gate) (40 same-strand reads) | - |
| OR-T | LOC101151653 | - | - | - | - | - | - | - | none > none > none/- | (e) < 2 primary reads on the territory (1 same-strand reads) | - |
| OR-T | LOC129533458 | - | - | - | - | - | - | - | none > none > none/- | (e) >= 2 reads, none assembled (floor / gate) (2 same-strand reads) | - |
| OR-T | LOC101144080 | - | - | - | - | - | - | - | none > none > none/- | (e) >= 2 reads, none assembled (floor / gate) (6 same-strand reads) | - |
| OR-T | LOC115933306 | - | DN_NC_073228.2_54243365_14 | NC_073228.2:54243366-54254354 | --- | - | MCL11 (7) | 0.94 | right > right > right/- | - | acc 0 bp; partner core 0 / acc 0 |
| OR-T | LOC129533792 | - | - | - | - | - | - | - | none > none > none/- | (e) >= 2 reads, none assembled (floor / gate) (6 same-strand reads) | - |
| OR-T | LOC101125558 | - | DN_NC_073228.2_54400031_14 | NC_073228.2:54400032-54411025 | --- | - | MCL11 (7) | 1.0 | right > right > right/- | - | acc 0 bp; partner core 0 / acc 0 |
| OR-T | LOC129533808 | - | DN_NC_073228.2_54880850_14 | NC_073228.2:54880851-54891837 | --- | - | MCL11 (7) | 1.0 | right > right > right/- | - | acc 0 bp; partner core 0 / acc 0 |
| OR-T | LOC129533806 | - | DN_NC_073228.2_54930301_14 | NC_073228.2:54930302-54941287 | --- | - | MCL11 (7) | 1.0 | right > right > right/- | - | acc 0 bp; partner core 0 / acc 0 |
| OR-T | LOC115934662 | - | DN_NC_073228.2_56462364_14 | NC_073228.2:56462365-56473358 | --- | - | MCL11 (7) | 1.0 | right > right > right/- | - | acc 0 bp; partner core 0 / acc 0 |
| OR-T | LOC129533797 | - | DN_NC_073228.2_56531256_14 | NC_073228.2:56531257-56618886 | PLR | LOC115933339:dominant L15/C7/P0 [no-row] | MCL11 (7) | 1.0 | right > right > right/- | (c) in family; partner bases carried: core 0 / container 885 | acc 1603 bp; partner core 0 / acc 885; rel MCL45 |
| OR-T | LOC129533813 | - | DN_NC_073228.2_56834958_14 | NC_073228.2:56834959-56845968 | --- | - | MCL11 (7) | 1.0 | right > right > right/- | - | acc 0 bp; partner core 0 / acc 0 |
| OR-T | LOC109026840 | - | - | - | - | - | - | - | none > none > none/- | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| OR-T | LOC115931404 | - | - | - | - | - | - | - | none > none > none/- | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| OR-T | LOC134759231 | - | - | - | - | - | - | - | none > none > none/- | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| KB-N | NPIPB10P | - | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (4 same-strand reads) | - |
| KB-N | NPIPA1 | - | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (9 same-strand reads) | - |
| KB-N | NPIPB1P | - | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (8 same-strand reads) | - |
| KB-N | NPIPB4 | - | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (63 same-strand reads) | - |
| KB-N | NPIPB11 | - | DN_NC_073242.2_21767099_2 | NC_073242.2:21767100-21769374 | --- | - | MCL3 (10) | 0.115 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| KB-N | NPIPB13 | - | DN_NC_073242.2_22132712_2 | NC_073242.2:22132713-22136157 | --- | - | MCL3 (10) | 0.25 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| KB-N | NPIPA8 | - | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (4 same-strand reads) | - |
| KB-N | NPIPA7 | - | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (1 same-strand reads) | - |
| KB-N | NPIPB8 | - | DN_NC_073242.2_28370348_3 | NC_073242.2:28370349-28375211 | --- | - | MCL3 (10) | 0.497 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| KB-N | LOC124907807 | - | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (10 same-strand reads) | - |
| KB-N | NPIPA5 | - | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (17 same-strand reads) | - |
| KB-N | NPIPA9 | - | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| KB-N | NPIPA6 | - | DN_NC_073242.2_30592214_1 | NC_073242.2:30592215-30598727 | --- | - | MCL3 (10) | 0.603 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| KB-N | LOC124907834 | - | DN_NC_073242.2_31453972_1 | NC_073242.2:31453973-31462773 | --- | - | MCL3 (10) | 0.786 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| KB-N | LOC124907808 | - | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (10 same-strand reads) | - |
| KB-N | NPIPA2 | - | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (138 same-strand reads) | - |
| KB-N | NPIPB2 | - | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (43 same-strand reads) | - |
| KB-N | LOC128966608 | - | DN_NC_073242.2_35555501_2 | NC_073242.2:35555502-35559153 | --- | - | MCL3 (10) | 0.264 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| KB-N | NPIPB3 | - | DN_NC_073242.2_35962552_2 | NC_073242.2:35962553-35966128 | --- | - | MCL3 (10) | 0.274 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| KB-N | NPIPB7 | - | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (67 same-strand reads) | - |
| KB-N | NPIPB6 | - | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (13 same-strand reads) | - |
| KB-N | NPIPB14P | - | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (22 same-strand reads) | - |
| KB-N | NPIPB12 | - | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (8 same-strand reads) | - |
| KB-N | NPIPB5 | - | - | - | - | - | - | - | none > none > none/none | (e) >= 2 reads, none assembled (floor / gate) (11 same-strand reads) | - |
| KB-N | NPIPB15 | - | DN_NC_073244.2_21077241_2 | NC_073244.2:21077242-21082460 | --- | - | MCL3 (10) | 0.491 | right > right > right/right | - | acc 0 bp; partner core 0 / acc 0 |
| KB-T | LOC101151653 | - | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| KB-T | LOC129533458 | - | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (1 same-strand reads) | - |
| KB-T | LOC101144080 | - | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (1 same-strand reads) | - |
| KB-T | LOC115933306 | - | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| KB-T | LOC129533792 | - | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| KB-T | LOC101125558 | - | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (1 same-strand reads) | - |
| KB-T | LOC129533808 | - | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| KB-T | LOC129533806 | - | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| KB-T | LOC115934662 | - | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| KB-T | LOC129533797 | - | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| KB-T | LOC129533813 | - | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| KB-T | LOC109026840 | - | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (1 same-strand reads) | - |
| KB-T | LOC115931404 | - | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
| KB-T | LOC134759231 | - | - | - | - | - | - | - | none > none > none/none | (e) < 2 primary reads on the territory (0 same-strand reads) | - |
