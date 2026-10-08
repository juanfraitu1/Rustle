# Pre-registration: candidate DEFINITIONS of a unit boundary inside fused transcripts (definition study; no pipeline change)

**Written 2026-09-30 (KEY=cdef) before any definition's call was scored on any substrate.** Nothing in `src/`, `tools/`,
`bench/` is edited; nothing is committed or pushed; no scoring of families is done; species are never pooled; no chimp (PTR) or
orangutan (PPY) product is opened. Scratch `/mnt/linuxdisk/tmp/rustle_figures_dev/container_defs/`; frozen instruments
`/mnt/linuxdisk/tmp/rustle_figures/container_defs_frozen/` (sha1 §6). The orchestrating session's task is the mandate.
This binds once Amendment 1 records this file's sha1 (before any DEV scoring).

## 0. The question, and what was seen before this file

Container v1 (`bench/family_container.py`, PREREG_fusion_container_sim §1) calls a block CORE iff it aligns to another member
of the SAME family, and failed for three reasons: it runs after the families stage and never changes a family; for CO-DUPLICATED
partners (PKD1P pieces of PKD1P3-NPIPA1, PKD1P1-NPIPA6, ...) every block is core and the container is empty; and it depends on the
family the locus was placed in, the very decision a fusion corrupts. F1 / F1v2 (`bridge_regroup.rs`, the assemble default) define a
bridge from READS (a junction that splits the other transcripts into an upstream and a downstream group, with a PAS-proven
unprimed 3' end upstream inside the bridged intron and an own-promoter start downstream) and keep it only when the link carries
fewer reads than each side (F1v2), so DOMINANT links (PKD1P3-NPIPA1, PDXDC2P-NPIPB14P) are never split, and a side with no
standalone reads has no read proof at all (r1134). Node cut at the annotated boundary (r846), truth-side chimera policy (r845),
read-level chimeric-bridge split (r967) and a cover for F1v2 bridges (r1184-r1188) are closed and are not re-proposed.

**Question.** Which DEFINITION of a unit boundary (a junction of a transcript that separates two units, each of which may be a
member of a different family) is clean (no new free constant beyond the shipped edge rule, combinatorial), can be computed at
the families stage without the family it corrupts, and finds boundaries where container v1 and F1 / F1v2 cannot: dominant links,
co-duplicated partners, and sides with no standalone reads? Candidates: V1 (container v1 as a transcript call), R (F1's bridge
junction), R2 (F1v2's), W (an alignment-witness rule), and their combinations.

**Seen before this file (plainly).** Read: PREREG_fusion_container_sim (incl. Outcome and human descriptive), PREREG_f1_bridge_locus,
PREREG_f1v2_readshare (incl. verification), PREREG_o1_cover_growth (incl. Outcome); register rows 845, 846, 967, 1013-1018, 1134,
1145-1151, 1184-1188; `family_container.py`; `bridge_regroup.rs` (lines 1-420) and `f1_bridge.py`; the seven memory files named in the
task. The parallel sibling pre-registration `PREREG_container_units_mechanism_2026-09-30.md` (a known-junction execution test) was
found and its §0-§5 skimmed; none of its results exist and none is used. Truth-side only (annotation x default GTF x F1 tables; no
definition's calls): the gates of §6 (G-T: the four per-transcript FUSED readings of r1147 reproduced exactly), a census of
positives, windows and strata on every substrate INCLUDING the held-out ones (§5 states the counts), four example transcripts with
multi-intron windows (testis), six with no window (OR6737 NC_073224.2, a non-dev contig), and the list of chr16 (dev) truth windows
whose parents are NPIP / PKD1P / SMG1P / PDXDC genes. **Seen on the PAF side only (no call):** the testis and chr16 PAF scans
(counts of records and of exon-exon records) and the frozen container's summary on testis and chr16. One plumbing spot check of W on
the named NPIP readthrough fusions of chr16 printed nothing (no exact-window positive has both parents among the names) and no W
output was seen. **Not seen:** any V1, R, R2 or W call of any product, any score, any NULL.

## 1. The definitions (binding)

### 1.1 Objects

- **Transcript:** a spliced (>= 2 exons) transcript of the sample's DEFAULT (BASE) GTF, the 2026-09-25 default-pipeline product
  (A119b `631c9f11`, testis `50f239d6`, OR6737 `3a8c6410`, KB3781 `f1022f4b`): the GTF BEFORE F1v2 regrouping, where fused loci exist.
- **Locus:** a `gene_id` of that GTF, keyed `CONTIG:START-END` over the exons of ALL its transcripts (the families stage's key; two
  gene_ids with one key are one locus). Its **exon bases** = the union of all its transcripts' exons (merged where they share
  >= 1 base).
- **Call:** `(transcript, intron j)`, j = 0.. in genomic order (the intron between exon j and j+1). A definition returns, per
  transcript, the set of introns it calls a boundary. No definition reads the annotation.
- **Product:** one families-stage run on the BASE GTF (frozen `mcl_families`, driver flags `--min-exonic-bp 1 --min-shared-exon-frac 0.60
  --emit-units`; the 09-29 default `--min-cov-shorter 0.70` changes edge WEIGHTS only, admission is identical, r1019): genome-wide for
  testis / OR6737 / KB3781; one chromosome at a time for A119b (as the O1-cover test; cross-chromosome pairs are absent). Its
  admitted-edge graph = `emu.py` (fac9a560) R0, gated byte for byte against the product's `clusters.tsv` (G-E).

### 1.2 V1 (container v1 as a transcript call)

```
V1 (inputs: the product's clusters.tsv, loci.tsv, loci.gff3, loci.paf and the GTF; instrument: the frozen family_container.py e197ccb3, unchanged)
  Every exon of a transcript of a CLUSTERED locus (or of a record that loci.tsv folds into one) lies in one exon block of that locus;
  the block is CORE (an aligned column joins an exon base of the block to an exon base of another member of the SAME family) or
  ACCESSORY; an accessory block is RELATED to family Y when the same column test joins it to an exon base of a member of Y != own.
  A transcript is called iff its exons read, in genomic order, as ONE contiguous run of core exons and ONE contiguous run of accessory
  exons (exactly one core/accessory transition) and at least one accessory exon's block is related to another family.
  The call is the intron at the transition.
```
Constants: none beyond the shipped edge rule and MCL that define the family. Needs at run time: the families clusters (hence a family
decision) and the PAF with CIGARs.

### 1.3 R and R2 (F1 and F1v2 bridge junctions)

```
R   (inputs: F1's junction table, BAM-derived: f1_bridge.py 37ee8e77 --mode full, run on the BASE GTF; the shipped Rust port writes the same
     table byte for byte on the held-out inputs, bench/ASSEMBLY_POLISH.md addendum 3)
  Every transcript of T_J for every junction J with BRIDGE(J) = STRUCTURAL(J) and UP-proof(J) and DOWN-proof(J) (the RULE block of
  PREREG_f1_bridge_locus §1.1), WITHOUT any read-share condition. The call is J.
R2  (f1v2.py b4e788ad --rule min)
  R restricted to the junctions with MINORITY(J): reads(T_J) < reads(UP_J) and reads(T_J) < reads(DOWN_J). These are exactly the bridges of the
  shipped default (`--bridge-regroup f1v2`).
```
Constants: F1's (21 bp 3' cluster gap, >= 2 reads, PAS hexamers AATAAA / ATTAAA in oriented [mode-35, mode-10], >= 60% A in 20 bp or A6
priming, start-cluster gap 100 bp and >= 3 reads, V1 >= 1 -- PREREG_f1_bridge_locus §1.2); R2 adds only the strict ordering
(1/2 is the definition of "minority"). Needs: the BAM (proofs) and the GTF; no PAF, no families.

### 1.4 W (alignment-witness rule)

```
W   (inputs: the GTF, the product's loci.paf with CIGARs, the product's admitted-edge graph)
  WITNESS SET of an exon e of a transcript of locus m: the loci m' != m, with spans that do NOT overlap m's span (same contig and
  overlapping spans = the same DNA, not another copy), such that some aligned CIGAR column (M, =, X) of a PAF record between m and m'
  joins a base of e to an exon base of m' (the container's column test).
  For intron j of the transcript: WL = the union of the witness sets of exons 0..j, WR = the union over exons j+1..;
    LO = WL \ WR ("occurs without the other side"), RO = WR \ WL, S = WL ∩ WR (the loci that SPAN the boundary: co-duplicated fused loci).
  Intron j is a SUPPORTED BOUNDARY iff
    (a) LO != {} and RO != {} ("A occurs without B and B occurs without A"), and
    (b) in the edge graph G with the nodes {m} ∪ S deleted, no connected component contains both a member of LO and a member of RO
        (every LO node and every RO node lie in different components; an LO or RO node with no edge is its own component).
  The call is j. No constant: a witness needs one column; the graph is an existing object.
```
Why this form. (i) The witness test is the container's column test, so W and V1 read the same evidence, but W reads it per EXON of the
transcript, not per block of a locus already placed in a family: the boundary is judged before, and independently of, the family decision
(failure iii). (ii) Loci that SPAN the boundary are removed before connectivity, because a co-duplicated fused locus aligns to both
halves and would otherwise connect the two families (failure ii); a real single gene whose paralogs carry both halves has LO or RO empty
and is never called. (iii) m is removed because every witness aligns to m by construction. (iv) Overlapping loci are not witnesses
because an alignment of two loci that share genome bases is the shared stretch aligned to itself; this is a property of the
witness, not a number. (v) "no component contains both" (every-pair reading) is registered over "some pair in different components"
(which would fire whenever one isolated witness exists). (vi) The transcript's own exons, not the locus's union, define the sides,
so the X-only and Y-only transcripts of a fused locus are not called and the X-Y transcript is.

### 1.5 W variants and the DEV rule that fixes the registered W\*

| name | witnesses | graph G of (b) |
|---|---|---|
| **WA** (registered default) | every non-overlapping locus with an exon-exon column | the shipped ADMITTED edges (the graph the families are made from) |
| WN | only loci with an admitted edge to m | admitted edges |
| WC | as WA | the witness relation itself (pairs with an exon-exon column) |
| WP | as WA | every PAF record (any identity, length, exon content): "all PAF hits" |

All four are computed on every substrate. **W\*** (used by the combinations and by the headline statements) is fixed on DEV ONLY:
F_w(V, s) = 2 TP / (2 TP + FP + FN) in the WINDOW reading (§2), all strata, on the DEV substrate of each dev species s (A119b chr16 +
chr17; OR6737 NC_073244.2); W\* = WA unless another variant V has F_w(V, s) >= F_w(WA, s) on BOTH dev species and > on at least one;
among several, the largest sum, ties in the order WN, WC, WP (`select_w.py`). After Amendment 2 W\* never changes.

### 1.6 Combinations (rule frozen)

Set operations on the `(transcript, intron)` call sets, with W = W\*: **R2 ∪ W** (headline combination), **R ∪ W**, **R2 ∩ W**,
**V1 ∪ R2 ∪ W**. No weights, no voting, no extra condition.

### 1.7 Every constant, by definition

| definition | constants | source |
|---|---|---|
| V1 | none of its own | the shipped edge rule (identity >= 0.7, cov_longer >= 0.3, >= 300 bp, >= 1 exonic base, shared exon fraction >= 0.60), MCL I 2.8 prune 1e-9, size >= 2, fold within clusters, all of which define the family |
| R | F1's (see §1.3) | PREREG_f1_bridge_locus §1.2, all inherited from the readthrough filter and `--polish-tes` |
| R2 | R's; the strict ordering | PREREG_f1v2_readshare §1 |
| WA, WN | none | the shipped edge rule defines the graph |
| WC | none | |
| WP | none | |

The only numbers of this file that are not inherited are in the JUDGING of §7 (the power floor n >= 10), never in a definition.

## 2. Truth (annotation only)

Annotation: human = `/mnt/linuxdisk/tmp/regress/chm13.gff` (the CHM13 v2.0 RefSeq full GFF; never `HSA_genomic.gff`); gorilla =
`winloci_data/GGO_genomic.gff`. Gene = a `gene` / `pseudogene` record not described as a readthrough; its exon bases = the merged exon rows
of its transcripts (CDS rows if none; the span if none) -- the gene set and exon model of `fused_pt.py` (c42a4605), gated (G-T). For a
spliced transcript t (exons on t's strand; RefSeq readthrough-described records are NOT in the lookup):
- G_j = the genes whose exon bases share >= 1 base with exon j; H = the union (the hit genes).
- A pair (a, b) of hit genes QUALIFIES iff their spans are disjoint, or they are the two parents named by a RefSeq readthrough record
  (a record named `A-B` whose parents resolve among the annotated genes of the contig and strand: human 178 of 209 resolve).
- **P (positive):** some pair qualifies. **N (negative):** exactly one hit gene (t lies inside one annotated gene). **U0:** no hit gene.
  **U2:** >= 2 hit genes, no qualifying pair.
- A qualifying pair (a, b) with every exon over a before every exon over b gives the **WINDOW** [last exon over a, first exon over b - 1]
  of introns. A window of ONE intron is **EXACT** (the intron between the last exon over gene 1 and the first over gene 2); a wider
  window (exons over no gene, or over a third gene, sit between) is **WIDE**. A wide window that contains an exact one is dropped.
  A positive with no window (an exon over both genes, interleaved exons) is **PX**: counted, never judged.

**Two readings, both reported everywhere.** STRICT (the task's "exact junction"): units = exact windows; a call inside a wide window is
AMBIGUOUS (neither TP nor FP). WINDOW (primary, because the motivating cases are wide windows: on chr16 the 19 PDXDC2P-NPIPB14P and the
PKD1P3-NPIPA1 transcripts have wide windows): units = all windows, a call anywhere inside a window is a TP. For every call:
**TP** (inside a counted window), **AMB** (inside a wide window, STRICT only), **FP_wrong** (a positive, outside all its windows),
**FP_N** (a negative: a single-gene false split), **UNJ** (U0, U2, PX: counted, never in a precision).
Readthrough-described records are positives "too" through their named parents (the qualifying rule above); a record whose parents do not
resolve (31 human) contributes nothing.

## 3. Strata of the positives (windows; nothing here enters a definition)

1. **Multi-copy parents.** A gene is multi-copy iff in a Compara Primates family with >= 2 genes (human), or Liftoff multi-copy (>= 1 extra
   copy at sequence_ID >= 0.95, both exon unions >= 200 bp; `figures/_liftoff.copy_pairs`), or in a RefSeq NAME FAMILY (>= 2 annotated genes
   with the same normalised description: lower case, trailing "-like" removed, placeholders "uncharacterized", "LOC<n>", "long intergenic"
   excluded). `1a` = the window has >= 1 multi-copy parent (left or right gene set); `1b` both sides; `1s` strict (>= 1 parent in Compara or
   Liftoff); `1n` no parent multi-copy. The family-relevant precision uses the calls on transcripts that touch a multi-copy gene.
2. **Dominant link vs bridge (F1v2's share rule).** For the window's introns, F1's decomposition of the locus's transcripts (recomputed from
   the GTF; gated on the bridge rows against F1v2's table): `2_bridge` = STRUCTURAL and MINORITY; `2_dominant` = STRUCTURAL and the link carries
   at least as many reads as a side; `2_nonstruct` = not STRUCTURAL (no standalone upstream or downstream component in the gene_id).
   (A bridge that also passes F1's read proofs is an R2 bridge.)
3. **Co-duplicated partner.** Under the WA witness relation, S of the window's first intron: `3_codup_fused` = S contains a locus that holds a
   positive transcript (another fused locus spans the boundary); `3_codup_any` = S non-empty; `3_no_span` = S empty.
4. **Standalone support.** `4_standalone_both` = at least one GTF transcript (any locus, same strand, other than t) overlaps the left parents'
   exon bases without the right parents', and one the right's without the left's; `4_standalone_not` otherwise.
5. `rt_named` = the parents are the two named parents of a RefSeq readthrough record (human only).

## 4. Scoring, matched NULL, single-gene cost

Per definition, per substrate, per stratum: truth units (windows), TP units (windows with >= 1 call inside), FN = units - TP,
recall = TP / units; FP = FP_wrong + FP_N (calls); precision = TP_calls / (TP_calls + FP) where TP_calls counts calls inside counted windows
(STRICT: exact windows only; WINDOW: exact and wide); **family-relevant precision** (calls on transcripts touching a multi-copy gene);
FP_N split by whether the gene is multi-copy; **single genes cut** = the distinct single genes with >= 1 FP_N call, and the calls on the
twelve recurrent F1 cases (CALD1, COL12A1, ARHGEF3, NCAPD3, BAZ2B, DLC1, LARP1B, MCM9, MIA3, PGBD1, PLEKHA5, RABGAP1L); loci (gene_ids)
with a TP and with an FP; UNJ counts by class. **Positives found by one definition and missed by the others:** for every window, which of
{V1, R, R2, W\*} finds it (patterns and unique counts). **Matched NULL:** for each definition, each call is replaced by ONE uniformly random
intron of the same transcript (5 draws, seeds `cdnull:<sample>:<substrate>:<definition>:<k>`, duplicates collapse) and scored identically.
The unit of a call is the transcript (the per-transcript reading of r1147); a transcript with several isoforms counts once per isoform, so
loci counts are reported beside.

## 5. Substrates and exposure

| substrate | contigs | role | products |
|---|---|---|---|
| **A119b DEV** | chr16, chr17 | dev (human) | chr16: the F1v2 dev GTF `hsa16.BASE.gtf` (same transcripts and gene_ids as the default GTF's chr16; only TPM differs), BASE.fam of `o1_cover/dev/hsa16`; chr17: `o1_cover/held/human_A119b/chr/chr17` |
| **OR6737 DEV** | NC_073244.2 | dev (gorilla) | `rt_arms/gorilla_OR6737` BASE.fam (genome-wide), scored on the contig |
| **A119b HELD** | chr1-12, 14, 15, 19, X, Y, M | held-out (human) | `o1_cover/held/human_A119b/chr/*` (per chromosome) |
| **testis HELD** | every contig | held-out (human, other library) | `rt_arms/human_testis` BASE.fam (genome-wide) |
| **KB3781 HELD** | every contig | held-out (gorilla, other individual) | `rt_arms/gorilla_KB3781` BASE.fam (genome-wide) |
| OR6737 REST | every contig but NC_073244.2 | DESCRIPTIVE (F1's held-out contigs; never in a statement) | as OR6737 DEV |

Not scored: A119b chr13 (its BASE families run is unfinished), chr18, chr20, chr21, chr22 and unplaced contigs (no BASE families product;
chr20-22 and chr18 were outside the O1-cover substrate). R and R2 could be scored on them from F1's tables; that is not done here.

**Exposure (plainly).** The held-out substrates were NOT fresh for the read rule: A119b (minus chr16 / 20 / 21 / 22) and testis were the
held-out substrates of the F1v2 test (verdict-bearing; its fifth use of each), the gorilla samples the held-out substrates of the F1 test
and the DEV substrates of the F1v2 design; chr17 was a held-out chromosome of both the F1v2 and the O1-cover tests and is DEV here by the
task; chr16 and NC_073244.2 were dev of F1 / F1v2. **R and R2 therefore carry no fresh evidence here; V1 and W are new** (V1 was run only on
the simulation and described on chr16). Truth-side counts seen before this file (exact windows / wide windows / PX transcripts):
A119b DEV 368 / 147 / 116, HELD 2,754 / 1,145 / 502; testis 303 / 37 / 30; OR6737 NC_073244.2 65 / 6 / 21, rest 581 / 164 / 306; KB3781
653 / 245 / 402; window-strata counts were seen for the exact windows (multi-copy parents: A119b HELD 382, DEV 93, testis 47, OR dev 9,
KB 212; F1 classes BRIDGE / DOMINANT / NONSTRUCT: A119b HELD 309 / 149 / 2,296, testis 18 / 8 / 277, KB 84 / 15 / 554, OR dev 14 / 0 / 51).

## 6. Instruments, inputs and gates

Frozen, `/mnt/linuxdisk/tmp/rustle_figures/container_defs_frozen/SHA1SUMS` (sha1 `12af2d41`): `cd_lib.py` 6bfd1b6e, `cd_defs.py` 5acdb674,
`cd_scan.py` 4ffaeadc, `score.py` 9d37f045, `run_score.py` 707b3157, `run_defs.py` ccb97deb, `truth_table.py` 3f62f199, `products.py`
4ead1f0e, `run_emu.py` 35acdbfa, `select_w.py` 192e4129, `test_cd.py` e5e3b83a (23 tests pass), `census.py` a23fa6cc, `gate_truth.py`
ae30deb4, `build_ann.py` 0f70529f. Reused frozen: `f1_bridge.py` 37ee8e77 (parser, components, side), `family_container.py` e197ccb3,
`emu.py` fac9a560, `f1v2.py` b4e788ad and `fused_pt.py` c42a4605 (G-T), `figures/_liftoff.py` 952fc546 (read-only). Interpreter
`/home/juanfra/miniforge3/bin/python3` 3.13.12 + pysam 0.23.3. Inputs: GTFs as §1.1 (+ `hsa16.BASE.gtf` 46a24de7); F1 tables A119b `4b42546d`
/ `2104f07f` (junctions / F1v2 bridges), testis `2f5d8307` / `a4a361d7`, OR6737 `dcb25438` / `95a2579c`, KB3781 `09a8de18` / `8bba2b49`; Compara
`8e8affdd`; Liftoff loci human `ba639599`, gorilla `d48bb51a`; PAF and graph sha1s in
`container_defs/tmp/input_sha1.txt` (testis PAF `9b0e1e68`, graph `68f28590`; OR6737 `6ead9746` / `b33a4685`; KB3781 `45b2e683` / `1a1a85cf`;
chr16 `bbedf886` / `fcdabbb3`; chr17 `d9a07dda` / `c0006d75`).

| gate | what must hold | status |
|---|---|---|
| G0 | `sha1sum -c SHA1SUMS` and the reused frozen files match, before every run | each run |
| G-T | `cd_lib`'s annotation records equal `fused_pt`'s (human) and its per-transcript FUSED count equals the frozen reading on the same transcripts (the published BASE numbers 1,754 / 194 / 502 / 567) | **passed before the freeze** (A119b, testis, OR6737, KB3781: exact) |
| G-E | `emu.py` R0 reproduces the product's `clusters.tsv` byte for byte (so its graph is the shipped graph) | **passed** for chr16, chr17..chrM, testis, OR6737, KB3781 |
| G-K | the PAF scan finds every PAF name among the GTF loci (`unknown_locus` = 0) and the loci.gff3 keys equal the GTF keys | each product |
| G-X | `cd_scan.exon_runs` equals the frozen container's `exon_columns` base for base on 400 random CIGARs, both strands | unit test, passed |
| G-R | the F1 / F1v2 tables of a sample have the bridge-junction counts of the frozen stats (A119b 887 / 425, testis 22 / 14, OR6737 97 / 64, KB3781 92 / 64) | each sample |
| G-F | for bridge rows, F1's decomposition recomputed by `cd_lib.f1_sides` equals F1v2's table (n_up, n_down, reads) | each sample |

A failed gate is a bug: it stops that product / substrate and is fixed in the instrument only (recorded as an amendment), never in a
definition, a stratum, a null or a criterion.

## 7. Judging (descriptive; no default is decided here)

A statement about a definition D and a stratum X on a held-out substrate s is made only when X has >= 10 window units on s (a power guard
of the judging, not of any definition). D **detects** X on s iff TP(D, X) exceeds the maximum over the five NULL draws AND D's precision
exceeds the maximum NULL precision (both readings reported). D **covers** X for a species iff it detects X on every judged held-out
substrate of the species (human: A119b HELD and testis HELD separately; gorilla: KB3781). The **false-split cost** of D is reported
beside it: FP_N calls, single genes cut, FP_N per 1,000 negative transcripts, the same inside multi-copy genes, the twelve named genes.
**Run-time needs** are stated per definition (reads, PAF, family clusters) with the pipeline stage it could occupy. No clause of this file
is a verdict on adopting a definition.

## 8. Predictions (this author's probabilities, before any DEV or held-out score)

1. All gates pass (0.90); G-R and G-F pass on every sample (0.85).
2. **V1** has window recall < 0.05 on every held-out substrate (0.85) and recall 0 in `3_codup_fused` (0.9: a co-duplicated partner is core).
3. **R2** has window recall 0 in `2_dominant` and `2_nonstruct` by construction (certain), recall in `2_bridge` >= 0.6 (0.7), and precision > R's on
   every held-out substrate (0.9).
4. **W\*** has window recall above R2's in `1a` on every held-out substrate (0.80) and precision below R2's on every held-out substrate (0.80).
5. W\* has recall > 0 in `2_dominant` and in `4_standalone_not` on A119b HELD (0.80).
6. W\* **covers** `3_codup_fused` on A119b HELD when it has >= 10 units (0.45).
7. The single-gene cost of W\* exceeds R2's (single genes cut) on every held-out substrate (0.75); W\* cuts >= 3 of the twelve named genes on KB3781 (0.40).
8. Matched NULL: the maximum NULL precision is below D's for D = R2 (0.9), R (0.8), V1 (0.8), W\* (0.60).
9. No definition or combination reaches precision >= 0.9 and recall >= 0.5 (window reading) on any held-out substrate (0.90).
10. On A119b HELD, >= half of W\*'s TP windows are found by neither R2 nor V1 (0.80); R2 ∪ W\* has higher recall than either and lower precision than R2 (0.95).
11. DEV rule: WA retained (0.60); WN (0.25), WC (0.10), WP (0.05).

## 9. Falsifiers of the design reasoning (reported whatever happens)

- **Z1 "Removing the spanning loci separates the two halves of a co-duplicated fusion."** Falsified on a substrate where, among windows with (a)
  true and `3_codup_fused`, (b) holds for fewer than half.
- **Z2 "A single gene's paralogs carry both halves and are removed as spanning, so W does not cut single genes inside multi-copy families."**
  Falsified on a substrate where W\*'s family-relevant window precision is below 0.5.
- **Z3 "One aligned column is a specific witness."** Falsified where W\*'s precision does not exceed its maximum NULL precision.
- **Z4 "The boundary can be judged without the family decision."** Reported: the share of W\*'s TP windows whose transcript sits in an unclustered
  locus (V1 cannot see them).

## 10. Order, stop rules, machine rules

1. Amendment 1 (this file's sha1, before any DEV scoring). 2. DEV products: per product `cd_scan` (+ `--merge`), the frozen container, `run_defs`;
`truth_table` per sample (done); `run_score A119b:DEV`, `run_score OR6737:DEV`; `select_w.py` -> W\*. 3. Amendment 2: DEV numbers (brief), W\*,
before any held-out product is scored. 4. Held-out products and scores in the order A119b HELD, testis, KB3781, then OR6737 REST.
5. Outcome. After Amendment 2 nothing changes in a definition, a stratum, a null, a criterion, a truth or a substrate. A failed gate stops its
product; a step that hits the time cap is re-run once.
Heavy steps (emulator and PAF scans of genome-wide products) via `bash tools/rlock.sh heavy`, light steps via `light`; foreground;
`TMPDIR` under `/mnt/linuxdisk`; never `pkill -f`; one heavy job at a time (another session shares the lock; waits, never bypassed).

## 11. Not in this test

Families scoring; any change to `src/`, the assembler, the families stage or a default; the sibling's known-junction execution; W / V1 on the F1v2
default's own families products (BASE is where fused loci exist); cross-chromosome witnesses on A119b; fusions between an annotated and an
unannotated gene (the truth cannot see them: they count as FP_N when the transcript hits one annotated gene, and UNJ when it hits none);
chimp and orangutan.

## 12. Hostile self-review (fixes applied above)

1. **"The exact-junction reading drops 29% of the windows, including the PDXDC2P-NPIPB14P fusion."** Yes; hence two readings, WINDOW primary.
2. **"A call inside a wide window is a cheap TP, and the null inherits it."** The NULL draws a random intron of the same transcript under the same
   truth; a window of k introns of n makes the NULL's TP probability k / n, which is the baseline a definition must beat.
3. **"The truth cannot see unannotated parents."** True; such calls are FP_N (one annotated gene) or UNJ (none): an underestimate of precision for
   every definition alike, not a ranking change.
4. **"W's witness relation (one column) is permissive."** It is the container's test and the user's; WN restricts to admitted neighbours, WC / WP vary
   the graph, and Z3 / the NULL ask whether it is specific.
5. **"Excluding overlapping loci is a choice."** It is a property of the witness (same DNA); the effect of admitting them is not measured here.
6. **"A119b per chromosome hides cross-chromosome copies."** Declared; testis and gorilla are genome-wide; A119b's numbers for V1 and W are lower bounds.
7. **"Four W variants are selection freedom."** One is registered (WA); the DEV rule is conservative (dominance on both dev species under F_w, which
   penalises calling nothing) and the gorilla dev has 65 exact windows, so a switch is unlikely; all four are reported on the held-out substrates.
8. **"The F1 class strata use the GTF, not the BAM."** The decomposition is F1's own; proofs (BAM) enter only R / R2, and the bridge rows are gated.
9. **"Standalone support is measured by assembled transcripts, not reads."** The assembler's floor is 2 reads; `reads` are summed and reported.
10. **"Multi-copy by RefSeq description is broad (33% of gorilla genes)."** Stratum `1s` (Compara / Liftoff only) is reported beside `1a`.
11. **"chr17 is a held-out contig of earlier tests."** Disclosed in §5; it is DEV here by the task and scores nothing in a held-out statement.
12. **"Metric traps (feedback_metric_traps)."** The universe is truth-side and fixed per substrate (annotation x GTF; no denominator is conditioned on a
    prediction); UNJ is reported, never silently dropped; the NULL is matched on transcripts and call counts; species are never pooled.

## Amendments

### Amendment 1 -- the freeze (2026-09-30 10:36, written BEFORE any DEV scoring and before `run_defs` was run on any product)

**This file's sha1 before this amendment:** `0ca94ec2922374e281d050dde4e38d0fdb8eb834` (29,235 bytes; a byte copy is kept at
`/mnt/linuxdisk/tmp/rustle_figures/container_defs_frozen/PREREG_container_units_definition_2026-09-30.pre_amendment1.md`, born 10:35:27).
The frozen text equals this file with the Amendment and Outcome bodies removed. Acceptance: the orchestrating session's task is the mandate.
Instruments as §6 (`SHA1SUMS` `12af2d41`); every run is preceded by `sha1sum -c` (G0). Run so far, all truth-side or PAF-side, none scored: the
annotation caches, `gate_truth.py` (G-T), `truth_table.py` on the four samples, the frozen container on testis and chr16, `cd_scan.py` on testis
and chr16, `run_emu.py` (G-E) on OR6737, KB3781, chr16. Products still to be made for DEV: chr17 scan + container, OR6737 scan + container
(+ the OR6737 and KB3781 scans for the held-out step). `run_defs.py` and `run_score.py` have NOT been run on any real product.

### Amendment 2 -- DEV, one instrument fix, and the registered W (2026-09-30 10:53, written BEFORE any held-out product is scanned, run or scored)

(File sha1 at 10:53:12, the last edit before the first held-out command (the chr1 scan, 10:53:47): `13e2a15b6c7b44ad2bd42700d29b4e94297d8246`; a first draft of this heading carried the wrong time 11:40, corrected afterwards.)

**DEV products run** (frozen instruments, G0 before each; G-K: `unknown_locus` = 0 in every scan; G-R: F1 / F1v2 bridge junctions 887 / 425 on A119b,
97 / 64 on OR6737): A119b chr16 (hsa16 dev GTF) and chr17, OR6737 (genome-wide product, scored on NC_073244.2). Truth tables rebuilt by the
frozen `truth_table.py` (content sha1 of the tx tables: A119b `8f0d84af`, testis `ec7790c4`, OR6737 `bfa9048d`, KB3781 `271fb641`).
**Instrument fix (one, after the first DEV scoring; definitions, strata, nulls and criteria unchanged).** `score.strata_of` assigned no stratum 3
label to a window whose locus has NO witness at all (`run_defs` stores diagnostics only for loci with PAF hits, and S is empty for them by
definition): those windows are now `3_no_span`. Effect on the A119b DEV units: `3_codup_fused` / `3_codup_any` / `3_no_span` = 169 / 217 / 298
(before the fix 81 was shown as `3_no_span`). `score.py` 9d37f045 -> d4731117, `test_cd.py` e5e3b83a -> 0b63978a (24 tests pass), `SHA1SUMS`
12af2d41 -> `ea1555be`; no other file changed. DEV was re-scored after the fix; every number below is the re-scored one.

**DEV results (window reading; units = windows: A119b 515 = 368 exact + 147 wide; OR6737 71 = 65 + 6).**

| A119b chr16+chr17 | calls | TP+AMB | FP (FP_N) | precision (strict / window) | recall (strict / window) | NULL max precision (window) | single genes cut |
|---|---|---|---|---|---|---|---|
| V1 | 151 | 54 | 84 (69) | .27 / .39 | .084 / .105 | .087 | 25 |
| R | 518 | 51 | 461 (448) | .087 / .100 | .120 / .099 | .016 | 34 |
| R2 | 90 | 32 | 57 (57) | .305 / .360 | .068 / .062 | .067 | 7 |
| **WA** | 7,392 | 237 | 6,212 (5,113) | .018 / .037 | .318 / .359 | .028 | 158 |
| WN / WC / WP | 109 / 3,571 / 1,012 | 0 / 95 / 20 | 83 / 3,103 / 979 | 0 / .011 / .003 (strict) | 0 / .095 / .008 | .043 / .024 / .014 | 3 / 103 / 26 |
| R2 ∪ WA | 7,475 | 267 | 6,264 | .022 / .041 | .386 / .417 | .028 | 163 |

| OR6737 NC_073244.2 | calls | TP+AMB | FP | precision (window) | recall (window) | NULL max precision | single genes cut |
|---|---|---|---|---|---|---|---|
| V1 | 59 | 0 | 59 | 0 | 0 | .017 | 12 |
| R = R2 | 10 | 9 | 1 | .90 | .127 | .30 | 1 |
| **WA** | 438 | 17 | 288 | .056 | .099 | .060 | 39 |
| WN / WC / WP | 0 / 92 / 22 | 0 / 5 / 0 | 0 / 58 / 22 | - / .079 / 0 | 0 / .042 / 0 | - / .077 / 0 | 0 / 14 / 5 |

**The §1.5 rule:** F_w(WA, WN, WC, WP) = .0535 / 0 / .039 / .0199 (A119b) and .0383 / 0 / .0455 / 0 (OR6737); WC beats WA on OR6737 and loses on A119b, so no
variant dominates: **W\* = WA** (`select_w.json`). **Read plainly before the held-out step:** on DEV the registered W floods (precision .037, at chance:
the matched NULL reaches .028), cutting 158 single genes on A119b; R2 is the only definition whose precision clears its NULL by a wide margin.

**Post-hoc DEV diagnostics (exploratory; nothing of it is a candidate; the registered set is unchanged).** Among WA's A119b DEV calls, 69% are on
transcripts inside one annotated gene (FP_N), 15% are wrong junctions of positives, 13% are UNJ; FP_N calls sit in loci with many spanning loci (median |S| 6 vs 2
for the TP windows) and 47% of them touch a multi-copy gene. Constant-free features tried on the calls (LO witness at exon j and RO witness at exon
j+1; exactly one supported intron per transcript; S empty or not) raise the precision at most to .17 at a small fraction of the recall, so no W2 is registered.
No held-out product has been read.

### Amendment 3 -- execution notes (2026-09-30, after the held-out scoring; no definition, stratum, null, criterion, truth or substrate changed)

1. **Gate scripts added after the scoring (frozen, `SHA1SUMS` `7cf79b24`):** `gate_f1.py` c7f13924 (G-R, G-F), `gate_keys.py` f8a11db2 (G-K), `report.py`
   b6fce667 (formatting only). The first `gate_f1.py` compared F1v2's `n_up` (transcripts in UP components) with `cd_lib`'s `n_up` (components) and
   reported false mismatches; that comparison was removed (the gate compares n_TJ, reads_TJ, reads_up, reads_down and keep); every row then agreed.
2. **Post hoc scripts** (after the held-out scoring, descriptive only, archived in `container_defs_frozen/posthoc/`): `post_tx.py` (transcript-level
   view), `top_fp.py` (genes carrying FP_N), `named16.py` (chr16 dev case study), `falsifiers.py` (Z1, Z4), `diag_dev*.py` (DEV features of WA calls), `example_nbpf1.py` (one NBPF1 transcript), `lead_dev.py` (the DEV-only lead: R in W-flagged transcripts).
3. **Timing.** Amendment 2's heading first carried 11:40; the real time was 10:53 and the order (file edit 10:53:12, first held-out artifact 10:53:47) is in its text.
4. **Provenance of the products (BASE arm only).** All-vs-all `minimap2 2.30-r1287 -x asm20 -c -X -N 50 -p 0.1 --secondary=yes` (sharded by `tools/mm2_shard.sh`), `mcl_families --from-gtf <default GTF>
   --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units`: A119b chromosomes by the O1-cover run of 09-29 (`fj_bin_frozen` 91ef2e1c; chr13 unfinished, never re-run here), testis / OR6737 / KB3781 by the 09-25 / 09-27
   `rt_arms` runs (rt-era frozen binary). The ALL / CORE arms of the O1-cover test (F1v2 GTFs) were not used. G-E shows that one emulator reproduces every product.
5. **Machine.** Every run of mine used `tools/rlock.sh` (light, heavy for the genome-wide scans of OR6737 / KB3781 and the emulator); no background job of mine; an independent
   verification agent ran light jobs of its own (a few short exploratory scripts of its own, < 15 s and < 2 GB each, bypassed `rlock.sh light`; disclosed in its report and in the Outcome).

## Outcome (2026-09-30)

**Gates.** G0 (every run), G-T (exact), G-E (all products), G-K (23 products), G-X (unit test), G-R (887 / 425, 22 / 14, 97 / 64, 92 / 64), G-F (0 mismatches in 1,098 bridge
rows) passed. Two instrument fixes (Amendments 2 and 3), both before any held-out statement. Held-out products: A119b 18 chromosomes, testis, KB3781, OR6737 (rest, descriptive).

**Held-out, WINDOW reading** (windows found / NULL maximum over 5 random-junction draws; `tables/score.*.json`, full tables in `tables/report.md`):

| substrate (windows) | definition | calls | recall | precision | NULL max precision | found / NULL max | single genes cut |
|---|---|---|---|---|---|---|---|
| A119b HELD (3,899) | V1 | 807 | .018 | .100 | .054 | 69 / 37 | 174 |
| | R | 5,755 | .095 | .067 | .014 | 371 / 78 | 433 |
| | **R2** | 1,146 | .067 | **.238** | .049 | 262 / 54 | 143 |
| | WA | 32,217 | .097 | .027 | .026 | 378 / 353 | 813 |
| | R2 ∪ WA | 33,278 | .161 | .035 | .026 | 629 / 379 | 939 |
| testis HELD (340) | V1 | 23 | .012 | .235 | .118 | 4 / 2 | 11 |
| | R | 35 | .056 | .576 | .212 | 19 / 7 | 6 |
| | **R2** | 17 | .038 | **.867** | .333 | 13 / 5 | 1 |
| | WA | 593 | .079 | .060 | .045 | 27 / 16 | 119 |
| | R2 ∪ WA | 609 | .115 | .080 | .046 | 39 / 18 | 120 |
| KB3781 HELD (898) | V1 | 134 | .033 | .250 | .075 | 30 / 9 | 35 |
| | R | 236 | .089 | .342 | .051 | 80 / 12 | 36 |
| | **R2** | 91 | .077 | **.767** | .078 | 69 / 7 | 12 |
| | WA | 6,629 | .126 | .042 | .037 | 113 / 89 | 269 |
| | R2 ∪ WA | 6,717 | .199 | .052 | .038 | 179 / 98 | 280 |

STRICT (exact windows): R2 precision .215 / .867 / .767, V1 .070 / .235 / .062, WA .007 / .036 / .008. FP_N per 1,000 negatives: R2 5.7 / 0.1 / 0.3, WA 183 / 21 / 76.

**Strata (window reading, found / units [NULL max]).**

| stratum | A119b HELD V1 / R / R2 / WA | testis V1 / R / R2 / WA | KB3781 V1 / R / R2 / WA |
|---|---|---|---|
| dominant link | 8/313 [4] / 109/313 [22] / 0 / 15/313 [12] | (8 units) | 4/23 [1] / 11/23 [4] / 0 / 11/23 [8] |
| co-duplicated partner | 13/530 [20] / 51/530 [15] / 7/530 [4] / 203/530 [199] | 0/26 / 0 / 0 / 15/26 [7] | 4/89 [1] / 6/89 [2] / 1/89 [0] / 28/89 [19] |
| no standalone side | 25/1,356 [13] / 9 [3] / 5 [1] / 120 [121] | 3/252 [2] / 0 / 0 / 23 [12] | 26/540 [8] / 0 / 0 / 72 [65] |

**Plain answer.** No registered definition covers the dominant, co-duplicated or no-standalone-read cases at a usable precision. R2 (the shipped F1v2 bridge) has the only usable precision
(.24 / .87 / .77) and covers only its designed stratum (minority links with standalone reads on both sides: recall .61 / .68 / .77 of those windows); by construction it finds no dominant link and no side without standalone reads.
R finds the most dominant links (A119b 109 of 313, KB3781 11 of 23) at precision .067 / .342 and 433 / 36 single genes cut (V1 finds 8 of 313 and 4 of 23, W 15 and 11). V1 finds 1-3% of windows and, with the container's empty co-duplicated blocks, at most chance on co-duplicated partners.
W (WA) reaches co-duplicated (38% / 58% / 31% of windows) and no-standalone windows (9-13%) and, on the chr16 dev region, all 19 NPIPB14P | PDXDC2P and 3 PKD1P3 | NPIPA1 transcripts where V1, R and R2 call none, but its precision is at chance
(.027 / .060 / .042 against .026 / .045 / .037) with 813 / 119 / 269 single genes cut; its false splits are duplication mosaics (NBPF1, GTF2I, ANKRD36 / ANKRD20A, GOLGA, DPY19L2P2, PMS2P). Combinations add recall and lose precision (R2 ∪ W .035 / .080 / .052).
Run-time needs: V1 families clusters + PAF (after `families`, blind to unclustered loci: 8% / 9% / 5% of transcripts are in clustered loci); R / R2 the BAM and the GTF (`assemble`, shipped); W the GTF, the all-vs-all PAF and the admitted edge graph (families stage, before MCL; no family decision, no reads).

**Predictions (§8).** Hits: 1, 2a, 3, 4, 5, 6 (literal; x1.02), 7a, 8, 9, 10, 11. Misses: 2b (V1 finds 13 of 530 and 4 of 89 co-duplicated windows, not 0), 7b (W cuts 1 of the named twelve on KB3781, not >= 3).
**Falsifiers (§9).** Z1 holds (with LO and RO non-empty, removing the spanning loci separates the halves in 88% of the A119b windows with a fused co-duplicated locus and in 100% of testis and KB3781; 100% in OR6737-rest, descriptive). Z2 fires on every substrate
(W's family-relevant window precision .028 / .113 / .060). Z3 is not falsified literally (W above its NULL by .001-.015 of precision; x1.07, x1.7, x1.3 in windows found). Z4: 13% / 22% / 32% of W's true windows are in loci in no family.

**Read this before quoting it.**
- The registered §7 "detects" clause is weak: WA passes it in most strata with margins of x1.0-1.3 (A119b, KB3781); quote ratios, not stars.
- The truth is annotation-only: a call inside one annotated gene is FP_N even when it marks a duplication-module boundary or a fusion with an unannotated gene. W's false splits are of that kind (post hoc, §5 of the report); the annotation cannot grade them.
- The exact-junction reading excludes 29% of the windows (wide windows, incl. the NPIP-region fusions) and the positives with no junction (PX: 502 / 30 / 402); the window reading is primary for that reason.
- R and R2 carry no fresh evidence (spent for bridge work); V1 and W are new. A119b's V1 and W are lower bounds (per-chromosome families and alignments).
- Post hoc, DEV only (held-out not used): R restricted to transcripts that W flags has precision .257 against R's .100 on A119b DEV (79 calls vs 518, 17 of 84 dominant windows kept, 5 genes cut vs 34); a lead for a new pre-registration, not a result.
- Post hoc (after the held-out scoring, not candidates): W flags transcripts 4-10x enriched for positives (precision over judged transcripts .092 / .112 / .125 against base rates .025 / .013 / .013) but calls 5.0 / 1.8 / 4.5 introns per flagged transcript and places none; constant-free features tried on DEV calls lift junction precision to at most .17.

**Kept products.** `container_defs/tables/` (`score.human_A119b.{DEV,HELD}.json`, `score.human_testis.HELD.json`, `score.gorilla_OR6737.{DEV,REST}.json`, `score.gorilla_KB3781.HELD.json`, `select_w.json`, `report.md`, `report.json`, `truth.*.pkl`), product dirs `dev/`, `held/` (hits, calls, container blocks),
frozen instruments `/mnt/linuxdisk/tmp/rustle_figures/container_defs_frozen/` (`SHA1SUMS` `7cf79b24`). Report `scratchpad/figs/container_defs.md`.

### Independent verification (2026-09-30)

**CONFIRMED EXACTLY on human testis (genome-wide; one substrate).** A separate agent re-implemented, from §1.1, §1.4, §2 and §4 alone, the truth, V1, R, R2, WA and the scoring without opening any `cd_*.py`, `truth_table.py`,
`score.py`, `run_defs.py`, `run_score.py` or result file before it had written its own numbers (`container_defs/verify/mine.json`, `compare.json`).
- **Truth:** labels and windows identical on all 24,802 spliced transcripts (P 337, N 22,683, U0 1,094, U2 688; 303 exact + 37 wide windows; 30 positives without a window; 178 of 209 readthrough parents resolved).
- **WA:** 593 calls on 328 transcripts, the `(transcript, intron)` sets identical (640 of 115,817 introns had LO and RO non-empty; 47 of those were connected and not called).
- **R** 35 calls (22 bridge rows) and **R2** 17 calls (14 kept) identical. **V1** 23 calls identical when folded records are mapped to their representative, as §1.2 says (its first, key-only reading gave 21: 125 transcripts and 2 genuine calls fewer).
- **Scoring:** all 14 compared fields equal for V1, R, R2 and WA under both readings (WA, WINDOW reading: 32 TP [19 exact + 13 AMB], 19 FP_wrong, 483 FP_N, 59 UNJ, 27 of 340 windows found).
- **Its own checks:** locus coordinate convention confirmed from splice motifs (GT-AG at 99.0% of 8,411 introns, 0 at +-1), CIGAR walk identity .92 (+) and .915 (-), a per-column brute-force expansion equals the interval witness sets on 5,784 exons of 804 loci, five called transcripts hand-traced to raw PAF columns.
- **Scope and deviation:** only testis was recomputed independently (the other substrates run the same code paths, gated per product); a handful of the agent's short exploratory scripts (< 15 s, < 2 GB each) ran without `rlock.sh light`; nothing in the repo was edited.

