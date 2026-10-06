# Pre-registration: family-relative CORE-SUPPORT MEMBERSHIP (container v2), a route that needs no fusion detector (dev only)

**Written 2026-09-30 (KEY=csm), before any CSM class, support, core, prune, recruit or null draw was computed on a product and
before any score of a cover arm exists.** Dev only: gorilla fusion simulation S, human A119b chr16 (H) and chr17 (C17), and a
gorilla OR6737 NC_073244.2 census (G). Python prototype in scratch `/mnt/linuxdisk/tmp/rustle_figures_dev/container_csm/`;
nothing in `src/`, `tools/`, `bench/` is edited (`bench/family_container.py` is IMPORTED, frozen sha1 e197ccb3); nothing is
committed or pushed. Species are never pooled. No chimp (PTR) or orangutan (PPY) product is opened and no held-out substrate is
touched (A119b outside chr16 / chr17, testis, KB3781, OR6737 outside NC_073244.2). This stage is DESIGN on dev.

## 0. What was seen before this file (and what was not)

**Seen.** The records of the whole container line: `PREREG_fusion_container_sim_2026-09-28.md` (container v1, Outcome),
`CONTAINER_HEADROOM_2026-09-30.md` (a unit-aware container has +2 copies of headroom: NPIPB5 on chr16 and LOC100420311 on chr17; container v1
is EMPTY for co-duplicated partners: NPIPA1, A6, A7, A9 and TBC1D3G), `figs/container_defs.md` (no junction detector is precise
everywhere), `figs/container_units_mech.md` and `PREREG_container_units_mechanism_2026-09-30.md` (an ORACLE unit split: S M1 7 / 7 / 8 / 9,
H NPIP 25 of 26, Compara F .699, 0 of 49 correct members lost; its bars and scorers), the header of `PREREG_container_units_v2_dev_2026-09-30.md`
(a sibling run, still in progress; none of its products is used), the register rows on closure / core-majority (725, 728, 732, 735-739:
"a family is what half of its members share": inert on NPIP, returns amplicons on chrY, start-dependent, empty denominator at one member), the frozen
container (`bench/family_container.py`), the scorers `score_s.py`, `score_h.py`, `o1_cover_frozen/{score,lo_score,sim_score}.py`, the
`npf_audit` `audit.py` / `npf.py`, the headroom `hl_core.py`.
From products: the stored BASE numbers of the mechanism test (S M1 2 / 0 / 0 / 3 at f = .1 / .5 / .9 / 1, copies in NPIP 24 / 15 / 13 / 13 / 18, M2 0 / 0 / 0 / 3;
H 24 of 26, 30 NPIP members, precision .900, 3 partners, Compara .500 / 1.0 / .667 with 99 of 99 pairs, Liftoff 17 of 36, referee F .236) and
the stored v1 container tables (S `*.container.blocks.tsv`, chr16 / chr17 `D.fam.container.tsv`). The PARTITION only: family sizes
(chr16 111 families, 372 members, sizes 2 x 57, 3 x 27, 4 x 13, 5 x 9, then 6, 7, 9, 28, 30; chr17 78 families, 226 members, 2 x 50, 3 x 15, largest 12 / 11;
S 34-36 families, 13-25 in the largest; gorilla NC_073244.2: 25 families lie wholly on it, 87 members, sizes 2 x 15, 3 x 5, 4 x 2, 5, 8, 21, and 11 families
also have members on other contigs). **Instrument work done before this file, none of it a CSM result:** `csm.py` (prototype, 11 unit tests on a hand-made
fixture, 7 of 7 mutants caught) and gate G-V1 (below), which compares only the V1 labels (support >= 2) with the stored v1 container and prints nothing about a majority.
**Not seen:** any class support, core class, prune, recruit, relation or null of the majority rule on any real product; any cover-aware score of any arm;
the TBC1D3 or NPIP class supports; PTR / PPY / held-out substrates.

## 1. Question and scope

Container v1 (core = one aligned exon column to ANY other member) calls co-duplicated partner pieces core; unit splitting needs a junction detector and no
detector is precise (register 1166D; `container_defs`: R2 .24 / .87 / .77 at recall .07 / .04 / .08, W at chance); the headroom study bounds a unit-aware container at +2
copies. **Question:** can a family-RELATIVE definition of core, the strict-majority consensus of a family's own members, (a) make the container informative for
co-duplicated partners, (b) recover the fused copies that the partition misplaces (NPIPB5, LOC100420311, the S fused copies) by RECRUITING a locus into every family
to whose consensus it contributes a block, and (c) PRUNE members that share nothing with the consensus (partner-only loci dragged in), all WITHOUT a fusion detector
and WITHOUT changing the MCL partition (membership is a cover derived post hoc)? **What it cannot settle:** any default (dev, oracle bars); held-out behaviour;
families of 1-3 members, where "majority" is "all" or "two of three".

## 2. The rule (binding)

Inputs are the products of the families stage (`clusters.tsv`, `loci.tsv`, `loci.gff3`, `loci.paf` with CIGARs, the families-input GTF). Nothing is fitted.

1. **Loci and blocks (v1, extended).** LOCUS = a representative key `CONTIG:START-END` plus every annotation record `loci.tsv` folds into it; a record that is not folded
   is its own representative. A locus is CLUSTERED iff its representative is a `clusters.tsv` member, else UNCLUSTERED (unlike v1, unclustered loci get blocks: a
   recruit comes from outside the family). BLOCKS = the union of the exons of all transcripts of all gene_ids of all records of the locus, merged where they overlap
   (abutting exons stay separate).
2. **Joins (v1's column test).** Block b of locus m and block b' of locus m' (m != m') are JOINED iff one aligned CIGAR column (M / = / X) of one PAF record maps an exon
   base of b onto an exon base of b' (strand-blind; every PAF record counts: no identity, length or primary filter; same-locus records are skipped). Joins with both
   sides unclustered are not needed and are not stored.
3. **Classes and support.** Inside a family F (MEMBERS = its clustered loci, |F| = their number in the evaluated product): nodes = the blocks of its members; classes =
   connected components under the joins between two members of F. SUPPORT s_F(c) = the number of distinct members with a block in c (a locus counts once however many
   blocks it has in c).
4. **Core.** MAJ: class c is CORE iff 2 s_F(c) > |F| (strict majority). V1 (the comparator): c is core iff s_F(c) >= 2, which is container v1's test exactly (a block is joined to a block
   of another member iff its class has two members). A family with NO core class is CORELESS: membership is left as clustered; under MAJ its blocks are UNDETERMINED
   (neither core nor accessory, no accessory call is made); under V1 they are accessory, as in v1.
5. **Cover membership.** In a non-coreless family F: a member is CORE if >= 1 of its blocks is in a core class; otherwise it is PRUNED (PERIPHERAL: reported, not counted as a
   member of F). A locus r that is not a member of F is RECRUITED into F iff >= 1 block of r is joined to a block of a core class of F. A locus may be recruited into several
   families and keeps its partition family (unless pruned from it). Families stay the MCL partition; the cover is derived post hoc and changes no family.
6. **Container v2 output** (per clustered locus and block): class id, support, |F|, class in {core, accessory, undetermined}; for every non-core block, RELATIONS = every other family
   F' such that the block is joined to a block of a CORE class of F'. (v1 related a block to ANY member of F'.) For a cover membership (r, F) the CORE blocks of r for F are the blocks of r joined to a
   core class of F (for a member: its blocks in core classes).

## 3. Variants, comparators, null (declared now; no tuning)

| arm | prune | recruit | core test | role |
|---|---|---|---|---|
| **BASE** | no | no | v1 container (V1 labels) | the current defaults' partition and v1 container |
| **V-maj** | yes | yes | MAJ | the registered rule |
| **PRUNE** | yes | no | MAJ | prune only |
| **RECRUIT** | no | yes | MAJ | recruit only |
| **V1CORE** | yes | yes | V1 (support >= 2) | the known-failure comparator: no majority |
| **NULL-1..5** | same counts as V-maj | same counts as V-maj | - | per family F: as many members pruned as V-maj prunes from F (uniform from F's members), as many loci recruited as V-maj recruits into F (uniform from the pool of loci outside F joined, in any block and class, to a block of a member of F); `random.Random(20260930 + k)`, families and pools in sorted key order; the same draws give NULL-P (prunes only) and NULL-R (recruits only) |
| ORACLE (reference, not re-run) | - | - | - | the mechanism test's oracle unit split (arm A): S M1 7 / 7 / 8 / 9 at f = .1 / .5 / .9 / 1, H NPIP 25 of 26, Compara F .699, 109 pairs, 0 of 49 lost |

PRUNE and RECRUIT decisions are independent (a pruned member has no core block, hence no recruiting join); V-maj = PRUNE + RECRUIT. V1CORE applies the same machinery with s >= 2.

## 4. Inputs, gates

- **S** (gorilla simulation, 10 fusions, f = 0, .1, .5, .9, 1): the BASE products of `container_units_mech/S/f<f>/BASE/` (`f<f>.fam.*`, `f<f>.families.gtf`, `f<f>.container.*`), produced by the e163d955 defaults
  (`copy_assign --assemble-only` with `--bridge-regroup f1v2`, `mcl_families --min-cov-shorter 0.70 --emit-container`). **H** (human A119b chr16): `container_units_mech/H/BASE/hsa16.*`
  (byte-identical to `container_headroom/data/human_A119b/fam/chr16/D.*`, checked). **C17** (chr17): `container_headroom/data/human_A119b/fam/chr17/D.*` with `asm_f1v2/chr17/*.families.gtf`.
  **G** (gorilla OR6737): `container_headroom/data/gorilla_OR6737/fam/c5/D.*` FILTERED AT FILE LEVEL (awk, before Python) to NC_073244.2: GTF and `loci.gff3` lines on the contig, `loci.tsv` and PAF records whose names are both on
  it, `clusters.tsv` rows of the 25 families that lie WHOLLY on it (the 11 cross-contig families are not evaluated: their majority would involve other contigs; their loci are read as unclustered). Their replays are
  byte-identical per the headroom NOTES.
- **G0** `sha1sum` of `bench/family_container.py` (e197ccb3), `tools/rlock.sh` (30f424a9), the copied instruments. **G1** `test_csm.py` passes (11 tests) and the 7-mutant check catches 7 of 7.
  **G-V1** (DONE before this file): `csm.py`'s V1 labels equal the stored v1 container block for block, class and geometry: S f0.5 1,320 of 1,320, chr16 1,675 of 1,675, chr17 1,219 of 1,219 (0 mismatches); run on the other S f before scoring.
  **G-S** the S scorer with an EMPTY membership reproduces `results/S_scores.json` BASE (copies in NPIP, M1, M1u, M2, M3v1 precision / recall, at every f). **G-H** the H scorer with an empty membership reproduces
  the stored BASE row (copies 24, fused 10, 30 members, precision .900, 3 partners, Compara .500 / 1.0 / .667, pairs 99 / 99, Liftoff 17 of 36) and its Python referee equals the frozen `family_score` binary (sens .134, prec .976, F .236, pairs 139 / 140 / 1474).
  **G-C17** chr17 BASE: TBC1D3 11 of 16, Compara .320 / .923 / .475, pairs 44 / 47. **G-N** the NULL draws are deterministic and have V-maj's per-family counts. **G-F** the G files contain no other contig.
  A failed gate stops that substrate: nothing is scored through it.

## 5. Measures (dev; per substrate and arm; species never pooled)

### 5.1 S (scorer logic of `score_s.py`, holders / NPIP family / placement from `npf_audit` exactly as the mechanism test; cover adaptation below)

- Cover adaptation (documented): the HOLDER of a copy / partner is unchanged (the locus with the most same-strand exon bp over the copy territory / partner exon union); the copy is IN NPIP iff the NPIP cluster id (BASE's) is in the cover families of
  the holder's locus (partition family unless pruned, plus recruited).
- **M1** fused copies in NPIP /10; **M1u** unfused /15; copies in NPIP /25; copies LOST = copies in NPIP under BASE not under the arm (by name).
- **M2 (literal)** S partners /20 (and the 10 fused partners) whose holder locus is in the NPIP cover. **M2c (primary partner-leak measure for a cover arm)** = partners whose holder locus is in the NPIP cover AND whose exon union has >= 1 base in a block of that
  locus that is CORE for the (locus, NPIP) membership. Reason: a cover puts a fused locus in both families by design, so the literal M2 counts every recruited fused locus that absorbed its partner and is >= the number of such recruits by construction; M2c asks whether the PARTNER's sequence is carried as NPIP core. Both are reported; the bars use M2c.
- **M3v1 / M4** as `fusion_container_sim` / `score_s_m45.py` with the arm's core definition: accessory blocks of the fused member (the carrier of the fusion junction, else the copy's holder) wrt its partition family; precision = |A n P_all| / |A|, recall = |A n P_fus| / |P_fus|; M4 = whether its accessory blocks relate (v2 relation) to the partition family of the partner's holder. **M3n**
  (NPIP-relative, cover arms): the same with A = the blocks of the member that are not CORE for (member, NPIP), for members in the NPIP cover. M5 = accessory bp on unfused holders.
- **Relation P / R (lenient)** as the mechanism test's M3u with the cover: one record per fused pair that has a carrier locus; CORRECT iff the carrier's cover families contain a cluster that MATCHES the copy's reference family F*_c (BASE f = 0; threshold-free co-membership, the test's `match()`), contain a cluster that matches the partner's reference family F*_p (an empty
  F*_p is matched iff no NPIP-core block of the carrier holds the partner's bases), and the partner is not carried as NPIP core (M2c = 0 for the pair). Precision = correct / emitted; recall = correct pairs / 10.
- Per-f census: families, coreless, pruned, recruited (by family), what the recruits are (copies, partner loci, other).

### 5.2 H (human chr16; scorer logic of `score_h.py` and `o1_cover_frozen/score.py`; cover adaptations documented and gated)

- **NPIP copies /26 and which** (the 26 Dishuck copies, `audit.human_truth()`; holder rule as `npip_dev.py`), by cover; **NPIP members** = BASE members minus pruned plus recruits; **NPIP precision (unit = locus)** = members whose same-strand exons overlap >= 1 truth copy / members; **non-copy members** (partners) listed with `gene_at` labels; the 11 FUSED_AUDIT copies, the four
  readthrough copies (NPIPA1, NPIPA6, NPIPA9, NPIPB14P), NPIPB5, PKD1P6-NPIPP1 (its NPIPP1 half).
- **Compara** (Primates families, `--chrom ALL` semantics with contigs other than chr16 dropped; `o1_cover_frozen/score.py` gene space, gate G-H against the binary): bipartite sens / prec / F (scipy assignment; pairwise counts beside) and pairwise TP / predicted / truth.
  Cover gene sets: a family's set = `gene_at` labels of its non-pruned members' representatives (as `core_sets`) plus one label per recruit. **PRIMARY label (the o1_cover convention for attached members, "a peripheral member is correct iff the truth lists one of its parents in that family")**: the recruit's PARENTS =
  the annotated genes whose span overlaps >= 1 base of a block of the recruit that is joined to a core class of the family; label = (1) a parent already in the family's gene set, else (2) a parent in the truth universe whose truth family is the family's plurality truth family, else (3) the locus's own `gene_at`. **STRICT label** = the locus's own `gene_at` (the BASE convention, no benefit of the doubt): reported beside every
  PRIMARY number. The attachment judge of `score.py` (own / truth / correct, judged iff the family has a labelled core member) gives the recruit precision. Pairs from recruited members are counted as pairs of their label with the family's other labels (set semantics, as every pair in the scorer).
  Predictions are intersected with the truth universe (register 991 / 992 trap): **unjudgeable recruits** (label outside the universe) are counted and reported, since they cannot hurt a universe-limited precision.
- **Liftoff** copy-pair recall: rows = non-pruned members (representative exon union) plus recruits (once per family), `lo_score.py`. **Referee** (protein-homology referee, `fig7/current/human_chr16_ref.*`): F by a Python re-implementation of the same gene-space scorer (PRIMARY / STRICT), gated to the frozen binary on BASE.
- **Correct members lost** (mechanism §6.2): a correct member of a BASE cluster = a gene labelled by `gene_at` that shares a Compara truth family with >= 1 other gene of that cluster; LOST iff no arm family holds it together with any of those mates (cover sets). Also the NPIP copies lost and the BASE NPIP copy-members pruned.
- Container v2 on chr16: blocks core -> accessory / undetermined relative to v1, the relations, the containers of NPIPA1 / A6 / A7 / A9 / B14P / B3 / LOC128966608 / B4, NPIPB5.

### 5.3 C17 (human chr17, TBC1D3)
TBC1D3 copies /16 and LOC100420311 (headroom `hl_core` fate rows, holder rule unchanged, cover placement), TBC1D3 members and precision (units overlapping a truth copy / members; BASE 11 of 11), Compara chr17 sens / prec / F and pairs (gate 44 / 47) under PRIMARY and STRICT, correct members lost.

### 5.4 G (gorilla NC_073244.2; census only)
Families wholly on the contig: coreless, pruned, recruited, what they are (annotation from `GGO_genomic.gff`, RepeatMasker). No family-level truth exists for gorilla except the copy sets, which are off this contig (NPIP T_member: 1 copy here, on a cross-contig family; TBC1D3: none): the NPIP / TBC1D3 membership of gorilla is measured by S only, and the G census is not a score.

### 5.5 Safety census (every dev contig: chr16, chr17, NC_073244.2, S per f)
Families coreless (by |F|); members pruned and recruited (counts, per family, by |F| stratum 2 / 3 / 4-5 / >= 6); recruits per family; **what the recruited and pruned loci ARE**: Compara truth family and whether it is the target family's plurality truth family; single-exon stubs (every transcript of the locus has one exon); repeat-derived (>= 50% of the
joined-block bases (recruits) or of all block bases (prunes) covered by RepeatMasker interspersed repeats; Alu separately; `hs1.repeatMasker.out.gz`, `GCF_029281585.2.repeatMasker.out`; the 50% is a description, not a rule); pseudogene fragments (overlap an annotated pseudogene record); fused / readthrough loci (exons over >= 2 annotated records with disjoint spans);
recruits that OVERLAP a member of the target family on the genome (same-place joins); pruned GROUPS (>= 2 pruned members of one family sharing a class: a divergent subfamily signature, e.g. NPIPA vs NPIPB). The **support histogram** of classes in NPIP and in the two other largest chr16 families; the number of classes whose blocks include two blocks of one locus (class chaining). **A per-locus table for every recruit and prune on chr16** with its truth status
(`tables/chr16_recruits_prunes.tsv`).

## 6. Decision rule and verdict mapping (integer comparisons; fixed now)

**Lexicographic rule over admissible arms.** (1) ADMISSIBLE iff zero correct members are lost on S (copies in NPIP under BASE and not under the arm, at every f), on H (Dishuck copies in NPIP, Compara-correct members by PRIMARY cover sets) and on C17 (TBC1D3 copies, Compara-correct members). An arm with a lost member is ranked below every admissible arm and its verdict is NOT. (2) Among admissible arms the larger H NPIP copies, then C17 TBC1D3 copies,
then H Compara true pairs (PRIMARY), then S sum over f = .1 / .5 / .9 / 1 of M1 (BASE: 2 + 0 + 0 + 3 = 5). (3) Then fewer recruits plus prunes (all dev contigs, S excluded).

**Clauses (per arm; S, H, C17 rows separate; nothing pooled).**
- **A (admissible)** as (1).
- **D (guards).** S: M1u >= BASE's at every f. H: non-copy NPIP members <= BASE's (3) and Compara pairwise precision (PRIMARY) >= BASE's - 0.01. C17: non-copy TBC1D3 members <= BASE's (0).
- **E (gain over BASE).** strictly above BASE in at least one of: S sum of M1 (5), H NPIP copies (24), C17 TBC1D3 copies (11), H Compara true pairs (99); OR a CLEANING gain: H non-copy NPIP members < 3 or C17 non-copy TBC1D3 members < 0 (impossible) with the headline measures not below BASE.
- **B (oracle bars, S).** M1 >= 9 at every f in {.1, .5, .9, 1} AND M2c = 0 at every f (the oracle bar; the ceiling was 9 for a partition because NPIPA7 is outside NPIP at f = 0, a cover could exceed it). **C (oracle bars, H and C17).** H NPIP copies >= 25 AND Compara F (PRIMARY) >= .699; C17 TBC1D3 copies >= 12.
- **WORKS** iff A and D and B and C. **PARTIAL** iff A and D and E and not WORKS (the failing bar clauses are named). **NOT** otherwise (a correct member lost, a guard failed, or no gain).
- **Attribution clause.** A gain (clause E item) is attributed to the core-support criterion only if no NULL seed also reaches it; otherwise it is reported as a recruiting / pruning effect, not as core support. NULL arms are expected to clear nothing; they are scored with the same clauses.
- Reported beside, never decisive: STRICT labels, literal M2, M3v1 / M3n / M4, relation P / R, Liftoff, referee F, NPIP precision, census.

## 7. Predictions (prior probabilities, stated before any CSM product) and falsifiers

P1. Gates G-S, G-H, G-C17 pass on the first attempt (0.85); the referee re-implementation matches the binary (0.85).
P2. **V-maj S**: M1 >= 9 at all four f (0.55); sum of M1 above BASE's 5 (0.92); M2c = 0 at every f (0.75); literal M2 > 0 at some f (0.95, by construction); M1u >= BASE's at every f (0.75).
P3. **V-maj H**: NPIPB5 is in NPIP (0.75); NPIP copies = 25 (0.6), >= 26 (0.05); the four readthrough copies all stay in NPIP (0.8); PKD1P6-NPIPP1 in NPIP (0.10); zero correct members lost (0.6; PRUNE alone 0.6, RECRUIT alone 0.97).
P4. **V-maj C17**: LOC100420311 recruited into TBC1D3 (0.65); TBC1D3 copies >= 12 (0.6).
P5. **Census scale.** Coreless families on chr16 between 3 and 20 of 111 (0.75); members pruned on chr16 between 1 and 25 (0.75); recruits (locus, family) on chr16 >= 20 (0.6), >= 100 (0.25); V1CORE recruits >= 3 x V-maj's (0.85) and V1CORE prunes <= V-maj's (0.9); NPIP non-copy members after V-maj exceed 3 (0.55); recruits into NPIP on chr16 number 2-8 (0.6).
P6. **What they are.** >= 50% of pruned loci on chr16 are single-exon stubs or pseudogene fragments (0.6); >= 1 chr16 recruit is repeat-derived by the description above (0.6); >= 1 pruned group of >= 2 members in one family (0.5); no NPIPA clade member is pruned (0.8).
P7. **Verdicts.** V-maj: WORKS 0.10, PARTIAL 0.45, NOT 0.45. PRUNE: NOT 0.65 (no gain: it cannot add a copy or a pair), PARTIAL 0.30 (cleaning). RECRUIT: PARTIAL 0.45, NOT 0.45, WORKS 0.10. V1CORE: NOT 0.90. NULL: no seed clears clause B or C (0.95); some seed shows a gain item of clause E (0.25).
P8. The container v2 fills v1's empty containers: NPIPA1, A6, A9 (and A7) have a non-empty accessory set under MAJ (0.75); the accessory set of NPIPA1 / A6 / A9 relates to the PKD1P family (0.55).
**Falsifiers.** Z1 V-maj or PRUNE loses a correct member (the majority core removes real divergent members: heterogeneous families). Z2 a NULL seed reaches V-maj's gain (the gain is recruiting anything joined, not core support). Z3 V1CORE is indistinguishable from V-maj (the majority adds nothing). Z4 M2c > 0 at any f for V-maj (the container does not separate partner sequence). Z5 V-maj
recruits >= 10% of all loci on a contig (over-recruiting through repeat-derived or same-place joins). Z6 a pruned group is a real subfamily (NPIPA vs NPIPB type) in any family.

## 8. Hostile self-review

1. **Repeat classes.** Alu- or L1-derived exons form giant classes with high support inside families of repeat-rich loci; a repeat class can be the majority "core" and recruit every locus with a repeat-derived exon that the aligner reports (the PAF keeps <= 50 targets per query, `-N 50 -p 0.1`). The census reports repeat fraction per recruit; no filter is part of the rule, so a repeat-driven result is a finding about the rule.
2. **Heterogeneous families.** A strict majority prunes the minority clade of a two-clade family (NPIPA vs NPIPB, ABCC-type paralogs); the rule has no clade notion (register 735-739: the closure is start-dependent and non-unique on palindromes). The pruned-group census is the detector; Z1 and Z6 are its falsifiers.
3. **Tiny families.** 57 of 111 chr16 families and 50 of 78 chr17 families have two members, where "majority" is "both": the rule can neither prune nor rank there, only recruit, and a pair joined at interval level but not at column level is CORELESS (container v1's smoke: 21 of 338 testis families had no core block). Families of three need two members to agree. The census is stratified by |F|.
4. **Class chaining.** Blocks are unions of the exons of ALL transcripts of a locus, so one long block (a readthrough locus) can connect classes that are separate elsewhere; connected components then inflate support. Reported as classes holding two blocks of one locus and as support > |F| / 2 classes per family; not corrected.
5. **Alignment coverage.** 81% of referee same-family pairs have no alignment at all (r1026): a member without an alignment to the consensus is pruned although real; the PAF caps targets per query at 50, so a locus in a 30-member family competes for slots with non-members. Counted inside "correct members lost".
6. **Same-place joins.** Strand-blind joins between loci that overlap on the genome (an antisense or alternative-isoform locus) align to themselves; such pairs recruit each other. Reported (recruits overlapping a member on the genome); not excluded.
7. **Recruit != member.** A recruit is evidence of a shared exon module, not of homology at gene level: duplicated exon modules shared by two families (segmental duplications) produce recruits in both directions by design; the cover counts them.
8. **The scorers.** PRIMARY labels give the benefit of the doubt to a recruit whose parent is in the target truth family (o1_cover's convention); STRICT does not; both are printed. Compara is universe-limited, so recruits of genes outside the universe are free to it: the NPIP / TBC1D3 non-copy counts are the counterweight (copy truth is not universe-limited).
9. **M2.** The literal M2 cannot be 0 for a cover that recruits fused loci; the bar uses M2c. A reader who insists on the literal M2 gets NOT for every recruiting arm by construction (stated so before the product).
10. **Dev only.** chr16 is the development block of every fusion rule; S has 10 fusions, annotation-built error-free reads (circular); the gorilla real-data piece is a census of families wholly on one contig (NPIP and TBC1D3 of gorilla are off it); no held-out number exists, so nothing here supports a default.
11. **The bar numbers are mine** (clause constants 9, 25, .699, 12, 3, .01 are the oracle's and BASE's own numbers, nothing is fitted); the rule's only constants are the majority 1/2 (a definition) and, in V1CORE, 2 (the comparator).

## 9. Order, machine rules, files

1. This file; then (Amendment 1) the frozen instruments with sha1, the gate results G-S / G-H / G-C17 / G-N / G-F, and every implementation ambiguity resolved before scoring; then the runs (S, H, C17, G; census) in that order; then Outcome. Light jobs only (`tools/rlock.sh light`; every step is < 2 GB and seconds to minutes), foreground, never `pkill -f`, TMPDIR and scratch under
   `/mnt/linuxdisk/tmp/rustle_figures_dev/container_csm/`. Register rows, suffix L, NOT appended. Report `figs/container_csm.md` in the session scratchpad. Nothing in `src/`, `tools/`, `bench/` is edited; nothing is committed or pushed.

## Amendments

**Amendment 1 (2026-09-30 14:45, after the instruments and the gates, before ANY CSM outcome was computed or printed on a real product).** Scratch
`/mnt/linuxdisk/tmp/rustle_figures_dev/container_csm/` (`lib/`, `out/`, `tables/`, `G/`, `logs/`); frozen copy of the instruments `/mnt/linuxdisk/tmp/rustle_figures/container_csm_frozen/` with `SHA1SUMS`. Prereg sha1 before this
amendment `a95db421`. Instrument sha1 (first 8): `csm.py` c3d7141e, `common.py` 4e1a8f8b, `score_s_csm.py` 2a40b17b, `score_h_csm.py` 9ddb91bf, `score_c17_csm.py` 7f306456, `census.py` 2f8d9f58, `verdict.py` 9fa5d501, `run_all.py` 8a6b4cd8,
`gate_v1.py` 556110b4, `gate_s.py` 3abbda5c, `gate_n.py` 33f16942, `test_csm.py` 8fb56971, `test_pipeline.py` 7076483f, `mutants.py` e46d08a8, `smoke.py` da2a3533. Imported unchanged: `bench/family_container.py` e197ccb3, `container_units_mech/lib/{score_s,score_h,score_s_all}.py`,
`o1_cover_frozen/{score,lo_score,sim_score}.py`, `npf_audit/{audit,npf}.py` (npf fdc1a7d3, audit fd0305ca), `container_headroom/code/{hl_core,fate,cells}.py`, `tools/rlock.sh` 30f424a9.

*Gates, all passed (outputs in `out/gates_g0_g1.txt`, `out/S_gate_base.json`, `out/H_gate_base.json`, `out/C17_gate_base.json`):*
- **G0** family_container.py e197ccb3, rlock.sh 30f424a9. **G1** `test_csm.py` 11 tests and `test_pipeline.py` 15 tests (census, tables, verdict logic on synthetic dicts) pass; **7 of 7** mutants of `csm.py` are caught.
- **G-V1** V1 labels == stored container v1, class and block geometry, block for block: S f0.0 1,331 / f0.1 1,322 / f0.5 1,320 / f0.9 1,321 / f1.0 1,342, chr16 1,675, chr17 1,219 (0 mismatches everywhere). So V1CORE's class labels ARE container v1's.
- **G-S** the S scorer with the EMPTY membership reproduces `S_scores.json` BASE at all five f: copies in NPIP 24 / 15 / 13 / 13 / 18, M1 9 / 2 / 0 / 0 / 3, M1u 15 / 13 / 13 / 13 / 15, M2 0 / 0 / 0 / 0 / 3 (fused 3), M3v1 precision / recall and accessory bp (.2439 / .3984, .215 / .4174, .2295 / .4174, .4092 / .6212), M4 and M5 (`S_m45.json`).
  **G-S2** (smoke) M3v1 through the V1-mode `blocks_dict` path equals the stored one (f0.5: .215 / .4174).
- **G-H** the H scorer with the empty membership reproduces the stored BASE row under BOTH label conventions: Compara 17 truth families / 54 genes / 20 clusters, sens .500 prec 1.000 F .6667, pairs 99 / 99 of 193, matched 27, exact 5; copies 24, fused 10, 30 members, 27 copy members, precision .900, 3 partners; Liftoff 17 of 36; and the Python referee equals the frozen `family_score` binary
  (80 families / 306 genes, sens .134 prec .976 F .236, pairs 139 / 140 of 1,474).
- **G-C17** chr17 BASE: TBC1D3 11 of 16, 11 members (11 copy members), LOC100420311 outside; Compara 25 families / 75 genes, sens .320 prec .923 F .475, pairs 44 / 47 of 109; Liftoff 1 of 14 (all equal the stored D row).
- **G-N** for every dev product the five NULL draws are deterministic, have V-maj's per-family prune and recruit counts, prune only within the family and recruit only from its touch pool; PRUNE and RECRUIT are the two halves of V-maj (asserted, no counts printed). **G-F** the G files (GTF, loci.gff3, loci.tsv, PAF, clusters.tsv, RepeatMasker rows) contain only NC_073244.2 (0 lines on another contig; 25 whole-contig families, 87 members;
  the file-level filter is awk, before Python).
- **Smoke** (`smoke.py`): every scoring and census path ran on the real products with a SYNTHETIC membership (random prunes / recruits, `random.Random(7)`, not a CSM rule, not one of the registered seeds) printing only structural checks; it is the only use of a non-registered membership and no outcome was printed.

*Implementation ambiguities resolved (binding from here on):*
1. `|F|` = the number of `clusters.tsv` rows of the family in the evaluated product (G: after the whole-contig filter). A locus with several blocks counts once per class. The class root is the smallest (locus key, block) node.
2. A fold target that is not a `clusters.tsv` member is an UNCLUSTERED locus and carries its folded records' exons; fold chains are resolved (none occur). The 11 cross-contig gorilla families (G) are unclustered in the G product.
3. Coreless handling: under MAJ, `class` = `undetermined` for every non-core block of a coreless family (no accessory call, no M3 accessory bases, no prune); relations are still computed for every non-core block. Under V1 the label is v1's (accessory).
4. M3v1 / M4 use the container of the arm's MODE (MAJ for V-maj, PRUNE, RECRUIT, NULL; V1 for V1CORE; the stored v1 container for BASE); prune / recruit do not change a label. The V1CORE relation is "joined to a block of a core class (support >= 2) of another family", which is not v1's "any block" (BASE uses the stored v1 container).
5. A recruit's CORE blocks for (r, F) are the blocks of r joined to a block of a core class of F, recomputed from the joins for NULL recruits too (a NULL recruit with no such join has none, so its M2c is 0 by definition).
6. M2c: a partner counts iff the holder locus is in the NPIP cover AND >= 1 base of the partner's all-transcript exon union lies in a core block of (holder locus, NPIP). Relation P / R (lenient): as §5.1, `match()` of the mechanism test on the BASE f = 0 reference families, an empty F*_p is satisfied iff the partner is not carried as NPIP core; reported, not in the verdict.
7. PRIMARY parents = every gene of the labeller with overlap > 0 (family_score's `min(e, ge) - max(s, gs) > 0`) with a block of the recruit joined to a core class of the target family; label rules (1)-(3) as §5.2 with the family's pruned core gene set; the plurality truth family is taken over the core genes in the universe after pruning.
   STRICT = the representative key's own `gene_at`. The judge's `jall` = labels of every record of the non-pruned members (folded records included).
8. Lost members (H, C17): the mechanism rule on the cover gene sets of the arm vs BASE's; H and C17 clause A also require no Dishuck / TBC1D3 copy lost. NPIP family id = BASE's (most copies); TBC1D3 family id = BASE's `fam_cl` of the headroom fate table.
9. NPIP / TBC1D3 members = BASE members minus pruned plus recruits; copy members = members whose same-strand exons (gene_ids of the representative key) overlap >= 1 truth copy; non-copy members = the rest (listed with `gene_at` labels).
10. Clause thresholds: C uses `round(F, 3) >= .699` (the oracle's .6988), H NPIP copies >= 25, C17 >= 12; E compares integers with BASE's gated values (S sum M1 = 5, H c1 = 24, C17 = 11, H pairs = 99). The NULL arms are the five seeds 20260931-35 (`random.Random(20260930 + k)`), and NULL-P / NULL-R are the same draws split.
11. Census repeat classes: RepeatMasker rows of class SINE, LINE, LTR, DNA, RC, Retroposon are "interspersed"; `SINE/Alu` separately; `repeat_derived` = >= 50% of the joined (recruit) or all (prune) block bases covered (a description). Pseudogene = an annotated record with `pseudogene` in its biotype overlapping an exon; fused = exons over >= 2 annotated records with disjoint spans.
12. The G census scores nothing against a truth (no family truth for gorilla off the copy sets): counts, what the loci are (annotation, repeats), and the support histograms of the three largest families; recruits there can only come from loci on the contig.

*Order.* `run_all.py s`, `h`, `c17`, `g`, then `verdict`; the first outputs are final. A defect found after the first scoring will be reported as a deviation with both numbers (as the mechanism test did), never by changing a bar, arm or threshold.

## Outcome (2026-09-30 15:20; scored with the Amendment 1 definitions, no deviation after it; full report scratchpad `figs/container_csm.md`, data `/mnt/linuxdisk/tmp/rustle_figures_dev/container_csm/`)

**Registered verdicts (§6, integer clauses, `out/verdict.json`, `tables/verdict.md`): V-maj NOT, RECRUIT NOT, V1CORE NOT, PRUNE PARTIAL, NULL 1-5 NOT.** V-maj is ADMISSIBLE (clause A holds on S, H and C17: no copy lost, Compara lost 0 of 49) and exceeds the oracle unit split's headline on every dev cell, but fails clause D (guards) on H and C17:

| measure (dev; species separate) | BASE | V-maj | oracle (unit split) | NULL 1-5 |
|---|---|---|---|---|
| S M1 /10 at f = 0 / .1 / .5 / .9 / 1 (sum over f > 0) | 9 / 2 / 0 / 0 / 3 (5) | **10 / 10 / 10 / 10 / 10 (40)** | - / 7 / 7 / 8 / 9 | sum 35 / 31 / 27 / 40 / 37 |
| S M2c (partner carried as NPIP core) / literal M2 at f >= .1 | 0 / 0 0 0 3 | **0 / 6 9 9 10** | 0 | 0 / 2-13 |
| S unfused copies in NPIP /15, copies lost | 15 13 13 13 15, - | 15 at every f, none | - | 14-15, 1 in 4 of 5 seeds at f = 0 |
| H NPIP copies /26; members; non-copy members; precision | 24; 30; 3; .900 | **26**; 44; **15**; .659 | 25 | 22-23; 44; 17-19; .57-.61 |
| H Compara PRIMARY F, pairs (STRICT F, pairs); lost members | .667, 99 (.667, 99); - | **.759, 168 (.736, 131)**; 0 | .699, 109 | .70-.76, 161-165 (.64-.70, 110-127); 0 |
| H Liftoff /36; referee F | 17; .236 | 20; .295 | - | 20; .270-.290 |
| C17 TBC1D3 /16; non-copy members; Compara PRIMARY pairs | 11; 0; 44 / 47 | **13**; **6**; 48 / 71 | 12 | 12-13; 7-10; 48 / 73-82 |
| recruit pairs (loci, % of loci) chr16 / chr17 / gorilla NC_073244.2 | 0 | 1,076 (355, 15.1%) / 537 (225, 9.3%) / 206 (68, 6.4%) | - | same counts |

- **Clauses.** D_H fails (non-copy NPIP members 15 > 3; pairwise precision holds, 1.000), D_C fails (6 > 0); B (S oracle bars: M1 >= 9 at f > 0 and M2c = 0) and C (H NPIP >= 25, F >= .699; C17 >= 12) hold. V1CORE also fails B (M2c 1 / 1 / 1 / 5 at f >= .1). PRUNE: A and D hold, E only through the cleaning item (NPIP non-copy 3 -> 2, precision .900 -> .929), B and C fail; its 2 prunes are NPIP fragment loci (`chr16:22353087-22784890`, overlapping NPIPB5's territory,
  and `chr16:30536712-30638496`), 0 on chr17 and gorilla. **Attribution:** S sum M1, C17 copies and H pairs are reached by all five NULL seeds (Z2); the H NPIP copy gain (26) by none (22-23).
- **Why (mechanism, read from the products).** The majority keeps 93% of the support >= 2 recruits (1,076 of 1,154 pairs on chr16) and 78% of the pool of loci joined to a family at all (1,382 pairs), so the null and V1CORE are nearly the same set. 885 of the 1,076 chr16 recruits enter through ONE joined block; 569 are repeat-derived, 417 single-exon stubs, 234 overlap a pseudogene record, 366 are unclustered loci (1 truth-consistent), 226 overlap a member on the genome.
  Family sizes: 57 of 111 chr16 families have two members (majority = both); 6 coreless families, all pairs; 25 families are all-single-block and take 24% of the recruits. Support histogram of NPIP (MCL1, 30): {1: 45, 2: 31, 3: 10, 5: 3, 6: 1, 10: 1, 24: 3, 27: 1, 28: 1}, five core classes (support 24-28); TBC1D3 (MCL3, 11): {1: 5, 2: 7, 5: 1, 11: 9}.
- **Predictions (§7).** Hit: P1, P2 (all five), P3 NPIPB5 / readthrough copies / zero lost, P4, P5 coreless 6 / pruned 2 / recruits >= 100 / V1CORE prunes <= V-maj's / NPIP non-copy > 3, P6 repeat-derived recruits / a pruned group of 2 / no NPIPA clade pruned, P7 V-maj, RECRUIT and V1CORE NOT, P8 (v1's empty containers filled, related to the PKD1P family). Miss: NPIP copies = 25 (26: PKD1P6-NPIPP1 too, p = .10), V1CORE recruits >= 3x V-maj's (1.07x),
  recruits into NPIP 2-8 (16), pruned loci stubs / pseudogenes (0 of 2), PRUNE NOT (PARTIAL), "no NULL seed clears B" (2 of 5 do on S), "no NULL gain item" (all five have one). **Falsifiers fired:** Z2, Z3 (V1CORE ~ V-maj), Z5 (15.1% of chr16 loci recruited). Not fired: Z1, Z4, Z6.
- **Independent verification.** A numpy per-column re-implementation written from §2 without importing `csm.py` or `family_container.py` reproduces the joins, class supports, coreless families, prunes and recruits exactly on chr16 (9,058 joins, 1,076 recruits), chr17 (5,648, 537), gorilla NC_073244.2 (1,448, 206) and S f = 0 / .5 / 1.
- **Post hoc (labelled, not registered, dev only; do not cite as tests).** (i) Partner bases (RefSeq partner records of the fused copies that sit in their family) called accessory: chr16 8 copies 22.7% (v1) -> **82.6%** (v2), the copy's own bases 0% -> 6.4%; chr17 2 copies 20.9% -> 100%, 0% -> 2.2%; NPIPA1 / A6 / A9 go from an empty container to 13-15 accessory blocks related to the PKD1P family (MCL25).
  (ii) The relation channel (an accessory block of a CLUSTERED locus joined to a CORE class of another family) relates exactly 2 chr16 loci to NPIP, PKD1P6-NPIPP1 and NPIPB5, both true copies, 0 false; on S it finds 7-11 of the copies BASE misses (2 false loci); on chr17 it gives TBC1D3G / TBC1D3D fragment loci and USP6 and misses LOC100420311 / TBC1D3P5 (a 2-member family whose core is the module). Scored as a membership list with the registered clauses (REL): S PARTIAL, H PARTIAL, C17 NOT (one non-copy, USP6).
  (iii) The registered same-strand guard counts 9 antisense / fragment loci of the copies among the 13 non-copy NPIP recruits; strand-blind NPIP precision is .967 (BASE) -> .909 (V-maj), 4 loci overlap no copy (PKD1, ACSM1, two repeat-derived stubs). (iv) The joined fraction of a recruit's own exonic bp does not separate true from false recruits (NPIPB5 .22, PKD1P6-NPIPP1 .16, PKD1 .02, stubs 1.0).
- **Recommended frozen rule (words).** Container v2: within each family, classes of exon blocks joined by aligned exon-on-both-sides columns; a class is core iff a strict majority of the family's members carry it; blocks core / accessory / undetermined (coreless family); relations of a non-core block to the CORE classes of other families. Membership stays the MCL partition: no recruiting (NOT), no pruning (a `peripheral` column only). Opt-in next to v1, byte-identical v1 default.
  Register next, on the reserved substrates: the container v2 partner-base measure and the relation channel as a candidate-membership rule, with a null that does not draw from a pool the rule exhausts. **Implementation note:** `figs/container_csm.md` §11 (extend `family_container.rs` `project_paf` to keep the partner block index, `classify` two-pass union-find / majority, columns `class_id support fam_size peripheral`, `--container-core support2|majority`, driver `RUSTLE_FAMILY_CONTAINER=2`, tests from the `test_csm.py` fixture and v1 byte-identity).
- **Deviations:** none after Amendment 1 (no bar, arm, threshold or truth changed; no scorer defect found). `smoke.py` exercised the code paths with a synthetic membership before the runs. **Limits:** dev only; S has 10 fusions and error-free reads; the gorilla NPIP / TBC1D3 families are off the dev contig (census only); Compara is universe-limited (936 of 977 judged chr16 recruits unjudgeable); PRIMARY labels favour recruits (STRICT beside); the null pool is 78% exhausted by the rule.
- **Draft register rows (suffix L, not appended):** 1194L-1200L in `figs/container_csm.md` §12. Prereg sha1 before this Outcome `60caa524`.
