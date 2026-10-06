# Pre-registration: container v2 (majority core) on RESERVED substrates: (a) partner-base accessory calls, (b) the RELATION CHANNEL as a candidate-membership rule

**Written 2026-09-30, before any reserved substrate was scored** (only structural gates G-V1 / G-R1 touched them, §0). KEY=v2ho. Python over EXISTING families products:
no aligner was run and none will be; nothing in `src/`, `tools/`, `bench/` is edited (`bench/family_container.py` e197ccb3 and the frozen CSM instruments
`/mnt/linuxdisk/tmp/rustle_figures/container_csm_frozen/` are IMPORTED unchanged); nothing is committed, staged or pushed. Species are never pooled.
Scratch `/mnt/linuxdisk/tmp/rustle_figures_dev/container_v2_heldout/` (`lib/` instruments, `dev/`, `out/`, `tables/`, `logs/`).
Predecessors: `PREREG_container_csm_2026-09-30.md` (container v2 = CSM majority core: DEV only; its Outcome recommended this test), `CONTAINER_HEADROOM_2026-09-30.md`.

## 0. What was seen before this file (and what was not)

**Seen (dev, as in the CSM record).** Container v2 on human chr16 / chr17, the gorilla simulation S and the gorilla OR6737 NC_073244.2 census (`figs/container_csm.md` §0-§13; register drafts 1194L-1200L):
V-maj recruits 1,076 / 537 / 206 pairs, the post hoc observations (partner bases accessory v1 22.7% -> v2 82.6% on chr16, 20.9% -> 100% on chr17; the relation channel relates exactly NPIPB5 and
PKD1P6-NPIPP1 to NPIP on chr16). **Dev looks made for THIS file (all on chr16, chr17, S f0.5, G = NC_073244.2 products only):** (i) the structure of the family-level core-to-core relation RC2 (§2.3: it is registered NOT);
(ii) a dry run of the new instruments on chr16, chr17 and G (RC / NULL counts, pair judge states, fused units): RC strict = **119 pairs / 32 loci (1.4% of 2,351 loci) on chr16, 51 / 28 (1.2% of 2,430) on chr17, 14 / 10 (0.9% of 1,062) on G**;
RC into NPIP on chr16 = 2 candidates, both truth copies (NPIPB5, PKD1P6-NPIPP1; same-strand non-copy 0; NULL-A draws into NPIP gain no copy in any of 5 seeds, NULL-B seeds gain 1-2 copies: 25, 25, 26, 25, 25 of 26); chr17 RC into the TBC1D3 family = 3 candidates, gains none (MCL62 is a 2-member family:
its members are core in their own family, not accessory), NULL-B seeds gain 0-2 copies; pair judge on chr16 (PRIMARY, same-strand): RC supported 21 of 119 (yield .231, judged precision 1.000) vs NULL-A .135-.187 and NULL-B .221-.244 (NULL-B seeds are at or above RC), STRICT yield .093 vs NULL-A .062-.113 and NULL-B .044-.074;
chr17 yield .087 vs NULL-A .043-.114, NULL-B .083-.175; gorilla G: no pair judged (13 of 14 unjudgeable, 1 own; Liftoff is the only gorilla family truth: 408 labelled genes); part (a) units on dev: chr16 6 fused NPIP members, partner bases accessory v1 17.6% -> v2 97.2%, the copy's own bases accessory 5.5%; chr17 2 members, 20.9% -> 100%, 2.2%.
These numbers are DEV, they fixed the statistics of §6 (see the iteration log, §8.12) and **no dev number enters a verdict.**
**Seen on reserved substrates in earlier work (disclosed exposure).** All reserved products served the F1 / F1v2 / COVER / containment / headroom studies (families, Compara F, cluster sets); the headroom study counted fused loci on the A119b chromosomes
(41 present-unplaced Compara multi-copy genes in fused loci, 24 on chr1 / chr19; 15 loci joining two multi-copy families), read the copy fate tables of the seven cells (`fate.*.D.json`) and the container-v1 partner-base sums of twelve fused copies (chr16, chr17, testis: 846 core / 61 accessory bases), and the CSM G census opened **OR6737 NC_073244.2** (25 whole-contig families: recruit counts, support histograms; nothing scored).
NC_073244.2 is also the gorilla development contig of the O1 work. **Not seen:** any container-v2 class, support, core, RC candidate, NULL draw, fused-unit count or truth-support state on any reserved product.
**What the two structural gates touched:** G-V1 (csm's V1 labels equal the stored container-v1 blocks, class and geometry) and G-R1 (the stored v1 `rel_families` of accessory blocks equal the relations read from the joins) were run on every reserved product and printed only pass / mismatch counts and structural counts (loci, clustered loci, families, joins, PAF accounting, §5).

## 1. Question and scope

Two things looked good on dev without being validated. **(a) The container.** Does the majority-core container call the partner bases of fused family members accessory (v1 calls them core whenever the partner is co-duplicated), without calling the copy's own bases accessory?
**(b) The relation channel.** Is "an accessory block of a clustered locus joined to a core class of another family" a precise candidate-membership signal, better than matched random draws, harmless to the partition, cheap in volume, and does it bring NPIP / TBC1D3 copies into their family?
**Cannot settle:** any default; the CSM recruit rule (registered NOT); families of 1-3 members (v2 = v1 there unless coreless); gorilla family-level truth (none exists beyond Liftoff pairs at identity >= .95 and the two copy truths).

## 2. Definitions (binding; frozen in words before any reserved score)

### 2.1 Container v2 (the CSM majority core, unchanged)
Loci, blocks, joins, classes and support exactly as `PREREG_container_csm_2026-09-30.md` §2 (`csm.py` c3d7141e): a block is joined to a block of another locus iff one aligned CIGAR column (M / = / X) of one PAF record maps an exon base onto an exon base (strand-blind, every record);
inside a family F the classes are the connected components of its blocks under joins between two members; support = distinct members; **a class is CORE iff 2 s > |F|**; a family with no core class is CORELESS (its blocks are UNDETERMINED: no accessory call); a block is ACCESSORY iff its class is not core in a non-coreless family.
Membership of the partition is the MCL partition, unchanged. Container v1 (core iff support >= 2) is the comparator; for a non-coreless family of 2 or 3 members v1 = v2 block for block.

### 2.2 The relation channel RC (the registered lead)
A CLUSTERED locus L of family F is a **candidate member of G != F** iff >= 1 ACCESSORY block of L is joined (column test) to a block lying in a CORE class of G (G must have a core class). The output is a list of (locus, family) candidate pairs; **the partition is unchanged, RC only adds pairs.**
**Unclustered loci: RC does nothing** (no container row, no candidate pair; the rule reads the container of a member). Coreless source families: nothing (no accessory call). **RC-u** (sensitivity, beside): coreless blocks counted as accessory (the dev post hoc REL).
n_G = the number of candidate loci of target family G; a **candidate locus** is a locus in >= 1 candidate pair.

### 2.3 RC2 (family-level core-to-core relation): REGISTERED NOT, before the reserved run
Definition tried: two families F, G where a CORE class of F is joined to a CORE class of G; F's members are candidate members of G. Direction variants: symmetric; subsumption (every core class of F joined to a core class of G); smaller-into-larger. **DEV behaviour (chr16, chr17, G = NC_073244.2, S f0.5):**

| product | families (core-core to another family) | ordered family pairs | candidate (locus, family) pairs: symmetric / subsumption / smaller-into-larger | RC pairs |
|---|---|---|---|---|
| chr16 | 111 (60) | 154 | 794 / 409 / 174 | 119 |
| chr17 | 78 (31) | 64 | 282 / 152 / 79 | 51 |
| G | 25 (12) | 44 | 121 / 37 / 34 | 14 |
| S f0.5 | 35 (22) | 54 | 266 / 63 / 52 | 132 |

It merges NPIP fragment families wholesale (chr16 MCL5, 9 NPIP fragment loci, is subsumed by NPIP MCL1 and by MCL6, MCL64, MCL75) and repeat families wholesale (chr16 MCL2, 28 loci in ONE class, is subsumed by seven families; chr17 MCL0, 12 loci, by five; gorilla ZNF families MCL35 / MCL62, 3 loci each, into five families each);
the motivating pair chr17 MCL62 -> MCL3 holds, but MCL62 is also core-core joined to MCL0, MCL5, MCL16, MCL25, MCL2 and MCL28, and "mutual subsumption" (MCL3 is subsumed by MCL62) means the direction needs a size rule; the only constant-free tie-break (argmax target) does not stop the wholesale merges.
**RC2 is therefore NOT: it is not run, not scored and not offered on any reserved substrate.** (Its dev structure is the evidence; it is a finding about the rule, not a measurement on reserved data.)

### 2.4 Pool and the two matched NULLs (none draws from a pool the rule exhausts)
**POOL(G)** = every CLUSTERED locus L outside G with >= 1 block joined to a block of ANY class of ANY member of G (core or not, accessory or not). RC is a subset of the pool; the census reports RC / pool (dev .12-.14).
**NULL-A (5 seeds, `random.Random(20261000 + k)`, k = 1..5):** for every target family G (sorted), n_G loci drawn uniformly without replacement from POOL(G).
**NULL-B (5 seeds, `random.Random(20262000 + k)`):** random ACCESSORY blocks: each clustered locus of a non-coreless family keeps its NUMBER of accessory blocks but the blocks are redrawn uniformly among its blocks; the RC test (joined to a core class of another family) is applied to the redrawn blocks; per target family min(n_G, available) loci are then drawn uniformly (the shortfall is reported). NULL-B tests whether the ACCESSORY designation matters beyond "joined to a core class".
One rng per (product, seed); families, loci and pools in sorted key order.

## 3. Truth tiers and the pair judge (per species; nothing is pooled across species)

Gene = (contig, Name) of a RefSeq gene / pseudogene / ncRNA_gene record (`families_gw/species/<sp>/genes_only.gff`, exon unions `genes.tsv`, frozen `npf.load_genes_tsv`).
**Tiers** (a tier maps a gene to a family label or to None = outside that truth): **compara** (human only: Ensembl Compara Primates, whole-genome families of >= 2 genes present in the annotation, `SC.load_truth`, first family per gene); **liftoff** (connected components of genes linked by Liftoff copy pairs, sequence_ID >= .95, both exon unions >= 200 bp, no read-support filter: a source gene and every non-readthrough gene covering >= 50% of the copy's exon bases; human 660 labelled genes, gorilla 408);
**copy** (only the target family G* holding the most truth copies of NPIP / TBC1D3: `truth.{hsa,ggo}.json` territories); **names** (REPORTED BESIDE, never in a decision: the alphabetic prefix of the RefSeq symbol before its first digit, >= 3 letters, LOC / LINC / MIR / SNOR excluded: old name-based truths were retracted, and the rule is not constant-free).
**Readings.** STRICT: the locus's own `gene_at` label (the BASE convention). PRIMARY: every gene whose span overlaps a block of L joined to a block of a member of G (any class: the same blocks for RC and for both NULLs). The family's gene set = the `gene_at` labels of every record of its members. Copy tier: PRIMARY = the joined blocks overlap a truth copy (>= 1 bp), STRICT = a strict majority of the exon bases of L lies in truth-copy territory; same strand (registered) or strand-blind (beside).
**Judge states** of a pair (L, G), first match: *unlabelled* (no gene) | *own* (every labelled gene of L is already a gene of G: a same-place locus, neutral, excluded from every fraction) | *supported* / *unsupported* (>= 1 tier judges: both sides labelled in it; supported iff the labels meet in some tier of the decision union compara + liftoff + copy) | *unjudgeable* (new genes exist, no tier has an opinion).
Reported separately and never pooled away: judged, unjudgeable, own, unlabelled, per tier. **Decision statistics** per (arm, reading): **yield = supported / (pairs - own)**, **precision = supported / (supported + unsupported)** (judged). (Dev showed precision saturates at 1.0 for every arm, so yield is the informative statistic; both are required, §6.)

## 4. Part (a): partner-base accessory calls of fused family members

**Unit = a FUSED MEMBER**: a clustered locus L of a family F with |F| >= 4 (for |F| <= 3 and non-coreless, v1 = v2 block for block: reported beside) whose exon blocks overlap (>= 1 exon base, same strand as L when L has one) >= 2 annotated gene records (not readthrough-described) with disjoint spans,
>= 1 COPY-SIDE and >= 1 PARTNER-SIDE. **copy-side** = the record carries the MAJORITY LABEL of F in the species' family-level truth (Compara for human, Liftoff for gorilla: the label carried by >= 2 members and a strict majority of F's LABELLED members), or, for G*, a strict majority of its exon bases lies in a truth copy (same strand); **partner-side** = every other record (another family's label or none: a co-duplicated PKD1P beside NPIPA1 is a partner).
Bases: the exon bases of L's blocks in the merged exon unions of the copy-side records (copy bases) and of the partner-side records minus the copy bases (partner bases); blocks are classed by v1 (core / accessory) and v2 (core / accessory / undetermined); **decision bases = the bases in blocks determined under v2** (the same bases for v1 and v2); units of coreless families are listed and excluded.
Per group: A1 / A2 = partner bases in accessory blocks under v1 / v2; K1 / K2 = copy bases accessory under v1 / v2; B = determined partner bases; Cd = determined copy bases; sign = units whose partner accessory fraction is higher / lower under v2. Beside: partner records with a label other than the family's (positive non-membership), unlabelled-partner share, |F| <= 3 units.
**Evaluable group:** >= 5 units with determined partner bases (the smallest n at which an all-one-direction sign test reaches one-sided p < .05).

## 5. Substrates, exposure, evaluation filter

| group | species | products (arm D = current defaults: f1v2 regroup + `--min-cov-shorter 0.70`) | gates |
|---|---|---|---|
| **A119b** | human | 18 per-chromosome products: chr1-12, 14, 15, 19, X, Y, M (`repfam/chrN/D.fam.*`, families GTF and PAF `rep/chrN/`); chrM has no family | G-V1 and G-R1 PASS on all 18 (blocks 0 mismatches) |
| **testis** | human | genome-wide `defaults_flip/runs/new/human_testis.*` (identical to `repfam/D.fam.*`, cmp) | PASS (4,888 blocks) |
| **OR6737** | gorilla | five contigs NC_073224.2 / 228.2 / 241.2 / 242.2 / 244.2, `fam/c5/D.fam.*` | PASS (3,669 blocks) |
| **KB3781** | gorilla | same five contigs | **FAIL: its `D.fam.loci.paf` is absent from disk; the only KB PAF is the pre-flip genome-wide BASE run: replay gives 62 of 2,324 class mismatches and 2 relation mismatches. Registered NOT SCORED (no aligner run).** A replay arm **KB3781r** is run beside, labelled inexact, never in a verdict |

**Evaluation filter (dev contigs).** Decisions in (a) and (b2) exclude pairs / units whose LOCUS lies on a development contig: testis chr16, chr17, chr18 (human NPIP / TBC1D3 blocks and the development chromosomes); OR6737 NC_073244.2 (O1 development contig, opened by the CSM census); the same filter is applied to RC and every NULL seed alike (draws are made on the whole product, the filter is applied to the pairs). The NPIP / TBC1D3 copy blocks (b3) are scored as blocks on their copy contigs. A119b contains no development chromosome (chr16, chr17 are excluded from the set).
**Truth copies for (b3):** testis NPIP (26 Dishuck, chr16) and TBC1D3 (16 RefSeq, chr17); OR6737 NPIP (25 T_member) and TBC1D3 (14); holder = most same-strand exon bp over the copy (headroom rule); G* = the cluster holding most copies.

## 6. Clauses, evaluability, verdict mapping (integer / exact-fraction comparisons; fixed now)

**Part (a), per evaluable group:** (a1) A2 > A1 (the partner bases in accessory blocks under v2 strictly exceed v1's, on the same determined bases) and sign v2_better > v2_worse; (a2) 10 x K2 <= Cd, where K2 = copy bases in accessory blocks under v2 and Cd = determined copy bases (the copy's own bases accessory <= 10%). A group passes iff (a1) and (a2).
**(a) verdict:** WORKS iff >= 1 human and >= 1 gorilla group evaluable and every evaluable group passes; PARTIAL iff >= 1 evaluable group passes and not WORKS; NOT otherwise (including "no evaluable group": not shown, labelled unpowered).
**Part (b), clauses:**
(b1) harmless: in every group the partition equals the clusters file, every candidate L is outside its target, and **zero truth copies are lost** (asserted; RC only adds).
(b2) selective: a group is evaluable iff RC has >= 10 judged pairs under BOTH readings; an evaluable group holds iff, for EVERY one of the 10 NULL seeds (A1-A5, B1-B5) and in BOTH readings (PRIMARY and STRICT), yield_RC > yield_seed (strict) AND precision_RC >= precision_seed. (b2) holds iff >= 1 human group is evaluable and every evaluable group holds.
(b3) gain: a copy block (group, family) is evaluable iff its REACHABLE headroom >= 1 (truth copies whose holder is clustered outside G* and lies in POOL(G*); unclustered or pool-less holders are unreachable by construction); it holds iff RC copies in G* >= BASE, > every NULL seed's (10 seeds) and same-strand non-copy candidates <= copies gained (strand-blind reading beside). (b3) holds iff >= 1 block is evaluable and every evaluable block holds.
(b5) volume: per group, candidate loci <= 5% of all loci of the product (A119b summed over chromosomes; the maximum chromosome reported). The 5% is half of the CSM falsifier Z5 (10%), not fitted.
(b6) no dev number is used: nothing from `dev/` enters `out/verdict.json`.
**(b) verdict:** WORKS iff b1 and b2 and b3 and b5; PARTIAL iff b1 and b5 and (b2 or b3) and not WORKS; NOT otherwise.
**Expected from the headroom study (stated now):** the headroom oracle gains no copy on the reserved copy cells (testis 5 of 26, OR NPIP 10 of 25, OR TBC1D3 7 of 14 in family; wrong-family copies testis NPIPB4 / NPIPB5 and OR NPIPB7 have no alignment to their family), so (b3) is expected NOT EVALUABLE everywhere and **(b) can then reach PARTIAL at most**; gorilla has no family-level truth beyond Liftoff and the copy truths, so the gorilla groups are expected NOT EVALUABLE for (a) and (b2). The test that can decide (a) and (b2) is A119b (and testis).

## 7. Predictions (prior probabilities) and falsifiers

P1. Structural gates hold on the first attempt for every reserved product except KB3781 (met: 20 of 21; KB3781 failed as stated in §5).
P2 (a). A119b evaluable (>= 5 units) 0.80; testis 0.25; OR6737 0.10. Where evaluable: A2 > A1 0.90 and sign v2 better > worse 0.85 (A119b); copy bases accessory <= 10% 0.70 (dev 5.5% / 2.2%, but sub-clade exons of large families become accessory). Verdict (a): WORKS 0.03, PARTIAL 0.65, NOT 0.32.
P3 (b). Volume <= 5% in every group 0.90. (b2) evaluable in A119b 0.95, testis 0.50, OR6737 0.10; **RC beats all ten NULL seeds on yield in both readings on A119b 0.15 (dev: RC beat NULL-A, not NULL-B)**; precision >= all seeds 0.85. (b3) evaluable anywhere 0.08. Verdict (b): WORKS 0.01, PARTIAL 0.15, NOT 0.84.
P4. RC candidates are dominated by annotation-fused or pseudogene loci (>= 40% of candidate loci) 0.60; >= 25% of RC pairs are same-place joins (overlap a member on the genome) 0.60; >= 25% of RC pairs are repeat-derived (>= 50% of joined bases in interspersed repeats) 0.40.
**Falsifiers.** Z1 a truth copy lost or a candidate inside its own family (instrument defect). Z2 a NULL seed reaches RC's yield (the accessory designation or the core-class join adds nothing over the pool). Z3 volume > 5%. Z4 copy bases accessory > 10% (container over-calls the family's own module). Z5 v2 partner accessory <= v1 (no container gain). Z6 RC pairs concentrated: >= 50% of pairs in the 3 largest target families of a group (a few repeat / SD families carry the rule).

## 8. Hostile self-review

1. **Repeat-derived exons.** Candidate pairs whose joined bases are >= 50% in RepeatMasker interspersed repeats are counted per group (description, not a rule); a repeat-driven result is a finding about RC.
2. **Single-exon stubs.** The container is informative only for multi-block loci; single-block candidates are counted (dev 0 of 32 / 28 / 10).
3. **Segmental-duplication mosaics (NBPF, GOLGA, ANKRD, ZNF).** Class chaining through one shared module can make a repeat-like class the majority core of several families; concentration (top-3 target share) and classes holding two blocks of one locus are reported.
4. **Two-member families.** 50-57% of families are pairs, where "majority" is "both" and v2 = v1; 2-member TARGET families take the largest RC share on dev (36 of 119 pairs); stratified.
5. **The USP6 case.** A real chimeric gene (USP6 beside TBC1D3) is a non-copy under the copy truth and a chimera by biology; annotation-fused candidate loci are counted and their support is reported apart.
6. **Unjudgeable labels.** Compara is universe-limited and Liftoff at identity >= .95 covers recent duplicates only: most pairs are unjudgeable (dev 50-70%); supported / unsupported are reported beside unjudgeable / own / unlabelled in every table, and the gorilla groups may have no judged pair at all.
7. **Same-place joins.** Strand-blind joins between loci that overlap on the genome (antisense or alternative loci) are counted; `own` pairs are neutral in the fractions.
8. **Saturating precision.** Precision among judged pairs ties at 1.0 for nearly every arm on dev; requiring it beside a strict yield comparison keeps the clause from being decided by ties; yield counts unjudgeable pairs as unsupported and can reflect annotation coverage, not quality (both arms pay the same price).
9. **NULL-B shortfall.** The redraw may offer fewer loci than RC's n_G (dev 6-15 of 119 / 51); fractions, not counts, are compared and the shortfall is printed.
10. **Species and samples.** A119b and testis are different samples of one species and share the annotation and the Compara truth; OR6737 is one gorilla sample (KB3781 is not scored); no pooling across species; the A119b chromosomes are summed within the sample only.
11. **Dev exposure of the truths.** NPIP / TBC1D3 truths were developed on chr16 / chr17; the testis NPIP cell uses the same truth copies on a different sample (and is filtered from the (b2) decision); OR6737 NC_073244.2 is census-exposed and filtered.
12. **Iteration log of the dev-made choices (honest list).** The copy-side rule of (a) went through four dev versions before this file: (1) "label shared with another member" (counts co-duplicated partners as copy: rejected on chr16, copy bases accessory 36%); (2) plurality label per tier (Liftoff's PKD1P component outvoted NPIP: rejected); (3) strict majority of ALL members (NPIP MCL1 has exactly 15 of 30 labelled: no label: rejected as knife-edge); (4) the final rule, a label carried by >= 2 members and a strict majority of the LABELLED members, one tier per species. The selectivity statistics went from precision alone (saturated) to yield plus precision. The 5-unit and 10-pair floors and the 5% bound are round numbers, not fitted.

## 9. Order, machine rules, files

1. This file. 2. **Amendment 1**: instrument sha1 (`lib/`), unit tests, mutants, G-N (null determinism and counts), independent numpy check of RC pairs and pools on the dev products, the structural gates on every reserved product, every ambiguity resolved, all BEFORE any reserved score. 3. Runs: A119b (all 18), testis, OR6737, then KB3781r (beside), then `verdict` and `tables`; the first outputs are final. 4. **Outcome** appended here; register rows with suffix N drafted in the report, not appended; report to `figs/container_v2_heldout.md`.
Light jobs only (`tools/rlock.sh light`; every product is seconds and < 0.3 GB), foreground, TMPDIR and scratch under `/mnt/linuxdisk`, never `pkill -f`; another agent holds the heavy lock at times and nothing here needs it.
A defect found after the first reserved score is reported as a deviation with both numbers, never by changing a clause, statistic, seed or filter.

## Amendments

**Amendment 1 (2026-09-30, about 15:40, after the instruments and the gates, before ANY reserved outcome was computed or printed; the frozen `SHA1SUMS` is stamped 15:39, the first reserved product file 15:41).** Frozen instruments (read-only copies, `sha1sum -c` passes):
`/mnt/linuxdisk/tmp/rustle_figures/container_v2_heldout_frozen/` (`SHA1SUMS`, `SHA1SUMS.imported`, `SHA1SUMS.inputs`, `SHA1SUMS.products`); scratch `/mnt/linuxdisk/tmp/rustle_figures_dev/container_v2_heldout/` (`lib/` = the same bytes). Prereg sha1 before this amendment `bfa37dbe`.
Instrument sha1 (first 8): `csm.py` c3d7141e, `common.py` 4e1a8f8b, `test_csm.py` 8fb56971 (the three are the CSM instruments, unchanged), `rc.py` cd33eeac, `truths.py` c75c2c05, `fusedm.py` d402dbc5, `product_run.py` 401d9c57, `verdict_v2.py` c3778058,
`gate_only.py` a07400e2, `gate_n_rc.py` 9359a877, `verify_rc.py` 038753c4, `test_rc.py` 9366574f, `test_verdict.py` 761215d3, `mutants_v2.py` 158a1431, `specs_reserved.json` 16ee1944, `specs_dev.json` 4d1862d2, `explore_rc2.py` d885884e.
Imported unchanged: `bench/family_container.py` e197ccb3, `tools/rlock.sh` 30f424a9, `o1_cover_frozen/score.py` 25242385, `npf.py` fdc1a7d3, `figures/_liftoff.py` 952fc546, `figures/figlib.py` 3bfd0414, `container_headroom/code/hl_core.py` 1bc19291 and `oracle.py` 6975625e.
Annotation and truth inputs (sha1 first 8 in `SHA1SUMS.inputs`): human / gorilla `genes_only.gff` 3274ab7a / e6c79d53 and `genes.tsv` 8c55e66a / 04045ac0, Compara Primates 8e8affdd, Liftoff loci human / gorilla ba639599 / d48bb51a, truth copies human / gorilla 644f89cc / 3a7d6a5d, RepeatMasker hs1 0df06713 and GCF_029281585.2 0eeef79b. The 126 product files (GTF, PAF, clusters, loci, container) are listed with sha1 and size in `SHA1SUMS.products`.

*Gates, all passed before any outcome (outputs in `out/gates_g0_g1.txt`, `out/gates/*.gate.json`, `out/verify_rc.*.json`):*
- **G0** sha1 of `family_container.py` and `rlock.sh`; **G1** `test_csm.py` 11 tests, `test_rc.py` 22 tests (RC, pools, both NULLs, the pair judge, the fused-unit finder, interval helpers), `test_verdict.py` 10 tests (every clause, floor and verdict branch) pass from the frozen copies; **25 of 25 mutants** of `rc.py`, `fusedm.py`, `truths.py`, `verdict_v2.py` are caught (one more, dropping the redundant guard `g != fid` of `rc_pairs`, is EQUIVALENT: an accessory block joined to a block of its own family lies in the same class, so it cannot be in a core class it is not in).
- **G-N** (dev products chr16, chr17, G) the RC / NULL-A / NULL-B sets are deterministic, every pair is clustered and outside its family, every arm pair lies in its target's pool, NULL-A has RC's per-target counts, NULL-B at most RC's and offered + shortfall = RC's count.
- **G-I independent verification** (`verify_rc.py`, numpy per-column expansion and union-find, no import of `csm.py`, `rc.py` or `family_container.py`): joins, RC pairs, RC-u pairs and pools are IDENTICAL to `rc.py` on chr16 (9,058 joins; RC 119, RC-u 120, pool 965 pairs), chr17 (5,648; 51, 51, 365), G (1,448; 14, 14, 113), S f0.5 (5,000; 132, 157, 648); every NULL-A draw lies in the independent pool and every NULL-B draw within the independent offer bound.
- **G-C** the V-maj recruit recomputation in `product_run.py` equals the CSM census (1,076 / 537 / 206 pairs) and RC is a subset of the clustered V-maj recruits (119 of 710, 51 of 263, 14 of 111).
- **G-V1 and G-R1 on the reserved products (frozen copies, `gate_only.py`):** PASS on the 18 A119b chromosomes (chr1 2,719 blocks ... chrY 1,187; chrM no family), testis (4,888 blocks, 342 families, 18,304 joins) and OR6737 (3,669 blocks, 136 families, 13,843 joins), 0 class / geometry / relation mismatches everywhere. **KB3781 and KB3781r FAIL** (2,324 blocks: 62 class mismatches, 2 relation mismatches): KB3781 is **NOT SCORED** as registered (§5); KB3781r (the same replay, `beside_inexact`) is run and reported beside, never in a verdict.
- A dev dry run of the final driver reproduces the pre-refactor numbers byte for byte (chr16).

*Implementation ambiguities resolved (binding from here on):*
1. RC is computed from the joins directly (`rc.rc_pairs`: accessory = block not in a core class of a non-coreless family; target block in a core class of another family), which equals the container-v2 relation rows of accessory blocks of non-coreless families (the independent numpy check agrees); `joined_blocks_any` (blocks of L joined to a block of any member of G) defines the PRIMARY blocks for RC and both NULLs alike.
2. `own` pairs (every labelled gene of L is already a gene of G) and `unlabelled` pairs are separate states; yield = supported / (pairs - own), precision = supported / (supported + unsupported); unlabelled and unjudgeable pairs stay in the yield denominator; the decision uses the same-strand copy tier (strand-blind reading beside).
3. The copy tier's strand is `Locus.strand` (the `loci.gff3` gene strand of the representative key); the holder rule and G* are the headroom rule (`hl_core.holder_of`, cluster holding most copies).
4. Evaluation filter by the LOCUS contig (dev contigs: testis chr16 / chr17 / chr18, OR6737 and KB3781r NC_073244.2); applied to pairs and units of RC and every NULL; the copy blocks (b3) and the census are unfiltered; volume (b5) uses all loci and all candidate loci of the product.
5. Part (a): |F| >= 4 decides (|F| <= 3 reported), coreless-family units are listed and excluded, the family label tier is Compara for human and Liftoff for gorilla, `family tau` needs >= 2 carrying members AND a strict majority of the LABELLED members; copy-tier copy-side needs a strict majority of the record's exon bases in a truth copy; overlap with the locus blocks is on exon bases of records on the locus strand.
6. A119b is decided as ONE group (the 18 per-chromosome products summed, counts only); testis, OR6737 and KB3781r are separate groups; nothing is pooled across species or samples.
7. NULL seeds 20261001-5 (A) and 20262001-5 (B), one `random.Random(seed)` per (product, seed); the pair tables `tables/rc_pairs/<product>.tsv` list every RC candidate with its judge states (hostile review).
8. The three tiers of the decision union are compara (human), liftoff, copy; `names` is reported beside per tier only.

*Order.* `product_run.py` on A119b_chr1-4, chr5-9, chr10-15, chr19 / X / Y / M, then testis, OR6737, KB3781r (light, foreground, each a few seconds to a minute), then `verdict_v2.py`; the first outputs are final. A defect found after the first scoring will be reported as a deviation with both numbers, never by changing a clause, statistic, seed or filter.

## Outcome (2026-09-30, about 15:50; the products ran 15:41-15:45, `verdict.json` 15:45; scored with the Amendment 1 definitions, no deviation after it; full report scratchpad `figs/container_v2_heldout.md`, data `/mnt/linuxdisk/tmp/rustle_figures_dev/container_v2_heldout/` `out/verdict.json`, `tables/verdict.md`, `tables/rc_pairs/`)

**Registered verdicts (§6, exact-fraction clauses): (a) NOT, (b) NOT; RC2 NOT (registered before the run); KB3781 not scored (gate failure, as registered).** Every product was run once, the first outputs are final.

| measure (reserved; species and samples separate) | A119b (18 chromosomes, human) | testis (human) | OR6737 (gorilla) | KB3781r (replay, beside) |
|---|---|---|---|---|
| G-V1 / G-R1 | PASS x 18 | PASS | PASS | FAIL (62 of 2,324 class, 2 relation mismatches) |
| fused members, |F| >= 4 (decision / with partner bases) | 7 / 7 | 1 / 1 | 1 / 1 | 0 |
| partner bases accessory v1 -> v2 (bp of determined bases) | 480 -> 757 of 7,990 (6.0% -> 9.5%) | 0 -> 0 of 9 | 885 -> 885 of 885 | - |
| units with the partner fraction higher / lower / equal under v2 | 2 / 0 / 5 | 0 / 0 / 1 | 0 / 0 / 1 | - |
| copy's own bases accessory under v2 | **4,459 of 35,166 (12.7%)** (v1 621, 1.8%) | 0 of 1,535 | 491 of 2,622 | - |
| (a) clause | evaluable; a1 holds, **a2 fails** | not evaluable | not evaluable | not evaluable |
| RC candidate pairs / loci (% of all loci) / RC-u pairs | 788 / 323 (0.50%) / 960 | 18 / 16 (0.13%) / 19 | 45 / 36 (0.69%) / 47 | 26 / 20 (0.43%) / 29 |
| pair judge, PRIMARY: supported / unsupported / unjudgeable / own / unlabelled | 68 / 33 / 516 / 96 / 75 | 0 / 5 / 9 / 2 / 1 | 0 / 0 / 23 / 2 / 1 | 0 / 0 / 14 / 0 / 0 |
| yield RC vs NULL-A (5) vs NULL-B (5), PRIMARY | **.098** vs .067-.080 vs .074-**.100** | not evaluable (5 judged) | not evaluable | - |
| yield RC vs NULL-A vs NULL-B, STRICT | .082 vs .052-.062 vs .056-.074 | not evaluable | not evaluable | - |
| precision RC vs max seed (PRIMARY / STRICT) | .673 vs .659 / .627 vs .594 | - | - | - |
| (b2) | evaluable; **fails: NULL-B3 yield 59/588 = .1003 >= RC 68/692 = .0983 (PRIMARY)**; the other 19 yield and all 20 precision comparisons hold | not evaluable | not evaluable | - |
| (b3) reachable headroom / RC candidates into the copy family | - | NPIP 0 of 2 / 0 | NPIP 0 of 1 / 0; TBC1D3 0 of 0 / 0 | NPIP 0 / 0 |
| (b1) copies lost / candidates inside their family | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| (b5) volume <= 5% | yes (max chromosome chr15 1.5%) | yes | yes | yes |

- **Clauses.** (a): only A119b is evaluable (7 units >= 5); the pooled partner accessory bases rise (757 > 480, sign 2 / 0) so (a1) holds, but the copy's own bases are 12.7% accessory under v2 (> 10%): (a2) fails, so the group and (a) fail; no gorilla group is evaluable (1 unit OR6737, 0 KB3781r). (b): b1 and b5 hold everywhere; b2 is evaluable only on A119b (judged 101 PRIMARY / 83 STRICT) and fails on one of 20 yield comparisons by .002; b3 is not evaluable on any block (every wrong-family copy is outside POOL(G*): testis NPIPB4 / NPIPB5, OR6737 NPIPB7; TBC1D3 none, KB3781r none), exactly as the headroom study predicted; so (b) = NOT although RC beats all five NULL-A seeds on yield in both readings and all ten seeds in STRICT.
- **Predictions (§7).** Hit: P1, P2 (A119b evaluable, A2 > A1, sign better > worse; testis and OR6737 not evaluable), (a) NOT (p .32), P3 volume, (b2) evaluable on A119b, RC NOT beating all ten seeds (p .85 against), precision >= all seeds, (b3) not evaluable, (b) NOT (p .84), P4 same-place >= 25% (25.5%), repeat-derived >= 25% (40.2%), fused-or-pseudogene loci >= 40% (49.8%). Miss: copy bases accessory <= 10% (12.7%, p .70), (a) PARTIAL favoured (p .65), testis (b2) evaluable (p .50: 5 judged).
- **Falsifiers fired.** Z2 (a NULL seed reaches RC's yield: NULL-B3, PRIMARY), Z4 (copy bases accessory 12.7% > 10%). Not fired: Z1, Z3, Z5, Z6 (group-level top-3 target families carry 4.9% of the pairs; six chromosomes carry 78%).
- **Post hoc (labelled, not registered).** (i) One locus decides (a2): NBPF20 (chr1, MCL10, 9 members) has 3,666 of its 20,101 determined copy bases (18%) in support-2 classes that v2 calls accessory; without it the copy bases are 793 of 15,065 (5.3%) and the partner fraction 1.4% -> 8.0%. (ii) Only 2 of the 7 units change class at all (chr15 GOLGA8N / ARHGAP11A: 110 bp core -> accessory in a 30-member family; chr19 ZNF69;ZNF763 / ZNF700: 167 bp); the other five have the same classes under v1 and v2; units of families of 2-3 members (12, non-coreless 11) are identical under v1 and v2 (341 accessory / 23,935 core bases each), as constructed. The dev effect (NPIP: partner bases accessory 17.6% -> 97.2%, 6 units in one 30-member family) is the extreme case of one co-duplicated partner family; the reserved chromosomes hold no such family except GOLGA8 / ARHGAP11.
  (iii) RC is tiny and SD-dominated: 788 pairs = 4.5% of the clustered V-maj recruit pairs (17,512) and 3.9% of the pool (20,146); 78% of the pairs sit on chr15, chr1, chr5, chr9, chr7, chr2; 29.7% of the candidate loci are annotation-fused, 37.2% overlap a pseudogene (union 49.8%), 25.5% of the pairs are same-place joins, 40.2% repeat-derived (judged precision .625 vs .688), 0 candidate loci are single-block; pair parents named GUSBP 66, GOLGA 60, ZNF 50, ANKRD 39, NBPF 34; 65% of the pairs are unjudgeable (PRIMARY). Judged precision by target size: 2 members .897 (26 of 29), 3 .647, 4-5 .538, >= 6 .520; fused .675 vs not fused .672 (the USP6-like chimeric loci are not worse); by chromosome the judged pairs are chr1 21 / 2, chr2 10 / 2, chr7 11 / 0, chr9 3 / 3, chr10 0 / 3, chr15 6 / 13, chr19 17 / 8, chrX 0 / 2 (supported / unsupported).
  (iv) `names` tier (beside): RC 212 supported / 179 unsupported vs NULL-A 134-160 / 128-152 and NULL-B 173-181 / 123-143: RC is not better than NULL-B on the weak names tier. (v) Compara tier alone: RC 64 / 28 vs NULL-A 43-53 / 24-33 and NULL-B 41-55 / 26-31; Liftoff tier: RC 4 / 6 vs NULL-A 0-2 / 2-8 and NULL-B 3-4 / 5-6. (vi) RC-u (coreless blocks counted) adds 172 pairs on A119b (960 vs 788).
- **Deviations:** none after Amendment 1 (no clause, statistic, seed, filter or arm changed; no instrument defect found). Machine use: every step was a light job; §9 expected < 0.3 GB per product, the testis / OR6737 / KB3781r batch peaked at 1.05 GB (46 s), the A119b batches at 0.32 GB (56-84 s). KB3781 was not scored as registered; KB3781r ran beside after its gates failed, labelled inexact. **Limits:** one human sample for the decisive groups (A119b; testis unpowered); gorilla has no family-level truth beyond Liftoff (408 genes) and the two copy truths, so both gorilla groups are unevaluable for (a) and (b2); (a) rests on 7 units; the headroom of the copy clause is zero on every reserved substrate; KB3781 needs one heavy minimap2 run to be scored.
- **Draft register rows (suffix N, not appended):** 1194N-1200N in `figs/container_v2_heldout.md` §12. Prereg sha1 before this Outcome `7e3bd986`.
