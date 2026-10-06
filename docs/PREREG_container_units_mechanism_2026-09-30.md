# Pre-registration: a UNIT split of fusion transcripts for the families input only — the ceiling of the execution under an ORACLE junction

**Written 2026-09-30, before any product of this test exists.** Dev only (gorilla fusion simulation; human A119b chr16).
Python prototype in scratch `/mnt/linuxdisk/tmp/rustle_figures_dev/container_units_mech/`; nothing in `src/`, `tools/`,
`bench/`; nothing committed or pushed. Species are never pooled. No chimp (PTR) or orangutan (PPY) product is opened.

## 0. What was seen before this file (and what was not)

Seen: the Outcome of `PREREG_fusion_container_sim_2026-09-28.md` (frozen pre-flip binaries: fused copies in NPIP
M1 = 9 / 1 / 0 / 0 / 4 at f = 0 / .1 / .5 / .9 / 1; partners in NPIP M2 = 0 / 1 / 0 / 0 / 4; every fusion assembled with its
exact intron chain from 3 reads up; the copy, the fusion and the partner come out as ONE gene_id); the dev tables of
`PREREG_o1_cover_growth_2026-09-29.md` §12 (Python F1v2 on the same simulation: fused copies in NPIP BASE / F1v2 =
9/9, 7/9, 1/5, 0/0, 1/1; human chr16 NPIP copies in NPIP 22 BASE, 24 ALL/F1v2, Compara F .615 -> .667); register r846 (a
node cut DOUBLES rather than separates, short pieces become hubs), r1013-r1018 (the node-replacement oracle of r1017: referee
F +0.021, NPIP sensitivity .600 -> .550, "the pieces fail the coverage gate that the fused locus passed"), r1121-r1127, r1145-r1151,
1165D-1169D; the source of `bridge_regroup.rs` (`components`, `rg3_pieces`, `regroup`), `mcl_families.rs` (`gtf_loci`, the
deferred edge rule), `tools/rustle_pipeline.sh`, the frozen scorers `score.py` (sim, sha1 8f512fb8) and `o1_cover_frozen/*`.
**Not seen:** any product of the e163d955 build on these BAMs, any unit-split product, any score of this test. The numbers
above were measured with OTHER binaries (pre-flip, frozen, or the Python F1v2); every BASE number of this test is re-measured.

## 1. Question and scope

If fusion transcripts were split into UNITS for the FAMILIES INPUT ONLY (the assembled GTF and its chains stay as they
are), can the existing families stage use the units productively, when the split junction is KNOWN (ORACLE)? This
isolates EXECUTION from DETECTION: a real detector would add its own errors; if even the oracle cannot make the
execution productive, no detector can, and if it can, a detector is worth a separate pre-registration. Nothing here can
justify a default: both substrates are dev and every oracle reads truth. The mechanism it isolates is the one named by the
session's results: a fused locus is placed by its most-read representative, and r1017/r1018 showed that replacing nodes
by per-gene nodes cost NPIP sensitivity because the pieces failed a gate the fused locus passed.

## 2. The execution under test ("unit split"; binding)

Input: a GTF (the driver's `PREFIX.gtf`, or `PREFIX.families.gtf`), and a set J(T) of ORACLE junctions per transcript T
(§3; introns of T as 1-based closed (s, e)). For a transcript T with exons e_1..e_n in TRANSCRIPTION order (genomic order
reversed on `-`) and m = |J(T)| >= 1 junctions:

1. **Units.** T is REPLACED in the families input by m + 1 unit transcripts: U_1 = e_1..e_{k_1}, U_2 = e_{k_1+1}..e_{k_2}, ...,
   cut at the oracle introns, in transcription order. Same contig and strand. Ids `<T>.U1 .. <T>.U<m+1>` (never a trailing
   `.<digits>`, which `npf.base_tid` would strip). Each unit keeps every attribute of T except `transcript_id`, `gene_id`
   (set by step 3) and `exon_number` (renumbered 1..), gains `fusion_of "<T>"`, `fusion_unit "<i>/<m+1>"` and
   `fusion_junction "<chrom>:<s>-<e>:<strand>,..."` (all cuts of T). Transcript line start/end = the unit's first/last exon.
2. **Reads of units** = the `reads` of T, for every unit (given). No other attribute is changed.
3. **Regroup, exactly as RG3/F1 do.** Every gene_id g of the input (not only those holding an oracle transcript) is re-derived
   from its surviving transcripts INCLUDING the units: `components` = same-contig, same-strand exon-overlap components
   (>= 1 shared exonic base; NOT junction sharing: that is what `bridge_regroup.rs::components` does, and this
   document follows the code), `rg3_pieces` names them (the piece whose representative max(reads, span, -line index) is best
   keeps `g`; the others become `<g>.rg<k>`, k = 2.. in the order of their representatives' line index). Units take T's
   line position (U_1, U_2, ... in place of T), so line-index ties resolve as they would for T. A Python port
   (`lib/units.py`) of `components` / `rg3_pieces`; **gate G1** (§5) proves the port.
4. **Fused LOCUS** = the input gene_id g(T) of T. Its **unit loci** = the new gene_ids of T's units. The families stage
   clusters the new loci (`mcl_families --from-gtf`, the driver's command, defaults incl. `--min-cov-shorter 0.70`, one locus per
   new gene_id, representative = its most-read transcript); families stay a STRICT PARTITION of unit loci. Each fused locus
   then **inherits the families of its unit loci**: fam(L) = the union over the units of its oracle transcripts. A
   **relation** is one record per oracle transcript T: (T, g(T), the unit loci in transcription order, their families).
5. **Consolidation ("one member per locus per family")** is a COUNTING rule applied after clustering, never inside it:
   a locus L is one member of family F iff >= 1 of its unit loci is in F; a locus whose units sit in the same family is
   counted once there. It changes no family; its effect (double-counted loci, family sizes by unit vs by locus, NPIP
   precision by locus) is reported (§9).
6. **Design choices (binding, made before any product):**
   (a) a unit with a single exon is KEPT as an ordinary single-exon transcript (no filter, no merge into a neighbour); its
   count, exonic length and degree are reported (§9); (b) a transcript with several oracle junctions gives m + 1 units;
   (c) a transcript without an oracle junction is untouched (only its gene_id may change by step 3); (d) an oracle
   transcript is never kept whole; (e) gene_ids of the input that hold no oracle transcript are regrouped by step 3 too, so
   the arms differ from BASE in regrouping as well as in splitting, which is why arm A0 (§4) exists.
7. **The families stage** is the driver's `families` stage, `tools/rustle_pipeline.sh families` (sha1 2c431091), verbatim:
   `mcl_families --from-gtf <families input> --fasta <genome> --threads 4 --min-exonic-bp 1 --min-shared-exon-frac 0.60
   --emit-units --out PREFIX.fam` with the e163d955 defaults (`--min-cov-shorter 0.70`); the aligner is the one inside
   `mcl_families` (`minimap2 -x asm20 -c -X -N 50 -p 0.1 --secondary=yes`, minimap2 2.30-r1287 on PATH, `RUSTLE_MINIMAP2`
   unset); its command line is read from each arm's own `PREFIX.cache/paf/*/key.tsv` and must be IDENTICAL across arms
   (register r1018's lesson: an arm is not comparable until its aligner invocation is copied from the baseline's own log).
   A units input is given to the driver as `PREFIX.gtf` with `RUSTLE_BRIDGE_REGROUP=off` (FAM_GTF = `PREFIX.gtf`).
   Each arm is its own all-vs-all (r1126: minimap2's output depends on the whole target set, so no PAF is reused across arms).

## 3. The oracles

**S (gorilla simulation).** J(T) = {the simulation's fusion intron `intron` of a FUSED pair of `ann/fusions.json`} (10 fused
pairs, sha1 9ba2fdf5; every fusion is assembled with its exact intron chain). A transcript carries the oracle junction iff
one of its introns equals (i0, i1) exactly on the pair's contig and strand; if a pair has no exact carrier in an arm with
fusion reads, the scorer's +-10 bp carrier rule is used for that pair and the fact is reported. The oracle set is defined on
the transcripts of the DEFAULT GTF and applied to whichever of them an arm's input contains. The 10 unfused pairs have no
fusion reads: no oracle junction.

**H (human A119b chr16).** Annotated genes = the RefSeq CHM13 full GFF's gene / pseudogene / ncRNA_gene records (the
species `genes.tsv` + `genes_only.gff` of `families_gw/species/human`, built from `/mnt/linuxdisk/tmp/regress/chm13.gff`,
1,689,421,265 bytes = `Reference/chm13v2.0_RefSeq_full.gff.gz` decompressed; never `HSA_genomic.gff`), exon union per
record, EXCLUDING records whose description contains `readthrough` (they are the fused models themselves); plus the 26
Dishuck copy records of `audit.human_truth()` (which supply exons for NPIPB14P, which has none in RefSeq, and the
NPIPP1 half of PKD1P6-NPIPP1), a same-named record replaced by the copy record. For each SPLICED transcript T (>= 2
exons) of the default GTF: G(T) = the genes on T's strand whose exon union shares >= 1 base with an exon of T; E_g =
the exon indices of T overlapping g; the intervals [min E_g, max E_g] over G(T) are merged where they overlap (a gene
nested in another's index interval, or two models over the same exons, never cut); between two consecutive merged
blocks the oracle junction is the WIDEST intron between the last exon of the upstream block and the first of the
downstream block (ties: the first). T is an oracle transcript iff it has >= 1 such cut (>= 2 blocks). Same contig
filter: chr16.

## 4. Arms and substrates

| arm | families input | what it isolates |
|---|---|---|
| **BASE** | the driver's `PREFIX.families.gtf` of the DEFAULT `assemble` (`--bridge-regroup f1v2`) = F1v2 alone | the current default |
| **A0** | plain `PREFIX.gtf` (`RUSTLE_BRIDGE_REGROUP=off`) regrouped by §2.3 with NO split | the regroup alone (RG3-like; natural bridges kept) |
| **A** (UNITS-A) | plain `PREFIX.gtf` + oracle unit split + §2.3 | oracle units on the plain GTF |
| **B** (UNITS-B) | `PREFIX.families.gtf` (F1v2's output: bridges already removed as relations) + oracle unit split of the oracle transcripts still present + §2.3 | oracle units on top of F1v2 |
| **N** (NULL) | as A, but each oracle transcript is cut at m uniformly random introns (seeded `stable_seed(20260930, tid)`; S: 1 draw, H: 3 draws) instead of its oracle introns | whether knowing the junction matters, or only shrinking nodes (r846's hub mechanism) |

F1v2 alone = BASE (the default). At f = 0 no transcript carries an oracle junction, so A = N = A0 there (asserted by
identical input GTFs; A and N are not re-run). S: f in {0, .1, .5, .9, 1} (nested subset BAMs `aln/f*.bam`, the simulation's own
pool, sha1s in §11); H: A119b chr16 (`A119b.t2t.bam`, `--region chr16:0-96330374`, the stored genome-wide best-AS table,
as `fusion_container_sim/human/run_hsa16.sh`).
Assembly = the driver's `assemble` command (strict junctions, the shipped polish, `--gtf-tpm`, secondary seeding at 0.98
with the genome-wide AS table; H: `--region` instead of `--genome-wide`, the only change), with `--bridge-regroup` default
(BASE) or `off` (A0/A/N source). Binaries: copies in scratch `bin/` of `/mnt/linuxdisk/home/juanfraitu/rustle_target/release`
(main e163d955; sha1s §11).

## 5. Gates (a failed gate stops that substrate; nothing is scored through it)

- **G0** the copied binaries' help states `--bridge-regroup` default f1v2 and `--min-cov-shorter` default 0.70.
- **G1** (regroup port) on the plain GTF, A0's gene_id of every transcript whose gene_id holds no F1v2 bridge junction equals
  the gene_id that `--bridge-regroup f1v2` wrote for it, on every S arm and on H (RG3's names are F1v2's for genes without a bridge).
- **G2** (S scorer fidelity) `lib/score_s.py` on the ORIGINAL frozen-binary run dirs `fusion_container_sim/runs/f*/` reproduces
  `score/score.json` (sha1 44459a90) for M1, M2 and the container M3 precision/recall at all five arms.
- **G3** (H machinery fidelity) `lib/score_h.py` on the dev BASE products (`o1_cover/dev/hsa16/BASE.*`, `hdev/hsa16.BASE.gtf`)
  reproduces `score.compara.json`'s BASE core (sens .444, prec 1.0, F .615, pairs 66/66) and `npip_dev.json`'s BASE row (22 copies in NPIP, 8 fused copies in NPIP).
- **G4** (assembly) on H: the new `off` GTF vs `hdev/hsa16.BASE.gtf` and the default GTF vs `hdev/hsa16.F1v2.gtf`: cmp result reported (byte-identical expected; a difference is reported, not fatal).
- **G5** (determinism) re-running `mcl_families` on an arm's input with the PAF replayed from its cache reproduces `clusters.tsv` byte for byte (BASE and A of every substrate).
- **G6** (every unit input) the units GTF parses as a GTF, every transcript id is unique, exons are strictly ordered and non-overlapping, the set of exon intervals of the units of T equals T's exons, and away from oracle transcripts the exon sets and `reads` equal the source's.

## 6. Measurements

### 6.1 S (gorilla; per arm BASE, A0, A, B, N; per f)

Scorer logic copied from `score.py` (sha1 8f512fb8) and gated by G2; holders, NPIP family and placement from `npf_audit/audit.py
run()` (config override as `ggo_npip_sim`), the loci of the ARM'S families input and its `clusters.tsv` / `loci.tsv` fold map.
- **M1** fused copies (of the 10) placed in the NPIP family (the cluster holding the most T_member copies). Also the copies in NPIP /25 and the unfused copies in NPIP /15 (**M1u**).
- **M2** S partners (of all 20) whose HOLDER (the locus with the most same-strand exon bp over the partner's all-transcript exon union, ties by reads) is in the NPIP family; also the count over the 10 fused partners.
- **M3u (relation precision / recall).** Reference families are taken from the BASE arm at f = 0 of THIS build: F*_c(pair) = the cluster of the copy's holder, F*_p(pair) = the cluster of the partner's holder (empty if unclustered). Family identity across arms is by
  co-membership, threshold-free: a cluster C (of arm X) "matches" F* iff C has a member locus, other than the pair's own copy and partner loci, whose same-strand exons share >= 1 base with a member locus of F* (other than the pair's own); an empty F* matches an unclustered unit only. For each oracle transcript T of a fused pair, the copy unit = the unit whose exons overlap `copy_fusion_ex` most, the partner unit likewise with `partner_fusion_ex`; the relation record is CORRECT iff the copy unit's cluster matches F*_c AND the partner unit's cluster matches F*_p (and is not the NPIP family unless F*_p is). **Precision** = correct relation records / emitted relation records (oracle transcripts of fused pairs present in the arm's input; one record per T); **recall** = fused pairs with >= 1 correct record / 10. BASE, A0 and any arm where a fused pair has no unit split emit no record for it: recall counts it as missed, precision's denominator does not contain it (reported as emitted / 10).
- **M3v1** the v1 container's M3 (accessory precision and recall against the partner's bases, `score.py`'s definitions, frozen `family_container.py` sha1 e197ccb3) for EVERY arm; for units arms the "fused member" rule of Amendment 2 §C falls to the copy's holder unit (no unit carries the fusion junction), so A = 0 and recall = 0 by construction: reported as such.
- **M4u** per fused pair: copy side / partner side outcomes (SAME / LEAK-into-NPIP / OTHER / UNCLUSTERED) for the relation; the partner "absorbed" flag of BASE (the partner's holder IS the fused locus).
- Identity counters: clusters unchanged vs BASE (by member-key sets, for clusters touching no oracle locus), number of clusters, largest cluster, members.

### 6.2 H (human chr16; per arm BASE, A0, A, B, N1-N3)

- **NPIP (26 Dishuck copies, `audit.human_truth()`; `npip_dev.py`'s holder rule: the locus with the most same-strand exon bp over the copy, ties by reads; NPIP family = the cluster holding most copies):** c1 = copies placed in NPIP (sensitivity = c1 / 26); distinct holder loci (collisions); **NPIP precision** = NPIP-family member loci whose exons overlap >= 1 truth copy / member loci; **partners dragged in** = NPIP-family member loci overlapping no truth copy, listed with `gene_at` labels; bipartite F of (copies, NPIP family) = the harmonic mean of c1 / 26 and the precision (one truth family: the optimal one-to-one match is the NPIP cluster); by unit and by consolidated locus; the per-copy fate table of the 11 copies in fused loci (`FUSED_AUDIT` of `npip_dev.py`: NPIPB2, NPIPA1, NPIPA6, NPIPA7, NPIPA9, NPIPB3, LOC128966608, NPIPB4, NPIPB5, NPIPB6, NPIPB14P): holder unit, its representative (tid, reads, exons), family, in NPIP?, and for a copy not recovered the reason (§9).
- **Genome-wide chr16 (family level; `score.py`'s gene space, gated equal to the binary by its own G2):** Compara Primates (`--chrom ALL` semantics, contigs other than chr16 dropped): bipartite sens / prec / F (scipy assignment; tie policy named, pairwise counts beside), pairwise TP / predicted; Liftoff copy-pair recall (`lo_score.py`, sample human_A119b); and, for the like-for-like comparison with r1017, the protein-homology referee (`fig7/current/human_chr16_ref.families.tsv`, `--chrom chr16`) and Soto / U2 NPIP families (`family_score --family NPIP`, chr16), by the frozen `family_score` (sha1 7723029b), also on r1017's own stored cluster files (`npf_audit/r1017/base_rerun` and `oracle_a`, a declared re-derivation of r1017's oracle).
- **Correct members lost** (for the no-loss bar): a correct member of a BASE cluster = a gene labelled by `gene_at` that shares a Compara truth family with >= 1 other gene of that cluster; it is LOST in an arm iff no arm cluster holds it together with any of those truth-mates. Counted over chr16 Compara and, separately, over the NPIP family.

## 7. Bars and verdicts (integer comparisons; fixed now)

**S**, per UNITS arm X in {A, B}, at EVERY f in {0.1, 0.5, 0.9, 1.0} (f = 0 is the control, reported):
- **S1** M1 >= 9 (and >= the arm's own f = 0 value). **S2** M2 = 0 (S partners in NPIP). **S3** at f >= 0.5: relation precision >= 0.90 and recall >= 0.90.
  **S4** M1u (unfused copies in NPIP) >= BASE's at the same f.
- **UNITS-X WORKS on S** iff S1-S4 all hold; **PARTIAL** iff S2 and S4 hold and S1 or S3 fails; **FAILS** iff S2 or S4 fails.

**H**, per UNITS arm X in {A, B}, against BASE:
- **H1** c1(X) >= c1(BASE), and no copy placed in NPIP by BASE is lost (set inclusion; gained reported). **H2** partners dragged into NPIP <= BASE's.
  **H3** Compara chr16 bipartite F >= BASE - 0.01 and pairwise precision >= BASE - 0.01 (0.01 = half the house gain bar, as `PREREG_o1_cover_growth` C4).
  **H4** correct members lost <= floor(0.02 x the number of BASE correct members) on chr16 Compara, and 0 in the NPIP family.
- **UNITS-X PASSES H** iff H1-H4 hold. Liftoff recall and the r1017-style numbers are reported beside; a Liftoff decrease of more than one pair is named.

**Overall.** **PRODUCTIVE** iff some arm WORKS on S and PASSES H (then a junction detector is worth its own pre-registration; this
says nothing about one). **S-ONLY** iff an arm WORKS on S and no arm PASSES H (the execution is substrate dependent). **NOT PRODUCTIVE**
iff no arm WORKS on S and no arm PASSES H: the units direction is closed irrespective of detection. Anything else is PARTIAL
and the clauses that failed are named. **Attribution clause:** a WORKS / PASSES is attributed to the JUNCTION only if the
matched null N does not also clear the same clause; if N clears it, the result is reported as a node-shrinking effect
(r846 mechanism b), not as recovery. Nothing is pooled across species; the two substrates are separate rows.

## 8. Predictions (prior probabilities stated before any product) and falsifiers

P1. BASE (this build) M1 at f = 0/.1/.5/.9/1 is within 2 of the Python-F1v2 dev row 9/9/5/0/1 and M2 = 0 at f <= .9 (0.65).
P2. A and B keep M2 = 0 at every f: the partner's own 30-read transcript decides its unit's representative (0.75).
P3. A keeps M1 >= 9 at f <= 0.5 (0.80) and at f = 0.9, 1.0 (0.50): at f = 1 the copy unit is one truncated half (missing its first or last exon) and must still pass `--min-shared-exon-frac 0.60` against full copies and the containment escape.
P4. B equals BASE at f <= 0.5 wherever F1v2 already removed the fusion (bridge), and equals A at f >= 0.9 (0.70).
P5. S3 passes at f >= 0.5 for A (0.40): 10 pairs, co-duplicated SMG1-like partners (4) and SNX29 / CNOT3 / PDXDC / NSMCE paralogs must land in their own reference families with the partner unit's co-members unchanged.
P6. A is WORKS on S (0.40); PARTIAL (0.45); FAILS (0.15).
P7. H: c1(A) >= c1(BASE) (0.55); Compara F(A) >= BASE - 0.01 (0.60); H4 holds (0.35); A PASSES H (0.25). The null N clears H1 in 0.30: if random cuts also raise NPIP sensitivity, the gain is node shrinking.
P8. "Same ceiling as r1017": NPIP sensitivity does NOT fall as it did in r1017 (0.600 -> 0.550), because the units keep read-derived representatives and the `--min-cov-shorter 0.70` escape admits contained pieces (0.55); the price is partners/hubs, not sensitivity.
P9. Single-exon and tiny units exist in H and form hub-like nodes (degree >= 10) at a rate above BASE (0.7); none in S (all S halves have >= 3 exons, read from `fusions.json`).
Falsifiers: Z1 M2 > 0 at any f for A or B (the units leak partners even with the oracle); Z2 c1(A) < c1(BASE) at H (r1017's mechanism repeats); Z3 N clears the clauses A clears (the gain is not the junction); Z4 S3 fails at f = 0.5 with M1/M2 passing (families right, relations wrong: the relation definition, not the families, is the problem); Z5 the families stage fails or takes > 2x BASE's time on units.

## 9. Safety readouts (reported, no bars)

(i) units that are single-exon or tiny (exonic length < 600 bp = 2 x the aligner's 300 bp floor) and their degree in the pre-MCL graph (`mcl_families --dump-graph`, no other effect; its `clusters.tsv` must equal the driver run's); hub-like = degree >= 10, with BASE's degree distribution (maximum, 99th percentile) and the largest connected component / largest cluster of every arm beside (the "one component swallowing everything" check); (ii) units failing the gate that the fused locus passed (r1018): for each oracle locus clustered in BASE (same f), every unit locus that is unclustered in the arm, with the reason from a Python port of the deferred edge rule (validated equal to the dumped graph): no PAF record >= 300 bp at identity >= 0.70, exonic conjunct, shared-exon fraction < 0.60, or coverage (cov_longer < 0.30 and cov_shorter < 0.70); (iii) loci whose unit loci land in the SAME family (double count), the consolidation rule's effect on family sizes and NPIP precision; (iv) `/usr/bin/time -v` wall time and peak RSS of every families run, loci / edges / clusters per arm; (v) whether downstream stages run: `mcl_families`' own copy table (`fam.copies.tsv/.fa`, the contract `copy_assign --families` reads) for the S f = 1 BASE and B / A arms, `copy_assign --families` on them if cheap, and `tools/rustle_pipeline.sh catalog` on the S f = 1 BAM if cheap (the legacy catalog reads the BAM, not the families GTF: what changes for the copy of a unit is reported from the copy table: its `tid` is a unit id, its sequence the half's spliced exon sum).

## 10. Hostile self-review

1. **The oracle reads the answer.** Any S/H pass is an upper bound; it never validates a detector, and the S fusion junctions are the simulation's own.
2. **S is circular and small.** Error-free full-length reads, 30 reads per partner transcript, 10 fusions, 4 of them SMG1-like partners in one co-duplicated family; every half has >= 3 exons, so S cannot test single-exon units (H does); dev substrate, n = 10.
3. **Holder-based "copy in NPIP" counts a fragment as a copy.** A unit that is one half of a truncated transcript placed in NPIP passes M1 while the copy table would hold a truncated copy: §9 (v) reports it, the bars do not.
4. **Relation truth is this pipeline's own f = 0 partition** (family identity by co-membership, no external truth): a partner whose reference family is wrong is "correct" if the units reproduce it. The Compara / referee numbers of H are the external check; S has none.
5. **Regroup confound.** Arms A/B/N also re-derive every gene_id; A0 (regroup only) and BASE (F1v2, which regroups too) separate the regroup from the split, and N separates the junction from the shrinking. A0 keeps natural bridges, BASE removes them: A vs BASE differs in that too.
6. **Node-shrinking trap** (metric trap of 09-18): a rise in sensitivity can be one component swallowing everything; precision, partners dragged in, the largest cluster and N are reported beside every sensitivity.
7. **Bipartite precision is not tie-invariant** (r1045): scipy's assignment is the named policy and pairwise counts are reported beside; a gene-space scorer that intersects with the truth universe hides unlabelled predictions: the unjudgeable predicted pairs are reported per arm.
8. **Dev only.** chr16 is the human dev contig of every earlier test (22 / 24 / 26 of 26 NPIP copies are known on it); no held-out chromosome, sample or species is touched, so nothing here supports a default.
9. **The bar numbers are mine** (M1 >= 9, tolerance 0.01, 2% of correct members): each is stated before the product; the integer comparisons leave no room after it.
10. **What this cannot say:** whether a real detector finds the junctions; register 1166D found 0 of 348 / 12 / 52 / 56 real F1v2 bridges joining two multi-copy families, so a positive result is a statement about an execution without a known real-data gain case.

## 11. Order, machine rules, files

1. This file, then `lib/units.py` (+ tests on hand-made GTFs), `lib/score_s.py`, `lib/score_h.py`; gates G0-G3 on old products and binaries; then the runs in the order S assemble -> S families -> S score -> H -> §9.
2. Heavy (`mcl_families`, minimap2, `gw_family_catalog`) under `bash tools/rlock.sh heavy`, light (assemble < 2 GB and < 3 min, Python scoring) under `light`; foreground; never `pkill -f`; `TMPDIR` under `/mnt/linuxdisk`.
3. Scratch `/mnt/linuxdisk/tmp/rustle_figures_dev/container_units_mech/` (< 40 GB). Register rows are drafted with suffix J in the report (`figs/container_units_mech.md`), never appended.
4. Binaries (scratch `bin/SHA1SUMS`): copy_assign 1325e9d1, mcl_families 67b2d40f, as_table 0452b1b9, gw_family_catalog d9978567, family_score 27aa9445 (e163d955 build); frozen scorer binary `fj_bin_frozen/family_score` 7723029b. Driver `tools/rustle_pipeline.sh` 2c431091, `tools/rlock.sh` 30f424a9. Sim pool: `ann/fusions.json` 9ba2fdf5, `ann/arms.json` 0d0fd0b4, `ann/pool.tsv` 3e8d427d, `score.py` 8f512fb8.

## Amendments

**Amendment 1 (2026-09-30 10:58, after the S assemblies, oracle tables and S families runs, BEFORE any units arm of S was
scored or opened; only the ORIGINAL frozen-binary runs were scored, as gate G2).** Prereg sha1 before this amendment
`e5949d143e3955b05dc283694f4b96bc9782b022`. Tools (scratch `lib/`): `units.py` e5c45213, `test_units.py` 2c6842b6 (7/7 pass),
`score_s.py` 5c436014 (the version used for the first S scoring; a later change is a further amendment with its diff),
`oracle_s.py` c184a036, `g1.py` 2737aa19, `run_families.sh` e02c726d, `run_container.sh` dbbfc074, `run_s_assemble.sh`
9eb44591. Binaries and driver as §11 (verified).
- **Gates.** G0 pass. **G1 pass on all five S arms** (0 mismatches of 626 / 620 / 629 / 629 / 626 compared transcripts; same
  transcript sets in plain, regrouped and default GTFs). **G2 pass**: `score_s.py` on `fusion_container_sim/runs/f*/` gives M1 =
  9 / 1 / 0 / 0 / 4, M2 = 0 / 1 / 0 / 0 / 4, copies in NPIP 23 / 15 / 14 / 14 / 18, NPIP units 27 / 20 / 20 / 20 / 24, container
  M3 = 0 / 0.317 0.533 / 0.214 0.417 / 0.214 0.417 / 0.389 0.621, equal to the Outcome of `PREREG_fusion_container_sim` and to
  `score/score.json`. G5, G6 are run at scoring time (§ order below); G3, G4 belong to H.
- **What the S products are.** The default assembly (F1v2) finds 0 / 4 / 4 / 5 / 0 F1 bridge junctions and KEEPS (minority, share
  < 1/2) 0 / 4 / 1 / 1 / 0 at f = 0 / .1 / .5 / .9 / 1 (most fusion links at f = .5 tie at share 0.5 and abstain). Every one of the
  10 fused pairs has exactly ONE exact carrier of its fusion intron at f = .1 .. 1 (0 at f = 0, all by exact equality; no +-10 bp
  fallback was needed). A splits all 10 oracle transcripts (20 units) in every arm with fusion reads; B splits those still in
  `families.gtf` (6 / 9 / 9 / 10 of 10 at f = .1 / .5 / .9 / 1). Every arm's aligner command, read from its own
  `PREFIX.cache/paf/*/key.tsv`, is `minimap2 -x asm20 -c -X -N 50 -p 0.1 --secondary=yes -t 4` (2.30-r1287).
- **Identical inputs are not re-run** (md5 of the families input, `alias_of.tsv`): at f = 0 every arm has the same input (A0, B
  and BASE's `families.gtf` are byte-identical, A and N have an empty oracle), so one run stands for all; at f = 1 A0 == BASE
  (no bridge) and B == A (no bridge to remove). Runtime of each run (58-121 s of families, 2.0-2.6 GB) is in the report.
- **Prior caveat for P1.** The "dev row 9/9/5/0/1" of §0/P1 is the locus_fix_design extended simulation (803 transcripts at
  f = .5, Python F1v2), not this pool (636 transcripts): P1 is a weaker prior than §8 reads. No bar uses it.
- **S4 is read as:** unfused copies in NPIP of the arm >= unfused copies in NPIP of BASE at the same f (BASE's f = 0 value is 14 of 15).
- **Order still to run:** S scoring (REF = BASE f = 0), G5 on BASE / A of S, then H (its scorers, G3, G4 and their own amendment
  before any H arm is scored), then §9.

**Amendment 2 (2026-09-30 11:15; after the first S scoring, the human assemblies and all seven H families runs; BEFORE any H
product was scored or opened beyond its transcript/loci counts).** (The file's sha1 between Amendments 1 and 2 was not
recorded; Amendment 1 holds the sha1 of the text before it, which is the frozen design.)
- **A defect in the S scorer's relation definition, found after the first S scoring and fixed; both are reported.**
  `score_s.py` v1 (sha1 5c436014; first scoring, `results/S_scores.v1.json`) built the reference families F*_c / F*_p and tested
  the family match (§6.1 M3u) on SAME-STRAND member exons, as the prereg text says ("whose same-strand exons share >= 1
  base"). Families are strand-blind (an inverted paralog is a member), so for partners with an inverted paralog
  (LOC129527626 CNOT3-like, whose family at f = 0 holds the NPIPB1P copy on `-`) F* came out EMPTY and a correct relation was
  scored wrong. v2 (sha1 9708f3f9, `results/S_scores.json`) uses each member locus's own exons on any strand, and adds the STRICT
  variant reported beside the bar (the unit's cluster must be the PLURALITY cluster of the reference family's members in the arm,
  ties to the smaller id; the barred definition is unchanged and lenient: >= 1 shared member). M1, M2, M1u and the container M3
  do not use it and are identical in v1 and v2. The bar S3 is sensitive to this fix (v1: A 8/10 at f = .1 / .5 / .9, 9/10 at
  f = 1; v2: 10/10 at every f; strict 8 / 8 / 9 / 10 of 10); both are in the report. No bar, arm or threshold changed.
- **`units.py` piece naming made collision-free** (sha1 3f66eab0; v1 e5c45213): a non-keeper piece name `<g>.rg<k>` that is already
  the gene_id of ANOTHER input gene (possible only on an input F1v2 already regrouped, i.e. arm B) moves to the next free k. The
  Rust pass never meets this (it runs once on names without `.rg`). Found when the H arm B asserted; **every S units GTF and the H
  A0 / A GTFs are byte-identical under v2** (cmp), so nothing scored so far changed. H arm B and N were built with v2.
- **Gates.** G5 pass: replaying each of the 19 S arms' families run from its own cache with `--dump-graph` reproduces
  `clusters.tsv` and `loci.tsv` byte for byte (`lib/run_dump.sh`; the dump is `PREFIX.dump.graph.tsv`, one row per edge). G6 pass on
  every S and H units GTF (`lib/g6.py`). **H: G1 pass** (9,115 transcripts compared, 0 mismatches; 12 bridged gene_ids skipped);
  **G4 pass, byte-identical**: the new `off` GTF == `hdev/hsa16.BASE.gtf`, the default GTF and `families.gtf` == the Python F1v2's
  `hdev/hsa16.F1v2.gtf` / `.families.gtf`. G3 (H machinery on the dev BASE products) runs with `score_h.py`.
- **The H oracle** (`lib/oracle_h.py` e893a607, `H/oracle.tsv` sha1 0371dd99): 8,782 spliced chr16 transcripts; 582 overlap >= 2
  annotated same-strand genes (gene set: 1,791 records - 12 `readthrough` + 26 Dishuck copy records, 1,781 after replacing 24
  same-named records); 329 of those are UNCUTTABLE (their genes' exon-index intervals overlap: shared exons or nested genes) and
  are left whole; **253 oracle transcripts, 256 cuts** (250 give 2 units, 3 give 3), 14 of them F1v2 bridges (so not in arm B's
  input). A splits 253 -> 509 units; B splits 239 -> 480. Null draws for H use `stable_seed(S, tid)` with S = 20260930, 20260931, 20260932
  (N1-N3; the S null uses 20260930). The oracle cuts the copies in fused loci only where the transcript overlaps two genes
  cleanly; which of the 11 fused copies it reaches is reported per copy.
- **H runs done before scoring:** families of BASE, A0, A, B, N1-N3 (45-57 s each, 2.3-2.5 GB; aligner command identical, from each
  cache key). S arms aliased: f0.0 all arms -> BASE; f1.0 A0 -> BASE, B -> A.
- **Still to do before any H score:** `score_h.py` and its gate G3; its sha1 goes into Amendment 3.

**Amendment 3 (2026-09-30 11:40; BEFORE any H arm was scored).** `lib/score_h.py` sha1 a286d14d (imports read-only the frozen
`o1_cover_frozen/score.py` 25242385, `sim_score.py` 18781a28, `lo_score.py` ebfe6f57, `emu.py` fac9a560 and the frozen
`family_score` 7723029b). **G3 pass** on the dev BASE products (`hdev/hsa16.BASE.gtf`, `o1_cover/dev/hsa16/BASE.*`): Compara chr16
sens .4444 / prec 1.000 / F .6154, pairs 66 / 66 of 193, the binary-gate G2 true; NPIP 22 of 26 copies in NPIP, 8 of the 11 fused
copies; unchanged from `score.compara.json` / `npip_dev.json`. Same instrument, same session: r1017's stored cluster files
re-score to its register values exactly: `base_rerun` referee F .214, NPIP-Soto F .727 at sens .600, U2 F .610; `oracle_a` (a
declared re-derivation of r1017's oracle) NPIP-Soto F .687 at sens .550, U2 F .632, referee F .250 (r1017's own .235 was not
reproduced by the re-derivation; the register value is quoted beside it). The "NPIP-Soto" number of r1017 is
`family_score --soto soto.families.tsv --chrom chr16 --family NPIP`; "referee" is `fig7/current/human_chr16_ref.families.tsv` with
`fig7/current/gff/human_chr16.gff`, `--chrom chr16`. H is now scored exactly as §6.2 and §7 say, arms BASE, A0, A, B, N1-N3.

## Outcome

**Outcome (2026-09-30 12:30; S scored with `score_s.py` 9708f3f9, H with `score_h.py` a286d14d, bars by `lib/verdict.py`; Amendments 1-3
are the only deviations; no bar, arm, threshold or truth changed after a product was seen).** Scratch
`/mnt/linuxdisk/tmp/rustle_figures_dev/container_units_mech/` (`results/S_scores.json`, `H_scores.json`, `verdict.json`, `tables.md`); report
`scratchpad/figs/container_units_mech.md`.

**Verdict, literally as §7 states it: NOT PRODUCTIVE** (no arm WORKS on S: A and B are both PARTIAL; no arm PASSES H: H2 fails 4 > 3).
Read the effect sizes before the word: the units arms are far better than the default on both dev substrates, and the clauses
they miss are missed by 1-2 members for three identified reasons, none of them a cost of the execution.

| S (gorilla, 10 fusions) | f = .1 | .5 | .9 | 1.0 | bar |
|---|---|---|---|---|---|
| BASE (e163d955 defaults) M1 / M2 | 2 / 0 | 0 / 0 | 0 / 0 | 3 / 3 | |
| A0 (regroup only) M1 / M2 | 1 / 1 | 0 / 0 | 0 / 0 | 3 / 3 | |
| **A** (oracle units on plain) M1 / M2 / M1u | **7** / 0 / 14 | **7** / 0 / 14 | **8** / 0 / 15 | **9** / 0 / 15 | S1 M1 >= 9: **fails f <= .9**; S2 M2 = 0: holds; S4 M1u >= BASE: holds |
| **B** (oracle units on F1v2) M1 / M2 | 7 / 0 | 7 / 0 | 8 / 0 | 9 / 0 | same as A |
| N (random cuts) M1 / M2 / M1u | 1 / 1 / 14 | 0 / 0 / 14 | 1 / 1 / 15 | 5 / 4 / **11** | fails S2 and S4 |
| relation P / R, A (lenient, the bar) | 1.0 / 1.0 | 1.0 / 1.0 | 1.0 / 1.0 | 1.0 / 1.0 | S3 >= .90 at f >= .5: holds (B: P 1.0, R .9 .9 1.0) |
| relation P / R, A (strict plurality, reported) | .8 / .8 | .8 / .8 | .9 / .9 | 1.0 / 1.0 | tracks M1 |

S verdicts: A PARTIAL, B PARTIAL (S2, S4 hold; S1 fails). The v1 container's M3 is 0 / 0 for every units arm by construction (no unit holds a
partner base); M4 "yes" 0 and M5 10.6-11.9 kb everywhere (the container's M4/M5 are not changed by units). Gates G0-G2, G5, G6 pass (H: G1, G3, G4 in Amendments 2-3; the secondary instruments were validated too: GD, the edge-rule port admits exactly the
dumped graph in 5 arms of S and H, and GM, the emulated MCL reproduces the partition).

| H (human chr16, dev) | BASE | A0 | **A** | **B** | N1-3 |
|---|---|---|---|---|---|
| copies in NPIP /26 (fused /11) | 24 (10) | 24 (10) | **25 (11)** | 25 (11) | 24 (10) |
| NPIP members / precision (unit) / partners dragged in | 30 / .900 / 3 | 29 / .897 / 3 | 33 / .879 / **4** | 33 / .879 / 4 | 29 / .931 / 2 |
| by consolidated locus: members / precision | 30 / .900 | 29 / .897 | 30 / .933 / partners 2 | 31 / .936 | 29 / .931 |
| Compara chr16 sens / prec / F (pairs) | .500 / 1 / .667 (99) | .482 / 1 / .650 (87) | **.537 / 1 / .699 (109)** | same as A | .482 / 1 / .650 (87) |
| Liftoff copy-pair recall /36 | 17 | 15 | 16 | 17 | 15 |
| protein referee F | .236 | .231 | .254 | .254 | .231 |
| Soto-NPIP F (sens) / U2 F | .800 (.700) / .645 | .765 (.650) / .623 | .800 (.700) / .635 | same as A | .765 (.650) / .623 |
| correct Compara members lost vs BASE | - | 1 | **0 of 49** | 0 | 1 |

H bars: A and B: H1 holds (25 >= 24, no copy lost, NPIPB5 gained), **H2 fails by unit (4 > 3; by consolidated locus 2 <= 3 it would hold)**, H3
holds (F +.032, pairwise precision 1.0), H4 holds (0 lost of 49). **The three NPIP members A changes are single-exon unit
representatives** (NPIPB5's new holder `DN_chr16_22797602_2.U2`, 8 reads, 1,472 bp; two 478-511 bp partner units): r846's mechanism
(b), at the scale of one copy and two partners, not a component swallowing (largest connected component 120 -> 131, largest cluster 30 -> 33,
nodes of degree >= 10 218 -> 226 of ~1,100). The null N does not reproduce any gain (24 / 26, F .650).

**Why the bars are missed (all read from products).** (1) **S1, NPIPB12 (f <= .9): the fused locus is not split** - the copy's standalone
transcript overlaps a native neighbour transcript that overlaps the partner's transcripts, so the exon-overlap regroup keeps one locus (9 of 10 oracle transcripts separate at
f < 1, 10 of 10 at f = 1; H: 215 of 253 = 85%, 78 of 89 fused loci): the F1v2 "only link" condition returns as the limit of unit splitting.
(2) **S1, NPIPB13 (f <= .5): one knife-edge edge** - the pair (NPIPB13, the NPIPB3-like locus) is admitted at f = 0 (w .836) and rejected at
f = .5 by `low_shared_exon` (< 0.60) because minimap2's alignment moved with the target set; re-admitting that single edge returns NPIPB13 to
NPIP (emulated MCL on the validated graph); M1 is invariant under 30 random node relabellings (MCL ties are not the cause). (3) **H2: single-exon
units** (106 of 253 oracle transcripts give >= 1; 15 are two-exon transcripts = two single-exon units; the annotation's own overlaps, e.g.
LOC124907830 over NPIPB5's first exon, are cut as if they were fusions). NPIPA7 is outside NPIP at f = 0 already (M1 <= 9 by construction).

**Safety (reported, §9).** Units per arm: S 20 (A) / 12-20 (B), none single-exon; H 509 (A), 122 single-exon, 144 < 600 bp exonic, 216 either, 133 with
an edge, 98 in a family. Hub-likeness: 94 of 509 units sit in nodes of degree >= 10 (18%) against 20% of BASE's nodes. **Units failing where the fused
locus passed (r1018): S 0 at f <= .9 and 1 at f = 1 (the ATXN2L partner unit, no PAF record: correct); H 6 (A) / 5 (B) unit loci, all partner pieces,
none a copy; reasons by the validated edge-rule port (gate GD: admitted pairs == dumped graph in every arm checked): no_exonic, low_shared_exon,
or edges admitted but left unclustered; `coverage` fired for none.** Consolidation: 0 (locus, family) double counts in every arm; 8-9 of 10 S loci and 14 of
89 H loci inherit >= 2 families; by locus the NPIP precision is .933 (A). Cost: families wall 47-61 s for A / B (BASE 49-68 s), 2.0-2.6 GB, PAF 10.3 vs 10.9 MB:
no cost. Downstream: `copy_assign --families` runs on the units copy tables (S f = 1; 4.4-4.7 GB, 70-76 s; rows 2,692 BASE vs 2,168 A); the copy of a unit is its
half transcript (`<T>.U1` / `.U2`: NPIPB13's NPIP copy = 7 exons without the copy's first exon); the legacy catalog reads the BAM and is unchanged (38 s, 1.8 GB, 104 copies / 21 families).

**Comparison with r1017's node-replacement oracle.** Same small ceiling on the referee (+0.018 here, +0.021 there), a different NPIP: r1017's oracle fell
(Soto-NPIP .727 / .600 -> .687 / .550, reproduced by re-scoring its stored files); the units oracle does not (.800 / .700 unchanged, +1 copy). The
cause is NOT `--min-cov-shorter` (a control with `--min-cov-shorter 0`: BASE 24 / 26 and A 25 / 26, Compara F .667 / .699, unchanged) and NOT the
coverage gate (it fired for no unit): the execution differs (read-derived transcript units regrouped by exon overlap against gene spans with clipped
exons) and the baseline differs (BASE already has RG3 and F1v2, which recover NPIPB2 / NPIPB6).

**Scorecard.** P1 missed (BASE M1 9/2/0/0/3, not 9/9/5/0/1; weak prior), P2 hit, P3 missed (A 7/7/8/9), P4 half (B = A in M1 at every f), P5 hit
(lenient; strict misses at f <= .5), P6 hit (PARTIAL), P7 hit (A does not pass H), P8 hit for the outcome and **refuted for its stated cause**, P9 half
(single-exon and tiny units exist in H, none in S; the hub rate is NOT above BASE's). Falsifiers Z1, Z2, Z5 not fired; Z3 not fired on S, and on H the null clears H1 / H3 only by
equality (no gain); Z4 not fired.

**Post hoc, not barred.** Probe Q1 (drop the 122 single-exon units from A): NPIPB5's gain disappears (24 / 26), the 4 partners remain, Compara F .699:
a single-exon guard does not fix H2. Raw relation matching (strand-specific, v1) and the strict variant are in the report; the barred lenient
definition passes S3 while M1 fails, which is a weakness of the definition, not evidence for the arm.

**Hostile review of this Outcome.** (a) The verdict word rests on two integer bars (M1 >= 9 at every f; partners <= BASE); a bar of 8 would have
WORKED on S at f >= .9 only. (b) S is a designed simulation (n = 10, NPIPA7 outside NPIP at f = 0); H is dev chr16. No held-out substrate was touched, so
nothing here supports a default. (c) The H oracle is annotation overlap, not a read-proven fusion: 329 of 582 two-gene transcripts are uncuttable and
some cuts split one gene. (d) The relation bar is lenient; the strict variant tracks M1. (e) The S3 number depends on one scorer fix (Amendment 2). (f) F1v2
(B) adds nothing on top of the oracle units here (B = A on the NPIP unit-level and Compara numbers; they differ by one Liftoff pair, 17 vs 16, and one by-locus NPIP member).
