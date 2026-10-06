# Pre-registration: the UNIT split with the assembler's OWN locus rule, and real detectors in the loop (dev only)

**Written 2026-09-30 (KEY=cuv2), before any product of this test exists.** Dev only (gorilla fusion simulation S; human A119b
chr16 H). Python prototype in scratch `/mnt/linuxdisk/tmp/rustle_figures_dev/container_units_v2/`; nothing in `src/`, `tools/`,
`bench/`; nothing committed or pushed. Species are never pooled. No chimp (PTR) or orangutan (PPY) product is opened and no held-out
substrate is touched (A119b outside chr16, testis, KB3781, OR6737 outside its dev contig: none). This stage is DESIGN on dev.
The orchestrating session's task is the mandate; the two sibling pre-registrations `PREREG_container_units_mechanism_2026-09-30.md`
(the oracle execution test, its Outcome) and `PREREG_container_units_definition_2026-09-30.md` (the detector definitions, its Outcome)
are the parents of this one.

## 0. What was seen before this file (and what was not)

**Seen.** The two parents in full, including Outcomes and Amendments; the mechanism report (`figs/container_units_mech.md`) and the
definition report (`figs/container_defs.md`); the scripts of the mechanism test (`lib/units.py`, `score_s.py`, `score_h.py`,
`oracle_s.py`, `oracle_h.py`, `verdict.py`, the run scripts) and the frozen definition instruments (`cd_defs.py`, `cd_scan.py`,
`cd_lib.py`, `run_defs.py`, `products.py`); the Rust source of `collapse_loci_groups` (`family_detect.rs`), its caller in
`copy_assign.rs` (the GTF emitter: `gene_id` = the raw base tid of the component representative), `regroup_gtf_lines` (RG3), all of
`bridge_regroup.rs` up to `rg3_pieces`, and `gtf_loci` in `mcl_families.rs`. From the products of the PARENT tests (not of this
test): the per-arm tables of the mechanism Outcome (S M1 / M2 per arm, H per arm), the stored F1 bridge tables of S
(`S/f*/BASE/f*.bridges.tsv`: 0 / 4 / 4 / 5 / 0 F1 bridge junctions, 0 / 4 / 1 / 1 / 0 kept by F1v2 at f = 0 / .1 / .5 / .9 / 1), the
call counts of the frozen definitions on chr16 (`calls.pkl`: R 179 cuts on 170 transcripts, R2 38 on 38, WA 4,974 on 929), and
ONE exploratory look at the NPIPB12 locus of the plain S f = 0.5 GTF (transcripts 549-554 of the plain GTF: the fusion transcript
550, the partner transcripts 549 / 551 / 552 / 553 and the copy's standalone transcript 554): 552 / 553 share the junction
(104292398-104292786) with the partner transcripts 549 / 551 and their last exon overlaps the first exon of 554, which shares
junctions only with the copy side of 550; hence my prediction P-A2 below is informed by that look. **Not seen:** any product of a
native-regroup arm (no unit arm of this test has been built), any detector-in-the-loop product, any score of this test, any W call
on the S substrate (no PAF of a plain-GTF S families run exists yet).

## 1. Question and scope

The mechanism test (oracle junctions, exon-overlap regroup) lifted the gorilla simulation from fused copies in NPIP 2 / 0 / 0 / 3
to 7 / 7 / 8 / 9 of 10 (f = .1 / .5 / .9 / 1) and human chr16 from 24 to 25 of 26 NPIP copies, at no cost, and missed its bars
for three reasons: (1) the exon-overlap regroup keeps a locus fused when a THIRD transcript overlaps both units (NPIPB12); (2) one
knife-edge gate edge (NPIPB13, unrelated to fusion; reported, not chased); (3) single-exon units. The definition study found no
detector of unit boundaries that is both precise and reaches the hard strata (R2 precision .24 / .87 / .77 at recall .07 / .04 / .08,
R reaches 35-48% of the dominant links at precision .07 / .58 / .34, W at chance).

**Part A (execution, oracle junctions).** Does re-deriving gene_ids with the assembler's OWN locus rule (junction-sharing
union-find; single-exon units attached by exon overlap) separate what the exon-overlap regroup left fused and lift M1 without
losing a member? **Part B (detector in the loop).** With each real detector's junction list in place of the oracle, what is kept
of the oracle's gain, what is paid in false splits, and which detector, if any, is worth a Rust port? **Part C.** The frozen rule
in words and an implementation note (no Rust is written). Nothing here can justify a default: every oracle reads truth and both
substrates are dev. Register rows to respect (not re-proposed): r846 (node cut doubles and makes hubs), r845, r967, r1017 / r1018,
r1053 / r1054 (node exon-model changes evict or make hubs), r1184-r1188 (a cover for F1v2 bridges).

## 2. The execution (binding)

### 2.1 Units and the families input
Exactly the mechanism test's §2.1-2.2: a transcript T with oracle / detector junctions J(T) (introns as 1-based closed (s, e)) is
REPLACED in the FAMILIES INPUT ONLY by m + 1 units `<T>.U1 .. <T>.U<m+1>` (transcription order; same contig, strand, `reads` of T;
attributes `fusion_of`, `fusion_unit`, `fusion_junction`), in T's line position; the assembled GTF itself is untouched; the families
stage is `mcl_families --from-gtf <input> --fasta <genome> --threads 4 --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units --out
PREFIX.fam` (the driver's families command; e163d955 binaries; defaults incl. `--min-cov-shorter 0.70`; the aligner inside
`mcl_families`, identical in every arm; `--dump-graph` is appended, which changes no product: gate G5'). The input is the plain GTF
(`--bridge-regroup off` product) whatever the detector.

### 2.2 E_nat: native regroup (arm family A1; the registered execution)
After the split, every gene_id of the families input is re-derived (the input gene_ids are discarded):
1. **Union-find on junctions.** Over ALL families-input transcripts in line order, keyed `(contig, donor, acceptor)` with donor = the
   exon's last base and acceptor = the next exon's first base (GTF 1-based closed), strand-blind: transcripts sharing ONE such key are
   in one component. This is an exact port of `family_detect::collapse_loci_groups` (union on first owner of each key, components by
   root, members in ascending line index).
2. **Representative** of a component = max over members of (reads, span, -line index), `span` = last exon end - first exon start + 1
   (the assembler's own tie-break: most reads, then longest span, then earliest index).
3. **Attachment (new; the assembler leaves a junction-less transcript as its own locus).** A single-exon UNIT (a unit with exactly one
   exon; original single-exon transcripts are NOT touched) joins the component of the transcript with >= 1 intron, same contig and
   same strand, with which it shares the most exonic bases (ties: the lower families-input line index); with no such overlap it is
   its own component. Attachment is computed after step 1 (a single-exon unit never joins another junction-less transcript and
   never links two components).
4. **Names.** The gene_id of the representative's INPUT gene_id is the component's name source o(G) (a unit's input gene_id is its
   parent's). Among the components with the same o(G) = g, the one whose representative is best by (reads, span, -index) keeps g and
   the others become `<g>.nat<k>`, k = 2.. in the order of their representatives' line index, skipping any name already taken.
   Nothing else about a transcript changes.
**No constant.** Every rule is a port of the assembler's own or a tie-break by line index.

### 2.3 Variants and comparators
- **A1s** (registered variant, H and S): step 3 attaches EVERY single-exon transcript of the input, original ones included.
- **A (RG3 execution, the comparator)**: the mechanism test's arm A (exon-overlap regroup, `units.py` 3f66eab0): products are REUSED
  (same binaries and inputs), re-scored with this test's scorers.
- **NAT0**: E_nat with an empty cut list (the native regroup alone): the control that separates the regroup from the split.
- **N1 (null under E_nat)**: the oracle transcripts cut at m uniformly random introns each (seeds as the mechanism test: S
  `stable_seed(20260930, tid)`, one draw; H 20260930 / 31 / 32).
- **A2 (conditional, §2.4).**

### 2.4 A2: only if E_nat leaves a third-transcript link
A residual link = a split transcript T (S oracle, H oracle) whose consecutive units U_i, U_{i+1} are still in ONE component after
§2.2. The residual links are listed per f (S) and for H with the linking transcripts, found as: the articulation transcripts X (X not
a unit of T) such that deleting X from the junction-sharing graph of the component and re-running step 1 puts U_i and U_{i+1} in
different components. If there is >= 1 residual link, **A2** is one more arm of the oracle execution: each articulation transcript
X keeps its junction edges only toward the side (component of the graph without X) with which it shares the most distinct junctions
(ties: the side containing the lower-indexed unit) and its edges toward the other side are dropped; the union-find is re-run once
with that edge set; representative, attachment and names as §2.2. A2 is constant-free (graph structure plus line order), is applied
only around split transcripts, and is reported with the number of links it separates and what the families do; it is run iff
there is >= 1 residual link after E_nat on S or H, and never tuned.

### 2.5 E*: the execution Part B uses (selection rule, fixed now)
Part A scores E_nat (A1), A1s, A2 (if any) and the comparator A with the SAME scorers. **E\*** = the execution with (1) the fewest
members lost (S: §6.1 L_S summed over f = 0 .. 1 incl. copies and partners; H: correct Compara members lost + NPIP copies lost),
then (2) the larger sum over f in {.1, .5, .9, 1} of S M1 plus H copies in NPIP c1, then (3) the larger H Compara TP pairs, then
(4) the order A1, A1s, A2, A. Part B is run only after E* is known; it is recorded in an amendment before any Part B arm is built.

## 3. Detectors (binding; junction lists on the PLAIN transcripts)
A cut = (T, j): transcript id of the plain GTF and intron index j in genomic order, converted to (s, e) by the transcript's exons
(asserted to be an intron of T). Per substrate:
- **ORACLE**: S = the exact carriers of each fused pair's fusion intron (`oracle_s.py`; 10 per f > 0, none at f = 0); H = the 253
  annotation-overlap oracle transcripts, 256 cuts (`oracle_h.py`, `H/oracle.tsv`). Reused from the mechanism test.
- **R2** = the rows of the assembler's `PREFIX.bridges.tsv` (default `--bridge-regroup f1v2` assemble) with keep = True: every
  transcript of T_J cut at J. **R** = every F1 bridge junction (`bridge` = True of `bridge_junctions.tsv` = all rows of `bridges.tsv`
  whatever the share decision): F1 without the share rule. S: from the assemble products of the e163d955 binary run with
  `--bridge-regroup f1` on each arm's BAM (gate G-R); H: the frozen `container_defs` R / R2 calls of chr16 (`calls.pkl`), gated
  against the H `bridges.tsv`.
- **W** = the frozen WA calls (`cd_defs.w_calls`, variant WA: witnesses = every non-overlapping locus with an exon-exon column,
  graph = the shipped admitted edges) on the plain GTF and its families product (PAF with CIGARs, admitted graph of
  `--dump-graph`). H: the stored chr16 calls (929 transcripts, 4,974 cuts). S: computed here with the frozen instruments unchanged
  on new plain-GTF families runs (f = 0 uses the BASE f = 0 product, whose input equals the plain GTF byte for byte).
- **RuW** = R union W (as `(T, j)` sets). **R|W** = R's cuts on transcripts with >= 1 W call (the definition study's post hoc DEV
  lead, **labelled post hoc**).
- **NULL_W** = W's calls with each call replaced by one uniformly random intron of the same transcript (seeded
  `stable_seed(20260930, "W", tid)`; duplicates collapse). **NULL_X** for every detector X whose S or H verdict is not NOT, built
  the same way (conditional stage B2).
A detector whose list is empty for an arm gives that arm the NAT0 input (aliased, not re-run).

## 4. Arms and substrates
- **S** (gorilla fusion simulation, `fusion_container_sim`, 10 fusions; nested subset BAMs f = 0, .1, .5, .9, 1; sha1s as the
  mechanism §11): BASE (default assemble, `families.gtf`), A0, A (mechanism products, reused and re-scored), **NAT0, A1, A1s, N1,
  (A2)**, and in Part B per detector **R2, R, W, RuW, R|W, NULL_W** (+ conditional nulls), f = 0 .. 1 (f = 0 is the false-split
  control: no fusion reads, every cut is false).
- **H** (human A119b chr16, dev; `--region chr16:0-96330374`; the plain GTF `H/PLAIN/hsa16.gtf`, byte-identical to the dev
  `hsa16.BASE.gtf`): BASE, A0, A, N1-N3 (mechanism products, reused), **NAT0, A1, A1s, N1-N3 under E_nat, (A2)**, and per detector
  the same Part B arms.
Identical families inputs (md5) are run once (alias table). Every arm is its own all-vs-all (no PAF is shared across arms).

## 5. Gates (a failed gate stops that substrate; nothing is scored through it)
- **G0** binaries `bin/SHA1SUMS` (copy_assign 1325e9d1, mcl_families 67b2d40f, as_table 0452b1b9, gw_family_catalog d9978567,
  family_score 27aa9445); driver 2c431091; `rlock.sh` 30f424a9; `sha1sum -c` of the frozen definition instruments
  (`container_defs_frozen/SHA1SUMS`) and `test_cd.py` passes; help states `--bridge-regroup` default f1v2, `--min-cov-shorter` 0.70.
- **G-A1a (the port is the assembler's rule).** On UNPOLISHED plain GTFs (`copy_assign --assemble-only --assembly-junctions strict
  --bridge-regroup off --gtf-tpm`, `--assembly-polish none`; S f = 0 .. 1 and H) the §2.2 step 1-2 port gives, for EVERY
  transcript, `gene_id == base_tid(transcript_id of its component representative)` (base tid = the id without a trailing `.<digits>`).
  Zero mismatches required; the number of gene_ids holding >= 2 components (raw-tid collisions) is reported.
- **G-A1b (on the shipped plain GTFs).** No component of the port spans two gene_ids (0 merges); the gene_ids the port splits into
  >= 2 components are counted and compared with RG3's exon-overlap pieces of the same GTF (both split / port only / RG3 only), each
  class with its cause (ghost: the polish dropped a bridge; collision; strand-crossing junction).
- **G-U** every units GTF parses, ids are unique, exons strictly ordered, the exon intervals of the units of T equal T's exons,
  every unsplit transcript is present with equal exons / reads / strand / contig (only its gene_id may differ; the mechanism's
  `g6.py`), and the builder is deterministic (two builds byte-identical).
- **G-S / G-H (scorer fidelity).** The scorers of this test reproduce, on the stored mechanism products, every S field of
  `S_scores.json` for BASE and A at the five f (M1, M1u, M2, copies in NPIP, relation counts, clusters) and every H field of
  `H_scores.json` for BASE, A and N1 (c1, fused, NPIP members / precision / partners by unit, Compara core and pairs, Liftoff, referee,
  Soto-NPIP, U2, lost vs BASE). The NEW readings (S lost members, H partners by consolidated locus) are added to, never substituted
  for, those fields.
- **G-W** the glue that computes W on a new product (graph reader, `w_calls`) reproduces the stored chr16 `calls['WA']` exactly from
  the stored chr16 scan and graph; `cd_scan` runs on each new S PAF with `unknown_locus = 0`.
- **G-R** S: the `bridge_junctions.tsv` of `--bridge-regroup f1` equals that of `f1v2` (same rows) and the f1v2 `bridges.tsv` rows are
  its bridge = True rows; H: the frozen R / R2 calls equal the rows of the H `bridges.tsv` (all / keep).
- **G5'** for the S f = .5 and H A1 arms (and the first detector arm of each substrate): the driver's own `families` command
  (no `--dump-graph`) replayed from the arm's cache gives `clusters.tsv` and `loci.tsv` byte-identical to the `--dump-graph` run.
- **G-D** the edge-rule port (`gate_diag.py`) admits exactly the dumped graph in the arms it is used on.

## 6. Measurements
### 6.1 S (per arm, per f; scorer = the mechanism's `score_s.py` 9708f3f9 + additions)
M1 (fused copies in NPIP /10), M1u (unfused /15), copies in NPIP /25, M2 (partners /20 whose holder is in NPIP), relation P / R
(lenient and strict; one record per CUT transcript, correct iff the transcript is an oracle carrier and the copy unit and partner unit
match F*_c / F*_p of BASE at f = 0; precision = correct / emitted over ALL cut transcripts of the detector, recall = pairs with >= 1
correct record / 10), clusters, largest cluster. New:
- **Separation.** T is separated iff every pair of consecutive units lies in distinct new gene_ids; counts of the 10 oracle
  transcripts (by f), NPIPB12 by name.
- **L_S, S lost members** (binding definition). Gene universe = the gorilla RefSeq genes of `ggo_npip_sim/ann/ggo3.gff` on the three
  sim contigs. A locus labels the genes sharing >= 1 same-strand exon base with its transcripts (each transcript's own strand); a
  family holds the genes its member loci label (a folded locus counts through its unit); two genes are CO-CLUSTERED in an arm iff a
  family holds both. P_REF = co-clustered pairs of BASE at f = 0. For BASE(f) the MATES of g are the genes h with (g, h) in P_REF
  and co-clustered in BASE(f); g is a correct member of BASE(f) iff it has a mate, and is LOST in X(f) iff no mate is co-clustered
  with g in X(f). Reported: lost genes, the copies and partners among them, the distinct P_REF families containing a lost gene.
- **False-split census** (detector arms; also f = 0): cuts, cut transcripts, cuts in the ORACLE (exact), cut transcripts lying inside
  exactly one annotated gene (FP_N) and the distinct genes they sit in; the family outcome of every cut transcript: **SAME** (its
  units sit in one family, >= 1 clustered, none unclustered: consolidation makes the split invisible), **DIFF** (>= 2 different
  families: a new relation / cover record), **ONE_UNCL** (>= 1 clustered and >= 1 unclustered unit), **ALL_UNCL**.

### 6.2 H (per arm; scorer = the mechanism's `score_h.py` a286d14d + additions)
NPIP (26 Dishuck copies): c1 (copies in the NPIP family), fused /11, members, precision, partners by unit, **partners by consolidated
locus**, per-copy fate of the 11 fused copies, the four readthrough copies (NPIPA1, NPIPA6, NPIPA9, NPIPB14P) and PKD1P6-NPIPP1
explicitly; Compara Primates chr16 bipartite sens / prec / F and pairwise TP / predicted / truth; Liftoff copy-pair recall /36;
protein referee F; Soto-NPIP F / sens; U2 F; clusters, largest; correct members lost vs BASE (Compara; the mechanism's function) and
NPIP copies lost (set inclusion). New:
- **Consolidated locus** (binding): a new gene_id holding >= 1 unit is counted as its INPUT gene_id (the pre-split locus); every other
  new gene_id as itself. "Partners by locus" = NPIP-family consolidated loci none of whose member loci overlaps a truth copy. (The
  mechanism test counted by input gene for every new gene_id, merging the pieces of loci that held no unit; the new rule keeps them.)
- **Separation**: oracle transcripts separated (of 253; the mechanism test's A: 215) and fused loci separated (input gene_ids holding
  >= 1 oracle transcript all of whose oracle transcripts are separated; the mechanism's A: 78 of 89), and which loci stay fused.
- **False-split census** of each detector arm: cuts, cut transcripts, cuts in the oracle (exact) and in a truth WINDOW of the
  definition study's truth table, FP_N cut transcripts and the genes they sit in; the SAME / DIFF / ONE_UNCL / ALL_UNCL outcome of every
  cut transcript; the DIFF records listed with the genes of each unit and their Compara families (judged: annotation-supported =
  the units overlap disjoint non-empty annotated gene sets; single-gene = all units overlap one gene); Compara families whose
  correct-member set shrank (from the lost list).
- **Safety readouts** (mechanism §9): units single-exon / < 600 bp, degrees, hubs (degree >= 10), largest connected component /
  cluster, units failing where the fused locus passed (`safety1.py` with the dumped graph), cost (wall, RSS).

### 6.3 A3: the container OUTPUT spec v2 (written, and emitted for the A1 oracle arms)
The spec's column layout is in the report (§A3); a prototype emitter (`relations.py`) writes `relations.tsv`, `members_by_locus.tsv`
and the v1 residual for the A1 oracle arms of S (f = .5, 1) and H, and the invariants are asserted (every relation row's unit
families equal `clusters.tsv`; one member per locus per family; every unit in exactly one row; the copy table rows of unit loci
unchanged in format).

## 7. Bars and verdicts (integer comparisons; fixed now)
### 7.1 Part A (execution)
The mechanism bars are RE-RUN for A1 (and A1s, A2) as information: S1 (M1 >= 9 and >= the arm's f = 0 value at f = .1 .. 1), S2 (M2 =
0), S3 (relation P and R >= .90 at f >= .5), S4 (M1u >= BASE's); H1 (c1 >= BASE, no BASE copy lost), H2 (partners <= BASE's, by unit
AND by consolidated locus: both stated), H3 (Compara F and pairwise precision >= BASE - .01), H4 (correct members lost <= floor(.02 x
49) = 0, none in NPIP). The absolute S1 cannot hold at f <= .5 while NPIPB13's edge is knife-edge (reported, not chased; NPIPA7 is outside
NPIP at f = 0 in BASE, ceiling M1 <= 9). **The Part A questions (no verdict word; they feed E\*):**
- **QA1** separation: A1 separates >= as many S oracle transcripts as A at every f (A: 9 / 9 / 9 / 10 at f = .1 / .5 / .9 / 1 of 10) and
  >= as many H fused loci as A (78 of 89); the loci still fused are listed.
- **QA2** no-regret: M1(A1, f) >= M1(A, f) at every f, M2(A1) = 0 at every f, no member lost (S and H).
- **QA3** H: c1 >= 25, partners by locus <= BASE's, Compara F >= A's - .01.
- The attribution clause of the mechanism test is kept: N1 (random cuts under E_nat) must not clear the same clause.

### 7.2 Part B (detectors, bars RELATIVE to the oracle arm of E*)
gain_O(f) = M1(oracle, E*, f) - M1(BASE, f) on S; gain_X likewise. On H: gain in copies c1 and in Compara TP pairs vs BASE.
- **S verdict.** NOT iff at any f in {0, .1, .5, .9, 1} a member is lost (L_S > 0 vs BASE(f)) or M2(X) > M2(BASE) or M1u(X) < M1u(BASE).
  Else WORKS iff at every f in {.1, .5, .9, 1} with gain_O(f) > 0: 4 x gain_X(f) >= 3 x gain_O(f). Else PARTIAL iff gain_X(f) >= 1 at some f;
  else NOT (no gain).
- **H verdict.** NOT iff a BASE copy is lost from NPIP, or a correct Compara member is lost (> 0), or partners by consolidated locus
  exceed BASE's, or Compara F or pairwise precision < BASE - .01. Else WORKS iff both 4 x gain_X >= 3 x gain_O hold (copies c1 and Compara TP
  pairs, each where gain_O > 0); PARTIAL iff gain_X >= 1 in c1 or pairs; else NOT.
- **Overall.** WORKS iff S and H are WORKS; NOT iff either is NOT; else PARTIAL.
- **Attribution.** A WORKS / PARTIAL is attributed to the JUNCTIONS only if the matched null NULL_X (run for every detector that is not
  NOT) does not clear the same clause; otherwise it is reported as a node-shrinking effect (r846 mechanism b).
### 7.3 Decision rule among detectors (lexicographic, no fitted constant)
(1) zero correct members lost on S (all f) and on H; (2) then the larger S sum over f in {.1, .5, .9, 1} of copies in NPIP /25, then
the larger H c1, then the larger H Compara TP pairs; (3) then fewer cuts (H cuts first, then S cuts summed over f). A detector with a
lost member is ranked only among those with lost members (by the number lost, fewer first). The detectors ranked are R2, R, W, RuW and
R|W; ORACLE is the reference, NULL_X the control.

## 8. Predictions (this author's probabilities, before any product) and falsifiers
**Part A.**
- P-A1 gates G-A1a (exact names on every unpolished GTF) 0.93; G-A1b (0 merges) 0.97; G-W, G-R, G-S, G-H 0.85 jointly.
- P-A2 A1 separates NPIPB12 at f <= .9 (0.80) and all 10 oracle transcripts at every f > 0 (0.70).
- P-A3 M1(A1, f) >= M1(A, f) at every f (0.80); M1(A1) = 8 / 8 / 9 / 9 at f = .1 / .5 / .9 / 1 (0.35), >= 8 / 8 / 8 / 9 (0.60).
- P-A4 M2(A1) = 0 at every f (0.80); M1u(A1) >= BASE at every f (0.90).
- P-A5 H: A1 separates >= 85 of the 89 fused loci (0.55) and >= 240 of the 253 oracle transcripts (0.60).
- P-A6 H: c1(A1) >= 25 (0.55), no BASE copy lost (0.75), partners by consolidated locus <= 3 (0.50), Compara F >= .699 (0.45), 0 lost
  of 49 (0.55); A1 passes H1-H4 by locus (0.30).
- P-A7 attachment matters: A1 and A1s differ in >= 1 NPIP placement on H (0.40); identical S products (0.85).
- P-A8 A2 is needed (>= 1 residual link after E_nat on S or H) (0.55); E* = A1 (0.55), A1s (0.15), A2 (0.10), A (0.20).
- P-A9 N1 clears no Part A clause that A1 clears (0.85).
**Part B.**
- P-B1 S: R2 and R are NOT WORKS (F1 needs a standalone side, R finds <= 5 of 10 fusions, none at f = 1): PARTIAL at best (0.90); R recovers
  >= 75% of the oracle gain at f = .9 (0.15), R2 at f = .1 (0.40).
- P-B2 S: W cuts >= 50 transcripts at f = 0 (0.60) and loses >= 1 member (0.75): W is NOT on S (0.75).
- P-B3 H: W and RuW are NOT (a lost Compara member or a lost NPIP copy) (0.90); W cuts >= 10x the oracle's transcripts' units (certain: 4,974 cuts).
- P-B4 H: R is NOT (>= 1 lost member) (0.55); R2 is not NOT (0.65) and keeps < 75% of the oracle's gain (0.85); R|W is not NOT (0.35).
- P-B5 No detector WORKS on both substrates (0.85); the detector ranked first by §7.3 is R2 (0.40), R|W (0.15), R (0.15), W / RuW (0.05), none
  without a lost member (0.25).
- P-B6 NULL_W clears no bar (0.85).
**Falsifiers.** Z1: M2(A1) > 0 at any f (the native rule leaks partners even with the oracle). Z2: A1 loses a member that A keeps. Z3: N1 clears
the clauses A1 clears (the gain is not the junction). Z4: a detector WORKS on both substrates (the definition study's pessimism is
overturned in the loop). Z5: RuW loses fewer members than R on H (W is not a flood when its cuts go through the families stage).
Z6: A1 needs > 2x BASE's families time.

## 9. Hostile self-review
1. **Oracles read the answer.** ORACLE arms are upper bounds; S is circular (error-free reads, 10 fusions, 4 SMG1-like partners in one co-duplicated
   family); the H oracle is annotation overlap, cuts some single genes and cannot cut 329 of 582 two-gene transcripts.
2. **Part A was designed after the mechanism Outcome and after one look at NPIPB12's transcripts.** The native rule is the assembler's
   own and was not tuned, but the choice to attach only single-exon UNITS (A1) versus every single-exon transcript (A1s) was made to keep
   non-oracle loci identical to the assembler's; both are run and E* is selected by a rule fixed here.
3. **The by-locus reading of partners** (new consolidated-locus rule) was chosen after seeing that the mechanism test failed H2 by unit (4 > 3)
   and passed by locus (2 <= 3); both are stated and the barred reading is the by-locus one because the A3 output (members by locus) is what a
   consumer reads. The alternative reading would have failed the mechanism's arms.
4. **S lost members are measured against this pipeline's own BASE at f = 0** (co-membership, no external truth); H uses Compara. The
   annotation cannot grade duplication-module boundaries: W's "false splits" may be real modules (definition Outcome); in this test they are
   judged by what they do to families, not by the annotation alone.
5. **Dev only.** chr16 is the human dev contig of every test; R and R2 carry no fresh evidence (spent on the held-out substrates of the bridge tests); nothing
   here supports a default. The detectors were NOT tuned on S or H; W\* was frozen by the definition study.
6. **Node-shrinking trap** (r846): every gain is read against NAT0 (regroup alone) and against the nulls; partners, the largest cluster and lost
   members beside every sensitivity.
7. **Bars are mine** (75% retention, 0 lost members from the task, floor(.02 x 49) = 0): stated before the product; the integer comparisons leave
   no room afterwards; a single lost member turns a detector into NOT, which is harsh on purpose (the task's rule (1)).
8. **Two-stage design (E\*).** Part B's execution is chosen from Part A by a fixed rule; the detectors are therefore compared under the execution
   that best served the ORACLE, which favours no detector in particular.
9. **The S W calls are computed on the plain-GTF families product of THIS build** (definition study: the default GTF of 09-25 with older binaries):
   G-W covers the glue, not the alignment; the PAF of a plain-GTF S product differs from the BASE's (bridges are loci there).
10. **What this cannot say:** that a detector exists for the dominant or standalone-free fusions (register 1166D: 0 of 348 / 12 / 52 / 56 real F1v2
    bridges joined two multi-copy families); anything about held-out substrates; a default.

## 10. Order, machine rules, files
1. This file; then the tools (`lib/native.py`, `units2.py`, tests), the gate runs (G-A1a needs unpolished assemblies; G-A1b, G-U, G-S, G-H, G-W, G-R),
   Amendment 1 (tool sha1s and gate results, before any unit arm of this test is built or scored). 2. Part A runs and scores (S, H), E*; Amendment 2
   (E*, before any Part B arm is built). 3. Part B (detector lists, arms, scores, census); conditional nulls. 4. A3 tables, Part C, Outcome.
2. Heavy (`mcl_families`, minimap2) under `bash tools/rlock.sh heavy` with `RLOCK_WAIT` / `RLOCK_TIMEOUT` raised for long waits and runs (the
   lock is shared with another session: wait, never bypass); light (assemblies < 2 GB and < 3 min, Python scoring) under `light`; foreground;
   never `pkill -f`; `TMPDIR` under `/mnt/linuxdisk`; one heavy run at a time. Scratch `/mnt/linuxdisk/tmp/rustle_figures_dev/container_units_v2/`
   (< 40 GB). No PTR / PPY, no held-out substrate.
3. Register rows are drafted with suffix K (first number 1194K) in the report `scratchpad/figs/container_units_v2.md`, never appended.

## Amendments

**Amendment 1 (2026-09-30 12:13; BEFORE any native-regroup unit arm was built or scored; only gate artefacts exist).** Prereg sha1 before this
amendment `6bc1013be8bd7bd695206c68ebba888a249cbf05` (written 12:05:46; a byte copy is `notes/PREREG.frozen_v0.md` in scratch); this is the frozen
design. Tools (scratch `lib/`): `units.py` 3f66eab0 (the mechanism's, unchanged), **`units2.py` 682fadeb** (E_nat / A1 / A1s / NAT0 / N1 / A2 and the
rg3 comparator mode), `gate_native.py` 2bfb6768, `test_units2.py` 0403bbf2 (9 tests pass), `run_s_gate_asm.sh` bad5ebd4, `run_s_f1.sh` 076dc6a9,
`run_h_gate_asm.sh` 8dcac76f. Binaries `bin/SHA1SUMS` verified (copy_assign 1325e9d1, mcl_families 67b2d40f, as_table 0452b1b9, gw_family_catalog
d9978567, family_score 27aa9445); driver 2c431091; `rlock.sh` 30f424a9; frozen definition instruments `sha1sum -c` OK (`SHA1SUMS` 7cf79b24) and
`test_cd.py` 24 / 24 pass.
- **G0 pass. G-A1a pass on all six unpolished plain GTFs** (S f = 0 / .1 / .5 / .9 / 1: 648 / 658 / 658 / 658 / 648 transcripts; H: 28,270): 0 name
  mismatches, 0 components spanning two gene_ids; gene_ids holding >= 2 components (raw-tid collisions) 2 per S arm and 56 in H. The S / H unpolished
  assemblies are `S/*/NOPOLISH`, `H/NOPOLISH` (9-15 s each, < 0.85 GB).
- **G-A1b pass on the shipped plain GTFs** (S 626 / 636 / 636 / 636 / 626 transcripts; H 9,473): 0 components span two gene_ids, 0 strand-crossing
  components; the port splits 2 gene_ids per S arm (both raw-tid collisions; RG3 splits none) and 84 gene_ids of H (RG3: 19): 19 split by both (18
  ghost candidates, 1 collision), **65 by the port only (55 collisions, 10 ghost candidates)**, 0 by RG3 only. The native regroup is therefore FINER than
  RG3's exon-overlap pieces on non-fusion loci too (two transcripts with no shared junction are two loci), which is why NAT0 (native regroup without
  cuts) is a separate arm and every gain is read against it.
- **Builder fidelity (G-U-A) pass:** `units2.py --regroup rg3` reproduces the mechanism's arm A (oracle) and A0 (regroup only) and N (random cuts) GTFs
  BYTE FOR BYTE on S f = .1 / .5 and on H (A, N1, N2, N3), so the new builder / writer is the mechanism's.
- **G-R pass:** S: the `bridge_junctions.tsv` of `--bridge-regroup f1` is byte-identical to the f1v2 run's on all five f (R = the bridge = True rows:
  0 / 4 / 4 / 5 / 0 junctions; R2 = keep: 0 / 4 / 1 / 1 / 0) and its bridge rows equal the f1v2 `bridges.tsv` rows; H: the frozen R (179 cuts) and R2 (38
  cuts) equal the rows of the H `bridges.tsv` (all / keep) exactly.
- **Still to run before any score:** G-S / G-H (scorers), G-W, G5', G-D, G-U on the new units GTFs. **Order still to run:** scorers (S lost members,
  consolidated-locus partners, relation records for detector arms), NAT0 / A1 / A1s / N1 builds, the S plain-GTF families products for W, then Part A runs.

**Amendment 2 (2026-09-30 13:16; Part A is run and scored; E\* is fixed here; BEFORE any Part B arm was built or scored).** File sha1 before this
amendment `2e72e87458dfbd909d0aeceae77beda344abbb2c` (byte copy `notes/PREREG.pre_amendment2.md`; it already contains Amendment 1). Tools in
scratch `lib/`: `units2.py` **8109249e** (one change since Amendment 1: `read_cuts` defaults `oracle_label` to `label` when the column is absent, which
touches only the last column of `units.tsv`; no cut, regroup or name changes; 9 tests pass), `score_s2.py` d3b15444, `score_h2.py` 2eacc2fe, `gate_h.py`
a0b50695, `verdictA.py` ef375cdb, `arms_spec.py` dafed6d0, `build_partA.py` e345bc4f, `run_arm.sh` a342640f, `run_queue.sh` 8610788a, `run_g5.sh` 75476569,
`w_calls.py` bdba582a; imported UNCHANGED from the mechanism test: `score_s.py` 9708f3f9, `score_h.py` a286d14d, `units.py` 3f66eab0, `g6.py`.
- **Gates.** G-S pass (0 mismatches of every compared field for BASE and A at the five f, 9 arms); G-H pass (0 mismatches of 110 fields for BASE, A, N1 of
  the mechanism products); G-U pass (`g6.py` bad 0 on every new units GTF of S and H; rebuilds byte-identical); G5' pass (S f = .5 A1 and H A1: replay
  without `--dump-graph`, clusters.tsv and loci.tsv byte-identical); G-W pass (the glue reproduces the stored chr16 `calls['WA']` exactly: 929 transcripts, 4,974
  cuts); `unknown_locus` = 0 on every S scan. The S products of arms PLAIN / NAT0 / A1 / N1 and H NAT0 / A1 / A2 / N1-N3 are in `S/`, `H/`; aliases in
  `alias_of.tsv`.
- **What the execution arms are (facts, not yet interpreted).** A1s is byte-identical to A1 on every S f and on H (attaching ORIGINAL single-exon
  transcripts moves none: the polished plain GTFs hold no single-exon transcript that overlaps a same-strand multi-exon one, the mono-shadow polish having
  removed them; 93 of H's 122 single-exon UNITS are attached): A1s is a no-op and appears in the tables as an alias. A2 on S is byte-identical to A1 (no residual
  link: 10 of 10 cuts separated at every f); on H the residual links after E_nat are 21 cut pairs (5 single-exon units attached to the neighbouring unit's
  component, 16 in components of 8-81 transcripts, 4 of them with an articulation transcript), A2 separates 5 more, and the family products are unchanged.
  **Separation definitions.** The mechanism report's "215 of 253 oracle transcripts, 78 of 89 fused loci" are "units in >= 2 distinct gene_ids" and "a locus with
  >= 1 separated transcript"; this test's strict versions (every adjacent pair / every transcript) give 214 and 76 for A. Both are tabled for every arm and QA1 is
  evaluated under the mechanism's definition (the one the prereg text quoted). **E\* selection (§2.5): A1** (lost members A1 0, A1s 0, A2 0, A 2 (A loses one
  gene at S f = .5 and f = .9 by L_S); sum M1(S) + c1(H) 61 / 61 / 61 / 56; Compara TP pairs 109 each; order A1 first). **Part B runs under A1 (E_nat).**
- **Part B detector lists so far:** R and R2 as gated (S: R 0 / 4 / 4 / 5 / 0 junctions, R2 0 / 4 / 1 / 1 / 0; H: R 179 cuts on 170 transcripts, R2 38 on 38, frozen
  calls equal to `bridges.tsv`); W on S computed with the frozen instruments on the plain-GTF families products: 132 / 137 / 137 / 133 / 108 transcripts and 601 / 657 /
  655 / 635 / 585 cuts at f = 0 / .1 / .5 / .9 / 1 (H: 929 transcripts, 4,974 cuts). No Part B arm exists yet.

**Amendment 3 (2026-09-30 13:44; Part B is RUN and SCORED; this amendment adds two POST HOC, UNBARRED extensions designed AFTER seeing the Part B scores;
no bar, arm, threshold or verdict of §7 changes).** File sha1 before this amendment `5962f137be188490b091dc22360ba32c67d13d13`. Part B arms ran as registered: detector lists by
`detectors.py` c973e3ee, arms by `build_partB.py` d28f410b (E_nat, arm directories `d<NAME>`), scoring by `score_s2.py` / `score_h2.py` (Amendment 2), verdicts by
`verdictB.py` 26d0ebb7 (§7.2 / §7.3 literally), hub / single-exon / gate statistics by `safety2.py` a8e46bc5 (the mechanism's `safety1.py` with this test's arm paths). What
was seen when Amendment 3 was written: every Part B row (S and H), the census tables, the safety tables; in words: R2, R and R|W are at BASE level on H and gain only on S at
f = .5 / .9 (R, R|W), W and R union W lose members on S and H, and the literal S verdict of R2 / R / R|W is NOT because M2 at f = 1 is 4 against BASE's 3, which is the
native regroup alone (NAT0 has the same 4 with no cut). **Two post hoc extensions (reported, never barred; `units2.py` 9fea69a1 adds `--regroup scoped`; `build_partX.py` d0616844):**
- **X1 scoped execution.** The native rule is applied ONLY inside the input gene_ids that hold a unit; every other gene_id keeps its RG3 exon-overlap pieces (the shipped
  F1v2 / RG3 partition); a single-exon unit attached to a transcript of an untouched gene takes that piece's name. Arms: ORACLE (S f = .1 .. 1, H), R (S f = .1 / .5 / .9,
  H), R2 (S f = .5 / .9, H), R|W (H); at f = 1 R and R2 have no cut and their scoped input is the RG3-regrouped plain GTF (= BASE's input there). The question: what should a port do
  with the loci it does not cut (the native regroup alone, NAT0, is neutral-to-harmful against BASE: S lost members 3 / 1 / 1 at f = .1 / .5 / .9, H Compara pairs 87 against 99)?
  Read as information: SAFE iff no member lost (S L_S, H Compara and NPIP copies), M2 <= BASE's at every f and H c1 / pairs >= BASE's.
- **X2 NULL_R under E_nat.** R with each call replaced by one random intron of the same transcript (seed `20260930\tR`, draws with replacement), S f = .1 / .5 / .9 and H: whether
  R's S gains at f = .5 / .9 are the junctions (they are not attributable to junctions if this null gains as much).

**Amendment 4 (2026-09-30 14:21; written AFTER the Outcome below was drafted and placed here with the Amendments; it changes no result, bar or verdict).** File sha1 before this
amendment `ed81e748d0faee37210d0d82e929480c233840c5`.
- **Independent verification of E_nat** (scratch `verify/`; a separate agent re-implemented §2.1-§2.2 from the prereg text alone, reusing only the mechanism test's `parse`, `introns`, `split_units`, and never opening
  `units2.py`, `gate_native.py` or `test_units2.py`): the families input equals A1's (S f = .5: 646 transcripts; H: 9,729); the **PARTITION into gene_ids is identical on S f = .5 (194 classes) and on H (3,003 classes)**; oracle
  transcripts fully separated 10 of 10 (S) and 232 of 253 (H), with >= 2 distinct gene_ids 10 and 233, as the implementation reports. Names are identical on S; on H 143 of 9,729 differ under the literal step order (the representative
  fixed at step 2, BEFORE attachment) and 0 under the reading `units2.py` implements (the representative is taken over the component AFTER attachment: an attached single-exon unit carries its parent's reads and outranks the step-1
  representative in 15 of the 40 components that receive units). **§2.2 is clarified: the representative used for the names (steps 2 and 4) is computed over the component after attachment.** No partition, family, score or verdict
  depends on names.
- **Other ambiguities the verifier named, all resolved as implemented:** attachment targets include multi-exon units; "shares the most exonic bases" is per transcript (per-component totals would change the partition); ties by the lower
  families-input index (35 of 93 attachments tie at the maximum, 14 of them across components: breaking the ties the other way changes 17 classes, so the rule is load-bearing; it is exercised on H only, S having no single-exon unit); the first sentence of step 4 is garbled ("the gene_id of the representative's INPUT gene_id": read as the representative's input gene_id, a
  unit's being its parent's; naming by the base tid instead differs in 724-756 names on H); "(the input gene_ids are discarded)" means not used for grouping; junction keys take donor and acceptor in genomic order on both strands.
- **Consolidated locus for groups that mix input gene_ids** (24 of H's 3,003 A1 genes, through attachment): mapped to the input gene_id of the group's representative; the alternative (the units' parent gene) changes no NPIP count
  (partners by locus and members by locus identical for A1, A2, R, R2, R|W, W, the scoped oracle and scoped R).
- **Order.** Amendment 3 was written after the post hoc X-arm INPUTS (GTFs of the scoped and NULL_R arms) had been built and before any X arm was run through the families stage or scored.
- **Tools added after Amendment 3:** `verdictX.py` 0e31e8bc, `relations.py` d90a33da (emitted for A1 S f = .5 / f = 1, A1 H and the scoped R arms S f = .5, H); `run_g5.sh` was also run on S f = .1 R2 and H R2 (byte-identical).

## Outcome

**Outcome (2026-09-30 14:15; Amendments 1-4 are the only deviations; no bar, arm, threshold or truth changed after a product was seen; scored with `score_s2.py` d3b15444, `score_h2.py` 2eacc2fe,
bars by `verdictA.py` ef375cdb / `verdictB.py` 26d0ebb7 / `verdictX.py`).** Scratch `/mnt/linuxdisk/tmp/rustle_figures_dev/container_units_v2/` (`results/*.json`, `partA_tables.md`, `partB_tables.md`, `a3/`); report
`scratchpad/figs/container_units_v2.md`.

**Gates.** G0, G-A1a (exact names on six unpolished GTFs), G-A1b (0 merges), G-U (+ byte-identical rebuilds and `--regroup rg3` reproducing the mechanism's A / A0 / N GTFs), G5' (S f = .5 A1, H A1, S f = .1 R2, H R2), G-S (0 mismatches), G-H (0 of 110
fields), G-W (929 transcripts, 4,974 cuts reproduced), G-R, GD (8,147 = 8,147) all pass. **Part A (oracle junctions, execution E_nat = A1).**

| S f = .1 / .5 / .9 / 1 | BASE | A (RG3 units) | NAT0 | **A1** | N1 (null) |
|---|---|---|---|---|---|
| fused copies in NPIP /10 | 2 / 0 / 0 / 3 | 7 / 7 / 8 / 9 | 1 / 0 / 0 / 4 | **9 / 9 / 9 / 9** | 1 / 1 / 1 / 5 |
| partners in NPIP /20 | 0 / 0 / 0 / 3 | 0 / 0 / 0 / 0 | 1 / 0 / 0 / 4 | **0 / 0 / 0 / 0** | 1 / 0 / 0 / 3 |
| unfused copies /15 | 13 / 13 / 13 / 15 | 14 / 14 / 15 / 15 | 14 / 14 / 14 / 15 | 15 / 15 / 14 / 15 | 14 / 15 / 15 / 11 |
| oracle transcripts separated /10 | | 9 / 9 / 9 / 10 | | **10 / 10 / 10 / 10** | 0 / 1 / 1 / 5 |
| relation P / R lenient (strict) | | 1 / 1 (.8 .8 .9 1) | | 1 / 1 (1 / 1 at every f) | .4 .5 .5 .7 |
| members lost (L_S genes) | | 0 / 1 / 1 / 0 | 3 / 1 / 1 / 0 | **0 / 0 / 0 / 0** | 3 / 2 / 2 / 0 |

S: A1 **WORKS** by the mechanism's absolute bars (S1-S4 at every f; A: PARTIAL). NPIPB12 (third transcript) separates at every f, NPIPB13's knife-edge edge lands in NPIP at f = .1 / .5, NPIPA7 is the one fused copy outside NPIP (outside at f = 0 already). At f >= .9 a fused copy's locus is represented by its
unit (truncated by a terminal exon). H (chr16): A1 = A on every family metric (25 / 26 copies incl. NPIPB5, F .667 -> .699, 109 pairs, Liftoff 16, referee .250, Soto-NPIP .800 / .700, 0 of 49 lost, none of the BASE copies lost), with 233 of 253 oracle transcripts and 83 of 89 fused loci separated (A: 215 / 78); H2 **fails by unit (4 > 3) and holds
by consolidated locus (2 <= 3)**; six fused loci stay fused (NPIPA1 | PKD1P3, MOSMO | VWA3A, FRG2KP | YBX3P1, LOC105379535 | LOC107987382, HAS3 | TANGO6, KLHDC4 | LOC100129215); A2 (articulation transcripts; needed on H only: 21 residual pairs) separates 5 more pairs and changes no NPIP / Compara metric; A1s == A1 (a no-op on the polished GTFs); the
null under the same execution gains nothing (copies 24, F .650-.683, 1 member lost). **E\* = A1.** QA1, QA2, QA3 hold. Families time 43-61 s (BASE 49-68 s), 2.0-2.5 GB. **A3:** the output spec (relations + members-by-locus + v1 residual) is in the report; `relations.py` emits it for S f = .5, S f = 1 and H A1 with every invariant asserted (S f = .5: 10 relation rows, 9 `cover`, 138 members, 0 double counts; H: 253 rows, 100 `cover`, 608 members, 8 double counts; v1 residual 179 / 177 / 190 accessory blocks).

**Part B (detectors through E_nat; literal §7.2 / §7.3).**

| detector | cuts S (f = 0 .. 1) / H | S fused copies f = .1 / .5 / .9 / 1 (BASE 2 / 0 / 0 / 3; oracle 9) | S members lost | H copies (BASE 24; oracle 25) / pairs (99; 109) / Compara members lost | S verdict | H verdict | overall |
|---|---|---|---|---|---|---|---|
| R2 | 0 / 4 / 1 / 1 / 0; 38 | 2 / 0 / 0 / 4 | 0 | 24 / 99 / 0 | NOT (M2 4 > 3 at f = 1 = NAT0) | NOT (no gain) | NOT |
| R | 0 / 4 / 4 / 5 / 0; 179 | 2 / 2 / 4 / 4 | 0 | 24 / 99 / 0 | NOT (same) | NOT (no gain) | NOT |
| R\|W (post hoc) | = R on S; 29 | 2 / 2 / 4 / 4 | 0 | 24 / 99 / 0 | NOT (same) | NOT (no gain) | NOT |
| R union W | 601 / 657 / 655 / 636 / 585; 5,139 | 4 / 3 / 4 / 4 | 2 / 4 / 5 / 5 / 5 | 17 (7 copies lost) / 96 / 1 | NOT | NOT | NOT |
| W | 601 / 657 / 655 / 635 / 585; 4,974 | 4 / 3 / 3 / 4 | 2 / 4 / 5 / 5 / 5 | 17 (same seven lost) / 96 / 1 | NOT | NOT | NOT |
| NULL_W (control) | 480-526; 3,762 | 2 / 2 / 5 / 6 | 2-4 | 20 / 102 / 1 | | | clears no bar |

**Decision rule (§7.3): 1. R|W, 2. R, 3. R2, 4. R union W, 5. W;** the first three are equal on every outcome (c1 24, pairs 99, 0 lost) and ordered by S copies sum (67 / 67 / 62) and by H cuts (29 / 179 / 38); no detector WORKS on both substrates; the conditional nulls were not triggered (no verdict above NOT). The matched null (NULL_W) clears no bar.
**False-split census (H):** R 160 of 179 cuts outside the oracle, 134 transcripts inside one gene (13 genes), families ALL_UNCL 131 / DIFF 30 (6 annotation-supported, 12 single-gene, 12 other: 24 are the isoforms of one gene, SMG1P1); R|W 21 of 29 outside, 9 transcripts in 1 gene, DIFF 6 (6 of 6 supported); W 4,882 of 4,974 outside, 629 transcripts in 90 genes, DIFF 352
(4 supported), 83% of its 5,903 units single-exon; on S W cuts 103 single-gene transcripts at f = 0. **The four readthrough copies:** unchanged for R2 / R / R|W; NPIPA9 leaves NPIP for W and R union W; PKD1P6-NPIPP1 stays outside NPIP in every arm (MCL24-MCL40).

**Post hoc, not barred (Amendment 3).** X1 scoped execution: oracle = A1; R, R2, R|W **SAFE** on S and H (0 lost, M2 <= BASE at every f incl. f = 1, M1u >= BASE, H identical to BASE); R gains +2 / +4 fused copies at f = .5 / .9 (M1 2 / 2 / 4 / 3). X2 NULL_R: copies 1 / 0 / 0, 3 / 1 / 1 genes lost, H pairs 87 and one member lost: R's gains are the junctions.
Reading: the F1 share rule is not needed once bridges are units (R precise on S, member-safe everywhere); the native regroup must be SCOPED to the gene_ids it cuts (NAT0 alone: S M2 4 at f = 1, H pairs 87 = A0, i.e. it lacks F1v2's bridge removal).

**Scorecard.** Part A: P-A1, P-A2, P-A4, P-A8, P-A9 hit; P-A3 hit except the exact vector (9 / 9 / 9 / 9, not 8 / 8 / 9 / 9); P-A5 miss (83 / 89, 233 / 253); P-A6 hit (Compara F equals A's .6988, a literal miss by .0002 of the stated .699); P-A7 half (A1s identical on S as predicted, on H too). Part B: P-B1 hit as a claim, both sub-claims missed (R keeps 44% at f = .9; R2 gains 0 at f = .1); P-B2, P-B3, P-B6 hit; P-B4 miss
(R loses no member; R, R2, R|W are NOT by no gain); P-B5 half (none WORKS on both; the first-ranked detector is R|W; three detectors lose no member). Falsifiers Z1-Z6 not fired.

**Hostile review of this Outcome.** (a) The S gain of A1 over RG3 rests on NPIPB12 and one flipped knife-edge edge; on H it is nil. (b) The literal Part B verdicts are partly rule artefacts (S M2 at f = 1 is NAT0's; H "no gain = NOT"): quoted literally, with the post hoc neutral readings labelled. (c) L_S is co-membership against this pipeline's own BASE at f = 0. (d) R and R2 carry no fresh evidence; R|W is post hoc and dev-selected; W was not tuned here. (e) M1 is holder based: at f >= .9 the member is a unit. (f) The census counts per transcript: 24 of R's 30 DIFF records are one gene.
(g) Dev only (S one simulation, H chr16); no held-out substrate, PTR or PPY product was opened; nothing here supports a default. (h) Disclosure: `truth.human_A119b.pkl` (the definition study's annotation-only truth table) holds every A119b contig; this test's census reads only its chr16 entries, but an exploratory look at its structure printed one chr1 record and the class counts over all contigs (annotation labels, no read, no product, no score; the parent Outcome had printed those counts).

**Recommendation (Part C, detail in the report).** Port the EXECUTION and the output spec behind an opt-in flag (`--bridge-regroup f1units`: R's bridge junctions without the share rule, units in the families input, scoped native regroup, relations + members-by-locus); do not port W; R|W (W-gated R: 6 of 6 relation records supported, 3 junctions cut) is the next pre-registration on the reserved substrates. Expected dev gain of the port: +2 / +4 fused copies on S at f = .5 / .9, none on chr16; the oracle shows the headroom (S 9 / 9 / 9 / 9, chr16 +1 copy and +10 Compara pairs).
