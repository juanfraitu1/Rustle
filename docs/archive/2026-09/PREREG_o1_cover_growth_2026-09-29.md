# Pre-registration: COVER-BY-GROWTH families for O1 (fusion loci may belong to two families), held out genome-wide

**Written 2026-09-29 (KEY=o1_cover) before any held-out COVER product exists.** The rule was designed and every
instrument was exercised on development data only (§12): human A119b chr16 (the F1v2 dev GTFs) and the gorilla fusion
simulation (the F1v2 simulation arms f0.0-f1.0). This file binds once Amendment 1 records its sha1 and the frozen
instruments' sha1s. Until then nothing held-out runs. Nothing in `src/` is edited; nothing is committed. The rule is a
Python post-processor over the families stage; a Rust port is a later decision, and so is any default change.

## 0. The question, and what is already known before the freeze

O1 emits a strict partition (0 / 2,670 loci in more than one cluster, ledger §6s9), so a real fusion such as
PKD1P6-NPIPP1 (110 MAPQ-60 reads over a canonical junction; memory `project_pkd1p6_npipp1_is_real`) cannot be a member of
both parent families. Node cutting (r846) and multi-label truths (r845) are closed and are not re-proposed. Soto 2025's
released algorithm, reproduced at ARI 0.97 today (`PREREG_soto_reconciliation_2026-09-29.md`), grows families ONLY
through core genes over pair-level edges; peripheral genes attach to every family they touch without extending it,
which yields their cover (481 / 491 exact cover-aware). F1v2 (`PREREG_f1v2_readshare_2026-09-29.md`, EFFECTIVE on human)
already turns minority readthrough links into explicit `fusion_of` relation records, and its verification showed that
its fused-locus gain over RG3 exists only if those bridges are kept as relations. **Question:** if the F1v2 bridges are
made PERIPHERAL nodes (families grow through the core only; a bridge joins every family it touches), do fusion loci land
in both parent families, without a partition regression and without spurious multi-memberships, on data where the
rule was never examined?

**Known before the freeze (disclosed; it caps this test).** A truth-side power count (`lib/power.py`, `power_lo.py`: the
F1v2 bridges, their parents = the `gene_at` labels of their `fusion_of` pieces, the Compara universe and Liftoff
multi-copy genes; NO family product of any held-out arm read) gives:

| substrate | bridges on the scored contigs | ≥ 2 distinct parents | ≥ 1 parent in Compara | parents in ≥ 2 Compara families | ≥ 1 / ≥ 2 Liftoff multi-copy parents |
|---|---|---|---|---|---|
| human A119b (scored contigs §3) | 348 | 193 | 10 | **0** | 0 / **0** |
| human testis (minus chr16, chr18) | 12 | 7 | 2 | **0** | 0 / **0** |
| gorilla KB3781 (minus the 3 NPIP contigs) | 52 | 41 | – | – | 0 / **0** |
| gorilla OR6737 (minus the 3 NPIP contigs) | 56 | 48 | – | – | 0 / **0** |
| human A119b chr16 (dev) | 12 | 9 | 2 | 0 | – |

(Liftoff multi-copy = the record has ≥ 1 extra copy at sequence_ID ≥ 0.95 with both exon unions ≥ 200 bp: 722 human,
392 gorilla genes.)

So **no held-out bridge joins two multi-copy families**: F1v2's minority bridges are readthroughs between single-copy
neighbours (the 10 A119b bridges with a Compara parent are ZNF / NBPF / GAGE-type fusions whose second parent is a
lncRNA, a single-copy gene or the same family). The motivating fusions (PKD1P6-NPIPP1, PDXDC2P-NPIPB14P, PKD1P3-NPIPA1)
are DOMINANT links that F1v2 by design never makes bridges (and on A119b chr16 the PKD1P loci sit in the NPIP family
itself, §12). **The two-family gain (C1) therefore cannot be judged on any held-out substrate; it is shown only on the
simulation (development).** The held-out part of this test is a SAFETY test: does the cover damage the partition,
over-cover, or attach bridges to families of neither parent? EFFECTIVE is unreachable held-out by construction; the best
possible held-out outcome is SAFE-INERT (§7).

## 1. The rule COVER (binding; `cover.py`, sha1 in Amendment 1)

### 1.1 Verbatim from the frozen file's docstring (the RULE block)

```
RULE (binding once pre-registered; docs/archive/2026-09/PREREG_o1_cover_growth_2026-09-29.md §1)
  Input: an F1v2 GTF in which every bridge transcript group is its own gene_id and carries `fusion_of`
  (F1v2.gtf), and the shipped families stage run on that GTF (`mcl_families --from-gtf`, driver flags), i.e. its
  `loci.paf`. The graph G = the shipped edge graph of that run (emu.py fac9a560 = the binary, gated byte for byte).
  PERIPHERAL set P (arm COVER) = the loci (gene_id groups) with >= 1 transcript carrying `fusion_of` (F1v2 bridge
  loci). CORE = every other locus.
  1. Core families = the shipped clustering (emu.mcl + fold_and_cluster = MCL I 2.8, prune 1e-9, fold-within-
     clusters, families of >= 2) run on G[CORE], the subgraph induced by the core loci: the PAF records with a
     peripheral end are left out before the edge step, nothing else changes.
  2. Each peripheral locus p attaches to EVERY core family that contains a core locus q with an admitted edge p-q
     in G (the same edge admission, emu.build_graph on the same PAF records). Attachment never merges families,
     never changes a core family, never makes a new family; a peripheral locus with no such q is in no family.
     Peripheral-peripheral edges are ignored. Peripheral loci are never folded.
  No constant beyond the shipped ones.
```

"Contains a core locus q" counts q whether it is a kept representative of the family or a record folded into one
(fold-within-clusters). A peripheral span key never equals a core span key (asserted). The bridges are exactly F1v2's
(`f1v2.py` b4e788ad `--rule min` over the frozen F1 `f1_bridge.py` 37ee8e77): the `<g>.fus<k>` gene_ids of `F1v2.gtf`.

### 1.2 Constants and sources

| element | value | source |
|---|---|---|
| edge admission, MCL, fold, family size | the shipped families stage (`--min-exonic-bp 1 --min-shared-exon-frac 0.60`, I 2.8, prune 1e-9, ≥ 2) | `docs/seeded_family_definition.md` §0★; emulated by `emu.py` fac9a560 (byte-exact, gate G1) |
| peripheral set | F1v2 bridge loci | F1v2 (EFFECTIVE, r1148-1151 drafts) |
| attachment | every family with ≥ 1 admitted core neighbour | Soto's released rule (a non-coding gene joins every family it has a kept pair with), transposed to our edge and our peripheral class |
| NULL matching strata | (scored contig or not) × degree in G (0, 1, 2, ≥ 3) | fixed here; degree governs how many families a node can touch |

Why no threshold: attachment is existential over the shipped edges, and the core/peripheral split is F1v2's own
bridge predicate. The design mirrors Soto (growth through core, attachment of periphery); only the peripheral class
differs (Soto: non-coding / processed genes by biotype; here: read-proven minority fusion links).

## 2. Arms (all from ONE families run per substrate on `F1v2.gtf`, plus the stored BASE run)

| arm | what | role |
|---|---|---|
| **BASE** | the shipped families on the BASE GTF (stored run; emu-gated, G1b) | comparator (the shipped default) |
| **ALL** | the binary's partition of `F1v2.gtf` (bridges as ordinary loci) | reported (the partition alternative for bridges) |
| **F1v2** | = COVER's core partition (bridges removed from the graph). On dev the binary run on `F1v2.families.gtf` equals it byte for byte (§12) | comparator (partition) |
| **COVER** | core partition + attachments of the bridges | treated |
| **NULL0-4** | the same construction with P = \|P\| random non-bridge loci drawn per (scored flag, degree) stratum, seeds `o1cover:<label>:s<k>` (per chromosome on A119b, label `human_A119b:<chrom>`); bridges are then ordinary core loci | matched null |
| **SE** (declared a priori; REPORTED ONLY) | P = bridges ∪ loci whose representative has ONE exon (the de novo analogue of Soto's processed-pseudogene / non-coding periphery) | reported, never judged. Prediction: it removes intronless coding families (ORs, histones) from the core, so its core sensitivity falls; dev chr16: 703 of 2,844 loci peripheral, Compara pairs 99 → 94 |

## 3. Substrates and exposure

Held-out = data where the COVER rule was never examined; each substrate is scored WITHOUT its NPIP-containing contigs
(the NPIP block is where every earlier fusion rule was designed).

- **V_A = human A119b**, per-chromosome families on chr1-15, chr17, chr19, chrX, chrY, chrM (excluded: chr16, chr18 =
  NPIP-containing; chr20, chr21, chr22 = F1v2 dev / excluded). **Per-chromosome, not genome-wide (declared deviation,
  cost):** a genome-wide A119b families run is ≈ 4-6 h per arm (93,532 loci, 1.98 Gb of spans; gorilla OR's 1.01 Gb took
  100 min). Per chromosome, BASE_c and ALL_c are the shipped binary on the chromosome's GTF lines (the product type of the
  chr16 dev runs); cross-chromosome edges are absent in every arm alike. G1 per chromosome.
- **V_T = human testis**, genome-wide families (BASE stored: `rt_arms/human_testis/human_testis.BASE.fam`), scored minus
  chr16, chr18.
- **V_K = gorilla KB3781**, **V_O = gorilla OR6737**, genome-wide families (BASE stored in `rt_arms`), scored minus
  NC_073242.2, NC_073241.2, NC_073244.2 (the three contigs holding a T_member NPIP copy; NC_073244.2 is also F1's dev
  contig).
- GTFs: BASE = `runs/<s>/<s>.gtf` (A119b 631c9f11, testis 50f239d6, the F1v2 test's BASE); F1v2 human = the F1v2 test's
  `F1v2.gtf` (A119b 976b93dc, testis 4afc60a2); F1v2 gorilla = frozen `f1v2.py --rule min` on (BASE, the F1 test's
  `F1.junctions.tsv`), regenerated today and byte-identical to the F1v2 design's `dev/<s>/min.gtf` (OR 40b1e5d7, KB
  ad177315).
- Species are never pooled; each substrate is reported alone.

**Exposure (plainly).** Every held-out sample is spent for BRIDGE work: the gorilla samples were F1's held-out and
F1v2's development substrates (the minority condition was designed on them); A119b and testis were F1v2's held-out.
The family-level COVER rule was never examined on any of them. What this author read of held-out data before this file:
the power counts of §0 (bridge parents × truth tables; no family product); the F1v2 Outcome; the npf_audit ape rows
(gorilla NPIP holders; the NPIP contigs are excluded from every family metric here and enter only the non-blind C6
control); file sizes, loci counts and span totals of the BASE GTFs; the stored BASE families logs of OR / KB (graph
size, cluster count). No held-out F1v2 families run, COVER, NULL or score was computed or opened.

## 4. Scoring (instruments `score.py`, `run_score.py`, `lo_score.py`, `npip_ctrl.py`, `clauses.py`)

**4.1 Gene space = family_score `--chrom ALL` semantics, re-derived in Python and gated (G2) against the frozen
`family_score` 7723029b binary on every core arm.** Gene = (contig, Name) of gene / pseudogene / ncRNA_gene records of
`families_gw/species/<sp>/genes_only.gff`; locus → one gene by `gene_at` (largest span overlap, first maximum); truth
families with ≥ 2 genes present; universe U; bipartite (scipy tie policy), pairwise over U; contigs outside the substrate
removed from clusters and truth first (the figures' `fs_score`). **Truth (human): Ensembl Compara Primates**
(`compara.Primates.families.tsv`, 426 families / 1,313 genes genome-wide; the figures' main human family reference).
Soto is not used (not independent; CHM13 v1.0 genes).

**4.2 Liftoff copy pairs (both species)**, the figures' `copy_pairs` / `pair_families` unchanged: (record, extra copy)
pairs with sequence_ID ≥ 0.95, both exon unions ≥ 200 bp, both read-supported (≥ 2 reads) in the sample, both on scored
contigs; recovered when one family's loci cover ≥ 50% of each locus's exon union. Family loci = each clusters.tsv
member's representative exon union; cover arms add each attachment. Pairwise recall only (a relation, not a partition).

**4.3 Cover-aware gene space** (arms with peripheral members). Items are genes, so a peripheral p contributes ONE label
to each family X it joins: (1) a parent of p already in X's core genes, else (2) a parent in U whose truth family is X's
plurality truth family, else (3) `gene_at(p)` (a wrong member when in U). **A peripheral member counts as correct iff the
truth lists one of its parents in that family.** Parents(p) = the distinct `gene_at` labels of p's `fusion_of` pieces
(a bridge) or of p itself (a NULL / SE locus). Core-only metrics are always reported beside, so the cover cannot hide a
partition regression.

**4.4 Attachment judge** (every (p, X) with X holding ≥ 1 member on a scored contig; *judged* iff X's core loci carry ≥ 1
gene label): **own** = a parent's own locus is in X (a core locus of X, kept or folded, is labelled with a parent);
**truth** (human) = X's core genes include a Compara family-mate of a parent; **lo** = one of X's members covers ≥ 50% of
a Liftoff copy-pair locus (record or extra copy, sequence_ID ≥ 0.95) of a parent. *Supported* = own ∨ truth ∨ lo.
**SELF** attachment = X holds one of p's own `fusion_of` pieces; **HOMOLOGY** = every other attachment. A SELF
attachment is supported by construction (the piece is labelled with the parent): disclosed, which is why C3 is judged on
HOMOLOGY attachments only.

**4.5 "Fusion member of both" (C1).** Eligible (truth side, fixed by the truth, not by a prediction): a bridge with two
parents q1, q2 in U whose Compara families differ. both(b) ⟺ b is attached to two distinct families X1, X2 with q1 in
X1 and q2 in X2, "q in X" meaning q's own locus is in X or X's core genes include a Compara family-mate of q.
**A fusion is correct in family X iff one of its parents is in X.**

**4.6 NPIP control (C6)**: copies = the 26 Dishuck chr16 copies (human) / the 25 T_member copies (gorilla); holder = the
arm's core locus with the most same-strand exon bp over the copy (ties: reads); NPIP family = the family holding the
most copies; under COVER a copy also counts when a bridge overlapping it is attached to the NPIP family.

## 5. Gates (a failed gate stops its substrate; fixed in the instrument only)

| gate | what must hold |
|---|---|
| G0 | the frozen `SHA1SUMS` hold; input GTF sha1s as §3 |
| G1 | `emu.py` R0 on the ALL run's PAF reproduces the binary's `clusters.tsv` byte for byte (per chromosome on A119b) |
| G1b | `emu.py` R0 on the stored BASE PAF reproduces the stored BASE `clusters.tsv` byte for byte |
| G2 | the Python gene-space core numbers equal the `family_score` binary's printed pooled and pairwise numbers, every human core arm |

## 6. Clauses (per held-out substrate; integers unless stated)

| clause | passes iff | states |
|---|---|---|
| **C1** gain: fusions in both parent families (human) | eligible E ≥ 1 (§4.5); PASS iff 2·both > E | NJ if E = 0 (**known: E = 0 on V_A and V_T, §0**); NA on gorilla (no family truth) |
| **C2** specificity of multi-membership | judged iff COVER has ≥ 1 peripheral locus attached to ≥ 2 families (attachment rows, §4.4). (a) multi-rate(COVER) > multi-rate(NULL_k) for every k (rates = multi / peripheral loci on the substrate, integer cross-multiplication); (b) among the judged attachments of COVER's multi-family loci, supported > unsupported | PASS = a∧b; FAIL = judged and ¬(a∧b); NJ |
| **C3** homology attachments correct (human) | judged iff ≥ 1 judged HOMOLOGY attachment; supported > unsupported | PASS / FAIL / NJ; on gorilla REPORTED only (no family truth; own ∨ lo is conservative: the sim's correct NPIPA7 → NPIP attachment is unsupported, §12) |
| **C4** no partition regression vs BASE | COVER's core (= F1v2) ≥ BASE − 0.01 on: human — Compara bipartite sensitivity, precision, F, pairwise sensitivity, pairwise precision, and Liftoff pair recall; gorilla — Liftoff pair recall | PASS / FAIL. 0.01 = half the house gain bar ΔF ≥ 0.02 (r845 CP-3), fixed here |
| **C5** the cover never lowers the score | human: cover-aware F ≥ core F and cover-aware pairwise precision ≥ core pairwise precision (zero tolerance); all: cover-aware Liftoff recall ≥ core recall | PASS / FAIL |
| **C6** NPIP no-regression (non-blind control) | copies in the NPIP family under COVER ≥ under BASE (testis: Dishuck 26; gorilla: T_member 25) | PASS / FAIL; NA on A119b (chr16 not computed) |

Reported, never judged: every arm's core and cover-aware numbers (ALL, F1v2, SE, NULL0-4), attachments per bridge with
SELF / HOMOLOGY class and judge flags, the attachment support of the NULL loci, multi-family bridges by name, the SE arm.

## 7. Verdict (nothing pooled)

- **REFUTE** iff on any held-out substrate C2, C3, C4, C5 or C6 = FAIL (the cover over-covers, attaches to wrong
  families, hides or causes a partition regression, lowers the score, or loses an NPIP copy).
- **EFFECTIVE** iff C1 = PASS on ≥ 1 held-out substrate and no refute trigger (**unreachable here: C1 is NJ / NA
  everywhere, §0**).
- **SAFE-INERT (KEEP OPT-IN; gain untested held-out)** iff no refute trigger and C1 is NJ / NA on every substrate.
- KEEP OPT-IN otherwise (C1 judged and failing without a refute trigger cannot occur here).

The gain is claimed only from the simulation (§12), labelled development evidence. A default change or a Rust port is
the user's call and would need a held-out gain substrate (§10).

## 8. Predictions (this author's probabilities, before any held-out number)

1. Gates G0, G1, G1b, G2 pass on every substrate (0.9).
2. COVER has ≥ 1 multi-family bridge: A119b 0.4, testis 0.1, KB 0.25, OR 0.25. C2 judged on ≥ 1 substrate (0.5); PASS
   where judged (0.6).
3. HOMOLOGY attachments exist on A119b (0.5); C3 PASS where judged (0.6).
4. C4 PASS: A119b 0.8, testis 0.85, KB 0.85, OR 0.85; COVER core Compara F ≥ BASE on A119b (0.6).
5. C5 PASS 0.9 per substrate; C6 PASS 0.9 per substrate.
6. Bridges attached to ≥ 1 family: A119b 10-40% (0.6); SELF ≥ HOMOLOGY attachments (0.8).
7. **Verdict:** SAFE-INERT 0.55, REFUTE 0.40, other 0.05.

## 9. Falsifiers of the design reasoning (reported whatever the verdict)

- **Z1 "bridges attach through their own pieces"**: falsified on a substrate where HOMOLOGY attachments ≥ SELF.
- **Z2 "minority bridges almost never join two families"**: falsified where multi-family bridges ≥ ⅒ of the bridges.
- **Z3 "removing bridges is no worse for the partition than removing random loci of the same degree"**: falsified where
  COVER-core Compara F (human) or Liftoff recall (gorilla) is below every NULL_k core.
- **Z4 "random peripheralization over-covers"**: reported: NULL multi-rate vs COVER (dev chr16: NULL 0-1 of 12, COVER 0).

## 10. Order, stop rules, machine rules, cost; not in this test

**Order** (driver `ho.py`, run from the frozen copy): G0 → testis (`fam human_testis ALL` until done → `gate1b` →
`cover` → `score`) → KB → OR → A119b (`split` → per chromosome `fam BASE`, `fam ALL`, `gate1b`, `cover` → `merge` →
`score`) → `clauses.py`. Tables to `o1_cover/held/tables/` before the next substrate starts.
**Stop rules.** After the freeze nothing changes in the rule, an arm, a null seed, a clause, a bar, a truth or a
substrate. A failed gate stops its substrate (fixed in the instrument only, as an amendment). A step that hits its time
cap is resumed (the shard wrapper's exit 75) or re-run once; a second failure makes that substrate's clauses not
measured. No variant (other peripheral classes, attachment thresholds, SE as a verdict arm) is substituted.
**Machine rules.** Families runs via `tools/rlock.sh heavy` (RLOCK_TIMEOUT 590; genome-wide runs through
`tools/mm2_shard.sh`, budget 420 s per call), Python steps via `light` (heavy when > 2 GB); foreground; `TMPDIR` under
`/mnt/linuxdisk`; never `pkill -f`; scratch `/mnt/linuxdisk/tmp/rustle_figures_dev/o1_cover/` (< 40 GB new data).
**Cost.** ALL families: testis ≈ 20 min, KB ≈ 65 min, OR ≈ 100 min (the BASE runs' times); A119b per chromosome ≈ 20
chromosomes × 2 arms × 1-5 min; cover + scoring ≈ 5-15 min per substrate. **≈ 4-5 h wall** plus lock waits.
**Not in this test:** a held-out GAIN substrate (needs a new simulation on a second genome, or a sample whose minority
bridges join two multi-copy families); dominant fusions (F1v2 never makes them bridges); SE as a rule; any Rust port.

## 11. Hostile self-review (fixes applied above)

1. **"The held-out test cannot show the gain."** Correct and stated first (§0): C1 is NJ on every held-out substrate by a
   truth-side count made before the freeze. The verdict is capped at SAFE-INERT; the gain is dev-only (simulation).
2. **"Self-overlap edges make bridge attachments 'correct' by construction."** Yes: a bridge's span contains its pieces,
   so the all-vs-all aligns them (identity ≈ 1.0) and the admitted SELF edge attaches it to its own pieces' families.
   The judge labels SELF attachments and C3 is judged on HOMOLOGY attachments only.
3. **"The NULL's multi-rate is biased low: random loci are not two-piece objects."** Intended: C2 asks whether
   multi-membership concentrates at bridges (two-gene objects) rather than arising wherever a node is made peripheral.
   Degree matching removes the first-order confound (a node with more neighbours touches more families).
4. **"F1v2 arm = COVER core is emulated, not the binary on `F1v2.families.gtf`."** On dev they are byte-identical (§12);
   the only possible difference is minimap2 target competition from the bridge sequences. Stated; the held-out F1v2 arm
   is the emulated core by definition.
5. **"Per-chromosome A119b is not the shipped genome-wide product."** Declared (§3), applied to every arm alike, and the
   other three substrates are genome-wide.
6. **"The gorilla judge is too weak."** It has no family truth; C3 is not a clause on gorilla (reported), and C2 (b) there
   relies on own ∨ Liftoff. Stated.
7. **"C4's tolerance favours the rule."** 0.01 is fixed here from the house gain bar (half of ΔF 0.02); every metric is
   reported, and C4 compares with BASE (the default), not with ALL.
8. **"Compara is protein-coding and bridges are mostly lncRNA-linked."** Yes; most attachments are unjudgeable by
   Compara (hence Liftoff and own). Reported as counts.
9. **"The cover-aware label choice is lenient."** It is the stated definition (a peripheral member is correct iff the
   truth lists one of its parents in that family) and it can only ADD one label per (p, X); rule (3) makes a wrong
   attachment a false member whenever `gene_at(p)` is in U. Core-only metrics are always beside it.
10. **"Metric traps."** Fixed universe per substrate (truth-side); no denominator conditioned on a prediction (C1's
    eligibility is truth-side; C2/C3 rates are over all attachments / peripheral loci); bipartite precision is
    tie-dependent (scipy policy, as family_score) and pairwise numbers are reported beside it; the SE arm's node
    removal is reported with its universe, never judged.

## 12. Development evidence (design phase; every number below was seen before this file)

Scratch `/mnt/linuxdisk/tmp/rustle_figures_dev/o1_cover/` (`dev/`, `lib/`); report `scratchpad/figs/o1_cover.md`.

**12.1 Human A119b chr16 (F1v2 dev GTFs `hdev/hsa16.{BASE,F1v2,F1v2.families}.gtf`; families = frozen
`fj_bin_frozen/mcl_families` 91ef2e1c, driver flags; 58-68 s, ≤ 2.5 GB each).** G1 (emu = binary on the ALL run) PASS;
the binary on `F1v2.families.gtf` = COVER's core, byte for byte; G2 PASS on all 10 arms.

| arm | Compara (chr16): bipartite sens / prec / F | pairs TP / pred (of 193) | NPIP copies in NPIP /26 | fused copies (npf_audit's 11) in NPIP |
|---|---|---|---|---|
| BASE | .444 / 1.0 / .615 | 66 / 66 | 22 | 8 |
| ALL | .500 / 1.0 / .667 | 99 / 99 | 24 | 10 |
| F1v2 = COVER core | .500 / 1.0 / .667 | 99 / 99 | 24 | 10 |
| COVER (cover-aware) | .500 / 1.0 / .667 | 99 / 99 | 24 (bridges add none) | 10 |
| NULL0-4 cores | .500 / 1.0 / .667 each | 99 / 99 | – | – |
| SE (reported) | .500 / 1.0 / .667 core; cover-aware 98 / 98 pairs | 94 / 94 | – | – |

12 bridge loci; 3 attached (1 family each: NPIPA7|LOC131696449 → NPIP MCL1, CLEC18A-region → MCL17, one at 63.9 Mb →
MCL108), all SELF and supported; **0 multi-family** (NULL 0-1 of 12). **PKD1P6-NPIPP1 is not in both families in any
arm:** it is not a bridge (dominant link), both halves are held by the NPIPP1-fragment cluster MCL25, and the PKD1P
loci's family on chr16 IS the NPIP family (MCL1) — "both" cannot exist here. NPIPB5 stays in MCL26 (SMG1P family), as
F1v2 / RG3 leave it.

**12.2 Gorilla fusion simulation** (locus_fix_design arms; 10 fused NPIP copies, 10 controls, 25 T_member copies; ALL
families on each arm's `F1v2.gtf`, 60-99 s; G1 PASS on all 5).

| f | bridges | fused copies in NPIP /10: BASE / F1v2 / COVER | copies /25 (COVER counting bridges) | **both families / eligible** | COVER multi / NULL multi (5 seeds, label `sim_<f>`) | controls /10, partners in NPIP |
|---|---|---|---|---|---|---|
| 0.0 | 1 | 9 / 9 / 9 | 23 | 0 / 0 | 0 / 0 | 10, 0 |
| 0.1 | 4 | 7 / 9 / 9 | 23 | **3 / 3** | 3 / 0-2 | 10, 0 |
| 0.5 | 6 | 1 / 5 / 5 | 19 (20) | **5 / 5** | 5 / 0-3 | 10, 0 |
| 0.9 | 1 | 0 / 0 / 0 | 15 | 0 / 0 | 0 / 0 | 10, 0 |
| 1.0 | 1 | 1 / 1 / 1 | 16 | 0 / 0 | 0 / 0 | 10, 1 |

Where F1v2 makes the fusion a bridge (minority, f ≤ 0.5) COVER puts every eligible fusion in the NPIP family AND its
partner's family (SMG1-like, NSMCE-like, SAGA29-like, CNOT3-like, SNX29), leaving the partition untouched; where the fusion
dominates (f ≥ 0.9) there is no bridge and COVER = F1v2 = BASE. At f = 0.5, 10 of 12 COVER attachments are supported
(own or Liftoff); the 2 unsupported are HOMOLOGY attachments: NPIPA7's bridge → the NPIP family (correct: NPIPA7's own
copy piece sits in MCL5, 8 NPIP-copy loci outside the main NPIP family MCL0, and gorilla NPIP paralogs are below
Liftoff's 0.95, so the judge cannot see it), and NPIPB2~SNX29's bridge → a third family MCL2 (wrong). NPIPA7's bridge
joins three families (MCL0, its own piece's MCL5, and its NSMCE-like partner's MCL25). Dry run of the full driver + `clauses.py` on chr16 and f0.5 (code paths only):
chr16 C1 NJ, C2 NJ, C3 NJ, C4-C6 PASS; f0.5 C2 PASS (5 vs 1-2 multi), C3 reported, C4-C6 PASS.

## Amendments

### Amendment 1 — the freeze (2026-09-29 16:12, written BEFORE any held-out command)

**This file's sha1 before this amendment:** `027d4c0250e57a7c23cca90f2a142264f7d2331a` (26,057 bytes; byte copy
`/mnt/linuxdisk/tmp/rustle_figures/o1_cover_heldout/PREREG_o1_cover_growth_2026-09-29.pre_amendment1.md`). The frozen
text equals this file with the Amendment and Outcome bodies removed. Acceptance: the orchestrating session's task
(design on dev → pre-register → run as pre-registered) is the mandate; no separate user acceptance was sought.

**Frozen instruments** `/mnt/linuxdisk/tmp/rustle_figures/o1_cover_frozen/SHA1SUMS` (sha1 `a00d1774`): `cover.py`
e07d849f, `emu.py` fac9a560 (the RG3 / F1 tests' emulator, unchanged), `score.py` 25242385, `run_score.py` 2412c303,
`lo_score.py` ebfe6f57, `npip_ctrl.py` 527322d6, `clauses.py` 635ff4c2, `ho.py` 4143e38f (driver; it also carries two
DEV pseudo-samples used only for the dry run), and the dev scripts `sim_score.py` 18781a28, `npip_dev.py` d5030d6e,
`power.py` 5047e478, `power_lo.py` 3fdf1006. Binaries `fj_bin_frozen/mcl_families` 91ef2e1c, `family_score` 7723029b.
Python: `python3` (linuxbrew 3.14.4, scipy) for every step; `figures/_liftoff.py` imported read-only from the repo.
**Commands:** `python3 <frozen>/ho.py <step> <sample> [...]` in the §10 order; heavy steps under
`RLOCK_TIMEOUT=590 bash tools/rlock.sh heavy`, the others under `light` (heavy when > 2 GB).

### Amendment 2 — execution only (2026-09-29 17:35, after testis was scored, before any A119b, KB or OR product)

1. **Per-chromosome A119b families also go through `tools/mm2_shard.sh`** (`ho.py` 4143e38f → 2c3b20b8, `SHA1SUMS`
   a00d1774 → 9a23456f). Reason: the testis run showed ≈ 6.6 s of all-vs-all per Mb of spans genome-wide, so chr1 / chr2
   (164 / 165 Mb of spans) could exceed the 590 s call cap without resumability. The wrapper's concatenated PAF equals a
   single minimap2 run byte for byte (its documented, cmp-checked property); no rule, arm, seed, clause, bar, truth or
   substrate changes. The genome-wide path is unchanged.
2. **Lock contention (machine only).** Other sessions hold the heavy lock and the CPUs intermittently (`window.py`,
   `as_table`, `_lrc.py`); the testis ALL run took 5 calls (16:14-17:17). KB shards run at 60-125 s each (88 shards).
   The resumable calls for KB / OR are issued by one serial loop (`for i ...; rlock heavy ho.py fam ...; done`), each
   call under `flock`, so at most one of my heavy jobs runs at a time.
3. **Order (17:50).** After 11 of KB's 88 shards, the KB loop was stopped between calls (PID of the loop shell, cwd
   checked; the running call finished) and the A119b per-chromosome families were run next, then KB and OR resume. The
   §10 order is changed only in time: no substrate's product is scored or opened before its own `score` step, and no
   clause depends on another substrate.
4. **chr13 (A119b) shard size (18:50).** chr13 holds 8,311 loci, 5,681 of them (59 Mb of spans) in its first 20 Mb (the
   acrocentric arm); its first 10-Mb query shard did not finish inside a 454 s budget, twice. The call chain was killed
   by PID (cwd checked) and chr13 is run last with `MM2_SHARD_BP=1000000` (1-Mb query shards; same index, byte-identical
   concatenation). If a single 1-Mb shard still cannot finish in a call, chr13 is dropped from V_A before any A119b
   product is scored, and the drop is reported with the counts it removes.
5. **chr13 dropped from V_A (22:22, before any A119b product was scored or merged).** At 1-Mb shards chr13 has 114
   query shards; shard 1 hit the 453 s deadline and needed a second call, shard 2 was estimated at 281 s: the acrocentric
   arm alone (≈ 60 shards) would take ≈ 5 h per arm. Per item 4, chr13 leaves V_A (`ho.py` 2c3b20b8 → bb2a0fd2: chr13
   removed from the per-chromosome list and added to the dropped contigs; `SHA1SUMS` → 1a589052). What it removes: 12
   of A119b's 348 substrate bridges and 12 of the 1,313 Compara genes. V_A is chr1-12, 14, 15, 17, 19, X, Y, M. The
   chr13 all-vs-all shards made so far are left in the shard cache, unused.
6. **KB and OR not measured (22:30).** A119b was scored at 22:24 and fixed the verdict at REFUTE (C2 and C5 FAIL on V_A,
   §7: any refute trigger on any substrate). No gorilla outcome can change a REFUTE, and finishing them would have taken
   ≈ 6-7 h more under the lock contention (KB: 16 of 88 shards after 4 calls; OR: 101 shards, not started). The KB loop
   was stopped between calls (loop shell killed by PID, cwd checked; the running call finished); OR was never started.
   V_K and V_O are reported as NOT MEASURED; no gorilla held-out COVER product exists.

## Outcome (2026-09-29 22:35)

**Gates.** G0 (sha1s, input GTFs), G1 (emu = binary on the ALL run: testis, and every one of A119b's 19 chromosomes),
G1b (emu = the stored BASE run: testis, 19 A119b chromosomes) and G2 (Python gene space = `family_score` binary: every
human core arm, 9 per substrate) all passed. NULL draws had no shortfall.

| | **human A119b** (per chromosome; chr1-12, 14, 15, 17, 19, X, Y, M) | **human testis** (genome-wide, minus chr16, chr18) |
|---|---|---|
| bridge loci on the substrate; attached; multi-family | 336; 33; **9** | 12; 2; 0 |
| attachments: SELF / HOMOLOGY (supported) | 36 (36) / 9 (7) | 2 (2) / 0 |
| Compara bipartite sens / prec / F, core: BASE → **COVER core (= F1v2)** | .2915 / .9365 / .4446 → **.2966 / .9375 / .4507** | .1532 / .9500 / .2639 → **.1540 / .9502 / .2651** |
| Compara pairs TP / pred: BASE → COVER core → **COVER cover-aware** | 551 / 643 → 554 / 646 → **558 / 659** (+4 true, +9 false) | 210 / 227 → 212 / 229 → 212 / 229 |
| cover-aware F; pairwise precision (core → cover-aware) | .4507 → **.4504**; .8576 → **.8467** | .2651 → .2651; .9258 → .9258 |
| ALL (bridges as loci) F; NULL0-4 core F | .4504; .4443-.4494 | .2651; .2603-.2639 |
| Liftoff pair recall: BASE / COVER core / COVER | .0745 / .0745 / .0745 (322 pairs) | .303 / .333 / .333 (33 pairs) |
| multi-family loci: COVER vs NULL0-4 (of 336 / 12) | **9 vs 5, 5, 9, 12, 9** | 0 vs 0 × 5 |
| judged attachments of multi-family loci supported: COVER vs NULL | 19 / 21 vs 2/15, 2/12, 5/19, 5/30, 5/26 | – |
| NPIP copies in the NPIP family: BASE / COVER (C6) | NA (chr16 not computed) | 5 / 5 of 26 |
| **C1** | NJ (0 eligible, as known) | NJ |
| **C2** | **FAIL** (a) 9 is not above 9 and 12; (b) 19 > 2 holds | NJ |
| **C3** | PASS (7 of 9 HOMOLOGY attachments supported) | NJ |
| **C4** | PASS (every Compara metric and Liftoff recall ≥ BASE) | PASS |
| **C5** | **FAIL** (cover-aware F .4504 < .4507; pair precision .8467 < .8576) | PASS |
| **C6** | NA | PASS |

gorilla KB3781 and OR6737: **NOT MEASURED** (Amendment 2, item 6).

**VERDICT (§7): REFUTE** (`clauses.py`: human_A119b C2 and C5). COVER stays a prototype; nothing is proposed for a
default or a Rust port.

**Read this before quoting it.**
- **C5 is one attachment.** All 9 false pairs come from the LINC01859|ZNF728 bridge attaching to `chr19_MCL1`, a
  12-gene chr19 cluster, mostly KZFPs, that already mixes 7 Compara families (ZNF100 / ZNF430 CF217, ZNF66 / ZNF675 CF223, ZNF85 / ZNF724
  CF42, ...). The attachment is truth-consistent (ZNF728's Compara family CF224 is MCL1's plurality: ZNF254, ZNF708,
  ZNF714) and adds 3 of the 4 true pairs; the pairwise metric charges it for the family's impurity. The other +1 true pair
  is NBPF1 joining an NBPF sub-family. The margin is .0003 F; the clause had zero tolerance and it fired.
- **C2: multi-membership is correct but not specific.** COVER's 9 multi-family bridges are mostly right (19 of 21
  judged attachments supported), random degree-matched loci's are mostly wrong (2-5 of 12-30), but random loci become
  multi-family just as often (5-12 of 336). The pre-registered clause asked for frequency, and frequency did not
  separate them. Post hoc: where MCL splits one large SD / KZFP / NBPF / ANKRD20A family into sub-families, ANY node
  touching two sub-families is attached to both.
- **What the 9 COVER multi-family bridges are (post hoc, not judged):** two KZFP tandem readthroughs placed in both
  parents' families (ZNF91|ZNF724 → MCL1 + MCL10; ZNF808|ZNF701 → MCL20 + MCL5, each family holding one parent's own
  locus); LINC01859|ZNF728 (MCL11 + MCL1); an NBPF1 fusion in 4 NBPF sub-families; single-gene loci whose pieces and
  paralogs sit in 2-3 MCL sub-families (CNTNAP3B, ZNG1C, LOC128966611, ANKRD20A2P region, LOC124904395). None has two
  parents in two Compara families, so C1 could not judge any of them.
- **The partition part is a gain:** COVER's core (= F1v2) beats BASE on every Compara metric on both human substrates
  and on Liftoff recall on testis (C4). This is the first genome-wide family-level measurement of F1v2 (A119b per
  chromosome); it is a property of F1v2, not of the cover.
- **SE (reported only) over-covers massively:** 1,263 of 29,847 single-exon loci attach, 579 to ≥ 2 families, 74 of
  1,928 judged attachments supported; core Compara F .4507 → .4275. Soto's periphery does not transfer as "one exon".

**Predictions (§8).** P1 gates: hit. P2 COVER multi on A119b (0.4): hit (9); testis none: hit; C2 judged (0.5): hit; C2
PASS (0.6): **missed**. P3 HOMOLOGY on A119b (0.5): hit; C3 PASS (0.6): hit. P4 C4 PASS A119b / testis: hit; COVER core
F ≥ BASE on A119b (0.6): hit. P5 C5 PASS (0.9): **missed on A119b**, hit on testis; C6: hit. P6 attached 10-40% on A119b
(0.6): **missed** (9.8%); SELF ≥ HOMOLOGY (0.8): hit. P7 verdict REFUTE (prior 0.40): hit. KB / OR predictions not
scored.

**Falsifiers (§9).** Z1 holds (SELF 36 ≥ HOMOLOGY 9; 2 ≥ 0). Z2 holds (9 of 336 < ⅒). Z3 holds (COVER-core F above
every NULL core: .4507 vs ≤ .4494, .2651 vs ≤ .2639; Liftoff .0745 vs .0683-.0745). Z4: random peripheralization makes
as many multi-family loci as the bridges (5-12 vs 9), mostly unsupported.

**What this does and does not show.**
- On a held-out human library, making F1v2's bridges peripheral and attaching them Soto-style puts 33 of 336 in a
  family and 9 in several, mostly through their own pieces (SELF 36 of 45) and mostly correctly, and never merges or
  changes a family.
- It is not specific (random loci of the same degree become multi-family as often) and it can lower a pairwise score by
  attaching a correct member to an impure family. Both pre-registered safety clauses fired, so the rule is refuted as
  specified.
- The gain the rule was built for — a fusion in both parent families — was shown only on the simulation (dev: 3/3 and
  5/5 eligible fusions at f = 0.1 / 0.5) and could not be judged on any held-out substrate: no F1v2 bridge joins two
  multi-copy (Compara or Liftoff) families. The motivating fusions (PKD1P6-NPIPP1 and the other PKD1P / PDXDC2P
  readthroughs) are dominant and never become bridges.
- Not measured: gorilla KB3781 / OR6737; chr13 of A119b.

**Kept products.** Tables `/mnt/linuxdisk/tmp/rustle_figures/o1_cover_heldout/tables/` (`verdict.json` e0c94404,
`human_A119b.{gene_space,att,liftoff}.json` b0d637c2 / 8794e6b8 / e1505ee6, `human_testis.{gene_space,att,liftoff,npip}
.json` e07a7a5b / 8000ff98 / 625ad793 / aae7e7e5, dev `sim_score.json`, `npip_dev.json`). Scratch
`/mnt/linuxdisk/tmp/rustle_figures_dev/o1_cover/` (8.1 GB: dev, held-out families runs, COVER / NULL / SE products);
shard cache `/mnt/linuxdisk/tmp/mm2_shard_cache/` (21 GB, includes the unfinished KB and chr13 shards). Report
`scratchpad/figs/o1_cover.md`.
