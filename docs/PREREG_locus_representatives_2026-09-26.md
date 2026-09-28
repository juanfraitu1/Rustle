# Pre-registration: locus representatives, one arm: regroup after polish (RG, rule RG3)

> **STATUS: PARKED (user, 2026-09-27) — no held-out run of any arm.** All four arms (A, B, C, RG) are answered
> descriptively on the development contigs only; see register rows 1121-1127. RG was parked after four hostile
> reviews for cost (~18-24 h full families, ~55-60 GB disk) rather than for a failed clause. The design below is kept
> as the starting point if RG is resumed; it would need a fresh acceptance and Amendment 1.

**Written 2026-09-26 as a four-arm file. Narrowed to arms C and RG on 2026-09-27, then to arm RG alone on 2026-09-27
(this version, KEY=revise4), before any held-out number of any arm exists.**
- This version replaces the C+RG version of this same untracked file: sha1 b9d72f0c6e091d5be936ce30fb8cd435a37855d1,
  archived read-only at
  `/mnt/linuxdisk/tmp/rustle_figures_dev/rep4_revise4/PREREG_locus_representatives_2026-09-26.narrowed_CRG.md`.
- That version replaced the four-arm version: sha1 037cdf27249ab637178bd5f408cbb584f5cc866a,
  `rep3_revise3/PREREG_locus_representatives_2026-09-26.fourarm.md`.
- Neither earlier version bound anything.

This file becomes binding when Amendment 1 records two things:
1. the user's explicit acceptance of this file, of the third reuse (§7.3) and of the lock time (§13.1: about 7-14 h
   expected, about 24 h at most);
2. this file's sha1.

Until then it binds nothing, and nothing held-out runs except the BASE-only prerequisites of §13 step 3, which produce no
arm number.

**User request (2026-09-26):** "see if we can make a better representative or set of representatives, and avoid
abstentions for O2, and a rule for copies with no reads". It answers the advisor's criticism: "our locus
representative might lose a lot of info; Liftoff projects all exons, isoforms and transcripts, whereas we just project
loci."

**User decisions, as relayed by the orchestrating session; Amendment 1 records the user's own confirmation:**
1. **NARROW (2026-09-26, night).** Arms **A** (rule R7, isoform evidence at the family edge) and **B** (the O2 tie-set
   answer TS0) are **dropped from held-out**. §3 reports their development evidence descriptively.
2. **C DESCRIPTIVE ONLY (2026-09-27).** Arm **C** (the semi-guided mode's copies with no reads) is reported as development
   evidence in §3.4, **with no held-out verdict; to be revisited later**. §3.4 also lists what a later test of C must
   fix. No held-out number of A, B or C is produced under this file. A later held-out test of any of them on these six
   samples is a fourth use (§7.3).
3. **Only RG goes to held-out,** with rule **RG3**: split a `gene_id` only where its surviving transcripts share no
   exonic base on the same strand (§2.1; `rep3_RG.md`). RG is judged on its own clauses (§8).
4. **RG's families cost is cut by an incremental re-alignment instrument, validated on dev before any held-out run**
   (§2.2, gate IRG7; `rep4_incr.md`). If the instrument does not pass its gates, the full families run is the method,
   and §13.1 states that cost.
5. **Substrates:** the six held-out samples, reused a **third** time (§7.3). Human chr16 and chr20 and gorilla
   NC_073244.2 are removed from scoring in every sample (§7.2).
6. **Prerequisites.**
   - **Queue 1**, `/mnt/linuxdisk/tmp/rustle_figures/rep_prereq_queue.sh` (log `rep_prereq_queue.log`), has been running
     since before acceptance.
   - **Queue 2**, `/mnt/linuxdisk/tmp/rustle_figures/rep_prereq2_queue.sh` (sha1 638f48fb), was **written on 2026-09-27
     and has not been started.** It builds what this file adds:
     - the genome-wide protein referee for four species;
     - the chimpanzee and orangutan genes tables, after gate IRG8.
   - Neither queue produces an arm number (§0, §13). Queue 2 starts after queue 1 finishes, on the user's go, before or
     after Amendment 1.

RG is not a default flip. A flip or a ship is the user's call, whatever the outcome.

## 0. What was seen before this file, and the procedural record

- **Seen: development contigs only (design evidence).** Human A119b chr16 and chr20, and gorilla OR6737 NC_073244.2.
  - Scout reports (`scratchpad/figs/`): `rep_edges.md` (A), `rep_o2.md` (B), `rep_noreads.md` (C), `rt4_ghost.md` (RG).
  - Hostile reviews:
    - first: `rep_critique.md`;
    - second: `rep2_review.md` (sha1 7cdb8187);
    - third: `rep3_review.md` (sha1 d86bfb7d; NOT READY: H1 costs and queueing, H2 C's vote floor).
    - §14 lists every item of each review and its fix, or the arm it went with.
  - Dev rebuilds and redesigns, with scratch in `/mnt/linuxdisk/tmp/rustle_figures_dev/rep*/`:
    - `rep2_B.md`, `rep2_C.md`, `rep2_RG.md`;
    - `rep3_C.md` (18f4c7ea), `rep3_RG.md` (d950ad1e);
    - `rep4_incr.md` (a8784fc2): the incremental families instrument.
  - Isoforms per locus in the three dev BASE GTFs (`rt3_impl/runs/new/<sub>/<sub>.seed.unset.gtf`): chr16 3.38
    (9,473 / 2,802), chr20 2.99 (5,631 / 1,882), NC_073244.2 3.42 (3,911 / 1,142).
- **What this revision's author (KEY=revise4) read.** The reports above and nothing held-out:
  - the scripts `rep_prereq_queue.sh` and `queue_lib.sh`;
  - the progress lines of `rep_prereq_queue.log`, with each line's result tail stripped (`sed 's/ :: .*$//'`), which
    leaves task names, call counts and seconds;
  - the process table;
  - the four referee pilot logs in the dev truth directory: contig and protein counts, shard timing, no family;
  - the source of `bench/truth.py` (protein-homology), `figures/_o1.py` (`annotation_cache`, `PH_DIR_DEFAULT`),
    `fam_call.sh`, `incr_mm2.py`'s docstring and `rg3_null.py`'s protocol;
  - the contig order of the chimpanzee and orangutan RefSeq GFFs (26 contigs each, one block per contig), which are
    annotation inputs.
- **Seen: held-out BASE rows, already published by the readthrough study** (review2 M1).
  - The six files `rt_arms/tables_v3/readthrough_eval.<sample>.annotated*.{tsv,json}`. Scope:
    - A119b rows exclude chr16 and chr20 (`annotated_minus_chr16_chr20`);
    - OR6737 rows exclude NC_073244.2 (`annotated_minus_NC_073244.2`);
    - testis, KB3781, PTR and PPY rows are whole-genome, so testis's include chr16 and chr20, and KB3781's include
      NC_073244.2.
  - For all six samples: `a.fused`, `a.absorbed_genes`, `b.tes_recovered`, `b.tss_recovered`, `d.found_annotated` and
    `d.reciprocal_one_to_one`, for the BASE, R and R3 arms. The v3 outcome and register rows r1117-r1120 quote them.
  - For testis and PTR, also `e.largest_family`, `e.liftoff_pair_recall` and the Compara and Soto rows:
    - human_testis BASE: Compara pairwise TP 222 of 3,053 truth pairs, 239 predicted; Soto 73 / 1,826, 91 predicted;
      338 families, largest 47;
    - chimp_PTR BASE: Liftoff pair recall 17 / 145; 379 families, largest 40.
  - **These are the baselines of RG's T, L, H and P2 clauses and of FU's neighbourhood.** RG predictions 1, 3, 4, 5 and
    6 (§10) were therefore set with that knowledge and are marked post-exposure. Prediction 2 (S) is not: no split class
    of any held-out sample was seen.
- **Seen: held-out, descriptive of BASE GTFs, not an arm.**
  - The main session's `rep_loss.py`:
    - gorilla_KB3781 has 4.2 isoforms per locus, and its representative carries a median 98.7% of its locus's exon bp;
      pooled, 20% of exon bp and 34% of junctions lie outside it;
    - human_testis has 1.9 isoforms per locus, with 11% and 22% outside.
    - It informed arm A's dropped predictions and RG's FU failure mode (§8.3).
  - The RG author (`rep2_RG` §1) read only the dev-contig rows of the held-out BASE GTFs of A119b (chr16, chr20) and
    OR6737 (NC_073244.2). They equal rt3's dev GTFs except `TPM`. The author also recorded the six BASE GTF sha1s.
- **Seen: earlier figure work.** Fig. 7 per-contig tables for A119b chr2, chr6, chr8 and chr10 and OR6737 NC_073234.2
  (de novo families of BASE type against Compara, Soto and the protein referee). These contigs lie inside V1 and V3, so
  V1 and V3 are also reported without them (§8.2).
- **Seen: queue 1 costs (timing only).** BASE families wall time:
  - gorilla_OR6737: 13 calls, about 104 min;
  - gorilla_KB3781: 9 calls, about 72 min;
  - human_testis: about 21 min (from earlier);
  - chimp_PTR: about 70 min (from earlier);
  - orangutan_PPY: running at the time of writing;
  - human_A119b: queued.
- **Does not exist, or not seen.**
  - No held-out number of RG3, NULL_RG3 or the old `rg.py` rule; no held-out number of any C rule.
  - **Queue 1's products:** Liftoff self-lifts for gorilla, orangutan and human; read-support tables for five samples;
    BASE families for human_A119b, gorilla_OR6737, gorilla_KB3781 and orangutan_PPY. **No author of this file or of
    `rep3_*` / `rep4_*` read them.**
  - **The genome-wide protein referee** (review3 H1, point 4) is only a partial pilot, with no families:
    `rustle_figures_dev/truth/gw/<species>/ph.*`, 1 of 21 to 23 blastp shards per species, over 24-26 contigs. Queue 2
    completes it.
  - No genes tables for chimpanzee or orangutan: `families_gw/species/` holds human and gorilla only. Queue 2 builds
    them after gate IRG8.
  - Not yet written: the budget-sharded passes of the incremental instrument and `fam_call_incr.sh` (§6, IRG7).
- **Procedural record (plain).**
  1. **Every rule was selected on the three development contigs.** Their dev numbers are the optimistic edge.
     - RG's adjacency was changed after review2 X2. Four adjacencies were measured on dev:
       - exact junction key: the old rule;
       - exon overlap;
       - exon-or-junction: identical to exon on 6/6 dev GTFs;
       - donor-or-acceptor: fails S on human dev.
     - Exon overlap was named before measuring, and kept. The piece names `.rg<k>` replaced `.<k>` (review2 L1).
     - The incremental instrument was designed on the same six dev families runs it is gated on (`rep4_incr` §3: two
       earlier predictors were rejected there). Gate IRG7-H (§6) is its held-out check.
     - C's selection history is in §3.4.
  2. **Machine-rule record.**
     - `rep3_RG` ran its six families runs under the lock.
     - `rep4_incr` ran its full, chaining-only, build-fidelity and incremental runs under the lock, at ≤ 2.5 GB. Three
       of its lock waits went through the tool's background mode, and the work inside the lock ran in the foreground
       (disclosed there).
     - This revision ran no heavy step. It ran one foreground read of two GFFs (contig order), one `paths` dry run of
       the genes-table helper (a registry lookup), and file copies (the frozen referee code, §6).
  3. **After Amendment 1, dev results cannot change a rule, a NULL, a clause, a floor or a tolerance** (review2 M6).
     - The dev tables that Amendment 1 records (§6) are reproductions of the numbers in `rep3_RG` and `rep4_incr` by the
       frozen instruments. The user accepts with them in view, and a dev-gate failure is an instrument bug, fixed in
       the instrument, never in the rule.
     - **Exception (review3 M6): IRG4's complete-universe dev rows are new numbers.** Nobody has seen them. They are
       recorded as found and cannot change anything.
     - The held-out stage runs whatever dev shows, unless the user stops the work.

## 1. The advisor's point, and where the representative actually enters

**The alignment is already Liftoff-like.**
- The families stage (`mcl_families --from-gtf`) aligns whole locus bodies all-vs-all. Each body spans every transcript
  of the `gene_id`, introns included, and Liftoff aligns gene bodies the same way.
- Rewriting the GTF so that a different transcript is the representative leaves `loci.paf` byte-identical on all three
  dev contigs (`rep_edges` §0).
- So "we only project loci" is false for the alignment. It is true for the four places below.

| where the locus model enters | what it sets | arm |
|---|---|---|
| the family edge filter | the exon model of the admission clauses, and the edge weight identity × cov_longer | A: **answered descriptively (§3.2)** |
| O2's copy table (`--emit-units`, one row per member locus) | the unit span that labels a read's placements | B: **answered descriptively (§3.3)** |
| nothing today | copies with no reads are never nodes, so nothing searches for them | C: **answered descriptively (§3.4)** |
| node identity: which transcripts share one `gene_id` | a polish can drop the only transcript bridging two pieces (a "ghost"), and a non-unique tid can join two components (a "collision"). Either way one representative can stand for two genes | **RG (held-out)** |

On dev, what one representative "loses" for families is mostly **other genes held in the same node**. Regrouping moves
Compara pairwise sensitivity 0.342 → 0.451 on chr16, about 20× any representative effect (`rep_edges` §3). RG tests that
directly.

## 2. The binding rule and its families method

### 2.1 Rule RG3: regroup after polish (same-strand exon-overlap components)

Verbatim from `rep3_RG` §1, realised by `rg3.py` ec17e5408dedc82e0d7201e62f77bb181a45ff95:

```
Input: a frozen BASE GTF (rt_arms/<s>/<s>.BASE.gtf = runs/<s>/<s>.gtf). Transcript = a `transcript` line; index i =
its order among transcript lines; strand = its column 7; exons = its `exon` lines (1-based closed, sorted); reads = the
`reads` attribute (absent -> 0); span = end - start + 1 of the transcript line.
Pieces: connected components of ONE gene_id's transcripts, two transcripts adjacent iff same chrom, same strand and
>= 1 shared exonic base (max(a1,a2) <= min(b1,b2) for some exon pair).
Piece representative = max (reads, span, -i).
Naming: a gene_id that is one piece is untouched. A gene_id with m >= 2 pieces: the piece whose representative is
max (reads, span, -i) keeps the gene_id; the others become "<gene_id>.rg<k>", k = 2..m in representative-index order.
Output: the input with `gene_id "<old>"` (first occurrence) replaced on every line of a renamed transcript; no line
added, removed or reordered.
Asserted on every run (exit 3): (1) identical to the input once gene_id is removed; (2) one piece per output gene_id;
(3) no new name equals an input gene_id OR an input transcript_id; (4) split-only; (5) no two output gene_ids of one
input gene_id share an exonic base on one strand.
```

- **Disclosure: a dev selection made after review2 X2.**
  - The four-arm version registered `rg.py` 92892fb3 (exact junction key). It fails S on 3/3 dev contigs (§4.1,
    r1125).
  - "Exon overlap" was named as the candidate before it was measured.
  - "Or shared junction" adds nothing: a shared junction implies a shared donor base on one strand.
  - Donor-or-acceptor fails S on human dev.
  - RG3's splits are a subset of the old rule's on 6/6 dev GTFs (0 new split `gene_id`s).
- **Properties.**
  - RG3 is threshold-free ("≥ 1 shared base" means "overlaps") and split-only.
  - It changes `gene_id` only, so intron chains and every `c.*` row are identical by construction.
  - RG3 runs on the BASE GTF only. R3 (EFFECTIVE, r1119) is not the default, so RG3 on R3 is dev-descriptive only (§4.1).
- **It is not a re-proposal of a closed result.**
  - It is not r846's node cut: it cuts no transcript, and every piece is an assembled exon-overlap component.
  - r1013 and r1017 closed trims and splits of fused loci at coordinates, and r1017's +0.021 is a ceiling, not a rule.
  - It merges nothing. The 1,818 / 588 / 143 BASE pairs of distinct spliced `gene_id`s that already share an exonic base
    on one strand (chr16 / chr20 / NC_073244.2) are the assembler's own grouping, out of scope.

### 2.2 The families stage for RG3 and NULL_RG3 (user decision 4; gates IRG7 and IRG7-H, §6)

The families stage is the shipped one on the arm GTF: `tools/rustle_pipeline.sh families` (d1816c35) with the frozen
`mcl_families` 6cf8183f and `--threads 4`. Only the all-vs-all alignment program behind `RUSTLE_MINIMAP2` differs
between the two methods.

```
Method F_INCR (used iff IRG7 passes on dev and IRG7-H passes on V2):
  RUSTLE_MINIMAP2 = the incremental instrument SHA1_INCR (incr_mm2.py fed8d80e plus the budget-sharded passes of
  IRG7; incr_core.py 5fbfa317; bin/minimap2-tmask 450a7316 = minimap2 v2.30 tag 79c9cc18 + mm2tmask.patch 4ed52d4f;
  stock minimap2 2.30-r1287, 0c244670), called by fam_call_incr.sh <s> RG|NULLRG (SHA1_FAM_INCR), with
  INCR_BASE_PREFIX = rt_arms/<s>/<s>.BASE.fam and INCR_BASE_MMI = BASE's mm2_shard-cached idx.mmi (or
  INCR_BASE_MIDOCC from BASE's shard log). Kept = arm records whose name, sequence and multiplicity equal BASE's.
  1. Index the arm FASTA (the full run's index). Probe mid_occ in it and in BASE's index. If they DIFFER: the full
     all-vs-all on the arm index (exact).
  2. Otherwise: a chaining-only pass (the families flags minus -c) over the kept records defines A_chain = the kept
     queries with a chain to a new record; the new and duplicate-name records get full blocks; A_chain is aligned to
     the new records only (target-mask patch: chains to unlisted targets are dropped after chaining, before
     alignment); every other kept pair keeps BASE's lines, minus the lines to split parents.
Method F_FULL (otherwise): fam_call.sh <s> RG|NULLRG (799f3e71): RUSTLE_MINIMAP2 = tools/mm2_shard.sh (04c2c639),
  the full all-vs-all (-x asm20 -c -X -N 50 -p 0.1 --secondary=yes -t 4).
```

- **What F_INCR is not.** Its PAF is not byte-identical to F_FULL's, and cannot be made so with any useful speedup
  (`rep4_incr` §2). A kept query's block depends on the whole target set in three ways:
  - seed use: occurrence counts, `mid_occ` and the high-occurrence rescue ranking;
  - anchor tie order: an unstable radix sort;
  - line order: a tie hash that includes the target's rank in the index.

  The oracle recompute set for a byte-identical PAF carries 88-98% of the alignment work.
- **Its acceptance object is `clusters.tsv`.** That is the only family product any RG scorer reads (besides the GTF and
  the truths).
- **Its semantics, stated.**
  - Under F_INCR a kept pair's lines are BASE's. So an RG3 − BASE graph difference comes only from changed records
    (exact lines) and from removed parents.
  - A NULL_RG3 run that falls back is exact.
  - Each run's `report.json` (mode, `mid_occ` of the arm and of BASE, N, A_chain, seconds) is reported, and so is each
    substrate's pairing: both incremental, one fallback, or both fallback.
- **Expected held-out regime (untested).**
  - The asm presets clamp `mid_occ` to [50, 500]. The dev loci sets of 25-55 Mb have raw values of 93-141.
  - Genome-wide, `mid_occ` should be clamped at 500 in BASE and in both arms. The fallback would then not fire, and
    NULL_RG3 would take the incremental path.
  - Dev showed that path on NULL_RG3 only in forced mode, and only with the superseded v2 wrapper (§4.2). That is why
    IRG7 requires forced mode on dev with the final wrapper, and IRG7-H cross-checks on held-out.

### 2.3 Where every binding number comes from (none is tuned here)

| constant | value | source |
|---|---|---|
| O1 edge constants | identity 0.70, block 300, cov_longer 0.30, exonic bp 1, shared fraction 0.60, inflation 2.8, prune 1e-9 | shipped `mcl_families` defaults and driver flags, used unchanged by RG's families stage |
| RG3 adjacency | same chrom, same strand, ≥ 1 shared exonic base | overlap; not a threshold |
| RG tolerances (P1, P2, H, T, L) | **0** (exact integer inequalities) | feasible because RG3 is split-only; every dev contig passes them exactly |
| "more than" / majority (S) | strict sign count | no constant; a tie fails |
| power floors | 30 truth families, pairs or **judged splits** (S); **fused gene PAIRS ≥ 50** (FU); universe ≥ 1,000 genes; ≥ 1,000 reference loci | R prereg (`PREREG_readthrough_ends_representatives_2026-09-25.md`) §4 floors (e), (a), (b), (d). **Two unit changes are stated** (review2 L2): floor (a) was 50 fused LOCI and is applied here to fused gene PAIRS; floor (e) of 30 is applied to judged splits (SEP + FRAG) |
| verdict-substrate rule | ≥ 1 human and ≥ 1 ape | R prereg §4 |
| seeds | NULL_RG3: string `"20260926:<sample>:<contig>"` | fixed in the draft; inputs sorted before any draw |
| NULL_RG3 depth bin | floor(log2 reads) of the piece representative | `rg_null.py` 08b5c317, carried verbatim into `rg3_null.py` (a matching granularity, not a rule threshold) |
| Liftoff pair truth (P2) | Fig. 6s rule: `sequence_ID` ≥ 0.95, exon unions ≥ 200 bp, ≥ 2 reads, ≥ 50% covering | `PREREG_genome_wide_families_2026-09-25.md` Amendment 1 item 3 |
| ends windows (T) | TES ±25 bp, TSS ±250 bp | R prereg §3 (b), from §6w3's dispersion |
| protein referee | longest CDS per gene, pseudogenes and V(D)J excluded, blastp e ≤ 1e-5, coverage ≥ 0.30 of the longer protein, MCL I = 2.8 | the §6ko rule in `bench/truth.py` (2c59a2d8), not re-tuned |
| F_INCR's fallback | full all-vs-all iff `mid_occ` (arm) ≠ `mid_occ` (BASE) | an equality, not a threshold (`rep4_incr` §4) |
| F_INCR's "kept" | name, sequence and multiplicity equal | an equality |
| call budgets | families 480 s per shard call; referee 540 s; `timeout 600` | machine rules; they change no output (IRG7 (c) checks this for the instrument) |

The C constants of the previous version (projection search, chain geometry, RepeatMasker, truth match, reference floor,
Liftoff for C) went with the arm (§3.4).

## 3. Answered descriptively on dev (no held-out verdict): arms A, B and C

Everything in this section is **dev design evidence**, from human A119b chr16/chr20 and gorilla OR6737 NC_073244.2.
- The rules were selected on these contigs.
- **No held-out verdict exists or will be produced under this file.**
- Register rows: r1121 (A), r1122 (B), r1123 and r1124 (C).

### 3.1 The answer to the advisor's criticism

The advisor's criticism was that one locus representative "might lose a lot of info", because Liftoff projects all
exons, isoforms and transcripts whereas we project loci. On the development contigs the answer has five parts.
1. **For families, the alignment is already Liftoff's.**
   - The families stage aligns every locus's whole body (the span of all its isoforms, introns included) all-vs-all.
   - Choosing a different representative leaves the alignment byte-identical.
   - The representative enters only the exon filter of an edge, its weight, the within-family fold and O2's unit.
2. **Isoforms add little at the family edge.**
   - Letting every isoform vote adds 13% edges on chr16 but at most +0.005 Compara pair sensitivity. It costs Compara
     pair precision 1.000 → 0.753, through a 48-locus NPIP + SMG1P hub (the r676 shape).
   - The representative-anchored variant R7 adds 4.9-11.2% edges. It leaves every judged pairwise truth unchanged
     (chr16, chr20) or better (gorilla: +15 referee pairs, 0 false), at the price of one referee gene on chr16.
3. **What one representative really hides is other genes in the same node.** Regrouping the node moves chr16 Compara
   pair sensitivity .342 → .451, about 20× the largest isoform effect. That is why RG is carried.
4. **For O2, the abstentions are not the representative's.**
   - 96-100% of AS-tied molecules have an NM-identical twin that no unit model can break.
   - The representative causes 285 of 14,706 real chr16 abstentions (1.9%).
   - Where it does lose information is the unit's *label*. In fused NPIP-region loci a quarter of simulated tie sets
     miss the source copy. All-isoform units fix 156 of those out of sample, but name 42 foreign copies.
5. **For copies with no reads, projecting full locus models the way Liftoff does finds more of them than literal
   Liftoff** on the gorilla dev contig (copy level b = 8, c = 0). But no acceptance rule survived dev, and C is reported
   only (§3.4).

In one line: the representative loses little that isoform evidence recovers at the family level; it mislabels O2 units
in fused loci; and the dominant loss is node identity.

### 3.2 Arm A: isoform evidence at the family edge (`rep_edges.md`; `emu.py` fac9a560 byte-exact to the frozen binary)

**R7** (representative-anchored isoform evidence): a pair is an edge iff both of these hold:
- (i) some single retained PAF record carries ≥ 1 exon base of both representatives;
- (ii) some pair of read-supported isoforms in the representatives' exon-overlap components passes O1's quantitative
  clauses.

The weight is the max over passing pairs. Nodes, fold, MCL and the O2 copy table are unchanged. Full text: the archived
four-arm version, §2.A.

| contig | arm | edges (+ / −) | families | loci in families | largest | Compara pairs TP / pred / truth | protein referee TP / pred | gained edges labelled TP / FP |
|---|---|---|---|---|---|---|---|---|
| chr16 | BASE | 5,822 | 111 | 353 | 27 | 66 / 66 / 193 | 100 / 101 | – |
| chr16 | R7 | 6,109 (+287 / −0) | 110 | 362 | 25 | 66 / 66 / 193 | 100 / 101 | Compara 12 / 0; referee 20 / 1; Soto 22 / 3 |
| chr16 | R1 (every isoform votes; the ceiling) | 6,602 (+780 / −0) | 138 | 454 | **48** | pS .342 → .347; pP 1.000 → **.753** | pP .990 → .707 | hidden pieces: 81 edges, **0 TP** |
| chr20 | BASE / R7 | 684 / 729 (+45) | 26 / 31 | 97 / 108 | 14 / 14 | (4 Compara families) | 3 / 4 both | – |
| NC_073244.2 | BASE / R7 | 430 / 478 (+48) | 29 / 32 | 93 / 105 | 21 / 21 | – | 275 / 275 → 290 / 290 | referee 38 / 0 |

- Relative edge gain: +4.9% (chr16), +6.6% (chr20), +11.2% (NC_073244.2). chr20 has the fewest isoforms per locus
  (2.99) but a larger gain than chr16 (3.38), so "the gain comes from isoforms" is not supported on dev.
- R7's one cost: −1 referee gene on chr16, bipartite F 0.225 → 0.220; pairwise unchanged. Liftoff within-contig pairs
  go 15 → 16 of 19.
- The relabel null is zero-width on dev, so every family change is the variant's.
- **Never measured:** NULL_A and NULL_Ar (matched and rank-matched relaxations), and the complete-universe referee.
  Without them, R7's gained edges cannot be attributed to isoform evidence rather than to relaxing the representative's
  own clause.
- The unanchored variants R1, R2, R2m, R4 and R6 rebuild the 46-48-locus NPIP + SMG1P hub (r676). R3 is a trade; R3c
  is worse than BASE.

### 3.3 Arm B: the O2 answer on AS-tied molecules (`rep_o2.md`, `rep2_B.md`, review2 X1/X3/H1; `copy_assign` f480847a)

| chr16 simulation, 2,127 AS-tied molecules (all NM-identical) | answered | single | wrong single (placement level) | set misses origin | Σ\|A\| |
|---|---|---|---|---|---|
| PRIMARY | 2,127 | 2,127 | 1,056 | 0 | 2,127 |
| PF (shipped default per-family table) | 993 | 560 | 548 | 147 | 1,699 |
| UC (union certificate, opt-in since r1103) | 0 | 0 | 0 | 0 | 0 |
| TS0 (the aligner's tie set as copies + `outside`) | 2,127 | 47 | 0 | 1 | 7,548 |
| TS1 (all-isoform units) | 2,127 | 25 | 0 | 1 | 7,949 |
| ISOPRIOR (expression prior) | 439 | 439 | 229 | 0 | 439 |

- **No unit model removes an abstention.**
  - Simulated ties are 2,127/2,127 (chr16) and 82/82 (gorilla) NM-identical; real chr16 ties 16,083 of 16,717 (96.2%).
  - The representative causes 285 of 14,706 real chr16 abstentions (1.9%, 258 of them twins).
  - At most 551 (3.7%) are resolvable by read evidence at all (`rep_o2` §4).
- **At the copy level TS0 is wrong where the representative's span is the unit** (review2 X1).
  - All 47 TS0 single answers are `{outside}`: two placements at chr16:21 Mb, neither overlapping a T0 unit.
  - All 47 are wrong by the simulator's source-copy id, and 22 by the locus-extent rule.
  - 547 of the 2,126 placement-valid tie sets (25.7%) do not contain the source copy. On gorilla, 82/82 name the source.
- **All-isoform units (T1)** name the source 156 more times. Leave-one-isoform-out gives the same 156 with 0 loss, so the
  gain is not self-inclusion. The cost is 42 foreign names (41 in MCL3).
- **The union certificate contradicts the read's own alignment** (review2 X3).
  - It certifies 0 of 1,599 simulated tied molecules in scope.
  - On real chr16, 244 of its 298 certified tied molecules go to a copy whose extent holds none of the read's best-AS
    placements (241 in MCL1, 196 to the fused NPIP unit MCL1:1).
  - It stays opt-in. This is a recorded warning, not a verdict.
- F-exact fails: 1 of 2,127 tie sets misses its origin (origin AS 3,371, twins 3,373 at the same NM).
- Scope never covered (review2 M7): molecules were simulated only from assembled copies, so reads from copies the
  catalog lacks never supplied an origin.

### 3.4 Arm C: copies with no reads, the semi-guided mode (`rep_noreads.md`, `rep2_C.md`, `rep3_C.md`; user decision 2)

**The semi-guided mode, as defined on dev.** It is the guided mode's "find new candidate loci" step, with the de novo
read-supported family models as the seed annotation (memory `project_two_modes_scope`). Every candidate descends from an
expressed member, so it is never genome-only discovery.
- **Seeds.** The members of one O1 family F: `mcl_families --from-gtf <BASE GTF of the contig> --min-exonic-bp 1
  --min-shared-exon-frac 0.60 --emit-units` (rt3_bin_frozen 0c4639e2).
  - A member's locus model is its full locus: every assembled isoform of its `gene_id`.
  - `mcl_families` emits no singleton families. So only genes with **≥ 2 expressed loci on the contig** are searched: on
    NC_073244.2, 29 families (sizes 2-21) hold 93 of the 1,142 BASE loci (review3 M3).
- **Projection, as Liftoff does it.** Each member's whole gene body (introns included) is aligned to the contig with
  `minimap2 -x asm20 -c --eqx -N 500 -p 0.01 --secondary=yes`.
  - Collinear hits are chained with `bench/guided_pipeline.py gene_body_chains` geometry: target gap ≤ Lq, span ≤ 2·Lq.
  - Every isoform's exon endpoints are lifted through the chain CIGARs.
  - `-p 0.1`, or Liftoff's `-p 0.5`, hides paralogues behind the self hit: the 30-kb KZFP bodies returned only it.
- **Candidates.**
  - A projection that touches a member's exons is a member hit. Anything else is a candidate.
  - The test must use exons: a 93-kb readthrough-widened member body swallowed an unexpressed neighbour.
  - Same-family candidates whose projected exons overlap form one candidate locus. Its model is the best chain's
    projection of the source member's full isoform set.
- **Read classes** (unique placements, UNIQ: primary + secondaries on the contig, best-AS set B per read), a partition
  applied in this order:
  - **X:** ≥ 1 sole read (expressed but unjoined, i.e. an O1 edge false negative);
  - **R:** more than half of the exon bp is RepeatMasker (the species RM track, a declared annotation input);
  - **(ii):** 0 sole reads but ≥ 1 tied read (an O2 twin);
  - **(i):** 0 sole and 0 tied reads, **the "no-read copy"**.
- **Status of a no-read copy.** It is a flagged class: never an O1 node, never a family member, never in O2's copy
  table. It is reference-present, so it is disjoint from O3's reference-absent class.

**Truth.** NEW = NAME ∪ LIFT95 ∪ LIFT01. It never uses O1's clause E (`rep2_C` §2):
- NAME is the RefSeq records whose NCBI description equals a seed record's (trailing "-like" / " pseudogene" removed).
- LIFT95 and LIFT01 are Liftoff 1.6.3 within-contig self-lifts with `-copies -sc 0.95`: literal, and with
  `-N 500 -p 0.01`.
- On NC_073244.2, NEW has 55 distinct copies (NAME 21, NAME and PROT 5, LIFT95 = LIFT01 29).
  - 46 are unread (42 of them in the CGB/NTF4-like tandem array MCL0).
  - 44 are class (i) (40 in MCL0).
- Human: NEW has 1 unread class-(i) copy on chr16 (MCL12) and 0 on chr20, and it confirms no flag on either contig. (The
  C+RG version's "human dev had no unread truth copies" was false; review3 L3.)

**What it found (NC_073244.2; conf/flags of class-(i) flags, then unread class-(i) copies found of 44).**

| arm | acceptance of a candidate | conf / flags | found | per-family C1 (good vs bad) | C2 vs literal Liftoff (win vs loss) | human chr16 bad families (chr20) |
|---|---|---|---|---|---|---|
| NR | family bar β(F) = the bottleneck of the maximum spanning tree of masked member-to-member projection weight | 45 / 48 | 42 | 3 > 1 | 2 > 1 | 0 (0) |
| O1C | O1's edge predicate verbatim on (member, projected copy) | 46 / 62 | 43 | 2 vs 8 | 2 vs 7 | 5 (1) |
| **O1C_L** | O1C with the copy's exon length set to its source member's | **45 / 49** | **42** | 2 > 1 | **1 = 1 (fails)** | **4 (0)** |
| O1C_L_FULL | the L form on full-model unions | 43 / 46 | 41 | – | – | – |
| LODN_lit | literal Liftoff seeded with the same member models | 35 / 36 | 34 | 1 > 0 | (reference) | 1 (0) |
| LODN_p01 | Liftoff with the projection's search | 41 / 42 | 39 | – | – | 1 (0) |

- **Every family-seeded search out-recalls literal Liftoff at the copy level, with c = 0.**
  - O1C_L vs LODN_lit: b = 8. That is 6 MCL0 array copies past Liftoff's `-N 50`, plus LOC101123748, LOC101130248 and
    LOC101152160.
  - O1C_L vs LODN_p01: b = 3. O1C_L vs NR: b = c = 0, with identical found sets.
  - Outside MCL0 there are 4 unread copies. NR confirms 3 of 5 flags and finds 2; both Liftoff runs flag and find 0.
  - The family calibration reaches divergent copies (the LILR-like MCL28 at 0.92-0.94 identity) that `-sc 0.95` cannot.
- **What the unconfirmed flags are.**
  - MCL20: an unannotated copy at 45.6 Mb (identity 0.92, coverage 0.84, 453 secondaries).
  - A family-attribution error: MCL3's extra flag sits on MCL5's truth copy.
  - On human chr16, O1C_L's 4 bad families are unannotated pericentromeric segmental-duplication copies at
    34.88-39.32 Mb, at 0.91-0.98 identity, each with 173-1,829 secondaries and no best-AS placement. They are
    unconfirmable by annotation, not shown to be false.
- **Class R is load-bearing on human.** It holds 115 (chr16) and 206 (chr20) O1C_L objects, and 0 confirmed objects or
  truth copies on gorilla in any arm. Without it, every unmasked arm, literal Liftoff included, would fail the human
  control by hundreds of flags.
- **No-read copies are not where O2's abstentions come from (dev).**
  - Of 1,325 reads with a secondary on a gorilla NR flag and their primary on a member, 0 are AS-tied (human: 0 of 99).
  - On chr16, 11,793 of the 11,961 member-touching AS-tied molecules (98.6%) have every twin on catalog members.
  - O1C_L's single class-(ii) object (MCL93, the HERC2P5 region) carries 42 of them.
- **Processed pseudogenes.** 109 intronless no-read candidates on gorilla: 103 unannotated, 1 confirmed. Genuine short
  fragments fail O1's 0.30 coverage clause.

**Why no rule is carried** (user decision 2, 2026-09-27).
1. **Every acceptance rule failed on dev.**
   - O1C is tautological on a projected copy (r1124). The copy's exons are the image of the chain being tested: the
     shared-exon ratio has median 1.000 (gorilla) and 0.978 (chr16), so only cov_longer is left, and chains reproducing
     2-5% of a member's exons pass it.
   - NR's family bar does no better than a permuted bar: 455 of 1,000 draws pooled, 405 of 1,000 per family (r1123).
   - O1C_L, the repair, was chosen after O1C failed. It fails C2 (1 = 1) and the human control on chr16 (4 > 1).
2. **Power.** 5-6 contributing families on dev, against a floor of 30. Only 2 of the 6 cast a C2 vote, and every
   verdict turns on MCL3, MCL20 and MCL28 (review3 H2).
3. **The judged truth leans toward Liftoff-like objects.** It confirms ≥ 0.95 identity or an exact NAME. So a divergent
   unannotated copy can never be confirmed, and C2 penalises the rule exactly where it is meant to add value
   (review3 M5). LIFT01 also shares the rule's own search (review3 M1).
4. **Cost.** Genome-wide C would take about 18-25 h of the lock (review3 H1). The 600 s call cap would also exclude the
   largest, array-richest contigs, which is non-random (review3 M4).

**What a later test of C must fix** (so the revisit starts here; review3 items):
- judge C1 and C2 on **decisive** families only (floors: #good + #bad ≥ 30; #dominates + #dominated ≥ 30). The dev
  decisive counts are C1 3 and C2 2 (H2);
- judge under NEW′ = NAME ∪ LIFT95, which equals NEW on dev (M1);
- state the scope as "within-contig, families of ≥ 2 expressed loci", and report the unsearched loci (M3);
- name the at-risk large contigs, and choose between exclusion and a longer bounded call (M4);
- report flags stratified by identity (≥ 0.95 vs < 0.95) and by NAME status (M5);
- apply floor (e) to the falsifiers C-F1 and C-F2. On gorilla dev their populations were 1 molecule and 1 object (L6);
- make the BAM secondary check operational: ≥ 1 FLAG-256 record, N recorded (L7);
- report O1C_L's side-flip count, which is 0 in 3,040 / 18,580 / 27,898 dev evaluations (L8);
- run gorilla (V3, V4) first, and skip the human controls if C is underpowered there (L9).

A held-out C test on V1-V4 would be a fourth use of those samples (§7.3).

**Frozen for the revisit:**
- `rep3_C/o1c.py` ce141fd4, `rep2_C/nr2.py` 100a9c00, `rep_noreads/noreads.py` 670087be, `rep2_C/run_lo.sh` cac2a7eb;
- outputs in `rep3_C/work/`, `rep2_C/work/` and `rep2_C/lo/`;
- `o1c.py score` reproduces `rep3_C/work/*` byte for byte (review3 check 2).

### 3.5 What is not claimed

- No held-out statement about R7, TS0, TS1, UC, NR, O1C, O1C_L or any other C arm. The A, B and C clauses, nulls,
  predictions and falsifiers of the earlier versions are withdrawn, not tested.
- The A and B instruments stay frozen on disk: `rep_edges/emu.py` fac9a560, `rep2_B/setans.py` d222e5a2,
  `summarize.py` 2e46738c, `summarize_real.py` 1918765f, `rep_o2/build_tables.py` 59788d72, `sim_iso.py` d5e7ce02.
  `emu.py` is still used by this file, for RG's freeze checks and reported rows only.

## 4. Development evidence for the carried arm (design only; selected on these contigs)

### 4.1 Arm RG (`rep3_RG.md`; frozen `mcl_families` 6cf8183f, `family_score` 27aa9445, `readthrough_eval.py` f6b99dcc)

**Split correctness S** (§8.1 definitions; every split `gene_id`; SEP / FRAG (PURE) / UNJ).

| input · rule | chr16 | chr20 | NC_073244.2 |
|---|---|---|---|
| BASE · old `rg.py` (junction key) | 84: 15 / **38** (32) / 31 → **fail** | 37: 4 / **22** (12) / 11 → **fail** | 8: 3 / **4** (3) / 1 → **fail** |
| BASE · **RG3** | 19: **16** / 2 (1) / 1 → **pass** | 10: **5** / 3 (1) / 2 → **pass** | 4: **4** / 0 / 0 → **pass** |
| BASE · donor\|acceptor | 51: 16 / 17 (15) / 18 → fail | 18: 5 / 7 (3) / 6 → fail | 4: 4 / 0 / 0 → pass |
| BASE · NULL_RG3 (matched to 19 / 10 / 4) | 19: 0 / 19 / 0 → fail | 10: 0 / 9 / 1 → fail | 4: 0 / 4 / 0 → fail |
| R3 · RG3 (descriptive) | 6: 4 / 1 / 1 → pass | 3: **1 / 1** / 1 → **fail (tie)** | 1: 1 / 0 / 0 → pass |

- **All dev values are below the 30-judged-split floor** (RG3: 18 / 8 / 4), so they are design evidence, not a pretest.
- **RG3's 5 dev fragmentations name the held-out failure modes.**
  - (a) The polish dropped the only bridge of one long gene: SRRM2 (520 / 44 reads), CEP250 (210 / 23), ANKRD11 (5′-UTR
    exons vs body).
  - (b) An annotated RefSeq readthrough record spans both parents, so a correct separation scores as a cut:
    FKBP1A-SDCBP2, PEDS1-UBE2V1. Dropping readthrough records (reported) makes chr20 7 vs 1.
- The old rule's 56 chr16 name-collision splits are 25 FRAG + 31 UNJ + 0 SEP (their representatives start at the same
  base). RG3 keeps 55 of the 56 together. All of RG3's SEPs are ghosts.

**Families, fused pairs, loci and ends.**

| contig | arm | fused PAIRS (cleared / new) | Compara TP/pred | referee TP/pred | largest (π) | pieces / in families | absorbed | TES / TSS genes | found_annotated |
|---|---|---|---|---|---|---|---|---|---|
| chr16 | BASE | 158 | 66/66 | 100/101 | 27 | – | 120 | 517 / 558 | 963 |
| | **RG3** | **138 (20/0)** | **87/87** | **125/126** | 29 (2) | 19 / 2 | 111 | 528 / 568 | 963 |
| | NULL_RG3 | 158 (0/0) | 66/66 | 100/101 | 27 | 19 / 1 | 119 | 520 / 558 | – |
| chr20 | BASE | 77 | – | 3/4 | 14 | – | 74 | 319 / 373 | 624 |
| | **RG3** | **69 (8/0)** | – | 3/4 | 14 | 10 / 0 | 68 | 324 / 380 | 624 |
| | NULL_RG3 | 77 (0/0) | – | 3/4 | 14 | 10 / 0 | 74 | 322 / 374 | – |
| NC_073244.2 | BASE | 63 | – | 275/275 | 21 | – | 65 | 574 / 680 | 869 |
| | **RG3** | **56 (7/0)** | – | 275/275 | 21 | 4 / 0 | 61 | 576 / 685 | 869 |
| | NULL_RG3 | 63 (0/0) | – | 275/275 | 21 | 4 / 0 | 65 | 576 / 682 | – |

- **The whole dev family gain is one truth family on chr16, which held-out excludes.**
  - NPIPB2 and NPIPB6, hidden inside the ghost-fused loci GSPT1~NPIPB2 and NPIPB6~EIF3CL, join NPIP.
  - Every gained TP pair is carried by an RG3 piece: Compara +21, referee +25.
  - Family sign b = 1, c = 0; the NULL gives 0 / 0.
  - RG3 matches the old rule on every truth: 0 TP and 0 FP pair difference on 3/3 contigs.
- **Fused pairs.** RG3 clears 35 of the old rule's 36 fused pairs. The one it keeps, PLK1~UBFD1, the old rule cleared
  only by cutting DCTN5, which lies between them. There are 0 new fused pairs on 3/3.
- **Loci.** +19 / +10 / +4 (+0.68 / +0.53 / +0.35%). The old rule added +4.0 / +2.3 / +0.8%.
- **Ends are now mostly attributable.** RG3 − NULL_RG3 is TES +8 / +2 / 0 and TSS +10 / +6 / +3. The old rule was
  within ±1 of its NULL on human. T stays safety-only.
- **Other rows.** Liftoff within-contig recall is unchanged (15/23, 1/1, 0/0). The identity gates pass on 3/3:
  - asserts 1-5;
  - every `c.*` row equal;
  - loci lost / gained 0;
  - 0 new fused pairs;
  - emu R0 byte-equal on 6/6 families runs;
  - relabel null zero-width 10/10.
- **Failure modes the families clauses can still show held-out:** a ghost piece joining a large family wrongly (P1, H),
  or a true pair lost through MCL re-flow (P2). Dev has 0 of either for RG3. Their dev negative controls are borrowed
  from arm A: R1 (Compara pP 1.000 → .753, largest 27 → 48) and R3c (pS .342 → .332).

### 4.2 The incremental families instrument (`rep4_incr.md`; scratch `rustle_figures_dev/rep4_incr/`)

**Acceptance on dev** (`clusters.tsv` byte-identical to the full `fam_call.sh`-style run; every run under the lock):

| dev arm | mode | `clusters.tsv` vs full | `copies.tsv` | family graph vs full | wall: full → incr |
|---|---|---|---|---|---|
| chr16 RG3 | incremental | **identical** | identical | same 5,877 edges; 2 reweighted (max Δ 0.0042) | 52.8 → 14.5 s |
| chr20 RG3 | incremental | **identical** | identical | identical (688 edges) | 42.9 → 7.0 s |
| NC_073244.2 RG3 | incremental | **identical** | identical | identical (430 edges) | 42.7 → 9.0 s |
| chr16 NULL_RG3 | fallback (`mid_occ` 139 → 141) | **identical** (PAF byte-identical) | identical | identical | 56.8 → 52.3 s |
| chr20 NULL_RG3 | fallback (93 → 94) | **identical** (PAF byte-identical) | identical | identical | 44.0 → 42.9 s |
| NC_073244.2 NULL_RG3 | fallback (139 → 140) | **identical** (PAF byte-identical) | identical | identical | 50.5 → 31.3 s |

- **Speed.** RG3 on dev: 30.5 s vs 138.4 s (4.5×). NULL_RG3 gets no speedup (fallback).
- **The PAF is not byte-identical on the RG3 arms.**
  - It differs on 42,613 / 13,315 / 4,497 lines, almost all in the `rl:i`, `cm`, `s1` and `de` tags or in line order.
  - The family view (non-self lines with ≥ 300 bp at ≥ 0.70 identity) differs on 45 / 4 / 43 lines of 83,265 / 35,545
    / 13,023.
- **Why byte-identity is impossible with any useful speedup.**
  - Share of the aligned work needed to recompute, by oracle: 0.88-0.98 for a byte-identical PAF; 0.25-0.92 for an
    identical line multiset; 0.02-0.70 for the floor of any reuse scheme.
  - Fixing `--mid-occ` does not help: the rescue ranking, tie order and tie hash remain, and it would change the
    reference run.
- **Soundness of A_chain.** 0 missed pieces-touching queries on 6/6, at 4-7% of the full run's cost.
  - The predictor "BASE line to the split parent" misses 222 on chr16 NULL_RG3, where the clusters then differ.
  - A seed-anchor predictor is sound but costs 0.85-0.89.
- **Target-mask patch fidelity.**
  - The masked lines equal the full run's A_chain → piece lines on 6/6 arms: byte-identical on 4/6, and on the two
    chr16 arms identical except column 12 (MAPQ), which families do not read.
  - Unmasked, the patched binary reproduces the full PAF byte for byte on 2 arms.
- **Forced incremental on NULL_RG3** (superseded v2 wrapper 96df0dfa, `incr_v2forced/`): `clusters.tsv` identical 3/3.
  But on chr16, 24 of 5,869 edges were lost and 2,422 reweighted (max Δ 0.395), and `copies.tsv` `max_family_identity`
  differs for 4 copies.
  - Cause: the NULL's pieces overlap in span, duplicating body sequence (chr16 +1.0 Mb, +1.8%). That shifts every
    query's seed use, which is why the frozen wrapper falls back when `mid_occ` moves.
- **Not done** (IRG7): budget-sharded passes (genome-wide, a chaining or masked pass is about 5-12 min on
  gorilla_OR6737, so one pass can exceed a 600 s call), and `fam_call_incr.sh`.

## 5. Arms

| arm | what | status |
|---|---|---|
| **BASE** | the shipped families stage on the BASE GTF, `rt_bin_frozen/mcl_families` 6cf8183f via `fam_call.sh` (799f3e71) | reused for human_testis and chimp_PTR; built by queue 1 for the other four |
| **RG3** | §2.1, then the families stage by F_INCR or F_FULL (§2.2) | **new, judged** |
| **NULL_RG3** | §5.1, then the same families method as RG3 | **new, judged in FU, S (as a control) and G** |
| old rule `rg.py` 92892fb3 | the GTF rewrite and S only, no families | reported (RG-F2) |

### 5.1 NULL_RG3

`rg3_null.py` 6ece9ae088766f8c22ecc251c40a883d1451c426 is `rg_null.py` 08b5c317 verbatim, with RG3's pieces and
`.rg<k>` names.

- **Seed.** Per contig, Python `random.Random` seeded with the string `"20260926:<sample>:<contig>"`, so the draw does
  not depend on PYTHONHASHSEED. Every list is sorted by (first exon start, `gene_id`) before a draw.
- **What it splits.** NULL_RG3 leaves RG3's split `gene_id`s intact. For each of them (g, in order) it splits one donor
  `gene_id` h that RG3 leaves as ONE piece and that holds ≥ 2 transcripts. Donors are drawn uniformly without
  replacement: one split per donor.
- **T_g (review3 L1).** For a split g, its non-keeper pieces in RG3's naming order carry n_p transcripts each.
  **T_g = Σ n_p**, the number of transcripts RG3 moves out of g. b_g = the maximum over those pieces of
  floor(log2(max(1, reads of the piece representative))).
- **Donor.** cand(h) = h's spliced transcripts other than rep(h) = max (reads, span, −index). The donor must have
  |cand(h)| ≥ T_g, including ≥ 1 transcript in depth bin b.
  - The bin search is b = b_g, b_g − 1, b_g + 1, b_g − 2, …, with the lower bin first on a tie.
  - If no bin has a donor, any donor with |cand(h)| ≥ T_g is used, and the draw is counted as bin-unmatched.
  - If there is none, a shortfall is counted, not replaced.
- **Pieces.** moved = one seed transcript in bin b, plus T_g − 1 others drawn from cand(h). They are filled into pieces
  of sizes n_1, n_2, … named `<h>.rg2`, `<h>.rg3`, …
- **Matched:** split `gene_id`s per contig, pieces per split, transcripts per piece, depth bin, spliced. On dev: 19 / 10
  / 4 splits and 171 / 30 / 12 transcripts, all at the exact bin.
- **Not matched, reported:** the donor's own depth; exon disjointness. NULL pieces overlap their keeper by construction,
  which is what makes the NULL a predicate control for S and FU.

## 6. Implementation gates (development contigs only; all before any held-out arm run, except IRG7-H)

**No tracked source is edited for this test.** Every step runs through frozen instruments on top of frozen binaries:
- `rt_bin_frozen/`: `mcl_families` 6cf8183f51ce6d56790dfd2f666823287e63c7c1 and `family_score`
  27aa94454638718ded0cd35fe31a87334ec7a372. On dev the three `mcl_families` builds (6cf8183f, 0c4639e2, 4d0b0817)
  give byte-identical outputs (`rep2_RG` §5).
- Scripts:
  - `bench/mechanism/readthrough_eval.py` f6b99dcc1a3f972427cf67bc3de3b846dd5c5f30;
  - `fam_call.sh` 799f3e71;
  - `tools/rustle_pipeline.sh` d1816c35;
  - `tools/mm2_shard.sh` 04c2c639. It is tracked and live, so the arm runner checks its sha1 (and
    `tools/rustle_pipeline.sh`'s) before every families call.
- The instrument of §2.2 (IRG7).
- The referee: `bench/truth.py` 2c59a2d8. Queue 2 runs frozen copies in `rustle_figures/rep_prereq2/frozen/`:
  `truth.py`, `lib.py` e7765bc0, `mcl_port.py` 8d0e9423 and the `mcl_port` binary 99749679.
- The genes-table builder (IRG8).

**A Rust port of RG happens only after its verdict**, behind an opt-in switch.
- Its unset output must be `cmp`-identical to the binary it replaces, and its set output `cmp`-identical to the
  instrument's, on dev and on every held-out substrate.
- The port plan is `rt4_ghost.md`, "Proposed Rust change", with RG3's adjacency.
- A port would run the full all-vs-all (F_FULL). That is why IRG7-H ties F_INCR to F_FULL.

**Any failure is fixed in the instrument, never in the rule.** A fixed instrument re-passes every dev gate of its arm,
byte for byte, before it runs again.

- **IRG1 (RG3 and its NULL).**
  - `rg3.py` ec17e540 must reproduce `rep3_RG/rg/*.RG3.gtf` (BASE and R3 × 3 contigs) byte for byte, and be idempotent.
  - `test_rg3.py` 03e10572b144bf8f05fc26048cfe9a871aefae3f must pass its 19 fixtures.
  - `rg3_null.py` 6ece9ae0 must reproduce `rep3_RG/null/*.NULL3.gtf`, and be `PYTHONHASHSEED`-independent (two seeds,
    `cmp`).
- **IRG2 (split classes).** SHA1_CLS is generalised from `rep3_RG/cls.py` 413390944e4a7e54213256ea649a31802283c135,
  with the species genes table and GFF as arguments. It must reproduce:
  - review2's table for the old rule: 84 = 15 / 32 pure / 6 / 31; 37 = 4 / 12 / 10 / 11; 8 = 3 / 3 / 1 / 1; NULL
    clean / pure 1 / 68, 0 / 31, 0 / 7;
  - §4.1's S table (`rep3_RG/cls/cls.json`).
- **IRG3 (RG scorer).** SHA1_RG_EVAL is generalised from `rep2_RG/rg_score.py` fea2167bfad967624f5f9a9b79a904a7721a72cb
  and `rep3_RG/score3.py` 792e8fa0abcd21f1641aa2f8c7cb92b548403284. It must reproduce `rep3_RG/scores.json`: FU, the
  referee rows, family signs, H, carriers and revealed copies.
- **IRG4 (complete-universe referee).**
  - SHA1_REP_PAIRS, restricted to the truth universe, must reproduce `family_score` 27aa9445's pairwise TP and predicted
    counts exactly on chr16, chr20 and NC_073244.2, for BASE, RG3 and NULL_RG3.
  - Its complete-universe dev rows, and the P1_all rows (§8.1), go into Amendment 1 as found (§0 item 3).
- **IRG5 (emulator).**
  - `emu.py` fac9a560 as R0 must reproduce byte for byte the six `clusters.tsv` of `rep3_RG/fam/` (RG3 and NULL_RG3 × 3
    contigs) and the three BASE ones of `rep2_RG/fam/`.
  - 10 random node relabellings must give the same partition.
- **IRG6 (the dev-contig restriction path, review3 L5).**
  - Pool the chr16 and chr20 dev products of BASE, RG3 and NULL_RG3: the GTFs concatenated, and the `clusters.tsv`
    files merged with family ids prefixed by contig.
  - `readthrough_eval.py` on the pooled GTF, with `--drop-contigs` = every annotated contig except chr16, must give a/b/d
    rows equal to its `--contigs chr16` rows.
  - SHA1_RG_EVAL and SHA1_REP_PAIRS, with the §7.2 restriction dropping chr20, must give pair rows equal to their
    chr16-only rows.
  - On a fixture copy in which one chr20 locus is appended to a chr16 family, exactly that locus's pairs must be
    removed.
- **IRG7 (the families instrument; user decision 4).** Before Amendment 1 the instrument gets the two pieces
  `rep4_incr` §6 lists:
  - every minimap2 pass routed through `tools/mm2_shard.sh paf` (with `MM2_SHARD_MINIMAP2 = bin/minimap2-tmask` and
    `MM2_TARGET_MASK` exported for the masked pass), so that each pass fits one 600 s call and resumes;
  - `fam_call_incr.sh <s> RG|NULLRG`, with `fam_call.sh`'s contract: 0 = the clusters exist, 75 = a budget stop,
    resumable. `INCR_WORK` is never deleted, and `INCR_BASE_MMI` is resolved from the md5 of BASE's `loci.fa`.

  The resulting SHA1_INCR and SHA1_FAM_INCR must pass all of the following, `cmp` byte for byte, on the six dev arms
  {RG3, NULL_RG3} × {chr16, chr20, NC_073244.2}:
  - **(a)** `clusters.tsv` and `copies.tsv` equal to the full runs in `rep4_incr/full/`, on 6/6. (`rep4_incr` met this
    with fed8d80e, unsharded.)
  - **(b)** with `INCR_FORCE_INCR=1` on the three NULL_RG3 arms, `clusters.tsv` equal to the full run, on 3/3. So far
    this has been shown only with the superseded v2 wrapper. It is required because genome-wide the fallback is not
    expected to fire (§2.2).
  - **(c)** the sharded passes' PAF equal to the unsharded fed8d80e PAF on 6/6. This includes, for each arm type, one
    run whose budget is small enough that ≥ 1 pass stops with 75 and resumes.
  - **(d)** the fallback PAF equal to the full run's on the three NULL_RG3 arms. (`rep4_incr` met this.)
  - **(e)** emu R0 on the instrument's PAF reproducing its own `clusters.tsv`, on 6/6.

  **If any of (a)-(e) fails, or if the pieces are not written when Amendment 1 is due, the method is F_FULL on every
  substrate.** Amendment 1 then says so, and the lock time is §13.1's F_FULL line. The PAF is not a gate object: it is
  not byte-identical, and cannot be (§2.2).
- **IRG7-H (the held-out cross-check; the first arm step, §13 step 4a).**
  - On V2 (human_testis, the cheapest full run: about 21 min per arm), RG3 and NULL_RG3 are run by both F_INCR and
    F_FULL. For each arm, the two `clusters.tsv` are compared with `cmp`. No score, no row and no other product is read.
  - **Identical for both arms:** F_INCR is the method on V1 and V3-V6.
  - **Any difference:** F_INCR is deemed failed on held-out, and every substrate uses F_FULL. The F_INCR products are
    moved aside unscored. The only thing reported is the number of differing lines.
  - **Disclosed weakness:** V2 has the fewest isoforms per locus (1.9, §0), so it is probably the instrument's easiest
    case. Every substrate's `report.json` is reported.
- **IRG8 (genes-table builder, review3 L11).**
  - The builder is `figures/_o1.py annotation_cache` (e3994c68; with `figlib.py` 3bfd0414, `samples.py` 52e87f3e,
    `assembly.py` 94861378), driven by `rustle_figures/rep_prereq2_genes.py` 52d3f70b.
  - It must rebuild the human and gorilla tables (`genes.tsv`, `genes_only.gff`, `mrna_models.tsv`; for human also
    `exons.gtf`) byte-identical to `families_gw/species/{human,gorilla}/`, which are the tables S was developed with.
  - Queue 2 runs this gate first, and its helper refuses to build the chimpanzee and orangutan tables otherwise. On
    failure, V5 and V6 are undecided (§9).
- **Freeze (Amendment 1).** Amendment 1 records:
  - SHA1_CLS, SHA1_RG_EVAL, SHA1_REP_PAIRS, SHA1_INCR and SHA1_FAM_INCR (or "F_FULL", with the reason);
  - every command line: `rg3.py`, `rg3_null.py`, `fam_call_incr.sh` / `fam_call.sh`, `readthrough_eval.py score`,
    `truth.py protein-homology`;
  - the dev tables of IRG2, IRG4 and IRG7;
  - the user's acceptance (header items 1-6);
  - this file's sha1.

**Freeze checks on held-out BASE** (after Amendment 1, before any arm). They read only BASE products:
- **(a)** emu R0 reproduces the substrate's BASE `fam.clusters.tsv` byte for byte. If it fails, the reported rows with
  dev-contig nodes deleted (§7.2) are not produced on that substrate; nothing else changes.
- **(b)** 10 random node relabellings of BASE give BASE's partition. If this fails, RG's family rows (P1, P2, H, G) are
  not measured on that substrate, because MCL tie order would confound them.
- **(c)** BASE's restricted rows (§7.2) are computed and written to Amendment 2. Unrestricted, they must reproduce the
  already-published `rt_arms/tables_v3/` rows (§0).
- **(d)** The instrument's BASE inputs exist: `<s>.BASE.fam.loci.{fa,paf}` and BASE's index or `mid_occ`. BASE's
  `mid_occ` is recorded. (The C+RG version's check (d), on the BAMs' `@PG`, served arm C only and is removed.)

## 7. Substrates, truths, exposure

### 7.1 Substrates (each judged on its own; species never pooled)

| id | sample | contigs scored |
|---|---|---|
| V1 | human_A119b | genome minus chr16, chr20 |
| V2 | human_testis | genome minus chr16, chr20 |
| V3 | gorilla_OR6737 | genome minus NC_073244.2 |
| V4 | gorilla_KB3781 | genome minus NC_073244.2 |
| V5 | chimp_PTR | whole genome |
| V6 | orangutan_PPY | whole genome |

- RG runs on all six.
- V3 and V4 share one genome and its truths, and so do V1 and V2. Only the libraries differ.
- **V5 and V6 carry no dev-contig name.** Their orthologues of human chr16/chr20 are scored. The design never looked at
  those genomes, but it did look at the human NPIP region, so this is disclosed.

### 7.2 The dev-contig rule (every sample)

- chr16 and chr20 (human samples) and NC_073244.2 (gorilla samples) are removed from every scored object:
  - truth pairs and predicted pairs with a locus on them;
  - fused pairs, splits (S), universe genes and reference loci on them.
- **How.**
  - `readthrough_eval.py score` gets `--drop-contigs chr16,chr20` on V1 and V2, and `--drop-contigs NC_073244.2` on V3
    and V4. It gets no flag on V5 and V6.
  - SHA1_CLS, SHA1_RG_EVAL and SHA1_REP_PAIRS apply the same sets.
  - Gate IRG6 checks the path.
- **Inputs are not restricted.** RG's families are built genome-wide, so a piece on a removed contig can still join a
  family; its pairs are then removed from scoring.
- **Reported, not judged:** RG's family rows with dev-contig nodes deleted before MCL (emu). A held-out pair merged
  through a chr16 hub shows up there.
- **What is lost.** NPIP-U2 lies on chr16 and is scored nowhere. So is RG's entire dev family gain (§4.1).

### 7.3 Exposure

- **This is the third verdict on the six samples**, after the R prereg (09-25) and v3 (09-26). It is the first on node
  regrouping.
  - v3 §0 item 5 requires a new library or the user's explicit acceptance for a third reuse.
  - The user's decision is recorded in the header. Amendment 1 records the user's confirmation before §13 step 4's arm
    runs.
- **After this study, V1-V6 are SPENT for node-regrouping work.** Any fourth verdict on them needs a new library or the
  user's explicit acceptance. That includes a held-out test of the descriptive arms A, B or C.
- **What was seen in advance:** listed in §0, including the published BASE rows of `tables_v3`.

### 7.4 Truths

| truth | substrates | role |
|---|---|---|
| **Protein referee, complete universe.** See the definition below the table | V1-V6 | **judged:** pairwise tp and precision |
| Protein referee, `family_score` form (27aa9445, `--chrom ALL --pairwise`) | V1-V6 | reported (Fig. 7 continuity) |
| **Compara Primates** pairs (`families_gw/species/human/compara.Primates.families.tsv`; `family_score` 27aa9445 via `readthrough_eval.py` f6b99dcc, pairwise, tie-free, r1045) | V1, V2 | **judged:** tp and precision (an upper bound) |
| **Liftoff copy pairs**, Fig. 6s rule (§2.3). The denominator is fixed | V1-V6 | **judged**, recall only |
| **Fused gene pairs, ends, loci** (`readthrough_eval.py` f6b99dcc: `fused_class` in pair form; TES/TSS; `d.found_annotated`; absorbed genes) | V1-V6 | **judged** |
| **Species genes table** (`genes.tsv` records `Name\|start1`, strand from `genes_only.gff`; for chimpanzee and orangutan built by queue 2 after IRG8) | V1-V6 (S; `gene_at` for the referee) | **judged** |
| Soto 2025 families; NPIP-U2 | V1, V2 | descriptive. Soto is "not independent" (house rule); NPIP-U2 is not scorable (§7.2) |

**The complete-universe protein referee.**
- It is `truth.py protein-homology --chrom ALL` (the §6ko rule, not re-tuned), run on each species' registry RefSeq GFF
  and genome. Human uses `/mnt/linuxdisk/tmp/regress/chm13.gff`, the full CHM13 RefSeq GFF, never `HSA_genomic.gff`.
- It gives `PF` families over (contig, Name) genes. Loci resolve to genes by `gene_at` on the species' genes table.
- **Universe** = every gene of the referee's `<prefix>.proteins.faa`: every non-excluded protein-coding gene, singletons
  included.
- A predicted pair is **judgeable** iff both loci resolve to universe genes, and **TP** iff both genes are in one `PF`.
- tp equals `family_score`'s. pred also counts the pairs with a singleton or cross-family gene, which `family_score`'s
  universe intersection deletes (r770/r991).
- **Where it lives.** The prefixes are `rustle_figures_dev/truth/gw/{human,gorilla,chimp,orangutan}/ph`
  (`figures/_o1.py` PH_DIR_DEFAULT; queue 2). The directory's name is historical, but its content is genome-wide
  truth: no dev-only author reads `ph.families.tsv`, `ph.edges.tsv` or the shards.

**Known weaknesses, stated now:**
- Compara and Liftoff pairs are positives-only or incomplete. Compara precision is an upper bound.
- The complete-universe referee cannot see non-coding or unannotated loci, so an arm that adds unlabelled members is
  flattered (r991).
  - §8.1 reports every arm's unjudgeable predicted pairs, and P1_all.
  - On dev, RG3 pieces joined families 2 times out of 33, and UNJ was 1 / 2 / 0 (review3 M7).
- `gene_at` under-credits passengers in fused loci (09-21 trap). The bias is identical across RG's arms.
- **S is judged against RefSeq.**
  - Annotated readthrough records (human) count against a correct separation. They are kept, which is conservative.
  - Unannotated pieces fall in UNJ, which is reported as a count at a stated truth completeness (09-22 trap).
- **F_INCR's kept pairs are BASE's** (§2.2). A family effect that a full run would create only through the global
  seeding shift is absent from F_INCR. IRG7-H checks on V2 whether that changes any cluster.

## 8. Metrics and clauses (integer arithmetic; per substrate; dev-contig rule applied)

### 8.1 Arm RG (zero tolerance; no new constant)

Notation: B = BASE, N = NULL_RG3, fu = fused gene pairs.

**Split classes** (SHA1_CLS).
- A **split** is an input `gene_id` whose transcripts carry ≥ 2 output `gene_id`s (pieces).
- A(P) = the annotated genes of the species genes table with an exon that shares ≥ 1 base, on the same strand, with an
  exon of some transcript of piece P.

Every split falls in exactly one class:
- **FRAG** (fragmentation): some annotated gene lies in A(P) for ≥ 2 pieces (one gene given two `gene_id`s). PURE (every
  piece has A(P) = {x}, review2's "pure fragmentation") is a subset, reported.
- **SEP** (clean separation): no gene is cut, and every A(P) is non-empty, so the pieces own pairwise-disjoint gene
  sets.
- **UNJ**: no gene is cut, and some piece touches no annotated exon on its strand. Reported, never judged.

| clause | measure | RG passes on s iff | floor |
|---|---|---|---|
| **FU** existence of separations | fused gene PAIRS: distinct unordered annotated gene pairs (same strand as the locus representative, non-overlapping spans) whose exon unions one spliced locus overlaps (`fused_class`). Pairs, not loci, so a split cannot raise the count (the `a.fused` trap) | fu(RG) < fu(B) **and** fu(RG) < fu(N) | fu(B) ≥ 50 |
| **S** split correctness | SEP and FRAG over RG3's splits | SEP > FRAG (strict; a tie fails) | SEP + FRAG ≥ 30 |
| **P1** precision | the complete-universe referee (V1-V6); Compara (V1, V2) | tp(RG) · pred(B) ≥ tp(B) · pred(RG) | ≥ 30 truth families |
| **P2** recall | the same truths' tp; Liftoff pair recall (V1-V6) | tp(RG) ≥ tp(B); rec(RG) ≥ rec(B) | ≥ 30 families; ≥ 30 read-supported pairs |
| **H** hub | L = the largest family's loci; π = RG3 pieces among its members | L(RG) ≤ L(B) + π | – |
| **T** ends | TES- and TSS-recovered genes (fixed universe) | tes(RG) ≥ tes(B) **and** tss(RG) ≥ tss(B) | universe ≥ 1,000 |
| **L** loci | `d.found_annotated`; absorbed genes | found(RG) ≥ found(B); absorbed(RG) ≤ absorbed(B) | ≥ 1,000 reference loci |
| **G** attributable family gain (label only) | on a judged recall truth: the pair gain, and the family-level sign (b = truth families whose TP pairs rise, c = those that fall) | tp(RG) − tp(B) > tp(N) − tp(B) **and** b > c | as P2 |

- **Identity gates.** These are not clauses: a failure is a bug, and it stops the substrate, whose RG clauses become not
  measured.
  - `rg3.py` asserts 1-5;
  - every `c.*` row equal to B's;
  - `a.loci_lost_vs_base` = `a.loci_gained_vs_base` = 0;
  - **0 new fused pairs vs BASE**. A split-only rewrite cannot fuse (review2 X2 fix 4);
  - 0 junction keys shared by two output `gene_id`s. This is structural: 0 on the BASE and RG3 dev GTFs, 3/3
    (review3 check 5);
  - the relabel null on RG3 zero-width;
  - emu R0 reproducing each arm's `clusters.tsv` from its own PAF.
- **Why S and FU together (rewritten after review3 M2).**
  - FU is an existence check: it passes iff at least one fused pair is cleared beyond the NULL, and NULL_RG3 cleared 0
    on 3/3 dev contigs.
  - A clean separation between span-disjoint genes clears a pair. So FU is implied by S whenever S passes, except when
    every SEP is between span-overlapping genes. It adds a failure path only where S is not judged.
  - FU is not a recall: clearing 1 of 500 pairs passes it. The recall side is carried by S's floor (≥ 30 judged splits)
    and by EFFECTIVE's requirement that S be judged on ≥ 1 human and ≥ 1 ape.
  - **So RG's EFFECTIVE rests on one informative clause, S, plus the non-inferiority and safety clauses P1, P2, H, T
    and L.** That is stated, not hidden.
  - S is the only clause that sees fragmentation. FU, P1, P2, T and L cannot: both fragments resolve to one gene, and
    extra ends only add.
  - **S counts any cut gene, not only PURE.** A split that separates two genes but cuts a third (the old PLK1 / DCTN5 /
    UBFD1 split) is an error.
  - A zero-tolerance form ("0 PURE fragmentations of genes with ≥ k reads") fails RG3 on dev for every k ≤ 44, so any
    passing k would be a free constant.
- **T is near-identity for a split-only rewrite.** Every transcript end survives in some piece. It is a safety clause,
  and no dev arm fails it (review3 §1).
- **Why G is only a label.**
  - r1017 bounds a perfect over-merge repair at +0.021 referee F.
  - The whole dev family gain is on chr16, which held-out excludes (§4.1).
  - FP pairs scale with family size: one wrong gene in a 40-member family adds up to 39 FP pairs. So P1 is RG's likeliest
    real failure.
- **Reported, never judged:**
  - PURE; UNJ; S with RefSeq readthrough records dropped; S of the old rule `rg.py` (RG-F2);
  - **the cleared fraction** (fu(B) − fu(RG)) / fu(B), with no bar (review3 M2);
  - **P1_all** (review3 M7): pred_all counts every predicted pair of the arm, with unjudgeable pairs counted as not TP;
    tp(RG) · pred_all(B) ≥ tp(B) · pred_all(RG) is reported beside P1, with each arm's unjudgeable pairs;
  - **revealed copies:** truth-labelled genes that are family members only through an RG3 piece and share a truth family
    with another member of their catalog family. They are counted per piece, so family size does not multiply them.
    Dev: {NPIPB2, NPIPB6} on chr16, 0 elsewhere and under both NULLs;
  - `a.fused` and `a.loci` (node-level, traps named); fused genes; reciprocal one-to-one loci; bipartite F (scipy tie
    policy, r1045);
  - Soto and NPIP (unrestricted, marked "includes dev contigs");
  - pieces per family; gained TP pairs by carrier (piece, keeper, unchanged); the name-based ghost / collision split;
    RG − N for T and L; the FRAG mechanisms (a) and (b) where computable;
  - **the families method per substrate:** F_INCR or F_FULL; each run's mode, `mid_occ` (arm / BASE), N, A_chain and
    seconds; the pairing of RG3 and NULL_RG3 modes.

### 8.2 Reported for every substrate (never judged)

- V1 and V3 recomputed without the Fig. 7 contigs: V1 without chr2, chr6, chr8 and chr10; V3 without NC_073234.2.
- RG's rows with dev-contig nodes deleted before MCL (emu).
- Every arm's unjudgeable predicted pairs, per truth.

### 8.3 Can every judged clause fail? (the negative control of each)

| clause | the outcome that fails it | shown on dev, or plausible held-out |
|---|---|---|
| S | SEP ≤ FRAG | **dev:** the old rule 3/3, donor\|acceptor on human, NULL_RG3 3/3, **RG3 on the chr20 R3 GTF (1 = 1)**. 5 of RG3's 33 dev splits are FRAG. Held-out mechanisms (a) and (b) of §4.1 |
| FU | S not judged, and no span-disjoint separation beyond the NULL | **dev:** NULL_RG3 as the tested arm (0 cleared). Held-out: a library whose polish drops few bridges (few ghosts), e.g. testis at 1.9 isoforms per locus (§0, seen). Implied by S wherever S passes and some SEP is span-disjoint (§8.1) |
| P1 | a piece joins a large family wrongly | dev control borrowed from arm A: R1 (host-attached pieces) at pair level on chr16 (Compara pP 1.000 → .753). Plausible on any NPIP/TBC1D3-like array |
| P2 | a true pair lost through MCL re-flow | dev control from arm A: R3c (Compara pS .342 → .332). Plausible |
| H | the largest family grows beyond its RG pieces | dev control from arm A: R1, 27 → 48 with π = 0 |
| T | a gene's end carried only by a fragment owned by another gene | near-impossible for a split-only rewrite; no dev arm fails it; a safety clause |
| L | a FRAG divides a gene so that no piece covers ≥ 50%, or absorbed genes rise | **dev:** the old NULL_RG (found_annotated 963 → 962 on chr16) |
| G (label) | no gain beyond the NULL, or b ≤ c | **dev:** NULL_RG3 (b = c = 0); chr20 and gorilla for RG3 |
| identity gates | a new fused pair, a differing `c.*` row, a failed assert | cannot fail for a correct split-only rewrite; each would be a bug |

## 9. Verdict (dev never enters; species never pooled)

**Common rules** (review2 M2).
- **Not judged:** a truth or clause below its power floor on s. It neither passes nor fails on s.
- **Not measured:** a clause that cannot be computed within the machine rules, or that a freeze check or identity gate
  disables. It caps the arm at KEEP OPT-IN, unless REFUTE already holds.
- **Undecided:** a BASE prerequisite of the substrate does not exist (review3 L4). These are:
  - its BASE families (queue 1);
  - the genome-wide protein referee of its species (queue 2);
  - for V5 and V6, the species genes table, which S and the referee's `gene_at` need (queue 2, after IRG8).

  Amendment 2 records any undecided substrate before any arm command. The user may then accept, by amendment and before
  any arm number exists, that the arm proceeds with the affected clauses not measured (capped at KEEP OPT-IN).
  Otherwise, with an undecided substrate the verdict is **undecided**, unless REFUTE already holds on the other
  substrates.
- **The families method never changes a clause.** Under F_FULL (IRG7 or IRG7-H failed), a families run that cannot
  finish within the machine rules makes P1, P2, H and G not measured on that substrate.

**Arm RG.**
- P(s) = P1 and P2 on every judged truth of s, and H, T and L where judged.
- P is **measured** on s iff ≥ 1 precision truth and ≥ 1 recall truth are judged on s.
- F_RG = the set of substrates on which P, FU or S fails.

| verdict | condition |
|---|---|
| **EFFECTIVE** (default candidate; the flip is the user's call) | F_RG empty; P measured on 6/6; FU and S each judged on ≥ 1 human and ≥ 1 ape. Labelled **"with family gain"** iff G holds on ≥ 1 human and ≥ 1 ape, else **"no attributable family gain"** |
| **REFUTE** | \|F_RG\| ≥ 2 |
| **KEEP OPT-IN** | \|F_RG\| = 1; or F_RG empty but P is not measured on some substrate, or FU or S is judged on no human or on no ape |

- Failures of different clauses on different substrates add up in F_RG.
- An identity-gate failure makes every RG clause not measured on that substrate.

## 10. Predictions (before any held-out number; set by this file's authors after `rep3_RG` and `rep4_incr`)

**Arm RG.** "Post-exposure" marks the predictions whose baselines were seen in the `tables_v3` BASE rows (§0).
1. The `rg3.py` asserts pass on 6/6 (p 0.95). Pieces add 0.2-1% loci, with A119b the most and testis the fewest
   (post-exposure).
2. S is judged on ≥ 5 of 6 (p 0.7) and passes wherever judged (p 0.65). Its likeliest failure is a human sample with
   many annotated readthrough records.
3. FU on 6/6 (p 0.8): a fused-pair reduction of 5-12%, against a NULL_RG3 reduction of ≤ 1% (post-exposure).
4. P1 on 6/6 (p 0.55); the likeliest failure is a ghost piece joining a large family. P2, H, T and L on 6/6 (p 0.8)
   (post-exposure).
5. G on ≥ 1 human and ≥ 1 ape (p 0.25). NULL_RG3's family rows equal BASE's within ±2 pairs per truth (p 0.8). Revealed
   copies ≥ 1 on ≥ 1 substrate (p 0.5). (Post-exposure: the testis and PTR BASE family rows were seen.)
6. RG3 − NULL_RG3 ≥ 0 for both TES and TSS on ≥ 4 of 6 (p 0.6) (post-exposure).
7. **Verdict:** EFFECTIVE with family gain ~0.05; EFFECTIVE with no attributable family gain ~0.20; KEEP OPT-IN ~0.40;
   REFUTE ~0.25; undecided ~0.10.
   - **Calibration (review3 L10).** The product of the component predictions under independence is
     0.55 × 0.8 × 0.65 × 0.8 ≈ 0.23, before the "S judged on a human and an ape" and "P measured on 6/6" requirements.
   - The components are positively correlated: P1, P2 and H fail together through one wrong join, and FU follows S. So
     EFFECTIVE (0.25 in total) sits a little above the product, not at the C+RG version's 0.38.
   - Undecided rose from 0.05 because the referee and the chimpanzee/orangutan tables do not exist yet.

**The families instrument.**
8. IRG7 passes on dev with the sharded passes, in forced mode included (p 0.85).
9. `mid_occ` is 500 (the asm clamp) in BASE and in both arms on 6/6 substrates, so NULL_RG3 takes the incremental path
   (p 0.8).
10. IRG7-H: `clusters.tsv` identical on V2 for both arms (p 0.5).
11. The lock time falls within §13.1's expected line (p 0.35), i.e. F_INCR is used and no NULL_RG3 falls back.

## 11. What would falsify the design reasoning (reported whatever the verdict)

- **RG-F1 "RG's family gain is ghost pieces joining their true family."**
  - It is measured on a substrate where RG3 gains TP pairs on a judged recall truth.
  - It is falsified if fewer than half of the gained TP pairs are carried by an RG3 piece.
  - n (gained pairs, and distinct carrier pieces) is reported beside it. Below 30 gained pairs (floor (e); review3 L6)
    it is reported as underpowered and cannot falsify.
  - Dev chr16: all of them (21 Compara / 25 referee pairs from 2 pieces, below the floor).
- **RG-F2 "The exact-junction adjacency is what fragmented genes."**
  - It is falsified if the old rule `rg.py` 92892fb3 passes S on at least half of the substrates where its S is judged.
  - Dev: it fails 3/3.
  - It costs a GTF rewrite and SHA1_CLS, with no families run.
- **The instrument, reported:** on every substrate that used F_INCR, the family-view line difference is not measurable
  without a full run. Only V2 (IRG7-H) says whether the carried BASE lines changed any cluster.

The C falsifiers C-F1 to C-F3 of the earlier versions are withdrawn with the arm (§3.4).

## 12. Not in this test (and why)

- **Arms A and B on held-out** (user decision 1; §3.2, §3.3). With them go:
  - R7's nulls and the complete-universe referee for A;
  - TS0 vs TS10, UC restricted to the tie set, T1 units, and the individual's copy number as a tie-breaker (first
    review L3).
- **Arm C on held-out** (user decision 2, 2026-09-27: descriptive only, to be revisited later). With it go:
  - O1C_L, NR, the LODN comparators and the NEW truth on held-out;
  - no-read copies on another contig;
  - adding no-read copies to O2's copy table;
  - chimpanzee and orangutan RepeatMasker.

  §3.4 lists what a later C test must fix.
- **RG + R7** (A is dropped). **RG3 on the R3 GTF**: R3 is not the default, so this stays dev-descriptive (§4.1).
- **Readthrough R3** stays EFFECTIVE and not flipped. BASE stays the input.
- **A merge rule** for the BASE `gene_id` pairs that already share an exonic base (1,818 / 588 / 143 on dev): a
  different arm, not proposed.
- **The annotation-only opportunity count** for RG (RefSeq readthrough records whose parents include a multi-member
  family outside chr16/chr20). It reads held-out truth tables, so it is the user's call and is not run.
- **A byte-identical incremental PAF.** It is impossible with any useful speedup (`rep4_incr` §2; draft r1126).
- **Dropped at dev (RG):** the exact-junction adjacency `rg.py` (fails S 3/3, r1125); donor-or-acceptor (fails S on
  human dev); exon-or-junction (identical to RG3); the `_rg`, `_rgu` and `.<k>` namings.
- **The Rust port of RG** (`--gtf-regroup`, `rt4_ghost.md`) happens only after RG's verdict (§6).

## 13. Order, lock time, stop rules, machine rules

1. **Accept.** The user accepts this file, arm C's move to descriptive, the third reuse and the lock time of §13.1.
   Amendment 1 records that acceptance and this file's sha1.
2. **Dev implementation and gates** (§6, dev only: IRG1-IRG7). Write Amendment 1: the SHA1_* values, the families
   method (F_INCR or F_FULL), the command lines and the dev tables. **No arm command runs on held-out data before
   Amendment 1 exists.**
3. **BASE-only held-out prerequisites.** They produce no arm number.
   - **Queue 1, running since before acceptance** (`rep_prereq_queue.sh`):
     - Liftoff self-lifts (gorilla, orangutan, human) and their read-support tables (five samples);
     - BASE families for human_A119b, gorilla_OR6737, gorilla_KB3781 and orangutan_PPY (`fam_call.sh`,
       `rt_bin_frozen`).
     - Calls are bounded, and exit 75 means run the same call again.
   - **Queue 2, written and not started** (`rep_prereq2_queue.sh` 638f48fb, helper `rep_prereq2_genes.py` 52d3f70b).
     It starts after queue 1 logs `REP PREREQ COMPLETE`, on the user's go, and refuses to start while queue 1 runs.
     - IRG8 on gorilla and human, then the chimpanzee and orangutan genes tables (`families_gw/species/{chimpanzee,
       orangutan}/`).
     - The genome-wide referee for gorilla, human, chimpanzee and orangutan, resuming the pilot prefixes (§7.4), from
       the frozen copies (§6), `--budget-s 540`, each call under `timeout 600`.
   - **After Amendment 1:** freeze checks (a)-(d) of §6. The restricted BASE rows, the queue-2 products' sha1s and any
     undecided substrate go into Amendment 2.
4. **Arm RG.**
   - **a. V2 first (IRG7-H).**
     - `rg3.py` and `rg3_null.py` on V2's BASE GTF.
     - RG3 and NULL_RG3 families by both F_INCR (`fam_call_incr.sh`) and F_FULL (`fam_call.sh`), then `cmp` of the
       `clusters.tsv`. The outcome fixes the method for the other five substrates (§6).
     - Skipped if Amendment 1 already set F_FULL.
   - **b. V1 and V3-V6.**
     - `rg3.py` and `rg3_null.py` on the BASE GTF. Input sha1s: A119b 631c9f11, testis 50f239d6, OR6737 3a8c6410,
       KB3781 f1022f4b, PTR dead746f, PPY b013a74a. The output sha1s are recorded.
     - Then `fam_call_incr.sh <s> RG` and `<s> NULLRG`, or `fam_call.sh` under F_FULL.
   - **c.** `readthrough_eval.py score --heldout --families`, with arms BASE, RG and NULLRG and the §7.2
     `--drop-contigs`. NULLRG is a plain arm, never `--null`, which would suppress its `e.*` rows.
   - **d.** SHA1_CLS (RG3 and NULL_RG3; the old `rg.py` rewrite, reported), SHA1_RG_EVAL and SHA1_REP_PAIRS; emu for the
     dev-contig-deleted rows and the relabel null.
5. **Verdict.** Compute RG's verdict (§9), append the Outcome, and give every number a new register row.

### 13.1 Lock time (review3 H1; what the user accepts)

| step | lock time | basis |
|---|---|---|
| queue 1 remainder (BASE families of PPY and A119b) | about 2-3 h left; already running, not part of this acceptance | queue timing (§0) |
| queue 2: IRG8 (human, gorilla) and the chimpanzee and orangutan tables | ≈ 0.2-0.5 h | unmeasured: four one-pass reads of 0.66-1.7 GB GFFs |
| queue 2: the referee, 86 remaining shards over 4 species | ≈ 3-4 h | pilot: 115-146 s per 1,000-protein shard (2.7-3.5 h), plus a 18-33 s GFF pass per call; edges and MCL unmeasured |
| dev gates IRG1-IRG7 | ≈ 0.3-0.5 h | dev families runs take 43-57 s (full) and 7-15 s (incremental) |
| freeze checks (a)-(d) on six BASE | ≈ 0.3-1 h | emu R0 takes 1.8 s on chr16; genome-wide unmeasured; 11 emu runs per sample |
| RG3 and NULL_RG3 GTF rewrites | < 0.1 h | 0.2-0.5 s per contig |
| **RG families, 2 arms × 6 samples: F_INCR, both arms incremental** | **≈ 2.2-5.3 h** | `rep4_incr` §5: 0.12-0.27 of the full run plus 18-30 min fixed per arm; full BASE runs 21-105 min per sample |
| … F_INCR, every NULL_RG3 falls back | ≈ 7.9-10.7 h | the same |
| … F_FULL (IRG7 or IRG7-H failed) | ≈ 13.6-15.9 h | the BASE family times (407-477 min per arm) |
| IRG7-H (V2: F_FULL for both arms, besides F_INCR) | ≈ 0.7 h | 2 × 21 min (F_INCR only) |
| scoring (`readthrough_eval` on 6 samples × 3 arms, SHA1_CLS, SHA1_RG_EVAL, SHA1_REP_PAIRS, emu rows) | ≈ 0.5-1.5 h | unmeasured genome-wide |
| **total: F_INCR, both arms incremental (expected)** | **≈ 7-14 h** | |
| total: F_INCR, every NULL_RG3 falls back | ≈ 13-19 h | |
| **total: F_FULL (the most the user accepts)** | **≈ 18-24 h** | includes up to 0.3 h of discarded V2 F_INCR runs |

For comparison, review3 costed the C+RG version at about 36-48 h, of which C was about 18-25 h.

**Stop rules.**
- **No change after Amendment 1** to a rule, a NULL, a clause, a floor or a tolerance. A bug fix re-runs every run of
  the arm and is recorded as an amendment.
- **A failed dev gate stops the arm, except IRG7 and IRG8.** An IRG7 failure switches the method to F_FULL. An IRG8
  failure leaves V5 and V6 undecided (§9). A failed freeze check or identity gate acts only as §6 and §8.1 state; it
  never touches a rule.
- **A step that cannot be computed within the machine rules makes its clause "not measured"**, which caps the verdict
  (§9).
- **No variant is ever substituted for a §2 rule.** Examples: `rg.py` or donor-or-acceptor for RG3; F_INCR on a
  substrate after IRG7-H has failed.

**Machine rules.**
- **One heavy process at a time, in the foreground:** `flock -w 900 /mnt/linuxdisk/tmp/rustle_heavy.lock timeout 600
  <cmd>`. The queues use `queue_lib.sh`'s serial flock with bounded calls. Genome-wide steps run as bounded, resumable
  calls.
- **Never `pkill -f`**; kill by PID. Binaries run only from their frozen copies.
- **Outputs and `TMPDIR` go under `/mnt/linuxdisk`.** Arm outputs: `/mnt/linuxdisk/tmp/rustle_figures/rt_arms/<sample>/`
  (next to BASE, as `fam_call.sh` writes them), plus `/mnt/linuxdisk/tmp/rustle_figures/rep_arms/<sample>/` for scores.
  Dev work: `/mnt/linuxdisk/tmp/rustle_figures_dev/rep*`.
- **Until Amendment 1, development contigs only.**

## 14. The hostile reviews, item by item

### 14.1 Third review (`rep3_review.md`, sha1 d86bfb7d; verdict NOT READY)

| item | arm | the defect | fix here, or where it went |
|---|---|---|---|
| **H1** | shared | the lock time is understated about 2×; the referee and the chimpanzee/orangutan genes tables are not queued; §0 misstates the referee | **Fixed.** C (18-25 h of the 36-48 h) left held-out with its verdict (user decision 2). §13.1 costs every remaining step: RG's families at the instrument's measured dev cost, with the F_FULL line if it fails. The header and §13 step 1 ask acceptance of that total. Queue 2 (`rep_prereq2_queue.sh`, written, not started) builds the referee for 4 species and the chimpanzee/orangutan tables. §0 now says "partial pilot, 1 of 21-23 shards per species, no families" |
| **H2** | C | C1/C2 floors count families that cast no vote | **Moot: C is descriptive.** Recorded in §3.4 as a requirement for any later C test (dev decisive counts C1 3, C2 2) |
| **M1** | C | the judged truth contains LIFT01, which shares O1C_L's search | Moot; recorded in §3.4 (judge under NEW′, which equals NEW on dev) |
| **M2** | RG | FU is implied by S; the "recall-type" rationale is false | **Fixed.** §8.1 rewritten; cleared fraction reported; §8.3 FU failing outcome ("S not judged and no span-disjoint separation") |
| **M3** | C | seeds are only families of ≥ 2 expressed loci | Moot; stated as C's scope in §3.4 (93 of 1,142 NC_073244.2 loci are family members) |
| **M4** | C | the 600 s cap excludes the largest contigs | Moot; recorded in §3.4 |
| **M5** | C | no identity-stratified row | Moot; recorded in §3.4 |
| **M6** | RG | IRG4's complete-universe rows are new numbers | **Fixed** (§0 item 3: recorded as found, cannot change anything) |
| **M7** | RG | P1 drops pairs whose loci resolve to no universe gene | **Fixed** (P1_all reported beside P1, §8.1; dev effect stated in §7.4) |
| L1 | RG | T_g undefined | **Fixed** (§5.1) |
| L2 | RG | post-exposure marking inconsistent | **Fixed** (§0 and §10 both mark predictions 1, 3, 4, 5 and 6; prediction 2 is not) |
| L3 | C | "Human dev had no unread truth copies" is false | Fixed in §3.4 (chr16 has 1 unread class-(i) NEW copy; NEW confirms no human flag) |
| L4 | RG | "undecided" omits the chimpanzee/orangutan genes tables | **Fixed** (§9) |
| L5 | RG | the dev-contig restriction has no gate; `--drop-contigs` not named | **Fixed** (IRG6; §7.2 and §13 step 4c name the flags) |
| L6 | shared | falsifiers have no floor | **Fixed** for RG-F1 (n reported; below 30 gained pairs it is underpowered). C-F1/C-F2 moot, recorded in §3.4 |
| L7 | C | freeze check (d) is not operational | Moot: the BAM check served C only and is removed; its operational form is recorded in §3.4 |
| L8 | C | O1C_L's side choice depends on genomic order | Moot; recorded in §3.4 (0 flips on dev) |
| L9 | C | run C on gorilla first | Moot; recorded in §3.4 |
| L10 | RG | EFFECTIVE above the independence product | **Fixed** (§10 item 7: 0.25, correlation stated) |
| L11 | RG | the genes-table builder is unfrozen | **Fixed** (IRG8; sha1 pins in the helper and in queue 2) |

### 14.2 Second review (`rep2_review.md`)

| item | arm | the defect | status in this file |
|---|---|---|---|
| **X1** | B | TS0 judged at placement level; at copy level it fails B1 on dev | Dropped with arm B. The copy-level numbers are in §3.3 |
| **X2** | RG | FU near-definitional; the rule mostly fragments genes | **Fixed.** Rule RG3 (exon overlap; disclosed dev selection, §0, §2.1); judged clause S required for EFFECTIVE and counted in F_RG; "0 new fused pairs" moved to the identity gates (§8.1); FU's remaining role stated honestly (review3 M2) |
| **X3** | B (UC) | B4 unmeasurable; the real failure only reported | Dropped with arm B. The 244-of-298 finding is in §3.3 and r1122 |
| **H1** | B | B5's third condition inverted | Dropped with arm B. The corrected reading is in §3.3 |
| **H2** | A | A1's null floor asymmetric | Dropped with arm A |
| **H3** | C | C4 pooled, so the array decides | Moot: C is descriptive. The per-family β test was run on dev and failed (405 of 1,000; r1123) |
| **H4** | C | EM misdescribed; C3's EM half automatic | Moot: C is descriptive |
| **M1** | shared | §0 missed the published `tables_v3` BASE rows | **Fixed** (§0, with scope per sample); post-exposure marking made consistent (review3 L2) |
| **M2** | shared | verdict rules contradict each other | **Fixed.** §9's common rules; F_RG; undecided includes the genes tables |
| **M3** | C | class (ii) is O2, not "no reads" | Moot: C is descriptive (§3.4 keeps the class partition) |
| **M4** | C | C-F2 cannot fire | Moot: C is descriptive |
| **M5** | A | A-F2 fires on ties | Dropped with arm A |
| **M6** | shared | could dev tables still change the rules? | **Fixed.** §0 item 3 and §13 step 1, plus the IRG4 exception (review3 M6) |
| **M7** | B | B's truth never covers reads from uncatalogued copies | Dropped with arm B; stated in §3.3 |
| **L1** | RG | piece names collide with transcript ids | **Fixed.** `.rg<k>` names; assert 3 also checks transcript ids (§2.1) |
| **L2** | RG | the fused-pair floor changed units | **Fixed.** Both unit changes stated in §2.3 |
| **L3** | C | 929 / 1,000 is a PRIM figure | Moot |
| **L4** | C | C judged on one genome | Moot |
| **L5** | shared | §2.E's Liftoff row read as universal | **Fixed.** §2.3 has only RG's Liftoff pair truth |
| **L6** | B | `setans.py` docstring | Dropped with arm B |
| **L7** | C | LIFT01 uses the rule's own search | Moot; carried into §3.4 (review3 M1) |
| clause table | shared | several clauses could not fail | **Fixed.** §8.3 names a failing outcome for every judged clause; RG3 fails S on the R3 chr20 GTF on dev; FU and T are described as what they are |

### 14.3 First review (`rep_critique.md`): status in this file

| item | status |
|---|---|
| K1, K1b (B's circular truth; NM filter) | dropped with arm B |
| K2 (C's circular truth, handicapped comparator) | moot: C is descriptive. NEW never used clause E, and the literal Liftoff was the comparator (§3.4) |
| H1 (A's safety clauses, referee dropped) | dropped with arm A; the complete-universe referee is kept for RG (P1, P2) |
| H2, H3 (C's primary-based flag; one array carries C2) | moot: C is descriptive |
| H4 (β never tested) | tested on dev and failed (r1123) |
| M1 (free constants) | holds: every RG clause is an exact inequality; S is a strict sign count |
| M2 (NULL_A too easy) | dropped with arm A |
| M3 (rule text vs instrument) | holds: §2.1 is written from `rg3.py`, and §2.2 from `incr_mm2.py` |
| M4 (exposure, bridges) | holds: V1/V3 without the Fig. 7 contigs; dev-contig-deleted rows; §7.3 |
| L1-L3 (B's cap, quantification, assumptions) | dropped with arm B |
| L4 (scope wording) | C moot; RG is seeded from assembled loci and does no genome-only discovery |
| L5 (A pre-announced as a gain test) | dropped with arm A |
| L6 (no O1/O2 conflation) | holds: RG changes only nodes |

## Amendments

(none yet)

**Amendment 1 must contain, before any held-out arm command:**
- the user's explicit acceptance of this file, of arm C's move to descriptive (user decision 2), of the third reuse and
  of the lock time of §13.1;
- this file's sha1;
- every SHA1_* value: SHA1_CLS, SHA1_RG_EVAL, SHA1_REP_PAIRS, SHA1_INCR and SHA1_FAM_INCR, or "F_FULL" with the failing
  IRG7 item;
- the command lines;
- the dev tables of §6: IRG2, IRG4 (complete-universe rows and P1_all, recorded as found) and IRG7.

**Amendment 2 (after the freeze checks, before any arm command) must contain:** the restricted BASE rows; the sha1s of
the four referee `ph.families.tsv` and of the chimpanzee and orangutan genes tables; BASE's `mid_occ` per substrate;
any undecided substrate, with the user's decision on it.
