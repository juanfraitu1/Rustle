# Pre-registration: held-out completeness of `--assemble-only` — `--polish-tes pas-end` judged; tools, `--polish-subchain drop` and `--polish-tss rescue` descriptive

**Written 2026-09-27, before any held-out run of `--polish-tes`, `--polish-subchain` or `--polish-tss`, and before any
SIRV run of our assembler.** This is the binding text. A default flip is the user's call whatever the outcome.

**User goal** (`/goal`): *"ensure the assembler part of the pipeline tries as best as possible to emit complete
transcripts, meaning that when running gffcompare afterwards there are few partial transcripts like m categories"*.
The goal's Stop hook asks for held-out verification of what shipped on dev (register r1130, r1131, r1133).

**What is tested.**
- **Judged:** `--polish-tes pas-end` (r1133) against BASE and against a matched NULL (§3).
- **Descriptive only:** BASE against StringTie / FLAIR / IsoSeq; `--polish-subchain drop` (r1130); `--polish-tss
  rescue` (r1131). Their held-out rows carry no verdict (§5, exposure).

---

## 0. What was seen before this file

- **Dev (design evidence; the contigs every option was selected on):** human A119b chr20 (and chr16) and gorilla OR6737
  NC_073244.2 — `tss2_impl.md`, `ct_levers.md`, `ct3_review.md` (session `scratchpad/figs/`), `bench/ASSEMBLY_POLISH.md`
  addenda 1-2, and this author's dev runs of the frozen binary and scripts (§2; scratch
  `/mnt/linuxdisk/tmp/rustle_figures_dev/complete_ho_dev/`).
- **SIRV (truth side and tools only):** the `ct_sirv` README facts (69 isoforms, 61 multi-exon, 9 true ISMs, 93.3%
  coordinate duplicates, SIRV107 strand artefact) and the three tools' SIRV rows, re-scored by the frozen `ho_eval.py`
  on the frozen tool GTFs: multi-exon 44 / 64 / 109, `=` 39 / 49 / 62, truth `=` isoforms 39 / 49 / 47 of 61 (StringTie
  / FLAIR / IsoSeq; identical to `ct_sirv`). **No SIRV output of our assembler exists.** The 380-read smoke output stays
  sealed in `ct_frozen/sealed/ct_sirv_smoke/`; this author only ran `stat`/`find -printf` on it (modes, ctimes).
- **Held-out:** this author read no class row, chain count or TES/TSS number of any held-out substrate. Read: file
  sizes, sha1s, the 09-25 BASE wall times (`runs/<s>/assemble.time`) and the rt3 arm call durations (≤ 452 s).
  **Seen by the user and earlier sessions:** BASE class shares of human_A119b (figures 1-2, chr20-22 and
  `eval_annotated`) and of gorilla_OR6737 (genome-wide, figures 1-2), and the BASE rows of the readthrough preregs.
- **Exposure.** V1-V6 carry the R/RQ1 (r1117) and R3 (r1119) verdicts. This is their **fourth** reuse (the count of
  `PREREG_complete_transcripts_2026-09-27.md` §0, whose own reuse was accepted by the user and then **not spent**:
  withdrawn on dev). See gate G0.

## 1. Arms (all genome-wide, the driver's `assemble` stage, frozen binaries `tss3_bin_frozen/`)

| arm | how | role |
|---|---|---|
| **BASE** | `rt_arms/<s>/<s>.BASE.gtf` → `runs/<s>/<s>.gtf` (09-25, `copy_assign` 753a3b4d; sha1s in `ct_frozen/INPUTS_SHA1.tsv`); SIRV: the driver, options unset | comparator (reused; G5 + IV_pre) |
| **PASEND** | `RUSTLE_POLISH_TES=pas-end` + driver | **judged** |
| **NULL** | `tes_null.py` on BASE + PASEND + BAM + FASTA (§3) | **judged comparator** of J2 |
| DROP | `RUSTLE_POLISH_SUBCHAIN=drop` + driver | descriptive |
| RESCUE | `RUSTLE_POLISH_TSS=rescue` + driver | descriptive |
| stringtie / flair / isoseq | the lab's de novo GTFs (`assembly/human/eval_annotated/`, `assembly/gorilla/eval_all/`; SIRV `ct_sirv/tools/eval/`), sha1s in `frozen/INPUTS.tsv` | descriptive; human_A119b, gorilla_OR6737 and SIRV only |
| SIRV kcd | BASE / PASEND / NULL with `--keep-coordinate-duplicates` (the driver's command expanded; IV1-SIRV checks the expansion) | S clauses in both count modes |

Every option other than the one named is unset (`env -u RUSTLE_POLISH_* -u RUSTLE_READTHROUGH_JUNCTIONS*`). The exact
commands are `frozen/run_ho.sh` (§9).

## 2. Dev design evidence (in-sample; frozen binary and scripts; chr20 / NC_073244.2; never in the verdict)

- **Reproduction.** tss3 (`copy_assign` b1709a96) dev BASE / PASEND / DROP GTFs are byte-identical to the tss2
  products; dev BASE sha1s equal `ct_frozen` DEV_BASE (6fa4dbaa / 41815166). `tes_null.py` re-derives **106/106 and
  118/118** of the binary's moves from the BAM (PARITY) and passes IV_pre; the scorer reproduces tss2_impl's labelled
  numbers exactly (moved 106 / 118, labelled 73 / 98, near 6→32 / 8→66, genes 28 recovered 2 lost / 43 and 5).
- **Why the judged metric is the moved set (decided here, on dev, before any held-out number).** The substrate-wide
  count of reference genes with a TES-within-50-bp query is **saturated**: BASE → PASEND 422 → 426 (chr20) and
  750 → 750 (gorilla), NULL 426 / 750. The genes whose moved transcript gains an annotated TES already have another
  transcript there. On the moved set M the same label moves 2 → 28 and 6 → 44 genes. A substrate-wide clause would
  have refuted on dev a rule whose every moved end is measured to improve; it is kept as a reported row.

| dev | \|M\| (labelled) | moved-set TES50 genes B / P / N | TES50 transcripts of M B / P / N | \|F\| (free) | substrate-wide TES50 genes B / P / N |
|---|---|---|---|---|---|
| chr20 | 106 (73) | 2 / **28** / 20 | 6 / 32 / 21 | 33 | 422 / 426 / 426 |
| NC_073244.2 | 118 (98) | 6 / **44** / 41 | 8 / 66 / 62 | 23 | 750 / 750 / 750 |

- **Descriptive arms on dev.** DROP: `c` 5.7 → 2.8% / 5.4 → 3.1%, chains −2 / −3, genes losing every `=`/`c` query 7 /
  7. RESCUE: chr20 +47 transcripts, +3 chains, +1 TSS250 gene; gorilla has no cap signal, so RESCUE = tag (0 added).
- **Tools on dev** (multi-exon shares `=` / `c` / `k` / `m` / `n` / `j` / other; chains):

| chr20 | = | c | k | m | n | j | other | chains | NC_073244.2 | = | c | m | chains |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BASE | .210 | .057 | .033 | .023 | .033 | .528 | .116 | 1,059 | BASE | .409 | .054 | .017 | 1,595 |
| StringTie | .168 | .008 | .055 | .069 | .068 | .412 | .221 | 861 | StringTie | .370 | .027 | .051 | 1,374 |
| FLAIR | .077 | .012 | .034 | .041 | .056 | .558 | .222 | 1,026 | FLAIR | .238 | .019 | .039 | 1,393 |
| IsoSeq | .048 | .031 | .022 | .042 | .072 | .540 | .244 | 1,253 | IsoSeq | .142 | .043 | .057 | 1,655 |

## 3. The NULL (binding; `tes_null.py`, verbatim in its docstring)

- **M** = the multi-exon transcripts PASEND moved (`tes_end_moved_from`).
- **Evidence**, re-derived per t ∈ M exactly as the binary builds it: primary spliced records (not flag 2820),
  deduplicated on (FLAG 0x10 strand, start, end, intron chain); t's own 3′ ends = records with t's strand and genomic
  intron chain, inside t's BASE last exon; single-linkage clusters (gap ≤ 21), ≥ 2 ends, mode = most ends (ties
  3′-most). `pas(c)` = AATAAA/ATTAAA wholly in oriented [c−35, c−10]; `primed(c)` = ≥ 60% A in [c+1, c+20] or an A6
  run. The rule's target is the most-3′ cluster with `pas ∧ ¬primed`.
- **NULL target** = the most-3′ cluster mode c < lb (t's emitted end) with **¬primed(c) ∧ ¬pas(c)**: the same move to an
  unprimed own-read cluster, with the PAS conjunct negated. Decision `free` (end moved there) or, when no such cluster
  exists, `tie` (NULL keeps PASEND's move, so t contributes equally to both arms). **F** = the `free` transcripts.
- **NULL.gtf** = PASEND.gtf byte for byte, except the free transcripts' transcript line and 3′-most exon line, and a
  `tes_null "free|tie"` attribute. `cov` is not recomputed (gffcompare ignores it).
- **Built-in stops (exit 3):** IV_pre — PASEND with `tes_*` stripped, `cov` masked and each move undone must equal BASE
  line for line; PARITY — for every t ∈ M, `primed(lb)`, #proven = `tes_clusters`, and the most-3′ proven mode = PASEND's
  new end. Both run on every substrate, so they also verify the never-rerun `--genome-wide` TES path (tss2_impl caveat)
  and the equivalence of the 09-25 BASE with tss3 beyond G5's one sample.
- Deterministic (a re-run is byte-identical; checked on dev).

## 4. Gates (`run_ho.sh gates`; recorded in Amendment 1 before step 2; any failure stops the work)

- **G0 reuse.** The main session records, in Amendment 1, that the user accepts spending the fourth reuse of V1-V6
  (accepted 09-27 for this goal's held-out set in the complete-transcripts prereg, unspent) on pas-end. Without it,
  J1/J2 are reported as descriptive and the verdict is "not tested".
- **G1 frozen code.** `sha1sum -c` passes in `tss3_bin_frozen/` (4/4), `complete_ho/frozen/` and `ct_frozen/` (7/7, the
  imported scorer core and REF_SHA1.tsv); driver `tools/rustle_pipeline.sh` = 2ec71247; `git diff --quiet HEAD --
  figures/ tools/rustle_pipeline.sh`.
- **G2 BASE sha1s** of the six = `ct_frozen/INPUTS_SHA1.tsv` (re-checked per sample before its arms).
- **G3 inputs** (`frozen/INPUTS.tsv`): tool and SIRV sha1s, gffcompare 0.12.10 (32cbdd64), `figures/inputs.local.tsv`
  and `samples.tsv` sha1s; size + mtime of every BAM and FASTA.
- **G4 SIRV seal unread:** `ct_sirv_smoke` is mode 500 / ctime 1790533721 with 12 files of mode 000, ctime unchanged.
- **G5 BASE re-run byte identity:** human_testis, frozen tss3, options unset → GTF `cmp`-identical to
  `runs/human_testis/human_testis.gtf` (params/quant reported).
- **G6** `tes_null.py selftest` and `ho_eval.py selftest` pass, and the first N lines of this file (the frozen text,
  everything above `## Amendments`) hash to `frozen/PREREG.sha1` (sha1 and N); amendments are appended below it.
- **Per substrate:** IV_pre and PARITY (§3); IVT and IVN (§6). A failure is fixed in the code (amendment; every
  PASEND/NULL arm re-run), never in a clause. If IV_pre fails because the 09-25 BASE differs from tss3 unset, BASE of
  that substrate is rebuilt with tss3 unset (amendment) and used for all its arms.

## 5. Substrates (species never pooled; dev contigs excluded everywhere)

| id | scorer key (restricted reference sha1, `ct_frozen/REF_SHA1.tsv`) | tools | prior exposure |
|---|---|---|---|
| **S0** | SIRV E0 of human_testis, `--sirv` (81051062); dedup **and** kcd | yes | fresh for our output (first use) |
| V1 | human_A119b `annotated_minus_chr16_chr20_chr21_chr22` (c89dfe1a) | yes | R/RQ1, R3; BASE classes in figs 1-2 |
| V2 | human_testis `annotated` (8ed998a1) | – | R/RQ1, R3 |
| V3 | gorilla_OR6737 `annotated_minus_NC_073244.2` (9271749f) | yes | R/RQ1, R3; r1074-r1076; figs 1-2 |
| V4 | gorilla_KB3781 `annotated` (a7e0e06e) | – | R/RQ1, R3 |
| V5 | chimp_PTR `annotated` (f4b58b50) | – | R/RQ1, R3 |
| V6 | orangutan_PPY `annotated` (35a55052) | – | R/RQ1, R3 |

- chr21/chr22 stay out of V1 (dev set of the 09-23 polish search); they are not scored here.
- S0 is disjoint from V2: the SIRV reads are V2's unmapped reads, re-aligned to SIRV1-7 only (`ct_sirv` README).
- **Circularity.** V3-V6 labels are Gnomon models (KB3781 is the mGorGor1 individual); whether their TESs followed
  the same reads' internally primed ends is unknown. EFFECTIVE therefore also needs J1 and J2 on V1 and V2 (§7).
- **BASE vs tools is descriptive** everywhere: BASE's held-out class shares were seen (§0).

## 6. Metrics and clauses (`ho_eval.py`; gffcompare 0.12.10; integers only)

- **Class shares** (tmap multi-exon queries): `=`, `c`, `k`, `m`, `n`, `j`, other; partial = c+k+m+n; matched chains
  (`.stats`); chain precision = `=` multi-exon / multi-exon; gffcompare intron-chain Sn/Pr.
- **TES50 / TSS250 label:** a query q and a multi-exon reference transcript r of gene g on the same contig and strand
  share the **last** (first) intron and |TES(q) − TES(r)| ≤ 50 bp (|TSS| ≤ 250 bp).
  - substrate-wide: `t.tes50_genes`, `t.tss250_genes`, and the transcript counts (reported);
  - **moved set** (judged): for X ∈ {BASE, PASEND, NULL}, `genes_X` = genes g hit by X's own 3′ end of some t ∈ M.
- **Descriptive per arm vs BASE:** `=`-pair losses/gains, chains Δ, added/dropped, **genes losing every multi-exon `=`/`c`
  query** (`v.genes_eqc_lost`), TES50/TSS250 genes recovered/lost, 3′/5′ moves (toward/away from the nearest reference
  TES of the same last intron).

| clause | passes iff | judged on |
|---|---|---|
| **J1** | genes_P > genes_B (NA if \|M\| < 10) | V1-V6 |
| **J2** | genes_P > genes_N (NA if \|F\| < 10) | V1-V6 |
| **J3** | PASEND's set of (query, ref) `=` pairs = BASE's, and matched chains equal | V1-V6, S0 dedup + kcd |
| **J4** | `=` multi-exon count and multi-exon count equal to BASE's (precision identical) | V1-V6, S0 dedup + kcd |
| IVT (gate) | PASEND vs BASE: 0 5′ moves, 0 added, 0 dropped, TSS250 gene set equal | all |
| IVN (gate) | NULL vs BASE: J3 and J4 hold, 0 5′ moves, 0 added/dropped | all |
| **S_EQ** | truth isoforms with a multi-exon `=` query: none lost vs BASE, and J3 | S0 dedup + kcd |
| S_ACTIVE | PASEND moves ≥ 1 SIRV end (else "inert") | S0 |
| **S_TES** | moved ends labelled by a truth isoform of the same last intron: away ≤ toward | S0 dedup + kcd |

NA and a missing substrate are not passes. J3/J4 are true by construction; a failure means the construction is broken.

## 7. Verdict (`ho_eval.py verdict --out complete_ho/tables [--final]`)

| verdict | condition |
|---|---|
| **EFFECTIVE (default candidate)** | gates pass; J3 and J4 everywhere; S not failed; **J1 ≥ 5/6 and J2 ≥ 5/6, with J1 and J2 on V1 and V2** |
| **KEEP OPT-IN** | not EFFECTIVE and not REFUTE (e.g. J1 holds but the PAS does not beat the matched NULL) |
| **REFUTE** | J3 or J4 fails anywhere; or S_EQ fails; or S_TES fails while active; or J1 passes on ≤ 4/6 |

- SIRV inert (S_ACTIVE false in both modes) is expected (r1132: a canonical PAS sits at 1 of 69 truth TESs). S0 can
  therefore only falsify pas-end, never support it; EFFECTIVE does not need S_ACTIVE.
- `--final` counts a substrate still missing after the pre-registered retries as J1/J2 not passed.
- Reported beside it: every §6 row per substrate; \|M\|, \|F\|, `null_var_pas` counts; the substrate-wide rows.

## 8. Predictions (before any number)

1. **BASE vs tools (V1, V3):** BASE has the lowest `m`, `n` and "other" shares and the highest `c` and `=` shares of the
   four arms; matched chains IsoSeq > BASE > FLAIR, StringTie; substrate-wide TES50 genes IsoSeq > BASE > FLAIR >
   StringTie. **S0:** BASE `m` ≤ every tool's (tools 0); BASE `c` below IsoSeq's .349; BASE dedup truth `=` isoforms
   ≥ 39 (StringTie), kcd ≥ dedup.
2. **pas-end:** J1 6/6; J3, J4, IVT, IVN, IV_pre, PARITY pass everywhere; \|M\| 1-4% of multi-exon transcripts;
   \|F\|/\|M\| .15-.40.
3. **J2** passes on V1 and V2 with the largest relative margins; on V3-V6 it passes with margins ≤ 10% of genes_P
   (dev gorilla 44 vs 41). Predicted verdict: **EFFECTIVE**; the named risk is J2 on the Gnomon substrates (→ KEEP
   OPT-IN).
4. **S0:** inert in both modes; S_EQ passes.
5. **DROP:** the `c` share falls 40-55% (relative) on every substrate; chains −0.1 to −0.3%; genes losing every `=`/`c`
   query 0.5-1.5% of BASE's; on S0 it loses at most one truth `=` isoform (SIRV303 is the only structurally exposed
   true ISM).
6. **RESCUE:** both gorilla libraries have no cap signal ⇒ RESCUE = BASE classes; on V1/V2 +0.5-1.5% multi-exon
   transcripts, chains +0.1-0.5%, chain precision down. No prediction for V5, V6 and S0 (cap signal unknown).

## 9. Order, stop rules, machine rules

`F=/mnt/linuxdisk/tmp/rustle_figures/complete_ho/frozen`; outputs under `complete_ho/` (tables in `complete_ho/tables`).
0. Freeze: done (Freeze record). G0 in Amendment 1.
1. `$F/run_ho.sh gates` → Amendment 1 copies `complete_ho/gates.txt`.
2. `$F/run_ho.sh sirv` → S0 rows recorded (Amendment 2) before step 3; nothing changes after them.
3. `$F/run_ho.sh sample <s>`: human_testis, chimp_PTR, orangutan_PPY, gorilla_KB3781, gorilla_OR6737, human_A119b.
4. `$F/run_ho.sh verdict` (then `--final` if a substrate cannot be built) → Outcome; register rows.

- **Heavy calls:** `flock -w 900 /mnt/linuxdisk/tmp/rustle_heavy.lock timeout 600`, foreground, serial, `TMPDIR`
  under `/mnt/linuxdisk`. An assembly that exits 124 gets **one** retry under `timeout 1200` (logged; the A119b
  arms took 405-452 s on 09-25/26); a second failure leaves that substrate missing. The scorer is resumable (exit 75).
  Never `pkill -f`.
- **Stops:** any sha1 mismatch; tes_null exit 3 (IV_pre / PARITY); IVT or IVN failing; an unexplained difference in
  IV1-SIRV. Nothing is re-tuned: not the rule, the NULL, a tolerance, a clause or a floor. The scorer refuses every
  non-dev substrate (S0 included) unless `--prereg` names a file containing its own sha1 (the Freeze record).
- **Disk (~75 GB free):** keep GTFs, `NULL.tsv`, `.stats/.tmap/.refmap`, tables, logs; delete `work/*/arms/*.pkl`
  after the Outcome. No BAM is written.

## 10. Not in this test

O1 reach (pas-end changes families, r1133) · `--polish-tss split` · `--polish-subchain tag` · any `m`-class lever
(closed: r1075 adopted, r1080 refuted) · chr21/chr22 · SQANTI3.

## 11. Hostile review of this file (by its author, before the freeze record) and fixes

| # | finding | fix (applied above) |
|---|---|---|
| H1 | The obvious clause (substrate-wide TES50 genes) is saturated: on dev it refutes pas-end (gorilla 750 = 750) although every moved end is measured to improve. | Judged metric = moved-set genes (§2, §6); substrate-wide kept as a row. The choice is dev-informed and disclosed; it is the gene count tss2_impl already reported, not a new statistic. |
| H2 | J1 is a low bar: M is selected on primed ends, which rarely sit at annotated TESs (dev 2 and 6 genes). | J1 alone never gives EFFECTIVE; J2 (PAS vs a PAS-free unprimed move of the same transcripts) carries the claim. |
| H3 | The NULL's `tie` fallback makes J2 depend on \|F\| only (dev 33 / 23). | \|F\| ≥ 10 floor (NA is not a pass); \|F\|/\|M\| predicted and reported. Ties are shared, so they cannot favour either arm. |
| H4 | The NULL's "PAS-free" is canonical-only; a variant hexamer can still sit at its target. | Kept (the literal complement of the rule's conjunct; it can only strengthen the NULL). `null_var_pas` is reported. |
| H5 | Gnomon labels on V3-V6 may be circular with the same reads. | EFFECTIVE also requires J1 and J2 on V1 and V2 (§7). |
| H6 | BASE comes from another binary (753a3b4d), and the `--genome-wide` TES path was never run. | G5 on human_testis plus IV_pre and PARITY on every substrate (§3). |
| H7 | S0 cannot support a PAS rule (r1132), so it could pass vacuously. | S0 is stated as falsification-only (S_EQ, S_TES while active); inertness is predicted, not counted. |
| H8 | The user's fourth-reuse acceptance was given for arm A, not for pas-end. | G0: recorded explicitly before step 2, else "not tested". |
| H9 | A119b arms are near the 600-s bound; a missing substrate could silently vanish from a 5/6 count. | One logged retry at 1200 s; `--final` counts a missing substrate as not passed. |
| H10 | "Precision identical" could be read on the rounded `.stats`. | J4 uses the integer `=` and multi-exon counts of the tmap. |
| H11 | The SIRV `kcd` arm is not a driver mode. | The expanded command must reproduce the driver's dedup GTF byte for byte first (IV1-SIRV). |
| H12 | Descriptive rows could be read as verdicts (DROP's `c` halving is the goal's headline). | §5/§7 state they are descriptive; DROP stays opt-in whatever they show (its dev trade-off, r1130, is the record). |
| H13 | The first draft only logged this file's sha1 at the gates, so the text could drift before the runs. | G6 now stops unless the frozen text hashes to `frozen/PREREG.sha1`; `run_ho.sh` refrozen (the table below). |
| H14 | V2 (human_testis) and S0 come from one library; were they independent? | Yes: disjoint reads (S0 = V2's unmapped reads); stated in §5. |
| H15 | The NULL takes the most-3′ PAS-free cluster, not a random one; is it a straw man? | It mirrors the rule's own choice with only the PAS conjunct negated, and on dev it recovers 20 / 41 genes (vs 2 / 6 for BASE): a strong comparator, not a weak one. |

## Freeze record

Frozen at `/mnt/linuxdisk/tmp/rustle_figures/complete_ho/frozen/` (`sha1sum -c SHA1SUMS` passes; files read-only):

| file | sha1 |
|---|---|
| `tes_null.py` (the NULL, IV_pre, PARITY) | bdcb847857131cae4b4149773f81d82ea5f1a19d |
| **`ho_eval.py` (the scorer)** | **597fcd08e4e3e0298946878667f88687df549fb6** |
| `run_ho.sh` (gates and steps) | baca765f09bc004fae9d2bc3ead3bf2da83f80b8 |
| `INPUTS.tsv` | 00f401cca102ea657cba598070bfa761f719b5cb |
| `SHA1SUMS` | be71242183dc867e5db3cc91d3928bb9928f92d6 |

Binaries `tss3_bin_frozen/`: `copy_assign` b1709a96ed1cf04bbb065cc6817343d73b225e69, `as_table` 0136ecac, `family_score`
27aa9445, `mcl_families` f155f2c4. Imported scorer core `ct_frozen/complete_eval.py` 6461144397294350fe7727e2bf1a1f626032ed03.
This file's own sha1 (all lines down to and including the `## Amendments` heading) is recorded outside it, in
`frozen/PREREG.sha1` (format: sha1, line count), and checked by gate G6.

## Amendments

**Amendment 1 (2026-09-27 22:10, before any held-out or SIRV number of any arm was scored or read).**
- **G0, user acceptance: GIVEN.** Asked in the main session ("Accept reusing the six held-out samples (dev contigs
  excluded) for the pas-end verification, recorded as a further reuse in the prereg?"), the user answered "Yes, accept
  and score". This is a further reuse of the six samples, after the readthrough R, R3 and the complete-transcripts
  arm-A acceptance; the exposure is as §0 states.
- **Gates G1-G6: PASS** (`/mnt/linuxdisk/tmp/rustle_figures/complete_ho/gates.txt`): frozen binaries/scripts/ct_frozen/
  driver/registry sha1s OK; BASE GTF sha1s of the six samples OK; tool GTFs OK; SIRV smoke seal unchanged (not opened);
  G5 human_testis BASE re-run with `tss3_bin_frozen` (options unset) byte-identical to `runs/human_testis/human_testis.gtf`
  (GTF, params.tsv, quant.tsv); scorer self-tests OK; the frozen prereg text (first 254 lines) sha1 a046be8c OK.
- **Deviations (recorded, none changes a rule, clause or arm):**
  - D1: every arm (PASEND, NULL, DROP, RESCUE on the six samples; the SIRV arms) was BUILT before G0 by
    `complete_ho/prebuild.sh`, which ran the frozen `run_ho.sh` functions with only `score` stubbed. No score was
    computed; `tables/` and `work/` were empty at 22:10.
  - D2: the human_A119b calls ran under an extra 588-s guard that never fired.
  - D3: the SIRV build logs showed |M| = 0 in both count modes (pas-end moves no SIRV end, as §2 predicted), and one
    log line showed the human_testis BASE transcript count (25,246). No other held-out number was read.
  - D4: the independent checker's author saw the header and 2 per-transcript rows of `chimp_PTR.NULL.tsv` (one `tie`,
    one `free` decision) while checking a file format; no count or class row.
- **Independent checker** (for the Outcome): `scratchpad/figs/ho_check.py` sha1 690e35cf, own parsers, validated on
  dev with 0 mismatches against the frozen scorer.

## Outcome (2026-09-27 22:25; scored after Amendment 1; independent recompute `ho_check.py` 690e35cf: 0 mismatches)

**Verdict, pas-end (judged): KEEP OPT-IN.** J1 (moved-set TES50 genes > BASE) passes 6/6; J2 (> matched NULL) passes
5/6; J3/J4 (`=` pairs, precision) pass everywhere; SIRV inert as predicted (|M| = 0: no SIRV end is internally primed),
S_EQ and S_TES pass, S_ACTIVE fails by construction. EFFECTIVE also required J1 and J2 on BOTH human samples;
human_testis fails J2 on a tie (13 vs 13 moved-set genes; only 65 transcripts are moved there, 0.3% of its multi-exon
output, the testis library has few internally primed ends). Moved-set TES50 genes BASE / pas-end / NULL: chimp 125 / 478
/ 419; KB3781 149 / 867 / 692; OR6737 114 / 674 / 588; A119b 249 / 883 / 711; testis 5 / 13 / 13; PPY 159 / 507 / 432.
Whole-substrate TES50 genes rise +101 / +205 / +155 / +128 / +1 / +67 (NULL +79 / +158 / +128 / +87 / +2 / +48).

**Descriptive (held-out, dev contigs and chr21/chr22 excluded):**
- **Retained intron `m` is the lowest of all arms wherever tools exist:** A119b ours 0.019 vs StringTie 0.053 / FLAIR
  0.037 / IsoSeq 0.034; OR6737 0.014 vs 0.040 / 0.029 / 0.042. SIRV: `m` = 0 for every arm.
- **Contained fragments `c` are the highest:** A119b 0.074 vs 0.011 / 0.015 / 0.036; OR6737 0.084 vs 0.027 / 0.022 /
  0.061. Partial share (c+k+m+n) A119b 0.148 vs 0.177 / 0.134 / 0.149; OR6737 0.138 vs 0.141 / 0.109 / 0.175.
- **`--polish-subchain drop` halves `c` on every substrate** (A119b 0.074 → 0.037; testis 0.060 → 0.051; OR6737 0.084 →
  0.047; KB3781 0.058 → 0.034; chimp 0.103 → 0.068; PPY 0.089 → 0.052) and brings the partial share to 0.115 / 0.104 on
  A119b / OR6737, the lowest of all arms. Cost: matched chains −98 / −17 / −26 / −25 / −22 / −37 (≤ 0.3%); genes losing
  every `=`/`c` query 134 / 28 / 60 / 41 / 74 / 189 (0.3–1.4%); TES50 genes −12 / −2 / −25 / −13 / −10 / −18.
  **SIRV (known complete isoforms, dedup): `c` 0.164 → 0.038 with all 48 truth `=` isoforms kept (48 → 48).**
- **`--polish-tss rescue`:** human A119b +1,711 transcripts, +94 chains, +13 TSS250 genes; chimp +245 / +19; the gorilla
  samples and testis unchanged (no cap signal, or nothing proven in scope).
- **SIRV truth isoforms recovered (`=`, of 61, dedup):** ours 48 in every arm; StringTie 39, FLAIR 49, IsoSeq 47.

**Reading.** On held-out data the default already emits the fewest retained-intron transcripts; its remaining partials
are contained 5′-truncated fragments, which the opt-in drop halves at ≤ 0.3% of chains and 0.3–1.4% of genes' only
compatible query, losing no SIRV truth isoform. pas-end improves 3′ ends selectively on 5/6 samples but misses the
pre-registered both-human requirement on a 65-transcript tie; it stays opt-in. Full tables:
`/mnt/linuxdisk/tmp/rustle_figures/complete_ho/check/outcome_rows.md` and `complete_ho/tables/`.
