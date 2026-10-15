# PREREG — `candidate_copies` few-read false flags: does a clean simulation produce them, and by which mechanism? (2026-10-14)

**Written before any `candidate_copies` run on the worlds below.** One generation-only probe (`synth_world.py build
--spec-world --seed 7`, no alignment, no stage run) was made to read the world's file layout; seed 7 is used nowhere below.

## The open item

Register 1247 / `docs/O3_CANDIDATES_AMENDMENT15_2026-10-05.md`: after the Amendment 15 consensus fix the opt-in stage
still fails both no-deletion controls — A14 C1' 10/53 families with a b/c/pri flag (bar ≤ 8), held-out H1 6/30 (bar ≤ 4).
The remaining false flags are described as "few-read (median 11.5 reads), mismatch-dominated consensus (median whole-length
d 0.036, median id x cov 0.974) — a future amendment's subject". Two mechanisms are suggested by the code and have never
been measured:
- **M1 mixing**: single linkage at de ≤ delta = 0.00958 chains reads of two copies; the column-wise vote gives a mosaic
  more than delta from every locus.
- **M2 consensus error**: reads of ONE present copy, but columns covered by < 3 members are never voted (`VOTE_MIN_COVER`)
  and ends covered by 2 survive the trim, so a 3-12-read cluster keeps its template read's errors and end artifacts.

A third explanation is not a stage defect: **M3 foreign transcripts**. Both controls use the KB3781 fibroblast library, and
the later unmapped-rescue study (`docs/PREREG_unmapped_rescue_2026-10-08.md`, Amendments 41-46) found about 4 in 5 of
that library's reference-absent flags to be other apes' transcripts. A median whole-length d of 0.036 is in the range of a
gorilla-orangutan transcript difference. Such reads ARE absent from the gorilla reference, so no RNA-internal rule can
separate them from a real absent copy at equal support; only a foreign-genome screen can. A simulation contains no M3 by
construction — which is what makes it the right instrument to ask whether M1/M2 exist at all.

## Phase 1 (this file): a diagnostic, not a rule

**Worlds.** `bench/o3_fewread/world.py` builds, per seed, one synthetic world with `bench/famsim` through
`synth_world.build_world` (random 60-kb background per family, synthetic 5-10-exon genes, the shipped HiFi read model
err 0.001 / indel 0.0003), 48 families:

| type | n | copies (in the reference unless marked) | expression (reads per copy) |
|---|---|---|---|
| NULL3 | 8 | A, B (2% SNPs), C (3%) | 80 each |
| NULLclose | 8 | A, B (1%), C (1.5%) — inside single-linkage chaining reach | 80 each |
| NULLlow | 8 | A, B (2%), C (3%) | A, B 80; C one of 4, 6, 8, 12, 16, 24, 32, 48 (one per family) |
| ALLELE@d | 10 | A, B (3%), H = A + d SNPs, **H absent** (d = 0.25, 0.5, 1, 2, 4%; 2 each) | 80 each |
| ERASED | 14 | A (2%), B (3%), **E absent** | A, B 80; E one of 4, 6, 8, 12, 16, 40, 80 (2 each) |

Read variants (both run): **E0** (famsim reads, 5' jitter 30) and **E3** (jitter 3 + the untemplated 5' G runs measured on
the real reads, Amendment 13: `synth_world.add_g5`). Seeds: **1 and 2** for Phase 1. Seeds 3, 4, 5 are reserved, unbuilt,
for a Phase 2 held-out test.

**Pipeline per world × variant.** Reads mapped with the shipped `sim.MM2` command to `genome.ref.fa`; a copies table
(`P.fam.copies.tsv` / `.fa` in `mcl_families`' layout: one family per synthetic family, one copy per IN-REFERENCE copy,
its truth exons, its spliced exon sum) built from the world's truth; `candidate_copies` at its defaults
(`--delta 0.00958 --min-cluster 3 --min-support 6`); minimap2 2.30, binary built `--release` from this file's commit.

**Scoring (`bench/o3_fewread/score.py`), per flagged candidate:** its reads (`reads.tsv`) mapped to their true copy (read
names carry it); `top_frac` = the share of the most common source copy; that copy's kind; the candidate's union sequence
(`contigs.fa`) aligned to every true transcript (`minimap2 -c -x splice:hq -uf`), best hit by identity x coverage. Class:
- **true**: ERASED family, `top_frac` ≥ 0.9 from E (the planted absent copy).
- **allele**: ALLELE family, `top_frac` ≥ 0.9 from H (correctly detected, but class b for C1'/H1).
- **M1 mixed**: `top_frac` < 0.9.
- **M2 pure-present**: `top_frac` ≥ 0.9 from an IN-REFERENCE copy.

A flag is **few-read** iff `n_reads` ≤ 20 (the real false flags' median is 11.5).

**Reported:** per type × variant × seed, families with ≥ 1 flag by class; every M1/M2 flag with n_reads, n_clusters, d,
top_frac, best truth hit; ERASED recall (families whose E is flagged) by E's expression; the n_reads of every true flag (the
cost curve any read floor would impose).

## Branching rule, fixed now

Let **F** = the number of (family, variant, seed) units of the three NULL types (8 × 3 × 2 variants × 2 seeds = 96 units)
carrying ≥ 1 few-read M1 or M2 flag.

- **F ≥ 5 (≥ 5.2%; the real H1 bar is 4/30 = 13%)** → the simulation reproduces the regime. Phase 2 = a NEW amendment to
  this file, written before seeds 3-5 are built, fixing (a) the rule chosen on seeds 1-2 only, (b) its held-out bars on
  seeds 3-5 (NULL-type false-flag units and ERASED recall per expression level), and (c) the real-data acceptance on the
  registered A14 / H substrates with the unchanged C1' / H1 / A13-1 bars, to run on the main machine.
- **F < 5** → the clean simulation does not reproduce the few-read false flags; no rule is developed on simulation (it
  would be fitted to nothing). The next step is the **M3 cause check on the real flags** (main machine, files
  `/mnt/linuxdisk/tmp/rna_allele/a14a15/`, `a15h/`): align the 16 + 10 false-flag unions to CHM13, mPanTro3 and mPonPyg2
  with the Addendum-1 foreign rule (coverage ≥ 0.8, identity ≥ 0.995 and above the gorilla host identity by ≥ d/2), pre-
  registered separately before that alignment is made.
- Whatever F is, any M1 or M2 flag is diagnosed individually (consensus vs template, per-column cover) and reported; a
  defect found that way gets its own fix and regression test, not a threshold.

Nothing in `candidate_copies` is changed before Phase 2 is registered.

# OUTCOME (2026-10-14) — F = 0: the clean simulation does not reproduce the few-read false flags

`candidate_copies` built `--release` at 49e62be3; `bench/o3_fewread/world.py build|run` and `score.py` under `pixi run`
(minimap2 2.30, samtools 1.22.1). Each world × variant: 48 families, about 130 s, peak RSS 0.24-0.34 GB. Run products:
`~/rustle_data/o3_fewread/w_s{1,2}/{E0,E3}/` on the laptop (`flags.tsv` = one row per flagged candidate).

| type (per world) | E0 s1 | E3 s1 | E0 s2 | E3 s2 |
|---|---|---|---|---|
| NULL3 / NULLclose / NULLlow: families with any flag | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| ALLELE (10): allele flags | 4 | 4 (+1 M1) | 6 | 6 |
| ERASED (14): E flagged | 12 | 12 | 12 | 12 |

- **F = 0 of 96 NULL-type units** (no M1 or M2 flag of any read count in a NULL family, including NULLclose's 1-1.5%
  paralogues inside chaining reach and NULLlow's present copies at 4-48 reads). **Branch: F < 5.** No rule is developed
  on simulation; the next step is the M3 cause check on the real flags (`docs/PREREG_o3_candidates_foreign_check_2026-10-14.md`).
- **The one M1 flag** is not few-read: `f29` (ALLELE@0.01, seed 1, E3), one cluster of 131 reads = 80 of the absent allele H
  + 51 of the present copy A (H is 1% from A, just above delta = 0.00958), best hit H at id x cov 0.9979. A chained
  allele+copy cluster, reported, not a few-read flag.
- **ERASED recall by E's reads (pooled over 4 units):** 4 reads 0/8 (below `--min-support 6`, by design), 6-80 reads 48/48.
  All 48 true flags are read-pure (top_frac 1.0), best hit their own E, id x cov ≥ 0.9969 (32 at ≥ 0.999). The 32 true flags
  at ≤ 20 reads have id x cov ≥ 0.9978: **on clean reads a 6-16-read cluster gives an accurate consensus**, so a read floor
  above 6 would remove real low-expression copies here without removing any false flag.
- **Allele flags:** 20, every one at d ≥ 0.01 (above delta), none at d ≤ 0.005 (linked, as designed); best hit H in 20/20.

**Reading.** On famsim reads (uniform HiFi errors, synthetic genes, one or two paralogues), neither mechanism M1 nor M2
produces few-read false flags. The real ones therefore come from something the simulation lacks: foreign or contaminant
transcripts (M3), real-read error structure (homopolymer indels, chimeras, incomplete processing), or real paralogue
complexity (many copies, gene conversion). M3 is the one with independent evidence on the same library, and it is directly
testable on the real flags; it goes first.
