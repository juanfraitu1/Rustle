# Pre-registration: readthrough junction filter v3 (arm R3): v2 with a reads-only majority guard on tier B

**Written 2026-09-26, before any held-out number of arm R2, arm R3 or their NULLs exists.** User goal (`/goal`):
*"develop effective readthrough reads filter"*. This file is the successor of `docs/archive/2026-09/PREREG_readthrough_v2_2026-09-26.md`
(the "v2 prereg"), which it **supersedes before v2's held-out stage**. It fixes one rule (R3) and its matched NULL
(NULL3). Everything else is v2's: the R prereg's metrics, substrates, power floors and scorer, v2's clause D and v2's
verdict logic, with R3 in place of R2. Where a v2 section still applies it is cited, not copied; every binding
constant and clause formula is copied verbatim. A default flip is the user's call whatever the outcome.

## 0. What was seen before this file, and the procedural record

- **Seen (development contigs only, design evidence).** Human A119b chr16 and chr20 and gorilla OR6737 NC_073244.2,
  and nothing else:
  - everything v2 §0 lists;
  - v2's own dev arms (`rt3_impl.md`): R2 assembled with the frozen `copy_assign` cdfad023. Clause D fails on both
    human dev contigs on chains alone: chr16 1,595 vs R 1,609 (bar 1,606), chr20 1,052 vs 1,060 (bar 1,058);
  - `rt4_tierb_guard.md`: what tier B costs on human dev, and the guard of §1 (session `scratchpad/figs/`, scripts
    and lists under `/mnt/linuxdisk/tmp/rustle_figures_dev/rt4_tierb_guard/`);
  - `rt4_ghost.md`: the regroup-after-polish lever (§10; not part of this test).
- **Seen (held-out, already published):** only what v2 §0 lists, i.e. the R prereg's Outcome, its addendum and the
  BASE / R rows of `rt_arms/tables/`. This author read the `a.fused` and `c.matching_intron_chains` rows of BASE, R and
  R's NULL there, for the projection of §2.5. They are the numbers v2 §6.4 is derived from.
- **Does not exist / not seen:**
  - no held-out number of R2, R3, NULL2 or NULL3: no held-out assembly of either rule has been run;
  - no assembly by a binary that holds the `r3` switch. The dev numbers of R3 in §2 are the list arm `GKmaj` of
    `rt4_tierb_guard`. That arm is `RUSTLE_READTHROUGH_JUNCTIONS=list:<file>` on cdfad023, where the file is R2's
    flagged set minus the junctions the guard protects. §2.4 shows that the scorer's port of the §1 rule reproduces
    that set exactly;
  - the parallel implementation step (`rt5_impl`) was editing `src/` while this file was written. Its in-progress
    code was read only to check that its dump format matches the scorer's port (column `N`, `NA` off tier B, tier
    `B_guarded`); none of its outputs was read. §4 binds whatever it produces.
- **Procedural record (plain).**
  1. **v2 §4 bound the author** with: "Dev results cannot change the rule, the NULL, a clause or a tolerance. The
     held-out stage runs whatever dev shows, unless the user stops the work. A stop would be recorded, and nothing
     would be re-tuned." v3 is a re-tune, and this file is its record:
     - after R2's dev arms showed clause D failing on both human dev contigs, a guard on tier B was **selected on the
       same three development contigs**;
     - it was chosen from **29 cost junctions** (tier-B-only junctions carried by the 29 matched reference chains R2
       loses against R) and **17 gain junctions** (the tier-B junctions that buy R2's extra fused-loci reduction);
     - it was compared with about ten alternative guards on those same contigs (`rt4_tierb_guard` §2).
  2. **v2 is superseded before its held-out stage.** Its held-out test (R2 judged against BASE and R) is never run,
     v2 never gets a verdict, and its §8 predictions are never scored. v2 stays on record as a dev arm (its dev table
     is in `rt3_impl.md`). R2's held-out rows produced beside R3 (§3) are **descriptive only**: they never enter a
     verdict and never become a v2 verdict.
  3. **Every dev number in this file is design evidence only.** It was selected on the contigs that measure it, so it
     is the optimistic edge. The guard's dev pass of clause D is by construction, not evidence.
  4. **What keeps the held-out stage a test:**
     - no held-out number of R2 or R3 exists when this file is frozen;
     - the guard adds no constant: it is a majority test of tier B's own premise (§2.2);
     - every clause, tolerance, floor and substrate is v2's, fixed before the guard existed.
  5. **Exposure.** The six held-out substrates served R's verdict (09-25/26). v2 never used them, so R3 is their
     second readthrough verdict. v2 §5's "third reuse" clause covers a rule chosen after seeing R2's held-out results,
     and that is not the case here. **A v4 chosen after R3's held-out results (or after R2's descriptive rows) needs a
     new library, or the user's explicit acceptance of a third reuse, written into that prereg.**

## 1. The v3 rule (binding, verbatim)

The rule is stated in integer arithmetic. It applies per canonical junction J with S ≥ 2, counted from PRIMARY spliced
reads. It is v2's rule with one extra conjunct in tier B:

```
S, U, V1, L exactly as the v2 prereg §1 defines and realises them (same population, same frozen statistics).
N  = N_span = same-strand SPLICED primaries (the population S, U and V1 are counted from) whose 5' end is upstream of
     J's donor (U's convention) and whose 3' end is beyond J's acceptor (V1's "beyond" convention), J reads included
     (N = S + K, K = those that do not use J).
tierA (the shipped R)  : U >= 20*S
exempt (ALE protection): 2*L >= S AND 3*V1 < 5*S
tierB (own promoter)   : U >= S AND V1 >= 4*S AND V1 > N      (the tie V1 = N is NOT flagged)
flag(J) = (tierA AND NOT exempt) OR tierB;  the chains (reads) carrying a flagged J are dropped before locus formation,
exactly as the r and r2 modes do. Name: RUSTLE_READTHROUGH_JUNCTIONS=r3.
```

**How the text is realised (no freedom left to the implementer):**
1. **S, U, V1, L, scope.** These are v2 §1 items 1-3 unchanged, as computed by the frozen cdfad023: canonical on its
   strand, S ≥ 2, primary = not secondary and not supplementary, spliced = ≥ 1 intron from the pool's own CIGAR
   parser. On every junction, R3's binary must report the same S, U, V1, L and the same tier A, exempt and tier-B-
   candidate bits as R2's (§4 IV2 (i)).
2. **N.**
   - Coordinates: d = J's donor (0-based first intron base) and a = J's acceptor (0-based first base after the intron).
     p5 and p3 are the read's 5′ and 3′ aligned bases in transcript orientation (0-based), as the rows S, U and V1 are
     built from.
   - Test, on `+`: p5 < d AND p3 ≥ a. On `−`: p5 ≥ a AND p3 < d. These are U's "upstream of the donor" test and V1's
     "beyond the acceptor" test. On both strands they mean: the genomic first aligned base < d and the genomic last
     aligned base ≥ a.
   - Population: the same per-region, same-strand spliced primaries as S, U and V1 (the window rule included). A read
     carrying J always counts (so N ≥ S), and so does a spliced read with any other intron chain across J; an
     unspliced read never counts (it is not in the population), nor does a read whose last aligned base is the
     intron's last base.
3. **N is needed only on a tier-B candidate** (in scope, U ≥ S and V1 ≥ 4S). The dumps write `NA` elsewhere.
4. **Identities** (a unit test checks each):
   - (i) tierB implies NOT exempt (v2 §1 item 4), so flag(J) = (tierA OR tierB) AND NOT exempt;
   - (ii) flag_v3 implies flag_v2: R3's flagged set is R2's flagged set minus the tier-B-only candidates with V1 ≤ N;
   - (iii) V1 = N is not flagged, and V1 = N + 1 is.
5. **Action.** v2 §1 item 5 unchanged. The removal code path is R2's.
6. **Outputs.**
   - v2 §1 item 6, plus a column `N` in `<out>.readthrough_junctions.tsv` and in the ALL dump
     (`RUSTLE_READTHROUGH_JUNCTIONS_ALL=1`).
   - `tier` is `A` (tier A, not exempt; takes precedence), `B` (tier B, guard passed), `B_guarded` (a tier-B candidate
     the guard protects, NOT flagged: exactly the rows R3 removes from R2), `A_exempt` or `-`. The `B_guarded` rows
     must be recoverable for the census (§7.1); the ALL dump suffices.
   - Under unset, `off`, `r`, `rq1`, **`r2`** and `list:`, every product stays byte-identical to cdfad023 (§4 IV1).
     r2 stays exactly as it is: it is a dev arm on record.

## 2. Why the guard (development diagnosis only; design evidence, selected on these contigs)

### 2.1 What tier B cost on human dev (`rt4_tierb_guard` answers 1-2)
- **Direct, and never a readthrough.** R2 loses 29 matched reference chains against R (19 / 10 / 0). Each holds
  exactly ONE R2-flagged junction, always tier-B-only and always an annotated intron of one gene.
- **J sits inside ONE transcription unit.** Many molecules join J's donor side to its acceptor side without using J:
  K = 15-1,319 on the 29 cost junctions, against 0-30 on the 17 gain junctions. On 18/29 the bypass skips an exon
  inside J's intron; the rest are alternative first exons or other splice-site pairs that span J. In 23/29, the V1
  reads that tier B reads as "B's own promoter" start on an internal exon of the bypass molecules: 5′-truncated
  through-molecules, or an alternative promoter of the same unit.

### 2.2 The guard: a majority test of tier B's own premise
- Tier B assumes that B is an independently initiated gene (v2 §2.4). `V1 > N_span` reads: "B's own starts outnumber
  all molecules that reach B's side from A's side".
- No constant is added. The weaker and stronger variants (c = 1/2: 2·V1 > N; c = 2: V1 > 2·N) give the same picture on
  dev (§2.3), and neither is an arm here.
- The bypass ratios K/V1, K/S and K/U separate cost from gain junctions at AUC up to .99 (n = 29 / 17). The other
  listed read features are weak to moderate (.12-.85, `rt4_tierb_guard` answer 3).

### 2.3 Dev table (full-BAM setting of v2 §2.1, seeded; R3 = the `GKmaj` list arm on cdfad023; scorer 2de04fd1)

| arm | chr16: fused / matched chains / lost vs BASE / precision | chr20 | gorilla NC_073244.2 |
|---|---|---|---|
| BASE | 101 / 1,605 / – / .1828 | 60 / 1,059 / – / .2098 | 53 / 1,595 / – / .4091 |
| R | 77 / 1,609 / 8 / .1849 | 47 / 1,060 / 2 / .2134 | 39 / 1,588 / 7 / .4130 |
| R2 | 69 / 1,595 / 23 / .1843 (**D fail**) | 41 / 1,052 / 10 / .2142 (**D fail**) | 37 / 1,597 / 0 / .4139 (D pass) |
| **R3** | **69 / 1,612 / 6 / .1836 (D pass)** | **41 / 1,062 / 0 / .2126 (D pass)** | **37 / 1,596 / 0 / .4130 (D pass)** |
| NULL of R3 | 98 / 1,553 / 53 | 59 / 1,035 / 24 | 50 / 1,580 / 16 |

- **A1** (fused reduction vs BASE): R3 −31.7% / −31.7% / −30.2%, the same as R2. The matched NULL gives −3.0% /
  −1.7% / −5.7%.
- **G1-G4** pass on all three:
  - G3 TES genes: 547 vs 517, .2890 vs .2736, 579 vs 574;
  - G4.annotated: 964 vs 963, pass, 867 vs 869.
- **Against R2:** the fused-loci set is identical on all three (0 pairs differ); chains +17 / +10 / −1; flagged
  junctions 1,762 / 877 / 122 → 1,024 / 459 / 106 (`B_guarded` 738 / 418 / 16); removed pool alignments 5,943 /
  2,774 / 397 → 3,751 / 1,564 / 363.
- **Junction precision** RT/(RT+ANN), labels only: R3 .794 / .783 / .885, against R .714 / .719 / .778 and R2 .609 /
  .579 / .836. The protected junctions are RT 3 / 5 / 2 against ANN 67 / 36 / 4.
- **Chain precision** is above BASE and below R on human (.1836 vs .1849, .2126 vs .2134). The drop comes from the
  fixed exemption (R2 without tier B: .1830 / .2116), not from the guard.
- **Specificity controls.** Scaling the ratio (c = 1/2: 2·V1 > N; c = 2: V1 > 2·N) gives 69-70 / 41 / 37 fused and
  1,609-1,613 / 1,060-1,062 / 1,596 chains; D passes at both. A random protection of the same number of tier-B
  junctions (matched by log2 S) keeps only 1 / 2 / 3 of the 9 / 7 / 3 fused loci tier B resolves.
- ⚠ **Thin margins:** closest gain ZNF205-AS1|ZNF213-AS1 (N 9 vs V1 10); closest protected cost CARHSP1 (N 237 vs
  V1 207).
- **What the guard cannot fix.** 2 of the 29 costs remain, TUBB3 (V1 92 vs N 17) and LOC124907834 (126 vs 70): both
  are two-promoter units whose internal promoter dominates. The guard gives back RT-labelled "dominant fusions" with
  N ≥ V1 (3 / 5 / 2 junctions), consistent with R's policy of keeping them, and re-loses one gorilla chain R2 had
  gained (LOC101145690, K 95 vs V1 12).
- **The chr16 loci-count jump (2,729 → 2,837) is not the guard** but the contig-level mono-exon polish-quantile cliff
  (mono-exon loci with ≤ 10 reads: 83 in BASE, 0 in R and R2, 84 in R3).

### 2.4 Port check (this file's author; dev only)
- **The scorer's v3 port reproduces the `GKmaj` lists exactly.** The port (`readthrough_eval.py` sha1 f6b99dcc,
  `v3_flag` / `tiers`) was run on synthetic R3 ALL dumps: `rt3_impl`'s R2 ALL dumps plus N = S + K on the tier-B
  candidates, K from `rt4_tierb_guard/sig.<sub>.tsv`. Flagged 1,024 / 459 / 106 with 0 port-only, 0 list-only and 0
  mismatches; `B_guarded` 738 / 418 / 16. `null` on those dumps reproduces `rt4_tierb_guard`'s NULL list byte for
  byte (body and summary) on NC_073244.2 and chr20.
- ⚠ **This checks the port, not the implementation.** `rt4`'s K counted all primaries, while N counts spliced
  primaries only. `rt4_tierb_guard` §3 reports identical guard decisions for the two on all three contigs. §4 IV3
  checks the real binary.

### 2.5 Projection to held-out (arithmetic on dev; NOT a measurement, NOT in any verdict)
- R3's dev gain over R is +7.9 / +10.0 points of fused reduction on human and +3.8 on gorilla; with R's dev→held-out
  shrinkage of about 0.73 (v2 §2.1) that is about +6 / +7 and +2.8. R3 − R matched chains on dev: +3 (+0.19%) / +2
  (+0.19%) / +8 (+0.50%).

| held-out substrate | FUSED BASE → R (published) | R3 projected | chains BASE → R (published) | R3 projected vs R |
|---|---|---|---|---|
| human_A119b −chr16/chr20 | 2,050 → 1,626 (−20.7%) | ~−27% (~1,500) | 35,956 → 35,796 | ≥ R − 0.1% |
| human_testis | 200 → 165 (−17.5%) | ~−24% (~152) | 11,399 → 11,357 | ≥ R − 0.1% |
| gorilla_OR6737 −NC_073244.2 | 549 → 454 (−17.3%) | ~−20% (~439) | 24,277 → 24,180 | ≥ R − 0.1% |
| gorilla_KB3781 | 662 → 466 (−29.6%) | ~−32% (~448) | 25,911 → 25,709 | ≥ R − 0.1% |
| chimp_PTR | 420 → 380 (−9.5%) | ~−12% (~368; A1 needs ≤ 378) | 21,031 → 20,952 | ≥ R − 0.1% |
| orangutan_PPY | 687 → 579 (−15.7%) | ~−18% (~560) | 20,243 → 20,175 | ≥ R − 0.1% |

- **This replaces v2 §2.6's warning.** Against D (§6.2), the projection puts both parts of D inside tolerance on all
  six substrates. That is expected, since the guard was selected for it (§0 item 3).
- **The projection assumes** that R3's dev chain gain over R shrinks toward zero rather than turning into a loss, that
  the guard keeps all of tier B's locus gain (as on dev), and that chimp's residual after R holds own-promoter
  junctions in gorilla's proportion.

## 3. Arms

| arm | what | status |
|---|---|---|
| **BASE** | current default; `rt_arms/<sample>/<sample>.BASE.gtf` (= `runs/<id>/<id>.gtf`), restricted | **reused** (v2 §3): the held-out arm GTFs and families of the R prereg; no re-assembly |
| **R** | `RUSTLE_READTHROUGH_JUNCTIONS=r`, frozen `copy_assign` 452454c6 | **reused** (v2 §3): `rt_arms/<sample>/<sample>.R.gtf` and, on human_testis and chimp_PTR, its families |
| **R3** | `RUSTLE_READTHROUGH_JUNCTIONS=r3` (§1), frozen `/mnt/linuxdisk/tmp/rustle_figures/rt3_bin_frozen/`, `copy_assign` sha1 f480847abefefc7ff7997c528c2fd4dce5bda57b | **new, judged**; genome-wide on the six samples: `rt_arms/<sample>/<sample>.R3.*` |
| **NULL3** | matched random removal to R3's flagged reads (§3.1), run as `list:` with R3's binary | **new**: `rt_arms/<sample>/<sample>.null3.tsv` and `<sample>.NULL3.*` |
| **R2** | `RUSTLE_READTHROUGH_JUNCTIONS=r2`, frozen `rt2_bin_frozen/` `copy_assign` cdfad023 | **new, DESCRIPTIVE ONLY** (never in a verdict): `rt_arms/<sample>/<sample>.R2.*` |

- **The reuse of BASE and R** is valid only under IV1 (§4), exactly as v2 §3 says. R's NULL (`<sample>.NULL.gtf`) may
  sit in the score call as the shared NULL of arm R's rows. It is not a comparator of any clause of this file.
- **R2's only question is: does the guard transfer?** It is answered by D(R2 vs R) beside D(R3 vs R), and by the
  R3 − R2 chains and fused loci per substrate (§7.1). R2 gets no NULL2 and no families. Its A1-G5 rows are printed
  and never read as a verdict.
- **Not run on held-out (by design).**
  - No component arms (v2 §3): no guard variant (c = 1/2, 2, K ≥ U, K ≥ S, the order test) and no "R3 minus the
    exemption".
  - The held-out stage answers "R3 vs BASE" and "R3 vs R" only. It must not become a design tool for a v4.

### 3.1 NULL3
v2 §3.1 unchanged, with F = R3's flagged set:
- **Draw.** `readthrough_eval.py null --flags <R3's <out>.readthrough_junctions.tsv>` per sample, with the default seed
  20260925 and one RNG per `sample:contig`.
- **Candidates** are the canonical S ≥ 2 junctions that R3 does not flag. That includes the exempted and the
  `B_guarded` junctions: both are readthrough-like, so drawing them can only make NULL3 stronger. This is conservative
  against R3 and accepted.
- **Target.** The draw targets the distinct primaries that carry an R3-flagged junction, matched by floor(log2 S).
- **Arm.** NULL3 gets metrics (a)-(d). A1's null part for R3 is judged against NULL3 only.

## 4. Implementation gates (dev only; all before any held-out R3 run)

Any failure is fixed in the code, never in the rule.
- **IV1, byte identity.** Unset, `off`, `r`, `rq1`, `r2` (with and without the ALL dump) and `list:`: every product
  that v2 §4 IV1 lists is `cmp`-identical to the frozen cdfad023 outputs, in the full-BAM setting on the three dev
  contigs, seeded and not seeded. `cargo test --release` passes.
- **IV2, port = spec.** On the three dev contigs:
  - (i) R3's ALL dump has the same S, U, V1, L and canonical values as R2's cdfad023 ALL dump on every row. Its tier
    equals R2's on every row, except that R2's `B` becomes `B` or `B_guarded`.
  - (ii) N equals an independent pysam recount of §1 item 2 on every tier-B candidate. The recount uses the full BAM,
    `-F 2308`, spliced primaries, the same strand rule, and the pool's CIGAR parse emulated as in `rt3_xcheck.py`. The
    only allowed differences are rows explained by the R prereg's Amendment 1 item 5, and each one is listed.
  - (iii) `readthrough_eval.py tiers` on R3's flags dump and ALL dump reports 0 mismatches and 0 `N_missing`. The
    `N_lt_S` count is reported.
  - (iv) Unit tests: V1 = N not flagged and V1 = N + 1 flagged; a spanning read with another intron chain counts; a J
    read counts; an unspliced read does not; a read ending on the intron's last base does not, one ending on the
    acceptor exon's first base does; the `−` strand; streaming = buffered; the identities of §1 item 4; tier A takes
    precedence over a guarded candidate.
- **IV3, action and the selected set.**
  - R3's flagged set on the three dev contigs equals `rt4_tierb_guard/lists/<sub>.GKmaj.tsv` (1,024 / 459 / 106). Any
    difference must come from the K-over-all-primaries vs N-over-spliced-primaries population, and each one is
    listed. An unexplained difference stops the work: the binary would not realise the rule that was selected.
  - 0 transcripts of R3's GTF (and of NULL3's) use a flagged (listed) junction.
  - `list:` of R3's own dump reproduces R3's GTF, quant, families, assignments and famcn byte for byte.
  - The `[readthrough]` log line names the arm, and the run cache is keyed on the switch.
- **Dev arms (reported, never in the verdict).** R3 and NULL3 on the three dev contigs, scored with
  `readthrough_eval.py score` (not `--heldout`) with `--null R3=NULL3 --versus R`; the table goes into Amendment 1
  beside §2.3. Dev results cannot change v3's rule, NULL, clauses or tolerances: the held-out stage runs whatever dev
  shows, unless the user stops the work. A stop is recorded, nothing is re-tuned, and a v4 is a new prereg (§0 item 5).
- **Freeze.** Amendment 1 records the following before any held-out R3 run:
  - the `copy_assign` sha1 serving R3 and NULL3. It replaces the placeholder token in §3's R3 row (the only
    placeholder in this file); the binary is copied to `rt3_bin_frozen/` with a `SHA1SUMS`;
  - the scorer sha1;
  - check (a) of v2 §4: re-scoring the reused BASE and R GTFs of the six held-out substrates reproduces every recorded
    row of `rt_arms/tables/`;
  - check (b): the scorer's selftest passes.
- **The scorer at the time of writing** is sha1 **f6b99dcc1a3f972427cf67bc3de3b846dd5c5f30**. Relative to v2's scorer
  2de04fd1 it adds only `v3_flag` / `v3_tier_label` / `v3_tier_ok`, the v3 check in `read_junction_table` and
  `tiers` (dispatched on an `N` column), `verdict --prereg v3` and one selftest group. `clauses` (A1-G5) is unchanged
  from 0a50ab11; `clause_d`, `clause_rows`, `score`, `null`, `verdict` and `verdict_v2` are unchanged from 2de04fd1,
  and their outputs on every existing held-out and dev table and R2 dump are byte-identical (report `rt5_prereg.md`).

## 5. Substrates

The six held-out verdict substrates are v2's (§5), each judged on its own; species are never pooled:

| id | substrate |
|---|---|
| V1 | human_A119b genome minus chr16, chr20 |
| V2 | human_testis whole genome |
| V3 | gorilla_OR6737 genome minus NC_073244.2 |
| V4 | gorilla_KB3781 whole genome |
| V5 | chimp_PTR whole genome |
| V6 | orangutan_PPY whole genome |

**Power floors** are the R prereg's §4: (a) BASE FUSED ≥ 50; (b) universe ≥ 1,000 genes; (c) BASE matched chains ≥
500; (d) ≥ 1,000 reference loci; (e) ≥ 30 families or pairs. They are properties of the reused BASE. As v2 §5
records, A1, G1, G2, G3 and G4.annotated qualify on all six, and G5 qualifies on human_testis and chimp_PTR.
G4.extra_copy is the exception (§6.3). v2 §5's honest note on exposure applies, with §0 item 5 above.

## 6. Metrics and clauses

### 6.1 The R prereg's clauses, unchanged (metrics (a)-(e) of its §3, same scorer code, same annotations, gffcompare 0.12.10)

| clause | measure | R3 passes iff |
|---|---|---|
| **A1** efficacy | FUSED loci (`a.fused`) | reduction vs BASE ≥ 10% **and** NULL3's reduction < half of R3's |
| **G1** precision | intron-chain precision (raw counts) | R3 ≥ BASE |
| **G2** sensitivity | matching reference intron chains (`c.matching_intron_chains`) | 100·R3 ≥ 99·BASE |
| **G3** ends | TES-recovered genes (fixed universe) | R3 ≥ BASE |
| **G4** loci kept | found read-supported reference loci, annotated (and extra copies, §6.3) | 100·R3 ≥ 99·BASE |
| **G5** families | human: Compara Primates bipartite F, sensitivity, precision; apes: Liftoff pair recall | each ≥ BASE − 0.005 |

X = 10%, Y = 1% and 0.005 are the R prereg's §5. They are not re-argued and not moved.

### 6.2 Clause D: R3 head to head with R (v2 §6.2, R3 in place of R2)

On every one of the six substrates:
- **D.fused:** FUSED(R3) ≤ FUSED(R). The tolerance is 0 loci.
- **D.chains:** 1000 · matched chains(R3) ≥ 998 · matched chains(R). That is 0.2% of R's count, on G2's measure.
- **D passes on a substrate iff both parts pass.**
- **A D failure is "larger than twice its tolerance"** iff FUSED(R3) > FUSED(R), or 1000 · chains(R3) < 996 ·
  chains(R).

Why D exists, why D.fused has no tolerance and why 0.2%: v2 §6.2, unchanged.

### 6.3 Measurement scope for G4 extra copies and G5
v2 §6.3 applies, with R3 in place of R2:
- **G5** is measured on human_testis and chimp_PTR only: human_testis with Compara Primates (426 families;
  `family_score` bipartite F, sensitivity and precision), chimp_PTR with Liftoff pair recall. Both use
  `tools/rustle_pipeline.sh families` on R3's GTF with the same `mcl_families` binaries and shard wrapper as the R
  prereg's addendum. BASE's and R's families there are reused.
- **G4.extra_copy** is reported, not judged, on chimp_PTR (below floor (d)) and not measured elsewhere.
- **"Not measured" elsewhere does NOT cap the verdict.** Any other judged clause that is not measured is not a pass.

### 6.4 The numbers the clauses compare against (v2 §6.4, verbatim; BASE and R are reused)

| substrate | A1: FUSED(R3) ≤ | D.fused: ≤ R | G2: chains ≥ | D.chains ≥ (2× tol ≥) | G1: precision ≥ | G3: TES ≥ | G4.annotated ≥ |
|---|---|---|---|---|---|---|---|
| human_testis | 180 | 165 | 11,286 | 11,335 (11,312) | .4596 | 5,069 | 5,652 |
| human_A119b −chr16/chr20 | 1,845 | 1,626 | 35,597 | 35,725 (35,653) | .1765 | 9,792 | 21,332 |
| gorilla_OR6737 −NC_073244.2 | 494 | 454 | 24,035 | 24,132 (24,084) | .3496 | 6,985 | 13,047 |
| gorilla_KB3781 | 595 | 466 | 25,652 | 25,658 (25,607) | .3594 | 6,177 | 12,275 |
| chimp_PTR | 378 | 380 | 20,821 | 20,911 (20,869) | .3915 | 8,378 | 12,345 |
| orangutan_PPY | 618 | 579 | 20,041 | 20,135 (20,095) | .1931 | 9,053 | 13,613 |

- **G1** is decided on raw counts (k·n′ ≥ k′·n), not on the 4-dp values printed here.
- **G5 floors:**
  - human_testis: F ≥ .268973, sensitivity ≥ .154939, precision ≥ .949545;
  - chimp_PTR: pair recall ≥ .112241.
- If the freeze re-score (§4) does not reproduce the source rows, the work stops.
- **Chimp's A1 needs FUSED(R3) ≤ 378**, i.e. two fewer fused loci than R's 380.

## 7. Verdicts (arm R3 only)

v2 §7 applies verbatim with R3 in place of R2 and NULL3 in place of NULL2.
- **Judged clauses.** P(s) = "every judged clause passes on substrate s". The judged clauses are A1 (with NULL3), G1,
  G2, G3 and G4.annotated, plus G5 on human_testis and chimp_PTR.
- **EFFECTIVE: adopt R3 as the default candidate. The flip is the user's call.** This needs all three:
  - P(s) on all six;
  - G5 measured and passing on both human_testis (F, sensitivity, precision) and chimp_PTR;
  - D passing on ≥ 5 of 6 with no D failure larger than twice its tolerance.
- **KEEP OPT-IN:** not EFFECTIVE, and every A1-G5 failure falls on at most one substrate. A judged clause that is not
  measured is not a pass: it blocks EFFECTIVE but is not counted as a failure.
- **REFUTE:** A1-G5 failures on ≥ 2 of the six.
- **Missing substrates.** A substrate missing from the tables makes R3 "undecided", unless ≥ 2 present substrates
  already fail (REFUTE).
- **Never pooled.** Species are never pooled, and dev never enters the verdict.
- **R2 is never judged.** Whatever its rows say, they neither make nor block R3's verdict. They are not a v2 verdict.

**Scorer realisation (binding on `verdict --prereg v3`).** It is v2's verdict logic (`verdict_v2`, unchanged) applied
to arm `R3` alone, after two name checks on R3's rows:
- R3's A1 row must name `NULL arm NULL3`, which is what `score --null R3=NULL3` writes. If it does not, A1 is **not a
  pass**, except a failure on the reduction alone (< 10%), which needs no NULL.
- R3's D row must read `vs R:`. If it does not, D is **not a pass**.

Every other arm with D rows (R2) prints as `descriptive` with its D per substrate. R prints as `not judged`, and an
absent R3 as `undecided`.

### 7.1 Reported beside the verdict (not clauses)
- **Everything in v2 §7.1**, for R3.
- **Guard census per substrate:** `B_guarded` junctions and the alignments they keep in the pool, by RT / ANN / novel
  label (labels only, as `rt2_protect` §1); the 20 most-read `B_guarded` junctions with N, V1 and gene pairs.
- **Guard transfer (R2, descriptive), per substrate:** D(R2 vs R) beside D(R3 vs R) (the scorer's "D per substrate"
  lines), chains(R3) − chains(R2) and FUSED(R3) − FUSED(R2).

## 8. Predictions (before any R3 number, dev or held-out; the dev numbers of §2 are the list-arm equivalent)

1. **A1.** R3 passes on 6/6. The median reduction is 20-26%, as §2.5 projects: the guard kept all of tier B's locus
   gain on dev, so this is v2's prediction.
   - **chimp_PTR:** 11-15%, a pass, but the least certain. If chimp's residual after R is mostly exon bridges or
     dominant fusions (v2 §2.5), it misses again.
   - **NULL3:** 1-8% (dev 1.7-5.7%).
2. **G1** passes on 6/6: above BASE, possibly below R on human, where the exemption returns unmatched chains (dev).
   **G3** passes on 6/6, within ±0.5% of R. **G4.annotated** passes on 6/6, with a loss ≤ 0.7%.
3. **G2** passes on 6/6. R3's chain loss vs BASE stays within R's held-out range (0.34-0.78%) ± 0.2 points. On
   gorilla_KB3781, D.chains (25,658) is stricter than G2 (25,652), so a D pass there implies a G2 pass.
4. **D.** D.fused passes on 6/6: R3 is below R by several loci on each (dev: 8 / 6 / 2). D.chains passes on ≥ 5 of 6,
   with R3 − R between −0.2% and +0.5%.
5. **G5:** |Δ| < 0.002 on both measured substrates; a pass. The largest family is unchanged.
6. **Guard transfer (R2, descriptive).**
   - D(R2 vs R) fails on both human substrates by more than twice its tolerance, while D(R3 vs R) passes on both.
   - chains(R3) − chains(R2) ≥ +0.5% of R2's on both human substrates (dev +1.07% / +0.95%).
   - |FUSED(R3) − FUSED(R2)| ≤ 1% of FUSED(BASE) on ≥ 5 of 6 (dev 0 / 0 / 0).
7. **Verdict:**

   | verdict | probability | what it takes |
   |---|---|---|
   | EFFECTIVE | ~0.40 | A1 holds on chimp, D holds on ≥ 5 |
   | KEEP OPT-IN | ~0.45 | typically chimp's A1 (R missed it at −9.5%) or one D.chains shortfall from the thin guard margins (§2.3) |
   | REFUTE | ~0.15 | failures on ≥ 2 substrates, e.g. chimp A1 plus G2 on a human substrate if tier B's human cost does not stay guarded |

   The dev evidence for the guard is selected (§0 item 3). A result worse than this is the expected direction of
   error, not a surprise.

## 9. What would falsify the design reasoning (reported whatever the verdict)

Each item names a design claim and the held-out observation that contradicts it. They are read from the §6-§7 measures
and the census; no extra held-out arm is run.

1. **"The guard keeps tier B's locus gain."** Falsified if the median over the six of (FUSED(R3) − FUSED(R2)) /
   FUSED(BASE) is ≥ 0.01. Dev: 0 / 0 / 0.
2. **"The guard removes tier B's chain cost where that cost exists."** On every substrate where D.chains(R2 vs R)
   fails, R3 must recover at least half of R2's shortfall: chains(R3) − chains(R2) ≥ ½ (chains(R) − chains(R2)). Dev:
   17 ≥ 7 on chr16 and 10 ≥ 4 on chr20.
   - Falsified if this fails on any such substrate.
   - If D.chains(R2 vs R) passes on all six, tier B's cost did not transfer. The guard's premise is then untested,
     not confirmed, and it is reported so.
3. **"The guard protects same-gene splices, not readthroughs."** Falsified if, among `B_guarded` junctions labelled RT
   or ANN, the RT ones are at least as many as the ANN ones on ≥ 2 substrates. Dev: RT 3 / 5 / 2 against ANN 67 /
   36 / 4.
4. **v2 §9 items 1, 2 and 4-7, with R3 in place of R2:**
   - (1) the median over the six of (FUSED(R) − FUSED(R3)) / FUSED(BASE) < 0.03;
   - (2) among R-flagged junctions that R3 exempts, those labelled RT are at least as many as those labelled ANN on
     ≥ 2 substrates; or R3 recovers fewer than half of R's lost reference chains, as the median over the six (R3's
     exemption is R2's);
   - (4) D.fused fails anywhere;
   - (5) NULL3's reduction ≥ half of R3's on any substrate;
   - (6) R3-fused loci whose gene pair is not fused in BASE exceed 5% of R3's fused loci on any substrate;
   - (7) R3's median reduction ≤ 17.4%.
5. **"The guarded tier B costs no chains beyond R"** (replaces v2 §9 item 3; dev R3 ≥ R on all three). Falsified if
   1000 · chains(R3) < 995 · chains(R) on ≥ 2 substrates.
6. **Not falsifiers:**
   - v2 §9 item 8;
   - a residual TUBB3-type cost: two-promoter units whose internal promoter dominates are outside any reads-only
     guard (§2.3).

## 10. Not in this test: the regroup-after-polish knob (`rt4_ghost.md`)

- **What it is.** Re-deriving `gene_id` after the polish splits a gene_id whose surviving transcripts form more than
  one shared-junction component. Such a gene_id is either a "ghost" (the polish dropped the only bridge) or a "tid
  collision" (two components got the same `DN_<chrom>_<start>_<n_exon>` string).
- **Dev numbers** (stable naming `_rgs`; fused loci before → after, chr16 / chr20 / gorilla): BASE 101 → 88, 60 → 55,
  53 → 49; R 77 → 74, 47 → 46, 39 → 37; R2 69 → 65, 41 → 41, 37 → 36. Every intron-chain row and G4 are identical by
  construction; TES genes rise slightly; chr16 families rise (R2 Compara F .633 → .667, on 17 families, below the 30
  floor). R3 was not regrouped on dev.
- **Why it is excluded here.** It cannot move clause D (chains are invariant). It is threshold-free and orthogonal to
  the junction rule: it acts on BASE too. Applied to one arm only it inflates A1 (R2 alone: 35.6 / 31.7 / 32.1%,
  against 26.1 / 25.5 / 26.5% when it is applied to every arm).
- **It gets its own knob and its own prereg later:** `--gtf-regroup`, opt-in, with clauses fused, G3 and families F,
  and the `a.fused` locus-count trap. No arm of this test is regrouped, and no regrouped GTF is scored here.

## 11. Order, stop rules, machine rules

1. **Implement.** Implement `r3` behind the switch. Run IV1-IV3 and the dev arms, freeze the binaries to
   `/mnt/linuxdisk/tmp/rustle_figures/rt3_bin_frozen/` with `SHA1SUMS`, and write Amendment 1: the sha1s (replacing
   §3's placeholder token), the IV results and the dev table. No held-out command runs before Amendment 1 exists.
2. **Held-out assemblies.** Run R3 genome-wide on the six samples (human_A119b, human_testis, gorilla_OR6737,
   gorilla_KB3781, chimp_PTR, orangutan_PPY) with `rt3_bin_frozen`, and R2 with `rt2_bin_frozen` (cdfad023). Both use
   the driver's `tools/rustle_pipeline.sh assemble` with the sample's BAM and FASTA, as BASE and R did. Outputs go to
   `/mnt/linuxdisk/tmp/rustle_figures/rt_arms/<sample>/<sample>.{R3,R2}.*`.
3. **NULL3.** Draw NULL3 from each sample's own R3 dump (§3.1), then run `NULL3` = `list:` with `rt3_bin_frozen`.
4. **Families.** Run the families stage on R3 for human_testis and chimp_PTR (§6.3). No other families are run.
5. **Score.** Run `readthrough_eval.py score --heldout` per substrate with the arms BASE, R, NULL (R's), R2, R3 and
   NULL3, and `--base BASE --null NULL --null R3=NULL3 --versus R`. The tables go to `rt_arms/tables_v3/`; the tables
   of R (`rt_arms/tables/`) are never overwritten.
6. **Verdict.** `readthrough_eval.py verdict --prereg v3 rt_arms/tables_v3/*.tsv`. The Outcome is appended here, and
   every number becomes a new register row. R2's rows are reported under "guard transfer", never as a v2 outcome.

**Stop rules.**
- No change to the code, the rule, the NULL or a clause after step 1's freeze. The only exception is a bug fix that
  re-runs every R3 and NULL3 arm and is recorded as an amendment.
- No guard variant is ever substituted for the §1 rule.
- A step that cannot be computed makes its clause "not measured", which caps the verdict (§6.3, §7) unless the clause
  is outside §6.3's scope.
- An unexplained IV3 difference from the `GKmaj` lists stops the work before step 2.

**Machine rules.**
- One heavy process at a time, in the foreground: `flock -w 900 /mnt/linuxdisk/tmp/rustle_heavy.lock timeout 600
  <cmd>`. Genome-wide steps are split into bounded, resumable calls (`--budget-s`, exit 75 = run the same command
  again).
- Outputs and `TMPDIR` go under `/mnt/linuxdisk`.
- Build only with `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release`, output
  captured to a file. Touch the edited `.rs` files first (/mnt/c timestamps).
- Take every binary from its frozen copy, never from what `cargo test` leaves in `target/release`.
- Never `pkill -f`: kill by PID.
- Until the freeze, development contigs only.

## Amendments

**Amendment 1 (2026-09-26, the freeze; written before any held-out R2, R3 or NULL3 run).**

§0 says no r3 binary had run when this file was written. The implementation and every check below happened after
that, all on the three dev contigs. The only held-out data read was check (a), a re-score of the already-published
BASE / R / RQ1 / NULL GTFs. Reports are in the session scratchpad: `rt5_impl.md`, `rt5_prereg.md`, `rt5_review.md`
(an independent review) and `rt6_freeze.md`.

1. **Frozen binaries** (`/mnt/linuxdisk/tmp/rustle_figures/rt3_bin_frozen/`, `SHA1SUMS` verifies 4/4).
   - `copy_assign` **f480847abefefc7ff7997c528c2fd4dce5bda57b**. It replaces §3's placeholder.
   - `as_table` f9556f36; `family_score` 7723029b (both unchanged).
   - `mcl_families` 0c4639e2: its bytes changed, its outputs are identical on one dev GTF. The families step (§11
     item 4) uses `rt_bin_frozen` via `fam_call.sh`, as BASE and R did.
   - Scorer: `readthrough_eval.py` **f6b99dcc1a3f972427cf67bc3de3b846dd5c5f30**, checked before and after these steps.
2. **IV1, byte identity: pass.**
   - The modes unset, off, r, rq1, r2, list:, rALL and r2ALL were run seeded and unseeded (45 runs). 312/312 products
     are `cmp`-identical to cdfad023's, with identical file sets and log lines.
   - `cargo test --release`: 909 passed, 0 failed (902 before, plus 7 new r3 tests).
3. **IV2, port = spec: pass.**
   - (i) r2 and r3 ALL dumps differ only where `B` became `B_guarded`.
   - (ii) An independent pysam recount (reviewer) matches the dump on every tier-B candidate: N, S and U on 676/676
     (chr20), 1,335 (chr16) and 76 (NC_073244.2). No guard decision differs.
   - (iii) `tiers` on the R3 dumps reports 0 mismatches, 0 `N_missing` and 0 `N_lt_S`.
   - (iv) The unit tests of §4 IV2(iv) are in the 7 new tests.
4. **IV3, action and the selected set: pass.**
   - R3's flagged set equals the GKmaj lists exactly: 1,024 / 459 / 106, of which 738 / 418 / 16 are `B_guarded`.
   - R3's GTF, quant, families, assignments and famcn are identical to the list:GKmaj arm's.
   - Transcripts that use a flagged junction: 0 in R3 and 0 in NULL3 (BASE, same lists: 174 / 78 / 44 and
     315 / 137 / 48).
5. **Streaming and genome-wide (not a §4 gate; checked because the held-out driver runs `--genome-wide`).**
   - Streaming = `--materialize-reads` for r3 on chr20 and NC_073244.2, except the O2-only `matched_reads`
     attribute, which also differs for unset and r2.
   - On a slice BAM holding only chr16 and chr20, the driver's `--genome-wide` mode matches per-contig `--region` runs
     for r3, r2 and unset:
     - the flags TSV and ALL dump are byte-identical;
     - the GTF differs only in the run-normalised TPM attribute, which the scorer never reads.
   - With the full-BAM best-AS table, the driver-mode R3 and R2 GTFs equal the dev GTFs and score identically.
6. **Freeze check (a): pass.** Re-scoring the six held-out substrates' BASE / R / RQ1 / NULL with f6b99dcc into
   `rt_arms/tables_check/` reproduces every value of `rt_arms/tables/`. The only differences are 4 BASE `e.status`
   notes: a file path inside "families not measured" changed because BASE is now given via a symlink. The R-prereg
   verdict is identical. **Check (b): pass** (selftest, 14 groups).
7. **Dev arms (dev only; never in the verdict).**
   - NULL3 was drawn from R3's own dumps, seed 20260925. Its list and GTF are byte-identical to rt4's `NULL_GKmaj`.
   - Every score row equals the §2.3 list-arm rows (228 / 228 / 222), so §2.3 now stands as measured on f480847a.
   - R's A1 is not measured on dev: no dev NULL of R exists.

   | contig | fused BASE / R / R2 / **R3** / NULL3 | matched chains BASE / R / R2 / **R3** (lost vs BASE) | R3 A1 (NULL3) | R3 G1-G4 | R3 D vs R | R2 D vs R (descriptive) |
   |---|---|---|---|---|---|---|
   | chr16 | 101 / 77 / 69 / **69** / 98 | 1,605 / 1,609 / 1,595 / **1,612** (6) | −31.7% (−3.0%) | pass | pass | fail (1,595 < 1,606) |
   | chr20 | 60 / 47 / 41 / **41** / 59 | 1,059 / 1,060 / 1,052 / **1,062** (0) | −31.7% (−1.7%) | pass | pass | fail (1,052 < 1,058) |
   | NC_073244.2 | 53 / 39 / 37 / **37** / 50 | 1,595 / 1,588 / 1,597 / **1,596** (0) | −30.2% (−5.7%) | pass | pass | pass |

   - R3 chain precision: .1836 / .2126 / .4130; BASE .1828 / .2098 / .4091; R .1849 / .2134 / .4130.
8. **The held-out stage (§11 items 2-6) now runs unchanged**, via `/mnt/linuxdisk/tmp/rustle_figures/rt3_arms_queue.sh`.

## Outcome (2026-09-26 15:17; R3 and NULL3 on `copy_assign` f480847a, R2 on cdfad023; scorer f6b99dcc; held-out substrates only)

`bench/mechanism/readthrough_eval.py verdict --prereg v3` over `/mnt/linuxdisk/tmp/rustle_figures/rt_arms/tables_v3/`:
- **R3 = EFFECTIVE** (D_pass 6/6). P(s) holds on all six substrates. G5 is measured and passes on human_testis and
  chimp_PTR. D passes on 6/6, with no failure. By §7 R3 becomes the default candidate. **The flip is the user's call.**
- **R2 (descriptive, never judged):** D 5/6. It fails on human_A119b by more than twice its tolerance.
- **Two independent audits** (session scratchpad `figs/rt7_recompute.md`, `figs/rt7_skeptic.md`) found no discrepancy:
  - provenance of binaries, modes and NULL3 lists;
  - every clause recomputed in integer form from the table rows;
  - the §7 mapping;
  - an own recount of FUSED and matched chains (own GFF/GTF code, own gffcompare 0.12.10 runs) on human_testis and
    orangutan_PPY, equal to the tables.
- **The caveats below qualify what "effective" means.** They do not change the verdict.

Substrate order everywhere: human_A119b −chr16/chr20 / human_testis / gorilla_OR6737 −NC_073244.2 / gorilla_KB3781 /
chimp_PTR / orangutan_PPY.

| substrate | FUSED BASE / R / **R3** / NULL3 | A1: R3 reduction (NULL3) | matched chains BASE / R / **R3** (R3 lost / gained vs BASE) | chain precision BASE / R / **R3** | TES genes BASE / R / **R3** | G4.annotated BASE / **R3** | D vs R: R3 − R fused, chains |
|---|---|---|---|---|---|---|---|
| human_A119b | 2,050 / 1,626 / **1,465** / 1,953 | −28.5% (−4.7%) | 35,956 / 35,796 / **35,843** (310 / 197) | .1765 / .1795 / **.1786** | 9,792 / 10,216 / **10,154** | 21,547 / **21,499** | −161, +47: pass |
| human_testis | 200 / 165 / **163** / 199 | −18.5% (−0.5%) | 11,399 / 11,357 / **11,381** (20 / 2) | .4596 / .4616 / **.4609** | 5,069 / 5,112 / **5,103** | 5,709 / **5,695** | −2, +24: pass |
| gorilla_OR6737 | 549 / 454 / **405** / 543 | −26.2% (−1.1%) | 24,277 / 24,180 / **24,273** (37 / 33) | .3496 / .3511 / **.3507** | 6,985 / 7,077 / **7,076** | 13,178 / **13,158** | −49, +93: pass |
| gorilla_KB3781 | 662 / 466 / **437** / 646 | −34.0% (−2.4%) | 25,911 / 25,709 / **25,878** (85 / 52) | .3594 / .3632 / **.3616** | 6,177 / 6,295 / **6,285** | 12,398 / **12,364** | −29, +169: pass |
| chimp_PTR | 420 / 380 / **353** / 417 | −16.0% (−0.7%) | 21,031 / 20,952 / **20,997** (38 / 4) | .3915 / .3927 / **.3925** | 8,378 / 8,419 / **8,430** | 12,469 / **12,448** | −27, +45: pass |
| orangutan_PPY | 687 / 579 / **497** / 678 | −27.7% (−1.3%) | 20,243 / 20,175 / **20,197** (66 / 20) | .1931 / .1942 / **.1940** | 9,053 / 9,221 / **9,199** | 13,750 / **13,718** | −82, +22: pass |

- **A1, G1, G2, G3 and G4.annotated pass on 6/6.** Median A1 reduction .2694; chimp_PTR clears the 10% bar that R
  missed (−16.0% vs R's −9.5%).
- **Tightest judged margins:**
  - G1 on gorilla_OR6737 (.35068 vs .34958);
  - chimp G5 (equal to BASE);
  - D.fused on human_testis (2 loci).
- **G5.**
  - human_testis, Compara Primates (426 families):
    - BASE F / sensitivity / precision: .273973 / .159939 / .954545;
    - R: .275098 / .160701 / .954751;
    - **R3: .275098 / .160701 / .954751** (floors .268973 / .154939 / .949545): pass.
    - Families 338 / 337 / 338; largest 47 in every arm.
  - chimp_PTR, Liftoff pair recall: **.117241 (17/145)** for BASE, R and R3 (floor .112241): pass.
    - Families 379 / 378 / 380; largest 40 in every arm.
- **Not measured, outside §6.3's judged scope:**
  - G4.extra_copy is reported on chimp only, below floor (d): 41/212 in every arm;
  - G5 is not measured on the other four substrates.
  - Neither caps the verdict.

### Guard transfer (R2 vs R3; descriptive, §7.1)
- **D(R2 vs R)** passes on 5/6. It fails on human_A119b by more than twice its tolerance: chains 35,526 against R's
  35,796, below the 2× bar of 35,653. FUSED is 1,461 ≤ 1,626. R2 also fails G2 there (35,526 < 35,597; not a verdict
  input).
- **R3 − R2**, per substrate:
  - chains +317 / +1 / +12 / +12 / +14 / +29, i.e. +0.89% / +0.009% / +0.05% / +0.05% / +0.07% / +0.14% of R2;
  - FUSED +4 / +1 / +1 / +2 / 0 / +1, at most 0.5% of FUSED(BASE).
- **Reading.** Tier B's chain cost reappeared on one held-out substrate: human_A119b, the library that supplied both
  human dev contigs.
  - There the guard turns R2's 270-chain shortfall against R into a 47-chain surplus.
  - On the other five substrates R2 already passes D. By §9 item 2's own wording, the guard's premise is **untested
    there, not confirmed**.
  - The guard costs 0-4 fused loci.
- **Census** (substrate-restricted):
  - `B_guarded` junctions: 21,180 / 35 / 875 / 2,558 / 383 / 2,884, i.e. 81 / 35 / 70 / 85 / 58 / 70% of R2's tier B.
  - Pool alignments they keep (R2 minus R3 "pool alignments removed", from the `[readthrough]` lines): 92,441 / 92 /
    2,519 / 8,341 / 978 / 8,598.
  - Labels RT / ANN / other: 70 / 1,677 / 19,433; 2 / 4 / 29; 7 / 34 / 834; 8 / 52 / 2,498; 3 / 66 / 314; 13 / 130 /
    2,741.
  - ⚠ **The list of the 20 most-read `B_guarded` junctions with N is not reported.**
    - The held-out runs wrote no ALL dump (`rt3_arms_queue.sh` did not set `RUSTLE_READTHROUGH_JUNCTIONS_ALL`).
    - R3's flags dump holds only flagged rows, and R2's dump has no N column.
    - So N of a guarded row is not on disk. §4 item 6's "the ALL dump suffices" was not realised on held-out.

### §8 predictions vs outcome
1. **A1 on 6/6: yes.**
   - Median 26.9% against the predicted 20-26%: just above the range.
   - chimp 16.0% against the predicted 11-15%: above the range, a pass.
   - NULL3 0.5-4.7% against the predicted 1-8%: testis (0.5%) and chimp (0.7%) fall below the range.
2. **G1 on 6/6: yes.**
   - "Possibly below R on human": R3's precision is below R's on 6/6, not only on human.
   - **G3** passes on 6/6, but lies within ±0.5% of R on 5/6 only. A119b is −0.61% (10,154 vs 10,216).
   - **G4.annotated** passes on 6/6, with a loss of 0.15-0.27% (≤ 0.7%): yes.
3. **G2 on 6/6: yes.**
   - R3's loss against BASE is 0.31 / 0.16 / 0.02 / 0.13 / 0.16 / 0.23%. The predicted band was 0.14-0.98%, and 4/6
     fall inside it.
   - OR6737 and KB3781 lose less than the band's floor, which is the favourable direction.
   - KB3781's D pass implies its G2 pass: yes.
4. **D.fused on 6/6: yes.** "Several loci below R" holds on 5/6; testis is only 2 below.
   - **D.chains on 6/6: yes.** R3 − R is +0.13 / +0.21 / +0.38 / +0.66 / +0.21 / +0.11%. That is inside the predicted
     −0.2% to +0.5% on 5/6; KB3781 is above it.
5. **G5 |Δ| < 0.002 on both, largest family unchanged: yes.** testis ΔF vs BASE is +.0011 (0 vs R); chimp Δ is 0.
6. **Guard transfer: failed on human_testis.**
   - (a) D(R2 vs R) was predicted to fail on both human substrates by more than 2×. It did on A119b; on testis it
     passes (162 / 11,380).
   - (b) R3 − R2 was predicted to be ≥ +0.5% of R2 on both human substrates. A119b gives +0.89%; testis gives
     +0.009%.
   - (c) |FUSED(R3) − FUSED(R2)| ≤ 1% of BASE on ≥ 5 of 6: 6/6, yes.
7. **Verdict: EFFECTIVE**, the outcome given ~0.40.

### §9 falsifiers vs outcome (none fires)
1. **Median of (FUSED(R3) − FUSED(R2)) / FUSED(BASE)** = .0019 < .01: not falsified.
2. **D.chains(R2 vs R)** fails on human_A119b only.
   - There R3 − R2 = 317 ≥ ½ · 270 = 135: not falsified.
   - On the other five substrates the guard's premise is **untested, not confirmed**.
3. **`B_guarded`, RT vs ANN:** RT < ANN on 6/6 (70 vs 1,677; 2 vs 4; 7 vs 34; 8 vs 52; 3 vs 66; 13 vs 130). Not
   falsified.
4. **v2 §9 items, with R3 in place of R2:**
   - (1) median of (FUSED(R) − FUSED(R3)) / FUSED(BASE) = .0714 ≥ .03;
   - (2) exempted junctions, RT < ANN on 6/6 (205 vs 543; 12 vs 49; 37 vs 204; 67 vs 334; 24 vs 136; 46 vs 131).
     Recovery of R's lost chains, median over the six = 53.6% ≥ 50% (per substrate 0.3 / 53.5 / 68.4 / 64.1 / 53.7 /
     12.0%). ⚠ This passes by 3.6 points; on A119b R3 recovers 1 of R's 311;
   - (4) D.fused fails nowhere;
   - (5) NULL3 is < half of R3 on 6/6;
   - (6) R3-fused gene pairs absent from BASE: 9 of 1,465 on A119b (0.6%), 0 elsewhere;
   - (7) median reduction 26.9% > 17.4%.
5. **1000 · chains(R3) < 995 · chains(R):** on 0 substrates. R3 is above R on 6/6.

### Caveats (audit 2, `rt7_skeptic.md`; none is a clause, none changes the verdict)
1. **The reduction does not reach the representative, which is the families' exonic input.**
   - `a.rep_fused`, R3 vs BASE: −3.2 / −5.0 / −0.8 / −5.1 / −0.6 / −3.6% (BASE 688 / 100 / 236 / 272 / 180 / 248).
     That is below 10% on 6/6. NULL3 is ≥ half of R3 on 3/6, and R3 is above R on 3/6 (+10 / +1 / +1 on A119b /
     OR6737 / chimp).
   - Only 34-50% of BASE-fused loci have a fused representative. The representative is the most-read transcript,
     usually the host gene's own isoform.
   - Other levels, R3 vs BASE:
     - locus span (the families' aligned sequence): −19.4 / −7.4 / −16.3 / −23.6 / −9.9 / −18.1%;
     - absorbed genes: −16.0 / −7.5 / −15.3 / −21.2 / −8.6 / −17.6% (testis: R3 is 3 worse than R);
     - reads in fused transcripts: −3.4 to −8.7%.
   - G5 equals R to 6 decimals. 333 of 338 (testis) and 370 of 380 (chimp) families are member-identical to R's.
   - **"Effective" is an assembly-level claim, not an O1 (family) claim.**
2. **Chains and precision.**
   - R3 matches fewer reference chains than BASE on 6/6: net −113 / −18 / −4 / −33 / −34 / −46.
   - The precision gain over BASE is all denominator: numerator −0.02 to −0.31%, denominator −0.33 to −1.48%.
   - The transcripts R3 removes match at .037 / .168 / .017 / .062 / .149 / .060, against BASE precision .177 / .460 /
     .350 / .359 / .392 / .193. NULL3's removed transcripts match at .141 / .447 / .287 / .300 / .356 / .151. So the
     removal is targeted.
   - NULL3 also passes G1 on 6/6, so G1 does not discriminate.
   - R3's precision is below R's on 6/6.
3. **Where R3's gain over R comes from: all of it is tier B.**
   - R3 resolves 215 / 12 / 58 / 46 / 33 / 92 R-fused loci; 178 / 11 / 49 / 39 / 28 / 79 of them carry a tier-B
     junction.
   - The fixed ALE exemption re-fuses 52 / 11 / 10 / 19 / 6 / 13 loci, i.e. 24 / 92 / 17 / 41 / 18 / 14% of that gain.
     That is why testis nets −2.
   - The guard gives back 3 / 1 / 1 / 1 / 0 / 1 loci.
   - 23-44% of the labelable tier-B flags are exact annotated introns (549 / 19 / 56 / 51 / 44 / 103), and 52-84% of
     R3's flags cannot be labelled.
4. **The flag set is knife-edge; the verdict is not.**
   - 20.3 / 42.0 / 34.4 / 24.3 / 40.6 / 32.7% of flagged junctions (11-28% of their reads) flip with ±1 read. The
     driver is the exemption's L test (20-47% of tier A). 50-55% of the flagged junctions sit at S = 2.
   - The guard itself is not knife-edge: V1 = N + 1 on ≤ 1.1% of tier B. V1 > 2N would drop 12-25% of tier B.
   - One step stricter: 20 → 21 drops 1.8-4.8% of tier A; 4 → 5 drops 2.0-6.6% of tier B.
   - Clause margins:
     - A1: 380 / 17 / 89 / 158 / 25 / 121 loci;
     - D.chains: +118 / +46 / +141 / +220 / +86 / +62 chains;
     - D.fused: 161 / **2** / 49 / 29 / 27 / 82 loci.
   - Only testis D.fused is fragile. 2 of its 11 tier-B resolutions hinge on knife-edge junctions, and a static proxy of
     U ≥ 2S would cost 3 loci against its 2-locus margin.
   - EFFECTIVE tolerates one D failure; a second would need ≥ 27 loci to move.
   - The unflagged side cannot be measured on held-out: there is no ALL dump.
5. **NULL3 is a random-reads null, not a fused-locus null.**
   - It is matched on alignments (0.91-1.09× R3's), yet it loses 3.1-4.7× R3's chains.
   - Per spliced locus, its fused change is −0.8 to +0.2%, so its raw reduction is locus loss.
   - No locus-matched null was run. A1's specificity half is therefore easy to pass.
6. **Correct readthrough counts as FUSED.** The gene set drops RefSeq's readthrough-described records.
   - 13 of the 35 BASE-fused loci that R3 resolves on testis, and 40 of 575 on A119b, overlap an annotated readthrough
     gene.
   - The ape annotations carry 0-1 such records, so this cannot be measured on 4/6.
7. **G4.annotated is below BASE on 6/6** (−48 / −14 / −20 / −34 / −21 / −32 loci), within the 1% tolerance.
8. **FUSED is the locus exon union (the R prereg §3 text: "exon union (all its transcripts)").**
   - 3.0-14.2% of BASE-fused loci are fused only by the union of different transcripts (224 / 6 / 47 / 94 / 12 / 48).
   - The scorer's metric note reads as per-transcript.
   - Counting only loci with a single fused transcript, R3 gives −22.2 / −16.5 / −21.1 / −25.2 / −14.5 / −24.1%: still
     ≥ 10% on 6/6, with NULL3 ≤ 5%. The scored number is 3-9 points more favourable.
9. **Exposure.**
   - This is the second readthrough verdict on these six substrates.
   - v3 is a re-tune: its guard was selected on the dev contigs after v2's dev D failure (§0 item 1).

### Procedural record (audit 1, `rt7_recompute.md`; no numeric effect)
1. **An out-of-queue held-out score.** `tables_v3/score_gorilla_KB3781.log` (14:16) was run outside the queue, during
   its STOP window (14:14-14:22).
   - R3's KB3781 numbers therefore existed before the families stage and the verdict.
   - The queue's 15:14 score reused that cache. Its key covers the GTF fingerprints and the scorer sha1, so the inputs
     were identical.
   - Every binary and list predates it, and nothing was re-tuned.
2. **Nothing is committed.** This file is untracked, and scorer f6b99dcc is an uncommitted working-tree change.
   Integrity rests on the sha1s recorded in the queue log and Amendment 1: 0bd2765a and f6b99dcc, both re-verified.
3. **G5's `family_score` binary.** G5 was scored with `rt_bin_frozen/family_score` 27aa9445, not the 7723029b that
   Amendment 1 item 1 lists. 27aa9445 is the binary behind BASE's and R's G5 rows and the §6.4 floors, so the
   comparison is like for like. R3's families used `rt_bin_frozen/mcl_families` (6cf8183f), as BASE and R did.
4. **NULL3 re-draws.** The lists were re-drawn byte-identical (list and summary) on testis, chimp, PPY and OR6737.
   A119b's was not re-drawn: its draw takes 627 s, over the 600 s cap. Its header and per-contig counts match R3's dump.
5. **Dev contigs inside held-out substrates.** human_testis and gorilla_KB3781 include chr16 / chr20 and NC_073244.2,
   as §5 specifies (other libraries and individuals).
6. **Dump coordinate convention.** The flags dump writes intron bounds 1-based inclusive, while §1 item 2 states d
   and a 0-based.
   - Under the dump's convention, 0 transcripts of R3's GTF use an R3-flagged junction, on all six.
   - IV2(ii) had already checked N against a pysam recount.

### What is NOT claimed
- **Not a family-level (O1) gain.** G5 equals R to 6 decimals, the representative barely moves, and G5 is measured on
  2 of 6 substrates (not on A119b, where the effect is largest).
- **Not more reference transcripts.** R3 matches fewer chains than BASE on 6/6. It matches more than R.
- **Not that the guard transfers in general.** Its premise was tested on one held-out substrate, the dev library.
- **Not that the constants 20, 4, 2L ≥ S, 3V1 < 5S are robust junction by junction.** They are robust only at the
  clause level. Stricter variants that add flags cannot be evaluated here.
- **Not that every FUSED locus removed was an error** (caveat 6), and not specificity against a locus-matched null
  (caveat 5).
- **Not a default flip.** That is the user's call.
- **Not a licence for a v4 on these substrates.** A v4 chosen after these results needs a new library, or the user's
  explicit acceptance of a third reuse, written into its prereg (§0 item 5).
