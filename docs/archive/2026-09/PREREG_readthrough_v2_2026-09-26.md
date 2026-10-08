# Pre-registration: readthrough junction filter v2 (arm R2): an alternative-last-exon exemption and an own-promoter tier

**Written 2026-09-26, before any held-out number of arm R2 or its NULL exists.** User goal (`/goal`): *"develop
effective readthrough reads filter"*. This file is the successor of `docs/PREREG_readthrough_ends_representatives_
2026-09-25.md` (the "R prereg"; register row 1117). It fixes one new rule (R2) and its matched NULL. It reuses the R
prereg's metrics, substrates, power floors and scorer. It adds a head-to-head clause D against the shipped opt-in arm
R, and it replaces the R prereg's §5 verdicts for this arm only. A default flip is the user's call whatever the outcome.

## 0. What existed and what was seen before this file

- **Seen (development contigs only, design evidence):**
  - `rt2_mech.md`: after arm R, which mechanisms still fuse loci.
  - `rt2_curves.md`: locus-level trade-off curves of ~20 junction rules.
  - `rt2_protect.md`: what R removes wrongly, and whether reads alone can protect it.

  All three are in the session's `scratchpad/figs/`, with scripts under `/mnt/linuxdisk/tmp/rustle_figures_dev/rt2_*`.
  They cover human A119b chr16 and chr20 and gorilla OR6737 NC_073244.2, and nothing else.
- **Seen (held-out, already published):** the R prereg's Outcome table and its addendum, i.e. six aggregate rows of BASE,
  R, RQ1 and NULL, plus the family clause on human_testis and chimp_PTR. To fix the measurement scope of §6.3, this
  author also listed `rt_arms/tables/` and read the G4.extra_copy and G5 rows of those tables for BASE and R. Their
  status: measured, below floor or not measured. The author also read the timestamps of the four `*.families.log` for
  the cost estimate. No junction-, locus- or read-level held-out quantity was read by this author.
- **Code read:** the R prereg's Amendments 1-3, which define S, U, V1, the population, the integer thresholds and the
  removal. Also `bench/mechanism/readthrough_eval.py`: the clause code A1-G5, `null`, and the constants block.
- **Does not exist / not seen:**
  - no R2, tierB, exemption or NULL2 number on any held-out contig;
  - **no assembly of R2 itself on any contig.** The design files measured the parts as separate arms: R; P = R minus
    the exemption (slice setting, see §2.3); and ORq80_50_9524 = R ∪ tierB (full-BAM setting).

  While this file was written, a parallel implementation step (`rt3_impl`) had uncommitted edits in
  `src/bin/copy_assign.rs`, `src/rustle/vg_family/denovo_assemble.rs` and `bench/mechanism/readthrough_eval.py`, the
  last with a `--prereg v2` reading of §7. None of its outputs was read. §4 binds whatever it produces, and §7 below
  binds the scorer, not the reverse. If the scorer's reading differs, the scorer is fixed (see §6.3 on G4.extra_copy).

## 1. The v2 rule (binding, verbatim)

The rule is stated in integer arithmetic. It applies per canonical junction J with S ≥ 2, counted from PRIMARY spliced
reads. The statistics are exactly those the current `RUSTLE_READTHROUGH_JUNCTIONS=r` code computes, plus L:

```
S  = primary reads using J;  U = spliced primaries starting upstream of J's donor and ending inside J's intron (R_up);
V1 = the existing Q1 count (transcripts from >= 3-read start clusters inside J's intron that splice out of their own
     first exon into the downstream exons);
L  = primary reads using J for which J is the LAST splice in transcript orientation (the acceptor exon is the read's
     terminal exon), i.e. j_last = L / S.
tierA (the shipped R) : U >= 20*S
exempt (ALE protection): 2*L >= S AND 3*V1 < 5*S      (J mostly leads into a terminal exon, no independent promoter)
tierB (own promoter)   : U >= S AND V1 >= 4*S          (R >= .5 and Q1 >= .8)
flag(J) = (tierA AND NOT exempt) OR tierB;  the chains (reads) carrying a flagged J are dropped before locus formation,
exactly as the r mode does. Name: RUSTLE_READTHROUGH_JUNCTIONS=r2.
```

**How the text is realised (no freedom left to the implementer):**
1. **Population and statistics.** S, U and V1 are the frozen binary's (`copy_assign` sha1 452454c6), as the R prereg's
   Amendment 1 defines them:
   - item 2: the integer thresholds;
   - item 3: V1's first-donor lookup, as in the script;
   - item 5: per-region population, primary = not secondary and not supplementary, spliced = ≥ 1 intron from the pool's
     own CIGAR parser.

   Scope as in R: the junction is canonical on its strand (GT-AG, GC-AG, AT-AC) and S ≥ 2. A junction outside the scope
   is never flagged.
2. **L** counts the same S reads, i.e. the same (donor, acceptor, read transcript strand) key and the same intron list.
   A read counts toward L when J is its last intron in transcript orientation: the rightmost intron for `+`, the
   leftmost for `−`. So a one-intron read using J counts. 0 ≤ L ≤ S.
3. **Integer equivalences (exact, no rounding):**
   - `2L ≥ S` ⇔ j_last ≥ 1/2;
   - `3V1 < 5S` ⇔ Q1 = V1/(V1+S) < 5/8 = .625, the doc's frozen Q1 cut;
   - `U ≥ S` ⇔ R = U/(U+S) ≥ 1/2;
   - `V1 ≥ 4S` ⇔ Q1 ≥ 4/5;
   - `U ≥ 20S` ⇔ R ≥ 20/21, which the R prereg prints as .9524.
4. **Identity.** tierB implies NOT exempt, because V1 ≥ 4S gives 3V1 ≥ 12S > 5S. So
   flag(J) = (tierA OR tierB) AND NOT exempt. The implementation must satisfy this identity; a unit test checks it.
5. **Action.** This is Amendment 1 item 4 of the R prereg, unchanged. Every alignment the assembler admits (primary or
   seeded good secondary) is removed from pass 1, widening, seeding and polish when its intron chain holds a flagged
   (donor, acceptor). Flags are computed once, from primaries.
6. **Outputs.**
   - `<out>.readthrough_junctions.tsv` has one row per flagged junction, with S, U, V1 and L and the tier(s) that
     flagged it (tierA only, tierB only, or both). The exact column and value names go in Amendment 1.
   - The ALL dump (`RUSTLE_READTHROUGH_JUNCTIONS_ALL=1`) under `r2` also carries L and the exempt bit. Under unset,
     `off`, `r`, `rq1` and `list:` every product stays byte-identical (§4 IV1).
   - The exempted junctions must be recoverable for the census of §7.1 (the ALL dump's exempt bit suffices).

## 2. Why each part (development diagnosis only)

### 2.1 The development contigs and how often they were read
- **Development = human A119b chr16 and chr20, and gorilla OR6737 NC_073244.2. They were read many times:**
  - chr16 chose R's θ;
  - all three were read by the doc's five refinements, by the R prereg's dev arms, by `rt2_mech`, by `rt2_curves`
    (about 25 list-switch arms per contig) and by `rt2_protect` (slice arms, rules P and PX, a 13-signal logistic).
- **v2's constants were picked from those grids on those contigs.** Dev numbers are optimistic by selection.
- **Dev→held-out shrinkage, the reference point.** R's own dev reduction was −21.7% to −26.4%; its held-out median was
  −17.4%, about 0.73 of the dev value.
- **Correct dev setting.** Dev arms must be restrictions of full-BAM runs: `--region <contig>` with the sample's
  genome-wide best-AS table (`rt2_curves` "Setting").
  - Slice BAMs admit secondaries whose real best alignment is off-slice. They give BASE fused 119 / 62 / 60 instead of
    101 / 60 / 53.
  - Primary-read junction statistics are the same in both settings: S and U recounts agree on 2,977 / 1,608 / 293
    junctions.

### 2.2 tierA: the shipped R, unchanged
- R passed every measured clause on 5 of 6 held-out substrates (register 1117). Its held-out results:

  | | range |
  |---|---|
  | fused loci | −9.5% to −29.6% |
  | intron-chain precision | up on 6/6 |
  | TES-recovered genes | up on 6/6 |
  | matched chains lost | ≤ 0.8% |
  | loci lost | ≤ 0.6% |
- On human, θ = 20/21 is the knee of R(θ) (`rt2_curves` answer 1).
  - BASE→R costs 0.2-0.3 reference chains per fused locus removed.
  - Lowering θ to .75 costs 2.3-4.5 chains per extra locus on human, and below .67 it costs 6-9.5 on every contig.
- Lowering θ is therefore not the lever. v2 keeps tierA and changes only its false removals and its reach.

### 2.3 exempt: alternative last exons (rule P of `rt2_protect`)
- **What R removes wrongly** (`rt2_protect` §1).
  - Of R's flags, 48 / 25 / 12 are exact annotated introns of the donor gene (ANN). Against the flags labelled RT, that
    gives junction precision .714 / .719 / .778.
  - They are mainly **alternative last exons** (21 / 16 / 8 annotated ALEs), intronic polyA (5 / 1 / 1) and unannotated
    polyA sites inside the intron (22 / 8 / 3).
  - They cost 15 real matched reference chains on dev, **every one directly**: the lost chain contains the flagged ANN
    junction. Two more lost chains are the curated readthrough CKLF-CMTM1, and losing those is correct.
- **Why reads can separate them** (§2). Two families of read signals live.
  - The acceptor side is independent: Q1 AUC .82 / .72 / .77.
  - J leads into a terminal exon: j_last AUC .68 / .74 / .87.
  - A readthrough enters gene B's multi-exon structure, while an ALE ends in the exon right after J.
  - Continuity, PAS, motif, PWM, strand and end distance are dead signals (`rt2_protect` §2).
- **Rule P = the exemption** (§3). Keep a flagged J when j_last ≥ .5 and Q1 < .625. Both constants are fixed points, not
  fits: .5 is "majority" and .625 is the doc's frozen Q1.
  - The effect is flat over j_last .34-.67 × Q1 cut .3-.8: ANN protected 14-17 / 6-13 / 9, RT protected 1-4.
  - Junction level:

    | | chr16 | chr20 | gorilla |
    |---|---|---|---|
    | ANN false removals, R → P | 48 → 33 | 25 → 14 | 12 → 3 |
    | RT recall, R → P | .504 → .496 | .520 → .488 | .420 → .400 |
    | precision, R → P | .714 → .781 | .719 → .811 | .778 → .930 |
  - Assembly level: P recovers 13 of the 15 real lost chains. It costs 1 fused locus per contig, about 1-1.7 points
    of reduction. ⚠ These assembly numbers are from the **slice setting** (BASE fused 119 / 62 / 60), which §2.1 calls
    wrong for assembly-level numbers.
  - The 2 chains still lost are FAM234A (multi-exon ALE, j_last .13) and CMTR2 (internal promoter, Q1 .99). Reads cannot
    separate them from a readthrough.
- **Why the exemption is not RQ1** (refuted, register 1117). RQ1 *requires* Q1 on every flag, which kills recall (A1
  failed on 5/6). The exemption only *removes* the flags that are terminal AND lack Q1.

### 2.4 tierB: the downstream gene's own promoter (ORq80_50_9524 of `rt2_curves`)
- **The dominant residual after R is class (1)** (`rt2_mech`): canonical, S ≥ 2 junctions that R scored and let through.
  - It is present in 75% / 77% / 56% of R's fused loci, and it is the only mechanism in 45 / 29 / 22 of them.
  - R misses these junctions because they are nowhere near its cut: their R medians are .50 / .57 / .75.
  - An oracle θ-sweep shows no clean cut: θ = .5 flags +4,471 / +2,775 / +407 junctions contig-wide.
- **The sub-class with the downstream gene's own promoter** (`rt2_curves` item 6).
  - Among the junction-bridged residual loci with R in [.5, .9524), 12 / 6 / 2 have Q1 ≥ .8 on every bridging
    transcript: gene B has its own start cluster and its own first exon, e.g. DN_chr16_21842574_5 (R .79, Q1 .98).
  - tierB targets exactly this sub-class: R ≥ .5 and Q1 ≥ .8. The .5 and .8 come from a grid of 5 Q1-gated variants
    and 7 plain OR variants, chosen on dev as the best human rule inside G2's 1% budget.
- **Evidence** (`rt2_curves` items 2-5).

  | ORq80_50_9524 vs R | chr16 | chr20 | gorilla |
  |---|---|---|---|
  | fused loci | −32.7% vs −23.8% | −33.3% vs −21.7% | −32.1% vs −26.4% |
  | net matched chains | −14 / 1,605 | −9 / 1,059 | −5 / 1,595 |
  | junction recall | .735 vs .504 | .659 vs .520 | .580 vs .420 |
  | junction FPR | .0115 vs .0044 | .0106 vs .0040 | .0020 vs .0012 |
  | margin over R's θ-curve at equal chain loss | +2.8 / +2.1 loci | +5.0 / +5.0 | +2.0 / +6.3 |

  - It is ROC-better than the θ-sweep at the junction level on all three contigs.
  - A matched NULL moves fused loci by −9.9% / −5.0% / −3.8% for the similar OR50_9524, and costs 2.2-3.5× more chains.
  - ⚠ `rt2_curves` itself says the margin over R's curve (0 to +5 loci per contig on 53-101 fused loci) is **within
    dev noise**. That is why this file exists.
- **tierB's cost is human-heavy.** Chains lost vs BASE:

  | | chr16 | chr20 | gorilla |
  |---|---|---|---|
  | R | 8 | 2 | 7 |
  | ORq80_50_9524 | 27 | 12 | 7 |

  The dev files did not examine why tierB's cost falls on human. A likely reason is that the curated human annotation
  holds more internal-promoter isoforms.

### 2.5 What v2 does not target (by design; the ceiling of any junction rule)
After R, some residual fused loci are out of reach of every junction rule (`rt2_curves` item 6, `rt2_mech`):
- 17 / 11 / 17 are exon bridges or grouping.
- 21 / 10 / 3 are dominant fusions, with R < .2 on every bridging transcript. They are kept on purpose (e.g.
  PKD1P6-NPIPP1, SLX1B-SULT1A4).
- The ghost-link regroup (class 6, 3 / 2 / 2 loci) is a separate lever and is not part of v2.

### 2.6 Dev projection of R2 (arithmetic on the separate arms; NOT a measurement, NOT in any verdict)
The projection is additive: ORq80_50_9524 (full-BAM setting) plus P's changes relative to R (slice setting). It assumes
tierB and the exemption act on disjoint junctions, which holds because tierB needs Q1 ≥ .8 and the exemption needs
Q1 < .625. It ignores interactions within a chain.

| dev contig | fused BASE → R → R2 (reduction) | matched chains BASE → R → R2 | R2 vs BASE | R2 vs R |
|---|---|---|---|---|
| chr16 | 101 → 77 → ~69 (−31.7%) | 1,605 → 1,609 → ~1,595 | −0.62% | −0.87% |
| chr20 | 60 → 47 → ~41 (−31.7%) | 1,059 → 1,060 → ~1,052 | −0.66% | −0.75% |
| NC_073244.2 | 53 → 39 → ~37 (−30.2%) | 1,595 → 1,588 → ~1,597 | +0.13% | +0.57% |

⚠ Read against D (§6.2), this projection puts the chains part of D **beyond twice its tolerance on both human dev
contigs**. The exemption's recovered chains (+4 / +2) do not pay for tierB's human cost (−18 / −10 net vs R). This is
recorded before any held-out number. It moves neither the rule nor D; §8 predicts from it. R2 itself must still be
assembled on dev (§4) and reported as dev.

## 3. Arms

| arm | what | status |
|---|---|---|
| **BASE** | current default; `rt_arms/<sample>/<sample>.BASE.gtf` (= `runs/<id>/<id>.gtf`), restricted | **reused**: the held-out arm GTFs and families of the R prereg; no re-assembly |
| **R** | `RUSTLE_READTHROUGH_JUNCTIONS=r`, frozen `copy_assign` 452454c6 | **reused**: `rt_arms/<sample>/<sample>.R.gtf` and, on human_testis and chimp_PTR, its families |
| **R2** | `RUSTLE_READTHROUGH_JUNCTIONS=r2` (§1) | new, genome-wide on the six samples |
| **NULL2** | matched random removal to R2's flagged reads | new; the scorer's null (§3.1) |

The reuse of BASE and R is valid only if IV1 (§4) holds: the new binary reproduces them byte for byte where it is
checked. If the implementation changes any code path shared with unset, `off`, `r` or `list:` (anything beyond a new
branch taken only under `r2`), BASE and R are re-assembled on all six samples with the new binary and the reused arms
are retired. That would be recorded as an amendment before any R2 number is scored. RQ1 and the old NULL are not arms
here; R's NULL is not reused for R2.

### 3.1 NULL2
- **Draw.** `readthrough_eval.py null --flags <R2's <out>.readthrough_junctions.tsv>` per sample, with the default seed
  20260925 and one RNG per `sample:contig`. This is the R prereg's Amendment 2 procedure unchanged, with F = R2's
  flagged set.
- **Candidates** are the canonical S ≥ 2 junctions that R2 does not flag. That includes R's exempted junctions. They
  are R-pass readthrough-like junctions, so drawing them can only make NULL2 stronger. This is conservative against R2
  and accepted.
- **Target.** The draw targets the distinct primaries that carry an R2-flagged junction, matched by floor(log2 S).
- **Arm.** NULL2 runs as `list:<NULL2 list>` with R2's binary. It gets metrics (a)-(d), as NULL did.
- **A1's null part** for R2 is judged against NULL2 only.

**Not run on held-out (by design).** There are no component arms on the held-out substrates: no "R minus exemption"
and no "R ∪ tierB". The held-out stage answers only "R2 vs BASE" and "R2 vs R". It must not become a design tool for a
v3. The attribution to parts comes only from the flag census (§7.1).

## 4. Implementation gates (dev only; all before any held-out R2 run)

Any failure is fixed in the code, never in the rule.
- **IV1, byte identity.**
  - Unset, `off`, `r`, `rq1` and `list:` of R's dump: every product that the R prereg's Amendment 3 item 4 lists is
    `cmp`-identical to the frozen 452454c6 outputs.
  - The products are compared in the full-BAM setting (`--region` with the genome-wide best-AS table, §2.1) on the
    three dev contigs, and on the dev slices already used by Amendment 3.
  - `cargo test --release` passes.
- **IV2, port = spec.** On the three dev contigs:
  - (i) The dump's S, U and V1 equal the frozen binary's ALL dump on every row.
  - (ii) L equals an independent pysam recount (the `rt2_protect` j_last instrument) on every junction. The only
    allowed differences are rows explained by the R prereg's Amendment 1 item 5 (`N I N` parse, leading or trailing
    `N`), and each one is listed.
  - (iii) Applying §1's integer rule in Python to the ALL dump's S, U, V1 and L reproduces the flagged set exactly:
    0 mismatches.
  - (iv) Unit tests cover:
    - the identity of §1 item 4;
    - the boundaries `2L = S`, `3V1 = 5S − 1`, `3V1 = 5S`, `U = S`, `V1 = 4S` and `U = 20S`;
    - a one-intron read counting toward L;
    - `−`-strand orientation.
- **IV3, action.**
  - 0 transcripts of R2's GTF (and of NULL2's) use a flagged (listed) junction.
  - `list:` of R2's own dump reproduces R2's GTF, quant, families, assignments and famcn byte for byte.
  - The `[readthrough]` log line names the arm, and the run cache is keyed on the switch, so no arm is served another
    arm's cached products.
- **Dev arms (reported, never in the verdict).**
  - R2 and NULL2 on the three dev contigs in the full-BAM setting, scored with `readthrough_eval.py score` (not
    `--heldout`), with clause D against R.
  - The table goes into Amendment 1 beside §2.6's projection.
  - Dev results cannot change the rule, the NULL, a clause or a tolerance. The held-out stage runs whatever dev shows,
    unless the user stops the work. A stop would be recorded, and nothing would be re-tuned.
- **Freeze.** Amendment 1 records the following before any held-out R2 run:
  - the `copy_assign` sha1 serving R2 and NULL2;
  - the scorer sha1;
  - that the scorer's A1-G5 clause code is unchanged from sha1 0a50ab11, checked by (a) re-scoring the reused BASE and R
    GTFs of the six held-out substrates, which must reproduce every recorded row of `rt_arms/tables/`, and (b) its
    selftest;
  - that the scorer's v2 verdict implements §7 exactly.

  Check (a) re-scores arms whose numbers are already published, and it is part of the freeze. A mismatch stops the
  work until explained.

## 5. Substrates

The **six held-out verdict substrates are the R prereg's**, each judged on its own; species are never pooled:

| id | substrate |
|---|---|
| V1 | human_A119b genome minus chr16, chr20 |
| V2 | human_testis whole genome |
| V3 | gorilla_OR6737 genome minus NC_073244.2 |
| V4 | gorilla_KB3781 whole genome |
| V5 | chimp_PTR whole genome |
| V6 | orangutan_PPY whole genome |

⚠ **Honest note on exposure.**
- **They are no longer untouched.** These six substrates served R's verdict (Outcome, 09-25) and its family addendum
  (09-26).
- **R's chimp miss motivated v2.** On chimp_PTR, A1 was −9.5% against a 10% bar, the only failed clause.
- **v2 was designed on the three development contigs only.** No candidate rule, statistic or locus was computed on or
  read from any held-out contig. What the designers knew of the held-out substrates is the six aggregate rows and the
  family addendum, not where R failed within them.
- **The residual risk is selection on a known failure.** A rule built to add reduction exists because one held-out
  substrate fell short. The bar therefore does not weaken any clause of R's. It adds D, and it keeps "all six" for
  EFFECTIVE.
- **After this prereg, these six substrates have served two readthrough verdicts.** A v3 rule chosen after seeing R2's
  held-out results needs a new library, or the user's explicit acceptance of a third reuse. That acceptance would be
  written into that v3 prereg.

**Power floors** are the R prereg's §4: (a) BASE FUSED ≥ 50; (b) universe ≥ 1,000 genes; (c) BASE matched chains ≥
500; (d) ≥ 1,000 reference loci; (e) ≥ 30 families or pairs. They are properties of the reused BASE. In R's held-out
tables A1, G1, G2, G3 and G4.annotated qualify on all six, and G5 qualifies on human_testis and chimp_PTR. So no
verdict clause of §7 is below its floor. G4.extra_copy is the exception (§6.3).

## 6. Metrics and clauses

### 6.1 The R prereg's clauses, unchanged (metrics (a)-(e) of its §3, same scorer code, same annotations, gffcompare 0.12.10)

| clause | measure | R2 passes iff |
|---|---|---|
| **A1** efficacy | FUSED loci (`a.fused`) | reduction vs BASE ≥ 10% **and** NULL2's reduction < half of R2's |
| **G1** precision | intron-chain precision (raw counts) | R2 ≥ BASE |
| **G2** sensitivity | matching reference intron chains (`c.matching_intron_chains`) | 100·R2 ≥ 99·BASE |
| **G3** ends | TES-recovered genes (fixed universe) | R2 ≥ BASE |
| **G4** loci kept | found read-supported reference loci, annotated (and extra copies, §6.3) | 100·R2 ≥ 99·BASE |
| **G5** families | human: Compara Primates bipartite F, sensitivity, precision; apes: Liftoff pair recall | each ≥ BASE − 0.005 |

The reasons for X = 10%, Y = 1% and 0.005 are the R prereg's §5. They are not re-argued, and they are not moved.

### 6.2 Clause D: R2 head to head with R (new)

On every one of the six substrates:
- **D.fused:** FUSED(R2) ≤ FUSED(R). The tolerance is 0 loci.
- **D.chains:** 1000 · matched chains(R2) ≥ 998 · matched chains(R). That is 0.2% of R's count, on G2's measure.
- **D passes on a substrate iff both parts pass.**
- **A D failure is "larger than twice its tolerance"** iff FUSED(R2) > FUSED(R), or 1000 · chains(R2) < 996 · chains(R).
  - This is the literal reading of a zero tolerance: any D.fused failure counts as larger than twice it.
  - So the one D failure that EFFECTIVE tolerates can only be a D.chains shortfall in (0.2%, 0.4%].
  - No fused tolerance was invented.

**Why D exists.** A1-G5 ask whether R2 beats BASE. Adopting R2 in place of the existing opt-in R also needs R2 to beat
R. Two failure modes would otherwise pass everything:
- **R2 is only a move down R's own θ-curve**, buying loci with chains, which `rt2_curves` shows costs 3-4.5 chains per
  locus on human.
- **The exemption costs more loci than tierB gains.**

**Why D.fused has no tolerance.** The design's claim (§2.3, §2.4) is that tierB's gain (+6 to +11 points on dev) exceeds
the exemption's cost (about 1 locus per contig). A substrate where R2 keeps more fused loci than R contradicts that
claim directly.

**Why 0.2% for D.chains.**
- (i) **It is G2's headroom above R on R's costliest held-out substrate.** On gorilla_KB3781, R lost 0.78% of BASE's
  matched chains and G2 allows 1.0%. An R2 that spends more than about 0.2% beyond R there is already at G2's edge. The
  two thresholds are 25,658 (D) and 25,652 (G2).
- (ii) **It is about half of R's own median held-out chain cost.** That median is 0.39%, range 0.34-0.78%. To be called
  better than R, R2 may spend at most half again what R spent.
- (iii) **It is one fifth of G2's Y**, so D cannot be met by a rule that uses most of G2's budget.

⚠ §2.6 projects D.chains at −0.87% and −0.75% on the two human dev contigs. D is therefore expected to bind on human
(§8). This was known when 0.2% was written down, and it is kept.

### 6.3 Pre-committed measurement scope for G4 extra copies and G5

- **G5 is measured on human_testis and chimp_PTR only.**
  - human_testis uses Compara families at Primates (426 families): `family_score` bipartite F, sensitivity and precision.
  - chimp_PTR uses Liftoff (record, extra copy) pairs: pair recall.
  - Both use the R prereg addendum's recipe: `tools/rustle_pipeline.sh families` on R2's GTF, with the same
    `mcl_families` binaries and shard wrapper. BASE's and R's families there are reused.
- **G4 extra copies can only be measured on chimp_PTR.** Only chimp has a merged Liftoff extra-copy table; R's tables
  say "Liftoff table not merged" on the other five. On chimp the extra-copy universe is below power floor (d): BASE
  found 41. So G4.extra_copy is **reported, not judged**, on chimp, and not measured elsewhere. G4.annotated is measured
  and judged on all six.
- **"Not measured" elsewhere does NOT cap the verdict.** This departs on purpose from the R prereg's §5/§8, where a
  clause that was not measured capped the verdict at keep opt-in. Two reasons:
  - **Cost.** No BASE or R families exist on the other four samples. Measuring G5 there would need 12 genome-wide
    families runs (BASE, R and R2 × 4). The addendum's runs took about 20 min (human_testis) to about 70 min (chimp)
    each, and the human_A119b assembly is the largest. That is roughly 4-14 h of heavy, serial machine time for a
    clause that the R result shows barely moves.
  - **The R result.** Under R the family measures moved by < 0.002: F +.0011, sensitivity +.0008, precision +.0003,
    pair recall 0. That is at most a fifth of G5's 0.005 tolerance, with no hub (largest family unchanged). R2 changes
    the read set by a similar order (tierB adds, the exemption removes). A G5 failure on an unmeasured substrate would
    need a family movement about 5× R's, in the harmful direction.
- **Binding on the scorer.** Its v2 verdict must not cap on G4.extra_copy or on G5 outside human_testis and chimp_PTR.
  Any other A1-G4.annotated clause that is "not measured" (for example a missing NULL2) is not a pass.

### 6.4 The numbers the clauses compare against, fixed now (BASE and R are reused, so these are known)

| substrate | A1: FUSED(R2) ≤ | D.fused: ≤ R | G2: chains ≥ | D.chains ≥ (2× tol ≥) | G1: precision ≥ | G3: TES ≥ | G4.annotated ≥ |
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
- **The table's source.** It is derived from the R prereg's Outcome table and the recorded clause rows. If the freeze
  re-score (§4) does not reproduce them, the work stops.
- **Chimp's A1 needs FUSED(R2) ≤ 378**, i.e. two fewer fused loci than R's 380.

## 7. Verdicts (arm R2)

Write P(s) for "every judged clause passes on substrate s". The judged clauses are A1 (with NULL2), G1, G2, G3 and
G4.annotated, plus G5 on human_testis and chimp_PTR.

- **EFFECTIVE: adopt R2 as the default candidate. The flip is the user's call.** This needs all three:
  - P(s) holds on all six substrates;
  - G5 is measured and passes on both human_testis (F, sensitivity, precision) and chimp_PTR;
  - D passes on ≥ 5 of 6, and no D failure is larger than twice its tolerance (§6.2).

  A flip changes the O1 default loci, so every downstream family and copy-assignment figure must be re-run.
- **KEEP OPT-IN:** R2 is not EFFECTIVE, and every A1-G5 failure falls on at most one substrate.
  - This covers the case "A1-G5 pass everywhere but D fails".
  - It also covers a judged clause that is not measured on a substrate: that is not a pass, so it blocks EFFECTIVE, but
    it is not counted as a failure.
  - Then R2 and R both stay opt-in, and the recorded numbers are what a user who opts in chooses between.
- **REFUTE:** A1-G5 clauses fail on ≥ 2 of the six substrates.
  - This gives a register row. The switch value stays opt-in until a cleanup wave removes it with the user's agreement.
- **A substrate missing from the tables** makes R2 "undecided", unless ≥ 2 present substrates already fail, which is
  REFUTE.
- **Never pooled.** Species are never pooled, and dev never enters the verdict.

### 7.1 Reported beside the verdict (not clauses)
- **Flag census per substrate:**
  - flagged junctions by tier: tierA only, tierB only, both;
  - exempted junctions;
  - removed alignments;
  - each of these by annotation label: RT, ANN, novel classes. Labels only, as `rt2_protect` §1.
- **The 20 most-read tierB-only flagged junctions**, with their gene pairs.
- **Chain turnover:**
  - reference chains R lost vs BASE that R2 matches ("recovered");
  - chains R2 lost that R did not.
- **Fusions R2 creates:** R2-fused loci whose annotated gene pair is not fused in BASE. On dev, R created none.
- **Loci and transcripts per arm;** loci lost and gained vs BASE and vs R.
- **Largest family** on human_testis and chimp_PTR (hub check).
- **NPIP on human_testis:** the Compara Primates family holding NPIP-named genes, with members found and bipartite F
  for BASE, R and R2. This is descriptive. The R prereg's NPIP cap needed both human samples, and families are measured
  here on human_testis only, so it cannot be evaluated as written. It does not change this verdict, and the user sees
  it before any flip.
- **Cost:** wall time and peak RSS per arm vs BASE.

## 8. Predictions (before any R2 number, dev or held-out)

1. **A1.** R2 passes on 6/6. Median reduction is 20-26%: R's held-out 17.4% plus tierB's dev gain net of the exemption,
   shrunk by R's dev→held-out ratio of about 0.73.
   - **chimp_PTR:** 11-15%, a pass, but the least certain. If chimp's residual is mostly exon bridges or dominant fusions
     (§2.5), it misses again.
   - **NULL2:** 1-8%.
2. **G1** passes on 6/6: above BASE, possibly below R on human where the exemption returns unmatched chains.
   **G3** passes on 6/6, within ±0.5% of R. **G4.annotated** passes on 6/6, with a loss ≤ 0.7%.
3. **G2 is the clause most at risk, on human.**
   - R already spent 0.37% (testis) and 0.45% (A119b) held-out. Adding tierB's dev human cost net of the exemption
     (about 0.75-0.9% vs R) projects −1.1% to −1.3% vs BASE: a **G2 failure on one or both human substrates** if dev
     transfers.
   - On the four ape substrates, whose dev analogue is gorilla, G2 passes (−0.2% to −0.8%).
4. **D.** D.fused passes on 6/6. D.chains fails on both human substrates by more than twice its tolerance, and passes
   on the four ape substrates, where the exemption's recoveries outweigh tierB's cost.
5. **G5:** |Δ| < 0.002 on both measured substrates; a pass. The largest family is unchanged.
6. **Verdict:**

   | verdict | probability | what it takes |
   |---|---|---|
   | KEEP OPT-IN | ~0.50 | D fails on human, G2 holds on at least one human substrate |
   | REFUTE | ~0.35 | G2 fails on both human substrates |
   | EFFECTIVE | ~0.15 | |

   A result better than this is evidence that the human dev chain cost was dev-specific (chr16 is the NPIP/SMG1P
   contig).

## 9. What would falsify the design reasoning (reported whatever the verdict)

Each item names the design claim and the held-out observation that contradicts it. They are read from the §6-§7
measures and the census; no extra held-out arm is run.

1. **"tierB's locus gain transfers."** Falsified if the median over the six of (FUSED(R) − FUSED(R2)) / FUSED(BASE) is
   < 0.03. Dev implies +5 to +10 points net. Below 3 points, the own-promoter class is a dev artefact, or it is rare
   outside chr16, chr20 and NC_073244.2.
2. **"The exemption protects alternative last exons, not readthroughs."**
   - Falsified if, among R-flagged junctions that R2 exempts, those labelled RT are at least as many as those labelled
     ANN on ≥ 2 substrates. On dev it was ANN 14-17 / 6-13 / 9 against RT 1-4.
   - Also falsified if R2 recovers fewer than half of R's lost reference chains, taken as the median over the six. On
     dev it was 13 of the 17 lost chains (15 real losses plus the 2 CKLF-CMTM1).
3. **"tierB's chain cost is human-specific"** (§2.4: 0 extra losses on gorilla). Falsified if R2 − R matched chains is
   ≤ −0.5% on ≥ 2 of the four ape substrates. The claim ties the cost to the curated human annotation. If the apes
   (Gnomon annotations) pay it too, the cost is a property of the rule.
4. **"The exemption's fused cost is small"** (about 1 locus per contig). Falsified if D.fused fails anywhere, i.e.
   FUSED(R2) > FUSED(R).
5. **"The gains are specific to the flagged junctions."** Falsified if NULL2's reduction is ≥ half of R2's on any
   substrate, which is A1's null part.
6. **"R2, like R, removes only existing fusions."** Falsified if R2-fused loci whose gene pair is not fused in BASE are
   > 5% of R2's fused loci on any substrate.
7. **"R2 improves on R at all."** Falsified if R2's median reduction is ≤ R's held-out median (17.4%).
8. **Not falsifiers.**
   - A D.chains failure alone falsifies only "R2 dominates R". It says nothing about the direction of the design.
   - A1 failing on chimp alone, with items 1-7 holding, says that chimp's residual is outside every junction rule's
     reach (§2.5). A readthrough rule is then the wrong lever for it.

## 10. Order, stop rules, machine rules

1. **Implement, then run IV1-IV3.** Implement behind the switch, then run IV1-IV3 and the dev arms. Amendment 1 records
   the sha1s, the IV results and the dev table.
2. **Held-out assemblies.** Run R2 genome-wide on the six samples (human_A119b, human_testis, gorilla_OR6737,
   gorilla_KB3781, chimp_PTR, orangutan_PPY), writing to `/mnt/linuxdisk/tmp/rustle_figures/rt_arms/<sample>/
   <sample>.R2.*`.
3. **NULL2.** Draw NULL2 from each sample's own R2 dump, then run `NULL2` = `list:`.
4. **Families.** Run the families stage on R2 for human_testis and chimp_PTR.
5. **Score.** `readthrough_eval.py score --heldout` on the six substrates, with BASE, R, R2 and NULL2 and `--null
   R2=NULL2 --versus R`. The new tables go to a new directory (`rt_arms/tables_v2/`); R's tables are never
   overwritten.
6. **Verdict.** `readthrough_eval.py verdict --prereg v2`. The Outcome is appended here, and every number becomes a new
   register row.

**No change to the code, the rule, the NULL or a clause after step 1's freeze.** The only exception is a bug fix that
re-runs every R2 and NULL2 arm and is recorded as an amendment. A step that cannot be computed makes its clause "not
measured". That caps the verdict (§6.3, §7) unless the clause is outside §6.3's scope.

**Machine rules:**
- One heavy process at a time, in the foreground: `flock -w 900 /mnt/linuxdisk/tmp/rustle_heavy.lock timeout 600
  <cmd>`, split per sample where needed.
- Outputs and `TMPDIR` go under `/mnt/linuxdisk`.
- Build only into `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target_dev`.
- Take the binary from the frozen copy, never from what `cargo test` leaves in `target/release` (R prereg Amendment 3,
  item 5).

## Amendments

**Amendment 1 (2026-09-26, before any held-out R2 run): superseded; the held-out stage of this file is never run.**
- R2 was implemented and frozen (`copy_assign` cdfad023; `rt3_impl.md`). On the dev arms, clause D failed on both
  human dev contigs on chains alone: chr16 1,595 vs R 1,609 (bar 1,606) and chr20 1,052 vs 1,060 (bar 1,058). Gorilla
  NC_073244.2 passed.
- §4 says dev results cannot change the rule. The rule was nevertheless changed: a guard on tier B was selected on the
  same dev contigs (`rt4_tierb_guard.md`). That re-tune is recorded as a new pre-registration,
  `docs/archive/2026-09/PREREG_readthrough_v3_2026-09-26.md` (arm R3), and not as an amendment of this file.
- R2 runs on the held-out substrates only as v3's descriptive "guard transfer" arm. None of its numbers is a v2
  outcome.
