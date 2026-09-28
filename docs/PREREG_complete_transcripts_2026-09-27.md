# Pre-registration: compatible-containment collapse (arm A), fewer partial transcripts from `--assemble-only`

**Written 2026-09-27, before any held-out or SIRV number of arm A or of its NULLs exists.** This is the binding text.
It replaces the session draft `PREREG_complete_transcripts_DRAFT.md` after the hostile review (`ct_critique.md`,
verdict NOT READY) and the build and freeze (`ct2_build.md`). All three reports are in the session scratchpad
(`…/scratchpad/figs/`). The frozen scripts are in `/mnt/linuxdisk/tmp/rustle_figures/ct_frozen/` (Amendment 1).
A default flip is the user's call, whatever the outcome.

**User goal** (`/goal`, 2026-09-27): *"ensure the assembler part of the pipeline tries as best as possible to emit
complete transcripts, meaning that when running gffcompare afterwards there are few partial transcripts like m
categories"*.

**User decisions (2026-09-27; also the acceptance line of Amendment 1).**
1. **The target is gffcompare class `c`** (partial fragments), **not `m`**.
   - We already emit the fewest `m` of the four assemblers: 2.3 / 1.7% of multi-exon queries, against 3.9-6.9% for
     StringTie, FLAIR and IsoSeq on the same dev contigs.
   - The `m` rule space is closed: r1075 was adopted and r1080 refuted.
   - `k`, `m` and `n` stay inside C1's partial share, so any change to them is measured.
2. **Held-out = the SIRV E0 spike-ins of human_testis** (a fresh truth) **plus the six genome samples minus their dev
   contigs.** This is recorded as the **fourth reuse** of the six samples: the R, RQ1 and R3 verdicts came before
   (r1117, r1119).
3. **chr21 and chr22 are dropped from human_A119b's held-out (V1)** because they were the dev set of the 09-23 polish
   rule search (r1062-r1085), and r1077's attribution is the premise of §2.1 (critique H3). They are scored and
   reported but never judged.

**What A does and does not do** (reworded per critique B2).
- **What it removes.** A removes every emitted transcript that is an **end-compatible contiguous sub-chain of another
  emitted transcript** at the same locus and strand.
- **What survives in its place.** On dev, for only 14-22% of the removed `c` is the surviving container the
  reference's complete form (`=` to the same reference). The other containers are:
  - another `c`;
  - a `j`;
  - a `k`, `m` or `n`;
  - an `=` to another reference.
- **What it does not do.** A completes nothing. It keeps every locus and every junction (identities (i) and (ii) of §1).
  Completion itself is at the read ceiling (r1077; §2.1).
- **Consequences.**
  - A partial-share reduction can be bought by trading an intron-correct `c` for a `j` container, which lies outside
    the partial share (c+k+m+n). Clause C5 (§6.3) and the `j`-drop rows (§7.1) exist for that reason.
  - Removing any set whose `=` rate is below base precision raises precision, so C3 is not evidence of selectivity.
    C3′ and the NULL_S conjunct of C1 are that evidence.

---

## 0. What was seen before this file, and the procedural record

- **Development contigs, design evidence.** Human A119b chr20 and chr16, and gorilla OR6737 NC_073244.2: full BAM,
  `--region`, seeded, the driver's flags, frozen `copy_assign` f480847a, dedup count mode.
  - `ct_attr.md`: where `c`, `k`, `m` and `n` come from.
  - `ct_levers.md`: about 40 lever variants. **Arm A (`L1e_compat10_cont_ge0.5`) was selected on these three
    contigs.**
  - The draft author's checks (`ct_prereg/`): the old NULL (5 seeds), the invariants and the SQANTI3 projection.
  - `ct_critique.md` (its scripts are in `ct_critic/`).
  - `ct2_build.md`, `ct2_build/summary.md`, `grid_table.md` and `exposure_dev.txt`:
    - NULL_S over 5 seeds;
    - **every clause in its critique-literal form and in the form fixed below;**
    - the TOL × ρ grid;
    - the truth-side exposure of the dev truths.
- **Truth side only (no read, no assembler output).**
  - The SIRV E0 truth (`ct_sirv`: 69 transcripts, 61 multi-exon, 9 true ISMs), its BAM, and the three tools' SIRV
    GTFs, re-scored by the frozen scorer. The re-score reproduces `ct_sirv`'s tool table exactly.
  - A's structural conditions (1)-(3) applied to the dev truths and to the SIRV truth.
- **Held-out, already published.** The BASE rows (matched chains and chain precision per substrate) of
  `PREREG_readthrough_v3_2026-09-26.md` Outcome, carried over from the draft (§2.9, §6.5). The held-out BASE/R3 GTFs and
  annotations were only hashed or stream-restricted to compute sha1s (`INPUTS_SHA1.tsv`, `REF_SHA1.tsv`).
- **Exists, and not read by this author or the build.** BASE's gffcompare class distribution on the held-out
  substrates:
  - `rustle_figures/assembly/<sample>/eval_all/`;
  - `rt_arms/work/score/<sample>/*/gffcompare/`;
  - `rt_arms/tables_v3` class rows;
  - publication figures 1-2 (OR6737 genome-wide, human chr20-22).

  They are properties of BASE, and the user and earlier sessions have seen them.
- **Does not exist.** Any held-out output of A, NULL_S or NULL, and any SIRV assembly of ours, other than the unscored
  380-read smoke test.
  - The smoke outputs were moved **unread** to `ct_frozen/sealed/ct_sirv_smoke/` (files mode 000, directory mode 500;
    Amendment 1).
  - A pointer file `ct_sirv/smoke/SEALED.txt` remains in the old location. Nobody opens the sealed directory.
- **Procedural record (plain).**
  1. **A was selected on the three dev contigs.** Its constants are plateau picks: TOL is flat over 5-50 bp and ρ
     over 0.25-0.5. TOL was **not** derived from an annotation-free measure (critique M1). Every dev number here is
     the optimistic edge.
  2. **The NULLs were designed after the lever.** NULL (draft §3.1) and NULL_S (critique B1) had their pools fixed
     before scoring and never changed after.
  3. **Decisions D1 and D2 (§6.2, §6.3) were made by this author after seeing the dev clause results of
     `ct2_build`,** which the build asked the author to decide.
     - Under the critique-literal clauses, A fails on dev:
       - C1's "NULL_S < ½ of A" conjunct fails 14/15;
       - C3′ "4·e_A ≤ e_NS" fails 6/15;
       - C5 fails 14/15;
       - H_NS fails 15/15.
     - **D1** makes every comparison against NULL_S **direction-only** (no constant). The reason is structural (§6.2):
       NULL_S's pool contains D by construction, so any constant factor against it measures the pool's size as
       much as condition (3). The ½ and ¼ constants are kept where they were designed, against a comparator
       disjoint from D (the old NULL, reported, and falsifier F1).
     - **D2** keeps C5 exactly as the critique wrote it (A ≤ NULL_S; it already has no constant). It makes C5 decide
       the **default-candidate tier** in place of H, which critique B1 demotes. C5 is not part of P(s).
     - **Both decisions are informed by dev numbers, and this is recorded as such.** D1 turns two dev failures into
       passes. D2 keeps a clause that fails on dev and makes it bind the default flip, so it does not rescue the
       lever's default candidacy.
     - A change to D1 or D2 is allowed only before any held-out or SIRV number of A, NULL_S or NULL exists, as a
       dated amendment.
  4. **What keeps the held-out stage a test.**
     - No held-out or SIRV number of A, NULL_S or NULL exists.
     - The rule has no constant beyond TOL and ρ, both frozen.
     - Every clause tolerance comes from an earlier prereg (§6), or the clause is direction-only.
     - The scorer is frozen and refuses held-out or SIRV scoring unless this file has an "Amendment 1" heading and its
       sha1.
  5. **Exposure** (§5). The six genome substrates carry three earlier verdicts (R, RQ1, R3). The user accepted their
     fourth reuse on 2026-09-27 (decision 2; Amendment 1). SIRV is fresh.

### 0.1 Review items and their fixes (`ct_critique.md`; build questions of `ct2_build.md`)

| item | problem | fix in this file | where |
|---|---|---|---|
| **B1** | The NULL was a random-drop sanity check (pool = condition (1)), so efficacy passed on removal alone; no arm isolated condition (3). | **NULL_S** (binding): the pool is "(1)+(2)+(4) against some x, D included", so NULL_S differs from A only in (3). C1's null conjunct and Q point to NULL_S; **C3′ added** (judged); **H demoted** to reported; H_NS reported. Constants against NULL_S: see D1. The old NULL is a reported row. | §3, §6.2, §6.4 |
| B1 (build) | Against a pool containing D, the ½ (C1) and ¼ (C3′) factors fail on dev 14/15 and 6/15, because 35-64% of NULL_S's draws are D itself. | **D1**: direction-only against NULL_S; in expectation this is exactly "what (3) accepts is more partial and less `=` than what (3) rejects" (§6.2). ½ and ¼ kept against the D-disjoint old NULL (reported row C1_NULL, falsifier F1). | §6.2, §9 |
| **B2** | Gene-level loss was hidden: `c` swapped for `j`/`k` counts as success in P. | **C5** (judged; A ≤ NULL_S genes losing their only `=`/`c` multi-exon query). Reported: the dropped-`c` × container table, BASE→A gene transitions, `j` drops, and C5 per removed `c`. Opening claim reworded. | §6.3, §7.1, header |
| B2 (build) | Is C5 the intended bar? It fails on dev 14/15, and every loss is a `c` whose container is `j`/`k`/`m`/`n`. | **D2**: C5 as written, binding the **default-candidate tier** (≥ 5 of 6, plus C5 on the curated stratum of V1 and V2). Definitions confirmed: gene_id through the refmap ref_id, multi-exon queries; the name variant is reported. | §6.3, §7 |
| **B3(a)** | 4/6 truths are Gnomon models that collapse 5′-degraded reads the way A does. | Step 1e records, per V3-V6, whether long-read RNA of the substrate's individual was Gnomon evidence (NCBI annotation report). Listed or undeterminable ⇒ **truth-circular**. Circular substrates never carry EFFECTIVE alone, and the efficacy claim cites V1, V2 and S0. | §5.3, §11 |
| **B3(b)** | Safety measured mostly on model-built references. | Strata curated (NM_/NR_) vs model (XM_/XR_). **C2 and C4 on the curated stratum are judged on V1 and V2**; C5.curated binds the tier on V1 and V2; `=` drops, C3′, C2, C4 and C5 are reported per stratum everywhere. | §5.2, §6.1, §6.3 |
| **B3(c)** | "The gate must pass the annotated target first": A applied to a human truth removes about 1% of it. | Per-truth **exposure bound** (A's (1)-(3) on the truth, TOL 0/10/50, per stratum) computed in step 1d and printed beside C2, with a pre-registered reading. | §5.2, §11 |
| **H1** | SIRV could pass vacuously, and the §9.5 falsifier could not fire. | **S0_active** (A drops ≥ 1 in dedup), otherwise S0 is "not measured" and blocks EFFECTIVE. An S1 failure in either mode counts. S2 out of scope ⇒ S0 validates safety only. F5 counts **every** true-ISM `=` removed; the end cause is descriptive. §8 states that S0 has about one isoform of safety power (SIRV303). | §6.4, §8, §9 |
| **H2** | The smoke output was readable, and SIRV was spent on Lever B. | Smoke sealed unread, with sha1s (Amendment 1). **Lever B is not run on SIRV.** After this test, SIRV is a second use for any later prereg. | §3, §10 |
| **H3** | V1 contained r1062-r1085's dev set; V3's exposure missed r1074-r1076. | V1 = A119b minus chr16, chr20, chr21 and chr22; chr21-chr22 reported only; the V3 exposure entry adds r1074-r1076. | §5 |
| **H4** | The user named `m`; the draft re-scoped to `c`. | Put to the user and answered: target = `c` (decision 1). | header |
| **M1** | Dev-plateau thresholds; TOL not derived. | TOL × ρ grid {5,10,20,50} × {0,0.5,1} as **descriptive** held-out rows; plateau falsifier F7 (> 25%). The missing annotation-free derivation of TOL is stated as a limitation, not repaired. | §3, §9, §0 |
| **M2** | Truth files not frozen. | `REF_SHA1.tsv` holds 12 restricted references, and the scorer stops on any mismatch. | Amendment 1 |
| **M3** | `j` drops invisible to P. | `d.dropped_j` is reported beside every class table, together with predicted against observed precision gain. | §7.1 |
| **L1** | Consumers of a default flip. | Flip prerequisites: a Rust port byte-identical to this post-processor, byte-identical when unset, and a family-catalog diff on one dev contig. | §7 |
| **L2** | BASE provenance. | Stated in §3.3. | §3.3 |
| **L3** | The old NULL was called "conservative against A". | Corrected with the measured numbers. | §3.2 |
| **L4** | Register hygiene. | If A is EFFECTIVE: correct §6q6/r866 ("isoseq's 5′-shorter collapse ≡ our ISM collapse"; end-compatibility was the missing part) and amend r1066/r1078 "rule search closed". | §7 |
| build | SQANTI3 is not run. | Q = in-house SQANTI-style FSM, which agrees id for id with SQANTI3 5.5.4 on dev (1,059 = 1,059; 1,585 = 1,585). Its denominator is multi-exon queries, **not** the house "all classified isoforms". | §6.1 |
| build | V1's new key no longer matches the published v3 BASE rows. | Freeze check (a) on the old key, plus an exact **additivity** check: V1 + chr21-chr22 = old key. | §4, §11 |

---

## 1. The rule (binding, verbatim)

```
Input: one assembled GTF (copy_assign --assemble-only), transcripts t with gene_id g(t), contig, strand s(t),
exons sorted by start as 0-based half-open [a, b) (GTF start-1, end), and reads(t) = the integer `reads` attribute
of t's transcript line (the count mode the GTF was assembled in; dedup by default).
chain(t) = ((e_0.b, e_1.a), (e_1.b, e_2.a), ...), m(t) = |chain(t)|.
A transcript y with m(y) >= 1 is DROPPED iff some transcript x != y has
  (1) g(x) = g(y), same contig, s(x) = s(y), and m(x) > m(y);
  (2) chain(y) = chain(x)[k .. k+m(y)-1] for some k      (exact coordinates, contiguous block);
  (3) y.e_0.a >= x.e_k.a - 10  AND  y.e_last.b <= x.e_(k+m(y)).b + 10
      (genomic coordinates: y's two terminal exons lie inside x's corresponding exons, up to TOL = 10 bp);
  (4) 2 * reads(x) >= reads(y)                            (support guard, rho = 0.5).
Every decision is taken against the whole input simultaneously (an x may itself be dropped). The transcript and
exon lines of a dropped y are removed; every other line is written byte for byte and in order.
Name: arm A = `compat_collapse.py --tol 10 --rho 0.5 --mode drop`
      (sha1 e14c5646b133acff731d1335f877d9ad2046fc3a).
```

**How the text is realised.** The implementer has no freedom left.
1. **Scope.** Mono-exon transcripts are never dropped and never serve as x. The lookup key is (gene_id, contig,
   strand), so a gene_id string that repeats on another contig or strand is a different locus.
2. **(3) is symmetric in genomic coordinates** and therefore strand-independent. (4) is implemented as "reject x when
   `reads(x) < 0.5 * reads(y)`", which is exact for integers.
3. **`reads` parse.** The pattern is anchored: `(?:^|; )reads "(\d+)"`, so `matched_reads` cannot shadow it. A
   missing `reads` would silently make ρ 0, so IV6 stops the work if any multi-exon transcript lacks it.
4. **Identities.**
   - (i) A transcript with the largest m in its (gene_id, contig, strand) is never dropped ⇒ **0 loci lost**.
   - (ii) (2) is transitive ⇒ every dropped chain is a contiguous block of a surviving chain ⇒ **0 junctions lost**
     per (gene_id, contig, strand).
   - (iii) *Not guaranteed:* that a surviving x satisfies all of (1)-(4) for y, because TOL and ρ are not transitive.
     On dev it holds for 243/243, 112/121 and 490/500 (chr20 / NC_073244.2 / chr16). This is reported on dev only;
     the frozen scorer does not compute it.
5. **Annotate mode** (`--mode annotate`) keeps every line. It appends `completeness
   "5p_fragment"|"3p_fragment"|"internal_fragment"; fragment_of "<x>";` to y's transcript line, where x is the
   qualifying container with the most reads (ties go to the earlier transcript_id). It is a product, not an arm, and
   the NULLs read D from it. Stripping the two attributes returns the input byte for byte (IV1).

---

## 2. Why this rule (development diagnosis only; design evidence, selected on these contigs)

### 2.1 What the partial classes are (`ct_attr`; dedup; chr20 / NC_073244.2 / chr16)
- **`c` fragments are 5′-truncated** (85 / 92 / 77%). Only 2 / 1 / 8 of them have a ≥ 2-read full chain that we fail
  to emit. For 74-81%, no full chain has ≥ 2 reads (r1077: the read ceiling).
- **Our excess `c` is co-emission.** 206 / 119 / 406 of the 286 / 212 / 552 `c` have an emitted container at the same
  locus. We emit 29.1 / 21.0 contained transcripts per 1,000 queries, against StringTie's 5.5 / 10.2.
- **`m` is already the lowest of the four assemblers.** The retained-intron rule (ratio 10) removes the lopsided
  cases. `k` is 65-82% 5′-extended, and its trimming was refuted (r1078/r1079).

### 2.2 Why end-compatibility separates fragments from real shorter isoforms (`ct_levers`)
- **Co-emitted `c`.** Median terminal-exon overhang past the container's exon is o5 = o3 = 0 bp. 171 / 108 / 333 of
  them lie within 10 bp at both ends.
- **`=` transcripts that also have an emitted container.** Median o5 is 80 / 33 / 55 bp, and only 6 / 10 / 8 lie
  within 10 bp at both ends.
- **The contrast.** A real shorter isoform starts in the container's intron (an alternative first exon). A truncated
  read starts inside the container's exon.
- **The guard (4)** protects annotated short forms whose container is a low-read over-extension. Examples: CD320, 50
  reads against a container of 3; TRMT1, 43 against 5; RPL18, 45 against 3.

### 2.3 Prior work (never re-proposed here)
- r1066 / r863 swept the ISM support ratio only.
- r1077 replaced fragments with raw containers.
- r1072 / r1079 drop or merge the 5′-extended form.
- r861 / r1064 restricted the ISM collapse to 3′-anchored or unprimed chains.
- **None of them tested an end-compatibility condition.** Plain variants fail on dev:
  - any sub-chain: −5.5 to −11.1% chains;
  - 5′-side only: −2.1 to −4.3% chains;
  - ratio 1.0: −1.04 / −1.43% chains on human.
- **Container classes of the dropped `c`** (critique B2): same-reference `=` 31/153, 21/96 and 40/295 (20 / 22 /
  14%). The rest are `c` 65/49/147, `j` 29/14/59, `k`/`m`/`n` 5/0/14, and other-reference `=` 23/12/35. The frozen
  scorer's container ("best surviving sub-chain container, most reads") gives, on chr20: same-ref `=` 35, `c` 55,
  `kmn` 5, `j` 32, other-ref `=` 26.

### 2.4 Dev table (binding seed 20260927; gffcompare 0.12.10 vs the rt6 contig-restricted RefSeq; dedup; frozen scorer)

P = (c+k+m+n) / multi-exon queries (tmap). Precision = `=` multi-exon queries / multi-exon queries
(`c.chain_precision`); on dev it equals matched chains / multi-exon queries.

| contig | arm | multi-exon q | = / c / k / m / n / j | P (rel. vs BASE) | precision % | genes `=` / genes `=`/`c` | drops: = / c / j |
|---|---|---|---|---|---|---|---|
| chr20 | BASE | 5,048 | 1059 / 286 / 169 / 115 / 167 / 2666 | 14.60% | 20.98 | 441 / 480 | |
| | **A** | 4,805 | 1057 / **133** / 166 / 111 / 164 / 2592 | **11.95% (−18.2%)** | 22.00 | 441 / 473 | 2 / 153 / 74 |
| | NULL_S | 4,805 | 1048 / 213 / 167 / 105 / 143 / 2553 | 13.07% (−10.5%) | 21.81 | 440 / 477 | 11 / 73 / 113 |
| | NULL | 4,805 | 1010 / 275 / 158 / 112 / 161 / 2508 | 14.69% (+0.6%) | 21.02 | 433 / 473 | 49 / 11 / 158 |
| NC_073244.2 | BASE | 3,899 | 1595 / 212 / 184 / 66 / 31 / 1740 | 12.64% | 40.91 | 823 / 885 | |
| | **A** | 3,778 | 1592 / **116** / 184 / 66 / 31 / 1719 | **10.51% (−16.9%)** | 42.14 | 821 / 878 | 3 / 96 / 21 |
| | NULL_S | 3,778 | 1586 / 144 / 183 / 63 / 30 / 1702 | 11.12% (−12.1%) | 41.98 | 822 / 880 | 9 / 68 / 38 |
| | NULL | 3,778 | 1565 / 207 / 175 / 61 / 29 / 1672 | 12.49% (−1.2%) | 41.42 | 821 / 884 | 30 / 5 / 68 |
| chr16 | BASE | 8,782 | 1605 / 552 / 269 / 261 / 290 / 4716 | 15.62% | 18.28 | 707 / 780 | |
| | **A** | 8,282 | 1602 / **257** / 266 / 253 / 282 / 4545 | **12.77% (−18.2%)** | 19.34 | 705 / 771 | 3 / 295 / 171 |
| | NULL_S | 8,282 | 1594 / 421 / 256 / 232 / 258 / 4456 | 14.09% (−9.8%) | 19.25 | 704 / 773 | 11 / 131 / 260 |
| | NULL | 8,282 | 1519 / 533 / 247 / 245 / 268 / 4388 | 15.61% (−0.1%) | 18.34 | 696 / 767 | 86 / 19 / 328 |

- Matched chains: 1,059 → 1,057, 1,595 → 1,592 and 1,605 → 1,602 under A (−0.19% each).
- IV3 holds on all three contigs: 0 loci, 0 junctions and 0 mono-exon transcripts lost.
- The in-house FSM share (Q) moves 20.98 → 22.00, 40.65 → 41.87 and 18.26 → 19.33 under A. NULL_S gives 21.81,
  41.72 and 19.23.
- The precision gain is almost arithmetic. With d = drops / multi-exon queries and base precision p, the gain is
  about p·d/(1−d): predicted +1.06 / +1.31 / +1.10 pt against observed +1.02 / +1.23 / +1.07.

### 2.5 NULL_S: why constants fail against it, and what (3) does per transcript (dev; `summary.md`)

- **Pool and overlap.**
  - The NULL_S pool has 692 / 238 / 1,392 transcripts, so f = |D| / |pool| = 35.1 / 50.8 / 35.9%.
  - Because of the strata, the realised share of NULL_S draws that fall in D is **35-64% over the 5 seeds** (binding
    seed: 98/243, 77/121, 201/500).
- **Decomposition, binding seed.** "Partial" means c+k+m+n by BASE class.

  | contig | D = A's drops: partial / `=` | NULL_S draws outside D: partial / `=` | NULL_S draws inside D: partial / `=` |
  |---|---|---|---|
  | chr20 | 163/243 (67%) / 2 (0.8%) | 41/145 (28%) / 10 (6.9%) | 68/98 (69%) / 1 |
  | NC_073244.2 | 96/121 (79%) / 3 (2.5%) | 10/44 (23%) / 8 (18.2%) | 63/77 (82%) / 1 |
  | chr16 | 314/500 (63%) / 3 (0.6%) | 83/299 (28%) / 10 (3.3%) | 122/201 (61%) / 1 |

- **Condition (3) is selective per transcript.**
  - What it accepts is 2.3-3.4× more partial than what it rejects within the same pool.
  - Its `=` rate relative to the rejected transcripts' is .12 / .14 / .18, below ¼.
  - Its partial rate passes the ½ factor too: rejected 28 / 23 / 28% against half of accepted, 33 / 40 / 31%.
- **The constant-factor clauses against NULL_S nevertheless fail** (5 seeds × 3 contigs):
  - "NULL_S reduction < ½ of A's" passes 1/15;
  - "4·e_A ≤ e_NS" passes 9/15 (chr20 5/5, NC_073244.2 0/5, chr16 4/5);
  - H_NS passes 0/15: NULL_S obtains 78-91% of A's precision gain.

  §6.2 derives why: a pool that contains D hands NULL_S a share of A's own drops.
- **Direction-only against NULL_S** (D1) passes 15/15 on each of C1s and C3′:
  - NULL_S relative partial-share change is −8.5…−10.5 / −10.0…−12.1 / −9.6…−11.2%, against A's −18.2 / −16.9 / −18.2%;
  - e_A = 2 / 3 / 3 against e_NS = 9-13 / 8-10 / 11-17.
- **C5 (genes losing their only multi-exon `=`/`c` query)** fails 14/15. A loses 7 / 7 / 9 genes; NULL_S loses 2-4 /
  2-5 / 4-9. The single pass is a tie, 9 vs 9, on chr16 seed 2.
  - On the curated stratum at the binding seed, A vs NULL_S is 6 vs 3 (chr20), 0 vs 0 (NC_073244.2) and 9 vs 8
    (chr16).
  - Every A loss is a gene whose only `=`/`c` query was a `c` whose surviving container is `j`, `k`, `m` or `n`.
    On chr20 the BASE→A transitions are c>j 5, c>k 2 and m>n 1.
  - **Per removed `c`**, A's gene losses are 7/153, 7/96 and 9/295 (4.6 / 7.3 / 3.1%), against NULL_S's 3/73, 5/68
    and 7/131 (4.1 / 7.4 / 5.3%). The excess tracks A's `c` removal. This is descriptive and not a clause.
- **Old NULL** (condition (1) pool, D excluded).
  - The C1 conjunct with ½ passes 15/15.
  - The `=` share among A's drops against the old NULL's is .04 / .10 / .035 (F1 bar ¼).
  - Its `c` share among drops is 4.5 / 4.1 / 3.8%, at or below BASE's `c` share of 5.7 / 5.4 / 6.3%. The old NULL
    pool is **not** partial-enriched, contrary to the draft's "conservative against A" (critique L3).

### 2.6 TOL × ρ grid (descriptive; `grid_table.md`)
- **At ρ 0.5, `c` removal varies by 6-8% relative across TOL 5-50:** 149-158, 90-97 and 287-310 transcripts. This is
  below the 25% plateau falsifier (F7).
- **TOL 0 halves the effect.** P changes by −8.6 / −10.6 / −8.7%, which would fail C1e on chr20 and chr16.
- **C5 losses come with the end tolerance.** At ρ 0.5 they go from 0 / 3 / 5 at TOL 0 to 7 / 7 / 9 at TOL 5.
- **ρ 1** gives P −11.7 / −9.9 / −10.5% and C5 losses of 5 / 2 / 6.
- **ρ 0** adds `=` drops (6 / 10 / 8) and C5 losses (8 / 12 / 10).
- No variant is substituted for §1 (§11 stop rules).

### 2.7 Truth-side exposure: A's (1)-(3) applied to the truth itself (`exposure_dev.txt`)

| truth | multi-exon refs | annotated ISMs | end-compatible at TOL 0 / 10 / 50 | removed share at TOL 10: all / curated / model |
|---|---|---|---|---|
| chr20 RefSeq | 4,295 | 274 | 34 / 42 / 58 | 0.98% / 1.72% (38/2,204) / 0.19% (4/2,079) |
| chr16 RefSeq | 6,690 | 489 | 54 / 73 / 90 | 1.09% / 1.54% (52/3,380) / 0.65% (21/3,235) |
| NC_073244.2 Gnomon | 5,568 | 348 | 3 / 14 / 36 | 0.25% / 1 of 6 curated / 0.24% (13/5,468) |
| SIRV E0 | 61 | 9 | 1 / 1 / 1 | 1.6% (SIRV303 only) |

- The curated human stratum exposes 6-7× the share of end-compatible real isoforms that the Gnomon gorilla truth does
  (1.54-1.72% against 0.24%). Safety is therefore judged on the curated stratum of V1 and V2 (§6.1).

### 2.8 SIRV truth structure (truth annotation only)
- Of the 9 true ISMs, only **SIRV303** satisfies (1)-(3) against a truth container (SIRV301, as a `3p_fragment`). This
  holds at TOL 0, 10 and 50.
- A 5′-truncated emitted copy of a true ISM whose TSS lies in the container's intron can still fall inside the
  container's exon and be removed. F5 counts every such removal.

### 2.9 Projection to held-out (arithmetic on dev; NOT a measurement, NOT in any verdict)

| substrate | BASE matched chains | BASE precision | Δprecision if d = 3% / 5% | the (reported) house bar H needs d ≳ |
|---|---|---|---|---|
| human_A119b −chr16/chr20 (published key; V1 now also drops chr21/chr22) | 35,956 | .1765 | +0.55 / +0.93 | 5.4% |
| human_testis | 11,399 | .4596 | +1.42 / +2.42 | 2.1% |
| gorilla_OR6737 −NC_073244.2 | 24,277 | .3496 | +1.08 / +1.84 | 2.8% |
| gorilla_KB3781 | 25,911 | .3594 | +1.11 / +1.89 | 2.7% |
| chimp_PTR | 21,031 | .3915 | +1.21 / +2.06 | 2.5% |
| orangutan_PPY | 20,243 | .1931 | +0.60 / +1.02 | 4.9% |

H is mostly a statement about BASE precision (critique B1). That is why it is reported and does not decide a tier.

---

## 3. Arms

| arm | what | status |
|---|---|---|
| **BASE** (genome) | the current default, readthrough unset: `rt_arms/<s>/<s>.BASE.gtf` (→ `runs/<s>/<s>.gtf`), dedup; sha1s in Amendment 1 | **reused** (§3.3) |
| **A** (genome) | §1 on BASE (`run_arms.sh`) | **new, judged** |
| **NULL_S** (genome) | §3.1 on BASE, binding seed 20260927 | **new, binding comparator** of C1s, C3′, Q and C5 |
| NULL (genome) | the draft's §3.1 (§3.2), binding seed | new, **reported** (row C1_NULL; falsifier F1) |
| A∘R3 (genome) | §1 on `rt_arms/<s>/<s>.R3.gtf`, scored against R3 | new, **descriptive only** (the user may flip R3) |
| grid (genome) | `compat_collapse.py --tol T --rho R --mode drop` on BASE, T ∈ {5, 10, 20, 50}, R ∈ {0, 0.5, 1}, minus (10, 0.5) | new, **descriptive only** (F7) |
| **BASE_dedup / BASE_kcd** (SIRV) | frozen `copy_assign` f480847a on `ct_sirv/bam/sirv_testis.bam`. dedup = the driver `tools/rustle_pipeline.sh assemble` (sha1 d1816c35); kcd = the driver's `stage_assemble` command expanded, plus `--keep-coordinate-duplicates` | **new** (fresh substrate) |
| **A_dedup / A_kcd** (SIRV) | §1 on each | **new, judged** (§6.4) |
| NULL_S / NULL (SIRV) | §3.1 / §3.2 on each | new, **descriptive only** (too few draws for a null test) |

- **Not run anywhere in this test:**
  - Lever B (`--polish-ism-ratio 0.5` + A), including on SIRV (critique H2);
  - extra NULL seeds on held-out;
  - any TOL/ρ variant as a judged arm;
  - the L3b orphan filter.
- The held-out stage answers "A vs BASE, with NULL_S as the (3)-isolating comparator" only.

### 3.1 NULL_S (binding; `complete_null.py --kind NULL_S`, sha1 9740f092eb2dec93f5e31683a7f3bed7bc6e8438)
- **Dropped set.** D = the transcripts whose transcript line carries `completeness` in A's annotate-mode output.
- **Pool.** Multi-exon transcripts y of BASE, **D included**, for which some x ≠ y of the same (gene_id, contig,
  strand) satisfies A's (1), (2) and (4): m(x) > m(y); chain(y) is a contiguous block of chain(x) with exact
  coordinates; 2·reads(x) ≥ reads(y).
  - NULL_S therefore differs from A only in condition (3).
  - The script asserts D ⊂ pool.
- **Strata.**
  - Locus class = the number of multi-exon transcripts of the (contig, strand, gene_id) in BASE, binned 2 | 3 | 4-5 |
    6-9 | ≥ 10.
  - Read bin = floor(log2(max(reads, 1))).
- **Draw.**
  - Per contig, the strata are taken in sorted order. From each, the script samples |D ∩ stratum| without replacement
    from the pool's candidates, sorted by transcript_id.
  - RNG: `random.Random(f"20260927:{sample}:{contig}:S")`.
  - A shortfall is filled from the same locus class at the nearest read bin (the lower bin on a tie), then from any
    pool transcript of the contig. Both counts are logged. "Unfilled" is reported, and the arm stands as drawn.
  - The log line gives |pool| and the number of draws inside D. These are reported as f (§7.1).
- **Action.** A's: the drawn ids' transcript and exon lines are removed; everything else is kept byte for byte.
- **Scope.** Drawn on every contig of the GTF; the scorer restricts to the substrate's contigs. RNGs are per contig,
  so dev contigs do not affect held-out draws.
- `--sample` is the sample id (`human_A119b`, …); on SIRV it is `sirv_E0_human_testis`.

### 3.2 NULL (reported; `complete_null.py --kind NULL`)
- The pool is multi-exon transcripts **not in D** that have a multi-exon sibling of the same (contig, strand, gene_id)
  with strictly more exons: condition (1) alone.
- Strata, draw, action and scope are as in §3.1. The RNG is `random.Random(f"20260927:{sample}:{contig}")`.
- It is byte-identical to the draft prototype bcfc5a50.
- It is **not** partial-enriched (§2.5), so against it the ½ and ¼ constants test only "A is not random".

### 3.3 Provenance of the reused BASE (critique L2)
- The held-out BASE GTFs were built on 09-25 by `copy_assign` 753a3b4d… (the driver key), not by f480847a. The
  v2/v3 IV1 chain (`PREREG_readthrough_v3_2026-09-26.md` Amendment 1) makes the two byte-equivalent in their outputs.
- The chr20 rows of the genome-wide A119b BASE equal the dev BASE: 5,631 = 5,631 transcripts, identical coordinates
  and `reads`. Only TPM differs, and A never reads it.
- **The reuse is valid only if** each BASE sha1 at scoring equals Amendment 1's record and freeze check (a) passes
  (§4, §11 step 1b).

---

## 4. Gates

Any failure is fixed in the code, never in the rule, NULL_S or a clause. A fix after this file is an amendment that
re-runs every A, NULL_S and NULL arm.

**Passed before this file (development contigs, truths and tool GTFs only; Amendment 1 records them):**
- **IV1, byte identity.**
  - `compat_collapse.py --mode drop` on a GTF with no qualifying pair returns it unchanged.
  - Drop mode = the input minus exactly the listed ids' lines, in order.
  - Annotate mode with the two attributes stripped = the input.
  - Annotate-mode gffcompare `.tmap`/`.stats` = BASE's (3/3).
  - Nothing is edited in Rust, so "unset" means the post-processor is not applied.
- **IV2, port = selection.** The frozen `compat_collapse.py` reproduces `ct_levers` `L1e_compat10_cont_ge0.5`'s dropped
  sets exactly: 243 / 121 / 500 ids. Drop and annotate outputs are byte-identical to `ct_prereg` e14c5646.
- **IV3, identities.** 0 loci, 0 junctions and 0 mono-exon transcripts are lost under A and NULL_S. The NULLs are
  deterministic, fill their targets exactly (unfilled 0, drawn from another stratum 0) and never draw a locus's
  longest transcript.
- **IV4 + scorer fixtures.** `complete_eval.py selftest` passes 29/29:
  - TOL 10 in / 11 out at both ends and on both strands;
  - the ρ edge;
  - another gene, the other strand, mono-exon;
  - y ⊂ x ⊂ z;
  - the `matched_reads` parse;
  - NULL_S determinism and pool;
  - class, gene, stratum, container and clause logic.
- **IV5, scorer = spec.** The frozen scorer reproduces the draft's dev table exactly (`ct_prereg/dev_table.tsv`
  707226c6: 11 integer columns × 7 arms × 3 contigs). It also reproduces `readthrough_eval.py`'s (f6b99dcc)
  `c.chain_precision` and matched chains on the same BASE.
- **IV6, inputs.** Every multi-exon transcript of the dev BASEs carries an integer `reads`. gffcompare is 0.12.10
  (binary sha1 32cbdd64).
- **SIRV technical gate (tools only).**
  - StringTie / FLAIR / IsoSeq reproduce `ct_sirv`'s table: n 44 / 64 / 109; `=` 39 / 49 / 62; sensitivity 39 / 49 /
    47 of 61.
  - The truth scored as a query gives 61/61 `=`.
  - A 9-drop line-subset of StringTie exercises every S-row.
- **Guards.** The scorer refuses:
  - a non-dev substrate without `--heldout --prereg`;
  - a prereg without an "Amendment 1" heading and the scorer's own sha1;
  - a SIRV arm that is not a frozen tool GTF or a line-subset of one;
  - a reference whose sha1 differs from `REF_SHA1.tsv`;
  - a gffcompare other than 0.12.10.

**Held-out gates, run in step 1 after this file (BASE and truth only; recorded in Amendment 2 before any held-out or
SIRV score of A, NULL_S or NULL):**
- **IV5(a), freeze check.** BASE alone is re-scored on the six keys of the published v3 tables:
  - `human_A119b --drop-contigs chr16,chr20` (reference fe75a714);
  - `human_testis`;
  - `gorilla_OR6737 --drop-contigs NC_073244.2`;
  - `gorilla_KB3781`, `chimp_PTR`, `orangutan_PPY`.

  It must reproduce `c.matching_intron_chains` and the `c.chain_precision` k/n of `rt_arms/tables_v3/` exactly. Only
  these rows are read.
- **IV5(a′), V1 additivity.** BASE on V1 (`--drop-contigs chr16,chr20,chr21,chr22`, reference c89dfe1a) plus BASE on
  `--contigs chr21,chr22` (reference 5d37e1d7) must equal the fe75a714 key exactly in matched chains, `=` multi-exon
  queries and multi-exon queries.
- **IV6** on each BASE: the scorer row `iv6.multi_without_reads` = 0.
- **IV7, exposure.** `complete_eval.py exposure --heldout` on V1-V6 (§5.2).
- **IV8, circularity record** (§5.3).
- **Any mismatch stops the work.**

**Per arm, before its clause rows are used:** `iv3.loci_lost` = `iv3.junctions_lost` = `iv3.dropped_mono` = 0 for A,
NULL_S and NULL, and `d.added` = 0. Otherwise stop.

---

## 5. Substrates, strata and exposure

Each substrate is judged on its own, and species are never pooled. The dev contigs (A119b chr16 and chr20; OR6737
NC_073244.2) are excluded everywhere.

| id | substrate (scorer key) | reads / count mode | prior exposure |
|---|---|---|---|
| **S0** | SIRV E0 of human_testis, SIRV1-7 (`--sirv`; truth `sirv_C.gtf` 81051062) | 11,103 primaries, 0 secondaries (seeding is a no-op); **dedup and kcd** | **fresh**: truth built 09-27; tool baselines seen; our output never scored; smoke sealed unread |
| V1 | human_A119b `annotated_minus_chr16_chr20_chr21_chr22` (c89dfe1a) | dedup | R/RQ1, R3 verdicts; BASE classes in `eval_all` and figs 1-2 (chr20-22) |
| V2 | human_testis `annotated` (8ed998a1); disjoint from S0, because the SIRV reads are this BAM's unmapped reads | dedup | R/RQ1, R3 verdicts; BASE classes in `eval_all` |
| V3 | gorilla_OR6737 `annotated_minus_NC_073244.2` (9271749f) | dedup | R/RQ1, R3 verdicts; **held-out substrate of r1074-r1076** (the current polish defaults); r1111; BASE classes and SQANTI3 in figs 1-2 |
| V4 | gorilla_KB3781 `annotated` (a7e0e06e) | dedup | R/RQ1, R3 verdicts; BASE classes in `eval_all` |
| V5 | chimp_PTR `annotated` (f4b58b50) | dedup | R/RQ1, R3 verdicts; BASE classes in `eval_all` |
| V6 | orangutan_PPY `annotated` (35a55052) | dedup | R/RQ1, R3 verdicts; BASE classes in `eval_all` |
| (rep.) | human_A119b `chr21-chr22` (5d37e1d7) | dedup | dev set of r1062-r1085; **scored, never judged** |

### 5.1 Exposure
- No number of A or NULL_S was ever computed on V1-V6 or S0, and A was designed without them.
- What was exposed is BASE. The fourth reuse of the six libraries weakens each claim more than the last. S0 is the
  one fresh test, and it is small (§6.4).
- **Power floors.**
  - The R prereg's floors hold on V1-V6: v3 §5 records ≥ 1,000 genes and ≥ 500 BASE matched chains on the same BASE
    files, and V1 minus chr21/22 keeps both.
  - C1e needs ≥ 10 BASE partials (1/X). This is met genome-wide.
  - S2 has its own floor (§6.4).

### 5.2 Strata and the exposure bound (critique B3(b), B3(c))
- **Strata** use the accession of the matched reference transcript (`rna-` stripped): **curated** = NM_/NR_,
  **model** = XM_/XR_, other, none.
  - The scorer reports C2, C4, C5, C3′ and the `=` drops per stratum.
  - **Judged on the curated stratum of V1 and V2:** C2.cur and C4.cur (§6.1), and C5.cur, which binds the tier
    (§6.3).
  - On V3-V6 the curated stratum is too small to judge (dev gorilla: 11 of 3,899 queries). It is reported.
- **Exposure bound X(s)** = A's (1)-(3) at TOL 10 applied to the restricted truth, as a share of its multi-exon
  references, per stratum. It is computed by `complete_eval.py exposure` in step 1d, which uses no reads and no
  assembler output.
- **Pre-registered reading:**
  - **(i)** X(s) is printed beside C2 and C2.cur in the Outcome.
  - **(ii)** If X_all(s) or X_cur(s) ≥ 1% (C2's tolerance), the substrate is labelled "C2 tolerance ≤ truth exposure".
    A perfect assembler would then lose about that share of `=` to A. C2 is still judged there: a real isoform
    lost is a real loss.
  - **(iii)** A's realised loss `g.refs_eq_lost` is reported as a fraction of the exposure count.

### 5.3 Truth circularity (critique B3(a))
- **Question.** Were long reads of the substrate's own individual Gnomon evidence for its truth?
- **Procedure (step 1e, recorded in Amendment 2 before any held-out score of A).** For each of V3-V6, read the NCBI
  Eukaryotic Genome Annotation report of its assembly:
  - gorilla GCF_029281585.2, for both gorilla samples;
  - chimp GCF_028858775.2;
  - orangutan GCF_028885625.2.

  Record whether any long-read RNA run (PacBio Iso-Seq or ONT) is listed as transcript evidence that comes from the
  substrate's individual:
  - KB3781 is the mGorGor1 individual (SRR27438212);
  - OR6737 is a different gorilla;
  - PTR and PPY are of unknown individual (`figures/samples.tsv`).
- **Rule.** Listed ⇒ **truth-circular**. Not determinable ⇒ **truth-circular (unknown)**, treated the same.
- **Effect.** A circular substrate is judged like any other: its failures count and its passes are needed for 6/6.
  EFFECTIVE, however, always also requires V1 and V2 (non-circular human RefSeq), so it never rests on a circular
  substrate alone. The Outcome states the efficacy claim on V1, V2 and S0, and on any Gnomon substrate shown
  non-circular.

---

## 6. Metrics and clauses

- All gffcompare numbers come from 0.12.10, with `-r` = the sample's annotation restricted to the substrate contigs
  (`assembly.restrict_gtf` rule; sha1 checked against `REF_SHA1.tsv`).
- "Multi-exon queries" = `.tmap` rows with `num_exons ≥ 2`, each emitted transcript counted once.
- Every clause is an integer test on the frozen scorer's rows. Subscripts: B = BASE, A = arm A, S = NULL_S, N = NULL.

### 6.1 Genome substrates V1-V6: P(s) = every clause below passes

| clause | measure (scorer row) | passes iff | tolerance source |
|---|---|---|---|
| **C1** efficacy + (3)-selectivity | p = c+k+m+n multi-exon queries, q = multi-exon queries (`c.partial` k, n) | **C1e**: 10·p_A·q_B ≤ 9·p_B·q_A (≥ 10% relative reduction) **and C1s**: p_A·q_S < p_S·q_A (A's partial share strictly below NULL_S's) | X = 10% = A1 of the R prereg (`PREREG_readthrough_ends_representatives_2026-09-25.md` §5); C1s is direction-only (D1) |
| **C2** chains | `c.matching_intron_chains` (M) | 100·M_A ≥ 99·M_B | Y = 1% = G2 of the R prereg §5 and the polish-lever bar (`PREREG_assembly_precision_levers_2026-09-23.md`) |
| **C3** precision | `c.chain_precision` k (`=` multi-exon) / n | k_A·n_B ≥ k_B·n_A | G1 of the R prereg (non-inferior, no constant) |
| **C3′** selectivity | e = `d.dropped_eq` (`=` multi-exon among the dropped), t = `d.dropped_multi` | e_A·t_S < e_S·t_A, **or** e_A = e_S = 0 | direction-only (D1); critique B1 item 3, with the ¼ moved to F1 |
| **C4** genes with `=` | `g.genes_eq` (G) | 100·G_A ≥ 99·G_B | Y = 1% = G4 of the R prereg |
| **Q** FSM share | `q.fsm` k (F) / n (multi-exon queries; in-house SQANTI-style FSM, not the house all-isoform denominator) | F_A·q_B > F_B·q_A **and** F_A·q_S > F_S·q_A (scorer row `Q`) | direction only |
| **C2.cur** (V1, V2 only) | `s.curated.refs_eq` (R_cur) | 100·R_cur,A ≥ 99·R_cur,B (scorer row `C2.curated`) | Y = 1%; critique B3(b) |
| **C4.cur** (V1, V2 only) | `s.curated.genes_eq` (G_cur) | 100·G_cur,A ≥ 99·G_cur,B (scorer row `C4.curated`) | Y = 1%; critique B3(b) |

- The scorer writes C2, C3, C4, Q, C2.curated and C4.curated as clause rows. C1e, C1s and C3′ are computed from
  the metric rows named in the table, and the note of the scorer's `C1` row repeats C1e ("efficacy True/False").
- **A clause that cannot be computed is "not measured":** it is not a pass, it blocks EFFECTIVE, and it is not a
  failure. Example: a stratum clause whose BASE count is 0. C3′ uses the rate form, so unequal drop counts (an
  unfilled NULL_S) do not make it unmeasurable.

### 6.2 D1: why every comparison against NULL_S is direction-only
- **Notation.** Let f be the share of NULL_S's draws that fall in D. Let π and ρ be the partial and `=` rates of the
  transcripts (3) **accepts** (D) and **rejects** (R = pool ∖ D). Let P_B be BASE's partial share. The algebra below
  is for a uniform draw; strata change the constants, not the argument.
- **Direction.** In expectation, "A's partial share < NULL_S's" ⟺ π_D > π_R, and "e_A < e_S" ⟺ ρ_D < ρ_R. These are
  exactly the claims that condition (3) selects partials and avoids complete transcripts. They have no constant and
  do not depend on f.
- **The ½ conjunct against NULL_S** ⟺ (f − ½)·π_D + (1 − f)·π_R < ½·P_B. At f = ½ this requires π_R < P_B: the
  sub-chains that (3) rejects would have to be less partial than the whole assembly. On dev π_R is 23-28% and P_B
  12.6-15.6%, and the realised f is 35-64%.
- **The ¼ conjunct against NULL_S** ⟺ ρ_D/ρ_R ≤ (1 − f)/(4 − f). The bar is ¼ only at f = 0 and 0.14 at f = ½.
- **Result.** Both constants measure the pool's size (f) as much as (3)'s selectivity. Against a D-disjoint
  comparator they mean what they were designed to mean. That comparator is the old NULL, reported as row
  `C1_NULL` and falsifier F1. Descriptively it is also D against NULL_S's own draws outside D (§2.5, §7.1).
- **Reported, not judged:** the scorer's critique-literal rows `C1` (½ against NULL_S) and `C3'` (¼ against NULL_S),
  labelled C1½NS and C3′¼NS in the Outcome, and H_NS.

### 6.3 D2: the default-candidate tier T (C5; critique B2)

| clause | measure (scorer row) | passes iff | where |
|---|---|---|---|
| **C5** | L = `g.genes_eqc_lost`: reference genes (gene_id of the refmap ref_id) with ≥ 1 multi-exon `=`/`c` query in BASE and none in the arm | L_A ≤ L_S (scorer row `C5`) | V1-V6 |
| **C5.cur** | the same on curated-stratum references (`s.curated.genes_eqc_lost`) | L_cur,A ≤ L_cur,S (scorer row `C5.curated`) | V1, V2 |

**T holds iff C5 passes on ≥ 5 of the 6 genome substrates and C5.cur passes on both V1 and V2.** The "≥ 5 of 6"
follows the v3 §7 precedent for a head-to-head bar.

- **Why a tier and not P(s).** Every dev C5 loss is a `c` whose surviving container is `j`, `k`, `m` or `n` and holds
  the removed chain as a contiguous block (identity (ii)). The gene keeps the junctions but loses gffcompare's
  intron-correct label. That is a real cost for a default, since every user would see it, but not a loss of
  assembled evidence. On an opt-in, and with the annotate product, the cost is the user's choice.
- **C5 is kept exactly as the critique wrote it**, with no constant and no normalisation. The per-removed-`c` rate
  (§2.5) is reported only.

### 6.4 SIRV (S0): exact counts against the complete truth

| clause | measure (scorer row) | passes iff | judged in |
|---|---|---|---|
| **S0_active** | A's drops (`d.dropped`) | ≥ 1 | dedup; if 0, S0 is **"not measured"** |
| **S1** safety | truth multi-exon isoforms that are the tmap ref of ≥ 1 multi-exon `=` query (`g.refs_eq_tmap`) | A ≥ B (0 tolerance: 1% of 61 < 1 isoform) | **dedup and kcd**; a failure in either mode is an S0 failure even when S0 is inactive |
| **S2** efficacy | partial share, as C1e | 10·p_A·q_B ≤ 9·p_B·q_A | dedup, only if p_B ≥ 10 (= 1/X); otherwise "outside scope", and S0 then validates **safety only** |
| **S3** precision | `=` multi-exon / multi-exon | k_A·n_B ≥ k_B·n_A | dedup |
| **Q0** FSM | in-house FSM share | F_A·n_B ≥ F_B·n_A (non-inferior) | dedup |

- **S0 passes iff** S0_active holds and S1 (dedup), S1 (kcd), S3, Q0 and, when in scope, S2 pass.
- **S0 fails iff** any of these that is measured fails, including S1 in either mode.
- **Why S0 cannot pass vacuously.** With S0_active false, S0 is "not measured", which is not a pass and blocks
  EFFECTIVE. S3 and Q0 are near-tautological whenever A acts on a complete truth, so the judged content of S0 is S1
  and, when in scope, S2.
- **Pre-registered SIRV facts** (`ct_sirv`):
  - 9 true ISMs; only SIRV303 is structurally at risk (§2.8);
  - SIRV701/705 share a chain, so the `=` ceiling is 60/61;
  - SIRV107 is expected unreachable (strand without `-uf`, and non-canonical introns under strict junctions);
  - 93.3% of primaries are coordinate duplicates, so dedup feeds about 15× fewer reads than kcd;
  - in kcd, (4) sees true counts.
- **Reported:** sensitivity /61 and /57 (`sirv.sens61`, `sirv.sens_ge2reads`) and precision beside StringTie / FLAIR
  / IsoSeq; the class table per mode; every `sirv.eq_lost.<isoform>` row with the removed query's end offsets
  (feeds F5).

### 6.5 The numbers the clauses compare against

All bars are relative to BASE's counts computed at scoring. The published rows below serve freeze check (a) and
orientation only.

| substrate | published BASE matched chains → C2 needs ≥ | published BASE precision |
|---|---|---|
| human_A119b −chr16/chr20 (freeze-check key; V1's own bar is computed at scoring) | 35,956 → (35,597 on the old key) | .1765 |
| human_testis | 11,399 → 11,286 | .4596 |
| gorilla_OR6737 −NC_073244.2 | 24,277 → 24,035 | .3496 |
| gorilla_KB3781 | 25,911 → 25,652 | .3594 |
| chimp_PTR | 21,031 → 20,821 | .3915 |
| orangutan_PPY | 20,243 → 20,041 | .1931 |

---

## 7. Verdicts (arm A only)

**Failing substrate.** A genome substrate s fails when P(s) is false through a clause that was measured. S0 fails as
defined in §6.4. S0 counts as one substrate, so there are seven.

| verdict | condition | meaning |
|---|---|---|
| **EFFECTIVE, default candidate** | P(s) on 6/6 (with C2.cur and C4.cur on V1 and V2); S0 measured and passing; **T** holds | A becomes the default candidate: a polish pass after the retained-intron filter. **The flip is the user's call** and needs the L1 prerequisites below. |
| **EFFECTIVE, opt-in** | as above, but T fails | Ship A as an opt-in flag (drop mode) plus the annotate product. The partial share falls ≥ 10%, (3) is selective, and no `=` chain or gene is lost beyond 1%, but A costs more intron-correct gene labels than NULL_S. The Outcome states that cost as numbers. |
| **KEEP OPT-IN** | not EFFECTIVE, and failing substrates ≤ 1 of 7; or only "not measured" blocks EFFECTIVE (e.g. S0 inactive) | opt-in, with a weaker claim |
| **REFUTE** | failing substrates ≥ 2 of 7 | not shipped; the annotate product (§1 item 5) stays a separate user decision |

- A missing substrate makes A "undecided", unless ≥ 2 present substrates already fail.
- Species are never pooled, and dev never enters the verdict.
- A∘R3, the grid, the chr21-chr22 rows, the old NULL, the SIRV NULLs, H, H_NS, C1½NS and C3′¼NS are never judged.
- The verdict is computed by the integer tests of §6 from the score TSVs. A helper written for this is allowed only
  if it first reproduces the dev clause table of Amendment 1 item 6 exactly; its sha1 goes into the Outcome.
- **Flip prerequisites (critique L1; not part of the verdict).** A Rust port must:
  - reproduce this post-processor's output byte for byte on the three dev contigs and the six BASE GTFs;
  - be byte-identical to f480847a when unset;
  - come with a family-catalog diff (`mcl_families` stage) on one dev contig, reported.
- **Register hygiene (critique L4), if EFFECTIVE:**
  - correct §6q6/r866's "isoseq's 5′-shorter collapse ≡ our ISM collapse" (end-compatibility was the missing part);
  - amend r1066's and r1078's "rule search closed".

### 7.1 Reported beside the verdict (not clauses)
- **Per substrate, for A, NULL_S and NULL:**
  - the full class table and the `c` share;
  - the ISM subtypes (`q.ism_3prime/5prime/internal`);
  - d = drops / multi-exon queries;
  - `=`, `c` and **`j` among the drops** (critique M3);
  - `=` drops per stratum;
  - the dropped-`c` × container table (`d.c_container.*`);
  - the BASE→arm gene transitions (`d.gene_transition.*`);
  - `g.genes_eqc_name_lost` (the critique's name variant);
  - C5 per removed `c`;
  - C2, C4 and C5 per stratum;
  - the reference chains lost and gained.
- **Selectivity context:**
  - f (draws of NULL_S inside D, and |pool|, from `null.log`);
  - C1½NS, C3′¼NS, C1_NULL, H and H_NS;
  - the precision gain, predicted p·d/(1−d) against observed.
  - If a helper that reproduces §2.5's dev decomposition exactly (`ct2_build/summarize.py` d5569e18) is run, also D
    against NULL_S's own draws outside D (π and ρ).
- **Exposure:** X(s) all / curated / model beside C2 and C2.cur (§5.2), and the circularity label (§5.3).
- **Descriptive arms:**
  - A∘R3 against R3 (C1e, C2-C4, Q);
  - chr21-chr22 (all rows);
  - the TOL × ρ grid (c removal, P, `=` drops, C5 losses).
- **SIRV:**
  - every S-row in both modes;
  - class tables;
  - sensitivity and precision beside the tools;
  - `sirv.eq_lost.*` with end offsets.

---

## 8. Predictions (before any A or NULL_S number on a held-out or SIRV substrate)

1. **C1e** passes on 6/6 (p ≈ 0.7). The median relative reduction is 12-18%. The weakest substrates are the two
   gorillas (dev co-emission 56% against 72-74% human) and human_testis (another library type).
2. **C1s** passes on 6/6 (p ≈ 0.8). The realised f is 35-65%, and NULL_S reduces P by about half of A's reduction.
3. **C2, C3 and C4** pass on 6/6. The chain loss is 0.1-0.4% and the gene loss ≤ 0.3%. **C2.cur and C4.cur** pass on
   V1 and V2 (p ≈ 0.85).
4. **C3′** passes on 6/6 (p ≈ 0.8). **Q** passes on 6/6 (p ≈ 0.8). Q nearly tracks C3, because in-house FSM ≈ `=`.
5. **C5 fails on ≥ 2 of 6 (p ≈ 0.85), so T fails.**
   - A loses about 0.8-1.5% of `=`/`c` genes. NULL_S loses about half that.
   - C5.cur fails on at least one of V1/V2 (p ≈ 0.8).
6. **Reported rows.**
   - H passes where BASE precision ≥ .35 (V2-V5) and fails on V1, with V6 borderline.
   - H_NS fails on ≥ 5 of 6.
   - C1½NS fails on ≥ 4 of 6, and C3′¼NS on ≥ 2 of 6.
7. **SIRV.**
   - S0_active in dedup: p ≈ 0.7. S2 in scope: p ≈ 0.5.
   - S1 passes in dedup (p ≈ 0.75) and kcd (p ≈ 0.7). A loss, if any, is SIRV303 (§2.8).
   - **S0 validates safety for about one isoform (SIRV303), not §2.2's premise in general** (critique H1).
8. **A∘R3** behaves like A on BASE: the same C1e-C4 outcome on ≥ 5 of 6.
9. **The TOL plateau holds** (F7 does not fire) on ≥ 5 of 6.
10. **Verdict.**

    | verdict | probability | typical path |
    |---|---|---|
    | EFFECTIVE, default candidate | ~0.05 | C5 passes on ≥ 5 of 6 and C5.cur on both human substrates |
    | EFFECTIVE, opt-in | ~0.45 | P(s) 6/6, S0 active and safe, C5 fails as on dev |
    | KEEP OPT-IN | ~0.30 | C1e < 10% on one gorilla substrate, or S0 inactive, or SIRV303 lost |
    | REFUTE | ~0.20 | C1e misses on ≥ 2 substrates, or an S0 failure plus one genome failure |

    Dev was selected (§0 item 1), so a worse result is the expected direction of error.

## 9. What would falsify the design reasoning (reported whatever the verdict)

1. **F1, "A removes fragments, not complete transcripts."** Falsified if the `=` share among A's drops is ≥ ¼ of the
   `=` share among the **old NULL's** drops (the D-disjoint comparator; §9.1 of the draft, constant unchanged) on ≥ 2
   substrates. Dev ratios: .04 / .10 / .035.
2. **F2, "most of our excess `c` is co-emitted with an end-compatible container."** Falsified if A removes < 30% of
   BASE's `c` on ≥ 2 substrates. Dev: 53 / 45 / 53%.
3. **F3, "the dev effect transfers."** Falsified if the median C1 relative reduction over the six is < 10%. Dev:
   16.9-18.2%.
4. **F4, "A leaves k, m and n alone."** Falsified if |Δ(k+m+n)| > 5% of BASE's on any substrate. Dev: −2.2 / 0 /
   −2.3%.
5. **F5, "a real shorter isoform starts in its container's intron."** Falsified on SIRV if A removes the `=` of **any**
   true ISM other than SIRV303, in either count mode (row `S9.5_true_ism_eq_removed` and `sirv.eq_lost.*`). Every
   removal counts. The emitted end offset against the truth end is reported as the cause, descriptively, and never
   used as an exemption (critique H1). SIRV303's loss is the pre-registered exposure (§2.8) and does not falsify.
6. **F6, "(3) is selective per transcript."** Falsified if C1s or C3′ fails on ≥ 2 substrates.
7. **F7, "the TOL plateau."** Falsified if `c` removal at ρ 0.5 varies by > 25% relative across TOL 5-50
   ((max − min)/min of `d.dropped_class.c` over TOL 5, 10, 20, 50) on ≥ 2 substrates. Dev: 6-8%. A substrate whose
   grid was not run counts as not measured.
8. **F8, "the precision gain is p·d/(1−d)."** Predicted against observed per substrate. This is arithmetic, not a
   test.
9. **Not falsifiers:**
   - `m` not moving (by design and by the user's decision);
   - orphan `c` remaining (no emitted container; r1077's read ceiling);
   - H failing where §2.9 predicts it;
   - C5 failing (predicted, §8 item 5).

## 10. Not in this test

- **Lever B** (`--polish-ism-ratio 0.5` + A). It is not run anywhere here, including SIRV. SIRV's BASE and A rows
  become visible in this test, so any later prereg, Lever B's included, counts SIRV as a **second use**, not a fresh
  truth.
- ρ 0, other TOLs (descriptive grid only), the 3′-flush ISM ratio 1.0 plus retained 8, and the L3b orphan
  read-container test.
- Any `m`/`n` rule (closed: r1075, r1080; user decision 1).
- Any completion or extension of fragments (r1072, r1077). Single-read or StringTie-style admission (r1067-r1073,
  r1082, r1083). 5′ trimming (r1078, r1079).
- An annotation-free derivation of TOL (critique M1). It is recorded as a limitation. If A is ever flipped to
  default, the derivation must come first.

## 11. Order, stop rules, machine rules

Shorthand used in the commands:

```
F=/mnt/linuxdisk/tmp/rustle_figures/ct_frozen      BIN=/mnt/linuxdisk/tmp/rustle_figures/rt3_bin_frozen
C=/mnt/linuxdisk/tmp/rustle_figures/complete_arms  D=/mnt/linuxdisk/tmp/rustle_figures_dev/ct_sirv
P=/mnt/c/Users/jfris/Desktop/Rustle/docs/PREREG_complete_transcripts_2026-09-27.md
L="flock -w 900 /mnt/linuxdisk/tmp/rustle_heavy.lock timeout 600"
SC="python3 $F/complete_eval.py score --heldout --prereg $P --budget-s 540"   # exit 75 = run again
```

Output directories (one stem per substrate key and count mode, so they are kept apart):
- step 1: `--out $C/tables_step1`;
- SIRV: `--out $C/sirv/tables`;
- genome: `--out $C/tables`;
- descriptive: `--out $C/tables_desc/<chr21_22|r3|grid>`.

`run_arms.sh` takes the heavy lock itself for each step. **Never wrap it in `$L`**, because the inner flock would wait
on the outer one.

0. **Freeze: this file and Amendment 1.** Done when this file is saved.
1. **Held-out BASE and truth checks** (no A, NULL_S or NULL arm is scored).
   - (a) `cd $F && sha1sum -c SHA1SUMS`, and `sha1sum -c` in `$BIN`. `sha1sum tools/rustle_pipeline.sh` must equal
     d1816c35. Every `readlink -f rt_arms/<s>/<s>.BASE.gtf` and `.R3.gtf` must match `INPUTS_SHA1.tsv`.
     `git diff --quiet HEAD -- figures/` must hold, because the scorer imports the substrate registry from there. The
     REF sha1 check stops any change that reaches a reference anyway.
   - (b) IV5(a) and IV5(a′):
     `$L $SC --out $C/tables_step1 --sample <s> [--drop-contigs … | --contigs chr21,chr22] --arm BASE=<BASE.gtf>
     --base BASE`. The keys are those of §4.
   - (c) IV6: `iv6.multi_without_reads` = 0 on each BASE, read from the same runs.
   - (d) IV7:
     `$L python3 $F/complete_eval.py exposure --heldout --sample <s> [--drop-contigs …]` for V1-V6.
   - (e) IV8: the circularity record (§5.3).
   - Write **Amendment 2** with (a)-(e). Only the rows named in (b)-(c) are read from these runs. Any mismatch stops
     the work.
2. **SIRV (fresh).** Outputs go to `$C/sirv/`, with `TMPDIR=$C/sirv/tmp`.
   - `$L bash tools/rustle_pipeline.sh assemble --bam $D/bam/sirv_testis.bam --fasta $D/ref/sirv_E0.fa
     --out $C/sirv/dedup/BASE --bin $BIN --threads 4`
   - **IV1-SIRV.** Run the driver's `stage_assemble` expanded:
     `env RUSTLE_GTF_SECONDARY=1 RUSTLE_GTF_SECONDARY_AS_RATIO=0.98 RUSTLE_GTF_SECONDARY_AS_TABLE=$C/sirv/dedup/BASE.molecules.tsv
     $BIN/copy_assign --assemble-only --genome-wide --assembly-junctions strict --assembly-polish full
     --polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7
     --polish-retained-ratio 10 --gtf-tpm --bam … --fasta … --out $C/sirv/iv1/BASE`, under `$L`.
     - Its GTF must be `cmp`-identical to `dedup/BASE.gtf`.
     - Only then is the same command run with `--keep-coordinate-duplicates` and `--out $C/sirv/kcd/BASE`.
     - A difference stops SIRV until explained. An explanation that leaves transcript and exon lines untouched is
       recorded as an amendment, and the check is repeated on those lines.
   - For m in dedup, kcd: `$F/run_arms.sh sirv_E0_human_testis $C/sirv/$m/BASE.gtf $C/sirv/$m/arms` (binding seed).
   - For m in dedup, kcd: `$L $SC --out $C/sirv/tables --sirv --arm BASE=… --arm A=…/A.gtf --arm NULL_S=…/NULL_S.gtf --arm NULL=…/NULL.gtf
     --base BASE --judged A --null-s NULL_S --null NULL --count-mode $m`.
   - Check the IV3/IV6 rows first.
3. **Record S0's rows** (Amendment 3, or Outcome part 1). Nothing may change after them. An S0 failure does not stop
   the genome stage; it enters §7.
4. **Genome substrates, smallest first:** human_testis, chimp_PTR, orangutan_PPY, gorilla_KB3781, gorilla_OR6737,
   human_A119b. Per sample:
   - re-check the BASE sha1;
   - `$F/run_arms.sh <s> rt_arms/<s>/<s>.BASE.gtf $C/<s>/arms`;
   - `$L $SC --out $C/tables --sample <s> [--drop-contigs chr16,chr20,chr21,chr22 | NC_073244.2] --arm BASE=… --arm A=…
     --arm NULL_S=… --arm NULL=… --base BASE --judged A --null-s NULL_S --null NULL --count-mode dedup`;
   - check the IV3 rows before any clause row is used.
5. **Verdict** (§7) from the TSVs. The Outcome is appended here, and every number becomes a register row
   (`docs/NEGATIVE_RESULTS_REGISTER.md` if refuted).
6. **Descriptive, after the Outcome's verdict line is written:**
   - chr21-chr22 (`--contigs chr21,chr22`, same arms);
   - A∘R3 (`$L python3 $F/compat_collapse.py R3.gtf …/A_R3.gtf --tol 10 --rho 0.5 --mode drop`, then
     `$L $SC --out $C/tables_desc/r3 … --arm BASE=R3.gtf --arm A=…/A_R3.gtf --base BASE --judged A`);
   - the grid: `compat_collapse.py --mode drop` per (TOL, ρ), all arms of a substrate in one
     `$L $SC --out $C/tables_desc/grid … --arm BASE=… --arm T<tol>_R<rho>=… --base BASE` call (resumable). TOL 5, 20
     and 50 at ρ 0.5 come first (F7), then the rest.

   Anything not run is reported as not measured.

**Stop rules.**
- After this file, nothing changes: not the rule, TOL, ρ, NULL_S, a clause, D1, D2 or the scorer. The one exception
  is a bug fix, recorded as an amendment, which re-runs every A, NULL_S and NULL arm.
- D1 and D2 may be changed only as §0 item 3 says.
- No variant (ρ 0, another TOL, Lever B, another NULL seed) is ever substituted for §1 or §3.1.
- These stop the work:
  - a sha1 mismatch (BASE, R3, binaries, driver, frozen scripts, references);
  - a failed freeze check (a) or (a′);
  - a non-zero IV3 or IV6 row;
  - any unexplained IV difference.
- Nobody opens `ct_frozen/sealed/`. It may be deleted after the Outcome, because it is this test's own scratch.

**Machine rules.**
- One heavy process at a time, in the foreground, under `$L`. This covers each gffcompare, each post-processor call
  on a genome-wide GTF (the A119b BASE is 334 MB), each NULL draw and each assembly.
- Any step that can exceed 600 s is split: the scorer caches gffcompare per arm (sha1 key) and exits 75 on budget,
  and the same command is run again.
- `TMPDIR` and all outputs go under `/mnt/linuxdisk`. Never `pkill -f`; kill by PID after
  `readlink /proc/<pid>/cwd`.
- No Rust edits and no cargo build. Binaries come only from `rt3_bin_frozen/`.
- **Disk (~84 GB free).**
  - Keep the drop and NULL lists, the logs, the score TSV/JSON, and the gffcompare `.tmap`, `.refmap` and `.stats`.
  - Delete the arm GTF copies once the Outcome is written. They regenerate byte-identically from BASE plus the frozen
    scripts.
  - Delete only this test's scratch.

---

## Amendments

### Amendment 1 (2026-09-27): the freeze; written before any held-out or SIRV number of A, NULL_S or NULL

1. **USER_REUSE_ACCEPTANCE: accepted 2026-09-27, via the user's two answers.**
   - (1) The target is class `c`, not `m`.
   - (2) The held-out set is the SIRV E0 spike-ins of human_testis plus the six samples minus dev contigs, recorded as
     the fourth reuse of the six.

   chr21 and chr22 leave V1 (critique H3), in line with (2).
2. **Frozen scripts.** `/mnt/linuxdisk/tmp/rustle_figures/ct_frozen/`: `sha1sum -c SHA1SUMS` passes 7/7, and the
   `SHA1SUMS` file itself is 9c547f5d62b544bcd474b0f3bf073c100efb3a56.

   | file | sha1 |
   |---|---|
   | compat_collapse.py (arm A; = `ct_prereg` e14c5646) | e14c5646b133acff731d1335f877d9ad2046fc3a |
   | complete_null.py (NULL_S binding, NULL reported) | 9740f092eb2dec93f5e31683a7f3bed7bc6e8438 |
   | **complete_eval.py (the scorer)** | **6461144397294350fe7727e2bf1a1f626032ed03** |
   | run_arms.sh | bcf50eb7d784d76dc41fe33b4b1c6a60585084de |
   | REF_SHA1.tsv | 06de71133ad0cadfc17d9c0641c67a248ba3192a |
   | INPUTS_SHA1.tsv | 1c79080836b14c50ac0a33ba4e034dfb3bbba8cf |
   | FREEZE_LOG.md | 6d5e3f9dde8b4faee267aabcc593b060a743c63d |

   Two scorer versions were superseded **before any scoring of our held-out or SIRV output** (`FREEZE_LOG.md`):
   - v0 c894943ac55615b2a449f147ab1d8a510af2c3a3, which had a typo in the StringTie SIRV sha1 constant;
   - v1 a97ecbfb92f526d25a8d8b6c8c1b210bc1bfe621, which lacked the SIRV clause rows and the tmap-based S1.

   Every dev table was re-scored with v2 (6461144…).
3. **Binaries, tools and the driver.**
   - `rt3_bin_frozen/copy_assign` f480847abefefc7ff7997c528c2fd4dce5bda57b;
   - `as_table` f9556f36eb76da65a9643e93912ac72c889df216;
   - gffcompare 0.12.10 `/home/juanfra/miniforge3/bin/gffcompare` 32cbdd6491c7f36dfcd8495cfc0e0fef7cdbeca4;
   - `bench/mechanism/readthrough_eval.py` f6b99dcc1a3f972427cf67bc3de3b846dd5c5f30 (IV5 reference);
   - driver `tools/rustle_pipeline.sh` d1816c3571c53d9fc3df6a63be47cc39a09dca60 (git 634c9551, unmodified).
4. **Inputs** (`INPUTS_SHA1.tsv`; hashes only).

   | sample | SHA1_BASE (`runs/<s>/<s>.gtf`) | SHA1_R3 (`rt_arms/<s>/<s>.R3.gtf`) |
   |---|---|---|
   | human_A119b | 631c9f114b3728c51405a6e4f019f20cee43f12e | d6afe1ea187dccc320925d0c99c1add302f3d7e9 |
   | human_testis | 50f239d6da21355f3d9723c495117980e1e9c25f | 60d4161fcd0e857f9f1b2eb26440a74a04089ed4 |
   | gorilla_OR6737 | 3a8c6410d8176d682d8112a8a32a19a7de5b1983 | 3e9735325437cb18edade20fb3acc94fc69d05cc |
   | gorilla_KB3781 | f1022f4ba65448c8f4c525c7e100e206b0174c57 | 6f13236a531a7aacbb8f74ad11c860290952b670 |
   | chimp_PTR | dead746f1aa3f7a4d79df5327c1884bbe3cf5873 | c1e36fc99803211da8141a1611d890f6a23443d1 |
   | orangutan_PPY | b013a74a22f923fb1e0aff81643a6265e7316bf5 | 16cc4c0f43f1ea996b0c93be2869e4348222f52d |

   - **Dev BASEs:** hsa20 6fa4dbaa6601c1cdcbf45ec9332a8bb1f1cf7443, ggo44 41815166d888e63d397fa8135e264dfb88ec1f1a,
     hsa16 46a24de7f6efdd8759acfd4596166eeabbef2d62.
   - **SIRV:**
     - truth `sirv_C.gtf` 81051062f744db0de12ca301ac4cc3fbb9841390;
     - `sirv_E0.fa` 4d7fbf49936b3281d854f0fe247f829698f55559;
     - `sirv_testis.bam` 1c5a0b51c5e3c22a6856b18c67117ead2d44baa5;
     - `molecules.tsv` bc1d0c9e4c26abd07a36bbaf45ac3ac4256aca8c;
     - `truth_support.tsv` 9c4fee0a8526e4208f23cb6f35a5726896ffc3bf;
     - tool GTFs: StringTie ece3289d586cb512103ecf3afcfbcb24b1cd4c52, FLAIR
       74fb52e388cd558416d33e7a0e116c30b6ffdec9, IsoSeq 9cdff97941831ba013545769ebac972df7d6a366.
5. **Restricted references** (`REF_SHA1.tsv`; the scorer stops on any mismatch; critique M2).

   | key | sha1 | role |
   |---|---|---|
   | human_A119b:chr20 / :chr16 | f839174dbf4a1cdcd08bdb0aaeb24085ba4480fd / 89e140c26665314c24e2df3627fbb8df09589ea9 | dev |
   | gorilla_OR6737:NC_073244.2 | 8c27a18a89bf3cfc3490f2d50efa15aca28d455c | dev |
   | sirv_E0_human_testis:SIRV1-7 | 81051062f744db0de12ca301ac4cc3fbb9841390 | S0 |
   | human_A119b:annotated_minus_chr16_chr20_chr21_chr22 | c89dfe1addd9a0c3417a8e8ffa37a25dbc7c41c7 | V1 |
   | human_A119b:chr21-chr22 | 5d37e1d7b6246c3ea82b508755d9a50aebfe46d3 | reported; IV5(a′) |
   | human_A119b:annotated_minus_chr16_chr20 | fe75a71496e518bf947b7bdd2ac51a8abd5b3503 | freeze check (a) only |
   | human_testis:annotated | 8ed998a1fb8cd52f14b1ad048c7c0db6d11a2854 | V2 |
   | gorilla_OR6737:annotated_minus_NC_073244.2 | 9271749f72dd362cba91907f0eba10e0bf039ab8 | V3 |
   | gorilla_KB3781:annotated | a7e0e06e8ea4f09a07bcc7d68031bdc1926b67ab | V4 |
   | chimp_PTR:annotated | f4b58b508d2d360fea5462b4852fad7e0ff884b0 | V5 |
   | orangutan_PPY:annotated | 35a550520091602e7d174500932a8dd84bdc7819 | V6 |

   Sources:
   - human `A119b.chm13v2.0_RefSeq.gtf.gz` 41f7f281…;
   - gorilla `GGO.GCF_029281585.2_RefSeq.gtf.gz` dd680eae…;
   - chimp and orangutan `assembly/<s>/annotation.gff_to_gtf.gtf` (sha1 as in the table).
6. **Dev gates and the dev clause table under this file's clauses** (never in the verdict; 5 NULL_S seeds × 3 contigs;
   §2.4-§2.5):
   - IV1, IV2, IV3, IV4 (selftest 29/29), IV5, IV6, the SIRV technical gate and the guards pass (§4).

   | contig | C1e | C1s | C2 / C3 / C4 | C3′ | Q | C2.cur / C4.cur | **C5** | C5.cur (binding) | H | H_NS | C1½NS / C3′¼NS / C1_NULL |
   |---|---|---|---|---|---|---|---|---|---|---|---|
   | chr20 | −18.2% ✓ | 5/5 | ✓ | 5/5 | 5/5 | ✓ | **0/5** (7 vs 2-4) | 6 vs 3 ✗ | ✓ | 0/5 | 1/5 / 5/5 / 5/5 |
   | NC_073244.2 | −16.9% ✓ | 5/5 | ✓ | 5/5 | 5/5 | ✓ | **0/5** (7 vs 2-5) | 0 vs 0 ✓ | ✓ | 0/5 | 0/5 / 0/5 / 5/5 |
   | chr16 | −18.2% ✓ | 5/5 | ✓ | 5/5 | 5/5 | ✓ | **1/5** (9 vs 4-9) | 9 vs 8 ✗ | ✓ | 0/5 | 0/5 / 4/5 / 5/5 |

   P(s) passes on all three contigs at every seed, and T fails. The dev pattern is therefore the
   "EFFECTIVE, opt-in" shape. Dev was selected, and dev never enters the verdict.
7. **Seal** (critique H2). All 12 files of `ct_sirv/smoke/` were moved unread to
   `ct_frozen/sealed/ct_sirv_smoke/`.
   - Manifest `sealed/ct_sirv_smoke.SHA1SUMS` = 414f98dbe4e9197111f43265b273f3ba6ae92d0d.
   - `smoke.gtf` = 34d222e4a049a4fb47dff7bfa5348d05ed129c13.
   - No other output of ours on SIRV exists, and `ct_sirv/tmp` is empty.
8. **Not yet run, by design (step 1, recorded in Amendment 2):**
   - freeze check (a) and V1 additivity (a′);
   - IV6 on the held-out BASEs;
   - held-out exposure (IV7);
   - the circularity record (IV8).
9. **Scorer realisation of this file's names.** The scorer's rows keep their frozen names:
   - `C1` = C1½NS (reported);
   - `C3'` = C3′¼NS (reported);
   - `C1_NULL` (reported);
   - `C2`, `C3`, `C4`, `Q`, `C5`, `C2.curated`, `C4.curated` and `C5.curated` are this file's clauses of the same
     name.

   C1e, C1s and C3′ are computed from the metric rows as §6.1 states.

---

## Outcome (2026-09-27): withdrawn on dev; held-out not spent

**Verdict: withdrawn before step 1.** Arm A was not run on any held-out or SIRV substrate.
- No held-out or SIRV number of A, NULL_S or NULL exists. Amendment 2 was never written, step 1 was never run, and
  `complete_arms/` was never created.
- Nobody opened `ct_frozen/sealed/`. **The SIRV E0 truth remains sealed and unspent** for our output: a later
  prereg may use it as a first use. The truth-side facts of §2.8 and the three tools' SIRV rows were seen before
  this file, as §0 records.
- The six genome substrates carry no verdict from this file. The fourth reuse accepted in Amendment 1 was not
  spent.
- The text above this section is sha1 76791696365211d5f583584da161e631ae14c212 (948 lines). None of it was changed.
- The §8 predictions stay unscored.
- All numbers below are **dev design evidence**, on the contigs A was selected on (§0 item 1): human A119b chr20 and
  chr16, gorilla OR6737 NC_073244.2. Species are never pooled.

### O.1 Clause placement: the critique's, not D1/D2
The main session **rejected D1 and D2** (§0 item 3), following ct2_review R1 and R3. The binding placement is the
critique's:
- **C5 is judged inside P(s)** (critique B2), not in a separate tier T.
- **C1s and C3′ take the R3 form**: A's drops D against NULL_S's stratum-matched draws **outside** D (`in_D` = 0),
  with classes from BASE's tmap. The constants are kept:
  - C1s ⇔ 2·part_R·n_D < part_D·n_R;
  - C3′ ⇔ 4·eq_D·n_R ≤ eq_R·n_D;
  - "not measured" if n_R = 0.
- The direction-only forms of D1 are reported only.

### O.2 A on dev under this placement: REFUTE-shaped on 3 of 3 contigs
**C5 fails 14 of 15** (5 NULL_S seeds × 3 contigs). The only pass is a tie, 9 = 9 on chr16 seed 2.

| contig | genes losing their only multi-exon `=`/`c` query: A vs NULL_S (binding seed) | NULL_S, 5-seed range / mean | C5 |
|---|---|---|---|
| chr20 | 7 (1.46% of 480) vs 3 | 2-4 / 3.2 | ✗ (0/5) |
| NC_073244.2 | 7 (0.79% of 885) vs 5 | 2-5 / 3.4 | ✗ (0/5) |
| chr16 | 9 (1.15% of 780) vs 7 | 4-9 / 6.6 | ✗ (1/5, a tie) |

Every other judged clause passes at the binding seed (chr20 / NC_073244.2 / chr16):
- **C1e**: the partial share falls by −18.2 / −16.9 / −18.2%;
- **C1s-R3**: the partial share of D against R is .671 vs .283, .793 vs .227 and .628 vs .278 (margins 1.19 / 1.75 /
  1.13). It passes on 5/5, 5/5 and 4/5 seeds;
- **C3′-R3**: `=` among the drops is 2/243 vs 10/145, 3/121 vs 8/44 and 3/500 vs 10/299 (5/5 seeds each);
- **C2, C3, C4 and Q** pass (matched chains −0.19% on each contig).

**Reading.** A is REFUTE-shaped on dev: P(s) fails on 3/3 contigs, and §7 refutes at 2 failing substrates. Dev is the
optimistic edge. Running the held-out stage would have spent the six substrates and SIRV on a predicted failure, so
the test was withdrawn.

**The `k`-container caveat does not rescue A.** In 4 of the 23 lost genes the container ⊇ the reference (CFAP61,
MATN4, LOC101151087, LOC124907834). Excluding them from A's count alone still leaves 5 / 6 / 8 losses, against the
NULL_S means of 3.2 / 3.4 / 6.6.

### O.3 What the lost genes are (`ct3_anatomy`; the annotation is used only as a label)
- **Pairs.** There are 23 genes and 36 F→C pairs (11 / 9 / 16).
  - 35 of the 36 F are **5′ fragments**. C is F's chain plus 1-28 upstream exons (median 3).
  - C never differs from F inside F's span.
  - For F's reference, C is `j` in 30 pairs, `k` in 4 and `m` in 2. The gene keeps its junctions but loses the
    intron-correct label.
- **Depth.** Lost genes are shallow: median gene depth 31 reads, against 178 for SAFE pairs (AUC .75). NULL_S's
  read-bin strata absorb most of that signal.
- **Attach point.** Every attach-point feature has AUC ≤ .67: link share, truncation share, 5′ sharpness and the
  overhang.
- **Shared with NULL_S.** 16 of the 23 genes are also lost under at least one NULL_S seed. C5 therefore asks A to be
  selective **within its own D**.

### O.4 The guard search: 14 reads-only guards, none passes; CLOSED on these three contigs
Sources: `ct3_anatomy`, `ct3_variants` and `ct3_review`.

**Guards.** Each guard is a reads-only conjunct added to A's (1)-(4). A keeps a drop only when some container passes
(1)-(4) and the guard:
- `dom`, `dom_lenient`, `plur`, `dom_attach`, `tail`, `dom_tail`, `dom_stop1`, `dom_tail_stop1` (ct3_anatomy);
- `jsup`, `domC`, `tss`, `corrob`, `dom|corrob`, `dom+stop1+domC` (ct3_variants).

**Inputs.** Every guard reads only BASE.gtf and the dedup primary reads. The rule code reads no annotation
(ct3_review §2).

**The search was in-sample and label-supervised.**
- Every guard was designed after the 23 lost genes had been seen, on the same three contigs.
- The selection used annotation labels on 36 pairs.
- 14 guards × 3 contigs × 5 frozen NULL_S seeds = 210 clause evaluations. Each guard has its own NULL_S, matched to
  its own D′.

**Result.** None passes every judged clause on all three contigs at the binding seed. Seed-pooled (in expectation),
none passes more than 2 of 3.
- **Closest: `dom+stop1+domC` (V6).** It passes chr20 and gorilla on 5/5 seeds. It fails chr16 on C1e (−9.1%) and on
  C1s-R3 (0/5). It keeps only 58 / 67 / 52% of A's `c` removal.
- **`dom_tail_stop1`** passes all three contigs only at non-binding seed 3.

**Why they fail.** The two walls sit on different contigs.
- **C5 binds on chr20**, where NULL_S loses only 0-3 genes.
- **C1s-R3 binds on chr16**, where A's own margin is 1.13.
- A guard that protects enough fragments to reach C5 parity returns those partials to NULL_S's pool outside D. That
  raises part_R and breaks C1s.

**Annotation oracle (not a rule).**
- The oracle drops A minus the 36 LOST F: 11 / 9 / 16 transcripts, 5-7% of A's drops.
- It passes every judged clause on all three contigs at the binding seed, and in 14 of 15 seed × contig cells.
- It keeps **93 / 93 / 95%** of A's `c` removal (P −16.8 / −15.6 / −17.3%).
- Pair protection:
  - the best guard, `dom|corrob`, protects 28% of LOST pairs at a cost of 2.7% of SAFE ones;
  - the oracle protects 100% at 0%.
- **The clause set is satisfiable by a subset of A. The limit is the reads' selectivity for the at-risk 5′
  fragments**, not a conflict between C5 and C1.
- chr16's C1s-R3 is tight even for the oracle: 1.02 at the binding seed, failing at seed 2. A later human prereg for
  any rule near A should expect C1s-R3 to be close to a coin flip. This is not a reason to relax the ½.

**Closed.** The guard search on human A119b chr16/chr20 and gorilla OR6737 NC_073244.2 is **closed**. A new guard
idea needs a fresh dev substrate outside V1-V6 (hold a substrate back).

### O.5 What shipped instead (opt-in; this file's verdict was not needed for it)
**`copy_assign --polish-subchain off|tag|drop`**, default `off`.
- **Status.** It was implemented after this Outcome's decision and was uncommitted when this was written. The tested
  binary is `ct_bin_frozen/copy_assign` 6bbd6442.
- **Driver.** `RUSTLE_POLISH_SUBCHAIN=tag|drop` adds the flag to the `assemble` stage. Unset gives the same command
  as before.
- **Rule.** It is an exact port of §1 at TOL 10 and ρ ½, with none of the O.4 guards. The ½-reads support condition
  (4) is kept.
- **Definition, in the shipped wording.** It flags *"an end-compatible contiguous sub-chain of a longer emitted
  transcript of the same locus with >= 1/2 its reads"*.
  - The tag does not say "incomplete". On dev, 63% of tagged transcripts are `c`, 30% are `j`, and 8 of 864 are `=`.
- **Vocabulary changed from §1 item 5** (ct2_review L2).
  - Old: `completeness "5p_fragment"|"3p_fragment"|"internal_fragment"; fragment_of "<x>";`.
  - New: `subchain_of "<x>"; subchain_missing "5p"|"3p"|"both";`. Here `subchain_missing` names the ends of x's
    chain that y lacks, in transcript orientation (a `-` strand swaps them).
  - The prototype's `5p_fragment` (m3 = 0: 5′-truncated, which is SQANTI3's `3prime_fragment`) becomes
    `subchain_missing "5p"`.
  - Dev counts, seeded, 864 tagged: 5p 804, 3p 54, both 6.
- **Verification** (`ct4_impl`, independently re-checked by `ct4_review`):
  - `off`, whether the flag is omitted or explicit, is byte-identical to the f480847a products: 394/394 over 56 dev
    runs, plus 120/120;
  - `tag` equals e14c5646 `--mode annotate` with the vocabulary mapped, 18/18 cells;
  - `drop` equals e14c5646 `--mode drop` in 18/18 cells. TPM is renormalised over the kept set;
  - streaming equals `--materialize-reads` once `matched_reads` is masked;
  - `--genome-wide` equals the per-contig runs concatenated;
  - `cargo test --release`: 917 passed, 0 failed, 13 ignored (8 new);
  - no GTF consumer parses the new attributes.

**Drop trade-off (dev, in-sample: A was selected on these contigs; no held-out test; species never pooled).**

| contig | `c` share of multi-exon queries | partial share (rel.) | `=` precision | `=` queries lost (matched chains) | genes losing their only `=`/`c`: drop vs NULL_S (binding; 5-seed mean) | `j` dropped | loci / mono-exon lost |
|---|---|---|---|---|---|---|---|
| human chr20 | 5.7 → 2.8% | −18.2% | +1.02 pt | 2 (1,059 → 1,057) | 7 (1.46%) vs 3; 3.2 | 74 | 0 / 0 |
| human chr16 | 6.3 → 3.1% | −18.2% | +1.07 pt | 3 (1,605 → 1,602) | 9 (1.15%) vs 7; 6.6 | 171 | 0 / 0 |
| gorilla NC_073244.2 | 5.4 → 3.1% | −16.9% | +1.23 pt | 3 (1,595 → 1,592) | 7 (0.79%) vs 5; 3.4 | 21 | 0 / 0 |

- **The precision gain is mostly the smaller denominator.** A matched random drop (NULL_S) gets 78-91% of it, and
  p·d/(1−d) predicts +1.06 / +1.10 / +1.31 pt (chr20 / chr16 / gorilla).
- **Unmeasured reach of `drop`.**
  - It removes the families-stage representative of 15-19 loci per contig.
  - It shrinks 8-11 gene spans per contig by at most 10 bp, and the driver's `flag` stage reads those spans.
  - `tag` is byte-neutral for every consumer.
- **Not a default candidate.** A default flip needs its own held-out prereg on untouched substrates. SIRV E0 is still
  available for it.

### O.6 Bookkeeping
- **Register.** Draft rows r1128 (the lever on dev), r1129 (the guard search, closed) and r1130 (what ships).
  Numbering is as drafted, after r1127.
- **L4 hygiene (§7) does not apply.** A was not found EFFECTIVE, so r866, r1066 and r1078 stand as written.
- **Reports** (session scratchpad `figs/`): `ct2_review.md`, `ct3_anatomy.md`, `ct3_variants.md`, `ct3_review.md`,
  `ct4_impl.md` and `ct4_review.md`.
- **Scripts and outputs.** `/mnt/linuxdisk/tmp/rustle_figures_dev/`:
  - `ct3_anatomy/`, `ct3_variants/`, `ct3_review/` (oracle: `oracle.sh` a6bedc5f, `oracle_eval.py` a0c236e4);
  - `ct4_impl/` (patches `subchain_copy_assign.patch` 866a0f16 and `subchain_driver.patch` fc4b9787);
  - `ct4_review/`.
- **Frozen.** `ct_frozen/` is kept unchanged. `sealed/` stays unread and may be deleted unread (§11); SIRV stays
  fresh either way.
