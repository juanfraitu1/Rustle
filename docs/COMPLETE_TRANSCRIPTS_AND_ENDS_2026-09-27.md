# Complete transcripts and read-proven transcript ends: findings (2026-09-27)

This document consolidates one day's work for the user and the advisor. It adds no new measurement.
- **Sources.** Every number is quoted from a source (key in brackets, listed in §10), or is marked "computed here"
  from one.
- **Species** are never pooled.
- **Dev** means the contigs the rules were designed on: human A119b chr20 and chr16, and gorilla OR6737 NC_073244.2.
  Dev numbers are in-sample.
- **Held-out** means the six genome samples minus the dev contigs and chr21/chr22, plus the human_testis SIRV E0
  spike-ins.
- **Scoring** is gffcompare 0.12.10 on multi-exon queries.

## 1. The question

The user's goal: *"ensure the assembler part of the pipeline tries as best as possible to emit complete transcripts,
meaning that when running gffcompare afterwards there are few partial transcripts like m categories"* [CT header].

| class | meaning (query vs its best reference transcript) | partial? |
|---|---|---|
| `=` | identical intron chain | no |
| `c` | contained: the query's chain is a compatible part of the reference's, usually a fragment missing exons at one end | yes |
| `k` | reverse containment: the query adds exons around the reference | yes |
| `m` / `n` | retained intron(s); every other intron matches (`m`) or not all do (`n`) | yes |
| `j` | shares at least one junction, otherwise differs (a candidate novel isoform) | no |

**Definitions** [HO §6]:
- class shares are over multi-exon queries;
- partial = `c`+`k`+`m`+`n`;
- chains = reference intron chains matched (`.stats`);
- chain precision = `=` queries / multi-exon queries.

**Why the target became `c`, not `m`** [CT header, decision 1]. The user made this call after seeing the dev data:
- we already emitted the fewest `m` of the four assemblers (2.3 / 1.7%, against 3.9-6.9% for StringTie, FLAIR and
  IsoSeq);
- the `m` rule space is closed (r1075 adopted, r1080 refuted);
- `k`, `m` and `n` stay inside the partial share, so they are still measured.

## 2. Where we stand

| dev (in-sample) [HO §2] | chr20 `=` | `c` | `m` | chains | NC_073244.2 `=` | `c` | `m` | chains |
|---|---|---|---|---|---|---|---|---|
| ours (BASE) | .210 | .057 | .023 | 1,059 | .409 | .054 | .017 | 1,595 |
| StringTie | .168 | .008 | .069 | 861 | .370 | .027 | .051 | 1,374 |
| FLAIR | .077 | .012 | .041 | 1,026 | .238 | .019 | .039 | 1,393 |
| IsoSeq | .048 | .031 | .042 | 1,253 | .142 | .043 | .057 | 1,655 |

**Held-out, with the lab's tool outputs** [HO-rows]. Every row here is **descriptive**: our BASE class shares had been
seen in figures 1-2 before the prereg was written [HO §0, §5].

| sample | arm | multi-exon | `=` | `c` | `k` | `m` | `n` | `j` | other | partial | chains |
|---|---|---|---|---|---|---|---|---|---|---|---|
| human A119b | ours | 195,666 | .178 | .074 | .029 | .019 | .026 | .512 | .162 | .148 | 34,752 |
| | StringTie | 190,865 | .150 | .011 | .047 | .053 | .067 | .427 | .246 | .177 | 28,616 |
| | FLAIR | 567,538 | .063 | .015 | .030 | .037 | .052 | .551 | .253 | .134 | 35,367 |
| | IsoSeq | 1,768,660 | .043 | .036 | .019 | .034 | .060 | .528 | .279 | .149 | 42,611 |
| gorilla OR6737 | ours | 69,446 | .350 | .084 | .032 | .014 | .009 | .476 | .036 | .138 | 24,277 |
| | StringTie | 63,777 | .344 | .027 | .043 | .040 | .030 | .462 | .054 | .141 | 21,946 |
| | FLAIR | 132,646 | .175 | .022 | .038 | .029 | .020 | .686 | .030 | .109 | 22,958 |
| | IsoSeq | 423,072 | .120 | .061 | .026 | .042 | .046 | .615 | .091 | .175 | 27,771 |

**What the held-out table shows** [r1137 for `m` and `c`; the other bullets are read here from HO-rows]:
- **`m`** is the lowest of the four arms on both samples, and every arm has `m` = 0 on SIRV.
- **`c`** is the highest of the four arms; it is our one excess partial class.
- **Total partial share** is middling: FLAIR's is lower on both samples.
- **Chain precision** (the `=` column) is the highest of the four arms.
- **Matched chains** trail FLAIR and IsoSeq on A119b, and IsoSeq on OR6737.
- **Output size differs.** FLAIR and IsoSeq emit 2.9× and 9.0× our multi-exon count on A119b, and 1.9× and 6.1× on
  OR6737 (computed here).

**SIRV E0** (dedup; 61 multi-exon truth isoforms) [HO-rows]:

| arm | multi-exon queries | `=` | `c` | `m` | truth isoforms recovered as `=` |
|---|---|---|---|---|---|
| ours | 61 | .787 | .164 | 0 | 48 |
| StringTie | 44 | .886 | .023 | 0 | 39 |
| FLAIR | 64 | .766 | 0 | 0 | 49 |
| IsoSeq | 109 | .569 | .349 | 0 | 47 |

**Pre-registered predictions that missed** (checked here against [HO-rows]; the Outcome did not score them, and none
enters a verdict) [HO §8]:
- "lowest `other` share" missed on OR6737 (FLAIR .030 vs .036);
- "chains IsoSeq > ours > FLAIR" missed on A119b (FLAIR 35,367 vs 34,752);
- "TES50 genes IsoSeq > ours > FLAIR" missed on OR6737 (FLAIR 10,722 vs 10,291);
- SIRV "kcd ≥ dedup" missed (45 < 48);
- "drop cuts `c` 40-55% everywhere" missed on chimp (34%), testis (15%) and SIRV (77%);
- "drop costs chains −0.1 to −0.3%" missed narrowly on KB3781 (−0.096%) and on SIRV (0);
- "drop costs 0.5-1.5% of genes their every `=`/`c` query" missed on testis (0.4%) and KB3781 (0.3%);
- "rescue adds 0.5-1.5% on both human samples" missed on testis;
- "pas-end moves 1-4%" missed on testis (0.3%);
- "J2 passes on both human samples with the largest relative margins" missed: testis ties, and A119b's margin
  (172 / 883 = 19.5%) is below KB3781's (175 / 867 = 20.2%); the predicted verdict EFFECTIVE did not happen;
- "J2 margins ≤ 10% off the human samples" missed: they are 12-20%.

## 3. What was tried for `c`, and why the default did not change

**Where `c` comes from** (dev, chr20 / gorilla / chr16) [ct_attr; CT §2.1].
- 85 / 92 / 77% of `c` are 5′-truncated.
- For 74-81% of them, no full chain has ≥ 2 reads: this is r1077's read ceiling. Only 2 / 1 / 8 have a ≥ 2-read full
  chain that we fail to emit.
- **The excess is co-emission.** 206 / 119 / 406 of the 286 / 212 / 552 `c` have an emitted container at the same
  locus. We emit 29.1 / 21.0 contained transcripts per 1,000 queries, against StringTie's 5.5 / 10.2 (chr20 / gorilla).

**End compatibility is what separates a fragment from a real shorter isoform** [ct_levers; CT §2.2].
- Co-emitted `c` overhang the container's exons by a median of 0 bp.
- `=` transcripts that have a container overhang it at the 5′ end by a median of 80 / 33 / 55 bp; only 6 / 10 / 8 lie
  within 10 bp at both ends.
- In other words, a real shorter isoform starts in the container's intron, while a truncated read starts inside its
  exon.
- Of about 40 lever variants, the existing dials do not help. A lower ISM ratio adds `c`, and a ratio ≥ 1.5 costs
  0.8-12% of chains.

**Arm A, the compatible-containment collapse** [CT §1]. It drops y when:
- y is an exact contiguous block of a longer x at the same locus;
- y's ends lie within 10 bp of x's corresponding exons;
- reads(x) ≥ ½·reads(y).

It keeps every locus and junction, and completes nothing: for only 14-22% of the removed `c` is the surviving container
an `=` of the same reference.

**It failed a gene-level clause** [CT O.2, O.3; r1128].
- **The null.** NULL_S drops at random from the same sub-chain pool, stratum-matched, without the end condition.
- **The clause.** C5 requires that A cost no more reference genes their only multi-exon `=`/`c` query than NULL_S does.
- **The failure.** C5 fails in **14 of 15** seed × contig cells: 7 / 7 / 9 genes against 3 / 5 / 7 at the binding
  seed (5-seed means 3.2 / 3.4 / 6.6). The one pass is a tie.
- **Everything else passes at the binding seed** (C1s-R3 passes 4 of 5 seeds on chr16). The partial share falls
  18.2 / 16.9 / 18.2% (relative), and chains fall 0.19%.
- **The precision gain is mostly the smaller denominator.** p·d/(1−d) predicts +1.06 / +1.31 / +1.10 pt against the
  observed +1.02 / +1.23 / +1.07, and NULL_S gets 78-91% of it.
- **What is lost.** 35 of the 36 lost pairs drop a 5′ fragment. For the dropped transcript's reference, the surviving
  container is `j` in 30 pairs, `k` in 4 and `m` in 2. The gene keeps its junctions but loses its intron-correct label.
  The lost genes are shallow: median depth 31 reads, against 178 for SAFE pairs.
- **Status.** The held-out test was withdrawn on dev and never spent.

**Fourteen reads-only guards, none passes** [CT O.4; r1129; ct3_anatomy; ct3_variants; ct3_review].
- None passes all three contigs at the binding seed, and seed-pooled none passes more than 2 of 3.
- An annotation oracle (A minus the 36 lost fragments; not a rule) passes all three and keeps 93 / 93 / 95% of A's `c`
  removal.
- The best guard protects 28% of lost pairs at 2.7% of safe ones; the oracle protects 100% at 0%.
- **So the clauses can be satisfied, and the limit is that the reads do not single out the at-risk 5′ fragments.**
- The search is closed on these contigs. A new guard needs a fresh dev substrate.

## 4. What shipped (all opt-in; off by default, off byte-identical) and the held-out numbers

```
RUSTLE_POLISH_SUBCHAIN=off|tag|drop  RUSTLE_POLISH_TSS=off|tag|rescue|split  RUSTLE_POLISH_TES=off|tag|pas-end \
  tools/rustle_pipeline.sh assemble --bam B --fasta G --out PREFIX
```

- **Driver.** An unset variable leaves the command unchanged, and an invalid value exits 2
  (`tools/rustle_pipeline.sh:94-114`). [AP]'s remark that the sub-chain value is not validated is stale.
- **Direct flags.** `copy_assign --polish-subchain|--polish-tss|--polish-tes`.
- **Order.** The sub-chain decision is taken on the final set, after any rescue, split or 3′ move [AP].
- **Verdicts.** Only pas-end was judged held-out; the drop and rescue rows are descriptive [HO header, §7].

| `--polish-subchain drop` [HO-rows; r1138] | `c` BASE → drop | partial | chains Δ (% of BASE) | genes losing every `=`/`c` (share) | TES50 genes lost |
|---|---|---|---|---|---|
| human A119b | .074 → .037 | .148 → .115 | −98 (0.28%) | 134 (0.8%) | 12 |
| human testis | .060 → .051 | .084 → .076 | −17 (0.15%) | 28 (0.4%) | 2 |
| gorilla OR6737 | .084 → .047 | .138 → .104 | −26 (0.11%) | 60 (0.5%) | 25 |
| gorilla KB3781 | .058 → .034 | .123 → .100 | −25 (0.10%) | 41 (0.3%) | 13 |
| chimp PTR | .103 → .068 | .164 → .131 | −22 (0.10%) | 74 (0.6%) | 10 |
| orangutan PPY | .089 → .052 | .154 → .119 | −37 (0.18%) | 189 (1.4%) | 18 |
| SIRV (dedup) | .164 → .038 | .164 → .038 | 0 | 0; all 48 truth `=` kept | 0 |

- **How much `c` falls** (computed here from the rounded shares): 50% on A119b, 34-44% on the other four non-testis
  samples, 15% on testis and 77% on SIRV. The Outcome and r1138 first said "halves `c` on every substrate", which holds
  on A119b and SIRV only; both were corrected on 2026-09-27 to "cuts `c` by 15-50%".
- **No matched null.** The held-out rows have none, so they do not revisit the dev C5 failure.
- **Unmeasured reach** [CT O.5]. `drop` changes the families-stage representative of 15-19 loci per dev contig, and
  shrinks 8-11 gene spans per contig by ≤ 10 bp; the `flag` stage reads those spans.
- **`tag`** only adds `subchain_of` and `subchain_missing "5p"|"3p"|"both"`. It never says "incomplete": on dev, 63% of
  tagged transcripts are `c`, 30% `j` and 8 of 864 `=` [r1130].

**`--polish-tes pas-end`: judged, KEEP OPT-IN** [HO Outcome; HO-rows; r1136].
- "Moved-set genes" of an arm are the genes whose annotated TES (same last intron, ≤ 50 bp) is hit by that arm's 3′
  end of some transcript that pas-end moved (for BASE, the unmoved end) [HO §6].
- The NULL moves the same transcripts to their most-3′ unprimed own cluster **without** a PAS, so it isolates the PAS.
  Where no such cluster exists it keeps pas-end's move (a tie); \|F\| counts the transcripts where it differs [HO §3].

| sample | moved \|M\| (% multi-exon) | \|F\| | moved-set genes BASE / pas-end / NULL | J1 | J2 | moves toward / away | whole-sample TES50 genes pas-end / NULL |
|---|---|---|---|---|---|---|---|
| human A119b | 5,734 (2.9%) | 1,746 | 249 / 883 / 711 | pass | pass | 1,690 / 2,288 | +128 / +87 |
| human testis | 65 (0.3%) | 13 | 5 / 13 / 13 | pass | **fail (tie)** | 11 / 36 | +1 / +2 |
| gorilla OR6737 | 2,191 (3.2%) | 658 | 114 / 674 / 588 | pass | pass | 1,393 / 480 | +155 / +128 |
| gorilla KB3781 | 2,874 (4.0%) | 1,098 | 149 / 867 / 692 | pass | pass | 1,802 / 586 | +205 / +158 |
| chimp PTR | 1,219 (2.3%) | 313 | 125 / 478 / 419 | pass | pass | 728 / 317 | +101 / +79 |
| orangutan PPY | 2,025 (1.9%) | 463 | 159 / 507 / 432 | pass | pass | 781 / 640 | +67 / +48 |

- **Clauses.** J1 = beats BASE; J2 = beats the NULL. J3/J4 (identical `=` pairs and precision) pass everywhere. SIRV is
  inert: no SIRV end is internally primed.
- **Why it stays opt-in.** EFFECTIVE needed J2 on both human samples, and testis ties 13 vs 13 on 65 moved transcripts.
- **Toward / away** is relative to the nearest annotated TES of the same last intron (read here from HO-rows; the
  Outcome does not discuss it). On both human samples more moves go away than toward, and so do the NULL's (A119b
  1,573 / 2,404; testis 11 / 36), yet the moved-set gene count rises. This is stated, not explained.

**`--polish-tss rescue`** (descriptive) [HO Outcome; HO-rows; r1139]. By design it acts only on contigs with a cap
signal (§5).
- A119b: +1,711 transcripts, +94 chains, +13 TSS250 genes; chain precision .178 → .177.
- Chimp: +245 transcripts, +19 chains.
- Both gorillas, testis and orangutan are unchanged (HO-rows; the Outcome and r1139 omit orangutan). The Outcome does
  not separate "no cap signal" from "nothing proven in scope".

## 5. The both-forms requirement, and transcript ends

**The requirement** [mem: read_proven_ends]:
- emit both the short form and its 5′-extended form when the reads prove both;
- do the same for 3′ ends;
- keep NPIP fusions from distorting the graphs.

**What the assembler did before** (dev, chr20 / gorilla / chr16) [tss_audit; tss_measure].
- **Pass 1 emits one transcript per intron chain.** It keys reads on the exact (contig, intron chain), and each end is
  the single outermost read (k = 1).
  - So a chain with two start clusters (case B) is never emitted twice: 0 of 11,468 / 6,630 / 20,910 raw chains are.
  - The TSS rests on one read within 10 bp for 51-55% of transcripts with ≥ 5 reads.
- **For a short chain nested in a 5′-extended one (case A), the ISM ratio (0.7) is the main loss.**
  - It removes 3,669 / 1,862 / 7,094 of 4,498 / 2,145 / 8,813 raw short forms.
  - The sub-chain drop then removes 227 / 116 / 474 more, only 2 / 2 / 6 of them at an annotated TSS.

**TSS proof (`--polish-tss`, reads only)** [AP addendum 2; tss_critique; r1131].
- **The scan.** A negative-binomial scan tests 5′ ends against a constant truncation hazard, fitted per contig on
  internal exons.
- **The acceptor stratum** is needed because truncated molecules pile at acceptor −2, at 97.6 / 83.1 / 114.5 × the
  hazard (chr20 / chr16 / gorilla). Without it, 66 / 63 / 33% of the first design's rescues were such piles.
- **The cap signature** is a 1-3 bp soft clip of only G at the read's 5′ end: the template-switch G opposite the m7G
  cap.
- **The cap signal.** A contig has one when capped reads are a majority at proven first-exon clusters **and** a
  minority in internal-exon bodies (exact binomial tests against ½).
  - chr20: 70% vs 4%. chr16: 69% vs 5%.
  - Gorilla: 4% vs 1% (capped share of all reads .034), so **no signal**.
  - Without a signal, `rescue` and `split` output exactly `tag`. Before this gate, gorilla rescued 38 forms for +1
    chain against a random-rescue null maximum of 1.
- **G1 fix** [main session 2026-09-27; AP]. A proven window counts for a form only when the form has its own start in
  it.
  - chr20 rescues go from 59 to 47.
  - `off` stays byte-identical (checked on chr20 [AP]), and `cargo test` passes 933.
  - Binary: `tss3_bin_frozen/copy_assign` b1709a96.
- **`split` (case B).** A chain with ≥ 2 cap-proven clusters becomes `<tid>` plus `<tid>_tss<i>`. On chr20 it splits
  49 chains (+54 records).
  - Only 4 have pieces at two distinct annotated TSSs. RefSeq rarely lists same-chain TSS variants, so this is a lower
    bound.
  - Its precision "gain" is a metric trap: an `=` twin counts twice.
  - It changes O1: chr20 goes from 26 to 27 families.
- **Dev rescue, chr20, after G1** [HO §2]. +47 transcripts and +3 chains. A depth-matched random rescue of the pre-fix
  59 forms gives a mean of 0.4 and a maximum of 2 [AP].
  - The rescued forms match a reference chain at .05 / .06 (chr20 / chr16), below the output's precision.
  - **It is a recall-for-precision trade.**

**TES proof (`--polish-tes`, reads plus genome)** [AP addendum 2; r1132; r1135]. ⚠ The triples in the next two bullets
are in AP's order, **chr20 / chr16 / gorilla**.
- **Why the genome.** The reads carry no poly(A) tail: they are FLNC-trimmed, and 78% have a 3′ clip of 0.
- **Why no scan.** The 3′ background is too clumped: a = 114 / 254 / 368, against 1.2 / 1.4 / 6.6 at the 5′ end.
- **The rule.** A cluster of the transcript's own 3′ ends (single linkage, ≥ 2 reads) is PAS-proven when AATAAA or
  ATTAAA lies wholly in [mode−35, mode−10] **and** the mode is not internally primed. Primed means ≥ 60% A in the
  20 bp downstream, or an A6 run (r1064's instrument).
- **Label-free separation.** The rule holds at:
  - .585 / .543 / .650 of emitted 3′ ends;
  - .022 / .074 / .125 of internal-exon 3′ piles;
  - .007 / .007 / .004 of random positions.

  That is a likelihood ratio of 27 / 7.3 / 5.2 (chr20 / chr16 / gorilla): the rule works in both species, but on
  gorilla about 5× less well than on chr20 (1.4× less than on chr16; computed here).
- **The boundary defect** (r1135). The k = 1 end is internally primed on 20.3 / 11.7 / 20.1% of multi-exon transcripts
  (chr20 / gorilla / chr16). It lies > 21 bp beyond the most distal PAS-proven own cluster on 15.6 / 16.0 / 17.5% of
  the chains that have one.
- **`pas-end`.** It moves a primed end upstream to the most-3′ PAS-proven own cluster. No transcript is added or
  removed.
  - Dev (chr20 / gorilla): 106 / 118 moves, a median of 868 / 728 bp. The annotated TES within 50 bp goes 6 → 32 and
    8 → 66. `=` is unchanged.
  - The PAS is what selects: an "any unprimed cluster" target gives a new-near rate of .317 / .574, against .438 /
    .673. The unprimed target recovers slightly more in absolute count (37 vs 32; 77 vs 65) but also loses more (9 vs
    6; 13 vs 7).
  - It changes O1 (chr20 26 → 27 families; gorilla 93 → 95 copies), and no truth scores that.
  - It ignores read support (G2): 7 of 106 chr20 moves go to a smaller cluster.

**Constants and where they come from** (the advisor asked for no arbitrary thresholds) [CT §0, §1; AP]:

| constant | value | status |
|---|---|---|
| sub-chain end tolerance TOL | 10 bp | plateau pick on dev (flat over 5-50 bp; 0 bp halves the effect); **not derived** from an annotation-free measure (critique M1, a stated limitation) |
| sub-chain support ρ | reads(x) ≥ ½·reads(y) | plateau pick (flat over 0.25-0.5); a support guard, not a 1/k division of reads |
| TSS window, level | W = 2·TOL+1 = 21 bp, α = 0.05 | α conventional; hazard, dispersion and acceptor excess r(d) are fitted per contig |
| TES hexamer window, priming | [−35, −10]; ≥ 60% A in 20 bp or A6 | standard PAS position; r1064's priming instrument, unchanged |
| TES cluster | ≥ 2 reads, linkage gap 21 bp | the pass-1 floor and the TSS window |
| label tolerances | TES ± 50 bp, TSS ± 250 bp | **labels only**, never in rule code; 250 bp is the §6w3 5′ dispersion |

## 6. What was decided not to build, and why

**TES rescue of PAS-proven 3′-shorter forms** [r1132].
- It rescued 187 / 19 / 411 forms (chr20 / gorilla / chr16) for +5 / +1 / +11 chains. The matched null gives a mean of
  +1.9 / +1.2 / +4.8 and a maximum of +5 / +2 / +8.
- Precision falls −0.65 / −0.17 / −0.70 pt, about the null's own cost (−0.71 / −0.17 / −0.76).
- At 2.7 / 5.3 / 2.7 chains per 100 added transcripts, this is r1064's regime (≈ 3.3). It is selective on chr16 only.

**Tandem APA split** [r1132].
- 395 / 271 / 490 chains carry ≥ 2 PAS-proven clusters, a median of 528 / 672 / 640 bp apart. RefSeq and Gnomon confirm
  1 / 1 / 1.
- The labels are blind to tandem APA, so the split is **unvalidated, not refuted**.
- SIRV cannot score it either: a canonical PAS sits at 1 of 69 truth TESs.

**Fusion as a relation between two parent loci** [r1134; tesfus_measure §3; tesfus_design §2].
- **Most fusions never become a node.** Of 42 / 31 / 61 events (chr20 / gorilla / chr16) with a PAS-proven parent-A
  end, 23 / 18 / 44 are not assembled or are dropped by the polish.
- **The chr16 family gain (Compara F .615 → .633) is one event, NPIPB2 → GSPT1.** It is a **ghost link**: the polish
  drops the fusion transcript, but `gene_id`, computed before the polish, still joins the two loci.
- **The gain is mostly not its own.** A cut-matched null matches it on Compara F and protein F in 2 of 5 seeds (it
  stays above every seed on Soto, NPIP-U2 and true pairs). Its structural arm is also a subset of the parked regroup
  RG3 (r1127: Compara true pairs 66 → 87, protein referee 100 → 125).
- **PKD1P6-NPIPP1 fails the parent proof.** NPIPP1 has neither end proven: its standalone 3′ ends sit on an Alu A-tract
  (internal priming, no hexamer), and its 5′ ends are uncapped (2%, 0%). PKD1P6 has both (84% capped; AATAAA,
  unprimed).
- **What fragments PKD1P6** is its non-canonical terminal junctions (CG-AG ×2 and CT-AG, on 119 of 130 primary reads)
  under strict canonicity (r1074), which is a closed lever. The fusion does not disturb NPIPP1's family in the dev
  assembly.

**Also not done.** The sub-chain guards are closed (§3). `--polish-tss split` and `tag` were not tested held-out
[HO §10].

## 7. SIRV as a truth (first use for our output)

**What it is.** Only human_testis is spiked: 11,103 SIRV E0 reads plus 17,488 ERCC reads, 2.3% of 1,233,001
[mem: spikein_standards].
- All spike-in reads are unmapped in the genome BAM [mem: spikein_standards]. The SIRV reads were re-aligned to
  SIRV1-7 only, so they are disjoint from the testis genome reads [HO §5].
- The truth has 69 isoforms: 61 multi-exon and 9 true ISMs.
- 93.3% of SIRV primaries are coordinate duplicates, so a `kcd` count mode is scored beside dedup [HO §0; CT §6.4].

**What it can test.**
- **Exact recovery of known complete isoforms.** We recover 48 of 61. The `=` ceiling is 60 (SIRV701 and SIRV705
  share a chain), and SIRV107 is expected to be unreachable.
- **Fragment emission beside known complete forms.** `c` is .164, and .038 under drop.
- **The safety of drop** (descriptive here; no clause judged it). SIRV303 is the only true ISM that is structurally
  exposed, and nothing was lost.

**What it cannot test** [tss_design; tesfus_design; r1132; mem: spikein_standards; CT §8].
- PAS-gated rules, which are inert on it by construction.
- Same-chain TSS or TES variants: it has no case-B truth.
- Paralogs: it has none, so it is not an O1 or O2 truth.
- More than one library. It carries about one isoform of safety power.

Any later prereg that uses SIRV is its second use.

## 8. Caches and file formats [main session 2026-09-27]

- **The cache.** `run_cache.rs` stores TSV, FASTA and PAF, keyed so that stale hits cannot happen. Exporting a
  `RUSTLE_POLISH_*` variable costs a needless catalog-cache miss, never a stale hit [AP].
- **Dependencies.** serde was removed from `Cargo.toml` on 2026-09-24 (`serde_json` stays as a test-only
  dev-dependency; `Cargo.toml:22-24, 92`).
- **Where a binary format would pay.** The only clear win is the ~1.4 GB genome-wide best-AS table `molecules.tsv`,
  which every seeded assembly re-reads. Parquet could help the Python analysis; HDF5 is not worth it.
- **Next step.** Measure the load time first.

## 9. Open items and the user's calls

1. **Should `--polish-subchain drop` become the default?**
   - For: held-out, it cuts `c` by 34-50% on five samples (15% on testis) at ≤ 0.3% of chains, and it keeps all 48
     SIRV isoforms.
   - Against: on dev it costs more intron-correct gene labels than a matched random drop, and the held-out rows have no
     such null.
   - A flip needs its own held-out prereg [CT O.5]. Its reach into the family representatives is unmeasured.
2. **Should `--polish-tes pas-end` become the default?** The verdict is KEEP OPT-IN, on a 65-transcript testis tie. The
   prereg leaves a flip to the user whatever the outcome [HO header]. A flip changes O1 families, with no truth to
   score them.
3. **Regroup.** RG3 (r1127) and fusion `relation` are parked together; one held-out families run costs about 18-24 h
   and 55-60 GB. Ghost links (16 / 4 / 28 dev loci, chr20 / gorilla / chr16) suggest recomputing `gene_id` after the
   polish, which is unmeasured [r1134].
4. **A TES/APA truth.** Only a poly(A)-site atlas (PolyASite or PolyA_DB lifted to CHM13; not on disk) could score case
   B or a TES rescue [r1132].
5. **A fresh substrate.** Any new guard or end rule needs a dev substrate outside V1-V6 [r1129; AP "Why both stay
   opt-in"].
   - The six samples now carry the readthrough R/RQ1 and R3 verdicts and this pas-end verdict [HO Amendment 1].
   - SIRV has been used once.
6. **Small items.**
   - pas-end ignores read support (G2).
   - The 21-bp linkage boundary is untested (G3).
   - A contig missing from `--fasta` fails silently (G4).
   - [AP]'s stale driver remark.
   - Check the cap signal on any new library before claiming anything for `--polish-tss`.

## 10. Provenance and source keys

**Commits.**
- `a03b27ab`: code (its commit message: off byte-identical 818/818 + 150/150, `cargo test` 933 passed; per [AP] the
  818/818 + 150/150 matrix was run on the tss2 build, and the G1 fix's `off` check on chr20).
- `2ac18d2b`: the complete-transcripts prereg and register rows 1118-1135.
- `a3cd54b4`: the held-out Outcome and rows 1136-1139.

**Preregs.**
- [CT] = `docs/PREREG_complete_transcripts_2026-09-27.md`: text sha1 76791696 (948 lines); withdrawn on dev.
- [HO] = `docs/PREREG_completeness_heldout_2026-09-27.md`: frozen text a046be8c (254 lines;
  `complete_ho/frozen/PREREG.sha1`).

**Frozen scripts.** `compat_collapse.py` e14c5646, `complete_null.py` 9740f092, `ho_eval.py` 597fcd08, `tes_null.py`
bdcb8478 and `run_ho.sh` baca765f. The independent checker is `ho_check.py` 690e35cf, with 0 mismatches.

**`copy_assign` binaries.**
- f480847a: dev BASE (`rt3_bin_frozen`).
- 753a3b4d: the 09-25 held-out BASE.
- 6bbd6442: sub-chain (`ct_bin_frozen`).
- fbe4fa86: TSS (`tss_bin_frozen`).
- aebfcf96: TSS fixes and TES (`tss2_bin_frozen`).
- **b1709a96: the G1 fix, used for the held-out runs (`tss3_bin_frozen`).**

**Register** (`docs/NEGATIVE_RESULTS_REGISTER.md`).
- r1062-r1085: prior polish work. r1074-r1076 adopted; r1077 is the read ceiling; r1078-r1083 refuted.
- r1127: RG3.
- r1128-r1130: sub-chain.
- r1131-r1135: ends and fusions.
- r1136-r1139: held-out.

**Other keys.**
- [HO-rows] = `/mnt/linuxdisk/tmp/rustle_figures/complete_ho/check/outcome_rows.md`.
- [AP] = `bench/ASSEMBLY_POLISH.md`, 2026-09-27 addenda 1-2.
- [mem: …] = `~/.claude/projects/-mnt-c-Users-jfris-Desktop/memory/project_*.md`.
- Session reports ([ct_attr], [ct_levers], [ct3_*], [tss_*], [tss2_*], [tesfus_*]) are in
  `/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/931c208e-8acb-4dd2-aacb-cf92d5ad051f/scratchpad/figs/`. That directory is
  **not durable**; the preregs, [AP] and the register hold the durable record.
