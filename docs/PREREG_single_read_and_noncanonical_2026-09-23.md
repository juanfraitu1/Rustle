# Pre-registration — can single-read chains and non-canonical junctions be admitted by a rule that holds precision?

**Written 2026-09-23 (§6z9), before any rate is computed.** User: *"can we try to learn how to emit the single
read chains and acceptable non-canonical junctions?"* — the two classes r1062/r1065 found the tools emit and
we refuse by construction (pass-1 floor of 2 reads; canonical-motif gate with the §6m8 majority tolerance).

## What is already known (not predictions)

- r1062: ≥ 95% of missed annotated transcripts have ≤ 1 exact-chain read; U1 − U2 = 707 reference chains on
  chr20-22 carried by exactly one read. isoseq's sensitivity lead at ≥ 1 read (80.6% vs our 68.4%) is this class.
- §6m8 / r866: relaxing canonicity itself recovered 0/40 NPIP junctions, and 36 of those 40 are non-canonical
  — "very likely alignment artifacts". §6au: of 18 deep non-canonical junctions, 4 were jitter off a canonical
  site and 3 recurred at exact coordinates in an independent library (real). The majority mode already
  tolerates a minority of non-canonical junctions ≤ 10 kb when a canonical majority fixes the strand.
- r1065: consensus chains we lack include 388 (human) / 239 (gorilla) with non-canonical junctions shared by
  all three tools — shared alignments, not independent evidence.
- `--rescue-singletons` (B2) already emits low-support chains, but only for certificate-ASSIGNED reads under
  `--gtf`; it is O2-scoped and not a general rule.

## Substrates

- **Development: human A119b chr20/21/22** (primary reads; RefSeq CHM13 as the label source). An independent
  human library (`human_testis`, chr20 slice in `bakeoff/human_chr20`) gives cross-library recurrence where it exists.
- **Held out: gorilla `GGO_mm.bam` NC_073244.2**, RefSeq GCF_029281585.2. No re-tuning.

## Part A — single-read spliced chains (all junctions canonical)

Population: every pass-1 group with exactly one primary read (chain not emitted by the polished output).
**Label:** the chain is exactly a RefSeq intron chain. **Base rate** = label rate over the whole population;
**reference rate** = the polished output's own annotation-exact rate on the same contigs (human 14.5%,
gorilla 32.9%).

Candidate rules, fixed now (each a predicate; reported alone and in the listed conjunctions):

| rule | predicate |
|---|---|
| A1 known-junctions | every junction of the chain is carried by ≥ 2 OTHER reads at the locus |
| A2 novel-isoform | the chain is neither a contiguous sub-chain nor a super-chain of an emitted transcript |
| A3 read quality | MAPQ 60, `de ≤ 0.01`, soft clips ≤ 20 bp at both ends |
| A4 not internally primed | < 60% A (T on −) in the 20 bp downstream of the 3′ end and no 6-A run |
| A5 anchored | the read's 5′ and 3′ ends lie within 50 bp of an emitted transcript's ends at the locus |
| A1∧A3, A1∧A2∧A3, A1∧A2∧A3∧A4, A1∧A3∧A5 | conjunctions |

Also a depth-3 decision tree on the same features, 5-fold CV on human, reported for information only; the
rule that goes to gorilla is the best **named predicate**, not the tree.

## Part B — non-canonical junctions

Population: every junction carried by ≥ 2 primary reads whose motif is not GT-AG / GC-AG / AT-AC on either
strand. **Label:** the junction is exactly a RefSeq intron. Candidate rules:

| rule | predicate |
|---|---|
| B1 no canonical neighbour | no canonical junction within ± 10 bp of either end carries more reads (not jitter) |
| B2 support | ≥ 5 reads |
| B3 majority at site | it is the most-supported junction sharing its donor or its acceptor |
| B4 recurrence | present in the independent library (chr20 only) |
| B5 short | intron ≤ 10 kb (the majority-mode tolerance bound) |
| B1∧B2, B1∧B2∧B3, B1∧B2∧B3∧B5 | conjunctions |

## Selection and bar — committed now

For each part, the rule taken to gorilla is the named predicate with the highest **precision** (label rate)
among those admitting ≥ 50 chains/junctions on human. On gorilla, the rule's admitted set is scored the same
way, and the rule is:

| held-out outcome | verdict |
|---|---|
| admitted-set precision ≥ the polished output's own annotation-exact rate on that contig, and ≥ 100 chains admitted | ⭐ **ADMISSIBLE** — pre-register an implementation arm (gffcompare 15-cell rerun) |
| precision ≥ half that rate | ⚠ **A DIAL** — recall at a stated precision cost, off by default |
| below half | ⛔ **NO** |

**Predicted, before looking:** Part A ⛔/⚠ — single-read chains are dominated by degraded/mis-spliced molecules;
A1∧A2∧A3 should enrich several-fold over the base rate but stay well under the polished rate, because the
polished rate is itself set by 2+-read evidence. Part B ⛔ — the deep non-canonical junctions are mostly
alignment jitter (§6au: 4 of 7 sites), and the annotation labels very few non-canonical introns, so even a
correct rule cannot show a high label rate; B1∧B2∧B3 will be the best and still below half.

I will not change the rules, the label, the substrates or the bar after seeing any number.

---

# OUTCOME (2026-09-23) — ⛔ **NO on both parts, by the bar, on both substrates.**

`/mnt/linuxdisk/tmp/gw22/ext/rules_spike.py` (pysam), primary reads only.

## Part A — single-read all-canonical chains

68% of all spliced chains in the BAM are single-read (human 94,415 of 138,099; gorilla 16,001 of 23,097).

| substrate | population | RefSeq-exact | **base rate** | polished output's own rate | best single rule | best conjunction | depth-3 tree (CV) |
|---|---|---|---|---|---|---|---|
| human chr20-22 (dev) | 58,492 | 704 | **1.20%** | 16.0% | A5 anchored 2.35% (n 12,695) | **A1∧A3∧A5 2.73%** (n 6,811) | 1.72% (admits 23,013) |
| gorilla NC_073244.2 (held out) | 12,251 | 525 | **4.29%** | 36.6% | A5 anchored 5.02% (n 4,087) | **A1∧A3∧A5 9.57%** (n 1,296) | 5.63% (admits 6,550) |

The same predicate is best on both substrates — a chain whose junctions are all seen in other reads, on a
clean read, whose ends sit within 50 bp of an emitted transcript's ends — and it enriches 2.3× over the base
rate. **On held-out gorilla it reaches 9.6% against a bar of 36.6% (⭐) or 18.3% (⚠): ⛔.** No predicate
comes within a factor of 4 of the polished output's precision. Read quality (A3), novelty (A2) and the
internal-priming test (A4) barely move the rate; the ends (A5) are the only feature with signal, which is the
§6w3 story again (5′ dispersion, 3′ tightness), not a rule.

## Part B — non-canonical junctions with ≥ 2 reads

| substrate | population | RefSeq introns among them | base rate | RefSeq non-canonical introns on the contigs | best rule |
|---|---|---|---|---|---|
| human chr20-22 | 11,349 | 12 | 0.11% | 121 of 19,776 | B4 recurrence in the independent testis library: **15 junctions recur, 3 are RefSeq (20%)** |
| gorilla NC_073244.2 | 435 | 0 | 0.00% | 133 of 15,135 | none (every rule 0%) |

Jitter (B1), support (B2), site majority (B3) and length (B5) all stay ≤ 0.4%; the dominant motifs among the
"clean" survivors (GTCT, GTAA, AGAC, TTAC) are the §6m7 artifact classes. The only signal is cross-library
recurrence (§6au's instrument), and it admits 15 junctions on a contig with 4,000 candidates. **⛔ on both
substrates.** The annotation itself carries 0.6-0.9% non-canonical introns, so even a perfect rule has almost
nothing to recover here.

## What this settles

isoseq collapse's sensitivity lead at ≥ 1 read (80.6% vs our 68.4%, r1065) is bought at a per-chain precision
of 1-4% on the single-read class; no structural, quality, priming or end-anchoring predicate lifts that class to
within a factor of 4 of the polished output. Single-read chains stay out of the default; if ever wanted, the
dial is the named predicate A1∧A3∧A5 at a stated ~3% (human) / ~10% (gorilla) precision. Non-canonical
junctions stay gated by the canonical majority rule; cross-library recurrence is the only admissible
evidence, and it needs a second library of the same tissue.

---

# ADDENDUM 1 (2026-09-23) — "how do the other tools do it?" — read from their outputs

Each tool's multi-exon chains, stratified by how many primary reads carry the exact chain (BAM tally) and by
motif, with the RefSeq-exact rate of each stratum (`tool_strata.py`):

| substrate | tool | chains | 0 exact reads | 1 read | 2 reads | ≥ 3 reads | ≥ 1 non-canonical junction | overall |
|---|---|---|---|---|---|---|---|---|
| human chr20-22 | **ours** | 14,370 | 0 | 0 | 4,115 (5.4%) | 10,255 (**20.2%**) | 1,391 (1.2%) | **16.0%** |
| | StringTie 3.0.1 | 12,627 | 1,791 (2.4%) | 2,321 (4.8%) | 3,058 (5.2%) | 5,457 (28.4%) | 2,103 (0.3%) | 14.8% |
| | FLAIR 3.0.1 | 35,791 | 0 | 20,124 (1.2%) | 6,245 (3.0%) | 9,422 (19.5%) | 11,504 (0.2%) | 6.3% |
| | isoseq collapse | 103,283 | 13,629 (0.0%) | 59,780 (0.8%) | 13,238 (2.4%) | 16,636 (11.5%) | 39,598 (0.0%) | 2.6% |
| gorilla NC_073244.2 | **ours** | 4,269 | 0 | 0 | 1,226 (17.0%) | 3,043 (**44.5%**) | 217 (0.0%) | **36.6%** |
| | StringTie 3.0.1 | 3,714 | 219 (6.8%) | 795 (14.3%) | 642 (22.1%) | 2,058 (53.2%) | 167 (0.0%) | 36.8% |
| | FLAIR 3.0.1 | 5,889 | 0 | 2,521 (4.0%) | 755 (12.1%) | 2,613 (45.7%) | 1,012 (0.0%) | 23.5% |
| | isoseq collapse | 16,382 | 1,445 (0.0%) | 10,267 (3.7%) | 1,655 (11.1%) | 3,015 (36.0%) | 3,692 (0.0%) | 10.0% |

**They do not do it; they have no floor.** isoseq emits 59,780 / 10,267 single-read chains at 0.8% / 3.7%
precision and 13,629 / 1,445 chains that NO read carries (its 5-bp fuzzy-junction merge manufactures them) at
0.0%; FLAIR 20,124 / 2,521 singletons at 1.2% / 4.0%; StringTie 2,321 / 795 at 4.8% / 14.3% plus zero-read
chains from flow decomposition. None checks motifs: every tool's non-canonical stratum sits at 0.0-0.3%.
Their singleton precision is the base rate Part A measured (1.2% / 4.3%) — the class is what it is; the tools
pay for it. ⚠ Two levers of OUR default the table exposes: the **majority-tolerated non-canonical chains
(§6m8) are 1,391 / 217 = 9.7% / 5.1% of our output at 1.2% / 0.0% precision**, and our **2-read stratum
(4,115 / 1,226) sits at 5.4% / 17.0% against 20.2% / 44.5% for ≥ 3 reads**.

# ADDENDUM 2 (2026-09-23, written before running) — derive the rule from simulated reads?

A simulation cannot contain the noise class that dominates real singletons (degraded, mis-spliced, pre-mRNA
and library-artefact molecules), so it cannot fit a rule for them. What it CAN give, with truth known, is the
part of the singleton and non-canonical population the ALIGNER manufactures from correct molecules:

- **S1**: every chr20 RefSeq multi-exon transcript × 1 read, HiFi-like errors (`bench/sim_reads.py`, err
  0.001, indel 0.0003), ends jittered ± 30 bp; **S2**: same at err 0.01 (a noisy-read arm); **S3**: 3 reads per
  transcript at err 0.001 (junction context for A1/A5). Aligned with the shipped minimap2 settings to chr20.
- Measured: per read, aligned chain vs truth (exact / junction shifted ≤ 10 bp / missing / extra / non-canonical
  introduced); the **motif spectrum of aligner-introduced non-canonical junctions** against the real data's
  (GTCT, GTAA, AGAC, TTAC, CCAG, CTAG …); and, in S1/S2, which read-level predicates (A3, A4) separate exact from
  artefact chains.

**Prediction:** the aligner corrupts < 5% of HiFi-like reads and < 15% of noisy ones; its non-canonical
junctions reproduce the real data's motif spectrum (proving the real ones are alignment artefacts, not
biology); and read-quality predicates separate artefact from exact chains well in the sim — which does NOT
transfer, because the real base rate (1.2%) is ~20× below what the sim's aligner-only error can explain.
The sim therefore bounds the rule: no read-level predicate derived from it can lift real singletons above
the aligner-artefact ceiling, and it will not change the ⛔.

## ADDENDUM 2 — OUTCOME (2026-09-23)

chr20, 4,295 multi-exon RefSeq transcripts, aligned to chr20 with the shipped minimap2 settings
(`/mnt/linuxdisk/tmp/gw22/sim/`):

| arm | reads | aligned chain = truth | missing junction(s) | other | shifted ≤ 10 bp | reads with an aligner-made non-canonical junction | A3 pass → exact / A3 fail → exact |
|---|---|---|---|---|---|---|---|
| S1 err 0.001, 1 read/tx | 4,295 | **94.1%** | 3.6% | 1.9% | 0.4% | 1.6% | 95.5% (n 4,232) / 3.2% (n 63) |
| S2 err 0.01, 1 read/tx | 4,295 | 93.0% | 3.9% | 2.0% | 1.0% | 2.0% | 99.0% (n 103) / 92.9% (n 4,191) |
| S3 err 0.001, 3 reads/tx | 12,885 | 94.3% | 3.4% | 1.9% | 0.3% | 1.5% | 95.7% / 1.6% |

**The aligner-only ceiling on a singleton's precision is ~94%; the real rate is 1.2% (human) / 4.3%
(gorilla).** So ≥ 95% of real single-read chains are wrong for reasons the aligner did not cause — molecules
that are not the annotated transcript (mis-spliced, degraded, pre-mRNA, unannotated, library artefacts) — and
no read-level rule learned on a simulation can see them. The read-quality predicate A3 separates aligner
artefacts almost perfectly in the sim (95.5% vs 3.2%) and moved the real rate from 1.20% to 1.31%: confirmed
non-transferable. At err 0.01 the same predicate fails 97.6% of reads (its `de ≤ 0.01` term), so it is an
error-rate gate, not a biology gate.

**Non-canonical junctions: the sim does NOT reproduce the real spectrum.** The aligner manufactures them at
1.5-2.0% of reads with GGAG, ATAG, AGGC, GTTC, TTAG, GTAC; the real ≥ 2-read non-canonical junctions are
GTCT, GTAA, AGAC, TTAC, CCAG, CTAG. The prediction ("shared alignment artefacts") was wrong on the motif
spectrum and right on the rate order: the real load (38% of isoseq's chains carry one) is far above what
correct molecules produce through this aligner, so those junctions come from molecules that are not annotated
transcripts, and the annotation calls them wrong at 99.9%. Either way they stay out.

**One aligner artefact that DOES pass our floor:** in S3, 375 distinct wrong chains arise from 730 wrong-chain
reads, and **161 of them are carried by ≥ 2 reads (3.75% of transcripts), 114 of those "missing junction"** — the
aligner reading through a short exon, consistently across reads of the same molecule. That is a precision lever on the ≥ 2-read side (retained-intron chains spanned by a dominant
junction), not a recall rule for singletons; it is noted for a separate pre-registration.

# ADDENDUM 3 (2026-09-23, written before running) — are our floors too demanding? Containment and collinear support

User: *"maybe these tools allow for a 'worse' prediction if for example an entire read supports an exon without
a dip in precision or if there is some colinearity."* Our pass-1 counts only reads carrying the EXACT chain.
isoseq collapse folds a 5′-shorter read into the longer chain that contains it (its FL count pools); StringTie
builds a chain from overlapping fragments. Three predicates on the same Part A population (single-read,
all-canonical chains), same label and bar:

| rule | predicate |
|---|---|
| A6 contained support | ≥ 1 OTHER read whose chain is a contiguous sub-chain of this chain (so pooled support ≥ 2, the isoseq collapse rule) |
| A6b | ≥ 2 such reads (pooled ≥ 3) |
| A7 collinear cover | every junction of the chain is carried by some other read whose chain is a contiguous sub-chain of it, and consecutive such fragments overlap (share a junction) — a StringTie-style path with no full-length witness |
| A6∧A3, A6∧A3∧A5, A7∧A3 | conjunctions with read quality / end anchoring |

Also reported: the same predicates on the tools' own zero- and one-read chains are already in Addendum 1
(StringTie's zero-read chains 2.4% / 6.8%, isoseq's 0.0%), which is what those rules yield when a tool applies them.

**Prediction:** A6 enriches over A1 (a contained read is stronger evidence than a shared junction) but stays
far below the polished rate, because 5′-truncated fragments of a WRONG chain are as common as of a right one;
A7 lands near StringTie's zero-read precision. ⛔ expected on both substrates.

## ADDENDUM 3 — OUTCOME (2026-09-23)

| predicate (single-read, all-canonical chains) | human chr20-22: n / precision | gorilla NC_073244.2: n / precision |
|---|---|---|
| base rate | 58,492 / 1.20% | 12,251 / 4.29% |
| A6 ≥ 1 contained read (isoseq-collapse pooling) | 44,166 / 1.35% | 8,181 / 4.36% |
| A6b ≥ 2 contained reads | 41,722 / 1.31% | 6,694 / 4.18% |
| **A7 collinear cover** (every junction carried by contained fragments that chain together) | 6,011 / 1.90% | 208 / **12.02%** |
| A6 ∧ A3 ∧ A5 | 11,234 / 2.39% | 2,622 / 5.45% |
| A7 ∧ A3 | 5,476 / 2.03% | 170 / 12.35% |
| **A7 ∧ A3 ∧ A5** | 1,459 / **4.18%** | 83 / **15.66%** |
| (previous best) A1 ∧ A3 ∧ A5 | 6,811 / 2.73% | 1,296 / 9.57% |
| polished output's own rate | 16.0% | 36.6% |

**Containment pooling is worthless as evidence**: 75% of human singletons (67% gorilla) have a contained
fragment, and the rate does not move (1.35% / 4.36% vs 1.20% / 4.29%) — a 5′-truncated fragment of a wrong
chain is exactly as common as of a right one, so isoseq's FL-count pooling adds count without adding truth.
**A collinear cover with no full-length witness is the strongest predicate found** — 2.8× base on gorilla —
which is StringTie's flow model, and it lands where StringTie's own zero-read chains land (Addendum 1: 2.4% /
6.8%) or a little above once quality and end-anchoring are added. On held-out gorilla the best conjunction
reaches **15.7% against the ⚠ bar of 18.3% and the ⭐ bar of 36.6%: ⛔ NO**, on 83 chains (13 true).

So the floors are not too demanding in the sense that matters: relaxing "one read carries the chain" to
"fragments chain across it" buys 61 true chains on chr20-22 for 1,398 false ones. If a recall dial is ever
wanted here, **A7 ∧ A3 ∧ A5 is the one to ship** — it is a structural rule (collinear path + clean read +
anchored ends), not a fitted threshold — at a stated ~4% / ~16% precision, off by default.

# ADDENDUM 4 (2026-09-23) — "fold in the 5′ end": let a single-read 5′ extension represent the emitted chain

The isoseq-collapse convention: reads that differ only by 5′ truncation are one isoform, represented by the
longest. Test: single-read canonical chains that contain an emitted transcript's chain with extra junctions at
ONE end only; label both forms against RefSeq (`fold5.py`; the same bar as Part A applies to the added set).

| substrate | extension side | n | singleton (long) is RefSeq | emitted (short) is RefSeq | fold-in FIXES (long true, short false) | fold-in BREAKS if it replaces (short true, long false) |
|---|---|---|---|---|---|---|
| human chr20-22 | **5′** | 8,938 | 199 (2.2%) | 455 (5.1%) | 174 (1.9%) | 430 (4.8%) |
| | 3′ | 4,045 | 63 (1.6%) | 750 (18.5%) | 49 (1.2%) | 736 (18.2%) |
| gorilla NC_073244.2 (held out) | **5′** | 716 | 60 (8.4%) | 115 (16.1%) | 54 (7.5%) | 109 (15.2%) |
| | 3′ | 628 | 24 (3.8%) | 351 (55.9%) | 8 (1.3%) | 335 (53.3%) |

⛔ **The shorter form we already emit is the annotated one 2.5× more often than its single-read 5′ extension**
(human 455 vs 199; gorilla 115 vs 60), so folding the 5′ end in as a replacement breaks 2.5 hits for every
one it fixes, and adding the long form instead admits 8,938 / 716 chains at 2.2% / 8.4% — below half the
polished rate (16.0 / 36.6) on both substrates. On the 3′ side the shorter form is right 12-14× more often
(polyA-anchored ends; a single-read 3′ read-through is almost never the transcript). This is the junction-level
face of §6w3 (5′ ends are dispersion) and of r861/r1064 (the ISM side): the 5′-longer read is not the better
witness. isoseq's convention is wrong on this library at the annotation level, and it is part of why its
≥ 3-read stratum sits at 11.5% (r1069).

# ADDENDUM 5 (2026-09-23, written before running) — floor of ONE read plus a learned noise cut

User: *"can we lower the floor to just 1 long read as long as we apply any other rule that cuts through the
noise?"* Addenda 1-4 tested hand-named predicates. This arm lets a model find the cut on a rich feature set,
with the operating point fixed by rule:

**Population and label** as Part A (single-read all-canonical chains; RefSeq-exact). **Features** (per chain,
from its read and its locus): MAPQ, `de`, identity (1 − NM/aligned length), 5′/3′ soft clips, read length,
aligned length, number of introns, min/median/max intron length, min internal exon length, first/last exon
length, strand agreement between junction motif and `ts` tag, min and median junction support from OTHER
reads, number and fraction of junctions carried by no other read, whether a novel junction's donor or acceptor
is reused by another junction, distance of any novel junction to the nearest other junction (jitter), number
of emitted transcripts at the locus, relation to them (exact / sub / super / other), distance of the read's 5′
and 3′ ends to the nearest emitted transcript end, locus read depth, downstream A-fraction and longest A-run
(priming), contained-fragment count (A6), collinear cover (A7).

**Model:** gradient-boosted trees (sklearn `HistGradientBoostingClassifier`, default depth, 300 iterations,
class-balanced), 5-fold CV on human chr20-22 → out-of-fold scores; then fit on all of human and score gorilla.

**Operating point, fixed now:** the human out-of-fold precision-recall curve is read at the ⭐ bar (16.0%) and
the ⚠ bar (8.0%); the score threshold at each is the smallest score whose admitted set has that precision with
≥ 100 chains. Those two thresholds are carried unchanged to gorilla; the verdict is the gorilla precision of
the admitted set at the ⭐ threshold (bar: ≥ 36.6% ⭐, ≥ 18.3% ⚠, else ⛔), with ≥ 100 chains admitted. Also
reported: gorilla's own out-of-fold curve (information only), and the model's top features, to see whether a
named rule approximates it.

**Prediction:** on human the curve reaches 16% only at a few hundred chains (the ends and the collinear cover
carry it), and the carried threshold reaches ⚠ at best on gorilla. If it reaches ⭐, the recall it buys will be
a few hundred chains per contig-set against ~700 single-read annotated chains — worth having as a dial, not a
change of floor.

## ADDENDUM 5 — OUTCOME (2026-09-23)

`single_feats.py` (33 features) → `HistGradientBoostingClassifier` (300 iterations, lr 0.05, balanced), 5-fold
out-of-fold on human chr20-22 (58,492 single-read chains, 704 annotated), then fit on all human and carried to
gorilla NC_073244.2 (12,251 / 525).

| operating point | human out-of-fold | gorilla, carried threshold (held out) |
|---|---|---|
| ⭐ threshold (score ≥ 0.8007, the smallest reaching 16.0% with ≥ 100 admitted) | 293 admitted, **16.0%**, 47 true (6.7% of the annotated singletons) | 220 admitted, **11.8%**, 26 true |
| ⚠ threshold (score ≥ 0.7009, 8.0%) | 2,087 admitted, 8.0%, 167 true (23.7%) | 1,126 admitted, 10.4%, 117 true |
| top 100 / 500 / 2,000 by score | 16.0% / 12.2% / 8.1% | 12.0% / 11.0% / 8.9% |

⛔ **At the carried ⭐ threshold gorilla sits at 11.8% against the ⚠ bar of 18.3% and the ⭐ bar of 36.6%.**
Information only — gorilla's OWN 5-fold model reaches 36.7% at 158 chains (top-100 46%), i.e. a per-library
self-trained cut can meet the bar for ~150–300 chains, which is 7–11% of that library's annotated single-read
chains; the cross-library model does not transfer the operating point (base rates 1.2% vs 4.3%, different
depth). Top features by permutation importance on gorilla: number of emitted transcripts at the locus, min and
median junction support from other reads, 5′ and 3′ end distance to an emitted transcript, fraction of novel
junctions, read length — the A1/A5 story with continuous resolution.

**Verdict on the floor:** a floor of one read is not rescued by a learned noise cut. The best any cut can do,
even trained on the same library, is a few hundred chains per contig-set at the polished precision, ~10% of the
single-read truths; carried across libraries it lands at the tools' level (≈ 12%). The floor stays at two
reads. The learned score is the natural recall DIAL if one is ever wanted (admit the top-k by score at a stated
precision), but it is a per-library fit, not a rule.
