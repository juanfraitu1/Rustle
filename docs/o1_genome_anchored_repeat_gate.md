# A P1-safe repeat gate: fix the universe, not the statistic

**Status 2026-08-19.** Candidate derived and measured offline. **Not pipeline-confirmed (T8).**
Extends [`o1_coverage_repair.md`](o1_coverage_repair.md) §5 and
[`o1_false_positive_rules.md`](o1_false_positive_rules.md); it does not replace their results.

## 1. The defect being fixed

Candidate **R5** — block promiscuity — was the strongest discriminator ever found for O1's false
merges, and it was rejected because it **breaks P1**: MRPS17's block scores **50** partners over the
whole catalog and **1** from a 4-node seed, so the same pair is rejected run-whole and accepted
run-from-seed. `E_r` would stop being a relation between two sequences and become a function of the
node set.

The shipped repeat-hub gate has the same disease. `bench/vg_repeat_catalog.py` states its universe
outright: *"Universe = the 3462 gene-assigned loci that participate in E_r families."* It counts
**catalog members**. (It is confined to `family_define`/`driver.rs`; `gw_family_catalog` never calls
it, so the shipped O1 definition is not currently exposed.)

**The diagnosis was wrong in one word.** R5's disease is not *"it is a repeat count"* — it is
*"it is counted over the NODE SET."* Change the universe and the statistic becomes pair-local.

## 2. The statistic

Count occurrences in the **fixed reference assembly** instead of in the catalog. `GGO.fasta` is not
a function of the seed, so the value depends only on `(a, b, genome)`.

```text
block_a, block_b       the aligned interval on each rep, from the one passing E_r record
S                      canonical k-mers (k=21) present in BOTH blocks -- their shared anchors
g(x)                   occurrences of anchor x in GGO.fasta, both strands
min_shared_gmult(a,b)  = min over S of g(x)
```

Reading: *does the sequence these two loci share contain **any** anchor that is rare in the genome?*
A real paralogue pair shares private sequence, so some anchor is near-unique. A pair bridged only by
a mobile element shares nothing but high-multiplicity anchors.

Library-free: no RepeatMasker, Dfam or RepBase enters the definition. Strand-agnostic: the query set
carries each anchor and its reverse complement, so scanning the assembly's forward strand yields the
two-strand count. `k = 21` is odd, so no anchor is its own reverse complement and none is double
counted.

**Substrate:** the exon-sum spliced reps in `GGO_gwcat.copies.fa` (verified seqlen/genomic-span
0.076–0.447), never genomic spans.

## 3. It passes the test R5 failed

Both statistics computed on the **same pairs** under two universes — the whole 1,415-rep catalog, and
a seed containing only the pair's own family. R5-analogue = number of distinct reps in the universe
carrying any of the pair's shared anchors.

| statistic | pairs whose value CHANGES between whole-catalog and seed |
|---|---:|
| R5-analogue (counted over the node set) | **94 / 147 = 0.6395** |
| **genome-anchored `min_shared_gmult`** | **0 / 147 = 0.0000** |

Largest R5 swings, with the genome-anchored value beside them:

| pair | arm | R5 whole → seed | g-mult (both universes) |
|---|---|---:|---:|
| TP18481 | TP | 228 → 8 | 5 |
| FP00050 | FP | 189 → 7 | 12,663 |
| FP00051 | FP | 178 → 7 | 12,973 |
| FP00052 | FP | 136 → 7 | 8,952 |
| FP00049 | FP | 114 → 6 | 1,347 |

The FP swings are the GWFAM210 MRPS17 **AluY** hub — the mechanism the anti-FP characterisation
identified as Group 1. Under the genome universe those pairs score 1,347–12,973 regardless of what
else is in the catalog.

## 4. Discrimination on the frozen arms

Unit = **pair**. GGO only. FP arm = the 14 gorilla catalog false merges; TP arm = the 150 true pairs.

| arm | n scored | median `min_shared_gmult` | range |
|---|---:|---:|---|
| FP | 12 | **182** | 1 – 12,973 |
| TP | 135 | **2** | 0 – 44 |

**AUC (FP scores higher than TP) = 0.9429.**

| cut `M` | FP rejected | TP lost | TP cost |
|---:|---:|---:|---:|
| 10 | 11/12 | 7/135 | 0.0519 |
| 20 | 10/12 | 2/135 | 0.0148 |
| **50** | **10/12** | **0/135** | **0.0000** |
| 100 | 9/12 | 0/135 | 0.0000 |
| 500 | 4/12 | 0/135 | 0.0000 |

At `M = 50` the rule rejects **10 of 12 scored false merges at zero cost to 135 true pairs**. For
comparison on the same FP arm: coverage-of-longer at 0.20 rejected 1/14, R2@0.05 rejected 2/74
overall, and the transcript-orientation guard rejects 6/14 GGO catalog FPs.

## 5. Where it abstains, and why that is honest

The statistic requires at least one **exact** shared 21-mer. It has none on **11 of 158 pairs = 7.0%**
— 2 FP and 9 TP — and every one of them sits at identity **0.6927–0.8031**, where exact 21-mers are
expected to be rare (0.78²¹ ≈ 0.005 per position). Abstention is not biased toward one arm.

Those pairs must be reported `GMULT_UNMEASURED` and fall through to the incumbent rule. Absence of a
shared anchor is missing data, never evidence of a repeat bridge.

## 6. Is it just softmask restated?

Largely the same signal, but not identical, and library-free.

| | count |
|---|---:|
| FPs with softmask ≥ 0.70 on BOTH sides | 10/12 |
| FPs rejected by `min_shared_gmult ≥ 50` | 10/12 |
| rejected by g-mult, MISSED by softmask | **1** — FP00048, softmask 0.689 / 0.643, g-mult 76 |
| caught by softmask, missed by g-mult | **1** — FP00055, softmask 0.929 / 0.735, g-mult 19 |

So it is neither a superset nor a restatement of the softmask gate. Its value over softmask is that
it takes no repeat library into the definition.

**FP00058 scores 1 and is correctly NOT rejected.** That case is the LAGE3 processed pseudogene
against its own parent — a **truth-label failure**, not a false merge. A statistic that rejected it
would be wrong.

## 7. What is NOT established

1. **T8.** Offline re-derivation. Nothing here has been through the shipped binary. The E_r rule is
   unchanged and no default moves on this evidence.
2. **The AUC is not held out.** The mechanism of these same 14 FPs was characterised as
   repeat-driven (Group 1 AluY hub, Groups 2–3 low-copy elements), so scoring them on repeat content
   is partly circular. **The TP half — 0/135 at `M = 50` — is the load-bearing number**, because the
   TP arm was not selected for repeat content. A held-out FP set is required before quoting the FP
   rejection rate as a rate.
3. **Coverage is 12/14 FP and 135/150 TP.** Six TP pairs have no single passing record under this
   rule and eleven pairs abstain; both are stated above rather than dropped silently.
4. `M = 50` was read off this table. It must be fixed on a held-out set, or the rule quoted as an
   ordering rather than a threshold.
5. **Gate or flag is undecided.** As a gate it changes `E_r` and every downstream certificate; as a
   flag it costs nothing and P1 is untouched either way.

## 8. Can it REPLACE coverage? No — refuted

The tempting move is to drop the coverage clause entirely. Coverage is scale-free (a ~1 kb Alu is
≥0.50 of any node under 2 kb), which is the whole named hole; a rare-anchor test is structurally
immune to that, because a repeat has no rare anchor at **any** node size. So:

```text
E_r_free(a,b)  iff  exists a record with identity >= tier floor
                    AND min_shared_gmult(record) < M
```

**This cannot be judged on the FP/TP arms** — both were built from pairs that already pass
coverage ≥ 0.50. A rule that drops coverage must be judged on what coverage is currently holding
back. All-vs-all over the 1,415 shipped GGO reps, both tiers, every identity-passing pair split by
whether it clears coverage:

| | pairs | within shipped family | cross-family |
|---|---:|---:|---:|
| `COV_PASS` (the shipped edge set) | 2,727 | 2,474 | 253 |
| `COV_FAIL` (what coverage rejects) | 14,111 | 830 | 13,281 |

(The 253 cross-family `COV_PASS` pairs are not an anomaly: families are γ-quasi-cliques *of* E_r
components, so a cross-family E_r edge is by definition an edge γ cut.)

| `M` | shipped edges kept | NEW edges from `COV_FAIL` | of which cross-family |
|---:|---:|---:|---:|
| 2 | 835/2,727 = 0.306 | 131 | 23 |
| 3 | 1,415/2,727 = 0.519 | 437 | 214 |
| 5 | 1,822/2,727 = 0.668 | 878 | 528 |
| 10 | 2,173/2,727 = 0.797 | 1,463 | 967 |
| 20 | 2,363/2,727 = 0.867 | 1,915 | 1,310 |
| 50 | 2,459/2,727 = 0.902 | 2,239 | 1,567 |

**There is no operating point.** Holding new cross-family edges at parity with γ's existing 253
puts `M` near 3, which discards **48% of the shipped edge set**. Recovering 90% of the shipped
edges costs **1,567 new cross-family edges, 6.2× the 253 that exist**. Recall loss starts before
merge suppression does — the same shape as the coverage-repair impossibility argument.

**Why, mechanistically.** The TP distribution has median `min_shared_gmult` = 2, so an admission
criterion strict enough to exclude repeats (`< 2`, i.e. a genuinely unique anchor) also excludes
most real paralogue pairs. The statistic's discriminative power lives at the **top** of its range —
it separates "definitely a mobile element" (100–13,000) from everything else, and does not separate
within the low end at all.

⟹ **The genome-anchored statistic is a VETO, never an admission criterion.** It belongs on top of
the coverage clause, not in place of it. Coverage stays.

## 9. Reproduction

```bash
cd /mnt/linuxdisk/home/juanfraitu/o1_gmult
python3 blocks.py     # recover each pair's aligned block from the shipped reps
python3 gmult.py      # one streaming pass over GGO.fasta, ~9 min, ~3 GB RSS
python3 eval.py       # seed-invariance, discrimination, softmask additivity
python3 covfree.py    # §8: can it replace coverage? (all-vs-all + one genome pass)
```

Outputs: `pair_blocks.tsv`/`.fa`, `gmult.tsv`, `seed_invariance.tsv`, `covfree.tsv`.
