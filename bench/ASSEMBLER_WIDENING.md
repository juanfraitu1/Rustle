# Assembler read-isoform widening + single-exon strand — vs `docs/PREREG_assembler_widening_2026-09-18.md` (md5 6d586b2d5ec6d7a6b15e014ef7fa5349)

**Verdict: the port is CORRECT and SHIPPED (opt-in), the strand fix WORKS, and neither closes the
reconstruction gap. W-2 FAILED — the widening is worth +2 junctions on the cluster, not the +21 §6m5
projected.** That projection was wrong and is corrected here.

Lib suite **880 passed / 0 failed** (3 new tests).

## What shipped

`denovo_assemble::pass1_skeletons_widened(reads, min_reads, min_terminal_support, snap, isoform_k)` —
a chain is admitted if `n >= min_reads` **OR** every one of its junctions has ≥ `isoform_k` reads, counted
per `(chrom, junction)` over all spliced reads. Chains are never concatenated; the pass can only ADD
skeletons. `copy_assign --read-isoform-k`, **default 0 = off**. `pass1_skeletons_robust_with` delegates at
`isoform_k = 0`, so every existing caller is unchanged.

## Results — 9-copy cluster (106 junctions) and all 26 copies (249)

| arm | transcripts | single-exon | **'+' fraction** | 9-copy | all 26 | complete |
|---|---|---|---|---|---|---|
| k=0 baseline | 750 | 97 | **1.000** | 65/106 (61%) | 151/249 (61%) | 3/26 |
| k=2 | 1201 | 97 | 1.000 | 67/106 (63%) | 162/249 (65%) | 3/26 |
| **k=3** | 1124 | 97 | 1.000 | **67/106 (63%)** | **161/249 (65%)** | 3/26 |
| k=5 | 1052 | 97 | 1.000 | 66/106 (62%) | 156/249 (63%) | 3/26 |
| k=8 | 997 | 97 | 1.000 | 65/106 (61%) | 153/249 (61%) | 3/26 |
| k=3 + `RUSTLE_GATE_MIN_READS=2` | 1124 | 97 | 1.000 | 67/106 (63%) | 161/249 (65%) | 3/26 |
| **`RUSTLE_READ_STRAND=1`** | 750 | 97 | **0.753** | 65/106 (61%) | 151/249 (61%) | 3/26 |
| READ_STRAND + k=3 + gate 2 | 1124 | 97 | **0.753** | 67/106 (63%) | 161/249 (65%) | 3/26 |

Targets: primary-read union **85/106 (80%)**, §6m4 ceiling **91/106 (86%)**.

| rule | bar | result |
|---|---|---|
| **W-1** off ⇒ byte-identical | exact | **PASSED** — `k=0` GTF is byte-identical to the §6m5 baseline |
| **W-2** ≥ 74/106 at some k | +10 | **FAILED** — best 67/106 (+2). All-26 gains +10 (151 → 161) |
| **W-3** ≤ 91/106, report transcript count | — | **PASSED** on coverage, but the trade is poor: **+50% transcripts (750 → 1124) for +2 junctions** |
| **S-1** single-exon `'+'` < 0.90 | 0.90 | **PASSED** — 1.000 → **0.753** |
| **S-2** no coverage loss | — | **PASSED** — identical (65/106, 151/249), exactly as pre-declared |

## Correction to §6m5

§6m5 wrote *"port read-isoform widening into the assembler (worth up to +21 junctions, 60.4% → 80.2%)"*.
**Measured: +2 on the cluster, +10 across 26 copies.** The projection assumed every junction visible in
reads was recoverable by a per-chain admission rule. It is not.

## Where the remaining 18 junctions actually go — still open

Not exact-chain collapse alone, and not the gate read floor:
- Of the 18 annotated junctions still missing at k=3: median depth **19**, median **18** distinct carrier
  chains, and a median of **9** carrier chains in which EVERY junction already has ≥ 3 reads. Only **4 of
  18** have no qualifying chain at all.
- So for 14 of 18 a qualifying skeleton IS admitted at pass 1 and the junction still never reaches the GTF.
- `RUSTLE_GATE_MIN_READS=2` changed **nothing** (identical transcript count and coverage), so the gate's
  read floor is not the filter either.

**The loss is between pass-1 admission and GTF emission, and it is not the read floor.** Remaining
candidates, untested: `assemble_gate`'s `max_span` / `min_spliced` / `max_spliced` / `pool_locus_support`,
the `collapse_loci_groups` gene grouping, or the `--gtf-copy-set` evidence rule that drops transcripts at
copies without evidence. That is the next thing to instrument — not more reads and not a bigger k.

## Change 3 (all-alignment pool) — NOT RUN

Feeding secondary alignments into assembly needs a code change (`reads_in_region` takes primaries), and the
entire remaining headroom to the ceiling is **+6 junctions**. Given W-2 failed at the much larger lever,
this was not spent. It stays available and is recorded as not done, not as done-and-negative.

## Recommendation

Keep `--read-isoform-k` at its **default 0**. It is correct, tested and byte-identical off, but at +2
junctions for +50% transcripts it is not a default. `RUSTLE_READ_STRAND=1` is a different matter — it fixes
a documented placeholder at zero measured cost and is worth proposing as a default on its own evidence
(386/400 vs 0.4867), independently of this experiment.
