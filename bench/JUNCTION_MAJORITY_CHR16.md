# Splitting the `build_spliced_seq` bucket — and the chr16 arm the code asked for

§6m7 localised the whole pass-1 → GTF loss to `build_spliced_seq` returning `None` (343 of 1,093 skeletons,
31.4%, at the default). This splits that bucket the three ways §6m7 named and measures the one recoverable
part.

## The default is STRICT: one bad junction kills the chain

`build_spliced_seq_with` has two modes (`denovo_assemble.rs`):
- **strict (the default)** — `junction_strand` returning `None` for ANY junction rejects the whole
  transcript, as does any strand disagreement. This is failure mode **(c)**: *one bad junction in an
  otherwise canonical chain*, and it is the shipped behaviour.
- **majority (`RUSTLE_JUNCTION_MAJORITY=1`)** — strand by canonical majority; a non-canonical junction is
  tolerated when its intron is ≤ `RUSTLE_JUNCTION_NC_MAX_BP` (10 kb); still rejects a spliced model with NO
  canonical junction, and still rejects a real strand conflict.

That is why §6m7 found `RUSTLE_JUNCTION_NC_MAX_BP` to be a no-op: `nc_max` is only consulted in majority
mode, which is off.

## The three-way split (k=0, 27 NPIP windows)

| bucket | skeletons | share of the 343 |
|---|---|---|
| **(c) one/few non-canonical junctions in an otherwise canonical chain** — recovered by majority | **238** | **69.4%** |
| (a)+(b) no canonical junction anywhere, strand conflict, or fetch failure | 105 | 30.6% |

`rej_seq` 343 → **105**; kept 750 → **988**. At k=3: `rej_seq` 485 → 152, kept 1,124 → **1,457**.

## What it recovers (canonical-only truth, per §6m7)

| arm | tx | canonical 9-copy (94) | canonical all-26 (209) | complete chains |
|---|---|---|---|---|
| k=0 strict (baseline) | 750 | 65 (69%) | 151 (72%) | 8/26 |
| k=3 strict | 1,124 | 67 (71%) | 161 (77%) | 9/26 |
| **k=0 MAJORITY** | 987 | **72 (77%)** | **162 (78%)** | **9/26** |
| **k=3 MAJORITY** | 1,456 | **74 (79%)** | **170 (81%)** | **10/26** |

Majority alone is worth **+7 cluster / +11 family junctions** — more than the widening delivered — and the
two compose: k=3 + majority is **+9 / +19** over the baseline, with complete chains 8 → 10.

## The cost, measured

| | strict | majority |
|---|---|---|
| junctions emitted | 613 | 759 (**+146, −0**) |
| of the 146 added: **read-supported** | — | **140 (95.9%)** |
| of the 146 added: canonical motif | — | 78 (53.4%) |
| median transcript span | 14,480 bp | 14,602 bp |
| transcripts > 100 kb | 74 | 76 |
| max transcript span | 593,043 | 593,043 |

**Nothing is removed and nothing runs away.** 96% of the added junctions are carried by real reads; the
span distribution barely moves, so the engulfment/over-merge failure the flag's comment warns about does
not appear here.

## ⚠ This is NOT the full chr16 arm the comment demands

`build_spliced_seq_with`'s comment says the known harm *"was measured ON chr16"* and that a **chr16 arm is
required before flipping**. This arm is **27 NPIP windows on chr16**, not a genome-wide chr16 run, and the
documented harm was genome-wide-shaped (families 121 → 117, copies 678 → 700, strictly-engulfed 60 → 63).

So this is **necessary but not sufficient**. It establishes that at NPIP the flag is a clear, cheap win
with no measurable span inflation. It does **not** clear the flag for a default flip — that still needs a
genome-wide chr16 catalog arm priced on families/copies/engulfment, which is the natural next experiment
and now has a concrete reason to be run.

**Recommendation:** use `RUSTLE_JUNCTION_MAJORITY=1` for pseudogene-family reconstruction work today; do not
change the default until the genome-wide chr16 arm exists.
