# Pre-registration: the Pareto set against a null arm, on the near-tie band alone, fresh one-genome simulations (2026-10-20)

**User (2026-10-20):** after register 1348's audit (a PSV-free null passes the same bars): "ok, lets try both, in order" —
first this test, then the FCGR2A lift check.

## Why

Register 1348 showed that the bars of `PREREG_o2_window_set` certify the 30-point window, not the set: a rule using no read
bases (minimap2's primary when no rival is within 10 points, every copy otherwise) passes them too. The set can show value
only where that null gives nothing: the **near-tie band**, reads whose closest rival is within 10 intron-free points. There
the null keeps every copy, so the question is whether the PSV columns narrow the pool, and at what cost in dropped truths.

## The rule (unchanged from `PREREG_o2_window_set`, code b4594ea2; `bench/o3_graph/o2window.py`)

- Window: alignments within 30 points of the best `ms` on two or more copies of one family, none off every copy, the
  primary among them; stitched alignments excluded.
- **Band under test:** reads first tied within 10 points (`window_copies` at W = 10, i.e. `ms > best − 10`).
- The set: `pareto_set` over the substitution columns of the common window.
- Single pick (fewest edits): descriptive.

## Arms, all scored on the band's reads

1. **The set.**
2. **Null, keep-all:** every copy in the 30-point window. It holds the truth whenever the truth is in the window and
   never narrows.
3. **Null, primary:** minimap2's primary alone (a one-copy set on every read).

## Substrate: round 26, fresh one-genome simulations

- `families26.json` (set `r26`, sha1 1826986addb5, 90 copies): 30 families drawn with `random.seed(2026)` from the 130 qualifying non-PAR families of the
  pool (`r19seq/families19cal_full2_all.json`) that no earlier round (families2-25) used. 2-9 copies each.
- Worlds INTACT, CROSS and SILENT under `profiles/{A119b,public}_testis_v4_ref.json` (no donor alleles), seeds 1601-1603:
  6 cells, 18 worlds, built by `fastset.py run` with round 25's `$V17` (`fast/heldout25/run_sim.sh`).
- Truth: each read's simulated copy. Every copy is in the reference.
- Scored read: determined, truth among the window's copies, and in the band. Pooled over the 18 worlds, reads keyed by
  world.

## Bars (the band only)

- **P:** at least 300 scored reads from at least 11 families.
- **H:** the set holds the truth on at least 99% of scored reads.
- **N (beats keep-all):** the set is one copy for at least 25% of scored reads. (Gorilla development: 42%; S1's band,
  seen after the fact: 38.6%. The keep-all null scores 0.)
- **Q (beats the primary):** of the set's singletons, at least 98% are right, and more are right than minimap2's primary
  on those same reads.

**Verdict:** PASS = P, H, N and Q hold; FAIL = P holds and another bar fails; UNDERPOWERED = P fails.

**A PASS would claim:** in one genome, among reads whose rival copy is within 10 intron-free points, the PSV columns
narrow at least a quarter of them to a single copy that is right at least 98% of the time and more often than minimap2's
coin toss, while keeping the true copy on at least 99% of all such reads. It would not claim anything about donors whose
copies differ from the reference (register 1348's S2), copies absent from the reference, or exact ties.

## Reported alongside, no bars

- the three arms on the whole 30-point window and on the 10-20 and 20-30 bands;
- the single pick's fixes and breaks;
- mean pool and set sizes;
- per family: scored, held, singletons, singletons right.

## Rules

- Round 26 is built and scored once with the code at this file's commit; bugs go into an amendment with the rerun.
- Heavy jobs foreground, one at a time, under 10 minutes each (one `fastset.py run` per command).
- The metric-auditor checks the verdict on a frozen copy before it is reported.
- Human numbers only; nothing pooled with gorilla.

## Amendments

(none yet)
