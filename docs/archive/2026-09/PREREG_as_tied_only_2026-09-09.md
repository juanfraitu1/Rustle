# PREREG — O2 admits AS-TIED MULTIMAPPERS ONLY (user correction, 2026-09-09)

**Written before the code change and before any re-run.**

## The correction (user, verbatim intent)
*"Copy assignment should only work for tied-AS multi-mapping reads, nothing else, so we need to ensure
anything else does not enter O2 at all."*

## Why it matters — measured BEFORE the change, on the two arms already in hand
`as_best` / `as_second` in `assignments.tsv` are **minimap2's own AS across the read's placements**
(`as_evidence_per_read`, from `BamRead::as_score`); the code notes they never feed the decision, so they are an
independent description of what O2 was handed.

| arm | assigned | AS exact tie | near-tie (≤5) | clear best (>5) | single placement |
|---|---|---|---|---|---|
| gorilla MCL1 (80 copies) | 6,142 | **6** | 14 | 1,152 | 4,970 |
| human MCL0 (26 copies) | 9,715 | 4,958 | 413 | 4,287 | 57 |

⛔ **On gorilla, 99.7 % of the molecules O2 "assigned" were never ambiguous** — one placement, or a clear
winner. They inflate every assignment rate the project has quoted. ⚠ On human it is ~half. The two substrates
disagree sharply, which is itself a reason to gate rather than to keep quoting a pooled rate.

⚠ **Scope of that measurement**: AS evidence is computed over the reads in the swept REGION, so placements on
other contigs (gorilla) or outside the copy intervals (the human subset BAM) are not counted. "Single
placement" therefore means *single placement seen in this region*, not genome-wide uniqueness.

## The rule to ship
A molecule is **O2-eligible** iff, among its placements inside the family's candidate loci, it has
**≥ 2 placements whose AS ≥ `--as-ratio` × best AS** (existing flag, default 0.98; 1.0 = exact tie).
Everything else does not enter O2: no certificate, no assignment, and — critically — **not in the rate
denominator**.

- new flag **`--as-tied-only`**, default **OFF** ⟹ every existing output stays byte-identical (default flips
  are the user's decision, and this one changes every headline number).
- ineligible molecules are **still emitted** with status **`unambiguous`** and their placement copy, so nothing
  is silently dropped and the count stays auditable. They are excluded from `assigned/tied/ambiguous` rates.
- the stderr rate line gains the eligible count as its own denominator.

## Predictions (write down before running)
| # | prediction | refuted by |
|---|---|---|
| **P1** | The gorilla assigned count collapses by **> 90 %** (6,142 → order 10–500), because only 20 of its assigned rows are ties or near-ties | a drop < 50 % |
| **P2** | The human assigned count falls by roughly **half**, not by 90 % | a fall outside 30–70 % |
| **P3** | Agreement on the surviving assignments **holds or improves** — the certificate was never using AS, so removing easy reads should not break it | agreement falls |
| **P4** | The abstention rate **rises** on both arms: the eligible set is exactly the hard set | abstention falls |
| **P5** ⚠ self-check | Some copies will lose ALL their support and become unsupported, since 4,970 gorilla assignments were single-placement reads | no copy loses support (would mean the gate is not biting) |

## ⛔ What would make this the wrong change
If, after gating, **abstention is so high that no copy retains a certificate**, then O2 as specified has no
addressable population on this substrate and the honest output is aggregate, not per-copy — the K = 0 frontier
argument, applied family-wide. That outcome must be reported, not tuned away by loosening `--as-ratio`.

## Rules held
⚠ Human and gorilla are never pooled. ⚠ The gate is measured on both arms and on the **held-back**
`fam_MCL2_073244` before any headline is restated. ⚠ Byte-identical escape verified by a no-flag re-run.
