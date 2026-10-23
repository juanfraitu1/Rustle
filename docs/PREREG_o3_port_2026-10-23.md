# Pre-registration: O3 port (Strainy's phasing + IsoCon's test), gorilla KB3781 (2026-10-23)

Committed before the port has been run on any panel. Design: `docs/superpowers/specs/2026-10-23-o3-strainy-isocon-design.md`.
The bars are the shipped stage's own bars (`candidate_copies`, Amendments 13-15: register rows 1221-1231 and 1247), so that the
port can replace the stage on equal terms.

## Arms, the same reads and the same scorer
- **P:** the port, run as designed.
- **K:** `candidate_copies` as it ran under Amendment 15. Its products are already on disk: a13a15 / a14a15 / a15h. It is not
  re-run.

## Decision rules
- **P1, development (deletion panel, a13a15).**
  - **Pass:** D right >= 0.80 x C = 5,842 (C = 7,303, Amendment 13b's comparator) AND false moves <= 5% of S reads.
  - **K on the same panel:** 7,882 and 0.77%, which passes.
- **P2, development (no-deletion control, a14a15).**
  - **Pass:** families with >= 1 false flag (class b / c / pri against mat / pat, Amendment 9's rule) <= 8 of 53.
  - **K:** 10 of 53, which fails.
- **P3, held out (a15h; run ONCE, after P1 and P2, with frozen code).**
  - **Pass:** families with >= 1 false flag <= 4 of 30 (H1) AND false moves <= 5% (H3).
  - **K:** 6 of 30 and 0.26%; H1 fails.
- **Adopt the port as the O3 stage iff P1, P2 and P3 all pass.**
  - P passes P1 but fails P2: "detects but not specific", not adopted.
  - P fails P1: "not sensitive", not adopted.

## Rules for development
- P1 and P2 are development. A change after looking is allowed only if it adds no new constant, is recorded in the doc with its
  reason, and is re-run on BOTH P1 and P2.
- Then the code is frozen (commit hash written here) and P3 is run once.

## Reported, not decided
- The per-copy haplotype count in the control against mat / pat: allele-level truth for every copy.
- How many certified flags depend on an editing-type column.
- Runtime.
- HG002 SPDYE1 (development case, human).
