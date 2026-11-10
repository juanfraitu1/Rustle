# PREREG r1606: tie groups — a floor on tied molecules per locus (2026-11-10)

Follows r1605 (`docs/PREREG_tie_groups_link_2026-11-10.md`, FAIL: purity 0.859, largest group 831 loci; link fraction inert).
Post hoc there: the largest group is 818 loci over every chromosome at ~1.6 molecules per locus. Written and committed
before any arm runs.

## Rule (`tie_groups.py build PREFIX as98 0 FLOOR`)
A locus exists only when >= FLOOR distinct tied molecules have a record in it; smaller loci are dropped before linking and
are not emitted. Arms: FLOOR 1 (control = r1605 fraction 0), 3, 5. Ties as98, link fraction 0, scorer `full` (r1605).

## Bars (per arm, as r1605; M = the 72 indistinguishable members)
- B1g grouped reach >= 58 of 72; B2f purity >= 0.90; B3 TAF11L models <= 5 kb reach >= 20 of 31; B5 largest group <= 60
  loci; B4 the job <= 1 h wall.
- An arm PASSES on all five. If 3 and 5 both pass, the higher grouped reach is kept; a tie goes to the higher purity.
- The control must reproduce r1605 (70 / 0.859 / 31 / 831).
- Descriptive: loci and groups left, and which indistinguishable members lose grouped reach.

## Outcome
(empty until the run)
