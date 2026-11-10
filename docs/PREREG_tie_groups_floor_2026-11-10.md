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

## Outcome (2026-11-10, job 4456389, `/anvme/workspace/v138da13-rustle/r1603/floor{1,3,5}.*`): all arms FAIL as registered
| floor | loci / groups | B1g grouped reach (>= 58) | B2f purity (>= 0.90) | B3 TAF11L (>= 20) | B5 largest group (<= 60) |
|---|---|---|---|---|---|
| 1 (control) | 5,541 / 2,599 | 70 | 201/234 = 0.859 | 31 | 831 |
| 3 | 2,838 / 450 | 68 | 113/145 = 0.779 | 31 | 719 |
| 5 | 2,154 / 308 | 62 | 95/125 = 0.760 | 31 | 543 |

- B4: 1,359 s - pass. The control reproduces r1605.
- The floor removes small pure groups (judged 234 -> 145 -> 125) and leaves the impure ones (33 -> 32 -> 30), so the purity
  ratio falls. The r1605 post hoc reading (a group of tiny loci) was wrong: the largest group survives the floor.
- Lost grouped reach: at 3 DEFB103A, OR4F5, PPIAL4E, PRAMEF7; at 5 also CFC1, CFC1B, DEFA1, DEFA1B, DEFB105A, MBD3L4.
- Post hoc (not registered): the largest group at floor 5 has 543 loci and 5,172 molecules, median 21 molecules and ~2 kb
  per locus; 377 loci overlap no CAT gene and the rest sit inside unrelated long genes (PRKCA, CADM2, LSAMP, AUTS2, 124
  protein-coding hits). That is a transcribed interspersed repeat, reached from the panel through one POTEF locus. The
  about 30 impure groups are stable under every rule tried (tie width, link fraction, floor).
