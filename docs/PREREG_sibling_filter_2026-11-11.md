# PREREG r1629: sibling concordance as a filter on the weakest default chains at family members (2026-11-11)

Goal (user): the family-wide assembly should share information across members. r1627/r1628 showed a sibling's model
is a valid second witness for a one-read chain but adds few chains. The same witness can be asked of chains the
assembler already emits on its floor. Written before any run.

## Rule D
At a family member B onto which at least one sibling model projects (r1519 parameters; B is "informed"), a default
multi-exon chain carried by exactly 2 reads is dropped when no sibling model projects onto B with exactly its intron
chain. Chains with 3 or more reads, single-exon models and members with no projection are untouched.

## Measure (correct = exact RefSeq intron chain; censat-excluded families out)
At informed members: D = the dropped chains, K = the 2-read chains kept (a sibling agrees).
- B1: share correct in D <= 0.5 x share correct among all 2-read chains at informed members (the filter separates).
- B2: correct chains lost (in D) <= 2% of all correct default chains at informed members.
- B3: precision of all default chains at informed members rises by >= 0.01 after the drop.
Step 0 on gorilla KB3781 (stop if any bar fails); then A119b and OR6737 held out, same bars. PASS on all = worth an
opt-in filter on `PREFIX.final.gtf`; FAIL = recorded, not ported.

## Tools
`bench/project/project_models.py` writes `concordant.tsv` (default chains a sibling projects to) and `informed.tsv`;
`bench/famcmp/sibling_filter.py` scores. No Rust.

## Outcome, step 0 (2026-11-11, job 4465569, `r1629/kb`): STOP as registered
Gorilla KB3781, default chains at informed members:

| set | chains | correct | share |
|---|---|---|---|
| all default chains | 3,123 | 614 | 0.197 |
| 2-read chains | 945 | 65 | 0.069 |
| D: 2-read, no sibling agrees (dropped) | 581 | 33 | 0.057 |
| K: 2-read, a sibling agrees (kept) | 364 | 32 | 0.088 |
| 3+ reads, a sibling agrees | 921 | 268 | 0.291 |
| 3+ reads, no sibling | 1,257 | 281 | 0.224 |

- B1 FAIL (0.057 against a bar of 0.034), B2 FAIL (33 correct lost against 12), B3 PASS (+0.032, which any cut of
  2-read chains would give). A119b and OR6737 not run.
- A sibling's agreement moves a 2-read chain from 0.057 to 0.088 and a 3+ read chain from 0.224 to 0.291: a real but
  weak signal. Half of the correct 2-read chains have no agreeing sibling (member-specific isoforms), so absence of
  agreement cannot be used to drop a chain.
