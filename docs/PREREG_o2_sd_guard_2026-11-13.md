# Row 1685 — O2 measurement: SD identity between copies as a guard on assignments (2026-11-13)

Goal (user, 2026-11-13): "explore what other ways can we include the segmental duplications so they help the
pipeline"; user's choice "Run the O2 measurement". **Read-assignment layer, which the thesis advisor has not approved:
opt-in measurement only, nothing wired, defaults unchanged, not for `advisor_share`.** O2 applies to AS-tied reads
only; no claim against the aligner at scale is made here.

**Mechanism tested.** Two copies held by a very young segmental duplication differ at few sites, so an assignment
between them may rest on one or two differences and be wrong more often. If so, abstaining on such copies removes
wrong assignments at little cost.

**Data (existing, KB3781 on its own T2T, strict haplotype truth, slack 0).** Row 1514's `--primary-on-family` arm:
development `o2mstie/pof_panel.strict3_dump.tsv` (393 right / 26 wrong), test `o2mstie/pof_heldout.strict3_dump.tsv`
(129 right / 22 wrong). Catalog `~/runs/gorilla_KB3781/kb.fam.copies.tsv`, gorilla SEDEF `sedef_dl/gorilla/final.bed`.
I know these right/wrong counts and that a few families hold most errors; no SD attribute of these rows has been seen.

**Attribute (blind, frozen, `bench/o2_contested/sd_guard.py`).** For an assigned row: the highest SEDEF fracMatch of
an SD pair holding the assigned copy and any other copy of its family (300 bp over each; none = 0).
**Guard arms:** abstain when the attribute >= T, T in {0.98, 0.99, 0.995}. Also reported, no bar: the reverse group
(no SD pair holds the copy with any sibling), in case errors sit in families that are not duplications at all.

**Bars.** Development: an arm that removes >= 50% of wrong rows and <= 20% of right rows; none -> STOP. Test: the
chosen arm meets the same two numbers on the held-out dump. Reported: the families behind the removed wrong rows,
since rows of one family are not independent.
**If it passes:** I report and ask; nothing is wired. **If it stops or fails:** the row closes.

## Result (2026-11-13, workspace `r1685_o2sd/sd_guard.out`) — STOP

Development dump (393 right / 26 wrong); every assigned copy shares an SD pair with a sibling.

| T | wrong removed | right removed | wrong rows by family |
|---|---|---|---|
| 0.98 | 25 (0.96) | 253 (0.64) | MCL5 12, MCL3 10, MCL48 2, MCL234 1 |
| 0.99 | 25 (0.96) | 237 (0.60) | same |
| 0.995 | 11 (0.42) | 138 (0.35) | MCL3 10, MCL234 1 |

No arm removes half of the wrong rows at a cost of at most 20% of the right ones; the held-out dump was not scored.
Reading: the mechanism is half true. Wrong assignments do sit on copies held by young duplications (25 of 26 at
>= 0.99), but so do 60% of the right ones, because that is where AS-tied reads and O2's assignments are in the first
place. Two families (MCL5, MCL3) hold 22 of the 26 wrong rows, so the errors are a family matter (row 1656's
sole-candidate finding), not an SD-identity matter. Row closed; nothing wired; defaults unchanged.
