# Row 1811: second ablation at the toned-down defaults: primate-only levers and untested levers (pre-registered 2026-11-19)

User (2026-11-19): try to tone down other overfit levers. Row 1793 called two levers overfit (both toned down, row 1810). Its
tables also showed levers that add correct chains on the two tuning samples and almost none elsewhere (transfrag extension: OR6737
+656 / -26, A119b +1,666 / -295; fly +2, mouse simulated +18 / -16, mouse ES +30 / -14) while adding wrong chains everywhere. That
pattern had no name in row 1793's reading. Registered before any arm ran. Nothing here changes a default.

**Arms** (one switch against `base` = today's defaults: junction snap on, sub-chain ratio 0.3, chain score cut 0.02; `RUSTLE_PTC=0`):
`ext0` (`RUSTLE_EXTEND_TRANSFRAGS=0`), `seed_primary` (`RUSTLE_SEED_POOL=primary`), `ext0_seed` (both), `bridge_off`
(`RUSTLE_BRIDGE_REGROUP=off`, untested so far), `layers_off` (`RUSTLE_ASSEMBLY_LAYERS=off`, untested; acts on `final.gtf` only).
**Data.** Tuning: OR6737, A119b. Non-primate: fly, mouse, zebrafish simulated; mouse ES, mouse heart. SIRV shown.
**Unit.** Exact intron chains; `PREFIX.gtf` (`layers_off`: `PREFIX.final.gtf`). Lever effect = base against the arm without it.

**Reading, fixed now.** On a dataset the lever *gives* when it adds (gained minus lost) at least 1% of the arm's correct chains.
- **PRIMATE-ONLY** = gives on both tuning samples; on at least 4 of the 5 non-primate datasets it does not give AND adds more
  wrong chains than correct ones.
- **GENERAL GAIN** = gives on both tuning samples and on at least 3 of 5 non-primate datasets.
- Otherwise described per dataset. The row 1793 reading (pays = higher precision at <= 5% cost) is shown beside it.
A PRIMATE-ONLY lever is put to the user as a candidate to switch off or restrict; the user rules.

**Limits.** Levers without a dial can only be on or off. RefSeq / GENCODE incomplete; simulated reads lack truncated molecules, and
extension exists to join truncated pieces, so a null effect on simulated reads is expected by design: the two real non-primate
libraries (ES, heart) carry the weight. Row 1801 (another session) tests duplicate counting, spliced share and the shadow ratio;
not repeated here.

## Result (2026-11-19): extension PRIMATE-ONLY; seed pool gives nowhere; two arms without effect

Jobs 4482800-4482806, every arm exit 0. Tables: workspace `r1811_ablation2` (`score.*.txt`, `wrong.*.txt`).
What the lever adds to `PREFIX.gtf` (base against the arm without it; net correct chains, share of the arm's correct chains, net
wrong chains):

| dataset | transfrag extension | seed pool (good against primary) |
|---|---|---|
| A119b | +1,435 (3.6%) / +16,507 | +276 (0.7%) / +12,446 |
| OR6737 | +644 (2.5%) / +2,526 | +162 (0.6%) / +1,093 |
| mouse heart | +47 (0.4%) / +311 | +24 (0.2%) / +198 |
| mouse ES | +19 (0.1%) / +1,009 | +54 (0.3%) / +1,834 |
| zebrafish simulated | +13 (0.1%) / +661 | +30 (0.2%) / +498 |
| mouse simulated | +8 (0.0%) / +987 | +38 (0.1%) / +221 |
| fly simulated | +3 (0.0%) / +40 | +4 (0.1%) / +784 |

Gained and lost separately (arm against base; correct +/-, wrong +/-): `score.*.txt`, e.g. A119b ext0 +244 / -1,679, +8,690 / -25,197.
- **Transfrag extension: PRIMATE-ONLY.** Gives on both tuning samples; on 5 of 5 non-primate datasets it does not give and adds
  more wrong chains than correct ones.
- **Seed pool:** gives on no dataset, the tuning samples included, so not PRIMATE-ONLY by the reading; it adds 4 to 45 wrong
  chains per correct one. Exact chains genome-wide do not measure what the pool was made for (copies of multi-copy families).
- **`bridge_off`:** the same chains as base on all seven. **`layers_off`:** writes no `final.gtf`; the arm cannot be scored as set up.
- SIRV (mouse ES): base 110 chains / 52 correct; without extension 90 / 53.

Wrong chains of base `PREFIX.gtf` by class (`wrong.py`, first class that fits; zebrafish simulated 4,034 / mouse ES 39,287 / SIRV 58):
piece of a true chain 1,951 / 3,637 / 6; no true intron 691 / 4,315 / 4; true introns in a combination not in the truth 669 /
9,249 / 17; new intron elsewhere 403 / 16,848 / 1; new intron within 10 bp of a true one 320 / 5,238 / 30. Without extension 500
of the 669 and 16 of the 17 false combinations go; with primary seeds 410 of the 691 "no true intron" chains go (zebrafish).

**User decision (2026-11-19, register 1812):** extension is opt-in again (`RUSTLE_EXTEND_TRANSFRAGS=1`), and the stage `assemble`
on its own seeds with primary alignments (`--seed-pool good` or `RUSTLE_SEED_POOL=good` = before; `all` and the family stages
keep good). Both together = arm `ext0_seed`.
