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
