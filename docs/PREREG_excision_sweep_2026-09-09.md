# PREREG — excision sweep over every human MCL0 copy (2026-09-09, before any run)

**Question.** §6hc's excision tested copies 2 and 22 only (97.4 % / 100 % abstention). The final pipeline (§6hf)
assigns 230 contested human molecules over 12 copies (2: 153, 22: 42, 12: 16, 9: 5, 10: 4, 19/18/11: 2, 8/7/21/16: 1).
For EVERY copy k ∈ 0..25: remove k from `copies16.tsv/.fa` (rows renumbered contiguously, `remap_no{k}.tsv`
new→old), re-run `copy_assign` with the shipped defaults, and follow (a) the molecules assigned to k in
`ours_final` and (b) everything else. Control (§6ff pattern): a molecule assigned to a copy that no longer exists
must abstain; assigning it elsewhere is a silent wrong call.

Scoring: `bench/o2_excision_sweep.py`. "Abstain" = not `assigned` in the excised run (tied / ambiguous /
origin-rejected / absent). "Moved" = `assigned` to another copy (mapped back through the remap). Sibling
identity = 1 − X/aligned over the genomic unit spans (minimap2 asm20 `--eqx`, all fragments summed, from
`scratchpad/human_gspans.paf`).

| # | prediction | refuted by |
|---|---|---|
| P1 | pooled over the 230: **≥ 95 % abstain** when their copy is removed | < 90 % |
| P2 | every copy with ≥ 5 assigned reads (2, 22, 12, 9) abstains **≥ 80 %** | any < 70 % |
| P3 | moved reads go to a near-identical sibling (**identity ≥ 0.97**) at **margin median < 25**; **0** move with margin ≥ 40 ("silent-confident") | > 30 % of moved to a copy < 0.97, or ≥ 3 silent-confident |
| P4 | in each of the 26 runs, contested molecules whose copy ≠ k keep status AND copy **≥ 99 %** (no copy is load-bearing for the others' verdicts) | any run < 97 % |
| P5 | removing a copy with 0 assignments changes ≤ 1 % of contested statuses | any such run > 3 % |
Human only (gorilla MCL1 has 4 assignments on 2 copies — reported, not predicted). Defaults as shipped
(§6hf); no flags.

## Outcome (2026-09-09) — `bakeoff/human/excise_all/no{0..25}`, scorer `bench/o2_excision_sweep.py`
| # | verdict |
|---|---|
| P1 | ✓ **219/230 = 95.2 %** abstain when their copy is removed (at the boundary) |
| P2 | ⚠ missed, not refuted: copy 12 abstains **12/16 = 75 %** (copies 2: 149/153, 22: 42/42, 9: 5/5) |
| P3 | ⛔ 11 moved, all to a sibling at identity ≥ 0.988 (median 0.992), margin median 20.7 — but **5 silent-confident** (margin ≥ 40): 12→19 ×2 (186, 180), 19→12 (110), 10→11 ×2 (76). Their BASE assignments were strong (n_decisive 11–17, margin 62–104): the read matches copy 12 at ≥ 16 columns where 12 ≠ 19, so against copy 19 it carries ≥ 16 substitutions over ~3 kb — and the origin certificate (Binomial, 0.3 %, alpha 0.001) needs ≈ 19 to reject. **A sibling ≲ 0.6 % diverged over the read footprint is invisible to the certificate** — the same limit that sets O3's ≥ 0.7 % flag |
| P4 | ⛔ refuted in letter: 24 of 26 runs keep < 99 % of the other copies' assignments (copy 6 removed: 84/230). Mechanism (`scratchpad` breakdown over all 26 runs): **272 designed abstentions** (the removed copy was a TIE PARTNER: the tie is now outside the catalog ⟹ forced `tied`, §6gz), **137 rows vanish** (92 molecules; open — see below), **38 events / 4 molecules assigned to a different copy** (copies 10↔11 at identity 0.992, base margin 6.9 = ONE column: they flip whenever the column set changes), 35 ambiguous. A tie partner IS load-bearing, by design; the genuinely unstable class is 4 one-column assignments |
| P5 | ⛔ copies 6 / 20 / 23 / 24 (0 assignments) change 44 % / 10 % / 11 % / 12 % of contested statuses — the same tie-partner mechanism (6/7/8 are the copy-2 population's tie set) |
**The fragile class, quantified:** 19/230 assignments rest on ≤ 2 columns (margin < 14), 6 of them at copy 12, 4 at 9, 3 at 10 — all inside near-identical sibling groups (2/8 at 0.996, 10/11 at 0.992, 12/19 at 0.988). **Open:** 92 molecules (primaries at copies 2/22/12) lose their row entirely when copies 6, 24, 23, 9, 11 or 10 are removed — not a wrong call (absent = abstained) but unexplained; to be traced in the gate before the sweep is quoted as a mechanism table.
