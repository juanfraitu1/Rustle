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
