# Pre-registration: O3 putative copies with the soft-clip signal on a fresh gorilla fibroblast held-out (P10), 2026-10-23

**Status:** committed before the panel is drawn; none of its families has been aligned or scored.
**Goal milestone 2** (`ThesisVault/Topics/Goal_O2_O3.md`). The user said "work on 2" after P6 failed recall (register 1380).

## What changed since P6
The putative tier now has a third signal, `foreign_clipping` (`o3port.py` e3597fbc, default on): >= half of a slot holder's reads
with a near-best record on the slot copy are soft-clipped there by >= 50 bp that are not a poly-A / poly-T tail, and >= 6 of them
clip at one boundary (same side, same 20-bp bin). The deleted copy's reads stop matching the survivor where its UTR or an exon
diverges; a copy's own alleles align end to end.

**Development (P5 and P6, both spent; looked at after P6's verdict):**
- P6 deletion arm: free-slot recall 2 / 13 -> 5 / 13; 7 / 7 putative flags D-majority.
- P6 control arm: 1 flag (GWFAM405, as before; nothing new). P5 control: 0 haplotypes meet the rule.
- The thresholds (50 bp, half, 6 at one boundary, 70% tail) were set once before looking at P6's numbers and not tuned.

## Panel (fixed before drawing)
Every family eligible by P4's rule (`o3port_p4.py eligible_families`) that is in none of P4's, P5's or P6's panels: 58 two-copy
families, sorted by name. Everything else as in P6: Amendment 7's deleted copy, `link_test.py reads` (<= 500 primaries per copy,
seed 1), the deletion arm masks only these 58 copies, the control arm is the unmasked `_pri`, the stage `p4/bin/candidate_copies`
(sha1 2e3cd592) gives the nets and the IsoCon side. Work dir `/mnt/linuxdisk/tmp/rna_allele/p10`.

**Port:** `bench/o3_graph/o3port.py` as of e3597fbc, default options.

## Decision rules (as P6)
- **HP1, precision (deletion arm).** True putative flags / all putative flags >= 0.50 (true = D-majority, or class a by Amendment
  9's rule).
- **HP2, specificity (control arm).** Families with >= 1 putative flag in class b, c or pri <= 15% of the panel (8 of 58).
- **HP3, recall (deletion arm).** Of the deleted copies whose largest D-majority haplotype (>= 6 reads) holds a survivor's allele
  slot uncertified, the fraction flagged putative >= 0.25.
- **Verdict.** Pass iff HP1, HP2 and HP3. HP1 or HP3 with nothing to count: underpowered (no pass). Run once; nothing changed after
  looking.

## Reported, not decided
Each signal's share of the flags (tie, splicing, clipping); certified flags; runtime.

## Result (run once, 2026-10-23; port e3597fbc; work dir `/mnt/linuxdisk/tmp/rna_allele/p10`)
- Reads 38,675 (19,571 D); the stage replicated on both arms (58 / 58); certified flags 38 (deletion arm), 3 (control).

| rule | result | bar | verdict |
|---|---|---|---|
| HP1 precision (deletion arm) | 6 / 6 putative copies D-majority (GWFAM119 11 reads, 139 88, 171 6, 349 95, 363 313, 74 142) | >= 0.50 | PASS |
| HP2 control families with a false putative flag | 0 (no putative flag at all on the control) | <= 8 | PASS |
| HP3 recall over free-slot deleted copies | 5 / 8 = 0.625 | >= 0.25 | PASS |

**Verdict: the putative tier PASSES.**
- Reported: every one of the 6 flags carries the clip signature (>= 59 % of its reads clipped at one boundary), so the soft-clip
  signal is what fires. Misses: GWFAM380 (454 reads, all clipped when every record is counted, but not on the near-best records
  the port uses), GWFAM259 and GWFAM409 (no clipping: the deleted copy aligns end to end on the survivor).
