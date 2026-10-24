# Pre-registration: O3 putative copies on a fresh gorilla fibroblast held-out (P6), 2026-10-23

**Status:** committed before the panel is drawn; no data of its families has been aligned or scored.
**Goal milestone 2** (`ThesisVault/Topics/Goal_O2_O3.md`): "copies that fit a free allele slot are caught at a pre-registered
precision".

## What is tested
**Putative copies** (`o3port.py`, `putative_copies`) are the uncertified tier: haplotypes that hold an allele slot, so RNA cannot
prove them extra (Theorem A), but where either:
- **(tie)** >= half of the haplotype's reads tie (within 10 AS of their best) to another copy too, while < half of the reads of
  the slot copy's other holders do; or
- **(splicing)** >= 2 of its majority junctions are non-canonical on the slot copy's DNA.

A putative copy needs >= 6 reads and a co-holder on its slot copy.

**Development** (P1 / P2; P5 post hoc for the splicing threshold):
- P1:
  - 11 deleted copies end up in a survivor's free slot; 5 of them are flagged;
  - the flags are 5 D-majority, 2 on GWFAM175's real paternal-only copy, and 1 false (GWFAM415).
- P2: putative flags only on GWFAM175's real copy.

## Panel (fixed before drawing)
80 fresh two-copy fibroblast families, all outside the 53 development families, refabsent's 34, and every family ever aligned or
scored:
- **the 40 families of P4's draw that P5 did not use** (P4's draw: `docs/PREREG_o3_two_tier_p4_2026-10-23.md`, never scored);
- **40 more**, drawn with `random.Random(1).sample` from the eligible families P4 did not draw:
  - eligibility is P4's rule: `intervals.tsv`, every copy listed, >= 2 copies, >= 20 clean reads per copy, clean span <= 200 kb,
    G3;
  - the list is sorted by family.

The rest is as in P5:
- the deleted copy is Amendment 7's;
- the reads come from `link_test.py reads` on this panel (<= 500 primaries per copy, seed 1);
- the deletion arm has only these 80 copies masked; the control arm is the unmasked `_pri`;
- the stage is `p4/bin/candidate_copies` (sha1 2e3cd592), and the port's nets and IsoCon side come from it.

**Port:** `bench/o3_graph/o3port.py` as of 29d63d5b, default options (correlation clustering, putative copies with both signals).

## Decision rules (putative copies)
- **HP1, precision (deletion arm).** True putative flags / all putative flags >= 0.50.
  - A flag is **true** when most of its reads come from the deleted copy (D-majority), or when its consensus classifies as class a
    (a haplotype-only locus: a real reference-absent copy) by Amendment 9's rule.
  - Development P1: 7 / 8.
- **HP2, specificity (control arm).** Families with >= 1 putative flag in class b, c or pri <= 15% of the panel (12 of 80).
  - Development P2: 0.
- **HP3, recall (deletion arm).** Of the deleted copies whose largest D-majority haplotype (>= 6 reads) holds a survivor's
  allele slot uncertified, the fraction flagged putative >= 0.25.
  - Development P1: 5 / 11, with both signals. Two-copy families have one survivor, so only the splicing signal can fire.
- **Verdict.** The putative tier passes iff HP1, HP2 and HP3 pass. HP1 or HP3 with nothing to count (no flag, or no free-slot
  copy) does not pass: underpowered.
- The run is made once; nothing is changed after looking.

## Reported, not decided
- The certified tier on the same panel: certified flags, their classes on the control, and each deleted copy's outcome. No
  realignment scoring.
- Runtime.

## Result (run once, 2026-10-23; port 29d63d5b; work dir `/mnt/linuxdisk/tmp/rna_allele/p6`)
- **Reads and nets.** 49,778 reads (23,011 D). The stage's nets were replicated on both arms (80 / 80 families).
- **Certified flags** (reported): 40 on the deletion arm, 4 on the control arm.

| rule | result | bar | verdict |
|---|---|---|---|
| HP1 precision (deletion arm) | 2 / 2 putative copies are D-majority (GWFAM155 266 reads, GWFAM416 485) | >= 0.50 | PASS |
| HP2 control families with a false putative flag | 1 (GWFAM405, c_unmatched) | <= 12 | PASS |
| HP3 recall over free-slot deleted copies | 2 / 13 = 0.15 | >= 0.25 | **FAIL** |

**Verdict: the putative tier DOES NOT PASS** (not sensitive enough). What it flags is right (2 / 2) and it rarely flags wrongly
(1 / 80 families on the control). But in two-copy families only the splicing signal can fire, and 11 of 13 deleted copies that
sit in a survivor's free slot splice canonically there. Milestone 2 stays open.
