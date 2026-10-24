# Pre-registration: O2 against minimap2's trusted calls on fresh A119b reads (2026-10-20)

**User (2026-10-20):** "no lets focus on redoing it on CHM13 with A119b". That replaces the cerebellum donor of
`PREREG_o2_fair_heldout2_2026-10-19.md`: its whole-genome alignment was stopped by the session for low memory.

**What is held out:** the reads, not the donor. A119b was O2's development library:
- version 2's validated-only rule was chosen after register 1339's errors on A119b fragments;
- k = 3 was chosen on the same fragments;
- the fair scorer's development numbers came from them (232 placed, 1 wrong, against the baseline's 103 wrong of 773).

This draw takes **none of those 15,573 source reads**. Every fragment, and every validating read's role, is new.

## Design (register 1341's, unchanged except for the draw)

- **Library:** A119b testis, `winloci_data/A119b.t2t.bam` (CHM13 v2.0; minimap2 splice:hq -uf --eqx -Y -N 50 -p 0.1).
- **Sources:** up to 60 per pool family, MAPQ-60 primaries of >= 1 kb, **seed 2**, excluding `o2frag_a119b/sources.tsv`
  (`o2frag.py sources --exclude`, 4a57339f).
- **Fragments:** 3 per source, 300-1,200 bp, seed 2.
- **Realignment:** with A119b's parameters plus `-K 100M` (batch size: memory only). Written unsorted, then sorted.
- **Truth:** the source's copy when its full read is clear by 30 AS; scored only when that copy is one of the fragment's
  tied copies.
- **O2:** `o2psv.py svp --k 3 --validated-only` (7bf68347). The validators are A119b's own reads, a fragment's own source
  left out.
- **Scorer:** `o2perfamily.py fair` (fad7e43b) on the 281 qualifying non-pseudoautosomal families, by source read, against
  minimap2's primary at MAPQ >= 1.

## Bars (unchanged)

- **F1:** at least 100 source reads placed.
- **F2:** at most 2% wrong.
- **F3:** precision at least 5 points above the baseline's.
- **F4:** more than 10 families with at least 5 placed and none wrong.

**Verdict:**
- **PASS:** F1 to F4 hold.
- **FAIL:** F1 holds and another bar fails.
- **UNDERPOWERED:** F1 fails.

**A PASS would claim:** on A119b reads O2 has never seen, O2 places almost-tied reads more precisely than minimap2's own
trusted calls, in more than 10 near-identical families.

**It would not claim:** another donor (registers 1341 and the cerebellum attempt), or exact ties. Exact ties carry no
separating PSV; register 1340.

**Rules:** code frozen at 4a57339f. Run script: `o2frag_a119b2/run.sh`. A bug is recorded as an amendment before any rerun.

## Results (`o2frag_a119b2/run.log`, `fair.json`)

**PASS: F1 to F4 hold.**
- **F1:** 197 source reads placed.
- **F2:** 0 wrong. With 0 of 197, the 95% upper bound on the error rate is 1.9%.
- **F3:** precision 1.000, against the MAPQ >= 1 baseline's 0.845 (671 placed, 104 wrong) on the same fragments.
- **F4:** 18 families with at least 5 placed and none wrong.

| fresh A119b draw | value |
|---|---|
| sources drawn / kept (clear by 30 AS) | 11,421 / 11,280, from 212 families |
| fragments with a truth | 33,840 |
| O2's decisions on almost-tied fragments | validated 277, sets 656, unexplained 859 |

**The 18 families (all members at 90-100% identity):**

| family | member identity | O2 right / wrong | minimap2 at MAPQ >= 1: right / wrong |
|---|---|---|---|
| AC1143242 | 0.991 | 6 / 0 | 8 / 2 |
| AC1378001 | 0.998 | 8 / 0 | 12 / 0 |
| BMS1P1 | 0.985-0.998 | 5 / 0 | 9 / 5 |
| CCT6P3 | 0.988 | 6 / 0 | 9 / 3 |
| CLEC18A | 0.971-0.9996 | 9 / 0 | 13 / 0 |
| GPAT2 | 0.992 | 10 / 0 | 18 / 0 |
| GPR89A | 0.971 | 11 / 0 | 16 / 0 |
| GTF2IRD2 | 0.988-0.993 | 8 / 0 | 23 / 5 |
| HSPA1A | 0.996 | 10 / 0 | 36 / 2 |
| HYDIN | 0.985 | 8 / 0 | 14 / 0 |
| MSTRG596 | 0.986 | 6 / 0 | 8 / 0 |
| NOMO1 | 0.996-0.998 | 5 / 0 | 13 / 1 |
| NOS2P2 | 0.997 | 5 / 0 | 9 / 8 |
| PDE4DIP | 0.911-0.998 | 10 / 0 | 14 / 1 |
| SPATA31A1 | 0.929-0.996 | 10 / 0 | 22 / 0 |
| SPATA31D4 | 0.906-0.977 | 6 / 0 | 11 / 1 |
| UGT2B17 | 0.966 | 6 / 0 | 7 / 0 |
| ZNF658 | 0.988 | 8 / 0 | 19 / 1 |

No qualifying family has an O2 error.

## Audit (metric-auditor, same day): holds with caveats

**Recomputed independently:**
- Truth and O2's rows match exactly: 197 placed, 0 wrong; F4's 18 families.
- The baseline is 670 placed, 103 wrong (0.846). That is one source off, because the row is filed under the primary's
  family. It does not flip.
- By fragment: O2 220 / 220 right; minimap2 at MAPQ >= 1 689 / 801 (0.860).
- The seed-1 and seed-2 sources share no read. The code was unchanged and the pre-registration preceded every output.
  No validation leakage.

**Caveats:**
1. **O2 certifies; it does not overrule.** Its call equals minimap2's primary on all 220 placed fragments, and every one of
   the 112 fragments minimap2 got wrong is one O2 abstained on. So O2 certifies about 29% of minimap2's trusted calls, at
   no error. Nothing simpler does that:
   - no MAPQ cut matches it (MAPQ >= 30: 91 placed, 0.879);
   - neither does a substitution count without the sample's validation (360 placed, 0.908).
2. **Held out by reads only.** The donor, families, PSV sites and validating reads are A119b's, where k = 3 and
   validated-only were chosen. 14 of the 18 F4 families also passed in development. This is a replicate. The other-donor
   test (register 1341) is still UNDERPOWERED.
3. **F4 does not single out O2.** minimap2's MAPQ >= 1 calls meet the same per-family rule in 19 families, and a floor of
   5 cannot bound a family's error rate.
4. **Proxy truth.** The truth is minimap2's call on the full read. A wrong full read would be invisible here; on HG002
   this truth agrees with the donor 98.4-99.8% at this margin.
