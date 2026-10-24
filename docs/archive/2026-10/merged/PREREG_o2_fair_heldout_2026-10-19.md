# Pre-registration: O2 against minimap2's own trusted calls, held out on a new donor (2026-10-19)

**Why:** register 1340's audit voided the per-family claim on two counts.
- The data was not held out.
- O2's precision was set against minimap2's on all reads, and minimap2's primary at MAPQ >= 1 passed that rule too.

O2's measured value is precision on almost-tied reads. On HG002 and A119b (development), O2's placements were 1.1% and 0.6%
wrong, while minimap2's MAPQ >= 1 calls on the same reads were 27% and 14% wrong. This test states that claim against the
equally selective baseline, on a donor never used in O2's development.

## Design

- **The regimen:** `o2psv.py svp --k 3 --validated-only` (7bf68347), unchanged.
- **The donor:** the public testis library (ENA ERR13885926, GENCODE CapTrap; `_from_wsl/human_val/human_testis.t2t.bam`).
- **Sources:** per pool family (`r19seq`), up to 500 MAPQ-60 primaries of at least 1 kb on a copy (`o2frag.py sources
  --per-family 500 --seed 1`).
- **Fragments:** 3 per source, 300-1,200 bp. Realigned with the library's own parameters:
  `minimap2 -ax splice:hq --eqx -Y -N 50 -p 0.1 --secondary=yes`.
- **Truth:** the source's copy when the realigned full read is at least 30 AS clear of any other locus. A fragment is
  scored only when that copy is one of its tied copies.
- **Validation:** the library's own reads. A fragment's own source never validates it.
- **Families:** the 281 qualifying families (every member pair at identity [0.90, 1.00)), with the 17 pseudoautosomal
  families left out (`famident/qualify_nonpar.json`).
- **Counting:** by source read. A source counts as right for a method when every fragment the method places is on its
  copy (`o2perfamily.py fair`, e1acb42d).
- **The baseline:** minimap2's primary for the same fragments, taken only at MAPQ >= 1.

## Bars

- **F1, powered:** O2 places at least 100 source reads.
- **F2, precise:** at most 2% of them are wrong.
- **F3, adds precision:** O2's precision is at least 5 points above the MAPQ >= 1 baseline's on the same fragments.
- **F4, the user's handful:** more than 10 families where O2 places at least 5 source reads and none wrong.

**Verdict:**
- **PASS:** F1 to F4 hold.
- **FAIL:** F1 holds and another bar fails.
- **UNDERPOWERED:** F1 fails. The library is shallow (28,662 MAPQ-60 reads on the pool's copies), so this is possible.

**Development, for reference only (A119b fragments, already seen):** O2 placed 232 source reads, 1 wrong (99.6%); the
MAPQ >= 1 baseline placed 773, 103 wrong (86.7%); 18 families pass F4's rule.

**Simulation (one genome):** register 1339's test A stands. On 18 fresh worlds, version 1 placed 349 reads, 0 wrong.
Version 2 placed 1,104, 0 wrong, over rounds 23-25. In one genome minimap2's MAPQ >= 1 calls are already 99.86% right,
so the simulation claim is no harm, not a precision gain.

**Rules:** code frozen at e1acb42d. A bug is recorded as an amendment before any rerun.

## Results (`o2frag_public/run_real.log`, `fair.json`)

**UNDERPOWERED: F1 fails, with 30 source reads placed against 100 needed.**
- F2 holds: 0 wrong.
- F3 holds: precision 1.000, against the MAPQ >= 1 baseline's 0.894 (160 placed, 17 wrong) on the same fragments.
- F4 fails: 1 family (ROPN1, 6 of 6).

| public testis library | value |
|---|---|
| sources drawn / kept (clear by 30 AS) | 6,674 from 215 families / 6,070 |
| fragments with a truth | 18,210 |
| O2's decisions on almost-tied fragments | validated 49, sets 153, unexplained 178 |

The library is about 45 times shallower than A119b. Few PSV positions reach 3 validating reads, so O2 places little.

**Where minimap2's trusted calls fail, O2 abstains or is right:**

| family | minimap2 at MAPQ >= 1: right / wrong | O2: right / wrong |
|---|---|---|
| TPSB2 | 13 / 7 | abstained |
| DHRS4 | 5 / 4 | 2 / 0 |
| HSPA1A | 19 / 2 | 2 / 0 |
| THOC3 | 9 / 1 | 2 / 0 |

The direction matches development (A119b: 99.6% against 86.7%), but this donor cannot carry the claim.
