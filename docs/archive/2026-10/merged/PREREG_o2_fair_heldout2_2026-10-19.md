# Pre-registration: O2 against minimap2's trusted calls on a deep new donor (2026-10-19)

**Why:** register 1341 held O2 out on the public testis donor and was UNDERPOWERED: 30 source reads placed, 0 wrong,
against minimap2's MAPQ >= 1 calls at 0.894. That library is shallow, so a PSV rarely reaches 3 validating reads.

This repeats the same test on a deep donor never used in O2's work.

## The donor

PacBio's public Kinnex full-length RNA release, `DATA-Revio-SCRI-Sample5-Cerebellum`
(https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/; released with the collaborator's permission).
- **Sample:** one prenatal cerebellum from Seattle Children's Research Institute (one donor). Revio, Kinnex full-length
  RNA kit.
- **Reads:** 44,962,471 FLNC reads (`2-FLNC/flnc.bam`).

## Design

- **Steps 1-2, as HG002's proof (`hg002proof/screen.sh`):**
  - the FLNC reads are screened against the pool's transcripts (`family_tx.fa`, `minimap2 -x map-hifi`); a hit at
    >= 0.90 identity over >= 200 bp makes a candidate;
  - the candidates are aligned to whole CHM13 v2.0 with A119b's parameters:
    `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`.

  That BAM is this donor's library. Its reads are the sources and the validators.
- **Step 3, exactly as register 1341** (`docs/PREREG_o2_fair_heldout_2026-10-19.md`):
  - **sources:** up to 500 per family, seed 1, MAPQ-60 primaries of >= 1 kb;
  - **fragments:** 3 per source, 300-1,200 bp, realigned with the same minimap2 parameters;
  - **truth:** the source's copy when its full read is clear by 30 AS, the fragment scored only when that copy is one of
    its tied copies;
  - **O2:** `o2psv.py svp --k 3 --validated-only`, a fragment's own source left out;
  - **scorer:** `o2perfamily.py fair` on the 281 qualifying non-pseudoautosomal families, counted by source read,
    against minimap2's primary at MAPQ >= 1.

## Bars (unchanged)

- **F1:** at least 100 source reads placed.
- **F2:** at most 2% of them wrong.
- **F3:** O2's precision at least 5 points above the baseline's.
- **F4:** more than 10 families with at least 5 placed and none wrong.

**Verdict:**
- **PASS:** F1 to F4 hold.
- **FAIL:** F1 holds and another bar fails.
- **UNDERPOWERED:** F1 fails.

**Rules:** code frozen at f2d19688 (regimen 7bf68347, scorer fad7e43b). Run script: `scri5/run.sh`. A bug is recorded as an
amendment before any rerun.
