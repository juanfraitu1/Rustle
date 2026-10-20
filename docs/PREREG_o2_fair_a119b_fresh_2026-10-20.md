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
