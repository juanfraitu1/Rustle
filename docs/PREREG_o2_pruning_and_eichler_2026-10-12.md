# O2 on the coin-toss genes, O2-pruned good pool, and Eichler's AS-margin rule, 2026-10-12

**Pre-registered before any S2 read is tested and before any pruned arm exists.** Follows register 1269 (secondary-only genes S1:
1 O2-supported gene per species; most are contradicted — paralog shadows). User (2026-10-12): "lets try on coin toss genes, and also
drop the transcripts in contradicted ones; also can we compare this with Eichler's rule that the copy belongs to the alignment that
is within 10 units of the AS?"

## 1. The certificate on S2 (coin-toss-starved genes)

`bench/tool_difference/s1_o2.py --stratum S2`, otherwise exactly as `docs/PREREG_s1_o2_certificate_2026-10-12.md` (same pairs: each
good secondary record overlapping the gene's exons vs the molecule's primary; same certificate; same gene classes). The S1 numbers
must reproduce exactly with the extended script (`--stratum S1`, the default).

## 2. Eichler's AS-margin rule beside O2 (per record, both strata)

As implemented 2026-09-21 (`docs/archive/2026-09/EICHLER_COMPARISON_2026-09-21.md`, register 935-936): a molecule is assigned to its
best-AS alignment iff no other alignment of it scores within 10 AS units (best - second >= 10, or no second alignment); otherwise it
is discarded. Per tested record (secondary in the stratum gene = copy B, primary = copy A): **B** if B's AS is the molecule's best and
the margin >= 10; **A** if the primary's AS is the best and the margin >= 10; **other** if the best is a third placement with margin
>= 10; **discard** if the margin < 10. Eichler gene classes, parallel to O2's: supported >= 2 records to B; contradicted 0 to B and
>= 2 to A or other; undecidable otherwise. Reported: per-record agreement where both rules call (O2 A/B vs Eichler A/B), and the
gene-class cross-table.

## 3. Pruned good-pool arms

For a rule R (O2 or Eichler), `OP_good_R` = `OP_good` minus every transcript that (a) touches an exon of a gene R calls contradicted
(S1 or S2) and (b) is GOOD-POOL-ONLY: its intron chain is absent from `OP_primary` on that chromosome (mono-exon: no `OP_primary`
transcript overlaps it). Transcripts built from the gene's own primary reads stay. Scored with `bench/tool_difference/score_op.py`
(same reference, strata and development-chromosome exclusion as register 1268), beside the existing arms.

## Claims and expectations (both species)

- **E1 (S2 certificate).** Written now: tied secondaries in S2 are certified to the primary copy more often than to the S2 copy.
- **E2 (pruning, annotation metrics).** Expected: genome-wide metrics move by < 0.5 point; stratum precision moves little, because a
  shadow can match the annotation as well as a real transcript (register 1269: reference-exact S1 transcripts sit mostly in
  contradicted genes). An annotation-only truth cannot score copy-of-origin; the number of shadow transcripts removed is the effect.
- **E3 (O2 vs Eichler).** Expected: where both call, they agree mostly (09-21: 97.1%); Eichler discards most of these records
  (margins are small among good secondaries), so it contradicts fewer genes and prunes fewer shadows.
Nothing is adopted: the outcome is descriptive, for the advisor's figure discussion. A test of copy-of-origin needs a truth with known
source copies (the family simulator, `bench/famsim`), named here as the follow-up.

**Instruments** (committed before running): `s1_o2.py` sha1 a79218c0bc6aa231f40cb7eaca7bca5da1ad1d12 (adds `--stratum` and the Eichler columns; 9 tests), `prune_good.py` sha1 4d17bbe843d8e99c4f406de683e9019abdfe37ca (2 tests).
