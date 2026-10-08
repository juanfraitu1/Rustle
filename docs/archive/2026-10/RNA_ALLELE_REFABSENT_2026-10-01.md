# Real reference-absent copies of KB3781: does the RNA-only chain flag them? (Amendment 10), 2026-10-01

Prereg: `docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md` Amendment 10 (commit ee32021d; truth built first, rules written before the
chain ran on these families). Scripts `bench/rna_allele/refabsent_truth.py` (truth), `control_test.py` (the chain, `--l`/`--w` on the
work dir), `refabsent_score.py` (rules); work dir `/mnt/linuxdisk/tmp/rna_allele/refabsent/`; outputs copied to
`docs/RNA_ALLELE_REFABSENT_score.out.txt`.

## The truth (independent of the chain)

- 378 families / 915 copies of the 2026-08-14 interval table; every copy lifted to its B haplotype (913/915 >= 50%) and aligned to both
  haplotype assemblies. **Haplotype-only (B-only) loci** = hits (identity >= 0.90, coverage >= 0.80) on the haplotype `_pri` did not take
  that chromosome from, overlapping the lifted B interval of no copy of the family: **127 loci in 34 families**.
- Back-aligned to `_pri`: **13 loci (6 families) are absent beyond delta** (best `_pri` identity < 0.9904, median 0.986 — detectable in
  principle); **114 loci (33 families) are within delta** (median 0.9986: alleles of `_pri` loci the table does not list, or near-identical
  duplicates — indistinguishable from an allele by construction).
- **Expressed** (>= 3 fibroblast reads best-placed there over both haplotypes, untied): **11 loci in 8 families**. Beyond delta:
  GWFAM175_B0 (281 reads), GWFAM26_B2 (13), GWFAM26_B3 (9), GWFAM175_B1 (6); within delta: 7 loci (4-77 reads). Two more beyond-delta
  loci have 2 reads. **One locus has >= 20 reads**: the biology of one individual and one tissue makes this a demonstration, not a rate.

## The chain on the 34 B-only families

Reads: every baseline record on a copy of the family, <= 2,000 per family (32,219 reads); IsoCon input <= 1,000 (20,844). 1,060 IsoCon
outputs -> 647 flagged (not in `_pri` at 0.999) -> 582 linked back (d <= delta) -> 65 new-copy contigs in 11 families -> **18 candidates;
7 flags (>= 2 transcripts, the floor pre-registered from Amendment 9) in 6 families**.

| locus (truth) | class | reads | flag matching it (best haplotype hit >= 0.999 on the locus) |
|---|---|---|---|
| **GWFAM175_B0** | beyond delta | 281 | **yes — 4 transcripts, paternal chr5 40,028,599-40,031,207, identity 0.9991** |
| GWFAM26_B2 | beyond delta | 13 | no (no candidate, no single transcript) |
| GWFAM26_B3 | beyond delta | 9 | no |
| GWFAM175_B1 | beyond delta | 6 | no |
| 7 within-delta loci | within delta | 4-77 | none flagged (their transcripts link to a `_pri` locus) |

- **D1 (detection) PASSES**: the only expressed beyond-delta locus with >= 20 reads is flagged. The three with 6-13 reads are not: with
  <= 1,000 of the family's reads in IsoCon's input and >= 3 reads per transcript, >= 2 transcripts, a locus needs ~6+ reads IN THE INPUT;
  these had ~3-7 after the sampling.
- **D2 HOLDS**: 0 of the 7 expressed within-delta loci is flagged — delta does what it says.
- **D3 HOLDS**: families with a false flag (class b or c, >= 2 transcripts) 6/34 = 17.6% (bar 20%).
- **The O2 side of the real missing copy (GWFAM175_B0, 281 reads):** in arm R0 (`_pri` only) 262 reads are placed on `_pri` copy
  GWFAM175:2 and 17 on copy 1 at a median `de` of **0.066** — a reference-absent copy's reads forced onto a paralog at 6.6% mismatch, the
  O3 signature; in arm C (`_pri` + candidates) **279/281** have their primary on the candidate (273 tie only among its own 4 transcripts,
  which the component rule places).

## The false flags, read closely (not the verdict)

- 4 of the 6 false-flag families carry a class-(b) "allele" flag with many transcripts — GWFAM256 (16 transcripts, maternal identity
  0.9996), GWFAM318 (11, paternal 0.9995), GWFAM70 (7, maternal 0.9996), GWFAM175:6 (5, paternal 1.0000, the paternal counterpart of
  `_pri` copy 5 whose transcripts are **2.4-4.2% from every `_pri` copy**). These are haplotype loci at the syntenic position of a `_pri`
  copy yet beyond allele divergence from it (otherwise they would have linked): divergent copies that occupy the homologous position
  (tandem-array shuffling), reference-absent in sequence, "alleles" only by position. The lift-based truth classes them as false; by the
  O3 question ("expressed sequence not in the reference") they are true. Amendment 10's rules stand as written; the reading is noted.
- The 2 class-(c) flags: GWFAM169 (2 transcripts, best haplotype hit 0.956, source outside the family) and GWFAM382 (9 transcripts, 0.926
  over a 376 kb spliced hit): transcripts in neither haplotype at consensus accuracy — the same kind as Amendment 9's unmatched flags.

## Reading

- On the real biology of the reference animal the chain finds the one well-expressed reference-absent copy, gives its reads a home at
  0.9991 identity instead of 6.6% mismatch, does not flag the alleles, and raises false flags in 17.6% of these families — half of which
  are divergent syntenic copies that a positional truth calls alleles.
- Power is the limit, not the method: fibroblasts express 4 beyond-delta loci in this individual, one of them well. The matched-truth
  design transfers to any individual with a diploid assembly (HG002 Kinnex + v1.1 is the human instance).
- Caveats: the 1,000-read cap costs the 6-13-read loci; `_pri` chrX/chrY have no B; unplaced haplotype contigs are not searched; the
  truth's lift-based allele/copy split is conservative for tandem arrays.
