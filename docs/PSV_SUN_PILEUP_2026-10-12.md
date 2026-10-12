# SUNs vs our PSVs (M1) + per-copy pileup check (M2) — results, gorilla KB3781 (2026-10-12)

Prereg `docs/PREREG_psv_sun_pileup_2026-10-12.md` (3be23dc3) + Amendment 1 (b658a58e, before any M1/M2 output).
Code `bench/psv_sun/` (6 unit tests); outputs `/mnt/linuxdisk/tmp/psv_sun/` (`sanity.json`, `gate_v.json`, `m1.json`,
`m1_families.tsv`, `m2_pairs.tsv`, `m2.json`, `truth/`). Subset S = 56 families / 208 copies of the 10-10 de novo catalog.
copy_assign snapshot (md5 8b0026b9) `--dump-psv`, 290 s / 9.5 GB under the heavy lock; every other step a light job.

## Finding 0 (gate V, unplanned) — many of O2's PSV columns are not homologous positions

- S1 PASS: 99.6% of the 10-10 rows reproduced.
- Gate V FAIL: only 54.6% of 216,593 (copy, column) pairs verify by independent local alignment (minimap2 asm20,
  either strand, all pieces) + genome base check.
- **10 / 56 families contain copies whose spliced transcripts are reverse-complement homologs** (strand calls on
  1-2-exon de novo copies). O2's own dump shows it: at O2's columns these copy pairs agree at median 0.263 (20 pairs,
  chance level) vs 0.508 for same-strand pairs (95 pairs); decisive-column density 723 / kb vs 9.8 / kb. O2's global
  copy-vs-copy alignment pairs NON-HOMOLOGOUS bases there, so those "PSVs" are alignment artifacts (they inflate
  n_decisive and certificate confidence).
- In strand-consistent families 37% of pairs are not covered by any local alignment (forced global alignment of
  non-homologous segments, or divergence beyond asm20 — not separated here).
- Only **5,235 of 47,134 (11.1%)** decisive columns are homologous in the strict sense (every carrying copy
  verified); strictness is part of that number (one unverified copy excludes the column).

## M1 — PSVs vs SUNs (homologous columns)

| quantity | value |
|---|---|
| homologous decisive columns | 5,235 |
| ... with >= 1 copy-private allele | 5,069 |
| ... with >= 1 SUN (private allele with a genome-unique 30-mer) | 1,594 (30% of homologous PSV columns) |
| private alleles / SUNs | 6,122 / 2,000 (33%) |
| G = private alleles with no genome-unique 30-mer | pooled 0.673; **per-family median 0.125** (43 families with >= 5 private) |

- G is concentrated: MCL116 (1,131/1,162), MCL316 (860/860), MCL38 (654/678), MCL90 (508/609) = 76% of the
  non-unique ones; mostly NOT soft-masked (3,471 unmasked vs 651 masked) — the allele + context exists elsewhere in
  the genome outside the expressed catalog (unexpressed copies the de novo catalog cannot contain).
- Prereg reading: G >= 0.05 -> recommend annotating O2's PSVs with genome-uniqueness. Holds per family (median .125)
  and pooled.
- **R not measurable:** none of the 529 O2-assigned rows covers a homologous column (reported as undefined, not 0).

## M2 — per-copy pileup at PSV columns vs KB3781's other haplotype

- 18,168 verified (copy, column) pairs; 5,370 covered by >= 10 reads (MAPQ>=20 primaries; O2-assigned reads add
  none here). **66 flagged**, 5,270 control.
- Truth: flagged 61/62 with a counterpart = B_ALT (**0.984**); control B differs 10/4,734 (**0.0021**); fold 469.
  **Bar PASS.** NO_COUNTERPART 4 flagged / 456 control.
- 58 of the 61 confirmed flags: the alternative base IS ANOTHER COPY'S PSV ALLELE (3 novel) — on one haplotype copy X
  carries copy Y's allele (median alt fraction 0.47 = heterozygous). Reads from that haplotype give O2 evidence for
  the wrong copy at that column. 30 copies in 17 families.
- No-second-haplotype chromosomes (chrX 4 pairs, chrY 2, mito 74): 0 flags / 80 (chrX alone underpowered).
- Flag transitions 68% (germline-like Ts/Tv; the genome confirms them, so not editing).

## Answer to the advisor

- Pileups do not replace alignment for FINDING PSVs, but a per-copy pileup after placement is a real, cheap check:
  1.2% of covered PSV columns in this animal are heterozygous (one haplotype carries another copy's allele), and the
  check is right 98% of the time against the animal's own haplotype assembly.
- SUN != old name for PSV: SUNs are the genome-unique, copy-private subset (~30% of homologous PSV columns here).
- The bigger problem the check exposed is upstream: O2's PSV columns include non-homologous artifacts (reverse-strand
  copies, forced global alignment). Proposed register rows (not added — register shared with another session):
  (1) "O2 PSV columns are homologous positions" — 54.6% of pairs verify; 10/56 families reverse-complement artifacts.
  (2) "A per-copy pileup flags real deviations from the reference PSV allele" — 61/62 vs 0.2% control, PASS.
