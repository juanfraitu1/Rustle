# PREREG — PSV columns only inside significant local blocks of the copy-vs-copy alignment (2026-10-12)

**Problem** (`docs/PSV_SUN_PILEUP_2026-10-12.md`; register 683). `discover_psvs` aligns each copy's spliced sequence to
copy 0's GLOBALLY (banded Gotoh; a 150-bp gap costs 182, 150 unrelated bases as mismatches ~112), so stretches with no
counterpart are forced into columns. Diagnostic `bench/psv_sun/gap_diag.py` on the 56-family KB3781 subset (HEAD
a37399a5 dump, 126,499 decisive (copy, column) pairs): 54.5% lie in stretches a chained `minimap2 asm20` alignment
supports ("aligned"); 31.6% in stretches asm20 misses but Smith-Waterman beats all 20 shuffles (divergent homologs:
`both_homologous` 21.3%, `no_alignment_homologous` 10.3%); 14.0% in stretches SW does not (`both_unrelated` 7.8%,
`terminal_overhang` 4.3%, `no_alignment_unrelated` 1.9%). Labels per pair: `/mnt/linuxdisk/tmp/psv_sun/gap_diag/pair_labels.tsv`
(frozen before the fix is written). Chained local alignment alone (register 577) would lose the 31.6%.

**Fix (decided before implementation).** Keep the global alignment; a column becomes a PSV column only inside a block
of the alignment path chosen by
  maximize  sum over kept blocks (S_block - T),   T = ln(L / alpha) / lambda
- column scores BLASTN-style: match +2, mismatch (or non-ACGT) -3, gap run of k columns -(5 + 2k);
- L = alignment length (columns), alpha = 1e-3 (the certificate's alpha), lambda solves
  0.25 e^{2 lambda} + 0.75 e^{-3 lambda} = 1 (ungapped Karlin-Altschul, uniform bases; K taken as 1);
- exact O(L) DP over prefix sums; columns outside kept blocks contribute no alignment (the copy is a gap there).
A block must pay its own significance cost, so an unrelated exon between two homologous blocks is cut out when its
negative score exceeds T, while a few PSV mismatches inside a block (a small dip) stay. Applies to every
`discover_psvs` engine (banded DP, poasta, minimap2).

**Acceptance (56-family subset, binary with vs without the fix, identical arguments):**
- **B1** assignments, quant, posterior/abundance byte-identical (the default read-star certificate does not read these
  columns). Otherwise stop and explain.
- **B2** pairs labelled `aligned`: >= 99% retained (same family, same copy-0 genome position, same copy, same allele).
- **B3** pairs labelled `both_unrelated` + `terminal_overhang` + `no_alignment_unrelated`: >= 90% removed.
- **B4** reported, no bar: retention of `both_homologous` and `no_alignment_homologous`.
- **B5** library tests pass (the 3 known `tests/pipeline_driver_dispatch.rs` failures are pre-existing and unrelated).
Fail of B2 or B3 -> the fix is not adopted as is.
