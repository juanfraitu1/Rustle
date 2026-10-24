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

## Result (2026-10-12) — B2 FAILED: not adopted as default; kept OPT-IN (`RUSTLE_PSV_BLOCKS=1`)

Run: binary with the mask vs a37399a5 behaviour, same 56-family subset and arguments (`/mnt/linuxdisk/tmp/psv_sun/
{run_fix2,run_blk}/`, scorer `bench/psv_sun/block_eval.py` -> `run_blk/block_eval.json`). PSV columns 33,015 -> 4,377.
- **B1 PASS**: assignments, quant byte-identical (copy_conversions: 5 -> 0 calls; 4 are MCL58, where each of its 4 copies
  is "converted" from the others — circular; not tested further).
- **B2 FAIL**: `aligned` pairs retained **0.257** (bar 0.99).
- **B3 PASS**: non-homologous pairs removed **0.9987** (both_unrelated 9,807/9,810, terminal_overhang 5,380/5,399,
  no_alignment_unrelated 2,386/2,386).
- **B4**: both_homologous retained 0.068, no_alignment_homologous 0.014.
- Library tests 1,023 pass (3 new: unrelated middle exon / unrelated terminal overhang -> no columns; 80%-identity
  homolog keeps all 80 columns).

**Why B2 failed (post hoc, not a re-scoring of the bar).** (1) The pair yardstick counts every copy with an allele at a
decisive column, including copies equal to copy 0, whose column disappears when the only differing copy was junk.
Restricted to pairs that DIFFER from copy 0: `aligned` retained 0.488. (2) The `aligned` label is about the STRETCH
(copy 0's offset inside a minimap2 block), not about O2's pairing: by the independent position verifier, 9,703 of
16,183 such differing pairs were the global DP's own MISALIGNMENTS (O2's column != minimap2's position) — retained
0.212 — while the 5,697 position-verified ones were retained **0.956**. So the mask mostly removes misalignment and
junk, but it also loses 4.4% of correct columns and cannot RECOVER the misaligned ones.

**Root cause exposed:** `banded_msa_pair` costs (mismatch 1, gap open 32, extend 1) make a ~30-base misalignment
cheaper than one gap, so inside homologous regions 60% of O2's apparent copy differences are misaligned columns.
The next step is the DP's scoring (match reward, log-odds mismatch, a cheap long-gap piece, free end gaps), with the
mask as a final guard — pre-registered with a position-verified yardstick on families held out from this subset.
Default behaviour: byte-identical to a37399a5 on all 7 outputs (`run_optin/`).
