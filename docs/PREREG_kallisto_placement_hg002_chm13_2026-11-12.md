# PREREG lr-kallisto PLACEMENT AWAY FROM THE REFERENCE, WITH TRUTH: HG002 reads on CHM13 (register 1651, 2026-11-12)

Why: r1647's exact-k-mer rule was 1-2% wrong where reads and reference are one individual; r1650 (OR6737) showed reach
and O2 agreement hold on another individual but could not measure the wrong rate, and 13.5% of its picks left the
primary's locus. This row measures the wrong rate away from the reference.

Sample: HG002 Kinnex aligned to CHM13 (`HG002_kinnex.t2t.bam`, existing). Catalog and transcripts: the A119b base run on
CHM13 (`~/runs/human_A119b/a1.fam.copies.*`, `a1.gtf`): loci and sequences are CHM13's; the catalog was built from
another sample's reads, which only limits which families exist. Truth: HG002's own haplotypes (r1644 `mat.mmi`,
`pat.mmi`), `haplotype_truth.py` + `strict_truth.py` slack 0, OUTSIDE rows dumped.
Truth caveat, stated first: copy loci are CHM13 copy sequences placed on HG002's haplotypes; where paralogs differ less
than the two individuals do, a copy can be placed on its paralog's locus. The primary is scored by the same truth, so
the paired comparison is fair; absolute wrong rates carry this noise for both.

Set: extended `as_table` rebuilt (the stored table has 15 columns); `panel.py` top 80, 2 Mb cost guard, first 40;
`copy_assign --tie-ms 10` (binary `target_r17`); contested rows with truth, per molecule.
Arm: r1647 unchanged (reference = `multi_iso_ref.py ... other` on `a1.*`; kallisto 0.52.0 `-k 63`,
`bus --long --threshold 0.8 -n`; one-locus rule). Control: minimap2 `map-hifi` on the same reference (no bar).
Scorer: pick right = span overlaps a truth copy; primary right = its record overlaps a truth copy.

Bars:
0. Power: kallisto assigns >= 100 decidable molecules, else UNDERPOWERED.
1. Wrong rate on its assigned decidable molecules <= 5%.
2. Tool right / primary wrong >= 2 x (tool wrong / primary right), first count >= 10.
3. Wall <= 3 h (limit = FAIL).
Readout, no bar: reach; negatives assigned and how many of those picks hold the primary; right / wrong of the picks
that leave the primary's locus (the r1650 risk); kallisto vs O2 on and off the shared set; the control.

Expected: no confident prediction. The HG002 maternal split (16 of 16 right, reach 0.11) says accuracy holds and reach
falls; r1650's moved reads say the wrong rate may rise.
Verdicts: 0-2 pass = PASS: the rule is a candidate certificate beside O2 (opt-in design row next, no default change).
1 or 2 fails = FAIL: the rule is reference-individual only.

Outputs: `/anvme/workspace/v138da13-rustle/r1651_kallisto_hg13/`.

## RESULT r1651 (2026-11-12, job 4467055, 51 min; `r1651_kallisto_hg13/{result,control,strict}.txt`)
Panel 40 families, 37 regions, 10.2 Mb; 52,511 contested rows, O2 assigns 2,289; truth 277,686 hits; 4,201 decidable
molecules, 4,152 negative. Reference 259,778 sequences. kallisto pseudoaligned 1,470 of 8,353 reads (0.176; 0.70 on the
same individual's own haplotype in r1647).
- Bar 0 PASS: 181 decidable molecules assigned (reach 0.043; r1647 0.28 and 0.60).
- Bar 1 FAIL: 166 right / 15 wrong = 8.3% (bar 5%). The primary on the same molecules: 163 / 18.
- Bar 2 FAIL: tool right / primary wrong 3, reverse 0 (first count under 10).
- Bar 3 PASS.
**FAIL as registered.** Away from the reference individual the exact-k-mer rule keeps 4% reach and its picks are the
primary's: 175 of 181 sit on the primary's locus (161 right / 14 wrong), 6 move (5 right / 1 wrong). The r1647
advantage over the primary does not transfer; r1650's reach on OR6737 did not predict this (gorilla to gorilla, 0.38).
Readout: negatives assigned 535 of 4,152 (0.129), 439 of them on the primary's locus. minimap2 control on the same
reference: 2,415 / 517 vs the primary 2,434 / 498 (4 vs 23). Where both kallisto and O2 assign (46): 44 right each.
Side finding, same truth and caveat: O2 (`--tie-ms 10`) on this cross-individual set is 170 right / 175 wrong on its
decidable assigned rows (strict_truth, leftmost-position primary 189 / 156); where kallisto abstains O2 is 121 right /
159 wrong. Either O2's certificate is unreliable when the copy sequences are another individual's, or the truth
(CHM13 copies placed on HG002 haplotypes) mislabels paralogs; MCL601 is 0 right / 9 wrong for kallisto AND the primary.
Not separated here; needs its own row before any cross-individual O2 figure is quoted.

## O2 SIDE FINDING, REGISTERED AS AN O2 ERROR (user decision 2026-11-12; `r1651_kallisto_hg13/o2_breakdown.txt`)
The user registers the cross-individual result as O2 being wrong, to drive improvements. The alternative (the truth
mislabels paralogs when CHM13 copies are placed on HG002 haplotypes) is not excluded by any run and stays noted here.
O2 `--tie-ms 10`, HG002 reads on CHM13, 345 decidable assigned rows, primary judged by record overlap:
| O2 | primary | O2's copy vs the primary's locus | rows |
|---|---|---|---|
| right | right | kept | 168 |
| right | wrong | moved | 2 |
| wrong | right | kept | 5 |
| wrong | right | moved | 65 |
| wrong | wrong | kept | 49 |
| wrong | wrong | moved | 56 |
- When O2 moves a read off the primary's locus it is right 2 times and wrong 121. When it keeps the locus: 168 / 54.
- 171 of the 175 wrong rows are three families: MCL189 (123 wrong, 0 right), MCL601 (27, 0), MCL223 (21 wrong, 7 right).
  The other nine families: 138 right / 4 wrong.
Pointers for O2 (each needs its own registered row): (1) a cross-individual guard on moves off the primary's locus;
(2) a per-family check that flags families where the sample's reads fit no catalog copy cleanly (the three above);
(3) the same breakdown on KB3781 and HG002-on-pat, where reads and reference are one individual, to see whether moves
are safe there (r1514: O2 right where the primary is wrong 5 vs 1).
