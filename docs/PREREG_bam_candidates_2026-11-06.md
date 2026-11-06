# PREREG: read-star only on BAM-reached copies (row 1557, 2026-11-07 00:10)

## Finding that motivates it (MCL1, 15,828 reads with primary in the block; gw_off baseline)
The nine chr13 hotspot families evaluate the same 30,125 molecules (212,723 read x family rows; most reads in
6-9 of the 9 families) while their copies share exons in only 14 copy pairs: at the read level they are one
super-family (a families-stage observation, recorded here, independent of this lever). Per molecule the BAM
holds a median 41 records reaching a median 8 (mean 8.9) of MCL1's 31 copies; all 48 assigned MCL1 reads
(13 with primary outside the block) have their copy reached. Read-star aligns every molecule to all 31.
Related: row 1514 (--primary-on-family, primary-location filter on ties, FAIL as registered) is a different
rule: it drops molecules; this restricts copies per molecule.

## Rule
A molecule's read-star candidates = the copies any of its BAM records (primary, secondary, supplementary)
overlaps; it is aligned to those only. Arm 3 (CLI, this row): per copy, `minimap2 -c --eqx -x splice
--secondary=yes -p 0.3 -N 31 -t 16` of that copy's genomic window against the reads reaching it; time the whole
loop; concatenate PAFs; compare with ctl_genomic.paf (469 s, all 31 copies).
## Bars (fixed before any number)
- Time: loop wall <= 160 s (>= 3x) at 16 threads.
- Observations: for (read, copy) pairs present in both PAFs, AS identical in >= 99% (the -p 0.3 pruning relative
  to the best over all targets is the known difference).
- Assigned: all 35 assigned reads (in the PAF) keep their copy as the unique best-AS copy among reached copies.
- Guard (the registered failure mode): of the 41 tied reads, the number whose best-AS copy becomes unique by
  >= 10 AS once unreached copies are removed (a tie that resolves only because its co-leader was pruned) <= 5.
## Prediction
Time passes (single-target chaining has no cross-copy seed hits). Observations pass. Assigned pass.
Guard FAILS: 28 of 35 assigned reads and 39 of 41 tied reads have a competitor within 10 AS that the BAM does
not reach, so ~30-39 ties are predicted to resolve falsely. If so the lever is closed at the CLI stage and the
Rust hook is not written. A pass would send it to the hot9 panel in Rust with the same guard.
