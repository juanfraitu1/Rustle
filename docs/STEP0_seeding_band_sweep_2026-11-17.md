# Step 0, row 1767: what do secondary alignments beyond the 0.98 band add to the assembly? (2026-11-16 / 17)

User question (2026-11-16): the assembler is already good on primary alignments and the near-tied secondaries add a little;
is there a wider band of secondaries that still gives good seeding or splicing information for the loci and isoforms that
are still missing? Descriptive, no bar set before the run, nothing here changes a default. Not for `advisor_share` until the
user says so.

Runs: `assemble` with today's defaults and the driver snapshot of the assembler artifact, seeding pool = primary only, then
`--seed-pool good --seed-as-ratio R` for R = 1 (exact ties), 0.99, 0.98 (the default), 0.97, 0.95, 0.9, 0.8, 0.65, 0.5, then
`--seed-pool all`. Samples: gorilla OR6737 (`GGO_mm.bam`), gorilla KB3781 (`GGO_KB3781_fib_mm.bam`), human A119b. Multi-exon exact intron chains against RefSeq. Workspace `bandsweep_2026-11-16` (= `r1767_bandsweep`),
outside the repo: `run.sbatch` (one job per sample, about 100 s per band on gorilla), `band.py`, `surv.py`, `complete.py`,
`sec_only_band.py` (row 1671's `bench/famcmp/sec_only_loci.py` with the band as `SEC_ONLY_RATIO`), `hap_join.py`, `hap.sbatch`.

## 1. Each band against the primary pool
Every band's assembly is the primary-pool assembly plus additions; losses are small (OR6737: 93 chains, 3 correct, at 0.98;
477 / 24 with every secondary).

OR6737, primary pool 75,018 chains / 26,263 correct / 0.350:
| band | chains | correct | precision | gained / correct | this step alone: new / correct / precision | new introns / in RefSeq | genes with a first exact chain |
|---|---|---|---|---|---|---|---|
| 1 | 75,284 | 26,334 | 0.350 | 293 / 72 | 293 / 72 / 0.246 | 483 / 0.63 | 53 |
| 0.99 | 75,735 | 26,392 | 0.349 | 777 / 131 | 496 / 61 / 0.123 | 1,418 / 0.62 | 81 |
| 0.98 | 76,139 | 26,424 | 0.347 | 1,214 / 164 | 468 / 33 / 0.071 | 2,178 / 0.59 | 94 |
| 0.97 | 76,486 | 26,448 | 0.346 | 1,588 / 189 | 420 / 26 / 0.062 | 2,740 / 0.56 | 107 |
| 0.95 | 77,275 | 26,480 | 0.343 | 2,424 / 224 | 889 / 34 / 0.038 | 3,978 / 0.50 | 123 |
| 0.9 | 78,978 | 26,562 | 0.336 | 4,176 / 307 | 1,889 / 85 / 0.045 | 6,509 / 0.47 | 180 |
| 0.8 | 80,978 | 26,618 | 0.329 | 6,235 / 365 | 2,220 / 62 / 0.028 | 9,018 / 0.42 | 220 |
| 0.65 | 83,034 | 26,680 | 0.321 | 8,341 / 432 | 2,386 / 68 / 0.029 | 11,562 / 0.39 | 254 |
| 0.5 | 85,164 | 26,724 | 0.314 | 10,539 / 477 | 2,547 / 49 / 0.019 | 14,089 / 0.36 | 277 |
| all | 98,220 | 26,882 | 0.274 | 23,679 / 643 | 13,769 / 168 / 0.012 | 28,519 / 0.24 | 382 |

KB3781, primary pool 73,652 / 26,240 / 0.356:
| band | gained / correct | this step alone: new / correct / precision | new introns / in RefSeq | genes with a first exact chain |
|---|---|---|---|---|
| 1 | 211 / 13 | 211 / 13 / 0.062 | 242 / 0.23 | 7 |
| 0.99 | 426 / 47 | 217 / 34 / 0.157 | 595 / 0.49 | 27 |
| 0.98 | 695 / 77 | 282 / 32 / 0.114 | 1,063 / 0.50 | 38 |
| 0.97 | 973 / 91 | 297 / 15 / 0.051 | 1,385 / 0.51 | 42 |
| 0.95 | 1,609 / 131 | 686 / 40 / 0.058 | 2,200 / 0.51 | 57 |
| 0.9 | 2,805 / 191 | 1,323 / 60 / 0.045 | 3,651 / 0.47 | 83 |
| 0.8 | 4,677 / 251 | 2,031 / 63 / 0.031 | 6,280 / 0.42 | 121 |
| 0.65 | 7,301 / 327 | 2,882 / 78 / 0.027 | 9,610 / 0.37 | 177 |
| 0.5 | 9,843 / 360 | 2,759 / 36 / 0.013 | 12,389 / 0.32 | 191 |
| all | 22,041 / 482 | 12,602 / 127 / 0.010 | 26,021 / 0.21 | 265 |

Reading: correct chains and genes reached rise with every band and precision falls steadily, with no cliff down to 0.9. The
chains a step adds are RefSeq-correct 1 to 16% of the time past exact ties; about half of the new introns are RefSeq introns
down to 0.9. The isoforms other tools find and the primary pool lacks are not here: of 4,714 such RefSeq chains on OR6737
(StringTie, FLAIR, IsoQuant, bambu) the 0.98 band recovers 22 and every secondary 61; KB3781 7 and 21 of 3,897.
Limit: RefSeq counts a real unannotated copy as wrong, and a RefSeq-correct chain at a paralog can be a shadow (section 4).

## 2. Does the chain score filter keep the good additions? (`surv.py`, `PREFIX.final.gtf` of each band)
On the primary pool the filter keeps 0.949 of correct and 0.728 of wrong chains (OR6737), 0.927 and 0.772 (KB3781). On the
chains a band adds:
| sample | band | added / correct | share of correct kept | share of wrong kept | precision of the kept additions |
|---|---|---|---|---|---|
| OR6737 | 1 | 293 / 72 | 0.40 | 0.31 | 0.296 |
| OR6737 | 0.98 | 1,214 / 164 | 0.41 | 0.38 | 0.143 |
| OR6737 | 0.9 | 4,176 / 307 | 0.44 | 0.33 | 0.097 |
| OR6737 | all | 23,679 / 643 | 0.83 | 0.49 | 0.045 |
| KB3781 | 1 | 211 / 13 | 0.54 | 0.28 | 0.111 |
| KB3781 | 0.98 | 695 / 77 | 0.49 | 0.33 | 0.157 |
| KB3781 | 0.9 | 2,805 / 191 | 0.55 | 0.35 | 0.103 |
| KB3781 | all | 22,041 / 482 | 0.80 | 0.51 | 0.034 |
The filter halves the additions without separating right from wrong. `final.gtf` as a whole still gains with width (OR6737
24,930 @0.413 primary pool, 25,200 @0.399 at 0.9; KB3781 24,327 @0.399 and 24,502 @0.391).

## 3. Are loci the primary pool already finds made more complete? (`complete.py`)
Found = RefSeq gene with >= 1 of its introns in the primary-pool assembly (15,205 OR6737; 13,351 KB3781).
| sample | measure | primary | 0.98 | 0.9 | all |
|---|---|---|---|---|---|
| OR6737 | mean share of a found gene's RefSeq introns covered | 0.8278 | 0.8289 | 0.8303 | 0.8326 |
| OR6737 | found genes with every RefSeq intron covered | 5,884 | 5,895 | 5,918 | 5,951 |
| OR6737 | found genes gaining an intron | | 82 | 172 | 327 |
| OR6737 | RefSeq chains with only a truncated model, completed (of 12,587) | | 10 | 16 | 38 |
| KB3781 | mean share of a found gene's RefSeq introns covered | 0.8515 | 0.8520 | 0.8529 | 0.8546 |
| KB3781 | found genes with every RefSeq intron covered | 5,594 | 5,602 | 5,615 | 5,632 |
| KB3781 | found genes gaining an intron | | 46 | 104 | 220 |
| KB3781 | RefSeq chains with only a truncated model, completed (of 9,291) | | 8 | 14 | 20 |
Where the gained correct chains sit at 0.98: OR6737 99 in genes the primary pool does not find, 17 first exact chain of a
found gene, 48 a further isoform of a gene with one; KB3781 37 / 5 / 35. Added chains inside found genes: OR6737 670 with
65 correct (0.097), KB3781 292 / 40 (0.137); about 0.05 at wide bands. Secondaries do not complete found loci.

## 4. Are the new loci real? KB3781 against its own haplotypes (`sec_only_band.py`, slack 2; method of row 1671)
Population per band: loci of that band's assembly with no primary alignment on an exon and >= 1 molecule whose secondary
is inside the band. Up to 30 molecules per locus aligned to both haplotype assemblies. L = molecules match this locus's
haplotype copy and not their primary's; P = the reverse (shadow); both; other; mixed; thin = fewer than 3 judged.
| band | loci in the assembly | only on secondaries | L | both | P | other | mixed | thin |
|---|---|---|---|---|---|---|---|---|
| 0.98 | 14,801 | 201 | 2 | 24 | 143 | 10 | 1 | 21 |
| 0.95 | 15,040 | 368 | 1 | 20 | 277 | 54 | 0 | 16 |
| 0.9 | 15,503 | 727 | 3 | 34 | 511 | 151 | 5 | 23 |
| 0.8 | 16,384 | 1,454 | 3 | 72 | 1,031 | 305 | 7 | 36 |
| 0.5 | 18,768 | 3,512 | 4 | 174 | 2,574 | 675 | 2 | 83 |
The 0.98 row reproduces row 1671 (194 loci, about 80% shadows) on today's defaults: 71% P. Widening multiplies the loci and
the share that is real does not grow: L stays at 1 to 4 loci, P stays at 70 to 75%, and `other` (molecules whose perfect
hit is on neither copy) grows from 5% to about 20%.

Added multi-exon chains that lie in a judged locus (`hap_join.py`): 127 of 695 at 0.98, 804 of 2,805 at 0.9, 4,072 of 9,843
at 0.5; the rest are in loci a primary alignment also touches and are not judged here. RefSeq-correct added chains by the
class of their locus:
| band | correct added chains in judged loci | in L | in both | in P | in other | in thin |
|---|---|---|---|---|---|---|
| 0.98 | 19 | 0 | 7 | 6 | 0 | 6 |
| 0.95 | 36 | 0 | 6 | 25 | 1 | 4 |
| 0.9 | 63 | 1 | 7 | 49 | 1 | 5 |
| 0.8 | 97 | 1 | 8 | 82 | 2 | 4 |
| 0.5 | 161 | 0 | 9 | 116 | 32 | 4 |
So RefSeq does not undercount these additions; it overcounts them. Past 0.98 most RefSeq-correct chains on new loci are
correct structures built from molecules of another copy.

Strict re-read (slack 0, same loci): no locus is L below the 0.5 band, and the correct chains move further into P.
| band | only on secondaries | L | both | P | other | mixed | thin | correct added chains in judged loci: L / both / P |
|---|---|---|---|---|---|---|---|---|
| 0.98 | 201 | 0 | 8 | 159 | 13 | 0 | 21 | 0 / 5 / 8 |
| 0.95 | 368 | 0 | 10 | 285 | 57 | 0 | 16 | 0 / 5 / 25 |
| 0.9 | 727 | 0 | 19 | 520 | 163 | 2 | 23 | 0 / 6 / 49 |
| 0.8 | 1,454 | 0 | 45 | 1,044 | 325 | 4 | 36 | 0 / 7 / 82 |
| 0.5 | 3,512 | 2 | 129 | 2,614 | 682 | 2 | 83 | 0 / 2 / 122 |

## 5. Human A119b (sections 1 only; `A119b/band.tsv`)
Primary pool 224,494 chains / 39,823 correct / 0.177; 10,586 RefSeq chains it lacks are found by another tool.
| band | chains | correct | precision | gained / correct | this step alone: new / correct / precision | new introns / in RefSeq | other tools' chains recovered | genes with a first exact chain |
|---|---|---|---|---|---|---|---|---|
| 1 | 226,619 | 39,904 | 0.176 | 2,444 / 90 | 2,444 / 90 / 0.037 | 2,047 / 0.20 | 24 | 38 |
| 0.99 | 231,103 | 40,024 | 0.173 | 7,602 / 218 | 5,437 / 129 / 0.024 | 6,373 / 0.16 | 40 | 105 |
| 0.98 | 236,301 | 40,086 | 0.170 | 13,044 / 286 | 5,957 / 72 / 0.012 | 11,332 / 0.11 | 47 | 140 |
| 0.97 | 241,044 | 40,120 | 0.166 | 17,986 / 324 | 5,614 / 44 / 0.008 | 15,956 / 0.09 | 53 | 154 |
| 0.95 | 247,680 | 40,158 | 0.162 | 24,929 / 366 | 7,823 / 44 / 0.006 | 22,349 / 0.08 | 63 | 172 |
| 0.9 | 256,321 | 40,214 | 0.157 | 34,254 / 440 | 10,779 / 83 / 0.008 | 30,344 / 0.07 | 79 | 208 |
| 0.8 | 268,252 | 40,301 | 0.150 | 46,852 / 537 | 15,131 / 109 / 0.007 | 42,925 / 0.06 | 90 | 277 |
| 0.65 | 285,024 | 40,330 | 0.142 | 64,186 / 596 | 18,753 / 88 / 0.005 | 60,120 / 0.06 | 103 | 290 |
| 0.5 | 303,310 | 40,364 | 0.133 | 82,882 / 655 | 20,036 / 67 / 0.003 | 78,905 / 0.05 | 112 | 318 |
| all | 436,392 | 40,529 | 0.093 | 216,968 / 870 | 136,560 / 232 / 0.002 | 217,717 / 0.02 | 133 | 441 |
The human sample is harsher than the gorillas on every column: a step's additions are 0.2 to 4% correct, and only 5 to 20%
of the new introns are RefSeq introns (about 50% on the gorillas down to 0.9). The gorilla reading "half of the new introns
are annotated" does not transfer. The filter, completeness and haplotype sections were not run on A119b (no haplotype truth).

## Summary
- A wider band is a net gain in RefSeq-correct chains and genes reached and a steady loss of precision; no band adds chains at
  the primary pool's precision, and the chain score filter does not select the good additions.
- Loci the primary pool finds are not made more complete by secondaries.
- On KB3781 the new loci of every band are mostly shadows by the haplotype truth, including most of their RefSeq-correct
  chains. The gain in section 1 is therefore largely not a gain in loci that are expressed.
- Not tested: a wider band as junction support only (no new seeds). Half of the new introns are RefSeq introns down to 0.9
  on both gorilla samples (5 to 20% on A119b), but section 4 says an annotated structure is not evidence of origin, so that test needs a truth of
  origin (KB3781), not RefSeq.
