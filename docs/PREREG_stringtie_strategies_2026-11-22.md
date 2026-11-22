# Row 1893: what StringTie's coverage trimming and its coverage-driven assembly contribute, by its own switches (pre-registered 2026-11-22)

User (2026-11-22): "can we try any of stringtie's strategies for coverage-based trimming or coverage-based extension?"
Record: the whole StringTie long-read flow was ported and compared in row 1658 (fewer correct chains at lower precision than
Rustle's default on gorilla and A119b; not kept), and the reads' own StringTie transfrag chains in row 1639. Both were judged
against ape / human RefSeq only, never on complete truth, and neither separated the two strategies named here. Rows 1882-1892:
Rustle's own extension (junction-pair links) pays only against the non-human ape RefSeq sets and is wrong on complete truth.

Before any port this row lets StringTie 3.0.3 isolate each strategy with its own switches:
- **Trimming:** `-L` against `-L -t` (`-t` = no trimming of predicted transcripts based on coverage).
- **Coverage-driven assembly / extension:** `-L` against `-L -R` (`-R` = the long reads cleaned and collapsed, not assembled).
  What `-L` reports and `-R` does not is what the splice graph and its coverage flow BUILD beyond the collapsed reads.
Nothing in Rustle changes. Workspace `r1893_stringtie_strategies`.

## Data
Nine datasets: fly, mouse, zebrafish simulated (complete truth); mouse ES (mouse contigs), mouse heart; WTC11 r2 (human iPS);
chimpanzee testis, Sumatran orangutan testis, KB3781 gorilla fibroblast. BAMs, truths and the existing `-L` outputs as in row
1880 (`r1880_ranked_vs_rivals/sets.tsv`); `-L -t` and `-L -R` are run here with `-p 8`.

## Measures
- Chains (multi-exon intron chains, `curve.py load`): chains / correct / precision of the three outputs.
- **BUILT** = chains of `-L` that `-R` lacks; **GONE** = chains of `-R` that `-L` lacks; each with its correct count.
- Trimming: chains that differ between `-L` and `-L -t`; among correct-chain transcripts the share with the 5' end and the 3'
  end within 50 bp of the truth (`r1831_mono_end_rules/score.py`), and exonic base precision, for both.
- Rustle beside them: today's `final.gtf` (row 1880 runs) on the same measures.

## Reading, fixed now
- **Extension is WORTH A PORT** when BUILT chains are at least half correct on at least 2 of the 3 complete-truth sets AND are
  at least as often correct as StringTie's whole `-L` output on at least 2 of the 3 real libraries with curated or mouse
  annotation (mouse ES, heart, WTC11). Then the next row specifies the port (which part of the flow, judged on which data).
- **APE-ANNOTATION EFFECT** when BUILT chains are under a quarter correct on at least 2 of the 3 complete-truth sets and above
  StringTie's own precision on at least 2 of the 3 ape libraries: the same pattern as Rustle's extension (row 1883); no port.
- Anything else: reported as it is, no port proposed from this row.
- **Trimming is WORTH A PORT** when `-L` beats `-L -t` by at least 2 points of ends within 50 bp (5' or 3') on at least 6 of 9
  datasets, at least 4 of the 5 non-ape ones, without losing more than 0.1% of the correct chains. Otherwise no port.
- A port is a build and a default question, both the user's call.

## Limits, stated in advance
- `-R` is StringTie's collapsed reads with a one-read floor: its chain count and precision are those of raw reads; only the
  set difference with `-L` is read.
- BUILT can include chains that exist in a read `-R` dropped while cleaning; the row cannot tell those from built ones.
- The switches isolate StringTie's strategies inside StringTie's pipeline (its junction correction, its 1.5 coverage floor);
  a port into Rustle would meet different inputs.
- Annotated ends are weak truth on real libraries (rows 1487 / 1831); the simulated truths are exact at both ends.

## Result (2026-11-22): neither strategy is worth a port; in long-read mode StringTie selects among read chains, it builds almost none

Job 4494858, eighteen runs exit 0; `r1893_stringtie_strategies/strat.NAME.txt` (`strat.py`). Chains / correct / precision.
| dataset | StringTie `-L` | `-L -R` (collapsed reads) | BUILT (in `-L`, not in `-R`) | GONE (in `-R`, not in `-L`) | Rustle `final.gtf` |
|---|---|---|---|---|---|
| fly simulated | 6,801 / 6,549 / 0.963 | 8,248 / 6,843 / 0.830 | 2 / 0 | 1,449 / 294 / 0.203 | 6,940 / 6,753 / 0.973 |
| mouse simulated | 24,730 / 23,705 / 0.959 | 38,157 / 31,280 / 0.820 | 69 / 2 / 0.029 | 13,496 / 7,577 / 0.561 | 25,864 / 25,077 / 0.970 |
| zebrafish simulated | 13,270 / 11,826 / 0.891 | 18,219 / 12,877 / 0.707 | 40 / 0 | 4,989 / 1,051 / 0.211 | 13,300 / 12,129 / 0.912 |
| mouse ES | 44,591 / 16,401 / 0.368 | 224,420 / 23,367 / 0.104 | 398 / 31 / 0.078 | 180,227 / 6,997 / 0.039 | 41,641 / 17,073 / 0.410 |
| mouse heart | 22,478 / 12,577 / 0.560 | 52,800 / 14,897 / 0.282 | 109 / 10 / 0.092 | 30,431 / 2,330 / 0.077 | 21,160 / 12,922 / 0.611 |
| WTC11 r2 | 36,963 / 20,064 / 0.543 | 101,098 / 22,261 / 0.220 | 367 / 71 / 0.193 | 64,502 / 2,268 / 0.035 | 31,336 / 20,016 / 0.639 |
| chimpanzee testis | 61,509 / 21,532 / 0.350 | 178,894 / 24,596 / 0.138 | 2,148 / 406 / 0.189 | 119,533 / 3,470 / 0.029 | 46,969 / 20,309 / 0.432 |
| Sumatran testis | 72,495 / 16,981 / 0.234 | 193,721 / 18,760 / 0.097 | 3,244 / 388 / 0.120 | 124,470 / 2,167 / 0.017 | 54,748 / 15,552 / 0.284 |
| KB3781 gorilla fibroblast | 58,164 / 22,087 / 0.380 | 251,259 / 26,454 / 0.105 | 1,241 / 181 / 0.146 | 194,336 / 4,548 / 0.023 | 53,153 / 23,407 / 0.440 |

**Coverage trimming (`-t`): no port.** With `-L`, switching trimming off changes nothing at all on nine of nine datasets: the
same chains, and the same share of ends within 50 bp to three decimals. In long-read mode the transcript ends do not come
from coverage trimming. (Ends within 50 bp, 5' / 3', StringTie against Rustle: simulated 0.90-0.99 / 0.94-0.97 against 0.99-1.00
/ 1.00; real libraries 0.48-0.61 / 0.60-0.69 against 0.50-0.66 / 0.56-0.69.)

**Coverage-driven assembly: no port.**
- What the graph and flow BUILD beyond the collapsed reads is small and mostly wrong: 2, 69, 40 chains on the simulated sets
  with 0, 2, 0 true; 8-9% annotated on the mouse libraries, 19% on WTC11, 12-19% on the three ape libraries. That is under
  StringTie's own precision on every one of the nine datasets, so it is not WORTH A PORT (complete truth 0 of 3) and also not
  the APE-ANNOTATION EFFECT as defined (BUILT is below StringTie's precision on the apes too). Registered verdict: "anything
  else", no port.
- **What StringTie does in long-read mode is select.** Its output is its collapsed read chains minus the ones its coverage
  and fraction floors drop: 97-99.9% of its chains are read chains. On Sumatran testis the collapsed reads hold 18,760
  annotated chains in 193,721; StringTie keeps 16,593 of them in 69,251. The BUILT chains are 388 of its 16,981 correct ones
  there (373 not in Rustle), about half of the 686-chain gap of row 1880; the other half is which one-read chains it keeps.
- So the premise in rows 1883 / 1889 that StringTie's lead on short-read testis libraries comes from chains longer than its
  reads is only half right. I wrote that reading; this table corrects it.

**What follows (not tested here).** The StringTie strategy that matters for the lost comparison is its SELECTION of low-count
read chains (coverage per base at least 1, isoform fraction 0.01 of the locus), not trimming and not extension. Rustle's
counterpart is the order of its one-read tier. Row 1890 tried one order (the model's `cov`) and it was worse everywhere; a
locus-share order in StringTie's sense has not been tried. It would be a new registered row, judged on all twenty datasets,
and it would re-open a question row 1892 called closed: the user's call.
