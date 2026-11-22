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
