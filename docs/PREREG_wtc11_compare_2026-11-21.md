# Row 1862: a blind comparison on WTC11, then WTC11 as development data for new rules (pre-registered 2026-11-21)

User (2026-11-21): run the comparison on the new library and use it to derive more rules for precision and recall. WTC11
replicate 2 (ENCODE ENCSR507JOF / ENCFF370NFS, row 1860) has been used for nothing but the polyA question: Rustle's assembly
at today's defaults exists (20,016 correct of 31,336, 0.639, seen), StringTie and IsoQuant have not been run on it.
Workspace `r1862_wtc11_compare`. Committed before either tool ran.

## Part 1, blind: predictions for the comparison (the page's four statements)
StringTie (`stringtie -L -p 8`) and IsoQuant 4.1.0 (`-d pacbio_ccs`, defaults), de novo, the commands of `apetools/job.sbatch`,
on `reads/HSA/WTC11_r2_mm.bam`; exact intron chains against RefSeq on CHM13 (`ptr_ppa_rules/gffcmp/hsa.ref.gtf`).
- P1: Rustle has more correct chains than StringTie and than IsoQuant.
- P2: Rustle's precision is not below StringTie's.
- P3: Rustle cut by `reads` to IsoQuant's correct count is at least as precise as IsoQuant minus 0.02.
- P4: no tool has both more correct chains and higher precision than Rustle.
Reported as predicted / not, each. This part is the only blind use of the library; after it, replicate 2 is development data.

## Part 2, development (descriptive, fixed tables; nothing built in this row)
On replicate 2, with the scripts of rows 1851 and 1853 unchanged (`why.py`'s classes and agreement, `iq.py`'s sets and intron
bins): (a) the classes of Rustle's wrong chains and how many another tool reports; (b) Rustle-only, IsoQuant-only and
StringTie-only sets; (c) the correct chains the other tools have and Rustle lacks: absent before the filter, or removed by
it; by intron count and by `reads` of the nearest Rustle chain where there is one.
A pattern is taken forward as a RULE CANDIDATE only when it is at least 5% of Rustle's wrong chains (precision) or at least
2% of the other tools' correct chains that Rustle lacks (recall) AND shows on mouse ES as well (rows 1851 / 1853 tables).

## Check data, set aside now
Replicates 1 and 3 of the same experiment (ENCFF563QZR, ENCFF245IPA) are downloaded and mapped with the same command and
NOT assembled or looked at in this row. Every rule derived in part 2 is its own registered row, built as an opt-in switch,
and judged on replicates 1 and 3 and on the twelve datasets. Limit, stated now: the replicates share the cell line and the
protocol with replicate 2; they guard against fitting one run, not against fitting the cell line.

## Expectation
P1-P4 hold (they hold on every dataset of the page). Part 2: the same picture as mouse ES, a Rustle-only low-count tail and
no class of its own; I do not expect a new rule from the precision side. The recall side (correct chains only the other
tools have, 4-8% of theirs on the other libraries) has not been examined before and is where a rule could come from.

## Result (2026-11-21)

Job 4491326 (StringTie exit 0, 151 s; IsoQuant exit 0, 205 s). `r1862_wtc11_compare/{cmp.wtc.txt,why.wtc.txt,cmp.es.txt}`
(`compare.py` of row 1838 and `why.py` of row 1851, unchanged apart from paths).

**Part 1, blind (175,219 truth chains):**
| tool | reported | correct | precision | recall |
|---|---|---|---|---|
| Rustle filtered | 31,336 | 20,016 | 0.639 | 0.1142 |
| Rustle unfiltered | 33,769 | 20,076 | 0.595 | 0.1146 |
| StringTie | 36,963 | 20,064 | 0.543 | 0.1145 |
| IsoQuant | 26,191 | 16,609 | 0.634 | 0.0948 |
- P1 (more correct chains than both): **NOT as predicted against StringTie** (20,016 against 20,064, 48 fewer); as predicted
  against IsoQuant (+3,407).
- P2 (precision not below StringTie's): as predicted (0.639 against 0.543).
- P3 (at IsoQuant's correct count within 0.02): as predicted, by a wide margin (reads >= 3: 16,823 / 23,110 = 0.728 against 0.634).
- P4 (not dominated): as predicted. Rustle is the most precise of the three as reported, on a library it had never seen.

**Part 2, development.**
(a) Wrong chains (11,320): 79% are also reported by another tool; all classes 100% canonical. Novel site beside an annotated
one 45.8%, piece of an annotated chain 17.1% (mouse ES 6.4%), combination 15.2%, no annotated intron 12.4%, shifted site
3.6% (28% shared: Rustle's own, as on ES). Precision on chains another tool shares: 0.672 / 0.677 / 0.674 for the three tools.
(b) Tool-only sets: Rustle 3,890 at 0.401; StringTie 10,578 at 0.209; IsoQuant 1,918 at 0.124. Chains all three report:
15,001 / 20,324 = 0.738. Rustle's own tail is the best of the three here (mouse ES: 0.135).
(c) **Recall side: true chains another tool has and Rustle's `final.gtf` lacks**, by where Rustle lost them (row 1838's classes):
| | total | removed by the filter | piece of a Rustle chain | Rustle has only a shorter one | near (<= 10 bp) | other form | an intron Rustle has nowhere |
|---|---|---|---|---|---|---|---|
| WTC11, in StringTie | 2,594 | 21 | 168 (6.5%) | 611 (23.6%) | 37 | 214 (8.2%) | 1,543 (59.5%) |
| WTC11, in IsoQuant | 621 | 24 | 192 (30.9%) | 116 (18.7%) | 29 | 41 | 219 (35.3%) |
| mouse ES, in StringTie | 1,140 | 54 | 278 (24.4%) | 152 (13.3%) | 15 | 108 (9.5%) | 533 (46.8%) |
| mouse ES, in IsoQuant | 795 | 68 | 502 (63.1%) | 48 (6.0%) | 10 | 101 | 66 (8.3%) |
StringTie's coverage of the chains Rustle misses: WTC11 under 2 for 70.5%, 2-5 for 24.7%, over 5 for 4.9%; mouse ES 36.0% /
36.9% / 27.1%. In the "intron Rustle has nowhere" class: WTC11 67.9% under 2; mouse ES 37.3% under 2, 40.3% at 2-5, 22.3% over 5.

**Patterns that pass the bar for a rule candidate (size, and on both libraries):**
1. **A true chain with an intron Rustle reports in no transcript**: 1,543 on WTC11 (7.7% of StringTie's correct chains) and
   533 on mouse ES (3.2%). On WTC11 two thirds are single-read chains (under Rustle's floor by design; the one-read tier with
   a partner library, row 1808, is the existing answer and the replicates now provide partners). On mouse ES 63% have
   coverage of 2 or more: Rustle drops real introns that have reads. Not examined before; the next row.
2. **Rustle has only a shorter chain** (611 and 152): the class transfrag extension was built for; extension is off by
   default since row 1811 because it paid on primates only. WTC11 is human: testable with the existing switch.
3. **Piece of a Rustle chain** (IsoQuant's side: 192 and 502): the sub-chain collapse, closed by rows 1839-1843.
The precision side gives no candidate: no class of Rustle's own beyond the shifted sites (3.6%), already measured (row 1842).
Expectation: P1 wrong by 48 chains; the rest held, including "the recall side is where a rule could come from".
