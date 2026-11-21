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
