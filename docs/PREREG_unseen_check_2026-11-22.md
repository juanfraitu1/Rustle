# Row 1887: the two defaults of row 1886 on libraries neither rule has seen (pre-registered 2026-11-22)

Row 1886 made the short-exon repair a default of `assemble` and writes `PREFIX.ranked.gtf`. Both were judged on datasets that
were in hand; the user's standing rule (2026-11-10 / reuse ruling) asks for one more unseen dataset per new rule. ENCODE
ENCSR507JOF WTC11 replicates 1 and 3 (`reads/HSA/WTC11_r1_mm.bam`, `_r3_mm.bam`) were mapped as check data in row 1862.
Replicate 3 was assembled once (row 1866, one-read promotion) and replicate 1 was a partner library there; neither has been run
with the repair, and no other assembler has been run on either. Human, so this checks transfer to unseen libraries, not to
another species (the non-primate side of both rules is rows 1878 and 1880: 7 of 7 datasets each).

## Runs (driver of commit 1886, nothing set unless said; `RUSTLE_PTC=0`)
Per replicate: `assemble` with `RUSTLE_SHORTEXONS=0` (before), `assemble` (default), `oneread` on the default; StringTie 3.0.3
`-L` and IsoQuant 4.1.0 with the commands of `apetools/job.sbatch`. Truth: RefSeq on CHM13 (`ptr_ppa_rules/gffcmp/hsa.ref.gtf`).
Workspace `r1887_unseen_check`.

## Reading, fixed now
- **Repair CONFIRMED** on a replicate when the default `final.gtf` has at least as many correct chains and at least the
  precision of the run without it, one strictly better (row 1878's "pays"). Confirmed overall = both replicates.
- **List CONFIRMED** on a replicate when `PREFIX.ranked.gtf`, cut at StringTie's and at IsoQuant's number of multi-exon chains,
  holds at least as many correct chains as each, and its `list_sensitive` set has at least StringTie's correct chains and at
  least its precision. Confirmed overall = both replicates.
- A miss is reported as a miss; it does not by itself undo a default (the user's call).
