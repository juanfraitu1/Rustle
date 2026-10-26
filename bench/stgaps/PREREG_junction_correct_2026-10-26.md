# Pre-registration: read-level junction correction — StringTie's long-read rule (ST) and a sequence-checked form (SEQ)

**Written 2026-10-26, before any run of either rule.** Branch `assembler-stringtie-gaps`. A default flip is the user's call.

## Prior evidence (seen)
Ledger LEDGER_GENERAL.md:604-700: `--read-junction-snap 25 / 10` (whole-junction snap to a better-supported junction)
cost chain Sn 3.2 / 2.6 and 2.3 / 1.4 points (gorilla NC_073244.2 / human chr20) with no precision gain. Register 1140:
in near pairs (<= 10 bp) 98.5-99.6% of reads align strictly better at their own site.

## Rules (`bench/stgaps/junction_correct.py`, applied to the BAM; the assembler is the frozen binary, unchanged)
Reads: primary records (flag & 2308 == 0) with >= 1 N. Sites: a junction's LEFT site (intron start) and RIGHT site
(intron end), keyed by (position, motif strand: GT-AG = +, CT-AC = -, else 0). Site support = number of such reads with
an intron at that site (StringTie `leftsupport` / `rightsupport`).
- **ST** (rlink.cpp 15300-15530 in `-L` mode, where every junction is "mismatched" so every one is eligible): each
  site independently moves to the same-strand site within < 25 bp (`sserror`) with the largest support, provided
  support_target × 0.9 (`1 - ERROR_PERC`) > support_own; ties → nearer, then left. One step (no chains). Not applied
  when the flanking CIGAR ops are not M/=/X of length > |shift| (counted).
- **SEQ**: ST's candidate, accepted only if the read's bases over the shifted segment have <= mismatches at the new
  placement than at the old (equal = moved: a sequence-equivalent placement).
A moved junction keeps the other site; CIGAR lengths are adjusted, query sequence and alignment start unchanged.

## Substrates
Dev (seen, never in the verdict): gorilla NC_073244.2 (GGO_mm.bam), human A119b chr20.
**Held-out:** human A119b **chr11, chr12**; gorilla **NC_073230.2** (GGO_mm.bam; in no register row).
Assembler: frozen `copy_assign` sha1 685801dd `--assemble-only --assembly-polish full --polish-isoform-fraction 0.02
--polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7 --polish-retained-ratio 10`, defaults
otherwise, on BASE (original BAM) and on each corrected BAM. Score: gffcompare 0.12.10 vs RefSeq
(`benchmark_collapse/ref/*_RefSeq.gtf.gz`, that contig).

## Clauses (per held-out substrate; chain = intron-chain level)
- **G (gain)**: chain Pr(arm) − Pr(BASE) >= +0.5 AND chain Sn(arm) − Sn(BASE) >= −0.3.
- **N (no loss)**: matching intron chains(arm) >= matching intron chains(BASE) − 1%.
**ADOPT as opt-in** = G and N on all 3 held-out substrates. Otherwise a register row, not built.
Predicted: ST fails N (as the 25 bp snap did); SEQ moves few reads and is near-inert (fails G).
