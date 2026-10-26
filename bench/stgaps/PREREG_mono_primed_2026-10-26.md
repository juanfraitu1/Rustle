# Pre-registration: drop internally primed single-exon transcripts (StringTie rlink.cpp:486/605, genomic form)

**Written 2026-10-26, before any held-out run.** Branch `assembler-stringtie-gaps`. A default flip is the user's call.

## Dev evidence (seen; never in the verdict)
Human A119b chr20 (`/mnt/linuxdisk/tmp/asm_gap/h20_du_nopolish.gtf` polished by `polish_port.py`, = shipped):
458 single-exon transcripts, 1 `=`. StringTie's exact read-sequence predicate flags 4 (inert). The genomic predicate
below flags 246 at share ≥ 0.8 (0 `=`); read-level primed rate 0.739 vs 0.066-0.081 with ends shifted ±1 kb / +5 kb.
Gorilla NC_073244.2 dev: 10 single-exon transcripts — inert there.

## Rule (frozen: `bench/stgaps/polya_mono_genomic.py` predicate, transcript level)
For each single-exon transcript of the polished GTF: reads = unspliced primary records (flag & 2820 == 0, no N)
overlapping ≥ 50% of it; a read is primed when the genome 20 bp after `reference_end` is ≥ 60% A or holds A6, or the
20 bp before `reference_start` is ≥ 60% T or holds T6. Drop the transcript when ≥ 1 read and primed/reads ≥ 0.8.
Multi-exon transcripts are untouched, so intron-chain metrics are unchanged by construction (checked: IV).

## Held-out
Human A119b `chr11` and `chr12` (`A119b.t2t.bam`, never used for any polish/assembly dev; chr16/chr20-22 were).
Arms: BASE = frozen `copy_assign` (sha1 685801dd…, main 3c64bc2a build) `--assemble-only --assembly-polish full
--polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7
--polish-retained-ratio 10`; DROP = BASE minus the rule's transcripts. gffcompare 0.12.10 vs CHM13 RefSeq
(`benchmark_collapse/ref/A119b.chm13v2.0_RefSeq.gtf.gz`, that chromosome). NULL = the same rule with every read's
ends shifted +1 kb.

## Clauses (each held-out chromosome)
- **IV**: DROP's multi-exon transcripts == BASE's (byte-identical lines); else STOP.
- **H1 (signal)**: flagged count under the rule ≥ 5 × flagged under NULL.
- **H2 (no loss)**: single-exon `=` transcripts dropped ≤ max(1, 10% of BASE's single-exon `=`).
- **H3 (gain)**: gffcompare transcript-level precision DROP > BASE.
**PASS** = IV + H1 + H2 + H3 on both chromosomes → wire as opt-in `--polish-mono-primed` (default off; flip is the
user's call). Any failure → register row, not built.
