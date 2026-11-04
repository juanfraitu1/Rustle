# PREREG MISSING-COPY HALLMARKS ON A NON-T2T GENOME (register 1537, 2026-11-04, before any number)

Question (user 2026-11-04): "comparing T2T vs non-T2T for certain gene families, what are the hallmarks of possibly missing
copies?" r1504 showed the detectors (flag, candidates) find 1/32 and 0/37 missing copies; its verdict said a detector must start
from reads that fail to place. The KB3781 mat/pat study (2026-10-08) already contradicts the literal reading: reads of a missing
copy are essentially never unmapped; they are absorbed on a near-identical relative (de ~0.001) or form a divergence pile. This row
measures, per expressed missing copy, what its reads do on the incomplete genome, against the same measure on present copies.
Data: the r1504 products as they stand (no realignment): `/anvme/workspace/v138da13-rustle/nont2t/{human,gorilla}`: `nt.bam`
(panel reads on GRCh38 / Kamilah, `-N 50 -p 0.1 --secondary=yes`), `truth2.tsv` (copy status by sequence, present loci),
`panel.genes.tsv`. Read origin = the T2T BAM primary record (`-F 0x904`) overlapping the copy's exons; a read overlapping two
panel copies is dropped (ambiguous origin; human 1,751 of 39,808).
Per read on the non-T2T genome (primary record): unmapped; MAPQ 0; tie = at least one secondary with AS equal to the primary;
absorbed = primary overlaps the locus of a PRESENT panel copy (truth2 `locus`), else elsewhere; de gap = de(non-T2T) - de(T2T);
soft clip = max clipped bases / read length.
Per copy (median over its reads, shares as fractions): unmapped share, MAPQ0 share, tie share, absorbed share, median de gap,
median clip share; plus the absorbing locus (the present copy taking most of its reads) and depth excess there = reads absorbed /
own T2T reads of that present copy.
Null: the same per-copy numbers on the expressed PRESENT copies (human 113). Threshold per feature = the 95th percentile of the
present copies (de gap, tie, MAPQ0, unmapped, clip: above; absorbed share is descriptive, not a hallmark, since present copies'
reads also land on their own locus). A copy is HALLMARKED when at least one of the five features exceeds its threshold.
Scorable = copies with >= 10 origin reads (human: 23 of 32 missing); smaller ones listed, not scored.
Bars (human, development): hallmarked >= 50% of scorable missing copies, split collapsed vs absent (reported both), AND present
copies hallmarked <= 25% (five features at the 95th percentile give ~23% by construction if independent; above 25% the rule is
noise). Gorilla = held out: the human thresholds applied once, reported as a second table, carries no bar (37 missing, many with
no homolog at all: a different regime).
Stop rule: PASS -> the hallmark set is the registered input of a read-level detector (next row, not this one); FAIL -> missing
copies are invisible at the read level on this panel and O3 on non-T2T genomes closes on reads too.
Script: `bench/nont2t/hallmarks.py` (reads read_origin.tsv + nt.bam + truth2.tsv -> per-read and per-copy tables).
