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

## RESULT r1537 (2026-11-04): FAIL as registered on human; gorilla held out separates for a different reason
Products `/anvme/workspace/v138da13-rustle/nont2t/{human,gorilla}/hm.{reads,copies,thresholds}.tsv`, `read_origin.tsv`.
Human origin reads 37,014 (1,751 ambiguous dropped); thresholds (95th percentile of 96 scorable present copies):
unmapped 0, MAPQ0 1.00, tie 1.00, de gap 0.0014, clip 0.0008.
| pair | copies | scorable (>= 10 reads) | hallmarked | bar |
|---|---|---|---|---|
| A119b -> GRCh38 | missing | 21 | 5 (0.24): collapsed 4 of 14, absent 1 of 7 | >= 0.50 **FAIL** |
| A119b -> GRCh38 | present | 96 | 9 (0.09) | <= 0.25 pass |
| GGO -> Kamilah (held out, human thresholds) | missing | 17 | 17 (1.00): collapsed 4 of 4, absent 13 of 13 | no bar |
| GGO -> Kamilah | present | 7 | 1 (0.14) | no bar |
What the reads do (descriptive, as registered): no human missing-copy read is unmapped; 15 of 21 scorable missing copies have
most reads absorbed on a present panel copy (absorbed share 0.39-1.00), at identity indistinguishable from that copy's own reads
(median de gap <= 0.002 for 19 of 21; the null is wide because whole GOLGA8 / TBC1D3 present loci differ by 0.002-0.003 between
CHM13 and GRCh38). Depth excess at the absorbing locus: CDY1 6.0, TSPY3 2.5, GOLGA8DP 1.1, the rest < 0.6; a present copy
(TBC1D3P3) shows 2.8, so depth excess alone does not separate either. Tie and MAPQ0 cannot be hallmarks on this panel: the present
null is saturated (Y-ampliconic present copies have every read tied), threshold 1.00.
Gorilla: all 17 separate on de gap and clip, and still 17 of 17 with gorilla's own present null (hm_own: thresholds de gap
0.0010, clip 0; present hallmarked 3 of 7). Present Kamilah copies place at mGorGor1 identity (de gap -0.0007 to 0.0010, one at
0.0075); missing ones form a divergence pile (0.0027-0.097) or are fully unmapped (LOC129532042, LOC134758753, LOC129529226;
61 reads, all flag-4 records in nt.bam). That is the divergence-pile hallmark of the KB3781 mat/pat study for the ABSENT class.
Two regimes, then: a copy whose nearest relative is far (absent on a draft) separates on identity; a copy collapsed onto a
near-identical relative (human GRCh38) does not.
### After the numbers (labelled, no bar changed)
1. Three human "missing" copies with absorbed share 0 are truth mislabels: GOLGA6L9 (579 reads) hits GRCh38 chr15:82430017-82439153
   at identity 1.000 and its reads sit at 82436415-82438701 inside that hit; GOLGA8M (171) hits chr15:28701953-28738384 at 0.9991,
   reads at 28668250-28716368; GOLGA6L17P (80) likewise at chr15 82 Mb. The r1504 claim rule gave each locus to a sibling with an
   equal hit (GOLGA6L11P, the other GOLGA8M gene id) and the claimer's locus span is shorter than the hit, so the reads fall outside
   it. GRCh38 has these copies; 830 reads sit in the wrong denominator. Excluding them: missing scorable 18, hallmarked 5 (0.28;
   collapsed 4 of 12, absent 1 of 6): still FAIL, bar unchanged.
1b. The ambiguous-read drop (a read over two panel copies) removes reads at adjacent near-identical copies, the class this study
   targets: GOLGA6L16P 750 -> unscored, GOLGA6L9 1,329 -> 579, RBMY1B 221 -> 45, RBMY1A1 209 -> 33, GOLGA8M 303 -> 171; hence 21
   scorable, not the 23 the PREREG stated. Bias noted; verdict unchanged.
2. Depth excess is the only feature that singles out the Y copies (CDY1, TSPY3) and it is confounded by present TBC1D3P3 (2.8).
Verdict: on a same-species near-complete genome (GRCh38) the reads of a missing copy are absorbed at near-identical identity, with
no read-level signature that the present copies lack: the COLLAPSED regime closes on reads (stop rule). On a draft genome (Kamilah)
ABSENT copies separate on identity against both nulls (17/17, 3/7 present false alarms with the gorilla null): that regime has a
read-level hallmark (the divergence pile), already seen on KB3781 mat/pat; it carries no bar here and would need its own row.
