# PREREG ISOSEQ-CLUSTER2 (2026-11-16, register 1759, written before any run)

Question (user 2026-11-16): IsoCon is a poor "purely de novo" comparator (r1478 WORSE, r1515 closed). Is `isoseq cluster2`
(PacBio's reference-free clustering of FLNC reads, isoseq 26.2.0, already in the sqanti3 env) the better one, and how does it compare with
our default assembly? Descriptive comparison of tools, no design change, no held-out.

Sample: gorilla OR6737 testis (the r1478 sample). Input `~/reads/GGO/GGO_OR6737.IS.8a36218d3f23-filtered.bam` (4,409,946 FLNC, lima + refine
by the producer). Whole sample, no read pools (IsoCon needed alignment-defined pools; cluster2 does not).

Arms.
- cluster2: `isoseq cluster2 -j 72 FLNC.bam transcripts.bam`, defaults (no `--singletons`); consensus transcripts mapped with the project
  flags (`minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`) to GGO.fa; primary records -> GTF (exons from the CIGAR,
  `reads` = FLNC count of the cluster). `bench/denovo_eval/cluster2_cmp.py`.
- ours: today's default `assemble` on `GGO_mm.bam` (`PREFIX.gtf`), rerun in the same job.
- IsoCon column in part 2 from the r1478 run (`~/work/jfim/runs/isocon_OR6737/isocon.gtf`), unchanged.

Gate (else INVALID): >= 0.90 of cluster2 transcripts have a primary alignment.

Part 1, whole sample, exact intron chains vs RefSeq (`curve.py` definition; multi-exon chains only). Both sets cut to the same number of
chains K = min(ours, cluster2), each ranked by its own read support. cluster2 BETTER if its correct chains at K > 1.02 x ours; WORSE if
< 0.98 x ours; TIED otherwise. Full-set chains / correct / precision reported beside it.

Part 2, the r1478 targets (`isocon_OR6737/targets.tsv`: 4 Yoo families, 30 single-copy genes), `toolcompare.py --targets`, same verdict rule
as r1478 per class against our arm: BETTER if it models more targets with junction precision >= ours - 0.05; WORSE if fewer targets, or
junction precision below ours - 0.05 at equal targets; TIED otherwise.

Reading. cluster2 replaces IsoCon as the reference-free comparator in the tool table if it is not WORSE than IsoCon in part 2 (multi-copy
targets modeled >= IsoCon's 20). Whatever the verdicts, nothing in Rustle changes.

Prediction: part 1 WORSE (clusters keep truncated and intron-retaining molecules as separate transcripts, and near-identical copies
collapse into one consensus); part 2 multi-copy WORSE than ours, about level with IsoCon.

Workspace `r1759_cluster2/`.

## Result (2026-11-16, job 4473746)
See register row 1759: gate PASS; part 1 WORSE (0.745 x ours at matched count); part 2 WORSE than ours in both classes (multi-copy 23 vs 27, single-copy 29 vs 30), above IsoCon (20 / 27). cluster2 replaces IsoCon as the reference-free comparator. Output `r1759_cluster2/{job.4473746.out,score.txt,cluster2.gtf}`.

## Follow-up, post hoc (2026-11-16; user: can a cluster2 strategy be adopted, or cluster2 be the purely de novo assembler?)
The 1,321 correct chains only cluster2 emits, by the number of primary reads of `GGO_mm.bam` with that exact chain: >= 2 reads 1,114
(our filters dropped them: the threshold trade closed held out by rows 1616-1618), 1 read 205 (the consensus recruits a second, slightly
erroneous read: about 0.8% of our correct chains), none 2. Reading: clustering without a genome recovers nothing the alignment route
does not already see; no strategy to port. cluster2 is the external genome-free column of the tool table, not a Rustle mode. One sample;
its ~92,800 wrong chains were not examined.
