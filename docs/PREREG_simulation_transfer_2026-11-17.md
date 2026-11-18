# Rows 1779 and 1780: does the assembler transfer to data it was never tuned on? (pre-registered 2026-11-17)

Advisor (meeting 2026-11-16): test the assembler with simulation. User (2026-11-17): run simulated reads (row 1779) and the LRGASP
mouse ES library with its SIRV spike-ins (row 1780); new samples are allowed for this step only.
Registered before any read was simulated or mapped. Nothing here changes a default.

**Why these data.** Every threshold and both chain scores were set on 12 primate samples. These four datasets share no species,
library or laboratory with them:

| dataset | row | reads | truth |
|---|---|---|---|
| mouse, simulated | 1779 | IsoSeqSim as in LRGASP (Pardo-Palacios 2024), GRCm39 + GENCODE vM27 | transcripts with >= 1 simulated read |
| fly, simulated | 1779 | same, BDGP6.46 + Ensembl 112 | same |
| mouse ES, real | 1780 | LRGASP cDNA-PacBio, ENCODE ENCSR172GXL, three replicates pooled | GENCODE vM27 (incomplete) |
| SIRV-Set 4 | 1780 | the spike-in reads of the same library | SIRV annotation ENCFF873AAI (complete by construction) |

**Simulation recipe (fixed).** The IsoSeqSim copy and the call of `LRGASP/lrgasp-simulation` (`src/simulate_pacbio.py`): Sequel
5' / 3' completeness tables, `--es 0.004 --ei 0.006 --ed 0.006`, `--polya`. Expression, mouse: `quantify.py` of that repository
on the real mouse ES reads (minimap2 to the vM27 transcripts, primary alignments counted). Read number = the real pooled read
count. Fly has no real library here: the mouse per-transcript counts, zeros included, are drawn without replacement onto the fly
transcripts (seed 1779). Poly-A tails are trimmed (`cutadapt --poly-a`) before mapping, as an Iso-Seq pipeline does.
No novel-isoform insertion step: for a de novo run every transcript is unannotated.

**Runs.** The shipped mapping (`genomes/map_reads.sbatch`) and the shipped all-tools job (`apetools/job.sbatch`), unchanged,
at today's defaults: ours `PREFIX.gtf`, ours `PREFIX.final.gtf`, StringTie `-L`, FLAIR, IsoQuant, each de novo. The guided arms the
job also produces (guide = truth with 5/16 of the transcripts hidden) are a separate labelled block, not judged.
A crash may be fixed; no parameter of any tool is changed after a result is seen. These samples are never mined for a rule.

**Unit.** Exact intron chains of multi-exon transcripts (`bench/famcmp/curve.py`). SIRV = chains on the SIRV contigs only.

**Bars**, per dataset:
- **B1, not dominated.** No other de novo tool at its default has both more correct chains and higher precision than ours
  `PREFIX.gtf`.
- **B2, curve against curve.** For every point of StringTie's coverage curve (cov >= 0, 2, 3, 4, 6, 11, 21) ours has a point on
  its read-count curve (same cuts) with at least as many correct chains and at least the same precision.
- **B3, the trained filter transfers** (not SIRV: too few chains). `PREFIX.final.gtf` against `PREFIX.gtf`: precision higher AND
  at least 0.93 of the correct chains kept (0.932 is the lowest share on the 12 fit samples).

**Verdict.** TRANSFERS = B1 and B2 on all four datasets and B3 on the three large ones. Otherwise the failed bars are named per
dataset; a failed bar is reported as a limit of the assembler, not repaired on these data.
Reported without a bar: recall and precision of every arm; recall by simulated read count (1, 2, 3-5, 6-20, > 20); new and wrong
chains of ours by class on the simulated sets; SIRV sensitivity and precision per tool; the guided block.

**Limits, stated in advance.** Simulated reads carry no degradation, nascent RNA, retained introns or template switching, so
absolute precision flatters every tool; only the comparison on identical reads is evidence. The fly expression profile is
transplanted, so its isoform ratios inside a gene are arbitrary. The SIRV annotation is public and its loci are simple
(the paper's own caveat). GENCODE mouse is incomplete, so "wrong" on the real library includes real unannotated chains.
bambu is not in the all-tools job; it is added only if its existing script runs unchanged.

**Change before any result (2026-11-18).** The first simulation job would have passed its 8 h limit: IsoSeqSim gives a whole
transcript to one worker, and ENSMUST00000020488.9 (123 kb, 16,029 reads) alone needed about 8 CPU-hours. The local copy now hands a
transcript out in pieces of about 500 kb of reads and reseeds each worker (forked workers shared one random state). Read model,
error rates, counts per transcript unchanged (chr4 test: counts identical apart from the reads the simulator drops as empty).
No read had been mapped. Original file kept as `py_isoseqsim_simulate_reads_normal.py.orig`.

**Limit added before any result (2026-11-18): over-long simulated reads.** `quantify.py` counts a read as one copy of the
transcript it lands on, however little of it the read covers, and the simulator then writes that transcript nearly full length.
In 300,000 real reads (261,632 mapped) a read covers 0.74 of its transcript on average, but 0.11 for transcripts > 15 kb (0.41%
of reads) and 0.03 for > 30 kb (0.21%). Worst case ENSMUST00000020488.9 (*D10Wsu102e*, last exon 122,802 bp): 463 of its 512
reads cover one 4.9 kb stretch (transcript 34.4-39.3 kb), yet it is simulated as 16,029 reads of 123 kb, about 6.6% of all
simulated bases. Left as is: it is the LRGASP recipe and every tool gets the same reads. Consequence: for transcripts > 15 kb the
simulated mouse set is easier than real data, and their "correct" chains say little. Recall is also reported with transcripts
> 15 kb left out; the bars stay on the full set.
