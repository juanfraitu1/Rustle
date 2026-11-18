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

**Noted before any result (2026-11-18): read numbers, and a restart.** Simulated reads handed to the tools: mouse 8,409,652
(47,461 truth transcripts), fly 2,357,595 (13,397). The fly set is smaller than "the real pooled read count" above: the
simulator takes the per-transcript counts as they are, and the fly draw holds 41,610 of the 142,434 mouse values (29% of the
reads). Left as is; the fly set is the same recipe at about 0.29 of the depth. The preparation job was restarted once more
(job 4480691): its FASTA-to-FASTQ step would have passed the 8 h limit; the replacement writes the same bytes (md5 on a 400 MB
sample), and the finished mouse simulation was kept. No read had been mapped.

## Result (2026-11-18): does not transfer as registered

Jobs 4480434 / 4480436 / 4480438, every arm finished (`bad=0`). Scripts and tables: workspace `r1779_simulated_reads/score/`
(`bars.py`, `extras.py`; `dme.txt`, `mmu.txt`, `mmu_es.txt`, `sirv.txt`, `*.extras.txt`).

| dataset | B1 | B2 | B3 |
|---|---|---|---|
| fly, simulated | PASS | FAIL, every cut | FAIL: 0.833 -> 0.974, kept 0.842 |
| mouse, simulated | PASS | FAIL, every cut | FAIL: 0.821 -> 0.916, kept 0.842 |
| mouse ES, real | PASS | FAIL at cov >= 3, 4, 6, 11 | FAIL: 0.320 -> 0.411, kept 0.916 |
| SIRV-Set 4 | PASS | FAIL, every cut | not judged |

Correct chains / precision, de novo:

| arm | fly (truth 9,379) | mouse sim (43,423) | mouse ES (114,584) | SIRV (60) |
|---|---|---|---|---|
| ours `PREFIX.gtf` | 7,058 / 0.833 | 27,175 / 0.821 | 17,714 / 0.320 | 47 / 0.398 |
| ours `PREFIX.final.gtf` | 5,945 / 0.974 | 22,892 / 0.916 | 16,233 / 0.411 | 41 / 0.547 |
| StringTie | 6,549 / 0.963 | 23,705 / 0.959 | 16,401 / 0.368 | 35 / 0.814 |
| IsoQuant | 6,592 / 0.972 | 24,098 / 0.967 | 13,443 / 0.472 | 39 / 0.848 |
| FLAIR | 6,614 / 0.720 | 27,048 / 0.535 | 16,686 / 0.096 | 48 / 0.081 |

**The limit.** Ours has the most correct chains on every dataset but SIRV (one fewer than FLAIR), and is never dominated. Its
precision is below StringTie and IsoQuant on all four, at every read-count cut on the simulated sets and SIRV. The filter trained
on the 12 primate samples raises precision everywhere but drops 16% of the correct chains on simulated reads (8% on real mouse),
more than the 7% allowed.

Recall by simulated reads of the chain (mouse; 1 / 2 / 3-5 / 6-20 / > 20): ours 0.006 / 0.433 / 0.583 / 0.699 / 0.910;
StringTie 0.016 / 0.346 / 0.409 / 0.533 / 0.858; IsoQuant 0.003 / 0.376 / 0.457 / 0.579 / 0.844; FLAIR 0.023 / 0.060 / 0.715 /
0.772 / 0.907. Fly, same order: ours 0.006 / 0.711 / 0.850 / 0.934 / 0.977; StringTie 0.012 / 0.649 / 0.704 / 0.818 / 0.951;
IsoQuant 0.004 / 0.683 / 0.749 / 0.841 / 0.936; FLAIR 0.033 / 0.059 / 0.937 / 0.952 / 0.960.
Without transcripts > 15 kb recall moves by at most 0.001 for every arm (117 mouse chains, 51 fly chains leave).

Wrong chains of ours `PREFIX.gtf`, first class that fits. Mouse, 5,940: new intron within 10 bp of a true one 3,075; contiguous
piece of a true chain 1,275; true introns in a combination not simulated 1,036; new intron elsewhere 298; no true intron 256.
Fly, 1,418: no true intron 835; within 10 bp 278; piece 195; combination 56; new intron elsewhere 54.
Not examined further: these samples are not mined.

Guided block (not judged; guide = truth with 5/16 hidden; hidden chains found / precision of non-guide chains). Fly, 2,901 hidden:
ours 1,714 / 0.656; StringTie 1,743 / 0.907; IsoQuant 1,992 / 0.964; FLAIR 134 / 0.635. Mouse simulated, 13,430 hidden: ours
6,109 / 0.763; StringTie 5,749 / 0.900; IsoQuant 7,185 / 0.961; FLAIR 2,054 / 0.717. Mouse ES, 35,741 hidden: ours 4,274 / 0.230;
StringTie 3,847 / 0.299; IsoQuant 4,085 / 0.220; FLAIR 2,322 / 0.055.
