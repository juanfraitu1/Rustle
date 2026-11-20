# Row 1830: does the support profile, used as a FILTER, remove wrong transcripts? (pre-registered 2026-11-20)

User (2026-11-19, handoff `docs/HANDOFF_PSU_TO_FRITZ_2026-11-19.md`, section 8): row 1828
(`docs/SUPPORT_PROFILE_2026-11-19.md`) only detects; test the flags as a filter. Judging data and bars chosen by the user on
2026-11-20. Committed before the profile is run on any judging dataset. Workspace `r1830_support_filter`.

## The filter, fixed
- Profile: `bench/famwise_eval/support_profile.py` as committed (blob `12a78da4`), default parameters (`--jitter 3 --near 20
  --thin-share 0.1 --min-exact 0`, primary records), `--multi-exon-only`. No re-tuning.
- Drop a multi-exon transcript of `PREFIX.final.gtf` when it carries SHIFTED_DOMINANT (at one junction the reads shifted 1-3 bp
  outnumber the exact ones) or THIN_EXACT (the weakest junction has fewer exact reads than 10% of the strongest).
- NOT filtered: WEAK_JUNCTION, NEAR_MISS_DOMINANT, UNCOVERED, NO_FULL_READ (counted per label, descriptive only). Single-exon
  transcripts are neither filtered nor counted.
- The filter is a post-step on a COPY of `final.gtf` (`PREFIX.supportfilter.gtf` in the workspace). No product of the driver changes.

## Judging data (user, 2026-11-20)
| | dataset | reads | truth |
|---|---|---|---|
| J1 | zebrafish, simulated | `reads/DRE/DRE_sim_mm.bam` | `r1800_snap_ism/dre/truth.gtf` (transcripts with >= 1 simulated read): complete |
| J2 | mouse, simulated | `reads/MMU/MMU_sim_mm.bam` | `r1779_simulated_reads/mmu/truth.gtf`: complete |
| J3 | fly, simulated | `reads/DME/DME_sim_mm.bam` | `r1779_simulated_reads/dme/truth.gtf`: complete |
| J4 | mouse ES cells, real (ENCSR172GXL) | `reads/MMU/MMU_ES_mm.bam` | GENCODE vM27 (`genomes/MMU/MMU.ref.gtf`): incomplete |

J4 is judged on the mouse contigs only. The library's SIRV and ERCC contigs are left out of the bars: SIRV built the rule (row
1828). They are reported apart, labelled development. Not judging data: SIRV, HG002.
Use of these four before this row: all were diagnostic or check data for other rules (rows 1779-1819). The support profile has
not been run on any of them. Known already from those rows: chain precision of `final.gtf` (row 1817) and that 18-27% of its
wrong chains have an intron within 10 bp of a true one (row 1818).

## Assembly
A fresh `tools/rustle_pipeline.sh assemble --bam BAM --fasta GENOME --out PREFIX` per dataset at today's defaults, repo at the
commit of this document, binaries built on Fritz on 2026-11-20 from `4e9bb9be`. The row-1817 assemblies are not reused.

## Labels and bars
Unit: multi-exon transcripts of `final.gtf`. TRUE = contig, strand and intron chain equal a truth transcript's
(`bench/famwise_eval/label_by_chain.py`, blob `89b1e494`); FALSE otherwise.
Per dataset:
- **F1, removal:** the filter drops >= 50% of the FALSE transcripts.
- **F2, loss:** the filter drops <= 2% of the TRUE transcripts.
Reported with them: precision and sensitivity (TRUE chains / truth chains) before and after, dropped transcripts by flag.

**Verdict.** PASS = F1 and F2 on all four datasets. Anything else = FAIL as registered, with the per-dataset table; a filter that
meets F2 everywhere but not F1 is reported as "safe, small" with its removal shares, and is not a pass.
**Decision.** A PASS makes the filter a candidate opt-in post-step. No default changes in this row; the user decides.

## Secondary, no bar
- StringTie on the same BAM as the reference point (the de novo run of rows 1779 / 1817 where one exists; a new run with the
  same command for zebrafish): its precision and sensitivity beside ours before and after.
- GffCompare class codes of the dropped and kept FALSE transcripts on J4.
- Each of the two flags alone.

## Limits, stated in advance
- J4's truth is incomplete: real unannotated transcripts count as FALSE, so F1 on J4 is measured against a mixed set.
- J1-J3 are simulated: no truncated reads, precision is flattered (row 1793), and their junction errors come from the simulator's
  error model and the aligner, which may differ from a real library's 1-3 bp shifts.
- Expectation written before the run: F1 is the hard bar. Row 1818 puts near-junction errors at 18-27% of the wrong chains in
  `final.gtf`, below 50%, so F1 passes only if THIN_EXACT also removes wrong chains of other classes.
- THIN_EXACT is relative to the transcript's strongest junction; a true isoform sharing a constitutive junction with a much more
  abundant one can be flagged. F2 measures that cost.

## State at registration
Nothing run. Fritz's `/anvme` has one storage target (`nvme-OST000b`) unavailable since the maintenance of 2026-11-19; J1's BAM
and parts of the Python environment lie on it. The runs start when it answers.

## Result (2026-11-20): FAIL as registered

Jobs 4483990 (mouse simulated), 4483991 (fly), 4483992 (mouse ES), 4483993 + 4483994 (zebrafish); repo `ab761645`, binaries of
2026-11-20; per dataset `r1830_support_filter/NAME/score.txt`.

| dataset | FALSE dropped | F1 (>= 50%) | TRUE dropped | F2 (<= 2%) | precision | sensitivity | StringTie, precision / sensitivity |
|---|---|---|---|---|---|---|---|
| J1 zebrafish simulated | 137 / 1,514 (9.0%) | FAIL | 121 / 12,381 (0.98%) | PASS | 0.8910 -> 0.8990 | 0.6784 -> 0.6717 | not run |
| J2 mouse simulated | 168 / 1,327 (12.7%) | FAIL | 454 / 27,159 (1.67%) | PASS | 0.9534 -> 0.9584 | 0.6255 -> 0.6150 | 0.9586 / 0.5459 |
| J3 fly simulated | 38 / 261 (14.6%) | FAIL | 60 / 7,041 (0.85%) | PASS | 0.9643 -> 0.9690 | 0.7507 -> 0.7443 | 0.9629 / 0.6983 |
| J4 mouse ES, mouse contigs | 6,218 / 33,382 (18.6%) | FAIL | 690 / 17,731 (3.89%) | FAIL | 0.3469 -> 0.3855 | 0.1547 -> 0.1487 | 0.3673 / 0.1431 |

F1 fails on all four datasets; F2 fails on the real library. The filter is not a candidate.
- On the three simulated sets the filter removes 9-15% of the FALSE transcripts, and on two of them it drops MORE true transcripts
  than false ones in absolute numbers (mouse 454 vs 168, fly 60 vs 38; zebrafish 121 vs 137). "Safe, small" by the wording above
  only in the sense that F2 holds; it is not a useful clean-up there.
- On mouse ES it drops nine FALSE (not in GENCODE) transcripts per TRUE one and raises precision by 0.039, at 3.9% of the TRUE ones.
- By flag (transcripts carrying it, TRUE / FALSE): SHIFTED_DOMINANT zebrafish 35 / 115, mouse simulated 387 / 155, fly 31 / 36,
  mouse ES 505 / 3,671; THIN_EXACT 86 / 25, 68 / 17, 29 / 3, 193 / 2,988. THIN_EXACT costs more true than false transcripts on
  every simulated set; SHIFTED_DOMINANT alone would still miss F1 everywhere (7.6-13.8% of FALSE on the simulated sets, 11.0% on ES).
- Not filtered, FALSE transcripts carrying NEAR_MISS_DOMINANT: 26, 72, 15, 2,521 (TRUE: 5, 120, 12, 142).
- The expectation written before the run held: near-junction errors are a minority of the wrong transcripts in `final.gtf`.
  The SIRV development figures (23 of 28 false flagged) did not transfer: SIRV is a deep library with deliberately close
  alternative sites, and its false transcripts are almost all shifted junctions.

**Deviations and what was not done.**
- J1 was run on `reads/DRE/DRE_sim_mm_b.bam`, a new mapping of the same simulated reads (`r1800_snap_ism/dre/reads.fastq.gz`,
  the command of `genomes/map_reads.sbatch`), because the registered BAM lies on the unavailable storage target. Its `final.gtf`
  has the transcript and correct-chain counts of row 1817 on the original BAM (12,381 correct, precision 0.891).
- Tools came from a fresh install of the same `pixi.lock` (`~/.pixi/detached-envs/rustle-14466684255792235947`), same reason.
- StringTie on zebrafish and the GffCompare view of J4 were not run (both tools have files on the unavailable target). Neither
  enters a bar.
- The SIRV md5 check of the 2026-11-20 binaries is still open (its input BAM is on the same target). The "before" numbers of
  all four datasets equal row 1817's `final.gtf` numbers (ES 17,731 correct against 17,730).
- Wall time of the profile: 27-32 min on fly and zebrafish, 2 h on mouse simulated and mouse ES (single-threaded).
