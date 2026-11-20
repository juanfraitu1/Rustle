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
