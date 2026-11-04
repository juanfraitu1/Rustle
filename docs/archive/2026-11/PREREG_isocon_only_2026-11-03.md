# PREREG ISOCON-ONLY TRANSCRIPTS (register 1515, 2026-11-03)

**Written before any number.** User (2026-11-03): "is there any of the clustering procedures that we could use for transcripts
that only isocon can find?" -> "yes lets do that". Closes the IsoCon question: which reference transcripts did IsoCon match on the
r1482 panel that our assembly did not, and does OUR cluster -> consensus -> significance-merge chain (`candidate_copies`, the port of
IsoCon steps 1, 3, 4) recover them when run on IsoCon's pools with every consensus emitted.

## What is known (not predicted)
- r1482 (A119b chr20, 370 single-copy genes): distinct CAT reference transcripts matched by an FSM: ours 631, IsoCon 645 (a net,
  not the set difference); IsoCon artifact share 11.0 vs ours 5.0; 52 candidates per gene, 42% ISM, 5,484 mono-exon.
- r1381: IsoCon's nearest-neighbour read correction gave nothing over correlation clustering on real reads.
- `candidate_copies` differs from IsoCon upstream of the merge: the initial partition is at delta 0.00958 (a copy-divergence
  threshold), the net is subsampled (`--max-reads`), pairs are tested only at minimizer-sketch share >= 0.5, and in-reference
  consensuses are dropped as `InReference` (never written). This row names what it tests: OUR chain on IsoCon's pools, all
  consensuses kept; NOT IsoCon's correction step (refuted) and NOT IsoCon's merge-to-exhaustion loop.
- The r1482 `ours` arm is `a1.gtf` of 2026-10-29, before the end-polish defaults (r1485-1488). `~/sqanti_sc_A119b/gffc/` holds
  gffcompare tracking of `ours_tssr`, `ours_med`, `ours_ext`, `ours_extmed` against `ref370.gtf` (the same CAT transcript ids).

## Part A: the Only set (deterministic from the r1482 files, computed AFTER this file is written)
- `S_arm` = distinct `associated_transcript` of rows with `structural_category == full-splice_match` in the UNFILTERED
  `~/sqanti_sc_A119b/sq3/{arm}/{arm}_classification.txt`. **Only = S_isocon \ S_ours.** Reported in two strata: multi-exon and
  mono-exon reference transcripts (SQANTI's mono-exon FSM is overlap-based).
- Per Only transcript (table `only.tsv`): gene, reference exon count, IsoCon's best FL support on it (`FL` column), ours' best
  category on the same gene (`ISM` of this transcript / `ISM` of another / NIC / NNC / none), and the DRIFT column: whether
  `ours_tssr` / `ours_med` / `ours_ext` / `ours_extmed` already match it (gffcompare class `=` in `gffc/*.tracking`). The drift
  column is a report, not a bar: if the end-polish defaults already match most of Only, the clustering question is dead and the
  write-up says so.

## Part B: our chain on the 370 pools, every consensus kept
- Code: `candidate_copies --all-clusters PATH` (opt-in, this row): after phase 1 (clusters -> consensus -> refinement -> significance
  merge) every family's every consensus is written to PATH (`>family:cK n_reads=N`), BEFORE the genome fates drop the in-reference
  ones. Without the flag: byte-identical products. The flag is no stage-key input and the run never replays from the cache.
- Input: one synthetic single-copy family per universe gene (`bench/sqanti_single/isocon_only.py copies`): `family_id` = the CAT
  gene id, copy 0, exons = the UNION of the gene's CAT exons on chr20 (so exon-skipping reads still attribute: `ATTRIB_MIN_READ_COV`
  0.5), `copies.fa` = the spliced exon-sum sequence in transcription orientation (the catalog's own format). Net = every primary and
  secondary record on the copy (the default pool, as IsoCon's). `--max-reads 5000` (the universe's depth ceiling: no subsampling).
  Defaults otherwise (`--delta 0.00958 --min-cluster 3 --min-support 6`). BAM `A119b.t2t.bam`, genome `chm13v2.0.fa`, a fresh
  `splice:hq` index in the workspace.
- Scoring = the r1482 protocol unchanged: `all_clusters.fa` mapped with the `isocon_targets.py gtf` command (`minimap2 -ax splice:hq
  -uf --eqx -N 50 -p 0.1 --secondary=yes`, primary alignments only), `sc_panel.py cut` to the universe spans, SQANTI3 6.1 `qc` +
  `filter rules` against the same `chr20_ref_cat.gtf` / `chr20.fa`, `sc_panel.py table` with the arm as a fifth row over all 370
  genes. `isocon_only.py recover` then lists, per Only transcript, whether an arm FSM names the same `associated_transcript`.

## Decision rule (fixed now)
1. **Recovery**: >= 50% of the multi-exon Only transcripts are matched by an arm FSM (same `associated_transcript`).
2. **Not WORSE by the r1482 rule**: the arm's artifact share <= ours + 5 = 10.0 points (FSM+ISM share reported, no bar: a
   candidate set is not an assembly).
3. **Wall**: the `candidate_copies` run finishes within `timeout 7200` on one node (IsoCon: 370 pools in 22 min).
PASS all three -> nothing is adopted here; the result goes to a design decision on an opt-in isoform arm of O3 (consensus
candidates are an O3 product, not the assembly). Any bar FAILS -> the IsoCon question is closed in the register: our chain does
not find what IsoCon found, or finds it at IsoCon's artifact cost. Part A is reported either way.

Run dir `~/isocon_only_A119b/` (job outputs), big files in workspace `isocon_only/`. Code: `src/candidate_copies.rs` (the flag),
`bench/sqanti_single/isocon_only.py` (`set | copies | recover`), `~/isocon_only_A119b/run.sbatch`.

## RESULT (2026-11-03, job 4421576; `~/isocon_only_A119b/{only.tsv,table.tsv,recover.tsv,sq3/arm}`, consensuses workspace `isocon_only/all_clusters.fa`)
**Part A.** S_ours 631, S_isocon 645; **Only = 154** (125 multi-exon, 29 mono-exon); the reverse set (ours-only) is 140, so the
645 vs 631 was a net of two large differences. Ours' best model on the Only gene: NIC 68, ISM of another transcript 39, NNC 11,
FSM of another 4, ISM of this one 2, nothing 1. Drift: ours_tssr / med / ext / extmed / all already match 3 / 4 / 3 / 4 / 4 of the
125 (the end-polish defaults do not close the gap). IsoCon's own read support on the multi-exon Only transcripts: median 3 reads
(63 of 125 at <= 3, 7 at >= 10).
**Part B.** candidate_copies wall 2,605 s (bar 3 PASS); 6,808 consensuses, 1,804 inside the universe spans after mapping.
| arm | transcripts | genes | FSM | ISM | NIC | NNC | artifact | FSM+ISM % | artifact % | rules PASS % |
|---|---|---|---|---|---|---|---|---|---|---|
| arm (all consensuses) | 1,804 | 368 | 314 | 123 | 161 | 157 | 1,049 | 24.2 | **58.1** | 41.4 |
Artifacts: genic_intron 759, genic 195, antisense 53, fusion 31. **Recovery: 12 of 125 multi-exon Only (9.6%)** (and 9 of 29
mono-exon). **❌ FAIL on bars 1 and 2.** Our cluster -> consensus -> merge chain on IsoCon's pools finds one in ten of what IsoCon
found, at an artifact share ten times ours: the consensus of a delta-0.00958 cluster is a copy-level object (intronic and
partial pieces dominate), not an isoform. The IsoCon question is closed: its single-copy isoform gain comes from its own
exhaustive per-read candidate loop, with support floors of 2-3 reads, and costs 2x artifacts; nothing of its clustering remains
to port. `--all-clusters` stays as a diagnostic hatch.
