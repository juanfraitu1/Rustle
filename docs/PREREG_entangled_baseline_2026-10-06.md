# Pre-registration: what StringTie and FLAIR do at overlapping genes, on ideal reads, and what "better" means (2026-10-06)

**Written before any StringTie, FLAIR or f1units output exists, and before `arms_score.py` was run on any arm.** Status of what was already seen: the default assembler products of the four ideal-expression runs exist and their per-COPY chain recovery (E1) and locus scores (E2-E4) are known (`docs/IDEAL_EXPRESSION_DEFAULT_2026-10-06.md`); the window-wide, per-GENE numbers of this study are new. DEV, human CHM13 v2.0 / CAT-Liftoff only, simulation, circular by construction (the reads come from the annotation that scores them).
User request (2026-10-06): "we need to be able to model the cases in which there are overlapping genes, ideally we should be able to provide what stringtie and flair do for those loci and then do better than that." Tools: `bench/entangled/` (`arms_score.py`, `pool.py`, `run.sh`, `flair_shims/`).

## 1. Substrate and arms

Reads: the four ideal-expression read sets (NPIP and TBC1D3 windows, replicates 1 and 2, seeds 20261006 / 20261007; 20,660 and 16,680 reads of 2,066 and 1,668 transcripts of 547 and 519 genes; 10 full-length jittered reads per transcript, error .001), the BAM `reads.bam` already used by the pipeline arm (minimap2 splice:hq, -N 50 -p 0.1, secondaries kept). Every arm sees the same alignments.

| arm | what | recipe |
|---|---|---|
| **S** | StringTie, annotation-free | `stringtie -L -p 4` on `reads.bam`, no `-G` (the lab's `benchmark_collapse/run_stringtie.sbatch`) |
| **F** | FLAIR, annotation-free | FLAIR's own BAM-to-BED12 filtering (`dofiltering`, primaries only, supplementary removed, quality 0) then `flair collapse --trust_ends --generate_map`, no `flair correct`, no annotation (the lab's `run_flair.sbatch`) |
| **D_asm** | our assembler's transcripts, default | `asm.gtf` of the registered default arm (HEAD binaries, f1v2 regroup, strict junctions, GOOD secondaries) |
| **D_fam** | the same after the families stage's input regroup | `asm.families.gtf` (bridge transcripts removed, gene_ids split) |
| **U_asm / U_fam** | the paused container-v2 lever, opt-in code at HEAD | the driver with `RUSTLE_BRIDGE_REGROUP=f1units` (each bridge transcript replaced by its units), `units.gtf` / `units.families.gtf` |

Versions (drift from the lab runs stated, not corrected): StringTie 3.0.3 (lab 3.0.1 and 3.0.3), FLAIR 3.0.0 (lab 3.0.1), minimap2 2.30 for FLAIR's internal alignment (lab 2.31). isoseq collapse is not run (the request names StringTie and FLAIR).

## 2. Truth, strata, metrics (all from `reads.transcripts.tsv`; computed by `arms_score.py`)

Truth = every simulated transcript of every gene in the windows, canonical chains (introns >= 50 bp and canonical). A gene is ENTANGLED (stratum **E**) iff its exon union shares >= 100 bp with another gene's exon union on the same strand (the rule that defined the stratum for the 25 + 16 family copies); the other genes with >= 1 multi-exon chain are stratum **N**. Genes without a multi-exon chain are in neither.
Per arm and stratum:
- **M1 chains**: truth chains (chrom, strand, ordered intron list) equal to the chain of >= 1 arm transcript; also genes COMPLETE (every chain recovered).
- **M2 artifacts**: arm multi-exon transcripts whose chain is no truth chain: fragment (contiguous part of a truth chain), fusion (exons overlap >= 2 truth genes on the strand), other. Reported as a count and per recovered chain.
- **M3 resolution**: a truth gene is RESOLVED iff some arm gene_id carries an exact chain of it and an exact chain of no other truth gene. `merged_ids` = arm gene_ids carrying exact chains of >= 2 truth genes.
All four runs are pooled (per family and together); replicates are never averaged away: the per-run rows are kept.

## 3. Predictions (falsifiable, fixed now)

- **P1** Both S and F recover a smaller share of chains at entangled genes (E) than at the others (N).
- **P2** D_asm has fewer artifact transcripts per recovered chain than S and than F.
- **P3 (parity)** D_asm recovers at least as many chains at entangled genes as S and as F (pooled). If this fails, the gap to the tools is at the assembler and is the first thing to repair.
- **P4 (units)** U_asm recovers at least as many chains as D_asm, and U_fam resolves strictly more entangled genes than D_fam, with no loss at stratum N (chains and resolved).

## 4. Definition of "better", fixed before any improvement arm is built

An arm is BETTER than the tools at overlapping genes iff, pooled over the four runs, ALL of: (a) its chains recovered at stratum E >= the larger of S and F; (b) its artifact transcripts at least as few as the smaller of S and F; (c) its resolved genes at E >= the larger of S and F AND strictly above D_fam; (d) no loss at stratum N relative to D_asm / D_fam (chains, complete genes, resolved). Family-level consequences (E2-E4 per copy, `bench/ideal_expression`) are reported beside, not part of this bar.
A BETTER arm earns a registered real-data check (A119b chr16 and chr17, gorilla OR6737) with the copy-recovery instruments; nothing is claimed for real reads from this study.
If no arm is BETTER, the report lists where each loses (which of a-d) and the next lever is registered separately.

## 5. Limits declared in advance

Simulation with the annotation as truth; ideal reads (no truncation, no readthrough, no depth limit); human only; two replicates of one design; the FLAIR and StringTie versions differ from the lab's; U is an existing opt-in lever, not a new design; "resolved" is tool-agnostic (it asks for an own gene_id) and for tools whose gene_id is per transcript it equals "any chain recovered".

## 6. What was run before the freeze

One truth-vs-truth smoke test of `arms_score.py` (the NPIP replicate-1 truth transcripts scored as an arm, no tool output): 2,066 transcripts, 1,899 multi-exon; stratum E 86 genes with 568 chains, stratum N 264 genes with 1,295 chains; every chain recovered, 0 artifacts, 0 merged ids, every gene resolved. Nothing else of S, F or U existed, and D had not been scored by this script.
