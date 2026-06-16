# Contiguous-core gate: real-pipeline measurement (Phase C)

Date: 2026-06-16. Measures the (robustness-fixed) contiguous-core family-merge
gate inside rustle's ACTUAL `--vg --vg-layer2` family formation on a real
chromosome, gate OFF vs ON.

- Gate env: `RUSTLE_VG_FAMILY_MIN_CORE_COVERAGE` (default 0.0 = OFF). Operating
  value tested: **0.13**.
- Binary rebuilt from `src/rustle/vg_family/family_graph.rs` (the prior-phase
  Rust fix is in `contiguous_core_coverage` / `merge_singletons_by_sequence`).
- Contig: **NC_073235.2** (gorilla autosome, 150.2 Mb, 188,134 reads), full
  contig. Inputs: `/home/juanfra/winloci_scratch/GGO.{bam,fasta}`, guide
  `/tmp/gw/st_NC_073235.2.gtf`. RAYON_NUM_THREADS=4, serial, OOM-safe.

## TL;DR

YES — the gate changes real-data family formation, in the **expected direction**.
At threshold 0.13 it drops exactly **5 domain-sharing exon merges** (contiguous
core 0.004–0.011 of the shorter exon) while leaving all **62 true-copy merges**
(core 0.31–1.00) intact. Downstream the change is small and benign: one real
locus (RSTL.647) gets its copy/novelty attribution cleaned; two others only get
cosmetic family-id renumbering. No transcripts gained or lost (2623 = 2623).

## Where the gate actually lives (architecture)

The gate is NOT at the cross-locus family-GROUPING level (`discover_family_groups*`,
k-mer Jaccard 0.20 over bundles — that decides which loci join a `FamilyGroup`,
and the gate never touches it). It is at the **cross-copy singleton-exon merge**
inside `merge_singletons_by_sequence` (family_graph.rs:~416-445), invoked from
`assemble_family_graph_from_copy_exons` (family_graph.rs:538-540) only when
`n_singletons >= 2 && total_exons <= 2000`. It governs which homologous exons
from DIFFERENT copies fuse into one shared family-graph node (ExonClass / PSV
column).

Consequence: families whose copies overlap in genomic POSITION merge via Stage-1
`cluster_by_position` and never see the gate. Only the dispersed-paralog
sequence-merge path (Stage-1b) is gated. The `[DIAG] family ... member_spans`
dump (gated by `RUSTLE_DIAG_TARGET`) reports family→loci membership but is
therefore the WRONG instrument for this gate; the right instrument is the
per-pair core-coverage at the merge loop.

## Instrumentation added (additive, default-off)

Added `RUSTLE_VG_CORE_GATE_TRACE` in `merge_singletons_by_sequence`: for every
Jaccard-passing cross-copy pair it prints
`[CORE_TRACE] cid_i= cid_j= jacc= core_cov= len_i= len_j= would_gate=`.
Computes `core_cov` for the trace even when the gate is off (so the full
distribution is visible); never changes merge decisions. Family-graph unit tests:
17 pass / 0 fail.

## Result 1 — core-coverage distribution of Jaccard-passing pairs (gate OFF)

67 Jaccard-passing cross-copy singleton pairs on the contig. Bimodal:

  core_cov < 0.05 : 5     <- domain-sharers (0.004, 0.008, 0.011, 0.011, 0.011)
  core_cov < 0.13 : 5     (no pairs between 0.011 and 0.314 — clean empty gap)
  core_cov >=0.31 : 62    true copies (mean of all 67 = 0.739)

The threshold 0.13 sits squarely in the empty bimodal gap.

The 5 low-core pairs (all in ONE family, copies cid=0 vs cid=2):

  jacc=0.386 core_cov=0.004 len_i=1140 len_j=2846
  jacc=0.636 core_cov=0.008 len_i=808  len_j=662
  jacc=0.707 core_cov=0.011 len_i=180  len_j=233   (x2)
  jacc=0.605 core_cov=0.011 len_i=180  len_j=233

These are textbook domain-sharers: enough shared 15-mer minimizers to clear the
Jaccard bar, but the longest contiguous homologous run covers <1.1% of the
shorter exon, and lengths are mismatched (e.g. 1140 vs 2846 bp). Not paralog
copies.

## Result 2 — gate ON @0.13 fires on exactly those 5 pairs

With `RUSTLE_VG_FAMILY_MIN_CORE_COVERAGE=0.13`, `would_gate=true` count = 5 —
precisely the 5 sub-0.013 pairs. The 62 true-copy merges are untouched. Direction
is exactly as designed: drop short-block domain-sharers, keep real copies.

## Result 3 — downstream output delta (OFF vs ON), full contig

Transcript count identical: 2623 vs 2623. GTF differs (md5 differs); the diff is
3 genes:

- **RSTL.447** (FAM_13 -> FAM_11) and **RSTL.508** (FAM_0 -> FAM_4): ONLY the
  `family_id` integer label changed. All coordinates, copy_id, abundances,
  family_verdict, copy counts byte-identical. Cosmetic renumbering caused by
  removing 5 merges upstream (family enumeration order shifts). Not a real
  formation change at these loci.

- **RSTL.647** (NC_073235.2:112,808,333-112,845,893, the domain-sharer family,
  cid0/cid2): REAL effect. Same 4 transcripts, same exon coordinates, but
  transcripts are re-ranked and the `rescue_class "strand_pure_minority"` +
  `copy_status "novel"` tags present with the gate OFF are DROPPED with it ON.
  Removing the spurious cross-copy merge that was contaminating this family's
  graph cleaned up its novel-copy attribution. This is the gate's intended
  precision behaviour on a real locus.

## Verdict

- Did the gate change family formation on real data? **Yes** — 5 merges removed
  on this contig.
- Expected direction (drop domain-sharers, keep copies)? **Yes, exactly** — all
  5 removed pairs are <0.013 contiguous core (domain-sharers); all 62 true-copy
  pairs (>=0.31) are preserved. Bimodal gap at 0.13 is clean.
- Output cost? Minimal/benign: 0 transcripts gained or lost; 1 locus
  (RSTL.647) cleaned; 2 cosmetic family-id renumbers.
- The gate is narrow on this contig because (a) only the dispersed-paralog
  Stage-1b merge path is gated (co-located families bypass it), and (b)
  genuine domain-sharers passing Jaccard are rare (5 of 67 pairs here, all in
  one family).

## Reproduce

  # rebuild
  cargo build --release
  # contig slice
  samtools faidx GGO.fasta NC_073235.2 > NC_073235.2.fasta; samtools faidx NC_073235.2.fasta
  samtools view -b GGO.bam NC_073235.2 > NC_073235.2.bam; samtools index NC_073235.2.bam
  # OFF (trace shows full distribution)
  RAYON_NUM_THREADS=4 RUSTLE_VG_CORE_GATE_TRACE=1 rustle --vg --vg-layer2 \
    --genome-fasta NC_073235.2.fasta -G st_NC_073235.2.gtf -L NC_073235.2.bam -o off.gtf
  # ON @0.13
  RAYON_NUM_THREADS=4 RUSTLE_VG_CORE_GATE_TRACE=1 RUSTLE_VG_FAMILY_MIN_CORE_COVERAGE=0.13 \
    rustle --vg --vg-layer2 --genome-fasta NC_073235.2.fasta -G st_NC_073235.2.gtf \
    -L NC_073235.2.bam -o on.gtf
  # compare: grep -c would_gate=true; diff off.gtf on.gtf

## What a genome-wide run would need

Per-contig serial (JOBS=1, the OOM protocol), 24 contigs, ~25-35 min total,
mirroring `bench/copy_recovery_eval/results_genomewide/gw_run.sh`. Run each
contig twice (gate OFF/ON) with `RUSTLE_VG_CORE_GATE_TRACE=1`; aggregate
`would_gate=true` counts and the core_cov histogram across contigs, and diff the
two GTF sets per contig (separating cosmetic family-id renumbers from real
structural/attribution changes, as done here). Expected scale: a handful of
domain-sharer merges per contig, concentrated in a few families.
