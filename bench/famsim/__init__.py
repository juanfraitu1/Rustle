"""famsim — controlled gene-family simulations that prove the condition they test.

Spec: docs/FAMSIM_DESIGN_2026-10-02.md. Usage: bench/FAMSIM.md. Run as `python3 bench/famsim <command>`.

Modules
  model       GeneModel: a copy's genomic sequence (transcript orientation) + its exons; the RNA chain derives from it
  ops         the mutation operators (snp, indel, exon_delete, splice_kill, exon_insert, exon_shuffle, invert, truncate,
              convert, intron_resize), each recording its realised coordinates
  template    the gene A from an annotation + genome (any species), a synthetic gene, or a saved template file; decoys
  chromosome  the artificial genome: background + planted copies -> FASTA, truth GTF/GFF3, copies table, manifest
  reads       IsoSeq-like reads per copy/isoform with the truth in the names (sim.simulate_reads + jitter)
  verify      re-derives the condition from the products alone -> verify.tsv (PASS/FAIL per claim)
  pipeline    minimap2 + the Rustle stages (assemble, de novo families, guided families, assign, flag)
  evaluate    scores every stage against the truth -> score.tsv, summary.md
  scenarios   the built-in ladder of scenarios
Only the standard library and pysam are imported; `bench/sim.py` and `bench/lib.py` are reused for the read model and
the scorers.
"""
import os
import sys

BENCH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BENCH not in sys.path:
    sys.path.insert(0, BENCH)
