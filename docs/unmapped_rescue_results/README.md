# Backing files for docs/UNMAPPED_RESCUE_2026-10-08.md and the Amendments of docs/PREREG_unmapped_rescue_2026-10-08.md

- `*.log`: stdout of the real-bed analyses (`apply_trim_real.py`, `eval_clusters.py`, `consensus_support.py`, `lrpap1_gate_flag.py`, `o3_haplotype.py`, `o3_bed.py`, `discover.py`, `stress_part.py`, `run_partition_real.py`), first line = the command.
- `dna_verify.jim.{log,json}`, `candidates_vs_human.log`, `discovery_candidates.fa`: Amendments 28 to 32.
- `worlds/<world>/`: the result files (< 120 KB) of each synthetic world. The worlds themselves (genomes, reads, SAM) were deleted to free disk; they rebuild deterministically from the seed with
  `python3 bench/unmapped_rescue/synth_world.py build <dir> --seed <seed> [flag]`: `synth` 20261008 (default), `synth2` 20261009 (classes 0.005..0.08, 4 per class), `synth3` 20261010,
  `synth4` 20261011, `synth5` 20261012, `synth6` 20261013 `--lead-g`, `synth7` 20261014, `synth8` 20261015, `synth9` 20261016 (all three `--partition-world`), `synth10` 20261017,
  `synth11` 20261018 (`--chain-world`), `synth12` 20261301 `--lead-g`, `synth13` 20261401 `--partition-world`, `synth14` 20261402 `--spec-world`,
  `chainw_<seed>_<depth>` = `--chain-world --seed <seed> --per-copy <depth>`.
