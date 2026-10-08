# Locus Representative Rule Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** add `mcl_families --representative most-junctions` (+ driver knob), then run the pre-registered dev + held-out comparison and decide the default.

**Architecture:** one flag in `mcl_families --from-gtf`'s locus builder (`GtfLocus` representative choice), threaded through `tools/rustle_pipeline.sh families` as `RUSTLE_REPRESENTATIVE`; a bench runner that re-runs the families stage per contig under both rules on the Figure 7 de novo GTFs and scores with Figure 7's own scorers plus `bench/copy_support.py`.

**Spec:** `docs/archive/2026-10/PREREG_locus_representative_rule_2026-10-04.md` (binding: rules, substrates, readouts, decision).

## Global Constraints

- Build only with `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target_m2`, `--release`, cargo output to a file.
- Heavy commands foreground via `bash tools/rlock.sh heavy ...`, each call < 10 min; light via `rlock light`; no backgrounding, no `pkill -f`.
- Never touch `/mnt/c/Users/jfris/Desktop/Rustle`; big outputs under `/mnt/linuxdisk/tmp/rep_rule/`; the Figure caches under `/mnt/linuxdisk/tmp/rustle_figures/` are READ-ONLY.
- Registered rules are fixed; harness bugs are fixed, rules never. Commits end with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>` and `Claude-Session: https://claude.ai/code/session_01DAyQQ6R8drUxY5GsM5kNkb`; author from the canonical repo's git config; no push.

## Review Focus

1. `--representative` unset or `most-reads` must be byte-identical to today's outputs (`loci.gff3`, `copies.tsv`, `clusters.tsv`) — a fixture GTF run both ways.
2. The junction count ignores gaps < 50 bp and counts gaps between exons sorted by coordinate (not file order).
3. Tie order under `most-junctions`: most reads, then longer span, then last `transcript_id` — a unit test with ties at each level.
4. The driver refuses `RUSTLE_REPRESENTATIVE` with a binary that predates the flag (the `--help` probe pattern already used for `--min-cov-shorter`).
5. The runner must use the SAME binary for both arms and record its sha1 in every log.

---

### Task 1: the flag

**Files:** `src/bin/mcl_families.rs` (GtfLocus builder ~:2539-2625, Args, the `--from-gtf` doc), `tools/rustle_pipeline.sh` (families stage: env knob + help probe + header comment), `docs/MODULE_STATUS.md` / `README.md` one line each if they list the families flags.

- [ ] Add `--representative {most-reads,most-junctions}` (clap `ValueEnum`, default `most-reads`) to Args; document it in the `--from-gtf` doc comment and the `GtfLocus` doc.
- [ ] In the locus builder: compute per transcript `junctions = number of gaps >= 50 bp between consecutive coordinate-sorted exons`; `most-junctions` key = `(junctions, reads, span)` over the sorted transcript ids (last id wins ties, as today); `most-reads` key unchanged.
- [ ] Unit tests: (a) the existing test still passes; (b) a locus where a 1-junction 10-read transcript loses to a 5-junction 2-read transcript under `most-junctions` and wins under `most-reads`; (c) ties at junctions resolve by reads, then span, then id; (d) a 30-bp gap is not a junction.
- [ ] Driver: `RUSTLE_REPRESENTATIVE` → `--representative "$RUSTLE_REPRESENTATIVE"` on the families stage, with the help-probe refusal when the binary lacks the flag; header comment; byte-identical when unset.
- [ ] Build, run the mcl_families tests + the driver e2e fixture if one covers `families`, commit.

### Task 2: the runner and the registered comparison

**Files:** `bench/rep_rule/run.sh`, `bench/rep_rule/score.py`, `docs/archive/2026-10/LOCUS_REPRESENTATIVE_RULE_2026-10-04.md`, register rows.

- [ ] Per contig (dev: human chr16, gorilla NC_073244.2; held-out: human chr2, chr6, chr8, chr10; gorilla NC_073234.2): copy the Figure 7 `.denovo.gtf` into `/mnt/linuxdisk/tmp/rep_rule/<species>_<contig>/`, run the driver's `families` stage twice (R_M: env unset; R_J: `RUSTLE_REPRESENTATIVE=most-junctions`) with the Task 1 binary, `--bin` the m2 release dir, `/usr/bin/time -v`, logs with the binary sha1; one heavy call per (contig, arm).
- [ ] H1: score each arm's `fam.clusters.tsv` exactly as Figure 7 does — import `figures/_o1_recovery.py` (`score_arm` / `family_score` against Compara and Soto with the contig's GFF slice; the Liftoff pair sensitivity through its `copy_pairs` / `pair_families` with the existing fig8 read-support products; `npip_u2` on chr16). If a reference product is missing, say so; do not rebuild heavy products.
- [ ] H2: from each arm's `fam.copies.tsv`: median junctions per copy (`n_exon - 1`), fraction with >= 2, total exon bp; the number of loci whose representative differs between arms (join `loci.gff3` by gene id).
- [ ] H3: `bench/copy_support.py` per contig and arm: copies = the contig's annotated protein-coding genes (human CAT slim GFF `winloci_data/gencode_chm13/chm13v2.0_CAT_Liftoff.slim.gff3.gz`; gorilla `winloci_data/GGO_genomic.gff`) written as a copies table (cid, family=ALL, name, chrom, strand, terr_lo0, terr_hi, isoform_gene) + a truth GTF of their transcripts; BAMs `winloci_data/A119b.t2t.bam` / `winloci_data/GGO_mm.bam`; `--loci R_M=<loci.gff3>,<gtf> --loci R_J=...` (the same assembled GTF for both; the loci GFF3 differ). On chr16 also the 25 NPIP copies (the page's tables) with the locus-level reading beside. Extend `copy_support.py` only if a needed input form is missing; the rule is fixed.
- [ ] Decision per the prereg; write the doc (per-contig tables, the genes that change, run times), register rows, and — iff adopted — flip the default (`--representative most-junctions` default + driver default) in a SEPARATE commit with the e2e/byte-identity note; otherwise record "R_J stays opt-in".
