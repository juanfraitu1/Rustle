# PREREG — `missing_copy_flag` positive control, re-measured with stable seeds (2026-10-14)

**Written before any scored run.** One plumbing run (3 genes, d = 0.02, seed 999, below) was made to check the commands
on this machine; seed 999 is excluded from every number here.

## Why

`REPRODUCE.md` (missing-copy section) quotes the positive control as "40/40 at 2% divergence, all confirmed", with
"re-measure before quoting": it was measured with the out-of-repo `/mnt/linuxdisk/tmp/gw22/o3/simB.py`, whose mutations
were seeded by Python's per-process `hash()` (B2, `bench/sim.py` docstring), so no read of it can be regenerated. The
2026-09-23 outcome table (`docs/archive/2026-09/PREREG_o3_rna_only_2026-09-23.md`, OUTCOME, positive control):

| divergence | fired / 40 | `reference_absent_candidate` | confirmed (`--confirm` the mutated copies) |
|---|---|---|---|
| 2% | 40 / 40 | 40 | 40 / 40 |
| 1% | 15 / 40 | 15 | 15 / 15 |

This pre-registration re-measures both rows with `sim.py missing-copy genomic` (stable seeds) and the current
`missing_copy_flag`, over five seeds instead of one, and adds a negative control the original did not have.

## Substrate (every input from NCBI, on the 2026-10-14 laptop)

- Genome: CHM13 v2.0 chr20 = NC_060944.1 (efetch, renamed `chr20`; 66,210,255 bp = the T2T length).
- Annotation: RefSeq GCF_009914755.1 (annotation date 2025-08-01), chr20 rows; `utilities gff-to-gtf ... chr20` gives
  4,574 transcripts = the record's `chr20_ref.gtf` count (`docs/DATA.md`).
- Tools: minimap2 2.30-r1287, samtools 1.22.1 (`pixi.toml`); `missing_copy_flag` built `--release` from the commit
  this file is committed in.

**Deviations from 2026-09-23, fixed now:**
1. The home-search index is chr20 only (`minimap2 -x splice:hq -d chm13_chr20.splice.mmi`), not whole-genome CHM13
   (laptop disk and memory). It can only hide an `unannotated_paralogue` call whose better home is on another
   chromosome. The consensus is a 1-2% mutant of the host locus itself, so a better home elsewhere would need a paralogue
   closer to the mutant than its own source; this is not expected to occur, and every `unannotated_paralogue` call is
   listed.
2. Gene IDs are RefSeq `gene-…` IDs (the GTF's `gene_id`); the 40 genes differ from 2026-09-23's in any case (unknown
   seed).
3. Five seeds per divergence (1, 2, 3, 4, 5) instead of one; pooled n = 200 planted copies per divergence.

## Commands (per divergence d ∈ {0.02, 0.01}, per seed s ∈ {1..5})

```bash
python -I bench/sim.py missing-copy genomic chr20_ref.gtf chm13_chr20.fa chr20 40 $d sim_d${d}_s${s} $s --threads 4
missing_copy_flag --bam sim_d${d}_s${s}.bam --fasta chm13_chr20.fa --loci chr20_ref.gtf \
    --index chm13_chr20.splice.mmi --gff chm13_chr20.RefSeq.gff --confirm copies=sim_d${d}_s${s}.copies.mmi \
    --threads 4 --out run_d${d}_s${s}
python bench/score.py missing-copy-control --sim sim_d${d}_s${s} --run run_d${d}_s${s}
```

**Negative control (N), d = 0.02 seeds only:** the `template|…` reads of `sim_d0.02_s${s}.fq` alone (the extra copy's
reads removed), mapped with the shipped `sim.MM2` command, then the same `missing_copy_flag` call (same `--confirm`
index) and the same scorer. Every flag default is used (`--m-min 0.10 --delta-min 0.01 --min-reads 10 --min-sub 3`).

## Unit and scoring

The unit is a **planted gene** (a line of `<sim>.genes.txt`), scored on the table row whose `locus` equals that gene ID.
Overlapping loci (nested genes, readthroughs) share the reads and fire with it (plumbing: 8 rows fired for 3 genes);
they are counted in a separate column, never in the denominators.
- **fired**: the row's `status` starts with `fired`. A gene with no row (fewer than `--min-reads` reads) is not fired.
- **RAC**: fired and `verdict == reference_absent_candidate`.
- **confirmed**: RAC and the `copies` confirmation column is 1 (the confirmation rule unchanged:
  `conf_identity >= 0.99` and `conf_identity - host_identity >= delta/2`).

## Bars (pooled over the five seeds)

| id | divergence | measure | bar | 2026-09-23 value |
|---|---|---|---|---|
| **P1a** | 2% | fired / 200 | ≥ 190 (0.95) | 40 / 40 |
| **P1b** | 2% | RAC / fired | ≥ 0.95 | 40 / 40 |
| **P1c** | 2% | confirmed / RAC | ≥ 0.95 | 40 / 40 |
| **P2a** | 1% | fired / 200 | in [50, 100] (0.25-0.50) | 15 / 40 = 0.375 |
| **P2b** | 1% | confirmed / RAC | ≥ 0.95 | 15 / 15 |
| **N** | template reads only (2% seeds' genes) | fired / 200 | ≤ 2 (0.01) | not measured |

Per-seed values are reported next to the pooled ones; the bars apply to the pooled values only.

## Reading, fixed now

- **P1a-c and P2a-b and N pass** → the control reproduces: `REPRODUCE.md` replaces "40/40 … re-measure before quoting" by
  the pooled numbers and this file's commands.
- **P1 fails** → the 40/40 does not reproduce on stable seeds; `REPRODUCE.md` and `docs/THESIS_OBJECTIVES.md` quote the
  new numbers, and the failing step (mixture, consistency, verdict or confirmation) is named from the scorer's funnel.
  No threshold of `missing_copy_flag` is changed in response (that would need its own pre-registration).
- **P2a outside its interval** → reported as a change in the 1% sensitivity edge, with the direction; not a failure of
  the detector by itself (the 2026-09-23 edge rested on one seed of 40).
- **N fails** → the detector fires on error-only reads at a simulated locus; reported as a false-positive rate on
  simulated negatives and investigated before any O3 number from `missing_copy_flag` is quoted.

## Plumbing run (seed 999, excluded)

3 genes at d = 0.02 (gene-CST4, gene-STK35, gene-SYS1-DBNDD2): 80 reads, simulation 5 s, flag about 2 s; the three
planted genes' rows and five overlapping rows fired. This checks the commands only and is no part of any bar.

# OUTCOME (2026-10-14) — every bar passes; the control reproduces on stable seeds

`missing_copy_flag` built `--release` at 49e62be3 (binary sha1 159643eb), `bash bench/o3_posctl.sh REF OUT` under
`pixi run`; md5 `chm13_chr20.fa` 1d0b5621077818120e7d40dae01ad4d7, `chr20_ref.gtf` 84145d8099109ddb1b58afc059017bc2.
Run products: `~/rustle_data/o3_posctl/scored/` on the laptop (one `*.per_gene.tsv` per run).

| seed | 2% fired | RAC | confirmed | N fired | 1% fired | RAC | confirmed |
|---|---|---|---|---|---|---|---|
| 1 | 39 | 39 | 39 | 0 | 13 | 13 | 12 |
| 2 | 40 | 40 | 39 | 0 | 21 | 21 | 20 |
| 3 | 40 | 40 | 40 | 0 | 21 | 21 | 21 |
| 4 | 40 | 39 | 39 | 0 | 19 | 18 | 18 |
| 5 | 40 | 40 | 40 | 0 | 17 | 17 | 17 |
| **pooled** | **199 / 200** | **198** | **197** | **0 / 200** | **91 / 200** | **90** | **88** |

| id | measure | value | bar | verdict |
|---|---|---|---|---|
| P1a | 2% fired / 200 | 199 (0.995) | ≥ 190 | ✅ PASS |
| P1b | 2% RAC / fired | 198 / 199 (0.995) | ≥ 0.95 | ✅ PASS |
| P1c | 2% confirmed / RAC | 197 / 198 (0.995) | ≥ 0.95 | ✅ PASS |
| P2a | 1% fired / 200 | 91 (0.455) | 50-100 | ✅ PASS (above the 2026-09-23 0.375, inside the interval) |
| P2b | 1% confirmed / RAC | 88 / 90 (0.978) | ≥ 0.95 | ✅ PASS |
| N | template-only fired / 200 | 0 | ≤ 2 | ✅ PASS (all 200 `no_mixture`) |

**Reading (as fixed):** the control reproduces. The 2026-09-23 "40/40 at 2%" becomes **199/200 fired, 198
`reference_absent_candidate`, 197 confirmed (five seeds)**; the 1% edge is **91/200 fired, 88/90 confirmed**; the detector
fired on none of the 200 planted genes when the extra copy's reads were removed. No `unannotated_paralogue` call occurred,
so deviation 1 (chr20-only home index) hid nothing on this substrate.

**The misses, named (reported, nothing re-tuned):** every one is at a lncRNA / uncharacterized `LOC` gene.
- *Mixture* (1, at 2%): `gene-LOC105372546` (seed 1), `no_mixture`.
- *Consistency* (2): `gene-LOC107985409` (2%, seed 4: n_psv 2, shared 0.26) and `gene-LOC105372659` (1%, seed 4: n_psv 4,
  shared 0.37), `scattered`.
- *Confirmation* (3): `gene-LINC00659` (2%, seed 2: conf 0.9943 vs host 0.9915, margin below delta/2 = 0.0067),
  `gene-LOC105372727` (1%, seed 1: conf 0.9886 < 0.99, host 1.0000) and `gene-LOC107985409` (1%, seed 2: conf 0.9936 vs
  host 0.9893). A host identity of 1.0000 at delta 0.02 says the patched consensus covers almost none of the PSVs (it is
  built over the template read's span only); the confirmation rule is then right to withhold.
