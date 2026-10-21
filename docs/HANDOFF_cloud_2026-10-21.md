# Handoff — local WSL machine → cloud (2026-10-21)

## Branch

Clone the repository and work on **`o3-fast`**: `git clone git@github.com:juanfraitu1/Rustle.git && git checkout o3-fast`.
- It contains `main` (`github/main` is an ancestor) plus 283 commits: the O2 / O3 work of 10-16 to 10-21.
- It is not merged into `main`.
- Registers 1347-1367 and their pre-registrations are on it; the latest summary is
  [`MULTIMAPPER_PREMISE_2026-10-21.md`](MULTIMAPPER_PREMISE_2026-10-21.md).

## Environment

- **Tools:** `pixi install` (see `pixi.toml`) gives minimap2 2.30, samtools 1.22.1, pysam, numpy and scipy.
  - `pixi install -e bench` adds gffcompare, stringtie and liftoff.
  - `pixi install -e isocon` adds IsoCon 0.3.3, with its own Python 3.8.
- **FLAIR 3.0.0** is not in pixi. It is not needed to continue: its outputs already exist (below).
- **Rust:** `pixi run build` (`CARGO_TARGET_DIR` defaults to `$HOME/rustle_target`; always `--release`; capture the log to a file).
- **Python tests:**
  - `python3 -m pytest -q bench/o3_graph` (481 pass at 3290e8a7);
  - `bench/rna_allele/test_npip_o3.py` (8);
  - `bench/tool_difference/test_s1_o2.py` (10);
  - `bench/o1_sd/test_sd_pairs.py` (5).
- **Paths:** scripts carry the local paths below (`/mnt/linuxdisk/...`, `/mnt/c/Users/jfris/Desktop/...`). On the cloud, either
  recreate those paths with symlinks, or pass the new paths where the scripts take arguments. The few hard-coded constants are:
  - `bench/rna_allele/npip_o3.py`: `BAM`, `GWFAM`, `TRUTH`, `GGO`;
  - `bench/o1_sd/sd_pairs.py`: catalogs and SD tracks;
  - `bench/o3_graph/donorsim.py` / `o2idealsim.py`: none (arguments only).

## State at handoff

| register | result | document |
|---|---|---|
| 1360, 1361 | the sample-aware corrector on NPIP / LRPAP1 tied reads in KB3781: UNDERPOWERED (final key) | `PREREG_o2_corrector_families_2026-10-20.md` |
| 1362 | the same corrector on ideal-expression reads from both haplotypes: FAIL (reads carry no positional evidence) | `PREREG_o2_corrector_idealsim_2026-10-21.md` |
| 1363 | reference-free O3 chain on KB3781 NPIP: finds the maternal-only units; FAIL on read recovery | `PREREG_o3_npip_maternal_copy_2026-10-21.md` |
| 1364, 1365 | candidates added to the reference: reads return (33 / 39), 0 false moves | `PREREG_o3_npip_candidates_realign_2026-10-21.md` |
| 1366 | premise check on existing StringTie / FLAIR runs (descriptive) | `MULTIMAPPER_PREMISE_2026-10-21.md` |
| 1367 | multi-copy gains are not shown to be real copies; the multi-mapping loci are SDs | `PREREG_multicopy_gains_real_2026-10-21.md` |
| 1368 | O1 family copies vs external SD calls: bar FAILS (72.1% / 82.4%); gene-like copies 93.8% / 92.0% paired; shortfall = transposon / satellite families and self-overlap families (an O1 defect) | `PREREG_o1_families_sd_pairs_2026-10-21.md` |

## Open items

1. **O1 self-overlap families** (register 1368): fixed as the opt-in `mcl_families --reject-self-hits`
   (`docs/PREREG_o1_self_hits_2026-10-21.md`).
   - The pre-registered v1 FAILED: gorilla B2 97.4%.
   - v2 (Amendment 1) keeps same-strand exon-sharing self-hits and PASSES all three catalogs. It is a development result,
     and its audit was pending at push time.
   - **Next:** confirm v2 on a fresh catalog (another sample's families run, pre-registered) before proposing it as the
     default. Also add a CAT-annotation families arm on CHM13 as the reference point.
2. **Classify the O1 copies that no SD pairs:** a quarter of the human high-identity copies, an eighth in gorilla. Check
   whether they are retrocopies (a single-exon copy with a multi-exon partner), groupings on short shared segments, tandem arrays,
   or repeat-driven. Pre-register first.
3. **Rerun the human Figure 3 table on CAT/Liftoff v2.0.** Register 1366's human numbers are scored on NCBI RefSeq.
4. **Generalize O3 beyond one locus:** a held-out reference-absent copy, then candidates plus realignment as in registers
   1363-1365.
5. **Older open threads:** `docs/PENDING_2026-10-11.md`, `docs/PENDING_2026-10-04.md`.

## Working rules this project follows

- **Pre-register first.** Write the decision rule (`docs/PREREG_*.md`, committed) before looking at results; change it only by
  a dated amendment.
- **Audit before reporting.** A fresh auditor (the `metric-auditor` agent) checks each verdict on a frozen, read-only copy with
  SHA1SUMS. Record the outcome in `docs/NEGATIVE_RESULTS_REGISTER.md`, and check that register before proposing anything.
- **Data invariants:**
  - never pool human and gorilla numbers;
  - human annotation = CAT/Liftoff v2.0;
  - CHM13 v1.0 ≠ v2.0;
  - primaries via `-F 2308`;
  - SAM is 1-based, exon tables 0-based;
  - gorilla SDs come from SEDEF, not our own `sd-scan` (register 1278: its stage-2 identity cut is a no-op).
- **O2 scope:** O2 = AS-tied multimappers only (a rival within 10 AS), allele-blind. A donor's own assembly may serve only as the
  answer key.
- **Families:** the multigene family (default `--rule multigene`) is the unit.
- **Push** only when asked.
- **Local-machine rules, not needed on a normal cloud machine:** the WSL crash rules (foreground, one heavy job at a time, no
  `pkill -f`) and the sparse-vhdx disk warnings.
- **Agent memory is not in this public repository.** To carry the local agent's notes over, copy
  `~/.claude/projects/-mnt-c-Users-jfris-Desktop/memory/` (2.8 MB) privately. The repository's `docs/MEMORY_DIGEST.md` holds the
  older compacted parts.

## Data to copy (local path → size → use)

### A. Needed for almost everything (~14 GB)

| file | size | use |
|---|---|---|
| `/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta` (+ `.fai`) | 3.4 GB | gorilla reference (GGO T2T, NC_ names) |
| `/mnt/linuxdisk/home/juanfraitu/gorilla_haps/{mat,pat}.fa` (+ `.fai`) | 3.4 + 3.2 GB | KB3781's two haplotypes: the answer keys |
| `/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa` (+ `.fai`) | 3.0 GB | human reference (human work only) |
| `/mnt/linuxdisk/tmp/rna_allele/{chrmap.tsv,mat.len.tsv,pat.len.tsv,fibro_cat.copies.tsv}` | < 1 MB | which haplotype contig equals the reference; the KB3781 catalog |
| `/mnt/linuxdisk/tmp/rustle_figures/runs/*/*.fam.copies.tsv` | < 10 MB | O1 catalogs (human_A119b, gorilla_KB3781, ...) |
| `/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_sedef_final.bed` | 103 MB | gorilla SDs (SEDEF) |
| `/mnt/c/Users/jfris/Desktop/CHM13_SDs.SD90.bed`, `CHM13_SDs.SD90.nosat.bed` | 30 + 19 MB | human SDs |
| `/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff` | 662 MB (gzip it) | gorilla RefSeq annotation |

### B. Reads: copy by task, or a subset

| file | size | advice |
|---|---|---|
| `/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam` (KB3781 fibroblasts, all records) | 22 GB | copy a **family-region subset** (recipe below) unless whole-genome work is planned |
| `/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_mm.bam` (gorilla OR6737 testis; the assembly benchmark) | 11 GB | only for the Figure 3 / StringTie-FLAIR work; subset possible |
| `/mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam` (human testis) | 90 GB | **do not copy**; subset to family regions if human work is planned |
| `/mnt/linuxdisk/tmp/rustle_figures/runs/{gorilla_OR6737,human_A119b}/*.molecules.tsv` | 260 MB / 1.4 GB | per-read best / second AS (Figure 3, register 1366-style analyses) |
| `/mnt/c/Users/jfris/Desktop/benchmark_collapse/` | 13 GB | copy only the GTFs (`stringtie_*/`, `flair_*/*.isoforms.gtf`, `assembler_*/`, `ref/`); skip FLAIR's BED / FASTA (7 GB) |
| `/mnt/linuxdisk/tmp/rustle_figures/assembly/*/eval_*/` (reference GTFs, tool GTFs, tmaps) | 7.1 GB | only for rescoring Figure 3; the tables in `figures/data/` are in the repo |

**Family-region subset recipe.** It keeps every record of every read that touches a catalog copy (±50 kb). Secondaries elsewhere
in the genome stay, so the O2 analyses still work:

```bash
awk -F'\t' 'NR>1 {s=$5-50000; if (s<0) s=0; print $4"\t"s"\t"$6+50000}' fibro_cat.copies.tsv | sort -k1,1 -k2,2n > fam.bed
samtools view -M -L fam.bed -F 4 BAM | cut -f1 | sort -u > fam.names
samtools view -@ 4 -N fam.names -b -o subset.bam BAM && samtools index subset.bam
```

### C. Do not copy (rebuild or regenerate)

- **Splice indexes** (`*.splice.mmi`, 13 GB each): rebuild with `minimap2 -x splice -d X.mmi X.fa`. That took 3.5 min and 23 GB RAM
  here; on a smaller machine, align against the FASTA directly.
- **Work directories:**
  - `/mnt/linuxdisk/tmp/o3_npip/` (17 GB, mostly an index copy);
  - `o3_graph/idealsim/` (2.5 GB);
  - `multicopy_gains/`, `o1_sd/`.
  
  Their results are in the documents above. Each has a frozen `audit_frozen_*` copy if a number needs rechecking.
