# Active working set — what is actually in use (2026-09-19)

Companion to `docs/CLEANUP_CANDIDATES.md`, which marks what is *probably dead*. This file marks what is
**live**, so the audit stops treating it as a candidate: `tools/audit_cleanup_candidates.py` anchors any
file named by a top-level `docs/*.md`, so everything listed here is protected by being listed here.

Re-derive the tiers with `git log --name-only`, the ledger, and the audit TSV; re-run the audit after
editing this file.

## Tier 1 — the live pipeline

What the current §6q7 workflow actually invokes. Touch these with care; they are load-bearing.

### Rust entry points (`src/bin/`, PROTECTED by the audit)

| binary | role |
|---|---|
| `copy_assign` | ⭐ the main pipeline: loci → isoform assembly → GTF; `--assemble-only` is the assembler product, `--assembly-polish` the §6p8-§6q6 filters; O2 assignment lives here too |
| `mcl_families` | the DNA-level family definition (`--min-exonic-bp`, `--min-shared-exon-frac`) |
| `family_define`, `mcl_refine`, `gw_family_catalog` | family catalog construction and refinement |
| `gamma_refine`, `parcn`, `filter_bam_by_as`, `index_bam`, `bam_header` | supporting steps |
| `asj`, `asj_verify` | ⚠ ASJ is a DROPPED objective — kept for provenance, not in the live path |
| `debug_poa` | diagnostic only |

### Scripts

| script | role | last used |
|---|---|---|
| `bench/assembly_polish.py` | Python mirror of the Rust `--assembly-polish` passes; byte-identical parity oracle | §6q6, 2026-09-19 |
| `bench/ism_collapse.py` | standalone ISM collapse used for the chr20 precision levers | §6p7 |
| `bench/readthrough_secondary_filter.py` | secondary-dominated readthrough flagging (opt-in) | §6n9/§6o0 |
| `tools/refseq_gff_to_gtf.py` | RefSeq GFF3 → gffread-style GTF; validated at 4,574 = 4,574 vs `chr20_ref.gtf`. **Needed because `gffread` is not installed on this machine.** | §6p9-§6q7 |
| `tools/audit_cleanup_candidates.py` | this audit; read-only, re-runnable | §6q8 |
| `bench/layer_order/lattice_*.py` | the §0★★★ nested edge-test lattice (levels, edges, filtration, truth, report) | §6p1-§6p5 |

### Off-repo, but part of the live path

- `/mnt/linuxdisk/home/juanfraitu/lattice_rules/engine.py` — `L3_CUT`/`L4_CUT`, `l4_refine()` (§6p1/§6p5).
- `/mnt/linuxdisk/tmp/flair_shims/` — shims that make FLAIR 3.0.0 runnable here (§6q6); a second, broken
  FLAIR install on `linuxdisk` shadows the working one.
- Benchmark substrates under `/mnt/linuxdisk/home/juanfraitu/bakeoff/` — `human_chr{20,11,7,14,5,9}`
  (shallow `human_testis.t2t.bam`), `a119b_chr20` and `ggo_NC_073244.2` (the lab's deep libraries).
  ⚠ Register row 867: the two human libraries are **not** interchangeable.

## Tier 2 — reproduces a recorded result

**152 scripts are named directly in `docs/o1_ledger.md`, `docs/NEGATIVE_RESULTS_REGISTER.md` or a
`docs/PREREG_*.md`.** They are the provenance of published numbers: not live, but deleting one makes a
ledger claim unreproducible. Archive, never delete. Enumerate them with:

```sh
grep -ohE '(bench|tools|scripts|analysis)/[A-Za-z0-9_./-]+\.(py|sh)' \
     docs/o1_ledger.md docs/NEGATIVE_RESULTS_REGISTER.md docs/PREREG_*.md | sort -u
```

Five are already classed as supersedable and are the archive-first candidates within this tier:
`bench/denovo_shared_def.py`, `bench/o3_flag_pass.py`, `bench/vg_repeat_catalog.py` (SUPERSEDED-PORTED),
`bench/gw_rebuild.sh`, `bench/gw_rebuild_v2.sh` (SUPERSEDED-CITED).

## Tier 3 — anchored only transitively (the review pool)

Of 973 scripts in the repo, 575 are KEEP-CITED but only 152 are named directly; **the other 423 are
anchored only through another un-verified file.** ⚠ The 2026-09-16 round-1 verification found the
KEEP-CITED control itself was **4/12 dead** against a ≤25% bar, so transitive anchoring over-keeps. Tier 3
is where a purge should look after the candidate classes are dealt with — but only with per-file checks,
not in bulk.

## Numbers behind the tiers

| | files | note |
|---|---|---|
| scripts in repo | 973 | `.py` + `.sh` under bench/tools/scripts/analysis |
| named directly in ledger/register/prereg | 152 | Tier 2 |
| KEEP-CITED (anchored, mostly transitively) | 575 | 423 of them Tier 3 |
| candidates | 397 | see `docs/CLEANUP_CANDIDATES.md` |
| all files audited | 2,926 | 1,479 candidates, 343 MB |
