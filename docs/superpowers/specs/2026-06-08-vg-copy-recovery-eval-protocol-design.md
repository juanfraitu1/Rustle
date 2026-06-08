# Formal evaluation protocol: rustle (VG) recovers paralog copies StringTie misses

**Date:** 2026-06-08
**Status:** Design approved, pending implementation plan
**Type:** Evaluation protocol (no phasing — copy recovery via secondary alignments only)

## Claim under test

On the same genome-wide long-read BAM, **rustle `--vg`** (which retains secondary/multi-mapping
alignments) recovers annotated paralog-family transcript copies that **StringTie2 misses**, and the
recoveries are **authentic** (backed by the copy's own primary-alignment evidence) rather than
secondary-alignment **phantoms**.

The protocol must produce a reproducible, advisor-defensible measurement of this claim, using
**SQANTI3** as the structural classifier and a **tool-independent authenticity guard** for the
paralog-specific over-counting failure mode.

## Why this protocol exists (the lessons it formalizes)

Prior evaluation was ad-hoc (`winloci_eval.sh` and friends) and carried three validity threats this
protocol eliminates:

1. **gffcompare / structural matchers reward fabrication on near-identical paralogs.** A transcript
   rustle assembles at a *silent* copy A's locus out of an expressed sister copy B's spilled-over
   secondary alignments sits at A's coordinates with the shared splice structure → SQANTI3 calls it
   **FSM to A** ("copy A recovered"), even though A is not expressed. This is a *locus-of-origin
   attribution* error, not a structural error (the structure is genuinely a real annotated
   transcript). It is specific to multi-copy families with multi-mapping reads — exactly rustle's
   target regime. Demonstrated cases: DAZ3 (3 strict / 38 tied, **0 primary reads**), TSPY (prim=0).
2. **Window fragility.** Per-locus extraction windows reverse outcomes (a 20 kb shift took RBMY from
   +1 to −5). Any windowed protocol inherits this.
3. **Non-determinism / LLM-in-the-loop.** Wrapping deterministic evaluation in agents mis-reported
   results (a screen over-reported 30 wins vs a true 20/4). The protocol is a single deterministic
   driver.

## Design decisions (locked)

- **Structural matcher:** SQANTI3 (`sqanti3_qc.py`). Recovery = a reference transcript with an
  **FSM** query isoform (an **ISM** set is tracked separately, never merged into the headline).
- **Authenticity guard:** tool-independent. A recovery is authentic iff the copy's locus has
  ≥ `k` reads whose **primary** alignment falls within the copy's exons (read straight from SAM
  flags; not a rustle-internal score). Default `k = 1`. The protocol reports how many FSM copies
  fail the guard (the phantom count).
- **Run scope:** genome-wide, **no windows**. Each tool runs once on the full BAM → one
  transcriptome GTF → one SQANTI3 run. The paralog universe is a post-hoc filter on reference
  transcripts. This eliminates window fragility entirely.
- **Paralog universe U:** tool-independent. RefSeq gene families with ≥2 copies, confirmed by
  pairwise sequence identity from self-alignment (pinned identity/length thresholds). Computed once
  from annotation + genome.
- **Baseline:** **StringTie2 `-L` only** (the named competitor, clean headline). Disclosed
  limitation: this confounds "uses secondaries" with "different assembler." An optional
  **rustle-primary-only** arm (off by default) can be enabled to attribute the win to secondaries.
- **rustle `--vg` config:** headline = the **intended product config (decisive gate ON)**; an
  optional **raw / no-gate ablation** is reported to show what the gate + guard buy.
- **Dataset:** GGO genome-wide long-read BAM (`$SCRATCH/GGO.bam`, 1.6 G, 380,369 secondary
  alignments), genome FASTA (`$SCRATCH/GGO.fasta`), RefSeq annotation (`$SCRATCH/GGO_tx.gff`,
  403 M). `$SCRATCH = /home/juanfra/winloci_scratch`.

## Architecture — six isolated units + one deterministic driver

Each unit has one responsibility and a file-level interface (a written artifact consumed by the
next unit), so units are independently testable.

```
GGO.bam (secondaries retained) + GGO.fasta + GGO_tx.gff(→GTF)
   │
   ├─[U1] transcriptome generation (once each, full BAM, no windows)
   │        ├─ StringTie2 -L          → stringtie.gtf
   │        ├─ rustle --vg [gate ON]  → rustle_vg.gtf            (headline)
   │        ├─ rustle --vg [no gate]  → rustle_vg_raw.gtf        (optional ablation)
   │        └─ rustle -L              → rustle_primary.gtf       (optional attribution arm)
   │
   ├─[U2] SQANTI3 QC per arm (vs RefSeq + genome) → <arm>_classification.txt
   │
   ├─[U3] paralog universe U (tool-independent, computed once) → universe.tsv
   │        RefSeq families ≥2 copies, identity-confirmed; transcript→family map
   │
   ├─[U4] recovery-set extraction per arm → <arm>_recovery.tsv
   │        FSM (and ISM) reference transcripts via associated_transcript, ∩ U
   │
   ├─[U5] authenticity guard (tool-independent, on rustle-VG recoveries) → authenticity.tsv
   │        primary-alignment read support per recovered copy locus; authentic/phantom
   │
   └─[U6] head-to-head + report → headtohead.tsv + REPORT.md
            rustle-VG-only recoveries in U, FSM/ISM × authentic/phantom;
            genome-wide SQANTI3 category distributions per arm; per-family table;
            worked examples (RABL2, RBMY)
```

### U1 — Transcriptome generation
- StringTie2: `stringtie -L -G <none for de novo> ...` on the full BAM. (De novo assembly; the
  reference is used only by SQANTI3, not given to the assemblers — both tools run reference-free so
  the comparison is of *assembly* not *guided quantification*.)
- rustle: `rustle --vg [pinned flags incl. decisive gate]` on the full BAM. Pinned config recorded
  with a hash. Optional arms (`--vg` no-gate; `-L` primary-only) gated behind config switches.
- **Determinism:** record exact command line, tool version, and config hash per arm into
  `provenance.json`. rustle's GTF line-order has a documented rayon nondeterminism — the protocol
  must sort each GTF canonically before SQANTI3 so downstream is stable.

### U2 — SQANTI3 QC
- One `sqanti3_qc.py <arm>.gtf <RefSeq.gtf> <genome.fasta>` per arm. No CAGE/polyA/short-read SJ
  inputs in v1 (note as a possible future enrichment).
- Output of interest: the `*_classification.txt` columns `isoform`, `structural_category`,
  `associated_gene`, `associated_transcript`, plus FL-count and junction columns for cross-checks.
- SQANTI3 is **not currently installed**; the plan must install it in a dedicated conda env (it
  needs cDNA_Cupcake / specific dependency pins) and record the version.

### U3 — Paralog universe U (tool-independent)
- From the RefSeq annotation: group transcripts by gene; identify genes belonging to multi-copy
  families. Confirm paralogy by extracting per-gene representative sequences and self-aligning
  (minimap2 or equivalent), keeping pairs with identity ≥ `MIN_IDENTITY` over ≥ `MIN_COV_FRAC` of
  length. Output `universe.tsv`: `transcript_id, gene_id, family_id, n_family_copies`.
- Pinned thresholds; the criterion never references any tool's output.

### U4 — Recovery-set extraction
- Parse each arm's SQANTI3 classification. For each reference transcript, record whether the arm has
  ≥1 query isoform with `structural_category == FSM` (FSM set) and separately `== ISM` (ISM set),
  keyed on `associated_transcript`. Intersect with U. Output `<arm>_recovery.tsv`:
  `ref_transcript, family_id, fsm(bool), ism(bool)`.

### U5 — Authenticity guard
- For each rustle-VG FSM/ISM recovery of a copy in U: resolve the copy's exon intervals (from the
  RefSeq annotation of the `associated_transcript`), then count reads in the BAM over those
  intervals whose alignment is **primary** (SAM flag: not 0x100 secondary, not 0x800 supplementary)
  — `samtools view` with flag filters. Authentic iff count ≥ `k` (default 1). Output
  `authenticity.tsv`: `ref_transcript, family_id, primary_support, authentic(bool)`.
- This is intentionally independent of rustle's `compute_copy_ownership` so the headline claim does
  not rest on a rustle-internal metric. (rustle's certificate may be reported as a secondary,
  mechanism-explaining column, clearly labeled, but is not the adjudicator.)

### U6 — Head-to-head + report
- **Headline set:** `recovered_by_rustleVG ∩ U` minus `recovered_by_StringTie2 ∩ U`, split by
  FSM/ISM and by authentic/phantom. Net win = authentic FSM gains; phantom = disclosed cost.
- **Global context:** genome-wide SQANTI3 category counts/fractions per arm (FSM, ISM, NIC, NNC,
  genic, antisense, fusion, intergenic) — shows rustle-VG does not wreck overall precision.
- **Per-family table:** for each family in U, copies recovered by each arm + authenticity.
- **Worked examples:** RABL2 and RBMY rows surfaced explicitly.
- Output `headtohead.tsv` + a human-readable `REPORT.md`.

## Determinism & reproducibility

- One driver script `run_protocol.sh` (+ small Python helpers); **no agent/LLM in the deterministic
  path**.
- All thresholds and tool flags in one `config.sh`/`config.json`; emit a config hash.
- Pin: rustle commit SHA, StringTie2 version (3.0.1 present), SQANTI3 version, minimap2 version,
  RefSeq annotation file hash.
- Canonical-sort every GTF before SQANTI3. Re-running the protocol on the same inputs reproduces
  `headtohead.tsv` byte-for-byte (modulo any documented SQANTI3 nondeterminism, which must be
  checked and noted).

## Protocol self-validation (oracles)

The protocol's own correctness is checked by known cases before any genome-wide claim is trusted:

- **Known positive:** RABL2 must appear as an authentic rustle-VG-only recovery (your worked example
  VG 9 vs ST 4).
- **Known negative:** DAZ3 must be classified PHANTOM by the authenticity guard (prim=0).
- **Unit tests:** U3 (universe set logic), U4 (recovery extraction from a toy classification.txt),
  U5 (primary-support counting from a toy BAM) tested on small synthetic fixtures with hand-computed
  expected outputs.

## Outputs

A results directory containing: `provenance.json`, the per-arm GTFs, per-arm SQANTI3
classifications, `universe.tsv`, per-arm `*_recovery.tsv`, `authenticity.tsv`, `headtohead.tsv`, and
`REPORT.md`.

## Non-goals (v1)

- No phasing / haplotype resolution (separate paused project).
- No CAGE/polyA/short-read SJ enrichment of SQANTI3 (possible later).
- No second organism (GGO only; generalizability is future work).
- No new rustle features — this evaluates the existing `--vg` capability.

## Open implementation questions (resolve during planning)

- Exact SQANTI3 version + dependency pin set that installs cleanly in this conda setup.
- Whether StringTie2 should get any non-default long-read flags (e.g. `-c`, `-f`) or strict defaults
  (default: strict defaults, documented).
- The precise pinned rustle `--vg` flag set for the headline arm (reference the existing win-stack
  config; record the exact flags).
- `MIN_IDENTITY` / `MIN_COV_FRAC` defaults for U3 and `k` for U5 (start with documented defaults;
  report sensitivity in the plan's validation step).
