# Pre-registration — LOCUS ASSEMBLY: reconstruct full-length NPIP copies from reads, O2 for the ambiguous ones

**Written 2026-09-18, BEFORE running the pipeline.** User: *"the assembler part must be applied here, can we
try to get all the reads (isoforms) at a locus and reconstruct the full length locus? For ambiguous ones we
need O2 obviously."*

No new algorithm is written. The shipped `copy_assign` already does steps 1–4 of that request (de-novo
assemble, co-located family detection, PSV + copy-specific-junction assignment with assign-or-abstain on
AS-tied multimappers, GTF of the transcripts O2 believes). What is new here is **the deliverable and its
scoring**: per annotated copy, reconstructed intron chain vs annotated chain, plus an explicit
*unrecoverable* flag.

## 0. The ceiling, already measured and binding (§6m4)

| read set | junction coverage (of 249) | complete chains |
|---|---|---|
| primary, MAPQ ≥ 30 | 167 (67.1%) | 5/26 |
| primary, any MAPQ | 188 (75.5%) | **7/26** ← the baseline to beat |
| ALL alignments = perfect O2 | 209 (83.9%) | **10/26** ← the ceiling |

**40 of 249 junctions (16.1%) have no read anywhere in this library.** So 16 of 26 copies CANNOT be
completed by any method on this substrate.

## 1. Run

- Binary: `copy_assign` built `--release` with `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target`.
- BAM: `winloci_data/A119b.t2t.bam` (`-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`).
- Regions: one window per NPIP member, body ± 2 kb — the SAME windows the §6m4 ceiling used, so the two are
  comparable. 27 members (26 spliced + 1 unspliced); the 26 spliced ones are scored.
- Flags: `--regions --families --gtf --read-provenance --flag-missing-copies`, everything else at shipped
  defaults. **No annotation is an input.** The exact command line is recorded in the report.
- Truth (scoring only): the 26 annotated intron chains from
  `Reference/chm13v2.0_RefSeq_full.gff.gz` via `dna_cert.load_nodes()`.
  ⚠ exons are 0-based half-open, read junctions 1-based ⇒ junction = `(exon_end+1, next_exon_start+1)`.

## 2. Pre-registered decision rules

**AS-1 (INTEGRITY, checked first).** Reconstruction must NOT exceed the §6m4 ceiling: ≤ 209/249 junctions
and ≤ 10/26 complete chains. **Exceeding either is a bug or truth leakage, not a success** — it must be
diagnosed before any number is reported.

**AS-2 (does the shipped pipeline reach the ceiling?).** Bars, both required:
(a) junction coverage ≥ **188/249 (75.5%)** — it must beat simply taking every primary read;
(b) complete chains ≥ **7/26**.
**PRE-DECLARED EXPECTATION: (b) is likely to FAIL**, because even a perfect O2 reaches only 10/26 and the
assembler must additionally solve 5′ truncation, which §6m0 showed it does not fully.

**AS-3 (the flagging half — the honest deliverable).** Every copy the pipeline declares reconstructed must
actually be reconstructible. Bar: **zero false "complete" claims** — no copy among the 16 that cannot be
complete may be emitted as a complete chain. Identification must come from the output alone
(`--flag-missing-copies`, provenance, or absence of a transcript), never from the truth.

**AS-4 (no truth leakage).** No GTF/annotation is passed as input. `--gtf` is an OUTPUT flag.

**AS-5 (report the abstentions).** The count of AS-tied molecules abstained must be reported, not hidden.
A pipeline that reaches a bar by assigning everything has not demonstrated O2.

## 3. Out of scope

- No change to `copy_assign` or to any family definition.
- No claim that a second library would help — that is a separate experiment the 16.1% motivates.
