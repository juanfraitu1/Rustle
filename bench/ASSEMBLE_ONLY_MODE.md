# `--assemble-only`: the assembler product, with none of the all-vs-all work

**What it is.** A mode that does ONLY what an assembler does — define loci and cluster reads into
transcripts — and skips family detection, homology refinement and copy assignment entirely.

```
copy_assign --bam R.bam --fasta genome.fa --regions R.txt --assemble-only --out OUT
```

It implies `--gtf`; the GTF is the product.

## What it runs, and what it skips

| stage | full mode | `--assemble-only` |
|---|---|---|
| read collection (`reads_in_region`) | yes | **yes** |
| intron-chain skeletons (`pass1_skeletons_widened`) | yes | **yes** |
| `assemble_gate` (motif, span, length, read floor) | yes | **yes** |
| `collapse_loci_groups` (locus grouping) | yes | **yes** |
| GTF emission | yes | **yes** |
| AS-tied gate | yes | **skipped** |
| `detect_and_assign` — co-located family detection, POA homology, E_r refine | yes | **skipped** |
| copy assignment / certificates | yes | **skipped** |

`<out>.families.tsv` and `<out>.assignments.tsv` are written **empty by construction** — there is no
assignment in this mode, and an empty file is the honest record of that rather than a missing one.

## Measured (27 NPIP windows, `RUSTLE_JUNCTION_MAJORITY=1 --read-isoform-k 3`)

| | wall | max RSS | GTF rows |
|---|---|---|---|
| full mode | 53.3 s | 525,840 KB | 14,730 |
| **`--assemble-only`** | **17.3 s** | 520,640 KB | 14,730 |

⭐**The GTF is BYTE-IDENTICAL between the two modes here, at 3.1x the speed.** Lib suite **883 passed / 0 failed**.

⚠**CORRECTION — identity holds only where NO families are detected.** On a region where detection fires
(chr16:11.9-19.0 Mb, full mode finds 8 families / 1,081 assignments), the GTFs DIFFER: assemble-only emits
**2,791 transcripts vs full mode's 2,777** — **14 extra, 0 missing, and all 2,777 shared ones byte-identical**.
The cause is `--gtf-copy-set` (default ON), which uses the detected families to drop "phantom" transcripts at
copies with no evidence; with no families that rule cannot fire. The speedup there is far larger: **7.2 s vs
557 s (77x)**. So: assemble-only is a strict superset of the full mode's transcripts, never a subset.

## Composes with the assembly knobs

`--read-isoform-k` (§6m6), `RUSTLE_JUNCTION_MAJORITY` (§6m8), `RUSTLE_GTF_SECONDARY` (§6n2),
`RUSTLE_GATE_CENSUS` (§6m7) all apply unchanged — those are all assembly-path settings.

## Why this exists

The advisor has repeatedly framed the work as an assembler (`project_assembler_framing_0909`). This gives
that framing a first-class entry point instead of requiring the full multi-copy machinery to be run and its
family output ignored. It is also the honest scope statement: **on hard loci the assembly is already
competitive** — §6hz measured 0.846 against flair 0.574, StringTie 0.502 and isoseq3 0.755 — and this mode
is exactly the thing those tools are compared against.

⚠ It is NOT the multi-copy product. Tied-AS multimapper resolution, copy assignment and the family
definition are the project's actual contribution and all live outside this mode.
