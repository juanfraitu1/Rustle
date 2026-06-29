# Read-chain isoform recall: VALIDATED real win over StringTie (2026-06-08)

**Question.** rustle emits ~7% more isoforms than StringTie. Are the extras *real*
(annotated isoforms ST misses = recall win) or *over-collection* (FP)?

**Answer: a genuine, replicated recall win.** rustle recovers ~3–5.5× more exact
annotated (FSM) isoforms than ST misses-the-other-way, and they are *well-expressed*
(ST drops them by parsimony, not coverage). Method: gffcompare vs RefSeq, per-tool,
on the chromosomes assembled earlier.

## FSM (exact annotated isoform) recovery — rustle vs StringTie
| chroms | rustle FSM | ST FSM | **rustle-only (ST misses)** | ST-only | Sn (R/ST) | Pr (R/ST) |
|---|---|---|---|---|---|---|
| NC_073227.2 + NC_073246.2 | 1286 | 1246 | **66** | 26 | 23.0 / 22.3 | 34.5 / 35.6 |
| NC_073243.2 | 637 | 610 | **33** | 6 | 23.7 / 22.7 | 33.4 / 33.1 |
| **total (3 chroms)** | — | — | **99** | **32** | +0.7–1.0 pp | ≈neutral |

- **~99 exact annotated isoforms recovered by rustle that StringTie misses, vs 32 the
  other way** (~3:1 to 5.5:1). These are class-`=` (FSM) = exact matches to the RefSeq
  annotation → **real, not over-collection.**
- **Mechanism is parsimony, not abundance:** rustle-only isoforms have *higher* median
  coverage (4090) than the both-recovered set (3706); 0/66 are low-coverage. StringTie
  drops these well-expressed alternative isoforms via its minimal-isoform-set assembly;
  rustle's read-chain (IsoSeq per-molecule collapse) keeps them.
- **Orthogonal to multimapping/VG:** earlier test showed the isoform advantage is
  concentrated at *singleton* genes (+0.29 tx/gene) not paralogs (+0.07) — this is the
  read-chain lever, not the VG/secondary-alignment infrastructure.

## The precision side (the watch-item)
rustle also emits more novel-junction (`j`) isoforms: +164 (dok7), +30 (NC_073243.2).
These are read-supported but unannotated — either genuine novel isoforms (RefSeq is
incomplete) or IsoSeq over-collection. They cost a little precision on dok7 (34.5 vs
35.6) but precision was *higher* on NC_073243.2 (33.4 vs 33.1), so they are not
systematically wrecking precision. Net: a recall gain at small/neutral precision cost.

## Bottom line
This is the first clean, replicated, **non-artifactual, non-identifiability-bound**
"rustle finds more real transcripts than StringTie" result of the investigation. It is
the **read-chain (IsoSeq) isoform-recall lever**, not multimapping/VG. Genome-wide
extrapolation from ~99 rustle-only FSM on 3 chromosomes → order ~700–1000 annotated
isoforms recovered that StringTie misses.

**Next:** quantify genome-wide (run all chromosomes), and characterize the rustle-only
isoforms (alt-TSS/TES? skipped exon? intron retention?) to frame the recall story; decide
whether to gate the `j` novel calls for precision.

Reproduce: `bench/copy_recovery_eval/iso_validate/` (gffcompare) + `/tmp/iso_paralog.py`.
