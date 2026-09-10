# PREREG — the de novo exon-sum width deficit: how much, which end, and does it correct? (2026-09-09, before the measurement)

**Question (user):** the de novo route's read-only span under-represents a copy's true (annotated) width.
By how much? Is it a fixed offset (fixable by a learned correction, no annotation needed at inference) or
does it scale with the copy? Do isoseq/flair/StringTie show the same pattern?

**Truth.** `copies16.tsv` (human MCL0) / gorilla `fam_MCL{1,7,58}` catalogs: `start`,`end` = the annotation's
own TSS–TES span (§9-03 pivot: `mcl_families` nodes are RefSeq gene/pseudogene records, `exons` is the
GFF's own exon list) — the trusted width for this comparison.

**De novo estimate, two forms, both annotation-free in their OWN coordinates** (only the copy interval is
used to bucket them, exactly as `tool_bakeoff.py`'s `copy_of()` already does — max positional overlap):
1. **`gtf`** — our own `--gtf` assembly (`ours_final2_g2.gtf`; `pass1_skeletons`/`assemble_gate` take only
   reads + genome, confirmed no `--families`/`--gff` dependency in their signatures). Per copy: min(transcript
   start) .. max(transcript end) over every family transcript whose `copy_index` matches (the gate's own
   positional tag), **before** any evidence/phantom filtering — the min_reads=3 collapse gate is the only
   filter active.
2. **`raw`** — every primary (`-F 2308`) whose alignment overlaps the copy's window, bucketed to its
   single nearest copy by max positional overlap (no min_reads floor, no chain-collapse): min(read start)
   .. max(read end). Isolates whether `min_reads` (not read *coverage*) is doing the clipping.

**Deficit**, strand-aware (5' / 3' relative to the gene, not chromosome coordinate direction):
`deficit_5 = |denovo_5'_end − true_5'_end|`, `deficit_3` likewise; `deficit_total = true_width − denovo_width`.

**Competitors, same rule**, transcripts from `flair_family.gtf` / `stringtie_family.gtf` / `isoseq_family.gff`
bucketed to copies by the identical positional-overlap rule (`tool_bakeoff.py::copy_of`).

Human MCL0 only (the largest catalog); gorilla reported as a secondary table, never pooled.

| # | prediction | refuted by |
|---|---|---|
| P1 | the `gtf` form under-represents true width in **≥ 80 %** of copies with `n_reads ≥ 10`, median relative deficit **10–40 %** | < 60 % under, or median outside 5–60 % |
| P2 | the deficit is **NOT a fixed additive constant**: coefficient of variation of `deficit_total` (bp) across copies is **> 0.5** (i.e. it scales with something, not a flat correction) | CV < 0.3 (a single constant would then fix most of it) |
| P3 | the deficit correlates with **`n_reads`** (more reads → smaller relative deficit, heterogeneous TSS/TES gets covered) — Spearman ρ ≤ −0.3 | ρ > −0.1 |
| P4 | `raw` (no min_reads floor) recovers **≥ 30 %** of the `gtf` form's median deficit — i.e. `min_reads=3` is a real contributor, not the whole story | `raw` and `gtf` deficits agree within 10 % (min_reads is irrelevant) |
| P5 | the deficit is **asymmetric**: one end (5' or 3') carries **≥ 65 %** of the total deficit, consistently across copies (sign test) — matching a known biology (IsoSeq/HiFi cDNA 5'-truncation or 3' internal priming) rather than random noise on both ends | neither end dominates (45–55 % split) |
| P6 | **isoseq / flair / StringTie show the SAME direction and a comparable magnitude** of deficit (within 2× of ours on the median), because the deficit is a property of the READS, not of any one assembler's algorithm | any competitor's median deficit is < 1/2 ours, or in the opposite direction |
| P7 | among competitors, whichever pads its calls (if any) shows a **smaller or reversed** deficit vs the others — a soft check for whether "snap-and-extend" behaviour is visible in an existing tool | — (report only) |
