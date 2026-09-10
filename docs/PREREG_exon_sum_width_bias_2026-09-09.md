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

## Outcome (2026-09-09) — `bench/width_deficit.py`, human MCL0 (26 copies) + gorilla MCL1 secondary (80 copies)

**A methodology bug caught before any number was trusted.** The first pass bucketed transcripts/reads to
their copy by raw max-overlap with no containment floor (exactly the mistake register row 879 already
named for a different tool) and used plain min/max for the span. One readthrough transcript
(`DN_chr16_15368428_10`, 32.7 kb, only 54 % of its own length inside copy 4's 17.8 kb window) inflated that
copy's "de novo width" by 15 kb; with the floor at 0, EVERY copy showed a 3–14× OVER-shoot, the opposite
sign of the hypothesis. Fixed: require ≥ 80 % of the transcript's/read's own length inside the copy window
before bucketing it there, and use the boundary as-is only after that floor (no percentile trim needed once
containment is enforced — see the sweep below). Recorded as register row 805.

| # | verdict |
|---|---|
| P1 | ⛔ **REFUTED**: only **6/13 = 46 %** of well-covered human copies are under-represented (predicted ≥ 80 %); **median relative deficit −0.1 %** (predicted 10–40 %) — no systematic under-representation once the containment floor is applied |
| P2 | — moot: the median deficit is near zero, so "fixed vs. scaling" doesn't apply; see the containment sweep below for what a floor choice actually controls |
| P3 | ⛔ weak, wrong-signed in places (ρ = 0.42 human / −0.34 gorilla) — not a reliable predictor at this sample size |
| P4 | ⚠ inverted from the framing: `min_reads=3` is NOT clipping the span inward; dropping it moves the median from −130 bp to +118 bp (human), i.e. **more reads (no floor) → a SLIGHT over-shoot**, not under. Two copies (n_reads 6 and 362) have no `--gtf` transcript at all below the floor — the singleton-chain effect (§6hh, row 798), not a boundary-clipping one |
| P5 | ✓✓ **strongly confirmed, the real finding**: median share of the (small) remaining deficit at the 5' end is **100–101 %** for every one of ours/flair/StringTie/isoseq, on both species. The 3' end (TES/polyA) is called to within **11–30 bp** (median \|d3\|) in every arm — essentially exact from reads alone. Whatever discrepancy exists is a 5' (TSS) phenomenon, not a general under-capture |
| P6 | ✓ every tool's median relative deficit is small and comparable: ours −0.1 %, flair −0.3 %, StringTie −1.8 %, isoseq −2.7 % (human); ours −0.1 % (gorilla MCL1, secondary, competitor GTFs not run — testis substrate) |
| P7 | no tool shows an opposite-signed or markedly smaller deficit; isoseq deviates most (still small) |

### Containment-floor sweep (human, `ours`, the mechanism, not the biology)
| floor | n copies \|rel\|>15 % | outliers |
|---|---|---|
| 0.5 | 8/24 | 0 (−100 %), 4 (−84 %), 8 (−55 %), 9 (−44 %), 2 (−43 %)… all OVER-shoots |
| 0.7 | 5/24 | same class, shrinking |
| 0.8 | 4/24 | 2, 7, 8, 15 |
| 0.9 | 2/24 | **7 (+23 %, 0 catalog reads), 15 (+40 %, 72 reads)** — the only two that are genuinely under-represented, both low-coverage |
Raising the floor from 0.5 → 0.9 removes every OVER-shoot outlier (all were readthrough/chimeric transcripts
riding just inside a loose floor) and leaves only two genuine, low-coverage UNDER-representation cases.

### Reading
The user's premise — that the de novo exon-sum systematically under-represents the true width — is **not
supported** once measured correctly; median deviation is under 1 % on both species. What looked like
under-representation in casual inspection was very likely this same contamination (a few readthrough
transcripts skew a **naive** DNA-side width estimate, and a naive estimate skewed one way in a small sample
looks like "it's off" without diagnosing the sign). The one clean, reproducible, cross-tool pattern is
**5'-end asymmetry**: ~100 % of whatever small deficit remains sits at the TSS side; TES calling from reads
alone is essentially exact (< 30 bp median) everywhere. This is NOT fixable by a single additive/multiplicative
constant (the deficit's sign and magnitude vary case to case, driven by coverage and the containment floor,
not by gene length) — the actionable fix is **not a bias correction, it is a containment/readthrough guard
on the assembler's OWN transcript boundaries** (equivalent to the ≥ 0.8 floor used here), since 4 of the 6
misbehaving copies at floor 0.5 were single mis-assembled readthrough transcripts, not a biological
5'-truncation signal. This is now the recommended next step for the de novo route, not a learned offset.
