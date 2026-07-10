# DAZ2 recovered — by fixing the assembly gate, not by any of the three proposed routes

**Date:** 2026-07-10. Follow-up to `bench/COLLAPSE_GATE_VALIDATION.md`, `bench/PRIOR_ART_LEVERS.md`.
An 11-agent probe tested three candidate routes and adversarially attacked the two empirical ones from three
lenses each. **All six refutation attempts failed to refute. All three routes are dead.** The fix came from the
fourth probe, which nobody had proposed.

## The three routes, refuted

**SDA depth-excess, transplanted to RNA — DEAD.** RNA depth is copy number × expression, and expression spans a
measured **9,106-fold** range across single-copy genes. DAZ1 reaches only **3.01×** the single-copy median and is
out-depthed by **24 of 100** random single-copy genes; it sits ~50× below SDA's own `mean + 3 s.d.` threshold.
Meanwhile single-copy **EEF1A1 scores 61.4×**. Depth-excess would rank EEF1A1, TSPYL1, DERPC and ATXN7L3B as the
strongest "collapses" and miss the one real one. It does not reject the false positive — **it amplifies it.**

**Cluster → consensus → realign — REFUTED, on its own accept case.** The DAZ1 read pile has **exactly 1 PSV
column** out of 5,096 columns at depth ≥ 10 (robust at minor-allele support 5/10/20). DAZ1 and DAZ2 are exonically
near-identical: of DAZ2's own 20 primary reads, **16 of 20 realign to the two copies as exact ties** (median AS gap
**0**), and 11 of 20 actually prefer DAZ1. The pile never splits, so the discriminator calls DAZ1 single-copy. It
does correctly reject EEF1A1 (both het clusters realign best to EEF1A1 itself; the pseudogenes are distant
secondaries at AS 870–984 vs 3472) — but rejecting the controls is worthless when the accept side fails.

A verifier landed the sharpest blow: **EEF1A1 downsampled to DAZ1's 200 reads yields PSV columns = 1**, identical
to DAZ1, across three seeds. So EEF1A1's apparent "richer substructure" — the χ(H) = 7 that killed the collapse
gate — was **a depth confound all along**, not biology.

**IsoCon step 1 — conditional, and the condition fails.** Reference-free edit-distance clustering would defeat
DAZ2's fragmentation, but under K = 0 it merges DAZ1 and DAZ2 into one candidate; the authors state the floor
themselves. Worse, it clusters by *sequence* and therefore discards the one signal that does separate DAZ2 — its
**distinct genomic junction coordinates (0 shared with DAZ1)**. Strictly more expensive and strictly weaker.

## The actual defect: a threshold applied to the wrong object

`assemble_gate` tested `GATE_MIN_READS = 3` against a **single intron chain**, i.e. against one *isoform*.

At DAZ2, 17 primary reads land inside the annotated span (16 at MAPQ 0). Twelve are spliced, and they fragment
into **9 distinct intron chains whose best support is 2 reads**. Every chain died at the gate. Yet all 12 share
the terminal junction `(42939630, 42943604)`, and **DAZ2 shares zero junctions with DAZ1** (58 vs 19 junctions,
exact overlap 0). A junction-coherent locus with 12 reads was shattered into sub-threshold isoform fragments and
never assembled.

The gate's threshold describes a *locus*. It was being applied to an *isoform*.

## The fix

`denovo_assemble::locus_support` — connected components of the **junction-incidence graph** (two skeletons
adjacent iff they share an exact `(chrom, donor, acceptor)`). `assemble_gate` now tests `min_reads` against the
component's summed support. This replaces a per-isoform threshold with a graph operation; the threshold that
remains applies to the object it was always meant to describe.

Two properties make it safe, and each is a test:

- **Single-exon skeletons never pool.** They carry no junctions, so each is its own component. Keying on the
  empty intron chain would union every unspliced read on a chromosome (measured: 746 reads across 44.3 Mbp).
- **Chimeric bridges never pool.** A skeleton that shares junctions with two skeletons whose spans are
  **disjoint** splices across two loci and belongs to neither; it keeps only its own read count. Exact, not
  thresholded — two spans either intersect or they do not.

The second guard was not foreseen; **the data forced it.** With naive pooling, GSTM silently lost GSTM5. The
culprit was a **2-read spliced transcript spanning `129191743–129222260`** carrying junctions of both GSTM5
(`129191742-129197751`) and GSTM1 (`129216297-129222748`). It bridged them in the incidence graph, inherited
their combined support, cleared the gate, and span-aware collapse then merged the two real copies into one locus.
The synthesis had predicted exactly this attack ("its correctness rests entirely on the claim that shared-junction
pooling never bridges paralogs"). It does bridge them, and the disjoint-span test is what stops it.

Separately, the readthrough filter now runs on **transcripts, before locus collapse**, not on reps after. Filtering
reps let the 298 kb DAZ readthrough become the representative of every locus it spanned — it absorbed DAZ2's
transcripts into DAZ1's group, and dropping the rep deleted DAZ2 with it.

## Result

| | before | after |
|---|---|---|
| **DAZ** | 1 rep, **0 families** | **2 copies** — `42783128-42859657` (DAZ1) and `42899568-42945549` (**DAZ2**, inside its annotated span, 0 shared junctions) — 1 family, **2213 / 2353 assigned**, 139 tied, 1 ambiguous |
| GSTM | 3 copies, 2675 assignments | **3 copies, 2675 assignments** (unchanged) |
| MAGEA | 1 family, 931 | unchanged |
| RBMY / TSPY | 888 / 218 | unchanged |
| r1 | 1 family, 665 | unchanged |
| **TSPYL1, EEF1A1** | 0 families | **0 families** (controls hold) |
| planted K=0 sim | 360 assigned / 1184 tied | unchanged |
| r4 | 1 family, 220 | **2 families, 818** — one pre-existing family byte-identical, one NEW 4-copy family over annotated loci (2 protein-coding, 2 lncRNA), carrying one flagged `Containment` overlap |

873 tests, byte-parity suite intact. `--no-pool-locus-support` reproduces the pre-fix gate exactly.

## What is NOT fixed, and will not be

**BPY2 is unrecoverable from this RNA.** Its window yields the DAZ family, not BPY2: copy A has **77 alignments,
all secondary, ~0 primaries**; expression is essentially nil. Position is not handed to us, depth-excess is
refuted, and the sequence routes need reads that do not exist. This is a data limit, not a method limit.

**DAZ's 139 tied reads are the K = 0 wall**, and that number is honest: DAZ1 and DAZ2 are exonically
near-identical (1 PSV column; DAZ2's own reads tie at median AS gap 0). What resolves the 2213 assigned reads is
**copy-specific junction structure**, not exonic PSVs — the junction term of the assignment gate doing exactly
what it exists for.

⚠ **Correction.** An earlier draft attributed the junction difference to the reps carrying "31 and 16 introns,"
as though the copies differed structurally. They do not. The recovered DAZ2 model is **5′-truncated**: its rep
starts at 42899568 against an annotated start of 42879918, covering **70.1%** of the annotated span. The 5′ gap
has mean primary depth 0.17× and contains exactly one read with a real intron chain — one read, below
`GATE_MIN_READS`, so the 5′ exons never assemble. The 16-vs-31 intron count is truncation, not divergence.

DAZ2 is nonetheless a **genuine second copy**, not spillover: minimap2 `-x asm20` aligns DAZ2's genomic span to
DAZ1's as a single alignment at **85.9% identity over 99.9% of DAZ1** (inverted, consistent with DAZ1 `-` /
DAZ2 `+`), and the two reps share **0 reads and 0 junctions**.

## The honest attack surface

The fix still *contains* `GATE_MIN_READS = 3`; it moves a magic number rather than removing one. The defensible
claim is narrower and stronger: the number now bounds **read support for a locus**, which is what a support
threshold means, and the operator that defines a locus is a graph property (connected components of junction
incidence, minus chimeric bridges) with two exactly-stated exclusions and no tunable constant. The remaining
constants — `GATE_MIN_READS`, `MAX_SPAN` — are unchanged from before and are not what this fix rests on.

Verified on one recovered copy (DAZ2) and one preserved copy (GSTM5). Prevalence genome-wide is unmeasured.

Related: `bench/COLLAPSE_GATE_VALIDATION.md`, `bench/PRIOR_ART_LEVERS.md`, `bench/YAG_CHECK.md`,
`reference_sda_vollger`, `reference_isocon_sahlin`.
