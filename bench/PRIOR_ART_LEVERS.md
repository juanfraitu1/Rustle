# What the prior art actually gives us: three open problems, three answers

**Date:** 2026-07-09. Open problems after `YAG_CHECK.md`: (1) the readthrough rule R4 had **no retrocopy
positive control**; (2) **DAZ2 / BPY2** exist but are never assembled; (3) RBMY's **19% ambiguous** rate was
unexplained. Each is addressed below by something already in the reference set.

---

## 1. The retrocopy control — closed, using our own shipped retrocopy filter

`bench/family_def_retrocopy_filter.py` (shipped 2026-06-23) already identifies **one-side-intronless retrocopy
edges** genome-wide — an intronless copy whose parent is spliced. The bake-off
(`project_conflict_criterion_bakeoff`) named the hardest instance: **EEF1A1 → LOC109023808**, where 3347 reads
cross-map from the highly-expressed spliced parent onto the intronless retrocopy.

That is precisely the object R4 (*drop a single-exon transcript iff ≥5 distinct junctions, each with ≥2 reads,
lie entirely inside its span*) is most likely to destroy. Measured:

| gene | locus | reads | MAPQ-0 | distinct junctions | R4 |
|---|---|---|---|---|---|
| **LOC109023808** (EEF1A1 retrocopy, intronless) | `NC_073243.2:97380144-97381766` | 15 | — | **0** | **kept** |
| **LOC134758341** (intronless, **100% MAPQ-0** ⇒ a near-identical paralog exists) | — | 13 | 1.00 | **0** | **kept** |
| EEF1A1 (spliced parent) | `NC_073229.2:97608632-97613071` | 3516 | — | 17 | rule n/a (spliced) |

The parent's spliced reads cannot deposit junctions inside the retrocopy: the retrocopy has no introns, so a
read aligning there aligns contiguously. **The mechanism that makes R4 safe for retrocopies is structural, not
statistical.**

A systematic sweep of 260 expressed intronless genes found only **1** with ≥20% MAPQ-0 reads — retrocopy-like
intronless loci are rare among expressed genes — and R4 keeps it. Combined with the earlier controls:
**0 false positives across 260 expressed intronless genes; artifact minimum remains 14 distinct junctions.**

**The gap flagged in `READTHROUGH_RULE_VALIDATION.md` is closed.** R4 is safe to implement.

---

## 2. Why we lose DAZ2 / BPY2, and what IsoCon does differently

Two things from `reference_isocon_sahlin`, and they are the whole answer.

**(a) IsoCon never sees this artifact, because of its data, not its algorithm.** IsoCon is **targeted**:
RT-PCR primers per family, so its reads are full-length **mature cDNA amplicons**. There is no unspliced
pre-mRNA in the library. Our substrate is whole-transcriptome Iso-Seq, which *does* contain intronic and
readthrough molecules — and those are exactly the reads that build a 164 kb single-exon "transcript" over DAZ.
**The readthrough artifact is a consequence of the data modality, not a flaw in the copy-assignment concept.**
Worth saying plainly to the advisor before he asks why IsoCon has no such problem.

**(b) IsoCon's candidate generation is reference-free, and ours is not.** IsoCon step 1 builds a
**nearest-neighbour graph on reads in edit-distance space**, partitions it into clusters each with one center
read, error-corrects within a cluster, then repartitions and re-corrects until convergence. Candidates come out
of *read sequence space*, never from a reference alignment.

Our de-novo assembly is reference-guided. In the DAZ window every one of the 12 assembled transcripts starts at
~42783130 — all isoforms of **DAZ1**. DAZ2's reads are highly similar to DAZ1's, so minimap2 places them on
DAZ1, and DAZ2 never becomes a locus. The only object over DAZ2 is the readthrough. An edit-distance read
clustering would separate DAZ1 from DAZ2 by their PSVs regardless of where the aligner put them.

**IsoCon flags DAZ as the hardest of its nine families** ("highly repetitive exon structure"). Our DAZ failure
is the same failure, at the same family, for a reason the paper already names.

**(c) IsoCon's filter tells us what to do with the readthrough's reads.** Step 2 assigns each read to its
closest candidate, tests each candidate pairwise at its variant positions, takes the **least-significant
p-value**, iteratively **removes non-significant candidates and reassigns their reads**, and retests. We already
use that per-variant real-vs-error test at the *read* level (the significance gate). We do **not** use the
candidate-removal-and-reassign loop. Dropping the readthrough should re-route its 12–14 reads to the surviving
candidates rather than discard them.

---

## 3. RBMY's 19% ambiguous rate — Eichler's rule reproduces it

`reference_eichler_tbc1d3`: any Iso-Seq read whose primary alignment score beats every other paralog cluster by
**≥ 10** is assigned; the rest are "marked as ambiguous and ignored." Since `copy_assign` now emits `as_margin`
per read, this is directly checkable on RBMY's distal cluster:

| our status | n | `AS ≥ 10` would assign | `AS < 10` (Eichler discards) |
|---|---|---|---|
| assigned | 503 | 406 | **95** |
| tied | 215 | 96 | 118 |
| **ambiguous** | **170** | 49 | **121** |

Our ambiguous reads have **median `as_margin` = 5**, below Eichler's cutoff: **121 of 170 (71%) would be
discarded as ambiguous by the published rule too.** The 19% ambiguous rate is not a bug — RBMY's copies are so
similar that the field-standard criterion also refuses to call them. Our class ≈ his class.

The converse is the method's added value: **95 reads that we assign on PSV evidence would be thrown away by
`AS ≥ 10`**, because AS alone cannot see the distinguishing columns. On this family we are simultaneously more
conservative (49 reads AS would assign, we abstain) and more productive (95 reads AS would discard, we assign
with a certificate).

Eichler also reports that **Minigraph/VG failed** on near-identical paralogs (copies became isolated nodes), and
fell back to phylogenetic grouping. So a graph aligner is not the lever for DAZ/RBMY either — consistent with
our own `rcf611_graph_vs_linear` refutation.

---

## 4. Also relevant, not acted on

- **SDA (Vollger 2019)**: PSVs are separated from het alleles by a **read-depth gate** (total depth ≈ single-copy
  coverage). SDA hits the same K=0 floor and names the same escape — reads > 100 kb. Corroborates that TSPY's
  52% tied at 34 reads is *coverage*, not the identifiability wall.
- **Clair3-RNA**: uneven-coverage normalization for RNA; DP/AD-style reporting for honest evaluation. The
  A→I editing filter is already shipped.

## Recommended order

1. **Implement R4** — the control is now in hand, and DAZ/BPY2 show the artifact suppresses whole families.
2. **Reassign the dropped transcript's reads** (IsoCon step 2), rather than discarding them.
3. **Edit-distance read clustering for candidate generation** (IsoCon step 1) — the only thing that recovers
   DAZ2/BPY2. This is a real piece of work, not a filter.

Related: `bench/YAG_CHECK.md`, `bench/READTHROUGH_RULE_VALIDATION.md`, `bench/family_def_retrocopy_filter.py`.
