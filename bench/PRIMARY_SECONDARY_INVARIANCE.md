# Is the recovered locus set invariant under the arbitrary primary/secondary tie-break?

**Date:** 2026-07-10. **Substrate:** `GGO_mm.bam` vs `GGO.fasta`, `copy_assign` at `b55a30b`, λ=58.
**Motivation (advisor).** The primary/secondary alignment label is arbitrary for a multi-mapping read: at
MAPQ = 0 minimap2 cannot tell which locus is correct and breaks the tie deterministically-but-meaninglessly.
Defining a locus by "≥ 1 primary read" is therefore only sound if the recovered locus set does **not** depend on
that arbitrary choice. This measures whether it does.

**Two tests, agreeing.**
1. **Adversarial bound (all tie-breaks at once).** minimap2 sets MAPQ = 0 *exactly* when the primary is tied.
   So a copy's MAPQ > 0 primary reads are **anchored** — no relabeling can move them — and its MAPQ = 0 reads are
   **swappable**. A copy with ≥ 3 anchored reads clears the locus gate under *every* tie-break.
2. **End-to-end relabel (one reversed tie-break, whole pipeline).** Secondaries in this BAM carry full SEQ with
   soft-clips, so flipping the 0x100 flag between a primary and a genuinely AS-tied secondary yields a
   different-but-valid primary assignment with SEQ/CIGAR/POS untouched (`relabel_ties.py`, promotes the
   maximally-different tied placement). Re-run `copy_assign` on the flipped BAM and compare the recovered copies.

## Result

| family | copies | anchored / swappable (bound) | tie flips | copies after flip | χ_H o→f | **invariant?** |
|---|---|---|---|---|---|---|
| **GSTM** | 3 | 1345 / **0** | **0 (no-op)** | 3 | 3→3 | **YES** (trivial: unique mappers) |
| **MAGEA** | 2 | 621 / **0** | **0 (no-op)** | 2 | 2→2 | **YES** (trivial) |
| **PCDHB** | 5 | 462 / **0** | **0 (no-op)** | 5 | 5→5 | **YES** (trivial) |
| **DAZ** | 2 | copy1 **1 / 17** | 41 | **2** | 2→2 | **YES** (empirical) |
| **RBMY** | 6 | copy3 2/5, copy5 2/0 | 12 | **6** | 6→6 | **YES** (empirical) |
| **TSPY** | 5 | **all copies 0 anchored** | 30 | **2** | 5→**2** | **NO — 5→2** |

## Reading it

**The advisor is right about the premise, and the data locates exactly where it matters.**

- **Divergent families (GSTM, MAGEA, PCDHB) — objection is vacuous.** Every primary read is MAPQ = 60; the
  relabeling flips **zero** reads (a literal identity operation). These are unique-mappers whose "primary" is not
  arbitrary at all. (This is also why they need `--homology-primary`: E_c finds 0 conflict edges — the reads never
  compete.) The locus set is invariant by construction.

- **Moderately-collapsed families (DAZ, RBMY) — arbitrary reads exist, but the copies don't depend on them.**
  DAZ2 has only **1 of 18** anchored primaries (the bound flags it as sensitive), yet flipping all 41 tied reads
  still recovers **2 copies at the same loci** (DAZ2 start shifts ~1.4 kb, a reassembly boundary wobble, not a
  lost locus). RBMY reshuffles read counts across its 6 tandem units but recovers all **6**. The reason: these
  copies are defined by **relabeling-invariant junction structure** (DAZ2's ~31 introns vs DAZ1's ~16; RBMY's
  per-unit junctions), which is present in a read regardless of which locus holds its primary label. Reads move
  between copies; the copies stay.

- **TSPY — the objection lands, exactly as predicted.** Every one of its primaries is MAPQ = 0 (0 anchored on
  all 5 copies). Its ~2.7 kb copies are exonically near-identical with little distinguishing junction structure,
  so per-copy existence rests on *which* copy the arbitrary primary pileup happened to feed. Reversing the
  tie-break **collapses 5 copies to 2** (χ_H 5→2): three copies fall below the gate when their arbitrary primaries
  are reassigned elsewhere. TSPY's 5-copy call is **not** invariant — its copy count ranges 2–5 across valid
  relabelings.

## Conclusion — the invariance criterion, and where the method meets it

"Define a locus by ≥ 1 primary" is unsound *only* in the regime the advisor names — near-identical exonic copies
at MAPQ = 0 — and the experiment shows the pipeline is already invariant everywhere its locus rests on
**relabeling-invariant evidence** (unique placement, or junction/PSV structure): 6 of 7 named families reproduce
their exact copy set under a full adversarial tie-break reversal. It fails precisely where no such evidence
exists (TSPY), which is the same **K=0 frontier** the thesis already flags as genuinely RNA-unresolvable
(`project_k0_frontier_unresolvable`) — there the honest output is aggregate/DNA, not a per-copy RNA count.

The soundness fix is to state locus support as the **relabeling-invariant set** — reads with *any* qualifying
alignment at the locus (primary **or** tied-secondary), gated on junction-incidence — rather than the reads that
happened to land primary there. `locus_support` already moves this way; the remaining primary-only step is the
skeleton `min_reads` gate, and TSPY is the measured case that motivates switching it. Reframed for defense: **a
locus is admissible iff its existence is invariant under primary↔secondary relabeling of tied reads** — a
criterion with no arbitrary threshold, which the method satisfies for every family that is resolvable at the RNA
level at all.

## Now emitted by the tool (the invariance certificate)

`copy_assign` writes the per-copy certificate into `quant.tsv` (two new columns) plus a stderr summary, so the
bound above is produced in-binary rather than by an offline `samtools` pass:

- **anchored_reads** = assigned reads that map uniquely (`mapq > 0`) to the copy — support that survives any
  primary/secondary relabeling of tied reads.
- **tie_invariant** = `anchored_reads ≥ GATE_MIN_READS` (3): TRUE ⟹ the copy exists under EVERY tie-break via
  unique mappers.
- **junction_invariant** = the copy is pinned by ≥ 3 reads carrying a **copy-specific junction**
  (`copy_junction_support`, the `junction_only` signal per copy): its splice structure identifies it regardless
  of the primary label — the second relabeling-invariant mechanism, and the one that actually rescues DAZ2.

A copy is invariant **overall** iff `tie_invariant OR junction_invariant`, and the stderr summary reports that.

Measured on the diagnostic regions (matches the table above), `tie / junction`:

| family | per-copy `tie` / `junction` | overall |
|---|---|---|
| GSTM | all `true` / `false` (unique mappers, not junction-pinned) | 3/3 |
| DAZ | DAZ1 `true`/`true`; **DAZ2 `false`/`true`** (1 unique mapper, junction-defined) | **2/2** |
| TSPY | all `false` / `false` (exonically identical — no copy-specific junction) | 0/5 |
| RFPL | readthroughs `true`/**`false`** (unique mappers, but single-exon ⟹ no junctions) | 4/4 |

`junction_invariant` is precise: it rescues DAZ2 (`false→true` overall), and it does **not** bless RFPL's
readthrough artifacts (single-exon ⟹ no junctions ⟹ `false`) or TSPY's exonically-identical copies (no
*distinctive* junction ⟹ `false`). RFPL's `tie_invariant=true` reflects that its copies are relabeling-robust
(they have unique mappers); their being artifacts is orthogonal and flagged by the Containment warning.

## Reproduce

```
# bound: anchored (MAPQ>0) vs swappable (MAPQ=0) primaries per copy span
samtools view -c -F 2308 -q 1 GGO_mm.bam <copy_span>     # anchored
samtools view -c -F 2308    GGO_mm.bam <copy_span>       # all primary
# end-to-end: reverse the tie-break and re-run
samtools view -b GGO_mm.bam <family_region> > f_orig.bam; samtools index f_orig.bam
python3 relabel_ties.py f_orig.bam f_flip.unsorted.bam; samtools sort -o f_flip.bam f_flip.unsorted.bam; samtools index f_flip.bam
copy_assign --bam f_flip.bam --fasta GGO.fasta --region <family_region> --min-copies 2 \
    --skip-poa-diagnostic --homology-primary --lambda-file <λ> --out f_flip
```

Related: `bench/KNOWN_FAMILY_REGRESSION.md` (the families under test), `project_k0_frontier_unresolvable`,
`project_daz2_locus_support` (DAZ2 = junction-defined), `project_family_def_readconflict` (E_c = ambiguity oracle).
