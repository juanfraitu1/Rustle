# The identifiability boundaries — what is principled-limit, not loose-end

A central thesis claim is that copy-resolution from RNA has **identifiability boundaries** that the theory
*predicts* and the methods *certify*, rather than gaps to be closed by more engineering. Two objectives sit
at those boundaries by their nature. Framing them as boundaries (and reaching-but-not-crossing them, with
certified abstention) is the result — crossing them requires orthogonal data (DNA), explicitly outside the
RNA-only scope.

## K=0 boundary — exonically-identical co-located copies

**Statement.** Two co-located copies whose spliced (exonic) sequences are byte-identical share *no*
distinguishing feature in RNA: a read carries the same bases/junctions whichever copy it came from. By
construction there is no PSV and no copy-specific junction → the read is **unidentifiable**.

**Theorem-predicted, not a failure.** This is the K=0 vertex of the identifiability spectrum (copy_assignment
theory: a copy is resolvable iff it carries a distinguishing feature; the Strong-Separation / K-frontier
results characterize exactly when a unique cover exists). K=0 violates the precondition.

**Empirically confirmed.** `project_k0_frontier`: at exonically-identical co-located loci the locus
references are NM:i:0 (byte-identical) and **0/386 reads** are resolvable — a hard floor, not a tuning issue.

**Now CERTIFIED by the significance gate.** The IsoCon gate routes such reads to `Tied` because every
competitor has an empty distinguishing set → `min_p_value = 1.0 ≥ α`. So the method does not guess (no 1/k);
it *proves* unresolvability per read. Validated on the sim5x **K=0 rung → 100% Tied**.

**Escapes (require leaving RNA-only or aggregating).** (i) reference-vs-reference NM gate (use the assembly,
not the reads); (ii) aggregate per-family quantification (the count is identifiable even when the partition
is not — copy_assignment theory §6b: famCN is a sufficient statistic, parCN is not); (iii) DNA.

## O4 boundary — copy vs allele (reference-absent confirmation)

**Statement.** Distinguishing a genuine extra **copy** from an extra **allele** (heterozygous site) is
information-theoretically impossible from RNA alone: both add one haplotype's worth of sequence to the pool.
RNA carries expression, not copy number.

**So O4 is DETECT-and-FLAG, bounded by design.** We *can* flag candidate reference-absent copies from
RNA-visible signals — collapsed-CNV (≥12 balanced co-segregating alt columns → 18 strong candidates) and
divergent-unmapped — and quantify the flag precision (**FP bound 7.39%**, `project_objectives_status`). What
RNA cannot do is **confirm** copy-vs-allele; that needs **DNA copy-number (parCN, Soto-2025-style)**. The
divergent-unmapped route was additionally **dry on T2T** (the reference already contains them), a second,
data-specific boundary.

**Framing.** O4 attained = *the RNA-attainable half* (detect + flag + FP-bound). The unattainable half
(confirmation) is an information boundary, and the catalog flags candidates *for* a DNA follow-up rather than
overclaiming them — the honest, theorem-consistent position.

## Why this is a contribution, not a shortfall

The methods **assign-or-abstain up to the boundary** and **certify** the boundary per read (Tied via
`min_p ≥ α`; flag-not-confirm for O4). That is exactly the Canzar-aesthetic result: a clean combinatorial
identifiability criterion, provable boundaries, no arbitrary 1/k apportionment across the unidentifiable.
Reaching the boundary with a certificate is the deliverable; the orthogonal-data escapes are named, not
pretended.
