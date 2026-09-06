# PREREG — L3: single-candidate ties under the pairwise rule (`-p 0`), 2026-09-05

Context: `docs/O1_O2_LOOSE_ENDS.md` L3. Under the genomic read-star a molecule whose only reported chain
(`minimap2 -p 0.3`) is on one candidate is TIED by construction (§6fa: no competitor, no columns). Row 704's
`-p 0` result ("halves the assigned") was measured under the intersection rule, which §6fc replaced. L1 made this
pressing: the EIF3C member of the NPIP cluster (4,774 reads, identity 0.903 to NPIP) is now a candidate and its
reads end as 4,888 ties (arm A of PREREG L1/L2) — the reads with the MOST evidence against every NPIP copy are
labelled as if they had none.
Design: `RUSTLE_STAR_P=0` (every chain reported), nothing else changed, on `sweep_v11/fam_MCL1_073242` (NPIP) and
`fam_MCL7_073242` (LCR16u); the pairwise certificate then sees every competitor's alignment.
Predictions:
- P1 NPIP 62 audited anchors: 0 wrong (arm A: 16 right, 0 wrong).
- P2 reads placed at the dropped EIF3C unit (MAPQ 60): assigned to it rises from 77 to ≥ 4,000 of 4,774 (the
  competitor alignments carry thousands of differing columns), with 0 assigned elsewhere.
- P3 MAPQ-60 placement agreement stays ≥ 0.995 on both families (arm A: 387/390 and 291/291).
- P4 wall time ≤ 4 × arm A (NPIP arm A 234 s max; §6fa's `-p 0` on the unit form took 631 s).
Fail rules: P1 fails ⟹ `-p 0` is not adopted, the wrong calls go to the register with their columns. P3 falls
below 0.995 ⟹ same. P2 fails while P1/P3 hold ⟹ single-candidate ties are not a reporting-floor artefact; the
tie label is retained and the reason recorded. P4 fails ⟹ adopt only if P1–P3 hold and note the cost.
