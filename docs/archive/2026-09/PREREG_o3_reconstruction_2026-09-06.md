# PREREG — reconstruction over every flagged pair (O3 item 2), 2026-09-06

For every `missing_copy` pair of the O3 pass (after §6fm's reruns): patch the candidate Y's locus over the covered
stretch with the consistent alleles (§6ff recovery, verbatim: the longest run of read-covered positions, gaps
≤ 50 bp; majority allele at each consistent site) and align the reconstruction to the WHOLE assembly
(`minimap2 -x asm20 -c -N 20 -p 0`), one run per sweep. Per flag report: stretch length, sites patched, identity
of the reconstruction to Y over the stretch, best assembly hit outside Y's extent (position, identity, the
catalog unit or RefSeq gene there), and the verdict:
- `existing_locus` — a hit outside Y's extent with identity ≥ the reconstruction's identity to Y: the reads come
  from a locus the assembly has (an O1 roster gap or an MCL cut, cross-referenced to the admission loci);
- `reference_absent` — no hit outside Y's extent reaches the identity to Y: the assembly holds no better
  explanation, a copy (or haplotype) absent from the reference;
- `unresolved` — too few sites to patch (< 3) or a stretch < 300 bp.
Predictions: P1 on the five excision positives the best hit outside Y is the true excised locus X (identity to X
≥ the §6ff recovery values 0.900 / 0.925 / 0.884 / 0.813 / 0.994). P2 on the intact catalogs ≥ half of the
flags with a patched stretch resolve to `existing_locus` (the admission classes suggest most missing origins
exist in the assembly). P3 the `reference_absent` flags are reported with their stretch and identity; no number is
chosen after looking. Fail rule: P1 fails ⟹ the reconstruction does not implement §6ff, fix first.
