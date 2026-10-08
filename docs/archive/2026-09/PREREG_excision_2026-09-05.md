# PREREG — EXCISION: remove a copy, follow its reads, detect "a copy should be here" (O3 in O2's vocabulary), 2026-09-05
Question (advisor, Part 1.4 of ADVISOR_QUESTIONS): what happens to the reads of a copy that is ~98 % identical to
another when that copy is absent from the catalog — and is there a pattern that says a copy is missing?
Instrument: `copy_assign` genomic read-star (the shipped default) on `rna_units_v9`, with one copy X removed
from the family's candidate set (copies.tsv/copies.fa re-indexed; nothing else changes).
Targets: NPIP copies 7 (nearest 0.985, 102 MAPQ-60 reads), 11 (0.987, 77), 2 (0.971, 98), 13 (0.966, 32); and
the ZNF569-like copy of the ZNF875 family (nearest 0.991, the > 99 % case).
Measurements, on X's reads (primary placement inside X, MAPQ 60):
- fate after excision: assigned to another copy (SILENT misassignment) / origin-rejected / tied / ambiguous-other.
- detector: among the family's origin-rejected reads grouped by their best candidate Y, the mismatch positions
  (X ops vs Y's padded locus, `minimap2 -x splice --eqx`) shared by >= 3 reads and by >= 50 % of the reads
  covering the position = CONSISTENT sites; density per kb of covered locus. Controls: the same statistic on the
  un-excised run's rejected reads, and on Y's own reads.
- recovery: Y's locus with the consistent sites patched by the majority read allele = the reconstructed copy;
  its identity to the true X locus vs to Y (minimap2).
Predictions:
- P1 silent misassignment <= 5 % of X's MAPQ-60 reads for every target (the origin certificate holds at 1.3–3.4 %
  divergence; the > 99 % case may exceed it — report).
- P2 >= 70 % of X's MAPQ-60 reads are origin-rejected (the rest tied or ambiguous), for targets <= 0.99.
- P3 consistent-site density >= 5 per kb in the excised group vs < 1 per kb in both controls.
- P4 the reconstructed copy is closer to X than to Y (identity), for every target where P3 holds.
Fail rule: P1 fails => the certificate does not protect against a missing copy at that divergence — the
identity at which it fails is the result. P3 fails => rejected reads carry no consistent signature and O3
cannot be built on them.
