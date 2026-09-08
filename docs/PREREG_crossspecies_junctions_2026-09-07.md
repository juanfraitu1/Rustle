# PREREG — cross-species replication of the read-through junctions (2026-09-07)

**Written before the run.** md5 in `mcl_ann/adj/readthrough/PREREG_xspecies.md5`.

## Why
Every artefact route nameable inside one library is closed (§6gb addenda 3–4: FLNC, canonical sites, no
microhomology tail, non-duplicate flanks, full-length spanning molecules). What no within-dataset statistic can
exclude is a systematic of THIS library and THIS individual. Orthogonal replication is the only open route.

## Design — junction probes, no liftover
For each junction, a **300 bp probe** = the last 150 bp of the upstream exon + the first 150 bp of the
downstream exon, i.e. the sequence that exists **only if the junction is spliced**. Another species' reads are
streamed against the probe set (`samtools fasta | minimap2 -x map-hifi`), and a read **supports** a probe when
its alignment covers the midpoint with **≥ 50 bp on each side**. No coordinate liftover is involved.

Three probe classes, built identically so nothing distinguishes them but their origin:
- **TEST** — the 42 guarded gorilla read-through junctions.
- **POSITIVE CONTROL** — ordinary introns of the same units (real spliceosomal junctions).
- **NEGATIVE CONTROL** — scrambled junctions: the upstream half of one read-through with the downstream half
  of another. These exist in no genome and must not replicate.

Species: **chimpanzee** (`PTR_mm.bam`, 3.98 M reads), an independent individual, library and species.

## Predictions
| # | prediction |
|---|---|
| **P1** | positive controls replicate at **≥ 0.60** |
| **P2** | negative controls replicate at **≤ 0.05** — if they do not, the probe test is not specific and nothing else may be read from it |
| **P3** | test junctions replicate at a rate **above the negative controls** |

## ⚠ Interpretation fixed in advance — the test is ASYMMETRIC
- **Replication is strong evidence FOR biology**: an independent individual, library and species using the
  same junction cannot be a systematic of the gorilla preparation.
- **Non-replication is WEAK evidence**: read-throughs can be lineage-specific, and NPIP is a fast-evolving,
  copy-number-variable family, so a junction may be genuinely gorilla-specific. A low test rate therefore
  leaves the question open; it does **not** demonstrate an artefact, and must not be reported as if it did.
- Expression differences between the two libraries (tissue, depth) also produce non-replication with no
  bearing on whether the junction is real. The positive-control rate bounds this.
