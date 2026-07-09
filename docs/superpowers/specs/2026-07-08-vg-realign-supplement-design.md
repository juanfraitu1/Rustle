# VG re-align supplement — Design

**Date:** 2026-07-08. **Substrate:** gorilla (GGO) HiFi Iso-Seq; O1 `FamilyGraph`, O4 `absent_copy`, the EM +
reference-free copy number.

## Goal

Close the last linear-reference-bias gap for the families where it matters: reads that the linear aligner
mis-places or drops are re-aligned to O1's per-family variation graph (`FamilyGraph`), and a re-alignment is
accepted only when it is **significantly** better than the linear one. This (a) **corrects** within-family
mis-assignment and (b) **discovers** copies not in the linear genome (a group of reads that consistently
threads a novel path). Bounded per-family, opt-in (`--vg-realign`), default output byte-identical. This is NOT a
genome-wide graph aligner — only candidate reads, only the family's own VG.

## Where the non-linear reference comes from

**O1.** `family_graph::build_family_graph` builds the `FamilyGraph` from a family's copies: exon-equivalence
classes = nodes (a shared exon is ONE node), junctions = edges, each copy = a path. O2 (the EM) only *annotates*
it (`per_copy_cov`). So the supplement re-aligns candidates to an object O1 already produced; O2 and the
reference-free copy number consume the corrected/admitted result.

## Pipeline placement

A supplement that runs, per family, AFTER O1 family detection and the linear O2 evidence build, and feeds its
corrected assignments + admitted copies into the existing EM / `chi_H` / `famcn_readonly` path. Entirely behind
`--vg-realign`; when off, nothing runs and all outputs are byte-identical.

## 1. Candidate collection

A read is a candidate iff its linear fit is poor or absent:
- **Mis-mapped (locus known):** a read mapped at a family locus with **low MAPQ** OR **soft-clip fraction ≥
  `clip_frac`** OR **per-base divergence (NM/de) ≥ `div_max`** to its locus. Routed directly to that family's VG.
- **Unmapped (separable stage):** unmapped reads carry no position, so route them by a **minimizer pre-screen**
  against family consensuses — a read is a candidate for family F if it shares ≥ `min_shared` canonical
  minimizers (`minimizers`, MINIMIZER_K/W) with F's consensus. Then re-align to F's VG. This stage is
  independent: the mis-mapped path ships and is testable without it; the unmapped screen is additive on top.

Thresholds are pre-screen SELECTORS (cheap gates that decide *what to re-align*), NOT accept criteria — the
significance gate (§3) decides what is *admitted*, so a permissive pre-screen only costs compute, not precision.

## 2. Re-align to the FamilyGraph (POA)

POA-align each candidate's sequence to the family VG using the existing aligner in `family_graph.rs`
(`aligner.align(&graph, seq)`). Produce, per read: the best copy-**path** (the exon-class sequence it threads),
its alignment score, and the per-column observed bases along that path (so the significance gate can score it in
the same PSV/allele space the linear evidence uses).

## 3. Significance-gated accept (reuse existing certificates)

- **Per-read correction:** build the read's `ReadFeatures` from its VG-path observations and score it against the
  family's `CopyProfile`s with `read_copy_evidence` (Task-1 machinery). Re-assign the read to the VG copy-path
  **only if** its evidence for that path is significant (`min_p < alpha/(K−1)`, the gate/EM certificate) AND the
  path fit is significantly better than the linear alignment (the linear locus scored in the same space). Reject
  otherwise (leave the linear assignment). Corrected reads flow into the EM as that copy's evidence.
- **Novel-copy discovery:** reads that thread a consistent path NO existing copy explains are pooled by path.
  A novel path is admitted as a reference-absent copy **only via O4's discipline** — reuse
  `absent_copy::admit_candidate` / its gate: ≥ `min_reads` reads, strand symmetry, and non-remap ≥ 98 % to any
  existing locus (so it is genuinely not in the reference). Admitted copies enter the copy set → EM + `chi_H` +
  `famcn_readonly` (the H1 flow). No arbitrary margin anywhere; the certificate is the same one O2/O4 use.

## 4. Output (additive)

`<out>.vg_realign.tsv`, header `read_name  family_id  action  target_copy  path  min_p  linear_locus` — one row
per candidate; `action ∈ {reassigned, novel-copy, rejected}`. Corrected assignments + admitted copies feed the
existing `.assignments.tsv`/`.em*.tsv`/`.famcn_readonly.tsv` emit (their values shift only under `--vg-realign`).
Default (flag off): no `.vg_realign.tsv`, all other outputs byte-identical.

## 5. Validation

1. **Sim (the demonstration).** Extend `bench/sim_genome.py`: plant a copy divergent enough (e.g. 6–10 %) that
   its reads **mis-map to the wrong locus or fail to map** on the linear genome. Run the pipeline with
   `--vg-realign` and assert: (a) the mis-mapped reads are re-assigned to the correct copy-path; (b) the
   divergent copy the linear pipeline missed is admitted as a novel copy and counted by `chi_H`/`famcn_readonly`.
   This is the direct proof of "find copies not in the genome via VG re-align." Emit `bench/VG_REALIGN.md`.
2. **Real (GGO).** Run on `GGO_mm.bam`; report corrected reads + any admitted VG copies. Honest caveat: like O4,
   admission may find few/zero genuine novel copies on GGO (the divergent-absent frontier is data-limited); the
   correction leg should still move reads. State whichever is observed.

## 6. Reuse / non-goals

- Reuses `FamilyGraph` + its POA aligner, `read_copy_evidence`/`min_p`, `absent_copy::admit_candidate`, and the
  EM/`chi_H`/`famcn_readonly` consumers. No new likelihood, no new admission model, no new copy-number metric.
- **Non-goals:** genome-wide VG-alignment; a new graph aligner; re-aligning confidently-well-mapped reads (out
  of the candidate set by design). Copy admission is mechanism — real-data yield is reported, not assumed.

## Files

- **Create** `src/rustle/vg_family/vg_realign.rs`: `collect_candidates` (mis-mapped + optional minimizer
  pre-screen), `realign_to_family_graph` (POA → path + observations), `accept_realignment` (significance vs
  linear, reuse `read_copy_evidence`), novel-path pooling + `absent_copy::admit_candidate`.
- **Modify** the assign path (`copy_assign_pipeline.rs`/`denovo_pipeline.rs`) to invoke it behind a
  `DenovoConfig`/CLI `vg_realign` flag and route corrected assignments + admitted copies into the existing
  structures.
- **Modify** `src/bin/copy_assign.rs`: `--vg-realign` (+ `--vg-realign-clip-frac`, `--vg-realign-div-max`,
  `--vg-realign-min-reads`) and emit `<out>.vg_realign.tsv`.
- **Create** `bench/sim_vg_realign.py` (planted mis-mapping/absent copy), `bench/VG_REALIGN.md`.
- **Test:** in-crate `#[cfg(test)]` for `accept_realignment` (significant path→reassign; insignificant→reject;
  K=0 path→no admission) + candidate selection; integration/sim for the end-to-end recovery.

## Testing (TDD)

Unit: `collect_candidates` selects low-MAPQ/high-clip/high-div reads and skips well-mapped ones; the minimizer
pre-screen routes an unmapped read to the family sharing its minimizers. `accept_realignment`: a read whose
VG-path evidence is significant AND beats the linear locus is `reassigned`; an insignificant one is `rejected`;
a novel path with < `min_reads` is not admitted; with ≥ `min_reads` + strand-symmetry + non-remap it is admitted
(mock the remap). Additivity: `--vg-realign` off ⟹ byte-identical outputs.

## Reproduce

```
copy_assign --bam GGO_mm.bam --fasta GGO.fasta --regions <fam> --vg-realign --lambda-global <λ> --out ca
python bench/sim_vg_realign.py    # planted mis-mapped + genome-absent copy -> recovered via VG re-align
```
