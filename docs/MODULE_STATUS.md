# Module status — what is SHIPPED vs what was TRIED

Enforced by `module_status_tests`: every module must carry a `//! **STATUS:** <TAG>` line, and
this file must list exactly the module set.

## Consolidated module registry (2026-10-05)

| module | status | note |
|---|---|---|
| `arms.rs` | **INFRASTRUCTURE** | Bulk-merged genome-projection, copy-graph, split, cache, and utility modules. |
| `bridge_regroup.rs` | **SHIPPED-DEFAULT** | Bridge-aware regrouping of assembled GTF. |
| `copy_assign.rs` | **SHIPPED-DEFAULT** | Copy assignment (PSV + junction likelihood); includes `em_copy_assign` and `copy_assign_pipeline`. |
| `denovo_pipeline.rs` | **SHIPPED-DEFAULT** | De-novo family detection driver; includes `denovo_assemble`. |
| `fam_from_gtf.rs` | **OTHER-BINARY** | `mcl_families --from-gtf` library; includes `annotation_families`, `family_container`, `family_relations`. |
| `family_hierarchy.rs` | **OPT-IN** | `mcl_families --emit-hierarchy` (2026-10-09): superfamily / family / expansion / duplicon layers, SD components, fusion loci, external families placed by containment (`docs/archive/2026-10/merged/PREREG_family_hierarchy_flag_2026-10-09.md`). Annotation only. |
| `family_detect.rs` | **SHIPPED-DEFAULT** | Strand-aware de-novo family detection; includes `family_graph`, `family_split`, `mosaic`, `read_conflict`, `family_rescue`. |
| `missing_copy.rs` | **OTHER-BINARY** | O3 RNA-only chain used by `missing_copy_flag` binary. |
| `candidates.rs` | **OTHER-BINARY** | Candidate copies and reference-absent copy admission; used by `candidate_copies` binary. |

**Vendored crate, not a family module:** `src/sedef/` (sedef-rs, segmental duplications; vendored 2026-10-13 from
the user's sandbox) — OTHER-BINARY: `utilities sd-scan` and the driver's `sd` stage; an opt-in input of
`mcl_families --emit-hierarchy --sd-cache`. Not adopted as the SD source (V1 FAIL), the SD98 source (V1b FAIL) or the
duplicon source (V2 FAIL): see `docs/archive/2026-10/merged/PREREG_sedef_rs_layers_2026-10-13.md`.
