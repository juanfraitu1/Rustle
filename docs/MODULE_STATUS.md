# Module status — what is SHIPPED vs what was TRIED

Enforced by `family::module_status_tests`: every module must carry a `//! **STATUS:** <TAG>` line, and
this file must list exactly the module set.

## Consolidated module registry (2026-10-05)

| module | status | note |
|---|---|---|
| `arms.rs` | **INFRASTRUCTURE** | Bulk-merged genome-projection, copy-graph, split, cache, and utility modules. |
| `bridge_regroup.rs` | **SHIPPED-DEFAULT** | Bridge-aware regrouping of assembled GTF. |
| `copy_assign.rs` | **SHIPPED-DEFAULT** | Copy assignment (PSV + junction likelihood); includes `em_copy_assign` and `copy_assign_pipeline`. |
| `denovo_pipeline.rs` | **SHIPPED-DEFAULT** | De-novo family detection driver; includes `denovo_assemble`. |
| `fam_from_gtf.rs` | **OTHER-BINARY** | `mcl_families --from-gtf` library; includes `annotation_families`, `family_container`, `family_relations`. |
| `family_detect.rs` | **SHIPPED-DEFAULT** | Strand-aware de-novo family detection; includes `family_graph`, `family_split`, `mosaic`, `read_conflict`, `family_rescue`. |
| `missing_copy.rs` | **OTHER-BINARY** | O3 RNA-only chain used by `missing_copy_flag` binary. |
| `o3.rs` | **OTHER-BINARY** | O3 candidates and reference-absent copy admission; used by `o3_candidates` binary. |
