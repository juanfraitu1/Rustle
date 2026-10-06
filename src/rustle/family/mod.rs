//! Family analysis for novel gene-family copies.
//!
//! Family detection and definition (O1), copy assignment under MAPQ-0 ambiguity (O2) and
//! missing/reference-absent copies (O3), with the structural detectors (mosaic, hidden_copy)
//! and k-mer-based novel-copy rescue they use. See `docs/o3_missing_copy_evidence.md` for why
//! the aligner misses the reads this module rescues.
//!
//! **STATUS:** INFRASTRUCTURE  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

pub mod arms;
pub mod bridge_regroup; // OPT-IN `copy_assign --bridge-regroup f1|f1v2`: bridge-aware regrouping of the assembled GTF, bridges kept as fusion_of relations and out of the families input (port of bench/f1_bridge.py + f1v2.py).
pub mod copy_assign; // Copy ASSIGNMENT: resolve a read to a known copy via PSV + junction likelihood.
pub mod denovo_pipeline; // Integration: de-novo family DETECTION driver (pass1->gate->collapse->detect->split).
pub mod fam_from_gtf; // the `--from-gtf` family stage as a library (loci, all-vs-all, copy table), imported by mcl_families.
pub mod family_detect; // Strand-aware de-novo family DETECTION: loci collapse + kmer prefilter + POA edges.
pub mod missing_copy; // O3 RNA-only chain: divergence mixture -> PSV consistency -> patched consensus -> home search -> screens -> verdict (§6ze).
pub mod o3; // O3 candidate copies and reference-absent admission gate.

// Re-export merged submodules so existing `crate::family::<name>` paths keep compiling.
pub use arms::catalog_input as catalog_input;
pub use arms::collapse_enumerate as collapse_enumerate;
pub use arms::copy_graph as copy_graph;
pub use arms::copy_split as copy_split;
pub use arms::from_genome as from_genome;
pub use arms::genome_projection as genome_projection;
pub use arms::linearize as linearize;
pub use arms::parcn as parcn;
pub use arms::readonly_copy_number as readonly_copy_number;
pub use arms::run_cache as run_cache;
pub use arms::seq_utils as seq_utils;
pub use arms::shared_definition as shared_definition;
pub use arms::single_copy as single_copy;
pub use arms::vg_realign as vg_realign;

pub use denovo_pipeline::denovo_assemble as denovo_assemble;

pub use fam_from_gtf::annotation_families as annotation_families;

pub use family_detect::family_rescue as family_rescue;
pub use family_detect::mosaic as mosaic;
pub use family_detect::read_conflict as read_conflict;

pub use o3::absent_copy as absent_copy;
pub use o3::o3_candidates as o3_candidates;


#[cfg(test)]
mod module_status_tests {
    //! Enforcement for `docs/MODULE_STATUS.md`.
    //!
    //! ⚠⚠ **WHY THIS EXISTS.** A `//!` header is a CLAIM, and a reachability survey of all 53 modules
    //! (ledger §6dj) found **30 of them disagreeing with what their callers actually support** — several
    //! describing themselves in the present tense as live analysis stages while having ZERO production
    //! callers. Shipped docstrings have also carried findings that were later retracted. These tests
    //! cannot verify a header's prose, but they CAN force every module to declare what it is, and stop
    //! the registry drifting away from the module set.

    use std::collections::BTreeSet;

    const TAGS: &[&str] = &[
        "SHIPPED-DEFAULT",   // reachable with no env var and no non-default flag — this is the method
        "OPT-IN",            // built and wired, behind a flag that defaults OFF — an arm, not the method
        "OTHER-BINARY",      // live, but only from a binary other than gw_family_catalog / copy_assign
        "REFUTED",           // implemented, MEASURED, and the measurement went against it — kept as an instrument
        "TEST-ONLY",         // no non-test callers anywhere in src/ — dead in every shipped binary
        "INFRASTRUCTURE",    // shared utility with no independent objective claim
        "AMBIGUOUS",         // could not be determined — must not be the resting state of a module
    ];

    fn module_files() -> Vec<(String, String)> {
        let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/src/rustle/family");
        let mut out = Vec::new();
        for e in std::fs::read_dir(dir).expect("family dir") {
            let p = e.expect("dir entry").path();
            if p.extension().and_then(|x| x.to_str()) != Some("rs") {
                continue;
            }
            let name = p.file_name().unwrap().to_string_lossy().to_string();
            if name == "mod.rs" {
                continue;
            }
            out.push((name, std::fs::read_to_string(&p).expect("read module")));
        }
        out.sort();
        out
    }

    /// Every module must say what it is. A new module cannot be added without declaring whether it
    /// ships, is an opt-in arm, or is not reachable at all.
    #[test]
    fn every_module_declares_a_status() {
        let mut missing = Vec::new();
        let mut bad = Vec::new();
        for (name, src) in module_files() {
            match src.lines().find(|l| l.starts_with("//! **STATUS:**")) {
                None => missing.push(name),
                Some(l) => {
                    let rest = l.trim_start_matches("//! **STATUS:**").trim();
                    if !TAGS.iter().any(|t| rest.starts_with(t)) {
                        bad.push(format!("{name}: {rest}"));
                    }
                }
            }
        }
        assert!(
            missing.is_empty(),
            "these modules declare no `//! **STATUS:** <TAG>` line (see docs/MODULE_STATUS.md): {missing:#?}"
        );
        assert!(bad.is_empty(), "unknown status tag; allowed are {TAGS:?}: {bad:#?}");
    }

    /// The registry and the module set must not drift apart — a module added or removed without
    /// updating `docs/MODULE_STATUS.md` is exactly how "what is shipped" stops being knowable.
    #[test]
    fn module_status_registry_covers_exactly_the_module_set() {
        let doc = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/docs/MODULE_STATUS.md"))
            .expect("docs/MODULE_STATUS.md must exist");
        let listed: BTreeSet<String> = doc
            .lines()
            .filter_map(|l| l.split('`').nth(1).map(|s| s.to_string()))
            .filter(|s| s.ends_with(".rs"))
            .collect();
        let actual: BTreeSet<String> = module_files().into_iter().map(|(n, _)| n).collect();
        let unlisted: Vec<&String> = actual.difference(&listed).collect();
        let stale: Vec<&String> = listed.difference(&actual).collect();
        assert!(unlisted.is_empty(), "modules missing from docs/MODULE_STATUS.md: {unlisted:#?}");
        assert!(stale.is_empty(), "docs/MODULE_STATUS.md lists modules that no longer exist: {stale:#?}");
    }
}
