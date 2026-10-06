//! Rustle — multi-copy gene-family analysis over long-read RNA (great-ape pan-transcriptomics).
//!
//! The thesis lives in the family-analysis modules: O1 family definition, O2 copy assignment under
//! MAPQ-0 ambiguity, O3 reference-absent / missing copies (detect and flag) — over `genome` plus two
//! small foundational modules (`types` for the deterministic hash aliases, `bam` for BAM opening and
//! CIGAR exon blocks). Every binary imports only these.
//!
//! The legacy StringTie assembler + network-flow stack (~50k lines, ~55 modules) was RETIRED on
//! 2026-07-14; see `docs/RETIREMENT_AND_MIGRATION.md`. Its residue was removed in three steps:
//! 2026-09-20 (StringTie-exact presets, the VG-HMM rescue cluster, the dropped ASJ objective's modules;
//! tag `retired-modules-2026-09-20`) and 2026-09-24 (the bundle/read types, `RunConfig`, the dead half
//! of `bam`, `util::constants`, and the `minimizers`/`bridge_detector`/`repeat_catalog` modules; tag
//! `notebook-2026-09-24`). Build new work in `family`.

pub mod bam;
pub mod genome;
pub mod types;

pub mod family;

// Re-export the family modules themselves so `crate::<module>` paths keep compiling.
pub use family::arms;
pub use family::bridge_regroup;
pub use family::candidates;
pub use family::copy_assign;
pub use family::denovo_pipeline;
pub use family::fam_from_gtf;
pub use family::family_detect;
pub use family::missing_copy;

// Re-export merged submodules at the crate root.
pub use family::arms::catalog_input;
pub use family::arms::collapse_enumerate;
pub use family::arms::copy_graph;
pub use family::arms::copy_split;
pub use family::arms::from_genome;
pub use family::arms::genome_projection;
pub use family::arms::linearize;
pub use family::arms::parcn;
pub use family::arms::readonly_copy_number;
pub use family::arms::run_cache;
pub use family::arms::seq_utils;
pub use family::arms::shared_definition;
pub use family::arms::single_copy;
pub use family::arms::vg_realign;

pub use family::denovo_pipeline::denovo_assemble;

pub use family::fam_from_gtf::annotation_families;

pub use family::family_detect::family_rescue;
pub use family::family_detect::mosaic;
pub use family::family_detect::read_conflict;

pub use family::candidates::absent_copy;
pub use family::candidates::detect;

/// Backwards-compatible alias: the `candidate_copies` binary's library module.
pub use family::candidates as candidate_copies;

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
        "SHIPPED-DEFAULT", // reachable with no env var and no non-default flag — this is the method
        "OPT-IN", // built and wired, behind a flag that defaults OFF — an arm, not the method
        "OTHER-BINARY", // live, but only from a binary other than gw_family_catalog / copy_assign
        "REFUTED", // implemented, MEASURED, and the measurement went against it — kept as an instrument
        "TEST-ONLY", // no non-test callers anywhere in src/ — dead in every shipped binary
        "INFRASTRUCTURE", // shared utility with no independent objective claim
        "AMBIGUOUS", // could not be determined — must not be the resting state of a module
    ];

    fn module_files() -> Vec<(String, String)> {
        // The family modules are inline in src/family.rs. Top-level modules start at column 0:
        // `pub mod X {`. We extract each block up to the next top-level module or EOF.
        let src = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/src/family.rs"))
            .expect("read src/family.rs");
        let mut out = Vec::new();
        let lines: Vec<&str> = src.lines().collect();
        let mut starts: Vec<usize> = Vec::new();
        for (i, line) in lines.iter().enumerate() {
            if line.starts_with("pub mod ") && line.trim().ends_with('{') {
                starts.push(i);
            }
        }
        for (k, &s) in starts.iter().enumerate() {
            let line = lines[s];
            let name = line["pub mod ".len()..line.trim().len() - 1]
                .trim()
                .split_whitespace()
                .next()
                .unwrap_or("");
            let end = starts.get(k + 1).copied().unwrap_or(lines.len());
            let body = &lines[s + 1..end.saturating_sub(1)]; // exclude closing `}`
            out.push((format!("{name}.rs"), body.join("\n")));
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
            match src
                .lines()
                .find(|l| l.trim_start().starts_with("//! **STATUS:**"))
            {
                None => missing.push(name),
                Some(l) => {
                    let rest = l.trim_start().trim_start_matches("//! **STATUS:**").trim();
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
        assert!(
            bad.is_empty(),
            "unknown status tag; allowed are {TAGS:?}: {bad:#?}"
        );
    }

    /// The registry and the module set must not drift apart — a module added or removed without
    /// updating `docs/MODULE_STATUS.md` is exactly how "what is shipped" stops being knowable.
    #[test]
    fn module_status_registry_covers_exactly_the_module_set() {
        let doc = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/docs/MODULE_STATUS.md"
        ))
        .expect("docs/MODULE_STATUS.md must exist");
        let listed: BTreeSet<String> = doc
            .lines()
            .filter_map(|l| l.split('`').nth(1).map(|s| s.to_string()))
            .filter(|s| s.ends_with(".rs"))
            .collect();
        let actual: BTreeSet<String> = module_files().into_iter().map(|(n, _)| n).collect();
        let unlisted: Vec<&String> = actual.difference(&listed).collect();
        let stale: Vec<&String> = listed.difference(&actual).collect();
        assert!(
            unlisted.is_empty(),
            "modules missing from docs/MODULE_STATUS.md: {unlisted:#?}"
        );
        assert!(
            stale.is_empty(),
            "docs/MODULE_STATUS.md lists modules that no longer exist: {stale:#?}"
        );
    }
}
