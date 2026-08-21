# docs/ — what is where

25 documents. Start here rather than grepping; most questions are answered by one of the first five.

## Read these first

| document | answers |
|---|---|
| [`NUMBERS.md`](NUMBERS.md) | **Every quotable figure with its substrate.** Look a number up here *before* quoting it — the two headline O1 rates are on different species |
| [`ONE_METHOD.md`](ONE_METHOD.md) | The canonical one-sentence statement per objective. If the code or a slide disagrees, they are wrong |
| [`NEGATIVE_RESULTS_REGISTER.md`](NEGATIVE_RESULTS_REGISTER.md) | **639 things already tried and why they failed. Read before proposing an approach.** §0 lists the nine routes closed on 2026-08-19 |
| [`HANDOFF.md`](HANDOFF.md) | Current state, open items in priority order, scratch directories |
| [`o1_catalog_provenance.md`](o1_catalog_provenance.md) | ⚠ Which catalog a number came from. The shipped 494-family catalog is **not reproducible**; the current default emits 627 |

## The definition

| | |
|---|---|
| [`seeded_family_definition.md`](seeded_family_definition.md) | The O1 definition itself — the normative document |
| [`joint_family_definition.md`](joint_family_definition.md) | DNA/RNA jointness: a property, not the definition |
| [`copy_assignment_definition.md`](copy_assignment_definition.md) | O2's rule |
| [`THESIS_OBJECTIVES.md`](THESIS_OBJECTIVES.md) · [`OBJECTIVES_AND_VERIFICATION.md`](OBJECTIVES_AND_VERIFICATION.md) | Scope, and the per-row verification ledger |
| [`VG_FAMILY_TERMS.md`](VG_FAMILY_TERMS.md) | What "locus", "copy", "isoform", "family" each mean |

## O1 — investigations

| | |
|---|---|
| [`o1_false_positive_rules.md`](o1_false_positive_rules.md) | **All E_r precision guards in one place.** Orientation guard (SHIPPED as default), genome-anchored veto (flag), junction predicate (refuted as a gate), plus the human-panel confirmation |
| [`o1_error_case_census.md`](o1_error_case_census.md) | What is actually wrong with the bad families — 30 definitional, 47 node-construction, 26 not-an-error — plus pathology (a) dissected |
| [`o1_coverage_repair.md`](o1_coverage_repair.md) | Why no threshold fixes the coverage clause (impossibility argument) |
| [`o1_read_evidence_repair.md`](o1_read_evidence_repair.md) | Why reads cannot fix it either |
| [`o1_duplication_provenance_model.md`](o1_duplication_provenance_model.md) | Block-aware provenance graph — representation kept, rooting deferred; includes the hierarchy slice |

## O2 / O3

| | |
|---|---|
| [`o2_reassignment_result.md`](o2_reassignment_result.md) | Reassignment ground truth: route refuted, with its pre-registration appended |
| [`o3_missing_copy_evidence.md`](o3_missing_copy_evidence.md) | The excision positive control and the S2 detector |
| [`o3_unmapped_route.md`](o3_unmapped_route.md) | Unmapped-read route: works in unique sequence, vacuous in paralogous |
| [`o3_overcollapse.md`](o3_overcollapse.md) | Is over-collapse happening, could we see it, and a scoped simulation |
| [`o3_haplotype_cnv_result.md`](o3_haplotype_cnv_result.md) | Haplotype copy-number run — uninformative by its own pre-registration (appended) |

## Reference / infrastructure

| | |
|---|---|
| [`linuxdisk_data_access.md`](linuxdisk_data_access.md) | Where the data lives and how to mount it |
| [`RETIREMENT_AND_MIGRATION.md`](RETIREMENT_AND_MIGRATION.md) | Python→Rust migration map (referenced from 7 source files) |
| [`MECHANISM_WALKTHROUGH_DAZ.md`](MECHANISM_WALKTHROUGH_DAZ.md) | End-to-end worked example on real DAZ reads, for the advisor |
| [`UNMAPPED_FAMILY_RESCUE.md`](UNMAPPED_FAMILY_RESCUE.md) | Why family-similar reads can still be unmapped |
| [`ONE_METHOD.md`](ONE_METHOD.md) · [`stricter_options_catalog.tsv`](stricter_options_catalog.tsv) | — |

## Conventions enforced by tests

* **Catalog provenance.** Any doc quoting `/494` or `/1415` carries the `CATALOG PROVENANCE.` banner.
  Enforced by `tests/docs_catalog_provenance.rs`, which fails the suite on an unlabelled reference.
