//! Variation-graph family analysis for novel gene-family copies.
//!
//! Builds a per-exon family graph (`family_graph`) over paralog copies and
//! drives structural detectors (mosaic, segdup, hidden_copy, positional) plus
//! k-mer-based novel-copy rescue. See `docs/UNMAPPED_FAMILY_RESCUE.md` for why
//! the aligner misses the reads this module rescues.

pub mod family_graph;
pub mod rescue;
pub mod diagnostic;
pub mod positional;
pub mod tandem;
pub mod consensus; // cross-copy consensus error-correction (subtractive precision lever)
pub mod mosaic;
pub mod segdup;
pub mod hidden_copy;
pub mod phasing;
pub mod secondary_index;
pub mod layer2;
pub mod psv_linkage; // Layer-2 "C": within-molecule PSV->junction linkage (PSV-column extraction).
pub mod copy_split; // Joint read-coherence + PSV decomposition into (copy, isoform) units.
pub mod copy_assign; // Copy ASSIGNMENT: resolve a read to a known copy via PSV + junction likelihood.
pub mod family_rescue; // Family-aware copy RESCUE: borrow-strength POA confirm of under-assembled copies.
pub mod family_detect; // Strand-aware de-novo family DETECTION: loci collapse + kmer prefilter + POA edges.
pub mod family_split; // De-novo family DECOMPOSITION: connected components + weighted-modularity Louvain.
pub mod denovo_assemble; // Integration: Pass-1 read-coherence skeletons + general-purpose assemble gate.
pub mod denovo_pipeline; // Integration: de-novo family DETECTION driver (pass1->gate->collapse->detect->split).

pub use family_graph::{ExonClass, FamilyGraph, JunctionEdge};
pub use diagnostic::{RescueClass, classify_internal, classify_external, cigar_has_long_indel};  // Task 6.1
pub use secondary_index::{SecondaryAlignment, SecondaryIndex};
