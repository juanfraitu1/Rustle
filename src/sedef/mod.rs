//! **STATUS:** OTHER-BINARY — `utilities sd-scan` and the driver's `sd` stage; opt-in input of
//! `mcl_families --emit-hierarchy --sd-cache`. Adoption as the SD / duplicon source is decided by
//! `docs/PREREG_sedef_rs_layers_2026-10-13.md`.
//!
//! Segmental duplications (SDs) of a genome with sedef-rs, a pure-Rust re-implementation of the SEDEF/BISER cascade
//! (seed -> chain -> refine -> decompose), vendored on 2026-10-13 from the user's sandbox
//! `/mnt/c/Users/jfris/Desktop/sd/sedef-rs` (docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md). The
//! sandbox stays the place to experiment; changes are ported here deliberately and listed below.
//!
//! Vendored files (md5 of the sandbox copies): anchors.rs 4b31457182ff3b13547434ac20c2384d, chain.rs
//! eadd615dacc6da648ffc480114b8fc03, decompose.rs dcbf640590c583ae1d4ff6545230e0dd, fasta.rs
//! 7a364ba96b2c631ea224bd7a3d1e2a01, kmer.rs f56c0eedb15d36fe5455d4f19befb896, refine.rs
//! 67e9c83a01f7fc4763b204221083ff90, seed.rs 9d79d48c254691514d0672dd35711a19. The sandbox's main.rs (its four-stage
//! CLI) is not vendored.
//!
//! Changes against the sandbox, each with its test:
//! - every `crate::<module>` path reads `crate::sedef::<module>`.
//! - `chain::chain_and_merge` (pass 2) and `merge_blocks` join chains / blocks only within one chromosome pair: the
//!   sandbox merged by window coordinates alone, fusing blocks of different chromosome pairs under the names its
//!   HashMap met first (test `chains_and_blocks_of_different_chromosome_pairs_never_merge`).
//! - `seed::build_windows` keeps only each sketch's hashes (the sandbox kept the 2,000-hash buffer, 16 KB per window),
//!   and `seed::find_candidates` counts pairs in `SeedParams::count_shards` shards, releasing each sketch as its entries
//!   are written: the same candidates, bounded memory (tests `window_sketches_keep_only_their_hashes`,
//!   `sharded_counting_gives_the_one_map_candidates`, `counting_releases_every_sketch`).

pub mod anchors;
pub mod chain;
pub mod decompose;
pub mod fasta;
pub mod kmer;
pub mod refine;
pub mod seed;
