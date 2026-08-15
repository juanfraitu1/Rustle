//! `--seed`: a QUERY over the emitted catalog. **Not** a term in the family definition.
//!
//! ## Why this is a projection and not a pipeline
//!
//! The seeded probe in `bench/crossspecies/seed_family.sh` builds its node set *from the seed*
//! (`V(s)` = intervals receiving ≥ β aligned bp from `s`). That construction is measurably
//! seed-dependent and it is **not** what the shipped binary does:
//!
//! * strict seed-invariance FAILS on it — seeding independently from each of the 19 annotated human
//!   NPIP genes gives only **4 distinct `F(s)` as sets of loci**, agreeing on **64 of 171** seed
//!   pairs, because `V(s)` inherits the *seed's own length* (seed NPIPB8, 10.6 kb → 27 loci of
//!   ~9.9–10.6 kb; seed NPIPB11, 25.2 kb → loci of 8.8–25.6 kb). Membership at the level of
//!   annotated genes is invariant (19/19 seeds → one gene set); locus EXTENT is not;
//! * on the same probe a gorilla single-copy CONTROL (`SDHA`) returns a **14-locus "family"** at
//!   induced density 0.747 — legal at both γ = 0.20 and γ = 0.40 — where the shipped binary
//!   correctly emits nothing. 52 of its 68 accepted edges have "coverage" > 1.0 (up to 2.019),
//!   because `(qe-qs)/min(|u|,|v|)` is unbounded above;
//! * the seeded probe and the shipped catalog disagree wholesale on gorilla (MAGEA4 seed → 1 locus
//!   vs 9 shipped copies; GSTM3 seed, 4,057 bp < β → 0 loci vs 4 shipped copies; HERC2 → 164 vs 2).
//!
//! The shipped binary's node set is built from the reads (or, under `--from-genome`, from genome
//! self-alignment) and never consults a seed. `gamma_quasi_clique_partition` starts from
//! `all_components` INCLUDING singletons, so the blocks **partition** the whole node set. Over a
//! seed-free node set, "the block containing `s`" is therefore a fact about a partition, not a
//! parameterised computation: every node lands in exactly one block, and any two nodes in the same
//! block return the same block. That is what makes `--seed` a legitimate query — and it is why the
//! query must read the emitted catalog rather than re-run anything with `s` in scope.
//!
//! ## Scope limit, stated up front
//!
//! The projection can only see blocks the run actually EMITTED. At the default `--min-copies 2` a
//! singleton block is not written, so a seed landing on one is reported as `ABSTAIN_NO_OVERLAP` and
//! is indistinguishable here from a seed landing on no transcribed locus at all. Re-run with
//! `--min-copies 1` to separate the two cases. Abstention is the O1-safe answer: "no emitted family
//! at this seed" is a legitimate result and is never guessed around.

use anyhow::{anyhow, Result};

/// A seed interval, stored 0-based half-open to match `copies.tsv`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeedLocus {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    /// The spec exactly as the user typed it, echoed into the output so a row is traceable.
    pub spec: String,
}

/// Outcome of projecting one seed onto the emitted catalog.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SeedHit {
    /// The seed overlaps `copy_idx` of family `family_idx` by `overlap_bp` bases.
    Hit {
        family_idx: usize,
        copy_idx: usize,
        overlap_bp: u64,
    },
    /// No emitted copy overlaps the seed. Honest answer; not an error.
    Abstain,
}

/// Parse `chrom:start-end`.
///
/// Coordinates are **1-based inclusive** (the GFF/samtools region convention, which is what a seed
/// pasted out of an annotation carries) and are converted to 0-based half-open on the way in, so
/// they line up with `copies.tsv`. `,` and `_` are accepted as digit separators. The chromosome
/// name may itself contain `:` (none of ours do, but RefSeq-style names contain `.` and `_`), so the
/// split is on the LAST `:`.
pub fn parse_seed(spec: &str) -> Result<SeedLocus> {
    let s = spec.trim();
    let (chrom, range) = s
        .rsplit_once(':')
        .ok_or_else(|| anyhow!("--seed {spec:?}: expected chrom:start-end (1-based inclusive)"))?;
    if chrom.is_empty() {
        return Err(anyhow!("--seed {spec:?}: empty chromosome"));
    }
    let (a, b) = range
        .rsplit_once('-')
        .ok_or_else(|| anyhow!("--seed {spec:?}: expected chrom:start-end (1-based inclusive)"))?;
    let clean = |x: &str| x.replace([',', '_'], "");
    let a: u64 = clean(a)
        .parse()
        .map_err(|_| anyhow!("--seed {spec:?}: start is not an integer"))?;
    let b: u64 = clean(b)
        .parse()
        .map_err(|_| anyhow!("--seed {spec:?}: end is not an integer"))?;
    if a == 0 {
        return Err(anyhow!("--seed {spec:?}: coordinates are 1-based, so start must be >= 1"));
    }
    if b < a {
        return Err(anyhow!("--seed {spec:?}: end < start"));
    }
    Ok(SeedLocus {
        chrom: chrom.to_string(),
        start: a - 1, // 1-based inclusive -> 0-based half-open
        end: b,
        spec: s.to_string(),
    })
}

/// Overlap of two half-open intervals, in bp (0 when disjoint or touching).
pub fn overlap_bp(a0: u64, a1: u64, b0: u64, b1: u64) -> u64 {
    let lo = a0.max(b0);
    let hi = a1.min(b1);
    hi.saturating_sub(lo)
}

/// Project a seed onto the emitted catalog.
///
/// `fams[fi][ci]` must be `(chrom, start, end)` in **exactly the order `emit_catalog` writes them** —
/// families in emitted order (so `fi` ↔ `GWFAM{fi}`), copies sorted by `(chrom, start)` (so `ci` ↔
/// `copy_idx`) — otherwise the reported ids do not name the rows the user can look up.
///
/// Rule: maximum overlapping bases wins. This is deliberate and sufficient rather than clever — on
/// the gorilla RABL2 check the shipped copy boundaries sat within 3–46 bp of the independently
/// seeded interval, so no seed is anywhere near two copies at once. Ties (including the degenerate
/// zero-length seed) resolve to the lowest `(family_idx, copy_idx)`, which is deterministic because
/// the caller's order is deterministic.
pub fn project_seed(seed: &SeedLocus, fams: &[Vec<(String, u64, u64)>]) -> SeedHit {
    let mut best: Option<(u64, usize, usize)> = None;
    for (fi, copies) in fams.iter().enumerate() {
        for (ci, (chrom, start, end)) in copies.iter().enumerate() {
            if chrom != &seed.chrom {
                continue;
            }
            let ov = overlap_bp(seed.start, seed.end, *start, *end);
            if ov == 0 {
                continue;
            }
            if best.map_or(true, |(b, _, _)| ov > b) {
                best = Some((ov, fi, ci));
            }
        }
    }
    match best {
        Some((overlap_bp, family_idx, copy_idx)) => SeedHit::Hit {
            family_idx,
            copy_idx,
            overlap_bp,
        },
        None => SeedHit::Abstain,
    }
}

/// Header of `<out>.seed.tsv`.
pub const SEED_TSV_HEADER: &str =
    "seed\tstatus\tfamily_id\tn_copies\tmember_idx\tchrom\tstart\tend\tis_seed_locus\tseed_overlap_bp";

/// Render the projection of one seed as the rows of `<out>.seed.tsv`: one row per MEMBER of the
/// component containing the seed (that is the object the query returns — the component, not the
/// single locus the seed happened to land on), or exactly one `ABSTAIN_NO_OVERLAP` row.
pub fn format_seed_rows(seed: &SeedLocus, fams: &[Vec<(String, u64, u64)>], hit: &SeedHit) -> Vec<String> {
    match hit {
        SeedHit::Abstain => vec![format!(
            "{}\tABSTAIN_NO_OVERLAP\t.\t0\t.\t.\t.\t.\t.\t0",
            seed.spec
        )],
        SeedHit::Hit {
            family_idx,
            copy_idx,
            overlap_bp,
        } => {
            let copies = &fams[*family_idx];
            copies
                .iter()
                .enumerate()
                .map(|(ci, (chrom, start, end))| {
                    format!(
                        "{}\tHIT\tGWFAM{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        seed.spec,
                        family_idx,
                        copies.len(),
                        ci,
                        chrom,
                        start,
                        end,
                        if ci == *copy_idx { "true" } else { "false" },
                        if ci == *copy_idx { *overlap_bp } else { 0 },
                    )
                })
                .collect()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fams() -> Vec<Vec<(String, u64, u64)>> {
        vec![
            // GWFAM0: two copies on chr1
            vec![
                ("chr1".to_string(), 1000, 2000),
                ("chr1".to_string(), 5000, 6000),
            ],
            // GWFAM1: two copies, one of them adjacent-but-disjoint to GWFAM0's second copy
            vec![
                ("chr1".to_string(), 6000, 7000),
                ("chr2".to_string(), 100, 900),
            ],
        ]
    }

    // ---- parse_seed ----

    #[test]
    fn parse_converts_one_based_inclusive_to_half_open() {
        let s = parse_seed("chr1:1001-2000").unwrap();
        assert_eq!(s.chrom, "chr1");
        assert_eq!((s.start, s.end), (1000, 2000));
        assert_eq!(s.spec, "chr1:1001-2000");
    }

    #[test]
    fn parse_accepts_refseq_style_names_and_digit_separators() {
        let s = parse_seed("NC_073236.2:140,592,916-140_608_794").unwrap();
        assert_eq!(s.chrom, "NC_073236.2");
        assert_eq!((s.start, s.end), (140_592_915, 140_608_794));
    }

    #[test]
    fn parse_rejects_zero_start_because_coordinates_are_one_based() {
        assert!(parse_seed("chr1:0-100").is_err());
    }

    #[test]
    fn parse_rejects_malformed_specs() {
        assert!(parse_seed("chr1").is_err());
        assert!(parse_seed("chr1:100").is_err());
        assert!(parse_seed("chr1:abc-100").is_err());
        assert!(parse_seed("chr1:200-100").is_err());
        assert!(parse_seed(":100-200").is_err());
    }

    #[test]
    fn single_base_seed_is_legal_and_has_length_one() {
        let s = parse_seed("chr1:1500-1500").unwrap();
        assert_eq!((s.start, s.end), (1499, 1500));
    }

    // ---- overlap ----

    #[test]
    fn touching_intervals_do_not_overlap() {
        assert_eq!(overlap_bp(0, 100, 100, 200), 0);
        assert_eq!(overlap_bp(0, 101, 100, 200), 1);
    }

    // ---- project_seed ----

    #[test]
    fn seed_inside_a_copy_returns_that_copy() {
        let f = fams();
        let s = parse_seed("chr1:1201-1300").unwrap();
        assert_eq!(
            project_seed(&s, &f),
            SeedHit::Hit { family_idx: 0, copy_idx: 0, overlap_bp: 100 }
        );
    }

    #[test]
    fn seed_off_every_copy_abstains_rather_than_guessing() {
        let f = fams();
        let s = parse_seed("chr1:3001-4000").unwrap();
        assert_eq!(project_seed(&s, &f), SeedHit::Abstain);
    }

    #[test]
    fn seed_on_an_absent_chromosome_abstains() {
        let f = fams();
        let s = parse_seed("chrX:1001-2000").unwrap();
        assert_eq!(project_seed(&s, &f), SeedHit::Abstain);
    }

    #[test]
    fn straddling_seed_goes_to_the_copy_it_overlaps_most_not_the_first_one() {
        // 5900..6600 overlaps GWFAM0 copy1 (5000..6000) by 100 bp and GWFAM1 copy0 (6000..7000) by
        // 600 bp. Max-overlap must win over first-seen, or the answer would depend on emit order.
        let f = fams();
        let s = parse_seed("chr1:5901-6600").unwrap();
        assert_eq!(
            project_seed(&s, &f),
            SeedHit::Hit { family_idx: 1, copy_idx: 0, overlap_bp: 600 }
        );
    }

    #[test]
    fn exact_ties_resolve_to_the_lowest_family_then_copy_index() {
        let f = vec![
            vec![("chr1".to_string(), 0, 100)],
            vec![("chr1".to_string(), 0, 100)],
        ];
        let s = parse_seed("chr1:1-100").unwrap();
        assert_eq!(
            project_seed(&s, &f),
            SeedHit::Hit { family_idx: 0, copy_idx: 0, overlap_bp: 100 }
        );
    }

    #[test]
    fn empty_catalog_abstains() {
        let s = parse_seed("chr1:1-100").unwrap();
        assert_eq!(project_seed(&s, &[]), SeedHit::Abstain);
    }

    // ---- the returned object is the COMPONENT, not the locus ----

    #[test]
    fn hit_rows_carry_every_member_of_the_component() {
        let f = fams();
        let s = parse_seed("chr1:1201-1300").unwrap();
        let hit = project_seed(&s, &f);
        let rows = format_seed_rows(&s, &f, &hit);
        assert_eq!(rows.len(), 2, "the query returns the whole component: {rows:?}");
        assert!(rows.iter().all(|r| r.contains("\tGWFAM0\t")));
        assert_eq!(rows.iter().filter(|r| r.ends_with("\ttrue\t100")).count(), 1);
        assert!(rows[1].contains("\tfalse\t0"), "non-seed members carry overlap 0: {}", rows[1]);
    }

    #[test]
    fn abstain_emits_exactly_one_row_and_names_the_reason() {
        let f = fams();
        let s = parse_seed("chr9:1-100").unwrap();
        let rows = format_seed_rows(&s, &f, &project_seed(&s, &f));
        assert_eq!(rows.len(), 1);
        assert!(rows[0].contains("ABSTAIN_NO_OVERLAP"), "{}", rows[0]);
        assert!(rows[0].starts_with("chr9:1-100\t"));
    }

    #[test]
    fn every_row_has_the_same_field_count_as_the_header() {
        let f = fams();
        let ncol = SEED_TSV_HEADER.split('\t').count();
        for spec in ["chr1:1201-1300", "chr9:1-100"] {
            let s = parse_seed(spec).unwrap();
            for r in format_seed_rows(&s, &f, &project_seed(&s, &f)) {
                assert_eq!(r.split('\t').count(), ncol, "{r}");
            }
        }
    }

    // ---- the property the projection inherits from the partition ----

    #[test]
    fn any_member_of_a_component_returns_the_same_component() {
        // P1 is a THEOREM over a seed-free node set: blocks partition the nodes, so seeding from any
        // member returns the identical block. Asserted here so a future change that makes the answer
        // seed-dependent fails a test rather than a re-measurement.
        let f = fams();
        let seeds = ["chr1:1001-2000", "chr1:5001-6000"];
        let rows: Vec<Vec<String>> = seeds
            .iter()
            .map(|spec| {
                let s = parse_seed(spec).unwrap();
                // strip the echoed seed spec + per-seed overlap columns; compare the COMPONENT.
                format_seed_rows(&s, &f, &project_seed(&s, &f))
                    .iter()
                    .map(|r| {
                        let c: Vec<&str> = r.split('\t').collect();
                        c[2..8].join("\t")
                    })
                    .collect()
            })
            .collect();
        assert_eq!(rows[0], rows[1]);
    }
}
