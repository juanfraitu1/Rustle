//! `--seed` is a QUERY over the emitted catalog, not an input to the definition.
//!
//! The load-bearing test here is `seed_does_not_change_the_catalog`: if the seed could move any
//! other output file, then "the seed is a query" would be false and `F(s)` would be a
//! seed-parameterised computation again — which is exactly the defect measured on the
//! `bench/crossspecies` seeded probe, whose node set `V(s)` inherits the seed's own length.

use std::process::Command;

const BAM: &str = "tests/fixtures/same_chrom_supplement/reads.bam";
const FASTA: &str = "tests/fixtures/same_chrom_supplement/genome.fa";

fn run(out: &str, extra: &[&str]) -> std::process::Output {
    let bin = env!("CARGO_BIN_EXE_gw_family_catalog");
    let mut args: Vec<&str> = vec!["--bam", BAM, "--fasta", FASTA, "--out", out];
    args.extend_from_slice(extra);
    Command::new(bin).args(&args).output().unwrap()
}

fn tmp(name: &str) -> String {
    std::env::temp_dir().join(name).to_str().unwrap().to_string()
}

/// The fixture emits one family, GWFAM0, with four copies: c1:0-260, c1:250-460, c1:380-550,
/// c2:0-260 (0-based half-open, as in copies.tsv).
#[test]
fn seed_inside_a_copy_returns_the_whole_component() {
    let out = tmp("seedproj_hit");
    // 1-based inclusive c1:1-260 == 0-based half-open c1:0-260 == copy 0 exactly.
    let o = run(&out, &["--seed", "c1:1-260"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let tsv = std::fs::read_to_string(format!("{out}.seed.tsv")).unwrap();
    let rows: Vec<&str> = tsv.lines().skip(1).collect();
    assert_eq!(
        rows.len(),
        4,
        "the query returns the COMPONENT (all 4 copies), not just the locus the seed hit:\n{tsv}"
    );
    assert!(rows.iter().all(|r| r.contains("\tHIT\tGWFAM0\t")), "{tsv}");
    // exactly one member is flagged as the one the seed landed on
    let seeded: Vec<&&str> = rows.iter().filter(|r| r.split('\t').nth(8) == Some("true")).collect();
    assert_eq!(seeded.len(), 1, "{tsv}");
    let f: Vec<&str> = seeded[0].split('\t').collect();
    assert_eq!((f[5], f[6], f[7]), ("c1", "0", "260"), "{tsv}");
    assert_eq!(f[9], "260", "overlap should be the full 260 bp: {tsv}");
}

#[test]
fn straddling_seed_goes_to_the_copy_it_overlaps_most() {
    let out = tmp("seedproj_straddle");
    // 1-based 400-560 == 0-based 399..560: overlaps c1:380-550 by 151 and c1:250-460 by 61.
    let o = run(&out, &["--seed", "c1:400-560"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let tsv = std::fs::read_to_string(format!("{out}.seed.tsv")).unwrap();
    let seeded: Vec<&str> = tsv
        .lines()
        .skip(1)
        .filter(|r| r.split('\t').nth(8) == Some("true"))
        .collect();
    assert_eq!(seeded.len(), 1, "{tsv}");
    let f: Vec<&str> = seeded[0].split('\t').collect();
    assert_eq!((f[5], f[6], f[7], f[9]), ("c1", "380", "550", "151"), "{tsv}");
}

#[test]
fn seed_off_every_copy_abstains_and_does_not_snap_to_a_neighbour() {
    let out = tmp("seedproj_abstain");
    let o = run(&out, &["--seed", "c1:900-1000", "--seed", "cZ:1-100"]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let tsv = std::fs::read_to_string(format!("{out}.seed.tsv")).unwrap();
    let rows: Vec<&str> = tsv.lines().skip(1).collect();
    assert_eq!(rows.len(), 2, "one ABSTAIN row per seed, no guessing:\n{tsv}");
    assert!(rows.iter().all(|r| r.contains("ABSTAIN_NO_OVERLAP")), "{tsv}");
    assert!(rows[0].starts_with("c1:900-1000\t") && rows[1].starts_with("cZ:1-100\t"), "{tsv}");
    // the operator is told what an abstention does and does not mean
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("--min-copies 1"), "abstention must name the singleton caveat:\n{err}");
}

/// THE property. `--seed` may add `<out>.seed.tsv` and nothing else.
#[test]
fn seed_does_not_change_the_catalog() {
    let plain = tmp("seedproj_plain");
    let seeded = tmp("seedproj_seeded");
    assert!(run(&plain, &[]).status.success());
    assert!(run(&seeded, &["--seed", "c1:1-260"]).status.success());
    for ext in ["families.tsv", "copies.tsv", "copies.fa"] {
        let a = std::fs::read(format!("{plain}.{ext}")).unwrap();
        let b = std::fs::read(format!("{seeded}.{ext}")).unwrap();
        assert_eq!(a, b, "--seed changed {ext}; it must be a projection, not an input");
    }
    assert!(
        std::fs::metadata(format!("{plain}.seed.tsv")).is_err(),
        "no --seed must mean no .seed.tsv (byte-identical output set)"
    );
}

/// Seed-invariance is a THEOREM once the node set is seed-free (the blocks partition the nodes), and
/// this asserts it end-to-end: seeding from any member of a component returns that same component.
#[test]
fn any_member_of_a_component_returns_the_same_component() {
    let out = tmp("seedproj_invariance");
    let o = run(
        &out,
        &["--seed", "c1:1-260", "--seed", "c1:251-460", "--seed", "c2:1-260"],
    );
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let tsv = std::fs::read_to_string(format!("{out}.seed.tsv")).unwrap();
    // group rows by the seed spec, compare the member list (cols 2..8: family_id..end)
    use std::collections::BTreeMap;
    let mut by_seed: BTreeMap<&str, Vec<String>> = BTreeMap::new();
    for l in tsv.lines().skip(1) {
        let c: Vec<&str> = l.split('\t').collect();
        by_seed.entry(c[0]).or_default().push(c[2..8].join("\t"));
    }
    assert_eq!(by_seed.len(), 3, "{tsv}");
    let sets: Vec<&Vec<String>> = by_seed.values().collect();
    assert_eq!(sets[0], sets[1], "F(s) must not depend on which member seeded it:\n{tsv}");
    assert_eq!(sets[0], sets[2], "F(s) must not depend on which member seeded it:\n{tsv}");
}

#[test]
fn a_malformed_seed_fails_loudly_instead_of_silently_abstaining() {
    let out = tmp("seedproj_bad");
    let o = run(&out, &["--seed", "c1:not-a-number"]);
    assert!(!o.status.success(), "a typo must not be reported as 'no family here'");
    assert!(
        std::fs::metadata(format!("{out}.seed.tsv")).is_err(),
        "a rejected seed must not leave a partial .seed.tsv"
    );
}
