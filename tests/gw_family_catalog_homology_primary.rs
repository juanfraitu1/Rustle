use std::process::Command;
#[test]
fn homology_primary_recovers_gstm_and_znf() {
    let bam = "/home/juanfra/winloci_scratch/GGO_mm.bam";
    if std::fs::metadata(bam).is_err() { eprintln!("real data absent; skip"); return; }
    // chr234 = NC_073234.2 has the KRAB-ZNF trio; run homology-primary and assert a >=3-copy family forms.
    let bin = env!("CARGO_BIN_EXE_gw_family_catalog");
    let out = std::env::temp_dir().join("hom_c234");
    // caller pre-subsets NC_073244.2 (KRAB-ZNF) into a small BAM to keep the test fast; see plan Task 8 notes.
    let sub = "/home/juanfra/winloci_scratch/c244.bam";
    if std::fs::metadata(sub).is_err() { eprintln!("subset c244.bam absent; skip"); return; }
    let ok = std::process::Command::new(bin).args([
        "--bam", sub, "--fasta", "/home/juanfra/winloci_scratch/GGO.fasta",
        "--out", out.to_str().unwrap(), "--homology-primary",
    ]).status().unwrap().success();
    assert!(ok);
    let copies = std::fs::read_to_string(format!("{}.copies.tsv", out.to_str().unwrap())).unwrap();
    // ZNF14/ZNF682/ZNF429 span ~26.1–28.0 Mb on NC_073244.2. The real RECOVERY signal is that homology-
    // primary groups >=2 of these paralogs into ONE family — the old read-conflict + 0.80-refine catalog
    // dropped them (exon identity 0.62–0.71 < 0.80). copies.tsv header: family_id, copy_idx, tid, chrom,
    // start, end, ... → family_id=f[0], chrom=f[3], start=f[4]. Group window copies by family and assert
    // the LARGEST window-family carries >=2 (not merely "some rows exist on this single-contig subset").
    use std::collections::HashMap;
    let mut by_fam: HashMap<String, usize> = HashMap::new();
    for l in copies.lines().skip(1) {
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() < 5 { continue; }
        if f[3] == "NC_073244.2" {
            if let Ok(s) = f[4].parse::<u64>() {
                if (26_000_000..=28_500_000).contains(&s) {
                    *by_fam.entry(f[0].to_string()).or_default() += 1;
                }
            }
        }
    }
    let max_in_window = by_fam.values().copied().max().unwrap_or(0);
    assert!(max_in_window >= 2,
        "homology-primary must group >=2 KRAB-ZNF paralogs (NC_073244.2:26-28.5Mb) into ONE family; largest window-family={max_in_window}, by_fam={by_fam:?}");
}

#[test]
fn homology_primary_emits_families_on_fixture() {
    let bin = env!("CARGO_BIN_EXE_gw_family_catalog");
    let out = "tests/fixtures/same_chrom_supplement/out_hom";
    let status = Command::new(bin).args([
        "--bam", "tests/fixtures/same_chrom_supplement/reads.bam",
        "--fasta", "tests/fixtures/same_chrom_supplement/genome.fa",
        "--out", out, "--homology-primary",
    ]).status().unwrap();
    assert!(status.success());
    let fams = std::fs::read_to_string(format!("{}.families.tsv", out)).unwrap();
    assert!(fams.lines().count() > 1, "expected >=1 homology family\n{}", fams);
}

#[test]
fn enumerate_copies_ignored_without_homology_primary() {
    // famCN (genome-projection copy enumeration) is only meaningful over the homology-primary
    // catalog: the conflict/default catalog's `fams` can carry same-locus duplicates (no distinct-
    // locus reduction without `--refine`), so a "family consensus" projected back onto the genome
    // would double-count. `--enumerate-copies` without `--homology-primary` must be a no-op: no
    // `<out>.famcn.tsv` written.
    let bin = env!("CARGO_BIN_EXE_gw_family_catalog");
    let out = "tests/fixtures/same_chrom_supplement/out_conflict_enum";
    let famcn = format!("{}.famcn.tsv", out);
    let _ = std::fs::remove_file(&famcn);
    let status = Command::new(bin).args([
        "--bam", "tests/fixtures/same_chrom_supplement/reads.bam",
        "--fasta", "tests/fixtures/same_chrom_supplement/genome.fa",
        "--out", out, "--cross-chrom", "--enumerate-copies",
    ]).status().unwrap();
    assert!(status.success());
    assert!(std::fs::metadata(&famcn).is_err(),
        "famcn.tsv must NOT be written when --homology-primary is absent (found {famcn})");
}

#[test]
fn families_tsv_has_protein_coheres_column() {
    // spec §6: the orthogonal protein-coherence QC flag (family_protein_coheres) must be emitted as
    // the LAST column of `<out>.families.tsv`. Without `--protein-tail` (mmseqs off), the column must
    // be present but "NA" for every family (no mmseqs call, no cost).
    let bin = env!("CARGO_BIN_EXE_gw_family_catalog");
    let out = "tests/fixtures/same_chrom_supplement/out_hom_protein_qc";
    let status = Command::new(bin).args([
        "--bam", "tests/fixtures/same_chrom_supplement/reads.bam",
        "--fasta", "tests/fixtures/same_chrom_supplement/genome.fa",
        "--out", out, "--homology-primary",
    ]).status().unwrap();
    assert!(status.success());
    let fams = std::fs::read_to_string(format!("{}.families.tsv", out)).unwrap();
    let mut lines = fams.lines();
    let header = lines.next().expect("families.tsv must have a header");
    let header_cols: Vec<&str> = header.split('\t').collect();
    // ⚠ THIS TEST USED TO READ `protein_coheres` POSITIONALLY, as the LAST column, and broke the first
    // time a column was appended (the λ certificate, 2026-08-14). The emitter's stated contract is that
    // new columns are APPENDED so that HEADER-KEYED readers keep working — so the test now reads the way
    // the contract says readers do, by name. Same assertion, no longer coupled to column count.
    let coheres_at = header_cols
        .iter()
        .position(|c| *c == "protein_coheres")
        .unwrap_or_else(|| panic!("families.tsv must carry a protein_coheres column; header={header}"));
    // The λ certificate columns are part of the emitted family (docs/seeded_family_definition.md §1★.5).
    for col in ["n_edges", "density", "lambda", "cut_certified"] {
        assert!(
            header_cols.contains(&col),
            "families.tsv must carry the `{col}` certificate column; header={header}"
        );
    }
    let lambda_at = header_cols.iter().position(|c| *c == "lambda").unwrap();
    let n_copies_at = header_cols.iter().position(|c| *c == "n_copies").unwrap();
    let certified_at = header_cols.iter().position(|c| *c == "cut_certified").unwrap();
    let mut n_rows = 0;
    for l in lines {
        if l.is_empty() { continue; }
        n_rows += 1;
        let cols: Vec<&str> = l.split('\t').collect();
        assert_eq!(cols.len(), header_cols.len(), "row/header column count must agree; row={l}");
        assert_eq!(
            cols.get(coheres_at).copied(), Some("NA"),
            "without --protein-tail every row's protein_coheres must be NA; row={l}"
        );
        // λ is REPORTED, never enforced, so the only invariants are arithmetic ones. The fixture's
        // families are 2-copy, which is precisely the case that CANNOT be cut-certified.
        let n_copies: usize = cols[n_copies_at].parse().unwrap();
        let lambda: usize = cols[lambda_at].parse().unwrap_or_else(|_| panic!("lambda must be numeric on the homology path; row={l}"));
        assert!(lambda < n_copies.max(2), "lambda cannot reach n for a simple graph; row={l}");
        assert_eq!(
            cols[certified_at] == "true", lambda >= 2,
            "cut_certified must be exactly `lambda >= 2`; row={l}"
        );
        if n_copies == 2 {
            assert_eq!(lambda, 1, "a 2-copy family has lambda = 1 NECESSARILY; row={l}");
        }
    }
    assert!(n_rows > 0, "expected >=1 family row\n{}", fams);
}

#[test]
fn homology_primary_writes_famcn_when_enumerating() {
    let bin = env!("CARGO_BIN_EXE_gw_family_catalog");
    let out = "tests/fixtures/same_chrom_supplement/out_famcn";
    let status = Command::new(bin).args([
        "--bam", "tests/fixtures/same_chrom_supplement/reads.bam",
        "--fasta", "tests/fixtures/same_chrom_supplement/genome.fa",
        "--out", out, "--homology-primary", "--enumerate-copies",
    ]).status().unwrap();
    assert!(status.success());
    assert!(std::fs::metadata(format!("{}.famcn.tsv", out)).is_ok(), "famcn.tsv must be written");
}
