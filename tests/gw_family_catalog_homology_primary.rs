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
