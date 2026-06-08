use rustle::family_manifest::{parse_family_manifest, create_family_groups_from_manifest, FamilyLocus};
use rustle::types::{Bundle, JunctionStats};
use std::io::Write;

fn mk_bundle(chrom: &str, start: u64, end: u64) -> Bundle {
    Bundle {
        chrom: chrom.into(),
        start,
        end,
        strand: '+',
        reads: Vec::new(),
        junction_stats: JunctionStats::default(),
        junction_pair_stats: Default::default(),
        bundlenodes: None,
        read_bnodes: None,
        bnode_colors: None,
        synthetic: false,
        rescue_class: None,
        vg_family_id: None,
        hp_tag: None,
        ps_tag: None,
    }
}

#[test]
fn parse_manifest_two_families() {
    let content = "family_id\tgene_id\tchrom\tstart\tend\tstrand\n\
                   FAM1\tgeneA\tchr1\t1000\t2000\t+\n\
                   FAM1\tgeneB\tchr2\t5000\t6000\t+\n\
                   FAM2\tgeneC\tchr3\t9000\t10000\t-\n";
    let mut f = tempfile::NamedTempFile::new().unwrap();
    write!(f, "{}", content).unwrap();
    let loci = parse_family_manifest(f.path()).unwrap();
    assert_eq!(loci.len(), 3);
    assert_eq!(loci[0].family_id, "FAM1");
    assert_eq!(loci[0].gene_id, "geneA");
    assert_eq!(loci[0].chrom, "chr1");
    assert_eq!(loci[0].start, 1000);
    assert_eq!(loci[0].end, 2000);
    assert_eq!(loci[0].strand, '+');
    assert_eq!(loci[2].family_id, "FAM2");
    assert_eq!(loci[2].strand, '-');
}

#[test]
fn parse_manifest_skips_blank_and_comment_lines() {
    let content = "family_id\tgene_id\tchrom\tstart\tend\tstrand\n\
                   \n\
                   # comment\n\
                   FAM1\tgeneA\tchr1\t1000\t2000\t+\n";
    let mut f = tempfile::NamedTempFile::new().unwrap();
    write!(f, "{}", content).unwrap();
    let loci = parse_family_manifest(f.path()).unwrap();
    assert_eq!(loci.len(), 1);
}

#[test]
fn create_groups_links_overlapping_bundles() {
    let loci = vec![
        FamilyLocus { family_id: "FAM1".into(), gene_id: "gA".into(), chrom: "chr1".into(), start: 1000, end: 2000, strand: '+' },
        FamilyLocus { family_id: "FAM1".into(), gene_id: "gB".into(), chrom: "chr2".into(), start: 5000, end: 6000, strand: '+' },
    ];
    let bundles = vec![
        mk_bundle("chr1", 1200, 1800),   // overlaps gA locus
        mk_bundle("chr2", 5100, 5900),   // overlaps gB locus
        mk_bundle("chr3", 0, 1000),      // no match
    ];
    let groups = create_family_groups_from_manifest(&loci, &bundles);
    assert_eq!(groups.len(), 1, "should produce exactly one family group");
    let g = &groups[0];
    let mut idxs = g.bundle_indices.clone();
    idxs.sort();
    assert_eq!(idxs, vec![0, 1]);
}

#[test]
fn create_groups_two_families() {
    let loci = vec![
        FamilyLocus { family_id: "FAM1".into(), gene_id: "gA".into(), chrom: "chr1".into(), start: 1000, end: 2000, strand: '+' },
        FamilyLocus { family_id: "FAM2".into(), gene_id: "gC".into(), chrom: "chr3".into(), start: 9000, end: 10000, strand: '-' },
    ];
    let bundles = vec![
        mk_bundle("chr1", 1500, 1800),
        mk_bundle("chr3", 9500, 9900),
    ];
    let groups = create_family_groups_from_manifest(&loci, &bundles);
    assert_eq!(groups.len(), 2);
}

#[test]
fn create_groups_empty_when_no_overlap() {
    let loci = vec![
        FamilyLocus { family_id: "FAM1".into(), gene_id: "gA".into(), chrom: "chr1".into(), start: 1000, end: 2000, strand: '+' },
    ];
    let bundles = vec![mk_bundle("chr1", 5000, 6000)]; // no overlap
    let groups = create_family_groups_from_manifest(&loci, &bundles);
    assert!(groups.is_empty());
}
