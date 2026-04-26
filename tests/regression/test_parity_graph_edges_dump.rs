use std::process::Command;
use std::path::Path;

#[test]
fn graph_edges_dump_produces_expected_schema() {
    let bam = "test_data/synthetic_family/reads_sorted.bam";
    let out_gtf = "/tmp/test_graph_edges_dump.gtf";
    let edges_tsv = "/tmp/test_graph_edges_dump.tsv";
    let _ = std::fs::remove_file(edges_tsv);

    let status = Command::new(env!("CARGO_BIN_EXE_rustle"))
        .args(["-L", "-o", out_gtf, bam])
        .env("RUSTLE_PARITY_GRAPH_EDGES_TSV", edges_tsv)
        .status()
        .expect("rustle run failed");
    assert!(status.success(), "rustle exited non-zero");

    assert!(Path::new(edges_tsv).exists(), "edges TSV not created");
    let content = std::fs::read_to_string(edges_tsv).expect("read edges TSV");
    let mut lines = content.lines();

    let header = lines.next().expect("missing header");
    let expected_cols = "source\tchrom\tbdstart\tbdend\tstrand\tfrom_idx\tto_idx\tfrom_start\tfrom_end\tto_start\tto_end\tedge_kind\tfrom_cov\tto_cov\tbottleneck_cov";
    assert_eq!(header, expected_cols, "header mismatch");

    let body_count = lines.filter(|l| !l.trim().is_empty()).count();
    assert!(body_count > 0, "no edge rows written");
}
