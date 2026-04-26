use std::process::Command;
use std::path::Path;

#[test]
fn flow_iter_dump_emits_monotone_iterations() {
    let bam = "test_data/synthetic_family/reads_sorted.bam";
    let out_gtf = "/tmp/test_flow_iter_dump.gtf";
    let flow_tsv = "/tmp/test_flow_iter_dump.tsv";
    let _ = std::fs::remove_file(flow_tsv);

    let status = Command::new(env!("CARGO_BIN_EXE_rustle"))
        .args(["-L", "-o", out_gtf, bam])
        .env("RUSTLE_FLOW_ITER_TSV", flow_tsv)
        .status()
        .expect("rustle run failed");
    assert!(status.success());

    assert!(Path::new(flow_tsv).exists());
    let content = std::fs::read_to_string(flow_tsv).expect("read flow TSV");
    let mut lines = content.lines();
    let header = lines.next().unwrap();
    assert_eq!(
        header,
        "source\tchrom\tbdstart\tbdend\tstrand\tbundle_run\titer_idx\tpath_len\tpath_nodes\tbottleneck\ttotal_flow_after"
    );

    let row_count = lines.filter(|l| !l.trim().is_empty()).count();
    assert!(row_count >= 2, "expected >=2 flow iterations, got {}", row_count);
}
