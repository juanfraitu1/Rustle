//! r1600 task 4: `RUSTLE_ASSIGN_ROUTER=1` (one routed pass for every window's light records, the tied gate, one routed
//! pass for the tied molecules' full records) writes every output byte for byte as the per-window query path, at 1 and
//! 4 region threads, on the `--tied-prefilter` path the pipeline runs.
use std::path::PathBuf;
use std::process::Command;

fn fx(name: &str) -> String {
    format!("{}/tests/fixtures/{name}", env!("CARGO_MANIFEST_DIR"))
}

/// Every `OUT.*` file of a run, name -> bytes.
fn outputs(dir: &PathBuf) -> std::collections::BTreeMap<String, Vec<u8>> {
    std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.file_name().unwrap().to_string_lossy().starts_with("o."))
        .map(|p| (p.file_name().unwrap().to_string_lossy().into_owned(), std::fs::read(&p).unwrap()))
        .collect()
}

fn run(fixture: &str, regions: &str, router: bool, region_threads: &str) -> (std::collections::BTreeMap<String, Vec<u8>>, String) {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR"))
        .join("copy_assign_router")
        .join(format!("{fixture}_{router}_{region_threads}"));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    let rf = d.join("regions.txt");
    std::fs::write(&rf, regions).unwrap();
    let f = fx(fixture);
    let mut c = Command::new(env!("CARGO_BIN_EXE_copy_assign"));
    c.args(["--bam", &format!("{f}/reads.bam"), "--fasta", &format!("{f}/genome.fa")])
        .args(["--families", &format!("{f}/copies.tsv"), "--copies-fa", &format!("{f}/copies.fa")])
        .args(["--regions", rf.to_str().unwrap(), "--region-threads", region_threads, "--tied-prefilter"])
        .args(["--out", d.join("o").to_str().unwrap()])
        .env("RUSTLE_TIMING", "1")
        .env_remove("RUSTLE_WINDOW_SWEEP");
    if router {
        c.env("RUSTLE_ASSIGN_ROUTER", "1");
    } else {
        c.env_remove("RUSTLE_ASSIGN_ROUTER");
    }
    let o = c.output().expect("spawn copy_assign");
    let err = String::from_utf8_lossy(&o.stderr).into_owned();
    assert!(o.status.success(), "{fixture} router={router} rt={region_threads} failed:\n{err}");
    (outputs(&d), err)
}

/// `[timing] reads_in_region REGION (N reads)` and `... light ...` counts, in order: what each region gathered.
fn read_counts(err: &str) -> Vec<String> {
    // reads_in_region / gathered lines end in a wall time (": 0.1s"), cut off; a light-records line is all count
    err.lines()
        .filter(|l| l.contains("reads_in_region") || l.contains("gathered from") || l.contains("light records"))
        .map(|l| if l.contains("light records") { l.to_string() } else { l.split(": ").next().unwrap_or("").to_string() })
        .collect()
}

#[test]
fn router_arm_is_byte_identical() {
    for (fixture, regions) in [("o2_origin_resolve", "c1:0-3500\n"), ("xfam_conflict", "x1:0-20000\n")] {
        let (base, base_err) = run(fixture, regions, false, "1");
        let rows = String::from_utf8_lossy(&base["o.assignments.tsv"]).lines().count();
        assert!(rows > 1, "{fixture}: the base arm assigned nothing ({rows} lines)");
        for (router, rt) in [(false, "4"), (true, "1"), (true, "4")] {
            let (got, err) = run(fixture, regions, router, rt);
            assert_eq!(got.keys().collect::<Vec<_>>(), base.keys().collect::<Vec<_>>(), "{fixture} router={router} rt={rt}: file set");
            for (k, v) in &base {
                assert!(got[k] == *v, "{fixture} router={router} rt={rt}: {k} differs");
            }
            assert_eq!(err.contains("[timing] router pass2"), router, "{fixture} router={router} rt={rt}: router path taken?");
            assert_eq!(read_counts(&err), read_counts(&base_err), "{fixture} router={router} rt={rt}: per-region read counts");
        }
    }
}
