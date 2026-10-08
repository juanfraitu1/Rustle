//! `tools/rustle_pipeline.sh` with `RUSTLE_BRIDGE_REGROUP=f1units`, `RUSTLE_BRIDGE_UNITS_LIST` and
//! `RUSTLE_FAMILY_RELATIONS=1` (`docs/archive/2026-09/PREREG_container_units_v2_dev_2026-09-30.md` Part C).
//!
//! The driver's `families` stage runs against a STUB `mcl_families` (a shell script that records its command line and
//! writes the products the driver summarises), so nothing here needs minimap2 or a genome: what is tested is what the
//! driver decides. `f1units` reads `PREFIX.families.gtf` like `f1` / `f1v2`, so the guard that refuses a missing or
//! older `families.gtf` is the one the other arms have; unset, the command is the one the driver has always run.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const DRIVER: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tools/rustle_pipeline.sh");

fn scratch(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR"))
        .join("pipeline_driver_units")
        .join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(d.join("bin")).expect("create scratch dir");
    let stub = d.join("bin/mcl_families");
    std::fs::write(
        &stub,
        "#!/bin/bash\n\
         if [ \"$1\" = --help ]; then echo '--emit-units <out>.copies.tsv --min-cov-shorter --emit-container --emit-relations'; exit 0; fi\n\
         ARGS=\"$*\"; echo \"$ARGS\" > \"$(dirname \"$0\")/args.txt\"\n\
         while [ $# -gt 0 ]; do case \"$1\" in --out) OUT=$2; shift 2;; *) shift;; esac; done\n\
         printf 'cluster_id\\tsize\\nMCL0\\t2\\n' > $OUT.clusters.tsv\n\
         printf 'family_id\\ttid\\nMCL0\\tt\\n' > $OUT.copies.tsv\n\
         case \" $ARGS \" in *\" --emit-relations \"*)\n\
           printf 'relation_id\\tx\\nREL1\\tF\\n' > $OUT.relations.tsv\n\
           printf 'family\\tlocus\\nMCL0\\tl\\n' > $OUT.members_by_locus.tsv;; esac\n",
    )
    .unwrap();
    Command::new("chmod")
        .args(["+x", stub.to_str().unwrap()])
        .status()
        .unwrap();
    d
}

/// `families` on `<dir>/run.*` with the stub, under `env` (the driver's own variables are cleared first).
fn families(dir: &Path, env: &[(&str, &str)]) -> Output {
    let mut c = Command::new("bash");
    c.arg(DRIVER)
        .args(["families", "--bam", "x.bam", "--fasta", "g.fa", "--out"])
        .arg(dir.join("run"))
        .arg("--bin")
        .arg(dir.join("bin"))
        .env_remove("RUSTLE_BRIDGE_REGROUP")
        .env_remove("RUSTLE_BRIDGE_UNITS_LIST")
        .env_remove("RUSTLE_FAMILY_RELATIONS")
        .env_remove("RUSTLE_GTF_REGROUP")
        .env("TMPDIR", dir);
    for (k, v) in env {
        c.env(k, v);
    }
    c.output().expect("bash failed to spawn")
}

fn err(o: &Output) -> String {
    String::from_utf8_lossy(&o.stderr).to_string()
}

fn recorded(dir: &Path) -> String {
    std::fs::read_to_string(dir.join("bin/args.txt"))
        .unwrap_or_default()
        .trim()
        .to_string()
}

/// Write `run.gtf` and then `run.families.gtf` (newer), or the other way round.
fn gtfs(dir: &Path, families_newer: bool) {
    let (first, second) = if families_newer {
        ("run.gtf", "run.families.gtf")
    } else {
        ("run.families.gtf", "run.gtf")
    };
    std::fs::write(
        dir.join(first),
        "c1\trustle\ttranscript\t1\t2\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\n",
    )
    .unwrap();
    std::thread::sleep(std::time::Duration::from_millis(50));
    std::fs::write(
        dir.join(second),
        "c1\trustle\ttranscript\t1\t2\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\n",
    )
    .unwrap();
}

#[test]
fn f1units_needs_a_families_gtf_no_older_than_the_gtf() {
    let dir = scratch("guard");
    std::fs::write(dir.join("run.gtf"), "x\n").unwrap();
    let o = families(&dir, &[("RUSTLE_BRIDGE_REGROUP", "f1units")]);
    assert_eq!(o.status.code(), Some(2), "{}", err(&o));
    assert!(
        err(&o).contains("RUSTLE_BRIDGE_REGROUP=f1units (unset = f1v2) but")
            && err(&o).contains("families.gtf is missing or older"),
        "{}",
        err(&o)
    );
    assert_eq!(recorded(&dir), "", "the stage did not start");
    // a families.gtf older than the GTF is stale too
    gtfs(&dir, false);
    let o = families(&dir, &[("RUSTLE_BRIDGE_REGROUP", "f1units")]);
    assert_eq!(o.status.code(), Some(2), "{}", err(&o));
    // newer: it reads PREFIX.families.gtf and the stage runs
    gtfs(&dir, true);
    let o = families(&dir, &[("RUSTLE_BRIDGE_REGROUP", "f1units")]);
    assert!(o.status.success(), "{}", err(&o));
    assert!(
        recorded(&dir).contains(&format!("--from-gtf {}/run.families.gtf", dir.display())),
        "{}",
        recorded(&dir)
    );
    // off with a families.gtf newer than the GTF means the GTF was assembled with a bridge mode: refused, as for f1 / f1v2
    let o = families(&dir, &[("RUSTLE_BRIDGE_REGROUP", "off")]);
    assert_eq!(o.status.code(), Some(2), "{}", err(&o));
    assert!(
        err(&o).contains("was assembled with a bridge mode"),
        "{}",
        err(&o)
    );
}

#[test]
fn unset_the_families_command_is_the_one_the_driver_always_ran_and_relations_add_one_flag() {
    let dir = scratch("command");
    gtfs(&dir, true);
    let base = format!(
        "--from-gtf {0}/run.families.gtf --fasta g.fa --threads 4 --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units --out {0}/run.fam",
        dir.display()
    );
    let o = families(&dir, &[]);
    assert!(o.status.success(), "{}", err(&o));
    assert_eq!(recorded(&dir), base, "the default (f1v2) families command");
    assert!(
        !dir.join("run.fam.relations.tsv").exists(),
        "no relations unless asked"
    );
    let o = families(&dir, &[("RUSTLE_FAMILY_RELATIONS", "0")]);
    assert!(o.status.success() && recorded(&dir) == base, "0 is off");
    // RUSTLE_FAMILY_RELATIONS=1: exactly one more flag, and the driver reports the tables
    let o = families(
        &dir,
        &[
            ("RUSTLE_FAMILY_RELATIONS", "1"),
            ("RUSTLE_BRIDGE_REGROUP", "f1units"),
        ],
    );
    assert!(o.status.success(), "{}", err(&o));
    assert_eq!(
        recorded(&dir),
        base.replace("--emit-units", "--emit-units --emit-relations")
    );
    assert!(
        err(&o).contains("families: relations 1 split transcripts, 0 cover, 1 members by locus"),
        "{}",
        err(&o)
    );
    let o = families(&dir, &[("RUSTLE_FAMILY_RELATIONS", "yes")]);
    assert_eq!(o.status.code(), Some(2));
    assert!(
        err(&o).contains("RUSTLE_FAMILY_RELATIONS must be 0 or 1"),
        "{}",
        err(&o)
    );
}

#[test]
fn the_units_list_belongs_to_f1units_and_to_a_file() {
    let dir = scratch("list");
    gtfs(&dir, true);
    let list = dir.join("bridges.tsv");
    std::fs::write(&list, "tid\tjunctions\nt\t1-2\n").unwrap();
    let l = list.to_str().unwrap();
    for arm in ["f1", "f1v2", "off"] {
        let o = families(
            &dir,
            &[
                ("RUSTLE_BRIDGE_REGROUP", arm),
                ("RUSTLE_BRIDGE_UNITS_LIST", l),
            ],
        );
        assert_eq!(o.status.code(), Some(2), "{arm}: {}", err(&o));
        assert!(
            err(&o).contains(
                "RUSTLE_BRIDGE_UNITS_LIST names the units of RUSTLE_BRIDGE_REGROUP=f1units"
            ),
            "{arm}: {}",
            err(&o)
        );
    }
    // unset means f1v2: the list is refused there too
    let o = families(&dir, &[("RUSTLE_BRIDGE_UNITS_LIST", l)]);
    assert_eq!(o.status.code(), Some(2), "{}", err(&o));
    let o = families(
        &dir,
        &[
            ("RUSTLE_BRIDGE_REGROUP", "f1units"),
            ("RUSTLE_BRIDGE_UNITS_LIST", "/nonexistent/l.tsv"),
        ],
    );
    assert_eq!(o.status.code(), Some(2));
    assert!(err(&o).contains("is missing or empty"), "{}", err(&o));
    let o = families(
        &dir,
        &[
            ("RUSTLE_BRIDGE_REGROUP", "f1units"),
            ("RUSTLE_BRIDGE_UNITS_LIST", "a b.tsv"),
        ],
    );
    assert_eq!(o.status.code(), Some(2));
    assert!(
        err(&o).contains("must not contain whitespace"),
        "{}",
        err(&o)
    );
    let o = families(
        &dir,
        &[
            ("RUSTLE_BRIDGE_REGROUP", "f1units"),
            ("RUSTLE_BRIDGE_UNITS_LIST", l),
        ],
    );
    assert!(o.status.success(), "{}", err(&o));
    let o = families(&dir, &[("RUSTLE_BRIDGE_REGROUP", "f2")]);
    assert_eq!(o.status.code(), Some(2));
    assert!(
        err(&o).contains("must be off, f1, f1v2 or f1units"),
        "{}",
        err(&o)
    );
    // f1units is refused with RUSTLE_GTF_REGROUP, as f1 / f1v2 are
    let o = families(
        &dir,
        &[
            ("RUSTLE_BRIDGE_REGROUP", "f1units"),
            ("RUSTLE_GTF_REGROUP", "1"),
        ],
    );
    assert_eq!(o.status.code(), Some(2));
    assert!(
        err(&o).contains("already splits every gene_id"),
        "{}",
        err(&o)
    );
}
