//! `copy_assign --assemble-only --bridge-regroup` end to end on a synthetic readthrough locus (`vg_family::bridge_regroup`).
//!
//! Gene X's reads end at a canonical PAS inside the intron [1251, 2200] that a 3-read bridge B skips on its way to
//! gene Y, whose reads start inside that intron at their own promoter. The assembler emits X, B and Y under one
//! `gene_id`; `--bridge-regroup` splits X and Y and keeps B as a relation. The frozen `bench/f1_bridge.py` (37ee8e77)
//! and `f1v2.py` (b4e788ad), run on the `--bridge-regroup off` GTF of this fixture, write exactly the files asserted
//! here (checked 2026-09-29). `f1v2` is the default since 2026-09-29: `off` is passed explicitly wherever the fused
//! locus is wanted, and the default run is asserted equal to the explicit `f1v2` run.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use noodles_core::Position;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::cigar::{op::Kind, Op};
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record::{Flags, MappingQuality};
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::{Record as _, RecordBuf};
use noodles_sam::header::record::value::{map::ReferenceSequence, Map};

const L: usize = 3000;
/// The flag-off assembly's one gene_id and its three transcripts.
const G: &str = "DN_c1_1000_2";
const X: &str = "DN_c1_1000_2";
const B: &str = "DN_c1_1000_4";
const Y: &str = "DN_c1_1997_3";
/// The ghost fixture's fourth transcript Z (20 reads), whose id is then the locus's gene_id.
const Z: &str = "DN_c1_2600_3";

fn scratch(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("copy_assign_bridge_regroup").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).expect("create scratch dir");
    d
}

/// The reads as 1-based closed exons: X (3' ends 1296-1300), the bridge B, Y (starts 1998-2001). Ends vary so that the
/// assembler's coordinate dedupe keeps >= 2 reads per chain. With `ghost`, also a 2-read read-through C from Y's
/// last two exons over a 250 bp exon into Z's last exon, and Z (20 reads) whose junction 2701-2800 lies inside that
/// exon: C joins Z to the locus, then `--assembly-polish`'s retained-intron rule (support 20 >= 10 x 2) drops it,
/// which leaves the gene_id in two exon-disjoint pieces (X/B/Y and Z), the `--gtf-regroup` case.
fn reads(ghost: bool) -> Vec<Vec<(u64, u64)>> {
    let mut v: Vec<Vec<(u64, u64)>> = Vec::new();
    for end in [1296, 1297, 1298, 1299, 1300, 1300] {
        v.push(vec![(1001, 1100), (1201, end)]);
    }
    for end in [2498, 2499, 2500] {
        v.push(vec![(1001, 1100), (1201, 1250), (2201, 2300), (2401, end)]);
    }
    for start in [1998, 1999, 2000, 2001, 2001] {
        v.push(vec![(start, 2100), (2201, 2300), (2401, 2500)]);
    }
    if ghost {
        for end in [2950, 2949] {
            v.push(vec![(2201, 2300), (2401, 2500), (2601, 2850), (2901, end)]);
        }
        for k in 0..20 {
            v.push(vec![(2601, 2700), (2801, 2850), (2901, 2950 - k)]);
        }
    }
    v
}

/// `c1`: all C but the canonical splice sites (GT..AG) and, when `pas`, AATAAA at 1270-1275 (inside [1300-35,
/// 1300-10] of X's 3' mode 1300; the 20 bp downstream are C: not internally primed). `ghost` adds the splice sites
/// of C and Z.
fn write_fixture(dir: &Path, pas: bool, ghost: bool) -> (String, String) {
    let mut seq = vec![b'C'; L];
    let mut put = |pos1: usize, s: &[u8]| seq[pos1 - 1..pos1 - 1 + s.len()].copy_from_slice(s);
    for d in [1101, 1251, 2101, 2301] {
        put(d, b"GT");
    }
    for a in [1199, 2199, 2399] {
        put(a, b"AG");
    }
    if ghost {
        for d in [2501, 2701, 2851] {
            put(d, b"GT");
        }
        for a in [2599, 2799, 2899] {
            put(a, b"AG");
        }
    }
    if pas {
        put(1270, b"AATAAA");
    }
    let fasta = dir.join("genome.fa");
    std::fs::write(&fasta, format!(">c1\n{}\n", String::from_utf8(seq).unwrap())).unwrap();
    std::fs::write(dir.join("genome.fa.fai"), format!("c1\t{L}\t4\t{L}\t{}\n", L + 1)).unwrap();

    let header = noodles_sam::Header::builder()
        .add_reference_sequence("c1", Map::<ReferenceSequence>::new(std::num::NonZeroUsize::try_from(L).unwrap()))
        .build();
    let mut recs: Vec<RecordBuf> = reads(ghost)
        .iter()
        .enumerate()
        .map(|(i, ex)| {
            let mut ops = vec![Op::new(Kind::Match, (ex[0].1 - ex[0].0 + 1) as usize)];
            for w in ex.windows(2) {
                ops.push(Op::new(Kind::Skip, (w[1].0 - w[0].1 - 1) as usize));
                ops.push(Op::new(Kind::Match, (w[1].1 - w[1].0 + 1) as usize));
            }
            let data: noodles_sam::alignment::record_buf::Data =
                [(Tag::new(b't', b's'), Value::Character(b'+'))].into_iter().collect();
            RecordBuf::builder()
                .set_name(format!("r{i}"))
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(Position::try_from(ex[0].0 as usize).unwrap())
                .set_mapping_quality(MappingQuality::new(60).unwrap())
                .set_cigar(ops.into_iter().collect())
                .set_data(data)
                .build()
        })
        .collect();
    recs.sort_by_key(|r| r.alignment_start().map(|p| p.get()).unwrap_or(0));
    let bam = dir.join("reads.bam");
    {
        let mut w = noodles_bam::io::Writer::new(std::fs::File::create(&bam).unwrap());
        w.write_header(&header).unwrap();
        for r in &recs {
            w.write_alignment_record(&header, r).unwrap();
        }
        w.try_finish().unwrap();
    }
    let mut reader = noodles_bam::io::Reader::new(std::fs::File::open(&bam).unwrap());
    let header = reader.read_header().unwrap();
    let mut indexer = noodles_csi::binning_index::Indexer::default();
    let mut chunk_start = reader.get_ref().virtual_position();
    let mut record = noodles_bam::Record::default();
    while reader.read_record(&mut record).unwrap() != 0 {
        let chunk_end = reader.get_ref().virtual_position();
        let ctx = match (
            record.reference_sequence_id().transpose().unwrap(),
            record.alignment_start().transpose().unwrap(),
            record.alignment_end().transpose().unwrap(),
        ) {
            (Some(id), Some(start), Some(end)) => Some((id, start, end, !record.flags().is_unmapped())),
            _ => None,
        };
        indexer
            .add_record(ctx, noodles_csi::binning_index::index::reference_sequence::bin::Chunk::new(chunk_start, chunk_end))
            .unwrap();
        chunk_start = chunk_end;
    }
    let index: noodles_bam::bai::Index = indexer.build(header.reference_sequences().len());
    noodles_bam::bai::write(dir.join("reads.bam.bai"), &index).unwrap();
    (bam.display().to_string(), fasta.display().to_string())
}

/// `copy_assign --assemble-only` on the fixture's one region; never asserts success.
fn run(dir: &Path, fx: &(String, String), out: &str, extra: &[&str]) -> (Output, String) {
    let out = dir.join(out).display().to_string();
    let o = Command::new(env!("CARGO_BIN_EXE_copy_assign"))
        .args(["--bam", &fx.0, "--fasta", &fx.1, "--region", "c1:1-3000", "--assemble-only", "--out", &out])
        .args(extra)
        .output()
        .expect("copy_assign failed to spawn");
    (o, out)
}

fn ok(o: &Output) {
    assert!(o.status.success(), "copy_assign failed:\n{}", String::from_utf8_lossy(&o.stderr));
}

fn read(path: &str) -> String {
    std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {path}: {e}"))
}

fn attr<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let pat = format!("{key} \"");
    let i = line.find(&pat)? + pat.len();
    Some(&line[i..i + line[i..].find('"')?])
}

/// transcript_id -> gene_id over the `transcript` lines.
fn genes(gtf: &str) -> Vec<(String, String)> {
    gtf.lines()
        .filter(|l| l.split('\t').nth(2) == Some("transcript"))
        .map(|l| (attr(l, "transcript_id").unwrap().to_string(), attr(l, "gene_id").unwrap().to_string()))
        .collect()
}

fn pairs(v: &[(&str, &str)]) -> Vec<(String, String)> {
    v.iter().map(|(a, b)| (a.to_string(), b.to_string())).collect()
}

/// A line without its gene_id and relation attributes (what `--bridge-regroup` may not change).
fn strip(line: &str) -> String {
    let mut s = line.replacen(&format!("gene_id \"{}\"", attr(line, "gene_id").unwrap_or("")), "", 1);
    if let Some(i) = s.find(" fusion_of \"") {
        s.truncate(i);
    }
    s
}

#[test]
fn a_proven_bridge_is_split_off_and_kept_as_a_relation() {
    let dir = scratch("bridge");
    let fx = write_fixture(&dir, true, false);
    let (o, off) = run(&dir, &fx, "off", &["--bridge-regroup", "off"]);
    ok(&o);
    assert_eq!(genes(&read(&format!("{off}.gtf"))), pairs(&[(X, G), (B, G), (Y, G)]), "the fused locus");
    for ext in ["families.gtf", "bridge_junctions.tsv", "bridges.tsv"] {
        assert!(!Path::new(&format!("{off}.{ext}")).exists(), "off writes no {ext}");
    }

    let (o, f1) = run(&dir, &fx, "f1", &["--bridge-regroup", "f1"]);
    ok(&o);
    let gtf = read(&format!("{f1}.gtf"));
    let (g2, fus) = (format!("{G}.rg2"), format!("{G}.fus1"));
    assert_eq!(genes(&gtf), pairs(&[(X, G), (B, &fus), (Y, &g2)]), "X keeps the name (5 reads > 4)");
    let b_line = gtf.lines().find(|l| l.contains("\ttranscript\t") && attr(l, "transcript_id") == Some(B)).unwrap();
    assert!(
        b_line.ends_with(&format!("low_confidence_reason \"none\"; fusion_of \"{G},{g2}\"; fusion_junction \"c1:1251-2200:+\";")),
        "{b_line}"
    );
    let off_gtf = read(&format!("{off}.gtf"));
    assert_eq!(off_gtf.lines().map(strip).collect::<Vec<_>>(), gtf.lines().map(strip).collect::<Vec<_>>(), "only names move");
    let fam = read(&format!("{f1}.families.gtf"));
    assert_eq!(
        fam.lines().collect::<Vec<_>>(),
        gtf.lines().filter(|l| attr(l, "transcript_id") != Some(B)).collect::<Vec<_>>(),
        "the families input is the GTF without the bridge's 5 lines"
    );
    let junctions = "gene\tchrom\ts\te\tstrand\tn_TJ\treads_TJ\tup\tdown\tinside\tU\tclusters\tup_proof\tV1\tdown_proof\tbridge\tTJ\n\
                     DN_c1_1000_2\tc1\t1251\t2200\t+\t1\t3\t1\t1\t0\t5\t1300:5:P\tTrue\t5\tTrue\tTrue\tDN_c1_1000_4\n";
    assert_eq!(read(&format!("{f1}.bridge_junctions.tsv")), junctions);
    assert!(!Path::new(&format!("{f1}.bridges.tsv")).exists(), "bridges.tsv is F1v2's table");
    let (p_off, p_f1) = (read(&format!("{off}.params.tsv")), read(&format!("{f1}.params.tsv")));
    let added: Vec<&str> = p_f1.lines().filter(|l| !p_off.lines().any(|x| x == *l)).collect();
    assert_eq!(added.len(), 15, "{added:?}");
    assert!(added.contains(&"bridge_regroup\tf1") && added.contains(&"bridge_regroup_bridge_transcripts\t1"));

    // F1v2: 3 reads < 5 (X) and < 4 (Y): a minority link, so the same regrouping, plus its table
    let (o, v2) = run(&dir, &fx, "f1v2", &["--bridge-regroup", "f1v2"]);
    ok(&o);
    assert_eq!(read(&format!("{v2}.gtf")), gtf);
    assert_eq!(read(&format!("{v2}.families.gtf")), fam);
    assert_eq!(read(&format!("{v2}.bridge_junctions.tsv")), junctions);
    assert_eq!(
        read(&format!("{v2}.bridges.tsv")),
        "gene\tchrom\ts\te\tstrand\tn_TJ\treads_TJ\tmax_TJ\tn_up\treads_up\tmax_up\tn_down\treads_down\tmax_down\tshare\tkeep\tTJ\n\
         DN_c1_1000_2\tc1\t1251\t2200\t+\t1\t3\t3\t1\t5\t5\t1\t4\t4\t0.4286\tTrue\tDN_c1_1000_4\n"
    );

    // the default (no flag) is f1v2 since 2026-09-29: every product of the explicit arm, byte for byte
    let (o, def) = run(&dir, &fx, "default", &[]);
    ok(&o);
    for ext in ["gtf", "families.gtf", "bridge_junctions.tsv", "bridges.tsv", "params.tsv"] {
        assert_eq!(read(&format!("{def}.{ext}")), read(&format!("{v2}.{ext}")), "{ext}: the default is f1v2");
    }
    assert!(String::from_utf8_lossy(&o.stderr).contains("BRIDGE REGROUP (f1v2)"), "the default logs its arm");

    // the buffered reader (its own indexed evidence pass) decides the same
    let (o, buf) = run(&dir, &fx, "buf", &["--bridge-regroup", "f1", "--materialize-reads"]);
    ok(&o);
    for ext in ["gtf", "families.gtf", "bridge_junctions.tsv"] {
        assert_eq!(read(&format!("{buf}.{ext}")), read(&format!("{f1}.{ext}")), "{ext}: buffered vs streaming");
    }
}

/// Without the PAS, X's 3' cluster is not proven: no bridge, and the output is `--gtf-regroup`'s. The ghost fixture
/// makes that a real rename: the polish (`--assembly-polish mono`, whose retained-intron rule is on by default) drops
/// C after it joined Z to the locus, so the one gene_id (Z's, 20 reads) is two exon-disjoint pieces, and both flags
/// name X, B and Y `.rg2`.
#[test]
fn without_a_pas_there_is_no_bridge_and_the_names_are_rg3s() {
    let dir = scratch("no_pas");
    let fx = write_fixture(&dir, false, true);
    let (o, off) = run(&dir, &fx, "off", &["--assembly-polish", "mono", "--bridge-regroup", "off"]);
    ok(&o);
    let (o, rg3) = run(&dir, &fx, "rg3", &["--assembly-polish", "mono", "--gtf-regroup", "--bridge-regroup", "off"]);
    ok(&o);
    let (o, f1) = run(&dir, &fx, "f1", &["--assembly-polish", "mono", "--bridge-regroup", "f1"]);
    ok(&o);
    let off_gtf = read(&format!("{off}.gtf"));
    assert_eq!(genes(&off_gtf), pairs(&[(X, Z), (B, Z), (Y, Z), (Z, Z)]), "one gene_id; the ghost C was dropped");
    let gtf = read(&format!("{f1}.gtf"));
    let rg2 = format!("{Z}.rg2");
    assert_eq!(genes(&gtf), pairs(&[(X, &rg2), (B, &rg2), (Y, &rg2), (Z, Z)]), "Z (20 reads) keeps the name");
    assert_eq!(gtf, read(&format!("{rg3}.gtf")), "no bridge: the names are --gtf-regroup's");
    assert_ne!(gtf, off_gtf);
    assert_eq!(read(&format!("{f1}.families.gtf")), gtf);
    let table = read(&format!("{f1}.bridge_junctions.tsv"));
    assert!(table.lines().nth(1).unwrap().ends_with("\t1300:5:u\tFalse\t5\tTrue\tFalse\tDN_c1_1000_4"), "{table}");
}

#[test]
fn bridge_regroup_is_refused_where_it_cannot_run() {
    let dir = scratch("refused");
    let fx = write_fixture(&dir, true, false);
    let err = |o: &Output| String::from_utf8_lossy(&o.stderr).to_string();
    let (o, _) = run(&dir, &fx, "both", &["--bridge-regroup", "f1", "--gtf-regroup"]);
    assert!(!o.status.success() && err(&o).contains("pass one of them"), "{}", err(&o));
    // the default arm (f1v2) contains --gtf-regroup's split too: RG3 alone needs an explicit off
    let (o, _) = run(&dir, &fx, "rg3_alone", &["--gtf-regroup"]);
    assert!(
        !o.status.success() && err(&o).contains("f1v2 (the default)") && err(&o).contains("--bridge-regroup off"),
        "{}",
        err(&o)
    );
    let regions = dir.join("regions.txt");
    std::fs::write(&regions, "c1:1-1500\nc1:1501-3000\n").unwrap();
    let o = Command::new(env!("CARGO_BIN_EXE_copy_assign"))
        .args(["--bam", &fx.0, "--fasta", &fx.1, "--regions", regions.to_str().unwrap(), "--assemble-only"])
        .args(["--bridge-regroup", "f1", "--out", dir.join("two").to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!o.status.success() && err(&o).contains("one region per contig"), "{}", err(&o));
    let o = Command::new(env!("CARGO_BIN_EXE_copy_assign"))
        .args(["--bam", &fx.0, "--fasta", &fx.1, "--region", "c1:1-3000", "--gtf", "--bridge-regroup", "f1"])
        .args(["--out", dir.join("gtf").to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!o.status.success() && err(&o).contains("--assemble-only"), "{}", err(&o));
}
