"""Tests of exact_truth.py (bench/famwise_eval): tiny unsorted BAMs, hand-made origin tables and GTFs."""
import os, subprocess, sys
import pytest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import records as R
import origin_bam as OB
from test_records import bam
import exact_truth as E

G = "0-100,150-250"                    # the source gene: the copy at c1:0, strand +; records are 100=50N100=
RECS = [("home", "c1", 0, 0, 200, 200), ("home", "c1", 1000, 256, 180, 180),      # S1, primary at the origin
        ("off", "c1", 1000, 0, 200, 200), ("off", "c1", 0, 256, 195, 195),        # lead 5, the secondary is at the origin
        ("tie0", "c1", 1000, 0, 200, 200), ("tie0", "c1", 0, 256, 200, 200),      # lead 0
        ("ten", "c1", 0, 0, 200, 200), ("ten", "c1", 1000, 256, 190, 190),        # lead 10 is still near-tied
        ("out", "c1", 1000, 0, 200, 200),                                         # no record at the origin
        ("merge", "c1", 0, 0, 200, 200), ("merge", "c1", 50, 256, 190, 190),      # two records, one place
        ("strand", "c1", 0, 16, 200, 200),                                        # the source's contig, the other strand
        ("intron", "c1", 0, 0, 200, 200),                                         # inside the gene span, on no exon
        ("contig", "c2", 0, 0, 200, 200),
        ("two", "c1", 1000, 0, 180, 180), ("two", "c1", 0, 256, 200, 200)]        # two places on the source, primary at the worse
SRC = {q: ("c1", "+", G) for q in ("home", "off", "tie0", "ten", "out", "merge", "strand", "contig")}
SRC.update(intron=("c1", "+", "100-150"), two=("c1", "+", "0-100,1000-1100"))

def origins(tmp_path, src=SRC, name="origin.tsv"):
    (tmp_path / name).write_text("read\tchrom\tstrand\tblocks\tgene\n" + "".join(f"{q}\t{c}\t{s}\t{b}\tG1\n" for q, (c, s, b) in src.items()))
    return str(tmp_path / name)

def run(tmp_path, capsys, *args, recs=RECS):
    """-> ({read: row dict}, stdout tables as lists of fields, path of the table)"""
    out = str(tmp_path / "exact.tsv")
    E.main([bam(tmp_path, recs), out, *args])
    lines = open(out).read().splitlines()
    assert lines[0] == "read\tstratum\tlead\tstatus\tpicks\torigin_loci\tloci\tsource\thome"
    rows = {l.split("\t")[0]: dict(zip(lines[0].split("\t"), l.split("\t"))) for l in lines[1:]}
    return rows, [l.split("\t") for l in capsys.readouterr().out.splitlines() if l], out

def test_origin_rule(tmp_path, capsys):
    rows, _, _ = run(tmp_path, capsys, "--origins", origins(tmp_path))
    got = {q: (r["stratum"], r["lead"], r["status"], r["picks"], r["origin_loci"], r["home"]) for q, r in rows.items()}
    assert got == {"home": ("S1", "20", "decidable", "primary=0", "0", "1"), "off": ("S2", "5", "decidable", "primary=0", "1", "1"),
                   "tie0": ("S2", "0", "decidable", "primary=0", "1", "1"), "ten": ("S2", "10", "decidable", "primary=0", "0", "1"),
                   "out": ("S0", "None", "outside", "primary=0", "", "0"), "merge": ("S0", "None", "decidable", "primary=0", "0", "1"),
                   "strand": ("S0", "None", "outside", "primary=0", "", "0"), "intron": ("S0", "None", "outside", "primary=0", "", "0"),
                   "contig": ("S0", "None", "outside", "primary=0", "", "0"), "two": ("S1", "20", "decidable", "primary=1", "0", "2")}
    assert rows["off"]["loci"] == "0:c1:1000-1250:+:AS200;1:c1:0-250:+:AS195" and rows["off"]["source"] == "c1:+"
    assert rows["merge"]["loci"] == "0:c1:0-300:+:AS200" and rows["strand"]["loci"] == "0:c1:0-250:-:AS200"

def test_tables(tmp_path, capsys):
    _, out, _ = run(tmp_path, capsys, "--origins", origins(tmp_path))
    w = lambda k, n: f"{k / n:.4f}\t" + "-".join(f"{x:.4f}" for x in R.wilson(k, n))
    assert out == [l.split("\t") for l in (
        "stratum\treads\toutside\tdecidable\tprimary right\tprimary wrong", "S0\t5\t4\t1\t1\t0", "S1\t2\t0\t2\t1\t1", "S2\t3\t0\t3\t1\t2",
        "lead\treads\tprimary at origin\tshare\twilson95", f"clear\t7\t2\t{w(2, 7)}", f"1-10\t2\t1\t{w(1, 2)}", f"0\t1\t0\t{w(0, 1)}")]

def test_command_line_from_another_directory(tmp_path, capsys):
    _, out, _ = run(tmp_path, capsys, "--origins", origins(tmp_path))
    p = subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(E.__file__)), "exact_truth.py"), "sel.bam", "cli.tsv",
                        "--origins", "origin.tsv"], cwd=tmp_path, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    assert [l.split("\t") for l in p.stdout.splitlines() if l] == out and (tmp_path / "cli.tsv").read_text() == (tmp_path / "exact.tsv").read_text()

def test_cf_split(tmp_path, capsys):
    cf = tmp_path / "cf.tsv"
    cf.write_text("read\tstratum\tlead\tstatus\tpicks\torigin_loci\tloci\nhome\tS1\t20\tdecidable\tprimary=0\t0\t\n"
                  "off\tS2\t5\ttie\tprimary=0\t0,1\t\ntie0\tS2\t0\toutside\tprimary=0\t\t\n")
    _, out, _ = run(tmp_path, capsys, "--origins", origins(tmp_path), "--cf", str(cf))
    i = out.index(["lead", "catalogue-free", "reads", "primary at origin", "share", "wilson95"])
    assert out[i - 3][:3] == ["clear", "7", "2"]                    # the table without the split stays
    assert [r[:4] for r in out[i + 1:i + 6]] == [["clear", "decidable", "1", "1"], ["clear", "not in the table", "6", "1"],
                                            ["1-10", "tie", "1", "0"], ["1-10", "not in the table", "1", "1"], ["0", "outside", "1", "0"]]

def test_truth_check_reads_the_catalogue_free_calls_against_the_exact_origin(tmp_path, capsys):
    cf = tmp_path / "cf.tsv"                # home: right call; off: exact origin is locus 1; tie0: its tie holds locus 1; ten: a tie without the origin; out: both outside
    cf.write_text("read\tstratum\tlead\tstatus\tpicks\torigin_loci\tloci\nhome\tS1\t20\tdecidable\tprimary=0\t0\t\n"
                  "off\tS2\t5\tdecidable\tprimary=0\t0\t\ntie0\tS2\t0\ttie\tprimary=0\t0,1\t\nten\tS2\t10\ttie\tprimary=0\t1\t\n"
                  "out\tS0\tNone\toutside\tprimary=0\t\t\ncontig\tS0\tNone\tdecidable\tprimary=0\t0\t\n")
    _, out, _ = run(tmp_path, capsys, "--origins", origins(tmp_path), "--cf", str(cf))
    i = out.index(["catalogue-free", "reads", "exact origin is the call / in the tie", "share", "wilson95", "exact outside"])
    assert [r[:3] + r[5:] for r in out[i + 1:i + 4]] == [["decidable", "3", "1", "1"], ["tie", "2", "1", "0"], ["outside", "1", "0", "1"]]

def test_null_puts_the_origin_at_the_primary(tmp_path, capsys):
    rows, out, _ = run(tmp_path, capsys, "--origins", origins(tmp_path), "--null")
    assert {r["status"] for r in rows.values()} == {"decidable"}
    assert all(r["picks"] == "primary=" + r["origin_loci"] for r in rows.values()) and rows["two"]["origin_loci"] == "1"
    assert out[1:4] == [["S0", "5", "0", "5", "5", "0"], ["S1", "2", "0", "2", "2", "0"], ["S2", "3", "0", "3", "3", "0"]]
    assert [r[:3] for r in out[5:]] == [["clear", "7", "7"], ["1-10", "2", "2"], ["0", "1", "1"]]

def test_readers_downstream(tmp_path, capsys):
    """records.load_truth and origin_bam.plan take the table as it is"""
    rows, _, out = run(tmp_path, capsys, "--origins", origins(tmp_path))
    t, recs = R.load_truth(out), R.load_records(str(tmp_path / "sel.bam"))
    assert {q: v.kind for q, v in t.items()} == {q: "one" if r["status"] == "decidable" else "outside" for q, r in rows.items()}
    assert t["off"].origin == {1} and t["off"].loci == [("c1", 1000, 1250, "+", 200), ("c1", 0, 250, "+", 195)]
    case = {q: OB.plan(recs[q], t[q])[0] for q in t}
    assert case["off"] == case["tie0"] == case["two"] == "moved" and OB.plan(recs["off"], t["off"])[2:] == ([1], 1)
    assert {case[q] for q in t if t[q].kind == "outside"} == {"primary_only"}
    assert "no_home" not in case.values() and case["home"] == "restricted" and case["merge"] == "kept_as_is"

def test_missing_read_fails(tmp_path, capsys):
    with pytest.raises(SystemExit, match="'out'"):
        run(tmp_path, capsys, "--origins", origins(tmp_path, {q: v for q, v in SRC.items() if q != "out"}))

def test_origins_needs_its_header(tmp_path, capsys):
    (tmp_path / "bare.tsv").write_text("home\tc1\t+\t0-100\n")
    with pytest.raises(SystemExit, match="header"):
        run(tmp_path, capsys, "--origins", str(tmp_path / "bare.tsv"))
    (tmp_path / "twice.tsv").write_text("read\tchrom\tstrand\tblocks\nhome\tc1\t+\t0-100\nhome\tc1\t+\t1000-1100\n")
    with pytest.raises(SystemExit, match="'home' has two rows"):
        run(tmp_path, capsys, "--origins", str(tmp_path / "twice.tsv"))
    assert not os.path.exists(tmp_path / "exact.tsv")

GTF = ('c1\tsim\ttranscript\t1\t250\t.\t+\t.\tgene_id "g"; transcript_id "T_1";\n'
       'c1\tsim\texon\t1\t100\t.\t+\t.\tgene_id "g"; transcript_id "T_1";\nc1\tsim\texon\t151\t250\t.\t+\t.\tgene_id "g"; transcript_id "T_1";\n'
       'c1\tsim\texon\t1101\t1200\t.\t+\t.\tgene_id "h"; transcript_id "T2";\n')

def test_gtf_source_from_the_read_name(tmp_path, capsys):
    (tmp_path / "t.gtf").write_text(GTF)
    one = lambda q: [(q, "c1", 0, 0, 200, 200)]
    a, b = "PacBio_simulated_read_7_T_1", "PacBio_simulated_read_8_T2"
    rows, _, _ = run(tmp_path, capsys, "--gtf", str(tmp_path / "t.gtf"), recs=one(a) + one(b))
    assert (rows[a]["status"], rows[a]["source"], rows[a]["origin_loci"]) == ("decidable", "T_1", "0")
    assert (rows[b]["status"], rows[b]["source"]) == ("outside", "T2")           # T2's exon is past the read
    rows, _, _ = run(tmp_path, capsys, "--gtf", str(tmp_path / "t.gtf"), "--name-re", r"^[^|]+\|(.+)$", recs=one("m1|T_1"))
    assert rows["m1|T_1"]["status"] == "decidable"
    with pytest.raises(SystemExit, match="T9"):
        run(tmp_path, capsys, "--gtf", str(tmp_path / "t.gtf"), recs=one("PacBio_simulated_read_9_T9"))
    with pytest.raises(SystemExit, match="'odd'"):
        run(tmp_path, capsys, "--gtf", str(tmp_path / "t.gtf"), recs=one("odd"))
