import subprocess, sys, tempfile, os, pathlib
GEN = pathlib.Path(__file__).resolve().parents[2] / "bench/mechanism/gen_heuristics.py"
REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(GEN.parent))
from gen_heuristics import verify_entry

def _toml(body):
    f = tempfile.NamedTemporaryFile("w", suffix=".toml", delete=False)
    f.write(body); f.close(); return f.name

def test_verify_passes_on_real_line():
    # denovo_pipeline.rs:2370 really is `min_identity: 0.80,`
    t = _toml('[[heuristic]]\nname="edge_id"\nfile="src/rustle/vg_family/denovo_pipeline.rs"\n'
              'line=2370\nvalue="0.80"\nstage="O1-edge"\ntier="decision"\nkind="arbitrary"\nrationale="edge identity floor"\n')
    out = tempfile.mktemp(suffix=".tsv")
    r = subprocess.run([sys.executable, str(GEN), "--repo-root", str(REPO), "--toml", t, "--out", out],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    rows = open(out).read().splitlines()
    assert rows[0] == "stage\ttier\tkind\tname\tvalue\tfile\tline\trationale"
    assert any("edge_id" in row and "0.80" in row for row in rows[1:])

def test_verify_fails_on_drifted_value():
    # value that is NOT on that line -> non-zero exit, names the entry
    t = _toml('[[heuristic]]\nname="edge_id"\nfile="src/rustle/vg_family/denovo_pipeline.rs"\n'
              'line=2370\nvalue="0.99"\nstage="O1-edge"\ntier="decision"\nkind="arbitrary"\nrationale="x"\n')
    out = tempfile.mktemp(suffix=".tsv")
    r = subprocess.run([sys.executable, str(GEN), "--repo-root", str(REPO), "--toml", t, "--out", out],
                       capture_output=True, text=True)
    assert r.returncode != 0
    assert "edge_id" in (r.stderr + r.stdout)

def test_verify_fails_on_digit_superset_drift():
    # unit test against verify_entry directly, on a synthetic source file, to prove
    # the digit-boundary fix: "3" must NOT match inside "13" (a digit-superset drift),
    # while "13" must match exactly.
    src_dir = tempfile.mkdtemp()
    src_path = pathlib.Path(src_dir) / "fake.rs"
    src_path.write_text("pub const GATE_MIN_READS: u32 = 13;\n")

    entry_drifted = {"name": "gate_min_reads", "file": "fake.rs", "line": 1, "value": "3"}
    err = verify_entry(entry_drifted, src_dir)
    assert err is not None, "digit-superset drift (3 -> 13) must be caught, not silently pass"
    assert "gate_min_reads" in err

    entry_correct = {"name": "gate_min_reads", "file": "fake.rs", "line": 1, "value": "13"}
    assert verify_entry(entry_correct, src_dir) is None
