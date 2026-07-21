import subprocess, sys, tempfile, os, pathlib
GEN = pathlib.Path(__file__).resolve().parents[2] / "bench/mechanism/gen_heuristics.py"
REPO = pathlib.Path(__file__).resolve().parents[2]

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
