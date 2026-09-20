import subprocess, sys, tempfile, os, pathlib
GEN = pathlib.Path(__file__).resolve().parents[2] / "bench/mechanism/gen_heuristics.py"
REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(GEN.parent))
from gen_heuristics import verify_entry, relocate, anchor_of

def _toml(body):
    f = tempfile.NamedTemporaryFile("w", suffix=".toml", delete=False)
    f.write(body); f.close(); return f.name

def _fake_repo(text):
    d = tempfile.mkdtemp()
    (pathlib.Path(d) / "fake.rs").write_text(text)
    return d

def _entry(**kw):
    e = {"name": "edge_id", "file": "fake.rs", "line": 1, "value": "0.80",
         "stage": "O1-edge", "tier": "decision", "kind": "arbitrary", "rationale": "edge identity floor"}
    e.update(kw); return e

def _toml_of(e):
    return _toml("[[heuristic]]\n" + "".join(f'{k}={v!r}\n'.replace("'", '"') for k, v in e.items()))

# These two tests used to hardcode `denovo_pipeline.rs:2370` ("really is min_identity: 0.80"). That line
# drifts every time anything is inserted above it, so the positive test failed and -- worse -- the negative
# test started passing for the WRONG reason (it expected a non-zero exit from a drifted VALUE and got one
# from a drifted LINE). Both now run against a synthetic source, so they test the generator's logic rather
# than the current shape of a real file.

def test_verify_passes_on_real_line():
    d = _fake_repo("// header\nmin_identity: 0.80,\n")
    t = _toml_of(_entry(line=2))
    out = tempfile.mktemp(suffix=".tsv")
    r = subprocess.run([sys.executable, str(GEN), "--repo-root", d, "--toml", t, "--out", out],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    rows = open(out).read().splitlines()
    assert rows[0] == "stage\ttier\tkind\tname\tvalue\tfile\tline\trationale"
    assert any("edge_id" in row and "0.80" in row for row in rows[1:])

def test_verify_fails_on_drifted_value():
    # The recorded line EXISTS and holds a value -- just not this one. Distinct from a drifted line number.
    d = _fake_repo("// header\nmin_identity: 0.80,\n")
    t = _toml_of(_entry(line=2, value="0.99"))
    out = tempfile.mktemp(suffix=".tsv")
    r = subprocess.run([sys.executable, str(GEN), "--repo-root", d, "--toml", t, "--out", out],
                       capture_output=True, text=True)
    assert r.returncode != 0
    assert "edge_id" in (r.stderr + r.stdout)

def test_verify_fails_on_digit_superset_drift():
    # unit test against verify_entry directly, on a synthetic source file, to prove
    # the digit-boundary fix: "3" must NOT match inside "13" (a digit-superset drift),
    # while "13" must match exactly.
    src_dir = _fake_repo("pub const GATE_MIN_READS: u32 = 13;\n")
    entry_drifted = {"name": "gate_min_reads", "file": "fake.rs", "line": 1, "value": "3"}
    err = verify_entry(entry_drifted, src_dir)
    assert err is not None, "digit-superset drift (3 -> 13) must be caught, not silently pass"
    assert "gate_min_reads" in err

    entry_correct = {"name": "gate_min_reads", "file": "fake.rs", "line": 1, "value": "13"}
    assert verify_entry(entry_correct, src_dir) is None

def test_relocate_prefers_the_declaration_over_a_doc_comment():
    """A doc comment citing the default sits directly ABOVE the declaration, so "first candidate at or below
    the old line" lands on the comment. That bug relocated min_identity, READ_CAP, PAD and 5 others onto
    `///` lines -- each a false claim about where the constant lives."""
    d = _fake_repo("// pad\n// pad\n/// defaults to asm20 min_identity 0.80 for the edge test\n"
                   "min_identity: 0.80,\n")
    new_line, why = relocate(_entry(line=1, name="RefineParams::min_identity"), d)
    assert new_line == 4, f"must pick the declaration (4), not the doc comment (3); got {new_line} ({why})"

def test_relocate_refuses_when_the_value_actually_changed():
    """A moved constant is bookkeeping; a retuned one invalidates the published registry. They must not be
    conflated, so relocate must REFUSE rather than silently repoint at some other 0.80 in the file."""
    d = _fake_repo("// pad\nmin_identity: 0.95,\nsomething_else: 0.80,\n")
    new_line, why = relocate(_entry(line=2, name="RefineParams::min_identity"), d)
    assert new_line is None, f"must refuse, got line {new_line}"
    assert "REAL CHANGE" in why

def test_relocate_anchors_a_cli_flag_whose_value_is_on_the_attribute_line():
    """clap puts the default on `#[arg(default_value_t = ..)]` and the field name on the NEXT line, so the
    anchor has to be found in a window around the value, not on its own line."""
    d = _fake_repo("// pad\n#[arg(long, default_value_t = 6.9)]\nmargin: f64,\n")
    assert anchor_of({"name": "--margin (default)", "line": 1}) == "margin"
    new_line, _ = relocate(_entry(line=1, name="--margin (default)", value="6.9"), d)
    assert new_line == 2, f"expected the attribute line 2, got {new_line}"
