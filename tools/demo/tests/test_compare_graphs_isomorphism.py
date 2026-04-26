import subprocess
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
TOOL = REPO / "tools" / "demo" / "compare_graphs_isomorphism.py"
FIXTURES = Path(__file__).parent / "fixtures"


def test_isomorphic_three_copies_produces_bijection(tmp_path):
    # All 3 copies have the same shape; iso must be found.
    out_json = tmp_path / "iso.json"
    out_png = tmp_path / "panel.png"

    cp = subprocess.run([
        "python3", str(TOOL),
        "--nodes-tsv", str(FIXTURES / "iso_triple_nodes.tsv"),
        "--edges-tsv", str(FIXTURES / "iso_triple_edges.tsv"),
        "--locus", "copy1:0-10000",
        "--locus", "copy2:0-10000",
        "--locus", "copy3:0-10000",
        "--output-json", str(out_json),
        "--output-png", str(out_png),
    ], capture_output=True, text=True)

    assert cp.returncode == 0, cp.stderr
    data = json.loads(out_json.read_text())
    assert data["isomorphic"] is True
    # bijection: {pair_label: {a_node_idx: b_node_idx}}
    assert "bijection" in data
    assert all(len(v) >= 3 for v in data["bijection"].values())


def test_non_isomorphic_reports_partial(tmp_path):
    out_json = tmp_path / "iso.json"
    out_png = tmp_path / "panel.png"

    cp = subprocess.run([
        "python3", str(TOOL),
        "--nodes-tsv", str(FIXTURES / "non_iso_nodes.tsv"),
        "--edges-tsv", str(FIXTURES / "non_iso_edges.tsv"),
        "--locus", "copyA:0-10000",
        "--locus", "copyB:0-10000",
        "--output-json", str(out_json),
        "--output-png", str(out_png),
    ], capture_output=True, text=True)

    assert cp.returncode == 0, cp.stderr
    data = json.loads(out_json.read_text())
    assert data["isomorphic"] is False
    assert "diff" in data
