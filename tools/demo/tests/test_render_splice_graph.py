import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
TOOL = REPO / "tools" / "demo" / "render_splice_graph.py"
NODES_FIXTURE = Path(__file__).parent / "fixtures" / "small_nodes.tsv"
EDGES_FIXTURE = Path(__file__).parent / "fixtures" / "small_edges.tsv"


def test_render_produces_png_for_small_graph(tmp_path):
    out_png = tmp_path / "graph.png"
    out_dot = tmp_path / "graph.dot"

    cp = subprocess.run([
        "python3", str(TOOL),
        "--nodes-tsv", str(NODES_FIXTURE),
        "--edges-tsv", str(EDGES_FIXTURE),
        "--output-png", str(out_png),
        "--output-dot", str(out_dot),
        "--locus", "chr_test:0-20000",
    ], capture_output=True, text=True)

    assert cp.returncode == 0, cp.stderr
    assert out_png.exists() and out_png.stat().st_size > 100
    assert out_dot.exists()

    dot = out_dot.read_text()
    assert "digraph" in dot
    assert "label=" in dot
    # Three nodes + 2 edges in fixture, expect at least 3 node entries.
    assert dot.count("[shape=") >= 3
