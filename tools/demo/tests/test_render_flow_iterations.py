import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
TOOL = REPO / "tools" / "demo" / "render_flow_iterations.py"
FIXTURES = Path(__file__).parent / "fixtures"


def test_renders_iteration_grid(tmp_path):
    out_png = tmp_path / "flow.png"
    cp = subprocess.run([
        "python3", str(TOOL),
        "--flow-tsv", str(FIXTURES / "flow_3iter.tsv"),
        "--nodes-tsv", str(FIXTURES / "small_nodes.tsv"),
        "--edges-tsv", str(FIXTURES / "small_edges.tsv"),
        "--locus", "chr_test:0-20000",
        "--output-png", str(out_png),
    ], capture_output=True, text=True)
    assert cp.returncode == 0, cp.stderr
    assert out_png.exists() and out_png.stat().st_size > 1000
