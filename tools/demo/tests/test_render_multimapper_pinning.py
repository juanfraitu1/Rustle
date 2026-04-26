import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
TOOL = REPO / "tools" / "demo" / "render_multimapper_pinning.py"
FIXTURES = Path(__file__).parent / "fixtures"


def test_renders_pinning_for_one_read(tmp_path):
    out_png = tmp_path / "pinning.png"
    cp = subprocess.run([
        "python3", str(TOOL),
        "--vg-trace-log", str(FIXTURES / "vg_trace_sample.log"),
        "--read-id", "READ_42",
        "--output-png", str(out_png),
    ], capture_output=True, text=True)
    assert cp.returncode == 0, cp.stderr
    assert out_png.exists() and out_png.stat().st_size > 500
