"""Runs every tools/*.py --selftest (audit 2026-11-19: they existed but nothing ran them). `pixi run test-py`."""
import os, subprocess, sys
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = ["locus_reps", "piece_merge", "piece_attach", "reps_table", "gtf_ptc", "oneread", "check_reads", "gtf_exclude_bed",
         "copy_table_check", "multicopy_rescue", "o2_family_report", "read_provenance", "empty_place_bam"]

@pytest.mark.parametrize("tool", TOOLS)
def test_selftest(tool):
    r = subprocess.run([sys.executable, os.path.join(HERE, tool + ".py"), "--selftest"], capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, f"{tool}: rc {r.returncode}\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}"
