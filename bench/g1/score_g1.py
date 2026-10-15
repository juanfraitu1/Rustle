"""Usage: python3 bench/g1/score_g1.py human[,gorilla]   (light; after `figures/_liftoff.py finder --species S`)

G1 (figure 8 prereg, docs/archive/2026-09/PREREG_liftoff_loci_2026-09-25.md §5) for the named species, with the
registered code (fig_loci._write_guided) into a scratch table; prints the G1 rows and the other rows of the table."""
import sys
from pathlib import Path
sys.path.insert(0, "/mnt/linuxdisk/home/juanfraitu/rustle_soto_subset/figures")
import figlib
import fig_loci
cfg = figlib.load_inputs(None)
cfg["fig8_guided_finder"] = "1"
out = Path("/mnt/linuxdisk/tmp/o1_g1/fig8_data")
fig_loci._write_guided(cfg, out, sys.argv[1].split(","), ["scratch table for G1 (o1_g1)"])
for r in figlib.read_table("fig8_guided", out):
    print("\t".join(str(r[k]) for k in ("species", "item", "identity_min", "sc_min", "n", "k", "frac", "ci_lo", "ci_hi", "status")))
