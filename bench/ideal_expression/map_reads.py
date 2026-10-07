#!/usr/bin/env python3
"""Map OUT.fq with the shipped minimap2 command against the whole-genome splice index, in read-disjoint parts, one part per call (bench/sim.py `_map_in_parts`:
verified parts, merged into OUT.bam + .bai). Re-run until it prints "mapping complete".

    map_reads.py OUT INDEX [--parts 4] [--threads 4] [--max-parts-per-call 1]
"""
import argparse
import os
import sys
from types import SimpleNamespace

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import sim  # noqa: E402

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("out")
ap.add_argument("index")
ap.add_argument("--parts", type=int, default=4)
ap.add_argument("--threads", type=int, default=4)
ap.add_argument("--max-parts-per-call", type=int, default=1)
a = ap.parse_args()
done = sim._map_in_parts(a.out, a.index, SimpleNamespace(parts=a.parts, threads=a.threads, max_parts_per_call=a.max_parts_per_call))
print("mapping complete" if done else "mapping parts remain: re-run to continue")
