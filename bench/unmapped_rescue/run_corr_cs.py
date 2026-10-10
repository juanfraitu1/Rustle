#!/usr/bin/env python3
"""Amendment 46: re-classify the Amendment 45 clusters with the support test on the CORRECTED reads (consensus and alignment reused). Output OUT/net_run_corr_cs/."""
import csv
import collections
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_net as RN  # noqa: E402
import seeds as SD  # noqa: E402

OUT = "/mnt/linuxdisk/tmp/o3_rescue/mattruth"
SRC = f"{OUT}/net_run_corr"
R = f"{OUT}/net_run_corr_cs"


def main():
    os.makedirs(R, exist_ok=True)
    for f in ("cons.jsonl", "cons.fa", "cons.R.paf", "cons.R.paf.done", "clusters.tsv"):
        if not os.path.exists(f"{R}/{f}"):
            shutil.copy(f"{SRC}/{f}", f"{R}/{f}")
    corr = SD.read_fa(RN.NET)
    corr.update(SD.read_fa(f"{SRC}/corrected.fa"))
    cl = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{R}/clusters.tsv"), delimiter="\t"):
        cl[r["cluster"]].append(r["read"])
    RN.finish(R, dict(cl), corr, support_seqs=corr)


if __name__ == "__main__":
    main()
