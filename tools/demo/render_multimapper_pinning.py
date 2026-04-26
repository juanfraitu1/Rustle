#!/usr/bin/env python3
"""Multi-mapper read pinning visualizer.

Parses [VG] stderr lines, picks one specified read, plots its weight at
each copy: BEFORE EM (uniform 1/NH) vs AFTER EM convergence.

Usage:
    render_multimapper_pinning.py --vg-trace-log path.log --read-id READ_X --output-png out.png
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

LINE_RE = re.compile(
    r"^\[VG\]\s+em_iter=(\d+)\s+read_id=(\S+)\s+copy=(\S+)\s+weight=([0-9.eE+-]+)"
)


def parse_log(path: Path, read_id: str) -> dict[int, dict[str, float]]:
    by_iter: dict[int, dict[str, float]] = defaultdict(dict)
    with path.open() as f:
        for line in f:
            m = LINE_RE.match(line.strip())
            if not m:
                continue
            it, rid, copy, w = int(m[1]), m[2], m[3], float(m[4])
            if rid != read_id:
                continue
            by_iter[it][copy] = w
    return dict(by_iter)


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vg-trace-log", required=True, type=Path)
    ap.add_argument("--read-id", required=True)
    ap.add_argument("--output-png", required=True, type=Path)
    args = ap.parse_args(argv)

    iters = parse_log(args.vg_trace_log, args.read_id)
    if 0 not in iters:
        raise SystemExit(f"no em_iter=0 row for read {args.read_id}")
    last_iter = max(iters.keys())

    pre = iters[0]
    post = iters[last_iter]
    copies = sorted(set(pre) | set(post))

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, label, weights in [(axes[0], "Before EM\n(StringTie-equivalent: uniform 1/NH)", pre),
                                (axes[1], f"After EM (iter {last_iter})\nRustle reweighted", post)]:
        vals = [weights.get(c, 0.0) for c in copies]
        bars = ax.bar(copies, vals, color="#1F4E79")
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("read weight")
        ax.set_title(label, fontsize=10)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.02, f"{v:.2f}",
                    ha="center", fontsize=9)
    fig.suptitle(f"Multi-mapper pinning: read = {args.read_id}", fontsize=12)
    fig.tight_layout()
    fig.savefig(args.output_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
