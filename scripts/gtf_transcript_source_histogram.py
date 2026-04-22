#!/usr/bin/env python3
"""Histogram of rustle/StringTie GTF transcript lines by `source` attribute (if present).

Usage:
  python3 scripts/gtf_transcript_source_histogram.py assembled.gtf

Helps locate which pipeline stages dominate output when sources are tagged
(flow, junction_path, ref_chain, guide:*, ...).
"""
from __future__ import annotations

import collections
import re
import sys

SRC_RE = re.compile(r'\bsource\s+"([^"]*)"', re.I)


def main() -> None:
    path = sys.argv[1] if len(sys.argv) > 1 else None
    if not path:
        print("usage: gtf_transcript_source_histogram.py file.gtf", file=sys.stderr)
        sys.exit(2)
    counts: collections.Counter[str] = collections.Counter()
    n_tx = 0
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[2] != "transcript":
                continue
            n_tx += 1
            attrs = parts[8] if len(parts) > 8 else ""
            m = SRC_RE.search(attrs)
            key = m.group(1) if m else "(none)"
            counts[key] += 1
    if n_tx == 0:
        print("no transcript rows found", file=sys.stderr)
        sys.exit(1)
    print(f"# transcripts={n_tx}")
    for k, c in counts.most_common():
        print(f"{c}\t{100.0 * c / n_tx:.2f}%\t{k}")


if __name__ == "__main__":
    main()
