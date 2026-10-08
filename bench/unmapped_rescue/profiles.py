#!/usr/bin/env python3
"""Per-family HMM profile arm (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 5). Pure parsers here; the pipeline is `run_profiles.py`."""
import collections

try:
    from attribute import union_length
except ImportError:                      # pragma: no cover
    from .attribute import union_length


def parse_tblout(lines, max_evalue=1e-3):
    """nhmmer --tblout lines -> [(target, query family HMM, start, end, e-value, score)] with the target span normalised to start <= end"""
    out = []
    for ln in lines:
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.split()
        a, b = int(f[6]), int(f[7])
        ev, sc = float(f[12]), float(f[13])
        if ev <= max_evalue:
            out.append((f[0], f[2], min(a, b), max(a, b), ev, sc))
    return out


def cover_scores_from_rows(rows):
    """{target: {family: bases covered by the union of the family's hit spans}}"""
    iv = collections.defaultdict(lambda: collections.defaultdict(list))
    for t, fam, a, b, _ev, _sc in rows:
        iv[t][fam].append((a, b))
    return {t: {f: union_length(v) for f, v in fams.items()} for t, fams in iv.items()}


def bit_scores_from_rows(rows):
    """{target: {family: best hit bit score}}"""
    out = collections.defaultdict(dict)
    for t, fam, _a, _b, _ev, sc in rows:
        out[t][fam] = max(out[t].get(fam, 0.0), sc)
    return dict(out)
