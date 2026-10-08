#!/usr/bin/env python3
"""Where do a consensus and its true transcript differ at the ends? (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 10, reading S4.)

end_offsets(ops) -> (left, right): from a global alignment (ops = [(op, n)], consensus = query, transcript = target; op in = X I D, I = consensus has extra bases),
the net consensus-minus-transcript base count before the first run of >= ANCHOR matches at each end. Positive = the consensus carries extra bases there,
negative = it lacks transcript bases."""
import re

ANCHOR = 8


def parse_edlib_cigar(s):
    return [(op, int(n)) for n, op in re.findall(r"(\d+)([=XID])", s)]


def _offset(ops):
    q = t = run = 0
    start = (0, 0)
    for op, n in ops:
        if op == "=":
            if run == 0:
                start = (q, t)
            run += n
            q += n
            t += n
            if run >= ANCHOR:
                return start[0] - start[1]
            continue
        run = 0
        if op in "XI":
            q += n
        if op in "XD":
            t += n
    return None


def end_offsets(ops):
    return _offset(ops), _offset(list(reversed(ops)))


def _anchor_index(ops):
    """index of the first op of the first run of >= ANCHOR matches (adjacent '=' ops merged), or None"""
    run, first = 0, None
    for i, (op, n) in enumerate(ops):
        if op == "=":
            if run == 0:
                first = i
            run += n
            if run >= ANCHOR:
                return first
        else:
            run = 0
    return None


def core_identity(ops):
    """identity between the first and the last run of >= ANCHOR matches: 1 - (X + I + D) / columns, the end gaps excluded"""
    a = _anchor_index(ops)
    rev = list(reversed(ops))
    b = _anchor_index(rev)
    if a is None or b is None:
        return None
    core = ops[a:len(ops) - b]
    cols = sum(n for _, n in core)
    edits = sum(n for op, n in core if op != "=")
    return 1 - edits / cols if cols else None


def compare(cons, transcript):
    """global (edlib NW) comparison -> dict(identity, edits, left, right, qlen, tlen). Needs edlib (miniforge python); its path CIGAR is extended (=XID)."""
    import edlib
    r = edlib.align(cons, transcript, mode="NW", task="path")
    ops = parse_edlib_cigar(r["cigar"])
    left, right = end_offsets(ops)
    cols = sum(n for _, n in ops)
    return dict(identity=1 - r["editDistance"] / max(1, cols), core_identity=core_identity(ops), edits=r["editDistance"], left=left, right=right, qlen=len(cons), tlen=len(transcript))
