#!/usr/bin/env python3
"""Amendment 43: what a flagged consensus differs from the primary in, from minimap2's long cs (--cs=long). Pure: cs_features."""
import re

TOK = re.compile(r"(=[ACGTN]+|\*[acgtn][acgtn]|\+[acgtn]+|-[acgtn]+|~[acgtn]{2}\d+[acgtn]{2})")


def cs_features(cs):
    """counts of matches, substitutions (and the a>g / t>c editing signature), insertions, deletions and their bases, and indels that only change the length of a
    homopolymer run of the reference (all indel bases equal to the reference base right before or right after)"""
    out = dict(match=0, sub=0, ag=0, ins=0, dele=0, ins_bp=0, del_bp=0, hp_indel=0)
    toks = TOK.findall(cs)
    prev_ref = ""                                    # the reference base just before the current position
    for i, t in enumerate(toks):
        c = t[0]
        if c == "=":
            out["match"] += len(t) - 1
            prev_ref = t[-1]
        elif c == "*":
            out["sub"] += 1
            if t[1:] in ("ag", "tc"):
                out["ag"] += 1
            prev_ref = t[1].upper()
        elif c in "+-":
            b = t[1:].upper()
            nxt = ""
            for u in toks[i + 1:]:
                if u[0] == "=":
                    nxt = u[1]
                    break
                if u[0] == "*":
                    nxt = u[1].upper()
                    break
                if u[0] == "-":
                    nxt = u[1].upper()
                    break
            if len(set(b)) == 1 and (b[0] == prev_ref or b[0] == nxt):
                out["hp_indel"] += 1
            if c == "+":
                out["ins"] += 1
                out["ins_bp"] += len(b)
            else:
                out["dele"] += 1
                out["del_bp"] += len(b)
                prev_ref = b[-1]
        else:  # intron
            prev_ref = ""
    return out
