#!/usr/bin/env python3
"""The registered identity x coverage and its gate-aware version (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 17).

ident x (aligned query bases) / (consensus length). Gate-aware: when the library gate is open and the consensus' unaligned 5' end is 1 to 3 bases, all G, those bases
leave the denominator. Nothing else changes."""

MAX_G = 3


def registered(ident, qs, qe, qlen):
    return ident * (qe - qs) / qlen


def gate_aware(ident, qs, qe, qlen, prefix, gate_open, max_g=MAX_G):
    """prefix = the consensus' first qs bases (the unaligned 5' end)"""
    if gate_open and 0 < qs <= max_g and len(prefix) == qs and set(prefix) == {"G"}:
        return ident * (qe - qs) / (qlen - qs)
    return registered(ident, qs, qe, qlen)
