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


DELTA = 0.00958
BAR = 0.999


def o3_class(ref_best, truth_best, delta=DELTA, bar=BAR):
    """Amendment 25. R = identity x coverage of the consensus' best hit on the reference genome (None if no hit), T = the same on the truth genome (the other haplotype's
    assembly, or the genome that holds every copy). COPY: T >= bar and R more than the allele cutoff below 1; ALLELE: T >= bar and R within the cutoff but under the bar;
    PRESENT: R >= bar; None: the consensus is not reconstructed on the truth genome (T < bar) and not present either."""
    if ref_best is not None and ref_best >= bar:
        return "PRESENT"
    if truth_best is None or truth_best < bar:
        return None
    if ref_best is None or ref_best < 1 - delta - 1e-12:
        return "COPY"
    return "ALLELE"


def registered_flag(ref_best, truth_best, bar=BAR):
    """the registered O3 flag: reconstructed on the truth genome, not present on the reference"""
    return truth_best is not None and truth_best >= bar and (ref_best is None or ref_best < bar)


def discovery_class(ref_best, mat_best, pat_best, delta=DELTA, bar=BAR):
    """Amendment 28. CONFIRMED: more than the allele cutoff from the primary reference and >= bar on the mother's or the father's assembly; NOVEL: more than the cutoff
    from the primary and < bar on both; ALLELE-LIKE: within the cutoff of the primary but under the bar; PRESENT: >= bar on the primary."""
    if ref_best is not None and ref_best >= bar:
        return "PRESENT"
    if ref_best is not None and ref_best >= 1 - delta - 1e-12:
        return "ALLELE-LIKE"
    return "CONFIRMED" if max(mat_best or 0.0, pat_best or 0.0) >= bar else "NOVEL"


def gate_aware_total(ident, aligned, qlen, lead_unaligned, prefix, gate_open, max_g=MAX_G):
    """identity x (aligned bases / consensus length); when the gate is open and the leading unaligned end is 1 to 3 bases, all G, it leaves the denominator"""
    if gate_open and 0 < lead_unaligned <= max_g and len(prefix) == lead_unaligned and set(prefix) == {"G"}:
        return ident * aligned / (qlen - lead_unaligned)
    return ident * aligned / qlen
