#!/usr/bin/env python3
"""Tests for the non-coding-aware promotion track (bench/promote_noncoding.py).
Hermetic: the module must import without pysam/pyabpoa/minimap2.
Run: PYTHONHASHSEED=0 pytest bench/test_promote_noncoding.py -v"""
import os
import sys
sys.path.insert(0, "bench")
import promote_noncoding as PN


def _rec(**kw):
    base = dict(genome_id=0.93, genome_cov=1.02, own_locus=True,
                alt_cols=82, alt_read_fraction=0.318, alt_reads=61, n_primary=192)
    base.update(kw)
    return base


def test_flagship_promotes():
    promote, call, reason = PN.classify_noncoding(_rec())
    assert promote is True
    assert call == "noncoding-candidate"
    assert reason  # non-empty


def test_artifact_low_identity_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(genome_id=0.28))
    assert promote is False and call == "artifact"


def test_near_reference_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(genome_id=0.99))
    assert promote is False and call == "~REF"


def test_chimeric_low_coverage_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(genome_cov=0.5))
    assert promote is False and call == "chimera"


def test_not_own_locus_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(own_locus=False))
    assert promote is False and call == "not-own-locus"


def test_thin_columns_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(alt_cols=2))
    assert promote is False and call == "thin"


def test_unbalanced_fraction_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(alt_read_fraction=0.05))
    assert promote is False and call == "thin"


def test_low_depth_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(alt_reads=3))
    assert promote is False and call == "thin"


def test_boundaries_inclusive():
    # genome_id == GID_LO promotes; == GID_HI is ~REF (excluded)
    assert PN.classify_noncoding(_rec(genome_id=PN.GID_LO))[0] is True
    assert PN.classify_noncoding(_rec(genome_id=PN.GID_HI))[1] == "~REF"
    # genome_cov == COV_MIN promotes; alt_cols == MIN_COLS promotes
    assert PN.classify_noncoding(_rec(genome_cov=PN.COV_MIN))[0] is True
    assert PN.classify_noncoding(_rec(alt_cols=PN.MIN_COLS))[0] is True
    # alt_read_fraction at both band edges promotes
    assert PN.classify_noncoding(_rec(alt_read_fraction=PN.AF_LO))[0] is True
    assert PN.classify_noncoding(_rec(alt_read_fraction=PN.AF_HI))[0] is True
