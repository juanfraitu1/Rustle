"""Tests for paf_units.py (docs/archive/2026-10/merged/PREREG_paf_duplicon_units_2026-10-10.md)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paf_units as P  # noqa: E402


def aln(qc, qs, qe, tc, ts, te, strand, cg, ident=1.0, qn=None, tn=None):
    return P.Aln(qc, qs, qe, tc, ts, te, strand, P.cigar_ops(cg), ident, max(qe - qs, te - ts), qn or qc, tn or tc)


def test_parse_region_is_one_based_inclusive():
    assert P.parse_region("chr1:112072-113065") == ("chr1", 112071, 113065)


def test_fields_become_genomic():
    f = ["chr1:1001-5000", "4000", "10", "2010", "+", "chr2:501-3000", "2500", "100", "2100", "1990", "2000", "60", "cg:Z:2000M"]
    a = P.aln_from_fields(f)
    assert (a.qc, a.qs, a.qe, a.tc, a.ts, a.te, a.strand) == ("chr1", 1010, 3010, "chr2", 600, 2600, "+")
    assert abs(a.ident - 0.995) < 1e-9 and a.block == 2000


def test_project_plus_strand_through_indels():
    a = aln("q", 100, 132, "t", 1000, 1033, "+", "10M2I10M3D10M")
    assert P.project(a, "q", 105) == ("t", 1005)
    assert P.project(a, "q", 111) == ("t", 1010)      # inside the insertion: collapses
    assert P.project(a, "q", 115) == ("t", 1013)
    assert P.project(a, "t", 1021) == ("q", 122)      # inside the deletion: collapses
    assert P.project(a, "q", 132) == ("t", 1033)


def test_project_minus_strand_reverses_query():
    a = aln("q", 100, 120, "t", 1000, 1020, "-", "20M")
    assert P.project(a, "q", 105) == ("t", 1015)
    assert P.project(a, "t", 1015) == ("q", 105)
    assert P.project(a, "q", 100) == ("t", 1020)


def test_breakpoints_snap_within_tolerance():
    bps = {}
    assert P.add_bp(bps, "c", 1000, 1000)
    assert not P.add_bp(bps, "c", 1500, 1000)
    assert P.add_bp(bps, "c", 2000, 1000)
    assert bps["c"] == [1000, 2000]


def test_keep_filters():
    ex = {"a": [(0, 10_000)], "b": [(0, 10_000)], "L1": [(0, 10_000)], "L2": [(0, 10_000)]}
    ok = aln("a", 0, 2000, "b", 0, 2000, "+", "2000M", 0.95)
    assert P.keep(ok, ex)
    assert P.keep(aln("a", 0, 2000, "b", 0, 2000, "+", "2000M", 0.90), ex)                  # the bar is inclusive
    assert not P.keep(aln("a", 0, 2000, "b", 0, 2000, "+", "2000M", 0.89), ex)
    assert not P.keep(aln("a", 0, 999, "b", 0, 999, "+", "999M", 0.95), ex)
    assert not P.keep(aln("c", 100, 2100, "c", 1000, 3000, "+", "2000M", 0.99, "L1", "L2"), ex)   # same DNA
    assert not P.keep(aln("c", 0, 2000, "d", 0, 2000, "+", "2000M", 0.99, "L1", "L1"), ex)        # same locus
    assert not P.keep(ok, {"a": [(5000, 6000)], "b": [(0, 10_000)]})                              # no exon on one side


def test_chain_xyz_gives_two_units():
    """X-Y aligned over AB, Y-Z over B only: units {A}, {B}; not one unit for everything."""
    alns = [aln("X", 0, 4000, "Y", 0, 4000, "+", "4000M"), aln("Y", 2000, 4000, "Z", 0, 2000, "+", "2000M")]
    bps = P.initial_breakpoints(alns, 1000)
    rounds, fix = P.propagate(alns, bps, 1000, 20)
    assert fix and bps["X"] == [0, 2000, 4000]
    at = P.atoms(P.covered(alns), bps, 1000)
    u = dict(zip(at, P.join_units(alns, at)))
    assert u[("X", 0, 2000)] == u[("Y", 0, 2000)]
    assert u[("X", 2000, 4000)] == u[("Y", 2000, 4000)] == u[("Z", 0, 2000)]
    assert u[("X", 0, 2000)] != u[("X", 2000, 4000)]


def test_join_needs_reciprocal_majority():
    a = aln("A", 0, 3000, "B", 0, 3000, "+", "3000M")
    at = [("A", 0, 3000), ("B", 0, 1200), ("B", 1200, 3000)]
    u = P.join_units([a], at)
    assert u[0] == u[2] and u[0] != u[1]


def test_short_atoms_are_not_units():
    at = P.atoms({"c": [(0, 3500)]}, {"c": [0, 1000, 1500, 3500]}, 1000)
    assert at == [("c", 0, 1000), ("c", 1500, 3500)]


def test_bed_lines_number_units_by_first_atom():
    at = [("c", 5000, 7000), ("c", 0, 2000), ("d", 0, 2000)]
    lines = P.bed_lines(at, [7, 3, 7])
    assert lines == ["c\t0\t2000\tPU1", "c\t5000\t7000\tPU2", "d\t0\t2000\tPU2"]
