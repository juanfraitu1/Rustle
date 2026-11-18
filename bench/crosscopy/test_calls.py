"""bench/crosscopy/calls.py: one call per key from PREFIX.fam.psi.tsv; classes, the two-bunch dedup, units."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from calls import load_calls, unit, flag_miskeyed, HEADER

def write(tmp_path, rows, name="f.tsv"):
    """rows: (family, member, anchor, col_a, col_o, chrom, strand, start, end, placed_w, bunch_placed_w, lo, hi, psi_placed_w, silent)"""
    out = []
    for fam, mem, anc, ca, co, chrom, st, s, e, pw, b, lo, hi, pp, sil in rows:
        d = dict.fromkeys(HEADER, "0")
        d.update(family=fam, member=mem, locus="G", anchor=anc, col_anchor=ca, col_other=co, class_anchor="core", class_other="core", chrom=chrom, strand=st,
                 start=s, end=e, placed_w=pw, max_w=pw, psi_placed_w=pp, psi_lo_w=lo, psi_hi_w=hi, bunch_placed_w=b, silent=sil, range_lo=lo, range_hi=hi, gap="NA")
        out.append("\t".join(str(d[h]) for h in HEADER))
    p = tmp_path / name
    p.write_text("\t".join(HEADER) + "\n" + "\n".join(out) + "\n")
    return str(p)

def test_classes_dedup_and_units(tmp_path):
    p = write(tmp_path, [
        ("F", "c1:1-500", "donor", "C0", "C2", "c1", "+", 101, 400, 12, 30, 0.4, 0.4, 0.4, 0),
        ("F", "c2:1-500", "donor", "C0", "C2", "c2", "-", 101, 400, 0, 40, 0.0, 0.0, 0.0, 1),
        ("F", "c1:1-500", "acceptor", "C2", "C0", "c1", "+", 101, 400, 12, 25, 0.3, 0.3, 0.3, 0),   # the same intron in its acceptor bunch
        ("F", "c2:1-500", "acceptor", "C2", "C0", "c2", "-", 101, 400, 0, 40, 0.0, 0.0, 0.0, 1),
        ("G", "c3:1-9", "donor", "C5", "C6", "c3", "+", 11, 90, 50, 50, 0.50, 0.52, 0.51, 0),
        ("G", "c4:1-9", "donor", "C5", "C6", "c4", "+", 11, 90, 60, 60, 0.51, 0.55, 0.53, 0),
        ("H", "c5:1-9", "donor", "C1", "C2", "c5", "+", 11, 90, 9, 9, 0.9, 0.9, 0.9, 0),            # 9 reads: not eligible
        ("H", "c6:1-9", "donor", "C1", "C2", "c6", "+", 11, 90, 50, 50, 0.1, 0.1, 0.1, 0)])
    C = load_calls(p)
    assert [(c.family, c.cls) for c in C] == [("F", "LARGE"), ("G", "SAME")]
    f = C[0]
    assert (f.anchor, f.gap, f.u, f.s, f.silent, f.copy_specific) == ("donor", 0.4, "c1:1-500", "c2:1-500", True, True)
    assert f.u_intron == ("c1", "+", 101, 400) and f.s_intron == ("c2", "-", 101, 400)
    assert unit(f, "GGO") == ("GGO", "c1", "+", 101, 400, "c2", "-", 400)

def test_small_overlap_and_shifted_sites(tmp_path):
    p = write(tmp_path, [
        ("S", "c1:1-9", "donor", "C0", "C1", "c1", "+", 11, 90, 30, 60, 0.50, 0.50, 0.5, 0),
        ("S", "c2:1-9", "donor", "C0", "C1", "c2", "+", 11, 90, 21, 60, 0.35, 0.40, 0.35, 0),        # gap 0.10: SMALL
        ("O", "c3:1-9", "donor", "C0", "C1", "c3", "+", 11, 90, 30, 60, 0.50, 0.60, 0.5, 0),
        ("O", "c4:1-9", "donor", "C0", "C1", "c4", "+", 11, 90, 21, 60, 0.30, 0.55, 0.35, 0),        # overlap, points 0.15 apart: no class
        ("T", "c5:1-9", "donor", "C0", "C1", "c5", "+", 11, 90, 10, 60, 0.9, 0.9, 0.17, 0),          # two shifted sites under one key at U:
        ("T", "c5:1-9", "donor", "C0", "C1", "c5", "+", 11, 93, 44, 60, 0.9, 0.9, 0.73, 0),          # the copy's range is the sum; U's intron = most reads
        ("T", "c6:1-9", "donor", "C0", "C1", "c6", "+", 11, 90, 12, 60, 0.2, 0.2, 0.2, 0)])
    C = {c.family: c for c in load_calls(p)}
    assert set(C) == {"S", "T"} and C["S"].cls == "SMALL" and C["S"].gap == 0.1 and not C["S"].copy_specific
    t = C["T"]
    assert t.cls == "LARGE" and t.gap == 0.7 and t.u_intron == ("c5", "+", 11, 93) and (t.reads_u, t.reads_s, t.lo_u, t.hi_s) == (60, 60, 0.9, 0.2)

def test_one_event_seen_from_both_copies_is_one_call(tmp_path):
    # copy c1 includes exon C1 at 0.6 and skips it at 0.4; c2 always skips: the inclusion key (U = c1) and the skip key (U = c2)
    # are the same event at the same donor, read from the two sides
    p = write(tmp_path, [
        ("M", "c1:1-500", "donor", "C0", "C1", "c1", "+", 101, 200, 60, 100, 0.6, 0.6, 0.6, 0),
        ("M", "c2:1-500", "donor", "C0", "C1", "c2", "+", 101, 200, 0, 100, 0.0, 0.0, 0.0, 1),
        ("M", "c1:1-500", "donor", "C0", "C2", "c1", "+", 101, 400, 40, 100, 0.4, 0.4, 0.4, 0),
        ("M", "c2:1-500", "donor", "C0", "C2", "c2", "+", 101, 400, 100, 100, 1.0, 1.0, 1.0, 0)])
    C = load_calls(p)
    assert len(C) == 1 and (C[0].col_other, C[0].u, C[0].s, C[0].gap) == ("C1", "c1:1-500", "c2:1-500", 0.6)

def test_same_boundary_two_junctions_one_donor_and_key_sets(tmp_path):
    p = write(tmp_path, [
        ("B", "c1:1-9", "donor", "C0", "C1", "c1", "+", 11, 90, 50, 50, 1.0, 1.0, 1.0, 0),          # points 1.0000 and 0.9500: SAME, at the boundary
        ("B", "c2:1-9", "donor", "C0", "C1", "c2", "+", 11, 90, 38, 40, 0.90, 1.0, 0.95, 0),
        ("D", "c3:1-900", "donor", "C0", "C1", "c3", "+", 101, 200, 50, 100, 0.5, 0.5, 0.5, 0),     # two junctions of one donor, both higher at c3:
        ("D", "c4:1-900", "donor", "C0", "C1", "c4", "+", 101, 200, 10, 100, 0.1, 0.1, 0.1, 0),     # two events, not one seen twice
        ("D", "c3:1-900", "donor", "C0", "C2", "c3", "+", 101, 400, 45, 100, 0.45, 0.45, 0.45, 0),
        ("D", "c4:1-900", "donor", "C0", "C2", "c4", "+", 101, 400, 10, 100, 0.1, 0.1, 0.1, 0),
        ("T", "c5:1-9", "donor", "C0", "C1", "c5", "+", 11, 90, 10, 60, 0.9, 0.9, 0.17, 0),
        ("T", "c5:1-9", "donor", "C0", "C1", "c5", "+", 11, 93, 44, 60, 0.9, 0.9, 0.73, 0),
        ("T", "c6:1-9", "donor", "C0", "C1", "c6", "-", 11, 90, 12, 60, 0.2, 0.2, 0.2, 0)])
    C = {(c.family, c.col_other): c for c in load_calls(p)}
    assert C[("B", "C1")].cls == "SAME" and {k for k in C if k[0] == "D"} == {("D", "C1"), ("D", "C2")}
    t = C[("T", "C1")]
    assert t.u_set == {("c5", "+", 11, 90), ("c5", "+", 11, 93)} and t.s_set == {("c6", "-", 11, 90)} and t.miskeyed is False
    assert unit(t, "G") == ("G", "c5", "+", 11, 93, "c6", "-", 90)                                  # S by its contig, strand and anchor site, not its locus name

def test_miskeyed_when_the_image_is_a_used_junction_under_another_key(tmp_path):
    p = write(tmp_path, [("K", "c1:1-900", "donor", "C0", "C1", "c1", "+", 101, 200, 90, 100, 0.9, 0.9, 0.9, 0),
                         ("K", "c2:1-900", "donor", "C0", "C1", "c2", "+", 101, 230, 25, 150, 0.17, 0.17, 0.17, 0)])
    rec = {("c1:1-900", "c2:1-900"): [([(0, 900, 0, 1)], 1, 1)]}                 # identity: c1's 101-200 is c2's 101-200
    c = load_calls(p)
    assert flag_miskeyed(c, {("c2", "+", 101, 200): [129] * 4}, rec)[0].miskeyed is True      # S uses that junction; the columns keyed it elsewhere
    assert flag_miskeyed(c, {("c2", "+", 101, 200): [2] * 4}, rec)[0].miskeyed is False       # two reads: not a used junction
    assert flag_miskeyed(c, {}, {})[0].miskeyed is False                                       # no image: nothing known against the call
