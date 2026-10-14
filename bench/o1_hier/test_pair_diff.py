import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import pair_diff as P  # noqa: E402

V3, V2 = ("3", "1.0000", "0.1000", "NA"), ("2", "1.0000", "0.0000", "NA")


def test_two_joined_singletons_are_the_only_allowed_change():
    before = {frozenset("abc"): V3, frozenset("gh"): V2}
    after = {frozenset("abc"): V3, frozenset("de"): V2, frozenset("gh"): V2}
    r = P.diff(before, after)
    assert (r["pairs_joined"], r["pairs_split"], r["other_changes"], r["unchanged"]) == (1, 0, 0, 2)
    assert (r["clusters_before"], r["clusters_after"], r["members_before"], r["members_after"]) == (2, 3, 5, 7)


def test_a_changed_member_set_is_an_other_change():
    r = P.diff({frozenset("abc"): V3}, {frozenset("ab"): V2, frozenset("cd"): V2})
    assert r["other_changes"] == 3 and r["pairs_joined"] == 0


def test_the_same_set_with_other_values_is_an_other_change():
    assert P.diff({frozenset("ab"): V2}, {frozenset("ab"): ("2", "0.5000", "0.0000", "NA")})["other_changes"] == 1


def test_a_pair_that_disappears_is_a_split():
    r = P.diff({frozenset("ab"): V2, frozenset("abc"[2:] + "xy"): V3}, {frozenset("cxy"): V3})
    assert (r["pairs_split"], r["other_changes"]) == (1, 0)


def test_load_groups_rows_by_cluster_and_ignores_ids(tmp_path):
    h = "cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n"
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    a.write_text(h + "MCL0\t2\t1.0000\t0.0000\tNA\tchr1\t1\t10\nMCL0\t2\t1.0000\t0.0000\tNA\tchr1\t21\t30\n")
    b.write_text(h + "MCL7\t2\t1.0000\t0.0000\tNA\tchr1\t21\t30\nMCL7\t2\t1.0000\t0.0000\tNA\tchr1\t1\t10\n")
    assert P.load(str(a)) == P.load(str(b)) == {frozenset({"chr1:1-10", "chr1:21-30"}): V2}


def test_joined_pairs_are_split_by_edge_weight(tmp_path):
    g = tmp_path / "g.tsv"
    g.write_text("a\tb\t0.300000\nd\te\t0.990000\nx\ty\t0.500000\n")
    before = {frozenset("xyz"): V3}
    after = {frozenset("xyz"): V3, frozenset("ab"): V2, frozenset("de"): V2}
    assert P.joined_bands(before, after, P.edge_weights(str(g))) == {"below_w_star": 1, "not_isolated": 0, "tie_regime": 1}


def test_uniform_cliques_below_one_are_counted(tmp_path):
    g = tmp_path / "g.tsv"
    g.write_text("a\tb\t0.900000\na\tc\t0.900000\nb\tc\t0.900000\n"   # uniform K3 below 1: counted
                 "d\te\t1.000000\nd\tf\t1.000000\ne\tf\t1.000000\n"   # uniform K3 at 1.0: exactly uniform, safe
                 "g\th\t0.900000\nh\ti\t0.900000\n"                   # an equal-weight path is not a clique
                 "j\tk\t0.900000\nj\tl\t0.800000\nk\tl\t0.900000\n")  # a clique with unequal weights
    assert P.symmetric_cliques(P.edge_weights(str(g))) == {"components": 1, "members": 3}


HDR = "cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n"
A, B, C = "chr1:1-10", "chr1:21-30", "chr1:41-50"


def write_clusters(path, clusters):
    """clusters: a list of member lists ('chrN:start-end') -> a clusters.tsv."""
    rows = []
    for i, members in enumerate(clusters):
        for m in members:
            chrom, rest = m.split(":")
            start, end = rest.split("-")
            rows.append(f"MCL{i}\t{len(members)}\t1.0000\t0.0000\tNA\t{chrom}\t{start}\t{end}\n")
    path.write_text(HDR + "".join(rows))


def test_a_joined_pair_with_another_edge_is_not_isolated(tmp_path):
    g = tmp_path / "g.tsv"
    g.write_text("a\tb\t0.300000\nb\tc\t0.300000\nx\ty\t0.300000\n")  # b has a second edge; x-y is alone
    before = {frozenset("zzz"): V3}
    after = {frozenset("zzz"): V3, frozenset("ab"): V2, frozenset("xy"): V2}
    # not_isolated is a flag ON TOP of the band counts: both joined pairs stay in their band
    assert P.joined_bands(before, after, P.edge_weights(str(g))) == {"below_w_star": 2, "not_isolated": 1}


def test_main_returns_the_verdict_as_a_status(tmp_path):
    before, same, gone, out = (tmp_path / n for n in ("before.tsv", "same.tsv", "gone.tsv", "out.json"))
    write_clusters(before, [[A, B]])
    write_clusters(same, [[A, B]])
    write_clusters(gone, [])  # the pair has disappeared: a split
    assert P.main([str(before), str(same), str(out)]) == 0
    assert json.loads(out.read_text())["verdict"] == "PASS"
    assert P.main([str(before), str(gone), str(out)]) == 1
    res = json.loads(out.read_text())
    assert (res["verdict"], res["pairs_split"]) == ("FAIL", 1)


def test_main_fails_a_joined_pair_that_is_not_an_isolated_component(tmp_path):
    before, after, g, out = (tmp_path / n for n in ("before.tsv", "after.tsv", "g.tsv", "out.json"))
    write_clusters(before, [])
    write_clusters(after, [[A, B]])  # the fix joined A-B ...
    g.write_text(f"{A}\t{B}\t0.300000\n{B}\t{C}\t0.300000\n")  # ... but B has a second edge in the dump
    assert P.main([str(before), str(after), str(out), str(g)]) == 1
    res = json.loads(out.read_text())
    assert res["joined_by_weight"] == {"below_w_star": 1, "not_isolated": 1}
    assert (res["pairs_split"], res["other_changes"], res["verdict"]) == (0, 0, "FAIL")  # only isolation fails it
    g.write_text(f"{A}\t{B}\t0.300000\n")  # the same join as an isolated two-node component
    assert P.main([str(before), str(after), str(out), str(g)]) == 0
    res = json.loads(out.read_text())
    assert (res["joined_by_weight"], res["verdict"]) == ({"below_w_star": 1, "not_isolated": 0}, "PASS")


def test_main_fails_a_joined_pair_with_no_edge_in_the_dump(tmp_path):
    before, after, g, out = (tmp_path / n for n in ("before.tsv", "after.tsv", "g.tsv", "out.json"))
    write_clusters(before, [])
    write_clusters(after, [[A, B]])  # the fix joined A-B ...
    g.write_text(f"{B}\t{C}\t0.300000\n")  # ... but the dump has no A-B edge (A has none at all)
    assert P.main([str(before), str(after), str(out), str(g)]) == 1
    res = json.loads(out.read_text())
    assert res["joined_by_weight"] == {"no_weight": 1, "not_isolated": 0}
    assert (res["pairs_split"], res["other_changes"], res["verdict"]) == (0, 0, "FAIL")  # only no_weight fails it
