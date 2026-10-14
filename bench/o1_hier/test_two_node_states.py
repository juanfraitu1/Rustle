import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import pair_diff as P  # noqa: E402
import two_node_states as T  # noqa: E402

V2, V3 = ("2", "1.0000", "0.0000", "NA"), ("3", "1.0000", "0.1000", "NA")


def fs(*names):
    return frozenset(names)


def synthetic():
    """A graph with every fate; -> (clusters before, clusters after, loci before, loci after, edge weights)."""
    w = {fs(u, v): x for u, v, x in [
        ("a1", "a2", 0.90),                                                       # tie regime, a 2-member family both times
        ("b1", "b2", 0.95),                                                       # tie regime, one locus both times
        ("c1", "c2", 0.30),                                                       # below w*: split -> a 2-member family
        ("d1", "d2", 0.40),                                                       # below w*: split -> one locus
        ("e1", "e2", 0.35),                                                       # below w*: split -> split (the fix failed)
        ("f1", "f2", 0.80), ("f2", "f3", 0.80),                                   # a path of 3 nodes: not a two-node component
        ("g1", "g2", 0.90), ("g1", "g3", 0.90), ("g2", "g3", 0.90),               # uniform K3 below 1: a 3-member family
        ("h1", "h2", 0.95), ("h1", "h3", 0.95), ("h2", "h3", 0.95),               # uniform K3 below 1: one locus
        ("i1", "i2", 0.80), ("i1", "i3", 0.80), ("i2", "i3", 0.80),               # uniform K3 below 1: neither
        ("j1", "j2", 1.00), ("j1", "j3", 1.00), ("j2", "j3", 1.00),               # uniform K3 at 1.0: exactly uniform, not counted
        ("k1", "k2", 0.90), ("k1", "k3", 0.80), ("k2", "k3", 0.90),               # K3 with unequal weights: not counted
    ]}
    before = {fs("a1", "a2"): V2, fs("g1", "g2", "g3"): V3, fs("i1", "i2"): V2}
    after = {fs("a1", "a2"): V2, fs("c1", "c2"): V2, fs("g1", "g2", "g3"): V3, fs("i1", "i2"): V2}
    loci_before = {"b2": "b1", "h2": "h1", "h3": "h1", "z2": "z1"}                # z2 -> z1 disappears after the fix
    # d2 -> d1 is the fold of a joined pair; f3 -> f2 is not
    loci_after = {"b2": "b1", "h2": "h1", "h3": "h1", "d2": "d1", "f3": "f2"}
    return before, after, loci_before, loci_after, w


def test_every_two_node_component_gets_its_band_and_fate():
    r = T.states(*synthetic())
    assert r["two_node_components"] == 5  # a, b, c, d, e; the path f and every K3 are not two-node components
    assert r["by_band"] == {
        "tie_regime": {"total": 2, "fates": {"folded->folded": 1, "listed->listed": 1}},
        "below_w_star": {"total": 3, "fates": {"split->folded": 1, "split->listed": 1, "split->split": 1}},
    }
    assert r["still_split_after"] == 1  # e: a pair the fix left split


def test_loci_line_changes_are_counted_and_placed():
    r = T.states(*synthetic())
    assert r["loci_tsv"] == {
        "lines_before": 4, "lines_after": 5, "added": 2, "removed_or_changed": 1,
        "added_in_two_node_component": {"tie_regime": 0, "below_w_star": 1},  # d2 -> d1
        "added_elsewhere": 1,                                                 # f3 -> f2 lies in a 3-node component
    }


def test_uniform_cliques_below_one_get_their_fates():
    before, after, loci_before, loci_after, w = synthetic()
    cl = T.states(before, after, loci_before, loci_after, w)["uniform_cliques_below_1"]
    assert cl == [
        {"members": ["g1", "g2", "g3"], "weight": 0.9, "before": "listed", "after": "listed"},
        {"members": ["h1", "h2", "h3"], "weight": 0.95, "before": "folded", "after": "folded"},
        {"members": ["i1", "i2", "i3"], "weight": 0.8, "before": "other", "after": "other"},
    ]
    # the same definition as pair_diff's count (K3 at 1.0, unequal weights and the path are all excluded)
    assert P.symmetric_cliques(w) == {"components": 3, "members": 9}


def write_tsvs(root, sample, clusters, loci, graph=None):
    """One sample's `<S>.fam.clusters.tsv` / `<S>.fam.loci.tsv` (or `<S>.graph.tsv`) under root."""
    root.mkdir(exist_ok=True)
    rows = ["cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n"]
    for i, members in enumerate(clusters):
        for m in members:
            chrom, rest = m.split(":")
            start, end = rest.split("-")
            rows.append(f"MCL{i}\t{len(members)}\t1.0000\t0.0000\tNA\t{chrom}\t{start}\t{end}\n")
    (root / f"{sample}.fam.clusters.tsv").write_text("".join(rows))
    loci_rows = "".join(f"{k}\t{v}\n" for k, v in loci.items())
    (root / f"{sample}.fam.loci.tsv").write_text("annotation\trepresentative\n" + loci_rows)
    if graph is not None:
        (root / f"{sample}.graph.tsv").write_text("".join(f"{u}\t{v}\t{x:.6f}\n" for u, v, x in graph))


def test_main_reads_the_files_per_sample_and_totals_them(tmp_path):
    A, B, C, D = "chr1:1-10", "chr1:21-30", "chr1:41-50", "chr1:61-70"
    graph = [(A, B, 0.9), (C, D, 0.3)]  # one tie-regime pair and one below w*
    write_tsvs(tmp_path / "before", "s1", [[A, B]], {})
    write_tsvs(tmp_path / "after", "s1", [[A, B], [C, D]], {})                      # s1: C-D joined and listed
    write_tsvs(tmp_path / "graph", "s1", [], {}, graph)
    write_tsvs(tmp_path / "before", "s2", [[A, B]], {})
    write_tsvs(tmp_path / "after", "s2", [[A, B]], {D: C})                          # s2: C-D joined and folded into one locus
    write_tsvs(tmp_path / "graph", "s2", [], {}, graph)
    out = tmp_path / "states.json"
    assert T.main([str(tmp_path / "before"), str(tmp_path / "after"), str(tmp_path / "graph"), str(out), "s1", "s2"]) == 0
    res = json.loads(out.read_text())
    assert list(res["samples"]) == ["s1", "s2"]
    assert res["samples"]["s1"]["by_band"]["below_w_star"]["fates"] == {"split->listed": 1}
    assert res["samples"]["s2"]["by_band"]["below_w_star"]["fates"] == {"split->folded": 1}
    assert res["samples"]["s2"]["loci_tsv"]["added_in_two_node_component"] == {"tie_regime": 0, "below_w_star": 1}
    assert res["total"] == {
        "two_node_components": 4,
        "by_band": {
            "tie_regime": {"total": 2, "fates": {"listed->listed": 2}},
            "below_w_star": {"total": 2, "fates": {"split->folded": 1, "split->listed": 1}},
        },
        "still_split_after": 0,
        "loci_tsv": {"added": 1, "removed_or_changed": 0, "added_elsewhere": 0,
                     "added_in_two_node_component": {"tie_regime": 0, "below_w_star": 1}},
        "uniform_cliques_below_1": {"components": 0, "members": 0, "fates_before": {}, "fates_after": {}},
    }
