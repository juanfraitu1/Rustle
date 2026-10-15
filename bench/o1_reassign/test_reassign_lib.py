import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reassign_lib as R  # noqa: E402


def adj(*edges):
    out = {}
    for u, v, w in edges:
        out.setdefault(u, []).append((v, w))
        out.setdefault(v, []).append((u, w))
    return out


def test_read_dump_and_clusters(tmp_path):
    (tmp_path / "g.tsv").write_text("a:1-2\tb:3-4\t0.500000\nb:3-4\tc:5-6\t0.900000\n")
    a = R.read_dump(str(tmp_path / "g.tsv"))
    assert sorted(a["b:3-4"]) == [("a:1-2", 0.5), ("c:5-6", 0.9)]
    (tmp_path / "c.tsv").write_text("cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n"
                                    "MCL0\t2\t1\t0\tNA\ta\t1\t2\nMCL0\t2\t1\t0\tNA\tb\t3\t4\n")
    assert R.read_clusters(str(tmp_path / "c.tsv")) == {"a:1-2": "MCL0", "b:3-4": "MCL0"}


def test_node_follows_a_stronger_edge_out_of_its_family():
    # x sits in F with in = 0.96 but has 0.99 into G (the CYKILR shape)
    fam = {"x": "F", "y": "F", "g1": "G", "g2": "G"}
    a = adj(("x", "y", 0.96), ("x", "g1", 0.99), ("g1", "g2", 0.99), ("y", "g2", 0.50))
    new, log = R.reassign(fam, a)
    assert new["x"] == "G"
    assert log["moves"][0][1:5] == ("x", "F", "G", 0.96)
    assert log["converged"] is True


def test_ties_stay_and_weaker_outside_edges_do_nothing():
    fam = {"x": "F", "y": "F", "g1": "G", "g2": "G"}
    a = adj(("x", "y", 0.9), ("x", "g1", 0.9), ("g1", "g2", 0.95), ("y", "g2", 0.3))
    new, log = R.reassign(fam, a)
    assert new == fam and log["moves"] == [] and log["passes"] == 1


def test_a_member_left_alone_follows_its_strongest_edge():
    # x and z leave F for G; y is then alone in F (in = 0) and follows its own strongest edge into G
    fam = {"x": "F", "z": "F", "y": "F", "g1": "G", "g2": "G"}
    a = adj(("x", "z", 0.5), ("y", "x", 0.4), ("x", "g1", 0.99), ("z", "g2", 0.98), ("g1", "g2", 0.99),
            ("y", "g1", 0.3))
    new, log = R.reassign(fam, a)
    assert new["x"] == new["z"] == new["y"] == "G"
    assert log["lost"] == [] and log["converged"] is True


def test_a_member_alone_without_any_family_edge_is_reported_lost():
    # MCL put y in F without a direct edge to x; once x leaves, y is alone with edges only to unclustered nodes
    fam = {"x": "F", "y": "F", "g1": "G", "g2": "G"}
    a = adj(("x", "g1", 0.99), ("g1", "g2", 0.99), ("y", "lonely", 0.8))
    new, log = R.reassign(fam, a)
    assert new["x"] == "G" and new["y"] == "F" and log["lost"] == ["y"]


def test_edges_to_unclustered_nodes_are_ignored():
    fam = {"x": "F", "y": "F"}
    a = adj(("x", "y", 0.5), ("x", "lonely", 0.99))
    new, log = R.reassign(fam, a)
    assert new == fam and log["moves"] == []


def test_outside_tie_between_families_is_broken_by_family_name():
    fam = {"x": "F", "y": "F", "a1": "MCL9", "a2": "MCL9", "b1": "MCL10", "b2": "MCL10"}
    a = adj(("x", "y", 0.2), ("x", "a1", 0.9), ("x", "b1", 0.9), ("a1", "a2", 0.9), ("b1", "b2", 0.9))
    new, _ = R.reassign(fam, a)
    assert new["x"] == "MCL10"          # min by name: "MCL10" < "MCL9"


def test_reassign_is_deterministic_and_does_not_mutate_its_input():
    fam = {"x": "F", "y": "F", "g1": "G", "g2": "G"}
    a = adj(("x", "y", 0.96), ("x", "g1", 0.99), ("g1", "g2", 0.99))
    before = dict(fam)
    assert R.reassign(fam, a) == R.reassign(fam, a)
    assert fam == before


def test_pass_cap_reports_non_convergence():
    # pass 1 moves x (and y after it); a second pass is needed to see that nothing moves, so one pass is not enough
    fam = {"x": "F", "y": "F", "g1": "G", "g2": "G"}
    a = adj(("x", "y", 0.5), ("x", "g1", 0.99), ("g1", "g2", 0.99))
    _, log = R.reassign(fam, a, max_passes=1)
    assert log["converged"] is False


def test_label_gains_and_losses():
    old = {"x": "F", "y": "F", "g1": "G", "g2": "G"}
    new = {"x": "G", "y": "F", "g1": "G", "g2": "G"}
    labels = {"x": {"S1"}, "g1": {"S1"}, "g2": {"S1"}, "y": {"S2"}}
    assert R.label_moves(["x"], old, new, labels) == dict(gain=1, loss=0, same=0, unlabelled=0)
    labels2 = {"x": {"S2"}, "y": {"S2"}, "g1": {"S1"}}
    assert R.label_moves(["x"], old, new, labels2) == dict(gain=0, loss=1, same=0, unlabelled=0)
    assert R.label_moves(["x"], old, new, {}) == dict(gain=0, loss=0, same=0, unlabelled=1)


def test_pair_changes():
    before = [("s1", True, {"F"}), ("s2", True, set()), ("s3", True, {"G"}), ("s4", False, set())]
    after = [("s1", True, {"G"}), ("s2", True, {"G"}), ("s3", True, set()), ("s4", True, set())]
    assert R.pair_changes(before, after) == dict(joined=1, split=1, kept=1, apart=0, covered_before=3, covered_after=4)


def test_an_unclustered_neighbour_does_not_hide_a_real_family():
    fam = {"x": "F", "y": "F", "g1": "G", "g2": "G"}
    a = adj(("x", "y", 0.5), ("x", "lonely", 0.99), ("x", "g1", 0.8), ("g1", "g2", 0.9))
    new, _ = R.reassign(fam, a)
    assert new["x"] == "G"
