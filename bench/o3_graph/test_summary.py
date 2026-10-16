import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S = _load("summary")


def test_a_removed_copy_is_detectable_when_its_reads_median_de_exceeds_delta():
    assert S.detectable([0.02, 0.01, 0.03]) is True
    assert S.detectable([0.002, 0.009, 0.0095]) is False
    assert S.detectable([]) is False


def test_the_o2_bars_compare_padded_union_with_the_current_and_locus_o2_per_world():
    worlds = {"INTACT": dict(cur_correct=100, graph_correct=99, graph_assigned=99, locus_correct=90),
              "REMOVED": dict(cur_correct=50, graph_correct=50, graph_assigned=50, locus_correct=45,
                              cur_assigned_noncat=20, graph_assigned_noncat=5),
              "BRANCH": dict(cur_correct=100, graph_correct=97, graph_assigned=97, locus_correct=90,
                             cur_assigned_noncat=30, graph_assigned_noncat=1)}
    b = S.o2_bars(worlds)
    assert b["A1"] == 0.97 and b["A2"] == 1.0
    assert b["A3_fewer"] is True and abs(b["A3_pooled"] - 6 / 50) < 1e-9
    assert abs(b["A4"] - 246 / 225) < 1e-9


def test_a_flagged_cluster_is_a_copy_an_allele_or_other():
    truth = {**{f"x{i}": dict(copy="F.9", hap="1") for i in range(20)},
             **{f"a{i}": dict(copy="F.1", hap="2") for i in range(19)}, "a19": dict(copy="F.1", hap="1"),
             **{f"m{i}": dict(copy="F.2", hap=str(1 + i % 2)) for i in range(20)}}
    nonc = {"F.9"}
    assert S.classify_cluster({f"x{i}" for i in range(20)}, truth, nonc)[0] == "copy"
    assert S.classify_cluster({f"a{i}" for i in range(20)}, truth, nonc)[0] == "allele"      # one haplotype of a catalogued copy
    assert S.classify_cluster({f"m{i}" for i in range(20)}, truth, nonc)[0] == "other"       # both haplotypes: not an allele split
