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


def test_o1_drops_are_scored_against_which_copies_were_silent():
    seeded = {"F": dict(copies=[dict(idx=0, origin="curated"), dict(idx=2, origin="curated")], dropped=[1, 3]),
              "G": dict(copies=[dict(idx=0, origin="curated")], dropped=[])}
    cat_copy = {("F", 0): "F.0", ("F", 1): "F.1", ("F", 2): "F.2", ("F", 3): "F.3", ("G", 0): "G.0"}
    kinds = {"F.0": "present", "F.1": "silent", "F.2": "present", "F.3": "present", "G.0": "present"}
    assert S.drop_outcome(seeded, cat_copy, kinds) == dict(silent=1, silent_dropped=1, expressed=4, expressed_dropped=1)


def test_a_read_s_fate_is_judged_by_what_kind_of_copy_it_came_from():
    # a catalogued (or O1-seeded uncatalogued) copy's read: right only on its own copy; an absent-copy call on it is an error
    assert S.read_fate("catalogued", assigned="F.0", own="F.0", in_o3=False) == "correct"
    assert S.read_fate("catalogued", assigned="F.1", own="F.0", in_o3=False) == "error"
    assert S.read_fate("catalogued", assigned=None, own="F.0", in_o3=True) == "error"
    assert S.read_fate("uncatalogued", assigned="F.5", own="F.5", in_o3=False) == "correct"
    assert S.read_fate("catalogued", assigned=None, own="F.0", in_o3=False) == "abstain"
    # an absent copy the genome does not explain: right only inside an absent-copy call
    assert S.read_fate("absent_detectable", assigned=None, own="F.9", in_o3=True) == "correct"
    assert S.read_fate("absent_detectable", assigned="F.1", own="F.9", in_o3=False) == "error"
    # an absent copy within delta of the genome cannot be told from the locus that explains it
    assert S.read_fate("absent_undetectable", assigned="F.1", own="F.9", in_o3=False) == "unknowable"
    assert S.read_fate("absent_undetectable", assigned=None, own="F.9", in_o3=True) == "correct"
