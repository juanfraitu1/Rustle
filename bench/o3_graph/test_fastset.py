import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


F = _load("fastset")


def test_reads_shared_by_worlds_are_mapped_once_by_sequence():
    worlds = {"INTACT": {"s1": "ACGT", "s2": "GGGG"}, "OUTSIDE": {"s9": "ACGT", "s3": "GGGG"}, "CROSS": {"s1": "ACGT", "s7": "TTTT"}}
    union, owners = F.union_reads(worlds)
    assert sorted(union.values()) == ["ACGT", "GGGG", "TTTT"]
    by_seq = {union[u]: sorted(v) for u, v in owners.items()}
    assert by_seq["ACGT"] == [("CROSS", "s1"), ("INTACT", "s1"), ("OUTSIDE", "s9")]
    assert by_seq["TTTT"] == [("CROSS", "s7")]


def test_a_union_alignment_is_handed_to_every_world_read_with_that_sequence_under_its_own_name():
    owners = {"u0": [("INTACT", "s1"), ("OUTSIDE", "s9")], "u1": [("CROSS", "s7")]}
    lines = ["u0\t0\tchr1\t11\t60\t4M\t*\t0\t0\tACGT\t*\tNM:i:0\n", "u0\t256\tchr2\t5\t0\t4M\t*\t0\t0\t*\t*\tNM:i:1\n",
             "u1\t16\tchr3\t7\t60\t4M\t*\t0\t0\tAAAA\t*\n"]
    out = F.split_union(lines, owners)
    assert [ln.split("\t")[:3] for ln in out["INTACT"]] == [["s1", "0", "chr1"], ["s1", "256", "chr2"]]
    assert [ln.split("\t")[0] for ln in out["OUTSIDE"]] == ["s9", "s9"] and out["OUTSIDE"][0].endswith("NM:i:0\n")
    assert [ln.split("\t")[0] for ln in out["CROSS"]] == ["s7"]


def test_a_local_reference_is_reused_only_for_the_family_set_it_was_built_for(tmp_path):
    fj = tmp_path / "f.json"
    fj.write_text('{"sets": {"s1": [1], "s2": [2]}}')
    F.write_local_source(str(tmp_path), str(fj), "s1")
    assert F.local_matches(str(tmp_path), str(fj), "s1") and not F.local_matches(str(tmp_path), str(fj), "s2")
    fj.write_text('{"sets": {"s1": [9], "s2": [2]}}')
    assert not F.local_matches(str(tmp_path), str(fj), "s1")


def _opts(**kw):
    import argparse
    base = dict(set_threads=4, threads=2, heavy_current=False, o2_distance="gapcomp", member_genome_check=True, all_arms=False,
                o3_procs=2, o2_genome="off", per_read_detectability=False, member_attribution=False, union_exact=False,
                o2_procs=2, member_attribution_delta=False, twins=False, error_floor=False, o1_hosts=False)
    base.update(kw)
    return argparse.Namespace(**base)


def test_a_world_runs_current_then_o1_o2_and_o3_and_the_genome_step_only_when_asked():
    steps = [s for s, _, _ in F.chain_steps(_opts(), "/R", "INTACT")]
    assert steps == ["current", "O1", "O2", "O3"]
    for mode, method in (("fallback", "integl"), ("combined", "integc")):
        steps = F.chain_steps(_opts(o2_genome=mode), "/R", "INTACT")
        assert [s for s, _, _ in steps] == ["current", "O1", "O2", "O2_genome", "O3"]
        assert steps[3][1] == f"python3 bench/o3_graph/o2graph.py run /R/INTACT {method} 2" and steps[3][2] == "light"


def test_the_set_is_scored_on_the_genome_step_s_table_and_per_read_when_asked():
    assert F.e2e_command(_opts(), "/R", ["A", "B"]) == "python3 bench/o3_graph/summary.py e2e /R/e2e.json /R/A /R/B"
    assert F.e2e_command(_opts(o2_genome="fallback", per_read_detectability=True), "/R", ["A"]) == (
        "python3 bench/o3_graph/summary.py e2e /R/e2e.json /R/A --o2 o2.integl.m2.tsv --per-read-detectability")
    assert F.e2e_command(_opts(o2_genome="combined"), "/R", ["A"]) == (
        "python3 bench/o3_graph/summary.py e2e /R/e2e.json /R/A --o2 o2.integc.m2.tsv")


def test_o3_attributes_its_call_members_only_when_asked():
    assert "--member-attribution" not in F.chain_steps(_opts(), "/R", "INTACT")[-1][1]
    assert F.chain_steps(_opts(member_attribution=True), "/R", "INTACT")[-1][1].endswith(
        "detect.py /R/INTACT --member-genome-check --member-attribution --arms I")


def test_o2_takes_the_exact_union_only_when_asked():
    assert F.chain_steps(_opts(), "/R", "INTACT")[2][1] == ("env O2_PROCS=2 python3 bench/o3_graph/o2graph.py run /R/INTACT integ 2 "
                                                         "--distance gapcomp")
    assert F.chain_steps(_opts(union_exact=True), "/R", "INTACT")[2][1].endswith("integ 2 --distance gapcomp --union-exact")


def test_o3_scales_attribution_by_delta_only_when_asked():
    assert "--member-attribution-delta" not in F.chain_steps(_opts(), "/R", "INTACT")[-1][1]
    assert "--member-attribution-delta --arms I" in F.chain_steps(_opts(member_attribution_delta=True), "/R", "INTACT")[-1][1]


def test_o1_marks_twins_only_when_asked():
    assert F.chain_steps(_opts(), "/R", "INTACT")[1][1] == "python3 bench/o3_graph/o1seed.py /R/INTACT"
    assert F.chain_steps(_opts(twins=True), "/R", "INTACT")[1][1] == "python3 bench/o3_graph/o1seed.py /R/INTACT --twins"


def test_o3_subtracts_the_error_floor_only_when_asked():
    assert "--error-floor" not in F.chain_steps(_opts(), "/R", "INTACT")[-1][1]
    assert "--error-floor --arms I" in F.chain_steps(_opts(error_floor=True), "/R", "INTACT")[-1][1]


def test_o3_takes_its_hosts_from_o1_only_when_asked():
    assert "--o1-hosts" not in F.chain_steps(_opts(), "/R", "INTACT")[-1][1]
    assert "--o1-hosts --arms I" in F.chain_steps(_opts(o1_hosts=True), "/R", "INTACT")[-1][1]


def test_the_error_floor_reaches_o2_o3_and_the_scorer():
    steps = dict((n, c) for n, c, _ in F.chain_steps(_opts(error_floor=True, o2_genome="fallback"), "/R", "INTACT"))
    assert steps["O2"].endswith("--error-floor") and steps["O2_genome"].endswith("integl 2 --error-floor")
    assert "--error-floor" in steps["O3"]
    assert F.e2e_command(_opts(error_floor=True), "/R", ["A"]).endswith("/R/A --error-floor")
