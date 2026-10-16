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
