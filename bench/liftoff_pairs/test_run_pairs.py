import json
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_pairs as RP  # noqa: E402

REAL = all(os.path.exists(p) for p in (RP.LIFTOFF, RP.RUNS, RP.FIG8))


def write(path, header, rows, comments=()):
    with open(path, "w") as fh:
        fh.writelines(f"# {c}\n" for c in comments)
        fh.write("\t".join(header) + "\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)


def toy(tmp_path, f1_k=1):
    """S1's placement and extra copy share family F1; S2 (family F2) and S3 (no family) are the decoys."""
    root, runs = tmp_path / "liftoff", tmp_path / "runs"
    (root / "human").mkdir(parents=True)
    (root / "support").mkdir()
    (runs / "human_testis").mkdir(parents=True)
    loci = [("S1", "in_place", 1000, 2000), ("S1", "extra_copy", 50000, 51000), ("S2", "in_place", 60000, 61000),
            ("S3", "in_place", 500000, 501000)]
    write(root / "human" / "liftoff_loci.tsv",
          ["species", "source_id", "cls", "contig", "start0", "end", "strand", "sequence_id", "exons", "exonic_len", "short"],
          [["human", s, c, "chr1", a, b, "+", "1.0", f"{a}-{b}", b - a, 0] for s, c, a, b in loci])
    write(root / "support" / "human_testis.tsv", ["source_id", "cls", "contig", "start0", "end", "reads_capped"],
          [[s, c, "chr1", a, b, 5] for s, c, a, b in loci])
    write(runs / "human_testis" / "human_testis.fam.copies.tsv",
          ["family_id", "copy_idx", "chrom", "start", "end", "strand", "exons"],
          [["F1", 0, "chr1", 1000, 2000, "+", "1000-2000"], ["F1", 1, "chr1", 50000, 51000, "+", "50000-51000"],
           ["F2", 0, "chr1", 60000, 61000, "+", "60000-61000"], ["F2", 1, "chr5", 1, 1001, "+", "1-1001"]])
    fig8 = tmp_path / "fig8_samples.tsv"
    write(fig8, ["sample", "species", "rustle_set", "scope", "measure", "stratum", "sc_min", "n", "k", "frac", "ci_lo",
                 "ci_hi", "status"],
          [["human_testis", "human", "families", scope, "pair_cofamily", "all", sc, 1, f1_k, "", "", "", "ok"]
           for scope in ("genome", "minus_dev") for sc in ("0.95", "0.98", "0.99", "1.00")],
          comments=["table: fig8_samples"])
    return ["--liftoff", str(root), "--runs", str(runs), "--fig8", str(fig8)]


def load(path):
    with open(path) as fh:
        return json.load(fh)


def test_inputs_only(tmp_path):
    out = tmp_path / "out"
    assert RP.main(["human_testis", "--out", str(out), "--inputs-only", *toy(tmp_path)]) == 0
    r = load(out / "inputs.json")
    assert r["f1"]["genome"]["0.95"] == [1, 1] and r["pairs"]["genome"]["0.95"] == 1 and r["pool"]["genome"] == 4
    assert sorted(os.listdir(out)) == ["inputs.json"]


def test_f1_mismatch_stops_the_run(tmp_path):
    with pytest.raises(SystemExit, match="F1 does not reproduce"):
        RP.main(["human_testis", "--out", str(tmp_path / "out"), "--inputs-only", *toy(tmp_path, f1_k=0)])


def test_full_toy_run(tmp_path):
    args = toy(tmp_path)
    for run in ("a", "b"):
        assert RP.main(["human_testis", "--out", str(tmp_path / run), "--null", "20", *args]) == 0
    t = load(tmp_path / "a" / "test.json")
    g = t["genome"]["0.95"]
    assert (g["n"], g["k"], g["null_max"], g["levels"]["0"]) == (1, 1, 0, 1) and g["p"] == 1 / 21
    assert t["primary"] == {"scope": "genome", "sc_min": "0.95"} and t["passes"] is False  # p 1/21 is not < 0.001
    with open(tmp_path / "a" / "test.json") as x, open(tmp_path / "b" / "test.json") as y:
        assert x.read() == y.read()


def test_unknown_sample_is_refused(tmp_path):
    with pytest.raises(SystemExit):
        RP.main(["chimp_PTR", "--out", str(tmp_path / "out")])


@pytest.mark.skipif(not REAL, reason="the figure 8 inputs are not on this machine")
def test_real_inputs_reproduce_f1(tmp_path):
    expect = {"human_A119b": ([148, 131], 514), "human_testis": ([34, 34], 38), "gorilla_KB3781": ([12, 11], 78)}
    for sid, (f1, n_pairs) in expect.items():
        assert RP.main([sid, "--out", str(tmp_path / sid), "--inputs-only"]) == 0
        r = load(tmp_path / sid / "inputs.json")
        assert r["f1"]["genome"]["0.95"] == f1 and r["pairs"]["genome"]["0.95"] == n_pairs, sid
