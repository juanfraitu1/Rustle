import json
import os
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import prepare as P  # noqa: E402
import score as S  # noqa: E402


def test_merge_regions_joins_overlapping_and_touching():
    lines = ["A\tc1:100-200\n", "B\tc1:150-300\n", "C\tc1:300-400\n", "D\tc2:5-10\n", "E\tc1:500-600\n"]
    assert P.merge_regions(lines) == ["c1:100-400", "c1:500-600", "c2:5-10"]


def test_original_seq_restores_reverse_records():
    fwd = types.SimpleNamespace(query_sequence="AACG", is_reverse=False)
    rev = types.SimpleNamespace(query_sequence="AACG", is_reverse=True)
    assert P.original_seq(fwd) == "AACG" and P.original_seq(rev) == "CGTT"


def test_aligned_blocks_skip_introns_deletions_and_insertions():
    rec = types.SimpleNamespace(reference_start=100, cigartuples=[(4, 5), (0, 10), (1, 3), (0, 5), (2, 2), (0, 5), (3, 50), (0, 10)])
    assert P.aligned_blocks(rec) == [(100, 110), (110, 115), (117, 122), (172, 182)]


def write(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)


def toy_arm(d, prefix, clusters, members, placed):
    write(f"{d}/{prefix}.assignments.tsv", ["read_name", "family_id", "status"], [(r, f, "assigned") for r, f in placed])
    write(f"{d}/{prefix}.psv_clusters.tsv", ["family", "host", "cluster", "n_reads", "n_sites", "n_alt", "n_unexplained_alt",
                                             "status", "explained_by", "est_div", "linkage", "editing_frac"],
          [(f, h, c, n, 3, 2, 2, st, "-", div, 1.0, 0.0) for f, h, c, n, st, div in clusters])
    write(f"{d}/{prefix}.psv_cluster_reads.tsv", ["read_name", "family", "host", "cluster"], members)


def test_score_removed_addedback_and_verdict(tmp_path):
    d = str(tmp_path)
    labels = [(f"a{i}", "A", "D", "A:2") for i in range(10)] + [(f"b{i}", "B", "D", "B:3") for i in range(10)] + \
             [("s1", "A", "S", "A:0")]
    write(f"{d}/labels.tsv", ["read", "family", "role", "copy"], labels)
    # masked arm: A's 10 D reads all placed on A, 9 in one unexplained cluster; B's D reads: 4 on B, 3 on A, 3 nowhere
    placed = [(f"a{i}", "A") for i in range(10)] + [("s1", "A")] + [(f"b{i}", "B") for i in range(4)] + \
             [(f"b{i}", "A") for i in range(4, 7)]
    toy_arm(d, "R", [("A", "0", "1", 9, "unexplained", 0.02), ("B", "0", "1", 4, "unexplained", 0.02)],
            [(f"a{i}", "A", "0", "1") for i in range(9)] + [(f"b{i}", "B", "0", "1") for i in range(4)], placed)
    assert S.main(["removed", f"{d}/labels.tsv", f"{d}/R", f"{d}/removed.tsv"]) == 0
    rows = {r["family"]: r for r in S.csv.DictReader(open(f"{d}/removed.tsv"), delimiter="\t")}
    assert (rows["A"]["recovered"], rows["A"]["label"], rows["A"]["n_placed"]) == ("1", "copy", "10")
    assert (rows["B"]["n_placed"], rows["B"]["on_other_family"], rows["B"]["not_placed"]) == ("4", "3", "3")
    # added back: A's reads explained (cleared); B still has an unexplained cluster with its reads
    toy_arm(d, "R0", [("A", "2", "1", 10, "explained", 0.0), ("B", "0", "1", 4, "unexplained", 0.02)],
            [(f"a{i}", "A", "2", "1") for i in range(10)] + [(f"b{i}", "B", "0", "1") for i in range(4)],
            [(f"a{i}", "A") for i in range(10)] + [(f"b{i}", "B") for i in range(4)])
    assert S.main(["addedback", f"{d}/labels.tsv", f"{d}/R0", f"{d}/addedback.tsv"]) == 0
    ab = {r["family"]: r["cleared"] for r in S.csv.DictReader(open(f"{d}/addedback.tsv"), delimiter="\t")}
    assert ab == {"A": "1", "B": "0"}
    rd = lambda p: list(S.csv.DictReader(open(p), delimiter="\t"))  # noqa: E731
    v = S.verdict(rd(f"{d}/removed.tsv"), rd(f"{d}/addedback.tsv"), None)
    # B has 4 placed reads (< 6): not eligible; A eligible, recovered and cleared
    assert (v["eligible"], v["recovered"], v["r1"], v["cleared"], v["r2"]) == (1, 1, 1.0, 1, 1.0)


def test_verdict_control_counts_random_stratum_only():
    ctl = [dict(stratum="random", bin="0", flag="-"), dict(stratum="random", bin=">delta", flag="copy"),
           dict(stratum="high_het", bin=">delta", flag="copy")] + [dict(stratum="random", bin="0", flag="-")] * 47
    v = S.verdict([], [], ctl)
    assert (v["control_genes"], v["false_copy"], v["r3"], v["r3_pass"]) == (49, 1, 1 / 49, False)


def test_verdict_by_copy_number(tmp_path):
    removed = [dict(family="A", n_placed="10", recovered="1", label="copy"), dict(family="B", n_placed="10", recovered="0", label="-"),
               dict(family="C", n_placed="10", recovered="1", label="copy")]
    added = [dict(family="A", cleared="1"), dict(family="C", cleared="0")]
    picked = {"A": 2, "B": 3, "C": 3}
    v = S.verdict(removed, added, None, copies_of=picked)
    assert v["by_copies"] == {"2": dict(eligible=1, recovered=1, cleared=1), ">=3": dict(eligible=2, recovered=1, cleared=0)}
