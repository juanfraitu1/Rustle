import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


O = _load("o2graph")
PATHS = {"copy0": {"1", "2", "3", "5"}, "copy1": {"1", "2", "4", "6"}, "copy2": {"1", "7", "4", "6"}}


def test_a_read_goes_to_the_copy_whose_path_holds_its_allele_nodes():
    # informative nodes: on some copies' paths but not all (node 1 is shared by all three)
    assert O.assign(["1", "2", "3", "5"], PATHS, margin=2) == ("copy0", {"copy0": 0, "copy1": 2, "copy2": 3})


def test_a_tie_or_a_thin_margin_abstains():
    assert O.assign(["1", "4", "6"], PATHS, margin=2)[0] is None           # copy1 and copy2 both hold every allele node
    assert O.assign(["1", "2", "4"], PATHS, margin=2)[0] is None           # copy1 misses 0, copy0 / copy2 miss 1: margin 1
    assert O.assign(["1"], PATHS, margin=2)[0] is None                     # no informative node at all


def test_paths_assign_the_copy_with_fewest_edits_by_the_margin():
    assert O.decide({"copy0": 1, "copy1": 4, "copy2": 6}, qlen=1000, margin=2) == ("copy0", "assigned")
    assert O.decide({"copy0": 1, "copy1": 2}, qlen=1000, margin=2) == (None, "thin")
    assert O.decide({"copy0": 3, "copy1": 3}, qlen=1000, margin=2) == (None, "tied")
    assert O.decide({"copy0": 1}, qlen=1000, margin=2) == ("copy0", "assigned")       # the only copy it aligns to


def test_a_read_no_copy_explains_within_delta_abstains():
    # 30 edits over 1,000 bases is 3% away from the best copy: an absent copy's read, O3's business
    assert O.decide({"copy0": 30, "copy1": 45}, qlen=1000, margin=2) == (None, "unexplained")
    assert O.decide({}, qlen=1000, margin=2) == (None, "unexplained")


def _world(tmp_path, cur_rows, o2_rows):
    (tmp_path / "cur").mkdir()
    (tmp_path / "graph").mkdir()
    (tmp_path / "P.fam.copies.tsv").write_text("family_id\tcopy_idx\ttid\nA\t0\tA.0\nA\t1\tA.1\nB\t0\tB.0\nB\t1\tB.1\n")
    (tmp_path / "reads.truth.E0.tsv").write_text("read\tcopy\tkind\tfamily\tisoform\tnovel\n"
                                                 "r1\tA.0\tpresent\tA\tt\t0\nr2\tB.0\tpresent\tB\tt\t0\nr3\tA.X\toutside\tA\tt\t0\n")
    (tmp_path / "cur" / "A.assignments.tsv").write_text(
        "read_name\tfamily_id\tstatus\tcatalog_copy_idx\n" + "".join("\t".join(r) + "\n" for r in cur_rows))
    (tmp_path / "graph" / "o2.t.tsv").write_text("read\tfamily\tcall\treason\n" + "".join("\t".join(r) + "\n" for r in o2_rows))
    return str(tmp_path)


def test_an_assignment_in_another_family_is_scored_wrong(tmp_path):
    # r1 (A.0): the current method names copy 0 in A (right) and copy 1 in B (a wrong-family call); r2 (B.0): path O2 names
    # copy 1 of A; r3 (an uncatalogued copy of A) assigned in B by the current method: wrong as well
    w = _world(tmp_path, [("r1", "A", "assigned", "0"), ("r1", "B", "assigned", "1"), ("r2", "B", "ambiguous", "-"),
                          ("r3", "B", "assigned", "0")],
               [("r1", "A", "0", "assigned"), ("r2", "A", "1", "assigned"), ("r2", "B", "-", "tied")])
    res = O.score_world(w, "o2.t.tsv")
    assert res["A"]["cur_correct"] == 1 and res["A"]["graph_correct"] == 1
    assert res["B"]["cur_wrong_family"] == 2 and res["A"].get("cur_wrong_family", 0) == 0
    assert res["A"]["graph_wrong_family"] == 1 and res["B"].get("graph_wrong_family", 0) == 0


def test_a_gain_counts_as_sequence_supported_only_when_the_genome_separates_the_copies():
    assert O.seq_supported({"copy0": 3, "copy1": 5}, "copy0")             # another copy's locus is 2 edits worse
    assert not O.seq_supported({"copy0": 3, "copy1": 4}, "copy0")         # within 1 edit: only the annotation separates them
    assert O.seq_supported({"copy0": 3}, "copy0")                         # no other locus explains the read at all
    assert not O.seq_supported({"copy1": 3}, "copy0")                     # its own locus does not explain it


def test_the_score_reports_sequence_supported_correct_calls(tmp_path):
    w = _world(tmp_path, [("r1", "A", "assigned", "0")], [("r1", "A", "0", "assigned")])
    res = O.score_world(w, "o2.t.tsv", locus={"A": {"r1": {"copy0": 2, "copy1": 3}}})
    assert res["A"]["cur_correct"] == 1 and res["A"].get("cur_correct_seq", 0) == 0
    res = O.score_world(w, "o2.t.tsv", locus={"A": {"r1": {"copy0": 2, "copy1": 9}}, "B": {"r1": {"copy0": 0}}})
    assert res["A"]["cur_correct_seq"] == 1 and res["A"]["graph_correct_seq"] == 1


def test_reads_of_a_paralog_annotated_structure_are_scored_on_their_own(tmp_path):
    w = _world(tmp_path, [("r1", "A", "assigned", "0"), ("r4", "A", "assigned", "0")],
               [("r1", "A", "0", "assigned"), ("r4", "A", "1", "assigned")])
    with open(tmp_path / "reads.truth.E0.tsv", "a") as fh:
        fh.write("r4\tA.0\tpresent\tA\tproj:A.1|t9\t0\n")              # copy A.0 expressing a structure annotated on A.1
    res = O.score_world(w, "o2.t.tsv")["A"]
    assert res["proj_reads"] == 1 and res["cur_correct_proj"] == 1
    assert res["graph_assigned_proj"] == 1 and res.get("graph_correct_proj", 0) == 0     # path O2 names the annotated paralog
    assert res["reads"] == 2 and res["graph_correct"] == 1


def test_union_paths_give_every_copy_each_paralog_structure_it_can_carry():
    copies = [dict(idx=0, isoforms={"a": [[0, 10], [20, 30]], "b": [[0, 10], [25, 30]]}),
              dict(idx=1, isoforms={"c": [[100, 110], [120, 130]]})]

    def projector(c, queries):            # copy 1 can carry a's structure; b does not project onto it; c projects onto 0
        out = {}
        for name in queries:
            if c["idx"] == 1:
                out[name] = [(100, 110), (120, 130)] if name.endswith("|a") else None
            else:
                out[name] = [(0, 10), (40, 50)]
        return out
    u = O.union_structures(copies, projector)
    assert sorted(u[0]) == ["a", "b", "proj:copy1|c"]          # c's intron (10, 40) is new to copy 0
    assert sorted(u[1]) == ["c"]                                # a projects onto copy 1 as copy 1's own chain: nothing new


def test_locus_o2_decides_on_genomic_locus_distances(tmp_path):
    import json
    w = _world(tmp_path, [("r1", "A", "assigned", "0"), ("r2", "B", "assigned", "0")], [])
    json.dump({"A": {"r1": {"copy0": 2, "copy1": 9}}, "B": {"r2": {"copy0": 3, "copy1": 4}}},
              open(tmp_path / "graph" / "locus.dist.json", "w"))
    (tmp_path / "reads.E0.fq").write_text("@r1\n" + "A" * 1000 + "\n+\n" + "I" * 1000 + "\n@r2\n" + "A" * 1000 + "\n+\n" + "I" * 1000 + "\n")
    O.run_locus(w, margin=2)
    rows = {(r.split("\t")[0], r.split("\t")[1]): r.split("\t")[2:] for r in (tmp_path / "graph" / "o2.locus.m2.tsv").read_text().splitlines()[1:]}
    assert rows[("r1", "A")] == ["0", "assigned"] and rows[("r2", "B")] == ["-", "thin"]


def _paf(q, qlen, qs, qe, strand, t, tlen, ts, te, nm):
    return "\t".join(map(str, [q, qlen, qs, qe, strand, t, tlen, ts, te, 0, 0, 60, f"NM:i:{nm}", "tp:A:P"]))


def test_read_bases_past_a_path_end_are_annotation_not_distance():
    paf = "\n".join([
        _paf("r", 1000, 0, 984, "+", "copy0|t", 984, 0, 984, 2),       # 16 bases past copy 0's path end: not counted
        _paf("r", 1000, 0, 1000, "+", "copy1|t", 1100, 50, 1050, 5),
        _paf("s", 1000, 0, 900, "+", "copy0|t", 1100, 0, 900, 1),       # stops 200 bp before the path end: a real difference
        _paf("u", 1000, 16, 1000, "-", "copy0|t", 984, 0, 984, 2)])     # minus strand: the read's start lies past the path end
    assert O._paf_dist(paf) == {"r": {"copy0": 18, "copy1": 5}, "s": {"copy0": 101}, "u": {"copy0": 18}}
    assert O._paf_dist(paf, ends=True) == {"r": {"copy0": 2, "copy1": 5}, "s": {"copy0": 101}, "u": {"copy0": 2}}


def test_padded_paths_extend_both_transcript_ends_with_the_copy_s_genome():
    assert O.padded([(100, 200), (300, 400)], 50) == [(50, 200), (300, 450)]
    assert O.padded([(30, 200)], 50) == [(0, 250)]                       # not before the contig start
    assert O.padded([(100, 200), (300, 400)], 0) == [(100, 200), (300, 400)]


def test_a_gain_is_psv_separated_when_the_runner_up_carries_the_read_s_structure():
    true = [(0, 100), (200, 300), (400, 500)]
    assert O.gain_kind(true, [[(0, 100), (205, 300), (400, 500)]]) == "psv"          # same chain within 10 bp
    assert O.gain_kind(true, [[(0, 100), (400, 500)], [(0, 100), (250, 300), (400, 500)]]) == "structure"
    assert O.gain_kind(true, []) == "structure"


def test_a_read_union_cannot_explain_falls_back_to_the_loci_unless_another_family_explains_it():
    loc = {"copy0": 2, "copy1": 9}
    assert O.fallback_call("1", "assigned", loc, 1000, other_explained=False) == ("1", "assigned")      # union decides
    assert O.fallback_call("-", "tied", loc, 1000, other_explained=False) == ("-", "tied")
    assert O.fallback_call("-", "unexplained", loc, 1000, other_explained=False) == ("0", "locus")     # structure-free
    assert O.fallback_call("-", "unexplained", loc, 1000, other_explained=True) == ("-", "other_family")
    assert O.fallback_call("-", "unexplained", {"copy0": 40, "copy1": 45}, 1000, other_explained=False) == ("-", "unexplained")
    assert O.fallback_call("-", "unexplained", {}, 1000, other_explained=False) == ("-", "unexplained")


def test_reads_of_a_cluster_with_linked_private_alleles_are_flagged_whatever_its_divergence():
    clusters = [dict(family="F", host="0", cluster="0", n_reads="25", status="unexplained", est_div="0.0012"),
                dict(family="F", host="0", cluster="1", n_reads="40", status="explained", est_div="0.0"),
                dict(family="F", host="1", cluster="0", n_reads="12", status="unexplained", est_div="0.02"),
                dict(family="G", host="0", cluster="0", n_reads="30", status="unresolved", est_div="0.001")]
    members = {("F", "0", "0"): {"a", "b"}, ("F", "0", "1"): {"c"}, ("F", "1", "0"): {"d"}, ("G", "0", "0"): {"e"}}
    assert O.flag_reads(clusters, members) == {("F", "a"), ("F", "b")}
    assert O.flag_reads(clusters, members, min_reads=10) == {("F", "a"), ("F", "b"), ("F", "d")}


def test_a_flagged_read_is_reported_as_a_candidate_not_assigned():
    rows = [dict(read="a", family="F", call="1", reason="assigned"), dict(read="b", family="F", call="-", reason="tied"),
            dict(read="c", family="F", call="0", reason="locus")]
    out = O.apply_flag(rows, {("F", "a"), ("F", "b")})
    assert [(r["call"], r["reason"]) for r in out] == [("-", "flagged"), ("-", "tied"), ("0", "locus")]


def test_a_new_seeded_copy_maps_to_the_truth_copy_whose_locus_it_overlaps():
    truth_copies = [dict(copy="F.0", chrom="c1", start="1000", end="3000"), dict(copy="F.5", chrom="c2", start="5000", end="6500"),
                    dict(copy="F.6", chrom="c2", start="9000", end="9500")]
    new = [dict(idx=7, contig="c2", start=5100, end=6300), dict(idx=8, contig="c3", start=0, end=900)]
    assert O.map_new_copies(new, truth_copies) == {7: "F.5", 8: None}


def test_calls_are_scored_against_each_read_s_true_copy_whatever_catalog_named_it():
    truth = {"r1": dict(copy="F.0", family="F", kind="present"), "r2": dict(copy="F.5", family="F", kind="outside"),
             "r3": dict(copy="F.9", family="F", kind="absent"), "r4": dict(copy="F.0", family="F", kind="present")}
    idx2copy = {("F", 0): "F.0", ("F", 7): "F.5", ("F", 8): None}
    calls = {("r1", "F"): 0, ("r2", "F"): 7, ("r3", "F"): 0, ("r4", "F"): 8}
    s = O.score_assign(calls, idx2copy, truth)
    assert s == dict(correct=2, wrong=2, outside_own=1)


def test_placements_of_an_undetectable_absent_copy_are_counted_apart_from_errors():
    truth = {"r1": dict(copy="F.0", family="F", kind="present"), "r2": dict(copy="F.9", family="F", kind="absent"),
             "r3": dict(copy="F.8", family="F", kind="absent")}
    idx2copy = {("F", 0): "F.0", ("F", 1): "F.1"}
    calls = {("r1", "F"): 0, ("r2", "F"): 1, ("r3", "F"): 1}
    # F.9 is within delta of the genome (undetectable): its reads cannot be told from that locus's own
    s = O.score_assign(calls, idx2copy, truth, unknowable={"F.9"})
    assert s == dict(correct=1, wrong=1, outside_own=0, unknowable_placed=1)
