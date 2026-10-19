import json
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


def test_round8_bars_hold_per_world_only_when_integrated_has_no_fewer_correct_and_no_more_errors():
    tab = {"INTACT": {"current": {"catalogued": {"correct": 90, "error": 5, "abstain": 5}},
                      "integrated": {"catalogued": {"correct": 90, "error": 5, "abstain": 5}}},
           "REMOVED": {"current": {"catalogued": {"correct": 80, "error": 2},
                                   "absent_undetectable": {"unknowable": 9, "abstain": 1}},
                       "integrated": {"catalogued": {"correct": 79, "error": 0},
                                      "absent_undetectable": {"unknowable": 3, "abstain": 7}}},
           "OUTSIDE": {"current": {"uncatalogued": {"correct": 1, "error": 10}},
                       "integrated": {"uncatalogued": {"correct": 9, "error": 11}}}}
    b = S.round8_bars(tab)
    assert b["INTACT"] == dict(current_correct=90, integrated_correct=90, current_error=5, integrated_error=5, B1=True, B2=True)
    assert b["REMOVED"]["B1"] is False and b["REMOVED"]["B2"] is True       # unknowable and abstain count for neither
    assert b["OUTSIDE"]["B1"] is True and b["OUTSIDE"]["B2"] is False
    assert b["verdict"] == "FAIL" and S.round8_bars({"INTACT": tab["INTACT"]})["verdict"] == "PASS"


def test_a_call_on_a_seeded_copy_that_is_no_truth_copy_is_a_placement_not_an_abstention():
    # round 8 review I3: O1 seeded OR4F17 (no truth copy); O2's calls there were scored as abstentions
    rows = [dict(read="a", family="F", call="5"), dict(read="b", family="F", call="0"), dict(read="c", family="F", call="-")]
    calls = S.integrated_calls(rows, {("F", 0): "F.0", ("F", 5): None})
    assert calls["b"] == {"F": "F.0"} and "c" not in calls
    placed = calls["a"]["F"]
    assert S.read_fate("catalogued", placed, "F.1", False) == "error"
    assert S.read_fate("absent_undetectable", placed, "F.9", False) == "unknowable"


def test_an_absent_copy_s_read_is_detectable_by_its_own_divergence_when_asked():
    # copy-level (rounds 7-8): the copy's median de decides every read; per read (round 9): the read's own primary de
    assert S.absent_class(True, 0.002, per_read=False) == "absent_detectable"
    assert S.absent_class(True, 0.002, per_read=True) == "absent_undetectable"
    assert S.absent_class(False, 0.05, per_read=True) == "absent_detectable"
    assert S.absent_class(False, 0.05, per_read=False) == "absent_undetectable"
    assert S.absent_class(True, None, per_read=True) == "absent_detectable"      # no primary alignment: the genome explains nothing


def test_e2e_scores_the_integrated_o2_table_it_is_given():
    assert S.e2e_args(["e2e", "o.json", "W1", "W2"]) == ("o.json", ["W1", "W2"], False, "o2.integ.m2.tsv", False, None, False)
    assert S.e2e_args(["e2e", "o.json", "W1", "--o2", "o2.integl.m2.tsv", "W2", "--per-read-detectability"]) == (
        "o.json", ["W1", "W2"], True, "o2.integl.m2.tsv", False, None, False)


def test_per_read_detectability_subtracts_the_world_s_error_floor_when_given():
    # delta (0.00958) is the p99 exonic divergence between KB3781's haplotype ASSEMBLIES: an error-free divergence; a read's
    # `de` adds its sequencing errors (round 11)
    assert S.absent_class(True, 0.0118, per_read=True) == "absent_detectable"
    assert S.absent_class(True, 0.0118, per_read=True, floor=0.0029) == "absent_undetectable"
    assert S.absent_class(True, 0.0211, per_read=True, floor=0.0029) == "absent_detectable"
    assert S.absent_class(True, None, per_read=True, floor=0.0029) == "absent_detectable"     # no primary: unexplained


def test_e2e_takes_the_error_floor_option():
    assert S.e2e_args(["e2e", "o.json", "W1", "--error-floor"]) == ("o.json", ["W1"], False, "o2.integ.m2.tsv", True, None, False)
    assert S.e2e_args(["e2e", "o.json", "W1"])[4] is False


def test_a_world_s_calls_are_each_system_s_placements_and_the_reads_inside_its_absent_copy_calls(tmp_path):
    w = tmp_path / "W"
    for d in ("cur", "graph"):
        (w / d).mkdir(parents=True)
    tsv = lambda p, rows: open(w / p, "w").write("\n".join("\t".join(r) for r in rows) + "\n")  # noqa: E731
    json.dump({"F": [{"idx": 0, "copy": "F.0"}, {"idx": 1, "copy": "F.1"}]}, open(w / "isoforms.json", "w"))
    json.dump({"F": {"copies": [{"origin": "catalog", "idx": 0}, {"origin": "new", "idx": 2, "contig": "c1", "start": 900,
                                                                  "end": 1500}]}}, open(w / "graph/seeded.json", "w"))
    tsv("copies.tsv", [("copy", "family", "chrom", "start", "end"), ("F.9", "F", "c1", "1000", "1400")])
    tsv("cur/A.assignments.tsv", [("read_name", "family_id", "status", "catalog_copy_idx"), ("r1", "F", "assigned", "1"),
                                  ("r2", "F", "ambiguous", "0")])
    tsv("cur/A.psv_cluster_reads.tsv", [("read_name", "family", "host", "cluster"), ("r3", "F", "0", "0"), ("r7", "F", "0", "1")])
    tsv("cur/A.psv_clusters.tsv", [("family", "host", "cluster", "status", "est_div"), ("F", "0", "0", "unexplained", "0.05"),
                                   ("F", "0", "1", "unexplained", "0.001")])        # within DELTA: an allele, no call
    tsv("graph/o2.integv.m2.tsv", [("read", "family", "call", "reason"), ("r1", "F", "0", "ok"), ("r2", "F", "-", "tied"),
                                   ("r6", "F", "2", "ok")])
    tsv("graph/calls.tsv", [("arm", "members"), ("I", "r4,r5"), ("II", "r8")])
    got = S.world_calls(str(w), "o2.integv.m2.tsv")
    assert got == dict(current={"r1": {"F": "F.1"}}, integrated={"r1": {"F": "F.0"}, "r6": {"F": "F.9"}},
                       current_o3={"r3"}, integrated_o3={"r4", "r5"})


def test_a_read_tied_within_the_slack_with_a_locus_on_no_copy_of_the_world_is_scored_apart(tmp_path):
    import pysam
    w = tmp_path / "W"
    w.mkdir()
    open(w / "copies.tsv", "w").write("copy\tfamily\tgene_id\tgene_name\tchrom\tstart\tend\tstrand\tkind\tin_reference\tcatalogued\n"
                                      "F.0\tF\tg0\tg0\tc1\t1000\t2000\t+\tpresent\t1\t1\n"
                                      "G.0\tG\tg1\tg1\tc1\t5000\t6000\t+\tpresent\t1\t1\n"
                                      "F.1\tF\tg2\tg2\tc1\t9000\t9900\t+\tabsent\t0\t1\n")
    hdr = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": "c1", "LN": 20000}]})
    with pysam.AlignmentFile(str(w / "E0.bam"), "wb", header=hdr) as out:
        for name, pos, flag, score in [("r1", 1100, 0, 300), ("r1", 15000, 256, 295),     # within 9 of a never-expressed locus
                                       ("r2", 1100, 0, 300), ("r2", 15000, 256, 290),     # 10 below: resolved
                                       ("r3", 1100, 0, 300), ("r3", 5100, 256, 300),      # tied with another family's copy
                                       ("r4", 9100, 0, 300)]:                             # the masked copy's locus: no copy there
            a = pysam.AlignedSegment(hdr)
            a.query_name, a.reference_id, a.reference_start, a.flag, a.mapping_quality = name, 0, pos, flag, 0
            a.cigarstring, a.query_sequence = "150=", "A" * 150
            a.set_tag("AS", score)
            out.write(a)
    assert S.off_family_ties(str(w), 9) == {"r1", "r4"}
    assert S.off_family_ties(str(w), 0) == {"r4"}


def test_the_bars_leave_the_off_family_ties_out_and_e2e_takes_their_slack():
    tab = {"W": {"current": {"catalogued": {"correct": 10, "error": 0}, "off_family_tie": {"correct": 5}},
                 "integrated": {"catalogued": {"correct": 10, "error": 0}, "off_family_tie": {"abstain": 5}}}}
    b = S.round8_bars(tab)["W"]
    assert b["current_correct"] == 10 and b["B1"] and b["B2"]
    assert S.e2e_args(["e2e", "o.json", "W1", "--o2", "x.tsv", "--off-family-ties", "9"]) == ("o.json", ["W1"], False, "x.tsv", False, 9, False)
    assert S.e2e_args(["e2e", "o.json", "W1"])[-2] is None


def test_off_family_ties_take_the_scope_s_sequence_ties_when_asked(tmp_path):
    import pysam
    w = tmp_path / "W"
    w.mkdir()
    open(w / "copies.tsv", "w").write("copy\tfamily\tgene_id\tgene_name\tchrom\tstart\tend\tstrand\tkind\tin_reference\tcatalogued\n"
                                      "F.0\tF\tg0\tg0\tc1\t1000\t2000\t+\tpresent\t1\t1\n")
    hdr = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": "c1", "LN": 20000}]})
    with pysam.AlignmentFile(str(w / "E0.bam"), "wb", header=hdr) as out:
        for name, pos, flag, score, cig in [("r1", 1100, 0, 300, "150="), ("r1", 15000, 256, 260, "150=")]:   # off locus, AS -40
            a = pysam.AlignedSegment(hdr)
            a.query_name, a.reference_id, a.reference_start, a.flag, a.mapping_quality = name, 0, pos, flag, 0
            a.cigarstring, a.query_sequence = cig, "A" * 150
            a.set_tag("AS", score)
            out.write(a)
    assert S.off_family_ties(str(w), 9) == set() and S.off_family_ties(str(w), 9, seq_tie=True) == {"r1"}
    assert S.e2e_args(["e2e", "o.json", "W1", "--off-family-ties", "9", "--seq-tie"])[-2:] == (9, True)
    assert S.e2e_args(["e2e", "o.json", "W1", "--off-family-ties", "9"])[-2:] == (9, False)


def test_set_fates_count_placements_against_the_truth_and_whether_each_tied_read_s_set_holds_its_copy():
    rows = [dict(read="a", family="F", call="0", reason="fewest_edits", base_call="-", candidates="0,1", set="0"),
            dict(read="b", family="F", call="-", reason="fewest_edit_set", base_call="-", candidates="0,1,2", set="0,1"),
            dict(read="c", family="F", call="-", reason="fewest_edit_set", base_call="-", candidates="0,dropped", set="dropped"),
            dict(read="d", family="F", call="1", reason="as_resolved", base_call="1", candidates="-", set="-"),
            dict(read="e", family="F", call="1", reason="assigned", base_call="1", candidates="0,1", set="1")]
    truth = {"a": "F.0", "b": "F.1", "c": "F.3", "d": "F.1", "e": "F.0"}
    seeded = {("F", 0): "F.0", ("F", 1): "F.1", ("F", 2): "F.2"}
    primary = {"a": "F.1", "b": "F.1", "c": "F.0", "e": "F.1"}
    got = S.set_fates(rows, truth, seeded, {"F": {"F.3"}}, primary)
    assert got == {"ambiguous": 4, "placed right": 1, "placed wrong": 1, "abstain": 2, "v17 placed right": 0,
                   "v17 placed wrong": 1, "v17 abstain": 3, "set holds truth": 3, "candidates": 9, "set size": 5,
                   "minimap2 right": 1, "minimap2 wrong": 3}
    assert S.set_fates(rows, truth, seeded, {"F": {"F.3"}}, primary, skip={"b"})["ambiguous"] == 3   # off-family ties apart


def test_an_alignment_lies_on_the_copy_of_its_strand_whose_exons_it_overlaps_most():
    cands = [("F.0", "+", [(100, 200), (300, 1000)], (100, 1000)),     # CYB5B: F.1 is one exon inside F.0's last exon
             ("F.1", "+", [(600, 1000)], (600, 1000)),
             ("F.2", "-", [(900, 1500)], (900, 1500))]                  # RPL38: an antisense neighbour over F.0's end
    assert S.attribute([(150, 200), (300, 700)], "+", cands) == "F.0"           # spliced: F.0's exons hold more of it
    assert S.attribute([(650, 990)], "+", cands) == "F.0"                       # inside both equally: the first
    assert S.attribute([(950, 1400)], "-", cands) == "F.2"                      # the antisense copy by its strand
    assert S.attribute([(5000, 5100)], "+", cands) is None


def test_a_read_with_one_tied_alignment_is_no_multimapper():
    assert S.tied_records([(300, 0), (280, 0)], slack=9) == 1                   # the rival 20 AS down: one placement
    assert S.tied_records([(300, 0), (280, 0)], slack=9, seq_tie=True) == 2     # ... unless it has no more edits
    assert S.tied_records([(300, 1), (295, 0)], slack=9) == 2
