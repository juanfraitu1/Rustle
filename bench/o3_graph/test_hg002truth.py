import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_hg002truth", os.path.join(HERE, "hg002truth.py"))
H = importlib.util.module_from_spec(spec)
spec.loader.exec_module(H)

IDX = {"c1": [(100, 300, "F.0"), (500, 700, "F.0"), (5100, 5300, "F.1"), (5500, 5700, "F.1")], "c2": [(100, 400, "G.0")]}


def test_a_read_is_tied_when_its_best_alignment_score_is_shared_by_two_copies():
    alns = [("c1", [(150, 300), (500, 650)], 900), ("c1", [(5150, 5300), (5500, 5650)], 900), ("c2", [(100, 300)], 400)]
    assert H.tied_copies(alns, IDX) == (900, {"F.0", "F.1"})
    alns[1] = ("c1", [(5150, 5300), (5500, 5650)], 899)                  # one point lower: a unique best
    assert H.tied_copies(alns, IDX) == (900, {"F.0"})
    assert H.tied_copies([("c1", [(310, 490)], 900)], IDX) == (900, set())   # on no copy's exons


def test_the_donor_s_own_genome_settles_a_tie_only_when_its_best_placements_lift_to_one_copy():
    lift = {("chr1_MATERNAL", 1000): ("c1", 150), ("chr1_PATERNAL", 2000): ("c1", 160),   # both alleles of copy F.0
            ("chr1_MATERNAL", 9000): ("c1", 5150), ("chr7_MATERNAL", 500): None}
    lf = lambda contig, s, e: [(lift[(contig, s)][0], lift[(contig, s)][1], lift[(contig, s)][1] + (e - s))] if lift.get((contig, s)) else []  # noqa: E731
    own = [("chr1_MATERNAL", [(1000, 1150)], 950), ("chr1_PATERNAL", [(2000, 2150)], 950), ("chr1_MATERNAL", [(9000, 9150)], 940)]
    assert H.own_truth(own, lf, IDX) == ("unique", "F.0")                  # one copy, both haplotypes
    own[2] = ("chr1_MATERNAL", [(9000, 9150)], 950)
    assert H.own_truth(own, lf, IDX) == ("own_tie", None)                  # the donor's copies tie too: undecidable
    assert H.own_truth([("chr7_MATERNAL", [(500, 650)], 950)], lf, IDX) == ("no_counterpart", None)   # a donor copy CHM13 lacks
    lift[("chr9_MATERNAL", 700)] = ("c2", 5000)                           # lifts onto CHM13, but onto no copy of any family
    assert H.own_truth([("chr9_MATERNAL", [(700, 850)], 950)], lf, IDX) == ("elsewhere", None)
    assert H.own_truth([], lf, IDX) == ("unmapped", None)


def test_the_tied_table_names_each_tied_read_s_copies_and_minimap2_s_coin_toss_primary(tmp_path):
    import csv
    import json
    import pysam
    cp = lambda gid, s: dict(gene_id=gid, chrom="c1", start=s, end=s + 600,  # noqa: E731
                             transcripts={"t": [[s, s + 200], [s + 400, s + 600]]})
    json.dump({"sets": {"S": [dict(family="F", copies=[cp("a", 100), cp("b", 5100)])]}}, open(tmp_path / "f.json", "w"))
    hdr = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": "c1", "LN": 10000}]})
    with pysam.AlignmentFile(str(tmp_path / "x.bam"), "wb", header=hdr) as out:
        for name, pos, flag, score, mapq in [("r1", 5150, 0, 280, 0), ("r1", 150, 256, 280, 0),     # tied; primary on b
                                             ("r2", 150, 0, 280, 7), ("r2", 5150, 256, 270, 0)]:    # a unique best
            a = pysam.AlignedSegment(hdr)
            a.query_name, a.reference_id, a.reference_start, a.flag, a.mapping_quality = name, 0, pos, flag, mapq
            a.cigarstring, a.query_sequence = "150M", "A" * 150
            a.set_tag("AS", score)
            out.write(a)
    assert H.main(["tied", str(tmp_path / "x.bam"), str(tmp_path / "f.json"), "S", str(tmp_path / "t.tsv")]) == 0
    rows = list(csv.DictReader(open(tmp_path / "t.tsv"), delimiter="\t"))
    assert [(r["read"], r["copies"], r["primary"], r["mapq"]) for r in rows] == [("r1", "F.0,F.1", "F.1", "0")]
    assert H.main(["untied", str(tmp_path / "x.bam"), str(tmp_path / "f.json"), "S", str(tmp_path / "u.tsv")]) == 0
    rows = list(csv.DictReader(open(tmp_path / "u.tsv"), delimiter="\t"))
    assert [(r["read"], r["copies"], r["primary"], r["mapq"]) for r in rows] == [("r2", "F.0", "F.0", "7")]


def test_a_tied_read_s_fate_is_judged_by_what_the_donor_s_own_genome_says():
    F, tied = H.tied_fate, {"F.0", "F.1"}
    assert F("unique", "F.0", "F.0", False, tied) == "correct"
    assert F("unique", "F.0", "F.1", False, tied) == "error"
    assert F("unique_elsewhere", "F.2", "F.0", False, tied) == "error"   # the donor's copy is one CHM13 did not tie it with
    assert F("unique", "F.0", None, False, tied) == "abstain"
    assert F("unique", "F.0", None, True, tied) == "error"               # an absent-copy call on a present copy's read
    assert F("own_tie", None, "F.1", False, tied) == "guess"             # undecidable even in the donor: placing it is a guess
    assert F("own_tie", None, "wrong-family", False, tied) == "error"    # ... but not off the copies it ties on
    assert F("own_tie", None, None, True, tied) == "error"               # the donor's tied copies are both in CHM13
    assert F("no_counterpart", None, "F.0", False, tied) == "error"      # a donor copy CHM13 lacks: no present copy is right
    assert F("no_counterpart", None, None, True, tied) == "correct"
    assert F("no_counterpart", None, None, False, tied) == "abstain"
    assert F("elsewhere", None, "F.0", False, tied) == "error"
    assert F("elsewhere", None, None, True, tied) == "o3"                # another CHM13 locus: reported apart


def test_the_score_tallies_each_method_s_fates_by_truth_class_over_the_world_s_tied_reads(tmp_path):
    import json
    w = tmp_path / "W"
    for d in ("cur", "graph"):
        (w / d).mkdir(parents=True)
    tsv = lambda p, rows: open(p, "w").write("\n".join("\t".join(r) for r in rows) + "\n")  # noqa: E731
    json.dump({"F": [{"idx": 0, "copy": "F.0"}, {"idx": 1, "copy": "F.1"}, {"idx": 2, "copy": "F.2"}]}, open(w / "isoforms.json", "w"))
    json.dump({"F": {"copies": [{"origin": "catalog", "idx": i} for i in range(3)]}}, open(w / "graph/seeded.json", "w"))
    tsv(w / "copies.tsv", [("copy", "family", "chrom", "start", "end")])
    tsv(w / "cur/A.assignments.tsv", [("read_name", "family_id", "status", "catalog_copy_idx"), ("r1", "F", "assigned", "0"),
                                      ("r2", "F", "ambiguous", "1")])
    tsv(w / "cur/A.psv_cluster_reads.tsv", [("read_name", "family", "host", "cluster"), ("r3", "F", "0", "0")])
    tsv(w / "cur/A.psv_clusters.tsv", [("family", "host", "cluster", "status", "est_div"), ("F", "0", "0", "unexplained", "0.05")])
    tsv(w / "graph/o2.integv.m2.tsv", [("read", "family", "call", "reason"), ("r1", "F", "0", "ok"), ("r2", "F", "2", "ok")])
    tsv(w / "graph/calls.tsv", [("arm", "members"), ("I", "r3")])
    open(w / "reads.E0.fq", "w").write("".join(f"@{q}\nACGT\n+\nIIII\n" for q in ("r1", "r2", "r3", "r4")))
    tsv(tmp_path / "tied.tsv", [("read", "as", "family", "copies", "primary", "mapq"), ("r1", "900", "F", "F.0,F.1", "F.1", "0"),
                                ("r2", "900", "F", "F.0,F.1", "F.0", "3"), ("r3", "900", "F", "F.0,F.1", "F.1", "0"),
                                ("r4", "900", "F", "F.0,F.1", "F.0", "0"), ("r5", "900", "F", "F.0,F.1", "F.0", "0")])
    tsv(tmp_path / "truth.tsv", [("read", "family", "tied_copies", "class", "truth"), ("r1", "F", "F.0,F.1", "unique", "F.0"),
                                 ("r2", "F", "F.0,F.1", "own_tie", ""), ("r3", "F", "F.0,F.1", "no_counterpart", ""),
                                 ("r4", "F", "F.0,F.1", "unmapped", ""), ("r5", "F", "F.0,F.1", "unique", "F.1")])
    out = tmp_path / "score.json"
    assert H.main(["score", str(tmp_path / "truth.tsv"), str(tmp_path / "tied.tsv"), str(w), str(out)]) == 0
    got = json.load(open(out))
    assert got["e2e"] == {"mm2": {"unique": {"error": 1}, "own_tie": {"guess": 1}, "no_counterpart": {"error": 1}},
                          "mm2_q1": {"unique": {"abstain": 1}, "own_tie": {"guess": 1}, "no_counterpart": {"abstain": 1}},
                          "current": {"unique": {"correct": 1}, "own_tie": {"abstain": 1}, "no_counterpart": {"correct": 1}},
                          "integrated": {"unique": {"correct": 1}, "own_tie": {"error": 1}, "no_counterpart": {"correct": 1}},
                          "current_tied": {"unique": {"correct": 1}, "own_tie": {"abstain": 1}, "no_counterpart": {"correct": 1}},
                          "integrated_tied": {"unique": {"correct": 1}, "own_tie": {"abstain": 1}, "no_counterpart": {"correct": 1}}}
    assert got["o2"]["integrated"] == {"unique": {"correct": 1}, "own_tie": {"error": 1}, "no_counterpart": {"abstain": 1}}
    assert got["o2"]["integrated_tied"]["own_tie"] == {"abstain": 1}     # F.2 is no copy r2 ties on: scoped out
    assert got["o3"]["current"] == {"no_counterpart": 1} and got["o3"]["mm2"] == {}
    assert got["skipped"] == {"unmapped": 1, "not_in_world": 1}


def test_the_donor_s_best_placement_must_match_near_perfectly_for_a_truth():
    assert H.near_perfect(2, 500) and not H.near_perfect(3, 500)          # at least 2 edits allowed
    assert H.near_perfect(6, 2000) and not H.near_perfect(7, 2000)        # 0.3% of the aligned bases


def test_the_donor_margin_is_how_far_the_best_placement_off_the_truth_copy_falls_behind():
    lift = {("chr1_MATERNAL", 1000): ("c1", 150), ("chr1_PATERNAL", 2000): ("c1", 160), ("chr1_MATERNAL", 9000): ("c1", 5150)}
    lf = lambda contig, s, e: [(lift[(contig, s)][0], lift[(contig, s)][1], lift[(contig, s)][1] + (e - s))] if lift.get((contig, s)) else []  # noqa: E731
    own = [("chr1_MATERNAL", [(1000, 1150)], 950, 1, 150), ("chr1_PATERNAL", [(2000, 2150)], 948, 2, 150),
           ("chr1_MATERNAL", [(9000, 9150)], 940, 3, 150), ("chr7_MATERNAL", [(500, 650)], 930, 4, 150)]
    assert H.own_truth(own, lf, IDX) == ("unique", "F.0")                  # five-field placements (NM, aligned bases) too
    assert H.donor_margin(own, lf, IDX, "F.0") == 10                       # F.1's 940; the paternal allele of F.0 is no rival
    assert H.donor_margin(own[:2], lf, IDX, "F.0") is None                 # nothing else placed it


def test_the_truth_table_classes_each_tied_read_by_the_donor_s_haplotypes_lifted_through_the_chain(tmp_path):
    import csv
    import json
    import pysam
    cp = lambda gid, s: dict(gene_id=gid, chrom="c1", start=s, end=s + 600,  # noqa: E731
                             transcripts={"t": [[s, s + 200], [s + 400, s + 600]]})
    json.dump({"sets": {"S": [dict(family="F", copies=[cp("a", 100), cp("b", 5100)])]}}, open(tmp_path / "f.json", "w"))
    open(tmp_path / "x.chain", "w").write("chain 1 h1_MATERNAL 10000 + 0 10000 c1 10000 + 0 10000 1\n10000\n\n"
                                          "chain 1 h1_PATERNAL 10000 + 0 10000 c1 10000 + 0 10000 2\n10000\n\n")

    def bam(path, contigs, recs):
        hdr = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": c, "LN": 10000} for c in contigs]})
        with pysam.AlignmentFile(str(path), "wb", header=hdr) as out:
            for i, (name, contig, pos, score, nm) in enumerate(recs):
                a = pysam.AlignedSegment(hdr)
                a.query_name, a.reference_id, a.reference_start = name, contigs.index(contig), pos
                a.flag, a.mapping_quality, a.cigarstring, a.query_sequence = 0 if i == 0 or recs[i - 1][0] != name else 256, 0, "150M", "A" * 150
                a.set_tag("AS", score)
                a.set_tag("NM", nm)
                out.write(a)
    bam(tmp_path / "mat.bam", ["h1_MATERNAL", "h2_MATERNAL"],
        [("r1", "h1_MATERNAL", 150, 300, 0), ("r1", "h1_MATERNAL", 5150, 290, 2),     # unique on F.0, margin 10
         ("r2", "h1_MATERNAL", 150, 300, 0),                                         # ... and on F.1 paternally: own_tie
         ("r3", "h2_MATERNAL", 150, 300, 0),                                         # no chain: a donor copy CHM13 lacks
         ("r4", "h1_MATERNAL", 150, 300, 9),                                         # 9 edits in 150 bases: poor_match
         ("r5", "h1_MATERNAL", 5150, 300, 0)])                                       # F.1, which CHM13 did not tie it on
    bam(tmp_path / "pat.bam", ["h1_PATERNAL"], [("r2", "h1_PATERNAL", 5150, 300, 0)])
    open(tmp_path / "tied.tsv", "w").write("read\tas\tfamily\tcopies\tprimary\tmapq\n" + "".join(
        f"r{i}\t900\tF\t{'F.0,F.2' if i == 5 else 'F.0,F.1'}\tF.0\t0\n" for i in range(1, 7)))
    assert H.main(["truth", str(tmp_path / "tied.tsv"), str(tmp_path / "mat.bam"), str(tmp_path / "pat.bam"),
                   str(tmp_path / "x.chain"), str(tmp_path / "f.json"), "S", str(tmp_path / "truth.tsv")]) == 0
    rows = [(r["read"], r["class"], r["truth"], r["donor_nm"], r["donor_margin"])
            for r in csv.DictReader(open(tmp_path / "truth.tsv"), delimiter="\t")]
    assert rows == [("r1", "unique", "F.0", "0", "10"), ("r2", "own_tie", "", "0", ""), ("r3", "no_counterpart", "", "0", ""),
                    ("r4", "poor_match", "", "9", ""), ("r5", "unique_elsewhere", "F.1", "0", ""), ("r6", "unmapped", "", "", "")]


def test_the_world_holds_the_families_with_scored_tied_reads_and_its_truth_the_reads_with_a_donor_copy():
    fams = [dict(family="F", copies=[]), dict(family="G", copies=[]), dict(family="K", copies=[])]
    rows = [dict(read="r1", family="F", **{"class": "unique"}, truth="F.0"),
            dict(read="r2", family="F", **{"class": "own_tie"}, truth=""),
            dict(read="r3", family="G", **{"class": "poor_match"}, truth=""),        # unscored: G has no scored read
            dict(read="r4", family="K", **{"class": "unique_elsewhere"}, truth="K.2"),
            dict(read="r5", family="K", **{"class": "no_counterpart"}, truth="")]
    assert [f["family"] for f in H.select_families(fams, rows, 1)] == ["F", "K"]
    assert [f["family"] for f in H.select_families(fams, rows, 2)] == ["F", "K"]
    assert H.select_families(fams, rows, 3) == []
    assert H.world_truth(rows) == [("r1", "F.0", "F"), ("r4", "K.2", "K")]


def _cards(placed_wrong, placed_right, own_tie_placed, decidable_right, mm2_wrong, o3_present=0, n_dec=200, n_tie=200):
    """An o2proof pooled card and an hg002truth score for one system ('integrated') and the coin toss."""
    def cls(correct, wrong, abstain):
        return {"all": dict(correct=correct, wrong=wrong, abstain=abstain)}
    o2p = {"pooled": {"card": {
        "integrated": {"decidable": cls(decidable_right, 0, n_dec - decidable_right),
                       "undecidable": cls(placed_right - decidable_right, placed_wrong, 300)},
        "mm2": {"decidable": cls(n_dec, 0, 0), "undecidable": cls(300 - mm2_wrong, mm2_wrong, 0)}}}}
    score = {"o2": {"integrated": {"own_tie": {"guess": own_tie_placed, "abstain": n_tie - own_tie_placed}}},
             "o3": {"integrated": {"unique": o3_present}},
             "e2e": {"integrated": {"unique": {"correct": 1000}, "own_tie": {"abstain": 200}}}}
    return o2p, score


def test_the_bars_pass_a_precise_calibrated_powerful_system_and_fail_or_withhold_otherwise():
    o2p, sc = _cards(placed_wrong=2, placed_right=300, own_tie_placed=4, decidable_right=190, mm2_wrong=150)
    b = H.bars(o2p, sc, systems=("integrated",))["integrated"]
    assert b["T1"]["holds"] and b["T2"]["holds"] and b["T3"]["holds"] and b["T4"]["holds"] and b["verdict"] == "PASS"
    o2p, sc = _cards(placed_wrong=20, placed_right=300, own_tie_placed=4, decidable_right=190, mm2_wrong=150)
    assert H.bars(o2p, sc, systems=("integrated",))["integrated"]["verdict"] == "FAIL"           # 20/320 wrong: T1
    o2p, sc = _cards(placed_wrong=0, placed_right=300, own_tie_placed=40, decidable_right=190, mm2_wrong=150)
    assert H.bars(o2p, sc, systems=("integrated",))["integrated"]["verdict"] == "FAIL"           # 40/200 guesses: T3
    o2p, sc = _cards(placed_wrong=0, placed_right=60, own_tie_placed=0, decidable_right=50, mm2_wrong=150, n_dec=50)
    b = H.bars(o2p, sc, systems=("integrated",))["integrated"]
    assert b["T1"]["powered"] is False and b["verdict"] == "UNDERPOWERED"                         # 60 placed < 100
    o2p, sc = _cards(placed_wrong=0, placed_right=300, own_tie_placed=0, decidable_right=190, mm2_wrong=10)
    b = H.bars(o2p, sc, systems=("integrated",))
    assert b["valid"] is False and b["integrated"]["T2"]["holds"] is None and b["integrated"]["verdict"] == "UNDERPOWERED"


def test_the_o3_safety_bar_counts_absent_copy_calls_on_reads_whose_copy_is_present():
    o2p, sc = _cards(placed_wrong=0, placed_right=300, own_tie_placed=0, decidable_right=190, mm2_wrong=150, o3_present=80)
    assert H.bars(o2p, sc, systems=("integrated",))["integrated"]["T5"]["holds"] is False         # 80 of 1,200 present
    assert H.wilson(0, 100)[0] == 0.0 and 0.03 < H.wilson(0, 100)[1] < 0.04


def test_minimap2_resolves_a_read_only_when_its_best_as_alignments_all_lie_on_one_copy():
    one = [("c1", [(150, 300)], 900), ("c1", [(5150, 5300)], 880)]
    assert H.resolved_copy(one, IDX) == "F.0"
    assert H.resolved_copy(one + [("c1", [(160, 300)], 900)], IDX) == "F.0"              # two best alignments, one copy
    assert H.resolved_copy(one + [("c1", [(5150, 5300)], 900)], IDX) is None             # an AS tie between copies: O2's
    assert H.resolved_copy(one + [("c9", [(150, 300)], 900)], IDX) is None               # an AS tie with a locus off the copies
    assert H.resolved_copy([("c9", [(150, 300)], 900)], IDX) is None and H.resolved_copy([], IDX) is None


def test_the_scoped_o2_keeps_a_call_only_on_one_of_the_copies_the_read_ties_on():
    assert H.scoped("F.1", {"F.0", "F.1"}) == "F.1"
    assert H.scoped("F.2", {"F.0", "F.1"}) is None and H.scoped("wrong-family", {"F.0", "F.1"}) is None
    assert H.scoped(None, {"F.0", "F.1"}) is None


def test_the_premise_bar_is_minimap2_s_agreement_with_the_donor_on_reads_it_resolves():
    rows = [{"class": "unique"}] * 995 + [{"class": "unique_elsewhere"}] * 5 + [{"class": "own_tie"}] * 50
    b = H.premise_bar(rows)
    assert b["value"] == 0.995 and b["n"] == 1000 and b["powered"] and b["holds"]
    assert H.premise_bar(rows[:990] + [{"class": "unique_elsewhere"}] * 20)["holds"] is False


def test_a_family_is_pseudoautosomal_when_every_copy_lies_in_a_par():
    cp = lambda ch, s, e: dict(chrom=ch, start=s, end=e)  # noqa: E731
    assert H.pseudoautosomal(dict(family="S", copies=[cp("chrX", 1229797, 1235760), cp("chrY", 1249920, 1256001)]))
    assert H.pseudoautosomal(dict(family="V", copies=[cp("chrX", 154109336, 154171772), cp("chrY", 62306324, 62368760)]))
    assert not H.pseudoautosomal(dict(family="R", copies=[cp("chrY", 2518432, 2609182), cp("chrY", 21631329, 21656255)]))
    assert not H.pseudoautosomal(dict(family="N", copies=[cp("chrX", 5392902, 5781161), cp("chrY", 15429200, 15752213)]))
    fams = [dict(family="S", copies=[cp("chrX", 1229797, 1235760), cp("chrY", 1249920, 1256001)]),
            dict(family="E", copies=[cp("chr16", 28000000, 28010000), cp("chr16", 28700000, 28710000)])]
    rows = [dict(read="a", family="S", **{"class": "own_tie"}, truth=""), dict(read="b", family="E", **{"class": "unique"}, truth="E.0")]
    assert [f["family"] for f in H.select_families(fams, rows, 1, par="exclude")] == ["E"]
    assert [f["family"] for f in H.select_families(fams, rows, 1, par="only")] == ["S"]
    assert [f["family"] for f in H.select_families(fams, rows, 1)] == ["S", "E"]


def test_the_premise_command_writes_the_bar_from_the_untied_truth_table(tmp_path):
    import json
    open(tmp_path / "t.tsv", "w").write("read\tfamily\ttied_copies\tclass\ttruth\tdonor_nm\tdonor_margin\n"
                                        + "".join(f"r{i}\tF\tF.0\tunique\tF.0\t0\t\n" for i in range(150))
                                        + "q\tF\tF.0\tunique_elsewhere\tF.1\t0\t\n")
    assert H.main(["premise", str(tmp_path / "t.tsv"), str(tmp_path / "p.json")]) == 0
    got = json.load(open(tmp_path / "p.json"))
    assert got["n"] == 151 and got["holds"] and got["powered"] and got["classes"] == {"unique": 150, "unique_elsewhere": 1}
