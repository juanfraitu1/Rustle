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
    json.dump({"F": [{"idx": 0, "copy": "F.0"}, {"idx": 1, "copy": "F.1"}]}, open(w / "isoforms.json", "w"))
    json.dump({"F": {"copies": [{"origin": "catalog", "idx": 0}, {"origin": "catalog", "idx": 1}]}}, open(w / "graph/seeded.json", "w"))
    tsv(w / "copies.tsv", [("copy", "family", "chrom", "start", "end")])
    tsv(w / "cur/A.assignments.tsv", [("read_name", "family_id", "status", "catalog_copy_idx"), ("r1", "F", "assigned", "0"),
                                      ("r2", "F", "ambiguous", "1")])
    tsv(w / "cur/A.psv_cluster_reads.tsv", [("read_name", "family", "host", "cluster"), ("r3", "F", "0", "0")])
    tsv(w / "cur/A.psv_clusters.tsv", [("family", "host", "cluster", "status", "est_div"), ("F", "0", "0", "unexplained", "0.05")])
    tsv(w / "graph/o2.integv.m2.tsv", [("read", "family", "call", "reason"), ("r1", "F", "0", "ok"), ("r2", "F", "-", "tied")])
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
    assert got["card"] == {"mm2": {"unique": {"error": 1}, "own_tie": {"guess": 1}, "no_counterpart": {"error": 1}},
                           "mm2_q1": {"unique": {"abstain": 1}, "own_tie": {"guess": 1}, "no_counterpart": {"abstain": 1}},
                           "current": {"unique": {"correct": 1}, "own_tie": {"abstain": 1}, "no_counterpart": {"correct": 1}},
                           "integrated": {"unique": {"correct": 1}, "own_tie": {"abstain": 1}, "no_counterpart": {"correct": 1}}}
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
