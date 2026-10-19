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
