import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import junction_truth as J


def test_read_junctions_forward_record():
    # 5S 10M 100N 20M : junction after 5 + 10 query bases
    assert J.read_junctions("5S10M100N20M", 0, 35) == [15]


def test_read_junctions_reverse_record_are_in_original_orientation():
    # same record on the reverse strand: original read = revcomp of SEQ, so position q becomes qlen - q
    assert J.read_junctions("5S10M100N20M", 16, 35) == [20]


def test_short_gaps_are_not_junctions():
    assert J.read_junctions("10M10N20M", 0, 30) == []


def test_paf_introns_plus_strand():
    # query 0-30 on target from 100: 10M 50N 20M -> intron target [110, 160), read junction 10
    assert J.paf_introns("10M50N20M", "+", 0, 30, 100) == [(10, 110, 160)]


def test_paf_introns_minus_strand_walk_query_backwards():
    # '-' hit, query interval [0, 30): the cigar walks the target forward and the query from qe down
    # 20M (query 30->10) 50N 10M (query 10->0): junction between query 9 and 10 -> 10
    assert J.paf_introns("20M50N10M", "-", 0, 30, 100) == [(10, 120, 170)]


def test_motif_in_transcript_orientation():
    seq = "AAGTAAAAAAAGCC"  # intron = [2, 12): GT ... AG on the forward strand
    assert J.motif(seq, 2, 12, "+") == ("GT", "AG")
    # on the reverse strand the transcript-oriented intron is revcomp(seq[2:12]) = CT...AC -> ("CT", "AC")
    assert J.motif(seq, 2, 12, "-") == ("CT", "AC")
    assert J.canonical(("GT", "AG")) and J.canonical(("GC", "AG")) and J.canonical(("AT", "AC"))
    assert not J.canonical(("CT", "AC")) and not J.canonical(("GA", "AG"))


def test_classify_excluded_produces_undetermined():
    juncs = [100, 300]
    good = [(101, 0, 0, True), (299, 0, 0, True)]          # (read pos, t0, t1, canonical)
    broken = [(100, 0, 0, True), (303, 0, 0, False)]
    partial = [(100, 0, 0, True)]                           # no intron near 300
    assert J.classify(juncs, good) == "PRODUCES"
    assert J.classify(juncs, broken) == "EXCLUDED"
    assert J.classify(juncs, partial) == "UNDETERMINED"
    assert J.classify(juncs, [(320, 0, 0, False)]) == "UNDETERMINED", "a broken intron 20 bp away is not this junction"


def test_one_locus_means_every_pair_overlaps():
    assert J.one_locus([("c", 0, 10), ("c", 5, 20)])
    assert not J.one_locus([("c", 0, 10), ("c", 15, 20)])
    assert not J.one_locus([("c", 0, 10), ("d", 0, 10)])
    assert not J.one_locus([])
