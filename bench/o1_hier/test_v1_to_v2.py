import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import v1_to_v2 as L  # noqa: E402

# v1 chr16 1000 bp, v2 995 bp: 500 aligned, 5 bp only in v1, 495 aligned
CHAIN = """chain 995 chr16 1000 + 0 1000 chr16 995 + 0 995 1
500 5 0
495

chain 1000 chr1 1000 + 0 1000 chr1 1000 + 0 1000 2
1000
"""


def lifter():
    return L.Lifter(L.read_chain(CHAIN.splitlines()))


def test_chain_blocks():
    assert L.read_chain(CHAIN.splitlines()) == {"chr16": [(0, 0, 500), (505, 500, 495)], "chr1": [(0, 0, 1000)]}


def test_identity_chromosome_lifts_unchanged():
    assert lifter().lift("chr1", 100, 200) == (100, 200, False)


def test_shift_after_a_deletion():
    assert lifter().lift("chr16", 600, 700) == (595, 695, False)


def test_an_interval_spanning_the_indel_lifts_with_it_and_is_flagged():
    assert lifter().lift("chr16", 400, 600) == (400, 595, True)


def test_an_end_inside_a_gap_is_unliftable():
    assert lifter().lift("chr16", 502, 510) is None
    assert lifter().lift("chr16", 490, 503) is None


def test_unknown_chromosome_is_unliftable():
    assert lifter().lift("chrY", 1, 5) is None


def test_rows_pairs_and_one_based():
    lines = ["chr16\t601\t700\tID_1\n", "chrY\t1\t5\tID_2\n"]
    out, c, _ = L.lift_rows(lifter(), lines, one_based=True)
    assert out == ["chr16\t596\t695\tID_1\n"] and c["rows_dropped"] == 1 and c["unliftable|chrY"] == 1
    out, c, _ = L.lift_rows(lifter(), ["chr1\t10\t20\tchr16\t600\t700\tS\n", "chr1\t10\t20\tchr16\t502\t510\tS\n"],
                            pairs=True)
    assert out == ["chr1\t10\t20\tchr16\t595\t695\tS\n"] and c["rows_dropped"] == 1
    out, _, _ = L.lift_rows(lifter(), ["chr16\t600\t700\tSD1\t0\t+\t600\t700\t1,1,1\n"])
    assert out == ["chr16\t595\t695\tSD1\t0\t+\t595\t695\t1,1,1\n"]
