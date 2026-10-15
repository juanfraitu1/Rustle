import os
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build as B  # noqa: E402
import evidence as E  # noqa: E402


# -------------------------------------------------------------------------------------------------------- evidence

def test_blocks_from_cigar():
    # 5S 10M 2I 5M 3D 4M 100N 6M 2S, from 1000
    al, ex = E.blocks(1000, [(4, 5), (0, 10), (1, 2), (0, 5), (2, 3), (0, 4), (3, 100), (0, 6), (4, 2)])
    assert al == [(1000, 1015), (1018, 1022), (1122, 1128)]
    assert ex == [(1000, 1022), (1122, 1128)]


def write_bam(path, reads):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 100000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, flag, start, cigar in sorted(reads, key=lambda r: r[2]):
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.flag, a.reference_id, a.reference_start, a.mapping_quality = name, flag, 0, start, 60
            a.cigartuples = cigar
            a.query_sequence = "A" * sum(n for op, n in cigar if op in (0, 1, 4, 7, 8))
            out.write(a)
    pysam.index(str(path))


def test_own_reads(tmp_path):
    p = tmp_path / "r.bam"
    write_bam(p, [("prim", 0, 1000, [(0, 100), (3, 500), (0, 100)]),
                  ("second", 256, 1000, [(0, 100)]),
                  ("suppl", 2048, 1000, [(0, 100)]),
                  ("intron_only", 0, 1200, [(0, 100)]),
                  ("prim2", 16, 1550, [(0, 100)])])
    with pysam.AlignmentFile(str(p)) as bam:
        got = E.own_reads(bam, "chr1", [(1000, 1100), (1600, 1700)])
    assert got == [("prim", [(1000, 1100), (1600, 1700)]), ("prim2", [(1550, 1650)])]


# ----------------------------------------------------------------------------------------------------------- build

PROJ = [(100, 200), (300, 400), (500, 600)]


def test_simple_path_is_one_transcript():
    loci = B.build([("r1", [(110, 200), (300, 390)]), ("r2", [(320, 400), (500, 580)])], PROJ)
    assert loci == [{"transcripts": [[(110, 200), (300, 400), (500, 580)]], "reads": 2,
                     "exons": [(110, 200), (300, 400), (500, 580)]}]


def test_branching_keeps_distinct_read_chains():
    reads = [("r1", [(100, 200), (300, 400)]), ("r2", [(100, 200), (500, 600)]), ("r3", [(310, 390)])]
    loci = B.build(reads, PROJ)
    assert loci[0]["transcripts"] == [[(100, 200), (300, 400)], [(100, 200), (500, 600)]] and loci[0]["reads"] == 3


def test_a_dropped_block_breaks_the_link():
    # r1's middle block (700-750) lies outside the projection, so its first and last blocks are not linked.
    reads = [("r1", [(100, 200), (700, 750), (500, 600)]), ("r2", [(120, 200)]), ("r3", [(500, 580)])]
    loci = B.build(reads, PROJ)
    assert [lc["exons"] for lc in loci] == [[(100, 200)], [(500, 600)]] and [lc["reads"] for lc in loci] == [2, 2]


def test_a_component_with_one_read_is_not_emitted():
    loci = B.build([("r1", [(100, 200), (300, 400)]), ("r2", [(500, 560)])], PROJ)
    assert loci == []


def test_single_exon_copy():
    loci = B.build([("r1", [(100, 150)]), ("r2", [(140, 200)])], PROJ)
    assert loci == [{"transcripts": [[(100, 200)]], "reads": 2, "exons": [(100, 200)]}]
