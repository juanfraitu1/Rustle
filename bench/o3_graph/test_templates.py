import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import templates as T  # noqa: E402

CANON = [(1, 100), (200, 300), (400, 500), (600, 700)]


def rec(strand, **tx):
    return {"chrom": "c", "strand": strand, "transcripts": dict(canon=CANON, **tx)}


def test_an_exon_skipping_isoform_is_a_skip_set_of_the_canonical():
    r = rec("+", skip2=[(1, 100), (400, 500), (600, 700)])
    assert T.cassette_isoforms(r) == ("canon", [("skip2", (2,))])


def test_alternative_ends_map_to_the_same_exons_and_a_novel_exon_is_rejected():
    r = rec("+", alt_tss=[(50, 100), (200, 300), (400, 500), (600, 700)],      # same chain: dropped
            novel=[(1, 100), (250, 300), (400, 500), (600, 650)],               # internal exon not in the canonical
            short_end=[(1, 100), (200, 300), (400, 450)])                       # ends inside exon 3: last exon starts at 400
    assert T.cassette_isoforms(r) == ("canon", [("short_end", (4,))])


def test_minus_strand_indices_follow_the_transcript_orientation():
    r = {"chrom": "c", "strand": "-", "transcripts": {"canon": [(1, 100), (200, 300), (400, 500)],
                                                      "skip_mid": [(1, 100), (400, 500)]}}
    assert T.cassette_isoforms(r) == ("canon", [("skip_mid", (2,))])


def test_internal_cassettes_are_skips_flanked_by_kept_exons():
    assert T.internal_cassettes([(2,), (4,), (2, 3), (1, 2)], 4) == [2, 3]


def test_assembled_genes_are_not_curated():
    assert not T.curated_name("MSTRG.3427") and T.curated_name("ZNF217") and T.curated_name("AL118508.1")
