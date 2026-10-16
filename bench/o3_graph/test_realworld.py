import importlib.util
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R = _load("realworld")
GENOME = "".join(random.Random(3).choice("ACGT") for _ in range(2000))


def test_snps_sit_on_genomic_positions_so_every_isoform_agrees():
    snps = {105: "A", 410: "T"}
    iso1 = R.splice(GENOME, [(100, 200), (400, 500)], "+", snps)
    iso2 = R.splice(GENOME, [(400, 500)], "+", snps)
    assert iso1[5] == "A" and iso1[100 + 10] == "T" and iso2[10] == "T"


def test_minus_strand_transcripts_are_reverse_complemented():
    s = R.splice(GENOME, [(100, 110), (200, 210)], "-", {})
    comp = str.maketrans("ACGT", "TGCA")
    assert s == (GENOME[100:110] + GENOME[200:210])[::-1].translate(comp)


def test_the_novel_exon_goes_into_every_isoform_whose_intron_spans_it():
    novel = "N" * 7
    p = 250                                                     # inside the main isoform's intron 200-300
    with_j = R.splice(GENOME, [(100, 200), (300, 400), (500, 600)], "+", {}, insert_at=p, novel=novel)
    assert with_j == GENOME[100:200] + novel + GENOME[300:400] + GENOME[500:600]
    skipping = R.splice(GENOME, [(100, 200), (500, 600)], "+", {}, insert_at=p, novel=novel)
    assert skipping == GENOME[100:200] + novel + GENOME[500:600]          # a real new exon is spliced in here too
    upstream = R.splice(GENOME, [(500, 600), (700, 800)], "+", {}, insert_at=p, novel=novel)
    assert novel not in upstream                                          # no intron of this isoform spans it


def test_anonymous_names_are_a_seeded_shuffle_carrying_no_truth():
    rows = [("FAM.0", "t1"), ("FAM.0", "t2"), ("FAM.1", "t1")]
    a, b = R.anon_names(rows, seed=4), R.anon_names(rows, seed=4)
    assert a == b and sorted(a) == ["s0000001", "s0000002", "s0000003"]


def test_the_novel_flag_sees_the_exon_on_either_strand():
    novel = "ACGTTGCAAC"
    assert R.carries(novel, "TT" + novel + "GG") and R.carries(novel, "TT" + novel[::-1].translate(R.COMP) + "GG")
    assert not R.carries(novel, "TTGG") and not R.carries("", "ACGT")
