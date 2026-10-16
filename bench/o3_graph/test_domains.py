import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


D = _load("domains")


def test_coverage_is_the_union_of_aligned_blocks_over_the_length():
    assert D.union_cover([(0, 50), (40, 100), (200, 250)], 500) == 150 / 500
    assert D.union_cover([], 100) == 0.0


def test_a_gene_pair_is_graded_by_how_much_of_each_gene_s_transcripts_its_partner_s_locus_holds():
    assert D.grade(0.95, 0.9) == "whole"            # each locus holds the other's transcript
    assert D.grade(0.97, 0.4) == "fragment"         # one is a truncated copy of the other
    assert D.grade(0.6, 0.3) == "partial"
    assert D.grade(0.3, 0.2) == "segment"           # they share only a segment: a domain sharer


def test_transcript_to_locus_coverage_takes_the_best_isoform_and_the_identity_floor():
    # PAF rows (query isoform, qlen, qs, qe, target locus, identity): gene of an isoform = the part before '|'
    rows = [("A|t1", 1000, 0, 900, "B", 0.98), ("A|t2", 800, 0, 400, "B", 0.99), ("B|u", 1200, 0, 1200, "A", 0.85),
            ("B|u", 1200, 0, 600, "A", 0.95)]
    cov = D.tx_locus_cover(rows, min_ident=0.9)
    assert cov[("A", "B")] == 0.9 and cov[("B", "A")] == 0.5
