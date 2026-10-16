import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


C = _load("calibrate")
REF = "ACGTAAAAACGTACGTCCCGTACGTACGTACGT" * 4             # an A homopolymer at 4-8, a C run at 16-18


def test_cigar_events_are_counted_with_their_homopolymer_context():
    # 4= 1D(A, inside the AAAAA run) 5= 1X 2= 1I(T, not a run) 10= 50N 3= 5S
    cig = [(7, 4), (2, 1), (7, 5), (8, 1), (7, 2), (1, 1), (7, 10), (3, 50), (7, 3), (4, 5)]
    s = C.cigar_events(cig, 0, "A" * 200, lambda p: REF[p % len(REF)])
    assert s["aligned"] == 25 and s["mismatch"] == [(10, "A")]
    assert s["del"] == [(4, 1, True)] and s["ins"] == [(13, 1, False)] and s["clip"] == 5


def test_a_mismatch_shared_by_enough_reads_is_a_donor_variant_not_an_error():
    counts = {("c1", 100, "T"): 6, ("c1", 200, "G"): 1, ("c1", 300, "A"): 9}
    cov = {("c1", 100): 12, ("c1", 200): 15, ("c1", 300): 10}
    v = C.classify_mismatches(counts, cov)
    assert v[("c1", 100, "T")] == "het" and v[("c1", 300, "A")] == "hom" and v[("c1", 200, "G")] == "error"


def test_truncation_is_measured_in_transcript_bases_on_the_transcript_s_strand():
    exons = [(100, 200), (300, 400), (500, 600)]                 # 300 transcript bases
    # a read covering 150-200 and 300-350: 50 bases missing at the plus-strand 5' end, 150 at the 3' end
    assert C.transcript_ends([(150, 200), (300, 350)], exons, "+") == (50, 150)
    assert C.transcript_ends([(150, 200), (300, 350)], exons, "-") == (150, 50)


def test_bases_aligned_inside_an_annotated_intron_are_retention():
    exons = [(100, 200), (300, 400)]
    assert C.intron_bases([(150, 260), (300, 350)], exons) == 60
    assert C.intron_bases([(150, 200), (300, 350)], exons) == 0


def test_an_intron_is_retained_only_when_one_aligned_block_spans_all_of_it():
    exons = [(100, 200), (300, 400), (500, 600)]
    assert C.retained_introns([(150, 450)], exons) == 1                  # spans 200-300 whole, 400-500 only partly
    assert C.retained_introns([(150, 260), (300, 350)], exons) == 0      # an alternative exon inside the intron
