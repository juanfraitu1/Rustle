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


def test_tag_summary_reports_per_kb_rates_and_tag_medians():
    rows = [dict(de=0.002, nm=3, aligned=1000, x=1, ins=1, dele=1, hp=1, indels=2, clip=0, qlen=1000),
            dict(de=0.004, nm=6, aligned=2000, x=2, ins=2, dele=2, hp=3, indels=4, clip=20, qlen=2020)]
    s = C.tag_summary(rows)
    assert s["reads"] == 2 and s["de_median"] == 0.003 and s["x_per_kb"] == 1.0 and s["ins_per_kb"] == 1.0
    assert s["hp_indel_frac"] == 4 / 6 and abs(s["clip_frac"] - 20 / 3020) < 1e-12


def test_the_expression_spread_leaves_out_unread_copies_and_poisson_counting_noise():
    import math
    sd = C.expression_log_sd([100, 100, 100, 0])           # equal expression; the unread copy is silent, not spread
    assert sd == 0.0
    counts = [50, 200]                                      # observed log-sd 0.693; Poisson noise of log n ~ 1/n
    raw = math.log(4) / 2
    assert abs(C.expression_log_sd(counts) - math.sqrt(raw ** 2 - (1 / 50 + 1 / 200) / 2)) < 1e-9


def test_error_stats_separate_a_shared_donor_allele_from_a_read_s_own_error():
    def rd(mm):
        return dict(ch="c1", qlen=100, ev=dict(aligned=100, mismatch=mm, ins=[], clip=0, blocks=[(0, 100)], **{"del": []}))
    reads = [rd([(50, "T"), (10, "G")]), rd([(50, "T")]), rd([(50, "T")])]
    s = C.error_stats(reads)
    assert s["hom_per_base"] == 3 / 300 and s["sub_error"] == 1 / 300 and s["het_per_base"] == 0


def test_the_adapter_consensus_is_taken_at_the_most_common_clip_length():
    clips = ["ACGTA", "ACGTA", "ACCTA", "ACGT", "TTTTTTT"]
    assert C.adapter_consensus(clips) == "ACGTA"


def test_a_read_counts_for_a_copy_only_when_an_aligned_block_overlaps_one_of_its_exons():
    exons = [(100, 200), (500, 600)]
    assert C.overlaps_exons([(150, 180)], exons)
    assert C.overlaps_exons([(10, 50), (590, 700)], exons)
    assert not C.overlaps_exons([(10, 100), (600, 900)], exons)          # half-open: touching is no overlap
    assert not C.overlaps_exons([(10, 90), (210, 490)], exons)           # a spliced read whose intron spans the exons


def test_paralog_expression_is_the_silent_fraction_and_the_median_ratio_to_the_family_s_coding_copies():
    rows = [("F1", "protein_coding", 40), ("F1", "protein_coding", 60), ("F1", "unprocessed_pseudogene", 0),
            ("F1", "unprocessed_pseudogene", 25), ("F1", "lncRNA", 10),
            ("F2", "protein_coding", 20), ("F2", "processed_pseudogene", 5),
            ("F3", "protein_coding", 3), ("F3", "unprocessed_pseudogene", 9)]       # coding copies too thin to compare: out
    got = C.paralog_expression(rows, min_coding=10)
    assert got["members"] == 4 and got["families"] == 2
    assert got["silent_frac"] == 0.25                                              # F1's first pseudogene
    assert got["ratio_median"] == 0.25                                             # 25/50, 10/50, 5/20 -> median 0.25
    assert got["by_class"]["unprocessed_pseudogene"] == dict(members=2, silent=1, ratio_median=0.5)
    assert got["by_class"]["lncRNA"] == dict(members=1, silent=0, ratio_median=0.2)
