import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_o2personal", os.path.join(HERE, "o2personal.py"))
PC = importlib.util.module_from_spec(spec)
spec.loader.exec_module(PC)


def test_a_gwfam_row_becomes_a_copy_named_family_dot_index(tmp_path):
    p = tmp_path / "c.tsv"
    p.write_text("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n"
                 "SM12\t3\tt\tNC_1\t100\t400\t2\t-\t0\t100-150,300-400\n")
    (r,) = PC.gwfam_rows(str(p))
    assert r == dict(copy="SM12.3", family="SM12", chrom="NC_1", start=100, end=400, strand="-",
                     exons=[(100, 150), (300, 400)])


def test_the_main_transcript_is_its_exons_joined_and_read_on_its_strand():
    ref = "AAAACCCCGGGGTTTT"
    fetch = lambda a, b: ref[a:b]
    assert PC.spliced(fetch, [(0, 2), (8, 10)], "+") == "AAGG"
    assert PC.spliced(fetch, [(0, 2), (8, 10)], "-") == "CCTT"


def test_the_family_split_is_fixed_by_the_name_alone():
    assert PC.half("NPIP") == PC.half("NPIP")
    halves = {PC.half(f"SM{i}") for i in range(20)}
    assert halves == {"dev", "heldout"}


def _read(events, lo=0, hi=100):
    return (set(events), [(lo, hi)])


def test_variant_events_follow_the_corrector_s_explained_and_carried_rules():
    a, b = ("S", 10, "A"), ("S", 20, "G")
    reads = [_read([a, b])] * 9 + [_read([a])] * 1 + [_read([])] * 10     # a in 10/20 reads, b in 9/20
    hom, het = PC.variant_events(reads)
    assert hom == set() and het == {a, b}
    hom, het = PC.variant_events([_read([a])] * 9 + [_read([])])          # 9/10 = carried
    assert hom == {a} and het == set()
    hom, het = PC.variant_events([_read([a])] + [_read([])] * 9)          # one read never counts
    assert hom == set() and het == set()


def test_two_haplotypes_separate_alleles_that_never_share_a_read():
    a, b, c = ("S", 10, "A"), ("S", 50, "T"), ("S", 80, "C")
    reads = [_read([a, b])] * 6 + [_read([c])] * 5
    haps = PC.haplotypes(reads)
    assert sorted(map(sorted, haps)) == sorted([sorted([a, b]), [c]])


def test_one_haplotype_when_the_reads_carry_no_heterozygous_event():
    a = ("S", 10, "A")
    assert PC.haplotypes([_read([a])] * 8) == [{a}]
    assert PC.haplotypes([_read([])] * 8) == [set()]


def test_reads_that_cover_no_het_site_do_not_split_the_copy():
    a = ("S", 10, "A")
    reads = [_read([a], 0, 50)] * 4 + [_read([], 0, 50)] * 4 + [_read([], 60, 100)] * 20
    haps = PC.haplotypes(reads)
    assert sorted(map(sorted, haps)) == sorted([[a], []])


def test_the_consensus_takes_events_most_reads_carry():
    a, b = ("S", 10, "A"), ("S", 20, "G")
    reads = [_read([a, b])] * 6 + [_read([a])] * 4
    assert PC.consensus(reads) == {a, b}
    assert PC.consensus([_read([a])] * 4 + [_read([])] * 6) == set()


def test_events_apply_to_the_reference_in_its_own_coordinates():
    ref = "ACGTACGTAC"                                                     # positions 100..109
    assert PC.apply_events(ref, 100, {("S", 101, "T")}) == "ATGTACGTAC"
    assert PC.apply_events(ref, 100, {("I", 102, "GG")}) == "ACGGGTACGTAC"
    assert PC.apply_events(ref, 100, {("D", 104, 2)}) == "ACGTGTAC"
    assert PC.apply_events(ref, 100, {("S", 101, "T"), ("D", 104, 2), ("I", 108, "A")}) == "ATGTGTAAC"
    assert PC.apply_events(ref, 100, {("S", 50, "T"), ("S", 200, "T")}) == ref             # outside: ignored


def test_conflicting_events_keep_the_more_carried_one():
    s1, s2, d = ("S", 10, "A"), ("S", 10, "G"), ("D", 9, 3)
    assert PC.compatible({s1: 5, s2: 3}) == {s1}
    assert PC.compatible({s1: 2, d: 4}) == {d}
    i1, i2 = ("I", 10, "A"), ("I", 10, "AA")
    assert PC.compatible({i1: 3, i2: 1, s1: 2}) == {i1, s1}               # an insertion before a base, a substitution of it


def test_the_truth_is_the_one_copy_within_one_substitution_of_the_best_placement():
    ref, alt = "ref", "alt"
    assert PC.truth_of([(100, "c", ref), (90, "d", ref), (98, "c", alt)]) == ("determined", "c", "refcons")
    assert PC.truth_of([(100, "c", alt), (95, "c", ref), (90, "d", ref)]) == ("determined", "c", "altonly")
    assert PC.truth_of([(100, "c", alt), (96, "c", ref), (90, "d", ref)]) == ("determined", "c", "refcons")
    assert PC.truth_of([(100, "c", ref), (97, "d", alt)]) == ("undetermined", None, None)
    assert PC.truth_of([(100, None, alt), (90, "c", ref)]) == ("off-copy", None, None)
    assert PC.truth_of([(100, "c", ref), (99, None, alt)]) == ("undetermined", None, None)
    assert PC.truth_of([]) == ("unmapped", None, None)


def test_index_names_map_to_assembly_contigs_by_unique_length():
    fai = {"CM1": 100, "CM2": 200, "JA1": 50, "JA2": 50}
    sq = {"chr1_mat": 100, "chr2_mat": 200, "u1": 50, "CM9": 7}
    assert PC.names_by_length(sq, fai) == {"chr1_mat": "CM1", "chr2_mat": "CM2"}
    assert PC.names_by_length({"CM1": 100}, fai) == {"CM1": "CM1"}           # already the assembly's name


def test_an_alt_ortholog_becomes_a_copy_row_on_its_contig_with_its_strand_turned():
    o = dict(copy="SM1.0", contig="CM5", start=900, end=1300, strand="-", exons=[(900, 1000), (1200, 1300)])
    rows = PC.alt_rows([o], {"SM1.0": "+"})
    assert rows == {"CM5": [("SM1.0", "-", [(900, 1000), (1200, 1300)], (900, 1300))]}
    assert PC.alt_rows([dict(o, strand="+")], {"SM1.0": "-"})["CM5"][0][1] == "-"


def _s(read, fam, truth, mm, corr, origin="altonly", margin=3, **arms):
    return dict(read=read, family=fam, truth=truth, minimap2=mm, corrector=corr, origin=origin, margin=margin, **arms)


def test_an_arm_is_tallied_against_minimap2_on_the_same_reads():
    rows = [_s("r1", "F", "c", "d", "c"), _s("r2", "F", "c", "c", "d"), _s("r3", "F", "c", "c", "c"),
            _s("r4", "F", "c", "d", "e")]
    n = PC.tally(rows, "corrector")
    assert (n["reads"], n["right"], n["fixes"], n["breaks"], n["both wrong"]) == (4, 2, 1, 1, 1)
    assert PC.tally(rows, "minimap2")["right"] == 2


def test_families_count_where_one_arm_is_right_more_often_than_another():
    rows = [_s("a", "F1", "c", "d", "c", personal="c"), _s("b", "F2", "c", "d", "c", personal="d"),
            _s("c", "F3", "c", "c", "c", personal="c")]
    assert PC.families_vs(rows, "personal", "corrector") == (0, 1)
    assert PC.families_fix_break(rows, "corrector") == (2, 0)


def test_a_verdict_needs_every_bar():
    good = dict(gain=1.5, fixes=50, breaks=5, better=12, worse=2, refcons_drop=0.0)
    assert PC.verdict(good, min_gain=1.0) == "PASS"
    assert PC.verdict(dict(good, worse=3), min_gain=1.0) == "PASS"           # 12 / 5 rounded up is 3
    assert PC.verdict(dict(good, worse=4), min_gain=1.0) == "FAIL"
    assert PC.verdict(dict(good, breaks=11), min_gain=1.0) == "FAIL"
    assert PC.verdict(dict(good, refcons_drop=0.3), min_gain=1.0) == "FAIL"
    assert PC.verdict(dict(good, better=10, worse=0), min_gain=1.0) == "FAIL"


def test_a_realignment_decides_only_a_unique_best_over_every_tied_copy():
    assert PC.decide({"c": 120, "d": 110}, {"c", "d"}) == ("c", "decided")
    assert PC.decide({"c": 120, "d": 120}, {"c", "d"}) == (None, "tie")
    assert PC.decide({"c": 120}, {"c", "d"}) == (None, "incomplete")
    assert PC.decide({"c": 120, "d": 110, "e": 130}, {"c", "d"}) == ("c", "decided")   # e is not tied: not a candidate
