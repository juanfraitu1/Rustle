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


def test_an_alignment_s_edits_are_its_nm_plus_the_read_bases_it_leaves_unaligned():
    assert PC.edits_of([(4, 3), (0, 50), (3, 1000), (0, 40), (5, 2)], 7) == 12      # 7 + 3 soft + 2 hard; the intron is free
    assert PC.edits_of([(0, 90)], 0) == 0


def test_the_truth_is_the_one_copy_with_strictly_fewest_edits_in_the_animal_s_genome():
    ref, alt = "ref", "alt"
    assert PC.truth_of([(0, "c", ref), (24, "d", ref), (2, "c", alt)]) == ("determined", "c", "refcons")
    assert PC.truth_of([(3, "c", alt), (4, "c", ref), (5, "d", ref)]) == ("determined", "c", "altonly")
    assert PC.truth_of([(3, "c", alt), (3, "c", ref), (5, "d", ref)]) == ("determined", "c", "refcons")
    assert PC.truth_of([(2, "c", ref), (2, "d", alt)]) == ("undetermined", None, None)
    assert PC.truth_of([(1, None, alt), (5, "c", ref)]) == ("off-copy", None, None)
    assert PC.truth_of([(1, "c", ref), (1, None, alt)]) == ("undetermined", None, None)
    assert PC.truth_of([]) == ("unmapped", None, None)


def test_window_edits_count_bases_substituted_inserted_and_deleted():
    assert PC.window_edits({("S", 10, "A"), ("I", 20, "GG"), ("D", 30, 3)}) == 6
    assert PC.window_edits(set()) == 0


def test_the_fewest_and_the_highest_take_a_read_only_when_unique():
    assert PC.fewest({"c": 1, "d": 3}, "d") == "c"
    assert PC.fewest({"c": 1, "d": 1}, "d") == "d"                                 # shared: minimap2's pick
    assert PC.highest({"c": 1204, "d": 1196}, "d") == "c"
    assert PC.highest({"c": 1204, "d": 1204}, "d") == "d"


def test_a_tie_s_structure_is_whether_its_copies_differ_in_intron_count():
    assert PC.structure({"c": 3, "d": 0}) == "intron-differ"
    assert PC.structure({"c": 3, "d": 3}) == "same"


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


def test_scored_reads_are_the_half_s_determined_reads_whose_truth_is_tied():
    dom = {"a": ("F", {"c", "d"}, "d", "d", 3), "b": ("F", {"c", "d"}, "d", "d", 0), "x": ("G", {"c", "d"}, "d", "d", 0),
           "u": ("F", {"c", "d"}, "d", "d", 0), "o": ("F", {"c", "d"}, "d", "d", 0)}
    t = lambda st, tr="", inn="", o="": dict(status=st, truth=tr, in_tie=inn, origin=o)
    truth = {"a": t("determined", "c", "1", "refcons"), "b": t("determined", "d", "1", "altonly"), "x": t("determined", "c", "1"),
             "u": t("undetermined"), "o": t("determined", "e", "0", "altonly")}
    qual = {"F": dict(half="dev", min_identity=0.97), "G": dict(half="heldout", min_identity=0.97)}
    arms = {q: dict(fewest_edits="c", max_as="d", structure="same") for q in dom}
    status, rows = PC.score_rows(dom, truth, qual, arms, {"a": dict(corrector="c")}, "dev")
    assert [r["read"] for r in rows] == ["a", "b"]
    assert rows[0]["fewest_subs"] == "c" and rows[1]["fewest_subs"] == "d"          # no run row: minimap2's pick
    assert status == {"determined refcons": 1, "determined altonly": 1, "undetermined": 1,
                      "determined altonly, truth outside the tie": 1}


def test_the_bars_of_an_arm_count_its_fixes_breaks_families_and_refcons_harm():
    rows = [dict(family=f"F{i}", truth="c", minimap2="d", o2="c", origin="altonly") for i in range(12)]
    rows += [dict(family="F0", truth="c", minimap2="c", o2="c", origin="refcons")] * 50
    b = PC.bars_for(rows, "o2")
    assert (b["fixes"], b["breaks"], b["better"], b["worse"], b["refcons_drop"]) == (12, 0, 12, 0, 0.0)
    assert b["gain"] == round(100 * 12 / 62, 2) and b["verdict"] == "PASS"
    rows[-1] = dict(rows[-1], o2="d")                                               # one refcons read broken
    assert PC.bars_for(rows, "o2")["refcons_drop"] == 2.0 and PC.bars_for(rows, "o2")["verdict"] == "FAIL"


def test_match_columns_read_the_eqx_cigar_in_the_primary_s_orientation():
    cig = [(4, 2), (7, 3), (8, 1), (2, 5), (7, 2), (1, 1), (8, 1)]                 # 2S 3= 1X 5D 2= 1I 1X  (10 read bases)
    eq, x = PC.match_columns(cig, False, 10)
    assert eq == {2, 3, 4, 6, 7} and x == {5, 9}
    eq, x = PC.match_columns(cig, True, 10)                                        # the alignment runs on the reverse strand
    assert eq == {7, 6, 5, 3, 2} and x == {4, 0}


def test_the_fewest_edit_copy_is_contradicted_where_its_mismatch_is_a_rival_s_match():
    cols = {"w": ({1, 2, 3}, {10}), "r": ({1, 2, 3, 10}, {20})}
    assert PC.contradicted("w", cols, (0, 30)) is True
    assert PC.contradicted("w", cols, (0, 5)) is False                             # outside the common window
    cols["r"] = ({1, 2, 3}, {10, 20})                                              # a shared error is no contradiction
    assert PC.contradicted("w", cols, (0, 30)) is False


def test_the_primary_oriented_common_window_is_where_every_tied_alignment_aligns():
    assert PC.window_primary({"a": ((0, 100), False), "b": ((10, 90), True)}, 100) == (10, 90)
    assert PC.window_primary({"a": ((5, 100), False), "b": ((0, 80), True)}, 100) == (20, 100)


def test_overlapping_counterparts_go_to_the_copy_with_the_lower_divergence():
    o = lambda s, e, de: dict(tname="m12", start=s, end=e, de=de)
    orth = {("A", "mat"): o(100, 200, 0.0086), ("B", "mat"): o(150, 260, 0.0017), ("C", "mat"): o(500, 600, 0.01)}
    kept = PC.resolve_by_de(orth)
    assert set(kept) == {("B", "mat"), ("C", "mat")}
