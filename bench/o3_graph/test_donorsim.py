import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_donorsim", os.path.join(HERE, "donorsim.py"))
DS = importlib.util.module_from_spec(spec)
spec.loader.exec_module(DS)


def _row(strand, mapq=60, tp="P", qs=0, qe=100, qlen=100, ts=1000, te=1100, cg="100M", t="ctg1"):
    return dict(qlen=qlen, qs=qs, qe=qe, strand=strand, tname=t, ts=ts, te=te, mapq=mapq, tp=tp, cg=cg)


def test_a_copy_s_exons_lift_through_its_flanked_alignment_on_either_strand():
    assert DS.lift_exons(_row("+"), [(10, 20), (40, 60)]) == [(1010, 1020), (1040, 1060)]
    assert DS.lift_exons(_row("-"), [(10, 20)]) == [(1080, 1090)]          # the query's reverse complement runs forward
    assert DS.lift_exons(_row("+", cg="50M10D50M"), [(10, 20), (60, 70)]) == [(1010, 1020), (1070, 1080)]
    assert DS.lift_exons(_row("+", qs=30), [(10, 20)]) is None                # an exon outside the alignment


def test_the_ortholog_is_the_longest_confident_primary_through_which_every_isoform_lifts():
    rows = [_row("+", mapq=60, ts=1000, te=1100), _row("+", mapq=60, tp="S", ts=5000, te=5100),
            _row("+", mapq=60, qs=0, qe=40, ts=9000, te=9040, cg="40M")]
    got = DS.pick_ortholog(rows, {"t1": [(10, 20)], "t2": [(30, 35)]})
    assert got["tname"] == "ctg1" and got["isoforms"] == {"t1": [(1010, 1020)], "t2": [(1030, 1035)]} and got["strand"] == "+"
    assert DS.pick_ortholog([_row("+", mapq=3)], {"t1": [(10, 20)]}) is None          # not confidently placed
    assert DS.pick_ortholog(rows, {"t1": [(90, 150)]}) is None                          # an isoform does not lift


def test_two_copies_whose_orthologs_overlap_on_a_haplotype_are_both_dropped():
    orth = {("F.0", "mat"): dict(tname="c1", start=100, end=200), ("F.1", "mat"): dict(tname="c1", start=150, end=250),
            ("F.2", "mat"): dict(tname="c1", start=300, end=400), ("F.0", "pat"): dict(tname="c9", start=100, end=200)}
    assert sorted(DS.drop_ambiguous(orth)) == [("F.0", "pat"), ("F.2", "mat")]


def test_an_exon_whose_lifted_length_balloons_through_a_large_target_gap_does_not_lift():
    # 50M 40000D 50M: an exon spanning the deletion would lift as 40,020 bp; the isoform is not projectable
    assert DS.lift_exons(_row("+", cg="50M40000D50M", te=41100), [(40, 60)]) is None
    assert DS.lift_exons(_row("+", cg="50M3D50M", te=1103), [(40, 60)]) == [(1040, 1063)]    # a small indel is fine


def test_the_exon_growth_guard_can_be_lifted_for_a_confident_primary():
    row = _row("+", cg="10M200D90M")                                         # an exon end lifts across a 200-bp gap
    assert DS.lift_exons(row, [(5, 20)]) is None
    assert DS.lift_exons(row, [(5, 20)], max_growth=None) == [(1005, 1220)]


def test_an_isoform_s_reads_split_over_the_copy_s_haplotypes_unless_each_instance_gets_its_own():
    assert DS.isoform_reads(100, 1.0, 2.0, 2.0, 2) == 50
    assert DS.isoform_reads(100, 1.0, 2.0, 2.0, 2, per_instance=True) == 100
    assert DS.isoform_reads(150, 1.0, 0.5, 2.5, 3) == 10
