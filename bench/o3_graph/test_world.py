import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import world as W  # noqa: E402

TPL = dict(gene="G1", chrom="chr20", strand="+", canonical="tx0", n_exons=6, cassettes=[3],
           isoforms=[dict(tid="t1", skip=[3]), dict(tid="t2", skip=[2, 3]), dict(tid="t3", skip=[6]), dict(tid="t4", skip=[5]),
                     dict(tid="t5", skip=[4])])


def test_the_main_isoform_skips_the_first_cassette_and_outweighs_the_canonical():
    iso = W.curated_isoforms(TPL, cap=3)
    assert iso[0] == {"skip": [3], "weight": W.MAIN_WEIGHT} and W.MAIN_WEIGHT > 1.0
    assert len(iso) == 3 and all(x["weight"] == W.OTHER_WEIGHT for x in iso[1:])


def test_an_inserted_exon_shifts_the_skips_after_it():
    assert W.shift_skips([{"skip": [2, 5], "weight": 2.0}], after=3) == [{"skip": [2, 6], "weight": 2.0}]


def test_rep_is_the_most_expressed_isoform():
    assert W.rep_isoform([{"skip": [3], "weight": 2.0}, {"skip": [5], "weight": 0.5}]) == "iso1"
    assert W.rep_isoform([{"skip": [3], "weight": 0.5}]) == "iso0"


def test_specs_cover_every_scenario_with_its_copies():
    sp = W.specs([TPL], seed=7, genome="g.fa", annotation="a.gtf")
    assert [s["scenario"] for s in sp] == list(W.SCENARIOS)
    by = {s["scenario"]: s for s in sp}
    assert sorted(by["OUTSIDE@0.03"]["kinds"].values()) == ["outside", "present", "present"]
    assert sorted(by["BRANCH"]["kinds"].values()) == ["absent", "present", "present"]
    twin = by["TWIN_ALT"]["copies"][1]
    assert twin["ops"] == [{"op": "snp", "rate": W.TWIN_EXON_RATE, "region": "exon", "exon": 3}]
    e = by["BRANCH"]["copies"][-1]
    assert e["ops"][0]["op"] == "exon_insert" and e["in_reference"] is False
    assert all(s["template"]["transcript"] == "tx0" and s["template"]["intron_cap"] == W.INTRON_CAP for s in sp)


def test_one_world_never_holds_two_families_of_one_template():
    sp = W.specs([TPL, dict(TPL, gene="G2")], seed=7, genome="g.fa", annotation="a.gtf")
    for sc in W.SCENARIOS:
        genes = [s["template"]["gene"] for s in sp if s["scenario"] == sc]
        assert sorted(genes) == ["G1", "G2"]
