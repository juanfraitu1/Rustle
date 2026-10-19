import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


P = _load("o2proof")

LOCI = [dict(copy="F.0", family="F", chrom="c1", start=1000, end=3000),
        dict(copy="F.1", family="F", chrom="c1", start=5000, end=7000),
        dict(copy="F.2", family="F", chrom="c2", start=100, end=2100),
        dict(copy="G.0", family="G", chrom="c1", start=9000, end=9500)]


def test_an_alignment_belongs_to_the_copy_whose_locus_it_overlaps_most():
    idx = P.locus_index(LOCI)
    assert P.copy_of(idx, "c1", 1200, 2900) == "F.0"
    assert P.copy_of(idx, "c1", 2900, 5200) == "F.1"          # 100 bp in F.0, 200 bp in F.1
    assert P.copy_of(idx, "c2", 50, 300) == "F.2"
    assert P.copy_of(idx, "c1", 3500, 4500) is None          # between loci
    assert P.copy_of(idx, "c3", 0, 100) is None
    # a locus that is no copy is keyed by its 10-kb bin, so its alignments compete as one locus
    assert P.locus_key(idx, "c1", 3500, 4500) == ("c1", 0)
    assert P.locus_key(idx, "c3", 25000, 25100) == ("c3", 2)
    assert P.locus_key(idx, "c1", 1200, 2900) == "F.0"


def test_a_copy_s_cost_is_its_best_alignment_s_edits_plus_unaligned_bases():
    # (copy, NM, aligned query bases) per alignment of a 1,000-base read
    alns = [("F.0", 3, 1000), ("F.1", 3, 990), ("F.1", 2, 900), ("G.0", 0, 400)]
    assert P.copy_costs(alns, 1000) == {"F.0": 3, "F.1": 13, "G.0": 600}


def test_the_margin_of_the_best_rival_classes_the_read():
    # O2 decides among the family's copies: undecidable = a family rival as good or better; weak = one edit; decidable =
    # two or more (or no family rival aligned); true_unaligned = no alignment on the true copy
    fam = {"F.0", "F.1", "F.2"}
    assert P.decidability({"F.0": 3, "F.1": 3}, "F.0", fam) == ("undecidable", 0)
    assert P.decidability({"F.0": 4, "F.1": 3}, "F.0", fam) == ("undecidable", -1)
    assert P.decidability({"F.0": 3, "F.1": 4}, "F.0", fam) == ("weak", 1)
    assert P.decidability({"F.0": 3, "F.1": 5, "F.2": 9}, "F.0", fam) == ("decidable", 2)
    assert P.decidability({"F.0": 3}, "F.0", fam) == ("decidable", None)
    assert P.decidability({"F.1": 3}, "F.0", fam) == ("true_unaligned", None)
    # a locus OFF the family that explains the read as well (another family's copy, a parent gene, an unlisted paralog:
    # round 13's UBE2Q2P, set 9's 49 USP17L genes outside the family) is O1's question, classed apart
    assert P.decidability({"F.0": 3, "G.0": 0}, "F.0", fam) == ("off_family", -3)
    assert P.decidability({"F.0": 3, ("c9", 12): 3, "F.1": 9}, "F.0", fam) == ("off_family", 0)
    assert P.decidability({"F.0": 3, ("c9", 12): 4, "F.1": 5}, "F.0", fam) == ("decidable", 2)


def test_an_outcome_is_correct_wrong_or_abstain():
    assert P.outcome(None, "F.0") == "abstain"
    assert P.outcome("F.0", "F.0") == "correct"
    assert P.outcome("F.1", "F.0") == "wrong"
    assert P.outcome("wrong-family", "F.0") == "wrong"


def test_the_scorecard_tallies_outcomes_by_method_class_and_mapq():
    rows = [dict(cls="undecidable", mapq0=True, calls={"m": "F.1", "o": None}, copy="F.0"),
            dict(cls="decidable", mapq0=True, calls={"m": "F.0", "o": "F.0"}, copy="F.0"),
            dict(cls="decidable", mapq0=False, calls={"m": "F.0", "o": "F.0"}, copy="F.0")]
    card = P.scorecard(rows, ["m", "o"])
    assert card["m"]["undecidable"]["mapq0"] == {"correct": 0, "wrong": 1, "abstain": 0}
    assert card["o"]["undecidable"]["mapq0"] == {"correct": 0, "wrong": 0, "abstain": 1}
    assert card["o"]["decidable"]["all"] == {"correct": 2, "wrong": 0, "abstain": 0}
    assert card["o"]["decidable"]["mapq0"] == {"correct": 1, "wrong": 0, "abstain": 0}


def test_the_scoped_systems_leave_untied_reads_to_minimap2_and_keep_a_tied_read_s_call_only_on_a_tied_copy():
    rows = [dict(read="t1", calls=dict(mm2="F.1", current="F.0", integrated="F.2")),     # tied on F.0 and F.1
            dict(read="u1", calls=dict(mm2="F.2", current=None, integrated="F.0"))]     # untied: minimap2 resolves it
    P.scope_rows(rows, {"t1": {"F.0", "F.1"}})
    assert rows[0]["calls"]["current_tied"] == "F.0" and rows[0]["calls"]["integrated_tied"] is None
    assert rows[1]["calls"]["current_tied"] == "F.2" and rows[1]["calls"]["integrated_tied"] == "F.2"
    assert P.score_args(["out.json", "W1", "--o2", "x.tsv", "W2", "--tied", "t.tsv"]) == ("out.json", ["W1", "W2"], "x.tsv", "t.tsv")
    assert P.score_args(["out.json", "W1"]) == ("out.json", ["W1"], "o2.integv.m2.tsv", None)
