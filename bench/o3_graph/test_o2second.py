import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_o2second", os.path.join(HERE, "o2second.py"))
S2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S2)

CPS = {"chr5": [("F.0", "+", [(100, 200), (300, 400)], (100, 400)), ("F.1", "+", [(900, 1000)], (900, 1000))]}


def test_the_total_variation_between_two_copy_shares():
    assert S2.tvd({"F.0": 3, "F.1": 1}, {"F.0": 1, "F.1": 1}) == 0.25
    assert S2.tvd({}, {"F.0": 1}) == 1.0 and S2.tvd({"F.0": 2}, {"F.0": 5}) == 0.0


def test_assembled_transcripts_weigh_on_the_family_copy_their_exons_lie_on():
    tx = {"t1": dict(chrom="chr5", blocks=[(150, 200), (300, 350)], strand="+"),
          "t2": dict(chrom="chr5", blocks=[(910, 990)], strand="+"),
          "t3": dict(chrom="chr9", blocks=[(0, 10)], strand="+")}
    assert S2.copy_weights(tx, {"t1": 5.0, "t2": 2.0, "t3": 9.0}, CPS, {"F.0", "F.1"}) == {"F.0": 5.0, "F.1": 2.0}


def test_an_attribute_is_read_by_its_own_name_not_inside_a_longer_one(tmp_path):
    g = tmp_path / "a.gtf"
    g.write_text('chr5\tX\ttranscript\t1\t9\t.\t+\t.\tgene_id "g"; transcript_id "t1"; matched_reads "7"; longcov "3.0"; '
                 'reads "18"; cov "30.0";\n')
    assert S2._weights(str(g), "reads", None) == {"t1": 18.0} and S2._weights(str(g), "cov", None) == {"t1": 30.0}


def test_a_read_is_credited_to_the_copy_of_the_transcript_covering_at_least_half_its_aligned_bases():
    idx = S2.transcript_index({"t1": dict(chrom="chr5", blocks=[(100, 200), (300, 400)], strand="+", copy="F.0"),
                               "t2": dict(chrom="chr5", blocks=[(150, 160)], strand="+", copy="F.1")})
    assert S2.credit("chr5", [(120, 200), (300, 330)], idx) == "F.0"                   # 110 of 110 bases in t1
    assert S2.credit("chr5", [(180, 260)], idx) is S2.NOWHERE                           # 20 of 80: no locus covers it
    assert S2.credit("chr9", [(0, 50)], idx) is S2.NOWHERE


def test_each_arm_admits_the_records_its_read_pool_takes():
    recs = [dict(chrom="chr5", blocks=[(100, 200)], AS=1000, prim=True),
            dict(chrom="chr5", blocks=[(900, 1000)], AS=990, prim=False),
            dict(chrom="chr5", blocks=[(5000, 5100)], AS=900, prim=False)]
    assert [r["AS"] for r in S2.admitted(recs, "primary")] == [1000]
    assert [r["AS"] for r in S2.admitted(recs, "good")] == [1000, 990]                 # AS >= 0.98 x best
    assert [r["AS"] for r in S2.admitted(recs, "home", ("chr5", 900, 1000))] == [990]  # only its records at home
    assert [r["AS"] for r in S2.admitted(recs, "home", None)] == [1000, 990]           # unhomed: the good pool


def test_a_tied_read_is_right_only_when_credited_to_its_true_copy_alone():
    assert S2.verdict({"F.1"}, "F.1") == "right" and S2.verdict({"F.0"}, "F.1") == "wrong"
    assert S2.verdict({"F.0", "F.1"}, "F.1") == "echo" and S2.verdict(set(), "F.1") == "uncredited"
    assert S2.verdict({S2.NOWHERE}, "F.1") == "uncredited" and S2.verdict({S2.NOWHERE, "F.1"}, "F.1") == "right"


def test_families_where_one_arm_credits_a_larger_share_of_reads_right_than_another():
    a = {"F": {"right": 8, "wrong": 2}, "G": {"right": 3, "echo": 7}, "H": {"right": 2}}
    b = {"F": {"right": 5, "wrong": 5}, "G": {"right": 4, "wrong": 6}, "H": {"right": 1, "wrong": 1}}
    assert S2.families_better(a, b, min_reads=10) == (["F"], ["G"])          # H has 2 reads
