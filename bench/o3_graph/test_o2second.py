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
