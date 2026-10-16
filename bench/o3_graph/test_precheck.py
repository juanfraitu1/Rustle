import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


P = _load("precheck")


def test_the_seed_is_the_smallest_with_enough_detectable_removed_copies():
    assert P.pick_seed({81: 1, 82: 3, 83: 5}, need=3) == 82
    assert P.pick_seed({81: 2, 82: 2}, need=3) is None


def test_a_copy_s_divergence_is_its_best_alignment_away_from_its_own_locus():
    alns = [dict(chrom="c1", start=1000, end=3000, score=900, de=0.0), dict(chrom="c1", start=50000, end=52000, score=850, de=0.012),
            dict(chrom="c2", start=0, end=2000, score=700, de=0.03)]
    assert P.away_divergence(alns, ("c1", 1000, 3000)) == 0.012
    assert P.away_divergence(alns[:1], ("c1", 1000, 3000)) == 1.0          # nothing else in the genome explains it


def test_sam_lines_parse_into_alignments_with_their_divergence():
    sam = ("@SQ\tSN:c1\tLN:100000\n"
           "q1\t0\tc1\t1001\t60\t100M300N100M\t*\t0\t0\t*\t*\tNM:i:2\tAS:i:180\tde:f:0.0123\n"
           "q1\t256\tc2\t51\t0\t200M\t*\t0\t0\t*\t*\tAS:i:150\tde:f:0.03\n"
           "q2\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n")
    alns = P.parse_sam(sam)
    assert set(alns) == {"q1"}
    assert alns["q1"][0] == dict(chrom="c1", start=1000, end=1500, score=180, de=0.0123)
    assert alns["q1"][1]["de"] == 0.03
