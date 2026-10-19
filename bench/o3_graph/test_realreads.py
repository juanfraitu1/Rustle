import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_realreads", os.path.join(HERE, "realreads.py"))
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)


def _fams():
    cp = lambda gid, ch, s, ex, bt="protein_coding": dict(gene_id=gid, gene_name=gid, chrom=ch, start=s, end=ex[-1][1],  # noqa: E731
                                                          strand="+", main="t", transcripts={"t": ex}, biotype=bt)
    return [dict(family="F", copies=[cp("a", "c1", 100, [[100, 300], [500, 700]]), cp("b", "c1", 5000, [[5000, 5200], [5400, 5600]]),
                                     cp("p", "c2", 100, [[100, 400]], "unprocessed_pseudogene")])]


def test_a_read_s_truth_is_the_copy_whose_exons_its_unique_primary_overlaps_most():
    idx = R.copy_exons(_fams())
    assert R.truth_copy("c1", [(150, 300), (500, 650)], 60, idx) == "F.0"
    assert R.truth_copy("c1", [(5100, 5200), (5400, 5420)], 60, idx) == "F.1"
    assert R.truth_copy("c1", [(150, 300)], 59, idx) is None            # not uniquely placed: no truth
    assert R.truth_copy("c1", [(310, 490)], 60, idx) is None            # intronic: on no exon
    assert R.truth_copy("c3", [(150, 300)], 60, idx) is None


def test_the_masked_copy_is_a_seeded_choice_among_expressed_coding_copies():
    counts = {"F.0": 50, "F.1": 25, "F.2": 400}                         # F.2 is a pseudogene: never masked
    got = {R.choose_masked(_fams(), counts, seed=s, min_reads=20)["F"] for s in range(30)}
    assert got == {"a", "b"}
    assert R.choose_masked(_fams(), dict(counts, **{"F.1": 5}), seed=3, min_reads=20) == {"F": "a"}
    assert R.choose_masked(_fams(), {"F.0": 1, "F.1": 2}, seed=3, min_reads=20) == {}   # no copy expressed enough: no test
