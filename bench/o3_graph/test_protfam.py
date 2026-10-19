import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_protfam", os.path.join(HERE, "protfam.py"))
P = importlib.util.module_from_spec(spec)
spec.loader.exec_module(P)

COMP = str.maketrans("ACGT", "TGCA")


def _fetch(seq):
    return lambda s, e: seq[s:e]


def test_a_cds_is_translated_across_its_blocks_on_either_strand_from_its_phase():
    g = "NN" + "ATGAA" + "CCCC" + "ATTTTAA" + "NN"                    # ATG AA|A TTT TAA split over two blocks
    assert P.translate_cds([(2, 7), (11, 18)], "+", _fetch(g)) == "MKF"  # the trailing stop is dropped
    rc = g[::-1].translate(COMP)                                    # the same gene on the minus strand
    n = len(g)
    assert P.translate_cds([(n - 18, n - 11), (n - 7, n - 2)], "-", _fetch(rc)) == "MKF"
    assert P.translate_cds([(1, 7), (11, 18)], "+", _fetch(g), phase=1) == "MKF"     # one base to skip first
    assert P.translate_cds([(2, 7), (11, 18)], "+", _fetch(g.replace("ATTTTAA", "ATAGTAA"))) == "MKX"   # internal stop: X


def test_each_gene_keeps_its_longest_cds():
    recs = [("G1", "T1", (0, 90)), ("G1", "T2", (0, 60)), ("G1", "T2", (100, 160)), ("G2", "T3", (5, 50))]
    assert P.longest_cds(recs) == {"G1": "T2", "G2": "T3"}


def test_two_proteins_are_homologous_at_half_identity_over_a_third_of_each():
    lens = {"a": 400, "b": 300, "c": 1000, "d": 200}
    hits = [("a", "b", 55.0, 1, 150, 1, 120),                        # 150/400 = 0.375 and 120/300 = 0.40
            ("a", "c", 80.0, 1, 380, 1, 300),                        # covers a, but only 0.30 of c: a domain, no family
            ("b", "d", 49.0, 1, 200, 1, 200),                        # identity under 50%
            ("a", "a", 100.0, 1, 400, 1, 400)]                       # self
    assert P.protein_links(hits, lens) == {("a", "b")}


def test_families_are_components_of_protein_and_near_identical_links_with_subfamilies_nested():
    prot = {("A", "B")}
    near = [("B", "B2"), ("C", "Cp")]                                 # B2: a near-identical copy; Cp: a pseudogene of C
    fams = P.families(prot, near, level2={"S1": ["B", "B2"], "S2": ["C", "Cp"]})
    assert fams == [dict(genes=["A", "B", "B2"], subfamilies=["S1"]), dict(genes=["C", "Cp"], subfamilies=["S2"])]


def test_the_member_check_reports_each_known_family_at_both_levels_and_every_miss():
    genes = {"g1": "DAZ1", "g2": "DAZ2", "g3": "DAZ3", "g4": "DAZ4", "x1": "OTHER1"}
    l1 = [dict(genes=["g1", "g2", "g3", "x1"], subfamilies=["S1"])]
    l2 = {"S1": ["g2", "g3"]}
    rep = P.member_report({"DAZ": r"^DAZ[1-4]$"}, genes, l1, l2)
    d = rep["DAZ"]
    assert d["members"] == 4 and d["level1_families"] == 1 and d["level1_largest"] == 3
    assert d["level1_missing"] == ["DAZ4"] and d["level1_other_genes"] == 1
    assert d["level2_subfamilies"] == 1 and d["level2_largest"] == 2 and d["level2_missing"] == ["DAZ1", "DAZ4"]


def test_a_build_can_reuse_the_protein_links_of_an_earlier_one():
    assert P.build_args(["out.json", "pool.json"]) == ("out.json", "pool.json", 4, None)
    assert P.build_args(["out.json", "pool.json", "--threads", "2", "--protein-links", "l.tsv"]) == ("out.json", "pool.json", 2, "l.tsv")
