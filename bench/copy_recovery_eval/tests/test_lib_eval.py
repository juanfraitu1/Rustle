import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import lib_eval


def test_build_universe_groups_paralog_families():
    gene_tx = {
        "RABL2A": ["rna-A1", "rna-A2"],
        "RABL2B": ["rna-B1"],
        "SOLO":   ["rna-S1"],
    }
    paralog_pairs = [("RABL2A", "RABL2B")]   # SOLO has no paralog
    uni = lib_eval.build_universe(gene_tx, paralog_pairs)
    tx_in = {r["transcript_id"] for r in uni}
    assert tx_in == {"rna-A1", "rna-A2", "rna-B1"}
    fam = {r["transcript_id"]: r["family_id"] for r in uni}
    assert fam["rna-A1"] == fam["rna-A2"] == fam["rna-B1"]
    assert "rna-S1" not in tx_in
    assert all(r["n_family_copies"] == 2 for r in uni)


def test_build_universe_transitive_closure():
    gene_tx = {"A": ["a"], "B": ["b"], "C": ["c"], "D": ["d"]}
    paralog_pairs = [("A", "B"), ("B", "C")]
    uni = lib_eval.build_universe(gene_tx, paralog_pairs)
    fam = {r["transcript_id"]: r["family_id"] for r in uni}
    assert fam["a"] == fam["b"] == fam["c"]
    assert "d" not in fam
    assert all(r["n_family_copies"] == 3 for r in uni)
