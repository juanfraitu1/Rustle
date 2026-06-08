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


def test_recovery_set_fsm_ism_within_universe():
    classif = [
        {"isoform": "q1", "structural_category": "full-splice_match",       "associated_transcript": "rna-A1"},
        {"isoform": "q2", "structural_category": "incomplete-splice_match", "associated_transcript": "rna-A2"},
        {"isoform": "q3", "structural_category": "novel_not_in_catalog",    "associated_transcript": "novel"},
        {"isoform": "q4", "structural_category": "full-splice_match",       "associated_transcript": "rna-S1"},  # not in U
    ]
    universe_tx = {"rna-A1", "rna-A2", "rna-B1"}
    rec = lib_eval.recovery_set(classif, universe_tx)
    assert rec["rna-A1"] == {"fsm": True, "ism": False}
    assert rec["rna-A2"] == {"fsm": False, "ism": True}
    assert "rna-S1" not in rec
    assert "novel" not in rec


def test_head_to_head_splits_authentic_and_phantom():
    rustle = {
        "rna-A1": {"fsm": True, "ism": False},   # ST misses, authentic -> WIN
        "rna-B1": {"fsm": True, "ism": False},   # ST misses, phantom   -> phantom
        "rna-C1": {"fsm": True, "ism": False},   # ST also has it       -> not rustle-only
    }
    stringtie = {
        "rna-C1": {"fsm": True, "ism": False},
    }
    authentic = {"rna-A1": True, "rna-B1": False, "rna-C1": True}
    fam = {"rna-A1": "FAM1", "rna-B1": "FAM2", "rna-C1": "FAM3"}
    res = lib_eval.head_to_head(rustle, stringtie, authentic, fam)
    assert res["rustle_only_fsm_authentic"] == ["rna-A1"]
    assert res["rustle_only_fsm_phantom"] == ["rna-B1"]
    assert res["n_win"] == 1
    assert res["n_phantom"] == 1
