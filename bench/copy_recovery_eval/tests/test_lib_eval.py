import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import lib_eval


def test_recovery_set_joins_across_rna_prefix():
    # SQANTI3 associated_transcript carries the gffread/GFF3 'rna-' prefix; the
    # universe uses the clean accession. They must still join.
    classif = [
        {"isoform": "q1", "structural_category": "full-splice_match", "associated_transcript": "rna-XM_055380435.2"},
        {"isoform": "q2", "structural_category": "full-splice_match", "associated_transcript": "rna-XR_999.1"},  # not in U
    ]
    universe_tx = {"XM_055380435.2", "XM_031015533.3"}  # clean accessions
    rec = lib_eval.recovery_set(classif, universe_tx)
    assert rec == {"XM_055380435.2": {"fsm": True, "ism": False}}


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
    # Clean (un-prefixed) ids isolate the FSM/ISM logic from prefix normalization
    # (which has its own test, test_recovery_set_joins_across_rna_prefix).
    classif = [
        {"isoform": "q1", "structural_category": "full-splice_match",       "associated_transcript": "A1"},
        {"isoform": "q2", "structural_category": "incomplete-splice_match", "associated_transcript": "A2"},
        {"isoform": "q3", "structural_category": "novel_not_in_catalog",    "associated_transcript": "novel"},
        {"isoform": "q4", "structural_category": "full-splice_match",       "associated_transcript": "S1"},  # not in U
    ]
    universe_tx = {"A1", "A2", "B1"}
    rec = lib_eval.recovery_set(classif, universe_tx)
    assert rec["A1"] == {"fsm": True, "ism": False}
    assert rec["A2"] == {"fsm": False, "ism": True}
    assert "S1" not in rec
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


def test_classify_authenticity_three_buckets():
    # decisive own-copy reads -> authentic
    assert lib_eval.classify_authenticity(5, 0, k_decisive=2) == lib_eval.AUTHENTIC
    assert lib_eval.classify_authenticity(2, 0, k_decisive=2) == lib_eval.AUTHENTIC
    # not enough decisive but tied reads exist -> unresolvable (NOT phantom)
    assert lib_eval.classify_authenticity(1, 30, k_decisive=2) == lib_eval.UNRESOLVABLE
    assert lib_eval.classify_authenticity(0, 1, k_decisive=2) == lib_eval.UNRESOLVABLE
    # no own evidence and no ties (sister spillover only) -> phantom
    assert lib_eval.classify_authenticity(0, 0, k_decisive=2) == lib_eval.PHANTOM
    assert lib_eval.classify_authenticity(1, 0, k_decisive=2) == lib_eval.PHANTOM


def test_head_to_head_three_way_status():
    rustle = {
        "A1": {"fsm": True, "ism": False},   # authentic -> win
        "B1": {"fsm": True, "ism": False},   # phantom
        "U1": {"fsm": True, "ism": False},   # unresolvable
    }
    stringtie = {}
    status = {"A1": lib_eval.AUTHENTIC, "B1": lib_eval.PHANTOM, "U1": lib_eval.UNRESOLVABLE}
    fam = {"A1": "FAM1", "B1": "FAM2", "U1": "FAM3"}
    res = lib_eval.head_to_head(rustle, stringtie, status, fam)
    assert res["rustle_only_fsm_authentic"] == ["A1"]
    assert res["rustle_only_fsm_phantom"] == ["B1"]
    assert res["rustle_only_fsm_unresolvable"] == ["U1"]
    assert (res["n_win"], res["n_phantom"], res["n_unresolvable"]) == (1, 1, 1)
    assert res["families_won"] == ["FAM1"]
