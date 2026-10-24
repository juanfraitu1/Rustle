"""Tests for check_e.py, arm (e) of docs/archive/2026-10/merged/PREREG_popfix_2026-10-12.md (run: python3 -m pytest -q bench/popfix)."""
import contextlib, io, json, os, sys

sys.path.insert(0, os.path.dirname(__file__))
import check_e as E
from runfixture import write_run, write_pop_rows

JOINS = [(f, c, "chrY") for f in ("F1", "F9") for c in (0, 1)]


def _pop(path):
    write_pop_rows(path, [("chrY", p, "F1", "FIXED", "K") for p in (100, 101, 200, 201)] +
                         [("chrY", 300, "F1", "POLYMORPHIC", "K"), ("chrY", 400, "F1", "NOT_ASSESSED", "-")])


FIXED_COLS = [(10, "A", "AC", [100, 200]), (20, "A", "AC", [101, 201])]


def _mol(name, cols=FIXED_COLS, fam="F1", contested=1, **kw):
    return dict(read=name, fam=fam, status="assigned", copy=0, cand=[0, 1], cols=cols, contested=contested, **kw)


def _tier(name, all_=("assigned", 0), fixed=("assigned", 0), fam="F1"):
    return (name, fam, all_[0], all_[1], fixed[0], fixed[1], 2, "PORTABLE")


def _run(tmp_path, stars, tiers, env=None, with_err=False):
    pop = str(tmp_path / "pop.tsv")
    _pop(pop)
    p, _ = write_run(str(tmp_path), stars, tiers=tiers, joins=JOINS)
    saved = os.environ.pop("YAG_FAMILIES", None)
    os.environ.update(env or {})
    buf, err = io.StringIO(), io.StringIO()
    try:
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(err):
            code = E.main(p, pop)
    finally:
        os.environ.pop("YAG_FAMILIES", None)
        if saved is not None:
            os.environ["YAG_FAMILIES"] = saved
    return (code, json.loads(buf.getvalue()), err.getvalue()) if with_err else (code, json.loads(buf.getvalue()))


def test_agreeing_eligible_molecules_pass(tmp_path):
    stars = [_mol("a"), _mol("t")]
    tiers = [_tier("a"), _tier("t", ("tied", 0), ("tied", 0))]
    code, out = _run(tmp_path, stars, tiers)
    assert code == 0
    assert (out["n_eligible"], out["n_assigned_any"], out["n_mismatch"]) == (2, 1, 0)


def test_status_and_copy_mismatches_exit_one_and_are_listed(tmp_path):
    stars = [_mol("copy"), _mol("all_only"), _mol("fixed_only"), _mol("ok")]
    tiers = [_tier("copy", ("assigned", 0), ("assigned", 1)),
             _tier("all_only", ("assigned", 0), ("ambiguous", 0)),
             _tier("fixed_only", ("tied", 0), ("assigned", 0)),
             _tier("ok")]
    code, out = _run(tmp_path, stars, tiers)
    assert code == 1
    assert (out["n_eligible"], out["n_assigned_any"], out["n_mismatch"]) == (4, 4, 3)
    kinds = {m["read_name"]: m["kind"] for m in out["first_mismatches"]}
    assert kinds == {"copy": "both_assigned_different_copy", "all_only": "only_status_all_assigned",
                     "fixed_only": "only_status_fixed_assigned"}


def test_first_ten_mismatches_only(tmp_path):
    stars = [_mol(f"m{i}") for i in range(12)]
    tiers = [_tier(f"m{i}", ("assigned", 0), ("assigned", 1)) for i in range(12)]
    code, out = _run(tmp_path, stars, tiers)
    assert code == 1 and out["n_mismatch"] == 12 and len(out["first_mismatches"]) == 10


def test_a_molecule_is_eligible_only_if_every_column_is_fixed_for_every_candidate(tmp_path):
    poly = FIXED_COLS + [(30, "A", "AC", [300, 200])]          # 300 is POLYMORPHIC
    nopos = FIXED_COLS + [(30, "A", "AC", [None, 200])]         # a candidate with no position
    other_chrom = FIXED_COLS + [(30, "A", "AC", [100, 201])]    # fine: still FIXED (control)
    stars = [_mol("poly", poly), _mol("nopos", nopos), _mol("ctl", other_chrom)]
    bad = ("assigned", 0), ("ambiguous", 0)
    tiers = [_tier("poly", *bad), _tier("nopos", *bad), _tier("ctl")]
    code, out = _run(tmp_path, stars, tiers)
    assert code == 0 and out["n_eligible"] == 1                 # only the control counts


def test_chrom_of_each_candidate_matters(tmp_path):
    pop = str(tmp_path / "pop.tsv")
    _pop(pop)
    p, _ = write_run(str(tmp_path), [_mol("x")], tiers=[_tier("x", ("assigned", 0), ("ambiguous", 0))],
                     joins=[("F1", 0, "chrY"), ("F1", 1, "chrX")])   # (chrX, 200) is not a FIXED site
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        code = E.main(p, pop)
    assert code == 0 and json.loads(buf.getvalue())["n_eligible"] == 0


def test_later_gates_uncontested_no_tier_zero_columns_and_non_yag_are_excluded(tmp_path):
    bad = ("assigned", 0), ("ambiguous", 0)
    stars = [_mol("origin", origin_rejected=1), _mol("tie_out", tie_outside=1), _mol("unc", contested=0),
             _mol("notier"), _mol("nocols", cols=[]), _mol("other_fam", fam="F9"), _mol("ok")]
    tiers = [_tier("origin", *bad), _tier("tie_out", *bad), _tier("unc", *bad),
             ("notier", "F1", "tied", 0, "NA", "NA", 0, "UNASSIGNED"), _tier("nocols", ("tied", 0), ("tied", 0)),
             _tier("other_fam", *bad, fam="F9"), _tier("ok")]
    fam = tmp_path / "yag.json"
    fam.write_text(json.dumps({"families": {"F1": ["A"]}}))
    code, out = _run(tmp_path, stars, tiers, {"YAG_FAMILIES": str(fam)})
    assert code == 0
    assert (out["n_eligible"], out["n_assigned_any"], out["n_mismatch"], out["n_eligible_no_columns"]) == (1, 1, 0, 1)
    # without the YAG filter the F9 molecule is read too (S is positional, not per family), and its disagreement counts
    code, out = _run(tmp_path, stars, tiers)
    assert code == 1 and [m["read_name"] for m in out["first_mismatches"]] == ["other_fam"]


def test_no_assigned_molecule_is_a_vacuous_pass_flagged_in_output_and_on_stderr(tmp_path):
    stars = [_mol("t1"), _mol("t2")]
    tiers = [_tier("t1", ("tied", 0), ("tied", 0)), _tier("t2", ("ambiguous", 0), ("ambiguous", 1))]
    code, out, err = _run(tmp_path, stars, tiers, with_err=True)
    assert code == 0                                    # exit code unchanged: not a mismatch
    assert (out["n_eligible"], out["n_assigned_any"], out["n_mismatch"]) == (2, 0, 0)
    assert out["vacuous"] is True
    assert "vacuous" in err.lower() and "n_assigned_any" in err


def test_an_assigned_molecule_makes_the_check_non_vacuous(tmp_path):
    code, out, err = _run(tmp_path, [_mol("a")], [_tier("a")], with_err=True)
    assert code == 0 and out["vacuous"] is False and err == ""
