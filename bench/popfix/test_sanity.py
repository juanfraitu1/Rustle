"""Tests for sanity.py, arm (a) of docs/archive/2026-10/merged/PREREG_popfix_2026-10-12.md (run: python3 -m pytest -q bench/popfix)."""
import io, json, os, sys, contextlib

sys.path.insert(0, os.path.dirname(__file__))
import sanity as Z
from runfixture import write_run, write_pop_rows


def _pop(path):
    # (chrom, pos, family_id, status, method)
    write_pop_rows(path, [
        ("chrY", 10, "F1", "FIXED", "K"), ("chrY", 11, "F1", "FIXED", "A"), ("chrY", 12, "F1", "POLYMORPHIC", "K"),
        ("chrY", 13, "F1", "UNCOVERED", "K"), ("chrY", 14, "F1", "NOT_ASSESSED", "-"),
        ("chrY", 20, "F2", "FIXED", "K"), ("chr1", 30, "F9", "NOT_ASSESSED", "-"),
    ])


def test_status_method_counts_per_yag_family_only(tmp_path):
    pop = str(tmp_path / "pop.tsv")
    _pop(pop)
    res = Z.status_method_counts(pop, ["F1", "F2", "F3"])
    assert res["F1"]["n_sites"] == 5
    assert res["F1"]["status_x_method"] == {"FIXED": {"K": 1, "A": 1}, "POLYMORPHIC": {"K": 1}, "UNCOVERED": {"K": 1},
                                            "NOT_ASSESSED": {"-": 1}}
    assert res["F2"]["status_x_method"] == {"FIXED": {"K": 1}}
    assert res["F3"] == {"n_sites": 0, "status_x_method": {}}   # a YAG family with no rows is listed, with zeros
    assert "F9" not in res                                       # a non-YAG family is not


def _star(name, cols, cand=(0, 1), contested=1, fam="F1"):
    return dict(read=name, fam=fam, status="tied", copy=0, cand=list(cand), cols=cols, contested=contested)


def test_column_fractions_count_contested_molecules_only(tmp_path):
    pop = str(tmp_path / "pop.tsv")
    _pop(pop)
    cols = [
        (1, "A", "AC", [10, 11]),      # both positions FIXED
        (2, "A", "AC", [10, 12]),      # both in the table, one POLYMORPHIC: ANY yes, FIXED no
        (3, "A", "AC", [10, 99]),      # 99 is not a catalog site: neither
        (4, "A", "AC", [None, None]),  # no position at all: not a numerator, counted as no-position
        (5, "A", "AC", [None, 11]),    # one positioned candidate and it is FIXED: counts for both
    ]
    stars = [_star("c1", cols), _star("u1", cols, contested=0), _star("other", cols, fam="F9")]
    p, _ = write_run(str(tmp_path), stars, joins=[("F1", 0, "chrY"), ("F1", 1, "chrY"), ("F9", 0, "chrY"), ("F9", 1, "chrY")])
    res = Z.column_site_fractions(p, pop, ["F1"])
    f = res["F1"]
    assert (f["n_molecules"], f["n_cols"], f["n_cols_no_position"]) == (1, 5, 1)
    assert (f["n_cols_all_any"], f["n_cols_all_fixed"]) == (3, 2)
    assert f["frac_any"] == 0.6 and f["frac_fixed"] == 0.4
    assert f["frac_any_of_positioned"] == 0.75 and f["frac_fixed_of_positioned"] == 0.5


def test_column_fractions_use_each_candidates_own_chrom(tmp_path):
    pop = str(tmp_path / "pop.tsv")
    _pop(pop)
    cols = [(1, "A", "AC", [10, 10])]
    p, _ = write_run(str(tmp_path), [_star("c1", cols)], joins=[("F1", 0, "chrY"), ("F1", 1, "chrX")])
    f = Z.column_site_fractions(p, pop, ["F1"])["F1"]   # (chrX, 10) is not a site: candidate 1 is outside both sets
    assert (f["n_cols_all_any"], f["n_cols_all_fixed"]) == (0, 0)


def test_main_writes_json_and_works_without_a_run(tmp_path):
    pop = str(tmp_path / "pop.tsv")
    _pop(pop)
    fam = tmp_path / "yag.json"
    fam.write_text(json.dumps({"families": {"F1": ["A"], "F2": ["B"]}}))
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        Z.main(pop, str(fam))
    out = json.loads(buf.getvalue())
    assert sorted(out["status_x_method"]) == ["F1", "F2"] and "column_fractions" not in out
    p, _ = write_run(str(tmp_path), [_star("c1", [(1, "A", "AC", [10, 11])])], joins=[("F1", 0, "chrY"), ("F1", 1, "chrY")])
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        Z.main(pop, str(fam), p)
    out = json.loads(buf.getvalue())
    assert out["column_fractions"]["F1"]["frac_fixed"] == 1.0
