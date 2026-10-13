"""Tests for split_half.py (run: python3 -m pytest -q bench/yag_o2).

The population mask is PER PAIR (ruling R27, mirroring Rust `fixed_only_certificate`): a column counts for the pair
(a, b) iff both candidates have a genome position there and both (chrom, pos) are in the site set S; candidate a's
PSV-score tie key counts only columns where a's own position is in S.
"""
import io, json, os, random, sys
import contextlib

sys.path.insert(0, os.path.dirname(__file__))
import split_half as S

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "popfix"))
from runfixture import write_run, write_pop  # noqa: E402  (synthetic run directories: star dump, assignments, tiers)


def run_main(prefix, lengths, out, env=None):
    saved = {k: os.environ.get(k) for k in ("POPULATION", "POPULATION_ANY", "YAG_FAMILIES")}
    for k in saved:
        os.environ.pop(k, None)
    os.environ.update(env or {})
    buf, err = io.StringIO(), io.StringIO()
    try:
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(err):
            S.main(prefix, lengths, out)
    finally:
        for k, v in saved.items():
            os.environ.pop(k, None)
            if v is not None:
                os.environ[k] = v
    return buf.getvalue(), err.getvalue()


# ---- the certificate under the per-pair mask -----------------------------------------------------------------------

SITES = {("chrY", 1000), ("chrY", 2000), ("chrY", 1010), ("chrY", 2010)}


def test_plan_case_per_pair_one_fixed_column_never_assigns():
    # the plan's case: column 1 is fixed for the pair (both positions in S), column 2 is not (1010 / 2010 absent)
    sites = {("chrY", 1000), ("chrY", 2000)}
    cols = [(10, "A", "AC"), (20, "G", "GT")]
    gpos = [[1000, 2000], [1010, 2010]]
    mcols = S.mask_columns(cols, gpos, ["chrY", "chrY"], sites)
    st, cp = S.certify([0, 1], mcols)
    assert st == "ambiguous"  # one fixed column (m = 1): p = alpha exactly, never assigns
    # both columns fixed -> two decisive columns -> assigned: the mask, not the data, made the difference
    mcols2 = S.mask_columns(cols, gpos, ["chrY", "chrY"], sites | {("chrY", 1010), ("chrY", 2010)})
    assert S.certify([0, 1], mcols2) == ("assigned", 0)
    # no mask at all (3-tuple columns): the original certificate
    assert S.certify([0, 1], cols) == ("assigned", 0)


def _three_copy_cols(z_b_in_sites):
    """3 candidates a=0, b=1, c=2; the read matches a. Columns 1-2: a vs b decisive, all fixed. Columns 3-4 ('Z'): the
    read matches b, a=A b=C c=G, a's and c's positions fixed and b's fixed only when z_b_in_sites. Columns 5-6: a vs c
    decisive, all fixed (b has a's allele)."""
    cols, gpos, sites = [], [], set()
    def add(rp, obs, al, positions, in_sites):
        cols.append((rp, obs, al))
        gpos.append(positions)
        for k, ok in zip(positions, in_sites):
            if ok:
                sites.add(("chrY", k))
    add(10, "A", "ACA", [100, 200, 300], [1, 1, 1])
    add(20, "A", "ACA", [101, 201, 301], [1, 1, 1])
    add(30, "C", "ACG", [102, 202, 302], [1, int(z_b_in_sites), 1])
    add(40, "C", "ACG", [103, 203, 303], [1, int(z_b_in_sites), 1])
    add(50, "A", "AAG", [104, 204, 304], [1, 1, 1])
    add(60, "A", "AAG", [105, 205, 305], [1, 1, 1])
    return cols, gpos, sites


def test_column_fixed_for_a_c_but_not_b_does_not_count_for_a_b():
    cols, gpos, sites = _three_copy_cols(z_b_in_sites=False)
    mcols = S.mask_columns(cols, gpos, ["chrY"] * 3, sites)
    # Z columns (read matches b) are dropped from (a, b) but kept for (a, c): a vs b = 2:0, a vs c = 2:0 -> assigned a
    assert S.certify([0, 1, 2], mcols) == ("assigned", 0)
    # make b's Z positions fixed: (a, b) now sees the two Z columns that favour b -> 2:2, margin 0 -> not assigned
    cols, gpos, sites = _three_copy_cols(z_b_in_sites=True)
    mcols = S.mask_columns(cols, gpos, ["chrY"] * 3, sites)
    assert S.certify([0, 1, 2], mcols)[0] == "ambiguous"


def test_per_pair_mask_is_not_a_per_column_mask():
    # 3 candidates a=0, b=1, c=2, the read matches a. X columns: a vs b decisive, a's and b's positions in S, c's not.
    # Z columns: a vs c decisive, a's and c's positions in S, b's not. No column is fixed for ALL three candidates, so
    # a per-column rule has no evidence at all (tied); the per-pair rule sees (a, b) at X and (a, c) at Z -> assigned.
    cols = [(10, "A", "ACA"), (20, "A", "ACA"), (30, "A", "AAG"), (40, "A", "AAG")]
    gpos = [[100, 200, 300], [101, 201, 301], [102, 202, 302], [103, 203, 303]]
    sites = {("chrY", 100), ("chrY", 200), ("chrY", 101), ("chrY", 201),   # X: a, b
             ("chrY", 102), ("chrY", 302), ("chrY", 103), ("chrY", 303)}   # Z: a, c
    mcols = S.mask_columns(cols, gpos, ["chrY"] * 3, sites)
    assert [c[3] for c in mcols] == [[1, 1, 0], [1, 1, 0], [1, 0, 1], [1, 0, 1]]
    assert S.certify([0, 1, 2], mcols) == ("assigned", 0)


def test_no_fixed_column_for_a_pair_is_a_tie_not_a_skip():
    cols, gpos, sites = _three_copy_cols(z_b_in_sites=False)
    # drop c's positions from S everywhere: (a, c) shares no fixed column -> Tied
    sites = {s for s in sites if s[1] < 300 or s[1] >= 400}
    mcols = S.mask_columns(cols, gpos, ["chrY"] * 3, sites)
    assert S.certify([0, 1, 2], mcols)[0] == "tied"


def test_chrom_matters_a_chrX_site_does_not_count_for_a_chrY_candidate():
    sites = {("chrX", 1000), ("chrY", 2000)}
    cols = [(10, "A", "AC"), (20, "A", "AC")]
    gpos = [[1000, 2000], [1000, 2000]]
    # candidate 0 on chrX: (chrX, 1000) and (chrY, 2000) are both in S -> both columns count -> assigned
    assert S.certify([0, 1], S.mask_columns(cols, gpos, ["chrX", "chrY"], sites)) == ("assigned", 0)
    # candidate 0 on chrY: (chrY, 1000) is NOT in S (only chrX's is) -> the columns never count -> tied
    assert S.certify([0, 1], S.mask_columns(cols, gpos, ["chrY", "chrY"], sites))[0] == "tied"


def test_candidate_without_a_position_is_never_in_s():
    sites = {("chrY", 2000)}
    cols = [(10, "A", "AC"), (20, "A", "AC")]
    gpos = [[None, 2000], [None, 2000]]
    assert S.certify([0, 1], S.mask_columns(cols, gpos, ["chrY", "chrY"], sites))[0] == "tied"


def test_missing_positions_mean_no_fixed_column():
    # R20: a dump with no positions for the columns (off the two-form path) -> tied, never the full-column call
    cols = [(10, "A", "AC"), (20, "A", "AC")]
    assert S.certify([0, 1], S.mask_columns(cols, None, ["chrY", "chrY"], {("chrY", 1)}))[0] == "tied"
    assert S.certify([0, 1], cols) == ("assigned", 0)


def test_score_tie_key_counts_only_columns_where_the_candidates_own_position_is_in_s():
    # the pair (0, 1) has no column in S for both, so the LLR key ties and the PSV score decides. Only candidate 1's
    # positions are in S: its score is +2 (it matches the read there) against 0 for candidate 0 -> best = 1. Were the
    # PAIR mask applied to the score, both scores would be 0 and candidate 0 (lower index) would win.
    cols = [(10, "A", "TA"), (20, "A", "TA")]
    gpos = [[500, 600], [501, 601]]
    sites = {("chrY", 600), ("chrY", 601)}
    mcols = S.mask_columns(cols, gpos, ["chrY", "chrY"], sites)
    assert [c[3] for c in mcols] == [[False, True], [False, True]]
    st, cp = S.certify([0, 1], mcols)
    assert (st, cp) == ("tied", 1)


# ---- site sets ------------------------------------------------------------------------------------------------------

def test_load_sites_fixed_only_vs_any_status(tmp_path):
    pop = tmp_path / "pop.tsv"
    write_pop(str(pop), [("chrY", 10, "FIXED"), ("chrY", 11, "POLYMORPHIC"), ("chrY", 12, "NOT_ASSESSED"),
                         ("chrY", 13, "UNCOVERED"), ("chr1", 10, "NOT_ASSESSED")])
    assert S.load_sites(str(pop)) == {("chrY", 10)}
    assert S.load_sites(str(pop), any_status=True) == {("chrY", 10), ("chrY", 11), ("chrY", 12), ("chrY", 13), ("chr1", 10)}


def test_population_any_admits_not_assessed_and_polymorphic_sites_end_to_end(tmp_path):
    pop = tmp_path / "pop.tsv"
    write_pop(str(pop), [("chrY", 100, "NOT_ASSESSED"), ("chrY", 200, "POLYMORPHIC"),
                         ("chrY", 101, "NOT_ASSESSED"), ("chrY", 201, "POLYMORPHIC")])
    cols = [(10, "A", "AC", [100, 200]), (20, "A", "AC", [101, 201])]
    star = dict(read="r1", fam="F", status="assigned", copy=0, cand=[0, 1], cols=cols)
    p, lp = write_run(str(tmp_path), [star], joins=[("F", 0, "chrY"), ("F", 1, "chrY")], lengths={"r1": 100})
    out = str(tmp_path / "o.json")
    run_main(p, lp, out, {"POPULATION": str(pop)})
    fixed_only = json.load(open(out))
    run_main(p, lp, out, {"POPULATION": str(pop), "POPULATION_ANY": "1"})
    any_mode = json.load(open(out))
    # FIXED-only: no site is FIXED -> tied (no half assigns); ANY: both sites present -> the certificate sees 2 columns
    assert fixed_only["counts"].get("some_half_assigns", 0) == 0
    assert any_mode["population"]["mode"] == "any"
    assert any_mode["population"]["full_read_status"]["assigned"] == 1
    assert fixed_only["population"]["mode"] == "fixed"
    assert fixed_only["population"]["full_read_status"]["tied"] == 1


# ---- main(): contested rows only, P0-fixed ---------------------------------------------------------------------------

def _mol(name, status, copy, cols, contested=1, cand=(0, 1), fam="F", **kw):
    return dict(read=name, fam=fam, status=status, copy=copy, cand=list(cand), cols=cols, contested=contested, **kw)


def _fixed_cols(n=4, off=0):
    return [(10 * (i + 1), "A", "AC", [100 + off + i, 200 + off + i]) for i in range(n)]


def _fixed_sites(path, n=4, off=0):
    rows = []
    for i in range(n):
        rows += [("chrY", 100 + off + i, "FIXED"), ("chrY", 200 + off + i, "FIXED")]
    write_pop(path, rows)


def test_p0_fixed_reproduces_tier_and_reports_mismatches(tmp_path):
    pop = str(tmp_path / "pop.tsv")
    _fixed_sites(pop)
    stars = [_mol("good", "assigned", 0, _fixed_cols()),
             _mol("tie_row", "tied", 0, [(10, "A", "AA", [100, 200])]),
             _mol("bad", "assigned", 0, _fixed_cols())]
    tiers = [("good", "F", "assigned", 0, "assigned", 0, 4, "PORTABLE"),
             ("tie_row", "F", "tied", 0, "tied", 0, 0, "UNASSIGNED"),
             ("bad", "F", "assigned", 0, "assigned", 1, 4, "PORTABLE"),   # wrong copy -> mismatch
             ("na_row", "F", "tied", 0, "NA", "NA", 0, "UNASSIGNED")]     # no tier: excluded from the gate
    p, lp = write_run(str(tmp_path), stars, tiers=tiers, joins=[("F", 0, "chrY"), ("F", 1, "chrY")],
                      lengths={"good": 100, "tie_row": 100, "bad": 100})
    out = str(tmp_path / "o.json")
    run_main(p, lp, out, {"POPULATION": pop})
    res = json.load(open(out))
    p0 = res["population"]["P0_fixed"]
    assert (p0["n_tier_rows"], p0["n_compared"], p0["n_same"]) == (3, 3, 2)
    assert res["P0_fixed_reproduction"] == round(2 / 3, 4)
    assert res["P0_fixed_pass"] is False
    assert [m["read_name"] for m in p0["first_mismatches"]] == ["bad"]
    assert p0["first_mismatches"][0]["tier"] == ["assigned", "1"]
    assert p0["first_mismatches"][0]["python"] == ["assigned", "0"]


def test_p0_fixed_counts_a_tier_row_without_a_star_row_as_a_failure(tmp_path):
    pop = str(tmp_path / "pop.tsv")
    _fixed_sites(pop)
    stars = [_mol("good", "assigned", 0, _fixed_cols())]
    tiers = [("good", "F", "assigned", 0, "assigned", 0, 4, "PORTABLE"),
             ("ghost", "F", "assigned", 0, "assigned", 0, 4, "PORTABLE")]
    p, lp = write_run(str(tmp_path), stars, tiers=tiers, joins=[("F", 0, "chrY"), ("F", 1, "chrY")],
                      lengths={"good": 100})
    out = str(tmp_path / "o.json")
    run_main(p, lp, out, {"POPULATION": pop})
    res = json.load(open(out))
    p0 = res["population"]["P0_fixed"]
    assert (p0["n_tier_rows"], p0["n_compared"], p0["n_no_star_row"], p0["n_same"]) == (2, 1, 1, 1)
    assert res["P0_fixed_reproduction"] == 0.5


def test_population_mode_reads_contested_rows_only(tmp_path):
    pop = str(tmp_path / "pop.tsv")
    _fixed_sites(pop)
    stars = [_mol("c1", "assigned", 0, _fixed_cols(), contested=1), _mol("u1", "assigned", 0, _fixed_cols(), contested=0)]
    tiers = [("c1", "F", "assigned", 0, "assigned", 0, 4, "PORTABLE")]
    p, lp = write_run(str(tmp_path), stars, tiers=tiers, joins=[("F", 0, "chrY"), ("F", 1, "chrY")],
                      lengths={"c1": 100, "u1": 100})
    out = str(tmp_path / "o.json")
    run_main(p, lp, out, {"POPULATION": pop})
    res = json.load(open(out))
    assert res["counts"]["uncontested_skipped"] == 1
    assert res["per_family"]["F"]["n"] == 1
    assert res["per_family"]["F"]["full|assigned"] == 1


def test_split_half_under_the_mask(tmp_path):
    # 4 fixed columns in the 5' half, 4 in the 3' half; read length 1000 -> half = 500
    pop = str(tmp_path / "pop.tsv")
    rows = []
    cols = []
    for i in range(4):
        cols.append((10 + i, "A", "AC", [100 + i, 200 + i]))
        cols.append((900 + i, "A", "AC", [150 + i, 250 + i]))
        rows += [("chrY", 100 + i, "FIXED"), ("chrY", 200 + i, "FIXED"), ("chrY", 150 + i, "FIXED"), ("chrY", 250 + i, "FIXED")]
    write_pop(pop, rows)
    stars = [_mol("m", "assigned", 0, cols)]
    p, lp = write_run(str(tmp_path), stars, joins=[("F", 0, "chrY"), ("F", 1, "chrY")], lengths={"m": 1000})
    out = str(tmp_path / "o.json")
    run_main(p, lp, out, {"POPULATION": pop})
    res = json.load(open(out))
    assert res["counts"]["both_halves_assign"] == 1 and res["counts"]["both_halves_agree"] == 1
    # drop the 3' half's sites from S: only the 5' half assigns
    write_pop(pop, [r for r in rows if r[1] < 150 or 200 <= r[1] < 250])
    run_main(p, lp, out, {"POPULATION": pop})
    res = json.load(open(out))
    assert res["counts"].get("both_halves_assign", 0) == 0 and res["counts"]["some_half_assigns"] == 1


def test_population_mode_needs_a_positions_column_and_a_family_join(tmp_path):
    pop = str(tmp_path / "pop.tsv")
    _fixed_sites(pop)
    stars = [_mol("m", "assigned", 0, _fixed_cols())]
    p, lp = write_run(str(tmp_path), stars, joins=[("F", 0, "chrY"), ("F", 1, "chrY")], positions=False, lengths={"m": 100})
    try:
        run_main(p, lp, str(tmp_path / "o.json"), {"POPULATION": pop})
        raise AssertionError("expected SystemExit for a dump without positions")
    except SystemExit as e:
        assert "positions" in str(e)
    p, lp = write_run(str(tmp_path), stars, joins=None, lengths={"m": 100}, prefix="nojoin")
    try:
        run_main(p, lp, str(tmp_path / "o.json"), {"POPULATION": pop})
        raise AssertionError("expected SystemExit for a run without family_join")
    except SystemExit as e:
        assert "family_join" in str(e)


def test_population_any_without_population_is_an_error(tmp_path):
    stars = [_mol("m", "assigned", 0, _fixed_cols())]
    p, lp = write_run(str(tmp_path), stars, joins=[("F", 0, "chrY"), ("F", 1, "chrY")], lengths={"m": 100})
    out = str(tmp_path / "a.json")
    try:
        run_main(p, lp, out, {"POPULATION_ANY": "1"})
        raise AssertionError("expected SystemExit: POPULATION_ANY=1 needs POPULATION")
    except SystemExit as e:
        assert "POPULATION_ANY" in str(e) and "POPULATION" in str(e).replace("POPULATION_ANY", "")
    assert not os.path.exists(out)              # nothing written: a silent FIXED/no-mask run would be the wrong arm
    # an unrecognised value is an error too (a stray "yes" must not silently select the FIXED set)
    pop = str(tmp_path / "pop.tsv")
    _fixed_sites(pop)
    try:
        run_main(p, lp, out, {"POPULATION": pop, "POPULATION_ANY": "yes"})
        raise AssertionError("expected SystemExit for POPULATION_ANY=yes")
    except SystemExit as e:
        assert "POPULATION_ANY" in str(e)
    # 0 / empty mean off: the plain no-POPULATION run
    run_main(p, lp, out, {"POPULATION_ANY": "0"})
    assert os.path.exists(out)


# ---- the prereg's denominator fields -----------------------------------------------------------------------------------

def _two_half_cols(a_half, b_half):
    """5' half columns (read pos < 500) and 3' half columns, 4 each if present, all FIXED for the pair."""
    cols = []
    for i in range(4):
        if a_half:
            cols.append((10 + i, "A", "AC", [100 + i, 200 + i]))
        if b_half:
            cols.append((900 + i, "A", "AC", [150 + i, 250 + i]))
    return cols


def _all_fixed_sites(path):
    rows = [("chrY", 100 + i, "FIXED") for i in range(4)] + [("chrY", 200 + i, "FIXED") for i in range(4)] + \
           [("chrY", 150 + i, "FIXED") for i in range(4)] + [("chrY", 250 + i, "FIXED") for i in range(4)]
    write_pop(path, rows)


def test_denominator_fields_population_block_and_per_family(tmp_path):
    pop = str(tmp_path / "pop.tsv")
    _all_fixed_sites(pop)
    stars = [
        _mol("both", "assigned", 0, _two_half_cols(True, True)),      # both halves assign
        _mol("one", "assigned", 0, _two_half_cols(True, False)),      # only the 5' half assigns
        _mol("none", "tied", 0, [(10, "A", "AC", [100, 200])]),        # no half assigns (1 column)
        _mol("nolen", "assigned", 0, _two_half_cols(True, True)),     # contested, read length unknown
        _mol("unc", "assigned", 0, _two_half_cols(True, True), contested=0),   # not contested: outside the denominator
        _mol("g1", "assigned", 0, _two_half_cols(True, True), fam="G"),  # a family whose only molecule has no length
    ]
    joins = [(f, c, "chrY") for f in ("F", "G") for c in (0, 1)]
    p, lp = write_run(str(tmp_path), stars, joins=joins, lengths={"both": 1000, "one": 1000, "none": 1000, "unc": 1000})
    out = str(tmp_path / "o.json")
    stdout, _ = run_main(p, lp, out, {"POPULATION": pop})
    res = json.load(open(out))
    pb = res["population"]
    assert (pb["n_contested"], pb["n_contested_with_length"]) == (5, 3)          # both, one, none (+ nolen, g1 w/o length)
    assert (pb["some_half_assigns"], pb["some_half_assigns_frac"]) == (2, round(2 / 3, 4))
    assert (pb["both_halves_assign"], pb["both_halves_assign_frac"]) == (1, round(1 / 3, 4))
    f = res["per_family"]["F"]
    assert (f["n_contested"], f["n_contested_with_length"]) == (4, 3)
    assert (f["some_half_assigns"], f["some_half_assigns_frac"]) == (2, round(2 / 3, 4))
    assert (f["both_halves_assign"], f["both_halves_assign_frac"]) == (1, round(1 / 3, 4))
    g = res["per_family"]["G"]
    assert (g["n_contested"], g["n_contested_with_length"]) == (1, 0)
    assert g["some_half_assigns"] == 0 and g["some_half_assigns_frac"] is None and g["both_halves_assign_frac"] is None
    # numerator and denominator come from the same length check: B_P1's n is the population block's numerator
    assert res["B_P1_split_half_concordance"]["n"] == pb["both_halves_assign"]
    # per_family is printed to stdout under POPULATION
    assert json.loads(stdout)["per_family"] == res["per_family"]


def test_denominator_fields_are_none_without_any_length(tmp_path):
    pop = str(tmp_path / "pop.tsv")
    _all_fixed_sites(pop)
    p, lp = write_run(str(tmp_path), [_mol("m", "assigned", 0, _two_half_cols(True, True))],
                      joins=[("F", 0, "chrY"), ("F", 1, "chrY")], lengths={})
    out = str(tmp_path / "o.json")
    run_main(p, lp, out, {"POPULATION": pop})
    pb = json.load(open(out))["population"]
    assert (pb["n_contested"], pb["n_contested_with_length"]) == (1, 0)
    assert pb["some_half_assigns_frac"] is None and pb["both_halves_assign_frac"] is None


def test_per_family_is_not_truncated_under_population_but_is_without(tmp_path):
    pop = str(tmp_path / "pop.tsv")
    _all_fixed_sites(pop)
    fams = [f"F{i:02d}" for i in range(20)]
    stars = [_mol(f"m{i}", "assigned", 0, _two_half_cols(True, True), fam=f) for i, f in enumerate(fams)]
    joins = [(f, c, "chrY") for f in fams for c in (0, 1)]
    p, lp = write_run(str(tmp_path), stars, joins=joins, lengths={f"m{i}": 1000 for i in range(20)})
    out = str(tmp_path / "o.json")
    run_main(p, lp, out, {"POPULATION": pop})
    assert sorted(json.load(open(out))["per_family"]) == fams
    run_main(p, lp, out)                                    # no POPULATION: the original top-15 cut
    assert len(json.load(open(out))["per_family"]) == 15


# ---- no POPULATION: output byte-identical to the pre-change code ------------------------------------------------------

def synth_run(d, positions):
    """40 deterministic molecules over two families; golden output below was written by the pre-change split_half.py."""
    rnd = random.Random(20261013)
    stars, lengths = [], {}
    for i in range(40):
        fam = "F1" if i % 3 else "F2"
        nc = rnd.choice([2, 3, 4])
        cand = sorted(rnd.sample(range(6), nc))
        truth = rnd.randrange(nc)
        ncols = rnd.choice([0, 3, 8, 14, 20])
        cols = []
        for j in range(ncols):
            al = "".join(rnd.choice("ACGT.") if rnd.random() < 0.1 else rnd.choice("ACGT") for _ in range(nc))
            obs = al[truth] if (al[truth] not in ".-" and rnd.random() < 0.95) else rnd.choice("ACGT")
            rp = rnd.randrange(0, 1000)
            cols.append((rp, obs, al, [1000 * (c + 1) + rp if rnd.random() < 0.9 else None for c in range(nc)]))
        cols.sort(key=lambda c: c[0])
        st = rnd.choice(["assigned", "assigned", "ambiguous", "tied"])
        stars.append(dict(read=f"read{i}", fam=fam, status=st, copy=cand[truth], cand=cand, cols=cols,
                          contested=rnd.choice([0, 1]), origin_rejected=int(i % 11 == 0), tie_outside=int(i % 13 == 0)))
        if i % 17 != 5:
            lengths[f"read{i}"] = 1000
    return write_run(d, stars, positions=positions, joins=[(f, c, "chrY") for f in ("F1", "F2") for c in range(6)],
                     lengths=lengths)


# written by the pre-change split_half.py (git d4823af8) on synth_run(); main() output file and stdout
GOLDEN_OUT = '{\n "counts": {\n  "some_half_assigns": 32,\n  "both_halves_assign": 22,\n  "both_halves_agree": 22,\n  "p0_n": 33,\n  "p0_same": 16,\n  "p0_diff|tied->assigned": 8,\n  "p0_diff|tied->ambiguous": 2,\n  "o2_and_half_assign": 25,\n  "o2_and_half_agree": 25,\n  "p0_diff|ambiguous->assigned": 5,\n  "p0_diff|assigned->tied": 2,\n  "no_length": 3\n },\n "B_P0_reproduction": 0.4848,\n "B_P0_pass": false,\n "B_P1_split_half_concordance": {\n  "n": 22,\n  "agree": 22,\n  "rate": 1.0,\n  "verdict": "UNDERPOWERED"\n },\n "B_P2_full_vs_half": {\n  "n": 25,\n  "agree": 25,\n  "rate": 1.0,\n  "verdict": "UNDERPOWERED"\n },\n "per_family": {\n  "F1": {\n   "n": 24,\n   "some_half_assigns": 19,\n   "both": 13,\n   "agree": 13\n  },\n  "F2": {\n   "n": 13,\n   "some_half_assigns": 13,\n   "both": 9,\n   "agree": 9\n  }\n }\n}'
GOLDEN_STDOUT = '{\n "counts": {\n  "some_half_assigns": 32,\n  "both_halves_assign": 22,\n  "both_halves_agree": 22,\n  "p0_n": 33,\n  "p0_same": 16,\n  "p0_diff|tied->assigned": 8,\n  "p0_diff|tied->ambiguous": 2,\n  "o2_and_half_assign": 25,\n  "o2_and_half_agree": 25,\n  "p0_diff|ambiguous->assigned": 5,\n  "p0_diff|assigned->tied": 2,\n  "no_length": 3\n },\n "B_P0_reproduction": 0.4848,\n "B_P0_pass": false,\n "B_P1_split_half_concordance": {\n  "n": 22,\n  "agree": 22,\n  "rate": 1.0,\n  "verdict": "UNDERPOWERED"\n },\n "B_P2_full_vs_half": {\n  "n": 25,\n  "agree": 25,\n  "rate": 1.0,\n  "verdict": "UNDERPOWERED"\n }\n}\n'


def test_without_population_output_is_byte_identical_to_the_pre_change_code(tmp_path):
    for positions in (False, True):  # a dump with a positions column must not change a no-POPULATION run
        d = tmp_path / f"pos{positions}"
        d.mkdir()
        p, lp = synth_run(str(d), positions)
        out = str(d / "o.json")
        stdout, _ = run_main(p, lp, out)
        assert open(out).read() == GOLDEN_OUT
        assert stdout == GOLDEN_STDOUT
    assert json.loads(GOLDEN_OUT)["counts"]["both_halves_assign"] >= 3  # the fixture is not trivial
