#!/usr/bin/env python
"""Self-check for bench/family_rna_refine.py (the DEFAULT-ON RNA-only refinement stage).

Asserts:
  (a) LEGACY opt-out (--legacy / RUSTLE_RNA_ORACLE=0) writes nothing; shipped path untouched.
  (b) DEFAULT (no flag) writes the refined catalog and is deterministic (== --rna-oracle == env=1).
  (c) The allele demotion removes DHRSX + LOC129530050 (DNA-confirmed).
  (d) The residual-removed count matches RNA_ONLY_EDGE_ORACLE.md's recall-preserving row:
      residual remaining {allele 0, oversize 3, multifam 4}; shipped total 12; 6/12 named FP removed.

Run: /home/juanfra/miniforge3/bin/python bench/test_family_rna_refine.py
"""
import hashlib
import json
import os
import subprocess
import sys

BENCH = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(BENCH, "family_rna_refine.py")
OUT_TSV = os.path.join(BENCH, "family_rna_refine.tsv")
OUT_JSON = os.path.join(BENCH, "family_rna_refine.json")
SHIPPED = os.path.join(BENCH, "denovo_families.py")
PY = sys.executable


def _run(args, env_extra=None):
    env = dict(os.environ)
    env.pop("RUSTLE_RNA_ORACLE", None)
    if env_extra:
        env.update(env_extra)
    return subprocess.run([PY, SCRIPT] + args, cwd=BENCH, env=env,
                          capture_output=True, text=True)


def _md5(path):
    with open(path, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


def test_legacy_opt_out_writes_nothing():
    for p in (OUT_TSV, OUT_JSON):
        if os.path.exists(p):
            os.remove(p)
    shipped_before = _md5(SHIPPED)
    r = _run(["--legacy"])                          # legacy opt-out (flag form)
    assert r.returncode == 0, f"legacy exit={r.returncode}"
    assert r.stdout.startswith("legacy:"), repr(r.stdout)
    assert not os.path.exists(OUT_TSV) and not os.path.exists(OUT_JSON), "legacy wrote outputs"
    assert _md5(SHIPPED) == shipped_before, "denovo_families.py (shipped path) was modified"
    r2 = _run([], env_extra={"RUSTLE_RNA_ORACLE": "0"})   # legacy opt-out (env form)
    assert r2.returncode == 0 and r2.stdout.startswith("legacy:"), repr(r2.stdout)
    assert not os.path.exists(OUT_TSV), "env-legacy wrote outputs"
    print("(a) legacy opt-out (--legacy / env=0) writes nothing / shipped untouched : OK")


def test_default_deterministic():
    r1 = _run([])                                  # DEFAULT: no flag -> refined catalog
    assert r1.returncode == 0, r1.stderr
    assert os.path.exists(OUT_TSV) and os.path.exists(OUT_JSON), "default did not write the catalog"
    tsv1, json1 = _md5(OUT_TSV), _md5(OUT_JSON)
    r2 = _run(["--rna-oracle"])                     # deprecated no-op == default
    assert r2.returncode == 0, r2.stderr
    assert _md5(OUT_TSV) == tsv1 and _md5(OUT_JSON) == json1, "default != --rna-oracle"
    r3 = _run([], env_extra={"RUSTLE_RNA_ORACLE": "1"})   # explicit-enable == default
    assert r3.returncode == 0, r3.stderr
    assert _md5(OUT_TSV) == tsv1, "TSV not byte-identical across runs"
    assert _md5(OUT_JSON) == json1, "JSON not byte-identical across runs"
    print("(b) default deterministic (default == --rna-oracle == env=1; byte-identical) : OK")


def test_allele_demote_removes_known_fp():
    d = json.load(open(OUT_JSON))
    demoted = {x["gene"] for x in d["alleles_demoted"]}
    conf = {x["gene"] for x in d["alleles_demoted"] if x["dna_confirmed"]}
    assert "DHRSX" in demoted and "LOC129530050" in demoted, f"missing allele demotions: {demoted}"
    assert {"DHRSX", "LOC129530050"} <= conf, f"DHRSX/LOC129530050 not flagged DNA-confirmed: {conf}"
    print(f"(c) allele demotion removes DHRSX + LOC129530050 : OK  (demoted={sorted(demoted)})")


def test_residual_removed_matches_recall_preserving_row():
    d = json.load(open(OUT_JSON))
    rem = d["residual_fp"]["residual_remaining"]
    assert rem == {"allele": 0, "oversize": 3, "multifam": 4}, f"residual remaining drift: {rem}"
    assert d["residual_fp"]["shipped_residual_total"] == 12, d["residual_fp"]["shipped_residual_total"]
    assert d["residual_fp"]["named_removed"] == 6, d["residual_fp"]["named_removed"]
    bd = d["residual_fp"]["named_removed_breakdown"]
    assert bd == {"allele": 2, "oversize": 2, "multifam": 2}, f"removed breakdown drift: {bd}"
    print("(d) residual remaining {allele 0, oversize 3, multifam 4}; 6/12 removed : OK")


def test_rna_only_guard():
    d = json.load(open(OUT_JSON))
    g = d["guards"]
    assert g["edge_decision_features"] == ["core_recip", "aln_frac"], g["edge_decision_features"]
    assert g["demote_features"] == ["balanced_frac", "copy_like"], g["demote_features"]
    assert g["no_dna_in_inference"] is True and g["gamma"] == 0.2 and g["seed"] == 0
    print("(+) RNA-only inference guard (no DNA in decision) : OK")


if __name__ == "__main__":
    test_legacy_opt_out_writes_nothing()
    test_default_deterministic()
    test_allele_demote_removes_known_fp()
    test_residual_removed_matches_recall_preserving_row()
    test_rna_only_guard()
    print("\nALL TESTS PASSED")
