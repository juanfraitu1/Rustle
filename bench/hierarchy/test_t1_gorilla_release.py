#!/usr/bin/env python3
"""Tests of the release path of the Phase R2 report (docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md section 5 and Amendment 1), written after the
independent code review (two reviewers, 2026-10-07) and BEFORE the changes they asked for:

  * `report` itself refuses without the release token (also when compute_report is called directly), runs the positive control and the Gate 0 fingerprints
    of its own inputs, and prints INVALID and withholds every rate when either fails (exit 3);
  * the printout carries the Wilson and bootstrap intervals, the path and short-aligned pair counts, the B_viol caveat and its distance-class split, and the
    registered pair/family counts of Gate 6 next to the observed ones;
  * Gate 0 records the interpreter.

The second half is the reviewer's supplementary suite (it kills 29 mutants that the builder's 74 tests let survive), adapted to the new `report` arguments.
Run: python3 -B -m unittest bench/hierarchy/test_t1_gorilla_release.py   (made-up worlds and scratch directories only)
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import t1_gorilla_pairs as T  # noqa: E402
import test_t1_gorilla_pairs as B  # noqa: E402  (builder helpers: sedef_row, make_world, add_control, release_files, report_argv, run_main ...)

sedef_row = B.sedef_row


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


class ReleaseUnitTests(unittest.TestCase):
    def test_require_release_reads_the_first_line_only(self):
        with tempfile.TemporaryDirectory() as d:
            ok = os.path.join(d, "ok")
            with open(ok, "w") as fh:
                fh.write(T.RELEASE_TOKEN + "\nanything\n")
            T.require_release(ok)                                         # no exception
            late = os.path.join(d, "late")
            with open(late, "w") as fh:
                fh.write("\n" + T.RELEASE_TOKEN + "\n")
            for bad in (late, os.path.join(d, "missing"), None):
                with self.assertRaises(T.ReleaseRefused):
                    T.require_release(bad)

    def test_compute_report_refuses_by_itself(self):
        a = argparse.Namespace(release_file=None, sedef="/nonexistent/s", genes="/nonexistent/g", truth="/nonexistent/t", pairs="/nonexistent/p",
                               liftoff="/nonexistent/l", registry=None, outdir="/nonexistent/o", rebuild=False)
        with self.assertRaises(T.ReleaseRefused):
            T.compute_report(a)

    def test_b_viol_by_class_splits_the_pooled_count(self):
        U = {1: frozenset({1}), 2: frozenset({1}), 3: frozenset({2}), 4: frozenset({2})}
        control = [("A", "B", 1, 2, "same<100kb"), ("A", "B", 1, 3, "same<100kb"), ("A", "C", 2, 4, "cross-contig"), ("B", "C", 3, 4, "cross-contig"),
                   ("A", "C", 1, 4, "cross-contig")]
        by = T.b_viol_by_class(control, U)
        self.assertEqual(by, {"cross-contig": (1, 3), "same<100kb": (1, 2)})
        self.assertEqual(sum(k for k, _n in by.values()), T.b_viol(control, U)[0])
        self.assertEqual(sum(n for _k, n in by.values()), T.b_viol(control, U)[1])

    def rows(self, **override):
        out = []
        for n in T.REQUIRED_INPUTS:
            r = {"name": n, "ok": True, "sha256": "a" * 64, "expected_prefix": "a" * 16}
            r.update(override.get(n, {}))
            out.append(r)
        return out

    def test_validity_names_every_reason(self):
        passed = {"passed": True, "without_class": [], "connected": True}
        failed = {"passed": False, "without_class": ["L7"], "connected": True}
        self.assertEqual(T.validity(self.rows(), passed), (True, []))
        ok, why = T.validity(self.rows(genes={"ok": False, "sha256": "b" * 64, "expected_prefix": "c" * 16}), passed)
        self.assertFalse(ok)
        self.assertTrue(any("genes" in r for r in why))
        ok, why = T.validity(self.rows(), failed)
        self.assertFalse(ok)
        self.assertTrue(any("L7" in r for r in why))
        ok, why = T.validity(self.rows(), None)
        self.assertFalse(ok)
        self.assertTrue(any("not run" in r for r in why))
        ok, why = T.validity(self.rows(genes={"ok": False, "sha256": "b" * 64, "expected_prefix": "c" * 16}), failed)
        self.assertEqual(len(why), 2)

    def test_validity_demands_all_six_registered_inputs(self):
        passed = {"passed": True, "without_class": [], "connected": True}
        self.assertEqual(T.REQUIRED_INPUTS, ("sedef", "genes", "truth", "pairs", "liftoff", "dna_sd_atoms.py"))
        ok, why = T.validity([], passed)
        self.assertFalse(ok)
        self.assertEqual(len(why), 6)
        for drop in T.REQUIRED_INPUTS:
            ok, why = T.validity([r for r in self.rows() if r["name"] != drop], passed)
            self.assertFalse(ok, drop)
            self.assertTrue(any(drop in r for r in why), drop)
        ok, why = T.validity(self.rows(truth={"expected_prefix": None}), passed)
        self.assertFalse(ok)
        self.assertTrue(any("truth" in r and "registered" in r for r in why))

    def test_validity_reports_a_control_that_could_not_be_computed(self):
        ok, why = T.validity(self.rows(), None, control_error="ValueError: 3 gene-LRPAP1 rows")
        self.assertFalse(ok)
        self.assertEqual(len(why), 1)
        self.assertIn("could not be computed", why[0])
        self.assertIn("3 gene-LRPAP1 rows", why[0])

    def test_the_b_viol_caveat_does_not_call_the_pool_small(self):
        self.assertIn("ceiling-type", T.B_VIOL_CAVEAT)
        self.assertIn("not identity-matched", T.B_VIOL_CAVEAT)
        self.assertNotIn("small", T.B_VIOL_CAVEAT)

    def test_environment_record_names_the_interpreter_and_libraries(self):
        e = T.environment_record()
        self.assertEqual(e["python_version"], sys.version.split()[0])
        self.assertEqual(e["executable"], sys.executable)
        self.assertEqual(sorted(e["modules"]), ["numpy", "pysam", "scipy"])

    def test_read_registry(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "r.tsv")
            with open(p, "w") as fh:
                fh.write("sedef\t41aa1c9c53a18d68\n\ngenes\t6f9bfbd776e8ecc8\n")
            self.assertEqual(T.read_registry(p), {"sedef": "41aa1c9c53a18d68", "genes": "6f9bfbd776e8ecc8"})

    def test_the_default_registry_is_the_registered_constants(self):
        r = T.default_registry()
        self.assertEqual(r["sedef"], "41aa1c9c53a18d68")
        self.assertEqual(r["dna_sd_atoms.py"], T.ATOMS_SHA256[:16])
        self.assertEqual(sorted(r), ["dna_sd_atoms.py", "genes", "liftoff", "pairs", "sedef", "truth"])

    def table(self, n_pairs=38, n_families=33):
        pw = 30 / n_pairs if n_pairs else None
        return {"tau": 0.90, "n_pairs": n_pairs, "n_families": n_families, "shared": 30, "pair_weighted": pw, "family_weighted": 0.75, "a_viol": (1 - pw) if pw is not None else None,
                "wilson": [0.62, 0.88], "bootstrap": {"B": 2000, "seed": 20260930, "n_families": n_families, "pair_weighted": [0.6, 0.9], "family_weighted": [0.55, 0.92]},
                "b_viol": [5, 100], "b_viol_by_class": {"cross-contig": [4, 90], "same<100kb": [1, 10]}, "underpowered": False,
                "short_aligned_pairs": 4, "path_pairs": 23}

    def test_format_report_prints_every_number_the_reviewers_asked_for(self):
        rep = {"valid": True, "tau": {"0.90": self.table(), "0.95": self.table(15, 15), "0.98": self.table(5, 5)}}
        text = "\n".join(T.format_report(rep))
        for needle in ("Wilson", "0.62", "0.88", "bootstrap", "2000", "20260930", "0.600", "0.900", "path pairs", "23 of 38", "aligned < 0.5", "4",
                       "B_viol", "5 of 100", "same<100kb 1/10", "cross-contig 4/90", "ceiling-type", "registered 38/33: match", "registered 15/15: match",
                       "registered 5/5: match"):
            self.assertIn(needle, text, needle)
        text = "\n".join(T.format_report({"valid": True, "tau": {"0.90": self.table(37, 33)}}))
        self.assertIn("registered 38/33: DIFFERS", text)

    def test_format_report_handles_an_empty_table(self):
        t = self.table(0, 0)
        t.update(pair_weighted=None, family_weighted=None, a_viol=None, wilson=[None, None], bootstrap=None, b_viol=[0, 0], b_viol_by_class={}, underpowered=True,
                 short_aligned_pairs=0, path_pairs=0, shared=0)
        text = "\n".join(T.format_report({"valid": True, "tau": {"0.98": t}}))
        self.assertIn("NA", text)
        self.assertIn("UNDERPOWERED True", text)


class ReleaseEndToEndTests(unittest.TestCase):
    def test_a_valid_world_releases_the_tables_and_records_the_gates(self):
        d, paths = B.world_with_control()
        reg, rel = B.release_files(d, paths)
        out = os.path.join(d, "out")
        rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
        self.assertEqual(rc, 0, text)
        rep = load_json(os.path.join(out, "report.json"))
        self.assertTrue(rep["valid"])
        g = rep["gates"]
        self.assertTrue(g["positive_control"]["passed"])
        self.assertEqual(sorted(r["name"] for r in g["inputs"]), ["dna_sd_atoms.py", "genes", "liftoff", "pairs", "sedef", "truth"])
        self.assertTrue(all(r["ok"] for r in g["inputs"]))
        self.assertEqual(g["environment"]["python_version"], sys.version.split()[0])
        self.assertEqual(sum(k for k, _n in rep["tau"]["0.90"]["b_viol_by_class"].values()), rep["tau"]["0.90"]["b_viol"][0])
        self.assertEqual(sum(n for _k, n in rep["tau"]["0.90"]["b_viol_by_class"].values()), rep["tau"]["0.90"]["b_viol"][1])
        for k in ("0.90", "0.95", "0.98"):
            self.assertTrue(os.path.exists(os.path.join(out, f"report_pairs_{k}.tsv")))
        self.assertIn("Wilson", text)
        self.assertIn("ceiling-type", text)
        self.assertIn("registered 38/33: DIFFERS", text)       # the made-up world is not the real one

    def test_a_failing_positive_control_makes_the_report_invalid_and_withholds_the_rates(self):
        d, paths = B.world_with_control(scatter_one=True)
        reg, rel = B.release_files(d, paths)
        out = os.path.join(d, "out")
        rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
        self.assertEqual(rc, 3, text)
        self.assertIn("INVALID", text)
        self.assertIn("L7", text)
        self.assertNotIn("S-rate", text)
        rep = load_json(os.path.join(out, "report.json"))
        self.assertFalse(rep["valid"])
        self.assertNotIn("tau", rep)
        self.assertFalse(rep["gates"]["positive_control"]["passed"])
        self.assertEqual([f for f in os.listdir(out) if f.startswith("report_pairs_")], [])

    def test_a_registry_mismatch_makes_the_report_invalid_and_withholds_the_rates(self):
        d, paths = B.world_with_control()
        reg, rel = B.release_files(d, paths)
        lines = open(reg).read().splitlines()
        with open(reg, "w") as fh:
            fh.write("\n".join(("pairs\t" + "0" * 16) if ln.startswith("pairs\t") else ln for ln in lines) + "\n")
        out = os.path.join(d, "out")
        rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
        self.assertEqual(rc, 3, text)
        self.assertIn("INVALID", text)
        self.assertIn("pairs", text)
        rep = load_json(os.path.join(out, "report.json"))
        self.assertFalse(rep["valid"])
        self.assertNotIn("tau", rep)
        self.assertEqual([f for f in os.listdir(out) if f.startswith("report_pairs_")], [])

    def test_a_missing_self_lift_table_is_refused_before_anything_is_computed(self):
        d, paths = B.world_with_control()
        reg, rel = B.release_files(d, paths)
        paths = dict(paths, liftoff=os.path.join(d, "no_such_liftoff.tsv"))
        out = os.path.join(d, "out")
        rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
        self.assertEqual(rc, 2, text)
        self.assertIn("liftoff", text.lower())
        self.assertFalse(os.path.exists(out))

    def test_a_registry_that_omits_an_input_is_refused_not_silently_passed(self):
        d, paths = B.world_with_control()
        reg, rel = B.release_files(d, paths)
        lines = [ln for ln in open(reg).read().splitlines() if not ln.startswith("pairs\t")]
        with open(reg, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        out = os.path.join(d, "out")
        rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
        self.assertEqual(rc, 2, text)
        self.assertIn("registry", text.lower())
        self.assertFalse(os.path.exists(os.path.join(out, "report.json")))

    def test_counts_print_the_registered_recon_values_and_explain_the_difference(self):
        d, paths = B.world_with_control()
        out = os.path.join(d, "out")
        rc, text = B.run_main(["counts", "--sedef", paths["sedef"], "--genes", paths["genes"], "--truth", paths["truth"], "--pairs", paths["pairs"],
                               "--liftoff", paths["liftoff"], "--outdir", out])
        self.assertEqual(rc, 0, text)
        self.assertIn("registered recon", text)
        self.assertIn("40631", text.replace(",", ""))
        self.assertIn("12483", text.replace(",", ""))
        self.assertIn("372793", text.replace(",", ""))
        self.assertIn("un-deduplicated", text)
        self.assertIn("coverage > 0.5", text)

    def test_a_registry_with_a_short_or_empty_prefix_is_refused(self):
        for bad in ("", "abc", "0" * 15):
            d, paths = B.world_with_control()
            reg, rel = B.release_files(d, paths)
            lines = open(reg).read().splitlines()
            with open(reg, "w") as fh:
                fh.write("\n".join(("pairs\t" + bad) if ln.startswith("pairs\t") else ln for ln in lines) + "\n")
            out = os.path.join(d, "out")
            rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
            self.assertEqual(rc, 2, (bad, text))
            self.assertIn("prefix", text.lower())
            self.assertFalse(os.path.exists(os.path.join(out, "report.json")))

    def test_a_wrong_self_lift_table_is_invalid_with_a_reason_not_a_traceback(self):
        d, paths = B.world_with_control()
        rows = open(paths["liftoff"]).read().splitlines()
        with open(paths["liftoff"], "w") as fh:
            fh.write("\n".join(rows[:4]) + "\n")                    # header + 3 copies: lrpap_copies expects 8
        reg, rel = B.release_files(d, paths)
        out = os.path.join(d, "out")
        rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
        self.assertEqual(rc, 3, text)
        self.assertIn("INVALID", text)
        self.assertIn("could not be computed", text)
        rep = load_json(os.path.join(out, "report.json"))
        self.assertFalse(rep["valid"])
        self.assertNotIn("tau", rep)
        self.assertIsNone(rep["gates"]["positive_control"])
        self.assertIn("gene-LRPAP1", rep["gates"]["positive_control_error"])

    def test_a_liftoff_file_with_other_columns_is_invalid_too(self):
        d, paths = B.world_with_control()
        with open(paths["liftoff"], "w") as fh:
            fh.write("a\tb\tc\n1\t2\t3\n")
        reg, rel = B.release_files(d, paths)
        out = os.path.join(d, "out")
        rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
        self.assertEqual(rc, 3, text)
        self.assertIn("could not be computed", text)

    def test_a_bad_registry_is_refused_before_any_atoms_are_built(self):
        d, paths = B.world_with_control()
        reg, rel = B.release_files(d, paths)
        lines = [ln for ln in open(reg).read().splitlines() if not ln.startswith("genes\t")]
        with open(reg, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        out = os.path.join(d, "out")
        rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
        self.assertEqual(rc, 2, text)
        self.assertFalse(os.path.exists(out))                       # not even the atoms directory

    def test_the_report_records_its_registry_and_its_own_hash(self):
        import hashlib
        d, paths = B.world_with_control()
        reg, rel = B.release_files(d, paths)
        out = os.path.join(d, "out")
        rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
        self.assertEqual(rc, 0, text)
        g = load_json(os.path.join(out, "report.json"))["gates"]
        self.assertEqual(g["registry"], reg)
        self.assertEqual(g["scorer_sha256"], hashlib.sha256(open(T.__file__, "rb").read()).hexdigest())

    def test_an_invalid_rerun_removes_the_stale_pair_tables_of_an_earlier_valid_run(self):
        d, paths = B.world_with_control()
        reg, rel = B.release_files(d, paths)
        out = os.path.join(d, "out")
        rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
        self.assertEqual(rc, 0, text)
        self.assertTrue(os.path.exists(os.path.join(out, "report_pairs_0.90.tsv")))
        lines = open(reg).read().splitlines()
        with open(reg, "w") as fh:
            fh.write("\n".join(("pairs\t" + "0" * 16) if ln.startswith("pairs\t") else ln for ln in lines) + "\n")
        rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
        self.assertEqual(rc, 3, text)
        self.assertEqual([f for f in os.listdir(out) if f.startswith("report_pairs_")], [])
        self.assertFalse(load_json(os.path.join(out, "report.json"))["valid"])

    def test_the_released_report_is_identical_under_two_hash_seeds(self):
        d, paths = B.world_with_control()
        reg, rel = B.release_files(d, paths)
        outs = []
        for seed in ("0", "1"):
            out = os.path.join(d, f"out{seed}")
            p = subprocess.run([sys.executable, "-B", os.path.join(HERE, "t1_gorilla_pairs.py")] + B.report_argv(paths, reg, rel, out),
                               capture_output=True, text=True, env=dict(os.environ, PYTHONHASHSEED=seed))
            self.assertEqual(p.returncode, 0, p.stderr)
            outs.append(out)
        names = sorted(f for f in os.listdir(outs[0]) if f.startswith("report"))
        self.assertEqual(names, sorted(f for f in os.listdir(outs[1]) if f.startswith("report")))
        for n in names:
            a, b = (open(os.path.join(o, n)).read() for o in outs)
            if n == "report.json":
                a, b = (json.loads(x) for x in (a, b))
                for x in (a, b):
                    x["gates"].pop("environment")
                    for r in x["gates"]["inputs"]:
                        r.pop("path")
            self.assertEqual(a, b, n)

    def test_gate0_records_the_interpreter_and_passes_against_a_matching_registry(self):
        d, paths = B.world_with_control()
        reg, _rel = B.release_files(d, paths)
        out = os.path.join(d, "out")
        rc, text = B.run_main(["gate0", "--sedef", paths["sedef"], "--genes", paths["genes"], "--truth", paths["truth"], "--pairs", paths["pairs"],
                               "--liftoff", paths["liftoff"], "--registry", reg, "--outdir", out])
        self.assertEqual(rc, 0, text)
        env = open(os.path.join(out, "gate0_r2_env.tsv")).read()
        self.assertIn("python_version\t" + sys.version.split()[0], env)
        self.assertIn("executable\t" + sys.executable, env)
        self.assertIn("numpy\t", env)
        self.assertIn("python " + sys.version.split()[0], text)


# ---------------------------------------------------------------------------------------------------------------------------
# the reviewer's supplementary tests (each kills a mutant that the builder's 74 tests let survive), adapted to the new report arguments
# ---------------------------------------------------------------------------------------------------------------------------
def row_with(c1, s1, e1, c2, s2, e2, frac=0.95, aln_len=None, cigar=None, col22=None, strand2="+"):
    line = sedef_row(c1, s1, e1, c2, s2, e2, strand2=strand2, cigar=cigar, frac=frac, aln_len=aln_len)
    f = line.split("\t")
    if col22 is not None:
        f[21] = f"{col22:.6f}"
    return "\t".join(f)


class SedefColumnTests(unittest.TestCase):
    def test_s02_the_tau_column_is_21_not_22(self):
        rows, _ = T.read_sedef([row_with("A", 0, 5000, "B", 0, 5000, frac=0.95, col22=0.50)])
        self.assertAlmostEqual(rows[0].frac, 0.95, places=6)
        self.assertEqual(len(T.retain(rows, 0.90)), 1)

    def test_s03_dedupe_priority_is_fracmatch_before_aln_len(self):
        def tagged(frac, aln_len, tag):
            f = row_with("A", 0, 5000, "B", 10, 5010, frac=frac, aln_len=aln_len).split("\t")
            f[12] = f[12] + tag
            return "\t".join(f)
        rows, _ = T.read_sedef([tagged(0.97, 100, ";keepme"), tagged(0.95, 300, ";dropme")])
        kept, removed = T.dedupe(rows)
        self.assertEqual(removed, 1)
        self.assertIn(";keepme", kept[0].line)

    def test_s05_s09_s10_tuples_that_differ_in_one_coordinate_are_different_tuples(self):
        base = ("A", 0, 5000, "B", 10, 5010)
        variants = [("A", 0, 5000, "B", 10, 5011), ("A", 1, 5000, "B", 10, 5010), ("A", 0, 5000, "C", 10, 5010), ("A", 0, 5001, "B", 10, 5010), ("A", 0, 5000, "B", 11, 5010)]
        for v in variants:
            rows, _ = T.read_sedef([sedef_row(*base), sedef_row(*v)])
            kept, removed = T.dedupe(rows)
            self.assertEqual((len(kept), removed), (2, 0), v)


class ClassEdgeReadingTests(unittest.TestCase):
    """C01/C02/C03: an A-B alignment that covers 30% of the 10 kb atom A and 75% of the 4 kb atom B1."""

    def make(self):
        d = tempfile.mkdtemp(prefix="supp_cls_")
        lines = [
            # A [0,10000) aligned to B [100000,110000) over the first 3000 bp only: 3000M7000D7000I  (A: M+D = 10000, B: M+I = 10000)
            row_with("A", 0, 10000, "B", 100000, 110000, frac=0.95, cigar="3000M7000D7000I"),
            # a second row cuts B at 104000: B [100000,104000) vs C [0,4000), fully aligned
            row_with("B", 100000, 104000, "C", 0, 4000, frac=0.95, cigar="4000M"),
        ]
        rows, _ = T.read_sedef(lines)
        rows, _ = T.dedupe(rows)
        return T.prepare_tau(rows, 0.90, d)

    def test_the_two_sided_minimum_coverage_decides_the_class_edge(self):
        info = self.make()
        edges = T.read_edges(info["prefix"])
        cov = {(i, j): c for i, j, _id, c in edges}
        self.assertTrue(any(abs(c - 0.3) < 1e-6 for c in cov.values()), cov)
        # classes: {B1, C} (full coverage); A and B2 stay alone -> 4 atoms, 3 classes
        self.assertEqual((info["atoms"], info["classes"]), (4, 3))


class DevelopmentSliceTests(unittest.TestCase):
    def test_g02_to_g05_all_four_development_contigs_are_registered(self):
        self.assertEqual(set(T.DEV_CONTIGS), {"NC_073241.2", "NC_073242.2", "NC_073244.2", "NC_073234.2"})

    def test_a_family_with_one_member_on_any_development_contig_leaves_the_verdict_set(self):
        for c in T.DEV_CONTIGS:
            genes_txt = "\n".join(["contig\tname\ttype\tbiotype\tstart1\tend\texons", f"{c}\tx\tgene\tpc\t1001\t2000\t1000-2000", "C1\ty\tgene\tpc\t1001\t2000\t1000-2000",
                                   "C1\tz\tgene\tpc\t5001\t6000\t5000-6000"])
            genes, index, _ = T.read_genes(genes_txt.splitlines())
            truth = "\n".join([B.TRUTH_HEAD, f"F1\twhole\t1to1\tgene-x\t{c}\t1001\t2000", "F1\twhole\t1to1\tgene-y\tC1\t1001\t2000", "F2\twhole\t1to1\tgene-z\tC1\t5001\t6000"])
            layer = T.read_layer(truth.splitlines(), index)
            self.assertEqual(sorted(T.verdict_families(layer, genes)), ["F2"], c)


class EligibilityBoundaryTests(unittest.TestCase):
    def test_p01_the_first_gene_must_have_a_class_too(self):
        genes, _i, _d = T.read_genes(["contig\tname\ttype\tbiotype\tstart1\tend\texons", "C1\tx\tgene\tpc\t1001\t2000\t1000-2000", "C1\ty\tgene\tpc\t5001\t6000\t5000-6000"])
        self.assertFalse(T.pair_eligible(0, 1, genes, {0: frozenset(), 1: frozenset({3})}))
        self.assertFalse(T.pair_eligible(1, 0, genes, {0: frozenset(), 1: frozenset({3})}))
        self.assertTrue(T.pair_eligible(0, 1, genes, {0: frozenset({3}), 1: frozenset({3})}))

    def test_p03_a_single_shared_exonic_bp_makes_a_pair_ineligible(self):
        genes, _i, _d = T.read_genes(["contig\tname\ttype\tbiotype\tstart1\tend\texons", "C1\tx\tgene\tpc\t1001\t2000\t1000-2000", "C1\ty\tgene\tpc\t2000\t3000\t1999-3000"])
        self.assertEqual(T.shared_exonic_bp(genes[0], genes[1]), 1)
        self.assertFalse(T.pair_eligible(0, 1, genes, {0: frozenset({1}), 1: frozenset({1})}))


def world_with_dev_family(d):
    """make_world plus the control, plus: family F6 with a member on NC_073241.2 (development slice) and one on NC_073224.2, joined by a 0.99 SD row; family F7
    with two genes at identity 0.92 whose genes sit in a 0.99 SD pair (so they have atoms at tau 0.98 although their pairs.tsv identity is 0.92)."""
    paths = B.make_world(d)
    sedef = [x for x in open(paths["sedef"]).read().splitlines() if not x.startswith("#")]
    sedef += [
        row_with("NC_073224.2", 120000, 126000, "NC_073241.2", 10000, 16000, frac=0.99),
        row_with("NC_073224.2", 150000, 156000, "NC_073227.2", 700000, 706000, frac=0.99),
    ]
    sedef.append(B.SEDEF_HEADER)
    genes = open(paths["genes"]).read().splitlines()
    genes += ["NC_073224.2\tgE1\tgene\tprotein_coding\t120501\t125000\t121000-124000",
              "NC_073241.2\tgE2\tgene\tprotein_coding\t10501\t15000\t11000-14000",
              "NC_073224.2\tgF1\tgene\tprotein_coding\t150501\t155000\t151000-154000",
              "NC_073227.2\tgF2\tgene\tprotein_coding\t700501\t705000\t701000-704000"]
    truth = open(paths["truth"]).read().splitlines()
    truth += ["F6\twhole\t1to1\tgene-gE1\tNC_073224.2\t120501\t125000", "F6\twhole\t1to1\tgene-gE2\tNC_073241.2\t10501\t15000",
              "F7\twhole\t1to1\tgene-gF1\tNC_073224.2\t150501\t155000", "F7\twhole\t1to1\tgene-gF2\tNC_073227.2\t700501\t705000"]
    pairs = open(paths["pairs"]).read().splitlines()
    pairs += [B.pairs_row("F6", "gene-gE1", "gene-gE2", "0.97"), B.pairs_row("F7", "gene-gF1", "gene-gF2", "0.92")]
    for name, lines in (("sedef", sedef), ("genes", genes), ("truth", truth), ("pairs", pairs)):
        with open(paths[name], "w") as fh:
            fh.write("\n".join(lines) + "\n")
    return B.add_control(d, paths)


class WiringTests(unittest.TestCase):
    def test_w04_to_w06_counts_keep_verdict_and_all_apart(self):
        d = tempfile.mkdtemp(prefix="supp_wire_")
        paths = world_with_dev_family(d)
        out = os.path.join(d, "out")
        rc = T.main(["counts", "--sedef", paths["sedef"], "--genes", paths["genes"], "--truth", paths["truth"], "--pairs", paths["pairs"], "--outdir", out])
        self.assertEqual(rc, 0)
        t = load_json(os.path.join(out, "counts.json"))["tau"]["0.90"]
        # F6 has a member on the development slice: it is in 'all' but not in the verdict set
        self.assertEqual(sum(t["same_family_pairs"]["all"].values()) - sum(t["same_family_pairs"]["verdict"].values()), 1)
        self.assertEqual(t["depth_matched"]["all"]["pairs"] - t["depth_matched"]["verdict"]["pairs"], 1)
        self.assertGreater(sum(t["control_pairs"]["all"].values()), sum(t["control_pairs"]["verdict"].values()))

    def test_w01_w02_w08_the_report_uses_the_verdict_set_and_its_own_tau(self):
        d = tempfile.mkdtemp(prefix="supp_rep_")
        paths = world_with_dev_family(d)
        reg, rel = B.release_files(d, paths)
        out = os.path.join(d, "out")
        rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
        self.assertEqual(rc, 0, text)
        rep = load_json(os.path.join(out, "report.json"))
        fams90 = {ln.split("\t")[0] for ln in open(os.path.join(out, "report_pairs_0.90.tsv")).read().splitlines()[1:]}
        self.assertNotIn("F6", fams90)                       # development-slice family out of the verdict set (W01)
        self.assertEqual(rep["tau"]["0.98"]["n_pairs"], 0)   # F7 has atoms at 0.98 but identity 0.92 (W08)
        self.assertEqual(rep["tau"]["0.90"]["b_viol"][1], len(T.control_pairs(*self._universe(paths, 0.90))))   # W02

    def _universe(self, paths, tau):
        genes, index, _ = T.read_genes(open(paths["genes"]))
        layer = T.read_layer(open(paths["truth"]), index)
        verdict = T.verdict_families(layer, genes)
        d = tempfile.mkdtemp(prefix="supp_u_")
        rows, _ = T.read_sedef(open(paths["sedef"]))
        rows, _ = T.dedupe(rows)
        info = T.prepare_tau(rows, tau, d)
        U = T.universe_for_tau(info, genes, layer)
        return layer, genes, U, verdict

    def dup_world(self):
        d = tempfile.mkdtemp(prefix="supp_dd_")
        paths = B.make_world(d)
        lines = [x for x in open(paths["sedef"]).read().splitlines() if not x.startswith("#")]
        # the same tuple twice: the 0.97 row aligns 10% of the atom pair (kept), the 0.91 row aligns all of it (removed by the identical-tuple rule)
        lines += [row_with("NC_073224.2", 200000, 205000, "NC_073227.2", 800000, 805000, frac=0.97, cigar="500M4500D4500I"),
                  row_with("NC_073224.2", 200000, 205000, "NC_073227.2", 800000, 805000, frac=0.91, cigar="5000M"), B.SEDEF_HEADER]
        genes = open(paths["genes"]).read().splitlines() + ["NC_073224.2\tgG1\tgene\tprotein_coding\t200501\t204500\t201000-204000",
                                                          "NC_073227.2\tgG2\tgene\tprotein_coding\t800501\t804500\t801000-804000"]
        truth = open(paths["truth"]).read().splitlines() + ["F8\twhole\t1to1\tgene-gG1\tNC_073224.2\t200501\t204500", "F8\twhole\t1to1\tgene-gG2\tNC_073227.2\t800501\t804500"]
        pairs = open(paths["pairs"]).read().splitlines() + [B.pairs_row("F8", "gene-gG1", "gene-gG2", "0.95")]
        for name, ls in (("sedef", lines), ("genes", genes), ("truth", truth), ("pairs", pairs)):
            with open(paths[name], "w") as fh:
                fh.write("\n".join(ls) + "\n")
        return d, B.add_control(d, paths)

    def test_w07_counts_build_atoms_from_the_deduplicated_rows(self):
        d, paths = self.dup_world()
        out = os.path.join(d, "out")
        self.assertEqual(T.main(["counts", "--sedef", paths["sedef"], "--genes", paths["genes"], "--truth", paths["truth"], "--pairs", paths["pairs"], "--outdir", out]), 0)
        t = load_json(os.path.join(out, "counts.json"))["tau"]["0.90"]
        self.assertEqual(t["rows_retained"], 3)   # make_world's two 0.9+ rows and ONE of the two duplicates

    def test_w03_the_report_builds_atoms_from_the_deduplicated_rows(self):
        d, paths = self.dup_world()
        reg, rel = B.release_files(d, paths)
        out = os.path.join(d, "out")
        rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
        self.assertEqual(rc, 0, text)
        rows = [ln.split("\t") for ln in open(os.path.join(out, "report_pairs_0.90.tsv")).read().splitlines()]
        head, body = rows[0], rows[1:]
        f8 = [r for r in body if r[0] == "F8"]
        self.assertEqual(len(f8), 1)
        self.assertEqual(f8[0][head.index("S")], "0")      # the kept duplicate covers 10% of the pair: no class edge


class BootstrapValueTests(unittest.TestCase):
    def test_r05_r06_r07_bootstrap_values_equal_a_literal_reimplementation(self):
        import random
        by_family = {"B": [0], "A": [1, 1], "C": [1, 0, 0]}
        pw, fw = T.bootstrap_values(by_family, B=50, seed=7)
        fams = ["A", "B", "C"]
        rng = random.Random(7)
        epw, efw = [], []
        for _ in range(50):
            draw = [fams[rng.randrange(3)] for _ in range(3)]
            n = sum(len(by_family[f]) for f in draw)
            epw.append(sum(sum(by_family[f]) for f in draw) / n)
            efw.append(sum(sum(by_family[f]) / len(by_family[f]) for f in draw) / 3)
        self.assertEqual((pw, fw), (epw, efw))


class MoreWiringTests(unittest.TestCase):
    def test_m02_depth_matching_uses_identity_not_identity_cs(self):
        row = B.pairs_row("F1", "gene-g0", "gene-g1", "0.95").split("\t")
        row[12] = "0.50"                                   # identity_cs
        rows = T.read_pairs([B.PAIRS_HEAD, "\t".join(row)])
        self.assertAlmostEqual(rows[0]["identity"], 0.95)

    def test_c10_run_atoms_verifies_the_script_hash_before_running_it(self):
        import unittest.mock as mock
        d = tempfile.mkdtemp(prefix="supp_hash_")
        bed = os.path.join(d, "x.bed")
        with open(bed, "w") as fh:
            fh.write(sedef_row("A", 0, 5000, "A", 20000, 25000) + "\n")
        cf = os.path.join(d, "c.txt")
        with open(cf, "w") as fh:
            fh.write("A\n")
        with mock.patch.object(T, "check_atoms_script", side_effect=RuntimeError("hash mismatch")) as chk:
            with self.assertRaises(RuntimeError):
                T.run_atoms(bed, cf, os.path.join(d, "o"))
            self.assertEqual(chk.call_count, 1)

    def test_w09_the_report_bootstrap_uses_the_registered_seed_and_b(self):
        d = tempfile.mkdtemp(prefix="supp_seed_")
        paths = world_with_dev_family(d)
        reg, rel = B.release_files(d, paths)
        out = os.path.join(d, "out")
        rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
        self.assertEqual(rc, 0, text)
        boot = load_json(os.path.join(out, "report.json"))["tau"]["0.90"]["bootstrap"]
        self.assertEqual((boot["B"], boot["seed"]), (2000, 20260930))


class RegisteredConstantsTests(unittest.TestCase):
    def test_registered_constants_and_positive_control_tau(self):
        self.assertEqual((T.TAU_PRIMARY, T.TAUS), (0.90, (0.90, 0.95, 0.98)))
        self.assertEqual((T.TOUCH_BP, T.CLASS_COVERAGE, T.MIN_PAIRS, T.MIN_FAMILIES, T.BOOTSTRAP_B, T.SEED), (100, 0.5, 20, 8, 2000, 20260930))
        self.assertEqual(T.REGISTERED_DEPTH_MATCHED, {"0.90": (38, 33), "0.95": (15, 15), "0.98": (5, 5)})

    def test_x05_the_positive_control_uses_the_tau_090_atoms(self):
        # copies are placed in an atom that exists only at tau 0.90 (SD row at fracMatch 0.92): the control must pass at 0.90 and would fail at 0.95
        d = tempfile.mkdtemp(prefix="supp_pc_")
        paths = B.make_world(d)
        names = ["LRPAP1", "L1", "L2", "L3", "L4", "L5", "L6", "L7"]
        with open(paths["genes"], "a") as fh:
            for i, n in enumerate(names):
                fh.write(f"NC_073224.2\t{n}\tgene\tprotein_coding\t{50201 + 300 * i}\t{50500 + 300 * i}\t{50200 + 300 * i}-{50500 + 300 * i}\n")   # inside the 0.92 atom [50000,56000)
        rows = [B.LIFTOFF_HEAD]
        for i, n in enumerate(names):
            if i == 0:
                rows.append(B.liftoff_row("gene-LRPAP1", "LRPAP1", "in_place", "NC_073224.2", 50200, 50500, "1.0", "-"))
            else:
                rows.append(B.liftoff_row("gene-LRPAP1", "LRPAP1", "dropped_M1", "NC_073224.2", 50200 + 300 * i, 50500 + 300 * i, "0.97", f"overlaps gene-{n} (NC_073224.2)"))
        liftoff = os.path.join(d, "liftoff.tsv")
        with open(liftoff, "w") as fh:
            fh.write("\n".join(rows) + "\n")
        out = os.path.join(d, "out")
        self.assertEqual(T.main(["counts", "--sedef", paths["sedef"], "--genes", paths["genes"], "--truth", paths["truth"], "--pairs", paths["pairs"], "--liftoff", liftoff,
                                 "--outdir", out]), 0)
        self.assertTrue(load_json(os.path.join(out, "counts.json"))["positive_control"]["passed"])


if __name__ == "__main__":
    unittest.main()
