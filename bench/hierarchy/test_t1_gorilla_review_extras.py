#!/usr/bin/env python3
"""Reviewer's extra tests for the survivors of the mutation run on the NEW release code (all pass on the frozen code; hermetic: made-up worlds under a temp
directory, `report` always with explicit made-up inputs and a scratch --outdir)."""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import t1_gorilla_pairs as T  # noqa: E402
import test_t1_gorilla_pairs as B  # noqa: E402

NAMES = ["sedef", "genes", "truth", "pairs", "liftoff", "dna_sd_atoms.py"]


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def f3(x):
    return "NA" if x is None else f"{x:.3f}"


class NearTokenTests(unittest.TestCase):
    def put(self, d, text):
        f = os.path.join(d, "rel")
        with open(f, "w", encoding="utf-8") as fh:
            fh.write(text)
        return f

    def test_only_the_exact_token_line_is_accepted(self):
        with tempfile.TemporaryDirectory() as d:
            T.require_release(self.put(d, T.RELEASE_TOKEN + "\n"))
            T.require_release(self.put(d, "  " + T.RELEASE_TOKEN + " \r\n"))
            for bad in (T.RELEASE_TOKEN[:-1] + "X\n", T.RELEASE_TOKEN[:-1] + "\n", T.RELEASE_TOKEN + "X\n", T.RELEASE_TOKEN.lower() + "\n", T.RELEASE_TOKEN.title() + "\n",
                        T.RELEASE_TOKEN.replace(" ", "  ", 1) + "\n", T.RELEASE_TOKEN + " " + T.RELEASE_TOKEN + "\n", "", "\n", T.RELEASE_TOKEN[:10] + "\n"):
                with self.assertRaises(T.ReleaseRefused, msg=repr(bad)):
                    T.require_release(self.put(d, bad))

    def test_a_near_token_is_refused_with_the_release_message_and_not_by_a_later_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            f = self.put(d, T.RELEASE_TOKEN[:-1] + "X\n")
            p = subprocess.run([sys.executable, "-B", os.path.join(HERE, "t1_gorilla_pairs.py"), "report", "--sedef", "/nonexistent/s", "--genes", "/nonexistent/g",
                                "--truth", "/nonexistent/t", "--pairs", "/nonexistent/p", "--liftoff", "/nonexistent/l", "--outdir", os.path.join(d, "out"),
                                "--release-file", f], capture_output=True, text=True)
            self.assertEqual(p.returncode, 2, p.stderr)
            self.assertIn("a release file with the registered token is required", p.stderr)
            self.assertFalse(os.path.exists(os.path.join(d, "out")))


class RegistryNameTests(unittest.TestCase):
    def reg(self, d, drop=None, short=None):
        path = os.path.join(d, "r.tsv")
        with open(path, "w") as fh:
            for m in NAMES:
                if m != drop:
                    fh.write(f"{m}\t{'a' * (15 if m == short else 16)}\n")
        return path

    def test_a_registry_lacking_any_one_input_is_refused_by_name(self):
        for n in NAMES:
            with tempfile.TemporaryDirectory() as d:
                with self.assertRaises(ValueError, msg=n) as cm:
                    T._registry_for(argparse.Namespace(registry=self.reg(d, drop=n)))
                self.assertTrue(str(cm.exception).startswith("registry "), str(cm.exception))
                self.assertIn(n, str(cm.exception))

    def test_a_short_prefix_for_any_one_input_is_refused_by_name(self):
        for n in NAMES:
            with tempfile.TemporaryDirectory() as d:
                with self.assertRaises(ValueError, msg=n) as cm:
                    T._registry_for(argparse.Namespace(registry=self.reg(d, short=n)))
                self.assertTrue(str(cm.exception).startswith("registry "), str(cm.exception))
                self.assertIn(n, str(cm.exception))

    def test_the_default_registry_passes_the_completeness_check(self):
        self.assertEqual(sorted(T._registry_for(argparse.Namespace(registry=None))), sorted(NAMES))

    def test_gate0_refuses_an_incomplete_registry_and_writes_nothing(self):
        d, paths = B.world_with_control()
        reg, _rel = B.release_files(d, paths)
        with open(reg, "w") as fh:
            fh.write("".join(f"{m}\t{'a' * 16}\n" for m in NAMES if m != "liftoff"))
        out = os.path.join(d, "out")
        rc, text = B.run_main(["gate0", "--sedef", paths["sedef"], "--genes", paths["genes"], "--truth", paths["truth"], "--pairs", paths["pairs"],
                               "--liftoff", paths["liftoff"], "--registry", reg, "--outdir", out])
        self.assertEqual(rc, 2, text)
        self.assertIn("registry", text.lower())
        self.assertFalse(os.path.exists(os.path.join(out, "gate0_r2.tsv")))

    def test_every_input_is_fingerprinted_against_its_own_prefix(self):
        for n in NAMES:
            d, paths = B.world_with_control()
            reg, rel = B.release_files(d, paths)
            with open(reg, "w") as fh:
                for m in NAMES:
                    pre = T.fingerprint(paths[m] if m in paths else T.ATOMS_SCRIPT)["sha256"][:16]
                    fh.write(f"{m}\t{'f' * 16 if m == n else pre}\n")
            out = os.path.join(d, "out")
            rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
            self.assertEqual(rc, 3, (n, text))
            self.assertIn(f"input {n}:", text)
            rep = load_json(os.path.join(out, "report.json"))
            self.assertFalse(rep["valid"])
            self.assertNotIn("tau", rep)


class PrintedReportTests(unittest.TestCase):
    def table(self, n_pairs, n_families, shared, **kw):
        t = {"tau": 0.9, "n_pairs": n_pairs, "n_families": n_families, "shared": shared, "pair_weighted": shared / n_pairs, "family_weighted": 0.7,
             "a_viol": 1 - shared / n_pairs, "wilson": [0.61, 0.87],
             "bootstrap": {"B": 2000, "seed": 20260930, "n_families": n_families, "pair_weighted": [0.58, 0.91], "family_weighted": [0.52, 0.93]},
             "b_viol": [5, 100], "b_viol_by_class": {"cross-contig": [3, 90], "same<100kb": [2, 10]}, "underpowered": False, "short_aligned_pairs": 7, "path_pairs": 23}
        t.update(kw)
        return t

    def test_golden_lines_of_one_block(self):
        lines = T.format_report({"valid": True, "tau": {"0.90": self.table(38, 33, 29)}})
        self.assertEqual(lines[:5], [
            "tau 0.90: 38 depth-matched pairs in 33 families (registered 38/33: match); UNDERPOWERED False",
            "  S-rate pair-weighted 0.763 (Wilson 0.610-0.870), family-weighted 0.700; A_viol 0.237",
            "  bootstrap (2000 family resamples, seed 20260930): pair-weighted 0.580-0.910, family-weighted 0.520-0.930",
            "  path pairs (a gene touches two or more classes): 23 of 38; pairs aligned < 0.5 of the shorter transcript: 7",
            "  B_viol: 5 of 100 eligible different-family pairs share a class (by distance class: cross-contig 3/90, same<100kb 2/10)"])
        self.assertEqual(len(lines), 6)
        self.assertIn("ceiling-type quantity", lines[5])

    def test_an_invalid_report_prints_no_rate_even_if_a_table_is_attached(self):
        lines = T.format_report({"valid": False, "reasons": ["why"], "tau": {"0.90": self.table(38, 33, 29)}})
        self.assertEqual(len(lines), 2)
        self.assertTrue(lines[0].startswith("INVALID"))
        self.assertEqual(lines[1], "  why")
        self.assertNotIn("0.763", "\n".join(lines))

    def test_the_printed_block_of_a_run_equals_report_json_line_for_line(self):
        d, paths = B.world_with_control()
        reg, rel = B.release_files(d, paths)
        out = os.path.join(d, "out")
        rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
        self.assertEqual(rc, 0, text)
        rep = load_json(os.path.join(out, "report.json"))
        printed = text.splitlines()
        self.assertTrue(printed[0].startswith("interpreter: python "), printed[0])
        self.assertNotIn("pair_rows", rep)
        for k, t in sorted(rep["tau"].items()):
            reg_tag = T.REGISTERED_DEPTH_MATCHED[k]
            exp = [f"tau {k}: {t['n_pairs']} depth-matched pairs in {t['n_families']} families (registered {reg_tag[0]}/{reg_tag[1]}: DIFFERS); UNDERPOWERED {t['underpowered']}",
                   f"  S-rate pair-weighted {f3(t['pair_weighted'])} (Wilson {f3(t['wilson'][0])}-{f3(t['wilson'][1])}), family-weighted {f3(t['family_weighted'])}; A_viol {f3(t['a_viol'])}"]
            b = t["bootstrap"]
            exp.append(f"  bootstrap ({b['B']} family resamples, seed {b['seed']}): pair-weighted {f3(b['pair_weighted'][0])}-{f3(b['pair_weighted'][1])}, "
                       f"family-weighted {f3(b['family_weighted'][0])}-{f3(b['family_weighted'][1])}" if b else "  bootstrap: NA (no families)")
            exp.append(f"  path pairs (a gene touches two or more classes): {t['path_pairs']} of {t['n_pairs']}; pairs aligned < 0.5 of the shorter transcript: {t['short_aligned_pairs']}")
            bc = ", ".join(f"{c} {v[0]}/{v[1]}" for c, v in sorted(t["b_viol_by_class"].items()))
            exp.append(f"  B_viol: {t['b_viol'][0]} of {t['b_viol'][1]} eligible different-family pairs share a class (by distance class: {bc or 'none'})")
            i = printed.index(exp[0])
            self.assertEqual(printed[i:i + 5], exp, k)


class ControlAndUniverseTauTests(unittest.TestCase):
    def world(self, d):
        """make_world with the eight LRPAP1 copies inside the atom [50000,56000) of the 0.92 SD row (an atom that exists at tau 0.90 only), plus family F20
        (identity 0.99) whose genes sit in that row's atoms: depth-matched at every tau by identity, eligible only at tau 0.90 by its atoms."""
        paths = B.make_world(d)
        names = ["LRPAP1", "L1", "L2", "L3", "L4", "L5", "L6", "L7"]
        extra_genes = ["NC_073224.2\tgX1\tgene\tprotein_coding\t54001\t54500\t54000-54500", "NC_073227.2\tgX2\tgene\tprotein_coding\t305001\t305500\t305000-305500"]
        with open(paths["genes"], "a") as fh:
            for i, n in enumerate(names):
                fh.write(f"NC_073224.2\t{n}\tgene\tprotein_coding\t{50201 + 300 * i}\t{50500 + 300 * i}\t{50200 + 300 * i}-{50500 + 300 * i}\n")
            fh.write("\n".join(extra_genes) + "\n")
        with open(paths["truth"], "a") as fh:
            fh.write("F20\twhole\t1to1\tgene-gX1\tNC_073224.2\t54001\t54500\nF20\twhole\t1to1\tgene-gX2\tNC_073227.2\t305001\t305500\n")
        with open(paths["pairs"], "a") as fh:
            fh.write(B.pairs_row("F20", "gene-gX1", "gene-gX2", "0.99") + "\n")
        rows = [B.LIFTOFF_HEAD]
        for i, n in enumerate(names):
            if i == 0:
                rows.append(B.liftoff_row("gene-LRPAP1", "LRPAP1", "in_place", "NC_073224.2", 50200, 50500, "1.0", "-"))
            else:
                rows.append(B.liftoff_row("gene-LRPAP1", "LRPAP1", "dropped_M1", "NC_073224.2", 50200 + 300 * i, 50500 + 300 * i, "0.97", f"overlaps gene-{n} (NC_073224.2)"))
        paths["liftoff"] = os.path.join(d, "liftoff.tsv")
        with open(paths["liftoff"], "w") as fh:
            fh.write("\n".join(rows) + "\n")
        return paths

    def test_the_positive_control_is_evaluated_at_tau_090_and_every_table_at_its_own_tau(self):
        d = tempfile.mkdtemp(prefix="extras_ctl_")
        paths = self.world(d)
        reg, rel = B.release_files(d, paths)
        out = os.path.join(d, "out")
        rc, text = B.run_main(B.report_argv(paths, reg, rel, out))
        self.assertEqual(rc, 0, text)
        rep = load_json(os.path.join(out, "report.json"))
        self.assertTrue(rep["valid"])
        self.assertTrue(rep["gates"]["positive_control"]["passed"])
        n = {k: rep["tau"][k]["n_pairs"] for k in ("0.90", "0.95", "0.98")}
        self.assertEqual(n, {"0.90": 3, "0.95": 1, "0.98": 0})       # F1 .97, F2 .92, F20 .99 at 0.90; F20's atoms do not exist at 0.95 and 0.98


class SurfaceTests(unittest.TestCase):
    def test_parser_defaults(self):
        ap = T.build_parser()
        self.assertEqual(ap.parse_args(["report"]).liftoff, T.DEFAULT_PATHS["liftoff"])
        self.assertIsNone(ap.parse_args(["counts"]).liftoff)
        self.assertIsNone(ap.parse_args(["gate0"]).liftoff)
        self.assertIsNone(ap.parse_args(["report"]).registry)
        self.assertIsNone(ap.parse_args(["gate0"]).registry)

    def test_counts_print_the_recon_line_once_and_its_differences_are_those_of_tau_090(self):
        d, paths = B.world_with_control()
        out = os.path.join(d, "out")
        rc, text = B.run_main(["counts", "--sedef", paths["sedef"], "--genes", paths["genes"], "--truth", paths["truth"], "--pairs", paths["pairs"],
                               "--liftoff", paths["liftoff"], "--outdir", out])
        self.assertEqual(rc, 0, text)
        t = load_json(os.path.join(out, "counts.json"))["tau"]["0.90"]
        self.assertEqual(text.count("registered recon values"), 1)
        self.assertIn(f"this run (identical tuples collapsed, coverage >= 0.5, as registered): {t['atoms']} atoms ({t['atoms'] - 40631:+d}), "
                      f"{t['classes']} classes ({t['classes'] - 12483:+d}), {t['edges_total']} atom edges ({t['edges_total'] - 372793:+d})", text)
        line = [ln for ln in text.splitlines() if "registered recon values" in ln][0]
        prev = text.splitlines()[text.splitlines().index(line) - 1]
        self.assertTrue(prev.startswith("tau 0.90:"), prev)

    def test_a_malformed_self_lift_table_gives_an_invalid_report_without_tables(self):
        # the reviewer's first version asserted the traceback of the frozen code; after its own recommendation 3 the run is INVALID with a recorded reason
        d, paths = B.world_with_control()
        reg, rel = B.release_files(d, paths)
        bad = os.path.join(d, "bad_liftoff.tsv")
        with open(bad, "w") as fh:
            fh.write("a\tb\n1\t2\n")
        out = os.path.join(d, "out")
        rc, text = B.run_main(B.report_argv(dict(paths, liftoff=bad), reg, rel, out))
        self.assertEqual(rc, 3, text)
        rep = load_json(os.path.join(out, "report.json"))
        self.assertFalse(rep["valid"])
        self.assertNotIn("tau", rep)
        self.assertEqual([f for f in os.listdir(out) if f.startswith("report_pairs")], [])


if __name__ == "__main__":
    unittest.main()
