#!/usr/bin/env python3
"""Tests of the loaders and of the gates of t1_refine.py on a tiny fixture (no real data).
Run: python3 -B -m unittest bench/hierarchy/test_t1_refine.py"""
import json
import os
import random
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import t1_lib as T  # noqa: E402
import t1_refine as R  # noqa: E402

S1C_HEAD = "Gene Name\tENSEMBL\tGene ID\tBiotype\tIn Table S1 (SD98 gene set)\tSD Unit\tMedian famCN\tNo. Assigned Families\tFamily ID\n"


def build_fixture(d):
    """Six genes on chr1: G1 and G2 clean (ID_1), G3 family-less, G4 multi-family, G5 clean without any duplicon, G6 clean and a folded annotation."""
    s1c = [
        ("A1", "G1", "protein_coding", "Yes", "10.0", "1", "ID_1"),
        ("A2", "G2", "protein_coding", "Yes", "12.0", "1", "ID_1"),
        ("A3", "G3", "lncRNA", "No", "N/A", "1", "Unassigned_3"),
        ("A4", "G4", "lncRNA", "Yes", "N/A", "2", "ID_1"),
        ("A4", "G4", "lncRNA", "Yes", "N/A", "2", "ID_2"),
        ("A5", "G5", "protein_coding", "Yes", "N/A", "1", "ID_2"),
        ("A6", "G6", "protein_coding", "Yes", "6.0", "1", "ID_3"),
    ]
    with open(f"{d}/s1c.tsv", "w") as fh:
        fh.write(S1C_HEAD)
        for name, g, bt, el, cn, nfam, fid in s1c:
            fh.write(f"{name}\tENS\t{g}\t{bt}\t{el}\tSD0\t{cn}\t{nfam}\t{fid}\n")
    genes = [
        dict(id="G1", n="A1", bt="protein_coding", k="p", b=0, c="chr1", s=100, e=400, st="+", x=[[100, 200], [300, 400]], sf=["ID_1"], cs=10.0, co=9.0),
        dict(id="G2", n="A2", bt="protein_coding", k="p", b=0, c="chr1", s=1000, e=1100, st="+", x=[[1000, 1100]], sf=["ID_1"], cs=12.0, co=11.0),
        dict(id="G3", n="A3", bt="lncRNA", k="n", b=0, c="chr1", s=5000, e=5100, st="-", x=[[5000, 5100]], sf=["Unassigned_3"], cs=None, co=None),
        dict(id="G4", n="A4", bt="lncRNA", k="n", b=0, c="chr1", s=7000, e=7100, st="+", x=[[7000, 7100]], sf=["ID_1", "ID_2"], cs=None, co=None),
        dict(id="G5", n="A5", bt="protein_coding", k="p", b=0, c="chr1", s=9000, e=9100, st="+", x=[[9000, 9100]], sf=["ID_2"], cs=None, co=None),
        dict(id="G6", n="A6", bt="protein_coding", k="p", b=0, c="chr1", s=11000, e=11100, st="+", x=[[11000, 11100]], sf=["ID_3"], cs=6.0, co=6.5),
    ]
    with open(f"{d}/families_cn.json", "w") as fh:
        json.dump(dict(genes=genes, links=[], sedef_links=[], families=[], combos=[]), fh)
    with open(f"{d}/exons.tsv", "w") as fh:
        fh.write("gene_id\tchrom\tstart\tend\n")
        for g in genes:
            for a, b in g["x"]:
                fh.write(f"{g['id']}\t{g['c']}\t{a}\t{b}\n")
        fh.write("OTHER\tchr1\t1\t2\n")                      # a gene that is not in S1C is ignored
    with open(f"{d}/colors.bed", "w") as fh:
        for a, b, lab in ((90, 210, "SD1"), (250, 420, "SD2"), (900, 1200, "SD3"), (4900, 5200, "SD4"), (6900, 7200, "SD5"), (10900, 11200, "SD6")):
            fh.write(f"chr1\t{a}\t{b}\t{lab}\t0\t+\t{a}\t{b}\t190,190,190\n")
    with open(f"{d}/split.tsv", "w") as fh:
        fh.write("gene_id\tunit_id\thalf\tn_unit_genes\tn_unit_families\tcross_half\n")
        for g in genes:
            fh.write(f"{g['id']}\tU0\t{'dev' if g['id'] in ('G1', 'G3') else 'heldout'}\t1\t1\t0\n")
    with open(f"{d}/rt.tsv", "w") as fh:
        fh.write("gene_id\tname\tbiotype\n")
        fh.write("G2\tA2\tprotein_coding\n")
    with open(f"{d}/regions.bed", "w") as fh:
        fh.write("chr1\t0\t1500\nchr1\t4000\t8000\n")
    fo = f"{d}/fo"
    os.makedirs(fo)
    with open(f"{fo}/names.tsv", "w") as fh:
        for g, loc in (("G1", "chr1:100-400"), ("G2", "chr1:1000-1100"), ("G3", "chr1:5000-5100"), ("G4", "chr1:7000-7100"),
                       ("G5", "chr1:9000-9100"), ("G6", "chr1:11000-11100")):
            fh.write(f"{g}\t{loc}\n")
    with open(f"{fo}/ours.clusters.tsv", "w") as fh:
        fh.write("cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n")
        for cid, loc in (("MCL0", (100, 400)), ("MCL0", (1000, 1100)), ("MCL1", (5000, 5100)), ("MCL1", (7000, 7100))):
            fh.write(f"{cid}\t2\t1.0\t0.5\tNA\tchr1\t{loc[0]}\t{loc[1]}\n")
    with open(f"{fo}/ours.loci.tsv", "w") as fh:
        fh.write("annotation\trepresentative\nchr1:11000-11100\tchr1:100-400\n")
    with open(f"{fo}/ours.graph.tsv", "w") as fh:
        fh.write("chr1:100-400\tchr1:1000-1100\t0.9\nchr1:5000-5100\tchr1:7000-7100\t0.8\n")
    return dict(s1c=f"{d}/s1c.tsv", families_cn=f"{d}/families_cn.json", sd98_exons=f"{d}/exons.tsv", colors=f"{d}/colors.bed",
                split=f"{d}/split.tsv", rt_genes=f"{d}/rt.tsv", sd98_regions=f"{d}/regions.bed", fo_dir=fo)


class UniverseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.paths = build_fixture(cls.tmp.name)
        cls.U = T.build_universe(cls.paths)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_classes_follow_the_number_of_assigned_soto_families(self):
        cls = {g: self.U.info[g]["cls"] for g in self.U.genes}
        self.assertEqual(cls, {"G1": "clean", "G2": "clean", "G3": "famless", "G4": "multi", "G5": "clean", "G6": "clean"})

    def test_V_drops_multi_family_genes_and_genes_without_an_exonic_duplicon(self):
        self.assertEqual(self.U.V, ["G1", "G2", "G3", "G6"])
        self.assertEqual(self.U.V_clean, ["G1", "G2", "G6"])
        self.assertEqual(self.U.V_famless, ["G3"])

    def test_strings_are_built_from_the_registered_exons(self):
        self.assertEqual(self.U.strings["G1"].u1, frozenset({"SD1", "SD2"}))
        self.assertEqual(self.U.strings["G1"].ordered, ("SD1", "SD2"))
        self.assertEqual(self.U.strings["G2"].u0, "SD3")

    def test_the_readthrough_list_and_eligibility_are_read(self):
        self.assertEqual(self.U.readthrough, {"G2"})
        self.assertTrue(self.U.info["G1"]["eligible"])
        self.assertFalse(self.U.info["G3"]["eligible"])

    def test_the_families_of_families_cn_are_kept_as_listed(self):
        self.assertEqual(self.U.info["G1"]["sf"], ("ID_1",))
        self.assertEqual(self.U.info["G4"]["sf"], ("ID_1", "ID_2"))
        self.assertEqual(self.U.info["G3"]["sf"], ("Unassigned_3",))

    def test_the_half_and_the_copy_numbers_are_read(self):
        self.assertEqual(self.U.info["G1"]["half"], "dev")
        self.assertEqual(self.U.info["G2"]["half"], "heldout")
        self.assertEqual((self.U.info["G1"]["cs"], self.U.info["G1"]["co"]), (10.0, 9.0))
        self.assertIsNone(self.U.info["G3"]["cs"])

    def test_the_family_of_a_gene_follows_rules_one_to_three(self):
        lab = {g: (self.U.fo_label[g], self.U.fo_rule[g]) for g in self.U.genes}
        self.assertEqual(lab["G1"], ("MCL0", 1))
        self.assertEqual(lab["G6"], ("MCL0", 2))          # a folded annotation takes the cluster of its representative
        self.assertEqual(lab["G3"], ("MCL1", 1))
        self.assertEqual(lab["G5"], (None, None))

    def test_V_F_is_V_with_a_cluster_and_nodes_follow_the_fold(self):
        self.assertEqual(self.U.VF, ["G1", "G2", "G3", "G6"])
        self.assertEqual(self.U.gene_node["G6"], "chr1:100-400")
        self.assertIsNone(self.U.gene_node["G5"])

    def test_strata_are_the_components_of_the_dumped_graph(self):
        self.assertEqual(self.U.strata, [("G1", ["G1", "G2", "G6"]), ("G3", ["G3"])])

    def test_the_region_of_a_gene_is_the_one_of_largest_exonic_overlap(self):
        self.assertEqual(self.U.region["G1"], ("chr1", 0, 1500))
        self.assertEqual(self.U.region["G3"], ("chr1", 4000, 8000))
        self.assertIsNone(self.U.region["G6"])

    def test_a_gene_missing_from_the_exon_table_is_an_error(self):
        with tempfile.TemporaryDirectory() as d:
            paths = build_fixture(d)
            with open(paths["sd98_exons"]) as fh:
                lines = fh.read().splitlines()
            with open(paths["sd98_exons"], "w") as fh:
                fh.write("\n".join(l for l in lines if not l.startswith("G3")) + "\n")
            with self.assertRaises(ValueError):
                T.build_universe(paths)


class ManifestTests(unittest.TestCase):
    def test_a_manifest_row_has_the_hash_size_and_line_count(self):
        with tempfile.TemporaryDirectory() as d:
            p = f"{d}/x.txt"
            with open(p, "w") as fh:
                fh.write("a\nb\nc\n")
            row = R.manifest_row(p)
            self.assertEqual((row["bytes"], row["lines"]), (6, 3))
            self.assertEqual(row["sha256"], T.sha256_file(p))
            self.assertEqual(len(row["sha256"]), 16)

    def test_a_big_file_gets_size_and_mtime_only_when_hashing_is_off(self):
        with tempfile.TemporaryDirectory() as d:
            p = f"{d}/x.txt"
            with open(p, "w") as fh:
                fh.write("abc")
            row = R.manifest_row(p, hash_it=False)
            self.assertIsNone(row["sha256"])
            self.assertEqual(row["bytes"], 3)
            self.assertIn("mtime", row)


class TripwireTests(unittest.TestCase):
    def test_a_count_that_differs_from_its_registered_value_is_a_mismatch(self):
        rows = R.compare_registered({"a": 5, "b": 3}, {"a": 5, "b": 4})
        self.assertEqual({r["name"]: r["ok"] for r in rows}, {"a": True, "b": False})
        self.assertFalse(R.all_ok(rows))
        self.assertTrue(R.all_ok(rows[:1]))

    def test_a_missing_obtained_value_is_a_mismatch(self):
        rows = R.compare_registered({"a": 5}, {})
        self.assertFalse(rows[0]["ok"])


class DenominatorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.paths = build_fixture(cls.tmp.name)
        cls.U = T.build_universe(cls.paths)
        cls.den = R.denominators(cls.U)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_universe_sizes(self):
        d = self.den
        self.assertEqual((d["genes_total"], d["class_clean"], d["class_multi"], d["class_famless"]), (6, 4, 1, 1))
        self.assertEqual((d["genes_with_exonic_duplicon"], d["V"], d["V_clean"], d["V_famless"]), (5, 4, 3, 1))
        self.assertEqual(d["non_multi_without_duplicon"], 1)

    def test_exon_counts_over_all_genes_and_over_V(self):
        d = self.den
        self.assertEqual((d["merged_exons_all"], d["labelled_exons_all"], d["merged_exons_V"], d["labelled_exons_V"]), (7, 6, 5, 5))
        self.assertEqual((d["tie_exons_all"], d["tie_exons_V"], d["tie_genes_V"]), (0, 0, 0))

    def test_string_marginals(self):
        d = self.den
        self.assertEqual((d["U1_distinct_V"], d["U1_distinct_V_union"], d["U1_distinct_V_num"]), (4, 4, 4))
        self.assertEqual((d["U0_distinct_V"], d["U1_distinct_clean"]), (4, 3))

    def test_the_f_o_layer_counts(self):
        d = self.den
        self.assertEqual((d["fo_clusters"], d["fo_rows"], d["fo_size2_rows"]), (2, 4, 2))
        self.assertEqual((d["VF"], d["VF_rule1"], d["VF_rule2"], d["V_unclustered"]), (4, 3, 1, 0))
        self.assertEqual((d["graph_components"], d["strata_with_VF"], d["strata_variation"], d["strata_single"]), (2, 2, 0, 2))

    def test_exclusion_arms(self):
        d = self.den
        self.assertEqual((d["readthrough_V"], d["readthrough_all"]), (1, 1))
        self.assertEqual((d["non_eligible_V"], d["non_eligible_all"]), (1, 1))
        self.assertEqual((d["span_pairs_V"], d["span_genes_V"], d["span_pairs_all"], d["span_genes_all"]), (0, 0, 0, 0))
        self.assertEqual((d["short_V"], d["short_all"]), (4, 6))

    def test_the_multi_family_genes_are_counted_outside_V(self):
        d = self.den
        self.assertEqual(d["multi_family_genes"], 1)
        self.assertEqual(d["fo_clusters_with_multi_family_genes"], 1)
        self.assertEqual((d["fo_members_not_in_V_multi"], d["fo_members_not_in_V_no_duplicon"]), (1, 0))

    def test_size_two_and_small_clusters_are_counted_in_the_layer_universe(self):
        d = self.den
        self.assertEqual((d["fo_clusters_touching_VF"], d["fo_clusters_ge2_VF"], d["fo_clusters_lt2_VF"], d["fo_clusters_empty_VF"]), (2, 1, 1, 0))
        self.assertEqual((d["fo_size2_VF"], d["fo_size2_both"], d["fo_size2_rows"]), (0, 0, 2))

    def test_folded_annotations_and_unconnected_genes(self):
        d = self.den
        self.assertEqual((d["fo_folded_records"], d["fo_folded_representatives"]), (1, 1))
        self.assertEqual((d["V_folded"], d["V_folded_rep_not_row"]), (1, 0))
        self.assertEqual((d["V_never_connected"], d["V_singleton_nodes"]), (0, 0))
        self.assertEqual((d["graph_nodes"], d["graph_edges"]), (4, 2))

    def test_halves_regions_and_copy_numbers(self):
        d = self.den
        self.assertEqual((d["V_dev"], d["V_heldout"]), (2, 2))
        self.assertEqual((d["fo_clusters_straddling_halves"], d["fo_genes_in_straddling_clusters"], d["fo_genes_in_clusters_ge2"]), (1, 3, 3))
        self.assertEqual(d["regions_V"], 2)
        self.assertEqual((d["V_no_famcn_cs"], d["V_no_famcn_co"], d["V_no_famcn_either"]), (1, 1, 1))

    def test_soto_families_in_the_clean_universe(self):
        d = self.den
        self.assertEqual((d["fs_families_V_clean"], d["fs_families_V_clean_ge2"]), (2, 1))
        self.assertEqual((d["fs_families_all_clean"], d["fs_families_all_clean_ge2"]), (3, 1))

    def test_genes_with_two_or_more_dominant_duplicons(self):
        d = self.den
        self.assertEqual((d["genes_u1_ge2_all"], d["genes_u1_ge2_V"]), (1, 1))

    def test_tie_and_sensitivity_counts_are_zero_on_a_fixture_without_ties(self):
        d = self.den
        self.assertEqual((d["num_diff_exons_V"], d["num_diff_genes_V"], d["num_u1_changed_V"]), (0, 0, 0))

    def test_every_value_is_a_plain_number_so_the_table_is_deterministic(self):
        for k, v in self.den.items():
            self.assertIsInstance(v, (int, float, str, list, dict, tuple, type(None)), k)


class StringsTableTests(unittest.TestCase):
    def test_the_table_has_strings_only_and_is_sorted_by_gene(self):
        with tempfile.TemporaryDirectory() as d:
            U = T.build_universe(build_fixture(d))
            text = R.strings_table(U)
        lines = text.strip().split("\n")
        self.assertEqual(lines[0], "gene\tu0\tu1\tordered\tn_exons")
        self.assertEqual([l.split("\t")[0] for l in lines[1:]], ["G1", "G2", "G3", "G6"])
        self.assertNotIn("MCL", text)
        self.assertNotIn("ID_", text)
        self.assertEqual(lines[1], "G1\tSD1\tSD1,SD2\tSD1>SD2\t2")


class DeterminismTests(unittest.TestCase):
    def test_identical_directories_agree_and_a_changed_file_does_not(self):
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            for d in (a, b):
                with open(f"{d}/x.tsv", "w") as fh:
                    fh.write("1\n")
            self.assertEqual(R.compare_dirs(a, b, ["x.tsv"]), [])
            with open(f"{b}/x.tsv", "w") as fh:
                fh.write("2\n")
            self.assertEqual(R.compare_dirs(a, b, ["x.tsv"]), ["x.tsv"])
            self.assertEqual(R.compare_dirs(a, b, ["missing.tsv"]), ["missing.tsv"])


class GuardTests(unittest.TestCase):
    def test_the_decisive_cells_refuse_real_data_without_the_go_file(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(SystemExit):
                R.require_go(d)
            with open(f"{d}/GO_R1", "w") as fh:
                fh.write("go\n")
            R.require_go(d)                                      # no exception now

    def test_the_decisive_cells_run_on_a_synthetic_world(self):
        out = R.run_decisive_synthetic(seed=0, B=100)
        self.assertEqual(set(out), {"t1a", "t1b", "t1c", "t1e"})
        self.assertEqual(out["t1a"]["verdict"], "REFUTED")
        self.assertIn(out["t1b"]["class"], ("CUT-ALIGNED", "PARTIAL", "NOT-DISTINGUISHED", "UNDERPOWERED", "UNINFORMATIVE"))
        self.assertGreater(out["t1e"]["n_pairs"], 0)


class GateFourTests(unittest.TestCase):
    def test_every_deterministic_synthetic_clause_passes(self):
        rows = R.gate4_det(B=200)
        bad = [r for r in rows if not r["ok"]]
        self.assertEqual(bad, [], bad)
        self.assertGreaterEqual(len(rows), 12)

    def test_the_calibration_block_reports_the_share_of_small_p(self):
        res = R.gate4_calib_block(0, 2, B=100, nrep=2)
        self.assertEqual(len(res["m1c"]), 2)
        self.assertTrue(all(0 < p <= 1 for p in res["m1c"] + res["m2"]))


# ---------------------------------------------------------------- stage 2: worlds, cells, Gate 1 helpers, Gate 5, Gate 3
class WorldAndCellTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.U = T.build_universe(build_fixture(cls.tmp.name))
        cls.w = R.world_of_universe(cls.U)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_the_world_is_V_F_with_its_labels_strings_strata_and_halves(self):
        w = self.w
        self.assertEqual(w["genes"], ["G1", "G2", "G3", "G6"])
        self.assertEqual(w["fam"], {"G1": "MCL0", "G2": "MCL0", "G3": "MCL1", "G6": "MCL0"})
        self.assertEqual(w["u1"]["G1"], frozenset({"SD1", "SD2"}))
        self.assertEqual(w["ordered"]["G1"], ("SD1", "SD2"))
        self.assertEqual(w["strata"], [("G1", ["G1", "G2", "G6"]), ("G3", ["G3"])])
        self.assertEqual(w["half"]["G2"], "heldout")

    def test_an_arm_deletes_genes_from_the_universe_and_the_strata(self):
        w = R.world_of_universe(self.U, delete={"G2"})
        self.assertEqual(w["genes"], ["G1", "G3", "G6"])
        self.assertEqual(w["strata"], [("G1", ["G1", "G6"]), ("G3", ["G3"])])

    def test_a_half_keeps_only_the_genes_of_that_half(self):
        w = R.world_of_universe(self.U, half="dev")
        self.assertEqual(w["genes"], ["G1", "G3"])

    def test_t1a_counts_the_exceptions_and_refutes_on_a_single_one(self):
        r = R.t1a_report(self.w)
        self.assertEqual((r["n"], r["m1c"], r["m2"], r["m1"]), (4, 2, 0, 2))
        self.assertEqual(r["verdict"], "REFUTED")
        self.assertEqual(len(r["exceptions_d1"]), r["m1c"])
        self.assertEqual(len(r["exceptions_d2"]), r["m2"])
        self.assertEqual(r["verdict_scope"], "primary")

    def test_t1a_is_not_refuted_only_when_both_counts_are_zero(self):
        w = dict(self.w, genes=["G3"], fam={"G3": "MCL1"}, strata=[("G3", ["G3"])])
        r = R.t1a_report(w)
        self.assertEqual((r["m1c"], r["m2"], r["verdict"]), (0, 0, "NOT-REFUTED"))

    def test_t1a_reports_the_size_two_stratum_and_the_macro_versions(self):
        w = dict(self.w, genes=["G1", "G2"], fam={"G1": "A", "G2": "A"}, strata=[("G1", ["G1", "G2"])])
        r = R.t1a_report(w)
        self.assertEqual(r["size2"]["n"], 2)
        self.assertEqual(r["size2"]["m1c"], 1)
        self.assertAlmostEqual(r["macro_d1"], 0.5)

    def test_the_cell_over_a_world_without_variation_bearing_strata_is_underpowered_and_prints_the_counts_first(self):
        cell = R.t1bc_report(self.w, "m1c", B=50)
        self.assertEqual((cell["class"], cell["n_var"], cell["n_live"]), ("UNDERPOWERED", 0, 0))
        text = R.format_cell(cell)
        self.assertLess(text.index("variation-bearing strata"), text.index("p = "))
        self.assertLess(text.index("live strata"), text.index("p = "))

    def test_t1e_over_a_world_without_pairs_reports_zero_pairs(self):
        self.assertEqual(R.t1e_report(self.w, B=20)["n_pairs"], 0)


class ArmTests(unittest.TestCase):
    def test_the_arms_delete_genes_without_looking_at_any_family_label(self):
        with tempfile.TemporaryDirectory() as d:
            U = T.build_universe(build_fixture(d))
        arms = R.arm_deletions(U)
        self.assertEqual(set(arms), {"none", "readthrough", "non_eligible", "same_locus", "all"})
        self.assertEqual(arms["none"], set())
        self.assertEqual(arms["readthrough"], {"G2"})
        self.assertEqual(arms["non_eligible"], {"G3"})
        self.assertEqual(arms["all"], {"G2", "G3"})

    def test_the_denominator_guard_flags_an_arm_that_removes_more_than_ten_percent_of_V(self):
        self.assertTrue(R.denominator_guard(215, 2142))
        self.assertFalse(R.denominator_guard(214, 2142))


class GateFiveTests(unittest.TestCase):
    def test_the_scorer_reproduces_the_hand_worked_partition_deficiency_with_singleton_and_ignored_ungrouped_genes(self):
        clean = {"a": "F1", "b": "F1", "c": "F1", "d": "F2", "e": "F2", "f": "F3"}
        Q = {"a": "Q1", "b": "Q1", "c": "Q2", "d": "Q3", "e": "Q3"}              # c is split off F1, f is ungrouped
        v = R.gate5b_values(clean, Q)
        self.assertEqual((v["D_S_to_Q_singleton"], v["S_split_singleton"], v["S_scored_singleton"]), (1, 1, 2))
        self.assertEqual((v["D_S_to_Q_ignored"], v["S_split_ignored"]), (1, 1))

    def test_monotone_u1_to_u2_holds_on_random_partitions_of_made_up_strings(self):
        rng = random.Random(3)
        genes = [f"g{i:03d}" for i in range(120)]
        u1 = {g: frozenset(rng.sample("ABCDEF", rng.randint(1, 3))) for g in genes}
        cn = {g: float(rng.randint(2, 30)) for g in genes}
        rows = R.gate5a_surrogate(u1, cn, seeds=(0, 1), source="toy")
        self.assertGreaterEqual(len(rows), 12)
        self.assertTrue(all(r["ok"] for r in rows), [r for r in rows if not r["ok"]])
        self.assertTrue(all(r["u2_refines_u1"] for r in rows))

    def test_a_broken_scorer_would_be_caught(self):
        rng = random.Random(5)
        genes = [f"g{i:03d}" for i in range(60)]
        u1 = {g: frozenset(rng.sample("ABCD", rng.randint(1, 2))) for g in genes}
        cn = {g: float(rng.randint(2, 30)) for g in genes}
        rows = R.gate5a_surrogate(u1, cn, seeds=(0, 1), source="toy", _break=True)
        self.assertFalse(all(r["ok"] for r in rows))


class UnitCoverTableTests(unittest.TestCase):
    TABLE = """# KEY=unitcovercn
| arm | famCN | units | Jaccard dev / held-out | p dev / held-out | exact | recall | precision | clean: own family / 2+ families / empty |
|---|---|---|---|---|---|---|---|---|
| A | - | 1546 | 0.439 / 0.473 | 0.1667 / 0.1667 | 0.107 | 0.555 | 0.718 | 0.546 / 0.303 / 0.020 |
| C | ours | 2400 | 0.492 / 0.601 | 0.1667 / 0.1667 | 0.161 | 0.817 | 0.608 | 0.626 / 0.225 / 0.053 |
"""

    def test_the_frozen_table_is_parsed_per_arm_and_famcn_source(self):
        t = R.parse_unitcover_table(self.TABLE)
        self.assertEqual(set(t), {("A", "-"), ("C", "ours")})
        self.assertEqual(t[("A", "-")], dict(units=1546, dev=0.439, heldout=0.473, exact=0.107, recall=0.555, precision=0.718,
                                              own=0.546, multi=0.303, empty=0.020))

    def test_numbers_are_compared_at_the_three_decimals_the_frozen_script_prints(self):
        mine = dict(n_units=1546, obs={"dev": (0.43949, 10), "heldout": (0.47349, 10)}, exact=0.1074, recall=0.5549, precision=0.7181,
                    clean_ok=0.5462, clean_multi=0.3031, clean_empty=0.0199)
        row = R.unit_cover_row(mine)
        self.assertEqual(row, dict(units=1546, dev=0.439, heldout=0.473, exact=0.107, recall=0.555, precision=0.718, own=0.546, multi=0.303, empty=0.020))


class GateThreeTests(unittest.TestCase):
    def test_the_deterministic_tables_do_not_depend_on_the_call_and_hold_no_family_label_next_to_a_string(self):
        with tempfile.TemporaryDirectory() as d:
            U = T.build_universe(build_fixture(d))
        a = R.gate3_tables(U, sizes=R.default_sizes(8), B=60)
        b = R.gate3_tables(U, sizes=R.default_sizes(8), B=60)
        self.assertEqual(a, b)
        self.assertEqual(set(a), {"strings_table.tsv", "denominators.json", "gate4_det.json", "synthetic_cells.json", "synthetic_t1e.json", "exceptions.json"})
        for line in a["strings_table.tsv"].splitlines():
            self.assertNotIn("MCL", line)

    def test_the_canary_lists_a_set_of_strings_in_hash_order_so_the_two_seeds_can_be_told_apart(self):
        c = R.hash_canary()
        self.assertEqual(sorted(c.split()), sorted(R.CANARY_WORDS))
        self.assertEqual(len(R.CANARY_WORDS), len(set(R.CANARY_WORDS)))
        self.assertGreaterEqual(len(R.CANARY_WORDS), 20)

    def test_the_comparison_is_by_file_content(self):
        with tempfile.TemporaryDirectory() as x, tempfile.TemporaryDirectory() as y:
            for d in (x, y):
                for n, t in (("a.json", "1"), ("b.tsv", "2")):
                    with open(f"{d}/{n}", "w") as fh:
                        fh.write(t)
            self.assertEqual(R.compare_dirs(x, y, ["a.json", "b.tsv"]), [])


class GateSwitchTests(unittest.TestCase):
    def test_a_cell_needs_every_covering_gate_to_pass(self):
        with tempfile.TemporaryDirectory() as d:
            for g in R.GATES_FOR_R1:
                with open(f"{d}/{g}.json", "w") as fh:
                    json.dump({"status": "PASS"}, fh)
            R.require_gates(d)
            with open(f"{d}/gate4.json", "w") as fh:
                json.dump({"status": "FAIL"}, fh)
            with self.assertRaises(SystemExit):
                R.require_gates(d)
            os.remove(f"{d}/gate4.json")
            with self.assertRaises(SystemExit):
                R.require_gates(d)

    def test_the_r0_summary_lists_every_gate_with_its_status(self):
        with tempfile.TemporaryDirectory() as d:
            for g, st in (("gate0", "PASS"), ("gate1", "PASS"), ("gate2", "SKIPPED"), ("gate3", "PASS"), ("gate4", "FAIL"), ("gate5", "PASS")):
                with open(f"{d}/{g}.json", "w") as fh:
                    json.dump({"status": st}, fh)
            text, ok = R.r0_summary(d)
        self.assertFalse(ok)
        self.assertIn("gate4: FAIL", text)
        self.assertIn("gate2: SKIPPED", text)
        self.assertIn("gate6: MISSING", text)

    def test_the_parser_knows_every_subcommand_of_the_directive(self):
        ap = R.build_parser()
        for cmd in ("inputs", "hash-paf", "gate0", "gate1", "gate2", "gate3", "gate4", "gate5", "gate6", "r0", "t1a", "t1b", "t1c", "t1e"):
            ap.parse_args([cmd])
        a = ap.parse_args(["gate4", "--part", "calib", "--start", "0", "--end", "20"])
        self.assertEqual((a.part, a.start, a.end), ("calib", 0, 20))
        a = ap.parse_args(["t1b", "--synthetic", "--B", "100"])
        self.assertTrue(a.synthetic)
        self.assertEqual(a.B, 100)

    def test_a_real_cell_is_refused_without_the_go_file(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(SystemExit):
                R.guard_real_cell(run_dir=d, out_dir=d)


# ---------------------------------------------------------------- stage 3: layers, sensitivities, sensitivity nulls, the level sweep
class LayerAndNullTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.U = T.build_universe(build_fixture(cls.tmp.name))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_strata_by_a_key_are_ordered_by_their_smallest_gene_id(self):
        strata = R.strata_by(["g3", "g1", "g2", "g4"], {"g1": "x", "g2": "y", "g3": "x", "g4": "y"}.get)
        self.assertEqual(strata, [("g1", ["g1", "g3"]), ("g2", ["g2", "g4"])])

    def test_F_S_is_the_clean_universe_with_the_soto_family_and_the_F_Q_groups_as_strata(self):
        w = R.layer_world(self.U, "F_S", {"G1": "q1", "G2": "q1"})
        self.assertEqual(w["genes"], ["G1", "G2", "G6"])
        self.assertEqual(w["fam"], {"G1": "ID_1", "G2": "ID_1", "G6": "ID_3"})
        self.assertEqual(w["strata"], [("G1", ["G1", "G2"]), ("G6", ["G6"])])           # G6 is in no F_Q group: the stratum 'no group'

    def test_genes_in_no_F_Q_group_share_one_stratum(self):
        w = R.layer_world(self.U, "F_S", {"G1": "q1"})
        self.assertEqual(w["strata"], [("G1", ["G1"]), ("G2", ["G2", "G6"])])

    def test_F_Q_is_V_with_ungrouped_genes_as_singleton_families(self):
        w = R.layer_world(self.U, "F_Q", {"G1": "q1", "G2": "q1"})
        self.assertEqual(w["genes"], ["G1", "G2", "G3", "G6"])
        self.assertEqual(w["fam"]["G1"], "q1")
        self.assertEqual(w["fam"]["G3"], "single:G3")

    def test_chromosome_and_region_strata_are_taken_from_the_universe(self):
        w = R.world_of_universe(self.U)
        self.assertEqual(R.null_strata(self.U, w, "chromosome"), [("G1", ["G1", "G2", "G3", "G6"])])
        reg = R.null_strata(self.U, w, "region")
        self.assertEqual(reg, [("G1", ["G1", "G2"]), ("G3", ["G3"]), ("G6", ["G6"])])           # G6 has no region and is alone in the stratum 'no region'
        self.assertEqual(R.null_strata(self.U, w, "governing"), w["strata"])

    def test_the_sensitivity_strings_are_rebuilt_from_the_exons(self):
        sens = R.sensitivity_u1(self.U, ["G1", "G2"])
        self.assertEqual(set(sens), {"union_per_label", "numeric_ties", "tie_set", "near_tie_band", "any_overlap"})
        self.assertEqual(sens["any_overlap"]["G1"], frozenset({"SD1", "SD2"}))
        self.assertEqual(sens["union_per_label"]["G2"], frozenset({"SD3"}))

    def test_the_rungs_u0_and_u2_are_added_to_a_world_without_touching_the_family_labels(self):
        w = R.world_of_universe(self.U)
        r = R.rung_strings(self.U, w)
        self.assertEqual(r["u0"]["G1"], frozenset({"SD1"}))
        self.assertEqual(set(r["u2co"]), {"G1", "G2", "G6"})              # G3 has no copy number
        self.assertEqual(r["u2co"]["G1"], frozenset({("SD1", 0), ("SD2", 0)}))
        self.assertEqual(r["ordered_u2co"]["G1"], (("SD1", 0), ("SD2", 0)))

    def test_the_sweep_directory_is_named_after_the_inflation(self):
        self.assertTrue(R.sweep_dir(2.0).endswith("/sweep_I2.0"))
        self.assertTrue(R.sweep_dir(4.0).endswith("/sweep_I4.0"))

    def test_t1a_adds_the_audit_records_of_every_family_and_string(self):
        w = R.world_of_universe(self.U)
        r = R.t1a_report(w)
        self.assertEqual([a["family"] for a in r["audit_d1"]], ["MCL0"])
        self.assertEqual(r["audit_d1"][0]["n"], 3)
        self.assertEqual(r["audit_d1"][0]["W"], 1)
        self.assertEqual(r["audit_d2"], [])

    def test_the_singleton_arm_counts_unclustered_genes_sharing_a_string_with_a_clustered_gene(self):
        rep = R.singleton_arm(self.U)
        self.assertEqual((rep["unclustered"], rep["unclustered_sharing_a_string_with_a_clustered_gene"], rep["m2"]), (0, 0, 0))

    def test_a_sensitivity_cell_has_no_class_and_runs_on_a_world(self):
        w = R.world_of_universe(self.U)
        r = R.sensitivity_null_report(self.U, w, "chromosome", B=30)
        self.assertEqual(r["null"], "chromosome")
        self.assertNotIn("class", r["m1c"])


class BuilderTests(unittest.TestCase):
    """The report builders behind t1a, t1b, t1c and t1e, run on the made-up fixture (never on real labels)."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.U = T.build_universe(build_fixture(cls.tmp.name))
        cls.Q = {"G1": "q1", "G2": "q1"}

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_t1a_builds_the_arms_halves_rungs_sensitivities_and_layers(self):
        res, text = R.build_t1a(self.U, self.Q)
        self.assertEqual(set(res["arms"]), {"none", "readthrough", "non_eligible", "same_locus", "all"})
        self.assertEqual(set(res["halves"]), {"dev", "heldout"})
        self.assertEqual(res["arms"]["none"]["verdict"], "REFUTED")
        self.assertEqual(res["arms"]["none"]["verdict_scope"], "primary")
        self.assertEqual(res["arms"]["readthrough"]["verdict_scope"], "descriptive")
        self.assertEqual(res["arms"]["readthrough"]["deleted"], 1)
        self.assertEqual(set(res["rungs"]), {"U0", "U2CO", "U2CS"})
        self.assertEqual(res["rungs"]["U2CO"]["genes_without_famcn"], 1)
        self.assertEqual(set(res["sensitivities"]), {"union_per_label", "numeric_ties", "tie_set", "near_tie_band", "any_overlap"})
        self.assertEqual(set(res["layers"]), {"F_Q", "F_S"})
        self.assertIn("m1c (D1) = 2", text)

    def test_the_non_primary_arms_never_issue_not_refuted(self):
        res, _ = R.build_t1a(self.U, self.Q)
        for arm in ("readthrough", "non_eligible", "same_locus", "all"):
            self.assertNotEqual(res["arms"][arm]["verdict"], "NOT-REFUTED")

    def test_t1b_and_t1c_run_on_the_fixture_and_print_the_counts_first(self):
        for stat in ("m1c", "m2"):
            cell, text = R.build_t1bc(self.U, stat, B=20)
            self.assertEqual(cell["class"], "UNDERPOWERED")
            self.assertLess(text.index("live strata"), text.index("p = "))

    def test_a_sensitivity_null_is_reported_without_a_class(self):
        for kind in ("chromosome", "region"):
            res, text = R.build_t1bc(self.U, "m1c", B=20, null=kind)
            self.assertEqual(res["null"], kind)
            self.assertNotIn("class", res["m1c"])
            self.assertIn("no class", text)

    def test_the_descriptive_F_S_cell_uses_the_F_Q_groups_as_strata_and_has_no_class(self):
        res, text = R.build_t1bc(self.U, "m1c", B=20, layer="F_S", Q=self.Q)
        self.assertEqual(res["layer"], "F_S")
        self.assertNotIn("class", res["m1c"])
        self.assertIn("F_Q", text)

    def test_t1e_builds_on_the_fixture(self):
        res, text = R.build_t1e(self.U, B=20)
        self.assertEqual(res["n_pairs"], 0)
        self.assertIn("T1-e", text)

    def test_the_level_specific_tag_compares_content_classes_only(self):
        tag = R.level_specific_tag
        self.assertEqual(tag("CUT-ALIGNED", "CUT-ALIGNED"), "")
        self.assertEqual(tag("CUT-ALIGNED", "PARTIAL"), "LEVEL-SPECIFIC")
        self.assertEqual(tag("PARTIAL", "NOT-DISTINGUISHED"), "LEVEL-SPECIFIC")
        self.assertEqual(tag("CUT-ALIGNED", "UNDERPOWERED"), "UNDERPOWERED (no LEVEL-SPECIFIC comparison)")
        self.assertEqual(tag("CUT-ALIGNED", "UNINFORMATIVE"), "UNINFORMATIVE (no LEVEL-SPECIFIC comparison)")


if __name__ == "__main__":
    unittest.main()
