#!/usr/bin/env python3
"""Tests of the gorilla SEDEF arm (Phase R2, T1-d) of docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md section 5 and Amendment 1.

Run: python3 -B -m unittest bench/hierarchy/test_t1_gorilla_pairs.py

Written BEFORE the implementation (each group was run and seen to fail first). The synthetic cases here are the Gate 4 style controls of
the arm: they exercise the scorers on made-up tables and never touch the real verdict set.
"""
import collections
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)

ATOMS = os.path.join(REPO, "bench", "dna_sd_atoms.py")
ATOMS_SHA256 = "9d96d976eb52dca6e056b1bc8256f006835d72b265a8697c925d15cfcb562bce"


def sedef_row(c1, s1, e1, c2, s2, e2, strand2="+", cigar=None, frac=0.95, aln_len=None, extra=""):
    """A 34-column SEDEF line as the gorilla table has it (0-based BED coordinates). Columns used by the code under test:
    1-6 coordinates, 9-10 strands, 12 aln_len, 17-18 matchB / mismatchB, 21 fracMatch, 33 cigar."""
    n = e1 - s1
    cigar = cigar if cigar is not None else f"{n}M"
    aln_len = aln_len if aln_len is not None else n
    matches = int(round(frac * 1000))
    f = [""] * 34
    f[0], f[1], f[2], f[3], f[4], f[5] = c1, str(s1), str(e1), c2, str(s2), str(e2)
    f[6], f[7], f[8], f[9] = "S", "5.2", "+", strand2
    f[10], f[11], f[12] = str(n), str(aln_len), "m=5.2;g=0.1" + extra
    for i in range(13, 32):
        f[i] = "0"
    f[16], f[17] = str(matches), str(1000 - matches)
    f[20], f[21] = f"{frac:.6f}", f"{frac:.6f}"
    f[32], f[33] = cigar, f"{frac:.6f}"
    return "\t".join(f)


def read_text(path):
    with open(path) as fh:
        return fh.read()


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


SEDEF_HEADER = "#chr1\tstart1\tend1\tchr2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\tmax_len\taln_len\tcomment\t" + "\t".join(f"c{i}" for i in range(14, 35))


class RestoredAtomsScriptTests(unittest.TestCase):
    """bench/dna_sd_atoms.py is restored unchanged from the attic (git b509675c)."""

    def test_the_restored_script_is_byte_identical_to_the_attic_copy(self):
        self.assertTrue(os.path.exists(ATOMS), "bench/dna_sd_atoms.py is not restored")
        with open(ATOMS, "rb") as fh:
            self.assertEqual(hashlib.sha256(fh.read()).hexdigest(), ATOMS_SHA256)

    def run_atoms(self, rows, contigs, mode="cigar"):
        d = tempfile.mkdtemp(prefix="r2_atoms_")
        sd = os.path.join(d, "sd.bed")
        with open(sd, "w") as fh:
            fh.write("\n".join(rows) + "\n")
        cf = os.path.join(d, "contigs.txt")
        with open(cf, "w") as fh:
            fh.write("\n".join(contigs) + "\n")
        out = os.path.join(d, "out")
        p = subprocess.run([sys.executable, "-B", ATOMS, sd, "gorilla", "@" + cf, out, mode], capture_output=True, text=True)
        return p, out

    def read(self, path):
        with open(path) as fh:
            return [ln.rstrip("\n").split("\t") for ln in fh][1:]

    def test_cigar_mode_makes_two_atoms_and_one_full_coverage_edge_for_one_pair_and_skips_the_header_line(self):
        rows = [sedef_row("A", 0, 5000, "A", 20000, 25000), SEDEF_HEADER]
        p, out = self.run_atoms(rows, ["A"])
        self.assertEqual(p.returncode, 0, p.stderr)
        nodes = self.read(out + ".nodes.tsv")
        self.assertEqual([(n[1], n[2], n[3]) for n in nodes], [("A", "0", "5000"), ("A", "20000", "25000")])
        edges = self.read(out + ".edges.tsv")
        self.assertEqual(len(edges), 1)
        self.assertEqual((edges[0][0], edges[0][1]), ("0", "1"))
        self.assertAlmostEqual(float(edges[0][2]), 0.95, places=6)
        self.assertAlmostEqual(float(edges[0][3]), 1.0, places=6)

    def test_a_pair_on_a_contig_that_is_not_listed_is_ignored(self):
        p, out = self.run_atoms([sedef_row("A", 0, 5000, "B", 20000, 25000)], ["A"])
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(self.read(out + ".nodes.tsv"), [])

    def test_cigar_mode_aborts_on_a_row_without_a_cigar(self):
        bad = sedef_row("A", 0, 5000, "A", 20000, 25000, cigar="")
        p, _ = self.run_atoms([bad], ["A"])
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("ABORT", p.stderr)

    def test_a_blocked_alignment_gives_a_partial_coverage_edge(self):
        # side A 0-6000 aligned to B 20000-26000 as 2000M2000D2000M ... : the gap removes the middle of A from the alignment
        row = sedef_row("A", 0, 6000, "A", 20000, 24000, cigar="2000M2000D2000M")
        p, out = self.run_atoms([row], ["A"])
        self.assertEqual(p.returncode, 0, p.stderr)
        edges = self.read(out + ".edges.tsv")
        self.assertEqual(len(edges), 1)
        cov = float(edges[0][3])
        self.assertLess(cov, 1.0)
        self.assertGreater(cov, 0.0)


# ---------------------------------------------------------------------------------------------------------------------------
# group B: the SEDEF table (header line, mitochondrial rows, the identical-tuple rule, the tau filter)
# ---------------------------------------------------------------------------------------------------------------------------
import t1_gorilla_pairs as T  # noqa: E402  (imported late so that group A fails and passes on its own)


class SedefReadingTests(unittest.TestCase):
    def test_the_header_line_and_the_mitochondrial_rows_are_dropped_and_counted(self):
        lines = [
            sedef_row("NC_011120.1", 100, 900, "NC_073224.2", 1000, 1800),
            sedef_row("NC_073224.2", 0, 5000, "NC_073227.2", 100000, 105000),
            sedef_row("NC_073224.2", 9000, 9800, "NC_011120.1", 100, 900),
            SEDEF_HEADER,
        ]
        rows, stats = T.read_sedef(lines)
        self.assertEqual(len(rows), 1)
        self.assertEqual(stats["data_rows"], 3)
        self.assertEqual(stats["header_lines"], 1)
        self.assertEqual(stats["mito_rows_dropped"], 2)
        self.assertEqual(rows[0].key, ("NC_073224.2", 0, 5000, "NC_073227.2", 100000, 105000))

    def test_columns_are_read_by_the_registered_one_based_numbers(self):
        # column 21 = fracMatch, column 12 = aln_len (1-based)
        r = sedef_row("A", 0, 5000, "B", 0, 5000, frac=0.9123, aln_len=4321)
        rows, _ = T.read_sedef([r])
        self.assertAlmostEqual(rows[0].frac, 0.9123, places=6)
        self.assertEqual(rows[0].aln_len, 4321)

    def test_a_row_that_does_not_have_34_columns_stops_the_run(self):
        with self.assertRaises(ValueError):
            T.read_sedef(["A\t0\t5000\tB\t0\t5000\tS"])

    def test_blank_lines_are_ignored(self):
        rows, stats = T.read_sedef(["", sedef_row("A", 0, 5000, "B", 0, 5000), ""])
        self.assertEqual((len(rows), stats["data_rows"]), (1, 1))


class DedupeTests(unittest.TestCase):
    def rows(self, specs):
        lines = [sedef_row("A", 0, 5000, "B", 10, 5010, frac=f, aln_len=n, extra=tag) for f, n, tag in specs]
        rows, _ = T.read_sedef(lines)
        return rows

    def test_an_identical_tuple_keeps_the_row_with_the_largest_fracmatch(self):
        rows = self.rows([(0.91, 100, ";x1"), (0.95, 100, ";x2"), (0.93, 100, ";x3")])
        kept, removed = T.dedupe(rows)
        self.assertEqual((len(kept), removed), (1, 2))
        self.assertIn(";x2", kept[0].line)

    def test_ties_go_to_the_larger_aln_len_then_to_the_first_in_file(self):
        rows = self.rows([(0.95, 100, ";a"), (0.95, 300, ";b"), (0.95, 300, ";c")])
        kept, _ = T.dedupe(rows)
        self.assertIn(";b", kept[0].line)

    def test_the_reverse_orientation_of_a_tuple_is_a_different_tuple(self):
        a = sedef_row("A", 0, 5000, "B", 10, 5010)
        b = sedef_row("B", 10, 5010, "A", 0, 5000)
        rows, _ = T.read_sedef([a, b])
        kept, removed = T.dedupe(rows)
        self.assertEqual((len(kept), removed), (2, 0))

    def test_survivors_keep_file_order(self):
        lines = [sedef_row("A", 0, 5000, "B", 0, 5000, frac=0.91), sedef_row("A", 9000, 9500, "B", 9000, 9500), sedef_row("A", 0, 5000, "B", 0, 5000, frac=0.97)]
        rows, _ = T.read_sedef(lines)
        kept, _ = T.dedupe(rows)
        self.assertEqual([r.key[1] for r in kept], [9000, 0])


class TauFilterTests(unittest.TestCase):
    def test_a_row_is_retained_iff_fracmatch_is_at_least_tau(self):
        lines = [sedef_row("A", i * 10000, i * 10000 + 5000, "B", i * 10000, i * 10000 + 5000, frac=f) for i, f in enumerate([0.8999, 0.9, 0.95, 0.98, 0.9801])]
        rows, _ = T.read_sedef(lines)
        self.assertEqual(len(T.retain(rows, 0.90)), 4)
        self.assertEqual(len(T.retain(rows, 0.95)), 3)
        self.assertEqual(len(T.retain(rows, 0.98)), 2)

    def test_dedupe_then_filter_equals_filter_then_dedupe(self):
        lines = [sedef_row("A", 0, 5000, "B", 0, 5000, frac=f) for f in (0.85, 0.96, 0.91)] + [sedef_row("A", 9000, 9500, "B", 9000, 9500, frac=0.88)]
        rows, _ = T.read_sedef(lines)
        for tau in (0.90, 0.95, 0.98):
            one = [r.key for r in T.retain(T.dedupe(rows)[0], tau)]
            two = [r.key for r in T.dedupe(T.retain(rows, tau))[0]]
            self.assertEqual(one, two)

    def test_the_filtered_table_is_written_with_the_original_lines(self):
        lines = [sedef_row("A", 0, 5000, "B", 0, 5000, frac=0.96), sedef_row("A", 9000, 9500, "B", 9000, 9500, frac=0.8)]
        rows, _ = T.read_sedef(lines)
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "o.bed")
            T.write_sedef(T.retain(rows, 0.9), path)
            with open(path) as fh:
                self.assertEqual(fh.read(), lines[0] + "\n")


# ---------------------------------------------------------------------------------------------------------------------------
# group C: atoms, atom classes (union-find over atom edges), the 100 bp touch rule
# ---------------------------------------------------------------------------------------------------------------------------
class ClassBuilderTests(unittest.TestCase):
    def test_an_edge_needs_two_sided_coverage_of_at_least_one_half_inclusive(self):
        edges = [(0, 1, 0.9, 0.5), (1, 2, 0.9, 0.4999), (3, 4, 0.9, 1.0)]
        class_of, used = T.build_classes(5, edges)
        self.assertEqual(used, 2)
        self.assertEqual(class_of[0], class_of[1])
        self.assertNotEqual(class_of[1], class_of[2])
        self.assertEqual(class_of[3], class_of[4])

    def test_the_class_id_is_the_smallest_atom_index_of_the_class(self):
        edges = [(4, 2, 0.9, 1.0), (2, 7, 0.9, 1.0), (9, 8, 0.9, 1.0)]
        class_of, _ = T.build_classes(10, edges)
        self.assertEqual(class_of[7], 2)
        self.assertEqual(class_of[4], 2)
        self.assertEqual(class_of[9], 8)
        self.assertEqual(class_of[0], 0)

    def test_classes_are_transitive_over_chains_and_independent_of_edge_order(self):
        edges = [(0, 1, 0.9, 1.0), (1, 2, 0.9, 1.0), (2, 3, 0.9, 1.0)]
        a, _ = T.build_classes(5, edges)
        b, _ = T.build_classes(5, list(reversed(edges)))
        self.assertEqual(a, b)
        self.assertEqual(len(set(a[:4])), 1)
        self.assertEqual(a[4], 4)

    def test_class_statistics(self):
        class_of, _ = T.build_classes(6, [(0, 1, 0.9, 1.0), (1, 2, 0.9, 1.0), (3, 4, 0.9, 1.0)])
        st = T.class_stats(class_of)
        self.assertEqual(st["classes"], 3)
        self.assertEqual(st["largest"], 3)
        self.assertEqual(st["singletons"], 1)


class TouchRuleTests(unittest.TestCase):
    def setUp(self):
        # contig A: atoms 0 [100,200) class 0, 1 [300,400) class 0 (same class), 2 [1000,1500) class 2; contig B: atom 3 [0,1000) class 3
        atoms = [("A", 100, 200), ("A", 300, 400), ("A", 1000, 1500), ("B", 0, 1000)]
        self.atoms = atoms
        self.index = T.AtomIndex(atoms)
        self.class_of = [0, 0, 2, 3]

    def test_the_index_returns_overlapping_atoms_with_their_overlap_in_bp(self):
        self.assertEqual(self.index.overlaps("A", 150, 350), [(0, 50), (1, 50)])
        self.assertEqual(self.index.overlaps("A", 200, 300), [])
        self.assertEqual(self.index.overlaps("A", 0, 100), [])
        self.assertEqual(self.index.overlaps("C", 0, 10 ** 6), [])
        self.assertEqual(self.index.overlaps("B", 500, 600), [(3, 100)])

    def test_bp_are_summed_per_class_over_the_gene_exon_blocks(self):
        t = T.touch_classes([(150, 200), (300, 350), (1200, 1210)], "A", self.index, self.class_of)
        self.assertEqual(t, {0: 100, 2: 10})

    def test_a_gene_touches_a_class_iff_it_has_at_least_100_exonic_bp_in_its_atoms(self):
        self.assertEqual(T.gene_U({0: 100, 2: 99}), frozenset({0}))
        self.assertEqual(T.gene_U({0: 99}), frozenset())
        self.assertEqual(T.gene_U({}), frozenset())

    def test_60_and_40_bp_in_two_atoms_of_one_class_touch_the_class_but_60_and_60_in_two_classes_touch_none(self):
        one = T.gene_U(T.touch_classes([(140, 200), (300, 340)], "A", self.index, self.class_of))
        self.assertEqual(one, frozenset({0}))
        # 60 bp in class 0 (atom 0) and 60 bp in class 2 (atom 2)
        two = T.gene_U(T.touch_classes([(140, 200), (1000, 1060)], "A", self.index, self.class_of))
        self.assertEqual(two, frozenset())


class AtomsRunnerTests(unittest.TestCase):
    def table(self, rows):
        d = tempfile.mkdtemp(prefix="r2_prepare_")
        return d

    def test_prepare_tau_builds_atoms_and_classes_for_a_small_table(self):
        lines = [
            sedef_row("NC_073224.2", 0, 5000, "NC_073227.2", 100000, 105000, frac=0.96),
            sedef_row("NC_073224.2", 50000, 56000, "NC_073227.2", 300000, 306000, frac=0.92),
            sedef_row("NC_073224.2", 70000, 75000, "NC_073227.2", 500000, 505000, frac=0.80),   # below every tau
            SEDEF_HEADER,
        ]
        rows, _ = T.read_sedef(lines)
        rows, _ = T.dedupe(rows)
        d = tempfile.mkdtemp(prefix="r2_prepare_")
        info = T.prepare_tau(rows, 0.90, d)
        self.assertEqual(info["rows_retained"], 2)
        self.assertEqual(info["atoms"], 4)
        self.assertEqual(info["edges_total"], 2)
        self.assertEqual(info["edges_used"], 2)
        self.assertEqual(info["classes"], 2)
        info95 = T.prepare_tau(rows, 0.95, d)
        self.assertEqual((info95["rows_retained"], info95["atoms"], info95["classes"]), (1, 2, 1))
        self.assertTrue(os.path.exists(os.path.join(d, "sd_ge090.bed")))
        atoms = T.read_nodes(os.path.join(d, "out090"))
        self.assertEqual(atoms[0], ("NC_073224.2", 0, 5000))

    def test_the_runner_refuses_a_script_whose_hash_is_not_the_registered_one(self):
        d = tempfile.mkdtemp(prefix="r2_prepare_")
        with self.assertRaises(RuntimeError):
            T.check_atoms_script(expected_sha256="0" * 64)
        T.check_atoms_script()


# ---------------------------------------------------------------------------------------------------------------------------
# group D: genes, the Compara layer, pair eligibility, distance classes, the verdict set, depth-matched pairs
# ---------------------------------------------------------------------------------------------------------------------------
GENES_TXT = "\n".join([
    "contig\tname\ttype\tbiotype\tstart1\tend\texons",
    "C1\tg0\tgene\tprotein_coding\t1001\t3000\t1000-1500,1400-1600,2500-3000",
    "C1\tg1\tgene\tprotein_coding\t5001\t7000\t5000-5500",
    "C1\tg2\tgene\tprotein_coding\t500001\t502000\t500000-500400",
    "C1\tg3\tgene\tprotein_coding\t2000001\t2001000\t2000000-2000300",
    "C2\tg4\tpseudogene\tpseudogene\t101\t900\t100-900",
    "C2\tg5\tgene\tlncRNA\t1001\t1300\t",
    "NC_073241.2\tg6\tgene\tprotein_coding\t1001\t2000\t1000-1800",
    "C1\tg7\tgene\tprotein_coding\t1501\t2800\t1500-2700",
]) + "\n"

TRUTH_HEAD = "family_id\tfamily_status\tmember_status\tgorilla_gff_id\tgorilla_contig\tstart1\tend1"
TRUTH_TXT = "\n".join([
    TRUTH_HEAD,
    "F1\twhole\t1to1\tgene-g0\tC1\t1001\t3000",
    "F1\twhole\t1to1\tgene-g1\tC1\t5001\t7000",
    "F2\tpartial\t1to1\tgene-g2\tC1\t500001\t502000",
    "F2\tpartial\tabsent\t\t\t\t",
    "F3\tnone\t1to1\tgene-g4\tC2\t101\t900",
    "F4\twhole\t1to1\tgene-g3\tC1\t2000001\t2001000",
    "F4\twhole\t1to1\tgene-g6\tNC_073241.2\t1001\t2000",
    "F5\twhole\t1to1\tgene-g7\tC1\t1501\t2800",
]) + "\n"


class GeneAndLayerTests(unittest.TestCase):
    def setUp(self):
        self.genes, self.index, self.dups = T.read_genes(GENES_TXT.splitlines())

    def test_exon_blocks_are_merged_and_the_join_key_is_contig_start1_end(self):
        g0 = self.genes[self.index[("C1", 1001, 3000)]]
        self.assertEqual(g0.blocks, ((1000, 1600), (2500, 3000)))
        self.assertEqual(g0.name, "g0")
        self.assertEqual(self.genes[self.index[("C2", 1001, 1300)]].blocks, ())
        self.assertEqual(self.dups, [])

    def test_a_duplicated_join_key_is_reported(self):
        lines = GENES_TXT.splitlines() + ["C1\tdup\tgene\tprotein_coding\t1001\t3000\t1000-1100"]
        _g, _i, dups = T.read_genes(lines)
        self.assertEqual(dups, [("C1", 1001, 3000)])

    def test_the_layer_keeps_1to1_rows_only_and_counts_all_families(self):
        layer = T.read_layer(TRUTH_TXT.splitlines(), self.index)
        self.assertEqual(sorted(layer.families), ["F1", "F2", "F3", "F4", "F5"])
        self.assertEqual(layer.n_families_all, 5)
        self.assertEqual(len(layer.families["F1"]), 2)
        self.assertEqual(layer.status["F2"], "partial")
        g = layer.gff_to_gene["gene-g0"]
        self.assertEqual(self.genes[g].name, "g0")
        self.assertEqual(layer.gene_family[g], "F1")

    def test_a_1to1_row_whose_key_is_not_in_the_gene_table_stops_the_run(self):
        bad = TRUTH_TXT + "F9\twhole\t1to1\tgene-zz\tC1\t9\t99\n"
        with self.assertRaises(ValueError):
            T.read_layer(bad.splitlines(), self.index)


class PairUniverseTests(unittest.TestCase):
    def setUp(self):
        self.genes, self.index, _ = T.read_genes(GENES_TXT.splitlines())
        self.layer = T.read_layer(TRUTH_TXT.splitlines(), self.index)
        self.byname = {g.name: i for i, g in enumerate(self.genes)}

    def test_exon_overlap_counts_shared_bp_of_the_exon_unions(self):
        self.assertEqual(T.exon_overlap_bp(((0, 100), (200, 300)), ((50, 250),)), 100)
        self.assertEqual(T.exon_overlap_bp(((0, 100),), ((100, 200),)), 0)
        g0, g7 = self.genes[self.byname["g0"]], self.genes[self.byname["g7"]]
        self.assertEqual(T.exon_overlap_bp(g0.blocks, g7.blocks), 100 + 200)   # (1500,1600) and (2500,2700)

    def test_genes_on_different_contigs_never_share_exonic_bp_even_at_identical_coordinates(self):
        genes, _i, _d = T.read_genes(["contig\tname\ttype\tbiotype\tstart1\tend\texons",
                                      "C1\tx\tgene\tprotein_coding\t1001\t3000\t1000-2000",
                                      "C2\ty\tgene\tprotein_coding\t1001\t3000\t1000-2000",
                                      "C1\tz\tgene\tprotein_coding\t1501\t3000\t1500-2000"])
        U = {0: frozenset({1}), 1: frozenset({1}), 2: frozenset({1})}
        self.assertEqual(T.shared_exonic_bp(genes[0], genes[1]), 0)
        self.assertTrue(T.pair_eligible(0, 1, genes, U))
        self.assertEqual(T.shared_exonic_bp(genes[0], genes[2]), 500)
        self.assertFalse(T.pair_eligible(0, 2, genes, U))

    def test_distance_classes_and_their_boundaries(self):
        def G(c, s, e):
            return T.Gene(c, "x", "gene", "pc", s + 1, e, ((s, e),))
        a = G("C1", 0, 1000)
        self.assertEqual(T.distance_class(a, G("C1", 500, 1500)), "same<100kb")           # overlapping spans: gap 0
        self.assertEqual(T.distance_class(a, G("C1", 100999, 102000)), "same<100kb")       # gap 99,999
        self.assertEqual(T.distance_class(a, G("C1", 101000, 102000)), "same100kb-1Mb")    # gap 100,000
        self.assertEqual(T.distance_class(a, G("C1", 1000999, 1002000)), "same100kb-1Mb")  # gap 999,999
        self.assertEqual(T.distance_class(a, G("C1", 1001000, 1002000)), "same>1Mb")       # gap 1,000,000
        self.assertEqual(T.distance_class(a, G("C2", 0, 1000)), "cross-contig")

    def test_the_verdict_set_excludes_every_family_with_a_member_on_the_four_development_contigs(self):
        v = T.verdict_families(self.layer, self.genes)
        self.assertEqual(sorted(v), ["F1", "F2", "F3", "F5"])      # F4 has a member on NC_073241.2

    def test_a_pair_is_eligible_iff_both_genes_have_a_class_differ_and_share_no_exon_bp(self):
        U = {self.byname["g0"]: frozenset({1}), self.byname["g1"]: frozenset({1}), self.byname["g7"]: frozenset({1}), self.byname["g2"]: frozenset()}
        ok = lambda a, b: T.pair_eligible(self.byname[a], self.byname[b], self.genes, U)   # noqa: E731
        self.assertTrue(ok("g0", "g1"))
        self.assertFalse(ok("g0", "g2"))          # g2 has no class
        self.assertFalse(ok("g0", "g0"))          # same gene
        self.assertFalse(ok("g0", "g7"))          # exon unions share 300 bp
        self.assertFalse(ok("g0", "g3"))          # g3 not in U at all

    def test_same_family_and_different_family_pools_are_counted_by_distance_class(self):
        names = ["g0", "g1", "g2", "g3", "g4", "g7"]
        U = {self.byname[n]: frozenset({0}) for n in names}
        verdict = T.verdict_families(self.layer, self.genes)
        same = T.same_family_pairs(self.layer, self.genes, U, verdict)
        self.assertEqual([(f, self.genes[a].name, self.genes[b].name, c) for f, a, b, c in same], [("F1", "g0", "g1", "same<100kb")])
        diff = T.control_pairs(self.layer, self.genes, U, verdict)
        counts = collections.Counter(c for *_x, c in diff)
        # verdict pool genes: g0, g1 (F1), g2 (F2), g4 (F3), g7 (F5); g0-g7 share exon bp and are dropped; g3 (F4) is outside the verdict set
        self.assertEqual(sum(counts.values()), 8)       # 10 pairs of 5 genes, minus the same-family pair g0-g1, minus g0-g7
        self.assertEqual(counts["cross-contig"], 4)
        self.assertEqual(counts["same<100kb"], 1)       # g1-g7
        self.assertEqual(counts["same100kb-1Mb"], 3)    # g0-g2, g1-g2, g2-g7

    def test_partial_family_pairs_are_counted_apart(self):
        # make F2 have two eligible members: add a second projected gene to F2 through a custom layer
        txt = TRUTH_TXT + "F2\tpartial\t1to1\tgene-g5\tC2\t1001\t1300\n"
        layer = T.read_layer(txt.splitlines(), self.index)
        U = {self.byname["g2"]: frozenset({0}), self.byname["g5"]: frozenset({0})}
        same = T.same_family_pairs(layer, self.genes, U, T.verdict_families(layer, self.genes))
        self.assertEqual([(f, layer.status[f]) for f, *_x in same], [("F2", "partial")])


PAIRS_HEAD = "family_id\ta\tb\ta_gene\tb_gene\ttarget\ttlen\tqlen\tn_records\tcov_short\tcov_long\tidentity\tidentity_cs\tmatches\tblock\tpass_C4"


def pairs_row(fam, ag, bg, ident, cov_short=0.8):
    return "\t".join([fam, "a", "b", ag, bg, ag, "1000", "1000", "1", str(cov_short), "0.5", ident, ident, "10", "10", "1"])


class DepthMatchedTests(unittest.TestCase):
    def setUp(self):
        self.genes, self.index, _ = T.read_genes(GENES_TXT.splitlines())
        self.layer = T.read_layer(TRUTH_TXT.splitlines(), self.index)
        self.byname = {g.name: i for i, g in enumerate(self.genes)}
        self.U = {self.byname[n]: frozenset({0}) for n in ("g0", "g1", "g2", "g3", "g4", "g6", "g7")}
        self.verdict = T.verdict_families(self.layer, self.genes)

    def rows(self, lines):
        return T.read_pairs([PAIRS_HEAD] + lines)

    def test_identity_na_is_excluded_and_identity_at_least_tau_is_depth_matched(self):
        rows = self.rows([pairs_row("F1", "gene-g0", "gene-g1", "0.9"), pairs_row("F1", "gene-g0", "gene-g1", "NA"), pairs_row("F1", "gene-g0", "gene-g1", "0.89")])
        self.assertEqual(rows[1]["identity"], None)
        dm = T.depth_matched(rows, self.layer, self.genes, self.U, 0.90, self.verdict)
        self.assertEqual(len(dm), 1)
        self.assertEqual(dm[0]["family_id"], "F1")
        self.assertEqual(len(T.depth_matched(rows, self.layer, self.genes, self.U, 0.95, self.verdict)), 0)

    def test_a_pair_with_shared_exon_bp_or_a_gene_without_a_class_is_not_eligible(self):
        txt = TRUTH_TXT + "F1\twhole\t1to1\tgene-g7\tC1\t1501\t2800\n"
        txt = txt.replace("F5\twhole\t1to1\tgene-g7\tC1\t1501\t2800\n", "")
        layer = T.read_layer(txt.splitlines(), self.index)
        verdict = T.verdict_families(layer, self.genes)
        rows = self.rows([pairs_row("F1", "gene-g0", "gene-g7", "0.95"), pairs_row("F1", "gene-g0", "gene-g1", "0.95")])
        U = dict(self.U)
        dm = T.depth_matched(rows, layer, self.genes, U, 0.90, verdict)
        self.assertEqual([(r["a_gene"], r["b_gene"]) for r in dm], [("gene-g0", "gene-g1")])     # g7 and g0 share exonic bp
        U[self.byname["g1"]] = frozenset()
        self.assertEqual(T.depth_matched(rows, layer, self.genes, U, 0.90, verdict), [])

    def test_a_row_whose_family_differs_from_the_family_of_its_genes_stops_the_run(self):
        rows = self.rows([pairs_row("F5", "gene-g0", "gene-g1", "0.95")])
        with self.assertRaises(ValueError):
            T.depth_matched(rows, self.layer, self.genes, self.U, 0.90, self.verdict)

    def test_the_verdict_set_filter_and_the_whole_family_count(self):
        rows = self.rows([pairs_row("F4", "gene-g3", "gene-g6", "0.97"), pairs_row("F1", "gene-g0", "gene-g1", "0.97")])
        in_verdict = T.depth_matched(rows, self.layer, self.genes, self.U, 0.90, self.verdict)
        everywhere = T.depth_matched(rows, self.layer, self.genes, self.U, 0.90, None)
        self.assertEqual(([r["family_id"] for r in in_verdict], [r["family_id"] for r in everywhere]), (["F1"], ["F1", "F4"]))
        self.assertEqual(len(T.eligible_rows(rows, self.layer, self.genes, self.U, self.verdict)), 1)

    def test_a_pair_row_that_names_a_gene_outside_the_layer_stops_the_run(self):
        rows = self.rows([pairs_row("F1", "gene-g0", "gene-nope", "0.95")])
        with self.assertRaises(ValueError):
            T.depth_matched(rows, self.layer, self.genes, self.U, 0.90, self.verdict)

    def test_the_aligned_fraction_is_carried_beside_the_identity(self):
        rows = self.rows([pairs_row("F1", "gene-g0", "gene-g1", "0.95", cov_short=0.31)])
        dm = T.depth_matched(rows, self.layer, self.genes, self.U, 0.90, self.verdict)
        self.assertAlmostEqual(dm[0]["cov_short"], 0.31)


# ---------------------------------------------------------------------------------------------------------------------------
# group E: the report statistics (exercised on made-up tables only)
# ---------------------------------------------------------------------------------------------------------------------------
class WilsonTests(unittest.TestCase):
    def test_known_wilson_intervals(self):
        lo, hi = T.wilson(10, 20)
        self.assertAlmostEqual(lo, 0.2993, places=3)
        self.assertAlmostEqual(hi, 0.7007, places=3)
        lo, hi = T.wilson(0, 20)
        self.assertAlmostEqual(lo, 0.0, places=6)
        self.assertAlmostEqual(hi, 0.1611, places=3)
        lo, hi = T.wilson(38, 38)
        self.assertAlmostEqual(lo, 0.9082, places=3)
        self.assertAlmostEqual(hi, 1.0, places=6)

    def test_an_empty_table_has_no_interval(self):
        self.assertEqual(T.wilson(0, 0), (None, None))


class RateTests(unittest.TestCase):
    def setUp(self):
        # family A: three pairs, all sharing a class; family B: one pair sharing none
        self.U = {1: frozenset({7}), 2: frozenset({7}), 3: frozenset({7, 8}), 4: frozenset({8}), 5: frozenset({9}), 6: frozenset({10})}
        self.pairs = [("A", 1, 2), ("A", 2, 3), ("A", 3, 4), ("B", 5, 6)]

    def test_pair_weighted_and_family_weighted_rates(self):
        r = T.s_rate(self.pairs, self.U)
        self.assertEqual((r["n_pairs"], r["n_families"], r["shared"]), (4, 2, 3))
        self.assertAlmostEqual(r["pair_weighted"], 0.75)
        self.assertAlmostEqual(r["family_weighted"], 0.5)
        self.assertAlmostEqual(r["a_viol"], 0.25)
        self.assertEqual(len(r["wilson"]), 2)

    def test_b_viol_counts_eligible_pairs_of_different_families_that_share_a_class(self):
        control = [("A", "B", 1, 5, "cross-contig"), ("A", "B", 3, 4, "same<100kb"), ("A", "B", 2, 6, "same>1Mb")]
        self.assertEqual(T.b_viol(control, self.U), (1, 3))

    def test_the_underpowered_flag(self):
        self.assertTrue(T.underpowered(19, 10))
        self.assertTrue(T.underpowered(20, 7))
        self.assertFalse(T.underpowered(20, 8))


class BootstrapTests(unittest.TestCase):
    def test_percentile_interval_uses_floor_and_ceil_positions(self):
        vals = list(range(2000))
        self.assertEqual(T.percentile_interval(vals), (50, 1949))
        self.assertEqual(T.percentile_interval(list(reversed(vals))), (50, 1949))

    def test_family_cluster_bootstrap_is_reproducible_and_bounded(self):
        by_family = {"A": [1, 1, 1], "B": [0], "C": [1, 0], "D": [1]}
        a = T.bootstrap_rates(by_family)
        b = T.bootstrap_rates(by_family)
        self.assertEqual(a, b)
        self.assertEqual(a["B"], 2000)
        self.assertEqual(a["seed"], 20260930)
        lo, hi = a["pair_weighted"]
        self.assertTrue(0.0 <= lo <= hi <= 1.0)
        v1 = T.bootstrap_values(by_family)
        self.assertEqual(v1, T.bootstrap_values(by_family))
        self.assertNotEqual(v1, T.bootstrap_values(by_family, seed=1))

    def test_all_shared_gives_a_degenerate_interval_at_one(self):
        r = T.bootstrap_rates({"A": [1, 1], "B": [1], "C": [1, 1, 1]})
        self.assertEqual(r["pair_weighted"], (1.0, 1.0))
        self.assertEqual(r["family_weighted"], (1.0, 1.0))

    def test_the_pair_weighted_and_the_family_weighted_resamples_use_the_same_family_draws(self):
        r = T.bootstrap_rates({"A": [1, 1, 1, 1], "B": [0]}, B=200, seed=5)
        self.assertEqual(r["B"], 200)
        self.assertTrue(r["pair_weighted"][0] <= r["pair_weighted"][1])


class ReportTableTests(unittest.TestCase):
    def test_a_made_up_table_is_assembled_with_every_registered_quantity(self):
        U = {1: frozenset({7}), 2: frozenset({7}), 3: frozenset({8}), 4: frozenset({8}), 5: frozenset({9})}
        dm = [{"family_id": "A", "g": 1, "h": 2, "identity": 0.95, "cov_short": 0.9}, {"family_id": "B", "g": 3, "h": 5, "identity": 0.92, "cov_short": 0.3}]
        control = [("A", "B", 1, 3, "cross-contig"), ("A", "B", 2, 4, "cross-contig")]
        t = T.report_table(dm, control, U, tau=0.90)
        self.assertEqual((t["n_pairs"], t["n_families"], t["shared"]), (2, 2, 1))
        self.assertTrue(t["underpowered"])
        self.assertEqual(t["b_viol"], [0, 2])
        self.assertEqual(t["short_aligned_pairs"], 1)        # cov_short < 0.5
        self.assertEqual(t["bootstrap"]["B"], 2000)
        self.assertEqual(t["tau"], 0.90)


# ---------------------------------------------------------------------------------------------------------------------------
# group F: the positive control, Gate 0 fingerprints, the counts command on a made-up world, determinism, the report guard
# ---------------------------------------------------------------------------------------------------------------------------
LIFTOFF_HEAD = "species\tshard\tsource_id\tname\trtype\tbiotype\tstratum\tshort\tsrc_contig\tsrc_start0\tsrc_end\tsrc_strand\tcls\ttag\tcontig\tstart0\tend\tstrand\tcoverage\tsequence_id\texons\texonic_len\tnote"


def liftoff_row(source, name, cls, contig, start0, end, seqid, note):
    return "\t".join(["gorilla", "x", source, name, "gene", "pc", "other", "0", contig, "0", "1", "+", cls, "0", contig, str(start0), str(end), "+", "1.0", seqid, "0-1", "1", note])


class PositiveControlTests(unittest.TestCase):
    def setUp(self):
        names = ["LRPAP1", "L1", "L2", "L3", "L4", "L5", "L6", "L7"]
        lines = ["contig\tname\ttype\tbiotype\tstart1\tend\texons"]
        for i, n in enumerate(names):
            lines.append(f"CT{i}\t{n}\tgene\tprotein_coding\t1001\t3000\t1000-2200")
        self.genes, self.index, _ = T.read_genes(lines)
        rows = [LIFTOFF_HEAD, liftoff_row("gene-LRPAP1", "LRPAP1", "in_place", "CT0", 1000, 3000, "1.0", "-")]
        for i, n in enumerate(names[1:], start=1):
            rows.append(liftoff_row("gene-LRPAP1", "LRPAP1", "dropped_M1", f"CT{i}", 1000, 3000, "0.97", f"overlaps gene-{n} (CT{i})"))
        rows.append(liftoff_row("gene-OTHER", "OTHER", "in_place", "CT0", 1, 2, "1.0", "-"))
        self.lines = rows

    def test_the_eight_copies_are_the_gene_lrpap1_rows_of_the_self_lift(self):
        c = T.lrpap_copies(self.lines, self.genes)
        self.assertEqual(len(c), 8)
        self.assertEqual([x["name"] for x in c][:3], ["LRPAP1", "L1", "L2"])
        self.assertEqual(c[0]["cls"], "in_place")

    def test_a_table_without_eight_copies_stops_the_run(self):
        with self.assertRaises(ValueError):
            T.lrpap_copies(self.lines[:-2], self.genes)

    def test_the_control_passes_when_every_copy_has_a_class_and_the_copy_graph_is_connected(self):
        c = T.lrpap_copies(self.lines, self.genes)
        U = {x["gene"]: frozenset({1}) for x in c}
        U[c[3]["gene"]] = frozenset({2, 1})
        r = T.positive_control(c, U)
        self.assertTrue(r["passed"])
        self.assertTrue(r["connected"])
        self.assertEqual(r["without_class"], [])

    def test_the_control_fails_for_a_copy_without_a_class_and_for_a_disconnected_graph(self):
        c = T.lrpap_copies(self.lines, self.genes)
        U = {x["gene"]: frozenset({1}) for x in c}
        U[c[2]["gene"]] = frozenset()
        r = T.positive_control(c, U)
        self.assertFalse(r["passed"])
        self.assertEqual(r["without_class"], ["L2"])
        U = {x["gene"]: frozenset({1}) for x in c}
        U[c[5]["gene"]] = frozenset({9})
        r = T.positive_control(c, U)
        self.assertFalse(r["passed"])
        self.assertFalse(r["connected"])
        self.assertEqual(len(r["components"]), 2)

    def test_the_identity_range_of_the_non_reference_copies_is_reported(self):
        c = T.lrpap_copies(self.lines, self.genes)
        U = {x["gene"]: frozenset({1}) for x in c}
        self.assertEqual(T.positive_control(c, U)["identity_range"], (0.97, 0.97))


class FingerprintTests(unittest.TestCase):
    def test_size_lines_and_sha256(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.txt")
            with open(p, "w") as fh:
                fh.write("abc\ndef\n")
            f = T.fingerprint(p)
            self.assertEqual((f["size"], f["lines"]), (8, 2))
            self.assertEqual(f["sha256"], hashlib.sha256(b"abc\ndef\n").hexdigest())

    def test_gate0_marks_a_prefix_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.txt")
            with open(p, "w") as fh:
                fh.write("abc\n")
            good = hashlib.sha256(b"abc\n").hexdigest()[:16]
            rows = T.gate0_rows({"sedef": p}, {"sedef": good})
            self.assertTrue(rows[0]["ok"])
            rows = T.gate0_rows({"sedef": p}, {"sedef": "0" * 16})
            self.assertFalse(rows[0]["ok"])


def make_world(d):
    """A made-up gorilla: two contigs, two SD pairs (two atom classes at tau 0.90), seven genes, four families."""
    sedef = [
        sedef_row("NC_073224.2", 0, 5000, "NC_073227.2", 100000, 105000, frac=0.96),
        sedef_row("NC_073224.2", 0, 5000, "NC_073227.2", 100000, 105000, frac=0.93),     # identical tuple, lower fracMatch: removed
        sedef_row("NC_073224.2", 50000, 56000, "NC_073227.2", 300000, 306000, frac=0.92),
        sedef_row("NC_073224.2", 70000, 75000, "NC_073227.2", 500000, 505000, frac=0.80),
        sedef_row("NC_011120.1", 100, 900, "NC_073224.2", 9000, 9800, frac=0.97),
        SEDEF_HEADER,
    ]
    genes = ["contig\tname\ttype\tbiotype\tstart1\tend\texons",
             "NC_073224.2\tgA1\tgene\tprotein_coding\t1001\t3000\t1000-2200",
             "NC_073227.2\tgA2\tgene\tprotein_coding\t100501\t103000\t101000-102000",
             "NC_073224.2\tgC1\tgene\tprotein_coding\t3001\t4500\t3500-4200",
             "NC_073224.2\tgB1\tgene\tprotein_coding\t51001\t53000\t51000-52500",
             "NC_073227.2\tgB2\tgene\tprotein_coding\t301001\t303000\t301000-302000",
             "NC_073227.2\tgD1\tgene\tprotein_coding\t304001\t305000\t304000-304500",
             "NC_073224.2\tgN\tgene\tprotein_coding\t80001\t82000\t80000-81000"]
    truth = [TRUTH_HEAD,
             "F1\twhole\t1to1\tgene-gA1\tNC_073224.2\t1001\t3000", "F1\twhole\t1to1\tgene-gA2\tNC_073227.2\t100501\t103000",
             "F2\twhole\t1to1\tgene-gB1\tNC_073224.2\t51001\t53000", "F2\twhole\t1to1\tgene-gB2\tNC_073227.2\t301001\t303000",
             "F3\tpartial\t1to1\tgene-gC1\tNC_073224.2\t3001\t4500", "F3\tpartial\t1to1\tgene-gD1\tNC_073227.2\t304001\t305000",
             "F4\twhole\t1to1\tgene-gN\tNC_073224.2\t80001\t82000", "F5\tnone\tabsent\t\t\t\t"]
    pairs = [PAIRS_HEAD, pairs_row("F1", "gene-gA1", "gene-gA2", "0.97"), pairs_row("F2", "gene-gB1", "gene-gB2", "0.92")]
    paths = {}
    for name, lines in (("sedef", sedef), ("genes", genes), ("truth", truth), ("pairs", pairs)):
        paths[name] = os.path.join(d, name + ".tsv")
        with open(paths[name], "w") as fh:
            fh.write("\n".join(lines) + "\n")
    return paths


def add_control(d, paths, scatter_one=False):
    """Append the eight LRPAP1 copies (genes + a self-lift table) to a made-up world; all sit in atom 0 unless scatter_one puts the last outside every atom."""
    names = ["LRPAP1", "L1", "L2", "L3", "L4", "L5", "L6", "L7"]
    with open(paths["genes"], "a") as fh:
        for i, n in enumerate(names):
            if scatter_one and i == 7:
                fh.write(f"NC_073224.2\t{n}\tgene\tprotein_coding\t{90001 + 300 * i}\t{90400 + 300 * i}\t{90000 + 300 * i}-{90300 + 300 * i}\n")   # outside every atom
            else:
                fh.write(f"NC_073224.2\t{n}\tgene\tprotein_coding\t{201 + 300 * i}\t{500 + 300 * i}\t{200 + 300 * i}-{500 + 300 * i}\n")     # inside atom 0
    rows = [LIFTOFF_HEAD]
    for i, n in enumerate(names):
        start0 = 200 + 300 * i
        if i == 0:
            rows.append(liftoff_row("gene-LRPAP1", "LRPAP1", "in_place", "NC_073224.2", start0, start0 + 300, "1.0", "-"))
        else:
            rows.append(liftoff_row("gene-LRPAP1", "LRPAP1", "dropped_M1", "NC_073224.2", start0, start0 + 300, "0.97", f"overlaps gene-{n} (NC_073224.2)"))
    paths["liftoff"] = os.path.join(d, "liftoff.tsv")
    with open(paths["liftoff"], "w") as fh:
        fh.write("\n".join(rows) + "\n")
    return paths


def world_with_control(scatter_one=False):
    d = tempfile.mkdtemp(prefix="r2_ctl_")
    return d, add_control(d, make_world(d), scatter_one)


def write_registry(path, paths):
    """A registry (name, sha256 prefix) of a made-up world's inputs: what the registered constants are for the real inputs."""
    with open(path, "w") as fh:
        for name in ("sedef", "genes", "truth", "pairs", "liftoff"):
            fh.write(f"{name}\t{T.fingerprint(paths[name])['sha256'][:16]}\n")
        fh.write(f"dna_sd_atoms.py\t{T.ATOMS_SHA256[:16]}\n")


def release_files(d, paths):
    """(registry path, release-file path) for a made-up world."""
    reg, rel = os.path.join(d, "registry.tsv"), os.path.join(d, "RELEASE")
    write_registry(reg, paths)
    with open(rel, "w") as fh:
        fh.write(T.RELEASE_TOKEN + "\n")
    return reg, rel


def report_argv(paths, reg, rel, out, *extra):
    return ["report", "--sedef", paths["sedef"], "--genes", paths["genes"], "--truth", paths["truth"], "--pairs", paths["pairs"],
            "--liftoff", paths["liftoff"], "--registry", reg, "--outdir", out, "--release-file", rel] + list(extra)


def run_main(argv):
    """T.main(argv) with stdout and stderr captured: (exit code, text)."""
    import contextlib
    import io
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        rc = T.main(argv)
    return rc, buf.getvalue()


class CountsCommandTests(unittest.TestCase):
    def run_counts(self, extra_env=None):
        d = tempfile.mkdtemp(prefix="r2_world_")
        paths = make_world(d)
        out = os.path.join(d, "out")
        cmd = [sys.executable, "-B", os.path.join(HERE, "t1_gorilla_pairs.py"), "counts", "--sedef", paths["sedef"], "--genes", paths["genes"],
               "--truth", paths["truth"], "--pairs", paths["pairs"], "--outdir", out]
        env = dict(os.environ, **(extra_env or {}))
        p = subprocess.run(cmd, capture_output=True, text=True, env=env)
        return p, out

    def test_the_counts_of_the_made_up_world(self):
        p, out = self.run_counts()
        self.assertEqual(p.returncode, 0, p.stderr)
        c = load_json(os.path.join(out, "counts.json"))
        s = c["sedef"]
        self.assertEqual((s["data_rows"], s["header_lines"], s["mito_rows_dropped"], s["duplicates_removed"]), (5, 1, 1, 1))
        t = c["tau"]["0.90"]
        self.assertEqual((t["rows_retained"], t["atoms"], t["classes"]), (2, 4, 2))
        self.assertEqual((t["labelled_genes"], t["eligible_labelled_genes"], t["path_genes"]), (7, 6, 0))
        self.assertEqual(t["same_family_pairs"]["verdict"], {"cross-contig": 3})
        self.assertEqual(t["same_family_pairs"]["partial_verdict"], 1)
        self.assertEqual(t["control_pairs"]["verdict"], {"cross-contig": 6, "same<100kb": 4, "same100kb-1Mb": 2})
        self.assertEqual(t["depth_matched"]["verdict"], {"pairs": 2, "families": 2})
        self.assertEqual(t["pairs_tsv_eligible_verdict"], 2)
        self.assertTrue(t["underpowered_verdict"])
        self.assertEqual(c["tau"]["0.95"]["depth_matched"]["verdict"], {"pairs": 1, "families": 1})
        self.assertEqual(c["tau"]["0.98"]["depth_matched"]["verdict"], {"pairs": 0, "families": 0})
        self.assertEqual(c["layer"]["genes_1to1"], 7)
        self.assertEqual(c["layer"]["families_all"], 5)
        self.assertEqual(c["layer"]["families_with_two_or_more_genes"], 3)
        self.assertEqual(c["pairs_tsv"]["rows"], 2)
        self.assertIsNone(c["positive_control"])

    def test_denominators_are_printed_and_no_s_rate_is_computed(self):
        p, out = self.run_counts()
        self.assertEqual(p.returncode, 0, p.stderr)
        for word in ("S-rate", "A_viol", "B_viol", "pair_weighted"):
            self.assertNotIn(word, p.stdout)
        self.assertNotIn("shared", read_text(os.path.join(out, "counts.json")))
        self.assertIn("denominators", p.stdout.lower())

    def test_the_counts_are_identical_under_two_hash_seeds(self):
        a, oa = self.run_counts({"PYTHONHASHSEED": "0"})
        b, ob = self.run_counts({"PYTHONHASHSEED": "1"})
        self.assertEqual(a.returncode, 0, a.stderr)
        self.assertEqual(b.returncode, 0, b.stderr)
        ja = load_json(os.path.join(oa, "counts.json"))
        jb = load_json(os.path.join(ob, "counts.json"))
        for j in (ja, jb):
            j.pop("inputs", None)
        self.assertEqual(ja, jb)
        self.assertEqual(a.stdout, b.stdout)


class ReportGuardTests(unittest.TestCase):
    """Hermetic: explicit nonexistent inputs and a scratch output directory, so that a broken guard can neither compute nor write a real result."""

    def run_report(self, *extra):
        d = tempfile.mkdtemp(prefix="r2_guard_")
        out = os.path.join(d, "out")
        p = subprocess.run([sys.executable, "-B", os.path.join(HERE, "t1_gorilla_pairs.py"), "report", "--sedef", "/nonexistent/s", "--genes", "/nonexistent/g",
                            "--truth", "/nonexistent/t", "--pairs", "/nonexistent/p", "--liftoff", "/nonexistent/l", "--outdir", out] + list(extra),
                           capture_output=True, text=True)
        return p, out

    def test_the_report_command_refuses_without_a_release_file(self):
        p, out = self.run_report()
        self.assertEqual(p.returncode, 2, p.stderr)
        self.assertIn("release", (p.stdout + p.stderr).lower())
        self.assertFalse(os.path.exists(out))

    def test_a_wrong_or_near_release_file_is_refused_too(self):
        with tempfile.TemporaryDirectory() as d:
            for text in ("not the token\n", T.RELEASE_TOKEN[:-1] + "X\n", " " + "\n" + T.RELEASE_TOKEN + "\n"):
                f = os.path.join(d, "rel")
                with open(f, "w") as fh:
                    fh.write(text)
                p, out = self.run_report("--release-file", f)
                self.assertEqual(p.returncode, 2, (text, p.stderr))
                self.assertFalse(os.path.exists(out))


class ReuseAndReportTests(unittest.TestCase):
    def test_prepare_tau_reuses_finished_outputs_when_the_retained_table_is_unchanged(self):
        lines = [sedef_row("NC_073224.2", 0, 5000, "NC_073227.2", 100000, 105000, frac=0.96)]
        rows, _ = T.read_sedef(lines)
        d = tempfile.mkdtemp(prefix="r2_reuse_")
        first = T.prepare_tau(rows, 0.90, d, reuse=True)
        self.assertFalse(first["reused"])
        again = T.prepare_tau(rows, 0.90, d, reuse=True)
        self.assertTrue(again["reused"])
        self.assertEqual((first["atoms"], first["classes"]), (again["atoms"], again["classes"]))
        rows2, _ = T.read_sedef(lines + [sedef_row("NC_073224.2", 20000, 26000, "NC_073227.2", 200000, 206000, frac=0.95)])
        changed = T.prepare_tau(rows2, 0.90, d, reuse=True)
        self.assertFalse(changed["reused"])
        self.assertEqual(changed["atoms"], 4)

    def test_the_report_runs_on_the_made_up_world_once_it_is_released(self):
        d, paths = world_with_control()
        reg, rel = release_files(d, paths)
        out = os.path.join(d, "out")
        rc, _text = run_main(report_argv(paths, reg, rel, out))
        self.assertEqual(rc, 0)
        rep = load_json(os.path.join(out, "report.json"))
        self.assertTrue(rep["valid"])
        t = rep["tau"]["0.90"]
        self.assertEqual((t["n_pairs"], t["n_families"], t["shared"]), (2, 2, 2))
        self.assertEqual(t["pair_weighted"], 1.0)
        self.assertEqual(t["a_viol"], 0.0)
        self.assertEqual(t["b_viol"], [4, 12])
        self.assertTrue(t["underpowered"])
        self.assertEqual(rep["tau"]["0.98"]["n_pairs"], 0)
        rows = read_text(os.path.join(out, "report_pairs_0.90.tsv")).splitlines()
        self.assertEqual(len(rows), 3)          # header + the two depth-matched pairs
        self.assertIn("path", rows[0])


class ControlAndGateCommandTests(unittest.TestCase):
    """counts with a self-lift table (the positive control branch) and gate0 on a made-up world."""

    def world_with_control(self, scatter_one=False):
        return world_with_control(scatter_one)

    def counts(self, d, paths):
        out = os.path.join(d, "out")
        rc = T.main(["counts", "--sedef", paths["sedef"], "--genes", paths["genes"], "--truth", paths["truth"], "--pairs", paths["pairs"],
                     "--liftoff", paths["liftoff"], "--outdir", out])
        return rc, load_json(os.path.join(out, "counts.json"))

    def test_the_positive_control_passes_when_all_eight_copies_sit_in_one_class(self):
        d, paths = self.world_with_control()
        rc, c = self.counts(d, paths)
        self.assertEqual(rc, 0)
        pc = c["positive_control"]
        self.assertTrue(pc["passed"])
        self.assertTrue(pc["connected"])
        self.assertEqual(pc["without_class"], [])
        self.assertEqual(pc["identity_range"], [0.97, 0.97])

    def test_the_positive_control_fails_when_a_copy_touches_no_class(self):
        d, paths = self.world_with_control(scatter_one=True)
        rc, c = self.counts(d, paths)
        self.assertEqual(rc, 0)
        self.assertFalse(c["positive_control"]["passed"])
        self.assertEqual(c["positive_control"]["without_class"], ["L7"])

    def test_gate0_reports_a_mismatch_against_the_registered_prefixes_and_fails(self):
        d, paths = self.world_with_control()
        out = os.path.join(d, "out")
        rc = T.main(["gate0", "--sedef", paths["sedef"], "--genes", paths["genes"], "--truth", paths["truth"], "--pairs", paths["pairs"],
                     "--liftoff", paths["liftoff"], "--outdir", out])
        self.assertEqual(rc, 1)
        tsv = read_text(os.path.join(out, "gate0_r2.tsv"))
        self.assertIn("dna_sd_atoms.py", tsv)
        self.assertIn("False", tsv)
        self.assertIn("True", tsv)       # the restored script itself matches its registered hash


class PrepareCommandAndPlantedRateTests(unittest.TestCase):
    def test_the_prepare_command_writes_one_summary_per_tau(self):
        d = tempfile.mkdtemp(prefix="r2_prep_")
        paths = make_world(d)
        out = os.path.join(d, "out")
        rc = T.main(["prepare", "--sedef", paths["sedef"], "--genes", paths["genes"], "--truth", paths["truth"], "--pairs", paths["pairs"], "--outdir", out])
        self.assertEqual(rc, 0)
        j = load_json(os.path.join(out, "prepare.json"))
        self.assertEqual(sorted(j["tau"]), ["0.90", "0.95", "0.98"])
        self.assertEqual((j["tau"]["0.90"]["atoms"], j["tau"]["0.90"]["classes"]), (4, 2))
        self.assertEqual(j["sedef"]["duplicates_removed"], 1)
        again = T.main(["prepare", "--sedef", paths["sedef"], "--genes", paths["genes"], "--truth", paths["truth"], "--pairs", paths["pairs"], "--outdir", out])
        self.assertEqual(again, 0)
        self.assertTrue(load_json(os.path.join(out, "prepare.json"))["tau"]["0.90"]["reused"])

    def test_planted_rates_all_shared_and_none_shared(self):
        U = {1: frozenset({1}), 2: frozenset({1}), 3: frozenset({2}), 4: frozenset({3})}
        r = T.s_rate([("A", 1, 2), ("B", 1, 2)], U)
        self.assertEqual((r["pair_weighted"], r["family_weighted"], r["a_viol"]), (1.0, 1.0, 0.0))
        r = T.s_rate([("A", 1, 3), ("B", 3, 4)], U)
        self.assertEqual((r["pair_weighted"], r["family_weighted"], r["a_viol"]), (0.0, 0.0, 1.0))
        b = T.bootstrap_rates(r["by_family"])
        self.assertEqual(b["pair_weighted"], (0.0, 0.0))

    def test_an_empty_table_reports_no_rate_and_is_underpowered(self):
        t = T.report_table([], [], {}, 0.98)
        self.assertEqual((t["n_pairs"], t["n_families"], t["pair_weighted"], t["a_viol"]), (0, 0, None, None))
        self.assertIsNone(t["bootstrap"])
        self.assertTrue(t["underpowered"])


class PathFlagAndPairTableTests(unittest.TestCase):
    def setUp(self):
        self.genes, self.index, _ = T.read_genes(GENES_TXT.splitlines())
        self.byname = {g.name: i for i, g in enumerate(self.genes)}

    def test_pairs_with_a_gene_that_touches_two_or_more_classes_are_flagged_path(self):
        g0, g1, g2 = self.byname["g0"], self.byname["g1"], self.byname["g2"]
        U = {g0: frozenset({7, 8}), g1: frozenset({7}), g2: frozenset({9})}
        dm = [{"family_id": "A", "g": g0, "h": g1, "identity": 0.95, "cov_short": 0.9, "a_gene": "gene-g0", "b_gene": "gene-g1", "dist_class": "same<100kb"},
              {"family_id": "B", "g": g1, "h": g2, "identity": 0.92, "cov_short": 0.3, "a_gene": "gene-g1", "b_gene": "gene-g2", "dist_class": "same100kb-1Mb"}]
        rows = T.pair_rows(dm, U, self.genes)
        self.assertEqual([r["S"] for r in rows], [1, 0])
        self.assertEqual([r["path"] for r in rows], [1, 0])
        self.assertEqual((rows[0]["n_classes_a"], rows[0]["n_classes_b"]), (2, 1))
        self.assertEqual(rows[1]["dist_class"], "same100kb-1Mb")
        t = T.report_table(dm, [], U, 0.90)
        self.assertEqual(t["path_pairs"], 1)


if __name__ == "__main__":
    unittest.main()
