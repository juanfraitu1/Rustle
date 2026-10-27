#!/usr/bin/env python3
"""Unit tests of t1_lib.py (docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md, KEY=dupfamboundary).
Run: python3 -B -m unittest bench/hierarchy/test_t1_lib.py
Chunks: 1 labels and unit strings, 2 statistics, 3 strata and the null, 4 classes and T1-e."""
import os
import random
import sys
import unittest
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import t1_lib as T  # noqa: E402


# ---------------------------------------------------------------- chunk 1: labels and unit strings
def colors(rows):
    return T.Colors.from_rows(rows)


class IntervalTests(unittest.TestCase):
    def test_touching_and_overlapping_intervals_merge(self):
        self.assertEqual(T.merge_intervals([[5, 10], [1, 3], [3, 4], [9, 12]]), [[1, 4], [5, 12]])

    def test_a_gap_keeps_intervals_apart(self):
        self.assertEqual(T.merge_intervals([[0, 5], [6, 9]]), [[0, 5], [6, 9]])

    def test_the_input_is_not_changed(self):
        iv = [[3, 4], [1, 3]]
        T.merge_intervals(iv)
        self.assertEqual(iv, [[3, 4], [1, 3]])


class CompositionTests(unittest.TestCase):
    def test_bases_are_the_overlap_of_every_segment_with_the_exon(self):
        c = colors([("c", 0, 60, "SD1"), ("c", 40, 200, "SD2")])
        self.assertEqual(T.composition(c, "c", 0, 100), Counter({"SD1": 60, "SD2": 60}))

    def test_same_label_overlaps_add(self):
        # two overlapping segments of one label both count (the comp() of soto_m2_unit_cover_cn.py)
        c = colors([("c", 0, 100, "X"), ("c", 50, 150, "X")])
        self.assertEqual(T.composition(c, "c", 0, 100)["X"], 150)

    def test_union_per_label_counts_a_base_once(self):
        c = colors([("c", 0, 100, "X"), ("c", 50, 150, "X")])
        self.assertEqual(T.composition(c, "c", 0, 100, union=True)["X"], 100)

    def test_another_chromosome_and_no_overlap_give_nothing(self):
        c = colors([("c", 0, 100, "X"), ("d", 0, 100, "Y")])
        self.assertEqual(T.composition(c, "c", 100, 200), Counter())
        self.assertEqual(T.composition(c, "e", 0, 100), Counter())

    def test_a_long_early_segment_is_not_missed(self):
        # the lookup must start maxlen before the exon: a segment that begins far upstream still overlaps
        c = colors([("c", 0, 10000, "LONG"), ("c", 9000, 9100, "S")])
        self.assertEqual(T.composition(c, "c", 9500, 9600), Counter({"LONG": 100}))


class DominantTests(unittest.TestCase):
    def test_most_bases_wins(self):
        self.assertEqual(T.pick_dominant(Counter({"A": 10, "B": 30}), "lex"), frozenset({"B"}))

    def test_a_tie_goes_to_the_lexicographically_smallest_label_string(self):
        comp = Counter({"SD10": 5, "SD9": 5, "SD3": 4})
        self.assertEqual(T.pick_dominant(comp, "lex"), frozenset({"SD10"}))   # 'SD10' < 'SD9' as strings

    def test_numeric_id_order_is_a_sensitivity(self):
        comp = Counter({"SD10": 5, "SD9": 5, "SD3": 4})
        self.assertEqual(T.pick_dominant(comp, "num"), frozenset({"SD9"}))

    def test_a_tie_can_yield_the_set_of_tied_labels(self):
        comp = Counter({"SD10": 5, "SD9": 5, "SD3": 4})
        self.assertEqual(T.pick_dominant(comp, "set"), frozenset({"SD10", "SD9"}))

    def test_the_near_tie_band_keeps_labels_within_five_percent_of_the_top(self):
        comp = Counter({"A": 100, "B": 95, "C": 94})
        self.assertEqual(T.pick_dominant(comp, "band"), frozenset({"A", "B"}))

    def test_no_composition_gives_no_label(self):
        for tie in ("lex", "num", "set", "band"):
            self.assertEqual(T.pick_dominant(Counter(), tie), frozenset())


class GeneStringTests(unittest.TestCase):
    def setUp(self):
        self.c = colors([
            ("c", 0, 100, "A"),                      # exon 1: A 100
            ("c", 200, 260, "B"), ("c", 260, 300, "A"),   # exon 2: B 60, A 40
        ])
        self.ex = [[0, 100], [200, 300], [400, 500]]  # exon 3 overlaps nothing

    def test_u1_is_the_set_of_dominant_labels_and_unlabelled_exons_drop_out(self):
        s = T.gene_strings(self.ex, "c", self.c)
        self.assertEqual(s.u1, frozenset({"A", "B"}))
        self.assertEqual(s.exon_dom, [frozenset({"A"}), frozenset({"B"}), frozenset()])

    def test_u0_is_the_label_with_most_bases_pooled_over_the_gene(self):
        s = T.gene_strings(self.ex, "c", self.c)
        self.assertEqual(s.pooled, Counter({"A": 140, "B": 60}))
        self.assertEqual(s.u0, "A")

    def test_the_ordered_string_is_genomic_order_with_repeats_collapsed(self):
        s = T.gene_strings(self.ex, "c", self.c)
        self.assertEqual(s.ordered, ("A", "B"))
        c2 = colors([("c", 0, 100, "A"), ("c", 200, 300, "A")])
        self.assertEqual(T.gene_strings([[0, 100], [200, 300]], "c", c2).ordered, ("A",))

    def test_exons_are_merged_before_labelling(self):
        c = colors([("c", 0, 100, "A")])
        s = T.gene_strings([[0, 50], [50, 100]], "c", c)   # touching exons form one merged exon
        self.assertEqual(s.n_exons, 1)
        self.assertEqual(s.pooled, Counter({"A": 100}))

    def test_a_gene_without_a_duplicon_has_empty_strings(self):
        s = T.gene_strings([[1000, 1100]], "c", self.c)
        self.assertIsNone(s.u0)
        self.assertEqual(s.u1, frozenset())
        self.assertEqual(s.ordered, ())

    def test_the_principal_duplicon_need_not_be_dominant_on_any_exon(self):
        # B is never the top label of an exon but has the most bases pooled: U0 not in U1 (10 such genes in V)
        c = colors([
            ("c", 0, 60, "A"), ("c", 0, 59, "B"),
            ("c", 200, 270, "C"), ("c", 200, 259, "B"),
            ("c", 400, 470, "D"), ("c", 400, 459, "B"),
        ])
        s = T.gene_strings([[0, 100], [200, 300], [400, 500]], "c", c)
        self.assertEqual(s.u1, frozenset({"A", "C", "D"}))
        self.assertEqual(s.u0, "B")

    def test_any_overlap_lists_every_label_with_a_base_on_an_exon(self):
        s = T.gene_strings(self.ex, "c", self.c)
        self.assertEqual(s.any_overlap, frozenset({"A", "B"}))

    def test_sensitivity_modes_change_only_what_they_name(self):
        c = colors([("c", 0, 50, "SD10"), ("c", 50, 100, "SD9")])
        lex = T.gene_strings([[0, 100]], "c", c, tie="lex")
        num = T.gene_strings([[0, 100]], "c", c, tie="num")
        st = T.gene_strings([[0, 100]], "c", c, tie="set")
        self.assertEqual((lex.u1, num.u1, st.u1), (frozenset({"SD10"}), frozenset({"SD9"}), frozenset({"SD10", "SD9"})))
        self.assertEqual((lex.u0, num.u0), ("SD10", "SD9"))


class CopyNumberClassTests(unittest.TestCase):
    def test_single_linkage_cuts_where_consecutive_values_differ_by_two_or_more(self):
        self.assertEqual(T.cn_classes([5.0, 1.0, 2.9, 3.0]), [[1.0, 2.9, 3.0], [5.0]])

    def test_a_gap_of_exactly_two_starts_a_new_class_like_the_frozen_comparator(self):
        self.assertEqual(T.cn_classes([1.0, 3.0]), [[1.0], [3.0]])

    def test_empty_input(self):
        self.assertEqual(T.cn_classes([]), [])

    def test_u2_pairs_each_label_with_the_class_of_the_genes_own_copy_number(self):
        u1 = {"g1": frozenset({"A"}), "g2": frozenset({"A"}), "g3": frozenset({"A", "B"}), "g4": frozenset({"B"})}
        cn = {"g1": 1.0, "g2": 6.0, "g3": 1.5, "g4": None}
        u2 = T.u2_strings(u1, cn)
        self.assertEqual(u2["g1"], frozenset({("A", 0)}))
        self.assertEqual(u2["g2"], frozenset({("A", 1)}))
        self.assertEqual(u2["g3"], frozenset({("A", 0), ("B", 0)}))
        self.assertNotIn("g4", u2)          # a gene without a copy number is outside U2

    def test_u2_refines_u1_by_dropping_the_class(self):
        rng = random.Random(7)
        u1 = {f"g{i}": frozenset(rng.sample("ABCDE", rng.randint(1, 3))) for i in range(60)}
        cn = {g: round(rng.uniform(0, 12), 1) for g in u1}
        u2 = T.u2_strings(u1, cn)
        for g, s in u2.items():
            self.assertEqual(frozenset(d for d, _ in s), u1[g])


# ---------------------------------------------------------------- chunk 2: statistics
def fs(*x):
    return frozenset(x)


class DeficiencyTests(unittest.TestCase):
    def setUp(self):
        self.fam = {"g1": "f1", "g2": "f1", "g3": "f1", "g4": "f2", "g5": "f2"}
        self.str = {"g1": fs("a"), "g2": fs("a"), "g3": fs("b"), "g4": fs("b"), "g5": fs("b")}

    def test_m1_counts_genes_outside_the_plurality_string_of_each_family(self):
        self.assertEqual(T.m1(self.fam, self.str), 1)

    def test_m2_counts_genes_outside_the_best_family_of_each_string(self):
        self.assertEqual(T.m2(self.fam, self.str), 1)

    def test_m1c_equals_m1_when_no_string_is_nested(self):
        self.assertEqual(T.m1c(self.fam, self.str), 1)

    def test_nesting_costs_nothing_in_m1c_but_counts_in_m1(self):
        fam = {"g1": "f", "g2": "f", "g3": "f"}
        st = {"g1": fs("a"), "g2": fs("a", "b"), "g3": fs("a", "b", "c")}
        self.assertEqual(T.m1(fam, st), 2)
        self.assertEqual(T.m1c(fam, st), 0)

    def test_a_branch_costs_one_gene(self):
        fam = {"g1": "f", "g2": "f", "g3": "f"}
        st = {"g1": fs("a"), "g2": fs("a", "b"), "g3": fs("a", "c")}
        self.assertEqual(T.m1c(fam, st), 1)

    def test_the_chain_with_most_genes_wins(self):
        fam = {f"g{i}": "f" for i in range(6)}
        st = {"g0": fs("a"), "g1": fs("a", "b"), "g2": fs("a", "b"), "g3": fs("a", "b"), "g4": fs("a", "c"), "g5": fs("a", "c")}
        self.assertEqual(T.m1c(fam, st), 2)          # chain a < ab holds 4 genes; a < ac holds 3

    def test_families_and_strings_of_one_gene_are_not_scored(self):
        fam = {"g1": "f1", "g2": "f2", "g3": "f3"}
        st = {"g1": fs("a"), "g2": fs("b"), "g3": fs("c")}
        self.assertEqual((T.m1(fam, st), T.m1c(fam, st), T.m2(fam, st)), (0, 0, 0))

    def test_the_universe_argument_restricts_every_count(self):
        self.assertEqual(T.m1(self.fam, self.str, genes=["g1", "g2", "g4", "g5"]), 0)
        self.assertEqual(T.m2(self.fam, self.str, genes=["g1", "g2", "g4", "g5"]), 0)

    def test_the_result_does_not_depend_on_the_input_order_of_the_genes(self):
        rng = random.Random(3)
        genes = sorted(self.fam)
        for _ in range(20):
            rng.shuffle(genes)
            self.assertEqual((T.m1(self.fam, self.str, genes), T.m1c(self.fam, self.str, genes), T.m2(self.fam, self.str, genes)), (1, 1, 1))

    def test_the_chain_dynamic_programme_equals_brute_force_on_random_families(self):
        rng = random.Random(11)
        for _ in range(300):
            n = rng.randint(2, 9)
            fam = {f"g{i}": "f" for i in range(n)}
            st = {g: frozenset(rng.sample("ABCD", rng.randint(1, 3))) for g in fam}
            cnt = Counter(st.values())
            self.assertEqual(T.chain_weight(cnt), brute_chain_weight(cnt))
            self.assertEqual(T.m1c(fam, st), n - brute_chain_weight(cnt))

    def test_families_equal_to_the_string_partition_give_zero(self):
        rng = random.Random(5)
        st = {f"g{i}": fs(rng.choice("ABCDEF")) for i in range(50)}
        fam = {g: "".join(s) for g, s in st.items()}
        self.assertEqual((T.m1(fam, st), T.m1c(fam, st), T.m2(fam, st)), (0, 0, 0))

    def test_relabelling_k_genes_costs_exactly_k_when_the_native_string_stays_the_plurality(self):
        rng = random.Random(9)
        for _ in range(100):
            nfam = rng.randint(3, 6)
            sizes = [rng.randint(5, 9) for _ in range(nfam)]
            fam, st = {}, {}
            for i, n in enumerate(sizes):
                for j in range(n):
                    fam[f"g{i}_{j}"] = f"f{i}"
                    st[f"g{i}_{j}"] = fs(f"s{i}")
            k = rng.randint(1, 3)
            for g in rng.sample(sorted(fam), k):
                own = int(fam[g][1:])
                st[g] = fs(f"s{rng.choice([x for x in range(nfam) if x != own])}")
            # at most 3 genes move, every family has >= 5 genes, so each native string keeps the plurality of its family and of its string
            self.assertEqual((T.m1(fam, st), T.m2(fam, st)), (k, k))

    def test_refining_the_strings_never_lowers_m1_and_never_raises_m2(self):
        rng = random.Random(21)
        for _ in range(200):
            n = rng.randint(6, 30)
            genes = [f"g{i}" for i in range(n)]
            fam = {g: f"f{rng.randint(0, 4)}" for g in genes}
            coarse = {g: fs(rng.choice("ABC")) for g in genes}
            fine = {g: fs((next(iter(coarse[g])), rng.randint(0, 2))) for g in genes}
            self.assertLessEqual(T.m1(fam, coarse), T.m1(fam, fine))
            self.assertGreaterEqual(T.m2(fam, coarse), T.m2(fam, fine))

    def test_synthetic_families_as_one_string_nested_truncations_and_branches(self):
        one = {f"g{i}": fs("a", "b") for i in range(5)}
        nested = {f"g{i}": fs(*"abcdefgh"[: i + 1]) for i in range(5)}
        branch = {"g0": fs("a"), "g1": fs("a", "b"), "g2": fs("a", "c"), "g3": fs("a", "d")}
        for st, expect_zero in ((one, True), (nested, True), (branch, False)):
            fam = {g: "f" for g in st}
            self.assertEqual(T.m1c(fam, st) == 0, expect_zero)

    def test_macro_versions_average_the_per_family_and_per_string_shares(self):
        # f1 has 3 genes with deficit 1 (share 1/3), f2 has 2 genes and no deficit
        self.assertAlmostEqual(T.macro_d1(self.fam, self.str), (1 / 3 + 0) / 2)
        # strings: a (2 genes, share 0), b (3 genes, share 1/3)
        self.assertAlmostEqual(T.macro_d2(self.fam, self.str), (0 + 1 / 3) / 2)

    def test_macro_of_an_empty_scope_is_zero(self):
        self.assertEqual(T.macro_d1({}, {}), 0.0)
        self.assertEqual(T.macro_d2({}, {}), 0.0)


def brute_chain_weight(cnt):
    """Maximum total count over chains of distinct strings under proper-subset inclusion, by enumeration."""
    items = sorted(cnt, key=lambda s: (len(s), sorted(s)))
    best = 0
    n = len(items)
    for mask in range(1 << n):
        chosen = [items[i] for i in range(n) if mask >> i & 1]
        ok = all(a < b for a, b in zip(chosen, chosen[1:]))      # sorted by size, each a proper subset of the next
        if ok:
            best = max(best, sum(cnt[s] for s in chosen))
    return best


class EndTruncationTests(unittest.TestCase):
    def test_m1t_counts_genes_outside_the_largest_pairwise_comparable_set(self):
        fam = {f"g{i}": "f" for i in range(5)}
        od = {"g0": ("A", "B", "C"), "g1": ("A", "B"), "g2": ("B", "C"), "g3": ("A",), "g4": ("C",)}
        self.assertEqual(T.m1t(fam, od), 2)

    def test_comparability_is_pairwise_not_chained(self):
        # BC is a prefix of BCD and BCD is a suffix of ABCD, but BC is neither a prefix nor a suffix of ABCD
        fam = {"g0": "f", "g1": "f", "g2": "f"}
        od = {"g0": ("A", "B", "C", "D"), "g1": ("B", "C"), "g2": ("B", "C", "D")}
        self.assertEqual(T.m1t(fam, od), 1)

    def test_equal_strings_are_comparable(self):
        fam = {"g0": "f", "g1": "f", "g2": "f"}
        od = {"g0": ("A", "B"), "g1": ("A", "B"), "g2": ("A", "B")}
        self.assertEqual(T.m1t(fam, od), 0)


class PathDefectTests(unittest.TestCase):
    def defect(self, strings):
        fam = {f"g{i}": "f" for i in range(len(strings))}
        od = {f"g{i}": s for i, s in enumerate(strings)}
        return T.path_defect(fam, od, {g: frozenset(s) for g, s in od.items()})

    def test_a_path_has_defect_zero(self):
        self.assertEqual(self.defect([("A", "B", "C")])[0], 0)
        self.assertEqual(self.defect([("A", "B"), ("B", "C")])[0], 0)

    def test_two_disjoint_units_have_defect_one(self):
        self.assertEqual(self.defect([("A",), ("B",)])[0], 1)

    def test_a_three_leaf_star_has_defect_one(self):
        self.assertEqual(self.defect([("A", "B"), ("A", "C"), ("A", "D")])[0], 1)

    def test_a_triangle_has_defect_one(self):
        self.assertEqual(self.defect([("A", "B"), ("B", "C"), ("C", "A")])[0], 1)

    def test_repeated_pairs_do_not_make_multi_edges(self):
        self.assertEqual(self.defect([("A", "B"), ("A", "B"), ("B", "A")])[0], 0)

    def test_ownership_multiplicity_counts_the_extra_families_that_hold_a_unit(self):
        fam = {"g0": "f1", "g1": "f1", "g2": "f2", "g3": "f2", "g4": "f3"}
        od = {"g0": ("A", "B"), "g1": ("A", "B"), "g2": ("B", "C"), "g3": ("B", "C"), "g4": ("Z",)}
        strings = {g: frozenset(s) for g, s in od.items()}
        delta, mult = T.path_defect(fam, od, strings)
        # unit B is held by f1 and f2 (mu = 2); f3 has one gene and is not scored
        self.assertEqual(mult, 1)


class AuditTests(unittest.TestCase):
    def setUp(self):
        self.fam = {"g1": "f", "g2": "f", "g3": "f"}
        self.str = {"g1": fs("a"), "g2": fs("a", "b"), "g3": fs("a", "c")}
        self.info = {"g1": dict(nex=1, interior=fs()), "g2": dict(nex=2, interior=fs()), "g3": dict(nex=5, interior=fs("x"))}

    def test_ties_between_equal_weight_chains_go_to_the_smallest_key(self):
        rec = T.audit_d1(self.fam, self.str, self.info)["f"]
        self.assertEqual(rec["kept_chain"], [("a",), ("a", "b")])
        self.assertEqual(rec["deficit_genes"], ["g3"])
        self.assertEqual(len(rec["tied_chains"]), 2)
        self.assertEqual((rec["n_tied_chains"], rec["tied_truncated"]), (2, False))
        self.assertEqual(sorted(len(c["deficit_genes"]) for c in rec["tied_chains"]), [1, 1])

    def brute_chains(self, cnt):
        items = sorted(cnt, key=lambda x: (len(x), tuple(sorted(x))))
        out = []

        def grow(chain):
            out.append(list(chain))
            for t in items:
                if t > chain[-1]:
                    grow(chain + [t])
        for s0 in items:
            grow([s0])
        w = lambda c: sum(cnt[x] for x in c)
        top = max(w(c) for c in out)
        best = [c for c in out if w(c) == top]
        best.sort(key=lambda c: [tuple(sorted(x)) for x in c])
        return top, best

    def test_the_chain_search_equals_a_brute_force_enumeration_on_random_families(self):
        rng = random.Random(21)
        for _ in range(60):
            cnt = {}
            for _ in range(rng.randint(1, 9)):
                cnt[frozenset(rng.sample("ABCD", rng.randint(1, 3)))] = rng.randint(1, 3)
            top, chains, n_tied = T._best_chains(cnt)
            btop, bchains = self.brute_chains(cnt)
            self.assertEqual(top, btop)
            self.assertEqual(chains, bchains)
            self.assertEqual(n_tied, len(bchains))

    def test_a_family_with_hundreds_of_tied_chains_is_enumerated_only_up_to_the_cap_and_the_kept_chain_is_still_the_smallest_key(self):
        import itertools
        labels = "ABCDE"
        cnt = {frozenset(c): 1 for k in range(1, 6) for c in itertools.combinations(labels, k)}
        top, chains, n_tied = T._best_chains(cnt, cap=10)
        self.assertEqual((top, n_tied, len(chains)), (5, 120, 10))                 # 5! maximal chains from a singleton to the full set
        self.assertEqual(chains[0], [frozenset("A"), frozenset("AB"), frozenset("ABC"), frozenset("ABCD"), frozenset("ABCDE")])
        self.assertEqual(chains, sorted(chains, key=lambda c: [tuple(sorted(x)) for x in c]))
        fam = {f"g{i}": "F" for i in range(len(cnt))}
        st = dict(zip(fam, cnt))
        rec = T.audit_d1(fam, st, {g: dict(nex=3, interior=fs()) for g in fam}, cap=10)["F"]
        self.assertEqual((len(rec["tied_chains"]), rec["n_tied_chains"], rec["tied_truncated"]), (10, 120, True))
        self.assertEqual(rec["kept_chain"], [("A",), ("A", "B"), ("A", "B", "C"), ("A", "B", "C", "D"), ("A", "B", "C", "D", "E")])

    def test_label_short_for_a_gene_with_at_most_two_merged_exons(self):
        info = dict(self.info)
        info["g3"] = dict(nex=2, interior=fs("x"))
        self.assertEqual(T.audit_d1(self.fam, self.str, info)["f"]["labels"]["g3"], "short")

    def test_label_ends_only_when_the_interior_equals_that_of_a_kept_gene_with_three_exons(self):
        # g1 and g2 are kept; g1 has three or more exons and the same interior set as g3
        info = {"g1": dict(nex=4, interior=fs("x")), "g2": dict(nex=4, interior=fs("y")), "g3": dict(nex=5, interior=fs("x"))}
        self.assertEqual(T.audit_d1(self.fam, self.str, info)["f"]["labels"]["g3"], "ends-only")

    def test_label_unexplained_is_the_rest(self):
        info = {"g1": dict(nex=4, interior=fs("z")), "g2": dict(nex=4, interior=fs("y")), "g3": dict(nex=5, interior=fs("x"))}
        self.assertEqual(T.audit_d1(self.fam, self.str, info)["f"]["labels"]["g3"], "unexplained")

    def test_the_audit_never_changes_the_deficiency(self):
        total = sum(len(r["deficit_genes"]) for r in T.audit_d1(self.fam, self.str, self.info).values())
        self.assertEqual(total, T.m1c(self.fam, self.str))

    def test_d2_deficit_genes_are_those_outside_the_best_family_with_ties_to_the_smaller_id(self):
        fam = {"g1": "f2", "g2": "f1", "g3": "f1", "g4": "f2"}
        st = {g: fs("s") for g in fam}
        rec = T.audit_d2(fam, st, {g: dict(nex=3, interior=fs()) for g in fam})[("s",)]
        self.assertEqual(rec["best_family"], "f1")
        self.assertEqual(rec["deficit_genes"], ["g1", "g4"])
        self.assertEqual(rec["tied_best_families"], ["f1", "f2"])           # the tie is listed so that its break can be audited

    def test_d2_lists_a_single_best_family_when_there_is_no_tie(self):
        fam = {"g1": "f2", "g2": "f1", "g3": "f1"}
        rec = T.audit_d2(fam, {g: fs("s") for g in fam}, {g: dict(nex=3, interior=fs()) for g in fam})[("s",)]
        self.assertEqual((rec["best_family"], rec["tied_best_families"]), ("f1", ["f1"]))

    def test_exception_lists_are_identical_under_any_input_order(self):
        rng = random.Random(4)
        genes = sorted(self.fam)
        ref = T.exceptions_d1(self.fam, self.str, self.info), T.exceptions_d2(self.fam, self.str, self.info)
        for _ in range(10):
            rng.shuffle(genes)
            fam = {g: self.fam[g] for g in genes}
            st = {g: self.str[g] for g in genes}
            info = {g: self.info[g] for g in genes}
            self.assertEqual((T.exceptions_d1(fam, st, info), T.exceptions_d2(fam, st, info)), ref)

    def test_the_exception_counts_equal_the_statistics(self):
        rng = random.Random(8)
        for _ in range(50):
            n = rng.randint(6, 30)
            genes = [f"g{i:02d}" for i in range(n)]
            fam = {g: f"f{rng.randint(0, 3)}" for g in genes}
            st = {g: frozenset(rng.sample("ABC", rng.randint(1, 2))) for g in genes}
            info = {g: dict(nex=rng.randint(1, 6), interior=fs()) for g in genes}
            self.assertEqual(len(T.exceptions_d1(fam, st, info)), T.m1c(fam, st))
            self.assertEqual(len(T.exceptions_d2(fam, st, info)), T.m2(fam, st))


# ---------------------------------------------------------------- chunk 3: strata, the null, classes
def world(seed, n_strata=12, mode="independent", sizes=None, nclu=(2, 4)):
    return T.synthetic_world(seed, n_strata=n_strata, mode=mode, sizes=sizes, nclu=nclu)


def merge_worlds(a, b):
    """Two synthetic worlds as one: the genes, labels and strata of b get the prefix 'x' so that nothing collides."""
    fam, strings, info = dict(a["fam"]), dict(a["strings"]), dict(a["info"])
    strata = list(a["strata"])
    for g, f in b["fam"].items():
        fam["x" + g] = "x" + f
        strings["x" + g] = b["strings"][g]
        info["x" + g] = b["info"][g]
    for sid, genes in b["strata"]:
        strata.append(("x" + sid, ["x" + g for g in genes]))
    strata.sort()
    return dict(fam=fam, strings=strings, strata=strata, info=info)


def anti_aligned_world(n):
    """n strata of four genes: clusters {g1,g2} and {g3,g4}, strings x, y, x, y. The real cut is worse than the mean random re-cut, so no
    z is negative in either cell."""
    fam, strings, strata, info = {}, {}, [], {}
    for i in range(n):
        genes = [f"a{i:02d}_{j}" for j in range(4)]
        for j, g in enumerate(genes):
            fam[g] = f"A{i:02d}_{j // 2}"
            strings[g] = frozenset({"x" if j % 2 == 0 else "y"})
            info[g] = dict(nex=3, interior=frozenset())
        strata.append((genes[0], genes))
    return dict(fam=fam, strings=strings, strata=strata, info=info)


def chain_world(n):
    """n strata of six genes in two clusters of three whose strings {a} < {a,b} < {a,b,c} form one inclusion chain: m1c cannot move,
    m2 can."""
    fam, strings, strata, info = {}, {}, [], {}
    chain = [frozenset("a"), frozenset("a"), frozenset("ab"), frozenset("ab"), frozenset("abc"), frozenset("abc")]
    for i in range(n):
        genes = [f"c{i:02d}_{j}" for j in range(6)]
        for j, g in enumerate(genes):
            fam[g] = f"C{i:02d}_{j % 2}"
            strings[g] = chain[j]
            info[g] = dict(nex=3, interior=frozenset())
        strata.append((genes[0], genes))
    return dict(fam=fam, strings=strings, strata=strata, info=info)


class StrataTests(unittest.TestCase):
    def test_components_are_labelled_by_their_smallest_node_string(self):
        comp = T.graph_components([("b", "c", 1.0), ("a", "b", 1.0), ("e", "d", 1.0)])
        self.assertEqual(comp, {"a": "a", "b": "a", "c": "a", "d": "d", "e": "d"})

    def test_strata_are_ordered_by_their_smallest_gene_id_and_named_after_it(self):
        comp = {"n1": "A", "n2": "A", "n3": "B"}
        gene_node = {"g3": "n3", "g1": "n1", "g2": "n2", "gX": "n3"}
        strata = T.build_strata(["g3", "g1", "g2", "gX"], gene_node, comp)
        self.assertEqual(strata, [("g1", ["g1", "g2"]), ("g3", ["g3", "gX"])])

    def test_genes_without_a_node_are_outside_every_stratum(self):
        strata = T.build_strata(["g1", "g2"], {"g1": "n1", "g2": None}, {"n1": "A"})
        self.assertEqual(strata, [("g1", ["g1"])])

    def test_variation_bearing_strata_hold_at_least_two_labels(self):
        fam = {"g1": "f1", "g2": "f1", "g3": "f1", "g4": "f2"}
        strata = [("g1", ["g1", "g2"]), ("g3", ["g3", "g4"])]
        self.assertEqual(T.variation_strata(fam, strata), [("g3", ["g3", "g4"])])


class PValueAndClassTests(unittest.TestCase):
    def test_p_is_the_lower_tail_with_the_plus_one_rule(self):
        self.assertAlmostEqual(T.p_lower(3, [5, 4, 3, 2, 1]), (1 + 3) / 6)
        self.assertAlmostEqual(T.p_lower(0, [5, 4, 3]), 1 / 4)

    def test_class_thresholds_are_inclusive_at_001_and_005(self):
        self.assertEqual(T.content_class(0.01), "CUT-ALIGNED")
        self.assertEqual(T.content_class(0.0100001), "PARTIAL")
        self.assertEqual(T.content_class(0.05), "PARTIAL")
        self.assertEqual(T.content_class(0.0500001), "NOT-DISTINGUISHED")

    def test_precedence_invalid_underpowered_uninformative_underpowered_by_live_strata_then_content(self):
        c = T.classify_cut                                       # (p_all, p_removed, n_var, n_live)
        self.assertEqual(c(0.001, 0.001, 30, 30, valid=False)["class"], "INVALID")
        self.assertEqual(c(0.001, 0.001, 7, 7)["class"], "UNDERPOWERED")
        self.assertEqual(c(0.001, 0.001, 7, 0)["class"], "UNDERPOWERED")          # fewer than 8 strata comes before no live stratum
        self.assertEqual(c(0.001, 0.001, 30, 0)["class"], "UNINFORMATIVE")
        self.assertEqual(c(0.001, 0.001, 8, 0)["class"], "UNINFORMATIVE")
        self.assertEqual(c(0.001, 0.001, 30, 5)["class"], "UNDERPOWERED")         # at least 8 strata, fewer than 8 live
        self.assertEqual(c(0.001, 0.001, 30, 7)["class"], "UNDERPOWERED")
        self.assertEqual(c(0.001, 0.001, 30, 8)["class"], "CUT-ALIGNED")
        self.assertEqual(c(0.001, 0.001, 30, 30)["class"], "CUT-ALIGNED")

    def test_the_issued_class_is_the_lower_of_the_two_p_classes(self):
        self.assertEqual(T.classify_cut(0.0005, 0.04, 30, 30)["class"], "PARTIAL")
        self.assertEqual(T.classify_cut(0.0005, 0.4, 30, 30)["class"], "NOT-DISTINGUISHED")
        self.assertEqual(T.classify_cut(0.03, 0.0005, 30, 30)["class"], "PARTIAL")

    def test_without_a_removal_the_class_is_that_of_all_strata(self):
        r = T.classify_cut(0.0005, None, 30, 30)
        self.assertEqual((r["class"], r["class_all"], r["class_removed"]), ("CUT-ALIGNED", "CUT-ALIGNED", "CUT-ALIGNED"))

    def test_exactly_eight_independent_units_are_enough(self):
        self.assertEqual(T.classify_cut(0.5, 0.5, 8, 8)["class"], "NOT-DISTINGUISHED")
        self.assertEqual(T.classify_cut(0.5, 0.5, 7, 7)["class"], "UNDERPOWERED")
        self.assertEqual(T.classify_cut(0.5, 0.5, 8, 7)["class"], "UNDERPOWERED")


class NullEngineTests(unittest.TestCase):
    def test_the_observed_statistics_equal_the_direct_ones(self):
        for seed in range(5):
            w = world(seed)
            eng = T.NullEngine(w["fam"], w["strings"], w["strata"])
            obs = eng.observed()
            self.assertEqual(obs["m1c"], T.m1c(w["fam"], w["strings"]))
            self.assertEqual(obs["m2"], T.m2(w["fam"], w["strings"]))

    def test_every_draw_equals_a_direct_computation_on_the_permuted_labels(self):
        for seed in range(4):
            w = world(seed)
            eng = T.NullEngine(w["fam"], w["strings"], w["strata"])
            res = eng.draws(25)
            for b, fam_b in enumerate(eng.naive_draws(25)):
                self.assertEqual(res["tot_m1c"][b], T.m1c(fam_b, w["strings"]), (seed, b))
                self.assertEqual(res["tot_m2"][b], T.m2(fam_b, w["strings"]), (seed, b))

    def test_stratum_statistics_are_the_statistic_on_the_genes_of_that_stratum_alone(self):
        w = world(3)
        eng = T.NullEngine(w["fam"], w["strings"], w["strata"])
        res = eng.draws(10)
        for b, fam_b in enumerate(eng.naive_draws(10)):
            for j, (sid, genes) in enumerate(eng.var):
                self.assertEqual(res["st_m1c"][b][j], T.m1c(fam_b, w["strings"], genes))
                self.assertEqual(res["st_m2"][b][j], T.m2(fam_b, w["strings"], genes))

    def test_blocks_of_draws_continue_one_random_stream(self):
        w = world(2)
        eng = T.NullEngine(w["fam"], w["strings"], w["strata"])
        whole = eng.draws(40)
        a, b = eng.draws(40, start=0, end=15), eng.draws(40, start=15, end=40)
        self.assertEqual(whole["tot_m1c"], a["tot_m1c"] + b["tot_m1c"])
        self.assertEqual(whole["tot_m2"], a["tot_m2"] + b["tot_m2"])

    def test_the_registered_procedure_shuffles_the_observed_label_list_of_each_variation_stratum_in_order(self):
        w = world(1)
        eng = T.NullEngine(w["fam"], w["strings"], w["strata"])
        rng = random.Random(T.SEED_T1)
        expect = {}
        for sid, genes in eng.var:
            lab = [w["fam"][g] for g in genes]
            perm = lab[:]
            rng.shuffle(perm)
            expect.update(zip(genes, perm))
        first = next(eng.naive_draws(1))
        for g, lab in expect.items():
            self.assertEqual(first[g], lab)

    def test_cluster_sizes_are_kept_inside_each_stratum(self):
        w = world(6)
        eng = T.NullEngine(w["fam"], w["strings"], w["strata"])
        for fam_b in eng.naive_draws(10):
            for sid, genes in eng.var:
                self.assertEqual(Counter(fam_b[g] for g in genes), Counter(w["fam"][g] for g in genes))

    def test_a_family_that_spans_two_strata_is_refused(self):
        fam = {"g1": "f", "g2": "f", "g3": "h", "g4": "h"}
        st = {g: frozenset({"a"}) for g in fam}
        with self.assertRaises(ValueError):
            T.NullEngine(fam, st, [("g1", ["g1", "g3"]), ("g2", ["g2", "g4"])])

    def test_the_same_seed_gives_the_same_draws_and_another_seed_does_not(self):
        w = world(0)
        a = T.NullEngine(w["fam"], w["strings"], w["strata"]).draws(30)
        b = T.NullEngine(w["fam"], w["strings"], w["strata"]).draws(30)
        c = T.NullEngine(w["fam"], w["strings"], w["strata"], seed=1).draws(30)
        self.assertEqual(a["tot_m2"], b["tot_m2"])
        self.assertNotEqual(a["tot_m2"], c["tot_m2"])


class CellOutcomeTests(unittest.TestCase):
    """The synthetic controls of Gate 4 that concern the null (B is small here; the gate itself uses 10,000)."""

    B = 600

    def cell(self, w, stat):
        return T.run_cell(w["fam"], w["strings"], w["strata"], stat=stat, B=self.B)

    def test_a_family_partition_equal_to_the_strata_is_underpowered(self):
        w = world(0, mode="one_cluster_per_stratum")
        self.assertEqual(self.cell(w, "m1c")["class"], "UNDERPOWERED")

    def test_eight_strata_with_two_clusters_and_one_string_each_are_uninformative(self):
        w = world(0, n_strata=8, mode="one_string_per_stratum", nclu=(2, 2))
        for stat in ("m1c", "m2"):
            r = self.cell(w, stat)
            self.assertEqual(r["class"], "UNINFORMATIVE", stat)
            self.assertEqual((r["n_var"], r["n_live"]), (8, 0), stat)

    def test_ten_such_strata_are_uninformative_too(self):
        w = world(2, n_strata=10, mode="one_string_per_stratum", nclu=(2, 2))
        for stat in ("m1c", "m2"):
            self.assertEqual(self.cell(w, stat)["class"], "UNINFORMATIVE", stat)

    def test_eight_variation_bearing_strata_of_which_only_two_are_live_are_underpowered(self):
        # a small vocabulary, so that strings repeat inside a stratum and m2 can move there
        w = merge_worlds(T.synthetic_world(3, n_strata=2, mode="independent", sizes=[12, 12], nclu=(3, 3), vocab="ABC"),
                         world(4, n_strata=6, mode="one_string_per_stratum", nclu=(2, 2)))
        for stat in ("m1c", "m2"):
            r = self.cell(w, stat)
            self.assertEqual(r["n_var"], 8, stat)
            self.assertEqual(r["n_live"], 2, stat)
            self.assertEqual(r["class"], "UNDERPOWERED", stat)

    def test_a_stratum_whose_strings_form_one_inclusion_chain_is_dead_for_m1c_but_not_for_m2(self):
        w = chain_world(10)
        r1 = self.cell(w, "m1c")
        r2 = self.cell(w, "m2")
        self.assertEqual((r1["n_var"], r1["n_live"], r1["class"]), (10, 0, "UNINFORMATIVE"))
        self.assertEqual((r2["n_var"], r2["n_live"]), (10, 10))
        self.assertNotIn(r2["class"], ("UNINFORMATIVE", "UNDERPOWERED"))

    def test_the_live_strata_are_those_whose_own_statistic_takes_two_values_over_the_draws(self):
        w = world(7, n_strata=12, mode="independent", sizes=[10] * 12, nclu=(2, 3))
        eng = T.NullEngine(w["fam"], w["strings"], w["strata"])
        dr = eng.draws(200)
        for stat in ("m1c", "m2"):
            cols = [[row[j] for row in dr["st_" + stat]] for j in range(len(eng.var))]
            expect = sum(1 for c in cols if len(set(c)) >= 2)
            self.assertEqual(T.run_cell(w["fam"], w["strings"], w["strata"], stat=stat, B=200)["n_live"], expect)

    def test_seven_strata_of_that_kind_are_underpowered_first(self):
        w = world(0, n_strata=7, mode="one_string_per_stratum", nclu=(2, 2))
        self.assertEqual(self.cell(w, "m1c")["class"], "UNDERPOWERED")

    def test_a_partition_that_follows_the_strings_is_cut_aligned_in_both_cells(self):
        w = world(0, n_strata=14, mode="follow", sizes=[10] * 14)
        for stat in ("m1c", "m2"):
            r = self.cell(w, stat)
            self.assertLessEqual(r["p"], 0.01, stat)
            self.assertEqual(r["class"], "CUT-ALIGNED", stat)

    def test_strings_nested_with_the_cut_give_the_lowest_p_in_both_cells(self):
        for mode in ("merge", "split"):
            w = world(1, n_strata=14, mode=mode, sizes=[12] * 14, nclu=(3, 4))
            for stat in ("m1c", "m2"):
                self.assertLessEqual(self.cell(w, stat)["p"], 0.01, (mode, stat))

    def test_strings_independent_of_the_cut_are_rarely_significant(self):
        hits = 0
        for seed in range(30):
            w = world(100 + seed, n_strata=10, mode="independent", sizes=[10] * 10, nclu=(2, 3))
            r = T.run_cell(w["fam"], w["strings"], w["strata"], stat="m1c", B=200)
            hits += r["p"] <= 0.05
        self.assertLessEqual(hits, 8)       # about 1.5 expected out of 30


class RobustnessTests(unittest.TestCase):
    def one_aligned_stratum(self, seed=5):
        # one large variation-bearing stratum follows its strings; thirteen small ones carry strings independent of the cut
        w = world(seed, n_strata=14, mode="independent", sizes=[40] + [6] * 13, nclu=(2, 3))
        sid, genes = T.variation_strata(w["fam"], w["strata"])[0]
        for g in genes:
            w["strings"][g] = frozenset({"own_" + w["fam"][g]})
        return w, sid

    def test_one_aligned_stratum_is_found_as_the_one_that_carries_the_class(self):
        for seed in range(5, 11):
            w, sid = self.one_aligned_stratum(seed)
            r = T.run_cell(w["fam"], w["strings"], w["strata"], stat="m1c", B=400)
            self.assertEqual(r["removed_stratum"], sid, seed)
            self.assertLessEqual(r["p"], 0.05, seed)
            self.assertGreater(r["p_removed"], r["p"], seed)

    def test_the_class_does_not_survive_the_removal_of_the_only_aligned_stratum(self):
        w, sid = self.one_aligned_stratum(5)
        r = T.run_cell(w["fam"], w["strings"], w["strata"], stat="m1c", B=800)
        self.assertGreater(r["p_removed"], 0.05)
        self.assertNotEqual(r["class"], "CUT-ALIGNED")
        self.assertEqual(r["class_removed"], "NOT-DISTINGUISHED")

    def test_z_is_zero_when_a_stratum_cannot_move(self):
        w = world(0, n_strata=9, mode="one_string_per_stratum", nclu=(2, 2))
        r = T.run_cell(w["fam"], w["strings"], w["strata"], stat="m1c", B=50)
        self.assertTrue(all(z == 0 for z in r["z"].values()))

    def test_the_most_negative_z_is_removed_and_ties_go_to_the_smaller_stratum_id(self):
        self.assertEqual(T.pick_removed({"s2": -1.5, "s1": -1.5, "s3": 0.2}), "s1")
        self.assertEqual(T.pick_removed({"s2": -2.0, "s1": -1.5}), "s2")

    def test_the_removal_step_is_skipped_when_no_z_is_negative(self):
        self.assertIsNone(T.pick_removed({"s1": 0.0, "s2": 0.3}))
        self.assertIsNone(T.pick_removed({"s1": 0.0}))
        self.assertIsNone(T.pick_removed({}))
        self.assertEqual(T.pick_removed({"s1": 0.0, "s2": -1e-9}), "s2")

    def test_a_cell_whose_every_z_is_non_negative_skips_the_removal_and_keeps_the_class_of_all_strata(self):
        w = anti_aligned_world(10)
        for stat in ("m1c", "m2"):
            r = T.run_cell(w["fam"], w["strings"], w["strata"], stat=stat, B=300)
            self.assertTrue(all(z >= 0 for z in r["z"].values()), stat)
            self.assertIsNone(r["removed_stratum"], stat)
            self.assertIsNone(r["p_removed"], stat)
            self.assertEqual(r["class_removed"], r["class_all"], stat)
            self.assertEqual(r["class"], r["class_all"], stat)

    def test_when_every_remaining_stratum_is_dead_p_is_one_and_the_class_of_the_remainder_is_not_distinguished(self):
        # one aligned live stratum, eleven strata whose genes all carry one string (dead in both cells)
        w = merge_worlds(world(5, n_strata=1, mode="independent", sizes=[40], nclu=(3, 3)),
                         world(6, n_strata=11, mode="one_string_per_stratum", nclu=(2, 2)))
        sid, genes = T.variation_strata(w["fam"], w["strata"])[0]
        for g in genes:
            w["strings"][g] = frozenset({"own_" + w["fam"][g]})
        for stat in ("m1c", "m2"):
            r = T.run_cell(w["fam"], w["strings"], w["strata"], stat=stat, B=300)
            self.assertEqual(r["n_live"], 1, stat)
            self.assertEqual(r["removed_stratum"], sid, stat)
            self.assertEqual(r["p_removed"], 1.0, stat)
            self.assertEqual(r["class_removed"], "NOT-DISTINGUISHED", stat)
            self.assertEqual(r["class"], "UNDERPOWERED", stat)                       # one live stratum is fewer than eight

    def test_the_removal_for_the_two_cells_is_made_separately(self):
        w, sid = self.one_aligned_stratum(5)
        a = T.run_cell(w["fam"], w["strings"], w["strata"], stat="m1c", B=300)
        b = T.run_cell(w["fam"], w["strings"], w["strata"], stat="m2", B=300)
        self.assertIn("z", a)
        self.assertIn("z", b)
        self.assertNotEqual(a["z"], b["z"])

    def test_the_stratum_rank_is_the_one_based_rank_by_smallest_gene_id_among_all_strata(self):
        w, sid = self.one_aligned_stratum(5)
        r = T.run_cell(w["fam"], w["strings"], w["strata"], stat="m1c", B=200)
        ids = [s for s, _ in w["strata"]]
        self.assertEqual(r["removed_rank"], ids.index(r["removed_stratum"]) + 1)
        self.assertEqual(r["ranks"], {s: ids.index(s) + 1 for s in r["z"]})

    def test_both_cells_from_shared_draws_equal_the_single_cells(self):
        w, sid = self.one_aligned_stratum(5)
        both = T.run_cells(w["fam"], w["strings"], w["strata"], B=300)
        for stat in ("m1c", "m2"):
            one = T.run_cell(w["fam"], w["strings"], w["strata"], stat=stat, B=300)
            for k in ("obs", "p", "p_upper", "p_removed", "removed_stratum", "n_live", "n_var", "class", "class_all", "class_removed", "z"):
                self.assertEqual(both[stat][k], one[k], (stat, k))

    def test_removal_for_m2_recomputes_the_statistic_without_the_stratum(self):
        w, sid = self.one_aligned_stratum(5)
        eng = T.NullEngine(w["fam"], w["strings"], w["strata"])
        rem = eng.draws_without(sid, 12)
        genes_out = {g for s, gs in w["strata"] if s == sid for g in gs}
        keep = [g for g in w["fam"] if g not in genes_out]
        for b, fam_b in enumerate(eng.naive_draws(12)):
            self.assertEqual(rem["m2"][b], T.m2(fam_b, w["strings"], keep))
            self.assertEqual(rem["m1c"][b], T.m1c(fam_b, w["strings"], keep))


# ---------------------------------------------------------------- chunk 4: T1-e (descriptive, independent evidence)
class RegionTests(unittest.TestCase):
    def setUp(self):
        self.regions = T.Regions([("chr1", 0, 100), ("chr1", 150, 300), ("chr2", 0, 500)])

    def test_the_region_with_the_most_exonic_bases_wins(self):
        self.assertEqual(T.best_region(self.regions, "chr1", [[0, 10], [160, 200]]), ("chr1", 150, 300))
        self.assertEqual(T.best_region(self.regions, "chr1", [[50, 100], [160, 200]]), ("chr1", 0, 100))

    def test_a_tie_goes_to_the_smaller_region_id(self):
        self.assertEqual(T.best_region(self.regions, "chr1", [[90, 100], [160, 170]]), ("chr1", 0, 100))
        self.assertEqual(T.best_region(self.regions, "chr1", [[80, 100], [150, 170]]), ("chr1", 0, 100))

    def test_no_overlap_or_another_chromosome_gives_no_region(self):
        self.assertIsNone(T.best_region(self.regions, "chr3", [[0, 10]]))
        self.assertIsNone(T.best_region(self.regions, "chr1", [[100, 150]]))

    def test_a_long_region_that_starts_far_upstream_is_found(self):
        regions = T.Regions([("chr1", 0, 100000), ("chr1", 90000, 90100)])
        self.assertEqual(T.best_region(regions, "chr1", [[95000, 95100]]), ("chr1", 0, 100000))


class RankAndFitTests(unittest.TestCase):
    def test_average_ranks_for_ties_are_one_based(self):
        self.assertEqual(T.average_ranks([3, 1, 3, 2]), [3.5, 1.0, 3.5, 2.0])

    def pairs(self, seed=0, n_strata=30, per=6, effect=-2.0):
        rng = random.Random(seed)
        rows = []
        for s in range(n_strata):
            off = rng.uniform(-5, 5)            # a stratum level that fixed effects must remove
            for _ in range(per):
                x_str = rng.randint(0, 1)
                x_clu = rng.randint(0, 1)
                w = rng.uniform(0.5, 1.0)
                y = off + effect * x_str + 0.5 * x_clu + 1.0 * w + rng.gauss(0, 0.3)
                rows.append(dict(stratum=f"s{s:03d}", w=w, x_clu=x_clu, x_str=x_str, x_jac=0.5 * x_str, y=y))
        return rows

    def test_the_planted_effect_is_recovered_with_stratum_fixed_effects(self):
        fit = T.fe_ols(self.pairs(), "y", ["w", "x_clu", "x_str"])
        self.assertAlmostEqual(fit["coef"]["x_str"], -2.0, delta=0.3)
        self.assertAlmostEqual(fit["coef"]["x_clu"], 0.5, delta=0.3)

    def test_a_pure_stratum_level_has_no_effect_after_demeaning(self):
        rng = random.Random(1)
        rows = []
        for s in range(10):
            off = rng.uniform(-9, 9)
            for _ in range(6):
                rows.append(dict(stratum=f"s{s}", w=rng.random(), x_clu=rng.randint(0, 1), x_str=rng.randint(0, 1), y=off))
        fit = T.fe_ols(rows, "y", ["w", "x_clu", "x_str"])
        for k in ("w", "x_clu", "x_str"):
            self.assertAlmostEqual(fit["coef"][k], 0.0, places=9)

    def test_the_fit_equals_numpy_least_squares_on_the_stacked_demeaned_data(self):
        import numpy as np
        rows = self.pairs(seed=2, n_strata=8)
        fit = T.fe_ols(rows, "y", ["w", "x_clu", "x_str"])
        by = {}
        for r in rows:
            by.setdefault(r["stratum"], []).append(r)
        X, y = [], []
        for rs in by.values():
            m = {k: sum(r[k] for r in rs) / len(rs) for k in ("w", "x_clu", "x_str", "y")}
            for r in rs:
                X.append([r["w"] - m["w"], r["x_clu"] - m["x_clu"], r["x_str"] - m["x_str"]])
                y.append(r["y"] - m["y"])
        beta = np.linalg.lstsq(np.array(X), np.array(y), rcond=None)[0]
        for k, b in zip(("w", "x_clu", "x_str"), beta):
            self.assertAlmostEqual(fit["coef"][k], b, places=8)

    def test_collinear_regressors_do_not_raise_and_are_flagged(self):
        rows = self.pairs(seed=3)
        for r in rows:
            r["x_clu"] = r["x_str"]
        fit = T.fe_ols(rows, "y", ["w", "x_clu", "x_str"])
        self.assertTrue(fit["rank_deficient"])
        self.assertTrue(all(abs(v) < 1e6 for v in fit["coef"].values()))

    def test_the_bootstrap_resamples_strata_and_is_reproducible(self):
        rows = self.pairs(seed=4)
        a = T.bootstrap_fe(rows, "y", ["w", "x_clu", "x_str"], B=300)
        b = T.bootstrap_fe(rows, "y", ["w", "x_clu", "x_str"], B=300)
        self.assertEqual(a["interval"], b["interval"])
        lo, hi = a["interval"]["x_str"]
        self.assertLess(lo, -2.0 + 0.5)
        self.assertGreater(hi, -2.0 - 0.5)
        self.assertLess(hi, 0.0)            # the planted effect excludes 0 with 30 strata

    def test_the_bootstrap_draws_among_all_listed_strata_including_those_without_pairs(self):
        import numpy as np
        xs = ["w", "x_clu", "x_str"]
        rows = self.pairs(seed=6, n_strata=2, per=8)
        ids = ["s000", "s001", "empty"]
        out = T.bootstrap_fe(rows, "y", xs, B=1, seed=5, strata=ids)
        blocks = T._stratum_blocks(rows, "y", xs, strata=ids)
        self.assertEqual([b[0] for b in blocks], ids)
        self.assertEqual(blocks[2][1].shape, (0, 3))
        rng = random.Random(5)
        idx = [rng.randrange(3) for _ in range(3)]                         # 3 draws among 3 strata, not 2 among 2
        a = sum(blocks[i][1].T @ blocks[i][1] for i in idx)
        b = sum(blocks[i][1].T @ blocks[i][2] for i in idx)
        beta = np.linalg.lstsq(a, b, rcond=None)[0]
        for name, v in zip(xs, beta):
            self.assertAlmostEqual(out["replicates"][name][0], float(v), places=9)
        plain = T.bootstrap_fe(rows, "y", xs, B=1, seed=5)
        self.assertEqual(len(T._stratum_blocks(rows, "y", xs)), 2)
        self.assertEqual(set(plain["point"]), set(xs))

    def test_the_point_estimate_does_not_depend_on_empty_strata(self):
        xs = ["w", "x_clu", "x_str"]
        rows = self.pairs(seed=7, n_strata=4, per=6)
        a = T.bootstrap_fe(rows, "y", xs, B=5, seed=1)["point"]
        b = T.bootstrap_fe(rows, "y", xs, B=5, seed=1, strata=["a_empty", "s000", "s001", "s002", "s003"])["point"]
        for k in xs:
            self.assertAlmostEqual(a[k], b[k], places=10)

    def test_the_percentile_interval_uses_the_nearest_rank_of_the_sorted_replicates(self):
        self.assertEqual(T.percentile_interval(list(range(1, 2001)), 2000), (50, 1950))

    def test_an_outcome_unrelated_to_the_regressors_gives_an_interval_around_zero(self):
        rng = random.Random(8)
        rows = [dict(stratum=f"s{s}", w=rng.random(), x_clu=rng.randint(0, 1), x_str=rng.randint(0, 1), y=rng.gauss(0, 1))
                for s in range(30) for _ in range(6)]
        lo, hi = T.bootstrap_fe(rows, "y", ["w", "x_clu", "x_str"], B=400)["interval"]["x_str"]
        self.assertLess(lo, 0.0)
        self.assertGreater(hi, 0.0)


class PairBuilderTests(unittest.TestCase):
    def test_pairs_join_genes_of_one_stratum_whose_nodes_share_an_edge(self):
        fam = {"g1": "f", "g2": "f", "g3": "h", "g4": "h"}
        u1 = {"g1": fs("a"), "g2": fs("a"), "g3": fs("b"), "g4": fs("a", "b")}
        co = {"g1": 4.0, "g2": 5.0, "g3": 9.0, "g4": None}
        node = {"g1": "n1", "g2": "n2", "g3": "n3", "g4": "n4"}
        edges = {("n1", "n2"): 0.9, ("n2", "n3"): 0.8, ("n3", "n4"): 0.7}
        region = {"g1": "R1", "g2": "R1", "g3": "R2", "g4": "R2"}
        rows = T.t1e_pairs([("g1", ["g1", "g2", "g3", "g4"])], fam, u1, co, node, edges, region)
        by = {(r["g"], r["h"]): r for r in rows}
        self.assertEqual(sorted(by), [("g1", "g2"), ("g2", "g3")])          # g4 has no copy number; g1-g3 share no edge
        r = by[("g1", "g2")]
        self.assertEqual((r["w"], r["x_clu"], r["x_str"], r["y_loc"], r["y_cn_raw"]), (0.9, 1, 1, 1, 1.0))
        r = by[("g2", "g3")]
        self.assertEqual((r["x_clu"], r["x_str"], r["y_loc"], r["x_jac"]), (0, 0, 0, 0.0))

    def test_two_genes_on_one_node_are_not_joined_and_an_edge_may_be_listed_either_way(self):
        fam = {"g1": "f", "g2": "f", "g3": "f"}
        u1 = {g: fs("a") for g in fam}
        co = {"g1": 1.0, "g2": 2.0, "g3": 3.0}
        node = {"g1": "n1", "g2": "n1", "g3": "n2"}
        rows = T.t1e_pairs([("g1", ["g1", "g2", "g3"])], fam, u1, co, node, {("n2", "n1"): 0.5}, {})
        self.assertEqual(sorted((r["g"], r["h"]) for r in rows), [("g1", "g3"), ("g2", "g3")])

    def test_the_outcome_is_ranked_over_all_pairs_of_all_strata(self):
        fam = {"g1": "f", "g2": "f", "g3": "h", "g4": "h"}
        u1 = {g: fs("a") for g in fam}
        co = {"g1": 1.0, "g2": 4.0, "g3": 1.0, "g4": 2.0}
        node = {g: "n" + g[1] for g in fam}
        rows = T.t1e_pairs([("g1", ["g1", "g2"]), ("g3", ["g3", "g4"])], fam, u1, co, node, {("n1", "n2"): 0.5, ("n3", "n4"): 0.6}, {})
        T.add_cn_rank(rows)
        self.assertEqual([r["y_cn"] for r in rows], [2.0, 1.0])


# ---------------------------------------------------------------- chunk 5: gate helpers (partition deficiency, ties, sized worlds, arms)
class PartitionDeficiencyTests(unittest.TestCase):
    def test_blocks_that_nest_in_the_cells_cost_nothing(self):
        blocks = {"S1": ["a", "b"], "S2": ["c", "d", "e"]}
        cell = {"a": "Q1", "b": "Q1", "c": "Q2", "d": "Q2", "e": "Q2"}
        r = T.d_partition(blocks, cell.get)
        self.assertEqual((r["moved"], r["split"], r["blocks"]), (0, 0, 2))

    def test_a_split_block_costs_the_genes_outside_its_largest_cell(self):
        blocks = {"S": ["a", "b", "c", "d"]}
        cell = {"a": "Q1", "b": "Q1", "c": "Q1", "d": "Q2"}
        r = T.d_partition(blocks, cell.get)
        self.assertEqual((r["moved"], r["split"]), (1, 1))

    def test_ungrouped_genes_are_singleton_cells_unless_they_are_ignored(self):
        blocks = {"S": ["a", "b", "c"]}
        cell = {"a": "Q1", "b": "Q1", "c": None}
        sing = T.d_partition(blocks, cell.get)
        self.assertEqual((sing["moved"], sing["split"], sing["blocks"]), (1, 1, 1))
        ign = T.d_partition(blocks, cell.get, ignore_ungrouped=True)
        self.assertEqual((ign["moved"], ign["split"], ign["blocks"]), (0, 0, 1))

    def test_two_ungrouped_genes_are_two_different_singleton_cells(self):
        r = T.d_partition({"S": ["a", "b"]}, {"a": None, "b": None}.get)
        self.assertEqual((r["moved"], r["split"]), (1, 1))

    def test_blocks_with_fewer_than_two_genes_are_not_scored_and_the_test_follows_the_ignore_rule(self):
        blocks = {"S1": ["a"], "S2": ["b", "c"]}
        cell = {"a": "Q1", "b": "Q2", "c": None}
        self.assertEqual(T.d_partition(blocks, cell.get)["blocks"], 1)
        self.assertEqual(T.d_partition(blocks, cell.get, ignore_ungrouped=True)["blocks"], 0)


class ExonLabelTests(unittest.TestCase):
    def test_an_exon_with_two_equal_top_labels_is_a_tie_and_the_lex_and_numeric_picks_are_reported(self):
        c = colors([("chr1", 0, 100, "SD10"), ("chr1", 50, 150, "SD9")])
        rows = T.exon_label_rows(c, "chr1", [[50, 100]])
        self.assertEqual(len(rows), 1)
        r = rows[0]
        self.assertEqual((r["labelled"], r["tie"], r["lex"], r["num"]), (True, True, "SD10", "SD9"))

    def test_no_tie_when_the_top_label_is_alone_and_no_label_when_nothing_overlaps(self):
        c = colors([("chr1", 0, 100, "SD1"), ("chr1", 50, 150, "SD2")])
        rows = T.exon_label_rows(c, "chr1", [[0, 100], [500, 600]])
        self.assertEqual((rows[0]["labelled"], rows[0]["tie"], rows[0]["lex"]), (True, False, "SD1"))
        self.assertEqual((rows[1]["labelled"], rows[1]["tie"], rows[1]["lex"]), (False, False, None))

    def test_a_single_overlapping_label_is_never_a_tie(self):
        c = colors([("chr1", 0, 100, "SD1")])
        self.assertFalse(T.exon_label_rows(c, "chr1", [[0, 50]])[0]["tie"])


class SizedWorldTests(unittest.TestCase):
    SIZES = [[5, 3, 2], [4, 4], [6, 2, 2, 1], [3, 3, 3]]

    def test_cluster_sizes_per_variation_bearing_stratum_are_read_back_in_descending_order(self):
        w = T.world_from_sizes(self.SIZES, "independent", seed=0)
        self.assertEqual(T.cluster_sizes(w["fam"], w["strata"]), [[5, 3, 2], [4, 4], [6, 2, 2, 1], [3, 3, 3]])

    def test_identity_gives_every_cluster_its_own_string(self):
        w = T.world_from_sizes(self.SIZES, "identity", seed=0)
        for sid, genes in w["strata"]:
            by = {}
            for g in genes:
                by.setdefault(w["fam"][g], set()).add(w["strings"][g])
            self.assertTrue(all(len(v) == 1 for v in by.values()))
            self.assertEqual(len({next(iter(v)) for v in by.values()}), len(by))

    def test_union_of_adjacent_pairs_acts_only_in_strata_with_at_least_three_clusters(self):
        w = T.world_from_sizes(self.SIZES, "pairs", seed=0)
        n_strings = [len({w["strings"][g] for g in genes}) for _, genes in w["strata"]]
        self.assertEqual(n_strings, [2, 2, 2, 2])          # 3 clusters -> 2 strings, 2 clusters -> identity (2), 4 -> 2, 3 -> 2

    def test_the_split_gives_two_alternating_types_to_every_cluster_of_at_least_three_genes(self):
        w = T.world_from_sizes(self.SIZES, "split", seed=0)
        for sid, genes in w["strata"]:
            by = {}
            for g in genes:
                by.setdefault(w["fam"][g], []).append(w["strings"][g])
            for lab, strs in by.items():
                self.assertEqual(len(set(strs)), 2 if len(strs) >= 3 else 1, lab)
                if len(strs) >= 3:
                    self.assertEqual(strs[0], strs[2])
                    self.assertNotEqual(strs[0], strs[1])

    def test_only_the_largest_stratum_follows_its_strings_in_the_planted_set(self):
        w = T.world_from_sizes(self.SIZES, "largest_follows", seed=3)

        def follows(genes):
            by = {}
            for g in genes:
                by.setdefault(w["fam"][g], set()).add(w["strings"][g])
            return all(len(v) == 1 for v in by.values()) and len({next(iter(v)) for v in by.values()}) == len(by)

        sizes = {sid: len(genes) for sid, genes in w["strata"]}
        largest = min((sid for sid in sizes if sizes[sid] == max(sizes.values())))
        self.assertEqual(sizes[largest], 11)
        self.assertTrue(follows(dict(w["strata"])[largest]))
        self.assertFalse(all(follows(genes) for sid, genes in w["strata"] if sid != largest))

    def test_gene_ids_sort_in_stratum_order_and_every_world_is_reproducible(self):
        a = T.world_from_sizes(self.SIZES, "independent", seed=4)
        b = T.world_from_sizes(self.SIZES, "independent", seed=4)
        self.assertEqual(a["strings"], b["strings"])
        ids = [sid for sid, _ in a["strata"]]
        self.assertEqual(ids, sorted(ids))


class SameLocusTests(unittest.TestCase):
    def info(self, **kw):
        base = dict(chrom="chr1", strand="+", start=0, end=100, exons=[[0, 100]])
        base.update(kw)
        return base

    def test_the_gene_with_the_smaller_exon_union_of_an_overlapping_same_strand_pair_is_deleted(self):
        info = {"g1": self.info(start=0, end=100, exons=[[0, 100]]), "g2": self.info(start=50, end=200, exons=[[50, 80]])}
        pairs, losers = T.same_locus_pairs(info, ["g1", "g2"])
        self.assertEqual(pairs, [("g1", "g2")])
        self.assertEqual(losers, {"g2"})

    def test_equal_exon_unions_delete_the_larger_gene_id(self):
        info = {"g1": self.info(exons=[[0, 50]]), "g2": self.info(start=10, end=60, exons=[[10, 60]])}
        _, losers = T.same_locus_pairs(info, ["g1", "g2"])
        self.assertEqual(losers, {"g2"})

    def test_other_strands_chromosomes_and_touching_spans_do_not_pair(self):
        info = {"a": self.info(), "b": self.info(strand="-"), "c": self.info(chrom="chr2"), "d": self.info(start=100, end=200, exons=[[100, 200]])}
        pairs, losers = T.same_locus_pairs(info, ["a", "b", "c", "d"])
        self.assertEqual((pairs, losers), ([], set()))

    def test_a_gene_in_two_pairs_is_deleted_once_and_the_pairs_are_independent(self):
        info = {"a": self.info(start=0, end=100, exons=[[0, 10]]), "b": self.info(start=50, end=150, exons=[[50, 100]]),
                "c": self.info(start=120, end=300, exons=[[120, 300]])}
        pairs, losers = T.same_locus_pairs(info, ["a", "b", "c"])
        self.assertEqual(pairs, [("a", "b"), ("b", "c")])
        self.assertEqual(losers, {"a", "b"})

    def test_only_the_listed_genes_are_paired(self):
        info = {"a": self.info(), "b": self.info(start=10, end=50)}
        self.assertEqual(T.same_locus_pairs(info, ["a"]), ([], set()))


class UnitCoverPortTests(unittest.TestCase):
    """The port of the arms of soto_m2_unit_cover_cn.py on a toy table whose answer is worked out by hand."""

    def genes(self):
        def g(i, S, comp, dom, cs=None, co=None, half="dev"):
            return dict(id=i, S=set(S), comp=Counter(comp), dom=set(dom), cn={"s1c": cs, "ours": co}, half=half, c="chr1",
                        s=int(i[1:]) * 1000, e=int(i[1:]) * 1000 + 500)
        return [g("g1", ["ID_1"], {"D1": 100}, ["D1"]),
                g("g2", ["ID_1"], {"D1": 80, "D9": 5}, ["D1"]),
                g("g3", ["ID_2"], {"D2": 90}, ["D2"]),
                g("g4", ["ID_1", "ID_2"], {"D1": 50, "D2": 60}, ["D1", "D2"], half="dev"),
                g("g5", ["ID_1", "ID_3"], {"D9": 70}, ["D9"], half="heldout")]

    def test_arm_a_predicts_the_owner_of_every_duplicon_a_multi_gene_overlaps(self):
        r = T.unit_cover(self.genes(), "A", "s1c")
        self.assertEqual(r["pred"]["g4"], {"ID_1", "ID_2"})
        self.assertEqual(r["pred"]["g5"], {"ID_1"})          # D9 is owned by ID_1: g2 votes with 5 bp, nobody else overlaps it
        self.assertEqual(r["n_units"], 3)                      # D1, D2, D9

    def test_arm_b_uses_only_dominant_duplicons_so_the_sliver_unit_disappears(self):
        r = T.unit_cover(self.genes(), "B", "s1c")
        self.assertEqual(r["n_units"], 2)
        self.assertEqual(r["pred"]["g5"], set())

    def test_exact_recall_and_precision_are_means_over_the_multi_family_genes(self):
        r = T.unit_cover(self.genes(), "A", "s1c")
        self.assertAlmostEqual(r["exact"], 0.5)                # g4 exact, g5 ({ID_1,ID_3} against {ID_1}) not
        self.assertAlmostEqual(r["recall"], (1.0 + 0.5) / 2)
        self.assertAlmostEqual(r["precision"], (1.0 + 1.0) / 2)

    def test_the_halves_are_scored_separately_by_the_mean_jaccard_of_the_multi_family_genes(self):
        r = T.unit_cover(self.genes(), "A", "s1c")
        self.assertEqual(r["obs"]["dev"], (1.0, 1))
        self.assertEqual(r["obs"]["heldout"], (0.5, 1))

    def test_the_clean_gene_check_leaves_the_gene_out_of_its_own_vote(self):
        r = T.unit_cover(self.genes(), "A", "s1c")
        # g1 and g2 vote for each other's duplicon D1 (own family); g3's only duplicon D2 has no other voter -> empty
        self.assertAlmostEqual(r["clean_ok"], 2 / 3)
        self.assertAlmostEqual(r["clean_empty"], 1 / 3)
        self.assertAlmostEqual(r["clean_multi"], 0.0)

    def test_copy_number_classes_split_a_duplicon_into_units_by_single_linkage(self):
        gs = self.genes()
        gs[0]["cn"]["s1c"], gs[1]["cn"]["s1c"], gs[2]["cn"]["s1c"] = 2.0, 9.0, 3.0
        r = T.unit_cover(gs, "C", "s1c")
        self.assertEqual(r["n_units"], 4)                      # D1 -> two classes (2.0 | 9.0), D2 -> one, D9 -> one (g2 has a copy number)

    def test_a_gene_with_a_copy_number_uses_only_the_classes_within_less_than_two_of_it(self):
        gs = self.genes()
        gs[0]["cn"]["s1c"], gs[1]["cn"]["s1c"] = 2.0, 9.0
        gs[1]["S"] = {"ID_5"}                                  # the two classes of D1 now have different owners
        gs[3]["cn"]["s1c"] = 8.0                               # g4 lies within < 2 of the class {9.0} only
        gs[3]["comp"], gs[3]["dom"] = Counter({"D1": 50}), {"D1"}
        r = T.unit_cover(gs, "C", "s1c")
        self.assertEqual(r["pred"]["g4"], {"ID_5"})
        gs[3]["cn"]["s1c"] = None                              # without a copy number the gene uses every class of the duplicon
        self.assertEqual(T.unit_cover(gs, "C", "s1c")["pred"]["g4"], {"ID_1", "ID_5"})

    def test_an_empty_half_scores_zero_instead_of_failing(self):
        gs = self.genes()
        gs[3]["half"] = "heldout"
        r = T.unit_cover(gs, "A", "s1c")
        self.assertEqual(r["obs"]["dev"], (0.0, 0))


# ---------------------------------------------------------------- chunk 6: strata that families may span, and the U2 class index
def spanning_world(seed, n=40, nstrata=5, nfam=7, vocab="ABC"):
    """Genes in strata that families are free to span (the shape of the chromosome, SD98-region and F_Q-group nulls)."""
    rng = random.Random(seed)
    genes = [f"g{i:03d}" for i in range(n)]
    fam = {g: f"F{rng.randrange(nfam)}" for g in genes}
    strings = {g: frozenset(rng.sample(vocab, rng.randint(1, 3))) for g in genes}
    cut = sorted(rng.sample(range(1, n), nstrata - 1))
    strata, lo = [], 0
    for hi in cut + [n]:
        chunk = genes[lo:hi]
        strata.append((chunk[0], chunk))
        lo = hi
    return fam, strings, strata


class GeneralNullTests(unittest.TestCase):
    def test_the_observed_statistics_equal_the_direct_ones(self):
        for seed in range(5):
            fam, st, strata = spanning_world(seed)
            eng = T.GeneralNull(fam, st, strata)
            obs = eng.observed()
            self.assertEqual((obs["m1c"], obs["m2"]), (T.m1c(fam, st), T.m2(fam, st)))

    def test_every_draw_equals_a_direct_computation_on_the_permuted_labels(self):
        for seed in range(4):
            fam, st, strata = spanning_world(seed)
            eng = T.GeneralNull(fam, st, strata)
            res = eng.draws(25)
            for b, fam_b in enumerate(eng.naive_draws(25)):
                self.assertEqual(res["tot_m1c"][b], T.m1c(fam_b, st), (seed, b))
                self.assertEqual(res["tot_m2"][b], T.m2(fam_b, st), (seed, b))

    def test_it_agrees_draw_for_draw_with_the_confined_engine_when_families_do_not_span_strata(self):
        for seed in range(4):
            w = world(seed)
            a = T.NullEngine(w["fam"], w["strings"], w["strata"]).draws(30)
            b = T.GeneralNull(w["fam"], w["strings"], w["strata"]).draws(30)
            self.assertEqual(a["tot_m1c"], b["tot_m1c"], seed)
            self.assertEqual(a["tot_m2"], b["tot_m2"], seed)

    def test_the_label_multiset_of_every_stratum_is_kept(self):
        fam, st, strata = spanning_world(2)
        eng = T.GeneralNull(fam, st, strata)
        for fam_b in eng.naive_draws(10):
            for sid, genes in strata:
                self.assertEqual(Counter(fam_b[g] for g in genes), Counter(fam[g] for g in genes))

    def test_blocks_of_draws_continue_one_random_stream(self):
        fam, st, strata = spanning_world(3)
        eng = T.GeneralNull(fam, st, strata)
        whole = eng.draws(30)
        a, b = eng.draws(30, start=0, end=11), eng.draws(30, start=11, end=30)
        self.assertEqual(whole["tot_m1c"], a["tot_m1c"] + b["tot_m1c"])
        self.assertEqual(whole["tot_m2"], a["tot_m2"] + b["tot_m2"])

    def test_a_gene_in_two_strata_is_refused(self):
        with self.assertRaises(ValueError):
            T.GeneralNull({"a": "F", "b": "F"}, {"a": fs("x"), "b": fs("x")}, [("a", ["a", "b"]), ("b", ["b"])])

    def test_the_cell_reports_p_and_the_tails_without_a_class(self):
        fam, st, strata = spanning_world(1)
        r = T.run_general_cells(fam, st, strata, B=200)
        for stat in ("m1c", "m2"):
            c = r[stat]
            self.assertNotIn("class", c)
            self.assertGreater(c["p"], 0.0)
            self.assertLessEqual(c["p"], 1.0)
            self.assertLessEqual(c["p_upper"], 1.0)
            self.assertEqual(c["B"], 200)

    def test_strata_that_cannot_move_leave_the_statistic_fixed(self):
        genes = [f"g{i}" for i in range(6)]
        fam = {g: "F" for g in genes}
        st = {g: fs("a") for g in genes}
        r = T.run_general_cells(fam, st, [("g0", genes)], B=20)
        self.assertEqual((r["m1c"]["draw_min"], r["m1c"]["draw_max"]), (0, 0))


class U2IndexTests(unittest.TestCase):
    def test_the_index_maps_each_label_and_copy_number_to_its_class(self):
        u1 = {"g1": fs("a"), "g2": fs("a"), "g3": fs("a", "b")}
        cn = {"g1": 2.0, "g2": 9.0, "g3": 3.0}
        idx = T.u2_index(u1, cn)
        self.assertEqual(idx["a"], {2.0: 0, 3.0: 0, 9.0: 1})
        self.assertEqual(idx["b"], {3.0: 0})

    def test_u2_strings_are_built_from_the_index(self):
        u1 = {"g1": fs("a"), "g2": fs("a"), "g3": fs("a", "b")}
        cn = {"g1": 2.0, "g2": 9.0, "g3": 3.0}
        self.assertEqual(T.u2_strings(u1, cn), {"g1": fs(("a", 0)), "g2": fs(("a", 1)), "g3": fs(("a", 0), ("b", 0))})

    def test_the_ordered_u2_string_carries_the_class_of_each_label(self):
        u1 = {"g1": fs("a"), "g3": fs("a", "b")}
        cn = {"g1": 2.0, "g3": 9.0}
        idx = T.u2_index(u1, cn)
        out = T.u2_ordered({"g1": ("a",), "g3": ("a", "b")}, cn, idx)
        self.assertEqual(out, {"g1": (("a", 0),), "g3": (("a", 1), ("b", 0))})

    def test_a_gene_without_a_copy_number_has_no_u2_string(self):
        u1 = {"g1": fs("a"), "g2": fs("a")}
        self.assertEqual(set(T.u2_ordered({"g1": ("a",), "g2": ("a",)}, {"g1": 1.0, "g2": None}, T.u2_index(u1, {"g1": 1.0, "g2": None}))), {"g1"})


if __name__ == "__main__":
    unittest.main()
