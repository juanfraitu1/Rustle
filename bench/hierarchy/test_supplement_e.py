#!/usr/bin/env python3
"""test_supplement_e.py -- supplementary tests for the Phase E code (expansion_table.py, liftoff_cross.py).

Written by the mutation review ('mut'): every test here PASSES on the frozen code and was written to kill a mutant that the builder's 171 tests let survive
(see the review report).  Hermetic: made-up tables, temp directories (tempfile.mkdtemp under $TMPDIR), stub scripts; no Liftoff, no minimap2, no real input and
nothing under the real run directory is read or written (every class derives from T.HermeticCase, which replaces E_ROOT, the self-lift root, the Liftoff binary, the gorilla
FASTA and GFF, the minimap2 of the Liftoff environment and the outgroups' files by paths that cannot exist, before any setUp runs; run the module under harness/guard.py to prove it).
Fix round of 2026-10-07: the tests that pinned a behaviour the builder changed on purpose were updated (listed in fix_round/CHANGES.md).  Fixture builders (rec, row, lrow,
outgroup_rows, make_indexes, WorldCase, make_run_dir, ...) are imported from the builder's test module; none of its tests is re-run here.

Run: python3 -B -m unittest bench/hierarchy/test_supplement_e.py     (from the repository root)
"""
import collections
import contextlib
import csv
import hashlib
import inspect
import io
import json
import os
import random
import re
import signal
import stat
import subprocess
import sys
import tempfile
import time
import types
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import expansion_table as X  # noqa: E402
import liftoff_cross as LC  # noqa: E402
import test_expansion_table as T  # noqa: E402  (fixture builders only)

rec, row, lrow = T.rec, T.row, T.lrow
REAL = T.REAL        # the registered paths as the modules define them: every test of this module runs with them replaced by unreachable ones (T.HermeticCase)


def rd(path):
    with open(path) as fh:
        return fh.read()


def write_exec(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)
    return path


def alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    # a zombie still answers kill(0); read its state
    try:
        with open(f"/proc/{pid}/stat") as fh:
            return fh.read().rsplit(")", 1)[1].split()[0] != "Z"
    except OSError:
        return False



# ===================================================================================================== A. registered parameters
class RegisteredParameterTests(T.HermeticCase):
    """the registered values themselves (Amendment 2 sections 3, 6-8; prereg section 7): nothing in the builder's tests names them"""

    def test_thresholds_of_the_prereg(self):
        self.assertEqual((X.SHORT_BP, X.HUB_MAX, X.MATCH_FRAC, X.MIN_COVERAGE, X.MIN_SEQ_ID, X.SEED_T2, X.RANK_CAP), (200, 51, 0.5, 0.5, 0.5, 20261007, 2500))
        self.assertEqual(X.COPY_CLASSES, ("extra_copy", "dropped_M1", "dropped_M2"))
        self.assertEqual(X.PLACED_CLASSES, ("in_place", "moved"))
        self.assertEqual(X.PRIMARY_CLASSES, ("in_place", "moved", "partial", "unmapped"))
        self.assertEqual(X.REASON_ORDER, ("outgroup_copy_row", "no_placement", "no_outgroup_locus"))
        self.assertEqual(X.CONTROL_REASONS, ("outgroup_copy_row", "different_triple", "no_placement"))
        self.assertEqual(X.NOT_ELIGIBLE, ("no_placement",))
        self.assertEqual(X.OUTGROUPS, (("H", "human"), ("C", "chimpanzee"), ("O", "orangutan")))

    def test_the_lrpap1_control_specification(self):
        # fix round, addendum 2 (D-m2): + expected_triple, the gene-LRPAP1 anchor of each outgroup (in_place, single-copy in all three self-lifts)
        self.assertEqual(X.LRPAP1_CONTROL, dict(reference="gene-LRPAP1", not_required=["gene-LOC129523503", "gene-LOC129526389"],
                                                fragments=["gene-LOC134756368", "gene-LOC115932954", "gene-LOC115932756"], min_share=0.5,
                                                expected_triple=("gene-LRPAP1", "gene-LRPAP1", "gene-LRPAP1")))

    def test_the_driver_parameters_of_amendment_2(self):
        self.assertEqual((LC.CALL_BUDGET_S, LC.CALL_LIMIT_S, LC.HARD_LIMIT_S, LC.RSS_LIMIT_GB, LC.PER_RECORD_S, LC.STARTUP_S), (530, 555, 590, 16.0, 0.09, 25.0))
        self.assertEqual((LC.TERM_GRACE_S, LC.KILL_WAIT_S, LC.INSESSION_MARGIN_S, LC.POLL_S, LC.WALL_FACTOR), (10.0, 10.0, 5.0, 0.5, 2.0))
        self.assertEqual((LC.INPUT_BIG_BYTES, LC.SHALLOW_SHA_BYTES, LC.N_PLANTED), (10 ** 9, 100 * 10 ** 6, 100))
        self.assertEqual(LC.predict_wall(1000), 25.0 + 1000 * 0.09 * 1.25)
        self.assertEqual(LC.predict_wall(0), 25.0)

    def test_function_defaults_that_main_relies_on(self):
        sig = lambda f: {k: v.default for k, v in inspect.signature(f).parameters.items()}  # noqa: E731
        d = sig(LC.run_budgeted)
        self.assertEqual((d["budget_s"], d["hard_s"], d["per_record_s"], d["startup_s"], d["rss_limit_gb"], d["call_limit_s"], d["ensure_links"]),
                         (530, 590, 0.09, 25.0, 16.0, 555, True))
        d = sig(LC.run_shard)
        self.assertEqual((d["hard_s"], d["rss_limit_gb"], d["poll_s"]), (590, 16.0, 0.5))
        d = sig(LC.timed_watch)
        self.assertEqual((d["rss_limit_gb"], d["poll_s"]), (16.0, 0.5))
        self.assertEqual(sig(LC.check_rss)["limit_gb"], 16.0)
        self.assertEqual(sig(LC.plan_calls)["budget_s"], 530)
        self.assertEqual(sig(LC.predict_wall)["margin"], 1.25)
        self.assertEqual(sig(LC.liftoff_command)["threads"], 4)
        self.assertEqual(sig(X.plan_run_shards)["cap"], 2500)
        self.assertEqual(sig(X.rank_plan)["cap"], 2500)
        self.assertEqual((sig(X.plant_negative)["n_records"], sig(X.plant_negative)["seed"]), (100, 20261007))
        self.assertEqual(sig(X.check_negative)["bar"], 1)
        self.assertEqual((sig(LC.build_all)["n_neg"], sig(LC.build_all)["hub_mode"], sig(LC.build_all)["strict"]), (100, "nodes", False))
        self.assertEqual(sig(X.called_components)["hub_mode"], "nodes")
        self.assertEqual(sig(X.write_outputs)["hub_mode"], "nodes")
        self.assertEqual(sig(X.Component.size)["mode"], "nodes")
        self.assertEqual(sig(X.Component.hub)["mode"], "nodes")

    def test_the_cross_species_inputs_are_those_of_amendment_2(self):
        h = "/mnt/linuxdisk/home/juanfraitu"
        og = REAL.OUTGROUPS
        self.assertEqual((og["H"]["species"], og["H"]["genome"], og["H"]["mmi"]), ("human", h + "/winloci_data/chm13v2.0.fa", h + "/winloci_data/chm13v2.0.fa.mmi"))
        self.assertEqual((og["C"]["species"], og["C"]["genome"], og["C"]["mmi"]),
                         ("chimpanzee", h + "/winloci_data/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna", h + "/ggo_npip/ref/PTR.fa.mmi"))
        self.assertEqual((og["O"]["species"], og["O"]["genome"], og["O"]["mmi"]),
                         ("orangutan", h + "/winloci_data/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna", h + "/o1_falsemerge/lit/evid/PPY/genome/PPY.fa.mmi"))
        for k, sp in (("H", "human"), ("C", "chimpanzee"), ("O", "orangutan")):
            self.assertEqual(og[k]["selflift"], f"/mnt/linuxdisk/tmp/rustle_figures/liftoff/{sp}/liftoff_loci.tsv")
        self.assertEqual(REAL.E_ROOT, "/mnt/linuxdisk/tmp/hier_run/e")
        self.assertEqual(REAL.GORILLA_FASTA, h + "/_from_wsl/winloci_scratch/GGO.fasta")
        self.assertEqual(REAL.GORILLA_GFF, h + "/winloci_data/GGO_genomic.gff")
        self.assertEqual(REAL.LIFTOFF_BIN, "/home/juanfra/miniforge3/envs/liftoff/bin/liftoff")
        self.assertEqual(REAL.MINIMAP2_ENV, "/home/juanfra/miniforge3/envs/liftoff/bin/minimap2")
        self.assertEqual(REAL.LIFT_ROOT, "/mnt/linuxdisk/tmp/rustle_figures/liftoff")

    def test_the_files_whose_identity_is_frozen(self):
        # fix round 2: + bench/hierarchy/test_supplement_e2.py (the second mutation review's tests): 8 files
        self.assertEqual(LC.CODE_FILES, ("bench/hierarchy/expansion_table.py", "bench/hierarchy/liftoff_cross.py", "bench/hierarchy/test_expansion_table.py",
                                         "bench/hierarchy/test_supplement_e.py", "bench/hierarchy/test_supplement_e2.py", "figures/_liftoff.py", "figures/figlib.py", "tools/rlock.sh"))

    def test_the_registered_counts_of_the_input_report(self):
        e = X.EXPECTED
        self.assertEqual((e["components_size_ge2"], e["annotated_records_in_components"], e["annotated_ge200bp"], e["components_ge2_such"]), (609, 2976, 2249, 307))
        self.assertEqual((e["hubs"], e["called"], e["records"], e["J"], e["runs_with_at_most_10"]), (7, 300, 1403, 44, 26))
        self.assertEqual(e["by_rtype"], {"gene": 907, "pseudogene": 496})
        self.assertEqual(e["records_per_run_first5"], [300, 300, 165, 116, 81])
        self.assertEqual(e["supported"], {"gorilla_KB3781": 80, "gorilla_OR6737": 139, "either": 157})
        self.assertEqual(e["supported_hubs"], {"gorilla_KB3781": 4, "gorilla_OR6737": 5})
        r = X.EXPECTED_RECORDS_READING
        self.assertEqual((r["hubs"], r["called"], r["records"], r["by_rtype"], r["records_per_run_first5"]),
                         (5, 302, 1446, {"gene": 944, "pseudogene": 502}, [302, 302, 167, 118, 83]))


# ===================================================================================================== B. interval helpers
class IntervalTests(T.HermeticCase):
    def bases(self, iv):
        return {b for s, e in iv for b in range(s, e)}

    def test_merge_iv_contained_adjacent_and_unsorted_blocks(self):
        self.assertEqual(X.merge_iv([(0, 100), (10, 20)]), [(0, 100)])
        self.assertEqual(X.merge_iv([(0, 10), (10, 20)]), [(0, 20)])
        self.assertEqual(X.merge_iv([(20, 30), (0, 10), (5, 12)]), [(0, 12), (20, 30)])
        self.assertEqual(X.merge_iv([(5, 50), (0, 10), (7, 20)]), [(0, 50)])
        self.assertEqual(X.merge_iv([]), [])

    def test_merge_iv_agrees_with_a_base_count(self):
        rng = random.Random(11)
        for _ in range(300):
            iv = []
            for _ in range(rng.randrange(0, 7)):
                s = rng.randrange(0, 50)
                iv.append((s, s + rng.randrange(1, 15)))
            m = X.merge_iv(iv)
            self.assertEqual(self.bases(m), self.bases(iv))
            self.assertTrue(all(a[1] < b[0] for a, b in zip(m, m[1:])), (iv, m))
            self.assertEqual(X.iv_len(m), len(self.bases(iv)))

    def test_iv_inter_agrees_with_a_base_count(self):
        rng = random.Random(5)

        def rand_iv():
            pts = sorted(rng.sample(range(0, 80), rng.randrange(0, 6) * 2))
            return [(pts[i], pts[i + 1]) for i in range(0, len(pts), 2)]
        for _ in range(500):
            a, b = rand_iv(), rand_iv()
            self.assertEqual(X.iv_inter(a, b), len(self.bases(a) & self.bases(b)), (a, b))

    def test_iv_inter_named_cases(self):
        self.assertEqual(X.iv_inter([(0, 10), (20, 30)], [(5, 25)]), 10)
        self.assertEqual(X.iv_inter([(0, 100)], [(10, 20), (30, 40), (90, 150)]), 30)
        self.assertEqual(X.iv_inter([(10, 20), (30, 40)], [(0, 100)]), 20)
        self.assertEqual(X.iv_inter([(0, 10)], [(10, 20)]), 0)
        self.assertEqual(X.iv_inter([], [(0, 10)]), 0)

    def test_parse_iv(self):
        self.assertEqual(X.parse_iv("5-10,20-30"), [(5, 10), (20, 30)])
        for s in ("", "-", "NA"):
            self.assertEqual(X.parse_iv(s), [])
        self.assertEqual(X.iv_str([(5, 10), (20, 30)]), "5-10,20-30")


# ===================================================================================================== C. pool, such records, called set
class PoolRuleTests(T.HermeticCase):
    def pool_with(self, specs):
        """specs: (id, cls, contig, start0, end, rtype); every record carries one dropped_M1 copy row on the first record, so all join one component"""
        recs, rows = {}, []
        first = specs[0][0]
        for rid, cls, contig, s0, e, rtype in specs:
            recs[rid] = rec(rid, contig, s0, e, rtype=rtype)
            rows.append(row(rid, cls=cls, contig=contig, start0=s0, end=e, rtype=rtype))
            if rid != first:
                rows.append(row(rid, cls="dropped_M1", contig="cx", start0=0, end=100, note=f"overlaps {first} (c1)"))
        for r in rows:
            if r["cls"] in X.COPY_CLASSES:
                r["tag"] = "1"
        return X.build_pool(rows, recs)

    def test_such_records_include_moved_and_pseudogene_records_and_the_200_bp_edge(self):
        pool = self.pool_with([("g-in", "in_place", "c1", 0, 2000, "gene"), ("g-moved", "moved", "c2", 0, 2000, "gene"),
                               ("g-pseudo", "in_place", "c3", 0, 2000, "pseudogene"), ("g-199", "in_place", "c4", 0, 199, "gene"),
                               ("g-200", "in_place", "c5", 0, 200, "gene"), ("g-part", "partial", "c6", 0, 2000, "gene"),
                               ("g-unm", "unmapped", "c7", 0, 2000, "gene")])
        comp = next(iter(pool.components.values()))
        self.assertEqual([r["id"] for r in comp.such()], ["g-in", "g-moved", "g-pseudo", "g-200"])
        self.assertEqual(sorted(r["rtype"] for r in comp.such()), ["gene", "gene", "gene", "pseudogene"])
        self.assertEqual({r["id"]: r["cls"] for r in comp.recs.values()}["g-moved"], "moved")
        self.assertEqual({r["id"]: r["cls"] for r in comp.recs.values()}["g-unm"], "unmapped")

    def test_the_called_set_counts_such_records_not_annotated_records(self):
        pool = self.pool_with([("g-a", "in_place", "c1", 0, 2000, "gene"), ("g-short", "in_place", "c2", 0, 150, "gene"), ("g-part", "partial", "c3", 0, 2000, "gene")])
        comp = next(iter(pool.components.values()))
        self.assertEqual((len(comp.record_ids), len(comp.such())), (3, 1))
        self.assertEqual(X.called_components(pool, "nodes"), [])
        pool = self.pool_with([("g-a", "in_place", "c1", 0, 2000, "gene"), ("g-b", "moved", "c2", 0, 2000, "gene")])
        self.assertEqual(len(X.called_components(pool, "nodes")), 1)

    def test_an_unknown_hub_mode_is_an_error(self):
        pool = self.pool_with([("g-a", "in_place", "c1", 0, 2000, "gene"), ("g-b", "in_place", "c2", 0, 2000, "gene")])
        comp = next(iter(pool.components.values()))
        with self.assertRaises(ValueError):
            comp.size("bogus")
        with self.assertRaises(ValueError):
            comp.hub("bogus")

    def test_a_copy_site_id_keeps_the_strand_and_two_strands_are_two_sites(self):
        rows = [row("g-a", contig="c1"), row("g-a", cls="extra_copy", contig="c9", start0=100, end=900, strand="+"),
                row("g-a", cls="extra_copy", contig="c9", start0=100, end=900, strand="-")]
        for r in rows[1:]:
            r["tag"] = "1"
        pool = X.build_pool(rows, {"g-a": rec("g-a", "c1", 0, 1000)})
        comp = next(iter(pool.components.values()))
        self.assertEqual(comp.sites, ["copy@c9:100-900+", "copy@c9:100-900-"])

    def test_the_m2_note_joins_the_record_after_the_words_copy_of(self):
        rows = [row("g-a", contig="c1"), row("g-b", contig="c2"), row("g-a", cls="dropped_M2", contig="c2", start0=0, end=500, note="overlaps copy of g-b (c2)")]
        rows[-1]["tag"] = "1"
        pool = X.build_pool(rows, {"g-a": rec("g-a", "c1", 0, 1000), "g-b": rec("g-b", "c2", 0, 1000)})
        self.assertEqual([sorted(c.record_ids) for c in pool.components.values()], [["g-a", "g-b"]])
        self.assertEqual(pool.components["PC_g-a"].copy_rows[0]["partner"], "g-b")


# ===================================================================================================== D. rank plan and shards
def shard_recs(n, step=100, length=50, contig="c1", prefix="g"):
    return [dict(id=f"{prefix}{i:03d}", contig=contig, start0=i * step, end=i * step + length, idx=i, component=f"P{i}") for i in range(n)]


class ShardTests(T.HermeticCase):
    def sizes(self, shards):
        return [len(ids) for _, ids in shards]

    def test_exactly_cap_records_is_one_shard(self):
        self.assertEqual(self.sizes(X.plan_run_shards("r01", shard_recs(4), cap=4)), [4])
        self.assertEqual(self.sizes(X.plan_run_shards("r01", shard_recs(5), cap=4)), [3, 2])

    def test_blocks_are_balanced_and_numbered_from_one(self):
        s = X.plan_run_shards("r07", shard_recs(9), cap=4)
        self.assertEqual((self.sizes(s), [lab for lab, _ in s]), ([3, 3, 3], ["r07.b1", "r07.b2", "r07.b3"]))
        s = X.plan_run_shards("r07", shard_recs(10), cap=4)
        self.assertEqual(self.sizes(s), [4, 4, 2])

    def test_a_block_ends_only_where_no_record_spans(self):
        recs = shard_recs(10)
        recs[3]["end"] = 480                     # g003 reaches over g004 (start 400): no cut at position 4
        s = X.plan_run_shards("r01", recs, cap=4)
        self.assertEqual(self.sizes(s), [5, 4, 1])
        self.assertEqual(s[0][1], [f"g{i:03d}" for i in range(5)])
        self.assertEqual([i for _, ids in s for i in ids], [f"g{i:03d}" for i in range(10)])

    def test_adjacent_records_may_be_cut_between(self):
        recs = [dict(id=f"a{i}", contig="c1", start0=i * 100, end=(i + 1) * 100, idx=i, component=f"P{i}") for i in range(4)]       # half-open: touching, not overlapping
        self.assertEqual(X.record_free_cuts(recs), [1, 2, 3])

    def test_a_long_record_keeps_its_span_until_it_ends(self):
        recs = [dict(id="L", contig="c1", start0=0, end=1000, idx=0, component="A"), dict(id="s1", contig="c1", start0=100, end=150, idx=1, component="B"),
                dict(id="s2", contig="c1", start0=200, end=250, idx=2, component="C"), dict(id="t", contig="c1", start0=1000, end=1100, idx=3, component="D")]
        self.assertEqual(X.record_free_cuts(recs), [3])

    def test_a_contig_change_is_a_cut_even_inside_the_span_of_a_long_record_of_the_other_contig(self):
        recs = [dict(id="L", contig="c1", start0=0, end=10 ** 6, idx=0, component="A"), dict(id="m", contig="c2", start0=5, end=9, idx=1, component="B")]
        self.assertEqual(X.record_free_cuts(recs), [1])

    def test_records_are_sorted_by_contig_start_end_then_id_before_cutting(self):
        recs = [dict(id="z", contig="c2", start0=0, end=10, idx=0, component="A"), dict(id="b", contig="c1", start0=50, end=60, idx=1, component="B"),
                dict(id="a", contig="c1", start0=50, end=60, idx=2, component="C"), dict(id="c", contig="c1", start0=50, end=55, idx=3, component="D")]
        s = X.plan_run_shards("r01", recs, cap=10)
        self.assertEqual(s[0][1], ["c", "a", "b", "z"])
        by_id = {r["id"]: r for r in recs}
        halves = X.split_shard("r01", [r["id"] for r in recs], by_id)
        self.assertEqual([i for _, ids in halves for i in ids], ["c", "a", "b", "z"])
        self.assertEqual([ids for _, ids in halves], [["c", "a", "b"], ["z"]])

    def test_interleaved_contigs_in_the_input_do_not_make_cuts_inside_a_contig(self):
        recs = [dict(id=f"{c}{i}", contig=c, start0=i * 10, end=i * 10 + 500, idx=0, component="x") for c in ("c1", "c2") for i in range(3)]
        recs.sort(key=lambda r: (r["start0"], r["contig"]))             # start0 first: contigs interleave
        by_id = {r["id"]: r for r in recs}
        halves = X.split_shard("r01", [r["id"] for r in recs], by_id)
        self.assertEqual([ids for _, ids in halves], [["c10", "c11", "c12"], ["c20", "c21", "c22"]])

    def test_split_ties_go_to_the_earlier_cut_and_the_middle_is_not_floored(self):
        def halves(ends):
            recs = shard_recs(5)
            for i, e in ends.items():
                recs[i]["end"] = e
            by_id = {r["id"]: r for r in recs}
            return X.record_free_cuts(recs), self.sizes(X.split_shard("r01", [r["id"] for r in recs], by_id))
        self.assertEqual(halves({}), ([1, 2, 3, 4], [2, 3]))                       # cuts 2 and 3 are equally far from 2.5: the earlier one
        self.assertEqual(halves({1: 250, 2: 350}), ([1, 4], [1, 4]))               # 1 and 4 equally far: the earlier one
        self.assertEqual(halves({1: 250, 3: 450}), ([1, 3], [3, 2]))               # 3 is nearer to 2.5 (not to floor(2.5) = 2, which ties)

    def test_split_prefers_the_cut_nearest_the_middle_not_the_first_or_last(self):
        recs = shard_recs(7)
        recs[1]["end"] = 300                                             # no cut at 2
        by_id = {r["id"]: r for r in recs}
        cuts = X.record_free_cuts(recs)
        self.assertEqual(cuts, [1, 3, 4, 5, 6])
        self.assertEqual(self.sizes(X.split_shard("r01", [r["id"] for r in recs], by_id)), [3, 4])

    def test_the_halves_are_named_s1_s2_and_partition_the_shard(self):
        recs = shard_recs(6)
        by_id = {r["id"]: r for r in recs}
        halves = X.split_shard("r01.b2", [r["id"] for r in recs], by_id)
        self.assertEqual([lab for lab, _ in halves], ["r01.b2.s1", "r01.b2.s2"])
        self.assertEqual([i for _, ids in halves for i in ids], [r["id"] for r in recs])

    def test_a_single_record_shard_cannot_be_split(self):
        recs = shard_recs(1)
        with self.assertRaises(X.NoSplit):
            X.split_shard("r01", [recs[0]["id"]], {recs[0]["id"]: recs[0]})

    def comps(self, n, per=2):
        out = []
        for i in range(n):
            c = X.Component(f"PC_{i:05d}")
            for k in range(per):
                r = rec(f"g{i:05d}-{k}", f"ctg{i:05d}", 0, 2000)
                r["cls"] = "in_place"
                c.recs[r["id"]] = r
                c.record_ids.append(r["id"])
            out.append(c)
        return out

    def test_the_default_cap_is_2500_records_per_shard(self):
        for n, want in ((2200, 1), (2500, 1), (2501, 2)):
            plan = X.rank_plan(self.comps(n))
            self.assertEqual(len(plan[0]["records"]), n)
            self.assertEqual(len(plan[0]["shards"]), want, n)
        plan = X.rank_plan(self.comps(2501))
        self.assertEqual([lab for lab, _ in plan[0]["shards"]], ["r01.b1", "r01.b2"])
        self.assertEqual(sorted(len(ids) for _, ids in plan[0]["shards"]), [1250, 1251])

    def test_an_empty_called_set_has_no_runs(self):
        self.assertEqual(X.rank_plan([]), [])
        self.assertEqual(X.plan_summary(X.rank_plan([]))["J"], 0)

    def test_rank_labels_of_the_control_plans(self):
        plan = X.rank_plan(self.comps(3, per=12), label_prefix="pos")
        self.assertEqual([r["label"] for r in plan][:3], ["pos1", "pos2", "pos3"])
        self.assertEqual(plan[-1]["label"], "pos12")
        plan = X.rank_plan(self.comps(3, per=12))
        self.assertEqual([r["label"] for r in plan][:2] + [plan[-1]["label"]], ["r01", "r02", "r12"])
        plan = X.rank_plan(self.comps(1, per=100))
        self.assertEqual((plan[0]["label"], plan[-1]["label"]), ("r001", "r100"))


# ===================================================================================================== E. outgroup locus and record calls
class OutgroupIndexRuleTests(T.HermeticCase):
    def index(self, specs, **kw):
        """specs: (id, contig, start0, end, exons or None, [copy classes])"""
        rows = []
        for rid, contig, s0, e, exons, copies in specs:
            rows.append(row(rid, contig=contig, start0=s0, end=e, exons=exons))
            for cls in copies:
                rows.append(row(rid, cls=cls, contig="x9", start0=0, end=500, note="-" if cls == "extra_copy" else "overlaps gene-Q (x)"))
        return X.OutgroupIndex(rows, "H", **kw)

    def test_a_long_early_record_is_found_by_a_late_query(self):
        specs = [("long", "h1", 0, 100000, "90000-100000", [])] + [(f"s{i}", "h1", 200 + i * 300, 300 + i * 300, None, []) for i in range(10)]
        m = self.index(specs).match("h1", [(95000, 99000)])
        self.assertEqual((m["status"], m["locus"]), ("ok", "long"))

    def test_a_record_that_starts_inside_the_query_is_a_candidate(self):
        idx = self.index([("a", "h1", 1500, 3500, None, [])])
        self.assertEqual(idx.match("h1", [(1000, 3000)])["locus"], "a")

    def test_only_the_last_exon_block_of_the_query_may_overlap(self):
        idx = self.index([("a", "h1", 9000, 10000, None, [])])
        self.assertEqual(idx.match("h1", [(0, 500), (9000, 10000)])["locus"], "a")

    def test_overlapping_exon_blocks_of_an_outgroup_record_are_merged(self):
        idx = self.index([("a", "h1", 0, 150, "0-100,50-150", [])])
        m = idx.match("h1", [(0, 150)])
        self.assertEqual((m["status"], m["locus"]), ("ok", "a"))
        idx = self.index([("a", "h1", 0, 1000, "0-100,50-150", [])])
        self.assertEqual(idx.match("h1", [(120, 1000)])["status"], "no_outgroup_locus")       # 30 of 150 bases of the union

    def test_every_copy_class_makes_the_record_multi_copy(self):
        for cls in ("extra_copy", "dropped_M1", "dropped_M2"):
            idx = self.index([("a", "h1", 1000, 3000, None, [cls])])
            self.assertEqual(idx.match("h1", [(1000, 3000)])["status"], "outgroup_copy_row", cls)

    def test_a_copy_row_candidate_beside_a_single_copy_candidate_is_a_copy_row_not_ambiguous(self):
        idx = self.index([("single", "h1", 1000, 3000, None, []), ("multi", "h1", 1200, 3200, None, ["extra_copy"])])
        m = idx.match("h1", [(1200, 3000)])
        self.assertEqual((m["status"], m["locus"]), ("outgroup_copy_row", "multi"))
        self.assertEqual(m["candidates"], ["multi", "single"])

    def test_two_copy_row_candidates_are_a_copy_row(self):
        idx = self.index([("m1", "h1", 1000, 3000, None, ["extra_copy"]), ("m2", "h1", 1200, 3200, None, ["dropped_M2"])])
        self.assertEqual(idx.match("h1", [(1200, 3000)])["status"], "outgroup_copy_row")

    def test_candidates_are_sorted_by_id(self):
        idx = self.index([("zz", "h1", 1000, 3000, None, []), ("aa", "h1", 1100, 3100, None, [])])
        self.assertEqual(idx.match("h1", [(1100, 3000)])["candidates"], ["aa", "zz"])

    def test_strict_mode_counts_a_landing_site_of_an_m2_copy_row(self):
        rows = [row("tgt", contig="h1", start0=1000, end=3000), row("src", contig="h2", start0=1000, end=3000),
                row("src", cls="dropped_M2", contig="h1", start0=1000, end=3000, note="overlaps copy of tgt (h1)")]
        self.assertEqual(X.OutgroupIndex(rows, "H").match("h1", [(1000, 3000)])["status"], "ok")
        self.assertEqual(X.OutgroupIndex(rows, "H", strict_single_copy=True).match("h1", [(1000, 3000)])["status"], "outgroup_copy_row")


class RecordCallRuleTests(T.HermeticCase):
    def setUp(self):
        self.idx = T.make_indexes()
        self.none = {"H": set(), "C": set(), "O": set()}

    def lifted(self, h, c, o):
        return {k: {"g1": v} for k, v in T.lifted_for("g1", h, c, o).items()}

    def test_an_ambiguous_outgroup_locus_is_no_outgroup_locus_with_detail_ambiguous(self):
        idx = T.make_indexes(h_extra=[("gene-Hdup", "h1", 1500, 3500, [])])
        c = X.classify_record("g1", self.lifted(("h1", 1500, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000)), self.none, idx)
        self.assertEqual((c["triple"], c["reason"]), (None, "no_outgroup_locus"))
        self.assertEqual((c["per"]["H"]["status"], c["per"]["H"]["detail"], c["per"]["H"]["candidates"]), ("ambiguous", "ambiguous", ["gene-H1", "gene-Hdup"]))

    def test_a_placement_on_a_copy_row_locus_has_the_detail_outgroup_locus_has_copy_row(self):
        c = X.classify_record("g1", self.lifted(("h3", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000)), self.none, self.idx)
        self.assertEqual((c["per"]["H"]["status"], c["per"]["H"]["detail"]), ("outgroup_copy_row", "outgroup_locus_has_copy_row"))
        self.assertEqual(c["per"]["H"]["locus"], "gene-H3")          # the locus is named even though the status is not ok

    def test_a_record_with_only_copy_rows_has_no_primary_row(self):
        lifted = self.lifted(("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000))
        lifted["H"]["g1"] = [lrow("g1_1", 1, "h1", 1000, 3000, cov=0.9, sid=0.97)]
        c = X.classify_record("g1", lifted, self.none, self.idx)
        self.assertEqual((c["per"]["H"]["status"], c["per"]["H"]["detail"]), ("no_placement", "no_primary_row"))

    def test_partial_and_low_identity_together_on_a_single_row_is_partial_only(self):
        lifted = self.lifted(("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000))
        lifted["H"]["g1"] = [lrow("g1", 0, "h1", 1000, 3000, cov=0.4, sid=0.4, partial=1, low=1)]
        c = X.classify_record("g1", lifted, self.none, self.idx)
        self.assertEqual(c["per"]["H"]["detail"], "partial_only")

    def test_coverage_of_exactly_one_half_counts_as_a_qualifying_placement(self):
        lifted = self.lifted(("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000))
        lifted["H"]["g1"] = [lrow("g1", 0, "h1", 1000, 3000, cov=0.5, sid=0.9)]
        c = X.classify_record("g1", lifted, self.none, self.idx)
        self.assertEqual((c["per"]["H"]["n_place"], c["per"]["H"]["status"]), (1, "ok"))
        lifted["H"]["g1"] = [lrow("g1", 0, "h1", 1000, 3000, cov=0.4999, sid=0.9)]
        self.assertEqual(X.classify_record("g1", lifted, self.none, self.idx)["per"]["H"]["n_place"], 0)

    def test_an_unreadable_coverage_is_not_a_qualifying_placement_and_does_not_crash(self):
        lifted = self.lifted(("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000))
        lifted["H"]["g1"] = [lrow("g1", 0, "h1", 1000, 3000, cov="NA", sid="NA")]
        c = X.classify_record("g1", lifted, self.none, self.idx)
        self.assertEqual((c["per"]["H"]["n_place"], c["per"]["H"]["status"]), (0, "no_placement"))

    def test_a_copy_row_of_any_tag_flags_the_record_and_keeps_the_locus_of_its_primary_placement(self):
        # (was test_a_copy_row_of_any_tag_disqualifies: the member-level exclusion is replaced by the registered class-level drop, D2 revised)
        lifted = self.lifted(("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000))
        for tag in (1, 2):
            lifted["H"]["g1"] = [lrow("g1", 0, "h1", 1000, 3000), lrow(f"g1_{tag}", tag, "h2", 1000, 3000, cov=0.9)]
            c = X.classify_record("g1", lifted, self.none, self.idx)
            self.assertEqual((c["triple"], c["reason"], c["multiplicity"]), (("gene-H1", "gene-C1", "gene-O1"), None, ["H"]), tag)
            self.assertEqual((c["per"]["H"]["copy_rows"], c["per"]["H"]["detail"]), (1, "reverse_copy_row"))

    def test_the_human_only_table_does_not_join_records_on_a_multi_copy_human_locus(self):
        comp = T.recs_comp("PC_a", ["a1", "a2", "a3"])
        calls = {}
        for rid, h in (("a1", ("h3", 1000, 3000)), ("a2", ("h3", 1000, 3000)), ("a3", ("h1", 1000, 3000))):
            lifted = {k: {rid: v} for k, v in T.lifted_for(rid, h, ("c1", 1000, 3000), ("o1", 1000, 3000)).items()}
            calls[rid] = X.classify_record(rid, lifted, self.none, self.idx)
        res = X.call_expansions([comp], calls, use=("H",), prefix="EH")
        self.assertEqual(res["expansions"], [])


# ===================================================================================================== F. expansions, siblings, control pool, outputs
class ExpansionRuleTests(T.HermeticCase):
    def test_expansion_ids_follow_the_first_member_not_the_triple(self):
        comp = T.recs_comp("PC_a", ["a1", "a2", "a3", "a4"])
        T1b = ("gene-H9", "gene-C9", "gene-O9")                 # sorts after T2 but its members come first
        res = X.call_expansions([comp], T.calls_for({"a1": T1b, "a2": T1b, "a3": T.T2, "a4": T.T2}))
        self.assertEqual([(e["id"], e["triple"], [m["id"] for m in e["members"]]) for e in res["expansions"]],
                         [("EX0001", T1b, ["a1", "a2"]), ("EX0002", T.T2, ["a3", "a4"])])

    def test_a_triple_with_a_single_expansion_is_no_sibling_group(self):
        res = X.call_expansions([T.recs_comp("PC_a", ["a1", "a2"]), T.recs_comp("PC_b", ["b1", "b2"])], T.calls_for({"a1": T.T1, "a2": T.T1, "b1": T.T2, "b2": T.T2}))
        self.assertEqual(len(res["expansions"]), 2)
        self.assertEqual(res["siblings"], [])

    def test_sibling_groups_are_sorted_by_triple(self):
        comps = [T.recs_comp("PC_a", ["a1", "a2"]), T.recs_comp("PC_b", ["b1", "b2"]), T.recs_comp("PC_c", ["c1", "c2"]), T.recs_comp("PC_d", ["d1", "d2"])]
        calls = T.calls_for({"a1": T.T2, "a2": T.T2, "b1": T.T1, "b2": T.T1, "c1": T.T2, "c2": T.T2, "d1": T.T1, "d2": T.T1})
        res = X.call_expansions(comps, calls)
        self.assertEqual(res["siblings"], [dict(triple=T.T1, expansions=["EX0002", "EX0004"]), dict(triple=T.T2, expansions=["EX0001", "EX0003"])])

    def test_the_guard_looks_at_the_members_of_the_class_only(self):
        comp = T.recs_comp("PC_a", ["a1", "a2", "a3"])
        calls = T.calls_for({"a1": T.T1, "a2": T.T1, "a3": T.T2})
        calls["a3"]["per"]["H"]["n_place"] = 3                   # a record of ANOTHER class has several placements
        res = X.call_expansions([comp], calls)
        self.assertEqual([e["triple"] for e in res["expansions"]], [T.T1])
        self.assertEqual(res["dropped_classes"], [])

    def test_the_guard_looks_at_every_member_not_only_the_first(self):
        comp = T.recs_comp("PC_a", ["a1", "a2", "a3"])
        calls = T.calls_for({"a1": T.T1, "a2": T.T1, "a3": T.T1})
        calls["a3"]["per"]["C"]["n_place"] = 2
        res = X.call_expansions([comp], calls)
        self.assertEqual((res["expansions"], len(res["dropped_classes"])), ([], 1))
        calls = T.calls_for({"a1": T.T1, "a2": T.T1, "a3": T.T1})
        calls["a1"]["per"]["O"]["n_place"] = 2
        self.assertEqual(res["dropped_classes"][0]["members"], ["a1", "a2", "a3"])
        self.assertEqual(len(X.call_expansions([comp], calls)["dropped_classes"]), 1)

    def test_a_component_with_one_triple_record_and_one_unplaced_record_is_no_control_component(self):
        comp = T.recs_comp("PC_a", ["a1", "a2"])
        calls = T.calls_for({"a1": T.T1, "a2": None})
        res = X.call_expansions([comp], calls)
        self.assertEqual(res["components"]["PC_a"]["different_loci"], 1)
        self.assertEqual(X.control_pool_rows([comp], calls, res, joining={}), [])

    def test_a_class_of_one_record_in_a_component_with_one_locus_is_not_an_expansion(self):
        comp = T.recs_comp("PC_a", ["a1", "a2", "a3"])
        calls = T.calls_for({"a1": T.T1, "a2": T.T2, "a3": T.T2})
        res = X.call_expansions([comp], calls)
        self.assertEqual([[m["id"] for m in e["members"]] for e in res["expansions"]], [["a2", "a3"]])
        self.assertEqual(res["components"]["PC_a"]["different_loci"], 2)

    def test_the_classes_of_a_component_with_unrelated_reasons_are_counted_once_each(self):
        comp = T.recs_comp("PC_a", ["a1", "a2", "a3", "a4"])
        calls = T.calls_for({"a1": T.T1, "a2": T.T1, "a3": None, "a4": None})
        calls["a3"]["reason"], calls["a4"]["reason"] = "no_outgroup_locus", "no_placement"
        res = X.call_expansions([comp], calls)
        self.assertEqual(res["reasons"], {"no_outgroup_locus": 1, "no_placement": 1})
        self.assertEqual(res["records_in_expansions"], 2)


class JoiningRowsAndOutputsTests(T.HermeticCase):
    def pool(self):
        rows, recs = T.mini_selflift()
        return X.build_pool(rows, recs), rows

    def test_joining_copy_rows_list_both_ends_sorted_and_skip_unannotated_sites(self):
        pool, _ = self.pool()
        j = X.joining_copy_rows(pool)
        self.assertEqual(j["gene-A"], [("gene-B", "dropped_M1", "0.97"), ("gene-C", "dropped_M1", "0.99")])
        self.assertEqual(j["gene-B"], [("gene-A", "dropped_M1", "0.97")])
        self.assertEqual(j["gene-C"], [("gene-A", "dropped_M1", "0.99")])
        self.assertEqual(j["gene-D"], [("gene-E", "dropped_M2", "0.96")])
        self.assertEqual(j["gene-E"], [("gene-D", "dropped_M2", "0.96")])
        self.assertNotIn("gene-F", j)                      # its only copy row lands on an unannotated site

    def test_two_joining_rows_of_one_record_are_joined_by_semicolons_in_sorted_order(self):
        comp = T.recs_comp("PC_a", ["a1", "a2"])
        calls = T.calls_for({"a1": T.T1, "a2": T.T2})
        res = X.call_expansions([comp], calls)
        joining = {"a1": [("a3", "dropped_M1", "0.97"), ("a2", "extra_copy", "0.99")]}
        rows = X.control_pool_rows([comp], calls, res, joining=joining)
        self.assertEqual(rows[0]["joining_copy_rows"], "a3:dropped_M1:0.97;a2:extra_copy:0.99")

    def write(self, calls=None, hub_mode="nodes", control=(), called=None):
        pool, rows = self.pool()
        called = called if called is not None else X.called_components(pool, hub_mode)
        if calls is None:
            calls = {}
            for c in called:
                for r in c.such():
                    calls[r["id"]] = T.calls_for({r["id"]: T.T1})[r["id"]]
        res = X.call_expansions(called, calls)
        tmp = tempfile.mkdtemp()
        X.write_outputs(res, pool, called, calls, tmp, rows=rows, hub_mode=hub_mode, control=control, extra=dict(note_from_extra=7))
        return tmp, res, pool

    def table(self, tmp, name):
        with open(os.path.join(tmp, name)) as fh:
            return list(csv.reader(fh, delimiter="\t"))

    def test_the_e_table_and_member_exons_are_in_genome_order_not_id_order(self):
        comp = T.recs_comp("PC_a", ["zz", "yy"])           # zz is at 0, yy at 10000
        pool = X.Pool()
        pool.components["PC_a"] = comp
        calls = T.calls_for({"zz": T.T1, "yy": T.T1})
        res = X.call_expansions([comp], calls)
        tmp = tempfile.mkdtemp()
        X.write_outputs(res, pool, [comp], calls, tmp)
        self.assertEqual([r[2] for r in self.table(tmp, "etable.tsv")[1:]], ["c1:0-2000", "c1:10000-12000"])
        self.assertEqual([r[0] for r in self.table(tmp, "member_exons.tsv")[1:]], ["c1:0-2000", "c1:10000-12000"])

    def test_member_exons_are_written_as_merged_blocks(self):
        comp = X.Component("PC_a")
        for rid, c, ex in (("a1", "c1", [(0, 100), (50, 150), (150, 300)]), ("a2", "c2", [(10, 20)])):
            r = rec(rid, c, 0, 1000, exons=ex)
            r["cls"] = "in_place"
            comp.recs[rid] = r
            comp.record_ids.append(rid)
        pool = X.Pool()
        pool.components["PC_a"] = comp
        tmp = tempfile.mkdtemp()
        X.write_outputs(X.call_expansions([], {}), pool, [], {}, tmp)
        got = {r[2]: r[3] for r in self.table(tmp, "member_exons.tsv")[1:]}
        self.assertEqual(got, {"a1": "0-300", "a2": "10-20"})

    def test_extension_sites_are_sorted_and_carry_the_sequence_id(self):
        rows = [row("g-a", contig="c1"), row("g-b", contig="c2")]
        for site, src, sid in (("c9", "g-b", "0.91"), ("c3", "g-a", "0.95")):
            rows.append(row(src, cls="extra_copy", contig=site, start0=100, end=900, seqid=sid, cov="0.99"))
            rows[-1]["tag"] = "1"
        rows.append(row("g-a", cls="dropped_M1", contig="c2", start0=0, end=500, note="overlaps g-b (c2)"))
        rows[-1]["tag"] = "1"
        recs = {"g-a": rec("g-a", "c1", 0, 1000), "g-b": rec("g-b", "c2", 0, 1000)}
        pool = X.build_pool(rows, recs)
        tmp = tempfile.mkdtemp()
        X.write_outputs(X.call_expansions([], {}), pool, [], {}, tmp, rows=rows)
        ext = self.table(tmp, "extension.tsv")
        self.assertEqual([r[0] for r in ext[1:]], ["copy@c3:100-900+", "copy@c9:100-900+"])
        self.assertEqual([r[7] for r in ext[1:]], ["0.95", "0.91"])          # the sequence_id column, not the coverage

    def test_siblings_dropped_classes_and_record_calls_tables(self):
        comps = [T.recs_comp("PC_a", ["a1", "a2"]), T.recs_comp("PC_b", ["b1", "b2"]), T.recs_comp("PC_c", ["c1", "c2", "c3"])]
        calls = T.calls_for({"a1": T.T1, "a2": T.T1, "b1": T.T1, "b2": T.T1, "c1": T.T2, "c2": T.T2, "c3": None})
        calls["c2"]["per"]["C"]["n_place"] = 2                 # the class of c1 c2 is dropped
        calls["c3"]["reason"] = "outgroup_copy_row"
        calls["a1"]["per"]["H"].update(status="ok", detail="", locus="gene-H1")
        calls["c3"]["per"]["H"].update(status="outgroup_copy_row", detail="outgroup_locus_has_copy_row", locus="gene-H3")
        pool = X.Pool()
        for c in comps:
            pool.components[c.id] = c
        res = X.call_expansions(comps, calls)
        tmp = tempfile.mkdtemp()
        X.write_outputs(res, pool, comps, calls, tmp)
        self.assertEqual(self.table(tmp, "siblings.tsv"), [["h", "c", "o", "expansions"], ["gene-H1", "gene-C1", "gene-O1", "EX0001,EX0002"]])
        self.assertEqual(self.table(tmp, "dropped_classes.tsv"),
                         [["pool_component_id", "h", "c", "o", "members", "reason"], ["PC_c", "gene-H2", "gene-C2", "gene-O2", "c1,c2", "multiplicity"]])
        rc = self.table(tmp, "record_calls.tsv")
        self.assertEqual(rc[1][:5], ["c1:0-2000", "a1", "PC_a", "", "gene-H1,gene-C1,gene-O1"])
        hdr = rc[0]
        c3 = next(r for r in rc[1:] if r[1] == "c3")
        self.assertEqual((c3[3], c3[4]), ("outgroup_copy_row", ""))
        self.assertEqual((c3[hdr.index("H_status")], c3[hdr.index("H_detail")], c3[hdr.index("H_locus")]), ("outgroup_copy_row", "outgroup_locus_has_copy_row", ""))
        a1 = next(r for r in rc[1:] if r[1] == "a1")
        self.assertEqual((a1[hdr.index("H_status")], a1[hdr.index("H_n_place")], a1[hdr.index("H_copy_rows")], a1[hdr.index("H_locus")]), ("ok", "1", "0", "gene-H1"))
        self.assertEqual(hdr[5:10] + hdr[10:15], ["H_status", "H_detail", "H_n_place", "H_copy_rows", "H_locus", "C_status", "C_detail", "C_n_place", "C_copy_rows", "C_locus"])

    def test_summary_counts_and_the_hash_list(self):
        tmp, res, pool = self.write()
        s = json.loads(rd(os.path.join(tmp, "summary.json")))
        self.assertEqual(s["note_from_extra"], 7)
        self.assertEqual(sorted(s["sha256"]), sorted(["etable.tsv", "member_exons.tsv", "control_pool.tsv", "extension.tsv", "siblings.tsv", "record_calls.tsv",
                                                      "dropped_classes.tsv", "components.tsv", "etable_human_only.tsv"]))
        for n, h in s["sha256"].items():
            self.assertEqual(h, hashlib.sha256(open(os.path.join(tmp, n), "rb").read()).hexdigest(), n)
        keys = [m.group(1) for m in re.finditer(r'^ "([^"]+)":', rd(os.path.join(tmp, "summary.json")), re.M)]
        self.assertEqual(keys, sorted(keys))

    def test_summary_counts_with_several_components(self):
        comps = [T.recs_comp("PC_a", ["a1", "a2"]), T.recs_comp("PC_b", ["b1", "b2", "b3"]), T.recs_comp("PC_c", ["c1", "c2"]), T.recs_comp("PC_d", ["d1", "d2"])]
        calls = T.calls_for({"a1": T.T1, "a2": T.T1,                           # a: one locus, an expansion
                             "b1": T.T1, "b2": T.T2, "b3": T.T2,               # b: two loci
                             "c1": T.T1, "c2": None,                           # c: one locus, one unplaced
                             "d1": None, "d2": None})                          # d: nothing
        pool = X.Pool()
        for c in comps:
            pool.components[c.id] = c
        pool.components["PC_hub"] = T.recs_comp("PC_hub", [f"h{i}" for i in range(60)])
        res = X.call_expansions(comps, calls)
        tmp = tempfile.mkdtemp()
        X.write_outputs(res, pool, comps, calls, tmp)
        s = json.loads(rd(os.path.join(tmp, "summary.json")))
        self.assertEqual((s["called_components"], s["expansions"], s["records_in_expansions"]), (4, 2, 4))
        self.assertEqual((s["components_with_a_triple"], s["components_with_different_loci"], s["classes_in_components_with_different_loci"]), (3, 1, 2))
        self.assertEqual(s["hubs"], 1)
        self.assertEqual(s["records_without_triple_by_reason"], {"no_placement": 3})
        cp = self.table(tmp, "control_pool.tsv")
        self.assertEqual((s["control_pool_rows"], s["control_pool_rows_eligible"]), (len(cp) - 1, sum(int(r[4]) for r in cp[1:])))
        by = {r[0]: r for r in self.table(tmp, "components.tsv")[1:]}
        self.assertEqual(by["PC_hub"][5], "hub")
        self.assertEqual(by["PC_a"][5], "verdict")

    def test_control_pool_eligibility_is_counted_in_the_summary(self):
        comp = T.recs_comp("PC_a", ["a1", "a2", "a3"])
        calls = T.calls_for({"a1": T.T1, "a2": T.T2, "a3": None})
        pool = X.Pool()
        pool.components["PC_a"] = comp
        res = X.call_expansions([comp], calls)
        tmp = tempfile.mkdtemp()
        X.write_outputs(res, pool, [comp], calls, tmp)
        s = json.loads(rd(os.path.join(tmp, "summary.json")))
        self.assertEqual((s["control_pool_rows"], s["control_pool_rows_eligible"]), (3, 2))


# ===================================================================================================== G. Gate 7 helpers and the planted negative control
class PositiveControlRuleTests(T.HermeticCase):
    def setUp(self):
        self.pool, self.rows, self.spec = T.lrpap_fixture()

    def test_exactly_half_of_the_exon_bases_qualifies(self):
        recs = T.lrpap_fixture_records()
        recs["gene-HALF"] = rec("gene-HALF", "c2", 5000, 6000)             # the landing interval of the reference on c2 is 5000-7000: 1000 of 2000 would be 50%
        recs["gene-HALF"] = rec("gene-HALF", "c2", 6000, 8000)             # 6000-7000 inside the 5000-7000 interval = 1000 of 2000
        rows = list(self.rows) + [row("gene-HALF", contig="c2", start0=6000, end=8000),
                                  row("gene-HALF", cls="dropped_M1", contig="c1", start0=0, end=100, note="overlaps gene-REF (c1)")]
        rows[-1]["tag"] = "1"
        pool = X.build_pool(rows, recs)
        ctrl = X.lrpap_control(pool, rows, self.spec)
        self.assertEqual(ctrl["shares"]["gene-HALF"], 0.5)
        self.assertIn("gene-HALF", ctrl["qualifying"])
        self.assertAlmostEqual(X.lrpap_control(pool, rows, dict(self.spec, min_share=0.5000001))["shares"]["gene-HALF"], 0.5)
        self.assertNotIn("gene-HALF", X.lrpap_control(pool, rows, dict(self.spec, min_share=0.5000001))["qualifying"])

    def test_the_reference_must_be_in_the_pool(self):
        with self.assertRaises(ValueError):
            X.lrpap_control(self.pool, self.rows, dict(self.spec, reference="gene-NOPE"))

    def test_fragments_that_are_all_unknown_are_an_error(self):
        with self.assertRaises(ValueError):
            X.lrpap_control(self.pool, self.rows, dict(self.spec, fragments=["gene-NOPE1", "gene-NOPE2"]))

    def test_the_named_not_required_records_are_listed_only_when_they_do_not_qualify(self):
        ctrl = X.lrpap_control(self.pool, self.rows, dict(self.spec, not_required=["gene-X2", "gene-X1"]))
        self.assertEqual(ctrl["not_required"], ["gene-X2"])
        self.assertEqual(ctrl["spec_problems"], ["gene-X1"])
        self.assertEqual(ctrl["other_non_qualifying"], ["gene-N"])

    def test_the_copy_interval_of_the_in_place_row_counts(self):
        rows = [r for r in self.rows if not (r["source_id"] == "gene-REF" and r["cls"] == "in_place")]
        ctrl = X.lrpap_control(X.build_pool(self.rows, T.lrpap_fixture_records()), rows, self.spec)
        self.assertNotIn("gene-REF", ctrl["qualifying"])
        rows = [r for r in self.rows if not (r["source_id"] == "gene-REF" and r["cls"] == "dropped_M1")]
        ctrl = X.lrpap_control(X.build_pool(self.rows, T.lrpap_fixture_records()), rows, self.spec)
        self.assertEqual(ctrl["qualifying"], ["gene-REF"])

    def called(self):
        return X.called_components(self.pool, "nodes")

    def calls(self, **triples):
        calls = {r["id"]: dict(triple=None, reason="no_placement", per={k: dict(n_place=1, detail="") for k in "HCO"}) for c in self.called() for r in c.such()}
        for rid, t in triples.items():
            calls[rid.replace("_", "-")]["triple"], calls[rid.replace("_", "-")]["reason"] = t, None
        return calls

    def test_a_big_class_that_misses_a_qualifying_record_fails(self):
        ctrl = X.lrpap_control(self.pool, self.rows, self.spec)
        calls = self.calls(**{"gene_REF": T.T1, "gene_X1": T.T1, "gene_X2": T.T1, "gene_N": T.T1, "gene_F1": T.T1, "gene_F2": T.T1})      # 6 members, gene-X3 missing
        res = X.call_expansions(self.called(), calls)
        out = X.check_positive_control(ctrl, res, calls)
        self.assertFalse(out["pass"])
        self.assertEqual(out["missing"], [dict(record="gene-X3", reason="no_placement")])

    def test_the_reasons_for_a_missing_record(self):
        ctrl = X.lrpap_control(self.pool, self.rows, self.spec)
        calls = self.calls(**{"gene_REF": T.T1, "gene_X1": T.T1, "gene_X3": T.T2})
        res = X.call_expansions(self.called(), calls)
        out = X.check_positive_control(ctrl, res, calls)
        self.assertFalse(out["pass"])
        self.assertEqual(out["missing"], [dict(record="gene-X3", reason="no partner with the same triple")])
        self.assertEqual(out["expansion"], res["expansions"][0]["id"])

    def test_the_best_partial_class_is_the_earliest_among_equals(self):
        ctrl = X.lrpap_control(self.pool, self.rows, self.spec)
        calls = self.calls(**{"gene_REF": T.T1, "gene_X1": T.T1, "gene_X3": T.T2, "gene_N": T.T2})            # two classes with two members each; one of them holds REF
        res = X.call_expansions(self.called(), calls)
        out = X.check_positive_control(ctrl, res, calls)
        self.assertEqual(out["expansion"], "EX0001")
        self.assertEqual(sorted(out["class_members"]), ["gene-REF", "gene-X1"])
        calls = self.calls(**{"gene_REF": T.T1, "gene_X1": T.T1, "gene_N": T.T1, "gene_X3": T.T2, "gene_X2": T.T2})
        res = X.call_expansions(self.called(), calls)
        out = X.check_positive_control(ctrl, res, calls)
        self.assertEqual(sorted(out["class_members"]), ["gene-N", "gene-REF", "gene-X1"])          # the larger overlap with Q wins, not the larger class
        self.assertEqual(out["missing"], [dict(record="gene-X3", reason="in another class")])

    def test_not_required_records_report_their_reason_or_that_they_were_not_called(self):
        ctrl = X.lrpap_control(self.pool, self.rows, self.spec)
        calls = self.calls(**{"gene_REF": T.T1, "gene_X1": T.T1, "gene_X3": T.T1})
        calls["gene-X2"]["reason"] = None
        res = X.call_expansions(self.called(), calls)
        out = X.check_positive_control(ctrl, res, calls)
        self.assertEqual(out["not_required"], [dict(record="gene-X2", in_class=False, reason="no_partner")])
        self.assertTrue(out["pass"])
        del calls["gene-X2"]
        out = X.check_positive_control(ctrl, res, calls)
        self.assertEqual(out["not_required"], [dict(record="gene-X2", in_class=False, reason="not called")])

    def test_split_control_removes_each_component_once(self):
        ctrl = X.lrpap_control(self.pool, self.rows, self.spec)
        same = dict(ctrl, fragments=dict(ctrl["fragments"], component=ctrl["component"]))
        verdict, control = X.split_control(self.called(), same)
        self.assertEqual([c.id for c in control], [ctrl["component"]])
        verdict, control = X.split_control(self.called(), ctrl)
        self.assertEqual(len(control), 2)

    def test_the_fragments_leave_the_verdict_set_with_the_control(self):
        ctrl = X.lrpap_control(self.pool, self.rows, self.spec)
        verdict, control = X.split_control(self.called(), ctrl)
        self.assertFalse({c.id for c in verdict} & {c.id for c in control})
        self.assertIn(ctrl["fragments"]["component"], {c.id for c in control})


class PlantedNegativeRuleTests(T.HermeticCase):
    def world(self, per_contig=12, n_contigs=26):
        return T.NegativeControlTests().world(n_contigs=n_contigs, per_contig=per_contig)

    def test_the_pairs_follow_the_registered_greedy_rule(self):
        rows, recs = self.world()
        pool = X.build_pool(rows, recs)
        el = X.eligible_for_planting(rows, recs, pool)
        chosen = random.Random(20261007).sample(el, 100)
        left, want = list(chosen), []
        while left:
            a = left.pop(0)
            j = next(k for k, b in enumerate(left) if recs[b]["contig"] != recs[a]["contig"])
            want.append(sorted([a, left.pop(j)]))
        comps = X.plant_negative(rows, recs, pool)
        self.assertEqual([sorted(r["id"] for r in c.such()) for c in comps], want)
        self.assertEqual([c.id for c in comps][:2] + [comps[-1].id], ["NEG01", "NEG02", "NEG50"])

    def test_exactly_the_requested_number_of_eligible_records_is_enough(self):
        rows, recs = [], {}
        for c in range(4):                                                      # 4 contigs x 2 long records, none in the pool
            for i in range(2):
                rid = f"gene-e{c}x{i}"
                recs[rid] = rec(rid, f"ctg{c}", i * 5000, i * 5000 + 2000)
                rows.append(row(rid, contig=f"ctg{c}", start0=i * 5000, end=i * 5000 + 2000))
        pool = X.build_pool(rows, recs)
        self.assertEqual(len(X.eligible_for_planting(rows, recs, pool)), 8)
        ok = []
        for seed in range(1, 60):
            try:
                ok.append((seed, X.plant_negative(rows, recs, pool, seed=seed, n_records=8)))
            except ValueError as e:
                self.assertNotIn("eligible", str(e))                           # a pairing failure, never the size check
        self.assertTrue(ok)
        self.assertEqual(len(ok[0][1]), 4)
        with self.assertRaises(ValueError) as cm:
            X.plant_negative(rows, recs, pool, n_records=9)
        self.assertIn("only 8 eligible records for 9", str(cm.exception))

    def test_records_on_one_contig_cannot_be_paired(self):
        rows, recs = [], {}
        for i in range(4):
            rid = f"gene-s{i}"
            recs[rid] = rec(rid, "ctg0", i * 5000, i * 5000 + 2000)
            rows.append(row(rid, contig="ctg0", start0=i * 5000, end=i * 5000 + 2000))
        with self.assertRaises(ValueError) as cm:
            X.plant_negative(rows, recs, X.build_pool(rows, recs), n_records=4)
        self.assertIn("share a contig", str(cm.exception))

    def test_the_200_bp_edge_and_copy_classes_of_the_planting_pool(self):
        rows = [row("g199", contig="c1", start0=0, end=199), row("g200", contig="c2", start0=0, end=200), row("gmoved", cls="moved", contig="c3", start0=0, end=2000),
                row("gm2", contig="c4", start0=0, end=2000), row("gz", contig="c5", start0=0, end=2000),
                row("gm2", cls="dropped_M2", contig="c5", start0=0, end=500, note="overlaps copy of gz (c5)"), row("gok", contig="c6", start0=0, end=2000)]
        rows[-2]["tag"] = "1"
        recs = {r["source_id"]: rec(r["source_id"], r["contig"], 0, int(r["end"])) for r in rows}
        pool = X.build_pool(rows, recs)
        self.assertEqual(sorted(pool.node_component), ["gm2", "gz"])
        self.assertEqual(X.eligible_for_planting(rows, recs, pool), ["g200", "gok"])
        self.assertEqual(X.planting_counts(rows, pool), dict(in_place_ge200=4, no_copy_row_as_source=3, target_of_a_copy_row=1, eligible=2))

    def test_the_planting_pool_is_sorted_ids(self):
        rows, recs = self.world(per_contig=4, n_contigs=3)
        rows.reverse()
        el = X.eligible_for_planting(rows, recs, X.build_pool(rows, recs))
        self.assertEqual(el, sorted(el))

    def test_the_negative_control_bar(self):
        res = dict(expansions=[1])
        self.assertEqual(X.check_negative(res, 50), dict(called=1, of=50, bar=1, **{"pass": True}))
        self.assertTrue(X.check_negative(dict(expansions=[]), 50)["pass"])
        self.assertFalse(X.check_negative(dict(expansions=[1, 2]), 50)["pass"])
        self.assertTrue(X.check_negative(dict(expansions=[1, 2]), 50, bar=2)["pass"])


class CompareExpectedTests(T.HermeticCase):
    def report(self, runs=(300, 300, 165, 116, 81, 70, 60)):
        return dict(components_size_ge2=609, annotated_records_in_components=2976, annotated_ge200bp=2249, components_ge2_such=307,
                    nodes=dict(hubs=7, called=300, records=1403, by_rtype={"gene": 907, "pseudogene": 496}, J=44,
                               plan=dict(records_per_run=list(runs), runs_with_at_most_10=26)),
                    supported={"gorilla_KB3781": dict(components=80, hubs_nodes=4, hubs_records=4), "gorilla_OR6737": dict(components=139, hubs_nodes=5, hubs_records=5),
                               "either": dict(components=157, hubs_nodes=5, hubs_records=5)})

    def test_a_long_run_list_matches_on_its_first_five(self):
        self.assertEqual(X.compare_expected(self.report(), "nodes"), [])
        rep = self.report(runs=(300, 300, 165, 116, 80, 70))
        self.assertEqual(X.compare_expected(rep, "nodes"), ["records_per_run_first5: expected [300, 300, 165, 116, 81], got [300, 300, 165, 116, 80]"])
        rep = self.report(runs=(300, 300, 165, 116, 81))
        rep["nodes"]["plan"]["records_per_run"] = [300, 300, 165, 116, 81]
        self.assertEqual(X.compare_expected(rep, "nodes"), [])

    def test_every_registered_number_is_compared_in_both_directions(self):
        for k, v in (("runs_with_at_most_10", 25), ("runs_with_at_most_10", 27)):
            rep = self.report()
            rep["nodes"]["plan"][k] = v
            self.assertEqual(len(X.compare_expected(rep, "nodes")), 1, (k, v))
        rep = self.report()
        rep["supported"]["either"]["components"] = 156
        self.assertEqual(X.compare_expected(rep, "nodes"), ["supported either: expected 157 components, got 156"])
        rep = self.report()
        del rep["supported"]["either"]
        self.assertEqual(X.compare_expected(rep, "nodes"), [])                     # a library that was not read is not a difference
        rep = self.report()
        rep["supported"]["gorilla_KB3781"]["hubs_nodes"] = 3
        self.assertEqual(X.compare_expected(rep, "nodes"), ["supported gorilla_KB3781: expected 4 hubs, got 3"])



# ===================================================================================================== H. the watchdog: environment, log, kill, memory poll
class WatchBase(T.HermeticCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.started = []

    def tearDown(self):
        for pid in self.started:
            try:
                os.killpg(pid, signal.SIGKILL)
            except OSError:
                pass

    def stub(self, name, body):
        return write_exec(os.path.join(self.tmp, name), "#!/bin/bash\n" + body)

    def watch(self, cmd, **kw):
        kw.setdefault("hard_s", 30)
        kw.setdefault("poll_s", 0.1)
        kw.setdefault("env", {})
        t0 = time.time()
        rc, peak = LC.timed_watch([cmd] if isinstance(cmd, str) else cmd, os.path.join(self.tmp, "time.txt"), os.path.join(self.tmp, "log.txt"), cwd=self.tmp, **kw)
        return rc, peak, time.time() - t0

    def pid_of(self, name="pid.txt"):
        return int(rd(os.path.join(self.tmp, name)).strip())


class TimedWatchTests(WatchBase):
    def test_the_child_sees_the_liftoff_environment_first_on_the_path_and_the_registered_tmpdir(self):
        stub = self.stub("env.sh", 'printf "%s" "$PATH" > path.txt\nprintf "%s" "${TMPDIR-UNSET}" > tmpdir.txt\n')
        with mock.patch.dict(os.environ, {"PATH": "/usr/bin:/bin"}):
            os.environ.pop("TMPDIR", None)
            self.assertEqual(self.watch(stub)[0], 0)
            self.assertEqual(rd(os.path.join(self.tmp, "path.txt")), f"{LC.L.LIFTOFF_ENV}/bin:/usr/bin:/bin")
            self.assertEqual(rd(os.path.join(self.tmp, "tmpdir.txt")), "/mnt/linuxdisk/tmp")
            os.environ["TMPDIR"] = "/somewhere/else"
            self.watch(stub)
            self.assertEqual(rd(os.path.join(self.tmp, "tmpdir.txt")), "/somewhere/else")
            self.watch(stub, env={"TMPDIR": "/from/caller", "EXTRA": "1"})
            self.assertEqual(rd(os.path.join(self.tmp, "tmpdir.txt")), "/from/caller")

    def test_the_log_holds_both_streams_and_the_time_record_the_exit_status(self):
        stub = self.stub("out.sh", 'echo to-stdout\necho to-stderr 1>&2\nexit 0\n')
        self.assertEqual(self.watch(stub)[0], 0)
        log = rd(os.path.join(self.tmp, "log.txt"))
        self.assertIn("to-stdout", log)
        self.assertIn("to-stderr", log)
        self.assertEqual(LC.read_time(os.path.join(self.tmp, "time.txt"))["exit"], 0.0)

    def test_a_run_past_the_hard_limit_is_killed_on_time_with_nothing_left_running(self):
        stub = self.stub("long.sh", 'echo $$ > pid.txt\nexec sleep 30\n')
        RecordingPopen.waits = []
        t_start = time.monotonic()
        with mock.patch.object(LC.subprocess, "Popen", RecordingPopen):
            rc, peak, elapsed = self.watch(stub, hard_s=1, poll_s=8)
        self.assertEqual(rc, LC.RC_TOO_LONG)
        loop = RecordingPopen.waits[:-1]                                # the last wait of the supervised process reaps it after the kill
        polls = [t for t, at in loop]
        self.assertTrue(polls and max(polls) <= 1.0, polls)          # the poll interval must not delay the kill beyond the hard limit
        self.assertGreaterEqual(max(at + t for t, at in loop) - t_start, 1.0)       # and the kill is not earlier than the hard limit (recorded monotonic stamps, no race)
        self.assertFalse(alive(self.pid_of()))                        # the session is verified empty before timed_watch returns: no waiting needed

    def test_a_run_over_the_memory_limit_is_killed_with_nothing_left_running(self):
        stub = self.stub("long.sh", 'echo $$ > pid.txt\nexec sleep 30\n')
        with mock.patch.object(LC, "group_rss_kb", return_value=1_050_000):
            rc, peak, elapsed = self.watch(stub, rss_limit_gb=1.0)
        self.assertEqual((rc, peak), (LC.RC_TOO_BIG, 1_050_000))
        self.assertFalse(alive(self.pid_of()))                        # timed_watch verifies the session is empty before it returns

    def test_the_memory_limit_is_in_decimal_gigabytes_and_exclusive(self):
        stub = self.stub("short.sh", 'exec sleep 0.6\n')
        for kb, want in ((900_000, 0), (950_000, 0), (1_000_000, 0), (1_000_001, LC.RC_TOO_BIG), (1_050_000, LC.RC_TOO_BIG), (20_000_000, LC.RC_TOO_BIG)):
            with mock.patch.object(LC, "group_rss_kb", return_value=kb):
                rc, peak, _ = self.watch(stub, rss_limit_gb=1.0)
            self.assertEqual(rc, want, kb)
            self.assertEqual(peak, kb)

    def test_a_run_longer_than_a_poll_but_inside_the_limits_completes_and_reports_the_highest_rss(self):
        stub = self.stub("short.sh", 'exec sleep 0.8\n')
        vals = iter([100, 700])
        with mock.patch.object(LC, "group_rss_kb", side_effect=lambda pid: next(vals, 300)):
            rc, peak, elapsed = self.watch(stub, poll_s=0.15)
        self.assertEqual((rc, peak), (0, 700))

    def test_a_bad_ps_line_is_skipped_and_the_poll_follows_the_interval(self):
        stub = self.stub("short.sh", 'exec sleep 0.8\n')
        seen = []

        def rss(pid):
            seen.append(pid)
            if len(seen) == 1:
                raise ValueError("bad ps line")
            return 200
        with mock.patch.object(LC, "group_rss_kb", side_effect=rss):
            rc, peak, _ = self.watch(stub, poll_s=0.2)
        self.assertEqual((rc, peak), (0, 200))
        self.assertTrue(2 <= len(seen) <= 6, len(seen))
        seen.clear()
        with mock.patch.object(LC, "group_rss_kb", side_effect=lambda pid: seen.append(pid) or 5):
            self.watch(self.stub("one.sh", 'exec sleep 1.0\n'), poll_s=0.4)
        self.assertTrue(1 <= len(seen) <= 4, len(seen))


class KillGroupTests(WatchBase):
    """kill_session (replaces _kill_group): SIGTERM first, SIGKILL for what ignores it, and the session is empty when it returns"""

    def group(self, body):
        p = subprocess.Popen(["/bin/bash", "-c", body], cwd=self.tmp, start_new_session=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.started.append(p.pid)
        pids = []
        for _ in range(400):                                             # until the leader and its two children are up (readiness, not a fixed delay)
            pids = [int(x) for x in subprocess.run(["ps", "-o", "pid=", "--sid", str(p.pid)], capture_output=True, text=True).stdout.split()]
            if len(pids) >= 3:
                break
            time.sleep(0.05)
        self.assertGreaterEqual(len(pids), 3, pids)                      # the leader and its two children
        return p, pids

    def signals(self):
        """patch os.kill / os.killpg to record the signals they send (and still send them)"""
        sent, real_killpg, real_kill = [], os.killpg, os.kill
        return (sent, mock.patch.object(LC.os, "killpg", lambda pg, sig: (sent.append(sig), real_killpg(pg, sig))[1]),
                mock.patch.object(LC.os, "kill", lambda pid, sig: (sent.append(sig), real_kill(pid, sig))[1]))

    def test_the_whole_session_dies_on_the_first_signal(self):
        p, pids = self.group("sleep 60 & sleep 60 & wait")
        sent, a, b = self.signals()
        with a, b:
            left = LC.kill_session(p.pid, term_grace_s=5, kill_wait_s=5)
        self.assertEqual(left, [])
        self.assertEqual(set(sent), {signal.SIGTERM})                    # nothing needed more than SIGTERM
        p.wait()
        for pid in pids:
            self.assertFalse(alive(pid), pid)

    def test_a_session_that_ignores_sigterm_is_killed_for_good_after_the_grace_period(self):
        p, pids = self.group('trap "" TERM; sleep 60 & sleep 60 & wait')
        sent, a, b = self.signals()
        with a, b:
            left = LC.kill_session(p.pid, term_grace_s=0.3, kill_wait_s=5)
        self.assertEqual(left, [])
        self.assertEqual(sent[0], signal.SIGTERM)
        self.assertIn(signal.SIGKILL, sent)
        self.assertLess(sent.index(signal.SIGTERM), sent.index(signal.SIGKILL))
        p.wait()
        for pid in pids:
            self.assertFalse(alive(pid), pid)

    def test_an_empty_session_is_left_alone_and_pid_one_is_refused(self):
        p = subprocess.Popen(["/bin/true"], start_new_session=True)
        p.wait()
        sent, a, b = self.signals()
        with a, b:
            self.assertEqual(LC.kill_session(p.pid), [])
            # fix round 2 (S1): the group gets its SIGTERM before any listing, so for a session that is already empty one harmless SIGTERM goes to a group that no
            # longer exists (ESRCH); no SIGKILL, and nothing is sent to any single process
            self.assertEqual(sent, [signal.SIGTERM])
            sent.clear()
            for bad in (0, 1):
                with self.assertRaises(ValueError):
                    LC.kill_session(bad)
            self.assertEqual(sent, [])


class GroupRssTests(T.HermeticCase):
    def test_the_resident_sets_of_the_whole_session_are_summed(self):
        # fix round 2 (S1): the sizes are read from /proc (LC._session_members), not from a forked `ps`; a zombie is not counted
        members = [(11, "S", 100), (12, "R", 200), (13, "S", 300), (14, "Z", 0)]
        with mock.patch.object(LC, "_session_members", return_value=members) as listing, mock.patch.object(LC.subprocess, "run", side_effect=AssertionError("ps was forked")):
            self.assertEqual(LC.group_rss_kb(4242), 600)
            self.assertEqual(LC.session_pids(4242), [11, 12, 13])
        self.assertEqual(listing.call_args_list[0][0], (4242,))

    def test_a_real_session_with_two_children_is_summed_not_maximised(self):
        p = subprocess.Popen(["/bin/bash", "-c", "sleep 5 & sleep 5 & wait"], start_new_session=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        try:
            time.sleep(0.4)
            each = [int(x) for x in subprocess.run(["ps", "-o", "rss=", "--sid", str(p.pid)], capture_output=True, text=True).stdout.split()]
            self.assertEqual(len(each), 3)
            got = LC.group_rss_kb(p.pid)
            self.assertGreaterEqual(got, sum(each) - 64)
            self.assertGreater(got, max(each))
        finally:
            os.killpg(p.pid, signal.SIGKILL)
            p.wait()


class CheckRssTests(T.HermeticCase):
    def test_the_default_limit_is_16_gb_decimal_and_exclusive(self):
        LC.check_rss({"peak_rss_kb": 16_000_000})
        LC.check_rss({"peak_rss_kb": 15_900_000})
        LC.check_rss({})
        for kb in (16_000_001, 16_100_000, 16_500_000, 30_000_000):
            with self.assertRaises(LC.StopRun):
                LC.check_rss({"peak_rss_kb": kb})
        with self.assertRaises(LC.StopRun):
            LC.check_rss({"peak_rss_kb": 2_100_000}, limit_gb=2.0)
        LC.check_rss({"peak_rss_kb": 2_000_000}, limit_gb=2.0)


# ===================================================================================================== I. one shard: stale files, skip rule, results
class RunShardRuleTests(T.HermeticCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.lifted = os.path.join(self.tmp, "stub_lifted.gff3")
        with open(self.lifted, "w") as fh:
            fh.write(T.lifted_gff([("gene-A", "h1", 1001, 2000, "+", "0.98", "0.97", 0, "False", "False", [(1001, 1400), (1701, 2000)], "gene"),
                                   ("gene-A_1", "h2", 501, 1500, "-", "0.96", "0.95", 1, "False", "False", [(501, 1500)], "gene")]))
        self.d = os.path.join(self.tmp, "run")
        os.makedirs(self.d)
        with open(os.path.join(self.d, "in.gff"), "w") as fh:
            fh.write(T.REF_GFF)
        self.env = {"STUB_LIFTED": self.lifted}
        self.flags = write_exec(os.path.join(self.tmp, "flags.sh"),
                                "#!/bin/bash\nout=lifted.gff3; while [ $# -gt 0 ]; do case \"$1\" in -o) out=$2; shift 2;; *) shift;; esac; done\n"
                                "[ -e in.gff_db ] && : > seen_db\n[ -e lifted.gff3 ] && : > seen_lifted\n[ -e done.key ] && : > seen_donekey\n"
                                "mkdir -p inter; echo x > inter/f; : > in.gff_db\ncp \"$STUB_LIFTED\" \"$out\"; printf 'gene-B\\ngene-P\\ngene-Z\\n' > unmapped.txt\n")
        self.noout = write_exec(os.path.join(self.tmp, "noout.sh"), "#!/bin/bash\nexit 0\n")

    def run_it(self, stub=None, **kw):
        kw.setdefault("hard_s", 60)
        kw.setdefault("env", self.env)
        kw.setdefault("target", "/t.fa")
        return LC.run_shard(self.d, ["gene-A", "gene-B"], liftoff_bin=stub or T.stub_liftoff(self.tmp), reference="/r.fa", minimap2="/m", **kw)

    def p(self, name):
        return os.path.join(self.d, name)

    def test_a_missing_in_gff_is_a_stop(self):
        os.remove(self.p("in.gff"))
        with self.assertRaises(LC.StopRun):
            self.run_it()

    def test_the_result_reports_features_time_and_memory(self):
        res = self.run_it()
        self.assertEqual((res["status"], res["features"]), ("done", 2))
        self.assertEqual(res["wall_s"], LC.read_time(self.p("time.txt"))["wall_s"])
        self.assertIsNotNone(res["wall_s"])
        self.assertGreater(res["peak_rss_kb"], 0)
        self.assertEqual(res["peak_rss_kb"], LC.read_time(self.p("time.txt"))["peak_rss_kb"])

    def test_a_run_is_skipped_only_when_loci_and_a_matching_done_key_are_both_there(self):
        self.assertEqual(self.run_it()["status"], "done")
        self.assertEqual(self.run_it()["status"], "skipped")
        os.remove(self.p("loci.tsv"))
        self.assertEqual(self.run_it()["status"], "done")
        os.remove(self.p("done.key"))
        self.assertEqual(self.run_it()["status"], "done")
        os.remove(self.p("done.key"))
        os.remove(self.p("loci.tsv"))
        self.assertEqual(self.run_it()["status"], "done")

    def test_a_loci_table_without_a_done_key_is_not_a_finished_run_and_does_not_crash(self):
        with open(self.p("loci.tsv"), "w") as fh:
            fh.write("half a run\n")
        self.assertEqual(self.run_it()["status"], "done")
        self.assertNotIn("half a run", rd(self.p("loci.tsv")))

    def test_another_input_or_other_parameters_redo_a_finished_run(self):
        self.assertEqual(self.run_it()["status"], "done")
        with open(self.p("in.gff"), "a") as fh:
            fh.write("NC_1\tRefSeq\tgene\t20001\t21000\t.\t+\t.\tID=gene-Q;Name=Q\n")
        self.assertEqual(self.run_it()["status"], "done")
        self.assertEqual(self.run_it()["status"], "skipped")
        with mock.patch.object(LC.L, "params_key", return_value={"liftoff": "other"}):
            self.assertEqual(self.run_it()["status"], "done")
        self.assertEqual(self.run_it()["status"], "done")                       # and back: the key records the parameters

    def test_nothing_of_an_earlier_attempt_is_left_when_liftoff_starts(self):
        for name in ("in.gff_db", "lifted.gff3", "lifted.gff3.gz", "unmapped.txt", "done.key", "toolong"):
            with open(self.p(name), "w") as fh:
                fh.write("stale\n")
        os.makedirs(self.p("inter"))
        self.run_it(stub=self.flags)
        for name in ("seen_db", "seen_lifted", "seen_donekey"):
            self.assertFalse(os.path.exists(self.p(name)), name)
        self.assertFalse(os.path.exists(self.p("toolong")))

    def test_scratch_files_of_a_finished_run_are_removed(self):
        self.run_it(stub=self.flags)
        for name in ("inter", "in.gff_db", "lifted.gff3"):
            self.assertFalse(os.path.exists(self.p(name)), name)
        self.assertTrue(os.path.exists(self.p("lifted.gff3.gz")))
        self.assertTrue(os.path.exists(self.p("done.key")))
        self.assertEqual(rd(self.p("types.txt")), "pseudogene\n")

    def test_a_liftoff_that_exits_without_an_output_is_a_stop_even_beside_a_stale_output(self):
        with open(self.p("lifted.gff3"), "w") as fh:
            fh.write(rd(self.lifted))
        with self.assertRaises(LC.StopRun) as cm:
            self.run_it(stub=self.noout)
        self.assertIn("liftoff failed", str(cm.exception))
        self.assertFalse(os.path.exists(self.p("done.key")))

    def test_a_failed_rerun_leaves_no_done_key_of_the_earlier_run(self):
        self.run_it()
        self.assertTrue(LC._is_done(__import__("pathlib").Path(self.d)))
        with self.assertRaises(LC.StopRun):
            self.run_it(stub=T.stub_liftoff(self.tmp, exit_code=3), target="/other.fa")
        self.assertFalse(os.path.exists(self.p("done.key")))
        self.assertFalse(LC._is_done(__import__("pathlib").Path(self.d)))

    def test_a_toolong_run_names_its_records_and_the_marker_is_cleared_by_a_later_success(self):
        res = self.run_it(stub=T.stub_liftoff(self.tmp, sleep=30), hard_s=1)
        self.assertEqual((res["status"], res["records"]), ("toolong", 2))
        self.assertTrue(os.path.exists(self.p("toolong")))
        self.assertIn("killed after 1 s at 2 records", rd(self.p("toolong")))
        self.assertEqual(self.run_it()["status"], "done")
        self.assertFalse(os.path.exists(self.p("toolong")))

    def test_the_post_run_memory_check_uses_the_limit_that_was_passed(self):
        with self.assertRaises(LC.StopRun) as cm:
            self.run_it(rss_limit_gb=0.0005)                     # the stub shell peaks at a few MB
        self.assertIn("peak RSS", str(cm.exception))
        self.assertFalse(os.path.exists(self.p("done.key")))
        self.assertEqual(self.run_it(rss_limit_gb=2.0)["status"], "done")

    def test_the_watchdog_poll_interval_is_the_one_that_was_passed(self):
        long = write_exec(os.path.join(self.tmp, "long.sh"), "#!/bin/bash\nexec sleep 30\n")
        RecordingPopen.waits = []
        with mock.patch.object(LC, "group_rss_kb", return_value=2_000_000), mock.patch.object(LC.subprocess, "Popen", RecordingPopen):
            with self.assertRaises(LC.StopRun) as cm:
                self.run_it(stub=long, hard_s=120, rss_limit_gb=1.0, poll_s=0.05)
        first = [t for t, at in RecordingPopen.waits if t not in (None, 15)][0]
        self.assertEqual(first, 0.05)                           # the default poll is 2 s
        self.assertIn("RSS", str(cm.exception))
        self.assertTrue(os.path.exists(self.p("rss_stop")))
        self.assertIn("2 records", rd(self.p("rss_stop")))


# ===================================================================================================== J. unit packing, run GFFs, plan and identity freezes
class PackingTests(T.HermeticCase):
    def test_the_default_call_budget_is_530_seconds(self):
        units = [("H", "a"), ("H", "b")]
        self.assertEqual(LC.plan_calls(units, {("H", "a"): 300, ("H", "b"): 230}), [units])
        self.assertEqual(LC.plan_calls(units, {("H", "a"): 300, ("H", "b"): 231}), [[units[0]], [units[1]]])

    def test_a_second_call_accumulates_from_zero_and_counts_each_unit_once(self):
        units = [("H", "a"), ("H", "b"), ("H", "c"), ("H", "d")]
        pred = {("H", "a"): 600, ("H", "b"): 300, ("H", "c"): 230, ("H", "d"): 1}
        self.assertEqual(LC.plan_calls(units, pred, budget_s=530), [[units[0]], [units[1], units[2]], [units[3]]])
        pred = {("H", "a"): 200, ("H", "b"): 200, ("H", "c"): 130, ("H", "d"): 1}
        self.assertEqual(LC.plan_calls(units[:3], pred, budget_s=530), [units[:3]])
        self.assertEqual(LC.plan_calls(units, pred, budget_s=530), [units[:3], units[3:]])

    def test_predict_wall_has_the_registered_defaults(self):
        self.assertEqual(LC.predict_wall(100), 25.0 + 100 * 0.09 * 1.25)
        self.assertEqual(LC.predict_wall(300), 25.0 + 33.75)


class RunGffAndFreezeTests(T.HermeticCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.gff = os.path.join(self.tmp, "ref.gff")
        with open(self.gff, "w") as fh:
            fh.write(T.REF_GFF)

    def test_a_changed_assignment_with_the_same_labels_rewrites_the_run_gffs(self):
        out = os.path.join(self.tmp, "gff")
        self.assertTrue(LC.ensure_run_gffs(self.gff, {"gene-A": "r01", "gene-B": "r01"}, out))
        self.assertTrue(LC.ensure_run_gffs(self.gff, {"gene-A": "r01", "gene-P": "r01"}, out))
        text = rd(os.path.join(out, "r01", "in.gff"))
        self.assertIn("gene-P", text)
        self.assertNotIn("gene-B", text)

    def test_a_missing_run_gff_is_written_again(self):
        out = os.path.join(self.tmp, "gff")
        assign = {"gene-A": "r01", "gene-B": "r02"}
        self.assertTrue(LC.ensure_run_gffs(self.gff, assign, out))
        os.remove(os.path.join(out, "r02", "in.gff"))
        self.assertTrue(LC.ensure_run_gffs(self.gff, assign, out))
        self.assertTrue(os.path.exists(os.path.join(out, "r02", "in.gff")))
        self.assertFalse(LC.ensure_run_gffs(self.gff, assign, out))

    def test_ids_are_read_from_the_ID_attribute_only(self):
        with open(self.gff, "w") as fh:
            fh.write("##gff-version 3\nNC_1\tRefSeq\tgene\t101\t1100\t.\t+\t.\tAlt_ID=wrong;Name=A;ID=gene-A\n"
                     "NC_1\tRefSeq\tmRNA\t101\t1100\t.\t+\t.\tID=rna-A1;Alt_Parent=wrong;Parent=gene-A\n"
                     "NC_1\tRefSeq\texon\t101\t500\t.\t+\t.\tID=exon-A1-1;Parent=rna-A1\n")
        counts = LC.write_run_gffs(self.gff, {"gene-A": "r01"}, os.path.join(self.tmp, "gff"))
        self.assertEqual(counts, {"r01": 1})
        self.assertEqual(len(rd(os.path.join(self.tmp, "gff", "r01", "in.gff")).splitlines()), 4)

    def test_a_line_with_neither_id_nor_parent_is_skipped(self):
        with open(self.gff, "a") as fh:
            fh.write("NC_1\tRefSeq\tbiological_region\t1\t10\t.\t+\t.\tNote=no id here\n")
        LC.write_run_gffs(self.gff, {"gene-A": "r01"}, os.path.join(self.tmp, "gff"))
        self.assertNotIn("biological_region", rd(os.path.join(self.tmp, "gff", "r01", "in.gff")))

    def test_children_are_routed_by_their_first_parent(self):
        with open(self.gff, "w") as fh:
            fh.write("##gff-version 3\nNC_1\tRefSeq\tgene\t101\t1100\t.\t+\t.\tID=gene-A\nNC_1\tRefSeq\tgene\t2001\t3000\t.\t+\t.\tID=gene-B\n"
                     "NC_1\tRefSeq\tmRNA\t101\t1100\t.\t+\t.\tID=rna-A\tParent=gene-A\n".replace("\tParent", ";Parent") +
                     "NC_1\tRefSeq\tmRNA\t2001\t3000\t.\t+\t.\tID=rna-B;Parent=gene-B\n"
                     "NC_1\tRefSeq\texon\t101\t500\t.\t+\t.\tID=exon-shared;Parent=rna-A,rna-B\n")
        counts = LC.write_run_gffs(self.gff, {"gene-A": "r01", "gene-B": "r02"}, os.path.join(self.tmp, "gff"))
        self.assertIn("exon-shared", rd(os.path.join(self.tmp, "gff", "r01", "in.gff")))
        self.assertNotIn("exon-shared", rd(os.path.join(self.tmp, "gff", "r02", "in.gff")))
        self.assertEqual(counts, {"r01": 1, "r02": 1})

    def test_the_plan_digest_format(self):
        plan = [dict(label="r02", records=[dict(id="c")]), dict(label="r01", records=[dict(id="b"), dict(id="a")])]
        want = hashlib.sha256(b"nodes" + b"|r01:a,b" + b"|r02:c").hexdigest()
        self.assertEqual(LC.plan_digest(plan, "nodes"), want)
        self.assertEqual(LC.plan_digest([dict(label="r01", records=[dict(id="a"), dict(id="b")]), plan[0]], "nodes"), want)
        self.assertNotEqual(LC.plan_digest(plan, "records"), want)
        self.assertEqual(LC.plan_digest([], "nodes"), hashlib.sha256(b"nodes").hexdigest())

    def test_the_freezes_create_their_directory_and_their_files(self):
        deep = os.path.join(self.tmp, "a", "b")
        LC.freeze_plan(deep, "abc")
        self.assertEqual(rd(os.path.join(deep, "plan.digest")), "abc\n")
        deep2 = os.path.join(self.tmp, "c", "d")
        LC.freeze_identity(deep2, {"f": "1" * 64}, "abcdef1")
        self.assertEqual(json.loads(rd(os.path.join(deep2, "run_identity.json"))), {"amendment_commit": "abcdef1", "code_sha256": {"f": "1" * 64}})
        LC.freeze_identity(deep2, {"f": "1" * 64}, "abcdef1")
        with self.assertRaises(LC.StopRun):
            LC.freeze_identity(deep2, {"f": "2" * 64}, "abcdef1")
        with self.assertRaises(LC.StopRun):
            LC.freeze_identity(deep2, {"f": "1" * 64}, "abcdef2")

    def test_the_build_check_stops_on_another_plan_and_accepts_the_same(self):
        LC.freeze_plan(self.tmp, "abc")
        LC.check_plan_frozen(self.tmp, "abc")
        with self.assertRaises(LC.StopRun):
            LC.check_plan_frozen(self.tmp, "abd")


class CommittedCodeTests(T.HermeticCase):
    """the default CODE_FILES, in a throw-away repository with the same layout"""

    def setUp(self):
        self.repo = tempfile.mkdtemp()
        T.git(self.repo, "init", "-q")
        for f in LC.CODE_FILES + ("docs/AMENDMENT.md",):
            os.makedirs(os.path.dirname(os.path.join(self.repo, f)), exist_ok=True)
            with open(os.path.join(self.repo, f), "w") as fh:
                fh.write(f"# {f}\n")
        T.git(self.repo, "add", "-A")
        T.git(self.repo, "commit", "-q", "-m", "amendment and code")
        self.amend = T.git(self.repo, "rev-parse", "HEAD")

    def edit(self, f):
        with open(os.path.join(self.repo, f), "a") as fh:
            fh.write("# edit\n")

    def test_clean_code_is_accepted_with_the_default_files_and_the_head_is_returned(self):
        self.assertEqual(LC.check_committed_code(self.amend, self.repo), self.amend)
        self.edit("docs/AMENDMENT.md")                                      # a file that is not code
        T.git(self.repo, "commit", "-q", "-am", "later")
        head = T.git(self.repo, "rev-parse", "HEAD")
        self.assertNotEqual(head, self.amend)
        self.assertEqual(LC.check_committed_code(self.amend, self.repo), head)

    def test_each_default_code_file_is_checked(self):
        for f in LC.CODE_FILES:
            self.edit(f)
            with self.assertRaises(LC.StopRun, msg=f) as cm:
                LC.check_committed_code(self.amend, self.repo)
            self.assertIn(f, str(cm.exception))
            T.git(self.repo, "checkout", "-q", "--", f)
        os.remove(os.path.join(self.repo, LC.CODE_FILES[2]))
        with self.assertRaises(LC.StopRun):
            LC.check_committed_code(self.amend, self.repo)

    def test_a_staged_but_uncommitted_change_is_refused(self):
        self.edit(LC.CODE_FILES[0])
        T.git(self.repo, "add", LC.CODE_FILES[0])
        with self.assertRaises(LC.StopRun) as cm:
            LC.check_committed_code(self.amend, self.repo)
        self.assertIn("differs from HEAD", str(cm.exception))

    def test_a_modified_file_that_is_not_code_does_not_stop_the_run(self):
        self.edit("docs/AMENDMENT.md")
        self.assertEqual(LC.check_committed_code(self.amend, self.repo), self.amend)

    def test_the_code_identity_covers_the_default_files_with_full_hashes(self):
        ident = LC.code_identity(self.repo)
        self.assertEqual(sorted(ident), sorted(LC.CODE_FILES))
        for f, h in ident.items():
            self.assertEqual(h, hashlib.sha256(open(os.path.join(self.repo, f), "rb").read()).hexdigest())
            self.assertEqual(len(h), 64)



# ===================================================================================================== K. the budgeted driver
class DriverBase(T.HermeticCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.gff = os.path.join(self.tmp, "ref.gff")
        with open(self.gff, "w") as fh:
            fh.write(T.REF_GFF)
        self.e_root = os.path.join(self.tmp, "e")
        gffdir = os.path.join(self.e_root, "gff")
        LC.write_run_gffs(self.gff, {"gene-A": "r01", "gene-B": "r01", "gene-P": "r02", "gene-Z": "ctl1"}, gffdir)
        LC.write_run_gffs(self.gff, {"gene-A": "r04", "gene-B": "r04", "gene-P": "r04", "gene-Z": "r04"}, gffdir)
        mk = lambda i, c, s0, e: dict(id=i, component="PC", contig=c, start0=s0, end=e, strand="+", idx=0, rtype="gene", biotype="x", exonic_len=e - s0)  # noqa: E731
        a, b, p_, z = mk("gene-A", "NC_1", 100, 1100), mk("gene-B", "NC_1", 2000, 3000), mk("gene-P", "NC_1", 5000, 6000), mk("gene-Z", "NC_1", 8000, 9000)
        self.ids4 = ["gene-A", "gene-B", "gene-P", "gene-Z"]
        self.plan = {"r01": dict(run=1, label="r01", records=[a, b], shards=[("r01", ["gene-A", "gene-B"])]),
                     "r02": dict(run=2, label="r02", records=[p_], shards=[("r02", ["gene-P"])]),
                     "ctl1": dict(run=1, label="ctl1", records=[z], shards=[("ctl1", ["gene-Z"])]),
                     "r04": dict(run=4, label="r04", records=[a, b, p_, z], shards=[("r04", list(self.ids4))])}
        self.lifted = os.path.join(self.tmp, "stub_lifted.gff3")
        with open(self.lifted, "w") as fh:
            fh.write(T.lifted_gff([("gene-A", "h1", 1001, 2000, "+", "0.98", "0.97", 0, "False", "False", [(1001, 2000)], "gene")]))
        self.kw = dict(e_root=self.e_root, outgroups={"H": dict(species="human", genome="/t/h.fa", mmi="/t/h.mmi")}, reference="/r/g.fa",
                       liftoff_bin=T.stub_liftoff_counting(self.tmp), minimap2="/m/mm2", hard_s=60, env={"STUB_LIFTED": self.lifted}, ensure_links=False)

    def run_dir(self, label, species="human"):
        return os.path.join(self.e_root, "lift", species, "runs", label)

    def done(self, label):
        return os.path.exists(os.path.join(self.run_dir(label), "done.key"))


class DriverRuleTests(DriverBase):
    def slow(self, n="1", **extra):
        return dict(self.kw, env={"STUB_LIFTED": self.lifted, "STUB_MAX_GENES": n}, **extra)

    def test_a_unit_stopped_by_the_end_of_the_call_runs_again_in_the_next_call(self):
        units = [("H", "r02"), ("H", "r01")]
        res = LC.run_budgeted(units, self.plan, budget_s=500, call_limit_s=2, **self.slow(hard_s=60))
        self.assertEqual((res["exit_code"], res["deadline"], res["split"]), (75, ("H", "r01"), None))
        self.assertEqual(res["pending"], [("H", "r01")])
        res = LC.run_budgeted(units, self.plan, budget_s=500, **self.kw)
        self.assertEqual((res["exit_code"], res["done"]), (0, units))
        self.assertTrue(self.done("r01") and self.done("r02"))

    def test_the_first_unit_stopped_by_a_short_call_limit_is_too_long_for_a_call_and_is_split(self):
        res = LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, call_limit_s=2, **self.slow(hard_s=60))
        self.assertEqual((res["exit_code"], res["split"], res["deadline"]), (75, ("H", "r01"), None))
        self.assertEqual(LC.load_splits(self.e_root, "human"), ["r01"])

    def test_a_later_unit_stopped_at_exactly_the_hard_limit_is_split_not_retried(self):
        res = LC.run_budgeted([("H", "r02"), ("H", "r01")], self.plan, budget_s=500, **self.slow(hard_s=1))
        self.assertEqual((res["exit_code"], res["split"], res["deadline"]), (75, ("H", "r01"), None))
        self.assertEqual((res["done"], res["pending"]), ([("H", "r02")], [("H", "r01")]))
        self.assertTrue(os.path.exists(os.path.join(self.run_dir("r01"), "toolong")))

    def test_splits_accumulate_and_a_half_is_recorded_under_its_own_label(self):
        kw = self.slow(hard_s=1)
        seen = []
        for _ in range(8):
            res = LC.run_budgeted([("H", "r04")], self.plan, budget_s=500, **kw)
            seen.append((res["exit_code"], res["split"]))
            if res["exit_code"] == 0:
                break
        self.assertEqual(seen, [(75, ("H", "r04")), (75, ("H", "r04")), (75, ("H", "r04")), (0, None)])
        self.assertEqual(LC.load_splits(self.e_root, "human"), ["r04", "r04.s1", "r04.s2"])
        for lab in ("r04.s1.s1", "r04.s1.s2", "r04.s2.s1", "r04.s2.s2"):
            self.assertTrue(self.done(lab), lab)
        self.assertFalse(self.done("r04"))
        with open(os.path.join(self.run_dir("r04.s1.s1"), "in.gff")) as fh:
            self.assertEqual(sum(1 for l in fh if "\tgene\t" in l or "\tpseudogene\t" in l), 1)

    def test_the_directory_of_a_unit_follows_the_species_of_the_outgroup(self):
        kw = dict(self.kw, outgroups={"H": dict(species="zz", genome="/t/h.fa", mmi="/t/h.mmi")})
        res = LC.run_budgeted([("H", "r02")], self.plan, budget_s=500, **kw)
        self.assertEqual(res["exit_code"], 0)
        self.assertTrue(os.path.exists(os.path.join(self.e_root, "lift", "zz", "runs", "r02", "done.key")))
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "lift", "human")))
        self.assertTrue(str(LC.run_dir_of("H", "r02", self.e_root, kw["outgroups"])).endswith("/lift/zz/runs/r02"))

    def test_a_run_gff_that_is_not_there_is_a_stop(self):
        plan = dict(self.plan, r09=dict(run=9, label="r09", records=self.plan["r02"]["records"], shards=[("r09", ["gene-P"])]))
        with self.assertRaises(LC.StopRun) as cm:
            LC.run_budgeted([("H", "r09")], plan, budget_s=500, **self.kw)
        self.assertIn("write the run GFFs first", str(cm.exception))

    def test_a_changed_run_gff_replaces_the_in_gff_of_an_unfinished_run(self):
        seen = write_exec(os.path.join(self.tmp, "seen.sh"), "#!/bin/bash\ncp in.gff seen.gff\nexit 3\n")
        with self.assertRaises(LC.StopRun):
            LC.run_budgeted([("H", "r02")], self.plan, budget_s=500, **dict(self.kw, liftoff_bin=seen))
        old = rd(os.path.join(self.run_dir("r02"), "seen.gff"))
        with open(os.path.join(self.e_root, "gff", "r02", "in.gff"), "a") as fh:
            fh.write("NC_1\tRefSeq\tgene\t20001\t21000\t.\t+\t.\tID=gene-Q;Name=Q\n")
        with self.assertRaises(LC.StopRun):
            LC.run_budgeted([("H", "r02")], self.plan, budget_s=500, **dict(self.kw, liftoff_bin=seen))
        new = rd(os.path.join(self.run_dir("r02"), "seen.gff"))
        self.assertNotIn("gene-Q", old)
        self.assertIn("gene-Q", new)
        self.assertEqual(new, rd(os.path.join(self.e_root, "gff", "r02", "in.gff")))

    def test_effective_shards_put_the_halves_where_the_shard_was_and_split_again_when_told(self):
        recs = shard_recs(4)
        by_id = {r["id"]: r for r in recs}
        shards = [("A", ["g000", "g001"]), ("B", ["g002", "g003"])]
        self.assertEqual(LC.effective_shards(shards, by_id, []), shards)
        self.assertEqual(LC.effective_shards(shards, by_id, ["A"]), [("A.s1", ["g000"]), ("A.s2", ["g001"]), ("B", ["g002", "g003"])])
        self.assertEqual(LC.effective_shards([("A", [r["id"] for r in recs])], by_id, ["A", "A.s2"]),
                         [("A.s1", ["g000", "g001"]), ("A.s2.s1", ["g002"]), ("A.s2.s2", ["g003"])])
        self.assertEqual(LC.shard_labels([dict(shards=shards, records=recs)], self.e_root, {"H": dict(species="human")}, "H"), ["A", "B"])

    def test_a_recorded_split_of_a_shard_without_a_cut_is_a_stop(self):
        recs = shard_recs(2)
        by_id = {r["id"]: r for r in recs}
        with self.assertRaises(LC.StopRun):
            LC.effective_shards([("A", ["g000"])], by_id, ["A"])

    def test_splits_are_stored_sorted_and_unique_per_outgroup(self):
        deep = os.path.join(self.tmp, "x", "y")
        LC.save_splits(deep, "human", ["b", "a", "a"])
        self.assertEqual(json.loads(rd(os.path.join(deep, "lift", "human", "splits.json"))), ["a", "b"])
        self.assertEqual(LC.load_splits(deep, "human"), ["a", "b"])
        self.assertEqual(LC.load_splits(deep, "chimpanzee"), [])

    def test_a_shard_is_done_only_with_both_its_done_key_and_its_loci(self):
        d = __import__("pathlib").Path(self.tmp) / "shard"
        d.mkdir()
        self.assertFalse(LC._is_done(d))
        (d / "done.key").write_text("{}")
        self.assertFalse(LC._is_done(d))
        (d / "loci.tsv").write_text("x")
        self.assertTrue(LC._is_done(d))
        (d / "done.key").unlink()
        self.assertFalse(LC._is_done(d))


class DriverClockTests(DriverBase):
    """run_budgeted with a fake clock and a fake run_shard: the numbers of the call limit and the packing"""

    def world(self, durations):
        self.now, self.seen, durations = [1000.0], [], list(durations)

        def fake(d, ids, **kw):
            self.seen.append(dict(kw, label=d.name, n=len(ids)))
            self.now[0] += durations.pop(0)
            (d / "done.key").write_text("{}")
            (d / "loci.tsv").write_text("x")
            return dict(status="done", wall_s=1.0, peak_rss_kb=1)
        return mock.patch.object(LC, "run_shard", fake)

    def go(self, units, plan=None, **extra):
        kw = {k: v for k, v in self.kw.items() if k not in ("hard_s", "liftoff_bin")}
        kw.update(liftoff_bin="/l", clock=lambda: self.now[0])
        kw.update(extra)
        return LC.run_budgeted(units, plan or self.plan, **kw)

    def test_the_time_left_to_a_run_is_the_call_limit_less_what_the_call_has_used(self):
        with self.world([100, 479.5, 0]):
            res = self.go([("H", "r01"), ("H", "r02"), ("H", "ctl1")], budget_s=1e9)
        self.assertEqual(res["exit_code"], 0)
        self.assertEqual([k["hard_s"] for k in self.seen], [555.0, 455.0, 1.0])
        self.assertEqual({k["rss_limit_gb"] for k in self.seen}, {16.0})
        self.assertEqual([k["label"] for k in self.seen], ["r01", "r02", "ctl1"])

    def test_without_a_call_limit_a_run_gets_the_hard_limit_of_590_seconds(self):
        with self.world([0]):
            self.go([("H", "r01")], budget_s=1e9, call_limit_s=10_000)
        self.assertEqual(self.seen[0]["hard_s"], 590)

    def test_a_unit_that_exactly_fills_the_budget_starts(self):
        with self.world([100, 100]):
            res = self.go([("H", "r01"), ("H", "r02")], budget_s=200, startup_s=100, per_record_s=0.0)
        self.assertEqual((res["exit_code"], len(res["done"])), (0, 2))

    def test_a_unit_that_does_not_fit_the_budget_by_a_second_waits_for_the_next_call(self):
        with self.world([100, 100]):
            res = self.go([("H", "r01"), ("H", "r02")], budget_s=199, startup_s=100, per_record_s=0.0)
        self.assertEqual((res["exit_code"], res["done"], res["pending"]), (75, [("H", "r01")], [("H", "r02")]))

    def test_the_first_unit_of_a_call_starts_even_when_it_alone_is_over_the_budget(self):
        with self.world([100]):
            res = self.go([("H", "r01"), ("H", "r02")], budget_s=50, startup_s=100, per_record_s=0.0)
        self.assertEqual((res["exit_code"], res["done"], res["pending"]), (75, [("H", "r01")], [("H", "r02")]))

    def test_the_prediction_of_a_unit_sums_its_unfinished_shards_only(self):
        plan = dict(self.plan)
        plan["r04"] = dict(plan["r04"], shards=[("r04.b1", ["gene-A", "gene-B"]), ("r04.b2", ["gene-P", "gene-Z"])])
        # one shard finished: only the other one is predicted (100 s), 100 + 100 fits 200
        d = os.path.join(self.run_dir("r04.b1"))
        os.makedirs(d)
        # a finished shard is one made by THIS command from THIS in.gff: its done.key holds the key of the call (the command go() gives run_budgeted; a blocked shard's
        # in.gff is its own file, absent here)
        key = LC._key(LC.liftoff_command("/t/h.fa", "/r/g.fa", "/m/mm2", "/l"), __import__("pathlib").Path(d) / "in.gff")
        for n, text in (("done.key", json.dumps(key)), ("loci.tsv", "x")):
            with open(os.path.join(d, n), "w") as fh:
                fh.write(text)
        with self.world([100, 100]):
            res = self.go([("H", "r01"), ("H", "r04")], plan, budget_s=200, startup_s=100, per_record_s=0.0)
        self.assertEqual((res["exit_code"], [k["label"] for k in self.seen]), (0, ["r01", "r04.b2"]))
        # nothing finished: two start-ups (200 s) are predicted, 100 + 200 does not fit 250
        os.remove(os.path.join(d, "done.key"))
        os.remove(os.path.join(d, "loci.tsv"))
        shutil_rm = [os.path.join(self.run_dir(x), n) for x in ("r01", "r04.b2") for n in ("done.key", "loci.tsv")]
        for f in shutil_rm:
            os.remove(f)
        with self.world([100, 100, 100]):
            res = self.go([("H", "r01"), ("H", "r04")], plan, budget_s=250, startup_s=100, per_record_s=0.0)
        self.assertEqual((res["exit_code"], res["done"], res["pending"]), (75, [("H", "r01")], [("H", "r04")]))
        self.assertEqual([k["label"] for k in self.seen], ["r01"])

    def test_a_unit_predicts_its_records_at_the_per_record_rate_with_the_margin(self):
        d1, p2 = LC.predict_wall(2), LC.predict_wall(1)                       # r01 holds two records, r02 one
        self.assertEqual((d1, p2), (25.0 + 2 * 0.09 * 1.25, 25.0 + 1 * 0.09 * 1.25))
        with self.world([d1, 0]):
            res = self.go([("H", "r01"), ("H", "r02")], budget_s=d1 + p2 + 0.001)
        self.assertEqual(res["exit_code"], 0)

    def test_a_unit_whose_prediction_is_a_hair_over_the_budget_waits(self):
        d1, p2 = LC.predict_wall(2), LC.predict_wall(1)
        with self.world([d1, 0]):
            res = self.go([("H", "r01"), ("H", "r02")], budget_s=d1 + p2 - 0.001)
        self.assertEqual((res["exit_code"], res["pending"]), (75, [("H", "r02")]))

    def test_the_target_is_the_genome_link_and_the_three_links_point_at_the_inputs(self):
        T.make_fake_genome(os.path.join(self.tmp, "g.fa"))                # a FASTA, its .fai and a k15 w10 index of it: the links are made to files that pass the preflight
        og = {"H": dict(species="human", genome=os.path.join(self.tmp, "g.fa"), mmi=os.path.join(self.tmp, "g.fa.mmi"))}
        with mock.patch.dict(LC.OUTGROUPS, og, clear=True), self.world([0]):
            self.go([("H", "r01")], budget_s=1e9, ensure_links=True, outgroups=og)
        d = os.path.join(self.e_root, "lift", "human")
        self.assertEqual(self.seen[0]["target"], os.path.join(d, "genome.fa"))
        self.assertEqual(os.readlink(os.path.join(d, "genome.fa")), og["H"]["genome"])
        self.assertEqual(os.readlink(os.path.join(d, "genome.fa.fai")), og["H"]["genome"] + ".fai")
        self.assertEqual(os.readlink(os.path.join(d, "genome.fa.mmi")), og["H"]["mmi"])
        self.assertTrue(os.path.samefile(os.path.join(d, "genome.fa.mmi"), og["H"]["mmi"]))


# ===================================================================================================== L. Gate 7 from finished runs
class GateRuleTests(T.WorldCase):
    def gate(self, group):
        return LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)(group)

    def test_a_failed_control_stops_every_later_call_not_only_the_first(self):
        self.finish_controls(neg_ok=False)
        for _ in range(3):
            with self.assertRaises(LC.StopRun) as cm:
                self.gate("verdict")
            self.assertIn("negative", str(cm.exception))
        self.finish_controls(neg_ok=True)
        self.place("pos1", {"gene-REF": None, "gene-F1": 2})
        for _ in range(2):
            with self.assertRaises(LC.StopRun):
                self.gate("neg")

    def test_the_saved_verdict_records_a_failure_as_a_failure(self):
        self.finish_controls(neg_ok=False)
        with self.assertRaises(LC.StopRun):
            self.gate("verdict")
        with open(os.path.join(self.e_root, "gate7.json")) as fh:
            g = json.load(fh)
        self.assertIs(g["valid"], False)
        self.assertIs(g["positive"]["pass"], True)
        self.assertIs(g["negative"]["pass"], False)
        self.assertEqual(sorted(g["control_runs"]), ["neg", "pos"])

    def test_the_neg_gate_records_no_overall_verdict(self):
        self.finish_controls()
        self.gate("neg")
        with open(os.path.join(self.e_root, "gate7.json")) as fh:
            g = json.load(fh)
        self.assertIsNone(g["valid"])
        self.assertEqual(sorted(g["control_runs"]), ["pos"])
        self.assertNotIn("negative", g)

    def test_a_negative_control_run_that_changes_after_a_pass_is_judged_again(self):
        self.finish_controls()
        self.gate("verdict")
        self.place("ctl2", {r["id"]: 1 for r in self.w["ctl"][1]["records"]})             # every planted pair now sits on one locus
        with self.assertRaises(LC.StopRun) as cm:
            self.gate("verdict")
        self.assertIn("negative", str(cm.exception))

    def test_the_digest_of_a_control_depends_on_its_group_runs_outgroups_and_every_file(self):
        self.finish_controls()
        w, e, og = self.w, self.e_root, self.og
        d_pos, d_neg = LC.control_digest(w, e, og, "pos"), LC.control_digest(w, e, og, "neg")
        self.assertNotEqual(d_pos, d_neg)
        self.assertEqual(d_pos, LC.control_digest(w, e, og, "pos"))
        base = os.path.join(e, "lift", "human", "runs", "pos1")
        for name in ("loci.tsv", "unmapped.txt", "done.key"):
            path = os.path.join(base, name)
            keep = open(path, "rb").read()
            with open(path, "ab") as fh:
                fh.write(b"\n")
            self.assertNotEqual(LC.control_digest(w, e, og, "pos"), d_pos, name)
            with open(path, "wb") as fh:
                fh.write(keep)
            self.assertEqual(LC.control_digest(w, e, og, "pos"), d_pos, name)
        # a missing file is not an empty file
        path = os.path.join(base, "unmapped.txt")
        self.assertEqual(os.path.getsize(path), 0)
        os.remove(path)
        self.assertNotEqual(LC.control_digest(w, e, og, "pos"), d_pos)
        open(path, "w").close()
        # the same files in another outgroup are another digest
        a, b = os.path.join(e, "lift", "human", "runs", "pos1"), os.path.join(e, "lift", "chimpanzee", "runs", "pos1")
        tmp = os.path.join(self.tmp, "swap")
        os.rename(a, tmp)
        os.rename(b, a)
        os.rename(tmp, b)
        self.assertNotEqual(LC.control_digest(w, e, og, "pos"), d_pos)

    def test_the_digest_of_the_negative_group_follows_the_planted_runs_only(self):
        self.finish_controls()
        w, e, og = self.w, self.e_root, self.og
        d_neg = LC.control_digest(w, e, og, "neg")
        with open(os.path.join(e, "lift", "human", "runs", "pos1", "loci.tsv"), "a") as fh:
            fh.write("\n")
        self.assertEqual(LC.control_digest(w, e, og, "neg"), d_neg)
        with open(os.path.join(e, "lift", "human", "runs", "ctl1", "loci.tsv"), "a") as fh:
            fh.write("\n")
        self.assertNotEqual(LC.control_digest(w, e, og, "neg"), d_neg)

    def test_the_positive_control_is_judged_on_the_whole_triple(self):
        self.finish_controls()
        self.place("pos3", {"gene-X1": 3, "gene-F3": 2}, outgroups="C")            # in the chimpanzee gene-X1 sits on another locus
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        self.assertFalse(res["positive"]["pass"])
        self.assertEqual([m["record"] for m in res["positive"]["missing"]], ["gene-X1"])

    def test_the_negative_control_is_judged_on_the_whole_triple(self):
        self.finish_controls()
        self.place("ctl2", {r["id"]: 1 for r in self.w["ctl"][1]["records"]}, outgroups="H")       # the human loci agree, the chimpanzee and orangutan ones do not
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        self.assertEqual((res["negative"]["called"], res["negative"]["pass"]), (0, True))
        self.assertTrue(res["valid"])

    def test_the_positive_control_reports_its_triple(self):
        self.finish_controls()
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        self.assertEqual(res["positive"]["triple"], ["gene-H1", "gene-C1", "gene-O1"])
        self.assertEqual([e["id"] for e in res["_control_expansions"]], ["EC0001", "EC0002"])
        self.assertEqual(res["negative"]["called_pairs"], [])

    def test_the_failure_text_names_only_the_control_that_failed(self):
        both = dict(positive=dict(**{"pass": False}, missing=[dict(record="gene-X1", reason="no_placement")]), negative=dict(called=2, of=50, bar=1, **{"pass": False}))
        text = LC.failure_text(both)
        self.assertIn("positive control", text)
        self.assertIn("negative control: 2 of 50", text)
        only_neg = dict(positive=dict(**{"pass": True}, missing=[]), negative=dict(called=2, of=50, bar=1, **{"pass": False}))
        self.assertNotIn("positive", LC.failure_text(only_neg))
        only_pos = dict(positive=dict(**{"pass": False}, missing=[]), negative=dict(called=0, of=50, bar=1, **{"pass": True}))
        self.assertNotIn("negative", LC.failure_text(only_pos))
        self.assertEqual(LC.failure_text(dict(positive=dict(**{"pass": True}, missing=[]))), "")
        self.assertEqual(LC.failure_text({}), "")


# ===================================================================================================== M. build: tables, status, exit codes
def rename_contigs(rows, recs, mapping):
    for r in rows:
        for k in ("src_contig", "contig"):
            r[k] = mapping.get(r[k], r[k])
    for r in recs.values():
        r["contig"] = mapping.get(r["contig"], r["contig"])


class HubWorld(T.WorldCase):
    """WorldCase plus one component of 4 long and 45 short annotated records and 5 unannotated sites: 54 nodes (a hub) but 49 records (not a hub)"""

    def setUp(self):
        super().setUp()
        rows, recs = T.HubAndCalledSetTests().big_pool(n_short=45, n_sites=5)
        rename_contigs(rows, recs, {"c1": "hb1", "c2": "hb2", "c9": "hb9"})
        self.rows, self.recs = self.rows + rows, dict(self.recs, **recs)
        T.write_selflift(self.root, self.rows, self.recs)
        self.w = LC.plan_from_real_inputs("nodes", lift_root=self.root, spec=self.spec, n_neg=4)


class HubReadingTests(HubWorld):
    def test_the_hub_reading_decides_the_called_set_of_the_plan(self):
        nodes = LC.plan_from_real_inputs("nodes", lift_root=self.root, spec=self.spec, n_neg=4)
        records = LC.plan_from_real_inputs("records", lift_root=self.root, spec=self.spec, n_neg=4)
        self.assertEqual((len(nodes["called"]), len(records["called"])), (4, 5))
        self.assertEqual(records["hub_mode"], "records")
        self.assertEqual(sum(len(r["records"]) for r in records["plan"]) - sum(len(r["records"]) for r in nodes["plan"]), 4)
        self.assertNotEqual(LC.plan_digest(nodes["plan"] + nodes["pos"] + nodes["ctl"], "nodes"), LC.plan_digest(records["plan"] + records["pos"] + records["ctl"], "records"))

    def test_the_build_uses_the_hub_reading_it_was_given_and_the_frozen_plan_of_that_reading(self):
        self.finish_controls()
        self.finish_verdict()
        LC.freeze_plan(self.e_root, LC.plan_digest(self.w["plan"] + self.w["pos"] + self.w["ctl"], "nodes"))
        res = LC.build_all(self.e_root, os.path.join(self.tmp, "out"), lift_root=self.root, outgroups=self.og, hub_mode="nodes", spec=self.spec, n_neg=4)
        self.assertTrue(res["gate7"]["valid"])
        with open(os.path.join(self.tmp, "out", "components.tsv")) as fh:
            by = {r[0]: r for r in list(csv.reader(fh, delimiter="\t"))[1:]}
        hub = next(c for c in by if by[c][1] == "54")
        self.assertEqual((by[hub][5], by[hub][4]), ("hub", "4"))
        self.assertEqual(by["PC_gene-V1"][5], "verdict")
        self.assertEqual(by["PC_gene-N"][5], "control")
        self.assertEqual(by["PC_gene-F1"][5], "control")
        with self.assertRaises(LC.StopRun):                                                  # the other reading is another plan
            LC.build_all(self.e_root, os.path.join(self.tmp, "out2"), lift_root=self.root, outgroups=self.og, hub_mode="records", spec=self.spec, n_neg=4)


class BuildRuleTests(T.WorldCase):
    def build(self, **kw):
        return LC.build_all(self.e_root, os.path.join(self.tmp, "out"), lift_root=self.root, outgroups=self.og, hub_mode="nodes", spec=self.spec, n_neg=4, **kw)

    def table(self, name):
        with open(os.path.join(self.tmp, "out", name)) as fh:
            return list(csv.reader(fh, delimiter="\t"))

    def ready(self):
        self.finish_controls()
        self.finish_verdict()

    def test_a_build_from_runs_of_another_plan_is_a_stop(self):
        self.ready()
        LC.freeze_plan(self.e_root, "0" * 64)
        with self.assertRaises(LC.StopRun):
            self.build()

    def test_the_summary_is_about_the_verdict_set_only(self):
        self.ready()
        res = self.build()
        with open(os.path.join(self.tmp, "out", "summary.json")) as fh:
            s = json.load(fh)
        self.assertEqual((s["called_components"], s["verdict_components"], s["verdict_records"]), (2, 2, 4))
        self.assertEqual((s["expansions"], s["records_in_expansions"], res["expansions"]), (1, 2, 1))
        self.assertEqual(s["records_without_triple_by_reason"], {})
        self.assertEqual(sorted(s["control_components"]), ["PC_gene-F1", "PC_gene-N"])
        self.assertEqual((s["hub_mode"], s["strict_single_copy"], s["gate7_valid"]), ("nodes", False, True))

    def test_the_verdict_set_has_its_own_status_and_the_controls_theirs(self):
        self.ready()
        self.build()
        by = {r[0]: r for r in self.table("components.tsv")[1:]}
        self.assertEqual({c: by[c][5] for c in by}, {"PC_gene-N": "control", "PC_gene-F1": "control", "PC_gene-V1": "verdict", "PC_gene-W1": "verdict"})

    def test_gate7_json_carries_the_digests_of_the_control_outputs_and_the_planted_count(self):
        self.ready()
        self.build()
        with open(os.path.join(self.tmp, "out", "gate7.json")) as fh:
            g = json.load(fh)
        self.assertEqual(g["control_runs"], {"pos": LC.control_digest(self.w, self.e_root, self.og, "pos"), "neg": LC.control_digest(self.w, self.e_root, self.og, "neg")})
        self.assertTrue(all(len(h) == 64 for h in g["control_runs"].values()))
        self.assertEqual((g["valid"], g["planted_components"], g["not_required"]), (True, 2, ["gene-X2"]))
        self.assertEqual(sorted(g["qualifying"]), ["gene-REF", "gene-X1", "gene-X3"])

    def test_sha256sums_list_the_files_in_name_order_whatever_the_directory_order(self):
        self.ready()
        real = os.listdir
        with mock.patch("os.listdir", side_effect=lambda p=".": sorted(real(p), reverse=True)):
            self.build()
        names = [line.partition("  ")[2] for line in rd(os.path.join(self.tmp, "out", "SHA256SUMS")).splitlines()]
        self.assertEqual(names, sorted(names))
        self.assertGreater(len(names), 10)

    def test_the_strict_reading_is_recorded_and_applied(self):
        self.ready()
        rows = T.outgroup_rows(*[(f"gene-H{n}", f"h{n}", 1000, 3000, []) for n in (1, 2, 3, 4)])
        rows.append(row("gene-Hsrc", contig="h8", start0=0, end=500))
        rows.append(row("gene-Hsrc", cls="dropped_M1", contig="h9", start0=0, end=500, note="overlaps gene-H1 (h1)"))
        rows[-1]["tag"] = "1"
        T.write_table(self.og["H"]["selflift"], rows)
        loose = self.build()
        self.assertTrue(loose["gate7"]["valid"])
        with open(os.path.join(self.tmp, "out", "summary.json")) as fh:
            self.assertIs(json.load(fh)["strict_single_copy"], False)
        with self.assertRaises(LC.StopRun):                              # never written over the table of the other reading
            self.build(strict=True)
        strict = LC.build_all(self.e_root, os.path.join(self.tmp, "out_strict"), lift_root=self.root, outgroups=self.og, hub_mode="nodes", spec=self.spec, n_neg=4, strict=True)
        self.assertFalse(strict["gate7"]["valid"])
        self.assertTrue(strict["status"].startswith("INVALID: positive control"))
        with open(os.path.join(self.tmp, "out_strict", "summary.json")) as fh:
            self.assertIs(json.load(fh)["strict_single_copy"], True)
        self.assertEqual(rd(os.path.join(self.tmp, "out_strict", "STATUS")).strip(), strict["status"])

    def test_status_names_the_control_that_failed(self):
        self.finish_controls(neg_ok=False)
        self.finish_verdict()
        res = self.build()
        status = rd(os.path.join(self.tmp, "out", "STATUS"))
        self.assertTrue(status.startswith("INVALID: "))
        self.assertIn("negative control: 2 of 2", status)
        self.assertNotIn("positive control", status)
        self.assertEqual(res["status"] + "\n", status)


# ===================================================================================================== N. the command line
STUB_LOGGING = T.STUB_CLI_LIFTOFF.replace("import json, os, re, sys\n", "import json, os, re, sys\nopen('argv.json', 'w').write(json.dumps(sys.argv[1:]))\n", 1)


class CliBase(T.WorldCase):
    """main() on the fixture world with a stub Liftoff (as the builder's end-to-end tests do), plus the pieces the builder's tests leave out"""

    def setUp(self):
        super().setUp()
        self.fasta = os.path.join(self.tmp, "gorilla", "gorilla.fa")
        T.make_fake_genome(self.fasta, mmi=False)                          # the reference must exist with a .fai that is not older (preflight)
        self.gff_path = os.path.join(self.tmp, "gorilla.gff")
        with open(self.gff_path, "w") as fh:
            fh.write("##gff-version 3\n")
            for rid, r in self.recs.items():
                fh.write(f"{r['contig']}\tRefSeq\tgene\t{r['start0'] + 1}\t{r['end']}\t.\t+\t.\tID={rid};Name={rid}\n")
        self.stub = write_exec(os.path.join(self.tmp, "stub_liftoff.py"), STUB_LOGGING)
        self.map_path = os.path.join(self.tmp, "stub_map.json")
        frag = {"gene-F1", "gene-F2", "gene-F3"}
        locus = {"gene-V1": 3, "gene-V2": 3, "gene-W1": 4, "gene-W2": 2}
        for run in self.w["pos"]:
            for r in run["records"]:
                locus[r["id"]] = 2 if r["id"] in frag else 1
        for j, run in enumerate(self.w["ctl"]):
            for r in run["records"]:
                locus[r["id"]] = 1 if j == 0 else 2
        with open(self.map_path, "w") as fh:
            json.dump({sp: {rid: [f"{loc}{n}", 1000, 3000] for rid, n in locus.items()} for _, sp, loc in T.OG_SPECIES}, fh)
        for pt in (mock.patch.dict(LC.OUTGROUPS, {k: dict(t, assembly="x") for k, t in self.og.items()}, clear=True), mock.patch.object(X, "LRPAP1_CONTROL", self.spec),
                   mock.patch.dict(os.environ, {"STUB_MAP": self.map_path})):
            pt.start()
            self.addCleanup(pt.stop)

    def cli(self, cmd, *extra, plan_patch=True, rehearsal=True):
        argv = [cmd, "--out", self.e_root, "--lift-root", self.root, "--n-planted", "4", *extra]
        if cmd == "run":
            argv += ["--amendment-commit", "abc1234", "--no-code-check", "--gorilla-gff", self.gff_path, "--gorilla-fasta", self.fasta, "--liftoff-bin", self.stub,
                     "--minimap2", "/m/mm2"]
        if cmd == "build" and rehearsal:                                   # the runs were made by a stub: build refuses them unless it is told this is a rehearsal
            argv.append("--rehearsal")
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = LC.main(argv)
        return rc, out.getvalue(), err.getvalue()

    def to_the_end(self, *extra):
        rcs = []
        for _ in range(8):
            rc, out, err = self.cli("run", *extra)
            rcs.append(rc)
            if rc != 75:
                break
        return rcs, err


class CliRunTests(CliBase):
    def test_liftoff_gets_the_registered_command_with_the_genome_link_the_gorilla_fasta_and_the_given_minimap2(self):
        rc, out, err = self.cli("run")
        self.assertEqual(rc, 75, err)
        argv = json.loads(rd(os.path.join(self.e_root, "lift", "human", "runs", "pos1", "argv.json")))
        self.assertEqual(argv, [os.path.join(self.e_root, "lift", "human", "genome.fa"), self.fasta, "-g", "in.gff", "-f", "types.txt", "-o", "lifted.gff3", "-u", "unmapped.txt",
                                "-dir", "inter", "-copies", "-sc", "0.95", "-a", "0.5", "-s", "0.5", "-overlap", "0.1", "-p", "4", "-m", "/m/mm2"])
        log = rd(os.path.join(self.e_root, "run.log"))
        self.assertIn("gate(pos) passed", log)
        self.assertIn("H pos1 [pos]: done", log)
        d = os.path.join(self.e_root, "lift", "human")
        self.assertEqual(sorted(os.listdir(d))[:4], ["genome.fa", "genome.fa.fai", "genome.fa.mmi", "runs"])

    def test_the_budget_of_the_call_is_the_one_on_the_command_line(self):
        self.assertEqual(self.cli("run")[0], 75)
        self.assertEqual(self.cli("run")[0], 75)
        rc, out, err = self.cli("run", "--budget-s", "1")
        self.assertEqual(rc, 75, err)
        res = json.loads(out)
        self.assertEqual((res["group"], len(res["done"]), len(res["pending"])), ("verdict", 22, 5))      # 21 control units were done before
        rc, out, err = self.cli("run")
        self.assertEqual(rc, 0, err)

    def test_the_strict_single_copy_reading_reaches_the_gate_and_the_build(self):
        rows = T.outgroup_rows(*[(f"gene-H{n}", f"h{n}", 1000, 3000, []) for n in (1, 2, 3, 4)])
        rows.append(row("gene-Hsrc", contig="h8", start0=0, end=500))
        rows.append(row("gene-Hsrc", cls="dropped_M1", contig="h9", start0=0, end=500, note="overlaps gene-H1 (h1)"))
        rows[-1]["tag"] = "1"
        T.write_table(self.og["H"]["selflift"], rows)
        self.assertEqual(self.cli("run")[0], 75)
        rc, out, err = self.cli("run", "--strict-single-copy")
        self.assertEqual(rc, 3)
        self.assertIn("Gate 7 FAILED", err)
        rcs, err = self.to_the_end()
        self.assertEqual(rcs[-1], 0, err)
        rc, out, err = self.cli("build", "--strict-single-copy")
        self.assertEqual(rc, 2, err)                                       # a strict build needs a table directory of its own
        rc, out, err = self.cli("build", "--strict-single-copy", "--table-dir", os.path.join(self.e_root, "table_strict"))
        self.assertEqual(rc, 4, err)
        self.assertEqual(json.loads(out)["status"][:8], "INVALID:")
        self.assertIs(json.load(open(os.path.join(self.e_root, "table_strict", "summary.json")))["strict_single_copy"], True)

    def test_the_build_prints_its_result_and_exits_zero_only_for_a_valid_table(self):
        rcs, err = self.to_the_end()
        self.assertEqual(rcs[-1], 0, err)
        rc, out, err = self.cli("build")
        self.assertEqual(rc, 0, err)
        printed = json.loads(out)
        self.assertEqual((printed["status"], printed["expansions"], printed["gate7"]), ("VALID", 1, True))

    def test_a_build_before_the_runs_are_finished_is_a_stop_with_exit_code_3(self):
        rc, out, err = self.cli("build")
        self.assertEqual(rc, 3)
        self.assertTrue(err.startswith("STOP:"))
        self.assertEqual(out, "")

    def test_a_build_with_another_hub_reading_than_the_runs_is_a_stop_with_exit_code_3(self):
        rcs, err = self.to_the_end()
        rc, out, err = self.cli("build", "--hub-mode", "records")
        self.assertEqual(rc, 3)
        self.assertIn("plan digest", err)

    def test_the_gffs_command_needs_no_amendment_and_writes_under_gff(self):
        with mock.patch.object(LC, "check_committed_code", side_effect=AssertionError("not asked")):
            rc, out, err = self.cli("gffs", "--gorilla-gff", self.gff_path)
        self.assertEqual(rc, 0, err)
        self.assertIn("run GFFs for", out)
        self.assertTrue(os.path.exists(os.path.join(self.e_root, "gff", "pos1", "in.gff")))
        self.assertTrue(os.path.exists(os.path.join(self.e_root, "gff", "r01", "in.gff")))
        self.assertTrue(os.path.exists(os.path.join(self.e_root, "gff", "ctl2", "in.gff")))

    def test_an_amendment_commit_is_seven_to_forty_lowercase_hex_digits_and_nothing_else(self):
        def not_reached(*a, **k):
            raise AssertionError("the plan must not be read for a refused command")
        with mock.patch.object(LC, "plan_from_real_inputs", not_reached), contextlib.redirect_stderr(io.StringIO()):
            for bad in ("", "abc12", "abc123", "ABCDEF1", "abcdef1x", "abcdef1 ", "g" * 10, "a" * 41):
                argv = ["run", "--out", self.e_root] + (["--amendment-commit", bad] if bad else [])
                self.assertEqual(LC.main(argv), 2, bad)


class CliWiringTests(CliBase):
    """the committed-code and identity checks of the run command, with the git check replaced by a mock"""

    def run_main(self, *extra, check=None, identity=None, tools=None):
        argv = ["run", "--out", self.e_root, "--lift-root", self.root, "--n-planted", "4", "--amendment-commit", "abcdef1", "--gorilla-gff", self.gff_path, "--gorilla-fasta", self.fasta,
                "--liftoff-bin", self.stub, "--minimap2", "/m/mm2", *extra]
        out, err = io.StringIO(), io.StringIO()
        calls = []
        self.needles = None

        def fake_check(commit, *a, **k):
            calls.append(commit)
            self.needles = k.get("needles")
            if check:
                raise LC.StopRun(check)
            return "HEAD"
        tools = {"liftoff": "1" * 64, "minimap2": "2" * 64, "gnu_time": "3" * 64} if tools is None else tools
        with mock.patch.object(LC, "check_committed_code", fake_check), mock.patch.object(LC, "code_identity", lambda *a, **k: identity or {"code.py": "1" * 64}), \
                mock.patch.object(LC, "resolve_commit", lambda arg, repo=None: arg), \
                mock.patch.object(LC, "tool_identity", lambda *a, **k: tools), contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = LC.main(argv)
        return rc, out.getvalue(), err.getvalue(), calls

    def test_the_code_is_checked_and_its_identity_frozen_before_anything_is_planned_or_run(self):
        rc, out, err, calls = self.run_main()
        self.assertEqual((rc, calls), (75, ["abcdef1"]))
        ident = json.loads(rd(os.path.join(self.e_root, "run_identity.json")))
        self.assertEqual((ident["amendment_commit"], ident["code_sha256"]), ("abcdef1", {"code.py": "1" * 64}))
        self.assertEqual(sorted(ident), ["amendment_commit", "code_sha256", "command", "inputs", "tools"])           # the identity also pins the tools, the command and the inputs (D8)
        self.assertTrue(os.path.exists(os.path.join(self.e_root, "plan.digest")))

    def test_the_amendment_must_list_the_code_files_and_the_tools(self):
        rc, out, err, calls = self.run_main()
        self.assertEqual(rc, 75, err)
        ident = json.loads(rd(os.path.join(self.e_root, "run_identity.json")))
        base = {"code.py": "1" * 64, "tool liftoff": "1" * 64, "tool minimap2": "2" * 64, "tool gnu_time": "3" * 64}
        # fix round 2 (G6): on the first call the Amendment must also list the plan digest and the sha256 of every pinned input
        self.assertEqual({k: v for k, v in self.needles.items() if k in base}, base)
        self.assertEqual(set(self.needles) - set(base), {"plan digest"} | {f"input {n}" for n, p in ident["inputs"].items() if p["sha256"]})
        rc, out, err, calls = self.run_main()                                                  # a later call: the identity froze the plan and the inputs
        self.assertEqual(self.needles, base)

    def test_a_tool_that_does_not_exist_stops_the_call_before_the_code_check_and_before_anything_is_recorded(self):
        rc, out, err, calls = self.run_main(tools={"liftoff": "1" * 64, "minimap2": "absent", "gnu_time": "3" * 64})
        self.assertEqual((rc, calls), (3, []))
        self.assertIn("minimap2", err)
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "run_identity.json")))
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "plan.digest")))

    def test_a_later_call_with_other_code_is_a_stop_and_runs_nothing(self):
        self.assertEqual(self.run_main()[0], 75)
        before = sorted(os.listdir(os.path.join(self.e_root, "lift", "human", "runs")))
        rc, out, err, calls = self.run_main(identity={"code.py": "2" * 64})
        self.assertEqual(rc, 3)
        self.assertIn("run_identity.json", err)
        self.assertEqual(sorted(os.listdir(os.path.join(self.e_root, "lift", "human", "runs"))), before)

    def test_uncommitted_code_stops_the_call_before_the_identity_or_the_plan_is_recorded(self):
        rc, out, err, calls = self.run_main(check="x.py differs from HEAD")
        self.assertEqual(rc, 3)
        self.assertIn("differs from HEAD", err)
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "run_identity.json")))
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "plan.digest")))
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "lift")))

    def test_the_no_code_check_flag_is_off_unless_given(self):
        self.assertEqual(self.run_main("--no-code-check")[3], [])
        os.remove(os.path.join(self.e_root, "run_identity.json")) if os.path.exists(os.path.join(self.e_root, "run_identity.json")) else None
        self.assertEqual(self.run_main()[3], ["abcdef1"])


class CliDefaultsTests(T.WorldCase):
    """what a command without options does: the configuration of the real run (every file access is replaced by a recorder)"""

    def setUp(self):
        super().setUp()
        self.seen = {}
        seen = self.seen

        def fake_plan(hub_mode, lift_root, spec=None, n_neg=100):
            seen["plan"] = (hub_mode, lift_root, n_neg)
            return self.w

        def fake_run(units, plan_by_label, **kw):
            seen["run"] = dict(kw, units=units, plan_by_label=plan_by_label)
            return dict(done=[], pending=[], split=None, deadline=None, group="pos", exit_code=0)
        rec_ = lambda name, ret=None: (lambda *a, **k: (seen.__setitem__(name, (a, k)), ret)[1])  # noqa: E731

        class FakeLock:
            """the driver lock without a file: main() would create E_ROOT/run.lock"""
            fd = 7

            def __init__(self, e_root):
                seen["lock"] = e_root

            def __enter__(self):
                return self

            def __exit__(self, *exc):
                return False
        pathlib = __import__("pathlib")
        patches = [mock.patch.object(LC, "plan_from_real_inputs", fake_plan), mock.patch.object(LC, "run_budgeted", fake_run),
                   # the steps the fix round added in front of the first launch, replaced by recorders so that main() touches no file (tools, inputs, lock, links):
                   mock.patch.object(LC, "DriverLock", FakeLock), mock.patch.object(LC, "guard_no_stray", rec_("stray")), mock.patch.object(LC, "check_liftoff_binary", rec_("binary")),
                   mock.patch.object(LC, "check_reference", rec_("reference")), mock.patch.object(LC, "check_outgroup_files", rec_("outgroup_files")),
                   mock.patch.object(LC, "tool_identity", lambda *a, **k: {"liftoff": "1" * 64, "minimap2": "2" * 64, "gnu_time": "3" * 64}),
                   mock.patch.object(LC, "pin_inputs", lambda paths, deep=True: {}), mock.patch.object(LC, "load_identity", lambda e: None),
                   mock.patch.object(LC, "verify_identity", rec_("verify_identity")),
                   mock.patch.object(LC, "make_lift_dir", lambda e: (seen.__setitem__("lift_dir", e), pathlib.Path(e) / "lift")[1]),
                   mock.patch.object(LC, "append_run_log", rec_("run_log")),
                   mock.patch.object(LC, "ensure_run_gffs", rec_("gffs", True)), mock.patch.object(LC, "freeze_plan", rec_("freeze_plan")),
                   mock.patch.object(LC, "freeze_identity", rec_("freeze_identity")), mock.patch.object(LC, "check_committed_code", rec_("check", "HEAD")),
                   mock.patch.object(LC, "resolve_commit", lambda arg, repo=None: arg),
                   mock.patch.object(LC, "code_identity", lambda *a, **k: {"f": "1" * 64}), mock.patch.object(LC, "write_run_gffs", rec_("write_run_gffs", {"r01": 1})),
                   mock.patch.object(LC, "make_gate", lambda w, e, og, fn, **kw: ("gate", e, og, fn)), mock.patch.object(LC.L, "mm2_wrapper", lambda root: seen.__setitem__("mm2_root", root) or "/sentinel/mm2.sh"),
                   mock.patch.object(LC, "build_all", lambda *a, **k: (seen.__setitem__("build", (a, k)), self.build_result)[1]),
                   mock.patch.object(LC, "describe_plans", lambda w, budget: (seen.__setitem__("describe_budget", budget), {"ok": True})[1])]
        for pt in patches:
            pt.start()
            self.addCleanup(pt.stop)
        self.build_result = dict(status="VALID", expansions=3, gate7=dict(valid=True), paths={}, valid=True, dropped_classes=0, records_in_dropped_classes=0)

    def main(self, *argv):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = LC.main(list(argv))
        return rc, out.getvalue(), err.getvalue()

    def test_run_defaults(self):
        rc, out, err = self.main("run", "--amendment-commit", "abcdef1")
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.seen["plan"], ("nodes", X.LIFT_ROOT, 100))
        kw = self.seen["run"]
        self.assertEqual((kw["e_root"], kw["reference"], kw["liftoff_bin"], kw["minimap2"], kw["budget_s"]),
                         (LC.E_ROOT, LC.GORILLA_FASTA, LC.LIFTOFF_BIN, "/sentinel/mm2.sh", 530.0))
        self.assertIs(kw["outgroups"], LC.OUTGROUPS)
        self.assertEqual(kw["group_of"], LC.group_map(self.w["plan"], self.w["pos"], self.w["ctl"]))
        self.assertEqual(kw["gate"][:3], ("gate", LC.E_ROOT, LC.OUTGROUPS))
        self.assertEqual(self.seen["mm2_root"], __import__("pathlib").Path(LC.E_ROOT) / "lift")
        self.assertEqual(self.seen["gffs"][0][0], LC.GORILLA_GFF)
        self.assertEqual(self.seen["gffs"][0][2], LC.E_ROOT + "/gff")
        self.assertEqual((kw["t0"], kw["pass_fds"], sorted(kw["wall_pred_s"])), (LC.CALL_T0, (7,), ["C", "H", "O"]))        # the call starts with the process; the lock is inherited; stop rule 6 is on
        self.assertTrue(callable(kw["on_first_launch"]))
        self.assertEqual(self.seen["lock"], LC.E_ROOT)
        self.assertEqual(self.seen["check"][0], ("abcdef1",))
        self.assertEqual(self.seen["freeze_identity"][0], (LC.E_ROOT, {"f": "1" * 64}, "abcdef1"))
        self.assertEqual(self.seen["freeze_plan"][0][0], LC.E_ROOT)
        self.assertEqual(self.seen["freeze_plan"][0][1], LC.plan_digest(self.w["plan"] + self.w["pos"] + self.w["ctl"], "nodes"))
        assigned = self.seen["gffs"][0][1]
        self.assertEqual(sorted(set(assigned.values())), sorted(r["label"] for r in self.w["plan"] + self.w["pos"] + self.w["ctl"]))
        self.assertEqual(sorted(kw["plan_by_label"]), sorted(r["label"] for r in self.w["plan"] + self.w["pos"] + self.w["ctl"]))

    def test_run_options_are_passed_on(self):
        rc, out, err = self.main("run", "--amendment-commit", "abcdef1", "--out", "/o", "--budget-s", "123", "--hub-mode", "records", "--liftoff-bin", "/b/l", "--gorilla-fasta", "/f.fa",
                                 "--gorilla-gff", "/a.gff", "--minimap2", "/m/mm2", "--lift-root", "/lr", "--n-planted", "6", "--strict-single-copy")
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.seen["plan"], ("records", "/lr", 6))
        kw = self.seen["run"]
        self.assertEqual((kw["e_root"], kw["reference"], kw["liftoff_bin"], kw["minimap2"], kw["budget_s"]), ("/o", "/f.fa", "/b/l", "/m/mm2", 123.0))
        self.assertEqual(self.seen["gffs"][0][0], "/a.gff")
        self.assertEqual(self.seen["gffs"][0][2], "/o/gff")
        self.assertEqual(self.seen["freeze_plan"][0][1], LC.plan_digest(self.w["plan"] + self.w["pos"] + self.w["ctl"], "records"))
        with mock.patch.object(LC, "load_indexes", lambda og, strict=False: ("indexes", strict)):
            self.assertEqual(kw["gate"][3](), ("indexes", True))
        self.assertNotIn("mm2_root", self.seen)

    def test_the_gate_reads_the_loose_reading_by_default(self):
        self.main("run", "--amendment-commit", "abcdef1")
        with mock.patch.object(LC, "load_indexes", lambda og, strict=False: ("indexes", strict)):
            self.assertEqual(self.seen["run"]["gate"][3](), ("indexes", False))

    def test_plan_and_gffs_defaults(self):
        rc, out, err = self.main("plan")
        self.assertEqual((rc, self.seen["plan"], self.seen["describe_budget"], json.loads(out)), (0, ("nodes", X.LIFT_ROOT, 100), 530.0, {"ok": True}))
        rc, out, err = self.main("gffs")
        self.assertEqual(rc, 0)
        args = self.seen["gffs"][0]                                        # gffs goes through ensure_run_gffs, which also writes the assignment marker
        self.assertEqual((args[0], args[2]), (LC.GORILLA_GFF, LC.E_ROOT + "/gff"))

    def test_build_defaults_and_exit_codes(self):
        rc, out, err = self.main("build")
        self.assertEqual(rc, 0)
        a, k = self.seen["build"]
        self.assertEqual((a[0], a[1]), (LC.E_ROOT, LC.E_ROOT + "/table"))
        self.assertEqual((k["lift_root"], k["hub_mode"], k["n_neg"], k["strict"]), (X.LIFT_ROOT, "nodes", 100, False))
        self.assertIs(k["outgroups"], LC.OUTGROUPS)
        self.assertEqual(json.loads(out), dict(status="VALID", expansions=3, dropped_classes=0, records_in_dropped_classes=0, gate7=True))
        self.build_result = dict(status="INVALID: x", expansions=0, gate7=dict(valid=False), paths={}, valid=False, dropped_classes=0, records_in_dropped_classes=0)
        self.assertEqual(self.main("build", "--table-dir", "/t", "--strict-single-copy")[0], 4)
        self.assertEqual(self.seen["build"][0][1], "/t")
        self.assertIs(self.seen["build"][1]["strict"], True)
        with mock.patch.object(LC, "build_all", side_effect=LC.StopRun("no")):
            rc, out, err = self.main("build")
        self.assertEqual((rc, err), (3, "STOP: no\n"))



# ===================================================================================================== O. second round: the finer rules
class RecordingPopen(subprocess.Popen):
    """records the timeouts of wait() on the supervised process (the one run under /usr/bin/time), not on helpers such as ps"""
    waits = []

    def __init__(self, args, *a, **k):
        self._main = bool(args) and args[0] in ("/usr/bin/timeout", "/usr/bin/time")
        super().__init__(args, *a, **k)

    def wait(self, timeout=None):
        if self._main:
            RecordingPopen.waits.append((timeout, time.monotonic()))
        return super().wait(timeout)


class WatchTimingTests(WatchBase):
    def test_no_wait_is_longer_than_the_time_left_to_the_hard_limit(self):
        stub = self.stub("long.sh", 'exec sleep 30\n')
        RecordingPopen.waits = []
        t_start = time.monotonic()
        with mock.patch.object(LC.subprocess, "Popen", RecordingPopen):
            rc, peak, elapsed = self.watch(stub, hard_s=1.0, poll_s=0.4)
        self.assertEqual(rc, LC.RC_TOO_LONG)
        polls = [(t, at) for t, at in RecordingPopen.waits[:-1] if t is not None]          # the last wait reaps the killed process
        self.assertGreaterEqual(len(polls), 3)
        for timeout, at in polls:
            self.assertLessEqual(timeout, 1.0 - (at - t_start) + 0.05, (timeout, at - t_start))


class PackingBoundaryTests(T.HermeticCase):
    def test_a_unit_a_hair_over_the_budget_starts_a_new_call_in_the_first_call_and_in_every_later_call(self):
        a, b, x = ("H", "a"), ("H", "b"), ("H", "x")
        self.assertEqual(LC.plan_calls([a, b], {a: 300, b: 230.05}), [[a], [b]])
        self.assertEqual(LC.plan_calls([x, a, b], {x: 600, a: 300, b: 230.05}), [[x], [a], [b]])
        self.assertEqual(LC.plan_calls([a, b], {a: 300, b: 230.0}), [[a, b]])
        self.assertEqual(LC.plan_calls([x, a, b], {x: 600, a: 300, b: 230.0}), [[x], [a, b]])


class RunGffDetailTests(T.HermeticCase):
    def test_the_run_gff_mark_ignores_the_insertion_order_of_the_assignment_and_the_return_is_a_bool(self):
        tmp = tempfile.mkdtemp()
        gff = os.path.join(tmp, "ref.gff")
        with open(gff, "w") as fh:
            fh.write(T.REF_GFF)
        out = os.path.join(tmp, "gff")
        self.assertIs(LC.ensure_run_gffs(gff, {"gene-A": "r01", "gene-B": "r02"}, out), True)
        self.assertIs(LC.ensure_run_gffs(gff, {"gene-B": "r02", "gene-A": "r01"}, out), False)

    def test_an_unknown_amendment_commit_is_named_as_such(self):
        repo = tempfile.mkdtemp()
        T.git(repo, "init", "-q")
        with open(os.path.join(repo, "a.py"), "w") as fh:
            fh.write("x\n")
        T.git(repo, "add", "-A")
        T.git(repo, "commit", "-q", "-m", "x")
        with self.assertRaises(LC.StopRun) as cm:
            LC.check_committed_code("0" * 40, repo, ("a.py",))
        self.assertIn("is not a commit", str(cm.exception))

    def test_the_ids_of_a_shard_keep_their_genome_order(self):
        recs = [dict(id="gene-z", contig="c1", start0=0, end=10), dict(id="gene-a", contig="c1", start0=100, end=110)]
        by_id = {r["id"]: r for r in recs}
        self.assertEqual(LC.effective_shards([("A", ["gene-z", "gene-a"])], by_id, []), [("A", ["gene-z", "gene-a"])])


class PoolDetailTests(T.HermeticCase):
    def test_the_sort_key_of_records_is_contig_start_end_id(self):
        a, b = rec("a", "c1", 0, 5000), rec("b", "c1", 100, 400)
        self.assertEqual([r["id"] for r in sorted([b, a], key=X.record_key)], ["a", "b"])
        self.assertEqual(X.record_key(a), ("c1", 0, 5000, "a"))
        comp = X.Component("PC_x")
        for r in (b, a):
            r["cls"] = "in_place"
            comp.recs[r["id"]] = r
        self.assertEqual([r["id"] for r in comp.such()], ["a", "b"])
        plan = X.rank_plan([comp])
        self.assertEqual([r["records"][0]["id"] for r in plan], ["a", "b"])

    def test_the_members_of_a_component_are_sorted_lists(self):
        rows, recs = [row("g-0", contig="c0")], {"g-0": rec("g-0", "c0", 0, 1000)}
        for i in range(1, 9):
            rid = f"g-{i}"
            recs[rid] = rec(rid, f"c{i}", 0, 1000)
            rows.append(row(rid, contig=f"c{i}"))
            rows.append(row(rid, cls="dropped_M1", contig="cx", start0=i, end=i + 100, note="overlaps g-0 (c0)"))
        for i in range(8):
            rows.append(row("g-0", cls="extra_copy", contig=f"s{(i * 5) % 8}", start0=0, end=500))
        for r in rows:
            if r["cls"] in X.COPY_CLASSES:
                r["tag"] = "1"
        comp = next(iter(X.build_pool(rows, recs).components.values()))
        self.assertEqual(comp.record_ids, sorted(comp.record_ids))
        self.assertEqual(len(comp.record_ids), 9)
        self.assertEqual(comp.sites, sorted(comp.sites))
        self.assertEqual(len(comp.sites), 8)

    def test_the_number_of_copy_rows_is_counted_and_reported(self):
        rows, recs = T.mini_selflift()
        self.assertEqual(X.build_pool(rows, recs).n_copy_rows, 4)
        root = tempfile.mkdtemp()
        T.write_selflift(root, rows, recs)
        self.assertEqual(X.pool_report(root, support_libs=())["copy_rows"], 4)
        self.assertEqual(X.Pool().n_copy_rows, 0)

    def test_an_annotated_record_may_start_with_the_letters_copy(self):
        rows = [row("g-a", contig="c1"), row("copy-gene", contig="c2"), row("copy-gene", cls="dropped_M1", contig="c1", start0=0, end=100, note="overlaps g-a (c1)")]
        rows[-1]["tag"] = "1"
        comp = next(iter(X.build_pool(rows, {"g-a": rec("g-a", "c1", 0, 1000), "copy-gene": rec("copy-gene", "c2", 0, 1000)}).components.values()))
        self.assertEqual((comp.record_ids, comp.sites), (["copy-gene", "g-a"], []))

    def test_a_record_without_a_primary_row_is_a_member_but_never_a_such_record(self):
        rows = [row("g-a", contig="c1"), row("g-a", cls="dropped_M1", contig="c2", start0=0, end=500, note="overlaps g-b (c2)")]
        rows[-1]["tag"] = "1"
        comp = next(iter(X.build_pool(rows, {"g-a": rec("g-a", "c1", 0, 1000), "g-b": rec("g-b", "c2", 0, 1000)}).components.values()))
        self.assertEqual(comp.recs["g-b"]["cls"], "unmapped")
        self.assertEqual([r["id"] for r in comp.such()], ["g-a"])

    def test_a_record_missing_from_the_record_table_is_rebuilt_from_its_first_row(self):
        rows = [row("g-a", contig="c1", start0=10, end=1010, exons="10-300,500-1010", strand="-", rtype="pseudogene", biotype="transcribed_pseudogene"), row("g-b", contig="c2")]
        rows.append(row("g-b", cls="dropped_M1", contig="c1", start0=0, end=500, note="overlaps g-a (c1)"))
        rows[-1]["tag"] = "1"
        comp = next(iter(X.build_pool(rows, {}).components.values()))
        a = comp.recs["g-a"]
        self.assertEqual((a["id"], a["name"], a["rtype"], a["biotype"], a["contig"], a["start0"], a["end"], a["strand"], a["exons"], a["exonic_len"], a["idx"], a["cls"]),
                         ("g-a", "g-a", "pseudogene", "transcribed_pseudogene", "c1", 10, 1010, "-", [(10, 300), (500, 1010)], 800, 0, "in_place"))
        row_without = dict(rows[0])
        del row_without["rtype"]
        del row_without["name"]
        got = X._record_from_row(row_without)
        self.assertEqual((got["rtype"], got["name"]), ("gene", "g-a"))

    def test_rank_runs_list_the_components_by_id_whatever_the_order_of_the_input(self):
        def comp(cid):
            c = X.Component(cid)
            for k in range(2):
                r = rec(f"{cid}-{k}", f"ctg-{cid}", k * 10000, k * 10000 + 2000)
                r["cls"] = "in_place"
                c.recs[r["id"]] = r
            return c
        plan = X.rank_plan([comp("PC_b"), comp("PC_c"), comp("PC_a")])
        self.assertEqual([[r["component"] for r in run["records"]] for run in plan], [["PC_a", "PC_b", "PC_c"]] * 2)

    def test_a_plan_of_exactly_cap_records_has_the_label_of_its_run(self):
        self.assertEqual(X.plan_run_shards("r01", shard_recs(4), cap=4), [("r01", [f"g{i:03d}" for i in range(4)])])

    def test_the_status_of_a_component_that_is_neither_verdict_nor_control_nor_hub(self):
        comps = []
        for cid, n in (("PC_two", 2), ("PC_one", 1), ("PC_three", 3)):
            c = X.Component(cid)
            for k in range(n):
                r = rec(f"{cid}-{k}", f"ctg-{cid}", k * 10000, k * 10000 + 2000)
                r["cls"] = "in_place"
                c.recs[r["id"]] = r
                c.record_ids.append(r["id"])
            comps.append(c)
        pool = X.Pool()
        for c in comps:
            pool.components[c.id] = c
        tmp = tempfile.mkdtemp()
        paths = X.write_outputs(X.call_expansions([], {}), pool, [], {}, tmp)
        with open(os.path.join(tmp, "components.tsv")) as fh:
            by = {r[0]: r[5] for r in list(csv.reader(fh, delimiter="\t"))[1:]}
        self.assertEqual(by, {"PC_two": "not_called", "PC_one": "too_few_such_records", "PC_three": "not_called"})
        self.assertEqual(sorted(paths), sorted(["etable.tsv", "member_exons.tsv", "control_pool.tsv", "extension.tsv", "siblings.tsv", "record_calls.tsv", "dropped_classes.tsv",
                                                "components.tsv", "etable_human_only.tsv", "summary.json"]))
        for n, path in paths.items():
            self.assertEqual(path, os.path.join(tmp, n))
            self.assertTrue(os.path.exists(path))

    def test_joining_copy_rows_keep_duplicates_and_are_sorted_not_in_file_order(self):
        rows = [row("g-a", contig="c1"), row("g-b", contig="c2"), row("g-c", contig="c3")]
        extra = [row("g-c", cls="dropped_M1", contig="c1", start0=0, end=100, note="overlaps g-a (c1)", seqid="0.99"),
                 row("g-a", cls="dropped_M1", contig="c2", start0=0, end=100, note="overlaps g-b (c2)", seqid="0.97"),
                 row("g-a", cls="dropped_M1", contig="c2", start0=500, end=600, note="overlaps g-b (c2)", seqid="0.97")]
        for r in extra:
            r["tag"] = "1"
        recs = {i: rec(i, c, 0, 1000) for i, c in (("g-a", "c1"), ("g-b", "c2"), ("g-c", "c3"))}
        j = X.joining_copy_rows(X.build_pool(rows + extra, recs))
        self.assertEqual(j["g-a"], [("g-b", "dropped_M1", "0.97"), ("g-b", "dropped_M1", "0.97"), ("g-c", "dropped_M1", "0.99")])
        self.assertEqual(j["g-b"], [("g-a", "dropped_M1", "0.97"), ("g-a", "dropped_M1", "0.97")])

    def test_the_index_names_the_first_copy_row_candidate_and_counts_copy_rows_of_a_record(self):
        rows = []
        for rid, s0 in (("m2", 1200), ("m1", 1000)):
            rows.append(row(rid, contig="h1", start0=s0, end=s0 + 2000))
            rows.append(row(rid, cls="extra_copy", contig="x9", start0=0, end=500))
        idx = X.OutgroupIndex(rows, "H")
        m = idx.match("h1", [(1200, 3000)])
        self.assertEqual((m["status"], m["locus"], m["candidates"]), ("outgroup_copy_row", "m1", ["m1", "m2"]))
        lifted = {k: {"g1": v} for k, v in T.lifted_for("g1", ("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000)).items()}
        lifted["C"]["g1"] = [lrow("g1", 0, "c1", 1000, 3000), lrow("g1_1", 1, "c9", 0, 100, cov=0.1), lrow("g1_2", 2, "c8", 0, 100, cov=0.1)]
        c = X.classify_record("g1", lifted, {"H": set(), "C": set(), "O": set()}, T.make_indexes())
        self.assertEqual((c["per"]["C"]["copy_rows"], c["per"]["C"]["n_place"]), (2, 1))


class ReadLiftedDetailTests(T.HermeticCase):
    def test_record_ids_with_underscores_keep_their_copy_rows(self):
        import pathlib
        import _liftoff
        tmp = tempfile.mkdtemp()
        d = os.path.join(tmp, "runs", "r01")
        os.makedirs(d)
        with open(os.path.join(d, "lifted.gff3"), "w") as fh:
            fh.write(T.lifted_gff([("gene-A_B", "h1", 1001, 2000, "+", "0.99", "0.98", 0, "False", "False", [(1001, 2000)], "gene"),
                                   ("gene-A_B_1", "h2", 501, 1500, "+", "0.9", "0.96", 1, "False", "False", [(501, 1500)], "gene"),
                                   ("gene-C", "h3", 1001, 2000, "+", "0.99", "0.98", 0, "False", "False", [(1001, 2000)], "gene"),
                                   ("gene-C_1", "h4", 1001, 2000, "+", "0.99", "0.98", 1, "False", "False", [(1001, 2000)], "gene")]))
        _liftoff.parse_lifted(pathlib.Path(d) / "lifted.gff3", pathlib.Path(d) / "loci.tsv")
        with open(os.path.join(d, "unmapped.txt"), "w") as fh:
            fh.write("gene-U_1\n\n")
        with open(os.path.join(d, "done.key"), "w") as fh:
            fh.write("{}")
        lifted, unmapped = X.read_lifted_runs(tmp, known_ids={"gene-A_B", "gene-C", "gene-U_1"})
        self.assertEqual(sorted(lifted), ["gene-A_B", "gene-C"])
        self.assertEqual(sorted(r["feature_id"] for r in lifted["gene-A_B"]), ["gene-A_B", "gene-A_B_1"])
        self.assertEqual(sorted(r["feature_id"] for r in lifted["gene-C"]), ["gene-C", "gene-C_1"])
        self.assertEqual(unmapped, {"gene-U_1"})


class PositiveControlDetailTests(T.HermeticCase):
    def setUp(self):
        self.pool, self.rows, self.spec = T.lrpap_fixture()
        self.ctrl = X.lrpap_control(self.pool, self.rows, self.spec)
        self.called = X.called_components(self.pool, "nodes")

    def calls(self, **triples):
        calls = {r["id"]: dict(triple=None, reason="no_placement", per={k: dict(n_place=1, detail="") for k in "HCO"}) for c in self.called for r in c.such()}
        for rid, t in triples.items():
            calls[rid.replace("_", "-")]["triple"], calls[rid.replace("_", "-")]["reason"] = t, None
        return calls

    def check(self, **triples):
        calls = self.calls(**triples)
        res = X.call_expansions(self.called, calls)
        return X.check_positive_control(self.ctrl, res, calls), res

    def test_the_target_class_is_the_one_with_the_most_qualifying_records_not_the_biggest(self):
        out, res = self.check(**{"gene_REF": T.T1, "gene_X1": T.T1, "gene_X3": T.T2, "gene_X2": T.T2, "gene_N": T.T2})
        self.assertFalse(out["pass"])
        self.assertEqual(sorted(out["class_members"]), ["gene-REF", "gene-X1"])
        self.assertEqual(out["missing"], [dict(record="gene-X3", reason="in another class")])

    def test_the_target_class_is_the_best_one_even_when_it_comes_later(self):
        out, res = self.check(**{"gene_REF": T.T1, "gene_N": T.T1, "gene_X1": T.T2, "gene_X3": T.T2, "gene_X2": T.T2})
        self.assertEqual([e["triple"] for e in res["expansions"]], [T.T1, T.T2])
        self.assertEqual(sorted(out["class_members"]), ["gene-X1", "gene-X2", "gene-X3"])
        self.assertEqual(out["expansion"], "EX0002")
        self.assertEqual(out["missing"], [dict(record="gene-REF", reason="in another class")])

    def test_equal_overlaps_go_to_the_earlier_class(self):
        out, res = self.check(**{"gene_REF": T.T1, "gene_N": T.T1, "gene_X1": T.T2, "gene_X2": T.T2})
        self.assertEqual(out["expansion"], "EX0001")
        self.assertEqual(sorted(out["class_members"]), ["gene-N", "gene-REF"])

    def test_a_class_with_no_qualifying_record_is_no_target(self):
        out, res = self.check(**{"gene_X2": T.T1, "gene_N": T.T1})
        self.assertEqual(len(res["expansions"]), 1)
        self.assertFalse(out["pass"])
        self.assertIsNone(out["expansion"])
        self.assertEqual(out["class_members"], [])
        self.assertEqual({m["record"] for m in out["missing"]}, {"gene-REF", "gene-X1", "gene-X3"})

    def test_a_class_with_a_single_qualifying_record_is_the_target(self):
        out, res = self.check(**{"gene_REF": T.T1, "gene_N": T.T1})
        self.assertEqual(out["expansion"], "EX0001")
        self.assertEqual({m["record"] for m in out["missing"]}, {"gene-X1", "gene-X3"})

    def test_class_members_are_listed_in_sorted_order(self):
        ids = [f"gene-m{k}" for k in (7, 2, 9, 4, 0, 8, 5, 1, 6, 3)]
        ctrl = dict(qualifying=list(ids), not_required=[], fragments=dict(component="PC_f", records=[]), expected_triple=T.T1)      # addendum 2: a control has its expected triple
        res = dict(expansions=[dict(id="EX0001", triple=T.T1, component="PC_a", members=[dict(id=i) for i in ids])])
        out = X.check_positive_control(ctrl, res, {})
        self.assertTrue(out["pass"])
        self.assertEqual(out["class_members"], sorted(ids))

    def test_the_share_of_a_record_counts_its_exon_union_once(self):
        recs = T.lrpap_fixture_records()
        recs["gene-DUP"] = rec("gene-DUP", "c2", 5000, 8000, exons=[(5000, 5100), (5050, 5150)])        # union 150 bp, 100 of them in the interval 5000-5100... see below
        rows = list(self.rows) + [row("gene-DUP", contig="c2", start0=5000, end=8000, exons="5000-5100,5050-5150"),
                                  row("gene-DUP", cls="dropped_M1", contig="c1", start0=0, end=100, note="overlaps gene-REF (c1)")]
        rows[-1]["tag"] = "1"
        rows = [r for r in rows if not (r["source_id"] == "gene-REF" and r["contig"] == "c2")]
        rows.append(row("gene-REF", cls="dropped_M1", contig="c2", start0=5000, end=5100, note="overlaps gene-X1 (c2)"))
        rows[-1]["tag"] = "1"
        ctrl = X.lrpap_control(X.build_pool(rows, recs), rows, dict(self.spec, min_share=0.6))
        self.assertAlmostEqual(ctrl["shares"]["gene-DUP"], 100 / 150)
        self.assertIn("gene-DUP", ctrl["qualifying"])


class PlantedDetailTests(T.HermeticCase):
    def test_planting_works_from_the_rows_alone_when_the_record_table_lacks_the_records(self):
        rows, recs = T.NegativeControlTests().world()
        pool = X.build_pool(rows, recs)
        full = X.plant_negative(rows, recs, pool)
        bare = X.plant_negative(rows, {}, pool)
        self.assertEqual([[r["id"] for r in c.such()] for c in bare], [[r["id"] for r in c.such()] for c in full])
        for c in bare:
            a, b = c.such()
            self.assertNotEqual(a["contig"], b["contig"])
            self.assertTrue(all(r["exonic_len"] >= 200 and r["cls"] == "in_place" for r in (a, b)))



# ===================================================================================================== P. third round: schedule numbers, report numbers, command constants
class PlanDescriptionTests(T.WorldCase):
    def test_the_schedule_numbers_of_the_plan_command(self):
        d = LC.describe_plans(self.w, 530)
        self.assertEqual((d["units"], d["calls"], d["predicted_work_s"], d["predicted_call_s"]), (27, 3, 680, [378, 151, 151]))
        self.assertEqual(d["groups"], dict(pos=dict(units=15, predicted_work_s=378, calls=1), neg=dict(units=6, predicted_work_s=151, calls=1),
                                           verdict=dict(units=6, predicted_work_s=151, calls=1)))
        self.assertEqual(d["plan_digest"], LC.plan_digest(self.w["plan"] + self.w["pos"] + self.w["ctl"], "nodes"))
        self.assertEqual((d["planted_components"], d["verdict"]["components"], d["verdict"]["records"]), (2, 2, 4))
        self.assertEqual(d["control_components"], {"PC_gene-N": 5, "PC_gene-F1": 3})
        self.assertEqual((d["positive_plan"]["J"], d["positive_plan"]["records"], d["planted_plan"]["records_per_run"]), (5, 8, [2, 2]))
        self.assertEqual(d["called_set"]["components"], 4)
        self.assertEqual(d["called_set"]["records"], 12)
        self.assertEqual(d["lrpap1"]["qualifying"], ["gene-REF", "gene-X1", "gene-X3"])

    def test_a_smaller_budget_makes_more_calls_and_each_group_keeps_its_own(self):
        d = LC.describe_plans(self.w, 100)
        self.assertGreater(d["calls"], 3)
        self.assertEqual(sum(g["calls"] for g in d["groups"].values()), d["calls"])
        self.assertEqual(sum(g["units"] for g in d["groups"].values()), d["units"])
        self.assertEqual(len(d["predicted_call_s"]), d["calls"])

    def test_a_budget_that_fits_no_second_unit_gives_calls_of_one_unit_each(self):
        d = LC.describe_plans(self.w, 1)
        self.assertEqual((d["calls"], d["units"]), (27, 27))
        self.assertEqual({g: v["calls"] for g, v in d["groups"].items()}, dict(pos=15, neg=6, verdict=6))
        self.assertEqual(len(d["predicted_call_s"]), 27)

    def test_the_plan_summary_counts_runs_of_at_most_ten_records(self):
        def plan(*sizes):
            return [dict(records=[dict(id=f"x{i}") for i in range(n)], shards=[("r", [])]) for n in sizes]
        s = X.plan_summary(plan(10, 11, 1, 0))
        self.assertEqual((s["J"], s["records"], s["records_per_run"], s["runs_with_at_most_10"], s["shards"]), (4, 22, [10, 11, 1, 0], 3, 4))


class ReportDetailTests(T.HermeticCase):
    def hub_pool(self, prefix, n_short, n_sites):
        """a component of 2 long records and n_short short ones joined to the first long one, plus n_sites unannotated sites"""
        recs, rows = {}, []
        for k in range(2):
            rid = f"gene-{prefix}L{k}"
            recs[rid] = rec(rid, f"{prefix}1", k * 10000, k * 10000 + 2000)
            rows.append(row(rid, contig=f"{prefix}1", start0=k * 10000, end=k * 10000 + 2000))
        rows.append(row(f"gene-{prefix}L1", cls="dropped_M1", contig=f"{prefix}1", start0=5, end=1900, note=f"overlaps gene-{prefix}L0 ({prefix}1)"))
        for i in range(n_short):
            rid = f"gene-{prefix}S{i}"
            recs[rid] = rec(rid, f"{prefix}2", i * 1000, i * 1000 + 100)
            rows.append(row(rid, contig=f"{prefix}2", start0=i * 1000, end=i * 1000 + 100))
            rows.append(row(rid, cls="dropped_M1", contig=f"{prefix}1", start0=100, end=200, note=f"overlaps gene-{prefix}L0 ({prefix}1)"))
        for i in range(n_sites):
            rows.append(row(f"gene-{prefix}L0", cls="extra_copy", contig=f"{prefix}9", start0=i * 3000, end=i * 3000 + 1000))
        for r in rows:
            if r["cls"] in X.COPY_CLASSES:
                r["tag"] = "1"
        return rows, recs

    def test_hubs_are_listed_by_size_and_supported_hubs_are_counted_under_both_readings(self):
        ra, ca = self.hub_pool("a", 60, 0)             # 62 records: a hub under both readings (the smaller id)
        rb, cb = self.hub_pool("b", 49, 5)             # 51 records + 5 sites = 56 nodes: a hub under 'nodes' only
        rc, cc = self.hub_pool("c", 3, 0)              # 5 records: not a hub
        rows, recs = ra + rb + rc, dict(ca, **cb, **cc)
        root = tempfile.mkdtemp()
        T.write_selflift(root, rows, recs, support={"lib": {"gene-aL0": 3, "gene-aL1": 2, "gene-bL0": 9, "gene-bL1": 9, "gene-cL0": 2, "gene-cL1": 2}})
        rep = X.pool_report(root, support_libs=("lib",))
        self.assertEqual(rep["nodes"]["hub_sizes"], [56, 62])
        self.assertEqual(rep["records"]["hub_sizes"], [62])
        self.assertEqual((rep["nodes"]["hubs"], rep["nodes"]["called"], rep["records"]["hubs"], rep["records"]["called"]), (2, 1, 1, 2))
        self.assertEqual(rep["supported"]["lib"], dict(components=3, hubs_nodes=2, hubs_records=1))
        self.assertEqual(rep["supported"]["either"], dict(components=3, hubs_nodes=2, hubs_records=1))
        self.assertEqual(rep["unannotated_sites"], 5)

    def test_the_200_bp_edge_of_the_report(self):
        rows = [row("g-a", contig="c1", start0=0, end=1000), row("g-b", contig="c2", start0=0, end=200), row("g-c", contig="c3", start0=0, end=199)]
        rows += [row("g-b", cls="dropped_M1", contig="c1", start0=0, end=100, note="overlaps g-a (c1)"), row("g-c", cls="dropped_M1", contig="c1", start0=200, end=300, note="overlaps g-a (c1)")]
        for r in rows[-2:]:
            r["tag"] = "1"
        recs = {"g-a": rec("g-a", "c1", 0, 1000), "g-b": rec("g-b", "c2", 0, 200), "g-c": rec("g-c", "c3", 0, 199)}
        root = tempfile.mkdtemp()
        T.write_selflift(root, rows, recs)
        rep = X.pool_report(root, support_libs=())
        self.assertEqual((rep["annotated_records_in_components"], rep["annotated_ge200bp"], rep["components_ge2_such"]), (3, 2, 1))


class CommandConstantTests(T.HermeticCase):
    def test_every_liftoff_option_reads_its_own_constant(self):
        base = LC.liftoff_command("/t", "/r", "/m", liftoff_bin="/l", threads=4)
        for const, flag, val in (("LIFT_SC", "-sc", 0.71), ("LIFT_A", "-a", 0.62), ("LIFT_S", "-s", 0.53), ("LIFT_OVERLAP", "-overlap", 0.34)):
            with mock.patch.object(LC, const, val):
                cmd = LC.liftoff_command("/t", "/r", "/m", liftoff_bin="/l", threads=4)
            self.assertEqual([c for c in cmd if c != base[cmd.index(flag) + 1] or True].count(flag), 1)
            self.assertEqual(cmd[cmd.index(flag) + 1], str(val), const)
            others = [(f, cmd[cmd.index(f) + 1]) for f in ("-sc", "-a", "-s", "-overlap") if f != flag]
            self.assertEqual(others, [(f, base[base.index(f) + 1]) for f in ("-sc", "-a", "-s", "-overlap") if f != flag])
        self.assertEqual(LC.liftoff_command("/t", "/r", "/m", threads=7)[-4:-2], ["-p", "7"])

    def test_the_default_of_load_indexes_is_the_loose_reading(self):
        self.assertIs(inspect.signature(LC.load_indexes).parameters["strict"].default, False)
        self.assertIs(inspect.signature(LC.plan_from_real_inputs).parameters["hub_mode"].default, "nodes")
        self.assertEqual(inspect.signature(LC.plan_from_real_inputs).parameters["n_neg"].default, 100)


class PositiveControlBigClassTests(T.HermeticCase):
    def test_a_class_of_five_records_that_misses_a_qualifying_one_fails(self):
        ctrl = dict(qualifying=["q1", "q2", "q3"], not_required=[], fragments=dict(component="PC_f", records=[]), expected_triple=T.T1)      # addendum 2
        members = ["q1", "q2", "o1", "o2", "o3"]
        res = dict(expansions=[dict(id="EX0001", triple=T.T1, component="PC_a", members=[dict(id=i) for i in members])])
        calls = {"q3": dict(triple=None, reason="no_placement")}
        out = X.check_positive_control(ctrl, res, calls)
        self.assertFalse(out["pass"])
        self.assertEqual(out["missing"], [dict(record="q3", reason="no_placement")])
        res["expansions"][0]["members"].append(dict(id="q3"))
        self.assertTrue(X.check_positive_control(ctrl, res, calls)["pass"])


if __name__ == "__main__":
    unittest.main()
