#!/usr/bin/env python3
"""Tests of expansion_table.py and liftoff_cross.py (docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md section 7, Gate 7, Amendment 1).

Written before the code (watch-fail-then-pass). Everything here uses made-up tables and stub Liftoff GFF3 text: no Liftoff, no minimap2, no real data.
Run: python3 -B -m unittest bench/hierarchy/test_expansion_table.py
"""
import collections
import contextlib
import copy
import csv
import fcntl
import io
import json
import os
import random
import shutil
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

# The registered paths as the modules define them, taken before any test patches them (the tests that pin the registered values read these)
REAL = types.SimpleNamespace(E_ROOT=LC.E_ROOT, LIFT_ROOT=X.LIFT_ROOT, LIFTOFF_BIN=LC.LIFTOFF_BIN, GORILLA_FASTA=LC.GORILLA_FASTA, GORILLA_GFF=LC.GORILLA_GFF,
                             MINIMAP2_ENV=LC.MINIMAP2_ENV, OUTGROUPS=copy.deepcopy(LC.OUTGROUPS))
VOID = "/nonexistent-hermetic-guard"                   # below a root that cannot be created: a read fails with FileNotFoundError, a write with PermissionError


class HermeticCase(unittest.TestCase):
    """Every test starts with the default paths of the run replaced by paths that cannot exist (E_ROOT, the self-lift root, the Liftoff binary, the gorilla FASTA and
    GFF, the minimap2 of the Liftoff environment, the three outgroups' files): a test that forgets to pass its own path, or a regression that ignores an option, fails
    loudly instead of reading or writing the real run directory, the real tables, the real annotation or the real Liftoff environment (finding 'hermeticity under
    regression' of the mutation review). The patches are applied in run(), so a subclass whose setUp does not call super().setUp() is covered too."""

    def run(self, result=None):
        void_og = {k: dict(t, genome=f"{VOID}/og/{k}/{os.path.basename(t['genome'])}", mmi=f"{VOID}/og/{k}/{os.path.basename(t['mmi'])}",
                           selflift=f"{VOID}/og/{k}/{os.path.basename(t['selflift'])}") for k, t in REAL.OUTGROUPS.items()}
        patches = [mock.patch.object(LC, "E_ROOT", f"{VOID}/E_ROOT"), mock.patch.object(X, "LIFT_ROOT", f"{VOID}/liftoff"),
                   mock.patch.object(LC, "LIFTOFF_BIN", f"{VOID}/liftoff_bin"), mock.patch.object(LC, "GORILLA_FASTA", f"{VOID}/GGO.fasta"),
                   mock.patch.object(LC, "GORILLA_GFF", f"{VOID}/GGO_genomic.gff"), mock.patch.object(LC, "MINIMAP2_ENV", f"{VOID}/minimap2"),
                   mock.patch.dict(LC.OUTGROUPS, void_og, clear=True)]
        for p in patches:
            p.start()
        try:
            return super().run(result)
        finally:
            for p in reversed(patches):
                p.stop()


def rd(path):
    with open(path) as fh:
        return fh.read()


# ----------------------------------------------------------------------------------------------- builders
def rec(rid, contig="c1", start0=0, end=1000, strand="+", exons=None, rtype="gene", biotype="protein_coding", name=None):
    """a record of records.tsv (the annotation) as read_records returns it"""
    exons = exons if exons is not None else [(start0, end)]
    return dict(id=rid, name=name or rid, rtype=rtype, biotype=biotype, contig=contig, start0=start0, end=end, strand=strand,
                exons=list(exons), exonic_len=sum(e - s for s, e in exons), idx=0)


def row(source, cls="in_place", contig=None, start0=0, end=1000, strand="+", note="-", seqid="1.0", cov="1.0", exons=None, short=None, rtype="gene",
        biotype="protein_coding"):
    """a row of liftoff_loci.tsv (the merged self-lift table)"""
    exons = exons if exons is not None else f"{start0}-{end}"
    el = sum(int(b) - int(a) for a, b in (x.split("-") for x in exons.split(",")))
    return dict(species="gorilla", shard="s", source_id=source, name=source, rtype=rtype, biotype=biotype, stratum="protein_coding",
                short="1" if (el < 200 if short is None else short) else "0", src_contig=contig or "c1", src_start0=str(start0), src_end=str(end),
                src_strand=strand, cls=cls, tag="0", contig=contig or "c1", start0=str(start0), end=str(end), strand=strand, coverage=cov,
                sequence_id=seqid, exons=exons, exonic_len=str(el), note=note)


def mini_selflift():
    """records A B C (a component through an M1 note), D E (one through an M2 note), F with an unannotated extra copy site, G alone, H short"""
    recs = {r["id"]: r for r in [
        rec("gene-A", "c1", 1000, 3000), rec("gene-B", "c2", 5000, 7000), rec("gene-C", "c1", 20000, 22000),
        rec("gene-D", "c3", 100, 2100), rec("gene-E", "c3", 9000, 11000), rec("gene-F", "c4", 0, 1500), rec("gene-G", "c4", 5000, 6500),
        rec("gene-H", "c5", 0, 100, rtype="gene", biotype="tRNA")]}
    rows = [row("gene-A", contig="c1", start0=1000, end=3000), row("gene-B", contig="c2", start0=5000, end=7000),
            row("gene-C", contig="c1", start0=20000, end=22000), row("gene-D", contig="c3", start0=100, end=2100),
            row("gene-E", contig="c3", start0=9000, end=11000), row("gene-F", contig="c4", start0=0, end=1500),
            row("gene-G", contig="c4", start0=5000, end=6500), row("gene-H", contig="c5", start0=0, end=100)]
    # A has a copy landing on B (M1) and C has a copy landing on A's placement? no: C's copy overlaps A (M1 'overlaps gene-A')
    rows.append(row("gene-A", cls="dropped_M1", contig="c2", start0=5100, end=6900, note="overlaps gene-B (c2)", seqid="0.97"))
    rows.append(row("gene-C", cls="dropped_M1", contig="c1", start0=1100, end=2900, note="overlaps gene-A (c1)", seqid="0.99"))
    rows.append(row("gene-D", cls="dropped_M2", contig="c3", start0=9100, end=10900, note="overlaps copy of gene-E (c3)", seqid="0.96"))
    rows.append(row("gene-F", cls="extra_copy", contig="c9", start0=300, end=1800, strand="-", seqid="0.98"))
    for r in rows:
        if r["cls"] in X.COPY_CLASSES:
            r["tag"] = "1"
    return rows, recs


class PoolTests(HermeticCase):
    def setUp(self):
        self.rows, self.recs = mini_selflift()

    def test_copy_rows_join_records_through_the_note_and_unannotated_sites_through_copy_nodes(self):
        pool = X.build_pool(self.rows, self.recs)
        comps = sorted(sorted(c.record_ids) for c in pool.components.values())
        self.assertEqual(comps, [["gene-A", "gene-B", "gene-C"], ["gene-D", "gene-E"], ["gene-F"]])
        f = next(c for c in pool.components.values() if "gene-F" in c.record_ids)
        self.assertEqual(f.sites, ["copy@c9:300-1800-"])

    def test_records_without_a_copy_row_are_not_in_the_pool(self):
        pool = X.build_pool(self.rows, self.recs)
        everything = {r for c in pool.components.values() for r in c.record_ids}
        self.assertNotIn("gene-G", everything)
        self.assertNotIn("gene-H", everything)

    def test_component_size_counts_annotated_records_of_any_length_and_copy_sites(self):
        pool = X.build_pool(self.rows, self.recs)
        f = next(c for c in pool.components.values() if "gene-F" in c.record_ids)
        self.assertEqual(f.size("nodes"), 2)          # gene-F and its copy@ site
        self.assertEqual(f.size("records"), 1)
        abc = next(c for c in pool.components.values() if "gene-A" in c.record_ids)
        self.assertEqual((abc.size("nodes"), abc.size("records")), (3, 3))

    def test_such_records_are_placed_annotated_records_of_at_least_200_bp(self):
        rows, recs = self.rows, self.recs
        recs["gene-S"] = rec("gene-S", "c1", 40000, 40100)              # 100 bp: short
        recs["gene-P"] = rec("gene-P", "c1", 50000, 52000)              # class partial: not placed
        rows.append(row("gene-S", contig="c1", start0=40000, end=40100))
        rows.append(row("gene-P", cls="partial", contig="c1", start0=50000, end=52000))
        rows.append(row("gene-S", cls="dropped_M1", contig="c1", start0=1200, end=1300, note="overlaps gene-A (c1)"))
        rows.append(row("gene-P", cls="dropped_M1", contig="c1", start0=1400, end=1500, note="overlaps gene-A (c1)"))
        for r in rows:
            if r["cls"] in X.COPY_CLASSES:
                r["tag"] = "1"
        pool = X.build_pool(rows, recs)
        abc = next(c for c in pool.components.values() if "gene-A" in c.record_ids)
        self.assertEqual(sorted(r["id"] for r in abc.such()), ["gene-A", "gene-B", "gene-C"])
        self.assertEqual(abc.size("records"), 5)        # the short and the partial record are members

    def test_such_records_are_sorted_by_contig_then_start(self):
        pool = X.build_pool(self.rows, self.recs)
        abc = next(c for c in pool.components.values() if "gene-A" in c.record_ids)
        self.assertEqual([r["id"] for r in abc.such()], ["gene-A", "gene-C", "gene-B"])     # c1:1000, c1:20000, c2:5000

    def test_rank_order_is_contig_as_a_python_string_then_start0_then_end_then_id(self):
        comp = X.Component("PC_t")
        for rid, c, s, e in (("t3", "c1", 100, 900), ("t2", "c1", 100, 500), ("t1", "c1", 100, 500), ("t0", "c10", 5, 900), ("t4", "c2", 1, 900)):
            r = rec(rid, c, s, e)
            r["cls"] = "in_place"
            comp.recs[rid] = r
            comp.record_ids.append(rid)
        # 'c1' < 'c10' < 'c2' as strings; the same contig and start0 are ordered by end, then by id
        self.assertEqual([r["id"] for r in comp.such()], ["t1", "t2", "t3", "t0", "t4"])
        plan = X.rank_plan([comp])
        self.assertEqual([r["records"][0]["id"] for r in plan], ["t1", "t2", "t3", "t0", "t4"])

    def test_a_component_id_is_stable_and_names_its_first_record(self):
        pool = X.build_pool(self.rows, self.recs)
        self.assertEqual(sorted(pool.components), ["PC_gene-A", "PC_gene-D", "PC_gene-F"])

    def test_every_copy_row_is_kept_with_its_partner_and_sequence_id(self):
        pool = X.build_pool(self.rows, self.recs)
        abc = pool.components["PC_gene-A"]
        pairs = sorted((c["source_id"], c["partner"], c["cls"], c["sequence_id"]) for c in abc.copy_rows)
        self.assertEqual(pairs, [("gene-A", "gene-B", "dropped_M1", "0.97"), ("gene-C", "gene-A", "dropped_M1", "0.99")])


class HubAndCalledSetTests(HermeticCase):
    def big_pool(self, n_short=40, n_sites=0):
        """a component of 4 long records plus n_short short records, each short record linked by an M1 copy of gene-L0"""
        recs, rows = {}, []
        for i in range(4):
            r = rec(f"gene-L{i}", "c1", i * 10000, i * 10000 + 2000)
            recs[r["id"]] = r
            rows.append(row(r["id"], contig="c1", start0=r["start0"], end=r["end"]))
        for i in range(n_short):
            r = rec(f"gene-S{i}", "c2", i * 1000, i * 1000 + 100)
            recs[r["id"]] = r
            rows.append(row(r["id"], contig="c2", start0=r["start0"], end=r["end"]))
            rows.append(row(f"gene-S{i}", cls="dropped_M1", contig="c1", start0=100, end=200, note="overlaps gene-L0 (c1)"))
        for i in range(1, 4):
            rows.append(row("gene-L0", cls="dropped_M1", contig="c1", start0=i * 10000 + 5, end=i * 10000 + 1900, note=f"overlaps gene-L{i} (c1)"))
        for i in range(n_sites):
            rows.append(row("gene-L0", cls="extra_copy", contig="c9", start0=i * 3000, end=i * 3000 + 1000))
        for r in rows:
            if r["cls"] in X.COPY_CLASSES:
                r["tag"] = "1"
        return rows, recs

    def test_the_hub_rule_is_more_than_51_members_counting_short_records(self):
        rows, recs = self.big_pool(n_short=48)           # 4 + 48 = 52 annotated records, only 4 of them long
        pool = X.build_pool(rows, recs)
        c = next(iter(pool.components.values()))
        self.assertEqual((c.size("records"), len(c.such())), (52, 4))
        self.assertTrue(c.hub("nodes") and c.hub("records"))
        self.assertEqual(X.called_components(pool, "nodes"), [])

    def test_exactly_51_members_is_not_a_hub(self):
        rows, recs = self.big_pool(n_short=47)           # 51
        pool = X.build_pool(rows, recs)
        c = next(iter(pool.components.values()))
        self.assertEqual(c.size("records"), 51)
        self.assertFalse(c.hub("records"))
        self.assertEqual(len(X.called_components(pool, "records")), 1)

    def test_copy_sites_count_as_members_in_the_literal_node_mode_only(self):
        rows, recs = self.big_pool(n_short=45, n_sites=5)   # 49 records + 5 sites = 54 nodes
        pool = X.build_pool(rows, recs)
        c = next(iter(pool.components.values()))
        self.assertEqual((c.size("records"), c.size("nodes")), (49, 54))
        self.assertTrue(c.hub("nodes"))
        self.assertFalse(c.hub("records"))
        self.assertEqual(len(X.called_components(pool, "nodes")), 0)
        self.assertEqual(len(X.called_components(pool, "records")), 1)

    def test_a_component_with_one_such_record_is_not_called(self):
        rows, recs = mini_selflift()
        pool = X.build_pool(rows, recs)
        called = X.called_components(pool, "nodes")
        self.assertEqual(sorted(c.id for c in called), ["PC_gene-A", "PC_gene-D"])     # F has one record
        self.assertEqual(len(pool.components), 3)


class RankPlanTests(HermeticCase):
    def comps(self):
        rows, recs = mini_selflift()
        recs["gene-A2"] = rec("gene-A2", "c0", 50, 2050)
        rows.append(row("gene-A2", contig="c0", start0=50, end=2050))
        rows.append(row("gene-A2", cls="dropped_M1", contig="c1", start0=1300, end=1600, note="overlaps gene-A (c1)"))
        for r in rows:
            if r["cls"] in X.COPY_CLASSES:
                r["tag"] = "1"
        return X.called_components(X.build_pool(rows, recs), "nodes")

    def test_run_j_holds_the_jth_record_of_every_called_component(self):
        plan = X.rank_plan(self.comps())
        self.assertEqual([len(r["records"]) for r in plan], [2, 2, 1, 1])      # ABC+A2 has 4 records, DE has 2
        self.assertEqual([r["run"] for r in plan], [1, 2, 3, 4])
        by_comp = collections.defaultdict(list)
        for r in plan:
            for rc in r["records"]:
                by_comp[rc["component"]].append((r["run"], rc["id"]))
        self.assertEqual(by_comp["PC_gene-A"], [(1, "gene-A2"), (2, "gene-A"), (3, "gene-C"), (4, "gene-B")])    # c0, c1:1000, c1:20000, c2
        self.assertEqual(by_comp["PC_gene-D"], [(1, "gene-D"), (2, "gene-E")])

    def test_no_run_holds_two_records_of_one_component(self):
        plan = X.rank_plan(self.comps())
        for r in plan:
            comps = [rc["component"] for rc in r["records"]]
            self.assertEqual(len(comps), len(set(comps)))

    def test_records_of_a_run_are_ordered_by_component_anchor(self):
        plan = X.rank_plan(self.comps())
        self.assertEqual([rc["component"] for rc in plan[0]["records"]], ["PC_gene-A", "PC_gene-D"])

    def test_the_number_of_runs_is_the_largest_such_count(self):
        plan = X.rank_plan(self.comps())
        self.assertEqual(len(plan), 4)

    def test_a_run_is_one_shard_when_it_has_at_most_cap_records(self):
        plan = X.rank_plan(self.comps(), cap=2500)
        self.assertTrue(all(len(r["shards"]) == 1 for r in plan))
        self.assertEqual(plan[0]["shards"][0][0], "r01")

    def test_a_run_with_more_than_cap_records_is_cut_at_record_free_positions(self):
        recs = [dict(id=f"g{i}", contig="c1", start0=i * 100, end=i * 100 + 50, idx=i, component=f"P{i}") for i in range(10)]
        recs[4]["end"] = 700                                   # g4 spans g5 and g6 (starts 500, 600): a cut inside is not record-free
        shards = X.plan_run_shards("r01", recs, cap=4)
        by_id = {r["id"]: r for r in recs}
        self.assertEqual(sorted(i for _, ids in shards for i in ids), sorted(by_id))
        self.assertTrue(all(len(ids) <= 4 for _, ids in shards))
        self.assertTrue(all(lab.startswith("r01.") for lab, _ in shards))
        for (_, a), (_, b) in zip(shards, shards[1:]):
            self.assertLessEqual(max(by_id[i]["end"] for i in a), min(by_id[i]["start0"] for i in b))      # nothing spans a cut

    def test_contigs_change_is_a_record_free_cut(self):
        recs = [dict(id=f"a{i}", contig="c1", start0=i * 10, end=i * 10 + 50, idx=i, component=f"P{i}") for i in range(3)]
        recs += [dict(id=f"b{i}", contig="c2", start0=i * 10, end=i * 10 + 50, idx=10 + i, component=f"Q{i}") for i in range(3)]
        cuts = X.record_free_cuts(sorted(recs, key=lambda r: (r["contig"], r["start0"], r["end"], r["id"])))
        self.assertEqual(cuts, [3])           # overlapping a0 a1 a2 on c1 give no cut; the contig change gives one

    def test_split_cuts_nearest_the_middle_and_names_the_halves(self):
        recs = [dict(id=f"g{i}", contig="c1", start0=i * 100, end=i * 100 + 50, idx=i, component=f"P{i}") for i in range(8)]
        by_id = {r["id"]: r for r in recs}
        halves = X.split_shard("r05", [r["id"] for r in recs], by_id)
        self.assertEqual([lab for lab, _ in halves], ["r05.s1", "r05.s2"])
        self.assertEqual([len(ids) for _, ids in halves], [4, 4])

    def test_a_shard_without_a_record_free_cut_cannot_be_split(self):
        recs = [dict(id=f"g{i}", contig="c1", start0=i, end=i + 1000, idx=i, component=f"P{i}") for i in range(5)]
        by_id = {r["id"]: r for r in recs}
        with self.assertRaises(X.NoSplit):
            X.split_shard("r01", [r["id"] for r in recs], by_id)

    def test_summary_of_a_plan(self):
        plan = X.rank_plan(self.comps())
        s = X.plan_summary(plan)
        self.assertEqual((s["J"], s["records"]), (4, 6))
        self.assertEqual(s["records_per_run"], [2, 2, 1, 1])
        self.assertEqual(s["runs_with_at_most_10"], 4)


# ----------------------------------------------------------------------------------------------- group B: the Liftoff command, the GFFs, the run function
REF_GFF = """##gff-version 3
NC_1\tRefSeq\tregion\t1\t100000\t.\t+\t.\tID=NC_1:1..100000;Name=1
NC_1\tRefSeq\tgene\t101\t1100\t.\t+\t.\tID=gene-A;Name=A;gene_biotype=protein_coding
NC_1\tRefSeq\tmRNA\t101\t1100\t.\t+\t.\tID=rna-A1;Parent=gene-A
NC_1\tRefSeq\texon\t101\t500\t.\t+\t.\tID=exon-A1-1;Parent=rna-A1
NC_1\tRefSeq\texon\t801\t1100\t.\t+\t.\tID=exon-A1-2;Parent=rna-A1
NC_1\tRefSeq\tCDS\t150\t500\t.\t+\t.\tID=cds-A1;Parent=rna-A1
NC_1\tRefSeq\tgene\t2001\t3000\t.\t-\t.\tID=gene-B;Name=B;gene_biotype=lncRNA
NC_1\tRefSeq\tlnc_RNA\t2001\t3000\t.\t-\t.\tID=rna-B1;Parent=gene-B
NC_1\tRefSeq\texon\t2001\t3000\t.\t-\t.\tID=exon-B1-1;Parent=rna-B1
NC_1\tRefSeq\tpseudogene\t5001\t6000\t.\t+\t.\tID=gene-P;Name=P;gene_biotype=pseudogene
NC_1\tRefSeq\texon\t5001\t6000\t.\t+\t.\tID=exon-P-1;Parent=gene-P
NC_1\tRefSeq\tcDNA_match\t5001\t6000\t.\t+\t.\tID=cdna-1;Target=x 1 1000
NC_1\tRefSeq\tgene\t8001\t9000\t.\t+\t.\tID=gene-Z;Name=Z;gene_biotype=protein_coding
NC_1\tRefSeq\texon\t8001\t9000\t.\t+\t.\tID=exon-Z-1;Parent=gene-Z
"""


def lifted_gff(rows) -> str:
    """a stub Liftoff output: rows = (id, contig, start1, end1, strand, coverage, seqid, extra_copy_number, partial, low_identity, [exon (s1, e1)...], rtype)"""
    out = ["##gff-version 3"]
    for fid, contig, s, e, st, cov, sid, tag, part, low, exs, rtype in rows:
        out.append(f"{contig}\tLiftoff\t{rtype}\t{s}\t{e}\t.\t{st}\t.\tID={fid};coverage={cov};sequence_ID={sid};extra_copy_number={tag};"
                   f"copy_num_ID={fid}_{tag};partial_mapping={part};low_identity={low}")
        for k, (a, b) in enumerate(exs, 1):
            out.append(f"{contig}\tLiftoff\texon\t{a}\t{b}\t.\t{st}\t.\tID=exon-{fid}-{k};Parent={fid}")
    return "\n".join(out) + "\n"


class CommandTests(HermeticCase):
    def test_the_registered_command_in_the_registered_order(self):
        cmd = LC.liftoff_command("/t/target.fa", "/r/ref.fa", "/m/minimap2", liftoff_bin="/l/liftoff")
        self.assertEqual(cmd, ["/l/liftoff", "/t/target.fa", "/r/ref.fa", "-g", "in.gff", "-f", "types.txt", "-o", "lifted.gff3", "-u", "unmapped.txt",
                               "-dir", "inter", "-copies", "-sc", "0.95", "-a", "0.5", "-s", "0.5", "-overlap", "0.1", "-p", "4", "-m", "/m/minimap2"])

    def test_the_outgroup_is_the_target_and_gorilla_the_reference(self):
        t = LC.OUTGROUPS["H"]
        cmd = LC.liftoff_command(t["genome"], LC.GORILLA_FASTA, "/m/minimap2")
        self.assertEqual(cmd[1], t["genome"])
        self.assertEqual(cmd[2], LC.GORILLA_FASTA)

    def test_types_file_holds_pseudogene_only(self):
        self.assertEqual(LC.types_text(), "pseudogene\n")

    def test_three_outgroups_with_their_existing_indexes(self):
        self.assertEqual(sorted(LC.OUTGROUPS), ["C", "H", "O"])
        for k, t in LC.OUTGROUPS.items():
            self.assertTrue(t["mmi"].endswith(".mmi"))
            self.assertIn(t["species"], ("human", "chimpanzee", "orangutan"))

    def test_registered_parameters_equal_the_self_lift_parameters(self):
        import _liftoff
        self.assertEqual((LC.LIFT_SC, LC.LIFT_A, LC.LIFT_S, LC.LIFT_OVERLAP, LC.LIFT_THREADS, LC.HARD_LIMIT_S),
                         (_liftoff.LIFT_SC, _liftoff.LIFT_A, _liftoff.LIFT_S, _liftoff.LIFT_OVERLAP, _liftoff.LIFT_THREADS, _liftoff.HARD_LIMIT_S))
        self.assertEqual((LC.LIFT_SC, LC.LIFT_OVERLAP, LC.HARD_LIMIT_S), (0.95, 0.1, 590))


class RunGffTests(HermeticCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.gff = os.path.join(self.tmp, "ref.gff")
        with open(self.gff, "w") as fh:
            fh.write(REF_GFF)

    def test_each_run_gets_the_lines_of_its_records_and_nothing_else(self):
        assign = {"gene-A": "r01", "gene-P": "r01", "gene-B": "r02"}
        counts = LC.write_run_gffs(self.gff, assign, os.path.join(self.tmp, "gff"))
        r1 = rd(os.path.join(self.tmp, "gff", "r01", "in.gff")).splitlines()
        r2 = rd(os.path.join(self.tmp, "gff", "r02", "in.gff")).splitlines()
        self.assertEqual(r1[0], "##gff-version 3")
        ids1 = [x.split("ID=")[1].split(";")[0] for x in r1[1:]]
        self.assertEqual(ids1, ["gene-A", "rna-A1", "exon-A1-1", "exon-A1-2", "cds-A1", "gene-P", "exon-P-1"])
        ids2 = [x.split("ID=")[1].split(";")[0] for x in r2[1:]]
        self.assertEqual(ids2, ["gene-B", "rna-B1", "exon-B1-1"])
        self.assertEqual(counts, {"r01": 2, "r02": 1})

    def test_region_and_match_lines_and_unselected_records_are_left_out(self):
        counts = LC.write_run_gffs(self.gff, {"gene-A": "r01"}, os.path.join(self.tmp, "gff"))
        text = rd(os.path.join(self.tmp, "gff", "r01", "in.gff"))
        self.assertNotIn("region", text)
        self.assertNotIn("cDNA_match", text)
        self.assertNotIn("gene-Z", text)

    def test_lines_are_copied_unchanged(self):
        LC.write_run_gffs(self.gff, {"gene-B": "r02"}, os.path.join(self.tmp, "gff"))
        got = rd(os.path.join(self.tmp, "gff", "r02", "in.gff")).splitlines()[1:]
        want = [l for l in REF_GFF.splitlines() if "gene-B" in l or "rna-B1" in l or "exon-B1" in l]
        self.assertEqual(got, want)

    def test_a_record_that_is_not_in_the_annotation_is_an_error(self):
        with self.assertRaises(KeyError):
            LC.write_run_gffs(self.gff, {"gene-A": "r01", "gene-NOPE": "r01"}, os.path.join(self.tmp, "gff"))


LIFTOFF_STUB = r"""#!/usr/bin/env python3
# a stand-in for `liftoff`: copies the prepared GFF3 to -o and, as Liftoff does, lists in -u every parent record of in.gff that it did not lift
import os, re, shutil, sys, time
args = sys.argv[1:]
out = args[args.index("-o") + 1] if "-o" in args else "lifted.gff3"
un = args[args.index("-u") + 1] if "-u" in args else "unmapped.txt"
parents = []
for line in open("in.gff"):
    f = line.rstrip("\n").split("\t")
    if len(f) == 9 and f[2] in ("gene", "pseudogene") and "Parent=" not in f[8]:
        parents.append(re.search(r"ID=([^;]+)", f[8]).group(1))
if len(parents) > int(os.environ.get("STUB_MAX_GENES", "%(max_genes)d")):
    time.sleep(30)
time.sleep(%(sleep)s)
shutil.copy(os.environ["STUB_LIFTED"], out)
lifted = set()
for line in open(out):
    f = line.rstrip("\n").split("\t")
    if len(f) == 9 and f[2] in ("gene", "pseudogene") and "Parent=" not in f[8]:
        lifted.add(re.sub(r"_[0-9]+$", "", re.search(r"ID=([^;]+)", f[8]).group(1)))
mode = os.environ.get("STUB_UNMAPPED", "rest")
if mode == "rest":
    open(un, "w").write("".join(p + "\n" for p in parents if p not in lifted))
elif mode != "none":
    open(un, "w").write("".join(x + "\n" for x in mode.split(",") if x))
sys.exit(%(exit_code)d)
"""


def stub_liftoff(directory, exit_code=0, sleep=0, max_genes=10 ** 6, name="liftoff_stub.py"):
    p = os.path.join(directory, name)
    with open(p, "w") as fh:
        fh.write(LIFTOFF_STUB % dict(exit_code=exit_code, sleep=sleep, max_genes=max_genes))
    os.chmod(p, os.stat(p).st_mode | stat.S_IXUSR)
    return p


class RunShardTests(HermeticCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.lifted = os.path.join(self.tmp, "stub_lifted.gff3")
        with open(self.lifted, "w") as fh:
            fh.write(lifted_gff([("gene-A", "h1", 1001, 2000, "+", "0.98", "0.97", 0, "False", "False", [(1001, 1400), (1701, 2000)], "gene"),
                                 ("gene-A_1", "h2", 501, 1500, "-", "0.96", "0.95", 1, "False", "False", [(501, 1500)], "gene")]))
        self.run_dir = os.path.join(self.tmp, "run")
        os.makedirs(self.run_dir)
        with open(os.path.join(self.run_dir, "in.gff"), "w") as fh:
            fh.write(REF_GFF)

    def test_a_finished_run_leaves_parsed_loci_a_done_key_and_a_time_record(self):
        env = {"STUB_LIFTED": self.lifted}
        res = LC.run_shard(self.run_dir, ["gene-A", "gene-B"], liftoff_bin=stub_liftoff(self.tmp), target="/t.fa", reference="/r.fa", minimap2="/m", hard_s=60,
                           env=env)
        self.assertEqual(res["status"], "done")
        loci = rd(os.path.join(self.run_dir, "loci.tsv")).splitlines()
        self.assertEqual(loci[0].split("\t")[:4], ["feature_id", "copy_num_id", "tag", "rtype"])
        self.assertEqual(len(loci), 3)
        self.assertTrue(os.path.exists(os.path.join(self.run_dir, "done.key")))
        self.assertTrue(os.path.exists(os.path.join(self.run_dir, "lifted.gff3.gz")))
        self.assertFalse(os.path.exists(os.path.join(self.run_dir, "lifted.gff3")))
        self.assertEqual(rd(os.path.join(self.run_dir, "types.txt")), "pseudogene\n")
        self.assertIn("wall_s", rd(os.path.join(self.run_dir, "time.txt")))

    def test_a_done_run_is_skipped_and_a_changed_command_is_redone(self):
        env = {"STUB_LIFTED": self.lifted}
        kw = dict(liftoff_bin=stub_liftoff(self.tmp), target="/t.fa", reference="/r.fa", minimap2="/m", hard_s=60, env=env)
        LC.run_shard(self.run_dir, ["gene-A", "gene-B"], **kw)
        self.assertEqual(LC.run_shard(self.run_dir, ["gene-A", "gene-B"], **kw)["status"], "skipped")
        kw2 = dict(kw, target="/other.fa")
        self.assertEqual(LC.run_shard(self.run_dir, ["gene-A", "gene-B"], **kw2)["status"], "done")

    def test_a_run_that_outlives_the_hard_limit_is_killed_and_recorded(self):
        env = {"STUB_LIFTED": self.lifted}
        res = LC.run_shard(self.run_dir, ["gene-A", "gene-B"], liftoff_bin=stub_liftoff(self.tmp, sleep=30), target="/t.fa", reference="/r.fa", minimap2="/m",
                           hard_s=1, env=env)
        self.assertEqual(res["status"], "toolong")
        self.assertTrue(os.path.exists(os.path.join(self.run_dir, "toolong")))
        self.assertFalse(os.path.exists(os.path.join(self.run_dir, "done.key")))

    def test_a_failing_liftoff_is_a_stop_not_a_retry(self):
        env = {"STUB_LIFTED": self.lifted}
        with self.assertRaises(LC.StopRun):
            LC.run_shard(self.run_dir, ["gene-A", "gene-B"], liftoff_bin=stub_liftoff(self.tmp, exit_code=3), target="/t.fa", reference="/r.fa", minimap2="/m",
                         hard_s=60, env=env)

    def test_peak_memory_over_the_limit_is_a_stop(self):
        LC.check_rss({"peak_rss_kb": 10.2e6}, limit_gb=16)
        with self.assertRaises(LC.StopRun):
            LC.check_rss({"peak_rss_kb": 17.5e6}, limit_gb=16)

    def test_a_run_that_grows_past_the_memory_limit_is_stopped_while_it_runs(self):
        stub = os.path.join(self.tmp, "hog.sh")
        with open(stub, "w") as fh:
            fh.write("#!/bin/bash\npython3 -c \"import time; b = b'x' * (300 * 1024 * 1024); time.sleep(60)\"\n")
        os.chmod(stub, os.stat(stub).st_mode | stat.S_IXUSR)
        with self.assertRaises(LC.StopRun) as cm:                  # hard_s is 120 s and the stub sleeps 60: only the watchdog can end this with a StopRun naming the RSS
            LC.run_shard(self.run_dir, ["gene-A", "gene-B"], liftoff_bin=stub, target="/t.fa", reference="/r.fa", minimap2="/m", hard_s=120, env={},
                         rss_limit_gb=0.15, poll_s=0.2)
        self.assertIn("RSS", str(cm.exception))
        self.assertFalse(os.path.exists(os.path.join(self.run_dir, "done.key")))

    def test_the_memory_watchdog_does_not_disturb_a_normal_run(self):
        env = {"STUB_LIFTED": self.lifted}
        res = LC.run_shard(self.run_dir, ["gene-A", "gene-B"], liftoff_bin=stub_liftoff(self.tmp), target="/t.fa", reference="/r.fa", minimap2="/m", hard_s=60,
                           env=env, rss_limit_gb=16, poll_s=0.2)
        self.assertEqual(res["status"], "done")

    # --- fix round, major 4: Liftoff ignores minimap2's exit status; a run is done only when every parent record of its in.gff is accounted for
    def shard(self, **env):
        return LC.run_shard(self.run_dir, ["gene-A", "gene-B", "gene-P", "gene-Z"], liftoff_bin=stub_liftoff(self.tmp), target="/t.fa", reference="/r.fa", minimap2="/m",
                            hard_s=60, env=dict({"STUB_LIFTED": self.lifted}, **env))

    def test_a_record_that_is_neither_lifted_nor_listed_unmapped_fails_the_run(self):
        with self.assertRaises(LC.StopRun) as cm:
            self.shard(STUB_UNMAPPED="gene-B")                  # gene-P and gene-Z vanished, as after a minimap2 that was killed half way
        msg = str(cm.exception)
        self.assertIn("2 of 4", msg)
        self.assertIn("gene-P", msg)
        self.assertIn("gene-Z", msg)
        self.assertFalse(os.path.exists(os.path.join(self.run_dir, "done.key")))
        self.assertIn("gene-P", rd(os.path.join(self.run_dir, "incomplete")))
        self.assertTrue(os.path.exists(os.path.join(self.run_dir, "lifted.gff3")))               # kept for the post-mortem

    def test_a_missing_unmapped_file_fails_the_run(self):
        with self.assertRaises(LC.StopRun) as cm:
            self.shard(STUB_UNMAPPED="none")
        self.assertIn("unmapped.txt", str(cm.exception))
        self.assertFalse(os.path.exists(os.path.join(self.run_dir, "done.key")))

    def test_a_record_with_only_a_copy_row_is_accounted_for(self):
        with open(self.lifted, "w") as fh:
            fh.write(lifted_gff([("gene-B_1", "h3", 100, 900, "+", "0.8", "0.96", 1, "False", "False", [(100, 900)], "gene")]))
        res = self.shard(STUB_UNMAPPED="gene-A,gene-P,gene-Z")
        self.assertEqual(res["status"], "done")

    def test_a_listed_unmapped_record_is_accounted_for(self):
        self.assertEqual(self.shard()["status"], "done")             # gene-A lifted (with a copy), gene-B, gene-P, gene-Z listed by the stub

    def test_the_run_reports_its_own_elapsed_time_and_the_summed_peak(self):
        res = self.shard()
        self.assertGreater(res["elapsed_s"], 0)
        self.assertIn("peak_summed_kb", res)

    # --- fix round: the files of an earlier attempt do not survive a new one
    def test_stale_files_of_an_earlier_attempt_do_not_survive_a_new_run(self):
        self.shard()
        with open(os.path.join(self.run_dir, "rss_stop"), "w") as fh:
            fh.write("an old stop\n")
        with open(os.path.join(self.run_dir, "time.txt"), "w") as fh:
            fh.write("wall_s\t123.0\npeak_rss_kb\t1\n")
        res = LC.run_shard(self.run_dir, ["gene-A", "gene-B", "gene-P", "gene-Z"], liftoff_bin=stub_liftoff(self.tmp, sleep=30, name="slow.py"), target="/t.fa",
                           reference="/r.fa", minimap2="/m", hard_s=1, env={"STUB_LIFTED": self.lifted})
        self.assertEqual(res["status"], "toolong")
        self.assertIsNone(res["wall_s"], "a killed run does not report the wall of the run before it")
        self.assertFalse(os.path.exists(os.path.join(self.run_dir, "rss_stop")))
        self.assertFalse(os.path.exists(os.path.join(self.run_dir, "loci.tsv")))

    def test_done_key_is_written_whole(self):
        self.shard()
        self.assertFalse(os.path.exists(os.path.join(self.run_dir, "done.key.tmp")))
        json.loads(rd(os.path.join(self.run_dir, "done.key")))


class FrozenPlanTests(HermeticCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        gff = os.path.join(self.tmp, "ref.gff")
        with open(gff, "w") as fh:
            fh.write(REF_GFF)
        self.gff, self.out = gff, os.path.join(self.tmp, "gff")

    def test_the_run_gffs_are_written_once_and_again_only_when_the_assignment_changes(self):
        assign = {"gene-A": "r01", "gene-B": "r02"}
        self.assertTrue(LC.ensure_run_gffs(self.gff, assign, self.out))
        first = os.path.getmtime(os.path.join(self.out, "r01", "in.gff"))
        self.assertFalse(LC.ensure_run_gffs(self.gff, assign, self.out))
        self.assertEqual(os.path.getmtime(os.path.join(self.out, "r01", "in.gff")), first)
        self.assertTrue(LC.ensure_run_gffs(self.gff, {"gene-A": "r01", "gene-B": "r01"}, self.out))
        self.assertIn("gene-B", rd(os.path.join(self.out, "r01", "in.gff")))

    def test_the_plan_digest_depends_on_the_labels_the_records_and_the_hub_reading(self):
        plan = [dict(label="r01", records=[dict(id="a"), dict(id="b")]), dict(label="r02", records=[dict(id="c")])]
        d = LC.plan_digest(plan, "nodes")
        self.assertEqual(d, LC.plan_digest(list(reversed(plan)), "nodes"))
        self.assertNotEqual(d, LC.plan_digest(plan, "records"))
        self.assertNotEqual(d, LC.plan_digest([dict(label="r01", records=[dict(id="a")]), plan[1]], "nodes"))

    def test_the_first_call_freezes_the_plan_and_a_different_plan_later_is_a_stop(self):
        LC.freeze_plan(self.tmp, "abc")
        LC.freeze_plan(self.tmp, "abc")
        with self.assertRaises(LC.StopRun):
            LC.freeze_plan(self.tmp, "abd")
        self.assertEqual(rd(os.path.join(self.tmp, "plan.digest")).strip(), "abc")

    def test_building_with_a_different_plan_than_the_runs_is_a_stop(self):
        LC.freeze_plan(self.tmp, "abc")
        LC.check_plan_frozen(self.tmp, "abc")
        with self.assertRaises(LC.StopRun):
            LC.check_plan_frozen(self.tmp, "other")
        LC.check_plan_frozen(os.path.join(self.tmp, "nowhere"), "anything")        # nothing was run: nothing to compare


def git(repo, *args):
    import subprocess
    return subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false", *args], capture_output=True, text=True,
                          check=True).stdout.strip()


class CodeIdentityTests(HermeticCase):
    """the heavy run may start only from code that is committed after the Amendment, and later calls must run the same code"""

    def setUp(self):
        self.repo = tempfile.mkdtemp()
        git(self.repo, "init", "-q")
        os.makedirs(os.path.join(self.repo, "code"))
        for f in ("a.py", "b.py"):
            with open(os.path.join(self.repo, "code", f), "w") as fh:
                fh.write(f"print('{f}')\n")
        with open(os.path.join(self.repo, "AMENDMENT.md"), "w") as fh:
            fh.write("amendment\n")
        git(self.repo, "add", "-A")
        git(self.repo, "commit", "-q", "-m", "amendment and code")
        self.amend = git(self.repo, "rev-parse", "HEAD")
        self.files = ("code/a.py", "code/b.py")

    def test_clean_tracked_code_after_the_amendment_commit_is_accepted(self):
        self.assertEqual(LC.check_committed_code(self.amend, self.repo, self.files), git(self.repo, "rev-parse", "HEAD"))

    def test_a_modified_file_is_refused(self):
        with open(os.path.join(self.repo, "code", "a.py"), "a") as fh:
            fh.write("# edit\n")
        with self.assertRaises(LC.StopRun) as cm:
            LC.check_committed_code(self.amend, self.repo, self.files)
        self.assertIn("differs from HEAD", str(cm.exception))

    def test_an_untracked_file_is_refused(self):
        with open(os.path.join(self.repo, "code", "c.py"), "w") as fh:
            fh.write("x = 1\n")
        with self.assertRaises(LC.StopRun) as cm:
            LC.check_committed_code(self.amend, self.repo, self.files + ("code/c.py",))
        self.assertIn("not tracked", str(cm.exception))

    def test_an_amendment_commit_that_is_not_in_the_history_is_refused(self):
        with self.assertRaises(LC.StopRun):
            LC.check_committed_code("0" * 40, self.repo, self.files)
        git(self.repo, "checkout", "-q", "-b", "other", "HEAD~0")
        with open(os.path.join(self.repo, "other.txt"), "w") as fh:
            fh.write("x\n")
        git(self.repo, "add", "-A")
        git(self.repo, "commit", "-q", "-m", "other")
        other = git(self.repo, "rev-parse", "HEAD")
        git(self.repo, "checkout", "-q", "-")
        with self.assertRaises(LC.StopRun) as cm:
            LC.check_committed_code(other, self.repo, self.files)
        self.assertIn("does not contain", str(cm.exception))

    def test_the_first_call_records_the_code_and_a_later_call_with_other_code_is_a_stop(self):
        e = tempfile.mkdtemp()
        ident = LC.code_identity(self.repo, self.files)
        self.assertEqual(sorted(ident), sorted(self.files))
        LC.freeze_identity(e, ident, self.amend)
        LC.freeze_identity(e, ident, self.amend)
        with open(os.path.join(self.repo, "code", "a.py"), "a") as fh:
            fh.write("# edit\n")
        with self.assertRaises(LC.StopRun):
            LC.freeze_identity(e, LC.code_identity(self.repo, self.files), self.amend)
        with self.assertRaises(LC.StopRun):
            LC.freeze_identity(e, ident, "f" * 40)

    # --- fix round: a weak code gate
    def test_a_file_hidden_by_assume_unchanged_is_still_refused(self):
        git(self.repo, "update-index", "--assume-unchanged", "code/a.py")
        with open(os.path.join(self.repo, "code", "a.py"), "a") as fh:
            fh.write("# edit that git diff no longer sees\n")
        with self.assertRaises(LC.StopRun) as cm:
            LC.check_committed_code(self.amend, self.repo, self.files)
        self.assertIn("differs from HEAD", str(cm.exception))

    def commit_amendment(self, text):
        with open(os.path.join(self.repo, "AMENDMENT.md"), "w") as fh:
            fh.write(text)
        git(self.repo, "add", "-A")
        git(self.repo, "commit", "-q", "-m", "amendment text")
        return git(self.repo, "rev-parse", "HEAD")

    def test_the_amendment_commit_must_hold_the_amendment_text(self):
        with self.assertRaises(LC.StopRun) as cm:
            LC.check_committed_code(self.amend, self.repo, self.files, amendment_file="AMENDMENT.md", needles={})
        self.assertIn("Amendment 2", str(cm.exception))
        with self.assertRaises(LC.StopRun):
            LC.check_committed_code(self.amend, self.repo, self.files, amendment_file="NO_SUCH_FILE.md", needles={})

    def test_the_amendment_text_must_list_the_sha256_of_every_code_file_and_tool(self):
        ident = LC.code_identity(self.repo, self.files)
        needles = dict(ident, **{"tool liftoff": "ab" * 32})
        head_only = "## Amendment 2\n" + "\n".join(f"{f} {h}" for f, h in ident.items()) + "\n"
        amend = self.commit_amendment(head_only)
        with self.assertRaises(LC.StopRun) as cm:
            LC.check_committed_code(amend, self.repo, self.files, amendment_file="AMENDMENT.md", needles=needles)
        self.assertIn("tool liftoff", str(cm.exception))
        # (fix round 2, G8: the Amendment also names the code commit, "Commit `<40-hex>` holds"; here the code commit is the first one)
        amend = self.commit_amendment(head_only + "liftoff " + "ab" * 32 + f"\nCommit `{self.amend}` holds\n")
        self.assertEqual(LC.check_committed_code(amend, self.repo, self.files, amendment_file="AMENDMENT.md", needles=needles), git(self.repo, "rev-parse", "HEAD"))

    def test_code_committed_after_the_amendment_that_the_amendment_does_not_list_is_refused(self):
        ident = LC.code_identity(self.repo, self.files)
        amend = self.commit_amendment("## Amendment 2\n" + "\n".join(f"{f} {h}" for f, h in ident.items()) + "\n")
        with open(os.path.join(self.repo, "code", "a.py"), "a") as fh:
            fh.write("# edited and committed after the Amendment\n")
        git(self.repo, "add", "-A")
        git(self.repo, "commit", "-q", "-m", "late code edit")
        with self.assertRaises(LC.StopRun) as cm:
            LC.check_committed_code(amend, self.repo, self.files, amendment_file="AMENDMENT.md", needles=LC.code_identity(self.repo, self.files))
        self.assertIn("code/a.py", str(cm.exception))

    # --- fix round 2, G8: the mark, the placeholders and the code commit named by the Amendment
    def amendment_text(self, mark="## Amendment 2", commit=None, extra=""):
        """an Amendment text that satisfies every rule: the mark, the hashes of the code files, 'Commit `<40-hex>` holds' (the first commit holds the code)"""
        ident = LC.code_identity(self.repo, self.files)
        return mark + "\n" + "\n".join(f"{f} {h}" for f, h in ident.items()) + f"\nCommit `{commit or self.amend}` holds\n" + extra

    def check(self, text, needles=None):
        amend = self.commit_amendment(text)
        return LC.check_committed_code(amend, self.repo, self.files, amendment_file="AMENDMENT.md", needles=LC.code_identity(self.repo, self.files) if needles is None else needles)

    def test_an_amendment_that_satisfies_every_rule_is_accepted(self):
        self.assertEqual(self.check(self.amendment_text()), git(self.repo, "rev-parse", "HEAD"))

    def test_the_mark_is_a_heading_at_the_start_of_a_line(self):
        for bad in ("### Amendment 2", "# Amendment 2", "## Amendment II", "## Amendment 20", "text ## Amendment 2", " ## Amendment 2", "##Amendment 2"):
            with self.assertRaises(LC.StopRun, msg=bad) as cm:
                self.check(self.amendment_text(mark=bad))
            self.assertIn("Amendment 2", str(cm.exception), bad)
        self.assertEqual(self.check(self.amendment_text(mark="Some preface\n\n## Amendment 2 (2026-10-07, committed before the first Phase E call)")), git(self.repo, "rev-parse", "HEAD"))

    def test_the_literal_code_commit_placeholder_is_refused_but_the_amendment_commit_cannot_be_filled_in(self):
        # (the coordinator's clarification of G8: refuse ONLY the literal <CODE_COMMIT>)
        with self.assertRaises(LC.StopRun) as cm:
            self.check(self.amendment_text(extra="see <CODE_COMMIT>\n"))
        self.assertIn("<CODE_COMMIT>", str(cm.exception))
        with self.assertRaises(LC.StopRun):
            self.check(self.amendment_text(extra="Commit `<CODE_COMMIT>` holds\n"))
        # a commit cannot contain its own hash: the command block of the committed Amendment keeps <AMENDMENT_COMMIT>; other angle-bracket words are not refused either
        ok = self.amendment_text(extra="run --amendment-commit <AMENDMENT_COMMIT> --budget-s 530\nE_ROOT/lift/<species>/runs/<label>/ <FOO_COMMIT> <reason>\n")
        self.assertEqual(self.check(ok), git(self.repo, "rev-parse", "HEAD"))

    def test_the_amendment_must_name_the_code_commit(self):
        ident = LC.code_identity(self.repo, self.files)
        text = "## Amendment 2\n" + "\n".join(f"{f} {h}" for f, h in ident.items()) + "\n"
        for body in (text, text + "Commit <CODE> holds\n", text + f"The commit {self.amend[:12]} holds\n", text + f"Commit `{self.amend[:12]}` holds\n",
                     text + f"Commit `{self.amend.upper()}` holds\n"):
            with self.assertRaises(LC.StopRun, msg=body[-60:]) as cm:
                self.check(body)
            self.assertIn("Commit `<40-hex>` holds", str(cm.exception))

    def test_the_code_commit_named_must_exist(self):
        with self.assertRaises(LC.StopRun) as cm:
            self.check(self.amendment_text(commit="0" * 40))
        self.assertIn("is not a commit", str(cm.exception))
        self.assertIn("0" * 12, str(cm.exception))

    def test_the_code_commit_named_must_hold_the_listed_code_files_with_the_listed_hashes(self):
        stale = self.amend                                                     # the first commit: code/a.py as it was before the edit below
        with open(os.path.join(self.repo, "code", "a.py"), "a") as fh:
            fh.write("# an edit made after the stale commit\n")
        git(self.repo, "add", "-A")
        git(self.repo, "commit", "-q", "-m", "code edit")
        code_commit = git(self.repo, "rev-parse", "HEAD")
        with self.assertRaises(LC.StopRun) as cm:
            self.check(self.amendment_text(commit=stale))                      # the hashes listed are those of the edited a.py, which the stale commit does not hold
        self.assertIn("code/a.py", str(cm.exception))
        self.assertIn("another sha256", str(cm.exception))
        self.assertEqual(self.check(self.amendment_text(commit=code_commit)), git(self.repo, "rev-parse", "HEAD"))

    def test_a_code_commit_that_lacks_a_listed_file_is_refused(self):
        empty_tree = git(self.repo, "hash-object", "-t", "tree", "-w", "/dev/null")
        bare = git(self.repo, "commit-tree", empty_tree, "-m", "a commit with no files")              # plumbing: HEAD and the work tree are not touched
        with self.assertRaises(LC.StopRun) as cm:
            self.check(self.amendment_text(commit=bare))
        self.assertIn("does not contain code/a.py", str(cm.exception))

    def test_a_code_file_without_a_needle_is_not_compared_with_the_code_commit(self):
        with open(os.path.join(self.repo, "code", "b.py"), "a") as fh:
            fh.write("# an edit the Amendment does not list\n")
        git(self.repo, "add", "-A")
        git(self.repo, "commit", "-q", "-m", "edit b")
        only_a = {"code/a.py": LC.code_identity(self.repo, self.files)["code/a.py"]}
        text = f"## Amendment 2\n{only_a['code/a.py']}\nCommit `{self.amend}` holds\n"        # the first commit holds an older b.py: there is no needle for it
        self.assertEqual(self.check(text, needles=only_a), git(self.repo, "rev-parse", "HEAD"))

    def test_the_plan_digest_and_the_inputs_are_needles_like_the_code_hashes(self):
        digest, small = "d" * 64, "5" * 64
        needles = dict(LC.code_identity(self.repo, self.files), **{"plan digest": digest, "input gorilla_fai": small})
        with self.assertRaises(LC.StopRun) as cm:
            self.check(self.amendment_text(extra=f"Plan digest `{digest}`\n"), needles=needles)
        self.assertIn("input gorilla_fai", str(cm.exception))
        self.assertEqual(self.check(self.amendment_text(extra=f"Plan digest `{digest}`\nsmall `{small}`\n"), needles=needles), git(self.repo, "rev-parse", "HEAD"))
        with self.assertRaises(LC.StopRun) as cm:
            self.check(self.amendment_text(extra=f"small `{small}`\n"), needles=needles)
        self.assertIn("plan digest", str(cm.exception))

    # --- fix round 2, G4: the commit is resolved to its 40-hex value
    def test_an_abbreviated_amendment_commit_is_resolved_to_40_hex(self):
        self.assertEqual(LC.resolve_commit(self.amend[:7], self.repo), self.amend)
        self.assertEqual(LC.resolve_commit(self.amend[:12], self.repo), self.amend)
        self.assertEqual(LC.resolve_commit(self.amend, self.repo), self.amend)

    def test_a_name_that_is_not_a_commit_is_a_stop(self):
        for bad in ("0" * 40, "abcdef1", "deadbeef"):
            with self.assertRaises(LC.StopRun, msg=bad) as cm:
                LC.resolve_commit(bad, self.repo)
            self.assertIn(bad, str(cm.exception))
        tree = git(self.repo, "rev-parse", "HEAD^{tree}")
        with self.assertRaises(LC.StopRun):
            LC.resolve_commit(tree[:12], self.repo)                              # a tree is not a commit

    def test_the_frozen_code_set_covers_every_file_the_run_imports_or_executes(self):
        for f in ("bench/hierarchy/expansion_table.py", "bench/hierarchy/liftoff_cross.py", "bench/hierarchy/test_expansion_table.py", "bench/hierarchy/test_supplement_e.py",
                  "bench/hierarchy/test_supplement_e2.py", "figures/_liftoff.py", "figures/figlib.py", "tools/rlock.sh"):
            self.assertIn(f, LC.CODE_FILES)
            self.assertTrue(os.path.isfile(os.path.join(str(LC.REPO), f)), f)
        self.assertEqual(len(LC.CODE_FILES), 8)                                  # fix round 2: the second mutation review's tests joined the frozen set

    def test_the_code_commit_check_covers_each_of_the_default_code_files(self):
        # a throw-away repository with the default layout; the code commit named by the Amendment is made stale in ONE file at a time (that file is edited and committed after it)
        repo = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, repo, True)
        git(repo, "init", "-q")
        for f in LC.CODE_FILES + ("AMENDMENT.md",):
            os.makedirs(os.path.dirname(os.path.join(repo, f)) or repo, exist_ok=True)
            with open(os.path.join(repo, f), "w") as fh:
                fh.write(f"# {f}\n")
        git(repo, "add", "-A")
        git(repo, "commit", "-q", "-m", "code")

        def amend(code_commit):
            with open(os.path.join(repo, "AMENDMENT.md"), "w") as fh:
                fh.write("## Amendment 2\n" + "\n".join(f"{f} {h}" for f, h in LC.code_identity(repo).items()) + f"\nCommit `{code_commit}` holds\n")
            git(repo, "add", "-A")
            git(repo, "commit", "-q", "-m", "amendment")
            return git(repo, "rev-parse", "HEAD")
        for f in LC.CODE_FILES:
            before = git(repo, "rev-parse", "HEAD")
            with open(os.path.join(repo, f), "a") as fh:
                fh.write("# edited after the commit the Amendment will name\n")
            git(repo, "add", "-A")
            git(repo, "commit", "-q", "-m", f"edit {f}")
            after = git(repo, "rev-parse", "HEAD")
            with self.assertRaises(LC.StopRun, msg=f) as cm:
                LC.check_committed_code(amend(before), repo, amendment_file="AMENDMENT.md", needles=LC.code_identity(repo))
            self.assertIn(f, str(cm.exception))
            self.assertIn("another sha256", str(cm.exception))
            self.assertEqual(LC.check_committed_code(amend(after), repo, amendment_file="AMENDMENT.md", needles=LC.code_identity(repo)), git(repo, "rev-parse", "HEAD"))


class ScheduleTests(HermeticCase):
    def plan(self):
        recs = {c: [dict(id=f"{c}{i}", component=c, contig="c1", start0=i * 100, end=i * 100 + 50, idx=0, rtype="gene", biotype="x", strand="+",
                         exonic_len=50, rank=i + 1) for i in range(n)] for c, n in (("A", 5), ("B", 3))}
        plan = []
        for j in range(1, 6):
            rs = [recs[c][j - 1] for c in ("A", "B") if len(recs[c]) >= j]
            plan.append(dict(run=j, label=f"r{j:02d}", records=rs, shards=[(f"r{j:02d}", [r["id"] for r in rs])]))
        return plan

    def small(self, prefix, n_runs):
        return [dict(run=j, label=f"{prefix}{j}", records=[dict(id=f"{prefix}{j}x", component=f"{prefix}", contig="c1", start0=j * 100, end=j * 100 + 50,
                                                                  idx=0, rtype="gene", biotype="x", strand="+", exonic_len=50, rank=j)],
                     shards=[(f"{prefix}{j}", [f"{prefix}{j}x"])]) for j in range(1, n_runs + 1)]

    def test_the_two_controls_come_first_then_the_verdict_set_each_outgroup_by_outgroup(self):
        units = LC.schedule_units(self.plan(), self.small("pos", 2), self.small("ctl", 2))
        self.assertEqual(units[:6], [("H", "pos1"), ("H", "pos2"), ("C", "pos1"), ("C", "pos2"), ("O", "pos1"), ("O", "pos2")])
        self.assertEqual(units[6:12], [("H", "ctl1"), ("H", "ctl2"), ("C", "ctl1"), ("C", "ctl2"), ("O", "ctl1"), ("O", "ctl2")])
        self.assertEqual(units[12:15], [("H", "r01"), ("H", "r02"), ("H", "r03")])
        self.assertEqual([u for u in units if u[1] == "r05"], [("H", "r05"), ("C", "r05"), ("O", "r05")])
        self.assertEqual(len(units), 3 * 2 + 3 * 2 + 3 * 5)

    def test_every_label_has_one_group_and_no_verdict_unit_precedes_a_control_unit(self):
        groups = LC.group_map(self.plan(), self.small("pos", 2), self.small("ctl", 2))
        self.assertEqual(groups["pos1"], "pos")
        self.assertEqual(groups["ctl2"], "neg")
        self.assertEqual(groups["r03"], "verdict")
        units = LC.schedule_units(self.plan(), self.small("pos", 2), self.small("ctl", 2))
        order = [groups[lab] for _, lab in units]
        self.assertEqual(order, sorted(order, key=("pos", "neg", "verdict").index))

    def test_a_call_never_mixes_groups_even_when_the_units_would_fit(self):
        units = [("H", "pos1"), ("H", "ctl1"), ("H", "r01")]
        groups = {"pos1": "pos", "ctl1": "neg", "r01": "verdict"}
        pred = {u: 10 for u in units}
        self.assertEqual(LC.plan_calls(units, pred, budget_s=500), [units])
        self.assertEqual(LC.plan_calls(units, pred, budget_s=500, group_of=groups), [[u] for u in units])

    def test_a_call_starts_a_unit_only_if_its_predicted_time_fits_the_budget(self):
        calls = LC.plan_calls([("H", "r01"), ("H", "r02"), ("H", "r03")], {("H", "r01"): 300, ("H", "r02"): 200, ("H", "r03"): 150}, budget_s=500)
        self.assertEqual(calls, [[("H", "r01"), ("H", "r02")], [("H", "r03")]])

    def test_the_first_unit_of_a_call_always_runs(self):
        calls = LC.plan_calls([("H", "r01"), ("H", "r02")], {("H", "r01"): 900, ("H", "r02"): 100}, budget_s=500)
        self.assertEqual(calls, [[("H", "r01")], [("H", "r02")]])

    def test_predicted_wall_is_start_up_plus_per_record_time(self):
        self.assertEqual(LC.predict_wall(100, per_record_s=0.1, startup_s=20, margin=1.0), 30.0)
        self.assertEqual(LC.predict_wall(100, per_record_s=0.1, startup_s=20, margin=1.5), 20 + 10 * 1.5)


# ----------------------------------------------------------------------------------------------- group C: the E-table
def lrow(fid, tag=0, contig="h1", s0=1000, e=3000, cov="1.0", sid="0.98", partial=0, low=0, exons=None, strand="+", rtype="gene"):
    """a parent-level row of a run's loci.tsv (figures/_liftoff.py::parse_lifted)"""
    exons = exons if exons is not None else f"{s0}-{e}"
    el = sum(int(b) - int(a) for a, b in (x.split("-") for x in exons.split(",")))
    return dict(feature_id=fid, copy_num_id=f"{fid}_{tag}", tag=str(tag), rtype=rtype, contig=contig, start0=str(s0), end=str(e), strand=strand,
                coverage=str(cov), sequence_id=str(sid), partial=str(partial), low_identity=str(low), exons=exons, exonic_len=str(el))


def outgroup_rows(*records):
    """records = (id, contig, start0, end, [copy_classes]) -> rows of an outgroup's own liftoff_loci.tsv"""
    rows = []
    for rid, contig, s0, e, copies in records:
        rows.append(row(rid, contig=contig, start0=s0, end=e))
        for cls in copies:
            rows.append(row(rid, cls=cls, contig="x9", start0=0, end=500, note="-" if cls == "extra_copy" else "overlaps gene-Q (x)"))
    return rows


def make_indexes(h_extra=(), c_extra=(), o_extra=()):
    base = {"H": [("gene-H1", "h1", 1000, 3000, []), ("gene-H2", "h2", 1000, 3000, []), ("gene-H3", "h3", 1000, 3000, ["extra_copy"])],
            "C": [("gene-C1", "c1", 1000, 3000, []), ("gene-C2", "c2", 1000, 3000, []), ("gene-C3", "c3", 1000, 3000, [])],
            "O": [("gene-O1", "o1", 1000, 3000, []), ("gene-O2", "o2", 1000, 3000, []), ("gene-O3", "o3", 1000, 3000, [])]}
    for k, extra in (("H", h_extra), ("C", c_extra), ("O", o_extra)):
        base[k] = base[k] + list(extra)
    return {k: X.OutgroupIndex(outgroup_rows(*v), k) for k, v in base.items()}


def lifted_for(rid, h, c, o):
    """lifted rows of one gorilla record in the three outgroups: h, c, o = (contig, s0, e) of its primary placement, or None for no row"""
    out = {}
    for k, loc in (("H", h), ("C", c), ("O", o)):
        out[k] = [] if loc is None else [lrow(rid, 0, *loc)]
    return out


class PlacementRuleTests(HermeticCase):
    def test_a_placement_needs_coverage_and_identity_of_one_half_and_no_flag(self):
        self.assertTrue(X.placement_ok(lrow("g", cov=0.5, sid=0.5)))
        self.assertFalse(X.placement_ok(lrow("g", cov=0.49, sid=0.9)))
        self.assertFalse(X.placement_ok(lrow("g", cov=0.9, sid=0.49)))
        self.assertFalse(X.placement_ok(lrow("g", cov=0.9, sid=0.9, partial=1)))
        self.assertFalse(X.placement_ok(lrow("g", cov=0.9, sid=0.9, low=1)))
        self.assertFalse(X.placement_ok(lrow("g", cov="NA", sid="NA")))

    def test_a_copy_row_is_never_the_placement(self):
        self.assertFalse(X.placement_ok(lrow("g_1", tag=1, cov=1.0, sid=1.0)))


class OutgroupLocusTests(HermeticCase):
    def setUp(self):
        self.idx = make_indexes()

    def test_exactly_half_of_the_shorter_exon_union_is_enough(self):
        m = self.idx["H"].match("h1", [(1000, 1500)])             # shorter union is the placement (500 bp), shares 500 of 500
        self.assertEqual((m["status"], m["locus"]), ("ok", "gene-H1"))
        m = self.idx["H"].match("h1", [(2000, 4000)])             # shares 1000 of min(2000, 2000) = 0.5
        self.assertEqual((m["status"], m["locus"]), ("ok", "gene-H1"))

    def test_a_hair_below_half_is_no_locus(self):
        m = self.idx["H"].match("h1", [(2001, 4001)])             # 999 of 2000
        self.assertEqual(m["status"], "no_outgroup_locus")

    def test_the_shorter_union_is_the_denominator_even_when_it_is_the_outgroup_record(self):
        idx = make_indexes(h_extra=[("gene-Hsmall", "h5", 1000, 1100, [])])
        m = idx["H"].match("h5", [(500, 1100)])                   # the 100 bp record lies wholly inside the 600 bp placement
        self.assertEqual((m["status"], m["locus"]), ("ok", "gene-Hsmall"))

    def test_another_contig_is_no_locus(self):
        self.assertEqual(self.idx["H"].match("h9", [(1000, 3000)])["status"], "no_outgroup_locus")

    def test_a_locus_that_has_a_copy_row_is_an_outgroup_copy_row(self):
        m = self.idx["H"].match("h3", [(1000, 3000)])
        self.assertEqual((m["status"], m["locus"]), ("outgroup_copy_row", "gene-H3"))

    def test_two_single_copy_candidates_are_ambiguous(self):
        idx = make_indexes(h_extra=[("gene-Hdup", "h1", 1500, 3500, [])])
        m = idx["H"].match("h1", [(1500, 3000)])
        self.assertEqual(m["status"], "ambiguous")
        self.assertEqual(sorted(m["candidates"]), ["gene-H1", "gene-Hdup"])

    def test_only_in_place_records_are_candidates(self):
        rows = outgroup_rows(("gene-M", "h1", 1000, 3000, []))
        rows[0]["cls"] = "moved"
        rows += [row("gene-P", cls="partial", contig="h1", start0=1000, end=3000)]
        idx = X.OutgroupIndex(rows, "H")
        self.assertEqual(idx.match("h1", [(1000, 3000)])["status"], "no_outgroup_locus")

    def test_exons_not_spans_decide_the_overlap(self):
        rows = outgroup_rows(("gene-X", "h1", 0, 10000, []))
        rows[0]["exons"] = "0-500,9500-10000"
        idx = X.OutgroupIndex(rows, "H")
        self.assertEqual(idx.match("h1", [(4000, 6000)])["status"], "no_outgroup_locus")      # inside the span, between the exons
        self.assertEqual(idx.match("h1", [(0, 500)])["locus"], "gene-X")

    def test_strict_single_copy_also_refuses_landing_sites_of_other_copies(self):
        rows = outgroup_rows(("gene-L", "h1", 1000, 3000, []), ("gene-S", "h2", 1000, 3000, ["dropped_M1"]))
        rows[-1]["note"] = "overlaps gene-L (x)"
        loose = X.OutgroupIndex(rows, "H")
        strict = X.OutgroupIndex(rows, "H", strict_single_copy=True)
        self.assertEqual(loose.match("h1", [(1000, 3000)])["status"], "ok")
        self.assertEqual(strict.match("h1", [(1000, 3000)])["status"], "outgroup_copy_row")


class RecordCallTests(HermeticCase):
    def setUp(self):
        self.idx = make_indexes()

    def call(self, rid, h, c, o, extra_rows=None, unmapped=None):
        lifted = {k: {rid: v} for k, v in lifted_for(rid, h, c, o).items()}
        if extra_rows:
            for k, rs in extra_rows.items():
                lifted[k][rid] = lifted[k].get(rid, []) + rs
        return X.classify_record(rid, lifted, unmapped or {"H": set(), "C": set(), "O": set()}, self.idx)

    def test_a_record_placed_on_one_single_copy_locus_in_each_outgroup_has_a_triple(self):
        c = self.call("g1", ("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000))
        self.assertEqual(c["triple"], ("gene-H1", "gene-C1", "gene-O1"))
        self.assertIsNone(c["reason"])
        self.assertEqual([c["per"][k]["n_place"] for k in "HCO"], [1, 1, 1])

    def test_no_row_and_an_unmapped_listing_are_no_placement(self):
        c = self.call("g1", None, ("c1", 1000, 3000), ("o1", 1000, 3000))
        self.assertEqual((c["triple"], c["reason"], c["per"]["H"]["detail"]), (None, "no_placement", "no_row"))
        c = self.call("g1", None, ("c1", 1000, 3000), ("o1", 1000, 3000), unmapped={"H": {"g1"}, "C": set(), "O": set()})
        self.assertEqual(c["per"]["H"]["detail"], "unmapped")

    def test_a_partial_only_row_is_no_placement_and_is_counted_as_partial_only(self):
        lifted = {k: {"g1": v} for k, v in lifted_for("g1", ("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000)).items()}
        lifted["H"]["g1"] = [lrow("g1", 0, "h1", 1000, 3000, cov=0.4, sid=0.9, partial=1)]
        c = X.classify_record("g1", lifted, {"H": set(), "C": set(), "O": set()}, self.idx)
        self.assertEqual((c["reason"], c["per"]["H"]["detail"]), ("no_placement", "partial_only"))

    def test_a_flagged_or_weak_primary_row_among_other_rows_is_not_a_placement_but_not_partial_only(self):
        lifted = {k: {"g1": v} for k, v in lifted_for("g1", ("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000)).items()}
        none = {"H": set(), "C": set(), "O": set()}
        lifted["H"]["g1"] = [lrow("g1", 0, "h1", 1000, 3000, cov=0.9, sid=0.9, low=1)]
        self.assertEqual(X.classify_record("g1", lifted, none, self.idx)["per"]["H"]["detail"], "not_a_placement")
        lifted["H"]["g1"] = [lrow("g1", 0, "h1", 1000, 3000, cov=0.4, sid=0.9, partial=1), lrow("g1_1", 1, "h9", 0, 500, cov=0.2, sid=0.96)]
        self.assertEqual(X.classify_record("g1", lifted, none, self.idx)["per"]["H"]["detail"], "not_a_placement")

    # --- fix round, D2 revised (the registered rule, prereg section 7 step 3: a class with a member that has more than one placement is dropped and counted).
    # (replaces test_a_copy_row_of_the_reverse_lift_disqualifies_whatever_its_coverage, which pinned the member-level exclusion: a reverse copy row gave the
    # record no triple at all, so the registered class-level drop could never fire)
    def test_a_record_with_a_reverse_copy_row_keeps_the_locus_of_its_primary_placement_flagged_multiplicity(self):
        extra = {"C": [lrow("g1_1", 1, "c9", 0, 800, cov=0.2, sid=0.96)]}                          # whatever its coverage
        c = self.call("g1", ("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000), extra_rows=extra)
        self.assertEqual((c["triple"], c["reason"]), (("gene-H1", "gene-C1", "gene-O1"), None))
        self.assertEqual((c["per"]["C"]["status"], c["per"]["C"]["locus"], c["per"]["C"]["detail"]), ("ok", "gene-C1", "reverse_copy_row"))
        self.assertEqual((c["per"]["C"]["n_place"], c["per"]["C"]["copy_rows"]), (1, 1))
        self.assertEqual([c["per"][k]["multiplicity"] for k in "HCO"], [False, True, False])
        self.assertEqual(c["multiplicity"], ["C"])

    def test_a_record_without_a_reverse_copy_row_is_not_flagged(self):
        c = self.call("g1", ("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000))
        self.assertEqual((c["multiplicity"], [c["per"][k]["multiplicity"] for k in "HCO"]), ([], [False, False, False]))

    def test_a_reverse_copy_row_does_not_hide_the_copy_row_of_the_outgroups_own_locus(self):
        extra = {"H": [lrow("g1_1", 1, "h9", 0, 800, cov=0.9, sid=0.97)]}
        c = self.call("g1", ("h3", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000), extra_rows=extra)        # h3 is the locus that has a copy row in the outgroup's own lift
        self.assertEqual((c["triple"], c["reason"], c["per"]["H"]["status"], c["per"]["H"]["detail"]), (None, "outgroup_copy_row", "outgroup_copy_row", "outgroup_locus_has_copy_row"))
        self.assertTrue(c["per"]["H"]["multiplicity"])

    def test_a_reverse_copy_row_of_a_record_with_no_placement_changes_nothing(self):
        lifted = {k: {"g1": v} for k, v in lifted_for("g1", ("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000)).items()}
        lifted["H"]["g1"] = [lrow("g1", 0, "h1", 1000, 3000, cov=0.4, sid=0.9, partial=1), lrow("g1_1", 1, "h9", 0, 800, cov=0.9, sid=0.97)]
        c = X.classify_record("g1", lifted, {"H": set(), "C": set(), "O": set()}, self.idx)
        self.assertEqual((c["triple"], c["reason"], c["per"]["H"]["status"], c["per"]["H"]["detail"]), (None, "no_placement", "no_placement", "not_a_placement"))
        self.assertFalse(c["per"]["H"]["multiplicity"])
        lifted["H"]["g1"] = [lrow("g1", 0, "h1", 1000, 3000, cov=0.4, sid=0.9, partial=1)]               # the only row is partial: partial_only, as before
        c = X.classify_record("g1", lifted, {"H": set(), "C": set(), "O": set()}, self.idx)
        self.assertEqual(c["per"]["H"]["detail"], "partial_only")

    def test_a_placement_on_an_outgroup_locus_with_a_copy_row_is_an_outgroup_copy_row(self):
        c = self.call("g1", ("h3", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000))
        self.assertEqual((c["triple"], c["reason"]), (None, "outgroup_copy_row"))

    def test_a_placement_on_no_annotated_locus_is_no_outgroup_locus(self):
        c = self.call("g1", ("h1", 1000, 3000), ("c1", 5000, 6000), ("o1", 1000, 3000))
        self.assertEqual((c["triple"], c["reason"]), (None, "no_outgroup_locus"))

    def test_the_reason_precedence_is_copy_row_then_no_placement_then_no_locus(self):
        c = self.call("g1", ("h3", 1000, 3000), None, ("o9", 5000, 6000))
        self.assertEqual(c["reason"], "outgroup_copy_row")
        c = self.call("g1", None, ("c9", 1000, 3000), ("o1", 1000, 3000))
        self.assertEqual(c["reason"], "no_placement")

    def test_qualifying_placements_are_counted_with_coverage_alone(self):
        extra = {"H": [lrow("g1_1", 1, "h2", 1000, 3000, cov=0.7, sid=0.5), lrow("g1_2", 2, "h2", 5000, 6000, cov=0.3, sid=0.99)]}
        c = self.call("g1", ("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000), extra_rows=extra)
        self.assertEqual(c["per"]["H"]["n_place"], 2)         # the primary and the copy with coverage 0.7; the 0.3 row does not count


def calls_for(spec):
    """spec: {record: triple or None, ...} -> record calls shaped like classify_record's (reason / n_place set by hand)"""
    out = {}
    for rid, t in spec.items():
        out[rid] = dict(triple=t, reason=None if t else "no_placement", per={k: dict(n_place=1, detail="") for k in "HCO"})
    return out


def recs_comp(cid, ids, contig="c1"):
    comp = X.Component(cid)
    for i, rid in enumerate(ids):
        r = rec(rid, contig, i * 10000, i * 10000 + 2000)
        r["cls"] = "in_place"
        comp.recs[rid] = r
        comp.record_ids.append(rid)
    return comp


T1 = ("gene-H1", "gene-C1", "gene-O1")
T2 = ("gene-H2", "gene-C2", "gene-O2")


class ExpansionTests(HermeticCase):
    def test_records_with_the_same_triple_are_one_expansion(self):
        comp = recs_comp("PC_a", ["a1", "a2", "a3"])
        res = X.call_expansions([comp], calls_for({"a1": T1, "a2": T1, "a3": None}))
        self.assertEqual(len(res["expansions"]), 1)
        ex = res["expansions"][0]
        self.assertEqual((ex["component"], ex["triple"], [m["id"] for m in ex["members"]]), ("PC_a", T1, ["a1", "a2"]))

    def test_one_record_per_triple_is_not_an_expansion(self):
        comp = recs_comp("PC_a", ["a1", "a2"])
        res = X.call_expansions([comp], calls_for({"a1": T1, "a2": T2}))
        self.assertEqual(res["expansions"], [])
        self.assertEqual(res["components"]["PC_a"]["different_loci"], 2)

    def test_a_component_with_two_triples_yields_two_expansions_without_chaining(self):
        comp = recs_comp("PC_a", ["a1", "a2", "a3", "a4"])
        res = X.call_expansions([comp], calls_for({"a1": T1, "a2": T1, "a3": T2, "a4": T2}))
        self.assertEqual(sorted(e["triple"] for e in res["expansions"]), sorted([T1, T2]))
        self.assertEqual(res["components"]["PC_a"]["different_loci"], 2)

    def test_expansion_ids_are_sequential_in_component_then_member_order(self):
        res = X.call_expansions([recs_comp("PC_b", ["b1", "b2"]), recs_comp("PC_a", ["a1", "a2"])], calls_for({"a1": T1, "a2": T1, "b1": T2, "b2": T2}))
        self.assertEqual([(e["id"], e["component"]) for e in res["expansions"]], [("EX0001", "PC_a"), ("EX0002", "PC_b")])

    def test_components_with_the_same_triple_are_siblings_and_stay_separate(self):
        res = X.call_expansions([recs_comp("PC_a", ["a1", "a2"]), recs_comp("PC_b", ["b1", "b2"])], calls_for({"a1": T1, "a2": T1, "b1": T1, "b2": T1}))
        self.assertEqual(len(res["expansions"]), 2)
        self.assertEqual(res["siblings"], [dict(triple=T1, expansions=["EX0001", "EX0002"])])

    def test_a_class_with_a_member_of_more_than_one_qualifying_placement_is_dropped_and_counted(self):
        comp = recs_comp("PC_a", ["a1", "a2"])
        calls = calls_for({"a1": T1, "a2": T1})
        calls["a2"]["per"]["O"]["n_place"] = 2
        res = X.call_expansions([comp], calls)
        self.assertEqual(res["expansions"], [])
        self.assertEqual(res["dropped_classes"], [dict(component="PC_a", triple=T1, members=["a1", "a2"], reason="multiplicity")])

    def test_three_records_on_one_locus_one_with_a_reverse_copy_row_the_class_is_dropped_and_counted(self):
        # the plan reviewer's scenario: three records of one component place on human H, the third also has a reverse copy row at H'; the registered rule drops the
        # CLASS (the old code kept the other two as an expansion and dropped only the member)
        comp = recs_comp("PC_a", ["a1", "a2", "a3"])
        idx = make_indexes()
        none = {"H": set(), "C": set(), "O": set()}
        calls = {}
        for rid in ("a1", "a2", "a3"):
            lifted = {k: {rid: v} for k, v in lifted_for(rid, ("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000)).items()}
            if rid == "a3":
                lifted["H"][rid] = lifted["H"][rid] + [lrow("a3_1", 1, "h9", 0, 800, cov=0.97, sid=0.97)]            # H' : a second place in the human genome
            calls[rid] = X.classify_record(rid, lifted, none, idx)
        self.assertEqual([calls[r]["triple"] for r in ("a1", "a2", "a3")], [T1, T1, T1])
        res = X.call_expansions([comp], calls)
        self.assertEqual(res["expansions"], [])
        self.assertEqual(res["dropped_classes"], [dict(component="PC_a", triple=T1, members=["a1", "a2", "a3"], reason="multiplicity")])
        self.assertEqual(res["reasons"], {})                                                    # nobody lacks a triple
        self.assertEqual(res["components"]["PC_a"], dict(n_such=3, different_loci=1, contributes=False))
        self.assertEqual(res["records_in_expansions"], 0)

    def test_the_same_three_records_without_the_copy_row_are_one_expansion(self):
        comp = recs_comp("PC_a", ["a1", "a2", "a3"])
        idx = make_indexes()
        none = {"H": set(), "C": set(), "O": set()}
        calls = {rid: X.classify_record(rid, {k: {rid: v} for k, v in lifted_for(rid, ("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000)).items()}, none, idx)
                 for rid in ("a1", "a2", "a3")}
        res = X.call_expansions([comp], calls)
        self.assertEqual([[m["id"] for m in e["members"]] for e in res["expansions"]], [["a1", "a2", "a3"]])
        self.assertEqual(res["dropped_classes"], [])

    def test_a_flagged_record_alone_in_its_class_is_not_a_dropped_class(self):
        comp = recs_comp("PC_a", ["a1", "a2"])
        calls = calls_for({"a1": T1, "a2": T2})
        calls["a2"]["per"]["H"]["multiplicity"] = True
        res = X.call_expansions([comp], calls)
        self.assertEqual((res["expansions"], res["dropped_classes"]), ([], []))                 # a class of one is no expansion candidate; nothing is dropped
        self.assertEqual(res["components"]["PC_a"]["different_loci"], 2)

    def test_the_flag_alone_drops_the_class_whatever_n_place_says(self):
        comp = recs_comp("PC_a", ["a1", "a2"])
        calls = calls_for({"a1": T1, "a2": T1})
        calls["a1"]["per"]["C"]["multiplicity"] = True
        res = X.call_expansions([comp], calls)
        self.assertEqual((res["expansions"], len(res["dropped_classes"])), ([], 1))
        calls = calls_for({"a1": T1, "a2": T1})
        calls["a1"]["per"]["C"]["copy_rows"] = 1                                                # a count alone (no flag set by hand) is enough as well
        self.assertEqual(len(X.call_expansions([comp], calls)["dropped_classes"]), 1)

    def test_records_without_a_triple_are_counted_by_reason(self):
        comp = recs_comp("PC_a", ["a1", "a2", "a3"])
        calls = calls_for({"a1": T1, "a2": T1, "a3": None})
        calls["a3"]["reason"] = "outgroup_copy_row"
        res = X.call_expansions([comp], calls)
        self.assertEqual(res["reasons"], {"outgroup_copy_row": 1})

    def test_only_such_records_take_part(self):
        comp = recs_comp("PC_a", ["a1", "a2", "a3"])
        comp.recs["a3"]["exonic_len"] = 150
        comp.recs["a3"]["exons"] = [(0, 150)]
        res = X.call_expansions([comp], calls_for({"a1": T1, "a2": T1, "a3": T1}))
        self.assertEqual([m["id"] for m in res["expansions"][0]["members"]], ["a1", "a2"])


class HumanOnlyTests(HermeticCase):
    """Amendment 3, Y4 sensitivity: expansions called relative to human only (the same lifts, a triple reduced to its human locus h)"""

    def setUp(self):
        self.idx = make_indexes()
        self.comp = recs_comp("PC_a", ["a1", "a2", "a3", "a4"])

    def calls(self):
        place = {"a1": (("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000)),
                 "a2": (("h1", 1000, 3000), ("c1", 1000, 3000), ("o2", 1000, 3000)),          # same h and c, another o
                 "a3": (("h1", 1000, 3000), None, ("o1", 1000, 3000)),                        # h ok, no placement in C
                 "a4": (("h2", 1000, 3000), ("c2", 1000, 3000), ("o2", 1000, 3000))}
        none = {"H": set(), "C": set(), "O": set()}
        out = {}
        for rid, (h, c, o) in place.items():
            lifted = {k: {rid: v} for k, v in lifted_for(rid, h, c, o).items()}
            out[rid] = X.classify_record(rid, lifted, none, self.idx)
        return out

    def test_the_unanimity_table_needs_the_same_triple_the_human_only_table_the_same_human_locus(self):
        calls = self.calls()
        self.assertEqual(X.call_expansions([self.comp], calls)["expansions"], [])
        res = X.call_expansions([self.comp], calls, use=("H",), prefix="EH")
        self.assertEqual([(e["id"], e["triple"], [m["id"] for m in e["members"]]) for e in res["expansions"]], [("EH0001", ("gene-H1",), ["a1", "a2", "a3"])])

    def test_a_record_needs_an_ok_human_locus_for_the_human_only_table(self):
        calls = self.calls()
        calls["a3"]["per"]["H"].update(status="outgroup_copy_row", locus=None)
        res = X.call_expansions([self.comp], calls, use=("H",), prefix="EH")
        self.assertEqual([m["id"] for m in res["expansions"][0]["members"]], ["a1", "a2"])

    def test_the_multiplicity_guard_of_the_human_only_table_looks_at_human_only(self):
        calls = self.calls()
        calls["a1"]["per"]["O"]["n_place"] = 2
        self.assertEqual(len(X.call_expansions([self.comp], calls, use=("H",), prefix="EH")["expansions"]), 1)
        calls["a1"]["per"]["H"]["n_place"] = 2
        res = X.call_expansions([self.comp], calls, use=("H",), prefix="EH")
        self.assertEqual((res["expansions"], len(res["dropped_classes"])), ([], 1))

    def test_a_reverse_copy_row_in_another_outgroup_does_not_drop_the_human_only_class(self):
        calls = self.calls()
        calls["a1"]["per"]["C"]["multiplicity"] = True
        self.assertEqual(len(X.call_expansions([self.comp], calls, use=("H",), prefix="EH")["expansions"]), 1)
        calls["a1"]["per"]["H"]["multiplicity"] = True
        res = X.call_expansions([self.comp], calls, use=("H",), prefix="EH")
        self.assertEqual((res["expansions"], len(res["dropped_classes"])), ([], 1))

    def test_the_tables_are_written_with_blank_c_and_o(self):
        calls = self.calls()
        pool = X.Pool()
        pool.components["PC_a"] = self.comp
        res = X.call_expansions([self.comp], calls)
        tmp = tempfile.mkdtemp()
        X.write_outputs(res, pool, [self.comp], calls, tmp)
        with open(os.path.join(tmp, "etable_human_only.tsv")) as fh:
            t = list(csv.reader(fh, delimiter="\t"))
        self.assertEqual(t[0], ["expansion_id", "pool_component_id", "record_id", "biotype", "exon_union_bp", "h", "c", "o", "kind"])
        self.assertEqual([(r[0], r[2], r[5], r[6], r[7]) for r in t[1:]],
                         [("EH0001", "c1:0-2000", "gene-H1", "", ""), ("EH0001", "c1:10000-12000", "gene-H1", "", ""), ("EH0001", "c1:20000-22000", "gene-H1", "", "")])
        with open(os.path.join(tmp, "record_calls.tsv")) as fh:
            rc = list(csv.reader(fh, delimiter="\t"))
        self.assertEqual(rc[0][5:10], ["H_status", "H_detail", "H_n_place", "H_copy_rows", "H_locus"])
        self.assertEqual([r[rc[0].index("H_locus")] for r in rc[1:]], ["gene-H1", "gene-H1", "gene-H1", "gene-H2"])
        self.assertEqual([r[rc[0].index("C_locus")] for r in rc[1:]], ["gene-C1", "gene-C1", "", "gene-C2"])


class ControlPoolTests(HermeticCase):
    def test_a_component_with_an_outgroup_copy_row_record_and_no_expansion_is_a_control_component(self):
        comp = recs_comp("PC_a", ["a1", "a2"])
        calls = calls_for({"a1": T1, "a2": None})
        calls["a2"]["reason"] = "outgroup_copy_row"
        res = X.call_expansions([comp], calls)
        rows = X.control_pool_rows([comp], calls, res, joining={})
        # only the three registered reasons are listed: the record with the one triple of its component has none of them
        self.assertEqual([(r["pool_component_id"], r["record_id"], r["reason"], r["eligible"]) for r in rows],
                         [("PC_a", "c1:10000-12000", "outgroup_copy_row", 1)])

    def test_different_triples_make_a_control_component(self):
        comp = recs_comp("PC_a", ["a1", "a2"])
        calls = calls_for({"a1": T1, "a2": T2})
        res = X.call_expansions([comp], calls)
        rows = X.control_pool_rows([comp], calls, res, joining={})
        self.assertEqual({r["reason"] for r in rows}, {"different_triple"})
        self.assertEqual({r["eligible"] for r in rows}, {1})

    def test_a_record_with_reason_no_placement_is_listed_but_not_eligible_as_a_member_of_a_control_set(self):
        comp = recs_comp("PC_a", ["a1", "a2", "a3"])
        calls = calls_for({"a1": T1, "a2": T2, "a3": None})            # a3: no placement
        res = X.call_expansions([comp], calls)
        rows = X.control_pool_rows([comp], calls, res, joining={})
        self.assertEqual([(r["reason"], r["eligible"]) for r in rows], [("different_triple", 1), ("different_triple", 1), ("no_placement", 0)])

    def test_a_record_placed_on_no_outgroup_locus_is_not_one_of_the_registered_reasons_and_is_not_listed(self):
        comp = recs_comp("PC_a", ["a1", "a2", "a3"])
        calls = calls_for({"a1": T1, "a2": T2, "a3": None})
        calls["a3"]["reason"] = "no_outgroup_locus"
        res = X.call_expansions([comp], calls)
        rows = X.control_pool_rows([comp], calls, res, joining={})
        self.assertEqual([r["record_id"] for r in rows], ["c1:0-2000", "c1:10000-12000"])

    def test_a_component_that_contributes_an_expansion_is_not_in_the_control_pool(self):
        comp = recs_comp("PC_a", ["a1", "a2", "a3"])
        calls = calls_for({"a1": T1, "a2": T1, "a3": T2})
        res = X.call_expansions([comp], calls)
        self.assertEqual(X.control_pool_rows([comp], calls, res, joining={}), [])

    def test_a_component_whose_records_all_lack_a_placement_is_not_a_control_component(self):
        comp = recs_comp("PC_a", ["a1", "a2"])
        calls = calls_for({"a1": None, "a2": None})
        res = X.call_expansions([comp], calls)
        self.assertEqual(X.control_pool_rows([comp], calls, res, joining={}), [])

    def test_a_component_with_a_dropped_class_may_still_be_a_control_component(self):
        comp = recs_comp("PC_a", ["a1", "a2", "a3"])
        calls = calls_for({"a1": T1, "a2": T1, "a3": None})
        calls["a3"]["reason"] = "outgroup_copy_row"
        calls["a2"]["per"]["O"]["n_place"] = 2
        res = X.call_expansions([comp], calls)
        self.assertEqual(res["expansions"], [])
        rows = X.control_pool_rows([comp], calls, res, joining={})
        # fix round, addendum 2: a2 (flagged, a second place in O) is listed as well, under outgroup_copy_row (before: only a3, whose reason was set by hand)
        self.assertEqual([(r["record_id"], r["reason"]) for r in rows], [("c1:10000-12000", "outgroup_copy_row"), ("c1:20000-22000", "outgroup_copy_row")])

    # --- fix round, addendum 2: the records flagged multiplicity are listed under outgroup_copy_row, so that the control pool does not shrink when a reverse copy
    # row no longer turns the record into one without a triple (D2 revised). Before D2 revised such a record WAS outgroup_copy_row, eligible 1.
    def test_a_flagged_record_with_a_triple_is_listed_as_outgroup_copy_row_and_makes_its_component_a_control_component(self):
        comp = recs_comp("PC_a", ["a1", "a2"])
        calls = calls_for({"a1": T1, "a2": T1})
        calls["a2"]["per"]["H"]["copy_rows"] = 1                              # the reverse lift gave a2 a second place in H
        res = X.call_expansions([comp], calls)
        self.assertEqual((res["expansions"], len(res["dropped_classes"])), ([], 1))
        rows = X.control_pool_rows([comp], calls, res, joining={})
        # a1 has the one triple of its component: none of the three registered reasons, not listed
        self.assertEqual([(r["pool_component_id"], r["record_id"], r["reason"], r["eligible"]) for r in rows], [("PC_a", "c1:10000-12000", "outgroup_copy_row", 1)])

    def test_the_flag_counts_whichever_outgroup_it_is_in_and_whichever_way_the_call_was_built(self):
        for how in (lambda p: p.update(copy_rows=1), lambda p: p.update(n_place=2), lambda p: p.update(multiplicity=True)):
            for o in "HCO":
                comp = recs_comp("PC_a", ["a1", "a2"])
                calls = calls_for({"a1": T1, "a2": T1})
                how(calls["a2"]["per"][o])
                rows = X.control_pool_rows([comp], calls, X.call_expansions([comp], calls), joining={})
                self.assertEqual([(r["record_id"], r["reason"], r["eligible"]) for r in rows], [("c1:10000-12000", "outgroup_copy_row", 1)], o)

    def test_the_pool_holds_what_it_held_when_the_reverse_copy_row_gave_the_record_the_reason_outgroup_copy_row_and_no_triple(self):
        comp = recs_comp("PC_a", ["a1", "a2", "a3"])
        old = calls_for({"a1": T1, "a2": None, "a3": T2})              # the reviewed representation: a2 has no triple, its reason is outgroup_copy_row
        old["a2"]["reason"] = "outgroup_copy_row"
        new = calls_for({"a1": T1, "a2": T1, "a3": T2})                # now: a2 keeps the triple of its primary placement and is flagged
        new["a2"]["per"]["H"]["copy_rows"] = 1

        def listing(calls):
            return [(r["pool_component_id"], r["record_id"], r["reason"], r["eligible"]) for r in X.control_pool_rows([comp], calls, X.call_expansions([comp], calls), joining={})]
        self.assertEqual(listing(new), listing(old))
        self.assertEqual([x[2] for x in listing(new)], ["different_triple", "outgroup_copy_row", "different_triple"])

    def test_a_flagged_record_whose_triple_is_missing_for_another_reason_is_listed_as_outgroup_copy_row_and_eligible(self):
        # the reason precedence of the reviewed version: outgroup_copy_row outranks no_placement (REASON_ORDER); the record is eligible although C has no placement for it
        idx = make_indexes()
        none = {"H": set(), "C": set(), "O": set()}

        def call(rid, h, c, o, extra=None):
            lifted = {k: {rid: v} for k, v in lifted_for(rid, h, c, o).items()}
            for k, rs in (extra or {}).items():
                lifted[k][rid] = lifted[k].get(rid, []) + rs
            return X.classify_record(rid, lifted, none, idx)
        comp = recs_comp("PC_a", ["a1", "a2"])
        calls = {"a1": call("a1", ("h1", 1000, 3000), ("c1", 1000, 3000), ("o1", 1000, 3000)),
                 "a2": call("a2", ("h1", 1000, 3000), None, ("o1", 1000, 3000), extra={"H": [lrow("a2_1", 1, "h9", 0, 800, cov=0.9, sid=0.97)]})}
        self.assertEqual((calls["a2"]["triple"], calls["a2"]["reason"], calls["a2"]["multiplicity"]), (None, "no_placement", ["H"]))
        rows = X.control_pool_rows([comp], calls, X.call_expansions([comp], calls), joining={})
        self.assertEqual([(r["record_id"], r["reason"], r["eligible"]) for r in rows], [("c1:10000-12000", "outgroup_copy_row", 1)])

    def test_a_record_without_a_triple_and_without_a_flag_keeps_its_own_reason(self):
        comp = recs_comp("PC_a", ["a1", "a2", "a3"])
        calls = calls_for({"a1": T1, "a2": T1, "a3": None})
        calls["a2"]["per"]["H"]["copy_rows"] = 1
        res = X.call_expansions([comp], calls)
        rows = X.control_pool_rows([comp], calls, res, joining={})
        self.assertEqual([(r["record_id"], r["reason"], r["eligible"]) for r in rows], [("c1:10000-12000", "outgroup_copy_row", 1), ("c1:20000-22000", "no_placement", 0)])

    def test_a_flagged_record_of_a_component_that_contributes_an_expansion_is_not_listed(self):
        comp = recs_comp("PC_a", ["a1", "a2", "a3", "a4"])
        calls = calls_for({"a1": T1, "a2": T1, "a3": T2, "a4": T2})
        calls["a2"]["per"]["H"]["copy_rows"] = 1                              # the class {a1, a2} is dropped, the class {a3, a4} is an expansion
        res = X.call_expansions([comp], calls)
        self.assertEqual(([e["triple"] for e in res["expansions"]], len(res["dropped_classes"])), ([T2], 1))
        self.assertEqual(X.control_pool_rows([comp], calls, res, joining={}), [])

    def test_two_dropped_classes_list_their_flagged_records_and_the_other_records_by_triple(self):
        comp = recs_comp("PC_a", ["a1", "a2", "a3", "a4"])
        calls = calls_for({"a1": T1, "a2": T1, "a3": T2, "a4": T2})
        calls["a2"]["per"]["H"]["copy_rows"] = 1
        calls["a4"]["per"]["O"]["n_place"] = 2
        res = X.call_expansions([comp], calls)
        self.assertEqual((res["expansions"], len(res["dropped_classes"])), ([], 2))
        rows = X.control_pool_rows([comp], calls, res, joining={})
        self.assertEqual([(r["record_id"], r["reason"], r["eligible"]) for r in rows],
                         [("c1:0-2000", "different_triple", 1), ("c1:10000-12000", "outgroup_copy_row", 1),
                          ("c1:20000-22000", "different_triple", 1), ("c1:30000-32000", "outgroup_copy_row", 1)])

    def test_a_flagged_record_keeps_its_joining_copy_rows(self):
        comp = recs_comp("PC_a", ["a1", "a2"])
        calls = calls_for({"a1": T1, "a2": T1})
        calls["a2"]["per"]["H"]["copy_rows"] = 1
        res = X.call_expansions([comp], calls)
        rows = X.control_pool_rows([comp], calls, res, joining={"a2": [("a1", "dropped_M1", "0.97")]})
        self.assertEqual([r["joining_copy_rows"] for r in rows], ["a1:dropped_M1:0.97"])

    def test_joining_copy_rows_are_listed_per_record(self):
        comp = recs_comp("PC_a", ["a1", "a2"])
        calls = calls_for({"a1": T1, "a2": T2})
        res = X.call_expansions([comp], calls)
        joining = {"a1": [("a2", "dropped_M1", "0.97")], "a2": [("a1", "dropped_M1", "0.97")]}
        rows = X.control_pool_rows([comp], calls, res, joining=joining)
        self.assertEqual(rows[0]["joining_copy_rows"], "a2:dropped_M1:0.97")


class OutputTests(HermeticCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        rows, recs = mini_selflift()
        self.pool = X.build_pool(rows, recs)
        self.called = X.called_components(self.pool, "nodes")
        calls = {}
        for c in self.called:
            for r in c.such():
                calls[r["id"]] = dict(triple=T1 if r["id"] in ("gene-A", "gene-C") else None, reason=None if r["id"] in ("gene-A", "gene-C") else "no_placement",
                                      per={k: dict(n_place=1, detail="") for k in "HCO"})
        self.calls = calls
        self.res = X.call_expansions(self.called, calls)
        self.paths = X.write_outputs(self.res, self.pool, self.called, self.calls, self.tmp, rows=rows)

    def table(self, name):
        with open(os.path.join(self.tmp, name)) as fh:
            return list(csv.reader(fh, delimiter="\t"))

    def test_the_e_table_has_the_registered_columns_and_record_ids(self):
        t = self.table("etable.tsv")
        self.assertEqual(t[0], ["expansion_id", "pool_component_id", "record_id", "biotype", "exon_union_bp", "h", "c", "o", "kind"])
        self.assertEqual(t[1], ["EX0001", "PC_gene-A", "c1:1000-3000", "protein_coding", "2000", "gene-H1", "gene-C1", "gene-O1", "annotated"])
        self.assertEqual(t[2][2], "c1:20000-22000")
        self.assertEqual(len(t), 3)

    def test_member_exons_lists_every_annotated_pool_record_0_based_half_open(self):
        t = self.table("member_exons.tsv")
        self.assertEqual(t[0], ["record_id", "pool_component_id", "gene_id", "exons"])
        ids = {r[2] for r in t[1:]}
        self.assertEqual(ids, {"gene-A", "gene-B", "gene-C", "gene-D", "gene-E", "gene-F"})
        a = next(r for r in t if r[2] == "gene-A")
        self.assertEqual((a[0], a[3]), ("c1:1000-3000", "1000-3000"))

    def test_extension_lists_the_unannotated_copy_sites(self):
        t = self.table("extension.tsv")
        self.assertEqual(t[0], ["site_id", "pool_component_id", "contig", "start0", "end", "strand", "source_id", "sequence_id", "cls"])
        self.assertEqual(t[1][:3], ["copy@c9:300-1800-", "PC_gene-F", "c9"])
        self.assertEqual(len(t), 2)

    def test_control_pool_and_dropped_tables_exist_with_their_headers(self):
        t = self.table("control_pool.tsv")
        self.assertEqual(t[0], ["pool_component_id", "record_id", "reason", "joining_copy_rows", "eligible"])
        t = self.table("record_calls.tsv")
        self.assertEqual(t[0][:5], ["record_id", "gene_id", "pool_component_id", "reason", "triple"])

    def test_statuses_hub_control_verdict_and_too_few_such_records(self):
        pool = X.build_pool(*mini_selflift())
        a, d, f = pool.components["PC_gene-A"], pool.components["PC_gene-D"], pool.components["PC_gene-F"]
        calls = {r["id"]: dict(triple=T1, reason=None, per={k: dict(n_place=1, detail="") for k in "HCO"}) for r in a.such()}
        res = X.call_expansions([a], calls)
        tmp = tempfile.mkdtemp()
        X.write_outputs(res, pool, [a], calls, tmp, hub_mode="nodes", control=[d])
        with open(os.path.join(tmp, "components.tsv")) as fh:
            by = {r[0]: r for r in list(csv.reader(fh, delimiter="\t"))[1:]}
        self.assertEqual(rd(os.path.join(tmp, "components.tsv")).splitlines()[0].split("\t"),
                         ["pool_component_id", "members_nodes", "annotated_records", "copy_sites", "such_records", "status", "expansions"])
        self.assertEqual([by[c][5] for c in ("PC_gene-A", "PC_gene-D", "PC_gene-F")], ["verdict", "control", "too_few_such_records"])
        self.assertEqual(by["PC_gene-A"][1:], ["3", "3", "0", "3", "verdict", "1"])
        self.assertEqual(by["PC_gene-F"][1:], ["2", "1", "1", "1", "too_few_such_records", "0"])
        pool.components["PC_gene-F"].sites += ["copy@x:%d-%d+" % (i, i + 9) for i in range(60)]                # 61 nodes: a hub under 'nodes'
        X.write_outputs(res, pool, [a], calls, tmp, hub_mode="nodes", control=[d])
        with open(os.path.join(tmp, "components.tsv")) as fh:
            by = {r[0]: r for r in list(csv.reader(fh, delimiter="\t"))[1:]}
        self.assertEqual(by["PC_gene-F"][5], "hub")

    def test_the_summary_counts_what_was_called_and_dropped(self):
        s = json.loads(rd(os.path.join(self.tmp, "summary.json")))
        self.assertEqual((s["expansions"], s["records_in_expansions"], s["called_components"]), (1, 2, 2))
        self.assertEqual(s["records_without_triple_by_reason"], {"no_placement": 3})
        # every outgroup: the 3 records without a triple have no row in the stub calls (status None) so nothing is counted under a detail
        self.assertEqual(s["no_placement_by_detail"], {"H": {}, "C": {}, "O": {}})
        self.assertEqual((s["components_with_a_triple"], s["components_with_different_loci"], s["classes_in_components_with_different_loci"]), (1, 0, 0))
        for name in ("etable.tsv", "member_exons.tsv", "control_pool.tsv", "extension.tsv"):
            self.assertIn(name, s["sha256"])
            self.assertEqual(len(s["sha256"][name]), 64)

    def test_unplaced_records_are_counted_by_the_detail_of_their_non_placement_in_each_outgroup(self):
        res = dict(self.res)
        calls = {k: dict(v, per={o: dict(v["per"][o]) for o in "HCO"}) for k, v in self.calls.items()}
        calls["gene-B"]["per"]["H"].update(status="no_placement", detail="partial_only")
        calls["gene-B"]["per"]["C"].update(status="no_placement", detail="unmapped")
        calls["gene-D"]["per"]["H"].update(status="no_placement", detail="partial_only")
        other = tempfile.mkdtemp()
        X.write_outputs(res, self.pool, self.called, calls, other)
        s = json.loads(rd(os.path.join(other, "summary.json")))
        self.assertEqual(s["no_placement_by_detail"], {"H": {"partial_only": 2}, "C": {"unmapped": 1}, "O": {}})

    def test_a_dropped_class_is_listed_and_counted_in_the_tables_and_the_summary(self):
        comp = recs_comp("PC_a", ["a1", "a2", "a3"])
        pool = X.Pool()
        pool.components["PC_a"] = comp
        calls = calls_for({"a1": T1, "a2": T1, "a3": T1})
        for rid in calls:
            calls[rid]["per"]["H"].update(status="ok", locus="gene-H1")
        calls["a3"]["per"]["H"].update(multiplicity=True, copy_rows=1, detail="reverse_copy_row")
        res = X.call_expansions([comp], calls)
        tmp = tempfile.mkdtemp()
        X.write_outputs(res, pool, [comp], calls, tmp)
        with open(os.path.join(tmp, "dropped_classes.tsv")) as fh:
            t = list(csv.reader(fh, delimiter="\t"))
        self.assertEqual(t[1:], [["PC_a", "gene-H1", "gene-C1", "gene-O1", "a1,a2,a3", "multiplicity"]])
        s = json.loads(rd(os.path.join(tmp, "summary.json")))
        self.assertEqual((s["dropped_classes"], s["records_in_dropped_classes"], s["expansions"]), (1, 3, 0))
        self.assertEqual(s["records_with_reverse_copy_rows"], {"H": 1, "C": 0, "O": 0})
        self.assertEqual((s["human_only_dropped_classes"], s["human_only_expansions"]), (1, 0))                    # the flag is in H: the human-only class is dropped too
        with open(os.path.join(tmp, "record_calls.tsv")) as fh:
            rc = list(csv.reader(fh, delimiter="\t"))
        self.assertEqual(rc[0][-1], "multiplicity")
        self.assertEqual([r[-1] for r in rc[1:]], ["", "", "H"])                                    # the flagged record names its outgroups

    def test_outputs_are_identical_on_a_second_write(self):
        other = tempfile.mkdtemp()
        rows, recs = mini_selflift()
        X.write_outputs(self.res, self.pool, self.called, self.calls, other, rows=rows)
        for name in ("etable.tsv", "member_exons.tsv", "control_pool.tsv", "extension.tsv", "siblings.tsv", "record_calls.tsv"):
            self.assertEqual(rd(os.path.join(self.tmp, name)), rd(os.path.join(other, name)), name)

    def test_the_summary_counts_hubs_among_all_components_and_among_those_with_two_such_records(self):
        pool = X.build_pool(*mini_selflift())
        pool.components["PC_gene-A"].sites += ["copy@x:%d-%d+" % (i, i + 9) for i in range(60)]          # 3 such records, 63 nodes: a hub
        pool.components["PC_gene-F"].sites += ["copy@y:%d-%d+" % (i, i + 9) for i in range(60)]          # 1 such record, 61 nodes: a hub nobody could call
        called = X.called_components(pool, "nodes")
        calls = {r["id"]: dict(triple=T1, reason=None, per={k: dict(n_place=1, detail="") for k in "HCO"}) for c in called for r in c.such()}
        tmp = tempfile.mkdtemp()
        X.write_outputs(X.call_expansions(called, calls), pool, called, calls, tmp, hub_mode="nodes")
        s = json.loads(rd(os.path.join(tmp, "summary.json")))
        self.assertEqual((s["hubs"], s["hubs_with_two_such_records"]), (2, 1))


class ReadLiftedTests(HermeticCase):
    def test_runs_of_one_outgroup_are_grouped_by_gorilla_record(self):
        tmp = tempfile.mkdtemp()
        import _liftoff
        for label, rows in (("r01", [("gene-A", "h1", 1001, 3000, "+", "0.99", "0.97", 0, "False", "False", [(1001, 3000)], "gene")]),
                            ("r02", [("gene-B", "h2", 5001, 7000, "-", "0.9", "0.9", 0, "False", "False", [(5001, 7000)], "gene"),
                                     ("gene-B_1", "h3", 100, 900, "+", "0.8", "0.96", 1, "False", "False", [(100, 900)], "gene")])):
            d = os.path.join(tmp, "runs", label)
            os.makedirs(d)
            with open(os.path.join(d, "lifted.gff3"), "w") as fh:
                fh.write(lifted_gff(rows))
            _liftoff.parse_lifted(__import__("pathlib").Path(d) / "lifted.gff3", __import__("pathlib").Path(d) / "loci.tsv")
            with open(os.path.join(d, "unmapped.txt"), "w") as fh:
                fh.write("gene-Z\n" if label == "r01" else "")
            with open(os.path.join(d, "done.key"), "w") as fh:
                fh.write("{}")
        lifted, unmapped = X.read_lifted_runs(tmp, known_ids={"gene-A", "gene-B", "gene-Z"})
        self.assertEqual(sorted(lifted), ["gene-A", "gene-B"])
        self.assertEqual(sorted(r["feature_id"] for r in lifted["gene-B"]), ["gene-B", "gene-B_1"])
        self.assertEqual(unmapped, {"gene-Z"})

    def test_only_the_named_runs_are_read_and_each_must_be_finished(self):
        tmp = tempfile.mkdtemp()
        import _liftoff
        import pathlib
        for label, rid, contig in (("r01", "gene-A", "h1"), ("r02", "gene-B", "h2")):
            d = os.path.join(tmp, "runs", label)
            os.makedirs(d)
            with open(os.path.join(d, "lifted.gff3"), "w") as fh:
                fh.write(lifted_gff([(rid, contig, 1001, 3000, "+", "0.99", "0.97", 0, "False", "False", [(1001, 3000)], "gene")]))
            _liftoff.parse_lifted(pathlib.Path(d) / "lifted.gff3", pathlib.Path(d) / "loci.tsv")
            with open(os.path.join(d, "unmapped.txt"), "w") as fh:               # (edited in the fix round: a finished run always has its unmapped.txt)
                fh.write("")
            with open(os.path.join(d, "done.key"), "w") as fh:
                fh.write("{}")
        lifted, _ = X.read_lifted_runs(tmp, known_ids={"gene-A", "gene-B"}, labels=["r02"])
        self.assertEqual(sorted(lifted), ["gene-B"])
        with self.assertRaises(RuntimeError):
            X.read_lifted_runs(tmp, known_ids={"gene-A", "gene-B"}, labels=["r02", "r03"])

    def test_an_unfinished_run_is_an_error(self):
        tmp = tempfile.mkdtemp()
        os.makedirs(os.path.join(tmp, "runs", "r01"))
        with self.assertRaises(RuntimeError):
            X.read_lifted_runs(tmp, known_ids=set(), labels=["r01"])

    def test_a_finished_run_without_unmapped_txt_is_an_error_not_an_empty_list(self):
        tmp = tempfile.mkdtemp()
        import _liftoff
        import pathlib
        d = os.path.join(tmp, "runs", "r01")
        os.makedirs(d)
        with open(os.path.join(d, "lifted.gff3"), "w") as fh:
            fh.write(lifted_gff([("gene-A", "h1", 1001, 3000, "+", "0.99", "0.97", 0, "False", "False", [(1001, 3000)], "gene")]))
        _liftoff.parse_lifted(pathlib.Path(d) / "lifted.gff3", pathlib.Path(d) / "loci.tsv")
        with open(os.path.join(d, "done.key"), "w") as fh:
            fh.write("{}")
        with self.assertRaises(RuntimeError) as cm:
            X.read_lifted_runs(tmp, known_ids={"gene-A"}, labels=["r01"])
        self.assertIn("unmapped.txt", str(cm.exception))

    def test_the_source_record_of_a_row_is_the_id_or_the_id_without_its_copy_tag(self):
        self.assertEqual(X.source_of("gene-A", {"gene-A"}), "gene-A")
        self.assertEqual(X.source_of("gene-A_3", {"gene-A"}), "gene-A")
        self.assertEqual(X.source_of("gene-KEG83_p11", {"gene-KEG83_p11"}), "gene-KEG83_p11")
        self.assertEqual(X.source_of("gene-KEG83_p11_1", {"gene-KEG83_p11"}), "gene-KEG83_p11")

    def test_ids_that_look_like_copy_tags_of_other_ids_are_refused(self):
        X.check_no_id_collisions(["gene-A", "gene-B", "gene-KEG83_p11", "gene-A_1x"])
        X.check_no_id_collisions(["gene-A_1"])                                                  # alone it is nobody's copy
        with self.assertRaises(ValueError) as cm:
            X.check_no_id_collisions(["gene-A", "gene-A_1", "gene-B"])
        self.assertIn("gene-A_1", str(cm.exception))


# ----------------------------------------------------------------------------------------------- group D: Gate 7
def lrpap_fixture():
    recs = {r["id"]: r for r in [
        rec("gene-REF", "c1", 1000, 3000), rec("gene-X1", "c2", 5000, 7000), rec("gene-X2", "c3", 8000, 30000), rec("gene-X3", "c4", 2500, 3500),
        rec("gene-N", "c1", 8000, 9000, exons=[(8000, 8100), (8900, 9000)], biotype="lncRNA"), rec("gene-F1", "c5", 0, 2000), rec("gene-F2", "c6", 0, 2000), rec("gene-F3", "c7", 0, 2000)]}
    rows = [row(r["id"], contig=r["contig"], start0=r["start0"], end=r["end"]) for r in recs.values()]
    rows += [row("gene-REF", cls="dropped_M1", contig="c2", start0=5000, end=7000, note="overlaps gene-X1 (c2)"),
             row("gene-REF", cls="dropped_M1", contig="c3", start0=9000, end=11000, note="overlaps gene-X2 (c3)"),
             row("gene-REF", cls="dropped_M1", contig="c4", start0=2000, end=4000, note="overlaps gene-X3 (c4)"),
             row("gene-REF", cls="dropped_M1", contig="c1", start0=8200, end=8800, note="overlaps gene-N (c1)"),
             row("gene-F1", cls="dropped_M1", contig="c6", start0=0, end=1500, note="overlaps gene-F2 (c6)"),
             row("gene-F2", cls="dropped_M1", contig="c7", start0=0, end=1500, note="overlaps gene-F3 (c7)")]
    for r in rows:
        if r["cls"] in X.COPY_CLASSES:
            r["tag"] = "1"
    spec = dict(reference="gene-REF", not_required=["gene-X2"], fragments=["gene-F1", "gene-F2", "gene-F3"], min_share=0.5, expected_triple=T1)      # addendum 2: the loci of the reference in H, C, O
    return X.build_pool(rows, recs), rows, spec


class PositiveControlTests(HermeticCase):
    def setUp(self):
        self.pool, self.rows, self.spec = lrpap_fixture()
        self.ctrl = X.lrpap_control(self.pool, self.rows, self.spec)

    def test_the_control_records_are_those_with_half_their_exon_bases_inside_a_copy_interval(self):
        self.assertEqual(self.ctrl["qualifying"], ["gene-REF", "gene-X1", "gene-X3"])
        self.assertAlmostEqual(self.ctrl["shares"]["gene-X2"], 2000 / 22000)
        self.assertEqual(self.ctrl["shares"]["gene-N"], 0.0)
        self.assertEqual(self.ctrl["shares"]["gene-X1"], 1.0)

    def test_the_not_required_records_are_checked_to_be_below_half(self):
        self.assertEqual(self.ctrl["not_required"], ["gene-X2"])
        self.assertEqual(self.ctrl["other_non_qualifying"], ["gene-N"])      # non-qualifying and not named: listed, not a failure
        self.assertEqual(self.ctrl["spec_problems"], [])
        spec = dict(self.spec, not_required=["gene-X2", "gene-X1"])          # X1 lies wholly inside a copy interval: the spec is wrong about it
        self.assertEqual(X.lrpap_control(self.pool, self.rows, spec)["spec_problems"], ["gene-X1"])

    def test_a_copy_interval_is_the_landing_interval_of_any_row_of_the_reference_in_the_self_lift(self):
        self.assertEqual(sorted(set(r["contig"] for r in self.rows if r["source_id"] == "gene-REF")), ["c1", "c2", "c3", "c4"])
        self.assertIn("gene-REF", self.ctrl["qualifying"])                   # through its own in_place row
        self.assertIn("gene-X3", self.ctrl["qualifying"])                    # through a dropped_M1 landing site
        rows = [r for r in self.rows if not (r["source_id"] == "gene-REF" and r["cls"] == "dropped_M1" and r["contig"] == "c4")]
        self.assertNotIn("gene-X3", X.lrpap_control(X.build_pool(self.rows, lrpap_fixture_records()), rows, self.spec)["qualifying"])

    def test_the_control_components_are_the_lrpap1_component_and_the_fragment_component(self):
        called = X.called_components(self.pool, "nodes")
        verdict, control = X.split_control(called, self.ctrl)
        self.assertEqual(sorted(c.id for c in control), sorted([self.ctrl["component"], self.ctrl["fragments"]["component"]]))
        self.assertEqual([c.id for c in verdict], [c.id for c in called if c.id not in {x.id for x in control}])

    def test_a_positive_control_component_that_is_not_called_is_an_error(self):
        with self.assertRaises(ValueError):
            X.split_control([], self.ctrl)

    def test_the_fragments_are_another_component(self):
        self.assertNotEqual(self.ctrl["component"], self.ctrl["fragments"]["component"])
        self.assertEqual(self.ctrl["fragments"]["records"], ["gene-F1", "gene-F2", "gene-F3"])

    def test_a_fragment_list_spread_over_two_components_is_an_error(self):
        spec = dict(self.spec, fragments=["gene-F1", "gene-X1"])
        with self.assertRaises(ValueError):
            X.lrpap_control(self.pool, self.rows, spec)

    def called(self):
        return X.called_components(self.pool, "nodes")

    def test_one_class_with_every_qualifying_record_passes(self):
        calls = {r["id"]: dict(triple=None, reason="no_placement", per={k: dict(n_place=1, detail="") for k in "HCO"})
                 for c in self.called() for r in c.such()}
        for rid in ("gene-REF", "gene-X1", "gene-X3"):
            calls[rid]["triple"], calls[rid]["reason"] = T1, None
        for rid in ("gene-F1", "gene-F2"):
            calls[rid]["triple"], calls[rid]["reason"] = T1, None
        res = X.call_expansions(self.called(), calls)
        out = X.check_positive_control(self.ctrl, res, calls)
        self.assertTrue(out["pass"], out)
        self.assertEqual(out["missing"], [])
        self.assertTrue(out["fragments_sibling"])
        self.assertEqual(out["not_required"], [dict(record="gene-X2", in_class=False, reason="no_placement")])

    def test_a_qualifying_record_that_is_missing_fails_the_control_and_names_the_reason(self):
        calls = {r["id"]: dict(triple=None, reason="no_placement", per={k: dict(n_place=1, detail="") for k in "HCO"})
                 for c in self.called() for r in c.such()}
        for rid in ("gene-REF", "gene-X1"):
            calls[rid]["triple"], calls[rid]["reason"] = T1, None
        calls["gene-X3"]["reason"] = "outgroup_copy_row"
        res = X.call_expansions(self.called(), calls)
        out = X.check_positive_control(self.ctrl, res, calls)
        self.assertFalse(out["pass"])
        self.assertEqual(out["missing"], [dict(record="gene-X3", reason="outgroup_copy_row")])

    def test_qualifying_records_split_over_two_classes_fail(self):
        calls = {r["id"]: dict(triple=None, reason="no_placement", per={k: dict(n_place=1, detail="") for k in "HCO"})
                 for c in self.called() for r in c.such()}
        for rid, t in (("gene-REF", T1), ("gene-X1", T1), ("gene-X3", T2), ("gene-N", T2)):
            calls[rid]["triple"], calls[rid]["reason"] = t, None
        res = X.call_expansions(self.called(), calls)
        out = X.check_positive_control(self.ctrl, res, calls)
        self.assertFalse(out["pass"])
        self.assertEqual(out["missing"], [dict(record="gene-X3", reason="in another class")])

    def test_the_not_required_records_may_join_the_class_without_failing_it(self):
        calls = {r["id"]: dict(triple=T1, reason=None, per={k: dict(n_place=1, detail="") for k in "HCO"}) for c in self.called() for r in c.such()}
        res = X.call_expansions(self.called(), calls)
        out = X.check_positive_control(self.ctrl, res, calls)
        self.assertTrue(out["pass"])
        self.assertEqual(out["not_required"], [dict(record="gene-X2", in_class=True, reason="")])

    # --- fix round, addendum 2 (D-m2): the class that holds every qualifying record must also sit on the EXPECTED triple (stricter than the registered text,
    # which asks only that the qualifying records form one class). For LRPAP1 the expected triple is the three gene-LRPAP1 anchors of the outgroups.
    def whole_class_calls(self, triple, records=("gene-REF", "gene-X1", "gene-X3")):
        calls = {r["id"]: dict(triple=None, reason="no_placement", per={k: dict(n_place=1, detail="") for k in "HCO"}) for c in self.called() for r in c.such()}
        for rid in records:
            calls[rid]["triple"], calls[rid]["reason"] = triple, None
        return calls

    def judge(self, calls, ctrl=None):
        return X.check_positive_control(ctrl or self.ctrl, X.call_expansions(self.called(), calls), calls)

    def test_the_control_carries_the_expected_triple_of_its_specification(self):
        self.assertEqual(self.ctrl["expected_triple"], T1)

    def test_a_specification_without_an_expected_triple_is_refused(self):
        spec = {k: v for k, v in self.spec.items() if k != "expected_triple"}
        with self.assertRaises(KeyError):
            X.lrpap_control(self.pool, self.rows, spec)

    def test_a_control_without_an_expected_triple_cannot_be_judged(self):
        ctrl = {k: v for k, v in self.ctrl.items() if k != "expected_triple"}
        with self.assertRaises(KeyError):
            self.judge(self.whole_class_calls(T1), ctrl)

    def test_the_class_that_holds_every_qualifying_record_passes_on_the_expected_triple_and_the_output_says_so(self):
        out = self.judge(self.whole_class_calls(T1))
        self.assertTrue(out["pass"], out)
        self.assertIs(out["triple_ok"], True)
        self.assertEqual(out["expected_triple"], list(T1))

    def test_one_class_with_every_qualifying_record_on_another_triple_fails_although_nothing_is_missing(self):
        # the adversarial case of the logic review (scenario P7): the qualifying records form ONE class, but not at the LRPAP1 loci
        for wrong in (T2, ("gene-H1", "gene-C1", "gene-O2"), ("gene-H1", "gene-C2", "gene-O1"), ("gene-H2", "gene-C1", "gene-O1"), ("gene-C1", "gene-H1", "gene-O1")):
            out = self.judge(self.whole_class_calls(wrong))
            self.assertFalse(out["pass"], wrong)
            self.assertEqual(out["missing"], [], wrong)                                         # the class is whole: the triple alone fails it
            self.assertIs(out["triple_ok"], False, wrong)
            self.assertEqual(out["expected_triple"], list(T1), wrong)
            self.assertEqual(sorted(out["class_members"]), ["gene-REF", "gene-X1", "gene-X3"], wrong)

    def test_the_class_triple_decides_when_a_larger_class_holds_the_qualifying_records_and_others(self):
        calls = self.whole_class_calls(T2, records=("gene-REF", "gene-X1", "gene-X3", "gene-X2", "gene-N"))     # the not-required records join the wrong class
        out = self.judge(calls)
        self.assertEqual((out["pass"], out["triple_ok"]), (False, False))
        out = self.judge(self.whole_class_calls(T1, records=("gene-REF", "gene-X1", "gene-X3", "gene-X2", "gene-N")))
        self.assertEqual((out["pass"], out["triple_ok"]), (True, True))

    def test_the_triple_is_not_judged_when_no_class_holds_every_qualifying_record(self):
        out = self.judge(self.whole_class_calls(T2, records=("gene-REF", "gene-X1")))               # gene-X3 is missing: the control fails on that, as before
        self.assertFalse(out["pass"])
        self.assertIsNone(out["triple_ok"])
        self.assertEqual(out["missing"], [dict(record="gene-X3", reason="no_placement")])
        out = self.judge(self.whole_class_calls(T1, records=("gene-REF", "gene-X1")))
        self.assertEqual((out["pass"], out["triple_ok"]), (False, None))


    # --- fix round 2, G1 (decision of the coordinator): the positive control is judged on the multiplicity flags of the QUALIFYING records only. A not-required record,
    # a lncRNA flank or any other record of the control components that carries a reverse copy row is left out of the control's class (it can neither fail the control
    # nor take the class down with it); it stays flagged, and the registered rule of call_expansions (the class is dropped as a whole) is unchanged.
    def flag(self, calls, record, outgroup="H"):
        calls[record]["per"][outgroup]["copy_rows"] = 1
        return calls

    def test_a_flagged_not_required_record_is_left_out_of_the_control_and_the_control_passes(self):
        members = ("gene-REF", "gene-X1", "gene-X3", "gene-X2", "gene-N")                      # the not-required records join the class
        calls = self.flag(self.whole_class_calls(T1, records=members), "gene-X2", "C")
        view, flagged = X.control_view(self.ctrl, self.called(), calls)
        self.assertEqual(flagged, [dict(record="gene-X2", outgroups=["C"], qualifying=False, left_out=True)])
        self.assertEqual((view["gene-X2"]["triple"], view["gene-X2"]["reason"]), (None, "multiplicity_left_out"))
        out = self.judge(view)
        self.assertTrue(out["pass"], out)
        self.assertEqual(sorted(out["class_members"]), ["gene-N", "gene-REF", "gene-X1", "gene-X3"])
        self.assertEqual(out["not_required"], [dict(record="gene-X2", in_class=False, reason="multiplicity_left_out")])

    def test_the_registered_rule_still_drops_the_whole_class_in_the_real_call_expansions(self):
        members = ("gene-REF", "gene-X1", "gene-X3", "gene-X2", "gene-N")
        calls = self.flag(self.whole_class_calls(T1, records=members), "gene-X2", "C")
        real = X.call_expansions(self.called(), calls)
        self.assertEqual(real["expansions"], [])
        self.assertEqual([(d["triple"], sorted(d["members"]), d["reason"]) for d in real["dropped_classes"]], [(T1, sorted(members), "multiplicity")])
        self.assertTrue(X.has_multiplicity(calls["gene-X2"]["per"]["C"]))                      # and the record is still flagged
        self.assertFalse(self.judge(calls)["pass"])                                            # judged on the registered classes, the control would fail: that is what the view prevents

    def test_a_flagged_qualifying_record_still_fails_the_control(self):
        members = ("gene-REF", "gene-X1", "gene-X3", "gene-X2", "gene-N")
        calls = self.flag(self.whole_class_calls(T1, records=members), "gene-X1", "H")
        view, flagged = X.control_view(self.ctrl, self.called(), calls)
        self.assertEqual(flagged, [dict(record="gene-X1", outgroups=["H"], qualifying=True, left_out=False)])
        self.assertEqual(view["gene-X1"]["triple"], T1)                                          # a qualifying record is never left out
        out = self.judge(view)
        self.assertFalse(out["pass"])
        self.assertEqual({m["reason"] for m in out["missing"]}, {"class dropped for multiplicity"})

    def test_both_kinds_flagged_the_control_fails_and_only_the_not_qualifying_one_is_left_out(self):
        members = ("gene-REF", "gene-X1", "gene-X3", "gene-X2", "gene-N")
        calls = self.flag(self.flag(self.whole_class_calls(T1, records=members), "gene-X1", "H"), "gene-N", "O")
        view, flagged = X.control_view(self.ctrl, self.called(), calls)
        self.assertEqual(flagged, [dict(record="gene-N", outgroups=["O"], qualifying=False, left_out=True), dict(record="gene-X1", outgroups=["H"], qualifying=True, left_out=False)])
        self.assertFalse(self.judge(view)["pass"])

    def test_flagged_records_are_listed_in_the_order_of_their_ids(self):
        # gene-REF comes before gene-N in the component (contig order) and after it by id: the list is sorted by id
        calls = self.flag(self.flag(self.whole_class_calls(T1, records=("gene-REF", "gene-X1", "gene-X3", "gene-N")), "gene-REF", "C"), "gene-N", "O")
        _, flagged = X.control_view(self.ctrl, self.called(), calls)
        self.assertEqual([f["record"] for f in flagged], ["gene-N", "gene-REF"])

    def test_a_flagged_record_of_the_fragment_component_is_left_out_too(self):
        calls = self.flag(self.whole_class_calls(T1, records=("gene-REF", "gene-X1", "gene-X3", "gene-F1", "gene-F2")), "gene-F2", "O")
        view, flagged = X.control_view(self.ctrl, self.called(), calls)
        self.assertEqual(flagged, [dict(record="gene-F2", outgroups=["O"], qualifying=False, left_out=True)])
        self.assertTrue(self.judge(view)["pass"])

    def test_a_flag_in_several_outgroups_names_them_all_in_order(self):
        calls = self.whole_class_calls(T1)
        for o in "OH":
            self.flag(calls, "gene-X1", o)
        _, flagged = X.control_view(self.ctrl, self.called(), calls)
        self.assertEqual(flagged, [dict(record="gene-X1", outgroups=["H", "O"], qualifying=True, left_out=False)])

    def test_without_a_flag_the_view_is_the_calls_and_nothing_is_reported(self):
        calls = self.whole_class_calls(T1)
        view, flagged = X.control_view(self.ctrl, self.called(), calls)
        self.assertEqual((view, flagged), (calls, []))

    def test_the_view_does_not_change_the_calls_it_is_made_from(self):
        calls = self.flag(self.whole_class_calls(T1, records=("gene-REF", "gene-X1", "gene-X3", "gene-X2")), "gene-X2")
        before = copy.deepcopy(calls)
        X.control_view(self.ctrl, self.called(), calls)
        self.assertEqual(calls, before)


class NegativeControlTests(HermeticCase):
    def world(self, n_contigs=26, per_contig=12):
        recs, rows = {}, []
        for c in range(n_contigs):
            for i in range(per_contig):
                rid = f"gene-n{c:02d}x{i:02d}"
                recs[rid] = rec(rid, f"ctg{c:02d}", i * 5000, i * 5000 + (150 if i == 0 else 2000))
                rows.append(row(rid, contig=f"ctg{c:02d}", start0=i * 5000, end=i * 5000 + (150 if i == 0 else 2000)))
        # a pool component: one copy row joins ctg00 records 1 and 2 (they are not eligible)
        rows.append(row("gene-n00x01", cls="dropped_M1", contig="ctg00", start0=10000, end=11000, note="overlaps gene-n00x02 (ctg00)"))
        rows[-1]["tag"] = "1"
        return rows, recs

    def test_one_hundred_records_in_fifty_pairs_on_different_contigs_deterministically(self):
        rows, recs = self.world()
        pool = X.build_pool(rows, recs)
        comps = X.plant_negative(rows, recs, pool)
        self.assertEqual(len(comps), 50)
        ids = [r["id"] for c in comps for r in c.such()]
        self.assertEqual((len(ids), len(set(ids))), (100, 100))
        for c in comps:
            a, b = c.such()
            self.assertNotEqual(a["contig"], b["contig"])
        again = X.plant_negative(rows, recs, pool)
        self.assertEqual([[r["id"] for r in c.such()] for c in comps], [[r["id"] for r in c.such()] for c in again])

    def test_only_long_in_place_records_with_no_copy_row_and_outside_the_pool_are_eligible(self):
        rows, recs = self.world()
        pool = X.build_pool(rows, recs)
        comps = X.plant_negative(rows, recs, pool)
        ids = {r["id"] for c in comps for r in c.such()}
        self.assertFalse({"gene-n00x01", "gene-n00x02"} & ids)                     # in the pool
        self.assertFalse({f"gene-n{c:02d}x00" for c in range(26)} & ids)           # 150 bp: short
        self.assertTrue(all(recs[i]["exonic_len"] >= 200 for i in ids))

    def test_a_record_with_a_copy_row_as_source_is_never_planted_even_when_it_is_in_no_pool_component(self):
        rows, recs = self.world()
        rows.append(row("gene-n05x03", cls="dropped_M1", contig="ctg05", start0=15000, end=17000, note="overlaps gene-n05x03 (ctg05)"))   # a copy on its own source
        rows[-1]["tag"] = "1"
        pool = X.build_pool(rows, recs)
        self.assertNotIn("gene-n05x03", pool.node_component)                 # a one-node component is not kept
        self.assertNotIn("gene-n05x03", X.eligible_for_planting(rows, recs, pool))

    def test_the_planting_pool_is_counted_by_the_three_exclusions(self):
        rows, recs = self.world()
        pool = X.build_pool(rows, recs)
        # 26 contigs x 11 long records; gene-n00x01 is the source of a copy row, gene-n00x02 its target
        self.assertEqual(X.planting_counts(rows, pool), dict(in_place_ge200=286, no_copy_row_as_source=285, target_of_a_copy_row=1, eligible=284))

    def test_the_sample_is_drawn_from_the_sorted_ids_with_the_t2_seed(self):
        rows, recs = self.world()
        pool = X.build_pool(rows, recs)
        eligible = X.eligible_for_planting(rows, recs, pool)
        self.assertEqual(eligible, sorted(eligible))
        want = random.Random(20261007).sample(eligible, 100)
        got = X.plant_negative(rows, recs, pool)
        self.assertEqual(sorted(r["id"] for c in got for r in c.such()), sorted(want))

    def test_too_few_eligible_records_is_an_error(self):
        rows, recs = self.world(n_contigs=2, per_contig=3)
        with self.assertRaises(ValueError):
            X.plant_negative(rows, recs, X.build_pool(rows, recs))

    def test_the_planted_components_get_two_rank_runs(self):
        rows, recs = self.world()
        comps = X.plant_negative(rows, recs, X.build_pool(rows, recs))
        plan = X.rank_plan(comps, label_prefix="ctl")
        self.assertEqual([r["label"] for r in plan], ["ctl1", "ctl2"])
        self.assertEqual([len(r["records"]) for r in plan], [50, 50])

    def test_at_most_one_called_expansion_of_fifty_passes(self):
        rows, recs = self.world()
        comps = X.plant_negative(rows, recs, X.build_pool(rows, recs))
        calls = {r["id"]: dict(triple=None, reason="no_outgroup_locus", per={k: dict(n_place=1, detail="") for k in "HCO"}) for c in comps for r in c.such()}
        for c in comps[:1]:
            for r in c.such():
                calls[r["id"]]["triple"], calls[r["id"]]["reason"] = T1, None
        res = X.call_expansions(comps, calls)
        self.assertEqual(X.check_negative(res, 50), dict(called=1, of=50, bar=1, **{"pass": True}))
        for c in comps[1:2]:
            for r in c.such():
                calls[r["id"]]["triple"], calls[r["id"]]["reason"] = T2, None
        res = X.call_expansions(comps, calls)
        self.assertFalse(X.check_negative(res, 50)["pass"])


# ----------------------------------------------------------------------------------------------- group B2: the budgeted driver
def stub_liftoff_counting(directory):
    """a stub that sleeps (past the hard limit) when its in.gff holds more than STUB_MAX_GENES records (default 99), else copies STUB_LIFTED"""
    return stub_liftoff(directory, max_genes=99, name="liftoff_count_stub.py")


class DriverTests(HermeticCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        gff = os.path.join(self.tmp, "ref.gff")
        with open(gff, "w") as fh:
            fh.write(REF_GFF)
        self.e_root = os.path.join(self.tmp, "e")
        LC.write_run_gffs(gff, {"gene-A": "r01", "gene-B": "r01", "gene-P": "r02", "gene-Z": "ctl1"}, os.path.join(self.e_root, "gff"))
        mk = lambda i, c, s0, e: dict(id=i, component="PC", contig=c, start0=s0, end=e, strand="+", idx=0, rtype="gene", biotype="x", exonic_len=e - s0)  # noqa: E731
        a, b, p_, z = mk("gene-A", "NC_1", 100, 1100), mk("gene-B", "NC_1", 2000, 3000), mk("gene-P", "NC_1", 5000, 6000), mk("gene-Z", "NC_1", 8000, 9000)
        self.plan = {"r01": dict(run=1, label="r01", records=[a, b], shards=[("r01", ["gene-A", "gene-B"])]),
                     "r02": dict(run=2, label="r02", records=[p_], shards=[("r02", ["gene-P"])]),
                     "ctl1": dict(run=1, label="ctl1", records=[z], shards=[("ctl1", ["gene-Z"])])}
        self.lifted = os.path.join(self.tmp, "stub_lifted.gff3")
        with open(self.lifted, "w") as fh:
            fh.write(lifted_gff([("gene-A", "h1", 1001, 2000, "+", "0.98", "0.97", 0, "False", "False", [(1001, 2000)], "gene")]))
        self.kw = dict(e_root=self.e_root, outgroups={"H": dict(species="human", genome="/t/h.fa", mmi="/t/h.mmi")}, reference="/r/g.fa",
                       liftoff_bin=stub_liftoff_counting(self.tmp), minimap2="/m/mm2", hard_s=60, env={"STUB_LIFTED": self.lifted}, ensure_links=False)

    def test_units_run_in_order_and_leave_done_markers(self):
        res = LC.run_budgeted([("H", "r01"), ("H", "r02"), ("H", "ctl1")], self.plan, budget_s=500, **self.kw)
        self.assertEqual((res["exit_code"], res["pending"], res["split"]), (0, [], None))
        self.assertEqual(res["done"], [("H", "r01"), ("H", "r02"), ("H", "ctl1")])
        for lab in ("r01", "r02", "ctl1"):
            self.assertTrue(os.path.exists(os.path.join(self.e_root, "lift", "human", "runs", lab, "done.key")))

    def test_a_unit_that_does_not_fit_the_call_is_left_for_the_next_call(self):
        def ticking():                          # a clock that advances 15 s per reading, so the stub's instant runs still consume the budget
            it = iter(range(0, 10000, 15))
            return lambda: next(it)
        res = LC.run_budgeted([("H", "r01"), ("H", "r02")], self.plan, budget_s=40, startup_s=30, per_record_s=0.0, clock=ticking(), **self.kw)
        self.assertEqual(res["done"], [("H", "r01")])
        self.assertEqual((res["exit_code"], res["pending"]), (75, [("H", "r02")]))
        res = LC.run_budgeted([("H", "r01"), ("H", "r02")], self.plan, budget_s=40, startup_s=30, per_record_s=0.0, clock=ticking(), **self.kw)
        self.assertEqual(res["done"], [("H", "r01"), ("H", "r02")])
        self.assertEqual(res["exit_code"], 0)

    def test_finished_units_are_not_run_again(self):
        LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, **self.kw)
        mark = os.path.join(self.e_root, "lift", "human", "runs", "r01", "time.txt")
        before = os.path.getmtime(mark)
        res = LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, **self.kw)
        self.assertEqual(res["done"], [("H", "r01")])
        self.assertEqual(os.path.getmtime(mark), before)

    def test_a_run_past_the_hard_limit_is_split_and_the_halves_run_in_the_next_call(self):
        kw = dict(self.kw, hard_s=1, env={"STUB_LIFTED": self.lifted, "STUB_MAX_GENES": "1"})
        res = LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, **kw)
        self.assertEqual((res["exit_code"], res["split"]), (75, ("H", "r01")))
        res = LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, **kw)
        self.assertEqual((res["exit_code"], res["done"]), (0, [("H", "r01")]))
        for half in ("r01.s1", "r01.s2"):
            self.assertTrue(os.path.exists(os.path.join(self.e_root, "lift", "human", "runs", half, "done.key")), half)
        with open(os.path.join(self.e_root, "lift", "human", "runs", "r01.s1", "in.gff")) as fh:
            self.assertEqual(sum(1 for l in fh if "\tgene\t" in l), 1)

    def test_a_later_unit_killed_by_the_end_of_the_call_is_run_again_not_split(self):
        kw = dict(self.kw, hard_s=60, env={"STUB_LIFTED": self.lifted, "STUB_MAX_GENES": "1"})        # r01 holds two genes: the stub sleeps 30 s
        res = LC.run_budgeted([("H", "r02"), ("H", "r01")], self.plan, budget_s=500, call_limit_s=2, **kw)
        self.assertEqual((res["exit_code"], res["split"], res["deadline"]), (75, None, ("H", "r01")))
        self.assertEqual(res["done"], [("H", "r02")])
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "lift", "human", "splits.json")))
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "lift", "human", "runs", "r01", "toolong")))

    def test_the_first_unit_of_a_call_has_the_whole_hard_limit_so_a_shard_that_long_is_split(self):
        kw = dict(self.kw, hard_s=1, env={"STUB_LIFTED": self.lifted, "STUB_MAX_GENES": "1"})
        res = LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, call_limit_s=1, **kw)
        self.assertEqual((res["exit_code"], res["split"], res["deadline"]), (75, ("H", "r01"), None))

    def test_a_shard_that_cannot_be_split_is_a_stop(self):
        self.plan["r01"]["records"][1]["start0"] = 500          # gene-B now overlaps gene-A: no record-free cut
        kw = dict(self.kw, hard_s=1, env={"STUB_LIFTED": self.lifted, "STUB_MAX_GENES": "1"})
        with self.assertRaises(LC.StopRun):
            LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, **kw)

    def test_a_failing_liftoff_stops_the_call_and_keeps_what_finished(self):
        kw = dict(self.kw, liftoff_bin=stub_liftoff(self.tmp, exit_code=3))
        with self.assertRaises(LC.StopRun):
            LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, **kw)

    # --- Gate 7 sequencing: the controls are lifted and checked first, in calls of their own
    def gated(self):
        gff = os.path.join(self.tmp, "ref.gff")
        LC.write_run_gffs(gff, {"gene-B": "pos1"}, os.path.join(self.e_root, "gff"))
        self.plan["pos1"] = dict(run=1, label="pos1", records=[dict(self.plan["r01"]["records"][1])], shards=[("pos1", ["gene-B"])])
        return [("H", "pos1"), ("H", "ctl1"), ("H", "r01"), ("H", "r02")], {"pos1": "pos", "ctl1": "neg", "r01": "verdict", "r02": "verdict"}

    def test_a_call_ends_where_the_group_changes_and_the_gate_is_asked_before_the_first_unit_of_each_group(self):
        units, groups = self.gated()
        asked, calls = [], []
        for _ in range(6):
            res = LC.run_budgeted(units, self.plan, budget_s=500, group_of=groups, gate=asked.append, **self.kw)
            calls.append((res["group"], [u[1] for u in res["done"]], res["exit_code"]))
            if res["exit_code"] == 0:
                break
        self.assertEqual(calls, [("pos", ["pos1"], 75), ("neg", ["pos1", "ctl1"], 75), ("verdict", ["pos1", "ctl1", "r01", "r02"], 0)])
        self.assertEqual(asked, ["pos", "neg", "verdict"])

    def test_a_refusal_of_the_gate_stops_the_call_before_any_unit_of_the_group_starts(self):
        units, groups = self.gated()

        def gate(group):
            if group == "verdict":
                raise LC.StopRun("Gate 7: the positive control failed")
        LC.run_budgeted(units, self.plan, budget_s=500, group_of=groups, gate=gate, **self.kw)       # pos
        LC.run_budgeted(units, self.plan, budget_s=500, group_of=groups, gate=gate, **self.kw)       # neg
        with self.assertRaises(LC.StopRun):
            LC.run_budgeted(units, self.plan, budget_s=500, group_of=groups, gate=gate, **self.kw)
        for lab in ("r01", "r02"):
            self.assertFalse(os.path.exists(os.path.join(self.e_root, "lift", "human", "runs", lab)), lab)

    def test_the_gate_is_not_asked_when_every_unit_of_the_call_is_already_done(self):
        units, groups = self.gated()
        for _ in range(4):
            LC.run_budgeted(units, self.plan, budget_s=500, group_of=groups, gate=lambda g: None, **self.kw)
        asked = []
        res = LC.run_budgeted(units, self.plan, budget_s=500, group_of=groups, gate=asked.append, **self.kw)
        self.assertEqual((res["exit_code"], asked), (0, []))

    def test_the_run_command_refuses_to_start_without_an_amendment_commit(self):
        self.assertEqual(LC.main(["run"]), 2)
        self.assertEqual(LC.main(["run", "--amendment-commit", "nothex"]), 2)

    # --- fix round, stop rule 6: the summed Liftoff wall of one outgroup may not pass twice its prediction
    def fake_runs(self, elapsed):
        """run_shard replaced by a function that finishes the run at once and reports the given elapsed seconds, one per run (no sleeping: the budget logic sees
        exactly these numbers)"""
        it = iter(elapsed)

        def fake(d, ids, **kw):
            e = next(it)
            kw["timing"]["elapsed_s"] = e
            (d / "done.key").write_text("{}")
            (d / "loci.tsv").write_text("x")
            return dict(status="done", wall_s=e, elapsed_s=e, peak_rss_kb=1, peak_summed_kb=1)
        return mock.patch.object(LC, "run_shard", fake)

    def test_the_summed_wall_of_an_outgroup_stops_the_call_when_it_passes_twice_the_prediction(self):
        units = [("H", "r01"), ("H", "r02"), ("H", "ctl1")]
        with self.fake_runs([1.1, 1.1, 1.1]), self.assertRaises(LC.StopRun) as cm:
            LC.run_budgeted(units, self.plan, budget_s=500, wall_pred_s={"H": 1.0}, **self.kw)                    # limit 2 x 1.0 s
        msg = str(cm.exception)
        self.assertIn("time budget", msg)
        self.assertIn("human", msg)
        self.assertIn("took 2 s in all, more than 2 x the predicted 1 s (2 s)", msg)
        runs = os.path.join(self.e_root, "lift", "human", "runs")
        for lab in ("r01", "r02"):
            self.assertTrue(os.path.exists(os.path.join(runs, lab, "done.key")), lab)                       # what finished stays valid
        self.assertFalse(os.path.exists(os.path.join(runs, "ctl1")), "no run starts after the limit was passed")
        led = LC.WallLedger(self.e_root)
        self.assertAlmostEqual(led.total("H"), 2.2)
        self.assertEqual([(e["label"], e["status"]) for e in led.entries("H")], [("r01", "done"), ("r02", "done")])

    def test_the_wall_just_at_the_limit_does_not_stop_and_a_hair_over_does(self):
        units = [("H", "r01"), ("H", "r02")]
        with self.fake_runs([1.0, 1.0]):
            self.assertEqual(LC.run_budgeted(units, self.plan, budget_s=500, wall_pred_s={"H": 1.0}, **self.kw)["exit_code"], 0)          # 2.0 s = the limit: not over
        self.assertEqual(LC.WallLedger(self.e_root).total("H"), 2.0)
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, True)
        kw = dict(self.kw, e_root=os.path.join(tmp, "e"))
        LC.write_run_gffs(os.path.join(self.tmp, "ref.gff"), {"gene-A": "r01", "gene-B": "r01", "gene-P": "r02"}, os.path.join(tmp, "e", "gff"))
        with self.fake_runs([1.0, 1.001]), self.assertRaises(LC.StopRun):
            LC.run_budgeted(units, self.plan, budget_s=500, wall_pred_s={"H": 1.0}, **kw)

    def test_the_stop_is_not_forgotten_by_the_next_call(self):
        units = [("H", "r01"), ("H", "r02"), ("H", "ctl1")]
        with self.fake_runs([1.1, 1.1]), self.assertRaises(LC.StopRun):
            LC.run_budgeted(units, self.plan, budget_s=500, wall_pred_s={"H": 1.0}, **self.kw)
        launched = []

        def must_not_run(*a, **k):
            launched.append(a)
            raise AssertionError("a run was started after the budget was passed")
        with mock.patch.object(LC, "run_shard", must_not_run), self.assertRaises(LC.StopRun) as cm:           # a new call: the sum is read from E_ROOT/wall.json
            LC.run_budgeted(units, self.plan, budget_s=500, wall_pred_s={"H": 1.0}, **self.kw)
        self.assertIn("time budget", str(cm.exception))
        self.assertEqual(launched, [])
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "lift", "human", "runs", "ctl1")))
        self.assertEqual(len(LC.WallLedger(self.e_root).entries("H")), 2)

    def test_the_wall_budget_is_per_outgroup(self):
        og = dict(self.kw["outgroups"], C=dict(species="chimpanzee", genome="/t/c.fa", mmi="/t/c.mmi"))
        kw = dict(self.kw, outgroups=og)
        pred = {"H": 100.0, "C": 1.0}
        with self.fake_runs([1.1, 1.1]), self.assertRaises(LC.StopRun) as cm:
            LC.run_budgeted([("C", "r01"), ("C", "r02")], self.plan, budget_s=500, wall_pred_s=pred, **kw)
        self.assertIn("chimpanzee", str(cm.exception))
        with self.fake_runs([5.0]):
            res = LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, wall_pred_s=pred, **kw)               # human is within its own budget
        self.assertEqual(res["exit_code"], 0)
        self.assertEqual(res["wall_s"], {"C": 2, "H": 5})

    def test_a_run_that_raises_a_stop_still_counts_its_wall(self):
        def failing(d, ids, **kw):
            kw["timing"]["elapsed_s"] = 3.5
            raise LC.StopRun("liftoff failed")
        with mock.patch.object(LC, "run_shard", failing), self.assertRaises(LC.StopRun):
            LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, wall_pred_s={"H": 100.0}, **self.kw)
        self.assertEqual([(e["label"], e["status"], e["elapsed_s"]) for e in LC.WallLedger(self.e_root).entries("H")], [("r01", "stopped", 3.5)])

    def test_a_run_killed_at_the_hard_limit_counts_its_wall_too(self):
        kw = dict(self.kw, hard_s=1, env={"STUB_LIFTED": self.lifted, "STUB_MAX_GENES": "1"})                # r01 holds two genes: the stub sleeps, killed at 1 s
        LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, wall_pred_s={"H": 100.0}, **kw)
        e = LC.WallLedger(self.e_root).entries("H")[0]
        self.assertEqual(e["status"], "toolong")
        self.assertGreaterEqual(e["elapsed_s"], 1.0)

    def test_the_predicted_wall_of_an_outgroup_is_the_sum_over_its_units(self):
        units = [("H", "r01"), ("H", "r02"), ("C", "r01")]
        got = LC.predicted_wall_by_outgroup(units, self.plan, per_record_s=0.1, startup_s=10.0)
        self.assertAlmostEqual(got["H"], (10 + 2 * 0.1 * 1.25) + (10 + 1 * 0.1 * 1.25))
        self.assertAlmostEqual(got["C"], 10 + 2 * 0.1 * 1.25)

    # --- fix round, major 5: a finished run is finished only if it was made by THIS command from THIS in.gff
    def test_a_finished_run_made_with_another_command_is_run_again(self):
        LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, **self.kw)                                  # made by one binary
        other = stub_liftoff(self.tmp, name="other_liftoff.py")
        res = LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, **dict(self.kw, liftoff_bin=other))
        self.assertEqual(res["exit_code"], 0)
        key = json.loads(rd(os.path.join(self.e_root, "lift", "human", "runs", "r01", "done.key")))
        self.assertEqual(key["cmd"][0], other, "the output of the first binary was not taken for a finished run")

    def test_a_finished_run_whose_in_gff_changed_is_run_again(self):
        LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, **self.kw)
        src = os.path.join(self.e_root, "gff", "r01", "in.gff")
        with open(src, "a") as fh:
            fh.write("# the reference GFF of this run was regenerated\n")
        LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, **self.kw)
        key = json.loads(rd(os.path.join(self.e_root, "lift", "human", "runs", "r01", "done.key")))
        import hashlib
        self.assertEqual(key["in_gff_sha256"], hashlib.sha256(open(src, "rb").read()).hexdigest()[:16])

    def test_a_done_key_that_is_not_json_is_not_a_finished_run(self):
        LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, **self.kw)
        kpath = os.path.join(self.e_root, "lift", "human", "runs", "r01", "done.key")
        with open(kpath, "w") as fh:
            fh.write("{\"cmd\": [")
        res = LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, **self.kw)
        self.assertEqual(res["exit_code"], 0)
        json.loads(rd(kpath))

    def test_is_done_compares_the_stored_key_with_the_expected_one(self):
        d = os.path.join(self.tmp, "run_x")
        os.makedirs(d)
        self.assertFalse(LC._is_done(__import__("pathlib").Path(d)))
        for fn, text in (("done.key", json.dumps({"a": 1})), ("loci.tsv", "x\n")):
            with open(os.path.join(d, fn), "w") as fh:
                fh.write(text)
        P = __import__("pathlib").Path(d)
        self.assertTrue(LC._is_done(P))                                   # no expectation given: the files decide (the readers of finished runs)
        self.assertTrue(LC._is_done(P, {"a": 1}))
        self.assertFalse(LC._is_done(P, {"a": 2}))

    # --- fix round: the start of the call is the start of the process
    def test_the_start_up_work_of_the_process_is_charged_to_the_call(self):
        clock = lambda: 520.0                                                  # noqa: E731  (520 s after the start of the process, when the first unit is reached)
        res = LC.run_budgeted([("H", "r01"), ("H", "r02")], self.plan, budget_s=530, startup_s=15, per_record_s=0.0, clock=clock, t0=0.0, **self.kw)
        self.assertEqual(res["done"], [("H", "r01")])
        self.assertEqual((res["exit_code"], res["pending"]), (75, [("H", "r02")]))
        res = LC.run_budgeted([("H", "r02"), ("H", "ctl1")], self.plan, budget_s=530, startup_s=15, per_record_s=0.0, clock=clock, **self.kw)    # no t0: from now
        self.assertEqual((res["exit_code"], res["done"]), (0, [("H", "r02"), ("H", "ctl1")]))

    def test_the_defaults_leave_room_for_the_kill_sequence_inside_the_lock_timeout(self):
        self.assertEqual(LC.CALL_LIMIT_S, 555)
        self.assertLess(LC.CALL_LIMIT_S + LC.TERM_GRACE_S + LC.KILL_WAIT_S + 5, 585)          # 585 = RLOCK_TIMEOUT, itself below the 590 s hard limit
        self.assertEqual(LC.POLL_S, 0.5)
        self.assertEqual(LC.WALL_FACTOR, 2.0)

    def test_the_run_line_reports_the_summed_peak_and_the_elapsed_time(self):
        lines = []
        LC.run_budgeted([("H", "r01")], self.plan, budget_s=500, log=lines.append, **self.kw)
        self.assertTrue(any("summed" in l and "elapsed" in l for l in lines), lines)


# ----------------------------------------------------------------------------------------------- group E: the input-only report
def write_selflift(root, rows, recs, support=None):
    os.makedirs(os.path.join(root, "gorilla"), exist_ok=True)
    os.makedirs(os.path.join(root, "support"), exist_ok=True)
    cols = ["species", "shard", "source_id", "name", "rtype", "biotype", "stratum", "short", "src_contig", "src_start0", "src_end", "src_strand", "cls", "tag",
            "contig", "start0", "end", "strand", "coverage", "sequence_id", "exons", "exonic_len", "note"]
    with open(os.path.join(root, "gorilla", "liftoff_loci.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(os.path.join(root, "gorilla", "records.tsv"), "w") as fh:
        fh.write("idx\tid\tname\trtype\tbiotype\tcontig\tstart0\tend\tstrand\texons\texonic_len\tshard\n")
        for i, r in enumerate(recs.values()):
            fh.write("\t".join(map(str, [i, r["id"], r["name"], r["rtype"], r["biotype"], r["contig"], r["start0"], r["end"], r["strand"],
                                         X.iv_str(r["exons"]), r["exonic_len"], "s"])) + "\n")
    for lib, sup in (support or {}).items():
        with open(os.path.join(root, "support", lib + ".tsv"), "w") as fh:
            fh.write("source_id\tcls\tcontig\tstart0\tend\treads_capped\n")
            for rid, n in sup.items():
                fh.write(f"{rid}\tin_place\tc1\t0\t10\t{n}\n")


class ReportTests(HermeticCase):
    def test_the_pool_report_counts_components_hubs_the_called_set_and_support(self):
        rows, recs = mini_selflift()
        root = tempfile.mkdtemp()
        write_selflift(root, rows, recs, support={"gorilla_KB3781": {"gene-A": 5, "gene-B": 2, "gene-C": 0, "gene-D": 9, "gene-E": 1},
                                                  "gorilla_OR6737": {"gene-A": 2, "gene-C": 2, "gene-D": 2, "gene-E": 2}})
        rep = X.pool_report(root, support_libs=("gorilla_KB3781", "gorilla_OR6737"))
        self.assertEqual((rep["components_size_ge2"], rep["annotated_records_in_components"], rep["annotated_ge200bp"]), (3, 6, 6))
        self.assertEqual(rep["components_ge2_such"], 2)
        self.assertEqual(rep["nodes"]["called"], 2)
        self.assertEqual(rep["nodes"]["records"], 5)
        self.assertEqual(rep["nodes"]["plan"]["records_per_run"], [2, 2, 1])
        self.assertEqual(rep["supported"]["gorilla_KB3781"]["components"], 1)          # A,B (5 and 2 reads): yes; D,E: 9 and 1: no
        self.assertEqual(rep["supported"]["gorilla_OR6737"]["components"], 2)
        self.assertEqual(rep["supported"]["either"]["components"], 2)

    def test_the_two_hub_readings_can_differ(self):
        rows, recs = HubAndCalledSetTests().big_pool(n_short=45, n_sites=5)
        root = tempfile.mkdtemp()
        write_selflift(root, rows, recs)
        rep = X.pool_report(root, support_libs=())
        self.assertEqual((rep["nodes"]["hubs"], rep["records"]["hubs"]), (1, 0))
        self.assertEqual((rep["nodes"]["called"], rep["records"]["called"]), (0, 1))

    def test_the_expected_numbers_are_compared_and_differences_listed(self):
        rep = dict(components_size_ge2=609, annotated_records_in_components=2976, annotated_ge200bp=2249, components_ge2_such=307,
                   records=dict(hubs=5, called=302, records=1446, by_rtype={"gene": 944, "pseudogene": 502}, J=44,
                                plan=dict(records_per_run=[302, 302, 167, 118, 83], runs_with_at_most_10=26)),
                   nodes=dict(hubs=7, called=300, records=1403, by_rtype={"gene": 907, "pseudogene": 496}, J=44,
                              plan=dict(records_per_run=[300, 300, 165, 116, 81], runs_with_at_most_10=26)),
                   supported={"gorilla_KB3781": dict(components=80, hubs_nodes=4, hubs_records=4),
                              "gorilla_OR6737": dict(components=139, hubs_nodes=5, hubs_records=5), "either": dict(components=157, hubs_nodes=5, hubs_records=5)})
        # the registered rule (hubs count the copy@ sites) gives 7 / 300 / 1,403; the older counting of annotated records only gives 5 / 302 / 1,446
        self.assertEqual(X.compare_expected(rep, "nodes"), [])
        self.assertEqual(X.compare_expected(rep, "records"), [])
        self.assertEqual(X.EXPECTED["hubs"], 7)
        self.assertEqual((X.EXPECTED["called"], X.EXPECTED["records"], X.EXPECTED["records_per_run_first5"]), (300, 1403, [300, 300, 165, 116, 81]))
        rep["nodes"]["hubs"] = 5
        rep["supported"]["gorilla_OR6737"]["hubs_nodes"] = 6
        diffs = X.compare_expected(rep, "nodes")
        self.assertIn("hubs: expected 7, got 5", diffs)
        self.assertIn("supported gorilla_OR6737: expected 5 hubs, got 6", diffs)


# ----------------------------------------------------------------------------------------------- group F: from finished runs to the E-table (stub run outputs)
def lrpap_fixture_records():
    pool, rows, spec = lrpap_fixture()
    return {rid: {k: v for k, v in r.items() if k != "cls"} for c in pool.components.values() for rid, r in c.recs.items()}


def make_run_dir(base, label, placements, unmapped=()):
    """a finished run: loci.tsv parsed from a stub Liftoff GFF3, unmapped.txt, done.key. placements = {record: (contig, start0, end)} (0-based half-open)"""
    import pathlib
    import _liftoff
    d = os.path.join(base, "runs", label)
    os.makedirs(d, exist_ok=True)
    rows = [(rid, c, s0 + 1, e, "+", "0.99", "0.98", 0, "False", "False", [(s0 + 1, e)], "gene") for rid, (c, s0, e) in placements.items()]
    with open(os.path.join(d, "lifted.gff3"), "w") as fh:
        fh.write(lifted_gff(rows))
    _liftoff.parse_lifted(pathlib.Path(d) / "lifted.gff3", pathlib.Path(d) / "loci.tsv")
    with open(os.path.join(d, "unmapped.txt"), "w") as fh:
        fh.write("".join(f"{u}\n" for u in unmapped))
    with open(os.path.join(d, "done.key"), "w") as fh:
        fh.write(json.dumps({"label": label, "placements": sorted(placements)}))


def write_table(path, rows):
    cols = list(rows[0].keys())
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")


OG_SPECIES = (("H", "human", "h"), ("C", "chimpanzee", "c"), ("O", "orangutan", "o"))


def make_fake_genome(fasta, names=("s1", "s2"), mmi=True):
    """a tiny FASTA with its .fai (newer than the FASTA) and a minimap2 index that has the header of the real ones (k15 w10, one sequence per .fai line): what the
    preflight of the links looks at. Nothing here is a real genome or index."""
    import struct
    os.makedirs(os.path.dirname(fasta), exist_ok=True)
    with open(fasta, "w") as fh:
        for n in names:
            fh.write(f">{n}\nACGTACGTAC\n")
    old = time.time() - 100
    os.utime(fasta, (old, old))
    with open(fasta + ".fai", "w") as fh:
        for k, n in enumerate(names):
            fh.write(f"{n}\t10\t{len(n) + 2 + k * (len(n) + 14)}\t10\t11\n")
    if mmi:
        with open(fasta + ".mmi", "wb") as fh:
            fh.write(b"MMI\x02" + struct.pack("<5I", 10, 15, 14, len(names), 0) + b"\0" * 16)


class WorldCase(HermeticCase):
    """a small gorilla pool and three outgroups. Control components: the LRPAP1-like component (REF X1 X2 X3 N; ranks follow contig then start) and the
    fragment component (F1 F2 F3); verdict components: V (V1 V2) and W (W1 W2); eight single planted-negative candidates on eight contigs.
    Outgroup loci: <k>1 .. <k>4 (single-copy) in every outgroup."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        pool, rows, self.spec = lrpap_fixture()
        recs = lrpap_fixture_records()
        for rid, c, s0 in (("gene-V1", "c8", 1000), ("gene-V2", "c9", 1000), ("gene-W1", "c10", 1000), ("gene-W2", "c11", 1000)):
            recs[rid] = rec(rid, c, s0, s0 + 2000)
            rows.append(row(rid, contig=c, start0=s0, end=s0 + 2000))
        rows.append(row("gene-V1", cls="dropped_M1", contig="c9", start0=1000, end=3000, note="overlaps gene-V2 (c9)"))
        rows.append(row("gene-W1", cls="dropped_M1", contig="c11", start0=1000, end=3000, note="overlaps gene-W2 (c11)"))
        for k in range(8):
            rid = f"gene-G{k}"
            recs[rid] = rec(rid, f"g{k}", 5000, 7000)
            rows.append(row(rid, contig=f"g{k}", start0=5000, end=7000))
        for r in rows:
            if r["cls"] in X.COPY_CLASSES:
                r["tag"] = "1"
        self.rows, self.recs = rows, recs
        self.root = os.path.join(self.tmp, "lift")
        write_selflift(self.root, rows, recs)
        self.og = {}
        for key, sp, loc in OG_SPECIES:
            os.makedirs(os.path.join(self.root, sp), exist_ok=True)
            write_table(os.path.join(self.root, sp, "liftoff_loci.tsv"),
                        outgroup_rows(*[(f"gene-{key}{n}", f"{loc}{n}", 1000, 3000, []) for n in (1, 2, 3, 4)]))
            fasta = os.path.join(self.tmp, "genomes", f"{sp}.fa")
            make_fake_genome(fasta)
            self.og[key] = dict(species=sp, genome=fasta, mmi=fasta + ".mmi", selflift=os.path.join(self.root, sp, "liftoff_loci.tsv"))
        self.e_root = os.path.join(self.tmp, "e")
        self.w = LC.plan_from_real_inputs("nodes", lift_root=self.root, spec=self.spec, n_neg=4)

    def load_indexes(self):
        return {k: X.OutgroupIndex(X.read_selflift(t["selflift"]), k) for k, t in self.og.items()}

    def place(self, label, loci, outgroups="HCO"):
        """finish run `label` in the outgroups: loci = {gorilla record: locus number 1..4 (or None = no row)}"""
        for key, sp, loc in OG_SPECIES:
            if key not in outgroups:
                continue
            pl = {rid: (f"{loc}{n}", 1000, 3000) for rid, n in loci.items() if n is not None}
            make_run_dir(os.path.join(self.e_root, "lift", sp), label, pl, unmapped=[rid for rid, n in loci.items() if n is None])

    def add_copy_row(self, species, label, rid, tag=1, contig="x9", s0=0, e=800, cov=0.97):
        """append a reverse copy row (extra_copy_number > 0) of gorilla record `rid` to the finished run `label` of `species`: a second place in that genome"""
        path = os.path.join(self.e_root, "lift", species, "runs", label, "loci.tsv")
        with open(path, "a") as fh:
            fh.write("\t".join(map(str, [f"{rid}_{tag}", f"{rid}_{tag}", tag, "gene", contig, s0, e, "+", cov, 0.97, 0, 0, f"{s0}-{e}", e - s0])) + "\n")

    def ranks(self, plan_key):
        return [[r["id"] for r in run["records"]] for run in self.w[plan_key]]

    def finish_controls(self, neg_ok=True):
        """positive control: the LRPAP1-like records all on locus 1, the fragments on locus 2; planted pairs on loci 1 and 2 (neg_ok=False: both members of
        every pair on locus 1)"""
        frag = {"gene-F1", "gene-F2", "gene-F3"}
        for run in self.w["pos"]:
            self.place(run["label"], {r["id"]: (2 if r["id"] in frag else 1) for r in run["records"]})
        for j, run in enumerate(self.w["ctl"]):
            self.place(run["label"], {r["id"]: (1 if (j == 0 or not neg_ok) else 2) for r in run["records"]})

    def finish_verdict(self):
        locus = {"gene-V1": 3, "gene-V2": 3, "gene-W1": 4, "gene-W2": 2}
        for run in self.w["plan"]:
            self.place(run["label"], {r["id"]: locus[r["id"]] for r in run["records"]})


class PlansTests(WorldCase):
    def test_the_control_components_leave_the_verdict_plan(self):
        ctl_ids = {r["id"] for run in self.w["pos"] for r in run["records"]}
        ver_ids = {r["id"] for run in self.w["plan"] for r in run["records"]}
        self.assertEqual(ctl_ids, {"gene-REF", "gene-N", "gene-X1", "gene-X2", "gene-X3", "gene-F1", "gene-F2", "gene-F3"})
        self.assertEqual(ver_ids, {"gene-V1", "gene-V2", "gene-W1", "gene-W2"})
        self.assertFalse(ctl_ids & ver_ids)

    def test_the_positive_runs_share_ranks_between_the_two_components(self):
        self.assertEqual([run["label"] for run in self.w["pos"]], ["pos1", "pos2", "pos3", "pos4", "pos5"])
        self.assertEqual([sorted(x) for x in self.ranks("pos")],
                         [["gene-F1", "gene-REF"], ["gene-F2", "gene-N"], ["gene-F3", "gene-X1"], ["gene-X2"], ["gene-X3"]])

    def test_the_planted_runs_hold_one_member_of_every_fake_pair(self):
        self.assertEqual([run["label"] for run in self.w["ctl"]], ["ctl1", "ctl2"])
        self.assertEqual([len(run["records"]) for run in self.w["ctl"]], [2, 2])
        planted = {r["id"] for run in self.w["ctl"] for r in run["records"]}
        self.assertTrue(planted <= {f"gene-G{k}" for k in range(8)})


class EvaluateControlsTests(WorldCase):
    def test_both_controls_pass_before_any_verdict_run_exists(self):
        self.finish_controls()
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        self.assertTrue(res["positive"]["pass"], res["positive"])
        self.assertEqual(res["negative"], dict(called=0, of=2, bar=1, called_pairs=[], records_with_triple=4, pairs_with_both_triples=2, dropped_pairs=0, **{"pass": True}))
        self.assertTrue(res["valid"])
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "lift", "human", "runs", "r01")))

    def test_a_qualifying_record_with_a_reverse_copy_row_drops_the_class_and_fails_the_positive_control_saying_so(self):
        self.finish_controls()
        self.add_copy_row("human", "pos3", "gene-X1")
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        self.assertFalse(res["positive"]["pass"])
        self.assertEqual({m["reason"] for m in res["positive"]["missing"]}, {"class dropped for multiplicity"})
        self.assertIn("class dropped for multiplicity", LC.failure_text(res))
        with self.assertRaises(LC.StopRun) as cm:
            LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("neg")
        self.assertIn("class dropped for multiplicity", str(cm.exception))

    def test_the_not_required_records_may_join_the_class(self):
        self.finish_controls()
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        self.assertEqual(sorted(res["positive"]["class_members"]), ["gene-N", "gene-REF", "gene-X1", "gene-X2", "gene-X3"])
        self.assertTrue(res["positive"]["fragments_sibling"] is False)

    # --- fix round 2, G1: the flags of records that are not qualifying do not decide the positive control
    def test_a_flagged_not_required_record_joins_the_class_and_the_control_passes_with_the_flag_reported(self):
        self.finish_controls()
        self.add_copy_row("human", "pos4", "gene-X2")                              # gene-X2 is not required; it has a second place in the human genome
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        p = res["positive"]
        self.assertTrue(p["pass"], p)
        self.assertEqual(p["flagged_members"], [dict(record="gene-X2", outgroups=["H"], qualifying=False, left_out=True)])
        self.assertEqual(p["left_out"], ["gene-X2"])
        self.assertEqual(sorted(p["class_members"]), ["gene-N", "gene-REF", "gene-X1", "gene-X3"])
        self.assertEqual(p["not_required"], [dict(record="gene-X2", in_class=False, reason="multiplicity_left_out")])
        self.assertTrue(res["valid"])
        self.assertEqual(LC.failure_text(res), "")
        # the registered rule is untouched: the same calls, through the real call_expansions, drop the class as a whole (and the control E-table follows it)
        calls = LC.classify_group(self.w["pos"], self.w["control"], self.e_root, self.og, self.load_indexes())
        real = X.call_expansions(self.w["control"], calls, prefix="EC")
        self.assertEqual([(sorted(d["members"]), d["reason"]) for d in real["dropped_classes"]], [(["gene-N", "gene-REF", "gene-X1", "gene-X2", "gene-X3"], "multiplicity")])
        self.assertEqual([e["triple"] for e in res["_control_expansions"]], [e["triple"] for e in real["expansions"]])
        self.assertTrue(X.has_multiplicity(calls["gene-X2"]["per"]["H"]))

    def test_a_flagged_qualifying_record_fails_the_control_and_the_text_names_the_record_and_the_outgroup(self):
        self.finish_controls()
        self.add_copy_row("human", "pos3", "gene-X1")
        self.add_copy_row("chimpanzee", "pos3", "gene-X1")
        self.add_copy_row("orangutan", "pos4", "gene-X2")                          # a not-required record flagged too: it is left out, not blamed
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        p = res["positive"]
        self.assertFalse(p["pass"])
        self.assertEqual(p["flagged_members"], [dict(record="gene-X1", outgroups=["H", "C"], qualifying=True, left_out=False),
                                                dict(record="gene-X2", outgroups=["O"], qualifying=False, left_out=True)])
        self.assertEqual(p["left_out"], ["gene-X2"])
        text = LC.failure_text(res)
        self.assertIn("positive control (LRPAP1)", text)
        self.assertIn("flagged multiplicity", text)
        self.assertIn("gene-X1: H,C", text)
        self.assertIn("gene-X2: O", text)
        self.assertIn("left out of the control", text)
        self.assertIn("class dropped for multiplicity", text)

    def test_the_gate_lets_the_next_group_start_when_only_a_not_required_record_is_flagged(self):
        self.finish_controls()
        self.add_copy_row("chimpanzee", "pos4", "gene-X2")
        LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("neg")                  # no StopRun
        with open(os.path.join(self.e_root, "gate7.json")) as fh:
            g = json.load(fh)
        self.assertEqual((g["positive"]["pass"], g["positive"]["left_out"]), (True, ["gene-X2"]))
        self.assertEqual(g["positive"]["flagged_members"], [dict(record="gene-X2", outgroups=["C"], qualifying=False, left_out=True)])

    def test_the_gate_stops_when_a_qualifying_record_is_flagged_and_says_which(self):
        self.finish_controls()
        self.add_copy_row("orangutan", "pos1", "gene-REF")
        with self.assertRaises(LC.StopRun) as cm:
            LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("neg")
        self.assertIn("gene-REF: O", str(cm.exception))

    def test_the_failure_text_of_a_hand_built_result_lists_the_flagged_members_and_the_left_out_ones(self):
        p = {"pass": False, "missing": [dict(record="gene-A", reason="class dropped for multiplicity")], "triple_ok": None,
             "flagged_members": [dict(record="gene-A", outgroups=["H", "O"], qualifying=True, left_out=False), dict(record="gene-B", outgroups=["C"], qualifying=False, left_out=True)],
             "left_out": ["gene-B"]}
        text = LC.failure_text(dict(positive=p))
        self.assertIn("flagged multiplicity (record: outgroups): gene-A: H,O, gene-B: C", text)
        self.assertIn("left out of the control (not qualifying records, never a reason to fail it): gene-B", text)
        p2 = dict(p, flagged_members=[], left_out=[])
        self.assertNotIn("flagged multiplicity", LC.failure_text(dict(positive=p2)))
        self.assertNotIn("flagged multiplicity", LC.failure_text(dict(positive={"pass": False, "missing": []})))                  # a result without the new keys

    # --- fix round 2, G2: what the negative control actually saw (no new threshold)
    def test_the_negative_control_reports_how_many_planted_records_got_a_triple_and_how_many_pairs_were_dropped(self):
        self.finish_controls(neg_ok=False)                                         # both members of both pairs on locus 1: two joined pairs
        planted = [r["id"] for r in self.w["ctl"][0]["records"]]
        self.add_copy_row("human", "ctl1", planted[0])                             # one member of the first pair has a second place: that pair's class is dropped
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        n = res["negative"]
        self.assertEqual((n["called"], n["of"], n["pass"]), (1, 2, True))           # one called pair is within the bar; the dropped one is not counted as called
        self.assertEqual((n["records_with_triple"], n["pairs_with_both_triples"], n["dropped_pairs"]), (4, 2, 1))

    def test_a_planted_pair_without_triples_is_not_counted_as_having_both(self):
        self.finish_controls()
        self.place("ctl2", {r["id"]: None for r in self.w["ctl"][1]["records"]})        # the second members are lifted nowhere
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        n = res["negative"]
        self.assertEqual((n["called"], n["pass"], n["records_with_triple"], n["pairs_with_both_triples"], n["dropped_pairs"]), (0, True, 2, 0, 0))

    def test_every_planted_pair_dropped_still_passes_but_says_so(self):
        self.finish_controls(neg_ok=False)
        for r in self.w["ctl"][0]["records"]:                                      # ctl1 holds one member of every pair
            self.add_copy_row("human", "ctl1", r["id"])
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        n = res["negative"]
        self.assertEqual((n["called"], n["pass"], n["dropped_pairs"], n["pairs_with_both_triples"]), (0, True, 2, 2))      # no threshold is invented: the numbers are reported

    def test_the_failure_text_of_the_negative_control_carries_the_counts(self):
        n = dict(called=2, of=50, bar=1, records_with_triple=97, pairs_with_both_triples=48, dropped_pairs=3, **{"pass": False})
        text = LC.failure_text(dict(negative=n))
        self.assertIn("negative control: 2 of 50 planted pairs came out as an expansion (bar 1)", text)
        self.assertIn("planted records with a triple 97", text)
        self.assertIn("pairs with both triples 48", text)
        self.assertIn("pairs dropped for multiplicity 3", text)
        self.assertEqual(LC.failure_text(dict(negative=dict(called=2, of=50, bar=1, **{"pass": False}))), "negative control: 2 of 50 planted pairs came out as an expansion (bar 1)")

    # --- fix round, addendum 2 (D-m2): one class that holds every qualifying record, on the wrong loci, fails the positive control and says which triple it has
    def place_lrpap_like_records(self, loci_by_outgroup):
        """the LRPAP1-like records of the positive runs on locus number loci_by_outgroup[key] in each outgroup, the fragments on locus 2; planted pairs as in finish_controls"""
        frag = {"gene-F1", "gene-F2", "gene-F3"}
        for run in self.w["pos"]:
            for key, n in loci_by_outgroup.items():
                self.place(run["label"], {r["id"]: (2 if r["id"] in frag else n) for r in run["records"]}, outgroups=key)
        for j, run in enumerate(self.w["ctl"]):
            self.place(run["label"], {r["id"]: (1 if j == 0 else 2) for r in run["records"]})

    def test_one_class_holding_every_qualifying_record_on_the_wrong_loci_fails_the_positive_control_and_names_both_triples(self):
        self.place_lrpap_like_records({"H": 3, "C": 3, "O": 3})
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        p = res["positive"]
        self.assertFalse(p["pass"])
        self.assertEqual((p["missing"], p["triple_ok"]), ([], False))                         # the class is whole; only its triple is wrong
        self.assertEqual((p["triple"], p["expected_triple"]), (["gene-H3", "gene-C3", "gene-O3"], ["gene-H1", "gene-C1", "gene-O1"]))
        self.assertEqual(sorted(p["class_members"]), ["gene-N", "gene-REF", "gene-X1", "gene-X2", "gene-X3"])
        self.assertTrue(res["negative"]["pass"])
        self.assertFalse(res["valid"])
        text = LC.failure_text(res)
        for needle in ("positive control", "gene-H3", "gene-O3", "gene-H1", "gene-O1", "expected"):
            self.assertIn(needle, text)
        self.assertNotIn("no single class holds every qualifying record", text)               # nothing is missing: the old message would mislead

    def test_a_triple_that_is_wrong_in_one_outgroup_only_fails_the_positive_control(self):
        self.place_lrpap_like_records({"H": 1, "C": 1, "O": 3})
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        self.assertEqual((res["positive"]["pass"], res["positive"]["triple_ok"], res["positive"]["triple"]), (False, False, ["gene-H1", "gene-C1", "gene-O3"]))
        self.assertFalse(res["valid"])

    def test_the_gate_stops_on_a_wrong_triple_and_the_saved_report_says_why(self):
        self.place_lrpap_like_records({"H": 3, "C": 3, "O": 3})
        with self.assertRaises(LC.StopRun) as cm:
            LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("neg")
        self.assertIn("Gate 7 FAILED before the neg group", str(cm.exception))
        self.assertIn("expected", str(cm.exception))
        with open(os.path.join(self.e_root, "gate7.json")) as fh:
            g = json.load(fh)
        self.assertEqual((g["positive"]["pass"], g["positive"]["triple_ok"], g["positive"]["expected_triple"]), (False, False, ["gene-H1", "gene-C1", "gene-O1"]))
        with self.assertRaises(LC.StopRun):                                                   # and it keeps stopping: the verdict group is not reachable either
            LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("verdict")

    def test_the_expected_triple_is_reported_when_the_positive_control_passes(self):
        self.finish_controls()
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        self.assertEqual((res["positive"]["pass"], res["positive"]["triple_ok"], res["positive"]["triple"], res["positive"]["expected_triple"]),
                         (True, True, ["gene-H1", "gene-C1", "gene-O1"], ["gene-H1", "gene-C1", "gene-O1"]))

    def test_the_failure_text_of_a_hand_built_result_keeps_the_missing_record_message(self):
        res = dict(positive={"pass": False, "missing": [dict(record="gene-X1", reason="no_placement")], "triple_ok": None})
        self.assertIn("missing [('gene-X1', 'no_placement')]", LC.failure_text(res))
        res = dict(positive={"pass": False, "missing": [], "triple_ok": False, "triple": ["a", "b", "c"], "expected_triple": ["x", "y", "z"]})
        self.assertIn("has the triple (a, b, c), not the expected (x, y, z)", LC.failure_text(res))

    def test_a_qualifying_record_without_a_placement_fails_the_positive_control(self):
        self.finish_controls()
        self.place("pos1", {"gene-REF": None, "gene-F1": 2})
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        self.assertFalse(res["positive"]["pass"])
        self.assertEqual(res["positive"]["missing"], [dict(record="gene-REF", reason="no_placement")])
        self.assertFalse(res["valid"])

    def test_a_fragment_that_takes_the_locus_of_a_full_length_copy_in_a_shared_run_fails_the_control(self):
        self.finish_controls()
        self.place("pos3", {"gene-X1": None, "gene-F3": 1})                     # the displaced copy is lifted elsewhere, i.e. not at all here
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        self.assertFalse(res["positive"]["pass"])
        self.assertEqual([m["record"] for m in res["positive"]["missing"]], ["gene-X1"])

    def test_two_planted_pairs_that_land_on_one_locus_fail_the_negative_control(self):
        self.finish_controls(neg_ok=False)
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        self.assertEqual((res["negative"]["called"], res["negative"]["pass"]), (2, False))
        self.assertTrue(res["positive"]["pass"])
        self.assertFalse(res["valid"])

    def test_one_group_can_be_evaluated_alone(self):
        self.finish_controls()
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes(), which=("pos",))
        self.assertIn("positive", res)
        self.assertNotIn("negative", res)

    def test_a_verdict_run_that_holds_a_control_record_does_not_rescue_a_failing_control(self):
        self.finish_controls()
        self.place("pos1", {"gene-REF": None, "gene-F1": 2})                    # REF is not placed in its control run ...
        self.place("r01", {"gene-V1": 3, "gene-REF": 1})                        # ... but a stray verdict run places it on the LRPAP1 locus
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        self.assertFalse(res["positive"]["pass"])
        self.assertEqual([m["record"] for m in res["positive"]["missing"]], ["gene-REF"])

    def test_a_control_run_that_was_split_is_read_through_its_halves(self):
        self.finish_controls()
        d = os.path.join(self.e_root, "lift", "human", "runs", "pos1")
        # pos1 holds gene-REF (c1) and gene-F1 (c5): two contigs, so a record-free cut exists; its halves replace it in the three outgroups
        for key, sp, loc in OG_SPECIES:
            LC.save_splits(self.e_root, sp, ["pos1"])
            make_run_dir(os.path.join(self.e_root, "lift", sp), "pos1.s1", {"gene-REF": (f"{loc}1", 1000, 3000)})
            make_run_dir(os.path.join(self.e_root, "lift", sp), "pos1.s2", {"gene-F1": (f"{loc}2", 1000, 3000)})
        import shutil
        for key, sp, loc in OG_SPECIES:
            shutil.rmtree(os.path.join(self.e_root, "lift", sp, "runs", "pos1"))        # created by this test
        res = LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())
        self.assertTrue(res["positive"]["pass"], res["positive"])

    def test_an_unfinished_control_run_is_a_stop(self):
        self.finish_controls()
        os.remove(os.path.join(self.e_root, "lift", "orangutan", "runs", "pos4", "done.key"))
        with self.assertRaises(LC.StopRun):
            LC.evaluate_controls(self.w, self.e_root, self.og, self.load_indexes())


class GateTests(WorldCase):
    def test_nothing_is_required_before_the_first_control(self):
        LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("pos")
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "gate7.json")))

    def test_an_unfinished_control_run_is_a_stop_not_a_verdict(self):
        self.finish_controls()
        os.remove(os.path.join(self.e_root, "lift", "chimpanzee", "runs", "ctl2", "done.key"))
        LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("neg")           # the positive control is complete
        with self.assertRaises(LC.StopRun) as cm:
            LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("verdict")
        self.assertIn("not finished", str(cm.exception))

    # --- fix round 2, G5: a hand-edited report with malformed keys is no saved report (it never was a decision), not a TypeError
    def test_a_hand_edited_report_with_malformed_keys_is_replaced_and_the_gate_still_decides(self):
        self.finish_controls()
        path = os.path.join(self.e_root, "gate7.json")
        for text in ('{"control_runs": "x"}', '{"control_runs": null}', '{"control_runs": [1, 2]}', '{"control_runs": 5, "positive": 5, "negative": [], "valid": "yes", "anomalies": 3}',
                     '{"positive": "p", "negative": 7}', '[]', '"s"', 'null', '{'):
            with open(path, "w") as fh:
                fh.write(text)
            LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("verdict")                    # the controls pass: no TypeError, no stop
            g = json.loads(rd(path))
            self.assertEqual(sorted(g["control_runs"]), ["neg", "pos"], text)
            self.assertIs(g["valid"], True, text)
            self.assertTrue(g["positive"]["pass"] and g["negative"]["pass"], text)

    def test_a_malformed_report_does_not_leave_junk_beside_the_group_that_was_evaluated(self):
        self.finish_controls()
        path = os.path.join(self.e_root, "gate7.json")
        with open(path, "w") as fh:
            fh.write('{"control_runs": 5, "positive": 5, "negative": []}')
        LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("neg")                          # only the positive control is evaluated by this gate
        g = json.loads(rd(path))
        self.assertEqual(sorted(g), ["anomalies", "control_runs", "positive", "valid"])
        self.assertTrue(g["positive"]["pass"])

    def test_a_malformed_report_does_not_hide_a_failed_control(self):
        self.finish_controls(neg_ok=False)
        with open(os.path.join(self.e_root, "gate7.json"), "w") as fh:
            fh.write('{"control_runs": "x"}')
        with self.assertRaises(LC.StopRun) as cm:
            LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("verdict")
        self.assertIn("negative", str(cm.exception))

    # --- fix round 2, G7: a finished run whose unmapped.txt is missing or damaged is a stop (exit 3), not a defect of the driver (exit 70)
    def test_a_finished_control_run_without_unmapped_txt_is_a_stop_of_the_gate(self):
        self.finish_controls()
        os.remove(os.path.join(self.e_root, "lift", "chimpanzee", "runs", "pos5", "unmapped.txt"))
        with self.assertRaises(LC.StopRun) as cm:
            LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("neg")
        self.assertIn("unmapped.txt", str(cm.exception))

    def test_a_finished_control_run_with_a_damaged_unmapped_txt_is_a_stop_of_the_gate(self):
        self.finish_controls()
        with open(os.path.join(self.e_root, "lift", "human", "runs", "pos5", "unmapped.txt"), "wb") as fh:
            fh.write(b"\xff\xfe\x00 not text\n")
        with self.assertRaises(LC.StopRun) as cm:
            LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("neg")
        self.assertIn("pos5", str(cm.exception))

    def test_a_finished_control_run_with_a_damaged_loci_tsv_is_a_stop_of_the_gate(self):
        self.finish_controls()
        with open(os.path.join(self.e_root, "lift", "orangutan", "runs", "pos2", "loci.tsv"), "w") as fh:
            fh.write("not\ta\tloci table\n1\t2\t3\n")
        with self.assertRaises(LC.StopRun) as cm:
            LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("neg")
        self.assertIn("cannot be read", str(cm.exception))

    def test_a_finished_control_run_without_loci_tsv_is_not_finished_and_stays_a_stop(self):
        self.finish_controls()
        os.remove(os.path.join(self.e_root, "lift", "orangutan", "runs", "pos2", "loci.tsv"))
        with self.assertRaises(LC.StopRun) as cm:
            LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("neg")
        self.assertIn("not finished", str(cm.exception))

    def test_the_negative_group_needs_the_positive_control_to_pass(self):
        self.finish_controls()
        gate = LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)
        gate("neg")
        self.place("pos1", {"gene-REF": None, "gene-F1": 2})
        with self.assertRaises(LC.StopRun):
            LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("neg")

    def test_the_verdict_group_needs_both_controls_to_pass(self):
        self.finish_controls(neg_ok=False)
        with self.assertRaises(LC.StopRun) as cm:
            LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("verdict")
        self.assertIn("negative", str(cm.exception))

    def test_the_gate_evaluates_the_controls_every_time_it_is_asked(self):
        # (replaces test_a_pass_is_remembered_until_a_control_run_changes: the saved verdict could be forged, and evaluating costs about 2 s)
        self.finish_controls()
        loads = []

        def loader():
            loads.append(1)
            return self.load_indexes()
        LC.make_gate(self.w, self.e_root, self.og, loader)("verdict")
        LC.make_gate(self.w, self.e_root, self.og, loader)("verdict")
        self.assertEqual(len(loads), 2)
        self.place("pos1", {"gene-REF": None, "gene-F1": 2})                      # an output that changed: the gate sees it
        with self.assertRaises(LC.StopRun):
            LC.make_gate(self.w, self.e_root, self.og, loader)("verdict")
        self.assertEqual(len(loads), 3)

    def test_a_forged_gate7_json_does_not_unlock_the_next_group(self):
        self.finish_controls()
        self.place("pos1", {"gene-REF": None, "gene-F1": 2})                      # the positive control failed
        with self.assertRaises(LC.StopRun):
            LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("neg")
        forged = dict(control_runs={"pos": LC.control_digest(self.w, self.e_root, self.og, "pos")}, positive={"pass": True}, valid=True)
        with open(os.path.join(self.e_root, "gate7.json"), "w") as fh:
            json.dump(forged, fh)
        with self.assertRaises(LC.StopRun):
            LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("neg")

    def test_the_gate_writes_the_verdict_it_used(self):
        self.finish_controls()
        LC.make_gate(self.w, self.e_root, self.og, self.load_indexes)("verdict")
        with open(os.path.join(self.e_root, "gate7.json")) as fh:
            g = json.load(fh)
        self.assertTrue(g["valid"])
        self.assertIn("control_runs", g)


class BuildAllTests(WorldCase):
    def build(self):
        return LC.build_all(self.e_root, os.path.join(self.tmp, "out"), lift_root=self.root, outgroups=self.og, hub_mode="nodes", spec=self.spec, n_neg=4)

    def table(self, name):
        with open(os.path.join(self.tmp, "out", name)) as fh:
            return list(csv.reader(fh, delimiter="\t"))

    def test_a_valid_table_holds_the_verdict_set_and_reports_lrpap1_separately(self):
        self.finish_controls()
        self.finish_verdict()
        res = self.build()
        self.assertTrue(res["gate7"]["valid"])
        self.assertEqual(rd(os.path.join(self.tmp, "out", "STATUS")).strip(), "VALID")
        et = self.table("etable.tsv")
        self.assertEqual([r[2] for r in et[1:]], ["c8:1000-3000", "c9:1000-3000"])                # V1, V2: one expansion at locus 3
        self.assertEqual({r[1] for r in et[1:]}, {"PC_gene-V1"})
        self.assertEqual({r[5:8] == ["gene-H3", "gene-C3", "gene-O3"] for r in et[1:]}, {True})
        with open(os.path.join(self.tmp, "out", "gate7.json")) as fh:
            g = json.load(fh)
        self.assertTrue(g["positive"]["pass"])
        self.assertEqual(sorted(g["positive"]["class_members"]), ["gene-N", "gene-REF", "gene-X1", "gene-X2", "gene-X3"])

    def test_sha256sums_cover_every_file_of_the_table_directory(self):
        import hashlib
        self.finish_controls()
        self.finish_verdict()
        self.build()
        out = os.path.join(self.tmp, "out")
        listed = {}
        for line in rd(os.path.join(out, "SHA256SUMS")).splitlines():
            h, _, name = line.partition("  ")
            listed[name] = h
        files = sorted(f for f in os.listdir(out) if f != "SHA256SUMS")
        self.assertEqual(sorted(listed), files)
        for name in files:
            with open(os.path.join(out, name), "rb") as fh:
                self.assertEqual(listed[name], hashlib.sha256(fh.read()).hexdigest(), name)
        self.assertIn("etable.tsv", listed)
        self.assertIn("gate7.json", listed)

    def test_the_control_expansions_are_written_apart_from_the_e_table(self):
        self.finish_controls()
        self.finish_verdict()
        self.build()
        ec = self.table("etable_control.tsv")
        self.assertEqual(ec[0], X.ETABLE_HEADER)
        self.assertEqual(sorted({r[0] for r in ec[1:]}), ["EC0001", "EC0002"])
        self.assertEqual(sorted(r[2] for r in ec[1:] if r[1] == "PC_gene-N"), ["c1:1000-3000", "c1:8000-9000", "c2:5000-7000", "c3:8000-30000", "c4:2500-3500"])
        self.assertEqual(sorted(r[2] for r in ec[1:] if r[1] == "PC_gene-F1"), ["c5:0-2000", "c6:0-2000", "c7:0-2000"])
        self.assertFalse({r[1] for r in self.table("etable.tsv")[1:]} & {"PC_gene-N", "PC_gene-F1"})

    def test_the_planted_pairs_that_were_called_are_listed_in_gate7(self):
        self.finish_controls(neg_ok=False)
        self.finish_verdict()
        res = self.build()
        called = res["gate7"]["negative"]["called_pairs"]
        self.assertEqual(len(called), 2)
        self.assertTrue(all(len(c["members"]) == 2 and c["component"].startswith("NEG") for c in called))

    def test_the_control_pool_comes_from_the_verdict_set_and_marks_eligibility(self):
        self.finish_controls()
        self.finish_verdict()
        self.build()
        cp = self.table("control_pool.tsv")
        self.assertEqual(cp[0], ["pool_component_id", "record_id", "reason", "joining_copy_rows", "eligible"])
        self.assertEqual([(r[0], r[1], r[2], r[4]) for r in cp[1:]],
                         [("PC_gene-W1", "c10:1000-3000", "different_triple", "1"), ("PC_gene-W1", "c11:1000-3000", "different_triple", "1")])

    def test_member_exons_and_extension_cover_the_whole_pool_including_the_control_components(self):
        self.finish_controls()
        self.finish_verdict()
        self.build()
        ids = {r[2] for r in self.table("member_exons.tsv")[1:]}
        self.assertTrue({"gene-REF", "gene-F1", "gene-V1", "gene-W2"} <= ids)

    def test_a_failing_positive_control_marks_the_table_invalid(self):
        self.finish_controls()
        self.place("pos1", {"gene-REF": None, "gene-F1": 2})
        self.finish_verdict()
        res = self.build()
        self.assertFalse(res["gate7"]["valid"])
        self.assertTrue(rd(os.path.join(self.tmp, "out", "STATUS")).startswith("INVALID"))
        self.assertIn("positive", rd(os.path.join(self.tmp, "out", "STATUS")))

    def test_a_failing_negative_control_marks_the_table_invalid(self):
        self.finish_controls(neg_ok=False)
        self.finish_verdict()
        res = self.build()
        self.assertFalse(res["gate7"]["valid"])
        self.assertIn("negative", rd(os.path.join(self.tmp, "out", "STATUS")))

    def test_unfinished_verdict_runs_are_an_error(self):
        self.finish_controls()
        self.place("r01", {"gene-V1": 3, "gene-W1": 4})
        with self.assertRaises(RuntimeError):
            self.build()

    def test_unfinished_control_runs_are_an_error(self):
        self.finish_verdict()
        with self.assertRaises(RuntimeError):
            self.build()

    # --- fix round: tripwires and the table directory
    def test_a_record_with_no_row_and_no_unmapped_listing_makes_the_table_invalid(self):
        self.finish_controls()
        self.finish_verdict()
        for key, sp, loc in OG_SPECIES:                                     # r01 of the world holds V1 and W1: V1 vanished from the human run (no row, not listed)
            if key == "H":
                make_run_dir(os.path.join(self.e_root, "lift", sp), "r01", {"gene-W1": (f"{loc}4", 1000, 3000)}, unmapped=[])
        res = self.build()
        self.assertFalse(res["valid"])
        status = rd(os.path.join(self.tmp, "out", "STATUS"))
        self.assertTrue(status.startswith("INVALID"), status)
        self.assertIn("H: no_row 1", status)                                # the anomaly names its outgroup and its kind

    def test_a_record_with_a_copy_row_and_no_primary_row_makes_the_table_invalid(self):
        self.finish_controls()
        self.finish_verdict()
        import pathlib
        import _liftoff
        d = os.path.join(self.e_root, "lift", "chimpanzee", "runs", "r01")
        with open(os.path.join(d, "lifted.gff3"), "w") as fh:
            fh.write(lifted_gff([("gene-V1_1", "c3", 1001, 3000, "+", "0.99", "0.98", 1, "False", "False", [(1001, 3000)], "gene"),
                                 ("gene-W1", "c4", 1001, 3000, "+", "0.99", "0.98", 0, "False", "False", [(1001, 3000)], "gene")]))
        _liftoff.parse_lifted(pathlib.Path(d) / "lifted.gff3", pathlib.Path(d) / "loci.tsv")
        res = self.build()
        self.assertFalse(res["valid"])
        self.assertIn("C: no_primary_row 1", rd(os.path.join(self.tmp, "out", "STATUS")))

    def test_a_verdict_class_with_a_reverse_copy_row_is_dropped_counted_and_reported(self):
        self.finish_controls()
        self.finish_verdict()
        self.add_copy_row("human", "r02", "gene-V2")                       # V2 also has a second place in the human genome: the class {V1, V2} is dropped as a whole
        res = self.build()
        self.assertTrue(res["valid"])                                       # a dropped class is a count, not an anomaly
        self.assertEqual((res["expansions"], res["dropped_classes"], res["records_in_dropped_classes"]), (0, 1, 2))
        self.assertEqual(self.table("etable.tsv")[1:], [])
        self.assertEqual(self.table("dropped_classes.tsv")[1:], [["PC_gene-V1", "gene-H3", "gene-C3", "gene-O3", "gene-V1,gene-V2", "multiplicity"]])
        s = json.loads(rd(os.path.join(self.tmp, "out", "summary.json")))
        self.assertEqual((s["dropped_classes"], s["records_in_dropped_classes"], s["records_with_reverse_copy_rows"]), (1, 2, {"H": 1, "C": 0, "O": 0}))
        calls = {r[1]: r for r in self.table("record_calls.tsv")[1:]}
        self.assertEqual(calls["gene-V2"][-1], "H")
        self.assertEqual(calls["gene-V1"][-1], "")

    def test_the_flagged_record_of_a_dropped_verdict_class_is_in_the_control_pool_as_outgroup_copy_row(self):
        # fix round, addendum 2: without the copy row the component V is a plain expansion and is not in the pool; with it the class is dropped and V2 is listed
        self.finish_controls()
        self.finish_verdict()
        self.build()
        pool_of = lambda: [(r[0], r[1], r[2], r[4]) for r in self.table("control_pool.tsv")[1:] if r[2] == "outgroup_copy_row"]    # noqa: E731
        self.assertEqual(pool_of(), [])
        self.add_copy_row("human", "r02", "gene-V2")
        self.build()                                                       # the same table directory: a build overwrites its own reading
        self.assertEqual(pool_of(), [("PC_gene-V1", "c9:1000-3000", "outgroup_copy_row", "1")])
        s = json.loads(rd(os.path.join(self.tmp, "out", "summary.json")))
        self.assertEqual(s["control_pool_rows"], len(self.table("control_pool.tsv")) - 1)
        self.assertEqual(s["control_pool_rows_eligible"], sum(int(r[4]) for r in self.table("control_pool.tsv")[1:]))

    def test_a_positive_control_class_on_another_triple_makes_the_table_invalid_and_says_which(self):
        frag = {"gene-F1", "gene-F2", "gene-F3"}
        for run in self.w["pos"]:
            self.place(run["label"], {r["id"]: (2 if r["id"] in frag else 3) for r in run["records"]})                # the LRPAP1-like records on locus 3: one class, wrong loci
        for j, run in enumerate(self.w["ctl"]):
            self.place(run["label"], {r["id"]: (1 if j == 0 else 2) for r in run["records"]})
        self.finish_verdict()
        res = self.build()
        self.assertFalse(res["valid"])
        status = rd(os.path.join(self.tmp, "out", "STATUS"))
        self.assertTrue(status.startswith("INVALID: positive control"), status)
        self.assertIn("gene-H3", status)
        self.assertIn("gene-H1", status)
        self.assertNotIn("no single class holds every qualifying record", status)
        g = json.loads(rd(os.path.join(self.tmp, "out", "gate7.json")))
        self.assertEqual((g["positive"]["pass"], g["positive"]["triple_ok"], g["positive"]["expected_triple"]), (False, False, ["gene-H1", "gene-C1", "gene-O1"]))
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "out", "etable.tsv")))                                  # the tables are written all the same

    # --- fix round 2, G3: after a gate stop the neg / verdict runs were never lifted; the build says which control failed instead of "N runs not finished"
    def finish_pos_only(self, broken=False):
        """the state a gate stop leaves: the pos group lifted (broken: gene-REF placed nowhere), neg and verdict never started"""
        frag = {"gene-F1", "gene-F2", "gene-F3"}
        for run in self.w["pos"]:
            self.place(run["label"], {r["id"]: (2 if r["id"] in frag else 1) for r in run["records"]})
        if broken:
            self.place("pos1", {"gene-REF": None, "gene-F1": 2})

    def finish_neg_only(self, neg_ok=True):
        for j, run in enumerate(self.w["ctl"]):
            self.place(run["label"], {r["id"]: (1 if (j == 0 or not neg_ok) else 2) for r in run["records"]})

    def partial_dir(self):
        return os.path.join(self.tmp, "out")

    def test_after_a_gate_stop_on_the_positive_control_the_build_names_the_control_and_writes_gate7_only(self):
        self.finish_pos_only(broken=True)
        res = self.build()
        self.assertFalse(res["valid"])
        self.assertTrue(res["status"].startswith("INVALID: positive control"), res["status"])
        self.assertNotIn("not finished", res["status"])
        out = self.partial_dir()
        self.assertEqual(sorted(os.listdir(out)), ["SHA256SUMS", "STATUS", "gate7.json"])
        self.assertEqual(rd(os.path.join(out, "STATUS")).strip(), res["status"])
        g = json.loads(rd(os.path.join(out, "gate7.json")))
        self.assertEqual((g["valid"], g["positive"]["pass"], "negative" in g), (False, False, False))
        self.assertEqual(sorted(g["control_runs"]), ["pos"])
        self.assertEqual((res["expansions"], res["dropped_classes"], res["records_in_dropped_classes"]), (0, 0, 0))
        import hashlib
        for line in rd(os.path.join(out, "SHA256SUMS")).splitlines():
            h, _, name = line.partition("  ")
            self.assertEqual(h, hashlib.sha256(open(os.path.join(out, name), "rb").read()).hexdigest(), name)

    def test_after_a_gate_stop_on_the_negative_control_the_build_names_the_control_and_writes_gate7_only(self):
        self.finish_pos_only()
        self.finish_neg_only(neg_ok=False)
        res = self.build()
        self.assertFalse(res["valid"])
        self.assertTrue(res["status"].startswith("INVALID: negative control: 2 of 2"), res["status"])
        self.assertNotIn("positive control", res["status"])
        self.assertEqual(sorted(os.listdir(self.partial_dir())), ["SHA256SUMS", "STATUS", "gate7.json"])
        g = json.loads(rd(os.path.join(self.partial_dir(), "gate7.json")))
        self.assertEqual((g["valid"], g["positive"]["pass"], g["negative"]["pass"], g["negative"]["called"]), (False, True, False, 2))
        self.assertEqual(sorted(g["control_runs"]), ["neg", "pos"])

    def test_a_failed_positive_control_and_finished_negative_runs_are_both_reported(self):
        self.finish_pos_only(broken=True)
        self.finish_neg_only(neg_ok=False)
        res = self.build()
        self.assertTrue(res["status"].startswith("INVALID: positive control"), res["status"])
        self.assertIn("negative control: 2 of 2", res["status"])
        self.assertEqual(sorted(os.listdir(self.partial_dir())), ["SHA256SUMS", "STATUS", "gate7.json"])

    def test_anomalies_in_the_only_finished_group_are_reported_as_invalid_too(self):
        self.finish_pos_only()
        for key, sp, loc in OG_SPECIES:                                            # pos5 holds gene-X3: it vanished from the human run (no row, not listed)
            if key == "H":
                make_run_dir(os.path.join(self.e_root, "lift", sp), "pos5", {}, unmapped=[])
        res = self.build()
        self.assertFalse(res["valid"])
        self.assertIn("H: no_row 1", res["status"])
        self.assertEqual(sorted(os.listdir(self.partial_dir())), ["SHA256SUMS", "STATUS", "gate7.json"])

    def test_unfinished_negative_runs_after_a_passing_positive_control_are_still_a_stop(self):
        self.finish_pos_only()
        with self.assertRaises(LC.StopRun) as cm:
            self.build()
        self.assertIn("not finished", str(cm.exception))
        self.assertFalse(os.path.exists(self.partial_dir()))

    def test_unfinished_verdict_runs_after_two_good_controls_are_still_a_stop(self):
        self.finish_controls()
        with self.assertRaises(LC.StopRun) as cm:
            self.build()
        self.assertIn("not finished", str(cm.exception))
        self.assertFalse(os.path.exists(self.partial_dir()))

    def test_unfinished_positive_runs_are_a_stop_even_when_nothing_else_exists(self):
        with self.assertRaises(LC.StopRun) as cm:
            self.build()
        self.assertIn("not finished", str(cm.exception))

    def test_a_failed_control_with_every_run_finished_still_writes_the_whole_invalid_table(self):
        self.finish_controls(neg_ok=False)
        self.finish_verdict()
        res = self.build()
        self.assertFalse(res["valid"])
        self.assertIn("etable.tsv", os.listdir(self.partial_dir()))
        self.assertTrue(res["status"].startswith("INVALID: negative control"))

    # --- fix round 2, G7
    def test_a_finished_run_without_unmapped_txt_is_a_stop_of_the_build_not_a_defect(self):
        self.finish_controls()
        self.finish_verdict()
        os.remove(os.path.join(self.e_root, "lift", "chimpanzee", "runs", "r01", "unmapped.txt"))
        with self.assertRaises(LC.StopRun) as cm:
            self.build()
        self.assertIn("unmapped.txt", str(cm.exception))

    def test_a_finished_run_with_a_damaged_unmapped_txt_is_a_stop_of_the_build_not_a_defect(self):
        self.finish_controls()
        self.finish_verdict()
        with open(os.path.join(self.e_root, "lift", "human", "runs", "r02", "unmapped.txt"), "wb") as fh:
            fh.write(b"\xff\xfe\x00 not text\n")
        with self.assertRaises(LC.StopRun):
            self.build()

    def test_gate7_json_and_the_summary_carry_the_diagnostics_of_the_two_controls(self):
        # fix round 2, G1 and G2
        self.finish_controls()
        self.finish_verdict()
        self.add_copy_row("human", "pos4", "gene-X2")                                  # a not-required record flagged: left out of the control, which passes
        res = self.build()
        self.assertTrue(res["valid"])
        g = json.loads(rd(os.path.join(self.tmp, "out", "gate7.json")))
        self.assertEqual(g["positive"]["flagged_members"], [dict(record="gene-X2", outgroups=["H"], qualifying=False, left_out=True)])
        self.assertEqual(g["positive"]["left_out"], ["gene-X2"])
        self.assertEqual((g["negative"]["records_with_triple"], g["negative"]["pairs_with_both_triples"], g["negative"]["dropped_pairs"]), (4, 2, 0))
        s = json.loads(rd(os.path.join(self.tmp, "out", "summary.json")))
        self.assertEqual((s["negative_records_with_triple"], s["negative_pairs_with_both_triples"], s["negative_dropped_pairs"]), (4, 2, 0))

    def test_a_valid_build_is_valid(self):
        self.finish_controls()
        self.finish_verdict()
        self.assertTrue(self.build()["valid"])

    def test_build_refuses_a_directory_that_holds_a_table_of_another_reading(self):
        self.finish_controls()
        self.finish_verdict()
        self.build()
        with self.assertRaises(LC.StopRun) as cm:
            LC.build_all(self.e_root, os.path.join(self.tmp, "out"), lift_root=self.root, outgroups=self.og, hub_mode="nodes", spec=self.spec, n_neg=4,
                         strict=True)
        self.assertIn("strict", str(cm.exception))
        with self.assertRaises(LC.StopRun):
            LC.build_all(self.e_root, os.path.join(self.tmp, "out"), lift_root=self.root, outgroups=self.og, hub_mode="records", spec=self.spec, n_neg=4)


STUB_CLI_LIFTOFF = r"""#!/usr/bin/env python3
import json, os, re, sys
args = sys.argv[1:]
target, out, un = args[0], args[args.index("-o") + 1], args[args.index("-u") + 1]
species = target.split("/lift/")[1].split("/")[0]
placed = json.load(open(os.environ["STUB_MAP"]))[species]
ids = []
for line in open("in.gff"):
    f = line.rstrip("\n").split("\t")
    if len(f) == 9 and f[2] in ("gene", "pseudogene") and "Parent=" not in f[8]:
        ids.append(re.search(r"ID=([^;]+)", f[8]).group(1))
lines, unmapped = ["##gff-version 3"], []
for rid in ids:
    loc = placed.get(rid)
    if loc is None:
        unmapped.append(rid)
        continue
    contig, s0, e = loc
    lines.append(f"{contig}\tLiftoff\tgene\t{s0 + 1}\t{e}\t.\t+\t.\tID={rid};coverage=0.99;sequence_ID=0.98;extra_copy_number=0;copy_num_ID={rid}_0;partial_mapping=False;low_identity=False")
    lines.append(f"{contig}\tLiftoff\texon\t{s0 + 1}\t{e}\t.\t+\t.\tID=exon-{rid}-1;Parent={rid}")
open(out, "w").write("\n".join(lines) + "\n")
open(un, "w").write("".join(f"{u}\n" for u in unmapped))
"""


class CliCase(WorldCase):
    """main() driven call after call with a stub Liftoff that places every record where STUB_MAP says (the harness of the CLI tests)"""

    def setUp(self):
        super().setUp()
        self.fasta = os.path.join(self.tmp, "gorilla", "gorilla.fa")
        make_fake_genome(self.fasta, mmi=False)
        self.gff = os.path.join(self.tmp, "gorilla.gff")
        with open(self.gff, "w") as fh:
            fh.write("##gff-version 3\n")
            for rid, r in self.recs.items():
                fh.write(f"{r['contig']}\tRefSeq\tgene\t{r['start0'] + 1}\t{r['end']}\t.\t+\t.\tID={rid};Name={rid}\n")
        self.stub = os.path.join(self.tmp, "stub_liftoff.py")
        with open(self.stub, "w") as fh:
            fh.write(STUB_CLI_LIFTOFF)
        os.chmod(self.stub, os.stat(self.stub).st_mode | stat.S_IXUSR)
        self.map_path = os.path.join(self.tmp, "stub_map.json")
        self.placed = self.good_map()
        patches = [mock.patch.dict(LC.OUTGROUPS, {k: dict(t, assembly="x") for k, t in self.og.items()}, clear=True),
                   mock.patch.object(X, "LRPAP1_CONTROL", self.spec),
                   mock.patch.dict(os.environ, {"STUB_MAP": self.map_path})]
        for pt in patches:
            pt.start()
            self.addCleanup(pt.stop)

    def good_map(self):
        frag = {"gene-F1", "gene-F2", "gene-F3"}
        locus = {"gene-V1": 3, "gene-V2": 3, "gene-W1": 4, "gene-W2": 2}
        for run in self.w["pos"]:
            for r in run["records"]:
                locus[r["id"]] = 2 if r["id"] in frag else 1
        for j, run in enumerate(self.w["ctl"]):
            for r in run["records"]:
                locus[r["id"]] = 1 if j == 0 else 2
        return {sp: {rid: [f"{loc}{n}", 1000, 3000] for rid, n in locus.items()} for _, sp, loc in OG_SPECIES}

    def write_map(self):
        with open(self.map_path, "w") as fh:
            json.dump(self.placed, fh)

    def cli(self, cmd, *extra, out=None, liftoff_bin=None, minimap2="/m/mm2", rehearsal=True):
        argv = [cmd, "--out", out or self.e_root, "--lift-root", self.root, "--n-planted", "4", *extra]
        if cmd == "run":
            argv += ["--amendment-commit", "abc1234", "--no-code-check", "--gorilla-gff", self.gff, "--gorilla-fasta", self.fasta,
                     "--liftoff-bin", liftoff_bin or self.stub]
            if minimap2:                            # None = the default per-call cache wrapper (the real default path)
                argv += ["--minimap2", minimap2]
        if cmd == "build" and rehearsal:            # the runs were made by a stub: `build` refuses them unless it is told this is a rehearsal
            argv.append("--rehearsal")
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = LC.main(argv)
        return rc, out.getvalue(), err.getvalue()

    def run_until_done(self, limit=8):
        seen = []
        for _ in range(limit):
            rc, out, err = self.cli("run")
            seen.append((rc, json.loads(out)["group"] if out.strip() else None))
            if rc != 75:
                return seen, err
        return seen, ""

    def run_dirs(self, species="human"):
        base = os.path.join(self.e_root, "lift", species, "runs")
        return sorted(os.listdir(base)) if os.path.isdir(base) else []


class CliEndToEndTests(CliCase):
    """the groups, the gate, the exit codes, the build"""

    def test_three_calls_controls_then_verdict_then_a_valid_table(self):
        self.write_map()
        seen, _ = self.run_until_done()
        self.assertEqual(seen, [(75, "pos"), (75, "neg"), (0, "verdict")])
        self.assertEqual(self.run_dirs(), ["ctl1", "ctl2", "pos1", "pos2", "pos3", "pos4", "pos5", "r01", "r02"])
        rc, out, err = self.cli("build")
        self.assertEqual(rc, 0, err)
        self.assertEqual(rd(os.path.join(self.e_root, "table", "STATUS")).strip(), "VALID")
        with open(os.path.join(self.e_root, "table", "etable.tsv")) as fh:
            self.assertEqual({r[2] for r in list(csv.reader(fh, delimiter="\t"))[1:]}, {"c8:1000-3000", "c9:1000-3000"})

    def test_a_failing_positive_control_stops_the_next_call_before_any_other_run(self):
        for sp in self.placed:
            self.placed[sp]["gene-REF"] = None                                  # the LRPAP1 reference is placed nowhere
        self.write_map()
        seen, err = self.run_until_done()
        self.assertEqual(seen, [(75, "pos"), (3, None)])
        self.assertIn("Gate 7 FAILED", err)
        self.assertIn("positive", err)
        self.assertEqual(self.run_dirs(), ["pos1", "pos2", "pos3", "pos4", "pos5"])

    def test_a_positive_control_class_on_the_wrong_loci_stops_the_next_call_naming_both_triples(self):
        # fix round, addendum 2 (D-m2): every LRPAP1-like record is placed, all on locus 3 of each outgroup: ONE class holds them all, but it is not the expected one
        frag = {"gene-F1", "gene-F2", "gene-F3"}
        for _, sp, loc in OG_SPECIES:
            for run in self.w["pos"]:
                for r in run["records"]:
                    if r["id"] not in frag:
                        self.placed[sp][r["id"]] = [f"{loc}3", 1000, 3000]
        self.write_map()
        seen, err = self.run_until_done()
        self.assertEqual(seen, [(75, "pos"), (3, None)])
        self.assertIn("Gate 7 FAILED", err)
        for needle in ("gene-H3", "gene-O3", "gene-H1", "gene-O1", "expected"):
            self.assertIn(needle, err)
        self.assertNotIn("no single class holds every qualifying record", err)
        self.assertEqual(self.run_dirs(), ["pos1", "pos2", "pos3", "pos4", "pos5"])
        with open(os.path.join(self.e_root, "gate7.json")) as fh:
            self.assertIs(json.load(fh)["positive"]["triple_ok"], False)

    # --- fix round 2, G3 and G7 on the command line
    def test_a_build_after_a_gate_stop_says_which_control_failed_and_exits_4(self):
        for sp in self.placed:
            self.placed[sp]["gene-REF"] = None                                  # the LRPAP1 reference is placed nowhere
        self.write_map()
        seen, err = self.run_until_done()
        self.assertEqual(seen, [(75, "pos"), (3, None)])
        rc, out, err = self.cli("build")
        self.assertEqual(rc, 4, err)
        self.assertTrue(json.loads(out)["status"].startswith("INVALID: positive control"), out)
        self.assertNotIn("not finished", out + err)
        self.assertEqual(sorted(os.listdir(os.path.join(self.e_root, "table"))), ["SHA256SUMS", "STATUS", "gate7.json"])
        self.assertTrue(rd(os.path.join(self.e_root, "table", "STATUS")).startswith("INVALID: positive control"))

    def test_a_build_after_a_gate_stop_on_the_negative_control_says_so(self):
        for sp in self.placed:
            for r in self.w["ctl"][1]["records"]:
                self.placed[sp][r["id"]] = self.placed[sp][self.w["ctl"][0]["records"][0]["id"]]       # every planted record on one locus
        self.write_map()
        seen, err = self.run_until_done()
        self.assertEqual(seen, [(75, "pos"), (75, "neg"), (3, None)])
        rc, out, err = self.cli("build")
        self.assertEqual(rc, 4, err)
        self.assertTrue(json.loads(out)["status"].startswith("INVALID: negative control"), out)

    def test_a_build_with_a_finished_run_that_lacks_unmapped_txt_is_exit_3_not_70(self):
        self.write_map()
        seen, _ = self.run_until_done()
        self.assertEqual(seen[-1][0], 0, seen)
        os.remove(os.path.join(self.e_root, "lift", "chimpanzee", "runs", "pos5", "unmapped.txt"))
        rc, out, err = self.cli("build")
        self.assertEqual(rc, 3, err)
        self.assertIn("STOP:", err)
        self.assertIn("unmapped.txt", err)
        self.assertNotIn("Traceback", err)

    def test_a_failing_negative_control_stops_the_call_that_would_start_the_verdict_set(self):
        for sp in self.placed:
            for r in self.w["ctl"][1]["records"]:
                self.placed[sp][r["id"]] = self.placed[sp][self.w["ctl"][0]["records"][0]["id"]]       # every planted record on one locus
        self.write_map()
        seen, err = self.run_until_done()
        self.assertEqual(seen, [(75, "pos"), (75, "neg"), (3, None)])
        self.assertIn("negative", err)
        self.assertEqual(self.run_dirs(), ["ctl1", "ctl2", "pos1", "pos2", "pos3", "pos4", "pos5"])
        with open(os.path.join(self.e_root, "gate7.json")) as fh:
            self.assertFalse(json.load(fh)["negative"]["pass"])

    def test_another_hub_reading_in_a_later_call_is_a_stop(self):
        self.write_map()
        self.cli("run")
        rc, out, err = self.cli("run", "--hub-mode", "records")
        self.assertEqual(rc, 3)
        self.assertIn("plan digest", err)

    def test_the_plan_command_prints_the_digest_and_the_groups(self):
        rc, out, err = self.cli("plan")
        self.assertEqual(rc, 0, err)
        d = json.loads(out)
        self.assertEqual(sorted(d["groups"]), ["neg", "pos", "verdict"])
        self.assertEqual(d["groups"]["verdict"]["units"], 6)
        self.assertEqual(len(d["plan_digest"]), 64)
        self.assertEqual(d["verdict"]["components"], 2)

    # --- fix round, blocker 1: the first call on a fresh root (nothing in it) with the DEFAULT minimap2 wrapper
    def test_the_first_call_on_a_fresh_root_with_the_default_wrapper_runs(self):
        self.write_map()
        fresh = os.path.join(self.tmp, "fresh_root")                      # no lift/, no plan.digest: like the real E_ROOT before the first call
        rc, out, err = self.cli("run", out=fresh, minimap2=None)          # no --minimap2: the per-call cache wrapper written into E_ROOT/lift
        self.assertEqual(rc, 75, err)
        wrapper = os.path.join(fresh, "lift", "mm2_cache.sh")
        self.assertTrue(os.access(wrapper, os.X_OK), "the driver creates E_ROOT/lift before it writes the wrapper there")
        self.assertEqual(json.loads(out)["group"], "pos")

    def test_a_call_that_cannot_start_leaves_no_frozen_state(self):
        self.write_map()
        fresh = os.path.join(self.tmp, "fresh_root2")
        rc, out, err = self.cli("run", out=fresh, liftoff_bin=os.path.join(self.tmp, "no_such_liftoff"))
        self.assertEqual(rc, 3, err)
        self.assertIn("liftoff", err)
        for name in ("plan.digest", "run_identity.json"):                  # a crash before the first Liftoff process must not trip the freeze of stop rule 5
            self.assertFalse(os.path.exists(os.path.join(fresh, name)), name)
        rc, out, err = self.cli("run", out=fresh)                          # the corrected call is then the first call
        self.assertEqual(rc, 75, err)
        self.assertTrue(os.path.exists(os.path.join(fresh, "plan.digest")))
        self.assertTrue(os.path.exists(os.path.join(fresh, "run_identity.json")))

    def test_a_missing_outgroup_file_is_a_stop_before_anything_is_frozen(self):
        self.write_map()
        fresh = os.path.join(self.tmp, "fresh_root3")
        os.remove(self.og["O"]["mmi"])                                      # the Liftoff index of the orangutan is gone
        rc, out, err = self.cli("run", out=fresh)
        self.assertEqual(rc, 3, err)
        self.assertIn(self.og["O"]["mmi"], err)
        for name in ("plan.digest", "run_identity.json"):
            self.assertFalse(os.path.exists(os.path.join(fresh, name)), name)
        self.assertFalse(os.path.isdir(os.path.join(fresh, "lift", "human", "runs")))

    def test_a_missing_timeout_binary_is_a_stop_before_anything_is_frozen(self):
        self.write_map()
        fresh = os.path.join(self.tmp, "fresh_root5")
        with mock.patch.object(LC, "TIMEOUT_BIN", os.path.join(self.tmp, "no_such_timeout")):
            rc, out, err = self.cli("run", out=fresh)
        self.assertEqual(rc, 3, err)
        self.assertIn("no_such_timeout", err)
        self.assertFalse(os.path.exists(os.path.join(fresh, "plan.digest")))

    def test_a_missing_reference_fasta_is_a_stop_before_anything_is_frozen(self):
        self.write_map()
        fresh = os.path.join(self.tmp, "fresh_root4")
        os.remove(self.fasta + ".fai")
        rc, out, err = self.cli("run", out=fresh)
        self.assertEqual(rc, 3, err)
        self.assertIn(".fai", err)
        self.assertFalse(os.path.exists(os.path.join(fresh, "plan.digest")))

    # --- fix round, stop rule 6 through the command line
    def test_the_command_line_stops_with_exit_3_when_the_time_budget_of_an_outgroup_is_passed(self):
        self.write_map()
        with mock.patch.object(LC, "WALL_FACTOR", 1e-9):
            rc, out, err = self.cli("run")
            self.assertEqual(rc, 3, err)
            self.assertIn("time budget", err)
            self.assertEqual(self.run_dirs(), ["pos1"])                     # the run that crossed the limit finished; no other started
            rc, out, err = self.cli("run")                                  # the sum is remembered in E_ROOT/wall.json: nothing new starts
            self.assertEqual(rc, 3, err)
            self.assertEqual(self.run_dirs(), ["pos1"])
        self.assertEqual(len(LC.WallLedger(self.e_root).entries("H")), 1)

    def test_the_plan_command_prints_the_predicted_wall_and_the_limit_of_every_outgroup(self):
        rc, out, err = self.cli("plan")
        d = json.loads(out)
        self.assertEqual(sorted(d["wall_predicted_s"]), ["C", "H", "O"])
        by_label = {r["label"]: r for r in self.w["plan"] + self.w["pos"] + self.w["ctl"]}
        want = sum(LC.predict_wall(len(by_label[lab]["records"])) for lab in by_label)
        self.assertAlmostEqual(d["wall_predicted_s"]["H"], want, delta=1)
        self.assertAlmostEqual(d["wall_limit_s"]["H"], 2 * want, delta=1)

    # --- fix round: gffs writes the marker so that the first heavy call does not rewrite the 60 GFFs
    def test_gffs_writes_the_assignment_marker_and_the_first_run_does_not_rewrite(self):
        rc, out, err = self.cli("gffs", "--gorilla-gff", self.gff)
        self.assertEqual(rc, 0, err)
        mark = os.path.join(self.e_root, "gff", "assign.sha256")
        self.assertTrue(os.path.exists(mark))
        first = os.path.getmtime(os.path.join(self.e_root, "gff", "pos1", "in.gff"))
        self.write_map()
        rc, out, err = self.cli("run")
        self.assertEqual(os.path.getmtime(os.path.join(self.e_root, "gff", "pos1", "in.gff")), first)

    def test_the_run_gffs_are_rewritten_when_the_annotation_changed(self):
        assign = {"gene-A": "r01"}
        gff = os.path.join(self.tmp, "ref.gff")
        with open(gff, "w") as fh:
            fh.write(REF_GFF)
        out = os.path.join(self.tmp, "gffout")
        self.assertTrue(LC.ensure_run_gffs(gff, assign, out))
        self.assertFalse(LC.ensure_run_gffs(gff, assign, out))
        with open(gff, "a") as fh:
            fh.write("NC_1\tRefSeq\tregion\t1\t10\t.\t+\t.\tID=x\n")            # the annotation file changed under the same assignment
        self.assertTrue(LC.ensure_run_gffs(gff, assign, out))

    # --- fix round: a build never overwrites a table, and a strict build has a table directory of its own
    def finish_all(self):
        self.write_map()
        seen, _ = self.run_until_done()
        self.assertEqual(seen[-1][0], 0, seen)

    def test_a_strict_build_needs_a_table_dir_of_its_own(self):
        self.finish_all()
        rc, out, err = self.cli("build", "--strict-single-copy")
        self.assertEqual(rc, 2, err)
        self.assertIn("table-dir", err)
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "table")))
        rc, out, err = self.cli("build", "--strict-single-copy", "--table-dir", os.path.join(self.e_root, "table"))
        self.assertEqual(rc, 2, err)
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "table")))
        rc, out, err = self.cli("build", "--strict-single-copy", "--table-dir", os.path.join(self.e_root, "table_strict"))
        self.assertEqual(rc, 0, err)
        self.assertTrue(json.loads(rd(os.path.join(self.e_root, "table_strict", "summary.json")))["strict_single_copy"])

    def test_a_build_never_writes_into_a_table_directory_that_is_not_empty(self):
        self.finish_all()
        self.assertEqual(self.cli("build")[0], 0)
        before = rd(os.path.join(self.e_root, "table", "SHA256SUMS"))
        rc, out, err = self.cli("build")
        self.assertEqual(rc, 2, err)
        self.assertIn("not empty", err)
        self.assertEqual(rd(os.path.join(self.e_root, "table", "SHA256SUMS")), before)

    def test_the_build_output_counts_the_classes_dropped_for_multiplicity(self):
        self.finish_all()
        path = os.path.join(self.e_root, "lift", "human", "runs", "r02", "loci.tsv")
        with open(path, "a") as fh:                                         # gene-V2 has a second place in the human genome
            fh.write("\t".join(map(str, ["gene-V2_1", "gene-V2_1", 1, "gene", "x9", 0, 800, "+", 0.97, 0.97, 0, 0, "0-800", 800])) + "\n")
        rc, out, err = self.cli("build")
        self.assertEqual(rc, 0, err)
        printed = json.loads(out)
        self.assertEqual((printed["expansions"], printed["dropped_classes"], printed["records_in_dropped_classes"]), (0, 1, 2))
        self.assertEqual(printed["status"], "VALID")

    def test_build_prints_only_the_status_when_a_control_failed(self):
        self.finish_all()
        # the controls of the finished runs are made to fail afterwards: the positive control (every LRPAP1-like record unplaced) -> the table is INVALID
        for sp in ("human", "chimpanzee", "orangutan"):
            d = os.path.join(self.e_root, "lift", sp, "runs", "pos1")
            with open(os.path.join(d, "loci.tsv"), "w") as fh:
                fh.write("\t".join(LC.L.LOCI_HEAD) + "\n")
            with open(os.path.join(d, "unmapped.txt"), "w") as fh:
                fh.write("gene-REF\ngene-F1\n")
        rc, out, err = self.cli("build")
        self.assertEqual(rc, 4, err)
        printed = json.loads(out)
        self.assertTrue(printed["status"].startswith("INVALID"))
        self.assertNotIn("expansions", printed)

    # --- fix round: the exit codes
    def test_a_defect_of_the_driver_is_exit_70_with_a_traceback_not_a_stop(self):
        with mock.patch.object(LC, "plan_from_real_inputs", side_effect=ValueError("boom")):
            rc, out, err = self.cli("plan")
        self.assertEqual(rc, 70)
        self.assertIn("Traceback", err)
        self.assertIn("boom", err)

    def test_a_stop_rule_in_plan_is_exit_3(self):
        with mock.patch.object(LC, "plan_from_real_inputs", side_effect=LC.StopRun("no")):
            rc, out, err = self.cli("plan")
        self.assertEqual(rc, 3)
        self.assertIn("STOP: no", err)

    def test_main_hands_the_start_of_the_process_to_the_driver(self):
        seen = {}

        def fake(units, plan_by_label, **kw):
            seen.update(kw)
            return {"done": [], "pending": [], "split": None, "deadline": None, "group": None, "exit_code": 0}
        with mock.patch.object(LC, "run_budgeted", fake):
            rc, out, err = self.cli("run")
        self.assertEqual(rc, 0, err)
        self.assertEqual(seen["t0"], LC.CALL_T0)
        self.assertIn("wall_pred_s", seen)
        self.assertEqual(len(seen["pass_fds"]), 1)


# ----------------------------------------------------------------------------------------------- group G: the fix round (review of 2026-10-07)
def ps_session(sid):
    """live (not zombie) pids of the session `sid`"""
    out = subprocess.run(["ps", "-o", "pid=,stat=", "--sid", str(sid)], capture_output=True, text=True).stdout
    return [int(line.split()[0]) for line in out.splitlines() if len(line.split()) >= 2 and not line.split()[1].startswith("Z")]


FAKE_LIFTOFF = r"""#!/bin/bash
# stands in for Liftoff's process tree: records its ids, then holds two sleeping children (the pool and minimap2)
echo "$$ $(ps -o pgid= -p $$ | tr -d ' ') $(ps -o sid= -p $$ | tr -d ' ')" > "$PWD/pids.txt"
%s
sleep 60 &
sleep 60 &
wait
"""

# fix round 2 (S1): the same tree, but it needs no `ps` to start (second-round safety review: the clean-up and the watchdog must not depend on a forked `ps`)
FAKE_LIFTOFF_NOPS = r"""#!/bin/bash
read -r line < /proc/$$/stat
rest=${line##*) }
set -- $rest
echo "$$ $3 $4" > "$PWD/pids.txt"
%s
/usr/bin/sleep 60 &
/usr/bin/sleep 60 &
wait
"""

DRIVER_SCRIPT = r'''
import sys
sys.path.insert(0, sys.argv[1])
import liftoff_cross as LC
from pathlib import Path
e_root, stub = sys.argv[2], sys.argv[3]
hard, margin, grace = float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
d = Path(e_root) / "run"
d.mkdir(parents=True, exist_ok=True)
lock = LC.DriverLock(e_root).__enter__() if sys.argv[7] == "lock" else None
rc, peak = LC.timed_watch([stub], d / "time.txt", d / "liftoff.log", cwd=d, hard_s=hard, insession_margin_s=margin, term_grace_s=grace, poll_s=0.2,
                          pass_fds=(lock.fd,) if lock else ())
print("RC", rc, flush=True)
'''


LEDGER_DRIVER_SCRIPT = r'''
import sys
sys.path.insert(0, sys.argv[1])
import liftoff_cross as LC
e_root, stub, gff = sys.argv[2], sys.argv[3], sys.argv[4]
LC.write_run_gffs(gff, {"gene-A": "r01"}, e_root + "/gff")
a = dict(id="gene-A", component="PC", contig="NC_1", start0=100, end=1100, strand="+", idx=0, rtype="gene", biotype="x", exonic_len=1000)
plan = {"r01": dict(run=1, label="r01", records=[a], shards=[("r01", ["gene-A"])])}
LC.run_budgeted([("H", "r01")], plan, e_root=e_root, outgroups={"H": dict(species="human", genome="/t/h.fa", mmi="/t/h.mmi")}, reference="/r/g.fa",
                liftoff_bin=stub, minimap2="/m/mm2", hard_s=60, ensure_links=False, wall_pred_s={"H": 1000.0}, term_grace_s=1)
'''


class ProcessSafetyTests(HermeticCase):
    """stand-ins made of `sleep` only (no Liftoff, no minimap2): the driver may die, be signalled or be killed while a run is in flight, and the Liftoff
    session (8 to 14 GB in the real run) must not outlive it"""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.e_root = os.path.join(self.tmp, "e")
        self.sids, self.procs = [], []
        self.addCleanup(self.reap)
        self.fake = self.script("fake_liftoff.sh", FAKE_LIFTOFF % "")
        self.deaf = self.script("deaf_liftoff.sh", FAKE_LIFTOFF % "trap '' TERM")
        self.nops = self.script("nops_liftoff.sh", FAKE_LIFTOFF_NOPS % "")
        self.driver = self.script("driver.py", DRIVER_SCRIPT)

    def script(self, name, text):
        p = os.path.join(self.tmp, name)
        with open(p, "w") as fh:
            fh.write(text)
        os.chmod(p, 0o755)
        return p

    def reap(self):
        """whatever a test left behind: processes of the sessions it started whose cwd is under its temp dir, then the driver processes"""
        for sid in self.sids:
            for pid in ps_session(sid):
                try:
                    if os.readlink(f"/proc/{pid}/cwd").startswith(self.tmp):
                        os.kill(pid, signal.SIGKILL)
                except (OSError, ProcessLookupError):
                    pass
        for p in self.procs:
            if p.poll() is None:
                p.kill()
            p.wait()
            for fh in (p.stdout, p.stderr):
                if fh is not None:
                    fh.close()

    def ids(self, run_dir):
        """(pid, pgid, sid) the stand-in wrote, once its two sleeping children are up"""
        pids = os.path.join(run_dir, "pids.txt")
        t0 = time.time()
        while time.time() - t0 < 20:
            words = rd(pids).split() if os.path.exists(pids) else []
            if len(words) == 3:
                pid, pgid, sid = map(int, words)
                if len(ps_session(sid)) >= 4:
                    self.sids.append(sid)
                    return pid, pgid, sid
            time.sleep(0.05)
        self.fail("the stand-in did not come up")

    def start(self, stub, hard, margin=5.0, grace=1.0, lock=False, env=None):
        p = subprocess.Popen([sys.executable, "-B", self.driver, HERE, self.e_root, stub, str(hard), str(margin), str(grace), "lock" if lock else "nolock"],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
        self.procs.append(p)
        try:
            pid, pgid, sid = self.ids(os.path.join(self.e_root, "run"))
        except AssertionError:
            p.kill()
            self.fail(f"driver failed: {p.communicate()[1][-600:]}")
        return p, sid

    def wait_empty(self, sid, seconds):
        end = time.time() + seconds
        while ps_session(sid) and time.time() < end:
            time.sleep(0.1)
        return ps_session(sid)

    # --- SIGTERM, SIGHUP and SIGINT to the driver (the `timeout 590` of tools/rlock.sh, a hang-up, a user interrupt)
    def test_sigterm_to_the_driver_ends_the_whole_liftoff_session(self):
        p, sid = self.start(self.fake, hard=60)
        p.send_signal(signal.SIGTERM)
        p.wait(timeout=30)
        self.assertEqual(p.returncode, 128 + signal.SIGTERM)
        self.assertEqual(ps_session(sid), [])

    def test_sighup_to_the_driver_ends_the_whole_liftoff_session(self):
        p, sid = self.start(self.fake, hard=60)
        p.send_signal(signal.SIGHUP)
        p.wait(timeout=30)
        self.assertEqual(p.returncode, 128 + signal.SIGHUP)
        self.assertEqual(ps_session(sid), [])

    def test_sigint_to_the_driver_ends_the_whole_liftoff_session(self):
        p, sid = self.start(self.fake, hard=60)
        p.send_signal(signal.SIGINT)
        p.wait(timeout=30)
        self.assertEqual(p.returncode, 128 + signal.SIGINT)
        self.assertEqual(ps_session(sid), [])

    # --- the driver cannot even run its finally blocks: SIGKILL (an OOM kill of the driver)
    def test_the_session_ends_by_itself_when_the_driver_is_killed(self):
        p, sid = self.start(self.fake, hard=3, margin=1.0, grace=1)             # the in-session limit is 4 s after the start: seconds of margin for the check below
        p.kill()
        p.wait()
        self.assertTrue(ps_session(sid), "right after the kill the stand-in is still alive")
        self.assertEqual(self.wait_empty(sid, 15), [], "the in-session `timeout -s KILL` ended it")

    def test_the_leftover_processes_of_a_killed_driver_keep_the_run_lock(self):
        p, sid = self.start(self.fake, hard=3, margin=1.0, grace=1, lock=True)
        p.kill()
        p.wait()
        lock_path = os.path.join(self.e_root, "run.lock")

        def held():
            fd = os.open(lock_path, os.O_RDWR)
            try:
                fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                return False
            except BlockingIOError:
                return True
            finally:
                os.close(fd)
        self.assertTrue(held(), "the orphans inherited the lock")
        with self.assertRaises(LC.StopRun):
            LC.DriverLock(self.e_root).__enter__()
        self.assertEqual(self.wait_empty(sid, 15), [])
        self.assertFalse(held())

    def test_the_flock_timeout_chain_of_the_heavy_call_leaves_nothing_behind(self):
        """tools/rlock.sh heavy is `flock -w W LOCK timeout T cmd`: when `timeout` ends the driver the Liftoff tree must be gone, so the heavy lock is not free
        while an 8 to 14 GB process of it is still alive"""
        lock = os.path.join(self.tmp, "heavy.lock")
        p = subprocess.Popen(["flock", "-w", "5", lock, "timeout", "3", sys.executable, "-B", self.driver, HERE, self.e_root, self.fake, "60", "5", "1", "nolock"],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.procs.append(p)
        pid, pgid, sid = self.ids(os.path.join(self.e_root, "run"))
        p.wait(timeout=30)
        self.assertEqual(p.returncode, 124)
        self.assertEqual(ps_session(sid), [])
        self.assertEqual(subprocess.run(["flock", "-n", lock, "true"]).returncode, 0, "the heavy lock is free again, and only now")

    def test_a_signal_during_a_run_is_recorded_in_the_wall_ledger_and_the_session_ends(self):
        gff = self.script("ref.gff", REF_GFF)
        drv = self.script("ledger_driver.py", LEDGER_DRIVER_SCRIPT)
        p = subprocess.Popen([sys.executable, "-B", drv, HERE, self.e_root, self.fake, gff], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.procs.append(p)
        pid, pgid, sid = self.ids(os.path.join(self.e_root, "lift", "human", "runs", "r01"))
        p.send_signal(signal.SIGTERM)
        p.wait(timeout=30)
        self.assertEqual(p.returncode, 128 + signal.SIGTERM)
        self.assertEqual(ps_session(sid), [])
        entries = LC.WallLedger(self.e_root).entries("H")
        self.assertEqual([(e["label"], e["status"]) for e in entries], [("r01", "interrupted")])
        self.assertGreater(entries[0]["elapsed_s"], 0)
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "lift", "human", "runs", "r01", "done.key")))

    # --- the kill sequence itself
    def test_children_that_ignore_sigterm_are_killed_and_the_session_is_empty_on_return(self):
        d = os.path.join(self.e_root, "run")
        os.makedirs(d)
        rc, _ = LC.timed_watch([self.deaf], os.path.join(d, "time.txt"), os.path.join(d, "liftoff.log"), cwd=d, hard_s=1, poll_s=0.2, term_grace_s=1)
        self.assertEqual(rc, LC.RC_TOO_LONG)
        _, _, sid = map(int, rd(os.path.join(d, "pids.txt")).split())
        self.sids.append(sid)
        self.assertEqual(ps_session(sid), [])

    def test_a_child_left_behind_by_a_finished_run_is_killed(self):
        stub = self.script("leaves.sh", "#!/bin/bash\necho \"$$ $(ps -o pgid= -p $$ | tr -d ' ') $(ps -o sid= -p $$ | tr -d ' ')\" > \"$PWD/pids.txt\"\nsleep 60 &\nexit 0\n")
        d = os.path.join(self.e_root, "run")
        os.makedirs(d)
        rc, _ = LC.timed_watch([stub], os.path.join(d, "time.txt"), os.path.join(d, "liftoff.log"), cwd=d, hard_s=30, poll_s=0.2, term_grace_s=1)
        self.assertEqual(rc, 0)
        _, _, sid = map(int, rd(os.path.join(d, "pids.txt")).split())
        self.sids.append(sid)
        self.assertEqual(ps_session(sid), [])

    def test_the_handlers_turn_a_signal_into_system_exit_and_leave_the_old_ones_in_place(self):
        script = self.script("sig.py", "import os, signal, sys\nsys.path.insert(0, sys.argv[1])\nimport liftoff_cross as LC\nold = signal.getsignal(signal.SIGTERM)\n"
                                       "try:\n    with LC.terminate_on_signals():\n        try:\n            os.kill(os.getpid(), signal.SIGTERM)\n"
                                       "            print('NOT REACHED')\n        finally:\n            print('finally ran')\n"
                                       "except SystemExit as e:\n    print('exit', e.code)\n"
                                       "print('restored', signal.getsignal(signal.SIGTERM) is old)\n")
        r = subprocess.run([sys.executable, "-B", script, HERE], capture_output=True, text=True, timeout=30)
        self.assertEqual(r.stdout.split("\n")[:3], ["finally ran", "exit 143", "restored True"], r.stdout + r.stderr)

    # --- fix round 2, S1: the clean-up and the RSS watchdog read /proc; a `ps` that cannot be started, prints nothing or stalls changes nothing
    def bin_dir(self, name, ps_text=None):
        """a directory to be the whole PATH of a driver: empty (no `ps` at all) or with a fake `ps` whose text is given"""
        d = os.path.join(self.tmp, name)
        os.makedirs(d)
        if ps_text is not None:
            self.script(os.path.join(name, "ps"), ps_text)
        return d

    def env_with_path(self, d):
        return dict(os.environ, PATH=d)

    def test_sigterm_ends_the_whole_session_when_ps_cannot_be_started(self):
        p, sid = self.start(self.nops, hard=60, env=self.env_with_path(self.bin_dir("noprog")))
        p.send_signal(signal.SIGTERM)
        p.wait(timeout=30)
        self.assertEqual(p.returncode, 128 + signal.SIGTERM, p.stderr.read()[-400:])
        self.assertEqual(ps_session(sid), [])

    def test_sigterm_ends_the_whole_session_when_ps_prints_nothing(self):
        p, sid = self.start(self.nops, hard=60, env=self.env_with_path(self.bin_dir("silentps", "#!/bin/sh\nexit 0\n")))
        p.send_signal(signal.SIGTERM)
        p.wait(timeout=30)
        self.assertEqual(p.returncode, 128 + signal.SIGTERM)
        self.assertEqual(ps_session(sid), [])

    def test_a_ps_that_stalls_is_never_started_by_the_driver(self):
        mark = os.path.join(self.tmp, "ps_was_called")
        stall = f'#!/bin/sh\necho called >> "{mark}"\nexec /usr/bin/sleep 60\n'
        p, sid = self.start(self.nops, hard=60, env=self.env_with_path(self.bin_dir("stallps", stall)))
        p.send_signal(signal.SIGTERM)
        p.wait(timeout=60)
        self.assertEqual(p.returncode, 128 + signal.SIGTERM, p.stderr.read()[-400:])
        self.assertEqual(ps_session(sid), [])
        self.assertFalse(os.path.exists(mark), "the watchdog or the clean-up forked ps")

    def test_the_rss_watchdog_does_not_need_ps(self):
        d = os.path.join(self.e_root, "run")
        os.makedirs(d)
        with mock.patch.dict(os.environ, {"PATH": self.bin_dir("noprog2")}):
            rc, peak = LC.timed_watch([self.nops], os.path.join(d, "time.txt"), os.path.join(d, "liftoff.log"), cwd=d, hard_s=30, poll_s=0.1, term_grace_s=1,
                                      rss_limit_gb=1e-6)                         # 1 kB: any session is over it
        self.assertEqual(rc, LC.RC_TOO_BIG)
        self.assertGreater(peak, 0)
        _, _, sid = map(int, rd(os.path.join(d, "pids.txt")).split())
        self.sids.append(sid)
        self.assertEqual(ps_session(sid), [])

    def test_an_os_error_of_one_poll_is_skipped_and_the_run_goes_on(self):
        d = os.path.join(self.e_root, "run")
        os.makedirs(d)
        quick = self.script("quick.sh", "#!/bin/bash\n/usr/bin/sleep 1.0\n")
        seen = []

        def rss(pid):
            seen.append(pid)
            if len(seen) == 1:
                raise OSError(24, "Too many open files")
            return 200
        with mock.patch.object(LC, "group_rss_kb", side_effect=rss):
            rc, peak = LC.timed_watch([quick], os.path.join(d, "time.txt"), os.path.join(d, "liftoff.log"), cwd=d, hard_s=30, poll_s=0.2)
        self.assertEqual((rc, peak), (0, 200))
        self.assertGreaterEqual(len(seen), 2)

    def live_session(self, body_ids=True):
        """a stand-in tree in a session of its own, started by the test itself: (Popen, sid) once the leader and its two sleepers are up"""
        d = os.path.join(self.tmp, "live")
        os.makedirs(d, exist_ok=True)
        p = subprocess.Popen([self.nops], cwd=d, start_new_session=True, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.procs.append(p)
        self.sids.append(p.pid)
        end = time.time() + 20
        while len(ps_session(p.pid)) < 3 and time.time() < end:
            time.sleep(0.05)
        self.assertGreaterEqual(len(ps_session(p.pid)), 3)
        return p, p.pid

    def test_the_session_is_listed_from_proc_and_agrees_with_ps(self):
        p, sid = self.live_session()
        with mock.patch.dict(os.environ, {"PATH": self.bin_dir("noprog3")}):          # no ps for the code under test
            got, rss = sorted(LC.session_pids(sid)), LC.group_rss_kb(sid)
        self.assertEqual(got, sorted(ps_session(sid)))
        self.assertGreater(rss, 0)
        self.assertLess(rss, 500_000)                                                # a bash and two sleeps are a few MB, not GB
        ps_rss = sum(int(x) for x in subprocess.run(["ps", "-o", "rss=", "--sid", str(sid)], capture_output=True, text=True).stdout.split())
        self.assertLessEqual(abs(rss - ps_rss), 256, (rss, ps_rss))                  # the resident set in kB, as `ps` reports it (not the virtual size, not a number of pages)
        self.assertEqual(LC.kill_session(sid, term_grace_s=5, kill_wait_s=5), [])
        p.wait(timeout=10)

    def test_a_zombie_is_neither_listed_nor_counted(self):
        d = os.path.join(self.tmp, "zombie")
        os.makedirs(d)
        # the python parent never waits for its child: the child stays a zombie of the same session
        script = self.script("zombie_parent.py", "import subprocess, time\nsubprocess.Popen(['/usr/bin/true'])\ntime.sleep(60)\n")
        p = subprocess.Popen([sys.executable, "-B", script], cwd=d, start_new_session=True, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.procs.append(p)
        self.sids.append(p.pid)
        end = time.time() + 20
        states = []
        while time.time() < end:
            states = [st for _, st, _ in LC._session_members(p.pid)]
            if "Z" in states:
                break
            time.sleep(0.05)
        self.assertIn("Z", states, states)
        members = LC._session_members(p.pid)
        zombie = next(pid for pid, st, _ in members if st == "Z")
        with mock.patch.object(LC, "_session_members", return_value=members):             # one snapshot for both sides: the sizes of live processes move
            self.assertNotIn(zombie, LC.session_pids(p.pid))
            self.assertEqual(LC.group_rss_kb(p.pid), sum(r for _, st, r in members if st != "Z"))
        self.assertEqual(LC.kill_session(p.pid, term_grace_s=5, kill_wait_s=5), [])
        p.wait(timeout=10)

    def test_the_clean_up_sends_sigterm_to_the_group_before_it_lists_anything(self):
        p, sid = self.live_session()
        order = []
        real_killpg, real_members = os.killpg, LC._session_members
        with mock.patch.object(LC.os, "killpg", lambda pg, sig: (order.append(("killpg", sig)), real_killpg(pg, sig))[1]), \
                mock.patch.object(LC, "_session_members", lambda s: (order.append(("list",)), real_members(s))[1]):
            self.assertEqual(LC.kill_session(sid, term_grace_s=5, kill_wait_s=5), [])
        self.assertEqual(order[0], ("killpg", signal.SIGTERM))
        p.wait(timeout=10)

    def test_a_failing_listing_does_not_keep_the_group_from_getting_its_sigterm(self):
        p, sid = self.live_session()
        with mock.patch.object(LC, "_session_members", side_effect=OSError(24, "Too many open files")):
            with self.assertRaises(OSError):
                LC.kill_session(sid, term_grace_s=0.3, kill_wait_s=0.3)
        p.wait(timeout=10)                                                           # the group was signalled before the listing failed: the leader died of it
        self.assertEqual(ps_session(sid), [])

    # --- fix round 2, S2: the in-session backstop is a SIGKILL (GNU timeout's -k never fires once `time` has died of the first SIGTERM)
    def test_the_in_session_backstop_kills_a_tree_that_ignores_sigterm_when_the_driver_is_killed(self):
        p, sid = self.start(self.deaf, hard=3, margin=1.0, grace=1)                  # the backstop fires 4 s after the start
        p.kill()
        p.wait()
        self.assertTrue(ps_session(sid), "right after the kill the stand-in is still alive")
        self.assertEqual(self.wait_empty(sid, 20), [], "`timeout -s KILL` ended a tree that ignores SIGTERM")

    def test_a_backstop_that_fires_before_the_drivers_own_limit_is_reported_as_too_long_and_leaves_nothing_running(self):
        # `timeout -s KILL` kills its whole group, itself included: the leader then looks SIGKILLed (rc -9) instead of timed out (124). When the driver was too slow to
        # notice its own limit, the run still exceeded it: the shard is classified as too long (and split), not as a Liftoff failure.
        d = os.path.join(self.e_root, "run")
        os.makedirs(d)
        rc, _ = LC.timed_watch([self.nops], os.path.join(d, "time.txt"), os.path.join(d, "liftoff.log"), cwd=d, hard_s=30, insession_margin_s=-29.0, poll_s=0.5, term_grace_s=1)
        self.assertEqual(rc, LC.RC_TOO_LONG)                                       # the backstop fired at 1 s, the driver's limit is 30 s
        _, _, sid = map(int, rd(os.path.join(d, "pids.txt")).split())
        self.sids.append(sid)
        self.assertEqual(ps_session(sid), [])

    def test_the_tree_runs_under_timeout_dash_s_kill(self):
        seen = []

        class Rec(subprocess.Popen):
            def __init__(self, args, *a, **k):
                seen.append(list(args))
                super().__init__(args, *a, **k)
        d = os.path.join(self.e_root, "run")
        os.makedirs(d)
        quick = self.script("quick2.sh", "#!/bin/bash\n/usr/bin/sleep 0.2\n")
        with mock.patch.object(LC.subprocess, "Popen", Rec):
            rc, _ = LC.timed_watch([quick], os.path.join(d, "time.txt"), os.path.join(d, "liftoff.log"), cwd=d, hard_s=20, insession_margin_s=3, term_grace_s=1, poll_s=0.1)
        self.assertEqual(rc, 0)
        argv = next(a for a in seen if a and a[0] == LC.TIMEOUT_BIN)
        self.assertEqual(argv[:4], [LC.TIMEOUT_BIN, "-s", "KILL", "23.000"])
        self.assertNotIn("-k", argv[:6])
        self.assertEqual(argv[4], LC.TIME_BIN)

    # --- fix round 2, S3: a signal that lands while the child is being created is remembered and acted on as soon as the session is known
    def test_a_signal_that_lands_while_the_child_is_being_created_still_ends_the_session(self):
        created = []

        class SigPopen(subprocess.Popen):
            def _execute_child(self, *a, **k):
                super()._execute_child(*a, **k)
                created.append(self.pid)
                os.kill(os.getpid(), signal.SIGTERM)                                # the child exists; Popen has not returned: `p` is not assigned yet
        d = os.path.join(self.e_root, "run")
        os.makedirs(d)
        with mock.patch.object(LC.subprocess, "Popen", SigPopen):
            with self.assertRaises(SystemExit) as cm:
                LC.timed_watch([self.nops], os.path.join(d, "time.txt"), os.path.join(d, "liftoff.log"), cwd=d, hard_s=60, poll_s=0.1, term_grace_s=1)
        self.assertEqual(cm.exception.code, 128 + signal.SIGTERM)
        self.assertEqual(len(created), 1)
        self.sids.append(created[0])
        self.assertEqual(ps_session(created[0]), [])

    def test_a_launch_that_fails_leaves_the_deferral_off(self):
        d = os.path.join(self.e_root, "run")
        os.makedirs(d)
        with mock.patch.object(LC.subprocess, "Popen", side_effect=OSError(11, "Resource temporarily unavailable")):
            with self.assertRaises(OSError):
                LC.timed_watch([self.nops], os.path.join(d, "time.txt"), os.path.join(d, "liftoff.log"), cwd=d, hard_s=30)
        self.assertFalse(LC._DEFER["on"])                                          # a later signal is raised again, not swallowed

    def test_a_deferred_signal_does_not_leak_into_the_next_run(self):
        d = os.path.join(self.e_root, "run")
        os.makedirs(d)
        quick = self.script("quick3.sh", "#!/bin/bash\n/usr/bin/sleep 0.2\n")
        LC._DEFER.update(on=False, sig=signal.SIGHUP)                                # a stale flag from a run that ended some other way
        rc, _ = LC.timed_watch([quick], os.path.join(d, "time.txt"), os.path.join(d, "liftoff.log"), cwd=d, hard_s=20, poll_s=0.1, term_grace_s=1)
        self.assertEqual(rc, 0)
        self.assertFalse(LC._DEFER["on"])

    def test_the_child_does_not_inherit_a_blocked_signal_mask(self):
        # the reason the deferral is a flag and not signals_held() around Popen: a blocked mask is inherited by the child (measured: SIGTERM then took 10 s)
        d = os.path.join(self.e_root, "run")
        os.makedirs(d)
        probe = self.script("mask.sh", "#!/bin/bash\ngrep '^SigBlk' /proc/$$/status > mask.txt\n")
        rc, _ = LC.timed_watch([probe], os.path.join(d, "time.txt"), os.path.join(d, "liftoff.log"), cwd=d, hard_s=20, poll_s=0.1, term_grace_s=1)
        self.assertEqual(rc, 0)
        blocked = int(rd(os.path.join(d, "mask.txt")).split()[1], 16)
        for sig in (signal.SIGTERM, signal.SIGHUP, signal.SIGINT):
            self.assertFalse(blocked & (1 << (sig - 1)), sig)


class StrayProcessTests(HermeticCase):
    """the pre-launch guard: no liftoff / minimap2 of this user may be running in the run root"""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.e_root = os.path.join(self.tmp, "e")
        self.run_dir = os.path.join(self.e_root, "lift", "human", "runs", "pos1")
        os.makedirs(self.run_dir)
        os.makedirs(os.path.join(self.tmp, "elsewhere"))
        self.procs = []
        self.addCleanup(self.reap)

    def reap(self):
        for p in self.procs:
            try:
                os.killpg(p.pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass
            p.wait()

    def spawn(self, name, cwd, *args):
        """a process whose argv names `name` (a script that sleeps) started in `cwd` with `args`"""
        script = os.path.join(self.tmp, "bin", name)
        os.makedirs(os.path.dirname(script), exist_ok=True)
        with open(script, "w") as fh:
            fh.write("#!/bin/bash\nsleep 60 &\nwait\n")
        os.chmod(script, 0o755)
        p = subprocess.Popen([script, *args], cwd=cwd, start_new_session=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.procs.append(p)
        t0 = time.time()
        while time.time() - t0 < 10:                                           # wait until bash has become the script (cmdline readable)
            try:
                if name in rd(f"/proc/{p.pid}/cmdline"):
                    break
            except OSError:
                pass
            time.sleep(0.02)
        return p

    def test_a_liftoff_process_whose_cwd_is_in_the_root_is_found(self):
        p = self.spawn("liftoff", self.run_dir)
        hits = LC.find_stray_processes(self.e_root)
        self.assertIn(p.pid, [pid for pid, _ in hits])

    def test_a_minimap2_process_with_the_root_in_its_command_line_is_found(self):
        p = self.spawn("minimap2", os.path.join(self.tmp, "elsewhere"), os.path.join(self.e_root, "lift", "human", "genome.fa.mmi"))
        self.assertIn(p.pid, [pid for pid, _ in LC.find_stray_processes(self.e_root)])

    def test_the_cache_wrapper_is_found_too(self):
        p = self.spawn("mm2_cache.sh", self.run_dir)
        self.assertIn(p.pid, [pid for pid, _ in LC.find_stray_processes(self.e_root)])

    def test_a_liftoff_process_of_another_root_is_not_ours(self):
        self.spawn("liftoff", os.path.join(self.tmp, "elsewhere"), "/somewhere/else/genome.fa")
        self.assertEqual(LC.find_stray_processes(self.e_root), [])

    def test_other_processes_in_the_root_are_not_strays(self):
        self.spawn("python-driver", self.run_dir, os.path.join(self.e_root, "x"))
        self.assertEqual(LC.find_stray_processes(self.e_root), [])

    @staticmethod
    def eventually(fn, cap=30.0, step=0.02):
        """poll `fn` until it returns something true or `cap` seconds have passed (a generous cap, never an upper bound on speed); the assertion that follows is exact"""
        end = time.time() + cap
        while True:
            v = fn()
            if v or time.time() >= end:
                return v
            time.sleep(step)

    def test_the_guard_refuses_and_names_the_process(self):
        # fix round 2 (F8 of the mutation review: 21 failures in 4,000 runs under CPU load): the bash stand-in forks before it execs, and a killed process is still listed (its
        # command line and working directory readable) for a moment after SIGKILL, so the guard was asked once too early. Both ends are polled now.
        p = self.spawn("liftoff", self.run_dir)
        self.assertTrue(self.eventually(lambda: p.pid in [pid for pid, _ in LC.find_stray_processes(self.e_root)]), "the stand-in never showed up as a stray")
        with self.assertRaises(LC.StopRun) as cm:
            LC.guard_no_stray(self.e_root)
        self.assertIn(str(p.pid), str(cm.exception))
        os.killpg(p.pid, signal.SIGKILL)
        p.wait()

        def no_refusal():
            try:
                LC.guard_no_stray(self.e_root)
            except LC.StopRun:
                return False
            return True
        self.assertTrue(self.eventually(no_refusal), "the killed stand-in (and the process it forked) is still listed as a stray after the cap")           # gone: no refusal


class DriverLockTests(HermeticCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.e_root = os.path.join(self.tmp, "e")

    def test_a_second_driver_is_refused_while_the_first_holds_the_lock(self):
        with LC.DriverLock(self.e_root):
            with self.assertRaises(LC.StopRun) as cm:
                LC.DriverLock(self.e_root).__enter__()
            self.assertIn("run.lock", str(cm.exception))
        with LC.DriverLock(self.e_root):                                        # released with the first driver
            pass

    def test_the_lock_file_names_its_holder(self):
        with LC.DriverLock(self.e_root):
            self.assertEqual(rd(os.path.join(self.e_root, "run.lock")).split()[0], str(os.getpid()))


class StrayGuardInTheDriverTests(HermeticCase):
    """run_budgeted refuses to start a run while a liftoff of this root is alive"""

    def test_a_stray_liftoff_stops_the_call_before_any_run_directory_is_made(self):
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, True)
        e_root = os.path.join(tmp, "e")
        gff = os.path.join(tmp, "ref.gff")
        with open(gff, "w") as fh:
            fh.write(REF_GFF)
        LC.write_run_gffs(gff, {"gene-A": "r01"}, os.path.join(e_root, "gff"))
        run_dir = os.path.join(e_root, "lift", "human", "runs", "elsewhere")
        os.makedirs(run_dir)
        script = os.path.join(tmp, "liftoff")
        with open(script, "w") as fh:
            fh.write("#!/bin/bash\nsleep 60 &\nwait\n")
        os.chmod(script, 0o755)
        stray = subprocess.Popen([script], cwd=run_dir, start_new_session=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.addCleanup(lambda: (os.killpg(stray.pid, signal.SIGKILL), stray.wait()))
        time.sleep(0.2)
        a = dict(id="gene-A", component="PC", contig="NC_1", start0=100, end=1100, strand="+", idx=0, rtype="gene", biotype="x", exonic_len=1000)
        plan = {"r01": dict(run=1, label="r01", records=[a], shards=[("r01", ["gene-A"])])}
        with self.assertRaises(LC.StopRun) as cm:
            LC.run_budgeted([("H", "r01")], plan, e_root=e_root, outgroups={"H": dict(species="human", genome="/t/h.fa", mmi="/t/h.mmi")}, reference="/r/g.fa",
                            liftoff_bin="/bin/false", minimap2="/m/mm2", hard_s=60, ensure_links=False)
        self.assertIn(str(stray.pid), str(cm.exception))
        self.assertFalse(os.path.exists(os.path.join(e_root, "lift", "human", "runs", "r01")))

    def test_the_guard_runs_before_the_freeze_so_a_refused_call_leaves_nothing_frozen(self):
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, True)
        e_root = os.path.join(tmp, "e")
        gff = os.path.join(tmp, "ref.gff")
        with open(gff, "w") as fh:
            fh.write(REF_GFF)
        LC.write_run_gffs(gff, {"gene-A": "r01"}, os.path.join(e_root, "gff"))
        run_dir = os.path.join(e_root, "lift", "human", "runs", "elsewhere")
        os.makedirs(run_dir)
        script = os.path.join(tmp, "liftoff")
        with open(script, "w") as fh:
            fh.write("#!/bin/bash\nsleep 60 &\nwait\n")
        os.chmod(script, 0o755)
        stray = subprocess.Popen([script], cwd=run_dir, start_new_session=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.addCleanup(lambda: (os.killpg(stray.pid, signal.SIGKILL), stray.wait()))
        time.sleep(0.2)
        a = dict(id="gene-A", component="PC", contig="NC_1", start0=100, end=1100, strand="+", idx=0, rtype="gene", biotype="x", exonic_len=1000)
        plan = {"r01": dict(run=1, label="r01", records=[a], shards=[("r01", ["gene-A"])])}
        froze = []
        with self.assertRaises(LC.StopRun):
            LC.run_budgeted([("H", "r01")], plan, e_root=e_root, outgroups={"H": dict(species="human", genome="/t/h.fa", mmi="/t/h.mmi")}, reference="/r/g.fa",
                            liftoff_bin="/bin/false", minimap2="/m/mm2", hard_s=60, ensure_links=False, on_first_launch=lambda: froze.append(1))
        self.assertEqual(froze, [])


class WallLedgerTests(HermeticCase):
    """stop rule 6 needs the Liftoff wall of every run of an outgroup, over all the calls of the phase"""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def test_wall_is_summed_per_outgroup_and_read_back_by_the_next_call(self):
        led = LC.WallLedger(self.tmp)
        led.add("H", "r01", 10.0, "done")
        led.add("H", "r02", 5.5, "toolong")
        led.add("C", "r01", 1.0, "done")
        again = LC.WallLedger(self.tmp)
        self.assertAlmostEqual(again.total("H"), 15.5)
        self.assertAlmostEqual(again.total("C"), 1.0)
        self.assertEqual(again.total("O"), 0.0)
        self.assertEqual([e["label"] for e in again.entries("H")], ["r01", "r02"])

    def test_the_file_is_written_whole(self):
        LC.WallLedger(self.tmp).add("H", "r01", 1.0, "done")
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "wall.json.tmp")))
        self.assertEqual(json.loads(rd(os.path.join(self.tmp, "wall.json")))["H"][0]["label"], "r01")

    def test_a_corrupt_ledger_is_a_stop_not_a_zero(self):
        with open(os.path.join(self.tmp, "wall.json"), "w") as fh:
            fh.write("{\"H\": [")
        with self.assertRaises(LC.StopRun) as cm:
            LC.WallLedger(self.tmp)
        self.assertIn("wall.json", str(cm.exception))


class InputPinTests(HermeticCase):
    """run_identity.json pins every input of Amendment 2 section 3: sha256 below 1 GB, size and mtime above"""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.small, self.big = os.path.join(self.tmp, "small.tsv"), os.path.join(self.tmp, "big.fa")
        with open(self.small, "w") as fh:
            fh.write("abc\n")
        with open(self.big, "w") as fh:
            fh.write("x" * 30)

    def test_small_files_carry_a_sha256_big_ones_only_size_and_mtime(self):
        import hashlib
        with mock.patch.object(LC, "INPUT_BIG_BYTES", 20):
            pins = LC.pin_inputs({"small": self.small, "big": self.big}, deep=True)
        self.assertEqual(pins["small"]["sha256"], hashlib.sha256(b"abc\n").hexdigest())
        self.assertIsNone(pins["big"]["sha256"])
        for name, path in (("small", self.small), ("big", self.big)):
            self.assertEqual((pins[name]["path"], pins[name]["size"], pins[name]["mtime_ns"]), (path, os.path.getsize(path), os.stat(path).st_mtime_ns))

    def test_a_shallow_pin_skips_the_hash_of_a_file_above_the_shallow_limit(self):
        with mock.patch.object(LC, "SHALLOW_SHA_BYTES", 2):
            self.assertIsNone(LC.pin_inputs({"small": self.small}, deep=False)["small"]["sha256"])
            self.assertEqual(len(LC.pin_inputs({"small": self.small}, deep=True)["small"]["sha256"]), 64)

    def test_nothing_changed_nothing_is_reported(self):
        a = LC.pin_inputs({"small": self.small, "big": self.big}, deep=True)
        self.assertEqual(LC.compare_pins(a, LC.pin_inputs({"small": self.small, "big": self.big}, deep=False)), [])

    def test_a_changed_content_size_or_mtime_is_reported_by_input_name(self):
        a = LC.pin_inputs({"small": self.small, "big": self.big}, deep=True)
        with open(self.small, "w") as fh:                                   # same size, other content
            fh.write("abd\n")
        st = os.stat(self.small)
        os.utime(self.small, ns=(st.st_atime_ns, a["small"]["mtime_ns"]))   # even with the old mtime put back
        problems = LC.compare_pins(a, LC.pin_inputs({"small": self.small, "big": self.big}, deep=True))
        self.assertEqual(len(problems), 1)
        self.assertIn("small", problems[0])
        self.assertIn("changed", problems[0])
        with open(self.big, "a") as fh:
            fh.write("y")
        problems = LC.compare_pins(a, LC.pin_inputs({"small": self.small, "big": self.big}, deep=False))
        self.assertTrue(any("big" in x and "size" in x for x in problems), problems)
        os.utime(self.big, ns=(1, 1))
        self.assertTrue(any("big" in x and "mtime" in x for x in LC.compare_pins(a, LC.pin_inputs({"big": self.big}, deep=False))))

    def test_an_input_that_is_not_pinned_or_not_there_is_reported(self):
        a = LC.pin_inputs({"small": self.small}, deep=True)
        b = LC.pin_inputs({"small": self.small, "big": self.big}, deep=True)
        self.assertTrue(any("big" in x for x in LC.compare_pins(a, b)))
        self.assertTrue(any("big" in x for x in LC.compare_pins(b, a)))

    def test_a_missing_input_cannot_be_pinned(self):
        with self.assertRaises(LC.StopRun) as cm:
            LC.pin_inputs({"gone": os.path.join(self.tmp, "gone.tsv")}, deep=True)
        self.assertIn("gone.tsv", str(cm.exception))

    def test_the_registered_configuration_has_the_seventeen_inputs_of_section_3(self):
        paths = LC.input_paths(gorilla_fasta="/g.fa", gorilla_gff="/g.gff", lift_root="/lift", outgroups=LC.OUTGROUPS)
        want = {"gorilla_fasta", "gorilla_fai", "gorilla_gff", "gorilla_selflift", "gorilla_records"}
        want |= {f"{k}_{x}" for k in "HCO" for x in ("fasta", "fai", "mmi", "selflift")}
        self.assertEqual(set(paths), want)
        self.assertEqual(len(paths), 17)
        self.assertEqual(paths["gorilla_fai"], "/g.fa.fai")
        self.assertEqual(paths["gorilla_selflift"], "/lift/gorilla/liftoff_loci.tsv")
        self.assertTrue(paths["H_fasta"].endswith("chm13v2.0.fa") and paths["H_fai"].endswith("chm13v2.0.fa.fai") and paths["H_mmi"].endswith("chm13v2.0.fa.mmi"))


class OutgroupPreflightTests(HermeticCase):
    """the three links Liftoff is given must point at files that exist and that no tool will rewrite"""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.fa = os.path.join(self.tmp, "g", "h.fa")
        make_fake_genome(self.fa)
        self.og = {"H": dict(species="human", genome=self.fa, mmi=self.fa + ".mmi")}
        self.root = os.path.join(self.tmp, "e")

    def link(self, name):
        return os.path.join(self.root, "lift", "human", name)

    def test_links_are_made_to_files_that_exist_and_made_again_without_complaint(self):
        d = LC.ensure_target_links("H", self.root, self.og)
        for name, src in (("genome.fa", self.fa), ("genome.fa.fai", self.fa + ".fai"), ("genome.fa.mmi", self.fa + ".mmi")):
            self.assertEqual(os.path.realpath(d / name), os.path.realpath(src))
        LC.ensure_target_links("H", self.root, self.og)

    def test_a_missing_source_is_a_stop_and_no_link_is_made(self):
        os.remove(self.fa + ".mmi")
        with self.assertRaises(LC.StopRun) as cm:
            LC.ensure_target_links("H", self.root, self.og)
        self.assertIn(self.fa + ".mmi", str(cm.exception))
        self.assertFalse(os.path.lexists(self.link("genome.fa")))

    def test_an_empty_source_is_a_stop(self):
        open(self.fa + ".fai", "w").close()
        with self.assertRaises(LC.StopRun):
            LC.ensure_target_links("H", self.root, self.og)

    def test_a_fai_older_than_its_fasta_is_a_stop_because_pyfaidx_would_rewrite_it_through_the_link(self):
        os.utime(self.fa + ".fai", (1, 1))
        with self.assertRaises(LC.StopRun) as cm:
            LC.ensure_target_links("H", self.root, self.og)
        self.assertIn("older", str(cm.exception))

    def test_an_index_of_another_k_or_w_or_sequence_count_is_a_stop(self):
        import struct
        for w, k, n in ((11, 15, 2), (10, 21, 2), (10, 15, 1)):
            with open(self.fa + ".mmi", "wb") as fh:
                fh.write(b"MMI\x02" + struct.pack("<5I", w, k, 14, n, 0) + b"\0" * 16)
            with self.assertRaises(LC.StopRun) as cm:
                LC.ensure_target_links("H", self.root, self.og)
            self.assertIn("index", str(cm.exception))

    def test_an_existing_link_to_another_file_or_to_nothing_is_a_stop(self):
        os.makedirs(os.path.dirname(self.link("genome.fa")))
        other = os.path.join(self.tmp, "other.fa")
        open(other, "w").write(">x\nA\n")
        os.symlink(other, self.link("genome.fa"))
        with self.assertRaises(LC.StopRun) as cm:
            LC.ensure_target_links("H", self.root, self.og)
        self.assertIn("not a link to", str(cm.exception))
        os.remove(self.link("genome.fa"))
        os.symlink(os.path.join(self.tmp, "gone.fa"), self.link("genome.fa"))                         # dangling: Liftoff would build an index through it
        with self.assertRaises(LC.StopRun):
            LC.ensure_target_links("H", self.root, self.og)

    def test_the_reference_needs_its_fai_and_the_fai_must_not_be_older(self):
        ref = os.path.join(self.tmp, "ref", "g.fa")
        make_fake_genome(ref, mmi=False)
        LC.check_reference(ref)
        os.utime(ref + ".fai", (1, 1))
        with self.assertRaises(LC.StopRun):
            LC.check_reference(ref)
        os.remove(ref + ".fai")
        with self.assertRaises(LC.StopRun):
            LC.check_reference(ref)


class RehearsalGuardTests(CliCase):
    """a stub rehearsal must never leave output that the real run (or `build` of the real root) takes for its own"""

    def setUp(self):
        super().setUp()
        self.real = os.path.join(self.tmp, "real_e")
        for pt in (mock.patch.object(LC, "E_ROOT", self.real), mock.patch.object(X, "LIFT_ROOT", self.root)):
            pt.start()
            self.addCleanup(pt.stop)

    def raw(self, *argv):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = LC.main(list(argv))
        return rc, out.getvalue(), err.getvalue()

    def test_flags_that_change_what_a_run_is_made_of_are_refused_on_the_real_root(self):
        cases = {"--liftoff-bin": [self.stub], "--minimap2": ["/m/mm2"], "--no-code-check": [], "--lift-root": ["/elsewhere/lift"], "--gorilla-gff": [self.gff],
                 "--gorilla-fasta": [self.fasta], "--n-planted": ["4"], "--strict-single-copy": []}
        for flag, val in cases.items():
            rc, out, err = self.raw("run", "--out", self.real, "--amendment-commit", "abc1234", flag, *val)
            self.assertEqual(rc, 2, (flag, err))
            self.assertIn(flag, err)
            self.assertIn("rehearsal", err)
            self.assertFalse(os.path.exists(self.real), f"{flag} created something in the real root")

    # --- fix round 2, G6: the other hub reading is a rehearsal flag too (a typo would freeze another plan on the real root)
    def test_a_hub_mode_other_than_the_registered_one_is_refused_on_the_real_root(self):
        for argv in (("run", "--out", self.real, "--amendment-commit", "abc1234", "--hub-mode", "records"), ("gffs", "--out", self.real, "--hub-mode", "records"),
                     ("build", "--out", self.real, "--hub-mode", "records")):
            rc, out, err = self.raw(*argv)
            self.assertEqual(rc, 2, (argv, err))
            self.assertIn("--hub-mode", err)
            self.assertFalse(os.path.exists(self.real), argv)

    def test_the_registered_hub_mode_is_not_a_rehearsal_flag_and_plan_is_never_refused(self):
        a = types.SimpleNamespace(cmd="run", liftoff_bin=LC.LIFTOFF_BIN, minimap2=None, no_code_check=False, lift_root=X.LIFT_ROOT, gorilla_gff=LC.GORILLA_GFF,
                                  gorilla_fasta=LC.GORILLA_FASTA, n_planted=LC.N_PLANTED, strict_single_copy=False, rehearsal=False, hub_mode="nodes")
        self.assertEqual(LC.rehearsal_flags(a), [])
        a.hub_mode = "records"
        self.assertEqual(LC.rehearsal_flags(a), ["--hub-mode"])
        rc, out, err = self.raw("plan", "--out", self.real, "--hub-mode", "records", "--lift-root", self.root, "--n-planted", "4")
        self.assertEqual(rc, 0, err)

    def test_a_rehearsal_with_the_other_hub_mode_leaves_a_marker_naming_it(self):
        self.write_map()
        rc, out, err = self.cli("run", "--hub-mode", "records")
        self.assertIn(rc, (75, 0), err)
        self.assertIn("--hub-mode", json.loads(rd(os.path.join(self.e_root, "REHEARSAL.json")))["flags"])

    def test_the_default_out_is_the_real_root_and_is_refused_too(self):
        rc, out, err = self.raw("run", "--amendment-commit", "abc1234", "--liftoff-bin", self.stub)             # no --out: E_ROOT
        self.assertEqual(rc, 2, err)
        self.assertFalse(os.path.exists(self.real))

    def test_gffs_and_build_flags_are_refused_on_the_real_root_too(self):
        for argv in (("gffs", "--out", self.real, "--gorilla-gff", self.gff), ("build", "--out", self.real, "--rehearsal"),
                     ("build", "--out", self.real, "--lift-root", "/elsewhere/lift"), ("build", "--out", self.real, "--n-planted", "4")):
            rc, out, err = self.raw(*argv)
            self.assertEqual(rc, 2, (argv, err))
            self.assertFalse(os.path.exists(self.real), argv)

    def test_a_directory_below_the_real_root_is_a_rehearsal_directory(self):
        below = os.path.join(self.real, "rehearsal9")                       # like the existing E_ROOT/rehearsal2
        self.write_map()
        rc, out, err = self.cli("run", out=below)
        self.assertEqual(rc, 75, err)
        self.assertTrue(os.path.exists(os.path.join(below, "REHEARSAL.json")))

    def test_a_stub_run_leaves_a_marker_and_build_refuses_it_unless_it_is_a_rehearsal(self):
        self.write_map()
        seen, _ = self.run_until_done()
        self.assertEqual(seen[-1][0], 0, seen)
        marker = os.path.join(self.e_root, "REHEARSAL.json")
        self.assertTrue(os.path.exists(marker))
        self.assertIn("--liftoff-bin", json.loads(rd(marker))["flags"])
        rc, out, err = self.cli("build", rehearsal=False)
        self.assertEqual(rc, 2, err)
        self.assertIn("rehearsal", err.lower())
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "table")))
        self.assertEqual(self.cli("build")[0], 0)

    def test_a_second_driver_is_refused_while_one_holds_the_run_lock(self):
        self.write_map()
        with LC.DriverLock(self.e_root):
            rc, out, err = self.cli("run")
        self.assertEqual(rc, 3, err)
        self.assertIn("run.lock", err)
        self.assertEqual(self.run_dirs(), [])


class PinnedInputsInTheRunTests(CliCase):
    def first_call(self):
        self.write_map()
        rc, out, err = self.cli("run")
        self.assertEqual(rc, 75, err)
        return json.loads(rd(os.path.join(self.e_root, "run_identity.json")))

    def test_the_first_call_records_the_inputs_the_tools_and_the_command(self):
        ident = self.first_call()
        want = {"gorilla_fasta", "gorilla_fai", "gorilla_gff", "gorilla_selflift", "gorilla_records"} | {f"{k}_{x}" for k in "HCO" for x in ("fasta", "fai", "mmi", "selflift")}
        self.assertEqual(set(ident["inputs"]), want)
        self.assertEqual(ident["inputs"]["H_selflift"]["path"], self.og["H"]["selflift"])
        self.assertEqual(len(ident["inputs"]["H_selflift"]["sha256"]), 64)
        self.assertEqual(sorted(ident["tools"]), ["gnu_time", "liftoff", "minimap2"])
        self.assertEqual(ident["command"]["liftoff_bin"], self.stub)
        self.assertEqual(ident["command"]["reference"], self.fasta)
        self.assertIn("figures/figlib.py", ident["code_sha256"])
        self.assertEqual(sorted(ident["code_sha256"]), sorted(LC.CODE_FILES))              # fix round 2: all 8 files of the frozen set, the second mutation review's tests included
        self.assertIn("bench/hierarchy/test_supplement_e2.py", ident["code_sha256"])

    def test_a_changed_input_is_a_stop_in_a_later_call_before_anything_starts(self):
        self.first_call()
        before = self.run_dirs()
        with open(self.og["H"]["selflift"], "a") as fh:
            fh.write("\n")                                                  # the self-lift table of the human outgroup was regenerated
        rc, out, err = self.cli("run")
        self.assertEqual(rc, 3, err)
        self.assertIn("H_selflift", err)
        self.assertIn("changed", err)
        self.assertEqual(self.run_dirs(), before)

    def test_a_touched_input_is_found_through_its_mtime(self):
        self.first_call()
        st = os.stat(self.og["C"]["mmi"])
        os.utime(self.og["C"]["mmi"], ns=(st.st_atime_ns, st.st_mtime_ns + 10 ** 9))
        rc, out, err = self.cli("run")
        self.assertEqual(rc, 3, err)
        self.assertIn("C_mmi", err)

    def test_other_tools_than_the_first_call_used_are_a_stop(self):
        self.first_call()
        other = os.path.join(self.tmp, "stub2.py")
        shutil.copy(self.stub, other)
        with open(other, "a") as fh:
            fh.write("\n# another liftoff\n")
        os.chmod(other, 0o755)
        rc, out, err = self.cli("run", liftoff_bin=other)
        self.assertEqual(rc, 3, err)
        self.assertIn("identity", err)


class RealModeFlowTests(CliCase):
    """the registered configuration with the registered paths replaced by temp files: no stub flag is on the command line, so the rules of the real root apply
    (the run identity, the keys of the runs and the pins are checked by `build`)"""

    def setUp(self):
        super().setUp()
        self.mm2 = os.path.join(self.tmp, "minimap2_stand_in")
        with open(self.mm2, "wb") as fh:
            fh.write(b"not a minimap2")
        # fix round 2 (G4, G6): the code check records what it is asked; "git rev-parse" is replaced by a resolver of two made-up commits
        self.checked = []
        self.full = "abc1234" + "0" * 33
        self.other = "def5678" + "0" * 33

        def fake_check(commit, *a, **k):
            self.checked.append((commit, k.get("needles")))
            return "HEAD"

        def fake_resolve(arg, repo=None):
            for full in (self.full, self.other):
                if full.startswith(arg):
                    return full
            raise LC.StopRun(f"{arg} is not a commit")
        for pt in (mock.patch.object(LC, "E_ROOT", self.e_root), mock.patch.object(LC, "LIFTOFF_BIN", self.stub), mock.patch.object(LC, "GORILLA_GFF", self.gff),
                   mock.patch.object(LC, "GORILLA_FASTA", self.fasta), mock.patch.object(LC, "N_PLANTED", 4), mock.patch.object(LC, "MINIMAP2_ENV", self.mm2),
                   mock.patch.object(X, "LIFT_ROOT", self.root), mock.patch.object(LC, "check_committed_code", fake_check), mock.patch.object(LC, "resolve_commit", fake_resolve)):
            pt.start()
            self.addCleanup(pt.stop)
        self.write_map()

    def real(self, cmd, *extra, commit="abc1234"):
        argv = [cmd, *extra] + (["--amendment-commit", commit] if cmd == "run" else [])
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = LC.main(argv)
        return rc, out.getvalue(), err.getvalue()

    def run_all(self):
        seen, err = [], ""
        for _ in range(8):
            rc, out, err = self.real("run")
            seen.append((rc, json.loads(out)["group"] if out.strip() else None))
            if rc != 75:
                break
        return seen, err

    # --- fix round 2, G4: the commit is resolved to 40 hex before it is checked, stored and compared
    def test_an_abbreviated_amendment_commit_is_stored_as_40_hex_and_any_spelling_of_it_is_accepted_later(self):
        rc, out, err = self.real("run", commit="abc1234")
        self.assertEqual(rc, 75, err)
        self.assertEqual(json.loads(rd(os.path.join(self.e_root, "run_identity.json")))["amendment_commit"], self.full)
        self.assertEqual(self.checked[0][0], self.full)
        for spelling in (self.full, "abc12340", "abc1234"):
            rc, out, err = self.real("run", commit=spelling)
            self.assertIn(rc, (75, 0), (spelling, err))
        self.assertEqual({c for c, _ in self.checked}, {self.full})

    def test_another_commit_than_the_first_call_recorded_is_still_a_stop(self):
        self.assertEqual(self.real("run", commit="abc1234")[0], 75)
        rc, out, err = self.real("run", commit="def5678")
        self.assertEqual(rc, 3, err)
        self.assertIn("Amendment commit", err)

    def test_a_commit_that_does_not_resolve_is_a_stop_before_anything_is_recorded(self):
        rc, out, err = self.real("run", commit="1111111")
        self.assertEqual(rc, 3, err)
        self.assertIn("is not a commit", err)
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "run_identity.json")))
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "plan.digest")))
        self.assertEqual(self.checked, [])

    # --- fix round 2, G6: on the first call the committed Amendment must also list the plan digest and the sha256 of every pinned input
    def test_the_first_call_asks_the_amendment_for_the_plan_digest_and_the_pinned_inputs_and_a_later_call_does_not(self):
        self.assertEqual(self.real("run")[0], 75)
        ident = json.loads(rd(os.path.join(self.e_root, "run_identity.json")))
        w = LC.plan_from_real_inputs("nodes", self.root, n_neg=4)
        digest = LC.plan_digest(w["plan"] + w["pos"] + w["ctl"], "nodes")
        base = {**ident["code_sha256"], **{f"tool {k}": v for k, v in ident["tools"].items()}}
        want = {**base, "plan digest": digest, **{f"input {n}": p["sha256"] for n, p in ident["inputs"].items() if p["sha256"]}}
        self.assertEqual(self.checked[0][1], want)
        self.assertEqual(len([k for k in want if k.startswith("input ")]), 17)                      # the fixture's inputs are all small; the real ones: 10 of 17
        self.real("run")
        self.assertEqual(self.checked[1][1], base)                                                  # identity exists: the plan and the inputs are compared with what it froze

    def test_only_the_inputs_that_carry_a_sha256_are_needles(self):
        paths = LC.input_paths(self.fasta, self.gff, self.root, LC.OUTGROUPS)
        sizes = sorted(os.path.getsize(p) for p in paths.values())
        limit = (sizes[0] + sizes[-1]) // 2                                       # some inputs below it (hashed), some above (size and mtime only: the 7 big files of the real run)
        self.assertTrue(sizes[0] < limit < sizes[-1], sizes)
        with mock.patch.object(LC, "INPUT_BIG_BYTES", limit):
            self.assertEqual(self.real("run")[0], 75)
        pins = json.loads(rd(os.path.join(self.e_root, "run_identity.json")))["inputs"]
        hashed = {n for n, p in pins.items() if p["sha256"]}
        self.assertTrue(hashed and hashed != set(pins), "the fixture must have both kinds")
        self.assertEqual({k[len("input "):] for k in self.checked[0][1] if k.startswith("input ")}, hashed)

    def test_a_refusal_of_the_code_check_on_the_first_call_leaves_nothing_frozen_with_the_inputs_already_hashed(self):
        def refuse(commit, *a, **k):
            raise LC.StopRun("the Amendment does not list the sha256 of input gorilla_fai")
        with mock.patch.object(LC, "check_committed_code", refuse):
            rc, out, err = self.real("run")
        self.assertEqual(rc, 3, err)
        self.assertIn("input gorilla_fai", err)
        for name in ("run_identity.json", "plan.digest", "lift"):
            self.assertFalse(os.path.exists(os.path.join(self.e_root, name)), name)

    def test_the_whole_flow_without_a_single_stub_flag(self):
        seen, err = self.run_all()
        self.assertEqual(seen, [(75, "pos"), (75, "neg"), (0, "verdict")], err)
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "REHEARSAL.json")))
        self.assertTrue(os.access(os.path.join(self.e_root, "lift", "mm2_cache.sh"), os.X_OK))
        rc, out, err = self.real("build")
        self.assertEqual(rc, 0, err)
        self.assertEqual(rd(os.path.join(self.e_root, "table", "STATUS")).strip(), "VALID")

    def test_a_stub_flag_is_refused_in_this_mode(self):
        rc, out, err = self.real("run", "--minimap2", "/m/mm2")
        self.assertEqual(rc, 2, err)

    def test_build_refuses_a_run_made_with_another_command(self):
        self.run_all()
        kpath = os.path.join(self.e_root, "lift", "orangutan", "runs", "r01", "done.key")
        key = json.loads(rd(kpath))
        key["cmd"][0] = "/some/other/liftoff"
        with open(kpath, "w") as fh:
            json.dump(key, fh)
        rc, out, err = self.real("build")
        self.assertEqual(rc, 3, err)
        self.assertIn("another command", err)
        self.assertIn("r01", err)
        self.assertFalse(os.path.exists(os.path.join(self.e_root, "table")))

    def test_build_refuses_runs_made_from_a_run_gff_that_has_changed_since(self):
        self.run_all()
        with open(os.path.join(self.e_root, "gff", "pos1", "in.gff"), "a") as fh:
            fh.write("# regenerated\n")
        rc, out, err = self.real("build")
        self.assertEqual(rc, 3, err)
        self.assertIn("another command or from another in.gff", err)
        self.assertIn("pos1", err)

    def test_build_refuses_other_code_than_the_runs_were_made_with(self):
        self.run_all()
        ip = os.path.join(self.e_root, "run_identity.json")
        ident = json.loads(rd(ip))
        ident["code_sha256"]["figures/figlib.py"] = "0" * 64
        with open(ip, "w") as fh:
            json.dump(ident, fh)
        rc, out, err = self.real("build")
        self.assertEqual(rc, 3, err)
        self.assertIn("other code", err)

    def test_build_refuses_a_changed_input(self):
        self.run_all()
        with open(self.og["H"]["selflift"], "a") as fh:
            fh.write("\n")
        rc, out, err = self.real("build")
        self.assertEqual(rc, 3, err)
        self.assertIn("H_selflift", err)

    def test_build_of_a_root_without_a_run_identity_is_refused(self):
        self.run_all()
        os.remove(os.path.join(self.e_root, "run_identity.json"))
        rc, out, err = self.real("build")
        self.assertEqual(rc, 3, err)
        self.assertIn("run_identity.json", err)

    def test_the_sensitivity_build_has_a_table_dir_of_its_own_and_checks_the_same_things(self):
        self.run_all()
        rc, out, err = self.real("build", "--strict-single-copy")
        self.assertEqual(rc, 2, err)
        rc, out, err = self.real("build", "--strict-single-copy", "--table-dir", os.path.join(self.e_root, "table_strict"))
        self.assertEqual(rc, 0, err)


class AtomicStateTests(HermeticCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def test_a_state_file_is_written_whole_or_not_at_all(self):
        p = os.path.join(self.tmp, "x", "state.json")
        LC._atomic_write(p, "one\n")
        LC._atomic_write(p, "two\n")
        self.assertEqual(rd(p), "two\n")
        self.assertEqual(os.listdir(os.path.dirname(p)), ["state.json"])

    def test_a_corrupt_splits_json_is_a_stop_with_a_clear_message(self):
        d = os.path.join(self.tmp, "lift", "human")
        os.makedirs(d)
        with open(os.path.join(d, "splits.json"), "w") as fh:
            fh.write("[\"r01\"")
        with self.assertRaises(LC.StopRun) as cm:
            LC.load_splits(self.tmp, "human")
        self.assertIn("splits.json", str(cm.exception))

    def test_save_splits_is_atomic_and_sorted(self):
        LC.save_splits(self.tmp, "human", ["r02", "r01", "r02"])
        self.assertEqual(json.loads(rd(os.path.join(self.tmp, "lift", "human", "splits.json"))), ["r01", "r02"])
        self.assertEqual(sorted(os.listdir(os.path.join(self.tmp, "lift", "human"))), ["splits.json"])

    def test_a_corrupt_run_identity_is_a_stop(self):
        with open(os.path.join(self.tmp, "run_identity.json"), "w") as fh:
            fh.write("{")
        with self.assertRaises(LC.StopRun):
            LC.freeze_identity(self.tmp, {"a.py": "0" * 64}, "abc1234")

    def test_comparing_does_not_write(self):
        LC.freeze_plan(self.tmp, "abc", write=False)
        LC.freeze_identity(self.tmp, {"a.py": "0" * 64}, "abc1234", write=False)
        self.assertEqual(os.listdir(self.tmp), [])


class PlanChecksTests(WorldCase):
    def test_a_plan_with_a_record_named_like_the_copy_of_another_is_refused(self):
        rows, recs = list(self.rows), dict(self.recs)
        recs["gene-V2_1"] = rec("gene-V2_1", "c12", 1000, 3000)
        rows.append(row("gene-V2_1", contig="c12", start0=1000, end=3000))
        cp = row("gene-V2_1", cls="dropped_M1", contig="c9", start0=1000, end=3000, note="overlaps gene-V2 (c9)")
        cp["tag"] = "1"
        rows.append(cp)
        root = os.path.join(self.tmp, "lift_collision")
        write_selflift(root, rows, recs)
        with self.assertRaises(ValueError) as cm:
            LC.plan_from_real_inputs("nodes", lift_root=root, spec=self.spec, n_neg=4)
        self.assertIn("gene-V2_1", str(cm.exception))


class HermeticBaseTests(HermeticCase):
    """the base itself: inside any test the defaults of the run point nowhere, and a command that relies on them fails loudly without touching anything"""

    def test_the_default_paths_of_the_run_cannot_exist(self):
        for path in (LC.E_ROOT, X.LIFT_ROOT, LC.LIFTOFF_BIN, LC.GORILLA_FASTA, LC.GORILLA_GFF, LC.MINIMAP2_ENV, *(v for t in LC.OUTGROUPS.values() for v in (t["genome"], t["mmi"], t["selflift"]))):
            self.assertTrue(path.startswith(VOID + "/"), path)
            self.assertFalse(os.path.lexists(path), path)
        self.assertEqual(sorted(LC.OUTGROUPS), ["C", "H", "O"])

    def test_the_registered_values_are_still_known_to_the_tests_that_pin_them(self):
        self.assertEqual(REAL.E_ROOT, "/mnt/linuxdisk/tmp/hier_run/e")
        self.assertEqual(REAL.LIFT_ROOT, "/mnt/linuxdisk/tmp/rustle_figures/liftoff")
        self.assertTrue(REAL.OUTGROUPS["H"]["genome"].endswith("/winloci_data/chm13v2.0.fa"))

    def test_a_command_that_relies_on_the_defaults_fails_loudly_and_creates_nothing(self):
        for argv in (["plan"], ["run", "--amendment-commit", "abc1234"], ["build"], ["gffs"]):
            out, err = io.StringIO(), io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                rc = LC.main(argv)
            self.assertNotEqual(rc, 0, argv)
            self.assertFalse(os.path.lexists(VOID), argv)                    # nothing could be created below the unreachable root
            self.assertEqual(out.getvalue(), "", argv)

    def test_the_directory_helpers_have_no_default_for_the_run_root(self):
        with self.assertRaises(TypeError):
            LC.ensure_target_links("H")
        with self.assertRaises(TypeError):
            LC.run_dir_of("H", "r01")

    def test_the_functions_that_read_the_self_lift_tables_resolve_their_default_when_they_are_called(self):
        with self.assertRaises(FileNotFoundError) as cm:
            LC.plan_from_real_inputs("nodes")
        self.assertIn(VOID, str(cm.exception))
        with self.assertRaises(FileNotFoundError):
            X.pool_report()
        with self.assertRaises(FileNotFoundError):
            LC.build_all(os.path.join(tempfile.mkdtemp(), "e"), os.path.join(tempfile.mkdtemp(), "out"))


if __name__ == "__main__":
    unittest.main()
