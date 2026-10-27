#!/usr/bin/env python3
"""Hermetic tests of the ideal-arm pilot harness (docs/PREREG_ideal_arm_pilot_2026-10-07.md): sim_genes.py, score_soto_nodes.py, the gates and run_soto_pilot.sh.

    python3 bench/ideal_expression/test_soto_pilot.py          (run it under PYTHONHASHSEED=0 and =1)

Everything is built in temporary directories from a TOY genome (three contigs, a few hundred kb) and a toy GFF with canonical introns; no real genome, no real index, no real
minimap2 (the runner smoke test puts stub scripts first on PATH), no pipeline binary. The only repository files read are the tracked Soto tables and the baseline document.
"""
import collections
import contextlib
import csv
import hashlib
import io
import json
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
import warnings

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, BENCH)
SIM_GENES = os.path.join(HERE, "sim_genes.py")
SCORER = os.path.join(HERE, "score_soto_nodes.py")
RUNNER = os.path.join(HERE, "run_soto_pilot.sh")
S1C = os.path.join(BENCH, "soto", "soto_famCN_S1C.tsv")
FIDELITY_DOC = os.path.join(BENCH, "SOTO_FAMILY_LOCUS_FIDELITY_2026-09-22.md")


def tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def jload(path):
    with open(path) as fh:
        return json.load(fh)


def slurp(path):
    with open(path) as fh:
        return fh.read()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------- the toy world
CONTIG_LEN = {"chr1": 120_000, "chr2": 120_000, "chr3": 450_000}
# gene id -> (contig, strand, name, biotype, {transcript id: [(start0, end) exon blocks]}); a gene whose id is in NONCANONICAL gets "AA...TT" introns
TOY_GENES = {
    "TG_A": ("chr1", "+", "AAA1", "protein_coding", {"TG_A.1": [(1000, 1400), (2000, 2300), (3000, 3500)], "TG_A.2": [(1000, 1400), (3000, 3500)]}),
    "TG_B": ("chr1", "-", "BBB1", "lncRNA", {"TG_B.1": [(5000, 5400), (5800, 6200)]}),
    "TG_N1": ("chr1", "+", "NEIGH1", "processed_pseudogene", {"TG_N1.1": [(8000, 8500)]}),
    "TG_S130": ("chr1", "+", "S130", "lncRNA", {"TG_S130.1": [(31000, 31130)]}),
    "TG_NC": ("chr1", "+", "NONCAN1", "unprocessed_pseudogene", {"TG_NC.1": [(20000, 20400), (20800, 21200)]}),
    "TG_SH": ("chr1", "+", "SHORT1", "transcribed_unprocessed_pseudogene", {"TG_SH.1": [(30000, 30300), (30330, 30600)]}),
    "TG_FAR": ("chr1", "+", "FAR1", "lncRNA", {"TG_FAR.1": [(70000, 70500), (70900, 71400)]}),
    "TG_P": ("chr2", "+", "PARA1", "protein_coding", {"TG_P.1": [(11000, 11400), (12000, 12300), (13000, 13500)], "TG_P.2": [(11000, 11400), (13000, 13500)]}),
    "TG_TINY": ("chr2", "+", "TINY1", "lncRNA", {"TG_TINY.1": [(30000, 30060)]}),
    "TG_BIG": ("chr3", "+", "BIG1", "protein_coding", {"TG_BIG.1": [(1000, 1300), (260000, 260300)], "TG_BIG.2": [(300000, 316000), (317000, 332000)], "TG_BIG.3": [(400000, 400400)]}),
}
NONCANONICAL = {"TG_NC"}
SOTO_LIST = ["TG_A", "TG_B", "TG_NC", "TG_SH", "TG_P", "TG_TINY", "TG_BIG", "TG_MISS"]   # TG_MISS is listed but absent from the annotation
FAMILIES = {"TG_A": ["ID_1"], "TG_P": ["ID_1", "ID_20"], "TG_B": ["ID_2"], "TG_NC": ["ID_3"], "TG_SH": ["ID_3"], "TG_TINY": ["ID_4"], "TG_BIG": ["ID_5"], "TG_MISS": ["ID_6"]}
HALF = {"TG_A": "dev", "TG_B": "heldout", "TG_NC": "dev", "TG_SH": "dev", "TG_P": "dev", "TG_TINY": "dev", "TG_BIG": "heldout", "TG_MISS": "dev"}
SD98 = {"TG_A": "Yes", "TG_P": "Yes", "TG_B": "Yes", "TG_NC": "No", "TG_SH": "No", "TG_TINY": "Yes", "TG_BIG": "No", "TG_MISS": "No"}
SOTO_BIOTYPE = {"TG_A": "protein_coding", "TG_B": "lncRNA", "TG_NC": "unprocessed_pseudogene", "TG_SH": "transcribed_unprocessed_pseudogene", "TG_P": "protein_coding",
                "TG_TINY": "lncRNA", "TG_BIG": "protein_coding", "TG_MISS": "protein_coding"}


def build_world(root):
    """Write the toy FASTA (+fai), the tabix-indexed GFF, a plain gene list, an S1C-like table, a split table; returns a dict of paths."""
    rng = random.Random(20261007)
    seq = {c: list("".join(rng.choices("ACGT", k=n))) for c, n in CONTIG_LEN.items()}
    for gid, (chrom, strand, _name, _bt, txs) in TOY_GENES.items():
        motif = ("AA", "TT") if gid in NONCANONICAL else (("GT", "AG") if strand == "+" else ("CT", "AC"))
        for blocks in txs.values():
            for (_s, pe), (s, _e) in zip(blocks, blocks[1:]):
                if s - pe >= 50:
                    seq[chrom][pe:pe + 2] = list(motif[0])
                    seq[chrom][s - 2:s] = list(motif[1])
    seq["chr2"][10500:14000] = seq["chr1"][500:4000]   # TG_P is a byte-identical paralog of TG_A (same strand, same relative coordinates + 10000)
    fa = os.path.join(root, "toy.fa")
    with open(fa, "w") as fh:
        for c in CONTIG_LEN:
            fh.write(f">{c}\n")
            s = "".join(seq[c])
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    pysam.faidx(fa)
    lines = []
    for gid, (chrom, strand, name, bt, txs) in TOY_GENES.items():
        allb = [b for bl in txs.values() for b in bl]
        lines.append((chrom, min(s for s, _ in allb), 0, f"{chrom}\tCAT\tgene\t{min(s for s, _ in allb) + 1}\t{max(e for _, e in allb)}\t.\t{strand}\t.\tID={gid};Name={gid};gene_id={gid};gene_name={name};gene_biotype={bt}"))
        for tid, bl in txs.items():
            lines.append((chrom, bl[0][0], 1, f"{chrom}\tCAT\ttranscript\t{bl[0][0] + 1}\t{bl[-1][1]}\t.\t{strand}\t.\tID={tid};Parent={gid};gene_id={gid};transcript_biotype={bt}"))
            for i, (s, e) in enumerate(bl, 1):
                lines.append((chrom, s, 2, f"{chrom}\tCAT\texon\t{s + 1}\t{e}\t.\t{strand}\t.\tID={tid}.exon{i};Parent={tid};gene={gid}"))
    order = {c: i for i, c in enumerate(CONTIG_LEN)}
    lines.sort(key=lambda t: (order[t[0]], t[1], t[2], t[3]))
    gff = os.path.join(root, "toy.gff3")
    with open(gff, "w") as fh:
        fh.write("##gff-version 3\n")
        for t in lines:
            fh.write(t[3] + "\n")
    gff_gz = pysam.tabix_index(gff, preset="gff", force=True)
    genes = os.path.join(root, "genes.txt")
    with open(genes, "w") as fh:
        fh.write("# the toy Soto list\n" + "\n".join(SOTO_LIST) + "\n\n")
    s1c = os.path.join(root, "toy_s1c.tsv")
    with open(s1c, "w") as fh:
        fh.write("Gene Name\tGene ID\tBiotype\tIn Table S1 (SD98 gene set)\tFamily ID\n")
        for g in SOTO_LIST:
            nm = TOY_GENES[g][2] if g in TOY_GENES else "MISSING1"
            for fam in FAMILIES[g]:
                fh.write(f"{nm}\t{g}\t{SOTO_BIOTYPE[g]}\t{SD98[g]}\t{fam}\n")
    split = os.path.join(root, "toy_split.tsv")
    with open(split, "w") as fh:
        fh.write("gene_id\tunit_id\thalf\tn_unit_genes\tn_unit_families\tcross_half\n")
        for g in SOTO_LIST:
            fh.write(f"{g}\tU0\t{HALF[g]}\t1\t1\t0\n")
    return dict(fa=fa, gff=gff_gz, genes=genes, s1c=s1c, split=split, root=root)


def sim_args(world, out, **kw):
    """argv for sim_genes.main on the toy world; kw overrides/extends (value None drops the flag, True adds a bare flag)."""
    opts = {"--out": out, "--genes": world["genes"], "--gff": world["gff"], "--genome": world["fa"], "--s1c": world["s1c"], "--split": world["split"],
            "--pad": "2000", "--reads-per-tx": "4", "--seed": "7"}
    for k, v in kw.items():
        opts["--" + k.replace("_", "-")] = v
    argv = []
    for k, v in opts.items():
        if v is None:
            continue
        argv.append(k)
        if v is not True:
            argv.append(str(v))
    return argv


def fastq_records(path):
    with open(path) as fh:
        recs = []
        while True:
            h = fh.readline()
            if not h:
                return recs
            recs.append((h.rstrip("\n")[1:], fh.readline().rstrip("\n"), fh.readline(), fh.readline().rstrip("\n")))


class WorldCase(unittest.TestCase):
    """One toy world per test class."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="soto_pilot_test_")
        cls.world = build_world(cls.tmp)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def run_sim(self, name, **kw):
        import sim_genes
        out = os.path.join(self.tmp, name)
        with contextlib.redirect_stdout(io.StringIO()):
            rc = sim_genes.main(sim_args(self.world, out, **kw))
        self.assertIn(rc, (0, None))
        return out


# ------------------------------------------------------------------------------------------------- sim_genes
class SimGenesTests(WorldCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.sg = __import__("sim_genes")
        cls.out = os.path.join(cls.tmp, "ideal")
        with contextlib.redirect_stdout(io.StringIO()):
            cls.sg.main(sim_args(cls.world, cls.out))

    def test_default_gene_list_is_the_2334_unique_gene_ids_of_s1c(self):
        ids = self.sg.read_gene_list(S1C)
        self.assertEqual(len(ids), 2334)
        self.assertEqual(len(set(ids)), 2334)
        self.assertTrue(all(re.match(r"^(CHM13|LOFF)_G\d+$", g) for g in ids))
        self.assertEqual(os.path.abspath(self.sg.DEFAULT_GENES), os.path.abspath(S1C))

    def test_plain_gene_list_keeps_first_seen_order_and_drops_comments_blanks_duplicates(self):
        p = os.path.join(self.tmp, "plain.txt")
        with open(p, "w") as fh:
            fh.write("# c\nB\n\nA\nB\n  C  \n")
        self.assertEqual(self.sg.read_gene_list(p), ["B", "A", "C"])

    def test_windows_are_gene_spans_plus_pad_merged_and_clipped(self):
        w = [(r["chrom"], int(r["start0"]), int(r["end"])) for r in tsv(self.out + ".windows.tsv")]
        self.assertEqual(w, [("chr1", 0, 8200), ("chr1", 18000, 23200), ("chr1", 28000, 32600), ("chr2", 9000, 15500), ("chr2", 28000, 32060), ("chr3", 0, 402400)])

    def test_pad_zero_windows_are_the_gene_spans(self):
        out = self.run_sim("pad0", pad=0)
        w = [(r["chrom"], int(r["start0"]), int(r["end"])) for r in tsv(out + ".windows.tsv")]
        self.assertEqual(w, [("chr1", 1000, 3500), ("chr1", 5000, 6200), ("chr1", 20000, 21200), ("chr1", 30000, 30600), ("chr2", 11000, 13500), ("chr2", 30000, 30060), ("chr3", 1000, 400400)])
        ids = {r["gene"] for r in tsv(out + ".transcripts.tsv")}
        self.assertNotIn("TG_N1", ids)   # 8000-8500 is outside every pad-0 window

    def test_neighbours_inside_windows_are_simulated_and_far_genes_are_not(self):
        rows = tsv(self.out + ".transcripts.tsv")
        genes = {r["gene"]: r["is_soto"] for r in rows}
        self.assertEqual(genes["TG_N1"], "0")
        self.assertEqual(genes["TG_A"], "1")
        self.assertNotIn("TG_FAR", genes)
        self.assertEqual(sorted(genes), ["TG_A", "TG_B", "TG_BIG", "TG_N1", "TG_NC", "TG_P", "TG_S130", "TG_SH", "TG_TINY"])

    def test_canonicalisation_skip_rules_and_reasons_are_those_of_sim_windows(self):
        by = {r["transcript"]: r for r in tsv(self.out + ".transcripts.tsv")}
        self.assertEqual((by["TG_NC.1"]["canon_blocks"], by["TG_NC.1"]["chain"], by["TG_NC.1"]["n_noncanonical"], by["TG_NC.1"]["spliced_len"]), ("20000-21200", "", "1", "1200"))
        self.assertEqual((by["TG_SH.1"]["canon_blocks"], by["TG_SH.1"]["n_short"]), ("30000-30600", "1"))
        self.assertEqual((by["TG_A.1"]["chain"], by["TG_A.1"]["spliced_len"]), ("1400-2000,2300-3000", "1200"))
        self.assertEqual(by["TG_B.1"]["chain"], "5400-5800")
        self.assertEqual({t: by[t]["skip_reason"] for t in by if by[t]["skip_reason"]},
                         {"TG_TINY.1": "under 120 bp", "TG_BIG.1": "kept intron > 200 kb", "TG_BIG.2": "molecule 31000 bp > 30000"})
        self.assertEqual({t for t in by if by[t]["simulated"] == "1"}, {"TG_A.1", "TG_A.2", "TG_B.1", "TG_N1.1", "TG_NC.1", "TG_SH.1", "TG_S130.1", "TG_P.1", "TG_P.2", "TG_BIG.3"})
        sk = {r["transcript"]: r["reason"] for r in tsv(self.out + ".skipped.tsv")}
        self.assertEqual(sk, {"TG_TINY.1": "under 120 bp", "TG_BIG.1": "kept intron > 200 kb", "TG_BIG.2": "molecule 31000 bp > 30000"})

    def test_ideal_arm_gives_reads_per_tx_to_every_simulated_transcript_and_none_to_the_rest(self):
        rows = tsv(self.out + ".transcripts.tsv")
        self.assertEqual({r["transcript"]: int(r["n_reads"]) for r in rows if r["simulated"] == "1"}, {t: 4 for t in
                         ("TG_A.1", "TG_A.2", "TG_B.1", "TG_N1.1", "TG_NC.1", "TG_SH.1", "TG_S130.1", "TG_P.1", "TG_P.2", "TG_BIG.3")})
        self.assertTrue(all(int(r["n_reads"]) == 0 for r in rows if r["simulated"] == "0"))
        recs = fastq_records(self.out + ".fq")
        self.assertEqual(len(recs), 40)
        names = collections.Counter(n.rsplit("|", 1)[0] for n, *_ in recs)
        self.assertEqual(names["TG_A.1|TG_A|fl"], 4)
        self.assertEqual(sorted(int(n.rsplit("|", 1)[1]) for n, *_ in recs if n.startswith("TG_A.1|")), [0, 1, 2, 3])
        self.assertTrue(all(len(s) == len(q) for _, s, _, q in recs))
        self.assertEqual(jload(self.out + ".sim.json")["reads"], 40)

    def test_targets_table_lists_every_listed_gene_with_status_and_reason(self):
        t = {r["gene"]: r for r in tsv(self.out + ".targets.tsv")}
        self.assertEqual(sorted(t), sorted(SOTO_LIST))
        self.assertEqual({g: t[g]["status"] for g in t}, {"TG_A": "simulated", "TG_B": "simulated", "TG_NC": "simulated", "TG_SH": "simulated", "TG_P": "simulated", "TG_BIG": "simulated",
                                                          "TG_TINY": "skipped", "TG_MISS": "skipped"})
        self.assertEqual(t["TG_MISS"]["reason"], "not in annotation")
        self.assertEqual(t["TG_TINY"]["reason"], "no simulable transcript (under 120 bp: 1)")
        self.assertEqual(t["TG_BIG"]["n_tx_annotated"], "3")
        self.assertEqual((t["TG_BIG"]["n_tx_simulated"], t["TG_BIG"]["n_reads"]), ("1", "4"))
        self.assertEqual((t["TG_A"]["s1c_families"], t["TG_P"]["s1c_families"], t["TG_P"]["n_families"]), ("ID_1", "ID_1;ID_20", "2"))
        self.assertEqual((t["TG_A"]["half"], t["TG_B"]["half"]), ("dev", "heldout"))
        self.assertEqual((t["TG_A"]["biotype"], t["TG_A"]["s1c_biotype"], t["TG_A"]["gene_name"], t["TG_A"]["s1c_name"]), ("protein_coding", "protein_coding", "AAA1", "AAA1"))
        self.assertEqual((t["TG_A"]["chrom"], t["TG_A"]["strand"], t["TG_A"]["start0"], t["TG_A"]["end"], t["TG_A"]["exon_bp"], t["TG_A"]["n_exons"]), ("chr1", "+", "1000", "3500", "1200", "3"))
        self.assertEqual(t["TG_NC"]["single_exon"], "1")
        self.assertEqual(t["TG_A"]["single_exon"], "0")
        self.assertEqual((t["TG_B"]["single_exon"], t["TG_B"]["n_exons"]), ("0", "2"))   # two exons are not a single exon
        self.assertEqual(t["TG_TINY"]["single_exon"], "")

    def test_sim_json_counts_and_fastq_checksum(self):
        m = jload(self.out + ".sim.json")
        self.assertEqual((m["genes_listed"], m["genes_in_annotation"], m["genes_missing"]), (8, 7, 1))
        self.assertEqual((m["transcripts"], m["simulated_transcripts"], m["reads"], m["windows"]), (13, 10, 40, 6))
        self.assertEqual((m["genes_in_windows"], m["neighbour_genes"]), (9, 2))
        self.assertEqual(m["fastq_sha256"], sha256(self.out + ".fq"))
        self.assertEqual(m["seed"], 7)
        self.assertEqual(m["skipped"], {"under 120 bp": 1, "kept intron > 200 kb": 1, "molecule > 30000 bp": 1})

    def test_read_ends_vary_and_stay_within_the_jitter(self):
        import lib
        out = self.run_sim("err0", err=0)   # error-free reads are exact substrings of the transcript, so the trims can be read off
        fa = pysam.FastaFile(self.world["fa"])
        rows = {r["transcript"]: r for r in tsv(out + ".transcripts.tsv") if r["simulated"] == "1"}
        recs = fastq_records(out + ".fq")
        five, three = collections.defaultdict(set), collections.defaultdict(set)
        for name, seq, _, _ in recs:
            t = name.split("|")[0]
            r = rows[t]
            cb = [tuple(map(int, b.split("-"))) for b in r["canon_blocks"].split(",")]
            body = lib.spliced1(fa, r["chrom"], [(s + 1, e) for s, e in cb], r["strand"])
            lo = body.find(seq)
            self.assertGreaterEqual(lo, 0, f"{name} is not a substring of its transcript")
            hi = len(body) - lo - len(seq)
            self.assertTrue(0 <= lo <= 30 and 0 <= hi <= 30, (name, lo, hi))
            five[t].add(lo)
            three[t].add(hi)
        self.assertTrue(any(len(v) > 1 for v in five.values()))
        self.assertTrue(any(len(v) > 1 for v in three.values()))
        self.assertGreater(len(set().union(*five.values())), 3)
        self.assertGreater(len(set().union(*three.values())), 3)

    def test_labels_stay_with_the_source_copy_when_paralogs_are_identical(self):
        import lib
        fa = pysam.FastaFile(self.world["fa"])
        rows = {r["transcript"]: r for r in tsv(self.out + ".transcripts.tsv")}
        body = {}
        for t in ("TG_A.1", "TG_P.1"):
            cb = [tuple(map(int, b.split("-"))) for b in rows[t]["canon_blocks"].split(",")]
            body[t] = lib.spliced1(fa, rows[t]["chrom"], [(s + 1, e) for s, e in cb], rows[t]["strand"])
        self.assertEqual(body["TG_A.1"], body["TG_P.1"], "the toy paralogs must be sequence-identical for this test to mean anything")
        for name, *_ in fastq_records(self.out + ".fq"):
            t, g, kind, _k = name.split("|")
            self.assertEqual(kind, "fl")
            self.assertEqual(g, rows[t]["gene"])
        self.assertEqual(sum(1 for n, *_ in fastq_records(self.out + ".fq") if n.split("|")[1] == "TG_P"), 8)
        self.assertEqual(sum(1 for n, *_ in fastq_records(self.out + ".fq") if n.split("|")[1] == "TG_A"), 8)

    def test_reads_of_a_minus_strand_gene_are_reverse_complemented_transcripts(self):
        import lib
        out = self.run_sim("err0b", err=0)
        fa = pysam.FastaFile(self.world["fa"])
        r = {x["transcript"]: x for x in tsv(out + ".transcripts.tsv")}["TG_B.1"]
        cb = [tuple(map(int, b.split("-"))) for b in r["canon_blocks"].split(",")]
        plus = "".join(fa.fetch(r["chrom"], s, e) for s, e in cb).upper()
        reads = [s for n, s, _, _ in fastq_records(out + ".fq") if n.startswith("TG_B.1|")]
        self.assertEqual(len(reads), 4)
        self.assertTrue(all(lib.rc(s) in plus for s in reads) and not any(s in plus for s in reads))

    def test_same_seed_same_bytes_different_seed_different_reads_same_tables(self):
        a, b = self.run_sim("s7a", seed=7), self.run_sim("s8", seed=8)
        self.assertEqual(sha256(a + ".fq"), sha256(self.out + ".fq"))
        self.assertNotEqual(sha256(b + ".fq"), sha256(a + ".fq"))
        for suf in (".transcripts.tsv", ".targets.tsv", ".windows.tsv", ".skipped.tsv"):
            self.assertEqual(sha256(a + suf), sha256(b + suf), suf)

    def test_no_reads_writes_the_tables_but_no_fastq(self):
        out = self.run_sim("noreads", no_reads=True)
        self.assertFalse(os.path.exists(out + ".fq"))
        for suf in (".transcripts.tsv", ".targets.tsv", ".windows.tsv", ".skipped.tsv", ".sim.json"):
            self.assertTrue(os.path.exists(out + suf), suf)
        self.assertEqual(sha256(out + ".transcripts.tsv"), sha256(self.out + ".transcripts.tsv"))

    def test_output_is_byte_identical_under_hash_seeds_0_and_1(self):
        outs = []
        for hs in ("0", "1"):
            out = os.path.join(self.tmp, f"hs{hs}")
            subprocess.run([sys.executable, SIM_GENES] + sim_args(self.world, out), check=True, capture_output=True, env={**os.environ, "PYTHONHASHSEED": hs})
            outs.append(out)
        for suf in (".fq", ".transcripts.tsv", ".targets.tsv", ".windows.tsv", ".skipped.tsv", ".sim.json"):
            self.assertEqual(sha256(outs[0] + suf), sha256(outs[1] + suf), suf)
        self.assertEqual(sha256(outs[0] + ".fq"), sha256(self.out + ".fq"))

    def test_processes_do_not_change_the_output(self):
        out = self.run_sim("procs2", procs=2)
        for suf in (".fq", ".transcripts.tsv", ".targets.tsv", ".sim.json"):
            self.assertEqual(sha256(out + suf), sha256(self.out + suf), suf)

    def test_first_restricts_the_gene_list_to_its_first_n_genes(self):
        out = self.run_sim("first3", first=3)
        self.assertEqual([r["gene"] for r in tsv(out + ".targets.tsv")], ["TG_A", "TG_B", "TG_NC"])
        self.assertEqual(jload(out + ".sim.json")["genes_listed"], 3)
        self.assertNotIn("TG_P", {r["gene"] for r in tsv(out + ".transcripts.tsv")})

    def test_the_registered_seed_is_the_default_and_differs_from_the_earlier_replicates(self):
        self.assertEqual(self.sg.SEED, 20261008)
        self.assertNotIn(self.sg.SEED, (20261006, 20261007))
        self.assertIn(f"[--seed {self.sg.SEED}]", self.sg.__doc__)      # the usage line (and --help) must not advertise another default

    def test_targets_carry_arm_independent_truth_columns(self):
        t = {r["gene"]: r for r in tsv(self.out + ".targets.tsv")}
        self.assertEqual({g: t[g]["simulable"] for g in t}, {"TG_A": "1", "TG_B": "1", "TG_NC": "1", "TG_SH": "1", "TG_P": "1", "TG_BIG": "1", "TG_TINY": "0", "TG_MISS": "0"})
        self.assertEqual({g: t[g]["sd98"] for g in ("TG_A", "TG_NC", "TG_TINY")}, {"TG_A": "Yes", "TG_NC": "No", "TG_TINY": "Yes"})
        self.assertEqual({g: t[g]["identical_seq_genes"] for g in t}, {"TG_A": "1", "TG_P": "1", "TG_B": "0", "TG_NC": "0", "TG_SH": "0", "TG_BIG": "0", "TG_TINY": "0", "TG_MISS": "0"})
        self.assertEqual(t["TG_A"]["source"], "")   # the toy ids are neither CHM13_ nor LOFF_

    def test_source_is_cat_or_liftoff_by_the_gene_id_prefix(self):
        self.assertEqual((self.sg.source_of("CHM13_G0000001"), self.sg.source_of("LOFF_G0000001"), self.sg.source_of("X")), ("CAT", "Liftoff", ""))

    def test_the_truth_columns_do_not_depend_on_the_reads_of_the_arm(self):
        import json as _json
        w = os.path.join(self.tmp, "w_truth.tsv")
        with open(w, "w") as fh:
            fh.write("TG_A\t3\nTG_B\t0\nTG_NC\t0\n")
        weighted = self.run_sim("w_truth", weights=w, total_reads=30)
        a = {r["gene"]: r for r in tsv(self.out + ".targets.tsv")}
        b = {r["gene"]: r for r in tsv(weighted + ".targets.tsv")}
        for g in ("TG_B", "TG_NC"):
            self.assertEqual((b[g]["status"], b[g]["n_reads"], b[g]["reason"]), ("skipped", "0", "weight 0"))
            for col in ("simulable", "exon_bp", "n_exons", "single_exon", "n_tx_eligible", "identical_seq_genes"):
                self.assertEqual(b[g][col], a[g][col], (g, col))
        self.assertEqual((b["TG_NC"]["exon_bp"], b["TG_NC"]["single_exon"]), ("1200", "1"))

    def test_single_exon_means_every_eligible_transcript_is_unspliced(self):
        t = {r["gene"]: r for r in tsv(self.out + ".targets.tsv")}
        self.assertEqual({g: t[g]["single_exon"] for g in ("TG_A", "TG_B", "TG_NC", "TG_SH", "TG_BIG", "TG_TINY")}, {"TG_A": "0", "TG_B": "0", "TG_NC": "1", "TG_SH": "1", "TG_BIG": "1", "TG_TINY": ""})

    def test_a_gene_whose_eligible_transcripts_are_all_unspliced_is_single_exon_even_with_two_disjoint_blocks(self):
        def row(t, g, blocks, chain):
            return dict(transcript=t, gene=g, canon_blocks=blocks, chain=chain, eligible=1, simulated=1, n_reads=4, skip_reason="", seq_sha1="")
        rows = [row("X.1", "X", "100-200", ""), row("X.2", "X", "500-700", ""), row("Y.1", "Y", "100-200", ""), row("Y.2", "Y", "100-200,500-700", "200-500")]
        info = {g: dict(chrom="chr1", strand="+", name=g, biotype="lncRNA", start0=100, end=700, n_tx=2) for g in "XY"}
        t = {r["gene"]: r for r in self.sg.target_rows(["X", "Y"], info, rows, {}, {})}
        self.assertEqual((t["X"]["single_exon"], t["X"]["n_exons"]), (1, 2))    # two blocks, but no transcript is spliced
        self.assertEqual(t["Y"]["single_exon"], 0)                              # one spliced isoform makes it a multi-exon gene

    def test_sim_json_documents_the_genes_that_cannot_be_simulated_and_the_identical_sequences(self):
        m = jload(self.out + ".sim.json")
        self.assertEqual(m["not_simulable"], {"genes": 2, "sd98_eligible": 1, "by_reason": {"under 120 bp": 1, "not in annotation": 1}})
        self.assertEqual(m["identical_sequences"], {"sequences_in_2plus_genes": 2, "genes": 2, "listed_genes": 2, "sequences_in_2plus_listed_genes": 2, "listed_genes_in_those": 2})
        rows = tsv(self.out + ".identical_seq.tsv")
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(r["genes"] == "TG_A;TG_P" and r["n_genes"] == "2" for r in rows))

    def test_read_bases_equal_the_bases_of_the_fastq(self):
        m = jload(self.out + ".sim.json")
        self.assertEqual(m["read_bases"], sum(len(s) for _, s, _, _ in fastq_records(self.out + ".fq")))

    def test_readnames_lists_every_read_in_file_order(self):
        names = os.path.join(self.tmp, "names.txt")
        r = subprocess.run([sys.executable, SIM_GENES, "readnames", "--out", self.out, "--to", names], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(slurp(names).split("\n")[:-1], [n for n, *_ in fastq_records(self.out + ".fq")])

    def test_transcripts_table_carries_the_sequence_hash_of_eligible_transcripts_only(self):
        rows = {r["transcript"]: r for r in tsv(self.out + ".transcripts.tsv")}
        self.assertEqual(rows["TG_A.1"]["seq_sha1"], rows["TG_P.1"]["seq_sha1"])
        self.assertNotEqual(rows["TG_A.1"]["seq_sha1"], rows["TG_A.2"]["seq_sha1"])
        self.assertEqual((rows["TG_TINY.1"]["seq_sha1"], rows["TG_BIG.1"]["seq_sha1"]), ("", ""))
        self.assertEqual(len(rows["TG_A.1"]["seq_sha1"]), 40)

    def test_more_than_four_processes_are_refused(self):
        import sim_genes
        with self.assertRaises(SystemExit):
            with contextlib.redirect_stderr(io.StringIO()):
                sim_genes.main(sim_args(self.world, os.path.join(self.tmp, "p9"), procs=9))


class SkipRuleTests(unittest.TestCase):
    """The simulate/skip rule of sim_windows.main, restated in sim_genes.skip_reason, at its boundaries."""

    def setUp(self):
        self.sg = __import__("sim_genes")

    def test_spliced_length_boundary_is_120_bp(self):
        self.assertEqual(self.sg.skip_reason(119, []), "under 120 bp")
        self.assertEqual(self.sg.skip_reason(120, []), "")

    def test_molecule_boundary_is_30_kb(self):
        self.assertEqual(self.sg.skip_reason(30_000, []), "")
        self.assertEqual(self.sg.skip_reason(30_001, []), "molecule 30001 bp > 30000")

    def test_kept_intron_boundary_is_200_kb(self):
        self.assertEqual(self.sg.skip_reason(500, [(1000, 201_000)]), "")
        self.assertEqual(self.sg.skip_reason(500, [(1000, 201_001)]), "kept intron > 200 kb")

    def test_the_constants_are_those_of_sim_windows(self):
        import sim_windows as sw
        self.assertEqual((sw.MAX_MOLECULE, sw.MAX_INTRON, sw.MIN_INTRON), (30_000, 200_000, 50))

    def test_reason_classes_group_the_molecule_lengths(self):
        self.assertEqual(self.sg.reason_class("molecule 31000 bp > 30000"), "molecule > 30000 bp")
        self.assertEqual(self.sg.reason_class("under 120 bp"), "under 120 bp")


class WindowArithmeticTests(unittest.TestCase):
    def setUp(self):
        self.sg = __import__("sim_genes")

    def test_touching_and_overlapping_windows_merge_and_separate_ones_do_not(self):
        self.assertEqual(self.sg.merge_windows([(0, 10), (10, 20)], 0), [(0, 20)])
        self.assertEqual(self.sg.merge_windows([(0, 10), (11, 20)], 0), [(0, 10), (11, 20)])
        self.assertEqual(self.sg.merge_windows([(100, 200), (150, 400)], 10), [(90, 410)])

    def test_windows_are_padded_clipped_at_zero_and_at_the_contig_end(self):
        self.assertEqual(self.sg.merge_windows([(100, 200)], 500), [(0, 700)])
        self.assertEqual(self.sg.merge_windows([(100, 200)], 500, limit=300), [(0, 300)])

    def test_input_order_is_irrelevant(self):
        self.assertEqual(self.sg.merge_windows([(50, 60), (0, 10), (5, 20)], 0), [(0, 20), (50, 60)])


class ParityWithSimWindowsTests(WorldCase):
    """sim_genes on the toy genes of chr1 must reproduce sim_windows.py (the earlier ideal harness) transcript row for row and read for read."""

    def test_rows_and_reads_equal_sim_windows(self):
        import sim_genes
        import sim_windows as sw
        copies = os.path.join(self.tmp, "ann")
        os.makedirs(copies)
        with open(os.path.join(copies, "copies.hsa.tsv"), "w") as fh:
            fh.write("family\tchrom\tterr_lo0\tterr_hi\tisoform_gene\tcid\tstrand\treadthrough\tname\tcat_name\n")
            for i, g in enumerate(("TG_A", "TG_B", "TG_SH")):
                chrom, strand, name, _bt, txs = TOY_GENES[g]
                blocks = [b for bl in txs.values() for b in bl]
                fh.write(f"NPIP\tchr1\t{min(s for s, _ in blocks)}\t{max(e for _, e in blocks)}\t{g}\th{i}\t{strand}\t0\t{name}\t{name}\n")
        old = (sw.CAT_GFF, sw.GENOME, sw.ANN, dict(sw.CHROM))
        try:
            sw.CAT_GFF, sw.GENOME, sw.ANN = self.world["gff"], self.world["fa"], copies
            sw.CHROM["NPIP"] = "chr1"
            ref = os.path.join(self.tmp, "ref")
            argv = sys.argv
            sys.argv = ["sim_windows.py", "--family", "NPIP", "--out", ref, "--pad", "2000", "--reads-per-tx", "4", "--seed", "7"]
            try:
                with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
                    warnings.simplefilter("ignore", ResourceWarning)   # sim_windows.main leaves files open
                    sw.main()
            finally:
                sys.argv = argv
        finally:
            sw.CAT_GFF, sw.GENOME, sw.ANN = old[0], old[1], old[2]
            sw.CHROM.clear(); sw.CHROM.update(old[3])
        genes = os.path.join(self.tmp, "chr1_genes.txt")
        with open(genes, "w") as fh:
            fh.write("TG_A\nTG_B\nTG_SH\n")
        out = os.path.join(self.tmp, "mine")
        with contextlib.redirect_stdout(io.StringIO()):
            sim_genes.main(sim_args(self.world, out, genes=genes))
        mine = {r["transcript"]: r for r in tsv(out + ".transcripts.tsv")}
        theirs = {r["transcript"]: r for r in tsv(ref + ".transcripts.tsv")}
        self.assertEqual(sorted(mine), sorted(theirs))
        self.assertIn("TG_S130.1", mine)   # the 130 bp neighbour: most of its reads fall back to the untrimmed transcript
        for t in theirs:
            for col in ("gene", "chrom", "strand", "raw_blocks", "canon_blocks", "chain", "n_raw_introns", "n_noncanonical", "n_short", "spliced_len", "simulated", "skip_reason", "biotype", "gene_name"):
                self.assertEqual(mine[t][col], theirs[t][col], (t, col))
        self.assertEqual(slurp(out + ".fq"), slurp(ref + ".fq"))


# ------------------------------------------------------------------------------------------------- weights
class AllocationTests(unittest.TestCase):
    def setUp(self):
        self.sg = __import__("sim_genes")

    def test_largest_remainder_sums_to_the_total_and_is_proportional(self):
        a = self.sg.allocate_weighted({"a": 3.0, "b": 1.0, "c": 0.0, "d": 2.0}, 10)
        self.assertEqual(a, {"a": 5, "b": 2, "c": 0, "d": 3})

    def test_remainder_ties_go_to_the_smaller_gene_id_and_input_order_is_irrelevant(self):
        self.assertEqual(self.sg.allocate_weighted({"c": 1.0, "a": 1.0, "b": 1.0}, 2), {"a": 1, "b": 1, "c": 0})
        w = {f"g{i}": float(i % 5 + 1) for i in range(40)}
        a = self.sg.allocate_weighted(w, 977)
        keys = list(w)
        random.Random(3).shuffle(keys)
        self.assertEqual(a, self.sg.allocate_weighted({k: w[k] for k in keys}, 977))
        self.assertEqual(sum(a.values()), 977)

    def test_all_zero_weights_cannot_be_allocated(self):
        with self.assertRaises(ValueError):
            self.sg.allocate_weighted({"a": 0.0, "b": 0.0}, 5)

    def test_negative_or_nan_weights_are_refused(self):
        for bad in (-1.0, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                self.sg.allocate_weighted({"a": 1.0, "b": bad}, 5)

    def test_split_equal_gives_the_remainder_to_the_first_transcripts(self):
        self.assertEqual(self.sg.split_equal(7, 3), [3, 2, 2])
        self.assertEqual(self.sg.split_equal(2, 3), [1, 1, 0])
        self.assertEqual(self.sg.split_equal(0, 3), [0, 0, 0])


class WeightsFileTests(unittest.TestCase):
    def setUp(self):
        self.sg = __import__("sim_genes")
        self.tmp = tempfile.mkdtemp(prefix="soto_weights_test_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def write(self, text):
        p = os.path.join(self.tmp, "w.tsv")
        with open(p, "w") as fh:
            fh.write(text)
        return p

    def test_header_comments_blank_lines_and_space_separated_lines_are_accepted(self):
        self.assertEqual(self.sg.read_weights(self.write("gene_id\tweight\n# c\n\nA\t1.5\nB 2\n")), {"A": 1.5, "B": 2.0})
        self.assertEqual(self.sg.read_weights(self.write("A\t3\nB\t0\n")), {"A": 3.0, "B": 0.0})

    def test_a_repeated_gene_a_negative_a_non_finite_or_a_non_numeric_weight_is_an_error(self):
        for text in ("A\t1\nA\t2\n", "A\t-1\n", "A\tnan\n", "A\tinf\n", "A\t1\nB\tabc\n", "A\n"):
            with self.assertRaises(ValueError, msg=text):
                self.sg.read_weights(self.write(text))


class WeightedArmTests(WorldCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.weights = os.path.join(cls.tmp, "w.tsv")
        with open(cls.weights, "w") as fh:
            fh.write("gene_id\tweight\nTG_A\t3\nTG_B\t1\nTG_NC\t0\nTG_SH\t2\nTG_P\t4\nTG_BIG\t1\n")

    def test_reads_follow_the_weights_with_largest_remainder_rounding_and_equal_split(self):
        out = self.run_sim("w100", weights=self.weights, total_reads=100)
        rows = tsv(out + ".transcripts.tsv")
        per_gene = collections.Counter()
        for r in rows:
            per_gene[r["gene"]] += int(r["n_reads"])
        # simulable genes and weights: A 3, B 1, N1 and S130 (neighbours: the median of the positive Soto weights = 2), SH 2, P 4, BIG 1 -> total 15; NC has weight 0.
        # quotas 20, 6.67, 13.33, 13.33, 13.33, 26.67, 6.67: floors sum to 97, the 3 left go to the three 2/3 remainders (B, BIG, P: tie, smaller gene id first, all three fit)
        self.assertEqual(dict(per_gene), {"TG_A": 20, "TG_B": 7, "TG_N1": 13, "TG_S130": 13, "TG_NC": 0, "TG_SH": 13, "TG_P": 27, "TG_BIG": 7, "TG_TINY": 0})
        self.assertEqual(sum(per_gene.values()), 100)
        by = {r["transcript"]: int(r["n_reads"]) for r in rows}
        self.assertEqual((by["TG_A.1"], by["TG_A.2"], by["TG_P.1"], by["TG_P.2"]), (10, 10, 14, 13))
        self.assertEqual(len(fastq_records(out + ".fq")), 100)
        self.assertEqual(jload(out + ".sim.json")["reads"], 100)

    def test_zero_weight_and_unreadable_genes_are_skipped_with_a_reason(self):
        out = self.run_sim("w100b", weights=self.weights, total_reads=100)
        t = {r["gene"]: r for r in tsv(out + ".targets.tsv")}
        self.assertEqual((t["TG_NC"]["status"], t["TG_NC"]["reason"]), ("skipped", "weight 0"))
        self.assertEqual(t["TG_TINY"]["status"], "skipped")
        self.assertEqual((t["TG_A"]["status"], t["TG_A"]["n_reads"], t["TG_A"]["n_tx_simulated"]), ("simulated", "20", "2"))
        rows = {r["transcript"]: r for r in tsv(out + ".transcripts.tsv")}
        self.assertEqual(rows["TG_NC.1"]["simulated"], "0")

    def test_neighbour_weight_is_documented_in_the_json_and_can_be_set(self):
        out = self.run_sim("w100c", weights=self.weights, total_reads=100)
        m = jload(out + ".sim.json")["weights"]
        self.assertEqual((m["neighbour_weight"], m["neighbour_weight_rule"], m["soto_with_weight"], m["soto_missing"]), (2.0, "median of the positive Soto weights", 6, 2))
        out2 = self.run_sim("w100d", weights=self.weights, total_reads=100, neighbour_weight="median-all")
        self.assertEqual(jload(out2 + ".sim.json")["weights"]["neighbour_weight"], 1.5)
        out3 = self.run_sim("w100e", weights=self.weights, total_reads=100, neighbour_weight="0")
        rows = tsv(out3 + ".transcripts.tsv")
        self.assertEqual(sum(int(r["n_reads"]) for r in rows if r["gene"] == "TG_N1"), 0)
        self.assertEqual(sum(int(r["n_reads"]) for r in rows), 100)

    def test_missing_weight_option_gives_listed_genes_without_a_line_a_weight(self):
        out = self.run_sim("w100f", weights=self.weights, total_reads=100, missing_weight="1")
        t = {r["gene"]: r for r in tsv(out + ".targets.tsv")}
        self.assertEqual(t["TG_TINY"]["status"], "skipped")   # still nothing to simulate
        self.assertEqual(t["TG_A"]["status"], "simulated")

    def test_a_bad_neighbour_weight_is_a_clean_error(self):
        import sim_genes
        for bad in ("abc", "-1", "nan"):
            with self.assertRaises(SystemExit) as cm:
                with contextlib.redirect_stderr(io.StringIO()):
                    sim_genes.main(sim_args(self.world, os.path.join(self.tmp, "badnw"), weights=self.weights, total_reads=100, neighbour_weight=bad))
            self.assertNotEqual(cm.exception.code, 0, bad)
            self.assertIn("--neighbour-weight", str(cm.exception.code), bad)

    def test_positive_weights_of_genes_without_a_simulable_transcript_are_counted(self):
        w = os.path.join(self.tmp, "w_tiny.tsv")
        with open(w, "w") as fh:
            fh.write("TG_A\t3\nTG_TINY\t5\nTG_MISS\t1\n")
        out = self.run_sim("w_tiny", weights=w, total_reads=50)
        m = jload(out + ".sim.json")["weights"]
        self.assertEqual(m["listed_with_positive_weight_and_no_simulable_transcript"], 2)   # TG_TINY (all transcripts rejected) and TG_MISS (not in the annotation)
        t = {r["gene"]: r for r in tsv(out + ".targets.tsv")}
        self.assertEqual(t["TG_TINY"]["status"], "skipped")
        self.assertEqual(sum(int(r["n_reads"]) for r in t.values()), 50 - sum(int(r["n_reads"]) for r in tsv(out + ".transcripts.tsv") if r["is_soto"] == "0"))

    def test_weights_need_a_total_and_a_total_needs_weights(self):
        import sim_genes
        for kw in ({"weights": self.weights}, {"total_reads": 10}):
            with self.assertRaises(SystemExit):
                with contextlib.redirect_stderr(io.StringIO()):
                    sim_genes.main(sim_args(self.world, os.path.join(self.tmp, "bad"), **kw))

    def test_weighted_output_is_byte_identical_under_hash_seeds_0_and_1(self):
        outs = []
        for hs in ("0", "1"):
            out = os.path.join(self.tmp, f"whs{hs}")
            subprocess.run([sys.executable, SIM_GENES] + sim_args(self.world, out, weights=self.weights, total_reads=100), check=True, capture_output=True, env={**os.environ, "PYTHONHASHSEED": hs})
            outs.append(out)
        for suf in (".fq", ".transcripts.tsv", ".targets.tsv", ".sim.json"):
            self.assertEqual(sha256(outs[0] + suf), sha256(outs[1] + suf), suf)


# ------------------------------------------------------------------------------------------------- gates
def copy_sim(src, dst, suffixes=(".fq", ".transcripts.tsv", ".targets.tsv", ".windows.tsv", ".skipped.tsv", ".sim.json")):
    for suf in suffixes:
        shutil.copyfile(src + suf, dst + suf)


class GateTests(WorldCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.sg = __import__("sim_genes")
        cls.out = os.path.join(cls.tmp, "g")
        with contextlib.redirect_stdout(io.StringIO()):
            cls.sg.main(sim_args(cls.world, cls.out, reads_per_tx=10))

    def fresh(self, name):
        dst = os.path.join(self.tmp, name)
        copy_sim(self.out, dst)
        return dst

    def g1(self, prefix, genes=True, genome=True):
        return self.sg.gate_manifest(prefix, genes=self.world["genes"] if genes else None, genome=self.world["fa"] if genome else None)

    def test_g1_passes_on_a_good_simulation(self):
        ok, lines = self.g1(self.out)
        self.assertTrue(ok, lines)
        self.assertTrue(any(ln.startswith("G1 PASS") for ln in lines), lines)

    def test_g1_fails_when_a_read_is_missing(self):
        d = self.fresh("g1_missing")
        recs = fastq_records(d + ".fq")
        with open(d + ".fq", "w") as fh:
            for n, s, plus, q in recs[1:]:
                fh.write(f"@{n}\n{s}\n{plus}{q}\n")
        ok, lines = self.g1(d)
        self.assertFalse(ok)
        self.assertTrue(any(ln.startswith("G1 FAIL") and "reads" in ln for ln in lines), lines)

    def test_g1_fails_when_sim_json_disagrees_with_the_reads(self):
        d = self.fresh("g1_json")
        m = jload(d + ".sim.json")
        m["reads"] += 1
        with open(d + ".sim.json", "w") as fh:
            json.dump(m, fh)
        ok, lines = self.g1(d)
        self.assertFalse(ok)
        self.assertTrue(any("sim.json" in ln for ln in lines), lines)

    def test_g1_fails_when_a_read_names_a_transcript_that_is_not_listed(self):
        d = self.fresh("g1_name")
        text = slurp(d + ".fq").replace("@TG_B.1|TG_B|fl|0", "@TG_ZZ.1|TG_B|fl|0", 1)
        with open(d + ".fq", "w") as fh:
            fh.write(text)
        ok, lines = self.g1(d)
        self.assertFalse(ok)
        self.assertTrue(any("not a listed transcript" in ln for ln in lines), lines)

    def test_g1_fails_when_a_read_carries_the_wrong_gene(self):
        d = self.fresh("g1_gene")
        text = slurp(d + ".fq").replace("@TG_B.1|TG_B|fl|0", "@TG_B.1|TG_A|fl|0", 1)
        with open(d + ".fq", "w") as fh:
            fh.write(text)
        ok, lines = self.g1(d)
        self.assertFalse(ok)

    def test_g1_fails_when_a_listed_gene_has_no_target_row(self):
        d = self.fresh("g1_target")
        rows = [r for r in tsv(d + ".targets.tsv") if r["gene"] != "TG_SH"]
        with open(d + ".targets.tsv", "w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        ok, lines = self.g1(d)
        self.assertFalse(ok)
        self.assertTrue(any("TG_SH" in ln for ln in lines), lines)

    def test_g1_fails_when_a_skipped_gene_has_no_reason(self):
        d = self.fresh("g1_reason")
        rows = tsv(d + ".targets.tsv")
        for r in rows:
            if r["gene"] == "TG_MISS":
                r["reason"] = ""
        with open(d + ".targets.tsv", "w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        ok, lines = self.g1(d)
        self.assertFalse(ok)
        self.assertTrue(any("TG_MISS" in ln and "reason" in ln for ln in lines), lines)

    def test_g1_fails_on_a_kept_intron_that_is_not_canonical(self):
        d = self.fresh("g1_intron")
        text = slurp(d + ".transcripts.tsv").replace("1400-2000,2300-3000", "1500-1900,2300-3000", 1)
        with open(d + ".transcripts.tsv", "w") as fh:
            fh.write(text)
        ok, lines = self.g1(d)
        self.assertFalse(ok)
        self.assertTrue(any("intron" in ln and ln.startswith("G1 FAIL") for ln in lines), lines)

    def test_g2_passes_when_read_ends_vary(self):
        ok, lines = self.sg.gate_read_ends(self.out, self.world["fa"])
        self.assertTrue(ok, lines)
        self.assertTrue(lines[0].startswith("G2 PASS"), lines)

    def test_g2_fails_when_no_end_is_trimmed(self):
        out = self.run_sim("g2_nojit", reads_per_tx=10, jitter=0)
        ok, lines = self.sg.gate_read_ends(out, self.world["fa"])
        self.assertFalse(ok)
        self.assertTrue(lines[0].startswith("G2 FAIL"), lines)

    def test_g2_fails_on_pooled_constant_ends_even_when_no_transcript_has_enough_reads_for_the_per_transcript_test(self):
        out = self.run_sim("g2_nojit4", reads_per_tx=4, jitter=0)
        ok, lines = self.sg.gate_read_ends(out, self.world["fa"])
        self.assertFalse(ok)
        self.assertTrue(any("distinct value" in ln for ln in lines), lines)

    def test_g2_fails_when_one_transcript_has_identical_ends_on_all_its_reads(self):
        d = self.fresh("g2_one")
        recs = fastq_records(d + ".fq")
        first = next(s for n, s, _, _ in recs if n.startswith("TG_N1.1|"))
        with open(d + ".fq", "w") as fh:
            for n, s, plus, q in recs:
                if n.startswith("TG_N1.1|"):
                    s, q = first, "~" * len(first)
                fh.write(f"@{n}\n{s}\n{plus}{q}\n")
        ok, lines = self.sg.gate_read_ends(d, self.world["fa"])
        self.assertFalse(ok)
        self.assertTrue(any("constant" in ln and "TG_N1.1" in ln for ln in lines), lines)

    def test_g3_compare_accepts_identical_simulations_and_names_the_file_that_differs(self):
        a, b = self.fresh("g3a"), self.fresh("g3b")
        ok, lines = self.sg.gate_compare(a, b)
        self.assertTrue(ok, lines)
        self.assertTrue(lines[0].startswith("G3 PASS"), lines)
        with open(b + ".targets.tsv", "a") as fh:
            fh.write("x\n")
        ok, lines = self.sg.gate_compare(a, b)
        self.assertFalse(ok)
        self.assertTrue(lines[0].startswith("G3 FAIL") and ".targets.tsv" in lines[0], lines)
        os.remove(b + ".skipped.tsv")
        ok, lines = self.sg.gate_compare(a, b)
        self.assertFalse(ok)
        self.assertIn(".skipped.tsv", lines[0])

    def test_g3_compare_can_be_restricted_to_some_files(self):
        a = self.fresh("g3c_a")
        b = os.path.join(self.tmp, "g3c_b")
        with contextlib.redirect_stdout(io.StringIO()):
            self.sg.main(sim_args(self.world, b, reads_per_tx=10, no_reads=True))
        ok, lines = self.sg.gate_compare(a, b)   # a has a FASTQ, b has none
        self.assertFalse(ok)
        ok, lines = self.sg.gate_compare(a, b, suffixes=(".transcripts.tsv", ".targets.tsv", ".windows.tsv", ".skipped.tsv"))
        self.assertTrue(ok, lines)
        r = subprocess.run([sys.executable, SIM_GENES, "compare", "--only", ".transcripts.tsv,.targets.tsv", a, b], capture_output=True, text=True)
        self.assertEqual((r.returncode, r.stdout.startswith("G3 PASS")), (0, True), r.stdout + r.stderr)

    def test_g3_compare_can_write_its_verdict_as_json(self):
        a, b = self.fresh("g3j_a"), self.fresh("g3j_b")
        j = os.path.join(self.tmp, "g3j.json")
        r = subprocess.run([sys.executable, SIM_GENES, "compare", "--json", j, a, b], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(jload(j)["ok"], True)
        self.assertTrue(jload(j)["lines"][0].startswith("G3 PASS"))
        with open(b + ".fq", "a") as fh:
            fh.write("\n")
        r = subprocess.run([sys.executable, SIM_GENES, "compare", "--json", j, a, b], capture_output=True, text=True)
        self.assertEqual(r.returncode, 1)
        self.assertEqual(jload(j)["ok"], False)

    def test_a_shallower_arm_is_exactly_the_first_reads_of_the_deeper_one(self):
        small = self.run_sim("depth4", reads_per_tx=4)
        large = self.run_sim("depth10", reads_per_tx=10)
        ok, lines = self.sg.gate_depth_subset(small, large)
        self.assertTrue(ok, lines)
        self.assertTrue(lines[0].startswith("G-depth PASS"), lines)
        # and the reads are what a direct 4-per-transcript simulation gives: the same file as the first four reads of each transcript
        big = {n: s for n, s, _, _ in fastq_records(large + ".fq")}
        for n, s, _, _ in fastq_records(small + ".fq"):
            self.assertEqual(big[n], s)

    def test_the_depth_gate_fails_when_a_shallow_read_differs_or_is_missing_from_the_deep_set(self):
        small = self.fresh("depth_tamper_small")
        large = self.run_sim("depth_tamper_large", reads_per_tx=10)
        text = slurp(small + ".fq")
        first = text.splitlines()[1]
        with open(small + ".fq", "w") as fh:
            fh.write(text.replace(first, first[::-1], 1))
        ok, lines = self.sg.gate_depth_subset(small, large)
        self.assertFalse(ok)
        self.assertTrue(lines[0].startswith("G-depth FAIL"), lines)

    def test_the_depth_gate_fails_when_the_tables_disagree(self):
        small = self.run_sim("depth_tab_small", reads_per_tx=4)
        large = self.run_sim("depth_tab_large", reads_per_tx=10)
        rows = tsv(large + ".transcripts.tsv")
        rows[0]["chain"] = "1-2"
        with open(large + ".transcripts.tsv", "w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        ok, lines = self.sg.gate_depth_subset(small, large)
        self.assertFalse(ok)

    def test_cli_depth_subset(self):
        small = self.run_sim("depth_cli_small", reads_per_tx=4)
        large = self.run_sim("depth_cli_large", reads_per_tx=10)
        j = os.path.join(self.tmp, "depth.json")
        r = subprocess.run([sys.executable, SIM_GENES, "depth-subset", "--json", j, small, large], capture_output=True, text=True)
        self.assertEqual((r.returncode, r.stdout.startswith("G-depth PASS")), (0, True), r.stdout + r.stderr)
        self.assertEqual(jload(j)["ok"], True)

    def test_control_writes_the_annotation_as_loci_and_as_one_node_per_gene(self):
        ctrl = os.path.join(self.tmp, "ctrl0")
        r = subprocess.run([sys.executable, SIM_GENES, "control", "--out", self.out, "--to", ctrl], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        gtf = [ln.split("\t") for ln in slurp(ctrl + ".gtf").splitlines()]
        tx = [f for f in gtf if f[2] == "transcript"]
        self.assertEqual(len(tx), 10)    # every eligible transcript, whatever the arm
        self.assertIn('gene_id "TG_A"; transcript_id "TG_A.1"', tx[0][8])
        nodes = [ln.split("\t") for ln in slurp(ctrl + ".loci.gff3").splitlines() if ln and not ln.startswith("#")]
        genes = {f[8].split(";")[1][5:]: f for f in nodes if f[2] == "gene"}
        self.assertEqual(sorted(genes), ["TG_A", "TG_B", "TG_BIG", "TG_N1", "TG_NC", "TG_P", "TG_S130", "TG_SH"])
        exons = collections.defaultdict(list)
        for f in nodes:
            if f[2] == "exon":
                exons[f[8].split("gene=")[1]].append((int(f[3]) - 1, int(f[4])))
        self.assertEqual(sorted(exons["TG_A"]), [(1000, 1400), (2000, 2300), (3000, 3500)])   # the longest isoform represents the gene
        self.assertEqual(sorted(exons["TG_BIG"]), [(400000, 400400)])

    def test_control_can_shift_every_coordinate(self):
        ctrl = os.path.join(self.tmp, "ctrl50")
        r = subprocess.run([sys.executable, SIM_GENES, "control", "--out", self.out, "--to", ctrl, "--shift", "50000"], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        nodes = [ln.split("\t") for ln in slurp(ctrl + ".loci.gff3").splitlines() if ln and not ln.startswith("#")]
        a_exons = sorted((int(f[3]) - 1, int(f[4])) for f in nodes if f[2] == "exon" and f[8].endswith("gene=TG_A"))
        self.assertEqual(a_exons, [(51000, 51400), (52000, 52300), (53000, 53500)])
        tx = [ln.split("\t") for ln in slurp(ctrl + ".gtf").splitlines() if ln.split("\t")[2] == "transcript"]
        self.assertTrue(any(f[3] == "51001" for f in tx))

    def test_cli_modes_print_pass_or_fail_and_set_the_exit_status(self):
        r = subprocess.run([sys.executable, SIM_GENES, "gates", "--out", self.out, "--genome", self.world["fa"], "--genes", self.world["genes"]], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("G1 PASS", r.stdout)
        self.assertIn("G2 PASS", r.stdout)
        a, b = self.fresh("cli_a"), self.fresh("cli_b")
        r = subprocess.run([sys.executable, SIM_GENES, "compare", a, b], capture_output=True, text=True)
        self.assertEqual((r.returncode, r.stdout.startswith("G3 PASS")), (0, True), r.stdout + r.stderr)
        with open(b + ".fq", "a") as fh:
            fh.write("\n")
        r = subprocess.run([sys.executable, SIM_GENES, "compare", a, b], capture_output=True, text=True)
        self.assertEqual((r.returncode, r.stdout.startswith("G3 FAIL")), (1, True), r.stdout + r.stderr)
        d = self.fresh("cli_bad")
        with open(d + ".fq", "w") as fh:
            fh.write("")
        r = subprocess.run([sys.executable, SIM_GENES, "gates", "--out", d, "--genome", self.world["fa"], "--genes", self.world["genes"]], capture_output=True, text=True)
        self.assertEqual(r.returncode, 1)
        self.assertIn("G1 FAIL", r.stdout)


def write_bam(path, world, fq_prefix, place, secondary=()):
    """BAM of the reads of a simulation: place(name) -> None (unmapped) or (contig, pos0); `secondary` = [(name, contig, pos0)] extra secondary records."""
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": c, "LN": n} for c, n in CONTIG_LEN.items()]}
    ids = {c: i for i, c in enumerate(CONTIG_LEN)}
    seqs = {n: s for n, s, _, _ in fastq_records(fq_prefix + ".fq")}
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        def rec(name, flag, contig, pos):
            a = pysam.AlignedSegment(bam.header)
            a.query_name, a.flag, a.query_sequence = name, flag, seqs[name]
            if contig is None:
                a.is_unmapped = True
            else:
                a.reference_id, a.reference_start, a.mapping_quality, a.cigarstring = ids[contig], pos, 60, f"{len(seqs[name])}M"
            return a
        for name in seqs:
            p = place(name)
            bam.write(rec(name, 4 if p is None else 0, *(p or (None, 0))))
        for name, contig, pos in secondary:
            bam.write(rec(name, 256, contig, pos))


class BamQcTests(WorldCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.sg = __import__("sim_genes")
        cls.out = os.path.join(cls.tmp, "q")
        with contextlib.redirect_stdout(io.StringIO()):
            cls.sg.main(sim_args(cls.world, cls.out))
        cls.spans = {}   # gene -> (contig, first simulated exon start): where a correctly placed read of the gene lands
        for r in tsv(cls.out + ".transcripts.tsv"):
            if r["simulated"] == "1":
                start = int(r["canon_blocks"].split(",")[0].split("-")[0])
                if r["gene"] not in cls.spans or start < cls.spans[r["gene"]][1]:
                    cls.spans[r["gene"]] = (r["chrom"], start)

    def place(self, name):
        g = name.split("|")[1]
        if name.startswith("TG_B.1|TG_B|fl|0") or name.startswith("TG_B.1|TG_B|fl|1"):
            return None                              # two unmapped reads
        if g == "TG_P":
            return ("chr1", 11000)                   # the paralog's reads land on another contig at the SAME coordinates as the source exon (chr2:11000)
        return self.spans[g]

    def test_g4_counts_records_and_prints_the_two_shares(self):
        bam = os.path.join(self.tmp, "q.bam")
        write_bam(bam, self.world, self.out, self.place, secondary=[("TG_A.1|TG_A|fl|0", "chr2", 11000), ("TG_P.1|TG_P|fl|0", "chr1", 1000), ("TG_P.1|TG_P|fl|1", "chr1", 1000)])
        ok, lines, st = self.sg.gate_bamqc(self.out, bam)
        self.assertTrue(ok, lines)
        self.assertTrue(lines[0].startswith("G4 PASS"), lines)
        self.assertEqual((st["reads"], st["primary_or_unmapped_records"], st["mapped"], st["unmapped"], st["secondary_records"]), (40, 40, 38, 2, 3))
        self.assertEqual(st["primary_on_source_gene"], 30)   # 38 mapped, minus the 8 reads of TG_P that landed on TG_A
        self.assertEqual(st["soto"]["reads"], 32)
        self.assertEqual(st["neighbour"]["reads"], 8)
        self.assertEqual(st["neighbour"]["primary_on_source_gene"], 8)

    def test_g4_fails_when_a_read_has_no_record(self):
        bam = os.path.join(self.tmp, "q_short.bam")
        recs = [r for r in fastq_records(self.out + ".fq") if not r[0].startswith("TG_SH.1|TG_SH|fl|3")]
        short = os.path.join(self.tmp, "q_short_prefix")
        copy_sim(self.out, short)
        with open(short + ".fq", "w") as fh:
            for n, s, plus, q in recs:
                fh.write(f"@{n}\n{s}\n{plus}{q}\n")
        write_bam(bam, self.world, short, self.place)
        ok, lines, _ = self.sg.gate_bamqc(self.out, bam)   # the FASTQ of `out` has one more read than the BAM
        self.assertFalse(ok)
        self.assertTrue(lines[0].startswith("G4 FAIL"), lines)

    def test_g4_fails_when_a_read_has_two_primary_records(self):
        bam = os.path.join(self.tmp, "q_dup.bam")
        header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": c, "LN": n} for c, n in CONTIG_LEN.items()]}
        seqs = {n: s for n, s, _, _ in fastq_records(self.out + ".fq")}
        with pysam.AlignmentFile(bam, "wb", header=header) as fh:
            for i, name in enumerate(list(seqs) + ["TG_A.1|TG_A|fl|0"]):
                a = pysam.AlignedSegment(fh.header)
                a.query_name, a.flag, a.query_sequence = name, 0, seqs[name]
                a.reference_id, a.reference_start, a.mapping_quality, a.cigarstring = 0, 1000, 60, f"{len(seqs[name])}M"
                fh.write(a)
        ok, lines, _ = self.sg.gate_bamqc(self.out, bam)
        self.assertFalse(ok)
        self.assertTrue(lines[0].startswith("G4 FAIL"), lines)

    def test_cli_bamqc(self):
        bam = os.path.join(self.tmp, "q_cli.bam")
        write_bam(bam, self.world, self.out, self.place)
        r = subprocess.run([sys.executable, SIM_GENES, "bamqc", "--out", self.out, "--bam", bam], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("G4 PASS", r.stdout)
        self.assertTrue(os.path.exists(self.out + ".g4.json"))


# ------------------------------------------------------------------------------------------------- score_soto_nodes
def write_truth(root, name, genes):
    """Truth tables as sim_genes.py writes them, from a spec {gene: dict(chrom, strand, blocks | isoforms, biotype, soto, fam, half, name, n_reads (per transcript, default 10), simulable (default True),
    sd98, reason)}. Eligible transcripts are `<gene>.t<i>`; a gene with simulable=False has none."""
    sg = __import__("sim_genes")
    prefix = os.path.join(root, name)
    with open(prefix + ".transcripts.tsv", "w") as ft, open(prefix + ".targets.tsv", "w") as fg:
        ft.write("\t".join(sg.TX_COLUMNS) + "\n")
        fg.write("\t".join(sg.TARGET_COLUMNS) + "\n")
        for g, d in genes.items():
            isos = d.get("isoforms") or [d["blocks"]]
            simulable = d.get("simulable", True)
            nr = d.get("n_reads", 10)
            union = sorted({b for iso in isos for b in iso}) if simulable else []
            merged = []
            for s, e in union:
                if merged and s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            if simulable:
                for i, iso in enumerate(isos, 1):
                    cb = ",".join(f"{x}-{y}" for x, y in iso)
                    chain = ",".join(f"{iso[k][1]}-{iso[k + 1][0]}" for k in range(len(iso) - 1))
                    row = dict(transcript=f"{g}.t{i}", gene=g, gene_name=d.get("name", g), biotype=d["biotype"], is_soto=int(d["soto"]), chrom=d["chrom"], strand=d["strand"], raw_blocks=cb,
                               canon_blocks=cb, chain=chain, n_raw_introns=len(iso) - 1, n_noncanonical=0, n_short=0, spliced_len=sum(y - x for x, y in iso), seq_sha1="x" * 40, eligible=1,
                               simulated=int(nr > 0), skip_reason="", n_reads=nr)
                    ft.write("\t".join(str(row[c]) for c in sg.TX_COLUMNS) + "\n")
            if d["soto"]:
                fam = d.get("fam", ["ID_1"])
                single = int(all(len(iso) == 1 for iso in isos)) if simulable else ""
                row = dict(gene=g, gene_name=d.get("name", g), s1c_name=d.get("name", g), biotype=d["biotype"], s1c_biotype=d["biotype"], s1c_families=";".join(fam), n_families=len(fam),
                           half=d.get("half", "dev"), source=sg.source_of(g), sd98=d.get("sd98", "Yes"), chrom=d["chrom"] if simulable else "", strand=d["strand"] if simulable else "",
                           start0=merged[0][0] if merged else "", end=merged[-1][1] if merged else "", n_tx_annotated=len(isos) if simulable else 0, n_tx_eligible=len(isos) if simulable else 0,
                           n_tx_simulated=len(isos) if simulable and nr > 0 else 0, n_reads=nr * len(isos) if simulable else 0, exon_bp=sum(y - x for x, y in merged), n_exons=len(merged),
                           single_exon=single, simulable=int(simulable), identical_seq_genes=0,
                           status="simulated" if simulable and nr > 0 else "skipped", reason=("" if simulable and nr > 0 else d.get("reason", "weight 0" if simulable else "no simulable transcript (under 120 bp: 1)")))
                fg.write("\t".join(str(row[c]) for c in sg.TARGET_COLUMNS) + "\n")
    return prefix


def write_nodes(path, nodes):
    """The families stage's PREFIX.fam.loci.gff3: a gene line and the exon lines of the representative, per node. nodes: {id: dict(chrom, strand, exons=[(start0, end)])}."""
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for nid, d in nodes.items():
            ex = sorted(d["exons"])
            fh.write(f'{d["chrom"]}\t.\tgene\t{ex[0][0] + 1}\t{ex[-1][1]}\t.\t{d["strand"]}\t.\tID=gene-{nid};Name={nid}\n')
            for a, b in ex:
                fh.write(f'{d["chrom"]}\t.\texon\t{a + 1}\t{b}\t.\t{d["strand"]}\t.\tParent=gene-{nid};gene={nid}\n')


def node(chrom, strand, exons):
    return dict(chrom=chrom, strand=strand, exons=exons)


def write_gtf(path, units):
    """units: [dict(chrom, strand, gene, tid, exons=[(start0, end)])] as the assembler's GTF (1-based closed coordinates; a transcript line, then its exons)."""
    with open(path, "w") as fh:
        for u in units:
            ex = sorted(u["exons"])
            fh.write(f'{u["chrom"]}\trustle\ttranscript\t{ex[0][0] + 1}\t{ex[-1][1]}\t.\t{u["strand"]}\t.\tgene_id "{u["gene"]}"; transcript_id "{u["tid"]}"; reads "5";\n')
            for i, (a, b) in enumerate(ex, 1):
                fh.write(f'{u["chrom"]}\trustle\texon\t{a + 1}\t{b}\t.\t{u["strand"]}\t.\tgene_id "{u["gene"]}"; transcript_id "{u["tid"]}"; exon_number "{i}";\n')


def unit(chrom, strand, gene, tid, exons):
    return dict(chrom=chrom, strand=strand, gene=gene, tid=tid, exons=exons)


def write_read_bam(path, recs):
    """recs: [(read name, flag, contig, pos0, length)]; flag 4 = unmapped record. A tiny BAM for the scorer's read statistics."""
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": c, "LN": n} for c, n in CONTIG_LEN.items()]}
    ids = {c: i for i, c in enumerate(CONTIG_LEN)}
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for name, flag, contig, pos, ln in recs:
            a = pysam.AlignedSegment(bam.header)
            a.query_name, a.flag, a.query_sequence = name, flag, "A" * ln
            if flag & 4:
                a.is_unmapped = True
            else:
                a.reference_id, a.reference_start, a.mapping_quality, a.cigarstring = ids[contig], pos, 60, f"{ln}M"
            bam.write(a)


def run_scorer(truth, out, *extra, hashseed=None):
    env = dict(os.environ)
    if hashseed is not None:
        env["PYTHONHASHSEED"] = str(hashseed)
    cut = any(x in extra for x in ("--single-cut-file", "--freeze-single-cut", "--no-single-cut"))
    args = [sys.executable, SCORER, "--truth", truth, "--out", out, *extra] + ([] if cut else ["--no-single-cut"])
    return subprocess.run(args, capture_output=True, text=True, env=env)


class ScorerBase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="soto_score_test_")
        cls.sn = __import__("score_soto_nodes")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def run_nodes(self, name, genes, nodes, *extra, ok=0):
        truth = write_truth(self.tmp, name, genes)
        loci = os.path.join(self.tmp, name + ".loci.gff3")
        write_nodes(loci, nodes)
        out = os.path.join(self.tmp, name + ".score")
        r = run_scorer(truth, out, "--unit", "node", "--loci-gff3", loci, *extra)
        self.assertEqual(r.returncode, ok, r.stdout + r.stderr)
        return r, out

    def run_gtf(self, name, genes, units, *extra, ok=0):
        truth = write_truth(self.tmp, name, genes)
        gtf = os.path.join(self.tmp, name + ".gtf")
        write_gtf(gtf, units)
        out = os.path.join(self.tmp, name + ".score")
        r = run_scorer(truth, out, "--gtf", gtf, *(extra if "--unit" in extra else ("--unit", "transcript", *extra)))
        self.assertEqual(r.returncode, ok, r.stdout + r.stderr)
        return r, out

    def result(self, r, out):
        return {x["gene"]: x for x in tsv(out + ".per_gene.tsv")}, jload(out + ".summary.json"), r.stdout


def soto(chrom, strand, blocks, biotype="lncRNA", **kw):
    return dict(chrom=chrom, strand=strand, blocks=blocks, biotype=biotype, soto=True, **kw)


def other(chrom, strand, blocks, biotype="lncRNA", **kw):
    return dict(chrom=chrom, strand=strand, blocks=blocks, biotype=biotype, soto=False, **kw)


class ScorerClassTests(ScorerBase):
    # --- a node that is right, and the sign conventions of the distances
    def test_a_node_with_the_gene_s_chain_is_node_ok_and_distances_follow_the_gene_strand(self):
        genes = {"P": soto("chr1", "+", [(1000, 1200), (2000, 2300)], "protein_coding"), "T": soto("chr1", "+", [(10000, 10400)]), "M": soto("chr1", "-", [(20000, 20400)])}
        nodes = {"N1": node("chr1", "+", [(1000, 1200), (2000, 2300)]), "N2": node("chr1", "+", [(10100, 10400)]), "N3": node("chr1", "-", [(20050, 20400)])}
        per, summ, _ = self.result(*self.run_nodes("basic", genes, nodes))
        self.assertEqual({g: per[g]["class"] for g in per}, {"P": "node_ok", "T": "node_ok", "M": "node_ok"})
        self.assertEqual((per["P"]["unit"], per["P"]["jaccard"], per["P"]["d5"], per["P"]["d3"], per["P"]["width_ratio"]), ("N1", "1.0000", "0", "0", "1.0000"))
        self.assertEqual((per["T"]["d5"], per["T"]["d3"], per["T"]["jaccard"], per["T"]["width_ratio"]), ("100", "0", "0.7500", "0.7500"))
        self.assertEqual((per["M"]["d5"], per["M"]["d3"]), ("0", "50"))
        self.assertEqual(summ["classes"]["node_ok"], 3)

    # --- a 1-bp overlap is not a match: two-sided, at least half of the SMALLER exon union
    def test_the_overlap_must_be_at_least_half_of_the_smaller_exon_union(self):
        genes = {"IN": soto("chr1", "+", [(1000, 2000)]), "TINY": soto("chr1", "+", [(5000, 6000)]), "HALF": soto("chr1", "+", [(9000, 10000)]), "JUST": soto("chr1", "+", [(13000, 14000)])}
        nodes = {"A": node("chr1", "+", [(1900, 2000)]),      # inside the gene: 100 of 100 bp of the smaller union (the gene is 10 times longer)
                 "B": node("chr1", "+", [(5950, 7000)]),      # 50 bp of a 1000 bp gene and a 1050 bp node: weak
                 "C": node("chr1", "+", [(9500, 10500)]),     # exactly 500 of 1000
                 "D": node("chr1", "+", [(13501, 14501)])}    # 499 of 1000
        per, summ, _ = self.result(*self.run_nodes("twosided", genes, nodes))
        self.assertEqual({g: per[g]["class"] for g in per}, {"IN": "node_ok", "TINY": "weak_overlap", "HALF": "node_ok", "JUST": "weak_overlap"})
        self.assertEqual((per["TINY"]["matched"], per["TINY"]["overlap_bp"]), ("1", "50"))
        self.assertEqual(summ["classes"]["weak_overlap"], 2)

    def test_a_gene_without_an_overlapping_node_is_no_locus_and_the_split_by_exon_count_is_kept(self):
        genes = {"A": soto("chr1", "+", [(1000, 1400)]), "B": soto("chr1", "+", [(5000, 5400)]), "C": soto("chr1", "+", [(8000, 8200), (8600, 8900)])}
        nodes = {"N1": node("chr1", "+", [(1000, 1400)]), "N9": node("chr2", "+", [(5000, 5400)]), "N8": node("chr1", "+", [(5400, 6000)])}   # B: other contig and abutting only
        per, summ, _ = self.result(*self.run_nodes("noloc", genes, nodes))
        self.assertEqual({g: per[g]["class"] for g in per}, {"A": "node_ok", "B": "no_locus", "C": "no_locus"})
        self.assertEqual((per["B"]["n_candidates"], per["B"]["matched"]), ("0", "0"))
        self.assertEqual(summ["classes_by_stratum"]["single_exon"]["no_locus"], 1)
        self.assertEqual(summ["classes_by_stratum"]["multi_exon"]["no_locus"], 1)

    def test_the_loser_of_the_one_to_one_competition_is_outcompeted(self):
        genes = {"F1": soto("chr1", "+", [(1000, 1600)]), "F2": soto("chr1", "+", [(2000, 2400)])}
        nodes = {"L": node("chr1", "+", [(1000, 1600), (2000, 2400)])}
        per, summ, _ = self.result(*self.run_nodes("outc", genes, nodes))
        self.assertEqual((per["F1"]["matched"], per["F2"]["matched"], per["F2"]["n_candidates"], per["F2"]["class"]), ("1", "0", "1", "outcompeted"))

    # --- fusion is the EXCESS over the gene's own annotation
    def test_a_node_that_covers_a_gene_the_gene_s_own_annotation_does_not_cover_is_fused_excess(self):
        genes = {"F1": soto("chr1", "+", [(1000, 1600)]), "F2": soto("chr1", "+", [(2000, 2400)])}
        nodes = {"L": node("chr1", "+", [(1000, 1600), (2000, 2400)])}
        per, summ, _ = self.result(*self.run_nodes("excess", genes, nodes))
        self.assertEqual((per["F1"]["class"], per["F1"]["cov_span"], per["F1"]["own_cov_span"], per["F1"]["excess_span"], per["F1"]["excess_exon"]), ("fused_excess", "1", "0", "1", "1"))
        self.assertEqual(summ["fusion"]["excess_span"], 1)

    def test_a_gene_nested_in_the_intron_of_the_matched_gene_is_the_annotation_s_own_overlap_not_excess(self):
        genes = {"H": soto("chr1", "+", [(1000, 1300), (5000, 5300)]), "G": other("chr1", "+", [(2000, 2400)])}
        nodes = {"NH": node("chr1", "+", [(1000, 1300), (5000, 5300)]), "NG": node("chr1", "+", [(2000, 2400)])}
        per, summ, _ = self.result(*self.run_nodes("nested", genes, nodes))
        self.assertEqual((per["H"]["cov_span"], per["H"]["own_cov_span"], per["H"]["excess_span"], per["H"]["class"]), ("1", "1", "0", "node_ok"))
        self.assertEqual(per["H"]["n_genes_span"], "2")          # the raw call (the baseline's) says fused; the excess is zero
        self.assertEqual(summ["fusion"]["raw_span"], 1)
        self.assertEqual(summ["fusion"]["excess_span"], 0)
        # the same node with an extra exon over the nested gene: now it covers exons the gene's own annotation does not
        nodes2 = {"NH": node("chr1", "+", [(1000, 1300), (2000, 2400), (5000, 5300)])}
        per, summ, _ = self.result(*self.run_nodes("nested2", genes, nodes2, "--fused-rule", "exon"))
        self.assertEqual((per["H"]["excess_exon"], per["H"]["class"]), ("1", "fused_excess"))
        self.assertEqual(per["H"]["class_span"], "extent_wrong")      # under the span rule the nested gene was already inside the gene's own span: only the chain is wrong

    def test_a_node_covers_a_gene_by_exons_when_it_overlaps_at_least_half_of_the_gene_s_exon_bp(self):
        genes = {"S": soto("chr1", "+", [(1000, 1200)]), "N": other("chr1", "+", [(1500, 1700)])}
        half = {"U": node("chr1", "+", [(1000, 1200), (1600, 1700)])}      # 100 of N's 200 exon bp: exactly half
        per, _, _ = self.result(*self.run_nodes("exon_half", genes, half))
        self.assertEqual((per["S"]["cov_exon"], per["S"]["excess_exon"]), ("1", "1"))
        less = {"U": node("chr1", "+", [(1000, 1200), (1601, 1700)])}      # 99 of 200
        per, _, _ = self.result(*self.run_nodes("exon_less", genes, less))
        self.assertEqual((per["S"]["cov_exon"], per["S"]["excess_exon"]), ("0", "0"))

    def test_both_fusion_rules_are_reported_in_every_run_and_the_class_follows_the_chosen_one(self):
        genes = {"H": soto("chr1", "+", [(1000, 1300), (5000, 5300)], "protein_coding"), "G": other("chr1", "+", [(2000, 2400)])}
        nodes = {"NH": node("chr1", "+", [(1000, 1300), (2000, 2400), (5000, 5300)])}
        r, out = self.run_nodes("rules", genes, nodes)
        per, summ, text = self.result(r, out)
        self.assertEqual(sorted(summ["classes_by_rule"]), ["exon", "span"])
        self.assertEqual(summ["classes_by_rule"]["exon"]["fused_excess"], 1)
        self.assertEqual(summ["classes_by_rule"]["span"]["extent_wrong"], 1)
        self.assertEqual(summ["classes"], summ["classes_by_rule"]["span"])   # the default class rule is the span rule
        self.assertIn("fused-rule=span", text)

    # --- extent: the chain of an annotated isoform; for single-exon genes a frozen cut
    def test_the_node_chain_must_equal_the_chain_of_an_annotated_isoform_of_the_gene(self):
        iso1 = [(1000, 1200), (2000, 2300), (3000, 3400)]
        iso2 = [(1000, 1200), (3000, 3400)]
        genes = {"X": soto("chr1", "+", None, "protein_coding", isoforms=[iso1, iso2]), "Y": soto("chr1", "+", None, "protein_coding", isoforms=[[(10000, 10300), (11000, 11300)]]),
                 "Z": soto("chr1", "+", None, "protein_coding", isoforms=[[(20000, 20300), (21000, 21300)]])}
        nodes = {"NX": node("chr1", "+", iso2),                                          # the second isoform: an annotated chain
                 "NY": node("chr1", "+", [(10000, 10300), (10500, 10700), (11000, 11300)]),   # an exon the gene does not have
                 "NZ": node("chr1", "+", [(20000, 21300)])}                              # unspliced over a spliced gene
        per, summ, _ = self.result(*self.run_nodes("chains", genes, nodes))
        self.assertEqual({g: per[g]["class"] for g in per}, {"X": "node_ok", "Y": "extent_wrong", "Z": "extent_wrong"})
        self.assertEqual((per["X"]["chain_equal"], per["Y"]["chain_equal"], per["Z"]["chain_equal"]), ("1", "0", "0"))

    def test_a_single_exon_gene_matched_by_a_spliced_node_has_the_wrong_extent_without_any_cut(self):
        genes = {"S": soto("chr1", "+", [(1000, 2000)])}
        nodes = {"N": node("chr1", "+", [(1000, 1400), (1600, 2000)])}
        per, _, _ = self.result(*self.run_nodes("spliced_single", genes, nodes))
        self.assertEqual(per["S"]["class"], "extent_wrong")

    def test_the_ends_are_not_compared_when_the_chain_is_right(self):
        genes = {"X": soto("chr1", "+", [(1000, 1200), (2000, 2300)], "protein_coding")}
        nodes = {"N": node("chr1", "+", [(900, 1200), (2000, 2400)])}   # both outer ends moved by 100 bp: the chain is the same
        per, _, _ = self.result(*self.run_nodes("ends", genes, nodes))
        self.assertEqual(per["X"]["class"], "node_ok")

    def _cut_world(self, jacs_dev, jacs_held):
        genes, nodes = {}, {}
        i = 0
        for half, jacs in (("dev", jacs_dev), ("heldout", jacs_held)):
            for j in jacs:
                s = 1000 + 3000 * i
                genes[f"G{i:03d}"] = soto("chr1", "+", [(s, s + 1000)], half=half)
                nodes[f"N{i:03d}"] = node("chr1", "+", [(s, s + int(round(j * 1000)))])
                i += 1
        return genes, nodes

    def test_the_single_exon_cut_is_frozen_from_the_dev_half_of_the_run_that_is_told_to_freeze_it(self):
        dev = [0.2, 0.3, 0.8] + [0.8 + 0.005 * k for k in range(1, 39)]   # 41 values: the 5th percentile is the third smallest, 0.8
        genes, nodes = self._cut_world(dev, [0.5, 0.9])
        cutfile = os.path.join(self.tmp, "cut1.json")
        r, out = self.run_nodes("freeze", genes, nodes, "--freeze-single-cut", cutfile)
        per, summ, text = self.result(r, out)
        cut = jload(cutfile)
        self.assertAlmostEqual(cut["value"], 0.8, places=9)
        self.assertEqual((cut["percentile"], cut["n"], cut["half"]), (5, 41, "dev"))
        self.assertIn("lowest 5% of the", cut["rule"])
        bad = {g for g in per if per[g]["class"] == "extent_wrong"}
        self.assertEqual(bad, {"G000", "G001", "G041"})   # strictly below 0.8 on the dev half, and the held-out gene at 0.5 (frozen cut applied as is)
        self.assertEqual(per["G042"]["class"], "node_ok")
        self.assertIn("FROZE", text)
        self.assertEqual(summ["single_cut"]["value"], 0.8)

    def test_a_frozen_cut_is_not_recomputed_in_another_arm(self):
        genes, nodes = self._cut_world([0.2, 0.3, 0.8] + [0.8 + 0.005 * k for k in range(1, 39)], [])
        cutfile = os.path.join(self.tmp, "cut2.json")
        self.run_nodes("freeze2", genes, nodes, "--freeze-single-cut", cutfile)
        genes2, nodes2 = self._cut_world([0.5] * 10, [0.5] * 10)     # a worse arm: every gene at 0.5
        r, out = self.run_nodes("arm2", genes2, nodes2, "--single-cut-file", cutfile)
        per, summ, text = self.result(r, out)
        self.assertEqual(summ["classes"]["extent_wrong"], 20)          # all of them: no per-arm tail
        self.assertEqual(summ["single_cut"]["value"], 0.8)
        self.assertIn("frozen", text)

    def test_the_freeze_refuses_to_overwrite_and_a_cut_option_is_required(self):
        genes, nodes = self._cut_world([0.5, 0.6, 0.7, 0.9], [])
        cutfile = os.path.join(self.tmp, "cut3.json")
        self.run_nodes("freeze3", genes, nodes, "--freeze-single-cut", cutfile)
        r, _ = self.run_nodes("freeze3b", genes, nodes, "--freeze-single-cut", cutfile, ok=2)
        self.assertIn("already exists", r.stderr)
        truth = write_truth(self.tmp, "nocut", genes)
        loci = os.path.join(self.tmp, "nocut.loci.gff3")
        write_nodes(loci, nodes)
        r = subprocess.run([sys.executable, SCORER, "--truth", truth, "--out", os.path.join(self.tmp, "nocut.score"), "--unit", "node", "--loci-gff3", loci], capture_output=True, text=True)
        self.assertEqual(r.returncode, 2)
        self.assertIn("--single-cut-file", r.stderr)

    def test_no_single_cut_leaves_the_extent_of_single_exon_genes_unjudged_and_says_so(self):
        genes, nodes = self._cut_world([0.1, 0.9], [])
        r, out = self.run_nodes("nocut2", genes, nodes)
        per, summ, text = self.result(r, out)
        self.assertEqual({per[g]["class"] for g in per}, {"node_ok"})
        self.assertEqual(summ["single_cut"], None)
        self.assertIn("not judged", text)

    # --- reads
    def test_a_gene_with_no_reads_in_the_arm_is_no_reads_even_if_a_node_sits_on_it_and_that_node_is_counted_apart(self):
        genes = {"Z": soto("chr1", "+", [(1000, 1400)], n_reads=0), "Y": soto("chr1", "+", [(5000, 5400)], n_reads=0), "OK": soto("chr1", "+", [(9000, 9400)])}
        nodes = {"N1": node("chr1", "+", [(1000, 1400)]), "N3": node("chr1", "+", [(9000, 9400)])}
        per, summ, _ = self.result(*self.run_nodes("noreads", genes, nodes))
        self.assertEqual({g: per[g]["class"] for g in per}, {"Z": "no_reads", "Y": "no_reads", "OK": "node_ok"})
        self.assertEqual(summ["no_reads"], {"genes": 2, "with_node": 1})
        self.assertEqual(summ["population"]["simulable"], 3)      # the zero-read genes stay in the denominator

    def test_reads_that_never_align_on_their_gene_make_it_reads_not_on_gene(self):
        genes = {"A": soto("chr1", "+", [(1000, 1400)]), "B": soto("chr1", "+", [(5000, 5400)]), "C": soto("chr1", "+", [(9000, 9400)])}
        nodes = {"N1": node("chr1", "+", [(1000, 1400)]), "N2": node("chr1", "+", [(5000, 5400)]), "N3": node("chr1", "+", [(9000, 9400)])}
        recs = []
        for k in range(4):
            recs.append((f"A.t1|A|fl|{k}", 0, "chr1", 1000, 300))
            recs.append((f"B.t1|B|fl|{k}", 0, "chr1", 9000, 300))      # B's reads all land on C
            recs.append((f"B.t1|B|fl|{k}", 256, "chr1", 5000, 300))    # with a secondary on B itself
            recs.append((f"C.t1|C|fl|{k}", 0 if k < 2 else 4, "chr1", 9000, 300))
        bam = os.path.join(self.tmp, "rnog.bam")
        write_read_bam(bam, recs)
        per, summ, text = self.result(*self.run_nodes("rnog", genes, nodes, "--bam", bam))
        self.assertEqual({g: per[g]["class"] for g in per}, {"A": "node_ok", "B": "reads_not_on_gene", "C": "node_ok"})
        self.assertEqual((per["A"]["reads_on_gene"], per["B"]["reads_on_gene"], per["C"]["reads_on_gene"]), ("4", "0", "2"))
        # B has a node although none of its reads has a primary alignment on it: a locus that exists through shared (secondary) alignments
        self.assertEqual(summ["recovered"], {"nodes_with_overlap": 3, "own_read": 2, "shared_only": 1, "unknown": 0})
        self.assertIn("shared", text)

    def test_without_a_bam_the_read_provenance_of_a_node_is_unknown_not_assumed(self):
        genes = {"A": soto("chr1", "+", [(1000, 1400)])}
        nodes = {"N1": node("chr1", "+", [(1000, 1400)])}
        per, summ, _ = self.result(*self.run_nodes("nobam", genes, nodes))
        self.assertEqual(summ["recovered"]["unknown"], 1)

    def test_a_real_library_has_no_simulated_reads_and_takes_its_read_counts_from_a_table(self):
        genes = {"A": soto("chr1", "+", [(1000, 1400)]), "B": soto("chr1", "+", [(5000, 5400)])}
        nodes = {"N1": node("chr1", "+", [(1000, 1400)])}
        reads = os.path.join(self.tmp, "reads.tsv")
        with open(reads, "w") as fh:
            fh.write("gene_id\tprim_on_exons\nA\t7\nB\t0\n")
        per, summ, text = self.result(*self.run_nodes("real", genes, nodes, "--real-library", "--reads-per-gene", reads, "--reads-column", "prim_on_exons"))
        self.assertEqual({g: per[g]["class"] for g in per}, {"A": "node_ok", "B": "no_reads"})
        self.assertEqual((per["A"]["n_reads"], per["B"]["n_reads"]), ("7", "0"))      # the real library's reads, not the ideal arm's 10 per transcript
        self.assertIn("arm=real", text)

    def test_a_real_library_without_a_read_table_does_not_judge_reads(self):
        genes = {"A": soto("chr1", "+", [(1000, 1400)], n_reads=0)}     # the ideal arm's n_reads (here 0) must not be read as the real library's
        nodes = {"N1": node("chr1", "+", [(1000, 1400)])}
        per, _, _ = self.result(*self.run_nodes("real2", genes, nodes, "--real-library"))
        self.assertEqual(per["A"]["class"], "node_ok")

    # --- population and universe
    def test_the_population_is_the_simulable_listed_genes_by_id_and_the_others_are_listed_apart(self):
        genes = {"A": soto("chr1", "+", [(1000, 1400)]), "TINY": soto("chr1", "+", [(1, 2)], "miRNA", simulable=False, sd98="Yes"),
                 "TINY2": soto("chr1", "+", [(1, 2)], "miRNA", simulable=False, sd98="No"), "N": other("chr1", "+", [(5000, 5400)])}
        nodes = {"N1": node("chr1", "+", [(1000, 1400)])}
        per, summ, text = self.result(*self.run_nodes("pop", genes, nodes))
        self.assertEqual(summ["population"]["listed"], 3)
        self.assertEqual(summ["population"]["simulable"], 1)
        self.assertEqual(summ["population"]["not_simulable"], {"genes": 2, "sd98_eligible": 1, "by_reason": {"no simulable transcript": 2}})
        self.assertEqual({g: per[g]["class"] for g in per}, {"A": "node_ok", "TINY": "not_simulable", "TINY2": "not_simulable"})
        self.assertNotIn("N", per)
        self.assertEqual(sum(summ["classes"].values()), 1)

    def test_the_gene_universe_comes_from_the_eligible_transcripts_not_from_the_arm_s_reads(self):
        genes = {"S": soto("chr1", "+", [(1000, 1200)]), "N": other("chr1", "+", [(1100, 1900)], n_reads=0)}   # a silent neighbour that overlaps more
        nodes = {"L": node("chr1", "+", [(1000, 1900)])}
        per, summ, _ = self.result(*self.run_nodes("universe", genes, nodes))
        self.assertEqual(per["S"]["class"], "outcompeted")           # the neighbour still takes the node: it is part of the ideal arm's universe
        self.assertEqual(summ["n_universe_genes"], 2)

    def test_competition_can_be_restricted_to_the_listed_genes(self):
        genes = {"S": soto("chr1", "+", [(1000, 1200)]), "N": other("chr1", "+", [(1100, 1900)])}
        nodes = {"L": node("chr1", "+", [(1000, 1900)])}
        per, summ, text = self.result(*self.run_nodes("compete_soto", genes, nodes, "--compete", "soto"))
        self.assertEqual(per["S"]["matched"], "1")
        self.assertIn("compete=soto", text)

    def test_a_subset_is_given_by_gene_id_and_a_name_column_is_refused(self):
        genes = {"A": soto("chr1", "+", [(1000, 1400)], name="NAMEA"), "B": soto("chr1", "+", [(5000, 5400)], name="NAMEB")}
        nodes = {"N1": node("chr1", "+", [(1000, 1400)])}
        ids = os.path.join(self.tmp, "ids.tsv")
        with open(ids, "w") as fh:
            fh.write("Gene ID\nA\nQ\n")
        per, summ, text = self.result(*self.run_nodes("subset", genes, nodes, "--subset", ids, "--subset-label", "ids"))
        self.assertEqual((summ["subsets"]["ids"]["n_ids"], summ["subsets"]["ids"]["n"]), (2, 1))
        self.assertEqual(summ["subsets"]["ids"]["classes"]["node_ok"], 1)
        names = os.path.join(self.tmp, "names.tsv")
        with open(names, "w") as fh:
            fh.write("family\tgene\nID_1\tNAMEA\n")
        r, _ = self.run_nodes("subset2", genes, nodes, "--subset", names, ok=2)
        self.assertIn("gene ids", r.stderr)

    # --- units
    def test_a_transcript_of_any_gtf_and_a_gene_id_locus_are_units_beside_the_node(self):
        genes = {"A": soto("chr1", "+", None, isoforms=[[(1000, 1200), (2000, 2300)]])}
        units = [unit("chr1", "+", "LOC", "T1", [(1000, 1200)]), unit("chr1", "+", "LOC", "T2", [(2000, 2300)])]
        per, summ, text = self.result(*self.run_gtf("unit_tx", genes, units))
        self.assertEqual((per["A"]["unit"], per["A"]["jaccard"], per["A"]["class"]), ("T2", "0.6000", "extent_wrong"))
        self.assertIn("unit=transcript", text)
        per, summ, text = self.result(*self.run_gtf("unit_loc", genes, units, "--unit", "locus"))
        self.assertEqual((per["A"]["unit"], per["A"]["jaccard"], per["A"]["class"]), ("LOC", "1.0000", "node_ok"))
        self.assertIn("unit=locus", text)

    def test_the_node_unit_needs_the_loci_file_and_the_other_units_need_a_gtf(self):
        genes = {"A": soto("chr1", "+", [(1000, 1400)])}
        truth = write_truth(self.tmp, "need", genes)
        r = run_scorer(truth, os.path.join(self.tmp, "need.score"), "--unit", "node")
        self.assertEqual(r.returncode, 2)
        self.assertIn("--loci-gff3", r.stderr)
        r = run_scorer(truth, os.path.join(self.tmp, "need2.score"), "--unit", "transcript")
        self.assertEqual(r.returncode, 2)
        self.assertIn("--gtf", r.stderr)

    # --- strata, tables, template
    def test_every_table_states_the_unit_the_competition_scope_and_the_fusion_rule(self):
        genes = {"A": soto("chr1", "+", [(1000, 1400)], "protein_coding")}
        nodes = {"N1": node("chr1", "+", [(1000, 1400)])}
        r, out = self.run_nodes("ctx", genes, nodes)
        text = r.stdout
        titles = [ln for ln in text.splitlines() if ln.startswith("TABLE")]
        self.assertGreaterEqual(len(titles), 5)
        for ln in titles:
            self.assertIn("unit=node", ln)
            self.assertIn("compete=all", ln)
            self.assertIn("fused-rule=span", ln)
            self.assertIn("match=two-sided>=0.5", ln)

    def test_classes_are_tabulated_by_biotype_source_half_and_exon_count(self):
        genes = {"A": soto("chr1", "+", [(1000, 1400)], "protein_coding", half="dev"), "CHM13_G1": soto("chr1", "+", [(5000, 5400)], "lncRNA", half="heldout"),
                 "LOFF_G2": soto("chr1", "+", [(9000, 9400), (9800, 10000)], "unprocessed_pseudogene", half="dev")}
        nodes = {"N1": node("chr1", "+", [(1000, 1400)]), "N3": node("chr1", "+", [(9000, 9400), (9800, 10000)])}
        per, summ, _ = self.result(*self.run_nodes("strata", genes, nodes))
        self.assertEqual(summ["classes_by_biotype"]["protein_coding"]["node_ok"], 1)
        self.assertEqual(summ["classes_by_biotype"]["lncRNA"]["no_locus"], 1)
        self.assertEqual(summ["classes_by_source"]["Liftoff"]["node_ok"], 1)
        self.assertEqual(summ["classes_by_source"]["CAT"]["no_locus"], 1)
        self.assertEqual(summ["classes_by_half"]["heldout"]["no_locus"], 1)
        self.assertEqual(summ["classes_by_stratum"]["multi_exon"]["node_ok"], 1)
        self.assertEqual(list(summ["classes_by_biotype"]), ["protein_coding", "transcribed_unprocessed_pseudogene", "unprocessed_pseudogene", "processed_pseudogene", "lncRNA", "other"])

    def test_descriptors_are_given_for_node_ok_genes_and_labelled_as_not_estimands(self):
        genes = {"A": soto("chr1", "+", [(1000, 2000)]), "B": soto("chr1", "+", [(5000, 6000)])}
        nodes = {"N1": node("chr1", "+", [(1000, 1800)]), "N2": node("chr1", "+", [(5000, 5500)])}
        _, summ, text = self.result(*self.run_nodes("desc", genes, nodes))
        self.assertAlmostEqual(summ["descriptors"]["node_ok"]["median_jaccard"], 0.65)
        self.assertIn("not estimands", text)

    def test_the_71_percent_appears_once_as_a_labelled_non_like_for_like_line_and_no_biotype_comparison_is_made_with_it(self):
        genes = {"A": soto("chr1", "+", [(1000, 1400)], "protein_coding")}
        _, summ, text = self.result(*self.run_nodes("ref", genes, {"N1": node("chr1", "+", [(1000, 1400)])}))
        self.assertEqual(text.count("71.0%"), 1)
        line = next(ln for ln in text.splitlines() if "71.0%" in ln)
        self.assertIn("NOT like-for-like", line)
        self.assertNotIn("79.3", text)
        self.assertNotIn("gain", text)
        self.assertEqual(summ["reference_line_not_like_for_like"]["matched_share"], 71.0)

    def test_the_report_prints_the_tie_policy_and_the_assembly_polish_counts(self):
        genes = {"A": soto("chr1", "+", [(1000, 1400)])}
        log = os.path.join(self.tmp, "assemble.log")
        with open(log, "w") as fh:
            fh.write("[copy_assign] x\n[copy_assign] ⭐ ASSEMBLY POLISH (full): 2593 transcripts -> ISM dropped 287 -> mono floor [0, 10, 9] (per contig) reads dropped 115 -> 2105 kept\nlast\n")
        _, summ, text = self.result(*self.run_nodes("polish", genes, {"N1": node("chr1", "+", [(1000, 1400)])}, "--assemble-log", log))
        self.assertIn("tie policy", text)
        self.assertIn("mono floor [0, 10, 9]", text)
        self.assertEqual(len(summ["assembly_polish"]), 1)

    def test_the_report_is_a_description_with_a_template_sentence_and_no_verdict(self):
        genes = {"A": soto("chr1", "+", [(1000, 1400)], "protein_coding")}
        _, summ, text = self.result(*self.run_nodes("verdict", genes, {"N1": node("chr1", "+", [(1000, 1400)])}))
        self.assertIn("INTERPRETATION TEMPLATE", text)
        for word in ("PASS", "FAIL", "VALID", "confirmed", "refuted", "improved", "worse", "gain"):
            self.assertNotIn(word, text)
        self.assertEqual(summ["interpretation_template"].count("\n"), 0)

    def test_comparators_are_read_from_their_summaries_and_printed_beside_the_arm(self):
        genes = {"A": soto("chr1", "+", [(1000, 1400)]), "B": soto("chr1", "+", [(5000, 5400)])}
        r1, out1 = self.run_nodes("cmp_a", genes, {"N1": node("chr1", "+", [(1000, 1400)]), "N2": node("chr1", "+", [(5000, 5400)])})
        r2, out2 = self.run_nodes("cmp_b", genes, {"N1": node("chr1", "+", [(1000, 1400)])}, "--compare", f"ideal={out1}.summary.json")
        text = r2.stdout
        self.assertIn("COMPARATORS", text)
        line = next(ln for ln in text.splitlines() if ln.startswith("ideal "))
        self.assertIn("unit=node", line)
        self.assertIn("100.0%", line)
        self.assertEqual(jload(out2 + ".summary.json")["comparators"][0]["label"], "ideal")

    # --- determinism
    def test_the_result_does_not_depend_on_the_order_of_the_inputs_or_on_tie_breaking(self):
        rng = random.Random(5)
        genes, nodes = {}, {}
        for i in range(30):
            s = 1000 + 2500 * i
            genes[f"G{i:02d}"] = soto("chr1", rng.choice("+-"), [(s, s + 800)], rng.choice(["protein_coding", "lncRNA"]))
            nodes[f"N{i:02d}"] = node("chr1", "+", [(s, s + 800)])
        genes["TA"] = soto("chr2", "+", [(1000, 1500)])
        genes["TB"] = soto("chr2", "+", [(1500, 2000)])
        nodes["NT"] = node("chr2", "+", [(1000, 2000)])           # an exact tie between two genes for one node
        r, out = self.run_nodes("order_a", genes, nodes)
        a, sa, _ = self.result(r, out)
        keys, nk = list(genes), list(nodes)
        rng.shuffle(keys); rng.shuffle(nk)
        r2, out2 = self.run_nodes("order_b", {k: genes[k] for k in keys}, {k: nodes[k] for k in nk})
        b, sb, _ = self.result(r2, out2)
        self.assertEqual(slurp(out + ".per_gene.tsv"), slurp(out2 + ".per_gene.tsv"))
        self.assertEqual({k: v for k, v in sa.items() if k != "inputs"}, {k: v for k, v in sb.items() if k != "inputs"})
        self.assertEqual((a["TA"]["class"], a["TB"]["class"]), ("fused_excess", "outcompeted"))   # the documented tie-break: the smaller gene id first
        self.assertEqual(sa["tie"]["sensitive_genes"], 2)
        self.assertIn("tie policy", r.stdout)

    def test_a_tie_between_nodes_of_identical_geometry_changes_no_outcome(self):
        genes = {"A": soto("chr1", "+", [(1000, 1400), (2000, 2300)])}
        nodes = {"N1": node("chr1", "+", [(1000, 1400), (2000, 2300)]), "N2": node("chr1", "+", [(1000, 1400), (2000, 2300)])}
        per, summ, _ = self.result(*self.run_nodes("harmless_tie", genes, nodes))
        self.assertEqual((per["A"]["matched"], summ["tie"]["sensitive_genes"]), ("1", 0))

    def test_tables_are_byte_identical_under_hash_seeds_0_and_1(self):
        genes = {f"G{i}": soto("chr1", "+", [(1000 + 3000 * i, 2000 + 3000 * i)], ("lncRNA", "protein_coding", "processed_pseudogene")[i % 3], fam=[f"ID_{i % 4}"]) for i in range(12)}
        nodes = {f"N{i}": node("chr1", "+", [(1000 + 3000 * i, 1900 + 3000 * i)]) for i in range(12)}
        truth = write_truth(self.tmp, "hs", genes)
        loci = os.path.join(self.tmp, "hs.loci.gff3")
        write_nodes(loci, nodes)
        outs = []
        for hs in (0, 1):
            out = os.path.join(self.tmp, f"hs{hs}.score")
            r = run_scorer(truth, out, "--unit", "node", "--loci-gff3", loci, hashseed=hs)
            self.assertEqual(r.returncode, 0, r.stderr)
            outs.append((out, r.stdout))
        for suf in (".per_gene.tsv", ".by_family.tsv", ".summary.json", ".report.txt"):
            self.assertEqual(slurp(outs[0][0] + suf), slurp(outs[1][0] + suf), suf)
        self.assertEqual(outs[0][1].replace(outs[0][0], "X"), outs[1][1].replace(outs[1][0], "X"))

    def test_family_table_counts_a_gene_in_every_family_it_belongs_to(self):
        genes = {"A1": soto("chr1", "+", [(1000, 1400)], fam=["ID_1"]), "A2": soto("chr1", "+", [(3000, 3400)], fam=["ID_1", "ID_2"]), "B1": soto("chr1", "+", [(5000, 5400)], fam=["ID_2"])}
        nodes = {"N1": node("chr1", "+", [(1000, 1400)]), "N2": node("chr1", "+", [(3000, 3400)])}
        r, out = self.run_nodes("fam", genes, nodes)
        fam = {x["family"]: x for x in tsv(out + ".by_family.tsv")}
        self.assertEqual((fam["ID_1"]["n_scored"], fam["ID_1"]["node_ok"], fam["ID_1"]["all_node_ok"]), ("2", "2", "1"))
        self.assertEqual((fam["ID_2"]["n_scored"], fam["ID_2"]["node_ok"], fam["ID_2"]["all_node_ok"]), ("2", "1", "0"))


class ScorerInputsTests(ScorerBase):
    def test_the_families_stage_nodes_are_read_from_the_loci_gff3_with_the_representative_exons(self):
        loci = os.path.join(self.tmp, "u.loci.gff3")
        write_nodes(loci, {"DN_chr1_100_1": node("chr1", "+", [(100, 200), (300, 400)]), "DN_chr2_5_1": node("chr2", "-", [(5, 50)])})
        u = self.sn.load_units(loci, "node", None)
        self.assertEqual(sorted(u), ["DN_chr1_100_1", "DN_chr2_5_1"])
        self.assertEqual((u["DN_chr1_100_1"]["exons"], u["DN_chr1_100_1"]["span"], u["DN_chr1_100_1"]["strand"]), ([(100, 200), (300, 400)], (100, 400), "+"))
        self.assertEqual(u["DN_chr2_5_1"]["chrom"], "chr2")

    def test_any_gtf_is_accepted_whatever_its_attributes(self):
        gtf = os.path.join(self.tmp, "any.gtf")
        with open(gtf, "w") as fh:
            fh.write('# comment\nchr1\tsomething\ttranscript\t101\t400\t.\t+\t.\ttranscript_id "t1"; gene_id "g1"; extra "x";\n'
                     'chr1\tsomething\texon\t101\t200\t.\t+\t.\ttranscript_id "t1"; gene_id "g1";\nchr1\tsomething\texon\t301\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
                     'chr1\tsomething\tCDS\t101\t200\t.\t+\t0\ttranscript_id "t1"; gene_id "g1";\n')
        u = self.sn.load_units(gtf, "transcript", None)
        self.assertEqual(u["t1"]["exons"], [(100, 200), (300, 400)])
        self.assertEqual(self.sn.load_units(gtf, "locus", None)["g1"]["exons"], [(100, 200), (300, 400)])

    def test_units_on_contigs_without_genes_are_dropped_on_load(self):
        gtf = os.path.join(self.tmp, "drop.gtf")
        write_gtf(gtf, [unit("chr1", "+", "a", "t1", [(10, 20)]), unit("chrZ", "+", "b", "t2", [(10, 20)])])
        self.assertEqual(sorted(self.sn.load_units(gtf, "transcript", {"chr1"})), ["t1"])

    def test_percentile_is_the_linear_interpolation_used_by_numpy(self):
        import statistics
        vals = [0.1, 0.5, 0.2, 0.9, 0.4, 0.3, 0.8]
        self.assertAlmostEqual(self.sn.percentile(vals, 5), statistics.quantiles(vals, n=20, method="inclusive")[0])
        self.assertEqual(self.sn.percentile([3.0], 5), 3.0)
        self.assertIsNone(self.sn.percentile([], 5))

    def test_the_reference_line_is_read_from_the_fidelity_document(self):
        b = self.sn.read_baseline(FIDELITY_DOC)
        self.assertEqual((b["n"], b["matched"], b["fused"]), (1078, 765, 258))
        self.assertAlmostEqual(b["matched_share"], 71.0)
        self.assertIsNone(self.sn.read_baseline(os.path.join(self.tmp, "absent.md")))


# ------------------------------------------------------------------------------------------------- equivalence with the original fidelity tool
GIT_REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
ORIGINAL_COMMIT = "8db314c7^"   # the last commit that still held bench/locus_probes.py (`soto-fidelity` = the former bench/soto_family_locus_fidelity.py)


class EquivalenceWithOriginalTests(unittest.TestCase):
    """The historical `locus_probes.py soto-fidelity` (git history) and score_soto_nodes.py make the same one-to-one assignments, d5, d3 and raw fused calls on a random fixture when the new
    scorer is asked for the old setting (transcripts as units, competition among the listed genes, the one-bp overlap rule is only a CLASS here, not a match condition): the new classes are
    built on the same assignment, so this is the scorer's link to the baseline tool."""

    @classmethod
    def setUpClass(cls):
        r = subprocess.run(["git", "-C", GIT_REPO, "show", f"{ORIGINAL_COMMIT}:bench/locus_probes.py"], capture_output=True, text=True)
        if r.returncode != 0:
            raise unittest.SkipTest("the historical bench/locus_probes.py is not in this checkout's git history")
        cls.tmp = tempfile.mkdtemp(prefix="soto_equiv_test_")
        cls.original = os.path.join(cls.tmp, "locus_probes.py")
        with open(cls.original, "w") as fh:
            fh.write(r.stdout)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_same_assignments_distances_and_raw_fusion_calls(self):
        rng = random.Random(11)
        genes, units, gff, soto_rows, last_unit = {}, [], ["##gff-version 3"], ["Gene Name\tFamily ID"], {}
        pos = {"chr1": 5000, "chr2": 5000}
        for i in range(300):
            chrom = "chr1" if i % 2 == 0 else "chr2"
            strand = rng.choice("+-")
            n_ex = rng.choice([1, 2, 3])
            start = pos[chrom] + rng.randint(300, 3000)
            blocks, cur = [], start
            for _ in range(n_ex):
                ln = rng.randint(120, 700)
                blocks.append((cur, cur + ln))
                cur += ln + rng.randint(80, 600)
            pos[chrom] = blocks[-1][1]
            listed = rng.random() < 0.7
            name = f"GENE{i:03d}"
            genes[name] = dict(chrom=chrom, strand=strand, blocks=blocks, biotype="lncRNA", soto=listed, name=name, fam=[f"ID_{i % 9}"])
            gff.append(f"{chrom}\tRefSeq\tgene\t{blocks[0][0] + 1}\t{blocks[-1][1]}\t.\t{strand}\t.\tID=gene-{name};Name={name}")
            gff.append(f"{chrom}\tRefSeq\tmRNA\t{blocks[0][0] + 1}\t{blocks[-1][1]}\t.\t{strand}\t.\tID=rna-{name};Parent=gene-{name}")
            for a, b in blocks:
                gff.append(f"{chrom}\tRefSeq\texon\t{a + 1}\t{b}\t.\t{strand}\t.\tID=exon-{name}-{a};Parent=rna-{name}")
            if listed:
                soto_rows.append(f"{name}\tID_{i % 9}")
            kind = rng.random()
            if kind < 0.15:
                continue
            if kind < 0.45 and chrom in last_unit:
                last_unit[chrom]["exons"] = sorted(last_unit[chrom]["exons"] + blocks)   # one transcript through this gene and the previous one on the contig
                continue
            ex = [(a + rng.randint(-40, 60), b + rng.randint(-60, 40)) for a, b in blocks]
            ex = [(a, b) for a, b in ex if b - a >= 60]
            if ex:
                last_unit[chrom] = unit(chrom, strand if rng.random() < 0.9 else "-", f"L{i}", f"T{i:03d}", ex)
                units.append(last_unit[chrom])
        gffp = os.path.join(self.tmp, "ref.gff3")
        with open(gffp, "w") as fh:
            fh.write("\n".join(gff) + "\n")
        sotop = os.path.join(self.tmp, "soto.tsv")
        with open(sotop, "w") as fh:
            fh.write("\n".join(soto_rows) + "\n")
        gtfp = os.path.join(self.tmp, "asm.gtf")
        write_gtf(gtfp, units)
        r = subprocess.run([sys.executable, self.original, "soto-fidelity", "--gtf", gtfp, "--gff", gffp, "--soto", sotop, "--out", os.path.join(self.tmp, "orig")], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        orig = {x["gene"]: x for x in tsv(os.path.join(self.tmp, "orig.per_gene.tsv"))}
        truth = write_truth(self.tmp, "mine", genes)
        out = os.path.join(self.tmp, "mine.score")
        r = run_scorer(truth, out, "--gtf", gtfp, "--unit", "transcript", "--compete", "soto")
        self.assertEqual(r.returncode, 0, r.stderr)
        mine = {x["gene"]: x for x in tsv(out + ".per_gene.tsv")}
        self.assertEqual(sorted(orig), sorted(mine))
        self.assertGreater(sum(1 for x in orig.values() if x["matched"] == "1"), 100)
        self.assertGreater(sum(1 for x in orig.values() if x["fused"] == "1"), 8)
        self.assertGreater(sum(1 for x in orig.values() if x["matched"] == "0"), 20)
        for g, o in orig.items():
            m = mine[g]
            self.assertEqual(o["matched"], m["matched"], g)
            if o["matched"] == "1":
                self.assertEqual((o["locus"], o["d5"], o["d3"], o["fused"]), (m["unit"], m["d5"], m["d3"], str(int(int(m["n_genes_span"]) >= 2))), g)
                self.assertAlmostEqual(float(o["exonic_jaccard"]), float(m["jaccard"]), delta=0.02, msg=g)
                self.assertAlmostEqual(float(o["width_ratio"]), float(m["width_ratio"]), delta=0.02, msg=g)


# ------------------------------------------------------------------------------------------------- run_soto_pilot.sh
STUB_RLOCK = """#!/bin/bash
# test stub of tools/rlock.sh: record the lock class and the command, then run the command (no lock, no timeout)
echo "$*" >> "$RS_LOCKLOG"
class=$1; shift
exec "$@"
"""
STUB_DRIVER = """#!/bin/bash
# test stub of tools/rustle_pipeline.sh. `assemble`: the annotation as an assembled GTF (the families input lacks TG_SH.1, as a bridge transcript would be); `families`: the nodes
# (fam.loci.gff3) and a clusters file, or a shard-budget stop / a failure when asked to.
stage=$1; shift
while [ $# -gt 0 ]; do case "$1" in --bam) BAM=$2; shift 2;; --out) OUT=$2; shift 2;; *) shift;; esac; done
echo "$stage $BAM $OUT" >> "$RS_LOCKLOG.driver"
HERE=$(dirname "$0")/../bench/ideal_expression
case $stage in
  assemble)
    python3 "$HERE/sim_genes.py" control --out "${BAM%.bam}" --to "$OUT.ctrl" > /dev/null
    mv "$OUT.ctrl.gtf" "$OUT.gtf"; rm -f "$OUT.ctrl.loci.gff3"
    grep -v 'transcript_id "TG_SH.1"' "$OUT.gtf" > "$OUT.families.gtf.tmp"; mv "$OUT.families.gtf.tmp" "$OUT.families.gtf"; touch "$OUT.families.gtf"
    echo "[copy_assign] ⭐ ASSEMBLY POLISH (full): 12 transcripts -> ISM dropped 1 -> mono floor [4, 3] (per contig) reads dropped 2 -> 9 kept" > "$OUT.assemble.log" ;;
  families)
    if [ -n "$STUB_FAMILIES_BUDGET_STOP" ]; then
      mkdir -p "$MM2_SHARD_DIR/key/q"; echo "[mm2_shard] budget: shard 2/4 not started (remaining 3s, estimate 40s); 1/4 shards done" >> "$MM2_SHARD_DIR/key/q/wrapper.log"
      if [ -n "$STUB_FAMILIES_PROGRESS" ]; then n=$(ls "$MM2_SHARD_DIR/key/q"/shard_*.ok 2>/dev/null | wc -l); touch "$MM2_SHARD_DIR/key/idx.mmi" "$MM2_SHARD_DIR/key/q/shard_$(printf %04d "$n").ok"; fi
      echo "Error: minimap2 all-vs-all failed" >&2; exit 1
    fi
    if [ -n "$STUB_FAMILIES_FAIL" ]; then echo "Error: something else" >&2; exit 1; fi
    python3 "$HERE/sim_genes.py" control --out "${BAM%.bam}" --to "$OUT.ctrl" > /dev/null
    mv "$OUT.ctrl.loci.gff3" "$OUT.fam.loci.gff3"; rm -f "$OUT.ctrl.gtf"
    printf 'cluster_id\\tsize\\nMCL0\\t2\\n' > "$OUT.fam.clusters.tsv" ;;
esac
echo "[rustle_pipeline] stub $stage done" >&2
"""
STUB_MINIMAP2 = """#!/usr/bin/env python3
# test stub of minimap2: --version, else a SAM with every read placed at the first exon of its source transcript (RS_STUB_TRUTH = a transcripts table); the FASTQ is the last argument
import csv, os, sys
if sys.argv[1:2] == ["--version"]:
    print("2.30-stub"); sys.exit(0)
sys.stderr.write("[M::main] Version: 2.30-stub\\n[M::main] CMD: minimap2 " + " ".join(sys.argv[1:]) + "\\n")
truth = {r["transcript"]: r for r in csv.DictReader(open(os.environ["RS_STUB_TRUTH"]), delimiter="\\t")}
print("@HD\\tVN:1.6\\tSO:unsorted")
for c, n in (("chr1", 120000), ("chr2", 120000), ("chr3", 450000)):
    print(f"@SQ\\tSN:{c}\\tLN:{n}")
with open(sys.argv[-1]) as fh:
    for h in fh:
        name = h[1:].strip(); seq = fh.readline().strip(); fh.readline(); q = fh.readline().strip()
        r = truth.get(name.split("|")[0])
        if r is None:
            print(f"{name}\\t4\\t*\\t0\\t0\\t*\\t*\\t0\\t0\\t{seq}\\t{q}")
        else:
            pos = int(r["canon_blocks"].split(",")[0].split("-")[0]) + 1
            print(f"{name}\\t0\\t{r['chrom']}\\t{pos}\\t60\\t{len(seq)}M\\t*\\t0\\t0\\t{seq}\\t{q}")
"""


def write_exec(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, 0o755)


def git(root, *args, date=None):
    env = dict(os.environ, GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@t", GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@t")
    if date:
        env["GIT_COMMITTER_DATE"] = env["GIT_AUTHOR_DATE"] = date
    return subprocess.run(["git", "-C", root, *args], capture_output=True, text=True, env=env, check=True).stdout.strip()


class RunnerTests(unittest.TestCase):
    STAGES = ["truth", "truthqc", "map", "cut", "bamqc", "assemble", "families", "anno", "score", "score-ctrl", "score-real", "report", "gates-sim", "gates", "cleanup", "status"]
    DEPTHS = ("4", "12")

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="soto_runner_test_")
        os.makedirs(os.path.join(cls.tmp, "world"))
        cls.world = build_world(os.path.join(cls.tmp, "world"))
        cls.repo = os.path.join(cls.tmp, "repo")
        ie = os.path.join(cls.repo, "bench", "ideal_expression")
        for d in (ie, os.path.join(cls.repo, "tools"), os.path.join(cls.repo, "src")):
            os.makedirs(d)
        for f in ("run_soto_pilot.sh", "sim_genes.py", "sim_windows.py", "score_soto_nodes.py", "map_reads.py"):
            shutil.copy(os.path.join(HERE, f), ie)
        for f in ("lib.py", "sim.py", "SOTO_FAMILY_LOCUS_FIDELITY_2026-09-22.md"):
            shutil.copy(os.path.join(BENCH, f), os.path.join(cls.repo, "bench"))
        write_exec(os.path.join(cls.repo, "tools", "rlock.sh"), STUB_RLOCK)
        write_exec(os.path.join(cls.repo, "tools", "rustle_pipeline.sh"), STUB_DRIVER)
        write_exec(os.path.join(cls.repo, "tools", "mm2_shard.sh"), "#!/bin/bash\nexit 0\n")
        with open(os.path.join(cls.repo, "src", "lib.rs"), "w") as fh:
            fh.write("// stub\n")
        with open(os.path.join(cls.repo, "Cargo.toml"), "w") as fh:
            fh.write("[package]\n")
        git(cls.repo, "init", "-q")
        git(cls.repo, "add", "-A")
        git(cls.repo, "commit", "-q", "-m", "init", date="2026-01-01T00:00:00")
        cls.bin = os.path.join(cls.tmp, "bin")
        os.makedirs(cls.bin)
        for b in ("copy_assign", "mcl_families", "as_table", "family_score"):
            write_exec(os.path.join(cls.bin, b), f"#!/bin/bash\necho {b} stub\n")
        cls.stubs = os.path.join(cls.tmp, "stubs")
        os.makedirs(cls.stubs)
        write_exec(os.path.join(cls.stubs, "minimap2"), STUB_MINIMAP2)
        cls.index = os.path.join(cls.tmp, "toy.mmi")
        with open(cls.index, "w") as fh:
            fh.write("not an index; the stub minimap2 ignores it\n")
        # toy real-library products for score-real: an old-style GTF, and a current-style run with nodes
        cls.real = os.path.join(cls.tmp, "real")
        os.makedirs(cls.real)
        write_gtf(os.path.join(cls.real, "old.gtf"), [unit("chr1", "+", "G1", "T1", [(1000, 1400), (2000, 2300), (3000, 3500)]), unit("chr1", "-", "G2", "T2", [(5000, 5400), (5800, 6200)])])
        write_gtf(os.path.join(cls.real, "f1v2.families.gtf"), [unit("chr1", "+", "G1", "T1", [(1000, 1400), (2000, 2300), (3000, 3500)])])
        write_gtf(os.path.join(cls.real, "testis.families.gtf"), [unit("chr1", "+", "G1", "T1", [(1000, 1400), (2000, 2300), (3000, 3500)])])
        write_nodes(os.path.join(cls.real, "testis.fam.loci.gff3"), {"N1": node("chr1", "+", [(1000, 1400), (2000, 2300), (3000, 3500)])})
        with open(os.path.join(cls.real, "testis.assemble.log"), "w") as fh:
            fh.write("[copy_assign] ⭐ ASSEMBLY POLISH (full): 100 transcripts -> ISM dropped 5 -> 90 kept\n")
        with open(os.path.join(cls.real, "reads.tsv"), "w") as fh:
            fh.write("gene_id\tprim_on_exons\nTG_A\t9\nTG_B\t0\n")
        cls.runner = os.path.join(ie, "run_soto_pilot.sh")
        cls.work = os.path.join(cls.tmp, "work")
        cls.lockfile = os.path.join(cls.tmp, "locks.log")
        cls._full = False

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def env(self, work=None, **extra):
        e = {k: v for k, v in os.environ.items() if not k.startswith("RUSTLE_") and not k.startswith("RS_")}
        w = work or self.work
        e.update(RS_WORK=w, RS_BIN=self.bin, RS_GENOME=self.world["fa"], RS_GFF=self.world["gff"], RS_INDEX=self.index, RS_GENES=self.world["genes"], RS_S1C=self.world["s1c"],
                 RS_SPLIT=self.world["split"], RS_PAD="2000", RS_SEED="7", RS_DEPTHS=" ".join(self.DEPTHS), RS_PARTS="3", RS_BASELINE_PER_GENE="", RS_LOCKLOG=self.lockfile,
                 RS_STUB_TRUTH=os.path.join(w, "d12", "reads.transcripts.tsv"), RS_REAL_A119B_OLD=os.path.join(self.real, "old.gtf"), RS_REAL_A119B_F1V2=os.path.join(self.real, "f1v2.families.gtf"),
                 RS_REAL_TESTIS=os.path.join(self.real, "testis"), RS_A119B_READS=os.path.join(self.real, "reads.tsv"), PATH=self.stubs + os.pathsep + e.get("PATH", ""))
        e.update(extra)
        return e

    def stage(self, *args, ok=0, **envx):
        r = subprocess.run(["bash", self.runner, *args], capture_output=True, text=True, env=self.env(**envx))
        if ok is not None:
            self.assertEqual(r.returncode, ok, f"{args}\n{r.stdout}\n{r.stderr}")
        return r

    def D(self, depth):
        return os.path.join(self.work, f"d{depth}")

    def locks(self):
        return slurp(self.lockfile).splitlines() if os.path.exists(self.lockfile) else []

    def full_run(self):
        cls = type(self)
        if cls._full:
            return
        self.stage("truth")
        self.stage("truthqc")
        cls._map_outputs = [self.stage("map").stdout for _ in range(3)]
        self.stage("cut")
        for d in self.DEPTHS:
            self.stage("bamqc", d)
        for d in self.DEPTHS:
            self.stage("assemble", d)
            self.stage("families", d)
            self.stage("anno", d)
        self.stage("score", self.DEPTHS[0])      # the calibration depth freezes the single-exon cut
        for d in self.DEPTHS[1:]:
            self.stage("score", d)
        for d in self.DEPTHS:
            self.stage("score-ctrl", d)
        self.stage("score-real")
        self.stage("gates-sim")
        cls._full = True

    # --- no products needed
    def test_help_lists_every_stage_the_depths_and_the_machine_rules(self):
        r = self.stage("--help")
        for st in self.STAGES:
            self.assertIn(st, r.stdout)
        for word in ("rlock.sh", "RS_WORK", "RS_DEPTHS", "20261008", "ceil(read Mb / 35)", "no progress"):
            self.assertIn(word, r.stdout)

    def test_help_and_status_only_read_and_create_no_directory(self):
        work = os.path.join(self.tmp, "readonly_work")
        for st in ("help", "--help", "-h", "status"):
            self.stage(st, work=work)
        self.assertFalse(os.path.exists(work), "help and status must not create the work directory")

    def test_refuses_to_run_with_a_rustle_variable_in_the_calling_environment(self):
        r = self.stage("status", ok=2, RUSTLE_JUNCTION_MAJORITY="1")
        self.assertIn("refusing", r.stderr)
        self.assertIn("RUSTLE_JUNCTION_MAJORITY", r.stderr)

    def test_dry_run_prints_the_commands_with_their_lock_class_and_writes_nothing(self):
        work = os.path.join(self.tmp, "dry_work")
        out = {st: self.stage("--dry-run", *st.split(), work=work).stdout for st in ("truth", "map", "cut", "assemble 4", "families 4", "score 4", "gates 4", "cleanup")}
        self.assertFalse(os.path.exists(work), "a dry run must create nothing")
        self.assertIn("DRY light", out["truth"]); self.assertIn("sim_genes.py", out["truth"]); self.assertIn("PYTHONHASHSEED=0", out["truth"]); self.assertIn("--reads-per-tx 12", out["truth"])
        self.assertIn("--reads-per-tx 4", out["truth"])
        self.assertIn("--seed 7", out["truth"])
        self.assertIn("DRY heavy", out["map"]); self.assertIn("--max-parts-per-call 1", out["map"]); self.assertIn("--parts 3", out["map"]); self.assertIn("/d12/reads", out["map"])
        self.assertIn("samtools view", out["cut"])
        self.assertIn("DRY heavy", out["assemble 4"]); self.assertIn("rustle_pipeline.sh assemble", out["assemble 4"]); self.assertIn("--threads 4", out["assemble 4"]); self.assertIn("/d4/asm", out["assemble 4"])
        fam = out["families 4"]
        self.assertIn("DRY heavy", fam); self.assertIn("rustle_pipeline.sh families", fam)
        self.assertRegex(fam, r"RUSTLE_MINIMAP2=\S+/d4/mm2_shard_exec\.sh")
        self.assertRegex(fam, r"MM2_SHARD_DIR=\S+/d4/mm2")
        self.assertIn("MM2_SHARD_BUDGET_S=480", fam)
        self.assertRegex(fam, r"MM2_SHARD_DEADLINE=\d{10}")
        self.assertIn("DRY light", out["score 4"]); self.assertIn("score_soto_nodes.py", out["score 4"]); self.assertNotIn("DRY heavy", out["score 4"])
        self.assertIn("--unit node", out["score 4"]); self.assertIn("--freeze-single-cut", out["score 4"])
        self.assertIn("DRY light", out["gates 4"])

    def test_the_default_process_count_is_one_and_a_single_process_is_a_light_job(self):
        work = os.path.join(self.tmp, "dry_work2")
        one = self.stage("--dry-run", "truth", work=work).stdout
        self.assertIn("--procs 1", one); self.assertNotIn("DRY heavy", one)
        many = self.stage("--dry-run", "truth", work=work, RS_PROCS="3").stdout
        self.assertIn("DRY heavy", many); self.assertIn("--procs 3", many)

    def test_the_number_of_mapping_parts_is_the_read_megabases_over_35_rounded_up(self):
        work = os.path.join(self.tmp, "parts_work")
        os.makedirs(os.path.join(work, "d12"))
        for bases, want in ((1, 1), (35_000_000, 1), (35_000_001, 2), (740_000_000, 22)):
            with open(os.path.join(work, "d12", "reads.sim.json"), "w") as fh:
                json.dump(dict(read_bases=bases), fh)
            env = {k: v for k, v in self.env(work=work).items() if k != "RS_PARTS"}
            r = subprocess.run(["bash", self.runner, "--dry-run", "map"], capture_output=True, text=True, env=env)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn(f"--parts {want} ", r.stdout + " ", (bases, r.stdout))

    def test_a_depth_argument_is_required_where_a_stage_is_per_arm_and_must_be_one_of_the_depths(self):
        work = os.path.join(self.tmp, "arg_work")
        r = self.stage("assemble", ok=2, work=work)
        self.assertIn("depth", r.stderr)
        r = self.stage("assemble", "7", ok=2, work=work)
        self.assertIn("one of", r.stderr)

    def test_stages_refuse_to_run_before_their_inputs_exist(self):
        work = os.path.join(self.tmp, "order_work")
        for st, msg in (("truthqc", "truth first"), ("map", "truth first"), ("cut", "map first"), ("bamqc 4", "map first"), ("assemble 4", "map first"), ("families 4", "assemble first"),
                        ("anno 4", "truth first"), ("score 4", "families first"), ("score-ctrl 4", "anno first"), ("score-real", "truth first"), ("gates-sim", "truth first"),
                        ("gates 4", "truth first"), ("cleanup", "map first"), ("report 4", "families first")):
            r = self.stage(*st.split(), ok=2, work=work)
            self.assertIn(msg, r.stderr, st)

    def test_assemble_needs_a_g4_that_passed_and_map_needs_g1_and_g2_that_passed(self):
        work = os.path.join(self.tmp, "gate_work")
        os.makedirs(os.path.join(work, "d4")); os.makedirs(os.path.join(work, "d12"))
        with open(os.path.join(work, "d12", "reads.fq"), "w") as fh:
            fh.write("@a\nAC\n+\n~~\n")
        with open(os.path.join(work, "d12", "reads.bam"), "w") as fh:
            fh.write("x")
        with open(os.path.join(work, "d4", "reads.bam"), "w") as fh:
            fh.write("x")
        self.assertIn("truthqc first", self.stage("map", ok=2, work=work).stderr)
        with open(os.path.join(work, "d4", "reads.gates.json"), "w") as fh:
            json.dump(dict(G1=dict(ok=True, lines=[]), G2=dict(ok=True, lines=[])), fh)
        with open(os.path.join(work, "d12", "reads.gates.json"), "w") as fh:
            json.dump(dict(G1=dict(ok=True, lines=[]), G2=dict(ok=False, lines=["G2 FAIL"])), fh)
        self.assertIn("G2 did not pass", self.stage("map", ok=2, work=work).stderr)
        self.assertIn("bamqc first", self.stage("assemble", "4", ok=2, work=work).stderr)
        with open(os.path.join(work, "d4", "reads.g4.json"), "w") as fh:
            json.dump(dict(ok=False, lines=["G4 FAIL"]), fh)
        self.assertIn("G4 did not pass", self.stage("assemble", "4", ok=2, work=work).stderr)

    def test_map_refuses_while_the_shallow_arms_are_not_proved_to_be_the_deep_pool_s_first_reads(self):
        work = os.path.join(self.tmp, "gdepth_work")
        for d in self.DEPTHS:
            os.makedirs(os.path.join(work, f"d{d}"))
            with open(os.path.join(work, f"d{d}", "reads.gates.json"), "w") as fh:
                json.dump(dict(G1=dict(ok=True, lines=[]), G2=dict(ok=True, lines=[])), fh)
        with open(os.path.join(work, "d12", "reads.fq"), "w") as fh:
            fh.write("@a\nAC\n+\n~~\n")
        self.assertIn("G-depth did not pass", self.stage("map", ok=2, work=work).stderr)          # G1 and G2 passed, the depth gate never ran
        with open(os.path.join(work, "depth_subset.json"), "w") as fh:
            json.dump(dict(ok=False, lines=["G-depth FAIL"]), fh)
        self.assertIn("G-depth did not pass", self.stage("map", ok=2, work=work).stderr)          # it ran and failed
        self.assertFalse(os.path.exists(os.path.join(work, "d12", "map.parts")), "nothing may be mapped behind a failed gate")

    # --- the whole toy run
    def test_full_toy_run_stage_by_stage(self):
        self.full_run()
        W = self.work
        for f in ("truth.run.log", "depth_subset.json", "single_cut.json", "g3_sim.json", "real/old_tx.summary.json", "real/testis_node.summary.json"):
            self.assertTrue(os.path.exists(os.path.join(W, f)), f)
        for d in self.DEPTHS:
            D = self.D(d)
            for f in ("reads.fq", "reads.transcripts.tsv", "reads.targets.tsv", "reads.sim.json", "reads.gates.json", "reads.bam", "reads.bam.bai", "reads.g4.json", "asm.gtf", "asm.families.gtf",
                      "asm.fam.loci.gff3", "asm.fam.clusters.tsv", "asm.assemble.log", "asm.run.log", "ctrl0.gtf", "ctrl0.loci.gff3", "ctrl50k.loci.gff3", "score.per_gene.tsv", "score.summary.json",
                      "score.report.txt", "score_tx.summary.json", "score_locus.summary.json", "score_soto.summary.json", "score_raw.summary.json", "score_ctrl0.summary.json",
                      "score_ctrl0_locus.summary.json", "score_ctrl50k.summary.json", "mm2_shard_exec.sh"):
                self.assertTrue(os.path.exists(os.path.join(D, f)), f"d{d}/{f}")
        self.assertEqual(jload(os.path.join(self.D("12"), "reads.sim.json"))["reads"], 120)
        self.assertEqual(jload(os.path.join(self.D("4"), "reads.sim.json"))["reads"], 40)
        self.assertIn("seed\t7", slurp(os.path.join(W, "truth.run.log")))
        outs = type(self)._map_outputs
        self.assertIn("mapping parts remain", outs[0]); self.assertIn("mapping parts remain", outs[1]); self.assertIn("mapping complete", outs[2])
        self.assertTrue(jload(os.path.join(W, "depth_subset.json"))["ok"])
        self.assertTrue(os.access(os.path.join(self.D("4"), "mm2_shard_exec.sh"), os.X_OK))

    def test_the_shallow_arm_s_bam_is_the_deep_pool_cut_by_read_name(self):
        self.full_run()
        deep, shallow = os.path.join(self.D("12"), "reads"), os.path.join(self.D("4"), "reads")
        names = {n for n, *_ in fastq_records(shallow + ".fq")}
        cut = subprocess.run(["samtools", "view", shallow + ".bam"], capture_output=True, text=True, check=True).stdout.splitlines()
        cut_names = {ln.split("\t")[0] for ln in cut}
        self.assertEqual(cut_names, names)
        self.assertEqual(len(cut), len(names))                      # one record per read (the stub maps every read once)
        full = subprocess.run(["samtools", "view", "-c", "-F", "2304", deep + ".bam"], capture_output=True, text=True, check=True).stdout.strip()
        self.assertEqual(full, "120")
        want = sorted(ln for ln in subprocess.run(["samtools", "view", deep + ".bam"], capture_output=True, text=True, check=True).stdout.splitlines() if ln.split("\t")[0] in names)
        self.assertEqual(sorted(cut), want)                         # the same records as the pool's, only those reads

    def test_assemble_prints_the_polish_counts_and_the_families_stage_uses_a_shim_and_a_shard_directory_per_arm(self):
        self.full_run()
        self.assertIn("ASSEMBLY POLISH", slurp(os.path.join(self.D("4"), "polish.txt")))
        locks = self.locks()
        fam = [ln for ln in locks if "rustle_pipeline.sh families" in ln]
        self.assertTrue(fam)
        for d in self.DEPTHS:
            line = next(ln for ln in fam if f"/d{d}/asm" in ln)
            self.assertIn(f"RUSTLE_MINIMAP2={self.D(d)}/mm2_shard_exec.sh", line)
            self.assertIn(f"MM2_SHARD_DIR={self.D(d)}/mm2", line)
            self.assertIn("MM2_SHARD_BUDGET_S=480", line)
        self.assertIn("exec bash", slurp(os.path.join(self.D("4"), "mm2_shard_exec.sh")))
        self.assertIn("tools/mm2_shard.sh", slurp(os.path.join(self.D("4"), "mm2_shard_exec.sh")))

    def test_assemble_removes_the_best_as_table_of_its_prefix_because_the_driver_reuses_any_table_that_names_the_same_bam_path(self):
        self.full_run()
        stale = [os.path.join(self.D("4"), "asm.molecules.tsv"), os.path.join(self.D("4"), "asm.molecules.tsv.asbin")]
        for f in stale:
            with open(f, "w") as fh:
                fh.write("#as_table\tbam=stale\n")
        self.stage("assemble", "4")
        for f in stale:
            self.assertFalse(os.path.exists(f), f"{f} survived: a re-cut BAM at the same path would be assembled with the old table")

    def test_lock_classes_follow_the_machine_rules(self):
        self.full_run()
        locks = self.locks()
        by = lambda needle: {ln.split()[0] for ln in locks if needle in ln}
        self.assertEqual(by("map_reads.py"), {"heavy"})
        self.assertEqual(by("rustle_pipeline.sh assemble"), {"heavy"})
        self.assertEqual(by("rustle_pipeline.sh families"), {"heavy"})
        self.assertEqual(by("score_soto_nodes.py"), {"light"})
        self.assertEqual(by("sim_genes.py bamqc"), {"light"})
        self.assertEqual(by("sim_genes.py gates"), {"light"})
        self.assertEqual(by("sim_genes.py control"), {"light"})
        self.assertEqual(by("samtools view"), {"light"})
        self.assertEqual(len([ln for ln in locks if "map_reads.py" in ln]), 3)   # one part per call
        self.assertTrue(all("pkill" not in ln and "nohup" not in ln for ln in locks))

    def test_families_reports_a_spent_shard_budget_as_call_again_and_any_other_failure_as_a_failure(self):
        self.full_run()
        clusters = os.path.join(self.D("4"), "asm.fam.clusters.tsv")
        try:
            r = self.stage("families", "4", ok=75, STUB_FAMILIES_BUDGET_STOP="1")
            self.assertIn("call again", r.stderr)
            self.assertFalse(os.path.exists(clusters), "a stale clusters file must not survive a budget stop")
            r = self.stage("families", "4", ok=1, STUB_FAMILIES_FAIL="1")
            self.assertIn("families failed", r.stderr)
        finally:
            self.stage("families", "4")
        self.assertTrue(os.path.exists(clusters))

    def test_a_budget_stop_that_made_no_progress_since_the_previous_call_is_a_failure_not_another_call_again(self):
        self.full_run()
        progress = os.path.join(self.D("4"), "families.progress")
        try:
            first = self.stage("families", "4", ok=75, STUB_FAMILIES_BUDGET_STOP="1", STUB_FAMILIES_PROGRESS="1")        # the first stop has nothing to be compared with
            self.assertIn("call again", first.stderr)
            second = self.stage("families", "4", ok=75, STUB_FAMILIES_BUDGET_STOP="1", STUB_FAMILIES_PROGRESS="1")       # one more shard is done
            self.assertIn("call again", second.stderr)
            stalled = self.stage("families", "4", ok=1, STUB_FAMILIES_BUDGET_STOP="1")                                  # nothing new: a single step does not fit the budget
            self.assertIn("no progress", stalled.stderr)
            self.assertNotIn("call again", stalled.stderr)
            shards = os.path.join(self.D("4"), "mm2", "key")
            for rel in ("idx.mmi", "q/shard_0000.ok", "q/shard_0001.ok"):
                os.remove(os.path.join(shards, rel))                                                                    # the shard cache was cleared by hand: that is a restart, not a stall
            self.assertIn("call again", self.stage("families", "4", ok=75, STUB_FAMILIES_BUDGET_STOP="1").stderr)
        finally:
            self.stage("families", "4")
        self.assertFalse(os.path.exists(progress), "a finished stage forgets its progress record")

    def test_the_single_exon_cut_is_frozen_by_the_calibration_depth_first_and_read_by_the_others(self):
        work = os.path.join(self.tmp, "cut_work")
        r = self.stage("score", "12", ok=2, work=work)
        self.assertIn("freeze", r.stderr)
        self.full_run()
        cut = jload(os.path.join(self.work, "single_cut.json"))
        self.assertEqual((cut["half"], cut["percentile"]), ("dev", 5))
        self.assertIn("FROZE", slurp(os.path.join(self.D("4"), "score.report.txt")) + slurp(os.path.join(self.D("4"), "score.log")))
        self.assertIn("frozen cut", slurp(os.path.join(self.D("12"), "score.report.txt")))
        before = slurp(os.path.join(self.work, "single_cut.json"))
        self.stage("score", "4")      # a second call reads the frozen cut, it does not compute another
        self.assertEqual(slurp(os.path.join(self.work, "single_cut.json")), before)
        self.assertIn("frozen cut", slurp(os.path.join(self.D("4"), "score.report.txt")))

    def test_only_the_calibration_depth_freezes_the_single_exon_cut(self):
        work = os.path.join(self.tmp, "dry_work3")
        shallow = self.stage("--dry-run", "score", "4", work=work).stdout
        deep = self.stage("--dry-run", "score", "12", work=work).stdout
        self.assertIn("--freeze-single-cut", shallow)
        self.assertNotIn("--freeze-single-cut", deep)
        self.assertIn("--freeze-single-cut", self.stage("--dry-run", "score", "12", work=work, RS_CUT_DEPTH="12").stdout)      # the calibration depth is a setting
        self.assertNotIn("--freeze-single-cut", self.stage("--dry-run", "score", "4", work=work, RS_CUT_DEPTH="12").stdout)

    def test_the_primary_score_uses_the_nodes_and_the_variants_the_transcripts_the_loci_and_the_raw_gtf(self):
        self.full_run()
        D = self.D("12")
        ctx = lambda f: jload(os.path.join(D, f + ".summary.json"))["context"]
        self.assertIn("unit=node", ctx("score")); self.assertIn("compete=all", ctx("score"))
        self.assertIn("unit=transcript", ctx("score_tx")); self.assertIn("unit=locus", ctx("score_locus")); self.assertIn("compete=soto", ctx("score_soto")); self.assertIn("unit=transcript", ctx("score_raw"))
        tx = {r["gene"]: r["class"] for r in tsv(os.path.join(D, "score_tx.per_gene.tsv"))}
        raw = {r["gene"]: r["class"] for r in tsv(os.path.join(D, "score_raw.per_gene.tsv"))}
        self.assertEqual(tx["TG_SH"], "no_locus")       # TG_SH.1 is not in asm.families.gtf
        self.assertEqual(raw["TG_SH"], "node_ok")       # but is in asm.gtf
        node_ok = {r["gene"]: r["class"] for r in tsv(os.path.join(D, "score.per_gene.tsv"))}
        self.assertEqual(node_ok["TG_A"], "node_ok")
        summ = jload(os.path.join(D, "score.summary.json"))
        self.assertGreater(len(summ["assembly_polish"]), 0)          # the polish counts of the assemble log travel into the report
        self.assertEqual(summ["population"]["simulable"], 6)         # TG_A, TG_B, TG_NC, TG_SH, TG_P, TG_BIG; TG_TINY and TG_MISS are listed apart
        self.assertEqual(summ["population"]["not_simulable"]["genes"], 2)

    def test_the_two_controls_are_scored_and_the_shifted_one_collapses(self):
        self.full_run()
        pos = jload(os.path.join(self.D("12"), "score_ctrl0.summary.json"))
        neg = jload(os.path.join(self.D("12"), "score_ctrl50k.summary.json"))
        self.assertEqual((pos["classes"]["fused_excess"], pos["classes"]["extent_wrong"], pos["classes"]["no_locus"]), (0, 0, 0))
        self.assertGreater(pos["classes"]["node_ok"], 0)
        self.assertEqual(neg["classes"]["node_ok"], 0)
        self.assertEqual(neg["classes"]["no_locus"], neg["population"]["simulable"])

    def test_the_real_library_comparators_are_scored_by_the_same_scorer_on_the_ideal_arm_s_truth(self):
        self.full_run()
        R = os.path.join(self.work, "real")
        for f in ("old_tx", "old_locus", "f1v2_tx", "f1v2_locus", "testis_tx", "testis_locus", "testis_node"):
            s = jload(os.path.join(R, f + ".summary.json"))
            self.assertIn("arm=real", s["context"])
            self.assertEqual(s["population"]["simulable"], 6)
        old = {r["gene"]: r["class"] for r in tsv(os.path.join(R, "old_tx.per_gene.tsv"))}
        self.assertEqual(old["TG_B"], "no_reads")          # the read table says TG_B has no primary reads in the real library
        self.assertEqual(old["TG_A"], "node_ok")
        self.assertIn("ASSEMBLY POLISH", slurp(os.path.join(R, "testis_node.report.txt")))

    def test_report_prints_the_comparators_beside_the_arm(self):
        self.full_run()
        r = self.stage("report", "12")
        self.assertIn("COMPARATORS", r.stdout)
        for label in ("this arm", "ctrl0 node", "ctrl50k node", "old A119b tx", "testis node"):
            self.assertIn(label, r.stdout)

    def test_gates_print_pass_or_fail_per_gate_and_a_table_is_valid_only_if_all_pass(self):
        self.full_run()
        for d in self.DEPTHS:
            r = self.stage("gates", d)
            for g in ("G1 PASS", "G2 PASS", "G3 PASS", "G4 PASS", "G5 PASS", "G6 PASS", "G7 PASS"):
                self.assertIn(g, r.stdout, (d, r.stdout))
            if d != self.DEPTHS[-1]:
                self.assertIn("G-depth PASS", r.stdout)
            self.assertRegex(r.stdout, r"(?m)^VALID\b")
            self.assertTrue(jload(os.path.join(self.D(d), "gates.json"))["valid"])
            for word in ("node_ok", "no_locus", "Jaccard"):
                self.assertNotIn(word, r.stdout)       # a gate prints PASS or FAIL, not an outcome statistic

    def test_g6_fails_when_the_positive_control_loses_a_gene_it_cannot_lose(self):
        self.full_run()
        f = os.path.join(self.D("4"), "score_ctrl0.summary.json")
        keep = slurp(f)
        try:
            s = json.loads(keep)
            s["classes"]["extent_wrong"] = 1
            with open(f, "w") as fh:
                json.dump(s, fh)
            r = self.stage("gates", "4", ok=1)
            self.assertIn("G6 FAIL", r.stdout)
            self.assertRegex(r.stdout, r"(?m)^INVALID\b")
        finally:
            with open(f, "w") as fh:
                fh.write(keep)

    def test_g7_fails_when_the_shifted_control_does_not_collapse(self):
        self.full_run()
        f = os.path.join(self.D("4"), "score_ctrl50k.summary.json")
        keep = slurp(f)
        try:
            s = json.loads(keep)
            s["classes"]["node_ok"] = jload(os.path.join(self.D("4"), "score_ctrl0.summary.json"))["classes"]["node_ok"]
            with open(f, "w") as fh:
                json.dump(s, fh)
            r = self.stage("gates", "4", ok=1)
            self.assertIn("G7 FAIL", r.stdout)
        finally:
            with open(f, "w") as fh:
                fh.write(keep)

    def test_gates_are_invalid_when_the_fastq_no_longer_matches_its_manifest(self):
        self.full_run()
        fq = os.path.join(self.D("12"), "reads.fq")
        keep = slurp(fq)
        try:
            with open(fq, "w") as fh:
                fh.write("".join(keep.splitlines(True)[:-4]))   # drop the last read
            r = self.stage("gates", "12", ok=1)
            self.assertIn("G1 FAIL", r.stdout)
            self.assertRegex(r.stdout, r"(?m)^INVALID\b")
        finally:
            with open(fq, "w") as fh:
                fh.write(keep)

    def test_a_crashed_recheck_cannot_pass_on_a_stale_result(self):
        self.full_run()
        r = self.stage("gates", "12", ok=1, RS_GENOME=os.path.join(self.tmp, "no_such_genome.fa"))
        self.assertIn("G1 FAIL", r.stdout)
        self.assertRegex(r.stdout, r"(?m)^INVALID\b")

    def test_g5_fails_when_a_binary_changed_after_assemble(self):
        self.full_run()
        b = os.path.join(self.bin, "copy_assign")
        keep = slurp(b)
        try:
            with open(b, "a") as fh:
                fh.write("# rebuilt\n")
            r = self.stage("gates", "12", ok=1)
            self.assertIn("G5 FAIL", r.stdout)
            self.assertIn("copy_assign", r.stdout)
        finally:
            with open(b, "w") as fh:
                fh.write(keep)
            os.chmod(b, 0o755)
            future = max(os.path.getmtime(os.path.join(self.bin, x)) for x in os.listdir(self.bin))
            for x in os.listdir(self.bin):
                os.utime(os.path.join(self.bin, x), (future, future))

    def test_g5_fails_when_a_script_changed_after_the_truth_stage(self):
        self.full_run()
        script = os.path.join(self.repo, "bench", "ideal_expression", "map_reads.py")
        keep = slurp(script)
        try:
            with open(script, "a") as fh:
                fh.write("# edited mid-run\n")
            r = self.stage("gates", "12", ok=1)
            self.assertIn("G5 FAIL", r.stdout)
            self.assertIn("map_reads.py", r.stdout)
        finally:
            with open(script, "w") as fh:
                fh.write(keep)
        self.stage("gates", "12")

    def test_g5_fails_for_a_source_commit_newer_than_the_binaries_unless_it_is_whitelisted(self):
        self.full_run()
        with open(os.path.join(self.repo, "src", "lib.rs"), "a") as fh:
            fh.write("// a later change\n")
        git(self.repo, "add", "-A")
        git(self.repo, "commit", "-q", "-m", "later", date="2099-01-01T00:00:00")
        sha = git(self.repo, "rev-parse", "--short", "HEAD")
        try:
            r = self.stage("gates", "12", ok=1)
            self.assertIn("G5 FAIL", r.stdout)
            self.assertIn(sha, r.stdout)
            r = self.stage("gates", "12", ok=0, RS_SRC_OK=sha)
            self.assertIn("G5 PASS", r.stdout)
            self.assertIn(sha, r.stdout)
        finally:
            git(self.repo, "reset", "-q", "--hard", "HEAD~1")

    def test_truth_refuses_to_overwrite_a_simulation_that_exists(self):
        self.full_run()
        before = sha256(os.path.join(self.D("12"), "reads.fq"))
        r = self.stage("truth", ok=2)
        self.assertIn("refusing to overwrite", r.stderr)
        self.assertEqual(sha256(os.path.join(self.D("12"), "reads.fq")), before)

    def test_cleanup_deletes_the_part_files_by_exact_name_and_only_after_the_bam_is_complete(self):
        self.full_run()
        pool = os.path.join(self.D("12"), "reads")
        parts = [f"{pool}.part{i}{ext}" for i in range(3) for ext in (".fq", ".bam")]
        self.assertTrue(all(os.path.exists(p) for p in parts))
        keep_other = os.path.join(self.D("12"), "reads.part_notes.txt")       # a look-alike that is not a part file must survive
        with open(keep_other, "w") as fh:
            fh.write("keep")
        self.stage("cleanup")
        self.assertFalse(any(os.path.exists(p) for p in parts))
        self.assertTrue(os.path.exists(keep_other))
        self.assertTrue(os.path.exists(pool + ".bam") and os.path.exists(pool + ".fq"))
        self.stage("cleanup")      # idempotent

    def test_status_lists_the_products_and_the_mapping_progress(self):
        work = os.path.join(self.tmp, "status_work")
        r = self.stage("status", work=work)
        self.assertIn("parts mapped", r.stdout)
        self.assertIn("d12", r.stdout)
        self.assertIn("single_cut.json", r.stdout)



if __name__ == "__main__":
    unittest.main()
