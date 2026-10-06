#!/usr/bin/env python3
"""Unit tests for famsim (stdlib unittest; no minimap2, no pysam needed except for the planting round-trip, which uses
samtools faidx + pysam if present). Run: python3 bench/famsim/test_famsim.py"""
import json
import os
import random
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from famsim.model import GeneModel, rc, seeded, synthetic_model  # noqa: E402
from famsim.ops import apply_ops  # noqa: E402

EXONS, INTRONS = [200, 150, 300, 120, 90, 400], [1000, 2000, 800, 1200, 900]


def fresh(seed=1):
    return synthetic_model(EXONS, INTRONS, random.Random(seed))


class TestModel(unittest.TestCase):
    def test_synthetic_is_canonical_and_spans_exons(self):
        m = fresh()
        self.assertTrue(m.canonical())
        self.assertEqual(len(m.seq), sum(EXONS) + sum(INTRONS))
        self.assertEqual([len(e) for e in m.exons], EXONS)
        self.assertEqual(len(m.chain_seq()), sum(EXONS))
        self.assertEqual(m.chain_junctions(), m.introns())

    def test_protected_positions_are_the_splice_dinucleotides(self):
        m = fresh()
        for s, e in m.introns():
            self.assertEqual(m.seq[s:s + 2], "GT"); self.assertEqual(m.seq[e - 2:e], "AG")
            for p in (s, s + 1, e - 2, e - 1):
                self.assertIn(p, m.protected())

    def test_edit_shifts_downstream_exons_and_refuses_boundaries(self):
        m = fresh()
        before = [(e.start, e.end) for e in m.exons]
        m.edit(m.exons[1].start + 10, 0, "ACGT")                       # insertion inside exon 2
        self.assertEqual(m.exons[1].end, before[1][1] + 4)
        self.assertEqual(m.exons[2].start, before[2][0] + 4)
        with self.assertRaises(ValueError):
            m.edit(m.exons[1].end, 2, "")                                  # the donor GT
        with self.assertRaises(ValueError):
            m.edit(m.exons[1].end - 5, 10, "")                             # straddles an exon end


class TestOps(unittest.TestCase):
    def test_snp_rate_and_protection(self):
        m = fresh(); t = fresh()
        apply_ops(m, [{"op": "snp", "rate": 1.0}], seeded(1, "x"))     # every unprotected base changes
        diff = sum(a != b for a, b in zip(m.seq, t.seq))
        self.assertEqual(diff, len(t.seq) - len(t.protected()))
        self.assertTrue(m.canonical())
        m2 = fresh(); apply_ops(m2, [{"op": "snp", "rate": 0.02}], seeded(1, "x"))
        self.assertEqual(m2.ops[0]["n"], round(0.02 * (len(t.seq) - len(t.protected()))))

    def test_snp_is_reproducible(self):
        a = fresh(); b = fresh()
        apply_ops(a, [{"op": "snp", "rate": 0.05}], seeded(7, "c")); apply_ops(b, [{"op": "snp", "rate": 0.05}], seeded(7, "c"))
        self.assertEqual(a.seq, b.seq)

    def test_exon_delete_internal_and_terminal(self):
        m = fresh(); apply_ops(m, [{"op": "exon_delete", "exon": 3}], seeded(1))
        self.assertEqual([e.label for e in m.exons], ["e1", "e2", "e4", "e5", "e6"])
        self.assertEqual(len(m.seq), sum(EXONS) + sum(INTRONS) - 300)
        self.assertTrue(m.canonical())
        m = fresh(); apply_ops(m, [{"op": "exon_delete", "exon": 1}], seeded(1))
        self.assertEqual(m.exons[0].label, "e2"); self.assertEqual(m.exons[0].start, 0)
        self.assertEqual(len(m.seq), sum(EXONS) + sum(INTRONS) - 200 - 1000); self.assertTrue(m.canonical())
        m = fresh(); apply_ops(m, [{"op": "exon_delete", "exon": 6}], seeded(1))
        self.assertEqual(m.exons[-1].label, "e5"); self.assertEqual(m.exons[-1].end, len(m.seq)); self.assertTrue(m.canonical())

    def test_splice_kill(self):
        m = fresh(); t = fresh(); apply_ops(m, [{"op": "splice_kill", "exon": 2}], seeded(1))
        self.assertEqual(len(m.seq), len(t.seq))
        self.assertFalse(m.exons[1].in_rna)
        self.assertEqual([e.label for e in m.rna_exons()], ["e1", "e3", "e4", "e5", "e6"])
        self.assertEqual(m.seq[m.exons[1].end:m.exons[1].end + 2], "CT")
        with self.assertRaises(ValueError):
            apply_ops(fresh(), [{"op": "splice_kill", "exon": 1}], seeded(1))

    def test_exon_insert_random_and_duplicate(self):
        m = fresh(); apply_ops(m, [{"op": "exon_insert", "after": 2, "length": 120}], seeded(1))
        self.assertEqual(len(m.exons), 7); self.assertEqual(m.exons[2].label, "ins1"); self.assertEqual(len(m.exons[2]), 120)
        self.assertTrue(m.canonical()); self.assertEqual(len(m.seq), sum(EXONS) + sum(INTRONS) + 124)
        m = fresh(); apply_ops(m, [{"op": "exon_insert", "after": 1, "source_exon": 3}], seeded(1))
        self.assertEqual(m.exons[1].label, "dupe3"); self.assertEqual(m.seq[m.exons[1].start:m.exons[1].end], m.seq[m.exons[3].start:m.exons[3].end])
        self.assertTrue(m.canonical())

    def test_exon_shuffle(self):
        m = fresh(); t = fresh(); apply_ops(m, [{"op": "exon_shuffle", "a": 2, "b": 4}], seeded(1))
        self.assertEqual([e.label for e in m.exons], ["e1", "e4", "e3", "e2", "e5", "e6"])
        self.assertEqual(m.seq[m.exons[1].start:m.exons[1].end], t.seq[t.exons[3].start:t.exons[3].end])
        self.assertEqual(m.seq[m.exons[3].start:m.exons[3].end], t.seq[t.exons[1].start:t.exons[1].end])
        self.assertEqual(len(m.seq), len(t.seq)); self.assertTrue(m.canonical())

    def test_invert_exon_intron_span(self):
        m = fresh(); t = fresh(); apply_ops(m, [{"op": "invert", "exon": 4}], seeded(1))
        e = m.exons[3]; self.assertTrue(e.inverted); self.assertFalse(e.in_rna)
        self.assertEqual(m.seq[e.start:e.end], rc(t.seq[t.exons[3].start:t.exons[3].end])); self.assertTrue(m.canonical())
        m = fresh(); apply_ops(m, [{"op": "invert", "intron": 2}], seeded(1))
        s, e_ = t.introns()[1]
        self.assertEqual(m.seq[s + 6:e_ - 6], rc(t.seq[s + 6:e_ - 6])); self.assertEqual(m.chain_seq(), t.chain_seq()); self.assertTrue(m.canonical())
        m = fresh(); a, b = t.exons[1].start - 50, t.exons[2].end + 50
        apply_ops(m, [{"op": "invert", "span": [a, b]}], seeded(1))
        self.assertEqual(m.seq[a:b], rc(t.seq[a:b]))
        self.assertEqual([e.label for e in m.exons], ["e1", "e3", "e2", "e4", "e5", "e6"])
        self.assertEqual([e.label for e in m.rna_exons()], ["e1", "e4", "e5", "e6"])
        with self.assertRaises(ValueError):
            apply_ops(fresh(), [{"op": "invert", "span": [10, 250]}], seeded(1))    # cuts exon 1

    def test_truncate(self):
        m = fresh(); apply_ops(m, [{"op": "truncate", "side": 5, "exons": 2}], seeded(1))
        self.assertEqual([e.label for e in m.exons], ["e3", "e4", "e5", "e6"]); self.assertEqual(m.exons[0].start, 0); self.assertTrue(m.canonical())
        m = fresh(); apply_ops(m, [{"op": "truncate", "side": 3, "exons": 1}], seeded(1))
        self.assertEqual(m.exons[-1].label, "e5"); self.assertEqual(m.exons[-1].end, len(m.seq)); self.assertTrue(m.canonical())
        m = fresh(); apply_ops(m, [{"op": "truncate", "side": 5, "bp": 250}], seeded(1))   # lands in intron 1 -> moves to exon 2
        self.assertEqual(m.exons[0].label, "e2"); self.assertEqual(m.exons[0].start, 0)
        m = fresh(); apply_ops(m, [{"op": "truncate", "side": 5, "bp": 50}], seeded(1))    # inside exon 1 -> exon shortened
        self.assertEqual(m.exons[0].label, "e1"); self.assertEqual(len(m.exons[0]), 150)

    def test_convert_uses_the_donor_s_current_sequence(self):
        t = fresh(); d = fresh(); apply_ops(d, [{"op": "snp", "rate": 0.1}], seeded(2))
        m = fresh(); apply_ops(m, [{"op": "convert", "from": "A2", "exon": 3}], seeded(3), donors={"A2": d})
        self.assertEqual(m.seq[m.exons[2].start:m.exons[2].end], d.seq[d.exons[2].start:d.exons[2].end])
        self.assertEqual(m.seq[:m.exons[2].start], t.seq[:t.exons[2].start])
        with self.assertRaises(ValueError):
            apply_ops(fresh(), [{"op": "convert", "from": "nope", "exon": 3}], seeded(3), donors={})

    def test_intron_resize(self):
        m = fresh(); apply_ops(m, [{"op": "intron_resize", "intron": 2, "length": 500}], seeded(1))
        self.assertEqual(m.introns()[1][1] - m.introns()[1][0], 500); self.assertTrue(m.canonical())
        self.assertEqual(len(m.seq), sum(EXONS) + sum(INTRONS) - 1500)

    def test_op_records_and_json_roundtrip(self):
        m = fresh(); apply_ops(m, [{"op": "snp", "rate": 0.01}, {"op": "exon_delete", "exon": 2}], seeded(1))
        self.assertEqual([o["op"] for o in m.ops], ["snp", "exon_delete"])
        m2 = GeneModel.from_json(json.loads(json.dumps(m.to_json())))
        self.assertEqual(m2.seq, m.seq); self.assertEqual([e.to_json() for e in m2.exons], [e.to_json() for e in m.exons])


@unittest.skipUnless(shutil.which("samtools"), "samtools not on PATH")
class TestPlantingAndReads(unittest.TestCase):
    SPEC = {"name": "t", "seed": 5, "background": {"source": "random", "length": 60000},
            "template": {"source": "synthetic", "exons": EXONS, "introns": INTRONS},
            "copies": [{"id": "A"}, {"id": "A2", "strand": "-", "ops": [{"op": "snp", "rate": 0.02}, {"op": "exon_delete", "exon": 3}]},
                       {"id": "A3", "in_reference": False, "expression": 7, "ops": [{"op": "invert", "exon": 2}]},
                       {"id": "A4", "expression": 0}],
            "layout": {"spacing": 5000, "start": 3000}, "reads": {"per_copy": 10, "jitter": 20}}

    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d)

    def test_truth_coordinates_read_back(self):
        import pysam
        from famsim import chromosome, reads
        planted, man = chromosome.build(self.SPEC, self.d, log=lambda *a: None)
        fa = pysam.FastaFile(os.path.join(self.d, "genome.truth.fa"))
        ref = pysam.FastaFile(os.path.join(self.d, "genome.ref.fa"))
        self.assertNotIn("absent_A3", ref.references); self.assertIn("absent_A3", fa.references)
        for p in planted:
            for s, e, label, in_rna, inv in p.dna_exon_intervals():
                g = fa.fetch(p.contig, s, e).upper()
                ex = next(x for x in p.model.exons if x.label == label)
                self.assertEqual(rc(g) if p.strand == "-" else g, p.model.seq[ex.start:ex.end], f"{p.id} {label}")
        n = reads.simulate(planted, self.SPEC["reads"], self.d, 5, log=lambda *a: None)
        self.assertEqual(n, {"A": 10, "A2": 10, "A3": 7, "A4": 0})
        truth = reads.load_truth(self.d)
        for name, t in truth.items():
            p = next(x for x in planted if x.id == t["copy"])
            self.assertTrue(t["junctions"] <= p.junctions(), name)
            self.assertTrue(name.startswith(t["copy"] + "|"))
        # minus-strand junctions are genomic intervals inside the planted copy
        p2 = next(x for x in planted if x.id == "A2")
        for s, e in p2.junctions():
            self.assertTrue(p2.pos <= s < e <= p2.end)
            self.assertEqual(rc(fa.fetch("sim", s, s + 2).upper()), "AG"); self.assertEqual(rc(fa.fetch("sim", e - 2, e).upper()), "GT")
        # determinism
        d2 = tempfile.mkdtemp()
        try:
            chromosome.build(self.SPEC, d2, log=lambda *a: None); reads.simulate(planted, self.SPEC["reads"], d2, 5, log=lambda *a: None)
            for fn in ("genome.truth.fa", "truth.gtf", "reads.fq"):
                self.assertEqual(open(os.path.join(self.d, fn)).read(), open(os.path.join(d2, fn)).read(), fn)
        finally:
            shutil.rmtree(d2)

    @unittest.skipUnless(shutil.which("minimap2"), "minimap2 not on PATH")
    def test_verify_passes_and_catches_a_corrupted_genome(self):
        from famsim import chromosome, reads, verify
        planted, _ = chromosome.build(self.SPEC, self.d, log=lambda *a: None)
        reads.simulate(planted, self.SPEC["reads"], self.d, 5, log=lambda *a: None)
        R = verify.run(self.d, log=lambda *a: None)
        self.assertEqual(R.failed, [], R.failed)
        # corrupt exon 2 of A in the truth genome: planted_in_genome for A must FAIL, nothing else about A2-A4 changes
        p = next(x for x in planted if x.id == "A"); s, e = p.genomic(p.model.exons[1].start, p.model.exons[1].end)
        contigs = verify.read_fasta(os.path.join(self.d, "genome.truth.fa"))
        contigs["sim"] = contigs["sim"][:s] + rc(contigs["sim"][s:e]) + contigs["sim"][e:]
        chromosome._write_fasta(os.path.join(self.d, "genome.truth.fa"), contigs)
        os.remove(os.path.join(self.d, "genome.truth.fa.fai"))
        R = verify.run(self.d, log=lambda *a: None)
        self.assertTrue(any(r[0] == "planted_in_genome" and r[1] == "A" for r in R.failed))


if __name__ == "__main__":
    unittest.main(verbosity=1)
