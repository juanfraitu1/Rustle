#!/usr/bin/env python3
"""Run: pixi run -e bench python3 -B -m unittest bench/famhmm/test_famhmm.py"""
import io
import os
import random
import shutil
import sys
import tempfile
import unittest
from contextlib import redirect_stderr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famhmm as F  # noqa: E402

ORFS = [
    ("orf1 source=c1 coords=5..124 length=40 frame=1 desc=x", "M" * 40),
    ("orf2 source=c1 coords=10..189 length=60 frame=2 desc=x", "M" * 60),
]


HAVE_HMMER = shutil.which("hmmbuild") is not None

CODONS = ["GCT", "GAA", "AAA", "CTG", "GAT", "TTC", "GGT", "CAC", "ATC", "AAG", "ATG", "AAC", "CCA", "CAG", "CGT", "TCT", "ACC", "GTG", "TGG", "TAC"]


def toy_family(tmp, n_orf=3, seed=1):
    """3-copy toy family: an ATG + 60 sense codons + TAA ORF at offset 0 (frame 1) in a 300-nt sequence; copies 2 and 3 carry
    2 synonymous-position-free mutations outside the ORF; copies beyond n_orf get the ORF scrambled to kill it.
    Returns (prefix, copies, loci) where loci keys are header field 3."""
    rng = random.Random(seed)
    orf = "ATG" + "".join(rng.choice(CODONS) for _ in range(60)) + "TAA"
    tail = "".join(rng.choice("ACGT") for _ in range(300 - len(orf)))
    base = orf + tail
    copies = []
    for i in range(3):
        s = list(base)
        for _ in range(2 * i):
            j = rng.randrange(len(orf) + 5, len(base))
            s[j] = "ACGT"[("ACGT".index(s[j]) + 1) % 4]
        if i >= n_orf:                       # kill the ORF: a stop every 10 codons in every frame
            for j in range(3, len(orf) - 3, 10):
                s[j:j + 3] = list("TAA")
        copies.append(("FAM1|%d|chr1:%d-%d|+|nexon=1" % (i, 1000 * i + 1, 1000 * i + 300), "".join(s)))
    prefix = os.path.join(tmp, "toy")
    with open(prefix + ".fam.copies.fa", "w") as fh:
        for h, s in copies:
            fh.write(">%s\n%s\n" % (h, s))
    with open(prefix + ".fam.copies.tsv", "w") as fh:
        fh.write("family_id\tcopy_idx\ttid\n" + "".join("FAM1\t%d\tDN%d\n" % (i, i) for i in range(3)))
    return prefix, copies


@unittest.skipUnless(HAVE_HMMER, "hmmbuild not on PATH")
class Build(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_build_toy_family(self):
        prefix, _ = toy_family(self.tmp)
        out = os.path.join(self.tmp, "out")
        rows = F.build(prefix, out)
        self.assertEqual(rows, [{"family": "FAM1", "n_copies": 3, "prot_status": "ok"}])
        for f in ("dna/FAM1.hmm", "prot/FAM1.hmm", "families.dna.hmm.h3m", "families.prot.hmm.h3m", "build.tsv"):
            self.assertTrue(os.path.exists(os.path.join(out, f)), f)
        orf = [ln.split("\t") for ln in open(os.path.join(out, "orf.tsv")).read().splitlines()[1:]]
        self.assertEqual([(r[2], r[3]) for r in orf], [("1", "0")] * 3)

    def test_build_prot_none_when_single_orf(self):
        prefix, _ = toy_family(self.tmp, n_orf=1)
        out = os.path.join(self.tmp, "out")
        rows = F.build(prefix, out)
        self.assertEqual(rows[0]["prot_status"], "none")
        self.assertTrue(os.path.exists(os.path.join(out, "dna/FAM1.hmm")))
        self.assertFalse(os.path.exists(os.path.join(out, "prot/FAM1.hmm")))

    def test_build_only_missing_family_warns(self):
        prefix, _ = toy_family(self.tmp)
        err = io.StringIO()
        with redirect_stderr(err):
            rows = F.build(prefix, os.path.join(self.tmp, "out"), only={"NOPE"})
        self.assertEqual(rows, [])
        self.assertIn("skip NOPE", err.getvalue())


DOMTBL = [
    "# target name accession tlen query name accession qlen E-value score bias # of c-Evalue i-Evalue score bias from to from to from to acc description",
    "Kinase   PF00069.1  260 FAM1  -  400  1e-50  170.0 0.1  1  2  1e-40  1e-38  80.0 0.0  1 100  120 220  118 222  0.9  -",
    "SH2      PF00017.1   80 FAM1  -  400  1e-20   70.0 0.0  2  2  1e-18  1e-17  60.0 0.0  1  80    5   85    3   87  0.9  -",
    "Junk     PF99999.1   50 FAM1  -  400  1e-02    5.0 0.0  1  1  1e-02  1e-02   5.0 0.0  1  50  300  350  298  352  0.5  -",
]


class Scan(unittest.TestCase):
    def test_parse_domtblout_orders_domains_by_position(self):
        d = F.parse_domtblout(DOMTBL)
        self.assertEqual([x[1] for x in d["FAM1"]], ["SH2", "Kinase"])
        self.assertEqual(d["FAM1"][0][0], 5)


class Helpers(unittest.TestCase):
    def test_copies_by_family_groups_on_first_field(self):
        recs = [("MCL0|0|chr1:1-9|+|nexon=1", "A"), ("MCL1|0|chr2:1-9|+|nexon=1", "C"), ("MCL0|1|chr1:20-29|-|nexon=1", "G")]
        fams = F.copies_by_family_records(recs)
        self.assertEqual(sorted(fams), ["MCL0", "MCL1"])
        self.assertEqual(len(fams["MCL0"]), 2)

    def test_longest_orf_picks_longest_and_parses_coords(self):
        hdr, frame, off, aa = F.longest_orf(ORFS)
        self.assertTrue(hdr.startswith("orf2"))
        self.assertEqual((frame, off, aa), (2, 9, "M" * 60))

    def test_longest_orf_none_when_empty(self):
        self.assertIsNone(F.longest_orf([]))


if __name__ == "__main__":
    unittest.main()
