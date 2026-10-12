"""Unit tests for the coordinate helpers of psv_sun.py (run: python3 -m pytest -q bench/psv_sun)."""
import os, random, subprocess, sys, tempfile
sys.path.insert(0, os.path.dirname(__file__))
import psv_sun as P

EX = [(100, 110), (200, 205), (300, 320)]  # 0-based half-open, L = 35


def test_tx_offset_genome_pos_roundtrip_both_strands():
    for strand in "+-":
        seen = set()
        for s, e in EX:
            for g in range(s, e):
                o = P.tx_offset(EX, strand, g)
                assert P.genome_pos(EX, strand, o) == g
                seen.add(o)
        assert seen == set(range(35))
    assert P.tx_offset(EX, "+", 150) is None
    assert P.tx_offset(EX, "+", 100) == 0 and P.tx_offset(EX, "-", 100) == 34


def test_cs_pairs_sub_ins_del():
    # target ACGTACGT, query ACcTAGGACGT? keep simple: 3 match, sub, ins 'gg', del 'a', 2 match
    pairs = list(P.cs_pairs(":3*gc+gg-a:2", 10, 0))
    assert pairs == [(0, 10), (1, 11), (2, 12), (3, 13), (6, 15), (7, 16)]


def test_canon():
    assert P.canon("TTT") == "AAA" and P.canon("ACG") == "ACG" and P.canon("CGT") == "ACG"


def _paf(t, q):
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/t.fa", "w").write(f">t\n{t}\n")
        open(f"{d}/q.fa", "w").write(f">q\n{q}\n")
        return P.run_paf(["minimap2", "-cx", "asm5", "--cs", f"{d}/t.fa", f"{d}/q.fa"])


def test_paf_qmap_plus_and_minus_strand_with_minimap2():
    rng = random.Random(3)
    t = "".join(rng.choice("ACGT") for _ in range(6000))
    q = t[1000:4000]
    q = q[:1500] + ("A" if q[1500] != "A" else "C") + q[1501:]  # one substitution
    for strand, query in (("+", q), ("-", P.rc(q))):
        r = [x for x in _paf(t, query) if x["tp"] == "P"][0]
        assert r["strand"] == strand
        m = P.paf_qmap(r)
        for qf in (0, 700, 1500, 2999):
            # forward-query offset qf of q corresponds to target 1000 + qf
            off = qf if strand == "+" else len(q) - 1 - qf
            assert m[off] == 1000 + qf


def test_fasta_fetch_against_samtools():
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/g.fa", "w").write(">c1\nACGTACGTAC\nGGTTAACC\n>c2\nTTTT\n")
        subprocess.run(["samtools", "faidx", f"{d}/g.fa"], check=True)
        F = P.Fasta(f"{d}/g.fa")
        assert F.fetch("c1", 8, 12) == "ACGG" and F.fetch("c2", 0, 4) == "TTTT" and F.fetch("c1", 17, 30) == "C"
