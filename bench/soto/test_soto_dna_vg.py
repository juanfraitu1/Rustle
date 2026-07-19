#!/usr/bin/env python3
"""Tests for the DNA variation-graph ceiling demo (bench/soto/soto_dna_vg.py).
Run: PYTHONHASHSEED=0 pytest bench/soto/test_soto_dna_vg.py -v"""
import os
import sys
sys.path.insert(0, "bench/soto")
import soto_dna_vg as V


def test_parse_family_members():
    bed = [
        "chr1\t100\t200\tSRGAP2C|ID_462\t0\t.",
        "chr7\t68507281\t68522154\tPMS2P4|ID_8\t0\t.",
        "chr1\t300\t400\tSRGAP2|ID_462\t0\t.",
    ]
    m = V.parse_family_members(bed, "ID_462")
    assert m == [("SRGAP2C", "chr1", 100, 200), ("SRGAP2", "chr1", 300, 400)]


def test_member_seq_plus_one_offset():
    fa = {"chr1:101-200": "ACGT", "chr7:1-9": "TTTT"}
    # BED start 100 (0-based) -> header start 101 (1-based)
    assert V.member_seq(fa, "chr1", 100, 200) == "ACGT"


def test_read_fasta(tmp_path):
    p = tmp_path / "m.fa"
    p.write_text(">chr1:101-200\nAC\nGT\n>chr7:1-9\nTTTT\n")
    assert V.read_fasta(str(p)) == {"chr1:101-200": "ACGT", "chr7:1-9": "TTTT"}
