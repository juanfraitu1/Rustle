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


def test_msa_to_gfa_snp_bubble():
    # col0 "A" invariant; col1 variant (C/G/C); col2-3 "GT" invariant
    rows = ["ACGT", "AGGT", "ACGT"]
    gfa, paths = V.msa_to_gfa(rows, ["m1", "m2", "m3"])
    assert paths["m1"] == ["1", "2", "4"]   # A , C , GT
    assert paths["m2"] == ["1", "3", "4"]   # A , G , GT
    assert paths["m3"] == ["1", "2", "4"]
    assert "S\t2\tC" in gfa and "S\t3\tG" in gfa
    assert "S\t1\tA" in gfa and "S\t4\tGT" in gfa
    assert "P\tm2\t1+,3+,4+\t*" in gfa
    # link m2 traverses 1->3 and 3->4
    links = {(l.split("\t")[1], l.split("\t")[3]) for l in gfa.splitlines() if l.startswith("L\t")}
    assert ("1", "3") in links and ("3", "4") in links


def test_msa_to_gfa_indel_skips_gap_member():
    # col1 gap in m2 -> m2 skips that allele node
    rows = ["ACGT", "A-GT", "ACGT"]
    gfa, paths = V.msa_to_gfa(rows, ["m1", "m2", "m3"])
    assert paths["m1"] == ["1", "2", "3"]
    assert paths["m2"] == ["1", "3"]        # skips the gap region
    assert paths["m3"] == ["1", "2", "3"]
    links = {(l.split("\t")[1], l.split("\t")[3]) for l in gfa.splitlines() if l.startswith("L\t")}
    assert ("1", "3") in links              # m2's skip link


def test_msa_to_gfa_all_identical_single_node():
    rows = ["ACGT", "ACGT"]
    gfa, paths = V.msa_to_gfa(rows, ["a", "b"])
    assert paths["a"] == ["1"] and paths["b"] == ["1"]
    assert "S\t1\tACGT" in gfa
    assert not any(l.startswith("L\t") for l in gfa.splitlines())  # no links, one node
