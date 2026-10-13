"""Unit tests for the k-mer layer (run: python3 -m pytest -q bench/popfix)."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import kmer_layer as K

G = {"chrY": "A" * 100 + "ACGTTGCATGCAAGTCCGATAGCTAGGCTTACGATCGTAGCTAGCTTAGCAT" + "C" * 100}


def test_allele_kmers_are_the_31_windows_containing_the_site():
    ks = K.allele_kmers(G, "chrY", 120)
    assert len(ks) == 31
    assert all(len(k) == 31 for k in ks)
    raw = G["chrY"][90:151]
    assert {K.canon(raw[i:i + 31]) for i in range(31)} == set(ks)


def test_flank_kmers_sit_40_bp_from_the_site():
    left, right = K.flank_kmers(G, "chrY", 120)
    assert left == K.canon(G["chrY"][50:81])      # ends at 120 - 40
    assert right == K.canon(G["chrY"][160:191])   # starts at 120 + 40


def test_fixed_needs_half_the_men_covering_and_every_allele_carried():
    sites = {"A": ["a1"], "C": ["c1"]}           # allele -> its sites' allele k-mer sets (ids)
    cover = [True, True, True, False]            # man covers the column
    carry = [{"A", "C"}, {"A", "C"}, {"A", "C"}, set()]
    assert K.decide(sites, cover, carry, n_men=4) == "FIXED"
    carry[1] = {"A"}                             # one covering man lost allele C (homogenised)
    assert K.decide(sites, cover, carry, n_men=4) == "POLYMORPHIC"
    assert K.decide(sites, [True, False, False, False], carry, n_men=4) == "UNCOVERED"
