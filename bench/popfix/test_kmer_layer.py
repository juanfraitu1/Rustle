"""Unit tests for the k-mer layer (run: python3 -m pytest -q bench/popfix)."""
import os, sys
import pytest
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


# ---- fix round 1 -------------------------------------------------------------------------------

def _write(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def _psv(tmp_path, copies, pos_rows):
    cat = str(tmp_path / "copies.tsv")
    _write(cat, ["tid", "chrom"], [("t1", "chrY"), ("t2", "chrY")])
    prefix = str(tmp_path / "p")
    _write(prefix + ".psv_copies.tsv", ["family_id", "copy_index", "copy_tid", "alleles", "strand"], copies)
    _write(prefix + ".psv_pos.tsv", ["family_id", "col_index", "copy_index", "genome_pos", "base"], pos_rows)
    return prefix, cat


def _mixed_genome():
    g = list("ACGT" * 60)            # 240 bp
    g[120] = "G"                     # copy 0, '+' strand: frame base G = forward G
    g[130] = "C"                     # copy 1, '-' strand: forward C, frame base complement = G
    return {"chrY": "".join(g)}


def test_mixed_strand_columns_carry_the_frame_allele(tmp_path):
    prefix, cat = _psv(tmp_path,
                       copies=[("F", 0, "t1", "G", "+"), ("F", 1, "t2", "G", "-")],
                       pos_rows=[("F", 0, 0, 120, "G"), ("F", 0, 1, 130, "G")])
    cols = K.columns(prefix, cat, _mixed_genome())
    assert cols[("F", 0)] == [(0, "chrY", 120, "G", "+"), (1, "chrY", 130, "G", "-")]


def test_mixed_strand_pair_with_shared_frame_allele_is_one_allele(tmp_path):
    # Two copies share the frame allele G but have complementary forward bases (G and C).
    # Grouping by the psv base must give ONE allele, so a man carrying the k-mers of copy 0 is FIXED.
    g = _mixed_genome()
    prefix, cat = _psv(tmp_path,
                       copies=[("F", 0, "t1", "G", "+"), ("F", 1, "t2", "G", "-")],
                       pos_rows=[("F", 0, 0, 120, "G"), ("F", 0, 1, 130, "G")])
    plan = K.build_plan(K.columns(prefix, cat, g), g)
    al, fl, lab = plan[("F", 0)]
    assert set(lab.values()) == {"G"}
    pm = set(al[0]) | {x for f in fl.values() for x in f}
    out = K.kmer_status(plan, iter([pm, pm]), 2)
    assert out[("F", 0)]["status"] == "FIXED"
    assert out[("F", 0)]["n_carry_all"] == 2


def test_frame_check_raises_naming_family_col_copy_on_mismatch(tmp_path):
    # forward base at 120 is G; psv base A on '+' matches neither G nor its complement
    prefix, cat = _psv(tmp_path, copies=[("F", 0, "t1", "A", "+")], pos_rows=[("F", 0, 0, 120, "A")])
    with pytest.raises(ValueError, match="family F col 0 copy 0"):
        K.columns(prefix, cat, _mixed_genome())


def test_frame_check_rejects_minus_strand_base_equal_to_forward(tmp_path):
    # '-' strand: psv base must be complement(forward G) = C; base G is wrong
    prefix, cat = _psv(tmp_path, copies=[("F", 0, "t1", "G", "-")], pos_rows=[("F", 0, 0, 120, "G")])
    with pytest.raises(ValueError, match="family F col 0 copy 0"):
        K.columns(prefix, cat, _mixed_genome())


def test_zero_assemblies_refuse_instead_of_fixing(tmp_path):
    with pytest.raises(SystemExit):
        K.find_men(str(tmp_path))


def test_find_men_returns_sorted_assemblies(tmp_path):
    for n in ("b", "a"):
        (tmp_path / f"{n}.chrY.fasta").write_text(">x\nA\n")
    assert K.find_men(str(tmp_path)) == [str(tmp_path / "a.chrY.fasta"), str(tmp_path / "b.chrY.fasta")]


def test_zero_men_is_uncovered_never_fixed():
    plan = {("F", 0): ({0: ["a"]}, {0: ("l", "r")}, {0: "G"})}
    out = K.kmer_status(plan, iter([]), 0)
    assert out[("F", 0)]["status"] == "UNCOVERED"
    assert K.decide({"G"}, [], [], 0) == "UNCOVERED"


def test_not_assessed_column_is_passed_through():
    out = K.kmer_status({("F", 1): None}, iter([]), 0)
    assert out[("F", 1)]["status"] == "NOT_ASSESSED"


PLAN = {("F", 0): ({0: ["a"]}, {0: ("l", "r")}, {0: "G"})}
COV = {"a", "l", "r"}


def test_boundary_exactly_ceil_half_covering_is_eligible():
    out = K.kmer_status(PLAN, iter([COV, COV, set(), set()]), 4)     # 2 of 4 = ceil(4/2)
    assert out[("F", 0)]["status"] == "FIXED"
    assert out[("F", 0)]["n_cover"] == 2


def test_fewer_than_ceil_half_covering_is_uncovered():
    out = K.kmer_status(PLAN, iter([COV, COV, set(), set(), set()]), 5)  # 2 of 5 < ceil(5/2)=3
    assert out[("F", 0)]["status"] == "UNCOVERED"


def test_presence_count_must_match_n_men():
    with pytest.raises(ValueError):
        K.kmer_status(PLAN, iter([COV]), 2)


# ---- fix round 2: pin the presence -> carry -> status path and n_carry_all -------------------------

TWO = {("F", 0): ({0: ["a0"], 1: ["b0"]}, {0: ("l0", "r0"), 1: ("l1", "r1")}, {0: "A", 1: "B"})}
FLANKS = {"l0", "r0", "l1", "r1"}
BOTH = {"a0", "b0"} | FLANKS


def test_two_allele_polymorphic_when_one_covering_man_lacks_allele_b():
    no_b = {"a0"} | FLANKS                      # covering man whose presence set lacks allele B's k-mers
    out = K.kmer_status(TWO, iter([BOTH, BOTH, no_b, BOTH]), 4)
    r = out[("F", 0)]
    assert r["status"] == "POLYMORPHIC"
    assert r["n_cover"] == 4
    assert r["n_carry_all"] == 3                # exactly the covering men that carry both alleles


def test_two_allele_fixed_when_every_covering_man_carries_both():
    out = K.kmer_status(TWO, iter([BOTH, BOTH, BOTH]), 3)
    r = out[("F", 0)]
    assert r["status"] == "FIXED"
    assert r["n_cover"] == 3
    assert r["n_carry_all"] == r["n_cover"] == 3


def test_man_with_only_one_of_two_flank_kmers_is_not_covering():
    half = {"a0", "b0", "l0", "l1", "r1"}       # site 0's right flank r0 missing
    out = K.kmer_status(TWO, iter([BOTH, BOTH, BOTH, half]), 4)
    r = out[("F", 0)]
    assert r["n_cover"] == 3                    # the half-flank man does not cover the column
    assert r["n_carry_all"] == 3                # and is excluded from carry counting
    assert r["status"] == "FIXED"


def _mixed_random_genome():
    import random
    rnd = random.Random(20261012)               # non-periodic, so the two sites' k-mer windows are distinct
    g = [rnd.choice("ACGT") for _ in range(240)]
    g[120] = "G"                                # copy 0, '+': frame base G
    g[130] = "C"                                # copy 1, '-': forward C, frame base complement = G
    return {"chrY": "".join(g)}


def test_mixed_strand_grouping_uses_psv_base_not_forward_base(tmp_path):
    # psv-base grouping: one allele G. Each man carries only one copy's exclusive k-mers, yet carries G -> FIXED.
    # Forward-base grouping would give alleles {G, C}, each man carries one -> POLYMORPHIC.
    g = _mixed_random_genome()
    prefix, cat = _psv(tmp_path,
                       copies=[("F", 0, "t1", "G", "+"), ("F", 1, "t2", "G", "-")],
                       pos_rows=[("F", 0, 0, 120, "G"), ("F", 0, 1, 130, "G")])
    plan = K.build_plan(K.columns(prefix, cat, g), g)
    al, fl, lab = plan[("F", 0)]
    assert set(lab.values()) == {"G"}
    flanks = {x for f in fl.values() for x in f}
    only0 = set(al[0]) - set(al[1])
    only1 = set(al[1]) - set(al[0])
    assert only0 and only1
    out = K.kmer_status(plan, iter([only0 | flanks, only1 | flanks]), 2)
    assert out[("F", 0)]["status"] == "FIXED"
    assert out[("F", 0)]["n_carry_all"] == 2
