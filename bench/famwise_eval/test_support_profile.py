"""tests of support_profile.py"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import support_profile as S

T = dict(chrom="c", strand="+", exons=[(0, 100), (200, 300), (400, 500)])
TINT = [(100, 200), (300, 400)]

def read(*blocks):
    b = list(blocks)
    return b, [(x[1], y[0]) for x, y in zip(b, b[1:])]

def test_classify():
    assert S.classify(*read((10, 100), (200, 300), (400, 480)), TINT, 3, 20)[0] == "compatible"
    assert S.classify(*read((10, 102), (201, 300)), TINT, 3, 20)[0] == "compatible"          # within jitter (profile then counts it as shifted)
    assert S.classify(*read((10, 108), (200, 300)), TINT, 3, 20)[0] == "near_miss"           # donor 8 bp off
    assert S.classify(*read((10, 100), (400, 480)), TINT, 3, 20)[0] == "skip"                # exon 2 skipped
    assert S.classify(*read((10, 60), (80, 100)), TINT, 3, 20)[0] == "other_junction"
    assert S.classify(*read((10, 250)), TINT, 3, 20)[0] == "retention"                       # block runs through intron 1
    assert S.classify(*read((10, 108), (200, 300)), TINT, 3, 20, frozenset({(108, 200)}))[0] == "alt_site"

def test_profile_flags_and_cover():
    reads = [read((0, 100), (200, 300), (400, 500))] * 3 + [read((10, 108), (200, 300))] * 5
    p = S.profile(T, reads, 3, 20, 50, 3)
    assert p["compatible"] == 3 and p["full_chain"] == 3 and p["five_complete"] == 3 and p["three_complete"] == 3
    assert p["cover1"] == 1.0 and p["coverk"] == 1.0 and p["j_near"] == [5, 0]
    assert "NEAR_MISS_DOMINANT" in p["flags"] and "WEAK_JUNCTION" not in p["flags"]
    p2 = S.profile(T, [read((0, 100), (200, 300))] * 2, 3, 20, 50, 3)
    assert p2["flags"] == ["WEAK_JUNCTION", "UNCOVERED", "NO_FULL_READ"] and p2["j_exact"] == [2, 0]

def test_shifted_reads_do_not_vouch():
    # 20 reads 2 bp off at junction 1, 3 exact: the 2 bp-shifted reads are not compatible; the junction is SHIFTED_DOMINANT and THIN_EXACT
    reads = [read((0, 102), (202, 300), (400, 500))] * 20 + [read((0, 100), (200, 300), (400, 500))] * 3
    p = S.profile(T, reads, 3, 20, 50, 3)
    assert p["compatible"] == 3 and p["why"]["shifted"] == 20 and p["j_exact"] == [3, 3] and p["j_jit"] == [20, 0]
    assert "SHIFTED_DOMINANT" in p["flags"] and "THIN_EXACT" not in p["flags"]          # 3 vs 3: not thin relative to the strongest
    p3 = S.profile(T, [read((0, 100), (200, 300), (400, 500))] * 2 + [read((200, 300), (400, 500))] * 40, 3, 20, 50, 3)
    assert p3["j_exact"] == [2, 42] and "THIN_EXACT" in p3["flags"]                         # 2 < 10% of 42

def test_minus_strand_ends():
    t = dict(T, strand="-")
    p = S.profile(t, [read((0, 100), (200, 300), (400, 420))], 3, 20, 50, 1)
    assert p["three_complete"] == 1 and p["five_complete"] == 0                              # 5' end of a minus transcript is on the right
