import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import copyclusters as CC  # noqa: E402


def hap(n, own=0, other=0, elsewhere=0, nomatch=0, nohit=0):
    return {"n_reads": str(n), "own_family": str(own), "other_family": str(other), "elsewhere": str(elsewhere),
            "no_match_0.995": str(nomatch), "no_hit": str(nohit)}


def site(unexpl, uncovered, share):
    return {"n_unexpl": str(unexpl), "unexpl_uncovered": str(uncovered), "median_read_share": str(share)}


def test_classify_rules_in_order():
    assert CC.classify(site(10, 0, 1.0), hap(10, own=2, other=7, nomatch=1)) == "other_family_copy"
    assert CC.classify(site(10, 0, 1.0), hap(10, own=2, nomatch=8)) == "not_in_genome"
    assert CC.classify(site(10, 0, 1.0), hap(10, own=2, elsewhere=6, nomatch=2)) == "absent_from_ref"
    assert CC.classify(site(10, 6, 1.0), hap(10, own=9, nomatch=1)) == "catalog_lacks_exons"
    assert CC.classify(site(10, 2, 0.05), hap(10, own=9, nomatch=1)) == "minority_segment"
    assert CC.classify(site(10, 2, 0.9), hap(10, own=9, nomatch=1)) == "own_copy_variant"
