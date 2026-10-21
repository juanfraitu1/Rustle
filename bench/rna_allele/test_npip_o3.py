"""PREREG_o3_npip_maternal_copy_2026-10-21: the driver's own pieces (the net's orientation, the target reads, the bars)."""
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("rna_allele_npip_o3", os.path.join(HERE, "npip_o3.py"))
N = importlib.util.module_from_spec(spec)
spec.loader.exec_module(N)

T = ("CM054600.2", 6720000, 6820000)


def test_a_net_read_is_written_in_sequencing_orientation():
    assert N.sequencing_orientation("AACG", False) == "AACG"
    assert N.sequencing_orientation("AACG", True) == "CGTT"


def test_a_target_read_fits_the_target_best_by_five_ms_and_a_long_intron_alignment_does_not_count():
    a = lambda acc, s, e, ms, intron=0: (acc, s, e, ms, intron)
    assert N.in_target([a("CM054600.2", 6750000, 6753000, 1500), a("CM054600.2", 17400000, 17403000, 1494)], T) is True
    assert N.in_target([a("CM054600.2", 6750000, 6753000, 1500), a("CM054570.2", 100, 3000, 1497)], T) is False
    assert N.in_target([a("CM054600.2", 6750000, 6753000, 1500), a("CM054570.2", 100, 300000, 1600, 150000)], T) is True
    assert N.in_target([a("CM054582.2", 1, 2, 10)], T) is False
    assert N.in_target([a("CM054600.2", 6819000, 6822000, 900)], T) is True               # overlapping the target's end


def test_a_haplotype_hit_lies_in_the_target_when_it_overlaps_it_on_the_same_contig():
    assert N.hit_in_target("mat:CM054600.2:6750123-6753000 0.9995", T) is True
    assert N.hit_in_target("mat:CM054600.2:17395900-17400000 1.0000", T) is False
    assert N.hit_in_target("pat:CM054576.2:6750123-6753000 1.0000", T) is False
    assert N.hit_in_target("", T) is False


def _c(name, n, cls, in_t, reads):
    return dict(candidate=name, n_contigs=n, cls=cls, in_T=in_t, reads=set(reads))


def test_the_bars_need_one_counted_target_candidate_holding_and_pure_and_no_counted_false_flag():
    good = [_c("NPIP:a", 3, "a_haplotype_only", True, "r1 r2 r3 r4 r5 r6 r7 r8 r9 r10".split()),
            _c("NPIP:b", 1, "c_unmatched", False, ["x"]),                                  # one contig: listed, not counted
            _c("NPIP:c", 2, "a_haplotype_only", False, ["y"])]                             # another true flag: reported only
    rt = set("r1 r2 r3 r4 r5 r6 r7 r8 r9 r10 r11 r12".split())
    b = N.bars(good, rt)
    assert {k: b[k] for k in ("E1", "E2", "E3", "E4", "E5", "verdict")} == dict(E1=True, E2=True, E3=True, E4=True, E5=True,
                                                                                 verdict="PASS")
    assert b["held"] == 10 and b["target_reads"] == 12 and b["other_true"] == ["NPIP:c"] and b["false"] == []
    assert N.bars(good + [_c("NPIP:d", 2, "b_allele", False, ["z"])], rt)["E5"] is False
    assert N.bars(good + [_c("NPIP:e", 2, "a_haplotype_only", True, ["r11"])], rt)["E2"] is False
    impure = [_c("NPIP:a", 3, "a_haplotype_only", True, "r1 r2 r3 r4 r5 r6 r7 x y z".split())]
    assert N.bars(impure, rt)["E3"] is True and N.bars(impure, rt)["E4"] is False                  # holds 7 of 12; 7 of its 10 reads
    assert N.bars([], rt)["verdict"] == "FAIL"
