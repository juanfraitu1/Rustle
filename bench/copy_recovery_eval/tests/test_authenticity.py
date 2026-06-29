import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import pysam
import importlib.util

def _load_guard():
    path = os.path.join(os.path.dirname(__file__), "..", "50_authenticity_guard.py")
    spec = importlib.util.spec_from_file_location("guard", path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m

def _make_bam(tmp):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 10000}]}
    bam = os.path.join(tmp, "toy.bam")
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        a = pysam.AlignedSegment(out.header)
        a.query_name = "r1"; a.flag = 0; a.reference_id = 0
        a.reference_start = 1000; a.mapping_quality = 60
        a.cigartuples = [(0, 100)]; a.query_sequence = "A"*100
        a.query_qualities = [30]*100
        out.write(a)
        b = pysam.AlignedSegment(out.header)
        b.query_name = "r2"; b.flag = 0x100; b.reference_id = 0
        b.reference_start = 1010; b.mapping_quality = 0
        b.cigartuples = [(0, 100)]; b.query_sequence = "A"*100
        b.query_qualities = [30]*100
        out.write(b)
    pysam.index(bam)
    return bam

def test_primary_support_counts_only_primary(tmp_path):
    guard = _load_guard()
    bam = _make_bam(str(tmp_path))
    n = guard.primary_support(bam, "chr1", [(1000, 1100)])
    assert n == 1
    assert guard.primary_support(bam, "chr1", [(5000, 5100)]) == 0


def _aln(out, name, flag, ref_id, start, AS, mapq):
    a = pysam.AlignedSegment(out.header)
    a.query_name = name; a.flag = flag; a.reference_id = ref_id
    a.reference_start = start; a.mapping_quality = mapq
    a.cigartuples = [(0, 100)]; a.query_sequence = "A" * 100
    a.query_qualities = [30] * 100
    a.set_tag("AS", AS, value_type="i")
    out.write(a)


def _make_decisive_bam(tmp):
    # two copies: copyA = chrA:1000-1100, sister copyB = chrB:1000-1100
    header = {"HD": {"VN": "1.6"},
              "SQ": [{"SN": "chrA", "LN": 5000}, {"SN": "chrB", "LN": 5000}]}
    bam = os.path.join(tmp, "dec.bam")
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        # write coordinate-sorted: all chrA (ref_id 0) then all chrB (ref_id 1).
        # chrA placements
        _aln(out, "r_own",     0,     0, 1000, 100, 60)  # better on A -> decisive for A
        _aln(out, "r_tied",    0,     0, 1010, 100, 0)   # equal on A and B -> tied
        _aln(out, "r_foreign", 0x100, 0, 1020, 80,  0)   # worse on A -> sister-owned
        # chrB (sister) placements
        _aln(out, "r_own",     0x100, 1, 1000, 90,  0)
        _aln(out, "r_tied",    0x100, 1, 1010, 100, 0)
        _aln(out, "r_foreign", 0,     1, 1020, 100, 60)
    pysam.index(bam)
    return bam


def test_decisive_support_splits_own_tied_foreign(tmp_path):
    guard = _load_guard()
    bam_path = _make_decisive_bam(str(tmp_path))
    bam = pysam.AlignmentFile(bam_path, "rb")
    d = guard.decisive_support(bam, "chrA", [(1000, 1100)], sisters=[("chrB", [(1000, 1100)])])
    assert d["n_decisive"] == 1   # r_own
    assert d["n_tied"] == 1       # r_tied
    assert d["n_foreign"] == 1    # r_foreign
    assert d["n_total"] == 3


def test_decisive_support_unique_copy_is_decisive(tmp_path):
    guard = _load_guard()
    bam_path = _make_decisive_bam(str(tmp_path))
    bam = pysam.AlignmentFile(bam_path, "rb")
    # no sisters supplied -> every read on the copy maps uniquely within the family
    d = guard.decisive_support(bam, "chrA", [(1000, 1100)], sisters=[])
    assert d["n_decisive"] == 3 and d["n_tied"] == 0 and d["n_foreign"] == 0
