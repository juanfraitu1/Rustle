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
