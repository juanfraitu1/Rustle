"""Tests of bench/stgaps/famasm.py (plan docs/superpowers/plans/2026-10-31-multicopy-assembly.md)."""
import os
import sys

import pysam
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famasm  # noqa: E402


def test_good_secondary_threshold():
    assert famasm.good_secondary(98, 100)
    assert not famasm.good_secondary(97, 100)


def _bam(path, header_len=10_000):
    hdr = {"HD": {"VN": "1.6"}, "SQ": [{"LN": header_len, "SN": "c1"}]}
    return pysam.AlignmentFile(path, "wb", header=hdr)


def _rec(name, start, cigar, flag, as_):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.flag = flag
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = 60
    a.cigarstring = cigar
    a.query_sequence = "A" * sum(n for op, n in a.cigartuples if op in (0, 1, 4, 7, 8))
    a.set_tag("AS", as_)
    return a


def test_pool_keeps_one_record_per_copy_and_drops_bad_secondaries(tmp_path):
    # copies: F1 copy 1 = c1:1000-2000, copy 2 = c1:5000-6000
    copies = tmp_path / "copies.tsv"
    copies.write_text("family_id\tcopy_idx\tchrom\tlocus_start\tlocus_end\tn_reads\n"
                      "F1\t1\tc1\t1000\t2000\t10\nF1\t2\tc1\t5000\t6000\t3\n")
    bam = tmp_path / "r.bam"
    with _bam(str(bam)) as w:
        w.write(_rec("r1", 1100, "100M300N100M", 0, 100))      # primary on copy 1
        w.write(_rec("r1", 5100, "100M300N100M", 256, 99))     # good secondary on copy 2
        w.write(_rec("r2", 5100, "200M", 256, 80))             # bad secondary (80 < 0.98 x 100)
        w.write(_rec("r2", 1100, "200M", 0, 100))              # r2 primary on copy 1, unspliced
        w.write(_rec("r3", 1100, "200M", 2048, 100))           # supplementary: dropped
    pysam.sort("-o", str(tmp_path / "s.bam"), str(bam))
    pysam.index(str(tmp_path / "s.bam"))
    mols = tmp_path / "molecules.tsv"
    mols.write_text("#as_table\tbam=x\nname\tbest_as\tsecond_as\tn_records\tprimary_chrom\tprimary_as\n"
                    "r1\t100\t99\t2\tc1\t100\nr2\t100\t80\t2\tc1\t100\nr3\t100\t0\t1\tc1\t100\n")
    out = tmp_path / "out"
    famasm.pool(str(copies), str(tmp_path / "s.bam"), str(mols), str(out))
    rows = [l.rstrip("\n").split("\t") for l in open(out / "pool" / "F1.tsv")]
    assert rows[0] == ["read", "copy_idx", "is_primary", "AS", "chain"]
    body = sorted(rows[1:])
    assert body == [["r1", "1", "1", "100", "1200-1500"], ["r1", "2", "0", "99", "5200-5500"], ["r2", "1", "1", "100", ""]]
    stats = [l.rstrip("\n").split("\t") for l in open(out / "pool.stats.tsv")]
    assert stats[0] == ["family", "molecules", "records", "secondary_kept"]
    assert stats[1] == ["F1", "2", "3", "1"]
