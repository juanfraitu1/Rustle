import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dup_evidence as de  # noqa: E402

SEDEF_GGO = "\t".join(["NC_011120.1", "1911", "7456", "NC_073224.2", "125568083", "125573628", "S", "5.2", "+", "-",
                       "5545", "5547", "m=5.2;g=0.1", "2", "2", "5543", "5256", "287", "259", "28", "0.948223",
                       "0.947539", "0.053651", "0.0543876", "4", "5545", "5301", "5023", "5256", "287", "4", "4",
                       "673M1I2734M1I328M1D1289M1D519M", "0.947539"])
BISER = "\t".join(["chrA", "99997", "120000", "chrB", "59997", "80000", "toy:toy", "3.7", "+", "+", "20003", "20003",
                   "20003M", "X=3.7;ID=0"])


def test_from_sedef_gorilla():
    p = de.from_sedef(SEDEF_GGO, "gorilla")
    assert p[:8] == ("NC_011120.1", 1911, 7456, "NC_073224.2", 125568083, 125573628, "+", "-")
    assert abs(p[8] - 0.948223) < 1e-9 and p[9] == "673M1I2734M1I328M1D1289M1D519M" and p[10] == "D1:sedef"


def test_from_biser_identity_from_error():
    p = de.from_biser(BISER)
    assert p[:8] == ("chrA", 99997, 120000, "chrB", 59997, 80000, "+", "+")
    assert abs(p[8] - 0.963) < 1e-9 and p[9] == "20003M" and p[10] == "D2:biser"


def test_normalize_cigar_merges_eq_x():
    assert de.normalize_cigar("10=2X5=3I4=1D") == "17M3I4M1D"


def paf(q, qlen, qs, qe, strand, t, tlen, ts, te, nm, bl, cg):
    return "\t".join(map(str, [q, qlen, qs, qe, strand, t, tlen, ts, te, nm, bl, 60, f"cg:Z:{cg}"]))


def test_selfpaf_side_a_is_target_and_filters():
    # chunk chr1@1000 query 0-2000 aligned to chr2 5000-7000, identity 0.95, canonical (chr1 > chr2 so target is side A)
    p = de.from_selfpaf(paf("chr2@1000", 5000, 0, 2000, "+", "chr1", 9000, 5000, 7000, 1900, 2000, "2000M"))
    assert p == ("chr1", 5000, 7000, "chr2", 1000, 3000, "+", "+", 0.95, "2000M", "D3:selfaln")
    # self hit (same chrom, overlapping) dropped; short dropped; low identity dropped
    assert de.from_selfpaf(paf("chr1@0", 9000, 5000, 7000, "+", "chr1", 9000, 5000, 7000, 2000, 2000, "2000M")) is None
    assert de.from_selfpaf(paf("chr2@0", 5000, 0, 900, "+", "chr1", 9000, 0, 900, 900, 900, "900M")) is None
    assert de.from_selfpaf(paf("chr2@0", 5000, 0, 2000, "+", "chr1", 9000, 0, 2000, 1700, 2000, "2000M")) is None
    # non-canonical orientation (query side sorts first) dropped: the mirror record carries the pair
    assert de.from_selfpaf(paf("chr1@0", 9000, 5000, 7000, "+", "chr2", 5000, 1000, 3000, 1900, 2000, "2000M")) is None


def test_to_sedef_gorilla_roundtrip_identity_check():
    p = ("chrA", 0, 1000, "chrB", 10, 1010, "+", "-", 0.9123, "1000M", "D2:biser")
    f = de.to_sedef_gorilla(p).split("\t")
    m, mm, frac = float(f[16]), float(f[17]), float(f[20])
    assert abs(m / (m + mm) - frac) <= 1e-4 and f[32] == "1000M" and f[8] == "+" and f[9] == "-"
    assert de.from_sedef("\t".join(f), "gorilla")[:8] == p[:8]


def test_pairs_file_roundtrip(tmp_path):
    p = ("chrA", 0, 1000, "chrB", 10, 1010, "+", "-", 0.9123, "1000M", "D2:biser")
    de.write_pairs([p], tmp_path / "x.pairs.tsv")
    assert de.read_pairs(tmp_path / "x.pairs.tsv") == [p]


def test_from_selfpaf_malformed_lines():
    # short line with fewer than 12 tab fields
    assert de.from_selfpaf("chr1\t100") is None
    # 13-field PAF line whose query name has no "@"
    assert de.from_selfpaf(paf("chr2", 5000, 0, 2000, "+", "chr1", 9000, 5000, 7000, 1900, 2000, "2000M")) is None


def test_normalize_cigar_rejects_invalid_ops():
    # normalize_cigar must raise ValueError for any op other than M, I, D, =, X
    with pytest.raises(ValueError):
        de.normalize_cigar("10M5S")


def test_valley_bimodal_and_monotone():
    hist = {1: 1000, 2: 400, 3: 120, 4: 30, 5: 12, 6: 20, 7: 45, 8: 60, 9: 40, 10: 10}
    assert de.valley(hist) == 5
    assert de.valley({1: 1000, 2: 500, 3: 250, 4: 100, 5: 50}) is None
    assert de.valley({1: 5}) is None


def test_run_atomic_failure(tmp_path):
    """Test that _run with a failing command leaves no output file and logs the error."""
    out_file = str(tmp_path / "output.txt")
    log_file = str(tmp_path / "run.log")
    with pytest.raises(RuntimeError) as exc_info:
        de._run([sys.executable, "-c", "import sys; sys.exit(3)"], out=out_file, log=log_file)
    assert "exit 3" in str(exc_info.value)
    assert not os.path.exists(out_file), f"Output file should not exist after failure, but found {out_file}"
    assert not os.path.exists(f"{out_file}.tmp"), f"Tmp file should not exist, but found {out_file}.tmp"
    assert os.path.exists(log_file), f"Log file should exist after failure"
    with open(log_file) as f:
        log_content = f.read()
    assert "sys.exit(3)" in log_content or "exit 3" in log_content


def test_run_atomic_success(tmp_path):
    """Test that _run with a succeeding command atomically writes to the output file."""
    out_file = str(tmp_path / "output.txt")
    log_file = str(tmp_path / "run.log")
    de._run([sys.executable, "-c", "import sys; sys.stdout.write('hello')"], out=out_file, log=log_file)
    assert os.path.exists(out_file), f"Output file should exist after success"
    assert not os.path.exists(f"{out_file}.tmp"), f"Tmp file should be cleaned up after success"
    with open(out_file) as f:
        assert f.read() == "hello"
    assert os.path.exists(log_file), f"Log file should exist after success"
    with open(log_file) as f:
        log_content = f.read()
    assert "sys.stdout.write" in log_content
