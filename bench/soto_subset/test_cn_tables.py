import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cn_tables as CN  # noqa: E402

S1E_HDR = ["rid", "kind", "s1e_row", "gene", "gene_id", "chrom", "start", "end", "fam_HG002"]
REAL = all(os.path.exists(p) for p in (CN.S1E, CN.READ_DEPTH, CN.S1C, CN.TRUTH_V2))


def write(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)


def test_num():
    assert CN.num("3.5") == 3.5 and CN.num("0") == 0.0
    assert CN.num("") is None and CN.num("nan") is None and CN.num("N/A") is None
    assert CN.num("abc") is None and CN.num(None) is None


def test_assembly_median_of_resolved_rows(tmp_path):
    p = tmp_path / "s1e.tsv"
    write(p, S1E_HDR, [["S0", "S1E", 0, "A", "g1", "chr1", 100, 200, "18.0"],
                       ["S1", "S1E", 1, "A", "g1", "chr1", 300, 400, "20.0"],
                       ["S2", "S1E", 2, "A", "g1", "chr1", 500, 600, ""],
                       ["S3", "S1E", 3, "B", "g2", "chr2", 100, 200, ""],
                       ["C0", "CTRL", -1, "C", "", "chr3", 1, 20001, "2.0"]])
    vals, unresolved = CN.assembly_cn(str(p))
    assert vals == {"g1": 19.0} and unresolved == {"g2"}


def test_s1c_repeats_and_conflicts(tmp_path):
    p = tmp_path / "s1c.tsv"
    hdr = ["Gene ID", "Family ID", "Median famCN"]
    write(p, hdr, [["g1", "ID_1", "4.5"], ["g1", "ID_2", "4.5"], ["g2", "ID_1", "N/A"]])
    assert CN.s1c_cn(str(p)) == {"g1": 4.5} and CN.s1c_genes(str(p)) == ["g1", "g2"]
    write(p, hdr, [["g1", "ID_1", "4.5"], ["g1", "ID_2", "5.0"]])
    with pytest.raises(ValueError, match="g1"):
        CN.s1c_cn(str(p))


def test_read_depth_column(tmp_path):
    p = tmp_path / "rd.tsv"
    write(p, ["gene_id", "famCN", "famCN_sotoiv"], [["g1", "3.0", "2.5"], ["g2", "4.0", "nan"]])
    assert CN.read_depth_cn(str(p)) == {"g1": 2.5}


def test_account():
    assert CN.account(["g1", "g2", "g3"], {"g1": 2.0}, {"g2"}) == {"valued": 1, "unresolved": 1, "absent": 1}


def test_check_coordinates(tmp_path):
    truth = tmp_path / "t.tsv"
    write(truth, ["label", "gene_id", "chrom", "start", "end", "strand", "exons"],
          [["ID_1", "g1", "chr1", 101, 200, "+", "101-200"], ["ID_1", "g2", "chr1", 301, 400, "+", "301-400"]])
    s1e = tmp_path / "s1e.tsv"
    write(s1e, S1E_HDR, [["S0", "S1E", 0, "A", "g1", "chr1", 90, 101, "2"],    # bases 91-101: shares 101
                         ["S1", "S1E", 1, "B", "g2", "chr1", 400, 410, "2"],   # bases 401-410: past the body (400)
                         ["S2", "S1E", 2, "C", "g9", "chr1", 1, 50, "2"]])     # gene not in the lifted table
    ok, bad, absent = CN.check_coordinates(str(s1e), str(truth))
    assert ok == 1 and absent == 1 and [b[0] for b in bad] == ["S1"]


def test_write_table_is_read_by_the_recipe(tmp_path):
    sys.path.insert(0, os.path.join(HERE, "..", "soto"))
    import soto_replication as R
    p = tmp_path / "asm.tsv"
    CN.write_table(str(p), {"g1": 18.0, "g2": 2.5})
    assert R.load_famcn(str(p), "famCN") == {"g1": 18.0, "g2": 2.5}


@pytest.mark.skipif(not REAL, reason="the frozen inputs are not on this machine")
def test_real_tables(capsys):
    assert CN.main([]) == 0
    genes = CN.s1c_genes(CN.S1C)
    asm, unresolved = CN.assembly_cn(CN.S1E)
    assert len(genes) == 2334
    assert CN.account(genes, asm, unresolved) == {"valued": 1534, "unresolved": 257, "absent": 543}
    assert CN.account(genes, CN.read_depth_cn(CN.READ_DEPTH)) == {"valued": 2334, "unresolved": 0, "absent": 0}
    assert CN.account(genes, CN.s1c_cn(CN.S1C)) == {"valued": 1793, "unresolved": 0, "absent": 541}
    ok, bad, absent = CN.check_coordinates(CN.S1E, CN.TRUTH_V2)
    assert (ok, len(bad), absent) == (1710, 0, 121)
