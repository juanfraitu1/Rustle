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
    assert rows[0] == ["read", "copy_idx", "is_primary", "AS", "chain", "ref_start", "ref_end"]
    body = sorted(rows[1:])
    assert body == [["r1", "1", "1", "100", "1200-1500", "1100", "1600"], ["r1", "2", "0", "99", "5200-5500", "5100", "5600"],
                    ["r2", "1", "1", "100", "", "1100", "1300"]]
    stats = [l.rstrip("\n").split("\t") for l in open(out / "pool.stats.tsv")]
    assert stats[0] == ["family", "molecules", "records", "secondary_kept"]
    assert stats[1] == ["F1", "2", "3", "1"]


# ---- lift ------------------------------------------------------------------------------------------------------
# F1: copy 1 = c1:1000-2000 (representative, 10 reads), copy 2 = c1:5000-6000 (3 reads).
_COPIES = {("F1", "1"): ("c1", 1000, 2000, 10), ("F1", "2"): ("c1", 5000, 6000, 3)}
_REPS = {"F1": "1"}
# copy 2 -> rep: first half on '+' (identity), second half on '-' (reversed).
_LIFTS = {("F1", "2"): ([0, 500], [(0, 500, 0, "+"), (500, 1000, 500, "-")])}


def test_private_junction_admits_only_at_own_copy():
    lifts = {("F1", "2"): ([0], [(0, 500, 0, "+")])}   # only the first 500 bp of copy 2 align
    rows = [("r1", "2", 1, ((5200, 5300),)), ("r2", "2", 1, ((5700, 5800),))]
    out = famasm.lift_rows(_COPIES, _REPS, lifts, "F1", rows)
    assert out[0] == ("r1", "2", 1, "199-300", 0)
    assert out[1][3].startswith("c2:") and out[1][4] == 1


def test_missing_lift_is_private(tmp_path):
    out = famasm.lift_rows(_COPIES, _REPS, {}, "F1", [("r1", "2", 1, ((5200, 5300),))])
    assert out == [("r1", "2", 1, "c2:5200-5300", 1)]


def test_inverse_lift_roundtrip():
    chain = ((5200, 5300), (5700, 5800))
    fam = famasm.fw_layer.fam_chain(_COPIES, _REPS, _LIFTS, "F1", "2", chain)
    assert fam == ((199, 300), (699, 800))
    assert famasm.inverse_lift(_LIFTS, _COPIES, _REPS, "F1", "2", fam) == chain
    # at the representative the inverse is a shift
    assert famasm.inverse_lift(_LIFTS, _COPIES, _REPS, "F1", "1", ((199, 300),)) == ((1200, 1300),)
    # a family junction outside every aligned column has no placement
    assert famasm.inverse_lift({("F1", "2"): ([0], [(0, 500, 0, "+")])}, _COPIES, _REPS, "F1", "2", ((699, 800),)) is None


# ---- admit -----------------------------------------------------------------------------------------------------
J1, J2, J3 = (1200, 1300), (1400, 1500), (1600, 1700)


def test_family_chain_floor_counts_molecules_over_copies():
    rows = [("r1", "1", 1, "199-300", 0), ("r2", "2", 1, "199-300", 0), ("r3", "1", 1, "399-500", 0),
            ("r4", "2", 1, "c2:5200-5300", 1), ("r5", "2", 0, "c2:5200-5300", 1)]
    assert famasm.family_chains(rows) == {((199, 300),)}


def test_admit_whole_needs_two_certified():
    chain = (J1, J2)
    two = [("r1", chain, 1100, 1800), ("r2", chain, 1150, 1750)]
    assert famasm.admit_at_copy(chain, two) == "whole"
    assert famasm.admit_at_copy(chain, two[:1]) is None


def test_admit_compat_covers_every_junction_twice():
    chain = (J1, J2, J3)
    left = [("a", (J1, J2), 1100, 1550), ("b", (J1, J2), 1150, 1580)]
    right = [("c", (J2, J3), 1350, 1800), ("d", (J2, J3), 1320, 1750)]
    assert famasm.admit_at_copy(chain, left + right) == "compat"
    assert famasm.admit_at_copy(chain, left + [("e", (J3,), 1550, 1800)]) is None
    # end-incompatible: the read's first exon runs into chain's first intron
    bad = [("c", (J2, J3), 1250, 1800), ("d", (J2, J3), 1320, 1750)]
    assert famasm.admit_at_copy(chain, left + bad) is None
    # a read whose junctions are not a consecutive run of chain's never counts
    assert famasm.admit_at_copy(chain, left + [("f", (J1, J3), 1100, 1800), ("g", (J1, J3), 1100, 1800)]) is None


def _admit_family(pool_rows, lift_rows, certified, lifts=_LIFTS):
    return famasm.admit_family("F1", pool_rows, lift_rows, certified, _COPIES, _REPS, lifts)


def test_chain_comes_from_certified_record():
    # r1 and r2 certified to copy 2; records on copies 1 and 2 (same family chain 199-300).
    pool = [("r1", "1", 1, ((1200, 1300),), 1100, 1400), ("r1", "2", 0, ((5200, 5300),), 5100, 5400),
            ("r2", "1", 1, ((1200, 1300),), 1100, 1400), ("r2", "2", 0, ((5200, 5300),), 5100, 5400),
            ("r3", "1", 1, ((1200, 1300),), 1100, 1400)]
    lift = [("r1", "1", 1, "199-300", 0), ("r1", "2", 0, "199-300", 0), ("r2", "1", 1, "199-300", 0),
            ("r2", "2", 0, "199-300", 0), ("r3", "1", 1, "199-300", 0)]
    rows, stats = _admit_family(pool, lift, {"r1": "2", "r2": "2", "r3": "1"})
    assert rows == [("2", "199-300", "5200-5300", "whole", 2, 5100, 5400)]
    assert stats["certified_without_record"] == 0 and stats["family_chains"] == 1
    # r3 certified to copy 2 but with a record on copy 1 only: counted, never used
    rows, stats = _admit_family(pool, lift, {"r1": "2", "r2": "2", "r3": "2"})
    assert rows[0][4] == 2 and stats["certified_without_record"] == 1


def test_inverse_lift_dedups():
    # synthetic lift where representative positions 199-300 and 1199-1300 both come from copy-2 positions 199-300
    lifts = {("F1", "2"): ([0, 0], [(0, 500, 0, "+"), (0, 500, 1000, "+")])}
    pool = [("r1", "2", 1, ((5200, 5300),), 5100, 5400), ("r2", "2", 1, ((5200, 5300),), 5150, 5450),
            ("r3", "1", 1, ((1200, 1300),), 1100, 1400), ("r4", "1", 1, ((2200, 2300),), 2100, 2400),
            ("r5", "1", 1, ((2200, 2300),), 2150, 2450)]
    lift = [("r1", "2", 1, "199-300", 0), ("r2", "2", 1, "199-300", 0), ("r3", "1", 1, "199-300", 0),
            ("r4", "1", 1, "1199-1300", 0), ("r5", "1", 1, "1199-1300", 0)]
    rows, stats = _admit_family(pool, lift, {"r1": "2", "r2": "2", "r3": "1", "r4": "1", "r5": "1"}, lifts)
    at2 = [r for r in rows if r[0] == "2"]
    assert len(at2) == 1 and at2[0][2] == "5200-5300" and at2[0][4] == 2 and at2[0][1] == "1199-1300;199-300"
    assert stats["uninvertible"] == 0


# ---- gtf -------------------------------------------------------------------------------------------------------
_COPIES_TSV = ("family_id\tcopy_idx\tchrom\tstrand\tlocus_start\tlocus_end\tn_reads\n"
               "F1\t1\tc1\t+\t1000\t2000\t10\nF1\t2\tc1\t+\t5000\t6000\t3\n")
_BASE = ('c1\trustle\ttranscript\t1101\t1400\t.\t+\t.\tgene_id "DN_1"; transcript_id "DN_1"; reads "3";\n'
         'c1\trustle\texon\t1101\t1200\t.\t+\t.\tgene_id "DN_1"; transcript_id "DN_1"; exon_number "1";\n'
         'c1\trustle\texon\t1301\t1400\t.\t+\t.\tgene_id "DN_1"; transcript_id "DN_1"; exon_number "2";\n')
_ADMIT_HDR = "copy_idx\tfam_chain\tgenomic_chain\tvia\tcertified_reads\tspan_start\tspan_end\n"


def _gtf_setup(tmp_path, admit_rows):
    (tmp_path / "copies.tsv").write_text(_COPIES_TSV)
    (tmp_path / "base.gtf").write_text(_BASE)
    (tmp_path / "out" / "admit").mkdir(parents=True)
    (tmp_path / "out" / "admit" / "F1.tsv").write_text(_ADMIT_HDR + admit_rows)
    famasm.gtf(str(tmp_path / "out"), str(tmp_path / "base.gtf"), str(tmp_path / "copies.tsv"), str(tmp_path / "famasm.gtf"))
    return (tmp_path / "famasm.gtf").read_text(), (tmp_path / "out" / "gtf.stats.tsv").read_text()


def test_no_duplicate_of_base(tmp_path):
    out, stats = _gtf_setup(tmp_path, "1\t199-300\t1200-1300\twhole\t2\t1100\t1400\n")
    assert out == _BASE
    assert stats == "appended\talready_in_base\n0\t1\n"


def test_base_lines_verbatim(tmp_path):
    out, stats = _gtf_setup(tmp_path, "2\t199-300\t5200-5300\twhole\t2\t5100\t5400\n")
    assert out.startswith(_BASE)
    tail = out[len(_BASE):].splitlines()
    assert [l.split("\t")[:8] for l in tail] == [["c1", "famasm", "transcript", "5101", "5400", ".", "+", "."],
                                                  ["c1", "famasm", "exon", "5101", "5200", ".", "+", "."],
                                                  ["c1", "famasm", "exon", "5301", "5400", ".", "+", "."]]
    assert 'gene_id "FAMASM_F1_2"; transcript_id "FAMASM_F1_2_1"; reads "2"; famasm "1"; family "F1"; via "whole"; certified_reads "2";' in tail[0]
    assert stats == "appended\talready_in_base\n1\t0\n"


def test_single_copy_untouched(tmp_path):
    (tmp_path / "copies.tsv").write_text("family_id\tcopy_idx\tchrom\tstrand\tlocus_start\tlocus_end\tn_reads\nF1\t1\tc1\t+\t1000\t2000\t10\n")
    (tmp_path / "base.gtf").write_bytes(_BASE.encode() + b"c1\trustle\ttranscript\t3000\t3100\t.\t-\t.\tgene_id \"DN_2\"; transcript_id \"DN_2\";\n")
    (tmp_path / "out").mkdir()
    famasm.gtf(str(tmp_path / "out"), str(tmp_path / "base.gtf"), str(tmp_path / "copies.tsv"), str(tmp_path / "famasm.gtf"))
    assert (tmp_path / "famasm.gtf").read_bytes() == (tmp_path / "base.gtf").read_bytes()


# ---- score -----------------------------------------------------------------------------------------------------
def _score_setup(tmp_path, famasm_multi="12.5\t2000\t300\t15.0", single="12.9\t1000\t129\t12.9", base_kept=True,
                 ms_run="110\t\t\t2", pairs_run="0.730\t0.745"):
    r = tmp_path
    (r / "copies.tsv").write_text(_COPIES_TSV)
    (r / "base.gtf").write_text(_BASE)
    extra = ('c1\tfamasm\ttranscript\t5101\t5400\t.\t+\t.\tgene_id "FAMASM_F1_2"; transcript_id "FAMASM_F1_2_1"; famasm "1";\n'
             'c1\tfamasm\texon\t5101\t5200\t.\t+\t.\tgene_id "FAMASM_F1_2"; transcript_id "FAMASM_F1_2_1"; famasm "1";\n')
    (r / "famasm.gtf").write_text((_BASE if base_kept else _BASE.replace("1101", "1102", 1)) + extra)
    (r / "famasm.env").write_text(f"SAMPLE=t\nSPECIES=human\nPREFIX=x\nBASE_GTF={r / 'base.gtf'}\nCOPIES={r / 'copies.tsv'}\n")
    (r / "strata.txt").write_text("reference multi-exon transcripts: {}\narm\tstratum\tref_tx\trecalled\trecall%\tqueries\t=\tprecision%\n"
                                  "base\tsingle\t1000\t129\t12.9\t1000\t129\t12.9\nbase\tmulti\t1000\t110\t11.0\t2000\t280\t14.0\n"
                                  f"famasm\tsingle\t1000\t129\t{single}\nfamasm\tmulti\t1000\t125\t{famasm_multi}\n")
    (r / "mainstays.txt").write_text("run\tfamily\texpressed\trecalled\tfamilies\trel_foreign\tunrel_foreign\tunrel_names\tannotated\trecalled_annot\n"
                                     "base.tsv\tALL\t123\t108\t\t\t2\t\t131\t108\n" + f"famasm.fam.copies.tsv\tALL\t123\t{ms_run}\t\t131\t108\n")
    (r / "pairs.txt").write_text("# level\nrun\tpred_pairs\tignored_noncoding\thit\tfalse\tof_which_related\tprecision\trecall\n"
                                 "base.tsv\t2228\t49502\t1617\t611\t524\t0.726\t0.741\n" + f"famasm.fam.copies.tsv\t2000\t1\t1\t1\t1\t{pairs_run}\n")
    (r / "gtf.stats.tsv").write_text("appended\talready_in_base\n1\t0\n")
    return famasm.score(str(r))


def test_score_pass(tmp_path):
    v = _score_setup(tmp_path)
    assert "M1 PASS" in v and "M2 PASS" in v and "M3 PASS" in v and "M4 PASS" in v and "verdict PASS" in v
    assert "single_copy_identical yes" in v and "base_lines_kept yes" in v


def test_score_fails_on_small_gain_and_lost_base_line(tmp_path):
    v = _score_setup(tmp_path, famasm_multi="11.9\t2000\t300\t15.0")
    assert "M1 FAIL" in v and "verdict FAIL" in v
    v = _score_setup(tmp_path, base_kept=False)
    assert "M3 FAIL" in v and "base_lines_kept no" in v
    v = _score_setup(tmp_path, ms_run="107\t\t\t2")
    assert "M4 FAIL" in v
    v = _score_setup(tmp_path, pairs_run="0.710\t0.745")
    assert "M4 FAIL" in v
