import gzip
import importlib.util
import os
import subprocess

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    s = importlib.util.spec_from_file_location(f"o3_graph_{name}_t", os.path.join(HERE, f"{name}.py"))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


D, R = _load("o2vc_data"), _load("o2vc_rule")
VCF = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
       "chr1\t101\t.\tA\tG\t30\tPASS\t.\tGT:GQ:PS:DP:AF:PQ\t0|1:30:101:12:0.50:20\n"
       "chr1\t201\t.\tC\tT,G\t20\tMultiallelic\t.\tGT:GQ:DP:AF\t1/2:20:10:0.40,0.60\n"
       "chr1\t301\t.\tT\tC\t5\tLowQual\t.\tGT:GQ:DP:AF\t0/1:5:8:0.30\n")


def test_parse_vcf_reads_every_record_0_based_with_gt_filter_and_depth(tmp_path):
    for p in (tmp_path / "x.vcf", tmp_path / "x.vcf.gz"):
        (gzip.open if p.suffix == ".gz" else open)(p, "wt").write(VCF)
        v = D.parse_vcf(str(p))
        assert v[("chr1", 100)] == R.Call("A", ("G",), (0, 1), "PASS", 12)
        assert v[("chr1", 200)] == R.Call("C", ("T", "G"), (1, 2), "Multiallelic", 10)
        assert v[("chr1", 300)] == R.Call("T", ("C",), (0, 1), "LowQual", 8)


def test_longcallr_runs_with_the_masseq_preset_min_mapq_0_and_no_bam(monkeypatch, tmp_path):
    seen = {}

    def fake(cmd, **kw):
        seen["cmd"] = cmd
        open(str(tmp_path / "o") + ".vcf", "w").write("#CHROM\n")       # longcallR 2.0.0 writes OUT_PREFIX.vcf, plain text (observed)
        return subprocess.CompletedProcess(cmd, 0, "", "")
    monkeypatch.setattr(D.subprocess, "run", fake)
    out = D.run_longcallr("t.bam", "g.fa", str(tmp_path / "o"), binary="LC", threads=2, region="chr1:1-100")
    assert seen["cmd"] == ["LC", "-b", "t.bam", "-f", "g.fa", "-p", "hifi-masseq", "--min-mapq", "0", "--no-bam-output",
                           "-t", "2", "-o", str(tmp_path / "o"), "-r", "chr1:1-100"]
    assert out == str(tmp_path / "o") + ".vcf"


def _fake_run(seen, rc=0, err="", out="", vcf=None):
    """A stand-in for subprocess.run: keeps the command, writes the VCF when given one, answers rc, out and err."""
    def fake(cmd, **kw):
        seen["cmd"] = cmd
        if vcf:
            open(vcf, "w").write("#CHROM\n")
        return subprocess.CompletedProcess(cmd, rc, out, err)
    return fake


def test_a_record_with_a_gt_that_calls_no_allele_is_kept_with_an_empty_gt(tmp_path):
    gts = ["./.", ".|.", ".", "./1"]                    # missing, phased missing, haploid missing, and one allele missing
    p = tmp_path / "m.vcf"
    p.write_text(VCF + "".join(f"chr1\t{401 + i}\t.\tG\tA\t3\tLowQual\t.\tGT:GQ:DP:AF\t{gt}:3:7:0.10\n" for i, gt in enumerate(gts)))
    v = D.parse_vcf(str(p))
    assert [v[("chr1", 400 + i)] for i in range(4)] == [R.Call("G", ("A",), (), "LowQual", 7)] * 4
    assert v[("chr1", 100)].gt == (0, 1)                # the called records are as before


def test_longcallr_gets_no_r_without_a_region(monkeypatch, tmp_path):
    seen = {}
    monkeypatch.setattr(D.subprocess, "run", _fake_run(seen, vcf=str(tmp_path / "o") + ".vcf"))
    assert D.run_longcallr("t.bam", "g.fa", str(tmp_path / "o"), binary="LC") == str(tmp_path / "o") + ".vcf"
    assert seen["cmd"] == ["LC", "-b", "t.bam", "-f", "g.fa", "-p", "hifi-masseq", "--min-mapq", "0", "--no-bam-output",
                           "-t", "2", "-o", str(tmp_path / "o")]


def test_a_failing_longcallr_raises_with_the_exit_code_and_the_last_twenty_stderr_lines(monkeypatch, tmp_path):
    err = "".join(f"line {i}\n" for i in range(30))
    monkeypatch.setattr(D.subprocess, "run", _fake_run({}, rc=101, err=err))
    with pytest.raises(RuntimeError) as e:
        D.run_longcallr("t.bam", "g.fa", str(tmp_path / "o"), binary="LC", region="chr1:0-100")
    msg = str(e.value)
    assert "101" in msg and "line 10" in msg and "line 29" in msg and "line 9" not in msg


def test_longcallr_exiting_0_without_its_vcf_is_an_error_not_a_missing_path(monkeypatch, tmp_path):
    monkeypatch.setattr(D.subprocess, "run", _fake_run({}))
    with pytest.raises(RuntimeError, match="wrote no VCF"):
        D.run_longcallr("t.bam", "g.fa", str(tmp_path / "o"), binary="LC")


def test_the_caller_s_version_is_the_line_its_version_flag_prints(monkeypatch):
    seen = {}
    monkeypatch.setattr(D.subprocess, "run", _fake_run(seen, out="longcallR 2.0.0\n"))
    assert D.longcallr_version("LC") == "longcallR 2.0.0"
    assert seen["cmd"] == ["LC", "--version"]


def test_longcallr_gets_the_directory_of_its_output_prefix_made_for_it_panics_without_one(monkeypatch, tmp_path):
    prefix = str(tmp_path / "vcf" / "F.0")

    def fake(cmd, **kw):
        assert os.path.isdir(tmp_path / "vcf")                                  # longcallR 2.0.0 does not make it (exit 101)
        open(prefix + ".vcf", "w").write("#CHROM\n")
        return subprocess.CompletedProcess(cmd, 0, "", "")
    monkeypatch.setattr(D.subprocess, "run", fake)
    assert D.run_longcallr("t.bam", "g.fa", prefix, binary="LC") == prefix + ".vcf"
