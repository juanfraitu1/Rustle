import json
import os
import random
import sys

import pysam
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_tg as R  # noqa: E402

E1, I1, E2 = 300, 500, 300          # exon, intron, exon lengths of gene A
A0, C0 = 1000, 20000                # gene A at 1000; its exact copy at 20000


def toy_genome(path):
    rng = random.Random(3)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    gene = rand(E1) + "GT" + rand(I1 - 4) + "AG" + rand(E2)
    seq = rand(A0) + gene + rand(C0 - A0 - len(gene)) + gene + rand(30000 - C0 - len(gene))
    path.write_text(">chr1\n" + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n")
    pysam.faidx(str(path))
    return seq


def toy_inputs(tmp_path):
    seq = toy_genome(tmp_path / "g.fa")
    a = [(A0, A0 + E1), (A0 + E1 + I1, A0 + E1 + I1 + E2)]
    (tmp_path / "base.gtf").write_text(
        f'chr1\tt\ttranscript\t{a[0][0] + 1}\t{a[1][1]}\t.\t+\t.\tgene_id "A"; transcript_id "A.1"; reads "5";\n'
        + "".join(f'chr1\tt\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "A"; transcript_id "A.1";\n' for s, e in a))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(seq)}]}
    reads = [("c1", 0, C0, [(0, E1), (3, I1), (0, E2)]), ("c2", 0, C0 + 50, [(0, E1 - 50), (3, I1), (0, 200)]),
             ("c3", 0, C0 + 100, [(0, E1 - 100), (3, I1), (0, E2)]), ("a1", 0, A0, [(0, E1), (3, I1), (0, E2)])]
    with pysam.AlignmentFile(str(tmp_path / "r.bam"), "wb", header=header) as out:
        for name, flag, start, cig in sorted(reads, key=lambda r: r[2]):
            x = pysam.AlignedSegment(out.header)
            x.query_name, x.flag, x.reference_id, x.reference_start, x.mapping_quality = name, flag, 0, start, 60
            x.cigartuples = cig
            x.query_sequence = "A" * sum(n for op, n in cig if op == 0)
            out.write(x)
    pysam.index(str(tmp_path / "r.bam"))
    return tmp_path


def test_end_to_end_on_a_toy_copy(tmp_path):
    d = toy_inputs(tmp_path)
    out = str(d / "tg")
    assert R.main(["templates", "--base", str(d / "base.gtf"), "--genome", str(d / "g.fa"), "--out", out]) == 0
    assert R.main(["project", "--out", out, "--index", str(d / "g.fa"), "--threads", "1"]) == 0
    assert R.main(["assemble", "--base", str(d / "base.gtf"), "--bam", str(d / "r.bam"), "--out", out,
                   "--label", "toy"]) == 0
    rep = json.load(open(os.path.join(out, "report.json")))
    assert rep["merged_candidates"] == 1 and rep["loci"] == 1 and rep["transcripts"] == 1
    exons = [ln.split("\t")[3:5] for ln in open(os.path.join(out, "tg.gtf")) if ln.split("\t")[2] == "exon"]
    assert exons == [[str(C0 + 1), str(C0 + E1)], [str(C0 + E1 + I1 + 1), str(C0 + E1 + I1 + E2)]]
    assert 'template "A.1"' in open(os.path.join(out, "tg.gtf")).read()
    merged = open(os.path.join(out, "merged.gtf")).read()
    assert merged.startswith(open(d / "base.gtf").read()) and "TG_toy_000000" in merged


def test_assemble_refuses_unprojected_batches(tmp_path):
    d = toy_inputs(tmp_path)
    out = str(d / "tg")
    R.main(["templates", "--base", str(d / "base.gtf"), "--genome", str(d / "g.fa"), "--out", out])
    with pytest.raises(SystemExit, match="missing"):
        R.main(["assemble", "--base", str(d / "base.gtf"), "--bam", str(d / "r.bam"), "--out", out, "--label", "x"])


def run_until_projected(d):
    out = str(d / "tg")
    R.main(["templates", "--base", str(d / "base.gtf"), "--genome", str(d / "g.fa"), "--out", out])
    R.main(["project", "--out", out, "--index", str(d / "g.fa"), "--threads", "1"])
    return out


def test_contig_names_that_disagree_fail_loudly(tmp_path):
    d = toy_inputs(tmp_path)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrZ", "LN": 1000}]}
    with pysam.AlignmentFile(str(d / "z.bam"), "wb", header=header):
        pass
    pysam.index(str(d / "z.bam"))
    out = run_until_projected(d)
    with pytest.raises(SystemExit, match="name contigs differently"):
        R.main(["assemble", "--base", str(d / "base.gtf"), "--bam", str(d / "z.bam"), "--out", out, "--label", "t"])


def test_assemble_is_deterministic(tmp_path):
    d = toy_inputs(tmp_path)
    out = run_until_projected(d)
    runs = []
    for _ in range(2):
        R.main(["assemble", "--base", str(d / "base.gtf"), "--bam", str(d / "r.bam"), "--out", out, "--label", "t"])
        runs.append((open(os.path.join(out, "tg.gtf")).read(), open(os.path.join(out, "report.json")).read()))
    assert runs[0] == runs[1] and runs[0][0]
