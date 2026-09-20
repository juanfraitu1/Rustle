import argparse
import os
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ape_substrate as asub  # noqa: E402
import dup_evidence as de  # noqa: E402

GFF = [
    "chr7\tRefSeq\tgene\t100\t900\t.\t+\t.\tID=gene-A;Name=A\n",
    "chr7\tRefSeq\tmRNA\t100\t900\t.\t+\t.\tID=rna-A1;Parent=gene-A\n",
    "chr7\tRefSeq\texon\t100\t200\t.\t+\t.\tID=exon-A1;Parent=rna-A1\n",
    "chr7\tRefSeq\tgene\t1000\t1900\t.\t+\t.\tID=gene-B;Name=B\n",
    "chr7\tRefSeq\texon\t1000\t1100\t.\t+\t.\tID=exon-B1;Parent=gene-B\n",
    "chr7\tRefSeq\tpseudogene\t2000\t2900\t.\t+\t.\tID=gene-C;Name=C\n",
]


def test_gff_batches_keep_children_with_gene():
    b = asub.gff_batches(GFF, max_genes=2)
    assert len(b) == 2
    assert [l.split("\t")[8].split(";")[0] for l in b[0]] == ["ID=gene-A", "ID=rna-A1", "ID=exon-A1", "ID=gene-B", "ID=exon-B1"]
    assert b[1] == [GFF[5]]


def test_select_substrate_handles_translocation():
    counts = {("chr7", "t7"): 900, ("chr7", "t3"): 20, ("chr17", "t17"): 600, ("chr17", "t5"): 300,
              ("chr15", "t15"): 700, ("chr16", "t16"): 800, ("chr16", "t1"): 50}
    assert asub.select_substrate(counts) == ["t15", "t16", "t17", "t5", "t7"]


def _lift_args(tmp_path, max_genes=1500, budget=10):
    return argparse.Namespace(species="X", target=str(tmp_path / "target.fa"), out=str(tmp_path / "out"),
                               max_genes=max_genes, threads=2, budget=budget)


def _seed_batches(tmp_path, max_genes=1500):
    bdir = tmp_path / "out" / "batches"
    bdir.mkdir(parents=True)
    (bdir / "chr7.000.gff").write_text(GFF[0])
    (bdir / "max_genes.txt").write_text(str(max_genes))
    return bdir


def test_cmd_lift_index_budget_expiry_no_liftoff_no_mmi(tmp_path, monkeypatch):
    _seed_batches(tmp_path)
    calls = []

    def fake_run_budget(cmd, budget, stdout=None, log=None):
        calls.append(cmd[0])
        return False  # simulate budget expiry

    monkeypatch.setattr(de, "run_budget", fake_run_budget)
    a = _lift_args(tmp_path)
    asub.cmd_lift(a)
    assert calls == ["minimap2"]  # Liftoff was never invoked
    assert not os.path.exists(a.target + ".mmi")
    assert not os.path.exists(a.target + ".mmi.tmp")


def test_cmd_lift_index_build_calledprocesserror_exits(tmp_path, monkeypatch):
    _seed_batches(tmp_path)

    def fake_run_budget(cmd, budget, stdout=None, log=None):
        open(cmd[2], "w").write("partial")  # cmd[2] is the .mmi.tmp path
        raise subprocess.CalledProcessError(1, cmd)

    monkeypatch.setattr(de, "run_budget", fake_run_budget)
    a = _lift_args(tmp_path)
    with pytest.raises(SystemExit) as exc:
        asub.cmd_lift(a)
    assert exc.value.code != 0
    assert not os.path.exists(a.target + ".mmi.tmp")
    assert not os.path.exists(a.target + ".mmi")


def test_cmd_lift_liftoff_calledprocesserror_removes_tmp_and_exits(tmp_path, monkeypatch):
    _seed_batches(tmp_path)
    a = _lift_args(tmp_path)
    open(a.target + ".mmi", "w").write("fake-index")  # skip the index-build branch

    def fake_run_budget(cmd, budget, stdout=None, log=None):
        out_tmp = cmd[cmd.index("-o") + 1]
        open(out_tmp, "w").write("partial")
        raise subprocess.CalledProcessError(1, cmd)

    monkeypatch.setattr(de, "run_budget", fake_run_budget)
    with pytest.raises(SystemExit) as exc:
        asub.cmd_lift(a)
    assert exc.value.code != 0
    bdir = f"{a.out}/batches"
    assert not os.path.exists(f"{bdir}/chr7.000.lifted.gff3.tmp")
    assert not os.path.exists(f"{bdir}/chr7.000.lifted.gff3")


def test_cmd_lift_rebatching_max_genes_mismatch_exits(tmp_path, monkeypatch):
    _seed_batches(tmp_path, max_genes=1500)

    def boom(cmd, budget, stdout=None, log=None):
        raise AssertionError("run_budget should not be called on a max-genes mismatch")

    monkeypatch.setattr(de, "run_budget", boom)
    a = _lift_args(tmp_path, max_genes=500)
    with pytest.raises(SystemExit) as exc:
        asub.cmd_lift(a)
    assert exc.value.code != 0


def test_cmd_lift_regeneration_deletes_stale_lifted_files(tmp_path, monkeypatch):
    bdir = tmp_path / "out" / "batches"
    bdir.mkdir(parents=True)
    stale = {
        "chr7.000.lifted.gff3": "stale-lifted",
        "chr7.000.unmapped.txt": "stale-unmapped",
        "chr7.000.liftoff.log": "stale-log",
    }
    for name, content in stale.items():
        (bdir / name).write_text(content)
    monkeypatch.setattr(asub, "HUMAN_CHROMS", ())  # skip tabix; no new .gff produced
    a = _lift_args(tmp_path, max_genes=500)
    open(a.target + ".mmi", "w").write("fake-index")  # skip the index-build branch too
    asub.cmd_lift(a)
    for name in stale:
        assert not os.path.exists(bdir / name), f"{name} should have been deleted on regeneration"
    assert (bdir / "max_genes.txt").read_text().strip() == "500"
