import json
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_reassign as RR  # noqa: E402

CL_HEAD = "cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n"
CP_HEAD = "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\tgene_id\n"


def toy(tmp_path, copies_family_y="F"):
    run = tmp_path / "runs" / "human_testis"
    run.mkdir(parents=True)
    s = "human_testis"
    (run / f"{s}.fam.clusters.tsv").write_text(
        CL_HEAD + "F\t2\t1\t0\tNA\tc\t1\t100\nF\t2\t1\t0\tNA\tc\t201\t300\nG\t2\t1\t0\tNA\tc\t401\t500\n"
        "G\t2\t1\t0\tNA\tc\t601\t700\n")
    (run / f"{s}.fam.loci.gff3").write_text(
        "##gff-version 3\n" + "".join(f"c\t.\tgene\t{a}\t{b}\t.\t+\t.\tID=gene-{n};Name={n}\n"
                                      for n, a, b in (("x", 1, 100), ("y", 201, 300), ("g1", 401, 500),
                                                      ("g2", 601, 700), ("x2", 50, 90))))
    (run / f"{s}.fam.loci.tsv").write_text("annotation\trepresentative\nc:50-90\tc:1-100\n")
    rows = [("F", "x2"), ("F" if copies_family_y == "F" else copies_family_y, "y"), ("G", "g1"), ("G", "g2")]
    (run / f"{s}.fam.copies.tsv").write_text(CP_HEAD + "".join(
        f"{f}\t{i}\tt{i}\tc\t0\t1\t1\t+\t3\t0-1\t{g}\n" for i, (f, g) in enumerate(rows)))
    out = tmp_path / "out" / s
    out.mkdir(parents=True)
    (out / "graph.tsv").write_text("c:1-100\tc:201-300\t0.960000\nc:1-100\tc:401-500\t0.990000\n"
                                   "c:401-500\tc:601-700\t0.990000\n")
    return tmp_path


def test_copy_nodes_map_genes_through_the_fold(tmp_path):
    t = toy(tmp_path)
    p = RR.paths(str(t / "runs"), str(t / "out"), "human_testis")
    shipped = RR.R.read_clusters(p["shipped_clusters"])
    rows, nodes = RR.copy_nodes(p["copies"], p["gff"], p["fold"], shipped)
    assert nodes == ["c:1-100", "c:201-300", "c:401-500", "c:601-700"]
    assert RR.gene_nodes(p["gff"])["x2"] == "c:50-90"          # folded into c:1-100 by the fold table


def test_copy_nodes_refuse_a_row_off_its_shipped_family(tmp_path):
    t = toy(tmp_path, copies_family_y="G")
    p = RR.paths(str(t / "runs"), str(t / "out"), "human_testis")
    with pytest.raises(SystemExit):
        RR.copy_nodes(p["copies"], p["gff"], p["fold"], RR.R.read_clusters(p["shipped_clusters"]))


def test_apply_relabels_the_copy_table(tmp_path):
    t = toy(tmp_path)
    assert RR.main(["apply", "human_testis", "--out", str(t / "out"), "--runs", str(t / "runs")]) == 0
    p = RR.paths(str(t / "runs"), str(t / "out"), "human_testis")
    fams = [l.split("\t")[0] for l in open(p["new_copies"]).read().splitlines()[1:]]
    assert fams == ["G", "G", "G", "G"]                      # x follows its 0.99 edge; y, left alone, follows x
    res = json.load(open(p["apply"]))
    assert (res["moves"], res["families_before"], res["families_after"], res["lost"]) == (2, 2, 1, [])
    assert (res["largest_before"], res["largest_after"], res["copies_relabelled"]) == (2, 4, 2)


def test_soto_pairs():
    fam = {"a": "F", "b": "F", "c": "G", "d": "G"}
    labels = {"a": {"S1"}, "b": {"S1"}, "c": {"S1"}, "d": {"S2"}}
    r = RR.soto_pairs(fam, labels)
    # family pairs: (a,b) S1=S1 yes, (c,d) no -> precision 1/2; Soto pairs: (a,b), (a,c), (b,c) -> recall 1/3
    assert (r["tp"], r["pred"], r["truth"]) == (1, 2, 3)
    assert r["precision"] == 0.5 and abs(r["recall"] - 1 / 3) < 1e-12
