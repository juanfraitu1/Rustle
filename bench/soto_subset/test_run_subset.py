import json
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_subset as RS  # noqa: E402

REAL = all(os.path.exists(p) for p in (RS.F_O, RS.F_O_GFF, RS.SHIPPED, RS.RUNS, RS.CN.S1E, RS.CN.READ_DEPTH,
                                       RS.CN.TRUTH_V2))


def write(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)


def toy(tmp_path, s1e_chrom="chr1"):
    """Seven genes: ID_1 = {g1, g2} and ID_2 = {g3, g4} inside our family F1 at copy numbers 2 and 10; ID_3 =
    {g6, g7} with g7 on an unclustered locus; g5 Unassigned. -> the path overrides for main()."""
    loc = {"g1": ("chr1", 100, 200), "g2": ("chr1", 300, 400), "g3": ("chr1", 500, 600), "g4": ("chr1", 700, 800),
           "g5": ("chr1", 900, 1000), "g6": ("chr2", 100, 200), "g7": ("chr2", 300, 400)}
    lab = {"g1": "ID_1", "g2": "ID_1", "g3": "ID_2", "g4": "ID_2", "g5": "Unassigned_1", "g6": "ID_3", "g7": "ID_3"}
    s1c_cn = {"g1": "2.0", "g2": "2.1", "g3": "10.0", "g4": "10.2"}
    write(tmp_path / "s1c.tsv", ["Gene ID", "Family ID", "Median famCN"],
          [[g, lab[g], s1c_cn.get(g, "N/A")] for g in sorted(lab)])
    with open(tmp_path / "loci.gff3", "w") as fh:
        for g, (c, s, e) in sorted(loc.items()):
            fh.write(f"{c}\t.\tgene\t{s}\t{e}\t.\t+\t.\tID=gene-{g};Name={g}\n")
    fo = tmp_path / "fo"
    fo.mkdir()
    write(fo / "ours.hierarchy.members.tsv", ["sd_complex", "superfamily", "family", "member", "expansions"],
          [["-", "-", "F1", f"chr1:{s}-{e}", "-"] for g, (c, s, e) in sorted(loc.items()) if g in ("g1", "g2", "g3", "g4")]
          + [["-", "-", "F2", "chr2:100-200", "-"], ["-", "-", "F2", "chr2:500-600", "-"]])
    write(fo / "ours.loci.tsv", ["annotation", "representative"], [])
    asm = {"g1": 2, "g2": 2, "g3": 10, "g4": 10}
    write(tmp_path / "s1e.tsv", ["rid", "kind", "s1e_row", "gene", "gene_id", "chrom", "start", "end", "fam_HG002"],
          [[f"S{i}", "S1E", i, g, g, s1e_chrom, loc[g][1] - 1, loc[g][2], asm[g]] for i, g in enumerate(sorted(asm))])
    rd = {"g1": 2.1, "g2": 2.3, "g3": 9.8, "g4": 10.4, "g5": 3.0, "g6": 4.0, "g7": 8.0}
    write(tmp_path / "rd.tsv", ["gene_id", "famCN_sotoiv"], [[g, v] for g, v in sorted(rd.items())])
    write(tmp_path / "truth.tsv", ["label", "gene_id", "chrom", "start", "end", "strand", "exons"],
          [[lab[g], g, c, s, e, "+", f"{s}-{e}"] for g, (c, s, e) in sorted(loc.items()) if g != "g5"])
    return ["--s1c", str(tmp_path / "s1c.tsv"), "--s1e", str(tmp_path / "s1e.tsv"), "--read-depth",
            str(tmp_path / "rd.tsv"), "--truth-v2", str(tmp_path / "truth.tsv"), "--fo", str(fo), "--fo-gff",
            str(tmp_path / "loci.gff3")]


def load(path):
    with open(path) as fh:
        return json.load(fh)


def test_inputs_only(tmp_path):
    out = tmp_path / "out"
    assert RS.main(["guided", "--out", str(out), "--inputs-only", *toy(tmp_path)]) == 0
    r = load(out / "inputs.json")
    assert r["members"] == 7 and r["coverage"]["genes_on_family_members"] == 5
    assert r["cn_genes"]["assembly"] == {"valued": 4, "unresolved": 0, "absent": 3}
    assert r["cn_genes"]["read_depth"]["valued"] == 7 and r["cn_genes"]["s1c"]["valued"] == 4
    assert r["common_members"] == 4 and r["families_ge2_in_common"] == 2
    assert r["coordinate_check"] == {"ok": 4, "disagree": 0, "gene_absent": 0}
    assert sorted(os.listdir(out)) == ["inputs.json"]


def test_full_toy_run_is_deterministic(tmp_path):
    args = toy(tmp_path)
    for run in ("a", "b"):
        assert RS.main(["guided", "--out", str(tmp_path / run), "--null", "20", "--boot", "20", "--no-recipe",
                        *args]) == 0
    t1, t2 = load(tmp_path / "a" / "t1.json"), load(tmp_path / "a" / "t2.json")
    assert t1["contained"] == 2 and t1["families_ge2"] == 3
    assert t1["violation_classes"] == {"coverage": 1, "cover": 0, "conflict": 0}
    assert t2["common"]["assembly"]["observed"]["exact"] == 2 and t2["common"]["assembly"]["observed"]["ari"] == 1.0
    assert t2["ours_alone"]["exact"] == 0 and t2["ceiling"]["exact"] == 2
    with open(tmp_path / "a" / "violations.tsv") as fh:
        assert fh.read() == "label\tsize\tblocks\tclass\nID_3\t2\t2\tcoverage\n"
    for f in ("t1.json", "t2.json", "violations.tsv"):
        with open(tmp_path / "a" / f) as x, open(tmp_path / "b" / f) as y:
            assert x.read() == y.read(), f


def test_coordinate_disagreement_stops_the_run(tmp_path):
    with pytest.raises(SystemExit, match="disagrees"):
        RS.main(["guided", "--out", str(tmp_path / "out"), "--inputs-only", *toy(tmp_path, s1e_chrom="chr9")])


def test_denovo_needs_a_sample(tmp_path):
    with pytest.raises(SystemExit):
        RS.main(["denovo", "--out", str(tmp_path / "out")])


@pytest.mark.skipif(not REAL, reason="the frozen inputs are not on this machine")
def test_real_inputs_only(tmp_path):
    assert RS.main(["guided", "--out", str(tmp_path / "g"), "--inputs-only"]) == 0
    r = load(tmp_path / "g" / "inputs.json")
    assert r["members"] == 2334 and r["coverage"]["genes_on_family_members"] == 1796
    assert r["cn_genes"]["assembly"] == {"valued": 1534, "unresolved": 257, "absent": 543}
    assert r["cn_genes"]["s1c"] == {"valued": 1793, "unresolved": 0, "absent": 541}
    assert r["cn_genes"]["read_depth"]["valued"] == 2334
    assert r["common_members"] == 1534 and r["families_ge2_in_common"] == 357
    assert r["coordinate_check"] == {"ok": 1710, "disagree": 0, "gene_absent": 121}
    for s, on in (("human_A119b", 1211), ("human_testis", 286)):
        assert RS.main(["denovo", "--sample", s, "--out", str(tmp_path / s), "--inputs-only"]) == 0
        assert load(tmp_path / s / "inputs.json")["coverage"]["records_on_members"] == on
