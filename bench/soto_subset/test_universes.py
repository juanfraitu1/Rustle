import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import universes as U  # noqa: E402

F_O = "/mnt/linuxdisk/tmp/o1_layers/v2/dupmasker"
F_O_GFF = "/mnt/linuxdisk/tmp/ourhom/loci.gff3"
SHIPPED = "/mnt/linuxdisk/tmp/o1_close/hier_shipped"
RUNS = "/mnt/linuxdisk/tmp/rustle_figures/runs"
TRUTH_V2 = "/mnt/linuxdisk/tmp/o1_close/lift/soto_truth.v2.tsv"
S1C = os.path.join(HERE, "..", "soto", "soto_famCN_S1C.tsv")
REAL = all(os.path.exists(p) for p in (F_O, F_O_GFF, SHIPPED, RUNS, TRUTH_V2, S1C))


def write(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)


def gff(path, genes):
    """genes: [(name, chrom, start, end, strand, [(s, e)] exons or None)], GFF 1-based closed."""
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for name, c, s, e, st, ex in genes:
            fh.write(f"{c}\t.\tgene\t{s}\t{e}\t.\t{st}\t.\tID=gene-{name};Name={name}\n")
            for i, (a, b) in enumerate(ex or []):
                fh.write(f"{c}\t.\texon\t{a}\t{b}\t.\t{st}\t.\tID=exon-{name}-{i};Parent=gene-{name}\n")


def members(path, fams):
    write(path, ["sd_complex", "superfamily", "family", "member", "expansions"],
          [["-", "-", fam, k, "-"] for fam, keys in fams.items() for k in keys])


def fold(path, pairs):
    write(path, ["annotation", "representative"], pairs)


def test_soto_views():
    labels = {"g1": {"ID_1"}, "g2": {"ID_1", "ID_2"}, "g3": {"Unassigned_9"}, "g4": {"ID_2"}}
    sets, soto_of, multi = U.soto_views(labels, ["ID_1", "ID_2", "ID_3", "Unassigned_9"])
    assert sets == {"ID_1": frozenset({"g1", "g2"}), "ID_2": frozenset({"g2", "g4"}), "ID_3": frozenset()}
    assert soto_of == {"g1": "ID_1", "g3": "Unassigned_9", "g4": "ID_2"} and multi == frozenset({"g2"})


def guided_fixture(tmp_path, extra_gene=False):
    s1c = tmp_path / "s1c.tsv"
    rows = [["g1", "ID_1"], ["g2", "ID_1"], ["g2", "ID_2"], ["g3", "ID_2"], ["g4", "Unassigned_7"]]
    write(s1c, ["Gene ID", "Family ID"], rows + ([["g9", "ID_3"]] if extra_gene else []))
    gff(tmp_path / "loci.gff3", [("g1", "chr1", 100, 200, "+", None), ("g2", "chr1", 300, 400, "+", None),
                                 ("g3", "chr1", 150, 180, "+", None), ("g4", "chr1", 900, 1000, "+", None)])
    members(tmp_path / "members.tsv", {"F1": ["chr1:100-200", "chr1:300-400"]})
    fold(tmp_path / "loci.tsv", [["chr1:150-180", "chr1:100-200"]])
    return [str(tmp_path / x) for x in ("s1c.tsv", "loci.gff3", "members.tsv", "loci.tsv")]


def test_guided_blocks(tmp_path):
    u = U.guided(*guided_fixture(tmp_path))
    assert u.block_of == {"g1": "F:F1", "g2": "F:F1", "g3": "F:F1", "g4": "L:chr1:900-1000"}
    assert u.clustered == frozenset({"g1", "g2", "g3"}) and u.multi == frozenset({"g2"})
    assert u.soto_sets == {"ID_1": frozenset({"g1", "g2"}), "ID_2": frozenset({"g2", "g3"})}
    assert u.soto_of == {"g1": "ID_1", "g3": "ID_2", "g4": "Unassigned_7"}
    assert u.coverage == dict(genes=4, genes_on_family_members=3, loci=3, loci_with_several_genes=1, families=1)
    assert u.genes_of["g3"] == frozenset({"g3"})


def test_guided_gene_without_gff_line_fails_loudly(tmp_path):
    with pytest.raises(ValueError, match="g9"):
        U.guided(*guided_fixture(tmp_path, extra_gene=True))


def test_denovo_reach(tmp_path):
    gff(tmp_path / "loci.gff3", [("L1", "chr1", 100, 500, "+", [(100, 150), (400, 500)]),
                                 ("L2", "chr1", 1000, 1200, "+", [(1000, 1200)]),
                                 ("L3", "chr1", 2000, 2100, "-", [(2000, 2100)]),
                                 ("L4", "chr1", 3000, 3100, "+", [(3000, 3100)])])
    members(tmp_path / "members.tsv", {"F1": ["chr1:100-500", "chr1:1000-1200"]})
    fold(tmp_path / "loci.tsv", [["chr1:3000-3100", "chr1:100-500"]])
    truth = tmp_path / "truth.tsv"
    write(truth, ["label", "gene_id", "chrom", "start", "end", "strand", "exons"],
          [["ID_1", "gA", "chr1", 120, 140, "+", "120-140"],      # exon inside L1's first exon
           ["ID_1", "gB", "chr1", 1100, 1150, "+", "1100-1150"],  # L2
           ["ID_2", "gC", "chr1", 2050, 2060, "+", "2050-2060"],  # L3 is on the other strand: no reach
           ["ID_2", "gD", "chr1", 3050, 3060, "+", "3050-3060"],  # L4, folded into L1
           ["ID_3", "gE", "chr1", 2010, 2020, "-", "2010-2020"],  # L3, a locus in no family
           ["ID_3", "gG", "chr1", 200, 300, "+", "200-300"],      # inside L1's intron: no reach
           ["ID_4", "gF", "chr1", 9000, 9100, "+", "9000-9100"]])  # nothing
    u = U.denovo(str(truth), str(tmp_path / "loci.gff3"), str(tmp_path / "members.tsv"), str(tmp_path / "loci.tsv"))
    assert u.block_of == {"chr1:100-500": "F:F1", "chr1:1000-1200": "F:F1", "chr1:2000-2100": "L:chr1:2000-2100"}
    assert u.soto_sets == {"ID_1": frozenset({"chr1:100-500", "chr1:1000-1200"}), "ID_2": frozenset({"chr1:100-500"}),
                           "ID_3": frozenset({"chr1:2000-2100"}), "ID_4": frozenset()}
    assert u.soto_of == {"chr1:1000-1200": "ID_1", "chr1:2000-2100": "ID_3"} and u.multi == frozenset({"chr1:100-500"})
    assert u.genes_of["chr1:100-500"] == frozenset({"gA", "gD"})
    assert u.coverage == dict(genes=7, records=7, records_on_members=3, records_on_unclustered_only=1,
                              records_without_locus=3, reached_loci=3, reached_clustered=2)
    cn = U.member_cn({"gA": 2.0, "gD": 4.0, "gB": 10.0, "gZ": 99.0}, u.genes_of)
    assert cn == {"chr1:100-500": 3.0, "chr1:1000-1200": 10.0}


@pytest.mark.skipif(not REAL, reason="the frozen inputs are not on this machine")
def test_real_universes_reproduce_plan1_mapping():
    g = U.guided(S1C, F_O_GFF, f"{F_O}/ours.hierarchy.members.tsv", f"{F_O}/ours.loci.tsv")
    assert g.coverage["genes"] == 2334 and g.coverage["genes_on_family_members"] == 1796 and g.coverage["loci"] == 2241
    expect = {"human_A119b": (1211, 155, 1092, 1396, 1115), "human_testis": (286, 100, 2072, 332, 238)}
    for s, (on, unc, none, reached, clustered) in expect.items():
        d = f"{SHIPPED}/{s}"
        u = U.denovo(TRUTH_V2, f"{RUNS}/{s}/{s}.fam.loci.gff3", f"{d}/ours.hierarchy.members.tsv", f"{d}/ours.loci.tsv")
        c = u.coverage
        assert (c["records_on_members"], c["records_on_unclustered_only"], c["records_without_locus"],
                c["reached_loci"], c["reached_clustered"]) == (on, unc, none, reached, clustered), s
        assert len(u.soto_sets) == 491, s
