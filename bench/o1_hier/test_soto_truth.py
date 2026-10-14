import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import soto_truth as S  # noqa: E402
import v1_to_v2 as L  # noqa: E402


def bed(chrom, s, e, strand, sizes, offs, gid):
    f = [chrom, str(s), str(e), "x", "0", strand, str(s), str(e), "0", str(len(sizes)),
         ",".join(map(str, sizes)) + ",", ",".join(map(str, offs)) + ","] + ["."] * 6 + [gid]
    return "\t".join(f) + "\n"


def test_exons_are_the_union_of_the_transcripts():
    lines = [bed("chr1", 100, 400, "+", [50, 100], [0, 200], "G1"), bed("chr1", 120, 500, "+", [60, 80], [0, 300], "G1")]
    assert S.cat_genes(lines) == {"G1": ("chr1", "+", [(100, 180), (300, 400), (420, 500)])}


def test_rows_join_by_gene_id_and_write_gff_coordinates():
    genes = {"G1": ("chr1", "-", [(100, 180), (300, 500)])}
    out, c = S.build([("ID_1", "G1"), ("ID_1", "G9")], genes)
    assert out == [("ID_1", "G1", "chr1", 101, 500, "-", "101-180,301-500")]
    assert c["genes_not_in_cat"] == 1 and c["labels"] == 1


def test_lift_drops_an_exon_in_a_gap_and_shifts_the_rest():
    chain = L.read_chain("chain 1 chr1 1000 + 0 1000 chr1 995 + 0 995 1\n500 5 0\n495\n".splitlines())
    genes = {"G1": ("chr1", "+", [(498, 503), (600, 700)])}
    out, c = S.build([("ID_1", "G1")], genes, L.Lifter(chain))
    assert out == [("ID_1", "G1", "chr1", 596, 695, "+", "596-695")] and c["exons_unliftable"] == 1


def test_s1c_drops_unassigned(tmp_path):
    p = tmp_path / "s1c.tsv"
    p.write_text("Gene Name\tGene ID\tFamily ID\nA\tG1\tID_1\nB\tG2\tUnassigned_3\n")
    assert S.s1c_rows(str(p)) == [("ID_1", "G1")]
