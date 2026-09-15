import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ape_substrate as asub  # noqa: E402

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
