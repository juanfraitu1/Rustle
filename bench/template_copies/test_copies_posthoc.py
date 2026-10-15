import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import copies_posthoc as P  # noqa: E402
from templates import Transcript  # noqa: E402


def test_template_place_four_classes():
    gene = {"source_id": "G", "cls": "in_place", "contig": "c", "start0": "1000", "end": "9000",
            "iv": [(1000, 2000), (8000, 9000)]}
    genes = P.K.LF.Index([gene])
    t = lambda iv: Transcript("t", "t", "c", "+", iv, 3)  # noqa: E731
    assert P.template_place(t([(1000, 2000)]), genes) == "gene_exons"
    assert P.template_place(t([(1900, 2100), (3000, 4000)]), genes) == "exon_overlap"
    assert P.template_place(t([(3000, 4000)]), genes) == "intronic"
    assert P.template_place(t([(20000, 21000)]), genes) == "intergenic"
