"""tools/annotated_families.py: copies.tsv is 0-based half-open, a BED element 0-based, a GFF3 gene 1-based closed."""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
HEAD = "family_id\tcopy_idx\ttid\tchrom\tstart\tend\n"


def run(copies, elements, *flags):
    with tempfile.TemporaryDirectory() as d:
        c, e = os.path.join(d, "c.tsv"), os.path.join(d, "e")
        open(c, "w").write(HEAD + copies)
        open(e, "w").write(elements)
        out = subprocess.run([sys.executable, "-I", os.path.join(HERE, "annotated_families.py"), c, e, *flags],
                             capture_output=True, text=True, check=True)
        return out.stdout.split()


def test_a_copy_that_only_touches_an_element_is_not_inside_it():
    # copy bases 1001-2000; BED element bases 901-1000: adjacent, no shared base
    assert run("F\t0\tt\tc1\t1000\t2000\n", "c1\t900\t1000\n", "--majority-in") == []
    assert run("F\t0\tt\tc1\t1000\t2000\n", "c1\t900\t1001\n", "--majority-in") == ["F"]


def test_a_copy_right_after_a_gene_is_unannotated():
    gff = "c1\tx\tgene\t901\t1000\t.\t+\t.\tID=g\n"
    assert run("F\t0\tt\tc1\t1000\t2000\n", gff) == ["F"]
    assert run("F\t0\tt\tc1\t999\t2000\n", gff) == []
