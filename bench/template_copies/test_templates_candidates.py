import os
import sys

import pysam
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import candidates as C  # noqa: E402
import templates as T  # noqa: E402

GTF = """chr1\tt\ttranscript\t101\t400\t.\t+\t.\tgene_id "g1"; transcript_id "g1.a"; reads "5";
chr1\tt\texon\t301\t400\t.\t+\t.\tgene_id "g1"; transcript_id "g1.a";
chr1\tt\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "g1.a";
chr1\tt\ttranscript\t101\t450\t.\t+\t.\tgene_id "g1"; transcript_id "g1.b"; reads "5";
chr1\tt\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "g1.b";
chr1\tt\texon\t301\t450\t.\t+\t.\tgene_id "g1"; transcript_id "g1.b";
chr1\tt\ttranscript\t1001\t1100\t.\t-\t.\tgene_id "g2"; transcript_id "g2.a"; reads "2";
chr1\tt\texon\t1001\t1100\t.\t-\t.\tgene_id "g2"; transcript_id "g2.a";
"""


# ------------------------------------------------------------------------------------------------------- templates

def test_read_gtf_and_representatives(tmp_path):
    p = tmp_path / "base.gtf"
    p.write_text(GTF)
    tx = T.read_gtf(str(p))
    assert tx["g1.a"].exons == [(100, 200), (300, 400)] and tx["g1.a"].reads == 5 and tx["g2.a"].strand == "-"
    reps = T.representatives(tx)
    assert reps["g1"].tid == "g1.b"   # tie on reads (5): the longer exon sum (250 > 200) wins
    assert reps["g2"].tid == "g2.a"
    assert T.loci(tx) == [{"contig": "chr1", "locus": "g1", "iv": [(100, 200), (300, 450)]},
                          {"contig": "chr1", "locus": "g2", "iv": [(1000, 1100)]}]


def test_representative_tie_on_length_takes_smaller_id():
    t1 = T.Transcript("b", "g", "chr1", "+", [(0, 10)], 3)
    t2 = T.Transcript("a", "g", "chr1", "+", [(20, 30)], 3)
    assert T.representatives({"b": t1, "a": t2})["g"].tid == "a"


def test_spliced_and_round_trip(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + "A" * 10 + "CCGG" + "T" * 6 + "GGTA" + "A" * 6 + "\n")
    pysam.faidx(str(fa))
    f = pysam.FastaFile(str(fa))
    plus = T.Transcript("p", "g", "chr1", "+", [(10, 14), (20, 24)], 1)
    minus = T.Transcript("m", "g", "chr1", "-", [(10, 14), (20, 24)], 1)
    assert T.spliced(f, plus) == "CCGGGGTA" and T.spliced(f, minus) == "TACCCCGG"
    T.write_templates({"g": plus}, f, str(tmp_path / "t.fa"), str(tmp_path / "t.tsv"))
    assert (tmp_path / "t.fa").read_text() == ">p\nCCGGGGTA\n"
    assert T.read_templates(str(tmp_path / "t.tsv")) == {"p": plus}


# ------------------------------------------------------------------------------------------------------ candidates

def paf_line(q, qlen, qs, qe, strand, ts, te, matches, block, cg):
    return "\t".join(map(str, [q, qlen, qs, qe, strand, "chr1", 100000, ts, te, matches, block, 60, f"cg:Z:{cg}"])) + "\n"


def test_hits_filters_and_projects(tmp_path):
    p = tmp_path / "h.paf"
    p.write_text(paf_line("ok", 200, 0, 200, "+", 1000, 1700, 180, 200, "100M500N100M")
                 + paf_line("low_id", 200, 0, 200, "+", 3000, 3200, 150, 200, "200M")
                 + paf_line("short", 200, 0, 90, "+", 5000, 5090, 90, 90, "90M")
                 + "\t".join(map(str, ["nocg", 200, 0, 200, "+", "chr1", 100000, 7000, 7200, 200, 200, 60])) + "\n")
    hs = C.hits(str(p))
    assert [h.template for h in hs] == ["ok"] and hs[0].exons == [(1000, 1100), (1600, 1700)]
    assert hs[0].identity == 0.9 and hs[0].aligned == 1.0


def hit(t, s, e, ident=0.9, strand="+", contig="chr1"):
    return C.Hit(t, contig, strand, [(s, e)], ident, 1.0)


def test_uncovered_uses_the_half_rule():
    base = [{"contig": "chr1", "locus": "L", "iv": [(0, 600)]}]
    h1, h2 = hit("a", 0, 1000), hit("b", 0, 2000)   # base covers 600/1000 of h1, 600/2000 of h2
    assert C.uncovered([h1, h2], base) == [h2]


def test_merge_candidates():
    a, b = hit("a", 0, 1000, 0.90), hit("b", 400, 1400, 0.95)   # overlap 600 >= 500: one candidate, b kept
    c = hit("c", 900, 1900, 0.99)                               # overlaps b by 500 of 1000: joins (transitive)
    d = hit("d", 0, 1000, 0.99, strand="-")                     # other strand: separate
    e = hit("e", 1800, 3000, 0.80)                               # overlaps c by 100: separate
    groups = C.merge_candidates([a, b, c, d, e])
    kept = [(k.template, sorted(o.template for o in others)) for k, others in groups]
    assert kept == [("c", ["a", "b"]), ("e", []), ("d", [])]


def test_gtf_without_reads_attribute_takes_the_longest(tmp_path):
    p = tmp_path / "st.gtf"
    p.write_text('chr1\tStringTie\ttranscript\t1\t100\t.\t+\t.\tgene_id "G"; transcript_id "G.1";\n'
                 'chr1\tStringTie\texon\t1\t100\t.\t+\t.\tgene_id "G"; transcript_id "G.1";\n'
                 'chr1\tStringTie\texon\t1\t300\t.\t+\t.\tgene_id "G"; transcript_id "G.2";\n')
    reps = T.representatives(T.read_gtf(str(p)))
    assert reps["G"].tid == "G.2" and reps["G"].reads == 0


def test_a_template_on_a_contig_missing_from_the_genome_fails_loudly(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + "A" * 100 + "\n")
    pysam.faidx(str(fa))
    with pysam.FastaFile(str(fa)) as f, pytest.raises(KeyError, match="chrZ"):
        T.spliced(f, T.Transcript("x", "g", "chrZ", "+", [(0, 10)], 1))
