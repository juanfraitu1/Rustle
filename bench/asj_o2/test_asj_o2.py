"""Tests for asj_o2.py (docs/PREREG_asj_o2_2026-10-10.md)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import asj_o2 as A  # noqa: E402


def test_cs_walk_plus_strand():
    # query = primary span starting at 1000; aligned from query offset 2: 3 matches, SNV (B base g, primary a), 4 matches, deletion of 2 B bases, 6 matches
    snvs, runs, indels = A.cs_walk(":3*ga:4-tt:6", "+", 2, 16, 1000)
    assert snvs == {1005: ("A", "G")}
    assert runs == [(1002, 1005), (1006, 1010), (1010, 1016)]
    assert indels == [1010]


def test_cs_walk_minus_strand_complements_and_reverses():
    # query (primary span at 0) aligned on '-': the alignment walks the query from qe - 1 downwards
    snvs, runs, indels = A.cs_walk(":2*ca:3", "-", 0, 6, 0)
    # walk: q 5,4 match; SNV at q 3 (B base c -> primary-strand G, primary base a -> T); q 2,1,0 match
    assert snvs == {3: ("T", "G")}
    assert sorted(runs) == [(0, 3), (4, 6)]
    assert indels == []


def test_cs_walk_insertion_marks_indel_and_advances_query():
    snvs, runs, indels = A.cs_walk(":2+ac:2", "+", 0, 6, 100)
    assert snvs == {} and runs == [(100, 102), (104, 106)] and indels == [102]


def test_het_sites_drop_snvs_near_indels():
    snvs = {100: ("A", "G"), 120: ("C", "T")}
    assert A.het_sites(snvs, [104], 5) == {120: ("C", "T")}
    assert A.het_sites(snvs, [106], 5) == snvs


def test_read_obs_walks_cigar():
    # 2S 3= 1X 2D 3= 100N 4= ; read starts at ref 10
    cig = [(4, 2), (7, 3), (8, 1), (2, 2), (7, 3), (3, 100), (7, 4)]
    seq = "nn" + "ACG" + "T" + "CCA" + "GGAA"
    obs, introns, span, xs = A.read_obs(cig, 10, seq, [10, 13, 14, 16, 120, 300])
    assert obs == {10: "A", 13: "T", 16: "C", 120: "G"}        # 14 deleted, 300 outside
    assert introns == [(19, 119)] and span == (10, 123) and xs == [(13, "T")]


def test_certificate_one_site_never_assigns():
    assert A.certify([("A", "A", "G")]) == (None, 0.001)


def test_certificate_two_sites_assign():
    hap, p = A.certify([("A", "A", "G"), ("C", "C", "T")])
    assert hap == "S" and p < 1e-3


def test_certificate_tie_abstains_and_neither_counts_as_column():
    assert A.certify([("A", "A", "G"), ("T", "C", "T")])[0] is None
    hap, p = A.certify([("G", "A", "G"), ("T", "C", "T"), ("N", "G", "A")])
    assert hap == "B" and abs(p - A.upper_tail(2, 3, A.EPS)) < 1e-18


def test_upper_tail_binomial():
    assert abs(A.upper_tail(1, 1, 0.001) - 0.001) < 1e-15
    assert abs(A.upper_tail(2, 2, 0.001) - 1e-6) < 1e-15
    assert A.upper_tail(0, 3, 0.001) == 1.0


def test_junction_tests_spanning_used_and_skips():
    j = (100, 200)
    reads = [(0, (j,), (50, 300))] * 4 + [(1, (), (50, 300))] * 4 + [(1, (), (150, 300))] * 5 + [(None, (j,), (50, 300))]
    t = A.junction_tests(reads)
    assert len(t) == 1
    r = t[0]
    assert (r["u0"], r["s0"], r["u1"], r["s1"]) == (4, 4, 0, 4) and r["dpsi"] == -1.0 and r["p"] < 0.05
    # constitutive in both groups -> not tested; carried by < 3 reads -> not tested
    assert A.junction_tests([(0, (j,), (0, 300))] * 3 + [(1, (j,), (0, 300))] * 3) == []
    assert A.junction_tests([(0, (j,), (0, 300))] * 2 + [(1, (), (0, 300))] * 5) == []


def test_bh():
    q = A.bh([0.01, 0.04, 0.03, 0.5])
    assert [round(x, 4) for x in q] == [0.04, 0.0533, 0.0533, 0.5]


def test_pick_anchor_rules():
    cand = {10: (10, "G", 8), 20: (20, "T", 2), 30: (6, "A", 6), 40: (9, "C", 9)}
    # 10: depth 18 frac .44 ok; 20: alt 2 < 5; 30: depth 12 ok; 40: depth 18 ok -> highest depth, smallest pos
    assert A.pick_anchor(cand) == 10
    assert A.pick_anchor({50: (30, "A", 5)}) is None              # minor fraction 5/35 < 0.35


def test_x_candidates_need_five_reads_of_one_base():
    xs = [[(5, "A")]] * 5 + [[(5, "C")]] * 4 + [[(9, "T")]] * 4
    assert A.x_candidates(xs) == {5: ("A", 5)}


def test_classify_anchor():
    snvs = {100: ("A", "G"), 150: ("C", "T")}
    het = {100: ("A", "G")}
    runs, indels = [(0, 100), (101, 150), (151, 400)], [146]
    assert A.classify(100, "A", "G", het, snvs, runs, indels) == "ALLELIC"
    assert A.classify(100, "A", "C", het, snvs, runs, indels) == "UNRESOLVED"     # third base
    assert A.classify(200, "T", "C", het, snvs, runs, indels) == "NON-ALLELIC"
    assert A.classify(149, "T", "C", het, snvs, runs, indels) == "UNRESOLVED"     # near an indel
    assert A.classify(150, "C", "T", het, snvs, runs, indels) == "UNRESOLVED"     # SNV dropped near an indel
    assert A.classify(500, "T", "C", het, snvs, runs, indels) == "UNRESOLVED"     # not covered


def test_concordance():
    pairs = [("ref", "S")] * 9 + [("alt", "B")] * 9 + [("ref", "B")] * 2
    assert abs(A.concordance(pairs) - 18 / 20) < 1e-12
    assert A.concordance([]) is None


def test_parity_halves():
    het = {10: ("A", "G"), 20: ("C", "T"), 30: ("G", "A")}
    assert A.halves(het) == ({10: ("A", "G"), 30: ("G", "A")}, {20: ("C", "T")})


def test_gene_table_counts_exons_per_transcript():
    gff = "\n".join([
        "c\t.\tgene\t101\t500\t.\t+\t.\tID=gene-A;Name=A;gene_biotype=protein_coding",
        "c\t.\tmRNA\t101\t500\t.\t+\t.\tID=rna-A1;Parent=gene-A",
        "c\t.\texon\t101\t150\t.\t+\t.\tID=e1;Parent=rna-A1",
        "c\t.\texon\t300\t500\t.\t+\t.\tID=e2;Parent=rna-A1",
        "c\t.\tgene\t601\t700\t.\t-\t.\tID=gene-B;Name=B;gene_biotype=pseudogene",
        "c\t.\texon\t601\t700\t.\t-\t.\tID=e3;Parent=gene-B",
        "#comment"])
    assert A.gene_table(gff.splitlines()) == [("c", 100, 500, "gene-A", "A", "+", 2), ("c", 600, 700, "gene-B", "B", "-", 1)]


def _gene_reads():
    s_seq = ["N"] * 350
    s_seq[50], s_seq[150] = "A", "C"
    b_seq = ["N"] * 450
    b_seq[50], b_seq[150] = "G", "T"
    s_read = ([(7, 250), (3, 100), (7, 100)], 50, "".join(s_seq), 60)
    b_read = ([(7, 50), (8, 1), (7, 99), (8, 1), (7, 299)], 50, "".join(b_seq), 0)
    return [s_read] * 6 + [b_read] * 6


def test_process_gene_both_arms():
    het = {100: ("A", "G"), 200: ("C", "T")}
    pri = {100: "A", 200: "C"}
    r = A.process_gene(_gene_reads(), het, dict(het), [(0, 100), (101, 200), (201, 600)], [], lambda p: pri.get(p, "N"))
    assert (r["n_reads"], r["n_het"]) == (12, 2)
    assert (r["armA"]["S"], r["armA"]["B"], r["armA"]["abstain"]) == (6, 6, 0)
    assert len(r["testsA"]) == 1 and r["testsA"][0]["intron"] == (300, 400) and r["testsA"][0]["dpsi"] == -1.0
    assert r["T1"] == {"both": 0, "discordant": 0}
    an = r["anchor"]
    assert (an["pos"], an["ref"], an["alt"], an["class"], an["conc"], an["n_conc"]) == (100, "A", "G", "ALLELIC", None, 0)
    assert (an["mapq0_ref"], an["mapq0_alt"]) == (0.0, 1.0)
    assert len(r["testsR"]) == 1 and r["testsR"][0]["u0"] == 6


def test_allelic_alignment_is_the_best_primary_covering_ninety_percent():
    def rec(q, ql, qs, qe, matches, tp):
        return "\t".join([q, str(ql), str(qs), str(qe), "+", "B", "1000", "0", "100", str(matches), "100", "60", f"tp:A:{tp}", "cs:Z::100"])
    lines = [rec("g1", 100, 0, 95, 90, "P"), rec("g1", 100, 0, 99, 99, "S"), rec("g2", 100, 0, 80, 80, "P"), rec("g3", 100, 5, 100, 80, "P"),
             rec("g3", 100, 0, 50, 50, "P")]
    best = A.allelic_alignments(lines)
    assert sorted(best) == ["g1", "g3"] and best["g3"][2] == "5" and best["g1"][9] == "90"
