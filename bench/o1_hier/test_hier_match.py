import os
import random
import sys

sys.path.insert(0, os.path.dirname(__file__))
import hier_match as H  # noqa: E402

NODES = {
    ("family", "F1"): frozenset({"a", "b", "c", "d"}),
    ("expansion", "F1.x1"): frozenset({"a", "b"}),
    ("expansion", "F1.x2"): frozenset({"c", "d"}),
    ("family", "F2"): frozenset({"e", "f"}),
    ("expansion", "F2.x1"): frozenset({"e", "f"}),  # identical to its family: the tie policy picks the family
    ("superfamily", "SF1"): frozenset({"a", "b", "c", "d", "e", "f"}),
}
MEMBERS = sorted("abcdef")
FAMILY_OF = {"a": "F1", "b": "F1", "c": "F1", "d": "F1", "e": "F2", "f": "F2"}


def test_truth_one_level_down_is_found_at_the_expansion_level():
    truth = {"T1": frozenset({"a", "b"}), "T2": frozenset({"c", "d"})}
    m = {x.label: x for x in H.Scorer(NODES, MEMBERS).run(truth)}
    assert (m["T1"].level, m["T1"].node_id, m["T1"].jaccard) == ("expansion", "F1.x1", 1.0)
    assert (m["T2"].level, m["T2"].node_id) == ("expansion", "F1.x2")


def test_the_flat_score_misses_what_the_hierarchy_finds():
    truth = {"T1": frozenset({"a", "b"}), "T2": frozenset({"c", "d"})}
    fam = {k: v for k, v in NODES.items() if k[0] == "family"}
    assert H.summary(H.Scorer(fam, MEMBERS).run(truth))["exact"] == 0
    assert H.summary(H.Scorer(NODES, MEMBERS).run(truth))["exact"] == 2


def test_identical_sets_at_two_levels_are_matched_once_with_the_named_tie_policy():
    (m,) = H.Scorer(NODES, MEMBERS).run({"T": frozenset({"e", "f"})})
    assert (m.level, m.node_id) == ("family", "F2")


def test_one_to_one_gives_a_contested_node_to_one_truth_family_only():
    truth = {"T1": frozenset({"e", "f"}), "T2": frozenset({"e", "f"})}
    ms = H.Scorer({("family", "F2"): NODES[("family", "F2")]}, MEMBERS).run(truth)
    assert sorted(m.jaccard for m in ms) == [0.0, 1.0]
    assert H.exact_any(truth, NODES) == {"any": 2, "expansion": 2, "family": 2}


def R(chrom, s, e, strand, exons):
    return H.Record(chrom, s, e, strand, exons)


LOCI = {  # key: (chrom, strands, exons) as load_loci_gff returns them
    "chr1:100-400": ("chr1", frozenset("+"), [(100, 150), (350, 400)]),
    "chr1:1000-1400": ("chr1", frozenset("+"), [(1000, 1100), (1300, 1400)]),
    "chr1:2000-2100": ("chr1", frozenset("-"), [(2000, 2100)]),
    "chr1:3000-3100": ("chr1", frozenset("+"), [(3000, 3100)]),   # folded into chr1:1000-1400
    "chr1:5000-5100": ("chr1", frozenset("+"), [(5000, 5100)]),   # a catalog locus outside every family
}
MEM = ["chr1:100-400", "chr1:1000-1400", "chr1:2000-2100"]
FOLD = {"chr1:3000-3100": "chr1:1000-1400"}


def test_a_shared_exon_on_the_same_strand_joins_and_an_intron_or_the_other_strand_does_not():
    truth = {"T": [R("chr1", 120, 130, "+", [(120, 130)]),      # exon on exon, same strand: joins chr1:100-400
                   R("chr1", 200, 300, "+", [(200, 300)]),      # inside the intron only: no member
                   R("chr1", 2050, 2060, "+", [(2050, 2060)]),  # exon on exon, other strand: no member
                   R("chr1", 3050, 3060, "+", [(3050, 3060)])]}  # a folded record's exon: its representative
    sets, per, cov = H.truth_sets(truth, MEM, fold=FOLD, loci=LOCI)
    assert sets["T"] == frozenset({"chr1:100-400", "chr1:1000-1400"})
    assert [len(x) for x in per["T"]] == [1, 0, 0, 1]
    assert (cov["records_on_members"], cov["records_without_locus"]) == (2, 2)


def test_coverage_separates_an_unclustered_locus_from_no_locus_and_counts_empty_families():
    truth = {"A": [R("chr1", 120, 130, "+", [(120, 130)]), R("chr1", 5010, 5020, "+", [(5010, 5020)])],
             "B": [R("chr1", 9000, 9100, "+", [(9000, 9100)])]}
    sets, per, cov = H.truth_sets(truth, MEM, fold=FOLD, loci=LOCI)
    assert (cov["records_on_unclustered_locus"], cov["records_without_locus"]) == (1, 1)
    assert cov["labels_without_member"] == ["B"] and cov["labels_one_member"] == 1 and cov["labels_no_member"] == 1
    node = {("family", "F"): frozenset(MEM)}
    (m,) = H.Scorer(node, MEM).run({"A": sets["A"]})
    assert H.gene_recall(m, node, per) == 0.5  # the record on the unclustered locus counts against the family


def test_key_mode_needs_identical_coordinates(tmp_path):
    t = tmp_path / "t.tsv"
    t.write_text("chr1\t100\t400\tT1\nchr1\t3000\t3100\tT1\nchr1\t101\t400\tT2\nchr3\t1\t9\tUnassigned_4\n")
    truth = H.load_truth(str(t))
    assert sorted(truth) == ["T1", "T2"] and truth["T2"][0].exons == [(101, 400)]
    sets, _, _ = H.truth_sets(truth, MEM, fold=FOLD, how="key")
    assert sets["T1"] == frozenset({"chr1:100-400", "chr1:1000-1400"}) and sets["T2"] == frozenset()


def test_the_soto_truth_table_is_read_with_its_strand_and_exons(tmp_path):
    t = tmp_path / "t.tsv"
    t.write_text("label\tgene_id\tchrom\tstart\tend\tstrand\texons\nID_1\tG1\tchr1\t101\t500\t-\t101-180,301-500\n")
    assert H.load_truth(str(t)) == {"ID_1": [H.Record("chr1", 101, 500, "-", [(101, 180), (301, 500)])]}


def test_loci_gff_pools_exons_and_uses_the_span_without_exon_lines(tmp_path):
    g = tmp_path / "l.gff3"
    g.write_text("##gff-version 3\n"
                 "chr1\t.\tgene\t10\t90\t.\t+\t.\tID=gene-A;Name=A\n"
                 "chr1\t.\texon\t10\t20\t.\t+\t.\tParent=gene-A;gene=A\n"
                 "chr1\t.\texon\t80\t90\t.\t+\t.\tParent=gene-A;gene=A\n"
                 "chr2\t.\tgene\t5\t50\t.\t-\t.\tID=gene-B;Name=B\n")
    assert H.load_loci_gff(str(g)) == {"chr1:10-90": ("chr1", frozenset("+"), [(10, 20), (80, 90)]),
                                        "chr2:5-50": ("chr2", frozenset("-"), [(5, 50)])}


def test_permutations_keep_sizes_and_within_family_keeps_families():
    rng = random.Random(7)
    for kind in ("global", "within_family"):
        p = H.permutation(MEMBERS, FAMILY_OF, kind, rng)
        assert sorted(p) == MEMBERS and sorted(p.values()) == MEMBERS
        moved = H.relabel(NODES, p)
        assert all(len(moved[k]) == len(v) for k, v in NODES.items())
    p = H.permutation(MEMBERS, FAMILY_OF, "within_family", rng)
    assert all(FAMILY_OF[p[m]] == FAMILY_OF[m] for m in MEMBERS)
    moved = H.relabel(NODES, p)
    assert all(moved[k] == v for k, v in NODES.items() if k[0] in ("family", "superfamily"))


def test_self_match_is_perfect_and_the_null_is_not():
    truth = {k[1]: v for k, v in NODES.items() if k[0] == "expansion"}
    sc = H.Scorer(NODES, MEMBERS)
    assert H.summary(sc.run(truth))["exact"] == 3
    ex, mj = H.null(sc, truth, MEMBERS, FAMILY_OF, "global", 200, 7)
    assert len(ex) == 200 and max(ex) <= 3 and sum(ex) / len(ex) < 3
    assert H.p_value(ex, 3) > 1 / 201  # six members: chance hits are common, as they should be


def test_cli_end_to_end(tmp_path):
    mem = tmp_path / "m.tsv"
    mem.write_text("sd_complex\tsuperfamily\tfamily\tmember\texpansions\n"
                   "CX0\tSF0\tMCL0\tchr1:1-10\tMCL0.x1\nCX0\tSF0\tMCL0\tchr1:21-30\tMCL0.x1\n"
                   "CX0\tSF0\tMCL0\tchr1:41-50\tMCL0.x2\nCX0\tSF0\tMCL0\tchr1:61-70\tMCL0.x2\n")
    tr = tmp_path / "t.tsv"
    tr.write_text("chr1\t2\t5\tA\nchr1\t25\t26\tA\nchr1\t45\t46\tB\nchr1\t65\t69\tB\n")
    gff = tmp_path / "l.gff3"
    gff.write_text("".join(f"chr1\t.\tgene\t{s}\t{s + 9}\t.\t+\t.\tID=gene-{s}\n" for s in (1, 21, 41, 61)))
    out = tmp_path / "o.json"
    H.main([str(mem), str(tr), str(out), "--loci-gff", str(gff), "--null", "20"])
    import json
    r = json.load(open(out))
    assert r["hierarchy"]["exact"] == 2 and r["flat"]["exact"] == 0
    assert r["hierarchy"]["exact_by_level"] == {"expansion": 2}
    assert set(r["null"]) == {"global", "within_family"}
