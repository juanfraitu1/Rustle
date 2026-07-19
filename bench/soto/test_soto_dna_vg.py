#!/usr/bin/env python3
"""Tests for the DNA variation-graph ceiling demo (bench/soto/soto_dna_vg.py).
Run: PYTHONHASHSEED=0 pytest bench/soto/test_soto_dna_vg.py -v"""
import os
import sys
sys.path.insert(0, "bench/soto")
import soto_dna_vg as V


def test_parse_family_members():
    bed = [
        "chr1\t100\t200\tSRGAP2C|ID_462\t0\t.",
        "chr7\t68507281\t68522154\tPMS2P4|ID_8\t0\t.",
        "chr1\t300\t400\tSRGAP2|ID_462\t0\t.",
    ]
    m = V.parse_family_members(bed, "ID_462")
    assert m == [("SRGAP2C", "chr1", 100, 200), ("SRGAP2", "chr1", 300, 400)]


def test_member_seq_plus_one_offset():
    fa = {"chr1:101-200": "ACGT", "chr7:1-9": "TTTT"}
    # BED start 100 (0-based) -> header start 101 (1-based)
    assert V.member_seq(fa, "chr1", 100, 200) == "ACGT"


def test_read_fasta(tmp_path):
    p = tmp_path / "m.fa"
    p.write_text(">chr1:101-200\nAC\nGT\n>chr7:1-9\nTTTT\n")
    assert V.read_fasta(str(p)) == {"chr1:101-200": "ACGT", "chr7:1-9": "TTTT"}


def test_msa_to_gfa_snp_bubble():
    # col0 "A" invariant; col1 variant (C/G/C); col2-3 "GT" invariant
    rows = ["ACGT", "AGGT", "ACGT"]
    gfa, paths = V.msa_to_gfa(rows, ["m1", "m2", "m3"])
    assert paths["m1"] == ["1", "2", "4"]   # A , C , GT
    assert paths["m2"] == ["1", "3", "4"]   # A , G , GT
    assert paths["m3"] == ["1", "2", "4"]
    assert "S\t2\tC" in gfa and "S\t3\tG" in gfa
    assert "S\t1\tA" in gfa and "S\t4\tGT" in gfa
    assert "P\tm2\t1+,3+,4+\t*" in gfa
    # link m2 traverses 1->3 and 3->4
    links = {(l.split("\t")[1], l.split("\t")[3]) for l in gfa.splitlines() if l.startswith("L\t")}
    assert ("1", "3") in links and ("3", "4") in links


def test_msa_to_gfa_indel_skips_gap_member():
    # col1 gap in m2 -> m2 skips that allele node
    rows = ["ACGT", "A-GT", "ACGT"]
    gfa, paths = V.msa_to_gfa(rows, ["m1", "m2", "m3"])
    assert paths["m1"] == ["1", "2", "3"]
    assert paths["m2"] == ["1", "3"]        # skips the gap region
    assert paths["m3"] == ["1", "2", "3"]
    links = {(l.split("\t")[1], l.split("\t")[3]) for l in gfa.splitlines() if l.startswith("L\t")}
    assert ("1", "3") in links              # m2's skip link


def test_msa_to_gfa_all_identical_single_node():
    rows = ["ACGT", "ACGT"]
    gfa, paths = V.msa_to_gfa(rows, ["a", "b"])
    assert paths["a"] == ["1"] and paths["b"] == ["1"]
    assert "S\t1\tACGT" in gfa
    assert not any(l.startswith("L\t") for l in gfa.splitlines())  # no links, one node


def test_load_detection():
    tsv = [
        "family_id\tgene\tchrom\tstart\tend\tdetected\trecovered_by",
        "ID_462\tSRGAP2C\tchr1\t121194173\t121402237\tY\tRNA-split",
        "ID_26\tSLC9B1P1\tchr16\t32804124\t32821138\tN\t",
    ]
    d = V.load_detection(tsv)
    assert d[("chr1", 121194173, 121402237)] == (True, "RNA-split")
    assert d[("chr16", 32804124, 32821138)] == (False, "")


def test_node_colour():
    det = {"a": True, "b": True, "c": False}
    assert V.node_colour({"a", "b"}, det) == V.GREEN     # all recovered
    assert V.node_colour({"c"}, det) == V.RED            # all DNA-only
    assert V.node_colour({"a", "c"}, det) == V.GREY      # mixed = shared/conserved
    assert V.node_colour(set(), det) == V.GREY           # empty


def test_colours_csv():
    paths = {"a": ["1", "2", "4"], "b": ["1", "3", "4"], "c": ["1", "2", "4"]}
    det = {"a": True, "b": False, "c": True}
    csv = V.colours_csv(paths, det)
    rows = dict(l.split(",") for l in csv.strip().splitlines()[1:])
    assert rows["1"] == V.GREY    # a,b,c -> mixed
    assert rows["2"] == V.GREEN   # a,c -> recovered
    assert rows["3"] == V.RED     # b -> DNA-only
    assert rows["4"] == V.GREY    # a,b,c -> mixed
    assert csv.splitlines()[0] == "Node,Colour"


def test_legend_tsv():
    members = [("AMY1A", "chr1", 100, 200), ("AC1", "chr1", 300, 400)]
    det = {"AMY1A": True, "AC1": False}
    rec = {"AMY1A": "RNA-split", "AC1": ""}
    cause = {"AC1": "expressed-K=0: exon-homogenized"}
    leg = V.legend_tsv(members, det, rec, cause)
    assert leg.splitlines()[0] == "gene\tlocus\tdetected\trecovered_by\tcause\tcolour"
    assert "AMY1A\tchr1:100-200\tY\tRNA-split\t\t#1e8e3e" in leg          # green: empty cause
    assert "AC1\tchr1:300-400\tN\t\texpressed-K=0: exon-homogenized\t#d93025" in leg  # red: cause shown


def test_load_causes():
    tsv = ["family_id\tgene\tchrom\tstart\tend\ta\tb\tc\td\te\tf\tcause",
           "ID_131\tAC105272.1\tchr1\t103516936\t103517135\t.\t.\t.\t.\t.\t.\texpressed-K=0: exon-homogenized"]
    d = V.load_causes(tsv)
    assert d[("chr1", 103516936, 103517135)] == "expressed-K=0: exon-homogenized"


def test_build_family_presence_and_colours(monkeypatch):
    # isolate build_family's orchestration from abpoa (keeps the test hermetic — no pyabpoa needed).
    # these 3 seqs are equal-length with one SNP, so the identity "MSA" is the correct alignment.
    monkeypatch.setattr(V, "abpoa_msa", lambda seqs: list(seqs))
    fa = {"chr1:101-200": "ACGTACGT", "chr1:301-400": "ACGTACGT", "chr1:501-600": "ACGTTCGT"}
    members = [("g1", "chr1", 100, 200), ("g2", "chr1", 300, 400), ("g3", "chr1", 500, 600)]
    detection = {("chr1", 100, 200): (True, "RNA-split"),
                 ("chr1", 300, 400): (True, "projection"),
                 ("chr1", 500, 600): (False, "")}
    r = V.build_family("ID_test", members, fa, detection, {})
    assert r["n_members"] == 3 and r["n_present"] == 3 and r["missing"] == []
    # every gene is a P-line
    plines = {l.split("\t")[1] for l in r["gfa"].splitlines() if l.startswith("P\t")}
    assert plines == {"g1", "g2", "g3"}
    # g3 differs (T at col4) -> its unique node is red; shared nodes grey; g1/g2 identical
    assert V.RED in r["colours"] and V.GREY in r["colours"]
    assert "g3\tchr1:500-600\tN" in r["legend"]


def test_build_family_missing_member_logged_not_counted():
    fa = {"chr1:101-200": "ACGT"}   # g2 absent from fa
    members = [("g1", "chr1", 100, 200), ("g2", "chr9", 100, 200)]
    detection = {("chr1", 100, 200): (True, "RNA-split")}
    r = V.build_family("ID_x", members, fa, detection, {})
    assert r["n_members"] == 2 and r["n_present"] == 1 and r["missing"] == ["g2"]


def test_apply_window_extracts_and_clamps():
    seqs = ["AAAACCCCGGGG", "TTTTAAAACCCC", "GGGG"]
    spans = [(0, 6), None, (0, 4)]        # m0 -> [0:6] clamped to 4; m1 -> None fallback first 4; m2 -> [0:4]
    out = V.apply_window(seqs, spans, window_bp=4)
    assert out == ["AAAA", "TTTT", "GGGG"]


def test_window_members_passthrough_when_all_fit():
    seqs = ["ACGT", "ACGA", "ACGT"]        # all <= window_bp -> unchanged, no minimap2 invoked
    out, windowed, span = V.window_members(seqs, window_bp=100)
    assert out == seqs and windowed is False and span is None


def test_build_family_windowed_flag_false_for_small(monkeypatch):
    monkeypatch.setattr(V, "abpoa_msa", lambda s: list(s))
    fa = {"chr1:101-108": "ACGTACGT", "chr1:201-208": "ACGTTCGT"}
    members = [("g1", "chr1", 100, 108), ("g2", "chr1", 200, 208)]
    detection = {("chr1", 100, 108): (True, "RNA-split"), ("chr1", 200, 208): (False, "")}
    r = V.build_family("ID_s", members, fa, detection, {})
    assert r["windowed"] is False and r["n_present"] == 2 and r["window_span"] is None
