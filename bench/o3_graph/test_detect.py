import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    """By file path under a unique name: bench has a `{name}.py` in many directories (score.py in 13), and a plain import would
    shadow another directory's module in a combined pytest run."""
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


D = _load("detect")

REF = "ACGTACGTAC"
SITES = {2: "T", 5: "G"}
CONS = {2: "alt", 5: "alt", 7: "ref"}


def test_a_copy_explains_with_every_alt_and_no_reference_conflict():
    assert D.explains(CONS, SITES, REF, {2: "T", 5: "G", 7: "T"})
    assert not D.explains(CONS, SITES, REF, {2: "T", 7: "T"})            # no base at an alt: a representative lacking the exon
    assert not D.explains(CONS, SITES, REF, {2: "T", 5: "G", 7: "A"})    # conflict where the cluster is reference


def test_a_copy_map_is_the_union_of_its_isoforms():
    assert D.union_map([{1: "A", 2: "C"}, {2: "C", 9: "G"}]) == {1: "A", 2: "C", 9: "G"}


def test_novel_bases_are_large_insertions_internal_gaps_and_unaligned_ends():
    # two aligned pieces of a 1,000-bp read with a 150-bp gap between them, one carrying a 40-bp insertion
    recs = [dict(qlen=1000, qs=0, qe=400, cigar="200=40I160="), dict(qlen=1000, qs=550, qe=1000, cigar="450=")]
    assert D.novel_bp(recs, min_len=30) == 40 + 150
    assert D.novel_bp([dict(qlen=1000, qs=0, qe=1000, cigar="500=20I480=")], min_len=30) == 0
    # the aligner stops at the novel exon: the rest of the read is an unaligned end; a short end (< min_len) is not novel
    assert D.novel_bp([dict(qlen=1000, qs=0, qe=700, cigar="700=")], min_len=30) == 300
    assert D.novel_bp([dict(qlen=1000, qs=10, qe=1000, cigar="990=")], min_len=30) == 0


def test_genome_novel_bases_come_from_large_insertions_and_soft_clips():
    assert D.cigar_novel([(0, 300), (1, 45), (0, 300), (4, 60)], min_len=30) == 105
    assert D.cigar_novel([(0, 300), (1, 10), (3, 5000), (0, 300)], min_len=30) == 0


def test_genome_check_keeps_only_what_the_genome_does_not_explain():
    assert D.genome_keeps("psv", median_de=0.02, median_novel=0)
    assert not D.genome_keeps("psv", median_de=0.003, median_novel=0)
    assert D.genome_keeps("branch", median_de=0.001, median_novel=150)
    assert not D.genome_keeps("branch", median_de=0.001, median_novel=0)


def test_novel_intervals_locate_the_insertion_the_gap_and_the_ends():
    recs = [dict(qlen=1000, qs=0, qe=400, cigar="200=40I160="), dict(qlen=1000, qs=550, qe=1000, cigar="450=")]
    assert D.novel_intervals(recs, min_len=30) == [(200, 240), (400, 550)]
    assert D.novel_intervals([dict(qlen=1000, qs=60, qe=700, cigar="300=5X335=")], min_len=30) == [(0, 60), (700, 1000)]


def test_reads_sharing_a_novel_segment_are_one_branch_cluster():
    hits = [("r1#0", "r2#0", 0.995, 140, 150, 150), ("r2#0", "r3#1", 0.99, 145, 150, 160),
            ("r4#0", "r1#0", 0.80, 140, 150, 150), ("r5#0", "r6#0", 0.999, 40, 150, 150)]
    groups = D.link_segments(["r1", "r2", "r3", "r4", "r5", "r6"], hits, min_ident=0.97)
    assert sorted(map(sorted, groups)) == [["r1", "r2", "r3"], ["r4"], ["r5"], ["r6"]]


def test_a_forced_alignment_of_novel_sequence_is_a_dense_error_window():
    # 200 matches, 60 bases aligned at ~50% errors (what an aligner does with a random insert), 200 matches
    dense = "1X1=" * 30
    rec = dict(qlen=460, qs=0, qe=460, cigar="200=" + dense + "200=")
    iv = D.novel_intervals([rec], min_len=30)
    assert len(iv) == 1 and 190 <= iv[0][0] <= 205 and 255 <= iv[0][1] <= 270
    # paralog-level divergence (a mismatch every 30 bases) is not novel sequence
    sparse = "29=1X" * 10
    assert D.novel_intervals([dict(qlen=500, qs=0, qe=500, cigar="100=" + sparse + "100=")], min_len=30) == []


def test_genome_novel_bases_include_a_forced_alignment_split_into_small_insertions():
    # minimap2's alignment of a 150-bp novel exon: several 20-29-bp insertions between tiny matches
    forced = [(7, 213), (1, 23), (7, 2), (1, 22), (7, 5), (1, 20), (8, 1), (7, 4), (1, 21), (8, 1), (7, 1), (8, 1), (7, 3),
              (8, 1), (7, 2), (1, 2), (7, 5), (8, 2), (3, 5178), (7, 4), (1, 29), (7, 5), (2, 1), (7, 300)]
    assert D.cigar_novel(forced, min_len=30) >= 100


def test_a_supplementary_alignment_explains_the_soft_clip_it_covers():
    prim = dict(cigar=[(4, 500), (7, 1000)], reverse=False)            # 500-bp leading clip
    supp = dict(cigar=[(7, 500), (5, 1000)], reverse=False)            # the clipped 500 bp aligned elsewhere (readthrough)
    assert D.genome_novel(prim, [], 1500) == 500
    assert D.genome_novel(prim, [supp], 1500) == 0
    # a reverse-strand supplementary: its stored interval [1000, 1500) is [0, 500) of the read as sequenced
    supp_rc = dict(cigar=[(5, 1000), (7, 500)], reverse=True)
    assert D.genome_novel(prim, [supp_rc], 1500) == 0


def test_a_member_without_a_genome_alignment_counts_as_wholly_novel():
    gst = {"r1": (0.002, 0), "r2": (0.003, 0)}
    de, nov = D.member_genome(["r1", "r2", "r3"], gst, {"r1": 1000, "r2": 1000, "r3": 900})
    assert de == [0.002, 0.003, 1.0] and nov == [0, 0, 900]


def test_the_genome_is_read_lazily_by_slice(tmp_path):
    import pysam
    fa = tmp_path / "g.fa"
    fa.write_text(">c1\nacgtACGTAA\n>c2\nGGGG\n")
    pysam.faidx(str(fa))
    g = D.Genome(str(fa))
    assert len(g["c1"]) == 10 and g["c1"][2:6] == "GTAC" and g["c2"][0:100] == "GGGG"


def test_minus_strand_loci_are_oriented_to_their_transcript():
    seq, ex = D.orient("AACCGGTTAC", [(2, 4), (6, 9)], "-", lo=0)
    assert seq == "GTAACCGGTT" and ex == [(1, 4), (6, 8)]
    assert D.orient("AACCGGTTAC", [(2, 4)], "+", lo=0) == ("AACCGGTTAC", [(2, 4)])


def test_s4_clusters_only_reads_no_catalogued_copy_explains():
    assert D.restrict(["r1", "r2", "r3", "r4"], {"r2", "r4", "r9"}) == ["r2", "r4"]
    assert D.restrict(["r1", "r2"], None) == []                      # no path-O2 table: the arm has no reads


def test_each_call_lands_in_its_arms():
    assert D.arms_of("R0", keep=True) == ["R0", "S3", "S23"] and D.arms_of("R0", keep=False) == ["R0"]
    assert D.arms_of("S1", keep=True) == ["S1", "S12", "S123"]
    assert D.arms_of("S2", keep=True) == ["S12", "S123", "S23", "S4n", "S4"]
    assert D.arms_of("S4p", keep=False) == ["S4n"]
