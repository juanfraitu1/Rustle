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
    keep = {("r2", "A"), ("r4", "A"), ("r9", "A")}
    assert D.restrict(["r1", "r2", "r3", "r4"], "A", keep) == ["r2", "r4"]
    assert D.restrict(["r1", "r2"], "A", None) == []                 # no path-O2 table: the arm has no reads


def test_s4_keys_unexplained_reads_by_family():
    # r2 is unexplained only in sister family B, where it was hosted by mistake; path O2 explains it in its own family A,
    # so it must not enter A's S4 pool
    keep = {("r2", "B"), ("r4", "A")}
    assert D.restrict(["r1", "r2", "r4"], "A", keep) == ["r4"]
    assert D.restrict(["r1", "r2", "r4"], "B", keep) == ["r2"]


def test_each_call_lands_in_its_arms():
    assert D.arms_of("R0", keep=True) == ["R0", "S3", "S23"] and D.arms_of("R0", keep=False) == ["R0"]
    assert D.arms_of("S1", keep=True) == ["S1", "S12", "S123"]
    assert D.arms_of("S2", keep=True)[:5] == ["S12", "S123", "S23", "S4n", "S4"]
    assert D.arms_of("S4p", keep=False) == ["S4n"]


def _mosaic_case(n_sites=30):
    ref = "A" * 1000
    pos = list(range(10, 10 + 30 * n_sites, 30))
    sites = {p: "C" for p in pos}
    cons = {p: "alt" for p in pos}
    half = set(pos[: n_sites // 2])
    copy0 = {p: ("C" if p in half else "A") for p in range(1000)}       # carries the first half of the alts
    copy1 = {p: ("C" if p in sites and p not in half else "A") for p in range(1000)}
    cov = {"m": {p: ref[p] for p in range(1000)}}
    return ref, sites, cons, {"copy0": copy0, "copy1": copy1}, cov


def test_a_cluster_far_from_every_single_copy_is_a_mosaic_copy():
    ref, sites, cons, maps, cov = _mosaic_case()
    # every alt is carried by some copy, so no alt is unexplained: the old rule leaves it `unresolved`
    assert D.status_of(["m"], cons, sites, ref, maps, cov, {})[0] == "unresolved"
    st, div = D.status_of(["m"], cons, sites, ref, maps, cov, {}, mosaic=True)
    assert st == "mosaic" and abs(div - 15 / 1000) < 1e-9          # the nearest copy differs at 15 of 1,000 bases


def test_a_mosaic_needs_the_nearest_copy_beyond_delta():
    ref, sites, cons, maps, cov = _mosaic_case(n_sites=12)          # nearest copy 6 / 1,000 = 0.6% < delta
    assert D.status_of(["m"], cons, sites, ref, maps, cov, {}, mosaic=True)[0] == "unresolved"


def test_a_copy_covering_too_little_of_the_cluster_is_no_candidate():
    ref, sites, cons, maps, cov = _mosaic_case()
    maps["copy2"] = {p: "A" for p in range(100)}                     # matches the consensus, but over 10% of it
    assert D.nearest_copy(cons, sites, ref, maps, set(range(1000)))[0] == "copy0"


def test_the_union_o2_arm_and_the_mosaic_arms_land_in_their_arms():
    assert D.arms_of("S1m", keep=True) == ["S123m"] and D.arms_of("S1m", keep=False) == []
    assert D.arms_of("Gp", keep=True) == ["Gn", "G"] and D.arms_of("Gp", keep=False) == ["Gn"]
    assert D.arms_of("S2", keep=True)[:8] == ["S12", "S123", "S23", "S4n", "S4", "S123m", "Gn", "G"]


def test_h_pools_the_family_s_unexplained_reads_across_its_hosts():
    by_host = {("F", 0): ["a", "b"], ("F", 1): ["c", "d"], ("G", 0): ["e"]}
    keep = {("a", "F"), ("c", "F"), ("d", "F"), ("e", "G"), ("b", "G")}
    assert D.family_pool(by_host, "F", keep) == ["a", "c", "d"]
    assert D.family_pool(by_host, "F", None) == []
    assert D.arms_of("Hp", keep=True) == ["Hn", "H"] and D.arms_of("Hp", keep=False) == ["Hn"]
    a = D.arms_of("S2", keep=True)
    assert "Hn" in a and a.index("Hn") + 1 == a.index("H")


def test_a_crashing_vg_msga_is_retried_in_fasta_order(tmp_path, monkeypatch):
    calls = []

    def fake_sh(cmd, binary_out=None):
        calls.append(cmd)
        if cmd.startswith("vg msga") and " -a" not in cmd:
            raise SystemExit("failed: vg msga ... Signal 6")
        if binary_out:
            open(binary_out, "wb").close()
    monkeypatch.setattr(D, "sh", fake_sh)
    g = {"c1": "ACGT" * 500}
    cps = [dict(idx=0, contig="c1", start=100, end=900, strand="+", isoforms={"t": [(100, 300), (500, 900)]}),
           dict(idx=1, contig="c1", start=1000, end=1800, strand="-", isoforms={"u": [(1000, 1800)]})]
    D.family_graph(str(tmp_path / "f"), g, cps)
    msga = [c for c in calls if c.startswith("vg msga")]
    assert len(msga) == 2 and msga[1].endswith(" -a")


def test_a_linked_cluster_with_few_private_alts_far_from_every_single_copy_is_a_mosaic():
    # KIR (round 3): 5 copy-private alts (est_div 0.005 < delta, so no PSV call), the other 35 alts split between two copies;
    # the nearest single copy is 2.2% away -- the documented mosaic rule applies to every cluster no copy explains
    ref = "A" * 1000
    pos = list(range(10, 10 + 25 * 40, 25))
    sites, cons = {p: "C" for p in pos}, {p: "alt" for p in pos}
    private, a, b = set(pos[:5]), set(pos[5:23]), set(pos[23:])
    maps = {"copy0": {p: ("C" if p in a else "A") for p in range(1000)},
            "copy1": {p: ("C" if p in b else "A") for p in range(1000)}}
    cov = {m: {p: (sites[p] if p in sites else ref[p]) for p in range(1000)} for m in ("m1", "m2", "m3")}
    vec = {m: {p: "alt" for p in pos} for m in cov}
    assert D.status_of(["m1", "m2", "m3"], cons, sites, ref, maps, cov, vec)[0] == "unexplained"     # est_div 0.005: no call
    st, div = D.status_of(["m1", "m2", "m3"], cons, sites, ref, maps, cov, vec, mosaic=True)
    assert st == "mosaic" and abs(div - 22 / 1000) < 1e-9


def test_the_same_pooled_member_set_is_called_once_per_family():
    calls = [["Hp", "F", "psv", "0", 0, ["a", "b", "c"], 0.02], ["Hp", "F", "psv", "1", 0, ["a", "b", "c"], 0.03],
             ["Hp", "F", "psv", "1", 1, ["d", "e", "f"], 0.02], ["Gp", "F", "psv", "1", 0, ["a", "b", "c"], 0.02]]
    out = D.dedupe(calls)
    assert [c[:4] for c in out] == [["Hp", "F", "psv", "0"], ["Hp", "F", "psv", "1"], ["Gp", "F", "psv", "1"]]
    assert sorted(out[1][5]) == ["d", "e", "f"]


def test_k_lands_in_its_arms():
    assert D.arms_of("Kp", keep=True) == ["Kn", "K"] and D.arms_of("Kp", keep=False) == ["Kn"]
    a = D.arms_of("S2", keep=True)
    assert "Kn" in a and a.index("Kn") + 1 == a.index("K")


def test_the_integrated_arm_lands_in_its_arms():
    assert D.arms_of("Ip", keep=True) == ["In", "I"] and D.arms_of("Ip", keep=False) == ["In"]
    assert D.arms_of("S2", keep=True)[-2:] == ["In", "I"]


def test_a_member_the_genome_explains_is_dropped_from_an_absent_copy_call():
    # psv call: keep members whose own placement diverges > DELTA; branch call: or carrying >= MIN_NOVEL unexplained bases
    assert D.member_filter("psv", ["a", "b", "c"], {"a": 0.02, "b": 0.003, "c": 0.0}, {"a": 0, "b": 0, "c": 0}) == ["a"]
    assert D.member_filter("branch", ["a", "b"], {"a": 0.001, "b": 0.001}, {"a": D.MIN_NOVEL, "b": 0}) == ["a"]


def test_vg_msga_retries_reverse_then_rotations_of_the_loci_order():
    # round 9 amendment 1: vg msga (v1.73) asserts on KRTAP4's loci in FASTA order (default and -a); other orders build
    assert D.msga_orders(["a", "b", "c"]) == [["c", "b", "a"], ["b", "c", "a"], ["c", "a", "b"]]
    assert D.msga_orders(["a"]) == []


def test_a_family_graph_cache_key_is_the_content_of_its_loci_and_isoform_files(tmp_path):
    for d in ("a", "b"):
        (tmp_path / d).mkdir()
        (tmp_path / d / "loci.fa").write_text(">copy0\nACGT\n")
        (tmp_path / d / "copy0.gtf").write_text("copy0\tsim\texon\t1\t4\n")
    assert D.graph_key(str(tmp_path / "a"), [0]) == D.graph_key(str(tmp_path / "b"), [0])
    (tmp_path / "b" / "copy0.gtf").write_text("copy0\tsim\texon\t1\t3\n")
    assert D.graph_key(str(tmp_path / "a"), [0]) != D.graph_key(str(tmp_path / "b"), [0])


def test_off_the_variable_sites_a_cluster_shows_what_its_reads_show_not_the_host_s_base():
    # round 10 development (set 9, LCE1): one copy's reads clustered on a paralog host; a position covered by a single read is
    # no variable site, and there the read shows its own copy's base -- agreement with that copy, not a difference
    ref = "A" * 1000
    single, shared, allele = set(range(100, 300, 10)), set(range(420, 520, 10)), 700
    maps = {"copy0": {p: "A" for p in range(1000)},
            "copy1": {p: ("G" if p in single | shared else "A") for p in range(1000)}}
    sites = {p: "G" for p in shared}
    sites[allele] = "T"                                               # a donor allele no copy carries
    cons = {p: "alt" for p in sites}

    def read(lo, hi):
        return {p: ("T" if p == allele else maps["copy1"][p]) for p in range(lo, hi)}
    cov = {"m1": read(0, 600), "m2": read(400, 1000), "m3": read(400, 1000)}
    assert D.nearest_copy(cons, sites, ref, maps, set(range(1000)), D.cluster_bases(["m1", "m2", "m3"], cov)) == ("copy1", 0.001)
    assert D.status_of(["m1", "m2", "m3"], cons, sites, ref, maps, cov, {}, mosaic=True)[0] == "unresolved"


def _attribution_case():
    # host = copy0 (all A); copy1 differs from the host at 10..90 step 10 (G); the absent copy at 200..290 step 10 (T)
    ref = "A" * 400
    c1 = set(range(10, 100, 10))
    maps = {"copy0": {p: "A" for p in range(400)}, "copy1": {p: ("G" if p in c1 else "A") for p in range(400)}}
    private = set(range(200, 300, 10))
    sites = {p: "T" for p in private}
    cons = {p: "alt" for p in private}
    return ref, maps, sites, cons, c1, private


def test_an_absent_copy_call_keeps_the_reads_it_explains_better_than_every_catalogued_copy():
    ref, maps, sites, cons, c1, private = _attribution_case()
    cov = {f"a{i}": {p: ("T" if p in private else "A") for p in range(150, 400)} for i in range(3)}   # the absent copy's reads
    cov["p0"] = {p: "A" for p in range(150, 400)}                    # copy0's read swept in: its own copy fits it better
    cov["p1"] = {p: ("G" if p in c1 else "A") for p in range(0, 120)}  # copy1's read: no position tells the call from copy1
    assert D.attribute(sorted(cov), cons, sites, ref, maps, cov) == ["a0", "a1", "a2"]


def test_a_cluster_of_several_present_copies_reads_keeps_none_of_them():
    ref, maps, sites, _, c1, _ = _attribution_case()
    cov = {f"x{i}": {p: "A" for p in range(0, 150)} for i in range(2)}
    cov.update({f"y{i}": {p: ("G" if p in c1 else "A") for p in range(0, 150)} for i in range(2)})
    sites, cons = {p: "G" for p in c1}, {p: "alt" for p in c1}      # the mixture's consensus: copy1's bases
    assert D.attribute(sorted(cov), cons, sites, ref, maps, cov) == []


def test_a_present_copy_s_donor_alleles_do_not_make_its_reads_an_absent_copy_s_when_scaled_by_delta():
    # round 10 (TAF11L.33): reads carrying their own copy's two donor alleles beat it by 2 -- O2's margin, but under delta of
    # a 400-bp read; an absent copy diverged beyond delta beats every catalogued copy by more
    ref, maps, sites, cons, c1, private = _attribution_case()
    two = {200, 210}
    sites2, cons2 = {p: "T" for p in two}, {p: "alt" for p in two}
    cov = {f"d{i}": {p: ("T" if p in two else "A") for p in range(0, 400)} for i in range(3)}
    assert D.attribute(sorted(cov), cons2, sites2, ref, maps, cov) == ["d0", "d1", "d2"]             # margin 2 alone keeps them
    assert D.attribute(sorted(cov), cons2, sites2, ref, maps, cov, delta=D.DELTA) == []
    cov = {f"a{i}": {p: ("T" if p in private else "A") for p in range(0, 400)} for i in range(3)}    # 10 private alts in 400 bp
    assert D.attribute(sorted(cov), cons, sites, ref, maps, cov, delta=D.DELTA) == ["a0", "a1", "a2"]


def test_the_genome_check_subtracts_the_world_s_error_floor_when_given():
    # round 10 (TAF11L.33 hap-2 reads, median de 0.0103): observed divergence = the copy's + sequencing error; uniquely placed
    # reads show the error floor (A119b median de 0.0027-0.0031)
    assert D.genome_keeps("psv", 0.0103, 0)
    assert not D.genome_keeps("psv", 0.0103, 0, floor=0.0029)
    assert D.genome_keeps("psv", 0.0211, 0, floor=0.0029)
    assert D.member_filter("psv", ["a", "b"], {"a": 0.0103, "b": 0.03}, {}, floor=0.0029) == ["b"]


def test_the_error_floor_is_the_median_de_of_uniquely_placed_primaries(tmp_path):
    import pysam
    hdr = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "c1", "LN": 10000}]}
    with pysam.AlignmentFile(str(tmp_path / "u.bam"), "wb", header=hdr) as out:
        for i, (flag, mapq, de) in enumerate([(0, 60, 0.001), (0, 60, 0.003), (0, 60, 0.002), (0, 0, 0.05), (256, 60, 0.09),
                                              (0, 60, 0.004)]):
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.flag, a.reference_id, a.reference_start, a.mapping_quality = f"r{i}", flag, 0, 100 * i, mapq
            a.cigarstring = "100M"
            a.set_tag("de", de, "f")
            out.write(a)
    pysam.sort("-o", str(tmp_path / "E0.bam"), str(tmp_path / "u.bam"))
    pysam.index(str(tmp_path / "E0.bam"))
    assert abs(D.error_floor(str(tmp_path / "E0.bam")) - 0.0025) < 1e-6     # median of 0.001, 0.002, 0.003, 0.004


def test_a_read_no_catalogued_copy_can_be_compared_with_is_not_attributed_to_an_absent_copy():
    # round 10 review (TAF11L.44's host window holds 225 bp of neighbour TAF11L.43; reads aligned there are covered by no copy's
    # map): with no competitor, attribution kept them by default
    ref, maps, sites, cons, c1, private = _attribution_case()
    maps = {n: {p: b for p, b in cc.items() if p < 150} for n, cc in maps.items()}     # the maps end before the call's sites
    cov = {f"a{i}": {p: ("T" if p in private else "A") for p in range(150, 400)} for i in range(3)}
    assert D.attribute(sorted(cov), cons, sites, ref, maps, cov) == []


def test_o1_hosting_places_each_family_read_on_the_seeded_copy_its_primary_overlaps_most(tmp_path):
    # round 10 review: arm I took its hosts and their reads from copy_assign's run; --o1-hosts takes them from O1 (seeded.json)
    import json
    import pysam
    w = tmp_path / "W"
    (w / "graph").mkdir(parents=True)
    seq = "".join("ACGT"[(i * 7 + i // 3) % 4] for i in range(20000))
    (w / "genome.ref.fa").write_text(f">c1\n{seq}\n")
    pysam.faidx(str(w / "genome.ref.fa"))
    copies = [dict(idx=0, contig="c1", start=1000, end=3000, strand="+", rep="t", isoforms={"t": [[1000, 1500], [2000, 3000]]},
                   origin="curated"),
              dict(idx=4, contig="c1", start=10000, end=11000, strand="-", rep="obs0", isoforms={"obs0": [[10000, 11000]]},
                   origin="new")]
    json.dump({"F": dict(copies=copies, dropped=[1], family_reads=["r1", "r2", "r3"])}, open(w / "graph" / "seeded.json", "w"))
    hdr = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "c1", "LN": 20000}]}
    with pysam.AlignmentFile(str(w / "u.bam"), "wb", header=hdr) as out:
        for name, flag, pos, cig in [("r1", 0, 1100, "400M500N800M"), ("r1", 256, 10100, "500M"), ("r2", 16, 10200, "600M"),
                                     ("r3", 0, 15000, "300M")]:
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.flag, a.reference_id, a.reference_start, a.cigarstring = name, flag, 0, pos, cig
            a.mapping_quality = 0 if flag & 256 else 60
            out.write(a)
    pysam.sort("-o", str(w / "E0.bam"), str(w / "u.bam"))
    pysam.index(str(w / "E0.bam"))
    fam, by_host = D.o1_hosts(str(w), D.Genome(str(w / "genome.ref.fa")))
    assert [(c["idx"], c["rs"], c["re"]) for c in fam["F"]] == [(0, 1000, 3000), (4, 10000, 11000)]
    assert fam["F"][0]["seq"] == seq[1000:1500] + seq[2000:3000]
    assert dict(by_host) == {("F", 0): ["r1"], ("F", 4): ["r2"]}                # r3's primary lies on no seeded copy
