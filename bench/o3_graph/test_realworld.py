import importlib.util
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R = _load("realworld")
GENOME = "".join(random.Random(3).choice("ACGT") for _ in range(2000))


def test_snps_sit_on_genomic_positions_so_every_isoform_agrees():
    snps = {105: "A", 410: "T"}
    iso1 = R.splice(GENOME, [(100, 200), (400, 500)], "+", snps)
    iso2 = R.splice(GENOME, [(400, 500)], "+", snps)
    assert iso1[5] == "A" and iso1[100 + 10] == "T" and iso2[10] == "T"


def test_minus_strand_transcripts_are_reverse_complemented():
    s = R.splice(GENOME, [(100, 110), (200, 210)], "-", {})
    comp = str.maketrans("ACGT", "TGCA")
    assert s == (GENOME[100:110] + GENOME[200:210])[::-1].translate(comp)


def test_the_novel_exon_goes_into_every_isoform_whose_intron_spans_it():
    novel = "N" * 7
    p = 250                                                     # inside the main isoform's intron 200-300
    with_j = R.splice(GENOME, [(100, 200), (300, 400), (500, 600)], "+", {}, insert_at=p, novel=novel)
    assert with_j == GENOME[100:200] + novel + GENOME[300:400] + GENOME[500:600]
    skipping = R.splice(GENOME, [(100, 200), (500, 600)], "+", {}, insert_at=p, novel=novel)
    assert skipping == GENOME[100:200] + novel + GENOME[500:600]          # a real new exon is spliced in here too
    upstream = R.splice(GENOME, [(500, 600), (700, 800)], "+", {}, insert_at=p, novel=novel)
    assert novel not in upstream                                          # no intron of this isoform spans it


def test_anonymous_names_are_a_seeded_shuffle_carrying_no_truth():
    rows = [("FAM.0", "t1"), ("FAM.0", "t2"), ("FAM.1", "t1")]
    a, b = R.anon_names(rows, seed=4), R.anon_names(rows, seed=4)
    assert a == b and sorted(a) == ["s0000001", "s0000002", "s0000003"]


def test_the_novel_flag_sees_the_exon_on_either_strand():
    novel = "ACGTTGCAAC"
    assert R.carries(novel, "TT" + novel + "GG") and R.carries(novel, "TT" + novel[::-1].translate(R.COMP) + "GG")
    assert not R.carries(novel, "TTGG") and not R.carries("", "ACGT")


def test_exon_blocks_follow_the_cigar_and_split_only_at_introns():
    assert R.exons_from_cigar(1000, "100M300N200M") == [(1000, 1100), (1400, 1600)]
    # a deletion stays inside its exon, an insertion takes no reference
    assert R.exons_from_cigar(1000, "50M2D48M300N10I190M") == [(1000, 1100), (1400, 1590)]
    assert R.exons_from_cigar(0, "10=1X9=") == [(0, 20)]


def test_junctions_are_canonical_on_the_transcript_strand():
    g = "A" * 100 + "GT" + "C" * 96 + "AG" + "A" * 100                  # an intron 100-200 with GT..AG on the plus strand
    assert R.canonical_junctions(g, [(0, 100), (200, 300)], "+")
    assert not R.canonical_junctions(g, [(0, 100), (200, 300)], "-")  # read as minus it is CT..AC-less
    m = "A" * 100 + "CT" + "G" * 96 + "AC" + "A" * 100                  # GT..AG on the minus strand
    assert R.canonical_junctions(m, [(0, 100), (200, 300)], "-")
    assert R.canonical_junctions(g, [(0, 300)], "+")                    # one exon: nothing to check


def test_a_copy_expresses_only_paralog_structures_its_own_catalog_entry_lacks():
    own = {"t1": [(0, 100), (200, 300)]}
    proj = {"p|a": [(0, 100), (200, 300)],                             # the copy's own intron chain: skipped
            "p|b": [(0, 100), (250, 300)], "q|c": [(0, 100), (250, 300)],  # one new chain, offered twice: kept once
            "q|d": None,                                               # did not project
            "q|e": [(0, 50), (200, 300)],
            "q|f": [(0, 103), (198, 300)]}                             # the own chain shifted 2-3 bp (a junction artifact)
    got = R.cross_structures(own, proj, max_proj=5, rng=random.Random(1))
    assert sorted(got) == ["p|b", "q|e"] or sorted(got) == ["q|c", "q|e"]
    assert len(R.cross_structures(own, proj, max_proj=1, rng=random.Random(1))) == 1


def test_a_paralog_isoform_projects_onto_a_copy_locus_through_its_splice_sites():
    rng = random.Random(7)
    rnd = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    e1, e2, e3 = rnd(300), rnd(250), rnd(300)
    locus = rnd(500) + e1 + "GT" + rnd(596) + "AG" + e2 + "GT" + rnd(796) + "AG" + e3 + rnd(500)
    para = list(e1 + e3)                                               # the paralog's isoform skips exon 2, 1% diverged
    for p in rng.sample(range(len(para)), 6):
        para[p] = "A" if para[p] != "A" else "C"
    ex = R.project("".join(para), locus, lo=10000, strand="+", n_exons=2)
    s3 = 500 + 300 + 600 + 250 + 800
    assert ex == [(10500, 10800), (10000 + s3, 10000 + s3 + 300)]
    assert R.project("".join(para), locus, lo=10000, strand="+", n_exons=3) is None   # not the paralog's structure


def test_a_minus_strand_copy_takes_the_projection_on_its_own_strand():
    rng = random.Random(7)
    rnd = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    e1, e2, e3 = rnd(300), rnd(250), rnd(300)
    plus = rnd(500) + e1 + "GT" + rnd(596) + "AG" + e2 + "GT" + rnd(796) + "AG" + e3 + rnd(500)
    minus = plus[::-1].translate(R.COMP)                               # the same gene on the minus strand
    ex = R.project(e1 + e3, minus, lo=0, strand="-", n_exons=2)
    L = len(plus)
    s3 = 500 + 300 + 600 + 250 + 800
    assert ex == [(L - s3 - 300, L - s3), (L - 800, L - 500)]
    assert R.project(e1 + e3, minus, lo=0, strand="+", n_exons=2) is None   # the copy's strand disagrees


def test_positions_lift_through_an_alignment_cigar():
    lift = R.lift_positions("100M5I50M3D100M", qs=0, ts=1000, positions=[0, 99, 102, 105, 154, 155, 254, 300])
    assert lift == {0: 1000, 99: 1099, 102: 1100, 105: 1100, 154: 1149, 155: 1153, 254: 1252, 300: None}


def test_a_paralog_isoform_lifts_over_through_the_locus_alignment():
    rng = random.Random(11)
    rnd = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    e1, e2, e3 = rnd(300), rnd(250), rnd(300)
    i1, i2 = "GT" + rnd(596) + "AG", "GT" + rnd(796) + "AG"
    c_locus = rnd(500) + e1 + i1 + e2 + i2 + e3 + rnd(500)
    p = list(c_locus)
    for k in rng.sample(range(len(p)), 30):                          # the paralog: ~1% diverged, a 4-bp deletion in intron 1
        p[k] = "A" if p[k] != "A" else "C"
    p = "".join(p)
    cut = 500 + 300 + 200
    p_locus = p[:cut] + p[cut + 4:]
    s3 = 500 + 300 + 600 + 250 + 800
    p_exons = [(500, 800), (s3 - 4, s3 - 4 + 300)]                  # exon 3 sits 4 bp earlier on the paralog
    assert R.liftover(p_locus, [p_exons], c_locus) == [[(500, 800), (s3, s3 + 300)]]


def test_heterozygous_sites_are_drawn_from_the_copy_s_exonic_bases_with_a_different_base():
    rng = random.Random(5)
    pos = list(range(100, 2100))
    het = R.het_sites(pos, 0.01, rng, lambda p: GENOME[p % len(GENOME)])
    assert 10 <= len(het) <= 30 and set(het) <= set(pos)
    assert all(b != GENOME[p % len(GENOME)] and b in "ACGT" for p, b in het.items())
    assert R.het_sites(pos, 0.0, rng, lambda p: "A") == {}


def test_a_second_haplotype_carries_the_alleles_on_every_isoform():
    het = {105: "A", 410: "T"}
    hap2 = R.splice(GENOME, [(100, 200), (400, 500)], "+", het)
    hap1 = R.splice(GENOME, [(100, 200), (400, 500)], "+", {})
    assert hap2[5] == "A" and hap2[110] == "T" and hap1[5] == GENOME[105] and hap1[110] == GENOME[410]


def test_a_world_with_heterozygous_alleles_builds_and_records_each_read_s_haplotype(tmp_path, monkeypatch):
    import json
    import pysam
    g = "".join(random.Random(9).choice("ACGT") for _ in range(30000))
    fa = tmp_path / "g.fa"
    fa.write_text(">c1\n" + g + "\n")
    pysam.faidx(str(fa))
    monkeypatch.setattr(R, "GENOME", str(fa))
    cp = lambda gid, lo: dict(gene_id=gid, gene_name=gid, chrom="c1", start=lo, end=lo + 2000, strand="+", main="t",  # noqa: E731
                              transcripts={"t": [[lo, lo + 600], [lo + 1200, lo + 2000]]})
    fams = {"sets": {"x": [dict(family="F", copies=[cp("g0", 1000), cp("g1", 10000), cp("g2", 20000)])]}}
    (tmp_path / "f.json").write_text(json.dumps(fams))
    R.main(["build", str(tmp_path / "w"), str(tmp_path / "f.json"), "--set", "x", "--scenario", "INTACT", "--seed", "3",
            "--het", "0.01", "--per-copy", "20"])
    rows = (tmp_path / "w" / "reads.truth.E0.tsv").read_text().splitlines()
    assert rows[0].split("\t")[-1] == "hap"
    haps = {r.split("\t")[-1] for r in rows[1:]}
    assert haps == {"1", "2"}
    assert sum(json.load(open(tmp_path / "w" / "het.json"))["sites"].values()) == 3 * 14   # 1% of each copy's 1,400 exonic bases


def test_a_silent_world_keeps_every_copy_in_the_catalog_but_a_third_express_nothing(tmp_path, monkeypatch):
    import json
    import pysam
    g = "".join(random.Random(9).choice("ACGT") for _ in range(70000))
    fa = tmp_path / "g.fa"
    fa.write_text(">c1\n" + g + "\n")
    pysam.faidx(str(fa))
    monkeypatch.setattr(R, "GENOME", str(fa))
    cp = lambda gid, lo: dict(gene_id=gid, gene_name=gid, chrom="c1", start=lo, end=lo + 2000, strand="+", main="t",  # noqa: E731
                              transcripts={"t": [[lo, lo + 600], [lo + 1200, lo + 2000]]})
    fams = {"sets": {"x": [dict(family="F", copies=[cp(f"g{k}", 1000 + 10000 * k) for k in range(6)])]}}
    (tmp_path / "f.json").write_text(json.dumps(fams))
    R.main(["build", str(tmp_path / "w"), str(tmp_path / "f.json"), "--set", "x", "--scenario", "SILENT", "--seed", "3",
            "--per-copy", "20"])
    kinds = json.load(open(tmp_path / "w" / "kinds.json"))
    silent = sorted(c for c, k in kinds.items() if k == "silent")
    assert len(silent) == 2                                                    # a third of 6
    copies = (tmp_path / "w" / "P.fam.copies.tsv").read_text().splitlines()[1:]
    assert len(copies) == 6                                                    # all still catalogued
    expressed = {r.split("\t")[1] for r in (tmp_path / "w" / "reads.truth.E0.tsv").read_text().splitlines()[1:]}
    assert expressed.isdisjoint(silent) and len(expressed) == 4


def test_homopolymer_runs_are_found():
    assert R.hp_runs("ACGTAAAAACCCG") == [(4, 9), (9, 12)]


def test_edits_apply_against_the_original_coordinates():
    assert R.apply_edits("ACGTACGT", [(2, "S", "T"), (5, "D", 1), (0, "I", "GG")]) == "GGACTTAGT"


def test_homopolymer_indels_land_inside_runs():
    rng = random.Random(5)
    seq = "".join(rng.choice("ACGT") for _ in range(300)) + "AAAAAA" + "".join(rng.choice("CGT") for _ in range(300))
    runs = R.hp_runs(seq)
    prof = dict(sub_error=0.0, ins_event=0.0, del_event=0.02, del_hp_frac=1.0, ins_hp_frac=1.0, del_len={"1": 1.0},
                ins_len={"1": 1.0})
    edits = R.draw_edits(seq, random.Random(1), prof)
    assert edits and all(k == "D" and any(s <= p < e for s, e in runs) for p, k, _ in edits)


def test_truncation_removes_transcript_fractions_from_each_end_but_keeps_a_minimum():
    assert len(R.truncate("A" * 1000, 0.1, 0.2, min_len=100)) == 700
    assert R.truncate("ACGT" * 50, 0.6, 0.3, min_len=100) == "ACGT" * 50


def test_a_retained_intron_merges_its_two_exons():
    assert R.retain([(0, 10), (20, 30), (40, 50)], {1}) == [(0, 10), (20, 50)]
    assert R.retain([(0, 10), (20, 30)], set()) == [(0, 10), (20, 30)]


def test_expression_multipliers_are_deterministic_and_flat_without_spread():
    assert R.expr_mult(3, "F.1", 0.0) == 1.0
    assert R.expr_mult(3, "F.1", 0.86) == R.expr_mult(3, "F.1", 0.86) != R.expr_mult(3, "F.2", 0.86)


def test_a_truncation_that_would_leave_too_little_is_redrawn_not_dropped():
    pairs = [[0.9, 0.0], [0.0, 0.95], [0.1, 0.1]]           # only the last leaves >= min_len of a 400-base transcript
    out = R.draw_truncation("A" * 400, pairs, random.Random(0), min_len=200)
    assert len(out) == 400 - 40 - 40
    assert R.draw_truncation("A" * 400, [[0.9, 0.0]], random.Random(0), min_len=200, tries=3) == "A" * 400


def test_adapters_are_attached_at_drawn_clip_lengths():
    out = R.attach_adapters("CCCCCC", "AAATTT", "GGGTTT", (2, 4))
    assert out == "TT" + "CCCCCC" + "GGGT"
    assert R.attach_adapters("CC", "AAA", "GGG", (5, 0)) == "AAA" + "CC"     # a clip longer than the adapter: all of it


def test_a_projection_may_lack_a_terminal_exon_too_short_for_the_aligner_when_asked():
    # round 10 development (WASH): a 6-bp last exon cannot be spliced-aligned, so the projection had one exon too few
    rng = random.Random(11)
    rnd = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    e1, e2, e3 = rnd(300), rnd(250), "TCCAGA"
    locus = rnd(500) + e1 + "GT" + rnd(596) + "AG" + e2 + "GT" + rnd(796) + "AG" + e3 + rnd(500)
    assert R.project(e1 + e2 + e3, locus, lo=0, strand="+", n_exons=3) is None
    ex = R.project(e1 + e2 + e3, locus, lo=0, strand="+", n_exons=[300, 250, 6], micro=30)
    assert ex is not None and len(ex) == 2 and ex[0] == (500, 800) and ex[1][0] == 1400
    assert R.project(e1 + e2 + e3, locus, lo=0, strand="+", n_exons=[300, 250, 6], micro=5) is None


def test_structures_a_few_bases_apart_at_a_junction_stay_apart_when_asked():
    own = {"t1": [(0, 100), (200, 300)]}
    proj = {"q|f": [(0, 103), (198, 300)]}                            # 2-3 bp from the own chain
    assert R.cross_structures(own, proj, max_proj=None) == {}
    assert R.cross_structures(own, proj, max_proj=None, tol=0) == {"q|f": [(0, 103), (198, 300)]}


def test_a_copy_s_heterozygous_sites_can_be_drawn_rather_than_rounded():
    # round 13 review I1: round(2 x het_per_base x exonic bp) gave 86 of 89 public copies no het site and none two -- an
    # effectively haploid donor; drawn per base, the same rate gives the copies their expected spread
    import random
    pos, rate = range(1500), 0.0002
    rounded = [len(R.het_sites(pos, rate, random.Random(s), lambda p: "A")) for s in range(300)]
    drawn = [len(R.het_sites(pos, rate, random.Random(s), lambda p: "A", draw=True)) for s in range(300)]
    assert set(rounded) == {0}
    assert 0.15 < sum(drawn) / len(drawn) < 0.45 and max(drawn) >= 2
    assert R.het_sites(pos, rate, random.Random(7), lambda p: "A", draw=True) == R.het_sites(pos, rate, random.Random(7),
                                                                                              lambda p: "A", draw=True)


def test_a_paralog_is_silent_at_the_measured_fraction_or_expressed_at_the_measured_ratio():
    prof = dict(paralog_silent_frac=0.25, paralog_expr_ratio=0.3)
    assert R.paralog_mult(3, "F.1", "protein_coding", prof) == 1.0                 # coding copies: unchanged
    assert R.paralog_mult(3, "F.1", "unprocessed_pseudogene", {}) == 1.0           # a profile without the measurement
    assert R.paralog_mult(3, "F.1", "unprocessed_pseudogene", dict(prof, paralog_silent_frac=0.0)) == 0.3
    assert R.paralog_mult(3, "F.1", "lncRNA", dict(prof, paralog_silent_frac=1.0)) == 0.0
    draws = [R.paralog_mult(3, f"F.{k}", "unprocessed_pseudogene", prof) for k in range(4000)]
    assert draws == [R.paralog_mult(3, f"F.{k}", "unprocessed_pseudogene", prof) for k in range(4000)]   # deterministic
    assert set(draws) == {0.0, 0.3} and abs(draws.count(0.0) / 4000 - 0.25) < 0.03


def _paralog_world(tmp_path, monkeypatch, scenario, prof=None):
    import json
    import pysam
    g = "".join(random.Random(9).choice("ACGT") for _ in range(70000))
    fa = tmp_path / "g.fa"
    fa.write_text(">c1\n" + g + "\n")
    pysam.faidx(str(fa))
    monkeypatch.setattr(R, "GENOME", str(fa))

    def cp(gid, lo, bt, n_ex):
        ex = [[lo + 400 * k, lo + 400 * k + 200] for k in range(n_ex)]
        return dict(gene_id=gid, gene_name=gid, chrom="c1", start=lo, end=ex[-1][1], strand="+", main="t", transcripts={"t": ex},
                    biotype=bt)
    copies = [cp(f"c{k}", 1000 + 10000 * k, "protein_coding", 2) for k in range(3)]
    copies += [cp(f"p{k}", 31000 + 10000 * k, "unprocessed_pseudogene", 4) for k in range(3)]   # more exons than the coding
    (tmp_path / "f.json").write_text(json.dumps({"sets": {"x": [dict(family="F", copies=copies)]}}))
    extra = []
    if prof is not None:
        (tmp_path / "p.json").write_text(json.dumps(prof))
        extra = ["--sim-profile", str(tmp_path / "p.json")]
    R.main(["build", str(tmp_path / "w"), str(tmp_path / "f.json"), "--set", "x", "--scenario", scenario, "--seed", "3",
            "--per-copy", "20"] + extra)
    kinds = json.load(open(tmp_path / "w" / "kinds.json"))
    reads = {r.split("\t")[1] for r in (tmp_path / "w" / "reads.truth.E0.tsv").read_text().splitlines()[1:]}
    return kinds, reads


def test_a_scenario_s_copy_is_drawn_among_the_coding_copies():
    cp = lambda gid, bt, n_ex: dict(gene_id=gid, biotype=bt, main="t", transcripts={"t": [[0, 1]] * n_ex})  # noqa: E731
    fam = dict(family="F", copies=[cp("c0", "protein_coding", 2), cp("c1", "protein_coding", 2),
                                   cp("p0", "unprocessed_pseudogene", 4), cp("p1", "lncRNA", 3)])
    assert {R.special_copy(fam, s)["gene_id"] for s in range(40)} == {"c0", "c1"}


def test_a_silent_world_silences_coding_copies_and_leaves_its_paralogs_to_the_measured_draw(tmp_path, monkeypatch):
    kinds, reads = _paralog_world(tmp_path, monkeypatch, "SILENT")
    assert [c for c, k in kinds.items() if k == "silent"] in (["F.0"], ["F.1"], ["F.2"])   # a third of the 3 coding copies
    assert {"F.3", "F.4", "F.5"} <= reads                                                # no profile: paralogs as coding


def test_a_paralog_drawn_silent_is_catalogued_and_has_no_read(tmp_path, monkeypatch):
    import json
    prof = json.load(open(os.path.join(HERE, "profiles", "A119b_testis_v2.json")))
    prof.update(paralog_silent_frac=1.0, paralog_expr_ratio=0.5)
    kinds, reads = _paralog_world(tmp_path, monkeypatch, "INTACT", prof)
    assert kinds == {"F.0": "present", "F.1": "present", "F.2": "present", "F.3": "silent", "F.4": "silent", "F.5": "silent"}
    assert reads == {"F.0", "F.1", "F.2"}
    copies = (tmp_path / "w" / "P.fam.copies.tsv").read_text().splitlines()[1:]
    assert len(copies) == 6


def test_an_expressed_paralog_has_its_measured_share_of_a_coding_copy_s_reads(tmp_path, monkeypatch):
    import collections
    import json
    prof = json.load(open(os.path.join(HERE, "profiles", "A119b_testis_v2.json")))
    prof.update(paralog_silent_frac=0.0, paralog_expr_ratio=0.5, expr_log_sd=0.0)
    kinds, _ = _paralog_world(tmp_path, monkeypatch, "INTACT", prof)
    assert set(kinds.values()) == {"present"}
    n = collections.Counter(r.split("\t")[1] for r in (tmp_path / "w" / "reads.truth.E0.tsv").read_text().splitlines()[1:])
    assert [n[f"F.{k}"] for k in range(6)] == [20, 20, 20, 10, 10, 10]


def test_a_branch_world_with_paralogs_builds_and_its_branch_copy_is_expressed_as_a_coding_copy(tmp_path, monkeypatch):
    import collections
    import json
    prof = json.load(open(os.path.join(HERE, "profiles", "A119b_testis_v2.json")))
    prof.update(paralog_silent_frac=0.0, paralog_expr_ratio=0.5, expr_log_sd=0.0)
    kinds, reads = _paralog_world(tmp_path, monkeypatch, "BRANCH", prof)
    assert "F.B" in kinds and "F.B" in reads                         # the extra absent copy, from a coding copy
    n = collections.Counter(r.split("\t")[1] for r in (tmp_path / "w" / "reads.truth.E0.tsv").read_text().splitlines()[1:])
    assert n["F.B"] >= n["F.3"]                                       # not scaled down as a paralog


def test_a_world_s_regions_end_at_their_contig_s_end(tmp_path, monkeypatch):
    # round 20's draw holds chrM genes (16,569 bp) with their nuclear copies: a copy's span +- REGION_PAD runs past the contig
    _paralog_world(tmp_path, monkeypatch, "INTACT")
    for ln in (tmp_path / "w" / "regions.txt").read_text().split():
        ch, rng = ln.split(":")
        s, e = map(int, rng.split("-"))
        assert 1 <= s <= e <= 70000, ln                              # the fixture's contig is 70 kb


def test_a_world_can_take_real_reads_and_a_named_copy_to_remove(tmp_path, monkeypatch):
    # the real-read masking test (2026-10-18): the reads and their truth come from a real library; the copy to remove is
    # named (realreads.py choose), not drawn; nothing is simulated
    import json
    import pysam
    g = "".join(random.Random(9).choice("ACGT") for _ in range(70000))
    fa = tmp_path / "g.fa"
    fa.write_text(">c1\n" + g + "\n")
    pysam.faidx(str(fa))
    monkeypatch.setattr(R, "GENOME", str(fa))
    cp = lambda gid, lo: dict(gene_id=gid, gene_name=gid, chrom="c1", start=lo, end=lo + 2000, strand="+", main="t",  # noqa: E731
                              transcripts={"t": [[lo, lo + 600], [lo + 1200, lo + 2000]]}, biotype="protein_coding")
    (tmp_path / "f.json").write_text(json.dumps({"sets": {"x": [dict(family="F", copies=[cp(f"g{k}", 1000 + 10000 * k) for k in range(3)])]}}))
    (tmp_path / "r.fq").write_text("".join(f"@m/{i}/ccs\n{g[1000 + i:1500 + i]}\n+\n{'~' * 500}\n" for i in range(4)))
    (tmp_path / "t.tsv").write_text("read\tcopy\tfamily\nm/0/ccs\tF.1\tF\nm/1/ccs\tF.1\tF\nm/2/ccs\tF.0\tF\n")
    (tmp_path / "s.json").write_text(json.dumps({"F": "g1"}))
    R.main(["build", str(tmp_path / "w"), str(tmp_path / "f.json"), "--set", "x", "--scenario", "REMOVED", "--seed", "3",
            "--real-reads", str(tmp_path / "r.fq"), "--real-truth", str(tmp_path / "t.tsv"), "--special-copies", str(tmp_path / "s.json")])
    kinds = json.load(open(tmp_path / "w" / "kinds.json"))
    assert kinds == {"F.0": "present", "F.1": "absent", "F.2": "present"}
    rows = [r.split("\t") for r in (tmp_path / "w" / "reads.truth.E0.tsv").read_text().splitlines()]
    assert rows[0] == ["read", "copy", "kind", "family", "isoform", "novel", "hap"]
    assert sorted(r[:3] for r in rows[1:]) == [["m/0/ccs", "F.1", "absent"], ["m/1/ccs", "F.1", "absent"], ["m/2/ccs", "F.0", "present"]]
    assert (tmp_path / "w" / "reads.E0.fq").read_text() == (tmp_path / "r.fq").read_text()          # every read, as given
    assert len((tmp_path / "w" / "P.fam.copies.tsv").read_text().splitlines()) == 1 + 2             # the removed copy: not catalogued
