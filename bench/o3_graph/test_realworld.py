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
