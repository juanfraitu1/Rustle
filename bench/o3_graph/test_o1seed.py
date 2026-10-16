import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S = _load("o1seed")


def _aln(read, chrom, start, blocks, strand="+"):
    """A primary alignment: reference blocks (exons) from `start`."""
    return dict(read=read, chrom=chrom, strand=strand, blocks=[(start + a, start + b) for a, b in blocks])


def test_a_catalogued_copy_is_kept_only_with_rna_support():
    loci = {0: ("c1", 1000, 3000), 1: ("c1", 10000, 12000)}
    prim = [_aln(f"r{i}", "c1", 1100, [(0, 300), (800, 1200)]) for i in range(6)] + [_aln("s0", "c1", 10100, [(0, 400)])]
    sup = S.support(prim, loci)
    assert sup == {0: 6, 1: 1}
    assert S.kept(sup, min_support=5) == [0]


def test_primary_alignments_off_the_catalog_form_new_loci():
    loci = {0: ("c1", 1000, 3000)}
    prim = ([_aln(f"r{i}", "c2", 5000 + 10 * i, [(0, 300), (900, 1300)]) for i in range(6)]
            + [_aln(f"x{i}", "c3", 100, [(0, 200)]) for i in range(2)]            # too few reads
            + [_aln("y", "c1", 1500, [(0, 200)])])                               # on the catalog
    new = S.new_loci(prim, loci, min_support=5)
    assert len(new) == 1 and new[0]["chrom"] == "c2" and new[0]["n_reads"] == 6
    assert new[0]["start"] == 5000 and new[0]["end"] == 5050 + 1300 and new[0]["strand"] == "+"


def test_intron_chains_seen_in_enough_reads_become_structures():
    prim = ([_aln(f"r{i}", "c1", 1000 + i, [(0, 300 - i), (800 - i, 1200 - i)]) for i in range(4)]
            + [_aln("z", "c1", 1000, [(0, 300), (600, 1200)])])                  # a chain seen once
    ch = S.observed_chains(prim, min_reads=3)
    assert len(ch) == 1
    (exons, n), = ch
    assert n == 4 and exons == [(1000, 1300), (1800, 2200)]


def test_exons_split_at_introns_and_deletions_stay_inside():
    # 100M 300N 50M 2D 48M: two exons, the deletion inside the second
    assert S.exons_from_cigartuples(1000, [(0, 100), (3, 300), (0, 50), (2, 2), (0, 48)]) == [(1000, 1100), (1400, 1500)]
    assert S.exons_from_cigartuples(0, [(4, 20), (0, 100), (1, 5), (0, 50)]) == [(0, 150)]


def test_a_family_keeps_supported_copies_adds_an_expressed_off_catalog_locus_and_never_seeds_another_family_s_copy():
    cps = [dict(idx=0, contig="c1", start=1000, end=3000, strand="+", rep="t", isoforms={"t": [[1000, 1300], [1800, 2200]]}),
           dict(idx=1, contig="c1", start=10000, end=12000, strand="+", rep="u", isoforms={"u": [[10000, 10400]]})]
    prim = ([_aln(f"a{i}", "c1", 1000, [(0, 300), (800, 1200)]) for i in range(6)]
            + [_aln(f"n{i}", "c2", 5000, [(0, 300), (900, 1300)]) for i in range(5)]       # an uncatalogued locus
            + [_aln(f"o{i}", "c9", 100, [(0, 300), (900, 1300)]) for i in range(5)])        # another family's copy
    all_loci = {("F", 0): ("c1", 0, 4000), ("F", 1): ("c1", 9000, 13000), ("G", 0): ("c9", 0, 3000)}
    copies, dropped = S.seed_family(cps, prim, all_loci)
    assert dropped == [1]                                                          # copy 1: no reads
    assert [c["origin"] for c in copies] == ["curated", "new"] and copies[1]["contig"] == "c2" and copies[1]["idx"] == 2
    assert copies[1]["isoforms"] == {"obs0": [[5000, 5300], [5900, 6300]]}


def test_a_new_locus_needs_reads_placed_there_uniquely():
    loci = {0: ("c1", 1000, 3000)}
    ties = [dict(_aln(f"t{i}", "c2", 5000, [(0, 300), (900, 1300)]), mapq=0) for i in range(8)]   # coin flips
    uniq = [dict(_aln(f"u{i}", "c3", 5000, [(0, 300), (900, 1300)]), mapq=60) for i in range(5)]
    new = S.new_loci(ties + uniq, loci, min_support=5)
    assert [g["chrom"] for g in new] == ["c3"]


def test_a_copy_s_observed_structures_come_from_reads_placed_there_uniquely():
    cps = [dict(idx=0, contig="c1", start=1000, end=3000, strand="+", rep="t", isoforms={"t": [[1000, 1300], [1800, 2200]]})]
    flips = [dict(_aln(f"f{i}", "c1", 1000, [(0, 300), (600, 1200)]), mapq=0) for i in range(4)]     # a paralog's structure
    own = [dict(_aln(f"o{i}", "c1", 1000, [(0, 300), (800, 1200)]), mapq=60) for i in range(5)]
    copies, _ = S.seed_family(cps, flips + own, {("F", 0): ("c1", 0, 4000)})
    assert list(copies[0]["isoforms"]) == ["t"]                      # the coin-flip structure is not the copy's


def test_adjacent_copies_do_not_take_each_other_s_structures():
    # copy 1 ends at 3000, copy 0... copy 2 starts at 3091 (CSAG): a read of copy 2 overlaps copy 1's padded window
    cps = [dict(idx=1, contig="c1", start=1000, end=3000, strand="+", rep="a", isoforms={"a": [[1000, 1300], [2500, 3000]]}),
           dict(idx=2, contig="c1", start=3091, end=6000, strand="+", rep="b", isoforms={"b": [[3091, 3400], [5500, 6000]]})]
    reads2 = [dict(_aln(f"r{i}", "c1", 3091, [(0, 309), (1000, 1300)]), mapq=60) for i in range(5)]   # copy 2's own new chain
    copies, _ = S.seed_family(cps, reads2, {("F", 1): ("c1", 0, 4000), ("F", 2): ("c1", 2091, 7000)})
    by = {c["idx"]: c for c in copies}
    assert list(by[1]["isoforms"]) == ["a"]                                    # nothing of copy 2's
    assert any(n.startswith("obs") for n in by[2]["isoforms"])                 # copy 2 keeps its own observed chain


def test_a_new_locus_counts_any_unique_placement_but_a_copy_s_structures_need_near_unique_reads():
    loci = {0: ("c1", 1000, 3000)}
    # a real uncatalogued copy near-identical to a paralog: unique placements at intermediate mapq (ANKRD20A.4, all 1-49)
    mid_new = [dict(_aln(f"w{i}", "c3", 5000, [(0, 300), (900, 1300)]), mapq=20) for i in range(6)]
    assert [g["chrom"] for g in S.new_loci(mid_new, loci, min_support=5)] == ["c3"]
    cps = [dict(idx=0, contig="c1", start=1000, end=3000, strand="+", rep="t", isoforms={"t": [[1000, 1300], [1800, 2200]]})]
    mid = [dict(_aln(f"m{i}", "c1", 1000, [(0, 300), (600, 1200)]), mapq=30) for i in range(5)]
    copies, _ = S.seed_family(cps, mid, {("F", 0): ("c1", 0, 4000)})
    assert list(copies[0]["isoforms"]) == ["t"]


def test_affinity_counts_aligned_bases_on_the_family_s_copies_not_on_their_padded_windows():
    cps = [dict(idx=0, contig="c1", start=1000, end=3000), dict(idx=1, contig="c2", start=500, end=900)]
    alns = [_aln("r", "c1", 1000, [(0, 300), (800, 1200)]),        # 300 + 400 aligned bases on copy 0 (the intron is no base)
            _aln("r", "c2", 100, [(0, 450)]),                       # 50 bases on copy 1 (500-550)
            _aln("s", "c1", 3005, [(0, 500)])]                      # inside copy 0's padded window only
    assert S.affinity(alns, cps) == {"r": 750}


def _new(idx, contig, start, end, reads):
    return dict(idx=idx, contig=contig, start=start, end=end, strand="+", rep="obs0", isoforms={"obs0": [[start, end]]},
                origin="new", support=len(reads), reads=sorted(reads))


def test_a_new_locus_seeded_in_two_families_stays_only_in_the_family_its_reads_align_to_most():
    # WASH.7 (round 7): an uncatalogued WASH copy 1 kb from a DDX11L copy; its reads also align near DDX11L copies
    reads = [f"w{i}" for i in range(6)]
    seeded = {"WASH": dict(copies=[_new(9, "chr3", 100, 900, reads)], dropped=[]),
              "DDX11L": dict(copies=[_new(10, "chr3", 100, 900, reads)], dropped=[])}
    aff = {"WASH": {q: 2000 for q in reads}, "DDX11L": {q: 40 for q in reads}}
    out = S.one_family(seeded, aff)
    assert [c["idx"] for c in out["WASH"]["copies"]] == [9]
    assert out["DDX11L"]["copies"] == []


def test_a_new_locus_its_reads_tie_between_families_is_seeded_in_none():
    reads = ["a", "b", "c", "d", "e"]
    seeded = {"F": dict(copies=[_new(3, "c1", 100, 900, reads)], dropped=[]),
              "G": dict(copies=[_new(4, "c1", 150, 950, reads)], dropped=[])}
    out = S.one_family(seeded, {"F": {}, "G": {}})
    assert out["F"]["copies"] == [] and out["G"]["copies"] == []


def test_a_new_locus_seeded_in_one_family_and_curated_copies_are_left_alone():
    cur = dict(idx=0, contig="c1", start=100, end=900, strand="+", rep="t", isoforms={"t": [[100, 900]]}, origin="curated",
               support=9)
    seeded = {"F": dict(copies=[cur, _new(1, "c5", 100, 900, ["a"])], dropped=[]),
              "G": dict(copies=[dict(cur, idx=0)], dropped=[])}
    out = S.one_family(seeded, {"F": {}, "G": {}})
    assert [c["idx"] for c in out["F"]["copies"]] == [0, 1] and [c["idx"] for c in out["G"]["copies"]] == [0]


def test_a_new_copy_records_the_reads_that_seeded_it():
    cps = [dict(idx=0, contig="c1", start=1000, end=3000, strand="+", rep="t", isoforms={"t": [[1000, 1300], [1800, 2200]]})]
    prim = [_aln(f"n{i}", "c2", 5000, [(0, 300), (900, 1300)]) for i in range(5)]
    copies, _ = S.seed_family(cps, prim, {("F", 0): ("c1", 0, 4000)})
    new = [c for c in copies if c["origin"] == "new"]
    assert len(new) == 1 and new[0]["reads"] == [f"n{i}" for i in range(5)]
