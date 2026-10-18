import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_seqfamilies", os.path.join(HERE, "seqfamilies.py"))
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)


def test_a_gene_s_model_is_carried_by_another_gene_s_genome_when_the_alignment_s_introns_are_its_junctions():
    # gene A: exons of 100, 50 and 80 bp (junctions at 100 and 150 of its transcript), spliced-aligned to B's locus from 1,000
    cg = [(100, "M"), (300, "N"), (50, "M"), (400, "N"), (80, "M")]
    assert S.structure_match([100, 150], 230, "+", 1000, cg)
    assert S.structure_match([95, 155], 230, "+", 1000, cg)                                     # within 10 bp
    assert not S.structure_match([100, 120], 230, "+", 1000, cg)                                # an intron at no junction
    assert not S.structure_match([100, 150, 200], 230, "+", 1000, cg)                           # a junction B's genome lacks
    # a processed retrocopy (one exon) spliced onto its parent's locus: introns the copy does not have
    assert not S.structure_match([], 230, "+", 1000, cg)
    # B's genome holds A's two exons contiguously (OR4F21 onto OR4F16: 398M1D596M): A's intron is not there
    assert not S.structure_match([398], 994, "+", 1000, [(398, "M"), (1, "D"), (596, "M")])
    # minus strand: the query runs backwards along the target, its junctions counted from its other end
    assert S.structure_match([80, 130], 230, "-", 1000, cg)
    # two single-exon genes
    assert S.structure_match([], 230, "+", 1000, [(230, "M")])


def test_families_are_the_components_of_links_with_both_sequence_and_structure_and_avoid_spent_genes():
    links = [("a", "b", True), ("b", "c", True), ("c", "d", False), ("e", "f", True), ("g", "h", True)]
    fams = S.components(links, spent={"h"})
    assert fams == [["a", "b", "c"], ["e", "f"]]                    # c-d fails structure; g-h touches a spent gene


def test_the_draw_takes_families_in_seeded_order_up_to_the_copy_budget():
    fams = [dict(family=f"F{i}", copies=[{}] * k) for i, k in enumerate([2, 3, 2, 40, 2, 5, 2, 2])]
    got = S.draw(fams, seed=17, budget=10, max_total=20)
    assert sum(len(f["copies"]) for f in got) >= 10 and sum(len(f["copies"]) for f in got) <= 20
    assert all(len(f["copies"]) < 40 for f in got)                  # the 40-copy family would exceed max_total
    assert got == S.draw(fams, seed=17, budget=10, max_total=20)    # deterministic


def test_a_junction_bordering_an_unaligned_end_of_the_transcript_needs_no_intron():
    # LGALS7 (round 17): its main transcript ends in a 6-bp exon minimap2 cannot align -- the twin's alignment covers 460 of
    # 466 bp with 2 introns against 3 junctions. A junction at the edge of an unaligned end is no structural difference.
    cg = [(169, "M"), (689, "N"), (202, "M"), (984, "N"), (89, "M")]
    assert S.structure_match([169, 371, 460], 466, "+", 1000, cg)              # junction 460: the aligned part ends there
    assert S.structure_match([6, 175, 377], 466, "+", 1000, cg, q0=6)           # an unaligned 5' piece of 6 bp
    assert not S.structure_match([169, 300, 371, 460], 466, "+", 1000, cg)     # a junction inside the aligned part: missing


def _g(gid, chrom, start, end, n=1, biotype="protein_coding", name=None):
    return dict(gene_id=gid, gene_name=name or gid, chrom=chrom, start=start, end=end, strand="+", n=n, biotype=biotype)


def test_paralogs_overlapping_a_coding_gene_are_dropped_and_among_overlapping_paralogs_a_pseudogene_wins():
    coding = [_g("OR4F16", "chr5", 1000, 2000)]
    nc = [_g("OR4F1P", "chr5", 1019, 1958, biotype="unprocessed_pseudogene"),         # annotated inside OR4F16: dropped
          _g("lnc1", "chr6", 100, 5000, n=4, biotype="lncRNA"),                       # holds OR4F7P: the pseudogene is the copy
          _g("OR4F7P", "chr6", 3000, 3900, n=1, biotype="unprocessed_pseudogene"),
          _g("OR4F2P", "chr11", 100, 1000, biotype="unprocessed_pseudogene"),
          _g("lnc2", "chr12", 100, 900, n=1, biotype="lncRNA"),                       # two overlapping lncRNAs: more transcripts
          _g("lnc3", "chr12", 500, 1500, n=3, biotype="lncRNA")]
    assert {g["gene_id"] for g in S.add_paralog_genes(coding, nc)} == {"OR4F7P", "OR4F2P", "lnc3"}


def test_a_paralog_links_to_a_coding_gene_only_when_its_locus_carries_the_coding_gene_s_model():
    coding = {"A", "B"}
    cover = {("A", "B"): 0.9, ("B", "A"): 0.3, ("A", "P"): 0.95, ("P", "A"): 1.0, ("Q", "A"): 1.0, ("A", "Q"): 0.2,
             ("P", "Q"): 1.0, ("Q", "P"): 1.0}
    struct = {("A", "B"): False, ("B", "A"): True, ("A", "P"): True, ("P", "A"): False, ("Q", "A"): True, ("P", "Q"): True,
              ("Q", "P"): True}
    assert S.pair_link("A", "B", cover, struct, coding) is True        # coding pair: either direction, as before
    assert S.pair_link("A", "P", cover, struct, coding) is True        # P's locus carries A's model
    assert S.pair_link("P", "A", cover, struct, coding) is True        # the same pair, given the other way round
    assert S.pair_link("A", "Q", cover, struct, coding) is None        # Q's own model sits in A's locus: not a copy of A
    assert S.pair_link("P", "Q", cover, struct, coding) is None        # two paralogs: no link of their own
    struct[("A", "P")] = False
    assert S.pair_link("A", "P", cover, struct, coding) is False       # a sequence link without A's exon chain


def test_a_paralog_family_is_named_after_its_first_coding_gene():
    gs = [_g("P1", "chr1", 10, 20, biotype="unprocessed_pseudogene", name="OR4F2P"), _g("C1", "chr1", 50, 60, name="OR4F21"),
          _g("C2", "chr5", 5, 9, name="OR4F16")]
    used = {}
    assert S.family_name(gs, used, "S19") == "S19OR4F21"
    assert S.family_name(gs, used, "S19") == "S19OR4F21x2"


def test_the_draw_reads_a_pool_of_any_set_name(tmp_path):
    import json
    pool = dict(genome="g", sets={"seq19": [dict(family=f"F{i}", copies=[{}] * 2) for i in range(60)]})
    json.dump(pool, open(tmp_path / "all.json", "w"))
    S.cmd_draw(str(tmp_path / "set.json"), str(tmp_path / "all.json"), seed=19, set_name="r19")
    got = json.load(open(tmp_path / "set.json"))
    assert list(got["sets"]) == ["r19"] and sum(len(f["copies"]) for f in got["sets"]["r19"]) == 90
