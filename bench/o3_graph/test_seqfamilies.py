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
