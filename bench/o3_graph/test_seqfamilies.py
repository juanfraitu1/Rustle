import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_seqfamilies", os.path.join(HERE, "seqfamilies.py"))
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)


def test_two_genes_share_a_model_when_the_alignment_s_introns_are_both_genes_introns():
    # gene A: exons of 100, 50 and 80 bp (junctions at 100 and 150 of its transcript), aligned to B's locus from 1,000 with
    # introns of 300 and 400 bp -- B's own introns
    cg = [(100, "M"), (300, "N"), (50, "M"), (400, "N"), (80, "M")]
    b_introns = [(1100, 1400), (1450, 1850)]
    assert S.structure_match([100, 150], 230, "+", 1000, cg, b_introns)
    assert S.structure_match([100, 150], 230, "+", 1000, cg, [(1105, 1400), (1450, 1845)])        # within 10 bp
    assert not S.structure_match([100, 150], 230, "+", 1000, cg, [(1130, 1400), (1450, 1850)])    # a moved junction
    assert not S.structure_match([100, 150], 230, "+", 1000, cg, [(1100, 1400)])                  # B lacks an intron
    # a processed retrocopy (one exon) spliced onto its parent's locus: introns the copy does not have
    assert not S.structure_match([], 230, "+", 1000, cg, b_introns)
    # minus strand: the query runs backwards along the target, its junctions counted from its other end
    assert S.structure_match([80, 130], 230, "-", 1000, cg, b_introns)
    # two single-exon genes
    assert S.structure_match([], 230, "+", 1000, [(230, "M")], [])


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
