import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


F = _load("families")


def test_the_largest_component_of_paralogs_at_the_identity_floor():
    hits = [("g1", "g2", 0.97, 900, 1000, 1000), ("g2", "g3", 0.93, 800, 1000, 1000), ("g4", "g5", 0.99, 900, 1000, 1000),
            ("g1", "g4", 0.85, 900, 1000, 1000), ("g3", "g6", 0.95, 300, 1000, 1000)]
    assert F.largest_component(["g1", "g2", "g3", "g4", "g5", "g6"], hits, min_ident=0.90) == ["g1", "g2", "g3"]


def test_isoforms_keep_the_canonical_first_and_distinct_chains_up_to_the_cap():
    tx = {"t1": [(0, 100), (200, 300)], "t2": [(0, 100), (200, 300), (400, 500)], "t3": [(0, 100), (400, 500)],
          "t4": [(0, 100), (200, 300)], "t5": [(50, 100), (400, 500)], "t6": [(0, 100), (250, 300), (400, 500)]}
    main, iso = F.pick_isoforms(tx, cap=3)
    assert main == "t2" and list(iso) == ["t2", "t1", "t3"]


def test_overlapping_genes_keep_the_one_with_more_transcripts():
    genes = [dict(gene_id="a", chrom="c", start=0, end=100, n=2), dict(gene_id="b", chrom="c", start=50, end=150, n=5),
             dict(gene_id="d", chrom="c", start=200, end=300, n=1)]
    assert [g["gene_id"] for g in F.drop_overlaps(genes)] == ["b", "d"]
