import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


E = _load("evidence")


def test_a_clipped_end_costs_its_edits_against_the_adjacent_path_not_its_length():
    # the read's first 4 bases carry one error; aligned against the 4 path bases before the core they cost 1 edit, not 4
    assert E.end_cost("TCGA", "GGGGACGA", anchored_at_end=True) == 1
    assert E.end_cost("ACGA", "GGGGACGA", anchored_at_end=True) == 0
    # a suffix: anchored at the path's start, free beyond
    assert E.end_cost("TTGC", "TTGCAAAA", anchored_at_end=False) == 0
    # past the path's end every remaining read base costs one
    assert E.end_cost("ACGTACGT", "ACG", anchored_at_end=False) == 5


def test_semiglobal_distance_replaces_clipped_bases_by_their_edits():
    core = "ACGTTGCAGGCTAGCTAGGATCCGATCGATTACG"
    target = "GGGG" + "ACGA" + core + "CCCC"
    read = "TCGA" + core                       # 4 clipped bases, one of them an error
    # a PAF line: read aligned from its base 4 to the end, target from 8 on, NM 0
    paf = "\t".join(["r", str(len(read)), "4", str(len(read)), "+", "copy1|t", str(len(target)), "8", str(8 + len(core)),
                     str(len(core)), str(len(core)), "60", "NM:i:0"])
    assert E.clip_distance(paf) == 4
    assert E.semiglobal_distance(paf, read, target) == 1


def test_a_column_is_a_sun_of_the_one_copy_that_differs_from_all_others():
    copies = {0: "AAC-", 1: "AGCT", 2: "AGTT", 3: "AGTA"}
    # col 0: no difference; col 1: copy 0 private (SUN of 0); col 2: two vs two (shared); col 3: '-' is no allele, A private to 3
    assert E.column_classes(copies) == [None, ("sun", 0), ("shared", None), ("sun", 3)]


def test_a_read_s_support_for_its_copy_is_split_into_sun_and_shared_columns():
    copies = {0: "AACT", 1: "AGCT", 2: "AGTA", 3: "AGTT"}
    cls = E.column_classes(copies)            # [None, sun of 0, shared (CC vs TT), sun of 2]
    assert cls == [None, ("sun", 0), ("shared", None), ("sun", 2)]
    # copy 1's read: rules out copy 0 at col 1 and copy 2 at col 3 (another copy's SUN), shares col 2 with copy 0
    assert E.support("?GCT", copies, 1, cls) == dict(sun=0, shared=1, against_sun=2, conflict=0)
    # copy 0's read: its own SUN at col 1, a shared column, and a conflict at col 3
    assert E.support("?ACA", copies, 0, cls) == dict(sun=1, shared=1, against_sun=0, conflict=1)


def test_two_copies_differ_by_pair_columns():
    assert E.column_classes({0: "AC", 1: "AG"}) == [None, ("pair", None)]


def test_a_read_star_proof_is_split_into_own_sun_shared_and_conflict_columns():
    # candidates 0,1,2,3 (positions in the proof), assigned candidate 1
    cols = "10:T:ATTT,20:C:CCGG,30:A:GAGG,40:T:.TC.,50:G:AAAA"
    s = E.star_support(cols, 1)
    # col 10: copy 0 private (A), read T rules it out -> against_sun; col 20: CC vs GG shared, read C = copy 1 -> shared;
    # col 30: copy 1 private (A), read A -> own sun; col 40: only two candidates have bases -> pair (read T = copy 1);
    # col 50: read G matches no candidate (all A) -> conflict
    assert s == dict(sun=1, shared=1, pair=1, against_sun=1, conflict=1, cols=5)


def test_gap_compressed_distance_counts_each_indel_once():
    core = "ACGTTGCAGGCTAGCTAGGATCCGATCGATTACG"
    target = "GGGG" + core + "CCCC"
    # aligned core with a 5-bp insertion and 1 mismatch: NM 6 (5 inserted bases + 1 mismatch), gap-compressed 2
    paf = "\t".join(["r", "39", "0", "39", "+", "copy1|t", str(len(target)), "4", str(4 + len(core)), "33", "39", "60",
                     "NM:i:6", "cg:Z:20M5I14M"])
    assert E.gap_compressed_distance(paf, "A" * 39, target) == 2
