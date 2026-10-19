import importlib.util
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_samplecopy", os.path.join(HERE, "samplecopy.py"))
C = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C)


def _reads(hap, n, cols, err, rng, cover=1.0):
    """n reads of a haplotype {col: base}, each column observed with probability `cover`, miscalled at rate err."""
    out = []
    for _ in range(n):
        r = {}
        for c in cols:
            if rng.random() < cover:
                b = hap[c]
                r[c] = rng.choice([x for x in "ACGT" if x != b]) if rng.random() < err else b
        out.append(r)
    return out


def test_a_column_is_informative_when_two_alleles_each_have_enough_reads():
    reads = [{1: "A", 2: "C", 3: "G"}] * 6 + [{1: "T", 2: "C", 3: "G"}] * 4 + [{1: "A", 2: "C", 3: "T"}]
    assert C.informative_columns(reads, min_reads=3, min_frac=0.1) == [1]          # column 3's T: one read, an error
    assert C.informative_columns(reads, min_reads=1, min_frac=0.05) == [1, 3]


def test_the_mixture_separates_two_haplotypes_that_differ_at_a_few_columns_despite_errors_and_gaps():
    rng = random.Random(1)
    cols = list(range(10))
    h1 = {c: "A" for c in cols}
    h2 = {**h1, 2: "C", 5: "G", 8: "T"}
    reads = _reads(h1, 40, cols, 0.01, rng, cover=0.7) + _reads(h2, 25, cols, 0.01, rng, cover=0.7)
    fit = C.fit_clusters(reads, cols, kmax=4, seed=0)
    assert len(fit["haplotypes"]) == 2                                              # 4 started, merged back to 2
    groups = [max(range(len(p)), key=p.__getitem__) for p in fit["posteriors"]]
    seen = [i for i, r in enumerate(reads) if any(c in r for c in (2, 5, 8))]    # a read with no distinguishing column
    g1 = {groups[i] for i in seen if i < 40}                                       # is undecidable: no claim on it
    g2 = {groups[i] for i in seen if i >= 40}
    assert len(g1) == 1 and len(g2) == 1 and g1 != g2
    assert len(seen) < len(reads)                                                  # the data does hold such reads


def test_clusters_closer_than_two_columns_merge_and_tiny_ones_dissolve():
    rng = random.Random(2)
    cols = list(range(8))
    h1 = {c: "A" for c in cols}
    reads = _reads(h1, 30, cols, 0.0, rng) + _reads({**h1, 4: "C"}, 20, cols, 0.0, rng)   # one column apart
    assert len(C.fit_clusters(reads, cols, kmax=3, seed=0)["haplotypes"]) == 1
    reads = _reads(h1, 30, cols, 0.0, rng) + _reads({**h1, 1: "G", 6: "T"}, 2, cols, 0.0, rng)   # two reads only
    assert len(C.fit_clusters(reads, cols, kmax=3, seed=0, min_size=3)["haplotypes"]) == 1


def test_a_cluster_takes_the_copy_its_untied_reads_resolve_to_by_a_clear_majority():
    member = [0, 0, 0, 0, 1, 1, 1, 2]
    resolved = ["F.0", "F.0", "F.0", None, "F.1", "F.2", None, "F.2"]                 # None: a tied read (no resolution)
    assert C.label_clusters(member, resolved, n_clusters=3, min_reads=2, min_share=0.8) == {0: "F.0", 1: None, 2: None}
    assert C.label_clusters(member, resolved, n_clusters=3, min_reads=1, min_share=0.5) == {0: "F.0", 1: "F.1", 2: "F.2"}


def test_a_tied_read_takes_its_cluster_s_copy_only_when_sure_and_the_copy_is_one_it_ties_on():
    labels = {0: "F.0", 1: "F.1", 2: None}
    assert C.assign_tied([0.995, 0.005, 0.0], labels, {"F.0", "F.1"}) == "F.0"
    assert C.assign_tied([0.9, 0.1, 0.0], labels, {"F.0", "F.1"}) is None             # not sure
    assert C.assign_tied([0.0, 0.0, 1.0], labels, {"F.0", "F.1"}) is None             # an unlabelled cluster
    assert C.assign_tied([0.0, 1.0, 0.0], labels, {"F.0", "F.2"}) is None             # not a copy the read ties on


def test_a_read_s_pileup_is_its_base_at_every_aligned_reference_position_indels_left_out():
    import pysam
    hdr = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": "a", "LN": 1000}]})
    a = pysam.AlignedSegment(hdr)
    a.query_name, a.reference_id, a.reference_start, a.flag = "r", 0, 100, 0
    a.cigarstring, a.query_sequence = "2S3=1X2=1I2=1D2=", "GGACGTACTTGCA"    # GG clipped; ACG, T (X), AC, T inserted, TG, -, CA
    p = C.read_pileup(a)
    assert p == {100: "A", 101: "C", 102: "G", 103: "T", 104: "A", 105: "C", 106: "T", 107: "G", 109: "C", 110: "A"}


def test_the_anchor_is_the_copy_most_untied_reads_resolve_to():
    assert C.choose_anchor({"F.2": 9, "F.0": 3}, ["F.0", "F.1", "F.2"]) == "F.2"
    assert C.choose_anchor({}, ["F.0", "F.1"]) == "F.0"


def test_a_read_s_evidence_ignores_the_cluster_sizes_so_a_read_with_no_informative_site_stays_undecided():
    haps = [{1: "A", 2: "C"}, {1: "T", 2: "C"}]
    ev = C.evidence_posteriors([{1: "A"}, {2: "C"}, {}], haps, 0.01)
    assert ev[0][0] > 0.99                                   # it carries the site that separates them
    assert ev[1] == [0.5, 0.5] and ev[2] == [0.5, 0.5]       # a shared or missing site: no evidence, whatever the weights
