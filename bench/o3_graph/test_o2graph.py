import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


O = _load("o2graph")
PATHS = {"copy0": {"1", "2", "3", "5"}, "copy1": {"1", "2", "4", "6"}, "copy2": {"1", "7", "4", "6"}}


def test_a_read_goes_to_the_copy_whose_path_holds_its_allele_nodes():
    # informative nodes: on some copies' paths but not all (node 1 is shared by all three)
    assert O.assign(["1", "2", "3", "5"], PATHS, margin=2) == ("copy0", {"copy0": 0, "copy1": 2, "copy2": 3})


def test_a_tie_or_a_thin_margin_abstains():
    assert O.assign(["1", "4", "6"], PATHS, margin=2)[0] is None           # copy1 and copy2 both hold every allele node
    assert O.assign(["1", "2", "4"], PATHS, margin=2)[0] is None           # copy1 misses 0, copy0 / copy2 miss 1: margin 1
    assert O.assign(["1"], PATHS, margin=2)[0] is None                     # no informative node at all


def test_paths_assign_the_copy_with_fewest_edits_by_the_margin():
    assert O.decide({"copy0": 1, "copy1": 4, "copy2": 6}, qlen=1000, margin=2) == ("copy0", "assigned")
    assert O.decide({"copy0": 1, "copy1": 2}, qlen=1000, margin=2) == (None, "thin")
    assert O.decide({"copy0": 3, "copy1": 3}, qlen=1000, margin=2) == (None, "tied")
    assert O.decide({"copy0": 1}, qlen=1000, margin=2) == ("copy0", "assigned")       # the only copy it aligns to


def test_a_read_no_copy_explains_within_delta_abstains():
    # 30 edits over 1,000 bases is 3% away from the best copy: an absent copy's read, O3's business
    assert O.decide({"copy0": 30, "copy1": 45}, qlen=1000, margin=2) == (None, "unexplained")
    assert O.decide({}, qlen=1000, margin=2) == (None, "unexplained")


def _world(tmp_path, cur_rows, o2_rows):
    (tmp_path / "cur").mkdir()
    (tmp_path / "graph").mkdir()
    (tmp_path / "P.fam.copies.tsv").write_text("family_id\tcopy_idx\ttid\nA\t0\tA.0\nA\t1\tA.1\nB\t0\tB.0\nB\t1\tB.1\n")
    (tmp_path / "reads.truth.E0.tsv").write_text("read\tcopy\tkind\tfamily\tisoform\tnovel\n"
                                                 "r1\tA.0\tpresent\tA\tt\t0\nr2\tB.0\tpresent\tB\tt\t0\nr3\tA.X\toutside\tA\tt\t0\n")
    (tmp_path / "cur" / "A.assignments.tsv").write_text(
        "read_name\tfamily_id\tstatus\tcatalog_copy_idx\n" + "".join("\t".join(r) + "\n" for r in cur_rows))
    (tmp_path / "graph" / "o2.t.tsv").write_text("read\tfamily\tcall\treason\n" + "".join("\t".join(r) + "\n" for r in o2_rows))
    return str(tmp_path)


def test_an_assignment_in_another_family_is_scored_wrong(tmp_path):
    # r1 (A.0): the current method names copy 0 in A (right) and copy 1 in B (a wrong-family call); r2 (B.0): path O2 names
    # copy 1 of A; r3 (an uncatalogued copy of A) assigned in B by the current method: wrong as well
    w = _world(tmp_path, [("r1", "A", "assigned", "0"), ("r1", "B", "assigned", "1"), ("r2", "B", "ambiguous", "-"),
                          ("r3", "B", "assigned", "0")],
               [("r1", "A", "0", "assigned"), ("r2", "A", "1", "assigned"), ("r2", "B", "-", "tied")])
    res = O.score_world(w, "o2.t.tsv")
    assert res["A"]["cur_correct"] == 1 and res["A"]["graph_correct"] == 1
    assert res["B"]["cur_wrong_family"] == 2 and res["A"].get("cur_wrong_family", 0) == 0
    assert res["A"]["graph_wrong_family"] == 1 and res["B"].get("graph_wrong_family", 0) == 0


def test_a_gain_counts_as_sequence_supported_only_when_the_genome_separates_the_copies():
    assert O.seq_supported({"copy0": 3, "copy1": 5}, "copy0")             # another copy's locus is 2 edits worse
    assert not O.seq_supported({"copy0": 3, "copy1": 4}, "copy0")         # within 1 edit: only the annotation separates them
    assert O.seq_supported({"copy0": 3}, "copy0")                         # no other locus explains the read at all
    assert not O.seq_supported({"copy1": 3}, "copy0")                     # its own locus does not explain it


def test_the_score_reports_sequence_supported_correct_calls(tmp_path):
    w = _world(tmp_path, [("r1", "A", "assigned", "0")], [("r1", "A", "0", "assigned")])
    res = O.score_world(w, "o2.t.tsv", locus={"A": {"r1": {"copy0": 2, "copy1": 3}}})
    assert res["A"]["cur_correct"] == 1 and res["A"].get("cur_correct_seq", 0) == 0
    res = O.score_world(w, "o2.t.tsv", locus={"A": {"r1": {"copy0": 2, "copy1": 9}}, "B": {"r1": {"copy0": 0}}})
    assert res["A"]["cur_correct_seq"] == 1 and res["A"]["graph_correct_seq"] == 1


def test_reads_of_a_paralog_annotated_structure_are_scored_on_their_own(tmp_path):
    w = _world(tmp_path, [("r1", "A", "assigned", "0"), ("r4", "A", "assigned", "0")],
               [("r1", "A", "0", "assigned"), ("r4", "A", "1", "assigned")])
    with open(tmp_path / "reads.truth.E0.tsv", "a") as fh:
        fh.write("r4\tA.0\tpresent\tA\tproj:A.1|t9\t0\n")              # copy A.0 expressing a structure annotated on A.1
    res = O.score_world(w, "o2.t.tsv")["A"]
    assert res["proj_reads"] == 1 and res["cur_correct_proj"] == 1
    assert res["graph_assigned_proj"] == 1 and res.get("graph_correct_proj", 0) == 0     # path O2 names the annotated paralog
    assert res["reads"] == 2 and res["graph_correct"] == 1


def test_union_paths_give_every_copy_each_paralog_structure_it_can_carry():
    copies = [dict(idx=0, isoforms={"a": [[0, 10], [20, 30]], "b": [[0, 10], [25, 30]]}),
              dict(idx=1, isoforms={"c": [[100, 110], [120, 130]]})]

    def projector(c, queries):            # copy 1 can carry a's structure; b does not project onto it; c projects onto 0
        out = {}
        for name in queries:
            if c["idx"] == 1:
                out[name] = [(100, 110), (120, 130)] if name.endswith("|a") else None
            else:
                out[name] = [(0, 10), (40, 50)]
        return out
    u = O.union_structures(copies, projector)
    assert sorted(u[0]) == ["a", "b", "proj:copy1|c"]          # c's intron (10, 40) is new to copy 0
    assert sorted(u[1]) == ["c"]                                # a projects onto copy 1 as copy 1's own chain: nothing new


def test_locus_o2_decides_on_genomic_locus_distances(tmp_path):
    import json
    w = _world(tmp_path, [("r1", "A", "assigned", "0"), ("r2", "B", "assigned", "0")], [])
    json.dump({"A": {"r1": {"copy0": 2, "copy1": 9}}, "B": {"r2": {"copy0": 3, "copy1": 4}}},
              open(tmp_path / "graph" / "locus.dist.json", "w"))
    (tmp_path / "reads.E0.fq").write_text("@r1\n" + "A" * 1000 + "\n+\n" + "I" * 1000 + "\n@r2\n" + "A" * 1000 + "\n+\n" + "I" * 1000 + "\n")
    O.run_locus(w, margin=2)
    rows = {(r.split("\t")[0], r.split("\t")[1]): r.split("\t")[2:] for r in (tmp_path / "graph" / "o2.locus.m2.tsv").read_text().splitlines()[1:]}
    assert rows[("r1", "A")] == ["0", "assigned"] and rows[("r2", "B")] == ["-", "thin"]


def _paf(q, qlen, qs, qe, strand, t, tlen, ts, te, nm):
    return "\t".join(map(str, [q, qlen, qs, qe, strand, t, tlen, ts, te, 0, 0, 60, f"NM:i:{nm}", "tp:A:P"]))


def test_read_bases_past_a_path_end_are_annotation_not_distance():
    paf = "\n".join([
        _paf("r", 1000, 0, 984, "+", "copy0|t", 984, 0, 984, 2),       # 16 bases past copy 0's path end: not counted
        _paf("r", 1000, 0, 1000, "+", "copy1|t", 1100, 50, 1050, 5),
        _paf("s", 1000, 0, 900, "+", "copy0|t", 1100, 0, 900, 1),       # stops 200 bp before the path end: a real difference
        _paf("u", 1000, 16, 1000, "-", "copy0|t", 984, 0, 984, 2)])     # minus strand: the read's start lies past the path end
    assert O._paf_dist(paf) == {"r": {"copy0": 18, "copy1": 5}, "s": {"copy0": 101}, "u": {"copy0": 18}}
    assert O._paf_dist(paf, ends=True) == {"r": {"copy0": 2, "copy1": 5}, "s": {"copy0": 101}, "u": {"copy0": 2}}


def test_padded_paths_extend_both_transcript_ends_with_the_copy_s_genome():
    assert O.padded([(100, 200), (300, 400)], 50) == [(50, 200), (300, 450)]
    assert O.padded([(30, 200)], 50) == [(0, 250)]                       # not before the contig start
    assert O.padded([(100, 200), (300, 400)], 0) == [(100, 200), (300, 400)]


def test_a_gain_is_psv_separated_when_the_runner_up_carries_the_read_s_structure():
    true = [(0, 100), (200, 300), (400, 500)]
    assert O.gain_kind(true, [[(0, 100), (205, 300), (400, 500)]]) == "psv"          # same chain within 10 bp
    assert O.gain_kind(true, [[(0, 100), (400, 500)], [(0, 100), (250, 300), (400, 500)]]) == "structure"
    assert O.gain_kind(true, []) == "structure"


def test_a_read_union_cannot_explain_falls_back_to_the_loci_unless_another_family_explains_it():
    loc = {"copy0": 2, "copy1": 9}
    assert O.fallback_call("1", "assigned", loc, 1000, other_explained=False) == ("1", "assigned")      # union decides
    assert O.fallback_call("-", "tied", loc, 1000, other_explained=False) == ("-", "tied")
    assert O.fallback_call("-", "unexplained", loc, 1000, other_explained=False) == ("0", "locus")     # structure-free
    assert O.fallback_call("-", "unexplained", loc, 1000, other_explained=True) == ("-", "other_family")
    assert O.fallback_call("-", "unexplained", {"copy0": 40, "copy1": 45}, 1000, other_explained=False) == ("-", "unexplained")
    assert O.fallback_call("-", "unexplained", {}, 1000, other_explained=False) == ("-", "unexplained")


def test_reads_of_a_cluster_with_linked_private_alleles_are_flagged_whatever_its_divergence():
    clusters = [dict(family="F", host="0", cluster="0", n_reads="25", status="unexplained", est_div="0.0012"),
                dict(family="F", host="0", cluster="1", n_reads="40", status="explained", est_div="0.0"),
                dict(family="F", host="1", cluster="0", n_reads="12", status="unexplained", est_div="0.02"),
                dict(family="G", host="0", cluster="0", n_reads="30", status="unresolved", est_div="0.001")]
    members = {("F", "0", "0"): {"a", "b"}, ("F", "0", "1"): {"c"}, ("F", "1", "0"): {"d"}, ("G", "0", "0"): {"e"}}
    assert O.flag_reads(clusters, members) == {("F", "a"), ("F", "b")}
    assert O.flag_reads(clusters, members, min_reads=10) == {("F", "a"), ("F", "b"), ("F", "d")}


def test_a_flagged_read_is_reported_as_a_candidate_not_assigned():
    rows = [dict(read="a", family="F", call="1", reason="assigned"), dict(read="b", family="F", call="-", reason="tied"),
            dict(read="c", family="F", call="0", reason="locus")]
    out = O.apply_flag(rows, {("F", "a"), ("F", "b")})
    assert [(r["call"], r["reason"]) for r in out] == [("-", "flagged"), ("-", "tied"), ("0", "locus")]


def test_a_new_seeded_copy_maps_to_the_truth_copy_whose_locus_it_overlaps():
    truth_copies = [dict(copy="F.0", chrom="c1", start="1000", end="3000"), dict(copy="F.5", chrom="c2", start="5000", end="6500"),
                    dict(copy="F.6", chrom="c2", start="9000", end="9500")]
    new = [dict(idx=7, contig="c2", start=5100, end=6300), dict(idx=8, contig="c3", start=0, end=900)]
    assert O.map_new_copies(new, truth_copies) == {7: "F.5", 8: None}


def test_calls_are_scored_against_each_read_s_true_copy_whatever_catalog_named_it():
    truth = {"r1": dict(copy="F.0", family="F", kind="present"), "r2": dict(copy="F.5", family="F", kind="outside"),
             "r3": dict(copy="F.9", family="F", kind="absent"), "r4": dict(copy="F.0", family="F", kind="present")}
    idx2copy = {("F", 0): "F.0", ("F", 7): "F.5", ("F", 8): None}
    calls = {("r1", "F"): 0, ("r2", "F"): 7, ("r3", "F"): 0, ("r4", "F"): 8}
    s = O.score_assign(calls, idx2copy, truth)
    assert s == dict(correct=2, wrong=2, outside_own=1)


def test_placements_of_an_undetectable_absent_copy_are_counted_apart_from_errors():
    truth = {"r1": dict(copy="F.0", family="F", kind="present"), "r2": dict(copy="F.9", family="F", kind="absent"),
             "r3": dict(copy="F.8", family="F", kind="absent")}
    idx2copy = {("F", 0): "F.0", ("F", 1): "F.1"}
    calls = {("r1", "F"): 0, ("r2", "F"): 1, ("r3", "F"): 1}
    # F.9 is within delta of the genome (undetectable): its reads cannot be told from that locus's own
    s = O.score_assign(calls, idx2copy, truth, unknowable={"F.9"})
    assert s == dict(correct=1, wrong=1, outside_own=0, unknowable_placed=1)


def test_gap_compressed_distance_counts_indel_events_and_aligns_clipped_ends():
    core = "ACGTTGCAGGCTAGCTAGGATCCGATCGATTACG"
    tx = {"copy1|t": "GGGG" + "ACGA" + core + "CCCC", "copy2|t": "GGGG" + "TTTT" + core + "CCCC"}
    read = "TCGA" + core                                    # 4 clipped bases; vs copy1's flank 1 edit, vs copy2's 3
    def paf(t, nm, cg):
        return "\t".join(["r", str(len(read)), "4", str(len(read)), "+", t, str(len(tx[t])), "8", str(8 + len(core)),
                          str(len(core)), str(len(core)), "60", f"NM:i:{nm}", f"cg:Z:{cg}"])
    lines = "\n".join([paf("copy1|t", 6, "10M5D24M"), paf("copy2|t", 1, "34M")])
    d = O._paf_dist_gc(lines, {"r": read}, tx)
    assert d == {"r": {"copy1": 1 + 1 + 1, "copy2": 1 + 3}}     # copy1: 1 mismatch + 1 deletion event + 1 end edit


def test_the_pruned_gap_compressed_distance_equals_scoring_every_line():
    import random
    rng = random.Random(4)
    core = "".join(rng.choice("ACGT") for _ in range(300))
    tx = {f"copy{i}|t{j}": "".join(rng.choice("ACGT") for _ in range(40)) + core + "".join(rng.choice("ACGT") for _ in range(40))
          for i in range(3) for j in range(3)}
    reads = {f"r{k}": "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 30))) + core for k in range(4)}
    lines = []
    for q, rd in reads.items():
        for t in tx:
            qs, nm = len(rd) - len(core), rng.randint(0, 8)
            lines.append("\t".join([q, str(len(rd)), str(qs), str(len(rd)), "+", t, str(len(tx[t])), "40", "340", "300", "300",
                                    "60", f"NM:i:{nm}", f"cg:Z:{150}M{rng.randint(0, 3)}D{150}M"]))
    rng.shuffle(lines)
    paf = "\n".join(lines)
    naive = {}
    for ln in lines:
        f = ln.split("\t")
        c = f[5].split("|")[0]
        naive.setdefault(f[0], {})
        naive[f[0]][c] = min(O.gap_compressed(ln, reads[f[0]], tx[f[5]]), naive[f[0]].get(c, 10 ** 9))
    assert O._paf_dist_gc(paf, reads, tx) == naive


def test_the_locus_fallback_distance_counts_indel_events_and_refuses_novel_sequence():
    # pysam op codes: 7 '=', 8 X, 1 I, 2 D, 3 N, 4 S
    assert O.locus_distance([(7, 300), (8, 1), (3, 2000), (7, 200), (2, 4), (7, 100), (1, 1), (7, 50)]) == 3
    assert O.locus_distance([(7, 300), (1, O.NOVEL_INS), (7, 300)]) is None        # a novel exon: O3's, not O2's
    assert O.locus_distance([(4, O.NOVEL_INS), (7, 600)]) is None                  # an unaligned end of that length too
    assert O.locus_distance([(4, 5), (7, 600)]) == 0                                # a short clip is not sequence


def test_the_locus_fallback_assigns_only_a_clear_winner_and_calls_a_genome_tie_tied():
    assert O.locus_decide({3: 1, 5: 4}, qlen=1000) == (3, "locus")
    assert O.locus_decide({3: 1, 5: 2}, qlen=1000) == (None, "locus_thin")
    assert O.locus_decide({3: 1, 5: 1}, qlen=1000) == (None, "locus_tied")
    assert O.locus_decide({3: 20}, qlen=1000) == (None, "unexplained")              # beyond delta: still a candidate for O3
    assert O.locus_decide({}, qlen=1000) == (None, "unexplained")


def test_the_locus_fallback_leaves_a_read_the_genome_places_off_the_seeded_copies_unexplained():
    # ELSEWHERE: the read's best alignment outside every seeded copy of the family (an unseeded or dropped copy, another family)
    assert O.locus_decide({3: 4, O.ELSEWHERE: 0}, qlen=1000) == (None, "unexplained")
    assert O.locus_decide({3: 1, O.ELSEWHERE: 2}, qlen=1000) == (None, "locus_thin")
    assert O.locus_decide({3: 0, O.ELSEWHERE: 5}, qlen=1000) == (3, "locus")


def _bam(path, reads):
    """reads: [(name, flag, contig, pos0, cigar)] on chrA/chrB, sorted and indexed."""
    import pysam
    hdr = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrA", "LN": 100000}, {"SN": "chrB", "LN": 100000}]}
    with pysam.AlignmentFile(path + ".u.bam", "wb", header=hdr) as out:
        for name, flag, contig, pos, cigar in reads:
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.flag, a.reference_name, a.reference_start, a.cigarstring = name, flag, contig, pos, cigar
            a.mapping_quality = 0 if flag & 256 else 60
            out.write(a)
    pysam.sort("-o", path, path + ".u.bam")
    pysam.index(path)


def test_the_locus_fallback_places_only_reads_integrated_o2_left_unexplained(tmp_path):
    import csv
    import json
    w = tmp_path / "W"
    (w / "graph").mkdir(parents=True)
    copies = [dict(idx=0, contig="chrA", start=1000, end=5000, origin="curated"),
              dict(idx=1, contig="chrB", start=1000, end=5000, origin="curated")]
    json.dump({"F": dict(copies=copies, dropped=[], family_reads=["r1", "r2", "r3", "r4", "r5"]),
               "G": dict(copies=[dict(idx=0, contig="chrA", start=50000, end=52000, origin="curated")], dropped=[],
                         family_reads=["r5"])}, open(w / "graph" / "seeded.json", "w"))
    rows = [("r1", "F", "0", "assigned"), ("r2", "F", "-", "unexplained"), ("r3", "F", "-", "unexplained"),
            ("r4", "F", "-", "unexplained"), ("r5", "F", "-", "unexplained"), ("r5", "G", "0", "assigned"),
            ("r6", "F", "-", "unexplained")]
    with open(w / "graph" / "o2.integ.m2.tsv", "w") as fh:
        fh.write("read\tfamily\tcall\treason\n" + "".join("\t".join(r) + "\n" for r in rows))
    _bam(str(w / "E0.bam"), [
        ("r2", 0, "chrB", 2000, "500=1X499="), ("r2", 256, "chrA", 2000, "500=4X496="),    # copy 1 by 3: placed
        ("r3", 0, "chrA", 2000, "500=1X499="), ("r3", 256, "chrB", 2000, "500=2X498="),    # margin 1: thin
        ("r4", 0, "chrA", 30000, "1000="), ("r4", 256, "chrA", 2000, "500=3X497="),        # best off the seeded copies
        ("r5", 0, "chrA", 2000, "1000="),                                                  # O2 explains it in family G
        ("r6", 0, "chrA", 2000, "400=60I540=")])                                           # a novel insertion
    assert O.run_locus_fallback(str(w)) == 0
    out = {(r["read"], r["family"]): (r["call"], r["reason"])
           for r in csv.DictReader(open(w / "graph" / "o2.integl.m2.tsv"), delimiter="\t")}
    assert out == {("r1", "F"): ("0", "assigned"), ("r2", "F"): ("1", "locus"), ("r3", "F"): ("-", "locus_thin"),
                   ("r4", "F"): ("-", "unexplained"), ("r5", "F"): ("-", "other_family"), ("r5", "G"): ("0", "assigned"),
                   ("r6", "F"): ("-", "unexplained")}
