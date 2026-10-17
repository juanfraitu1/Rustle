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


def test_the_combined_o2_lets_a_copy_explain_a_read_by_its_genome_alignment_where_no_path_carries_the_structure(tmp_path):
    import csv
    import json
    w = tmp_path / "W"
    d = w / "graph" / "o2paths.union.pad.integ"
    d.mkdir(parents=True)
    copies = [dict(idx=0, contig="chrA", start=1000, end=5000, origin="curated"),
              dict(idx=1, contig="chrB", start=1000, end=5000, origin="curated")]
    json.dump({"F": dict(copies=copies, dropped=[], family_reads=["r7", "r8", "r9", "r2"])}, open(w / "graph" / "seeded.json", "w"))
    rows = [("r7", "F", "0", "assigned"), ("r8", "F", "0", "assigned"), ("r9", "F", "0", "assigned"), ("r2", "F", "-", "unexplained")]
    with open(w / "graph" / "o2.integ.m2.tsv", "w") as fh:
        fh.write("read\tfamily\tcall\treason\n" + "".join("\t".join(r) + "\n" for r in rows))
    json.dump({"r7": {"0": 4, "1": 30}, "r8": {"0": 1, "1": 6}, "r9": {"0": 2, "1": 9}, "r2": {"0": 40}}, open(d / "F.dist.json", "w"))
    (w / "reads.E0.fq").write_text("".join(f"@{n}\n{'A' * 1000}\n+\n{'I' * 1000}\n" for n in ("r7", "r8", "r9", "r2")))
    _bam(str(w / "E0.bam"), [
        ("r7", 0, "chrB", 2000, "1000="), ("r7", 256, "chrA", 2000, "500=4X496="),     # copy 1 by its genome: no path of it fits
        ("r8", 0, "chrA", 2000, "500=1X499="), ("r8", 256, "chrB", 2000, "500=5X495="),  # paths and genome agree
        ("r9", 0, "chrA", 30000, "1000="), ("r9", 256, "chrA", 2000, "500=3X497="),      # the genome's best is off the copies
        ("r2", 0, "chrB", 2000, "500=1X499="), ("r2", 256, "chrA", 2000, "500=4X496=")])
    assert O.run_locus_fallback(str(w), mode="all") == 0
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integc.m2.tsv"), delimiter="\t")}
    assert out == {"r7": ("1", "genome_assigned"), "r8": ("0", "assigned"), "r9": ("-", "genome_unexplained"),
                   "r2": ("1", "genome_assigned")}
    # mode "seeded": the guard against the best alignment off the copies only for the reads the paths leave unexplained
    assert O.run_locus_fallback(str(w), mode="seeded") == 0
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integs.m2.tsv"), delimiter="\t")}
    assert out == {"r7": ("1", "genome_assigned"), "r8": ("0", "assigned"), "r9": ("0", "assigned"), "r2": ("1", "locus")}
    # mode "resolve": O2's own calls stand; only its unexplained (fallback) and tied / thin reads (genome) are decided here
    with open(w / "graph" / "o2.integ.m2.tsv", "a") as fh:
        fh.write("r10\tF\t-\ttied\nr11\tF\t-\tthin\n")
    d2 = json.load(open(d / "F.dist.json"))
    d2.update({"r10": {"0": 3, "1": 3}, "r11": {"0": 2, "1": 3}})
    json.dump(d2, open(d / "F.dist.json", "w"))
    (w / "reads.E0.fq").write_text("".join(f"@{n}\n{'A' * 1000}\n+\n{'I' * 1000}\n" for n in ("r7", "r8", "r9", "r2", "r10", "r11")))
    _bam(str(w / "E0.bam"), [
        ("r7", 0, "chrB", 2000, "1000="), ("r7", 256, "chrA", 2000, "500=4X496="),
        ("r8", 0, "chrA", 2000, "500=1X499="), ("r8", 256, "chrB", 2000, "500=5X495="),
        ("r9", 0, "chrA", 30000, "1000="), ("r9", 256, "chrA", 2000, "500=3X497="),
        ("r2", 0, "chrB", 2000, "500=1X499="), ("r2", 256, "chrA", 2000, "500=4X496="),
        ("r10", 0, "chrA", 2000, "1000="), ("r10", 256, "chrB", 2000, "500=3X497="),   # the genome breaks the paths' tie
        ("r11", 0, "chrA", 2000, "500=2X498="), ("r11", 256, "chrB", 2000, "500=3X497=")])  # still thin
    assert O.run_locus_fallback(str(w), mode="resolve") == 0
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integr.m2.tsv"), delimiter="\t")}
    assert out == {"r7": ("0", "assigned"), "r8": ("0", "assigned"), "r9": ("0", "assigned"), "r2": ("1", "locus"),
                   "r10": ("0", "genome_assigned"), "r11": ("-", "thin")}


def test_o2_paths_dump_each_read_s_distance_to_every_copy_for_the_combined_decision():
    assert O.path_dists({"r1": {"copy0": 3, "copy12": 7}}) == {"r1": {"0": 3, "12": 7}}


def test_exact_union_paths_keep_a_paralog_structure_a_few_bases_off_the_copy_s_own():
    copies = [dict(idx=0, isoforms={"a": [[0, 10], [20, 30]]}), dict(idx=1, isoforms={"c": [[100, 110], [120, 130]]})]

    def projector(c, queries):            # a lands on copy 1 three bases off c's junctions
        return {n: ([(100, 113), (123, 130)] if c["idx"] == 1 else None) for n in queries}
    assert sorted(O.union_structures(copies, projector)[1]) == ["c"]
    assert sorted(O.union_structures(copies, projector, tol=0)[1]) == ["c", "proj:copy0|a"]


def test_the_projector_offers_each_paralog_isoform_s_exon_lengths_only_for_the_micro_exon_rule():
    assert O.projection_spec([(0, 10), (20, 26)], micro=0) == 2
    assert O.projection_spec([(0, 10), (20, 26)], micro=30) == [10, 6]


def _end_cost_dp(p, t, anchored_at_end):
    if len(p) > O.EXT_CAP:
        return len(p)
    if anchored_at_end:
        p, t = p[::-1], t[::-1]
    prev = list(range(len(t) + 1))
    for i in range(1, len(p) + 1):
        cur = [i] + [0] * len(t)
        for j in range(1, len(t) + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (p[i - 1] != t[j - 1]))
        prev = cur
    return min(prev)


def test_the_bit_parallel_end_cost_equals_the_dynamic_programme():
    import random
    rng = random.Random(3)
    cases = [("", ""), ("", "ACGT"), ("ACGT", ""), ("A" * 51, "A" * 80), ("ACGTTG", "ACGTTGCA"), ("AAAA", "TTTT")]
    for _ in range(3000):
        m, n = rng.randint(0, 52), rng.randint(0, 85)
        base = "".join(rng.choice("ACGT") for _ in range(max(m, n)))
        p = "".join(c if rng.random() > 0.15 else rng.choice("ACGTN") for c in base[:m])
        cases.append((p, base[:n] if rng.random() < 0.7 else "".join(rng.choice("ACGT") for _ in range(n))))
    for p, t in cases:
        for end in (False, True):
            assert O.end_cost(p, t, end) == _end_cost_dp(p, t, end), (p, t, end)


def test_a_within_delta_twin_takes_a_read_from_its_curated_copy_only_by_the_margin():
    tw = {"copy3": "copy0"}                                      # O1's new copy 3 lies within delta of curated copy 0
    assert O.decide({"copy0": 2, "copy3": 1, "copy1": 9}, qlen=1000, twins=tw) == ("copy0", "assigned")   # thin: the catalog's
    assert O.decide({"copy0": 4, "copy3": 1, "copy1": 9}, qlen=1000, twins=tw) == ("copy3", "assigned")   # by the margin
    assert O.decide({"copy0": 1, "copy3": 1, "copy1": 2}, qlen=1000, twins=tw) == (None, "thin")          # 0 vs 1 still thin
    assert O.decide({"copy3": 1, "copy1": 9}, qlen=1000, twins=tw) == ("copy3", "assigned")               # no curated rival
    assert O.decide({"copy0": 2, "copy3": 1, "copy1": 9}, qlen=1000) == (None, "thin")                     # without twins


def test_the_unexplained_gate_widens_by_the_world_s_error_floor_when_given():
    assert O.decide({"copy0": 2, "copy1": 9}, qlen=184) == (None, "unexplained")                        # 2 > 0.00958 x 184
    assert O.decide({"copy0": 2, "copy1": 9}, qlen=184, max_div=O.DELTA + 0.0029) == ("copy0", "assigned")
    assert O.gate(0.0029) == O.DELTA + 0.0029 and O.gate(0.0) == O.DELTA


def test_structural_gap_compression_counts_a_long_indel_at_its_length():
    # round 12 review C1: a paralog's path absorbed a retained intron as ONE 142-bp insertion (one edit) while the read's own
    # copy's path clipped it at full length; gapstruct compresses only indels < NOVEL_INS bases (sequencing-type)
    f = ["r", "900", "0", "900", "+", "t", "900", "0", "758", "750", "900", "60", "NM:i:146", "cg:Z:274M142I365M1I116M"]
    assert O._gc_core(f) == 5                                    # 3 mismatches + 2 indel events
    assert O._gc_core(f, struct=True) == 3 + 142 + 1             # the 142-bp insertion at its length, the 1-bp one compressed
    g = ["r", "500", "0", "500", "+", "t", "500", "0", "514", "485", "500", "60", "NM:i:16", "cg:Z:306M15I141M"]
    assert O._gc_core(g) == O._gc_core(g, struct=True) == 1 + 1  # a 15-bp repeat-unit indel stays one event


def test_the_genome_vetoes_a_path_call_its_own_alignments_contradict_by_the_margin():
    # round 12 review C1: CDRT15.0 reads with a retained intron -- a paralog's path absorbs the intron as one gap-compressed
    # insertion while the own paths clip it; the genome puts them on CDRT15.0 (MAPQ 53-60, NM 0-4)
    assert O.veto_decide("1", {0: 0, 1: 5}, qlen=400, primary=(0, 60)) == ("0", "veto_locus")
    assert O.veto_decide("1", {0: 0}, qlen=400, primary=(0, 53)) == ("0", "veto_locus")   # no alignment at the called copy
    assert O.veto_decide("1", {0: 0, 1: 0}, qlen=400, primary=(0, 60)) is None   # DNA-identical: the paths' structure decides
    assert O.veto_decide("1", {0: 0, 1: 1}, qlen=400, primary=(0, 60)) is None   # under the margin
    assert O.veto_decide("1", {0: 20}, qlen=400, primary=(0, 60)) is None        # the genome explains it nowhere better
    assert O.veto_decide("1", {O.ELSEWHERE: 0, 1: 6}, qlen=400, primary=(O.ELSEWHERE, 60)) == ("-", "veto_unexplained")
    assert O.veto_decide("1", {0: 0, 2: 0, 1: 6}, qlen=400, primary=(0, 60)) is None     # a genome tie keeps the path call


def test_the_genome_vetoes_only_where_it_places_the_read_confidently():
    # development (set 8 public, DDX11L): the true copy's alignment soft-clips a 15-bp 5' exon the aligner missed, a paralog
    # aligns it (primary, MAPQ 1) -- the veto moved a correct call; every veto that fixed an error had MAPQ 53-60
    assert O.veto_decide("6", {5: 3, 6: 16}, qlen=330, primary=(5, 1)) is None
    assert O.veto_decide("6", {5: 3, 6: 16}, qlen=330, primary=(5, O.VETO_MAPQ)) == ("5", "veto_locus")
    assert O.veto_decide("6", {5: 3, 6: 16}, qlen=330, primary=(6, 60)) is None   # the primary is not on the genome's winner
    assert O.veto_decide("1", {O.ELSEWHERE: 0, 1: 6}, qlen=400, primary=(O.ELSEWHERE, 20)) is None


def test_the_veto_mode_rechecks_assigned_reads_and_places_unexplained_ones_as_the_fallback(tmp_path):
    import csv
    import json
    w = tmp_path / "W"
    (w / "graph").mkdir(parents=True)
    copies = [dict(idx=0, contig="chrA", start=1000, end=5000, origin="curated"),
              dict(idx=1, contig="chrB", start=1000, end=5000, origin="curated")]
    json.dump({"F": dict(copies=copies, dropped=[], family_reads=["r1", "r2", "r3"])},
              open(w / "graph" / "seeded.json", "w"))
    rows = [("r1", "F", "1", "assigned"), ("r2", "F", "0", "assigned"), ("r3", "F", "-", "unexplained")]
    with open(w / "graph" / "o2.integ.m2.tsv", "w") as fh:
        fh.write("read\tfamily\tcall\treason\n" + "".join("\t".join(r) + "\n" for r in rows))
    _bam(str(w / "E0.bam"), [
        ("r1", 0, "chrA", 2000, "1000="), ("r1", 256, "chrB", 2000, "500=5X495="),       # path says 1, genome says 0 by 5
        ("r2", 0, "chrA", 2000, "1000="), ("r2", 256, "chrB", 2000, "1000="),             # genome tie: the path call stands
        ("r3", 0, "chrB", 2000, "500=1X499="), ("r3", 256, "chrA", 2000, "500=4X496=")])  # the fallback, as integl
    assert O.run_locus_fallback(str(w), mode="veto") == 0
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integv.m2.tsv"), delimiter="\t")}
    assert out == {"r1": ("0", "veto_locus"), "r2": ("0", "assigned"), "r3": ("1", "locus")}


def test_the_veto_charges_clipped_and_long_inserted_bases_in_full():
    # round 12 development: KRTAP2.1 reads align to KRTAP2.0 as 27S46=105N5=1X82=1X53= -- 2 by locus_distance, a sub-NOVEL_INS clip
    # being free -- against 5 on their own locus, which aligns those bases; the veto then moved them to KRTAP2.0
    assert O.veto_distance([(4, 27), (7, 46), (3, 105), (7, 5), (8, 1), (7, 82), (8, 1), (7, 53)]) == 29
    assert O.veto_distance([(7, 300), (1, O.NOVEL_INS + 10), (7, 300)]) == O.NOVEL_INS + 10
    assert O.veto_distance([(7, 300), (1, 3), (8, 2), (2, 5), (7, 100), (5, 40)]) == 4      # hard clips are not the read's


def test_the_veto_keeps_a_path_call_whose_rival_alignment_only_clips_the_difference(tmp_path):
    import csv
    import json
    w = tmp_path / "W"
    (w / "graph").mkdir(parents=True)
    copies = [dict(idx=0, contig="chrA", start=1000, end=5000, origin="curated"),
              dict(idx=1, contig="chrB", start=1000, end=5000, origin="curated")]
    json.dump({"F": dict(copies=copies, dropped=[], family_reads=["r1"])}, open(w / "graph" / "seeded.json", "w"))
    with open(w / "graph" / "o2.integ.m2.tsv", "w") as fh:
        fh.write("read\tfamily\tcall\treason\nr1\tF\t1\tassigned\n")
    _bam(str(w / "E0.bam"), [("r1", 0, "chrA", 2000, "27S473="), ("r1", 256, "chrB", 2000, "24=3I1X3=1X3=1X465=")])
    assert O.run_locus_fallback(str(w), mode="veto") == 0
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integv.m2.tsv"), delimiter="\t")}
    assert out == {"r1": ("1", "assigned")}


def test_the_alignment_cache_returns_minimap2_s_own_lines_and_skips_reads_it_holds(tmp_path):
    # O2_ALN_CACHE (2026-10-17): keyed by the path set and minimap2's arguments, per read SEQUENCE (worlds of one family set
    # share their reads under other names); exact, because minimap2 aligns each read on its own
    import random
    import subprocess
    rng = random.Random(3)
    path = "".join(rng.choice("ACGT") for _ in range(3000))
    (tmp_path / "tx.fa").write_text(f">copy0|t1\n{path}\n>copy1|t2\n{path[:1500]}{'A' * 40}{path[1540:]}\n")
    reads = {"a": path[100:900], "b": path[1200:2100], "c": "".join(rng.choice("ACGT") for _ in range(600))}
    args = ["minimap2", "-c", "-x", "map-hifi", "-N", "200", "-p", "0.3", "--secondary=yes"]
    (tmp_path / "r.fa").write_text("".join(f">{n}\n{s}\n" for n, s in reads.items()))
    direct = subprocess.run(args + [str(tmp_path / "tx.fa"), str(tmp_path / "r.fa")], capture_output=True, text=True,
                            check=True).stdout
    cache = tmp_path / "cache"
    got = O.cached_minimap2(str(cache), args, str(tmp_path / "tx.fa"), reads, str(tmp_path / "w1"))
    # minimap2 breaks equal-score ties by a hash of the read's NAME, so a tied line's tp / s2 tags can differ; O2's distances
    # read NM, the CIGAR and the coordinates only, and those are the same
    tx = {"copy0|t1": path, "copy1|t2": f"{path[:1500]}{'A' * 40}{path[1540:]}"}
    strip = lambda paf: sorted("\t".join(x for x in ln.split("\t") if not x.startswith(("tp:", "s1:", "s2:")))  # noqa: E731
                               for ln in paf.splitlines())
    assert strip(got) == strip(direct)
    assert O._paf_dist(got) == O._paf_dist(direct)
    assert O._paf_dist_gc(got, reads, tx) == O._paf_dist_gc(direct, reads, tx)
    # another world: the same sequences under other names, one new read; only the new one is aligned
    reads2 = {"x": reads["a"], "y": reads["b"], "z": path[2000:2800]}
    got2 = O.cached_minimap2(str(cache), args, str(tmp_path / "tx.fa"), reads2, str(tmp_path / "w2"))
    assert O._last_aligned == 1
    assert sorted(ln.split("\t", 1)[1] for ln in got2.splitlines() if ln.split("\t")[0] in "xy") == sorted(
        ln.split("\t", 1)[1] for ln in got.splitlines() if ln.split("\t")[0] in "ab")
    assert any(ln.startswith("z\t") for ln in got2.splitlines())


def test_the_guard_abstains_when_an_unseeded_locus_explains_the_read_as_well():
    # --o2-genome guard (2026-10-17, the O2 proof's OUTSIDE worlds): O1 did not seed C4's / HSFY's uncatalogued near-identical
    # twin, so its reads had one candidate in the catalog; the genome still sees the twin as a locus off every seeded copy
    E = O.ELSEWHERE
    assert O.tie_guard("0", {0: 3, E: 3})                  # a tie with the unseeded locus
    assert O.tie_guard("0", {0: 4, 1: 9, E: 2})            # the unseeded locus is better
    assert not O.tie_guard("0", {0: 3, E: 4})              # the called copy is better
    assert not O.tie_guard("0", {0: 3, 1: 3})              # a seeded rival is path O2's business, not the guard's
    assert not O.tie_guard("0", {E: 3})                    # no alignment on the called copy: the veto's case


def test_twin_abstain_leaves_a_read_its_twins_cannot_separate():
    # --twin-abstain (2026-10-17, the O2 proof's OUTSIDE worlds: VCY's uncatalogued twin is real, and the annotation prior
    # gave its 94 reads to the curated copy): within the margin either way, the read is undecidable
    tw = {"copy1": "copy0"}
    assert O.decide({"copy0": 5, "copy1": 4}, 1000, twins=tw) == ("copy0", "assigned")              # the prior (default)
    assert O.decide({"copy0": 5, "copy1": 4}, 1000, twins=tw, twin_abstain=True) == (None, "tied_twin")
    assert O.decide({"copy0": 4, "copy1": 4}, 1000, twins=tw, twin_abstain=True) == (None, "tied_twin")
    assert O.decide({"copy0": 3, "copy1": 5}, 1000, twins=tw, twin_abstain=True) == ("copy0", "assigned")  # by the margin
    assert O.decide({"copy0": 7, "copy1": 4}, 1000, twins=tw, twin_abstain=True) == ("copy1", "assigned")
