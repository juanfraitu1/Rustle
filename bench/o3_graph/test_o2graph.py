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


def test_a_path_call_the_genome_contradicts_is_abstained_when_no_confident_winner_takes_it():
    # --veto-abstain (2026-10-17, round 15's FAM90A errors): the called copy 45 aligns 16 edits away on the genome, four
    # siblings 5 away and tied among themselves -- the veto finds no single confident winner and kept the wrong path call
    g = {45: 16, 20: 5, 29: 5, 41: 5, 42: 5}
    assert O.veto_decide("45", g, 1445, primary=(20, 0), ldist=g) is None                          # today: kept
    assert O.veto_decide("45", g, 1445, primary=(20, 0), ldist=g, contradict_abstain=True) == ("-", "genome_contradicts")
    # a confident unique winner is still a veto, a genome tie or a thin gap still keeps the path call
    g2 = {45: 16, 20: 5, 29: 9}
    assert O.veto_decide("45", g2, 1445, primary=(20, 60), ldist=g2, contradict_abstain=True) == ("20", "veto_locus")
    assert O.veto_decide("45", {45: 5, 20: 5}, 1445, primary=(20, 0), ldist={45: 5, 20: 5}, contradict_abstain=True) is None
    assert O.veto_decide("45", {45: 6, 20: 5}, 1445, primary=(20, 0), ldist={45: 6, 20: 5}, contradict_abstain=True) is None
    # no genome alignment on the called copy (minimap2's secondary cap in a 57-copy family): no evidence either way
    assert O.veto_decide("45", {20: 5, 29: 5}, 1445, primary=(20, 0), ldist={20: 5, 29: 5}, contradict_abstain=True) is None


def test_the_contradiction_must_hold_on_the_called_copy_s_aligned_bases():
    # A copy's genome cost lies between its locus_distance (clips under NOVEL_INS free) and its veto_distance (clips paid in
    # full). The call is contradicted only when a rival's worst case beats the call's best case by the margin.
    # Set 8 public DDX11L.6 (correct): its alignment soft-clips a 15-bp terminal exon the aligner missed -- 16 with the clip,
    # 1 on its aligned bases -- and a paralog aligns the whole read at 3 (development: the first rule abstained it)
    v, l = {5: 3, 6: 16, 0: 5}, {5: 3, 6: 1, 0: 5}
    assert O.veto_decide("6", v, 330, primary=(5, 1), ldist=l, contradict_abstain=True) is None
    # r10 public DAZ.2 (correct): a 4-bp end clip on its own copy, the rivals at 1 either way
    assert O.veto_decide("2", {-1: 1, 0: 1, 1: 1, 2: 5}, 805, primary=(-1, 0), ldist={-1: 1, 0: 1, 1: 1, 2: 1},
                         contradict_abstain=True) is None
    # set 8 public OUTSIDE AQP7.3 (an uncatalogued copy's read, path-called to copy 2): the rival beats it by the margin
    # with and without clips
    v = {-1: 1, 0: 31, 1: 3, 2: 2, 3: 0}
    assert O.veto_decide("2", v, 336, primary=(3, 5), ldist={-1: 1, 0: 14, 1: 3, 2: 2, 3: 0},
                         contradict_abstain=True) == ("-", "genome_contradicts")
    # the call's only alignment leaves >= NOVEL_INS bases unaligned (no locus_distance): no evidence on its aligned bases
    assert O.veto_decide("45", {45: 40, 20: 5}, 1445, primary=(20, 0), ldist={20: 5}, contradict_abstain=True) is None


def test_a_locus_call_stands_only_when_its_clipped_bases_cannot_hand_the_margin_to_a_rival():
    # --locus-clipped (not adopted), the locus fallback: a call made on the genome alone must win by the margin in its worst
    # case (its clipped bases paid) against every rival's best case. Round 15 A119b 502 CROSS s0000493 (FAM90A.51's read): copy 45's
    # alignment clips the read's last 11 bases -- 3 edits on the rest, 14 counting them -- and nine siblings align all of it
    # at 5; locus_distance alone gave copy 45 the margin
    ld = {45: 3, 46: 5, 47: 5, 51: 5, 54: 5}
    vd = {45: 14, 46: 5, 47: 5, 51: 5, 54: 5}
    assert O.locus_decide(ld, 361) == (45, "locus")                                                # today
    assert O.locus_decide(ld, 361, vdist=vd) == (None, "locus_clipped")
    assert O.locus_decide({3: 1, 5: 4}, 1000, vdist={3: 1, 5: 4}) == (3, "locus")                  # no clip: as today
    assert O.locus_decide({3: 1, 5: 6}, 1000, vdist={3: 3, 5: 6}) == (3, "locus")                  # the margin survives it
    assert O.locus_decide({3: 1, O.ELSEWHERE: 4}, 1000, vdist={3: 3, O.ELSEWHERE: 4}) == (None, "locus_clipped")
    # a twin folded into the call (O1 --twins) is no rival, as in decide()
    assert O.locus_decide({3: 1, 7: 1, 5: 6}, 1000, twins={7: 3}, vdist={3: 2, 7: 1, 5: 6}) == (3, "locus")


def test_the_veto_mode_abstains_a_clipped_locus_call_with_locus_clipped(tmp_path):
    import csv
    import json
    w = tmp_path / "W"
    (w / "graph").mkdir(parents=True)
    copies = [dict(idx=0, contig="chrA", start=1000, end=5000, origin="curated"),
              dict(idx=1, contig="chrB", start=1000, end=5000, origin="curated")]
    json.dump({"F": dict(copies=copies, dropped=[], family_reads=["r1"])}, open(w / "graph" / "seeded.json", "w"))
    with open(w / "graph" / "o2.integ.m2.tsv", "w") as fh:
        fh.write("read\tfamily\tcall\treason\nr1\tF\t-\tunexplained\n")
    # copy 0 clips the read's last 11 bases (3 edits on the rest); copy 1 aligns all of it with 5
    _bam(str(w / "E0.bam"), [("r1", 0, "chrA", 2000, "100=1X100=1X100=1X86=11S"), ("r1", 256, "chrB", 2000, "100=5X295=")])
    assert O.run_locus_fallback(str(w), mode="veto") == 0
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integv.m2.tsv"), delimiter="\t")}
    assert out == {"r1": ("0", "locus")}                                                          # today
    O._VETO_ABSTAIN = True                                                                       # the path rule alone: as today
    try:
        assert O.run_locus_fallback(str(w), mode="veto") == 0
    finally:
        O._VETO_ABSTAIN = False
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integv.m2.tsv"), delimiter="\t")}
    assert out == {"r1": ("0", "locus")}
    O._LOCUS_CLIPPED = True
    try:
        assert O.run_locus_fallback(str(w), mode="veto") == 0
    finally:
        O._LOCUS_CLIPPED = False
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integv.m2.tsv"), delimiter="\t")}
    assert out == {"r1": ("-", "locus_clipped")}


# --o2-significance (2026-10-17, the O2 margin): r12 A119b's DUX4 errors were path calls by 2 gap-compressed edits that were
# two 1-bp indels, against 63 rivals; the v2 profiles' homopolymer indels are ~20x likelier than any one substitution
_CORE = "ACGTAAAACGTCAGTCGATCGGCTAGCTAGCATCGATCGACTAGC"   # AAAA at 4-7 (a homopolymer run), GG at 20-21 (a run of 2)
_WTS = {"X": 8.0, "hp": 4.0, "indel": 7.0}


def _wpaf(read, t, nm, cg, qs=0, ts=0):
    tl = len(_CORE)
    return "\t".join(["r", str(len(read)), str(qs), str(len(read)), "+", t, str(tl), str(ts), str(tl), "40", str(tl), "60",
                      f"NM:i:{nm}", f"cg:Z:{cg}"])


def test_the_weighted_distance_prices_each_edit_by_its_type():
    tx = {"copy0|t": _CORE}
    n = len(_CORE)
    hp_del = _CORE[:4] + _CORE[5:]                         # one A of AAAA: a homopolymer deletion
    assert O.gap_weighted(_wpaf(hp_del, "copy0|t", 1, f"4M1D{n - 5}M"), hp_del, _CORE, _WTS) == 4.0
    hp_del_r = _CORE[:7] + _CORE[8:]                       # the same deletion placed at the run's other end
    assert O.gap_weighted(_wpaf(hp_del_r, "copy0|t", 1, f"7M1D{n - 8}M"), hp_del_r, _CORE, _WTS) == 4.0
    del2 = _CORE[:20] + _CORE[21:]                         # one G of GG: a run of 2 is no homopolymer
    assert O.gap_weighted(_wpaf(del2, "copy0|t", 1, f"20M1D{n - 21}M"), del2, _CORE, _WTS) == 7.0
    hp_ins = _CORE[:4] + "A" + _CORE[4:]                   # AAAA -> AAAAA
    assert O.gap_weighted(_wpaf(hp_ins, "copy0|t", 1, f"4M1I{n - 4}M"), hp_ins, _CORE, _WTS) == 4.0
    ins = _CORE[:12] + "G" + _CORE[12:]                    # a G between C and A
    assert O.gap_weighted(_wpaf(ins, "copy0|t", 1, f"12M1I{n - 12}M"), ins, _CORE, _WTS) == 7.0
    sub = _CORE[:30] + "T" + _CORE[31:]                    # C -> T
    assert O.gap_weighted(_wpaf(sub, "copy0|t", 1, f"{n}M"), sub, _CORE, _WTS) == 8.0
    two = _CORE[:4] + _CORE[5:30] + "T" + _CORE[31:]       # a homopolymer deletion and a substitution
    assert O.gap_weighted(_wpaf(two, "copy0|t", 2, f"4M1D{n - 5}M"), two, _CORE, _WTS) == 12.0


def test_the_weighted_distance_prices_clipped_end_edits_as_substitutions():
    tx_seq = "GGGG" + "ACGA" + _CORE
    read = "TCGA" + _CORE                                  # 4 clipped bases, 1 of them differs from the path's flank
    ln = "\t".join(["r", str(len(read)), "4", str(len(read)), "+", "copy0|t", str(len(tx_seq)), "8", str(len(tx_seq)),
                    str(len(_CORE)), str(len(_CORE)), "60", "NM:i:0", f"cg:Z:{len(_CORE)}M"])
    assert O.gap_weighted(ln, read, tx_seq, _WTS) == 8.0
    assert O.gap_compressed(ln, read, tx_seq) == 1


def test_the_pruned_weighted_distance_equals_scoring_every_line():
    import random
    rng = random.Random(5)
    core = "".join(rng.choice("ACGT") for _ in range(300))
    tx = {f"copy{i}|t{j}": "".join(rng.choice("ACGT") for _ in range(40)) + core + "".join(rng.choice("ACGT") for _ in range(40))
          for i in range(3) for j in range(3)}
    reads = {f"r{k}": "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 30))) + core for k in range(4)}
    lines = []
    for q, rd in reads.items():
        for t in tx:
            qs, nd = len(rd) - len(core), rng.randint(0, 3)
            nm = nd + rng.randint(0, 6)
            lines.append("\t".join([q, str(len(rd)), str(qs), str(len(rd)), "+", t, str(len(tx[t])), "40", str(340 + nd), "300",
                                    "300", "60", f"NM:i:{nm}", f"cg:Z:150M{nd}D150M" if nd else "cg:Z:300M"]))
    rng.shuffle(lines)
    naive = {}
    for ln in lines:
        f = ln.split("\t")
        c = f[5].split("|")[0]
        naive.setdefault(f[0], {})
        naive[f[0]][c] = min(O.gap_weighted(ln, reads[f[0]], tx[f[5]], _WTS), naive[f[0]].get(c, float("inf")))
    assert O._paf_dist_w("\n".join(lines), reads, tx, _WTS) == naive


def test_edit_weights_are_minus_log_rates_of_each_edit_type_on_uniquely_placed_reads(tmp_path):
    import math
    import pysam
    w = tmp_path / "W"
    w.mkdir()
    ref = "ACGT" * 50 + "AAAAA" + "CGTA" * 50 + "TTT" + "GCAT" * 50     # runs >= 3: AAAAA, TTT (the CGTA/ACGT AA, TT joins aside)
    (w / "genome.ref.fa").write_text(f">chrA\n{ref}\n")
    pysam.faidx(str(w / "genome.ref.fa"))
    hdr = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrA", "LN": len(ref)}]}
    recs = [  # (name, flag, mapq, pos, cigar, sequence)
        ("a", 0, 60, 0, "100=", ref[:100]),
        ("b", 0, 60, 150, "51=1D150=", ref[150:201] + ref[202:352]),     # deletion of ref[201] -- inside AAAAA (200-204)
        ("c", 0, 60, 0, "50=1X49=", ref[:50] + ("A" if ref[50] != "A" else "C") + ref[51:100]),
        ("d", 0, 60, 0, "30=1I70=", ref[:30] + "G" + ref[30:100]),      # a G before ref[30] = G: a run of 1, no homopolymer
        ("e", 0, 10, 0, "50=1X49=", ref[:50] + ("A" if ref[50] != "A" else "C") + ref[51:100]),  # MAPQ 10: not counted
        ("f", 256, 60, 0, "50=1X49=", ""),                                  # a secondary: not counted
    ]
    with pysam.AlignmentFile(str(w / "u.bam"), "wb", header=hdr) as out:
        for name, flag, mq, pos, cg, seq in recs:
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.flag, a.reference_name, a.reference_start, a.cigarstring = name, flag, "chrA", pos, cg
            a.mapping_quality = mq
            if seq:
                a.query_sequence = seq
            out.write(a)
    pysam.sort("-o", str(w / "E0.bam"), str(w / "u.bam"))
    pysam.index(str(w / "E0.bam"))
    c = O.edit_counts(str(w))
    assert c["aligned"] == 100 + 201 + 100 + 100 and c["X"] == 1
    assert c["hp"] == 1 and c["indel"] == 1                # b's deletion inside AAAAA; d's insertion of G at a C|G boundary
    assert c["hp_runs"] == 1                                # AAAAA, inside b's aligned blocks only (TTT is past every read)
    wts = O.weights_from_counts(c)
    assert math.isclose(wts["X"], math.log(3 * (c["aligned"] + 2) / (c["X"] + 1)))
    assert math.isclose(wts["hp"], math.log((c["hp_runs"] + 2) / (c["hp"] + 1)))
    assert math.isclose(wts["indel"], math.log((c["aligned"] + 2) / (c["indel"] + 1)))


def test_a_path_call_must_clear_alpha_over_its_rivals_number_in_weighted_edits():
    import math
    # r12 A119b DUX4.9's read on copy 19: two homopolymer indels closer than 63 siblings -- 8 < ln(63 / 1e-3) = 11.05
    wd = {"19": 12.0, **{str(i): 20.0 for i in range(64) if i != 19}}
    pd = {"19": 2, **{str(i): 4 for i in range(64) if i != 19}}
    assert not O.significant("19", pd, wd)
    assert O.significant("19", pd, {**wd, **{str(i): 28.0 for i in range(64) if i != 19}})   # two substitutions: 16
    assert O.significant("19", {"19": 2, "3": 4, "4": 4}, {"19": 12.0, "3": 20.0, "4": 20.0})  # 8 >= ln(2 / 1e-3) = 7.6
    assert math.isclose(math.log(63 / O.ALPHA), 11.05, abs_tol=0.01)
    assert O.significant("19", {"19": 2}, {"19": 12.0})                       # a sole candidate
    assert O.significant("19", pd, {})                                        # no weighted distances: no gate
    # a twin folded into the call (decide's rule: within the margin of it) is no rival
    assert O.significant("3", {"3": 1, "7": 1, "5": 6}, {"3": 8.0, "7": 8.0, "5": 48.0}, twins={7: 3})
    assert not O.significant("3", {"3": 1, "7": 1, "5": 6}, {"3": 8.0, "7": 8.0, "5": 48.0})


def test_the_veto_mode_abstains_an_insignificant_path_call_with_o2_significance(tmp_path):
    import csv
    import json
    w = tmp_path / "W"
    d = w / "graph" / "o2paths.union.pad.integ"
    d.mkdir(parents=True)
    copies = [dict(idx=i, contig="chrA", start=1000 + 10000 * i, end=5000 + 10000 * i, origin="curated") for i in range(3)]
    json.dump({"F": dict(copies=copies, dropped=[], family_reads=["r1", "r2"])}, open(w / "graph" / "seeded.json", "w"))
    with open(w / "graph" / "o2.integ.m2.tsv", "w") as fh:
        fh.write("read\tfamily\tcall\treason\nr1\tF\t0\tassigned\nr2\tF\t0\tassigned\n")
    json.dump({"r1": {"0": 0, "1": 2, "2": 2}, "r2": {"0": 0, "1": 2, "2": 2}}, open(d / "F.dist.json", "w"))
    json.dump({"r1": {"0": 0.0, "1": 6.0, "2": 6.0}, "r2": {"0": 0.0, "1": 16.0, "2": 16.0}}, open(d / "F.wdist.json", "w"))
    _bam(str(w / "E0.bam"), [("r1", 0, "chrA", 2000, "1000="), ("r2", 0, "chrA", 2000, "1000=")])
    assert O.run_locus_fallback(str(w), mode="veto") == 0
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integv.m2.tsv"), delimiter="\t")}
    assert out == {"r1": ("0", "assigned"), "r2": ("0", "assigned")}                             # today
    O._SIGNIFICANCE = True
    try:
        assert O.run_locus_fallback(str(w), mode="veto") == 0
    finally:
        O._SIGNIFICANCE = False
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integv.m2.tsv"), delimiter="\t")}
    assert out == {"r1": ("-", "insignificant"), "r2": ("0", "assigned")}      # 6 < ln(2 / 1e-3) = 7.6 <= 16


def test_a_locus_call_pays_only_for_read_bases_a_rival_aligns_and_its_own_alignment_clips():
    # --locus-pairwise (2026-10-17): against each rival a copy pays its alignment's edits plus the read bases the rival's
    # alignment explains and its own leaves out -- bases every copy clips (a short exon past a junction, a read end) are
    # nobody's evidence. hits: {copy: (edits, read start, read end of its alignment)}
    # round 15 A119b 502 CROSS s0000493 (FAM90A.51's read): copy 45 aligns read[0:350] with 3 edits, its siblings all 361
    # bases with 5 -- copy 45 pays 3 + 11 = 14 against 5
    ld = {45: 3, 46: 5, 51: 5}
    assert O.locus_decide(ld, 361) == (45, "locus")
    assert O.locus_decide(ld, 361, hits={45: (3, 0, 350), 46: (5, 0, 361), 51: (5, 0, 361)}) == (None, "locus_clipped")
    # r10 public HNRNPCL.0 (correct): copy 0 aligns read[0:179] exactly, copies 1 and 2 read[15:194] with 2 and 3 edits
    hits = {0: (0, 0, 179), 1: (2, 15, 194), 2: (3, 15, 194)}
    assert O.locus_decide({0: 0, 1: 2, 2: 3}, 194, hits=hits) == (0, "locus")
    # r10 public OUTSIDE FAM236.2 (correct): copy 1 aligns read[2:149] exactly, its folded twin 3 read[9:156], copies 0
    # and 2 read[9:153] and read[5:149] with 9 edits
    hits = {1: (0, 2, 149), 3: (0, 9, 156), 0: (9, 9, 153), 2: (9, 5, 149)}
    assert O.locus_decide({1: 0, 3: 0, 0: 9, 2: 9}, 158, twins={3: 1}, hits=hits) == (1, "locus")


def test_genome_hits_keep_each_copy_s_best_alignment_with_its_read_interval(tmp_path):
    import json
    w = tmp_path / "W"
    (w / "graph").mkdir(parents=True)
    spans = {"F": [("chrA", 1000, 5000, 0), ("chrB", 1000, 5000, 1)]}
    _bam(str(w / "E0.bam"), [("r1", 0, "chrA", 2000, "100=1X100=1X100=1X86=11S"), ("r1", 256, "chrB", 2000, "100=5X295="),
                             ("r1", 256, "chrB", 3000, "50S350="),                 # a worse alignment of copy 1: 50 bases out
                             ("r2", 16, "chrA", 2000, "10S390=")])                 # reverse strand: the clip is the read's end
    hits = O.genome_hits(str(w), {"r1": {"F"}, "r2": {"F"}}, spans)
    assert hits[("r1", "F")] == {0: (3, 0, 389), 1: (5, 0, 400)}
    assert hits[("r2", "F")] == {0: (0, 0, 390)}


def test_the_veto_mode_abstains_a_locus_call_its_rival_explains_better_pairwise_with_locus_pairwise(tmp_path):
    import csv
    import json
    w = tmp_path / "W"
    (w / "graph").mkdir(parents=True)
    copies = [dict(idx=0, contig="chrA", start=1000, end=5000, origin="curated"),
              dict(idx=1, contig="chrB", start=1000, end=5000, origin="curated")]
    json.dump({"F": dict(copies=copies, dropped=[], family_reads=["r1"])}, open(w / "graph" / "seeded.json", "w"))
    with open(w / "graph" / "o2.integ.m2.tsv", "w") as fh:
        fh.write("read\tfamily\tcall\treason\nr1\tF\t-\tunexplained\n")
    _bam(str(w / "E0.bam"), [("r1", 0, "chrA", 2000, "100=1X100=1X100=1X86=11S"), ("r1", 256, "chrB", 2000, "100=5X295=")])
    O._LOCUS_PAIRWISE = True
    try:
        assert O.run_locus_fallback(str(w), mode="veto") == 0
    finally:
        O._LOCUS_PAIRWISE = False
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integv.m2.tsv"), delimiter="\t")}
    assert out == {"r1": ("-", "locus_clipped")}


def test_a_path_margin_over_a_rival_whose_own_locus_explains_the_read_better_than_its_paths_is_structure_only():
    # --structure-check (2026-10-17): r12 A119b KRTAP2.0's read retains an intron near its start; copy 1's annotated long
    # first exon absorbs it (path 2), copy 0's paths clip it (85), and the genome holds the read identically at both (0)
    pd = {1: 2, 0: 85, 3: 85, 2: 86}
    g = {0: 0, 1: 0, 2: 2, 3: 2}
    assert O.structure_only("1", pd, g, g)
    # rivals whose paths fit the read as their loci do: the paths decide (a sequence-made margin)
    assert not O.structure_only("0", {0: 0, 1: 2}, {0: 0, 1: 2}, {0: 0, 1: 2})
    # a rival whose paths miss the read's structure but whose locus is 3 edits worse: the genome separates them
    assert not O.structure_only("0", {0: 0, 1: 85}, {0: 0, 1: 3}, {0: 0, 1: 3})
    # r15 public 503 CROSS FAM90A.12's 201-bp read: only copy 11's paths align it, eight loci hold it exactly
    gg = {i: 0 for i in range(5, 13)}
    assert O.structure_only("11", {11: 0}, gg, gg)
    # r10 public DAZ.2 (correct): its rivals' alignments run 40-58 kb past their loci (none contained): no evidence
    assert not O.structure_only("2", {2: 1}, {2: 5}, {2: 1})
    # a twin folded into the call is no rival
    assert not O.structure_only("3", {3: 0, 7: 9}, {3: 0, 7: 0}, {3: 0, 7: 0}, twins={7: 3})
    assert O.structure_only("3", {3: 0, 7: 9}, {3: 0, 7: 0}, {3: 0, 7: 0})


def test_genome_distances_can_keep_only_alignments_inside_each_copy_s_padded_locus(tmp_path):
    w = tmp_path / "W"
    w.mkdir()
    spans = {"F": [("chrA", 1000, 5000, 0), ("chrB", 1000, 5000, 1)]}
    _bam(str(w / "E0.bam"), [("r1", 0, "chrA", 2000, "1000="), ("r1", 256, "chrB", 4500, "400=3X597="),   # chrB: runs past
                             ("r1", 256, "chrB", 1000 - O.PATH_PAD, "500=1X499=")])                          # inside the pad
    d, _ = O.genome_distances(str(w), {"r1": {"F"}}, spans)
    assert d[("r1", "F")] == {0: 0, 1: 1}
    d, _ = O.genome_distances(str(w), {"r1": {"F"}}, spans, contained=O.PATH_PAD)
    assert d[("r1", "F")] == {0: 0, 1: 1}
    _bam(str(w / "E0.bam"), [("r1", 0, "chrA", 2000, "1000="), ("r1", 256, "chrB", 4500, "400=3X597=")])
    d, _ = O.genome_distances(str(w), {"r1": {"F"}}, spans, contained=O.PATH_PAD)
    assert d[("r1", "F")] == {0: 0}


def test_the_veto_mode_abstains_a_structure_only_path_call_with_structure_check(tmp_path):
    import csv
    import json
    w = tmp_path / "W"
    d = w / "graph" / "o2paths.union.pad.integ"
    d.mkdir(parents=True)
    copies = [dict(idx=0, contig="chrA", start=1000, end=5000, origin="curated"),
              dict(idx=1, contig="chrB", start=1000, end=5000, origin="curated")]
    json.dump({"F": dict(copies=copies, dropped=[], family_reads=["r1", "r2"])}, open(w / "graph" / "seeded.json", "w"))
    with open(w / "graph" / "o2.integ.m2.tsv", "w") as fh:
        fh.write("read\tfamily\tcall\treason\nr1\tF\t1\tassigned\nr2\tF\t1\tassigned\n")
    json.dump({"r1": {"1": 2, "0": 85}, "r2": {"1": 2, "0": 4}}, open(d / "F.dist.json", "w"))
    _bam(str(w / "E0.bam"), [("r1", 0, "chrA", 2000, "1000="), ("r1", 256, "chrB", 2000, "1000="),         # loci identical
                             ("r2", 0, "chrB", 2000, "1000="), ("r2", 256, "chrA", 2000, "500=4X496=")])    # paths = loci
    O._STRUCTURE_CHECK = True
    try:
        assert O.run_locus_fallback(str(w), mode="veto") == 0
    finally:
        O._STRUCTURE_CHECK = False
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integv.m2.tsv"), delimiter="\t")}
    assert out == {"r1": ("-", "structure_only"), "r2": ("1", "assigned")}


def test_the_poisson_quantile_is_the_smallest_count_whose_cdf_reaches_p():
    import math
    for lam in (0.0, 0.8, 3.2, 9.6, 40.0):
        q = O.poisson_quantile(lam, 1 - O.ALPHA)
        cdf = lambda n: sum(math.exp(-lam) * lam ** i / math.factorial(i) for i in range(n + 1))
        assert cdf(q) >= 1 - O.ALPHA and (q == 0 or cdf(q - 1) < 1 - O.ALPHA)


def test_a_short_read_s_few_sequencing_edits_are_within_reach_of_its_own_copy():
    # --length-gate (2026-10-17, round 16's B1 shortfall): the fixed gate (delta + floor) x length allowed a 281-bp read 3.6
    # edits; SLC35G.1's read carried 4 and was left unexplained, where the current method placed it. Within reach: allelic
    # divergence up to delta (p99 between haplotypes) plus sequencing errors up to Poisson(floor x length)'s 1 - alpha quantile
    floor = 0.0032
    assert 4 > O.gate(floor) * 281                                                  # today: unexplained
    assert O.within_reach(4, 281, floor) and O.within_reach(7, 281, floor)
    assert not O.within_reach(12, 281, floor)
    assert O.within_reach(40, 3000, floor) and not O.within_reach(60, 3000, floor)   # long reads: close to delta + floor


def test_the_genome_step_rescues_path_calls_the_fixed_gate_and_margin_left_out():
    floor = 0.0032
    # an unexplained sole-candidate read within reach (round 16 A119b 601 SLC35G.1, 281 bp, 4 edits)
    assert O.rescue("unexplained", {"1": 4}, {"1": 22.4}, 281, floor, length_gate=True) == (1, "length_gate")
    assert O.rescue("unexplained", {"1": 4}, {"1": 22.4}, 281, floor) is None                    # without the flag
    assert O.rescue("unexplained", {"1": 12}, {"1": 90.0}, 281, floor, length_gate=True) is None   # beyond reach
    # a margin-1 read whose one difference is a substitution, against 2 rivals (round 16 A119b 602 SPAG.1): 7.7 >= ln(2/1e-3)
    pd, wd = {"1": 3, "0": 4, "2": 5}, {"1": 14.8, "0": 22.5, "2": 30.2}
    assert O.rescue("thin", pd, wd, 789, floor, significance_margin=True) == (1, "significant")
    assert O.rescue("thin", pd, wd, 789, floor) is None
    # ... but not against 4 rivals (CLEC.2: 7.6 < ln(4/1e-3) = 8.3), nor a difference of one homopolymer indel
    assert O.rescue("thin", {"2": 9, "3": 10, "4": 12, "5": 13, "6": 13}, {"2": 52.0, "3": 59.6, "4": 74.9, "5": 80.0,
                                                                          "6": 80.0}, 3452, floor, significance_margin=True) is None
    assert O.rescue("tied", {"1": 3, "0": 3}, {"1": 18.9, "0": 22.8}, 700, floor, significance_margin=True) is None
    # an unexplained read the length gate leaves thin goes on to the significance margin
    assert O.rescue("unexplained", {"3": 4, "2": 5}, {"3": 22.8, "2": 30.5}, 312, floor, length_gate=True,
                    significance_margin=True) == (3, "significant")


def test_the_veto_mode_places_rescued_path_calls_through_the_same_checks(tmp_path):
    import csv
    import json
    w = tmp_path / "W"
    d = w / "graph" / "o2paths.union.pad.integ"
    d.mkdir(parents=True)
    copies = [dict(idx=i, contig="chrA", start=1000 + 10000 * i, end=5000 + 10000 * i, origin="curated") for i in range(3)]
    json.dump({"F": dict(copies=copies, dropped=[], family_reads=["r1", "r2", "r3"])}, open(w / "graph" / "seeded.json", "w"))
    with open(w / "graph" / "o2.integ.m2.tsv", "w") as fh:
        fh.write("read\tfamily\tcall\treason\nr1\tF\t-\tunexplained\nr2\tF\t-\tthin\nr3\tF\t-\tthin\n")
    json.dump({"r1": {"1": 4}, "r2": {"1": 3, "0": 4, "2": 5}, "r3": {"1": 3, "0": 4, "2": 5}}, open(d / "F.dist.json", "w"))
    json.dump({"r1": {"1": 22.4}, "r2": {"1": 14.8, "0": 22.5, "2": 30.2}, "r3": {"1": 14.8, "0": 22.5, "2": 30.2}},
              open(d / "F.wdist.json", "w"))
    (w / "reads.E0.fq").write_text("".join(f"@{n}\n{'A' * k}\n+\n{'I' * k}\n" for n, k in (("r1", 281), ("r2", 789), ("r3", 789))))
    _bam(str(w / "E0.bam"), [("r1", 0, "chrA", 12000, "281="),                          # the genome agrees with copy 1
                             ("r2", 0, "chrA", 12000, "789="),
                             ("r3", 0, "chrA", 2000, "789="), ("r3", 256, "chrA", 12000, "400=9X380=")])  # copy 0, confidently
    O._LENGTH_GATE = O._SIGNIFICANCE_MARGIN = True
    try:
        assert O.run_locus_fallback(str(w), mode="veto", floor=0.0032) == 0
    finally:
        O._LENGTH_GATE = O._SIGNIFICANCE_MARGIN = False
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integv.m2.tsv"), delimiter="\t")}
    assert out == {"r1": ("1", "length_gate"), "r2": ("1", "significant"), "r3": ("0", "veto_locus")}


def test_reads_too_short_for_map_hifi_s_score_floor_are_aligned_to_the_paths_again(tmp_path):
    # --short-path-align (2026-10-17, round 16): map-hifi's minimum DP score (-s 200) is beyond any read under ~200 bp, so
    # round 16's 150-200-bp reads had no path alignment at all (OR4A.5, GGTLC.1); with -s 40 they align to their own copy's
    # isoform. The genome step realigns only the reads the path step left without any alignment.
    import csv
    import json
    import random
    rng = random.Random(7)
    exon = "".join(rng.choice("ACGT") for _ in range(400))
    other = exon[:200] + "".join(rng.choice("ACGT") for _ in range(200))       # copy 1 differs over the second half
    w = tmp_path / "W"
    d = w / "graph" / "o2paths.union.pad.integ"
    d.mkdir(parents=True)
    (d / "F.tx.fa").write_text(f">copy0|t0\n{exon}\n>copy1|t1\n{other}\n")
    copies = [dict(idx=0, contig="chrA", start=1000, end=5000, origin="curated"),
              dict(idx=1, contig="chrB", start=1000, end=5000, origin="curated")]
    json.dump({"F": dict(copies=copies, dropped=[], family_reads=["r1"])}, open(w / "graph" / "seeded.json", "w"))
    with open(w / "graph" / "o2.integ.m2.tsv", "w") as fh:
        fh.write("read\tfamily\tcall\treason\nr1\tF\t-\tunexplained\n")
    json.dump({}, open(d / "F.dist.json", "w"))                                   # no path alignment
    read = exon[220:400]                                                          # 180 bp from copy 0's own half
    (w / "reads.E0.fq").write_text(f"@r1\n{read}\n+\n{'I' * len(read)}\n")
    _bam(str(w / "E0.bam"), [("r1", 0, "chrA", 2000, "180=")])
    assert O.run_locus_fallback(str(w), mode="veto") == 0
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integv.m2.tsv"), delimiter="\t")}
    assert out == {"r1": ("0", "locus")}                                          # today: the genome fallback
    O._SHORT_PATH_ALIGN = True
    try:
        assert O.run_locus_fallback(str(w), mode="veto") == 0
    finally:
        O._SHORT_PATH_ALIGN = False
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integv.m2.tsv"), delimiter="\t")}
    assert out == {"r1": ("0", "short_path")}


def _struct_world(tmp_path):
    import json
    import random
    rng = random.Random(11)
    e1, e2, e3 = ("".join(rng.choice("ACGT") for _ in range(n)) for n in (150, 76, 150))
    e2b = e2[:20] + ("A" if e2[20] != "A" else "C") + e2[21:40] + ("G" if e2[40] != "G" else "T") + e2[41:60] + (
        "C" if e2[60] != "C" else "A") + e2[61:]                                  # copy 1's exon 2: three PSVs
    w = tmp_path / "W"
    d = w / "graph" / "o2paths.union.pad.integ"
    d.mkdir(parents=True)
    (d / "F.tx.fa").write_text(f">copy0|full\n{e1 + e2 + e3}\n>copy0|skip\n{e1 + e3}\n>copy1|full\n{e1 + e2b + e3}\n"
                               f">copy1|skip\n{e1 + e3}\n")
    copies = [dict(idx=0, contig="chrA", start=1000, end=5000, origin="curated"),
              dict(idx=1, contig="chrB", start=1000, end=5000, origin="curated")]
    json.dump({"F": dict(copies=copies, dropped=[], family_reads=["r1"])}, open(w / "graph" / "seeded.json", "w"))
    read = e1[:70] + "T" + e1[70:] + e2 + e3                                      # copy 0's full isoform, one inserted base
    (w / "reads.E0.fq").write_text(f"@r1\n{read}\n+\n{'I' * len(read)}\n")
    return w, d, read


def test_structure_aware_distances_compare_copies_on_the_isoform_that_fits_the_read(tmp_path):
    # --struct-rescue (2026-10-17, round 16's thin calls): gap compression counts a whole missing exon as one insertion, so a
    # read's best path can be an isoform that lacks one of its exons -- and the paralogous bases in that exon go uncompared.
    # SPAG.1's reads: 3 vs 4 on copy 1's exon-skipping isoform (a 76-bp insertion), 4 vs 9 on their own isoform.
    w, d, read = _struct_world(tmp_path)
    gc = O.struct_distances(str(d), "F", {"r1": read}, struct=False)
    gs = O.struct_distances(str(d), "F", {"r1": read})
    assert gc["r1"]["1"] - gc["r1"]["0"] == 1                                    # thin: copy 1's skip isoform, one event
    assert gs["r1"] == {"0": 1, "1": 4}                                           # each copy on the read's own isoform


def test_the_veto_mode_rescues_a_thin_path_call_on_structure_aware_distances(tmp_path):
    import csv
    import json
    w, d, read = _struct_world(tmp_path)
    with open(w / "graph" / "o2.integ.m2.tsv", "w") as fh:
        fh.write("read\tfamily\tcall\treason\nr1\tF\t-\tthin\n")
    json.dump({"r1": {"0": 1, "1": 2}}, open(d / "F.dist.json", "w"))
    _bam(str(w / "E0.bam"), [("r1", 0, "chrA", 2000, f"{len(read)}=")])
    O._STRUCT_RESCUE = True
    try:
        assert O.run_locus_fallback(str(w), mode="veto") == 0
    finally:
        O._STRUCT_RESCUE = False
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integv.m2.tsv"), delimiter="\t")}
    assert out == {"r1": ("0", "struct_aware")}


def test_the_twin_guard_abstains_a_call_an_uncatalogued_locus_explains_as_well_but_not_a_dropped_silent_copy(tmp_path):
    # --twin-guard (2026-10-18, round 17: every cell failed OUTSIDE -- a pair's uncatalogued twin, unseeded by O1, had its
    # reads placed on the one remaining copy; LGALS7 223, CXorf49 197 errors). The guard of register 1315 abstains when a locus
    # off every seeded copy explains the read as well; it was rejected because a catalogued silent copy O1 dropped is such a
    # locus too (round 10: array units). A dropped catalogued copy's locus is no longer "off the catalogue".
    import csv
    import json
    w = tmp_path / "W"
    (w / "graph").mkdir(parents=True)
    json.dump({"F": dict(copies=[dict(idx=0, contig="chrA", start=1000, end=5000, origin="curated")], dropped=[1],
                         family_reads=["r1", "r2"])}, open(w / "graph" / "seeded.json", "w"))
    json.dump({"F": [dict(copy="F.0", idx=0, contig="chrA", start=1000, end=5000), dict(copy="F.1", idx=1, contig="chrB",
                                                                                       start=1000, end=5000)]},
              open(w / "isoforms.json", "w"))
    with open(w / "graph" / "o2.integ.m2.tsv", "w") as fh:
        fh.write("read\tfamily\tcall\treason\nr1\tF\t0\tassigned\nr2\tF\t0\tassigned\n")
    _bam(str(w / "E0.bam"), [("r1", 0, "chrA", 2000, "1000="), ("r1", 256, "chrB", 30000, "1000="),   # an uncatalogued twin
                             ("r2", 0, "chrA", 2000, "1000="), ("r2", 256, "chrB", 2000, "1000=")])   # the dropped silent copy
    assert O.run_locus_fallback(str(w), mode="veto") == 0
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integv.m2.tsv"), delimiter="\t")}
    assert out == {"r1": ("0", "assigned"), "r2": ("0", "assigned")}                              # today
    O._TWIN_GUARD = True
    try:
        assert O.run_locus_fallback(str(w), mode="veto") == 0
    finally:
        O._TWIN_GUARD = False
    out = {r["read"]: (r["call"], r["reason"]) for r in csv.DictReader(open(w / "graph" / "o2.integv.m2.tsv"), delimiter="\t")}
    assert out == {"r1": ("-", "uncatalogued_tie"), "r2": ("0", "assigned")}


def test_o2_acts_only_on_as_tied_reads_minimap2_resolves_the_rest_and_the_unexplained_are_o3_s():
    S, E = O.as_scope, O.ELSEWHERE
    assert S({}, "1", "assigned", 1000) == ("-", "as_unaligned")
    assert S({E: 0}, "-", "unexplained", 1000) == ("-", "as_elsewhere")                # best placement off the family's copies
    assert S({0: 0, E: 0}, "0", "assigned", 1000) == ("-", "as_elsewhere")             # tied with an unseeded locus
    assert S({1: 3}, "0", "assigned", 1000) == ("1", "as_resolved")                     # one best copy: minimap2's, whatever O2 said
    assert S({1: 20}, "1", "assigned", 1000) == ("-", "as_unexplained")                 # beyond the gate: nothing explains it
    assert S({1: None}, "1", "assigned", 1000) == ("-", "as_unexplained")               # the locus lacks the read's sequence
    assert S({1: 12}, "-", "thin", 1000, floor=0.005) == ("1", "as_resolved")           # the gate adds the error floor
    assert S({0: 2, 1: 2}, "1", "assigned", 1000) == ("1", "assigned")                  # an AS tie: O2's call among the tied copies
    assert S({0: 2, 1: 2}, "2", "veto_locus", 1000) == ("-", "as_tie_outside")          # ... never off them
    assert S({0: 2, 1: 2}, "-", "tied", 1000) == ("-", "tied")                          # ... or O2's abstention


def test_the_best_as_loci_of_a_read_are_its_best_scoring_alignments_seeded_copies_or_elsewhere(tmp_path):
    import pysam
    w = tmp_path / "W"
    w.mkdir()
    hdr = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": "c1", "LN": 20000}]})
    with pysam.AlignmentFile(str(w / "E0.bam"), "wb", header=hdr) as out:
        for name, pos, flag, score, cig in [("r1", 1100, 0, 300, "148=1X1="), ("r1", 5100, 256, 300, "150="),   # tie 0 / 1
                                            ("r1", 9100, 256, 250, "150="),
                                            ("r2", 1100, 0, 300, "150="), ("r2", 5100, 256, 290, "150="),        # resolved on 0
                                            ("r3", 15000, 0, 300, "150="), ("r3", 5100, 256, 300, "150=")]:      # tie with off
            a = pysam.AlignedSegment(hdr)
            a.query_name, a.reference_id, a.reference_start, a.flag, a.mapping_quality = name, 0, pos, flag, 0
            a.cigarstring, a.query_sequence = cig, "A" * 150
            a.set_tag("AS", score)
            out.write(a)
    spans = {"F": [("c1", 1000, 2000, 0), ("c1", 5000, 6000, 1), ("c1", 9000, 10000, 2)]}
    got, qlen = O.best_as_loci(str(w), {"r1": {"F"}, "r2": {"F"}, "r3": {"F"}}, spans)
    assert got[("r1", "F")] == {0: 1, 1: 0} and got[("r2", "F")] == {0: 0} and got[("r3", "F")] == {1: 0, O.ELSEWHERE: 0}
    assert qlen["r1"] == 150
