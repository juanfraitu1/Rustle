"""Tests of donor_reads.py (PREREG r1797 / r1798): counterparts, expression scenarios, read drawing."""
import collections, json, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import donor_reads as D
from records import Rec

def row(qs, qe, ts, te, tname="chr1_MATERNAL", strand="+", mapq=60, cg=None):
    return dict(qlen=10000, qs=qs, qe=qe, strand=strand, tname=tname, ts=ts, te=te, mapq=mapq, tp="P", cg=cg or f"{qe - qs}M", de=0.001)

def copy(gid, start, exons_by_tid, chrom="chr1", strand="+", family="F1"):
    end = max(e[-1][1] for e in exons_by_tid.values())
    return dict(chrom=chrom, start=start, end=end, strand=strand, main=D.canonical(exons_by_tid), transcripts=exons_by_tid, gene_id=gid, gene_name=gid, biotype="protein_coding")

def test_merge_abutting():
    assert D.merge_abutting([(0, 10), (10, 20), (30, 40)]) == [(0, 20), (30, 40)]

def test_pick_row_takes_the_row_covering_the_member_not_the_longest():
    flank, member = row(0, 4000, 100, 4100), row(5000, 6000, 9000, 10000)
    assert D.pick_row([flank, member], "chr1", 5100, 5900) is member
    assert D.pick_row([row(5000, 6000, 0, 1000, tname="chr2_MATERNAL")], "chr1", 5100, 5900) is None      # not the homologous chromosome

def test_counterpart_keeps_the_isoforms_that_lift_and_flips_the_strand_on_a_minus_row():
    c = copy("g1", 1100, {"t1": [(1100, 1200), (1300, 1400)], "t2": [(1100, 1200), (1900, 2100)]})
    cp, why, _ = D.counterpart([row(0, 1000, 5000, 6000)], c, lo=1000)          # the row ends before t2's last exon
    assert why == "counterpart" and set(cp["isoforms"]) == {"t1"} and cp["isoforms"]["t1"] == [(5100, 5200), (5300, 5400)] and cp["strand"] == "+"
    cp, _, _ = D.counterpart([row(0, 10000, 0, 10000, strand="-")], c, lo=1000)
    assert cp["strand"] == "-" and set(cp["isoforms"]) == {"t1", "t2"}
    assert D.counterpart([row(0, 1000, 5000, 6000, mapq=5)], c, lo=1000)[:2] == (None, "mapq < 20")
    assert D.counterpart([], c, lo=1000)[:2] == (None, "no row on the homologous chromosome")
    assert D.counterpart([row(100, 150, 5000, 5050)], c, lo=1000)[:2] == (None, "no isoform lifts")

def test_collapsed_drops_same_family_lifts_that_overlap_on_the_haplotype_only():
    cp = lambda s, e, strand="+": dict(tname="chr1_MATERNAL", strand=strand, start=s, end=e, isoforms={})
    ref = {"a": ("chr1", 0, 100), "b": ("chr1", 500, 600), "c": ("chr1", 550, 650), "d": ("chr1", 900, 1000)}
    fam = {"a": "F1", "b": "F1", "c": "F1", "d": "F2"}
    assert D.collapsed({"a": cp(0, 100), "b": cp(20, 120)}, fam, ref) == {"a", "b"}
    assert D.collapsed({"b": cp(0, 100), "c": cp(20, 120)}, fam, ref) == set()            # they overlap on the reference too
    assert D.collapsed({"a": cp(0, 100), "d": cp(20, 120)}, fam, ref) == set()            # two families
    assert D.collapsed({"a": cp(0, 100), "b": cp(60, 160)}, fam, ref) == set()            # 40% of the shorter
    assert D.collapsed({"a": cp(0, 100), "b": cp(20, 120, "-")}, fam, ref) == set()       # opposite strands

FAM = {"sets": {"truth": [{"family": "F1", "copies": [
    copy("gA", 100, {"a1": [(100, 200), (300, 400), (500, 600)], "a2": [(100, 200), (500, 600)]}),
    copy("gB", 2100, {"b1": [(2100, 2200), (2300, 2400), (2500, 2600)]})]}]}}

def rec(start, end, chain=(), ms=100, primary=True, chrom="chr1", strand="+"):
    return Rec(chrom, start, end, strand, tuple(chain), ms, ms, primary)

def test_member_of_takes_the_member_with_most_exon_bases_on_the_strand():
    idx = D.member_index(FAM)
    assert D.member_of(rec(150, 550, [(200, 301), (400, 501)]), idx) == "gA"
    assert D.member_of(rec(150, 550, strand="-"), idx) is None and D.member_of(rec(700, 900), idx) is None
    two = {"sets": {"t": [{"family": "F", "copies": [copy("gZ", 100, {"z": [(100, 200)]}), copy("gY", 100, {"y": [(100, 200)]}), copy("gX", 150, {"x": [(150, 400)]})]}]}}
    assert D.member_of(rec(100, 200), D.member_index(two)) == "gY" and D.member_of(rec(150, 400), D.member_index(two)) == "gX"     # tie: lower id; else most bases

def test_isoform_weights_split_a_read_over_the_isoforms_it_fits():
    tx = FAM["sets"]["truth"][0]["copies"][0]["transcripts"]
    full, part, skip, novel = ((200, 300), (400, 500)), ((200, 300),), ((200, 500),), ((200, 250),)
    assert D.isoform_weights([full, part, skip, (), novel], tx) == {"a1": 2.0, "a2": 1.0}      # single-exon and novel reads carry none
    assert D.isoform_weights([(), novel], tx) == {"a1": 1.0}                                    # no read fits: the canonical isoform

def test_scenarios_move_only_the_near_tied_reads_of_a_group():
    idx = D.member_index(FAM)
    reads = {f"u{i}": [rec(150, 550, [(200, 301), (400, 501)])] for i in range(5)}                                      # decided, on gA
    reads.update({f"t{i}": [rec(150, 550, ms=100), rec(2150, 2550, ms=95, primary=False)] for i in range(3)})         # near-tied, primary on gA
    reads.update({"v0": [rec(2150, 2550, ms=100), rec(150, 550, ms=100, primary=False)]})                             # near-tied, primary on gB
    reads.update({"w0": [rec(2150, 2550, ms=100), rec(9000, 9400, ms=99, primary=False)]})                            # a non-member place in the band: left alone
    reads.update({"x0": [rec(2150, 2550, ms=100), rec(150, 550, ms=50, primary=False)]})                              # clear lead
    c = D.scenario_counts(reads, idx)
    assert c["observed"] == {"gA": 8, "gB": 3}
    assert c["to-major"] == {"gA": 9, "gB": 2} and c["to-minor"] == {"gA": 5, "gB": 6}

def test_cells_split_evenly_over_haplotypes_and_renormalise_over_lifted_isoforms():
    lifted = {"mat": {"a1", "a2"}, "pat": {"a2"}}
    assert D.cells(40, {"a1": 3.0, "a2": 1.0}, lifted, "a1") == {("mat", "a1"): 15, ("mat", "a2"): 5, ("pat", "a2"): 20}
    assert D.cells(10, {"a1": 1.0}, {"pat": {"a2"}}, "a1") == {("pat", "a2"): 10}      # no weighted isoform lifts: the canonical of those that do

def test_cut_redraws_until_enough_is_left_and_full_needs_every_junction():
    rng = random.Random(1)
    for _ in range(50):
        n5, n3 = D.cut(1000, [[0.5, 0.45], [0.1, 0.0]], rng)
        assert 1000 - n5 - n3 >= 200
    assert D.cut(250, [[0.5, 0.4]], rng) == (0, 0)                                    # no draw fits: the whole transcript
    assert D.is_full([100, 200], 300, 80, 80, False) and not D.is_full([100, 200], 300, 81, 0, False)
    assert not D.is_full([100, 200], 300, 0, 0, True) and not D.is_full([], 300, 0, 0, False)

PROF = {"classes": [1000, 2000, 3000, 5000], "end_pairs": {"0": [[0.0, 0.0]], "1": [], "2": [], "3": [], "4": []}, "rates": [[0.0, 0.0, 0.0]],
        "ins_len": {"1": 1.0}, "del_len": {"1": 1.0}, "ins_hp_frac": 0.0, "del_hp_frac": 0.0, "kept_intron_per_intron": 0.0}

def test_draw_reads_is_seeded_sense_and_carries_its_origin(tmp_path):
    genome = "".join(random.Random(7).choice("ACGT") for _ in range(3000))
    fetch = {"ref": lambda contig, a, b: genome[a:b]}
    fam = json.loads(json.dumps(FAM)); fam["sets"]["truth"][0]["copies"][1]["strand"] = "-"
    cps = D.reference_counterparts(fam)
    counts = {"gA": (4, {"a1": 1.0}), "gB": (2, {"b1": 1.0})}
    a = D.draw_reads(fam, cps, counts, PROF, fetch, seed=2026)
    assert a == D.draw_reads(fam, cps, counts, PROF, fetch, seed=2026) and a != D.draw_reads(fam, cps, counts, PROF, fetch, seed=1)
    assert len(a) == 6 and all("|" not in r["read"] for r in a) and len({r["read"] for r in a}) == 6
    sense = genome[100:200] + genome[300:400] + genome[500:600]
    rc = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))
    by = collections.defaultdict(list)
    for r in a:
        by[r["member"]].append(r)
    assert all(r["seq"] == sense and r["full"] == 1 and r["hap"] == "ref" and r["isoform"] == "a1" for r in by["gA"])
    assert all(r["seq"] == rc(genome[2100:2200] + genome[2300:2400] + genome[2500:2600]) for r in by["gB"])
    assert by["gA"][0]["blocks"] == "100-200,300-400,500-600" and by["gA"][0]["chrom"] == "chr1" and by["gB"][0]["strand"] == "-"

def test_errors_follow_the_drawn_read_rates_without_a_second_poisson():
    prof = dict(PROF, rates=[[0.01, 0.0, 0.0]])
    genome = "".join(random.Random(7).choice("ACGT") for _ in range(3000))
    reads = D.draw_reads(FAM, D.reference_counterparts(FAM), {"gA": (30, {"a1": 1.0})}, prof, {"ref": lambda c, a, b: genome[a:b]}, seed=3)
    assert {r["edits"] for r in reads} == {3}                                           # 300 bases x 0.01, every read
    assert all(len(r["seq"]) == 300 for r in reads)

def test_families_json_merges_abutting_exons_and_names_the_canonical():
    out = D.family_json([("F1", "g1", "G1", "chr1", "+", "protein_coding")], {"g1": {"t1": [(0, 10), (10, 20), (40, 50)], "t2": [(0, 20)]}})
    c = out["sets"]["truth"][0]["copies"][0]
    assert c["transcripts"]["t1"] == [[0, 20], [40, 50]] and c["main"] == "t1" and (c["start"], c["end"]) == (0, 50)

# ---- r1809: reads that copy a real read's record

def test_scenario_members_agree_with_the_counts():
    idx = D.member_index(FAM)
    reads = {f"t{i}": [rec(150, 550, ms=100), rec(2150, 2550, ms=95, primary=False)] for i in range(3)}
    reads.update({"v0": [rec(2150, 2550, ms=100), rec(150, 550, ms=100, primary=False)], "u0": [rec(150, 550)], "z0": [rec(9000, 9400)]})
    mem = D.scenario_members(reads, idx)
    assert mem["observed"] == {"t0": "gA", "t1": "gA", "t2": "gA", "v0": "gB", "u0": "gA"}                # z0 is on no member
    assert mem["to-major"]["v0"] == "gA" and mem["to-minor"]["t0"] == "gB" and mem["to-minor"]["u0"] == "gA"
    assert {s: dict(collections.Counter(m.values())) for s, m in mem.items()} == D.scenario_counts(reads, idx)

def test_templates_copy_the_record_at_the_member_the_read_goes_to():
    idx = D.member_index(FAM)
    reads = {"t0": [rec(150, 550, [(200, 301), (400, 501)], ms=100), rec(2150, 2560, [(2200, 2301)], ms=95, primary=False), rec(2150, 2500, ms=90, primary=False)],
             "v0": [rec(2150, 2550, ms=100), rec(150, 550, ms=100, primary=False)]}
    reads["t1"] = reads["t0"]                                                               # gA holds the most near-tied primaries
    t = {(s, q): (m, b) for s, q, m, b in D.templates(reads, idx)}
    assert t["observed", "t0"] == ("gA", [(150, 200), (300, 400), (500, 550)])              # the primary's blocks, split at introns only
    assert t["to-minor", "t0"] == ("gB", [(2150, 2200), (2300, 2560)])                      # its best-ms record on gB
    assert t["to-major", "v0"] == ("gA", [(150, 550)]) and t["to-major", "t0"] == t["observed", "t0"]

def test_a_moved_read_without_a_record_at_its_destination_stays():
    three = {"sets": {"t": [{"family": "F", "copies": [copy("gA", 100, {"a": [(100, 600)]}), copy("gB", 2100, {"b": [(2100, 2600)]}), copy("gC", 4100, {"c": [(4100, 4600)]})]}]}}
    idx = D.member_index(three)
    reads = {f"ab{i}": [rec(150, 550), rec(2150, 2550, primary=False)] for i in range(3)}
    reads.update({f"bc{i}": [rec(2150, 2550), rec(4150, 4550, primary=False)] for i in range(2)})
    reads.update({"c0": [rec(4150, 4550), rec(2150, 2550, primary=False)]})                  # one group gA-gB-gC; gA has most, gC fewest
    t = {(s, q): m for s, q, m, _ in D.templates(reads, idx)}
    assert t["to-major", "bc0"] == "gB" and t["to-major", "c0"] == "gC"                      # no record on gA: they stay
    assert t["to-minor", "ab0"] == "gA" and t["to-minor", "bc0"] == "gC"
    assert D.templates.stayed == {"to-major": 3, "to-minor": 3}                              # counted, as registered

def test_lift_many_is_lift_exons_for_every_block_list():
    for strand in "+-":
        r = row(0, 3000, 5000, 7990, strand=strand, cg="1000M10D2000M")
        lists = [[(100, 200), (300, 400)], [(950, 1050)], [(2900, 3100)], [(10, 20), (40, 50), (2000, 2100)]]
        assert D.lift_many(r, lists) == [D.DS.lift_exons(r, b) and [tuple(x) for x in D.DS.lift_exons(r, b)] for b in lists]
        assert D.lift_many(r, lists)[2] is None and D.lift_many(r, lists)[0] is not None    # outside the row; inside it
    assert D.lift_many(row(0, 3000, 5000, 8100, cg="1000M100D2000M"), [[(950, 1050)]]) == [None]      # a block growing by more than 50

def test_template_reads_are_the_lifted_blocks_in_the_sense_of_the_transcript():
    genome = "".join(random.Random(7).choice("ACGT") for _ in range(3000))
    rc = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))
    fam = json.loads(json.dumps(FAM)); fam["sets"]["truth"][0]["copies"][1]["strand"] = "-"
    tpl = [("r1", "gA", [(100, 200), (300, 400), (500, 600)]), ("r2", "gA", [(150, 200), (450, 580)]), ("r3", "gA", [(105, 200), (300, 400), (500, 600)]),
           ("r4", "gA", [(190, 200), (300, 400), (500, 600)]), ("r5", "gB", [(2100, 2200), (2300, 2400), (2500, 2600)]), ("r6", "gB", [(2900, 2950)])]
    lifts = {"gA": {"h1": ("c", "+", [[(b[0] + 1000, b[1] + 1000) for b in bl] for _, g, bl in tpl if g == "gA"]),
                    "h2": ("c", "+", [None] * 4)},
             "gB": {"h1": ("c", "-", [[(3100, 3200), (3300, 3400), (3500, 3600)], None])}}
    hap = "N" * 1000 + genome
    fetch = {"h1": lambda c, a, b: hap[a:b] if a >= 1000 else genome[a:b], "h2": lambda c, a, b: genome[a:b]}
    a = D.draw_template_reads(fam, tpl, lifts, PROF, fetch, seed=2026)
    assert a == D.draw_template_reads(fam, tpl, lifts, PROF, fetch, seed=2026) and len(a) == 5 and all("|" not in r["read"] for r in a)
    by = {r["template"]: r for r in a}
    assert "r6" not in by and all(r["hap"] == "h1" for r in a)                               # lifts nowhere: not drawn; h2 never lifts
    assert by["r1"]["seq"] == genome[100:200] + genome[300:400] + genome[500:600] and by["r2"]["seq"] == genome[150:200] + genome[450:580]
    assert by["r5"]["seq"] == rc(genome[2100:2200] + genome[2300:2400] + genome[2500:2600]) and by["r5"]["strand"] == "-"
    assert (by["r1"]["full"], by["r1"]["isoform"]) == (1, "a1") and by["r3"]["full"] == 1 and (by["r5"]["full"], by["r5"]["isoform"]) == (1, "b1")
    assert by["r2"]["full"] == 0 and by["r4"]["full"] == 0                                    # a chain of no isoform; an end block under 20 bases
    assert by["r1"]["blocks"] == "100-200,300-400,500-600" and by["r2"]["blocks"] == by["r1"]["blocks"]      # the member's exon union, as before
