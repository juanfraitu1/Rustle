"""Tests of mimic_bar.py, and of the pure parts of mimic_stats.py and read_stats.py (bench/famwise_eval)."""
import collections, os, subprocess, sys
import pysam
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mimic_bar as B, mimic_stats as M, read_stats as RS

TEN = list(B.FAMILIES)
OTHERS = ["CHM13_G9000001", "CHM13_G9000002"]
FAM_HEAD = "family\treads\tnear_ties\tshare\ttie\tdecidable\tno near-perfect hit\toutside\tcompetitor unlocatable\tprimary_wrong\tSEVERAL COPIES\tONE PLACE\n"

def run(script, *args):
    p = subprocess.run([sys.executable, os.path.join(HERE, script), *map(str, args)], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    return p.stdout

def real():
    """a held-out quarter: the ten families at 1000 / 10000, two more at 100 / 10000"""
    return {**{g: (1000, 10000) for g in TEN}, **{g: (100, 10000) for g in OTHERS}}

def tables(tmp_path, name, fams, wrong=(47, 497)):
    """PREFIX.fam.tsv + PREFIX.stats.tsv; fams = {first gene id: (near_ties, reads)}; the value column is junk on purpose"""
    with open(tmp_path / f"{name}.fam.tsv", "w") as fh:
        fh.write(FAM_HEAD)
        for i, (g, (k, n)) in enumerate(fams.items()):
            fh.write(f"{g},CHM13_X{i:07d} (+3)\t{n}\t{k}\t{k / n:.4f}" + "\t0" * 8 + "\n")
    with open(tmp_path / f"{name}.stats.tsv", "w") as fh:
        fh.write("stat\tvalue\tk\tn\tlo\thi\nstratum_S2\t0.5\t1\t2\t0\t1\nfam_p90\t0.1510\t\t154\t0.06\t0.25\n")
        if wrong:
            fh.write(f"S2_decidable_primary_wrong\tNA\t{wrong[0]}\t{wrong[1]}\t0\t1\n")
    return str(tmp_path / name)

def bar(tmp_path, sim, *extra, held=None, **kw):
    out = run("mimic_bar.py", "--seen", tables(tmp_path, "seen", real()), "--held", tables(tmp_path, "held", held or real()),
              "--sim", tables(tmp_path, "sim", sim, **kw), *extra)
    lines = out.splitlines()
    assert lines[0].split("\t") == ["item", "detail", "seen", "held", "sim", "verdict"]
    return [l.split("\t") for l in lines[1:-1]], lines[-1]

def verdict(rows, item):
    return next(r[5] for r in rows if r[0] == item)

FULL = ("--ends", "0.709", "0.700", "--replay", "Rustle=1.04:1.00", "--replay", "StringTie=1.00:1.02")

# ---- the two-way test of a share
def test_overlap_passes():
    assert B.same(10, 100, 12, 100) == (True, "overlap")

def test_factor_passes_without_overlap():
    assert not B.overlap(1000, 10000, 1400, 10000)
    assert B.same(1000, 10000, 1400, 10000) == (True, "factor 1.40")
    assert B.same(1400, 10000, 1000, 10000) == (True, "factor 1.40")

def test_both_fail():
    assert B.same(1000, 10000, 2000, 10000) == (False, "disjoint, factor 2.00")
    assert not B.same(1000, 10000, 1500, 10000)[0]            # a factor of exactly 1.5 is not less than 1.5

def test_zero_shares():
    assert B.same(0, 10000, 0, 500)[0]                         # both zero
    assert B.same(0, 10, 1, 10) == (True, "overlap")           # one zero: the interval test alone
    assert B.same(0, 10000, 50, 10000) == (False, "disjoint, one share zero")

# ---- F
def test_seven_of_ten_fails_and_eight_passes(tmp_path):
    sim = real(); sim.update({g: (2000, 10000) for g in TEN[:3]})
    rows, last = bar(tmp_path, sim, *FULL)
    assert verdict(rows, "F") == "fail" and [r[5] for r in rows if r[0].startswith("F") and r[0] != "F"] == ["fail"] * 3 + ["pass"] * 7
    assert last == "BAR fail: F"
    sim[TEN[2]] = (1000, 10000)
    rows, last = bar(tmp_path, sim, *FULL)
    assert verdict(rows, "F") == "pass" and last == "BAR pass"

def test_missing_family_fails(tmp_path):
    sim = real(); del sim[TEN[5]]
    rows, _ = bar(tmp_path, sim)
    row = next(r for r in rows if r[1].startswith(TEN[5]))
    assert row[4] == "missing" and row[5] == "fail" and verdict(rows, "F") == "pass"
    for g in TEN[6:8]:
        del sim[g]
    assert verdict(bar(tmp_path, sim)[0], "F") == "fail"

def test_family_rows_show_seen_held_sim(tmp_path):
    rows, _ = bar(tmp_path, real())
    fam = [r for r in rows if r[0].startswith("F") and r[0] != "F"]
    assert [r[1].split()[0] for r in fam] == TEN and all(r[2].startswith("1000/10000 0.1000") and r[3] == r[2] == r[4] for r in fam)

# ---- S
def test_pooled_share_leaves_out_exactly_the_first_four(tmp_path):
    sim = real(); sim.update({g: (9000, 10000) for g in TEN[:4]})
    rows, _ = bar(tmp_path, sim)
    s = next(r for r in rows if r[0] == "S")
    assert s[3].startswith("6200/80000 ") and s[4].startswith("6200/80000 ") and s[5] == "pass"      # the four do not count
    sim[TEN[4]] = (5000, 10000)                                                                        # the fifth does
    s = next(r for r in bar(tmp_path, sim)[0] if r[0] == "S")
    assert s[4].startswith("10200/80000 ") and s[5] == "fail"

# ---- W
def test_w_reads_k_and_n(tmp_path):
    rows, _ = bar(tmp_path, real(), wrong=(47, 600))
    w = next(r for r in rows if r[0] == "W")
    assert w[3].startswith("47/497 0.0946") and w[4].startswith("47/600 0.0783") and w[5] == "pass"
    assert verdict(bar(tmp_path, real(), wrong=(200, 500))[0], "W") == "fail"

def test_w_without_decidable_reads_is_not_given(tmp_path):
    for wrong in (None, (0, 0)):                               # 0 of 0 has the interval 0-1: it would pass against anything
        rows, last = bar(tmp_path, real(), *FULL, wrong=wrong)
        assert verdict(rows, "W") == "not given" and last == "BAR incomplete: not given W"

# ---- E and A
def test_ends_at_the_edge(tmp_path):
    assert verdict(bar(tmp_path, real(), "--ends", "0.709", "0.809")[0], "E") == "pass"    # 0.809 - 0.709 > 0.1 in binary floats
    assert verdict(bar(tmp_path, real(), "--ends", "0.709", "0.609")[0], "E") == "pass"
    assert verdict(bar(tmp_path, real(), "--ends", "0.709", "0.8091")[0], "E") == "fail"

def test_ends_from_read_stats_tables(tmp_path):
    for name, v in (("h.tsv", "0.7090"), ("s.tsv", "0.9000")):
        (tmp_path / name).write_text(f"stat\tvalue\nreads\t10\nboth_ends_strict\t{v}\n")
    e = next(r for r in bar(tmp_path, real(), "--ends", tmp_path / "h.tsv", tmp_path / "s.tsv")[0] if r[0] == "E")
    assert e[3:] == ["0.7090", "0.9000", "fail"]

def test_replay_at_the_edge(tmp_path):
    rows, last = bar(tmp_path, real(), "--ends", "0.7", "0.7", "--replay", "Rustle=1.04:1.14", "--replay", "StringTie=1.00:0.90")
    assert verdict(rows, "A") == "pass" and last == "BAR pass"
    rows, last = bar(tmp_path, real(), "--ends", "0.7", "0.7", "--replay", "Rustle=1.04:1.15", "--replay", "StringTie=1.00:0.90")
    assert verdict(rows, "A") == "fail" and last == "BAR fail: A"
    assert [(r[0], r[1].split()[0], r[5]) for r in rows if r[0] in ("A1", "A2")] == [("A1", "Rustle", "fail"), ("A2", "StringTie", "pass")]

def test_replay_needs_both_registered_tools(tmp_path):
    rows, last = bar(tmp_path, real(), "--ends", "0.7", "0.7", "--replay", "Rustle=1.04:1.04")
    assert verdict(rows, "A") == "not given" and last == "BAR incomplete: not given A"

# ---- the last line
def test_not_given_reads_incomplete(tmp_path):
    rows, last = bar(tmp_path, real())
    assert verdict(rows, "E") == verdict(rows, "A") == "not given" and verdict(rows, "F") == verdict(rows, "S") == verdict(rows, "W") == "pass"
    assert last == "BAR incomplete: not given E, A"

def test_a_failed_item_decides_even_when_others_are_not_given(tmp_path):
    assert bar(tmp_path, real(), wrong=(200, 500))[1] == "BAR fail: W; not given E, A"

def test_fam_table_missing_is_not_given(tmp_path):
    sim = tables(tmp_path, "kb", real()); os.remove(sim + ".fam.tsv")
    out = run("mimic_bar.py", "--seen", tables(tmp_path, "seen", real()), "--held", tables(tmp_path, "held", real()), "--sim", sim)
    assert out.splitlines()[-1] == "BAR incomplete: not given F, S, E, A"

# ---- mimic_stats.py
def test_wilson():
    lo, hi = M.wilson(47, 497)
    assert (f"{lo:.4f}", f"{hi:.4f}") == ("0.0719", "0.1235") and M.wilson(0, 0) == (0.0, 1.0) and M.wilson(0, 50)[0] == 0.0

def test_primary_wrong():
    assert M.wrong({"picks": "primary=1;eichler-keep=0", "origin_loci": "0,2"}) and not M.wrong({"picks": "primary=2;eichler-keep=0", "origin_loci": "0,2"})

def test_famstats():
    F = [collections.Counter(n=100, S2=k) for k in (60, 20, 10, 6, 4, 0)]
    s = M.famstats(F)
    assert s["pooled"] == 100 / 600 and s["half_in"] == 1 and s["zero"] == 1 / 6 and s["ge50"] == 1 / 6
    assert s["pooled_minus_top1"] == 40 / 500 and s["pooled_minus_top4"] == 4 / 200 and s["top4_of_S2"] == 0.96

# ---- read_stats.py
def test_blocks_split_at_introns_only():
    assert RS.blocks(100, [(4, 5), (7, 50), (2, 3), (8, 1), (1, 4), (7, 46), (3, 900), (7, 100)]) == [(100, 200), (1100, 1200)]

def test_ends_of_a_read_on_a_transcript():
    ex = [(100, 200), (1100, 1200), (2100, 2300)]              # 400 bases
    assert RS.judge("+", 150, 2200, ex) == (True, True, 50, 100, 0, 400, 3)
    assert RS.judge("+", 1100, 2200, ex)[:5] == (False, True, 100, 100, 1)
    assert RS.judge("-", 150, 2200, ex)[:5] == (True, True, 100, 50, 0)
    assert RS.judge("-", 150, 1200, ex)[:5] == (False, True, 200, 50, 1)
    assert RS.judge("+", 150, 1200, ex)[:2] == (True, False)

def seg(header, name, start, cigar, flag=0):
    a = pysam.AlignedSegment(header)
    a.query_name, a.reference_id, a.reference_start, a.flag, a.mapping_quality, a.cigartuples = name, 0, start, flag, 60, cigar
    a.query_sequence = "A" * sum(l for op, l in cigar if op in (1, 4, 7, 8))
    a.set_tag("AS", 100, "i"); a.set_tag("ms", 100, "i")
    return a

def test_read_stats_table(tmp_path):
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": "c1", "LN": 10000}]})
    with pysam.AlignmentFile(str(tmp_path / "sel.bam"), "wb", header=header) as fh:
        fh.write(seg(header, "full", 150, [(7, 50), (3, 900), (7, 99), (8, 1), (3, 900), (7, 100)]))    # 250 bases, 1 mismatch
        fh.write(seg(header, "full", 5000, [(7, 250)], 256))
        fh.write(seg(header, "cut", 1100, [(7, 100), (3, 900), (7, 50), (2, 2), (7, 50)]))               # 200 bases, 1 deletion
    gtf = tmp_path / "ref.gtf"
    gtf.write_text("".join(f'c1\tx\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "g"; transcript_id "t1";\n' for s, e in ((100, 200), (1100, 1200), (2100, 2300))))
    out = run("read_stats.py", tmp_path / "sel.bam", gtf, "--out", tmp_path / "st.tsv", "--pct", 100)
    assert "strict (best-overlap transcript): 5' in first exon 0.500, 5' truncated 0.500; 3' in last exon 1.000; both 0.500" in out
    t = dict(l.rstrip("\n").split("\t") for l in open(tmp_path / "st.tsv"))
    assert t["stat"] == "value" and (t["reads"], t["read_length_p50"], t["secondaries_per_read"], t["share_without_secondaries"]) == ("2", "250", "0.5000", "0.5000")
    assert (t["first_exon_5p_strict"], t["last_exon_3p"], t["both_ends_strict"]) == ("0.5000", "1.0000", "0.5000")
    assert t["events_per_aligned_base"] == f"{2 / 450:.6f}"
