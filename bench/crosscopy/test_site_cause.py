"""bench/crosscopy/site_cause.py: the junction's site strength against its competitors', at the copy that uses it and the one that does not."""
import math, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from site_cause import advantage, verdict, main, pool
from test_calls import write

def build_two_copies(tmp_path, break_acceptor_on=None, truncate_s=None, slip_s=0, decoy=False):
    """two copies (c1, c2; plus strand) of one gene: exon 1 61-160, intron 161-360, exon 2 361-420, intron 421-620, exon 3 621-720.
    break_acceptor_on: that copy's exon-2 acceptor and the bases around it are replaced by a run of C.
    Default: c1 includes exon 2 in 60 of 100 reads and skips it in 40, c2 only skips. With truncate_s the competitor is an
    alternative acceptor at 330 instead of the skip, and contig c2 ends at that base (so its exon-2 acceptor has no window).
    With slip_s, c2 carries that many extra bases inside intron 1 which the pair alignment does not show (one `=` block):
    a junction carried over from c1 lands that many bases short of c2's real acceptor."""
    rnd = random.Random(7); R = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
    intron = lambda: "GTAAGT" + R(174) + "TTTTTTTTTTCTCTTTTCAG"
    seq = R(60) + R(100) + intron() + R(60) + intron() + R(100) + R(60)
    assert len(seq) == 780 and seq[160:162] == "GT" and seq[358:360] == "AG" and seq[618:620] == "AG"
    if decoy:      # both copies: no AG at the exon-2 acceptor, a full acceptor 4 bases downstream
        seq = seq[:344] + "TTTTTTTTTTCTCTTTTCAG" + seq[364:]
    c = {"c1": seq, "c2": seq}
    if break_acceptor_on:
        c[break_acceptor_on] = seq[:346] + "C" * 22 + seq[368:]      # no AG within 6 bases of the acceptor either
    if truncate_s:
        c["c2"] = seq[:truncate_s]
    if slip_s:
        c["c2"] = seq[:250] + "C" * slip_s + seq[250:]
    fa = tmp_path / "g.fa"; fa.write_text("".join(f">{k}\n{v}\n" for k, v in c.items()))
    n2 = len(c["c2"]) - slip_s; m1, m2 = "c1:1-780", f"c2:1-{n2}"
    psi = write(tmp_path, [("F", m1, "donor", "C0", "C1", "c1", "+", 161, 360, 60, 100, 0.6, 0.6, 0.6, 0),
                           ("F", m2, "donor", "C0", "C1", "c2", "+", 161, 360, 0, 100, 0.0, 0.0, 0.0, 1)])
    comp = (161, 330) if truncate_s else (161, 620)
    dump = tmp_path / "i.tsv"
    dump.write_text("chrom\tstrand\tstart\tend\tplaced\tmax\tplaced_w\tmax_w\n" + "".join(
        f"{ch}\t+\t{s}\t{e}\t{n}\t{n}\t{n}\t{n}\n" for ch, s, e, n in [("c1", 161, 360, 60), ("c1", *comp, 40), ("c1", 421, 620, 60), ("c2", comp[0], comp[1] + slip_s, 100)]))
    paf = tmp_path / "p.paf"; paf.write_text(f"{m1}\t780\t0\t{n2}\t+\t{m2}\t{n2}\t0\t{n2}\t{n2 - 2}\t{n2}\t0\tcg:Z:{n2}=\n")
    ex = lambda ch, t, s, e: f'{ch}\tr\texon\t{s}\t{e}\t.\t+\t.\tgene_id "G{ch}"; transcript_id "{t}";\n'
    gtf = tmp_path / "g.gtf"
    gtf.write_text(ex("c1", "T1", 61, 160) + ex("c1", "T1", 361, 420) + ex("c1", "T1", 621, 720) + ex("c1", "T2", 61, 160) + ex("c1", "T2", 621, 720)
                   + ("" if truncate_s else ex("c2", "T3", 61, 160) + ex("c2", "T3", 621, 720)))
    return str(fa), psi, str(dump), str(paf), str(gtf)

def read_rows(p):
    lines = open(p).read().splitlines()
    return [dict(zip(lines[0].split("\t"), l.split("\t"))) for l in lines[1:]]

def test_advantage_and_verdict():
    assert advantage(8.0, [3.0, 6.5]) == 1.5 and advantage(8.0, []) == 8.0
    nan = float("nan")       # a competitor without a window: not located, wherever it stands in the list
    assert all(math.isnan(advantage(*a)) for a in ((5.0, [3.0, nan]), (5.0, [nan, 3.0]), (nan, [3.0])))
    assert verdict(1.0) == "AGREES" and verdict(-1.0) == "OPPOSES" and verdict(0.99) == "NO SITE DIFFERENCE"
    assert verdict(float("nan")) == "UNLOCATED" and verdict(0.6, tau=0.5) == "AGREES"

def test_broken_acceptor_at_the_silent_copy(tmp_path):
    fa, psi, dump, paf, gtf = build_two_copies(tmp_path, break_acceptor_on="c2")
    main(["score", psi, dump, paf, gtf, fa, str(tmp_path / "o.tsv")])
    r = read_rows(tmp_path / "o.tsv")
    assert len(r) == 1 and r[0]["verdict_1"] == "AGREES" and float(r[0]["dA"]) >= 1 and r[0]["canonical_lost"] == "1"
    assert r[0]["n_competitors"] == "1" and r[0]["dropped_competitors"] == "0" and abs(float(r[0]["pair_identity"]) - 778 / 780) < 1e-4

def test_identical_copies_have_no_site_difference(tmp_path):
    fa, psi, dump, paf, gtf = build_two_copies(tmp_path)
    main(["score", psi, dump, paf, gtf, fa, str(tmp_path / "o.tsv")])
    r = read_rows(tmp_path / "o.tsv")[0]
    assert r["verdict_1"] == "NO SITE DIFFERENCE" and float(r["dA"]) == 0 and r["canonical_lost"] == "0"

def test_site_at_contig_end_is_unlocated(tmp_path):
    fa, psi, dump, paf, gtf = build_two_copies(tmp_path, truncate_s=365)       # c2 ends 5 bases after its exon-2 acceptor
    main(["score", psi, dump, paf, gtf, fa, str(tmp_path / "o.tsv")])
    assert read_rows(tmp_path / "o.tsv")[0]["verdict_1"] == "UNLOCATED"

def test_competitor_that_does_not_project(tmp_path):
    fa, psi, dump, paf, gtf = build_two_copies(tmp_path, break_acceptor_on="c2")
    open(paf, "w").write("c1:1-780\t780\t0\t500\t+\tc2:1-780\t780\t0\t500\t500\t500\t0\tcg:Z:500=\n")   # exon 3 not aligned: the skip (40% of U's reads) has no image
    main(["score", psi, dump, paf, gtf, fa, str(tmp_path / "o.tsv")])
    assert read_rows(tmp_path / "o.tsv")[0]["verdict_1"] == "UNLOCATED"
    d = open(dump).read().replace("c1\t+\t161\t620\t40\t40\t40\t40", "c1\t+\t161\t620\t5\t5\t5\t5").replace("c2\t+\t161\t620\t100\t100\t100\t100", "c2\t+\t161\t330\t100\t100\t100\t100")
    open(dump, "w").write(d)                                                    # now 5 of 65 reads at U (under 10%): dropped, the call is judged
    main(["score", psi, dump, paf, gtf, fa, str(tmp_path / "o.tsv")])
    r = read_rows(tmp_path / "o.tsv")[0]
    assert r["dropped_competitors"] == "1" and r["verdict_1"] in ("AGREES", "OPPOSES", "NO SITE DIFFERENCE")

def test_pool_counts_units_once_and_leaves_the_development_sample_out(tmp_path):
    head = "family cls gap u s unit reads_u dA verdict_0.5 verdict_1 verdict_2 pair_identity".split()
    def tab(name, rows):
        p = tmp_path / name
        p.write_text("\t".join(head) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows)); return str(p)
    # a member's name is its locus span in ONE sample's assembly: every sample names the same copies differently
    row = lambda fam, uc, sc, start, reads, v, cls="LARGE": [fam, cls, 0.5, f"{uc}:{reads}-9", f"{sc}:{reads}-8", f"{uc}:+:{start}-{start + 99}|{sc}:+:{start}", reads, 2.0, v, v, v, 0.99]
    dev = tab("dev.tsv", [row("F", "c1", "c2", 100, 50, "AGREES")])
    kb = tab("kb.tsv", [row("F", "c1", "c2", 100, 80, "AGREES"), row("G", "c1", "c2", 900, 30, "OPPOSES")])                 # first = a development unit
    t1 = tab("t1.tsv", [row("H", "c1", "c2", 100, 30, "AGREES"), row("H", "c2", "c1", 100, 30, "AGREES"),                   # U and S swapped: two units
                        row("I", "c1", "c2", 500, 30, "OPPOSES"), row("J", "c3", "c4", 100, 30, "NO SITE DIFFERENCE"), row("M", "c7", "c8", 100, 30, "MISKEYED")])
    t2 = tab("t2.tsv", [row("I", "c1", "c2", 500, 90, "AGREES"), row("K", "c5", "c6", 100, 30, "NO SITE DIFFERENCE", cls="SAME")])
    res = pool(("GGO", dev), [("GGO", kb), ("PTR", t1), ("PTR", t2)])
    assert res["units_large"] == 6 and res["agrees"] == 3 and res["opposes"] == 1          # G opposes; H twice; I once, from the sample with 90 reads; M counted, not judged
    assert res["site_diff_large"] == (4, 5) and res["site_diff_same"] == (0, 1) and res["miskeyed"] == 1

def test_one_rule_for_both_copies_of_a_pair(tmp_path):
    # identical copies whose junction has a poor acceptor and a full one 4 bases on: read at c1, carried over at c2. The 6-base
    # search must run on both or on neither, or identical sequences differ.
    fa, psi, dump, paf, gtf = build_two_copies(tmp_path, decoy=True)
    main(["score", psi, dump, paf, gtf, fa, str(tmp_path / "o.tsv")])
    r = read_rows(tmp_path / "o.tsv")[0]
    assert float(r["dA"]) == 0 and r["verdict_1"] == "NO SITE DIFFERENCE" and r["unit"] == "c1:+:161-360|c2:+:161"

def test_miskeyed_call_gets_no_verdict(tmp_path):
    fa, psi, dump, paf, gtf = build_two_copies(tmp_path, break_acceptor_on="c2")
    t = open(psi).read().replace("c2\t+\t161\t360", "c2\t+\t161\t330"); assert t != open(psi).read(); open(psi, "w").write(t)   # S's row under the key is another intron
    open(dump, "a").write("c2\t+\t161\t360\t50\t50\t50\t50\nc2\t+\t161\t330\t5\t5\t5\t5\n")                    # and S uses the image of U's junction
    main(["score", psi, dump, paf, gtf, fa, str(tmp_path / "o.tsv")])
    assert read_rows(tmp_path / "o.tsv")[0]["verdict_1"] == "MISKEYED"

def test_a_slip_of_the_carried_junction_is_not_a_lost_site(tmp_path):
    # c2 = c1 with 2 extra intronic bases the alignment hides: the image of c1's junction ends 2 bases before c2's intact acceptor
    fa, psi, dump, paf, gtf = build_two_copies(tmp_path, slip_s=2)
    main(["score", psi, dump, paf, gtf, fa, str(tmp_path / "o.tsv")])
    r = read_rows(tmp_path / "o.tsv")[0]
    assert r["verdict_1"] == "NO SITE DIFFERENCE" and r["canonical_lost"] == "0" and abs(float(r["dA"])) < 1

def test_events_beside_a1_count_linked_units_once(tmp_path):
    # a key gives one call, between its two extreme copies, so two samples of a genome can state one event as two units
    head = "family cls gap u s unit reads_u dA verdict_0.5 verdict_1 verdict_2 pair_identity".split()
    def tab(name, rows):
        p = tmp_path / name
        p.write_text("\t".join(head) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows)); return str(p)
    row = lambda unit, gap, v, cls="LARGE": ["F", cls, gap, "u", "s", unit, 30, 2.0, v, v, v, 0.99]
    dev = tab("dev.tsv", [])
    a = tab("a.tsv", [row("c1:+:101-200|c2:+:101", 0.5, "AGREES"), row("c9:+:101-200|c8:+:101", 0.4, "OPPOSES")])
    b = tab("b.tsv", [row("c1:+:101-200|c3:+:101", 0.9, "AGREES"),          # the same junction of U against another copy
                      row("c4:+:101-200|c3:+:101", 0.3, "OPPOSES"),         # another U against that copy: linked through S's site
                      row("c5:+:101-200|c6:+:101", 0.6, "NO SITE DIFFERENCE")])
    c = tab("c.tsv", [row("c1:+:101-200|c2:+:101", 0.5, "OPPOSES")])        # the same strings on another genome: another event
    res = pool(("GGO", dev), [("PTR", a), ("PTR", b), ("PPA", c)])
    assert (res["agrees"], res["opposes"]) == (2, 3)
    assert res["events"] == (4, 1, 2)       # 4 events (three PTR units are one); of those whose speaking unit (largest gap) has a site difference: AGREES 1, OPPOSES 2
