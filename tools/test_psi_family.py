"""tools/psi_family.py: introns keyed by the exon columns of one member; `differs` needs disjoint wide ranges and reads."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from psi_family import load_columns, key_rows, HEADER, load_pair_records, project_intron

def setup(tmp_path):
    # member A (+): exons 1-100 (C0), 201-300 (C1), 401-500 (C2); member B (-) on c2: 1-100 (C2), 201-300 (C1), 401-500 (C0)
    (tmp_path / "l.gff3").write_text("c1\t.\tgene\t1\t500\t.\t+\t.\tID=a\nc2\t.\tgene\t1\t500\t.\t-\t.\tID=b\n")
    (tmp_path / "c.tsv").write_text("family\tmembers\tcolumn\tclass\tn_loci\tn_exons\tcollision\tlen_ratio\texons\n" + "".join(
        f"F1\t2\tC{i}\tcore\t2\t2\t0\t1.00\tc1:1-500#{i}:{a}-{b},c2:1-500#{2 - i}:{c}-{d}\n"
        for i, (a, b, c, d) in enumerate([(0, 100, 400, 500), (200, 300, 200, 300), (400, 500, 0, 100)])))
    return load_columns(str(tmp_path / "c.tsv"), str(tmp_path / "l.gff3"))

def row(bunch, anchor, chrom, strand, s, e, placed, lo, hi):
    return dict(bunch=bunch, anchor=anchor, chrom=chrom, strand=strand, start=str(s), end=str(e), locus="G", in_gtf="1", placed=str(placed), max=str(placed),
                psi_placed="NA", psi_lo="NA", psi_hi="NA", placed_w=str(placed), max_w=str(placed), psi_placed_w="NA", psi_lo_w=str(lo), psi_hi_w=str(hi))

def test_keys_and_differs(tmp_path):
    cols = setup(tmp_path)
    assert cols[("c2", "-")][0][:2] == (1, 100) and cols[("c1", "+")][1][4] == "C1"
    rows = [row("a", "donor", "c1", "+", 101, 200, 18, 0.9, 0.95), row("a", "donor", "c1", "+", 101, 400, 2, 0.05, 0.1),      # A: C0->C1 0.9, C0->C2 skip
            row("b", "donor", "c2", "-", 301, 400, 4, 0.1, 0.3), row("b", "donor", "c2", "-", 101, 400, 26, 0.7, 0.9),        # B (minus): C0->C1 0.1-0.3, C0->C2
            row("x", "donor", "c1", "+", 101, 150, 9, 0.5, 0.5),                                                             # acceptor side in no exon
            row("y", "donor", "c1", "-", 101, 200, 9, 0.5, 0.5)]                                                             # antisense: no member on that strand
    out = key_rows(rows, cols, min_bunch=20)
    assert len(out) == 4 and all(len(r) == len(HEADER) for r in out)
    by = {(r[1], r[5]): r for r in out}
    assert by[("c1:1-500", "C1")][:8] == ["F1", "c1:1-500", "G", "donor", "C0", "C1", "core", "core"]
    assert by[("c2:1-500", "C1")][4] == "C0" and by[("c2:1-500", "C2")][10:12] == ["101", "400"]
    assert all(r[24] == 2 and r[25] == 1 and r[23] in (20, 30) for r in out)
    assert all(r[25] == 0 for r in key_rows(rows, cols, min_bunch=25))                 # A's bunch holds 20 certain reads only
    extra = rows + [row("a2", "donor", "c1", "+", 91, 200, 5, 1.0, 1.0)]               # a shallower bunch of A (shifted donor in C0): not summed
    e = key_rows(extra, cols)
    assert [r[25] for r in e] == [1, 1, 1, 1, 0] and e[4][23] == 5
    rows[2] = row("b", "donor", "c2", "-", 301, 400, 4, 0.1, 0.92)                     # ranges now overlap for C0->C1
    d = {(r[1], r[5]): r[25] for r in key_rows(rows, cols)}
    assert d[("c1:1-500", "C1")] == 0 and d[("c1:1-500", "C2")] == 1

def test_project_intron_plus_minus_and_gap(tmp_path):
    paf = tmp_path / "p.paf"   # c1:1-500 vs c2:1-500, 500 bases, minus strand, one 500= block; then a second pair with a gap
    paf.write_text("c1:1-500\t500\t0\t500\t-\tc2:1-500\t500\t0\t500\t500\t500\t0\tcg:Z:500=\n"
                   "c1:1-500\t500\t0\t500\t+\tc3:1-490\t490\t0\t490\t490\t500\t0\tcg:Z:95=10I395=\n")
    R = load_pair_records(str(paf), {"c1:1-500", "c2:1-500", "c3:1-490"})
    assert project_intron(R, "c1:1-500", "c2:1-500", 101, 200) == (301, 400)      # minus: base x maps to 501 - x
    assert project_intron(R, "c2:1-500", "c1:1-500", 301, 400) == (101, 200)      # and back
    assert project_intron(R, "c1:1-500", "c3:1-490", 201, 300) == (191, 290)      # past the 10-base insertion in c1
    assert project_intron(R, "c1:1-500", "c3:1-490", 101, 200) is None            # exonic base 100 sits in the insertion
    assert project_intron(R, "c1:1-500", "c9:1-9", 101, 200) is None              # no record for the pair

def test_old_columns_unchanged_and_new_ones(tmp_path):
    cols = setup(tmp_path)
    rows = [row("a", "donor", "c1", "+", 101, 200, 18, 0.9, 0.95), row("a", "donor", "c1", "+", 101, 400, 2, 0.05, 0.1),
            row("b", "donor", "c2", "-", 301, 400, 4, 0.1, 0.3), row("b", "donor", "c2", "-", 101, 400, 26, 0.7, 0.9)]
    out = key_rows(rows, cols, min_bunch=20)
    assert len(out) == 4 and all(len(r) == len(HEADER) == 31 for r in out)
    by = {(r[1], r[5]): r for r in out}
    a = by[("c1:1-500", "C1")]
    assert a[24:26] == [2, 1] and a[26:31] == [0, 2, 0.9, 0.95, 0.6]            # 0.9 - 0.3
    assert by[("c2:1-500", "C1")][30] == -0.85                                   # 0.1 - 0.95: below its partner

def test_silent_copy_gets_a_row(tmp_path):
    cols = setup(tmp_path)                                                       # B is A reversed: x on c1 = 501 - x on c2
    (tmp_path / "p.paf").write_text("c1:1-500\t500\t0\t500\t-\tc2:1-500\t500\t0\t500\t500\t500\t0\tcg:Z:500=\n")
    rec = load_pair_records(str(tmp_path / "p.paf"), {"c1:1-500", "c2:1-500"})
    rows = [row("a", "donor", "c1", "+", 101, 200, 18, 0.6, 0.6), row("a", "donor", "c1", "+", 101, 400, 12, 0.4, 0.4)]   # A skips C1 at 0.4
    intr = {("c2", "-", 301, 400): [40, 40, 40, 40]}                             # B: only C0->C1 from that donor, never the skip
    out = key_rows(rows, cols, min_bunch=20, introns=intr, records=rec)
    add = [r for r in out[2:]]
    skip = next(r for r in add if r[5] == "C2")
    assert skip[1] == "c2:1-500" and skip[10:12] == ["101", "400"] and skip[26] == 1 and skip[28:30] == [0.0, 0.0]
    incl = next(r for r in add if r[5] == "C1")                                  # B's constitutive junction: PSI 1, not silent
    assert incl[26] == 0 and incl[28:30] == [1.0, 1.0] and incl[18] == 40
    a_skip = next(r for r in out[:2] if r[5] == "C2")
    assert a_skip[27] == 2 and a_skip[30] == 0.4                                 # eligible with the added row; gap 0.4 - 0
    assert key_rows(rows, cols, min_bunch=20, introns={("c2", "-", 301, 400): [5, 5, 5, 5]}, records=rec)[2:] == []   # 5 reads: not eligible

def test_unlocated_and_existing_member(tmp_path):
    cols = setup(tmp_path)
    (tmp_path / "p.paf").write_text("c1:1-500\t500\t0\t90\t-\tc2:1-500\t500\t410\t500\t90\t90\t0\tcg:Z:90=\n")    # only exon C0 aligned
    rec = load_pair_records(str(tmp_path / "p.paf"), {"c1:1-500", "c2:1-500"})
    rows = [row("a", "donor", "c1", "+", 101, 200, 18, 0.6, 0.6), row("a", "donor", "c1", "+", 101, 400, 12, 0.4, 0.4)]
    st = {}
    out = key_rows(rows, cols, min_bunch=20, introns={("c2", "-", 301, 400): [40, 40, 40, 40]}, records=rec, stats=st)
    # B's own C0->C1 junction needs no image; the skip has neither an own intron nor an image: not located, A alone under that key
    assert len(out) == 3 and st == {"unlocated": 1} and out[2][5] == "C1" and out[2][26] == 0 and out[0][27] == 2 and out[0][30] == -0.4
    assert out[1][5] == "C2" and out[1][27] == 1 and out[1][30] == "NA"
    both = rows + [row("b", "donor", "c2", "-", 301, 400, 30, 0.9, 0.9), row("b", "donor", "c2", "-", 101, 400, 3, 0.1, 0.1)]
    (tmp_path / "q.paf").write_text("c1:1-500\t500\t0\t500\t-\tc2:1-500\t500\t0\t500\t500\t500\t0\tcg:Z:500=\n")
    rec = load_pair_records(str(tmp_path / "q.paf"), {"c1:1-500", "c2:1-500"})
    assert len(key_rows(both, cols, min_bunch=20, introns={}, records=rec)) == 4                   # B has rows under both keys: nothing added

FULL = "c1:1-500\t500\t0\t500\t-\tc2:1-500\t500\t0\t500\t500\t500\t0\tcg:Z:500=\n"     # x on c1 = 501 - x on c2

def test_own_shifted_junction_is_not_silent(tmp_path):
    # B joins C0 to C1 through an acceptor 7 bases inside its C1 exon (294-400, alone at its sites, so the psi stage kept no row);
    # the image of A's junction (301-400) has no read. B uses the junction between the two columns: its own intron is its row.
    cols = setup(tmp_path); (tmp_path / "p.paf").write_text(FULL)
    rec = load_pair_records(str(tmp_path / "p.paf"), {"c1:1-500": "+", "c2:1-500": "-"})
    rows = [row("a", "donor", "c1", "+", 101, 200, 30, 1.0, 1.0)]
    out = key_rows(rows, cols, min_bunch=20, introns={("c2", "-", 294, 400): [35, 35, 35, 35]}, records=rec)
    assert len(out) == 2 and out[1][1] == "c2:1-500" and out[1][10:12] == ["294", "400"] and out[1][18] == 35
    assert out[1][26] == 0 and out[1][28:30] == [1.0, 1.0] and out[0][27] == 2 and out[0][30] == 0.0

def test_image_off_the_copys_exon_gives_no_row(tmp_path):
    cols = setup(tmp_path)                                                       # this record puts A's 101-200 on c2's 251-350:
    (tmp_path / "p.paf").write_text("c1:1-500\t500\t0\t450\t-\tc2:1-500\t500\t0\t450\t450\t450\t0\tcg:Z:450=\n")   # 351 is in no exon of B
    rec = load_pair_records(str(tmp_path / "p.paf"), {"c1:1-500": "+", "c2:1-500": "-"})
    st = {}
    out = key_rows([row("a", "donor", "c1", "+", 101, 200, 30, 1.0, 1.0)], cols, min_bunch=20, introns={("c2", "-", 120, 350): [40, 40, 40, 40]}, records=rec, stats=st)
    assert len(out) == 1 and st == {"off_exon": 1}

def test_record_against_the_members_strands_is_not_used(tmp_path):
    (tmp_path / "p.paf").write_text(FULL.replace("\t-\t", "\t+\t") + FULL)       # A is +, B is -: only the - record relates them
    rec = load_pair_records(str(tmp_path / "p.paf"), {"c1:1-500": "+", "c2:1-500": "-"})
    assert len(rec[("c1:1-500", "c2:1-500")]) == 1 and project_intron(rec, "c1:1-500", "c2:1-500", 101, 200) == (301, 400)
