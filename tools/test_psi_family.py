"""tools/psi_family.py: introns keyed by the exon columns of one member; `differs` needs disjoint wide ranges and reads."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from psi_family import load_columns, key_rows, HEADER

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
