"""tools/psi_bunches.py: introns from CIGARs, placed <= max per band, PSI ranges hold psi_placed, bunches by shared site."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from psi_bunches import introns_of, gtf_introns, tally, bunch_rows, psi, summary, HEADER, intron_rows, INTRON_HEADER

A, B, C, X = ("c1", "+", 101, 200), ("c1", "+", 101, 300), ("c1", "+", 251, 300), ("c2", "+", 11, 90)

def test_introns_of_and_gtf(tmp_path):
    # 100M 100N 50M 2D 48M 50N 10M from 0-based 0: introns 101-200 and 301-350 (S / I move no reference base)
    assert introns_of(0, [(4, 5), (0, 100), (3, 100), (0, 50), (2, 2), (1, 3), (0, 48), (3, 50), (0, 10)]) == [(101, 200), (301, 350)]
    g = tmp_path / "g.gtf"
    ex = lambda s, e: f'c1\tr\texon\t{s}\t{e}\t.\t+\t.\tgene_id "G"; transcript_id "T";\n'
    g.write_text(ex(201, 250) + ex(1, 100) + ex(301, 400))
    assert gtf_introns(str(g)) == {("c1", "+", 101, 200): "G", ("c1", "+", 251, 300): "G"}

def test_tally_bands():
    solo = {A: 10, B: 2}
    tied = [("r1", True, (A,)), ("r1", True, (X,)),      # exact tie with another locus: tied in both bands
            ("r2", True, (A,)), ("r2", False, ())  ,      # near tie with an unspliced record: placed in exact, tied in wide
            ("r3", True, (A, C)), ("r3", True, (B,)),     # two records of one locus, different introns
            ("r4", True, (A,)), ("r4", True, (A,)),       # both records carry it: placed, counted once
            ("r5", False, (A,))]                          # primary on a contig not scanned: skipped
    c = tally(solo, tied)
    assert c[A] == [12, 14, 11, 14] and c[B] == [2, 3, 2, 3] and c[C] == [0, 1, 0, 1] and c[X] == [0, 1, 0, 1]
    assert all(v[0] <= v[1] and v[2] <= v[3] and v[2] <= v[0] and v[1] <= v[3] for v in c.values())

def test_bunches_and_ranges():
    cnt = {A: [12, 14, 11, 14], B: [2, 3, 2, 3], C: [0, 1, 0, 1], X: [5, 5, 5, 5]}
    rows = bunch_rows(cnt, {A: "G1"}, min_reads=3)
    assert [r[0] for r in rows] == ["c1:+:D:101"] * 2           # A, B share the donor; C is below 3 reads, X is alone
    a, b = rows
    assert a[6:9] == ["G1", 1, 2] and b[6:9] == ["G1", 0, 2]    # B takes the bunch's locus, not in the GTF
    assert len(a) == len(HEADER) and a[9:11] == [12, 14]
    for r in rows:
        for o in (9, 14):
            assert r[o + 3] <= r[o + 2] <= r[o + 4]
    assert abs(a[11] - 12 / 14) < 1e-9 and abs(a[12] - 12 / 15) < 1e-9 and abs(a[13] - 14 / 16) < 1e-9
    assert bunch_rows({("c1", "-", 101, 200): [3] * 4, ("c1", "-", 101, 300): [3] * 4}, {})[0][:2] == ["c1:-:A:101", "acceptor"]
    assert psi(0, 0, 0, 0) == (None, None, None) and psi(0, 4, 3, 7) == (0.0, 0.0, 4 / 7)
    assert summary(rows)[0].startswith("exact: 1 bunches, 0 ") and summary(bunch_rows({A: [0, 4, 0, 4], B: [3] * 4}, {}))[0].startswith("exact: 1 bunches, 1 ")

def test_intron_rows_keep_what_bunches_drop():
    cnt = {A: [12, 14, 11, 14], B: [2, 3, 2, 3], C: [0, 1, 0, 1], X: [5, 5, 5, 5], ("c1", "-", 5, 9): [0, 0, 0, 0]}
    rows = intron_rows(cnt)
    assert [tuple(r[:4]) for r in rows] == [("c1", "+", 101, 200), ("c1", "+", 101, 300), ("c1", "+", 251, 300), ("c2", "+", 11, 90)]
    assert rows[2][4:] == [0, 1, 0, 1] and len(INTRON_HEADER) == len(rows[0]) == 8   # C: under the 3-read floor, kept here

def test_gtf_introns_keep_contigs_apart(tmp_path):
    from psi_bunches import gtf_introns
    ex = lambda c, s, e: f'{c}\tx\texon\t{s}\t{e}\t.\t+\t.\tgene_id "g"; transcript_id "t1";\n'
    (tmp_path / "g.gtf").write_text(ex("c1", 1, 100) + ex("c1", 201, 300) + ex("c2", 151, 180) + ex("c2", 401, 500))   # one transcript id on two contigs
    assert set(gtf_introns(str(tmp_path / "g.gtf"))) == {("c1", "+", 101, 200), ("c2", "+", 181, 400)}
