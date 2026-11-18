"""Tests of sim_chains.py (bench/famwise_eval): the chain table of row 1798 on a two-family toy library."""
import csv, json, os, sys
import pytest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sim_chains as S

A1 = [[100, 200], [300, 400], [500, 600]]
A2 = [[100, 200], [300, 400], [700, 1000]]
A4 = [[100, 200], [250, 400]]
B1 = [[100, 200], [300, 400]]
B2 = [[100, 200], [500, 1000]]
SIL = [[5000, 5100], [5200, 5300]]

def copy(gene, chrom, strand, start, end, **tx):
    return dict(chrom=chrom, start=start, end=end, strand=strand, gene_id=gene, gene_name=gene.lower(), main=next(iter(tx)), transcripts=tx)

FAM = {"sets": {"toy": [
    {"family": "F1", "copies": [copy("A", "chr1", "+", 100, 1000, a1=A1, a2=A2, a3=[[100, 1000]], a4=A4),
                                copy("A2", "chr1", "+", 50, 600, x1=[[50, 200], [300, 400], [500, 600]]),   # a1's chain again
                                copy("S", "chr1", "+", 5000, 6000, s1=SIL + [[5900, 6000]]),
                                copy("T", "chr1", "+", 150, 900, t1=[[150, 200], [800, 900]])]},            # silent, on emitting A
    {"family": "F2", "copies": [copy("B", "chr2", "-", 100, 1000, b1=B1, b2=B2),
                                copy("U", "chr2", "-", 5000, 6000, u1=[[5000, 5100], [5900, 6000]]),
                                copy("V", "chr2", "+", 100, 1000, v1=[[100, 200], [900, 1000]])]}]}}        # silent: B is on the other strand
# read, member, isoform, full, exact row (status, picks, origin_loci) or None = not in the table
READS = [("r1", "A", "a1", 1, ("decidable", "primary=0", "0")), ("r2", "A", "a1", 1, ("decidable", "primary=1;eichler-keep=NA", "0,1")),
         ("r3", "A", "a1", 0, None),                                # not full: no count, and its lost primary affects nothing
         ("r4", "A2", "x1", 1, ("decidable", "primary=0", "0")),
         ("r5", "A", "a2", 1, ("decidable", "primary=1", "0")),      # wrong primary
         ("r6", "A", "a2", 1, ("decidable", "primary=0", "0")),
         ("r7", "A", "a3", 1, None), ("r8", "A", "a3", 1, None),      # single exon: no chain
         ("r9", "B", "b1", 1, ("decidable", "primary=0", "0")), ("r10", "B", "b1", 1, ("outside", "primary=0", "")),
         ("r11", "B", "b2", 1, None),                               # one full read is no reference chain
         ("r12", "A", "a4", 1, ("decidable", "primary=0", "0")), ("r13", "A", "a4", 1, None)]   # unmapped
GTFS = {"rustle:P": [("chr1", "+", A1), ("chr1", "+", A2), ("chr1", "+", SIL), ("chr1", "+", [[4000, 4100], [5200, 5300]]),
                     ("chr1", "+", [[5000, 6000]]), ("chr2", "+", [[5000, 5100], [5200, 5300]])],
        "rustle:O": [("chr1", "+", A1), ("chr1", "+", A2), ("chr1", "+", A4), ("chr2", "-", B1)],
        "stringtie:P": [("chr1", "+", A1), ("chr1", "+", SIL), ("chr1", "+", [[5000, 5100], [5200, 5300], [5400, 5500]])],
        "stringtie:O": [("chr2", "-", B1)]}

def gtf(path, txs):
    with open(path, "w") as fh:
        for i, (chrom, strand, exons) in enumerate(txs):
            fh.writelines(f'{chrom}\tx\texon\t{s + 1}\t{e}\t.\t{strand}\t.\tgene_id "g{i}"; transcript_id "t{i}";\n' for s, e in exons)
    return str(path)

def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def run(tmp_path, gtfs=GTFS, extra=()):
    fam, origin, exact, out = (str(tmp_path / n) for n in ("fam.json", "origin.tsv", "exact.tsv", "out"))
    json.dump(FAM, open(fam, "w"))
    member = {c["gene_id"]: (f["family"], c) for f in FAM["sets"]["toy"] for c in f["copies"]}
    with open(origin, "w") as o, open(exact, "w") as x:
        o.write("read\tchrom\tstrand\tblocks\tmember\tfamily\thap\tisoform\tfull\tmore\n")
        x.write("read\tstratum\tlead\tstatus\tpicks\torigin_loci\tloci\n")
        for q, m, t, full, ex in READS:
            o.write(f"{q}\t{member[m][1]['chrom']}\t{member[m][1]['strand']}\t.\t{m}\t{member[m][0]}\tmat\t{t}\t{full}\tz\n")
            if ex:
                x.write(f"{q}\tS1\t0\t{ex[0]}\t{ex[1]}\t{ex[2]}\t.\n")
    S.main([fam, origin, exact, out, *(f"{k}={gtf(tmp_path / (k.replace(':', '_') + '.gtf'), v)}" for k, v in gtfs.items()), *extra])
    by_intron = {r["introns"]: r for r in rows(out + ".chains.tsv")}
    return by_intron, {(r["tool"], r["subset"]): r for r in rows(out + ".summary.tsv")}, {r["member"]: r for r in rows(out + ".false.tsv")}

def test_key_is_the_one_chains_gives(tmp_path):
    c = FAM["sets"]["toy"][1]["copies"][0]
    assert {S.chain_key(c, B1)} == S.chains(gtf(tmp_path / "b.gtf", [("chr2", "-", B1)])) == {("chr2", "-", ((200, 300),))}
    assert S.chain_key(c, [[100, 150], [150, 200], [300, 400]]) == S.chain_key(c, B1)    # abutting exons are one exon
    assert S.chain_key(c, [[100, 1000]]) is None

@pytest.mark.parametrize("n,cls", [(2, "2-4"), (4, "2-4"), (5, "5-19"), (19, "5-19"), (20, "20+")])
def test_depth_classes(n, cls):
    assert S.depth(n) == cls

def test_at_origin(tmp_path):
    p = tmp_path / "e.tsv"
    p.write_text("read\tstratum\tlead\tstatus\tpicks\torigin_loci\tloci\n" + "".join(f"{q}\tS1\t0\t{st}\t{pk}\t{o}\t.\n" for q, st, pk, o in [
        ("ok", "decidable", "primary=1", "0,1"), ("wrong", "decidable", "primary=1", "0,10"), ("tie", "tie", "primary=0", "0"),
        ("out", "outside", "primary=0", ""), ("na", "decidable", "primary=NA", "0"), ("extra", "decidable", "primary=-1", "1")]))
    assert S.at_origin(str(p)) == {"ok"}

def test_chain_table(tmp_path):
    ch, _, _ = run(tmp_path)
    assert sorted(ch) == ["200-250", "200-300", "200-300,400-500", "200-300,400-700"]    # no single-exon a3, no one-read b2
    a1, a2, a4, b1 = ch["200-300,400-500"], ch["200-300,400-700"], ch["200-250"], ch["200-300"]
    # full reads only, summed over the two members that share the chain; every full read at home = not affected
    assert (a1["full_reads"], a1["off_reads"], a1["affected"], a1["members"], a1["families"], a1["depth"]) == ("3", "0", "0", "A,A2", "F1", "2-4")
    assert (a1["chrom"], a1["strand"], a1["first_intron_start"], a1["last_intron_end"]) == ("chr1", "+", "200", "500")
    for r in (a2, a4, b1):                                         # wrong primary, absent from the table, outside
        assert (r["full_reads"], r["off_reads"], r["affected"]) == ("2", "1", "1")
    assert [a2[k] for k in GTFS] == ["1", "1", "0", "0"] and [b1[k] for k in GTFS] == ["0", "1", "0", "1"]

def test_summary(tmp_path):
    _, sm, _ = run(tmp_path, {**GTFS, "rustle:D": GTFS["rustle:P"]})
    r = sm["rustle", "all"]
    assert list(r) == ["tool", "subset", "reference", "affected", "P", "O", "D", "C", "O_over_P", "D_over_P", "ref_P", "ref_O", "ref_D"]
    assert list(r.values())[2:] == ["4", "3", "1", "3", "1", "3.00", "3.00", "1.00", "2", "4", "2"]
    assert sm["rustle", "2-4"] == {**r, "subset": "2-4"} and sm["rustle", "20+"]["reference"] == "0"
    s = sm["stringtie", "all"]                                     # P recovers no affected chain; stringtie has no D
    assert list(s.values())[2:] == ["4", "3", "0", "1", "NA", "inf", "NA", "NA", "1", "1", "NA"]
    assert [k[1] for k in sm if k[0] == "rustle"] == ["all", "2-4", "5-19", "20+"]

def test_drop_families(tmp_path):
    _, sm, _ = run(tmp_path, extra=["--drop-families", "X9,F2"])
    r = sm["rustle", "without dropped families"]
    assert (r["reference"], r["affected"], r["P"], r["O"], r["C"], r["O_over_P"]) == ("3", "2", "1", "2", "2.00", "2.00")
    assert sm["rustle", "all"]["reference"] == "4"

def test_false_copies(tmp_path, capsys):
    _, _, fc = run(tmp_path)
    # S: the chain inside its span; not the one reaching outside it, not the single-exon one. T lies on emitting A.
    assert sorted(fc) == ["S", "TOTAL", "U", "V"] and (fc["S"]["gene_name"], fc["S"]["family"]) == ("s", "F1")
    assert [fc["S"][k] for k in GTFS] == ["1", "0", "2", "0"]
    assert [fc["U"][k] for k in GTFS] == ["0", "0", "0", "0"]      # the chr2 chain is on the other strand
    assert [fc["TOTAL"][k] for k in GTFS] == ["1", "0", "2", "0"]
    out = capsys.readouterr().out
    assert "TOTAL\t\t\t1\t0\t2\t0\n" in out and "silent members\t3\n" in out and "rustle\tall\t4\t3\t1\t3\t3.00\t3.00\t2\t4\n" in out

def test_missing_p_fails(tmp_path):
    with pytest.raises(SystemExit) as e:
        run(tmp_path, {"rustle:P": GTFS["rustle:P"], "stringtie:O": GTFS["stringtie:O"]})
    assert "stringtie" in str(e.value) and "P" in str(e.value)
