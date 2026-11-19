"""tests of undecidable.py (PREREG r1820); run: python3 -m pytest bench/famwise_eval/test_undecidable.py"""
import os, sys, gzip
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import records as R
import undecidable as U

def rec(chrom, s, e, ms, primary=False, strand="+", AS=None):
    return R.Rec(chrom, s, e, strand, (), ms, ms if AS is None else AS, primary)

def truth(kind, loci, origin):
    return R.Truth("S2", "tie" if kind != "one" else "decidable", kind, frozenset(origin), loci)

LOCI = [("c1", 100, 200, "+", 1000), ("c1", 5000, 5100, "+", 1000), ("c2", 0, 100, "+", 990)]

def test_u_star_takes_origin_records_off_the_primary_place():
    rs = [rec("c1", 100, 200, 500, True), rec("c1", 5000, 5100, 500), rec("c2", 0, 100, 495)]
    t = truth("several", LOCI, {0, 1})
    assert U.u_star(rs, t) == [1]                    # c2 is not an origin
    t2 = truth("several", LOCI, {0, 1, 2})
    assert U.u_star(rs, t2) == [1, 2]

def test_u_star_skips_a_second_record_of_the_primary_place():
    rs = [rec("c1", 100, 200, 500, True), rec("c1", 120, 190, 500), rec("c1", 5000, 5100, 500)]
    assert U.u_star(rs, truth("several", LOCI, {0, 1})) == [2]

def test_u_ms_only_on_exact_ties():
    rs = [rec("c1", 100, 200, 500, True), rec("c1", 5000, 5100, 500), rec("c2", 0, 100, 499)]
    assert U.u_ms(rs) == [1]
    rs2 = [rec("c1", 100, 200, 500, True), rec("c1", 5000, 5100, 499)]
    assert U.u_ms(rs2) == []

def test_merge_spans_by_strand_and_overlap():
    m = U.merge_spans([("c1", "+", 0, 10, "a"), ("c1", "+", 5, 20, "b"), ("c1", "-", 5, 20, "c"), ("c1", "+", 20, 30, "d")])
    assert [(v[2], v[3], v[4]) for v in m if v[1] == "+"] == [(0, 20, ["a", "b"]), (20, 30, ["d"])]

def test_landing_copies_and_classes():
    tr = {"u1": truth("several", LOCI, {0, 1}), "u2": truth("several", LOCI, {0, 1}),
          "d_home": truth("one", LOCI, {1}), "d_else": truth("one", LOCI, {0})}
    recs = {"u1": [rec("c1", 100, 200, 500, True), rec("c1", 5000, 5100, 500)],
            "u2": [rec("c1", 100, 200, 500, True), rec("c1", 5010, 5090, 500)],
            "d_home": [rec("c1", 5000, 5100, 500, True)],
            "d_else": [rec("c1", 100, 200, 500, True), rec("c1", 5000, 5100, 400)]}
    cps = U.landing_copies(tr, recs)
    assert len(cps) == 1 and cps[0]["start"] == 5000 and cps[0]["u_reads"] == 2 and cps[0]["u_records"] == 2
    st = U.classify(cps, tr, recs, own=2, shadow=1)
    v = st[cps[0]["id"]]
    assert (v["dd_origin"], v["dd_power_any"], v["dd_power_10"], v["cls"], v["cls_10"]) == (1, 1, 0, "own_1", "own_1")
    st = U.classify(cps, {k: x for k, x in tr.items() if k != "d_home"}, recs, own=2, shadow=1)
    assert st[cps[0]["id"]]["cls"] == "shadow" and st[cps[0]["id"]]["cls_10"] == "unjudgeable"

def test_contiguous_part():
    big = ((1, 2), (3, 4), (5, 6))
    assert U.contiguous_part(((3, 4), (5, 6)), big)
    assert not U.contiguous_part(((1, 2), (5, 6)), big)
    assert not U.contiguous_part(((3, 4),), big)     # one intron is not enough

def test_load_annot_gff3_and_ncbi_pseudogene(tmp_path):
    p = tmp_path / "a.gff3.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("c1\tCAT\tgene\t101\t400\t.\t+\t.\tID=G1;gene_name=AAA;gene_biotype=unprocessed_pseudogene\n"
                 "c1\tCAT\ttranscript\t101\t400\t.\t+\t.\tID=T1;Parent=G1;transcript_biotype=unprocessed_pseudogene\n"
                 "c1\tCAT\texon\t101\t150\t.\t+\t.\tID=T1.e1;Parent=T1\n"
                 "c1\tCAT\texon\t201\t250\t.\t+\t.\tID=T1.e2;Parent=T1\n"
                 "c1\tCAT\texon\t301\t400\t.\t+\t.\tID=T1.e3;Parent=T1\n"
                 "c1\tRefSeq\tpseudogene\t1001\t1200\t.\t-\t.\tID=gene-P;Name=P;gene_biotype=pseudogene\n"
                 "c1\tRefSeq\texon\t1001\t1050\t.\t-\t.\tID=id-P;Parent=gene-P\n"
                 "c1\tRefSeq\texon\t1101\t1200\t.\t-\t.\tID=id-P-2;Parent=gene-P\n")
    tx, genes = U.load_annot(str(p))
    assert U.chain(tx["T1"]["exons"]) == ((150, 200), (250, 300)) and tx["T1"]["gene"] == "G1"
    assert genes["G1"]["biotype"] == "unprocessed_pseudogene" and U.group(genes["G1"]["biotype"]) == "unprocessed_pseudogene"
    assert U.chain(tx["gene-P"]["exons"]) == ((1050, 1100),) and U.group("pseudogene") == "pseudogene_unsplit"
    cps = [dict(id="x", chrom="c1", strand="+", start=120, end=380), dict(id="y", chrom="c1", strand="-", start=990, end=1300)]
    ann = U.copy_annotation(cps, tx, genes)
    assert ann["x"]["gene"] == "AAA" and ann["x"]["chains"] == {((150, 200), (250, 300))} and ann["x"]["exonic"] == 200
    assert ann["y"]["group"] == "pseudogene_unsplit"

def test_sam_hits_coordinates_and_scores():
    l = "t1\t0\tcopy|c1|1000\t11\t60\t10=100N5=1X4=\t*\t0\t0\t*\t*"
    (q, tn, rev, exons, cov, ident), = U.sam_hits([l])
    assert exons == [(10, 20), (120, 130)] and cov == 1.0 and abs(ident - 0.95) < 1e-9 and not rev

def test_bars_c1_counts_and_cost():
    ann_chain = ((10, 20), (30, 40))
    cls = {f"s{i}": "shadow" for i in range(25)} | {"o": "own"}
    ann = {c: dict(chains={ann_chain}) for c in cls}
    A = {x: {c: (set(), []) for c in cls} for x in ("P", "PU", "PROJ", "D", "PUms")}
    for i in range(25):
        A["PU"][f"s{i}"][0].add(ann_chain)
        if i < 5:
            A["PROJ"][f"s{i}"][0].add(ann_chain); A["D"][f"s{i}"][0].add(ann_chain)
        if i < 20:
            A["PUms"][f"s{i}"][0].add(ann_chain)
    A["PU"]["o"][0].add(((11, 20), (30, 40)))
    L = dict(l.split("\t") for l in U.bars(cls, ann, A))
    assert L["C1 added depicted copies (PU not P)"] == "25" and L["C1"] == "HOLDS"
    assert L["C1b added depicted copies not depicted by PROJ"] == "20" and L["C1b"] == "HOLDS"
    assert L["C2"] == "NEW OVER DEFAULTS" and L["C3"] == "HOLDS"
    assert L["cost PU unmatched chains added on own copies"] == "1"
    A["PU"]["s0"][0].add(((99, 100), (200, 300)))
    L = dict(l.split("\t") for l in U.bars(cls, ann, A))
    assert L["C1 depicted share of chains PU adds"].startswith("25/26")

def test_stem_and_find_arrays():
    assert U.stem("TSPY10") == "TSPY" and U.stem("TSPY6P") == "TSPY" and U.stem("USP17L12") == "USP17L"
    assert U.stem("RBMY1A1") == "RBMY1A" and U.stem("AC114498.1") is None and U.stem("LOFF_G1") is None
    genes = {f"g{i}": dict(chrom="chrY", strand="+", start=i * 20000, end=i * 20000 + 3000, name=f"TSPY{i}", biotype="protein_coding")
             for i in range(1, 7)}
    genes["far"] = dict(chrom="chrY", strand="+", start=5_000_000, end=5_003_000, name="TSPY20", biotype="protein_coding")
    arr = U.find_arrays(genes, 5)
    assert len(arr) == 1 and len(arr[0][4]) == 6 and arr[0][0].startswith("TSPY@chrY")
