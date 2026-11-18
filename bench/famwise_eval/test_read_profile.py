"""Tests of read_profile.py (bench/famwise_eval): a tiny random genome, eight hand-made genes and hand-made =/X reads."""
import importlib.util, json, os, random, subprocess, sys
import pysam
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import read_profile as RP

FULL = "300=200N300="                 # a whole two-exon transcript of 600 bases
# gene id, biotype, contig, strand, transcripts as exon lists (0-based half-open)[, n_transcripts of the table]
GENES = [("G1", "protein_coding", "chr1", "+", [[(2000, 2300), (2500, 2800)]]),
         ("G2", "protein_coding", "chr1", "-", [[(5000, 5300), (5500, 5800)]]),
         ("M1", "protein_coding", "chr1", "+", [[(7000, 7300), (7500, 7800)]]),                                  # a family member
         ("L1", "lncRNA", "chr1", "+", [[(8500, 8800), (9000, 9300)]]),
         ("T2", "protein_coding", "chr1", "+", [[(10000, 10300), (10500, 10800)], [(10000, 10300), (10600, 10800)]], 1),   # two in the GFF3
         ("N2", "protein_coding", "chr1", "+", [[(13500, 13800), (14000, 14300)]], 2),                               # two in the table
         ("S1", "protein_coding", "chr1", "+", [[(12000, 12600)]]),                                              # one exon
         ("X1", "protein_coding", "chrX", "+", [[(1000, 1300), (1500, 1800)]])]                                  # not in --chroms
KEYS = {"source", "genes", "reads", "classes", "end_pairs", "rates", "ins_len", "del_len", "ins_hp_frac", "del_hp_frac",
        "kept_intron_per_intron", "sub_error", "ins_event", "del_event", "read_length", "full5_frac", "class_reads"}

def world(tmp_path, reads):
    """reads = (name, contig, start, cigar[, flag[, mapq]]) -> dict of the six paths the script takes"""
    p = {k: str(tmp_path / v) for k, v in dict(bam="lib.bam", genes="genes.tsv", gff="a.gff3", fa="g.fa", fam="fam.json", out="out.json").items()}
    rng, sizes = random.Random(0), {"chr1": 16000, "chrX": 4000}
    with open(p["fa"], "w") as fh:
        fh.writelines(f">{c}\n{''.join(rng.choice('ACGT') for _ in range(n))}\n" for c, n in sizes.items())
    pysam.faidx(p["fa"])
    rows, gff = ["gene_id\tgene_name\tgene_biotype\tsource\tchrom\tstart0\tend\tstrand\tn_transcripts\texon_blocks\n"], []
    for gid, bt, c, strand, txs, *n in GENES:
        s, e = min(x[0] for t in txs for x in t), max(x[1] for t in txs for x in t)
        rows.append(f"{gid}\t{gid}n\t{bt}\tCAT\t{c}\t{s}\t{e}\t{strand}\t{(n or [len(txs)])[0]}\t.\n")
        gff.append((c, s, f"{c}\tCAT\tgene\t{s + 1}\t{e}\t.\t{strand}\t.\tID={gid};gene_name={gid}n\n"))
        for i, t in enumerate(txs):
            gff.append((c, s, f"{c}\tCAT\ttranscript\t{t[0][0] + 1}\t{t[-1][1]}\t.\t{strand}\t.\tID={gid}.t{i};Parent={gid}\n"))
            gff += [(c, a, f"{c}\tCAT\texon\t{a + 1}\t{b}\t.\t{strand}\t.\tID={gid}.t{i}.e{a};Parent={gid}.t{i}\n") for a, b in t]
    open(p["genes"], "w").writelines(rows)
    open(p["gff"], "w").writelines(l for _, _, l in sorted(gff, key=lambda x: x[:2]))
    p["gff"] = pysam.tabix_index(p["gff"], preset="gff", force=True)
    json.dump({"sets": {"truth": [{"family": "F", "copies": [{"gene_id": "M1", "chrom": "chr1"}]}]}}, open(p["fam"], "w"))
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": c, "LN": n} for c, n in sizes.items()]}
    with pysam.AlignmentFile(p["bam"], "wb", header=hdr) as w:
        for q, c, start, cig, *rest in sorted(reads, key=lambda r: (r[1], r[2])):      # stable: ties keep the given order
            a = pysam.AlignedSegment()
            a.query_name, a.reference_id, a.reference_start, a.cigarstring = q, list(sizes).index(c), start, cig
            a.flag, a.mapping_quality = rest + [0, 60][len(rest):]
            a.query_sequence = "A" * a.infer_query_length()
            a.set_tag("AS", 100); a.set_tag("ms", 100)
            w.write(a)
    pysam.index(p["bam"])
    return p

def sampled(p, gene="G1", cap=50, seed=1):
    members = RP.load_members(p["fam"])
    with pysam.TabixFile(p["gff"]) as gff, pysam.AlignmentFile(p["bam"]) as bam, pysam.FastaFile(p["fa"]) as fa:
        g, exons = next(x for x in RP.gene_set(p["genes"], gff, members, ["chr1"]) if x[0]["gene_id"] == gene)
        return RP.gene_reads(bam, fa, g, exons, cap, seed)

def test_class_index_at_the_bounds():
    assert [RP.klass(n) for n in (1, 999, 1000, 1999, 2000, 2999, 3000, 4999, 5000, 90000)] == [0, 0, 1, 1, 2, 2, 3, 3, 4, 4]

def test_gene_set(tmp_path):
    p = world(tmp_path, [("r", "chr1", 2000, FULL)])
    with pysam.TabixFile(p["gff"]) as gff:
        ids = lambda chroms: [g["gene_id"] for g, _ in RP.gene_set(p["genes"], gff, RP.load_members(p["fam"]), chroms)]
        assert ids(["chr1"]) == ["G1", "G2"]              # member, lncRNA, two transcripts (GFF3 or table), one exon: out
        assert ids(["chr1", "chrX"]) == ["G1", "G2", "X1"]
        assert RP.gene_set(p["genes"], gff, set(), ["chr1"])[0][1] == [(2000, 2300), (2500, 2800)]
        assert "M1" in [g["gene_id"] for g, _ in RP.gene_set(p["genes"], gff, set(), ["chr1"])]

def test_qualifying_reads(tmp_path):
    p = world(tmp_path, [("ok", "chr1", 2000, FULL), ("mapq59", "chr1", 2000, FULL, 0, 59), ("secondary", "chr1", 2000, FULL, 256),
                         ("supplementary", "chr1", 2000, FULL, 2048), ("antisense", "chr1", 2000, FULL, 16),
                         ("at_flank", "chr1", 1000, "100=900N" + FULL), ("past_flank", "chr1", 999, "100=901N" + FULL),
                         ("at_end", "chr1", 2000, FULL + "900N100="), ("past_end", "chr1", 2000, FULL + "901N100="),
                         ("intronic", "chr1", 2350, "100=")])
    assert sorted(rd["name"] for rd in sampled(p)) == ["at_end", "at_flank", "ok"]      # intronic qualifies, has no transcript end
    assert sampled(p, "G2") == []

def test_alleles_come_from_the_full_pileup(tmp_path):
    """3 of 10 reads carry a substitution at 2100 (the individual's allele); `both` carries a second one alone (an error)"""
    carrier, both = "100=1X199=200N300=", "100=1X99=1X99=200N300="
    reads = [(f"p{i}", "chr1", 2000, FULL) for i in range(7)] + [("c1", "chr1", 2000, carrier), ("c2", "chr1", 2000, carrier), ("both", "chr1", 2000, both)]
    p = world(tmp_path, reads)
    assert {rd["name"]: rd["sub"] for rd in sampled(p) if rd["sub"]} == {"both": 1}
    one = {rd["name"]: rd["sub"] for seed in range(40) for rd in sampled(p, cap=1, seed=seed)}
    assert one["both"] == 1 and one["c1"] == 0 and one["c2"] == 0        # a sample of one read still knows the allele

def test_allele_needs_three_reads_and_a_fifth_of_the_coverage(tmp_path):
    two = [("a", "chr1", 2000, "100=1X199=200N300="), ("b", "chr1", 2000, "100=1X199=200N300="), ("c", "chr1", 2000, FULL)]
    assert sum(rd["sub"] for rd in sampled(world(tmp_path, two))) == 2                   # 2 reads of 3: errors
    thin = [(f"c{i}", "chr1", 2000, "100=1X199=200N300=") for i in range(3)] + [(f"p{i}", "chr1", 2000, FULL) for i in range(13)]
    assert sum(rd["sub"] for rd in sampled(world(tmp_path, thin))) == 3                  # 3 of 16 = 18.75%: errors
    fifth = thin[:-1]
    assert sum(rd["sub"] for rd in sampled(world(tmp_path, fifth))) == 0                 # 3 of 15 = 20%: the allele

def test_rates_ends_and_kept_introns(tmp_path):
    p = world(tmp_path, [("ev", "chr1", 2000, "100=1X99=2I100=200N3D297="),       # 597 aligned bases; one event of each kind
                         ("short", "chr1", 2030, "270=200N240="),                 # 30 bases missing at 5', 60 at 3'
                         ("kept", "chr1", 2000, "800="),                          # the intron retained
                         ("minus", "chr1", 5060, "240=200N300=", 16)])            # G2 is on -: its left end is the 3' end
    rows = sampled(p) + sampled(p, "G2")
    prof = RP.profile(rows, "lib.bam", 2)
    assert prof["genes"] == 2 and prof["reads"] == 4 and prof["classes"] == [1000, 2000, 3000, 5000]
    assert prof["rates"] == [[round(1 / 597, 6)] * 3, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    assert prof["end_pairs"] == {"0": [[0.0, 0.0], [0.0, 0.0], [0.05, 0.1], [0.0, 0.1]], "1": [], "2": [], "3": [], "4": []}
    assert prof["class_reads"] == {"0": 4, "1": 0, "2": 0, "3": 0, "4": 0}
    aligned = 597 + 800 + 510 + 540
    assert [round(prof[k] * aligned, 6) for k in ("sub_error", "ins_event", "del_event")] == [1, 1, 1]
    assert prof["kept_intron_per_intron"] == 1 / 4 and prof["full5_frac"] == 1.0
    assert prof["ins_len"] == {"2": 1.0} and prof["del_len"] == {"3": 1.0}
    assert prof["read_length"]["0.5"] == 599 and set(prof) == KEYS

def test_cap_and_seed(tmp_path):
    names = [f"r{i:02d}" for i in range(10)]
    p = world(tmp_path, [(q, "chr1", 2000, FULL) for q in names])
    got = lambda **kw: [rd["name"] for rd in sampled(p, **kw)]
    assert got() == names and got(cap=10, seed=7) == names                               # no more reads than the cap: all, BAM order
    assert got(cap=4, seed=7) == random.Random("7:G1").sample(names, 4) == got(cap=4, seed=7)
    assert len({tuple(got(cap=4, seed=s)) for s in range(5)}) > 1

def test_end_to_end(tmp_path):
    reads = [("ev", "chr1", 2000, "100=1X99=2I100=200N3D297="), ("ok", "chr1", 2000, FULL), ("minus", "chr1", 5060, "240=200N300=", 16)]
    reads += [(f"x{g}", c, t[0][0][0], FULL) for g, _, c, _, t, *_ in GENES if g in ("M1", "L1", "T2", "N2", "X1")] + [("xS1", "chr1", 12000, "600=")]
    p = world(tmp_path, reads)
    run = subprocess.run([sys.executable, os.path.join(HERE, "read_profile.py"), p["bam"], p["genes"], p["gff"], p["fa"], p["fam"],
                          p["out"], "--cap", "50", "--seed", "1"], capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    prof = json.load(open(p["out"]))
    assert set(prof) == KEYS and prof["source"] == p["bam"]
    assert (prof["genes"], prof["reads"], len(prof["rates"])) == (2, 3, 3)               # chrX is not in the default --chroms
    assert set(prof["end_pairs"]) == set(prof["class_reads"]) == {"0", "1", "2", "3", "4"} and prof["class_reads"]["0"] == 3
    assert "genes 2" in run.stdout and "reads 3" in run.stdout
    spec = importlib.util.spec_from_file_location("realworld", os.path.join(HERE, "..", "o3_graph", "realworld.py"))
    W = importlib.util.module_from_spec(spec); spec.loader.exec_module(W)
    seq, kinds = "".join(random.Random(3).choice("ACGT") for _ in range(6000)), set()
    for r in prof["rates"]:
        one = dict(sub_error=r[0], ins_event=r[1], del_event=r[2], **{k: prof[k] for k in ("ins_len", "del_len", "ins_hp_frac", "del_hp_frac")})
        kinds |= {k for _, k, _ in W.draw_edits(seq, random.Random(5), one)}
    assert kinds == {"S", "I", "D"}                                                      # the lengths and shares were consumed
