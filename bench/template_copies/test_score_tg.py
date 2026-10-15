import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import score_tg as S  # noqa: E402


def lrow(src, cls, s, e, short=False, seqid=1.0):
    return {"source_id": src, "cls": cls, "contig": "chr1", "start0": str(s), "end": str(e), "iv": [(s, e)],
            "short": short, "seqid_f": seqid}


def test_liftoff_scores_by_hand():
    loci = [lrow("S1", "in_place", 0, 1000), lrow("S1", "extra_copy", 5000, 6000), lrow("S1", "extra_copy", 9000, 10000),
            lrow("S2", "in_place", 20000, 21000), lrow("S2", "extra_copy", 30000, 31000, seqid=0.9)]
    sup = {S.LF.locus_key(r): 3 for r in loci}
    base = [{"contig": "chr1", "locus": "b1", "iv": [(0, 1000)]}, {"contig": "chr1", "locus": "b2", "iv": [(40000, 41000)]}]
    tg = [{"contig": "chr1", "locus": "t1", "iv": [(5000, 6000)]}]
    res = S.liftoff_scores(loci, sup, base, base + tg, tg, {"b1": 5, "b2": 9})
    # copies at >= 0.95: S1's two; base covers neither; S1's placement is assembled, so both are reachable; one reached
    assert res["recall"] == {"copies": 2, "base": 0, "merged": 1, "reachable": 2, "reached": 1}
    assert res["p1_located"] == {"base": 0.5, "merged": 2 / 3, "passes": True}
    assert res["p2"]["tg_located"] == 1.0 and res["p2"]["base_low_loci"] == 1 and res["p2"]["base_low_located"] == 1.0
    assert res["p2"]["passes"] is True


TABLE = """[human] genome-wide minus chr20 (sensitivity / precision)
arm                     base        intron         chain    transcript         locus
OP_primary_TG   50.6 /  18.2  66.7 /  56.4  23.9 /  15.5  22.0 /  14.0  36.1 /  23.5
StringTie       59.6 /  14.4  63.5 /  57.8  18.1 /  15.0  16.9 /  12.7  33.7 /  20.1
FLAIR           66.4 /   8.8  67.7 /  30.1  22.4 /   6.2  21.1 /   3.6  31.9 /   6.3
default_primary  48.1 /  17.4  60.8 /  64.4  21.8 /  18.5  20.3 /  16.1  34.5 /  22.8
stratum S1: 926 multi-exon reference transcripts in 304 genes
   OP_primary_TG  matched     0 / 926   precision 0 / 0 = 0.000
"""


def test_parse_score_op():
    arms = S.parse_score_op(TABLE)
    assert arms == {"OP_primary_TG": (23.9, 15.5), "StringTie": (18.1, 15.0), "FLAIR": (22.4, 6.2),
                    "default_primary": (21.8, 18.5)}


def test_parse_stats(tmp_path):
    p = tmp_path / "x.stats"
    p.write_text("#= Summary\n        Intron level:    66.7     |    56.4    |\nIntron chain level:    23.6     |    15.6    |\n")
    assert S.parse_stats(str(p)) == {"intron_sn": 66.7, "intron_pr": 56.4, "intron_chain_sn": 23.6, "intron_chain_pr": 15.6}


def test_interval_is_seeded():
    assert S.interval([1, 0, 1, 1]) == S.interval([1, 0, 1, 1]) and S.interval([]) != S.interval([])  # NaN != NaN


def test_lead():
    arms = S.parse_score_op(TABLE)
    assert S.lead(arms, "OP_primary_TG") == {"b1": True, "b2": True}     # 23.9 > 18.1 and > 22.4; 15.5 > 6.2
    assert S.lead({"OP_primary_TG": (23.9, 15.5)}, "OP_primary_TG") == {"b1": None, "b2": None}
    assert S.lead(dict(arms, FLAIR=(24.0, 6.2)), "OP_primary_TG")["b1"] is False


def test_chains_reads_a_gzipped_reference(tmp_path):
    import gzip
    gtf = ('chr1\tt\ttranscript\t101\t400\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
           'chr1\tt\texon\t101\t200\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
           'chr1\tt\texon\t301\t400\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
    ref = tmp_path / "ref.gtf.gz"
    with gzip.open(ref, "wt") as fh:
        fh.write(gtf)
    q = tmp_path / "q.gtf"
    q.write_text(gtf)
    assert S.chains(str(ref), str(tmp_path / "work"), [("q", str(q))]) == {"q": {"chain_sn": 100.0, "chain_pr": 100.0}}
