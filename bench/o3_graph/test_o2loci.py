import importlib.util
import os

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_o2loci", os.path.join(HERE, "o2loci.py"))
L = importlib.util.module_from_spec(spec)
spec.loader.exec_module(L)

SPEC = {"e6": (100, 200), "e7": (300, 354), "e8": (500, 600)}


def test_a_read_spanning_exons_6_to_8_is_full_length_with_20_bp_of_exon_7_else_it_skips_exon_7():
    assert L.classify([(150, 200), (300, 354), (500, 550)], SPEC) == "FL"
    assert L.classify([(150, 200), (500, 550)], SPEC) == "D7"
    assert L.classify([(150, 200), (340, 354), (500, 550)], SPEC) == "D7"                    # 14 bp of exon 7
    assert L.classify([(300, 354), (500, 550)], SPEC) is None                               # does not reach exon 6
    assert L.inclusion({"FL": 3, "D7": 1}) == 0.75 and L.inclusion({}) is None


def test_the_chosen_copy_s_best_alignment_becomes_the_read_s_one_primary_with_its_sequence():
    h = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr5", "LN": 1000}]})
    prim = pysam.AlignedSegment(h)
    prim.query_name, prim.flag, prim.reference_id, prim.reference_start = "r1", 0, 0, 10
    prim.cigarstring, prim.query_sequence = "4M", "ACGG"
    prim.query_qualities = pysam.qualitystring_to_array("IIII")
    sec = pysam.AlignedSegment(h)
    sec.query_name, sec.flag, sec.reference_id, sec.reference_start = "r1", 256 | 16, 0, 500
    sec.cigarstring = "4M"
    recs = [dict(copy="F.0", AS=8, rec=prim), dict(copy="F.1", AS=7, rec=sec), dict(copy="F.1", AS=5, rec=sec)]
    assert L.choose(recs, "F.1")["AS"] == 7 and L.choose(recs, "F.2") is None
    out = L.as_primary(sec, prim)
    assert out.flag == 16 and out.query_sequence == "CCGT" and out.reference_start == 500
    assert pysam.qualities_to_qualitystring(out.query_qualities) == "IIII"


def test_assembled_transcripts_are_weighted_by_stringtie_cov_or_flair_s_mapped_reads(tmp_path):
    gtf = tmp_path / "st.gtf"
    gtf.write_text(
        'chr5\tStringTie\ttranscript\t101\t600\t.\t+\t.\tgene_id "g"; transcript_id "t1"; longcov "3.0"; cov "30.0";\n'
        'chr5\tStringTie\texon\t151\t200\t.\t+\t.\tgene_id "g"; transcript_id "t1";\n'
        'chr5\tStringTie\texon\t301\t354\t.\t+\t.\tgene_id "g"; transcript_id "t1";\n'
        'chr5\tStringTie\texon\t501\t550\t.\t+\t.\tgene_id "g"; transcript_id "t1";\n'
        'chr5\tStringTie\ttranscript\t101\t600\t.\t+\t.\tgene_id "g"; transcript_id "t2"; cov "10.0";\n'
        'chr5\tStringTie\texon\t151\t200\t.\t+\t.\tgene_id "g"; transcript_id "t2";\n'
        'chr5\tStringTie\texon\t501\t550\t.\t+\t.\tgene_id "g"; transcript_id "t2";\n')
    tx = L.read_gtf(str(gtf))
    assert tx["t1"]["blocks"] == [(150, 200), (300, 354), (500, 550)] and tx["t1"]["cov"] == 30.0
    (tmp_path / "map.txt").write_text("t1_chr5:100\tra,rb\nt2_chr5:100\trc\n")                 # FLAIR: isoform_gene
    assert L.read_map(str(tmp_path / "map.txt")) == {"t1": 2, "t2": 1}
    copies = {"F.1": dict(chrom="chr5", span=(0, 1000), **SPEC)}
    assert L.tx_inclusion(tx, {t: v["cov"] for t, v in tx.items()}, copies) == {"F.1": {"FL": 30.0, "D7": 10.0}}


def test_read_gtf_keeps_rustle_s_reads_attribute(tmp_path):
    gtf = tmp_path / "r.gtf"
    gtf.write_text('chr5\trustle\ttranscript\t101\t300\t.\t+\t.\tgene_id "g"; transcript_id "t1"; reads "18"; cov "7.7";\n'
                   'chr5\trustle\texon\t101\t150\t.\t+\t.\tgene_id "g"; transcript_id "t1"; exon_number "1";\n')
    tx = L.read_gtf(str(gtf))
    assert tx["t1"]["reads"] == 18.0 and tx["t1"]["cov"] == 7.7 and tx["t1"]["blocks"] == [(100, 150)]


def test_an_abstained_read_is_left_out_of_the_o2_arm_and_kept_at_its_primary_in_o2coin():
    calls = {"r1": "F.1", "r2": "", "r3": None}
    assert L.arm_copy("o2", calls, "r1", "F.0") == "F.1"
    assert L.arm_copy("o2", calls, "r2", "F.0") is None                   # abstained: no record
    assert L.arm_copy("o2", calls, "r3", "F.0") == "F.0"                  # not tied: minimap2's primary
    assert L.arm_copy("o2coin", calls, "r2", "F.0") == "F.0"
    assert L.arm_copy("mm2", calls, "r1", "F.0") == "F.0"
