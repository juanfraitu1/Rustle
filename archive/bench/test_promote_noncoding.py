#!/usr/bin/env python3
"""Tests for the non-coding-aware promotion track (bench/promote_noncoding.py).
Hermetic: the module must import without pysam/pyabpoa/minimap2.
Run: PYTHONHASHSEED=0 pytest bench/test_promote_noncoding.py -v"""
import os
import sys
sys.path.insert(0, "bench")
import promote_noncoding as PN


def _rec(**kw):
    base = dict(genome_id=0.93, genome_cov=1.02, own_locus=True,
                alt_cols=82, alt_read_fraction=0.318, alt_reads=61, n_primary=192)
    base.update(kw)
    return base


def test_flagship_promotes():
    promote, call, reason = PN.classify_noncoding(_rec())
    assert promote is True
    assert call == "noncoding-candidate"
    assert reason  # non-empty


def test_artifact_low_identity_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(genome_id=0.28))
    assert promote is False and call == "artifact"


def test_near_reference_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(genome_id=0.99))
    assert promote is False and call == "~REF"


def test_chimeric_low_coverage_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(genome_cov=0.5))
    assert promote is False and call == "chimera"


def test_not_own_locus_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(own_locus=False))
    assert promote is False and call == "not-own-locus"


def test_thin_columns_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(alt_cols=2))
    assert promote is False and call == "thin"


def test_unbalanced_fraction_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(alt_read_fraction=0.05))
    assert promote is False and call == "thin"


def test_low_depth_rejected():
    promote, call, _ = PN.classify_noncoding(_rec(alt_reads=3))
    assert promote is False and call == "thin"


def test_boundaries_inclusive():
    # genome_id == GID_LO promotes; == GID_HI is ~REF (excluded)
    assert PN.classify_noncoding(_rec(genome_id=PN.GID_LO))[0] is True
    assert PN.classify_noncoding(_rec(genome_id=PN.GID_HI))[1] == "~REF"
    # genome_cov == COV_MIN promotes; alt_cols == MIN_COLS promotes
    assert PN.classify_noncoding(_rec(genome_cov=PN.COV_MIN))[0] is True
    assert PN.classify_noncoding(_rec(alt_cols=PN.MIN_COLS))[0] is True
    # alt_read_fraction at both band edges promotes
    assert PN.classify_noncoding(_rec(alt_read_fraction=PN.AF_LO))[0] is True
    assert PN.classify_noncoding(_rec(alt_read_fraction=PN.AF_HI))[0] is True


def test_coding_potential_label():
    assert PN.coding_potential(102, 3342) == "noncoding"            # flagship: 9.2%
    assert PN.coding_potential(1028, 4915) == "coding-capable/novel-protein"  # 62.7%
    assert PN.coding_potential(0, 0) == "noncoding"                 # guard div-by-zero


def test_orf_aa_finds_frame():
    # ATG GCA GCA TAA  -> M A A * -> ORF "MAA" length 3
    assert PN.orf_aa("ATGGCAGCATAA") == 3


def test_biotype_of_overlap():
    idx = PN.load_gff_from_lines([
        "NC_1\tX\tgene\t100\t500\t.\t+\t.\tgene_biotype=lncRNA",
        "NC_1\tX\tgene\t2000\t2500\t.\t+\t.\tgene_biotype=protein_coding",
    ])
    assert PN.biotype_of(idx, "NC_1", 150, 200) == "lncRNA"
    assert PN.biotype_of(idx, "NC_1", 2100, 2200) == "protein_coding-gene-body"
    assert PN.biotype_of(idx, "NC_1", 900, 950) == "intergenic"  # annotation loaded, no overlap
    assert PN.biotype_of({}, "NC_1", 150, 200) == "unknown"      # no annotation


def test_read_fasta_roundtrip(tmp_path):
    p = tmp_path / "c.fa"
    p.write_text(">a\nACGT\nACGT\n>b\nTTTT\n")
    d = PN.read_fasta(str(p))
    assert d == {"a": "ACGTACGT", "b": "TTTT"}


def test_build_record_honesty_rail():
    flag = dict(chrom="NC_073236.2", start=139051025, end=139225258,
                n_alt_positions=82, alt_read_fraction=0.3177, n_alt_reads=61, n_primary_reads=192)
    r = PN.build_record("NC_073236.2_139051025", flag, ("NC_073236.2", 139171134, 0.931, 1.02),
                        3342, 102, "lncRNA", None, "not-tested", "reason text")
    assert r["copy_vs_allele"] == "candidate-DNA-needed"
    assert r["status"] == "flagged-reference-divergent-candidate"
    assert r["track"] == "noncoding"
    assert r["divergence"] == 6.9 and r["genome_id"] == 0.931
    assert r["coding_potential"] == "noncoding" and r["biotype"] == "lncRNA"


def test_write_tsv_has_all_columns(tmp_path):
    flag = dict(chrom="NC_1", start=1, end=9, n_alt_positions=10, alt_read_fraction=0.3,
                n_alt_reads=6, n_primary_reads=20)
    r = PN.build_record("NC_1_1", flag, ("NC_1", 1, 0.9, 1.0), 300, 50, "intergenic", None, "not-tested", "x")
    p = tmp_path / "o.tsv"
    PN.write_tsv([r], str(p))
    header = p.read_text().splitlines()[0].split("\t")
    assert header[0] == "cid" and "copy_vs_allele" in header and "status" in header and len(header) == 21


def test_best_hits_from_paf_picks_best():
    # cid q1: two hits; the higher identity*coverage wins (line B)
    paf = "\n".join([
        # qname qlen qs qe strand tname tlen ts te nmatch blocklen mapq ...tags
        "q1\t1000\t0\t900\t+\tNC_1\t9\t500\t1400\t800\t900\t60\ttp:A:P",
        "q1\t1000\t0\t950\t+\tNC_2\t9\t10\t960\t930\t950\t40\ttp:A:S",
    ])
    hits = PN.best_hits_from_paf(paf)
    assert set(hits) == {"q1"}
    tname, pos, gid, gcov, mapq = hits["q1"]
    assert tname == "NC_2"          # 0.979*0.95 > 0.889*0.90
    assert pos == 10 and mapq == 40
    assert round(gid, 3) == 0.979 and round(gcov, 3) == 0.95


def test_best_hits_ignores_short_lines():
    assert PN.best_hits_from_paf("garbage\tline\n\n") == {}


def test_promote_all_flagship_and_exclusions():
    cons = {"C_flag": "ATG" + "GCA" * 200, "C_art": "AAAA", "C_ref": "CCCC"}
    flags = {
        "C_flag": dict(chrom="NC_1", start=1000, end=2000, n_alt_positions=40,
                       alt_read_fraction=0.32, n_alt_reads=30, n_primary_reads=90),
        "C_art": dict(chrom="NC_1", start=5000, end=6000, n_alt_positions=40,
                      alt_read_fraction=0.32, n_alt_reads=30, n_primary_reads=90),
        "C_ref": dict(chrom="NC_1", start=8000, end=9000, n_alt_positions=40,
                      alt_read_fraction=0.32, n_alt_reads=30, n_primary_reads=90),
    }
    hits = {
        "C_flag": ("NC_1", 1050, 0.90, 1.00, 60),   # divergent, own-locus -> promote
        "C_art": ("NC_1", 5050, 0.30, 1.00, 60),    # genome_id < GID_LO -> artifact
        "C_ref": ("NC_1", 8050, 0.99, 1.00, 60),    # genome_id >= GID_HI -> ~REF
    }
    gff = PN.load_gff_from_lines(["NC_1\tX\tgene\t900\t2100\t.\t+\t.\tgene_biotype=lncRNA"])
    promoted, tally = PN.promote_all(cons, flags, hits, gff, {}, {})
    assert [r["cid"] for r in promoted] == ["C_flag"]
    assert promoted[0]["biotype"] == "lncRNA"
    assert promoted[0]["copy_vs_allele"] == "candidate-DNA-needed"
    assert promoted[0]["protein"] == "not-tested"     # carried default
    assert tally["artifact"] == 1 and tally["~REF"] == 1 and tally["noncoding-candidate"] == 1


def test_promote_all_missing_flag_or_hit_tallied():
    cons = {"C_x": "ACGT"}
    promoted, tally = PN.promote_all(cons, {}, {}, {}, {}, {})
    assert promoted == [] and tally["no-flag"] == 1
    promoted, tally = PN.promote_all(cons, {"C_x": dict(chrom="NC_1", start=1, end=2,
        n_alt_positions=40, alt_read_fraction=0.3, n_alt_reads=30, n_primary_reads=90)}, {}, {}, {}, {})
    assert promoted == [] and tally["no-hit"] == 1


def test_promote_all_excludes_already_coding():
    cons = {"C_flag": "ATG" + "GCA" * 200}
    flags = {"C_flag": dict(chrom="NC_1", start=1000, end=2000, n_alt_positions=40,
                            alt_read_fraction=0.32, n_alt_reads=30, n_primary_reads=90)}
    hits = {"C_flag": ("NC_1", 1050, 0.90, 1.00, 60)}   # would promote if not excluded
    promoted, tally = PN.promote_all(cons, flags, hits, {}, {}, {}, exclude={"C_flag"})
    assert promoted == [] and tally["already-coding"] == 1


import json as _json

_CONS = "/home/juanfra/winloci_scratch/refabsent/gw_promoted/cons.fa"
_GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
_DATA_PRESENT = os.path.exists(_CONS) and os.path.exists(_GENOME) and os.path.exists(PN.MM2)

_STRONG_CREDIBLE = [
    "NC_073236.2_139051025",  # flagship lncRNA
    "NC_073236.2_44341131",   # lncRNA
    "NC_073234.2_27052843",   # ETV6 gene-body
    "NC_073230.2_167930813",  # VIPR2 gene-body
    "NC_073231.2_4861523",    # TDRP gene-body
    "NC_073242.2_18569043",   # intergenic, 1028-aa novel ORF
    "NC_073224.2_45327469",   # intergenic -> TRAF5
]
_ARTIFACTS = ["NC_073243.2_30120917", "NC_073231.2_39477379"]  # div ~70/74, genome_id ~0.29/0.26


import pytest


@pytest.mark.skipif(not _DATA_PRESENT, reason="cons.fa / GGO.fasta / minimap2 not present")
def test_integration_real_734(tmp_path):
    out = tmp_path / "gw_noncoding_copies.json"
    sys.argv = ["promote_noncoding", "--cons", _CONS, "--out", str(out), "--threads", "6"]
    PN.main()
    recs = _json.load(open(out))
    cids = {r["cid"] for r in recs}
    # flagship + all 7 strong-credible are promoted
    missing = [c for c in _STRONG_CREDIBLE if c not in cids]
    assert not missing, f"strong-credible missing from promotion: {missing}"
    # the 2 high-divergence repeat/chimera artifacts are excluded
    assert not (set(_ARTIFACTS) & cids), "artifact leaked into promotion"
    # sane yield band (looser bar than the workflow's ~19; report the count)
    print(f"\n[integration] promoted {len(recs)} non-coding candidates; "
          f"biotypes={sorted({r['biotype'] for r in recs})}")
    assert len(recs) >= 15
    # honesty rail holds on every record
    assert all(r["copy_vs_allele"] == "candidate-DNA-needed" for r in recs)
    assert all(r["status"] == "flagged-reference-divergent-candidate" for r in recs)
    assert all(r["track"] == "noncoding" for r in recs)
    # flagship specifics
    fs = next(r for r in recs if r["cid"] == "NC_073236.2_139051025")
    assert fs["coding_potential"] == "noncoding"
    if os.path.exists(PN.GFF):            # biotype needs the (linuxdisk-mounted) annotation
        assert fs["biotype"] == "lncRNA"
