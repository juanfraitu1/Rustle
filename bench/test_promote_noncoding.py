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
