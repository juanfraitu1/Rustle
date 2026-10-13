"""Unit test for the alignment layer's ortholog base read-through (needs minimap2)."""
import os, random, sys, tempfile
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "psv_sun"))
import pytest
import align_layer as A
import psv_sun as P


def test_ortholog_base_reads_the_mans_base_at_the_site():
    rng = random.Random(1)
    q = "".join(rng.choice("ACGT") for _ in range(30_000))      # locus +- 10 kb, site at 15,000
    man = list(rng.choice("ACGT") for _ in range(5_000)) + list(q) + list(rng.choice("ACGT") for _ in range(5_000))
    alt = "A" if q[15_000] != "A" else "C"
    man[5_000 + 15_000] = alt                                    # this man carries another base at the site
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/q.fa", "w").write(f">q\n{q}\n")
        open(f"{d}/m.fa", "w").write(">m\n" + "".join(man) + "\n")
        recs = P.run_paf(["minimap2", "-cx", "asm5", "--cs", f"{d}/m.fa", f"{d}/q.fa"])
        tseq = "".join(man)
    assert A.ortholog_base(recs, len(q), 10_000, 15_000, lambda name, s, e: tseq[s:e]) == alt


def test_no_single_flank_spanning_alignment_is_no_ortholog():
    assert A.ortholog_base([], 30_000, 10_000, 15_000, lambda n, s, e: "") is None


def test_reverse_strand_man_is_read_back_in_forward_reference_orientation():
    rng = random.Random(2)
    q = "".join(rng.choice("ACGT") for _ in range(30_000))
    man = list(rng.choice("ACGT") for _ in range(5_000)) + list(q) + list(rng.choice("ACGT") for _ in range(5_000))
    alt = "A" if q[15_000] != "A" else "C"
    man[5_000 + 15_000] = alt
    tseq = "".join(man)
    rc = tseq[::-1].translate(str.maketrans("ACGT", "TGCA"))     # the man as he sits on the reverse strand
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/q.fa", "w").write(f">q\n{q}\n")
        open(f"{d}/m.fa", "w").write(">m\n" + rc + "\n")
        recs = P.run_paf(["minimap2", "-cx", "asm5", "--cs", f"{d}/m.fa", f"{d}/q.fa"])
    assert recs and recs[0]["strand"] == "-"
    assert A.ortholog_base(recs, len(q), 10_000, 15_000, lambda n, s, e: rc[s:e]) == alt


def test_find_men_refuses_an_empty_assembly_dir(tmp_path):
    import kmer_layer as KL
    with pytest.raises(SystemExit):
        KL.find_men(str(tmp_path))


def test_a_man_with_two_copies_of_the_locus_has_no_unique_ortholog():
    # R8: every record spanning both flanks counts, primary or secondary. Two copies that differ at the site
    # make the base ambiguous, so the layer must abstain rather than pick one.
    rng = random.Random(3)
    q = "".join(rng.choice("ACGT") for _ in range(30_000))
    q2 = list(q)
    q2[15_000] = "A" if q[15_000] != "A" else "C"
    spacer = "".join(rng.choice("ACGT") for _ in range(5_000))
    tseq = q + spacer + "".join(q2)
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/q.fa", "w").write(f">q\n{q}\n")
        open(f"{d}/m.fa", "w").write(">m\n" + tseq + "\n")
        recs = P.run_paf(["minimap2", "-cx", "asm5", "--cs", "-N", "50", f"{d}/m.fa", f"{d}/q.fa"])
    spanning = [r for r in recs if r["qs"] < 10_000 and r["qe"] > 30_000 - 10_000]
    assert len(spanning) == 2, "minimap2 must report both copies (primary + secondary) for this test to bite"
    assert A.ortholog_base(recs, len(q), 10_000, 15_000, lambda n, s, e: tseq[s:e]) is None


def test_locus_window_refuses_sites_too_close_to_the_contig_start_or_end():
    # Site at p < FLANK: the window would be clamped at 0, shifting the site offset. Such a site must not qualify.
    seq = "ACGT" * 10_000                                   # 40 kb contig
    fetch = lambda c, s, e: seq[max(0, s):min(len(seq), e)]
    length = lambda c: len(seq)
    assert A.locus_window(fetch, length, "chrY", 5_000) is None           # p - FLANK < 0
    assert A.locus_window(fetch, length, "chrY", 35_000) is None          # p + FLANK + 1 > 40 kb
    w = A.locus_window(fetch, length, "chrY", 20_000)
    assert w == (seq[10_000:30_001], 10_000)                               # full window, site at offset FLANK


def test_site_at_contig_start_is_not_called_polymorphic_when_the_man_is_identical():
    # Regression for the clamp bug: an identical man must never yield a mismatch at a window-clamped site.
    seq = "".join(random.Random(4).choice("ACGT") for _ in range(40_000))
    fetch = lambda c, s, e: seq[max(0, s):min(len(seq), e)]
    assert A.locus_window(fetch, lambda c: len(seq), "chrY", 3_000) is None
