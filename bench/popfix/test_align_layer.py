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


# ---- resumable, cached alignment layer (task 4 operational fix) ---------------------------------

KSCAN = os.environ.get("POPFIX_KSCAN", "/mnt/linuxdisk/tmp/psv_sun/bin/kscan")
needs_kscan = pytest.mark.skipif(not os.path.exists(KSCAN), reason="compiled kscan not available")

PAF_LINE = ("q1\t20001\t0\t20001\t+\tchrY\t80000\t5000\t25001\t20000\t20001\t60\ttp:A:P\tcm:i:100\ts1:i:19000\t"
            "AS:i:39000\tcs:Z::10*ac:19989\n")


def test_parse_paf_reads_the_same_records_run_paf_returns_from_minimap2_stdout():
    (r,) = P.parse_paf(PAF_LINE)
    assert (r["qname"], r["qlen"], r["qs"], r["qe"], r["strand"]) == ("q1", 20001, 0, 20001, "+")
    assert (r["tname"], r["tlen"], r["ts"], r["te"], r["tp"], r["AS"]) == ("chrY", 80000, 5000, 25001, "P", 39000)
    assert r["cs"] == ":10*ac:19989"
    assert P.parse_paf("") == []


def test_cached_records_reads_a_complete_paf_and_ignores_a_leftover_part_file(tmp_path):
    cache, fa = str(tmp_path), "/some/where/man7.chrY.fasta"
    assert A.cached_records(cache, fa) is None
    (tmp_path / "man7.chrY.fasta.paf.part").write_text("garbage that is not a PAF record\n")
    assert A.cached_records(cache, fa) is None                       # a .part is never a cache hit
    (tmp_path / "man7.chrY.fasta.paf").write_text(PAF_LINE)
    assert A.cached_records(cache, fa) == P.parse_paf(PAF_LINE)
    (tmp_path / "man7.chrY.fasta.paf").write_text("")                # no alignments at all is a valid complete cache
    assert A.cached_records(cache, fa) == []


def _fasta(path, name, seq):
    with open(path, "w") as fh:
        fh.write(f">{name}\n{seq}\n")
    with open(path + ".fai", "w") as fh:
        fh.write(f"{name}\t{len(seq)}\t{len(name) + 2}\t{len(seq)}\t{len(seq) + 1}\n")


def _rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def _setup(tmp_path):
    """80 kb chrY, three PSV columns (one site each, all >= 10 kb from both ends), three men: man1 identical, man2
    carries another base at site 40,000, man3 is man1 on the reverse strand. -> main() args and the site list."""
    rng = random.Random(7)
    g = "".join(rng.choice("ACGT") for _ in range(80_000))
    sites = [12_000, 40_000, 68_000]
    gfa = str(tmp_path / "g.fa")
    _fasta(gfa, "chrY", g)
    asm = tmp_path / "asm"
    asm.mkdir()
    m2 = list(g)
    m2[40_000] = "A" if g[40_000] != "A" else "C"
    _fasta(str(asm / "man1.chrY.fasta"), "chrY", g)
    _fasta(str(asm / "man2.chrY.fasta"), "chrY", "".join(m2))
    _fasta(str(asm / "man3.chrY.fasta"), "chrY", _rc(g))
    cat = tmp_path / "copies.tsv"
    cat.write_text("tid\tchrom\nt0\tchrY\n")
    prefix = str(tmp_path / "p")
    (tmp_path / "p.psv_copies.tsv").write_text("family_id\tcopy_index\tcopy_tid\talleles\tstrand\nF1\t0\tt0\tAC\t+\n")
    (tmp_path / "p.psv_pos.tsv").write_text(
        "family_id\tcol_index\tcopy_index\tgenome_pos\tbase\n"
        + "".join(f"F1\t{i}\t0\t{p}\t{g[p]}\n" for i, p in enumerate(sites)))
    return prefix, str(cat), gfa, str(asm)


def _run(args, out):
    A.main(*args, str(out))
    return (out / "align_status.tsv").read_bytes()


def _paf_files(out):
    (cache,) = [p for p in out.iterdir() if p.is_dir() and p.name.startswith("paf_")]
    return cache, sorted(f for f in os.listdir(cache) if f.endswith(".paf"))


@needs_kscan
def test_main_runs_in_resumable_batches_and_matches_a_fresh_run(tmp_path, monkeypatch, capsys):
    prefix, cat, gfa, asm = _setup(tmp_path)
    args = (prefix, cat, gfa, asm, KSCAN)
    fresh = _run(args, tmp_path / "fresh")
    rows = [l.split("\t") for l in fresh.decode().splitlines()[1:]]
    assert [r[2] for r in rows] == ["FIXED", "POLYMORPHIC", "FIXED"], "the synthetic setup must exercise the verdicts"

    out = tmp_path / "batched"
    out.mkdir()
    (out / "align_status.tsv").write_text("stale: from an earlier, different run\n")
    monkeypatch.setenv("POPFIX_MAX_NEW_PAF", "1")
    with pytest.raises(SystemExit) as e:
        A.main(*args, str(out))
    assert e.value.code == 3
    assert "incomplete: 1/3 men cached; re-run" in capsys.readouterr().err
    cache, pafs = _paf_files(out)
    assert pafs == ["man1.chrY.fasta.paf"]                            # exactly one man computed, men in sorted order
    assert not any(f.endswith(".part") for f in os.listdir(cache))
    assert not (out / "align_status.tsv").exists()                    # neither stale nor partial

    monkeypatch.delenv("POPFIX_MAX_NEW_PAF")
    assert _run(args, out) == fresh
    assert _paf_files(out)[1] == ["man1.chrY.fasta.paf", "man2.chrY.fasta.paf", "man3.chrY.fasta.paf"]


@needs_kscan
def test_a_batch_that_computes_the_last_missing_man_finishes_instead_of_exiting_3(tmp_path, monkeypatch):
    prefix, cat, gfa, asm = _setup(tmp_path)
    args = (prefix, cat, gfa, asm, KSCAN)
    fresh = _run(args, tmp_path / "fresh")
    out = tmp_path / "o"
    monkeypatch.setenv("POPFIX_MAX_NEW_PAF", "1")
    for _ in range(2):
        with pytest.raises(SystemExit) as e:
            A.main(*args, str(out))
        assert e.value.code == 3
    assert len(_paf_files(out)[1]) == 2
    assert _run(args, out) == fresh            # one man left, budget of one: computing it leaves nobody missing


@needs_kscan
def test_a_fully_cached_run_never_calls_minimap2(tmp_path, monkeypatch):
    prefix, cat, gfa, asm = _setup(tmp_path)
    args = (prefix, cat, gfa, asm, KSCAN)
    first = _run(args, tmp_path / "o")

    def boom(*a, **k):
        raise AssertionError("minimap2 must not run when every man is cached")
    monkeypatch.setattr(A, "compute_paf", boom)
    assert _run(args, tmp_path / "o") == first


@needs_kscan
def test_a_leftover_part_file_is_recomputed_not_read(tmp_path):
    prefix, cat, gfa, asm = _setup(tmp_path)
    args = (prefix, cat, gfa, asm, KSCAN)
    fresh = _run(args, tmp_path / "fresh")
    out = tmp_path / "o"
    _run(args, out)
    cache, _ = _paf_files(out)
    os.remove(os.path.join(cache, "man2.chrY.fasta.paf"))             # as if the run died while computing man2 ...
    with open(os.path.join(cache, "man2.chrY.fasta.paf.part"), "w") as fh:
        fh.write("q\t1\t0\t1\t+\tchrY\t1\t0\t1\t1\t1\t60\tcs:Z::1\n")  # ... leaving a parsable but wrong fragment
    assert _run(args, out) == fresh
    assert not os.path.exists(os.path.join(cache, "man2.chrY.fasta.paf.part"))
    assert os.path.exists(os.path.join(cache, "man2.chrY.fasta.paf"))


def test_compute_paf_passes_the_thread_count_and_a_failed_minimap2_leaves_no_paf(tmp_path, monkeypatch):
    import subprocess
    seen = []

    def fake_ok(cmd, stdout=None, **kw):
        seen.append(cmd)
        stdout.write(PAF_LINE.encode())
        return subprocess.CompletedProcess(cmd, 0)
    monkeypatch.setattr(A.subprocess, "run", fake_ok)
    A.compute_paf(str(tmp_path), "/x/man1.chrY.fasta", "/x/loci.fa", 5)
    assert seen[0][seen[0].index("-t") + 1] == "5" and seen[0][-2:] == ["/x/man1.chrY.fasta", "/x/loci.fa"]
    assert (tmp_path / "man1.chrY.fasta.paf").read_text() == PAF_LINE
    assert not (tmp_path / "man1.chrY.fasta.paf.part").exists()

    def fake_crash(cmd, stdout=None, **kw):
        stdout.write(PAF_LINE.encode()[:30])                         # a truncated record, then minimap2 dies
        raise subprocess.CalledProcessError(137, cmd)
    monkeypatch.setattr(A.subprocess, "run", fake_crash)
    with pytest.raises(subprocess.CalledProcessError):
        A.compute_paf(str(tmp_path), "/x/man2.chrY.fasta", "/x/loci.fa", 2)
    assert not (tmp_path / "man2.chrY.fasta.paf").exists()           # only a .part remains, and that is never read
    assert A.cached_records(str(tmp_path), "/x/man2.chrY.fasta") is None
