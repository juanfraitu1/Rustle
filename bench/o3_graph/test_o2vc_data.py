import array
import gzip
import importlib.util
import os
import random
import subprocess
import types

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    s = importlib.util.spec_from_file_location(f"o3_graph_{name}_t", os.path.join(HERE, f"{name}.py"))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


D, R = _load("o2vc_data"), _load("o2vc_rule")
VCF = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
       "chr1\t101\t.\tA\tG\t30\tPASS\t.\tGT:GQ:PS:DP:AF:PQ\t0|1:30:101:12:0.50:20\n"
       "chr1\t201\t.\tC\tT,G\t20\tMultiallelic\t.\tGT:GQ:DP:AF\t1/2:20:10:0.40,0.60\n"
       "chr1\t301\t.\tT\tC\t5\tLowQual\t.\tGT:GQ:DP:AF\t0/1:5:8:0.30\n")


def test_parse_vcf_reads_every_record_0_based_with_gt_filter_and_depth(tmp_path):
    for p in (tmp_path / "x.vcf", tmp_path / "x.vcf.gz"):
        (gzip.open if p.suffix == ".gz" else open)(p, "wt").write(VCF)
        v = D.parse_vcf(str(p))
        assert v[("chr1", 100)] == R.Call("A", ("G",), (0, 1), "PASS", 12)
        assert v[("chr1", 200)] == R.Call("C", ("T", "G"), (1, 2), "Multiallelic", 10)
        assert v[("chr1", 300)] == R.Call("T", ("C",), (0, 1), "LowQual", 8)


def test_longcallr_runs_with_the_masseq_preset_min_mapq_0_and_no_bam(monkeypatch, tmp_path):
    seen = {}

    def fake(cmd, **kw):
        seen["cmd"] = cmd
        open(str(tmp_path / "o") + ".vcf", "w").write("#CHROM\n")       # longcallR 2.0.0 writes OUT_PREFIX.vcf, plain text (observed)
        return subprocess.CompletedProcess(cmd, 0, "", "")
    monkeypatch.setattr(D.subprocess, "run", fake)
    out = D.run_longcallr("t.bam", "g.fa", str(tmp_path / "o"), binary="LC", threads=2, region="chr1:1-100")
    assert seen["cmd"] == ["LC", "-b", "t.bam", "-f", "g.fa", "-p", "hifi-masseq", "--min-mapq", "0", "--no-bam-output",
                           "-t", "2", "-o", str(tmp_path / "o"), "-r", "chr1:1-100"]
    assert out == str(tmp_path / "o") + ".vcf"


def _fake_run(seen, rc=0, err="", out="", vcf=None):
    """A stand-in for subprocess.run: keeps the command, writes the VCF when given one, answers rc, out and err."""
    def fake(cmd, **kw):
        seen["cmd"] = cmd
        if vcf:
            open(vcf, "w").write("#CHROM\n")
        return subprocess.CompletedProcess(cmd, rc, out, err)
    return fake


def test_a_record_with_a_gt_that_calls_no_allele_is_kept_with_an_empty_gt(tmp_path):
    gts = ["./.", ".|.", ".", "./1"]                    # missing, phased missing, haploid missing, and one allele missing
    p = tmp_path / "m.vcf"
    p.write_text(VCF + "".join(f"chr1\t{401 + i}\t.\tG\tA\t3\tLowQual\t.\tGT:GQ:DP:AF\t{gt}:3:7:0.10\n" for i, gt in enumerate(gts)))
    v = D.parse_vcf(str(p))
    assert [v[("chr1", 400 + i)] for i in range(4)] == [R.Call("G", ("A",), (), "LowQual", 7)] * 4
    assert v[("chr1", 100)].gt == (0, 1)                # the called records are as before


def test_longcallr_gets_no_r_without_a_region(monkeypatch, tmp_path):
    seen = {}
    monkeypatch.setattr(D.subprocess, "run", _fake_run(seen, vcf=str(tmp_path / "o") + ".vcf"))
    assert D.run_longcallr("t.bam", "g.fa", str(tmp_path / "o"), binary="LC") == str(tmp_path / "o") + ".vcf"
    assert seen["cmd"] == ["LC", "-b", "t.bam", "-f", "g.fa", "-p", "hifi-masseq", "--min-mapq", "0", "--no-bam-output",
                           "-t", "2", "-o", str(tmp_path / "o")]


def test_a_failing_longcallr_raises_with_the_exit_code_and_the_last_twenty_stderr_lines(monkeypatch, tmp_path):
    err = "".join(f"line {i}\n" for i in range(30))
    monkeypatch.setattr(D.subprocess, "run", _fake_run({}, rc=101, err=err))
    with pytest.raises(RuntimeError) as e:
        D.run_longcallr("t.bam", "g.fa", str(tmp_path / "o"), binary="LC", region="chr1:0-100")
    msg = str(e.value)
    assert "101" in msg and "line 10" in msg and "line 29" in msg and "line 9" not in msg


def test_longcallr_exiting_0_without_its_vcf_is_an_error_not_a_missing_path(monkeypatch, tmp_path):
    monkeypatch.setattr(D.subprocess, "run", _fake_run({}))
    with pytest.raises(RuntimeError, match="wrote no VCF"):
        D.run_longcallr("t.bam", "g.fa", str(tmp_path / "o"), binary="LC")


def test_the_caller_s_version_is_the_line_its_version_flag_prints(monkeypatch):
    seen = {}
    monkeypatch.setattr(D.subprocess, "run", _fake_run(seen, out="longcallR 2.0.0\n"))
    assert D.longcallr_version("LC") == "longcallR 2.0.0"
    assert seen["cmd"] == ["LC", "--version"]


def test_longcallr_gets_the_directory_of_its_output_prefix_made_for_it_panics_without_one(monkeypatch, tmp_path):
    prefix = str(tmp_path / "vcf" / "F.0")

    def fake(cmd, **kw):
        assert os.path.isdir(tmp_path / "vcf")                                  # longcallR 2.0.0 does not make it (exit 101)
        open(prefix + ".vcf", "w").write("#CHROM\n")
        return subprocess.CompletedProcess(cmd, 0, "", "")
    monkeypatch.setattr(D.subprocess, "run", fake)
    assert D.run_longcallr("t.bam", "g.fa", prefix, binary="LC") == prefix + ".vcf"


# ---- Task 3: the BAM side (tied domain, teachers, read-frame columns, teacher BAM, depth) ----

def write_bam(path, refs, recs, copies_tsv):
    """A small world on disk in the directory path (made if missing): refs ({chrom: seq}) as genome.fa with its .fai,
    copies_tsv (rows (copy, chrom, start, end, strand); the family is the copy's name before its last dot) as a pool-format
    copies.tsv whose exons are the span, and recs as reads.bam, sorted and indexed. recs items: dict(name, chrom, start,
    cigar (a string), seq, qual, secondary, reverse, AS) with optional NM (0), de (0.0; None leaves the tag out), ms (AS),
    mapq (60, or 0 on a secondary) and flags (more FLAG bits, e.g. 0x800). seq is in the genome's strand, as a BAM holds it;
    seq None leaves SEQ unset (a secondary, as minimap2 writes one) and qual None leaves QUAL unset. Returns
    SimpleNamespace(bam, copies, genome) of paths."""
    import pysam
    os.makedirs(path, exist_ok=True)
    genome, copies, bam = (os.path.join(str(path), x) for x in ("genome.fa", "copies.tsv", "reads.bam"))
    with open(genome, "w") as fh:
        fh.write("".join(f">{ch}\n{seq}\n" for ch, seq in refs.items()))
    pysam.faidx(genome)
    with open(copies, "w") as fh:
        fh.write("copy\tfamily\tchrom\tstart\tend\tstrand\texons\n")
        for c, ch, s, e, strand in copies_tsv:
            fh.write(f"{c}\t{c.rsplit('.', 1)[0]}\t{ch}\t{s}\t{e}\t{strand}\t{s}-{e}\n")
    hdr = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": ch, "LN": len(s)} for ch, s in refs.items()]})
    tid = {ch: i for i, ch in enumerate(refs)}
    with pysam.AlignmentFile(bam + ".unsorted.bam", "wb", header=hdr) as out:
        for r in recs:
            a = pysam.AlignedSegment(hdr)
            a.query_name, a.reference_id, a.reference_start = r["name"], tid[r["chrom"]], r["start"]
            a.flag = (0x100 if r["secondary"] else 0) | (0x10 if r["reverse"] else 0) | r.get("flags", 0)
            a.mapping_quality = r.get("mapq", 0 if r["secondary"] else 60)
            a.cigarstring = r["cigar"]
            if r["seq"] is not None:
                a.query_sequence = r["seq"]
                if r["qual"] is not None:
                    a.query_qualities = array.array("B", r["qual"])
            a.set_tag("AS", r["AS"])
            a.set_tag("NM", r.get("NM", 0))
            if r.get("de", 0.0) is not None:
                a.set_tag("de", float(r.get("de", 0.0)), value_type="f")
            a.set_tag("ms", r.get("ms", r["AS"]))
            out.write(a)
    pysam.sort("-o", bam, bam + ".unsorted.bam")
    pysam.index(bam)
    os.remove(bam + ".unsorted.bam")
    return types.SimpleNamespace(bam=bam, copies=copies, genome=genome)


def _random_refs(seed, chroms=("chrA", "chrB"), n=2000):
    rng = random.Random(seed)
    return {ch: "".join(rng.choice("ACGT") for _ in range(n)) for ch in chroms}


def _rec(name, chrom, start, AS, cigar="600M", seq=None, qual=None, secondary=False, reverse=False, **kw):
    return dict(name=name, chrom=chrom, start=start, cigar=cigar, seq=seq, qual=qual, secondary=secondary, reverse=reverse,
                AS=AS, **kw)


def _no_runs(seed, chroms=("chrA", "chrB"), n=2000):
    """Random bases with no two equal neighbours, so no read of them holds a poly-run window, one sequence per chrom."""
    rng, out = random.Random(seed), {}
    for ch in chroms:
        s = []
        for _ in range(n):
            s.append(rng.choice([b for b in "ACGT" if not s or b != s[-1]]))
        out[ch] = "".join(s)
    return out


def _with_alts(seq, positions):
    """seq with an ALT at each position: the first base that is neither its own nor a neighbour's (it makes no run)."""
    s = list(seq)
    for p in positions:
        s[p] = next(b for b in "ACGT" if b not in (seq[p - 1], seq[p], seq[p + 1]))
    return "".join(s)


def make_world(tmp_path, alts=()):
    """The brief's world: chrA and chrB (2,000 random bases, seed 1); F.0 = chrA:100-700, F.1 = chrB:100-700, F.2 =
    chrB:1000-1600, G.0 = chrA:1000-1600, all +. Reads of 600 bp, 600M; the primaries carry F.0's reference bases, with an
    ALT at each of `alts` (chrA positions) when given. T1-T6: teachers of F.0 (T6 has base quality 5 at chrA:400); Y: F.0
    600 / F.1 590, a teacher (margin 10); W: F.0 600 / G.0 598, a teacher (off-family alignments do not count); X: F.0 600 /
    F.1 595, tied; Z: F.0 500 / F.1 595 / F.2 590 (reverse), tied with its primary outside the tied set; S: a 400-bp teacher
    (under the caller's read-length floor); H: a teacher with de 0.08. Returns (BAM, copies TSV); the genome is genome.fa
    beside the BAM."""
    refs = _random_refs(1)
    f0 = _with_alts(refs["chrA"], alts)[100:700]
    copies = [("F.0", "chrA", 100, 700, "+"), ("F.1", "chrB", 100, 700, "+"), ("F.2", "chrB", 1000, 1600, "+"),
              ("G.0", "chrA", 1000, 1600, "+")]

    def prim(name, AS=600, qual=None, **kw):
        return _rec(name, "chrA", 100, AS, seq=f0, qual=qual or [30] * 600, **kw)

    low = [30] * 600
    low[400 - 100] = 5                                                  # T6: base quality 5 at chrA:400
    recs = [prim(f"T{i}", qual=low if i == 6 else None) for i in range(1, 7)]
    recs += [prim("Y"), _rec("Y", "chrB", 100, 590, secondary=True)]
    recs += [prim("W"), _rec("W", "chrA", 1000, 598, secondary=True)]
    recs += [prim("X"), _rec("X", "chrB", 100, 595, secondary=True)]
    recs += [prim("Z", AS=500), _rec("Z", "chrB", 100, 595, secondary=True),
             _rec("Z", "chrB", 1000, 590, secondary=True, reverse=True)]
    recs += [_rec("S", "chrA", 100, 400, cigar="400M", seq=f0[:400], qual=[30] * 400), prim("H", de=0.08)]
    w = write_bam(tmp_path / "world", refs, recs, copies)
    return w.bam, w.copies


def test_tied_wide_keeps_an_as_tie_whose_primary_lies_below_it():
    alns = [dict(AS=1734, prim=True, copy="F.2"), dict(AS=1819, prim=False, copy="F.5"), dict(AS=1819, prim=False, copy="F.7")]
    per, prim = D.tied_wide(alns)
    assert prim == "F.2" and sorted(per) == ["F.5", "F.7"]
    assert D.tied_wide(alns + [dict(AS=1815, prim=False, copy=None)]) == (None, "F.2")
    assert D.tied_wide([dict(AS=900, prim=True, copy="F.0"), dict(AS=890, prim=False, copy="F.1")]) == (None, "F.0")
    assert D.tied_wide([dict(AS=900, prim=True, copy="F.0"), dict(AS=899, prim=False, copy="G.1")]) == (None, "F.0")


def test_a_teacher_is_a_primary_at_least_ten_as_above_every_other_copy_of_its_family():
    assert D.is_teacher(1000, [990]) and not D.is_teacher(1000, [991]) and D.is_teacher(1000, [])


def test_read_frame_positions_count_both_clip_kinds_skip_insertions_and_introns_and_flip_reverse_alignments():
    assert D.read_frame_positions(100, [(4, 2), (0, 3), (1, 1), (0, 2), (5, 2)], 10, False) == {2: 100, 3: 101, 4: 102, 6: 103, 7: 104}
    assert D.read_frame_positions(100, [(5, 2), (0, 3), (1, 1), (0, 2), (4, 2)], 10, True) == {7: 100, 6: 101, 5: 102, 3: 103, 2: 104}
    assert D.read_frame_positions(100, [(0, 2), (3, 50), (2, 1), (0, 2)], 4, False) == {0: 100, 1: 101, 2: 153, 3: 154}


def test_columns_are_the_read_bases_every_tied_alignment_aligns_with_each_copy_s_base_in_the_read_frame():
    ref = {"chrA": "AAAACGTACG", "chrB": "CGTGCGTTTT"}
    fetch = lambda ch, lo, hi: ref[ch][lo:hi]
    a0 = D.Aln("F.0", "chrA", 4, ((0, 5),), False, 5, 0, 0.0, 5)
    a1 = D.Aln("F.1", "chrB", 1, ((0, 5),), True, 3, 1, 0.2, 3)
    read = D.TiedRead("F", {"F.0": a0, "F.1": a1}, "F.0", 0, 2, "CGTAC", [30] * 5)
    cols = D.columns(read, fetch)
    assert [c["p"] for c in cols] == [0, 1, 2, 3, 4]
    assert cols[2] == dict(p=2, base="T", qual=30, at={"F.0": ("chrA", 6), "F.1": ("chrB", 3)}, ref={"F.0": "T", "F.1": "C"})
    ins = D.TiedRead("F", {"F.0": D.Aln("F.0", "chrA", 4, ((0, 2), (1, 1), (0, 2)), False, 3, 1, 0.2, 3), "F.1": a1},
                     "F.0", 0, 0, "CGTAC", [30] * 5)
    assert [c["p"] for c in D.columns(ins, fetch)] == [0, 1, 3, 4]


def test_scan_finds_the_tied_reads_and_the_teachers(tmp_path):
    bam, copies = make_world(tmp_path)                               # the world above
    s = D.scan(bam, copies)
    assert sorted(s.tied) == ["X", "Z"]
    assert s.tied["Z"].prim == "F.0" and sorted(s.tied["Z"].copies) == ["F.1", "F.2"] and s.tied["Z"].copies["F.2"].flip
    assert s.tied["X"].copies["F.1"].AS == 595 and s.tied["X"].margin == 5 and len(s.tied["X"].seq) == 600
    assert s.teachers["F.0"] == {"T1", "T2", "T3", "T4", "T5", "T6", "Y", "W", "S", "H"}


def test_the_teacher_bam_holds_only_the_teachers_primaries_and_depth_counts_reads_as_the_caller_does(tmp_path):
    import pysam
    bam, copies = make_world(tmp_path)
    s = D.scan(bam, copies)
    out = str(tmp_path / "teach.bam")
    assert D.write_teacher_bam(bam, {"F.0": ("chrA", 100, 700)}, s.teachers, out) == 10
    # the caller's DP (longcallR 2.0.0, hifi-masseq): T1-T6, Y, W and H; S is out (400 bp < 500), while H (de 0.08 < 0.5),
    # T6 (base quality 5: any quality counts) and position 110 (read base 10: HiFi trims no read end) are in
    fetch = pysam.FastaFile(os.path.join(os.path.dirname(bam), "genome.fa")).fetch      # required (Task 4 ruling 3)
    assert D.teacher_depth(out, "chrA", [110, 400], fetch) == {110: 9, 400: 9}


# ---- beyond the brief's tests: behaviour the brief states but does not pin ----

def test_scan_puts_a_reverse_primary_in_the_read_frame_and_skips_a_read_whose_cigar_and_seq_lengths_disagree(tmp_path, capsys):
    refs = _random_refs(2)
    copies = [("F.0", "chrA", 100, 700, "+"), ("F.1", "chrB", 100, 700, "+")]
    bam_seq, bam_qual = refs["chrA"][100:700], [10 + i % 30 for i in range(600)]    # as the BAM holds them (genome strand)
    recs = [_rec("R", "chrA", 100, 600, seq=bam_seq, qual=bam_qual, reverse=True),
            _rec("R", "chrB", 100, 600, secondary=True),
            # K: a hard-clipped primary, its SEQ 590 bases but its CIGAR 600 read bases
            _rec("K", "chrA", 100, 600, cigar="10H590M", seq=bam_seq[:590], qual=[30] * 590),
            _rec("K", "chrB", 100, 600, cigar="10H590M", secondary=True),
            # J: a tied secondary whose CIGAR holds 605 read bases against the primary's 600
            _rec("J", "chrA", 100, 600, seq=bam_seq, qual=[30] * 600),
            _rec("J", "chrB", 100, 600, cigar="5H600M", secondary=True)]
    w = write_bam(tmp_path / "w", refs, recs, copies)
    s = D.scan(w.bam, w.copies)
    assert sorted(s.tied) == ["R"]
    r = s.tied["R"]
    assert r.seq == bam_seq.translate(str.maketrans("ACGT", "TGCA"))[::-1] and list(r.qual) == bam_qual[::-1]
    assert r.copies["F.0"].flip and not r.copies["F.1"].flip and r.prim == "F.0" and r.mapq == 60 and r.margin == 0
    err = capsys.readouterr().err
    assert "; 2 skipped (a CIGAR's read bases" in err and "0 skipped (no primary SEQ)" in err


def test_scan_restricts_the_tied_reads_and_the_teachers_to_the_families_asked_for(tmp_path):
    bam, copies = make_world(tmp_path)
    f = D.scan(bam, copies, families={"F"})
    assert sorted(f.tied) == ["X", "Z"] and sorted(f.teachers) == ["F.0", "F.1", "F.2"] and f.teachers["F.1"] == set()
    g = D.scan(bam, copies, families={"G"})
    assert g.tied == {} and g.teachers == {"G.0": set()}


def test_columns_upper_case_a_soft_masked_reference_and_complement_it_for_a_flipped_copy():
    ref = {"chrA": "aaaacgtacg", "chrB": "cgtgcgtttt"}
    a0 = D.Aln("F.0", "chrA", 4, ((0, 5),), False, 5, 0, 0.0, 5)
    a1 = D.Aln("F.1", "chrB", 1, ((0, 5),), True, 3, 1, 0.2, 3)
    cols = D.columns(D.TiedRead("F", {"F.0": a0, "F.1": a1}, "F.0", 0, 2, "CGTAC", [30] * 5), lambda ch, lo, hi: ref[ch][lo:hi])
    assert [c["ref"] for c in cols] == [{"F.0": "C", "F.1": "C"}, {"F.0": "G", "F.1": "G"}, {"F.0": "T", "F.1": "C"},
                                        {"F.0": "A", "F.1": "A"}, {"F.0": "C", "F.1": "C"}]


def test_read_frame_positions_refuse_a_cigar_that_does_not_hold_the_whole_read():
    with pytest.raises(ValueError):
        D.read_frame_positions(100, [(0, 5)], 6, True)               # a hard clip left out would shift every flipped column


def make_depth_world(tmp_path):
    """Reads longcallR 2.0.0's hifi-masseq pileup treats differently, every one carrying an ALT at chrA:105, 400 and 695 and
    chrB:1000, so the caller writes a record, and its DP, there. chrA and chrB have no two equal neighbouring bases but for
    TTTTT planted at chrB:1001-1005. chrA, reads of 600M at 100 unless said: g1-g8 plain; de050 (de 0.5); de049 (de 0.49);
    nm (no de, NM 400); q2 (base quality 2: a Q0 at an ALT makes longcallR 2.0.0 panic, phase.rs:307); short (499M); len500
    (500M at 150); hard (100H450M: 450 bases of SEQ); del (a deletion at 400); nbase (N at 400); sec (a secondary with SEQ);
    supp (supplementary); dup (duplicate); mapq0 (MAPQ 0).
    chrB: f1-f7 (600M at 700: 1000 is read base 300; f1-f4 also carry an ALT at 1003, inside the run, so that the site gets
    a record), n1-n3 (600M at 990: read base 10), e1-e3 (600M at 410: read base 590), s1-s2 (50S600M at 990: read base 60,
    aligned base 10), h1 (10H50S600M at 990) and t1 (600M50S at 410: aligned base 590 of 600, 50 bases before the read's
    end). Returns write_bam's namespace."""
    refs = _no_runs(4)
    b = list(refs["chrB"])
    b[1001:1006] = "TTTTT"
    for k in (1000, 1006):                                          # the run is exactly five long
        if b[k] == "T":
            b[k] = next(x for x in "ACG" if x not in (b[k - 1], b[k + 1]))
    refs["chrB"] = "".join(b)
    A, B = _with_alts(refs["chrA"], (105, 400, 695)), _with_alts(refs["chrB"], (1000,))
    B2 = _with_alts(B, (1003,))

    def ra(name, start=100, cigar="600M", seq=None, qual=None, **kw):
        seq = A[start:start + 600] if seq is None else seq
        return _rec(name, "chrA", start, 600, cigar=cigar, seq=seq, qual=[30] * len(seq) if qual is None else qual, **kw)

    def rb(name, start, cigar="600M", seq=None):
        seq = B[start:start + 600] if seq is None else seq
        return _rec(name, "chrB", start, 600, cigar=cigar, seq=seq, qual=[30] * len(seq))

    recs = [ra(f"g{i}") for i in range(1, 9)]
    recs += [ra("de050", de=0.5), ra("de049", de=0.49), ra("nm", de=None, NM=400), ra("q2", qual=[2] * 600),
             ra("short", cigar="499M", seq=A[100:599]), ra("len500", start=150, cigar="500M", seq=A[150:650]),
             ra("hard", cigar="100H450M", seq=A[100:550]), ra("del", cigar="300M1D299M", seq=A[100:400] + A[401:700]),
             ra("nbase", seq=A[100:400] + "N" + A[401:700]), ra("sec", secondary=True), ra("supp", flags=0x800),
             ra("dup", flags=0x400), ra("mapq0", mapq=0)]
    recs += [rb(f"f{i}", 700, seq=B2[700:1300] if i <= 4 else None) for i in range(1, 8)]
    recs += [rb(f"n{i}", 990) for i in range(1, 4)]
    recs += [rb(f"e{i}", 410) for i in range(1, 4)]
    recs += [rb(f"s{i}", 990, cigar="50S600M", seq=refs["chrA"][1500:1550] + B[990:1590]) for i in range(1, 3)]
    recs += [rb("h1", 990, cigar="10H50S600M", seq=refs["chrA"][1500:1550] + B[990:1590]),
             rb("t1", 410, cigar="600M50S", seq=B[410:1010] + refs["chrA"][1500:1550])]
    return write_bam(tmp_path / "depth", refs, recs, [("F.0", "chrA", 100, 700, "+"), ("F.1", "chrB", 700, 1300, "+")])


def test_teacher_depth_counts_what_longcallr_hifi_masseq_counts(tmp_path):
    import pysam
    w = make_depth_world(tmp_path)
    fetch = pysam.FastaFile(w.genome).fetch
    # chrA (DP written by longcallR 2.0.0 on this world): g1-g8, de049, nm, q2, del (not at 400), nbase (not at 400), dup and
    # mapq0 count, len500 at 400 only; de050, short, hard, sec and supp never; read bases 5 (105) and 595 (695) count
    assert D.teacher_depth(w.bam, "chrA", [105, 400, 695], fetch) == {105: 15, 400: 14, 695: 15}
    # chrB:1000 is followed by TTTTT: it counts only where it lies >= 40 bases from both aligned ends (f1-f7), not 10 from
    # the start (n1-n3; s1-s2 and h1 after 50 soft-clipped bases) or from the end (e1-e3; t1 before 50 soft-clipped bases).
    # chrB:1003, inside the run, counts in all 17 reads: the run's base is the reference's there
    assert D.teacher_depth(w.bam, "chrB", [1000, 1003], fetch) == {1000: 7, 1003: 17}


@pytest.mark.skipif(not os.path.exists(D.LONGCALLR), reason="longcallR 2.0.0 is not installed")
def test_teacher_depth_is_the_dp_the_real_caller_writes_on_the_fixture(tmp_path):
    import pysam
    alts = [110, 125, 400, 690]                                     # read bases 10, 25, 300 and 590 of every 600-bp teacher
    bam, copies = make_world(tmp_path, alts=alts)
    teach, genome = str(tmp_path / "teach.bam"), os.path.join(os.path.dirname(bam), "genome.fa")
    D.write_teacher_bam(bam, {"F.0": ("chrA", 100, 700)}, D.scan(bam, copies).teachers, teach)
    calls = D.parse_vcf(D.run_longcallr(teach, genome, str(tmp_path / "vcf" / "F.0"), threads=1, region="chrA:1-2000"))
    assert all(("chrA", p) in calls for p in alts)
    assert D.teacher_depth(teach, "chrA", alts, pysam.FastaFile(genome).fetch) == {p: calls["chrA", p].dp for p in alts}
