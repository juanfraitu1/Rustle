import csv
import importlib.util
import json
import os
import random
import subprocess
import sys
import types

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    s = importlib.util.spec_from_file_location(f"o3_graph_{name}_t4", os.path.join(HERE, f"{name}.py"))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


O = _load("o2vc")
TD = _load("test_o2vc_data")                    # Task 3's write_bam, the one BAM writer of the o2vc tests


def family_refs(psv=()):
    """chrA (2,000 random bases, seed 11) with A at 300 and 350 and G at 420; chrB identical to chrA, so the reference has
    no PSV between F.0 and F.1, but for another base at each position of `psv`; chrM (1,000 random bases)."""
    rng = random.Random(11)
    a = [rng.choice("ACGT") for _ in range(2000)]
    a[300], a[350], a[420] = "A", "A", "G"
    b = [next(x for x in "ACGT" if x != a[i]) if i in psv else a[i] for i in range(len(a))]
    return {"chrA": "".join(a), "chrB": "".join(b), "chrM": "".join(rng.choice("ACGT") for _ in range(1000))}


def with_bases(seq, bases):
    """seq with bases ({position: base}) put in."""
    s = list(seq)
    for p, b in bases.items():
        s[p] = b
    return "".join(s)


def rec(name, chrom, start, AS, seq=None, qual=30, **kw):
    """A record for write_bam (Task 3's _rec): base quality `qual` throughout when seq is given (a list: per base; None:
    no QUAL)."""
    if seq is not None and isinstance(qual, int):
        qual = [qual] * len(seq)
    return TD._rec(name, chrom, start, AS, seq=seq, qual=qual if seq is not None else None, **kw)


def make_family_world(tmp_path, extra=(), psv=()):
    """The brief's world (family_refs). Copies, all +: F.0 = chrA:100-700, F.1 = chrB:100-700, F.2 = chrB:1000-1600 (no
    reads), M.0 = chrM:100-700, M.1 = chrA:1000-1600. Reads of 600 bp, 600M, base quality 30: t0-t7, teachers of F.0 with
    the reference; u0-u7, teachers of F.1 with chrB's reference but G at 300; X, F.0 594 (primary) / F.1 594, with chrA's
    reference but G at 300; Y1 and Y2, F.0 588 / F.1 588, with chrA's reference but C at 350 and T at 420; M1 and M2, M.0
    600 (primary) / M.1 600. `extra`: more records; `psv`: chrB positions given another base than chrA's. Returns
    write_bam's namespace (bam, copies, genome)."""
    refs = family_refs(psv)
    A = refs["chrA"]
    copies = [("F.0", "chrA", 100, 700, "+"), ("F.1", "chrB", 100, 700, "+"), ("F.2", "chrB", 1000, 1600, "+"),
              ("M.0", "chrM", 100, 700, "+"), ("M.1", "chrA", 1000, 1600, "+")]
    recs = [rec(f"t{i}", "chrA", 100, 600, seq=A[100:700]) for i in range(8)]
    recs += [rec(f"u{i}", "chrB", 100, 600, seq=with_bases(refs["chrB"], {300: "G"})[100:700]) for i in range(8)]
    recs += [rec("X", "chrA", 100, 594, seq=with_bases(A, {300: "G"})[100:700], NM=1),
             rec("X", "chrB", 100, 594, secondary=True, NM=1)]
    for q in ("Y1", "Y2"):
        recs += [rec(q, "chrA", 100, 588, seq=with_bases(A, {350: "C", 420: "T"})[100:700], NM=2),
                 rec(q, "chrB", 100, 588, secondary=True, NM=2)]
    for q in ("M1", "M2"):
        recs += [rec(q, "chrM", 100, 600, seq=refs["chrM"][100:700]), rec(q, "chrA", 1000, 600, secondary=True)]
    return TD.write_bam(tmp_path / "world", refs, recs + list(extra), copies)


VCF_HEAD = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n"


def fake_vcf_caller(prefixes, body):
    """A stand-in for o2vc_data.run_longcallr: keeps each out_prefix in prefixes, makes its directory and writes
    <out_prefix>.vcf (Task 2: where longcallR 2.0.0 writes it), holding body[chrom] (VCF lines) for a region on chrom."""
    def run(teacher_bam, genome, out_prefix, binary=None, threads=2, region=None):
        prefixes.append(out_prefix)
        os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
        with open(out_prefix + ".vcf", "w") as fh:
            fh.write(VCF_HEAD + body.get(region.split(":")[0], ""))
        return out_prefix + ".vcf"
    return run


def fake_caller(prefixes, ref="A"):
    """fake_vcf_caller holding one record, for a chrB region only: chrB 301 {ref}>G, 1/1, DP 8."""
    return fake_vcf_caller(prefixes, {"chrB": f"chrB\t301\t.\t{ref}\tG\t.\tPASS\t.\tGT:GQ:DP:AF\t1/1:30:8:1.00\n"})


def test_run_places_a_tied_read_by_the_sample_s_call_where_the_reference_alone_cannot_and_flags_a_shared_foreign_pair(tmp_path, monkeypatch):
    w, prefixes = make_family_world(tmp_path), []
    monkeypatch.setattr(O.D, "run_longcallr", fake_caller(prefixes))
    assert O.cmd_run(w.bam, w.copies, w.genome, str(tmp_path / "out"), exclude_chrM=True) == 0
    calls = {r["read"]: r for r in csv.DictReader(open(tmp_path / "out" / "calls.tsv"), delimiter="\t")}
    assert calls["X"]["vc"] == "F.1" and calls["X"]["vc_cols"] == "1" and calls["X"]["n3"] == ""
    assert calls["Y1"]["flag"] == calls["Y2"]["flag"] == "F#0" and calls["Y1"]["vc"] == ""
    assert not any(r["family"] == "M" for r in calls.values())                       # a chrM copy: family left out
    assert not any(p.endswith("F.2") or "/M." in p for p in prefixes)                # F.2 has no teachers: never called
    kinds = {r["kind"] for r in csv.DictReader(open(tmp_path / "out" / "F.columns.tsv"), delimiter="\t") if r["read"] == "X"}
    assert kinds == {"separating_call"}


def test_run_stops_when_a_record_s_ref_disagrees_with_the_genome(tmp_path, monkeypatch):
    w = make_family_world(tmp_path)
    monkeypatch.setattr(O.D, "run_longcallr", fake_caller([], ref="C"))
    with pytest.raises(ValueError, match="chrB:301"):
        O.cmd_run(w.bam, w.copies, w.genome, str(tmp_path / "out"))


# ---- beyond the brief's tests: the controller's rulings (task-4-rulings.md) and behaviour the brief states but does not pin ----

def _rows(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def test_teacher_depth_needs_the_genome_and_fails_at_once_without_it(tmp_path):
    w = make_family_world(tmp_path)                                 # ruling 3: fetch is required, never looked for late
    for args in ((w.bam, "chrA", [300]), (w.bam, "chrA", [300], None)):
        with pytest.raises(TypeError):
            O.D.teacher_depth(*args)


def test_run_raises_teacher_base_qualities_of_0_to_1_and_counts_them_and_the_tied_reads_without_qual(tmp_path, monkeypatch):
    import pysam
    A = family_refs()["chrA"]
    q0 = [30] * 600
    q0[10:15] = [0] * 5                                             # longcallR 2.0.0 panics on a Q0 at a variant site
    w = make_family_world(tmp_path, extra=[rec("t8", "chrA", 100, 600, seq=A[100:700], qual=q0),
                                           rec("N", "chrA", 100, 600, seq=A[100:700], qual=None),
                                           rec("N", "chrB", 100, 600, secondary=True)])
    monkeypatch.setattr(O.D, "run_longcallr", fake_caller([]))
    assert O.cmd_run(w.bam, w.copies, w.genome, str(tmp_path / "out"), exclude_chrM=True) == 0
    with pysam.AlignmentFile(str(tmp_path / "out" / "teach" / "F.bam")) as bam:
        quals = {r.query_name: list(r.query_qualities) for r in bam}
    assert quals["t8"][9:16] == [30, 1, 1, 1, 1, 1, 30] and min(min(v) for v in quals.values()) == 1 and len(quals) == 17
    run = json.load(open(tmp_path / "out" / "run.json"))
    assert run["q0_raised"] == dict(reads=1, bases=5) and run["tied_without_qual"] == 1     # N: QUAL 255 votes everywhere


def test_run_reuses_a_saved_vcf_and_never_calls_longcallr_again(tmp_path, monkeypatch):
    w, first, again, out = make_family_world(tmp_path), [], [], str(tmp_path / "out")
    monkeypatch.setattr(O.D, "run_longcallr", fake_caller(first))
    O.cmd_run(w.bam, w.copies, w.genome, out, exclude_chrM=True)
    with open(os.path.join(out, "vcf", "F.1.vcf"), "w") as fh:
        fh.write(VCF_HEAD)                                          # the saved VCF is the record, even saying less
    monkeypatch.setattr(O.D, "run_longcallr", fake_caller(again))
    assert O.cmd_run(w.bam, w.copies, w.genome, out, exclude_chrM=True) == 0
    calls = {r["read"]: r for r in _rows(os.path.join(out, "calls.tsv"))}
    assert len(first) == 2 and again == [] and calls["X"]["vc"] == "" and calls["Y1"]["flag"] == "F#0"
    assert json.load(open(os.path.join(out, "run.json")))["reused_vcfs"] == ["F.0", "F.1"]


def failing_on(chrom, prefixes):
    """fake_caller, but raising run_longcallr's RuntimeError for a region on chrom."""
    ok = fake_caller(prefixes)

    def run(teacher_bam, genome, out_prefix, binary=None, threads=2, region=None):
        if region.startswith(chrom + ":"):
            prefixes.append(out_prefix)
            raise RuntimeError(f"longcallR exited 101 ({region}):\nthread '<unnamed>' panicked at src/phase.rs:307:52")
        return ok(teacher_bam, genome, out_prefix, binary, threads, region)
    return run


def test_a_copy_whose_call_fails_is_unknown_everywhere_and_the_run_goes_on(tmp_path, monkeypatch):
    w, prefixes, out = make_family_world(tmp_path), [], tmp_path / "out"
    monkeypatch.setattr(O.D, "run_longcallr", failing_on("chrB", prefixes))
    assert O.cmd_run(w.bam, w.copies, w.genome, str(out), exclude_chrM=True) == 0
    calls = {r["read"]: r for r in _rows(out / "calls.tsv")}
    assert calls["X"]["vc"] == "" and calls["X"]["vc_set"] == "F.0,F.1" and calls["Y1"]["flag"] == calls["Y2"]["flag"] == ""
    assert _rows(out / "F.columns.tsv") == []                       # F.1 unknown: nothing separates, nothing is foreign
    run = json.load(open(out / "run.json"))
    assert list(run["failed_copies"]) == ["F.1"] and "exited 101" in run["failed_copies"]["F.1"]
    assert len(prefixes) == 2 and not os.path.exists(out / "vcf" / "F.1.vcf")         # nothing saved: a rerun calls again


def test_the_caller_region_is_the_span_padded_1000_with_start_1_based_and_end_exclusive_inside_the_contig():
    assert O.caller_region("chr1", 5000, 6000, 10 ** 6) == ("chr1:4001-7000", (4001, 6999))
    assert O.caller_region("chr1", 500, 6000, 6500) == ("chr1:1-6501", (1, 6500))     # START 0 and END > length + 1 panic


def test_longcallr_gets_each_copy_s_padded_span_and_that_contig_s_fasta_alone(tmp_path, monkeypatch):
    import pysam
    B = family_refs()["chrB"]
    w = make_family_world(tmp_path, extra=[rec(f"v{i}", "chrB", 1000, 600, seq=B[1000:1600]) for i in range(6)])
    seen, fake = [], fake_caller([])

    def caller(teacher_bam, genome, out_prefix, binary=None, threads=2, region=None):
        seen.append((os.path.basename(out_prefix), region, genome, binary, threads))
        return fake(teacher_bam, genome, out_prefix, binary, threads, region)
    monkeypatch.setattr(O.D, "run_longcallr", caller)
    assert O.cmd_run(w.bam, w.copies, w.genome, str(tmp_path / "out"), exclude_chrM=True, longcallr="LC", threads=3) == 0
    ref = str(tmp_path / "out" / "ref")
    assert seen == [("F.0", "chrA:1-1700", f"{ref}/chrA.fa", "LC", 3), ("F.1", "chrB:1-1700", f"{ref}/chrB.fa", "LC", 3),
                    ("F.2", "chrB:1-2001", f"{ref}/chrB.fa", "LC", 3)]      # F.2 (teachers v0-v5): 1600 + 1000 > 2,000 + 1
    fa = pysam.FastaFile(f"{ref}/chrB.fa")
    assert list(fa.references) == ["chrB"] and fa.fetch("chrB") == B


def test_a_tied_copy_without_teachers_is_unknown_so_it_is_never_beaten_and_blocks_foreign_alleles(tmp_path, monkeypatch):
    z = with_bases(family_refs()["chrB"], {350: "C", 420: "T"})[100:700]
    extra = []
    for q in ("Z1", "Z2"):                                          # Y's two alleles, tied on F.1 and teacherless F.2
        extra += [rec(q, "chrB", 100, 588, seq=z, NM=2), rec(q, "chrB", 1000, 588, secondary=True)]
    w = make_family_world(tmp_path, extra=extra)
    monkeypatch.setattr(O.D, "run_longcallr", fake_caller([]))
    O.cmd_run(w.bam, w.copies, w.genome, str(tmp_path / "out"), exclude_chrM=True)
    calls = {r["read"]: r for r in _rows(tmp_path / "out" / "calls.tsv")}
    assert [(calls[q]["vc"], calls[q]["vc_set"], calls[q]["flag"]) for q in ("Z1", "Z2")] == [("", "F.1,F.2", "")] * 2
    assert [calls[q]["n3"] for q in ("Z1", "Z2")] == ["", ""]          # N3 knows no more than VC: F.2 stays unknown
    assert calls["Y1"]["flag"] == "F#0"


def test_a_reference_psv_the_sample_confirms_places_the_read_under_both_vc_and_n3(tmp_path, monkeypatch):
    A = family_refs()["chrA"]
    w = make_family_world(tmp_path, psv=(500,), extra=[rec("P", "chrA", 100, 600, seq=A[100:700]),
                                                        rec("P", "chrB", 100, 594, secondary=True, NM=1)])
    monkeypatch.setattr(O.D, "run_longcallr", fake_caller([]))
    O.cmd_run(w.bam, w.copies, w.genome, str(tmp_path / "out"), exclude_chrM=True)
    p = {r["read"]: r for r in _rows(tmp_path / "out" / "calls.tsv")}["P"]
    assert (p["vc"], p["vc_cols"], p["n3"], p["margin"]) == ("F.0", "2", "F.0", "6")
    assert [(r["p"], r["kind"]) for r in _rows(tmp_path / "out" / "F.columns.tsv") if r["read"] == "P"] == [
        ("200", "separating_call"), ("400", "separating_ref")]     # chrB:300's call; chrA:500 vs chrB:500, confirmed
    f0 = {r["copy"]: r for r in _rows(tmp_path / "out" / "F.reads.tsv") if r["read"] == "P"}["F.0"]
    assert (f0["ref_for"], f0["ref_against"], f0["call_for"], f0["call_against"]) == ("1", "0", "2", "0")


def test_the_case_files_hold_tags_votes_and_columns_and_the_teacher_side_is_reported_apart(tmp_path, monkeypatch):
    odd = with_bases(family_refs()["chrA"], {350: "C", 420: "T"})[100:700]
    w = make_family_world(tmp_path, extra=[rec(q, "chrA", 100, 600, seq=odd) for q in ("ta", "tb")])     # teachers of F.0
    monkeypatch.setattr(O.D, "run_longcallr", fake_caller([]))
    O.cmd_run(w.bam, w.copies, w.genome, str(tmp_path / "out"), exclude_chrM=True)
    out = tmp_path / "out"
    x = {r["copy"]: r for r in _rows(out / "F.reads.tsv") if r["read"] == "X"}
    assert [x["F.0"][k] for k in ("primary", "AS", "NM", "subs", "indels", "call_for", "call_against")] == ["1", "594", "1", "1",
                                                                                                         "0", "0", "1"]
    assert [x["F.1"][k] for k in ("primary", "ref_for", "ref_against", "call_for", "call_against")] == ["0", "0", "0", "1", "0"]
    assert [r["copies"] for r in _rows(out / "F.columns.tsv") if r["read"] == "X"] == ["F.0|chrA|300|A|A|10;F.1|chrB|300|A|G|8"]
    assert _rows(out / "F.flags.tsv") == [dict(flag="F#0", n_reads="2", reads="Y1,Y2",
                                               alleles="F.0|chrA|350|C;F.0|chrA|420|T;F.1|chrB|350|C;F.1|chrB|420|T")]
    t = {r["copy"]: r for r in _rows(out / "F.teachers.tsv")}
    assert [t["F.0"]["teachers"], t["F.0"]["teacher_flag_reads"], t["F.1"]["records"], t["F.1"]["pass"]] == ["10", "2", "1", "1"]
    calls = {r["read"]: r for r in _rows(out / "calls.tsv")}
    assert calls["Y1"]["flag"] == "F#0" and calls["X"]["vc"] == "F.1"          # ta and tb's own flag changes no call
    run = json.load(open(out / "run.json"))
    assert (run["preset"], run["min_mapq"], run["families"], run["excluded_chrM"], run["depth_agreement"]) == (
        "hifi-masseq", 0, ["F"], ["M"], [1, 1])


def test_the_cli_runs_cmd_run_with_its_options(monkeypatch):
    seen = {}
    monkeypatch.setattr(O, "cmd_run", lambda *a, **k: seen.update(a=a, k=k) or 0)
    assert O.main(["run", "b.bam", "c.tsv", "g.fa", "out", "--families", "F,G", "--exclude-chrM", "--threads", "3"]) == 0
    assert seen == dict(a=("b.bam", "c.tsv", "g.fa", "out"), k=dict(families={"F", "G"}, exclude_chrM=True, threads=3,
                                                                     weak_het=O.WEAK_HET_DEFAULT))
    assert O.main(["run", "b.bam", "c.tsv", "g.fa", "out", "--weak-het", "unknown"]) == 0
    assert seen["k"]["weak_het"] == "unknown"
    with pytest.raises(SystemExit):
        O.main(["run", "b.bam", "c.tsv", "g.fa", "out", "--weak-het", "drop"])                # not a policy
    with pytest.raises(SystemExit):
        O.main(["run", "b.bam"])


def test_a_teacher_side_flag_needs_its_alleles_known(tmp_path, monkeypatch):
    B = family_refs()["chrB"]
    alt = {p: next(x for x in "ACGT" if x != B[p]) for p in (1100, 1200)}
    w = make_family_world(tmp_path, extra=[rec(q, "chrB", 1000, 600, seq=with_bases(B, alt)[1000:1600]) for q in ("v0", "v1")])
    monkeypatch.setattr(O.D, "run_longcallr", fake_caller([]))
    O.cmd_run(w.bam, w.copies, w.genome, str(tmp_path / "out"), exclude_chrM=True)
    t = {r["copy"]: r for r in _rows(tmp_path / "out" / "F.teachers.tsv")}
    assert (t["F.2"]["teachers"], t["F.2"]["teacher_flag_reads"]) == ("2", "0")       # depth 2 < 6: F.2 unknown there


def make_inverted_world(tmp_path):
    """Family F with an inverted copy (Review Focus 1): chrA as family_refs' and chrC its reverse complement (chrA:p is
    chrC:1999-p); F.0 = chrA:100-700 (+) and F.3 = chrC:1300-1900 (-), its image. Reads of 600 bp, 600M: t0-t7, teachers
    of F.0 with the reference; w0-w7, reverse teachers of F.3 with chrC's reference but C at 1699 (G at chrA:300 in the
    transcript's frame); X, F.0 594 (forward primary) / F.3 594 (reverse), with G at chrA:300; Y1 and Y2, F.0 588 / F.3 588
    (reverse), with C at chrA:350 and T at chrA:420."""
    A = family_refs()["chrA"]
    C = A.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    recs = [rec(f"t{i}", "chrA", 100, 600, seq=A[100:700]) for i in range(8)]
    recs += [rec(f"w{i}", "chrC", 1300, 600, seq=with_bases(C, {1699: "C"})[1300:1900], reverse=True) for i in range(8)]
    recs += [rec("X", "chrA", 100, 594, seq=with_bases(A, {300: "G"})[100:700], NM=1),
             rec("X", "chrC", 1300, 594, secondary=True, reverse=True, NM=1)]
    for q in ("Y1", "Y2"):
        recs += [rec(q, "chrA", 100, 588, seq=with_bases(A, {350: "C", 420: "T"})[100:700], NM=2),
                 rec(q, "chrC", 1300, 588, secondary=True, reverse=True, NM=2)]
    return TD.write_bam(tmp_path / "inv", {"chrA": A, "chrC": C}, recs,
                        [("F.0", "chrA", 100, 700, "+"), ("F.3", "chrC", 1300, 1900, "-")])


def test_a_copy_aligned_opposite_to_the_read_is_compared_in_the_read_s_frame_and_flagged_in_the_genome_s(tmp_path, monkeypatch):
    w = make_inverted_world(tmp_path)
    call = "chrC\t1700\t.\tT\tC\t.\tPASS\t.\tGT:GQ:DP:AF\t1/1:30:8:1.00\n"          # F.3's T>C: G at chrA:300 in X's frame
    monkeypatch.setattr(O.D, "run_longcallr", fake_vcf_caller([], {"chrC": call}))
    O.cmd_run(w.bam, w.copies, w.genome, str(tmp_path / "out"))
    calls = {r["read"]: r for r in _rows(tmp_path / "out" / "calls.tsv")}
    assert (calls["X"]["vc"], calls["X"]["n3"], calls["Y1"]["flag"], calls["Y2"]["flag"]) == ("F.3", "", "F#0", "F#0")
    assert [(r["kind"], r["copies"]) for r in _rows(tmp_path / "out" / "F.columns.tsv") if r["read"] == "X"] == [
        ("separating_call", "F.0|chrA|300|A|A|8;F.3|chrC|1699|A|G|8")]          # chrC's T and C, complemented
    assert [r["copies"] for r in _rows(tmp_path / "out" / "F.columns.tsv") if r["read"] == "Y1" and r["p"] == "250"] == [
        "F.0|chrA|350|A|A|8;F.3|chrC|1649|A|A|8"]                            # no record: chrC's T, in the read's frame
    assert _rows(tmp_path / "out" / "F.flags.tsv")[0]["alleles"] == "F.0|chrA|350|C;F.0|chrA|420|T;F.3|chrC|1579|A;F.3|chrC|1649|G"


def test_a_copy_is_unknown_outside_its_called_region_over_the_caller_s_max_depth_and_at_a_reference_n():
    call = O.R.Call("A", ("G",), (0, 1), "PASS", 7)
    cps = O.Copies({"F.0": ("chrA", 100, 700)}, {"F.0": {("chrA", 300): call}}, {"F.0": (1, 1699)}, None, None)
    cps.depth["F.0"].update({0: 9, 1: 9, 400: 50001, 500: 5, 600: 9, 1698: 9, 1699: 9})

    def g(p, ref="A"):
        return cps.genotype("F.0", "chrA", p, ref)
    assert g(300) == (frozenset("AG"), 7)                                       # a record: its GT and its DP
    assert [g(p) for p in (1, 600, 1698)] == [(frozenset("A"), 9)] * 3
    assert [g(p) for p in (0, 1699, 400, 500)] == [(None, None), (None, None), (None, 50001), (None, 5)]
    assert g(600, "N") == (None, None) and cps.genotype("F.1", "chrB", 600, "A") == (None, None)       # F.1: no call


def test_a_copy_nothing_beats_is_not_placed_unless_it_beats_every_other(tmp_path, monkeypatch):
    A = family_refs()["chrA"]                                       # ruling 1: a placement is assign()'s first element only
    recs = [rec(f"{ch}{i}", ch, 100, 600, seq=A[100:700]) for ch in ("chrA", "chrB", "chrD") for i in range(8)]
    recs += [rec("K", "chrA", 100, 600, seq=A[100:700]), rec("K", "chrB", 100, 600, secondary=True),
             rec("K", "chrD", 100, 600, secondary=True)]
    w = TD.write_bam(tmp_path / "chain", {"chrA": A, "chrB": A, "chrD": A}, recs,
                     [("F.0", "chrA", 100, 700, "+"), ("F.1", "chrB", 100, 700, "+"), ("F.2", "chrD", 100, 700, "+")])
    het, hom = "LowQual\t.\tGT:GQ:DP:AF\t0/1:5:8:0.50\n", "PASS\t.\tGT:GQ:DP:AF\t1/1:30:8:1.00\n"
    body = {"chrA": f"chrA\t351\t.\tA\tT\t.\t{het}",                                 # F.0 at 350: A/T
            "chrB": f"chrB\t301\t.\tA\tG\t.\t{hom}",                                 # F.1 at 300: G
            "chrD": f"chrD\t301\t.\tA\tG\t.\t{het}chrD\t351\t.\tA\tT\t.\t{hom}"}     # F.2: A/G at 300, T at 350
    monkeypatch.setattr(O.D, "run_longcallr", fake_vcf_caller([], body))
    O.cmd_run(w.bam, w.copies, w.genome, str(tmp_path / "out"))
    k = _rows(tmp_path / "out" / "calls.tsv")[0]                    # K (A at 300 and 350): F.0 beats F.1, F.1 beats F.2,
    assert (k["read"], k["vc"], k["vc_set"], k["vc_cols"], k["n3_set"]) == ("K", "", "F.0", "0", "F.0,F.1,F.2")  # F.0 ~ F.2


def test_a_reference_psv_the_sample_does_not_carry_is_erased_so_n3_and_vc_decide_apart(tmp_path, monkeypatch):
    refs = family_refs(psv=(500,))                                  # chrA:500 and chrB:500 differ in the reference
    A, B = refs["chrA"], refs["chrB"]
    a, b = A[500], B[500]
    c = next(x for x in "ACGT" if x not in (a, b) and not (x == "G" and a == "A"))     # foreign there, not editing
    recs = [rec(f"t{i}", "chrA", 100, 600, seq=A[100:700]) for i in range(8)]
    recs += [rec(f"u{i}", "chrB", 100, 600, seq=with_bases(B, {500: a})[100:700]) for i in range(8)]   # copy B carries A's base
    recs += [rec("P", "chrA", 100, 600, seq=A[100:700]), rec("P", "chrB", 100, 594, secondary=True, NM=1)]
    for q in ("Q1", "Q2"):                                          # a base neither copy holds at 500, a foreign C at 350
        recs += [rec(q, "chrA", 100, 588, seq=with_bases(A, {350: "C", 500: c})[100:700], NM=2),
                 rec(q, "chrB", 100, 588, secondary=True, NM=2)]
    w = TD.write_bam(tmp_path / "erased", {"chrA": A, "chrB": B}, recs,
                     [("F.0", "chrA", 100, 700, "+"), ("F.1", "chrB", 100, 700, "+")])
    call = f"chrB\t501\t.\t{b}\t{a}\t.\tPASS\t.\tGT:GQ:DP:AF\t1/1:30:8:1.00\n"
    monkeypatch.setattr(O.D, "run_longcallr", fake_vcf_caller([], {"chrB": call}))
    O.cmd_run(w.bam, w.copies, w.genome, str(tmp_path / "out"))
    calls = {r["read"]: r for r in _rows(tmp_path / "out" / "calls.tsv")}
    cols = {(r["read"], r["p"]): r for r in _rows(tmp_path / "out" / "F.columns.tsv")}
    assert (calls["P"]["n3"], calls["P"]["vc"], calls["P"]["vc_set"]) == ("F.0", "", "F.0,F.1")     # the reference's PSV
    assert [(p, r["kind"], r["copies"]) for (q, p), r in sorted(cols.items()) if q == "P"] == [     # places P, the sample's
        ("400", "erased_ref", f"F.0|chrA|500|{a}|{a}|8;F.1|chrB|500|{b}|{a}|8")]                    # calls do not
    assert (cols["Q1", "400"]["kind"], cols["Q1", "250"]["kind"]) == ("erased_ref", "foreign")       # erased_ref first,
    assert calls["Q1"]["flag"] == calls["Q2"]["flag"] == "F#0"      # and its foreign allele still counts toward the flag


# ---- Task 5: o2vc score. The brief's tests (all but score-worlds: rulings 1), the rulings' extras, then behaviour the brief
# ---- states but does not pin. A toy run directory is what `run` leaves: calls.tsv, <family>.reads.tsv, <family>.columns.tsv.

CALLS_HEAD = "read family candidates minimap2 mapq margin outside vc vc_set vc_cols n3 n3_set flag"
READS_HEAD = "read copy primary AS NM de ms subs indels ref_for ref_against call_for call_against"
COLUMNS_HEAD = "read p base qual kind copies"
TRUTH_HEAD = "read family tied_copies class truth donor_nm donor_margin"
CASE_ADDED = ["vc", "n3", "flag", "truth_class", "truth", "vc_outcome", "score_options"]


def _put(path, head, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(head.split()) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def tied_read(read, mm, truth, vc="", n3="", flag="", mapq=0, cols=0, kinds=(), fam="F", cands=None, tags=None, vc_set=None,
              n3_set=None):
    """One tied read of a toy run directory. mm: minimap2's primary copy ("" off every copy); vc, n3: the placed copy, "" when
    abstaining; flag: the flag ids, comma-joined; cols: the columns voting for VC; kinds: the kinds of its columns.tsv rows (a
    kind, or (kind, base quality) -- 30 otherwise); tags: {copy: (AS, NM, subs)}; vc_set, n3_set: the copies no copy beats (the
    placed copy when None); truth: (class, truth copy or None), None for a read without a truth row."""
    cands = cands or [f"{fam}.0", f"{fam}.1"]
    return dict(read=read, mm=mm, truth=truth, vc=vc, n3=n3, flag=flag, mapq=mapq, cols=cols, kinds=list(kinds), fam=fam,
                cands=cands, tags=tags or {c: (1000, 2, 1) for c in cands}, vc_set=vc if vc_set is None else vc_set,
                n3_set=n3 if n3_set is None else n3_set)


def _qual_kind(k):
    """(base quality, kind) of a `kinds` entry of tied_read: a kind alone has quality 30, (kind, quality) its own."""
    return (30, k) if isinstance(k, str) else (k[1], k[0])


def write_run(out, reads):
    """calls.tsv and, per family, reads.tsv and columns.tsv for `reads` (tied_read dicts), as o2vc run writes them."""
    os.makedirs(out, exist_ok=True)
    _put(f"{out}/calls.tsv", CALLS_HEAD, [
        [t["read"], t["fam"], ",".join(t["cands"]), t["mm"], t["mapq"], 0, int(t["mm"] not in t["cands"]), t["vc"], t["vc_set"],
         t["cols"], t["n3"], t["n3_set"], t["flag"]] for t in reads])
    for f in sorted({t["fam"] for t in reads}):
        mine = [t for t in reads if t["fam"] == f]
        _put(f"{out}/{f}.reads.tsv", READS_HEAD, [
            [t["read"], c, int(c == t["mm"]), t["tags"][c][0], "" if t["tags"][c][1] is None else t["tags"][c][1], 0.01, 900,
             "" if t["tags"][c][2] is None else t["tags"][c][2], 0, 0, 0, 0, 0] for t in mine for c in t["cands"]])
        _put(f"{out}/{f}.columns.tsv", COLUMNS_HEAD, [
            [t["read"], 100 + i, "A", *_qual_kind(k), f"{t['cands'][0]}|chrA|{300 + i}|A|A|10;{t['cands'][1]}|chrB|{300 + i}|A|G|8"]
            for t in mine for i, k in enumerate(t["kinds"])])


def write_truth(path, reads):
    rows = []
    for t in reads:
        if t["truth"] is not None:
            rows.append([t["read"], t["fam"], ",".join(t["cands"]), t["truth"][0], t["truth"][1] or "", 0, ""])
    _put(path, TRUTH_HEAD, rows)


def toy_dir(tmp_path, reads, par=()):
    """A run directory, its donor truth table and a PAR list for `reads`, as a namespace (out, truth, par, reads)."""
    out = str(tmp_path / "run")
    write_run(out, reads)
    truth, parp = str(tmp_path / "truth.tsv"), str(tmp_path / "par.json")
    write_truth(truth, reads)
    json.dump(list(par), open(parp, "w"))
    return types.SimpleNamespace(out=out, truth=truth, par=parp, reads=reads)


def make_score_toy(tmp_path):
    """The brief's toy: family F (copies F.0, F.1) with reads r1-r7, and family G (one decidable read, VC right) in par.json."""
    reads = [
        tied_read("r1", "F.0", ("unique", "F.1"), vc="F.1", cols=1, kinds=["separating_call"],
                  tags={"F.0": (1000, 3, 2), "F.1": (1000, 3, 1)}),             # AS and NM tie; the substitutions pick F.1
        tied_read("r2", "F.0", ("unique", "F.0"), vc="F.1", n3="F.1", mapq=3, cols=2, kinds=["separating_ref"] * 2),
        tied_read("r3", "F.1", ("unique", "F.1")),
        tied_read("r4", "F.0", ("own_tie", None), vc="F.0", cols=1, kinds=["separating_ref"]),
        tied_read("r5", "F.0", ("no_counterpart", None), flag="F#0", kinds=["foreign"]),
        tied_read("r6", "F.0", ("no_counterpart", None), flag="F#0", kinds=["foreign"]),
        tied_read("r7", "F.0", ("unique", "F.0"), flag="F#1", kinds=["foreign"]),
        tied_read("g1", "G.0", ("unique", "G.1"), vc="G.1", cols=1, kinds=["separating_call"], fam="G"),
    ]
    return toy_dir(tmp_path, reads, par=["G"])


def _score(d, name="score.A.F.json", **kw):
    """Score the toy `d` (donor truth; the PAR list; panel A on F unless kw says otherwise), return its JSON."""
    args = dict(panel="A", focus="F", exclude_p=d.par)
    args.update(kw)
    assert O.cmd_score(d.out, [d.truth], "donor", **args) == 0
    return json.load(open(os.path.join(d.out, name)))


def test_score_counts_each_arm_on_decidable_reads_and_the_calibration_bars(tmp_path):
    # The brief's test, its toy expectations recomputed by hand for the 2026-10-22 amendment (task-5 fix round 1): the arms are
    # over the SCORED reads -- r1, r2, r3, r7 (unique) and r5, r6 (no_counterpart: every placement wrong) -- not the four
    # decidable ones. VC places r1 (right) and r2 (wrong) and abstains on r3 and, flagged, on r5, r6, r7: abstained 2 -> 4.
    d = make_score_toy(tmp_path)
    assert O.cmd_score(d.out, [d.truth], "donor", panel="A", focus="F", exclude_p=d.par) == 0
    s = json.load(open(os.path.join(d.out, "score.A.F.json")))
    assert s["arms"]["VC"] == dict(placed=2, right=1, wrong=1, abstained=4, fixes=1, breaks=1, both_wrong=0)
    assert s["arms"]["N2"]["placed"] == 1 and s["arms"]["N2"]["right"] == 1
    assert s["K1"] == dict(reads=1, placed=1) and s["K2"] == dict(reads=5, flagged=1)
    assert [f["classes"] for f in s["flags"] if f["id"] == "F#0"] == [{"no_counterpart": 2}]
    assert s["bars"]["A1"]["pass"] is False and s["bars"]["A3"]["pass"] is True and s["bars"]["K1"]["pass"] is False
    assert "G" in s["excluded"] and "G" not in s["families"]


def test_evidence_names_a_sample_call_when_any_deciding_column_is_one():
    assert O.evidence({"separating_call", "separating_ref"}) == "sample call"
    assert O.evidence({"separating_ref"}) == "reference PSV"
    assert O.evidence({"separating_ref", "erased_ref"}) == "reference PSV corrected by the sample"      # rulings 8


def test_score_mode_none_reports_placed_abstained_and_flagged_per_arm_with_no_truth_and_no_right_or_wrong(tmp_path):
    d = make_score_toy(tmp_path)
    assert O.cmd_score(d.out, [], "none") == 0                       # the DAZ run: zero truth paths, so one pool of G and F
    s = json.load(open(os.path.join(d.out, "score.json")))
    assert (s["mode"], s["panel"], s["tied"], s["bars"], s["verdict"]) == ("none", None, 8, {}, None)
    assert s["arms"] == dict(N1=dict(placed=8, abstained=0, flagged=3), N2=dict(placed=1, abstained=7, flagged=3),
                             N3=dict(placed=1, abstained=7, flagged=3), VC=dict(placed=4, abstained=4, flagged=3))
    assert s["vc"]["placements"] == 4 and s["vc"]["single_column"] == dict(placements=3)         # r1, r4, g1 rest on one column
    assert not any(k in s for k in ("K1", "K2", "D2", "classes", "decidable", "scored", "precision", "fixes"))
    case = {r["read"]: r for r in _rows(os.path.join(d.out, "F.case.tsv"))}
    assert [case[q]["vc_outcome"] for q in ("r1", "r2", "r3", "r4", "r5", "r7")] == [
        "placed", "placed", "abstain", "placed", "flagged", "flagged"]
    assert {(r["truth_class"], r["truth"]) for r in case.values()} == {("", "")}               # no truth was read
    assert {r["score_options"] for r in case.values()} == {"mode=none truth=- n4=- panel=- focus=- exclude=-"}
    text = open(os.path.join(d.out, "F.summary.md")).read()
    assert "Score options: mode=none truth=- n4=- panel=- focus=- exclude=-" in text
    head = text.split("## Worked reads")[0].lower()
    assert "right" not in head and "wrong" not in head                        # nothing to be right or wrong about
    assert "### First VC placement: r1" in text and "### First abstention: r3" in text and "### First flagged read: r5" in text


def test_score_json_reports_every_arm_its_precision_with_counts_and_the_class_flag_and_fix_tables(tmp_path):
    d = make_score_toy(tmp_path)
    s = _score(d)
    assert (s["mode"], s["panel"], s["focus"], s["tied"], s["scored"]) == ("donor", "A", "F", 7, 6)   # r4 is own_tie: K1 only
    assert list(s["arms"]) == ["N1", "N2", "N3", "VC"]                                      # no --n4: no N4
    # the arms are over the six scored reads r1, r2, r3, r7 (unique) and r5, r6 (no_counterpart: every placement wrong)
    assert s["arms"]["N1"] == dict(placed=6, right=3, wrong=3, abstained=0, fixes=0, breaks=0, both_wrong=3)
    assert s["arms"]["N2"] == dict(placed=1, right=1, wrong=0, abstained=5, fixes=0, breaks=0, both_wrong=0)
    assert s["arms"]["N3"] == dict(placed=1, right=0, wrong=1, abstained=5, fixes=0, breaks=1, both_wrong=0)
    assert s["precision"]["VC"] == dict(wrong=1, placed=2, rate=0.5) and s["precision"]["N1"] == dict(wrong=3, placed=6, rate=0.5)
    # the decidable subset (truth copy among the tied copies: r1, r2, r3, r7), reported for comparison with earlier registers
    assert s["decidable"]["reads"] == 4
    assert s["decidable"]["arms"]["N1"] == dict(placed=4, right=3, wrong=1, abstained=0, fixes=0, breaks=0, both_wrong=1)
    assert s["decidable"]["arms"]["VC"] == dict(placed=2, right=1, wrong=1, abstained=2, fixes=1, breaks=1, both_wrong=0)
    assert s["decidable"]["precision"]["VC"] == dict(wrong=1, placed=2, rate=0.5)
    assert s["outside"]["tied"] == 0 and s["outside"]["scored"] == 0 and s["truth_outside_candidates"] == 0
    assert s["D2"] == dict(reads=3, moved=1, flagged=1)                                    # N1 right on r2, r3, r7
    assert s["classes"] == {"no_counterpart": dict(reads=2, placed=0, flagged=2), "own_tie": dict(reads=1, placed=1, flagged=0),
                            "unique": dict(reads=4, placed=2, flagged=1)}
    assert s["C1"] == dict(flags=2, counted=1)
    assert [(f["id"], f["family"], f["reads"], f["classes"]) for f in s["flags"]] == [
        ("F#0", "F", 2, {"no_counterpart": 2}), ("F#1", "F", 1, {"unique": 1})]
    assert s["fixes"] == dict(reads=1, evidence={"sample call": 1, "reference PSV corrected by the sample": 0, "reference PSV": 0},
                              tags_alone=dict(AS=0, NM=0, subs=1, none=0))
    assert s["vc"]["placements"] == 3 and s["vc"]["single_column"] == dict(placements=2, scored=1, right=1, wrong=0)
    assert s["vc"]["erased_ref"] == dict(placements=0, scored=0, right=0, wrong=0)
    assert s["vc"]["not_scored"] == dict(placements=1, classes={"own_tie": 1})                 # r4: placed, K1 only
    assert list(s["families"]) == ["F"] and s["families"]["F"]["arms"] == s["arms"]
    assert s["excluded"]["G"]["arms"]["VC"] == dict(placed=1, right=1, wrong=0, abstained=0, fixes=1, breaks=0, both_wrong=0)
    assert s["bars"]["A1"] == {"value": 0.5, "num": 1, "den": 2, "limit": 0.02, "op": "<=", "pass": False}
    assert set(s["bars"]) == {"A1", "A2", "A3", "K1", "K2"} and s["verdict"] == "FAIL"
    assert (s["bars"]["A2"]["value"], s["bars"]["A2"]["limit"], s["bars"]["A2"]["pass"]) == (1, 10, None)   # 3 N1-wrong < 10
    assert s["bars"]["A2"]["reason"] == "the focus family has 3 scored reads minimap2 places wrong, fewer than 10"
    assert (s["bars"]["A3"]["value"], s["bars"]["A3"]["limit"], s["bars"]["A3"]["pass"]) == (1, 0, True)       # N3 has no fix
    assert (s["bars"]["K1"]["num"], s["bars"]["K1"]["den"]) == (1, 1) and (s["bars"]["K2"]["num"], s["bars"]["K2"]["den"]) == (1, 5)


def test_score_prints_each_arm_with_its_counts_and_each_bar_with_its_verdict(tmp_path, capsys):
    d = make_score_toy(tmp_path)
    _score(d)
    out = capsys.readouterr().out
    assert "7 tied reads, 6 scored" in out
    assert "N1: placed 6, right 3, wrong 3, abstained 0; wrong/placed 3/6 = 0.5000" in out
    assert "VC: placed 2, right 1, wrong 1, abstained 4; wrong/placed 1/2 = 0.5000; fixes 1, breaks 1, both wrong 0" in out
    assert "decidable subset (4 reads): VC placed 2, right 1, wrong 1; N1 right 3, wrong 1" in out
    assert "K1: VC places 1 of 1 own_tie reads" in out and "K2: 1 of 5 reads with a present copy are flagged" in out
    assert "D2: 3 scored reads N1 places right; VC moves 1, flags 1" in out
    assert "C1: 1 of 2 flags hold >= 2 no_counterpart" in out
    assert "bar A1: wrong/placed 1/2 = 0.5000 <= 0.02: FAIL" in out and "bar A3: " in out and "verdict FAIL" in out
    assert "primary outside the tied copies: 0 tied reads, 0 scored" in out


def test_the_case_file_is_reads_tsv_plus_the_decision_the_truth_and_the_outcome(tmp_path):
    d = make_score_toy(tmp_path)
    _score(d)
    rows = _rows(os.path.join(d.out, "F.case.tsv"))
    assert list(rows[0]) == READS_HEAD.split() + CASE_ADDED
    assert [r["read"] for r in rows] == [q for q in ("r1", "r2", "r3", "r4", "r5", "r6", "r7") for _ in range(2)]     # reads.tsv's
    by = {r["read"]: r for r in rows}
    assert {q: by[q]["vc_outcome"] for q in by} == dict(r1="fix", r2="break", r3="abstain", r4="placed", r5="flagged",
                                                        r6="flagged", r7="flagged")
    assert [(by[q]["vc"], by[q]["n3"], by[q]["flag"]) for q in ("r1", "r2", "r5", "r7")] == [
        ("F.1", "", ""), ("F.1", "F.1", ""), ("", "", "F#0"), ("", "", "F#1")]
    assert [(by[q]["truth_class"], by[q]["truth"]) for q in ("r1", "r4", "r5")] == [
        ("unique", "F.1"), ("own_tie", ""), ("no_counterpart", "")]
    assert [r["AS"] for r in rows if r["read"] == "r1"] == ["1000", "1000"]                  # reads.tsv's own fields are kept
    assert os.path.exists(os.path.join(d.out, "G.case.tsv"))                                 # an excluded family is scored apart


def test_both_right_and_both_wrong_are_told_apart_and_a_primary_off_every_copy_is_a_wrong_placement(tmp_path):
    c3 = ["F.0", "F.1", "F.2"]
    d = toy_dir(tmp_path, [
        tied_read("b1", "F.0", ("unique", "F.0"), vc="F.0", cols=1, cands=c3),               # VC confirms a right primary
        tied_read("b2", "F.0", ("unique", "F.1"), vc="F.0", cols=1, cands=c3),               # VC confirms a wrong primary
        tied_read("b3", "F.0", ("unique", "F.1"), vc="F.2", cols=1, cands=c3),               # VC moves it to a third copy
        tied_read("b4", "", ("unique", "F.1"), vc="F.1", cols=1, cands=c3),                  # primary on no copy at all
        tied_read("b5", "H.0", ("unique", "F.0"), cands=c3),                                 # primary on another family's copy
    ])
    s = _score(d, "score.json", panel=None, focus=None, exclude_p=None)
    assert s["arms"]["N1"] == dict(placed=5, right=1, wrong=4, abstained=0, fixes=0, breaks=0, both_wrong=4)   # always placed
    assert s["arms"]["VC"] == dict(placed=4, right=2, wrong=2, abstained=1, fixes=1, breaks=0, both_wrong=2)
    assert s["outside"]["tied"] == 2 and s["outside"]["scored"] == 2                         # b4, b5: "primary outside" row
    assert s["outside"]["arms"]["N1"]["right"] == 0 and s["outside"]["arms"]["N1"]["wrong"] == 2
    assert s["outside"]["arms"]["VC"] == dict(placed=1, right=1, wrong=0, abstained=1, fixes=1, breaks=0, both_wrong=0)
    assert s["D2"] == dict(reads=1, moved=0, flagged=0)
    case = {r["read"]: r["vc_outcome"] for r in _rows(os.path.join(d.out, "F.case.tsv"))}
    assert case == dict(b1="both right", b2="both wrong", b3="both wrong", b4="fix", b5="abstain")
    text = open(os.path.join(d.out, "F.summary.md")).read()
    assert "Reads whose minimap2 primary lies outside the tied copies: 2 of 5 tied reads, 2 scored." in text
    assert text.count("| VC variant calls |") == 3                       # the arm table, the outside row's, the decidable subset's


def test_every_tied_read_with_a_definite_truth_is_scored_and_the_others_are_counted_apart(tmp_path):
    d = toy_dir(tmp_path, [
        tied_read("a1", "F.0", ("unique", "F.1"), vc="F.1", cols=1),                         # scored, a tied copy is the truth
        tied_read("a2", "F.0", ("unique", "F.2"), vc="F.1", cols=1),                         # scored, truth copy not tied
        tied_read("a3", "F.0", None, vc="F.0", cols=1),                                      # no truth row: counted only
        tied_read("a4", "F.0", ("poor_match", None), vc="F.0", cols=1),                      # counted only
        tied_read("a5", "F.0", ("elsewhere", None)),                                         # scored: every placement is wrong
        tied_read("a6", "F.0", ("unique_elsewhere", "F.1"), vc="F.1", cols=1),               # scored, a fix
        tied_read("a7", "F.0", ("no_counterpart", None), vc="F.1", cols=1),                  # scored: VC placed it, wrong
        tied_read("a8", "F.0", ("own_tie", None), vc="F.0", cols=1),                         # K1 only
        tied_read("a9", "F.0", ("unmapped", None)),                                          # counted only
    ])
    s = _score(d, "score.json", panel=None, focus=None, exclude_p=None)
    assert (s["tied"], s["scored"], s["decidable"]["reads"], s["truth_outside_candidates"]) == (9, 5, 2, 1)
    # scored: a1, a2, a5, a6, a7. N1 is right on none of them (a2's truth copy is F.2; a5, a7 have no truth copy)
    assert s["arms"]["N1"] == dict(placed=5, right=0, wrong=5, abstained=0, fixes=0, breaks=0, both_wrong=5)
    assert s["arms"]["VC"] == dict(placed=4, right=2, wrong=2, abstained=1, fixes=2, breaks=0, both_wrong=2)
    assert s["decidable"]["arms"]["VC"] == dict(placed=2, right=2, wrong=0, abstained=0, fixes=2, breaks=0, both_wrong=0)
    assert s["classes"] == {"elsewhere": dict(reads=1, placed=0, flagged=0), "no_counterpart": dict(reads=1, placed=1, flagged=0),
                            "no_truth_row": dict(reads=1, placed=1, flagged=0), "own_tie": dict(reads=1, placed=1, flagged=0),
                            "poor_match": dict(reads=1, placed=1, flagged=0), "unique": dict(reads=2, placed=2, flagged=0),
                            "unique_elsewhere": dict(reads=1, placed=1, flagged=0), "unmapped": dict(reads=1, placed=0, flagged=0)}
    assert s["K1"] == dict(reads=1, placed=1) and s["K2"] == dict(reads=4, flagged=0)       # a8; a1, a2, a6, a8: copy present
    assert s["vc"]["not_scored"] == dict(placements=3, classes={"no_truth_row": 1, "own_tie": 1, "poor_match": 1})
    case = {r["read"]: r["vc_outcome"] for r in _rows(os.path.join(d.out, "F.case.tsv"))}
    assert case == dict(a1="fix", a2="both wrong", a3="placed", a4="placed", a5="abstain", a6="fix", a7="both wrong", a8="placed",
                        a9="abstain")
    assert {r["read"]: r["truth_class"] for r in _rows(os.path.join(d.out, "F.case.tsv"))}["a3"] == "no_truth_row"


def test_a_certainly_wrong_vc_placement_reaches_the_bars_the_reviewers_probe(tmp_path):
    # Panel E on 79 VC placements, 30 of them wrong by construction: 20 reads whose minimap2 primary F.2 lies outside the tied
    # copies and IS the donor's copy (VC moves them: a break), and 10 no_counterpart reads VC places. Counting only the
    # decidable reads (truth among the tied copies) gave E1 0/49 and E2 2 > 0: PASS.
    d = toy_dir(tmp_path, (
        [tied_read(f"a{i:03d}", "F.0", ("unique", "F.1"), vc="F.1", cols=1) for i in range(2)]                    # fixes
        + [tied_read(f"b{i:03d}", "F.0", ("unique", "F.0"), vc="F.0", cols=1) for i in range(47)]                 # both right
        + [tied_read(f"c{i:03d}", "F.2", ("unique_elsewhere", "F.2"), vc="F.0", cols=1) for i in range(20)]       # breaks
        + [tied_read(f"d{i:03d}", "F.0", ("no_counterpart", None), vc="F.1", cols=1) for i in range(10)]))        # wrong
    s = _score(d, "score.E.json", panel="E", focus=None, exclude_p=None)
    assert s["arms"]["VC"] == dict(placed=79, right=49, wrong=30, abstained=0, fixes=2, breaks=20, both_wrong=10)
    assert s["arms"]["N1"] == dict(placed=79, right=67, wrong=12, abstained=0, fixes=0, breaks=0, both_wrong=12)
    assert (s["tied"], s["scored"], s["decidable"]["reads"], s["truth_outside_candidates"]) == (79, 79, 49, 20)
    assert s["decidable"]["arms"]["VC"] == dict(placed=49, right=49, wrong=0, abstained=0, fixes=2, breaks=0, both_wrong=0)
    assert s["outside"]["tied"] == 20 and s["outside"]["scored"] == 20                        # the primaries on F.2
    assert s["outside"]["arms"]["N1"]["right"] == 20 and s["outside"]["arms"]["VC"] == dict(
        placed=20, right=0, wrong=20, abstained=0, fixes=0, breaks=20, both_wrong=0)
    assert s["D2"] == dict(reads=67, moved=20, flagged=0)                                    # N1 right: b and the primaries on F.2
    assert s["bars"]["E1"]["num"] == 30 and s["bars"]["E1"]["den"] == 79 and s["bars"]["E1"]["pass"] is False
    assert (s["bars"]["E2"]["value"], s["bars"]["E2"]["limit"], s["bars"]["E2"]["pass"]) == (2, 20, False)
    assert s["bars"]["K2"]["pass"] is True and s["verdict"] == "FAIL"                         # not PASS


def test_an_excluded_family_is_scored_apart_and_never_enters_the_pooled_numbers(tmp_path):
    d = toy_dir(tmp_path, [
        tied_read("f1", "F.0", ("unique", "F.0"), vc="F.0", cols=1),
        tied_read("p1", "P.0", ("own_tie", None), vc="P.0", cols=1, fam="P"),                # K1 would count it
        tied_read("p2", "P.0", ("no_counterpart", None), flag="P#0", fam="P"),               # C1 would count the flag
        tied_read("p3", "P.0", ("no_counterpart", None), flag="P#0", fam="P"),
        tied_read("p4", "P.0", ("unique", "P.0"), flag="P#1", fam="P"),                      # K2 would count it
    ], par=["P"])
    s = _score(d, "score.C.json", panel="C", focus=None)
    assert list(s["families"]) == ["F"] and list(s["excluded"]) == ["P"] and s["tied"] == 1
    assert s["K1"] == dict(reads=0, placed=0) and s["K2"] == dict(reads=1, flagged=0)
    assert s["C1"] == dict(flags=0, counted=0) and s["flags"] == [] and s["arms"]["VC"]["placed"] == 1
    assert s["bars"]["C1"]["pass"] is None and s["verdict"] == "UNDERPOWERED"                # the pool holds no no_counterpart read
    p = s["excluded"]["P"]
    assert p["tied"] == 3 + 1 and p["K1"] == dict(reads=1, placed=1) and p["C1"] == dict(flags=2, counted=1)
    assert os.path.exists(os.path.join(d.out, "P.case.tsv")) and os.path.exists(os.path.join(d.out, "P.summary.md"))


def test_the_exclude_file_may_be_a_json_list_or_a_json_object_keyed_by_family(tmp_path):
    d = make_score_toy(tmp_path)
    keyed = str(tmp_path / "keyed.json")
    json.dump({"G": {"why": "PAR"}}, open(keyed, "w"))
    assert _score(d, exclude_p=keyed)["excluded"].keys() == _score(d)["excluded"].keys() == {"G"}


def test_n4_is_the_correctors_pick_and_minimap2s_for_a_read_it_does_not_hold(tmp_path):
    d = make_score_toy(tmp_path)
    n4 = str(tmp_path / "o2correct.tsv")
    _put(n4, "read family corrector minimap2 reason candidates margin scores",
         [["r1", "F", "F.1", "F.0", "sample", "F.0,F.1", 5, ""], ["r2", "F", "F.1", "F.0", "sample", "F.0,F.1", 4, ""]])
    s = _score(d, "score.json", panel=None, focus=None, n4_p=n4)
    assert list(s["arms"]) == ["N1", "N2", "N3", "N4", "VC"]
    # r1: F.1 right (N1 wrong) = fix; r2: F.1 wrong (N1 right) = break; r3, r7 absent: N1's, right; r5, r6 (no_counterpart)
    # absent: N1's, wrong
    assert s["arms"]["N4"] == dict(placed=6, right=3, wrong=3, abstained=0, fixes=1, breaks=1, both_wrong=2)
    assert "N4 corrector" in open(os.path.join(d.out, "F.summary.md")).read()


def n4_file(tmp_path, reads):
    """An o2correct-style --n4 file holding `reads` (each picking F.1)."""
    path = str(tmp_path / "o2correct.tsv")
    _put(path, "read family corrector minimap2 reason candidates margin scores",
         [[q, "F", "F.1", "F.0", "sample", "F.0,F.1", 5, ""] for q in reads])
    return path


def test_how_many_tied_reads_the_n4_file_holds_is_reported_so_a_silent_fallback_to_n1_is_visible(tmp_path, capsys):
    d = make_score_toy(tmp_path)
    s = _score(d, "score.json", panel=None, focus=None, n4_p=n4_file(tmp_path, ["r1", "r2", "zz"]))    # zz: not a tied read
    assert s["n4_held"] == 2 and s["families"]["F"]["n4_held"] == 2 and s["excluded"]["G"]["n4_held"] == 0
    assert "N4: the corrector holds 2 of 7 tied reads; the others are N1's" in capsys.readouterr().out
    assert "N4: the corrector holds 2 of 7 tied reads of this family; the others are N1's." in open(
        os.path.join(d.out, "F.summary.md")).read()
    assert "n4_held" not in _score(d, "score.json", panel=None, focus=None)                             # no --n4: nothing to report


def test_an_n4_file_that_holds_none_of_the_tied_reads_is_warned_about(tmp_path, capsys):
    d = make_score_toy(tmp_path)
    s = _score(d, "score.json", panel=None, focus=None, n4_p=n4_file(tmp_path, ["zz"]))
    assert s["n4_held"] == 0 and s["arms"]["N4"] == s["arms"]["N1"]                                     # N4 is N1 everywhere
    assert "holds none of the tied reads" in capsys.readouterr().err


def test_only_columns_that_can_vote_set_the_evidence_and_the_erased_counts(tmp_path):
    d = toy_dir(tmp_path, [                                       # a base votes from quality R.MIN_BASEQ = 10
        tied_read("m1", "H.0", ("unique", "H.1"), vc="H.1", cols=1, fam="H", kinds=[("separating_call", 9), "separating_ref"]),
        tied_read("m2", "H.0", ("unique", "H.1"), vc="H.1", cols=1, fam="H", kinds=[("separating_call", 10), "separating_ref"]),
        tied_read("m3", "H.0", ("unique", "H.0"), vc="H.0", cols=1, fam="H", kinds=[("erased_ref", 9), "separating_ref"]),
        tied_read("m4", "H.0", ("unique", "H.0"), vc="H.0", cols=1, fam="H", kinds=[("erased_ref", 10), "separating_ref"])])
    s = _score(d, "score.json", panel=None, focus=None, exclude_p=None)
    assert s["fixes"]["evidence"] == {"sample call": 1, "reference PSV corrected by the sample": 0, "reference PSV": 1}   # m2; m1
    assert s["vc"]["erased_ref"] == dict(placements=1, scored=1, right=1, wrong=0)                    # m4 only: m3's column is Q9
    assert s["vc"]["erased_dependent"] == dict(placements=1, scored=1, right=1, wrong=0)              # N3 abstains on m4
    text = open(os.path.join(d.out, "H.summary.md")).read()
    assert "| m1 | H.1 | H.0 | reference PSV | " in text and "| m2 | H.1 | H.0 | sample call | " in text
    assert "m1\t100\tA\t9\tseparating_call" in text                                    # the worked read lists every row


def test_the_score_options_are_recorded_in_the_summary_and_in_every_case_row(tmp_path):
    d = make_score_toy(tmp_path)
    n4 = n4_file(tmp_path, ["r1"])
    assert O.cmd_score(d.out, [d.truth, d.truth], "donor", n4_p=n4, panel="A", focus="F", exclude_p=d.par) == 0
    want = f"mode=donor truth={d.truth},{d.truth} n4={n4} panel=A focus=F exclude={d.par}"
    for fam in ("F", "G"):                                                                            # G is excluded: scored apart
        assert {r["score_options"] for r in _rows(os.path.join(d.out, f"{fam}.case.tsv"))} == {want}
        assert f"Score options: {want}" in open(os.path.join(d.out, f"{fam}.summary.md")).read()
    s = json.load(open(os.path.join(d.out, "score.A.F.json")))
    assert (s["truth_files"], s["n4_file"], s["exclude_file"]) == ([d.truth, d.truth], n4, d.par)


def test_each_fix_names_its_evidence_and_the_tags_that_alone_would_have_picked_it_and_erased_columns_are_counted(tmp_path):
    d = toy_dir(tmp_path, [
        tied_read("h1", "H.0", ("unique", "H.1"), vc="H.1", cols=1, kinds=["separating_call", "erased_ref"], fam="H",
                  tags={"H.0": (1000, 4, 1), "H.1": (1000, 3, 1)}),                          # a fix; only NM picks H.1
        tied_read("h2", "H.0", ("unique", "H.1"), vc="H.1", cols=2, kinds=["separating_ref"] * 2 + ["erased_ref"], fam="H",
                  tags={"H.0": (1005, 2, 1), "H.1": (1000, 2, 1)}),                          # a fix; AS picks H.0: wrong
        tied_read("h3", "H.0", ("unique", "H.1"), vc="H.1", n3="H.1", cols=1, kinds=["separating_ref"], fam="H",
                  tags={"H.0": (1000, 3, 3), "H.1": (1000, 3, 2)}),                          # a fix; only subs picks H.1
        tied_read("h4", "H.0", ("unique", "H.0"), vc="H.0", n3="H.1", cols=1, kinds=["separating_ref", "erased_ref"], fam="H"),
        tied_read("h5", "H.0", ("unique", "H.0"), vc="H.1", cols=1, kinds=["erased_ref", "separating_ref"], fam="H"),
        tied_read("h6", "H.0", ("unique", "H.0"), vc="H.0", n3="H.0", cols=3, fam="H",
                  kinds=["separating_ref"] * 3 + ["erased_ref"]),
    ])
    s = _score(d, "score.json", panel=None, focus=None)
    assert s["fixes"] == dict(reads=3, evidence={"sample call": 1, "reference PSV corrected by the sample": 1, "reference PSV": 1},
                              tags_alone=dict(AS=0, NM=1, subs=1, none=1))
    assert s["arms"]["VC"] == dict(placed=6, right=5, wrong=1, abstained=0, fixes=3, breaks=1, both_wrong=0)
    assert s["arms"]["N3"] == dict(placed=3, right=2, wrong=1, abstained=3, fixes=1, breaks=1, both_wrong=0)
    assert s["vc"] == dict(placements=6,
                           single_column=dict(placements=4, scored=4, right=3, wrong=1),             # h1 h3 h4 h5
                           erased_ref=dict(placements=5, scored=5, right=4, wrong=1),                # all but h3
                           erased_dependent=dict(placements=4, scored=4, right=3, wrong=1),          # h6: N3 agrees with VC
                           not_scored=dict(placements=0, classes={}))
    text = open(os.path.join(d.out, "H.summary.md")).read()
    assert "| h1 | H.1 | H.0 | sample call | no | yes | no |" in text
    assert "| h2 | H.1 | H.0 | reference PSV corrected by the sample | no | no | no |" in text
    assert "| h3 | H.1 | H.0 | reference PSV | no | no | yes |" in text
    assert "resting on a single column (`vc_cols == 1`): 4 of 6" in text
    assert "whose read has an `erased_ref` column: 5 of 6" in text
    assert "depending on erased columns" in text and ": 4 of 6 (scored 4: 3 right, 1 wrong)" in text


def test_the_summary_holds_the_arm_table_the_fixes_the_calibration_and_the_legend(tmp_path):
    d = make_score_toy(tmp_path)
    _score(d)
    text = open(os.path.join(d.out, "F.summary.md")).read()
    assert "Scored against the donor's truth. 7 tied reads (no_counterpart 2, own_tie 1, unique 4)." in text
    assert "Scored: 6 (no_counterpart 2, unique 4). Not scored: own_tie 1 (K1 only)." in text
    assert "| N1 minimap2 | 6 | 3 | 3 | 0 | 3/6 = 0.5000 | - | - | - |" in text
    assert "| N2 minimap2 at MAPQ >= 1 | 1 | 1 | 0 | 5 | 0/1 = 0.0000 | 0 | 0 | 0 |" in text
    assert "| N3 reference only | 1 | 0 | 1 | 5 | 1/1 = 1.0000 | 0 | 1 | 0 |" in text
    assert "| VC variant calls | 2 | 1 | 1 | 4 | 1/2 = 0.5000 | 1 | 1 | 0 |" in text
    assert "Decidable subset (the truth copy is one of the tied copies; the definition of earlier registers): 4 reads" in text
    assert "| N1 minimap2 | 4 | 3 | 1 | 0 | 1/4 = 0.2500 | - | - | - |" in text                  # the subset's own table
    assert "| VC variant calls | 2 | 1 | 1 | 2 | 1/2 = 0.5000 | 1 | 1 | 0 |" in text
    assert "compare arms by fixes and breaks" in text                                      # precision is not comparable across arms
    assert "which the arm table cannot see: 1 (own_tie 1)." in text                        # r4: placed, K1 only
    assert "| r1 | F.1 | F.0 | sample call | no | no | yes |" in text
    assert "- K1: VC places 1 of 1 own_tie reads." in text and "- K2: 1 of 5 reads with a present copy are flagged." in text
    assert "- D2: of 3 scored reads N1 places right, VC moves 1 and flags 1." in text
    assert "- flag F#0: reads 2 (no_counterpart 2)" in text and "- flag F#1: reads 1 (unique 1)" in text
    assert [p.split("\n")[0] for p in text.split("### ")[1:]] == ["First fix: r1", "First abstention: r3", "First flagged read: r5"]


def test_the_summary_shows_the_worked_reads_the_first_fix_abstention_and_flag_by_read_name(tmp_path):
    d = toy_dir(tmp_path, [                                                                  # file order is not name order
        tied_read("z9", "W.0", ("unique", "W.1"), vc="W.1", cols=1, kinds=["separating_ref"], fam="W"),
        tied_read("a3", "W.0", ("unique", "W.1"), vc="W.1", cols=1, kinds=["separating_call"], fam="W"),
        tied_read("m1", "W.0", ("unique", "W.0"), kinds=["separating_ref"] * 2, fam="W"),
        tied_read("c2", "W.1", ("unique", "W.1"), kinds=["separating_ref"] * 3, fam="W"),
        tied_read("q5", "W.0", ("unique", "W.0"), flag="W#0", kinds=["foreign"], fam="W"),
        tied_read("b4", "W.0", ("unique", "W.0"), flag="W#1", kinds=["foreign"] * 2, fam="W"),
    ])
    _score(d, "score.json", panel=None, focus=None)
    text = open(os.path.join(d.out, "W.summary.md")).read()
    parts = text.split("### ")[1:]
    assert [p.split("\n")[0] for p in parts] == ["First fix: a3", "First abstention: c2", "First flagged read: b4"]
    fix, abst, flag = parts
    assert fix.count("separating_call") == 1 and "separating_ref" not in fix and "VC W.1" in fix and "sample call" in fix
    assert abst.count("separating_ref") == 3 and "separating_call" not in abst and "VC abstains" in abst
    assert flag.count("foreign") == 2 and "W#1" in flag and "VC abstains" in flag
    assert "a3\tW.0\t1\t1000" in fix and "a3\tW.1\t0\t1000" in fix                           # its reads.tsv rows, verbatim


def test_n2_is_minimap2_only_from_mapq_one(tmp_path):
    d = toy_dir(tmp_path, [tied_read("m0", "F.0", ("unique", "F.0"), mapq=0), tied_read("m1", "F.0", ("unique", "F.0"), mapq=1),
                           tied_read("m2", "F.0", ("unique", "F.1"), mapq=60)])
    s = _score(d, "score.json", panel=None, focus=None, exclude_p=None)
    assert s["arms"]["N2"] == dict(placed=2, right=1, wrong=1, abstained=1, fixes=0, breaks=0, both_wrong=1)     # m1, m2 not m0


def test_d2_counts_only_the_reads_minimap2_places_right_that_vc_moves_or_flags(tmp_path):
    d = toy_dir(tmp_path, [
        tied_read("x1", "F.0", ("unique", "F.0"), flag="F#0"),                          # N1 right, flagged
        tied_read("x2", "F.0", ("unique", "F.1"), flag="F#0"),                          # N1 wrong, flagged: not in D2
        tied_read("x3", "F.0", ("unique", "F.0"), vc="F.1", cols=1),                    # N1 right, moved
        tied_read("x4", "F.0", ("unique", "F.0"), vc="F.0", cols=1),                    # N1 right, confirmed
        tied_read("x5", "F.0", ("unique", "F.1"), vc="F.1", cols=1),                    # N1 wrong, fixed: not in D2
        tied_read("x6", "F.0", ("unique", "F.0")),                                      # N1 right, VC abstains unflagged
        tied_read("x7", "F.1", ("unique", "F.1"), vc="F.0", cols=1)])                   # N1 right, moved
    s = _score(d, "score.json", panel=None, focus=None, exclude_p=None)
    assert s["D2"] == dict(reads=5, moved=2, flagged=1)
    assert "- D2: of 5 scored reads N1 places right, VC moves 2 and flags 1." in open(os.path.join(d.out, "F.summary.md")).read()


def test_a_placement_is_a_non_empty_vc_or_n3_never_inferred_from_the_set(tmp_path):
    c3 = ["F.0", "F.1", "F.2"]
    d = toy_dir(tmp_path, [
        tied_read("c1", "F.0", ("unique", "F.0"), cands=c3, vc_set="F.0", n3_set="F.0"),       # beat chain: kept {F.0}
        tied_read("c2", "F.0", ("unique", "F.0"), cands=c3, vc_set="", n3_set=""),             # beat cycle: empty set
        tied_read("c3", "F.0", ("unique", "F.0"), cands=c3, vc="F.0", n3="F.0", cols=1)])
    s = _score(d, "score.json", panel=None, focus=None, exclude_p=None)
    assert s["arms"]["VC"]["placed"] == 1 and s["arms"]["VC"]["abstained"] == 2
    assert s["arms"]["N3"]["placed"] == 1 and s["arms"]["N3"]["abstained"] == 2


def verdict_of(passes):
    """A panel's verdict from its bars' `pass` values: FAIL when any bar failed (evaluably), else UNDERPOWERED when any bar had
    nothing to measure (None), else PASS."""
    v = list(passes.values())
    return "FAIL" if any(p is False for p in v) else "UNDERPOWERED" if any(p is None for p in v) else "PASS"


def make_bulk(tmp_path, spec):
    """A run directory of family F from [(prefix, n, minimap2, truth, kwargs of tied_read)]: n reads each."""
    reads = [tied_read(f"{p}{i:03d}", mm, truth, **kw) for p, n, mm, truth, kw in spec for i in range(n)]
    return toy_dir(tmp_path, reads)


def panel_a_toy(tmp_path, fix=10, ok=38, brk=0, tie=50, tie_placed=1, flagged=2, n3_fix=False, miss=0):
    """`fix` reads VC places right where minimap2 is wrong (N3 too when n3_fix), `ok` both right, `brk` VC wrong where minimap2 is
    right, `tie` own_tie reads of which `tie_placed` VC places, `flagged` flagged reads whose copy is present, `miss` reads
    minimap2 places wrong and VC leaves."""
    return make_bulk(tmp_path, [
        ("a", fix, "F.0", ("unique", "F.1"), dict(vc="F.1", n3="F.1" if n3_fix else "", cols=1)),
        ("g", miss, "F.0", ("unique", "F.1"), {}),
        ("b", ok, "F.0", ("unique", "F.0"), dict(vc="F.0", cols=1)),
        ("c", brk, "F.0", ("unique", "F.0"), dict(vc="F.1", cols=1)),
        ("d", tie - tie_placed, "F.0", ("own_tie", None), {}),
        ("e", tie_placed, "F.0", ("own_tie", None), dict(vc="F.0", cols=1)),
        ("f", flagged, "F.0", ("unique", "F.0"), dict(flag="F#0"))])


@pytest.mark.parametrize("kw, expect", [
    ({}, dict(A1=True, A2=True, A3=True, K1=True, K2=True)),                      # 10 fixes, 0 of 48 wrong, 1 of 50, 2 of 100
    (dict(fix=9, ok=39), dict(A1=True, A2=None, A3=True, K1=True, K2=True)),      # 9 fixes of 9 minimap2-wrong: unreachable
    (dict(fix=9, ok=38, miss=1), dict(A1=True, A2=False, A3=True, K1=True, K2=True)),     # 9 fixes of 10 minimap2-wrong
    (dict(ok=39, brk=1), dict(A1=True, A2=True, A3=True, K1=True, K2=True)),      # 1 wrong of 50 placed = 0.02: on the bar
    (dict(brk=1), dict(A1=False, A2=True, A3=True, K1=True, K2=True)),            # 1 wrong of 49: over
    (dict(tie_placed=2), dict(A1=True, A2=True, A3=True, K1=False, K2=True)),     # 2 of 50 own_tie reads placed
    (dict(flagged=3), dict(A1=True, A2=True, A3=True, K1=True, K2=False)),        # 3 of 101 flagged
    (dict(n3_fix=True), dict(A1=True, A2=True, A3=False, K1=True, K2=True)),      # N3 fixes as many as VC: not more
    (dict(tie=0, tie_placed=0, flagged=0), dict(A1=True, A2=True, A3=True, K1=None, K2=True)),     # no own_tie read: K1 null
])
def test_panel_a_bars_hold_at_their_boundaries(tmp_path, kw, expect):
    s = _score(panel_a_toy(tmp_path, **kw), exclude_p=None)
    assert {b: v["pass"] for b, v in s["bars"].items()} == expect
    assert s["verdict"] == verdict_of(expect)                                       # K1 null (no own_tie read): UNDERPOWERED


def panel_e_toy(tmp_path, fix=2, ok=47, brk=1, flagged=0):
    return make_bulk(tmp_path, [
        ("a", fix, "F.0", ("unique", "F.1"), dict(vc="F.1", cols=1)),
        ("b", ok, "F.0", ("unique", "F.0"), dict(vc="F.0", cols=1)),
        ("c", brk, "F.0", ("unique", "F.0"), dict(vc="F.1", cols=1)),
        ("f", flagged, "F.0", ("unique", "F.0"), dict(flag="F#0"))])


@pytest.mark.parametrize("kw, expect", [
    ({}, dict(E1=True, E2=True, K2=True)),                                        # 1 wrong of 50; 2 fixes > 1 break; 0 flagged
    (dict(ok=46), dict(E1=False, E2=True, K2=True)),                              # 1 wrong of 49
    (dict(fix=1, ok=48), dict(E1=True, E2=False, K2=True)),                       # fixes = breaks
    (dict(flagged=2), dict(E1=True, E2=True, K2=False)),                          # 2 of 52 flagged
])
def test_panel_e_bars_hold_at_their_boundaries(tmp_path, kw, expect):
    s = _score(panel_e_toy(tmp_path, **kw), "score.E.json", panel="E", focus=None, exclude_p=None)
    assert {b: v["pass"] for b, v in s["bars"].items()} == expect and "K1" in s     # K1 is reported, not a bar of E
    assert s["verdict"] == verdict_of(expect)


def test_panel_c_needs_a_flag_with_two_no_counterpart_reads_and_an_excluded_family_s_flag_does_not_count(tmp_path):
    d = make_score_toy(tmp_path)
    s = _score(d, "score.C.json", panel="C", focus=None)
    assert s["bars"] == {"C1": {"value": 1, "limit": 1, "op": ">=", "pass": True}} and s["verdict"] == "PASS"
    par_f = str(tmp_path / "par_f.json")
    json.dump(["F"], open(par_f, "w"))
    s = _score(d, "score.C.json", panel="C", focus=None, exclude_p=par_f)                    # the pool, G: no no_counterpart
    assert s["bars"]["C1"]["pass"] is None and s["bars"]["C1"]["reason"] and s["verdict"] == "UNDERPOWERED"
    assert s["excluded"]["F"]["C1"] == dict(flags=2, counted=1)                              # F's flag is reported, apart


def test_c1_is_evaluable_and_fails_when_there_are_no_counterpart_reads_but_no_flag_holds_two(tmp_path):
    d = toy_dir(tmp_path, [tied_read("n1", "F.0", ("no_counterpart", None)),
                           tied_read("n2", "F.0", ("no_counterpart", None), flag="F#0")])
    s = _score(d, "score.C.json", panel="C", focus=None, exclude_p=None)
    assert s["bars"]["C1"] == {"value": 0, "limit": 1, "op": ">=", "pass": False} and s["verdict"] == "FAIL"


def test_c1_counts_the_flags_holding_at_least_two_no_counterpart_reads_in_flag_order(tmp_path):
    nc = ("no_counterpart", None)
    d = toy_dir(tmp_path, [
        tied_read("a1", "F.0", nc, flag="F#0"), tied_read("a2", "F.0", ("unique", "F.0"), flag="F#0"),      # one: not counted
        tied_read("b1", "F.0", nc, flag="F#1"), tied_read("b2", "F.0", nc, flag="F#1"),                      # two: counted
        tied_read("c1", "F.0", nc, flag="F#2"), tied_read("c2", "F.0", nc, flag="F#2"),
        tied_read("d1", "F.0", nc, flag="F#2,F#10"),                                                         # in two flags: both
        tied_read("e1", "F.0", ("unique", "F.0"), flag="F#10")])
    s = _score(d, "score.C.json", panel="C", focus=None, exclude_p=None)
    assert [(f["id"], f["reads"]) for f in s["flags"]] == [("F#0", 2), ("F#1", 2), ("F#2", 3), ("F#10", 2)]   # not F#10 before F#2
    assert [f["classes"] for f in s["flags"]] == [{"no_counterpart": 1, "unique": 1}, {"no_counterpart": 2},
                                                  {"no_counterpart": 3}, {"no_counterpart": 1, "unique": 1}]
    assert s["C1"] == dict(flags=4, counted=2) and s["bars"]["C1"]["value"] == 2


def test_panel_a_judges_the_focus_family_alone_and_pools_the_calibration_over_every_family_scored(tmp_path):
    d = toy_dir(tmp_path, (
        [tied_read(f"f{i:02d}", "F.0", ("unique", "F.0"), vc="F.0", cols=1) for i in range(5)]               # F: 5 placed, right
        + [tied_read(f"h{i:02d}", "H.0", ("unique", "H.0"), vc="H.1", cols=1, fam="H") for i in range(5)]      # H: 5 placed, wrong
        + [tied_read(f"t{i:02d}", "H.0", ("own_tie", None), vc="H.0", cols=1, fam="H") for i in range(2)]))   # H: own_tie placed
    s = _score(d, "score.A.F.json", exclude_p=None)
    assert s["arms"]["VC"]["wrong"] == 5 and s["families"]["F"]["arms"]["VC"]["wrong"] == 0
    assert (s["bars"]["A1"]["num"], s["bars"]["A1"]["den"], s["bars"]["A1"]["pass"]) == (0, 5, True)         # F alone
    assert (s["bars"]["K1"]["num"], s["bars"]["K1"]["den"]) == (2, 2)                                        # H's reads count


def test_a_tag_not_measured_on_a_copy_picks_nothing(tmp_path):
    d = toy_dir(tmp_path, [tied_read("n1", "F.0", ("unique", "F.1"), vc="F.1", cols=1, kinds=["separating_ref"],
                                     tags={"F.0": (1000, None, None), "F.1": (1000, 3, 1)})])
    s = _score(d, "score.json", panel=None, focus=None, exclude_p=None)
    assert s["fixes"]["reads"] == 1 and s["fixes"]["tags_alone"] == dict(AS=0, NM=0, subs=0, none=1)       # no NM on F.0: no pick
    assert "| n1 | F.1 | F.0 | reference PSV | no | no | no |" in open(os.path.join(d.out, "F.summary.md")).read()


def test_a_bar_with_nothing_to_measure_is_not_evaluable_and_the_panel_is_underpowered_not_failed(tmp_path):
    # This was the test that pinned E2 False at zero VC placements (a panel in which VC places nothing "failed"). Ruling (fix
    # round 1, finding 2): with no VC placement there is nothing to measure, so E1 and E2 are not evaluable.
    d = toy_dir(tmp_path, [tied_read("u1", "F.0", ("unique", "F.0")), tied_read("u2", "F.0", ("unique", "F.1"))])
    s = _score(d, "score.E.json", panel="E", focus=None, exclude_p=None)
    assert s["arms"]["VC"]["placed"] == 0 and s["precision"]["VC"] == dict(wrong=0, placed=0, rate=None)
    why = "VC placed no scored read"
    assert s["bars"]["E1"] == {"value": None, "num": 0, "den": 0, "limit": 0.02, "op": "<=", "pass": None, "reason": why}
    assert s["bars"]["E2"] == {"value": 0, "limit": 0, "op": ">", "pass": None, "reason": why}
    assert s["bars"]["K2"]["pass"] is True and s["verdict"] == "UNDERPOWERED"


def test_vc_placements_on_reads_that_are_not_scored_do_not_make_e1_and_e2_evaluable(tmp_path):
    d = toy_dir(tmp_path, [tied_read("u1", "F.0", ("unique", "F.0")), tied_read("t1", "F.0", ("own_tie", None), vc="F.1", cols=1)])
    s = _score(d, "score.E.json", panel="E", focus=None, exclude_p=None)
    assert s["vc"]["placements"] == 1 and s["arms"]["VC"]["placed"] == 0                     # t1: placed, but only K1 sees it
    assert s["bars"]["E1"]["pass"] is None and s["bars"]["E2"]["pass"] is None and s["verdict"] == "UNDERPOWERED"


def test_a2_and_a3_cannot_be_evaluated_when_the_focus_family_has_no_scored_read(tmp_path):
    d = toy_dir(tmp_path, [tied_read(f"t{i}", "F.0", ("own_tie", None)) for i in range(3)])      # own_tie only: K1 and nothing else
    s = _score(d, exclude_p=None)
    assert s["scored"] == 0 and s["K1"] == dict(reads=3, placed=0)
    why = "the focus family has no scored read"
    assert s["bars"]["A1"]["pass"] is None and s["bars"]["A1"]["reason"] == "VC placed no scored read of the focus family"
    assert (s["bars"]["A2"]["pass"], s["bars"]["A2"]["reason"]) == (None, why)
    assert (s["bars"]["A3"]["pass"], s["bars"]["A3"]["reason"]) == (None, why)
    assert s["bars"]["K1"]["pass"] is True and s["bars"]["K2"]["pass"] is True and s["verdict"] == "UNDERPOWERED"


def test_k2_cannot_be_evaluated_without_a_read_whose_copy_is_present_and_a_certain_wrong_placement_still_fails_e1(tmp_path):
    d = toy_dir(tmp_path, [tied_read("n1", "F.0", ("no_counterpart", None), vc="F.1", cols=1)])
    s = _score(d, "score.E.json", panel="E", focus=None, exclude_p=None)
    assert s["bars"]["K2"] == {"value": None, "num": 0, "den": 0, "limit": 0.02, "op": "<=", "pass": None,
                               "reason": "no read with a present copy"}
    assert s["K1"] == dict(reads=0, placed=0) and s["arms"]["VC"] == dict(placed=1, right=0, wrong=1, abstained=0, fixes=0,
                                                                          breaks=0, both_wrong=1)
    assert s["bars"]["E1"]["pass"] is False and s["verdict"] == "FAIL"                       # wrong 1 of 1: evaluable, failed


def test_a_bar_that_failed_still_fails_the_panel_when_another_bar_cannot_be_evaluated(tmp_path):
    # panel E: VC places nothing (E1 and E2 not evaluable), but half of the reads with a present copy are flagged: K2 failed
    d = toy_dir(tmp_path, [tied_read("u1", "F.0", ("unique", "F.0")), tied_read("u2", "F.0", ("unique", "F.1"), flag="F#0")])
    s = _score(d, "score.E.json", panel="E", focus=None, exclude_p=None)
    assert [s["bars"][b]["pass"] for b in ("E1", "E2", "K2")] == [None, None, False] and s["verdict"] == "FAIL"


def test_a_focus_family_with_scored_reads_where_vc_places_nothing_fails_a3_while_a1_and_a2_are_not_evaluable(tmp_path):
    # Task 6 ruling A3: none of these reads is placed wrong by minimap2, so no outcome could reach 10 fixes: A2 not evaluable
    d = toy_dir(tmp_path, [tied_read(f"u{i}", "F.0", ("unique", "F.0")) for i in range(3)])
    s = _score(d, exclude_p=None)
    assert [s["bars"][b]["pass"] for b in ("A1", "A2", "A3")] == [None, None, False] and s["verdict"] == "FAIL"


def _edit_calls(out, fn):
    rows = _rows(os.path.join(out, "calls.tsv"))
    fn(rows)
    _put(os.path.join(out, "calls.tsv"), CALLS_HEAD, [[r[k] for k in CALLS_HEAD.split()] for r in rows])


def test_inconsistent_or_unsupported_input_stops_the_scorer(tmp_path):
    d = make_score_toy(tmp_path)
    for mode in ("world", "units", "other"):
        with pytest.raises(ValueError, match="mode"):
            O.cmd_score(d.out, [d.truth], mode)
    with pytest.raises(ValueError, match="no truth"):
        O.cmd_score(d.out, [d.truth], "none")
    with pytest.raises(ValueError, match="truth"):
        O.cmd_score(d.out, [], "donor")
    world = str(tmp_path / "world_truth.tsv")                                                  # read, copy: --mode world's truth
    _put(world, "read copy", [["r1", "F.1"]])
    with pytest.raises(ValueError, match="class"):
        O.cmd_score(d.out, [world], "donor")
    with pytest.raises(ValueError, match="corrector"):
        O.cmd_score(d.out, [d.truth], "donor", n4_p=world)
    for panel in ("D", "Cg", "X"):
        with pytest.raises(ValueError, match="panel"):
            O.cmd_score(d.out, [d.truth], "donor", panel=panel)
    with pytest.raises(ValueError, match="focus"):
        O.cmd_score(d.out, [d.truth], "donor", panel="A")
    with pytest.raises(ValueError, match="focus"):
        O.cmd_score(d.out, [d.truth], "donor", panel="A", focus="G", exclude_p=d.par)
    with pytest.raises(ValueError, match="Z"):                                                  # a focus family absent from the run
        O.cmd_score(d.out, [d.truth], "donor", panel="A", focus="Z")
    for panel in (None, "C", "E"):                                                              # --focus belongs to --panel A
        with pytest.raises(ValueError, match="focus"):
            O.cmd_score(d.out, [d.truth], "donor", panel=panel, focus="F")
    with pytest.raises(ValueError, match="panel"):
        O.cmd_score(d.out, [], "none", panel="A", focus="F")

    def placed_and_flagged(rows):
        next(r for r in rows if r["read"] == "r5")["vc"] = "F.0"
    _edit_calls(d.out, placed_and_flagged)
    with pytest.raises(ValueError, match="r5.*flag"):
        O.cmd_score(d.out, [d.truth], "donor")
    d = make_score_toy(tmp_path)
    _edit_calls(d.out, lambda rows: rows.append(dict(rows[0])))
    with pytest.raises(ValueError, match="r1.*twice"):
        O.cmd_score(d.out, [d.truth], "donor")
    d = make_score_toy(tmp_path)
    _edit_calls(d.out, lambda rows: rows[0].update(outside="1"))
    with pytest.raises(ValueError, match="r1.*outside"):
        O.cmd_score(d.out, [d.truth], "donor")
    d = make_score_toy(tmp_path)
    with open(os.path.join(d.out, "F.reads.tsv"), "a") as fh:                                  # a read calls.tsv does not hold
        fh.write("\t".join(["zz", "F.0", "1", "1000", "2", "0.01", "900", "1", "0", "0", "0", "0", "0"]) + "\n")
    with pytest.raises(ValueError, match="zz.*calls.tsv"):
        O.cmd_score(d.out, [d.truth], "donor")
    d = make_score_toy(tmp_path)
    kept = [r for r in _rows(os.path.join(d.out, "F.reads.tsv")) if r["read"] != "r3"]         # a tied read without a reads.tsv row
    _put(os.path.join(d.out, "F.reads.tsv"), READS_HEAD, [[r[k] for k in READS_HEAD.split()] for r in kept])
    with pytest.raises(ValueError, match="r3.*no row"):
        O.cmd_score(d.out, [d.truth], "donor")
    d = make_score_toy(tmp_path)                                                                # unique, but no truth copy
    rows = _rows(d.truth)
    next(r for r in rows if r["read"] == "r1")["truth"] = ""
    _put(d.truth, TRUTH_HEAD, [[r[k] for k in TRUTH_HEAD.split()] for r in rows])
    with pytest.raises(ValueError, match="r1.*truth copy"):
        O.cmd_score(d.out, [d.truth], "donor")


def test_a_run_with_no_tied_read_scores_to_nothing_and_is_underpowered(tmp_path):
    d = toy_dir(tmp_path, [])
    s = _score(d, "score.E.json", panel="E", focus=None, exclude_p=None)
    assert s["tied"] == 0 and s["scored"] == 0 and s["families"] == {} and s["arms"]["VC"]["placed"] == 0 and s["flags"] == []
    assert {b: v["pass"] for b, v in s["bars"].items()} == dict(E1=None, E2=None, K2=None) and s["verdict"] == "UNDERPOWERED"
    with pytest.raises(ValueError, match="F"):                                                  # panel A: no family F in the run
        O.cmd_score(d.out, [d.truth], "donor", panel="A", focus="F")


def random_reads(rng, n, fams=("F", "G", "H", "P")):
    """n random tied reads over `fams` (2 or 3 copies each): every field the scorer reads drawn independently, a flagged read
    never placed (what run guarantees), minimap2's primary now and then outside the tied copies."""
    reads = []
    for i in range(n):
        f = rng.choice(fams)
        cands = [f"{f}.{k}" for k in range(rng.choice((2, 3)))]
        pool = cands + [""]
        flag = rng.choice(("", "", "", f"{f}#0", f"{f}#1", f"{f}#0,{f}#1"))
        cls = rng.choice(("unique", "unique", "unique", "unique_elsewhere", "own_tie", "no_counterpart", "elsewhere", "poor_match",
                          None))
        truth = None if cls is None else (cls, rng.choice(cands + [f"{f}.9"]) if cls in ("unique", "unique_elsewhere") else None)
        kinds = rng.sample(("separating_call", "separating_ref", "erased_ref", "foreign"), rng.randrange(0, 3))
        reads.append(tied_read(f"q{i:05d}", rng.choice(pool + [f"X.{rng.randrange(2)}"]), truth, fam=f, cands=cands,
                               vc="" if flag else rng.choice(pool), n3=rng.choice(pool), flag=flag,
                               mapq=rng.choice((0, 0, 1, 7, 60)), cols=rng.randrange(1, 4), kinds=kinds))
    return reads


def recount(reads, excluded):
    """The pooled counts of `reads` outside `excluded`, straight from the amended spec (a naive loop per read, no shared code):
    a tied read is SCORED when its class is unique, unique_elsewhere, elsewhere or no_counterpart; a placement is right only
    when the placed copy is the truth copy, which only unique and unique_elsewhere reads have."""
    zero = dict(placed=0, right=0, wrong=0, abstained=0, fixes=0, breaks=0, both_wrong=0)
    names = ("N1", "N2", "N3", "VC")
    arms, dec_arms = {a: dict(zero) for a in names}, {a: dict(zero) for a in names}
    k1, k2, d2, flags, classes = [0, 0], [0, 0], [0, 0, 0], {}, {}
    n_scored = n_dec = n_outside = n_not_scored_placed = 0
    not_scored = {}

    def tally(c, x, tc, n1):
        if x == "":
            c["abstained"] += 1
            return
        c["placed"] += 1
        if tc is not None and x == tc:
            c["right"] += 1
            c["fixes"] += 1 if n1 != tc else 0
        else:
            c["wrong"] += 1
            if tc is not None and n1 == tc:
                c["breaks"] += 1
            else:
                c["both_wrong"] += 1

    for t in reads:
        if t["fam"] in excluded:
            continue
        cls, tc = t["truth"] if t["truth"] is not None else ("no_truth_row", None)
        n1 = t["mm"] if t["mm"] else "(off)"
        got = {"N1": n1, "N2": n1 if t["mapq"] >= 1 else "", "N3": t["n3"], "VC": t["vc"]}
        c = classes.setdefault(cls, dict(reads=0, placed=0, flagged=0))
        c["reads"] += 1
        c["placed"] += 1 if t["vc"] else 0
        c["flagged"] += 1 if t["flag"] else 0
        if cls == "own_tie":
            k1[0] += 1
            k1[1] += 1 if t["vc"] else 0
        if cls in ("unique", "unique_elsewhere", "own_tie"):
            k2[0] += 1
            k2[1] += 1 if t["flag"] else 0
        for f in t["flag"].split(","):
            if f:
                flags.setdefault(f, {})
                flags[f][cls] = flags[f].get(cls, 0) + 1
        if cls not in ("unique", "unique_elsewhere", "elsewhere", "no_counterpart"):
            if t["vc"]:
                n_not_scored_placed += 1
                not_scored[cls] = not_scored.get(cls, 0) + 1
            continue
        n_scored += 1
        n_outside += 1 if t["mm"] not in t["cands"] else 0
        for a, x in got.items():
            tally(arms[a], x, tc, n1)
        if tc is not None and tc in t["cands"]:                              # the decidable subset
            n_dec += 1
            for a, x in got.items():
                tally(dec_arms[a], x, tc, n1)
        if tc is not None and n1 == tc:
            d2[0] += 1
            if got["VC"] != "" and got["VC"] != n1:
                d2[1] += 1
            if t["flag"]:
                d2[2] += 1
    return dict(arms=arms, dec_arms=dec_arms, K1=dict(reads=k1[0], placed=k1[1]), K2=dict(reads=k2[0], flagged=k2[1]),
                scored=n_scored, decidable=n_dec, outside_scored=n_outside, classes=dict(sorted(classes.items())),
                D2=dict(reads=d2[0], moved=d2[1], flagged=d2[2]),
                C1=dict(flags=len(flags), counted=sum(1 for c in flags.values() if c.get("no_counterpart", 0) >= 2)),
                flags={f: dict(sorted(c.items())) for f, c in sorted(flags.items())},
                not_scored=dict(placements=n_not_scored_placed, classes=dict(sorted(not_scored.items()))))


@pytest.mark.parametrize("seed", [1, 2, 3, 4, 5, 6])
def test_the_score_matches_a_direct_recount_on_random_runs(tmp_path, seed):
    rng = random.Random(seed)
    reads = random_reads(rng, 400)
    d = toy_dir(tmp_path, reads, par=["P"])
    s = _score(d, "score.json", panel=None, focus=None)
    want = recount(reads, {"P"})
    assert {a: s["arms"][a] for a in want["arms"]} == want["arms"]
    assert (s["K1"], s["K2"], s["D2"], s["C1"]) == (want["K1"], want["K2"], want["D2"], want["C1"])
    got = (s["scored"], s["decidable"]["reads"], s["outside"]["scored"])
    assert got == (want["scored"], want["decidable"], want["outside_scored"])
    assert {a: s["decidable"]["arms"][a] for a in want["dec_arms"]} == want["dec_arms"]
    assert s["classes"] == want["classes"] and s["vc"]["not_scored"] == want["not_scored"]
    assert {f["id"]: f["classes"] for f in s["flags"]} == want["flags"]
    assert sorted(s["families"]) == ["F", "G", "H"] and list(s["excluded"]) == ["P"]
    for fam in ("F", "G", "H"):                                    # a family alone is the recount of the reads outside the others
        alone = recount(reads, {x for x in ("F", "G", "H", "P") if x != fam})
        assert {a: s["families"][fam]["arms"][a] for a in alone["arms"]} == alone["arms"]
    p = recount(reads, {"F", "G", "H"})
    assert {a: s["excluded"]["P"]["arms"][a] for a in p["arms"]} == p["arms"] and s["excluded"]["P"]["K1"] == p["K1"]


def test_the_cli_runs_cmd_score_with_its_options(monkeypatch):
    seen = []
    monkeypatch.setattr(O, "cmd_score", lambda *a, **k: seen.append((a, k)) or 0)
    assert O.main(["score", "out", "t1.tsv", "t2.tsv", "--mode", "donor", "--n4", "n4.tsv", "--panel", "A", "--focus", "F",
                   "--exclude", "par.json"]) == 0
    assert O.main(["score", "out", "--mode", "none"]) == 0
    assert seen == [(("out", ["t1.tsv", "t2.tsv"], "donor"), dict(n4_p="n4.tsv", panel="A", focus="F", exclude_p="par.json")),
                    (("out", [], "none"), dict(n4_p=None, panel=None, focus=None, exclude_p=None))]
    with pytest.raises(SystemExit):
        O.main(["score", "out", "t.tsv"])                                                    # --mode is required
    with pytest.raises(SystemExit):
        O.main(["score"])
    with pytest.raises(SystemExit):
        O.main(["score", "out", "--mode"])                                                   # an option with no value
    with pytest.raises(SystemExit):
        O.main(["score", "out", "--mode", "none", "--threads", "3"])                         # run's option on score
    with pytest.raises(SystemExit):
        O.main(["run", "b.bam", "c.tsv", "g.fa", "out", "--mode", "donor"])                  # score's option on run


def test_the_scoring_is_the_same_under_any_hash_seed(tmp_path):
    d = make_score_toy(tmp_path)
    n4 = str(tmp_path / "o2correct.tsv")
    _put(n4, "read family corrector minimap2 reason candidates margin scores",
         [["r1", "F", "F.1", "F.0", "sample", "F.0,F.1", 5, ""]])
    seen = []
    for seed in ("0", "1", "12345"):
        r = subprocess.run([sys.executable, os.path.join(HERE, "o2vc.py"), "score", d.out, d.truth, "--mode", "donor", "--n4", n4,
                            "--panel", "A", "--focus", "F", "--exclude", d.par], env=dict(os.environ, PYTHONHASHSEED=seed),
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        seen.append([open(os.path.join(d.out, f), "rb").read() for f in ("score.A.F.json", "F.case.tsv", "F.summary.md",
                                                                           "G.case.tsv", "G.summary.md")] + [r.stdout])
    assert seen[0] == seen[1] == seen[2]


# ---- Task 6 rulings (task-6-rulings.md): the weak-het policy, the exit code, bar evaluability, the DAZ copy table ----

@pytest.mark.parametrize("af, policy, vc, n3, weak", [
    ("0.14", "keep", "", "F.0", 1),          # P1: F.1's weak het {C, A} shares F.0's A, the reference PSV is erased
    ("0.14", "unknown", "", "", 1),          # P2: F.1 unknown at the column, under N3 too (N3 knows only what VC knows)
    ("0.14", "ref", "F.0", "F.0", 1),        # P3: the weak alt dropped, F.1 is {C}: the PSV separates and Z's A votes F.0
    ("0.15", "ref", "", "F.0", 0),           # AF 0.15 is not weak: P3 changes nothing
])
def test_the_weak_het_policy_decides_a_column_whose_reference_psv_a_weak_het_erases(tmp_path, monkeypatch, af, policy, vc,
                                                                                     n3, weak):
    A = family_refs((300,))["chrA"]                                                   # chrB has C at 300, chrA A
    w = make_family_world(tmp_path, psv=(300,), extra=[rec("Z", "chrA", 100, 590, seq=A[100:700]),
                                                       rec("Z", "chrB", 100, 590, secondary=True)])
    body = {"chrB": f"chrB\t301\t.\tC\tA\t.\tLowQual\t.\tGT:GQ:DP:AF\t0/1:30:8:{af}\n"}
    monkeypatch.setattr(O.D, "run_longcallr", fake_vcf_caller([], body))
    out = tmp_path / "out"
    assert O.cmd_run(w.bam, w.copies, w.genome, str(out), exclude_chrM=True, weak_het=policy) == 0
    z = {r["read"]: r for r in _rows(out / "calls.tsv")}["Z"]
    assert (z["vc"], z["n3"]) == (vc, n3)
    assert json.load(open(out / "run.json"))["weak_het"] == dict(policy=policy, threshold=0.15, records=weak)


def test_run_exits_non_zero_when_every_attempted_longcallr_call_failed(tmp_path, monkeypatch):
    w, prefixes, out = make_family_world(tmp_path), [], tmp_path / "out"

    def fail(teacher_bam, genome, out_prefix, binary=None, threads=2, region=None):
        prefixes.append(out_prefix)
        raise RuntimeError(f"longcallR exited 101 ({region})")
    monkeypatch.setattr(O.D, "run_longcallr", fail)
    assert O.cmd_run(w.bam, w.copies, w.genome, str(out), exclude_chrM=True) == 1
    assert len(prefixes) == 2 and sorted(json.load(open(out / "run.json"))["failed_copies"]) == ["F.0", "F.1"]


def test_a2_is_not_evaluable_when_the_focus_family_has_fewer_than_ten_reads_minimap2_places_wrong(tmp_path):
    os.makedirs(tmp_path / "a")
    os.makedirs(tmp_path / "b")
    s = _score(panel_a_toy(tmp_path / "a", fix=9, ok=39), exclude_p=None)
    assert s["families"]["F"]["arms"]["N1"]["wrong"] == 9
    assert s["bars"]["A2"] == {"value": 9, "limit": 10, "op": ">=", "pass": None,
                               "reason": "the focus family has 9 scored reads minimap2 places wrong, fewer than 10"}
    assert s["verdict"] == "UNDERPOWERED"
    s = _score(panel_a_toy(tmp_path / "b", fix=9, ok=38, miss=1), exclude_p=None)       # reachable and missed: FAIL
    assert s["bars"]["A2"] == {"value": 9, "limit": 10, "op": ">=", "pass": False} and s["verdict"] == "FAIL"


def test_c1_is_not_evaluable_with_fewer_than_two_no_counterpart_tied_reads(tmp_path):
    d = toy_dir(tmp_path, [tied_read("n1", "F.0", ("no_counterpart", None), flag="F#0"),
                           tied_read("u1", "F.0", ("unique", "F.0"), flag="F#0")])
    s = _score(d, "score.C.json", panel="C", focus=None, exclude_p=None)
    assert s["bars"]["C1"] == {"value": 0, "limit": 1, "op": ">=", "pass": None,
                               "reason": "1 no_counterpart tied read in the pool, fewer than 2"}
    assert s["verdict"] == "UNDERPOWERED"


def test_copies_from_sets_writes_pool_rows_with_the_union_of_transcript_exons(tmp_path):
    p = tmp_path / "f.json"
    p.write_text(json.dumps({"sets": {"seq19": [{"family": "S19X", "copies": [
        {"chrom": "chrY", "start": 10, "end": 90, "strand": "-", "transcripts": {"t1": [[10, 20], [50, 60]], "t2": [[10, 20], [70, 90]]}}]}]}}))
    assert O.copies_from_sets(str(p), "seq19", "S19X") == [dict(copy="S19X.0", family="S19X", chrom="chrY", start=10, end=90,
                                                                strand="-", exons="10-20,50-60,70-90")]


def test_copies_from_sets_reads_the_family_file_s_own_forms_and_the_cli_writes_a_table_load_copies_reads(tmp_path):
    # families19cal_full2_all.json writes start and end as strings and the transcripts as a Python literal
    p, out = tmp_path / "f.json", tmp_path / "c.tsv"
    p.write_text(json.dumps({"sets": {"seq19": [
        {"family": "S19Y", "copies": [{"chrom": "chr1", "start": "5", "end": "40", "strand": "+", "transcripts": "{'t': [[5, 40]]}"}]},
        {"family": "S19X", "copies": [
            {"chrom": "chrY", "start": "100", "end": "300", "strand": "+", "transcripts": "{'a': [[200, 300], [100, 150]], 'b': [[100, 150]]}"},
            {"chrom": "chrY", "start": "1000", "end": "1100", "strand": "-", "transcripts": "{'c': [[1000, 1100]]}"}]}]}}))
    assert O.main(["copies-from-sets", str(p), "seq19", "S19X", str(out)]) == 0
    assert [(r["copy"], r["start"], r["end"], r["exons"]) for r in _rows(out)] == [
        ("S19X.0", "100", "300", "100-150,200-300"), ("S19X.1", "1000", "1100", "1000-1100")]
    assert open(out).readline().split() == "copy family chrom start end strand exons".split()
    assert sorted(c for v in O.P.load_copies(str(out)).values() for c, *_ in v) == ["S19X.0", "S19X.1"]
    with pytest.raises(ValueError, match="S19Z"):
        O.copies_from_sets(str(p), "seq19", "S19Z")
