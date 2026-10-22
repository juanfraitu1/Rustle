import csv
import importlib.util
import json
import os
import random

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
    assert seen == dict(a=("b.bam", "c.tsv", "g.fa", "out"), k=dict(families={"F", "G"}, exclude_chrM=True, threads=3))
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
