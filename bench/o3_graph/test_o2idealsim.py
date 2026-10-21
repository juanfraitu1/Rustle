"""PREREG_o2_corrector_idealsim_2026-10-21: the simulation's copies, instances, truth strata and bars."""
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_o2idealsim", os.path.join(HERE, "o2idealsim.py"))
IS = importlib.util.module_from_spec(spec)
spec.loader.exec_module(IS)


def _copy(name, chrom="NC_1", start=100, end=400, strand="+", exons=((100, 150), (300, 400))):
    return dict(copy=name, family=name.rsplit(".", 1)[0], chrom=chrom, start=start, end=end, strand=strand, exons=list(exons))


def test_the_family_set_lists_each_family_s_copies_by_catalog_index_with_one_coding_model():
    rows = [_copy("F.1", start=900, end=1000, exons=((900, 1000),)), _copy("F.0"), _copy("G.0", strand="-")]
    fs = IS.family_set(rows)
    (f, g) = fs["sets"]["ideal"]
    assert f["family"] == "F" and [c["start"] for c in f["copies"]] == [100, 900] and g["family"] == "G"
    assert f["copies"][0] == dict(chrom="NC_1", start=100, end=400, strand="+", transcripts={"main": [[100, 150], [300, 400]]},
                                  main="main", biotype="protein_coding")


def test_every_copy_has_a_reference_instance_and_its_counterpart_flagged_when_contested():
    rows = [_copy("F.0"), _copy("F.1", start=900, end=1000, exons=((900, 1000),))]
    orth = [dict(copy="F.0", hap="mat", contig="CM1", start=5100, end=5400, strand="-", mapq=60, exons=[[5100, 5200], [5300, 5400]]),
            dict(copy="F.1", hap="mat", contig="CM1", start=7000, end=7100, strand="+", mapq=60, exons=[[7000, 7100]])]
    inst = IS.instances(rows, orth, {("F.1", "mat")})
    assert [(i["copy"], i["hap"], i["contig"], i["strand"], i["contested"]) for i in inst] == [
        ("F.0", "ref", "NC_1", "+", 0), ("F.1", "ref", "NC_1", "+", 0), ("F.0", "mat", "CM1", "-", 0), ("F.1", "mat", "CM1", "+", 1)]
    assert inst[0]["isoforms"] == {"main": [[100, 150], [300, 400]]} and inst[2]["isoforms"] == {"main": [[5100, 5200], [5300, 5400]]}


def test_truth_strata_keep_reference_reads_and_uncontested_counterpart_reads():
    t = [dict(read="a", truth="F.0", hap="ref"), dict(read="b", truth="F.0", hap="mat"), dict(read="c", truth="F.1", hap="mat"),
         dict(read="d", truth="F.1", hap="ref")]
    s = IS.strata(t, {("F.1", "mat")})
    assert [r["read"] for r in s["ref"]] == ["a", "d"] and [r["read"] for r in s["alt"]] == ["b"]
    assert sorted(r["read"] for r in s["all"]) == ["a", "b", "d"]


def _score(reads, mm_right, fixes, breaks, better=(), worse=()):
    return dict(pooled={"reads": reads, "minimap2 right": mm_right, "fixes": fixes, "breaks": breaks},
                families_better=list(better), families_worse=list(worse))


def test_a_replicate_passes_only_when_every_bar_holds():
    ok = IS.replicate_bars(_score(200, 150, 30, 2, ("F", "G")), _score(120, 120, 0, 1))
    assert ok == dict(P=True, R=True, S=True, H=True, F=True, verdict="PASS")
    assert IS.replicate_bars(_score(200, 150, 30, 7, ("F", "G")), _score(120, 120, 0, 1))["verdict"] == "FAIL"   # S: 30 < 35
    assert IS.replicate_bars(_score(200, 150, 30, 2, ("F",)), _score(120, 120, 0, 1))["F"] is False
    assert IS.replicate_bars(_score(200, 150, 30, 2, ("F", "G"), ("H",)), _score(120, 120, 0, 1))["F"] is False
    assert IS.replicate_bars(_score(200, 150, 30, 2, ("F", "G")), _score(120, 117, 0, 3))["H"] is False
    assert IS.replicate_bars(_score(200, 192, 6, 0, ("F", "G")), _score(120, 120, 0, 0))["verdict"] == "UNDERPOWERED"
    assert IS.replicate_bars(_score(200, 150, 1, 0, ("F", "G")), _score(120, 120, 0, 0))["R"] is True          # 1 >= 2% of 50
    assert IS.replicate_bars(_score(200, 100, 1, 0, ("F", "G")), _score(120, 120, 0, 0))["R"] is False         # 1 < 2% of 100


def test_the_overall_verdict_needs_both_replicates_to_pass():
    assert IS.overall(["PASS", "PASS"]) == "PASS"
    assert IS.overall(["PASS", "FAIL"]) == "FAIL" and IS.overall(["UNDERPOWERED", "FAIL"]) == "FAIL"
    assert IS.overall(["PASS", "UNDERPOWERED"]) == "UNDERPOWERED"
