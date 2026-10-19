import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


F = _load("fastset")


def test_reads_shared_by_worlds_are_mapped_once_by_sequence():
    worlds = {"INTACT": {"s1": "ACGT", "s2": "GGGG"}, "OUTSIDE": {"s9": "ACGT", "s3": "GGGG"}, "CROSS": {"s1": "ACGT", "s7": "TTTT"}}
    union, owners = F.union_reads(worlds)
    assert sorted(union.values()) == ["ACGT", "GGGG", "TTTT"]
    by_seq = {union[u]: sorted(v) for u, v in owners.items()}
    assert by_seq["ACGT"] == [("CROSS", "s1"), ("INTACT", "s1"), ("OUTSIDE", "s9")]
    assert by_seq["TTTT"] == [("CROSS", "s7")]


def test_a_union_alignment_is_handed_to_every_world_read_with_that_sequence_under_its_own_name():
    owners = {"u0": [("INTACT", "s1"), ("OUTSIDE", "s9")], "u1": [("CROSS", "s7")]}
    lines = ["u0\t0\tchr1\t11\t60\t4M\t*\t0\t0\tACGT\t*\tNM:i:0\n", "u0\t256\tchr2\t5\t0\t4M\t*\t0\t0\t*\t*\tNM:i:1\n",
             "u1\t16\tchr3\t7\t60\t4M\t*\t0\t0\tAAAA\t*\n"]
    out = F.split_union(lines, owners)
    assert [ln.split("\t")[:3] for ln in out["INTACT"]] == [["s1", "0", "chr1"], ["s1", "256", "chr2"]]
    assert [ln.split("\t")[0] for ln in out["OUTSIDE"]] == ["s9", "s9"] and out["OUTSIDE"][0].endswith("NM:i:0\n")
    assert [ln.split("\t")[0] for ln in out["CROSS"]] == ["s7"]


def test_a_local_reference_is_reused_only_for_the_family_set_it_was_built_for(tmp_path):
    fj = tmp_path / "f.json"
    fj.write_text('{"sets": {"s1": [1], "s2": [2]}}')
    F.write_local_source(str(tmp_path), str(fj), "s1")
    assert F.local_matches(str(tmp_path), str(fj), "s1") and not F.local_matches(str(tmp_path), str(fj), "s2")
    fj.write_text('{"sets": {"s1": [9], "s2": [2]}}')
    assert not F.local_matches(str(tmp_path), str(fj), "s1")


def _opts(**kw):
    import argparse
    base = dict(set_threads=4, threads=2, heavy_current=False, o2_distance="gapcomp", member_genome_check=True, all_arms=False,
                o3_procs=2, o2_genome="off", per_read_detectability=False, member_attribution=False, union_exact=False,
                o2_procs=2, member_attribution_delta=False, twins=False, error_floor=False, o1_hosts=False,
                mapped_sites=False, branch_genome_novel=False, allele_discount=False,
                mosaic_shared=False, o3_pool="paths", haplotype_check=False,
                haplotype_check_genome=False, branch_genome_segments=False,
                branch_linear=False, aln_cache=False, twin_abstain=False, veto_abstain=False, o2_significance=False,
                locus_pairwise=False, structure_check=False, struct_rescue=False, short_path_align=False, twin_guard=False,
                ref_copies=False, local=None, ref_carried=False)
    base.update(kw)
    return argparse.Namespace(**base)


def test_a_world_runs_current_then_o1_o2_and_o3_and_the_genome_step_only_when_asked():
    steps = [s for s, _, _ in F.chain_steps(_opts(), "/R", "INTACT")]
    assert steps == ["current", "O1", "O2", "O3"]
    for mode, method in (("fallback", "integl"), ("combined", "integc")):
        steps = F.chain_steps(_opts(o2_genome=mode), "/R", "INTACT")
        assert [s for s, _, _ in steps] == ["current", "O1", "O2", "O2_genome", "O3"]
        assert steps[3][1] == f"python3 bench/o3_graph/o2graph.py run /R/INTACT {method} 2" and steps[3][2] == "light"


def test_the_set_is_scored_on_the_genome_step_s_table_and_per_read_when_asked():
    assert F.e2e_command(_opts(), "/R", ["A", "B"]) == "python3 bench/o3_graph/summary.py e2e /R/e2e.json /R/A /R/B"
    assert F.e2e_command(_opts(o2_genome="fallback", per_read_detectability=True), "/R", ["A"]) == (
        "python3 bench/o3_graph/summary.py e2e /R/e2e.json /R/A --o2 o2.integl.m2.tsv --per-read-detectability")
    assert F.e2e_command(_opts(o2_genome="combined"), "/R", ["A"]) == (
        "python3 bench/o3_graph/summary.py e2e /R/e2e.json /R/A --o2 o2.integc.m2.tsv")


def test_o3_attributes_its_call_members_only_when_asked():
    assert "--member-attribution" not in F.chain_steps(_opts(), "/R", "INTACT")[-1][1]
    assert F.chain_steps(_opts(member_attribution=True), "/R", "INTACT")[-1][1].endswith(
        "detect.py /R/INTACT --member-genome-check --member-attribution --arms I")


def test_o2_takes_the_exact_union_only_when_asked():
    assert F.chain_steps(_opts(), "/R", "INTACT")[2][1] == ("env O2_PROCS=2 O2_THREADS=2 python3 bench/o3_graph/o2graph.py run /R/INTACT integ 2 "
                                                         "--distance gapcomp")
    assert F.chain_steps(_opts(union_exact=True), "/R", "INTACT")[2][1].endswith("integ 2 --distance gapcomp --union-exact")


def test_o3_scales_attribution_by_delta_only_when_asked():
    assert "--member-attribution-delta" not in F.chain_steps(_opts(), "/R", "INTACT")[-1][1]
    assert "--member-attribution-delta --arms I" in F.chain_steps(_opts(member_attribution_delta=True), "/R", "INTACT")[-1][1]


def test_o1_marks_twins_only_when_asked():
    assert F.chain_steps(_opts(), "/R", "INTACT")[1][1] == "python3 bench/o3_graph/o1seed.py /R/INTACT"
    assert F.chain_steps(_opts(twins=True), "/R", "INTACT")[1][1] == "python3 bench/o3_graph/o1seed.py /R/INTACT --twins"


def test_o3_subtracts_the_error_floor_only_when_asked():
    assert "--error-floor" not in F.chain_steps(_opts(), "/R", "INTACT")[-1][1]
    assert "--error-floor --arms I" in F.chain_steps(_opts(error_floor=True), "/R", "INTACT")[-1][1]


def test_o3_takes_its_hosts_from_o1_only_when_asked():
    assert "--o1-hosts" not in F.chain_steps(_opts(), "/R", "INTACT")[-1][1]
    assert "--o1-hosts --arms I" in F.chain_steps(_opts(o1_hosts=True), "/R", "INTACT")[-1][1]


def test_the_error_floor_reaches_o2_o3_and_the_scorer():
    steps = dict((n, c) for n, c, _ in F.chain_steps(_opts(error_floor=True, o2_genome="fallback"), "/R", "INTACT"))
    assert steps["O2"].endswith("--error-floor") and steps["O2_genome"].endswith("integl 2 --error-floor")
    assert "--error-floor" in steps["O3"]
    assert F.e2e_command(_opts(error_floor=True), "/R", ["A"]).endswith("/R/A --error-floor")


def test_o3_clusters_on_mapped_sites_only_when_asked():
    assert "--mapped-sites" not in F.chain_steps(_opts(), "/R", "INTACT")[-1][1]
    assert "--mapped-sites --arms I" in F.chain_steps(_opts(mapped_sites=True), "/R", "INTACT")[-1][1]


def test_a_killed_step_s_empty_timing_file_reads_as_unknown():
    # round 12: rlock's timeout killed /usr/bin/time before it wrote; the empty file crashed the driver
    assert F.parse_time("12.30 4567\n") == (12.3, 4567)
    assert F.parse_time("Command terminated by signal 9\n590.01 2048\n") == (590.01, 2048)
    wall, rss = F.parse_time("")
    assert wall != wall and rss == 0                                   # nan: the step left no time


def test_a_step_s_timing_carries_its_cpu_seconds():
    # round 15's cost bar compared wall times at unequal threads (the current method's minimap2 at 8, each integrated step
    # at 2): every step now records its CPU seconds (user + system, its children's included) beside the wall time
    assert F.parse_time("12.30 4567 20.50 1.25\n") == (12.3, 4567)
    assert F.parse_cpu("12.30 4567 20.50 1.25\n") == 21.75
    assert F.parse_cpu("Command terminated by signal 9\n590.01 2048 900.00 3.00\n") == 903.0
    old, empty = F.parse_cpu("12.30 4567\n"), F.parse_cpu("")
    assert old != old and empty != empty                               # nan: no CPU time recorded


def test_the_branch_arm_filters_by_genome_novelty_only_when_asked():
    assert "--branch-genome-novel" not in F.chain_steps(_opts(), "/R", "INTACT")[-1][1]
    assert "--branch-genome-novel --arms I" in F.chain_steps(_opts(branch_genome_novel=True), "/R", "INTACT")[-1][1]


def test_o3_discounts_the_donor_s_alleles_only_when_asked():
    assert "--allele-discount" not in F.chain_steps(_opts(), "/R", "INTACT")[-1][1]
    assert "--allele-discount --arms I" in F.chain_steps(_opts(allele_discount=True), "/R", "INTACT")[-1][1]


def test_the_veto_genome_step_writes_and_scores_its_own_table():
    steps = dict((n, c) for n, c, _ in F.chain_steps(_opts(o2_genome="veto", error_floor=True), "/R", "INTACT"))
    assert steps["O2_genome"].endswith("integv 2 --error-floor")
    assert "--o2 o2.integv.m2.tsv" in F.e2e_command(_opts(o2_genome="veto"), "/R", ["A"])


def test_o3_takes_the_shared_position_mosaic_distance_and_the_final_o2_pool_only_when_asked():
    o3 = F.chain_steps(_opts(), "/R", "INTACT")[-1][1]
    assert "--mosaic-shared" not in o3 and "--pool-table" not in o3
    o3 = F.chain_steps(_opts(mosaic_shared=True, o3_pool="final", o2_genome="veto"), "/R", "INTACT")[-1][1]
    assert "--mosaic-shared" in o3 and "--pool-table o2.integv.m2.tsv" in o3


def test_o3_checks_a_call_s_haplotype_only_when_asked():
    assert "--haplotype-check" not in F.chain_steps(_opts(), "/R", "INTACT")[-1][1]
    assert "--haplotype-check" in F.chain_steps(_opts(haplotype_check=True), "/R", "INTACT")[-1][1]


def test_o2_and_o3_minimap2_runs_get_the_world_s_threads():
    steps = dict((n, c) for n, c, _ in F.chain_steps(_opts(threads=3), "/R", "INTACT"))
    assert "O2_THREADS=3" in steps["O2"] and "O3_THREADS=3" in steps["O3"]


def test_o3_checks_the_members_haplotype_on_their_own_loci_only_when_asked():
    assert "--haplotype-check-genome" not in F.chain_steps(_opts(), "/R", "INTACT")[-1][1]
    assert "--haplotype-check-genome" in F.chain_steps(_opts(haplotype_check_genome=True), "/R", "INTACT")[-1][1]
    assert "--branch-genome-segments" in F.chain_steps(_opts(branch_genome_segments=True), "/R", "INTACT")[-1][1]
    assert "--branch-genome-segments" not in F.chain_steps(_opts(), "/R", "INTACT")[-1][1]
    assert "--branch-linear" in F.chain_steps(_opts(branch_linear=True), "/R", "INTACT")[-1][1]
    st = dict((k, c) for k, c, _ in F.chain_steps(_opts(aln_cache=True), "/R", "INTACT"))
    assert "O2_ALN_CACHE=/R/.o2_aln_cache" in st["O2"] and "O3_ALN_CACHE=/R/.o3_aln_cache" in st["O3"]
    st0 = dict((k, c) for k, c, _ in F.chain_steps(_opts(), "/R", "INTACT"))
    assert "ALN_CACHE" not in st0["O2"] + st0["O3"]
    st = dict((k, c) for k, c, _ in F.chain_steps(_opts(twin_abstain=True, o2_genome="guard"), "/R", "INTACT"))
    assert st["O2"].endswith("--twin-abstain") and "integt 2" in st["O2_genome"] and "--twin-abstain" in st["O2_genome"]
    assert "--twin-abstain" not in st0["O2"]


def test_o2_s_genome_step_abstains_contradicted_and_insignificant_path_calls_only_when_asked():
    st0 = dict((k, c) for k, c, _ in F.chain_steps(_opts(o2_genome="veto"), "/R", "INTACT"))
    assert "--veto-abstain" not in st0["O2_genome"] and "--o2-significance" not in st0["O2"] + st0["O2_genome"]
    st = dict((k, c) for k, c, _ in F.chain_steps(_opts(o2_genome="veto", veto_abstain=True), "/R", "INTACT"))
    assert "integv 2" in st["O2_genome"] and "--veto-abstain" in st["O2_genome"] and "--veto-abstain" not in st["O2"]
    st = dict((k, c) for k, c, _ in F.chain_steps(_opts(o2_genome="veto", o2_significance=True), "/R", "INTACT"))
    assert "--o2-significance" in st["O2"] and "--o2-significance" in st["O2_genome"]    # integ weighs, the genome step gates


def test_o2_s_genome_step_takes_the_pairwise_locus_rule_and_the_structure_check_only_when_asked():
    st0 = dict((k, c) for k, c, _ in F.chain_steps(_opts(o2_genome="veto"), "/R", "INTACT"))
    assert "--locus-pairwise" not in st0["O2_genome"] and "--structure-check" not in st0["O2_genome"]
    st = dict((k, c) for k, c, _ in F.chain_steps(_opts(o2_genome="veto", locus_pairwise=True, structure_check=True), "/R",
                                                   "INTACT"))
    assert "--locus-pairwise" in st["O2_genome"] and "--structure-check" in st["O2_genome"]
    assert "--locus-pairwise" not in st["O2"] and "--structure-check" not in st["O2"]


def test_o2_s_genome_step_takes_the_structure_aware_rescue_and_short_read_realignment_only_when_asked():
    st0 = dict((k, c) for k, c, _ in F.chain_steps(_opts(o2_genome="veto"), "/R", "INTACT"))
    assert "--struct-rescue" not in st0["O2_genome"] and "--short-path-align" not in st0["O2_genome"]
    st = dict((k, c) for k, c, _ in F.chain_steps(_opts(o2_genome="veto", struct_rescue=True, short_path_align=True,
                                                         aln_cache=True), "/R", "INTACT"))
    assert "--struct-rescue" in st["O2_genome"] and "--short-path-align" in st["O2_genome"]
    assert "O2_ALN_CACHE=/R/.o2_aln_cache" in st["O2_genome"]          # the rescue reads the path alignments exactly


def test_o2_s_genome_step_takes_the_twin_guard_only_when_asked():
    st0 = dict((k, c) for k, c, _ in F.chain_steps(_opts(o2_genome="veto"), "/R", "INTACT"))
    st = dict((k, c) for k, c, _ in F.chain_steps(_opts(o2_genome="veto", twin_guard=True), "/R", "INTACT"))
    assert "--twin-guard" not in st0["O2_genome"] and "--twin-guard" in st["O2_genome"] and "integv 2" in st["O2_genome"]


def test_o1_searches_each_world_s_own_local_reference_for_copies_the_catalogue_lacks_only_when_asked():
    st0 = dict((k, c) for k, c, _ in F.chain_steps(_opts(), "/R", "INTACT"))
    assert "--ref-copies" not in st0["O1"]
    st = dict((k, c) for k, c, _ in F.chain_steps(_opts(ref_copies=True, local="/L"), "/R", "INTACT"))
    assert st["O1"].endswith("--ref-copies --local /L")
    st = dict((k, c) for k, c, _ in F.chain_steps(_opts(ref_copies=True, local="/L"), "/R", "REMOVED"))
    assert st["O1"].endswith("--ref-copies --local /R/local_masked")       # the removed copy is masked there
    st = dict((k, c) for k, c, _ in F.chain_steps(_opts(ref_copies=True), "/R", "OUTSIDE"))
    assert st["O1"].endswith("--ref-copies --local /R/local")
    st = dict((k, c) for k, c, _ in F.chain_steps(_opts(ref_copies=True, ref_carried=True, local="/L"), "/R", "INTACT"))
    assert st["O1"].endswith("--ref-copies --local /L --ref-carried")
    st = dict((k, c) for k, c, _ in F.chain_steps(_opts(ref_copies=True, aln_cache=True, local="/L"), "/R", "INTACT"))
    assert st["O1"].startswith("env O1_ALN_CACHE=/R/.o1_aln_cache python3 bench/o3_graph/o1seed.py /R/INTACT")
