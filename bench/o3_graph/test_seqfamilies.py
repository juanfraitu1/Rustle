import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_seqfamilies", os.path.join(HERE, "seqfamilies.py"))
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)


def test_a_gene_s_model_is_carried_by_another_gene_s_genome_when_the_alignment_s_introns_are_its_junctions():
    # gene A: exons of 100, 50 and 80 bp (junctions at 100 and 150 of its transcript), spliced-aligned to B's locus from 1,000
    cg = [(100, "M"), (300, "N"), (50, "M"), (400, "N"), (80, "M")]
    assert S.structure_match([100, 150], 230, "+", 1000, cg)
    assert S.structure_match([95, 155], 230, "+", 1000, cg)                                     # within 10 bp
    assert not S.structure_match([100, 120], 230, "+", 1000, cg)                                # an intron at no junction
    assert not S.structure_match([100, 150, 200], 230, "+", 1000, cg)                           # a junction B's genome lacks
    # a processed retrocopy (one exon) spliced onto its parent's locus: introns the copy does not have
    assert not S.structure_match([], 230, "+", 1000, cg)
    # B's genome holds A's two exons contiguously (OR4F21 onto OR4F16: 398M1D596M): A's intron is not there
    assert not S.structure_match([398], 994, "+", 1000, [(398, "M"), (1, "D"), (596, "M")])
    # minus strand: the query runs backwards along the target, its junctions counted from its other end
    assert S.structure_match([80, 130], 230, "-", 1000, cg)
    # two single-exon genes
    assert S.structure_match([], 230, "+", 1000, [(230, "M")])


def test_families_are_the_components_of_links_with_both_sequence_and_structure_and_avoid_spent_genes():
    links = [("a", "b", True), ("b", "c", True), ("c", "d", False), ("e", "f", True), ("g", "h", True)]
    fams = S.components(links, spent={"h"})
    assert fams == [["a", "b", "c"], ["e", "f"]]                    # c-d fails structure; g-h touches a spent gene


def test_the_draw_takes_families_in_seeded_order_up_to_the_copy_budget():
    fams = [dict(family=f"F{i}", copies=[{}] * k) for i, k in enumerate([2, 3, 2, 40, 2, 5, 2, 2])]
    got = S.draw(fams, seed=17, budget=10, max_total=20)
    assert sum(len(f["copies"]) for f in got) >= 10 and sum(len(f["copies"]) for f in got) <= 20
    assert all(len(f["copies"]) < 40 for f in got)                  # the 40-copy family would exceed max_total
    assert got == S.draw(fams, seed=17, budget=10, max_total=20)    # deterministic


def test_a_junction_bordering_an_unaligned_end_of_the_transcript_needs_no_intron():
    # LGALS7 (round 17): its main transcript ends in a 6-bp exon minimap2 cannot align -- the twin's alignment covers 460 of
    # 466 bp with 2 introns against 3 junctions. A junction at the edge of an unaligned end is no structural difference.
    cg = [(169, "M"), (689, "N"), (202, "M"), (984, "N"), (89, "M")]
    assert S.structure_match([169, 371, 460], 466, "+", 1000, cg)              # junction 460: the aligned part ends there
    assert S.structure_match([6, 175, 377], 466, "+", 1000, cg, q0=6)           # an unaligned 5' piece of 6 bp
    assert not S.structure_match([169, 300, 371, 460], 466, "+", 1000, cg)     # a junction inside the aligned part: missing


def _g(gid, chrom, start, end, n=1, biotype="protein_coding", name=None):
    return dict(gene_id=gid, gene_name=name or gid, chrom=chrom, start=start, end=end, strand="+", n=n, biotype=biotype)


def test_paralogs_overlapping_a_coding_gene_are_dropped_and_among_overlapping_paralogs_a_pseudogene_wins():
    coding = [_g("OR4F16", "chr5", 1000, 2000)]
    nc = [_g("OR4F1P", "chr5", 1019, 1958, biotype="unprocessed_pseudogene"),         # annotated inside OR4F16: dropped
          _g("lnc1", "chr6", 100, 5000, n=4, biotype="lncRNA"),                       # holds OR4F7P: the pseudogene is the copy
          _g("OR4F7P", "chr6", 3000, 3900, n=1, biotype="unprocessed_pseudogene"),
          _g("OR4F2P", "chr11", 100, 1000, biotype="unprocessed_pseudogene"),
          _g("lnc2", "chr12", 100, 900, n=1, biotype="lncRNA"),                       # two overlapping lncRNAs: more transcripts
          _g("lnc3", "chr12", 500, 1500, n=3, biotype="lncRNA")]
    assert {g["gene_id"] for g in S.add_paralog_genes(coding, nc)} == {"OR4F7P", "OR4F2P", "lnc3"}


def test_a_paralog_links_to_a_coding_gene_only_when_its_locus_carries_the_coding_gene_s_model():
    coding = {"A", "B"}
    cover = {("A", "B"): 0.9, ("B", "A"): 0.3, ("A", "P"): 0.95, ("P", "A"): 1.0, ("Q", "A"): 1.0, ("A", "Q"): 0.2,
             ("P", "Q"): 1.0, ("Q", "P"): 1.0}
    struct = {("A", "B"): False, ("B", "A"): True, ("A", "P"): True, ("P", "A"): False, ("Q", "A"): True, ("P", "Q"): True,
              ("Q", "P"): True}
    assert S.pair_link("A", "B", cover, struct, coding) is True        # coding pair: either direction, as before
    assert S.pair_link("A", "P", cover, struct, coding) is True        # P's locus carries A's model
    assert S.pair_link("P", "A", cover, struct, coding) is True        # the same pair, given the other way round
    assert S.pair_link("A", "Q", cover, struct, coding) is None        # Q's own model sits in A's locus: not a copy of A
    assert S.pair_link("P", "Q", cover, struct, coding) is None        # two paralogs: no link of their own
    struct[("A", "P")] = False
    assert S.pair_link("A", "P", cover, struct, coding) is False       # a sequence link without A's exon chain


def test_a_paralog_family_is_named_after_its_first_coding_gene():
    gs = [_g("P1", "chr1", 10, 20, biotype="unprocessed_pseudogene", name="OR4F2P"), _g("C1", "chr1", 50, 60, name="OR4F21"),
          _g("C2", "chr5", 5, 9, name="OR4F16")]
    used = {}
    assert S.family_name(gs, used, "S19") == "S19OR4F21"
    assert S.family_name(gs, used, "S19") == "S19OR4F21x2"


def test_the_draw_reads_a_pool_of_any_set_name(tmp_path):
    import json
    pool = dict(genome="g", sets={"seq19": [dict(family=f"F{i}", copies=[{}] * 2) for i in range(60)]})
    json.dump(pool, open(tmp_path / "all.json", "w"))
    S.cmd_draw(str(tmp_path / "set.json"), str(tmp_path / "all.json"), seed=19, set_name="r19")
    got = json.load(open(tmp_path / "set.json"))
    assert list(got["sets"]) == ["r19"] and sum(len(f["copies"]) for f in got["sets"]["r19"]) == 90


def _cigar_with_introns(at, qlen, intron=300):
    """A CIGAR aligning a qlen-bp transcript contiguously except for an intron after each query position in `at`."""
    out, q = [], 0
    for p in sorted(at):
        out += [(p - q, "M"), (intron, "N")]
        q = p
    return out + [(qlen - q, "M")]


def test_an_intron_chain_is_carried_when_each_junction_has_its_intron_even_if_the_other_genome_has_more():
    # GOLGA6L6 onto GOLGA6L22's locus (2026-10-18): all 8 junctions of GOLGA6L6 have their intron, and GOLGA6L22's genome has
    # two more inside GOLGA6L6's 1,431-bp exon (its own 765 + 327 + 297 split)
    js = [271, 400, 451, 545, 652, 711, 853, 2284]
    cg = _cigar_with_introns([271, 400, 451, 545, 652, 711, 853, 1639, 2095, 2284], 4031)
    assert S.structure_match(js, 4031, "+", 1000, cg, chain="carried")
    assert not S.structure_match(js, 4031, "+", 1000, cg)                     # the exact rule (default) asks for no more
    # a junction without its intron is still missing
    assert not S.structure_match(js + [3000], 4031, "+", 1000, cg, chain="carried")


def test_a_carried_intron_chain_still_keeps_retrocopies_out():
    cg = [(100, "M"), (300, "N"), (50, "M"), (400, "N"), (80, "M")]
    # a single-exon model (a retrocopy) spliced onto its parent's locus: introns it does not have
    assert not S.structure_match([], 230, "+", 1000, cg, chain="carried")
    # the parent's model onto the retrocopy's genome: its junctions find no intron
    assert not S.structure_match([100, 150], 230, "+", 1000, [(230, "M")], chain="carried")
    # two single-exon genes
    assert S.structure_match([], 230, "+", 1000, [(230, "M")], chain="carried")


def test_the_build_takes_its_intron_chain_rule_and_chain_retention_from_the_command_line():
    assert S.build_args(["out.json", "a.json", "b.json"]) == ("out.json", ["a.json", "b.json"], False, "exact", False)
    assert S.build_args(["out.json", "a.json", "--paralogs", "--intron-chain", "carried"]) == ("out.json", ["a.json"], True,
                                                                                                "carried", False)
    assert S.build_args(["out.json", "--retain-chains", "a.json"]) == ("out.json", ["a.json"], False, "exact", True)
    assert S.build_args(["out.json", "--transcribed-seeds"]) == ("out.json", [], False, "exact", False)
    assert S.build_seeds(["out.json", "--transcribed-seeds"]) and not S.build_seeds(["out.json"])


def test_transcribed_pseudogenes_seed_a_family_like_coding_genes_and_others_do_not():
    # XKRY (2026-10-18): both copies are transcribed pseudogenes in CAT and in RefSeq, so no coding gene seeds their family
    others = [_g("XKRY", "chrY", 100, 900, biotype="transcribed_unprocessed_pseudogene"),
              _g("XKRY2", "chrY", 5000, 5800, biotype="transcribed_unprocessed_pseudogene"),
              _g("P1", "chr1", 100, 900, biotype="unprocessed_pseudogene"),
              _g("L1", "chr2", 100, 900, biotype="lncRNA")]
    assert S.seed_genes([_g("C1", "chr3", 1, 2)], others) == {"C1", "XKRY", "XKRY2"}
    cover = {("XKRY", "XKRY2"): 1.0, ("XKRY2", "XKRY"): 1.0}
    struct = {("XKRY", "XKRY2"): True, ("XKRY2", "XKRY"): True}
    assert S.pair_link("XKRY", "XKRY2", cover, struct, {"XKRY", "XKRY2"}) is True       # seeds: they link
    assert S.pair_link("XKRY", "XKRY2", cover, struct, set()) is None                   # without the option: no link


def test_the_locus_alignment_keeps_every_chain_only_when_asked():
    # TSPY2 (2026-10-18): among ~45 near-identical TSPY loci, minimap2's default chain selection dropped the alignment of
    # 62 of 63 TSPY transcripts to the most divergent copy's locus; -P (retain all chains) kept 46
    assert "-P" not in S.splice_cmd(4)
    assert S.splice_cmd(4, retain=True)[-1] == "-P"
    assert S.splice_cmd(4)[:3] == ["minimap2", "-cx", "splice:hq"]


def test_the_build_takes_a_span_cap_or_none():
    assert S.build_span(["out.json"]) is None                        # families.MAX_SPAN, as before
    assert S.build_span(["out.json", "--max-span", "0"]) == 0        # no cap: the definition holds long genes (DAZ1, 61 kb)
    assert S.build_args(["out.json", "--max-span", "0", "a.json"]) == ("out.json", ["a.json"], False, "exact", False)


def test_genes_are_candidates_when_any_of_their_isoforms_align():
    # NPIPA1 (2026-10-18): its canonical is a 27-exon, 5.4-kb isoform at 0.87-0.90 identity to the other NPIP canonicals, so
    # the canonical prefilter never offered it to the locus rule -- though NPIPA2's canonical is carried by its locus
    paf = "\t".join(["NPIPA1|t3", "1084", "0", "1084", "+", "NPIPA2|t1", "757", "0", "757", "740", "757", "60"]) + "\n"
    paf += "\t".join(["NPIPA1|t3", "1084", "0", "1084", "+", "NPIPA1|t1", "5410", "0", "1084", "1084", "1084", "60"]) + "\n"
    paf += "\t".join(["NPIPA5|t1", "1084", "0", "500", "+", "NPIPA2|t1", "757", "0", "500", "400", "500", "60"]) + "\n"
    assert S.candidate_genes(paf) == {("NPIPA1", "NPIPA2")}                 # same-gene pairs out; 0.80 identity out
    assert S.build_args(["o.json", "--candidate-isoforms"]) == ("o.json", [], False, "exact", False)
    assert S.build_isoforms(["o.json", "--candidate-isoforms"]) and not S.build_isoforms(["o.json"])
