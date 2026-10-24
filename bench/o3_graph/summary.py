#!/usr/bin/env python3
"""Summaries the graph-O2/O3 verdicts are read from (round 4 review I4: these were scratch scripts before).

    summary.py o3 OUT.json WORLD ...      per scenario and arm: detections (REMOVED: of the DETECTABLE removed copies, median
                                          primary `de` on the world's genome > DELTA) and false-call families (score.py's rule)
    summary.py o2bars PREFIX WORLD ...    round 4's A1-A4: padded union O2 (o2.unionp.m2.tsv) against the current O2 and locus O2
                                          (o2.locus.m2.tsv), per world; writes PREFIX_<arm>.json (o2graph.py score)
    summary.py o2bars-unionf PREFIX WORLD ...   the same for round 5's graph O2 with the locus fallback (o2.unionf.m2.tsv)
    summary.py e2e OUT.json WORLD ... [--per-read-detectability] [--o2 O2_NAME] [--error-floor]   end to end: every read's fate (read_fate) by the kind of copy it came from, for the
                                          current method and the integrated O1 -> O2 -> O3
    summary.py bars E2E.json              round 8's pre-registered bars (docs/archive/2026-10/merged/PREREG_o3_graph_round8_2026-10-16.md) on an e2e
                                          table: per world, integrated correct >= current (B1) and errors <= current (B2)
    summary.py flags OUT.json WORLD ...   round 6: the flagged clusters (o2graph.flag_reads' rule) by truth class -- copy /
                                          allele (one haplotype of a catalogued copy) / other -- and the uncatalogued copies covered
"""
import collections
import csv
import importlib.util
import json
import os
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DELTA = 0.00958
FLAG_MIN_READS = 20
ARMS = ("CUR", "S23", "S3", "S123m", "G", "H", "Hn", "K", "Kn", "I", "In")
NONCAT_WORLDS = ("OUTSIDE", "REMOVED", "BRANCH")


def detectable(des, delta=DELTA):
    """O3's definition: nothing in the genome explains the reads within delta (their median primary `de` > delta)."""
    return bool(des) and statistics.median(des) > delta


def o2_bars(worlds):
    """worlds: {world name: counters (cur_correct, graph_correct, graph_assigned, graph_wrong_family, locus_correct,
    cur_assigned_noncat, graph_assigned_noncat)} -> A1 (lowest correct ratio), A2 (lowest precision), A3 (fewer in REMOVED and
    BRANCH; pooled ratio over OUTSIDE, REMOVED, BRANCH), A4 (padded union / locus correct, pooled)."""
    a1 = min(v["graph_correct"] / max(1, v["cur_correct"]) for v in worlds.values())
    a2 = min(v["graph_correct"] / max(1, v["graph_assigned"] + v.get("graph_wrong_family", 0)) for v in worlds.values())
    nc = {n: (v.get("graph_assigned_noncat", 0), v.get("cur_assigned_noncat", 0)) for n, v in worlds.items()
          if n.split("_")[0] in NONCAT_WORLDS}
    fewer = all(g < c for n, (g, c) in nc.items() if n.split("_")[0] in ("REMOVED", "BRANCH"))
    pooled = sum(g for g, _ in nc.values()) / max(1, sum(c for _, c in nc.values()))
    a4 = sum(v["graph_correct"] for v in worlds.values()) / max(1, sum(v["locus_correct"] for v in worlds.values()))
    return dict(A1=round(a1, 4), A2=round(a2, 4), A3_fewer=fewer, A3_pooled=pooled, A4=a4)


def classify_cluster(members, truth, nonc, purity=0.9):
    """A flagged cluster's truth (round 6): `copy` when >= purity of it is one uncatalogued copy's reads; `allele` when >= purity
    is one catalogued copy's reads AND >= purity of those come from one haplotype (a heterozygous split); `other` otherwise
    (another family's reads, an isoform artifact, a mixture). -> (class, majority copy, its share)."""
    comp = collections.Counter(truth[q]["copy"] for q in members if q in truth)
    if not comp:
        return "other", None, 0.0
    top, k = comp.most_common(1)[0]
    share = k / sum(comp.values())
    if share >= purity and top in nonc:
        return "copy", top, share
    if share >= purity:
        haps = collections.Counter(truth[q].get("hap", "1") for q in members if q in truth and truth[q]["copy"] == top)
        if haps and haps.most_common(1)[0][1] >= purity * k and len(haps) <= 2 and set(haps) != {"1"}:
            return "allele", top, share
    return "other", top, share


def cmd_flags(out, worlds):
    """Round 6: per world, the flagged clusters by truth class and the uncatalogued copies they cover (>= 80% of a copy's hosted
    reads in its `copy` clusters)."""
    res = {}
    for w in worlds:
        truth = {r["read"]: r for r in csv.DictReader(open(f"{w}/reads.truth.E0.tsv"), delimiter="\t")}
        kinds = json.load(open(f"{w}/kinds.json"))
        nonc = {c for c, k in kinds.items() if k in ("absent", "outside")}
        mem = collections.defaultdict(set)
        for r in csv.DictReader(open(f"{w}/cur/A.psv_cluster_reads.tsv"), delimiter="\t"):
            mem[(r["family"], r["host"], r["cluster"])].add(r["read_name"])
        hosted = collections.defaultdict(set)
        for r in csv.DictReader(open(f"{w}/cur/A.assignments.tsv"), delimiter="\t"):
            hosted[r["family_id"]].add(r["read_name"])
        cls = collections.Counter()
        false_fams = collections.defaultdict(set)
        covered = collections.defaultdict(set)
        for c in csv.DictReader(open(f"{w}/cur/A.psv_clusters.tsv"), delimiter="\t"):
            if c["status"] != "unexplained" or int(c["n_reads"]) < FLAG_MIN_READS:
                continue
            m = mem[(c["family"], c["host"], c["cluster"])]
            k, top, _ = classify_cluster(m, truth, nonc)
            cls[k] += 1
            if k == "copy":
                covered[top] |= m
            else:
                false_fams[k].add(c["family"])
        flagged_copies = []
        for cp in sorted(nonc):
            src = {q for q, t in truth.items() if t["copy"] == cp} & hosted[cp.rsplit(".", 1)[0]]
            if src and len(covered[cp] & src) >= 0.8 * len(src):
                flagged_copies.append(cp)
        name = os.path.basename(w.rstrip("/"))
        res[name] = dict(clusters=dict(cls), uncatalogued=len(nonc), flagged_copies=flagged_copies,
                         false_flag_families={k: sorted(v) for k, v in false_fams.items()},
                         n_false_flag_families=len(set().union(*false_fams.values())) if false_fams else 0)
        print(f"{name}: flagged clusters {dict(cls)} | uncatalogued copies flagged {len(flagged_copies)}/{len(nonc)} | "
              f"families with a false flag {res[name]['n_false_flag_families']} {res[name]['false_flag_families']}")
    json.dump(res, open(out, "w"), indent=1)
    return 0


def absent_class(copy_detectable, read_de, per_read=False, floor=0.0):
    """A REMOVED copy's read: detectable or not. Copy level (rounds 7-8): its copy's reads' median `de` on the masked genome
    > DELTA. Per read (round 9, `--per-read-detectability`): the read's OWN primary `de` > DELTA (a read without a primary
    alignment is unexplained by the genome, hence detectable) -- a read of a detectable copy that is near-identical to a
    present copy over its own span cannot be told from that copy's reads. floor (round 11, `--error-floor`): the world's
    sequencing-error floor (detect.error_floor) is subtracted from the read's `de` first -- DELTA is the p99 exonic divergence
    between KB3781's haplotype assemblies (register 1199), an error-free divergence."""
    if per_read:
        return "absent_detectable" if read_de is None or read_de - floor > DELTA else "absent_undetectable"
    return "absent_detectable" if copy_detectable else "absent_undetectable"


def read_fate(cls, assigned, own, in_o3):
    """End to end (O1 -> O2 -> O3): a read's fate under one system, by the kind of copy it came from. An absent-copy call
    (O3) takes precedence over an assignment (O2). catalogued / uncatalogued: right only on its own copy (an uncatalogued copy
    only once O1 seeded it); an O3 call on it is an error. absent_detectable: right only inside an O3 call. absent_undetectable
    (the genome explains it within delta): placing it is `unknowable`, an O3 call is right, abstaining is neutral."""
    if in_o3:
        return "error" if cls in ("catalogued", "uncatalogued") else "correct"
    if assigned is None:
        return "abstain"
    if cls in ("catalogued", "uncatalogued"):
        return "correct" if assigned == own else "error"
    return "error" if cls == "absent_detectable" else "unknowable"


def integrated_calls(rows, seeded_copy):
    """{read: {family: the truth copy of the seeded copy it was called to}} from o2.integ rows. A call on a seeded copy that
    matches no truth copy (map_new_copies: None) places the read on a locus that is no copy: a placement, not an abstention
    (round 8 review I3; the frozen round-8 scorer counted it as abstain)."""
    out = {}
    for r in rows:
        if r["call"] != "-":
            k = (r["family"], int(r["call"]))
            out.setdefault(r["read"], {})[r["family"]] = seeded_copy.get(k) or f"{k[0]}#{k[1]}: no truth copy"
    return out


def seeded_copies(w):
    """(cat_copy, seeded_copy) of a world: {(family, catalogue idx): truth copy} (isoforms.json), and the same for every seeded
    copy (graph/seeded.json), a new one mapped to the truth copy it overlaps most (o2graph.map_new_copies)."""
    O = _o2("o2graph")
    iso = json.load(open(f"{w}/isoforms.json"))
    cat_copy = {(f, c["idx"]): c["copy"] for f, cps in iso.items() for c in cps}
    sj = json.load(open(f"{w}/graph/seeded.json"))
    truth_copies = list(csv.DictReader(open(f"{w}/copies.tsv"), delimiter="\t"))
    seeded_copy = dict(cat_copy)
    for f, v in sj.items():
        for i, cp in O.map_new_copies([c for c in v["copies"] if c["origin"] == "new"], truth_copies).items():
            seeded_copy[(f, i)] = cp
    return cat_copy, seeded_copy


def world_calls(w, o2_name="o2.integ.m2.tsv"):
    """One world's calls under both systems: {current, integrated: {read: {family: truth copy}}, current_o3, integrated_o3: the
    reads inside an absent-copy call}. current: copy_assign's assigned reads and its unexplained PSV clusters beyond DELTA;
    integrated: the O2 table o2_name on graph/seeded.json (a new seeded copy maps to the truth copy it overlaps most) and O3's
    arm-I calls."""
    cat_copy, seeded_copy = seeded_copies(w)
    cur = {}
    for r in csv.DictReader(open(f"{w}/cur/A.assignments.tsv"), delimiter="\t"):
        if r["status"] == "assigned" and r["catalog_copy_idx"].isdigit():
            cur.setdefault(r["read_name"], {})[r["family_id"]] = cat_copy.get((r["family_id"], int(r["catalog_copy_idx"])))
    integ = integrated_calls(csv.DictReader(open(f"{w}/graph/{o2_name}"), delimiter="\t"), seeded_copy)
    mem = collections.defaultdict(set)
    for r in csv.DictReader(open(f"{w}/cur/A.psv_cluster_reads.tsv"), delimiter="\t"):
        mem[(r["family"], r["host"], r["cluster"])].add(r["read_name"])
    cur_o3 = set()
    for c in csv.DictReader(open(f"{w}/cur/A.psv_clusters.tsv"), delimiter="\t"):
        if c["status"] == "unexplained" and float(c["est_div"]) > DELTA:
            cur_o3 |= mem[(c["family"], c["host"], c["cluster"])]
    i_o3 = set()
    for c in csv.DictReader(open(f"{w}/graph/calls.tsv"), delimiter="\t"):
        if c["arm"] == "I":
            i_o3 |= set(c["members"].split(","))
    return dict(current=cur, integrated=integ, current_o3=cur_o3, integrated_o3=i_o3)


def off_family_ties(w, slack, seq_tie=False):
    """The reads whose alignments in WORLD/E0.bam within `slack` of their best AS include one on no copy of the world
    (copies.tsv, in the reference): a locus the simulation never expresses, so neither placing nor abstaining on the read can
    be judged against its truth (user 2026-10-19: such ties are scored apart; slack 9 = Eichler's within-10 rule). seq_tie: an
    alignment there with no more edits, counting its unaligned bases, than the best-AS one ties too (o2graph.best_as_loci)."""
    import pysam
    vd = _o2("o2graph").veto_distance
    loci = collections.defaultdict(list)
    for c in csv.DictReader(open(f"{w}/copies.tsv"), delimiter="\t"):
        if c["in_reference"] == "1":
            loci[c["chrom"]].append((int(c["start"]), int(c["end"])))
    al = collections.defaultdict(list)
    with pysam.AlignmentFile(f"{w}/E0.bam") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_supplementary or not rec.has_tag("AS"):
                continue
            on = any(rec.reference_start < e and rec.reference_end > s for s, e in loci.get(rec.reference_name, ()))
            al[rec.query_name].append((rec.get_tag("AS"), on, vd(rec.cigartuples)))
    out = set()
    for q, v in al.items():
        best = max(a for a, _, _ in v)
        ref = min(d for a, _, d in v if a == best)
        if any(not on for a, on, d in v if a >= best - slack or seq_tie and d <= ref):
            out.add(q)
    return out


SET_KEYS = ("ambiguous", "placed right", "placed wrong", "abstain", "v17 placed right", "v17 placed wrong", "v17 abstain",
            "set holds truth", "candidates", "set size", "minimap2 right", "minimap2 wrong")


def set_fates(rows, truth, seeded_copy, dropped, primary, skip=frozenset()):
    """The single-genome O2 (o2graph --fewest-edits) on the reads tied in their own family (rows with candidates):
    placements right / wrong / abstained and v17's (base_call) beside them; whether the fewest-edit set holds the read's
    copy; the summed sizes of the tied copies and of the set; minimap2's primary right / wrong. truth: {read: copy};
    seeded_copy: {(family, seeded idx): truth copy}; dropped: {family: the truth copies O1 dropped as silent}; primary:
    {read: the copy its primary alignment lies on}; skip: reads left out (the off-family ties)."""
    n = dict.fromkeys(SET_KEYS, 0)
    for r in rows:
        q, f = r["read"], r["family"]
        t = truth.get(q)
        if t is None or q in skip or r["candidates"] == "-" or f != t.rsplit(".", 1)[0]:
            continue
        held = set()
        for x in r["set"].split(","):
            held |= set(dropped.get(f, ())) if x == "dropped" else {seeded_copy.get((f, int(x)))}
        n["ambiguous"] += 1
        for key, c in (("", r["call"]), ("v17 ", r["base_call"])):
            n[key + ("abstain" if c == "-" else "placed right" if seeded_copy.get((f, int(c))) == t else "placed wrong")] += 1
        n["set holds truth"] += t in held
        n["candidates"] += len(r["candidates"].split(","))
        n["set size"] += len(r["set"].split(","))
        n["minimap2 right" if primary.get(q) == t else "minimap2 wrong"] += 1
    return n


def attribute(blocks, strand, cands):
    """The copy an alignment lies on: among cands [(copy, strand, exons, (span start, end))] whose span its aligned blocks
    touch, those on its strand (a sense read's), the one whose exons the blocks overlap most, then whose span, then the first;
    None off every copy. Amendment 1 of PREREG_o2_single_genome_2026-10-19 (review of d63df788): a nested copy (CYB5B's
    single-exon AC026464.2) or an antisense neighbour (RPL38's MGC16275) shares the span, so the first span hit was no
    attribution."""
    def ov(a, b):
        return sum(max(0, min(e, t) - max(s, u)) for s, e in a for u, t in b)
    hits = [h for h in cands if any(s < h[3][1] and e > h[3][0] for s, e in blocks)]
    if not hits:
        return None
    same = [h for h in hits if h[1] == strand] or hits
    return max(enumerate(same), key=lambda kv: (ov(blocks, kv[1][2]), ov(blocks, [kv[1][3]]), -kv[0]))[1][0]


def tied_records(alns, slack, seq_tie=False):
    """How many of a read's alignments [(AS, veto_distance)] are tied: within `slack` of the best AS, or (seq_tie) with no
    more edits than the best-AS one (o2graph.best_as_loci's filter)."""
    best = max(a for a, _ in alns)
    ref = min(d for a, d in alns if a == best)
    return sum(1 for a, d in alns if a >= best - slack or seq_tie and d <= ref)


def _copy_cands(w):
    """{chrom: [(copy, strand, exons, span)]} of the world's copies in the reference; a catalogued copy's exons are its
    isoforms' (isoforms.json), another's its span."""
    exons = collections.defaultdict(set)
    for cps in json.load(open(f"{w}/isoforms.json")).values():
        for c in cps:
            for ex in c["isoforms"].values():
                exons[c["copy"]] |= {tuple(x) for x in ex}
    out = collections.defaultdict(list)
    for c in csv.DictReader(open(f"{w}/copies.tsv"), delimiter="\t"):
        if c["in_reference"] == "1":
            sp = (int(c["start"]), int(c["end"]))
            out[c["chrom"]].append((c["copy"], c["strand"], sorted(exons.get(c["copy"], ())) or [sp], sp))
    return out


def primary_copies(w):
    """{read: the copy its primary alignment in WORLD/E0.bam lies on (attribute), or None}."""
    import pysam
    cands = _copy_cands(w)
    out = {}
    with pysam.AlignmentFile(f"{w}/E0.bam") as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue
            out[rec.query_name] = attribute(rec.get_blocks(), "-" if rec.is_reverse else "+", cands.get(rec.reference_name, ()))
    return out


def single_tied_alignment(w, slack, seq_tie=False):
    """The reads with ONE tied alignment (tied_records) in WORLD/E0.bam: minimap2 placed them once, so they are no
    multimappers even when that alignment spans two overlapping copies (Amendment 1: CYB5B, RPL38)."""
    import pysam
    vd = _o2("o2graph").veto_distance
    al = collections.defaultdict(list)
    with pysam.AlignmentFile(f"{w}/E0.bam") as bam:
        for rec in bam.fetch(until_eof=True):
            if not (rec.is_unmapped or rec.is_supplementary) and rec.has_tag("AS"):
                al[rec.query_name].append((rec.get_tag("AS"), vd(rec.cigartuples)))
    return {q for q, v in al.items() if tied_records(v, slack, seq_tie) == 1}


def cmd_sets(out, worlds, o2_name="o2.intege.m2.tsv", off_slack=None, seq_tie=False):
    """Per world (CELL/WORLD) and pooled: set_fates of the single-genome O2 table on every catalogued read, the off-family
    ties (off_family_ties) and the reads with one tied alignment (single_tied_alignment, Amendment 1) apart."""
    res, pooled = {}, dict.fromkeys(SET_KEYS, 0)
    for w in worlds:
        name = f"{os.path.basename(os.path.dirname(w.rstrip('/')))}/{os.path.basename(w.rstrip('/'))}"
        kinds = json.load(open(f"{w}/kinds.json"))
        truth = {r["read"]: r["copy"] for r in csv.DictReader(open(f"{w}/reads.truth.E0.tsv"), delimiter="\t")
                 if kinds.get(r["copy"], "present") in ("present", "silent")}
        cat_copy, seeded = seeded_copies(w)
        sj = json.load(open(f"{w}/graph/seeded.json"))
        dropped = {f: {cat_copy[(f, i)] for i in v.get("dropped", []) if (f, i) in cat_copy} for f, v in sj.items()}
        skip = (off_family_ties(w, off_slack, seq_tie) | single_tied_alignment(w, off_slack, seq_tie)
                if off_slack is not None else set())
        n = set_fates(csv.DictReader(open(f"{w}/graph/{o2_name}"), delimiter="\t"), truth, seeded, dropped,
                      primary_copies(w), skip)
        res[name] = n
        for k in SET_KEYS:
            pooled[k] += n[k]
        print(f"{name}: " + ", ".join(f"{k} {v}" for k, v in n.items()), flush=True)
    res["pooled"] = pooled
    print("POOLED: " + ", ".join(f"{k} {v}" for k, v in pooled.items()), flush=True)
    json.dump(res, open(out, "w"), indent=1)
    return 0


def cmd_e2e(out, worlds, per_read=False, o2_name="o2.integ.m2.tsv", error_floor=False, off_slack=None, seq_tie=False):
    """Per world and read class: the fates of every read under the current method (copy_assign O2 + CUR O3) and the integrated
    O1 -> O2 -> O3 (o2_name, by default o2.integ.m2.tsv, on graph/seeded.json + arm I)."""
    res = {}
    for w in worlds:
        name = os.path.basename(w.rstrip("/"))
        truth = {r["read"]: r for r in csv.DictReader(open(f"{w}/reads.truth.E0.tsv"), delimiter="\t")}
        kinds = json.load(open(f"{w}/kinds.json"))
        det = removed_detectability(w) if name.startswith("REMOVED") else {}
        rde = primary_de(w) if per_read and name.startswith("REMOVED") else {}
        floor = _o2("detect").error_floor(f"{w}/E0.bam") if error_floor and name.startswith("REMOVED") else 0.0

        def cls_of(t):
            k = kinds.get(t["copy"], "present")
            if k in ("present", "silent"):
                return "catalogued"
            if k == "outside":
                return "uncatalogued"
            if name.startswith("BRANCH"):                  # a novel-exon absent copy: detectable by construction
                return "absent_detectable"
            return absent_class(det.get(t["family"], (0, False))[1], rde.get(t["read"]), per_read, floor)
        wc = world_calls(w, o2_name)
        oft = off_family_ties(w, off_slack, seq_tie) if off_slack is not None else set()
        tab = {}
        for system in ("current", "integrated"):
            o2, o3 = wc[system], wc[f"{system}_o3"]
            cnt = collections.defaultdict(collections.Counter)
            for q, t in truth.items():
                calls = o2.get(q, {})
                assigned = calls.get(t["family"], next(iter(calls.values())) if calls else None)
                if t["family"] not in calls and calls:
                    assigned = "wrong-family"
                cl = cls_of(t)
                cnt["off_family_tie" if q in oft else cl][read_fate(cl, assigned, t["copy"], q in o3)] += 1
            tab[system] = {k: dict(v) for k, v in cnt.items()}
        res[name] = tab
        print(name)
        for system in ("current", "integrated"):
            print(f"   {system:10s} " + " | ".join(f"{cl}: " + " ".join(f"{k} {v}" for k, v in sorted(tab[system][cl].items()))
                                                for cl in sorted(tab[system])))
    json.dump(res, open(out, "w"), indent=1)
    return 0


def e2e_args(argv):
    """e2e OUT.json WORLD ... [--per-read-detectability] [--o2 O2_NAME] [--error-floor] [--off-family-ties SLACK] -> (out,
    worlds, per_read, o2_name, error_floor, slack, seq_tie); O2_NAME: the integrated system's O2 table in WORLD/graph (o2.integ.m2.tsv;
    o2.integl.m2.tsv: with the locus fallback); --error-floor: per-read detectability less each REMOVED world's error floor
    (absent_class); --off-family-ties: the reads off_family_ties finds are the class `off_family_tie`, which the bars leave out
    (--seq-tie: with the O2 scope's sequence ties)."""
    rest, o2_name, slack, i = [], "o2.integ.m2.tsv", None, 2
    while i < len(argv):
        if argv[i] == "--o2":
            o2_name, i = argv[i + 1], i + 2
            continue
        if argv[i] == "--off-family-ties":
            slack, i = int(argv[i + 1]), i + 2
            continue
        if not argv[i].startswith("--"):
            rest.append(argv[i])
        i += 1
    return argv[1], rest, "--per-read-detectability" in argv, o2_name, "--error-floor" in argv, slack, "--seq-tie" in argv


UNJUDGED = ("off_family_tie",)                  # read classes the bars leave out (summary.off_family_ties)


def round8_bars(tab):
    """Round 8 (pre-registered): per world of an e2e table, the reads judged correct and wrong under each system, summed over
    every read class but UNJUDGED (abstain and unknowable count for neither); B1 integrated correct >= current, B2 integrated errors <=
    current. verdict PASS only when both hold in every world."""
    out = {}
    for w, t in sorted(tab.items()):
        n = {(s, x): sum(c.get(x, 0) for k, c in t[s].items() if k not in UNJUDGED) for s in ("current", "integrated")
             for x in ("correct", "error")}
        out[w] = dict(current_correct=n[("current", "correct")], integrated_correct=n[("integrated", "correct")],
                      current_error=n[("current", "error")], integrated_error=n[("integrated", "error")],
                      B1=n[("integrated", "correct")] >= n[("current", "correct")],
                      B2=n[("integrated", "error")] <= n[("current", "error")])
    out["verdict"] = "PASS" if all(v["B1"] and v["B2"] for v in out.values()) else "FAIL"
    return out


def cmd_bars(path):
    b = round8_bars(json.load(open(path)))
    for w, v in sorted(b.items()):
        if w != "verdict":
            print(f"{w:8s} correct current {v['current_correct']:6d} integrated {v['integrated_correct']:6d} B1 {v['B1']} | "
                  f"errors current {v['current_error']:5d} integrated {v['integrated_error']:5d} B2 {v['B2']}")
    print("verdict", b["verdict"])
    return 0


def drop_outcome(seeded, cat_copy, kinds):
    """Round 7: O1's drops against the truth -- silent copies dropped (wanted) and expressed copies dropped (the absorption risk)."""
    out = collections.Counter(silent=0, silent_dropped=0, expressed=0, expressed_dropped=0)
    for f, v in seeded.items():
        dropped = set(v["dropped"])
        for (ff, i), cp in cat_copy.items():
            if ff != f:
                continue
            k = "silent" if kinds.get(cp) == "silent" else "expressed"
            out[k] += 1
            out[f"{k}_dropped"] += i in dropped
    return dict(out)


def _o2(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def cmd_integ(out, worlds):
    """Round 7, per world: correct / wrong calls of the current O2, padded union O2 on the curated catalog, and the integrated
    O1 -> O2 (o2.integ.m2.tsv on graph/seeded.json); reads of uncatalogued in-genome copies called to their own seeded copy; O1's
    drops (silent vs expressed copies) and new copies."""
    O = _o2("o2graph")
    res = {}
    for w in worlds:
        truth = {r["read"]: r for r in csv.DictReader(open(f"{w}/reads.truth.E0.tsv"), delimiter="\t")}
        iso = json.load(open(f"{w}/isoforms.json"))
        cat_copy = {(f, c["idx"]): c["copy"] for f, cps in iso.items() for c in cps}
        cur = {(r["read_name"], r["family_id"]): int(r["catalog_copy_idx"])
               for r in csv.DictReader(open(f"{w}/cur/A.assignments.tsv"), delimiter="\t")
               if r["status"] == "assigned" and r["catalog_copy_idx"].isdigit()}

        def calls_of(name):
            return {(r["read"], r["family"]): int(r["call"]) for r in csv.DictReader(open(f"{w}/graph/{name}"), delimiter="\t")
                    if r["call"] != "-"}
        sj = json.load(open(f"{w}/graph/seeded.json"))
        truth_copies = list(csv.DictReader(open(f"{w}/copies.tsv"), delimiter="\t"))
        seeded_copy = dict(cat_copy)
        for f, v in sj.items():
            for i, cp in O.map_new_copies([c for c in v["copies"] if c["origin"] == "new"], truth_copies).items():
                seeded_copy[(f, i)] = cp
        kinds = json.load(open(f"{w}/kinds.json"))
        unk = set()
        if os.path.basename(w.rstrip("/")) == "REMOVED":         # absent copies the genome explains within delta
            fam_abs = {c: k for c, k in kinds.items() if k == "absent"}
            det = removed_detectability(w)
            unk = {c for c in fam_abs if not det.get(c.rsplit(".", 1)[0], (0, False))[1]}
        r = dict(current=O.score_assign(cur, cat_copy, truth, unk),
                 padded_union=O.score_assign(calls_of("o2.unionp.m2.tsv"), cat_copy, truth, unk),
                 integrated=O.score_assign(calls_of("o2.integ.m2.tsv"), seeded_copy, truth, unk), unknowable=sorted(unk),
                 drops=drop_outcome(sj, cat_copy, kinds),
                 new_copies=[(f, c["idx"], seeded_copy.get((f, c["idx"]))) for f, v in sj.items() for c in v["copies"] if c["origin"] == "new"],
                 outside_reads=sum(1 for t in truth.values() if t["kind"] == "outside"))
        name = os.path.basename(w.rstrip("/"))
        res[name] = r
        print(f"{name}: current {r['current']} | padded union {r['padded_union']} | integrated {r['integrated']} | drops {r['drops']} "
              f"| new copies {len(r['new_copies'])} ({sum(1 for x in r['new_copies'] if x[2])} on a truth copy) | outside reads {r['outside_reads']}")
    json.dump(res, open(out, "w"), indent=1)
    return 0


def primary_de(w):
    """{read: its primary alignment's `de`} on the world's genome."""
    import pysam
    out = {}
    with pysam.AlignmentFile(f"{w}/E0.bam") as bam:
        for rec in bam.fetch(until_eof=True):
            if not (rec.is_unmapped or rec.is_secondary or rec.is_supplementary):
                out[rec.query_name] = rec.get_tag("de")
    return out


def removed_detectability(w):
    """{family: (median de, detectable)} for a REMOVED world's absent copies (primary alignments only)."""
    import pysam
    truth = {r["read"]: r for r in csv.DictReader(open(f"{w}/reads.truth.E0.tsv"), delimiter="\t")}
    de = collections.defaultdict(list)
    with pysam.AlignmentFile(f"{w}/E0.bam") as bam:
        for rec in bam.fetch(until_eof=True):
            t = truth.get(rec.query_name)
            if t and t["kind"] == "absent" and not (rec.is_unmapped or rec.is_secondary or rec.is_supplementary):
                de[t["family"]].append(rec.get_tag("de"))
    return {f: (statistics.median(v), detectable(v)) for f, v in de.items()}


def cmd_o3(out, worlds):
    subprocess.run([sys.executable, os.path.join(HERE, "score.py"), out] + worlds, check=True, capture_output=True)
    rows = json.load(open(out))["rows"]
    det = {}
    for w in worlds:
        if os.path.basename(w.rstrip("/")) == "REMOVED":
            for f, (_, ok) in removed_detectability(w).items():
                det[(w, f)] = ok
    by = collections.defaultdict(collections.Counter)
    for r in rows:
        if r["arm"] not in ARMS:
            continue
        k = by[(r["scenario"], r["arm"])]
        k["families"] += 1
        k["false"] += r["false_call"]
        if r["scenario"] == "REMOVED":
            d = det.get((r["world"], r["family"]), False)
            k["detectable"] += d
            k["detected"] += bool(d and r["detected"])
        elif r["detected"] is not None:
            k["sources"] += 1
            k["detected"] += bool(r["detected"])
    for sc in sorted({s for s, _ in by}):
        cells = []
        for a in ARMS:
            k = by.get((sc, a))
            if not k:
                continue
            of = k["detectable"] if sc == "REMOVED" else k["sources"]
            cells.append(f"{a}: " + (f"det {k['detected']}/{of} " if of else "") + f"false {k['false']}/{k['families']}")
        print(sc, " | ".join(cells))
    return 0


def cmd_o2bars(prefix, worlds, graph_arm="unionp.m2"):
    res = {}
    for arm in (graph_arm, "locus.m2"):
        p = f"{prefix}_{arm}.json"
        subprocess.run([sys.executable, os.path.join(HERE, "o2graph.py"), "score", p, f"o2.{arm}.tsv"] + worlds, check=True,
                       capture_output=True)
        res[arm] = json.load(open(p))["worlds"]
    per = {}
    for w in worlds:
        u, lo = collections.Counter(), collections.Counter()
        for v in res[graph_arm][w].values():
            u.update(v)
        for v in res["locus.m2"][w].values():
            lo.update(v)
        u["locus_correct"] = lo["graph_correct"]
        per[os.path.basename(w.rstrip("/"))] = u
        print(f"{os.path.basename(w)}: current {u['cur_correct']} padded union {u['graph_correct']}/{u['graph_assigned']} "
              f"locus {lo['graph_correct']} | non-catalogued placed: current {u.get('cur_assigned_noncat', 0)} "
              f"padded union {u.get('graph_assigned_noncat', 0)}")
    print(o2_bars(per))
    return 0


def main(argv):
    if argv[0] == "o3":
        return cmd_o3(argv[1], argv[2:])
    if argv[0] == "bars":                                # round 8: the pre-registered end-to-end bars
        return cmd_bars(argv[1])
    if argv[0] == "sets":                                # the single-genome O2 (o2.intege): placements and fewest-edit sets
        out, worlds, _, o2_name, _, slack, seq_tie = e2e_args(argv)
        return cmd_sets(out, worlds, o2_name if "--o2" in argv else "o2.intege.m2.tsv", slack, seq_tie)
    if argv[0] == "e2e":                                 # end to end: every read's fate, current vs integrated
        return cmd_e2e(*e2e_args(argv))
    if argv[0] == "integ":                               # round 7: the integrated O1 -> O2 against the current O2
        return cmd_integ(argv[1], argv[2:])
    if argv[0] == "flags":                               # round 6: the flag's clusters by truth class
        return cmd_flags(argv[1], argv[2:])
    if argv[0] == "o2bars-unionf":                       # round 5: the graph O2 with the locus fallback
        return cmd_o2bars(argv[1], argv[2:], graph_arm="unionf.m2")
    return cmd_o2bars(argv[1], argv[2:])


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
