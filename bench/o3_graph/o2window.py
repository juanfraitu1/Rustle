#!/usr/bin/env python3
"""Widening O2's window past Eichler's 10 AS with an intron-free score (user 2026-10-20: "we could also go past the 10 AS
difference from Eichler, if that is too close to the fully tied impossible cases"; "Sounds good, please proceed").

A development sweep on gorilla KB3781 (PREREG_o2_personal_copies' development families only; descriptive, no verdict). A read
is tied in window W when its alignments within W of its best `ms` (minimap2's score without the intron charges; on these reads
it tracks the edit count) lie on two or more copies of one family, none off them, minimap2's primary among them. Per window:
    - minimap2's pick (its primary's copy);
    - O2's pick: the tied copy with the unique fewest edits over the read bases every tied alignment aligns (o2personal);
    - O2's set: the tied copies no other tied copy beats at every substitution column where the two differ (pareto_set) --
      the pool reduced without dropping a copy the read does not contradict.
An alignment that lays bases on two copies' exons (stitched(): SM392's chimeras) is no placement for O2's pick or set.

    o2window.py prepare SUB.bam COPIES.tsv QUAL.json GENOME.fa HALF OUT_PREFIX [--windows 10,20,30,50]
        OUT.calls.tsv (per read and window: copies, minimap2, O2 pick, O2 set, stitched), OUT.domain.tsv (o2correct-donor
        columns, for o2personal.py truth), OUT.names (the reads to align to the animal's assemblies)
    o2window.py truth-ms DOMAIN.tsv COPIES.tsv ORTH.tsv CHRMAP.tsv OUT.tsv --bams mat=BAM:FAI pat=BAM:FAI
        the answer key by the intron-free score instead of edits (a clipped read end is cheap, a mismatch 5x dearer)
    o2window.py truth-sim OUT.tsv WORLD...               (read, status, truth) from simulated worlds' reads.truth.E0.tsv
    o2window.py truth-hg002 OUT.tsv TRUTH.tsv...         the same from hg002truth tables (class -> status)
    o2window.py bars OUT.json --calls CALLS... --truth TRUTH... [--window 30]     PREREG_o2_window_set's bars (bars())
    o2window.py band OUT.json --calls CALLS... --truth TRUTH...       PREREG_o2_set_null_band: set vs keep-all vs primary, by band
    o2window.py score CALLS.tsv TRUTH.tsv OUT.json      per window and band: the truth's statuses; minimap2 and O2's pick
                                                        against the truth; O2's set: holds the truth, size, singletons
"""
import collections
import csv
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
STITCH_BP = 50                                                  # bases on a second copy's exons that make a stitch


def _sibling(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}_ow", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


PC = _sibling("o2personal")
C = PC.C


def _ov(blocks, exons):
    return sum(max(0, min(e, t) - max(s, u)) for s, e in blocks for u, t in exons)


def stitched(blocks, exons):
    """Whether an alignment's aligned blocks lay >= STITCH_BP bases on the exons of two or more copies ({copy: exons} of
    the copies its span touches)."""
    return sum(1 for ex in exons.values() if _ov(blocks, ex) >= STITCH_BP) >= 2


def window_copies(alns, w):
    """({copies}, minimap2's primary copy) of a read tied in window w: its alignments [{ms, copy, prim}] whose ms is within
    w of its best (difference < w, as Eichler's 10 AS is SLACK 9) lie on >= 2 copies of one family, none off every copy,
    minimap2's primary among them; the set is None otherwise."""
    prim = next((a["copy"] for a in alns if a["prim"]), None)
    best = max(a["ms"] for a in alns)
    near = [a for a in alns if a["ms"] > best - w]
    cs = {a["copy"] for a in near}
    if None in cs or len(cs) < 2 or len({c.rsplit(".", 1)[0] for c in cs}) > 1 or not any(a["prim"] for a in near):
        return None, prim
    return cs, prim


def pareto_set(cols, window):
    """The copies no rival beats at every substitution column where the two differ inside the window: rival r beats c when
    the read matches r and not c (r '=', c 'X') somewhere and never the reverse. cols: {copy: (=, X) read positions}."""
    lo, hi = window

    def wins(r, c):                                             # columns where the read sides with r against c
        return sum(1 for i in cols[c][1] if lo <= i < hi and i in cols[r][0])

    return {c for c in cols if not any(r != c and wins(r, c) > 0 and wins(c, r) == 0 for r in cols)}


MARGIN_MS = 5                                                   # one substitution: a lost match (1) + the mismatch (4)


def truth_ms(alns, margin=MARGIN_MS):
    """A read's origin in the animal's genome by minimap2's intron-free score (ms: a mismatch costs 5, a clipped base 1)
    from its alignments [(ms, reference copy or None, 'ref'|'alt')]: the one copy every alignment within `margin` of the best
    lies on -> ('determined', copy, 'refcons' when a reference-haplotype alignment is among them, else 'altonly'); two
    copies, or a copy and an off-copy locus -> 'undetermined'; the best off every copy -> 'off-copy'."""
    if not alns:
        return "unmapped", None, None
    best = max(m for m, _, _ in alns)
    w = [(m, c, r) for m, c, r in alns if m > best - margin]
    cs = {c for _, c, _ in w}
    if len(cs) > 1:
        return "undetermined", None, None
    c = cs.pop()
    if c is None:
        return "off-copy", None, None
    return "determined", c, "refcons" if any(r == "ref" for _, _, r in w) else "altonly"


def match_columns_ref(cigar, flip, length, seq, start, fetch):
    """(read positions matching the reference, read positions mismatching it) in the primary's orientation, from an
    alignment with or without --eqx: an M op compares the read (seq, the alignment's orientation) with the reference at
    `start` (fetch(lo, hi)); = and X ops are trusted."""
    eq, x, q, r = set(), set(), 0, start
    for op, n in cigar:
        if op == 0:
            a, b = seq[q:q + n].upper(), fetch(r, r + n).upper()
            for j in range(n):
                (eq if a[j] == b[j] and a[j] != "N" else x).add(q + j)
        elif op == 7:
            eq.update(range(q, q + n))
        elif op == 8:
            x.update(range(q, q + n))
        if op in (0, 1, 4, 5, 7, 8):
            q += n
        if op in (0, 2, 3, 7, 8):
            r += n
    if flip:
        eq, x = {length - 1 - i for i in eq}, {length - 1 - i for i in x}
    return eq, x


def world_tag(w):
    """'<sample>_<world>' of a world directory .../<sample>/<world>."""
    w = w.rstrip("/")
    return f"{os.path.basename(os.path.dirname(w))}_{os.path.basename(w)}"


def calls_tag(p):
    return os.path.basename(p)[:-len(".calls.tsv")] if p.endswith(".calls.tsv") else os.path.basename(p)


def keyed(read, tag, pooled):
    """Simulated read names repeat across worlds (s0000001 in every one): pooled tables key a read by its world."""
    return f"{tag}/{read}" if pooled else read


def truth_sim(rows, in_reference):
    """{read: (status, copy)} from a simulated world's reads.truth: determined when the copy is in the reference."""
    return {r["read"]: ("determined", r["copy"]) if r["copy"] in in_reference else ("off-copy", "") for r in rows}


HG002_STATUS = {"unique": "determined", "unique_elsewhere": "determined", "own_tie": "undetermined",
                "no_counterpart": "off-copy", "elsewhere": "off-copy", "poor_match": "poor_match"}


def truth_hg002(rows):
    """{read: (status, copy)} from hg002truth's truth tables (class -> status; unique_elsewhere's truth is a CHM13 copy
    outside minimap2's tie, so it is determined too)."""
    return {r["read"]: (HG002_STATUS[r["class"]], r["truth"] if HG002_STATUS[r["class"]] == "determined" else "")
            for r in rows}


MIN_READS, MIN_FAMILIES = 300, 11
BAR_H, BAR_R, BAR_Q = 0.99, 0.50, 0.98


def bars(calls, truth):
    """The pre-registered bars (PREREG_o2_window_set_2026-10-20) of the set rule on one window's calls: P power (>= MIN_READS
    scored from >= MIN_FAMILIES families), H the set holds the truth on >= BAR_H of scored reads, R it is one copy for >=
    BAR_R of them, Q those singletons are right on >= BAR_Q and at least as often as minimap2's primary on the same reads.
    A read is scored when determined with the truth among the window's copies."""
    n = dict(reads=len(calls), scored=0, truth_held=0, singletons=0, singletons_right=0, minimap2_right=0,
             minimap2_right_on_singletons=0, o2_fixes=0, o2_breaks=0, primary_stitched=0, truth_outside=0)
    fams, status = set(), {}
    for c in calls:
        st, t = truth.get(c["read"], ("missing", ""))
        if st != "determined":
            status[st] = status.get(st, 0) + 1
            continue
        if t not in c["copies"].split(","):
            n["truth_outside"] += 1
            continue
        n["scored"] += 1
        fams.add(c["family"])
        sset = set(c["o2_set"].split(","))
        mm_ok = c["minimap2"] == t
        n["truth_held"] += t in sset
        n["minimap2_right"] += mm_ok
        n["primary_stitched"] += c["primary_stitched"] == "1"
        n["o2_fixes"] += c["o2"] == t and not mm_ok
        n["o2_breaks"] += mm_ok and c["o2"] != t
        if len(sset) == 1:
            n["singletons"] += 1
            n["singletons_right"] += sset == {t}
            n["minimap2_right_on_singletons"] += mm_ok
    k = n["scored"]
    n["families"] = len(fams)
    n["status"] = status
    n["P"] = k >= MIN_READS and len(fams) >= MIN_FAMILIES
    n["H"] = k > 0 and n["truth_held"] >= BAR_H * k
    n["R"] = k > 0 and n["singletons"] >= BAR_R * k
    n["Q"] = (n["singletons"] > 0 and n["singletons_right"] >= BAR_Q * n["singletons"]
              and n["singletons_right"] >= n["minimap2_right_on_singletons"])
    n["verdict"] = "UNDERPOWERED" if not n["P"] else "PASS" if n["H"] and n["R"] and n["Q"] else "FAIL"
    return n


BAND_N, BAND_Q = 0.25, 0.98


def _arm_counts(calls, truth, pick):
    n = dict(scored=0, truth_held=0, singletons=0, singletons_right=0, minimap2_right=0)
    fams = set()
    for c in calls:
        st, t = truth.get(c["read"], ("missing", ""))
        if st != "determined" or t not in c["copies"].split(","):
            continue
        s = pick(c)
        n["scored"] += 1
        fams.add(c["family"])
        n["truth_held"] += t in s
        n["minimap2_right"] += c["minimap2"] == t
        if len(s) == 1:
            n["singletons"] += 1
            n["singletons_right"] += s == {t}
    n["families"] = len(fams)
    return n


def band_bars(calls, truth):
    """PREREG_o2_set_null_band: on one band's calls, the set against two nulls (keep-all: every window copy; primary:
    minimap2's copy), and the bars P (MIN_READS / MIN_FAMILIES), H (set holds the truth >= BAR_H), N (set is one copy on
    >= BAND_N of scored reads), Q (singletons right >= BAND_Q and more than minimap2's primary on those reads)."""
    arms = {"set": _arm_counts(calls, truth, lambda c: set(c["o2_set"].split(","))),
            "keep_all": _arm_counts(calls, truth, lambda c: set(c["copies"].split(","))),
            "primary": _arm_counts(calls, truth, lambda c: {c["minimap2"]})}
    mm_on_sing = _arm_counts([c for c in calls if len(c["o2_set"].split(",")) == 1], truth, lambda c: {c["minimap2"]})
    st = arms["set"]
    k = st["scored"]
    P = k >= MIN_READS and st["families"] >= MIN_FAMILIES
    H = k > 0 and st["truth_held"] >= BAR_H * k
    N = k > 0 and st["singletons"] >= BAND_N * k
    Q = st["singletons"] > 0 and st["singletons_right"] >= BAND_Q * st["singletons"] and st["singletons_right"] > mm_on_sing["singletons_right"]
    arms["minimap2_right_on_set_singletons"] = mm_on_sing["singletons_right"]
    arms["bars"] = dict(P=P, H=H, N=N, Q=Q, verdict="UNDERPOWERED" if not P else "PASS" if H and N and Q else "FAIL")
    return arms


# ---------------------------------------------------------------------------------------------------------------- commands

def cmd_prepare(bam_p, copies_p, qual_p, genome_p, which, out, windows):
    """qual_p / which: a qual.json and its half ('dev' / 'heldout'), or '-' '-' for every family of the copies table."""
    import pysam
    P, S = PC._sibling("o2psv"), PC._sibling("summary")
    rows = PC.read_copies(copies_p)
    qual = json.load(open(qual_p)) if qual_p != "-" else {r["family"]: dict(half=which) for r in rows}
    exons = {r["copy"]: r["exons"] for r in rows}
    idx = C.span_index(P.load_copies(copies_p))
    fa = pysam.FastaFile(genome_p)
    wmax = max(windows)

    def info(rec):
        cands = C.overlapping(idx, rec.reference_name, rec.reference_start, rec.reference_end)
        bl = rec.get_blocks()
        c = S.attribute(bl, P._strand(rec), cands) if cands else None
        st = stitched(bl, {x[0]: exons[x[0]] for x in cands}) if cands else False
        return c, st

    alns = collections.defaultdict(list)
    with pysam.AlignmentFile(bam_p) as bam:                      # pass 1: every read's alignments, light
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_supplementary or not rec.has_tag("ms"):
                continue
            c, st = info(rec)
            alns[rec.query_name].append(dict(ms=rec.get_tag("ms"), copy=c, prim=not rec.is_secondary, st=st))
    dom = {}
    for q, al in alns.items():
        per = {w: window_copies(al, w)[0] for w in windows}
        hit = [w for w in windows if per[w]]
        if hit and qual.get(next(iter(per[hit[0]])).rsplit(".", 1)[0], {}).get("half") == which:
            dom[q] = per
    del alns
    print(f"[o2window] pass 1: {len(dom)} {which} reads tied in some window", flush=True)
    recs, seqs = collections.defaultdict(list), {}
    with pysam.AlignmentFile(bam_p) as bam:                      # pass 2: their alignments in full
        for rec in bam.fetch(until_eof=True):
            q = rec.query_name
            if q not in dom or rec.is_unmapped or rec.is_supplementary or not rec.has_tag("ms"):
                continue
            c, st = info(rec)
            recs[q].append(dict(ms=rec.get_tag("ms"), AS=rec.get_tag("AS"), copy=c, prim=not rec.is_secondary, st=st,
                                start=rec.reference_start, cig=tuple(rec.cigartuples), rev=rec.is_reverse,
                                ch=rec.reference_name))
            if not rec.is_secondary:
                seqs[q] = (rec.query_sequence, rec.is_reverse)

    def fetcher(chrom):
        return lambda lo, hi: fa.fetch(chrom, lo, hi)

    n = collections.Counter()
    with open(f"{out}.calls.tsv", "w") as fh, open(f"{out}.names", "w") as nh:
        fh.write("read\tfamily\twindow\tcopies\tminimap2\tprimary_stitched\to2\to2_set\tcandidates\n")
        for q in sorted(dom):
            nh.write(q + "\n")
            al = recs[q]
            seq, prev = seqs[q]
            cache = {}
            prim = next(a for a in al if a["prim"])
            best = max(a["ms"] for a in al)
            for w in windows:
                cs = dom[q][w]
                if not cs:
                    continue
                rep = {}                                        # per copy, its best non-stitched alignment in the window
                for a in al:
                    if a["ms"] > best - w and a["copy"] in cs and not a["st"]:
                        if a["copy"] not in rep or a["ms"] > rep[a["copy"]]["ms"]:
                            rep[a["copy"]] = a
                mm = prim["copy"]
                if len(rep) >= 2:
                    spans = {c: (C.read_span(a["cig"]), a["rev"] != prev) for c, a in rep.items()}
                    win = C.common_window(spans, len(seq))
                    ed = {}
                    for c, a in rep.items():
                        k = (a["start"], a["cig"], a["ch"], win[c])
                        if k not in cache:
                            cache[k] = PC.window_edits(C.alignment_events(a["start"], a["cig"], seq if a["rev"] == prev
                                                                          else PC.rc(seq), fetcher(a["ch"]), win[c])[0])
                        ed[c] = cache[k]
                    pick = PC.fewest(ed, mm)
                    cols = {c: match_columns_ref(a["cig"], a["rev"] != prev, len(seq), seq if a["rev"] == prev else PC.rc(seq),
                                                 a["start"], fetcher(a["ch"])) for c, a in rep.items()}
                    st_ = pareto_set(cols, PC.window_primary(spans, len(seq)))
                elif len(rep) == 1:
                    pick = next(iter(rep))
                    st_ = set(rep)
                else:
                    pick, st_ = mm, {mm}
                n[(w, "reads")] += 1
                n[(w, "O2 moves")] += pick != mm
                n[(w, "set of 1")] += len(st_) == 1
                fh.write(f"{q}\t{mm.rsplit('.', 1)[0]}\t{w}\t{','.join(sorted(cs))}\t{mm}\t{int(prim['st'])}\t{pick}\t"
                         f"{','.join(sorted(st_))}\t{','.join(sorted(rep))}\n")
    with open(f"{out}.domain.tsv", "w") as fh:                  # for o2personal.py truth (o2correct donor's columns)
        fh.write("read\tfamily\tcorrector\tminimap2\treason\tcandidates\tmargin\tscores\n")
        for q in sorted(dom):
            prim = next(a for a in recs[q] if a["prim"])
            cs = set().union(*(v for v in dom[q].values() if v))
            fh.write(f"{q}\t{prim['copy'].rsplit('.', 1)[0]}\t{prim['copy']}\t{prim['copy']}\twindow\t{','.join(sorted(cs))}\t0\t\n")
    print(f"[o2window] -> {out}.calls.tsv: " + ", ".join(f"W{w} {k} {v}" for (w, k), v in sorted(n.items())), flush=True)
    return 0


def cmd_score(calls_p, truth_p, out_p):
    truth = {r["read"]: r for r in csv.DictReader(open(truth_p), delimiter="\t")}
    by = collections.defaultdict(collections.Counter)
    sizes = collections.defaultdict(lambda: collections.defaultdict(list))
    prev_w = {}
    for r in csv.DictReader(open(calls_p), delimiter="\t"):
        t, w = truth[r["read"]], int(r["window"])
        cs = set(r["copies"].split(","))
        band = "first tied here" if r["read"] not in prev_w else "tied in a smaller window"
        prev_w.setdefault(r["read"], w)
        for key in (f"W{w}", f"W{w} {band}"):
            n = by[key]
            n["reads"] += 1
            if t["status"] != "determined":
                n[t["status"]] += 1
                continue
            if t["truth"] not in cs:
                n["truth outside the window"] += 1
                continue
            n["scored"] += 1
            n["primary stitched"] += r["primary_stitched"] == "1"
            mm_ok, o2_ok = r["minimap2"] == t["truth"], r["o2"] == t["truth"]
            n["minimap2 wrong"] += not mm_ok
            n["O2 fixes"] += o2_ok and not mm_ok
            n["O2 breaks"] += mm_ok and not o2_ok
            st = set(r["o2_set"].split(","))
            n["set holds the truth"] += t["truth"] in st
            n["set of 1"] += len(st) == 1
            n["set of 1, right"] += st == {t["truth"]}
            sizes[key]["pool"].append(len(cs))
            sizes[key]["set"].append(len(st))
    res = {}
    for key, n in sorted(by.items(), key=lambda kv: (int(kv[0].split()[0][1:]), kv[0])):
        d = dict(n)
        if sizes[key]["set"]:
            d["mean pool"] = round(sum(sizes[key]["pool"]) / len(sizes[key]["pool"]), 3)
            d["mean set"] = round(sum(sizes[key]["set"]) / len(sizes[key]["set"]), 3)
        res[key] = d
        print(f"[o2window] {key}: " + ", ".join(f"{k} {v}" for k, v in d.items()), flush=True)
    json.dump(res, open(out_p, "w"), indent=1)
    return 0


def cmd_truth_ms(domain_p, copies_p, orth_p, chrmap_p, out_p, bams):
    """truth_ms() for every read of a domain TSV (o2correct donor's columns) from its alignments to the animal's two
    assemblies, attributed as o2personal.cmd_truth attributes them (repeated here so the frozen o2personal stays as run)."""
    import pysam
    P, S = PC._sibling("o2psv"), PC._sibling("summary")
    dom = PC.read_domain(domain_p)
    refc = PC.reference_contigs(chrmap_p)
    rows = PC.read_copies(copies_p)
    idx_ref = C.span_index(P.load_copies(copies_p))
    idx_alt = C.span_index(PC.alt_rows(PC.read_orthologs(orth_p), {r["copy"]: r["strand"] for r in rows}))
    alns = collections.defaultdict(list)
    for h, (bam_p, fai_p) in bams.items():
        with pysam.AlignmentFile(bam_p) as bam:
            names = PC.names_by_length({x["SN"]: x["LN"] for x in bam.header.to_dict()["SQ"]}, PC._fai(fai_p))
            for rec in bam.fetch(until_eof=True):
                q = rec.query_name
                if q not in dom or rec.is_unmapped or rec.is_supplementary or not rec.has_tag("ms"):
                    continue
                cm = names.get(rec.reference_name)
                if cm in refc:
                    role, cands = "ref", C.overlapping(idx_ref, refc[cm][0], rec.reference_start, rec.reference_end)
                else:
                    role = "alt"
                    cands = C.overlapping(idx_alt, cm, rec.reference_start, rec.reference_end) if cm else []
                c = S.attribute(rec.get_blocks(), P._strand(rec), cands) if cands else None
                alns[q].append((rec.get_tag("ms"), c, role))
    n = collections.Counter()
    with open(out_p, "w") as fh:
        fh.write("read\tfamily\tstatus\ttruth\torigin\tbest_ms\trunner_up\n")
        for q, (f, *_rest) in sorted(dom.items()):
            st, t, o = truth_ms(alns.get(q, []))
            best = max((m for m, _, _ in alns[q]), default="")
            ru = max((m for m, c, _ in alns[q] if c != t), default="") if st == "determined" else ""
            fh.write(f"{q}\t{f}\t{st}\t{t or ''}\t{o or ''}\t{best}\t{ru}\n")
            n[st if st != "determined" else f"determined {o}"] += 1
    print(f"[o2window] ms truth -> {out_p}: " + ", ".join(f"{k} {v}" for k, v in sorted(n.items())), flush=True)
    return 0


def cmd_truth_table(kind, out_p, paths):
    """A truth table (read, status, truth) for the scorer: kind 'sim' from WORLD dirs (reads.truth.E0.tsv, copies.tsv),
    'hg002' from hg002truth tables."""
    out = {}
    for p in paths:
        if kind == "sim":
            inref = {r["copy"] for r in csv.DictReader(open(f"{p}/copies.tsv"), delimiter="\t") if r["in_reference"] == "1"}
            t = truth_sim(list(csv.DictReader(open(f"{p}/reads.truth.E0.tsv"), delimiter="\t")), inref)
            out.update({keyed(q, world_tag(p), True): v for q, v in t.items()})
        else:
            out.update(truth_hg002(list(csv.DictReader(open(p), delimiter="\t"))))
    with open(out_p, "w") as fh:
        fh.write("read\tstatus\ttruth\n")
        for q, (st, t) in sorted(out.items()):
            fh.write(f"{q}\t{st}\t{t}\n")
    print(f"[o2window] truth table ({kind}) -> {out_p}: {len(out)} reads", flush=True)
    return 0


def cmd_bars(out_p, calls_ps, truth_ps, window):
    """The pre-registered bars on one window, pooled over the calls tables (one per world or sample), plus the bands
    'first tied within 10 / 20 / 30' as descriptive strata."""
    truth = {}
    for p in truth_ps:
        for r in csv.DictReader(open(p), delimiter="\t"):
            truth[r["read"]] = (r["status"], r["truth"])
    calls, first = [], {}
    pooled = len(calls_ps) > 1                                  # several worlds: reads keyed by world (keyed())
    for p in calls_ps:
        for r in csv.DictReader(open(p), delimiter="\t"):
            r["read"] = keyed(r["read"], calls_tag(p), pooled)
            first.setdefault(r["read"], int(r["window"]))
            if int(r["window"]) == window:
                calls.append(r)
    res = {"all": bars(calls, truth)}
    for lo, hi in ((0, 10), (10, 20), (20, 30), (30, 50)):
        band = [c for c in calls if lo < first[c["read"]] <= hi or (lo == 0 and first[c["read"]] <= hi)]
        if band:
            res[f"first tied within {hi}"] = bars(band, truth)
    json.dump(res, open(out_p, "w"), indent=1)
    for k, b in res.items():
        print(f"[o2window] W{window} {k}: scored {b['scored']} / {b['families']} families, minimap2 right {b['minimap2_right']}, "
              f"held {b['truth_held']}, singletons {b['singletons']} (right {b['singletons_right']}, minimap2 right there "
              f"{b['minimap2_right_on_singletons']}), O2 pick +{b['o2_fixes']}/-{b['o2_breaks']}, status {b['status']}, "
              f"outside {b['truth_outside']} | P {b['P']} H {b['H']} R {b['R']} Q {b['Q']} -> {b['verdict']}", flush=True)
    return 0


def cmd_band(out_p, calls_ps, truth_ps):
    """band_bars on the near-tie band (reads first tied within 10 points), plus the arms on the 10-20 and 20-30 bands and
    the whole 30-point window, pooled over several calls tables (reads keyed by world)."""
    truth = {}
    for p in truth_ps:
        for r in csv.DictReader(open(p), delimiter="\t"):
            truth[r["read"]] = (r["status"], r["truth"])
    rows, first = [], {}
    pooled = len(calls_ps) > 1
    for p in calls_ps:
        for r in csv.DictReader(open(p), delimiter="\t"):
            r["read"] = keyed(r["read"], calls_tag(p), pooled)
            first.setdefault(r["read"], int(r["window"]))
            if r["window"] == "30":
                rows.append(r)
    res = {}
    for name, lo, hi in (("band <= 10 (under test)", 0, 10), ("band 10-20", 10, 20), ("band 20-30", 20, 30), ("window 30, all", 0, 30)):
        sub = [r for r in rows if lo < first[r["read"]] <= hi]
        res[name] = band_bars(sub, truth)
        b = res[name]
        print(f"[o2window] {name}: " + "; ".join(f"{a} scored {b[a]['scored']} held {b[a]['truth_held']} singletons {b[a]['singletons']} "
              f"right {b[a]['singletons_right']}" for a in ("set", "keep_all", "primary"))
              + f"; families {b['set']['families']}; minimap2 right on the set's singletons {b['minimap2_right_on_set_singletons']}; "
              f"bars {b['bars']}", flush=True)
    json.dump(res, open(out_p, "w"), indent=1)
    return 0


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    cmd, a = argv[0], argv[1:]
    if cmd == "prepare":
        ws = [int(x) for x in a[a.index("--windows") + 1].split(",")] if "--windows" in a else [10, 20, 30, 50]
        return cmd_prepare(*a[:6], ws)
    if cmd == "score":
        return cmd_score(*a[:3])
    if cmd == "truth-ms":
        return cmd_truth_ms(*a[:5], PC._pairs(a, "--bams"))
    if cmd in ("truth-sim", "truth-hg002"):
        return cmd_truth_table(cmd[6:], a[0], a[1:])
    if cmd == "band":
        i, j = a.index("--calls"), a.index("--truth")
        return cmd_band(a[0], a[i + 1:j], a[j + 1:])
    if cmd == "bars":
        win = int(a[a.index("--window") + 1]) if "--window" in a else 30
        if "--window" in a:
            k = a.index("--window")
            a = a[:k] + a[k + 2:]
        i, j = a.index("--calls"), a.index("--truth")
        return cmd_bars(a[0], a[i + 1:j], a[j + 1:], win)
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
