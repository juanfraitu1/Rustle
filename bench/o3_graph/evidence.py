#!/usr/bin/env python3
"""What O2's evidence rests on (2026-10-16 goal: can PSVs include SUNs; how much do the CIGAR, NM and de tags matter).

    evidence.py suns OUT.json WORLD ...     the current O2's PSV columns (copy_assign --dump-psv, WORLD/cur_dump/A.psv_*.tsv) classed
                                           as SUNs (one copy differs, all others agree), pair or shared columns; every assigned
                                           read's support split by class, and its accuracy against truth by class
    evidence.py star OUT.json WORLD ...     the current O2's read-star proofs (copy_assign --dump-star, WORLD/cur_dump/S.star_reads.tsv):
                                           every assigned read's decisive columns, classed among its candidates as its copy's own
                                           SUN, another copy's SUN it rules out, a two-copy pair or a shared column; accuracy by
                                           evidence and by the read's minimap2 `de`, against truth
    evidence.py o2ends OUT.json WORLD ...   the integrated O2's per-copy distance two ways -- minimap2's NM plus every soft-clipped
                                           read base (o2graph._paf_dist, what O2 runs) and a SEMI-GLOBAL distance (each clipped
                                           read end aligned against the path bases next to the aligned core; the read fully
                                           aligned, the path free beyond) -- and how the margin-2 decision changes, against truth

A clipped end is part of the read: counting its bases as edits (the round-8 error: 4 clipped bases, 1 of them a sequencing
error, made 1 vs 5) inflates the distance of exactly the copies whose path disagrees with an error near the read's end.
"""
import collections
import csv
import importlib.util
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXT_CAP, EXT_BAND = 50, 30
COMP = str.maketrans("ACGTN", "TGCAN")


def _load(name):
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def end_cost(p, t, anchored_at_end):
    """Edit distance of the whole read end `p` against path bases `t`, anchored where the end meets the aligned core (a prefix
    at t's end, a suffix at t's start) and free on the far side; read bases beyond the path cost one each. An end longer than
    EXT_CAP costs its length (what clipping counts)."""
    if len(p) > EXT_CAP:
        return len(p)
    if anchored_at_end:
        p, t = p[::-1], t[::-1]
    prev = list(range(len(t) + 1))
    for i in range(1, len(p) + 1):
        cur = [i] + [0] * len(t)
        a = p[i - 1]
        for j in range(1, len(t) + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a != t[j - 1]))
        prev = cur
    return min(prev)


def _paf(line):
    f = line.rstrip("\n").split("\t")
    nm = int(next(x for x in f[12:] if x.startswith("NM:i:"))[5:])
    return f, nm, int(f[1]), int(f[2]), int(f[3]), f[4], int(f[7]), int(f[8])


def clip_distance(line):
    """O2's distance (o2graph._paf_dist without `ends`): NM plus the read's unaligned bases."""
    _, nm, qlen, qs, qe, _, _, _ = _paf(line)
    return nm + qs + qlen - qe


def semiglobal_distance(line, read, target):
    """NM plus each clipped read end's edits against the path bases next to the aligned core."""
    _, nm, qlen, qs, qe, strand, ts, te = _paf(line)
    if strand == "-":
        read, qs, qe = read[::-1].translate(COMP), qlen - qe, qlen - qs
    pre, suf = read[:qs], read[qe:]
    return (nm + end_cost(pre, target[max(0, ts - len(pre) - EXT_BAND):ts], True)
            + end_cost(suf, target[te:te + len(suf) + EXT_BAND], False))


NO_ALLELE = set("-N.?")


def column_classes(copies):
    """Per PSV column of a family ({copy: allele string}, one character per column): None (no difference), ("sun", c) when
    copy c alone differs and every other copy agrees (c's singly unique nucleotide), ("pair", None) when only two copies
    have an allele, ("shared", None) otherwise (an allele several copies share)."""
    n = len(next(iter(copies.values())))
    out = []
    for j in range(n):
        al = {c: s[j] for c, s in copies.items() if s[j] not in NO_ALLELE}
        counts = collections.Counter(al.values())
        if len(counts) < 2:
            out.append(None)
        elif len(al) == 2:
            out.append(("pair", None))
        else:
            lone = [b for b, k in counts.items() if k == 1]
            if len(counts) == 2 and len(lone) == 1:
                out.append(("sun", next(c for c, b in al.items() if b == lone[0])))
            else:
                out.append(("shared", None))
    return out


def support(read, copies, a, cls):
    """A read's bases at the PSV columns (no allele: '.', '?', '-', 'N') against its called copy a: columns where it carries
    a's own SUN, columns that rule out another copy's SUN (it matches the shared majority), shared columns, and conflicts."""
    out = dict(sun=0, shared=0, against_sun=0, conflict=0)
    for j, b in enumerate(read):
        if b in NO_ALLELE or cls[j] is None or copies[a][j] in NO_ALLELE:
            continue
        if b != copies[a][j]:
            out["conflict"] += 1
        elif cls[j][0] == "sun":
            out["sun" if cls[j][1] == a else "against_sun"] += 1
        else:
            out["shared"] += 1
    return out


def star_support(columns, assigned):
    """One read-star proof (`readpos:readbase:candidate bases`, comma-joined; '.' = the candidate has no base) against the
    assigned candidate (its index in the proof's candidate list)."""
    out = dict(sun=0, shared=0, pair=0, against_sun=0, conflict=0, cols=0)
    for col in columns.split(","):
        if not col:
            continue
        _, rb, cands = col.split(":")
        out["cols"] += 1
        if cands[assigned] in NO_ALLELE:
            continue
        if rb != cands[assigned]:
            out["conflict"] += 1
            continue
        present = {i: b for i, b in enumerate(cands) if b not in NO_ALLELE}
        counts = collections.Counter(present.values())
        lone = [b for b, k in counts.items() if k == 1]
        if len(present) == 2:
            out["pair"] += 1
        elif len(counts) == 2 and len(lone) == 1:
            out["sun" if present[assigned] == lone[0] else "against_sun"] += 1
        else:
            out["shared"] += 1
    return out


def gap_compressed_distance(line, read, target):
    """Like minimap2's `de`, counted not divided: mismatches plus indel EVENTS (an indel of any length counts once) inside the
    aligned core, plus each clipped read end's semi-global edits (end_cost)."""
    f, nm, qlen, qs, qe, strand, ts, te = _paf(line)
    cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")
    num, events, indel_bases = "", 0, 0
    for ch in cg:
        if ch.isdigit():
            num += ch
        else:
            if ch in "ID":
                events += 1
                indel_bases += int(num)
            num = ""
    if strand == "-":
        read, qs, qe = read[::-1].translate(COMP), qlen - qe, qlen - qs
    pre, suf = read[:qs], read[qe:]
    return (nm - indel_bases + events + end_cost(pre, target[max(0, ts - len(pre) - EXT_BAND):ts], True)
            + end_cost(suf, target[te:te + len(suf) + EXT_BAND], False))


def _fasta(p):
    out, name = {}, None
    for ln in open(p):
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            out[name] = []
        else:
            out[name].append(ln.strip())
    return {k: "".join(v) for k, v in out.items()}


VARIANTS = ("clip", "semiglobal", "gapcomp")


def o2ends_world(w):
    """Re-run the integrated O2's path alignments; decide each read (o2graph.decide, margin 2, delta gate) on three distances --
    NM + clipped bases (what O2 runs), semi-global (clipped ends aligned), gap-compressed (indel events, not bases) -- and
    score each against truth: {variant: Counter(correct, wrong, abstain, unexplained)} over the reads of each family's own
    copies, plus any call a family makes on another family's read (wrong)."""
    O = _load("o2graph")
    d = f"{w}/graph/o2paths.union.pad.integ"
    truth = {r["read"]: r for r in csv.DictReader(open(f"{w}/reads.truth.E0.tsv"), delimiter="\t")}
    iso = json.load(open(f"{w}/isoforms.json"))
    seeded = {(f, c["idx"]): c["copy"] for f, cps in iso.items() for c in cps}
    sj = json.load(open(f"{w}/graph/seeded.json"))
    truth_copies = list(csv.DictReader(open(f"{w}/copies.tsv"), delimiter="\t"))
    for f, v in sj.items():
        for i, cp in O.map_new_copies([c for c in v["copies"] if c["origin"] == "new"], truth_copies).items():
            seeded[(f, i)] = cp
    tab = {v: collections.Counter() for v in VARIANTS}
    for fn in sorted(os.listdir(d)):
        if not fn.endswith(".tx.fa"):
            continue
        fam = fn[:-6]
        tx, reads = _fasta(f"{d}/{fam}.tx.fa"), _fasta(f"{d}/{fam}.reads.fa")
        res = subprocess.run(["minimap2", "-c", "-x", "map-hifi", "-N", str(max(200, 2 * len(tx))), "-p", "0.3", "--secondary=yes",
                              f"{d}/{fam}.tx.fa", f"{d}/{fam}.reads.fa"], capture_output=True, text=True, check=True).stdout
        dist = {v: collections.defaultdict(dict) for v in VARIANTS}
        for ln in res.splitlines():
            f = ln.split("\t")
            c, q = f[5].split("|")[0], f[0]
            for v, fn_ in (("clip", lambda: clip_distance(ln)), ("semiglobal", lambda: semiglobal_distance(ln, reads[q], tx[f[5]])),
                           ("gapcomp", lambda: gap_compressed_distance(ln, reads[q], tx[f[5]]))):
                x = fn_()
                dist[v][q][c] = min(x, dist[v][q].get(c, 10 ** 9))
        for q in reads:
            t = truth.get(q)
            if t is None:
                continue
            own_family = t["family"] == fam
            for v in VARIANTS:
                c, reason = O.decide(dist[v].get(q, {}), len(reads[q]))
                if c is None:
                    if own_family:
                        tab[v]["unexplained" if reason == "unexplained" else "abstain"] += 1
                else:
                    ok = own_family and seeded.get((fam, int(c[4:]))) == t["copy"]
                    tab[v]["correct" if ok else "wrong"] += 1
    return tab


def cmd_o2ends(out, worlds):
    res = {}
    for w in worlds:
        tab = o2ends_world(w)
        res[w] = {v: dict(c) for v, c in tab.items()}
        print(w.split("/o3_graph/")[-1] + ": " + " | ".join(f"{v} " + " ".join(f"{k} {c.get(k, 0)}" for k in
                                                       ("correct", "wrong", "abstain", "unexplained")) for v, c in tab.items()),
              flush=True)
    json.dump(res, open(out, "w"), indent=1)
    return 0


def suns_world(w):
    """Column classes per family and assigned reads' accuracy by the kind of PSV support they carry."""
    truth = {r["read"]: r["copy"] for r in csv.DictReader(open(f"{w}/reads.truth.E0.tsv"), delimiter="\t")}
    copies, tid = collections.defaultdict(dict), {}
    for r in csv.DictReader(open(f"{w}/cur_dump/A.psv_copies.tsv"), delimiter="\t"):
        copies[r["family_id"]][int(r["copy_index"])] = r["alleles"]
        tid[(r["family_id"], int(r["copy_index"]))] = r["copy_tid"]
    cls = {f: column_classes(c) for f, c in copies.items() if c}
    cols = collections.Counter(x[0] if x else "none" for v in cls.values() for x in v)
    groups = collections.defaultdict(collections.Counter)
    for r in csv.DictReader(open(f"{w}/cur_dump/A.psv_reads.tsv"), delimiter="\t"):
        if r["status"] != "assigned" or r["read_name"] not in truth or r["family_id"] not in cls:
            continue
        a = int(r["assigned_copy"])
        sp = support(r["alleles"], copies[r["family_id"]], a, cls[r["family_id"]])
        key = ("own SUN" if sp["sun"] else "rules out others' SUNs only" if sp["against_sun"] and not sp["shared"]
               else "shared/pair only" if sp["shared"] and not sp["against_sun"] else "mixed, no own SUN"
               if sp["shared"] else "none")
        ok = tid.get((r["family_id"], a)) == truth[r["read_name"]]
        groups[key]["correct" if ok else "wrong"] += 1
        groups[key]["conflict_reads"] += sp["conflict"] > 0
    return dict(columns=dict(cols), groups={k: dict(v) for k, v in groups.items()})


def cmd_suns(out, worlds):
    res = {w: suns_world(w) for w in worlds}
    for w, v in res.items():
        print(w, "columns", v["columns"], "| assigned reads by support:",
              {k: f"{g.get('correct', 0)}/{g.get('correct', 0) + g.get('wrong', 0)}" for k, g in v["groups"].items()}, flush=True)
    json.dump(res, open(out, "w"), indent=1)
    return 0


def star_world(w):
    import pysam
    truth = {r["read"]: r["copy"] for r in csv.DictReader(open(f"{w}/reads.truth.E0.tsv"), delimiter="\t")}
    iso = json.load(open(f"{w}/isoforms.json"))
    cat = {(f, c["idx"]): c["copy"] for f, cps in iso.items() for c in cps}
    de = {}
    with pysam.AlignmentFile(f"{w}/E0.bam") as bam:
        for r in bam.fetch(until_eof=True):
            if not (r.is_unmapped or r.is_secondary or r.is_supplementary):
                de[r.query_name] = r.get_tag("de")
    rows = []
    for r in csv.DictReader(open(f"{w}/cur_dump/S.star_reads.tsv"), delimiter="\t"):
        if r["status"] != "assigned" or r["read_name"] not in truth or not r["catalog_copy_idx"].isdigit():
            continue
        cands = r["candidates"].split(",")
        if r["assigned_copy"] not in cands:
            continue
        sp = star_support(r["columns"], cands.index(r["assigned_copy"]))
        sp["ok"] = cat.get((r["family_id"], int(r["catalog_copy_idx"]))) == truth[r["read_name"]]
        sp["de"] = de.get(r["read_name"], float("nan"))
        rows.append(sp)
    return rows


def _rate(rows):
    n = len(rows)
    return f"{sum(not r['ok'] for r in rows)}/{n}" if n else "0/0"


def cmd_star(out, worlds):
    res = {}
    for w in worlds:
        rows = star_world(w)
        sup = lambda r: r["sun"] + r["against_sun"] + r["pair"] + r["shared"]  # noqa: E731
        by_ev = {"own SUN": [r for r in rows if r["sun"]], "no own SUN": [r for r in rows if not r["sun"]]}
        by_n = {k: [r for r in rows if lo <= sup(r) <= hi] for k, (lo, hi) in
                {"1": (1, 1), "2": (2, 2), "3-5": (3, 5), "6+": (6, 10 ** 9), "0": (0, 0)}.items()}
        q = sorted(r["de"] for r in rows)
        cuts = [q[int(len(q) * f)] for f in (0.25, 0.5, 0.75)] if q else []
        by_de = {f"de q{i + 1}": [r for r in rows if (i == 0 or r["de"] > cuts[i - 1]) and (i == 3 or r["de"] <= cuts[i])]
                 for i in range(4)} if cuts else {}
        res[w] = dict(assigned=len(rows), wrong=sum(not r["ok"] for r in rows),
                      by_evidence={k: _rate(v) for k, v in by_ev.items()}, by_support_cols={k: _rate(v) for k, v in by_n.items()},
                      by_de={k: _rate(v) for k, v in by_de.items()}, de_cuts=cuts,
                      wrong_calls=[{k: v for k, v in r.items() if k != "ok"} for r in rows if not r["ok"]],
                      cols_classes={k: sum(r[k] for r in rows) for k in ("sun", "against_sun", "pair", "shared", "conflict")})
        v = res[w]
        print(f"{w}: wrong/assigned {v['wrong']}/{v['assigned']} | by evidence {v['by_evidence']} | by supporting columns "
              f"{v['by_support_cols']} | by de {v['by_de']} | supporting columns {v['cols_classes']}", flush=True)
    json.dump(res, open(out, "w"), indent=1)
    return 0


def main(argv):
    if argv[0] == "star":
        return cmd_star(argv[1], argv[2:])
    if argv[0] == "suns":
        return cmd_suns(argv[1], argv[2:])
    if argv[0] == "o2ends":
        return cmd_o2ends(argv[1], argv[2:])
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
