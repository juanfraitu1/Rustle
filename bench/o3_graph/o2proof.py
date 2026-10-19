#!/usr/bin/env python3
"""The O2 scorecard by DECIDABILITY (user, 2026-10-17: "can we use these to prove that O2 works?"): every catalogued read
is classed by what its own sequence can tell, then each assigner's outcome is tallied by that class and by MAPQ.

    o2proof.py score OUT.json WORLD ... [--o2 NAME] [--tied TIED.tsv]
        per world and pooled; NAME defaults to o2.integv.m2.tsv. With --tied (hg002truth.py tied: read, copies), also the
        user's scoped systems `current_tied` / `integrated_tied` (2026-10-18: "O2 can only work for tied AS reads,
        everything else is already resolved by minimap"): a tied read's call stands only on a copy it ties on, else it is
        an abstention; an untied read takes minimap2's primary.

The oracle:
- A read's alignments in WORLD/E0.bam (primary and secondary; supplementary ones are pieces) are placed on the copy whose
  locus each overlaps most, or, off every copy, on a 10-kb bin of the genome. A locus's COST = its best alignment's NM plus
  the read bases that alignment leaves unaligned (soft clips, partial alignments): a copy lacking an exon the read carries
  costs at least that exon.
- O2 decides among the copies of the read's family (its candidate set, in the reference genome, catalogued or not):
  margin = the best family rival's cost minus the true copy's.
- Classes:
  - `off_family`: a locus OFF the family explains the read at least as well (another family's copy, or a locus that is no
    copy: a parent gene, an unlisted paralog; round 13's UBE2Q2P reads, set 9's 49 USP17L genes outside the family). Which
    loci form the family is O1's question, so these reads are classed apart.
  - `undecidable`: margin <= 0. A family rival explains the read at least as well, so no sequence-based assigner can
    prefer the truth; placing the read is a guess.
  - `weak`: margin 1. One edit, as a single sequencing error makes.
  - `decidable`: margin >= 2 (O2's own margin), or no family rival aligned at all.
  - `true_unaligned`: the aligner never placed the read on its true copy.

The assigners:
- `mm2`: the primary alignment's copy, always.
- `mm2_q1`: the primary's copy when its MAPQ >= 1, else abstain.
- `current`: copy_assign's O2 (cur/A.assignments.tsv).
- `integrated`: the O1 -> O2 table (graph/NAME on graph/seeded.json, as summary.py e2e maps it).
An outcome is `correct` (the read's own copy), `wrong` (another copy, another family, or a locus that is no copy) or
`abstain`.

Scope: the reads of catalogued copies (truth kind present or silent), which is O2's domain. Absent and uncatalogued copies
are O3's.
"""
import bisect
import collections
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CLASSES = ("undecidable", "weak", "decidable", "off_family", "true_unaligned")
METHODS = ("mm2", "mm2_q1", "current", "integrated")
OUTCOMES = ("correct", "wrong", "abstain")


def _load(name):
    import importlib.util
    spec = importlib.util.spec_from_file_location(f"o3_graph_{name}", os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def locus_index(loci):
    """{chrom: sorted [(start, end, copy)]} of copy loci (0-based half-open)."""
    idx = collections.defaultdict(list)
    for c in loci:
        idx[c["chrom"]].append((int(c["start"]), int(c["end"]), c["copy"]))
    for v in idx.values():
        v.sort()
    return idx


def copy_of(idx, chrom, s, e):
    """The copy whose locus overlaps [s, e) most, or None."""
    v = idx.get(chrom, ())
    best, ov = None, 0
    for a, b, cp in v[:bisect.bisect_left(v, (e,))]:
        o = min(b, e) - max(a, s)
        if o > ov:
            best, ov = cp, o
    return best


def locus_key(idx, chrom, s, e, bin_bp=10000):
    """An alignment's locus: the copy it overlaps most, else its genome bin (chrom, start // bin_bp)."""
    return copy_of(idx, chrom, s, e) or (chrom, s // bin_bp)


def copy_costs(alns, read_len):
    """{copy: its best alignment's NM + the read bases it leaves unaligned} from [(copy, NM, aligned query bases)]."""
    out = {}
    for cp, nm, aligned in alns:
        if cp is None:
            continue
        c = nm + max(0, read_len - aligned)
        if cp not in out or c < out[cp]:
            out[cp] = c
    return out


def decidability(costs, true, family):
    """(class, margin): `off_family` when a locus off `family` costs no more than the true copy (margin = that locus's cost
    minus the true copy's), else the best family rival's cost minus the true copy's."""
    if true not in costs:
        return "true_unaligned", None
    off = [v for k, v in costs.items() if k not in family]
    if off and min(off) <= costs[true]:
        return "off_family", min(off) - costs[true]
    rival = [v for k, v in costs.items() if k != true and k in family]
    if not rival:
        return "decidable", None
    m = min(rival) - costs[true]
    return ("undecidable" if m <= 0 else "weak" if m == 1 else "decidable"), m


def outcome(call, true):
    if call is None:
        return "abstain"
    return "correct" if call == true else "wrong"


def scorecard(rows, methods):
    """{method: {class: {stratum: {correct, wrong, abstain}}}}, strata `all`, `mapq0` (the primary's MAPQ is 0) and `mapq1`."""
    card = {m: {c: {s: dict.fromkeys(OUTCOMES, 0) for s in ("all", "mapq0", "mapq1")} for c in CLASSES} for m in methods}
    for r in rows:
        for m in methods:
            o = outcome(r["calls"].get(m), r["copy"])
            for s in ("all", "mapq0" if r["mapq0"] else "mapq1"):
                card[m][r["cls"]][s][o] += 1
    return card


def read_lengths(fq):
    out = {}
    with open(fq) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().rstrip("\n")
            fh.readline()
            fh.readline()
            out[h[1:].split()[0]] = len(s)
    return out


def world_rows(w, o2_name="o2.integv.m2.tsv"):
    """One world's rows: every catalogued read with its class, margin, primary MAPQ and the four assigners' calls."""
    import pysam
    S = _load("summary")
    O = _load("o2graph")
    truth = {r["read"]: r for r in csv.DictReader(open(f"{w}/reads.truth.E0.tsv"), delimiter="\t")}
    scope = {q for q, t in truth.items() if t["kind"] in ("present", "silent")}
    copies = [c for c in csv.DictReader(open(f"{w}/copies.tsv"), delimiter="\t") if c["in_reference"] == "1"]
    fam_of = {c["copy"]: c["family"] for c in copies}
    family = collections.defaultdict(set)
    for c in copies:
        family[c["family"]].add(c["copy"])
    idx = locus_index(copies)
    lens = read_lengths(f"{w}/reads.E0.fq")
    alns, prim = collections.defaultdict(list), {}
    with pysam.AlignmentFile(f"{w}/E0.bam") as bam:
        for rec in bam.fetch(until_eof=True):
            q = rec.query_name
            if q not in scope or rec.is_unmapped or rec.is_supplementary or not rec.cigartuples:
                continue
            cp = locus_key(idx, rec.reference_name, rec.reference_start, rec.reference_end)
            aligned = sum(n for op, n in rec.cigartuples if op in (0, 1, 7, 8))
            alns[q].append((cp, rec.get_tag("NM") if rec.has_tag("NM") else 0, aligned))
            if not rec.is_secondary:
                prim[q] = (cp, rec.mapping_quality)
    iso = json.load(open(f"{w}/isoforms.json"))
    cat_copy = {(f, c["idx"]): c["copy"] for f, cps in iso.items() for c in cps}
    cur = collections.defaultdict(dict)
    for r in csv.DictReader(open(f"{w}/cur/A.assignments.tsv"), delimiter="\t"):
        if r["status"] == "assigned" and r["catalog_copy_idx"].isdigit():
            cur[r["read_name"]][r["family_id"]] = cat_copy.get((r["family_id"], int(r["catalog_copy_idx"])))
    sj = json.load(open(f"{w}/graph/seeded.json"))
    seeded_copy = dict(cat_copy)
    for f, v in sj.items():
        for i, cp in O.map_new_copies([c for c in v["copies"] if c["origin"] == "new"], copies).items():
            seeded_copy[(f, i)] = cp
    integ = S.integrated_calls(csv.DictReader(open(f"{w}/graph/{o2_name}"), delimiter="\t"), seeded_copy)

    def fam_call(calls, f):
        if not calls:
            return None
        if f in calls:
            return calls[f] or "no-truth-copy"
        return "wrong-family"

    rows = []
    for q in sorted(scope):
        t = truth[q]
        cls, m = decidability(copy_costs(alns.get(q, []), lens.get(q, 0)), t["copy"], family[t["family"]])
        pc, mq = prim.get(q, (None, 0))
        mm2 = None if q not in prim else (pc if isinstance(pc, str) and fam_of.get(pc) == t["family"] else "wrong-locus")
        rows.append(dict(read=q, family=t["family"], copy=t["copy"], cls=cls, margin=m, mapq0=mq == 0,
                         calls=dict(mm2=mm2, mm2_q1=mm2 if mq >= 1 else None,
                                    current=fam_call(cur.get(q), t["family"]),
                                    integrated=fam_call(integ.get(q), t["family"]))))
    return rows


def scope_rows(rows, tied):
    """Adds the scoped systems' calls to rows in place; tied = {read: the copies its best AS ties on}."""
    for r in rows:
        c = r["calls"]
        for m in ("current", "integrated"):
            c[f"{m}_tied"] = (c[m] if c[m] in tied[r["read"]] else None) if r["read"] in tied else c["mm2"]
    return rows


def summary_lines(card):
    """Per method and class: n, correct / wrong / abstain over all reads and over MAPQ-0 reads."""
    out = []
    for m in card:
        for c in CLASSES:
            a, z = card[m][c]["all"], card[m][c]["mapq0"]
            n = sum(a.values())
            if n:
                out.append(f"  {m:10s} {c:14s} n {n:6d} | correct {a['correct']:6d} wrong {a['wrong']:5d} abstain "
                           f"{a['abstain']:6d} | MAPQ 0: n {sum(z.values()):6d} correct {z['correct']:6d} wrong {z['wrong']:5d}")
    return out


def cmd_score(out, worlds, o2_name, tied_p=None):
    res, pooled = {}, []
    tied = ({r["read"]: set(r["copies"].split(",")) for r in csv.DictReader(open(tied_p), delimiter="\t")}
            if tied_p else None)
    methods = METHODS + (("current_tied", "integrated_tied") if tied is not None else ())
    for w in worlds:
        rows = world_rows(w, o2_name)
        if tied is not None:
            scope_rows(rows, tied)
        pooled += rows
        name = os.path.basename(w.rstrip("/"))
        res[name] = dict(card=scorecard(rows, methods),
                         classes=dict(collections.Counter(r["cls"] for r in rows)),
                         by_family={f: dict(collections.Counter(r["cls"] for r in rows if r["family"] == f))
                                    for f in sorted({r["family"] for r in rows})})
        print(f"{name}: " + ", ".join(f"{k} {v}" for k, v in sorted(res[name]["classes"].items())), flush=True)
    res["pooled"] = dict(card=scorecard(pooled, methods), classes=dict(collections.Counter(r["cls"] for r in pooled)))
    print("pooled:")
    print("\n".join(summary_lines(res["pooled"]["card"])))
    json.dump(res, open(out, "w"), indent=1)
    return 0


def score_args(argv):
    """(OUT, worlds, O2 table name, tied table or None) of `score`'s arguments."""
    opt = {k: argv[argv.index(k) + 1] for k in ("--o2", "--tied") if k in argv}
    worlds = [x for i, x in enumerate(argv[1:], 1) if x not in opt and argv[i - 1] not in ("--o2", "--tied")]
    return argv[0], worlds, opt.get("--o2", "o2.integv.m2.tsv"), opt.get("--tied")


def main(argv):
    if len(argv) >= 3 and argv[0] == "score":
        return cmd_score(*score_args(argv[1:]))
    print(__doc__, file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
