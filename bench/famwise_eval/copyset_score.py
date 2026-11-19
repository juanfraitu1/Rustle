#!/usr/bin/env python3
"""The measure of docs/GOAL_multicopy_families_2026-11-19.md: transcripts as claims about copies, scored against the
origin of their reads. A claim is a set of keys (contig, strand, intron chain): one transcript of a per-copy GTF, or the
transcripts of one family transcript (`fw_tx` attribute of the `famwise` GTF, `--fw NAME=FW.gtf`).

  copyset_score.py --truth READS.tsv --tie-split T.tsv --bam SEL.bam --out DIR
                   [--arm NAME=GTF ...] [--fw NAME=FW.gtf ...] [--min-reads 2]

Supporters of a claim: judged molecules (truth kinds one / several) with a record within 10 ms of their best (or the
primary) carrying exactly one of the claim's keys (origin_score.py's rule). Per supporter, the claim's copies are mapped
onto the read's truth loci; a claimed copy the read has no locus at counts as a copy the read did not come from.
  EQUAL   claimed copies = the read's origin set
  NARROW  claimed copies a proper subset (the claim states more resolution than the read has)
  WIDE    a proper superset (claims copies the read did not come from)
  OFF     disjoint;  MIXED otherwise
A claim's verdict is its supporters' strict plurality verdict (a tie for the top count = MIXED); no judged supporter =
UNJUDGED. Precision = EQUAL / judged claims; beside it, the molecule-weighted share = supporters of EQUAL claims / supporters
of judged claims (grouping copies into one claim shrinks the claim count; the weighted share does not move with it).
Recall units: a judged molecule whose records within 10 ms carry a chain at EVERY origin locus gives the unit = its keys
at its origin loci; units held by >= --min-reads molecules are true units. Recall per arm: structure (every key of the
unit is emitted by the arm) and resolution (one claim holds every key of the unit and claims exactly the unit's copies).
Regime of a unit or claim: `one` = a single copy, `set` = two or more.
Writes DIR/claims.tsv (arm x regime x verdict), DIR/recall.tsv, DIR/summary.tsv."""
import argparse, collections, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import records as R
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
from gtf_names import chain, opn

VERDICTS = ("EQUAL", "NARROW", "WIDE", "MIXED", "OFF")

def load_gtf(path, fw=False):
    """-> list of claims; a claim = dict(keys=frozenset((contig, strand, chain)), spans=[(contig, strand, s, e)]).
    fw: transcripts sharing a `fw_tx` value form one claim (transcripts without one stay alone)."""
    ex, at, grp = collections.defaultdict(list), {}, {}
    for l in opn(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8])
        if not t:
            continue
        t = t.group(1)
        if f[2] == "exon":
            ex[t].append((int(f[3]), int(f[4]))); at[t] = (f[0], f[6])
        elif f[2] == "transcript" and fw:
            m = re.search(r'fw_tx "([^"]+)"', f[8])
            if m:
                grp[t] = m.group(1)
    groups = collections.defaultdict(list)
    for t, e in ex.items():
        if len(e) < 2:
            continue
        e.sort()
        groups[grp.get(t, "tx:" + t)].append((at[t][0], at[t][1], chain(e), e[0][0] - 1, e[-1][1]))
    out = []
    for g, v in groups.items():
        out.append(dict(keys=frozenset((c, s, ch) for c, s, ch, _, _ in v), spans=[(c, s, a, b) for c, s, _, a, b in v]))
    return out

def overlaps(locus, span):
    c, s, e, st = locus[0], locus[1], locus[2], locus[3]
    return c == span[0] and st == span[1] and s < span[3] and span[2] < e

def claimed_local(claim, t):
    """the claim's copies in the read's locus indices; -1 = a claimed copy at none of the read's judged loci"""
    out = set()
    for sp in claim["spans"]:
        hit = [i for i, lo in enumerate(t.loci) if overlaps(lo, sp)]
        hit = [i for i in hit if R.locus_index(R.Rec(sp[0], sp[2], sp[3], sp[1], (), 0, 0, False), t) is not None] or [-1]
        out.add(hit[0])
    return out

def read_verdict(cl, origin):
    if cl == origin:
        return "EQUAL"
    if cl < origin:
        return "NARROW"
    if cl > origin:
        return "WIDE"
    return "OFF" if not (cl & origin) else "MIXED"

def plurality(vs):
    """the strict plurality verdict; a tie for the top count is MIXED"""
    c = collections.Counter(vs).most_common()
    return "MIXED" if len(c) > 1 and c[0][1] == c[1][1] else c[0][0]

class Containment:
    """read chain -> the arm keys that hold it as a contiguous run of introns (same contig and strand)"""
    def __init__(self, keys):
        self.by = collections.defaultdict(list)
        for k in keys:
            for i, j in enumerate(k[2]):
                self.by[(k[0], k[1], j)].append((k, i))
    def holders(self, c, st, ch):
        out = []
        for k, i in self.by.get((c, st, ch[0]), ()):
            if k[2][i:i + len(ch)] == ch:
                out.append(k)
        return out

def supporters(truth, recs, cont=None):
    """-> {key: [read]} over judged molecules (records within 10 ms of the best, or the primary); with a Containment
    index, a record supports every key holding its chain as a contiguous run, else only the key equal to its chain"""
    sup = collections.defaultdict(list)
    for q, rs in recs.items():
        t = truth.get(q)
        if t is None or t.kind not in R.JUDGED:
            continue
        for r in R.statement(rs, "ms", 10):
            if r.chain and R.locus_index(r, t) is not None:
                for k in (cont.holders(r.chrom, r.strand, r.chain) if cont else [(r.chrom, r.strand, r.chain)]):
                    sup[k].append(q)
    return sup

def units(truth, recs, min_reads):
    """-> {unit keys (frozenset): (n reads, regime, origin loci spans of one read)}"""
    n, info = collections.Counter(), {}
    for q, rs in recs.items():
        t = truth.get(q)
        if t is None or t.kind not in R.JUDGED:
            continue
        at = {}
        for r in R.statement(rs, "ms", 10):
            i = R.locus_index(r, t)
            if r.chain and i in t.origin and i not in at:
                at[i] = (r.chrom, r.strand, r.chain)
        if set(at) != set(t.origin):
            continue
        u = frozenset(at.values())
        if len(u) != len(t.origin):
            continue
        n[u] += 1
        info.setdefault(u, ("one" if len(t.origin) == 1 else "set", [t.loci[i] for i in sorted(t.origin)]))
    return {u: (k, *info[u]) for u, k in n.items() if k >= min_reads}

def same_copies(claim, loci):
    """the claim's spans and the unit's origin loci cover each other"""
    return all(any(overlaps(lo, sp) for lo in loci) for sp in claim["spans"]) and \
           all(any(overlaps(lo, sp) for sp in claim["spans"]) for lo in loci)

def n_places(spans):
    """overlapping spans on one contig and strand are one place"""
    out = []
    for c, st, s, e in sorted(spans):
        if out and out[-1][0] == c and out[-1][1] == st and s < out[-1][3]:
            out[-1][3] = max(out[-1][3], e)
        else:
            out.append([c, st, s, e])
    return len(out)

def score_arm(claims, sup, truth, dump=None):
    """-> Counter[(regime, verdict)]; dump: a list that receives (claim, regime, verdict, judged supporters)"""
    out = collections.Counter()
    for cl in claims:
        regime = "one" if n_places(cl["spans"]) == 1 else "set"
        reads = {q for k in cl["keys"] for q in sup.get(k, ())}
        vs = []
        for q in reads:
            t = truth[q]
            vs.append(read_verdict(claimed_local(cl, t), set(t.origin)))
        v = plurality(vs) if vs else "UNJUDGED"
        out[(regime, v)] += 1
        if dump is not None:
            dump.append((cl, regime, v, len(vs)))
        if vs:
            out[(regime, "w_judged")] += len(vs); out[(regime, "w_equal")] += len(vs) * (v == "EQUAL")
    return out

def recall_arm(claims, us, contained=False):
    keys = set().union(*(c["keys"] for c in claims)) if claims else set()
    cont = Containment(keys) if contained else None
    hold = (lambda k: set(cont.holders(*k))) if contained else (lambda k: {k} & keys)
    by_key = collections.defaultdict(list)
    for c in claims:
        for k in c["keys"]:
            by_key[k].append(c)
    out = collections.Counter()
    for u, (n, regime, loci) in us.items():
        out[(regime, "units")] += 1
        H = {k: hold(k) for k in u}
        if all(H.values()):
            out[(regime, "structure")] += 1
        cands = {id(c): c for h in H[next(iter(u))] for c in by_key.get(h, ())}
        if any(all(H[k] & c["keys"] for k in u) and same_copies(c, loci) for c in cands.values()):
            out[(regime, "resolution")] += 1
    return out

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True); ap.add_argument("--tie-split", required=True); ap.add_argument("--bam", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--arm", action="append", default=[]); ap.add_argument("--fw", action="append", default=[])
    ap.add_argument("--min-reads", type=int, default=2)
    ap.add_argument("--new-lost", help="arm name to compare every other arm with: claims by key set present in one and not the other, with verdicts (DIR/new_lost.tsv)")
    ap.add_argument("--match", choices=("exact", "contained"), default="exact",
                    help="exact: a read chain matches a transcript chain only when equal; contained: when it is a contiguous run of it")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    truth, recs = R.load_truth(a.truth, a.tie_split), R.load_records(a.bam)
    us = units(truth, recs, a.min_reads)
    arms = [(n, load_gtf(g)) for n, g in (x.split("=", 1) for x in a.arm)] + [(n, load_gtf(g, fw=True)) for n, g in (x.split("=", 1) for x in a.fw)]
    rows, rrows, summ = [], [], []
    sup = supporters(truth, recs) if a.match == "exact" else None
    dumps = {}
    for n, claims in arms:
        if a.match == "contained":
            sup = supporters(truth, recs, Containment(set().union(*(c["keys"] for c in claims)) if claims else set()))
        dumps[n] = []
        c, r = score_arm(claims, sup, truth, dumps[n]), recall_arm(claims, us, a.match == "contained")
        for reg in ("one", "set"):
            judged = sum(c[(reg, v)] for v in VERDICTS)
            rows.append((n, reg, sum(c[(reg, v)] for v in (*VERDICTS, "UNJUDGED")), judged, *(c[(reg, v)] for v in VERDICTS), c[(reg, "UNJUDGED")]))
            rrows.append((n, reg, r[(reg, "units")], r[(reg, "structure")], r[(reg, "resolution")]))
            lo, hi = R.wilson(c[(reg, "EQUAL")], judged)
            ru = r[(reg, "units")]
            wj = c[(reg, "w_judged")]
            summ.append(f"{n}\t{reg}\tprecision {c[(reg, 'EQUAL')]}/{judged} = {c[(reg, 'EQUAL')] / judged if judged else 0:.3f} [{lo:.3f}, {hi:.3f}]\t"
                        f"weighted {c[(reg, 'w_equal')] / wj if wj else 0:.3f}\t"
                        f"recall structure {r[(reg, 'structure')]}/{ru}, resolution {r[(reg, 'resolution')]}/{ru}")
    with open(f"{a.out}/claims.tsv", "w") as fh:
        fh.write("arm\tregime\tclaims\tjudged\t" + "\t".join(VERDICTS) + "\tUNJUDGED\n")
        fh.writelines("\t".join(map(str, x)) + "\n" for x in rows)
    with open(f"{a.out}/recall.tsv", "w") as fh:
        fh.write("arm\tregime\tunits\tstructure\tresolution\n")
        fh.writelines("\t".join(map(str, x)) + "\n" for x in rrows)
    if a.new_lost and a.new_lost in dumps:
        base = {d[0]["keys"]: d for d in dumps[a.new_lost]}
        with open(f"{a.out}/new_lost.tsv", "w") as fh:
            fh.write(f"arm\tside\tverdict\tclaims\tsupporters\n")
            for n, ds in dumps.items():
                if n == a.new_lost:
                    continue
                mine = {d[0]["keys"]: d for d in ds}
                for side, src, other in (("new", mine, base), ("lost", base, mine)):
                    cnt = collections.Counter(); sp = collections.Counter()
                    for k, d in src.items():
                        if k not in other:
                            cnt[d[2]] += 1; sp[d[2]] += d[3]
                    for v in (*VERDICTS, "UNJUDGED"):
                        fh.write(f"{n}\t{side}\t{v}\t{cnt[v]}\t{sp[v]}\n")
    with open(f"{a.out}/summary.tsv", "w") as fh:
        fh.write(f"# match {a.match}; true units (>= {a.min_reads} molecules): one {sum(v[1] == 'one' for v in us.values())}, set {sum(v[1] == 'set' for v in us.values())}\n")
        fh.write("\n".join(summ) + "\n")
    print(open(f"{a.out}/summary.tsv").read(), end="")

if __name__ == "__main__":
    main()
