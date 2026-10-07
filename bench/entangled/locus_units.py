#!/usr/bin/env python3
"""Rule 1 (primary-first representative) and Rule 2 (structural separators) of docs/PREREG_locus_units_2026-10-06.md, as a post-processor of a families-input GTF.

    locus_units.py --base IN.families.gtf --out OUT.families.gtf [--primary-gtf P.families.gtf] [--rule1] [--rule2 all|primary] [--side OUT.separators.tsv]

Rule 1 needs --primary-gtf (the primaries-only assembly of the same reads): a transcript whose chain (chrom, strand, ordered intron list) occurs there gets reads += PRIMARY_BONUS, so the
shipped representative key (reads, span, -index) prefers primary-supported chains; nothing else uses `reads` in the from-GTF families path (src/family.rs fam_from_gtf).
Rule 2: per locus (gene_id) the incidence graph has the multi-exon transcripts and the junctions carried by >= 2 of them; a transcript is a SEPARATOR (a bridge) iff it is an articulation point and its
chain runs from a junction of one of the groups its removal leaves, through a junction no group carries, to a junction of another group; it is removed from the output (listed in --side) and the other
transcripts are regrouped by the components of the graph without the separators. With `--rule2 primary` the graph holds the primary-supported transcripts only; the other
multi-exon transcripts of the locus attach to the component sharing most junctions with them (ties: the earliest), those sharing none form their own groups; single-exon transcripts attach to the
component they overlap most on their strand (else their own group). The component holding the best representative (reads, span, -index) keeps the gene_id, the others are `<gene_id>.s<k>`, k = 2..,
orphans `<gene_id>.t<m>` and single-exon `<gene_id>.u<m>`. Loci without a separator are unchanged byte for byte.
"""
import argparse
import collections
import re

import networkx as nx

PRIMARY_BONUS = 1_000_000
ATTR = re.compile(r'(\w+) "([^"]*)"')


class Rec:
    __slots__ = ("fields", "feature", "chrom", "strand", "start", "end", "attrs", "tid", "gid", "line")

    def __init__(self, fields, line):
        self.fields = fields
        self.chrom, self.feature, self.strand = fields[0], fields[2], fields[6]
        self.start, self.end = int(fields[3]), int(fields[4])
        self.attrs = dict(ATTR.findall(fields[8]))
        self.tid, self.gid = self.attrs.get("transcript_id", ""), self.attrs.get("gene_id", "")
        self.line = line

    def set_attr(self, key, value):
        pat = re.compile(r'(\b' + key + r' ")([^"]*)(")')
        self.fields[8] = pat.sub(lambda m: m.group(1) + value + m.group(3), self.fields[8], count=1)
        self.attrs[key] = value
        if key == "gene_id":
            self.gid = value

    def text(self):
        return "\t".join(self.fields)


def parse_gtf(src):
    text = src if "\n" in src else open(src).read()
    recs = []
    for i, ln in enumerate(text.splitlines()):
        if not ln or ln.startswith("#"):
            continue
        f = ln.split("\t")
        if len(f) >= 9 and f[2] in ("transcript", "exon"):
            recs.append(Rec(f, i))
    return recs


def write_gtf(recs, path):
    with open(path, "w") as fh:
        for r in recs:
            fh.write(r.text() + "\n")


def transcripts(recs):
    """tid -> dict(chrom, strand, gid, exons (0-based half-open, sorted), reads, line, rec)."""
    tx = {}
    for r in recs:
        if r.feature == "transcript":
            tx[r.tid] = dict(chrom=r.chrom, strand=r.strand, gid=r.gid, exons=[], reads=int(r.attrs.get("reads", "0") or 0), line=r.line, rec=r)
    for r in recs:
        if r.feature == "exon" and r.tid in tx:
            tx[r.tid]["exons"].append((r.start - 1, r.end))
    for d in tx.values():
        d["exons"].sort()
        d["chain"] = tuple((d["exons"][i][1], d["exons"][i + 1][0]) for i in range(len(d["exons"]) - 1))
        d["span"] = d["exons"][-1][1] - d["exons"][0][0] if d["exons"] else 0
    return tx


def _is_bridge(chain, group_junctions):
    """True iff the chain runs, in order, from a junction carried by one group, through a junction carried by no group, to a junction carried by another group."""
    lab = []
    for j in chain:
        hit = [i for i, s in enumerate(group_junctions) if j in s]
        lab.append(hit[0] if len(hit) == 1 else -1)
    for a in range(len(lab)):
        if lab[a] < 0:
            continue
        for p in range(a + 1, len(lab)):
            if lab[p] != -1:
                continue
            for b in range(p + 1, len(lab)):
                if lab[b] >= 0 and lab[b] != lab[a]:
                    return True
    return False


def separators_components(chains):
    """chains: tid -> tuple of junctions in chain order. A transcript t is a SEPARATOR (a bridge) iff it is an articulation point of the incidence graph (transcripts and the
    junctions carried by >= 2 transcripts) and its chain runs from a junction of one of the groups its removal leaves, through a junction that no group carries (private to t),
    to a junction of another group. Returns (separators, comp), comp = the component index of every other transcript in the graph without the separators."""
    count = collections.Counter(j for ch in chains.values() for j in set(ch))
    g = nx.Graph()
    for t, ch in chains.items():
        g.add_node(("T", t))
        for j in set(ch):
            if count[j] >= 2:
                g.add_edge(("T", t), ("J", j))
    sep = set()
    for v in nx.articulation_points(g):
        if v[0] != "T":
            continue
        t = v[1]
        h = g.copy()
        h.remove_node(v)
        groups = [[w[1] for w in cc if w[0] == "T"] for cc in nx.connected_components(h)]
        gj = [{j for u in grp for j in chains[u]} for grp in groups if grp]
        if _is_bridge(chains[t], gj):
            sep.add(t)
    h = g.copy()
    h.remove_nodes_from([("T", t) for t in sep])
    comp, k = {}, 0
    for c in sorted((sorted(w[1] for w in cc if w[0] == "T") for cc in nx.connected_components(h)), key=lambda ts: min(ts) if ts else ""):
        if not c:
            continue
        for t in c:
            comp[t] = k
        k += 1
    return sep, comp


def _overlap(a, b):
    i = j = tot = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi > lo:
            tot += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def _union(exons):
    out = []
    for s, e in sorted(exons):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def primary_chains(recs):
    tx = transcripts(recs)
    return {(d["chrom"], d["strand"], d["chain"]) for d in tx.values() if d["chain"]}


def primary_singles(recs):
    tx = transcripts(recs)
    return {(d["chrom"], d["strand"], tuple(d["exons"])) for d in tx.values() if len(d["exons"]) == 1}


def apply_rule1(recs, primary, singles=frozenset()):
    tx = transcripts(recs)
    for t, d in tx.items():
        key = (d["chrom"], d["strand"], d["chain"])
        marked = key in primary if d["chain"] else (d["chrom"], d["strand"], tuple(d["exons"])) in singles
        if marked:
            d["rec"].set_attr("reads", str(d["reads"] + PRIMARY_BONUS))
    return recs


def apply_rule2(recs, pool=None):
    """pool: None (all transcripts build the graph) or a set of transcript ids. Returns (records, side rows)."""
    tx = transcripts(recs)
    by_locus = collections.defaultdict(list)
    for t, d in tx.items():
        by_locus[d["gid"]].append(t)
    new_gid, drop, side = {}, set(), []
    best = lambda ts: max(ts, key=lambda t: (tx[t]["reads"], tx[t]["span"], -tx[t]["line"]))
    for gid, ts in by_locus.items():
        graph_ts = [t for t in ts if tx[t]["chain"] and (pool is None or t in pool)]
        if len(graph_ts) < 2:
            continue
        sep, comp = separators_components({t: tx[t]["chain"] for t in graph_ts})
        if not sep:
            continue
        groups = collections.defaultdict(list)
        for t, k in comp.items():
            groups[k].append(t)
        rest = [t for t in ts if t not in sep and t not in comp]
        for t in sorted(sep, key=lambda t: tx[t]["line"]):
            drop.add(t)
            side.append(dict(tid=t, gid=gid, groups=len(groups), reads=tx[t]["reads"]))
        multi_rest = [t for t in rest if tx[t]["chain"]]
        single_rest = [t for t in rest if not tx[t]["chain"]]
        # multi-exon transcripts outside the graph: the group sharing most junctions, else groups of their own
        orphan_pool = []
        for t in multi_rest:
            js = set(tx[t]["chain"])
            score = sorted(((len(js & {j for u in groups[k] for j in tx[u]["chain"]}), -k) for k in groups), reverse=True)
            if score and score[0][0] > 0:
                groups[-score[0][1]].append(t)
            else:
                orphan_pool.append(t)
        orphan_groups = []
        if orphan_pool:
            _, ocomp = separators_components({t: tx[t]["chain"] for t in orphan_pool})
            og = collections.defaultdict(list)
            for t in orphan_pool:
                og[ocomp.get(t, -1 - len(og))].append(t)
            orphan_groups = list(og.values())
        # naming: the group holding the best representative keeps gid
        glist = [groups[k] for k in sorted(groups)]
        glist.sort(key=lambda ts: min(tx[t]["line"] for t in ts))
        keeper = max(range(len(glist)), key=lambda i: (tx[best(glist[i])]["reads"], tx[best(glist[i])]["span"], -tx[best(glist[i])]["line"]))
        n = 2
        for i, g in enumerate(glist):
            name = gid if i == keeper else f"{gid}.s{n}"
            if i != keeper:
                n += 1
            for t in g:
                new_gid[t] = name
        for m, g in enumerate(orphan_groups, 1):
            for t in g:
                new_gid[t] = f"{gid}.t{m}"
        # single-exon: the group they overlap most on their strand, else their own
        un = 1
        gunion = {}
        for t, name in list(new_gid.items()):
            if t in sep:
                continue
        members = collections.defaultdict(list)
        for t, name in new_gid.items():
            members[name].append(t)
        for name, mem in members.items():
            gunion[name] = {}
            for t in mem:
                gunion[name].setdefault((tx[t]["chrom"], tx[t]["strand"]), []).extend(tx[t]["exons"])
        for t in sorted(single_rest, key=lambda t: tx[t]["line"]):
            key = (tx[t]["chrom"], tx[t]["strand"])
            sc = sorted(((_overlap(_union(gunion[name].get(key, [])), tx[t]["exons"]), name) for name in gunion), key=lambda x: (-x[0], x[1] != gid, x[1]))
            if sc and sc[0][0] > 0:
                new_gid[t] = sc[0][1]
            else:
                new_gid[t] = f"{gid}.u{un}"; un += 1
    out = []
    for r in recs:
        if r.tid in drop:
            continue
        if r.tid in new_gid and new_gid[r.tid] != r.gid:
            r.set_attr("gene_id", new_gid[r.tid])
        out.append(r)
    return out, side


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--primary-gtf")
    ap.add_argument("--rule1", action="store_true")
    ap.add_argument("--rule2", choices=["all", "primary"])
    ap.add_argument("--side")
    a = ap.parse_args()
    recs = parse_gtf(a.base)
    pchains = psingles = None
    if a.rule1 or a.rule2 == "primary":
        if not a.primary_gtf:
            raise SystemExit("--primary-gtf is required for --rule1 and --rule2 primary")
        prec = parse_gtf(a.primary_gtf)
        pchains, psingles = primary_chains(prec), primary_singles(prec)
    pool = None
    if a.rule2 == "primary":
        tx0 = transcripts(recs)
        pool = {t for t, d in tx0.items() if d["chain"] and (d["chrom"], d["strand"], d["chain"]) in pchains}
    if a.rule1:
        recs = apply_rule1(recs, pchains, psingles)
    side = []
    if a.rule2:
        recs, side = apply_rule2(recs, pool)
    write_gtf(recs, a.out)
    if a.side:
        with open(a.side, "w") as fh:
            fh.write("tid\tgid\tgroups\treads\n")
            for s in side:
                fh.write(f"{s['tid']}\t{s['gid']}\t{s['groups']}\t{s['reads']}\n")
    print(f"[locus_units] {len(transcripts(parse_gtf(a.base)))} -> {len(transcripts(recs))} transcripts, {len(side)} separators removed")


if __name__ == "__main__":
    main()
