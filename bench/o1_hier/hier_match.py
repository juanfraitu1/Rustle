#!/usr/bin/env python3
"""Hierarchy-aware O1 scorer (docs/superpowers/specs/2026-10-13-close-o1-design.md, section A).

    hier_match.py MEMBERS_TSV TRUTH_TSV OUT.json [--loci-gff LOCI_GFF3] [--fold LOCI_TSV] [--map exon|key]
                  [--null 1000] [--seed 7]        (--loci-gff is required by --map exon, the default)

Nodes: every family, expansion, superfamily and SD complex of `<out>.hierarchy.members.tsv`, each as its set of member
loci (`chrom:start-end`, GFF 1-based). Duplicons, SD blocks and external families are not nodes.
Truth: the `soto_truth.py` table (`label gene_id chrom start end strand exons`, GFF 1-based) or the 4-column
`chrom start end label` of `--hierarchy-external` (no strand; the span is its one exon); `Unassigned_*` labels are not
families. `--map exon` (default, de novo catalogs): a truth gene joins a member when one of its exons shares >= 1 bp,
on the same strand, with an exon of the member's locus or of a locus folded into it (`--fold`, the run's
`<out>.loci.tsv`; exons from `--loci-gff`, the catalog's `.fam.loci.gff3`). `--map key`: identical coordinates only,
directly or through the fold map (the 10-09 F_O universe, where the members ARE the annotation records).
A truth family is EVALUABLE with >= 2 members. Matching: one-to-one Hungarian (scipy linear_sum_assignment) on
-Jaccard; rows in label order, columns in (level rank, node id) order, rank family < expansion < superfamily <
sd_complex: the named tie policy. `exact_any` asks only whether an equal node exists (tie-free, not one-to-one).
Nulls keep every node's size and the nesting: `global` relabels all members by one random permutation per round;
`within_family` permutes members inside their own family, so families, superfamilies and SD complexes are unchanged
and only the expansions move.
"""
import argparse
import bisect
import collections
import json
import random
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix

LEVELS = ("family", "expansion", "superfamily", "sd_complex")
RANK = {level: i for i, level in enumerate(LEVELS)}
Match = collections.namedtuple("Match", "label level node_id jaccard recall precision inter t_size n_size")
Record = collections.namedtuple("Record", "chrom start end strand exons")


def parse_key(k):
    c, r = k.rsplit(":", 1)
    s, e = r.split("-")
    return c, int(s), int(e)


def load_nodes(members_tsv):
    """-> ({(level, node_id): frozenset(member)}, {member: family})."""
    nodes, family_of = collections.defaultdict(set), {}
    with open(members_tsv) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {k: hdr.index(k) for k in ("sd_complex", "superfamily", "family", "member", "expansions")}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            m = f[ix["member"]]
            family_of[m] = f[ix["family"]]
            for level in ("sd_complex", "superfamily", "family"):
                if f[ix[level]] not in ("", "-"):
                    nodes[(level, f[ix[level]])].add(m)
            for x in f[ix["expansions"]].split(","):
                if x not in ("", "-"):
                    nodes[("expansion", x)].add(m)
    return {k: frozenset(v) for k, v in nodes.items()}, family_of


def load_truth(path):
    """-> {label: [Record]} (GFF 1-based closed); `Unassigned_*` labels dropped."""
    out = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "label":
                continue
            if len(f) >= 7:
                label, chrom, s, e, strand = f[0], f[2], int(f[3]), int(f[4]), f[5]
                exons = [tuple(int(x) for x in b.split("-")) for b in f[6].split(",")]
            elif len(f) >= 4:
                label, chrom, s, e, strand = f[3], f[0], int(f[1]), int(f[2]), "."
                exons = [(s, e)]
            else:
                continue
            if not label.startswith("Unassigned_"):
                out[label].append(Record(chrom, s, e, strand, exons))
    return dict(out)


def load_fold(path):
    """`<out>.loci.tsv` (annotation -> representative) -> {record key: member key}."""
    with open(path) as fh:
        fh.readline()
        return dict(line.rstrip("\n").split("\t")[:2] for line in fh if line.strip())


def load_loci_gff(path):
    """Every locus of a catalog (`.fam.loci.gff3`: gene lines + exon lines by Parent) -> {key: (chrom, strands, exons)};
    a locus without exon lines is its span (as `exonless_span`); loci sharing a key pool their exons and strands."""
    genes, exons = {}, collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] == "gene":
                genes[at["ID"]] = (f"{f[0]}:{f[3]}-{f[4]}", f[0], f[6], int(f[3]), int(f[4]))
            elif f[2] == "exon":
                exons[at["Parent"]].append((int(f[3]), int(f[4])))
    out = {}
    for gid, (key, chrom, strand, s, e) in genes.items():
        _, strands, ex = out.get(key, (chrom, frozenset(), []))
        out[key] = (chrom, strands | {strand}, ex + (exons.get(gid) or [(s, e)]))
    return out


class SpanIndex:
    """Closed intervals per chromosome with a payload; overlap queries by binary search."""

    def __init__(self, items):
        by = collections.defaultdict(list)
        for c, s, e, p in items:
            by[c].append((s, e, p))
        self.by, self.starts, self.maxlen = by, {}, {}
        for c, v in by.items():
            v.sort(key=lambda x: (x[0], x[1]))
            self.starts[c] = [x[0] for x in v]
            self.maxlen[c] = max(x[1] - x[0] for x in v)

    def query(self, c, s, e):
        v = self.by.get(c)
        if not v:
            return []
        lo = bisect.bisect_left(self.starts[c], s - self.maxlen[c])
        hi = bisect.bisect_right(self.starts[c], e)
        return [p for a, b, p in v[lo:hi] if a <= e and s <= b]


def _strand_ok(a, b):
    return "." in a or "." in b or bool(a & b)


def truth_sets(truth, members, fold=None, loci=None, how="exon"):
    """truth {label: [Record]}; members = member keys; fold {record key: member key}; loci = load_loci_gff (exon mode).
    -> (sets {label: frozenset(member)}, per_record {label: [frozenset(member)]}, coverage)."""
    mset = set(members)
    fold = {r: m for r, m in (fold or {}).items() if m in mset}
    cov, all_idx = collections.Counter(), None
    if how == "exon":
        items = []
        for key, m in [(m, m) for m in sorted(mset)] + sorted(fold.items()):
            if key not in loci:
                cov["member_loci_missing_from_gff"] += 1
                c, s, e = parse_key(key)
                items.append((c, s, e, (frozenset("."), m)))
                continue
            chrom, strands, ex = loci[key]
            items.extend((chrom, s, e, (strands, m)) for s, e in ex)
        idx = SpanIndex(items)
        all_idx = SpanIndex([(chrom, s, e, (strands, 1)) for chrom, strands, ex in loci.values() for s, e in ex])

        def hit(r, index=idx):
            return {p for s, e in r.exons for st, p in index.query(r.chrom, s, e) if _strand_ok(st, {r.strand})}
    else:
        def hit(r, index=None):
            k = f"{r.chrom}:{r.start}-{r.end}"
            return {k} if k in mset else ({fold[k]} if k in fold else set())
    sets, per_record, without = {}, {}, []
    for label in sorted(truth):
        recs = [frozenset(hit(r)) for r in truth[label]]
        per_record[label] = recs
        sets[label] = frozenset().union(*recs)
        cov["records"] += len(recs)
        for r, h in zip(truth[label], recs):
            if h:
                cov["records_on_members"] += 1
                cov["records_on_several_members"] += len(h) > 1
            elif all_idx is not None and hit(r, all_idx):
                cov["records_on_unclustered_locus"] += 1
            else:
                cov["records_without_locus" if all_idx is not None else "records_without_member"] += 1
        n = len(sets[label])
        cov["labels"] += 1
        cov["labels_evaluable" if n >= 2 else ("labels_one_member" if n == 1 else "labels_no_member")] += 1
        if n == 0:
            without.append(label)
    shared = collections.Counter(m for s in sets.values() for m in s)
    out = dict(sorted(cov.items()))
    out["members_in_several_truth_families"] = sum(1 for v in shared.values() if v > 1)
    out["labels_without_member"] = without
    return sets, per_record, out


def _incidence(sets, index):
    rows = [index[m] for s in sets for m in s]
    cols = [j for j, s in enumerate(sets) for _ in s]
    return csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(index), len(sets)))


class Scorer:
    """One-to-one Hungarian between truth sets and nodes on -Jaccard, with the named tie policy (rows in label
    order, columns in (level rank, node id) order)."""

    def __init__(self, nodes, members):
        self.keys = sorted(nodes, key=lambda k: (RANK[k[0]], k[1]))
        self.index = {m: i for i, m in enumerate(sorted(members))}
        self.sizes = np.array([len(nodes[k]) for k in self.keys], float)
        self.N = _incidence([nodes[k] for k in self.keys], self.index)

    def run(self, truth):
        labels = sorted(truth)
        out = {}
        if labels and self.keys:
            T = _incidence([truth[x] for x in labels], self.index)
            inter = (T.T @ self.N).tocsc()
            used = np.unique(inter.nonzero()[1])
            if used.size:
                I = inter[:, used].toarray()
                ts = np.array([len(truth[x]) for x in labels], float)
                J = I / (ts[:, None] + self.sizes[used][None, :] - I)
                for i, j in zip(*linear_sum_assignment(-J)):
                    if I[i, j] > 0:
                        level, nid = self.keys[used[j]]
                        n = self.sizes[used[j]]
                        out[labels[i]] = Match(labels[i], level, nid, float(J[i, j]), float(I[i, j] / ts[i]),
                                               float(I[i, j] / n), int(I[i, j]), int(ts[i]), int(n))
        return [out.get(x) or Match(x, None, None, 0.0, 0.0, 0.0, 0, len(truth[x]), 0) for x in labels]


def summary(matches, nodes=None, per_record=None):
    n = len(matches)
    hit = [m for m in matches if m.level]
    exact = sum(m.jaccard == 1.0 for m in matches)
    out = {
        "n_truth": n,
        "n_matched": len(hit),
        "exact": exact,
        "exact_rate": exact / n if n else None,
        "mean_jaccard": sum(m.jaccard for m in matches) / n if n else None,
        "macro_recall": sum(m.recall for m in matches) / n if n else None,
        "macro_precision": sum(m.precision for m in hit) / len(hit) if hit else None,
        "micro_recall": sum(m.inter for m in matches) / sum(m.t_size for m in matches) if n else None,
        "micro_precision": sum(m.inter for m in hit) / sum(m.n_size for m in hit) if hit else None,
        "by_level": dict(sorted(collections.Counter(m.level for m in hit).items())),
        "exact_by_level": dict(sorted(collections.Counter(m.level for m in matches if m.jaccard == 1.0).items())),
    }
    if nodes is not None and per_record is not None and n:
        out["macro_gene_recall"] = sum(gene_recall(m, nodes, per_record) for m in matches) / n
    return out


def gene_recall(m, nodes, per_record):
    """Share of the truth family's RECORDS (with or without a member) that have a member in the matched node."""
    recs = per_record[m.label]
    if not m.level or not recs:
        return 0.0
    node = nodes[(m.level, m.node_id)]
    return sum(1 for r in recs if r & node) / len(recs)


def exact_any(truth, nodes):
    """Tie-free: per level, how many truth sets equal SOME node of that level (no one-to-one)."""
    levels_of = collections.defaultdict(set)
    for (level, _), s in nodes.items():
        levels_of[s].add(level)
    c = collections.Counter()
    for s in truth.values():
        ls = levels_of.get(s, set())
        for level in ls:
            c[level] += 1
        c["any"] += bool(ls)
    return dict(sorted(c.items()))


def permutation(members, family_of, kind, rng):
    """A bijection member -> member: `global` over all members, `within_family` inside each family."""
    groups = collections.defaultdict(list)
    for m in sorted(members):
        groups["*" if kind == "global" else family_of[m]].append(m)
    out = {}
    for g in sorted(groups):
        dst = groups[g][:]
        rng.shuffle(dst)
        out.update(zip(groups[g], dst))
    return out


def relabel(sets, mapping):
    return {k: frozenset(mapping[m] for m in v) for k, v in sets.items()}


def null(scorer, truth, members, family_of, kind, n, seed):
    """Relabelling the truth by sigma^-1 scores exactly as relabelling every node by sigma, so the nodes stay fixed.
    -> (exact counts, mean Jaccards), one per round."""
    rng = random.Random(f"{seed}:{kind}")
    ex, mj = [], []
    for _ in range(n):
        ms = scorer.run(relabel(truth, permutation(members, family_of, kind, rng)))
        ex.append(sum(m.jaccard == 1.0 for m in ms))
        mj.append(sum(m.jaccard for m in ms) / len(ms) if ms else 0.0)
    return ex, mj


def p_value(values, observed):
    return (1 + sum(v >= observed for v in values)) / (1 + len(values))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("members")
    ap.add_argument("truth")
    ap.add_argument("out")
    ap.add_argument("--fold")
    ap.add_argument("--loci-gff")
    ap.add_argument("--map", choices=("exon", "key"), default="exon")
    ap.add_argument("--null", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args(argv)
    nodes, family_of = load_nodes(a.members)
    members = sorted(family_of)
    if a.map == "exon" and not a.loci_gff:
        ap.error("--map exon needs --loci-gff (the catalog's .fam.loci.gff3)")
    sets, per_record, cov = truth_sets(load_truth(a.truth), members, load_fold(a.fold) if a.fold else None,
                                       load_loci_gff(a.loci_gff) if a.loci_gff else None, a.map)
    ev = {k: v for k, v in sets.items() if len(v) >= 2}
    hier_scorer = Scorer(nodes, members)
    fam_nodes = {k: v for k, v in nodes.items() if k[0] == "family"}
    hier, flat = hier_scorer.run(ev), Scorer(fam_nodes, members).run(ev)
    res = {"inputs": vars(a), "n_members": len(members),
           "n_nodes": dict(sorted(collections.Counter(k[0] for k in nodes).items())), "coverage": cov,
           "hierarchy": summary(hier, nodes, per_record), "flat": summary(flat, nodes, per_record),
           "exact_any": exact_any(ev, nodes), "null": {}}
    for kind in ("global", "within_family"):
        if a.null > 0:
            ex, mj = null(hier_scorer, ev, members, family_of, kind, a.null, a.seed)
            res["null"][kind] = {"rounds": a.null, "exact_mean": sum(ex) / len(ex), "exact_max": max(ex),
                                 "mean_jaccard_mean": sum(mj) / len(mj), "mean_jaccard_max": max(mj),
                                 "p_exact": p_value(ex, res["hierarchy"]["exact"]),
                                 "p_mean_jaccard": p_value(mj, res["hierarchy"]["mean_jaccard"])}
    fm = {m.label: m for m in flat}
    res["per_truth"] = [{"label": m.label, "n_records": len(per_record[m.label]), "n_members": m.t_size,
                         "level": m.level, "node": m.node_id, "jaccard": round(m.jaccard, 6),
                         "recall": round(m.recall, 6), "precision": round(m.precision, 6),
                         "gene_recall": round(gene_recall(m, nodes, per_record), 6),
                         "flat_node": fm[m.label].node_id, "flat_jaccard": round(fm[m.label].jaccard, 6)}
                        for m in hier]
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ("per_truth", "inputs")}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
