#!/usr/bin/env python3
"""Is a POSET a better family object than a partition of an undirected graph?
Per `docs/PREREG_containment_poset_2026-09-22.md` (§6y0, md5 `4c0be9cc`).

r1020 measured that containment between loci is a strict partial order (antisymmetry 0 violations,
transitivity 92.9%, longest chain 38, AUC 0.943 for same-family and NOT a size artefact). r1021 proposed
that §6s9's inexpressible fusion, §6u8's dropped cover and §6t7's asymmetry are one limitation: the output
is a partition of an undirected graph while the data is a poset. **This tests whether that scores.**

⚠ THE EDGE SET IS HELD FIXED. Both arms consume the shipped `--dump-graph` output, which already passed
every shipped conjunct; direction is then read off the PAF. No edge is added or removed by either arm, so a
difference is attributable to the structure alone.
⚠ The comparator is `mcl_port` MCL ON THE SAME GRAPH, never the shipped Rust F (register 917).
⚠ PRIMARY metric is PAIRWISE, because bipartite matching is not well defined for a cover and penalises one
by construction (§6u8's own warning).

Usage:
  containment_poset.py --graph X.graph.tsv --paf X.paf --gff chrN.genes.gff --truth T.tsv --chrom chrN
"""
import argparse
import bisect
import collections
import itertools
import re
import sys

sys.path.insert(0, 'bench')
import mcl_port


def merged_len(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return sum(e - s for s, e in out)


def read_graph(path):
    """the shipped pre-MCL graph: `u<TAB>v<TAB>w`, plus a self-row per node."""
    adj = collections.defaultdict(dict)
    nodes = set()
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) == 3 and f[0] != f[1]:
            a, b, w = f[0], f[1], float(f[2])
            adj[a][b] = max(adj[a].get(b, 0.0), w)
            adj[b][a] = max(adj[b].get(a, 0.0), w)
            nodes.add(a); nodes.add(b)
        elif f and f[0]:
            nodes.add(f[0])
    return adj, nodes


def directed(paf, keep, C):
    """x -> y ('x is contained in y') for edges already in the shipped graph."""
    qlen = {}
    acc = collections.defaultdict(lambda: [[], []])
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 11:
            continue
        q, t = f[0], f[5]
        qlen[q] = int(f[1]); qlen[t] = int(f[6])
        if q == t:
            continue
        key = (q, t) if q <= t else (t, q)
        if key not in keep:
            continue
        a = (int(f[2]), int(f[3])); b = (int(f[7]), int(f[8]))
        e = acc[key]
        if key[0] == q:
            e[0].append(a); e[1].append(b)
        else:
            e[0].append(b); e[1].append(a)
    below = collections.defaultdict(set)
    for (x, y), (ix, iy) in acc.items():
        if not ix:
            continue
        cx = merged_len(ix) / max(qlen[x], 1)
        cy = merged_len(iy) / max(qlen[y], 1)
        if cx >= C and cx > cy:
            below[x].add(y)
        elif cy >= C and cy > cx:
            below[y].add(x)
    return below


def transitive_closure(below):
    out = {k: set(v) for k, v in below.items()}
    changed = True
    while changed:
        changed = False
        for x in list(out):
            add = set()
            for y in out[x]:
                add |= out.get(y, set())
            add.discard(x)
            if not add <= out[x]:
                out[x] |= add; changed = True
    return out


def comparability_components(below, nodes):
    adj = collections.defaultdict(set)
    for x, ys in below.items():
        for y in ys:
            adj[x].add(y); adj[y].add(x)
    seen = set(); comps = []
    for n in nodes:
        if n in seen:
            continue
        stack = [n]; seen.add(n); c = []
        while stack:
            u = stack.pop(); c.append(u)
            for w in adj.get(u, ()):
                if w not in seen:
                    seen.add(w); stack.append(w)
        comps.append(c)
    return comps


def p1_downsets(below, nodes):
    """principal down-set of each MAXIMAL element -- the natural cover."""
    above = collections.defaultdict(set)
    for x, ys in below.items():
        for y in ys:
            above[y].add(x)
    maximal = [n for n in nodes if not below.get(n)]
    fams = {}
    for i, m in enumerate(maximal):
        seen = {m}; stack = [m]
        while stack:
            u = stack.pop()
            for w in above.get(u, ()):
                if w not in seen:
                    seen.add(w); stack.append(w)
        if len(seen) >= 2:
            fams[f'P1_{i}'] = sorted(seen)
    return fams


def p2_antichains(below, nodes):
    """maximal antichains, greedily, inside each comparability component."""
    comp = comparability_components(below, nodes)
    fams = {}; k = 0
    for c in comp:
        if len(c) < 2:
            continue
        rel = {x: below.get(x, set()) for x in c}
        rest = sorted(c, key=lambda n: -len(rel.get(n, ())))
        while rest:
            chain = []
            for n in rest:
                if all(n not in rel.get(m, ()) and m not in rel.get(n, ()) for m in chain):
                    chain.append(n)
            if len(chain) >= 2:
                fams[f'P2_{k}'] = sorted(chain); k += 1
            rest = [n for n in rest if n not in set(chain)]
            if not chain:
                break
    return fams


def pairwise(truth, pred):
    """cover-compatible: a pair is together iff it co-occurs in ANY family.

    ⚠ Restricted to the TRUTH'S OWN UNIVERSE (nodes the truth labels). Scoring every predicted pair
    against a truth that labels only part of the node set deflates precision by the universe mismatch,
    not by the method. Conditioning on the TRUTH is correct; conditioning on the PREDICTION is register
    770's trap and is not what this does.
    """
    universe = set().union(*truth.values()) if truth else set()

    def pairs(d):
        s = set()
        for members in d.values():
            m = sorted(set(members) & universe)
            for a, b in itertools.combinations(m, 2):
                s.add((a, b))
        return s
    T, P = pairs(truth), pairs(pred)
    if not T:
        return (float('nan'),) * 3
    tp = len(T & P)
    prec = tp / len(P) if P else 0.0
    rec = tp / len(T)
    f = 0.0 if prec + rec == 0 else 2 * prec * rec / (prec + rec)
    return prec, rec, f


def main():
    ap = argparse.ArgumentParser()
    for a in ('--graph', '--paf', '--gff', '--truth'):
        ap.add_argument(a, required=True)
    ap.add_argument('--chrom', default='chr16')
    ap.add_argument('--label', default='arm')
    ap.add_argument('--sweep', default='0.50,0.60,0.70,0.80,0.90')
    a = ap.parse_args()

    adj, nodes = read_graph(a.graph)
    keep = {(u, v) if u <= v else (v, u) for u in adj for v in adj[u]}

    # locus -> gene name, max overlap (the same resolver every other scorer uses)
    genes = []
    for ln in open(a.gff):
        if ln.startswith('#'):
            continue
        f = ln.split('\t')
        if len(f) < 9 or f[0] != a.chrom or f[2] not in ('gene', 'pseudogene', 'ncRNA_gene'):
            continue
        m = re.search(r'Name=([^;]+)', f[8])
        if m:
            genes.append((int(f[3]) - 1, int(f[4]), m.group(1)))
    genes.sort()
    gs_list = [g[0] for g in genes]

    def gname(node):
        m = re.match(r'(\S+):(\d+)-(\d+)$', node)
        if not m or m.group(1) != a.chrom:
            return None
        s, e = int(m.group(2)), int(m.group(3))
        best = None
        for idx in range(max(0, bisect.bisect_left(gs_list, s) - 40), len(genes)):
            g0, g1, g = genes[idx]
            if g0 > e:
                break
            ov = min(e, g1) - max(s, g0)
            if ov > 0 and (best is None or ov > best[0]):
                best = (ov, g)
        return best[1] if best else None

    lab = {}
    for ln in open(a.truth):
        f = ln.rstrip('\n').split('\t')
        if len(f) >= 2 and f[0] != 'Gene Name':
            lab.setdefault(f[0], set()).add(f[1])
    truth = collections.defaultdict(set)
    node_gene = {n: gname(n) for n in nodes}
    for n, g in node_gene.items():
        for fam in lab.get(g, ()):
            truth[fam].add(n)
    truth = {k: v for k, v in truth.items() if len(v) >= 2}

    mcl_f = mcl_port.mcl({(u, v): w for u in adj for v, w in adj[u].items() if u < v}, inflation=2.8)
    M = {f'M_{i}': c for i, c in enumerate(mcl_f) if len(set(c)) >= 2}
    big = lambda d: max((len(set(v)) for v in d.values()), default=0)
    print(f"{a.label}: nodes {len(nodes)} | truth {len(truth)} families | "
          f"M: {len(M)} fams, largest {big(M)}")
    p, r, f = pairwise(truth, M)
    print(f"  {'M (mcl_port)':<26} pairwise P {p:.3f} R {r:.3f} F {f:.3f}   fams {len(M):>4} largest {big(M):>4}")

    for C in [float(x) for x in a.sweep.split(',')]:
        raw = directed(a.paf, keep, C)
        for tag, rel in (('raw', raw), ('tclosed', transitive_closure(raw))):
            objs = {'P1 downsets': p1_downsets(rel, nodes), 'P2 antichains': p2_antichains(rel, nodes),
                    'P3 components (null)': {f'P3_{i}': c for i, c in
                                             enumerate(comparability_components(rel, nodes)) if len(c) >= 2}}
            for name, fams in objs.items():
                p, r, f = pairwise(truth, fams)
                print(f"  C={C:.2f} {tag:<8} {name:<21} P {p:.3f} R {r:.3f} F {f:.3f}   "
                      f"fams {len(fams):>4} largest {big(fams):>4}")


if __name__ == '__main__':
    main()
