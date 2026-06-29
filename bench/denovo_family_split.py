#!/usr/bin/env python3
"""Decompose the over-merged de-novo family components into real recent-duplicate families.

WHY: connected components of the POA-homology graph transitively close a SUPERFAMILY (e.g. KRAB zinc
fingers) into one blob, because genes that merely share a tandem DOMAIN pairwise-confirm. A real
recent-duplicate family is a DENSE, mutually-homologous subgraph; the over-merge is a SPARSE chain bridged
by lower-core_recip domain edges. So we keep the same POA edges but report each component's structure and
split large components by weighted-modularity (uses the core_recip signal; no arbitrary similarity cutoff).

Reads denovo_family_edges.tsv (a, b, core_recip) written by denovo_families.py -- NO POA re-run.
Run with /home/juanfra/miniforge3/bin/python (networkx). Deterministic (seeded).
"""
import os
import re
import sys
from collections import defaultdict, Counter

import networkx as nx

HERE = os.path.dirname(__file__)
EDGES = os.path.join(HERE, "denovo_family_edges.tsv")
META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
ANNOT = "/home/juanfra/winloci_scratch/annot_intervals.tsv"
OUT = os.path.join(HERE, "denovo_families_split.tsv")
MIN_DECOMP = 6        # only decompose components this size or larger (small families left intact)
RESOLUTION = 1.0      # modularity resolution (1.0 = standard, principled default)
SEED = 0

KNOWN = {"RABL2": [("NC_073235.2", 15131675, 15147533), ("NC_086018.1", 48818440, 48831690)],
         "APOBEC3": [("NC_086018.1", 37023493, 37058621)],
         "RFPL": [("NC_086018.1", 27445867, 30376055)]}


def load_annot():
    meta = {}
    for line in open(META):
        if line.startswith("id"):
            continue
        tid, c, s, e, st, ne, nr = line.rstrip("\n").split("\t")
        meta[tid] = (c, int(s), int(e))
    ann = defaultdict(list)
    for line in open(ANNOT):
        if line.startswith("chrom"):
            continue
        c, s, e, bt, g = line.rstrip("\n").split("\t")
        ann[c].append((int(s), int(e), g))
    for c in ann:
        ann[c].sort()
    return meta, ann


def gene_root(tid, meta, ann):
    c, s, e = meta[tid]
    best, bov = None, 0
    for as_, ae, g in ann.get(c, []):
        ov = min(e, ae) - max(s, as_)
        if ov > bov:
            bov, best = ov, g
    if best is None or best.startswith("LOC"):
        return None
    return re.sub(r"[0-9]+[A-Z]?$", "", best)   # TUBA1C -> TUBA, ZNF577 -> ZNF


def dedup_oversplit(members, meta, thr=0.50):
    """OVER-SPLIT GUARD (output-level). The intron-junction locus collapse leaves
    near-identical-junction variants at the SAME genomic position as separate "copies"
    (e.g. PRNP's five 14.60-14.615 Mb fragments). Members whose spans reciprocally
    overlap >= thr are one locus -> merge, keep the longest-span representative.
    Applied at OUTPUT only: family DETECTION/decomposition is untouched, so genuinely
    distinct paralog copies (non-overlapping spans) are preserved -- unlike merging
    BEFORE POA, which perturbs detection and drops real copies (MAGEA 11->6)."""
    ms = list(members)
    par = list(range(len(ms)))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    by_c = defaultdict(list)
    for i, t in enumerate(ms):
        by_c[meta[t][0]].append(i)
    for _c, idxs in by_c.items():
        idxs.sort(key=lambda i: meta[ms[i]][1])
        for a in range(len(idxs)):
            ia = idxs[a]; sa, ea = meta[ms[ia]][1], meta[ms[ia]][2]
            for b in range(a + 1, len(idxs)):
                ib = idxs[b]; sb, eb = meta[ms[ib]][1], meta[ms[ib]][2]
                if sb >= ea:                       # sorted by start -> no further overlap
                    break
                ov = min(ea, eb) - max(sa, sb)
                if ov > 0 and ov >= thr * min(ea - sa, eb - sb):
                    par[f(ia)] = f(ib)
    groups = defaultdict(list)
    for i in range(len(ms)):
        groups[f(i)].append(ms[i])
    return {max(g, key=lambda t: meta[t][2] - meta[t][1]) for g in groups.values()}


def concordance(members, meta, ann):
    roots = [gene_root(t, meta, ann) for t in members]
    roots = [r for r in roots if r]
    if len(roots) < 2:
        return None, None
    top, cnt = Counter(roots).most_common(1)[0]
    return top, cnt / len(roots)


def family_quality(fams, meta, ann, label):
    judged = conc = 0
    sizes = sorted((len(m) for m in fams), reverse=True)
    for m in fams:
        top, c = concordance(m, meta, ann)
        if c is None:
            continue
        judged += 1
        if c >= 0.5:
            conc += 1
    print(f"  [{label}] families={len(fams)} maxsize={sizes[0] if sizes else 0} "
          f">100={sum(1 for s in sizes if s > 100)} 21-100={sum(1 for s in sizes if 21 <= s <= 100)} "
          f"6-20={sum(1 for s in sizes if 6 <= s <= 20)} | name-concordant {conc}/{judged} "
          f"({conc/judged:.1%})" if judged else f"  [{label}] families={len(fams)}")


def main():
    meta, ann = load_annot()
    G = nx.Graph()
    for line in open(EDGES):
        if line.startswith("a\t"):
            continue
        a, b, cr = line.rstrip("\n").split("\t")
        G.add_edge(a, b, weight=float(cr))
    print(f"edge graph: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")

    raw = [c for c in nx.connected_components(G) if len(c) >= 2]
    print("\nBEFORE decomposition (connected components = current families):")
    family_quality(raw, meta, ann, "raw")

    # structural diagnostics on the largest components (dense superfamily vs sparse domain-chain?)
    print("\nstructure of the largest raw components:")
    for comp in sorted(raw, key=len, reverse=True)[:4]:
        sub = G.subgraph(comp)
        n, m = sub.number_of_nodes(), sub.number_of_edges()
        dens = 2 * m / (n * (n - 1)) if n > 1 else 0
        arts = len(list(nx.articulation_points(sub)))
        bicon = sum(1 for _ in nx.biconnected_components(sub))
        avgw = sum(d["weight"] for *_, d in sub.edges(data=True)) / m
        top, c = concordance(comp, meta, ann)
        print(f"  n={n:4d} edges={m:6d} density={dens:.3f} avg_core_recip={avgw:.3f} "
              f"articulation_pts={arts} biconnected_comps={bicon} dominant={top}({c:.0%})"
              if c else f"  n={n} density={dens:.3f}")

    # decompose components >= MIN_DECOMP by weighted modularity; keep smaller ones intact
    final = []
    n_decomposed = 0
    for comp in raw:
        if len(comp) < MIN_DECOMP:
            final.append(set(comp))
            continue
        sub = G.subgraph(comp)
        communities = nx.community.louvain_communities(sub, weight="weight", resolution=RESOLUTION, seed=SEED)
        communities = [set(c) for c in communities if len(c) >= 2]
        if len(communities) > 1:
            n_decomposed += 1
        final.extend(communities)
    print(f"\ndecomposed {n_decomposed} large components by weighted modularity (resolution={RESOLUTION})")
    print("AFTER decomposition:")
    family_quality(final, meta, ann, "split")

    # write split families with structural diagnostics. A discrete recent-duplicate family is DENSE and
    # mutually homologous; a non-discretizable homology web (e.g. the KRAB-ZNF continuum) is SPARSE with
    # many articulation points. We FLAG webs (size>=10 & density<0.30) rather than claim them as families.
    # NB: 0.30 ALIGNS with the validated DNA-manifest bar (make_dna_family_manifest.py overmerge_sparse)
    # and with src/rustle/vg_family/family_split.rs::WEB_MAX_DENSITY. Measured: the family-density
    # distribution is bimodal (median 1.0); the only families in [0.15,0.30) are multi-chromosome
    # domain-sharing over-merges (DSFAM0 = 164-member ZNF / 19 chroms). web_min_size stays 10 so a SMALL
    # sparse group (a real divergent family, e.g. MAGEB n=7 density 0.24) is NOT dropped. See L2.
    final.sort(key=len, reverse=True)
    memb2fam = {}
    n_web = 0
    n_written = 0
    n_oversplit_dropped = 0
    with open(OUT, "w") as fh:
        fh.write("family_id\tn_copies\tn_chroms\tdominant_gene_root\tname_concordance\t"
                 "density\tavg_core_recip\tn_articulation\tclass\tmembers\n")
        for i, m in enumerate(final):
            m = dedup_oversplit(m, meta)            # collapse over-split fragments (output-level)
            if len(m) < 2:                          # all members were one over-split locus (e.g. PRNP)
                n_oversplit_dropped += 1
                continue
            top, c = concordance(m, meta, ann)
            nchr = len({meta[t][0] for t in m})
            sub = G.subgraph(m)
            n, e = sub.number_of_nodes(), sub.number_of_edges()
            dens = 2 * e / (n * (n - 1)) if n > 1 else 1.0
            avgw = (sum(d["weight"] for *_, d in sub.edges(data=True)) / e) if e else 0.0
            arts = len(list(nx.articulation_points(sub))) if n > 2 else 0
            klass = "web" if (n >= 10 and dens < 0.30) else "family"  # 0.30 = aligned DNA-manifest bar (L2)
            if klass == "web":
                n_web += 1
            fid = f"DSFAM{i}"
            for t in m:
                memb2fam[t] = fid
            n_written += 1
            fh.write(f"{fid}\t{len(m)}\t{nchr}\t{top or '(novel)'}\t{c if c is not None else 'NA'}\t"
                     f"{dens:.3f}\t{avgw:.3f}\t{arts}\t{klass}\t{','.join(sorted(m))}\n")
    print(f"\n[wrote {OUT}]  ({n_written-n_web} discrete families + {n_web} flagged homology webs; "
          f"{n_oversplit_dropped} over-split single-locus 'families' dropped)")

    # known-family recovery after split (do they sit in clean, dedicated families now?)
    print("\nknown-family recovery after split:")
    for name, regions in KNOWN.items():
        hits = defaultdict(int)
        for tid, fid in memb2fam.items():
            c, s, e = meta[tid]
            for kc, ks, ke in regions:
                if c == kc and s < ke and e > ks:
                    hits[fid] += 1
        if not hits:
            print(f"  {name}: none in any family")
            continue
        parts = ", ".join(f"{fid}(n={n}, famsize={sum(1 for x in memb2fam.values() if x == fid)})"
                          for fid, n in sorted(hits.items(), key=lambda kv: -kv[1]))
        print(f"  {name}: {parts}")


if __name__ == "__main__":
    main()
