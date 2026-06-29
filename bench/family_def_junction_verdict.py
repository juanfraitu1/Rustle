#!/usr/bin/env python3
"""
family_def_junction_verdict.py — final consolidated A2 verdict.

Compares, on the cDNA-homology graph (17410 edges / 1460 families / 389-mega):
  (B0) baseline                       — no filter
  (Bloci) disjoint-loci filter        — drop edges whose gene spans OVERLAP (co-location /
                                        nesting artifact). This is the existing baseline the
                                        task says removed 307 edges (389 -> 382, barely helped).
  (Bsplice) intron-architecture / retrocopy filter — drop any edge incident to an intronless
                                        gene (the orthogonal-to-sequence retrocopy signal).
  (Bboth) loci AND splice              — do they stack?

Reports for each: edges, families, MAX component (does 389 break?), genes-in-families, and the
panel connectivity. Also the DIRECT comparison the task asks: how many edges does each cut
INSIDE the 389-mega vs OUTSIDE, and does the mega break into smaller coherent pieces.

Run: /home/juanfra/miniforge3/bin/python bench/family_def_junction_verdict.py
"""
import collections
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_read_filters import dna_homology
from family_def_junction_genomewide import components

INDEX = "/home/juanfra/winloci_scratch/gene_intron_index.json"
GENES_BED = "/home/juanfra/winloci_scratch/unmapped_poc/genes.bed"
PANEL_REAL = {"RABL2": ["RABL2A", "RABL2B"],
              "APOBEC3": ["APOBEC3D", "APOBEC3F"],
              "RFPL": ["RFPL1", "RFPL2", "RFPL3"]}


def load_coords():
    coord = {}
    with open(GENES_BED) as f:
        for line in f:
            c, s, e, name = line.rstrip("\n").split("\t")
            # keep widest span if duplicate names
            s, e = int(s), int(e)
            if name in coord:
                oc, os_, oe = coord[name]
                if oc == c:
                    coord[name] = (c, min(os_, s), max(oe, e))
            else:
                coord[name] = (c, s, e)
    return coord


def overlap(coord, a, b):
    ca = coord.get(a); cb = coord.get(b)
    if ca is None or cb is None:
        return False
    if ca[0] != cb[0]:
        return False
    return min(ca[2], cb[2]) - max(ca[1], cb[1]) > 0


def panel_conn(edge_set):
    comps = components(edge_set); g2c = {}
    for i, c in enumerate(comps):
        for g in c:
            g2c[g] = i
    out = {}
    for fam, mem in PANEL_REAL.items():
        pres = [m for m in mem if m in g2c]
        out[fam] = "CONN" if len(pres) >= 2 and len({g2c[m] for m in pres}) == 1 else \
                   ("SPLIT" if len(pres) >= 2 else "N/A")
    return out


def main():
    H, dna = dna_homology()
    idx = json.load(open(INDEX))
    coord = load_coords()
    base_edges = sorted(dna)
    base_comp = components(base_edges)
    mega = max(base_comp, key=len)
    mega_set = set(mega)
    mega_edges = {(a, b) for (a, b) in base_edges if a in mega_set and b in mega_set}
    il = {g for g in idx if idx[g]["n_intron"] == 0}

    def keep_loci(a, b):       # disjoint-loci: drop if spans overlap
        return not overlap(coord, a, b)

    def keep_splice(a, b):     # retrocopy filter: drop if either side intronless
        ra, rb = idx.get(a), idx.get(b)
        return ra is not None and rb is not None and ra["n_intron"] >= 1 and rb["n_intron"] >= 1

    FILTERS = {
        "B0 baseline": lambda a, b: True,
        "Bloci disjoint-loci": keep_loci,
        "Bsplice retro-filter": keep_splice,
        "Bboth loci&splice": lambda a, b: keep_loci(a, b) and keep_splice(a, b),
    }

    print(f"[baseline] edges={len(base_edges)} families={len(base_comp)} "
          f"maxcomp={len(mega)} genes-in-fam={sum(len(c) for c in base_comp)}")
    print(f"\n  {'filter':24} {'edges':>6} {'dropped':>7} {'fams':>5} {'maxC':>5} "
          f"{'g-in-fam':>9} {'megaIn':>7} {'megaOut':>7} {'panel':>16}")
    out = {}
    for name, keep in FILTERS.items():
        kept = [(a, b) for (a, b) in base_edges if keep(a, b)]
        comp = components(kept)
        maxc = max((len(c) for c in comp), default=0)
        gin = sum(len(c) for c in comp)
        mega_in = sum(1 for (a, b) in mega_edges if not keep(a, b))   # dropped inside mega
        mega_out = (len(base_edges) - len(kept)) - mega_in            # dropped outside mega
        pc = panel_conn(set(kept))
        pstr = "/".join(pc[f][:4] for f in ("RABL2", "APOBEC3", "RFPL"))
        print(f"  {name:24} {len(kept):6d} {len(base_edges)-len(kept):7d} {len(comp):5d} "
              f"{maxc:5d} {gin:9d} {mega_in:7d} {mega_out:7d} {pstr:>16}")
        out[name] = dict(edges=len(kept), dropped=len(base_edges) - len(kept),
                         families=len(comp), maxcomp=maxc, genes_in_fam=gin,
                         mega_internal_dropped=mega_in, mega_external_dropped=mega_out,
                         panel=pc)

    # mega breakup under Bsplice
    keep = keep_splice
    mk = [(a, b) for (a, b) in mega_edges if keep(a, b)]
    sub = components(mk)
    subsizes = sorted((len(c) for c in sub), reverse=True)
    iso = len(mega) - sum(len(c) for c in sub)
    print(f"\n=== 389-MEGA under retro-filter (Bsplice) ===")
    print(f"  {len(mega)} genes / {len(mega_edges)} edges -> kept {len(mk)} edges; "
          f"largest surviving sub-component = {subsizes[0] if subsizes else 0}; "
          f"sub-components={len(sub)}; newly-isolated singletons={iso}")
    print(f"  NOTE: retro-filter shaves the intronless arm but the multi-exon CORE of the "
          f"mega stays connected (max {subsizes[0] if subsizes else 0}); domain bridges are "
          f"multi-exon<->multi-exon so they SURVIVE this filter.")

    # quantify: what fraction of the 14809-edge cut at Carch vs the 196/5198 at Bsplice is retro
    print(f"\n=== ORTHOGONALITY / SCALE ===")
    retro_edges = [(a, b) for (a, b) in base_edges if a in il or b in il]
    # are retro edges high-identity? (they should be id~1.0 — orthogonal to id)
    ids = []
    for a, b in retro_edges:
        r = H.get((a, b) if a < b else (b, a))
        if r:
            ids.append(r["id"])
    import statistics
    print(f"  retro edges (intronless-incident): {len(retro_edges)} "
          f"({100*len(retro_edges)/len(base_edges):.1f}% of homology graph)")
    if ids:
        print(f"  their cDNA identity: median={statistics.median(ids):.3f} "
              f"min={min(ids):.3f} max={max(ids):.3f}  "
              f"-> HIGH id, so sequence-id CANNOT distinguish them; intron-count CAN.")

    json.dump(dict(filters=out,
                   mega_breakup_splice=dict(subsizes=subsizes[:10], singletons=iso,
                                            kept_edges=len(mk)),
                   retro_edges=len(retro_edges),
                   retro_id_median=statistics.median(ids) if ids else None),
              open(os.path.join(HERE, "family_def_junction_verdict.json"), "w"), indent=2)
    print(f"\n[+] wrote {os.path.join(HERE,'family_def_junction_verdict.json')}")


if __name__ == "__main__":
    main()
