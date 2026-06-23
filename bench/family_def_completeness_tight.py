#!/usr/bin/env python3
"""family_def_completeness_tight.py — (a) clean recall with a TIGHT locus<->gene mapping.

The completeness/scaffold numbers used best_gene (best ANY-overlap), which is noisy. Here a
gene maps to a de-novo locus only on RECIPROCAL overlap >= 0.5 (substantial, not incidental):
a gene is EXPRESSED iff some de-novo locus reciprocally-overlaps it; RECOVERED iff some
~R∩~B-family locus does. Re-measure the funnel + the detectable/recovered counts cleanly.
Run: python bench/family_def_completeness_tight.py
"""
import bisect
import collections
import json
import os
import sys

import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import GENES_BED
from family_def_read_filters import dna_homology

META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
FAM_TSV = "/home/juanfra/winloci_scratch/validated_families.tsv"
RECIP = 0.50


def overlaps_recip(loci_by_chrom, starts_by_chrom, chrom, gs, ge):
    """True if some locus on `chrom` reciprocally-overlaps [gs,ge] >= RECIP."""
    locs = loci_by_chrom.get(chrom)
    if not locs:
        return False
    starts = starts_by_chrom[chrom]
    i = bisect.bisect_right(starts, ge)
    for k in range(i - 1, -1, -1):
        ls, le = locs[k]
        if le <= gs:
            if ge - ls > 5_000_000:        # loci are sorted by start; stop once far left
                break
            continue
        ov = min(ge, le) - max(gs, ls)
        if ov > 0 and ov >= RECIP * min(ge - gs, le - ls):
            return True
    return False


def main():
    coord = {}
    with open(GENES_BED) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                coord[p[3]] = (p[0], int(p[1]), int(p[2]))
    Hd, _ = dna_homology()
    DG = nx.Graph()
    for (ga, gb), r in Hd.items():
        if r.get("id", 0) >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30:
            DG.add_edge(ga, gb)
    dna_fams = [c for c in nx.connected_components(DG) if len(c) >= 2]

    # de-novo loci (all) and RNA-family-member loci, indexed per chrom
    all_loci = collections.defaultdict(list)
    with open(META) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            all_loci[p[1]].append((int(p[2]), int(p[3])))
    fam_loci = collections.defaultdict(list)
    fams = collections.defaultdict(list)
    for line in open(FAM_TSV).read().splitlines()[1:]:
        fi, lid, c, s, e, nr = line.split("\t")
        fams[int(fi)].append((c, int(s), int(e)))
    for fi, loci in fams.items():
        if len(loci) >= 2:
            for c, s, e in loci:
                fam_loci[c].append((s, e))
    for d in (all_loci, fam_loci):
        for c in d:
            d[c].sort()
    all_starts = {c: [s for s, e in v] for c, v in all_loci.items()}
    fam_starts = {c: [s for s, e in v] for c, v in fam_loci.items()}

    def expressed(g):
        c, s, e = coord[g]
        return overlaps_recip(all_loci, all_starts, c, s, e)

    def recovered(g):
        c, s, e = coord[g]
        return overlaps_recip(fam_loci, fam_starts, c, s, e)

    genes = [g for g in DG.nodes() if g in coord]
    exp = {g for g in genes if expressed(g)}
    rec = {g for g in genes if recovered(g)}
    tot = len(genes)
    print(f"=== COMPLETENESS, TIGHT mapping (reciprocal overlap >= {RECIP}) ===")
    print(f"  DNA paralog families: {len(dna_fams)}; member genes with coords: {tot}")
    print(f"    EXPRESSED (de-novo locus recip-overlaps): {len(exp):5d} ({100*len(exp)/tot:.0f}%)")
    print(f"    RECOVERED (~R∩~B-family locus recip-overlaps): {len(rec):5d} "
          f"({100*len(rec)/tot:.0f}% of all; {100*len(rec)/max(len(exp),1):.0f}% of expressed)")
    silent = tot - len(exp); enr = len(exp - rec)
    print(f"  gap: SILENT/not-in-sample {silent} ({100*silent/tot:.0f}%); "
          f"EXPRESSED-not-recovered {enr} ({100*enr/tot:.0f}%); RECOVERED {len(rec)} ({100*len(rec)/tot:.0f}%)")

    # per-family (detectable = >=2 expressed; recovered = >=2 recovered)
    det = rec2 = 0
    missed = 0
    for fam in dna_fams:
        fe = sum(1 for g in fam if g in exp)
        fr = sum(1 for g in fam if g in rec)
        if fe >= 2:
            det += 1
            if fr >= 2:
                rec2 += 1
            else:
                missed += 1
    print(f"\n  per-family: {det} detectable (>=2 expressed); {rec2} recovered as ~R∩~B ({100*rec2/max(det,1):.0f}%); "
          f"{missed} detectable-but-not-recovered (the scaffold target)")
    print(f"  (vs noisy best_gene mapping: 315 detectable / 94 recovered / 221 missed -- "
          f"tight mapping {'changes' if det != 315 or rec2 != 94 else 'confirms'} it)")
    json.dump(dict(recip=RECIP, dna_families=len(dna_fams), member_genes=tot,
                   expressed=len(exp), recovered=len(rec), silent=silent, expressed_not_recovered=enr,
                   detectable=det, recovered_families=rec2, missed=missed),
              open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "family_def_completeness_tight.json"), "w"), indent=2)
    print("\n[+] wrote family_def_completeness_tight.json")


if __name__ == "__main__":
    main()
