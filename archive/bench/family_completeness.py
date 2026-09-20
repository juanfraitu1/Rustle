#!/usr/bin/env python3
"""Completeness sanity check for the annotation-based RNA-level family definition: build families as
connected components over the POA homology graph (bench/family_pairs_graded.tsv) and verify that KNOWN
multi-copy gene families are recovered (their members land in one component). Reports recall on a
hand-curated textbook set + the similarity-built universe, the component size distribution, and the
over-merge (mega-component / domain-hub chain) control. Annotated-only; extensible to unannotated.
"""
import csv
import os
from collections import defaultdict

BASE = os.path.dirname(__file__)
META = "/tmp/gene_reps_gw.meta.tsv"
GRADED = os.path.join(BASE, "family_pairs_graded.tsv")
UNIVERSE = os.path.join(BASE, "copy_recovery_eval/results/universe.tsv")
OUT_MD = os.path.join(BASE, "family_definition.md")
OUT_TSV = os.path.join(BASE, "families.tsv")
CORE_COV = 0.13       # homology edge: contiguous-core coverage
CORE_IDENT = 0.70     # homology edge: core-block identity (drops chance/very-distant)

# hand-curated high-confidence textbook multi-copy families (exact members where prefix is ambiguous;
# specific prefixes [TAS2R/DEFB/MAGEA/PRAMEF] expanded against the gene set).
CURATED_EXACT = {
    "RABL2": ["RABL2A", "RABL2B"],
    "APOBEC3": ["APOBEC3B", "APOBEC3C", "APOBEC3D", "APOBEC3F", "APOBEC3G", "APOBEC3H"],
    "RFPL": ["RFPL1", "RFPL2", "RFPL3", "RFPL4A", "RFPL4B"],
    "DAZ": ["DAZ1", "DAZL"],
    "GGT": ["GGT1", "GGT5", "GGT6", "GGT7", "GGTLC1", "GGTLC2"],
    "CRYBG": ["CRYBG1", "CRYBG2", "CRYBG3"],
}
CURATED_PREFIX = ["TAS2R", "DEFB", "MAGEA", "PRAMEF", "OR2", "OR4", "SIGLEC"]


def load_genes():
    g = set()
    with open(META) as fh:
        fh.readline()
        for line in fh:
            g.add(line.split("\t", 1)[0])
    return g


def main():
    genes = load_genes()
    # build known sets
    curated = {}
    for fam, members in CURATED_EXACT.items():
        present = [m for m in members if m in genes]
        if len(present) >= 2:
            curated[fam] = present
    for pre in CURATED_PREFIX:
        present = sorted(x for x in genes if x.startswith(pre))
        if len(present) >= 2:
            curated[pre + "*"] = present
    # universe families
    uni = defaultdict(set)
    with open(UNIVERSE) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if int(r["n_family_copies"]) >= 2:
                uni[r["family_id"]].add(r["gene_id"])
    uni = {f: sorted(m & genes) for f, m in uni.items() if len(m & genes) >= 2}

    # build family graph (connected components) at two edge stringencies
    def components(min_ident):
        parent = {}
        def find(x):
            parent.setdefault(x, x)
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        nedge = 0
        with open(GRADED) as fh:
            fh.readline()
            for line in fh:
                a, ca, b, cb, cc, ci, sc = line.rstrip("\n").split("\t")
                if float(cc) >= CORE_COV and float(ci) >= min_ident:
                    parent[find(a)] = find(b); nedge += 1
        comp = defaultdict(set)
        for g in list(parent):
            comp[find(g)].add(g)
        return [c for c in comp.values()], nedge

    comps, nedge = components(CORE_IDENT)
    g2c = {}
    for i, c in enumerate(comps):
        for g in c:
            g2c[g] = i

    def recall(known):
        rec = miss = 0; missed = []
        for fam, members in known.items():
            comps_hit = defaultdict(int)
            for m in members:
                if m in g2c:
                    comps_hit[g2c[m]] += 1
            best = max(comps_hit.values()) if comps_hit else 0
            if best >= 2:
                rec += 1
            else:
                miss += 1; missed.append(fam)
        return rec, miss, missed

    cr, cm, cmiss = recall(curated)
    ur, um, umiss = recall(uni)
    sizes = sorted((len(c) for c in comps), reverse=True)
    mega = [c for c in comps if len(c) >= 25]

    L = []
    def P(s=""): L.append(s)
    P("# Annotation-based RNA-level multi-copy gene family definition (completeness)")
    P()
    P("**Definition:** a multi-copy gene family = a maximal connected set of ANNOTATED genes whose "
      f"representative transcripts pairwise share a POA contiguous exon core ≥ {CORE_COV} at core "
      f"identity ≥ {CORE_IDENT}. Roster = the {len(genes):,} annotated gene models; grouping = the "
      "validated POA criterion; families = connected components. Annotated-only (extensible: the same "
      "criterion runs against de-novo/unannotated loci later to add unannotated copies).")
    P()
    P(f"- homology edges: {nedge:,} ; **families (components, size≥2): "
      f"{sum(1 for c in comps if len(c) >= 2):,}** ; genes in a family: "
      f"{sum(len(c) for c in comps if len(c) >= 2):,}")
    P(f"- component size distribution (top): {sizes[:12]}")
    P()
    P("## Completeness — are KNOWN multi-copy families recovered?")
    P(f"- **curated textbook families: {cr}/{cr+cm} recovered** "
      f"(members land in one component){' ; missed: ' + ', '.join(cmiss) if cmiss else ''}")
    P(f"- similarity-built universe families: **{ur}/{ur+um} recovered**"
      f"{' ; missed (sample): ' + ', '.join(umiss[:8]) if umiss else ''}")
    P()
    P("### curated families & their recovery")
    P("| family | members (annotated) | in one component? |")
    P("|---|---|---|")
    for fam, members in sorted(curated.items()):
        ch = defaultdict(int)
        for m in members:
            if m in g2c:
                ch[g2c[m]] += 1
        best = max(ch.values()) if ch else 0
        P(f"| {fam} | {len(members)} ({', '.join(members[:6])}{'…' if len(members) > 6 else ''}) | "
          f"{'YES (' + str(best) + '/' + str(len(members)) + ')' if best >= 2 else 'no'} |")
    P()
    P("## Over-merge control (the precision side of 'complete')")
    P(f"- mega-components (size≥25, likely domain-hub chains, NOT single families): {len(mega)} "
      f"(largest {sizes[0] if sizes else 0})")
    P("- these are the transitive-closure artifacts; a tighter clustering (mutual-core / community "
      "detection) would split them. Size-2/3 components are clean copy sets.")
    P()
    P("## Honest scope")
    P("- ANNOTATED genes only (this definition); unannotated/de-novo copies are added later by running "
      "the SAME POA criterion against de-novo loci (the cross-chrom discovery already does this).")
    P("- RNA-structural operational family (shared contiguous exon core), NOT the formal gene-tree/"
      "protein family — a parallel DNA tier would supply that (and the ground truth).")
    P("- universe is itself similarity-built (recovery there is partly internal consistency); the "
      "curated textbook set is the independent completeness check.")

    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(L) + "\n")
    with open(OUT_TSV, "w") as fh:
        fh.write("family_id\tn_genes\tgenes\n")
        for i, c in enumerate(sorted(comps, key=lambda c: -len(c))):
            if len(c) >= 2:
                fh.write(f"FAM{i}\t{len(c)}\t{','.join(sorted(c))}\n")
    print("\n".join(L))
    print(f"\n[wrote {OUT_MD} and {OUT_TSV}]")


if __name__ == "__main__":
    main()
