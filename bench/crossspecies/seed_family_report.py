#!/usr/bin/env python3
"""Report for the seeded-DNA family definition (bench/crossspecies/seed_family.sh).

⚠⚠ RULE CORRECTION 2026-08-11 (defect B1, report half). The docstring that used to sit here read
"Edge rule mirrors the shipped one: identity = nmatch/blocklen >= 0.80 ... E_r is the UNION of the
asm20 and the sensitive tiers." THAT WAS FALSE ON THREE AXES and every P2 density this file has ever
printed was computed under it:
  - identity metric: the binary uses `1 - de` (gap-compressed divergence), not nmatch/blocklen;
  - identity floor: the SENSITIVE tier floor is 0.60, and sensitive-only has been the default since
    2026-08-07, so 0.80 is the RETIRED asm20 leg's floor;
  - coverage form: `(qe-qs)/min(ql,tl)` is defect M1, a query-axis numerator over a possibly-target
    denominator; the shipped form's numerator axis follows the denominator.
`seed_family.sh` was corrected to issue the shipped minimap2 command on 2026-08-10, which left this
file applying a PANEL-era rule to a SHIPPED-tier PAF -- neither tier's rule.

The shipped rule now comes from ONE place, bench/soto/rustlib.py (the Python mirror of
`nucleotide_edges_scored`), so this file cannot drift from the binary again. The legacy rule is still
printed, labelled LEGACY, so the historical rows remain reproducible.
The max() denominator is reported alongside as the CERTIFICATE variant, never as the shipped rule.
"""
import sys, os
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "soto"))
import rustlib as er_tier  # noqa: E402 -- THE definition of the E_r rule

OUT = sys.argv[1]
ID_FLOOR, COV_FLOOR, GAMMA = 0.80, 0.50, 0.40


def parse_region(name):
    chrom, rng = name.rsplit(":", 1)
    a, b = rng.split("-")
    return chrom, int(a) - 1, int(b)


def shipped_edges(paf, nodes, denom="min"):
    """THE shipped E_r rule, from rustlib -- 1-de >= 0.60, M1-correct cov >= 0.50, ANY record."""
    return {frozenset(k) for k in
            er_tier.er_edges([(paf, er_tier.SENSITIVE_IDENTITY)], nodes=set(nodes),
                             min_coverage=er_tier.MIN_COVERAGE, denom=denom)}


def build_edges(paf, denom):
    """LEGACY rule (nmatch/blocklen >= 0.80, pre-M1 coverage). Retained for reproducibility ONLY."""
    edges = set()
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        q, ql, qs, qe, t, tl = f[0], int(f[1]), int(f[2]), int(f[3]), f[5], int(f[6])
        nm, bl = int(f[9]), int(f[10])
        if q == t or bl == 0:
            continue
        if nm / bl < ID_FLOOR:
            continue
        d = max(ql, tl) if denom == "max" else min(ql, tl)
        if d and (qe - qs) / d >= COV_FLOOR:
            edges.add(frozenset((q, t)))
    return edges


def components(nodes, edges):
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for e in edges:
        a, b = tuple(e)
        if a in parent and b in parent:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    groups = defaultdict(list)
    for n in nodes:
        groups[find(n)].append(n)
    return sorted(groups.values(), key=len, reverse=True)


def density(comp, edges):
    n = len(comp)
    if n < 2:
        return float("nan")
    s = set(comp)
    e = sum(1 for x in edges if set(x) <= s)
    return 2 * e / (n * (n - 1))


def load_bed(path):
    out = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3:
            out.append((f[0], int(f[1]), int(f[2]), f[3] if len(f) > 3 else ""))
    return out


def overlaps(chrom, a, b, bed):
    return [r for r in bed if r[0] == chrom and r[1] < b and a < r[2]]


def panel(tag):
    loci = [l.split()[0] for l in open(f"{OUT}/{tag}.loci.fa.fai")]
    print(f"\n{'='*66}\n{tag}: {len(loci)} candidate loci recovered from ONE seed\n{'='*66}")
    span = [parse_region(x) for x in loci]
    tot = sum(b - a for _, a, b in span)
    print(f"  total recovered sequence : {tot/1e6:.2f} Mb")
    print(f"  locus span median/min/max: "
          f"{sorted(b-a for _,a,b in span)[len(span)//2]:,} / "
          f"{min(b-a for _,a,b in span):,} / {max(b-a for _,a,b in span):,} bp")
    # PROVENANCE. Every density below is a statement about ONE tier AND ONE rule; print both.
    print(f"  E_r TIER: {er_tier.er_provenance()}")
    res = {}
    for rule in ("SHIPPED", "LEGACY"):
        for denom in ("min", "max"):
            if rule == "SHIPPED":
                edges = shipped_edges(f"{OUT}/{tag}.ava.paf", loci, denom)
                rulestr = (f"1-de >= {er_tier.SENSITIVE_IDENTITY:.2f}, "
                           f"cov-of-{denom}(M1) >= {er_tier.MIN_COVERAGE:.2f}, ANY record")
            else:
                edges = build_edges(f"{OUT}/{tag}.ava.paf", denom)
                rulestr = f"nmatch/blocklen >= {ID_FLOOR:.2f}, cov(pre-M1)/{denom} >= {COV_FLOOR:.2f}"
            edges = {e for e in edges if set(e) <= set(loci)}
            comps = components(loci, edges)
            big = comps[0]
            d = density(big, edges)
            print(f"\n  -- {rule} rule, coverage denominator = {denom}() --   [{rulestr}]")
            print(f"     edges {len(edges)}   components {len(comps)}   largest {len(big)}/{len(loci)}")
            print(f"     induced density of largest component = {d:.3f}"
                  f"   (gamma=0.40 margin {d/0.40:.2f}x; shipped gamma=0.20 margin {d/0.20:.2f}x)"
                  if d == d else "     density n/a")
            singles = sum(1 for c in comps if len(c) == 1)
            print(f"     singletons {singles}")
            res[(rule, denom)] = (comps, edges)
    # Downstream validation blocks keep reading res["min"]; bind it to the SHIPPED rule.
    res["min"] = res[("SHIPPED", "min")]
    res["max"] = res[("SHIPPED", "max")]
    return loci, res


# ------------------------------------------------------------------ HUMAN panel (control)
hs_loci, hs_res = panel("HS")
hs_bed = load_bed(f"{OUT}/hs_npip.bed")
big = hs_res["min"][0][0]
hit_genes, off_target = set(), []
for L in big:
    c, a, b = parse_region(L)
    ov = overlaps(c, a, b, hs_bed)
    if ov:
        hit_genes.update(g[3] for g in ov)
    else:
        off_target.append(L)
allg = {g[3] for g in hs_bed}
print(f"\n  VALIDATION vs CHM13 annotation (19 NPIP genes, none of which were given except the seed):")
print(f"     NPIP genes recovered in the largest component: {len(hit_genes)}/{len(allg)}")
if allg - hit_genes:
    print(f"     missed: {', '.join(sorted(allg - hit_genes))}")
print(f"     loci in the component NOT overlapping any annotated NPIP gene: "
      f"{len(off_target)}/{len(big)}")

# ------------------------------------------------------------------ GORILLA panel
ggo_loci, ggo_res = panel("GGO")
# ⚠ Counted by QUALIFYING EDGE, not by argmax best-hit.  Best-hit is an argmax and it DEFLATES
#   here: every gorilla locus's single best target is one of ~6 human genes, which reads as
#   "6/19 orthologs" when in fact all 19 clear the edge rule.
name_of = {}
bed19 = load_bed(f"{OUT}/hs_npip.bed")
for i, line in enumerate(open(f"{OUT}/hs_npip.regions")):
    name_of[line.strip()] = bed19[i][3]

qual, ids_all = defaultdict(set), []
for line in open(f"{OUT}/GGO_vs_hsnpip.paf"):
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        continue
    q, ql, qs, qe, t, tl = f[0], int(f[1]), int(f[2]), int(f[3]), f[5], int(f[6])
    nm, bl = int(f[9]), int(f[10])
    if bl == 0:
        continue
    idy, cov = nm / bl, (qe - qs) / max(min(ql, tl), 1)
    if idy >= ID_FLOOR and cov >= COV_FLOOR:
        qual[name_of.get(t, t)].add(q)
        ids_all.append(idy)

gbig = ggo_res["min"][0][0]
linked = set().union(*qual.values()) if qual else set()
print(f"\n  VALIDATION vs HUMAN NPIP genes (independent of the gorilla annotation we seeded from):")
print(f"     gorilla loci in the largest component with a qualifying human-NPIP edge: "
      f"{len(linked & set(gbig))}/{len(gbig)}")
if ids_all:
    ids_all.sort()
    print(f"     identity median {ids_all[len(ids_all)//2]:.3f}  "
          f"range {ids_all[0]:.3f}-{ids_all[-1]:.3f}")
allg19 = {r[3] for r in bed19}
print(f"     human NPIP genes with >=1 qualifying gorilla edge: {len(qual)}/{len(allg19)}")
if allg19 - set(qual):
    print(f"     none: {', '.join(sorted(allg19 - set(qual)))}")

named = defaultdict(set)
for line in open(f"{OUT}/GGO_loci_genes.tsv"):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 8:
        named[f"{f[0]}:{int(f[1])+1}-{f[2]}"].add(f[7])
allnames = set().union(*named.values()) if named else set()
real = {n for n in allnames if n and not n.startswith("LOC")}
print(f"\n  WHAT THE GORILLA ANNOTATION CALLS THEM:")
print(f"     loci overlapping >=1 annotated gene: {len(named)}/{len(ggo_loci)}")
print(f"     distinct gene symbols: {len(allnames)}   named (non-LOC): {len(real)}")
print(f"     named symbols: {', '.join(sorted(real)) if real else '(none)'}")
print(f"\n  => gorilla annotation supplied 1 NPIP name; the method returned "
      f"{len(gbig)} loci in one component.")
