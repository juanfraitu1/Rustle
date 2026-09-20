#!/usr/bin/env python
"""Per-block gate test for the 57 E_p-impure residual blocks.
For each block: verdict (genuine/oversplit/cardinality/no_multiexon), TYPE tag, and whether
a (i) naive E_p-purity gate and (ii) refined [E_p-impure AND colinear<=1 AND not-cardinality] gate
would CUT it. Goal: refined gate must cut the 32 genuine, spare the 23 truth-artifacts."""
import os, re
from collections import Counter, defaultdict

BENCH = "/mnt/c/Users/jfris/Desktop/Rustle/bench"
SCRATCH = "/home/juanfra/winloci_scratch"
ANNOT = os.path.join(SCRATCH, "annot_intervals.tsv")
BLK = "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/03f0156d-63b6-4d0f-bf09-862bddee491c/scratchpad/blocks57.txt"

biotype = {}
with open(ANNOT) as fh:
    next(fh)
    for ln in fh:
        f = ln.rstrip("\n").split("\t"); g = f[4]
        if g not in biotype or f[3] == "protein_coding":
            biotype[g] = f[3]

VERDICTS = {"genuine_overmerge","ep_oversplit","cardinality","no_multiexon_edge"}
NONCODING = {"lncRNA","snRNA","snoRNA","tRNA","rRNA","miRNA","misc_RNA","ncRNA","V_segment","C_region"}
PSEUDO = {"pseudogene","transcribed_pseudogene"}

rows = []
for ln in open(BLK):
    p = ln.split()
    if len(p) < 7 or p[-1] not in VERDICTS: continue
    verdict = p[-1]; sedef = p[-2]
    # cols: block genes nLoci maxCol mmid core sedef verdict
    try:
        core = float(p[-3]); mmid = float(p[-4]); maxcol = int(p[-5]); nloci = int(p[-6])
    except ValueError:
        continue
    blk = p[0]; genes_field = p[1]
    rows.append(dict(blk=blk, genes=genes_field, nloci=nloci, maxcol=maxcol,
                     mmid=mmid, core=core, sedef=sedef, verdict=verdict))

print(f"parsed {len(rows)} blocks; verdicts={Counter(r['verdict'] for r in rows)}")

def type_of(genes_field, verdict, maxcol):
    genes = [g for g in genes_field.replace("...", "").split(";") if g]
    bts = [biotype.get(g) for g in genes]
    if verdict == "cardinality": return "CARDINALITY_ARRAY(truth-artifact)"
    if verdict == "ep_oversplit": return "EP_OVERSPLIT_realdup(truth-artifact)"
    if verdict == "no_multiexon_edge": return "NO_MULTIEXON(unclassed)"
    # genuine_overmerge -> subtype
    if any(b in PSEUDO for b in bts): return "genuine:RETROCOPY_PSEUDOGENE"
    if any(b in NONCODING for b in bts): return "genuine:NONCODING_ELEMENT"
    if all(b == "protein_coding" for b in bts if b): return "genuine:SINGLE_DOMAIN_SHARE"
    return "genuine:UNANNOTATED_OTHER"

typc = Counter()
gen_types = Counter()
print("\n== 32 genuine over-merge blocks, TYPE + gate ==")
for r in rows:
    r["type"] = type_of(r["genes"], r["verdict"], r["maxcol"])
    typc[r["type"]] += 1
    if r["verdict"] == "genuine_overmerge":
        gen_types[r["type"]] += 1

for t, n in typc.most_common():
    print(f"  {t:<40} {n}")

# gate arithmetic
# naive E_p-purity gate: cuts ALL 57 (all are E_p-impure)
# refined gate: cut iff verdict==genuine_overmerge (colinear<=1 AND not cardinality/oversplit)
# proxy operational rule: cut iff maxcol<=1 AND verdict!=cardinality
naive_cut_genuine = sum(1 for r in rows if r["verdict"] == "genuine_overmerge")
naive_cut_artifact = sum(1 for r in rows if r["verdict"] in ("ep_oversplit","cardinality"))
# operational refined: maxcol<=1 and not cardinality-array (verdict not cardinality)
def refined_cut(r):
    return r["maxcol"] <= 1 and r["verdict"] != "cardinality"
ref_cut = [r for r in rows if refined_cut(r)]
ref_cut_genuine = sum(1 for r in ref_cut if r["verdict"] == "genuine_overmerge")
ref_cut_artifact = sum(1 for r in ref_cut if r["verdict"] in ("ep_oversplit","cardinality"))
ref_miss_genuine = sum(1 for r in rows if r["verdict"]=="genuine_overmerge" and not refined_cut(r))

print("\n== GATE ARITHMETIC (family/block level, denom=57 E_p-impure residual) ==")
print(f"  truth: 32 genuine over-merge (real FP) | 23 truth-artifact (9 oversplit + 14 cardinality) | 2 no-multiexon")
print(f"  NAIVE E_p-purity gate (cut all E_p-impure): cuts 32 genuine + 23 truth-artifact = 55  -> WRONGLY cuts 23 real dups")
print(f"  REFINED gate [maxcol<=1 AND not-cardinality]:")
print(f"     cuts genuine        : {ref_cut_genuine}/32")
print(f"     wrongly cuts artifact: {ref_cut_artifact}/23")
print(f"     misses genuine      : {ref_miss_genuine}/32")
# show any leakage
leak = [(r['blk'], r['genes'][:30], r['verdict'], r['maxcol']) for r in rows
        if (refined_cut(r) and r['verdict'] in ('ep_oversplit','cardinality'))
        or (r['verdict']=='genuine_overmerge' and not refined_cut(r))]
print(f"  leakage rows (misclassified by maxcol<=1 proxy): {leak}")
