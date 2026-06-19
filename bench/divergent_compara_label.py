#!/usr/bin/env python3
"""Build the GOLD divergent-paralog calibration set: take low-core (<0.20) both-named different-gene de-novo
locus pairs (the APOBEC3D/F pattern), fetch Ensembl Compara paralogy for their genes (independent truth,
resumable cache), and label each pair paralog/not. Output = divergent_calibration_pairs.tsv + the positive
rate. Reuses compara_fetch.py's REST/cache helpers (no re-fetch of already-cached genes)."""
import sys, csv, re, bisect, collections, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compara_fetch as cf

ANNOT = "/home/juanfra/winloci_scratch/annot_intervals.tsv"
META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
EDGES = "denovo_family_edges.tsv"
CORE_MAX = 0.20  # the divergent tail

# --- locus -> gene (max-overlap) ---
genes = collections.defaultdict(list)
for r in csv.DictReader(open(ANNOT), delimiter="\t"):
    genes[r["chrom"]].append((int(r["start"]), int(r["end"]), r["gene"]))
for c in genes:
    genes[c].sort()
starts = {c: [g[0] for g in v] for c, v in genes.items()}
def map_locus(ch, s, e):
    if ch not in genes:
        return None
    arr, st = genes[ch], starts[ch]
    i = bisect.bisect_right(st, e) - 1
    best, bov = None, 0
    while i >= 0 and arr[i][1] >= s - 200_000:
        gs, ge, nm = arr[i]
        ov = min(e, ge) - max(s, gs)
        if ov > bov:
            bov, best = ov, nm
        if gs < s - 200_000:
            break
        i -= 1
    return best if bov > 0 else None

coord = {r["id"]: (r["chrom"], int(r["start"]), int(r["end"]))
         for r in csv.DictReader(open(META), delimiter="\t")}

pairs = []          # (gene_a, gene_b, core_recip)
gene_set = set()
for r in csv.DictReader(open(EDGES), delimiter="\t"):
    cr = float(r["core_recip"])
    if cr >= CORE_MAX or r["a"] not in coord or r["b"] not in coord:
        continue
    ga = map_locus(*coord[r["a"]])
    gb = map_locus(*coord[r["b"]])
    if not ga or not gb or ga.startswith("LOC") or gb.startswith("LOC") or ga == gb:
        continue
    pairs.append((ga, gb, cr))
    gene_set.add(ga)
    gene_set.add(gb)
gene_set = sorted(gene_set)
print("divergent (<%.2f) both-named diff-gene pairs: %d  distinct genes: %d" % (CORE_MAX, len(pairs), len(gene_set)))

# --- fetch ENSG + paralogues per gene (resumable; reuse compara_fetch cache) ---
cache = cf.load_cache()
info = {}
for i, sym in enumerate(gene_set, 1):
    lk_url = "%s/lookup/symbol/homo_sapiens/%s?content-type=application/json" % (cf.REST, sym)
    lk, did = cf.fetch_cached(cache, "lookup", sym, lk_url)
    if did:
        time.sleep(cf.SPACING)
    ensg = None
    if lk.get("status") == "ok" and isinstance(lk.get("data"), dict):
        ensg = lk["data"].get("id")
    para = set()
    if lk.get("status") != "notfound":
        pa_url = ("%s/homology/symbol/homo_sapiens/%s?type=paralogues;compara=vertebrates;"
                  "content-type=application/json;format=condensed" % (cf.REST, sym))
        pa, did = cf.fetch_cached(cache, "homology_paralogues", sym, pa_url)
        if did:
            time.sleep(cf.SPACING)
        if pa.get("status") == "ok":
            data = pa.get("data") or {}
            d = data.get("data") if isinstance(data, dict) else None
            homs = d[0].get("homologies", []) if isinstance(d, list) and d else []
            para = {h.get("id") for h in homs if h.get("id")}
    info[sym] = {"ensg": ensg, "para": para}
    if i % 50 == 0:
        cf.save_cache(cache)
        print("  fetched %d/%d" % (i, len(gene_set)))
cf.save_cache(cache)

# --- label pairs (paralog iff one ENSG is in the other's paralogue set) ---
labelable = pos = 0
gp = collections.defaultdict(lambda: [0, 0])  # gene-pair -> [labelable, positive] (dedup view)
rows = []
for ga, gb, cr in pairs:
    ia, ib = info.get(ga), info.get(gb)
    if not ia or not ib or not ia["ensg"] or not ib["ensg"]:
        continue
    labelable += 1
    is_para = (ib["ensg"] in ia["para"]) or (ia["ensg"] in ib["para"])
    pos += is_para
    key = tuple(sorted((ga, gb)))
    gp[key][0] += 1
    gp[key][1] += int(is_para)
    rows.append((ga, gb, "%.3f" % cr, "paralog" if is_para else "not"))

gpos = sum(1 for k, v in gp.items() if v[1] > 0)
print("\n=== GOLD divergent calibration set ===")
print("locus pairs: %d | gold-labelable: %d | Compara-paralog: %d (%.0f%%)"
      % (len(pairs), labelable, pos, 100 * pos / labelable if labelable else 0))
print("distinct gene-pairs: %d | labelable: %d | with a Compara-paralog call: %d (%.0f%%)"
      % (len(gp), len(gp), gpos, 100 * gpos / len(gp) if gp else 0))
with open("divergent_calibration_pairs.tsv", "w") as f:
    f.write("gene_a\tgene_b\tcore_recip\tcompara\n")
    for r in rows:
        f.write("\t".join(r) + "\n")
print("wrote divergent_calibration_pairs.tsv")
