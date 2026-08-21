#!/usr/bin/env python3
"""Score the shipped HSA family catalog against ENSEMBL COMPARA paralogy — truth the annotation we
seed with did not generate.

WHY. The standing objection is that seeding O1 with an annotation makes any discovery circular. That
objection is precise and answerable, but it is about the EVALUATION, not the method: using annotation
as a seed is fine; scoring against truth DERIVED FROM THAT ANNOTATION is not. Compara paralogy is
independent of the CHM13 GFF used to seed and to name nodes — it comes from Ensembl's gene trees.

⚠ THE ANNOTATION IS STILL USED HERE, to map a node to a gene symbol. What it is NOT used for is
deciding whether two genes are paralogues. That is the whole distinction.

TWO DIRECTIONS, reported separately because they are different assertions:
  PRECISION  of same-family pairs whose two genes are both known to Compara, how many does Compara
             call paralogues?
  RECALL     of Compara paralogue pairs whose two genes both have a catalog node, how many does the
             catalog put in one family?

⚠ Compara paralogy and "same gene family" are not the same predicate — Compara includes ancient
paralogues that no RNA homology rule should join, and excludes recent copies it has not resolved. So
this BOUNDS agreement; it does not score correctness. Read the numbers as concordance.

Human/CHM13 only. Cache-first; a small capped number of live lookups is allowed for missing symbols.
"""
import csv, gzip, json, os, re, sys, time, urllib.request, collections

REPO  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D     = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
GFF   = "/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0_RefSeq_full.gff.gz"
CACHE = os.path.join(REPO, "bench", "compara_cache.json")
REST  = "https://rest.ensembl.org"
MAX_FETCH = int(os.environ.get("O1_COMPARA_MAX_FETCH", "0"))   # 0 = cache-only

cache = json.load(open(CACHE))
fetched = 0


def rest(key, url):
    """cache-first; only fetch when explicitly budgeted"""
    global fetched
    if key in cache:
        return cache[key].get("data")
    if fetched >= MAX_FETCH:
        return None
    try:
        req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
        d = json.loads(urllib.request.urlopen(req, timeout=20).read())
    except Exception:
        d = None
    cache[key] = {"data": d}
    fetched += 1
    time.sleep(0.12)
    return d


def ensg(sym):
    d = rest(f"lookup|{sym}", f"{REST}/lookup/symbol/homo_sapiens/{sym}?content-type=application/json")
    return d.get("id") if isinstance(d, dict) else None


def paralogues(sym):
    d = rest(f"homology_paralogues|{sym}",
             f"{REST}/homology/symbol/homo_sapiens/{sym}?type=paralogues;format=condensed;content-type=application/json")
    out = set()
    try:
        for blk in d["data"]:
            for h in blk.get("homologies", []):
                if h.get("id"):
                    out.add(h["id"])
    except Exception:
        return None
    return out


def load_genes():
    """protein-coding gene intervals from the APPROVED CHM13 GFF (never HSA_genomic.gff)"""
    g = collections.defaultdict(list)
    with gzip.open(GFF, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene" or "gene_biotype=protein_coding" not in f[8]:
                continue
            m = re.search(r"(?:^|;)gene=([^;]+)", f[8])
            if m:
                g[f[0]].append((int(f[3]) - 1, int(f[4]), m.group(1)))
    for c in g:
        g[c].sort()
    return g


def main():
    genes = load_genes()
    rows = list(csv.DictReader(open(f"{D}/HSA_gwcat.copies.tsv"), delimiter="\t"))
    fam = collections.defaultdict(list)
    sym = {}
    for r in rows:
        k = f"{r['family_id']}~{r['copy_idx']}"
        fam[r["family_id"]].append(k)
        c, s, e = r["chrom"], int(r["start"]), int(r["end"])
        best, bov = None, 0
        for g0, g1, name in genes.get(c, []):
            if g0 >= e:
                break
            ov = min(e, g1) - max(s, g0)
            if ov > bov:
                best, bov = name, ov
        if best:
            sym[k] = best
    print(f"HSA catalog: {len(fam)} families, {len(rows)} copies; mapped to a protein-coding gene: "
          f"{len(sym)}/{len(rows)} = {len(sym)/len(rows):.4f}", flush=True)

    # ---------- PRECISION ----------
    checkable = confirmed = 0
    same_gene = 0
    misses = []
    for f_, members in fam.items():
        named = [(m, sym[m]) for m in members if m in sym]
        for i in range(len(named)):
            for j in range(i + 1, len(named)):
                a, sa = named[i]; b, sb = named[j]
                if sa == sb:
                    same_gene += 1
                    continue
                pa, pb = paralogues(sa), paralogues(sb)
                ea, eb = ensg(sa), ensg(sb)
                if pa is None or pb is None or ea is None or eb is None:
                    continue
                checkable += 1
                if (eb in pa) or (ea in pb):
                    confirmed += 1
                else:
                    misses.append((f_, sa, sb))
    print(f"\n=== PRECISION vs INDEPENDENT Compara paralogy (unit = pair) ===")
    print(f"  same-family pairs mapping to the SAME gene (a split locus, not a paralogy claim): {same_gene}")
    print(f"  checkable pairs (both genes resolved in Compara):  {checkable}")
    if checkable:
        print(f"  Compara CONFIRMS paralogy:  {confirmed}/{checkable} = {confirmed/checkable:.4f}")
        print(f"  not confirmed: {len(misses)}   e.g. {misses[:6]}")

    # ---------- RECALL ----------
    bysym = collections.defaultdict(set)
    for k, s in sym.items():
        bysym[s].add(k.split("~")[0])
    tested = joined = 0
    for s in list(bysym):
        ps = paralogues(s)
        e_of = {}
        if not ps:
            continue
        for s2 in bysym:
            if s2 == s or s2 <= s:
                continue
            e2 = ensg(s2)
            if e2 and e2 in ps:
                tested += 1
                if bysym[s] & bysym[s2]:
                    joined += 1
    print(f"\n=== RECALL: Compara paralogue pairs whose BOTH genes have a catalog node ===")
    print(f"  testable pairs: {tested}")
    if tested:
        print(f"  catalog puts them in ONE family: {joined}/{tested} = {joined/tested:.4f}")
    print(f"\n  (live lookups used: {fetched}; cache entries: {len(cache)})")
    if fetched:
        json.dump(cache, open(CACHE, "w"))


if __name__ == "__main__":
    main()
