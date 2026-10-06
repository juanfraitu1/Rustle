#!/usr/bin/env python3
"""PSV ceiling: for every copy of every multi-copy family in a genome-wide catalog, the identity to its most similar sibling copy on
the SPLICED copy sequence (what a read sees), from one all-vs-all minimap2 asm20 run per species; plus, for gorilla, the outcome of the
real-read copy-assignment runs (rna_allele o2_*.assignments.tsv, 309 families) by identity band of the copy the read sits in.

    python3 bench/psv_ceiling/psv_ceiling.py          # writes /mnt/linuxdisk/tmp/psv_ceiling/{species}.nearest.tsv, o2_bands.tsv, data.json
Species are never pooled. Identity = matches / alignment block length of the best same-family hit whose aligned span covers >= 50 % of
the shorter copy; copies without such a hit are counted separately ("no alignment").
"""
import collections, csv, glob, json, os, statistics, subprocess
OUT = "/mnt/linuxdisk/tmp/psv_ceiling"
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
SPECIES = {
    "gorilla":  ("/mnt/linuxdisk/tmp/rna_allele/fibro_cat.copies.fa", "KB3781 fibroblast catalog (GWFAM, 667 families)"),
    "human":    ("/mnt/linuxdisk/tmp/rustle_figures/runs/human_testis/human_testis.cat.copies.fa", "testis legacy catalog"),
    "chimp":    ("/mnt/linuxdisk/tmp/rustle_figures/runs/chimp_PTR/chimp_PTR.cat.copies.fa", "chimpanzee legacy catalog"),
}
READ_LEN = 2162   # median gorilla KB3781 fibroblast read length (fibro.molecules.tsv)
BANDS = [("identical", 1.0, None), ("99.9-100%", 0.999, 1.0), ("99.5-99.9%", 0.995, 0.999), ("99-99.5%", 0.99, 0.995),
         ("98-99%", 0.98, 0.99), ("95-98%", 0.95, 0.98), ("<95%", None, 0.95)]


def band(x):
    if x is None: return "no alignment"
    if x >= 1.0: return "identical"
    for name, lo, hi in BANDS[1:]:
        if (lo is None or x >= lo) and (hi is None or x < hi): return name
    return "<95%"


def lengths(fa):
    L, cur = {}, None
    for ln in open(fa):
        if ln[0] == ">": cur = ln[1:].split()[0]; L[cur] = 0
        else: L[cur] += len(ln.strip())
    return L


def nearest(species, fa):
    paf = f"{OUT}/{species}.ava.paf"
    if not os.path.exists(paf):
        with open(paf + ".tmp", "w") as o, open(paf + ".log", "w") as e:
            subprocess.run([MM2, "-cx", "asm20", "-N", "100", "-p", "0.1", "-t", "2", fa, fa], stdout=o, stderr=e, check=True)
        os.replace(paf + ".tmp", paf)
    L = lengths(fa)
    best = {}
    for ln in open(paf):
        f = ln.split("\t"); q, t = f[0], f[5]
        if q == t or q.split("|")[0] != t.split("|")[0]: continue
        ql, tl = int(f[1]), int(f[6]); span = int(f[3]) - int(f[2])
        if span < 0.5 * min(ql, tl): continue
        ident = int(f[9]) / max(1, int(f[10]))
        if q not in best or ident > best[q][0]: best[q] = (ident, t, span, int(f[10]) - int(f[9]))
    rows = []
    for name in L:
        fam = name.split("|")[0]
        b = best.get(name)
        rows.append({"copy": name, "family": fam, "len": L[name], "nearest": b[1] if b else "", "identity": round(b[0], 6) if b else None,
                     "aligned_span": b[2] if b else 0, "diffs": b[3] if b else None, "band": band(b[0] if b else None)})
    with open(f"{OUT}/{species}.nearest.tsv", "w") as o:
        w = csv.DictWriter(o, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    return rows


def o2_bands(near):
    """Gorilla real reads: rna_allele o2_*.assignments.tsv (one row per contested read) x identity band of the copy the read sits in."""
    ident = {(r["family"], r["copy"].split("|")[1]): r["identity"] for r in near}
    cnt = collections.Counter(); fams = collections.defaultdict(set)
    for p in sorted(glob.glob("/mnt/linuxdisk/tmp/rna_allele/o2_*.assignments.tsv")):
        for x in csv.DictReader(open(p), delimiter="\t"):
            key = (x["family_id"], x["catalog_copy_idx"])
            b = band(ident.get(key)) if key in ident else "no alignment"
            if x["status"] == "assigned": cls = "assigned"
            elif x["origin_rejected"] == "1": cls = "matches no catalog copy"
            elif x["status"] == "tied" or x["n_decisive"] == "0": cls = "no distinguishing column"
            else: cls = "distinguishing columns, evidence insufficient"
            cnt[(b, cls)] += 1; fams[b].add(x["family_id"])
    rows = [{"band": b, "class": c, "n_reads": n, "n_families": len(fams[b])} for (b, c), n in sorted(cnt.items())]
    with open(f"{OUT}/o2_bands.tsv", "w") as o:
        w = csv.DictWriter(o, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    return rows


def main():
    data = {"read_len": READ_LEN, "species": {}, "bands": [b[0] for b in BANDS] + ["no alignment"]}
    for sp, (fa, label) in SPECIES.items():
        rows = nearest(sp, fa)
        ids = sorted(r["identity"] for r in rows if r["identity"] is not None)
        one_site = 1 - 1 / READ_LEN
        data["species"][sp] = {
            "label": label, "n_copies": len(rows), "n_families": len({r["family"] for r in rows}),
            "no_alignment": sum(r["identity"] is None for r in rows),
            "band_counts": collections.Counter(r["band"] for r in rows),
            "ge_one_site": sum(x >= one_site for x in ids), "ge_999": sum(x >= 0.999 for x in ids), "identical": sum(x >= 1.0 for x in ids),
            "median_identity": statistics.median(ids), "identities": ids,
        }
        print(sp, {k: v for k, v in data["species"][sp].items() if k != "identities"})
    data["o2"] = o2_bands(nearest("gorilla", SPECIES["gorilla"][0]))
    json.dump(data, open(f"{OUT}/data.json", "w"))
    print("o2 bands:", {(r["band"], r["class"]): r["n_reads"] for r in data["o2"]})


if __name__ == "__main__":
    main()
