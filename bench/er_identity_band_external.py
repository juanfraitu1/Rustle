#!/usr/bin/env python3
"""External corroboration of E_r edges, STRATIFIED BY IDENTITY BAND.

⭐ **The gap this closes.** The readiness audit (§6ax) records that the catalog's evidence covers
the wrong end of its own distribution: median edge identity is **0.8287** and **86.31% of edges sit
below 0.90**, while every external check to date (NPIP annotated median 0.9779, GOLGA6L7 0.9673)
lives at **>= 0.90**. Nothing external speaks to the ~0.83 band that is most of the catalog. This
instrument asks whether the gorilla-native RefSeq annotation does.

**Truth source.** `GGO_genomic.gff` (GCF_029281585.2-RS_2024_02, NHGRI_mGorGor1-v2.0_pri, taxon
9595) -- **gorilla-native**, human-curated relation, and produced with no knowledge of this
pipeline. Every `gene` and `pseudogene` feature carries a `description=` attribute.

⚠ **READ THE DESCRIPTION, NOT THE SYMBOL.** In this GFF RFPL is 9/9 and GOLGA6 14/14 LOC-named, so
a symbol grep sees essentially none of them -- the trap that produced the retracted §6b answer key.
The `description=` field survives LOC-naming ("beta-defensin 125-like" on a `LOC129533586` gene),
which is what makes it usable here.

⚠ **"uncharacterized LOCxxxxx" is NOT an annotation.** Such endpoints are counted UNINFORMATIVE,
never as disagreement -- never infer content from an absent name.

⚠ **Two endpoints hitting the SAME gene feature are excluded**, not counted as agreement: that is
a self-overlap, the W063 pathology, not external corroboration.

⚠ **The null matches the SIZE distribution, not the edge count.** An edge-count-matched null proves
nothing (register). Null pairs are drawn from the catalog's own nodes, matched on the exon-sum
length decile of BOTH endpoints and on same-chromosome status, since co-location inflates
description agreement on its own.

Usage:
  er_identity_band_external.py <gff> <nodes.tsv> <edges.tsv> [--null N] [--seed S]
"""
import sys
import re
import gzip
import random
import bisect
import collections

BANDS = [(0.60, 0.70), (0.70, 0.80), (0.80, 0.90), (0.90, 1.001)]
UNINFORMATIVE = re.compile(r"^uncharacterized\b", re.I)


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def norm_desc(d):
    """Normalise a RefSeq description to a family-level stem.

    Strips the curation prefix, isoform suffixes, the '-like' paralog marker and trailing
    member numerals, so 'beta-defensin 125-like' and 'beta-defensin 124' share a stem.
    """
    d = d.strip().lower()
    d = re.sub(r"^low quality protein:\s*", "", d)
    d = re.sub(r"\s+isoform\s+x?\d+$", "", d)
    d = d.replace("%2c", ",").replace("%3b", ";").replace("%25", "%")
    d = re.sub(r"[-\s]like$", "", d)
    d = re.sub(r"[\s,]+\d+$", "", d)          # trailing member number
    d = re.sub(r"\s+(alpha|beta|gamma|delta)$", "", d)
    d = re.sub(r"\s+", " ", d).strip()
    return d


def load_gff(path):
    """chrom -> (sorted starts, [(start, end, gene_id, biotype, stem)])"""
    per = collections.defaultdict(list)
    kept = skipped = 0
    for line in opener(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
            continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        desc = attrs.get("description", "")
        gid = attrs.get("ID", attrs.get("Name", ""))
        bio = attrs.get("gene_biotype", f[2])
        if not desc or UNINFORMATIVE.match(desc):
            skipped += 1
            continue
        stem = norm_desc(desc)
        if not stem:
            skipped += 1
            continue
        kept += 1
        per[f[0]].append((int(f[3]) - 1, int(f[4]), gid, bio, stem))
    idx = {}
    for c, v in per.items():
        v.sort()
        idx[c] = ([x[0] for x in v], v)
    sys.stderr.write(f"# GFF: {kept} informative gene/pseudogene features, "
                     f"{skipped} uninformative/undescribed, {len(idx)} contigs\n")
    return idx


def hits(idx, chrom, s, e):
    """Features overlapping [s,e), best-overlap first."""
    if chrom not in idx:
        return []
    starts, v = idx[chrom]
    i = bisect.bisect_left(starts, s)
    out = []
    j = max(0, i - 400)                       # features are <= a few Mb; walk back a window
    while j < len(v) and v[j][0] < e:
        fs, fe = v[j][0], v[j][1]
        ov = min(fe, e) - max(fs, s)
        if ov > 0:
            out.append((ov, v[j]))
        j += 1
    out.sort(key=lambda x: -x[0])
    return out


RECIP = 0.50


def annotate(idx, chrom, s, e):
    """-> (gene_id, stem) of the RECIPROCALLY-overlapping informative feature, or None.

    ⚠ **Reciprocal, not best-overlap-by-bp.** A first version took the feature with the largest
    absolute overlap; that labels any node lying inside a huge gene as that gene, and the symptom
    was unmistakable -- "titin" appeared on 26 endpoints and "myosin-11" on 27, i.e. small nodes
    inside long introns inheriting the host gene's name. The pair must cover >= RECIP of EACH
    OTHER, so the node has to be gene-sized and coincident, not merely contained. Endpoints with
    no reciprocal match are UNINFORMATIVE -- we do not know what gene they are -- never a
    disagreement.
    """
    best = None
    for ov, feat in hits(idx, chrom, s, e):
        fs, fe = feat[0], feat[1]
        f_node = ov / max(1, e - s)
        f_feat = ov / max(1, fe - fs)
        if f_node >= RECIP and f_feat >= RECIP:
            score = min(f_node, f_feat)
            if best is None or score > best[0]:
                best = (score, feat[2], feat[4])
    return (best[1], best[2]) if best else None


def verdict(a, b):
    """None = uninformative, 'same_gene' = excluded, True/False = corroborated or not."""
    if a is None or b is None:
        return None
    if a[0] == b[0]:
        return "same_gene"
    return a[1] == b[1]


def main():
    gff, nodes_f, edges_f = sys.argv[1], sys.argv[2], sys.argv[3]
    nnull = int(sys.argv[sys.argv.index("--null") + 1]) if "--null" in sys.argv else 20000
    seed = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 101
    rng = random.Random(seed)

    idx = load_gff(gff)

    nodes = []
    hdr = None
    for line in open(nodes_f):
        f = line.rstrip("\n").split("\t")
        if hdr is None:
            hdr = {n: i for i, n in enumerate(f)}
            continue
        nodes.append((f[hdr["chrom"]], int(f[hdr["start"]]), int(f[hdr["end"]]),
                      int(f[hdr["exon_sum_len"]])))
    ann = [annotate(idx, c, s, e) for c, s, e, _ in nodes]
    lens = sorted(n[3] for n in nodes)
    dec = lambda L: bisect.bisect_left(lens, L) * 10 // max(1, len(lens))

    # ---- real edges, stratified by co-location
    per_band = collections.defaultdict(lambda: collections.Counter())
    per_strat = collections.defaultdict(lambda: collections.Counter())
    hdr = None
    for line in open(edges_f):
        f = line.rstrip("\n").split("\t")
        if hdr is None:
            hdr = {n: i for i, n in enumerate(f)}
            continue
        ident = float(f[hdr["identity"]])
        a = annotate(idx, f[hdr["chrom_i"]], int(f[hdr["start_i"]]), int(f[hdr["end_i"]]))
        b = annotate(idx, f[hdr["chrom_j"]], int(f[hdr["start_j"]]), int(f[hdr["end_j"]]))
        v = verdict(a, b)
        same = f[hdr["chrom_i"]] == f[hdr["chrom_j"]]
        sep = (abs(int(f[hdr["start_i"]]) - int(f[hdr["start_j"]])) if same else -1)
        strat = "cross_chrom" if not same else ("same_chrom_far" if sep > 1_000_000
                                                else "same_chrom_near")
        for key in (strat, f"{strat}|{[b_ for b_ in BANDS if b_[0] <= ident < b_[1]]}"):
            pass
        c2 = per_strat[strat]
        c2["n"] += 1
        if v is None:
            c2["uninformative"] += 1
        elif v == "same_gene":
            c2["same_gene"] += 1
        else:
            c2["informative"] += 1
            c2["agree"] += 1 if v else 0
            for lo, hi in BANDS:
                if lo <= ident < hi:
                    k = per_strat[f"{strat} [{lo:.2f},{hi:.2f})"]
                    k["informative"] += 1
                    k["agree"] += 1 if v else 0
                    break
        for lo, hi in BANDS:
            if lo <= ident < hi:
                c = per_band[(lo, hi)]
                c["n"] += 1
                if v is None:
                    c["uninformative"] += 1
                elif v == "same_gene":
                    c["same_gene"] += 1
                else:
                    c["informative"] += 1
                    c["agree"] += 1 if v else 0
                break

    # ---- size-matched null, drawn from the catalog's own nodes
    null = collections.Counter()
    pool_by = collections.defaultdict(list)
    for i, (c, s, e, L) in enumerate(nodes):
        pool_by[dec(L)].append(i)
    tries = 0
    nulls = {k: collections.Counter() for k in
             ("cross_chrom", "same_chrom_far", "same_chrom_near")}
    by_chrom = collections.defaultdict(list)
    for i, (c, s_, e_, L) in enumerate(nodes):
        by_chrom[c].append(i)
    while (null["informative"] < nnull or
           min(v["informative"] for v in nulls.values()) < nnull // 10) and tries < nnull * 200:
        tries += 1
        d1, d2 = rng.randrange(10), rng.randrange(10)
        if not pool_by[d1] or not pool_by[d2]:
            continue
        i, j = rng.choice(pool_by[d1]), rng.choice(pool_by[d2])
        if i == j:
            continue
        v = verdict(ann[i], ann[j])
        if v is None or v == "same_gene":
            continue
        null["informative"] += 1
        null["agree"] += 1 if v else 0
        same = nodes[i][0] == nodes[j][0]
        sep = abs(nodes[i][1] - nodes[j][1]) if same else -1
        k = "cross_chrom" if not same else ("same_chrom_far" if sep > 1_000_000
                                            else "same_chrom_near")
        nulls[k]["informative"] += 1
        nulls[k]["agree"] += 1 if v else 0

    print("\n=== EXTERNAL CORROBORATION BY IDENTITY BAND ===")
    print("band\tedges\tuninformative\tsame_gene\tinformative\tagree\trate")
    tot = collections.Counter()
    for lo, hi in BANDS:
        c = per_band[(lo, hi)]
        r = c["agree"] / c["informative"] if c["informative"] else float("nan")
        print(f"[{lo:.2f},{hi:.2f})\t{c['n']}\t{c['uninformative']}\t{c['same_gene']}\t"
              f"{c['informative']}\t{c['agree']}\t{r:.4f}")
        tot.update(c)
    rt = tot["agree"] / tot["informative"] if tot["informative"] else float("nan")
    print(f"ALL\t{tot['n']}\t{tot['uninformative']}\t{tot['same_gene']}\t"
          f"{tot['informative']}\t{tot['agree']}\t{rt:.4f}")
    nr = null["agree"] / null["informative"] if null["informative"] else float("nan")
    print(f"\nSIZE-MATCHED NULL: {null['agree']}/{null['informative']} = {nr:.4f} "
          f"(seed {seed}, {tries} draws)")
    for lo, hi in BANDS:
        c = per_band[(lo, hi)]
        if c["informative"] and nr:
            r = c["agree"] / c["informative"]
            print(f"  enrichment [{lo:.2f},{hi:.2f}) = {r/nr:.2f}x  (vs pooled null)")

    print("\n=== PROXIMITY CONTROL: each stratum against ITS OWN geometry-matched null ===")
    print("stratum\tedges\tinformative\tagree\trate\tnull_rate\tenrichment")
    for k in ("cross_chrom", "same_chrom_far", "same_chrom_near"):
        c, nk = per_strat[k], nulls[k]
        if not c["informative"] or not nk["informative"]:
            print(f"{k}\t{c['n']}\t{c['informative']}\t{c['agree']}\t-\t-\t-")
            continue
        r = c["agree"] / c["informative"]
        n0 = nk["agree"] / nk["informative"]
        enr = (r / n0) if n0 else float("inf")
        print(f"{k}\t{c['n']}\t{c['informative']}\t{c['agree']}\t{r:.4f}\t"
              f"{n0:.4f} (n={nk['informative']})\t{enr:.2f}x")
    print("\n=== the decisive cell: CROSS-CHROMOSOME edges, by identity band ===")
    nk = nulls["cross_chrom"]
    n0 = nk["agree"] / nk["informative"] if nk["informative"] else float("nan")
    for lo, hi in BANDS:
        c = per_strat[f"cross_chrom [{lo:.2f},{hi:.2f})"]
        if c["informative"]:
            r = c["agree"] / c["informative"]
            print(f"  [{lo:.2f},{hi:.2f})  {c['agree']}/{c['informative']} = {r:.4f}   "
                  f"enrichment {r/n0:.2f}x" if n0 else "")


if __name__ == "__main__":
    main()
