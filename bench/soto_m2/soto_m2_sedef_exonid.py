#!/usr/bin/env python3
"""Pre-registered KEY=sedefexon (docs/archive/2026-09/PREREG_soto_sedef_exon_identity_2026-09-30.md): SEDEF-projected exon links kept only when
the exon pair itself is >= 98% identical.

Re-derives the native CHM13 v1.0 SEDEF links with soto_replication.py's own CIGAR walk (rows >= 0.98 identity, exons covered
>= 0.99 by the row's span, projected through the CIGAR), and for every projection measures the exon pair's identity:
1 - global edit distance (edlib) / max length, the projected side reverse-complemented when the row's second side is '-'.
Arms: A map-back alone; B map-back + every SEDEF link (must equal the frozen shared_exons_2334_finalv1_native.tsv); C map-back +
the SEDEF links with a projection at exon identity >= 0.98. Each arm classified exactly as in soto_m2_families.py.

    python3 bench/soto_m2/soto_m2_sedef_exonid.py --sedef final_v1_clean.bed --genome t2t-chm13-v1.0.fa.gz --cat cat_v4.bed \
        --geneset elig.tsv --full-geneset full.tsv --famcn-ours famcn_ours_allwssd.tsv --out-links links.tsv
Needs edlib and pysam (the miniforge python)."""
import argparse
import csv
import os
import sys
from collections import Counter, defaultdict

import edlib
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "soto"))
sys.path.insert(0, HERE)
import soto_replication as sr  # noqa: E402
from soto_m2_families import COMBOS, classify, is_pseudo  # noqa: E402

SOTO = os.path.join(HERE, "..", "soto")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def exon_identity(fa, c_e, e_lo, e_hi, c_p, p_lo, p_hi, flip, cache):
    key = (c_e, e_lo, e_hi, c_p, p_lo, p_hi, flip)
    if key not in cache:
        E = fa.fetch(c_e, e_lo, e_hi).upper()
        P = fa.fetch(c_p, p_lo, p_hi).upper()
        if flip:
            P = revcomp(P)
        L = max(len(E), len(P))
        cache[key] = 1.0 - edlib.align(E, P, mode="NW", task="distance")["editDistance"] / L if L else 0.0
    return cache[key]


def sedef_links(sedef_path, genome, cat_bed, geneset_path, min_identity=0.98, min_cov=0.99):
    """(g1, g2) sorted pair -> best exon-pair identity over the projections that create the link."""
    genes, _ = sr.load_geneset(geneset_path)
    exons, _ = sr.load_exons(cat_bed, genes)
    per_chrom = defaultdict(list)
    for g, evs in exons.items():
        for c, s, e in evs:
            per_chrom[c].append((s, e, g))
    for c in per_chrom:
        per_chrom[c].sort()
    fa = pysam.FastaFile(genome)
    best, cache = {}, {}
    n_rows = n_used = n_proj = 0

    def link(own_chrom, p_lo, p_hi, g_other, ident):
        for s2, e2, g_own in per_chrom.get(own_chrom, ()):
            if e2 <= p_lo:
                continue
            if s2 >= p_hi:
                break
            if g_own != g_other:
                k = tuple(sorted((g_own, g_other)))
                best[k] = max(best.get(k, -1.0), ident)

    for line in open(sedef_path):
        n_rows += 1
        f = line.rstrip("\n").split("\t")
        if len(f) < 33:
            continue
        try:
            identity = float(f[20])
        except ValueError:
            continue
        if identity < min_identity:
            continue
        c1, s1, e1, c2, s2, e2, strand2 = f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5]), f[9]
        blocks, p1_end, bw_end = sr.build_blocks(s1, sr.cigar_ops(f[32]))
        if p1_end != e1 or bw_end != (e2 - s2):
            continue
        n_used += 1
        flip = strand2 == "-"
        # other = side2 (bwalk axis), own = side1
        for s, e, g_other in per_chrom.get(c2, ()):
            if e <= s2:
                continue
            if s >= e2:
                break
            if (min(e, e2) - max(s, s2)) / (e - s) < min_cov:
                continue
            v_s, v_e = max(s, s2), min(e, e2)
            bw_lo, bw_hi = sorted((sr.genomic_to_bwalk(v_s, s2, e2, strand2), sr.genomic_to_bwalk(v_e, s2, e2, strand2)))
            proj = sr.project_bwalk_to_a(blocks, bw_lo, bw_hi)
            if proj is None:
                continue
            n_proj += 1
            ident = exon_identity(fa, c2, v_s, v_e, c1, proj[0], proj[1], flip, cache)
            link(c1, proj[0], proj[1], g_other, ident)
        # other = side1 (plain axis), own = side2
        for s, e, g_other in per_chrom.get(c1, ()):
            if e <= s1:
                continue
            if s >= e1:
                break
            if (min(e, e1) - max(s, s1)) / (e - s) < min_cov:
                continue
            v_s, v_e = max(s, s1), min(e, e1)
            proj = sr.project_a_to_bwalk(blocks, v_s, v_e)
            if proj is None:
                continue
            g_lo = sr.bwalk_to_genomic(proj[0], s2, e2, strand2)
            g_hi = sr.bwalk_to_genomic(proj[1], s2, e2, strand2)
            p_lo, p_hi = min(g_lo, g_hi), max(g_lo, g_hi)
            n_proj += 1
            ident = exon_identity(fa, c1, v_s, v_e, c2, p_lo, p_hi, flip, cache)
            link(c2, p_lo, p_hi, g_other, ident)
    print(f"[sedef] {n_rows} rows, {n_used} used (identity >= {min_identity}, CIGAR reconciles); {n_proj} exon projections; "
          f"{len(best)} links", file=sys.stderr)
    return best


def dedup(edges):
    seen, out = set(), []
    for x, y in edges:
        k = tuple(sorted((x, y)))
        if x != y and k not in seen:
            seen.add(k)
            out.append((x, y))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", default=os.path.join(SOTO, "soto_famCN_S1C.tsv"))
    ap.add_argument("--shared", default=os.path.join(SOTO, "shared_exons_5154_exon_mapback.tsv"))
    ap.add_argument("--frozen-native", default=os.path.join(SOTO, "shared_exons_2334_finalv1_native.tsv"))
    ap.add_argument("--split", default=os.path.join(SOTO, "soto_split_2026-09-29.tsv"))
    for k in ("sedef", "genome", "cat", "geneset", "full_geneset", "famcn_ours", "out_links"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    ap.add_argument("--exon-identity", type=float, default=0.98)
    a = ap.parse_args(argv)

    best = sedef_links(a.sedef, a.genome, a.cat, a.full_geneset)
    frozen = {tuple(sorted(e)) for e in sr.read_edges(a.frozen_native)}
    if set(best) != frozen:
        sys.exit(f"CORRECTNESS GATE FAILED: re-derived {len(best)} links vs frozen {len(frozen)} "
                 f"(only re-derived {len(set(best) - frozen)}, only frozen {len(frozen - set(best))})")
    print(f"correctness gate: the {len(best)} re-derived SEDEF links equal the frozen native-v1 edge file")

    mapback = sr.read_edges(a.shared)
    mb_keys = {tuple(sorted(e)) for e in mapback}
    with open(a.out_links, "w") as out:
        out.write("gene_a\tgene_b\tbest_exon_identity\tin_mapback\n")
        for (x, y), v in sorted(best.items()):
            out.write(f"{x}\t{y}\t{v:.4f}\t{int((x, y) in mb_keys)}\n")
    ids = sorted(best.values())
    only = sorted(v for k, v in best.items() if k not in mb_keys)
    q = lambda xs, p: xs[min(len(xs) - 1, int(p * len(xs)))] if xs else float("nan")
    keepC = [k for k, v in best.items() if v >= a.exon_identity]
    print(f"SEDEF links: {len(ids)}; exon identity quartiles {q(ids, .25):.3f} / {q(ids, .5):.3f} / {q(ids, .75):.3f}; "
          f">= {a.exon_identity}: {len(keepC)} ({len(keepC) / len(ids):.3f})")
    print(f"SEDEF-only links (not in map-back): {len(only)}; quartiles {q(only, .25):.3f} / {q(only, .5):.3f} / {q(only, .75):.3f}; "
          f">= {a.exon_identity}: {sum(v >= a.exon_identity for v in only)}\n")

    genes, _ = sr.load_geneset(a.geneset)
    full, _ = sr.load_geneset(a.full_geneset)
    clean, _amb = sr.load_truth(a.truth)
    info, fam_all, manual = {}, defaultdict(set), set()
    for r in csv.DictReader(open(a.truth), delimiter="\t"):
        g = r["Gene ID"]
        info.setdefault(g, dict(name=r["Gene Name"], biotype=r["Biotype"], fams=[],
                                backbone=r.get("In Table S1 (SD98 gene set)") == "Yes"))
        if r["Family ID"]:
            info[g]["fams"].append(r["Family ID"])
            if not r["Family ID"].startswith("Unassigned"):
                fam_all[r["Family ID"]].add(g)
        if r.get("Family MAD", "").strip().lower() == "manual merge":
            manual.add(r["Family ID"])
    cn_s1c = {g: v for g, v in sr.load_famcn(a.truth, "Median famCN").items() if g in genes}
    cn_ours = {g: v for g, v in sr.load_famcn(a.famcn_ours, "famCN_sotoiv").items() if g in genes}
    half = sr.load_split(a.split)
    every = set(info) | genes | full
    sedef_all = list(best)
    arms = {"A": mapback, "B": dedup(mapback + sedef_all), "C": dedup(mapback + keepC)}

    print("| copy numbers | filter | A exact (dev / held-out) | B exact | C exact (dev / held-out) | ARI A / B / C "
          "| C recovered vs A | C broken vs A |")
    print("|---|---|---|---|---|---|---|---|")
    verdict, details = None, []
    for src, cn_use in (("s1c", cn_s1c), ("ours", cn_ours)):
        for key0, label0, no_p, no_l in COMBOS:
            keep = {g for g in every if not (no_p and is_pseudo(info.get(g, {}).get("biotype", "")))
                    and not (no_l and info.get(g, {}).get("biotype", "") == "lncRNA")}
            res = {}
            for arm, edges in arms.items():
                rows, anc, _ps, _po = classify(keep, genes, full, edges, sedef_all, clean, fam_all, info, manual, cn_s1c,
                                               cn_use, filtered=bool(key0))
                cls = {r["family"]: r for r in rows}
                fam_half = {r["family"]: Counter(half.get(g) for g in r["clean"]).most_common(1)[0][0]
                            for r in rows if r["clean"]}
                ex = Counter(fam_half.get(f) for f, r in cls.items() if r["cls"] == "match")
                res[arm] = (cls, anc, ex, sum(r["cls"] == "match" for r in rows))
            A, B, C = res["A"], res["B"], res["C"]
            rec = sorted((f for f in C[0] if C[0][f]["cls"] == "match" and A[0][f]["cls"] != "match"), key=lambda x: int(x[3:]))
            brk = sorted((f for f in C[0] if C[0][f]["cls"] != "match" and A[0][f]["cls"] == "match"), key=lambda x: int(x[3:]))
            print(f"| {src} | {label0} | {A[3]} ({A[2]['dev']} / {A[2]['heldout']}) | {B[3]} | {C[3]} ({C[2]['dev']} / "
                  f"{C[2]['heldout']}) | {A[1]['ari']:.4f} / {B[1]['ari']:.4f} / {C[1]['ari']:.4f} | {len(rec)} | {len(brk)} |")
            details.append((src, label0, rec, brk, A[0], C[0]))
            if src == "s1c" and key0 == "":
                verdict = "CONVERGES" if C[3] > (A[3] + B[3]) / 2 else "DOES NOT CONVERGE"
    print(f"\nVERDICT (S1C copy numbers, all genes): {verdict}\n")
    for src, label0, rec, brk, A, C in details:
        if label0 != "all genes":
            continue
        r_txt = ", ".join(f"{f} (was {A[f]['cls']})" for f in rec) or "-"
        b_txt = ", ".join(f"{f} (+{len(C[f]['extra'])})" for f in brk) or "-"
        print(f"- {src}, {label0}: recovered {r_txt}; broken {b_txt}")

if __name__ == "__main__":
    main()
