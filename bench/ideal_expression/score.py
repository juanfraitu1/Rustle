#!/usr/bin/env python3
"""Score one ideal-expression run (docs/PREREG_ideal_expression_2026-10-06.md): can the default pipeline find the copies of a family?

    score.py --truth TRUTHPREFIX --strata STRATAPREFIX --asm ASMPREFIX --out OUT

TRUTHPREFIX.{transcripts,targets,named_family}.tsv come from sim_windows.py, STRATAPREFIX.{E0,R}.tsv from strata.py (run before the pipeline), ASMPREFIX = the driver's --out prefix
(`assemble` + `families`). For every family of targets.tsv (the family's copies, then the held-out family): holder, E2 own locus, E3 / E3' / E3c / E3p, E4 in the family, cluster
precision, IDEAL-FOUND = E2 and E3 and E4, the `nodes.py` own-node flag; the registered rule (YES / PARTLY / NO / CEILING-LIMITED) on R with every stratum printed beside.
Writes OUT.<family>.copies.tsv and OUT.summary.json. Coordinates 0-based half-open; a locus key is (chrom, 1-based gene-row start, end).
"""
import argparse
import collections
import csv
import json
import math
import re
import sys

PURITY = 0.5
COVER = 0.9
YES = 0.9
CP_MIN = 0.5
SHARED = 100


def rd(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def parse_chain(s):
    return tuple(tuple(map(int, x.split("-"))) for x in s.split(",")) if s else ()


def blocks(s):
    return [tuple(map(int, x.split("-"))) for x in s.split(",")] if s else []


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def bp(iv):
    return sum(e - s for s, e in iv)


def ov(a, b):
    i = j = tot = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi > lo:
            tot += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def chain_of(exons):
    ex = sorted(exons)
    return tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))


def span_key(s):
    m = re.match(r"^(.+):(\d+)-(\d+)$", s)
    return (m.group(1), int(m.group(2)), int(m.group(3)))


def load_loci(gff3):
    loci, by_key = {}, collections.defaultdict(list)
    for ln in open(gff3):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "gene":
            loci[at["Name"]] = dict(chrom=f[0], strand=f[6], s1=int(f[3]), e=int(f[4]), exons=[])
            by_key[(f[0], int(f[3]), int(f[4]))].append(at["Name"])
        elif f[2] == "exon":
            loci[at["gene"]]["exons"].append((int(f[3]) - 1, int(f[4])))
    for v in loci.values():
        v["exons"] = sorted(v["exons"])
        v["chain"] = chain_of(v["exons"])
        v["u"] = merge(v["exons"])
        v["bp"] = bp(v["u"])
    return loci, by_key


def load_gtf(path):
    """gene_id -> list of (transcript id, sorted exon blocks)."""
    ex = collections.defaultdict(list)
    gid = {}
    for ln in open(path):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        gid[t] = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        ex[t].append((int(f[3]) - 1, int(f[4])))
    out = collections.defaultdict(list)
    for t, e in ex.items():
        out[gid[t]].append((t, sorted(e)))
    return out


def lower(a, b):
    order = {"NO": 0, "PARTLY": 1, "YES": 2}
    return a if order.get(a, -1) <= order.get(b, -1) else b


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--strata", required=True)
    ap.add_argument("--asm", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--single", default=None, help="STRATAPREFIX.single_copy.tsv: the single-copy genes of the G7 negative control")
    a = ap.parse_args()
    trs, targets = rd(a.truth + ".transcripts.tsv"), rd(a.truth + ".targets.tsv")
    named = rd(a.truth + ".named_family.tsv")
    e0 = {r["cid"]: r for r in rd(a.strata + ".E0.tsv")}
    loci, by_key = load_loci(a.asm + ".fam.loci.gff3")
    fold = {span_key(r["annotation"]): span_key(r["representative"]) for r in rd(a.asm + ".fam.loci.tsv")}
    locus_tx = load_gtf(a.asm + ".families.gtf")
    cl_of_key, members = {}, collections.defaultdict(list)
    for r in rd(a.asm + ".fam.clusters.tsv"):
        k = (r["chrom"], int(r["start"]), int(r["end"]))
        if k not in by_key:
            sys.exit(f"INVALID: the cluster row {k} is on no locus")
        if len(by_key[k]) > 1:
            sys.exit(f"INVALID: the cluster row {k} is on a span shared by loci {by_key[k]}")
        cl_of_key[k] = r["cluster_id"]; members[r["cluster_id"]].append(by_key[k][0])
    name_cid = {}
    chains, union, tx_by_copy = collections.defaultdict(set), {}, collections.defaultdict(list)
    for r in trs:
        if r["copy"] and int(r["simulated"]):
            tx_by_copy[r["copy"]].append(r)
            if r["chain"]:
                chains[r["copy"]].add(parse_chain(r["chain"]))
    for cid, rows in tx_by_copy.items():
        union[cid] = merge([b for r in rows for b in blocks(r["canon_blocks"])])
    gene_union = collections.defaultdict(list)
    gene_info = {}
    for r in trs:
        if int(r["simulated"]):
            gene_union[r["gene"]].extend(blocks(r["canon_blocks"]))
            gene_info[r["gene"]] = (r["gene_name"], r["biotype"], r["strand"], r["chrom"])
    gene_union = {g: merge(b) for g, b in gene_union.items()}

    def resolve(name):
        L = loci[name]
        rk = fold.get((L["chrom"], L["s1"], L["e"]), (L["chrom"], L["s1"], L["e"]))
        names = by_key.get(rk)
        if names is None:
            return name
        if len(names) > 1 and rk != (L["chrom"], L["s1"], L["e"]):
            sys.exit(f"INVALID: the fold row to {rk} is on a span shared by loci {names}")
        return names[0] if rk != (L["chrom"], L["s1"], L["e"]) else name

    summaries = {}
    all_rows = {}
    for fam in sorted({t["family"] for t in targets}, key=lambda f: (any(t["heldout"] == "1" and t["family"] == f for t in targets), f)):
        tg = sorted([t for t in targets if t["family"] == fam], key=lambda t: t["cid"])
        rows, holder_of = [], {}
        for t in tg:
            cid, c_union, strand = t["cid"], union.get(t["cid"], []), t["strand"]
            mono = int(t["mono"])
            cand = []
            for name, L in loci.items():
                if L["chrom"] != t["chrom"] or not c_union:
                    continue
                o = ov(L["u"], c_union)
                if o > 0 and (mono or L["strand"] == strand):
                    cand.append((o / (L["bp"] + bp(c_union) - o), L["bp"], name, o))
            cand.sort(key=lambda x: (-x[0], x[1], x[2]))
            holder_of[cid] = resolve(cand[0][2]) if cand else None
            t["_folded"] = int(bool(cand) and holder_of[cid] != cand[0][2])
            top_any = max(((ov(L["u"], c_union), L["strand"]) for L in loci.values() if L["chrom"] == t["chrom"] and c_union), default=(0, strand))
            t["_split"] = len({resolve(x[2]) for x in cand if x[3] >= SHARED})
            t["_wrong_strand_top"] = int(top_any[0] > 0 and top_any[1] != strand and not mono)
            t["_cand"] = cand
        # K* : the cluster of the holder of the lowest-cid R copy that has a holder in a cluster
        kstar = None
        for t in tg:
            if int(e0[t["cid"]]["in_R"]) and holder_of[t["cid"]]:
                L = loci[holder_of[t["cid"]]]
                k = cl_of_key.get((L["chrom"], L["s1"], L["e"]))
                if k is not None:
                    kstar = k; break
        fam_union = merge([b for t in tg for b in union.get(t["cid"], [])])
        named_u = {n["gene"]: merge(blocks(n["canon_blocks"])) for n in named if (n["strand"] in "+-")}
        for t in tg:
            cid = t["cid"]; c_union = union.get(cid, []); H = holder_of[cid]
            L = loci.get(H) if H else None
            row = dict(cid=cid, name=t["name"], gene=t["gene"], family=fam, stratum_in=t["stratum"], in_R=int(e0[cid]["in_R"]), entangled_with=t["entangled_with"],
                       shared_bp=t["shared_bp"], shared_frac=t["shared_frac"], n_chains=t["n_multiexon_chains"], observable_P1=e0[cid]["observable_chains_P1"],
                       observable_P2=e0[cid]["observable_chains_P2"], observable_P3=e0[cid]["observable_chains_P3"], reads_back_share=e0[cid]["reads_back_share"],
                       holder=H or "", split_loci=t["_split"], wrong_strand_top=t["_wrong_strand_top"], folded_holder=t["_folded"])
            if L:
                o = ov(L["u"], c_union)
                rp = o / L["bp"]
                lu = merge([b for _, ex in locus_tx.get(H, []) for b in ex])
                lp = ov(lu, c_union) / bp(lu) if lu else 0.0
                row["rep_purity"], row["locus_purity"] = round(rp, 3), round(lp, 3)
                row["holder_shared_with"] = ";".join(sorted(x for x, h2 in holder_of.items() if h2 == H and x != cid))
                row["E2"] = int(not row["holder_shared_with"] and rp >= PURITY and lp >= PURITY)
                row["holder_cluster"] = cl_of_key.get((L["chrom"], L["s1"], L["e"]), "")
                multi_tx = [ex for _, ex in locus_tx.get(H, []) if len(ex) > 1]
                ch = sorted(chains.get(cid, []))
                if ch:
                    row["E3"] = int(L["chain"] in set(ch))
                    lo, hi = c_union[0][0], c_union[-1][1]
                    rep_in = tuple(j for j in L["chain"] if j[0] >= lo and j[1] <= hi)
                    ok = False
                    for c in ch:
                        need = max(2, math.ceil(len(c) / 2)); m = len(rep_in)
                        if m >= need and any(c[i:i + m] == rep_in for i in range(len(c) - m + 1)):
                            ok = True
                    row["E3r"] = int(row["E3"] or ok)
                    lt = {chain_of(ex) for ex in multi_tx}
                    row["E3c"] = round(sum(1 for c in ch if c in lt) / len(ch), 3)
                    all_tx = [ex for _, ex in locus_tx.get(H, [])]
                    row["E3p"] = round(sum(1 for ex in all_tx if chain_of(ex) not in set(ch)) / len(all_tx), 3) if all_tx else ""
                else:   # MONO: the representative covers >= 90% of the copy and has <= 110% of its bp
                    row["E3"] = row["E3r"] = int(c_union and o / bp(c_union) >= COVER and L["bp"] <= 1.1 * bp(c_union))
                    row["E3c"] = row["E3p"] = ""
                row["E4"] = int(kstar is not None and row["holder_cluster"] == kstar)
            else:
                row.update(rep_purity=0.0, locus_purity=0.0, holder_shared_with="", E2=0, holder_cluster="", E3=0, E3r=0, E3c="", E3p="", E4=0)
            row["IDEAL_FOUND"] = int(row["E2"] and row["E3"] and row["E4"])
            row["IDEAL_FOUND_relaxed"] = int(row["E2"] and row["E3r"] and row["E4"])
            # the registered own-node flag (nodes.py): a locus of a cluster that holds a locus overlapping a same-strand copy, overlapping this copy
            rows.append(row)
        # nodes.py own-node: clusters holding a locus that overlaps a same-strand copy of the family
        fam_clusters = set()
        for name, Lc in loci.items():
            k = cl_of_key.get((Lc["chrom"], Lc["s1"], Lc["e"]))
            if k is not None and any(Lc["chrom"] == t["chrom"] and Lc["strand"] == t["strand"] and ov(Lc["u"], union.get(t["cid"], [])) > 0 for t in tg):
                fam_clusters.add(k)
        for row, t in zip(rows, tg):
            row["own_node"] = int(any(cl_of_key.get((Lc["chrom"], Lc["s1"], Lc["e"])) in fam_clusters and Lc["strand"] == t["strand"] and ov(Lc["u"], union.get(t["cid"], [])) > 0
                                      for Lc in loci.values() if Lc["chrom"] == t["chrom"]))
            order = (("annotation:" + t["stratum"]) if t["stratum"] != "R_in" else "", "aligner" if t["stratum"] == "R_in" and not row["in_R"] else "",
                     "locus" if not row["E2"] else "", "family" if not row["E4"] else "", "representative" if not row["E3"] else "")
            row["miss_reason"] = next((x for x in order if x), "found") if not row["IDEAL_FOUND"] else "found"
        # cluster precision of K*
        other, cp_n, cp_d = [], 0, 0
        if kstar is not None:
            for m in members[kstar]:
                Lm = loci[m]; cp_d += 1
                if any(Lm["chrom"] == t["chrom"] and ov(Lm["u"], union.get(t["cid"], [])) > 0 for t in tg) or \
                   any(n["chrom"] == Lm["chrom"] and ov(Lm["u"], named_u[n["gene"]]) > 0 for n in named):
                    cp_n += 1
                else:
                    best = max(((ov(Lm["u"], u), g) for g, u in gene_union.items() if gene_info[g][3] == Lm["chrom"]), default=(0, ""))
                    other.append(f"{m}:{gene_info[best[1]][0] if best[0] else '-'}:{gene_info[best[1]][1] if best[0] else '-'}")
        cp = (cp_n / cp_d) if cp_d else 0.0
        R = [r for r in rows if r["in_R"]]; n = len(R); N = len(rows)
        s = lambda rs, k: sum(int(r[k]) for r in rs)
        need = math.ceil(YES * n) if n else 0
        if n < 0.5 * N:
            rule = "CEILING-LIMITED"
        elif s(R, "IDEAL_FOUND") >= need and cp >= CP_MIN:
            rule = "YES"
        elif s(R, "E2") >= need and s(R, "E4") >= need:
            rule = "PARTLY"
        else:
            rule = "NO"
        byst = {k: dict(copies=sum(1 for r in rows if r["stratum_in"] == k), IDEAL_FOUND=sum(int(r["IDEAL_FOUND"]) for r in rows if r["stratum_in"] == k)) for k in ("E", "C", "X", "R_in")}
        summaries[fam] = dict(family=fam, heldout=int(tg[0]["heldout"]), N=N, R=n, rule=(rule if not int(tg[0]["heldout"]) else "HELD-OUT (" + rule + ", no bar)"), n_folded_holders=sum(int(r["folded_holder"]) for r in rows), need=need, kstar=kstar, K_size=cp_d, CP=round(cp, 3), other_members=other,
                              R_counts=dict(E2=s(R, "E2"), E3=s(R, "E3"), E3r=s(R, "E3r"), E4=s(R, "E4"), IDEAL_FOUND=s(R, "IDEAL_FOUND"), IDEAL_FOUND_relaxed=s(R, "IDEAL_FOUND_relaxed")),
                              ALL_counts=dict(E2=s(rows, "E2"), E3=s(rows, "E3"), E4=s(rows, "E4"), IDEAL_FOUND=s(rows, "IDEAL_FOUND"), IDEAL_FOUND_relaxed=s(rows, "IDEAL_FOUND_relaxed"),
                                              own_node=s(rows, "own_node")),
                              by_stratum_in=byst, aligner_limited=sum(1 for r in rows if r["stratum_in"] == "R_in" and not r["in_R"]),
                              misses={k: sum(1 for r in rows if r["miss_reason"] == k) for k in sorted({r["miss_reason"] for r in rows})})
        all_rows[fam] = rows
        with open(f"{a.out}.{fam}.copies.tsv", "w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
        print(f"{fam}: N={N} R={n} rule={rule if not int(tg[0]['heldout']) else 'HELD-OUT/' + rule} IDEAL-FOUND on R {s(R, 'IDEAL_FOUND')}/{n} (need {need}), on ALL {s(rows, 'IDEAL_FOUND')}/{N}; E2 {s(R, 'E2')} E3 {s(R, 'E3')} E4 {s(R, 'E4')} on R; "
              f"CP {cp:.3f} (K*={kstar}, {cp_d} members); misses {summaries[fam]['misses']}")
    if a.single:   # G7: >= 95% of single-copy simulated genes get exactly one locus (same strand, >= 100 exonic bp) and none lies in a family's K*
        one = in_k = n_g = 0
        ks = {s["kstar"] for s in summaries.values() if not s["heldout"] and s["kstar"] is not None}
        for g in (r["gene"] for r in rd(a.single)):
            if g not in gene_union:
                continue
            n_g += 1
            hit = {resolve(nm) for nm, Lc in loci.items() if Lc["chrom"] == gene_info[g][3] and Lc["strand"] == gene_info[g][2] and ov(Lc["u"], gene_union[g]) >= SHARED}
            one += int(len(hit) == 1)
            in_k += int(any(cl_of_key.get((loci[h]["chrom"], loci[h]["s1"], loci[h]["e"])) in ks for h in hit))
        summaries["G7"] = dict(single_copy_genes=n_g, exactly_one_locus=one, share_one=round(one / n_g, 3) if n_g else None, in_Kstar=in_k, ok=bool(n_g and one / n_g >= 0.95 and in_k == 0))
        print("G7:", summaries["G7"])
    json.dump(summaries, open(a.out + ".summary.json", "w"), indent=1)


if __name__ == "__main__":
    main()
