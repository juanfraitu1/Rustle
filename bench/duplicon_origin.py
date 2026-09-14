#!/usr/bin/env python3
"""Prereg Addendum W: did a multi-copy gene family arise from duplicons? For each family, every member's region (gene
span +/- one gene-span length) is aligned to a reference member's region; a pair whose alignment chain continues
>= 1 kb past the gene on a matched side in BOTH members co-duplicated its flank (DUPLICON pair), otherwise the homology
stops at the gene (GENE-ONLY pair). A family is DUPLICON-BORNE when most of its aligned pairs co-duplicated flank.

Prereg Addendum X (refined, judged on held-out families): an aligned pair is GENOMIC when its chain shares >= 1 kb of
aligned bases that are non-exonic (introns or flanks) in BOTH member and reference, else RNA-LIKE; a family is
DUPLICON-DERIVED when >= 1 pair is aligned and >= 50% of aligned pairs are GENOMIC. Exons = union of the gene's
transcript exons, else the gene's own exon features, else the whole gene span.

usage: duplicon_origin.py --gff full.gff.gz --genome genome.fa --outdir DIR [--amy-truth truth.tsv] [--threads 4]
"""
import argparse
import bisect
import collections
import csv
import gzip
import math
import os
import pickle
import re
import statistics
import subprocess

import pysam

MIN_ID, EXT, CAP, CORE_FRAC, NONEXONIC = 0.80, 1000, 40, 0.90, 1000
FAMILIES = {
    "NPIP": (r"^NPIP[AB]\d+P?$", "core"), "TBC1D3": (r"^TBC1D3[A-Z]?$", "core"),
    "GOLGA8": (r"^GOLGA8[A-Z]*P?\d*$", "core"), "LRRC37A": (r"^LRRC37A\d*P?$", "core"), "RGPD": (r"^RGPD\d+$", "core"),
    "NBPF": (r"^NBPF\d+P?$", "core"), "SPATA31": (r"^SPATA31[A-Z]\d*P?\d*$", "core"), "PMS2": (r"^PMS2(P\d+)?$", "core"),
    "TRIM51": (r"^TRIM51[A-Z]*P?\d*$", "core"), "GUSB": (r"^GUSB(P\d+)?$", "core"),
    "GAPDH": (r"^GAPDH(P\d+)?$", "retro"), "PPIA": (r"^PPIA(P\d+)?$", "retro"), "EEF1A1": (r"^EEF1A1(P\d+)?$", "retro"),
    # Addendum X held-out families
    "NOTCH2NL": (r"^NOTCH2NL[A-Z]?$|^NOTCH2$", "sd"), "SRGAP2": (r"^SRGAP2[A-D]?$", "sd"),
    "ARHGAP11": (r"^ARHGAP11[AB]$", "sd"), "HYDIN": (r"^HYDIN\d?$", "sd"), "FAM72": (r"^FAM72[A-D]$", "sd"),
    "SMN": (r"^SMN[12]$", "sd"), "SERF1": (r"^SERF1[AB]$", "sd"), "ZNG1": (r"^ZNG1[A-F]$", "sd"),
    "NCF1": (r"^NCF1[BC]?$", "sd"),
    "RPL21": (r"^RPL21(P\d+)?$", "retro_ho"), "HMGB1": (r"^HMGB1(P\d+)?$", "retro_ho"),
    "NPM1": (r"^NPM1(P\d+)?$", "retro_ho"), "RPS2": (r"^RPS2(P\d+)?$", "retro_ho"), "TUBB": (r"^TUBB(P\d+)?$", "retro_ho"),
}


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def load_genes(gff):
    genes, tx_parent, exons = {}, {}, collections.defaultdict(int)
    feat_parent, exon_iv, cds_iv = {}, collections.defaultdict(list), collections.defaultdict(list)
    cds_tx = collections.defaultdict(list)
    with gzip.open(gff, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] in ("gene", "pseudogene"):
                genes[a["ID"]] = {"id": a["ID"], "name": a.get("Name", "?"), "chrom": f[0], "start0": int(f[3]) - 1,
                                  "end": int(f[4]), "strand": f[6], "biotype": a.get("gene_biotype", f[2])}
            elif f[2] in ("mRNA", "transcript", "lnc_RNA", "ncRNA", "primary_transcript") or f[2] == "exon":
                if f[2] == "exon":
                    exons[a.get("Parent", "")] += 1
                else:
                    tx_parent[a["ID"]] = a.get("Parent", "")
            if f[2] == "exon":
                exon_iv[a.get("Parent", "")].append((int(f[3]) - 1, int(f[4])))
            elif f[2] == "CDS":
                cds_iv[a.get("Parent", "")].append((int(f[3]) - 1, int(f[4])))
                cds_tx[a.get("Parent", "")].append((int(f[3]) - 1, int(f[4]), int(f[7]) if f[7] != "." else 0))
            elif "ID" in a and "Parent" in a:
                feat_parent[a["ID"]] = a["Parent"]

    def owner(p):
        up = p
        for _ in range(5):
            if up is None or up in genes:
                break
            up = feat_parent.get(up)
        return up if up in genes else None

    tx_ex, own_ex = collections.defaultdict(list), collections.defaultdict(list)
    tx_lists, cds = collections.defaultdict(list), collections.defaultdict(list)
    for p, iv in exon_iv.items():
        gid = owner(p)
        if gid is None:
            continue
        if p == gid:
            own_ex[gid] += iv
        else:
            tx_ex[gid] += iv
            tx_lists[gid].append(sorted(iv))
    for p, iv in cds_iv.items():
        gid = owner(p)
        if gid is not None:
            cds[gid] += iv
    tx_cds = collections.defaultdict(list)
    for p, segs in cds_tx.items():
        gid = owner(p)
        if gid is not None:
            tx_cds[gid].append(sorted(segs))
    for gid, g in genes.items():
        iv = tx_ex.get(gid) or own_ex.get(gid)
        g["exons"] = merge(iv) if iv else [(g["start0"], g["end"])]
        g["exon_source"] = "transcripts" if tx_ex.get(gid) else "gene" if own_ex.get(gid) else "span"
        g["transcripts"] = tx_lists.get(gid) or ([sorted(own_ex[gid])] if own_ex.get(gid) else [])
        g["cds"] = merge(cds.get(gid, []))
        g["tx_cds"] = tx_cds.get(gid, [])
    introns = collections.defaultdict(int)
    for tx, p in tx_parent.items():
        introns[p] = max(introns[p], exons.get(tx, 0) - 1)
    for gid, g in genes.items():
        g["introns"] = max(introns.get(gid, 0), exons.get(gid, 0) - 1, 0)
    return genes


def chains_from_paf(lines):
    by = collections.defaultdict(list)
    for line in lines:
        f = line.split("\t")
        by[(f[0], f[4])].append({"q": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4],
                                 "ts": int(f[7]), "te": int(f[8]), "nm": int(f[9]), "bl": int(f[10]),
                                 "cg": next((t[5:] for t in f[12:] if t.startswith("cg:Z:")), "")})
    chains = []
    for (q, strand), rs in by.items():
        rs.sort(key=lambda r: (r["ts"], r["te"]))
        Lq, cur = rs[0]["qlen"], []
        for r in rs + [None]:
            ok = False
            if r is not None and cur:
                gap = r["ts"] - max(x["te"] for x in cur)
                span = r["te"] - min(x["ts"] for x in cur)
                order = r["qs"] >= cur[-1]["qs"] if strand == "+" else r["qs"] <= cur[-1]["qs"]
                ok = gap <= Lq and span <= 2 * Lq and order
            if r is None or (cur and not ok):
                chains.append({"q": q, "strand": strand, "recs": cur, "nm": sum(x["nm"] for x in cur),
                               "bl": sum(x["bl"] for x in cur), "qs": min(x["qs"] for x in cur), "qe": max(x["qe"] for x in cur),
                               "ts": min(x["ts"] for x in cur), "te": max(x["te"] for x in cur),
                               "aligned": sum(e - s for s, e in merge([(x["qs"], x["qe"]) for x in cur]))})
                cur = []
            if r is not None:
                cur.append(r)
    return chains


def _hits(ex, ends, lo, hi):
    """Merged intervals of `ex` (sorted, `ends` = their ends) overlapping [lo, hi)."""
    i = bisect.bisect_right(ends, lo)
    while i < len(ex) and ex[i][0] < hi:
        yield ex[i]
        i += 1


def nonexonic_bp(recs, t_ex, q_ex):
    """Aligned bases (CIGAR M/=/X runs) outside the exons of BOTH sequences, counted once per reference position.
    t_ex / q_ex: merged exon intervals in reference / member region coordinates."""
    te, qe = [b for _, b in t_ex], [b for _, b in q_ex]
    out = []
    for r in recs:
        plus = r["strand"] == "+"
        t, q = r["ts"], (r["qs"] if plus else r["qe"] - 1)
        for n, op in re.findall(r"(\d+)([MIDNSHP=X])", r["cg"]):
            n = int(n)
            if op in "M=X":
                blocked = [(a - t, b - t) for a, b in _hits(t_ex, te, t, t + n)]
                if plus:
                    blocked += [(a - q, b - q) for a, b in _hits(q_ex, qe, q, q + n)]
                else:
                    blocked += [(q - b + 1, q - a + 1) for a, b in _hits(q_ex, qe, q - n + 1, q + 1)]
                cur = 0
                for s_, e_ in merge([(max(0, x), min(n, y)) for x, y in blocked if max(0, x) < min(n, y)]) + [(n, n)]:
                    if s_ > cur:
                        out.append((t + cur, t + s_))
                    cur = max(cur, e_)
                t += n
                q += n if plus else -n
            elif op in "DN":
                t += n
            elif op == "I":
                q += n if plus else -n
    return sum(e - s for s, e in merge(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--amy-truth")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--families", help="comma list (default all); PAFs and the parsed GFF are cached in --outdir")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    genome = pysam.FastaFile(a.genome)
    clen = dict(zip(genome.references, genome.lengths))
    gcache = f"{a.outdir}/genes.pkl"
    if os.path.exists(gcache):
        genes = pickle.load(open(gcache, "rb"))
    else:
        genes = load_genes(a.gff)
        pickle.dump(genes, open(gcache, "wb"))
    by_name = collections.defaultdict(list)
    for g in genes.values():
        by_name[g["name"]].append(g)

    fams = {}
    wanted = set(a.families.split(",")) if a.families else set(FAMILIES) | {"AMY"}
    for fam, (pat, kind) in FAMILIES.items():
        if fam not in wanted:
            continue
        members = [g for n, gs in by_name.items() if re.match(pat, n) for g in gs if g["chrom"] in clen]
        fams[fam] = (kind, members)
    if a.amy_truth and "AMY" in wanted:
        amy = []
        for r in csv.DictReader(open(a.amy_truth), delimiter="\t"):
            gs = [g for g in by_name.get(r["name"], []) if g["chrom"] == r["chrom"]]
            base = gs[0] if gs else {"introns": 0, "exons": [(int(r["start0"]), int(r["end"]))], "exon_source": "span"}
            amy.append({"name": r["name"], "chrom": r["chrom"], "start0": int(r["start0"]), "end": int(r["end"]),
                        "strand": r["strand"], "biotype": r["biotype"], "introns": base["introns"],
                        "exons": base["exons"], "exon_source": base["exon_source"]})
        fams["AMY"] = ("reported", amy)

    pair_rows, fam_rows = [], []
    for fam, (kind, members) in fams.items():
        pcs = [g for g in members if g["biotype"] == "protein_coding"] or members
        med = statistics.median(g["end"] - g["start0"] for g in pcs)
        ref = min(pcs, key=lambda g: (abs((g["end"] - g["start0"]) - med), g["name"]))
        others = sorted([g for g in members if g is not ref], key=lambda g: (g["name"], g["chrom"], g["start0"]))
        dropped = max(0, len(others) - CAP)
        others = others[:CAP]

        def region(g):
            L = g["end"] - g["start0"]
            s, e = max(0, g["start0"] - L), min(clen[g["chrom"]], g["end"] + L)
            return s, e

        rs, re_ = region(ref)
        tfa, qfa = f"{a.outdir}/{fam}.ref.fa", f"{a.outdir}/{fam}.members.fa"
        with open(tfa, "w") as fh:
            fh.write(f">ref\n{genome.fetch(ref['chrom'], rs, re_).upper()}\n")
        keys = {}
        with open(qfa, "w") as fh:
            for i, g in enumerate(others):
                s, e = region(g)
                keys[f"m{i}"] = (g, s, e)
                fh.write(f">m{i}\n{genome.fetch(g['chrom'], s, e).upper()}\n")
        pafp = f"{a.outdir}/{fam}.paf"
        if not os.path.exists(pafp):
            out = subprocess.run(["minimap2", "-c", "-x", "asm20", "-N", "50", "-p", "0.1", "-t", str(a.threads), tfa, qfa],
                                 capture_output=True, text=True, check=True).stdout
            open(pafp + ".tmp", "w").write(out)
            os.replace(pafp + ".tmp", pafp)
        paf = open(pafp).read().splitlines()
        chains = chains_from_paf(paf)
        gt0, gt1 = ref["start0"] - rs, ref["end"] - rs
        cov = collections.Counter()
        n_aligned = n_dup = n_gen = 0
        t_ex = merge([(x - rs, y - rs) for x, y in ref["exons"]])
        for key, (g, s, e) in keys.items():
            gq0, gq1 = g["start0"] - s, g["end"] - s
            cand = [c for c in chains if c["q"] == key and c["qs"] < gq1 and gq0 < c["qe"] and c["ts"] < gt1 and gt0 < c["te"]]
            row = {"family": fam, "kind": kind, "reference": ref["name"], "member": g["name"], "chrom": g["chrom"],
                   "member_introns": g["introns"], "ref_introns": ref["introns"], "identity": float("nan"),
                   "ext_q_left": 0, "ext_q_right": 0, "ext_t_left": 0, "ext_t_right": 0, "class": "UNALIGNED",
                   "member_exon_source": g["exon_source"], "ref_exon_source": ref["exon_source"], "shared_nonexonic": 0,
                   "class_x": "UNALIGNED"}
            if cand:
                c = max(cand, key=lambda c: c["aligned"])
                ident = c["nm"] / c["bl"]
                row["identity"] = ident
                if ident >= MIN_ID:
                    eql, eqr = gq0 - c["qs"], c["qe"] - gq1
                    etl, etr = gt0 - c["ts"], c["te"] - gt1
                    row.update(ext_q_left=eql, ext_q_right=eqr, ext_t_left=etl, ext_t_right=etr)
                    matched = [(eql, etl), (eqr, etr)] if c["strand"] == "+" else [(eql, etr), (eqr, etl)]
                    dup = any(x >= EXT and y >= EXT for x, y in matched)
                    row["class"] = "DUPLICON" if dup else "GENE-ONLY"
                    n_aligned += 1
                    n_dup += dup
                    ne = nonexonic_bp(c["recs"], t_ex, merge([(x - s, y - s) for x, y in g["exons"]]))
                    row["shared_nonexonic"] = ne
                    row["class_x"] = "GENOMIC" if ne >= NONEXONIC else "RNA-LIKE"
                    n_gen += ne >= NONEXONIC
                    for s_, e_ in merge([(r["ts"], r["te"]) for r in c["recs"]]):
                        for p in range(s_, e_):
                            cov[p] += 1
                else:
                    row["class"] = row["class_x"] = "LOW-IDENTITY"
            pair_rows.append(row)
        need = math.ceil(CORE_FRAC * n_aligned) if n_aligned else 0
        core = sum(1 for v in cov.values() if v >= need) if n_aligned else 0
        call = "DUPLICON-BORNE" if n_aligned >= 2 and n_dup / n_aligned >= 0.5 else "NOT"
        fam_rows.append({"family": fam, "kind": kind, "members": len(members), "dropped_by_cap": dropped,
                         "reference": ref["name"], "ref_span": ref["end"] - ref["start0"], "aligned": n_aligned,
                         "duplicon_pairs": n_dup, "dup_frac": n_dup / n_aligned if n_aligned else float("nan"),
                         "core_bp": core, "core_over_gene": core / (ref["end"] - ref["start0"]), "call": call,
                         "genomic_pairs": n_gen, "gen_frac": n_gen / n_aligned if n_aligned else float("nan"),
                         "call_x": "DUPLICON-DERIVED" if n_aligned >= 1 and n_gen / n_aligned >= 0.5 else "NOT"})

    tag = "" if not a.families else "." + "_".join(sorted(wanted))
    with open(f"{a.outdir}/pairs{tag}.tsv", "w") as fh:
        cols = list(pair_rows[0].keys())
        fh.write("\t".join(cols) + "\n")
        for r in pair_rows:
            fh.write("\t".join(f"{r[k]:.4f}" if isinstance(r[k], float) else str(r[k]) for k in cols) + "\n")
    print("family\tkind\tmembers\tdropped\treference\tref_span\taligned\tduplicon_pairs\tdup_frac\tcore_bp\tcore/gene\tcall"
          "\tgenomic_pairs\tgen_frac\tcall_x")
    for r in fam_rows:
        print(f"{r['family']}\t{r['kind']}\t{r['members']}\t{r['dropped_by_cap']}\t{r['reference']}\t{r['ref_span']}\t"
              f"{r['aligned']}\t{r['duplicon_pairs']}\t{r['dup_frac']:.3f}\t{r['core_bp']}\t{r['core_over_gene']:.2f}\t{r['call']}"
              f"\t{r['genomic_pairs']}\t{r['gen_frac']:.3f}\t{r['call_x']}")
    core_dup = sum(1 for r in fam_rows if r["kind"] == "core" and r["call"] == "DUPLICON-BORNE")
    retro_dup = sum(1 for r in fam_rows if r["kind"] == "retro" and r["call"] == "DUPLICON-BORNE")
    n_core = sum(1 for r in fam_rows if r["kind"] == "core")
    n_retro = sum(1 for r in fam_rows if r["kind"] == "retro")
    verdict = ("SUPPORTED" if core_dup >= 7 and retro_dup <= 1 else
               "UNINFORMATIVE (classifier calls retro families duplicon-borne)" if retro_dup >= 2 else "NOT SUPPORTED")
    print(f"\nREADING: core-duplicon families DUPLICON-BORNE {core_dup}/{n_core}; retrotransposition families DUPLICON-BORNE "
          f"{retro_dup}/{n_retro} -> {verdict}")
    nx = lambda kind: (sum(1 for r in fam_rows if r["kind"] == kind and r["call_x"] == "DUPLICON-DERIVED"),
                       sum(1 for r in fam_rows if r["kind"] == kind))
    print(f"DEVELOPMENT (Addendum X rule, not deciding): core DUPLICON-DERIVED {nx('core')[0]}/{nx('core')[1]}; "
          f"retro {nx('retro')[0]}/{nx('retro')[1]}; AMY {nx('reported')[0]}/{nx('reported')[1]}")
    sd, rho = nx("sd"), nx("retro_ho")
    if sd[1]:
        vx = ("SUPPORTED" if sd[0] >= 7 and rho[0] <= 1 else
              "DOES NOT DISCRIMINATE (>= 2 retro families DUPLICON-DERIVED)" if rho[0] >= 2 else "NOT SUPPORTED")
        print(f"READING X (held out): segmental-duplication families DUPLICON-DERIVED {sd[0]}/{sd[1]}; retrotransposition "
              f"families {rho[0]}/{rho[1]} -> {vx}")
    retro_sig = collections.Counter((r["family"], r["member_introns"] == 0 and r["ref_introns"] >= 1) for r in pair_rows)
    print("retro signature (member 0 introns, reference >= 1) per family:",
          {f: f"{retro_sig[(f, True)]}/{retro_sig[(f, True)] + retro_sig[(f, False)]}" for f in fams})


if __name__ == "__main__":
    main()
