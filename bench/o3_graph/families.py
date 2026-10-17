#!/usr/bin/env python3
"""Handpicked real gene families for the graph-O2/O3 simulation v2 (prereg docs/PREREG_o3_graph_real_2026-10-15.md): mainstay
multi-copy families (NPIP, TBC1D3, GOLGA, Soto's and others), their copies and curated isoforms from CHM13 v2.0's CAT/Liftoff
annotation (the project's human default).

    families.py build OUT.json [--cap 10] [--seed 1] [--min-cover 0.5]

Per handpicked group (a CAT gene-name pattern): the genes matching it, span <= MAX_SPAN, non-overlapping (an overlap keeps the
gene with more transcripts); their canonical transcripts (most exons, then longest span, then first id) spliced from the genome
on their strand; copies linked when their canonicals align at identity >= MIN_IDENT over >= half the shorter (minimap2 -cx
asm20); the largest linked component is the family, a seeded sample of `cap` copies when larger. Isoforms per copy: the
canonical (the MAIN isoform) plus up to ISO_CAP - 1 others with distinct intron chains, in id order.
"""
import argparse
import collections
import json
import os
import random
import re
import subprocess
import sys
import tempfile

CAT = "/mnt/linuxdisk/home/juanfraitu/winloci_data/gencode_chm13"
GENOME = "/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa"
MAX_SPAN, MIN_IDENT, ISO_CAP, MIN_COVER = 60000, 0.90, 4, 0.5
# handpicked before any run: development families, then held-out ones (NPIP, TBC1D3 and GOLGA8 are Soto's mainstays)
GROUPS = {
    "dev": {"SMN": r"^SMN[12]$", "NCF1": r"^NCF1[BC]?$", "PAGE": r"^PAGE[0-9]", "CT45A": r"^CT45A", "GAGE": r"^GAGE",
            "POTE": r"^POTE", "GSTA": r"^GSTA[0-9]", "TP53TG3": r"^TP53TG3"},
    "heldout": {"NPIPA": r"^NPIPA", "NPIPB": r"^NPIPB", "TBC1D3": r"^TBC1D3", "GOLGA8": r"^GOLGA8", "GOLGA6L": r"^GOLGA6L",
                "CT47A": r"^CT47A", "XAGE": r"^XAGE", "SPATA31": r"^SPATA31"},
    # held out for union O2 / G (docs/PREREG_o3_graph_union_2026-10-15.md), picked before any run: Soto's human-specific
    # duplications FAM72, ARHGAP11, TCAF2, and copy-number mainstays AMY, SULT1A, FCGR, APOBEC3, PSG; ARHGAP11 and APOBEC3 form
    # no family at MIN_IDENT (built, no world), so GTF2H2 (the SMA region) and LRRC37A (17q21.31) were added, still before any run
    "heldout2": {"FAM72": r"^FAM72[A-D]$", "ARHGAP11": r"^ARHGAP11[AB]$", "TCAF2": r"^TCAF2", "AMY": r"^AMY[12]",
                 "SULT1A": r"^SULT1A[1-4]$", "FCGR": r"^FCGR[23]", "APOBEC3": r"^APOBEC3[A-H]$", "PSG": r"^PSG[0-9]+$",
                 "GTF2H2": r"^GTF2H2", "LRRC37A": r"^LRRC37A"},
    # held out for round 3 (docs/PREREG_o3_graph_union2_2026-10-15.md), picked before any run: X-linked cancer-testis families
    # SSX, MAGEA, SPANX; UGT2B; the opsins OPN1LW/MW; CGB; KIR; CYP4F (spares, in order, should one form no family: TRIM,
    # LILR, DEFB10, CSAG)
    "heldout3": {"SSX": r"^SSX[0-9]+[AB]?$", "MAGEA": r"^MAGEA[0-9]+[AB]?$", "SPANX": r"^SPANX", "UGT2B": r"^UGT2B[0-9]+$",
                 "OPN1": r"^OPN1[LM]W", "CGB": r"^CGB[0-9]*$", "KIR": r"^KIR[23]D", "CYP4F": r"^CYP4F[0-9]+$"},
    # held out for round 4 (docs/PREREG_o3_graph_round4_2026-10-15.md), picked before any run, in this order; the first 8 with
    # >= 3 copies are the set (a removed copy must leave >= 2 catalogued hosts: round 3's CYP4F had none)
    "heldout4": {"TRIM": r"^TRIM(43|49|51|64)[A-Z]?$", "LILR": r"^LILR[AB][0-9]$", "DEFB10": r"^DEFB10[3-7][AB]$",
                 "CSAG": r"^CSAG[0-9]$", "MAGEB": r"^MAGEB[0-9]+$", "CTAG": r"^CTAG[12][AB]?$", "SERF1": r"^SERF1[AB]?$",
                 "CBWD": r"^CBWD[0-9]$", "MAGEC": r"^MAGEC[0-9]$", "CYP2A": r"^CYP2A[0-9]+$", "NBPF": r"^NBPF[0-9]+$",
                 "GTF2I": r"^GTF2I", "PRAMEF": r"^PRAMEF[0-9]+$", "VCX": r"^VCX[0-9]?[AB]?$", "SPDYE": r"^SPDYE[0-9]+$",
                 "FAM86": r"^FAM86[A-Z]$", "GPR89": r"^GPR89[ABC]$", "STAG3L": r"^STAG3L[0-9]$", "POM121": r"^POM121[A-Z]?$",
                 "CCDC144": r"^CCDC144[A-Z]+$"},
    # held out for round 5 (docs/PREREG_o3_graph_round5_2026-10-16.md), picked before any run, in this order; the first 8 with
    # >= 3 copies are the set; the TRIM and LILR subfamilies come first so that sister families (another family's reads hosted
    # on a family: round 4's dev-2 failure) can occur
    "heldout5": {"TRIM64": r"^TRIM64[A-Z]?$", "TRIM43": r"^TRIM43[A-Z]?$", "TRIM51": r"^TRIM51[A-Z]*$", "LILRA": r"^LILRA[0-9]$",
                 "LILRB": r"^LILRB[0-9]$", "SPDYE": r"^SPDYE[0-9]+[A-Z]?$", "STAG3L": r"^STAG3L[0-9]$", "MT1": r"^MT1[A-Z]$",
                 "GSTM": r"^GSTM[1-5]$", "CEACAM": r"^CEACAM[0-9]+$", "FAM230": r"^FAM230[A-Z]$", "HERC2P": r"^HERC2P[0-9]+$",
                 "CYP2D": r"^CYP2D[0-9]+P?$", "BTN3A": r"^BTN3A[1-3]$", "SPANXN": r"^SPANXN[0-9]$",
                 "GUSBP": r"^GUSBP[0-9]+$", "ANKRD20A": r"^ANKRD20A[0-9]+P?$", "FRG2": r"^FRG2[A-Z]*$", "GOLGA6": r"^GOLGA6[A-D]$",
                 "PMS2": r"^PMS2(CL|P[0-9]+)?$"},
    # held out for round 6 (docs/PREREG_o3_graph_round6_2026-10-16.md), picked before any run, in this order; the first 8 with
    # >= 3 copies are the set (FRG2, GOLGA6 and PMS2 were round-5 candidates past its first 8, never built into a world)
    "heldout6": {"FRG2": r"^FRG2[A-Z]*$", "GOLGA6": r"^GOLGA6[A-D]$", "PMS2": r"^PMS2(CL|P[0-9]+)?$", "CCL34": r"^CCL[34](L[0-9])?$",
                 "DEFA": r"^DEFA[0-9]+[A-Z]?$", "BCRP": r"^BCRP[0-9]+$", "GGT": r"^GGT[0-9]+[A-Z]*$", "CFHR": r"^CFHR[0-9]$",
                 "PRB": r"^PRB[0-9]$", "ANKRD30": r"^ANKRD30[A-Z]+$", "FRG1": r"^FRG1[A-Z]*$", "GSTT2": r"^GSTT2[A-Z]?$",
                 "CTAGE": r"^CTAGE[0-9]+$", "RIMBP3": r"^RIMBP3[A-Z]?$"},
    # held out for round 7 (docs/PREREG_o3_graph_round7_2026-10-16.md), picked before any run, in this order; the first 8 with
    # >= 3 copies are the set (CTAGE was a round-6 candidate past its first 8, never built into a world)
    "heldout7": {"CTAGE": r"^CTAGE[0-9]+$", "AQP7": r"^AQP7(P[0-9]+)?$", "PRSS": r"^PRSS[1-3](P[0-9])?$", "SPRR2": r"^SPRR2[A-G]$",
                 "FAM74A": r"^FAM74A[0-9]+$", "TUBB8": r"^TUBB8(P[0-9]+)?$", "DDX11L": r"^DDX11L[0-9]+$", "WASH": r"^WASH[0-9]P?$",
                 "SIGLEC": r"^SIGLEC[0-9]+[A-Z]?$", "IFITM": r"^IFITM[0-9]+$", "FAM138": r"^FAM138[A-Z]$", "FAM157": r"^FAM157[A-Z]$"},
    # held out for round 8 (docs/PREREG_o3_graph_round8_2026-10-16.md), picked before any run, in this order; the first 8 with
    # >= 3 copies are the set. Human-specific duplications (NOTCH2NL, SRGAP2) first, then copy-number mainstays (GH/CSH, USP17L,
    # the Y ampliconic TSPY, RBMY, BPY2, CDY), then others; none was in any earlier world (FAM86, GPR89 and POM121 were round-4
    # candidates past its first 8, never built into a world)
    "heldout8": {"NOTCH2NL": r"^NOTCH2NL[A-CR]?$", "SRGAP2": r"^SRGAP2[B-D]?$", "GH": r"^(GH[12]|CSH[12]|CSHL1)$",
                 "USP17L": r"^USP17L[0-9]+$", "TSPY": r"^TSPY[0-9]+$", "RBMY": r"^RBMY1[A-J][0-9]?$", "BPY2": r"^BPY2[BC]?$",
                 "CDY": r"^CDY[12][AB]?$", "OR4F": r"^OR4F[0-9]+$", "LCE1": r"^LCE1[A-F]$", "KRTAP5": r"^KRTAP5-[0-9]+$",
                 "PNMA6": r"^PNMA6[A-F]$", "PPIAL4": r"^PPIAL4[A-H]$", "FCGR1": r"^FCGR1[A-C]P?$", "NOMO": r"^NOMO[1-3]$",
                 "CROCC": r"^CROCC(P[0-9])?$", "RGPD": r"^RGPD[1-8]$", "FAM86": r"^FAM86[A-Z]$", "GPR89": r"^GPR89[ABC]$",
                 "POM121": r"^POM121[A-Z]?$"},
    # held out for round 9 (docs/PREREG_o3_graph_round9_2026-10-16.md), picked before any run, in this order; the first 8 with
    # >= 3 copies are the set, UNCAPPED (full families, --cap 1000). PPIAL4 and FCGR1 were round-8 candidates past its first 8,
    # never built into a world; none of the others was in any earlier list
    "heldout9": {"PPIAL4": r"^PPIAL4[A-H]$", "FCGR1": r"^FCGR1[A-C]P?$", "PGA": r"^PGA[3-5]$", "H4C": r"^H4C[0-9]+$",
                 "LCE2": r"^LCE2[A-D]$", "LCE3": r"^LCE3[A-E]$", "KRTAP4": r"^KRTAP4-[0-9]+$", "AGAP": r"^AGAP([4-9]|1[01])P?$",
                 "ELOA3": r"^ELOA3[B-D]?P?$", "FAM25": r"^FAM25[A-G]$", "CYP3A": r"^CYP3A[0-9]+$", "FAM47": r"^FAM47[A-C]$",
                 "FAM9": r"^FAM9[A-C]$", "KRTAP9": r"^KRTAP9-[0-9]+$", "DAZ": r"^DAZ[1-4]$"},
    # held out for round 10 (docs/PREREG_o3_graph_round10_2026-10-16.md), picked before any run, in this order; the first 8 with
    # >= 3 copies are the set, UNCAPPED, linked by the locus rule (--link locus --min-cover 0.8). CYP3A, FAM47, FAM9 and DAZ were
    # round-9 candidates past its first 8, never built into a world (KRTAP9 was, in the rebuilt set 10); none of the others was in
    # any earlier list
    "heldout10": {"CYP3A": r"^CYP3A[0-9]+$", "FAM47": r"^FAM47[A-C]$", "FAM9": r"^FAM9[A-C]$", "DAZ": r"^DAZ[1-4]$",
                  "HNRNPCL": r"^HNRNPCL[0-9]+$", "PRR20": r"^PRR20[A-G]$", "ZNF705": r"^ZNF705[A-G]P?$", "TAF11L": r"^TAF11L[0-9]+$",
                  "F8A": r"^F8A[1-3]$", "FAM236": r"^FAM236[A-D]$", "CENPVL": r"^CENPVL[0-9]$", "FOXD4L": r"^FOXD4L[0-9]$",
                  "OR2T": r"^OR2T[0-9]+$", "KRTAP10": r"^KRTAP10-[0-9]+$", "H2BC": r"^H2BC[0-9]+$", "FAM66": r"^FAM66[A-E]$",
                  "CXADRP": r"^CXADRP[0-9]$", "GOLGA2P": r"^GOLGA2P[0-9]+$"},
    # held out for round 11 (docs/PREREG_o3_graph_round11_2026-10-16.md), picked before any run, in this order; the first 8 with
    # >= 3 copies are the set, UNCAPPED, linked by the locus rule (--link locus --min-cover 0.8). FOXD4L, OR2T, KRTAP10, H2BC,
    # FAM66, CXADRP and GOLGA2P were round-10 candidates past its first 8, never built into a world; none of the others was in
    # any earlier list
    "heldout11": {"FOXD4L": r"^FOXD4L[0-9]$", "OR2T": r"^OR2T[0-9]+$", "KRTAP10": r"^KRTAP10-[0-9]+$", "H2BC": r"^H2BC[0-9]+$",
                  "FAM66": r"^FAM66[A-E]$", "CXADRP": r"^CXADRP[0-9]$", "GOLGA2P": r"^GOLGA2P[0-9]+$", "H3C": r"^H3C[0-9]+$",
                  "H2AC": r"^H2AC[0-9]+$", "KRTAP2": r"^KRTAP2-[0-9]+$", "KRTAP1": r"^KRTAP1-[0-9]+$", "PRR23": r"^PRR23[A-D]$",
                  "ETD": r"^ETD[ABC]$", "SULT1C": r"^SULT1C[0-9]$", "UGT2A": r"^UGT2A[1-3]$", "CDRT15": r"^CDRT15(L[0-9]+|P[0-9]+)?$",
                  "NAIP": r"^NAIP(P[0-9]+)?$", "HLA-DRB": r"^HLA-DRB[0-9]$", "TCEAL": r"^TCEAL[0-9]$", "ARMCX": r"^ARMCX[0-9]$",
                  "OR4C": r"^OR4C[0-9]+$", "OR4K": r"^OR4K[0-9]+$", "OR11H": r"^OR11H[0-9]+$", "SEC14L": r"^SEC14L[0-9]$",
                  "DUX4": r"^DUX4(L[0-9]+)?$"},
    # held out for round 12 (docs/PREREG_o3_graph_round12_2026-10-16.md), picked before any run, in this order; the first 8 with
    # >= 3 copies are the set, UNCAPPED, linked by the locus rule. KRTAP2, KRTAP1, ETD, CDRT15, NAIP, HLA-DRB, OR11H and DUX4
    # were round-11 candidates past its first 8, never built into a world; none of the others was in any earlier list
    "heldout12": {"KRTAP2": r"^KRTAP2-[0-9]+$", "KRTAP1": r"^KRTAP1-[0-9]+$", "ETD": r"^ETD[ABC]$",
                  "CDRT15": r"^CDRT15(L[0-9]+|P[0-9]+)?$", "NAIP": r"^NAIP(P[0-9]+)?$", "HLA-DRB": r"^HLA-DRB[0-9]$",
                  "OR11H": r"^OR11H[0-9]+$", "DUX4": r"^DUX4(L[0-9]+)?$", "POLR2J": r"^POLR2J[23]?$",
                  "KRTAP19": r"^KRTAP19-[0-9]+$", "KRTAP3": r"^KRTAP3-[0-9]+$", "MRGPRX": r"^MRGPRX[1-4]$"},
    # held out for round 13 (docs/PREREG_o3_graph_round13_2026-10-17.md), picked before any run, in this order; the first 8 with
    # >= 3 copies are the set, UNCAPPED, linked by the locus rule. MRGPRX was a round-12 candidate past its first 8 with >= 3
    # copies, never built into a world; none of the others was in any earlier list
    "heldout13": {"MRGPRX": r"^MRGPRX[1-4]$", "S100A7": r"^S100A7([AL]2?)?$", "KRTAP12": r"^KRTAP12-[0-9]+$",
                  "KRTAP13": r"^KRTAP13-[0-9]+$", "KRTAP20": r"^KRTAP20-[0-9]+$", "KRTAP21": r"^KRTAP21-[0-9]+$",
                  "KRTAP6": r"^KRTAP6-[0-9]+$", "HLA-ABC": r"^HLA-[ABC]$", "H1": r"^H1-[0-9]$", "OR2A": r"^OR2A[0-9]+$",
                  "SEPTIN7P": r"^SEPTIN7P[0-9]+$", "CTSLP": r"^CTSLP[0-9]+$", "UBE2Q2P": r"^UBE2Q2P[0-9]+$",
                  "PGAM1P": r"^PGAM1P[0-9]+$", "OR52": r"^OR52[A-Z][0-9]+$", "OR7E": r"^OR7E[0-9]+P?$"},
}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def canonical(tx):
    """Most exons, then the longest span, then the first id."""
    return min(tx, key=lambda t: (-len(tx[t]), -(max(e for _, e in tx[t]) - min(s for s, _ in tx[t])), t))


def chain(ex):
    ex = sorted(ex)
    return tuple((a[1], b[0]) for a, b in zip(ex, ex[1:]))


def pick_isoforms(tx, cap=ISO_CAP):
    """(main id, {id: exons}): the canonical first, then transcripts with a new intron chain (single-exon ones skipped)."""
    main = canonical(tx)
    out, seen = {main: sorted(tx[main])}, {chain(tx[main])}
    for t in sorted(tx):
        c = chain(tx[t])
        if len(out) >= cap:
            break
        if t != main and c and c not in seen:
            seen.add(c)
            out[t] = sorted(tx[t])
    return main, out


def drop_overlaps(genes):
    keep = []
    for g in sorted(genes, key=lambda g: (-g["n"], g["chrom"], g["start"])):
        if not any(k["chrom"] == g["chrom"] and k["start"] < g["end"] and g["start"] < k["end"] for k in keep):
            keep.append(g)
    return sorted(keep, key=lambda g: (g["chrom"], g["start"]))


def largest_component(names, hits, min_ident=MIN_IDENT, min_cover=MIN_COVER):
    """hits: (a, b, identity, block, len a, len b); linked at identity >= min_ident over >= min_cover of the shorter (default
    half; `--min-cover 0.8`: whole-transcript copies, not genes that share only a segment)."""
    par = {n: n for n in names}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for a, b, ident, blk, la, lb in hits:
        if a != b and a in par and b in par and ident >= min_ident and blk >= min_cover * min(la, lb):
            par[find(a)] = find(b)
    comp = collections.defaultdict(list)
    for n in names:
        comp[find(n)].append(n)
    return sorted(max(comp.values(), key=lambda c: (len(c), sorted(c)[0] < "~")))


def union_cover(blocks, length):
    """The union of [start, end) blocks over a length."""
    tot, cs, ce = 0, None, None
    for s0, e0 in sorted(blocks):
        if ce is None or s0 > ce:
            tot += (ce - cs) if ce is not None else 0
            cs, ce = s0, e0
        else:
            ce = max(ce, e0)
    tot += (ce - cs) if ce is not None else 0
    return tot / length if length else 0.0


def tx_locus_cover(rows, min_ident=MIN_IDENT):
    """rows (query isoform 'GENE|tx', query length, qs, qe, target locus GENE, identity) -> {(A, B): over A's isoforms, the best
    union of blocks aligned to B's locus at >= min_ident, over the isoform's length} (spliced to the partner's GENOME, so its
    intron divergence and its annotation's canonical choice do not enter)."""
    blocks, qlen = collections.defaultdict(list), {}
    for q, ql, qs, qe, t, ident in rows:
        g = q.split("|")[0]
        if g == t or ident < min_ident:
            continue
        blocks[(q, t)].append((qs, qe))
        qlen[q] = ql
    out = {}
    for (q, t), b in blocks.items():
        k = (q.split("|")[0], t)
        out[k] = max(out.get(k, 0.0), union_cover(b, qlen[q]))
    return out


def locus_cover(genes, fa, pad=1000):
    """tx_locus_cover for genes [{gene_id, chrom, start, end, strand, tx}]: every isoform, spliced (minimap2 -cx splice:hq), onto
    every gene's locus (+- pad)."""
    with tempfile.TemporaryDirectory() as d:
        with open(f"{d}/q.fa", "w") as fh:
            for g in genes:
                for t, ex in g["tx"].items():
                    fh.write(f">{g['gene_id']}|{t}\n{spliced(fa, g['chrom'], g['strand'], ex)}\n")
        with open(f"{d}/t.fa", "w") as fh:
            for g in genes:
                fh.write(f">{g['gene_id']}\n{fa.fetch(g['chrom'], max(0, g['start'] - pad), g['end'] + pad).upper()}\n")
        paf = subprocess.run(["minimap2", "-cx", "splice:hq", "-N", "1000", "-p", "0.01", "--secondary=yes", f"{d}/t.fa", f"{d}/q.fa"],
                             capture_output=True, text=True, check=True).stdout.splitlines()
    return tx_locus_cover([(x[0], int(x[1]), int(x[2]), int(x[3]), x[5], int(x[9]) / int(x[10]))
                           for x in (ln.split("\t") for ln in paf)])


def locus_component(names, cov, min_cover=0.8):
    """`--link locus` (2026-10-16): genes linked when one gene's genomic locus holds >= min_cover of the other's canonical
    transcript at >= MIN_IDENT (whole copies and fragments link; genes sharing only a segment -- domain sharers -- do not; the
    partner's locus, not its annotation, so its intron divergence and canonical choice do not enter); the largest component."""
    par = {n: n for n in names}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for (a, b), c in cov.items():
        if a in par and b in par and max(c, cov.get((b, a), 0.0)) >= min_cover:
            par[find(a)] = find(b)
    comp = collections.defaultdict(list)
    for n in names:
        comp[find(n)].append(n)
    return sorted(max(comp.values(), key=lambda c: (len(c), sorted(c)[0] < "~")))


def gene_transcripts(gff, gene):
    """{transcript id: exons (0-based half-open)} of one gene from the tabix-indexed CAT GFF3."""
    import pysam
    tx = collections.defaultdict(list)
    tids = set()
    for row in gff.fetch(gene["chrom"], gene["start"], gene["end"]):
        f = row.split("\t")
        attrs = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        if f[2] == "transcript" and attrs.get("Parent") == gene["gene_id"]:
            tids.add(attrs["ID"])
        elif f[2] == "exon":
            tx[attrs.get("Parent")].append((int(f[3]) - 1, int(f[4])))
    return {t: sorted(tx[t]) for t in tids if tx.get(t)}


def spliced(fa, chrom, strand, exons):
    s = "".join(fa.fetch(chrom, a, b) for a, b in exons).upper()
    return s[::-1].translate(COMP) if strand == "-" else s


def build_family(name, pattern, genes_tsv, gff, fa, cap, seed, min_cover=MIN_COVER, link="tx"):
    import csv
    rx = re.compile(pattern)
    genes = []
    for r in csv.DictReader(open(genes_tsv), delimiter="\t"):
        if rx.search(r["gene_name"]) and int(r["end"]) - int(r["start0"]) <= MAX_SPAN:
            genes.append(dict(gene_id=r["gene_id"], gene_name=r["gene_name"], chrom=r["chrom"], start=int(r["start0"]),
                              end=int(r["end"]), strand=r["strand"], n=int(r["n_transcripts"])))
    genes = drop_overlaps(genes)
    for g in genes:
        g["tx"] = gene_transcripts(gff, g)
    genes = [g for g in genes if g["tx"]]
    if len(genes) < 2:
        return None
    seqs = {g["gene_id"]: spliced(fa, g["chrom"], g["strand"], g["tx"][canonical(g["tx"])]) for g in genes}
    with tempfile.TemporaryDirectory() as d:
        with open(f"{d}/c.fa", "w") as fh:
            for k, v in seqs.items():
                fh.write(f">{k}\n{v}\n")
        out = subprocess.run(["minimap2", "-cx", "asm20", "-X", f"{d}/c.fa", f"{d}/c.fa"], capture_output=True, text=True).stdout
    hits = [(f[0], f[5], int(f[9]) / int(f[10]), int(f[10]), int(f[1]), int(f[6])) for f in (ln.split("\t") for ln in out.splitlines())]
    if link == "locus":                                   # canonical transcripts only: a short alternative isoform that is
        cov = locus_cover([dict(g, tx={canonical(g["tx"]): g["tx"][canonical(g["tx"])]}) for g in genes], fa)  # one domain
        # would link genes sharing only that domain (LILR 2 -> 9 with any isoform)
        comp = locus_component(list(seqs), cov, min_cover=min_cover)
        n_comp = len(comp)
    else:
        comp = largest_component(list(seqs), hits, min_cover=min_cover)
        n_comp = len(comp)
    if len(comp) < 2:
        return None
    if len(comp) > cap:
        comp = sorted(random.Random(f"{seed}:{name}").sample(comp, cap))
    by = {g["gene_id"]: g for g in genes}
    copies = []
    for gid in sorted(comp, key=lambda x: (by[x]["chrom"], by[x]["start"])):
        g = by[gid]
        main, iso = pick_isoforms(g["tx"])
        copies.append(dict(gene_id=gid, gene_name=g["gene_name"], chrom=g["chrom"], start=min(s for v in iso.values() for s, _ in v),
                           end=max(e for v in iso.values() for _, e in v), strand=g["strand"], main=main,
                           transcripts={t: [list(x) for x in v] for t, v in iso.items()}))
    return dict(family=name, pattern=pattern, candidates=len(genes), component=n_comp,
                copies=copies)


def cmd_build(a):
    import pysam
    gff = pysam.TabixFile(f"{CAT}/chm13v2.0_CAT_Liftoff.slim.gff3.gz")
    fa = pysam.FastaFile(GENOME)
    out = {}
    for st, groups in GROUPS.items():
        out[st] = []
        for name, pat in groups.items():
            fam = build_family(name, pat, f"{CAT}/chm13v2.0_CAT_Liftoff.genes.tsv", gff, fa, a.cap, a.seed, a.min_cover, a.link)
            if fam:
                out[st].append(fam)
            print(f"{st} {name}: " + (f"{len(fam['copies'])} copies (component {fam['component']} of {fam['candidates']}), "
                                      f"{sum(len(c['transcripts']) for c in fam['copies'])} isoforms" if fam else "no family"), flush=True)
    json.dump(dict(genome=GENOME, annotation=f"{CAT}/chm13v2.0_CAT_Liftoff.slim.gff3.gz", cap=a.cap, seed=a.seed,
                   **({"min_cover": a.min_cover} if a.min_cover != MIN_COVER else {}),
                   **({"link": a.link} if a.link != "tx" else {}), sets=out),
              open(a.out, "w"), indent=1)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("build")
    p.add_argument("out")
    p.add_argument("--cap", type=int, default=10)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--min-cover", type=float, default=MIN_COVER,
                   help="a link must align over this fraction of the shorter canonical transcript (default 0.5)")
    p.add_argument("--link", choices=("tx", "locus"), default="tx",
                   help="tx: canonical transcripts aligned to each other (default); locus: each gene's canonical transcript "
                        "spliced onto the other gene's genomic locus, >= --min-cover of it in one direction (copies and fragments "
                        "link, genes sharing only a segment do not)")
    p.set_defaults(func=cmd_build)
    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
