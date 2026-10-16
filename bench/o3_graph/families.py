#!/usr/bin/env python3
"""Handpicked real gene families for the graph-O2/O3 simulation v2 (prereg docs/PREREG_o3_graph_real_2026-10-15.md): mainstay
multi-copy families (NPIP, TBC1D3, GOLGA, Soto's and others), their copies and curated isoforms from CHM13 v2.0's CAT/Liftoff
annotation (the project's human default).

    families.py build OUT.json [--cap 10] [--seed 1]

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
MAX_SPAN, MIN_IDENT, ISO_CAP = 60000, 0.90, 4
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


def largest_component(names, hits, min_ident=MIN_IDENT):
    """hits: (a, b, identity, block, len a, len b); linked at identity >= min_ident over >= half the shorter."""
    par = {n: n for n in names}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for a, b, ident, blk, la, lb in hits:
        if a != b and a in par and b in par and ident >= min_ident and blk >= 0.5 * min(la, lb):
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


def build_family(name, pattern, genes_tsv, gff, fa, cap, seed):
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
    comp = largest_component(list(seqs), hits)
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
    return dict(family=name, pattern=pattern, candidates=len(genes), component=len(largest_component(list(seqs), hits)),
                copies=copies)


def cmd_build(a):
    import pysam
    gff = pysam.TabixFile(f"{CAT}/chm13v2.0_CAT_Liftoff.slim.gff3.gz")
    fa = pysam.FastaFile(GENOME)
    out = {}
    for st, groups in GROUPS.items():
        out[st] = []
        for name, pat in groups.items():
            fam = build_family(name, pat, f"{CAT}/chm13v2.0_CAT_Liftoff.genes.tsv", gff, fa, a.cap, a.seed)
            if fam:
                out[st].append(fam)
            print(f"{st} {name}: " + (f"{len(fam['copies'])} copies (component {fam['component']} of {fam['candidates']}), "
                                      f"{sum(len(c['transcripts']) for c in fam['copies'])} isoforms" if fam else "no family"), flush=True)
    json.dump(dict(genome=GENOME, annotation=f"{CAT}/chm13v2.0_CAT_Liftoff.slim.gff3.gz", cap=a.cap, seed=a.seed, sets=out),
              open(a.out, "w"), indent=1)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("build")
    p.add_argument("out")
    p.add_argument("--cap", type=int, default=10)
    p.add_argument("--seed", type=int, default=1)
    p.set_defaults(func=cmd_build)
    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
