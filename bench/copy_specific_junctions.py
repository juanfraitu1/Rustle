#!/usr/bin/env python3
"""COPY-specific junctions: do paralog COPIES splice differently? The unification of the three
threads — (1) discovered cross-chromosome copies, (2) reads assigned to each copy (here by genomic
locus; the PSV-based copy_split assignment applies in the collapsed regime, data-limited in GGO per
task c), (3) per-copy junction usage compared between copies (the ASJ machinery, allele->copy).

For each confirmed cross-chrom recent-dup pair with BOTH copies expressed: align their exon structures
(intron-chain NW), and for each HOMOLOGOUS junction compare PSI between the two copies (Fisher exact);
also flag COPY-PRIVATE exons (an exon one copy includes that the other lacks). BH-FDR.

Run with /home/juanfra/miniforge3/bin/python (scipy). Deterministic.
"""
import csv
import glob
import os
import re
from collections import defaultdict

import pysam
from scipy.stats import fisher_exact

GW = "/tmp/gw"
PAIRS = os.path.join(os.path.dirname(__file__), "crosschrom_graded.tsv")
BAM = "/home/juanfra/winloci_scratch/GGO.bam"
OUT_TSV = os.path.join(os.path.dirname(__file__), "copy_specific_junctions.tsv")
OUT_MD = os.path.join(os.path.dirname(__file__), "copy_specific_junctions.md")
CORE_IDENT_MIN = 0.90
MIN_COV = 15        # both copies must have >= this many reads
MIN_EXON = 3
EX_BP, EX_FRAC = 6, 0.10
GAP, MISS = -1.0, -2.0
MIN_J = 3          # junction observed >= this in a copy
MIN_SPAN = 5       # >= spanning reads per copy
DPSI = 0.30


def gff_attrs(f):
    d = {}
    for kv in f.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1); d[k.strip()] = v.strip()
    return d


def load_gene_exons():
    """gene -> (chrom, strand, exons_transcription_order[(gstart,gend)])  (rep = longest transcript)."""
    gmeta = {}; feat_gene = {}; tx_exons = defaultdict(list)
    for gff in sorted(glob.glob(os.path.join(GW, "ref_*.gff3"))):
        chrom = os.path.basename(gff)[len("ref_"):-len(".gff3")]
        for line in open(gff):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = gff_attrs(f[8])
            if f[2] == "gene":
                g = a.get("gene") or (a.get("ID", "")[5:] if a.get("ID", "").startswith("gene-") else a.get("ID"))
                if g:
                    gmeta.setdefault(g, (chrom, f[6]))
            elif f[2] in ("mRNA", "transcript", "primary_transcript"):
                if a.get("ID") and a.get("gene"):
                    feat_gene[(chrom, a["ID"])] = a["gene"]
            elif f[2] == "exon":
                p = a.get("Parent", "")
                if (chrom, p) in feat_gene:
                    tx_exons[(chrom, p)].append((int(f[3]) - 1, int(f[4])))   # 0-based half-open
    best = {}
    for (chrom, fid), exons in tx_exons.items():
        g = feat_gene[(chrom, fid)]
        tot = sum(e - s for s, e in exons)
        if g not in best or tot > best[g][2]:
            best[g] = (chrom, sorted(exons), tot)
    out = {}
    for g, (chrom, exons, _) in best.items():
        strand = gmeta.get(g, (chrom, "+"))[1]
        ordered = exons if strand == "+" else exons[::-1]   # transcription order
        out[g] = (chrom, strand, ordered)
    return out


def reads_at(bam, chrom, lo, hi):
    """(introns set, rstart, rend) per primary read overlapping [lo,hi]."""
    out = []
    for aln in bam.fetch(chrom, lo, hi):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.query_sequence is None:
            continue
        introns = []
        rpos = aln.reference_start
        for op, length in aln.cigartuples:
            if op == 3:
                introns.append((rpos, rpos + length)); rpos += length
            elif op in (0, 2, 7, 8):
                rpos += length
        out.append((set(introns), aln.reference_start, aln.reference_end))
    return out


def psi(reads, d, a):
    used = span = 0
    for js, rs, re in reads:
        if rs <= d and re >= a:
            span += 1
            if (d, a) in js:
                used += 1
    return used, span


def genomic_intron(ex_i, ex_j):
    """intron between two transcription-consecutive exons (genomic coords donor<acceptor)."""
    left, right = (ex_i, ex_j) if ex_i[0] < ex_j[0] else (ex_j, ex_i)
    return left[1], right[0]


def align_exons(A, B):
    """NW over exon-length vectors; return matched (iA,iB) index pairs (similar lengths)."""
    LA = [e - s for s, e in A]; LB = [e - s for s, e in B]
    na, nb = len(LA), len(LB)
    dp = [[0.0] * (nb + 1) for _ in range(na + 1)]
    bt = [[None] * (nb + 1) for _ in range(na + 1)]
    for i in range(1, na + 1):
        dp[i][0] = i * GAP; bt[i][0] = "u"
    for j in range(1, nb + 1):
        dp[0][j] = j * GAP; bt[0][j] = "l"
    for i in range(1, na + 1):
        for j in range(1, nb + 1):
            ok = abs(LA[i-1] - LB[j-1]) <= max(EX_BP, EX_FRAC * max(LA[i-1], LB[j-1]))
            d = dp[i-1][j-1] + (1.0 if ok else MISS)
            u = dp[i-1][j] + GAP; l = dp[i][j-1] + GAP
            best, w = d, ("d" if ok else "x")
            if u > best: best, w = u, "u"
            if l > best: best, w = l, "l"
            dp[i][j] = best; bt[i][j] = w
    i, j = na, nb; pairs = []
    while i > 0 or j > 0:
        w = bt[i][j]
        if w in ("d", "x"):
            if w == "d":
                pairs.append((i-1, j-1))
            i -= 1; j -= 1
        elif w == "u":
            i -= 1
        else:
            j -= 1
    pairs.reverse()
    return pairs


def main():
    genes = load_gene_exons()
    bam = pysam.AlignmentFile(BAM, "rb")
    pairs = [r for r in csv.DictReader(open(PAIRS), delimiter="\t")
             if float(r["core_ident"]) >= CORE_IDENT_MIN]
    rows = []
    n_pairs_both_expr = 0
    cov_cache = {}
    def ncov(g):
        if g in cov_cache:
            return cov_cache[g]
        if g not in genes:
            cov_cache[g] = 0; return 0
        c, st, ex = genes[g]
        lo = min(s for s, _ in ex); hi = max(e for _, e in ex)
        n = sum(1 for a in bam.fetch(c, lo, hi)
                if not (a.is_unmapped or a.is_secondary or a.is_supplementary)) if hi > lo else 0
        cov_cache[g] = n; return n

    for r in pairs:
        ga, gb = r["gene_a"], r["gene_b"]
        if ga not in genes or gb not in genes:
            continue
        if len(genes[ga][2]) < MIN_EXON or len(genes[gb][2]) < MIN_EXON:
            continue
        if ncov(ga) < MIN_COV or ncov(gb) < MIN_COV:
            continue
        n_pairs_both_expr += 1
        ca, sa, exA = genes[ga]; cb, sb, exB = genes[gb]
        ra = reads_at(bam, ca, min(s for s, _ in exA), max(e for _, e in exA))
        rb = reads_at(bam, cb, min(s for s, _ in exB), max(e for _, e in exB))
        matched = align_exons(exA, exB)
        mset = set(matched)
        # homologous junctions: consecutive matched exon pairs (i,j) and (i+1,j+1)
        pair_rows = []
        maxA = maxB = 0.0
        for (ia, ib) in matched:
            if (ia + 1, ib + 1) not in mset:
                continue
            dA, aA = genomic_intron(exA[ia], exA[ia + 1])
            dB, aB = genomic_intron(exB[ib], exB[ib + 1])
            uA, spA = psi(ra, dA, aA); uB, spB = psi(rb, dB, aB)
            if spA < MIN_SPAN or spB < MIN_SPAN or (uA + uB) < MIN_J:
                continue
            pA, pB = uA / spA, uB / spB
            maxA = max(maxA, pA); maxB = max(maxB, pB)
            if (pA >= 0.98 and pB >= 0.98) or (pA <= 0.02 and pB <= 0.02):
                continue
            _, pv = fisher_exact([[uA, spA - uA], [uB, spB - uB]])
            pair_rows.append(dict(gene_a=ga, gene_b=gb, chrom_a=ca, chrom_b=cb,
                             jA=f"{dA}-{aA}", jB=f"{dB}-{aB}", uA=uA, spA=spA, uB=uB, spB=spB,
                             psiA=round(pA, 3), psiB=round(pB, 3), dPSI=round(abs(pA - pB), 3),
                             kind="homologous", pval=pv))
        # confound control: a copy that splices NOTHING (max PSI<0.5 over its homologous junctions)
        # is a retrocopy/unspliced copy -> its dPSI=1 junctions are trivial. Genuine DIFFERENTIAL
        # splicing requires BOTH copies to actually splice.
        both_spliced = maxA >= 0.5 and maxB >= 0.5
        for r in pair_rows:
            r["both_spliced"] = int(both_spliced)
        rows.extend(pair_rows)
        # copy-private exons: a multi-flanked exon in one copy with no matched counterpart
        matchedA = {ia for ia, _ in matched}; matchedB = {ib for _, ib in matched}
        for (ex, rd, lab, chrom, gA, gB) in [(exA, ra, "A_private", ca, ga, gb),
                                             (exB, rb, "B_private", cb, gb, ga)]:
            ms = matchedA if lab == "A_private" else matchedB
            for k in range(1, len(ex) - 1):
                if k in ms:
                    continue
                d1, a1 = genomic_intron(ex[k - 1], ex[k]); d2, a2 = genomic_intron(ex[k], ex[k + 1])
                u1, s1 = psi(rd, d1, a1)
                if s1 >= MIN_SPAN and u1 >= MIN_J and u1 / s1 >= 0.5:
                    rows.append(dict(gene_a=gA, gene_b=gB, chrom_a=chrom, chrom_b="-",
                                     jA=f"{d1}-{a1}", jB="-", uA=u1, spA=s1, uB=0, spB=0,
                                     psiA=round(u1/s1, 3), psiB=0.0, dPSI=round(u1/s1, 3),
                                     kind=lab, pval=1.0))

    # BH-FDR over homologous-junction tests only
    hom = [r for r in rows if r["kind"] == "homologous"]
    hom.sort(key=lambda r: r["pval"]); m = len(hom); prev = 1.0
    for i in range(m - 1, -1, -1):
        prev = min(prev, hom[i]["pval"] * m / (i + 1)); hom[i]["q"] = prev
    for r in rows:
        r.setdefault("q", 1.0)

    csj_all = [r for r in hom if r["q"] < 0.05 and r["dPSI"] >= DPSI]
    csj = [r for r in csj_all if r.get("both_spliced")]            # genuine differential splicing
    unspliced = [r for r in csj_all if not r.get("both_spliced")]  # one copy is a retrocopy/unspliced
    priv = [r for r in rows if r["kind"].endswith("private")]
    csj_pairs = {(r["gene_a"], r["gene_b"]) for r in csj}

    for r in rows:
        r.setdefault("both_spliced", "")
    with open(OUT_TSV, "w") as fh:
        cols = ["kind", "both_spliced", "gene_a", "chrom_a", "jA", "psiA", "gene_b", "chrom_b", "jB",
                "psiB", "uA", "spA", "uB", "spB", "dPSI", "pval", "q"]
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda r: (-r["dPSI"], r.get("q", 1.0))):
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    L = []
    def P(s=""): L.append(s)
    P("# Copy-specific junctions — do paralog COPIES splice differently?")
    P()
    P("The unification: (1) discovered cross-chromosome copies → (2) reads assigned per copy (by locus "
      "here; PSV-based copy_split applies in the collapsed regime, data-limited in GGO) → (3) per-copy "
      "junction usage compared between copies (the ASJ machinery, allele→copy).")
    P()
    P(f"- cross-chrom recent-dup pairs (core_ident≥{CORE_IDENT_MIN}, both ≥{MIN_EXON} exons, BOTH "
      f"copies expressed ≥{MIN_COV} reads): **{n_pairs_both_expr}**")
    P(f"- homologous junctions compared between copies: **{len(hom)}**")
    P(f"- copy-specific junctions raw (q<0.05, |dPSI|≥{DPSI}): **{len(csj_all)}**, of which:")
    P(f"  - **DIFFERENTIAL splicing (both copies splice): {len(csj)}** across **{len(csj_pairs)}** copy pairs")
    P(f"  - one copy UNSPLICED/retrocopy (trivial dPSI): {len(unspliced)}")
    P(f"- copy-private exon junctions (one copy includes an exon the other lacks): **{len(priv)}**")
    P()
    P("## Top DIFFERENTIAL copy-specific junctions (both copies splice, differ at this junction)")
    P("| copy A | jA | PSI_A | copy B | jB | PSI_B | dPSI | q |")
    P("|---|---|---|---|---|---|---|---|")
    for r in sorted(csj, key=lambda r: (-r["dPSI"], r["q"]))[:20]:
        P(f"| {r['gene_a']} | {r['jA']} | {r['psiA']} | {r['gene_b']} | {r['jB']} | {r['psiB']} | "
          f"{r['dPSI']:.2f} | {r['q']:.1e} |")
    P()
    P("## Honest notes")
    P("- BOTH copies expressed is rare (~3% of recent-dup pairs) — most paralogs have one dominant "
      "copy (consistent with the headroom findings). So this is a small, real set, not a large catalog.")
    P("- read→copy assignment here is by genomic locus (cross-chrom copies map to distinct coords). The "
      "PSV-based copy_split assignment (interest #2) is exercised only in the COLLAPSED regime, which is "
      "essentially empty in GGO (task c) — it needs a deep co-located dataset (testis HiFi for DAZ/RBMY).")
    P("- homologous-junction mapping is via exon-length alignment; copy-private exons flag structure "
      "divergence (a copy gained/lost an exon).")
    P("- coverage-limited; deeper data finds more both-expressed pairs.")
    P()
    P("## Reproduce")
    P("- `MINIFORGE python bench/copy_specific_junctions.py`")
    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"\n[wrote {OUT_TSV} and {OUT_MD}]")


if __name__ == "__main__":
    main()
