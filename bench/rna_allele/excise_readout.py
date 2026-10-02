#!/usr/bin/env python3
"""Readouts R1-R3 of Amendment 4 (docs/PREREG_rna_allele_haplotype_count_2026-10-01.md): NPIPA2 hard-masked, the NPIP-net reads realigned
to the masked and the unmasked genome with the same minimap2.

R1 fate of NPIPA2's reads (baseline primary on its exons): landing copy, concentration, divergence (`de`) before and after.
R2 allele sites at the landing copy, baseline vs masked: the registered site test on its primaries of reads that are not AS-tied within
   this BAM (second AS < 0.98 x best over the read's records); PSV filter over the paralog hits that remain in the reference (hits on the
   masked span dropped). A site only in the masked arm whose minor base equals NPIPA2's base at that column = a fake allele.
   (O2-assigned tied reads are not added: O2 was not re-run on the subset; disclosed.)
R3 no-reference-match groups at the landing copy (the caller's novel-haplotype rule, NPIPA2 removed from the paralog list).

    python3 excise_readout.py --work /mnt/linuxdisk/tmp/rna_allele --x /mnt/linuxdisk/tmp/rna_allele/excise
"""
import argparse
import collections
import csv
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import caller as C  # noqa: E402  (Paralogs, homopolymer_near, phase, cigar_ops, blocks)

MASK = ("NC_073242.2", 32426793, 32456482)
COMP = str.maketrans("ACGTN", "TGCAN")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True)
    ap.add_argument("--x", required=True)
    a = ap.parse_args(argv)
    W, X = a.work, a.x
    genes = {r["gene_id"]: r for r in csv.DictReader(open(f"{W}/genes.tsv"), delimiter="\t") if r["biotype"] == "S_fam:NPIP"}
    ex = {g: C.blocks(r["exons"]) for g, r in genes.items()}
    name = {g: r["name"] for g, r in genes.items()}

    def on_copy(rd):
        out = []
        for g, r in genes.items():
            if r["chrom"] != rd.reference_name:
                continue
            bl = ex[g]
            if rd.reference_start < bl[-1][1] and bl[0][0] < rd.reference_end:
                if any(min(y, e) > max(x, s) for x, y in bl for s, e in rd.get_blocks()):
                    out.append(g)
        return out
    prim = {}
    for arm in ("base", "masked"):
        bam = pysam.AlignmentFile(f"{X}/{arm}.bam")
        d = {}
        for rd in bam.fetch(until_eof=True):
            if rd.is_unmapped:
                d.setdefault(rd.query_name, ("unmapped", [], None))
                continue
            if rd.is_secondary or rd.is_supplementary:
                continue
            d[rd.query_name] = (rd.reference_name, on_copy(rd), rd.get_tag("de") if rd.has_tag("de") else None)
        prim[arm] = d
    A = [q for q, (c, cps, de) in prim["base"].items() if "gN15" in cps]
    land = collections.Counter()
    de_b, de_m = [], []
    for q in A:
        c, cps, de = prim["masked"].get(q, ("missing", [], None))
        land[name[cps[0]] if cps else ("unmapped" if c == "unmapped" else f"elsewhere:{c}")] += 1
        if prim["base"][q][2] is not None:
            de_b.append(prim["base"][q][2])
        if de is not None:
            de_m.append(de)
    med = lambda v: sorted(v)[len(v) // 2] if v else float("nan")
    top, n_top = land.most_common(1)[0]
    print(f"R1: NPIPA2 reads (baseline primary on its exons) {len(A)}; after masking: {dict(land.most_common())}")
    print(f"    concentration on {top}: {n_top / len(A):.3f}; median de before {med(de_b):.4f} after {med(de_m):.4f}")
    L = next(g for g in genes if name[g] == top)
    r = genes[L]
    bl = ex[L]
    lo, hi = bl[0][0], bl[-1][1]
    fa = pysam.FastaFile("/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta")
    ref = fa.fetch(r["chrom"], max(0, lo - 20), hi + 20).upper()
    roff = max(0, lo - 20)
    # paralog hits of L's exon sum (registered definition), split into NPIPA2's hit(s) and the rest
    hits_all, hits_rest, hits_a2, qlen = [], [], [], None
    for ln in open(f"{W}/para.paf"):
        f = ln.rstrip("\n").split("\t")
        if f[0] != L:
            continue
        qlen = int(f[1])
        ident, cov = int(f[9]) / max(1, int(f[10])), (int(f[3]) - int(f[2])) / max(1, int(f[1]))
        own = f[5] == r["chrom"] and int(f[7]) < hi and lo < int(f[8])
        if own or ident < 0.90 or cov < 0.50:
            continue
        cg = next(x[5:] for x in f[12:] if x.startswith("cg:Z:"))
        h = dict(t=f[5], ts=int(f[7]), strand=f[4], qs=int(f[2]), qe=int(f[3]), ops=C.cigar_ops(cg))
        hits_all.append(h)
        masked_hit = f[5] == MASK[0] and int(f[7]) < MASK[2] and MASK[1] < int(f[8])
        (hits_a2 if masked_hit else hits_rest).append(h)
    P_rest, P_a2 = C.Paralogs(hits_rest, fa), C.Paralogs(hits_a2, fa)
    exlen = sum(y - x for x, y in bl)

    def esum_index(p):
        acc = 0
        for x, y in bl:
            if x <= p < y:
                i = acc + (p - x)
                return i if r["strand"] == "+" else exlen - 1 - i
            acc += y - x
        return None
    print(f"    landing copy {top} ({L}, {r['chrom']}:{lo}-{hi} {r['strand']}): paralog hits {len(hits_all)}, of which on NPIPA2 {len(hits_a2)}")
    for arm in ("base", "masked"):
        bam = pysam.AlignmentFile(f"{X}/{arm}.bam")
        # AS ties within this BAM
        best = collections.defaultdict(list)
        for rd in bam.fetch(until_eof=True):
            if not rd.is_unmapped and rd.has_tag("AS"):
                best[rd.query_name].append(rd.get_tag("AS"))
        tied = {q for q, v in best.items() if len(v) > 1 and sorted(v)[-2] >= 0.98 * max(v) and sorted(v)[-2] > 0}
        keep = lambda rd: not (rd.is_unmapped or rd.is_secondary or rd.is_supplementary) and rd.query_name not in tied
        cnt = bam.count_coverage(r["chrom"], lo, hi, quality_threshold=0, read_callback=keep)
        sites = []
        for x, y in bl:
            for p in range(x, y):
                c = [cnt[k][p - lo] for k in range(4)]
                tot = sum(c)
                if tot < 10:
                    continue
                order = sorted(range(4), key=lambda k: -c[k])
                if c[order[1]] < 3 or c[order[1]] < 0.20 * tot:
                    continue
                rb = ref[p - roff]
                b1, b2 = "ACGT"[order[0]], "ACGT"[order[1]]
                alt = ({b1, b2} - {rb}).pop() if rb in (b1, b2) else b2
                tr = (lambda b: b if r["strand"] == "+" else b.translate(COMP))
                if tr(rb) == "A" and tr(alt) == "G":
                    continue
                if C.homopolymer_near(ref, p, roff):
                    continue
                i = esum_index(p)
                gb = tr(rb)
                if any(b != "?" and b != gb for b in P_rest.bases_at(i, qlen)):
                    continue
                a2 = [b for b in P_a2.bases_at(i, qlen) if b not in ("?", None)]
                fake = bool(a2) and tr(alt) in a2 and a2[0] != gb
                sites.append((p, b1, b2, tot, c[order[1]], fake))
        nreads = sum(1 for rd in bam.fetch(r["chrom"], lo, hi) if keep(rd))
        print(f"R2 [{arm}] {top}: untied primaries {nreads}; allele sites {len(sites)}; of which NPIPA2's base (fake allele) "
              f"{sum(1 for s in sites if s[5])}; call {'2' if sites else '1+'}")
        for s in sites[:12]:
            print(f"      {r['chrom']}:{s[0]} {s[1]}/{s[2]} cov {s[3]} minor {s[4]} {'NPIPA2 base' if s[5] else ''}")
        # R3: no-reference-match groups at L, NPIPA2 removed from the paralog list
        psv = []
        for x, y in bl:
            for p in range(x, y):
                i = esum_index(p)
                rb = ref[p - roff]
                gb = rb if r["strand"] == "+" else rb.translate(COMP)
                pb = P_rest.bases_at(i, qlen)
                if any(b != "?" and b != gb for b in pb):
                    psv.append((p, gb, pb))
        pidx = {p: k for k, (p, _, _) in enumerate(psv)}
        groups = collections.Counter()
        members = collections.defaultdict(list)
        for rd in bam.fetch(r["chrom"], lo, hi):
            if rd.is_unmapped or rd.is_secondary or rd.is_supplementary:
                continue
            pat = {}
            for qp, rp in rd.get_aligned_pairs(matches_only=True):
                k = pidx.get(rp)
                if k is not None:
                    b = rd.query_sequence[qp]
                    pat[k] = b if r["strand"] == "+" else b.translate(COMP)
            if len(pat) < 2 or sum(1 for k, b in pat.items() if b != psv[k][1]) < 2:
                continue
            ok = True
            for hi_ in range(len(hits_rest)):
                cols = [(k, b) for k, b in pat.items() if psv[k][2][hi_] not in ("?", None)]
                if sum(1 for k, b in cols if b != psv[k][2][hi_]) < 2:
                    ok = False; break
            if ok:
                key = tuple(sorted(pat.items()))
                groups[key] += 1; members[key].append(rd.query_name)
        g2 = {k: v for k, v in groups.items() if v >= 2}
        print(f"R3 [{arm}] {top}: PSV columns {len(psv)}; no-reference-match groups (>= 2 reads) {len(g2)} holding "
              f"{sum(g2.values())} reads, of which NPIPA2 reads {sum(1 for k in g2 for q in members[k] if q in set(A))}")


if __name__ == "__main__":
    main()
