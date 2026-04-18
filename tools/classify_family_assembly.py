#!/usr/bin/env python3
"""
classify_family_assembly.py — classify an assembly's recovery of annotated
GOLGA6/GOLGA8 transcripts.

Classes per annotated locus:
  exact       : overlapping same-strand transcript with identical intron chain
  partial     : overlapping same-strand transcript whose chain is a contiguous
                subchain of the annotated chain
  wrong-strand: overlapping transcript on the opposite strand
  other       : overlapping same-strand transcript, chain doesn't match
  missing     : no overlapping transcript

Usage:
    python3 classify_family_assembly.py <GFF> <pattern> <ASSEMBLY.gtf> [--tol=10]

The pattern is a description substring such as
    'golgin subfamily A member 6' or 'golgin subfamily A member 8'
"""
import argparse
import sys
from collections import defaultdict


def parse_attrs_gff(s):
    out = {}
    for kv in s.strip().rstrip(";").split(";"):
        kv = kv.strip()
        if not kv:
            continue
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip().strip('"')
    return out


def parse_attrs_gtf(s):
    out = {}
    for kv in s.strip().rstrip(";").split(";"):
        kv = kv.strip()
        if not kv:
            continue
        if " " in kv:
            k, v = kv.split(" ", 1)
            out[k.strip()] = v.strip().strip('"')
    return out


def matches(gene_desc, pattern):
    clauses = [p.strip() for p in pattern.split("||") if p.strip()]
    for c in clauses:
        if c.lower() in gene_desc.lower():
            return True
    return False


def load_gff_loci(gff_path, pattern):
    """Return [(locus_name, chrom, strand, intron_chain_sorted, start, end)...]

    For each gene matching pattern, take the canonical (longest-intron-chain)
    mRNA as the reference transcript.
    """
    gene_by_id = {}
    mrna_by_gene = defaultdict(list)
    exon_by_mrna = defaultdict(list)
    mrna_parent = {}

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            chrom, src, feat, s, e, _, strand, _, attr = p
            a = parse_attrs_gff(attr)
            if feat == "gene":
                gene_by_id[a.get("ID", "")] = (
                    a.get("Name", a.get("gene", "?")),
                    chrom,
                    strand,
                    int(s),
                    int(e),
                    a.get("description", ""),
                )
            elif feat == "mRNA":
                gid = a.get("Parent", "")
                mid = a.get("ID", "")
                mrna_parent[mid] = gid
                mrna_by_gene[gid].append(mid)
            elif feat == "exon":
                mid = a.get("Parent", "")
                exon_by_mrna[mid].append((int(s), int(e)))

    loci = []
    for gid, (name, chrom, strand, gs, ge, desc) in gene_by_id.items():
        if not matches(desc, pattern):
            continue
        best_mid = None
        best_introns = -1
        best_chain = None
        best_span = None
        for mid in mrna_by_gene.get(gid, []):
            exs = sorted(exon_by_mrna.get(mid, []))
            if len(exs) < 2:
                continue
            chain = tuple((exs[i][1] + 1, exs[i + 1][0] - 1) for i in range(len(exs) - 1))
            if len(chain) > best_introns:
                best_introns = len(chain)
                best_mid = mid
                best_chain = chain
                best_span = (exs[0][0], exs[-1][1])
        if best_mid is None:
            # zero introns; still add as a locus with no chain
            loci.append((name, chrom, strand, (), (gs, ge)))
        else:
            loci.append((name, chrom, strand, best_chain, best_span))
    return sorted(loci, key=lambda x: (x[1], x[4][0]))


def load_gtf_transcripts(gtf_path):
    """Return list of (tid, chrom, strand, intron_chain, start, end)."""
    exons_by_tid = defaultdict(list)
    meta = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            chrom, src, feat, s, e, _, strand, _, attr = p
            a = parse_attrs_gtf(attr)
            tid = a.get("transcript_id", "")
            if not tid:
                continue
            if feat == "transcript":
                meta[tid] = (chrom, strand)
            elif feat == "exon":
                exons_by_tid[tid].append((int(s), int(e)))
                if tid not in meta:
                    meta[tid] = (chrom, strand)

    out = []
    for tid, exs in exons_by_tid.items():
        chrom, strand = meta[tid]
        exs = sorted(exs)
        chain = tuple((exs[i][1] + 1, exs[i + 1][0] - 1) for i in range(len(exs) - 1))
        out.append((tid, chrom, strand, chain, (exs[0][0], exs[-1][1])))
    return out


def chains_overlap(tol, a, b):
    if len(a) != len(b):
        return False
    for (ad, aa), (bd, ba) in zip(a, b):
        if abs(ad - bd) > tol or abs(aa - ba) > tol:
            return False
    return True


def is_contiguous_subchain(sub, full, tol):
    if len(sub) > len(full):
        return False
    n = len(sub)
    for off in range(len(full) - n + 1):
        ok = True
        for i in range(n):
            ad, aa = sub[i]
            bd, ba = full[off + i]
            if abs(ad - bd) > tol or abs(aa - ba) > tol:
                ok = False
                break
        if ok:
            return True
    return False


def spans_overlap(span_a, span_b):
    return not (span_a[1] < span_b[0] or span_b[1] < span_a[0])


def classify_locus(locus, tx_by_chrom, tol):
    name, chrom, strand, chain, span = locus
    cands_same = []
    cands_opp = []
    for (tid, tchrom, tstrand, tchain, tspan) in tx_by_chrom.get(chrom, []):
        if not spans_overlap(span, tspan):
            continue
        if tstrand == strand:
            cands_same.append((tid, tchain, tspan))
        else:
            cands_opp.append((tid, tchain, tspan))

    # check exact (same-strand, identical chain)
    for tid, tchain, _ in cands_same:
        if chain and chains_overlap(tol, chain, tchain):
            return ("exact", tid)

    # partial: same-strand, contiguous subchain (query's chain is subchain of annot's OR vice versa)
    for tid, tchain, _ in cands_same:
        if tchain and chain and is_contiguous_subchain(tchain, chain, tol):
            return ("partial", tid)

    # other: same strand, overlapping, chain differs
    if cands_same:
        best = max(cands_same, key=lambda x: len(x[1]))
        return ("other", best[0])

    # wrong-strand: opposite strand overlap
    if cands_opp:
        best = max(cands_opp, key=lambda x: len(x[1]))
        return ("wrong-strand", best[0])

    return ("missing", "—")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gff")
    ap.add_argument("pattern")
    ap.add_argument("assembly_gtf")
    ap.add_argument("--tol", type=int, default=10)
    ap.add_argument("--label", default=None)
    args = ap.parse_args()

    loci = load_gff_loci(args.gff, args.pattern)
    txs = load_gtf_transcripts(args.assembly_gtf)
    tx_by_chrom = defaultdict(list)
    for t in txs:
        tx_by_chrom[t[1]].append(t)

    counts = defaultdict(int)
    rows = []
    for locus in loci:
        cls, best_tid = classify_locus(locus, tx_by_chrom, args.tol)
        counts[cls] += 1
        name, chrom, strand, chain, span = locus
        rows.append((name, chrom, len(chain), best_tid, cls))

    label = args.label or args.assembly_gtf
    print(f"# {label}")
    print(f"# pattern: {args.pattern}")
    print(f"# tol={args.tol}  n_loci={len(loci)}")
    print()
    print("| Locus | Chrom | Introns | Best tx | Class |")
    print("|---|---|---:|---|---|")
    for (name, chrom, n_int, tid, cls) in rows:
        print(f"| {name} | {chrom} | {n_int} | {tid} | {cls} |")
    print()
    print("| Class | Count |")
    print("|---|---:|")
    for k in ("exact", "partial", "wrong-strand", "other", "missing"):
        print(f"| {k} | {counts.get(k, 0)} |")


if __name__ == "__main__":
    main()
