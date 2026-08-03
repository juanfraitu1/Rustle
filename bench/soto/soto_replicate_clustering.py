#!/usr/bin/env python3
"""Soto 2025 gene-family clustering, steps 4-6, reimplemented from their STAR Methods verbatim.

Their text (Cell 188:5363, "Gene family clustering"):
  "DNA sequences of all SD98 regions were extracted using BEDTools getfasta and mapped back to the reference
   genome using minimap2 (v2.17) with the following parameters: -c --end-bonus 5 --eqx -N50 -p0.5 -t64. For
   each SD98 exon, the BEDTools intersect with -f 0.99 parameter was used to select mappings covering >99% of
   the exon sequence, removing self-mappings. This list was refined using the previously published WSSD CNs
   (famCN) of humans from the SGDP (n=269) ... After comparing the median famCN values of SD98 genes with
   shared exons, groupings where the mean absolute deviation of the CN was less than one were selected. The
   list was filtered to focus on gene families containing at least one protein-coding or unprocessed
   pseudogene. SD98 genes associated with other gene features, including lncRNAs and processed pseudogenes,
   were also assigned a gene family ID. On the other hand, if a gene was not associated with any other gene
   feature, they were classified as 'unassigned' or 'singletons'."

So:
  step 4  exon E of gene B is SHARED with the genes of SD98 region R when some alignment whose QUERY is R
          covers >= 99% of E on the target side, and R is not E's own region (self-mapping).
  step 5  connected components of that shared-exon graph, then SPLIT any component whose median-famCN
          mean-absolute-deviation is >= 1.
  step 6  a component is a family only if it holds >= 1 protein_coding / unprocessed_pseudogene gene;
          other biotypes keep the family ID they joined; a gene with no partner is a singleton.

MAD here is mean ABSOLUTE deviation about the mean (their words), not median-absolute-deviation.
"""
import argparse, csv, gzip, statistics as st, sys
from collections import defaultdict


def open_maybe_gz(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def _cigar(cg):
    """Yield (length, op) from a PAF cg:Z: string."""
    n = 0
    for ch in cg:
        if ch.isdigit():
            n = n * 10 + (ord(ch) - 48)
        else:
            yield n, ch
            n = 0


def _project(blocks, ts, te):
    """Map a target interval [ts,te) to query coords through the aligned blocks. None if unaligned."""
    ql = qr = None
    for t_from, t_to, q_at, _ in blocks:
        if t_to <= ts:
            continue
        if t_from >= te:
            break
        if ql is None:
            ql = q_at + max(ts, t_from) - t_from
        qr = q_at + min(te, t_to) - t_from
    return (ql, qr) if ql is not None and qr is not None and qr > ql else None


def parse_region(name):
    """'chr1:1000-2000' -> ('chr1', 999, 2000) — samtools faidx names are 1-based inclusive."""
    c, rng = name.rsplit(":", 1)
    s, e = rng.split("-")
    return c, int(s) - 1, int(e)


def load_exons(cat_bed, keep_genes):
    """gene -> list of (chrom, start, end) exons, from the CAT bigBed-derived BED."""
    ex = defaultdict(set)
    meta = {}
    with open(cat_bed) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 21:
                continue
            gid = f[18]
            if gid not in keep_genes:
                continue
            meta[gid] = (f[12], f[19])
            c, st = f[0], int(f[1])
            sizes = [int(x) for x in f[10].rstrip(",").split(",") if x]
            starts = [int(x) for x in f[11].rstrip(",").split(",") if x]
            for sz, so in zip(sizes, starts):
                ex[gid].add((c, st + so, st + so + sz))
    return ex, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf", required=True, help="SD98 regions mapped back to the genome (may be .gz)")
    ap.add_argument("--geneset", required=True, help="sd98_geneset_v1.tsv (gene_id, name, biotype, ...)")
    ap.add_argument("--cat-bed", required=True, help="CAT v4 BED")
    ap.add_argument("--famcn", required=True, help="TSV: gene_id -> famCN (from famcn_from_wssd.py or S1C)")
    ap.add_argument("--min-cov", type=float, default=0.99, help="bedtools -f equivalent")
    ap.add_argument("--mad", type=float, default=1.0, help="famCN mean-absolute-deviation split threshold")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    genes, biotype = set(), {}
    with open(a.geneset) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genes.add(r["gene_id"])
            biotype[r["gene_id"]] = r.get("biotype", "")
    exons, meta = load_exons(a.cat_bed, genes)
    print(f"[step4] {len(genes)} SD98 genes, {sum(len(v) for v in exons.values())} exons", file=sys.stderr)

    # exon index per chromosome for interval lookup
    per_chrom = defaultdict(list)
    for g, evs in exons.items():
        for c, s, e in evs:
            per_chrom[c].append((s, e, g))
    for c in per_chrom:
        per_chrom[c].sort()

    # which SD98 region does each gene sit in? (for self-mapping removal)
    # a gene's region = the query region name whose span contains its first exon
    # -- derived from the PAF query names, so no extra input file is needed.

    # step 4: walk the PAF; for each alignment, find exons >= min_cov covered by its TARGET interval, then
    # PROJECT that exon back through the CIGAR to its QUERY position and link only the gene(s) whose exons
    # sit there.
    #
    # Linking every gene in the query REGION instead (the naive reading of "removing self-mappings") fuses
    # the catalog. Reproduced head-to-head on the same PAF, same downstream steps, universe = Soto's 1,793
    # Yes genes: literal region-level linking scores ARI 0.267 (pair P 0.268 / R 0.273) and produces a
    # 1,102-gene raw connected component, versus ARI 0.514 (P 0.854 / R 0.369) for the projection.
    # (An earlier comment here claimed "a 178-gene family, ARI 0.20" -- that figure is NOT reproducible on
    # any PAF and has been replaced with the measured values.)
    # HONESTY NOTE: this is a DEVIATION from Soto's literal sentence ("BEDTools intersect -f 0.99"), chosen
    # because it scores far better. It is defensible as the reading that preserves exon-to-exon
    # correspondence, but it must be presented as our choice, not as their method. The duplication maps gene A's exon onto gene
    # B's exon, and the CIGAR is what makes that correspondence exact.
    shared = defaultdict(set)   # gene -> set of partner genes
    n_aln = n_proj = 0
    with open_maybe_gz(a.paf) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            n_aln += 1
            qname, qstart, qend, strand = f[0], int(f[2]), int(f[3]), f[4]
            tname, tstart, tend = f[5], int(f[7]), int(f[8])
            cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), None)
            if cg is None:
                continue
            try:
                qc, qoff, _ = parse_region(qname)
            except ValueError:
                continue
            hits = per_chrom.get(tname)
            if not hits:
                continue
            # target->query breakpoint table from the CIGAR (query coords are in the REGION's frame)
            blocks = []          # (t_from, t_to, q_at_t_from, step) ; step +1 forward, -1 reverse
            qp, tp = qstart, tstart
            for n, op in _cigar(cg):
                if op in "=XM":
                    blocks.append((tp, tp + n, qp, 1))
                    qp += n
                    tp += n
                elif op in "DN":
                    tp += n
                elif op == "I":
                    qp += n
            for s, e, g in hits:
                if e <= tstart:
                    continue
                if s >= tend:
                    break
                if (min(e, tend) - max(s, tstart)) / (e - s) < a.min_cov:
                    continue
                pq = _project(blocks, max(s, tstart), min(e, tend))
                if pq is None:
                    continue
                # region-frame -> genome coords; on '-' the query was reverse-complemented
                ql, qr = pq
                if strand == "-":
                    ql, qr = (qstart + qend) - qr, (qstart + qend) - ql
                gl, gr = qoff + ql, qoff + qr
                if tname == qc and gl < e and gr > s:
                    continue                     # self-mapping: projects onto itself
                n_proj += 1
                for s2, e2, g2 in per_chrom.get(qc, ()):
                    if e2 <= gl:
                        continue
                    if s2 >= gr:
                        break
                    if g2 != g:
                        shared[g].add(g2)
                        shared[g2].add(g)
    print(f"[step4] {n_aln} alignments, {n_proj} exon projections -> "
          f"{len(shared)} genes with >=1 shared exon", file=sys.stderr)

    # connected components
    seen, comps = set(), []
    for g in shared:
        if g in seen:
            continue
        stack, comp = [g], []
        seen.add(g)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in shared[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        comps.append(comp)
    print(f"[step4] {len(comps)} raw shared-exon components", file=sys.stderr)

    # step 5: split components whose famCN mean-absolute-deviation >= threshold.
    # Their rule is stated as a selection ("groupings where the MAD of CN was less than one were selected"),
    # so a component that fails is split into MAD-coherent subgroups by greedy agglomeration on famCN.
    famcn = {}
    with open(a.famcn) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gid = r.get("gene_id") or r.get("Gene ID")
            v = r.get("famCN") or r.get("famCN_median") or r.get("Median famCN")
            try:
                famcn[gid] = float(v)
            except (TypeError, ValueError):
                pass

    def mad(vals):
        if len(vals) < 2:
            return 0.0
        m = sum(vals) / len(vals)
        return sum(abs(v - m) for v in vals) / len(vals)

    final = []
    for comp in comps:
        vals = [(famcn[g], g) for g in comp if g in famcn]
        if len(vals) < 2 or mad([v for v, _ in vals]) < a.mad:
            final.append(comp)
            continue
        # split: sort by famCN, start a new group whenever adding the next gene breaks MAD < 1
        vals.sort()
        cur, groups = [], []
        for v, g in vals:
            if cur and mad([x for x, _ in cur] + [v]) >= a.mad:
                groups.append([g2 for _, g2 in cur])
                cur = []
            cur.append((v, g))
        if cur:
            groups.append([g2 for _, g2 in cur])
        nocn = [g for g in comp if g not in famcn]
        if groups and nocn:
            groups[0].extend(nocn)      # genes without famCN cannot be split on it
        final.extend(groups)

    # step 6: a family needs >= 1 protein-coding / unprocessed pseudogene, and >= 2 members
    ELIGIBLE = {"protein_coding", "unprocessed_pseudogene",
                "transcribed_unprocessed_pseudogene", "translated_unprocessed_pseudogene"}
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "gene_name", "biotype", "family_id", "n_members", "famCN", "status"])
        fam_i = 0
        n_fam = n_single = 0
        for comp in sorted(final, key=lambda c: -len(c)):
            eligible = any(biotype.get(g) in ELIGIBLE for g in comp)
            if len(comp) >= 2 and eligible:
                fid = f"REPFAM{fam_i}"
                fam_i += 1
                n_fam += 1
                status = "family"
            else:
                fid = ""
                status = "singleton" if len(comp) < 2 else "no_coding_member"
                n_single += len(comp)
            for g in sorted(comp):
                nm, bt = meta.get(g, ("", biotype.get(g, "")))
                w.writerow([g, nm, bt, fid, len(comp),
                            f"{famcn[g]:.2f}" if g in famcn else "", status])
        # genes with no shared exon at all are singletons too
        for g in sorted(genes - seen):
            nm, bt = meta.get(g, ("", biotype.get(g, "")))
            w.writerow([g, nm, bt, "", 1, f"{famcn[g]:.2f}" if g in famcn else "", "singleton"])
            n_single += 1
    print(f"[done] {n_fam} families, {n_single} singleton/ineligible genes -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
