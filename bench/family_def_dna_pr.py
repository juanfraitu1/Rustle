#!/usr/bin/env python3
"""
family_def_dna_pr.py — precision / recall of the RNA read-conflict family definition
against an INDEPENDENT DNA-sequence ground truth.

DNA ground truth (the "biological" definition the RNA method is contrasted with):
  all-vs-all alignment of the longest cDNA per annotated gene (rep_cdna.fa, minimap2
  asm20).  For each unordered gene pair we keep the best identity and the aligned
  fraction in EACH direction (cov_a = how much of A aligns, cov_b = of B).
    - LOOSE paralog edge   : identity >= 0.90 AND max(cov_a,cov_b) >= 0.30   (the
      project's pinned bar, config.sh MIN_IDENTITY/MIN_COV_FRAC; UTRs of real paralogs
      diverge so a one-directional 30% bar keeps e.g. RABL2A/B at cov 0.34).
    - WHOLE-GENE paralog   : identity >= 0.90 AND min(cov_a,cov_b) >= 0.50   (reciprocal;
      excludes DOMAIN-SHARERS that match over only one shared domain/exon).
  Family = connected component (>=2 genes).  Independent of reads.

RNA definition (ours): the shipped de-tie read-conflict graph over the SAME 34k annotated
  gene vertices (family_def_genomewide.scan/pair_evidence/components, Δ=0.005).

Outputs precision/recall (edge level) for the loose, whole-gene, and expressed-restricted
DNA truths; decomposes the recall shortfall through the funnel expressed->cross-map->tied;
audits RNA-only edges by homology AND chromosome co-location; the panel confusion; and a
per-edge dump (family_def_dna_pr_edges.tsv) for independent verification.

Upstream (one-time, writes to /home/juanfra/winloci_scratch):
  1. cDNA per gene:   gffread -w GGO_cdna.fa -g GGO.fasta GGO_tx.gff
  2. longest per gene: python bench/build_rep_cdna.py            -> rep_cdna.fa
  3. all-vs-all:      minimap2 -cx asm20 -N 50 -p 0.1 -t 8 rep_cdna.fa rep_cdna.fa
                        | awk '...de:f...' -> rep_ava.tsv   (qgene tgene qlen tlen matches blocklen qstart qend strand mapq de)

Run: /home/juanfra/miniforge3/bin/python bench/family_def_dna_pr.py
"""
import collections
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_genomewide import (
    best_gene, ref_span, de_of, pair_evidence, components,
    GENES_BED, BAM, DELTA, DE_MAX, MIN_READS,
)

AVA = "/home/juanfra/winloci_scratch/rep_ava.tsv"
MIN_ID = 0.90
MIN_COV = 0.30       # loose (one-directional) bar
RECIP_COV = 0.50     # whole-gene (reciprocal) bar
EXPR_MIN = 1
EDGE_DUMP = os.path.join(HERE, "family_def_dna_pr_edges.tsv")


def load_genes():
    by_chrom = collections.defaultdict(list)
    g2pos = {}
    with open(GENES_BED) as f:
        for line in f:
            c, s, e, name = line.rstrip("\n").split("\t")
            by_chrom[c].append((int(s), int(e), name))
            if name not in g2pos:
                g2pos[name] = (c, int(s), int(e))
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom, g2pos


def parse_homology():
    """unordered gene-pair -> dict(id=best_identity, cov_a, cov_b) where cov_a is the max
    aligned fraction with the lexicographically-FIRST gene as query, cov_b the second."""
    H = {}
    n_rec = 0
    with open(AVA) as f:
        for line in f:
            q, t, ql, tl, m, bl, qs, qe, strand, mapq, de = line.rstrip("\n").split("\t")
            if q == t:
                continue
            n_rec += 1
            ql = int(ql); m = int(m); bl = int(bl); qs = int(qs); qe = int(qe)
            if bl == 0 or ql == 0:
                continue
            ident = m / bl
            cov = (qe - qs) / ql
            a, b = (q, t) if q < t else (t, q)
            key = (a, b)
            rec = H.get(key)
            if rec is None:
                rec = {"id": 0.0, "cov_a": 0.0, "cov_b": 0.0}
                H[key] = rec
            if ident > rec["id"]:
                rec["id"] = ident
            if q == a:
                rec["cov_a"] = max(rec["cov_a"], cov)
            else:
                rec["cov_b"] = max(rec["cov_b"], cov)
    return H, n_rec


def scan_rna(by_chrom):
    mm = collections.defaultdict(dict)
    gene_reads = collections.Counter()
    p = subprocess.Popen(["samtools", "view", "-f", "0x100", BAM],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        f = line.split("\t")
        if len(f) < 9:
            continue
        de = de_of(f[11:])
        if de is None or de > DE_MAX:
            continue
        g = best_gene(by_chrom, f[2], int(f[3]) - 1, int(f[3]) - 1 + ref_span(f[5]))
        if g is None:
            continue
        d = mm[f[0]]
        if g not in d or de < d[g]:
            d[g] = de
    p.wait()
    p = subprocess.Popen(["samtools", "view", "-F", "0x900", BAM],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        f = line.split("\t")
        if len(f) < 9:
            continue
        g = best_gene(by_chrom, f[2], int(f[3]) - 1, int(f[3]) - 1 + ref_span(f[5]))
        if g is None:
            continue
        gene_reads[g] += 1
        if f[0] in mm:
            de = de_of(f[11:])
            if de is not None and de <= DE_MAX:
                d = mm[f[0]]
                if g not in d or de < d[g]:
                    d[g] = de
    p.wait()
    return mm, gene_reads


def uf_components(edge_keys):
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        r = x
        while parent[r] != r:
            r = parent[r]
        while parent[x] != r:
            parent[x], x = r, parent[x]
        return r
    for a, b in edge_keys:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    comp = collections.defaultdict(set)
    for g in parent:
        comp[find(g)].add(g)
    return [c for c in comp.values() if len(c) >= 2]


def pr(rna, dna):
    inter = rna & dna
    p_ = len(inter) / max(len(rna), 1)
    r_ = len(inter) / max(len(dna), 1)
    return len(inter), len(rna - dna), len(dna - rna), p_, r_


def main():
    by_chrom, g2pos = load_genes()
    print(f"[vertices] {sum(len(v) for v in by_chrom.values()):,} gene intervals", flush=True)

    print("[DNA] parsing all-vs-all cDNA homology ...", flush=True)
    H, n_rec = parse_homology()
    dna_loose = {(a, b) for (a, b), r in H.items()
                 if r["id"] >= MIN_ID and max(r["cov_a"], r["cov_b"]) >= MIN_COV}
    dna_whole = {(a, b) for (a, b), r in H.items()
                 if r["id"] >= MIN_ID and min(r["cov_a"], r["cov_b"]) >= RECIP_COV}
    print(f"  {n_rec:,} cross-gene alignments", flush=True)
    print(f"  LOOSE  paralog edges (id>=.90, max-cov>=.30): {len(dna_loose):,} "
          f"-> {len(uf_components(dna_loose)):,} families", flush=True)
    print(f"  WHOLE  paralog edges (id>=.90, recip-cov>=.50): {len(dna_whole):,} "
          f"-> {len(uf_components(dna_whole)):,} families  (domain-sharers removed)", flush=True)

    print("[RNA] scanning GGO.bam (2 passes) ...", flush=True)
    mm, gene_reads = scan_rna(by_chrom)
    ev = pair_evidence(mm)
    rna_edge_list, rna_fams = components(ev, DELTA, DE_MAX, MIN_READS)
    rna = {((a, b) if a < b else (b, a)) for a, b, n in rna_edge_list}
    n_tied = {((a, b) if a < b else (b, a)): n for a, b, n in rna_edge_list}
    expressed = {g for g, n in gene_reads.items() if n >= EXPR_MIN}
    print(f"  {len(ev):,} cross-mapping pairs -> {len(rna):,} de-tie edges -> "
          f"{len(rna_fams):,} RNA families; expressed genes: {len(expressed):,}", flush=True)

    # expressed-restricted DNA truths (fair universe: RNA can only see expressed genes)
    dna_loose_expr = {(a, b) for (a, b) in dna_loose if a in expressed and b in expressed}
    dna_whole_expr = {(a, b) for (a, b) in dna_whole if a in expressed and b in expressed}

    def report(name, dna):
        tp, fp, fn, p_, r_ = pr(rna, dna)
        print(f"  {name:34} P={p_:.3f} R={r_:.3f}   TP={tp:5d} FP={fp:5d} FN={fn:6d}")
        return dict(name=name, tp=tp, fp=fp, fn=fn, precision=p_, recall=r_, dna=len(dna))

    print("\n" + "=" * 80)
    print("EDGE-LEVEL precision / recall   (RNA de-tie graph  vs  DNA paralogy graph)")
    print("=" * 80)
    R = {}
    R["loose"] = report("vs LOOSE DNA (all genes)", dna_loose)
    R["whole"] = report("vs WHOLE-GENE DNA (all genes)", dna_whole)
    R["loose_expr"] = report("vs LOOSE DNA (expressed only)", dna_loose_expr)
    R["whole_expr"] = report("vs WHOLE-GENE DNA (expressed only)", dna_whole_expr)

    # ---- PRECISION audit: RNA-only edges by homology AND co-location ----
    print("\n--- PRECISION audit: RNA-only edges (vs LOOSE DNA) by homology x co-location ---")
    rna_only = rna - dna_loose
    cat = collections.Counter()
    coloc = collections.Counter()
    for (a, b) in rna_only:
        r = H.get((a, b))
        if r is None:
            hc = "no_cDNA_homology"
        elif r["id"] >= 0.80:
            hc = "real_homology_id>=.80 (sub-bar)"
        else:
            hc = "weak_id<.80"
        cat[hc] += 1
        ca = g2pos.get(a, (None,))[0]; cb = g2pos.get(b, (None,))[0]
        if ca and cb and ca == cb:
            d = abs(g2pos[a][1] - g2pos[b][1])
            loc = "same_chrom_tandem<1Mb" if d < 1_000_000 else "same_chrom_far"
        else:
            loc = "cross_chrom (bridge-like)"
        coloc[(hc, loc)] += 1
    for hc, n in cat.most_common():
        print(f"  {n:5d}  {hc}")
    print("  -- no_cDNA_homology broken down by co-location (is it a missed local paralog or a bridge?) --")
    for (hc, loc), n in sorted(coloc.items()):
        if hc == "no_cDNA_homology":
            print(f"     {n:5d}  {loc}")
    same_chrom_nohom = sum(n for (hc, loc), n in coloc.items()
                           if hc == "no_cDNA_homology" and loc.startswith("same_chrom"))
    print(f"  => of {cat['no_cDNA_homology']} no-homology RNA edges, {same_chrom_nohom} are co-located "
          f"(likely real local paralogs the longest-cDNA truth missed); the rest are cross-chrom bridges.")

    # ---- RECALL decomposition funnel (vs LOOSE DNA) ----
    print("\n--- RECALL decomposition (vs LOOSE DNA): why DNA edges are 'missed' ---")
    funnel = collections.Counter()
    for (a, b) in (dna_loose - rna):
        ea, eb = a in expressed, b in expressed
        if not (ea and eb):
            funnel["1_unexpressed (>=1 copy no RNA)"] += 1
            continue
        crd = ev.get((a, b))
        if crd is None:
            funnel["2_no_cross_map (reads place uniquely)"] += 1
            continue
        tied = sum(1 for d, m in crd if d <= DELTA and m <= DE_MAX)
        if tied == 0:
            funnel["3_resolvable (cross-map, de-decidable)"] += 1
        elif tied < MIN_READS:
            funnel["4_sub_quorum (tied 1-2)"] += 1
        else:
            funnel["5_GENUINE_MISS (tied>=3, no edge)"] += 1
    for c in sorted(funnel):
        print(f"  {funnel[c]:6d}  {c}")
    tp = R["loose"]["tp"]
    crossmap_fn = funnel["3_resolvable (cross-map, de-decidable)"] + funnel["4_sub_quorum (tied 1-2)"] + funnel["5_GENUINE_MISS (tied>=3, no edge)"]
    tied_fn = funnel["4_sub_quorum (tied 1-2)"] + funnel["5_GENUINE_MISS (tied>=3, no edge)"]
    print(f"  recall | cross-mapping universe  = {tp}/({tp}+{crossmap_fn}) = {tp/max(tp+crossmap_fn,1):.3f}")
    print(f"  recall | TIED universe (drop resolvable) = {tp}/({tp}+{tied_fn}) = {tp/max(tp+tied_fn,1):.3f}")

    # ---- PANEL confusion vs DNA truth ----
    print("\n" + "=" * 80)
    print("PANEL (named pairs) confusion vs DNA truth")
    print("=" * 80)
    PANEL = [("RABL2A", "RABL2B"), ("APOBEC3D", "APOBEC3F"), ("RFPL1", "RFPL2"),
             ("RFPL2", "RFPL3"), ("CREB1", "METTL21A"), ("GCA", "KCNH7"),
             ("CASP8", "FLACC1"), ("ASDURF", "ASNSD1"), ("GPR39", "LYPD1")]
    for a, b in PANEL:
        key = (a, b) if a < b else (b, a)
        r = H.get(key)
        hs = (f"id={r['id']:.2f} covA={r['cov_a']:.2f} covB={r['cov_b']:.2f}") if r else "no aln"
        print(f"  {a:10}~{b:10} LOOSE-DNA={'Y' if key in dna_loose else 'n'}  "
              f"WHOLE-DNA={'Y' if key in dna_whole else 'n'}  RNA={'Y' if key in rna else 'n'}  ({hs})")

    # ---- per-edge dump for verification ----
    union = rna | dna_loose
    with open(EDGE_DUMP, "w") as o:
        o.write("gene_a\tgene_b\tin_rna\tin_dna_loose\tin_dna_whole\tn_tied\texpr_a\texpr_b\t"
                "cross_map\tchrom_a\tchrom_b\tsame_chrom\tdist_bp\tid\tcov_a\tcov_b\n")
        for (a, b) in sorted(union):
            r = H.get((a, b), {})
            ca = g2pos.get(a, ("?", -1, -1)); cb = g2pos.get(b, ("?", -1, -1))
            same = int(ca[0] == cb[0])
            dist = abs(ca[1] - cb[1]) if same else -1
            o.write(f"{a}\t{b}\t{int((a,b) in rna)}\t{int((a,b) in dna_loose)}\t{int((a,b) in dna_whole)}\t"
                    f"{n_tied.get((a,b),0)}\t{int(a in expressed)}\t{int(b in expressed)}\t"
                    f"{int((a,b) in ev)}\t{ca[0]}\t{cb[0]}\t{same}\t{dist}\t"
                    f"{r.get('id',0):.4f}\t{r.get('cov_a',0):.3f}\t{r.get('cov_b',0):.3f}\n")
    print(f"\n[+] wrote {EDGE_DUMP}  ({len(union):,} union edges)")

    out = dict(min_id=MIN_ID, min_cov=MIN_COV, recip_cov=RECIP_COV,
               rna_edges=len(rna), results=R, funnel=dict(funnel),
               rna_only_homology=dict(cat),
               recall_crossmap=tp/max(tp+crossmap_fn,1), recall_tied=tp/max(tp+tied_fn,1))
    with open(os.path.join(HERE, "family_def_dna_pr.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"[+] wrote {os.path.join(HERE, 'family_def_dna_pr.json')}")


if __name__ == "__main__":
    main()
