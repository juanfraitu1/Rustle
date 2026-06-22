#!/usr/bin/env python3
"""
family_def_vg_pipeline.py — genome-wide multi-copy gene family definition with VG
structural validation.

  DISCOVER : read-conflict (de-tie) graph over annotated-gene vertices -> candidate
             families (which loci are confusable).
  VALIDATE : for each gene in a candidate family, build its ALL-ISOFORM exon-union copy
             model from reads (isoforms aggregate, not pollute); align copies; an edge
             is KEPT only if the two copies share a real backbone (parallel paths /
             asymmetric real copy), REJECTED if no backbone (repeat-bridge).
  RESULT   : several families found, several edges rejected as bridges, only the real
             families kept. Audited against the DNA homology truth.

OOM-safe: copy models are built ONLY for genes in a candidate family (~1.7k, not 34k);
reads capped per locus; sequences streamed to disk; one small minimap2 pass.

Run (background): /home/juanfra/miniforge3/bin/python bench/family_def_vg_pipeline.py
"""
import collections
import json
import os
import subprocess
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_genomewide import (load_vertices, scan as conflict_scan, pair_evidence,
                                    components, DELTA, DE_MAX, MIN_READS)
from family_def_read_filters import dna_homology

BAM = "/home/juanfra/winloci_scratch/GGO.bam"
GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
GENES_BED = "/home/juanfra/winloci_scratch/unmapped_poc/genes.bed"
SEQS = "/home/juanfra/winloci_scratch/vg_pipeline_copies.fa"
PAF = "/home/juanfra/winloci_scratch/vg_pipeline_copies.paf"
PAD = 20000
MAX_READS = 3000
ID_MIN = 0.80
MAXCOV_MIN = 0.20     # below this in BOTH directions = no backbone = repeat-bridge


def gene_coords():
    d = {}
    with open(GENES_BED) as f:
        for line in f:
            c, s, e, name = line.rstrip("\n").split("\t")
            d.setdefault(name, (c, int(s), int(e)))
    return d


def merge(iv):
    iv.sort()
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def build_copy_model(bam, genome, chrom, start, end):
    blocks = []
    n = 0
    s0 = max(0, start - PAD)
    try:
        it = bam.fetch(chrom, s0, end + PAD)
    except (ValueError, KeyError):
        return ""
    for r in it:
        if r.is_unmapped or r.is_supplementary:
            continue
        if r.reference_end is None or r.reference_end < start or r.reference_start > end:
            continue
        b = r.get_blocks()
        if b:
            blocks.extend(b)
            n += 1
            if n >= MAX_READS:
                break
    if not blocks:
        return ""
    union = merge([list(x) for x in blocks])
    return "".join(genome.fetch(chrom, s, e).upper() for s, e in union)


def main():
    print("[1/4 discover] read-conflict graph over annotated-gene vertices ...", flush=True)
    by_chrom, v2chrom = load_vertices("genes")
    mm, n_sec, n_prim = conflict_scan(by_chrom)
    ev = pair_evidence(mm)
    edges, fams = components(ev, DELTA, DE_MAX, MIN_READS)
    edge_set = {(a, b) for a, b, n in edges}
    fam_genes = sorted({g for c in fams for g in c})
    print(f"  candidate: {len(fams)} families, {len(edge_set)} edges, {len(fam_genes)} genes", flush=True)

    print(f"[2/4 build] all-isoform exon-union copy models for {len(fam_genes)} family genes ...", flush=True)
    coords = gene_coords()
    bam = pysam.AlignmentFile(BAM, "rb")
    genome = pysam.FastaFile(GENOME)
    built = 0
    with open(SEQS, "w") as out:
        for i, g in enumerate(fam_genes):
            if g not in coords:
                continue
            c, s, e = coords[g]
            seq = build_copy_model(bam, genome, c, s, e)
            if len(seq) >= 200:
                out.write(f">{g}\n{seq}\n")
                built += 1
            if (i + 1) % 250 == 0:
                print(f"    ... {i+1}/{len(fam_genes)} ({built} models)", flush=True)
    bam.close(); genome.close()
    print(f"  built {built} copy models", flush=True)

    print("[3/4 validate] minimap2 backbone alignment of copy models ...", flush=True)
    with open(PAF, "w") as pf:
        subprocess.run(["minimap2", "-cx", "asm20", "-N", "20", "-p", "0.05", "-t", "6",
                        SEQS, SEQS], stdout=pf, stderr=subprocess.DEVNULL)
    H = {}
    with open(PAF) as f:
        for line in f:
            x = line.split("\t")
            if len(x) < 11:
                continue
            q, ql, qs, qe, t, m, bl = x[0], int(x[1]), int(x[2]), int(x[3]), x[5], int(x[9]), int(x[10])
            if q == t:
                continue
            ident = m / bl if bl else 0
            cov = (qe - qs) / ql if ql else 0
            key = (q, t) if q < t else (t, q)
            r = H.setdefault(key, {"id": 0.0, "ca": 0.0, "cb": 0.0})
            r["id"] = max(r["id"], ident)
            if q == key[0]:
                r["ca"] = max(r["ca"], cov)
            else:
                r["cb"] = max(r["cb"], cov)

    def has_backbone(a, b):
        r = H.get((a, b) if a < b else (b, a))
        if not r:
            return False
        return r["id"] >= ID_MIN and max(r["ca"], r["cb"]) >= MAXCOV_MIN

    kept = {(a, b) for (a, b) in edge_set if has_backbone(a, b)}
    rejected = edge_set - kept
    _, fams_kept = components({(a, b): ev[(a, b)] for (a, b) in kept}, DELTA, DE_MAX, MIN_READS)

    print("[4/4 result]", flush=True)
    print(f"  candidate families : {len(fams)}   edges : {len(edge_set)}")
    print(f"  REJECTED edges (no backbone = repeat-bridge): {len(rejected)} "
          f"({100*len(rejected)/max(len(edge_set),1):.0f}%)")
    print(f"  VALIDATED families : {len(fams_kept)}   kept edges : {len(kept)}")

    # audit the rejection against the DNA truth
    H_dna, dna = dna_homology()
    def cls(e):
        r = H_dna.get(e)
        if r is None or r["id"] == 0:
            return "no_homology"
        if r["id"] >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30:
            return "DNA_paralog"
        return "sub_bar" if r["id"] >= 0.80 else "weak"
    rc = collections.Counter(cls(e) for e in rejected)
    kc = collections.Counter(cls(e) for e in kept)
    print("\n  rejected edges by DNA class (want: mostly no_homology = correct rejection):")
    for k, v in rc.most_common():
        print(f"     {v:5d}  {k}")
    print("  kept edges by DNA class (want: enriched for DNA_paralog):")
    for k, v in kc.most_common():
        print(f"     {v:5d}  {k}")
    p_before = sum(1 for e in edge_set if e in dna) / max(len(edge_set), 1)
    p_after = sum(1 for e in kept if e in dna) / max(len(kept), 1)
    print(f"\n  precision vs DNA-paralog:  before={p_before:.3f}  ->  after VG validation={p_after:.3f}")

    # show the marquee bridges' fate
    print("\n  marquee bridges (should be REJECTED):")
    for a, b in [("OCLN", "SEPTIN7"), ("BCAS4", "CCDC30")]:
        e = (a, b) if a < b else (b, a)
        if e in edge_set:
            r = H.get(e, {})
            print(f"     {a}~{b}: {'REJECTED' if e in rejected else 'kept'}  "
                  f"(id={r.get('id',0):.2f} maxcov={max(r.get('ca',0),r.get('cb',0)):.2f})")

    json.dump(dict(candidate_families=len(fams), candidate_edges=len(edge_set),
                   rejected=len(rejected), validated_families=len(fams_kept), kept_edges=len(kept),
                   precision_before=p_before, precision_after=p_after,
                   rejected_classes=dict(rc), kept_classes=dict(kc)),
              open(os.path.join(HERE, "family_def_vg_pipeline.json"), "w"), indent=2)
    print(f"\n[+] wrote {os.path.join(HERE,'family_def_vg_pipeline.json')}")


if __name__ == "__main__":
    main()
