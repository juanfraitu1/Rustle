#!/usr/bin/env python3
"""family_def_newbam_validate.py — does the new multimapper BAM (-N 50 -p 0.1) strengthen
~R without breaking ~B's bridge-pruning?

The old GGO.bam used minimap2's default secondary cap, so ~R (read-confusability) was built
on a tiny fraction of the real cross-mapping. The new GGO_mm.bam (-N 50 -p 0.1 --secondary=yes)
surfaces ~65x more secondaries (de tags intact). Per paralog-rich chromosome, we re-run the
EXACT ~R scan (family_def_genomewide.scan logic, region-restricted) on OLD vs NEW and compare:
  (1) cross-mapping volume, candidate ~R pairs/families, quorum strength;
  (2) whether ~B (cached copy-model reciprocal coverage) still prunes the extra bridges.
Run: /home/juanfra/miniforge3/bin/python bench/family_def_newbam_validate.py NC_073244.2 [more chroms]
"""
import collections
import os
import subprocess
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import (de_of, best_gene, ref_span, pair_evidence, components,
                                    DELTA, DE_MAX, MIN_READS)

META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
SEQS = "/home/juanfra/winloci_scratch/vg_reinforce_copies.fa"
GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
OLD = "/home/juanfra/winloci_scratch/GGO.bam"
NEW = "/home/juanfra/winloci_scratch/GGO_mm.bam"
TMP = "/home/juanfra/winloci_scratch/newbamval"
COV_MIN = 0.30
PAD, MAX_READS, MIN_LEN = 5000, 3000, 300
SAM = "/home/juanfra/miniforge3/bin/samtools"


def merge(iv):
    iv.sort(); out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def build_model(bam, genome, chrom, start, end):
    """exon-union of primary reads at a locus (same as reinforce.build_model)."""
    blocks, n = [], 0
    try:
        it = bam.fetch(chrom, max(0, start - PAD), end + PAD)
    except (ValueError, KeyError):
        return ""
    for r in it:
        if r.is_unmapped or r.is_supplementary or r.is_secondary:
            continue
        if r.reference_end is None or r.reference_end < start or r.reference_start > end:
            continue
        b = r.get_blocks()
        if b:
            blocks.extend(b); n += 1
            if n >= MAX_READS:
                break
    if not blocks:
        return ""
    union = merge([list(x) for x in blocks])
    seq = "".join(genome.fetch(chrom, s, e).upper() for s, e in union)
    return seq if len(seq) >= MIN_LEN else ""


def load_denovo(chroms):
    by_chrom = collections.defaultdict(list)
    info = {}
    with open(META) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            c = p[1]
            if c not in chroms:
                continue
            vid, s, e = p[0], int(p[2]), int(p[3])
            by_chrom[c].append((s, e, vid))
            info[vid] = (c, s, e, int(p[6]))
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom, info


def scan_region(by_chrom, bam, region):
    """region-restricted version of family_def_genomewide.scan (same de/best_gene logic)."""
    mm = collections.defaultdict(dict)
    nsec = 0
    p = subprocess.Popen([SAM, "view", "-f", "0x100", bam, region],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        f = line.split("\t")
        if len(f) < 9:
            continue
        nsec += 1
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
    p = subprocess.Popen([SAM, "view", "-F", "0x900", bam, region],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        i = line.find("\t")
        if line[:i] not in mm:
            continue
        f = line.split("\t")
        de = de_of(f[11:])
        if de is None or de > DE_MAX:
            continue
        g = best_gene(by_chrom, f[2], int(f[3]) - 1, int(f[3]) - 1 + ref_span(f[5]))
        if g is None:
            continue
        d = mm[line[:i]]
        if g not in d or de < d[g]:
            d[g] = de
    p.wait()
    return mm, nsec


def load_models():
    models, name, buf = {}, None, []
    for line in open(SEQS):
        if line.startswith(">"):
            if name is not None:
                models[name] = "".join(buf)
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name is not None:
        models[name] = "".join(buf)
    return models


def recip_cov(a, b):
    os.makedirs(TMP, exist_ok=True)
    open(f"{TMP}/a.fa", "w").write(f">a\n{a}\n")
    open(f"{TMP}/b.fa", "w").write(f">b\n{b}\n")
    out = subprocess.run(f"minimap2 -c -x asm20 {TMP}/a.fa {TMP}/b.fa 2>/dev/null",
                         shell=True, capture_output=True, text=True).stdout
    best = 0.0
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 11:
            continue
        ql, qs, qe, tl, ts, te = int(f[1]), int(f[2]), int(f[3]), int(f[6]), int(f[7]), int(f[8])
        best = max(best, min((qe - qs) / max(ql, 1), (te - ts) / max(tl, 1)))
    return best


def run(label, bam, by_chrom, region):
    mm, nsec = scan_region(by_chrom, bam, region)
    multimappers = sum(1 for g in mm.values() if len(g) >= 2)
    ev = pair_evidence(mm)
    edges, fams = components(ev, DELTA, DE_MAX, MIN_READS)
    # quorum strength: conflicting-read count per candidate edge (n is in the edge tuple)
    counts = sorted((n for (ga, gb, n) in edges), reverse=True)
    print(f"  [{label}] secondaries_scanned={nsec}  multimapper_reads={multimappers}  "
          f"candidate_pairs(edges)={len(edges)}  candidate_families={len(fams)}")
    if counts:
        print(f"           quorum per edge: max={counts[0]} median={counts[len(counts)//2]} "
              f"edges>=10reads={sum(1 for c in counts if c>=10)}")
    return dict(mm=mm, ev=ev, edges=edges, fams=fams, multimappers=multimappers)


def prune_b(edges, models):
    """split candidate edges into ~B-pass (real copies) vs ~B-fail (bridges)."""
    kept = bridges = 0
    for (a, b, n) in edges:
        if not models.get(a) or not models.get(b):
            continue
        if recip_cov(models[a], models[b]) >= COV_MIN:
            kept += 1
        else:
            bridges += 1
    return kept, bridges


def main():
    chroms = sys.argv[1:] or ["NC_073244.2"]
    by_chrom, info = load_denovo(set(chroms))
    models = load_models()
    bam_old = pysam.AlignmentFile(OLD, "rb")
    genome = pysam.FastaFile(GENOME)
    print(f"=== ~R re-validation on new multimapper BAM ({','.join(chroms)}) ===")
    print(f"de-novo loci on these chroms: {sum(len(v) for v in by_chrom.values())}")
    for region in chroms:
        print(f"\n--- {region} ---")
        old = run("OLD", OLD, by_chrom, region)
        new = run("NEW", NEW, by_chrom, region)
        # resolve copy models for EVERY candidate locus (cache, else build from old BAM)
        loci = {g for (a, b, n) in list(old["edges"]) + list(new["edges"]) for g in (a, b)}
        built = 0
        for v in loci:
            if not models.get(v) and v in info:
                c, s, e, _ = info[v]
                models[v] = build_model(bam_old, genome, c, s, e)
                if models[v]:
                    built += 1
        print(f"  resolved models for {len(loci)} candidate loci ({built} built on-demand)")
        ok, ob = prune_b(old["edges"], models)
        nk, nb = prune_b(new["edges"], models)
        print(f"  [~B] OLD: {len(old['edges'])} edges -> {ok} real-copy | {ob} bridge")
        print(f"  [~B] NEW: {len(new['edges'])} edges -> {nk} real-copy | {nb} bridge")
        print(f"  NET: new BAM = {new['multimappers']-old['multimappers']:+d} multimapper reads, "
              f"{len(new['edges'])-len(old['edges']):+d} candidate edges, "
              f"{new['multimappers']/max(old['multimappers'],1):.1f}x cross-mapping; "
              f"~B keeps {nk} real copies (vs {ok}), prunes {nb} bridges (vs {ob}).")
    bam_old.close(); genome.close()


if __name__ == "__main__":
    main()
