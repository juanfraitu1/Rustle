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
import shutil
import subprocess
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import (de_of, best_gene, ref_span, pair_evidence, components,
                                    DELTA, DE_MAX, MIN_READS)

META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
SEQS = "/home/juanfra/winloci_scratch/vg_reinforce_copies.fa"
GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
# ledger §6am: OLD was "GGO.bam", which has been a SYMLINK to GGO_mm.bam since 2026-06-23 -> both arms
# were ONE inode and the OLD-vs-NEW comparison was against itself, so it could not fail. The
# default-secondary-cap alignment the recorded table (bench/FAMILY_DEF.md) was measured on survives
# under its original bytes as GGO_orig_defaultcap.bam; that file is the OLD arm. Both arms overridable.
OLD = os.environ.get("NEWBAM_OLD", "/home/juanfra/winloci_scratch/GGO_orig_defaultcap.bam")
NEW = os.environ.get("NEWBAM_NEW", "/home/juanfra/winloci_scratch/GGO_mm.bam")
TMP = "/home/juanfra/winloci_scratch/newbamval"
COV_MIN = 0.30
PAD, MAX_READS, MIN_LEN = 5000, 3000, 300
SAM = "/home/juanfra/miniforge3/bin/samtools"
MM2 = os.environ.get("MINIMAP2", "minimap2")


def _abort(msg):
    """Every guard exits nonzero with the reason (ledger §6am)."""
    sys.exit(f"ABORT family_def_newbam_validate: {msg}")


def check_inputs():
    """Pre-flight, before any scan. §6am: this instrument reported a comparison it never made."""
    for label, path in (("meta table", META), ("copy-model cache", SEQS), ("genome", GENOME),
                        ("OLD bam", OLD), ("NEW bam", NEW), ("samtools", SAM)):
        # §6am: a missing input let samtools/pysam yield nothing and the empty arms diffed as "no change".
        if not os.path.exists(path):
            _abort(f"{label} missing: {path}")
    for label, bam in (("OLD", OLD), ("NEW", NEW)):
        # §6am: without an index the region fetch returns 0 alignments, which scored as agreement.
        if not any(os.path.exists(bam + ext) for ext in (".bai", ".csi")):
            _abort(f"{label} bam has no .bai/.csi index: {bam} (region-restricted scan needs one)")
    # §6am, THE fix: if the arms are the same file the instrument is comparing a BAM to itself.
    so, sn = os.stat(OLD), os.stat(NEW)
    if (so.st_dev, so.st_ino) == (sn.st_dev, sn.st_ino):
        _abort(f"OLD and NEW are the SAME FILE (dev {so.st_dev}, inode {so.st_ino}): "
               f"{OLD} -> {os.path.realpath(OLD)} == {NEW} -> {os.path.realpath(NEW)}; "
               f"the comparison is against itself and cannot fail")
    # §6am: minimap2 absent => shell rc 127, coverage 0.0 for every pair, i.e. ~B "prunes" everything.
    if shutil.which(MM2) is None:
        _abort(f"minimap2 not found on PATH as {MM2!r} (~B reciprocal coverage needs it)")
    print(f"arms: OLD={os.path.realpath(OLD)}\n      NEW={os.path.realpath(NEW)}")


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
    rc = p.wait()
    # §6am: an unchecked samtools failure (bad region, missing index) returned 0 rows and the arm scored
    # as "no cross-mapping" instead of failing.
    if rc != 0:
        _abort(f"samtools view -f 0x100 {bam} {region} failed (exit {rc})")
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
    rc = p.wait()
    # §6am: same defect on the primary/supplementary pass - a silent failure dropped the quorum reads.
    if rc != 0:
        _abort(f"samtools view -F 0x900 {bam} {region} failed (exit {rc})")
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
    r = subprocess.run(f"{MM2} -c -x asm20 {TMP}/a.fa {TMP}/b.fa",
                       shell=True, capture_output=True, text=True)
    # §6am: stderr went to /dev/null and the exit code was ignored, so a failed aligner returned an empty
    # PAF, coverage 0.0, and every candidate edge was reported as a ~B-pruned "bridge".
    if r.returncode != 0:
        _abort(f"minimap2 failed (exit {r.returncode}): {r.stderr.strip()[:400]}")
    out = r.stdout
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
    return dict(mm=mm, ev=ev, edges=edges, fams=fams, multimappers=multimappers, nsec=nsec)


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
    check_inputs()
    by_chrom, info = load_denovo(set(chroms))
    # §6am: with no de-novo locus on these chroms best_gene() is None for every read, both arms score 0
    # edges, and the NET line prints "+0" as if the two BAMs agreed. An empty evidence set is not a result.
    if sum(len(v) for v in by_chrom.values()) == 0:
        _abort(f"no de-novo loci in {META} for {','.join(chroms)} (nothing to scan; check the chrom names)")
    models = load_models()
    # §6am: an empty copy-model cache makes ~B call every candidate edge a bridge, which reads as a pass.
    if not models:
        _abort(f"copy-model cache {SEQS} yielded 0 sequences")
    bam_old = pysam.AlignmentFile(OLD, "rb")
    genome = pysam.FastaFile(GENOME)
    print(f"=== ~R re-validation on new multimapper BAM ({','.join(chroms)}) ===")
    print(f"de-novo loci on these chroms: {sum(len(v) for v in by_chrom.values())}")
    for region in chroms:
        print(f"\n--- {region} ---")
        old = run("OLD", OLD, by_chrom, region)
        new = run("NEW", NEW, by_chrom, region)
        # §6am: an arm that scanned no secondary alignment contributes "0" to every column, and the
        # OLD-vs-NEW deltas then read as agreement. Refuse to report a comparison one arm never made.
        if old["nsec"] == 0 or new["nsec"] == 0:
            _abort(f"{region}: secondaries scanned OLD={old['nsec']} NEW={new['nsec']} - an arm returned "
                   f"no alignments on this region, so the comparison is empty")
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
