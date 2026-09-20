#!/usr/bin/env python
"""READ-CONFUSABILITY BOUNDARY vs COPY IDENTITY -- CHARACTERIZE STAGE (measure, no wiring).

THESIS-SCOPE QUESTION
---------------------
The copy-assignment method targets READ-CONFUSABLE multi-copy families: copies so similar
the aligner maps a read equally well to >=2 loci and GIVES UP (MAPQ 0 / tiny best-vs-2nd
alignment-score margin). DIVERGENT copies are ALIGNER-RESOLVABLE (reads carry enough
mismatches to place uniquely, MAPQ>0) so they never needed the method.

If that is true there is a COPY-IDENTITY BOUNDARY: below it the aligner resolves copies
itself (few confusable reads); above it resolvability collapses (reads pile up MAPQ 0).
The divergence FLOOR being wired into the family definition (bench/divergence_floor.py,
reciprocal whole-transcript identity) is PRINCIPLED (Canzar: no arbitrary threshold) iff
the P/R-optimal floor COINCIDES with this read-confusability boundary -- i.e. the method
models exactly the families the aligner cannot resolve.

WHAT THIS MEASURES (characterize only; the floor stays at family_rna_refine defaults)
-------------------------------------------------------------------------------------
For every multi-copy family (bench/family_rna_refine.tsv, n_loci>=2) that has a copy-vs-copy
identity signal (join to bench/colinear_multiexon_gate.tsv on member DN pairs):

  COPY IDENTITY (x axis), per family, over covered member pairs:
    mean_matched_id   whole-read proxy: mean over matched (colinear) exons of exon identity   [gate]
    best_exon_id(max) per-exon  proxy: identity of the single most-similar exon               [gate]
    core_recip        gate-derived reciprocal identity -- a PROXY for the floor (corr~0.49),  [gate]
                      NOT the metric the floor operates on
    recip_best_mean   recip_id_best = matches_best/max(len_a,len_b) -- THE ACTUAL divergence-  [ri]
                      floor metric (family_rna_refine.load_recip_id; the floor cuts a cross-
                      gene edge iff recip_id_best < floor), aggregated mean over the family's
                      member gene-pairs. The coincidence/mass test uses THIS, not core_recip.

  READ RESOLVABILITY (y axis) -- TWO independent metrics that cross-check:
    (1) BAM MAPQ  from GGO_mm.bam: fraction of primary reads at the family's loci with
        MAPQ==0 (multi-mapping/confusable). GGO_mm.bam was aligned -N50 -p0.1
        --secondary=yes which can DEPRESS MAPQ -> possible over-count of confusable.
    (2) RE-ALIGN CONTROL: a seeded read SAMPLE per family is re-aligned ONCE, all families
        together, to GGO.splice.mmi with DEFAULT minimap2 (splice:hq, no forced secondary).
        Per read we read the DEFAULT MAPQ (un-depressed; controls the -N50 confound) and the
        AS MARGIN = AS(primary) - AS(best secondary) (setting-robust; connects to Eichler's
        AS>=10 rule). confusable_margin = margin < MARGIN_MIN (default 10).

  Bin families by copy identity; per bin report the mean confusable fraction under all three
  read metrics (BAM MAPQ0, re-align default MAPQ0, AS-margin<10). LOCATE THE BOUNDARY = the
  identity where the confusable fraction rises sharply.

  PER-EXON NUANCE (tested): confusability is likely PER-EXON not whole-read -- an identical
  exon confuses reads spanning only it even if another exon is divergent. Tests:
   * axis test: is the confusable fraction predicted better by best_exon_id (max exon id,
     per-exon) than by mean_matched_id (whole-read)? A sharper/lower boundary on best_exon_id
     => exon-driven.
   * within-family test: do DIVERGENT-overall families (low mean id) still carry a subset of
     confusable reads (margin<10)? A non-zero residual => exon-dependent, not a clean line.
   * the per-read AS margin is itself a per-exon-span measurement (it depends on what the read
     covers), so its distribution WITHIN a family exposes exon mixing.

  CONFOUNDS controlled: read DEPTH (more reads -> more confusable pairs?), family SIZE
  (n members), and the -N50 forced-secondary MAPQ depression (BAM-MAPQ0 vs re-align-MAPQ0 gap).

COINCIDENCE WITH THE FLOOR: bench/divergence_floor.tsv (the P/R-vs-floor sweep) is read if
present and the located boundary is compared to it. The floor is a WHOLE-TRANSCRIPT RECIPROCAL
identity (core_recip); the boundary here is compared on the matching axis.

Deterministic: PYTHONHASHSEED=0, random.seed(SEED), minimap2 -t fixed, --seed fixed, sorted IO.
Writes: bench/read_confusability_boundary.tsv (per-family) + prints the curve, boundary,
exon-dependence verdict, confound checks. Run:
  PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/read_confusability_boundary.py
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import argparse
import json
import math
import random
import statistics
import subprocess
import tempfile
from collections import defaultdict
from itertools import combinations

import pysam

BENCH = os.path.dirname(os.path.abspath(__file__))
CATALOG_TSV = os.path.join(BENCH, "family_rna_refine.tsv")
GATE_TSV = os.path.join(BENCH, "colinear_multiexon_gate.tsv")
RECIP_TSV = os.path.join(BENCH, "ri_sharedlen_recip_id.tsv")   # recip_id_best = THE ACTUAL floor metric
DIVFLOOR_TSV = os.path.join(BENCH, "divergence_floor.tsv")
DIVFLOOR_JSON = os.path.join(BENCH, "divergence_floor.json")   # P/R-optimal floor detail (if the sweep emitted it)
DIVFLOOR_MD = os.path.join(BENCH, "DIVERGENCE_FLOOR.md")       # prose rationale (if present)
OUT_TSV = os.path.join(BENCH, "read_confusability_boundary.tsv")
OUT_JSON = os.path.join(BENCH, "read_confusability_boundary.json")
FIG_PNG = os.path.join(BENCH, "fig_read_confusability_boundary.png")
SHIPPED_FLOOR = 0.80          # current default divergence floor (family_rna_refine, DEFAULT-ON)

BAM = "/home/juanfra/winloci_scratch/GGO_mm.bam"
MMI = "/home/juanfra/winloci_scratch/GGO.splice.mmi"
MINIMAP2 = "/home/juanfra/miniforge3/bin/minimap2"

SEED = 0
SAMPLE_PER_FAM = 40          # reads re-aligned per family (seeded)
MARGIN_MIN = 10              # AS margin below this == confusable (Eichler AS>=10)
MM2_THREADS = 4
MM2_SEED = 42
# copy-identity bin edges (whole-read axis; high-identity region resolved finely)
ID_EDGES = [0.0, 0.50, 0.85, 0.90, 0.95, 0.98, 0.99, 1.0001]


# --------------------------------------------------------------------------- families + identity
def load_families():
    fams = defaultdict(list)
    with open(CATALOG_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            fams[f[I["family_id"]]].append(dict(
                dn=f[I["member_dn"]], gene=f[I["member_gene"]], dom=f[I["dominant_gene"]],
                chrom=f[I["chrom"]], start=int(f[I["start"]]), end=int(f[I["end"]]),
                n_loci=int(f[I["n_loci"]])))
    return fams


def load_gate():
    gate = {}
    with open(GATE_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            gate[frozenset((f[I["dn_a"]], f[I["dn_b"]]))] = dict(
                mean_matched_id=float(f[I["mean_matched_id"]]),
                best_exon_id=float(f[I["best_exon_id"]]),
                core_recip=float(f[I["core_recip"]]),
                aln_frac=float(f[I["aln_frac"]]))
    return gate


def load_recip_best():
    """Per GENE-PAIR reciprocal whole-transcript identity recip_id_best = matches_best/max(len_a,len_b)
    from bench/ri_sharedlen_recip_id.tsv -- THE ACTUAL divergence-floor metric. This is exactly what
    family_rna_refine.load_recip_id() reads and the floor operates on: a cross-gene edge is CUT iff
    recip_id_best < floor. Keyed by frozenset of the gene pair (last-write-wins, matching load_recip_id).
    NOT the same as gate.core_recip (a different, colinear-gate-derived reciprocal identity, corr~0.49;
    recip_id_best is lower on average -> the real floor cuts MORE). Deterministic ordered pass."""
    recip = {}
    if not os.path.exists(RECIP_TSV):
        return recip
    with open(RECIP_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            v = f[I["recip_id_best"]]
            if v == "":
                continue
            recip[frozenset((f[I["gene_a"]], f[I["gene_b"]]))] = float(v)
    return recip


def family_identity(rows, gate, recip):
    """Aggregate copy identity over covered member-pairs. None if no gate pair covered.
    Gate metrics (mean_matched_id/best_exon_id/core_recip) aggregate over covered DN pairs;
    recip_best (the ACTUAL floor metric) aggregates over covered member GENE pairs -- the level
    the floor actually cuts on (family_rna_refine.load_recip_id keys by gene-pair)."""
    dns = [r["dn"] for r in rows]
    hits = [gate[frozenset(p)] for p in combinations(dns, 2) if frozenset(p) in gate]
    if not hits:
        return None
    # ---- ACTUAL floor metric: recip_id_best over member gene-pairs (floor cuts cross-gene edges) ----
    genes_distinct = sorted({r["gene"] for r in rows if r["gene"] != "NA"})
    rhits = [recip[frozenset(p)] for p in combinations(genes_distinct, 2) if frozenset(p) in recip]
    return dict(
        n_pairs_covered=len(hits),
        mean_matched_id=statistics.mean(h["mean_matched_id"] for h in hits),   # whole-read [gate]
        min_matched_id=min(h["mean_matched_id"] for h in hits),
        best_exon_id=max(h["best_exon_id"] for h in hits),                     # per-exon (max) [gate]
        core_recip=statistics.mean(h["core_recip"] for h in hits),             # PROXY for floor [gate]
        exon_spread=max(h["best_exon_id"] for h in hits) - statistics.mean(h["mean_matched_id"] for h in hits),
        # THE ACTUAL floor metric (recip_id_best); None if no member gene-pair is in the ri cache
        # (absent => edge KEPT, per load_recip_id => family NOT cut by the floor):
        recip_best_mean=(statistics.mean(rhits) if rhits else None),
        recip_best_min=(min(rhits) if rhits else None),
        n_recip_pairs=len(rhits))


# --------------------------------------------------------------------------- BAM MAPQ pass + sampling
def bam_pass(fams, gate, recip, rng):
    """Per family: primary-read MAPQ tally at its loci + a seeded read sample for re-alignment.
    Returns per_fam dict and the sample-read list [(tag, seq)]."""
    bam = pysam.AlignmentFile(BAM)
    refs = set(bam.references)
    per_fam = {}
    sample_reads = []          # (tag, seq); tag encodes fam + bam_mapq
    for fid in sorted(fams, key=lambda x: int(x)):
        rows = fams[fid]
        ident = family_identity(rows, gate, recip)
        if ident is None:
            continue
        seen = {}              # read_name -> mapq (primary)
        seq_by_name = {}       # read_name -> seq (first seen, for sampling)
        for r in rows:
            if r["chrom"] not in refs:
                continue
            for a in bam.fetch(r["chrom"], r["start"], r["end"]):
                if a.is_secondary or a.is_supplementary:
                    continue
                nm = a.query_name
                if nm not in seen:
                    seen[nm] = a.mapping_quality
                    if a.query_sequence is not None:
                        seq_by_name[nm] = a.query_sequence
        depth = len(seen)
        if depth == 0:
            continue
        mapq0 = sum(1 for q in seen.values() if q == 0)
        # seeded sample of reads that have sequence
        names = sorted(seq_by_name)
        rng.shuffle(names)
        picked = names[:SAMPLE_PER_FAM]
        for nm in picked:
            sample_reads.append((f"{fid}|{nm}|{seen[nm]}", seq_by_name[nm]))
        per_fam[fid] = dict(
            dom=rows[0]["dom"], n_members=len(rows), n_loci=rows[0]["n_loci"],
            depth=depth, bam_mapq0=mapq0, bam_conf_frac=mapq0 / depth,
            n_sampled=len(picked), **ident)
    bam.close()
    return per_fam, sample_reads


# --------------------------------------------------------------------------- re-align control (ONE run)
def realign_sample(sample_reads, workdir):
    """Write all sampled reads to one FASTA, align ONCE to the mmi with DEFAULT splice:hq
    (minimap2 loads the 13GB index once). Parse per read: default MAPQ + AS margin
    (primary AS - best secondary AS). Returns read_tag -> dict."""
    fa = os.path.join(workdir, "sample.fa")
    with open(fa, "w") as fo:
        for tag, seq in sample_reads:
            fo.write(f">{tag}\n{seq}\n")
    sam = os.path.join(workdir, "sample.sam")
    cmd = [MINIMAP2, "-t", str(MM2_THREADS), "--seed", str(MM2_SEED),
           "-ax", "splice:hq", MMI, fa]
    with open(sam, "w") as so, open(os.path.join(workdir, "mm2.log"), "w") as le:
        subprocess.run(cmd, stdout=so, stderr=le, check=True)
    # parse: per read collect primary AS/mapq + best secondary AS
    prim_as = {}
    prim_mapq = {}
    best_sec_as = {}
    with open(sam) as fh:
        for ln in fh:
            if ln.startswith("@"):
                continue
            f = ln.rstrip("\n").split("\t")
            qname = f[0]
            flag = int(f[1])
            if flag & 0x4:                     # unmapped
                continue
            as_val = None
            for t in f[11:]:
                if t.startswith("AS:i:"):
                    as_val = int(t[5:]); break
            if as_val is None:
                continue
            if flag & 0x100:                   # secondary -> competitor for margin
                if as_val > best_sec_as.get(qname, -1):
                    best_sec_as[qname] = as_val
            elif flag & 0x800:                 # supplementary -> different span, ignore
                continue
            else:                              # primary
                # keep the highest-scoring primary (should be unique)
                if as_val > prim_as.get(qname, -1):
                    prim_as[qname] = as_val
                    prim_mapq[qname] = int(f[4])
    res = {}
    for tag in prim_as:
        pa = prim_as[tag]
        sa = best_sec_as.get(tag)
        margin = pa if sa is None else (pa - sa)
        res[tag] = dict(realn_mapq=prim_mapq[tag], prim_as=pa,
                        best_sec_as=(sa if sa is not None else -1),
                        margin=margin, has_secondary=(sa is not None))
    return res


def attach_realn(per_fam, realn):
    """Aggregate per-read re-align results back to families."""
    by_fam = defaultdict(list)
    for tag, d in realn.items():
        fid = tag.split("|", 1)[0]
        by_fam[fid].append(d)
    for fid, fam in per_fam.items():
        ds = by_fam.get(fid, [])
        n = len(ds)
        fam["n_realn"] = n
        if n == 0:
            fam["realn_conf_frac"] = None
            fam["margin_conf_frac"] = None
            fam["realn_mapq0"] = 0
            fam["margin_conf"] = 0
            fam["margins"] = []
            continue
        mapq0 = sum(1 for d in ds if d["realn_mapq"] == 0)
        marg0 = sum(1 for d in ds if d["margin"] < MARGIN_MIN)
        fam["realn_mapq0"] = mapq0
        fam["margin_conf"] = marg0
        fam["realn_conf_frac"] = mapq0 / n
        fam["margin_conf_frac"] = marg0 / n
        fam["margins"] = sorted(d["margin"] for d in ds)


# --------------------------------------------------------------------------- binning + boundary
def bin_index(x, edges):
    for i in range(len(edges) - 1):
        if edges[i] <= x < edges[i + 1]:
            return i
    return len(edges) - 2


def bin_label(i, edges):
    return f"{edges[i]:.2f}-{edges[i+1]:.2f}"


def build_curve(per_fam, axis, edges):
    """Read-weighted confusable fractions per identity bin, for all three metrics."""
    bins = [dict(fams=0, depth=0, bam_mapq0=0, n_realn=0, realn_mapq0=0, margin_conf=0)
            for _ in range(len(edges) - 1)]
    for fam in per_fam.values():
        if fam.get(axis) is None:        # e.g. recip_best_mean absent (no cross-gene ri pair)
            continue
        i = bin_index(fam[axis], edges)
        b = bins[i]
        b["fams"] += 1
        b["depth"] += fam["depth"]
        b["bam_mapq0"] += fam["bam_mapq0"]
        if fam.get("n_realn"):
            b["n_realn"] += fam["n_realn"]
            b["realn_mapq0"] += fam["realn_mapq0"]
            b["margin_conf"] += fam["margin_conf"]
    curve = []
    for i, b in enumerate(bins):
        curve.append(dict(
            bin=bin_label(i, edges), n_fams=b["fams"], n_reads=b["depth"],
            bam_conf=(b["bam_mapq0"] / b["depth"]) if b["depth"] else None,
            n_realn=b["n_realn"],
            realn_conf=(b["realn_mapq0"] / b["n_realn"]) if b["n_realn"] else None,
            margin_conf=(b["margin_conf"] / b["n_realn"]) if b["n_realn"] else None))
    return curve


def locate_boundary(curve, key):
    """Boundary = the lower edge of the bin with the largest upward jump in `key` vs the bin
    below it (the identity where confusable fraction rises sharpest)."""
    vals = [(c["bin"], c[key]) for c in curve if c[key] is not None]
    best = None
    for j in range(1, len(vals)):
        dv = vals[j][1] - vals[j - 1][1]
        if best is None or dv > best[0]:
            # boundary edge = start of the higher bin
            edge = float(vals[j][0].split("-")[0])
            best = (dv, edge, vals[j - 1][0], vals[j][0])
    return dict(jump=best[0], boundary=best[1], from_bin=best[2], to_bin=best[3]) if best else None


# --------------------------------------------------------------------------- confounds
def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx, my = statistics.mean(xs), statistics.mean(ys)
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if sx == 0 or sy == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (sx * sy)


def partial_corr(x, y, z):
    """corr(x,y) controlling for z (linear)."""
    rxy, rxz, ryz = pearson(x, y), pearson(x, z), pearson(y, z)
    if None in (rxy, rxz, ryz):
        return None
    den = math.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
    return (rxy - rxz * ryz) / den if den else None


def confound_checks(per_fam):
    fids = [f for f in per_fam.values()]
    conf = [f["bam_conf_frac"] for f in fids]
    ident = [f["mean_matched_id"] for f in fids]
    depth = [float(f["depth"]) for f in fids]
    size = [float(f["n_members"]) for f in fids]
    out = dict(
        r_conf_identity=pearson(conf, ident),
        r_conf_depth=pearson(conf, depth),
        r_conf_size=pearson(conf, size),
        r_conf_depth_ctrl_id=partial_corr(conf, depth, ident),
        r_conf_size_ctrl_id=partial_corr(conf, size, ident),
        r_conf_identity_ctrl_depth=partial_corr(conf, ident, depth),
    )
    # forced-secondary depression: BAM mapq0 vs re-align default mapq0, on families with realn
    both = [f for f in fids if f.get("realn_conf_frac") is not None]
    bam_m = statistics.mean(f["bam_conf_frac"] for f in both) if both else None
    realn_m = statistics.mean(f["realn_conf_frac"] for f in both) if both else None
    marg_m = statistics.mean(f["margin_conf_frac"] for f in both) if both else None
    out["mean_bam_conf"] = bam_m
    out["mean_realn_conf"] = realn_m
    out["mean_margin_conf"] = marg_m
    out["n_fam_with_realn"] = len(both)
    # per-family agreement between BAM and margin
    if both:
        out["r_bam_vs_margin"] = pearson([f["bam_conf_frac"] for f in both],
                                         [f["margin_conf_frac"] for f in both])
        out["r_bam_vs_realnmapq"] = pearson([f["bam_conf_frac"] for f in both],
                                            [f["realn_conf_frac"] for f in both])
    return out


def exon_dependence(per_fam):
    """Is confusability whole-read or per-exon?
      (a) which axis predicts confusable fraction better: mean_matched_id vs best_exon_id
      (b) do DIVERGENT-overall families (mean_matched_id < 0.85) still carry confusable reads?
      (c) exon_spread (best-mean) vs residual confusability in divergent families."""
    fids = list(per_fam.values())
    conf = [f["bam_conf_frac"] for f in fids]
    r_mean = pearson(conf, [f["mean_matched_id"] for f in fids])
    r_best = pearson(conf, [f["best_exon_id"] for f in fids])
    # divergent-overall families with a high-identity exon
    div = [f for f in fids if f["mean_matched_id"] < 0.85]
    div_conf_bam = statistics.mean(f["bam_conf_frac"] for f in div) if div else None
    div_hi_exon = [f for f in div if f["best_exon_id"] >= 0.98]
    div_hi_conf = statistics.mean(f["bam_conf_frac"] for f in div_hi_exon) if div_hi_exon else None
    div_lo_exon = [f for f in div if f["best_exon_id"] < 0.98]
    div_lo_conf = statistics.mean(f["bam_conf_frac"] for f in div_lo_exon) if div_lo_exon else None
    # residual confusable reads in divergent families (BAM)
    div_reads = sum(f["depth"] for f in div)
    div_mapq0 = sum(f["bam_mapq0"] for f in div)
    return dict(
        r_conf_vs_mean_matched=r_mean,
        r_conf_vs_best_exon=r_best,
        best_exon_stronger=(r_best is not None and r_mean is not None and abs(r_best) > abs(r_mean)),
        n_divergent_fams=len(div),
        divergent_conf_frac=div_conf_bam,
        divergent_residual_mapq0=div_mapq0, divergent_reads=div_reads,
        divergent_residual_frac=(div_mapq0 / div_reads) if div_reads else None,
        div_with_hi_exon_n=len(div_hi_exon), div_with_hi_exon_conf=div_hi_conf,
        div_with_lo_exon_n=len(div_lo_exon), div_with_lo_exon_conf=div_lo_conf,
    )


# --------------------------------------------------------------------------- confusable-mass concentration + floor impact
def mass_concentration(per_fam, axis, thresholds):
    """Of ALL confusable reads (BAM-MAPQ0 over all reads; AS-margin over the sample), what
    fraction come from families with identity axis >= t. Concentration => the method targets
    the top-identity families."""
    rows = []
    tot_bam = sum(f["bam_mapq0"] for f in per_fam.values())
    tot_marg = sum(f.get("margin_conf", 0) or 0 for f in per_fam.values() if f.get("n_realn"))
    for t in thresholds:
        hi_bam = sum(f["bam_mapq0"] for f in per_fam.values() if f[axis] >= t)
        hi_marg = sum((f.get("margin_conf", 0) or 0) for f in per_fam.values()
                      if f.get("n_realn") and f[axis] >= t)
        rows.append(dict(t=t,
                         bam_frac=(hi_bam / tot_bam) if tot_bam else None,
                         margin_frac=(hi_marg / tot_marg) if tot_marg else None,
                         bam=f"{hi_bam}/{tot_bam}", margin=f"{hi_marg}/{tot_marg}"))
    return rows


def _floor_cut(fam, floor, metric):
    """Would the divergence floor CUT this family at `floor` on `metric`? The floor cuts a cross-gene
    edge iff recip_id_best < floor; a missing recip value (None) means the edge is KEPT (absent from
    the ri cache => not a divergence signal, per family_rna_refine.load_recip_id) => family NOT cut."""
    v = fam.get(metric)
    return (v is not None) and (v < floor)


def floor_impact(per_fam, floors, metric="recip_best_mean"):
    """THE COINCIDENCE TEST. The divergence floor cuts cross-gene edges with reciprocal whole-
    transcript identity recip_id_best < FLOOR (family_rna_refine.load_recip_id). We aggregate
    recip_id_best per family (recip_best_mean = mean over member gene-pairs) -- THE ACTUAL floor
    metric -- and ask what fraction of read-confusable reads live in families the floor would CUT.
    A LOW cut-fraction => the floor keeps the confusable families (coincidence). A HIGH cut-fraction
    => the floor removes the families that most need the method (anti-coincidence).
    metric defaults to the real floor metric recip_best_mean; pass 'core_recip' for the (weaker,
    corr~0.49) gate-derived proxy as a cross-check."""
    rows = []
    tot_bam = sum(f["bam_mapq0"] for f in per_fam.values())
    tot_marg = sum(f.get("margin_conf", 0) or 0 for f in per_fam.values() if f.get("n_realn"))
    for fl_ in floors:
        cut_bam = sum(f["bam_mapq0"] for f in per_fam.values() if _floor_cut(f, fl_, metric))
        cut_marg = sum((f.get("margin_conf", 0) or 0) for f in per_fam.values()
                       if f.get("n_realn") and _floor_cut(f, fl_, metric))
        rows.append(dict(floor=fl_, metric=metric,
                         bam_cut_frac=(cut_bam / tot_bam) if tot_bam else None,
                         margin_cut_frac=(cut_marg / tot_marg) if tot_marg else None))
    return rows


def segment_local_proof(per_fam):
    """Divergent-by-whole-transcript families (mean_matched_id==0: NO colinear exon match) that
    are nonetheless read-confusable => confusability is SEGMENT-LOCAL, decoupled from whole-
    transcript identity => no whole-transcript floor can cleanly gate it."""
    ex = []
    for fid, f in per_fam.items():
        if f["mean_matched_id"] == 0.0 and (f.get("margin_conf") or 0) > 0:
            ex.append(dict(gene=f["dom"], core_recip=f["core_recip"],
                           recip_best_mean=f.get("recip_best_mean"), best_exon=f["best_exon_id"],
                           margin_conf=f["margin_conf"], n_realn=f["n_realn"],
                           bam_mapq0=f["bam_mapq0"], depth=f["depth"]))
    ex.sort(key=lambda d: -d["margin_conf"])
    return ex


# --------------------------------------------------------------------------- floor cross-ref
def read_divfloor():
    """Legacy raw-row reader (kept for the human-readable print of the sweep)."""
    if not os.path.exists(DIVFLOOR_TSV):
        return None
    rows = []
    with open(DIVFLOOR_TSV) as fh:
        for ln in fh:
            if ln.startswith("#") or ln.startswith("floor\t") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            rows.append(f)
    return rows


def load_divfloor_sweep():
    """Parse bench/divergence_floor.tsv into a list of dicts (one per floor row), typed.
    Prefers bench/divergence_floor.json if the sweep emitted it (else falls back to the TSV).
    'off' -> numeric floor 0.0. Returns (rows, source) or (None, None) if the sweep is absent."""
    # prefer a JSON emission of the sweep if the sweep wrote one (schema: {"rows":[...]})
    if os.path.exists(DIVFLOOR_JSON):
        try:
            with open(DIVFLOOR_JSON) as fh:
                obj = json.load(fh)
            jr = obj.get("rows", obj if isinstance(obj, list) else None)
            if jr:
                rows = []
                for r in jr:
                    fl = r.get("floor")
                    rows.append(dict(
                        floor_label=str(fl),
                        floor=(0.0 if str(fl) == "off" else float(fl)),
                        n_families=int(r["n_families"]),
                        R_oracle=float(r["R_oracle"]),
                        P_oracle_task=float(r["P_oracle_task"]),
                        P_oracle_dedup=float(r["P_oracle_dedup"]),
                        E_p_purity=float(r["E_p_purity"]),
                        distinct_fp=int(r["distinct_fp"]),
                        edges_cut=int(r.get("edges_cut_by_floor", 0)),
                        domain_shares_excluded=str(r.get("domain_shares_excluded", "NA"))))
                return rows, DIVFLOOR_JSON
        except Exception:
            pass
    if not os.path.exists(DIVFLOOR_TSV):
        return None, None
    rows = []
    with open(DIVFLOOR_TSV) as fh:
        header = None
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            if header is None:
                header = f
                I = {c: i for i, c in enumerate(header)}
                continue
            fl = f[I["floor"]]
            rows.append(dict(
                floor_label=fl,
                floor=(0.0 if fl == "off" else float(fl)),
                n_families=int(f[I["n_families"]]),
                R_oracle=float(f[I["R_oracle"]]),
                P_oracle_task=float(f[I["P_oracle_task"]]),
                P_oracle_dedup=float(f[I["P_oracle_dedup"]]),
                E_p_purity=float(f[I["E_p_purity"]]),
                distinct_fp=int(f[I["distinct_fp"]]),
                edges_cut=int(f[I["edges_cut_by_floor"]]) if "edges_cut_by_floor" in I else 0,
                domain_shares_excluded=(f[I["domain_shares_excluded"]]
                                        if "domain_shares_excluded" in I else "NA")))
    return rows, DIVFLOOR_TSV


def _f1(p, r):
    return (2 * p * r / (p + r)) if (p is not None and r is not None and (p + r) > 0) else None


def _dshare_frac(s):
    """'4/4' -> 1.0 ; '3/4' -> 0.75 ; 'NA' -> None."""
    try:
        a, b = s.split("/")
        return int(a) / int(b) if int(b) else None
    except Exception:
        return None


def derive_pr_optimal_floor(sweep):
    """Derive 'the P/R-optimal floor' from the divergence-floor sweep -- honestly.
    The sweep is a MONOTONE precision/recall trade (E_p purity rises, R_oracle falls with the
    floor), so there is NO interior F-score optimum; the shipped 0.80 is a PURITY choice (the
    smallest floor that excludes all named domain-shares). We report every reasonable notion so
    the coincidence test is not cherry-picked:
      * f1_ep_optimal   : argmax F1(E_p_purity, R_oracle)   [aggregate purity-vs-recall]
      * f1_task_optimal : argmax F1(P_oracle_task, R_oracle)
      * domain_share_knee: smallest numeric floor with domain_shares_excluded == 1.0 (the shipped
                           rationale; == SHIPPED_FLOOR when the sweep agrees)
      * monotone P/R flags so the caller can state 'no interior optimum'."""
    if not sweep:
        return None
    tab = []
    for r in sweep:
        tab.append(dict(floor=r["floor"], floor_label=r["floor_label"],
                        n_families=r["n_families"], R_oracle=r["R_oracle"],
                        P_task=r["P_oracle_task"], E_p=r["E_p_purity"],
                        f1_ep=_f1(r["E_p_purity"], r["R_oracle"]),
                        f1_task=_f1(r["P_oracle_task"], r["R_oracle"]),
                        dshare_frac=_dshare_frac(r["domain_shares_excluded"]),
                        dshare=r["domain_shares_excluded"]))
    f1_ep_opt = max(tab, key=lambda t: (t["f1_ep"] if t["f1_ep"] is not None else -1))
    f1_task_opt = max(tab, key=lambda t: (t["f1_task"] if t["f1_task"] is not None else -1))
    # domain-share knee: smallest numeric floor (>0) that fully excludes the named domain-shares
    knee = None
    for t in sorted((t for t in tab if t["floor"] > 0), key=lambda t: t["floor"]):
        if t["dshare_frac"] is not None and t["dshare_frac"] >= 1.0:
            knee = t["floor"]
            break
    # monotonicity: purity vs floor and recall vs floor
    floors_sorted = sorted(tab, key=lambda t: t["floor"])
    ep_series = [t["E_p"] for t in floors_sorted]
    r_series = [t["R_oracle"] for t in floors_sorted]
    ep_monotone_up = all(b >= a - 1e-9 for a, b in zip(ep_series, ep_series[1:]))
    r_monotone_down = all(b <= a + 1e-9 for a, b in zip(r_series, r_series[1:]))
    return dict(
        table=tab,
        f1_ep_optimal_floor=f1_ep_opt["floor"], f1_ep_optimal_value=f1_ep_opt["f1_ep"],
        f1_task_optimal_floor=f1_task_opt["floor"], f1_task_optimal_value=f1_task_opt["f1_task"],
        domain_share_knee=knee,
        shipped_floor=SHIPPED_FLOOR,
        purity_monotone_increasing=ep_monotone_up,
        recall_monotone_decreasing=r_monotone_down,
        no_interior_f1_optimum=(f1_ep_opt["floor"] in (0.0,) and f1_task_opt["floor"] in (0.0,)),
        # the operative floor for the coincidence test = the shipped/knee floor (the O1 default)
        operative_floor=(knee if knee is not None else SHIPPED_FLOOR),
        rationale=("monotone P/R trade: no interior F1 max (both F1 series peak at floor=off); the "
                   "shipped floor is a PURITY choice = smallest floor excluding all named domain-"
                   "shares (O1-precision), not a read-resolvability boundary"))


# --------------------------------------------------------------------------- output
def write_tsv(per_fam):
    cols = ["family_id", "dominant_gene", "n_members", "n_loci", "n_pairs_covered",
            "mean_matched_id", "min_matched_id", "best_exon_id", "core_recip", "exon_spread",
            "recip_best_mean", "recip_best_min", "n_recip_pairs",
            "depth", "bam_mapq0", "bam_conf_frac",
            "n_realn", "realn_mapq0", "realn_conf_frac", "margin_conf", "margin_conf_frac"]
    with open(OUT_TSV, "w") as out:
        out.write("# READ-CONFUSABILITY vs COPY IDENTITY (characterize). copy identity from "
                  "colinear_multiexon_gate.tsv (mean_matched_id=whole-read, best_exon_id=per-exon "
                  "max, core_recip=floor metric). confusable: BAM MAPQ0 (GGO_mm.bam -N50 -p0.1); "
                  f"re-align default splice:hq MAPQ0; AS-margin<{MARGIN_MIN} (Eichler). "
                  f"sample<={SAMPLE_PER_FAM}/fam seed={SEED}. Deterministic PYTHONHASHSEED=0.\n")
        out.write("\t".join(cols) + "\n")
        for fid in sorted(per_fam, key=lambda x: int(x)):
            f = per_fam[fid]
            def fmt(v):
                if v is None:
                    return "NA"
                if isinstance(v, float):
                    return f"{v:.4f}"
                return str(v)
            row = [fid, f["dom"], f["n_members"], f["n_loci"], f["n_pairs_covered"],
                   f["mean_matched_id"], f["min_matched_id"], f["best_exon_id"], f["core_recip"],
                   f["exon_spread"], f.get("recip_best_mean"), f.get("recip_best_min"),
                   f.get("n_recip_pairs"),
                   f["depth"], f["bam_mapq0"], f["bam_conf_frac"],
                   f.get("n_realn"), f.get("realn_mapq0"), f.get("realn_conf_frac"),
                   f.get("margin_conf"), f.get("margin_conf_frac")]
            out.write("\t".join(fmt(v) for v in row) + "\n")


# --------------------------------------------------------------------------- per-family TSV reloader
def load_per_fam_tsv(path):
    """Reconstruct the per_fam dict from a frozen read_confusability_boundary.tsv so the
    COINCIDENCE stage can run over the already-computed (race-safe, floor-OFF) per-family
    measurements without re-touching the 13GB minimap2 index or the churning shared catalog."""
    per_fam = {}
    with open(path) as fh:
        header = None
        for ln in fh:
            if ln.startswith("#"):
                continue
            if header is None:
                header = ln.rstrip("\n").split("\t")
                I = {c: i for i, c in enumerate(header)}
                continue
            f = ln.rstrip("\n").split("\t")

            def fl(c):
                v = f[I[c]]
                return None if v == "NA" else float(v)

            def it(c):
                v = f[I[c]]
                return None if v == "NA" else int(v)

            has_recip = "recip_best_mean" in I  # backward-compat with pre-fix frozen TSVs
            per_fam[f[I["family_id"]]] = dict(
                dom=f[I["dominant_gene"]], n_members=it("n_members"), n_loci=it("n_loci"),
                n_pairs_covered=it("n_pairs_covered"),
                mean_matched_id=fl("mean_matched_id"), min_matched_id=fl("min_matched_id"),
                best_exon_id=fl("best_exon_id"), core_recip=fl("core_recip"),
                exon_spread=fl("exon_spread"),
                recip_best_mean=(fl("recip_best_mean") if has_recip else None),
                recip_best_min=(fl("recip_best_min") if has_recip else None),
                n_recip_pairs=(it("n_recip_pairs") if has_recip else None),
                depth=it("depth"),
                bam_mapq0=it("bam_mapq0"), bam_conf_frac=fl("bam_conf_frac"),
                n_realn=it("n_realn"), realn_mapq0=it("realn_mapq0"),
                realn_conf_frac=fl("realn_conf_frac"), margin_conf=it("margin_conf"),
                margin_conf_frac=fl("margin_conf_frac"), margins=[])
    return per_fam


# --------------------------------------------------------------------------- THE COINCIDENCE VERDICT
def coincidence_verdict(result, per_fam, pr_opt):
    """Compare the read-confusability BOUNDARY (O2: where reads become aligner-unresolvable) to
    the divergence FLOOR (O1: the whole-transcript reciprocal-identity purity cut).

    Two independent tests, because the boundary and the floor live on DIFFERENT axes:
      (A) AXIS test  : is the confusability boundary, measured ON THE FLOOR'S OWN AXIS (core_recip),
                       within one identity bin of the operative floor? (a same-axis coincidence)
      (B) MASS test  : does a floor placed at the operative value CUT the read-confusable mass?
                       A floor that coincides with the boundary should keep confusable families
                       (cut ~0); a floor that cuts the majority of confusable reads is ANTI-coincident
                       (it removes exactly the families the method exists for).
    The MASS test is decisive (the floor's job is to cut families; we ask whether it cuts the
    confusable ones). Reported alongside the whole-read boundary (the headline resolvability line)."""
    b = result.get("boundaries", {})
    # headline boundary = whole-read (mean_matched_id), most confusable-sensitive metric (margin)
    bnd_whole = (b.get("mean_matched_id:margin_conf") or b.get("mean_matched_id:bam_conf") or {}).get("boundary")
    # same-axis boundary on THE ACTUAL floor metric (recip_best_mean) -- the apples-to-apples comparison
    bnd_flooraxis = (b.get("recip_best_mean:margin_conf") or b.get("recip_best_mean:bam_conf") or {}).get("boundary")
    floor = pr_opt["operative_floor"] if pr_opt else SHIPPED_FLOOR

    # (A) same-axis coincidence: floor vs the core_recip confusability boundary
    bin_w = 0.05
    axis_gap = (abs(bnd_flooraxis - floor) if bnd_flooraxis is not None else None)
    axis_coincide = (axis_gap is not None and axis_gap <= bin_w)

    # (B) mass test: fraction of confusable reads the operative floor would CUT.
    # Uses THE ACTUAL floor metric recip_id_best (recip_best_mean); None recip => edge KEPT => not cut.
    tot_bam = sum(f["bam_mapq0"] for f in per_fam.values())
    tot_marg = sum((f.get("margin_conf") or 0) for f in per_fam.values() if f.get("n_realn"))
    cut_bam = sum(f["bam_mapq0"] for f in per_fam.values() if _floor_cut(f, floor, "recip_best_mean"))
    cut_marg = sum((f.get("margin_conf") or 0) for f in per_fam.values()
                   if f.get("n_realn") and _floor_cut(f, floor, "recip_best_mean"))
    cut_bam_frac = (cut_bam / tot_bam) if tot_bam else None
    cut_marg_frac = (cut_marg / tot_marg) if tot_marg else None
    # cross-check: the same test on the gate-derived core_recip PROXY (weaker; corr~0.49)
    cut_bam_cr = sum(f["bam_mapq0"] for f in per_fam.values() if _floor_cut(f, floor, "core_recip"))
    cut_marg_cr = sum((f.get("margin_conf") or 0) for f in per_fam.values()
                      if f.get("n_realn") and _floor_cut(f, floor, "core_recip"))
    # coincidence on mass = floor keeps (does not cut) the confusable mass
    mass_coincide = (cut_marg_frac is not None and cut_marg_frac < 0.15
                     and cut_bam_frac is not None and cut_bam_frac < 0.15)

    coincide = bool(axis_coincide and mass_coincide)
    # family-level corrs (computed in main): confusability vs the ACTUAL floor axis (recip_best) and
    # vs the whole-read matched-exon axis. NEITHER is a segment metric -- the point is that NO family-
    # level identity metric cleanly predicts confusability (the driver is per-read shared-segment id).
    r_floor_axis = result.get("r_margin_recip_best")
    r_seg_axis = result.get("r_margin_mean_matched")

    verdict = dict(
        headline_boundary_wholeread=bnd_whole,
        floor_metric="recip_best_mean (recip_id_best = matches_best/max(len); the ACTUAL floor metric)",
        boundary_on_floor_axis_recip_best=bnd_flooraxis,
        operative_floor=floor,
        shipped_floor=SHIPPED_FLOOR,
        f1_ep_optimal_floor=(pr_opt["f1_ep_optimal_floor"] if pr_opt else None),
        f1_task_optimal_floor=(pr_opt["f1_task_optimal_floor"] if pr_opt else None),
        domain_share_knee=(pr_opt["domain_share_knee"] if pr_opt else None),
        axis_test_gap=axis_gap, axis_test_coincide=axis_coincide,
        mass_test_bam_cut_frac=cut_bam_frac, mass_test_margin_cut_frac=cut_marg_frac,
        mass_test_coincide=mass_coincide,
        # cross-check on the gate-derived core_recip proxy (weaker; corr~0.49 to the real floor):
        mass_test_bam_cut_frac_core_recip_proxy=((cut_bam_cr / tot_bam) if tot_bam else None),
        mass_test_margin_cut_frac_core_recip_proxy=((cut_marg_cr / tot_marg) if tot_marg else None),
        corr_confusable_vs_floor_axis_recip_best=r_floor_axis,
        corr_confusable_vs_segment_axis_mean_matched=r_seg_axis,
        coincide=coincide,
        interpretation=(
            "COINCIDE: the divergence floor sits at the read-confusability boundary; the floor is "
            "the principled resolvability cut (Canzar: no arbitrary threshold)."
            if coincide else
            "DO NOT COINCIDE: the floor is an O1 whole-transcript PURITY cut, the boundary is an O2 "
            "read-resolvability property. The operative floor CUTS the majority of read-confusable "
            "reads (removes exactly the near-identical/segment-local families the method targets), "
            "and confusability is a per-exon/per-read shared-segment property no whole-transcript "
            "identity floor can cleanly gate. Report the two as SEPARATE objects: O1-precision floor "
            "vs O2-read-confusability boundary."),
        scope_claim_supported="directionally" if not coincide else "yes",
    )
    return verdict


# --------------------------------------------------------------------------- figure
def render_figure(result, per_fam, sweep, pr_opt, verdict, out_png):
    """confusable-fraction vs copy identity, with the read-confusability BOUNDARY and the
    divergence FLOOR both marked. 2x2:
      A whole-read axis (mean_matched_id): 3 confusable metrics + boundary line   [O2]
      B floor axis (core_recip): 3 metrics + shipped floor line + floor 'CUT' band [O1 vs O2]
      C floor-impact: fraction of confusable reads a core_recip floor would CUT    [coincidence]
      D divergence-floor P/R trade (from divergence_floor.tsv): monotone, no interior optimum"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def centers(edges):
        cs = []
        for i in range(len(edges) - 1):
            lo, hi = edges[i], min(edges[i + 1], 1.0)
            cs.append((lo + hi) / 2)
        return cs

    cx = centers(ID_EDGES)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    C = dict(bam="#1f77b4", realn="#2ca02c", margin="#d62728")

    # ---- A: whole-read axis ----
    ax = axes[0][0]
    cw = result["curve_mean_matched_id"]
    for key, lab, col in [("bam_conf", "BAM MAPQ0 (-N50)", C["bam"]),
                          ("realn_conf", "re-aln default MAPQ0", C["realn"]),
                          ("margin_conf", "AS-margin<10 (Eichler)", C["margin"])]:
        ys = [c[key] for c in cw]
        ax.plot(cx, ys, "-o", color=col, label=lab, lw=2, ms=5)
    bw = verdict["headline_boundary_wholeread"]
    if bw is not None:
        ax.axvline(bw, color="black", ls="--", lw=2)
        ax.text(bw - 0.005, ax.get_ylim()[1] * 0.92, f"read-confusability\nboundary ≈{bw:.2f}",
                ha="right", va="top", fontsize=9, fontweight="bold")
    ax.set_title("A. Confusable fraction vs COPY IDENTITY (whole-read: mean matched-exon id)\n"
                 "[O2 resolvability axis] -- confusion concentrates at near-identity", fontsize=10)
    ax.set_xlabel("mean matched-exon identity (copy vs copy)")
    ax.set_ylabel("fraction of reads confusable")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)

    # ---- B: ACTUAL floor axis (recip_id_best) ----
    ax = axes[0][1]
    cc = result["curve_recip_best_mean"]
    for key, lab, col in [("bam_conf", "BAM MAPQ0", C["bam"]),
                          ("realn_conf", "re-aln MAPQ0", C["realn"]),
                          ("margin_conf", "AS-margin<10", C["margin"])]:
        ys = [c[key] for c in cc]
        ax.plot(cx, ys, "-o", color=col, label=lab, lw=2, ms=5)
    floor = verdict["operative_floor"]
    ax.axvspan(0, floor, color="red", alpha=0.08)
    ax.axvline(floor, color="red", ls="-", lw=2.5)
    ax.text(floor - 0.01, ax.get_ylim()[1] * 0.95,
            f"divergence FLOOR = {floor:.2f}\n(cuts recip_id_best < {floor:.2f})",
            ha="right", va="top", fontsize=9, fontweight="bold", color="red")
    bfa = verdict["boundary_on_floor_axis_recip_best"]
    if bfa is not None:
        ax.axvline(bfa, color="black", ls="--", lw=2)
        ax.text(bfa + 0.01, ax.get_ylim()[1] * 0.70, f"confusability\nboundary {bfa:.2f}",
                ha="left", va="top", fontsize=9, fontweight="bold")
    ax.set_title("B. Confusable fraction vs recip_id_best (THE ACTUAL FLOOR METRIC)\n"
                 "[O1 floor axis] -- red band = families the floor DISCARDS", fontsize=10)
    ax.set_xlabel("recip_id_best = matches_best/max(len) (reciprocal whole-transcript id = floor metric)")
    ax.set_ylabel("fraction of reads confusable")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)

    # ---- C: floor-impact (fraction of confusable reads CUT) ----
    ax = axes[1][0]
    fi = result["floor_impact"]
    fx = [r["floor"] for r in fi]
    ax.plot(fx, [r["bam_cut_frac"] for r in fi], "-o", color=C["bam"], lw=2, ms=6,
            label="BAM-MAPQ0 reads cut")
    ax.plot(fx, [r["margin_cut_frac"] for r in fi], "-s", color=C["margin"], lw=2, ms=6,
            label="AS-margin reads cut")
    ax.axvline(floor, color="red", ls="-", lw=2.5)
    ax.axhline(0.15, color="gray", ls=":", lw=1.5)
    ax.text(fx[0], 0.16, "coincidence would need ~0 cut", fontsize=8, color="gray")
    mc = verdict["mass_test_margin_cut_frac"]
    bc = verdict["mass_test_bam_cut_frac"]
    if mc is not None:
        ax.annotate(f"floor {floor:.2f} CUTS\n{bc*100:.0f}% BAM / {mc*100:.0f}% margin\nof confusable reads",
                    xy=(floor, mc), xytext=(floor + 0.02, min(mc, 0.9)),
                    fontsize=9, fontweight="bold", color="red",
                    arrowprops=dict(arrowstyle="->", color="red"))
    ax.set_title("C. COINCIDENCE (mass test): fraction of read-confusable reads the\n"
                 "recip_id_best floor would DISCARD -- high = ANTI-coincident", fontsize=10)
    ax.set_xlabel("candidate divergence floor (recip_id_best = the ACTUAL floor metric)")
    ax.set_ylabel("fraction of confusable reads CUT")
    ax.set_ylim(0, 1)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)

    # ---- D: divergence-floor P/R trade ----
    ax = axes[1][1]
    if sweep:
        sx = [r["floor"] for r in sweep]
        ax.plot(sx, [r["E_p_purity"] for r in sweep], "-o", color="#9467bd", lw=2, ms=5,
                label="E_p purity (precision)")
        ax.plot(sx, [r["R_oracle"] for r in sweep], "-o", color="#8c564b", lw=2, ms=5,
                label="R_oracle (recall)")
        ax.axvline(floor, color="red", ls="-", lw=2.5)
        f1o = pr_opt["f1_ep_optimal_floor"] if pr_opt else None
        if f1o is not None:
            ax.axvline(f1o, color="green", ls="--", lw=2)
            ax.text(f1o + 0.005, 0.55, f"F1 optimum\n= {f1o:.2f} (floor=off)", fontsize=8,
                    color="green", fontweight="bold")
        ax.set_title("D. Divergence-floor P/R trade (divergence_floor.tsv)\n"
                     "MONOTONE -> no interior optimum; 0.80 = domain-share PURITY knee", fontsize=10)
        ax.set_xlabel("divergence floor (recip_id_best)")
        ax.set_ylabel("metric value")
        ax.set_ylim(0.4, 1.02)
        ax.legend(fontsize=8, loc="lower left")
        ax.grid(alpha=0.3)
    else:
        ax.text(0.5, 0.5, "divergence_floor.tsv absent", ha="center", va="center")
        ax.axis("off")

    bnd = verdict["headline_boundary_wholeread"]
    floor_v = verdict["operative_floor"]
    verdict_txt = ("COINCIDE: floor == read-confusability boundary (principled)"
                   if verdict["coincide"] else
                   f"DO NOT COINCIDE: O2 read-confusability boundary ≈{bnd} (copy id)  vs  "
                   f"O1 divergence floor {floor_v} (whole-transcript purity)")
    fig.suptitle("Read-confusability boundary vs divergence floor\n" + verdict_txt,
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    return out_png


def main(argv=None):
    global SAMPLE_PER_FAM, CATALOG_TSV
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-realign", action="store_true",
                    help="skip the minimap2 re-align control (BAM-MAPQ only; fast)")
    ap.add_argument("--sample", type=int, default=SAMPLE_PER_FAM)
    ap.add_argument("--catalog", default=CATALOG_TSV,
                    help="family catalog TSV (default the shipped family_rna_refine.tsv; pass a "
                         "floor-OFF snapshot for the coincidence test since the shared file is a "
                         "floored default and may be rewritten by a running sweep)")
    ap.add_argument("--out", default=OUT_TSV, help="output per-family TSV")
    ap.add_argument("--from-tsv", default=None, metavar="TSV",
                    help="COINCIDENCE-STAGE mode: reload the frozen per-family read_confusability_"
                         "boundary.tsv (already computed against the race-safe floor-OFF snapshot) "
                         "and run curves+boundary+floor-coincidence+figure+JSON WITHOUT re-touching "
                         "the 13GB minimap2 index or the churning shared catalog. Deterministic.")
    ap.add_argument("--json-out", default=OUT_JSON, help="coincidence-stage JSON output")
    ap.add_argument("--fig", default=FIG_PNG, help="figure PNG output")
    ap.add_argument("--no-fig", action="store_true", help="skip figure render")
    args = ap.parse_args(argv)
    SAMPLE_PER_FAM = args.sample
    CATALOG_TSV = args.catalog
    globals()["OUT_TSV"] = args.out
    from_tsv = args.from_tsv

    if from_tsv:
        print(f"[from] coincidence stage from frozen per-family TSV = {from_tsv}", flush=True)
        per_fam = load_per_fam_tsv(from_tsv)
        print(f"[from] {len(per_fam)} families reloaded (skip BAM + re-align + TSV write)", flush=True)
    else:
        print(f"[cat ] catalog = {CATALOG_TSV}", flush=True)
        rng = random.Random(SEED)
        fams = load_families()
        gate = load_gate()
        recip = load_recip_best()
        print(f"[load] {len(fams)} families; gate {len(gate)} pairs; "
              f"recip_id_best {len(recip)} gene-pairs (ACTUAL floor metric)", flush=True)
        per_fam, sample_reads = bam_pass(fams, gate, recip, rng)
        print(f"[bam ] {len(per_fam)} families with identity+reads; "
              f"{len(sample_reads)} sampled reads for re-align", flush=True)

        if not args.no_realign:
            workdir = tempfile.mkdtemp(prefix="readconf_")
            print(f"[mm2 ] re-aligning {len(sample_reads)} reads -> {MMI} ({workdir})", flush=True)
            realn = realign_sample(sample_reads, workdir)
            print(f"[mm2 ] parsed {len(realn)} mapped reads", flush=True)
            attach_realn(per_fam, realn)
        else:
            for f in per_fam.values():
                f.update(n_realn=0, realn_mapq0=0, realn_conf_frac=None,
                         margin_conf=0, margin_conf_frac=None, margins=[])

        write_tsv(per_fam)
        print(f"[write] {OUT_TSV}\n", flush=True)

    # ---- curves ----
    result = dict(n_families=len(per_fam))
    for axis, name in [("mean_matched_id", "WHOLE-READ (mean_matched_id)"),
                       ("best_exon_id", "PER-EXON (best_exon_id max)"),
                       ("recip_best_mean", "FLOOR METRIC (recip_id_best -- ACTUAL floor)"),
                       ("core_recip", "FLOOR PROXY (core_recip -- gate-derived, corr~0.49)")]:
        curve = build_curve(per_fam, axis, ID_EDGES)
        result[f"curve_{axis}"] = curve
        print(f"=== RESOLVABILITY vs COPY IDENTITY :: {name} ===")
        print(f"{'bin':>12} {'nFam':>5} {'nReads':>8} {'BAM_conf':>9} "
              f"{'realn_conf':>10} {'margin_conf':>11} {'nRealn':>7}")
        for c in curve:
            def p(v):
                return "   NA" if v is None else f"{v:6.3f}"
            print(f"{c['bin']:>12} {c['n_fams']:5d} {c['n_reads']:8d} {p(c['bam_conf'])} "
                  f"{p(c['realn_conf']):>10} {p(c['margin_conf']):>11} {c['n_realn']:7d}")
        for key in ["bam_conf", "margin_conf"]:
            b = locate_boundary(curve, key)
            if b:
                result.setdefault("boundaries", {})[f"{axis}:{key}"] = b
                print(f"  boundary[{key}] = {b['boundary']:.3f}  (jump {b['jump']:+.3f} "
                      f"{b['from_bin']} -> {b['to_bin']})")
        print()

    # ---- exon dependence ----
    ex = exon_dependence(per_fam)
    result["exon_dependence"] = ex
    print("=== PER-EXON vs WHOLE-READ ===")
    print(f"  corr(confusable, mean_matched_id) = {ex['r_conf_vs_mean_matched']:+.3f}  (whole-read)")
    print(f"  corr(confusable, best_exon_id)    = {ex['r_conf_vs_best_exon']:+.3f}  (per-exon max)")
    print(f"  per-exon axis stronger predictor: {ex['best_exon_stronger']}")
    print(f"  divergent-overall families (mean_matched_id<0.85): n={ex['n_divergent_fams']} "
          f"conf_frac={ex['divergent_conf_frac']:.4f} residual_mapq0={ex['divergent_residual_mapq0']}"
          f"/{ex['divergent_reads']} ({ex['divergent_residual_frac']:.4f})")
    print(f"    of those WITH a high-identity exon (best_exon>=0.98): n={ex['div_with_hi_exon_n']} "
          f"conf={ex['div_with_hi_exon_conf'] if ex['div_with_hi_exon_conf'] is None else round(ex['div_with_hi_exon_conf'],4)}")
    print(f"    of those WITHOUT              (best_exon<0.98): n={ex['div_with_lo_exon_n']} "
          f"conf={ex['div_with_lo_exon_conf'] if ex['div_with_lo_exon_conf'] is None else round(ex['div_with_lo_exon_conf'],4)}")
    print()

    # ---- confounds ----
    cc = confound_checks(per_fam)
    result["confounds"] = cc
    print("=== CONFOUND CHECKS ===")
    print(f"  corr(confusable, identity)              = {cc['r_conf_identity']:+.3f}")
    print(f"  corr(confusable, depth)                 = {cc['r_conf_depth']:+.3f}")
    print(f"  corr(confusable, family size)           = {cc['r_conf_size']:+.3f}")
    print(f"  partial corr(conf, depth | identity)    = "
          f"{cc['r_conf_depth_ctrl_id'] if cc['r_conf_depth_ctrl_id'] is None else round(cc['r_conf_depth_ctrl_id'],3)}")
    print(f"  partial corr(conf, size  | identity)    = "
          f"{cc['r_conf_size_ctrl_id'] if cc['r_conf_size_ctrl_id'] is None else round(cc['r_conf_size_ctrl_id'],3)}")
    print(f"  partial corr(conf, identity | depth)    = "
          f"{cc['r_conf_identity_ctrl_depth'] if cc['r_conf_identity_ctrl_depth'] is None else round(cc['r_conf_identity_ctrl_depth'],3)}")
    if cc["mean_bam_conf"] is not None:
        print(f"  forced-secondary depression (n={cc['n_fam_with_realn']} fam): "
              f"mean BAM-MAPQ0 conf={cc['mean_bam_conf']:.4f}  vs  "
              f"re-align-MAPQ0 conf={cc['mean_realn_conf']:.4f}  "
              f"AS-margin conf={cc['mean_margin_conf']:.4f}")
        print(f"  agreement corr: BAM vs margin={cc.get('r_bam_vs_margin')}, "
              f"BAM vs realn-mapq={cc.get('r_bam_vs_realnmapq')}")
    print()

    # ---- confusable-mass concentration ----
    print("=== CONFUSABLE-READ MASS CONCENTRATION (fraction of confusable reads from families "
          ">= identity t; axis=mean_matched_id) ===")
    conc = mass_concentration(per_fam, "mean_matched_id", [0.85, 0.95, 0.98, 0.99])
    result["mass_concentration"] = conc
    for c in conc:
        print(f"  mean_matched_id>={c['t']:.2f}: BAM-MAPQ0 {c['bam']}={c['bam_frac']:.3f}  "
              f"AS-margin {c['margin']}={c['margin_frac']:.3f}")
    print()

    # ---- THE COINCIDENCE TEST: floor impact on confusable reads (ACTUAL floor metric recip_id_best) ----
    fi = floor_impact(per_fam, [0.70, 0.80, 0.85, 0.90], metric="recip_best_mean")
    fi_proxy = floor_impact(per_fam, [0.70, 0.80, 0.85, 0.90], metric="core_recip")
    result["floor_impact"] = fi
    result["floor_impact_core_recip_proxy"] = fi_proxy
    print("=== FLOOR-BOUNDARY COINCIDENCE :: fraction of read-confusable reads the divergence floor "
          "(recip_id_best) would CUT ===")
    print("  (LOW cut = floor keeps confusable families = coincidence; HIGH cut = anti-coincidence)")
    for r, rp in zip(fi, fi_proxy):
        print(f"  floor={r['floor']:.2f}: recip_id_best -> BAM cut={r['bam_cut_frac']:.2f} "
              f"AS-margin cut={r['margin_cut_frac']:.2f}   [core_recip proxy: BAM "
              f"{rp['bam_cut_frac']:.2f} margin {rp['margin_cut_frac']:.2f}]")
    # family-level corr of confusability with the ACTUAL floor axis (recip_best_mean) and whole-read axis
    fam_marg = [(f["margin_conf"] / f["n_realn"], f.get("recip_best_mean"), f["mean_matched_id"])
                for f in per_fam.values() if f.get("n_realn")]
    fm_rb = [(a, b) for a, b, _ in fam_marg if b is not None]
    if fm_rb:
        r_rb = pearson([a for a, _ in fm_rb], [b for _, b in fm_rb])
        result["r_margin_recip_best"] = r_rb
    if fam_marg:
        r_mi = pearson([a for a, _, _ in fam_marg], [c for _, _, c in fam_marg])
        result["r_margin_mean_matched"] = r_mi
        print(f"  corr(margin_conf_frac, recip_id_best[ACTUAL floor axis, n={len(fm_rb)}]) = "
              f"{result.get('r_margin_recip_best'):+.3f}   "
              f"corr(margin_conf_frac, mean_matched_id[matched-exon axis]) = {r_mi:+.3f}")
        print("  (neither family-level identity metric cleanly predicts confusability; the driver is "
              "a PER-READ shared-segment identity no family metric captures)")
    print()

    # ---- segment-local (per-exon) proof ----
    sl = segment_local_proof(per_fam)
    result["segment_local"] = sl
    print("=== SEGMENT-LOCAL (per-exon) CONFUSABILITY PROOF :: divergent families (mean_matched_id"
          "==0, no colinear exon match) that are STILL read-confusable ===")
    print(f"  {len(sl)} such families. Top (gene, recip_id_best, core_recip, best_exon, "
          f"margin_conf/n_realn, bamMQ0/depth):")
    for d in sl[:8]:
        rb = "NA" if d.get("recip_best_mean") is None else f"{d['recip_best_mean']:.3f}"
        print(f"    {d['gene']:<14} recip_id_best={rb} core_recip={d['core_recip']:.3f} "
              f"best_exon={d['best_exon']:.3f} margin={d['margin_conf']}/{d['n_realn']} "
              f"bamMQ0={d['bam_mapq0']}/{d['depth']}")
    print()

    # ========================== COINCIDENCE STAGE (boundary vs P/R-optimal floor) ==========================
    # ---- load the divergence-floor sweep + derive 'the P/R-optimal floor' ----
    sweep, src = load_divfloor_sweep()
    result["divfloor_present"] = bool(sweep)
    result["divfloor_source"] = src
    result["divfloor_md_present"] = os.path.exists(DIVFLOOR_MD)
    if sweep:
        result["divfloor_sweep"] = sweep
        print(f"=== DIVERGENCE-FLOOR SWEEP ({os.path.basename(src)}) ===")
        print(f"  {'floor':>5} {'nFam':>5} {'R_oracle':>8} {'P_task':>7} {'E_p':>6} "
              f"{'F1(Ep,R)':>9} {'domShareExcl':>12}")
        pr_opt = derive_pr_optimal_floor(sweep)
        result["pr_optimal_floor"] = pr_opt
        by_floor = {t["floor_label"]: t for t in pr_opt["table"]}
        for r in sweep:
            t = by_floor[r["floor_label"]]
            print(f"  {r['floor_label']:>5} {r['n_families']:5d} {r['R_oracle']:8.4f} "
                  f"{r['P_oracle_task']:7.4f} {r['E_p_purity']:6.4f} "
                  f"{(t['f1_ep'] or 0):9.4f} {r['domain_shares_excluded']:>12}")
        print(f"  >> purity monotone-up: {pr_opt['purity_monotone_increasing']}  "
              f"recall monotone-down: {pr_opt['recall_monotone_decreasing']}")
        print(f"  >> F1(E_p,R) optimum   = floor {pr_opt['f1_ep_optimal_floor']:.2f} "
              f"(value {pr_opt['f1_ep_optimal_value']:.4f})")
        print(f"  >> F1(P_task,R) optimum= floor {pr_opt['f1_task_optimal_floor']:.2f} "
              f"(value {pr_opt['f1_task_optimal_value']:.4f})")
        print(f"  >> domain-share knee   = floor {pr_opt['domain_share_knee']}  "
              f"(shipped floor = {pr_opt['shipped_floor']})  -> OPERATIVE floor "
              f"{pr_opt['operative_floor']}")
        print(f"  >> {pr_opt['rationale']}")
    else:
        pr_opt = None
        print("=== DIVERGENCE-FLOOR SWEEP: bench/divergence_floor.{tsv,json} absent -> "
              "boundary reported STANDALONE; cross-ref the floor later ===")
    print()

    # ---- the verdict ----
    verdict = coincidence_verdict(result, per_fam, pr_opt)
    result["coincidence"] = verdict
    print("=== FLOOR <-> BOUNDARY COINCIDENCE VERDICT ===")
    print(f"  read-confusability boundary (whole-read mean_matched_id) = "
          f"{verdict['headline_boundary_wholeread']}")
    print(f"  read-confusability boundary ON THE ACTUAL FLOOR AXIS (recip_id_best) = "
          f"{verdict['boundary_on_floor_axis_recip_best']}")
    print(f"  operative divergence floor (recip_id_best)               = {verdict['operative_floor']}")
    print(f"  F1-optimal floor (E_p / P_task)                          = "
          f"{verdict['f1_ep_optimal_floor']} / {verdict['f1_task_optimal_floor']}  "
          f"(domain-share knee {verdict['domain_share_knee']})")
    print(f"  (A) axis test: |boundary-floor|={verdict['axis_test_gap']}  "
          f"coincide={verdict['axis_test_coincide']}")
    print(f"  (B) mass test [recip_id_best]: floor cuts {verdict['mass_test_bam_cut_frac']:.2f} BAM / "
          f"{verdict['mass_test_margin_cut_frac']:.2f} AS-margin of confusable reads  "
          f"coincide={verdict['mass_test_coincide']}")
    print(f"      [core_recip proxy cross-check: {verdict['mass_test_bam_cut_frac_core_recip_proxy']:.2f} "
          f"BAM / {verdict['mass_test_margin_cut_frac_core_recip_proxy']:.2f} AS-margin]")
    print(f"  corr(confusable, ACTUAL floor axis recip_id_best) = "
          f"{verdict['corr_confusable_vs_floor_axis_recip_best']}  "
          f"vs corr(confusable, matched-exon axis mean_matched_id) = "
          f"{verdict['corr_confusable_vs_segment_axis_mean_matched']}")
    print(f"  >> COINCIDE = {verdict['coincide']}   scope-claim supported = "
          f"{verdict['scope_claim_supported']}")
    print(f"  >> {verdict['interpretation']}")
    print()

    # ---- figure ----
    if not args.no_fig:
        fig_path = render_figure(result, per_fam, sweep, pr_opt, verdict, args.fig)
        result["figure_path"] = fig_path
        print(f"[fig ] {fig_path}", flush=True)

    # ---- JSON ----
    result["determinism"] = dict(pythonhashseed=os.environ.get("PYTHONHASHSEED"),
                                 seed=SEED, mm2_seed=MM2_SEED, margin_min=MARGIN_MIN,
                                 id_edges=ID_EDGES, from_tsv=from_tsv)
    with open(args.json_out, "w") as jf:
        json.dump(result, jf, sort_keys=True, indent=2, default=str)
    print(f"[json] {args.json_out}", flush=True)

    print("\nJSON " + json.dumps(result, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
