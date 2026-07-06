#!/usr/bin/env python3
"""o2_vg_visualization.py -- MATERIALIZE the O2 copy-assignment VARIATION GRAPH and render it.

THE CANZAR-FLIP FRAMING (the thesis's distinctive claim, made concrete here):
  In copy-assignment (objective O2), a multi-copy gene family is ONE variation graph.
    * COPIES are PATHS through that graph.
    * PSVs (paralogous sequence variants) are the BUBBLES -- the positions where copies differ.
    * A SUN (Singly Unique Nucleotide, Sudmant 2010) is a bubble where ONE copy's allele is
      PRIVATE -> a single read over it pins that copy deterministically (Strong-Separation witness).
    * multi-mapping READS are SHARED EVIDENCE threaded through the graph: a read follows the
      partial path consistent with its observed alleles and is ASSIGNED to that copy-path via the
      significance gate (n_decisive >= 1 distinguishing bubble AND a log-likelihood-ratio margin,
      copy_assign.assign_read -- the per-read min_p certificate of THEORY.md Thm 4 / Thm 7).
    * a read that spans NO distinguishing bubble stays TIED -- the K=0 / K-frontier identifiability
      FLOOR (MCC = chi(H); copies with identical hap-rows collapse to one node -> no bubble -> no
      de-tie is possible even in principle).

This script does NOT re-derive the biology from scratch: it reuses the exact catalog machinery of
bench/psv_graph_genomewide.py (the 154 validated co-located families -> backbone, asm20 copy
alignment, splice:hq read alignment, read-supported PSV columns) and the REAL per-read gate
bench/copy_assign.assign_read. It then MATERIALIZES the VG as an explicit object (nodes / bubbles /
copy-paths / read-threadings / assignment) and RENDERS it (matplotlib backbone+bubbles+paths+reads,
plus an optional networkx graph view).

FLAGSHIP FAMILIES (validated against the real data, see O2_VG_VISUALIZATION.md):
  * fam 39  RABL2 (RABL2A + RABL2B + 3 LOC paralogs, 5 copies, cross-chromosomal)
            -> RESOLVABLE end: K = #copies, every copy Tier-1 SUN-identifiable, reads cleanly assigned.
  * fam 1   RGPD8 / RANBP2 paralog cluster (7 copies, one chromosome)
            -> K-FRONTIER FLOOR: 0 read-supported PSV bubbles, K=1, every read TIED.
  * fam 22  ANKRD18A / ANKRD18B cluster (6 copies)
            -> MID / within-graph frontier: some copies resolve (private bubbles), some collapse
               (shared hap-rows -> converging paths).

Outputs (per family <F> = family id):
  bench/fig_o2_vg_fam<F>_<gene>.png        matplotlib backbone+bubbles+paths+reads figure (REQUIRED)
  bench/fig_o2_vg_fam<F>_<gene>_graph.png  networkx bubble-graph view (optional, if networkx present)
  /home/juanfra/winloci_scratch/o2vg/fam<F>.json   the materialized VG object (reusable cache)

Reusable API:
  vg = materialize_family(fam_id)   -> dict: nodes, bubbles, copy_paths, reads, assignment, stats
  render_family(vg, out_png)        -> writes the figure

Deterministic (PYTHONHASHSEED=0 enforced in __main__; fixed layout, no RNG in layout).
Run:  PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/o2_vg_visualization.py [fam_ids...]
"""
import collections
import json
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import psv_graph_genomewide as pg   # noqa: E402  (reuse sh / column_alleles / dedup_copies / consts)
import copy_assign as ca            # noqa: E402  (reuse the REAL per-read gate assign_read)
import recombinant_abstain as ra    # noqa: E402  (DEFAULT-ON recombinant-read ABSTAIN leg of the gate)

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = "/home/juanfra/winloci_scratch/o2vg"
SUN_TSV = os.path.join(HERE, "sun_identifiability.tsv")

# ----- flagship roster (family id -> short gene label used in filenames/titles) -----
FLAGSHIP = {
    39: ("RABL2", "RESOLVABLE"),   # 5 copies, all Tier-1 SUN, cross-chrom
    1:  ("RGPD8", "FLOOR"),        # 7 copies, 0 PSV bubbles, all reads tied
    22: ("ANKRD18", "MID"),        # 6 copies, K=4/6, some collapse
}


# --------------------------------------------------------------------------- gene labels
def load_gene_labels():
    """(chrom,start) -> gene name, from sun_identifiability.tsv (same region frame as here)."""
    g = {}
    if os.path.exists(SUN_TSV):
        with open(SUN_TSV) as f:
            next(f)
            for line in f:
                p = line.rstrip("\n").split("\t")
                g[(p[2], int(p[3]))] = p[5]
    return g


# --------------------------------------------------------------------------- families table
def load_families():
    fams = collections.defaultdict(list)
    with open(pg.FAM_TSV) as f:
        next(f)
        for line in f:
            fi, lid, c, s, e, nr = line.rstrip("\n").split("\t")
            fams[int(fi)].append((lid, c, int(s), int(e), int(nr)))
    return fams


# --------------------------------------------------------------------------- MATERIALIZE the VG
def materialize_family(fam_id, use_cache=True):
    """Build the per-family variation-graph object from the REAL data. Returns a JSON-able dict:

      family, copies:[{name, chrom, start, end, len, gene, hap_group, tier, n_sun, sun_cols,
                       n_assigned_reads}],
      backbone_len,
      bubbles:[{col, pos, alleles:{allele:[copy,...]}, is_sun, sun_copy, n_read_support}],
      reads:[{name, best_copy|None, status(assigned|ambiguous|tied|no_cover), margin, n_span,
              n_decisive, span_cols:[first,last], obs:{col:base}}],
      assignment: Counter, per_copy: {copy: n}, K, cls, n_reads, tiers:Counter
    Assignment is the REAL gate copy_assign.assign_read (n_decisive>=1 AND margin>=MARGIN -> ASSIGNED,
    n_decisive>=1 but margin<MARGIN -> AMBIGUOUS, n_decisive==0 -> TIED, no PSV covered -> no_cover).
    """
    os.makedirs(CACHE, exist_ok=True)
    cpath = os.path.join(CACHE, f"fam{fam_id}.json")
    if use_cache and os.path.exists(cpath):
        vg = json.load(open(cpath))
        vg["tiers"] = {int(k): v for k, v in vg["tiers"].items()}   # JSON stringifies int keys
        ra.apply_abstain_to_vg(vg)      # DEFAULT-ON abstain leg (no-op under RUSTLE_NO_RECOMBINANT_ABSTAIN=1)
        return vg

    tmp = os.path.join(CACHE, f"tmp_fam{fam_id}")
    os.makedirs(tmp, exist_ok=True)
    genome = pysam.FastaFile(pg.GENOME)
    contig_len = dict(zip(genome.references, genome.lengths))
    gene_lbl = load_gene_labels()
    fams = load_families()

    regions = pg.dedup_copies(fams[fam_id])
    regions = sorted(regions, key=lambda r: -(r[2] - r[1]))
    names = [f"copy{i}" for i in range(len(regions))]
    bc, bs, be, _ = regions[0]
    bs, be = max(0, bs), min(contig_len.get(bc, be), be)
    refseq = genome.fetch(bc, bs, be).upper()

    open(f"{tmp}/ref.fa", "w").write(f">{names[0]}\n{refseq}\n")
    pg.sh(f"samtools faidx {tmp}/ref.fa")
    with open(f"{tmp}/copies.fa", "w") as o:
        for nm, (c, s, e, _) in list(zip(names, regions))[1:]:
            s2, e2 = max(0, s), min(contig_len.get(c, e), e)
            o.write(f">{nm}\n{genome.fetch(c, s2, e2).upper()}\n")
    pg.sh(f"minimap2 -ax asm20 {tmp}/ref.fa {tmp}/copies.fa 2>/dev/null "
          f"| samtools sort -o {tmp}/copies.bam - && samtools index {tmp}/copies.bam")

    # reads from all copy loci -> best-mapping copy (lowest de); the multimapper evidence pool
    bam = pysam.AlignmentFile(pg.BAM, "rb")
    best, seqs = {}, {}
    for i, (c, s, e, _) in enumerate(regions):
        s2, e2 = max(0, s - pg.PAD), e + pg.PAD
        for r in bam.fetch(c, s2, min(contig_len.get(c, e2), e2)):
            if r.is_unmapped or r.is_supplementary:
                continue
            if r.reference_end is None or r.reference_end < s or r.reference_start > e:
                continue
            de = dict(r.get_tags()).get("de")
            if de is None:
                continue
            if r.query_name not in best or de < best[r.query_name][1]:
                best[r.query_name] = (i, de)
            if not r.is_secondary and r.query_sequence and r.query_name not in seqs:
                seqs[r.query_name] = r.query_sequence
    bam.close()
    if len(seqs) > pg.READ_CAP:
        seqs = dict(list(seqs.items())[:pg.READ_CAP])
    with open(f"{tmp}/reads.fa", "w") as o:
        for n, sq in seqs.items():
            o.write(f">{n}\n{sq}\n")
    pg.sh(f"minimap2 -ax splice:hq {tmp}/ref.fa {tmp}/reads.fa 2>/dev/null "
          f"| samtools sort -o {tmp}/reads.bam - && samtools index {tmp}/reads.bam")

    copy_cols = pg.column_alleles(f"{tmp}/copies.bam", names[0])
    read_cols = pg.column_alleles(f"{tmp}/reads.bam", names[0])

    # PSV bubbles: all copies defined, >=2 distinct alleles, read-supported (>=MIN_READS_PSV)
    bubbles = []
    for p in range(len(refseq)):
        hap = {names[0]: refseq[p]}
        for nm in names[1:]:
            if nm in copy_cols.get(p, {}):
                hap[nm] = copy_cols[p][nm]
        if (len(hap) == len(names) and len(set(hap.values())) >= 2
                and len(read_cols.get(p, {})) >= pg.MIN_READS_PSV):
            alleles = collections.defaultdict(list)
            for nm in names:
                alleles[hap[nm]].append(nm)
            vals = list(hap.values())
            sun_copy = None
            for nm in names:
                if vals.count(hap[nm]) == 1:
                    sun_copy = nm  # a private allele at this column (report one representative)
                    break
            bubbles.append(dict(col=len(bubbles), pos=p, hap=hap,
                                alleles={a: v for a, v in alleles.items()},
                                is_sun=sun_copy is not None, sun_copy=sun_copy,
                                n_read_support=len(read_cols.get(p, {}))))
    hap_rows = {nm: tuple(b["hap"][nm] for b in bubbles) for nm in names}
    groups = collections.defaultdict(list)
    for nm in names:
        groups[hap_rows[nm]].append(nm)
    group_id = {}
    for gi, (hr, mem) in enumerate(sorted(groups.items(), key=lambda kv: -len(kv[1]))):
        for nm in mem:
            group_id[nm] = gi

    # per-copy SUN columns + tier (recomputed here, internally consistent -- NOT joined from tsv)
    sun_cols = collections.defaultdict(list)
    for j, b in enumerate(bubbles):
        vals = list(b["hap"].values())
        for nm in names:
            if vals.count(b["hap"][nm]) == 1:
                sun_cols[nm].append(j)
    tier = {}
    for nm in names:
        if len(groups[hap_rows[nm]]) > 1:
            tier[nm] = 3                       # hap-row collapses with another copy -> K-frontier
        elif sun_cols.get(nm):
            tier[nm] = 1                       # >=1 private column -> single-read SUN-identifiable
        else:
            tier[nm] = 2                       # unique hap-row but no SUN -> needs PSV linkage

    # copy allele vectors in COLUMN space (for the real gate)
    copy_vecs = {nm: {j: bubbles[j]["hap"][nm] for j in range(len(bubbles))} for nm in names}

    # thread every read through the graph via the REAL gate
    reads_out = []
    cats = collections.Counter()
    per_copy = collections.Counter()
    for r in sorted(best):
        obs = {}
        for j, b in enumerate(bubbles):
            base = read_cols.get(b["pos"], {}).get(r)
            if base is not None:
                obs[j] = base
        if not obs:
            cats["no_cover"] += 1
            reads_out.append(dict(name=r, best_copy=None, status="no_cover", margin=0.0,
                                  n_span=0, n_decisive=0, span_cols=None, obs={}))
            continue
        b_copy, margin, n_span, n_dec, _ = ca.assign_read(obs, copy_vecs)
        spanned = sorted(obs)
        span = [spanned[0], spanned[-1]]
        if n_dec >= 1 and margin >= ca.MARGIN:
            status, bc_ = "assigned", b_copy
            cats["assigned"] += 1
            per_copy[b_copy] += 1
        elif n_dec >= 1:
            status, bc_ = "ambiguous", None
            cats["ambiguous"] += 1
        else:
            status, bc_ = "tied", None
            cats["tied"] += 1
        reads_out.append(dict(name=r, best_copy=bc_, status=status, margin=round(float(margin), 2),
                              n_span=n_span, n_decisive=n_dec, span_cols=span,
                              obs={str(k): v for k, v in obs.items()}))

    K = len(groups)
    if len(bubbles) == 0:
        cls = "no_psv (K-FRONTIER floor)"
    elif K == 1:
        cls = "collapsed (K-FRONTIER floor)"
    elif K == len(names):
        cls = "fully_resolvable"
    else:
        cls = "partial (mixed resolvable / collapsed)"

    copies_out = []
    for i, (c, s, e, nr) in enumerate(regions):
        nm = names[i]
        copies_out.append(dict(name=nm, chrom=c, start=s, end=e, length=e - s,
                               gene=gene_lbl.get((c, s), "?"), hap_group=group_id[nm],
                               tier=tier[nm], n_sun=len(sun_cols.get(nm, [])),
                               sun_cols=sun_cols.get(nm, []), n_assigned_reads=per_copy.get(nm, 0)))

    vg = dict(
        family=fam_id, gene=FLAGSHIP.get(fam_id, ("?", "?"))[0],
        role=FLAGSHIP.get(fam_id, ("?", "?"))[1],
        n_copies=len(names), backbone_len=len(refseq), backbone_chrom=bc,
        backbone_start=bs, backbone_end=be,
        cross_chrom=len({r[0] for r in regions}) >= 2,
        chroms=sorted({r[0] for r in regions}),
        n_bubbles=len(bubbles), n_sun_bubbles=sum(1 for b in bubbles if b["is_sun"]),
        K=K, cls=cls, group_sizes=sorted((len(v) for v in groups.values()), reverse=True),
        tiers=dict(collections.Counter(tier.values())),
        copies=copies_out, bubbles=bubbles, reads=reads_out,
        n_reads=len(best), assignment=dict(cats), per_copy=dict(per_copy),
    )
    genome.close()
    json.dump(vg, open(cpath, "w"))     # CACHE stays the PURE-gate assignment (abstain leg applied on read)
    ra.apply_abstain_to_vg(vg)          # DEFAULT-ON abstain leg (no-op under RUSTLE_NO_RECOMBINANT_ABSTAIN=1)
    return vg


# --------------------------------------------------------------------------- RENDER (matplotlib)
COPY_COLORS = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e",
               "#17becf", "#8c564b", "#e377c2", "#bcbd22", "#7f7f7f"]
ALLELE_Y = {"A": 0.30, "C": 0.10, "G": -0.10, "T": -0.30, "N": 0.0}


def _pick_bubbles(vg, max_show=30):
    """Choose a legible, deterministic subset of bubbles: all SUN bubbles kept preferentially,
    then fill with evenly-spaced non-SUN bubbles across the backbone. Returns sorted col indices."""
    bubbles = vg["bubbles"]
    if len(bubbles) <= max_show:
        return list(range(len(bubbles)))
    sun = [b["col"] for b in bubbles if b["is_sun"]]
    non = [b["col"] for b in bubbles if not b["is_sun"]]
    keep = set()
    # evenly sample SUN bubbles first (they carry the identifiability story)
    if sun:
        step = max(1, len(sun) / min(len(sun), max_show))
        keep.update(sun[int(i * step)] for i in range(min(len(sun), max_show)))
    # fill remaining slots with evenly-spaced non-SUN bubbles
    rem = max_show - len(keep)
    if rem > 0 and non:
        step = max(1, len(non) / rem)
        keep.update(non[int(i * step)] for i in range(rem))
    return sorted(keep)[:max_show]


def render_family(vg, out_png, max_bubbles=30, max_reads=60):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyBboxPatch

    names = [c["name"] for c in vg["copies"]]
    gene_of = {c["name"]: c["gene"] for c in vg["copies"]}
    tier_of = {c["name"]: c["tier"] for c in vg["copies"]}
    color = {nm: COPY_COLORS[i % len(COPY_COLORS)] for i, nm in enumerate(names)}
    bubbles = vg["bubbles"]
    blen = max(vg["backbone_len"], 1)

    shown = _pick_bubbles(vg, max_bubbles)
    # lay shown bubbles out on an EVEN x-grid (readability), but keep a position ruler on top
    n_show = max(len(shown), 1)
    xgrid = {col: (0.06 + 0.90 * (k / max(n_show - 1, 1))) if n_show > 1 else 0.5
             for k, col in enumerate(shown)}
    col2x = lambda c: xgrid.get(c)

    fig = plt.figure(figsize=(16, 8.2), dpi=200)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.02], hspace=0.02,
                          top=0.90, bottom=0.10, left=0.015, right=0.985)
    axr = fig.add_subplot(gs[0])   # reads band (top)
    axg = fig.add_subplot(gs[1])   # variation graph (bottom)
    for ax in (axr, axg):
        ax.set_xlim(0, 1)
        ax.axis("off")

    # ---------------- bottom: the VARIATION GRAPH (backbone + bubbles + copy-paths) ----------
    BB_Y = 0.0
    axg.set_ylim(-0.78, 0.80)
    # backbone bar (shared consensus path)
    axg.add_patch(FancyBboxPatch((0.02, BB_Y - 0.045), 0.96, 0.09,
                                 boxstyle="round,pad=0.002", linewidth=0,
                                 facecolor="#dddddd", zorder=1))
    bb_label_y = 0.62 if bubbles else 0.20   # placed clear of bubbles / floor ylim

    if not bubbles:
        axg.set_ylim(-0.62, 0.28)                            # tighten (no bubbles above backbone)
        # ---- K-FRONTIER FLOOR: no bubble at all; every copy-path rides the backbone identically
        for i, nm in enumerate(names):
            y = BB_Y + (i - (len(names) - 1) / 2) * 0.018   # tiny fan just to show N overlapping paths
            axg.plot([0.03, 0.97], [y, y], color=color[nm], lw=1.8, alpha=0.9, zorder=3)
        axg.text(0.5, -0.42,
                 f"0 read-supported PSV bubbles  ->  all {len(names)} copy-paths COLLAPSE onto the\n"
                 f"backbone (identical hap-rows, MCC = chi(H) = 1)  ->  no distinguishing position\n"
                 f"exists  ->  every read is TIED. This is the K-frontier identifiability FLOOR.",
                 fontsize=13, ha="center", va="center", color="#333333",
                 bbox=dict(boxstyle="round,pad=0.6", fc="#f4f4f4", ec="#bbbbbb"))
    else:
        # draw shown bubbles: alleles as sublevel nodes; copy-paths thread their allele sublevel
        prev_y = {nm: BB_Y for nm in names}
        prev_x = {nm: 0.03 for nm in names}
        for k, col in enumerate(shown):
            b = bubbles[col]
            x = col2x(col)
            hap = b["hap"]
            alleles = b["alleles"]
            # allele nodes
            for a, members in alleles.items():
                ya = BB_Y + ALLELE_Y.get(a, 0.0)
                is_priv = len(members) == 1
                axg.add_patch(FancyBboxPatch((x - 0.007, ya - 0.03), 0.014, 0.06,
                              boxstyle="round,pad=0.002", linewidth=(1.4 if is_priv else 0.6),
                              edgecolor=("#000000" if is_priv else "#999999"),
                              facecolor=("#fff2b2" if is_priv else "#ffffff"), zorder=4))
                axg.text(x, ya, a, fontsize=6.5, ha="center", va="center", zorder=5,
                         fontweight=("bold" if is_priv else "normal"))
            # SUN marker
            if b["is_sun"]:
                axg.plot([x], [BB_Y + 0.44], marker="*", ms=8, color="#c8a200", zorder=6)
            # thread each copy-path to its allele sublevel at this bubble
            for nm in names:
                ya = BB_Y + ALLELE_Y.get(hap[nm], 0.0)
                axg.plot([prev_x[nm], x], [prev_y[nm], ya], color=color[nm], lw=1.3,
                         alpha=0.8, zorder=3)
                prev_y[nm], prev_x[nm] = ya, x
        for nm in names:  # exit to backbone end
            axg.plot([prev_x[nm], 0.97], [prev_y[nm], BB_Y], color=color[nm], lw=1.3,
                     alpha=0.8, zorder=3)
        axg.text(0.5, 0.62, f"PSV BUBBLES  (yellow/bold = SUN: one copy's PRIVATE allele "
                 f"-> single-read resolvable)   *  = SUN bubble", fontsize=8.5,
                 ha="center", va="bottom", color="#555555")
    axg.text(0.008, bb_label_y, "backbone = shared consensus path (gray bar)", fontsize=8,
             color="#888888", ha="left", va="center")

    # ---------------- top: READS threaded as arcs, colored by ASSIGNED copy / gray if unresolved --
    axr.set_ylim(0, 1)
    reads = vg["reads"]
    if not bubbles:
        # FLOOR family: no bubble exists -> no read can be threaded through a distinguishing
        # position; draw a stack of gray bars = every read TIED at the K-frontier floor.
        n_draw = min(len(reads), max_reads)
        step = max(1, len(reads) // max(n_draw, 1))
        sample = reads[::step][:n_draw]
        dy = 0.9 / max(len(sample), 1)
        for ri, _ in enumerate(sample):
            y = 0.95 - ri * dy
            # deterministic pseudo-span so bars vary in length (visual only; all are TIED)
            frac = 0.55 + 0.42 * (((ri * 2654435761) % 1000) / 1000.0)
            x0 = 0.04 + 0.06 * (((ri * 40503) % 1000) / 1000.0)
            axr.plot([x0, x0 + frac * 0.9], [y, y], color="#bdbdbd", lw=1.0,
                     solid_capstyle="round", alpha=0.55, zorder=2)
        axr.text(0.02, 0.995, "READS (shared multi-mapping evidence) -- ALL TIED: no distinguishing "
                 "bubble to thread through", fontsize=9, color="#444444", ha="left", va="top")
    else:
        # order: assigned (grouped by copy) first, then ambiguous, then tied; drop no_cover (partial
        # reads that reach no PSV column) from the drawn rows -- counted in the annotation instead.
        # 'recombinant' = the abstain leg's belongs-to-no-copy call; drawn with the abstainers.
        order = {"assigned": 0, "ambiguous": 1, "recombinant": 2, "tied": 3, "no_cover": 4}
        drawable = [r for r in reads if r["span_cols"] is not None and r["status"] != "no_cover"]
        drawable.sort(key=lambda r: (order[r["status"]],
                      names.index(r["best_copy"]) if r["best_copy"] else 99, r["span_cols"][0]))
        if len(drawable) > max_reads:                       # proportional, deterministic subsample
            by_status = collections.defaultdict(list)
            for r in drawable:
                by_status[r["status"]].append(r)
            picked = []
            for st in ("assigned", "ambiguous", "recombinant", "tied"):
                lst = by_status.get(st, [])
                if not lst:
                    continue
                q = max(1, round(max_reads * len(lst) / len(drawable)))
                picked.extend(lst[::max(1, len(lst) // max(q, 1))][:q])
            drawable = sorted(picked, key=lambda r: (order[r["status"]],
                              names.index(r["best_copy"]) if r["best_copy"] else 99,
                              r["span_cols"][0]))
        rows = max(len(drawable), 1)
        dy = 0.9 / rows
        for ri, r in enumerate(drawable):
            y = 0.95 - ri * dy
            c0, c1 = r["span_cols"]
            xs = [col2x(c) for c in shown if c0 <= c <= c1]
            if xs:
                x0, x1 = min(xs), max(xs)
            else:
                x0 = 0.06 + 0.90 * (c0 / max(len(bubbles) - 1, 1))
                x1 = 0.06 + 0.90 * (c1 / max(len(bubbles) - 1, 1))
            if x1 - x0 < 0.008:
                x1 = x0 + 0.008
            col_r = color[r["best_copy"]] if r["status"] == "assigned" else "#bdbdbd"
            lw = 1.8 if r["status"] == "assigned" else 1.0
            axr.plot([x0, x1], [y, y], color=col_r, lw=lw, solid_capstyle="round",
                     alpha=(0.95 if r["status"] == "assigned" else 0.5), zorder=2)
            for c in shown:                                  # tick where read crosses a shown bubble
                if c0 <= c <= c1:
                    axr.plot([col2x(c)], [y], marker="|", ms=3.0, color=col_r, zorder=3)
        axr.text(0.02, 0.995, "READS threaded through the graph (shared multi-mapping evidence) -- "
                 "colored by ASSIGNED copy, gray = unresolved", fontsize=9, color="#444444",
                 ha="left", va="top")

    # ---------------- title + annotation block ----------
    cp = vg["copies"]
    tier_str = " / ".join(f"T{c['tier']}:{c['gene']}({c['n_sun']}SUN)" for c in cp)
    asn = vg["assignment"]
    tot = max(vg["n_reads"], 1)
    pct = lambda k: 100 * asn.get(k, 0) / tot
    role = vg["role"]
    head = (f"O2 copy-assignment VARIATION GRAPH  |  family {vg['family']}  {vg['gene']}  "
            f"[{role}]  |  {vg['n_copies']} copies, "
            f"{'CROSS-CHROM ' if vg['cross_chrom'] else ''}{','.join(c[-3:] for c in vg['chroms'])}")
    sub = (f"{vg['n_bubbles']} PSV bubbles ({vg['n_sun_bubbles']} SUN)  |  K = {vg['K']} resolvable "
           f"haplotype-paths of {vg['n_copies']}  ({vg['cls']})  |  tiers "
           f"T1:{vg['tiers'].get(1,0)} T2:{vg['tiers'].get(2,0)} T3:{vg['tiers'].get(3,0)}")
    if not vg["bubbles"]:
        # floor: no bubble -> every read is TIED (the no_cover bookkeeping IS the tie here)
        n_tied_floor = asn.get("tied", 0) + asn.get("no_cover", 0)
        asg = (f"reads = {vg['n_reads']}   ASSIGNED 0 (0%)   "
               f"TIED (K-frontier floor) {n_tied_floor} ({100*n_tied_floor/tot:.0f}%)   "
               f"-- no distinguishing bubble exists, so NO read is resolvable, even in principle")
    else:
        # all five categories are listed so they SUM to n_reads (recombinant = abstain leg's
        # belongs-to-no-copy call; previously omitted, which made the footer appear not to total).
        asg = (f"reads = {vg['n_reads']}   ASSIGNED {asn.get('assigned',0)} ({pct('assigned'):.0f}%)   "
               f"AMBIGUOUS {asn.get('ambiguous',0)} ({pct('ambiguous'):.0f}%)   "
               f"RECOMBINANT/abstain {asn.get('recombinant',0)} ({pct('recombinant'):.0f}%)   "
               f"TIED {asn.get('tied',0)} ({pct('tied'):.0f}%)   "
               f"no-PSV-cover {asn.get('no_cover',0)} ({pct('no_cover'):.0f}%)")
    fig.suptitle(head, fontsize=13.5, fontweight="bold", y=0.985)
    fig.text(0.5, 0.925, sub, fontsize=10.5, ha="center", color="#333333")
    fig.text(0.5, 0.052, asg, fontsize=10.5, ha="center", color="#222222",
             bbox=dict(boxstyle="round,pad=0.4", fc="#eef4ff", ec="#9db8e0"))

    # copy legend (color = copy path; gene + tier + assigned reads)
    handles = []
    for c in cp:
        nm = c["name"]
        handles.append(Line2D([0], [0], color=color[nm], lw=3,
                       label=f"{nm} {c['gene']}  T{c['tier']}  ({c['n_assigned_reads']} reads)"))
    handles.append(Line2D([0], [0], color="#bdbdbd", lw=3, label="TIED / unresolved read (gray)"))
    fig.legend(handles=handles, loc="lower center", ncol=min(4, len(handles)),
               fontsize=8.5, frameon=False, bbox_to_anchor=(0.5, -0.005))

    fig.savefig(out_png, bbox_inches="tight", dpi=200)
    plt.close(fig)
    return out_png


# --------------------------------------------------------------------------- optional networkx view
def render_graph_view(vg, out_png, max_bubbles=24):
    try:
        import networkx as nx
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as ex:
        sys.stderr.write(f"[graph-view skipped: {ex}]\n")
        return None
    names = [c["name"] for c in vg["copies"]]
    color = {nm: COPY_COLORS[i % len(COPY_COLORS)] for i, nm in enumerate(names)}
    bubbles = vg["bubbles"]
    shown = _pick_bubbles(vg, max_bubbles)
    G = nx.DiGraph()
    pos = {}
    node_col, node_lbl = {}, {}
    G.add_node("S"); pos["S"] = (-1, 0); node_col["S"] = "#333333"; node_lbl["S"] = "src"
    G.add_node("T")
    for k, col in enumerate(shown):
        b = bubbles[col]
        for ai, (a, mem) in enumerate(sorted(b["alleles"].items())):
            nid = f"b{col}_{a}"
            G.add_node(nid)
            pos[nid] = (k, ALLELE_Y.get(a, 0.0) * 3 + (ai - len(b["alleles"]) / 2) * 0.15)
            priv = len(mem) == 1
            node_col[nid] = "#fff2b2" if priv else "#ffffff"
            node_lbl[nid] = a
    pos["T"] = (len(shown), 0); node_col["T"] = "#333333"; node_lbl["T"] = "snk"
    # copy paths
    for nm in names:
        prev = "S"
        for col in shown:
            a = bubbles[col]["hap"][nm]
            nid = f"b{col}_{a}"
            G.add_edge(prev, nid, copy=nm)
            prev = nid
        G.add_edge(prev, "T", copy=nm)
    fig, ax = plt.subplots(figsize=(16, 6), dpi=170)
    nx.draw_networkx_nodes(G, pos, node_size=260,
                           node_color=[node_col.get(n, "#ffffff") for n in G.nodes()],
                           edgecolors="#888888", ax=ax)
    nx.draw_networkx_labels(G, pos, labels={n: node_lbl.get(n, "") for n in G.nodes()},
                            font_size=7, ax=ax)
    for nm in names:
        edges = [(u, v) for u, v, d in G.edges(data=True) if d.get("copy") == nm]
        nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=color[nm], width=1.3,
                               arrows=False, ax=ax, alpha=0.7)
    ax.set_title(f"family {vg['family']} {vg['gene']} -- variation-graph bubble view "
                 f"(copy-paths thread allele nodes; yellow = SUN private allele)", fontsize=11)
    ax.axis("off")
    fig.savefig(out_png, bbox_inches="tight", dpi=170)
    plt.close(fig)
    return out_png


# --------------------------------------------------------------------------- driver
def run(fam_id, use_cache=True):
    vg = materialize_family(fam_id, use_cache=use_cache)
    gene = vg["gene"]
    base = os.path.join(HERE, f"fig_o2_vg_fam{fam_id}_{gene}")
    render_family(vg, base + ".png")
    render_graph_view(vg, base + "_graph.png")
    a = vg["assignment"]
    print(f"[fam {fam_id} {gene} {vg['role']}] copies={vg['n_copies']} bubbles={vg['n_bubbles']} "
          f"(SUN {vg['n_sun_bubbles']}) K={vg['K']} cls={vg['cls']}")
    print(f"   tiers T1/T2/T3 = {vg['tiers'].get(1,0)}/{vg['tiers'].get(2,0)}/{vg['tiers'].get(3,0)}"
          f"   reads={vg['n_reads']} assigned={a.get('assigned',0)} ambiguous={a.get('ambiguous',0)}"
          f" recombinant={a.get('recombinant',0)} tied={a.get('tied',0)} no_cover={a.get('no_cover',0)}")
    print(f"   per-copy assigned: {vg['per_copy']}")
    print(f"   -> {base}.png")
    return vg


if __name__ == "__main__":
    os.environ.setdefault("PYTHONHASHSEED", "0")
    args = sys.argv[1:]
    fresh = "--fresh" in args
    # DEFAULT-ON recombinant-abstain leg; --no-recombinant-abstain opts out (sets the env the leg reads),
    # mirroring --no-repeat-gate / --no-split-recombinants.
    if "--no-recombinant-abstain" in args:
        os.environ[ra.OPTOUT_ENV] = "1"
    nums = [int(a) for a in args if not a.startswith("--")]
    fam_ids = nums if nums else list(FLAGSHIP)
    for fid in fam_ids:
        run(fid, use_cache=not fresh)
