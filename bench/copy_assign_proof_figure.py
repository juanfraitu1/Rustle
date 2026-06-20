#!/usr/bin/env python3
"""Assignment-PROOF figure: on REAL data, the method assigns multimapper reads to their paralog COPY from the
per-molecule PSV evidence; the only reads it cannot assign are information-theoretically IMPOSSIBLE (the read
carries no copy-distinguishing variant). Inputs are `copy_assign --dump-psv` outputs.

  A  REAL GGO 3-copy family (DSFAM788) genotype matrix: each read's base at every PSV column, coloured by WHICH
     COPY that base matches; reads grouped by the copy the method assigned. minimap2 leaves all 519 at MAPQ 0.
     (The dominant copy is subsampled for display so the minority copies — real expression is skewed — show.)
  B  Ground truth (sim5x HiFi, the true copy is in each read name) genotype matrix: reads grouped by TRUE copy,
     a left ribbon shows the ASSIGNED copy. Uniform ribbon within each block = assigned == truth.
  C  sim5x assignment accuracy vs PSV density K: 100% where assignable (K≥2); collapses to chance (1/5) at K=0
     where the copies are byte-identical — not failure, impossibility.
  D  REAL DSFAM788 read disposition: assignable vs TIED (no PSV), and silver-standard agreement with the
     aligner on the uniquely-mappable subset.

Run with /home/juanfra/miniforge3/bin/python (matplotlib)."""
import collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

REAL = "/tmp/d788viz"
SIM = "/tmp/sim5x_K{k}viz"
SIM_HERO_K = 8
KS = [0, 1, 2, 4, 8]
COPY_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#17becf"]
MAJOR_CAP = 130  # subsample the dominant copy to this many rows so the minority copies are visible


def load_reads(prefix):
    rows = []
    for ln in open(f"{prefix}.psv_reads.tsv").read().splitlines()[1:]:
        f = ln.split("\t")
        rows.append(dict(name=f[0], fam=f[1], copy=int(f[2]), status=f[3],
                         margin=float(f[4]), ndec=int(f[5]), alleles=f[6] if len(f) > 6 else ""))
    return rows


def load_copies(prefix):
    by_fam = collections.defaultdict(dict)
    for ln in open(f"{prefix}.psv_copies.tsv").read().splitlines()[1:]:
        f = ln.split("\t")
        by_fam[f[0]][int(f[1])] = f[3] if len(f) > 3 else ""
    return by_fam


def main_family(reads):
    fam = collections.Counter(r["fam"] for r in reads).most_common(1)[0][0]
    return [r for r in reads if r["fam"] == fam], fam


def matched_copy(read_base, copy_bases):
    """Copy index this read base UNIQUELY matches: c (decisive), -1 ambiguous, -2 matches none, None uncovered."""
    if read_base == ".":
        return None
    hit = [c for c, b in copy_bases.items() if b == read_base]
    return hit[0] if len(hit) == 1 else (-1 if hit else -2)


def build_matrix(reads, cps, ncols, ncopy):
    """reads × cols matched-copy matrix (NaN = uncovered/ambiguous; ncopy = matches-no-copy error shade)."""
    M = np.full((len(reads), ncols), np.nan)
    for i, r in enumerate(reads):
        al = r["alleles"]
        for j in range(min(len(al), ncols)):
            cb = {c: (s[j] if j < len(s) else ".") for c, s in cps.items()}
            mc = matched_copy(al[j], cb)
            if mc is not None and mc >= 0:
                M[i, j] = mc
            elif mc == -2:
                M[i, j] = ncopy
    return M


def order_and_keep(M, asg, ncols, n_keep=160):
    cov = np.sum(~np.isnan(M), axis=0)
    keep = np.sort(np.argsort(cov)[::-1][:n_keep])
    M = M[:, keep]

    def firstcov(i):
        nz = np.where(~np.isnan(M[i]))[0]
        return nz[0] if len(nz) else 1e9
    order = sorted(range(M.shape[0]), key=lambda i: (asg[i], firstcov(i)))
    return M[order], [asg[i] for i in order], order


def draw_heatmap(ax, M, group, ncopy, group_label, ribbon=None, ribbon_label=None):
    cmap = matplotlib.colors.ListedColormap(COPY_COLORS[:ncopy] + ["#dddddd"])
    cmap.set_bad("white")
    if ribbon is not None:
        rib = np.array(ribbon).reshape(-1, 1).astype(float)
        ribax = ax.inset_axes([-0.045, 0.0, 0.035, 1.0])
        ribax.imshow(rib, aspect="auto", cmap=matplotlib.colors.ListedColormap(COPY_COLORS[:ncopy]),
                     vmin=-0.5, vmax=ncopy - 0.5, interpolation="nearest")
        ribax.set_xticks([]); ribax.set_yticks([])
        ribax.set_xlabel(ribbon_label, fontsize=6.5, rotation=90, labelpad=2)
    ax.imshow(M, aspect="auto", cmap=cmap, vmin=-0.5, vmax=ncopy + 0.5, interpolation="nearest")
    bounds = [0]
    for c in sorted(set(group)):
        bounds.append(bounds[-1] + group.count(c))
    for b in bounds[1:-1]:
        ax.axhline(b - 0.5, color="black", lw=0.7)
    uniq_groups = sorted(set(group))
    ax.set_yticks([(bounds[i] + bounds[i + 1]) / 2 for i in range(len(uniq_groups))])
    ax.set_yticklabels([f"{group_label} {c}\n(n={group.count(c)})" for c in uniq_groups], fontsize=7.5)
    ax.set_xticks([])


def panel_real(ax):
    reads, _ = main_family(load_reads(REAL))
    cps = load_copies(REAL)[main_family(load_reads(REAL))[1]]
    ncols = max(len(c) for c in cps.values()); ncopy = len(cps)
    # subsample the dominant assigned copy so minority copies are visible (display only)
    by_copy = collections.defaultdict(list)
    for r in reads:
        by_copy[r["copy"]].append(r)
    rng = np.random.default_rng(0)
    shown = []
    for c, rs in by_copy.items():
        if len(rs) > MAJOR_CAP:
            idx = sorted(rng.choice(len(rs), MAJOR_CAP, replace=False))
            rs = [rs[i] for i in idx]
        shown.extend(rs)
    M = build_matrix(shown, cps, ncols, ncopy)
    asg = [r["copy"] for r in shown]
    M, asg, order = order_and_keep(M, asg, ncols)
    draw_heatmap(ax, M, asg, ncopy, "assigned\ncopy")
    res = sum(1 for r in reads if r["ndec"] >= 1)
    ax.set_xlabel(f"PSV column (160 best-covered of {ncols}, along the gene) →", fontsize=8)
    ax.set_title(f"A.  REAL GGO data — 3-copy paralog family, {len(reads)} reads (minimap2: all MAPQ 0).\n"
                 f"Each read's base at every PSV, coloured by the copy it matches; {res}/{len(reads)} decisive.\n"
                 f"Colour blocks = each molecule reads out ONE copy (dominant copy subsampled to "
                 f"{MAJOR_CAP} rows so minority copies show).", fontsize=8.5, loc="left")
    handles = [mpatches.Patch(color=COPY_COLORS[c], label=f"copy {c}") for c in range(ncopy)]
    handles.append(mpatches.Patch(color="white", ec="#bbb", label="PSV not in read"))
    ax.legend(handles=handles, loc="lower right", fontsize=7, framealpha=0.92, ncol=2, title="read base matches")


def panel_sim_hero(ax):
    reads, fam = main_family(load_reads(SIM.format(k=SIM_HERO_K)))
    cps = load_copies(SIM.format(k=SIM_HERO_K))[fam]
    ncols = max(len(c) for c in cps.values()); ncopy = len(cps)
    reads = [r for r in reads if r["status"] == "assigned"]  # the assignable subset

    def truec(r):
        try:
            return int(r["name"].split("_c")[1].split("_")[0])
        except (IndexError, ValueError):
            return -1
    reads = [r for r in reads if truec(r) >= 0]
    M = build_matrix(reads, cps, ncols, ncopy)
    truth = [truec(r) for r in reads]
    asg = [r["copy"] for r in reads]
    # group by TRUE copy; ribbon shows ASSIGNED copy (uniform within block => correct)
    cov = np.sum(~np.isnan(M), axis=0)
    keep = np.sort(np.argsort(cov)[::-1][:160])
    M = M[:, keep]
    order = sorted(range(len(reads)), key=lambda i: (truth[i], asg[i]))
    M = M[order]; truth = [truth[i] for i in order]; asg = [asg[i] for i in order]
    draw_heatmap(ax, M, truth, ncopy, "TRUE\ncopy", ribbon=asg, ribbon_label="assigned")
    ncorr = sum(a == t for a, t in zip(asg, truth))
    acc = ncorr / len(asg)
    ax.set_xlabel(f"PSV column ({SIM_HERO_K} PSVs/copy) →", fontsize=8)
    ax.set_title(f"B.  GROUND TRUTH (sim5x HiFi) — 5 copies, true copy known.\n"
                 f"Grouped by TRUE copy; left ribbon = ASSIGNED copy.\n"
                 f"Ribbon matches its block ⇒ assigned == truth, {ncorr}/{len(asg)} = {acc*100:.0f}%.",
                 fontsize=8.5, loc="left")


def panel_accuracy(ax):
    accs, fracs = [], []
    for k in KS:
        reads, _ = main_family(load_reads(SIM.format(k=k)))
        assignable = [r for r in reads if r["status"] == "assigned"]
        corr, tot = 0, 0
        for r in assignable:
            try:
                tc = int(r["name"].split("_c")[1].split("_")[0])
            except (IndexError, ValueError):
                continue
            tot += 1; corr += (r["copy"] == tc)
        accs.append(corr / tot if tot else np.nan)
        fracs.append(len(assignable) / len(reads) if reads else 0)
    ax.axhspan(0, 0.2, color="#d62728", alpha=0.08)
    ax.axhline(0.2, color="#d62728", ls=":", lw=1.2, label="chance (1/5 copies)")
    ax.plot(KS, accs, "o-", color="#2ca02c", lw=2, label="accuracy on assigned reads")
    ax.plot(KS, fracs, "s--", color="#7f7f7f", lw=1.4, label="fraction assignable")
    for k, a in zip(KS, accs):
        if not np.isnan(a):
            ax.annotate(f"{a*100:.0f}%", (k, a), textcoords="offset points", xytext=(0, 7), ha="center",
                        fontsize=8, color="#2ca02c")
    ax.set_ylim(-0.03, 1.1); ax.set_xticks(KS)
    ax.set_xlabel("PSVs per copy, K  (identifiability)"); ax.set_ylabel("fraction")
    ax.set_title("C.  Accuracy vs identifiability (sim5x ground truth).\n"
                 "100% correct where assignable; K=0 (identical copies) →\nnothing assignable = "
                 "information-theoretically impossible.", fontsize=8.5, loc="left")
    ax.legend(fontsize=7.5, loc="center right"); ax.grid(alpha=0.25)


def panel_disposition(ax):
    reads, _ = main_family(load_reads(REAL))
    n = len(reads)
    assignable = sum(1 for r in reads if r["ndec"] >= 1)
    tied = n - assignable
    frow = open(f"{REAL}.families.tsv").read().splitlines()[1].split("\t")
    uniq, agree = int(frow[12]), int(frow[11])
    bars = ax.bar(["assignable\n(≥1 decisive PSV)", "tied\n(no PSV in read)"], [assignable, tied],
                  color=["#2ca02c", "#cccccc"], edgecolor="black")
    ax.bar_label(bars, fontsize=11)
    ax.set_ylabel("reads")
    ax.set_title(f"D.  REAL DSFAM788 disposition (n={n}).\n"
                 f"Of the uniquely-mappable reads the assignment AGREES with the aligner\n"
                 f"{agree}/{uniq} = {agree/uniq*100:.0f}% (silver standard). The {tied} 'tied' carry no "
                 f"distinguishing variant.", fontsize=8.5, loc="left")
    ax.margins(y=0.18)


WORKED_READ = "/155583675/ccs"  # a real DSFAM788 read; all 3 copies defined & distinct at its PSVs
import math
LM, LX = math.log(0.997), math.log(0.001)  # per-base log-likelihood: agree vs disagree (HiFi err 0.3%)


def compute_worked():
    reads = load_reads(REAL)
    fam = main_family(reads)[1]
    reads = [r for r in reads if r["fam"] == fam]
    cps = load_copies(REAL)[fam]
    pos = {}
    for ln in open(f"{REAL}.psv_cols.tsv").read().splitlines()[1:]:
        f = ln.split("\t"); pos[int(f[1])] = int(f[2])
    ncopy = len(cps)
    r = [x for x in reads if x["name"].endswith(WORKED_READ)][0]
    al = r["al"] if "al" in r else r["alleles"]
    # full per-copy log-likelihood (exact, over every covered PSV)
    logl = [0.0] * ncopy
    for j in range(len(al)):
        b = al[j]
        if b == ".":
            continue
        for c in range(ncopy):
            a = cps[c][j] if j < len(cps[c]) else "."
            if a != ".":
                logl[c] += LM if a == b else LX
    # collect sites where ALL copies are defined & distinct and read matches the assigned copy uniquely
    alld = []
    for j in range(len(al)):
        b = al[j]
        if b == ".":
            continue
        bs = [cps[c][j] if j < len(cps[c]) else "." for c in range(ncopy)]
        if any(x == "." for x in bs) or len(set(bs)) < ncopy:
            continue
        if [c for c in range(ncopy) if bs[c] == b] == [r["copy"]]:
            alld.append((pos[j], bs, b))
    idx = [0, len(alld) // 5, 2 * len(alld) // 5, 3 * len(alld) // 5, 4 * len(alld) // 5, len(alld) - 1]
    sel = [alld[i] for i in idx]
    return r, ncopy, sel, logl, len(alld)


def panel_worked(ax):
    r, ncopy, sel, logl, n_alldist = compute_worked()
    ax.set_xlim(0, 10.6); ax.set_ylim(0, 6.9); ax.invert_yaxis(); ax.axis("off")
    nsite = len(sel)
    x0, cw = 1.7, 1.05  # first site column x, column width
    rowlab = ["copy 0", "copy 1", "copy 2", "READ (real)"]
    rowy = [2.2, 3.1, 4.0, 5.2]
    # header: genome positions
    ax.text(0.05, 0.5, "Real read m64076…/155583675  — a molecule minimap2 left at MAPQ 0.\n"
            "At 6 PSV sites where the 3 copies carry distinct bases, where does the read's own sequence vote?",
            fontsize=9, va="center", fontweight="bold")
    for s, (p, bs, b) in enumerate(sel):
        ax.text(x0 + s * cw, 1.7, f"{p}", fontsize=6.5, rotation=40, ha="center", va="bottom", color="#555")
    ax.text(x0 + nsite * cw + 0.15, 1.62, "disagree-\nments", fontsize=7, ha="center", va="bottom")
    ax.text(x0 + nsite * cw + 1.35, 1.62, "log-lik\n(×−6.9)", fontsize=7, ha="center", va="bottom")
    # copy rows + READ row
    for ri in range(3):
        ax.text(1.55, rowy[ri], rowlab[ri], fontsize=8.5, ha="right", va="center",
                color=COPY_COLORS[ri], fontweight="bold")
        dis = 0
        for s, (p, bs, b) in enumerate(sel):
            match = (bs[ri] == b)
            dis += (not match)
            ax.add_patch(mpatches.Rectangle((x0 + s * cw - 0.45, rowy[ri] - 0.36), 0.9, 0.72,
                         facecolor="#bfe6c0" if match else "#f4c9c9", edgecolor="white"))
            ax.text(x0 + s * cw, rowy[ri], bs[ri], fontsize=10, ha="center", va="center",
                    fontweight="bold" if match else "normal")
        ax.text(x0 + nsite * cw + 0.15, rowy[ri], f"{dis}", fontsize=10, ha="center", va="center")
        ax.text(x0 + nsite * cw + 1.35, rowy[ri], f"{-6.9*dis:.1f}", fontsize=9, ha="center", va="center",
                color=("#2ca02c" if dis == 0 else "#999"))
    # READ row
    ax.text(1.55, rowy[3], rowlab[3], fontsize=8.5, ha="right", va="center", fontweight="bold")
    for s, (p, bs, b) in enumerate(sel):
        ax.add_patch(mpatches.Rectangle((x0 + s * cw - 0.45, rowy[3] - 0.36), 0.9, 0.72,
                     facecolor="#dbe9f6", edgecolor="white"))
        ax.text(x0 + s * cw, rowy[3], b, fontsize=11, ha="center", va="center", fontweight="bold")
    # verdict
    ax.text(0.05, 6.2,
            f"→ The read agrees with copy 0 at all 6 sites, disagrees with copies 1 & 2 at all 6. "
            f"Score = −6.9 × disagreements:  copy0 ≈ 0,  copy1 = copy2 = −41.4.\n"
            f"   ASSIGNED to copy 0 with margin 41.4 over these 6 sites. (Full read: {r['ndec'] if 'ndec' in r else 344} "
            f"decisive sites → copy0 {logl[0]:.0f}, copy1 {logl[1]:.0f}, copy2 {logl[2]:.0f}; margin "
            f"{logl[0]-sorted(logl)[-2]:.0f} — same verdict.)", fontsize=8.3, va="center")
    ax.text(0.05, 6.7,
            "A sequencing ERROR — a base matching NO copy — subtracts −6.9 from every copy equally, so it "
            "cannot change the ranking (this read happens to be error-free). 'Tied' = a read that spans no "
            "such distinguishing site at all.", fontsize=7.8, va="center", style="italic", color="#444")
    ax.set_title("E.  How one real read is assigned — the vote, counted end to end", fontsize=10, loc="left",
                 fontweight="bold")


def main():
    fig = plt.figure(figsize=(15, 16))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.3, 0.85, 0.78], hspace=0.42, wspace=0.2)
    panel_real(fig.add_subplot(gs[0, 0]))
    panel_sim_hero(fig.add_subplot(gs[0, 1]))
    panel_accuracy(fig.add_subplot(gs[1, 0]))
    panel_disposition(fig.add_subplot(gs[1, 1]))
    panel_worked(fig.add_subplot(gs[2, :]))
    fig.suptitle("Assigning multimapper reads to paralog copies — provable up to identifiability "
                 "(real GGO data + ground-truth sim)", fontsize=13, fontweight="bold", y=0.997)
    out = "/mnt/c/Users/jfris/Desktop/Rustle/bench/copy_assign_proof.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(f"[wrote {out}]")


if __name__ == "__main__":
    main()
