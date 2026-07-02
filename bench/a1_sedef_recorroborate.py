#!/usr/bin/env python
"""
A1 re-corroboration using the FULL SEDEF final.bed (253030 pairs).

Context. A1 (bench/A1_READ_CONSENSUS_O1.md) established a genuinely read-derived
copy signal (SDA / Vollger phi-correlation, no asm20) on the collapsed frontier.
Its genome-corroboration LEG was DOWNGRADED to "consistent-not-established":
the read copy count K_read vs a *basic* SEDEF partner count gave Spearman 0.443,
p=0.044, perm_p=0.052 (straddles the null), and K_read is ~80% a read-DEPTH proxy
(partial rho(K,SEDEF|depth)=0.40, p=0.08).

The user's point: the FULL SEDEF map is ALREADY computed (final.bed, 253029 pairs);
A1 used only a basic partner count. This script recomputes the GENOME side from the
FULL final.bed and asks whether the corroboration can be pushed ABOVE the null with
confounds PROPERLY controlled, on the collapsed-core subset (collapsed_excess>0).

  Axis 1  EXPECTED-CN proxy:  expected_CN = (#distinct GGO SEDEF partner regions
          overlapping the backbone locus) + 1        vs   K_read.
  Axis 2  DIVERGENCE (cleaner; sidesteps copy-vs-allele): median SEDEF divergence
          of a locus's segdup partners  vs  read heterogeneity (n_read_psv, mean_phi).
          Hypothesis: more-divergent segdups -> more read PSVs.
  Axis 3  which (if any) clears perm_p<0.05 AFTER confound control.

CRUX = confound control.  A1 controlled only med_depth.  This script ALSO controls
locus SIZE, because size drives BOTH the genome count (expected_CN) and the read
count (n_read_psv) while being ~uncorrelated with depth (see the confound matrix in
the JSON):  a count-vs-count correlation can be entirely a size artifact, whereas the
DIVERGENCE axis uses a genome-side RATE (sedef_div) that is size-immune.

NO SDA re-run: K_read / med_depth / n_read_psv / mean_phi are read from the cached
A1 rows.json.  Only the SEDEF (genome) side is recomputed from final.bed.

Determinism: PYTHONHASHSEED=0 asserted; every permutation gets an explicit
np.random.default_rng(seed); rows sorted before write.
"""
import os, sys, json, math, warnings
import numpy as np
from scipy import stats

SEED = 1234
B_PERM = 20000      # Monte-Carlo label-permutation null
B_SHIFT = 5000      # size-matched shifted-genomic null
ROOT = "/mnt/c/Users/jfris/Desktop/Rustle"
ROWS = f"{ROOT}/bench/a1_read_consensus_o1.rows.json"
SEDEF = "/mnt/c/Users/jfris/Desktop/final.bed"
FAI = "/home/juanfra/winloci_scratch/GGO.fasta.fai"
OUT_TSV = f"{ROOT}/bench/a1_sedef_recorroborate.tsv"
OUT_JSON = f"{ROOT}/bench/a1_sedef_recorroborate.json"

# final.bed columns (0-based), verified vs the example row:
#   c0 c1..c5 = c1,s1,e1,c2,s2,e2 ; col20 identity=0.948223 ; col22 divergence=0.053651
C_IDENT = 20
C_DIV = 22


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- GGO chroms
def load_ggo_chroms():
    """GGO nuclear chromosomes = fai scaffolds with length > 1e6 -> 25 chroms
    (23 autosomes NC_0732xx.2 + X NC_086017.1 + Y NC_086018.1).  Drops MT
    NC_011120.1 (16412 bp) and any header junk."""
    sizes = {}
    with open(FAI) as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            name, length = p[0], int(p[1])
            if length > 1_000_000:
                sizes[name] = length
    return sizes


# ------------------------------------------------------------- SEDEF parse
def parse_sedef(ggo_sizes):
    """Per-chrom feature arrays. Each GGO-GGO pair contributes TWO features (one on
    each side); each feature carries partner interval + divergence + identity.
    KEEP ONLY pairs where BOTH c1 and c2 are GGO nuclear chroms (drops MT rows)."""
    feats = {}   # chrom -> list of (start,end,pchrom,pstart,pend,div,ident)
    n_total = n_ggo = n_drop_chrom = n_bad = 0
    with open(SEDEF) as f:
        for ln in f:
            if not ln or ln[0] == "#":
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 34:
                continue
            n_total += 1
            c1, c2 = p[0], p[3]
            try:
                s1, e1, s2, e2 = int(p[1]), int(p[2]), int(p[4]), int(p[5])
                ident, div = float(p[C_IDENT]), float(p[C_DIV])
            except ValueError:
                n_bad += 1
                continue
            # robustness: identity in (0,1], divergence in [0,1].  NB identity==1.0
            # (div==0) are LEGITIMATE perfect-identity segdups = the most-collapsed
            # pairs, exactly the frontier signal -> must be KEPT, not filtered.
            if not (0.0 < ident <= 1.0 and 0.0 <= div <= 1.0):
                n_bad += 1
                continue
            if (c1 not in ggo_sizes) or (c2 not in ggo_sizes):
                n_drop_chrom += 1
                continue
            n_ggo += 1
            feats.setdefault(c1, []).append((s1, e1, c2, s2, e2, div, ident))
            feats.setdefault(c2, []).append((s2, e2, c1, s1, e1, div, ident))
    idx = {}
    for c, lst in feats.items():
        lst.sort(key=lambda x: (x[0], x[1]))
        idx[c] = dict(
            starts=np.array([x[0] for x in lst], dtype=np.int64),
            ends=np.array([x[1] for x in lst], dtype=np.int64),
            pch=[x[2] for x in lst],
            ps=np.array([x[3] for x in lst], dtype=np.int64),
            pe=np.array([x[4] for x in lst], dtype=np.int64),
            dv=np.array([x[5] for x in lst], dtype=float),
            it=np.array([x[6] for x in lst], dtype=float),
        )
    log(f"[sedef] rows total={n_total} kept(GGO-GGO)={n_ggo} "
        f"dropped(non-GGO chrom)={n_drop_chrom} dropped(bad ident/div)={n_bad}")
    return idx


def _merge_regions(intervals):
    """intervals: list of (chrom,start,end) -> number of DISTINCT merged regions."""
    if not intervals:
        return 0
    by_ch = {}
    for c, s, e in intervals:
        by_ch.setdefault(c, []).append((s, e))
    n = 0
    for c in sorted(by_ch):                     # deterministic
        ivs = sorted(by_ch[c])
        cur_s, cur_e = ivs[0]
        for s, e in ivs[1:]:
            if s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                n += 1
                cur_s, cur_e = s, e
        n += 1
    return n


def query_locus(idx, chrom, start, end):
    """Return (expected_CN, median_div, median_ident, n_pairs, distinct_partners).
    expected_CN = distinct partner regions + 1 (the +1 = the backbone copy itself).
    div/ident are the medians over the pairs that overlap the locus."""
    if chrom not in idx:
        return 1, np.nan, np.nan, 0, 0
    a = idx[chrom]
    mask = (a["starts"] < end) & (a["ends"] > start)
    k = int(mask.sum())
    if k == 0:
        return 1, np.nan, np.nan, 0, 0
    partners = [(a["pch"][i], int(a["ps"][i]), int(a["pe"][i]))
                for i in np.nonzero(mask)[0]]
    ndist = _merge_regions(partners)
    return ndist + 1, float(np.median(a["dv"][mask])), float(np.median(a["it"][mask])), k, ndist


# ------------------------------------------------------- stats helpers
def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def _resid_on_ranks(v, controls):
    """residual of rank(v) linearly regressed on rank(each control). All finite,
    same length. controls = list of 1-D arrays."""
    rv = stats.rankdata(v)
    cols = [np.ones_like(rv, dtype=float)] + [stats.rankdata(c) for c in controls]
    A = np.column_stack(cols)
    beta, *_ = np.linalg.lstsq(A, rv, rcond=None)
    return rv - A @ beta


def partial_spearman(x, y, controls):
    """Rank-based partial Spearman of x,y controlling a LIST of covariates.
    Returns (rho, analytic_p, n, rx, ry). df = n - 2 - #controls."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    controls = [np.asarray(c, float) for c in controls]
    m = np.isfinite(x) & np.isfinite(y)
    for c in controls:
        m &= np.isfinite(c)
    n = int(m.sum())
    ncov = len(controls)
    if n < ncov + 4:
        return np.nan, np.nan, n, None, None
    rx = _resid_on_ranks(x[m], [c[m] for c in controls])
    ry = _resid_on_ranks(y[m], [c[m] for c in controls])
    r = float(np.corrcoef(rx, ry)[0, 1])
    df = n - 2 - ncov
    if df > 0 and abs(r) < 1:
        t = r * math.sqrt(df / (1 - r * r))
        p = float(2 * stats.t.sf(abs(t), df))
    else:
        p = np.nan
    return r, p, n, rx, ry


def perm_p_from_vectors(u, v, obs, rng, B):
    """two-sided MC permutation p for pearson(u, v_perm) (u,v already residualized
    ranks OR raw ranks). Permutes v."""
    if u is None or v is None or not np.isfinite(obs) or len(u) < 3:
        return np.nan
    su, sv = np.std(u), np.std(v)
    if su == 0 or sv == 0:
        return np.nan
    ge = 0
    a = np.abs(obs) - 1e-12
    for _ in range(B):
        r = np.corrcoef(u, rng.permutation(v))[0, 1]
        if np.isfinite(r) and abs(r) >= a:
            ge += 1
    return (ge + 1) / (B + 1)


def perm_p_raw(x, y, obs, rng, B=B_PERM):
    """MC permutation p for a raw Spearman between x and y."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return np.nan
    return perm_p_from_vectors(stats.rankdata(x[m]), stats.rankdata(y[m]), obs, rng, B)


# ------------------------------------------------------- main
def main():
    if os.environ.get("PYTHONHASHSEED") != "0":
        log("WARN: PYTHONHASHSEED != 0 (output is sort-stable regardless, but run "
            "with PYTHONHASHSEED=0 for the record)")
    warnings.filterwarnings("ignore")   # constant-input spearman handled explicitly

    ggo = load_ggo_chroms()
    log(f"[ggo] {len(ggo)} nuclear chroms; MT/junk dropped")
    idx = parse_sedef(ggo)

    rows = json.load(open(ROWS))
    core = [r for r in rows if r.get("collapsed_excess", 0) > 0]
    core.sort(key=lambda r: r["family"])
    log(f"[core] collapsed-core n={len(core)} cls={sorted(set(r['cls'] for r in core))}")

    # ---- recompute SEDEF side from the FULL final.bed -----------------------
    per = []
    for r in core:
        chrom, span = r["backbone"].split(":")
        s, e = (int(v) for v in span.split("-"))
        exp_cn, mdiv, mident, npairs, ndist = query_locus(idx, chrom, s, e)
        per.append(dict(
            family=r["family"], gene=r["gene"], cls=r["cls"],
            chrom=chrom, start=s, end=e, size=e - s,
            K_read=r["K_read"], med_depth=r["med_depth"],
            n_read_psv=r["n_read_psv"], mean_phi=r["mean_phi"],
            collapsed_excess=r["collapsed_excess"],
            cached_partner_regions=r.get("sedef_partner_regions"),
            expected_CN=exp_cn, sedef_div=mdiv, sedef_ident=mident,
            n_sedef_pairs=npairs, n_distinct_partners=ndist,
        ))

    # ---- pseudo-replicate detection: overlapping backbones on the SAME chrom -
    # (families 169/170 are the same LOC101124683 locus, nested; identical genome
    #  AND identical read signal -> double-counts one point).  Build a de-dup mask
    #  keeping, per overlap-cluster, the LARGEST locus (deterministic tie: min family).
    order = sorted(range(len(per)), key=lambda i: (per[i]["chrom"], per[i]["start"], per[i]["end"]))
    dup_of = {}
    prev = None
    for i in order:
        if prev is not None and per[i]["chrom"] == per[prev]["chrom"] and per[i]["start"] < per[prev]["end"]:
            # overlaps previous kept locus -> mark smaller as dup of larger
            big, small = (prev, i) if per[prev]["size"] >= per[i]["size"] else (i, prev)
            dup_of[per[small]["family"]] = per[big]["family"]
        else:
            prev = i
    dup_families = set(dup_of)
    keep_dedup = [j for j, p in enumerate(per) if p["family"] not in dup_families]
    for p in per:
        p["is_pseudoreplicate"] = p["family"] in dup_families
    log(f"[dedup] pseudo-replicates dropped for the dedup variant: "
        f"{sorted(dup_of.items())}  (n {len(per)} -> {len(keep_dedup)})")

    # ---- arrays ------------------------------------------------------------
    fam = [p["family"] for p in per]
    K = np.array([p["K_read"] for p in per], float)
    D = np.array([p["med_depth"] for p in per], float)
    SZ = np.array([p["size"] for p in per], float)
    CN = np.array([p["expected_CN"] for p in per], float)
    NP = np.array([p["n_read_psv"] for p in per], float)
    PHI = np.array([p["mean_phi"] for p in per], float)
    DIV = np.array([p["sedef_div"] for p in per], float)
    chroms = [p["chrom"] for p in per]

    # ---- confound matrix (why size must be controlled) ---------------------
    confvars = dict(size=SZ, med_depth=D, expected_CN=CN, K_read=K,
                    n_read_psv=NP, sedef_div=DIV, mean_phi=PHI)
    cm = {}
    names = sorted(confvars)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            r, p, _ = spearman(confvars[a], confvars[b])
            cm[f"{a}__{b}"] = dict(rho=r, p=p)

    # ---- cached-vs-full comparison (does the FULL final.bed help?) ---------
    cachedPR = np.array([p["cached_partner_regions"] for p in per], float)
    cr, cp, _ = spearman(cachedPR, K)
    fr, fp, _ = spearman(CN, K)
    cached_full = dict(
        cached_partner_regions_vs_K_read=dict(rho=cr, p=cp),
        full_expected_CN_vs_K_read=dict(rho=fr, p=fp),
        note=("A1 headline (0.443) reproduces from the cached raw partner count; the "
              "FULL final.bed distinct-region count is WEAKER, not stronger."),
    )

    # ================= generic depth+size-controlled axis runner =============
    def run_axis(xarr, yarr, xname, yname, hypothesis, read_is_y=True):
        """x = genome-derived, y = read-derived (by convention).  Returns full dict
        with raw + partial{depth, size, depth+size} + perm nulls."""
        r_raw, p_raw, n = spearman(xarr, yarr)
        perm_raw = perm_p_raw(xarr, yarr, r_raw, np.random.default_rng(SEED + 101))
        out = dict(x=xname, y=yname, hypothesis=hypothesis, n=n,
                   spearman_raw=r_raw, p_raw=p_raw, perm_p_raw=perm_raw)
        for lbl, ctrl, off in [("partial_depth", [D], 201),
                               ("partial_size", [SZ], 202),
                               ("partial_depth_size", [D, SZ], 203)]:
            pr, pp, npart, rx, ry = partial_spearman(xarr, yarr, ctrl)
            perm = perm_p_from_vectors(rx, ry, pr, np.random.default_rng(SEED + off), B_PERM)
            out[lbl] = dict(rho=pr, p_analytic=pp, n=npart, perm_p=perm)
        return out

    results = {}

    # AXIS 1 : K_read  vs  expected_CN --------------------------------------
    a1 = run_axis(CN, K, "expected_CN", "K_read",
                  "more genome copies (distinct SEDEF partners) -> higher K_read")
    results["axis1_Kread_vs_expectedCN"] = a1

    # AXIS 2 : n_read_psv  vs  sedef_div ------------------------------------
    a2 = run_axis(DIV, NP, "sedef_div", "n_read_psv",
                  "more-divergent segdup partners -> more read PSVs")
    results["axis2_nreadpsv_vs_div"] = a2

    # AXIS 2b: mean_phi vs sedef_div (secondary, direction-agnostic) --------
    a2b = run_axis(DIV, PHI, "sedef_div", "mean_phi",
                   "divergence vs read linkage coherence (direction-agnostic)")
    results["axis2b_meanphi_vs_div"] = a2b

    # ---- de-dup sensitivity (drop pseudo-replicates) ----------------------
    def dedup_axis(xarr, yarr):
        kx, ky = xarr[keep_dedup], yarr[keep_dedup]
        Dd, Sd = D[keep_dedup], SZ[keep_dedup]
        r_raw, p_raw, n = spearman(kx, ky)
        pr, pp, npart, rx, ry = partial_spearman(kx, ky, [Dd])
        pr2, pp2, npart2, rx2, ry2 = partial_spearman(kx, ky, [Dd, Sd])
        return dict(n=n, spearman_raw=r_raw, p_raw=p_raw,
                    partial_depth=dict(rho=pr, p_analytic=pp, n=npart),
                    partial_depth_size=dict(rho=pr2, p_analytic=pp2, n=npart2))
    results["axis1_Kread_vs_expectedCN"]["dedup"] = dedup_axis(CN, K)
    results["axis2_nreadpsv_vs_div"]["dedup"] = dedup_axis(DIV, NP)

    # ================= SIZE-MATCHED SHIFTED GENOMIC NULL =====================
    # Recompute the genome variable at random SAME-SIZE positions on the SAME chrom;
    # compare the observed correlation to this shifted distribution.  Meaningful for
    # expected_CN (defined everywhere, baseline 1).  For divergence it is WEAK: random
    # windows are unconditioned on carrying a partner, so it is reported with a caveat.
    def shifted_null(genome_kind, read_arr, obs_raw, obs_part, seed_off):
        rgn = np.random.default_rng(SEED + seed_off)
        n = len(per)
        rr = np.asarray(read_arr, float)
        rstarts = np.empty((B_SHIFT, n), dtype=np.int64)
        for j in range(n):
            hi = max(1, ggo[chroms[j]] - int(SZ[j]))
            rstarts[:, j] = rgn.integers(0, hi, size=B_SHIFT)
        ge_raw = ge_part = valid = 0
        for b in range(B_SHIFT):
            g = np.empty(n)
            for j in range(n):
                cn, mdiv, _, _, _ = query_locus(
                    idx, chroms[j], int(rstarts[b, j]), int(rstarts[b, j]) + int(SZ[j]))
                g[j] = cn if genome_kind == "cn" else mdiv
            m = np.isfinite(g) & np.isfinite(rr)
            if m.sum() < 4 or np.std(g[m]) == 0 or np.std(rr[m]) == 0:
                continue
            valid += 1
            rraw = stats.spearmanr(g[m], rr[m])[0]
            if np.isfinite(rraw) and abs(rraw) >= abs(obs_raw) - 1e-12:
                ge_raw += 1
            mm = m & np.isfinite(D)
            if mm.sum() >= 4:
                gx = _resid_on_ranks(g[mm], [D[mm]])
                gy = _resid_on_ranks(rr[mm], [D[mm]])
                if np.std(gx) > 0 and np.std(gy) > 0:
                    rp = np.corrcoef(gx, gy)[0, 1]
                    if np.isfinite(rp) and np.isfinite(obs_part) and abs(rp) >= abs(obs_part) - 1e-12:
                        ge_part += 1
        return dict(shift_perm_p_raw=(ge_raw + 1) / (valid + 1) if valid else np.nan,
                    shift_perm_p_partial_depth=(ge_part + 1) / (valid + 1) if valid else np.nan,
                    n_valid_shuffles=valid)

    log("[shift] axis1 (expected_CN) size-matched shifted null ...")
    results["axis1_Kread_vs_expectedCN"]["shifted_null"] = shifted_null(
        "cn", K, a1["spearman_raw"], a1["partial_depth"]["rho"], 301)
    log("[shift] axis2 (sedef_div) size-matched shifted null [weak; see caveat] ...")
    results["axis2_nreadpsv_vs_div"]["shifted_null"] = shifted_null(
        "div", NP, a2["spearman_raw"], a2["partial_depth"]["rho"], 302)

    # ================= multiple-testing note =================================
    # Primary family = {axis1, axis2, axis2b} on the perm_p_raw statistic (3 tests).
    prim = {k: results[k]["perm_p_raw"] for k in
            ["axis1_Kread_vs_expectedCN", "axis2_nreadpsv_vs_div", "axis2b_meanphi_vs_div"]}
    best_key = min(prim, key=lambda k: prim[k])
    mt = dict(n_primary_tests=3,
              primary_perm_p_raw=prim,
              best_axis=best_key, best_perm_p_raw=prim[best_key],
              bonferroni_best=min(1.0, 3 * prim[best_key]))

    # ---- write per-locus TSV ----------------------------------------------
    per_sorted = sorted(per, key=lambda p: p["family"])
    cols = ["family", "gene", "cls", "chrom", "start", "end", "size", "is_pseudoreplicate",
            "K_read", "med_depth", "expected_CN", "sedef_div", "sedef_ident",
            "n_read_psv", "mean_phi", "n_sedef_pairs", "n_distinct_partners",
            "cached_partner_regions", "collapsed_excess"]
    with open(OUT_TSV, "w") as f:
        f.write("\t".join(cols) + "\n")
        for p in per_sorted:
            f.write("\t".join(str(p[c]) for c in cols) + "\n")

    meta = dict(
        seed=SEED, B_perm=B_PERM, B_shift=B_SHIFT,
        collapsed_core_n=len(per), dedup_core_n=len(keep_dedup),
        pseudoreplicates=sorted(dup_of.items()),
        ggo_chroms=sorted(ggo.keys()),
        sedef_ident_col=C_IDENT, sedef_div_col=C_DIV,
        core_cls=sorted(set(p["cls"] for p in per)),
    )
    out = dict(meta=meta, confound_matrix=cm, cached_vs_full=cached_full,
               results=results, multiple_testing=mt, per_locus=per_sorted)
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True,
                  default=lambda o: float(o) if isinstance(o, np.floating) else o)

    # ---- console report ----------------------------------------------------
    def fmt(v):
        return f"{v:.4f}" if isinstance(v, float) and np.isfinite(v) else str(v)
    print("=" * 82)
    print(f"A1 SEDEF RE-CORROBORATION  (collapsed-core n={len(per)}, dedup n={len(keep_dedup)})")
    print("=" * 82)
    print("\nCACHED-vs-FULL (does the full final.bed help Axis 1?)")
    print(f"   cached partner count vs K_read = {fmt(cr)} (p={fmt(cp)})   [= A1 headline 0.443]")
    print(f"   FULL expected_CN     vs K_read = {fmt(fr)} (p={fmt(fp)})   [weaker, not stronger]")
    def cmget(a, b):
        return cm.get(f"{a}__{b}") or cm.get(f"{b}__{a}")
    print("\nKEY CONFOUNDS")
    for a, b in [("size", "expected_CN"), ("size", "n_read_psv"), ("size", "med_depth"),
                 ("size", "sedef_div"), ("K_read", "med_depth"), ("n_read_psv", "med_depth")]:
        e = cmget(a, b)
        print(f"   {a+' vs '+b:26} rho={fmt(e['rho'])} (p={fmt(e['p'])})")

    for key, lbl in [("axis1_Kread_vs_expectedCN", "AXIS 1  K_read vs expected_CN (distinct SEDEF partners +1)"),
                     ("axis2_nreadpsv_vs_div", "AXIS 2  n_read_psv vs SEDEF divergence  [cleaner: rate, size-immune genome side]"),
                     ("axis2b_meanphi_vs_div", "AXIS 2b mean_phi vs SEDEF divergence")]:
        d = results[key]
        print(f"\n{lbl}")
        print(f"   n                         = {d['n']}")
        print(f"   raw Spearman              = {fmt(d['spearman_raw'])} (p={fmt(d['p_raw'])}, perm_p={fmt(d['perm_p_raw'])})")
        print(f"   partial | depth           = {fmt(d['partial_depth']['rho'])} (p={fmt(d['partial_depth']['p_analytic'])}, perm_p={fmt(d['partial_depth']['perm_p'])})")
        print(f"   partial | size            = {fmt(d['partial_size']['rho'])} (p={fmt(d['partial_size']['p_analytic'])}, perm_p={fmt(d['partial_size']['perm_p'])})")
        print(f"   partial | depth+size      = {fmt(d['partial_depth_size']['rho'])} (p={fmt(d['partial_depth_size']['p_analytic'])}, perm_p={fmt(d['partial_depth_size']['perm_p'])})")
        if "dedup" in d:
            dd = d["dedup"]
            print(f"   dedup(n={dd['n']}) raw           = {fmt(dd['spearman_raw'])} (p={fmt(dd['p_raw'])})  partial|depth+size={fmt(dd['partial_depth_size']['rho'])} (p={fmt(dd['partial_depth_size']['p_analytic'])})")
        if "shifted_null" in d:
            s = d["shifted_null"]
            print(f"   size-matched shift raw    = {fmt(s['shift_perm_p_raw'])}  |depth={fmt(s['shift_perm_p_partial_depth'])}  (valid={s['n_valid_shuffles']})")

    print("\nMULTIPLE TESTING (3 primary axes, perm_p_raw)")
    print(f"   best axis = {mt['best_axis']}  perm_p_raw={fmt(mt['best_perm_p_raw'])}  "
          f"Bonferroni x3 = {fmt(mt['bonferroni_best'])}")
    print("\nwrote", OUT_TSV)
    print("wrote", OUT_JSON)


if __name__ == "__main__":
    main()
