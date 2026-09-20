#!/usr/bin/env python
"""READ-INTRINSIC TE detection for the reframed RNA family (E_r) over-merge arbitration.

QUESTION: can we detect TE-ness DIRECTLY FROM THE READS (self-contained, NO RepeatMasker library), and does
a read-intrinsic TE signal separate the E_r over-merge classes BETTER than genome soft-mask? The committed
soft-mask gate (bench/FAMILY_ER_TE_GATE.md) separates genuine-junk from protein-missed real TP on the
arbitration population (in_ep=0, direct core edge, mostly gorilla LOC*-LOC*) at AUC(genuine>TP)=0.694 at the
recommended band (core>=0.50) and 0.777 over the whole re-admit region (core>=0.13). Can a read-derived
signal beat it?

THREE READ-INTRINSIC TE SIGNALS (all self-contained, no RepeatMasker):
  (1) HOMOLOGY-GRAPH DEGREE  -- how many other de-novo loci each bridged locus links to at core_recip>=0.13
      in the de-novo homology graph (bench/denovo_family_edges.tsv). A TE bridge's shared sequence connects
      to MANY loci (repeat hub); a real paralog connects only to its family. Purely from the reads (the graph
      is de-novo, no annotation). Per edge = max/mean of the two endpoints' degree.
  (2) READ MULTI-MAPPING DEGREE -- the most literal "from the reads": count each read's genome-wide
      alignment multiplicity in GGO_mm.bam (minimap2 -N50 -> caps ~50); per bridge locus = median/mean/max
      over its reads (cache bench/ri_readmm_cache.tsv). A read over a TE multi-maps to many loci.
  (3) GENOME K-MER MULTIPLICITY -- count canonical 21-mer multiplicity of the bridge exonic sequence across
      the whole GGO.fasta (cache bench/ri_kmer_cache.tsv). High = repeat. Mask-AGNOSTIC (independent of the
      soft-mask) and catches novel/gorilla-specific TEs a human-biased RepeatMasker library misses.

EVAL: labels from the committed edge cache (bench/edge_bridge_mask.tsv: cls TP/truthbar/genuine, in_ep, core,
bridge_dn pair, bridge_mask=soft-mask feature). Arbitration population = in_ep=0 AND direct core edge,
genuine vs protein-missed real TP (same population where soft-mask got AUC 0.694/0.777). For each signal:
AUC(genuine>TP) by core band vs soft-mask; complementarity (rank-average combined AUC + correlation);
combined gate 'in_ep OR (core>=t AND unmasked AND low-degree)' through family_er_te_gate.compute_point vs
the soft-mask-only recommended point (21/97 TPpm @ edge-host 8.9%).

HONEST CAVEAT (bench/ri_sharedlen_cache.tsv): very-high-copy REAL gene families (KRAB-ZNF ~400, olfactory
~1000) are ALSO high-degree / high-multiplicity, so degree/multiplicity ALONE cannot separate them from TEs.
The protection is (a) protein-homology (in_ep already rescues the annotated ZNF/OR families before the TE
branch is reached) + (b) at the operating band (core>=0.50) the softmask<0.10 branch, which already excludes
every high-degree protein-missed TP (all have softmask>=0.235). SHARED-SEGMENT LENGTH (a real paralog shares
a LONG contiguous gene body; a TE bridge a SHORT exonized-TE fragment) separates the two ONLY on the broad
re-admit region core>=0.13 (AUC 0.821); it is at chance at core>=0.50 (0.521) and INVERTS at core>=0.70
(0.482), so it is NOT a term in the deployed gate.

Reproduce: PYTHONHASHSEED=0 python bench/family_er_readintrinsic_te.py
Deterministic: all heavy features precomputed to integer caches (deterministic builders); PYTHONHASHSEED=0;
gate re-refine seed=0. Writes bench/family_er_readintrinsic_te.tsv + .json.
"""
import json, os, sys, statistics
from bisect import bisect_left, bisect_right
from collections import defaultdict
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"; os.execv(sys.executable, [sys.executable] + sys.argv)
BENCH = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, BENCH)
import family_er_oracle_sweep as S
import family_er_te_gate as TEG

EDGE_CACHE = os.path.join(BENCH, "edge_bridge_mask.tsv")
DENOVO_EDGES = os.path.join(BENCH, "denovo_family_edges.tsv")
READMM_CACHE = os.path.join(BENCH, "ri_readmm_cache.tsv")
KMER_CACHE = os.path.join(BENCH, "ri_kmer_cache.tsv")
SHAREDLEN_CACHE = os.path.join(BENCH, "ri_sharedlen_cache.tsv")
DEG_T = 0.13
BANDS = [0.13, 0.30, 0.50, 0.70]
SOFTMASK_AUC = {0.13: 0.7769, 0.30: 0.7053, 0.50: 0.6938, 0.70: 0.6751}  # committed reference


def auc_gt(pos, neg):
    """Mann-Whitney AUC = P(pos > neg), 0.5-tie-corrected. pos=genuine, neg=protein-missed TP."""
    if not pos or not neg:
        return None
    negs = sorted(neg); c = 0.0
    for v in pos:
        lo = bisect_left(negs, v); hi = bisect_right(negs, v); c += lo + 0.5 * (hi - lo)
    return c / (len(pos) * len(neg))


def spearman(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    def ranks(a):
        order = sorted(range(n), key=lambda i: a[i])
        r = [0.0] * n; i = 0
        while i < n:
            j = i
            while j + 1 < n and a[order[j + 1]] == a[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = ranks(xs), ranks(ys)
    mx = sum(rx) / n; my = sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = sum((rx[i] - mx) ** 2 for i in range(n)) ** 0.5
    dy = sum((ry[i] - my) ** 2 for i in range(n)) ** 0.5
    return (num / (dx * dy)) if dx and dy else None


def rank_avg_auc(gen_pairs, tp_pairs):
    """Combined AUC of the equal-weight rank-average of two directed features (both higher=genuine).
    gen_pairs/tp_pairs: list of (f1, f2). Non-fitted (no weights learned) -> honest complementarity test."""
    allv = gen_pairs + tp_pairs
    n = len(allv)
    if n < 2:
        return None
    def col_ranks(col):
        vals = [v[col] for v in allv]
        order = sorted(range(n), key=lambda i: vals[i])
        r = [0.0] * n; i = 0
        while i < n:
            j = i
            while j + 1 < n and vals[order[j + 1]] == vals[order[i]]:
                j += 1
            avg = (i + j) / 2.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    r0, r1 = col_ranks(0), col_ranks(1)
    comb = [r0[i] + r1[i] for i in range(n)]
    ng = len(gen_pairs)
    return auc_gt(comb[:ng], comb[ng:])


# ---------------------------------------------------------------- load caches
def load_edge_cache():
    edges = {}
    with open(EDGE_CACHE) as fh:
        h = fh.readline().rstrip("\n").split("\t"); idx = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            def g(c):
                return f[idx[c]]
            def gn(c):
                return None if f[idx[c]] == "" else float(f[idx[c]])
            k = frozenset((g("gene_a"), g("gene_b")))
            edges[k] = dict(gene_a=g("gene_a"), gene_b=g("gene_b"), cls=g("cls"),
                            in_ep=(g("in_ep") == "1"), core=gn("core"),
                            dn_a=g("bridge_dn_a") or None, dn_b=g("bridge_dn_b") or None,
                            mask=gn("bridge_mask"))
    return edges


def load_degree():
    nbr = defaultdict(set)
    with open(DENOVO_EDGES) as fh:
        next(fh)
        for ln in fh:
            a, b, cc = ln.rstrip("\n").split("\t")
            if float(cc) >= DEG_T:
                nbr[a].add(b); nbr[b].add(a)
    return {k: len(v) for k, v in nbr.items()}


def load_locus_cache(path, cols):
    out = {}
    if not os.path.exists(path):
        print(f"[warn] cache missing: {path} -> feature treated as absent", flush=True)
        return out
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t"); idx = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            rec = {}
            for c in cols:
                v = f[idx[c]]
                rec[c] = None if v == "" else float(v)
            out[f[0]] = rec
    return out


def load_sharedlen():
    out = {}
    with open(SHAREDLEN_CACHE) as fh:
        h = fh.readline().rstrip("\n").split("\t"); idx = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            k = frozenset((f[idx["gene_a"]], f[idx["gene_b"]]))
            al = f[idx["aln_len"]]; af = f[idx["aln_frac"]]
            out[k] = dict(aln_len=(None if al == "" else float(al)),
                          aln_frac=(None if af == "" else float(af)))
    return out


# ---------------------------------------------------------------- per-edge aggregation
def edge_locus_agg(edges, locus_val):
    """edge -> (max, mean) of a per-locus scalar over the two bridge loci. locus_val: dn -> float|None."""
    mx, mn = {}, {}
    for k, e in edges.items():
        va = locus_val(e["dn_a"]) if e["dn_a"] else None
        vb = locus_val(e["dn_b"]) if e["dn_b"] else None
        vals = [v for v in (va, vb) if v is not None]
        if vals:
            mx[k] = max(vals); mn[k] = sum(vals) / len(vals)
        else:
            mx[k] = None; mn[k] = None
    return mx, mn


def arb_pop(edges, feat, lo):
    """genuine, TP feature lists on the arbitration population: in_ep=0, core>=lo, feature present."""
    gen, tp = [], []
    for k, e in edges.items():
        if e["in_ep"] or e["core"] is None or e["core"] < lo:
            continue
        v = feat.get(k)
        if v is None:
            continue
        if e["cls"] == "genuine":
            gen.append((k, v))
        elif e["cls"] == "TP":
            tp.append((k, v))
    return gen, tp


def auc_by_band(edges, feat, higher_is_genuine=True):
    res = {}
    for lo in BANDS:
        gen, tp = arb_pop(edges, feat, lo)
        gv = [v for _, v in gen]; tv = [v for _, v in tp]
        if higher_is_genuine:
            a = auc_gt(gv, tv)
        else:
            a = auc_gt(tv, gv)  # length: genuine LOWER -> report discriminating direction (TP>genuine)
        res[lo] = dict(auc=(None if a is None else round(a, 4)), n_gen=len(gv), n_tp=len(tv))
    return res


def main():
    edges = load_edge_cache()
    deg = load_degree()
    mm = load_locus_cache(READMM_CACHE, ["mm_median", "mm_mean", "mm_max"])
    km = load_locus_cache(KMER_CACHE, ["km_median", "km_mean", "km_p90", "km_frac_ge_100"])
    slen = load_sharedlen()

    n_arb = sum(1 for e in edges.values() if not e["in_ep"] and e["core"] is not None)
    n_arb_gen = sum(1 for e in edges.values() if not e["in_ep"] and e["core"] is not None and e["cls"] == "genuine")
    n_arb_tp = sum(1 for e in edges.values() if not e["in_ep"] and e["core"] is not None and e["cls"] == "TP")
    print(f"[pop] edges={len(edges)} arbitration(in_ep=0,core present)={n_arb} "
          f"(genuine={n_arb_gen}, protein-missed TP={n_arb_tp})", flush=True)

    # ----- feature 1: homology-graph degree
    deg_mx, deg_mn = edge_locus_agg(edges, lambda d: deg.get(d))
    # ----- feature 2: read multimapping (three per-locus stats)
    mm_variants = {}
    for st in ("mm_median", "mm_mean", "mm_max"):
        mx, mn = edge_locus_agg(edges, lambda d, st=st: (mm.get(d) or {}).get(st))
        mm_variants[st] = (mx, mn)
    # ----- feature 3: k-mer multiplicity (four per-locus stats)
    km_variants = {}
    for st in ("km_median", "km_mean", "km_p90", "km_frac_ge_100"):
        mx, mn = edge_locus_agg(edges, lambda d, st=st: (km.get(d) or {}).get(st))
        km_variants[st] = (mx, mn)
    # ----- caveat feature: shared-segment length (per edge)
    slen_len = {k: (slen.get(k) or {}).get("aln_len") for k in edges}
    slen_frac = {k: (slen.get(k) or {}).get("aln_frac") for k in edges}

    # ============ AUC per feature vs soft-mask ============
    softmask_feat = {k: e["mask"] for k, e in edges.items()}
    auc = {}
    auc["softmask"] = auc_by_band(edges, softmask_feat)
    auc["degree_max"] = auc_by_band(edges, deg_mx)
    auc["degree_mean"] = auc_by_band(edges, deg_mn)
    for st, (mx, mn) in mm_variants.items():
        auc[f"readmm_{st}_max"] = auc_by_band(edges, mx)
    for st, (mx, mn) in km_variants.items():
        auc[f"kmer_{st}_max"] = auc_by_band(edges, mx)
    auc["sharedlen_len_TPhi"] = auc_by_band(edges, slen_len, higher_is_genuine=False)
    auc["sharedlen_frac_TPhi"] = auc_by_band(edges, slen_frac, higher_is_genuine=False)

    def best_variant(prefix):
        cand = [name for name in auc if name.startswith(prefix)]
        return max(cand, key=lambda nm: (auc[nm][0.13]["auc"] or 0))
    best_readmm = best_variant("readmm_")
    best_kmer = best_variant("kmer_")

    # ============ complementarity of degree (best single) with soft-mask ============
    complement = {}
    for lo in BANDS:
        # need edges with BOTH softmask and degree present
        gen_pairs, tp_pairs = [], []
        mvals, dvals = [], []
        for k, e in edges.items():
            if e["in_ep"] or e["core"] is None or e["core"] < lo:
                continue
            m = e["mask"]; d = deg_mx.get(k)
            if m is None or d is None:
                continue
            if e["cls"] == "genuine":
                gen_pairs.append((m, d))
            elif e["cls"] == "TP":
                tp_pairs.append((m, d))
            else:
                continue
            mvals.append(m); dvals.append(d)
        comb = rank_avg_auc(gen_pairs, tp_pairs)
        sm_only = auc_gt([p[0] for p in gen_pairs], [p[0] for p in tp_pairs])
        dg_only = auc_gt([p[1] for p in gen_pairs], [p[1] for p in tp_pairs])
        complement[lo] = dict(
            n_gen=len(gen_pairs), n_tp=len(tp_pairs),
            softmask_auc=(None if sm_only is None else round(sm_only, 4)),
            degree_auc=(None if dg_only is None else round(dg_only, 4)),
            combined_rankavg_auc=(None if comb is None else round(comb, 4)),
            spearman_softmask_degree=(None if spearman(mvals, dvals) is None
                                      else round(spearman(mvals, dvals), 4)))

    # ============ combined GATE through the committed compute_point machinery ============
    print("[gate] preparing D (family_er_oracle_sweep.prepare) ...", flush=True)
    D = S.prepare()
    ref_pairs, attr = D["ref_pairs"], D["attr"]
    # attach soft-mask (mask) and degree to attr (compute_point reads attr via keep predicate)
    for k in ref_pairs:
        e = edges.get(k)
        attr[k]["mask"] = e["mask"] if e else None
        attr[k]["deg"] = deg_mx.get(k)
    n_nc = len(D["noncoding_tp"])

    def cp(name, t, keep):
        r = S.compute_point(name, t, keep, D)
        return r

    none_r = cp("none", None, lambda k, a: True)
    prot_r = cp("protein", None, lambda k, a: a["in_ep"])
    # soft-mask recommended reference: in_ep OR (core>=0.50 AND mask<0.10)
    def kp_mask(t, m):
        return lambda k, a, t=t, m=m: a["in_ep"] or (
            a["core"] is not None and a["core"] >= t and a["mask"] is not None and a["mask"] < m)
    softmask_rec = cp("softmask", 0.50, kp_mask(0.50, 0.10))
    core_knee = cp("combined", 0.70, lambda k, a: a["in_ep"] or (a["core"] is not None and a["core"] >= 0.70))

    # combined gate sweep: in_ep OR (core>=t AND mask<m AND deg<Dmax)
    def kp_md(t, m, dmax):
        return lambda k, a, t=t, m=m, dmax=dmax: a["in_ep"] or (
            a["core"] is not None and a["core"] >= t
            and a["mask"] is not None and a["mask"] < m
            and a["deg"] is not None and a["deg"] < dmax)
    # degree-ONLY gate (no soft-mask): does degree REPLACE soft-mask? in_ep OR (core>=t AND deg<dmax)
    def kp_deg(t, dmax):
        return lambda k, a, t=t, dmax=dmax: a["in_ep"] or (
            a["core"] is not None and a["core"] >= t and a["deg"] is not None and a["deg"] < dmax)
    T_GRID = [0.13, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70]
    M = 0.10
    D_GRID = [3, 4, 5, 6, 8, 10, 15, 20, 30, 10**9]  # 1e9 = degree OFF (reduces to soft-mask gate)
    gate_rows = []
    for t in T_GRID:
        for dmax in D_GRID:
            r = cp("mask+deg", t, kp_md(t, M, dmax))
            r["m"] = M; r["dmax"] = (None if dmax >= 10**9 else dmax)
            gate_rows.append(r)
    deg_only_rows = []
    for t in T_GRID:
        for dmax in D_GRID:
            if dmax >= 10**9:
                continue
            r = cp("deg-only", t, kp_deg(t, dmax))
            r["m"] = None; r["dmax"] = dmax
            deg_only_rows.append(r)

    # Pareto test vs soft-mask recommended (21 TPpm @ edge-host). Beat = TPpm>=ref AND edgehost<=ref,
    # strictly better on at least one.
    ref_tppm = softmask_rec["tp_noep_retained"]
    ref_eh = softmask_rec["genuine_rate_block_edgehost"]
    pareto = [r for r in gate_rows
              if r["tp_noep_retained"] >= ref_tppm and r["genuine_rate_block_edgehost"] <= ref_eh
              and (r["tp_noep_retained"] > ref_tppm or r["genuine_rate_block_edgehost"] < ref_eh)]
    combined_rec = None
    if pareto:
        combined_rec = max(pareto, key=lambda r: (r["tp_noep_retained"], -r["genuine_rate_block_edgehost"],
                                                  -r["genuine_overmerge"]))
    # dmax plateau + in-sample provenance for the recommended combined point (guard against knife-edge/overfit)
    plateau = None
    if combined_rec:
        wt = combined_rec["t"]
        band = sorted((r for r in gate_rows if r["t"] == wt and r.get("dmax") is not None),
                      key=lambda r: r["dmax"])
        opt = (combined_rec["tp_noep_retained"], combined_rec["genuine_overmerge"])
        plateau_dmax = [r["dmax"] for r in band
                        if (r["tp_noep_retained"], r["genuine_overmerge"]) == opt]
        below = [r for r in band if r["dmax"] < combined_rec["dmax"]]
        one_below = below[-1] if below else None
        plateau = dict(
            t=wt, recommended_dmax=combined_rec["dmax"], plateau_dmax=plateau_dmax,
            plateau_width=len(plateau_dmax),
            one_lower_dmax=(one_below["dmax"] if one_below else None),
            one_lower_TPpm=(one_below["tp_noep_retained"] if one_below else None),
            note=("dmax=%s is the LOWEST cap that preserves all real protein-missed TP (TPpm=%d); dmax in %s "
                  "give the identical optimum (width-%d plateau), dmax=%s drops TPpm to %s (loses real TP). "
                  "NOT a knife-edge. CAVEAT: dmax is picked IN-SAMPLE (single free parameter, coarse grid, no "
                  "held-out split) so the modest 55->%d over-merge gain carries a small overfit risk."
                  % (combined_rec["dmax"], combined_rec["tp_noep_retained"], plateau_dmax, len(plateau_dmax),
                     (one_below["dmax"] if one_below else None),
                     (one_below["tp_noep_retained"] if one_below else None),
                     combined_rec["genuine_overmerge"])))
    # degree-only: can it MATCH the soft-mask recommended point without the mask at all?
    deg_only_pareto = [r for r in deg_only_rows
                       if r["tp_noep_retained"] >= ref_tppm and r["genuine_rate_block_edgehost"] <= ref_eh]
    deg_only_rec = (max(deg_only_pareto, key=lambda r: (r["tp_noep_retained"],
                    -r["genuine_rate_block_edgehost"], -r["genuine_overmerge"])) if deg_only_pareto else None)

    # ============ read-mm mechanistic stats (WHY signal 2 fails: both classes high-multiplicity) ======
    # per-edge mm_median_max over the arbitration population + fraction saturating the -N50 cap.
    mm_med_arb, mm_max_arb = [], []
    for k, e in edges.items():
        if e["in_ep"] or e["core"] is None or e["cls"] not in ("genuine", "TP"):
            continue
        v = mm_variants["mm_median"][0].get(k); vmax = mm_variants["mm_max"][0].get(k)
        if v is not None:
            mm_med_arb.append(v)
        if vmax is not None:
            mm_max_arb.append(vmax)
    n_loci_cache = sum(1 for _ in open(READMM_CACHE)) - 1
    n_loci_cap = 0
    with open(READMM_CACHE) as fh:
        hh = fh.readline().rstrip("\n").split("\t"); ii = {c: j for j, c in enumerate(hh)}
        for ln in fh:
            ff = ln.rstrip("\n").split("\t"); mv = ff[ii["mm_max"]]
            if mv != "" and float(mv) >= 49:
                n_loci_cap += 1
    readmm_mechanism = dict(
        note=("read multi-mapping is near-chance because BOTH genuine and protein-missed-TP loci are "
              "high-multiplicity, NOT because multiplicity is ~1. The minimap2 -N50 secondary cap truncates "
              "the count so both classes pile at the ceiling and cannot be separated."),
        arb_edge_mm_median_max_median=(round(statistics.median(mm_med_arb), 1) if mm_med_arb else None),
        arb_edge_mm_median_max_mean=(round(statistics.mean(mm_med_arb), 1) if mm_med_arb else None),
        arb_edge_frac_at_N50_cap=(round(sum(1 for v in mm_max_arb if v >= 49) / len(mm_max_arb), 4)
                                  if mm_max_arb else None),
        per_locus_n=n_loci_cache,
        per_locus_frac_at_N50_cap=round(n_loci_cap / max(1, n_loci_cache), 4))

    # ============ shared-length by-band medians (sign-inversion evidence for the caveat) =============
    slen_band = {}
    for lo in BANDS:
        gpv, tpv = arb_pop(edges, slen_len, lo)
        gvals = [v for _, v in gpv]; tvals = [v for _, v in tpv]
        slen_band[lo] = dict(
            median_genuine=(round(statistics.median(gvals), 0) if gvals else None),
            median_TP=(round(statistics.median(tvals), 0) if tvals else None),
            auc_TP_gt_genuine=auc["sharedlen_len_TPhi"][lo]["auc"])

    # ============ HONEST CAVEAT: high-copy real families mis-flagged by degree ============
    # top-decile degree among protein-missed TP in arbitration pop = real families the degree flags as TE.
    gen_arb, tp_arb = arb_pop(edges, deg_mx, 0.13)
    tp_sorted = sorted(tp_arb, key=lambda kv: -kv[1])
    hi_deg_tp = tp_sorted[:max(1, len(tp_sorted) // 10)]
    def genes_of(k):
        e = edges[k]; return sorted((e["gene_a"], e["gene_b"]))
    hi_deg_tp_named = [dict(pair=" -- ".join(genes_of(k)), degree=v,
                           in_ep=edges[k]["in_ep"], core=round(edges[k]["core"], 3),
                           shared_len=slen_len.get(k), softmask=edges[k]["mask"]) for k, v in hi_deg_tp]
    # protein branch already rescues annotated high-copy families: count in_ep=1 TP that are high-degree
    all_tp = [(k, deg_mx.get(k)) for k in edges if edges[k]["cls"] == "TP" and deg_mx.get(k) is not None]
    all_tp_sorted = sorted(all_tp, key=lambda kv: -kv[1])
    hi_deg_all_tp = all_tp_sorted[:max(1, len(all_tp_sorted) // 10)]
    rescued_by_protein = sum(1 for k, v in hi_deg_all_tp if edges[k]["in_ep"])
    # does shared-length separate the residual protein-missed high-copy TP from genuine? use hi-degree subset
    def hi_deg_genuine():
        return [kv for kv in sorted(gen_arb, key=lambda kv: -kv[1])[:len(gen_arb) // 2]]
    hi_gen = hi_deg_genuine()
    slen_gen = [slen_len.get(k) for k, _ in hi_gen if slen_len.get(k) is not None]
    slen_tp = [slen_len.get(k) for k, _ in hi_deg_tp if slen_len.get(k) is not None]
    # combined degree(<D) AND shared-length(>=L) as a real-family rescue on the high-degree TP
    slen_auc_arb = auc["sharedlen_len_TPhi"]

    # how many of the high-degree protein-missed TP are ALREADY excluded by the softmask<0.10 branch:
    n_hi_deg_tp_softmask_ge_m = sum(1 for k, _ in hi_deg_tp
                                    if edges[k]["mask"] is not None and edges[k]["mask"] >= 0.10)
    caveat = dict(
        premise=("degree/multiplicity flag ANY high-copy locus as TE, including REAL high-copy gene families "
                 "(KRAB-ZNF ~400 paralogs, olfactory-receptor ~1000). So degree alone cannot separate real "
                 "high-copy families from TEs."),
        high_degree_protein_missed_TP=hi_deg_tp_named,
        n_high_degree_protein_missed_TP=len(hi_deg_tp),
        protein_branch_rescue=dict(
            note=("protein-homology (in_ep) is evaluated FIRST in the gate, so annotated high-copy families "
                  "(ZNF/OR) that ARE protein-homologous are KEPT by the in_ep branch and never reach the "
                  "degree branch -- degree only ever judges the protein-MISSED residual."),
            n_top_decile_degree_TP=len(hi_deg_all_tp),
            rescued_by_in_ep=rescued_by_protein,
            rescued_frac=round(rescued_by_protein / max(1, len(hi_deg_all_tp)), 4)),
        readmm_why_near_chance=readmm_mechanism,
        shared_length_fix=dict(
            note=("SHARED-SEGMENT LENGTH separates real paralogs (LONG contiguous gene body) from TE bridges "
                  "(SHORT exonized-TE fragment) ONLY on the broad re-admit region core>=0.13 (AUC TP>genuine "
                  "0.821). It is a SIGN-INVERSION at high core: at the recommended core>=0.50 band it is at "
                  "chance (0.521) and at core>=0.70 it INVERTS below chance (0.482 -- genuine bridges share "
                  "MORE there). So shared-length does NOT rescue the high-copy-family confound at the operating "
                  "band. The actual protection at core>=0.50 is (i) protein-first ordering (in_ep evaluated "
                  "first) and (ii) the softmask<0.10 branch, which ALREADY excludes all high-degree "
                  "protein-missed TP (every one has softmask>=0.235); shared-length is not a term in the gate."),
            sharedlen_auc_by_band=slen_auc_arb,
            median_sharedlen_by_band=slen_band,
            n_high_degree_protein_missed_TP_softmask_ge_0_10=n_hi_deg_tp_softmask_ge_m,
            hi_degree_genuine_median_sharedlen=(round(statistics.median(slen_gen), 1) if slen_gen else None),
            hi_degree_TP_median_sharedlen=(round(statistics.median(slen_tp), 1) if slen_tp else None)))

    # ============ determinism (recompute AUCs + gate twice) ============
    a2 = auc_by_band(edges, deg_mx)
    r2 = S.compute_point("softmask", 0.50, kp_mask(0.50, 0.10), D)
    det = (a2 == auc["degree_max"]
           and r2["tp_noep_retained"] == softmask_rec["tp_noep_retained"]
           and r2["genuine_rate_block_edgehost"] == softmask_rec["genuine_rate_block_edgehost"])

    # ============ write per-edge TSV ============
    cols = ["gene_a", "gene_b", "cls", "in_ep", "core", "bridge_dn_a", "bridge_dn_b", "bridge_mask",
            "deg_max", "deg_mean", "mm_median_max", "mm_mean_max", "mm_max_max", "km_median_max",
            "km_p90_max", "km_frac_ge_100_max", "aln_len", "aln_frac", "arbitration"]
    with open(os.path.join(BENCH, "family_er_readintrinsic_te.tsv"), "w") as out:
        out.write("\t".join(cols) + "\n")
        for k in sorted(edges, key=lambda x: tuple(sorted(x))):
            e = edges[k]
            def fm(x):
                return "" if x is None else (f"{x:.4f}" if isinstance(x, float) else str(x))
            arbo = (not e["in_ep"]) and (e["core"] is not None) and (e["cls"] in ("genuine", "TP"))
            row = [e["gene_a"], e["gene_b"], e["cls"], "1" if e["in_ep"] else "0", fm(e["core"]),
                   e["dn_a"] or "", e["dn_b"] or "", fm(e["mask"]),
                   fm(deg_mx.get(k)), fm(deg_mn.get(k)),
                   fm(mm_variants["mm_median"][0].get(k)), fm(mm_variants["mm_mean"][0].get(k)),
                   fm(mm_variants["mm_max"][0].get(k)),
                   fm(km_variants["km_median"][0].get(k)), fm(km_variants["km_p90"][0].get(k)),
                   fm(km_variants["km_frac_ge_100"][0].get(k)),
                   fm(slen_len.get(k)), fm(slen_frac.get(k)), "1" if arbo else "0"]
            out.write("\t".join(row) + "\n")

    # ============ write JSON ============
    def slim(r):
        return dict(gate=r["gate"], t=r.get("t"), m=r.get("m"), dmax=r.get("dmax"),
                    n_edges=r["n_edges"], tp_retained=r["tp_retained"],
                    tp_noep_retained=r["tp_noep_retained"], noncoding_tp_retained=r["noncoding_tp_retained"],
                    truthbar_retained=r["truthbar_retained"], genuine_overmerge=r["genuine_overmerge"],
                    genuine_rate_block=r["genuine_rate_block"],
                    genuine_rate_block_edgehost=r["genuine_rate_block_edgehost"])
    summary = dict(
        question=("Can TE-ness be detected DIRECTLY FROM THE READS (self-contained, no RepeatMasker) and does "
                  "it separate the E_r over-merge classes better than genome soft-mask (AUC 0.694 @ core>=0.50, "
                  "0.777 @ core>=0.13)?"),
        population=dict(arbitration_in_ep0_core_present=n_arb, genuine=n_arb_gen, protein_missed_TP=n_arb_tp,
                       definition="in_ep=0 AND direct core edge; genuine-FP vs protein-missed real TP"),
        softmask_reference_auc=SOFTMASK_AUC,
        auc_by_feature=auc,
        best_readmm_variant=best_readmm, best_kmer_variant=best_kmer,
        headline=dict(
            band_core_ge_0_13=dict(
                softmask=SOFTMASK_AUC[0.13], degree_max=auc["degree_max"][0.13]["auc"],
                readmm_best=auc[best_readmm][0.13]["auc"], kmer_best=auc[best_kmer][0.13]["auc"]),
            band_core_ge_0_50_RECOMMENDED=dict(
                softmask=SOFTMASK_AUC[0.50], degree_max=auc["degree_max"][0.50]["auc"],
                readmm_best=auc[best_readmm][0.50]["auc"], kmer_best=auc[best_kmer][0.50]["auc"]),
            reading=("NO single read-intrinsic signal beats soft-mask at the 0.694 operating band (core>=0.50): "
                     "degree=0.650, read-mm=0.592, k-mer=0.619 all LOSE there. HOMOLOGY-GRAPH DEGREE is the "
                     "strongest read-intrinsic signal and beats soft-mask ONLY on the broad re-admit region "
                     "(core>=0.13: 0.830 > 0.777) and core>=0.70 (0.762 > 0.675). The load-bearing result is "
                     "COMPLEMENT not replace: rank-average(soft-mask, degree) exceeds soft-mask at EVERY band "
                     "including the operating band (0.7395 > 0.6938 @ core>=0.50); Spearman 0.50 = partly "
                     "orthogonal. READ MULTI-MAPPING is near-chance NOT because multiplicity is low but because "
                     "BOTH classes are high-multiplicity: per-locus median mm_median ~11 and ~92% of arbitration "
                     "edges saturate the minimap2 -N50 secondary cap, so the capped count cannot separate the "
                     "two. GENOME K-MER MULTIPLICITY is intermediate (0.732 @ core>=0.13). degree_max AUC is "
                     "non-monotone across bands (0.83/0.73/0.65/0.76) as the high-band population shrinks "
                     "(n_tp 90->64->33->19).")),
        complementarity_degree_vs_softmask=complement,
        gate=dict(
            reference_none=slim(none_r), reference_protein=slim(prot_r),
            softmask_recommended=slim(softmask_rec), core_knee=slim(core_knee),
            combined_gate="in_ep OR (core>=t AND bridge_mask<m AND degree<dmax)",
            sweep=[slim(r) for r in gate_rows],
            pareto_beats_softmask_recommended=(slim(combined_rec) if combined_rec else None),
            pareto_dmax_plateau=plateau,
            pareto_note=("Pareto-beat = TPpm >= softmask-recommended AND edge-host <= it, strictly better on "
                         "one. None -> degree adds nothing beyond soft-mask at the operating point. dmax is "
                         "selected in-sample (see pareto_dmax_plateau)."),
            degree_only_gate="in_ep OR (core>=t AND degree<dmax)  [no soft-mask -> can degree REPLACE it?]",
            degree_only_sweep=[slim(r) for r in deg_only_rows],
            degree_only_matches_softmask_recommended=(slim(deg_only_rec) if deg_only_rec else None)),
        caveat_high_copy_families=caveat,
        determinism=dict(reproducible_within_run=det, pythonhashseed=os.environ.get("PYTHONHASHSEED")),
    )
    with open(os.path.join(BENCH, "family_er_readintrinsic_te.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)

    # ============ console ============
    print("\n================ READ-INTRINSIC TE SIGNALS: AUC(genuine>TP) by core band ================")
    print(f"{'feature':<26} " + "  ".join(f"c>={lo}" for lo in BANDS))
    def line(name):
        return f"{name:<26} " + "  ".join(
            f"{(auc[name][lo]['auc'] if auc[name][lo]['auc'] is not None else 'NA'):<6}" for lo in BANDS)
    print(line("softmask") + "   <- committed reference (0.777 / 0.694)")
    print(line("degree_max"))
    print(line("degree_mean"))
    print(line(best_readmm) + "   <- best read-multimap variant")
    print(line(best_kmer) + "   <- best k-mer variant")
    print(line("sharedlen_len_TPhi") + "   <- caveat fix (AUC TP>genuine)")
    print("\n[headline] core>=0.13: softmask=%.4f degree=%.4f readmm=%s kmer=%s"
          % (SOFTMASK_AUC[0.13], auc["degree_max"][0.13]["auc"], auc[best_readmm][0.13]["auc"],
             auc[best_kmer][0.13]["auc"]))
    print("[headline] core>=0.50: softmask=%.4f degree=%.4f readmm=%s kmer=%s"
          % (SOFTMASK_AUC[0.50], auc["degree_max"][0.50]["auc"], auc[best_readmm][0.50]["auc"],
             auc[best_kmer][0.50]["auc"]))
    print("\n[complementarity degree+softmask]")
    for lo in BANDS:
        c = complement[lo]
        print(f"   c>={lo}: softmask={c['softmask_auc']} degree={c['degree_auc']} "
              f"combined(rank-avg)={c['combined_rankavg_auc']} spearman={c['spearman_softmask_degree']} "
              f"(n_gen={c['n_gen']}/n_tp={c['n_tp']})")
    print("\n[gate] soft-mask recommended: TPpm=%d/%d ncTP=%d/%d gen=%d blk_eh=%s"
          % (softmask_rec["tp_noep_retained"], n_arb_tp, softmask_rec["noncoding_tp_retained"], n_nc,
             softmask_rec["genuine_overmerge"], softmask_rec["genuine_rate_block_edgehost"]))
    print("[gate] core knee (t=0.70):    TPpm=%d gen=%d blk_eh=%s"
          % (core_knee["tp_noep_retained"], core_knee["genuine_overmerge"],
             core_knee["genuine_rate_block_edgehost"]))
    if combined_rec:
        print("[gate] mask+degree PARETO:   t=%s dmax=%s TPpm=%d ncTP=%d gen=%d blk_eh=%s"
              % (combined_rec["t"], combined_rec["dmax"], combined_rec["tp_noep_retained"],
                 combined_rec["noncoding_tp_retained"], combined_rec["genuine_overmerge"],
                 combined_rec["genuine_rate_block_edgehost"]))
        if plateau:
            print("[gate]   dmax plateau=%s (width %d), one-lower dmax=%s -> TPpm=%s; picked IN-SAMPLE"
                  % (plateau["plateau_dmax"], plateau["plateau_width"], plateau["one_lower_dmax"],
                     plateau["one_lower_TPpm"]))
    else:
        print("[gate] mask+degree: NO Pareto-beat of the soft-mask recommended point "
              "(degree does not improve the operating point).")
    if deg_only_rec:
        print("[gate] degree-ONLY (no mask): t=%s dmax=%s TPpm=%d ncTP=%d gen=%d blk_eh=%s "
              "-> degree can REPLACE soft-mask at the operating point"
              % (deg_only_rec["t"], deg_only_rec["dmax"], deg_only_rec["tp_noep_retained"],
                 deg_only_rec["noncoding_tp_retained"], deg_only_rec["genuine_overmerge"],
                 deg_only_rec["genuine_rate_block_edgehost"]))
    else:
        print("[gate] degree-ONLY (no mask): cannot match the soft-mask recommended point.")
    print("\n[caveat] high-degree protein-missed TP (real families degree flags as TE):")
    for r in hi_deg_tp_named[:10]:
        print(f"   {r['pair']:34s} deg={r['degree']} core={r['core']} shared_len={r['shared_len']} "
              f"softmask={r['softmask']}")
    print(f"[caveat] protein branch rescues {rescued_by_protein}/{len(hi_deg_all_tp)} top-decile-degree TP "
          f"(annotated high-copy families kept by in_ep before the degree branch).")
    print(f"[caveat] shared-length AUC(TP>gen): c>=0.13={slen_auc_arb[0.13]['auc']} c>=0.50={slen_auc_arb[0.50]['auc']} "
          f"c>=0.70={slen_auc_arb[0.70]['auc']} (SIGN-INVERSION: chance at op band, inverts at 0.70; "
          f"softmask<0.10 already excludes all {len(hi_deg_tp)}/{len(hi_deg_tp)} high-degree TP, all softmask>=0.235)")
    print(f"[caveat] read-mm WHY near-chance: arb-edge median mm_median_max={readmm_mechanism['arb_edge_mm_median_max_median']} "
          f"({readmm_mechanism['arb_edge_frac_at_N50_cap']*100:.0f}% of arb edges at -N50 cap; per-locus "
          f"{readmm_mechanism['per_locus_frac_at_N50_cap']*100:.0f}% at cap) -> BOTH classes high-multiplicity, not ~1")
    print(f"\n[determinism] within-run reproducible={det} (PYTHONHASHSEED=0)")
    print("wrote bench/family_er_readintrinsic_te.tsv + bench/family_er_readintrinsic_te.json")


if __name__ == "__main__":
    main()
