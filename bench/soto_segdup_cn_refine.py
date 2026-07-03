#!/usr/bin/env python3
"""bench/soto_segdup_cn_refine.py -- CHARACTERIZE stage.

Soto-2025-style DNA-CALIBRATED refine of the RNA family catalog, DNA-informed (NOT the RNA-only
gate). Soto's family definition = SD98 (segmental duplications at >= 98% identity) + keep exon
matches covering >99% of each exon + famCN refine (group only genes whose read-depth copy number
mean-abs-deviation < 1). The discriminator that kills a DOMAIN-SHARE (two unrelated genes sharing
ONE exon, e.g. RHD+SDHD) is NOT the exon-by-exon test (a domain-share shares a full exon) -- it is
(1) the SD98/segdup restriction (one shared exon is not an extended >=1kb >=90% segmental
DUPLICATION) and (2) famCN copy-number consistency. Both are DNA.

We have the gorilla analogs:
  * SEDEF segdup map  = /mnt/c/Users/jfris/Desktop/final.bed  (SD98 analog; final.bed = genome
    self-alignment, one line = one segdup PAIR, region A cols 0-2 / partner B cols 3-5, identity =
    fracMatch col 21 [1-indexed], fracMatchIndel col 22).
  * mGorGor1 diploid CN = bench/diploid_cn_oracle.tsv (famCN analog; per multi-copy family
    dominant-gene diploid/haploid CN).

TEST -- does the DNA duplication signal cleanly SPLIT the domain-share FP residual from real
paralog families, and at what recall cost on DIVERGENT real paralogs?
  * DOMAIN-SHARE FP (target)         = ep_impure genuine pairs (RHD+SDHD, RBBP4+SYNC, DEDD+NIT1,
    ...). EXPECT: their two loci are NOT a SEDEF segdup pair (one shared exon != extended dup).
  * REAL near-identical paralogs     = class TP (DNA-confirmed). EXPECT: their loci ARE a SEDEF pair.
  * REAL DIVERGENT paralogs (DECISIVE)= class truthbar (E_p-protein-confirmed, DNA did not confirm).
    The DNA-level version of the RNA erosion confound: what fraction is captured in SEDEF vs eroded
    below its ~90% identity floor?  Stratified by SEDEF identity floor 0.90 / 0.95 / 0.98.

A pair (dn_a, dn_b) -> two genome loci (chrom/start/end via denovo_transcripts.meta.tsv). SEDEF
segdup-support (primary) = the two loci sit on the two ends of ONE SEDEF segdup pair (locus A
duplicated to locus B), recorded at that pair's identity. Weaker signals also tracked: same-region
(both loci inside one segdup unit) and each-locus-in-any-segdup. CN-consistency = both genes resolve
to the multi-copy CN oracle and their famCN mean-abs-deviation < 1 (|asm_hapCN_a - asm_hapCN_b| <= 1).

HONESTY: DNA-informed (SEDEF + diploid CN), a Soto-style DNA-calibrated REFINE / validation, like
the diploid oracle -- NOT the RNA-only edge rule. Deterministic (PYTHONHASHSEED=0).
Pair substrate reused verbatim from bench/colinear_multiexon_gate.tsv (same 5571 direct multi-exon
E_r edges: 416 TP + 3382 truthbar + 1773 genuine).
Out: bench/soto_segdup_cn_refine.tsv (per pair) + .json (summary) + stdout.
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import bisect
import collections
import csv
import json

BENCH = os.path.dirname(os.path.abspath(__file__))

SEDEF_CANDIDATES = ["/mnt/c/Users/jfris/Desktop/final.bed",
                    "/mnt/c/Users/jfris/Desktop/Rustle/final.bed",
                    "/home/juanfra/winloci_scratch/final.bed"]
SEDEF = next((p for p in SEDEF_CANDIDATES if os.path.exists(p)), None)
META_TSV = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
PAIRS_TSV = os.path.join(BENCH, "colinear_multiexon_gate.tsv")
CN_TSV = os.path.join(BENCH, "diploid_cn_oracle.tsv")
OUT_TSV = os.path.join(BENCH, "soto_segdup_cn_refine.tsv")
OUT_JSON = os.path.join(BENCH, "soto_segdup_cn_refine.json")

FLOORS = [0.90, 0.95, 0.98]          # SEDEF identity floors toward SD98
BACKSCAN = 5_000_000                 # segdups are start-sorted; stop after a long non-overlapping run

# named domain-share exemplars (task roster; same as colinear_multiexon_gate)
DOMAIN_SHARE_EXEMPLARS = [("MOV10", "RHOC"), ("RBBP4", "SYNC"), ("RHD", "SDHD"),
                          ("DEDD", "NIT1"), ("KMT2C", "TMEM128"),
                          ("LOC129532935", "SORT1"), ("ABITRAM", "CTNNAL1"),
                          ("DIMT1", "KIF2A")]

# MAGE / GSTM2 control gene sets (cardinality controls; same as colinear_multiexon_gate)
GSTM2_GENES = {"GSTM2", "LOC115930164", "LOC115930576", "LOC101129940"}
MAGE_GENES = {"MAGEA9", "LOC129529978", "LOC129529986"}


# --------------------------------------------------------------------------- SEDEF segdup map
def load_sedef(path):
    """per-chrom start-sorted intervals (start, end, pid, side, id_fm, id_fi); each final.bed line =
       one segdup PAIR. id_fm = fracMatch (col 21 1-idx), id_fi = fracMatchIndel (col 22)."""
    by = collections.defaultdict(list)
    pid = 0
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 22:
                continue
            try:
                ca, sa, ea = f[0], int(f[1]), int(f[2])
                cb, sb, eb = f[3], int(f[4]), int(f[5])
                idfm = float(f[20]); idfi = float(f[21])
            except (ValueError, IndexError):
                continue  # header / malformed
            by[ca].append((sa, ea, pid, 0, idfm, idfi))
            by[cb].append((sb, eb, pid, 1, idfm, idfi))
            pid += 1
    starts = {}
    for c in by:
        by[c].sort()
        starts[c] = [iv[0] for iv in by[c]]
    return by, starts, pid


def overlaps(by, starts, chrom, s, e):
    """set of (pid, side, id_fm, id_fi) whose segdup interval overlaps [s,e) on chrom."""
    out = []
    arr = by.get(chrom)
    if not arr:
        return out
    st = starts[chrom]
    i = bisect.bisect_right(st, e)            # intervals starting at/before e
    j = i - 1
    while j >= 0:
        iv_s, iv_e, pid, side, idfm, idfi = arr[j]
        if iv_e <= s:
            if s - iv_s > BACKSCAN:
                break
            j -= 1
            continue
        if iv_s < e and iv_e > s:
            out.append((pid, side, idfm, idfi))
        j -= 1
    return out


def segdup_support(by, starts, locA, locB):
    """Test whether locus A is duplicated to locus B in SEDEF.
       Returns dict: cross_n, cross_best_idfm/idfi (A on one end, B on the other of ONE pair);
       same_region (both loci inside one segdup unit / same side); a_in_segdup / b_in_segdup."""
    hitA = overlaps(by, starts, *locA)
    hitB = overlaps(by, starts, *locB)
    a_sides = collections.defaultdict(set)   # pid -> {sides A hits}
    a_id = {}
    for pid, side, idfm, idfi in hitA:
        a_sides[pid].add(side); a_id[pid] = (idfm, idfi)
    b_sides = collections.defaultdict(set)
    for pid, side, idfm, idfi in hitB:
        b_sides[pid].add(side)
    cross_ids = []      # (idfm, idfi) of pids where A and B sit on OPPOSITE ends
    same_ids = []       # (idfm, idfi) of pids where A and B sit on the SAME end (one segdup unit)
    for pid in set(a_sides) & set(b_sides):
        sa = a_sides[pid]; sb = b_sides[pid]
        if (0 in sa and 1 in sb) or (1 in sa and 0 in sb):
            cross_ids.append(a_id[pid])
        if sa & sb:
            same_ids.append(a_id[pid])
    cross_best = max((x[0] for x in cross_ids), default=0.0)
    cross_best_fi = max((x[1] for x in cross_ids), default=0.0)
    same_best = max((x[0] for x in same_ids), default=0.0)
    return dict(
        cross_n=len(cross_ids),
        cross_best_idfm=round(cross_best, 4),
        cross_best_idfi=round(cross_best_fi, 4),
        same_region=int(bool(same_ids)),
        same_best_idfm=round(same_best, 4),
        a_in_segdup=int(bool(hitA)),
        b_in_segdup=int(bool(hitB)),
    )


# --------------------------------------------------------------------------- coords + CN oracle
def load_meta(path):
    meta = {}
    with open(path) as f:
        next(f)
        for line in f:
            t = line.rstrip("\n").split("\t")
            meta[t[0]] = (t[1], int(t[2]), int(t[3]))
    return meta


def load_cn(path):
    """gene -> best (max n_loci_ref) oracle entry: dipCN, asm_hapCN, n_loci_ref, chi_H, class."""
    gene = {}
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            g = row["gene"]
            if not g or g == "NA":
                continue
            try:
                nl = int(row["n_loci_ref"])
                dip = int(row["hap_CN_mat"]) + int(row["hap_CN_pat"])
                asm = int(row["asm_hapCN"])
                chi = int(row["chi_H"]) if row["chi_H"] not in ("NA", "") else -1
            except (ValueError, KeyError):
                continue
            rec = dict(dipCN=dip, asm_hapCN=asm, n_loci_ref=nl, chi_H=chi, cls=row["class"])
            if g not in gene or nl > gene[g]["n_loci_ref"]:
                gene[g] = rec
    return gene


def cn_consistency(cn, ga, gb):
    """Soto famCN refine: consistent iff both genes resolve to the CN oracle and their famCN
       (asm_hapCN) mean-abs-deviation < 1 (== |asm_hapCN_a - asm_hapCN_b| <= 1 for integers)."""
    a = cn.get(ga); b = cn.get(gb)
    if a is None and b is None:
        return "NA_both", None, None, None
    if a is None:
        return "NA_a", None, b["asm_hapCN"], None
    if b is None:
        return "NA_b", a["asm_hapCN"], None, None
    diff = abs(a["asm_hapCN"] - b["asm_hapCN"])
    status = "consistent" if diff <= 1 else "inconsistent"
    return status, a["asm_hapCN"], b["asm_hapCN"], diff


# --------------------------------------------------------------------------- summarise helpers
def pct(n, d):
    return round(100 * n / d, 2) if d else None


def capture_block(rows, floor_key):
    """for a list of pair-rows: SEDEF cross-support capture rates + weaker signals."""
    n = len(rows)
    if not n:
        return dict(n=0)
    cross_any = sum(1 for r in rows if r["cross_n"] >= 1)
    out = dict(n=n,
               cross_support_any_id_pct=pct(cross_any, n),
               both_in_segdup_pct=pct(sum(1 for r in rows if r["a_in_segdup"] and r["b_in_segdup"]), n),
               same_region_pct=pct(sum(1 for r in rows if r["same_region"]), n),
               no_segdup_at_all_pct=pct(sum(1 for r in rows
                                            if not (r["a_in_segdup"] or r["b_in_segdup"])), n))
    for fl in FLOORS:
        k = f"cross_support_ge_{int(fl*100)}_pct"
        out[k] = pct(sum(1 for r in rows if r["cross_best_idfm"] >= fl), n)
    return out


# =========================================================================== REFINE STAGE
# Apply the Soto-style DNA-calibrated refine to the SHIPPED catalog and re-measure the family-level
# truths. KEEP a cross-gene family edge iff the gene-pair is SEDEF-segdup-supported at identity
# >= floor (sweep 0.90/0.95/0.98 toward SD98); optional famCN (diploid-CN) leg also cuts CN-
# inconsistent edges. Mirrors bench/colinear_multiexon_gate.gate_stage EXACTLY (locus subgraph per
# catalog family + same-gene backbone guard; re-derive components; drop <2 distinct loci; re-eval
# via shipped eval_partition + oracle_residuals + gene-projection relabel) so the segdup refine and
# the RNA-only colinear gate are apples-to-apples. DNA-INFORMED, NOT the RNA-only edge rule.
def refine_stage(out_rows):
    import networkx as nx
    import family_level_pr_current as C
    import graph_def_refine_sweep as SW
    import rna_only_edge_oracle as RN
    import genome_family_def as G
    import colinear_multiexon_gate as CG   # reuse _build_block_graph / _partition / _gene_families / _eval_family

    # per-gene-pair segdup support = MAX cross_best_idfm over the gene-pair's locus-pairs (a gene-pair
    # is segdup-supported if ANY of its locus-pairs is a SEDEF segdup pair); + worst CN status.
    seg = {}
    cn_by_gp = {}   # gp -> (rank, status, cn_a, cn_b) worst-case (inconsistent dominates)
    for r in out_rows:
        gp = frozenset((r["gene_a"], r["gene_b"]))
        v = r["cross_best_idfm"]
        if gp not in seg or v > seg[gp]:
            seg[gp] = v
        rank = {"inconsistent": 2, "consistent": 1}.get(r["cn_status"], 0)
        if gp not in cn_by_gp or rank > cn_by_gp[gp][0]:
            cn_by_gp[gp] = (rank, r["cn_status"], r.get("cn_a", ""), r.get("cn_b", ""))

    ctx, raw_fams, edge_pairs = C.build_ctx()
    genes = ctx["genes"]; gene_of_dn = ctx["gene_of_dn"]
    Gr = ctx["G"]
    cat_blocks, fam_ids, fam_dom, member_gene = C.load_catalog_blocks()
    mc_oracle = C.oracle_multicopy_genes(ctx)
    n_mc = len(mc_oracle)

    sub = {frozenset((r["gene_a"], r["gene_b"])): r for r in out_rows}

    def partition_refine(floor, use_cn):
        """cut cross-gene edge iff (segdup support < floor) OR (use_cn and famCN-inconsistent).
           floor=None -> segdup leg inert (isolate CN leg). Gene-pairs not in substrate are NOT
           acted on (mirrors the colinear gate; substrate covers ~1812/1817 intra-fam cross edges)."""
        out = []
        for b in cat_blocks:
            H = CG._build_block_graph(b, gene_of_dn, Gr, nx)
            rm = []
            for u, v, d in H.edges(data=True):
                if d["kind"] != "cross":
                    continue
                gp = d["gp"]
                if gp not in seg:
                    continue
                cut = (floor is not None and seg[gp] < floor)
                if use_cn and cn_by_gp.get(gp, (0,))[0] == 2:   # famCN inconsistent
                    cut = True
                if cut:
                    rm.append((u, v))
            H.remove_edges_from(rm)
            for comp in nx.connected_components(H):
                if G.distinct_loci(comp, genes) >= 2:
                    out.append(set(comp))
        return out

    # ungated baseline (same operationalization as colinear gate -> isolates the refine delta)
    base_blocks = CG._partition(cat_blocks, gene_of_dn, Gr, genes, {}, None, nx, G)
    base_gf = CG._gene_families(base_blocks, gene_of_dn)
    base_ev = CG._eval_family(base_blocks, ctx, C, SW, RN, G, mc_oracle)
    base_R = base_ev["R_oracle"]; base_recov = base_ev["recovered_mc"]

    def comem(gf, ga, gb):
        return bool(gf.get(ga, set()) & gf.get(gb, set()))
    base_co = {gp: comem(base_gf, r["gene_a"], r["gene_b"]) for gp, r in sub.items()}
    n_base_co = sum(base_co.values())

    def ctrl_fams(gf_, gset):
        fams = set()
        for g in gset:
            fams |= gf_.get(g, set())
        return dict(n_families=len(fams), genes_present=sorted(g for g in gset if g in gf_))
    gstm2_base = ctrl_fams(base_gf, GSTM2_GENES)
    mage_base = ctrl_fams(base_gf, MAGE_GENES)

    def measure(blocks, label):
        gf = CG._gene_families(blocks, gene_of_dn)
        ev = CG._eval_family(blocks, ctx, C, SW, RN, G, mc_oracle)
        # pair-level SPLITS: co-membered in baseline, NOT co-membered after refine
        split = collections.defaultdict(list)
        for gp, r in sub.items():
            if base_co[gp] and not comem(gf, r["gene_a"], r["gene_b"]):
                split[r["cls"]].append(r)
        fp = split["genuine"]; tp = split["TP"]; tb = split["truthbar"]
        ep_fp = [r for r in fp if r["ep_impure"]]
        div_tb = [r for r in tb if r["mean_matched_id"] < 0.90]
        real_split = len(tp) + len(tb)
        total_split = len(fp) + real_split
        # named rosters
        named_fp_removed = []
        exemplar_status = {}
        for p in DOMAIN_SHARE_EXEMPLARS:
            gp = frozenset(p)
            if gp not in sub:
                exemplar_status[f"{p[0]}+{p[1]}"] = "absent"
            elif not base_co[gp]:
                exemplar_status[f"{p[0]}+{p[1]}"] = f"not-comembered(sd={seg.get(gp,0.0):.3f})"
            else:
                kept = comem(gf, *p)
                exemplar_status[f"{p[0]}+{p[1]}"] = (
                    f"{'kept' if kept else 'SPLIT'}(sd={seg.get(gp,0.0):.3f},{sub[gp]['cls']})")
                if not kept:
                    named_fp_removed.append(f"{p[0]}+{p[1]}(sd={seg.get(gp,0.0):.3f})")
        ep_fp_named = sorted(f"{r['gene_a']}+{r['gene_b']}(sd={r['cross_best_idfm']:.3f})" for r in ep_fp)
        div_named = [f"{r['gene_a']}+{r['gene_b']}(id{r['mean_matched_id']:.2f},sd{r['cross_best_idfm']:.3f})"
                     for r in sorted(tb, key=lambda r: (r["mean_matched_id"], r["gene_a"], r["gene_b"]))]
        tp_named = sorted(f"{r['gene_a']}+{r['gene_b']}(id{r['mean_matched_id']:.2f},sd{r['cross_best_idfm']:.3f})"
                          for r in tp)
        lost = sorted(base_recov - ev["recovered_mc"])
        gained = sorted(ev["recovered_mc"] - base_recov)
        gstm2 = ctrl_fams(gf, GSTM2_GENES); mage = ctrl_fams(gf, MAGE_GENES)
        gstm2_spared = (gstm2["genes_present"] == gstm2_base["genes_present"]
                        and gstm2["n_families"] == gstm2_base["n_families"])
        mage_spared = (mage["genes_present"] == mage_base["genes_present"]
                       and mage["n_families"] == mage_base["n_families"])
        return dict(
            refine=label,
            n_families=ev["n_families"], impure_blocks=ev["impure_blocks"],
            E_p=ev["E_p"], E_p_delta_vs_rederive_baseline=round(ev["E_p"] - base_ev["E_p"], 4),
            E_p_delta_vs_catalog_0_8918=round(ev["E_p"] - 0.8918, 4),
            R_oracle=f"{ev['R_oracle']}/{n_mc}", R_oracle_held_50_57=(ev["R_oracle"] >= base_R),
            oracle_genes_lost=lost, oracle_genes_gained=gained,
            # (a) domain-share FP removed
            FP_genuine_split=len(fp), FP_ep_impure_split=len(ep_fp),
            named_domain_share_FP_removed=named_fp_removed,
            exemplar_status=exemplar_status,
            ep_impure_FP_removed_named=ep_fp_named,
            # (b) real paralog recall cost, divergent broken out + named
            REAL_split_total=real_split, REAL_TP_split=len(tp),
            REAL_truthbar_divergent_split=len(tb),
            REAL_divergent_meanid_lt_0_90_split=len(div_tb),
            split_FP_fraction_pct=(round(100 * len(fp) / total_split, 1) if total_split else None),
            divergent_real_cut_named=div_named,
            TP_real_cut_named=tp_named,
            # (e) cardinality controls
            GSTM2_spared=gstm2_spared, MAGE_spared=mage_spared, GSTM2=gstm2, MAGE=mage,
        )

    per_floor = {}
    for fl in FLOORS:
        per_floor[f"floor_{int(fl*100)}"] = measure(
            partition_refine(fl, use_cn=False),
            f"keep cross-gene edge iff SEDEF cross-support fracMatch >= {fl:.2f} (SD{int(fl*100)})")

    # famCN (diploid-CN) leg in isolation: cut ONLY CN-inconsistent edges (segdup leg inert) --
    # does the famCN refine touch GSTM2 / MAGE cardinality?
    cn_only = measure(partition_refine(None, use_cn=True),
                      "famCN leg only: cut cross-gene edge iff diploid-CN inconsistent (|asm_hapCN "
                      "diff|>1); segdup leg inert")
    # best floor + famCN combined (the full Soto refine)
    best_fl = 0.98
    combined = measure(partition_refine(best_fl, use_cn=True),
                       f"full Soto refine: SEDEF cross-support >= {best_fl:.2f} AND famCN-consistent")

    # best floor: highest E_p among floors that HOLD R_oracle >= baseline.
    # HONEST NOTE (wording precision): removing the 5 sd=0.000 single-exon domain-share FP edges by
    # itself costs NO oracle genes -- but NO SEDEF floor is recall-safe as a whole. SD90's 49/57
    # (one oracle gene lost, LOC101127437) comes from DIVERGENT-REAL collateral, not from the clean
    # single-exon FP removals. "removes the clean FP at zero oracle cost" is true ONLY of those 5
    # edges, NOT of SD90 as an operating point.
    holding = [fl for fl in FLOORS if per_floor[f"floor_{int(fl*100)}"]["R_oracle_held_50_57"]]
    if holding:
        best_hold = max(holding, key=lambda fl: per_floor[f"floor_{int(fl*100)}"]["E_p"])
        best_note = (f"floor {best_hold:.2f}: highest E_p among SEDEF floors that HOLD R_oracle "
                     f">= baseline ({base_R}/{n_mc}). NOTE R_oracle measures the 57 near-identical "
                     "DNA-confirmed multi-copy oracle genes -- it can HOLD while the segdup refine "
                     "still splits divergent truthbar reals (that recall cost is REAL_truthbar_"
                     "divergent_split, measured separately).")
    else:
        best_hold = None
        best_note = f"no SEDEF floor holds R_oracle >= baseline ({base_R}/{n_mc})"

    # famCN cardinality verdict on GSTM2 / MAGE.
    # HONEST NOTE (counting precision): the famCN leg flags 6 CN-inconsistent gene-pairs, NOT 5 GSTM2:
    #   4 GSTM2-anchored (CN19 vs 3/13/16/4: LOC101129940/LOC115930164/LOC115930576/LOC115935025)
    #   + 1 GSTM-cluster TP FALSE-FLAG (LOC115930164+LOC115930576, CN13 vs 16 -- MAD<1 too strict for
    #     a large amplicon) = 5 GSTM-family pairs, PLUS FOXO1+LOC115933254 (CN2 vs 8) = 6 total.
    cn_pairs_inconsistent = sorted(
        f"{'+'.join(sorted(gp))}(cnA={cn_by_gp[gp][2]},cnB={cn_by_gp[gp][3]},cls={sub[gp]['cls']})"
        for gp in sub if cn_by_gp.get(gp, (0,))[0] == 2)
    gstm2_cn_effect = ("SPLIT" if not cn_only["GSTM2_spared"] else "unchanged")
    mage_cn_effect = ("SPLIT" if not cn_only["MAGE_spared"] else "unchanged")

    return dict(
        operationalization=("locus subgraph per catalog family + same-gene backbone (guard); KEEP "
                            "cross-gene edge iff gene-pair SEDEF-segdup-supported (cross_best_fracMatch"
                            " >= floor); re-derive components; drop <2 distinct loci; re-eval via "
                            "shipped eval_partition + oracle_residuals + gene-projection relabel. "
                            "DNA-INFORMED (SEDEF + diploid CN), NOT the RNA-only edge rule."),
        catalog_reference=dict(n_families=573, E_p=0.8918, R_oracle="50/57"),
        rederive_baseline_no_cut=dict(
            n_families=base_ev["n_families"], E_p=base_ev["E_p"], R_oracle=f"{base_R}/{n_mc}",
            impure_blocks=base_ev["impure_blocks"],
            note="same re-derivation baseline as the colinear gate; ALL refine deltas measured vs "
                 "THIS baseline."),
        substrate_pairs_comembered_in_baseline=n_base_co,
        per_floor=per_floor,
        famCN_leg_only=cn_only,
        full_soto_refine_sd98_plus_famCN=combined,
        best_floor=best_hold, best_floor_note=best_note,
        famCN_cardinality=dict(
            note="famCN (diploid-CN) refine is DNA-reachable only where BOTH genes resolve to the "
                 "51-gene multi-copy CN oracle; it adjudicates cardinality within amplicons, not the "
                 "domain-share residual (single-copy partners are absent from the oracle).",
            n_gene_pairs_CN_inconsistent=len(cn_pairs_inconsistent),
            CN_inconsistent_gene_pairs=cn_pairs_inconsistent,
            GSTM2_effect=gstm2_cn_effect, MAGE_effect=mage_cn_effect,
            GSTM2_baseline=gstm2_base, GSTM2_after_famCN=cn_only["GSTM2"],
            MAGE_baseline=mage_base, MAGE_after_famCN=cn_only["MAGE"]),
    )


def main():
    assert SEDEF, "final.bed (SEDEF) not found in " + str(SEDEF_CANDIDATES)
    print(f"[load] SEDEF segdup map {SEDEF} ...", flush=True)
    by, starts, n_pairs = load_sedef(SEDEF)
    print(f"  {sum(len(v) for v in by.values())} segdup intervals ({n_pairs} pairs) over "
          f"{len(by)} contigs", flush=True)
    meta = load_meta(META_TSV)
    cn = load_cn(CN_TSV)
    print(f"[load] {len(meta)} loci coords ; {len(cn)} CN-oracle genes", flush=True)

    # pair substrate (reuse colinear_multiexon_gate.tsv verbatim)
    pairs = []
    with open(PAIRS_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pairs.append(row)

    out_rows = []
    n_missing = 0
    for p in pairs:
        da, db = p["dn_a"], p["dn_b"]
        if da not in meta or db not in meta:
            n_missing += 1
            continue
        locA = meta[da]; locB = meta[db]
        sup = segdup_support(by, starts, locA, locB)
        cn_status, cn_a, cn_b, cn_diff = cn_consistency(cn, p["gene_a"], p["gene_b"])
        rec = dict(
            gene_a=p["gene_a"], gene_b=p["gene_b"], dn_a=da, dn_b=db,
            cls=p["cls"], ep_impure=int(p["ep_impure"]),
            gstm2=int(p["gstm2"]), mage=int(p["mage"]),
            in_dna_loose=int(p["in_dna_loose"]),
            mean_matched_id=float(p["mean_matched_id"]),
            div_bin=p["div_bin"], colinear=int(p["colinear"]),
            chromA=locA[0], startA=locA[1], endA=locA[2],
            chromB=locB[0], startB=locB[1], endB=locB[2],
            same_chrom=int(locA[0] == locB[0]),
            **sup,
            cn_status=cn_status,
            cn_a=("" if cn_a is None else cn_a), cn_b=("" if cn_b is None else cn_b),
            cn_diff=("" if cn_diff is None else cn_diff),
        )
        # Soto-style keep decisions at each SEDEF floor (keep edge iff cross-support >= floor)
        for fl in FLOORS:
            rec[f"segdup_keep_{int(fl*100)}"] = int(sup["cross_best_idfm"] >= fl)
        out_rows.append(rec)

    # ---------------------------------------------------------------- per-pair TSV
    cols = (["gene_a", "gene_b", "dn_a", "dn_b", "cls", "ep_impure", "gstm2", "mage",
             "in_dna_loose", "mean_matched_id", "div_bin", "colinear",
             "chromA", "startA", "endA", "chromB", "startB", "endB", "same_chrom",
             "cross_n", "cross_best_idfm", "cross_best_idfi", "same_region", "same_best_idfm",
             "a_in_segdup", "b_in_segdup", "cn_status", "cn_a", "cn_b", "cn_diff",
             "segdup_keep_90", "segdup_keep_95", "segdup_keep_98"])
    out_rows.sort(key=lambda r: (r["cls"], r["gene_a"], r["gene_b"], r["dn_a"], r["dn_b"]))
    with open(OUT_TSV, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in out_rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---------------------------------------------------------------- summaries
    by_cls = collections.defaultdict(list)
    for r in out_rows:
        by_cls[r["cls"]].append(r)

    # (target) DOMAIN-SHARE FP = ep_impure genuine  (+ colinear==1 subset)
    ds_fp = [r for r in by_cls["genuine"] if r["ep_impure"]]
    ds_fp_col1 = [r for r in ds_fp if r["colinear"] == 1]
    # broader over-merge FP = all genuine
    genuine = by_cls["genuine"]

    # REAL paralogs
    tp = by_cls["TP"]                       # near-identical (DNA-confirmed)
    tb = by_cls["truthbar"]                 # DIVERGENT (E_p-confirmed, DNA did not confirm)
    real = tp + tb

    # ---- (a) does segdup-support SEPARATE domain-share FP (OUT) from real (IN)? ----
    separation = dict(
        DOMAIN_SHARE_FP_ep_impure_genuine=capture_block(ds_fp, None),
        DOMAIN_SHARE_FP_ep_impure_genuine_colinear1=capture_block(ds_fp_col1, None),
        ALL_over_merge_FP_genuine=capture_block(genuine, None),
        REAL_near_identical_TP=capture_block(tp, None),
        REAL_divergent_truthbar=capture_block(tb, None),
        REAL_all=capture_block(real, None),
    )

    # ---- (b) DECISIVE: divergent real (truthbar) captured in SEDEF, by identity floor + RNA div ----
    tb_by_rnabin = collections.defaultdict(list)
    for r in tb:
        tb_by_rnabin[r["div_bin"]].append(r)
    rna_order = ["id>=0.95", "0.90-0.95", "0.80-0.90", "0.70-0.80", "<0.70"]
    divergent_stratified_by_rna = {}
    for b in rna_order:
        rs = tb_by_rnabin.get(b, [])
        if rs:
            divergent_stratified_by_rna[b] = capture_block(rs, None)

    # also stratify truthbar by their OWN best SEDEF identity (where the DNA erosion floor bites)
    def sedef_id_bin(r):
        if r["cross_n"] == 0:
            return "no_cross_pair"
        v = r["cross_best_idfm"]
        if v >= 0.98:
            return "sd>=0.98"
        if v >= 0.95:
            return "sd0.95-0.98"
        if v >= 0.90:
            return "sd0.90-0.95"
        return "sd<0.90"
    tb_sedef_id_hist = collections.Counter(sedef_id_bin(r) for r in tb)
    tp_sedef_id_hist = collections.Counter(sedef_id_bin(r) for r in tp)
    fp_sedef_id_hist = collections.Counter(sedef_id_bin(r) for r in ds_fp)

    # ---- Soto refine effect: FP removed vs real cut, at each floor (pair-level) ----
    refine_by_floor = {}
    for fl in FLOORS:
        kk = f"segdup_keep_{int(fl*100)}"
        def cut(rows):
            c = sum(1 for r in rows if not r[kk])
            return dict(n=len(rows), cut=c, cut_pct=pct(c, len(rows)),
                        kept=len(rows) - c, kept_pct=pct(len(rows) - c, len(rows)))
        refine_by_floor[f"floor_{int(fl*100)}"] = dict(
            DOMAIN_SHARE_FP_removed=cut(ds_fp),          # want HIGH cut
            ALL_genuine_FP_removed=cut(genuine),
            REAL_TP_cut=cut(tp),                         # want LOW cut
            REAL_divergent_truthbar_cut=cut(tb),         # the recall cost (want LOW, but expect high)
            REAL_all_cut=cut(real),
        )

    # ---- CN-consistency leg (sparse: only multi-copy CN-oracle genes) ----
    cn_status_hist = {cls: dict(collections.Counter(r["cn_status"] for r in rs))
                      for cls, rs in by_cls.items()}
    cn_both = [r for r in out_rows if r["cn_status"] in ("consistent", "inconsistent")]
    cn_leg = dict(
        note="famCN refine only adjudicates pairs where BOTH genes are in the multi-copy CN oracle; "
             "domain-share partners are typically single-copy genes ABSENT from the oracle (that "
             "absence is itself consistent with 'not a duplicated family').",
        n_pairs_with_both_genes_in_oracle=len(cn_both),
        by_class=cn_status_hist,
        consistent_pairs=sorted(f"{r['gene_a']}+{r['gene_b']}(cnA={r['cn_a']},cnB={r['cn_b']},"
                                f"cls={r['cls']})" for r in cn_both if r["cn_status"] == "consistent"),
        inconsistent_pairs=sorted(f"{r['gene_a']}+{r['gene_b']}(cnA={r['cn_a']},cnB={r['cn_b']},"
                                  f"cls={r['cls']})" for r in cn_both if r["cn_status"] == "inconsistent"),
    )

    # ---- named domain-share exemplars ----
    by_gp = collections.defaultdict(list)
    for r in out_rows:
        by_gp[frozenset((r["gene_a"], r["gene_b"]))].append(r)
    exemplar_out = []
    for ga, gb in DOMAIN_SHARE_EXEMPLARS:
        rs = by_gp.get(frozenset((ga, gb)), [])
        if rs:
            r = max(rs, key=lambda x: x["cross_best_idfm"])   # strongest-support representative
            exemplar_out.append(dict(
                pair=f"{ga}+{gb}", cls=r["cls"], ep_impure=r["ep_impure"], colinear=r["colinear"],
                same_chrom=r["same_chrom"], cross_n=r["cross_n"],
                cross_best_idfm=r["cross_best_idfm"], both_in_segdup=int(r["a_in_segdup"] and r["b_in_segdup"]),
                segdup_supported_90=r["segdup_keep_90"], cn_status=r["cn_status"]))
        else:
            exemplar_out.append(dict(pair=f"{ga}+{gb}", note="no direct multi-exon E_r edge in substrate"))

    # ---- GSTM2 / MAGE controls (cardinality: test the CN-consistency leg) ----
    def ctrl(flag):
        rs = [r for r in out_rows if r[flag]]
        return dict(
            n=len(rs),
            cls_mix=dict(collections.Counter(r["cls"] for r in rs)),
            cross_support_any_pct=pct(sum(1 for r in rs if r["cross_n"] >= 1), len(rs)),
            cross_support_ge_98_pct=pct(sum(1 for r in rs if r["cross_best_idfm"] >= 0.98), len(rs)),
            both_in_segdup_pct=pct(sum(1 for r in rs if r["a_in_segdup"] and r["b_in_segdup"]), len(rs)),
            cn_status_mix=dict(collections.Counter(r["cn_status"] for r in rs)),
            cn_consistent_pct=pct(sum(1 for r in rs if r["cn_status"] == "consistent"), len(rs)),
        )
    controls = dict(GSTM2=ctrl("gstm2"), MAGE=ctrl("mage"))

    summary = dict(
        substrate=dict(sedef_map=SEDEF, n_sedef_pairs=n_pairs, n_pair_substrate=len(out_rows),
                       n_pairs_dropped_no_coords=n_missing,
                       class_counts=dict(collections.Counter(r["cls"] for r in out_rows)),
                       note="DNA-INFORMED Soto-style refine (SEDEF segdup + diploid CN), NOT the "
                            "RNA-only edge rule. cross-support = the two loci are the two ends of one "
                            "SEDEF segdup pair (locus A dup to locus B), identity = fracMatch."),
        SEPARATION_segdup_support=separation,
        DECISIVE_divergent_real_capture=dict(
            definition="truthbar = E_p-protein-confirmed DIVERGENT real paralogs DNA did not confirm; "
                       "cross_support_ge_XX_pct = fraction whose loci are a SEDEF segdup pair at "
                       ">= that identity (the DNA recall the SD98 restriction retains).",
            truthbar_overall=capture_block(tb, None),
            truthbar_by_RNA_divergence_bin=divergent_stratified_by_rna,
            truthbar_SEDEF_identity_hist=dict(tb_sedef_id_hist),
            TP_SEDEF_identity_hist=dict(tp_sedef_id_hist),
            DOMAIN_SHARE_FP_SEDEF_identity_hist=dict(fp_sedef_id_hist),
        ),
        soto_refine_effect_by_floor=refine_by_floor,
        cn_consistency_leg=cn_leg,
        named_domain_share_exemplars=exemplar_out,
        controls_GSTM2_MAGE=controls,
        params=dict(FLOORS=FLOORS, identity="fracMatch (col21) primary; fracMatchIndel (col22) carried",
                    PYTHONHASHSEED=os.environ.get("PYTHONHASHSEED")),
        out_tsv=OUT_TSV,
    )

    # ---------------- REFINE STAGE: apply Soto-style segdup(+famCN) refine to the catalog ----------
    print("[refine] applying Soto-style SEDEF-segdup refine to the shipped catalog "
          "(floors 0.90/0.95/0.98 + famCN leg) ...", flush=True)
    summary["refine_stage"] = refine_stage(out_rows)

    with open(OUT_JSON, "w") as jf:
        json.dump(summary, jf, indent=1, sort_keys=True)
        jf.write("\n")
    print(json.dumps(summary, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
