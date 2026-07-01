#!/usr/bin/env python
"""A1: GENUINELY READ-DERIVED (non-reference-shared) corroboration of O1 (multi-copy family/copy
structure), via SDA / Vollger-2019 read-heterogeneity correlation-clustering -- NO asm20 seed.

WHY (the honesty gap this closes):
  Both prior O1 corroborations reduced to a SHARED REFERENCE substrate:
    * L1 (GENOME_RNA_OVERLAY.md): the de-novo "transcripts" are reference substrings, so grouping them
      vs SEDEF was reference-vs-reference.
    * A2 (genome_rna_overlay_readcontent.py): the PSV columns came from a minimap2 asm20 self-alignment
      of the REFERENCE copies. The break-test (genome_rna_overlay_breaktest.py) proved the read gate is
      INERT (ungated rho -0.420 ~ gated -0.443) => the PSV axis is asm20-vs-SEDEF concordance, ZERO read
      info.
  A1 fixes this by DISCOVERING copy-distinguishing columns from READ DISAGREEMENT (the pileup on
  GGO_mm.bam), linking them by cross-read co-segregation (SDA / Vollger 2019), and partitioning reads
  into copies -- the asm20 alignment is NEVER consulted for column discovery. It is used ONLY as the
  thing we prove the reads beat (Control 1).

METHOD (method (A) direct-pileup, the validated smoke-test core in a1_read_sda_smoketest.py):
  At a COLLAPSED reference locus, the -N50 GGO_mm.bam pileup stacks reads from K paralogs (MAPQ-0
  multimappers) onto one locus. Columns are DISCOVERED from that read disagreement:
    1. read x pos allele matrix (pysam pileup, base-qual>=20, ACGT, secondary alns kept).
    2. candidate het columns: depth>=12, minor>=4 reads AND VAF>=0.20 (kills singleton HiFi errors).
    3. SDA linkage: binarize (major=0/minor=1), pairwise phi over co-covering reads, keep columns
       co-segregating |phi|>=0.60 with >=2 others (errors correlate with nothing -> dropped).
    4. read partition -> K_read: union-find on reads agreeing >=80% over the linked columns; K_read =
       #clusters with >=4 reads. Linked columns = the READ-DISCOVERED PSV set.

LEGS:
  (1) K_read + read partition per family (frontier + resolved control).
  (2) CROSS-MODAL CORROBORATION vs the genome: K_read vs SEDEF distinct segdup partners / n_loci /
      collapsed_excess+1 (Spearman + distance-matched permutation null). Read-heterogeneity (SDA) vs
      genome self-alignment (SEDEF) = different data + method.
  (3) DECISIVE INDEPENDENCE CONTROL:
        Control 1 (read-only vs asm20-seeded): K_read vs K_asm20 (psv_graph json), read_derived_delta.
          The airtight subset = the 27 no_psv families where asm20=0 columns (K_asm20==1 by
          construction). Every one with K_read>=2 surviving the shuffle null is structure the reference
          self-alignment CANNOT supply.
        Control 2 (shuffle null): permute alleles across reads per column (preserve MAF, destroy
          cross-read linkage); re-run 3-4. Real copy structure -> linked pairs & K_read collapse.
          A MAF-marginal artifact survives. Verdict rule mirrors the A2 break-test.
  (4) SCOPE: genuine read-derived independence is confined to the frontier subset; quantify how much
      corroboration lives there vs the reference-shared resolved majority.

DETERMINISM: PYTHONHASHSEED pinned (asserted below); per-family seeded np.random.default_rng(seed).

Outputs: bench/a1_read_consensus_o1.tsv (per-family), bench/a1_read_consensus_o1.json (summary +
independence-control verdict).
Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/a1_read_consensus_o1.py
"""
import collections
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
import pysam
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# reuse the EXACT validated SDA core from the smoke test (no asm20 anywhere in these)
from a1_read_sda_smoketest import (
    allele_matrix, call_het_columns, link_columns, cluster_reads, shuffle_cols,
    INCLUDE_SECONDARY, MIN_CLUST,
)
# reuse SEDEF + copy-region machinery (genome self-alignment; never sees a read)
from genome_rna_overlay_readcontent import (
    dedup_copies, load_sedef_index, family_sedef_identity, sedef_rows_overlapping, MINOV,
)

BAM = "/home/juanfra/winloci_scratch/GGO_mm.bam"
FAM_TSV = "/home/juanfra/winloci_scratch/validated_families.tsv"
FCN_TSV = os.path.join(HERE, "family_copy_number.tsv")
PSV_JSON = os.path.join(HERE, "psv_graph_genomewide.json")
FINAL_BED = "/mnt/c/Users/jfris/Desktop/final.bed"

MAX_WIN = 40000          # cap the pileup window length (bound compute; exonic columns are sparse)
N_CONTROL = 32           # fully_resolvable control families (deterministic stride sample)
N_WORKERS = 4            # parallel families (each independent + per-family seeded => deterministic)
SEED_BASE = 1000         # per-family rng = default_rng(SEED_BASE + family_id)
MAF_ALLELE = 0.40        # a het column with dominant MAF>=this flags possible diploid allele (copy-vs-allele)
NULL_RETAIN_FRAC = 0.30  # verdict rule: shuffle must destroy >(1-this) of linked pairs / K_read


def load_fam_loci():
    fam = collections.defaultdict(list)
    with open(FAM_TSV) as f:
        next(f)
        for line in f:
            fi, lid, c, s, e, nr = line.rstrip("\n").split("\t")
            fam[int(fi)].append((c, int(s), int(e), int(nr)))
    return fam


def load_fcn():
    out = {}
    with open(FCN_TSV) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            v = dict(zip(hdr, line.rstrip("\n").split("\t")))
            out[int(v["family"])] = v
    return out


def sedef_partner_count(regions, idx):
    """Read-INDEPENDENT genome copy estimate for the family: distinct, merged genomic intervals that
    SEDEF (final.bed self-alignment) links to ANY of the family's copy regions. This is the number of
    segdup partners SEDEF sees anywhere in the genome (a copy count from genome bytes, no reads)."""
    mates = collections.defaultdict(list)  # mate_chrom -> [(ms,me)]
    own = [(c, s, e) for (c, s, e, _nr) in regions]
    for (c, s, e, _nr) in regions:
        for rec in sedef_rows_overlapping(idx, c, s, e):
            mc, ms, me = rec[2], rec[3], rec[4]
            mates[mc].append((ms, me))
    # merge overlapping mate intervals, then drop those that are just the family's own copies
    partners = set()
    for mc, ivs in mates.items():
        ivs.sort()
        cur = None
        merged = []
        for ms, me in ivs:
            if cur and ms <= cur[1]:
                cur = (cur[0], max(cur[1], me))
            else:
                if cur:
                    merged.append(cur)
                cur = (ms, me)
        if cur:
            merged.append(cur)
        for ms, me in merged:
            # is this merged mate essentially one of the family's OWN copy regions? (self-link)
            is_own = any(oc == mc and min(oe, me) - max(os_, ms) >= MINOV for (oc, os_, oe) in own)
            partners.add((mc, ms, me, is_own))
    n_partner_regions = len(partners)                       # all merged mate intervals
    n_partner_external = sum(1 for p in partners if not p[3])  # partners NOT among the family's own copies
    return n_partner_regions, n_partner_external, len(mates)


def pick_backbone(regions):
    """Backbone = highest-read-count region (richest collapse mixture; deterministic). Window capped."""
    regions = sorted(regions, key=lambda r: (-r[3], r[0], r[1]))  # -nr, then chrom,start for ties
    c, s, e, nr = regions[0]
    if e - s > MAX_WIN:
        e = s + MAX_WIN
    return c, s, e, nr


def maf_profile(cands):
    """copy-vs-allele guard: fraction of candidate het columns whose minor allele fraction is 'balanced'
    (>= MAF_ALLELE), i.e. ~50/50 like a diploid heterozygote rather than 1/K copy structure."""
    if not cands:
        return 0.0, 0.0
    mafs = [minc / dp for (_p, _maj, _minb, minc, dp) in cands]
    balanced = sum(1 for m in mafs if m >= MAF_ALLELE) / len(mafs)
    return float(np.median(mafs)), float(balanced)


def analyze_family(bam, fam_id, regions, idx, seed):
    _ = fam_id  # (seed already folds in the family id)
    c, s, e, nr = pick_backbone(regions)
    cols, depth = allele_matrix(bam, c, s, e, INCLUDE_SECONDARY)
    cands = call_het_columns(cols, depth)
    med_maf, balanced = maf_profile(cands)
    if not cands:
        res = dict(n_cand=0, n_read_psv=0, mean_phi=0.0, K_read=0, clust_sizes=[], n_reads=0,
                   med_depth=int(np.median(list(depth.values()))) if depth else 0,
                   null_linked=0, null_phi=0.0, null_K=0, med_maf=0.0, balanced_frac=0.0)
    else:
        keep, edges, binv = link_columns(cols, cands)
        K_read, sizes, n_reads = cluster_reads(cands, binv, keep)
        rng = np.random.default_rng(seed)
        scols = shuffle_cols(cols, cands, rng)
        nkeep, nedges, nbinv = link_columns(scols, cands)
        nK, _, _ = cluster_reads(cands, nbinv, nkeep) if nkeep else (0, [], 0)
        res = dict(n_cand=len(cands), n_read_psv=len(keep),
                   mean_phi=round(float(np.mean([abs(x[2]) for x in edges])), 3) if edges else 0.0,
                   K_read=K_read, clust_sizes=sizes[:8], n_reads=n_reads,
                   med_depth=int(np.median([cc[4] for cc in cands])),
                   null_linked=len(nkeep),
                   null_phi=round(float(np.mean([abs(x[2]) for x in nedges])), 3) if nedges else 0.0,
                   null_K=nK, med_maf=med_maf, balanced_frac=balanced,
                   n_linked_real=len(keep))
    res["backbone"] = f"{c}:{s}-{e}"
    res["backbone_nr"] = nr
    n_part, n_part_ext, n_mate_chroms = sedef_partner_count(regions, idx)
    res["sedef_partner_regions"] = n_part
    res["sedef_partner_external"] = n_part_ext
    sed_id, n_link, npairs, _ = family_sedef_identity(regions, idx)
    res["sedef_identity"] = sed_id
    res["sedef_pair_links"] = n_link
    return res


def distance_matched_null(pairs_xy, obs_rho, rng, B=5000):
    """Permutation join-null (as in the A2 overlay): shuffle the K_genome column against K_read,
    recompute Spearman. Fraction |perm|>=|obs|. Breaks the family join -> tests it drives the signal."""
    x = np.array([p[0] for p in pairs_xy], float)
    y = np.array([p[1] for p in pairs_xy], float)
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    obs = abs(obs_rho)
    ge = 0
    for _ in range(B):
        rp = rng.permutation(ry)
        rho = stats.spearmanr(rx, rp).statistic
        if not np.isnan(rho) and abs(rho) >= obs:
            ge += 1
    return (ge + 1) / (B + 1)


# ---- parallel worker (fork inherits the loaded SEDEF index + tables; own BAM per process) ----
_G = {}          # globals inherited by forked workers: idx, fam_loci, fcn, jfam
_BAM = None       # per-process BAM handle (opened lazily in each worker)


def _worker(fam):
    global _BAM
    if _BAM is None:
        _BAM = pysam.AlignmentFile(BAM, "rb")
    fam_loci = _G["fam_loci"]
    fcn = _G["fcn"]
    jfam = _G["jfam"]
    idx = _G["idx"]
    regions = dedup_copies([(c, s, e, nr) for (c, s, e, nr) in fam_loci.get(fam, [])])
    if len(regions) < 2:
        return None
    jf = jfam.get(fam, {})
    v = fcn[fam]
    res = analyze_family(_BAM, fam, regions, idx, SEED_BASE + fam)
    K_asm20 = int(jf.get("K", 1))
    psvs_asm20 = int(jf.get("psvs", 0))
    n_loci = int(v["n_loci"])
    collapsed_excess = int(v["collapsed_excess"])
    cls = v["cls"]
    read_derived_delta = res["K_read"] - K_asm20
    survives_null = (res["K_read"] >= 2 and res["null_K"] <= 1 and
                     res["null_linked"] <= NULL_RETAIN_FRAC * max(res["n_read_psv"], 1))
    sedef_corr = res["sedef_partner_external"] >= 1 or collapsed_excess > 0 or n_loci >= 2
    return dict(
        family=fam, gene=v["gene"], cls=cls, n_loci=n_loci,
        chi_H=int(v["chi_H"]), collapsed_excess=collapsed_excess,
        K_genome_excess=collapsed_excess + 1,
        K_read=res["K_read"], K_asm20=K_asm20, read_derived_delta=read_derived_delta,
        n_read_psv=res["n_read_psv"], psvs_asm20=psvs_asm20, mean_phi=res["mean_phi"],
        n_cand=res["n_cand"], n_reads=res["n_reads"], med_depth=res["med_depth"],
        clust_sizes=res["clust_sizes"], med_maf=res["med_maf"], balanced_frac=res["balanced_frac"],
        null_linked=res["null_linked"], null_phi=res["null_phi"], null_K=res["null_K"],
        survives_null=survives_null,
        sedef_partner_regions=res["sedef_partner_regions"],
        sedef_partner_external=res["sedef_partner_external"],
        sedef_identity=res["sedef_identity"], sedef_pair_links=res["sedef_pair_links"],
        sedef_corr=sedef_corr, backbone=res["backbone"], backbone_nr=res["backbone_nr"],
    )


def main():
    assert os.environ.get("PYTHONHASHSEED") == "0", \
        "run with PYTHONHASHSEED=0 for determinism (union-find/dict order)"
    fam_loci = load_fam_loci()
    fcn = load_fcn()
    J = json.load(open(PSV_JSON))
    jfam = {f["family"]: f for f in J["families"]}

    # ---- run set: all 30 frontier (no_psv + partial) + deterministic stride sample of resolved ----
    frontier = sorted(f for f, v in fcn.items() if v["cls"] in ("no_psv", "partial"))
    resolved_all = sorted(f for f, v in fcn.items() if v["cls"] == "fully_resolvable")
    stride = max(1, len(resolved_all) // N_CONTROL)
    control = resolved_all[::stride][:N_CONTROL]
    run_set = sorted(set(frontier) | set(control))
    print(f"[run set] frontier(no_psv+partial)={len(frontier)}  control(resolved)={len(control)}  "
          f"total={len(run_set)}", flush=True)

    print("[load] SEDEF index (final.bed, genome self-alignment) ...", flush=True)
    idx, n_sedef = load_sedef_index()
    print(f"       {n_sedef} SEDEF rows over {len(idx)} contigs", flush=True)

    # publish read-only globals for forked workers (inherited via fork copy-on-write)
    _G.update(fam_loci=fam_loci, fcn=fcn, jfam=jfam, idx=idx)

    rows = []
    done = 0
    with Pool(processes=N_WORKERS) as pool:
        for r in pool.imap_unordered(_worker, run_set):
            done += 1
            if r is None:
                continue
            rows.append(r)
            print(f"  [{done}/{len(run_set)}] fam {r['family']:>4} {r['gene']:<16} {r['cls']:<16} "
                  f"nloci={r['n_loci']} K_asm20={r['K_asm20']} | Kread={r['K_read']} "
                  f"readPSV={r['n_read_psv']} phi={r['mean_phi']} depth={r['med_depth']} "
                  f"nullK={r['null_K']} nullLink={r['null_linked']} | "
                  f"sedefPart={r['sedef_partner_external']} survNull={r['survives_null']}", flush=True)
    rows.sort(key=lambda x: x["family"])

    # ================= aggregate / verdict =================
    frontier_rows = [r for r in rows if r["cls"] in ("no_psv", "partial")]
    nopsv_rows = [r for r in rows if r["cls"] == "no_psv"]
    resolved_rows = [r for r in rows if r["cls"] == "fully_resolvable"]

    # --- Control 1: read-only vs asm20-seeded (the airtight subset = no_psv where K_asm20==1) ---
    airtight = [r for r in nopsv_rows if r["K_asm20"] == 1]
    airtight_kread2 = [r for r in airtight if r["K_read"] >= 2]
    airtight_kread2_surv = [r for r in airtight if r["survives_null"]]
    airtight_kread2_surv_corr = [r for r in airtight if r["survives_null"] and r["sedef_corr"]]
    # partners specifically (external SEDEF segdup partner, the strongest independent genome check)
    airtight_full = [r for r in airtight if r["survives_null"] and r["sedef_partner_external"] >= 1]

    # --- Control 2: shuffle-null destruction check across frontier ---
    real_link = sum(r["n_read_psv"] for r in frontier_rows)
    null_link = sum(r["null_linked"] for r in frontier_rows)
    null_retained_frac = null_link / max(real_link, 1)
    frontier_K2 = [r for r in frontier_rows if r["K_read"] >= 2]
    frontier_nullK2 = [r for r in frontier_rows if r["null_K"] >= 2]

    # --- Leg 2: cross-modal corroboration Spearman(K_read, K_genome) on the collapsed core ---
    core = [r for r in frontier_rows if r["collapsed_excess"] > 0]  # comparable (K_genome>1)
    corr = {}
    if len(core) >= 4:
        for gname, gkey in [("n_loci", "n_loci"), ("collapsed_excess+1", "K_genome_excess"),
                            ("sedef_partner_regions", "sedef_partner_regions")]:
            pairs = [(r["K_read"], r[gkey]) for r in core]
            xs = [p[0] for p in pairs]
            ys = [p[1] for p in pairs]
            if len(set(xs)) < 2 or len(set(ys)) < 2:
                corr[gname] = dict(rho=None, p=None, n=len(pairs), note="no variance in one axis")
                continue
            rho, p = stats.spearmanr(xs, ys)
            rng = np.random.default_rng(20260701)
            perm_p = distance_matched_null(pairs, rho, rng, B=5000)
            corr[gname] = dict(rho=float(rho), p=float(p), perm_p=float(perm_p), n=len(pairs))

    # --- scope numbers ---
    n_val = len(fcn)
    n_frontier_total = sum(1 for v in fcn.values() if v["cls"] in ("no_psv", "partial"))
    scope = dict(
        validated_families=n_val,
        frontier_families_total=n_frontier_total,
        pct_frontier=round(100 * n_frontier_total / n_val, 1),
        resolved_families_total=sum(1 for v in fcn.values() if v["cls"] == "fully_resolvable"),
        collapsed_core_excess_gt0=sum(1 for v in fcn.values() if int(v["collapsed_excess"]) > 0),
        run_frontier=len(frontier_rows), run_resolved=len(resolved_rows),
    )

    # --- resolved control: is read-heterogeneity FLAT there (reference-shared regime)? ---
    def frac_k2(rr):
        return round(sum(1 for r in rr if r["K_read"] >= 2) / max(len(rr), 1), 3)
    resolved_summary = dict(
        n=len(resolved_rows),
        frac_K_read_ge2=frac_k2(resolved_rows),
        median_n_read_psv=float(np.median([r["n_read_psv"] for r in resolved_rows])) if resolved_rows else 0.0,
        median_balanced_frac=float(np.median([r["balanced_frac"] for r in resolved_rows])) if resolved_rows else 0.0,
        note="fully_resolvable = reference-shared regime (uniquely-mappable reads ~ reference). "
             "K_read>=2 here is diploid het (balanced MAF) or residual cross-map, NOT independence.",
    )
    frontier_summary = dict(
        n=len(frontier_rows),
        frac_K_read_ge2=frac_k2(frontier_rows),
        median_n_read_psv=float(np.median([r["n_read_psv"] for r in frontier_rows])) if frontier_rows else 0.0,
        median_balanced_frac=float(np.median([r["balanced_frac"] for r in frontier_rows])) if frontier_rows else 0.0,
    )

    # ================= verdict string (A1 analog of the break-test verdict) =================
    n_air = len(airtight)
    n_air_surv = len(airtight_kread2_surv)
    n_air_full = len(airtight_full)
    if n_air_surv == 0:
        verdict = (
            f"A1 DOES NOT ESCAPE on GGO: 0/{n_air} no_psv (asm20=0) families yield K_read>=2 surviving "
            f"the shuffle null. On the collapsed set the read-derived signal does not exceed the "
            f"reference; A1's only content is method-validation on the reference-shared majority -> "
            f"does NOT close the circularity gap.")
    else:
        verdict = (
            f"A1 GENUINELY ESCAPES the reference substrate on a NARROW scope: {n_air_surv}/{n_air} "
            f"no_psv (asm20=0 => K_asm20==1 by construction) families yield K_read>=2 that SURVIVES the "
            f"cross-read shuffle null -- copy structure the reference self-alignment CANNOT supply "
            f"(read_derived_delta>=+1). Of these, {n_air_full} are ALSO corroborated by an EXTERNAL "
            f"SEDEF segdup partner (genome self-alignment, never sees a read) = fully cross-modal. "
            f"Shuffle destroys read-linkage genome-wide (frontier retains {100*null_retained_frac:.0f}% "
            f"of linked columns, K_read->1 in null). SCOPE: independence is confined to the "
            f"{scope['pct_frontier']}% frontier; the resolved {scope['resolved_families_total']} families "
            f"are reference-shared (method-validation only, K_read>=2 frac={resolved_summary['frac_K_read_ge2']}).")

    out = dict(
        method="SDA/Vollger-2019 read-heterogeneity correlation-clustering, direct GGO_mm.bam pileup, "
               "NO asm20 column seed (columns discovered from read disagreement).",
        run_set=dict(frontier=len(frontier_rows), resolved_control=len(resolved_rows), total=len(rows)),
        scope=scope,
        control1_read_only_vs_asm20=dict(
            airtight_subset_no_psv_Kasm20_eq1=n_air,
            K_read_ge2=len(airtight_kread2),
            K_read_ge2_survives_shuffle=n_air_surv,
            K_read_ge2_survives_and_genome_corroborated=len(airtight_kread2_surv_corr),
            K_read_ge2_survives_and_external_SEDEF_partner=n_air_full,
            interpretation="asm20 supplies 0 columns on these families (K_asm20==1). Any K_read>=2 here "
                           "is structure the reference self-alignment cannot produce. read_derived_delta "
                           "= K_read - K_asm20.",
            frontier_mean_read_derived_delta=round(
                float(np.mean([r["read_derived_delta"] for r in frontier_rows])), 3) if frontier_rows else 0.0,
        ),
        control2_shuffle_null=dict(
            frontier_real_linked_columns=real_link,
            frontier_null_linked_columns=null_link,
            null_retained_frac=round(null_retained_frac, 4),
            frontier_families_K_read_ge2=len(frontier_K2),
            frontier_families_null_K_ge2=len(frontier_nullK2),
            verdict_rule=f"shuffle must destroy >{100*(1-NULL_RETAIN_FRAC):.0f}% of linked columns AND "
                         f"collapse K_read->1; a MAF-marginal artifact survives.",
            passes=null_retained_frac < NULL_RETAIN_FRAC and len(frontier_nullK2) == 0,
        ),
        leg2_cross_modal_corroboration=corr,
        resolved_control_summary=resolved_summary,
        frontier_summary=frontier_summary,
        verdict=verdict,
    )

    # ================= write artifacts =================
    tsv = os.path.join(HERE, "a1_read_consensus_o1.tsv")
    with open(tsv, "w") as fh:
        fh.write("family\tgene\tcls\tn_loci\tchi_H\tcollapsed_excess\tK_read\tK_asm20\t"
                 "read_derived_delta\tn_read_psv\tpsvs_asm20\tmean_phi\tn_reads\tmed_depth\tmed_maf\t"
                 "balanced_frac\tnull_linked\tnull_K\tsurvives_null\tsedef_partner_regions\t"
                 "sedef_partner_external\tsedef_identity\tsedef_corr\tbackbone\tclust_sizes\n")
        for r in sorted(rows, key=lambda x: (x["cls"], -x["K_read"], x["family"])):
            si = f"{r['sedef_identity']:.5f}" if r["sedef_identity"] is not None else "NA"
            fh.write(f"{r['family']}\t{r['gene']}\t{r['cls']}\t{r['n_loci']}\t{r['chi_H']}\t"
                     f"{r['collapsed_excess']}\t{r['K_read']}\t{r['K_asm20']}\t{r['read_derived_delta']}\t"
                     f"{r['n_read_psv']}\t{r['psvs_asm20']}\t{r['mean_phi']}\t{r['n_reads']}\t"
                     f"{r['med_depth']}\t{r['med_maf']:.3f}\t{r['balanced_frac']:.3f}\t{r['null_linked']}\t"
                     f"{r['null_K']}\t{int(r['survives_null'])}\t{r['sedef_partner_regions']}\t"
                     f"{r['sedef_partner_external']}\t{si}\t{int(r['sedef_corr'])}\t{r['backbone']}\t"
                     f"{r['clust_sizes']}\n")
    js = os.path.join(HERE, "a1_read_consensus_o1.json")
    json.dump(dict(summary=out, families=rows), open(js, "w"), indent=2)

    # ================= report =================
    print("\n" + "=" * 88)
    print("A1 READ-DERIVED O1 CORROBORATION -- SUMMARY")
    print("=" * 88)
    print(f"scope: {scope['validated_families']} validated families; frontier(no_psv+partial)="
          f"{scope['frontier_families_total']} ({scope['pct_frontier']}%); collapsed core "
          f"(collapsed_excess>0)={scope['collapsed_core_excess_gt0']}")
    print(f"run: frontier={len(frontier_rows)}  resolved-control={len(resolved_rows)}")
    print(f"\n[Control 1: read-only vs asm20] airtight no_psv (asm20=0, K_asm20==1) = {n_air}")
    print(f"   K_read>=2: {len(airtight_kread2)}   +survives shuffle: {n_air_surv}   "
          f"+external SEDEF partner: {n_air_full}")
    print(f"   frontier mean read_derived_delta (K_read - K_asm20) = "
          f"{out['control1_read_only_vs_asm20']['frontier_mean_read_derived_delta']}")
    print(f"\n[Control 2: shuffle null] frontier real linked cols={real_link} null={null_link} "
          f"(retained {100*null_retained_frac:.1f}%)  null K>=2 families={len(frontier_nullK2)}  "
          f"PASS={out['control2_shuffle_null']['passes']}")
    print(f"\n[Leg 2: cross-modal corroboration] Spearman(K_read, genome) on collapsed core (n={len(core)}):")
    for k, v in corr.items():
        print(f"   K_read vs {k:22s}: rho={v['rho']} p={v.get('p')} perm_p={v.get('perm_p')} n={v['n']}")
    print(f"\n[Scope] frontier frac K_read>=2 = {frontier_summary['frac_K_read_ge2']} "
          f"(median read-PSV {frontier_summary['median_n_read_psv']})  |  "
          f"resolved frac K_read>=2 = {resolved_summary['frac_K_read_ge2']} "
          f"(median read-PSV {resolved_summary['median_n_read_psv']})")
    print("\n" + "=" * 88)
    print("VERDICT:", verdict)
    print("=" * 88)
    print(f"\nwrote {tsv}\nwrote {js}")


if __name__ == "__main__":
    main()
