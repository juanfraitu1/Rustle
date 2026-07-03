#!/usr/bin/env python
"""ABSOLUTE DIVERGENCE FLOOR -- MEASURE STAGE (sweep, no wiring).

WHAT THIS MEASURES
------------------
The shipped RNA family definition (bench/family_rna_refine.py) keeps a cross-gene
edge iff  core_recip>=0.19 AND aln_frac>=0.24  (+3 VG gates + gamma + demote). It has
NO identity floor, so it admits DIVERGENT paralogs and single-exon DOMAIN-SHARES whose
whole-transcript RECIPROCAL identity is low. This stage adds an ABSOLUTE RECIPROCAL-
IDENTITY FLOOR to the edge rule and SWEEPS it, to scope the family definition to gene
families whose members are >= FLOOR% similar (Soto's SD98 near-identical move, but as an
opt-out-able default).

THE METRIC (RNA-only; NOT DNA)
------------------------------
Per cross-gene edge, from the SAME spliced-exon minimap2 alignment that defines aln_frac
(bench/ri_sharedlen_recip_id.tsv; built by bench/ri_build_recip_id.py, aln_len reproduced
byte-for-byte vs ri_sharedlen_universal.tsv):
    recip_id_best = matches_best / max(len_a, len_b)
                  = min( matches/len_a , matches/len_b )
the RECIPROCAL whole-transcript identity = the MIN over the two members of
(aligned-AND-matching bases / member length). It is reciprocal: a small gene embedded in
one big exon of a large gene scores LOW on the large gene's side (max(len) big) -> cut.
aln_frac (block/min) does NOT have this property (a domain-share shares one exon at ~95%).

Class medians on the 1457 core+aln-passing edges (recip_id_best): TP(real cDNA) 0.615,
genuine(over-merge) 0.295, truthbar(divergent paralog) 0.320 -- the floor separates real
near-identical paralogs from divergent/domain-share over-merges by DESIGN.

THE FLOORED EDGE RULE
---------------------
  KEEP a cross-gene edge iff  core_recip>=0.19 AND aln_frac>=0.24
                              AND recip_id_best >= FLOOR
                              AND NOT(repeat-hub gate)          [everything else unchanged]
FLOOR=0.0 recovers the current divergent-inclusive catalog (md5 dca64cbd) byte-identical
(recip_id_best is always >= 0). within-gene / unannotated edges are NEVER floored.

SWEEP + EVAL (via the shipped bench/family_level_pr_current.py path)
-------------------------------------------------------------------
FLOOR in {off, 0.70, 0.75, 0.80, 0.85, 0.90}. For each: build the floored catalog, write
it to a temp TSV in the family_rna_refine.tsv format, and run family_level_pr_current.main()
on it (redirected outputs) to obtain -- from the SAME machinery the committed metrics use:
  R_oracle       : diploid-DNA-oracle multi-copy gene recovery (truth 3; with the shipped
                   gene-projection relabel credit)
  P_oracle       : 1-(allele+oversize+multifam)/oracle_mapped (taskformula) + dedup variant
  E_p purity     : protein-family block purity (truth 1)
  distinct_fp    : distinct diploid-oracle-offending blocks
  n_families     : multi-copy families
  FP composition : allele / oversize / multifam / E_p-impure-not-DNA-named
Plus the DOMAIN-SHARE exclusion probe (RHD+SDHD, MOV10+RHOC, RBBP4+SYNC, DEDD+NIT1 still
co-membered?). RNA-only (no DNA in any edge decision); the DNA/protein oracle only SCORES.

Deterministic (PYTHONHASHSEED=0, gamma=0.20 seed=0, sorted writes).
Writes: bench/divergence_floor.tsv (the P/R-vs-floor curve)
Run:    PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/divergence_floor.py
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import argparse
import contextlib
import io
import json
import statistics
import tempfile
from collections import Counter

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)

import family_er_pr as FP
import family_rna_refine as FR
import family_level_pr_current as PRC
import genome_family_def as G
import graph_def_refine_sweep as SW
import rna_only_edge_oracle as RO

RECIP_TSV = os.path.join(BENCH, "ri_sharedlen_recip_id.tsv")
OUT_TSV = os.path.join(BENCH, "divergence_floor.tsv")
BUILDER = os.path.join(BENCH, "ri_build_recip_id.py")

# swept floors (the P/R-vs-floor curve) + the off baseline (== current default dca64cbd)
FLOORS = [0.70, 0.75, 0.80, 0.85, 0.90]
NAMED_DOMAIN_SHARES = [("RHD", "SDHD"), ("MOV10", "RHOC"), ("RBBP4", "SYNC"), ("DEDD", "NIT1")]


# --------------------------------------------------------------------------- recip-id cache
def load_recip_id(metric="best"):
    """gene-pair (frozenset) -> reciprocal whole-transcript identity (recip_id_best by
    default; recip_id_sum available). Only cross-gene pairs passing core+aln are present
    (the only edges a floor can cut). Absent pair => not in dict (edge kept: absence is a
    technical missing-sequence, not a divergence signal -- and there are 0 such rows)."""
    if not os.path.exists(RECIP_TSV):
        import subprocess
        r = subprocess.run([sys.executable, BUILDER], env=dict(os.environ, PYTHONHASHSEED="0"))
        if r.returncode != 0:
            raise RuntimeError("ri_build_recip_id.py failed")
    col = "recip_id_sum" if metric == "sum" else "recip_id_best"
    d = {}
    with open(RECIP_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            v = f[I[col]]
            if v == "":
                continue
            d[frozenset((f[I["gene_a"]], f[I["gene_b"]]))] = float(v)
    return d


def recip_class_stats(metric="best"):
    """class medians/means + per-floor retained fractions on the core+aln-passing edges."""
    col = "recip_id_sum" if metric == "sum" else "recip_id_best"
    by = {"TP": [], "genuine": [], "truthbar": []}
    with open(RECIP_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[I[col]] == "":
                continue
            cls = f[I["cls"]]
            if cls in by:
                by[cls].append(float(f[I[col]]))
    stats = {}
    for cls, xs in by.items():
        stats[cls] = dict(n=len(xs),
                          median=round(statistics.median(xs), 4) if xs else None,
                          mean=round(statistics.mean(xs), 4) if xs else None,
                          frac_ge={f"{fl:.2f}": round(sum(1 for x in xs if x >= fl) / len(xs), 4)
                                   for fl in FLOORS} if xs else {})
    return stats


# --------------------------------------------------------------------------- floored build (COPY of FR.build_catalog + floor)
def build_catalog_floored(floor, recip_id, repeat_gate=True, gamma=FR.GAMMA,
                          split_recombinants=True, repeat_bridge_gate=True):
    """A COPY of family_rna_refine.build_catalog with ONE added edge conjunct: a cross-gene
    edge is kept only if its reciprocal whole-transcript identity recip_id[k] >= floor. floor
    is None or 0.0 => the conjunct is vacuous and the catalog is byte-identical to the shipped
    default. Everything else (core+aln, repeat-hub gate, gamma, recombinant-split, multi-repeat-
    bridge, allele-demote) is UNCHANGED. RNA-only; no DNA read here."""
    FR.rna_only_guard()
    floor = 0.0 if floor is None else float(floor)

    # ---- RNA features (exact oracle loaders, same as FR) ----
    gene_of_dn = RO.load_gene_of_dn()
    pair_core = RO.load_pair_core(gene_of_dn)
    univ_aln = RO.load_universal_aln()
    allele = RO.load_allele()
    pair_repeat_mult = FR.load_repeat_mult() if repeat_gate else {}

    # ---- shipped graph context ----
    meta = FP.load_meta(); annot = FP.load_annot(); gene_of = FP.gene_of_factory(annot)
    raw_fams = FP.load_raw_families(); edge_pairs = FP.load_edges()
    genes, gene_of_dn2, *_ = FP.build_genes_dict(raw_fams, meta, gene_of)
    all_nodes = set()
    for f in raw_fams:
        all_nodes.update(f)

    def core_aln_keep(k):
        c = pair_core.get(k); c = c if c is not None else 0.0
        a = univ_aln.get(k); a = a if a is not None else 0.0
        return (c >= FR.CORE_MIN) and (a >= FR.ALN_MIN)

    def repeat_hub(k):
        m = pair_repeat_mult.get(k)
        return (m is not None) and (m >= FR.REPEAT_MULT_MIN)

    def recip_ok(k):
        # DIVERGENCE FLOOR: keep only if reciprocal identity >= floor. Absent (missing-seq,
        # 0 rows) => keep (technical, not a divergence signal). floor==0 => always True.
        if floor <= 0.0:
            return True
        r = recip_id.get(k)
        return True if r is None else (r >= floor)

    import networkx as nx
    Gr = nx.Graph(); Gr.add_nodes_from(all_nodes)
    for a, b in edge_pairs:
        if a in all_nodes and b in all_nodes:
            Gr.add_edge(a, b)

    kept = set()
    n_cross_cut_floor = 0
    floor_cut_pairs = set()
    for u, v in Gr.edges():
        ga, gb = gene_of_dn2.get(u), gene_of_dn2.get(v)
        if ga is None or gb is None or ga == gb:
            kept.add(frozenset((u, v)))            # within-gene / unannotated: never floored
            continue
        k = frozenset((ga, gb))
        keep_ca = core_aln_keep(k)
        repeat_only = repeat_gate and repeat_hub(k)
        floor_ok = recip_ok(k)
        if keep_ca and not repeat_only and floor_ok:
            kept.add(frozenset((u, v)))
        elif keep_ca and not repeat_only and not floor_ok:
            n_cross_cut_floor += 1
            floor_cut_pairs.add(k)                 # cut SPECIFICALLY by the divergence floor

    # ---- shipped gamma-quasi-clique refinement ----
    comps = SW.components_from_edges(all_nodes, kept)
    refined = G.refine_families(comps, [tuple(e) for e in kept], genes, gamma, FR.SEED)

    # ---- recombinant-split gate ----
    if split_recombinants:
        import recombinant_split as RS
        refined, _ = RS.split_families(refined, genes, gene_of_dn2)
        refined = [b for b in refined if G.distinct_loci(b, genes) >= 2]

    # ---- multi-repeat-bridge gate ----
    if repeat_bridge_gate:
        import multi_repeat_bridge_gate as MRB
        refined, _ = MRB.split_families_repeat_bridge(
            refined, genes, gene_of_dn2, FR.REPEAT_BRIDGE_MULT_MIN, FR.REPEAT_BRIDGE_COUNT_MIN)
        refined = [b for b in refined if G.distinct_loci(b, genes) >= 2]

    # ---- allele DEMOTE ----
    def demote_gene(g):
        a = allele.get(g)
        return a is not None and a["balanced_frac"] >= FR.DEMOTE_BAL_MIN and a["copy_like"] <= FR.DEMOTE_COPY_MAX

    catalog = []
    for b in refined:
        gs = [gene_of_dn2[dn] for dn in b if dn in gene_of_dn2]
        if gs:
            dom, cnt = Counter(gs).most_common(1)[0]
            homog = cnt / len(gs)
            if demote_gene(dom) and homog >= 0.5 and G.distinct_loci(b, genes) >= 2:
                continue
        catalog.append(sorted(b))

    return dict(catalog=catalog, gene_of_dn=gene_of_dn2, genes=genes,
                n_cross_cut_floor=n_cross_cut_floor, floor_cut_pairs=floor_cut_pairs)


# --------------------------------------------------------------------------- catalog TSV (FR.write_outputs format)
def write_catalog_tsv(path, catalog, gene_of_dn, genes):
    """Write the floored catalog in the EXACT family_rna_refine.tsv format PRC reads."""
    fams_sorted = sorted(catalog, key=lambda b: tuple(sorted(b)))
    with open(path, "w") as out:
        out.write("family_id\tn_loci\tdominant_gene\tmember_dn\tmember_gene\tchrom\tstart\tend\n")
        for fid, b in enumerate(fams_sorted):
            gs = [gene_of_dn.get(dn, "NA") for dn in b]
            dom = Counter(g for g in gs if g != "NA").most_common(1)
            dom = dom[0][0] if dom else "NA"
            nl = G.distinct_loci(b, genes)
            for dn in sorted(b):
                g = genes[dn]
                out.write(f"{fid}\t{nl}\t{dom}\t{dn}\t{gene_of_dn.get(dn, 'NA')}\t"
                          f"{g['chrom']}\t{g['start']}\t{g['end']}\n")


# --------------------------------------------------------------------------- eval via shipped PRC path
def evaluate_via_prc(catalog_tsv):
    """Run family_level_pr_current.main() on a floored catalog TSV using the SHIPPED eval
    machinery, with all outputs redirected to a temp dir (never touches the committed
    family_level_pr_current.{tsv,json}). Returns the PRC `out` dict."""
    td = tempfile.mkdtemp(prefix="divfloor_prc_")
    old = (PRC.CATALOG_TSV, PRC.OUT_TSV, PRC.OUT_JSON)
    PRC.CATALOG_TSV = catalog_tsv
    PRC.OUT_TSV = os.path.join(td, "pr.tsv")
    PRC.OUT_JSON = os.path.join(td, "pr.json")
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            out = PRC.main()
    finally:
        PRC.CATALOG_TSV, PRC.OUT_TSV, PRC.OUT_JSON = old
    return out


def domain_share_comembership(catalog_tsv, pairs):
    """Which named pairs are STILL co-membered in the floored catalog (should be [] if excluded)."""
    import csv
    from collections import defaultdict
    fam = defaultdict(set)
    with open(catalog_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fam[row["family_id"]].add(row["member_gene"])
    res = {}
    for a, b in pairs:
        res[f"{a}+{b}"] = sorted(fid for fid, gs in fam.items() if a in gs and b in gs)
    return res


def metrics_from_prc(out):
    """Pull the requested scalars out of the PRC output dict."""
    orc = out["truth3_diploid_oracle"]["current"]
    ep = out["truth1_protein_Ep"]["current"]
    roster = out["residual_fp_roster"]
    return dict(
        n_families=out["n_families_current"],
        R_oracle=orc["recall_oracle"],
        R_oracle_frac=f"{orc['oracle_genes_recovered_multicopy']}/{orc['oracle_multicopy_genes']}",
        P_oracle_task=orc["precision_oracle_taskformula"],
        P_oracle_dedup=orc["precision_oracle_dedup"],
        oracle_mapped=orc["oracle_mapped_families"],
        distinct_fp=orc["distinct_fp_blocks"],
        fp_allele=orc["n_allele"], fp_oversize=orc["n_oversize"], fp_multifam=orc["n_multifam"],
        E_p_purity=ep["precision_Ep"], E_p_impure=ep["impure_blocks"], E_p_total=ep["total_blocks"],
        ep_impure_not_dna=len(roster["ep_impure_not_dna_named"]),
        fp_multifam_names=["+".join(e["genes"]) for e in roster["multifam"]],
        fp_oversize_names=["+".join(e["genes"]) for e in roster["oversize"]],
        fp_allele_names=["+".join(e["genes"]) for e in roster["allele"]],
    )


# --------------------------------------------------------------------------- sweep
def run_sweep(metric="best"):
    recip_id = load_recip_id(metric)
    cls_stats = recip_class_stats(metric)
    rows = []
    detail = {}
    # off baseline (floor disabled) + the swept floors
    for floor in [None] + FLOORS:
        built = build_catalog_floored(floor, recip_id)
        td = tempfile.mkdtemp(prefix="divfloor_cat_")
        cat_tsv = os.path.join(td, "family_rna_refine.tsv")
        write_catalog_tsv(cat_tsv, built["catalog"], built["gene_of_dn"], built["genes"])
        out = evaluate_via_prc(cat_tsv)
        m = metrics_from_prc(out)
        ds = domain_share_comembership(cat_tsv, NAMED_DOMAIN_SHARES)
        ds_excluded = [p for p, fids in ds.items() if not fids]
        floor_label = "off" if floor is None else f"{floor:.2f}"
        md5 = _md5(cat_tsv)
        row = dict(floor=floor_label, n_families=m["n_families"],
                   R_oracle=m["R_oracle"], R_oracle_frac=m["R_oracle_frac"],
                   P_oracle_task=m["P_oracle_task"], P_oracle_dedup=m["P_oracle_dedup"],
                   E_p_purity=m["E_p_purity"], distinct_fp=m["distinct_fp"],
                   oracle_mapped=m["oracle_mapped"],
                   fp_allele=m["fp_allele"], fp_oversize=m["fp_oversize"], fp_multifam=m["fp_multifam"],
                   ep_impure_not_dna=m["ep_impure_not_dna"],
                   edges_cut_by_floor=built["n_cross_cut_floor"],
                   domain_shares_excluded=f"{len(ds_excluded)}/{len(NAMED_DOMAIN_SHARES)}",
                   catalog_md5=md5[:8])
        rows.append(row)
        detail[floor_label] = dict(metrics=m, domain_shares=ds,
                                   domain_shares_excluded=ds_excluded,
                                   n_edges_cut_by_floor=built["n_cross_cut_floor"],
                                   catalog_md5=md5)
        print(f"[floor {floor_label:>4s}] nFam={m['n_families']:3d}  "
              f"R_oracle={m['R_oracle']} ({m['R_oracle_frac']})  "
              f"P_oracle={m['P_oracle_task']}/{m['P_oracle_dedup']}  "
              f"E_p={m['E_p_purity']}  distinctFP={m['distinct_fp']}  "
              f"FP(a/o/m)={m['fp_allele']}/{m['fp_oversize']}/{m['fp_multifam']}  "
              f"domShareExcl={row['domain_shares_excluded']}  "
              f"edgesCutByFloor={built['n_cross_cut_floor']}  md5={md5[:8]}", flush=True)
    return rows, detail, cls_stats


def _md5(path):
    import hashlib
    with open(path, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


def write_curve_tsv(rows, cls_stats):
    cols = ["floor", "n_families", "R_oracle", "R_oracle_frac", "P_oracle_task", "P_oracle_dedup",
            "E_p_purity", "distinct_fp", "oracle_mapped", "fp_allele", "fp_oversize", "fp_multifam",
            "ep_impure_not_dna", "edges_cut_by_floor", "domain_shares_excluded", "catalog_md5"]
    with open(OUT_TSV, "w") as out:
        out.write("# DIVERGENCE FLOOR SWEEP :: reciprocal whole-transcript identity (recip_id_best = "
                  "matches/max(len)); RNA-only edge gate ADDED to core>=0.19 AND aln>=0.24 (+VG gates + "
                  "gamma0.20 + demote). floor=off == current default (md5 dca64cbd). Eval via shipped "
                  "family_level_pr_current.py. Deterministic PYTHONHASHSEED=0.\n")
        out.write("# class medians (recip_id_best on 1457 core+aln-passing edges): "
                  + "  ".join(f"{c}(n={cls_stats[c]['n']},med={cls_stats[c]['median']})"
                              for c in ("TP", "genuine", "truthbar")) + "\n")
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"[write] {OUT_TSV}", flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description="divergence-floor sweep (measure stage).")
    ap.add_argument("--metric", choices=["best", "sum"], default="best",
                    help="reciprocal-identity column (default best = matches/max(len))")
    ap.add_argument("--json", action="store_true", help="also emit the full detail JSON to stdout")
    args = ap.parse_args(argv)
    rows, detail, cls_stats = run_sweep(args.metric)
    write_curve_tsv(rows, cls_stats)
    if args.json:
        print("JSON " + json.dumps(dict(rows=rows, detail=detail, class_stats=cls_stats),
                                   sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
