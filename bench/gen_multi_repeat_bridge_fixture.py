#!/usr/bin/env python3
"""BYTE-PARITY golden-fixture generator for the Rust port of the O1 multi-repeat-bridge
gate `bench/multi_repeat_bridge_gate.py` (MRB) -> `vg_family/multi_repeat_bridge.rs`
(the O1 over-merge-gate migration, RETIREMENT_AND_MIGRATION.md Tier-2 #7, full byte-parity).

It imports the SHIPPED `multi_repeat_bridge_gate` (MRB) + the shared exon-graph library
`recombination_bridge_detector` (RBD) + `vg_repeat_catalog` (V) and reconstructs the EXACT
(refined, genes, gene_of_dn) that `family_rna_refine.build_catalog` feeds to
`MRB.split_families_repeat_bridge` -- the SAME proven prefix as
`gen_recombinant_split_fixture.py`, THEN the recombinant-split stage + the >=2-distinct-loci
filter (build_catalog lines 388-392), i.e. the post-recombinant-split catalog the gate runs on.

For every selected block it runs the REAL per-block `MRB.split_families_repeat_bridge([b], ...)`
(so the (parts, split_info) are the shipped function's own output) and captures a SELF-CONTAINED
record: the block + the genes/gene_of_dn/strand + the RBD exon-list skeleton subset (family_exons)
+ the V introns-string skeleton + F.load_meta subset (locus_node_set) + the node mult/cyclic for
every involved VG node + the RAW case-carrying genome slices via a mock fetcher (so no 3.4G genome
is needed at test time) + the CHARACTERIZE plan.

Coverage (the task's required cases):
  * ALL 61 CUT (DISCONNECTED) families incl PDIA3/PRKAB2 (one block carries both);
  * the GSTM2 / MAGE SPARED negatives (connected single-component domain-share / cardinality
    families that STRUCTURALLY cannot fire -- the load-bearing no-full-shared-exon conjunct);
  * a >MAX_MEMBERS SKIP family (ZNF92, n=41 recs -> early return, never cut);
  * same-gene-guard-BLOCKED disconnected families (guard_ok=0 -> spared);
  * representative NON-cut families (connected + multi-component-but-not-bridged).

Aggregate: split_families_repeat_bridge over an ordered block list = the per-block results
concatenated (the function is a per-block loop), so the DEFAULT (cost-budgeted) and FULL
aggregates are derived from the per-block (parts, split_info) captures.

Two-phase (the reconstruct + per-family characterize is the slow part): phase 1 pickles the
per-block captures; phase 2 (cheap) assigns the default cost budget + writes the JSON, so the
budget can be tuned without re-running the heavy reconstruction.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/gen_multi_repeat_bridge_fixture.py
Deterministic (sorted selection, PYTHONHASHSEED=0).
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/..")  # repo root (relative paths)

import collections
import json
import pickle
import struct
import tempfile

sys.path.insert(0, "bench")
import family_er_pr as FP
import genome_family_def as G
import graph_def_refine_sweep as SW
import rna_only_edge_oracle as RO
import family_rna_refine as FR
import recombination_bridge_detector as RBD
import recombinant_split as RS
import multi_repeat_bridge_gate as MRB
import vg_repeat_catalog as V
import pysam
import networkx as nx

OUT = os.path.join("src", "rustle", "vg_family", "testdata", "multi_repeat_bridge_fixture.json")
# phase-1 (heavy reconstruct + capture) cache -- a pure speed optimisation; delete to rebuild from
# scratch (deterministic under PYTHONHASHSEED=0). Override with MRB_FIXTURE_CACHE.
CACHE = os.environ.get("MRB_FIXTURE_CACHE", os.path.join(tempfile.gettempdir(), "mrb_phase1.pkl"))

T, C = MRB.WIRED_MULT_MIN, MRB.WIRED_COUNT_MIN  # 8, 2


def fbits(x):
    return struct.unpack("<Q", struct.pack("<d", float(x)))[0]


# ------------------------------------------------------------------ reconstruction (== build_catalog prefix + recombinant split)
def reconstruct():
    FR.rna_only_guard()
    gene_of_dn = RO.load_gene_of_dn()
    pair_core = RO.load_pair_core(gene_of_dn)
    univ_aln = RO.load_universal_aln()
    pair_repeat_mult = FR.load_repeat_mult()
    meta = FP.load_meta()
    annot = FP.load_annot()
    gene_of = FP.gene_of_factory(annot)
    raw_fams = FP.load_raw_families()
    edge_pairs = FP.load_edges()
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

    Gr = nx.Graph(); Gr.add_nodes_from(all_nodes)
    for a, b in edge_pairs:
        if a in all_nodes and b in all_nodes:
            Gr.add_edge(a, b)
    kept = set()
    for u, v in Gr.edges():
        ga, gb = gene_of_dn2.get(u), gene_of_dn2.get(v)
        if ga is None or gb is None or ga == gb:
            kept.add(frozenset((u, v))); continue
        k = frozenset((ga, gb))
        if core_aln_keep(k) and not repeat_hub(k):
            kept.add(frozenset((u, v)))
    comps = SW.components_from_edges(all_nodes, kept)
    refined0 = G.refine_families(comps, [tuple(e) for e in kept], genes, FR.GAMMA, FR.SEED)
    # recombinant-split stage + >=2-distinct-loci filter (build_catalog 388-392)
    refined1, _ = RS.split_families(refined0, genes, gene_of_dn2)
    refined2 = [b for b in refined1 if G.distinct_loci(b, genes) >= 2]
    return refined2, genes, gene_of_dn2


# ------------------------------------------------------------------ per-block capture
def dp_cost(recs):
    n = len(recs); c = 0
    for a in range(n):
        for b in range(n):
            if a == b:
                continue
            for qa in recs[a]["exons"]:
                if len(qa) < RBD.MIN_EXON:
                    continue
                for tb in recs[b]["exons"]:
                    if len(tb) < RBD.MIN_EXON:
                        continue
                    lo, hi = (len(qa), len(tb)) if len(qa) <= len(tb) else (len(tb), len(qa))
                    if lo / hi < 0.5:
                        continue
                    c += 2 * len(qa) * len(tb)
    return c


def plan_json(pl):
    if pl is None:
        return None
    base = dict(
        skipped=bool(pl.get("skipped")),
        n_members=pl["n_members"],
        n_distinct_genes=pl["n_distinct_genes"],
        same_gene_family=bool(pl["same_gene_family"]),
    )
    if pl.get("skipped"):
        base.update(n_comp=0, connected=False, shares_full_exon=False, xcomp_bits=None,
                    no_full_shared_exon=False, n_cross_pairs=0, n_cross_pairs_with_node=0,
                    n_cross_shared_nodes=0, max_cross_mult=0, med_cross_mult=0, n_cross_cyclic=0,
                    rep_bridges=[], n_internal_shared_nodes=0, max_internal_mult=0,
                    n_internal_repeat_m8=0, min_comp_size=0, separates_same_gene=False,
                    guard_ok=False, comps_dn=[], all_dn=[])
        return base
    xid = pl["xcomp_best_exon_id"]
    base.update(
        n_comp=pl["n_comp"], connected=bool(pl["connected"]),
        shares_full_exon=bool(pl["shares_full_exon"]),
        xcomp_bits=(fbits(xid) if xid is not None else None),
        no_full_shared_exon=bool(pl["no_full_shared_exon"]),
        n_cross_pairs=pl["n_cross_pairs"], n_cross_pairs_with_node=pl["n_cross_pairs_with_node"],
        n_cross_shared_nodes=pl["n_cross_shared_nodes"], max_cross_mult=pl["max_cross_mult"],
        med_cross_mult=pl["med_cross_mult"], n_cross_cyclic=pl["n_cross_cyclic"],
        rep_bridges=[[t, pl["rep_bridges"][t]] for t in sorted(pl["rep_bridges"])],
        n_internal_shared_nodes=pl["n_internal_shared_nodes"],
        max_internal_mult=pl["max_internal_mult"], n_internal_repeat_m8=pl["n_internal_repeat_m8"],
        min_comp_size=pl["min_comp_size"], separates_same_gene=bool(pl["separates_same_gene"]),
        guard_ok=bool(pl["guard_ok"]),
        comps_dn=[sorted(c) for c in pl["comps_dn"]], all_dn=list(pl["all_dn"]),
    )
    return base


def capture(block, genes, gene_of_dn, res, fa):
    """Self-contained record: inputs + REAL characterize plan + REAL per-block split."""
    blk = sorted(block)
    members = MRB._block_members(blk, genes, gene_of_dn)
    R_skel, strand, vskel, meta = res["R_skel"], res["strand"], res["vskel"], res["meta"]
    mult, cyclic = res["mult"], res["cyclic"]

    pl = MRB.characterize(-1, "w", "w", members, R_skel, strand, fa, meta, vskel, fa, mult, cyclic)
    # REAL per-block split (the shipped function's own (parts, split_info))
    nr, si = MRB.split_families_repeat_bridge([blk], genes, gene_of_dn, T, C)
    parts = [list(p) for p in nr]
    split_info = si[0] if si else None

    recs = RBD.family_exons(members, R_skel, strand, fa)
    cost = dp_cost(recs) if (pl is not None and not pl.get("skipped")) else 0

    genes_out = {m["dn"]: [m["chrom"], m["start"], m["end"]] for m in members}
    gof_out = {dn: gene_of_dn[dn] for dn in blk if dn in gene_of_dn}
    strand_out = {m["dn"]: strand.get(m["dn"], "+") for m in members}

    # family_exons skeleton subset + slices (fa_up path)
    slices = {}
    r_skel_out = []
    seen_r = set()
    for m in members:
        key = (m["chrom"], m["start"], m["end"])
        exons = R_skel.get(key)
        if exons is not None and key not in seen_r:
            seen_r.add(key)
            r_skel_out.append(dict(chrom=key[0], start=key[1], end=key[2],
                                   exons=[list(e) for e in exons]))
            for (s, e) in exons:
                slk = (m["chrom"], s - 1, e)
                if slk not in slices:
                    slices[slk] = fa.fetch(m["chrom"], s - 1, e)

    # locus_node_set subset (only when NOT skipped: characterize reaches nodesets over RECS)
    meta_out = {}
    vskel_out = []
    seen_v = set()
    node_union = set()
    if pl is not None and not pl.get("skipped"):
        for r in recs:
            dn = r["dn"]
            mm = meta.get(dn)
            if mm is None:
                continue
            meta_out[dn] = [mm[0], mm[1], mm[2]]
            vkey = (mm[0], mm[1], mm[2])
            if vkey not in seen_v:
                seen_v.add(vkey)
                if vkey in vskel:
                    vskel_out.append(dict(chrom=vkey[0], start=vkey[1], end=vkey[2],
                                          introns=vskel[vkey]))
            chrom, exons = V.dn_exons(dn, meta, vskel)
            if exons:
                for (a, b) in exons:
                    slk = (chrom, a, b)
                    if slk not in slices:
                        slices[slk] = fa.fetch(chrom, a, b)
            ns = MRB.locus_node_set(dn, meta, vskel, fa)
            if ns:
                node_union |= ns
    node_mult_out = {str(c): [mult[c], cyclic[c]] for c in sorted(node_union) if c in mult}

    return dict(
        block=blk,
        genes=genes_out, gene_of_dn=gof_out, strand=strand_out,
        r_skel=r_skel_out, meta=meta_out, vskel=vskel_out,
        node_mult=node_mult_out,
        slices=[dict(chrom=k[0], start0=k[1], end0=k[2], bytes=v) for k, v in slices.items()],
        plan=plan_json(pl),
        parts=parts, split_info=split_info,
        cut=(len(parts) > 1), skipped=bool(pl.get("skipped")) if pl else False,
        n_recs=len(recs), cost=cost,
    )


# ------------------------------------------------------------------ phase 1 (heavy) -> pickle
def phase1():
    if os.path.exists(CACHE):
        print("[phase1] loading cache", CACHE, flush=True)
        with open(CACHE, "rb") as fh:
            return pickle.load(fh)

    print("[phase1] reconstruct refined2 ...", flush=True)
    refined2, genes, gene_of_dn = reconstruct()
    print("  refined2 blocks:", len(refined2))
    res = MRB._stage_loaders()
    fa = res["fa_up"]

    # scan every block: classify + dominant + member genes
    def domgenes(b):
        cnt = collections.Counter(gene_of_dn.get(dn) for dn in b if gene_of_dn.get(dn) not in (None, "NA"))
        dom = cnt.most_common(1)[0][0] if cnt else "NA"
        mg = sorted(set(gene_of_dn.get(dn) for dn in b if gene_of_dn.get(dn)))
        return dom, mg

    scan = []
    for b in refined2:
        members = MRB._block_members(b, genes, gene_of_dn)
        if len(members) < 2:
            continue
        pl = MRB.characterize(-1, "w", "w", members, res["R_skel"], res["strand"], fa,
                              res["meta"], res["vskel"], fa, res["mult"], res["cyclic"])
        if pl is None:
            continue
        parts, _ = MRB.split_family_repeat_bridge(b, genes, gene_of_dn, T, C, res=res)
        dom, mg = domgenes(b)
        scan.append(dict(block=sorted(b), dom=dom, mg=mg, pl_lite=dict(
            skipped=bool(pl.get("skipped")),
            n_comp=pl.get("n_comp"), connected=pl.get("connected"),
            no_full=pl.get("no_full_shared_exon"),
            rep8=(pl["rep_bridges"].get(8) if not pl.get("skipped") else None),
            guard=pl.get("guard_ok")), cut=(len(parts) > 1)))

    # ---- select candidate blocks ----
    def has(mg, pred):
        return any(g and pred(g) for g in mg)

    cut = [s for s in scan if s["cut"]]
    spared_gstm2 = [s for s in scan if not s["cut"] and not s["pl_lite"]["skipped"]
                    and s["pl_lite"]["connected"] and has(s["mg"], lambda g: g.startswith("GSTM2"))]
    spared_mage = [s for s in scan if not s["cut"] and not s["pl_lite"]["skipped"]
                   and s["pl_lite"]["connected"] and has(s["mg"], lambda g: g.startswith("MAGE"))]
    guard = [s for s in scan if not s["cut"] and not s["pl_lite"]["skipped"]
             and s["pl_lite"]["no_full"] and (s["pl_lite"]["rep8"] or 0) >= C and not s["pl_lite"]["guard"]]
    skip = [s for s in scan if s["pl_lite"]["skipped"]]
    noncut_conn = [s for s in scan if not s["cut"] and not s["pl_lite"]["skipped"]
                   and s["pl_lite"]["n_comp"] == 1]
    noncut_multi = [s for s in scan if not s["cut"] and not s["pl_lite"]["skipped"]
                    and (s["pl_lite"]["n_comp"] or 0) >= 2]

    print(f"  scan: cut={len(cut)} spared_gstm2={len(spared_gstm2)} spared_mage={len(spared_mage)} "
          f"guard={len(guard)} skip={len(skip)} noncut_conn={len(noncut_conn)} noncut_multi={len(noncut_multi)}")

    selected = {}  # key = tuple(block) -> (group, block)

    def sel(group, s):
        selected.setdefault(tuple(s["block"]), (group, s["block"]))

    for s in cut:
        sel("cut", s)
    # spared: keep the cheapest few of each (by n loci as a proxy; costs added later)
    for s in sorted(spared_gstm2, key=lambda s: len(s["block"])):
        sel("spared_gstm2", s)
    for s in sorted(spared_mage, key=lambda s: len(s["block"])):
        sel("spared_mage", s)
    for s in sorted(guard, key=lambda s: len(s["block"]))[:4]:
        sel("guard", s)
    for s in skip:
        sel("skip", s)
    for s in sorted(noncut_conn, key=lambda s: len(s["block"]))[:4]:
        sel("noncut", s)
    for s in sorted(noncut_multi, key=lambda s: len(s["block"]))[:4]:
        sel("noncut", s)

    # capture each selected block
    print(f"  capturing {len(selected)} blocks ...", flush=True)
    records = []
    for i, (key, (group, block)) in enumerate(sorted(selected.items())):
        rec = capture(block, genes, gene_of_dn, res, fa)
        rec["group"] = group
        records.append(rec)
        if (i + 1) % 10 == 0:
            print(f"    {i+1}/{len(selected)}", flush=True)

    out = dict(records=records, n_input=len(refined2), n_cut_full=len(cut))
    with open(CACHE, "wb") as fh:
        pickle.dump(out, fh)
    print("[phase1] cached ->", CACHE)
    return out


# ------------------------------------------------------------------ phase 2 (cheap) -> JSON
DEFAULT_BUDGET_TOTAL = 120_000_000  # ~ debug HW-DP budget for the default parity test (~11s)


def phase2(ph1):
    records = ph1["records"]
    # deterministic order = sorted by block first element (stable)
    records.sort(key=lambda r: r["block"][0])

    # assign `default`: a COST-BUDGETED subset that MUST contain the required cases (PDIA3/PRKAB2
    # cut, one GSTM2-spared, one MAGE-spared, the >MAX_MEMBERS skip, guard-blocked, non-cut),
    # then the cheapest remaining cuts under budget.  The FULL (all-record) set is the `#[ignore]`
    # release-run test.  Cheapest-by-cost picks keep the default `cargo test` snappy.
    def carries(r, gene):
        return any(g == gene for g in r["gene_of_dn"].values())

    for r in records:
        r["default"] = False
    used = 0

    def mark(r):
        nonlocal used
        if not r["default"]:
            r["default"] = True
            used += r["cost"]

    def cheapest_n(group, n):
        return sorted([r for r in records if r["group"] == group], key=lambda r: r["cost"])[:n]

    # required non-cut cases (cheapest representative of each)
    for r in cheapest_n("spared_gstm2", 1):   # ~55M (floor: GSTM2 domain-share is a big family)
        mark(r)
    for r in cheapest_n("spared_mage", 1):    # ~10M
        mark(r)
    for r in cheapest_n("guard", 2):          # HNRNPA3(0) + CXCL16(0.24M)
        mark(r)
    for r in cheapest_n("skip", 1):           # ZNF92 skip (no tensor -> ~0)
        mark(r)
    for r in cheapest_n("noncut", 2):         # a connected + a multi-comp non-cut (cheap)
        mark(r)
    # PDIA3/PRKAB2 cut block (required)
    for r in records:
        if r["group"] == "cut" and (carries(r, "PDIA3") or carries(r, "PRKAB2")):
            mark(r)
    # fill with cheapest cut blocks under budget
    for r in sorted([r for r in records if r["group"] == "cut" and not r["default"]],
                    key=lambda r: r["cost"]):
        if used + r["cost"] > DEFAULT_BUDGET_TOTAL:
            continue
        mark(r)

    # cheap aggregate subset (`agg`): the split_families wrapper + split_info bookkeeping is
    # order-independent per block, so a CHEAP subset (cost < AGG budget) suffices to test it
    # without re-running the expensive GSTM2/PDIA3 HW-DP a second time.
    AGG_BUDGET = 20_000_000
    for r in records:
        r["agg"] = bool(r["default"] and r["cost"] < AGG_BUDGET)

    # ---- aggregates (split_families over an ordered block list = per-block concat) ----
    def aggregate(recs):
        nr, si = [], []
        for r in recs:
            nr.extend([list(p) for p in r["parts"]])
            if r["split_info"] is not None:
                si.append(r["split_info"])
        return nr, si

    default_recs = [r for r in records if r["default"]]
    agg_recs = [r for r in records if r["agg"]]
    default_new_refined, default_split_info = aggregate(default_recs)
    agg_new_refined, agg_split_info = aggregate(agg_recs)
    full_new_refined, full_split_info = aggregate(records)

    n_cut_default = sum(1 for r in default_recs if r["cut"])
    n_cut_agg = sum(1 for r in agg_recs if r["cut"])
    n_cut_full = sum(1 for r in records if r["cut"])

    # strip internal helpers from records for the JSON (keep everything the Rust test needs)
    for r in records:
        r.pop("n_recs", None)

    fixture = dict(
        WIRED_MULT_MIN=T, WIRED_COUNT_MIN=C, MAX_MEMBERS=MRB.MAX_MEMBERS,
        MULT_GRID=list(MRB.MULT_GRID), ID_THRESH=MRB.ID_THRESH, K=V.K,
        n_input=ph1["n_input"], n_cut_full=ph1["n_cut_full"],
        n_cut_default=n_cut_default, n_cut_agg=n_cut_agg, n_cut_selected=n_cut_full,
        default_budget=used,
        records=records,
        default_new_refined=default_new_refined, default_split_info=default_split_info,
        agg_new_refined=agg_new_refined, agg_split_info=agg_split_info,
        full_new_refined=full_new_refined, full_split_info=full_split_info,
    )
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(fixture, fh, indent=1, sort_keys=True)
    sz = os.path.getsize(OUT)
    print(f"[phase2] wrote {os.path.normpath(OUT)} ({sz:,} bytes)")
    print(f"  records={len(records)} default={len(default_recs)} default_cost={used:,} "
          f"(~{used/11e6:.1f}s debug)")
    print(f"  cut: selected={n_cut_full} (full n_cut={ph1['n_cut_full']}) default_cut={n_cut_default}")
    groups = collections.Counter(r["group"] for r in records)
    print("  groups:", dict(groups))
    # sanity
    assert n_cut_full == ph1["n_cut_full"], "selected cut != full cut (should include ALL cut)"


def main():
    ph1 = phase1()
    phase2(ph1)


if __name__ == "__main__":
    main()
