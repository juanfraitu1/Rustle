#!/usr/bin/env python3
"""SECOND, ORTHOGONAL detection lever for currently-unplaced (singleton) genes: a sub-SD98-identity
shared-exon signal, corroborated by tight agreement in OUR OWN independently-computed famCN (via
famcn_from_wssd.py's WSSD read-depth computation -- NOT Soto's own published S1C truth table, so this is
NOT the same circularity as using S1C famCN to split an already-detected component).

WHY THIS EXISTS (docs/o1_ledger.md, follow-up to §6im). ID_328 (CU633904 paralogs, 8 members, all
unprocessed_pseudogene/eligible) was diagnosed (§6ik) as a `sedef_coverage_gap_below_98pct` case -- gene-
body identity ~0.9139-0.9151, under the SD98 gate. Investigating further (this script's origin) found TWO
compounding, more specific causes, not a fundamental ceiling:
  1. SEDEF calls one wide (~230kb) duplication UNIT spanning each gene pair, and reports ONE identity
     figure averaged over that whole unit -- but a direct realignment of just the ~14kb gene body itself
     (minimap2, both loci extracted with padding) shows up to 99.99% local identity. The unit-average
     dilutes a highly-conserved core with more divergent flanking sequence; gating on the row's declared
     identity is gating on the wrong statistic for this case.
  2. Independently, every one of these wide units' outer bounds ALSO straddles each acrocentric
     chromosome's own liftover-guard span (the §6il/§6im acrofix anchors don't reach this far into the
     rearranged region) -- confirmed via direct realignment of each side's FULL row span (or wider) against
     v1.0: chr14/15/21/22/13 all realign at 100% identity (or 99.999%+) end to end, at exactly their
     already-known "before"-regime offset, all the way past the row's own bounds with real margin. New,
     deeper anchors for these 5 chromosomes (in acro_extra_anchors.tsv) fix the liftover side.

Neither fix alone recovers ID_328 through the STANDARD >=98%-identity pipeline (the row's declared
identity is still ~0.914, so it never reaches the liftover check there at all). This script is a SEPARATE,
disclosed, opt-in lever, run ONLY over genes with no edge in the standard graph: lower the identity floor
(--min-identity, e.g. 0.85) to admit the same wide, real SEDEF calls SEDEF found for these loci, restricted
to the singleton population (so it can never perturb an already-correctly-formed family), then require CN
CORROBORATION (tight agreement in independently-computed famCN) before accepting a merge -- exactly the
"E_c earns an orthogonal role, used only where the primary graph is blind" precedent from the actual O1
side of this project (project_o1_shared_read_edges.md §6bi), not a threshold relaxation applied wholesale.

Two safeguards against the failure mode this mirrors (E_c-as-additive-merge-tier, refuted on the O1 side
for exactly this reason -- near-universal ambient repeat homology, invariant to any depth threshold):
  (a) scoped to the singleton population ONLY, never touching an already-placed gene or bridging two
      already-formed families (unlike O1's E_c, which was tested genome-wide against a moving target);
  (b) requires BOTH a real (if sub-98%) structural signal AND independent CN corroboration, not either
      alone -- CN values cluster near common integers by chance, and weak structural signal alone is
      exactly what the O1-side investigation showed is near-universal.

Usage: soto_weak_edge_cn_merge.py --weak-edges weak_edges_singletons2.tsv --geneset singleton_357_geneset.tsv
                                   --famcn famcn_ours_all.tsv --mad 1.0 --mad-statistic median
                                   --out cn_merge_families.tsv
"""
import argparse, csv, os, sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from soto_cluster_from_shared import ELIGIBLE, mad_mean, mad_median  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--weak-edges", required=True)
    ap.add_argument("--geneset", required=True, help="singleton-only geneset (gene_id, biotype)")
    ap.add_argument("--famcn", required=True, help="OUR OWN famCN, e.g. famcn_ours_all.tsv "
                     "(gene_id, famCN, famCN_mad, n_samples) -- NOT Soto's S1C table")
    ap.add_argument("--mad", type=float, default=1.0)
    ap.add_argument("--mad-statistic", choices=["mean", "median"], default="median")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    mad_fn = mad_median if a.mad_statistic == "median" else mad_mean

    biotype = {}
    with open(a.geneset) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            biotype[r["gene_id"]] = r.get("biotype", "")

    famcn = {}
    with open(a.famcn) as fh:
        for r in csv.reader(fh, delimiter="\t"):
            try:
                famcn[r[0]] = float(r[1])
            except (IndexError, ValueError):
                continue

    adj = defaultdict(set)
    with open(a.weak_edges) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ga, gb = r["gene_a"], r["gene_b"]
            if biotype.get(ga) in ELIGIBLE and biotype.get(gb) in ELIGIBLE:
                adj[ga].add(gb)
                adj[gb].add(ga)

    seen, comps = set(), []
    for g in adj:
        if g in seen:
            continue
        stack, comp = [g], []
        seen.add(g)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        comps.append(comp)

    accepted, rejected = [], []
    for comp in comps:
        vals = [famcn[g] for g in comp if g in famcn]
        if len(comp) >= 2 and len(vals) >= 2 and mad_fn(vals) < a.mad:
            accepted.append((comp, vals))
        else:
            rejected.append((comp, vals))

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "biotype", "family_id", "n_members", "our_famCN", "cn_mad", "status"])
        for i, (comp, vals) in enumerate(sorted(accepted, key=lambda cv: -len(cv[0]))):
            fid = f"CNWEAKFAM{i}"
            m = mad_fn(vals)
            for g in sorted(comp):
                w.writerow([g, biotype.get(g, ""), fid, len(comp),
                            f"{famcn[g]:.2f}" if g in famcn else "", f"{m:.3f}", "cn_weak_edge_family"])

    print(f"[weak-edge components] {len(comps)} raw eligible-only components", file=sys.stderr)
    print(f"[accepted] {len(accepted)} components pass CN agreement (mad<{a.mad}, {a.mad_statistic}) "
          f"-> {sum(len(c) for c, _ in accepted)} genes -> {a.out}", file=sys.stderr)
    print(f"[rejected] {len(rejected)} components fail CN agreement or lack CN data", file=sys.stderr)
    for comp, vals in sorted(rejected, key=lambda cv: -len(cv[0]))[:10]:
        print(f"  rejected: {len(comp)} genes, CN values={[round(v,2) for v in vals]}", file=sys.stderr)


if __name__ == "__main__":
    main()
