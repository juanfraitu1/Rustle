#!/usr/bin/env python
"""HEAD-TO-HEAD: DNA family catalog (E_a segdup, genome_families_refined.tsv) vs
RNA family catalog (E_r transcript-homology, family_rna_refine.tsv), scored against the
SAME non-circular truths:
  (a) PROTEIN E_p  : protein_families_refined.tsv co-membership (non-mega) + aln.m8 direct
                     reciprocal edges. Family-level purity (Ep-impure block rate) + edge-level
                     Ep-homologous pair precision.
  (b) cDNA 90%     : family_def_dna_pr_edges.tsv in_dna_loose (id>=90% cov>=30%). Pair precision
                     over ALL within-family pairs AND over EVALUABLE pairs (present as a truth row).
  (c) annotation   : gene-symbol-root dominant purity + single-root-family fraction.

Both catalogs are reduced to the SAME object: family -> set of distinct RefSeq gene symbols,
then clique-expanded to unordered within-family gene pairs. Identical scoring on both.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/family_dna_vs_rna_headtohead.py
"""
import json
import os
import re
import sys
from collections import defaultdict, Counter
from itertools import combinations

BENCH = os.path.dirname(os.path.abspath(__file__))
SCRATCH = "/home/juanfra/winloci_scratch"
ALN = os.path.join(SCRATCH, "aln.m8")
ANNOT = os.path.join(SCRATCH, "annot_intervals.tsv")
DNA_CAT = os.path.join(BENCH, "genome_families_refined.tsv")
RNA_CAT = os.path.join(BENCH, "family_rna_refine.tsv")
TRUTH = os.path.join(BENCH, "family_def_dna_pr_edges.tsv")
PRFAM = os.path.join(BENCH, "protein_families_refined.tsv")

ALPHA_P = 1e-5
C_P = 0.50
MEGA_MIN = 500


# ------------------------------------------------------------------ loaders
def load_dna_catalog():
    """DNA family -> set of distinct gene symbols (col 'genes')."""
    fams = {}
    with open(DNA_CAT) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        gi = hdr.index("genes")
        fi = hdr.index("family_id")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            genes = {g for g in f[gi].split(",") if g}
            fams[f[fi]] = genes
    return fams


def load_rna_catalog():
    """RNA family -> set of distinct member_gene symbols (already projected DN->gene)."""
    fams = defaultdict(set)
    with open(RNA_CAT) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        fi = hdr.index("family_id")
        mgi = hdr.index("member_gene")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            g = f[mgi]
            # keep only real gene symbols (drop unmapped DN ids / blanks)
            if g and not g.startswith("DN_"):
                fams[f[fi]].add(g)
    return dict(fams)


def load_prfam():
    """gene -> Ep family id; set of MEGA Ep ids (>=MEGA_MIN genes)."""
    g2f = {}
    sizes = {}
    with open(PRFAM) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            fid = f[0]
            genes = f[11].split(",")
            sizes[fid] = len(genes)
            for g in genes:
                g2f[g] = fid
    mega = {fid for fid, n in sizes.items() if n >= MEGA_MIN}
    return g2f, mega, sizes


def load_ep_direct():
    """reciprocal whole-protein edges from aln.m8 (evalue<=ALPHA_P, qcov,tcov>=C_P), keyed by gene."""
    fwd = {}
    with open(ALN) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            q, t = f[0], f[1]
            if q == t:
                continue
            try:
                ev, pid, qc, tc = float(f[2]), float(f[3]), float(f[4]), float(f[5])
            except ValueError:
                continue
            if ev > ALPHA_P or qc < C_P or tc < C_P:
                continue
            fwd[(q, t)] = pid
    edges = set()
    for (q, t) in fwd:
        if (t, q) in fwd:
            edges.add(frozenset((q, t)))
    return edges


def load_truth():
    """dna_loose set(frozenset), evaluable set(frozenset)=all rows present in truth file."""
    dna_loose = set()
    evaluable = set()
    with open(TRUTH) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            key = frozenset((f[0], f[1]))
            evaluable.add(key)
            if int(f[3]) == 1:      # in_dna_loose
                dna_loose.add(key)
    return dna_loose, evaluable


def load_biotype():
    """gene -> biotype (any interval)."""
    bt = {}
    with open(ANNOT) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            # contig start end biotype gene
            g = f[4]
            if g not in bt:
                bt[g] = f[3]
    return bt


# ------------------------------------------------------------------ symbol root
_root_re = re.compile(r'^([A-Za-z]+?[A-Za-z]*?)(\d+.*)?$')
def symbol_root(g):
    """Heuristic gene-symbol root: strip trailing digits+suffix. GOLGA6L7->GOLGA6L is kept as
    letters-before-first-digit-run? We use: letters up to the LAST alpha before a terminal digit
    block. Simpler robust rule: remove trailing run of digits and any trailing single letters that
    follow digits. Practically: take chars until the first digit that is followed only by digits/letters
    at the end. We approximate with: strip a terminal [0-9]+[A-Z]? and any preceding single-letter dup.
    For LOC ids there is no meaningful root -> return the whole id (each LOC becomes its own root)."""
    if g.startswith("LOC"):
        return g            # uncharacterized: no family root, self-root (conservative: won't inflate purity)
    # remove a terminal numeric+optional-letter suffix, iteratively for things like GOLGA6L7 -> GOLGA6L
    m = re.match(r'^(.*?)(\d+[A-Za-z]?)$', g)
    if m and m.group(1):
        return m.group(1).rstrip('-')
    return g


# ------------------------------------------------------------------ scoring
def within_family_pairs(fams):
    """catalog {fid:set(genes)} -> (multi-copy families list of gene-sets, all within-family pairs set)."""
    mc = [gs for gs in fams.values() if len(gs) >= 2]
    pairs = set()
    for gs in mc:
        for a, b in combinations(sorted(gs), 2):
            pairs.add(frozenset((a, b)))
    return mc, pairs


def ep_homologous(key, ep_direct, g2f, mega):
    if key in ep_direct:
        return True
    a, b = tuple(key)
    fa, fb = g2f.get(a), g2f.get(b)
    if fa is not None and fa == fb and fa not in mega:
        return True
    return False


def score_catalog(name, fams, ep_direct, g2f, mega, dna_loose, evaluable, biotype):
    mc, pairs = within_family_pairs(fams)
    n_fam_total = len(fams)
    n_fam_mc = len(mc)

    # ----- (a) Ep family-level purity (block over-merge) -----
    ep_evaluable_fams = 0      # families with >=2 genes carrying a non-mega Ep family
    ep_impure_fams = 0         # of those, spanning >=2 distinct non-mega Ep families
    ep_evaluable_incl_mega = 0
    ep_impure_incl_mega = 0
    for gs in mc:
        nonmega = {g2f[g] for g in gs if g in g2f and g2f[g] not in mega}
        allep = {g2f[g] for g in gs if g in g2f}
        # count genes mapping to a non-mega Ep fam
        n_nonmega_genes = sum(1 for g in gs if g in g2f and g2f[g] not in mega)
        n_ep_genes = sum(1 for g in gs if g in g2f)
        if n_nonmega_genes >= 2:
            ep_evaluable_fams += 1
            if len(nonmega) >= 2:
                ep_impure_fams += 1
        if n_ep_genes >= 2:
            ep_evaluable_incl_mega += 1
            if len(allep) >= 2:
                ep_impure_incl_mega += 1

    # ----- (a') Ep edge-level precision -----
    n_pairs = len(pairs)
    ep_hom = sum(1 for k in pairs if ep_homologous(k, ep_direct, g2f, mega))
    # restrict to pairs where BOTH genes carry an Ep family (evaluable-for-Ep)
    pairs_ep_eval = [k for k in pairs if all(g in g2f for g in k)]
    ep_hom_eval = sum(1 for k in pairs_ep_eval if ep_homologous(k, ep_direct, g2f, mega))

    # ----- (b) cDNA pair precision -----
    tp_cdna = sum(1 for k in pairs if k in dna_loose)
    pairs_in_truthfile = [k for k in pairs if k in evaluable]
    tp_cdna_eval = sum(1 for k in pairs_in_truthfile if k in dna_loose)

    # ----- (c) annotation symbol-root coherence -----
    single_root_fams = 0
    root_purity_sum = 0.0
    named_families = 0        # families with >=2 NAMED (non-LOC) genes
    named_single_root = 0
    named_purity_sum = 0.0
    for gs in mc:
        roots = [symbol_root(g) for g in gs]
        c = Counter(roots)
        dom = c.most_common(1)[0][1]
        purity = dom / len(gs)
        root_purity_sum += purity
        if len(c) == 1:
            single_root_fams += 1
        named = [g for g in gs if not g.startswith("LOC")]
        if len(named) >= 2:
            named_families += 1
            nr = Counter(symbol_root(g) for g in named)
            ndom = nr.most_common(1)[0][1]
            named_purity_sum += ndom / len(named)
            if len(nr) == 1:
                named_single_root += 1

    # ----- gene coverage -----
    all_genes = set()
    for gs in fams.values():
        all_genes |= gs

    return dict(
        name=name,
        n_families_total=n_fam_total,
        n_families_multicopy=n_fam_mc,
        n_distinct_genes=len(all_genes),
        n_within_family_pairs=n_pairs,
        # (a) Ep family purity
        ep_evaluable_families=ep_evaluable_fams,
        ep_impure_families=ep_impure_fams,
        ep_impure_rate=round(ep_impure_fams / max(1, ep_evaluable_fams), 4),
        ep_pure_rate=round(1 - ep_impure_fams / max(1, ep_evaluable_fams), 4),
        ep_evaluable_families_incl_mega=ep_evaluable_incl_mega,
        ep_impure_families_incl_mega=ep_impure_incl_mega,
        ep_impure_rate_incl_mega=round(ep_impure_incl_mega / max(1, ep_evaluable_incl_mega), 4),
        # (a') Ep edge precision
        ep_pair_prec_allpairs=round(ep_hom / max(1, n_pairs), 4),
        ep_hom_pairs=ep_hom,
        ep_evaluable_pairs=len(pairs_ep_eval),
        ep_pair_prec_evaluable=round(ep_hom_eval / max(1, len(pairs_ep_eval)), 4),
        ep_hom_pairs_evaluable=ep_hom_eval,
        # (b) cDNA
        cdna_tp_allpairs=tp_cdna,
        cdna_pair_prec_allpairs=round(tp_cdna / max(1, n_pairs), 4),
        cdna_pairs_in_truthfile=len(pairs_in_truthfile),
        cdna_tp_evaluable=tp_cdna_eval,
        cdna_pair_prec_evaluable=round(tp_cdna_eval / max(1, len(pairs_in_truthfile)), 4),
        # (c) annotation
        annot_single_root_rate=round(single_root_fams / max(1, n_fam_mc), 4),
        annot_root_purity_mean=round(root_purity_sum / max(1, n_fam_mc), 4),
        annot_named_families=named_families,
        annot_named_single_root_rate=round(named_single_root / max(1, named_families), 4),
        annot_named_root_purity_mean=round(named_purity_sum / max(1, named_families), 4),
        _pairs=pairs,          # for recall cross-analysis
        _genes=all_genes,
    )


def main():
    print("[load] catalogs / truths / Ep ...", flush=True)
    dna = load_dna_catalog()
    rna = load_rna_catalog()
    g2f, mega, sizes = load_prfam()
    ep_direct = load_ep_direct()
    dna_loose, evaluable = load_truth()
    biotype = load_biotype()
    print(f"       DNA_families={len(dna)}  RNA_families={len(rna)}  Ep_genes={len(g2f)}  "
          f"mega_Ep={sorted(mega)} (sizes={[sizes[m] for m in sorted(mega)]})  "
          f"ep_direct_edges={len(ep_direct)}  dna_loose_truth={len(dna_loose)}  "
          f"truth_rows={len(evaluable)}", flush=True)

    d = score_catalog("DNA(E_a segdup)", dna, ep_direct, g2f, mega, dna_loose, evaluable, biotype)
    r = score_catalog("RNA(E_r transcript)", rna, ep_direct, g2f, mega, dna_loose, evaluable, biotype)

    # ----- recall / completeness: paralog pairs (in_dna_loose) captured within-family -----
    dna_loose_genes = set()
    for k in dna_loose:
        dna_loose_genes |= set(k)
    def paralog_capture(res):
        pairs = res["_pairs"]
        genes = res["_genes"]
        captured = sum(1 for k in dna_loose if k in pairs)
        genes_cov = len(dna_loose_genes & genes)
        # of dna_loose pairs whose BOTH genes are present in this catalog's node set, how many co-grouped
        reachable = [k for k in dna_loose if all(g in genes for g in k)]
        cap_reach = sum(1 for k in reachable if k in pairs)
        return dict(paralog_pairs_captured=captured,
                    paralog_pair_recall_abs=round(captured / max(1, len(dna_loose)), 4),
                    dna_loose_genes_covered=genes_cov,
                    dna_loose_gene_coverage=round(genes_cov / max(1, len(dna_loose_genes)), 4),
                    reachable_pairs=len(reachable),
                    captured_of_reachable=cap_reach,
                    recall_on_reachable=round(cap_reach / max(1, len(reachable)), 4))
    d_rec = paralog_capture(d)
    r_rec = paralog_capture(r)

    # ----- Ep multi-copy family recall: of non-mega Ep fams with >=2 genes, how many have >=2 members grouped -----
    ep_fams = defaultdict(set)
    for g, fid in g2f.items():
        if fid not in mega:
            ep_fams[fid].add(g)
    ep_mc = {fid: gs for fid, gs in ep_fams.items() if len(gs) >= 2}
    def ep_family_recall(res):
        # for each catalog family, which Ep families does it place >=2 members into (co-group)?
        pairs = res["_pairs"]
        recovered = 0
        for fid, gs in ep_mc.items():
            # a catalog "recovers" this Ep family if >=2 of its genes share a within-family pair
            gl = sorted(gs)
            hit = any(frozenset((a, b)) in pairs for a, b in combinations(gl, 2))
            if hit:
                recovered += 1
        return dict(ep_mc_families=len(ep_mc), ep_mc_recovered=recovered,
                    ep_mc_recall=round(recovered / max(1, len(ep_mc)), 4))
    d_epr = ep_family_recall(d)
    r_epr = ep_family_recall(r)

    for res in (d, r):
        del res["_pairs"]; del res["_genes"]

    # ============================== report ==============================
    def row(label, key, fmt="{}"):
        return f"  {label:<46} DNA={fmt.format(d[key]):>12}   RNA={fmt.format(r[key]):>12}"
    print("\n============================================================")
    print("HEAD-TO-HEAD  DNA (E_a segdup)  vs  RNA (E_r transcript)")
    print("============================================================")
    print(row("families (total / multi-copy)",
              "n_families_total") + f"   [mc DNA={d['n_families_multicopy']} RNA={r['n_families_multicopy']}]")
    print(row("distinct genes", "n_distinct_genes"))
    print(row("within-family gene pairs", "n_within_family_pairs"))
    print("\n-- (a) PROTEIN E_p family purity (specificity) --")
    print(row("Ep-evaluable families (>=2 non-mega Ep genes)", "ep_evaluable_families"))
    print(row("Ep-IMPURE families (span >=2 Ep families)", "ep_impure_families"))
    print(row("Ep-impure RATE  (lower=more specific)", "ep_impure_rate"))
    print(row("Ep-PURE rate", "ep_pure_rate"))
    print("\n-- (a') PROTEIN E_p edge precision --")
    print(row("Ep-homolog pair precision (ALL pairs)", "ep_pair_prec_allpairs"))
    print(row("Ep-evaluable pairs (both genes in Ep)", "ep_evaluable_pairs"))
    print(row("Ep-homolog pair precision (EVALUABLE)", "ep_pair_prec_evaluable"))
    print("\n-- (b) cDNA 90% pair precision --")
    print(row("cDNA pair precision (ALL within-fam pairs)", "cdna_pair_prec_allpairs"))
    print(row("pairs present in truth file (evaluable)", "cdna_pairs_in_truthfile"))
    print(row("cDNA pair precision (EVALUABLE only)", "cdna_pair_prec_evaluable"))
    print("\n-- (c) annotation symbol-root coherence --")
    print(row("single-root family rate (all genes incl LOC)", "annot_single_root_rate"))
    print(row("mean dominant-root purity (all genes)", "annot_root_purity_mean"))
    print(row("named-gene families (>=2 non-LOC)", "annot_named_families"))
    print(row("single-root rate (NAMED only)", "annot_named_single_root_rate"))
    print(row("mean dominant-root purity (NAMED only)", "annot_named_root_purity_mean"))
    print("\n-- RECALL / COMPLETENESS (paralog copies from DNA/segdup cDNA-loose set) --")
    print(f"  cDNA-loose truth pairs={len(dna_loose)}  truth genes={len(dna_loose_genes)}")
    print(f"  {'paralog pairs captured within-family':<46} "
          f"DNA={d_rec['paralog_pairs_captured']:>12}   RNA={r_rec['paralog_pairs_captured']:>12}")
    print(f"  {'  = abs recall of cDNA-loose pairs':<46} "
          f"DNA={d_rec['paralog_pair_recall_abs']:>12}   RNA={r_rec['paralog_pair_recall_abs']:>12}")
    print(f"  {'cDNA-loose truth genes covered':<46} "
          f"DNA={d_rec['dna_loose_genes_covered']:>12}   RNA={r_rec['dna_loose_genes_covered']:>12}")
    print(f"  {'  = gene coverage':<46} "
          f"DNA={d_rec['dna_loose_gene_coverage']:>12}   RNA={r_rec['dna_loose_gene_coverage']:>12}")
    print(f"  {'reachable pairs (both genes in catalog)':<46} "
          f"DNA={d_rec['reachable_pairs']:>12}   RNA={r_rec['reachable_pairs']:>12}")
    print(f"  {'recall on reachable':<46} "
          f"DNA={d_rec['recall_on_reachable']:>12}   RNA={r_rec['recall_on_reachable']:>12}")
    print(f"  {'Ep multi-copy families recovered (>=2 grp)':<46} "
          f"DNA={d_epr['ep_mc_recovered']}/{d_epr['ep_mc_families']}={d_epr['ep_mc_recall']:>6}   "
          f"RNA={r_epr['ep_mc_recovered']}/{r_epr['ep_mc_families']}={r_epr['ep_mc_recall']:>6}")

    out = dict(dna=d, rna=r, dna_recall=d_rec, rna_recall=r_rec,
               dna_ep_family_recall=d_epr, rna_ep_family_recall=r_epr,
               truth=dict(dna_loose_pairs=len(dna_loose), dna_loose_genes=len(dna_loose_genes),
                          truth_rows=len(evaluable), ep_direct_edges=len(ep_direct),
                          mega_ep=sorted(mega)))
    with open(os.path.join(BENCH, "family_dna_vs_rna_headtohead.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print("\nwrote bench/family_dna_vs_rna_headtohead.json")


if __name__ == "__main__":
    main()
