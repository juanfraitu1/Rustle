#!/usr/bin/env python
"""TE/REPEAT-EXCLUSION gate for the reframed RNA family (E_r) edge oracle.

OPEN PROBLEM (bench/FAMILY_ER_ORACLE_SWEEP.md, commit 73a7120): the pure-protein gate cuts genuine
over-merges cleanly but drops ~97 protein-missed real TP (non-coding / LOC* paralogs); no core_recip
threshold separates real non-coding paralogs (KEEP) from TE-bridged junk (CUT) because BOTH are
non-protein-homologous with overlapping core_recip.

HYPOTHESIS: the soft-mask/TE status of the BRIDGING exonic sequence disambiguates them -- a TE-bridged
junk edge shares a REPEAT (masked) segment; a real non-coding paralog shares UNIQUE (unmasked) sequence.

METHOD (chosen in Verify): the production feature is the WHOLE-LOCUS MAX EXONIC SOFTMASK ("bridge_mask" in
the cache, but read it as a whole-locus proxy, NOT the shared-region mask) = max softmask (lowercase)
fraction over the two bridged de-novo loci's SPLICED (exon-only) sequence, taken over the argmax-core
DN--DN bridge edge, read from the soft-masked GGO.fasta. Spliced sequence drops intronic TEs (only
EXONIZED TEs matter at the RNA level). It is a whole-locus MAX (not the aligned/shared segment) because
divergent paralogs (ENO2-ENO3, OGDH-OGDHL) do NOT align at minimap2-asm20 at all, so an aligned-region
mask is UNDEFINED for exactly the pairs we care most about -> the whole-locus proxy is the justified
fallback, covers all genes (incl. LOC* with no annotated exon record), and needs no aligner -> the gate
is fully deterministic. An ALIGNED-REGION shared-mask CROSS-CHECK (minimap2 -cx asm20 -t1, deterministic)
is computed as a DIAGNOSTIC on named + asymmetric pairs to expose the MAX-inflation risk (a repeat-rich
partner inflating the edge mask even when the SHARED bridge is unique, e.g. FGF2-NUDT6 whole 0.223 but
shared region 0.000). It never touches the gate. Cached to bench/edge_bridge_mask.tsv.

COMBINED GATE:  E_r edge valid  <=>  in_ep  OR  (core_recip >= t  AND  bridge_mask < m)
  m = 1.0 reduces exactly to the existing 'combined' gate (in_ep OR core>=t) -> built-in sanity anchor
  (its rows must reproduce family_er_oracle_sweep 'combined' rows). Transitive-only pairs (no direct DN
  edge -> core=None -> bridge_mask=None) fail the CORE branch regardless -- the mask is never evaluated
  on them.

Sweeps t x m, and for every point emits (bench/family_er_te_gate.tsv): gate params, n_edges, TP retained,
tp_noep_retained (of the 97 protein-missed real TP -- the recovery target), noncoding_tp_retained (of 31),
truthbar_retained (the protein-homologous divergent-paralog proxy = the E_r WIN; TAUTOLOGICAL, rides the
in_ep branch, mask never touches it), genuine_overmerge + BOTH block faces (prfam-mixing + edge-host),
id-stratified reachable recall. JSON adds the curve + recommended point + mask_discrimination (class-mean
mask + AUC by core band) + an HONEST named confirmation built from the edges the MASK BRANCH ACTUALLY
DECIDES (genuine cut-by-mask, TP kept-by-mask, junk the mask can't catch) + the aligned-region cross-check
+ head-to-head vs pure-protein and combined-core(no-mask).

HONEST FRAMING (promoted from fine print to the stated conclusion, per adversarial review):
  (1) The hypothesis "real paralog = UNIQUE/unmasked, TE-junk = REPEAT/masked" is POPULATION-LEVEL FALSE
      for the real PROTEIN-HOMOLOGOUS paralogs: truthbar median mask 0.312 ~ genuine 0.346 (frac>=0.10:
      0.908 vs 0.926). Real paralogs are as repeat-rich as junk; they survive ONLY via the in_ep
      tautology, NOT the mask. The mask signal holds ONLY for the protein-MISSED TP subset (median 0.172).
  (2) The mask is a REAL but WEAK, NON-COMPLEMENTARY signal on the mask-decidable population (in_ep=0,
      core present): AUC(genuine>TP) = 0.777 over the whole re-admit region but DEGRADES to 0.694/0.675 at
      the operating band core>=0.50/0.70 -- it is NOT a clean TE disambiguator, it SHIFTS the operating
      point. At the recommended point it cuts 138 genuine per 12 protein-missed TP (11.5x, not the 19.2x
      sweep-aggregate).
  (3) The named showcase FPs AMY2A-ZNF91 / AMY2A-ZNF141 / RPL14-ZNF669 have core=None -> bridge_mask=None:
      they are cut by the CORE branch, NOT the mask. They are relabeled here as core-branch cuts and are
      NOT mask confirmations; the mask-decided demonstration uses edges that actually carry a bridge_mask.
  (4) The payoff is a MODEST Pareto SHIFT: +2/97 protein-missed TP (+1/31 non-coding) over the pure-core
      knee at LOWER edge-host over-merge; the open problem is shifted, not dissolved -- at m=0.10 the mask
      would still CUT 61% of the protein-missed TP it aims to recover. Both block faces are reported; the
      edge-host face is the arbiter.

Reproduce: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/family_er_te_gate.py
Deterministic: no aligner in the production path; block re-refine seed=0; hash seed pinned in-script.
"""
import json
import os
import statistics as _stat
import subprocess
import sys
import tempfile
from bisect import bisect_left, bisect_right
from collections import defaultdict, Counter

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)
import family_er_pr as F                 # loaders, projection
import family_er_oracle_sweep as S       # prepare(), compute_point(), predicate_for()
import genome_family_def as G            # refine_families (unused directly; via S)
import pysam

SCRATCH = "/home/juanfra/winloci_scratch"
GENOME = os.path.join(SCRATCH, "GGO.fasta")
SKEL = os.path.join(SCRATCH, "denovo_skeletons.tsv")
DENOVO_EDGES = os.path.join(BENCH, "denovo_family_edges.tsv")
MASK_CACHE = os.path.join(BENCH, "edge_bridge_mask.tsv")

T_GRID = [0.13, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80]
M_GRID = [0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 1.00]   # 1.00 = mask OFF (sanity: reduces to 'combined')
M_STAR = 0.10   # the Verify-chosen mask threshold (cut the core branch if bridge_mask >= 0.10)

# named confirmation pairs
NAMED_CUT = [("AMY2A", "ZNF91"), ("AMY2A", "ZNF141"), ("RPL14", "ZNF669")]      # TE-junk -> must be CUT
NAMED_KEEP_CODING = [("ENO2", "ENO3"), ("OGDH", "OGDHL"), ("GALNT14", "GALNT16")]  # real paralogs -> KEEP
NAMED_KEEP_NC = [("SURF2", "SURF4"), ("CTSA", "NEURL2")]                        # protein-missed real TP -> KEEP


# ---------------------------------------------------------------- exonic soft-mask (spliced sequence)
def load_skeletons():
    """(chrom,start,end) -> introns string. Keyed exactly as family_er_pr.load_meta values."""
    skel = {}
    with open(SKEL) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            skel[(f[0], int(f[1]), int(f[2]))] = f[5]
    return skel


def dn_exons(dn, meta, skel):
    """Spliced exon intervals for a DN locus = complement of the introns column within [start,end]."""
    mm = meta.get(dn)
    if mm is None:
        return None, None
    chrom, start, end = mm
    introns = skel.get((chrom, start, end))
    if introns is None:
        return chrom, [(start, end)]        # single-exon fallback (span)
    iv = []
    if introns.strip():
        for tok in introns.split(";"):
            a, b = tok.split("-")
            iv.append((int(a), int(b)))
    iv.sort()
    exons = []
    cur = start
    for a, b in iv:
        if a > cur:
            exons.append((cur, a))
        cur = b
    if end > cur:
        exons.append((cur, end))
    return chrom, exons


def build_bridge_mask(D, meta):
    """Per gene-pair edge -> (core, bridge_dn_a, bridge_dn_b, bridge_mask) where bridge = the argmax-core
    co-blocked DN--DN edge and bridge_mask = max exonic-softmask fraction of the two loci. This mirrors
    exactly how family_er_oracle_sweep builds gp_core (same co-blocked constraint), so the mask is
    attached to the SAME DN edge whose core the gate tests."""
    skel = load_skeletons()
    fa = pysam.FastaFile(GENOME)

    annot = F.load_annot()
    gene_of = F.gene_of_factory(annot)
    raw_fams = F.load_raw_families()
    genes, gene_of_dn, _, _, _ = F.build_genes_dict(raw_fams, meta, gene_of)
    dn2blk = {}
    for bi, blk in enumerate(D["blocks0"]):
        for dn in blk:
            dn2blk[dn] = bi

    # gene-pair -> (max_core, bridge_dn_a, bridge_dn_b) over direct co-blocked DN edges
    gp_bridge = {}
    with open(DENOVO_EDGES) as fh:
        next(fh)
        for ln in fh:
            a, b, cc = ln.rstrip("\n").split("\t")
            cc = float(cc)
            ga, gb = gene_of_dn.get(a), gene_of_dn.get(b)
            if ga is None or gb is None or ga == gb:
                continue
            if dn2blk.get(a) is None or dn2blk.get(a) != dn2blk.get(b):
                continue
            k = frozenset((ga, gb))
            if k not in gp_bridge or cc > gp_bridge[k][0]:
                gp_bridge[k] = (cc, a, b)

    _mc = {}
    def dn_mask(dn):
        if dn in _mc:
            return _mc[dn]
        chrom, exons = dn_exons(dn, meta, skel)
        if exons is None:
            _mc[dn] = None
            return None
        s = "".join(fa.fetch(chrom, a, b) for a, b in exons)
        v = None if not s else sum(1 for c in s if c.islower()) / len(s)
        _mc[dn] = v
        return v

    # bridge_mask per gene-pair (only pairs present in the refined set with a direct co-blocked DN edge)
    out = {}
    for k in D["ref_pairs"]:
        br = gp_bridge.get(k)
        if br is None:
            out[k] = dict(core=None, dn_a=None, dn_b=None, mask=None)
            continue
        cc, a, b = br
        ma, mb = dn_mask(a), dn_mask(b)
        mask = None
        if ma is not None and mb is not None:
            mask = max(ma, mb)
        elif ma is not None:
            mask = ma
        elif mb is not None:
            mask = mb
        out[k] = dict(core=cc, dn_a=a, dn_b=b, mask=mask, mask_a=ma, mask_b=mb)
    return out, gene_of_dn


def write_mask_cache(bridge, ref_pairs, attr):
    cols = ["gene_a", "gene_b", "cls", "in_ep", "core", "bridge_dn_a", "bridge_dn_b",
            "mask_a", "mask_b", "bridge_mask"]
    with open(MASK_CACHE, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for k in sorted(ref_pairs, key=lambda x: tuple(sorted(x))):
            a, b = sorted(k)
            br = bridge[k]
            def fm(x):
                return "" if x is None else f"{x:.4f}"
            fh.write("\t".join([a, b, attr[k]["cls"], str(int(attr[k]["in_ep"])),
                                fm(br["core"]), br["dn_a"] or "", br["dn_b"] or "",
                                fm(br.get("mask_a")), fm(br.get("mask_b")), fm(br["mask"])]) + "\n")


# ---------------------------------------------------------------- mask discrimination stats
def _median(xs):
    return round(_stat.median(xs), 4) if xs else None


def _auc_gt(pos, neg):
    """Mann-Whitney AUC = P(random pos mask > random neg mask), 0.5-tie-corrected. Here pos=genuine,
    neg=protein-missed TP, so AUC>0.5 means genuine tends to be MORE masked than real TP (the useful
    direction). 0.5 = no separation."""
    if not pos or not neg:
        return None
    negs = sorted(neg)
    c = 0.0
    for v in pos:
        lo = bisect_left(negs, v)
        hi = bisect_right(negs, v)
        c += lo + 0.5 * (hi - lo)
    return round(c / (len(pos) * len(neg)), 4)


def mask_discrimination(ref_pairs, attr):
    """Class-level mask distributions + AUC(genuine>TP) by core band on the MASK-DECIDABLE population
    (in_ep=0, direct core present). This is the honest evidence: (a) the hypothesis is population-FALSE
    for protein-homologous paralogs (truthbar mask ~ genuine mask), it holds ONLY for protein-missed TP;
    (b) AUC degrades from the whole-region 0.777 to ~0.68-0.69 at the operating band core>=0.5-0.7."""
    def stats(cls_filter):
        ms = [attr[k]["mask"] for k in ref_pairs if cls_filter(k) and attr[k]["mask"] is not None]
        if not ms:
            return dict(n=0, median=None, mean=None, frac_ge_0_10=None)
        return dict(n=len(ms), median=round(_stat.median(ms), 4), mean=round(_stat.mean(ms), 4),
                    frac_ge_0_10=round(sum(1 for m in ms if m >= 0.10) / len(ms), 4))
    by_class = dict(
        TP=stats(lambda k: attr[k]["cls"] == "TP"),
        TP_protein_missed=stats(lambda k: attr[k]["cls"] == "TP" and not attr[k]["in_ep"]),
        truthbar=stats(lambda k: attr[k]["cls"] == "truthbar"),
        genuine=stats(lambda k: attr[k]["cls"] == "genuine"))
    auc_band = {}
    for lo in (0.13, 0.30, 0.50, 0.70):
        gen = [attr[k]["mask"] for k in ref_pairs if attr[k]["cls"] == "genuine" and not attr[k]["in_ep"]
               and attr[k]["core"] is not None and attr[k]["core"] >= lo and attr[k]["mask"] is not None]
        tp = [attr[k]["mask"] for k in ref_pairs if attr[k]["cls"] == "TP" and not attr[k]["in_ep"]
              and attr[k]["core"] is not None and attr[k]["core"] >= lo and attr[k]["mask"] is not None]
        auc_band[f"core>={lo}"] = dict(n_genuine=len(gen), n_TP=len(tp), auc_genuine_gt_TP=_auc_gt(gen, tp),
                                       genuine_median=_median(gen), TP_median=_median(tp))
    return dict(
        by_class=by_class, auc_by_core_band=auc_band,
        conclusion=("HYPOTHESIS IS POPULATION-FALSE for real protein-homologous paralogs: truthbar median "
                    f"mask {by_class['truthbar']['median']} ~ genuine {by_class['genuine']['median']} "
                    f"(frac>=0.10 {by_class['truthbar']['frac_ge_0_10']} vs {by_class['genuine']['frac_ge_0_10']}); "
                    "they are saved ONLY by the in_ep tautology. The mask discriminates ONLY the "
                    f"protein-missed TP subset (median {by_class['TP_protein_missed']['median']}), and even "
                    "there AUC degrades from 0.777 (whole re-admit region) to ~0.68-0.69 at the operating "
                    "band core>=0.50-0.70 -- a WEAK, non-complementary signal, not a clean TE disambiguator."))


def mask_decided_examples(ref_pairs, attr, t, m, noncoding_tp):
    """The edges the MASK BRANCH ACTUALLY DECIDES at the recommended (t, m): in_ep=0 AND core>=t (so the
    core branch already admits them and ONLY the mask can still cut). This is the honest named
    confirmation -- NOT AMY2A-ZNF91 (core=None, cut by the core branch, mask never seen)."""
    pop = [k for k in ref_pairs if not attr[k]["in_ep"] and attr[k]["core"] is not None
           and attr[k]["core"] >= t and attr[k]["mask"] is not None]
    def nm(k):
        return " -- ".join(sorted(k))
    def rec(k):
        return dict(pair=nm(k), cls=attr[k]["cls"], core=round(attr[k]["core"], 3),
                    mask=round(attr[k]["mask"], 3), noncoding=(k in noncoding_tp))
    gen_cut = sorted([k for k in pop if attr[k]["cls"] == "genuine" and attr[k]["mask"] >= m],
                     key=lambda k: (-attr[k]["mask"], nm(k)))
    gen_surv = sorted([k for k in pop if attr[k]["cls"] == "genuine" and attr[k]["mask"] < m],
                      key=lambda k: (-attr[k]["core"], nm(k)))
    tp_keep = sorted([k for k in pop if attr[k]["cls"] == "TP" and attr[k]["mask"] < m],
                     key=lambda k: (-attr[k]["core"], nm(k)))
    tp_cut = sorted([k for k in pop if attr[k]["cls"] == "TP" and attr[k]["mask"] >= m],
                    key=lambda k: (-attr[k]["mask"], nm(k)))
    n_gen = len(gen_cut) + len(gen_surv)
    n_tp = len(tp_keep) + len(tp_cut)
    return dict(
        note=("edges with in_ep=0 AND core>=t -> admitted by the core branch, so ONLY the mask can cut "
              "them. THIS is where the mask acts; AMY2A-ZNF91/ZNF141/RPL14-ZNF669 are NOT here (core=None)."),
        band=f"in_ep=0, core>={t}", n_genuine=n_gen, n_TP=n_tp,
        genuine_cut_by_mask_n=len(gen_cut), genuine_surviving_mask_n=len(gen_surv),
        TP_kept_by_mask_n=len(tp_keep), TP_cut_by_mask_n=len(tp_cut),
        genuine_cut_rate=round(len(gen_cut) / n_gen, 4) if n_gen else None,
        TP_cut_rate=round(len(tp_cut) / n_tp, 4) if n_tp else None,
        per_capita_enrichment=(round((len(gen_cut) / n_gen) / (len(tp_cut) / n_tp), 2)
                               if n_gen and n_tp and len(tp_cut) else None),
        TE_junk_the_mask_CUTS=[rec(k) for k in gen_cut[:8]],       # high core, high mask -> genuine removed
        real_paralogs_the_mask_KEEPS=[rec(k) for k in tp_keep[:8]],  # high core, low mask -> real TP kept
        junk_the_mask_CANNOT_catch=[rec(k) for k in gen_surv[:8]],   # high core, low mask genuine -> residual
        real_TP_the_mask_WRONGLY_cuts=[rec(k) for k in tp_cut[:8]])  # high core, high mask real TP -> false cut


# ---------------------------------------------------------------- aligned-region shared-mask CROSS-CHECK
def _dn_spliced_seq(dn, meta, skel, fa):
    chrom, exons = dn_exons(dn, meta, skel)
    if exons is None:
        return None
    return "".join(fa.fetch(chrom, a, b) for a, b in exons)


def aligned_region_crosscheck(pairs, bridge, meta):
    """DIAGNOSTIC ONLY (never touches the gate): for each named/asymmetric pair align the two bridge DN
    transcripts (minimap2 -cx asm20 -t1, deterministic) and report the softmask fraction of the ALIGNED
    (shared) segment vs the whole-locus MAX proxy. Exposes MAX-inflation (repeat-rich partner inflating
    the edge mask when the SHARED region is unique) and confirms the proxy is UNDEFINED for divergent
    paralogs (they do not align at asm20 -> justifies the whole-locus fallback). Degrades gracefully if
    minimap2 is absent."""
    from shutil import which
    if which("minimap2") is None:
        return dict(available=False, note="minimap2 not on PATH; cross-check skipped")
    skel = load_skeletons()
    fa = pysam.FastaFile(GENOME)
    out = []
    for k in pairs:
        br = bridge.get(k)
        if br is None or br.get("dn_a") is None or br.get("dn_b") is None:
            out.append(dict(pair=" -- ".join(sorted(k)), whole_locus_max=None, aligned="no_bridge"))
            continue
        sa = _dn_spliced_seq(br["dn_a"], meta, skel, fa)
        sb = _dn_spliced_seq(br["dn_b"], meta, skel, fa)
        aln = None
        if sa and sb:
            with tempfile.TemporaryDirectory() as td:
                pa, pb = os.path.join(td, "a.fa"), os.path.join(td, "b.fa")
                with open(pa, "w") as fh:
                    fh.write(f">a\n{sa}\n")
                with open(pb, "w") as fh:
                    fh.write(f">b\n{sb}\n")
                r = subprocess.run(["minimap2", "-cx", "asm20", "-t", "1", pa, pb],
                                   capture_output=True, text=True)
            best = None
            for ln in r.stdout.splitlines():
                f = ln.split("\t")
                if len(f) < 12:
                    continue
                qs, qe, ts, te, ml = int(f[2]), int(f[3]), int(f[7]), int(f[8]), int(f[9])
                if best is None or ml > best[0]:
                    best = (ml, qs, qe, ts, te)
            if best is not None:
                ml, qs, qe, ts, te = best
                def lf(s):
                    return (round(sum(1 for c in s if c.islower()) / len(s), 4) if s else None)
                aln = dict(aln_len=ml, mask_query_side=lf(sb[qs:qe]), mask_target_side=lf(sa[ts:te]),
                           mask_shared=max(lf(sb[qs:qe]) or 0.0, lf(sa[ts:te]) or 0.0))
        rec = dict(pair=" -- ".join(sorted(k)), cls=None, whole_locus_max=br.get("mask"),
                   mask_a=br.get("mask_a"), mask_b=br.get("mask_b"),
                   aligned=(aln if aln is not None else "DID_NOT_ALIGN_asm20 (divergent -> proxy justified)"))
        out.append(rec)
    return dict(available=True,
                note=("aligned_region shared-mask cross-check (minimap2 -cx asm20 -t1). 'DID_NOT_ALIGN' = "
                      "divergent paralog, aligned mask UNDEFINED -> whole-locus MAX proxy is the fallback. "
                      "Where whole_locus_max >> mask_shared the proxy OVER-inflates via a repeat-rich "
                      "partner (MAX-inflation); where mask_shared >= whole_locus_max the shared bridge is "
                      "genuinely repeat."),
                pairs=out)


# ---------------------------------------------------------------- gate predicates
def kp_combined(t):
    return lambda k, a, t=t: a["in_ep"] or (a["core"] is not None and a["core"] >= t)


def kp_masked(t, m):
    # in_ep OR (core>=t AND bridge unmasked). m=1.0 -> mask<1.0 always true for real fracs -> == combined
    return lambda k, a, t=t, m=m: a["in_ep"] or (
        a["core"] is not None and a["core"] >= t
        and a["mask"] is not None and a["mask"] < m)


def name_status(D, kept_sets, pairs):
    """For each named pair: is it in the refined set, its cls/in_ep/core/mask, and kept at each gate."""
    attr = D["attr"]
    ref = D["ref_pairs"]
    res = []
    for pa, pb in pairs:
        k = frozenset((pa, pb))
        rec = dict(pair=f"{pa}--{pb}", in_refined=(k in ref))
        if k in ref:
            rec.update(cls=attr[k]["cls"], in_ep=bool(attr[k]["in_ep"]),
                       core=attr[k]["core"], mask=attr[k]["mask"])
            for gname, kept in kept_sets.items():
                rec[f"kept_{gname}"] = (k in kept)
        res.append(rec)
    return res


def main():
    print("[load] preparing data (reusing family_er_pr + family_er_oracle_sweep machinery) ...", flush=True)
    D = S.prepare()
    ref_pairs, attr = D["ref_pairs"], D["attr"]
    meta = F.load_meta()

    n_tp = sum(1 for k in ref_pairs if attr[k]["cls"] == "TP")
    n_tb = sum(1 for k in ref_pairs if attr[k]["cls"] == "truthbar")
    n_ge = sum(1 for k in ref_pairs if attr[k]["cls"] == "genuine")
    n_tp_noep = sum(1 for k in ref_pairs if attr[k]["cls"] == "TP" and not attr[k]["in_ep"])
    n_nc_tp = len(D["noncoding_tp"])
    print(f"[E_r] refined pairs={len(ref_pairs)} TP={n_tp} truthbar={n_tb} genuine={n_ge} "
          f"TP(protein-missed)={n_tp_noep} noncoding_TP={n_nc_tp}", flush=True)

    print("[mask] computing per-edge bridge_mask (exonic softmask MAX over argmax-core DN bridge) ...",
          flush=True)
    bridge, gene_of_dn = build_bridge_mask(D, meta)
    # attach mask to attr (compute_point predicates read attr)
    for k in ref_pairs:
        attr[k]["mask"] = bridge[k]["mask"]
    write_mask_cache(bridge, ref_pairs, attr)
    n_have_core = sum(1 for k in ref_pairs if attr[k]["core"] is not None)
    n_have_mask = sum(1 for k in ref_pairs if attr[k]["mask"] is not None)
    print(f"[mask] pairs with direct-core={n_have_core}  with bridge_mask={n_have_mask}  "
          f"(cached -> {MASK_CACHE})", flush=True)

    mask_disc = mask_discrimination(ref_pairs, attr)
    print("[mask-disc] " + mask_disc["conclusion"], flush=True)

    # ------------------------------------------------ operating points
    rows = []
    kept_sets = {}   # gate label -> kept set (for named confirmation)

    def add(gate, t, m, keep):
        r = S.compute_point(gate, t, keep, D)
        r["m"] = None if m is None else round(m, 3)
        r["gate"] = gate
        rows.append(r)
        kept_sets[f"{gate}_t{t if t is not None else 'NA'}_m{m if m is not None else 'NA'}"] = \
            {k for k in ref_pairs if keep(k, attr[k])}
        return r

    # references
    none_r = add("none", None, None, lambda k, a: True)
    prot_r = add("protein", None, None, lambda k, a: a["in_ep"])

    # sanity anchor vs the committed sweep
    assert (none_r["n_edges"], none_r["tp_retained"], none_r["truthbar_retained"],
            none_r["genuine_overmerge"]) == (10755, 450, 6142, 4163), f"none mismatch: {none_r}"
    assert (prot_r["n_edges"], prot_r["tp_retained"], prot_r["truthbar_retained"],
            prot_r["genuine_overmerge"]) == (6506, 353, 6142, 11), f"protein mismatch: {prot_r}"
    print("[sanity] none(10755/450/6142/4163) protein(6506/353/6142/11) reproduced  OK", flush=True)

    # combined (no mask) reference curve
    combined_rows = {}
    for t in T_GRID:
        combined_rows[t] = add("combined", t, None, kp_combined(t))

    # combined_mask sweep t x m  (m=1.0 must reproduce combined)
    mask_rows = {}
    for t in T_GRID:
        for m in M_GRID:
            r = add("combined_mask", t, m, kp_masked(t, m))
            mask_rows[(t, m)] = r

    # sanity: m=1.0 reduces to combined
    for t in T_GRID:
        c, mm = combined_rows[t], mask_rows[(t, 1.00)]
        assert (c["n_edges"], c["tp_retained"], c["genuine_overmerge"]) == \
               (mm["n_edges"], mm["tp_retained"], mm["genuine_overmerge"]), \
               f"m=1.0 != combined at t={t}: {c} vs {mm}"
    print("[sanity] combined_mask(m=1.0) == combined at all t  OK (mask OFF anchor)", flush=True)

    # ------------------------------------------------ recommended point
    # The core-only knee (family_er_oracle_sweep recommendation): combined t=0.70.
    core_knee = combined_rows[0.70]

    # PRIMARY RECOMMENDATION: hold the Verify-chosen mask threshold m*=0.10 (keeps the named non-coding
    # paralogs SURF2-SURF4 / CTSA-NEURL2, whose exonic mask ~0.07/0.01 < 0.10) and pick the t that
    # RECOVERS THE MOST protein-missed real TP subject to the arbiter over-merge face (edge-host) staying
    # AT OR BELOW the core knee. This is the honest "recover more at no-worse over-merge" point.
    cand = [r for r in rows if r["gate"] == "combined_mask" and r["m"] == round(M_STAR, 3)
            and r["genuine_rate_block_edgehost"] is not None
            and r["genuine_rate_block_edgehost"] <= core_knee["genuine_rate_block_edgehost"]
            and r["tp_noep_retained"] >= core_knee["tp_noep_retained"]]
    if cand:
        recommended = max(cand, key=lambda r: (r["tp_noep_retained"],
                                               -r["genuine_rate_block_edgehost"],
                                               -r["genuine_overmerge"]))
        rec_kind = "mstar0.10_maxrecovery_edgehost<=core_knee"
    else:
        recommended = mask_rows[(0.50, M_STAR)]
        rec_kind = "mstar0.10_fallback"

    # STRICT-PARETO ALTERNATIVE: the combined_mask point that dominates the core knee on protein-missed-TP
    # recovery AND on BOTH block faces (prfam + edge-host). More aggressive masking; may sacrifice a few
    # real non-coding TP with moderate exonic repeat (e.g. SURF2-SURF4 at mask 0.068 is cut for m<0.068).
    strict = [r for r in rows if r["gate"] == "combined_mask"
              and r["tp_noep_retained"] >= core_knee["tp_noep_retained"]
              and r["genuine_rate_block"] is not None
              and r["genuine_rate_block"] <= core_knee["genuine_rate_block"]
              and r["genuine_rate_block_edgehost"] is not None
              and r["genuine_rate_block_edgehost"] <= core_knee["genuine_rate_block_edgehost"]]
    strict_pareto = (max(strict, key=lambda r: (r["tp_noep_retained"], -r["genuine_overmerge"],
                                                -r["genuine_rate_block_edgehost"])) if strict else None)

    # ------------------------------------------------ does the mask point beat the ENTIRE no-mask curve?
    # The honest test of "is this a real gain from the mask": is there ANY no-mask combined point that
    # matches the recommended point's protein-missed-TP recovery AT NO WORSE edge-host over-merge? If NONE,
    # the recommended point lies strictly OFF the no-mask curve (a genuine Pareto shift, not reachable by
    # retuning t alone).
    rec_tppm = recommended["tp_noep_retained"]
    rec_eh = recommended["genuine_rate_block_edgehost"]
    no_mask_dominators = [dict(t=combined_rows[t]["t"], TPpm=combined_rows[t]["tp_noep_retained"],
                               edgehost=combined_rows[t]["genuine_rate_block_edgehost"])
                          for t in T_GRID
                          if combined_rows[t]["tp_noep_retained"] >= rec_tppm
                          and combined_rows[t]["genuine_rate_block_edgehost"] <= rec_eh]
    beats_no_mask = dict(
        recommended_TPpm=rec_tppm, recommended_edgehost=rec_eh,
        no_mask_curve=[dict(t=combined_rows[t]["t"], TPpm=combined_rows[t]["tp_noep_retained"],
                            edgehost=combined_rows[t]["genuine_rate_block_edgehost"]) for t in T_GRID],
        no_mask_points_matching_or_better=no_mask_dominators,
        recommended_lies_off_no_mask_curve=(len(no_mask_dominators) == 0),
        verdict=("NO no-mask combined point reaches TPpm>=%d at edge-host<=%.4f -> the recommended mask "
                 "point lies strictly OFF the no-mask curve (genuine, if modest, Pareto shift)."
                 % (rec_tppm, rec_eh) if not no_mask_dominators else
                 "a no-mask point already dominates -> the mask adds nothing here."))

    # ------------------------------------------------ head-to-head: mask vs no-mask at matched m*
    head2head = []
    for t in T_GRID:
        c = combined_rows[t]
        mm = mask_rows[(t, M_STAR)]
        head2head.append(dict(
            t=t,
            combined=dict(n=c["n_edges"], TP=c["tp_retained"], TPpm=c["tp_noep_retained"],
                          ncTP=c["noncoding_tp_retained"], genuine=c["genuine_overmerge"],
                          blk_prfam=c["genuine_rate_block"], blk_edgehost=c["genuine_rate_block_edgehost"]),
            combined_mask=dict(n=mm["n_edges"], TP=mm["tp_retained"], TPpm=mm["tp_noep_retained"],
                               ncTP=mm["noncoding_tp_retained"], genuine=mm["genuine_overmerge"],
                               blk_prfam=mm["genuine_rate_block"],
                               blk_edgehost=mm["genuine_rate_block_edgehost"]),
            d_genuine=mm["genuine_overmerge"] - c["genuine_overmerge"],
            d_TPpm=mm["tp_noep_retained"] - c["tp_noep_retained"],
            d_ncTP=mm["noncoding_tp_retained"] - c["noncoding_tp_retained"],
            d_truthbar=mm["truthbar_retained"] - c["truthbar_retained"],
            # the separation ratio: fraction of genuine cut vs fraction of protein-missed-TP cut by the mask
            genuine_cut=c["genuine_overmerge"] - mm["genuine_overmerge"],
            tppm_cut=c["tp_noep_retained"] - mm["tp_noep_retained"]))

    # separation summary at m* over the whole sweep: total genuine cut vs total TPpm cut by adding mask
    tot_gen_cut = sum(h["genuine_cut"] for h in head2head)
    tot_tppm_cut = sum(h["tppm_cut"] for h in head2head)
    sep_ratio = (tot_gen_cut / tot_tppm_cut) if tot_tppm_cut else None

    # separation AT THE RECOMMENDED point specifically (NOT the sweep aggregate) + genuine cut ATTRIBUTION
    # (which branch removes the over-merges: the CORE branch does the bulk, the mask only trims the residual).
    raw_t = next((tt for tt in T_GRID if round(tt, 3) == recommended["t"]), None)
    sep_at_rec = None
    if raw_t is not None:
        c_at = combined_rows[raw_t]                 # mask OFF at the recommended t
        m_at = mask_rows[(raw_t, M_STAR)]           # mask ON  at the recommended t
        base_gen = none_r["genuine_overmerge"]
        core_branch_cut = base_gen - c_at["genuine_overmerge"]      # genuine removed by core>=t alone
        mask_branch_cut = c_at["genuine_overmerge"] - m_at["genuine_overmerge"]  # additionally removed by mask
        tppm_cut_here = c_at["tp_noep_retained"] - m_at["tp_noep_retained"]
        sep_at_rec = dict(
            t=recommended["t"], m=recommended["m"],
            genuine_cut_by_mask=mask_branch_cut, tppm_cut_by_mask=tppm_cut_here,
            genuine_to_tppm_ratio_at_recommended=(round(mask_branch_cut / tppm_cut_here, 2)
                                                  if tppm_cut_here else None),
            attribution=dict(
                baseline_genuine=base_gen,
                cut_by_core_branch=core_branch_cut,
                cut_by_mask_branch=mask_branch_cut,
                surviving_genuine=m_at["genuine_overmerge"],
                note=("the CORE branch (in_ep OR core>=t) does the BULK of the over-merge removal "
                      f"({core_branch_cut} genuine edges, mostly transitive core=None or low-core); the "
                      f"MASK adds only {mask_branch_cut} more. Do NOT credit the block-rate headline to the "
                      "mask -- it is a residual trim on the core-branch survivors.")),
            note=("AT the recommended point (fixed t), turning the mask ON removes "
                  f"{mask_branch_cut} genuine per {tppm_cut_here} protein-missed TP; this is the honest "
                  "operating-point separation, NOT the 19.2x sweep-sum. It is a PRECISION KNOB: at fixed t "
                  "it REDUCES recall (drops TP) -- the +2 TPpm net gain vs the core knee comes from being "
                  "able to LOWER t (0.70->0.50), enabled by the mask's genuine-trim, not from mask "
                  "'recovering' any TP."))

    # ------------------------------------------------ named confirmation (HONEST: separated by what
    # actually decides each pair -- core branch vs mask branch vs protein branch)
    rec_keep = {k for k in ref_pairs
                if kp_masked(recommended["t"], recommended["m"])(k, attr[k])} if recommended["t"] is not None \
        else set(ref_pairs)
    named_sets = {"recommended": rec_keep,
                  "protein": kept_sets[f"protein_tNA_mNA"],
                  "combined_core_knee": kept_sets[f"combined_t0.7_mNA"]}

    def label_branch(k):
        """Which branch decides this pair at the recommended gate."""
        a = attr[k]
        if a["in_ep"]:
            return "PROTEIN (in_ep=1; mask never evaluated -- decorative)"
        if a["core"] is None:
            return "CORE=None (transitive-only, no bridge; cut by the core branch -- NOT a mask decision)"
        if a["core"] < recommended["t"]:
            return f"CORE<t (core {a['core']:.3f} < {recommended['t']}; cut by the core branch -- NOT a mask decision)"
        if a["mask"] is None:
            return f"MASK=None (core {a['core']:.3f}>=t but no mask -> CUT by the mask branch)"
        return f"MASK (core {a['core']:.3f}>=t; decided by mask {a['mask']:.3f} {'<' if a['mask'] < recommended['m'] else '>='} {recommended['m']})"

    def annotate(status_list):
        for rec in status_list:
            pa, pb = rec["pair"].split("--")          # name_status uses 'A--B' (no spaces)
            kk = frozenset((pa, pb))
            if kk in ref_pairs:
                rec["decided_by"] = label_branch(kk)
        return status_list

    named = dict(
        # RELABELED: these have core=None -> cut by the CORE branch, NOT the mask (NOT a mask confirmation).
        te_junk_cut_by_CORE_branch_not_mask=annotate(name_status(D, named_sets, NAMED_CUT)),
        # RELABELED: in_ep=1 -> kept by PROTEIN regardless of mask (mask decorative, does not test the gate).
        coding_paralog_kept_by_PROTEIN_not_mask=annotate(name_status(D, named_sets, NAMED_KEEP_CODING)),
        # the ONLY named pairs the mask actually decides (in_ep=0, core>=t, low mask -> mask KEEPS them):
        noncoding_tp_kept_by_mask=annotate(name_status(D, named_sets, NAMED_KEEP_NC)),
        note_named=("AMY2A-ZNF91/ZNF141/RPL14-ZNF669 do NOT verify the TE hypothesis -- they have "
                    "core=None/mask=None and are cut by the core branch. ENO2-ENO3/OGDH-OGDHL/GALNT14-"
                    "GALNT16 are in_ep=1, kept by protein regardless of mask. The only NAMED pairs the mask "
                    "decides are SURF2-SURF4 and CTSA-NEURL2 (both low-mask KEEP). For a genuine "
                    "cut-vs-keep-ON-THE-MASK demonstration see mask_decided_confirmation below."),
        # the HONEST confirmation: edges the mask branch actually decides at the recommended (t, m)
        mask_decided_confirmation=mask_decided_examples(ref_pairs, attr, recommended["t"],
                                                        recommended["m"], D["noncoding_tp"]))

    # hard confirmations (each labeled by WHICH branch does the cut/keep, to avoid the old
    # mask-misattribution: AMY2A-ZNF91 is cut by CORE=None, NOT the mask)
    amy_znf91 = frozenset(("AMY2A", "ZNF91"))
    eno = frozenset(("ENO2", "ENO3"))
    surf = frozenset(("SURF2", "SURF4"))
    conf = {}
    if amy_znf91 in ref_pairs:
        conf["AMY2A_ZNF91_cut_at_recommended"] = (amy_znf91 not in rec_keep)
        conf["AMY2A_ZNF91_cut_BY"] = ("CORE branch (core=None, no bridge_mask) -- NOT a mask cut"
                                      if attr[amy_znf91]["core"] is None else "core<t or mask")
        assert amy_znf91 not in rec_keep, "AMY2A-ZNF91 NOT cut at recommended point!"
        assert attr[amy_znf91]["core"] is None, "AMY2A-ZNF91 expected core=None (transitive-only)"
    else:
        conf["AMY2A_ZNF91_in_refined"] = False
    if eno in ref_pairs:
        conf["ENO2_ENO3_kept_at_recommended"] = (eno in rec_keep)
        conf["ENO2_ENO3_kept_BY"] = "PROTEIN branch (in_ep=1); mask decorative" if attr[eno]["in_ep"] else "core/mask"
    if surf in ref_pairs:
        # the ONE named pair the MASK actually keeps: in_ep=0, core>=t, mask<m
        conf["SURF2_SURF4_kept_at_recommended"] = (surf in rec_keep)
        conf["SURF2_SURF4_mask_decided"] = (not attr[surf]["in_ep"] and attr[surf]["core"] is not None
                                            and attr[surf]["core"] >= recommended["t"])
    # non-coding recovered = the noncoding_tp kept at recommended
    nc_kept = sorted(" -- ".join(sorted(k)) for k in D["noncoding_tp"] if k in rec_keep)
    nc_lost = sorted(" -- ".join(sorted(k)) for k in D["noncoding_tp"] if k not in rec_keep)

    # ------------------------------------------------ MAX-inflation set + aligned-region cross-check
    # protein-missed TP whose two loci disagree strongly in mask (|mask_a-mask_b|>0.15): the whole-locus
    # MAX may be driven by ONE repeat-rich partner, not a shared TE. Cross-check with the aligned (shared)
    # region mask. Also cross-check the named pairs.
    def _absdiff(k):
        ma, mb = bridge[k].get("mask_a"), bridge[k].get("mask_b")
        return abs(ma - mb) if (ma is not None and mb is not None) else 0.0
    max_infl_pairs = sorted([k for k in ref_pairs if attr[k]["cls"] == "TP" and not attr[k]["in_ep"]
                             and _absdiff(k) > 0.15], key=lambda k: (-_absdiff(k), tuple(sorted(k))))
    named_all = [frozenset(p) for p in (NAMED_CUT + NAMED_KEEP_CODING + NAMED_KEEP_NC)]
    crosscheck_pairs = [k for k in named_all if k in ref_pairs] + list(max_infl_pairs)
    aligned = aligned_region_crosscheck(crosscheck_pairs, bridge, meta)
    # attach cls to each cross-check record + summarize the MAX-inflation direction
    cls_by_label = {" -- ".join(sorted(k)): attr[k]["cls"] for k in crosscheck_pairs}
    for rec in aligned.get("pairs", []):
        rec["cls"] = cls_by_label.get(rec["pair"])
    max_inflation = dict(
        n_asymmetric_protein_missed_TP=len(max_infl_pairs),
        pairs=[dict(pair=" -- ".join(sorted(k)), mask_a=bridge[k]["mask_a"], mask_b=bridge[k]["mask_b"],
                    whole_locus_max=bridge[k]["mask"], core=round(attr[k]["core"], 3)) for k in max_infl_pairs],
        note=("for these the whole-locus MAX proxy is driven by ONE repeat-rich partner. The aligned "
              "cross-check resolves direction: FGF2-NUDT6 whole 0.223 but shared region 0.000 (unique -> "
              "proxy OVER-inflates, would be WRONGLY cut at m=0.10) vs ASB6-NTMT1 whole 0.344 shared 0.445 "
              "(genuinely repeat-bridged). A shared-region (or min) mask would fix the FGF2-NUDT6-type "
              "false cut, at the cost of being undefined for divergent paralogs that do not align."))

    # ------------------------------------------------ determinism (recompute recommended twice)
    r2 = S.compute_point(recommended["gate"], recommended["t"],
                         kp_masked(recommended["t"], recommended["m"]), D)
    det = (r2["n_edges"] == recommended["n_edges"]
           and r2["genuine_rate_block"] == recommended["genuine_rate_block"]
           and r2["genuine_rate_block_edgehost"] == recommended["genuine_rate_block_edgehost"]
           and r2["n_blocks"] == recommended["n_blocks"])

    # ============================== write TSV ==============================
    cols = ["gate", "t", "m", "n_edges", "tp_retained", "tp_noep_retained", "noncoding_tp_retained",
            "noncoding_tp_total", "truthbar_retained", "genuine_overmerge", "genuine_rate_edge",
            "genuine_prfam_invisible", "n_blocks", "genuine_rate_block", "genuine_rate_block_edgehost",
            "reach_recall_overall", "recall_id_0.90-0.95", "recall_id_0.95-0.99", "recall_id_>=0.99"]
    with open(os.path.join(BENCH, "family_er_te_gate.tsv"), "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")

    # ============================== write JSON ==============================
    summary = dict(
        method=dict(
            feature="WHOLE-LOCUS MAX exonic softmask (cache column 'bridge_mask', but it is NOT the "
                    "shared-region mask): max softmask fraction over the two loci of the argmax-core "
                    "co-blocked DN--DN bridge edge, on spliced exon-only seq from soft-masked GGO.fasta. "
                    "Whole-locus (not aligned/shared) because divergent paralogs (ENO2-ENO3, OGDH-OGDHL) do "
                    "NOT align at minimap2-asm20 -> aligned mask UNDEFINED for the pairs we care about; the "
                    "whole-locus proxy covers all genes and needs no aligner (deterministic gate). An "
                    "aligned-region shared-mask cross-check is computed as a DIAGNOSTIC (see "
                    "aligned_region_crosscheck) and exposes MAX-inflation on asymmetric pairs.",
            gate="in_ep OR (core_recip >= t AND bridge_mask < m)",
            mstar=M_STAR, t_grid=T_GRID, m_grid=M_GRID,
            deterministic="no aligner in the GATE path; block re-refine seed=0; PYTHONHASHSEED=0. The "
                          "aligned cross-check uses minimap2 -t1 (also deterministic) but never feeds the gate.",
            mask_off_anchor="m=1.0 reduces exactly to the 'combined' gate (verified per t)"),
        baseline=dict(refined_pairs=len(ref_pairs), TP=n_tp, truthbar=n_tb, genuine=n_ge,
                      TP_protein_missed=n_tp_noep, noncoding_TP=n_nc_tp,
                      pairs_with_direct_core=n_have_core, pairs_with_bridge_mask=n_have_mask),
        reference_points=dict(
            none=none_r, protein=prot_r, combined_core_knee_t070=core_knee),
        curve=rows,
        mask_discrimination=mask_disc,
        head_to_head_mstar=dict(
            note=f"combined (no mask) vs combined_mask at m={M_STAR}, per t. genuine_cut/tppm_cut = how many "
                 f"genuine vs protein-missed TP the mask branch removes at that t.",
            rows=head2head,
            separation_across_sweep=dict(
                total_genuine_cut_by_mask=tot_gen_cut,
                total_tppm_cut_by_mask=tot_tppm_cut,
                genuine_to_tppm_cut_ratio=(round(sep_ratio, 2) if sep_ratio else None),
                caveat="THIS IS A SUM OVER ALL t (double-counts nested cuts); it is NOT the operating-point "
                       "separation. Use separation_at_recommended for the honest figure."),
            separation_at_recommended=sep_at_rec),
        beats_no_mask_sweep=beats_no_mask,
        recommended=dict(kind=rec_kind, gate=recommended["gate"], t=recommended["t"], m=recommended["m"],
                         n_edges=recommended["n_edges"], tp_retained=recommended["tp_retained"],
                         tp_noep_retained=recommended["tp_noep_retained"],
                         noncoding_tp_retained=recommended["noncoding_tp_retained"],
                         noncoding_tp_total=n_nc_tp, truthbar_retained=recommended["truthbar_retained"],
                         genuine_overmerge=recommended["genuine_overmerge"],
                         genuine_rate_edge=recommended["genuine_rate_edge"],
                         genuine_rate_block=recommended["genuine_rate_block"],
                         genuine_rate_block_edgehost=recommended["genuine_rate_block_edgehost"],
                         vs_core_knee=dict(
                             core_knee_t=0.70,
                             d_tp_noep=recommended["tp_noep_retained"] - core_knee["tp_noep_retained"],
                             d_genuine=recommended["genuine_overmerge"] - core_knee["genuine_overmerge"],
                             d_block_prfam=round(recommended["genuine_rate_block"]
                                                 - core_knee["genuine_rate_block"], 4),
                             d_block_edgehost=round(recommended["genuine_rate_block_edgehost"]
                                                    - core_knee["genuine_rate_block_edgehost"], 4))),
        recommended_note=(
            "REFRAMED (per adversarial review): the mask is a PRECISION KNOB, not a recovery mechanism. At "
            "FIXED t it REDUCES recall (drops real TP alongside junk). The +2/97 protein-missed TP net gain "
            "over the pure-core knee (combined t=0.70) is a LOWER-t effect: because the mask trims genuine "
            "harder than TP, we can drop t from 0.70 to 0.50 and still hold the edge-host over-merge BELOW "
            "the core knee (8.9% vs 11.2%), and the lower t admits +2 more TP. Mask does NOT 'recover' any "
            "TP by itself. Truthbar (protein-homologous divergent-paralog proxy) is 100% retained but that "
            "is TAUTOLOGICAL -- it rides in_ep, the mask never touches it. Net vs core knee: +2 protein-"
            "missed TP (+1 non-coding), -31 genuine over-merge edges, -2.3pt edge-host block-rate."),
        strict_pareto_alternative=(None if strict_pareto is None else dict(
            note="dominates the core knee on protein-missed-TP recovery AND BOTH block faces; more "
                 "aggressive masking (m<=0.05) -> sacrifices some moderate-repeat real non-coding TP "
                 "(e.g. SURF2-SURF4, mask 0.068).",
            gate=strict_pareto["gate"], t=strict_pareto["t"], m=strict_pareto["m"],
            tp_retained=strict_pareto["tp_retained"], tp_noep_retained=strict_pareto["tp_noep_retained"],
            noncoding_tp_retained=strict_pareto["noncoding_tp_retained"],
            genuine_overmerge=strict_pareto["genuine_overmerge"],
            genuine_rate_block=strict_pareto["genuine_rate_block"],
            genuine_rate_block_edgehost=strict_pareto["genuine_rate_block_edgehost"])),
        named_confirmation=named,
        aligned_region_crosscheck=aligned,
        max_inflation=max_inflation,
        confirmations=conf,
        noncoding_recovered_at_recommended=nc_kept,
        noncoding_lost_at_recommended=nc_lost,
        honest_residual=(
            "PROMOTED FROM FINE PRINT (per adversarial review) -- three residual facts: "
            "(1) The hypothesis 'real paralog=unique/unmasked, TE-junk=repeat/masked' is POPULATION-FALSE "
            "for the real PROTEIN-HOMOLOGOUS paralogs: truthbar median mask 0.312 ~ genuine 0.346, so mask "
            "CANNOT separate them; they survive ONLY via the in_ep tautology. The signal exists ONLY on the "
            "protein-missed TP subset (median 0.172), and even there m=0.10 CUTS 61% of the very TP it aims "
            "to recover (real LOC paralogs are themselves repeat-derived). "
            "(2) On the mask-decidable population the signal is WEAK and DEGRADES at the operating band: "
            "AUC(genuine>TP) = 0.777 over the whole re-admit region but 0.694/0.675 at core>=0.50/0.70. So "
            "the mask SHIFTS the operating point (a modest Pareto move) but is NOT a clean TE disambiguator. "
            "(3) Residual TE-junk SURVIVES: near-identical unique-sequence LOC-family over-merges (e.g. "
            "LOC129530205--LOC129530238, core 1.0 mask 0.035) have low mask and the mask cannot touch them; "
            "and ~78% of the protein-missed TP stay dropped. INTRONIC-TE CAVEAT: the mask is computed on "
            "SPLICED exon sequence only, so intronic/unexonized TEs are (correctly) invisible -- a TE that "
            "bridges two loci only via a shared intron would not be caught, but such a bridge also would not "
            "survive the exonic core homology, so this is not a live failure mode here."),
        determinism=dict(recommended_reproducible=det, pythonhashseed=os.environ.get("PYTHONHASHSEED")),
    )
    with open(os.path.join(BENCH, "family_er_te_gate.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)

    # ============================== console ==============================
    def fmt(r):
        return (f"{r['gate']:<14} t={str(r['t'] or ''):>4} m={str(r['m'] if r['m'] is not None else ''):>4} "
                f"n={r['n_edges']:>5} TP={r['tp_retained']:>3} TPpm={r['tp_noep_retained']:>2}/{n_tp_noep} "
                f"ncTP={r['noncoding_tp_retained']:>2}/{n_nc_tp} tbar={r['truthbar_retained']:>4} "
                f"gen={r['genuine_overmerge']:>4} blk[prfam/eh]={r['genuine_rate_block']}/"
                f"{r['genuine_rate_block_edgehost']}")
    print("\n================ E_r TE-MASK GATE SWEEP ================")
    for r in rows:
        print("  " + fmt(r))
    print("\n[mask discrimination] class median mask: "
          f"TP={mask_disc['by_class']['TP']['median']} "
          f"TP(protein-missed)={mask_disc['by_class']['TP_protein_missed']['median']} "
          f"truthbar={mask_disc['by_class']['truthbar']['median']} genuine={mask_disc['by_class']['genuine']['median']} "
          f"-> HYPOTHESIS POPULATION-FALSE for protein paralogs (truthbar~genuine); holds only on protein-missed TP")
    print("[mask AUC(genuine>TP) by core band] " + ", ".join(
        f"{b}={mask_disc['auc_by_core_band'][b]['auc_genuine_gt_TP']}" for b in mask_disc['auc_by_core_band']))
    if sep_at_rec:
        print(f"[separation AT recommended] fixed t={sep_at_rec['t']}: mask cuts "
              f"{sep_at_rec['genuine_cut_by_mask']} genuine vs {sep_at_rec['tppm_cut_by_mask']} protein-missed TP "
              f"(ratio {sep_at_rec['genuine_to_tppm_ratio_at_recommended']}x) | attribution: core branch cuts "
              f"{sep_at_rec['attribution']['cut_by_core_branch']}, mask adds only "
              f"{sep_at_rec['attribution']['cut_by_mask_branch']}")
    print(f"[separation across sweep (SUM, double-counts)] {tot_gen_cut} genuine vs {tot_tppm_cut} "
          f"protein-missed TP  (ratio {sep_ratio if sep_ratio is None else round(sep_ratio,2)}x -- NOT the operating point)")
    print(f"\n[core-only knee]  combined t=0.70: TP={core_knee['tp_retained']} "
          f"TPpm={core_knee['tp_noep_retained']}/{n_tp_noep} ncTP={core_knee['noncoding_tp_retained']}/{n_nc_tp} "
          f"gen={core_knee['genuine_overmerge']} blk[prfam/eh]={core_knee['genuine_rate_block']}/"
          f"{core_knee['genuine_rate_block_edgehost']}")
    print(f"[RECOMMENDED]     {recommended['gate']} t={recommended['t']} m={recommended['m']} "
          f"({rec_kind}): TP={recommended['tp_retained']} TPpm={recommended['tp_noep_retained']}/{n_tp_noep} "
          f"ncTP={recommended['noncoding_tp_retained']}/{n_nc_tp} gen={recommended['genuine_overmerge']} "
          f"blk[prfam/eh]={recommended['genuine_rate_block']}/{recommended['genuine_rate_block_edgehost']}")
    if strict_pareto:
        print(f"[strict Pareto alt] {strict_pareto['gate']} t={strict_pareto['t']} m={strict_pareto['m']}: "
              f"TP={strict_pareto['tp_retained']} TPpm={strict_pareto['tp_noep_retained']}/{n_tp_noep} "
              f"ncTP={strict_pareto['noncoding_tp_retained']}/{n_nc_tp} gen={strict_pareto['genuine_overmerge']} "
              f"blk[prfam/eh]={strict_pareto['genuine_rate_block']}/{strict_pareto['genuine_rate_block_edgehost']}")
    print("\n[named pairs -- labeled by WHICH BRANCH decides them]")
    for grp, key in [("core-branch cut (NOT mask)", "te_junk_cut_by_CORE_branch_not_mask"),
                     ("protein-branch keep (NOT mask)", "coding_paralog_kept_by_PROTEIN_not_mask"),
                     ("MASK-branch keep", "noncoding_tp_kept_by_mask")]:
        for rec in named[key]:
            st = rec.get("kept_recommended")
            print(f"   {grp:30s} {rec['pair']:22s} cls={rec.get('cls')} in_ep={rec.get('in_ep')} "
                  f"core={rec.get('core')} mask={rec.get('mask')} kept@rec={st} | {rec.get('decided_by','')}")
    md = named["mask_decided_confirmation"]
    print(f"\n[MASK-DECIDED confirmation @ recommended] band {md['band']}: "
          f"genuine cut-by-mask {md['genuine_cut_by_mask_n']}/{md['n_genuine']} ({md['genuine_cut_rate']}) vs "
          f"real TP cut-by-mask {md['TP_cut_by_mask_n']}/{md['n_TP']} ({md['TP_cut_rate']}) "
          f"-> {md['per_capita_enrichment']}x per-capita")
    print("   TE-junk the mask CUTS (high core, high mask):")
    for e in md["TE_junk_the_mask_CUTS"][:4]:
        print(f"     {e['pair']:36s} core={e['core']} mask={e['mask']}")
    print("   real paralogs the mask KEEPS (high core, low mask):")
    for e in md["real_paralogs_the_mask_KEEPS"][:4]:
        print(f"     {e['pair']:36s} core={e['core']} mask={e['mask']} nc={e['noncoding']}")
    print("   junk the mask CANNOT catch (high core, low mask genuine -- residual):")
    for e in md["junk_the_mask_CANNOT_catch"][:3]:
        print(f"     {e['pair']:36s} core={e['core']} mask={e['mask']}")
    print(f"\n[aligned-region cross-check] available={aligned.get('available')}: "
          f"MAX-inflation on {max_inflation['n_asymmetric_protein_missed_TP']} asymmetric protein-missed TP")
    for rec in aligned.get("pairs", []):
        if rec.get("whole_locus_max") is not None:
            print(f"     {rec['pair']:28s} cls={rec.get('cls')} whole_MAX={rec['whole_locus_max']} "
                  f"maskA={rec.get('mask_a')} maskB={rec.get('mask_b')} aligned={rec['aligned']}")
    print(f"\n[confirmations] {conf}")
    print(f"[noncoding recovered @ recommended] {len(nc_kept)}/{n_nc_tp}: {nc_kept}")
    print(f"[determinism] recommended reproducible within-run={det} (PYTHONHASHSEED=0)")
    print("wrote bench/family_er_te_gate.tsv + bench/family_er_te_gate.json + bench/edge_bridge_mask.tsv")


if __name__ == "__main__":
    main()
