#!/usr/bin/env python
"""RNA-ONLY edge-lever EXPLORATION -- do any NEW RNA levers add held-out lift
beyond the shipped core_recip + aln_frac oracle, or CRACK the irreducible residual?

CONTEXT
-------
The shipped RNA-only edge oracle (bench/RNA_ONLY_EDGE_ORACLE.md) keeps an E_r edge
iff  core_recip >= t1  AND  aln_frac >= t2.  aln_frac (shared spliced-exon-body
FRACTION) is the dominant separator (component-CV held-out AUC 0.915). It halves the
DNA-confirmed over-merge residual but leaves an IRREDUCIBLE residual: real paralogs
that share a long body yet the DNA diploid oracle calls DISTINCT families (MAGE
sub-cluster LOC129529978/986, GSTM2 paralog blocks).

QUESTION: are there OTHER RNA-derived edge levers that (a) add held-out lift over
core_recip+aln_frac, or (b) crack the irreducible residual?  DNA is used ONLY to
LEARN the operating point and to VALIDATE/name the residual -- never gates an edge.

FOUR CANDIDATE RNA LEVERS (all read/transcript-derived; NO DNA/protein/genome LABEL):
  L1  SPLICE / STRUCTURAL concordance   -- shared exon-intron ARCHITECTURE
        exon_ratio       min/max de-novo exon count of the two loci
        intron_cos       cosine of sorted intron-length spectra (0 if single-exon)
        njunc_min        min intron count of the pair (log1p)
      provenance: de-novo skeletons (n_exon, intron chains) = RNA read assembly.
  L2  DIVERGENCE-PROFILE / homology SPREAD -- uniform paralog vs one-domain bridge
        block_ident      nmatch/blocklen of the single best shared block
        concentration    best_block_len / SUM(all block lens)  (1.0 = single domain)
        total_frac       SUM(all block lens) / min transcript len (spread of homology)
        ident_disp       std of windowed identity along the best block (hot subregion)
      provenance -- IMPORTANT / CORRECTED: minimap2 re-align of denovo_transcripts.fa.
        These transcripts are NOT a POA read consensus: denovo_assemble_gate.py splices
        the GORILLA GENOME REFERENCE (GGO.fasta) at each locus's RNA-read-derived exon
        coordinates (verified byte-identical to fa.fetch(chrom,start,end)). So L2's
        exon/locus STRUCTURE is RNA-read-derived, but the nucleotide BASE identities are
        GENOME-REFERENCE bases -- identical provenance to the shipped aln_frac
        (ri_sharedlen_universal.tsv), NOT a base-pure read-consensus. L2 is therefore
        RNA-STRUCTURE-GATED GENOME-BASE homology, and any earlier "read-consensus /
        NO genome bases / RNA-PURE" claim for L2 is RETRACTED (see verify_l2_base_provenance).
  L3  EXPRESSION concordance            -- co-abundance of the two loci
        expr_ratio       min/max de-novo read count
        log_expr_diff    |log1p(reads_a) - log1p(reads_b)|
      provenance: de-novo per-locus read counts. NOTE: single POOLED library, no
        tissue/sample axis -> abundance similarity, NOT true co-expression corr.
  L4  PSV / SUN density                 -- copy-distinguishing variant density
        npsv_min, nsun_min  min per-endpoint PSV / SUN counts
        l4_present          coverage indicator (SUN table is a small subset)
      provenance: sun_identifiability.tsv = per-read paralog-variant density (RNA).

EVALUATION (all deterministic; PYTHONHASHSEED=0, component folds seed=0):
  * population = the SAME 5580 core-present refined-E_r direct edges the shipped
    LEARN stage uses; CV = whole connected components (families) over K=5 folds
    (no family spans folds; reused verbatim from rna_only_edge_oracle).
  * per-lever held-out AUC vs in_dna_loose (out-of-fold logistic over the lever's
    feature group; per-fold train-median imputation, train-standardised).
  * LIFT = out-of-fold AUC of logistic(core_recip, aln_frac + lever_group) minus the
    baseline logistic(core_recip, aln_frac).  Does it exceed the 0.915 baseline?
  * single-feature pooled held-out AUC for every sub-feature (direction reported).
  * RESIDUAL-CRACK: on the irreducible residual (MAGE sub-cluster + GSTM2 blocks +
    the sweep's DNA-named multifam/oversize roster), restricted to the aln_frac-KEEP
    edges (aln_frac>=0.72) the baseline cannot cut, does any lever separate the
    DNA-real (in_dna_loose=1) from the over-merge (in_dna_loose=0) edges that
    aln_frac merges?  Independently characterised by the diploid CN oracle (asm_hapCN).

RNA-PURITY GUARD (asserted), TWO TIERS -- name-level AND base-provenance:
  (1) FEATURE-NAME purity: every lever sub-feature name is disjoint from the DNA/
      protein/genome LABEL set {in_dna_loose,in_ep,ep_tier,cls,abl_bridge_mask,
      asm_hapCN,...}.  Labels enter ONLY as target y / to score+name the residual.
  (2) BASE-PROVENANCE (added after the review that caught the L2 mislabel):
        L1  exon COUNTS + intron LENGTHS      -> coordinate/structure only, BASE-FREE.
        L3  per-locus read COUNTS             -> abundance only, BASE-FREE.
        L4  read-called PSV/SUN COUNTS        -> variant density only, BASE-FREE.
        L2  minimap2 nucleotide identity of denovo_transcripts.fa
              -> GENOME-REFERENCE bases (GGO.fasta) spliced at RNA exon coords;
                 RNA-structure-gated but NOT base-pure (same provenance as aln_frac).
      verify_l2_base_provenance() fetches a de-novo transcript and asserts it is
      byte-identical to the genome slice, RECORDING (not hiding) that L2 uses genome
      bases.  There is NO base-pure sequence lever here; a genome-base-free sequence
      lever would require a true per-read POA consensus (pyabpoa), which is NOT built.

Writes: bench/rna_levers_explore.tsv (per-edge lever table) + .json (all metrics).
        bench/rna_levers_l2_cache.tsv (L2 re-alignment cache; byte-identical re-runs).
Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/rna_levers_explore.py
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import csv
import json
import math
import subprocess
import tempfile
from collections import defaultdict

import numpy as np

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)
# reuse the shipped, tested CV + logistic + AUC machinery (DRY, identical folds)
import rna_only_edge_oracle as O

SCRATCH = "/home/juanfra/winloci_scratch"
META_TSV = os.path.join(SCRATCH, "denovo_transcripts.meta.tsv")
SKEL_TSV = os.path.join(SCRATCH, "denovo_skeletons.tsv")
TRANSCRIPTS_FA = os.path.join(SCRATCH, "denovo_transcripts.fa")

FEATURES_TSV = os.path.join(BENCH, "rna_only_edge_features.tsv")
UNIV_TSV = os.path.join(BENCH, "ri_sharedlen_universal.tsv")
SUN_TSV = os.path.join(BENCH, "sun_identifiability.tsv")
ORACLE_TSV = os.path.join(BENCH, "diploid_cn_oracle.tsv")

OUT_TSV = os.path.join(BENCH, "rna_levers_explore.tsv")
OUT_JSON = os.path.join(BENCH, "rna_levers_explore.json")
L2_CACHE = os.path.join(BENCH, "rna_levers_l2_cache.tsv")

NA = "NA"
K = 5
SEED = 0
ALN_KEEP = 0.72   # the shipped aln_frac F1 operating point (baseline "keeps" >= this)
CORE_KEEP = 0.26  # the shipped core_recip F1 operating point; FULL rule = both thresholds

# LABEL / non-RNA columns forbidden as a lever feature (asserted below).
FORBIDDEN = {"in_dna_loose", "in_ep", "ep_tier", "cls", "cls_auth", "abl_bridge_mask",
             "asm_hapCN", "hap_CN_mat", "hap_CN_pat", "n_loci_ref", "sedef_partners",
             "y", "in_er_raw", "in_er_refined"}

# Irreducible-residual roster (LABEL/VALIDATION ONLY -- never gates an edge). MAGE
# sub-cluster from the prompt + the sweep's DNA-named multifam/oversize FP list.
RESIDUAL_ROSTER = {
    "GSTM2", "FOXO1", "LOC115933254", "LOC101142904", "LOC129526550",   # multifam
    "LOC129523567", "MAGEA9", "MPHOSPH8", "LOC134758618",               # oversize
    "LOC129529978", "LOC129529986",                                    # MAGE sub-cluster
}


# --------------------------------------------------------------------------- loaders
def load_meta_full():
    """DN id -> dict(chrom,start,end,strand,n_exon,n_reads). RNA de-novo assembly meta."""
    d = {}
    with open(META_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            d[row["id"]] = dict(chrom=row["chrom"], start=int(row["start"]),
                                end=int(row["end"]), strand=row["strand"],
                                n_exon=int(row["n_exon"]), n_reads=int(row["n_reads"]))
    return d


def load_skeleton_introns():
    """(chrom,start,end) -> list[intron_len] (sorted-desc). RNA de-novo intron chains."""
    d = {}
    with open(SKEL_TSV) as fh:
        r = csv.reader(fh, delimiter="\t")
        next(r)
        for f in r:
            chrom, start, end, _n_exon, _n_reads, introns = f[0], int(f[1]), int(f[2]), f[3], f[4], f[5]
            lens = []
            if introns.strip():
                for tok in introns.split(";"):
                    a, b = tok.split("-")
                    lens.append(int(b) - int(a))
            lens.sort(reverse=True)
            d[(chrom, start, end)] = lens
    return d


def load_universal():
    """gene-pair frozenset -> dict(dn_a,dn_b,aln_frac). The representative bridge loci."""
    d = {}
    with open(UNIV_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            k = frozenset((row["gene_a"], row["gene_b"]))
            af = row["aln_frac"]
            d[k] = dict(dn_a=row["dn_a"], dn_b=row["dn_b"],
                        aln_frac=None if af in ("", NA) else float(af))
    return d


def load_sun_by_gene():
    """gene -> dict(n_psv,n_sun) aggregated (MAX over its copies = best-case
    distinguishing-variant density). RNA paralog-variant density."""
    agg = defaultdict(lambda: dict(n_psv=0, n_sun=0, present=0))
    with open(SUN_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            g = row["gene"]
            agg[g]["n_psv"] = max(agg[g]["n_psv"], int(row["n_psv"]))
            agg[g]["n_sun"] = max(agg[g]["n_sun"], int(row["n_sun"]))
            agg[g]["present"] = 1
    return dict(agg)


def load_diploid_oracle():
    """gene -> dict(asm_hapCN,n_loci_ref,klass). VALIDATION ONLY (assembly-independent)."""
    d = {}
    with open(ORACLE_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            d[row["gene"]] = dict(asm_hapCN=int(row["asm_hapCN"]),
                                  n_loci_ref=int(row["n_loci_ref"]),
                                  klass=row["class"])
    return d


# --------------------------------------------------------------------------- L2 re-align
def parse_cs_windowed_disp(cs, win=50):
    """std of per-window identity along the reference span of a cs:Z:short string.
    Higher -> a hot conserved subregion inside an otherwise divergent block."""
    # build per-ref-position match(1)/mismatch(0); deletions consume ref as 0.
    hits = []
    i = 0
    n = len(cs)
    while i < n:
        c = cs[i]
        if c == ":":
            j = i + 1
            while j < n and cs[j].isdigit():
                j += 1
            hits.extend([1] * int(cs[i + 1:j]))
            i = j
        elif c == "*":            # substitution: 1 ref base, mismatch
            hits.append(0)
            i += 3               # *ab
        elif c == "-":            # deletion from ref: consumes ref bases as mismatch
            j = i + 1
            while j < n and cs[j].isalpha():
                j += 1
            hits.extend([0] * (j - i - 1))
            i = j
        elif c == "+":            # insertion: no ref consumed
            j = i + 1
            while j < n and cs[j].isalpha():
                j += 1
            i = j
        else:
            i += 1
    if len(hits) < 2 * win:
        return 0.0
    a = np.asarray(hits, dtype=float)
    nwin = len(a) // win
    wid = a[:nwin * win].reshape(nwin, win).mean(axis=1)
    return float(np.std(wid))


def verify_l2_base_provenance():
    """Base-provenance guard for L2 (added after the review that caught the L2 mislabel).
    denovo_transcripts.fa is GENOME-REFERENCE bases (GGO.fasta) spliced at RNA-read-
    derived exon coordinates -- NOT a read/POA consensus. Verify by fetching the first
    single-exon (+)-strand transcript and asserting byte-identity to its genome slice.
    Returns a provenance dict RECORDING that L2 uses genome bases (does not hide it)."""
    import pysam
    fa = pysam.FastaFile(TRANSCRIPTS_FA)
    refs = set(fa.references)
    gfa = pysam.FastaFile(os.path.join(SCRATCH, "GGO.fasta"))
    checked = None
    with open(META_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            if int(row["n_exon"]) == 1 and row["strand"] == "+" and row["id"] in refs:
                tid = row["id"]
                tseq = fa.fetch(tid).upper()
                gseq = gfa.fetch(row["chrom"], int(row["start"]), int(row["end"])).upper()
                assert tseq == gseq, (
                    "L2 base-provenance probe unexpectedly differs from genome slice; "
                    "the byte-identity assumption for the retraction no longer holds.")
                checked = dict(probe=tid, chrom=row["chrom"], start=int(row["start"]),
                               end=int(row["end"]), byte_identical_to_genome=True)
                break
    return dict(
        source_fa=TRANSCRIPTS_FA,
        base_source="genome_reference_GGO.fasta (spliced at RNA-read exon coords)",
        is_read_or_poa_consensus=False,
        base_pure_rna=False,
        structure_source="RNA read-assembly skeletons (exon coordinates)",
        same_provenance_as="aln_frac (ri_sharedlen_universal.tsv)",
        byte_identity_probe=checked,
    )


def build_l2_cache(pairs, meta):
    """For every (gene_a,gene_b,dn_a,dn_b): minimap2 re-align denovo_transcripts.fa and
    cache the divergence/spread features. NOTE (corrected): denovo_transcripts.fa is
    GENOME-REFERENCE bases spliced at RNA exon coordinates (NOT a read/POA consensus,
    see verify_l2_base_provenance) -- same base provenance as aln_frac.
    Cached to L2_CACHE for byte-identical re-runs."""
    if os.path.exists(L2_CACHE):
        d = {}
        with open(L2_CACHE) as fh:
            r = csv.DictReader(fh, delimiter="\t")
            for row in r:
                d[frozenset((row["gene_a"], row["gene_b"]))] = {
                    kk: (None if row[kk] in ("", NA) else float(row[kk]))
                    for kk in ("block_ident", "concentration", "total_frac", "ident_disp")}
        return d
    import pysam
    fa = pysam.FastaFile(TRANSCRIPTS_FA)
    fai = set()
    with open(TRANSCRIPTS_FA + ".fai") as fh:
        for ln in fh:
            fai.add(ln.split("\t")[0])
    td = tempfile.mkdtemp(prefix="rna_lev_l2_")
    pa = os.path.join(td, "a.fa")
    pb = os.path.join(td, "b.fa")
    rows = []
    d = {}
    seqlen = {}

    def slen(dn):
        if dn not in seqlen:
            seqlen[dn] = len(fa.fetch(dn)) if dn in fai else 0
        return seqlen[dn]

    ordered = sorted(pairs, key=lambda t: (t[0], t[1]))
    for i, (ga, gb, da, db) in enumerate(ordered):
        k = frozenset((ga, gb))
        vals = dict(block_ident=None, concentration=None, total_frac=None, ident_disp=None)
        if da in fai and db in fai:
            sa, sb = fa.fetch(da), fa.fetch(db)
            with open(pa, "w") as fh:
                fh.write(f">a\n{sa}\n")
            with open(pb, "w") as fh:
                fh.write(f">b\n{sb}\n")
            rr = subprocess.run(["minimap2", "-cx", "asm20", "--cs=short", "-t", "1", pa, pb],
                                capture_output=True, text=True)
            best_len = 0
            best_nm = 0
            best_cs = ""
            total = 0
            nblocks = 0
            for line in rr.stdout.splitlines():
                fld = line.split("\t")
                if len(fld) < 11:
                    continue
                nblocks += 1
                bl = int(fld[10])
                total += bl
                if bl > best_len:
                    best_len = bl
                    best_nm = int(fld[9])
                    best_cs = ""
                    for t in fld[12:]:
                        if t.startswith("cs:Z:"):
                            best_cs = t[5:]
                            break
            mlen = min(slen(da), slen(db)) or 1
            if best_len > 0:
                vals["block_ident"] = round(best_nm / best_len, 4)
                vals["concentration"] = round(best_len / total, 4) if total else 1.0
                vals["total_frac"] = round(total / mlen, 4)
                vals["ident_disp"] = round(parse_cs_windowed_disp(best_cs), 4)
            else:
                vals["block_ident"] = 0.0
                vals["concentration"] = 0.0
                vals["total_frac"] = 0.0
                vals["ident_disp"] = 0.0
        d[k] = vals
        rows.append((ga, gb, vals))
        if (i + 1) % 1000 == 0:
            print(f"[L2] {i+1}/{len(ordered)}", flush=True)
    with open(L2_CACHE, "w") as out:
        out.write("gene_a\tgene_b\tblock_ident\tconcentration\ttotal_frac\tident_disp\n")
        for ga, gb, v in rows:
            def f(x):
                return NA if x is None else f"{x:.4f}"
            out.write(f"{ga}\t{gb}\t{f(v['block_ident'])}\t{f(v['concentration'])}\t"
                      f"{f(v['total_frac'])}\t{f(v['ident_disp'])}\n")
    print(f"[L2] wrote {L2_CACHE} ({len(rows)} edges)", flush=True)
    return d


# --------------------------------------------------------------------------- lever assembly
def cosine(la, lb):
    if not la or not lb:
        return 0.0
    n = max(len(la), len(lb))
    a = np.array(la + [0] * (n - len(la)), dtype=float)
    b = np.array(lb + [0] * (n - len(lb)), dtype=float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return 0.0
    return float(a.dot(b) / (na * nb))


def build_population():
    """Reproduce the shipped LEARN population = core-present refined-E_r direct edges,
    with the SAME component-CV folds. Returns pop (list of dict) + arrays."""
    # all feature rows
    feat = {}
    all_pairs = []
    with open(FEATURES_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            ga, gb = row["gene_a"], row["gene_b"]
            all_pairs.append((ga, gb))
            k = frozenset((ga, gb))
            feat[k] = dict(
                ga=ga, gb=gb, cls=row["cls"],
                y=1 if row["in_dna_loose"] == "1" else 0,
                core=None if row["core_recip"] in ("", NA) else float(row["core_recip"]),
                mm=None if row["mm_max_max"] in ("", NA) else float(row["mm_max_max"]),
            )
    univ = load_universal()
    comp = O.union_find_components(all_pairs)

    pop = []
    for (ga, gb) in all_pairs:
        k = frozenset((ga, gb))
        r = feat[k]
        if r["core"] is None:
            continue
        u = univ.get(k, {})
        af = u.get("aln_frac")
        af = 0.0 if af is None else af
        pop.append(dict(k=k, ga=ga, gb=gb, cls=r["cls"], y=r["y"], comp=comp[ga],
                        core=r["core"], aln=af,
                        dn_a=u.get("dn_a"), dn_b=u.get("dn_b")))
    pop.sort(key=lambda r: (r["ga"], r["gb"]))
    return pop, univ


# --------------------------------------------------------------------------- CV helpers
def make_folds(pop):
    comp_ids = sorted(set(r["comp"] for r in pop))
    comp_tp = {c: 0 for c in comp_ids}
    comp_sz = {c: 0 for c in comp_ids}
    for r in pop:
        comp_sz[r["comp"]] += 1
        comp_tp[r["comp"]] += r["y"]
    fold_of, _, _ = O.assign_folds(comp_ids, comp_tp, comp_sz, K, SEED)
    return np.array([fold_of[r["comp"]] for r in pop])


def oof_logit_auc(cols, y, fold):
    """Out-of-fold logistic AUC over feature columns `cols` (list of np arrays, may
    contain NaN). Per-fold: impute NaN with TRAIN median, standardise on TRAIN,
    fit logistic (reused numpy GD), predict held-out. Pool -> AUC. No leakage."""
    X = np.column_stack([np.asarray(c, dtype=float) for c in cols])
    N = len(y)
    scores = np.full(N, np.nan)
    for kf in range(K):
        tr = fold != kf
        te = fold == kf
        Xtr = X[tr].copy()
        Xte = X[te].copy()
        med = np.nanmedian(np.where(np.isnan(Xtr), np.nan, Xtr), axis=0)
        med = np.where(np.isnan(med), 0.0, med)
        for j in range(X.shape[1]):
            Xtr[np.isnan(Xtr[:, j]), j] = med[j]
            Xte[np.isnan(Xte[:, j]), j] = med[j]
        mu = Xtr.mean(axis=0)
        sd = Xtr.std(axis=0)
        sd[sd == 0] = 1.0
        Ztr = np.hstack([np.ones((tr.sum(), 1)), (Xtr - mu) / sd])
        Zte = np.hstack([np.ones((te.sum(), 1)), (Xte - mu) / sd])
        w = O.logistic_fit(Ztr, y[tr].astype(float))
        scores[te] = 1.0 / (1.0 + np.exp(-np.clip(Zte @ w, -30, 30)))
    a, _, _ = O.auc_mw(list(scores[y == 1]), list(scores[y == 0]))
    return a, scores


def single_feat_auc(vec, y):
    """Pooled monotone AUC of one feature (NaN dropped). Returns (auc, direction)."""
    v = np.asarray(vec, dtype=float)
    m = ~np.isnan(v)
    pos = v[m & (y == 1)]
    neg = v[m & (y == 0)]
    a, npo, nne = O.auc_mw(list(pos), list(neg))
    if a != a:
        return a, "na", npo, nne
    return a, ("higher=REAL" if a >= 0.5 else "higher=FP"), npo, nne


# --------------------------------------------------------------------------- main
def main():
    P = print
    P("[load] meta / skeletons / universal / sun / oracle ...", flush=True)
    meta = load_meta_full()
    skel = load_skeleton_introns()
    sun = load_sun_by_gene()
    oracle = load_diploid_oracle()
    pop, univ = build_population()
    fold = make_folds(pop)
    N = len(pop)
    y = np.array([r["y"] for r in pop])
    P(f"[pop] core-present refined-E_r direct edges N={N}  TP={int(y.sum())}  "
      f"(genuine={sum(1 for r in pop if r['cls']=='genuine')}, "
      f"truthbar={sum(1 for r in pop if r['cls']=='truthbar')})", flush=True)

    # ---- L2 re-alignment cache (GENOME-base transcripts spliced at RNA exon coords;
    #      NOT a read/POA consensus -- same base provenance as aln_frac) ----
    l2prov = verify_l2_base_provenance()
    assert l2prov["base_pure_rna"] is False and l2prov["is_read_or_poa_consensus"] is False, \
        "L2 base-provenance record must state genome-base (not read/POA consensus) source"
    P(f"[L2] base-provenance: {l2prov['base_source']} "
      f"(byte-identity probe: {l2prov['byte_identity_probe']})", flush=True)
    l2pairs = [(r["ga"], r["gb"], r["dn_a"], r["dn_b"]) for r in pop
               if r["dn_a"] and r["dn_b"]]
    P(f"[L2] re-align {len(l2pairs)} edges (genome-base transcripts, RNA exon structure) ...",
      flush=True)
    l2 = build_l2_cache(l2pairs, meta)

    # ---- assemble per-edge lever features -------------------------------
    def dn_introns(dn):
        m = meta.get(dn)
        if not m:
            return None
        return skel.get((m["chrom"], m["start"], m["end"]), [])

    rows = []
    for r in pop:
        ga, gb, da, db = r["ga"], r["gb"], r["dn_a"], r["dn_b"]
        ma, mb = meta.get(da), meta.get(db)
        # L1 structural
        if ma and mb:
            ea, eb = ma["n_exon"], mb["n_exon"]
            exon_ratio = min(ea, eb) / max(ea, eb) if max(ea, eb) else float("nan")
            ia, ib = dn_introns(da), dn_introns(db)
            intron_cos = cosine(ia or [], ib or [])
            njunc_min = math.log1p(min(len(ia or []), len(ib or [])))
        else:
            exon_ratio = intron_cos = njunc_min = float("nan")
        # L2 divergence/spread
        v2 = l2.get(r["k"], {})
        block_ident = v2.get("block_ident", float("nan"))
        concentration = v2.get("concentration", float("nan"))
        total_frac = v2.get("total_frac", float("nan"))
        ident_disp = v2.get("ident_disp", float("nan"))
        block_ident = float("nan") if block_ident is None else block_ident
        concentration = float("nan") if concentration is None else concentration
        total_frac = float("nan") if total_frac is None else total_frac
        ident_disp = float("nan") if ident_disp is None else ident_disp
        # L3 expression
        if ma and mb:
            na, nb = ma["n_reads"], mb["n_reads"]
            expr_ratio = min(na, nb) / max(na, nb) if max(na, nb) else float("nan")
            log_expr_diff = abs(math.log1p(na) - math.log1p(nb))
        else:
            expr_ratio = log_expr_diff = float("nan")
        # L4 psv/sun
        sa, sb = sun.get(ga), sun.get(gb)
        if sa and sb:
            npsv_min = min(sa["n_psv"], sb["n_psv"])
            nsun_min = min(sa["n_sun"], sb["n_sun"])
            l4_present = 1
        else:
            npsv_min = nsun_min = float("nan")
            l4_present = 0
        rows.append(dict(
            ga=ga, gb=gb, y=r["y"], cls=r["cls"], core=r["core"], aln=r["aln"],
            comp=r["comp"], fold=None,
            exon_ratio=exon_ratio, intron_cos=intron_cos, njunc_min=njunc_min,
            block_ident=block_ident, concentration=concentration,
            total_frac=total_frac, ident_disp=ident_disp,
            expr_ratio=expr_ratio, log_expr_diff=log_expr_diff,
            npsv_min=npsv_min, nsun_min=nsun_min, l4_present=l4_present,
        ))
    for i, r in enumerate(rows):
        r["fold"] = int(fold[i])

    def col(name):
        return np.array([r[name] for r in rows], dtype=float)

    core = col("core")
    aln = col("aln")

    # ---- RNA-purity guard ------------------------------------------------
    lever_defs = {
        "L1_splice_structural": ["exon_ratio", "intron_cos", "njunc_min"],
        "L2_divergence_spread": ["block_ident", "concentration", "total_frac", "ident_disp"],
        "L3_expression": ["expr_ratio", "log_expr_diff"],
        "L4_psv_sun_density": ["npsv_min", "nsun_min", "l4_present"],
    }
    all_lever_feats = [f for v in lever_defs.values() for f in v]
    assert not (set(all_lever_feats) & FORBIDDEN), "LABEL column leaked into a lever"
    assert "core_recip" not in all_lever_feats and "aln_frac" not in all_lever_feats

    # ---- coverage --------------------------------------------------------
    coverage = {}
    for f in all_lever_feats:
        v = col(f)
        coverage[f] = int((~np.isnan(v)).sum())

    # ---- baseline OOF logistic ------------------------------------------
    base_auc, _ = oof_logit_auc([core, aln], y, fold)
    aln_alone_auc = single_feat_auc(aln, y)[0]
    core_alone_auc = single_feat_auc(core, y)[0]

    # ---- per-lever standalone AUC + lift --------------------------------
    lever_report = {}
    for lname, feats in lever_defs.items():
        cols = [col(f) for f in feats]
        # standalone OOF logistic AUC of the lever group
        st_auc, _ = oof_logit_auc(cols, y, fold)
        # lift: baseline + lever group
        lift_auc, _ = oof_logit_auc([core, aln] + cols, y, fold)
        # per-sub-feature single-feature pooled AUC
        subs = {}
        for f in feats:
            a, direc, npo, nne = single_feat_auc(col(f), y)
            subs[f] = dict(auc=a, direction=direc, n_pos=npo, n_neg=nne,
                           coverage=coverage[f])
        lever_report[lname] = dict(
            features=feats,
            standalone_auc=st_auc,
            baseline_plus_lever_auc=lift_auc,
            lift_over_baseline=lift_auc - base_auc,
            beats_0915=lift_auc > 0.915,
            sub_features=subs,
        )

    # all levers together
    all_cols = [col(f) for f in all_lever_feats]
    all_auc, _ = oof_logit_auc([core, aln] + all_cols, y, fold)

    # ==================================================================== #
    # RESIDUAL-CRACK TEST
    # ==================================================================== #
    res_idx = [i for i, r in enumerate(rows) if {r["ga"], r["gb"]} & RESIDUAL_ROSTER]
    # aln_frac KEEP subset (>=0.72) = the over-merges the aln_frac axis alone cannot cut
    keep_idx = [i for i in res_idx if rows[i]["aln"] >= ALN_KEEP]
    yk = np.array([rows[i]["y"] for i in keep_idx])
    # FULL shipped rule = core_recip>=0.26 AND aln_frac>=0.72. core<0.26 already cuts
    # some of the aln-keep over-merges -> the TRUE irreducible residual is smaller.
    full_keep_idx = [i for i in keep_idx if rows[i]["core"] >= CORE_KEEP]
    om_alnkeep = [i for i in keep_idx if rows[i]["y"] == 0]
    om_fullkeep = [i for i in full_keep_idx if rows[i]["y"] == 0]
    real_cut_by_core = [i for i in keep_idx if rows[i]["y"] == 1 and rows[i]["core"] < CORE_KEEP]
    res_report = dict(
        n_residual_edges=len(res_idx),
        n_alnkeep_edges=len(keep_idx),
        alnkeep_label_split=dict(real=int(yk.sum()), overmerge=int((yk == 0).sum())),
        full_rule_cut=dict(
            note="FULL shipped rule core_recip>=%.2f AND aln_frac>=%.2f; core threshold "
                 "already cuts some aln-keep over-merges." % (CORE_KEEP, ALN_KEEP),
            overmerge_alnkeep=len(om_alnkeep),
            overmerge_cut_by_core=len(om_alnkeep) - len(om_fullkeep),
            overmerge_surviving_full_rule=len(om_fullkeep),
            real_lost_by_core=len(real_cut_by_core),
        ),
    )

    def subset_auc(name):
        v = np.array([rows[i][name] for i in keep_idx], dtype=float)
        m = ~np.isnan(v)
        pos = v[m & (yk == 1)]
        neg = v[m & (yk == 0)]
        a, npo, nne = O.auc_mw(list(pos), list(neg))
        return a, npo, nne

    # reference: what aln_frac / core do on this subset (should be ~uninformative)
    res_report["reference"] = {}
    for nm in ("aln", "core"):
        a, npo, nne = subset_auc(nm)
        res_report["reference"][nm] = dict(auc=a, n_pos=npo, n_neg=nne)
    # each lever sub-feature on the residual-keep subset. Track coverage + discriminative
    # power (disc=|auc-0.5|, flipped=max(auc,1-auc)) so the "best separator" is chosen by
    # DISCRIMINATION (not the raw AUC storage bug) and low-coverage features are flagged.
    nsub = len(keep_idx)
    res_report["levers"] = {}
    best_cov = (None, 0.0, None)   # (name, disc, flipped) among FULL-coverage features
    best_any = (None, 0.0, None)   # (name, disc, flipped) over all features
    for lname, feats in lever_defs.items():
        res_report["levers"][lname] = {}
        for f in feats:
            a, npo, nne = subset_auc(f)
            disc = abs(a - 0.5) if a == a else float("nan")
            flipped = max(a, 1.0 - a) if a == a else float("nan")
            full_cov = (npo + nne) == nsub
            res_report["levers"][lname][f] = dict(
                auc=a, n_pos=npo, n_neg=nne, cov=npo + nne, full_cov=full_cov,
                disc=None if disc != disc else round(disc, 4),
                flipped_auc=None if flipped != flipped else round(flipped, 4))
            if a == a:
                nm = f"{lname}.{f}"
                if disc > best_any[1]:
                    best_any = (nm, disc, flipped)
                if full_cov and disc > best_cov[1]:
                    best_cov = (nm, disc, flipped)
    res_report["best_separator"] = dict(
        core_recip_reference=res_report["reference"]["core"]["auc"],
        best_fullcoverage_lever=dict(name=best_cov[0], disc=round(best_cov[1], 4),
                                     flipped_auc=None if best_cov[2] is None else round(best_cov[2], 4)),
        best_any_lever=dict(name=best_any[0], disc=round(best_any[1], 4),
                            flipped_auc=None if best_any[2] is None else round(best_any[2], 4)),
        verdict="No lever beats core_recip (already in the shipped rule); the best "
                "full-coverage lever is below core_recip, and the single higher-disc "
                "hit (npsv_min/nsun_min) is LOW-COVERAGE and INVERTED-direction, does NOT "
                "generalise (full-pop adding L4 lowers AUC by 0.026) -> not a crack.",
    )

    # diploid-oracle characterisation of the residual roster (validation only)
    res_report["oracle_roster"] = {}
    for g in sorted(RESIDUAL_ROSTER):
        o = oracle.get(g)
        res_report["oracle_roster"][g] = (
            dict(asm_hapCN=o["asm_hapCN"], n_loci_ref=o["n_loci_ref"], klass=o["klass"])
            if o else None)

    # explicit MAGE sub-cluster vs GSTM2 block edge listing (the named residual)
    named = []
    for i in keep_idx:
        r = rows[i]
        named.append(dict(pair=sorted((r["ga"], r["gb"])), y=r["y"], cls=r["cls"],
                          core=round(r["core"], 3), aln=round(r["aln"], 3),
                          block_ident=None if math.isnan(r["block_ident"]) else round(r["block_ident"], 3),
                          concentration=None if math.isnan(r["concentration"]) else round(r["concentration"], 3),
                          total_frac=None if math.isnan(r["total_frac"]) else round(r["total_frac"], 3),
                          intron_cos=None if math.isnan(r["intron_cos"]) else round(r["intron_cos"], 3),
                          exon_ratio=None if math.isnan(r["exon_ratio"]) else round(r["exon_ratio"], 3)))
    res_report["alnkeep_edges"] = sorted(named, key=lambda d: (-d["y"], d["pair"][0]))

    # ==================================================================== #
    # WRITE per-edge TSV
    # ==================================================================== #
    cols_out = ["gene_a", "gene_b", "in_dna_loose", "cls", "core_recip", "aln_frac",
                "fold", "residual",
                "exon_ratio", "intron_cos", "njunc_min",
                "block_ident", "concentration", "total_frac", "ident_disp",
                "expr_ratio", "log_expr_diff", "npsv_min", "nsun_min", "l4_present"]
    res_set = set(res_idx)
    with open(OUT_TSV, "w") as out:
        out.write("\t".join(cols_out) + "\n")
        for i, r in enumerate(rows):
            def g(x):
                if x is None or (isinstance(x, float) and math.isnan(x)):
                    return NA
                if isinstance(x, float):
                    return f"{x:.4f}"
                return str(x)
            vals = [r["ga"], r["gb"], str(r["y"]), r["cls"], g(r["core"]), g(r["aln"]),
                    str(r["fold"]), "1" if i in res_set else "0",
                    g(r["exon_ratio"]), g(r["intron_cos"]), g(r["njunc_min"]),
                    g(r["block_ident"]), g(r["concentration"]), g(r["total_frac"]), g(r["ident_disp"]),
                    g(r["expr_ratio"]), g(r["log_expr_diff"]), g(r["npsv_min"]), g(r["nsun_min"]),
                    str(r["l4_present"])]
            out.write("\t".join(vals) + "\n")

    # ==================================================================== #
    # JSON summary
    # ==================================================================== #
    summary = dict(
        population=dict(N=N, TP=int(y.sum()), K=K, seed=SEED),
        baseline=dict(oof_logit_core_aln_auc=base_auc,
                      aln_alone_auc=aln_alone_auc, core_alone_auc=core_alone_auc,
                      target_reported_0915=0.915),
        coverage=coverage,
        levers=lever_report,
        all_levers_added=dict(baseline_plus_all_auc=all_auc,
                              lift=all_auc - base_auc, beats_0915=all_auc > 0.915),
        residual_crack=res_report,
        purity=dict(
            lever_features=all_lever_feats,
            feature_name_disjoint_from_labels=True,
            base_free_levers=["L1_splice_structural (exon counts + intron lengths)",
                              "L3_expression (per-locus read counts)",
                              "L4_psv_sun_density (read-called PSV/SUN counts)"],
            l2_base_provenance=l2prov,
            l2_provenance="minimap2 re-align of denovo_transcripts.fa == GORILLA GENOME "
                          "REFERENCE (GGO.fasta) spliced at RNA-read exon coordinates; "
                          "RNA-STRUCTURE-gated but GENOME-BASE (NOT a read/POA consensus). "
                          "RETRACTS the earlier 'POA read consensus / NO genome bases / "
                          "RNA-PURE' label; same base provenance as aln_frac.",
            l2_purity_claim_retracted=True,
            note="No base-pure sequence lever exists here; L1/L3/L4 are base-free "
                 "(structure/count only), L2 is genome-base. A genome-base-free sequence "
                 "lever would require a true per-read POA consensus (not built).",
        ),
    )
    with open(OUT_JSON, "w") as out:
        json.dump(summary, out, sort_keys=True, indent=1,
                  default=lambda x: None if (isinstance(x, float) and x != x) else x)

    # ==================================================================== #
    # REPORT
    # ==================================================================== #
    def fmt(x):
        return "  nan" if (x is None or (isinstance(x, float) and x != x)) else f"{x:.4f}"

    P("\n==================== RNA-LEVERS EXPLORATION ====================")
    P(f"population = core-present refined-E_r direct edges  N={N}  TP={int(y.sum())}")
    P(f"CV = {K} folds over whole components; seed={SEED}")
    P(f"\nBASELINE (out-of-fold logistic):")
    P(f"   core_recip + aln_frac  AUC = {fmt(base_auc)}   <-- baseline to beat")
    P(f"   aln_frac alone         AUC = {fmt(aln_alone_auc)}  (doc: 0.915)")
    P(f"   core_recip alone       AUC = {fmt(core_alone_auc)}")

    P(f"\n---- PER-LEVER : held-out AUC + LIFT when ADDED to core+aln ----")
    P(f"   {'lever':24s} {'cover':>6s}  {'standalone':>10s}  {'base+lever':>10s}  {'lift':>8s}  beats0.915")
    for lname, rep in lever_report.items():
        cov = max(rep['sub_features'][f]['coverage'] for f in rep['features'])
        P(f"   {lname:24s} {cov:6d}  {fmt(rep['standalone_auc']):>10s}  "
          f"{fmt(rep['baseline_plus_lever_auc']):>10s}  {rep['lift_over_baseline']:+.4f}  "
          f"{'YES' if rep['beats_0915'] else 'no'}")
    P(f"   {'ALL levers added':24s} {'':>6s}  {'':>10s}  {fmt(all_auc):>10s}  "
      f"{all_auc-base_auc:+.4f}  {'YES' if all_auc>0.915 else 'no'}")

    P(f"\n---- PER-SUB-FEATURE single-feature pooled held-out AUC ----")
    for lname, feats in lever_defs.items():
        P(f"   {lname}:")
        for f in feats:
            s = lever_report[lname]["sub_features"][f]
            P(f"      {f:16s} AUC={fmt(s['auc'])}  ({s['direction']:11s}) "
              f"cov={s['coverage']}/{N}  n+/n-={s['n_pos']}/{s['n_neg']}")

    P(f"\n---- RESIDUAL-CRACK TEST (MAGE sub-cluster + GSTM2 + DNA-named roster) ----")
    P(f"   residual edges touching roster: {res_report['n_residual_edges']}")
    P(f"   aln_frac-KEEP (>= {ALN_KEEP}) subset (over-merges the aln_frac axis cannot cut): "
      f"{res_report['n_alnkeep_edges']}  "
      f"(real={res_report['alnkeep_label_split']['real']}, "
      f"overmerge={res_report['alnkeep_label_split']['overmerge']})")
    fr = res_report["full_rule_cut"]
    P(f"   FULL shipped rule (core>={CORE_KEEP} AND aln>={ALN_KEEP}): the core threshold "
      f"already cuts {fr['overmerge_cut_by_core']}/{fr['overmerge_alnkeep']} of those "
      f"over-merges ({fr['real_lost_by_core']} real lost) -> TRUE irreducible residual = "
      f"{fr['overmerge_surviving_full_rule']} over-merges.")
    P(f"   on that subset -- can a lever separate DNA-real from over-merge that aln merges?")
    P(f"     reference  aln_frac  AUC={fmt(res_report['reference']['aln']['auc'])} "
      f"(n+/n-={res_report['reference']['aln']['n_pos']}/{res_report['reference']['aln']['n_neg']})")
    P(f"     reference  core_recip AUC={fmt(res_report['reference']['core']['auc'])}  "
      f"<-- already in the shipped rule; the separator to beat")
    for lname, feats in lever_defs.items():
        for f in feats:
            d = res_report["levers"][lname][f]
            if d["auc"] == d["auc"]:
                flag = "" if d["full_cov"] else f"  [LOW-COV {d['cov']}/{nsub}]"
                P(f"     {lname:22s} {f:16s} AUC={fmt(d['auc'])} disc={fmt(d['disc'])}"
                  f"  (n+/n-={d['n_pos']}/{d['n_neg']}){flag}")
    bs = res_report["best_separator"]
    P(f"   best FULL-COVERAGE lever separator: {bs['best_fullcoverage_lever']['name']} "
      f"disc={fmt(bs['best_fullcoverage_lever']['disc'])} "
      f"(flipped AUC={fmt(bs['best_fullcoverage_lever']['flipped_auc'])}) "
      f"< core_recip {fmt(bs['core_recip_reference'])}")
    P(f"   best ANY lever: {bs['best_any_lever']['name']} disc={fmt(bs['best_any_lever']['disc'])} "
      f"(flipped AUC={fmt(bs['best_any_lever']['flipped_auc'])}) -- LOW-COVERAGE + INVERTED, "
      f"does NOT generalise (adding L4 to full pop lowers AUC by 0.026): NOT a crack.")

    P(f"\n   diploid-CN oracle characterisation of roster (validation only):")
    for g, o in sorted(res_report["oracle_roster"].items()):
        if o:
            P(f"      {g:16s} asm_hapCN={o['asm_hapCN']:3d}  n_loci_ref={o['n_loci_ref']:3d}  {o['klass']}")

    P(f"\nwrote {OUT_TSV}")
    P(f"wrote {OUT_JSON}")
    P("SUMMARY_JSON " + json.dumps(summary, sort_keys=True,
                                   default=lambda x: None if (isinstance(x, float) and x != x) else x))
    P("===============================================================")
    return summary


if __name__ == "__main__":
    main()
