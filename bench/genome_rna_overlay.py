#!/usr/bin/env python
"""NON-CIRCULAR cross-modal corroboration of O1.

Question: do the GENOME-grounded multi-copy families (SEDEF segdup pairs + NCBI annotation,
built WITHOUT ever reading a read) and the RNA read-conflict families (built from the BAM only)
agree on which loci are paralogous copies?

Non-circularity: the RNA family catalog (denovo_families.tsv) is produced by twopass_denovo ->
denovo_assemble_gate -> denovo_families (POA/k-mer over transcript sequences). It NEVER reads
final.bed (SEDEF) or the GFF. SEDEF self-alignment never sees a read. The only shared substrate
is the bare reference genome (allowed). So a SEDEF-link agreement is two independent estimators of
the same biology converging, NOT a data leak.

METRICS (exactly as specified by the task methodology):
  1. HEADLINE  = annotation-free SEDEF-linkage precision. For each RNA family with >=2 mapped
                 distinct loci, for each within-family locus PAIR, test whether a SEDEF segdup pair
                 links the two loci's GENOMIC REGIONS (one side overlaps locus1, the other side
                 overlaps locus2). precision = #linked-pairs / #pairs-tested.  Capped at 100k pairs.
  2. NULL      = distance-stratified random pairs of de-novo loci (NOT co-members), matched on
                 chromosome and genomic separation. background = their SEDEF-linkage rate.
                 LIFT = precision/background, with Fisher + binomial significance + permutation p.
  3. RECALL    = conditioned on expression: of genome protein-coding families (>=2 genes) with >=2
                 RNA-EXPRESSED member genes, fraction where >=2 of those expressed genes land in ONE
                 RNA family.
  4. SECONDARY = annotation-based precision: within-RNA-family loci -> genes -> share a
                 genome_family_def family (GFAM).

Coordinates: meta and final.bed are BOTH 0-based half-open (no conversion). GFF is 1-based inclusive
(converted on load by genome_family_def.load_genes: start-1, end).

Run: /home/juanfra/miniforge3/bin/python bench/genome_rna_overlay.py
"""
import json
import os
import random
import re
import sys
from bisect import bisect_left, bisect_right
from collections import defaultdict

from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from genome_family_def import load_genes, load_sedef_pairs

BENCH = os.path.dirname(os.path.abspath(__file__))
META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
GFF = "/home/juanfra/winloci_scratch/GGO_genomic.gff"
SEDEF = next((p for p in ["/mnt/c/Users/jfris/Desktop/final.bed",
                          "/mnt/c/Users/jfris/Desktop/Rustle/final.bed"] if os.path.exists(p)), None)
RNA_FAM = os.path.join(BENCH, "denovo_families.tsv")
GFAM_ANNOT = os.path.join(BENCH, "genome_families_annotated.tsv")
GFAM_PC = os.path.join(BENCH, "genome_families_protein_coding.tsv")

LOCUS_COVER = 0.50      # a locus "overlaps" a segdup side iff >=50% of the locus is covered (pre-registered)
OVL_MODE = "cover50"    # primary overlap rule; sensitivity also runs "midpoint"
DN_GENE_COVER = 0.50    # a locus maps to a gene iff >=50% of the DN span is covered by the gene
PAIR_CAP = 100_000      # global cap on within-family pairs tested
PER_FAM_CAP = 3000      # per-family pair cap (so the DNFAM0 over-merge cannot dominate)
NULL_POOL = 100_000     # matched null pairs for the background rate estimate
PERM_B = 1000           # permutation replicates
PERM_SIZE = 2000        # matched null pairs per permutation replicate
DIST_EDGES = [0, 1e4, 1e5, 1e6, 1e7, float("inf")]   # absolute cis separation bins
RNG = random.Random(20260630)
MITO = {"NC_011120.1"}

ID_RE = re.compile(r"^DN_(N[CW]_\d+\.\d+)_(\d+)_(\d+)$")


# ----------------------------- inputs -----------------------------
def load_meta(path):
    """id -> (chrom, start, end) collapsed to UNION span over duplicate rows. Also per-contig universe."""
    span = {}     # id -> [chrom, min_start, max_end]
    with open(path) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            tid, chrom, s, e = f[0], f[1], int(f[2]), int(f[3])
            if chrom in MITO:
                continue
            if tid in span:
                r = span[tid]
                if s < r[1]:
                    r[1] = s
                if e > r[2]:
                    r[2] = e
            else:
                span[tid] = [chrom, s, e]
    out = {tid: (r[0], r[1], r[2]) for tid, r in span.items()}
    return out


def load_rna_families(path):
    fams = []
    with open(path) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            fams.append((f[0], f[2].split(",")))
    return fams


# ----------------------------- SEDEF interval index -----------------------------
def _ovl_ok(ss, se, s, e, mode):
    """does segdup side [ss,se] 'cover' locus [s,e] under `mode`?"""
    if mode == "midpoint":
        mid = (s + e) // 2
        return ss <= mid < se
    ov = min(se, e) - max(ss, s)            # cover50
    return ov > 0 and ov >= LOCUS_COVER * (e - s)


class SedefIndex:
    """per-contig sorted segdup SIDES, each carrying its partner side."""

    def __init__(self, pairs):
        sides = defaultdict(list)   # contig -> list of (start, end, partner_contig, partner_start, partner_end)
        for (cA, sA, eA, cB, sB, eB) in pairs:
            if cA in MITO or cB in MITO:
                continue
            sides[cA].append((sA, eA, cB, sB, eB))
            sides[cB].append((sB, eB, cA, sA, eA))
        self.starts = {}
        self.rows = {}
        for c, lst in sides.items():
            lst.sort()
            self.rows[c] = lst
            self.starts[c] = [r[0] for r in lst]

    def _overlapping_sides(self, chrom, s, e, mode):
        """sides on `chrom` that 'cover' locus [s,e] under `mode`."""
        rows = self.rows.get(chrom)
        if not rows:
            return []
        starts = self.starts[chrom]
        hi = bisect_left(starts, e)   # sides starting before locus end
        out = []
        # walk left; sides are start-sorted. A side can still overlap if its end > s.
        for j in range(hi - 1, -1, -1):
            ss, se = rows[j][0], rows[j][1]
            if se <= s:
                # start-sorted: once we are far enough left that even long sides can't reach, stop.
                if s - ss > 6_000_000:
                    break
                continue
            if _ovl_ok(ss, se, s, e, mode):
                out.append(rows[j])
        return out

    def partner_regions(self, chrom, s, e, mode):
        """merged partner intervals (the OPPOSITE segdup sides) for a locus, as dict contig->(starts[],rows[])."""
        by_c = defaultdict(list)
        for (_ss, _se, pc, ps, pe) in self._overlapping_sides(chrom, s, e, mode):
            by_c[pc].append((ps, pe))
        idx = {}
        for c, ivals in by_c.items():
            ivals.sort()
            merged = []
            cs = ce = None
            for a, b in ivals:
                if cs is None:
                    cs, ce = a, b
                elif a <= ce:
                    ce = max(ce, b)
                else:
                    merged.append((cs, ce)); cs, ce = a, b
            if cs is not None:
                merged.append((cs, ce))
            idx[c] = ([m[0] for m in merged], merged)
        return idx


# ----------------------------- generic interval index (genes) -----------------------------
def build_gene_index(genes, genes_by):
    starts = {c: [x[0] for x in v] for c, v in genes_by.items()}
    return genes_by, starts


def genes_covering_dn(genes_by, gstarts, chrom, s, e):
    """gene indices whose overlap covers >=DN_GENE_COVER of the DN span [s,e]; ranked desc by coverage."""
    ivs = genes_by.get(chrom)
    if not ivs:
        return []
    starts = gstarts[chrom]
    span = e - s
    if span <= 0:
        return []
    need = DN_GENE_COVER * span
    hi = bisect_left(starts, e)
    hits = []
    for j in range(hi - 1, -1, -1):
        gs, ge, gi = ivs[j]
        if ge <= s:
            if s - gs > 6_000_000:
                break
            continue
        ov = min(ge, e) - max(gs, s)
        if ov > 0 and ov >= need:
            hits.append((ov, gi))
    hits.sort(reverse=True)
    return [gi for _ov, gi in hits]


def genes_overlapping_interval(genes_by, gstarts, chrom, s, e):
    """ALL gene indices that overlap [s,e] at all (used to test gene expression by de-novo loci)."""
    ivs = genes_by.get(chrom)
    if not ivs:
        return []
    starts = gstarts[chrom]
    hi = bisect_left(starts, e)
    out = []
    for j in range(hi - 1, -1, -1):
        gs, ge, gi = ivs[j]
        if ge <= s:
            if s - gs > 6_000_000:
                break
            continue
        if min(ge, e) - max(gs, s) > 0:
            out.append(gi)
    return out


# ----------------------------- link test -----------------------------
def interval_covers_locus(idx_for_contig, s, e, mode):
    """does any merged partner interval 'cover' locus [s,e] under `mode` ?"""
    if idx_for_contig is None:
        return False
    starts, rows = idx_for_contig
    hi = bisect_right(starts, e)
    for j in range(hi - 1, -1, -1):
        ss, se = rows[j]
        if se <= s:
            if s - ss > 6_000_000:
                break
            continue
        if _ovl_ok(ss, se, s, e, mode):
            return True
    return False


def make_linker(sedef, mode):
    """link operating directly on intervals (c,s,e); cache keyed by the interval tuple."""
    cache = {}

    def partner(iv):
        pr = cache.get(iv)
        if pr is None:
            c, s, e = iv
            pr = sedef.partner_regions(c, s, e, mode)
            cache[iv] = pr
        return pr

    def link(iv_i, iv_j):
        ci, si, ei = iv_i
        cj, sj, ej = iv_j
        pr_i = partner(iv_i)
        if interval_covers_locus(pr_i.get(cj), sj, ej, mode):
            return True
        pr_j = partner(iv_j)            # symmetric safety (redundant in principle)
        return interval_covers_locus(pr_j.get(ci), si, ei, mode)

    return link, cache


# ----------------------------- strata helpers -----------------------------
def dist_bin(d):
    for k in range(len(DIST_EDGES) - 1):
        if DIST_EDGES[k] <= d < DIST_EDGES[k + 1]:
            return k
    return len(DIST_EDGES) - 2


def _strat_lift(o, n):
    """observed/null rate ratio; 'inf' when observed>0 but null rate is 0."""
    op = o[0] / o[1] if o[1] else 0.0
    npr = n[0] / n[1] if n[1] else 0.0
    if npr > 0:
        return round(op / npr, 2)
    return float("inf") if op > 0 else 0.0


def main():
    assert SEDEF, "final.bed not found"
    print(f"[load] meta   {META}")
    meta = load_meta(META)
    print(f"       {len(meta)} unique de-novo loci (mito excluded, union-span collapsed)")
    print(f"[load] sedef  {SEDEF}")
    pairs = load_sedef_pairs(SEDEF)
    print(f"       {len(pairs)} segdup pairs -> building side index")
    sedef = SedefIndex(pairs)
    print(f"[load] genes  {GFF}")
    genes, genes_by = load_genes(GFF)
    genes_by, gstarts = build_gene_index(genes, genes_by)
    print(f"       {len(genes)} genes")
    print(f"[load] rna families {RNA_FAM}")
    rna_fams = load_rna_families(RNA_FAM)
    print(f"       {len(rna_fams)} RNA families, {sum(len(m) for _, m in rna_fams)} members")

    # ---- universe of de-novo loci, per contig sorted (for null draws) ----
    uni_by_c = defaultdict(list)     # contig -> list of (start, id)
    for tid, (c, s, e) in meta.items():
        uni_by_c[c].append((s, tid))
    for c in uni_by_c:
        uni_by_c[c].sort()
    uni_starts = {c: [x[0] for x in v] for c, v in uni_by_c.items()}

    # locus -> rna family (for excluding co-members in null + recall)
    locus_fam = {}
    for fid, members in rna_fams:
        for m in members:
            locus_fam[m] = fid

    # =====================================================================
    # Build within-family DISTINCT LOCI (merge overlapping member intervals so
    # isoform-splits of the same locus are not counted as separate "copies").
    # Pairs / strata are MODE-INDEPENDENT; only the SEDEF link result depends on mode.
    # =====================================================================
    def merge_loci(intervals):
        by_c = defaultdict(list)
        for (c, s, e) in intervals:
            by_c[c].append((s, e))
        loci = []
        for c, ivl in by_c.items():
            ivl.sort()
            cs = ce = None
            for s, e in ivl:
                if cs is None:
                    cs, ce = s, e
                elif s <= ce:               # overlap -> same locus
                    ce = max(ce, e)
                else:
                    loci.append((c, cs, ce)); cs, ce = s, e
            if cs is not None:
                loci.append((c, cs, ce))
        return loci

    n_missing = 0
    per_fam = []          # (fid, n_loci_distinct, n_pairs_possible, n_pairs_sampled)
    fam_pairs = {}        # fid -> list of (iv_i, iv_j) interval pairs (capped)
    for fid, members in rna_fams:
        ivs = []
        for m in members:
            if m not in meta:
                n_missing += 1
                continue
            ivs.append(meta[m])
        loci = merge_loci(ivs)
        n_loci = len(loci)
        if n_loci < 2:
            per_fam.append((fid, n_loci, 0, 0))
            continue
        allp = [(loci[a], loci[b]) for a in range(n_loci) for b in range(a + 1, n_loci)]
        npos = len(allp)
        if npos > PER_FAM_CAP:
            allp = RNG.sample(allp, PER_FAM_CAP)
        fam_pairs[fid] = allp
        per_fam.append((fid, n_loci, npos, len(allp)))

    total_sampled = sum(len(p) for p in fam_pairs.values())
    if total_sampled > PAIR_CAP:
        scale = PAIR_CAP / total_sampled
        fam_pairs = {fid: RNG.sample(pl, min(max(1, int(round(len(pl) * scale))), len(pl)))
                     for fid, pl in fam_pairs.items()}

    # flat strata list (mode-independent) for matched-null sampling
    def stratum_of(iv_i, iv_j):
        ci, si, _ = iv_i; cj, sj, _ = iv_j
        if ci != cj:
            return ("trans", ci, cj)
        return ("cis", ci, dist_bin(abs(si - sj)))
    obs_strata = [stratum_of(a, b) for pl in fam_pairs.values() for (a, b) in pl]

    # ---- null draw functions (mode-independent; return interval pairs) ----
    def draw_cis(contig, bin_k):
        starts = uni_starts.get(contig)
        if not starts or len(starts) < 2:
            return None
        lo, hi = DIST_EDGES[bin_k], DIST_EDGES[bin_k + 1]
        for _ in range(20):
            ai = RNG.randrange(len(starts))
            sa = starts[ai]
            r0 = bisect_left(starts, sa + lo)
            r1 = bisect_left(starts, sa + hi) if hi != float("inf") else len(starts)
            l1 = bisect_right(starts, sa - lo)
            l0 = bisect_right(starts, sa - hi) if hi != float("inf") else 0
            cand = [c for c in (list(range(r0, r1)) + list(range(l0, l1))) if c != ai]
            if not cand:
                continue
            bi = RNG.choice(cand)
            ta, tb = uni_by_c[contig][ai][1], uni_by_c[contig][bi][1]
            if ta == tb:
                continue
            if locus_fam.get(ta) is not None and locus_fam.get(ta) == locus_fam.get(tb):
                continue
            return (meta[ta], meta[tb])
        return None

    def draw_trans(c1, c2):
        l1, l2 = uni_by_c.get(c1), uni_by_c.get(c2)
        if not l1 or not l2:
            return None
        for _ in range(10):
            ta = l1[RNG.randrange(len(l1))][1]
            tb = l2[RNG.randrange(len(l2))][1]
            if ta == tb:
                continue
            if locus_fam.get(ta) is not None and locus_fam.get(ta) == locus_fam.get(tb):
                continue
            return (meta[ta], meta[tb])
        return None

    def draw_matched(st):
        return draw_trans(st[1], st[2]) if st[0] == "trans" else draw_cis(st[1], st[2])

    # =====================================================================
    # eval_mode: METRIC 1 (headline) + METRIC 2 (null) for a given overlap mode
    # =====================================================================
    def eval_mode(mode, full=True):
        link, _ = make_linker(sedef, mode)
        fam_rows, strat_obs = [], defaultdict(lambda: [0, 0])
        n_linked = n_tested = nl_x = nt_x = fam_coherent = fam_with_pairs = 0
        for fid, plist in fam_pairs.items():
            nlinked_f = 0
            for (a, b) in plist:
                r = link(a, b)
                n_tested += 1
                kind = stratum_of(a, b)
                kk = kind[0] if kind[0] == "trans" else f"cis_bin{kind[2]}"
                strat_obs[kk][1] += 1
                if r:
                    strat_obs[kk][0] += 1
                    n_linked += 1; nlinked_f += 1
                if fid != "DNFAM0":
                    nt_x += 1
                    if r:
                        nl_x += 1
            fam_with_pairs += 1
            fam_coherent += int(nlinked_f > 0)
            nloci = next(p[1] for p in per_fam if p[0] == fid)
            fam_rows.append((fid, nloci, len(plist), nlinked_f, int(nlinked_f > 0)))
        for fid, nloci, npos, ns in per_fam:
            if fid not in fam_pairs:
                fam_rows.append((fid, nloci, 0, 0, 0))
        precision = n_linked / n_tested if n_tested else 0.0
        precision_x = nl_x / nt_x if nt_x else 0.0

        # null
        npool = NULL_POOL if full else NULL_POOL // 2
        null_linked = null_tested = 0
        strat_null = defaultdict(lambda: [0, 0])
        for _ in range(npool):
            st = obs_strata[RNG.randrange(len(obs_strata))]
            pr = draw_matched(st)
            if pr is None:
                continue
            r = link(pr[0], pr[1])
            null_tested += 1
            kk = st[0] if st[0] == "trans" else f"cis_bin{st[2]}"
            strat_null[kk][1] += 1
            if r:
                strat_null[kk][0] += 1; null_linked += 1
        background = null_linked / null_tested if null_tested else 0.0
        lift = (precision / background) if background > 0 else float("inf")

        res = dict(mode=mode, link=link, fam_rows=fam_rows,
                   n_linked=n_linked, n_tested=n_tested, precision=precision,
                   nl_x=nl_x, nt_x=nt_x, precision_x=precision_x,
                   fam_coherent=fam_coherent, fam_with_pairs=fam_with_pairs,
                   null_linked=null_linked, null_tested=null_tested,
                   background=background, lift=lift,
                   strat_obs=dict(strat_obs), strat_null=dict(strat_null))
        return res

    print(f"\n[1+2] running primary overlap mode = {OVL_MODE} ...")
    R = eval_mode(OVL_MODE, full=True)
    link = R["link"]
    fam_rows = R["fam_rows"]
    n_linked, n_tested, precision = R["n_linked"], R["n_tested"], R["precision"]
    nl_x, nt_x, precision_no_dnfam0 = R["nl_x"], R["nt_x"], R["precision_x"]
    fam_coherent, fam_with_pairs = R["fam_coherent"], R["fam_with_pairs"]
    fam_coherent_frac = fam_coherent / fam_with_pairs if fam_with_pairs else 0.0
    background, lift = R["background"], R["lift"]
    null_linked, null_tested = R["null_linked"], R["null_tested"]
    strat_obs, strat_null = R["strat_obs"], R["strat_null"]

    print(f"\n[1] HEADLINE precision (SEDEF-linked within-family DISTINCT-locus pairs):")
    print(f"    pairs possible (sum over fams, pre-cap)= {sum(f[2] for f in per_fam)}")
    print(f"    pairs tested (after caps)             = {n_tested}")
    print(f"    SEDEF-linked                          = {n_linked}")
    print(f"    precision                             = {precision:.4f}")
    print(f"    precision EXCLUDING DNFAM0 over-merge  = {precision_no_dnfam0:.4f}  ({nl_x}/{nt_x})")
    print(f"    family-equal coherent fraction (>=1 linked pair) = {fam_coherent_frac:.4f} "
          f"({fam_coherent}/{fam_with_pairs})")

    print(f"\n[2] NULL CONTROL (distance/chrom-matched, {null_tested} draws):")
    print(f"    background SEDEF-link rate = {background:.5f}  ({null_linked}/{null_tested})")
    print(f"    LIFT (precision/background) = {lift:.2f}")

    table = [[n_linked, n_tested - n_linked], [null_linked, null_tested - null_linked]]
    try:
        odds, fisher_p = stats.fisher_exact(table, alternative="greater")
    except Exception:
        fisher_p = float("nan")
    binom_p = stats.binomtest(n_linked, n_tested, background, alternative="greater").pvalue \
        if (background > 0 and n_tested) else float("nan")
    print(f"    Fisher exact (greater) p = {fisher_p:.3e};  binomial p = {binom_p:.3e}")

    # permutation p (primary mode)
    perm_ge = 0
    perm_rates = []
    sz = min(PERM_SIZE, len(obs_strata)) if obs_strata else 0
    if sz:
        for _ in range(PERM_B):
            kk = tt = 0
            for _ in range(sz):
                st = obs_strata[RNG.randrange(len(obs_strata))]
                pr = draw_matched(st)
                if pr is None:
                    continue
                tt += 1
                if link(pr[0], pr[1]):
                    kk += 1
            rate = kk / tt if tt else 0.0
            perm_rates.append(rate)
            if rate >= precision:
                perm_ge += 1
    perm_p = (perm_ge + 1) / (PERM_B + 1) if perm_rates else float("nan")
    print(f"    permutation p (B={PERM_B}, size={sz}) = {perm_p:.4g}"
          + (f"  (max null rate {max(perm_rates):.4f})" if perm_rates else ""))

    print(f"    per-stratum (observed | null):")
    for kind in sorted(set(list(strat_obs) + list(strat_null))):
        o = strat_obs.get(kind, [0, 0]); nn = strat_null.get(kind, [0, 0])
        op = o[0] / o[1] if o[1] else 0.0
        npr = nn[0] / nn[1] if nn[1] else 0.0
        lf = op / npr if npr > 0 else float("inf")
        print(f"      {kind:10s} obs {op:.3f} ({o[0]}/{o[1]})  null {npr:.4f} ({nn[0]}/{nn[1]})  lift {lf:.1f}")

    # ---- sensitivity: permissive midpoint-in-block overlap rule ----
    print(f"\n[1b] SENSITIVITY overlap rule = midpoint-in-block:")
    Rm = eval_mode("midpoint", full=False)
    print(f"    precision {Rm['precision']:.4f} (excl DNFAM0 {Rm['precision_x']:.4f}); "
          f"background {Rm['background']:.5f}; LIFT {Rm['lift']:.1f}")

    # =====================================================================
    # METRIC 3 : CONDITIONED RECALL (genome PC families recovered as RNA families)
    # =====================================================================
    # gene symbol -> intervals (protein_coding genes only)
    sym_to_iv = defaultdict(list)
    for gi, g in enumerate(genes):
        if g["biotype"] == "protein_coding":
            sym_to_iv[g["name"]].append((g["chrom"], g["start"], g["end"]))

    # universe de-novo loci -> for each gene, which de-novo loci overlap it (expression);
    # and which of those loci are RNA-family members (-> rna family).
    # Build gene-interval expression by scanning universe loci against gene index.
    gene_expr_loci = defaultdict(list)    # gene_idx -> list of de-novo locus ids overlapping
    for tid, (c, s, e) in meta.items():
        for gi in genes_overlapping_interval(genes_by, gstarts, c, s, e):
            gene_expr_loci[gi].append(tid)

    # map gene symbol -> set of rna families it lands in (via expressed member loci)
    def sym_rna_families(sym):
        fams = set()
        expressed = False
        for (c, s, e) in sym_to_iv.get(sym, []):
            # gene idx for this interval: find matching gene index
            pass
        return fams, expressed

    # build symbol -> gene indices
    sym_to_gidx = defaultdict(list)
    for gi, g in enumerate(genes):
        if g["biotype"] == "protein_coding":
            sym_to_gidx[g["name"]].append(gi)

    def sym_expressed_and_fams(sym):
        fams = set()
        expressed = False
        for gi in sym_to_gidx.get(sym, []):
            loci = gene_expr_loci.get(gi, [])
            if loci:
                expressed = True
            for tid in loci:
                f = locus_fam.get(tid)
                if f is not None:
                    fams.add(f)
        return expressed, fams

    # load genome PC families (member gene symbols)
    pc_fams = []
    with open(GFAM_PC) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            pc_fams.append((f[0], f[6].split(",")))

    n_pc = len(pc_fams)
    n_pc_expr2 = 0           # >=2 member genes expressed
    n_pc_recovered = 0       # >=2 expressed genes share ONE rna family
    recall_rows = []
    for gid, syms in pc_fams:
        expr_genes = 0
        fam_counter = defaultdict(int)
        for sym in set(syms):
            ex, fams = sym_expressed_and_fams(sym)
            if ex:
                expr_genes += 1
            for fam in fams:
                fam_counter[fam] += 1
        if expr_genes >= 2:
            n_pc_expr2 += 1
            recovered = any(v >= 2 for v in fam_counter.values())
            if recovered:
                n_pc_recovered += 1
            recall_rows.append((gid, expr_genes, recovered))
    cond_recall = n_pc_recovered / n_pc_expr2 if n_pc_expr2 else 0.0
    print(f"\n[3] CONDITIONED RECALL (genome protein-coding families):")
    print(f"    total PC families (>=2 genes)            = {n_pc}")
    print(f"    with >=2 RNA-EXPRESSED member genes       = {n_pc_expr2}")
    print(f"    of those, >=2 expressed genes in ONE RNA family = {n_pc_recovered}")
    print(f"    conditioned recall                        = {cond_recall:.4f}")

    # =====================================================================
    # METRIC 4 : SECONDARY annotation-based precision (loci -> genes -> share GFAM)
    # =====================================================================
    # gene symbol -> GFAM id (annotated catalog, broad)
    sym_gfam = {}
    with open(GFAM_ANNOT) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            fam = f[0]
            for sym in f[6].split(","):
                sym_gfam.setdefault(sym, set()).add(fam)

    # locus(interval) -> gene symbol set (>=50% DN span coverage)
    locus_genes = {}
    def genes_for_locus(iv):
        if iv in locus_genes:
            return locus_genes[iv]
        c, s, e = iv
        gidx = genes_covering_dn(genes_by, gstarts, c, s, e)
        syms = [genes[gi]["name"] for gi in gidx]
        locus_genes[iv] = syms
        return syms

    annot_tested = annot_linked = 0
    annot_both_annotated = 0
    for fid, plist in fam_pairs.items():
        for (i, j) in plist:
            gi = genes_for_locus(i); gj = genes_for_locus(j)
            if not gi or not gj:
                continue
            annot_both_annotated += 1
            fams_i = set()
            for s in gi:
                fams_i |= sym_gfam.get(s, set())
            shared = False
            for s in gj:
                if sym_gfam.get(s, set()) & fams_i:
                    # require distinct genes
                    if s not in gi or len(set(gi) | set(gj)) >= 2:
                        shared = True
                        break
            annot_tested += 1
            if shared:
                annot_linked += 1
    annot_precision = annot_linked / annot_both_annotated if annot_both_annotated else 0.0
    print(f"\n[4] SECONDARY annotation-based precision (loci->genes->share GFAM):")
    print(f"    pairs with genes on BOTH loci = {annot_both_annotated} / {n_tested} tested")
    print(f"    of those, share a GFAM family = {annot_linked}")
    print(f"    annotation precision          = {annot_precision:.4f}")

    # =====================================================================
    # write per-family TSV + JSON
    # =====================================================================
    fam_rows.sort(key=lambda r: (-r[2], r[0]))
    tsv_path = os.path.join(BENCH, "genome_rna_overlay_per_family.tsv")
    with open(tsv_path, "w") as fh:
        fh.write("family_id\tn_loci_mapped\tn_pairs\tn_sedef_linked\tcoherent_flag\n")
        for (fid, nloci, npairs, nlinked, coh) in fam_rows:
            fh.write(f"{fid}\t{nloci}\t{npairs}\t{nlinked}\t{coh}\n")

    summary = dict(
        non_circular="RNA families from BAM POA/k-mer; SEDEF from genome self-alignment; "
                     "only shared substrate = bare reference genome",
        n_rna_families=len(rna_fams),
        n_rna_members=sum(len(m) for _, m in rna_fams),
        n_members_missing_meta=n_missing,
        # metric 1
        headline_precision=precision,
        headline_n_linked=n_linked,
        headline_n_tested=n_tested,
        pairs_possible_precap=sum(f[2] for f in per_fam),
        pair_cap=PAIR_CAP,
        per_family_cap=PER_FAM_CAP,
        precision_excl_DNFAM0=precision_no_dnfam0,
        precision_excl_DNFAM0_counts=[nl_x, nt_x],
        family_equal_coherent_fraction=fam_coherent_frac,
        family_coherent_counts=[fam_coherent, fam_with_pairs],
        locus_cover_rule=LOCUS_COVER,
        per_stratum_precision={k: (v[0], v[1]) for k, v in strat_obs.items()},
        # metric 2
        null_background=background,
        null_counts=[null_linked, null_tested],
        lift=lift,
        fisher_p_greater=float(fisher_p),
        binomial_p_greater=float(binom_p),
        permutation_p=perm_p,
        permutation_B=PERM_B,
        permutation_size=sz,
        per_stratum_null={k: (v[0], v[1]) for k, v in strat_null.items()},
        per_stratum_lift={k: _strat_lift(strat_obs.get(k, [0, 0]), strat_null.get(k, [0, 0]))
                          for k in set(list(strat_obs) + list(strat_null))},
        sensitivity_midpoint_rule=dict(precision=Rm["precision"], precision_excl_DNFAM0=Rm["precision_x"],
                                       background=Rm["background"], lift=Rm["lift"],
                                       counts=[Rm["n_linked"], Rm["n_tested"]],
                                       null_counts=[Rm["null_linked"], Rm["null_tested"]]),
        # metric 3
        pc_families_total=n_pc,
        pc_families_expr_ge2=n_pc_expr2,
        pc_families_recovered=n_pc_recovered,
        conditioned_recall=cond_recall,
        # metric 4
        annot_precision=annot_precision,
        annot_counts=[annot_linked, annot_both_annotated],
        # reconcile (prior)
        prior_reconcile=dict(
            prior_precision_raw=0.644, prior_precision_eff=0.719, prior_recall=0.105,
            note="prior used cDNA-homology truth (partly circular, shares lineage with reads); "
                 "this overlay replaces it with read-independent SEDEF segdup linkage.",
        ),
    )
    json_path = os.path.join(BENCH, "genome_rna_overlay.json")
    json.dump(summary, open(json_path, "w"), indent=2)

    # ---- self-consistency checks ----
    assert n_linked <= n_tested
    assert null_linked <= null_tested
    assert 0.0 <= precision <= 1.0 and 0.0 <= background <= 1.0
    assert abs(sum(r[3] for r in fam_rows) - n_linked) == 0, "per-family linked sum != headline linked"
    assert sum(r[2] for r in fam_rows) == n_tested, "per-family pair sum != tested"

    print(f"\n[ok] wrote {tsv_path}")
    print(f"     wrote {json_path}")
    print(f"\n== SELF-CONSISTENCY: per-family linked sum {sum(r[3] for r in fam_rows)} == headline {n_linked}; "
          f"per-family pairs {sum(r[2] for r in fam_rows)} == tested {n_tested} ==")
    print(f"\n== SUMMARY ==")
    print(f"  HEADLINE precision      {precision:.3f}  ({n_linked}/{n_tested}); excl DNFAM0 {precision_no_dnfam0:.3f}")
    print(f"  NULL background         {background:.4f};  LIFT {lift:.1f}x  (Fisher p {fisher_p:.1e}, perm p {perm_p})")
    print(f"  CONDITIONED RECALL      {cond_recall:.3f}  ({n_pc_recovered}/{n_pc_expr2} expressed PC families)")
    print(f"  ANNOTATION precision    {annot_precision:.3f}  ({annot_linked}/{annot_both_annotated})")


if __name__ == "__main__":
    main()
