#!/usr/bin/env python3
"""Locus-construction probes (§6w1-§6w5), consolidated from four one-off scripts so the register rows they
produced (r969-r978) keep a runnable, ANCHORED generator without four files. Each subcommand is the original
script verbatim (only `main` renamed and colliding helpers prefixed per block); the original docstrings follow.

  locus_probes.py colinearity ...           was bench/colinearity_conjunct.py        (r969)
  locus_probes.py five-prime ...            was bench/five_prime_deficit.py          (r973-r975)
  locus_probes.py partner-discontinuity ... was bench/partner_discontinuity.py       (§6w2 finding 3)
  locus_probes.py soto-fidelity ...         was bench/soto_family_locus_fidelity.py  (r976-r978)
"""
import sys
import argparse
import collections
import csv
import gzip
import itertools
import re
import statistics


# ================================================================================================
# colinearity  (was bench/colinearity_conjunct.py)  renamed: main->main_colinearity_conjunct
# ================================================================================================
"""Does block colinearity separate real containments from repeat/domain ones?
Per `docs/PREREG_colinearity_conjunct_2026-09-21.md`.

r919 swept containment/identity from 0.30 to 0.99 on the 76 pairs the shipped rule rejects at
containment >= 0.90 (19 TRUE / 57 FALSE, held-out chr2/chr8/chr10) and found a FLAT precision curve
(~0.25) -- no pairwise scalar separates them. Every earlier pairwise script keeps only the BEST PAF
record per gene pair; 95.6% of chr2 gene pairs have more than one record (minimap2 run with -P), so
this is genuinely unused information.

    colinearity(pair) = Kendall-tau concordance of block order between A's and B's coordinates,
                         over every surviving PAF record for that pair (not just the best one)

Usage: colinearity_conjunct.py --gff G --soto S1C --pafs chr2=P1,chr8=P2,chr10=P3 --graphs DIR"""
MIN_BLOCK_NMATCH = 30

def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)

def gene_names(gff, chrom):
    n = {}
    for line in _open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] not in ('gene', 'pseudogene'):
            continue
        m = re.search(r'Name=([^;]+)', f[8])
        if m:
            n[f'{f[0]}:{f[3]}-{f[4]}'] = m.group(1)
    return n

def soto(s1c):
    ids = collections.defaultdict(set)
    for r in csv.DictReader(open(s1c), delimiter='\t'):
        fid = (r.get('Family ID') or '').strip(); nm = (r.get('Gene Name') or '').strip()
        if fid and fid != 'N/A' and nm:
            ids[nm].add(fid)
    return {nm: next(iter(v)) for nm, v in ids.items() if len(v) == 1}

def load_all_blocks(paf):
    """(a,b) sorted tuple -> list of (a_start,a_end,b_start,b_end,strand) in a's/b's own coords,
    plus the best (m, la, lb) per pair for the containment gate (matches vg_multiplicity_containment.py)."""
    blocks = collections.defaultdict(list)
    best = {}
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 12:
            continue
        a, b = f[0], f[5]
        if a == b:
            continue
        la, lb, m = int(f[1]), int(f[6]), int(f[9])
        qs, qe, ts, te = int(f[2]), int(f[3]), int(f[7]), int(f[8])
        strand = f[4]
        k = tuple(sorted((a, b)))
        if m >= MIN_BLOCK_NMATCH:
            if k[0] == a:
                blocks[k].append((qs, qe, ts, te, strand))
            else:
                blocks[k].append((ts, te, qs, qe, strand))
        rec = (m, la, lb, qs, qe, ts, te) if k[0] == a else (m, lb, la, ts, te, qs, qe)
        if k not in best or m > best[k][0]:
            best[k] = rec
    return best, blocks

def colinearity(pairblocks):
    """Kendall-tau concordance of block order; None if <2 blocks or mixed strand."""
    if len(pairblocks) < 2:
        return None, len(pairblocks)
    strands = {s for (_, _, _, _, s) in pairblocks}
    if len(strands) > 1:
        return 0.0, len(pairblocks)
    strand = strands.pop()
    ordered = sorted(pairblocks, key=lambda x: x[0])
    bstarts = [x[2] for x in ordered]
    pairs = list(itertools.combinations(bstarts, 2))
    if not pairs:
        return None, len(pairblocks)
    if strand == '+':
        concordant = sum(1 for x, y in pairs if x <= y)
    else:
        concordant = sum(1 for x, y in pairs if x >= y)
    return concordant / len(pairs), len(pairblocks)

def main_colinearity_conjunct():
    ap = argparse.ArgumentParser()
    for x in ('--gff', '--soto', '--pafs', '--graphs'):
        ap.add_argument(x, required=True)
    a = ap.parse_args()
    fam = soto(a.soto)

    rows = []
    for spec in a.pafs.split(','):
        chrom, paf = spec.split('=')
        names = gene_names(a.gff, chrom)
        shipped = set()
        for line in open(f'{a.graphs}/{chrom}.graph.tsv'):
            f = line.rstrip('\n').split('\t')
            if len(f) == 3 and f[0] != f[1]:
                shipped.add(tuple(sorted((f[0], f[1]))))
        best, blocks = load_all_blocks(paf)
        for k, (m, la, lb, qs, qe, ts, te) in best.items():
            if k in shipped:
                continue
            ga, gb = names.get(k[0]), names.get(k[1])
            if not ga or not gb:
                continue
            fa, fb = fam.get(ga), fam.get(gb)
            if not fa or not fb:
                continue
            lo = min(la, lb)
            if not lo or m / lo < 0.90:
                continue
            score, nblocks = colinearity(blocks.get(k, []))
            rows.append((fa == fb, score, nblocks, ga, gb))

    T = [r for r in rows if r[0]]
    F = [r for r in rows if not r[0]]
    print(f"population: {len(rows)} rejected containment>=0.90 pairs — {len(T)} TRUE, {len(F)} FALSE\n")

    Tb = [r[2] for r in T]; Fb = [r[2] for r in F]
    print(f"  block count   TRUE median {statistics.median(Tb) if Tb else 0:.1f}"
          f"   FALSE median {statistics.median(Fb) if Fb else 0:.1f}")
    Tunder = sum(1 for r in T if r[2] < 2); Funder = sum(1 for r in F if r[2] < 2)
    print(f"  <2 blocks (colinearity undefined): TRUE {Tunder}/{len(T)} ({100*Tunder/len(T) if T else 0:.1f}%)"
          f"   FALSE {Funder}/{len(F)} ({100*Funder/len(F) if F else 0:.1f}%)\n")

    Ts = [r[1] for r in T if r[1] is not None]
    Fs = [r[1] for r in F if r[1] is not None]
    print(f"  colinearity (n={len(Ts)+len(Fs)} scoreable)  TRUE median "
          f"{statistics.median(Ts) if Ts else float('nan'):.3f}   FALSE median "
          f"{statistics.median(Fs) if Fs else float('nan'):.3f}\n")

    print(f"  {'score >= t':>10} {'TRUE kept':>10} {'FALSE kept':>11} {'precision':>10} {'recall':>8}")
    best_p = 0.0
    for t in (0.50, 0.70, 0.80, 0.90, 0.95, 0.99, 1.00):
        t_kept = sum(1 for r in T if r[1] is not None and r[1] >= t)
        f_kept = sum(1 for r in F if r[1] is not None and r[1] >= t)
        p = t_kept / (t_kept + f_kept) if (t_kept + f_kept) else 0.0
        if t_kept >= 10:
            best_p = max(best_p, p)
        print(f"  {t:>10.2f} {t_kept:>10} {f_kept:>11} {p:>10.3f} {t_kept/len(T) if T else 0:>8.3f}")

    pairs = list(itertools.product(Ts, Fs))
    wins = sum(1 for x, y in pairs if x > y) + 0.5 * sum(1 for x, y in pairs if x == y)
    auc = wins / len(pairs) if pairs else float('nan')
    print(f"\n  AUC (colinearity separating TRUE from FALSE, scoreable only) = {auc:.3f}"
          f"   [r906's multiplicity AUC was 0.681]")
    print(f"  baseline precision (containment alone) = 0.250")
    print(f"  best precision at >=10 TRUE retained (scoreable-only threshold sweep) = {best_p:.3f}")

# ================================================================================================
# five-prime  (was bench/five_prime_deficit.py)  renamed: main->main_five_prime_deficit
# ================================================================================================
"""How much of a gene's 5' end does the constructed locus miss, and what predicts it?
Per `docs/PREREG_five_prime_deficit_2026-09-21.md` (§6w3).

r806 established that the de novo exon-sum does NOT systematically under-represent a copy's width
(median relative deficit -0.1%), but that "~100% of whatever deficit remains is at the 5' end; the 3'
end is exact (<30bp) on every tool tested". This measures that asymmetry directly, strand-aware, and
tests what predicts the 5' deficit.

⚠ Only loci containing EXACTLY ONE annotated gene are scored: a readthrough-fused locus's "5' deficit"
is undefined (§6v8), so including them would measure over-merge instead of truncation.

⚠ Truth at the 5' end is ambiguous -- RefSeq carries several transcripts per gene -- so BOTH are
reported: against the gene record's own terminus, and against the most extreme 5' end over all of that
gene's transcripts. Neither is picked after the fact.

Usage:
  five_prime_deficit.py --gtf dn16.gtf --gff chr16.genes.gff [--chrom chr16] [--label human-chr16]"""
def load_genes(gff, chrom):
    """name -> (start, end, strand, tx_min_start, tx_max_end) over gene + all its transcripts."""
    gene = {}
    gid_of = {}
    tx_extent = collections.defaultdict(lambda: [None, None])
    lines = [ln for ln in open(gff) if not ln.startswith('#')]
    for ln in lines:
        f = ln.rstrip('\n').split('\t')
        if len(f) < 9 or (chrom and f[0] != chrom):
            continue
        if f[2] in ('gene', 'pseudogene', 'ncRNA_gene'):
            n = re.search(r'Name=([^;]+)', f[8])
            i = re.search(r'ID=([^;]+)', f[8])
            if n and i:
                gene[n.group(1)] = [int(f[3]), int(f[4]), f[6]]
                gid_of[i.group(1)] = n.group(1)
    for ln in lines:
        f = ln.rstrip('\n').split('\t')
        if len(f) < 9 or (chrom and f[0] != chrom):
            continue
        p = re.search(r'Parent=([^;,]+)', f[8])
        if not p or p.group(1) not in gid_of:
            continue
        g = gid_of[p.group(1)]
        s, e = int(f[3]), int(f[4])
        cur = tx_extent[g]
        cur[0] = s if cur[0] is None else min(cur[0], s)
        cur[1] = e if cur[1] is None else max(cur[1], e)
    out = {}
    for g, (s, e, st) in gene.items():
        ts, te = tx_extent.get(g, [s, e])
        out[g] = (s, e, st, min(s, ts or s), max(e, te or e))
    return out

def load_loci(gtf, chrom):
    """locus id -> (start, end, strand, n_exons, reads, exon_list)."""
    tx = {}
    exons = collections.defaultdict(list)
    for ln in open(gtf):
        if ln.startswith('#'):
            continue
        f = ln.rstrip('\n').split('\t')
        if len(f) < 9 or (chrom and f[0] != chrom):
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m:
            continue
        tid = m.group(1)
        if f[2] == 'transcript':
            r = re.search(r'reads "(\d+)"', f[8])
            tx[tid] = [int(f[3]), int(f[4]), f[6], 0, int(r.group(1)) if r else 0]
        elif f[2] == 'exon':
            exons[tid].append((int(f[3]), int(f[4])))
    out = {}
    for tid, v in tx.items():
        ex = sorted(exons.get(tid, []))
        v[3] = len(ex)
        out[tid] = tuple(v) + (ex,)
    return out

def main_five_prime_deficit():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--gff', required=True)
    ap.add_argument('--chrom')
    ap.add_argument('--label', default='substrate')
    a = ap.parse_args()

    genes = load_genes(a.gff, a.chrom)
    loci = load_loci(a.gtf, a.chrom)
    gl = sorted((v[0], v[1], g) for g, v in genes.items())

    rows = []
    for lid, (ls, le, lstrand, nex, reads, ex) in loci.items():
        if nex < 2:
            continue
        inside = [g for (gs, ge, g) in gl
                  if not (ge < ls or gs > le)
                  and (min(le, ge) - max(ls, gs)) >= 0.5 * (ge - gs + 1)]
        if len(inside) != 1:
            continue
        g = inside[0]
        gs, ge, st, uts, ute = genes[g]
        if st == '+':
            d5, d3 = ls - gs, ge - le
            d5u = ls - uts
        else:
            d5, d3 = ge - le, ls - gs
            d5u = ute - le
        rows.append(dict(locus=lid, gene=g, strand=st, reads=reads, nex=nex,
                         span=le - ls, d5=d5, d3=d3, d5u=d5u))

    if not rows:
        print("no scoreable loci")
        return

    d5 = [r['d5'] for r in rows]
    d3 = [r['d3'] for r in rows]
    d5u = [r['d5u'] for r in rows]
    n = len(rows)
    print(f"=== {a.label} ===")
    print(f"scoreable loci (exactly 1 gene inside, >=2 exons): {n}\n")

    def line(tag, v):
        pos = sum(1 for x in v if x > 0)
        print(f"  {tag:26} median {statistics.median(v):>8.0f}   mean {statistics.mean(v):>9.0f}   "
              f"short>0 {100*pos/len(v):>5.1f}%   p75 {sorted(v)[int(.75*len(v))]:>8.0f}   "
              f"p90 {sorted(v)[int(.90*len(v))]:>9.0f}")
    line("d5 (vs gene record)", d5)
    line("d5 (vs transcript union)", d5u)
    line("d3 (vs gene record)", d3)
    print(f"\n  GATE 0  median d5 - median d3 = {statistics.median(d5)-statistics.median(d3):>.0f} bp"
          f"   (bar: >= 100 bp)")

    print("\n  |d5| vs |d3| (magnitude, ignoring direction):")
    print(f"    median |d5| {statistics.median([abs(x) for x in d5]):>8.0f}    "
          f"median |d3| {statistics.median([abs(x) for x in d3]):>8.0f}")
    within = lambda v, t: 100*sum(1 for x in v if abs(x) <= t)/len(v)
    for t in (30, 100, 500):
        print(f"    within +-{t:>4} bp:   5' {within(d5,t):>5.1f}%     3' {within(d3,t):>5.1f}%")

    print("\n  d5 by read depth (is the deficit a coverage effect?):")
    bands = [(2, 2), (3, 4), (5, 9), (10, 29), (30, 10**9)]
    print(f"    {'reads':>10} {'n':>5} {'median d5':>10} {'median d3':>10} {'%short':>7}")
    for lo, hi in bands:
        v = [r for r in rows if lo <= r['reads'] <= hi]
        if len(v) >= 5:
            m5 = statistics.median([r['d5'] for r in v])
            m3 = statistics.median([r['d3'] for r in v])
            sh = 100*sum(1 for r in v if r['d5'] > 0)/len(v)
            lbl = f"{lo}-{hi}" if hi < 10**9 else f">={lo}"
            print(f"    {lbl:>10} {len(v):>5} {m5:>10.0f} {m3:>10.0f} {sh:>6.1f}%")

# ================================================================================================
# partner-discontinuity  (was bench/partner_discontinuity.py)  renamed: main->main_partner_discontinuity
# ================================================================================================
"""Does a readthrough-fused de novo locus show a DISCONTINUITY in its PAF alignment-partner set
along its own coordinate axis?

Every node-split trigger refuted so far (registers 845/846/937-940/948-949/967/968) used READ-level
evidence inside the locus, and all failed the same way: a readthrough molecule is a genuine, abundant,
full-length transcript, so the "bridge" IS the population, not an outlier. Graph-structural splits were
tried on the graph's TOPOLOGY (r300 bridges 0.3%, r301 connectivity inverts, r827 lambda>=2) but never
on the POSITION of a partner's alignment along the locus's own axis, which is what this measures.

    For a cut c:  L = partners aligning left of c, R = partners aligning right of c
    score(c) = |L & R| / |L | R|          (Jaccard; low = the two halves have different relatives)
    discontinuity(locus) = min over c of score(c)

⚠ The naive form is DEGENERATE: with a weak floor the minimum is always found at an extreme cut where
one side has 2 partners. Both sides must carry real, EXCLUSIVE evidence, and the comparison must be
LENGTH-MATCHED -- a long locus has more partners and more chances to find a low-J cut, so an unmatched
fused-vs-clean comparison measures length, not fusion.

Usage:
  partner_discontinuity.py --paf dn16.paf --gff chr16.genes.gff [--locus X]..."""
MIN_EXCL = 5          # each side needs this many partners the other side does NOT have

CUT_LO, CUT_HI = 0.15, 0.85   # ignore cuts near the ends

NBINS = 100

MIN_FRAC = 0.05       # a partner counts for a side if it covers this much of that side

def load_paf(paf):
    hits = collections.defaultdict(list)
    qlen = {}
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 11:
            continue
        q, t = f[0], f[5]
        qlen[q] = int(f[1]); qlen[t] = int(f[6])
        if q == t:
            continue
        hits[q].append((int(f[2]), int(f[3]), t))
        hits[t].append((int(f[7]), int(f[8]), q))
    return hits, qlen

def genes_per_locus(gff, loci):
    """locus -> [gene names whose OWN span is >=50% inside the locus]."""
    spans = []
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] not in ('gene', 'pseudogene'):
            continue
        m = re.search(r'Name=([^;]+)', f[8])
        if m:
            spans.append((int(f[3]), int(f[4]), m.group(1)))
    spans.sort()
    out = {}
    for name in loci:
        m = re.match(r'^\S+:(\d+)-(\d+)$', name)
        if not m:
            continue
        s, e = int(m.group(1)), int(m.group(2))
        inside = []
        for (gs, ge, g) in spans:
            if ge < s:
                continue
            if gs > e:
                break
            ov = min(e, ge) - max(s, gs)
            if ov > 0 and ov >= 0.5 * (ge - gs + 1):
                inside.append(g)
        out[name] = inside
    return out

def discontinuity(records, n):
    if not records or n <= 0:
        return None
    if len({p for (_, _, p) in records}) < 2 * MIN_EXCL:
        return None
    best = None
    for b in range(1, NBINS):
        frac = b / NBINS
        if frac < CUT_LO or frac > CUT_HI:
            continue
        c = n * frac
        lenL, lenR = c, n - c
        L, R = set(), set()
        for (qs, qe, p) in records:
            if min(qe, c) - qs >= MIN_FRAC * lenL:
                L.add(p)
            if qe - max(qs, c) >= MIN_FRAC * lenR:
                R.add(p)
        if len(L - R) < MIN_EXCL or len(R - L) < MIN_EXCL:
            continue
        j = len(L & R) / len(L | R)
        if best is None or j < best[0]:
            best = (j, frac, len(L - R), len(R - L), len(L & R))
    return best

def band(n):
    for hi, lbl in ((20000, '<20kb'), (50000, '20-50kb'), (100000, '50-100kb')):
        if n < hi:
            return lbl
    return '>=100kb'

def main_partner_discontinuity():
    ap = argparse.ArgumentParser()
    ap.add_argument('--paf', required=True)
    ap.add_argument('--gff', required=True)
    ap.add_argument('--locus', action='append', default=[])
    a = ap.parse_args()

    hits, qlen = load_paf(a.paf)
    gpl = genes_per_locus(a.gff, list(qlen))

    rows = []
    for name, n in qlen.items():
        d = discontinuity(hits.get(name, []), n)
        if d is None:
            continue
        rows.append((name, n, len(gpl.get(name, [])), d))

    print(f"scoreable loci (>= {2*MIN_EXCL} partners, both sides >= {MIN_EXCL} exclusive): "
          f"{len(rows)} of {len(qlen)}\n")

    print("=== named loci of interest ===")
    print(f"{'locus':32} {'len':>8} {'genes':>6} {'minJ':>7} {'cut@':>6} {'exL':>5} {'exR':>5} {'shared':>7}")
    for name in a.locus:
        hit = [r for r in rows if r[0] == name]
        if not hit:
            print(f"{name:32} {qlen.get(name,0):>8} {'-':>6} {'not scoreable':>7}")
            continue
        _, n, ng, (j, frac, exL, exR, sh) = hit[0]
        print(f"{name:32} {n:>8} {ng:>6} {j:>7.3f} {frac:>6.2f} {exL:>5} {exR:>5} {sh:>7}")

    print("\n=== minJ by LENGTH BAND x GENES-INSIDE (the confound-controlled view) ===")
    print(f"{'band':>10} {'genes':>7} {'n':>5} {'median minJ':>12} {'mean':>7}")
    grp = collections.defaultdict(list)
    for (_, n, ng, d) in rows:
        cls = '1 gene' if ng == 1 else ('>=2 genes' if ng >= 2 else '0 genes')
        grp[(band(n), cls)].append(d[0])
    for b in ('<20kb', '20-50kb', '50-100kb', '>=100kb'):
        for cls in ('1 gene', '>=2 genes'):
            v = grp.get((b, cls), [])
            if len(v) >= 3:
                print(f"{b:>10} {cls:>7} {len(v):>5} {statistics.median(v):>12.3f} "
                      f"{statistics.mean(v):>7.3f}")

# ================================================================================================
# soto-fidelity  (was bench/soto_family_locus_fidelity.py)  renamed: _open->_open_soto_family_locus_fidelity, load_genes->load_genes_soto_family_locus_fidelity, load_loci->load_loci_soto_family_locus_fidelity, main->main_soto_family_locus_fidelity
# ================================================================================================
"""How faithfully does each SOTO family's genes get reproduced by our de novo loci?

§6w3 measured the 5'/3' boundary asymmetry over all single-gene loci on chr16 and found the error is
DISPERSION not bias (median |d5| 229 bp vs |d3| 22 bp). This asks the next question the user posed:
run that per SOTO FAMILY, genome-wide, and rank families by how badly their loci approximate the truth,
so the worst ones can be looked at individually.

Metric-trap guards that are NOT optional here:
  * strand-aware 5'/3' (a minus-strand gene's 5' is its END coordinate);
  * ONE-TO-ONE gene->locus matching, greedy by exonic overlap; a locus already taken by a better gene
    cannot be reused, and the loser counts as a MISS, never as a second match (feedback_metric_traps);
  * width ratios are reported as an IN-BAND FRACTION, never as a median ratio (same file);
  * a locus covering >1 Soto/annotated gene is FLAGGED and reported separately -- a readthrough-fused
    locus's boundary error is not a boundary error, it is over-merge (§6v8), and pooling the two
    reproduces exactly the confound §6w3's population filter was built to avoid.

Usage:
  soto_family_locus_fidelity.py --gtf ours_genome.gtf --gff RefSeq.gff.gz --soto soto_famCN_S1C.tsv       [--out prefix]"""
TX = {'mRNA', 'transcript', 'ncRNA', 'lnc_RNA', 'lncRNA', 'pseudogenic_transcript', 'primary_transcript',
      'tRNA', 'rRNA', 'snRNA', 'snoRNA', 'miRNA', 'misc_RNA', 'V_gene_segment', 'C_gene_segment',
      'J_gene_segment', 'ncRNA_gene'}

def _open_soto_family_locus_fidelity(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def merge(iv):
    if not iv:
        return []
    iv = sorted(iv)
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]

def ov_bp(a, b):
    i = j = t = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi > lo:
            t += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return t

def soto_families(path):
    fam = {}
    for r in csv.DictReader(open(path), delimiter='\t'):
        f = (r.get('Family ID') or '').strip()
        g = (r.get('Gene Name') or '').strip()
        if f and f != 'N/A' and g:
            fam.setdefault(g, f)
    return fam

def load_genes_soto_family_locus_fidelity(gff, wanted):
    """gene name -> (chrom, start, end, strand, merged exons) for every gene, two-pass (RefSeq's own
    file does not guarantee parent-before-child order -- register 965)."""
    lines = _open_soto_family_locus_fidelity(gff).read().splitlines()
    gene, gid = {}, {}
    for ln in lines:
        if not ln or ln[0] == '#':
            continue
        f = ln.split('\t')
        if len(f) < 9 or f[2] not in ('gene', 'pseudogene', 'ncRNA_gene'):
            continue
        n = re.search(r'Name=([^;]+)', f[8])
        i = re.search(r'ID=([^;]+)', f[8])
        if n and i:
            gene[n.group(1)] = [f[0], int(f[3]), int(f[4]), f[6]]
            gid[i.group(1)] = n.group(1)
    t2g = {}
    for ln in lines:
        if not ln or ln[0] == '#':
            continue
        f = ln.split('\t')
        if len(f) < 9 or f[2] not in TX:
            continue
        i = re.search(r'ID=([^;]+)', f[8])
        p = re.search(r'Parent=([^;,]+)', f[8])
        if i and p and p.group(1) in gid:
            t2g[i.group(1)] = gid[p.group(1)]
    ex = collections.defaultdict(list)
    for ln in lines:
        if not ln or ln[0] == '#':
            continue
        f = ln.split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        p = re.search(r'Parent=([^;,]+)', f[8])
        if not p:
            continue
        par = p.group(1)
        g = t2g.get(par) or gid.get(par)
        if g:
            ex[g].append((int(f[3]), int(f[4])))
    out = {}
    for g, v in gene.items():
        e = merge(ex.get(g, []))
        out[g] = (v[0], v[1], v[2], v[3], e or [(v[1], v[2])])
    return out

def load_loci_soto_family_locus_fidelity(gtf):
    tx, ex = {}, collections.defaultdict(list)
    for ln in open(gtf):
        if not ln or ln[0] == '#':
            continue
        f = ln.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m:
            continue
        t = m.group(1)
        if f[2] == 'transcript':
            tx[t] = (f[0], int(f[3]), int(f[4]), f[6])
        elif f[2] == 'exon':
            ex[t].append((int(f[3]), int(f[4])))
    return {t: (v[0], v[1], v[2], v[3], merge(ex.get(t, []))) for t, v in tx.items()}

def main_soto_family_locus_fidelity():
    ap = argparse.ArgumentParser()
    for a in ('--gtf', '--gff', '--soto'):
        ap.add_argument(a, required=True)
    ap.add_argument('--out', default='soto_fidelity')
    a = ap.parse_args()

    fam = soto_families(a.soto)
    genes = load_genes_soto_family_locus_fidelity(a.gff, set(fam))
    loci = load_loci_soto_family_locus_fidelity(a.gtf)

    # index loci by chrom for overlap search
    bychrom = collections.defaultdict(list)
    for lid, (c, s, e, st, ex) in loci.items():
        bychrom[c].append((s, e, lid))
    for c in bychrom:
        bychrom[c].sort()

    # how many annotated genes does each locus substantially cover? (fusion flag)
    gene_iv = collections.defaultdict(list)
    for g, (c, s, e, st, ex) in genes.items():
        gene_iv[c].append((s, e, g))
    for c in gene_iv:
        gene_iv[c].sort()
    ngenes_in = {}
    for lid, (c, s, e, st, ex) in loci.items():
        n = 0
        for (gs, ge, g) in gene_iv.get(c, []):
            if ge < s:
                continue
            if gs > e:
                break
            o = min(e, ge) - max(s, gs)
            if o > 0 and o >= 0.5 * (ge - gs + 1):
                n += 1
        ngenes_in[lid] = n

    # candidate (gene, locus) pairs by exonic overlap, then ONE-TO-ONE greedy
    cands = []
    for g, f in fam.items():
        gi = genes.get(g)
        if not gi:
            continue
        c, gs, ge, st, gex = gi
        for (ls, le, lid) in bychrom.get(c, []):
            if le < gs:
                continue
            if ls > ge:
                break
            o = ov_bp(gex, loci[lid][4])
            if o > 0:
                cands.append((o, g, lid))
    cands.sort(reverse=True)
    used_l, used_g, match = set(), set(), {}
    for o, g, lid in cands:
        if g in used_g or lid in used_l:
            continue
        used_g.add(g); used_l.add(lid); match[g] = (lid, o)

    rows = []
    for g, f in fam.items():
        gi = genes.get(g)
        if not gi:
            continue
        c, gs, ge, st, gex = gi
        glen = sum(e - s for s, e in gex)
        if g not in match:
            rows.append(dict(fam=f, gene=g, matched=False))
            continue
        lid, o = match[g]
        lc, ls, le, lst, lex = loci[lid]
        llen = sum(e - s for e, s in [(e, s) for s, e in lex])
        llen = sum(e - s for s, e in lex)
        d5 = (ls - gs) if st == '+' else (ge - le)
        d3 = (ge - le) if st == '+' else (ls - gs)
        jac = o / max(glen + llen - o, 1)
        rows.append(dict(fam=f, gene=g, matched=True, locus=lid, d5=d5, d3=d3,
                         jac=jac, ratio=llen / max(glen, 1), fused=ngenes_in[lid] >= 2))

    matched = [r for r in rows if r['matched']]
    clean = [r for r in matched if not r['fused']]
    fused = [r for r in matched if r['fused']]
    print(f"Soto genes with a family: {len(rows)}   matched to a de novo locus: {len(matched)} "
          f"({100*len(matched)/len(rows):.1f}%)")
    print(f"  of matched: {len(clean)} clean single-gene loci, {len(fused)} on a locus covering >=2 genes "
          f"({100*len(fused)/len(matched):.1f}% FUSED -- reported separately, never pooled)\n")

    def band(v, t):
        return 100 * sum(1 for x in v if abs(x) <= t) / len(v) if v else 0

    for lbl, grp in (("CLEAN (1 gene/locus)", clean), ("FUSED (>=2 genes/locus)", fused)):
        if not grp:
            continue
        d5 = [r['d5'] for r in grp]; d3 = [r['d3'] for r in grp]; jac = [r['jac'] for r in grp]
        print(f"{lbl}  n={len(grp)}")
        print(f"   median |d5| {statistics.median([abs(x) for x in d5]):>7.0f} bp   "
              f"median |d3| {statistics.median([abs(x) for x in d3]):>6.0f} bp   "
              f"median exonic Jaccard {statistics.median(jac):.3f}")
        print(f"   within +-100bp: 5' {band(d5,100):>5.1f}%   3' {band(d3,100):>5.1f}%   "
              f"| width ratio in [0.8,1.25]: "
              f"{100*sum(1 for r in grp if 0.8<=r['ratio']<=1.25)/len(grp):.1f}%\n")

    # per-family aggregate, clean genes only
    byfam = collections.defaultdict(list)
    for r in clean:
        byfam[r['fam']].append(r)
    agg = []
    for f, rs in byfam.items():
        if len(rs) < 2:
            continue
        agg.append((statistics.median([x['jac'] for x in rs]), len(rs), f,
                    statistics.median([abs(x['d5']) for x in rs]),
                    statistics.median([abs(x['d3']) for x in rs]),
                    statistics.median([x['ratio'] for x in rs])))
    agg.sort()
    print(f"=== WORST families by median exonic Jaccard (>=2 clean genes matched; n={len(agg)}) ===")
    print(f"{'family':>10} {'n':>3} {'medJac':>7} {'med|d5|':>8} {'med|d3|':>8} {'ratio':>7}  members")
    for jac, n, f, m5, m3, ratio in agg[:15]:
        mem = ','.join(sorted(x['gene'] for x in byfam[f])[:4])
        print(f"{f:>10} {n:>3} {jac:>7.3f} {m5:>8.0f} {m3:>8.0f} {ratio:>7.2f}  {mem}")
    print(f"\n=== BEST families (same cohort) ===")
    for jac, n, f, m5, m3, ratio in agg[-5:]:
        mem = ','.join(sorted(x['gene'] for x in byfam[f])[:4])
        print(f"{f:>10} {n:>3} {jac:>7.3f} {m5:>8.0f} {m3:>8.0f} {ratio:>7.2f}  {mem}")

    with open(f'{a.out}.per_gene.tsv', 'w') as fh:
        fh.write("family\tgene\tmatched\tlocus\td5\td3\texonic_jaccard\twidth_ratio\tfused\n")
        for r in sorted(rows, key=lambda x: (x['fam'], x['gene'])):
            if r['matched']:
                fh.write(f"{r['fam']}\t{r['gene']}\t1\t{r['locus']}\t{r['d5']}\t{r['d3']}"
                         f"\t{r['jac']:.4f}\t{r['ratio']:.4f}\t{int(r['fused'])}\n")
            else:
                fh.write(f"{r['fam']}\t{r['gene']}\t0\t\t\t\t\t\t\n")
    print(f"\nper-gene table -> {a.out}.per_gene.tsv")


SUBCOMMANDS = {'colinearity': main_colinearity_conjunct, 'five-prime': main_five_prime_deficit, 'partner-discontinuity': main_partner_discontinuity, 'soto-fidelity': main_soto_family_locus_fidelity}


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in SUBCOMMANDS:
        sys.exit('usage: locus_probes.py {' + '|'.join(SUBCOMMANDS) + '} [args]')
    SUBCOMMANDS[sys.argv.pop(1)]()


if __name__ == '__main__':
    main()
