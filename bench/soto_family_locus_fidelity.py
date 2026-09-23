#!/usr/bin/env python3
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
  soto_family_locus_fidelity.py --gtf ours_genome.gtf --gff RefSeq.gff.gz --soto soto_famCN_S1C.tsv \
      [--out prefix]
"""
import argparse
import collections
import csv
import gzip
import re
import statistics

TX = {'mRNA', 'transcript', 'ncRNA', 'lnc_RNA', 'lncRNA', 'pseudogenic_transcript', 'primary_transcript',
      'tRNA', 'rRNA', 'snRNA', 'snoRNA', 'miRNA', 'misc_RNA', 'V_gene_segment', 'C_gene_segment',
      'J_gene_segment', 'ncRNA_gene'}


def _open(p):
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


def load_genes(gff, wanted):
    """gene name -> (chrom, start, end, strand, merged exons) for every gene, two-pass (RefSeq's own
    file does not guarantee parent-before-child order -- register 965)."""
    lines = _open(gff).read().splitlines()
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


def load_loci(gtf):
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


def main():
    ap = argparse.ArgumentParser()
    for a in ('--gtf', '--gff', '--soto'):
        ap.add_argument(a, required=True)
    ap.add_argument('--out', default='soto_fidelity')
    a = ap.parse_args()

    fam = soto_families(a.soto)
    genes = load_genes(a.gff, set(fam))
    loci = load_loci(a.gtf)

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


if __name__ == '__main__':
    main()
