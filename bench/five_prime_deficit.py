#!/usr/bin/env python3
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
  five_prime_deficit.py --gtf dn16.gtf --gff chr16.genes.gff [--chrom chr16] [--label human-chr16]
"""
import argparse
import collections
import re
import statistics


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


def main():
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


if __name__ == '__main__':
    main()
