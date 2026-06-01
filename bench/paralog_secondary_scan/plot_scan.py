#!/usr/bin/env python3
"""Figure for the paralog secondary-dependence scan.
Panel A: funnel (expressed candidates -> secondary-dependent -> taxonomy).
Panel B: per-hit separation (G-anchored vs total chain reads, colored by taxonomy).
Panel C: per-read dNM (distinguishing-mismatch) distributions for 3 exemplars,
         showing the trimodal edit-distance signal (sibling | tie | own-copy).

Usage: plot_scan.py OUT.png
"""
import sys, re, json, bisect, collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def parse_cig(c): return re.findall(r'(\d+)([MIDNSH=X])', c)
def junctions(pos, cig):
    cur = pos; js = []
    for ln, op in parse_cig(cig):
        ln = int(ln)
        if op == 'N': js.append((cur - 1, cur + ln)); cur += ln
        elif op in 'MD=X': cur += ln
    return js
def alen(c): return sum(int(l) for l, o in parse_cig(c) if o in 'M=XD') or 1

hits = json.load(open('/tmp/hits_nm.json'))
genes = {}
for line in open('/tmp/gene_introns.tsv'):
    if line.startswith('gene_id'): continue
    gid, chrom, s, e, strand, ntx, chain = line.rstrip('\n').split('\t')
    ints = set()
    if chain:
        for pr in chain.split(','):
            d, a = pr.split(':'); ints.add((int(d), int(a)))
    genes[gid] = dict(chrom=chrom, start=int(s), end=int(e), strand=strand, introns=ints)

# per-read dNM for exemplars
reads = collections.defaultdict(list); pos_index = collections.defaultdict(list)
for line in open('/tmp/mm_placements.tsv'):
    f = line.rstrip('\n').split('\t')
    if len(f) < 6 or f[5] == '' or f[4] == '*': continue
    qn, flag, chrom, pos, cig, nm = f
    p = int(pos); reads[qn].append((chrom, p, cig, int(nm))); pos_index[chrom].append((p, qn))
for c in pos_index: pos_index[c].sort()

def dnm_list(gid):
    g = genes[gid]; K = max(2, min(8, (len(g['introns']) + 2) // 3))
    arr = pos_index.get(g['chrom'], [])
    lo = bisect.bisect_left(arr, (g['start'] - 10, '')); hi = bisect.bisect_right(arr, (g['end'] + 10, chr(0x10ffff)))
    qns = set(qn for _, qn in arr[lo:hi]); out = []
    for qn in qns:
        gnm = None; galen = None
        for (chrom, pos, cig, nm) in reads[qn]:
            if chrom != g['chrom'] or not (g['start'] - 10 <= pos <= g['end']): continue
            if sum(1 for j in junctions(pos, cig) if j in g['introns']) >= K:
                if gnm is None or nm < gnm: gnm = nm; galen = alen(cig)
        if gnm is None: continue
        best_other = None
        for (chrom, pos, cig, nm) in reads[qn]:
            if chrom == g['chrom'] and g['start'] - 10 <= pos <= g['end']: continue
            if alen(cig) < 0.7 * galen: continue
            if best_other is None or nm < best_other: best_other = nm
        if best_other is None: out.append(50)            # uniquely G (cap for display)
        else: out.append(max(-50, min(50, best_other - gnm)))
    return out

fig = plt.figure(figsize=(15, 4.6))

# Panel A: funnel
axA = fig.add_subplot(1, 3, 1)
tax = collections.Counter(h['nm_taxonomy'] for h in hits)
order = ['expressed_real_copy', 'pure_spillover', 'weak', 'ambiguous_dominated']
colors = {'expressed_real_copy': '#2c7fb8', 'pure_spillover': '#d95f02',
          'weak': '#999999', 'ambiguous_dominated': '#7570b3'}
vals = [tax.get(k, 0) for k in order]
axA.bar(['846\nexpressed\ncandidates'], [846], color='#cccccc', width=0.6)
axA.bar(['252\nsecondary-\ndependent'], [252], color='#bbbbbb', width=0.6)
bottom = 0
for k in order:
    axA.bar(['252 by\nedit-distance\ntaxonomy'], [tax.get(k, 0)], bottom=bottom,
            color=colors[k], width=0.6, label=f"{k} ({tax.get(k,0)})")
    bottom += tax.get(k, 0)
axA.set_ylabel('genes'); axA.set_title('A. Genome-wide funnel'); axA.legend(fontsize=7, loc='upper right')

# Panel B: separation
axB = fig.add_subplot(1, 3, 2)
for k in order:
    xs = [int(h['n_chain_reads']) for h in hits if h['nm_taxonomy'] == k]
    ys = [max(0.5, h['nm_anchor_G']) for h in hits if h['nm_taxonomy'] == k]
    axB.scatter(xs, ys, s=18, alpha=0.7, color=colors[k], label=k)
axB.set_xscale('log'); axB.set_yscale('log')
axB.set_xlabel('chain-supporting reads at locus'); axB.set_ylabel('reads decisively own (dNM>=2)')
axB.axhline(3, ls='--', c='k', lw=0.8); axB.set_title('B. Expression anchor vs coverage')
axB.legend(fontsize=7)

# Panel C: dNM exemplars
axC = fig.add_subplot(1, 3, 3)
exemplars = [('gene-LOC134756861', 'NBPF8-like (expressed)', '#2c7fb8'),
             ('gene-LOC129530216', 'DAZ3 (weak)', '#999999'),
             ('gene-LOC101136473', 'SORD-like (spillover)', '#d95f02')]
bins = np.arange(-30, 32, 2)
for gid, lab, col in exemplars:
    if gid in genes:
        d = dnm_list(gid)
        axC.hist(d, bins=bins, histtype='step', lw=1.8, color=col, label=f"{lab} n={len(d)}")
axC.axvspan(-2, 2, color='#eeeeee', zorder=0)
axC.axvline(0, c='k', lw=0.6)
axC.set_xlabel('dNM = NM(sibling) - NM(this gene)\n<-- read is sibling\'s | tie | read is this gene\'s -->')
axC.set_ylabel('reads'); axC.set_title('C. Per-read edit-distance signal'); axC.legend(fontsize=7)

plt.tight_layout()
plt.savefig(sys.argv[1], dpi=130)
print('wrote', sys.argv[1])
