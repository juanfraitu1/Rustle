#!/usr/bin/env python3
"""Chapter figures (S6): sweep status/agreement paired across unit catalogs, the ZNF569-like unit-extent
schematic, and the NPIP valid-anchor outcome table.  usage: chapter_figures.py <out_dir> <sweep_v2_dir>
<sweep_v3_dir> <sweep_v3_nofilter_dir>  (each sweep dir holds sweep_families.tsv from o2_sweep_analyse.py)."""
import sys, csv, json, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
out, d2, d3, d3n = sys.argv[1:5]
def fam(d): return {r['fam']: r for r in csv.DictReader(open(f'{d}/sweep_families.tsv'), delimiter='\t')}
F2, F3, F3n = fam(d2), fam(d3), fam(d3n)
fams = sorted(F3n, key=lambda f: -int(F3n[f]['reads_in_units']))
# --- Fig 1: status composition per family, three catalogs -------------------------------------------
fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
for ax, (F, title) in zip(axes, [(F2, 'units v2 (hull-clipped chain), PSV read-filter on'),
                                 (F3, 'units v3 (chain follows reads), read-filter on'),
                                 (F3n, 'units v3, read-filter off')]):
    tot = [int(F[f]['reads_in_units']) if f in F else 0 for f in fams]
    a = [int(F[f]['assigned']) / t if t else 0 for f, t in zip(fams, tot)]
    ti = [int(F[f]['tied']) / t if t else 0 for f, t in zip(fams, tot)]
    am = [int(F[f]['amb']) / t if t else 0 for f, t in zip(fams, tot)]
    x = range(len(fams))
    ax.bar(x, a, color='#2a6f97', label='assigned'); ax.bar(x, ti, bottom=a, color='#f6ae2d', label='tied')
    ax.bar(x, am, bottom=[i + j for i, j in zip(a, ti)], color='#c9c9c9', label='ambiguous')
    ax.set_ylim(0, 1); ax.set_ylabel('share of unit reads'); ax.set_title(title, fontsize=10, loc='left')
axes[0].legend(ncol=3, fontsize=8, loc='upper right')
axes[-1].set_xticks(range(len(fams))); axes[-1].set_xticklabels([f.replace('_0732', ' /') for f in fams], rotation=90, fontsize=7)
axes[-1].set_xlabel('sweep family (catalog cluster / contig), sorted by unit reads in v3')
fig.tight_layout(); fig.savefig(f'{out}/fig_sweep_status.png', dpi=160); plt.close(fig)
# --- Fig 2: MAPQ-60 placement agreement per family, v2 -> v3 ------------------------------------------
fig, ax = plt.subplots(figsize=(6.5, 7))
ys = []
for i, f in enumerate(fams):
    def ag(F):
        if f not in F or int(F[f]['q60_assigned']) == 0: return None
        return int(F[f]['q60_agree']) / int(F[f]['q60_assigned'])
    a2, a3 = ag(F2), ag(F3n)
    if a2 is not None and a3 is not None:
        ax.plot([a2, a3], [i, i], color='#bbbbbb', lw=1, zorder=1)
    if a2 is not None: ax.scatter(a2, i, color='#c9c9c9', edgecolor='k', zorder=2, s=28)
    if a3 is not None: ax.scatter(a3, i, color='#2a6f97', zorder=3, s=28)
ax.set_yticks(range(len(fams))); ax.set_yticklabels([f.replace('_0732', ' /') for f in fams], fontsize=7)
ax.set_xlabel('MAPQ-60 reads: O2 copy = aligner placement (fraction of assigned)'); ax.set_xlim(0, 1.02)
ax.scatter([], [], color='#c9c9c9', edgecolor='k', label='units v2, filter on'); ax.scatter([], [], color='#2a6f97', label='units v3, filter off')
ax.legend(loc='lower left', fontsize=8); ax.invert_yaxis(); fig.tight_layout(); fig.savefig(f'{out}/fig_sweep_agreement.png', dpi=160); plt.close(fig)
# --- Fig 3: the ZNF569-like locus: gene model, SD hull, v2 unit, v3 unit, read coverage -------------------
j = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'docs', 'figures', 'znf569_locus.json')))
lo, hi, bb = j['lo'], j['hi'], j['binbp']
fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(11, 4.6), sharex=True, gridspec_kw={'height_ratios': [2, 3]})
xs = [lo + i * bb for i in range(len(j['cov']))]
ax0.fill_between(xs, j['cov'], step='post', color='#2a6f97', alpha=0.8); ax0.set_ylabel('MAPQ-60 primary\nread blocks')
ax0.set_title(f"LOC101127631 (ZNF569-like), NC_073244.2 — {j['n_primary_q60']} MAPQ-60 reads", fontsize=10, loc='left')
def spans(s): return [tuple(map(int, p.split('-'))) for p in s.split(',')]
tracks = [('GFF exons (all isoforms)', j['gff_exons'], '#444444'), ('SD core hull (SEDEF, ≥ half the family)', spans(j['hull']), '#f6ae2d'),
          ('unit v2: chain clipped to the hull (809 bp)', spans(j['unit_v2']), '#b23a48'), ('unit v3: chain follows the reads', spans(j['unit_v3']), '#2a6f97')]
for k, (name, sp, col) in enumerate(tracks):
    y = len(tracks) - k
    for s, e in sp: ax1.add_patch(plt.Rectangle((s, y - 0.3), e - s, 0.6, color=col))
    if k >= 1: ax1.plot([min(s for s, _ in sp), max(e for _, e in sp)], [y, y], color=col, lw=0.8, zorder=0)
    ax1.text(lo + 200, y + 0.36, name, fontsize=8, va='bottom')
ax1.set_ylim(0.4, len(tracks) + 0.9); ax1.set_yticks([]); ax1.set_xlim(lo, hi); ax1.set_xlabel('genome position (bp)')
ax1.ticklabel_format(style='plain', axis='x'); fig.tight_layout(); fig.savefig(f'{out}/fig_znf569_locus.png', dpi=160); plt.close(fig)
# --- Fig 4: NPIP valid anchors (62), outcome per arm (o2scale/rescore_valid_truth.log, §6es) ----------------
arms = [('GFF-model cores, spliced', 5, 5, 38, 19), ('faithful units, spliced', 3, 3, 6, 53), ('units + genomic hulls (off)', 8, 3, 34, 20)]
fig, ax = plt.subplots(figsize=(7, 2.8))
for i, (n, asg, ok, tied, amb) in enumerate(arms):
    ax.barh(i, ok, color='#2a6f97'); ax.barh(i, asg - ok, left=ok, color='#b23a48'); ax.barh(i, tied, left=asg, color='#f6ae2d'); ax.barh(i, amb, left=asg + tied, color='#c9c9c9')
ax.set_yticks(range(len(arms))); ax.set_yticklabels([a[0] for a in arms], fontsize=9); ax.invert_yaxis(); ax.set_xlabel('62 junction-anchored NPIP reads (annotation gaps removed)')
for c, l in [('#2a6f97', 'assigned, agrees'), ('#b23a48', 'assigned, wrong'), ('#f6ae2d', 'tied'), ('#c9c9c9', 'ambiguous')]: ax.barh([], [], color=c, label=l)
ax.legend(ncol=4, fontsize=8, loc='lower right', bbox_to_anchor=(1, 1.02)); fig.tight_layout(); fig.savefig(f'{out}/fig_npip_anchors.png', dpi=160); plt.close(fig)
print('wrote 4 figures to', out)
