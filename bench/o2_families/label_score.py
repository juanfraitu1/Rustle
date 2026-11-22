#!/usr/bin/env python3
"""Score the o2_status labels of an extended assembly against planted truth, by unit (register 1911).

    label_score.py PLANT_DIR SCENARIO LABELLED.gtf [OTHER=GTF ...]

Needs PLANT_DIR/units.tsv (units.py). A transcript belongs to the planted copy it shares the most exon bases with (as
score.py). A unit's label = the best o2_status among the transcripts on its copies (expressed > shared >
homology_only); no transcript = absent. Prints, for the labelled assembly, expressed and silent units by label, and for
it and every OTHER assembly the units with a transcript (expressed, silent, total)."""
import sys, os, re, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score as S
RANK = {'expressed': 3, 'shared': 2, 'homology_only': 1}

def copy_labels(truth, tool, status):
    """-> {copy: best label or 'unlabelled'} for every copy with a transcript"""
    cex = collections.defaultdict(list); cloc = {}
    for t, (c, st, ex, g) in truth.items(): cex[g] += ex; cloc[g] = (c, st)
    by = collections.defaultdict(list)
    for g in cex: cex[g] = S.merge(cex[g]); by[cloc[g][0]].append(g)
    out = {}
    for t, (c, st, ex, _) in tool.items():
        best, bg = 0, None
        for g in by.get(c, ()):
            if st in '+-' and st != cloc[g][1]: continue
            if cex[g][-1][1] <= ex[0][0] or cex[g][0][0] >= ex[-1][1]: continue
            n = S.ov(ex, cex[g])
            if n >= 50 and (n > best or (n == best and g > (bg or ''))): best, bg = n, g
        if bg is None: continue
        l = status.get(t, 'unlabelled')
        if RANK.get(l, 0) >= RANK.get(out.get(bg, ''), -1): out[bg] = l
    return out

def main():
    d, sc, lab = sys.argv[1:4]; others = sys.argv[4:]
    truth, planted, fam = S.load_plan(d, sc)
    creads = collections.Counter()
    for t, n in planted.items(): creads[truth[t][3]] += n
    U = {l.split('\t')[1]: l.split('\t')[3] for l in list(open(f'{d}/units.tsv'))[1:]}
    units = collections.defaultdict(set)
    for g in creads: units[U[g]].add(g)
    uexp = {u for u, gs in units.items() if any(creads[g] > 0 for g in gs)}
    def unit_labels(path, with_status):
        st = {}
        if with_status:
            for l in open(path):
                f = l.split('\t')
                if len(f) > 8 and f[2] == 'transcript':
                    t = re.search(r'transcript_id "([^"]+)"', f[8]); s = re.search(r'o2_status "([^"]+)"', f[8])
                    if t and s: st[t.group(1)] = s.group(1)
        cl = copy_labels(truth, S.load_gtf(path), st); ul = {}
        for u, gs in units.items():
            ls = [cl[g] for g in gs if g in cl]
            if ls: ul[u] = max(ls, key=lambda x: RANK.get(x, 0))
        return ul
    ul = unit_labels(lab, True); c = collections.Counter()
    for u in units: c[('exp' if u in uexp else 'sil', ul.get(u, 'absent'))] += 1
    ne, ns = len(uexp), len(units) - len(uexp)
    print(f"scenario {sc}: expressed units {ne}: expressed {c[('exp','expressed')]}, shared {c[('exp','shared')]}, homology_only {c[('exp','homology_only')]}, absent {c[('exp','absent')]}")
    sw = c[('sil','expressed')] + c[('sil','shared')] + c[('sil','homology_only')]
    print(f"scenario {sc}: silent units {ns}: with a transcript {sw}: homology_only {c[('sil','homology_only')]} ({c[('sil','homology_only')] / sw:.3f} of them), shared {c[('sil','shared')]}, expressed {c[('sil','expressed')]}" if sw else f"scenario {sc}: silent units {ns}: none with a transcript")
    print('available units (expressed + silent = total):', end=' ')
    for x in [f'labelled={lab}'] + others:
        n, p = x.split('=', 1)
        if not os.path.exists(p): print(f'{n} MISSING;', end=' '); continue
        v = unit_labels(p, False); e = sum(1 for u in v if u in uexp); print(f'{n} {e} + {len(v) - e} = {len(v)};', end=' ')
    print()

if __name__ == '__main__':
    main()
