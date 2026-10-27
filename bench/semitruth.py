"""SEMITRUTH (docs/LEDGER_O1.md): human synteny truth for SEMISEC's unannotated silent candidates.

  semitruth.py GENOME_FA GGO_GFF CAT_GENES_TSV RMSK_OUT_GZ MMI WORK OUT_DIR...

Each OUT_DIR holds a `semisec` run already scored by `semisec.py report` (silent.json, cands.tsv). One minimap2 call for all of them.
"""
import bisect
import collections
import gzip
import json
import os
import random
import re
import subprocess
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import semisec as S  # noqa: E402

FLANK, MIN_PROJ = 2000, 0.5
CG = re.compile(r'(\d+)([MIDNSHP=X])')


def merged(ivs):
    out = []
    for s, e in sorted(ivs):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def covered(m, s, e):
    """bp of [s, e) inside the merged intervals m (sorted, disjoint)."""
    i = bisect.bisect_right([x[0] for x in m], s) - 1 if m else 0
    tot = 0
    for a, b in m[max(i, 0):]:
        if a >= e:
            break
        tot += max(0, min(b, e) - max(a, s))
    return tot


def rmsk(path):
    by = collections.defaultdict(list)
    for line in gzip.open(path, 'rt'):
        f = line.split()
        if len(f) < 15 or not f[0].isdigit():
            continue
        by[f[4]].append((int(f[5]) - 1, int(f[6])))
    return {c: merged(v) for c, v in by.items()}


def project(cig, strand, qs, qe, ts, a, b):
    """Target bases aligned to query bases in [a, b) (forward query coordinates): (min, max+1, n aligned)."""
    q = qs if strand == '+' else qe
    t, lo, hi, n = ts, None, None, 0
    for k, op in CG.findall(cig):
        k = int(k)
        if op in 'M=X':
            for j in range(k):
                qq = q + j if strand == '+' else q - 1 - j
                if a <= qq < b:
                    n += 1
                    lo = t + j if lo is None else min(lo, t + j)
                    hi = t + j + 1 if hi is None else max(hi, t + j + 1)
            q = q + k if strand == '+' else q - k
            t += k
        elif op == 'I':
            q = q + k if strand == '+' else q - k
        elif op in 'DN':
            t += k
    return lo, hi, n


def main(genome_fa, gff, cat, rm, mmi, work, *outs):
    os.makedirs(work, exist_ok=True)
    gi = S.biotypes(gff)
    fa = pysam.FastaFile(genome_fa)
    rnd = random.Random(1)
    items = []  # (sample, set, chrom, s, e, parent symbols)
    for out in outs:
        name = os.path.basename(out.rstrip('/'))
        sil = json.load(open(f'{out}/silent.json'))
        cands, _ = S.load(out)
        cidx = collections.defaultdict(list)
        for c in cands:
            cidx[c['chrom']].append((c['start'], c['end']))
        cidx = {k: S.Index(v) for k, v in cidx.items()}

        def parents(c):
            ch, iv = c['best_seed'].rsplit(':', 1)
            s0, e0 = map(int, iv.split('-'))
            return sorted({h[4][5:] for h in (gi[ch].over(s0, e0) if ch in gi else [])
                           if h[2] == 'protein_coding' and not h[4][5:].startswith('LOC')})

        for c in sil:
            st = 'U' if c['truth'] == 'none' else 'P' if c['truth'] == 'pseudogene' else None
            if st is None:
                continue
            items.append((name, st, c['chrom'], c['start'], c['end'], parents(c)))
            if st == 'U':
                L, clen = c['end'] - c['start'], fa.get_reference_length(c['chrom'])
                for _ in range(1000):
                    s = rnd.randrange(FLANK, max(FLANK + 1, clen - L - FLANK))
                    if not (gi.get(c['chrom']) and gi[c['chrom']].over(s, s + L)) and not (cidx.get(c['chrom']) and cidx[c['chrom']].over(s, s + L)):
                        items.append((name, 'N', c['chrom'], s, s + L, []))
                        break
    with open(f'{work}/q.fa', 'w') as w:
        for i, (_, _, ch, s, e, _) in enumerate(items):
            a = max(0, s - FLANK)
            b = min(fa.get_reference_length(ch), e + FLANK)
            w.write(f'>{i}|{s - a}|{e - a}\n{fa.fetch(ch, a, b)}\n')
    paf = f'{work}/q.paf'
    if not os.path.exists(paf):
        subprocess.run(f'minimap2 -x asm20 -c --secondary=no -t 5 {mmi} {work}/q.fa > {paf}.tmp 2> {work}/mm2.log && mv {paf}.tmp {paf}',
                       shell=True, check=True)
    hit = {}
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        tags = dict(x.split(':', 2)[::2] for x in f[12:])
        if tags.get('tp') != 'P' or int(f[11]) < 1:
            continue
        i, a, b = map(int, f[0].split('|'))
        lo, hi, n = project(tags['cg'], f[4], int(f[2]), int(f[3]), int(f[7]), a, b)
        if n >= MIN_PROJ * (b - a) and (i not in hit or n > hit[i][3]):
            hit[i] = (f[5], lo, hi, n)
    genes = collections.defaultdict(list)
    for line in open(cat):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'gene_id':
            continue
        genes[f[4]].append((int(f[5]), int(f[6]), f[2], f[1]))
    genes = {c: S.Index(v) for c, v in genes.items()}
    rep = rmsk(rm)
    hgnc = re.compile(r'^[A-Z][A-Z0-9-]*P\d+$')
    tab = collections.defaultdict(collections.Counter)
    rows = []
    for i, (smp, st, ch, s, e, par) in enumerate(items):
        rc = covered(rep.get(ch, []), s, e) / max(e - s, 1)
        if i not in hit:
            call, names = 'UNPLACED', []
        else:
            hc, lo, hi, _ = hit[i]
            hs = genes[hc].over(lo, hi) if hc in genes else []
            ps = [h for h in hs if 'pseudogene' in h[2]]
            call = 'PSEUDO' if ps else 'GENE' if hs else 'NONE'
            names = [h[3] for h in ps]
        k = tab[(smp, st)]
        k['n'] += 1
        k[call] += 1
        k['rep_bp'] += rc
        hn = [x for x in names if hgnc.match(x)]
        if call == 'PSEUDO' and hn and par:
            k['parent_den'] += 1
            k['parent_ok'] += any(re.match(rf'^{re.escape(p)}P\d+$', x) for x in hn for p in par)
        rows.append((smp, st, ch, s, e, call, ','.join(names), ','.join(par), f'{rc:.3f}', *(hit[i][:3] if i in hit else ('', '', ''))))
    with open(f'{work}/items.tsv', 'w') as w:
        w.write('sample\tset\tchrom\tstart\tend\tcall\thuman_pseudogenes\tparents\trepeat_frac\thchrom\thstart\thend\n')
        for r in rows:
            w.write('\t'.join(map(str, r)) + '\n')
    lines = ['| sample | set | n | PSEUDO | GENE | NONE | UNPLACED | repeat frac | parent match |', '|---|---|---|---|---|---|---|---|---|']
    for (smp, st), k in sorted(tab.items()):
        n = k['n']
        lines.append(f"| {smp} | {st} | {n} | " + ' | '.join(f"{k[x]} ({k[x] / n:.3f})" for x in ('PSEUDO', 'GENE', 'NONE', 'UNPLACED'))
                     + f" | {k['rep_bp'] / n:.3f} | {k['parent_ok']}/{k['parent_den']} |")
    for smp in sorted({x[0] for x in tab}):
        p, u, nn = (tab[(smp, x)] for x in 'PUN')
        rp, ru, rn = (t['PSEUDO'] / max(t['n'], 1) for t in (p, u, nn))
        t1, t2 = rp >= 0.50, ru >= 2 * rn
        lines.append(f"{smp}: T1 P PSEUDO {rp:.3f} >= 0.50 {'PASS' if t1 else 'FAIL'}; T2 U {ru:.3f} vs 2 x N {2 * rn:.3f} "
                     f"{'PASS' if t2 else 'FAIL'}; VERDICT {'PASS' if t1 and t2 else 'FAIL'}; est. real U copies "
                     f"{u['PSEUDO'] / max(rp, 1e-9):.0f} of {u['n']}")
    open(f'{work}/report.md', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main(*sys.argv[1:])
