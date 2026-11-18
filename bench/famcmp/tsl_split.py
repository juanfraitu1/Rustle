#!/usr/bin/env python3
"""PREREG r1788 (docs/PREREG_tsl_split_2026-11-18.md): recovery of reference intron chains by Ensembl transcript support level.

  tsl_split.py REF.gff3 OUT.json NAME=arm.gtf [NAME=arm.gtf ...]
REF = the CHM13 CAT/Liftoff annotation (transcript rows carry transcript_support_level, source_transcript, gene_id)."""
import collections, json, re, sys
LEVELS = ['1', '2', '3', '4', '5', 'NA', 'novel']


def chain(exons):
    e = sorted(exons); return tuple((a[1], b[0]) for a, b in zip(e, e[1:]))


assert chain([(10, 20), (1, 5)]) == ((5, 10),) and chain([(1, 5)]) == ()


def arm_chains(path):
    ex = collections.defaultdict(list)
    for l in open(path):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon':
            ex[(f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]), int(f[4])))
    return {(k[0], k[1], chain(v)) for k, v in ex.items() if len(v) > 1}


def main(ref, out, *arms):
    lvl, gene, ex = {}, {}, collections.defaultdict(list)
    for l in open(ref):
        f = l.split('\t')
        if len(f) < 9 or f[2] not in ('transcript', 'exon'): continue
        a = dict(x.split('=', 1) for x in f[8].strip().split(';') if '=' in x)
        if f[2] == 'transcript':
            t = a.get('transcript_support_level', 'NA')
            lvl[a['ID']] = 'novel' if a.get('source_transcript', 'N/A') == 'N/A' else (t if t in LEVELS else 'NA'); gene[a['ID']] = a['gene_id']
        else: ex[(f[0], f[6], a['Parent'])].append((int(f[3]), int(f[4])))
    R = {}                                                                 # chain -> (best level, gene)
    for (c, s, t), v in ex.items():
        if len(v) < 2: continue
        k = (c, s, chain(v))
        if k not in R or LEVELS.index(lvl[t]) < LEVELS.index(R[k][0]): R[k] = (lvl[t], gene[t])
    A = {n: arm_chains(p) for n, p in (x.split('=', 1) for x in arms)}
    expressed = {R[k][1] for a in A.values() for k in a if k in R}
    res = {'levels': LEVELS, 'ref': {}, 'ref_expressed': {}, 'arms': {}}
    for L in LEVELS:
        res['ref'][L] = sum(v[0] == L for v in R.values()); res['ref_expressed'][L] = sum(v[0] == L and v[1] in expressed for v in R.values())
    print('reference chains by level:', res['ref'], '\nin expressed genes:', res['ref_expressed'])
    for n, a in A.items():
        hit = collections.Counter(R[k][0] for k in a if k in R)
        res['arms'][n] = {'multi_exon_chains': len(a), 'matched': sum(hit.values()), 'recovered': {L: hit[L] for L in LEVELS}}
        print(f"\n{n}: {len(a)} multi-exon chains, {sum(hit.values())} match a reference chain\nlevel\tref\trecovered\tshare\tref_expressed\tshare_expressed\tshare_of_matches")
        for L in LEVELS:
            r, e = res['ref'][L], res['ref_expressed'][L]
            print(f"{L}\t{r}\t{hit[L]}\t{hit[L] / r:.4f}\t{e}\t{hit[L] / e if e else 0:.4f}\t{hit[L] / max(sum(hit.values()), 1):.4f}")
    json.dump(res, open(out, 'w'), indent=1)


if __name__ == '__main__': main(*sys.argv[1:])
