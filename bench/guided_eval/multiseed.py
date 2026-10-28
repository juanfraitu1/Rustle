"""PREREG GUIDED-S9 (docs/LEDGER_O1.md): many templates per family vs one.

  multiseed.py prep OUTDIR S8DIR1 S8DIR2 FASTA        templates = the S8 seed-1 and seed-2 genes of every family: gene body (aligned
                                                      already, names 'fam|gene') + the CDS of EVERY transcript ('fam|gene|tx' in cds_all.fa)
  multiseed.py score OUTDIR ARM PAF...                ARM M = all templates; U2 = bodies + each gene's longest CDS only
Body PAFs are named A*.paf (CDS-base mapping as in guided3). Both seeded genes are excluded from hits and denominators.
"""
import collections
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cdsseed as C  # noqa: E402
import famtruth as F  # noqa: E402
import guided3 as g3  # noqa: E402


def prep(outdir, d1, d2, fasta):
    T1, T2 = json.load(open(f'{d1}/truth.seed1.json')), json.load(open(f'{d2}/truth.seed2.json'))
    G = g3.load_genes(T1['species'])
    seeded = collections.defaultdict(dict)  # fam -> gene -> body meta
    for T in (T1, T2):
        for q, m in T['seeds'].items():
            seeded[m['fam']][m['gene']] = m
    genes = {g for v in seeded.values() for g in v}
    cds = C.cat_cds(genes) if T1['species'] == 'human' else g3.ggo_cds(genes)
    C.FA = fasta
    want = {}
    for k, gs in seeded.items():
        for g in gs:
            for tx, segs in cds.get(g, {}).items():
                want[f'{k}|{g}|{tx}'] = sorted(segs, key=lambda x: x[1])
    plus = [f'{c}:{s + 1}-{e}' for segs in want.values() for c, s, e, st in segs if st == '+']
    minus = [f'{c}:{s + 1}-{e}' for segs in want.values() for c, s, e, st in segs if st == '-']
    sp, sm = C.faidx(plus, f'{outdir}/plus'), C.faidx(minus, f'{outdir}/minus', rc=True)
    templates = {}
    seen = set()
    with open(f'{outdir}/cds_all.fa', 'w') as w:
        for name, segs in want.items():
            st = segs[0][3]
            seq = ''.join(sp[f'{c}:{s + 1}-{e}'] for c, s, e, _ in segs) if st == '+' else \
                ''.join(sm[f'{c}:{s + 1}-{e}'] for c, s, e, _ in reversed(segs))
            k, g, tx = name.split('|')
            if (g, seq) in seen:  # identical CDS from another transcript: one template is enough
                continue
            seen.add((g, seq))
            w.write(f'>{name}\n{seq}\n')
            templates[name] = dict(fam=k, gene=g, cds_len=len(seq), body_cds=None)
    for k, gs in seeded.items():
        for g, m in gs.items():
            templates[f'{k}|{g}'] = dict(fam=k, gene=g, cds_len=m['cds_len'], body_cds=m['body_cds'])
    for name, t in templates.items():
        t['longest'] = t['cds_len'] == max(x['cds_len'] for x in templates.values() if x['gene'] == t['gene'])
    json.dump(dict(species=T1['species'], families=T1['families'], pseudo=T2['pseudo'], extra=T2['extra'], judge=T2['judge'],
                   seeded={k: sorted(v) for k, v in seeded.items()}, templates=templates), open(f'{outdir}/truth.json', 'w'))
    n_cds = sum(1 for n in templates if n.count('|') == 2)
    print(f'families {len(seeded)}, seeded genes {len(genes)}, CDS templates {n_cds} (longest {sum(t["longest"] for n, t in templates.items() if n.count("|") == 2)})')


def score(outdir, arm, pafs):
    T = json.load(open(f'{outdir}/truth.json'))
    G = g3.load_genes(T['species'])
    idx = F.exon_index(G)
    fams = {k: set(v) for k, v in T['families'].items()}
    pseudo = {k: set(v) for k, v in T['pseudo'].items()}
    extra = {k: set(v) for k, v in T['extra'].items()}
    judge = set(T['judge'])
    tm = T['templates']
    use = {n for n, t in tm.items() if arm == 'M' or t['longest'] or n.count('|') == 1}
    body_cds = {n: [tuple(x) for x in t['body_cds']] for n, t in tm.items() if t['body_cds']}
    groups = collections.defaultdict(list)
    for q, t, st, bl, qi, ti, cdsb, ident in g3.records(pafs, body_cds):
        if q not in use or not bl or cdsb <= 0:
            continue
        g = C.gene_of(idx, G, t, st, bl)
        if g and g not in T['seeded'][tm[q]['fam']]:
            groups[(q, g, st)].append((qi, ti, cdsb, ident))
    found = collections.defaultdict(set)
    for (q, g, st), rs in groups.items():
        bases, ident = g3.chain(rs, st)
        if bases / tm[q]['cds_len'] >= g3.COV and ident >= g3.IDENT:
            found[tm[q]['fam']].add(g)
    c = collections.Counter()
    mm = mp = nf = nfp = 0.0
    per, named = {}, {}
    for k, sg in T['seeded'].items():
        others = fams[k] - set(sg)
        ps = pseudo.get(k, set())
        fd = found.get(k, set())
        c['judged'] += len(fd & judge); c['ok'] += len(fd & (others | ps | extra.get(k, set()))); c['unjudged'] += len(fd - judge)
        c['members'] += len(others); c['mf'] += len(fd & others); c['pseudo'] += len(ps); c['pf'] += len(fd & ps)
        if others:
            nf += 1; mm += len(fd & others) / len(others)
        if ps:
            nfp += 1; mp += len(fd & ps) / len(ps)
        r = dict(seeds=sorted(G[g]['name'] for g in sg), n_others=len(others), members=sorted(G[g]['name'] for g in fd & others),
                 pseudo=sorted(G[g]['name'] for g in fd & ps), n_pseudo=len(ps),
                 wrong=sorted(G[g]['name'] for g in (fd & judge) - others - ps - extra.get(k, set())))
        per[k] = r
        for lab, pat in g3.NAMED.items():
            if any(re.match(pat, G[g]['name']) for g in sg):
                named[lab] = r
    res = dict(macro_member=mm / max(nf, 1), macro_pseudo=mp / max(nfp, 1), precision=c['ok'] / max(c['judged'], 1),
               precision_strict=c['ok'] / max(c['judged'] + c['unjudged'], 1), families=nf, families_with_pseudo=nfp, **c)
    print(f"{T['species']} {arm}: MACRO member {res['macro_member']:.4f} pseudo {res['macro_pseudo']:.4f} | precision {res['precision']:.4f} "
          f"({c['ok']}/{c['judged']}) strict {res['precision_strict']:.4f} (+{c['unjudged']}) | micro member {c['mf']}/{c['members']} "
          f"pseudo {c['pf']}/{c['pseudo']} | {int(nf)} fams with others")
    for lab in sorted(named):
        r = named[lab]
        print(f"  {lab}: seeds {','.join(r['seeds'])} members {len(r['members'])}/{r['n_others']} pseudo {len(r['pseudo'])}/{r['n_pseudo']} wrong {len(r['wrong'])}")
    json.dump(dict(totals=res, per_family=per, named=named), open(f'{outdir}/score.{arm}.json', 'w'), indent=0)


if __name__ == '__main__':
    if sys.argv[1] == 'prep':
        os.makedirs(sys.argv[2], exist_ok=True)
        prep(*sys.argv[2:6])
    else:
        score(sys.argv[2], sys.argv[3], sys.argv[4:])
