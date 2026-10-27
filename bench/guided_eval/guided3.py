"""PREREG GUIDED-S3 (docs/LEDGER_O1.md): Liftoff-chained spliced seeds for the guided family search, human dev, gorilla held-out.

  guided3.py prep SPECIES OUTDIR [--seed 1|2]          truth (strict HGNC most-specific group; Compara Simiiformes for group-less genes),
                                                        seeds (cds.fa, body.fa, body CDS intervals)
  guided3.py score SPECIES OUTDIR ARM PAF... [--seed]  chained per (seed, target gene, strand); found at coverage >= 0.5, identity >= 0.5

SPECIES = human (CHM13 v2.0, CAT/Liftoff v2.0) or gorilla (mGorGor1 v2.0, RefSeq; truth by gene symbol, LOC genes outside).
"""
import argparse
import collections
import csv
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
import cdsseed as C  # noqa: E402
import famtruth as F  # noqa: E402

FASTA = {'human': '/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa',
         'gorilla': '/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta'}
COV, IDENT = 0.5, 0.5
NAMED = {**F.MAINSTAYS, **F.YAGS}


def load_genes(species):
    human = F.load_cat()
    if species == 'human':
        F.GENES.clear()
        F.GENES.update(human)
        return human
    return F.load_ggo(human)


def hgnc_tables():
    rows = list(csv.DictReader(open(C.HGNC), delimiter='\t'))
    size = collections.Counter(g for r in rows for g in r['gene_group_id'].split('|') if g)
    by_ens, by_sym = {}, {}
    for r in rows:
        gs = [g for g in r['gene_group_id'].split('|') if g]
        if not gs:
            continue
        best = min(gs, key=lambda g: (size[g], int(g)))
        if r['ensembl_gene_id']:
            by_ens[r['ensembl_gene_id']] = best
        by_sym[r['symbol']] = best
    return by_ens, by_sym, size


def group_of(species, d, by_ens, by_sym, size):
    if species == 'gorilla':
        return by_sym.get(d['name'])
    gs = {by_ens[e] for e in d['ens'] if e in by_ens}
    return min(gs, key=lambda g: (size[g], int(g))) if gs else None


def truth(species, G):
    by_ens, by_sym, size = hgnc_tables()
    pc = {g for g, d in G.items() if d['biotype'] == 'protein_coding' and (species == 'human' or not d['name'].startswith('LOC'))}
    grp = {g: group_of(species, G[g], by_ens, by_sym, size) for g in pc}
    fams = collections.defaultdict(set)
    for g in pc:
        if grp[g]:
            fams['H' + grp[g]].add(g)
    pairs, universe = F.truth_pairs(G, F.PARALOGS, F.LEVELS['simiiformes'])
    free = {g for g in pc if not grp[g]}
    for c in F.components({p for p in pairs if p[0] in free and p[1] in free}):
        fams['C' + min(c)] = set(c)
    fams = {k: v for k, v in fams.items() if len(v) >= 2}
    fam_of_group = {k[1:]: k for k in fams if k[0] == 'H'}
    key_fam = collections.defaultdict(set)
    for k, c in fams.items():
        for g in c:
            for s in C.stem_keys(G[g]['name']):
                key_fam[s].add(k)
    pseudo = collections.defaultdict(set)
    for g, d in G.items():
        if 'pseudogene' not in d['biotype'] or d['name'].startswith('LOC'):
            continue
        cand = set()
        gg = group_of(species, d, by_ens, by_sym, size)
        if gg in fam_of_group:
            cand.add(fam_of_group[gg])
        m = re.match(r'^(.+?)P[0-9]*$', d['name'])
        if m:
            cand |= key_fam.get(m.group(1), set())
        if len(cand) == 1:
            pseudo[next(iter(cand))].add(g)
    judge = set(universe) | set(pc) | {g for v in pseudo.values() for g in v}
    return fams, pseudo, judge


def ggo_cds(wanted):
    """RefSeq gene id -> {mRNA: [(chrom, s0, e, strand)]}."""
    tx_gene = {}
    out = collections.defaultdict(lambda: collections.defaultdict(list))
    for line in open(F.GGO_GFF):
        if line[0] == '#':
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        if f[2] == 'mRNA':
            tx_gene[a['ID']] = a.get('Parent')
        elif f[2] == 'CDS':
            g = tx_gene.get(a.get('Parent'))
            if g in wanted:
                out[g][a['Parent']].append((f[0], int(f[3]) - 1, int(f[4]), f[6]))
    return out


def order_key(species, G, g):
    return (min(G[g]['ens']) if G[g]['ens'] else 'ZZZ', g) if species == 'human' else (G[g]['name'], g)


def prep(a):
    G = load_genes(a.species)
    fams, pseudo, judge = truth(a.species, G)
    pcs = {g for c in fams.values() for g in c}
    cds = C.cat_cds(pcs) if a.species == 'human' else ggo_cds(pcs)
    C.FA = FASTA[a.species]
    seeds = {}
    for k, c in sorted(fams.items()):
        have = [g for g in c if cds.get(g)]
        label = [g for g in have if any(re.match(p, G[g]['name']) for p in NAMED.values())]
        order = sorted(label, key=lambda g: order_key(a.species, G, g)) + \
            sorted((g for g in have if g not in label), key=lambda g: order_key(a.species, G, g))
        if len(order) >= a.seed:
            g = order[a.seed - 1]
            tx = max(cds[g], key=lambda t: (sum(e - s for _, s, e, _ in cds[g][t]), t))
            seeds[k] = (g, sorted(cds[g][tx], key=lambda x: x[1]))
    plus = [f'{c}:{s + 1}-{e}' for g, segs in seeds.values() for c, s, e, st in segs if st == '+']
    minus = [f'{c}:{s + 1}-{e}' for g, segs in seeds.values() for c, s, e, st in segs if st == '-']
    sp = C.faidx(plus, os.path.join(a.outdir, 'plus'))
    sm = C.faidx(minus, os.path.join(a.outdir, 'minus'), rc=True)
    meta, body = {}, []
    with open(os.path.join(a.outdir, 'cds.fa'), 'w') as fc:
        for k, (g, segs) in seeds.items():
            st = segs[0][3]
            seq = ''.join(sp[f'{c}:{s + 1}-{e}'] for c, s, e, _ in segs) if st == '+' else \
                ''.join(sm[f'{c}:{s + 1}-{e}'] for c, s, e, _ in reversed(segs))
            name = f'{k}|{g}'
            fc.write(f'>{name}\n{seq}\n')
            d = G[g]
            b0, b1 = min(s for s, _ in d['exons']), max(e for _, e in d['exons'])
            b0, b1 = min(b0, segs[0][1]), max(b1, segs[-1][2])
            # seed CDS intervals in body-local coordinates (transcript orientation)
            loc = sorted((s - b0, e - b0) if st == '+' else (b1 - e, b1 - s) for _, s, e, _ in segs)
            body.append((f"{d['chrom']}:{b0 + 1}-{b1}", name, st))
            meta[name] = dict(fam=k, gene=g, cds_len=len(seq), body_cds=loc)
    bp = C.faidx([r for r, _, st in body if st != '-'], os.path.join(a.outdir, 'bodyp'))
    bm = C.faidx([r for r, _, st in body if st == '-'], os.path.join(a.outdir, 'bodym'), rc=True)
    with open(os.path.join(a.outdir, 'body.fa'), 'w') as w:
        for r, name, st in body:
            w.write(f'>{name}\n{(bm if st == "-" else bp)[r]}\n')
    json.dump(dict(species=a.species, seed=a.seed, families={k: sorted(v) for k, v in fams.items()},
                   pseudo={k: sorted(v) for k, v in pseudo.items()}, judge=sorted(judge), seeds=meta),
              open(os.path.join(a.outdir, f'truth.seed{a.seed}.json'), 'w'))
    sz = sorted(map(len, fams.values()), reverse=True)
    print(f'{a.species}: families {len(fams)} (largest {sz[:5]}), seeded {len(seeds)}, pseudogenes {sum(map(len, pseudo.values()))} '
          f'in {len(pseudo)} families, named seeds {sum(1 for k, (g, _) in seeds.items() if any(re.match(p, G[g]["name"]) for p in NAMED.values()))}')


def overlap(iv, s, e):
    return sum(max(0, min(b, e) - max(a, s)) for a, b in iv)


def records(paths, body_cds):
    """Per PAF record: (seed, chrom, strand, genomic blocks, seed interval, target interval, seed-CDS bases aligned, identity)."""
    for p in paths:
        bc = body_cds if os.path.basename(p).startswith('A') else {}  # gene-body PAFs are named A*.paf
        for line in open(p):
            f = line.rstrip('\n').split('\t')
            q, qs, qe, st, t, ts = f[0], int(f[2]), int(f[3]), f[4], f[5], int(f[7])
            cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), '')
            bl, pos, cur, qo, cdsb = [], ts, None, 0, 0
            for n, op in C.CIG.findall(cg):
                n = int(n)
                if op in 'M=X':
                    qa = qs + qo if st == '+' else qe - qo - n
                    cdsb += overlap(bc[q], qa, qa + n) if q in bc else n
                    qo += n
                if op in 'M=XD':
                    cur = [pos, pos + n] if cur is None else [cur[0], pos + n]
                    pos += n
                elif op == 'N':
                    if cur:
                        bl.append(tuple(cur))
                    cur = None
                    pos += n
                elif op == 'I':
                    qo += n
            if cur:
                bl.append(tuple(cur))
            yield q, t, st, bl, (qs, qe), (ts, int(f[8])), cdsb, int(f[9]) / int(f[10])


def chain(rs, strand):
    rs = sorted(rs, key=lambda r: (r[0][0], r[0][1]))
    best, dp = (0, 0.0), []
    for j, r in enumerate(rs):
        cur = (r[2], r[2] * r[3])
        for i in range(j):
            p = rs[i]
            if p[0][1] > r[0][0]:
                continue
            if (strand == '+' and p[1][1] > r[1][0]) or (strand == '-' and r[1][1] > p[1][0]):
                continue
            if dp[i][0] + r[2] > cur[0]:
                cur = (dp[i][0] + r[2], dp[i][1] + r[2] * r[3])
        dp.append(cur)
        if cur[0] > best[0]:
            best = cur
    return best[0], best[1] / best[0] if best[0] else 0.0


def score(a):
    T = json.load(open(os.path.join(a.outdir, f'truth.seed{a.seed}.json')))
    G = load_genes(a.species)
    idx = F.exon_index(G)
    fams = {k: set(v) for k, v in T['families'].items()}
    pseudo = {k: set(v) for k, v in T['pseudo'].items()}
    judge = set(T['judge'])
    body_cds = {q: [tuple(x) for x in m['body_cds']] for q, m in T['seeds'].items()}
    groups = collections.defaultdict(list)
    for q, t, st, bl, qi, ti, cdsb, ident in records(a.paf, body_cds):
        if q not in T['seeds'] or not bl or cdsb <= 0:
            continue
        g = C.gene_of(idx, G, t, st, bl)
        if g and g != T['seeds'][q]['gene']:
            groups[(q, g, st)].append((qi, ti, cdsb, ident))
    found = collections.defaultdict(set)
    for (q, g, st), rs in groups.items():
        bases, ident = chain(rs, st)
        if bases / T['seeds'][q]['cds_len'] >= COV and ident >= IDENT:
            found[q].add(g)
    c = collections.Counter()
    mm = mp = nf = nfp = 0.0
    per, named = {}, {}
    for q, m in T['seeds'].items():
        others = fams[m['fam']] - {m['gene']}
        ps = pseudo.get(m['fam'], set())
        fd = found.get(q, set())
        c['members'] += len(others); c['mf'] += len(fd & others); c['pseudo'] += len(ps); c['pf'] += len(fd & ps)
        c['judged'] += len(fd & judge); c['ok'] += len(fd & (others | ps)); c['unjudged'] += len(fd - judge)
        nf += 1; mm += len(fd & others) / len(others)
        if ps:
            nfp += 1; mp += len(fd & ps) / len(ps)
        r = dict(seed=G[m['gene']]['name'], n=len(others) + 1, members=sorted(G[g]['name'] for g in fd & others),
                 pseudo=sorted(G[g]['name'] for g in fd & ps), n_pseudo=len(ps), wrong=sorted(G[g]['name'] for g in (fd & judge) - others - ps))
        per[q] = r
        for lab, pat in NAMED.items():
            if re.match(pat, G[m['gene']]['name']):
                named[lab] = r
    res = dict(macro_member=mm / max(nf, 1), macro_pseudo=mp / max(nfp, 1), precision=c['ok'] / max(c['judged'], 1),
               precision_strict=c['ok'] / max(c['judged'] + c['unjudged'], 1), micro_member=c['mf'] / max(c['members'], 1),
               micro_pseudo=c['pf'] / max(c['pseudo'], 1), families=nf, families_with_pseudo=nfp, **c)
    print(f"{a.species} {a.arm} seed{a.seed}: MACRO member {res['macro_member']:.4f} pseudo {res['macro_pseudo']:.4f} | precision "
          f"{res['precision']:.3f} ({c['ok']}/{c['judged']}), unjudged-as-wrong {res['precision_strict']:.3f} (+{c['unjudged']}) | micro "
          f"member {c['mf']}/{c['members']} pseudo {c['pf']}/{c['pseudo']} | {int(nf)} fams")
    for lab in sorted(named):
        r = named[lab]
        print(f"  {lab}: seed {r['seed']} members {len(r['members'])}/{r['n'] - 1} pseudo {len(r['pseudo'])}/{r['n_pseudo']} wrong {len(r['wrong'])}")
    json.dump(dict(totals=res, per_family=per, named=named), open(os.path.join(a.outdir, f'score.{a.arm}.seed{a.seed}.json'), 'w'), indent=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd')
    ap.add_argument('species', choices=['human', 'gorilla'])
    ap.add_argument('outdir')
    ap.add_argument('arm', nargs='?')
    ap.add_argument('paf', nargs='*')
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    prep(a) if a.cmd == 'prep' else score(a)


if __name__ == '__main__':
    main()
