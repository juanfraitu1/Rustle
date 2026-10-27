"""PREREG GUIDED-S (docs/LEDGER_O1.md): spliced seeds for the guided family search.

  cdsseed.py prep OUTDIR [--level simiiformes] [--seed 1|2]   truth families, proxy pseudogenes, seeds (cds.fa, prot.fa, body.regions)
  cdsseed.py score OUTDIR ARM HITS [--level ...] [--seed 1|2]  ARM in A (PAF, gene body) / B (PAF, CDS) / C (miniprot GFF)

prep needs `samtools` for the sequences. Hits: a record covering >= 0.5 of the seed at identity >= 0.5 (Liftoff -a / -s defaults);
its gene = the CAT gene with the most exonic bp overlap, priority protein-coding > pseudogene > other; the seed's own gene excluded.
"""
import argparse
import bisect
import collections
import json
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
import famtruth as F  # noqa: E402

FA = '/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa'
COV, IDENT = 0.5, 0.5
CODON = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in 'TCAG' for y in 'TCAG' for z in 'TCAG'],
    'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')}
CIG = re.compile(r'(\d+)([MIDNSHP=X])')


HGNC = '/mnt/linuxdisk/tmp/hgnc/hgnc_complete_set.txt'


def hgnc_groups():
    """Ensembl gene id -> its most specific HGNC gene group id (fewest members; ties lowest id)."""
    import csv
    rows = list(csv.DictReader(open(HGNC), delimiter='\t'))
    size = collections.Counter(g for r in rows for g in r['gene_group_id'].split('|') if g)
    out = {}
    for r in rows:
        gs = [g for g in r['gene_group_id'].split('|') if g]
        if gs and r['ensembl_gene_id']:
            out[r['ensembl_gene_id']] = min(gs, key=lambda g: (size[g], int(g)))
    return out


def families(G, level, truth='compara'):
    pairs, universe = F.truth_pairs(G, F.PARALOGS, F.LEVELS[level])
    if truth == 'hgnc':
        grp = hgnc_groups()
        by = collections.defaultdict(list)
        for g, d in G.items():
            if d['biotype'] == 'protein_coding':
                for e in d['ens']:
                    if e in grp:
                        by[grp[e]].append(g)
        pairs = set(pairs)
        for gs in by.values():
            gs = sorted(set(gs))
            pairs.update((gs[0], x) for x in gs[1:])
        universe = set(universe) | {g for gs in by.values() for g in gs}
    fams = [c for c in F.components(pairs) if len(c) >= 2]
    fams.sort(key=lambda c: min(min(G[g]['ens']) for g in c))
    return fams, universe


def stem_keys(name):
    keys = {name}
    m = re.match(r'^(.*?)([A-Z]|[0-9]+)$', name)
    if m and m.group(1):
        keys.add(m.group(1))
    return keys


def hgnc_pseudogenes(G, fams):
    grp = hgnc_groups()
    fam_grp = collections.defaultdict(set)
    for i, c in enumerate(fams):
        for g in c:
            for e in G[g]['ens']:
                if e in grp:
                    fam_grp[grp[e]].add(i)
    out = collections.defaultdict(set)
    for g, d in G.items():
        if 'pseudogene' in d['biotype']:
            hit = {i for e in d['ens'] if e in grp for i in fam_grp.get(grp[e], ())}
            if len(hit) == 1:
                out[next(iter(hit))].add(g)
    return out


def proxy_pseudogenes(G, fams):
    key_fam = collections.defaultdict(set)
    for i, c in enumerate(fams):
        for g in c:
            for k in stem_keys(G[g]['name']):
                key_fam[k].add(i)
    out = collections.defaultdict(set)
    for g, d in G.items():
        if 'pseudogene' not in d['biotype']:
            continue
        m = re.match(r'^(.+?)P[0-9]*$', d['name'])
        if not m:
            continue
        hit = key_fam.get(m.group(1), set())
        if len(hit) == 1:
            out[next(iter(hit))].add(g)
    return out


def cat_cds(genes_wanted):
    """CAT gene id -> {transcript: [(chrom, s0, e, strand)]} for protein-coding CDS lines of the wanted genes."""
    out = collections.defaultdict(lambda: collections.defaultdict(list))
    for line in open(F.CAT_GFF):
        f = line.split('\t', 8)
        if len(f) < 9 or f[2] != 'CDS':
            continue
        a = dict(kv.split('=', 1) for kv in f[8].rstrip('\n;').split(';') if '=' in kv)
        g = a.get('gene_id')
        if g in genes_wanted:
            out[g][a['Parent']].append((f[0], int(f[3]) - 1, int(f[4]), f[6]))
    return out


def faidx(regions, path, rc=False):
    with open(path + '.regions', 'w') as w:
        w.write(''.join(r + '\n' for r in regions))
    cmd = ['samtools', 'faidx', FA, '-r', path + '.regions'] + (['-i', '--mark-strand', 'no'] if rc else [])
    seqs, name = {}, None
    for line in subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.splitlines():
        if line.startswith('>'):
            name = line[1:]
            seqs[name] = []
        else:
            seqs[name].append(line.upper())
    return {k: ''.join(v) for k, v in seqs.items()}


def prep(a):
    G = F.load_cat()
    F.GENES.update(G)
    fams, universe = families(G, a.level, a.truth)
    pseudo = proxy_pseudogenes(G, fams)
    if a.truth == 'hgnc':
        for i, v in hgnc_pseudogenes(G, fams).items():
            pseudo[i] |= v
    pc = {g for c in fams for g in c}
    cds = cat_cds(pc)
    seeds = {}
    for i, c in enumerate(fams):
        order = sorted((g for g in c if cds.get(g)), key=lambda g: (min(G[g]['ens']), g))
        if len(order) >= a.seed:
            g = order[a.seed - 1]
            tx = max(cds[g], key=lambda t: (sum(e - s for _, s, e, _ in cds[g][t]), t))
            seeds[i] = (g, tx, sorted(cds[g][tx], key=lambda x: x[1]))
    plus = [f'{c}:{s + 1}-{e}' for g, tx, segs in seeds.values() for c, s, e, st in segs if st == '+']
    minus = [f'{c}:{s + 1}-{e}' for g, tx, segs in seeds.values() for c, s, e, st in segs if st == '-']
    sp = faidx(plus, os.path.join(a.outdir, 'plus'))
    sm = faidx(minus, os.path.join(a.outdir, 'minus'), rc=True)
    meta = {}
    with open(os.path.join(a.outdir, 'cds.fa'), 'w') as fc, open(os.path.join(a.outdir, 'prot.fa'), 'w') as fp, \
            open(os.path.join(a.outdir, 'body.regions'), 'w') as fb:
        for i, (g, tx, segs) in seeds.items():
            st = segs[0][3]
            parts = [sp[f'{c}:{s + 1}-{e}'] for c, s, e, _ in segs] if st == '+' else \
                [sm[f'{c}:{s + 1}-{e}'] for c, s, e, _ in reversed(segs)]
            seq = ''.join(parts)
            prot = ''.join(CODON.get(seq[k:k + 3], 'X') for k in range(0, len(seq) - 2, 3)).rstrip('*').replace('*', 'X')
            name = f'F{i}|{g}'
            fc.write(f'>{name}\n{seq}\n')
            fp.write(f'>{name}\n{prot}\n')
            d = G[g]
            fb.write(f"{d['chrom']}:{min(s for s, _ in d['exons']) + 1}-{max(e for _, e in d['exons'])}\t{name}\t{d['strand']}\n")
            meta[name] = dict(fam=i, gene=g, cds_len=len(seq))
    rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(a.outdir, 'body.regions'))]
    bp = faidx([r for r, _, st in rows if st != '-'], os.path.join(a.outdir, 'bodyp'))
    bm = faidx([r for r, _, st in rows if st == '-'], os.path.join(a.outdir, 'bodym'), rc=True)
    with open(os.path.join(a.outdir, 'body.fa'), 'w') as w:  # transcript orientation, like the CDS seeds
        for r, name, st in rows:
            w.write(f'>{name}\n{(bm if st == "-" else bp)[r]}\n')
    json.dump(dict(level=a.level, seed=a.seed, truth=a.truth, universe=sorted(universe), families=[sorted(c) for c in fams], pseudo={str(k): sorted(v) for k, v in pseudo.items()},
                   seeds=meta), open(os.path.join(a.outdir, f'truth.{a.truth}.seed{a.seed}.json'), 'w'))
    print(f'families {len(fams)}, seeded {len(seeds)}, proxy pseudogenes {sum(map(len, pseudo.values()))} in {len(pseudo)} families')


def paf_records(path):
    """(seed name, chrom, strand, genomic blocks, query coverage, identity) per PAF record."""
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        q, ql, qs, qe, st, t, ts = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[7])
        cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), '')
        bl, pos, cur = [], ts, None
        for n, op in CIG.findall(cg):
            n = int(n)
            if op in 'M=XD':
                cur = [pos, pos + n] if cur is None else [cur[0], pos + n]
                pos += n
            elif op == 'N':
                if cur:
                    bl.append(tuple(cur))
                cur = None
                pos += n
        if cur:
            bl.append(tuple(cur))
        yield q, t, st, bl, (qe - qs) / ql, int(f[9]) / int(f[10])


def miniprot_records(path):
    cur = None
    for line in open(path):
        if line.startswith('##PAF'):
            if cur:
                yield cur
            f = line.rstrip('\n').split('\t')
            q, ql, qs, qe = f[1], int(f[2]), int(f[3]), int(f[4])
            cur = [q, f[6], f[5], [], (qe - qs) / ql, None]
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or cur is None:
            continue
        if f[2] == 'mRNA':
            m = re.search(r'Identity=([0-9.]+)', f[8])
            cur[5] = float(m.group(1)) if m else 0.0
        elif f[2] == 'CDS':
            cur[3].append((int(f[3]) - 1, int(f[4])))
    if cur:
        yield cur


def gene_of(idx, G, chrom, strand, bl):
    best = None
    for prio in ('protein_coding', 'pseudogene', None):
        share = collections.Counter()
        if (chrom, strand) not in idx:
            return None
        v, st, m = idx[(chrom, strand)]
        for s, e in bl:
            i = bisect.bisect_left(st, s - m)
            while i < len(v) and v[i][0] < e:
                x, y, g = v[i]
                if y > s and (prio is None or (prio == 'protein_coding' and G[g]['biotype'] == prio) or
                              (prio == 'pseudogene' and 'pseudogene' in G[g]['biotype'])):
                    share[g] += min(y, e) - max(x, s)
                i += 1
        if share:
            best = max(share, key=lambda g: (share[g], g))
            break
    return best


def score(a):
    T = json.load(open(os.path.join(a.outdir, f'truth.{a.truth}.seed{a.seed}.json')))
    G = F.load_cat()
    F.GENES.update(G)
    idx = F.exon_index(G)
    fams = [set(c) for c in T['families']]
    pseudo = {int(k): set(v) for k, v in T['pseudo'].items()}
    universe = set(T['universe']) if 'universe' in T else F.truth_pairs(G, F.PARALOGS, F.LEVELS[T['level']])[1]
    judge = set(universe) | {g for c in fams for g in c} | {g for v in pseudo.values() for g in v}
    recs = miniprot_records(a.hits) if a.arm == 'C' else paf_records(a.hits)
    found = collections.defaultdict(set)
    for q, chrom, strand, bl, cov, ident in recs:
        if q not in T['seeds'] or cov < COV or (ident or 0) < IDENT or not bl:
            continue
        g = gene_of(idx, G, chrom, strand, bl)
        if g and g != T['seeds'][q]['gene']:
            found[q].add(g)
    tot = collections.defaultdict(lambda: collections.Counter())
    per = {}
    for q, m in T['seeds'].items():
        i = m['fam']
        others = fams[i] - {m['gene']}
        ps = pseudo.get(i, set())
        fd = found.get(q, set())
        half = 'dev' if _half(G, fams[i]) == 0 else 'heldout'
        c = tot[half]
        c['members'] += len(others); c['members_found'] += len(fd & others)
        c['pseudo'] += len(ps); c['pseudo_found'] += len(fd & ps)
        c['judged'] += len(fd & judge); c['judged_ok'] += len(fd & (others | ps)); c['unjudged'] += len(fd - judge)
        c['fams'] += 1; c['macro_member'] += len(fd & others) / len(others) if others else 0
        if ps:
            c['fams_ps'] += 1; c['macro_pseudo'] += len(fd & ps) / len(ps)
        per[q] = dict(names=sorted(G[g]['name'] for g in fams[i]), found_members=sorted(G[g]['name'] for g in fd & others),
                      missed_members=sorted(G[g]['name'] for g in others - fd), found_pseudo=sorted(G[g]['name'] for g in fd & ps),
                      missed_pseudo=sorted(G[g]['name'] for g in ps - fd), wrong=sorted(G[g]['name'] for g in (fd & judge) - others - ps))
    res = {}
    for half, c in tot.items():
        res[half] = dict(c, member_recall=c['members_found'] / max(c['members'], 1), pseudo_recall=c['pseudo_found'] / max(c['pseudo'], 1),
                         precision=c['judged_ok'] / max(c['judged'], 1), macro_member=c['macro_member'] / max(c['fams'], 1),
                         macro_pseudo=c['macro_pseudo'] / max(c['fams_ps'], 1))
        if half not in a.show.split(','):
            continue
        print(f"{a.arm} seed{a.seed} {half}: member recall {res[half]['member_recall']:.3f} ({c['members_found']}/{c['members']}) "
              f"pseudo recall {res[half]['pseudo_recall']:.3f} ({c['pseudo_found']}/{c['pseudo']}) precision {res[half]['precision']:.3f} "
              f"({c['judged_ok']}/{c['judged']}) unjudged {c['unjudged']} | MACRO member {res[half]['macro_member']:.3f} "
              f"pseudo {res[half]['macro_pseudo']:.3f} ({c['fams']} fams, {c['fams_ps']} with pseudogenes)")
    named = {}
    for label, pat in {**F.MAINSTAYS, **F.YAGS}.items():
        rx = re.compile(pat)
        for q, r in per.items():
            if any(rx.match(n) for n in r['names']):
                named.setdefault(label, []).append(dict(seed=q, **{k: r[k] for k in ('found_members', 'missed_members', 'found_pseudo', 'missed_pseudo', 'wrong')}))
    json.dump(dict(totals=res, per_family=per, named=named), open(os.path.join(a.outdir, f'score.{a.arm}.{a.truth}.seed{a.seed}.json'), 'w'), indent=0)


def _half(G, fam):
    F.G = G
    return F.half(fam)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd')
    ap.add_argument('outdir')
    ap.add_argument('arm', nargs='?')
    ap.add_argument('hits', nargs='?')
    ap.add_argument('--level', default='simiiformes')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--truth', default='compara', choices=['compara', 'hgnc'])
    ap.add_argument('--show', default='dev', help='halves to PRINT (JSON holds both); held-out only after dev selection')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    prep(a) if a.cmd == 'prep' else score(a)


if __name__ == '__main__':
    main()
