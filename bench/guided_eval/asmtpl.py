"""PREREG ASMTPL (docs/LEDGER_O1.md): assembled transcripts as family templates, gorilla, 13 keyword families.

  asmtpl.py prep OUTDIR REGIONS_TSV GTF...     transcripts of the assembled GTFs at the family loci -> unique intron chains -> templates.fa
  asmtpl.py score OUTDIR ARM PAF...            ARM R1 = one template per family (most reads, then longest); RK = every unique chain
"""
import collections
import json
import os
import re
import sys
import urllib.parse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cdsseed as C  # noqa: E402
import famtruth as F  # noqa: E402
import guided3 as g3  # noqa: E402

KEY = {'NPIP': r'nuclear pore complex[- ]interacting protein', 'TBC1D3': r'TBC1 domain family member 3',
       'GOLGA6': r'golgin subfamily A member 6', 'GOLGA8': r'golgin subfamily A member 8',
       'TSPY': r'testis-specific Y-encoded protein', 'XAGE': r'X antigen family', 'DAZ': r'deleted in azoospermia protein',
       'CT47': r'cancer/testis antigen 47', 'SPATA31': r'spermatogenesis-associated protein 31',
       'RBMY': r'RNA-binding motif protein, Y chromosome', 'CDY': r'testis-specific chromodomain protein Y',
       'NBPF': r'neuroblastoma breakpoint family', 'USP17L': r'ubiquitin carboxyl-terminal hydrolase 17-like'}
SHARED = os.environ.get('ASMTPL_SHARED') == '1'
EXON_JUDGE = os.environ.get('ASMTPL_EXON') == '1'  # ASMTPL-3
MERGE = {'GOLGA6': 'GOLGA', 'GOLGA8': 'GOLGA'} if EXON_JUDGE else {}
SYM = {'NPIP': r'^NPIP', 'TBC1D3': r'^TBC1D3', 'GOLGA6': r'^GOLGA6', 'GOLGA8': r'^GOLGA8', 'TSPY': r'^TSPY', 'XAGE': r'^XAGE',
       'DAZ': r'^DAZ[0-9]', 'CT47': r'^CT47', 'SPATA31': r'^SPATA31', 'RBMY': r'^RBMY', 'CDY': r'^CDY[0-9]', 'NBPF': r'^NBPF', 'USP17L': r'^USP17L'}
DEV = {'NPIP', 'TBC1D3', 'GOLGA6', 'GOLGA8', 'USP17L', 'NBPF'}
FA = '/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta'


def keyword_of(G):
    """RefSeq gene id -> keyword family (by description) or None."""
    desc = {}
    for line in open(F.GGO_GFF):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] not in ('gene', 'pseudogene'):
            continue
        a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        desc[a['ID']] = urllib.parse.unquote(a.get('description', ''))
    out = {}
    for g in G:
        d = desc.get(g, '')
        out[g] = next((k for k, p in KEY.items() if re.search(p, d, re.I)), None) or \
            next((k for k, p in SYM.items() if re.match(p, G[g]['name'])), None)  # named members carry HGNC-style descriptions
    return out


def read_gtf(path):
    tx = {}
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] not in ('transcript', 'exon'):
            continue
        a = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
        t = tx.setdefault(a['transcript_id'], dict(chrom=f[0], strand=f[6], exons=[], reads=0))
        if f[2] == 'transcript':
            t['reads'] = int(a.get('reads', 0))
        else:
            t['exons'].append((int(f[3]) - 1, int(f[4])))
    for t in tx.values():
        t['exons'].sort()
    return tx


def prep(outdir, regions, gtfs):
    G = F.load_ggo(F.load_cat())
    kw = keyword_of(G)
    loci = collections.defaultdict(list)
    for line in open(regions):
        c, s, e, fam, gid, bt = line.rstrip('\n').split('\t')
        loci[c].append((int(s), int(e), fam, gid))
    chains = {}
    for p in gtfs:
        for tid, t in read_gtf(p).items():
            s, e = t['exons'][0][0], t['exons'][-1][1]
            hit = {(fam, gid) for ls, le, fam, gid in loci.get(t['chrom'], ()) if ls < e and le > s}
            fams = {fam for fam, _ in hit}
            if len(fams) != 1:
                continue
            key = (t['chrom'], t['strand'], tuple((a, b) for (_, a), (b, _) in zip(t['exons'], t['exons'][1:])) or tuple(t['exons']))
            cur = chains.get(key)
            if cur is None or (t['reads'], sum(b - a for a, b in t['exons'])) > (cur['reads'], cur['len']):
                chains[key] = dict(fam=next(iter(fams)), src=sorted(g for _, g in hit), chrom=t['chrom'], strand=t['strand'],
                                   exons=t['exons'], reads=t['reads'], len=sum(b - a for a, b in t['exons']), n_introns=len(t['exons']) - 1)
    span = collections.defaultdict(list)
    for g, d in G.items():
        span[d['chrom']].append((min(a for a, _ in d['exons']), max(b for _, b in d['exons']), g))
    for t in chains.values():
        s0, e0 = t['exons'][0][0], t['exons'][-1][1]
        over = {g for a, b, g in span[t['chrom']] if a < e0 and b > s0 and any(x < eb and y > ea for x, y in G[g]['exons'] for ea, eb in t['exons'])}
        t['chimeric'] = not SHARED and any(kw.get(g) != t['fam'] and G[g]['biotype'] == 'protein_coding' and not g3.unnamed(G[g]) for g in over)
        t['src'] = sorted(set(t['src']) | {g for g in over if kw.get(g) == t['fam']})
    if SHARED:  # ASMTPL-2: keep family-shared exons only
        sites = collections.Counter()
        for t in chains.values():
            for a, b in t['exons']:
                sites[(t['fam'], t['chrom'], t['strand'], a, b)] += 1
        inner = collections.Counter()
        for t in chains.values():
            ex = t['exons']
            if len(ex) > 1:
                inner[(t['fam'], t['chrom'], t['strand'], 'acc', ex[0][1])] += 1   # first exon: its donor side
                inner[(t['fam'], t['chrom'], t['strand'], 'don', ex[-1][0])] += 1  # last exon: its acceptor side
        kept = {}
        for key, t in chains.items():
            ex = t['exons']
            if len(ex) == 1:
                keep = ex
            else:
                keep = [e for i, e in enumerate(ex) if sites[(t['fam'], t['chrom'], t['strand'], e[0], e[1])] >= 2
                        or (i == 0 and inner[(t['fam'], t['chrom'], t['strand'], 'acc', e[1])] >= 2)
                        or (i == len(ex) - 1 and inner[(t['fam'], t['chrom'], t['strand'], 'don', e[0])] >= 2)]
            if keep:
                t = dict(t, exons=keep, len=sum(b - a for a, b in keep))
                kept.setdefault((t['chrom'], t['strand'], tuple(keep)), t)
        chains = kept
    C.FA = FA
    regs = [f"{t['chrom']}:{a + 1}-{b}" for t in chains.values() for a, b in t['exons']]
    sp = C.faidx(regs, f'{outdir}/ex')
    templates = {}
    with open(f'{outdir}/templates.fa', 'w') as w:
        for i, t in enumerate(chains.values()):
            seq = ''.join(sp[f"{t['chrom']}:{a + 1}-{b}"] for a, b in t['exons'])
            if t['strand'] == '-':
                seq = seq.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]
            name = f"{t['fam']}|T{i}"
            w.write(f'>{name}\n{seq}\n')
            templates[name] = dict(fam=t['fam'], src=t['src'], cds_len=len(seq), reads=t['reads'], n_introns=t['n_introns'], chimeric=t['chimeric'], body_cds=None)
    for fam in {t['fam'] for t in templates.values()}:
        best = max((n for n, t in templates.items() if t['fam'] == fam and not t['chimeric']), default=None, key=lambda n: (templates[n]['reads'], templates[n]['cds_len']))
        if best:
            templates[best]['r1'] = True
    fams = collections.defaultdict(set)
    pseudo = collections.defaultdict(set)
    for g, k in kw.items():
        if k and G[g]['biotype'] == 'protein_coding':
            fams[k].add(g)
        elif k and 'pseudogene' in G[g]['biotype']:
            pseudo[k].add(g)
    json.dump(dict(families={k: sorted(v) for k, v in fams.items()}, pseudo={k: sorted(v) for k, v in pseudo.items()},
                   keyword={g: k for g, k in kw.items() if k}, templates=templates), open(f'{outdir}/truth.json', 'w'))
    per = collections.Counter(t['fam'] for t in templates.values())
    chim = collections.Counter(t['fam'] for t in templates.values() if t['chimeric'])
    print('templates per family:', dict(sorted(per.items())), '| chimeric', dict(sorted(chim.items())), '| total', len(templates))


def score(outdir, arm, pafs):
    T = json.load(open(f'{outdir}/truth.json'))
    G = F.load_ggo(F.load_cat())
    idx = F.exon_index(G)
    tm = T['templates']
    for t in tm.values():
        t['fam'] = MERGE.get(t['fam'], t['fam'])
    kw = {g: MERGE.get(k, k) for g, k in T['keyword'].items()}
    fams, pseudo = collections.defaultdict(set), collections.defaultdict(set)
    for k, v in T['families'].items():
        fams[MERGE.get(k, k)] |= set(v)
    for k, v in T['pseudo'].items():
        pseudo[MERGE.get(k, k)] |= set(v)
    use = {n for n, t in tm.items() if (arm == 'RK' or t.get('r1')) and not t.get('chimeric')}
    src = collections.defaultdict(set)  # fam -> loci expressed in the libraries (reported, not excluded)
    for t in tm.values():
        src[t['fam']] |= set(t['src'])
    groups = collections.defaultdict(list)
    exbp = collections.defaultdict(lambda: [0, 0])
    for q, t, st, bl, qi, ti, cdsb, ident in g3.records(pafs, {}):
        if q not in use or not bl:
            continue
        g = C.gene_of(idx, G, t, st, bl)
        if g and g not in tm[q]['src']:  # a template never finds its own locus
            groups[(q, g, st)].append((qi, ti, cdsb, ident))
            exbp[(q, g, st)][0] += sum(max(0, min(b, e) - max(a, s)) for a, b in bl for s, e in G[g]['exons'])
            exbp[(q, g, st)][1] += sum(b - a for a, b in bl)
    found, putative = collections.defaultdict(set), collections.defaultdict(set)
    for (q, g, st), rs in groups.items():
        bases, ident = g3.chain(rs, st)
        if bases / tm[q]['cds_len'] >= g3.COV and ident >= g3.IDENT:
            if EXON_JUDGE and exbp[(q, g, st)][0] < 0.5 * exbp[(q, g, st)][1]:
                putative[tm[q]['fam']].add(g)  # blocks mostly inside the gene's introns: a copy the annotation does not name
            else:
                found[tm[q]['fam']].add(g)
    rows = {}
    for fam in sorted(src):
        others, ps, fd = fams.get(fam, set()), pseudo.get(fam, set()), found.get(fam, set())
        ok = fd & (others | ps)
        unj = {g for g in fd - ok if kw.get(g) is None and g3.unnamed(G[g])}
        wrong = fd - ok - unj
        rows[fam] = dict(split='dev' if fam in DEV or fam == 'GOLGA' else 'heldout', putative=len(putative.get(fam, set()) - fd), n_templates=sum(1 for n in use if tm[n]['fam'] == fam),
                         members=len(fd & others), n_members=len(others), pseudo=len(fd & ps), n_pseudo=len(ps),
                         ok=len(ok), wrong=len(wrong), unjudged=len(unj), wrong_names=sorted(G[g]['name'] + ':' + str(kw.get(g)) for g in wrong))
        r = rows[fam]
        print(f"{arm} {fam:8s} {r['split']:7s} tpl {r['n_templates']:3d} | members {r['members']}/{r['n_members']} pseudo {r['pseudo']}/{r['n_pseudo']} "
              f"| wrong {r['wrong']} unjudged {r['unjudged']} putative {r['putative']} {r['wrong_names'][:4]}")
    for split in ('dev', 'heldout'):
        rs = [r for r in rows.values() if r['split'] == split]
        ok, wr = sum(r['ok'] for r in rs), sum(r['wrong'] for r in rs)
        print(f"{arm} {split}: member {sum(r['members'] for r in rs)}/{sum(r['n_members'] for r in rs)} pseudo {sum(r['pseudo'] for r in rs)}/"
              f"{sum(r['n_pseudo'] for r in rs)} judged precision {ok / max(ok + wr, 1):.3f} ({ok}/{ok + wr}) unjudged {sum(r['unjudged'] for r in rs)} putative {sum(r['putative'] for r in rs)}")
    json.dump(rows, open(f'{outdir}/score.{arm}.json', 'w'), indent=0)


if __name__ == '__main__':
    if sys.argv[1] == 'prep':
        prep(sys.argv[2], sys.argv[3], sys.argv[4:])
    else:
        score(sys.argv[2], sys.argv[3], sys.argv[4:])
