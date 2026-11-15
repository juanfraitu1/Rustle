#!/usr/bin/env python3
"""Rows 1707-1709 (docs/PREREG_join_rules_2026-11-14.md): classes of gene-joining transcripts and what each rule does to them.
A join = a transcript with a crossing intron (left exon on one RefSeq gene of any biotype, right exon on another that shares
no exon with it). Class from the weakest crossing intron: a annotated (RefSeq or CAT intron), b1 / b2 / b0 unannotated with
>= 5 primary records and canonical sites (recurring in the witness BAM / not / witness not informative), w weak (2-4), c.
usage: join_truth.py X.gtf FINAL.gtf REFSEQ.gtf[.gz] CAT.gff3[.gz] BAM WITNESS.bam GENOME.fa OUTPREFIX   (pysam; famtruth env)
writes OUTPREFIX.joins.tsv, .tagged.tsv, .rtgenes.tsv, .tandem.tsv and prints the tables."""
import collections, csv, gzip, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
import fusion_audit as FA
from join_support import crossing
from cat_vs_refseq_lib import Idx, share

N_SUP, N_WIT, N_WEAK, N_FLANK = 5, 2, 2, 5
CANON = {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')}
RC = str.maketrans('ACGT', 'TGCA')


def klass(annot, prim, canon, wit, flank):
    """frozen class rule of the PREREG."""
    if annot: return 'a'
    if prim >= N_SUP and canon: return 'b1' if wit >= N_WIT else 'b2' if flank >= N_FLANK else 'b0'
    return 'w' if canon and prim >= N_WEAK else 'c'


def tx_table(gtf):
    T = {}
    for l in open(gtf):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9: continue
        t = FA.A('transcript_id', f[8])
        if f[2] == 'transcript': T[t] = dict(c=f[0], s=f[6], g=FA.A('gene_id', f[8]), reads=int(FA.A('reads', f[8]) or 0), tag='fusion_of "' in f[8], ex=[])
        elif f[2] == 'exon': T[t]['ex'].append((int(f[3]) - 1, int(f[4])))
    for v in T.values(): v['ex'].sort(); v['chain'] = (v['c'], v['s'], tuple((a[1], b[0]) for a, b in zip(v['ex'], v['ex'][1:])))
    return T


class Bam:
    def __init__(self, path):
        import pysam
        self.b = pysam.AlignmentFile(path); self.c = {}

    def intron(self, c, s, e):
        """(primary, primary MAPQ 0, secondary) records carrying exactly intron [s, e)."""
        k = (c, s, e)
        if k not in self.c:
            n = [0, 0, 0]
            for r in self.b.fetch(c, max(s - 1, 0), s):
                if r.is_supplementary or r.cigartuples is None: continue
                p = r.reference_start
                for op, ln in r.cigartuples:
                    if op == 3 and p == s and p + ln == e:
                        n[2 if r.is_secondary else 0] += 1; n[1] += (not r.is_secondary and r.mapping_quality == 0); break
                    if op in (0, 2, 3, 7, 8): p += ln
                    if p > s: break
            self.c[k] = n
        return self.c[k]

    def flank(self, c, s, e):
        """the smaller count of primary records over the base before the donor and the base after the acceptor."""
        f = lambda x: self.b.count(c, max(x, 0), max(x, 0) + 1, read_callback=lambda r: not (r.is_secondary or r.is_supplementary))
        return min(f(s - 1), f(e))


def main(xgtf, fgtf, rs, cat, bam, wit, fasta, out):
    import pysam
    R, rmeta = FA.load(rs, 'gene_id'); name = {g: rmeta[g][2].removeprefix('gene-') for g in R}; ri = Idx(R)
    tex = collections.defaultdict(list); ann = set(); gin = collections.defaultdict(set)
    for l in FA.op(rs):
        f = l.rstrip('\n').split('\t')
        if len(f) > 8 and f[2] == 'exon': tex[(FA.A('gene_id', f[8]), FA.A('transcript_id', f[8]), f[0], f[6])].append((int(f[3]) - 1, int(f[4])))
    for (g, _, c, _), e in tex.items():
        e.sort(); I = {(c, a[1], b[0]) for a, b in zip(e, e[1:])}; ann |= I; gin[g] |= I
    cex = collections.defaultdict(list)
    for l in FA.op(cat):
        f = l.rstrip('\n').split('\t')
        if len(f) > 8 and f[2] == 'exon': cex[(f[0], re.search(r'Parent=([^;]+)', f[8]).group(1))].append((int(f[3]) - 1, int(f[4])))
    ncat = 0
    for (c, _), e in cex.items():
        e.sort(); I = {(c, a[1], b[0]) for a, b in zip(e, e[1:])}; ncat += len(I - ann); ann |= I
    refseq_only = set().union(*gin.values())
    B, W, fa = Bam(bam), Bam(wit), pysam.FastaFile(fasta)

    def canon(c, s, e, strand):
        d, a = fa.fetch(c, s, s + 2).upper(), fa.fetch(c, e - 2, e).upper()
        return ((d, a) if strand == '+' else (a.translate(RC)[::-1], d.translate(RC)[::-1])) in CANON

    X, F = tx_table(xgtf), tx_table(fgtf); fchains = {v['chain'] for v in F.values()}; xchains = {v['chain'] for v in X.values()}
    own = set(); touch = collections.defaultdict(set)  # gene -> transcripts touching it without a disjoint second gene

    def judge(v):
        hits = ri.hits(v['c'], v['s'], [list(x) for x in v['ex']]); cr = crossing(v['ex'], hits, R)
        dis = any(not share(R[a][2], R[b][2]) for i, a in enumerate(hits) for b in hits[i + 1:])
        return hits, cr, dis

    n = collections.Counter(); rows = []; tagged = []; jgenes = collections.defaultdict(set)
    prev = {}  # shortcut: reuses any existing OUT.joins.tsv without checking it came from the same inputs, delete it to recount
    if os.path.exists(out + '.joins.tsv'):
        for l in list(open(out + '.joins.tsv'))[1:]:
            f = l.rstrip('\n').split('\t'); prev[f[0]] = tuple(int(x) if i in (3, 4) or 7 <= i <= 17 else x for i, x in enumerate(f))
    for t, v in X.items():
        hits, cr, dis = judge(v)
        if not dis: own.update(hits)
        if v['tag']:
            on = [{g for g in hits if share([list(e)], R[g][2])} for e in v['ex']]
            kind = 'join' if cr else 'inside one gene' if on and set.intersection(*on) else 'touches two genes, no crossing intron' if dis else 'exon outside any gene'
            tagged.append((t, v['c'], v['ex'][0][0], v['ex'][-1][1], v['s'], v['reads'], kind)); n[('tagged', kind)] += 1
        if not cr:
            n['touch two genes, no crossing intron'] += dis; continue
        for g in hits: jgenes[g].add(t)
        if t in prev: rows.append(prev[t]); continue
        best = None
        for s, e in cr:
            p, p0, sec = B.intron(v['c'], s, e)
            if best is None or p < best[2]: best = (s, e, p, p0, sec)
        s, e, p, p0, sec = best; an = (v['c'], s, e) in ann; cn = canon(v['c'], s, e, v['s']); w = W.intron(v['c'], s, e)[0]
        fl = W.flank(v['c'], s, e) if (not an and p >= N_SUP and cn and w < N_WIT) else -1
        k = klass(an, p, cn, w, fl); infin = v['chain'] in fchains; flag = (p + sec) > 0 and sec / (p + sec) >= 0.5
        on = [{g for g in hits if share([list(x)], R[g][2])} for x in v['ex']]; i = [x[1] for x in v['ex']].index(s)
        ga, gb = sorted(name[g] for g in on[i] - on[i + 1]), sorted(name[g] for g in on[i + 1] - on[i])
        rows.append((t, v['g'], v['c'], s, e, v['s'], k, int(an), int((v['c'], s, e) in refseq_only), p, p0, sec, int(cn), w, fl, int(v['tag']), int(infin), int(flag), v['reads'], ','.join(ga), ','.join(gb)))
    hdr = 'transcript\tlocus\tchrom\tintron_start\tintron_end\tstrand\tclass\tannotated\tin_refseq\tprimary\tprimary_mapq0\tsecondary\tcanonical\twitness_primary\twitness_flank\tf1_tagged\tin_final\tsec_flag\treads\tleft_genes\tright_genes'
    with open(out + '.joins.tsv', 'w') as fh:
        print(hdr, file=fh)
        for r in sorted(rows, key=lambda r: (r[6], -r[9])): print(*r, sep='\t', file=fh)
    with open(out + '.tagged.tsv', 'w') as fh:
        print('transcript\tchrom\tstart\tend\tstrand\treads\tkind', file=fh)
        for r in tagged: print(*r, sep='\t', file=fh)
    CL = ('a', 'b1', 'b2', 'b0', 'w', 'c'); col = hdr.split('\t').index; cnt = lambda f: [sum(1 for r in rows if r[col('class')] == k and f(r)) for k in CL]
    print(f'{xgtf}\ntranscripts {len(X)}; joins (crossing intron) {len(rows)}; touching two genes without one {n["touch two genes, no crossing intron"]}; CAT-only annotated introns {ncat}')
    print('row', *CL, 'all', sep='\t')
    for lab, f in (('joins', lambda r: 1), ('F1-tagged', lambda r: r[col('f1_tagged')]), ('absent from final.gtf', lambda r: not r[col('in_final')]),
                   ('secondary share >= 0.5', lambda r: r[col('sec_flag')]), ('primaries mostly MAPQ 0', lambda r: r[col('primary')] and 2 * r[col('primary_mapq0')] >= r[col('primary')])):
        x = cnt(f); print(lab, *x, sum(x), sep='\t')
    tg = collections.Counter(k for *_, k in tagged); nt = sum(tg.values())
    print('F1-tagged transcripts by content:', dict(tg), f'inside one gene {tg["inside one gene"] / max(nt, 1):.3f} (too strict if > 0.10)')
    ab = [r for r in rows if r[col('class')] in ('a', 'b1')]; cc = [r for r in rows if r[col('class')] == 'c']; fl = [r for r in rows if r[col('sec_flag')]]
    sh = lambda a, b: f'{a} of {b} = {a / max(b, 1):.3f}'
    print('F1 leaves untagged a + b1 joins:', sh(sum(not r[col('f1_tagged')] for r in ab), len(ab)), '(too loose if > 0.50)')
    print('final.gtf removes a + b1 joins:', sh(sum(not r[col('in_final')] for r in ab), len(ab)), '(too strict if > 0.05); class c joins passing:', sh(sum(r[col('in_final')] for r in cc), len(cc)), '(too loose if > 0.50)')
    print('secondary-share emulation flags a or b1:', sh(sum(r[col('class')] in ('a', 'b1') for r in fl), len(fl)), '(too strict if > 0.10)')
    fj = 0
    for v in F.values():
        if v['chain'] not in xchains and judge(v)[1]: fj += 1
    print('joins in final.gtf whose chain is not in x.gtf:', fj)
    for nm in ('PKD1P6-NPIPP1', 'PDXDC2P-NPIPB14P', 'PKD1P4-NPIPA8', 'PKD1P3-NPIPA1'):
        g = next((g for g in R if name[g] == nm), None)
        if g is None: print('control', nm, 'not in the annotation'); continue
        res = []
        for (gg, _, c, s_), e in tex.items():
            if gg != g: continue
            for s, e2 in crossing(e, [h for h in ri.hits(c, s_, [list(x) for x in e]) if h != g], R): res.append((s, e2, *B.intron(c, s, e2)))
        res = sorted(set(res)); print('control', nm, [(f'{s}-{e}', f'prim {p}', f'sec {q}', 'FLAG' if p + q and q / (p + q) >= 0.5 else 'keep') for s, e, p, _, q in res][:4] or 'no crossing intron')
    # row 1708: annotated readthrough genes
    xi, fi = ({(v['c'], *i) for v in T.values() for i in v['chain'][2]} for T in (X, F)); rt = []
    done = collections.defaultdict(set)
    for (g, _, c, s_), e in tex.items():
        for s, e2 in crossing(e, [h for h in ri.hits(c, s_, [list(x) for x in e]) if h != g], R): done[g].add((c, s, e2))
    old = out + '.rtgenes.tsv'
    if os.path.exists(old): rt = [(f[0], *map(int, f[1:])) for f in (l.split('\t') for l in list(open(old))[1:])]; done = {}
    for g, I in done.items():
        sup = {i: B.intron(*i)[0] for i in I}; ex_ = [i for i in I if sup[i] >= N_SUP]
        rt.append((name[g], len(I), max(sup.values()), int(bool(ex_)), int(any(i in xi for i in ex_)), int(any(i in fi for i in ex_)), int(any(i in xi for i in I))))
    with open(out + '.rtgenes.tsv', 'w') as fh:
        print('gene\tcrossing_introns\tmax_primary\texpressed\tin_x\tin_final\tany_intron_in_x', file=fh)
        for r in sorted(rt, key=lambda r: -r[2]): print(*r, sep='\t', file=fh)
    ne = sum(r[3] for r in rt)
    print(f'annotated readthrough genes (structural) {len(rt)}; expressed {ne}; recall x.gtf {sh(sum(r[4] for r in rt), ne)}; final.gtf {sh(sum(r[5] for r in rt), ne)} (bar 0.80); not expressed but modelled {sum(r[6] and not r[3] for r in rt)}')
    # row 1708: tandem truth pairs
    try:
        import famtruth as T_
        G = T_.load_cat(); pairs, _ = T_.truth_pairs(G, T_.PARALOGS, set(T_.LEVELS['primates']))
        byname = {name[g]: g for g in R}; tp = collections.Counter(); trows = []
        expr = lambda g: any(B.intron(*i)[0] >= N_SUP for i in sorted(gin[g])[:40])
        for a, b in pairs:
            a, b = byname.get(G[a]['name']), byname.get(G[b]['name'])
            if not a or not b or R[a][:2] != R[b][:2] or share(R[a][2], R[b][2]): continue
            gap = max(R[a][2][0][0], R[b][2][0][0]) - min(R[a][2][-1][1], R[b][2][-1][1])
            if gap > 1000000 or not gin[a] or not gin[b] or not (expr(a) and expr(b)): continue
            fused = bool(jgenes[a] & jgenes[b]); oa, ob = a in own, b in own
            o = 'two models' if oa and ob else 'fused only' if fused else 'both missing' if not (oa or ob) else 'one missing'
            tp[o] += 1; tp['two models and also a join'] += (oa and ob and fused); trows.append((name[a], name[b], R[a][0], gap, o, int(fused)))
        with open(out + '.tandem.tsv', 'w') as fh:
            print('gene_a\tgene_b\tchrom\tgap\toutcome\tjoined_by_a_transcript', file=fh)
            for r in sorted(trows): print(*r, sep='\t', file=fh)
        nt = len(trows); print(f'tandem truth pairs (same strand, <= 1 Mb, both spliced and expressed) {nt}: ' + ', '.join(f'{k} {tp[k]}' for k in ('two models', 'fused only', 'one missing', 'both missing', 'two models and also a join')) + f'; two models {tp["two models"] / max(nt, 1):.3f} (bar 0.80)')
    except (Exception, SystemExit) as e:
        print(f'tandem pairs: no paralog truth ({type(e).__name__}: {e})')


def _selfcheck():
    assert klass(True, 0, False, 0, 0) == 'a' and klass(False, 5, True, 2, -1) == 'b1' and klass(False, 5, True, 1, 5) == 'b2'
    assert klass(False, 5, True, 0, 4) == 'b0' and klass(False, 4, True, 9, 9) == 'w' and klass(False, 1, True, 0, 0) == 'c' and klass(False, 50, False, 9, 9) == 'c'


if __name__ == '__main__':
    _selfcheck(); main(*sys.argv[1:9])
