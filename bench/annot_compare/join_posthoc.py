#!/usr/bin/env python3
"""Post hoc tables for rows 1707-1709 (docs/PREREG_join_rules_2026-11-14.md; none of this decides a registered bar):
class (b) list per crossing intron, tandem pairs with a spliced model each in separate loci, cost of the F1 tag inside
one gene, chimeric pairs by span relation. Reads the outputs of join_truth.py / chimeric_report.py under PREFIX.
usage: join_posthoc.py X.gtf REFSEQ.gtf[.gz] CAT.genes.tsv GENOME.fa PREFIX"""
import collections, csv, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fusion_audit as FA
from join_truth import tx_table, RC
from cat_vs_refseq_lib import Idx, share

rel = lambda a, b: 'other chromosome' if a[0] != b[0] else 'overlapping spans' if a[1] < b[2] and b[1] < a[2] else 'distant >= 1 Mb' if max(a[1], b[1]) - min(a[2], b[2]) >= 1000000 else 'adjacent, same strand' if a[3] == b[3] else 'adjacent, opposite strands'


def main(xgtf, rs, cat, fasta, out):
    import pysam
    R, rmeta = FA.load(rs, 'gene_id'); name = {g: rmeta[g][2].removeprefix('gene-') for g in R}; ri = Idx(R); byname = {v: k for k, v in name.items()}
    gin = collections.defaultdict(set); tex = collections.defaultdict(list)
    for l in FA.op(rs):
        f = l.rstrip('\n').split('\t')
        if len(f) > 8 and f[2] == 'exon': tex[(FA.A('gene_id', f[8]), FA.A('transcript_id', f[8]))].append((int(f[3]) - 1, int(f[4])))
    for (g, _), e in tex.items(): e.sort(); gin[g].update((a[1], b[0]) for a, b in zip(e, e[1:]))
    X = tx_table(xgtf); own = collections.defaultdict(list)  # gene -> (locus, spliced on an annotated intron of the gene, tagged)
    for t, v in X.items():
        hits = ri.hits(v['c'], v['s'], [list(x) for x in v['ex']])
        if any(not share(R[a][2], R[b][2]) for i, a in enumerate(hits) for b in hits[i + 1:]): continue
        for g in hits: own[g].append((v['g'], bool(gin[g] & set(v['chain'][2])), v['tag']))
    # tandem pairs: a spliced own model each, in two loci
    rows = list(csv.DictReader(open(out + '.tandem.tsv'), delimiter='\t')); c = collections.Counter()
    with open(out + '.tandem.posthoc.tsv', 'w') as fh:
        print(*rows[0], 'spliced_model_a', 'spliced_model_b', 'separate_loci', sep='\t', file=fh)
        for r in rows:
            la, lb = ({l for l, sp, _ in own[byname[r[k]]] if sp} for k in ('gene_a', 'gene_b'))
            sep = bool(la and lb and not (la & lb))
            c['pairs'] += 1; c['two models'] += r['outcome'] == 'two models'; c['both with a spliced model'] += bool(la and lb); c['both spliced, separate loci'] += sep
            print(*r.values(), int(bool(la)), int(bool(lb)), int(sep), sep='\t', file=fh)
    print('tandem pairs (post hoc):', dict(c))
    # F1 tag inside one gene: genes left with no untagged model (the families input drops tagged transcripts)
    tg = collections.defaultdict(lambda: [0, 0])
    for g, ms in own.items():
        for _, _, tag in ms: tg[g][tag] += 1
    hit = [g for g, (u, t) in tg.items() if t]
    print(f'F1 (post hoc): genes with a tagged non-joining transcript {len(hit)}; of those with no untagged transcript {sum(tg[g][0] == 0 for g in hit)}')
    # class (b) list, one row per crossing intron
    fa = pysam.FastaFile(fasta); C = collections.defaultdict(list)
    for r in csv.DictReader(open(cat), delimiter='\t'): C[r['chrom']].append((int(r['start0']), int(r['end']), r['gene_name']))
    catat = lambda ch, p: ','.join(sorted({n for s, e, n in C[ch] if s <= p < e})[:4]) or '.'
    I = collections.defaultdict(list)
    for r in csv.DictReader(open(out + '.joins.tsv'), delimiter='\t'):
        if r['class'] in ('b1', 'b2', 'b0'): I[(r['chrom'], int(r['intron_start']), int(r['intron_end']), r['strand'])].append(r)
    n = collections.Counter()
    with open(out + '.class_b.tsv', 'w') as fh:
        print('chrom\tintron_start\tintron_end\tstrand\tclass\tprimary_reads\tprimary_mapq0\tsecondary_reads\tsplice_motif\tother_sample_primary_reads\tother_sample_flank_depth\ttranscripts\tF1_tagged\tin_final\trefseq_left\trefseq_right\tcat_left\tcat_right', file=fh)
        for (ch, s, e, st), rs_ in sorted(I.items(), key=lambda kv: -int(kv[1][0]['primary'])):
            r = rs_[0]; d, a = fa.fetch(ch, s, s + 2).upper(), fa.fetch(ch, e - 2, e).upper()
            if st == '-': d, a = a.translate(RC)[::-1], d.translate(RC)[::-1]
            n[r['class']] += 1
            print(ch, s, e, st, r['class'], r['primary'], r['primary_mapq0'], r['secondary'], d + '-' + a, r['witness_primary'], r['witness_flank'], len(rs_), sum(int(x['f1_tagged']) for x in rs_), sum(int(x['in_final']) for x in rs_),
                  r['left_genes'] or '.', r['right_genes'] or '.', catat(ch, s - 1), catat(ch, e), sep='\t', file=fh)
    print('class (b) crossing introns:', dict(n))
    # chimeric pairs by span relation
    sp = {name[g]: (R[g][0], R[g][2][0][0], R[g][2][-1][1], R[g][1]) for g in R}; k = collections.Counter()
    for r in csv.DictReader(open(out + '.chimeric.tsv'), delimiter='\t'):
        if int(r['molecules']) < 5 or r['gene_a'] not in sp or r['gene_b'] not in sp: continue
        x = rel(sp[r['gene_a']], sp[r['gene_b']]); k[('>= 5 molecules', x)] += 1; k[('candidate', x)] += r['candidate'] == '1'
    for x in sorted({b for _, b in k}): print('chimeric', x, '>= 5 molecules', k[('>= 5 molecules', x)], 'candidates', k[('candidate', x)], sep='\t')


if __name__ == '__main__':
    assert rel(('c', 0, 10, '+'), ('c', 5, 20, '-')) == 'overlapping spans' and rel(('c', 0, 10, '+'), ('c', 20, 30, '-')) == 'adjacent, opposite strands' and rel(('c', 0, 10, '+'), ('d', 0, 10, '+')) == 'other chromosome'
    main(*sys.argv[1:6])
