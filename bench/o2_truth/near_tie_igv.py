#!/usr/bin/env python3
"""r1829 side figure (descriptive; user 2026-11-19): IGV-like examples of undecidable and decidable near-tied reads, the
reads' records at BOTH copies drawn against the annotated transcripts of each copy's gene.

  near_tie_igv.py CLASSES.tsv READS.tsv TIE_SPLIT.tsv SEL.bam ANNOT.gff3.gz CAT_GENES.tsv OUT_PREFIX [N_UNDEC N_DEC]

Examples: the gene pairs with the most reads of class `2 undecidable` (several HG002 copies fit equally) and of the
decidable classes (3 reference lead, 4 own variant), among reads spliced at both places. Per example and copy: the CAT
transcripts of the copy's gene, and every record of the example's reads (at most 40 reads) overlapping the copy: aligned
blocks, mismatched bases (=/X CIGAR), primary or secondary, ms. A read's chain at a copy is KNOWN when its introns are a
contiguous run of the introns of one annotated transcript of that gene (exact coordinates).
Writes OUT_PREFIX.features.json (one row per drawn feature) and OUT_PREFIX.examples.json (one row per example and copy),
each {rows, meta} with this code."""
import collections, gzip, json, os, sys
import pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from near_tie_classes import truth, annot, gene_at

def main(classes, reads_tsv, split_tsv, sel, gff, genes_tsv, out, n_undec=2, n_dec=1):
    n_undec, n_dec = int(n_undec), int(n_dec)
    T = truth(reads_tsv)
    split = {}
    for l in open(split_tsv):
        f = l.rstrip('\n').split('\t')
        if f[0] != 'read':
            split[f[0]] = [int(x) for x in f[7].split(',') if x]
    cls = {}
    for l in open(classes):
        f = l.rstrip('\n').split('\t')
        if f[0] != 'read':
            cls[f[0]] = f[3]
    name = {}; biot = {}
    for l in open(genes_tsv):
        f = l.rstrip('\n').split('\t')
        if f[0] != 'gene_id':
            name[f[0]] = f[1]; biot[f[0]] = f[2]
    A = annot(gff)
    tx, ex, idx, longest = A
    R = collections.defaultdict(list)
    for r in pysam.AlignmentFile(sel, check_sq=False):
        if r.is_unmapped or r.is_supplementary or r.query_name not in cls:
            continue
        mm = []; p = r.reference_start; blocks = []; introns = []; bs = None
        for op, n in r.cigartuples:
            if op in (7, 8, 0):
                if bs is None: bs = p
                if op == 8: mm.extend(range(p, p + n))
                p += n
            elif op == 2:
                if bs is None: bs = p
                p += n
            elif op == 3:
                if bs is not None: blocks.append((bs, p)); bs = None
                introns.append((p, p + n)); p += n
        if bs is not None: blocks.append((bs, p))
        R[r.query_name].append(dict(chrom=r.reference_name, strand='-' if r.is_reverse else '+', s=r.reference_start, e=r.reference_end,
                                    primary=not r.is_secondary, ms=r.get_tag('ms'), AS=r.get_tag('AS'), blocks=blocks, introns=introns, mm=mm))
    gene_tx = collections.defaultdict(list)
    for t, (c, st, g) in tx.items():
        if ex.get(t): gene_tx[g].append(t)
    def tx_introns(t):
        e = sorted(ex[t]); return [(a[1], b[0]) for a, b in zip(e, e[1:])]
    def known(introns, g):
        if not introns: return False
        k = len(introns)
        for t in gene_tx[g]:
            ti = tx_introns(t)
            for i in range(len(ti) - k + 1):
                if ti[i:i + k] == introns: return True
        return False
    def rec_at(q, locus):
        c, s, e, st = locus
        rr = [x for x in R[q] if x['chrom'] == c and x['strand'] == st and x['s'] < e and s < x['e']]
        return max(rr, key=lambda x: x['ms']) if rr else None
    # candidate groups: (gene A, gene B) -> reads, for undecidable and decidable reads spliced at both places
    groups = {'undecidable': collections.defaultdict(list), 'decidable': collections.defaultdict(list)}
    for q, c in cls.items():
        d = T.get(q)
        if d is None: continue
        if c == '2 undecidable':
            places = [i for i in split.get(q, []) if 0 <= i < len(d['loci'])]; kind = 'undecidable'
        elif c in ('3 reference lead', '4 own variant'):
            o = next(iter(d['origin']))
            others = [i for i in range(len(d['loci'])) if i != o]
            recs = [(i, rec_at(q, d['loci'][i])) for i in others]
            recs = [(i, x) for i, x in recs if x]
            if not recs: continue
            places = [o, max(recs, key=lambda y: y[1]['ms'])[0]]; kind = 'decidable'
        else:
            continue
        recs = [(i, rec_at(q, d['loci'][i])) for i in places]
        recs = sorted([(i, x) for i, x in recs if x], key=lambda y: -y[1]['ms'])[:2]
        if len(recs) < 2 or not all(x['introns'] for _, x in recs): continue
        gs = [gene_at(A, *d['loci'][i]) for i, _ in recs]
        if None in gs or gs[0] == gs[1]: continue
        if kind == 'decidable' and recs[0][0] != places[0]:
            recs = recs[::-1]; gs = gs[::-1]        # origin first
        key = tuple(gs) if kind == 'decidable' else tuple(sorted(gs))
        groups[kind][key].append((q, [d['loci'][i] for i, _ in recs] if kind == 'decidable' or tuple(gs) == key else [d['loci'][i] for i, _ in recs][::-1]))
    feats, exs = [], []
    # undecidable examples: pairs of protein-coding genes (rDNA and other repeat units are excluded from assignment anyway)
    pc = lambda k: all(biot.get(g) == 'protein_coding' for g in k)
    picks = [('undecidable', k, v) for k, v in sorted(((k, v) for k, v in groups['undecidable'].items() if pc(k)), key=lambda kv: -len(kv[1]))[:n_undec]] + \
            [('decidable', k, v) for k, v in sorted(groups['decidable'].items(), key=lambda kv: -len(kv[1]))[:n_dec]]
    for n, (kind, gs, members) in enumerate(picks):
        exid = f"{kind}-{n + 1}"
        members = members[:40]
        # one lane per read, the same in both panels: ordered by the read's start at copy 0 (5' to 3' on its strand)
        first = {q: rec_at(q, loci[0]) for q, loci in members}
        order = sorted((q for q, _ in members if first[q]), key=lambda q: first[q]['s'])
        order += [q for q, _ in members if q not in order]
        rank = {q: i for i, q in enumerate(order)}
        for ci, g in enumerate(gs):
            c, st = tx[gene_tx[g][0]][0], tx[gene_tx[g][0]][1]
            # window: the gene's transcripts and the reads' records here
            spans = [(min(a for a, _ in ex[t]), max(b for _, b in ex[t])) for t in gene_tx[g]]
            recs = [(q, rec_at(q, loci[ci])) for q, loci in members]
            recs = [(q, x) for q, x in recs if x]
            lo = min([s for s, _ in spans] + [x['s'] for _, x in recs]) - 200
            hi = max([e for _, e in spans] + [x['e'] for _, x in recs]) + 200
            lane = 0
            ts = sorted(gene_tx[g], key=lambda t: -sum(b - a for a, b in ex[t]))[:8]
            for t in ts:
                e = sorted(ex[t])
                feats.append(dict(example=exid, copy=ci, lane=lane, kind='tx_line', start=e[0][0], end=e[-1][1], label=t, primary=None, read=None))
                for a, b in e:
                    feats.append(dict(example=exid, copy=ci, lane=lane, kind='tx_exon', start=a, end=b, label=t, primary=None, read=None))
                lane += 1
            lane = 8   # reads start on the same lane in both panels, whatever the number of transcripts
            nk = 0; nmm = 0; nprim = 0
            for q, x in recs:
                k = known(x['introns'], g); nk += k; nmm += len(x['mm']); nprim += x['primary']
                L = lane + rank[q]
                feats.append(dict(example=exid, copy=ci, lane=L, kind='read_line', start=x['s'], end=x['e'], label=f"{'known chain' if k else 'chain not annotated'}; ms {x['ms']}, AS {x['AS']}; {len(x['mm'])} mismatches",
                                  primary=x['primary'], read=q))
                for a, b in x['blocks']:
                    feats.append(dict(example=exid, copy=ci, lane=L, kind='read_block', start=a, end=b, label=None, primary=x['primary'], read=q))
                for p in x['mm']:
                    feats.append(dict(example=exid, copy=ci, lane=L, kind='mismatch', start=p, end=p + 1, label=None, primary=x['primary'], read=q))
            gid = g
            exs.append(dict(key=f"{exid}|{ci}", example=exid, kind=kind, copy=ci, gene=name.get(gid, gid), biotype=biot.get(gid, ''),
                            chrom=c, strand=st, start=lo, end=hi, transcripts_shown=len(ts), transcripts_total=len(gene_tx[g]),
                            reads=len(recs), reads_known_chain=nk, primary_here=nprim,
                            mismatches_per_read=round(nmm / max(1, len(recs)), 2),
                            role=('origin' if ci == 0 else 'competitor') if kind == 'decidable' else 'tied copy'))
    code = open(os.path.abspath(__file__)).read()
    meta = [{'label': 'near_tie_igv.py', 'value': code, 'type': 'code', 'lang': 'python'},
            {'label': 'Input', 'value': f"{classes}; {sel} (HG002 Kinnex on CHM13, 30,000 near-tied reads); {gff}"}]
    json.dump({'rows': feats, 'meta': meta}, open(out + '.features.json', 'w'))
    json.dump({'rows': exs, 'meta': meta}, open(out + '.examples.json', 'w'))
    for e in exs:
        print(e)

if __name__ == '__main__':
    main(*sys.argv[1:])
