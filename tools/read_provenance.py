#!/usr/bin/env python3
"""Where each read went: from its alignment to the transcript it supports, with what was changed on the way.

    read_provenance.py ASSEMBLY.gtf BAM OUT_PREFIX [--fasta GENOME.fa] [--as-table TAB] [--shortexons PREFIX.shortexons.reads.tsv]
                       [--sequences PREFIX.o2.sequences.tsv] [--copies PREFIX.fam.copies.tsv] [--window 10] [--threads N]

RECONSTRUCTED after the run from the alignments and the assembly; the assembler itself keeps no per-read record. What is
read off exactly and what is inferred is said per column. One row per PRIMARY alignment.

OUT_PREFIX.reads.tsv.gz
  read, contig, start, end, strand, introns          the alignment (1-based closed span; strand of the transcript)
  relation       exact                    its intron chain is a transcript's chain                                   [exact]
                 exact_after_replacement  after replacing junctions (below) its chain is a transcript's              [inferred]
                 partial                  its chain is a contiguous piece of exactly one transcript (a shorter read) [exact]
                 partial_after_replacement, partial_ambiguous (a piece of several transcripts: listed, up to 5)
                 unspliced_model          no intron; >= 80% of it inside a single-exon transcript                    [exact]
                 unspliced_in_exon        no intron; >= 80% of it inside one exon of a spliced transcript (not counted as support)
                 none                     no transcript carries its chain
  transcript, gene                        the transcript(s) and their gene
  replaced       junctions of the read that are in no transcript but lie within --window bases of one that is, old>new.
                 The assembler snaps and replaces weak junctions by a supported neighbour; this column infers that from
                 the outcome, it does not read the assembler's decision                                               [inferred]
  shortexon_rewrite   yes when the short-exon repair rewrote this read's alignment before assembly (--shortexons)     [exact]
  rival_within_10ms   yes when another alignment of the read is within 10 ms of its best (--as-table): a near-tied multimapper [exact]
  label          o2_status of the gene (expressed / shared / homology_only), when the assembly carries it            [exact]
  sequence       for a tied read, the sequence it supports in PREFIX.o2.sequences.tsv: locus read at / number (--sequences) [exact]
  reason         for relation none, why no transcript carries the chain:                                              [inferred]
                 non-canonical junction   an intron that is not GT-AG, GC-AG or AT-AC (the assembler's strict gate; --fasta)
                 under the read floor     fewer than 2 distinct alignments (start, end) share this chain on the contig
                 removed by the rule R    with --polish-reasons (the assembler's own table, PREFIX.polish_reasons.tsv): the chain
                                          was a transcript when the rules ran and R removed it                            [exact]
                 removed by the rule R (polish_rules.py)   with --unpolished and --reasons: the chain was a transcript before polishing
                                          and tools/polish_rules.py names the rule that removed it (ism, frac, ret, ...)  [exact]
                 removed by a later rule  the chain had N distinct alignments and passed the floor: isoform fraction,
                                          retained intron, spliced share or the final filter removed it (not resolved)
                 unspliced, no model here
  A tied read is not given a copy of origin: `sequence` and `label` are all that can be said (docs/O2_FAMILY_BENCHMARK.md).
OUT_PREFIX.transcripts.tsv   one row per transcript: the reads behind it by relation, the reads the assembler counted
  (attribute `reads`), reads rewritten by the short-exon repair, near-tied reads, the three most frequent junction replacements.
OUT_PREFIX.family_reads.tsv.gz (with --copies)   the rows of reads on a gene of a multi-copy family, or with a rival within 10 ms."""
import sys, os, re, gzip, argparse, collections, multiprocessing, bisect
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import eq_classes as E
import tie_table as TT

COLS = 'read\tcontig\tstart\tend\tstrand\tintrons\trelation\ttranscript\tgene\treplaced\tshortexon_rewrite\trival_within_10ms\tlabel\tsequence\treason'
CANON = {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')}
COMP = str.maketrans('ACGT', 'TGCA')
G = {}

def load_gtf(path):
    tx = {}; ex = collections.defaultdict(list)
    for l in open(path):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] not in ('transcript', 'exon'): continue
        t = re.search(r'transcript_id "([^"]+)"', f[8])
        if not t: continue
        t = t.group(1)
        if f[2] == 'exon': ex[t].append((int(f[3]), int(f[4])))
        else:
            g = re.search(r'gene_id "([^"]+)"', f[8]).group(1); rd = re.search(r'reads "([^"]+)"', f[8]); st = re.search(r'o2_status "([^"]+)"', f[8])
            tx[t] = [f[0], f[6], g, rd.group(1) if rd else '.', st.group(1) if st else '.']
    out = collections.defaultdict(list)
    for t, (c, s, g, rd, lab) in tx.items():
        e = sorted(ex.get(t, []))
        if not e: continue
        out[c].append((t, s, g, rd, lab, e, tuple((a[1], b[0] - 1) for a, b in zip(e, e[1:]))))
    return out

def introns(r):
    out = []; pos = r.reference_start
    for op, n in r.cigartuples:
        if op == 3: out.append((pos, pos + n)); pos += n
        elif op in (0, 2, 7, 8): pos += n
    return tuple(out)

def spans(v):
    """(start, end, id) sorted, with the starts and the longest span: so `inside` starts at the first span that can reach"""
    v = sorted(v)
    return v, [x[0] for x in v], max((x[1] - x[0] for x in v), default=0)

def contig_index(txs):
    """per strand: chain -> transcript, sorted junctions, junction -> transcripts, single-exon spans, exons of spliced transcripts"""
    ix = {}
    for s in '+-.':
        chain = {}; junc = set(); has = collections.defaultdict(set); mono = []; exons = []; info = {}
        for t, st, g, rd, lab, e, ch in txs:
            if st != s: continue
            info[t] = (g, lab, ch)
            if ch:
                chain.setdefault(ch, t); junc.update(ch)
                for j in ch: has[j].add(t)
                exons += [(a, b, t) for a, b in e]
            else: mono.append((e[0][0], e[0][1], t))
        js = sorted(junc)
        ix[s] = dict(chain=chain, junc=junc, js=js, jstart=[j[0] for j in js], has=has, mono=spans(mono), exons=spans(exons), info=info)
    return ix

def near_junction(j, x, w):
    i = bisect.bisect_left(x['jstart'], j[0] - w); best = None
    while i < len(x['js']) and x['js'][i][0] <= j[0] + w:
        k = x['js'][i]
        if abs(k[1] - j[1]) <= w:
            d = abs(k[0] - j[0]) + abs(k[1] - j[1])
            if best is None or d < best[0]: best = (d, k)
        i += 1
    return best[1] if best else None

def pieces(ch, x):
    """transcripts whose chain holds ch as a contiguous piece"""
    cand = x['has'].get(ch[0], ())
    out = []
    for t in cand:
        full = x['info'][t][2]; n = len(ch)
        if len(full) > n and any(full[i:i + n] == ch for i in range(len(full) - n + 1)): out.append(t)
    return sorted(out)

def inside(a, b, spans):
    """the span (start, end, id) that covers most of [a, b] (1-based closed), with the covered share"""
    spans, starts, longest = spans
    best = (0.0, None); L = b - a + 1
    for s, e, t in spans[bisect.bisect_left(starts, a - longest):]:
        if s > b: break
        if e < a: continue
        ov = min(b, e) - max(a, s) + 1
        if ov / L > best[0]: best = (ov / L, t)
    return best

def work(contig):
    import pysam
    txs = G['tx'].get(contig, []); ix = contig_index(txs); w = G['window']
    B = pysam.AlignmentFile(G['bam']); F = pysam.FastaFile(G['fasta']) if G['fasta'] else None
    rows = []; ends = collections.defaultdict(set); canon = {}
    for r in B.fetch(contig):
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        s = E.tx_strand(r); x = ix[s]; ch = introns(r); a, b = r.reference_start + 1, r.reference_end
        rel = 'none'; tids = []; rep = []; key = None
        if not ch:
            sh, t = inside(a, b, x['mono'])
            if t is not None and sh >= 0.8: rel = 'unspliced_model'; tids = [t]
            else:
                sh, t = inside(a, b, x['exons'])
                if t is not None and sh >= 0.8: rel = 'unspliced_in_exon'; tids = [t]
        else:
            cc = []
            for j in ch:
                if j in x['junc']: cc.append(j); continue
                k = near_junction(j, x, w)
                if k is None: cc.append(j)
                else: cc.append(k); rep.append((j, k))
            cc = tuple(cc); key = (s, cc); ends[key].add((a, b)); sfx = '_after_replacement' if rep else ''
            if cc in x['chain']: rel = 'exact' + sfx; tids = [x['chain'][cc]]
            else:
                p = pieces(cc, x)
                if len(p) == 1: rel = 'partial' + sfx; tids = p
                elif p: rel = 'partial_ambiguous'; tids = p[:5]
                elif F is not None and key not in canon:
                    ok = True
                    for i0, i1 in cc:
                        d = F.fetch(contig, i0, i0 + 2).upper(); ac = F.fetch(contig, i1 - 2, i1).upper()
                        if s == '-': d, ac = ac[::-1].translate(COMP), d[::-1].translate(COMP)
                        if (d, ac) not in CANON: ok = False; break
                    canon[key] = ok
        rows.append((r.query_name, a, b, s, len(ch), rel, tids, rep, key))
    part = f"{G['out']}.part.{G['order'][contig]:06d}.gz"; fam = G['family']; tie = G['ties']; se = G['shortexons']; sq = G['sequences']
    rl = G.get('rules', {}).get(contig, {})
    xl = G.get('exact', {}).get(contig, {})
    agg = collections.defaultdict(lambda: collections.Counter()); reps = collections.defaultdict(collections.Counter)
    with gzip.open(part, 'wt', compresslevel=3) as o, (gzip.open(part + '.fam', 'wt', compresslevel=3) if fam is not None else open(os.devnull, 'w')) as of:
        for q, a, b, s, n, rel, tids, rep, key in rows:
            info = ix[s]['info']; genes = sorted({info[t][0] for t in tids}); lab = info[tids[0]][1] if tids else '.'
            reason = '.'
            if rel == 'none':
                if key is None: reason = 'unspliced, no model here'
                elif canon.get(key) is False: reason = 'non-canonical junction'
                elif key in xl: reason = 'removed by the rule ' + xl[key]
                elif key in rl: reason = 'removed by the rule ' + rl[key] + ' (polish_rules.py)'
                elif len(ends[key]) < 2: reason = 'under the read floor (1 distinct alignment)'
                else: reason = f'removed by a later rule ({len(ends[key])} distinct alignments)'
            rv = '.' if tie is None else 'yes' if q in tie else 'no'
            line = '\t'.join((q, contig, str(a), str(b), s, str(n), rel, ','.join(tids) or '.', ','.join(genes) or '.',
                              ';'.join(f'{j[0] + 1}-{j[1]}>{k[0] + 1}-{k[1]}' for j, k in rep) or '.', 'yes' if q in se else 'no' if G['has_se'] else '.',
                              rv, lab, sq.get(q, '.'), reason)) + '\n'
            o.write(line)
            if fam is not None and (rv == 'yes' or any(g in fam for g in genes)): of.write(line)
            for t in tids:
                c = agg[t]; c[rel] += 1; c['shortexon'] += q in se; c['rival'] += q in tie if tie is not None else 0
                if rel != 'partial_ambiguous':
                    for j, k in rep: reps[t][f'{j[0] + 1}-{j[1]}>{k[0] + 1}-{k[1]}'] += 1
    trows = []
    for t, st, g, rd, lab, e, ch in txs:
        c = agg[t]
        trows.append('\t'.join(map(str, (t, g, contig, e[0][0], e[-1][1], st, len(ch), rd, c['exact'], c['exact_after_replacement'], c['partial'] + c['partial_after_replacement'],
                                          c['partial_ambiguous'], c['unspliced_model'], c['shortexon'], c['rival'],
                                          ';'.join(f'{k} x{v}' for k, v in reps[t].most_common(3)) or '.', lab))) + '\n')
    return contig, part, trows, collections.Counter(r[5] for r in rows)

def run(a):
    import pysam
    G.update(tx=load_gtf(a.gtf), bam=a.bam, fasta=a.fasta, window=a.window, out=a.out)
    refs = list(pysam.AlignmentFile(a.bam).references); G['order'] = {c: i for i, c in enumerate(refs)}
    if a.as_table: TT.check_table(a.as_table, a.bam)      # another BAM's scores would mark the wrong reads as near-tied, silently
    G['ties'] = E.load_ties(a.as_table, E.floor_fn(10), 'ms') if a.as_table else None
    se = set()
    if a.shortexons and os.path.exists(a.shortexons):
        with open(a.shortexons) as f:
            f.readline()
            for l in f: se.add(l.split('\t', 1)[0])
    G['shortexons'] = se; G['has_se'] = bool(a.shortexons)
    sq = {}
    if a.sequences and os.path.exists(a.sequences):
        with open(a.sequences) as f:
            h = f.readline().rstrip('\n').split('\t'); ia, ic, ir = h.index('read_at'), h.index('candidate'), h.index('read_names')
            for l in f:
                x = l.rstrip('\n').split('\t')
                for q in x[ir].split(','): sq[q] = f'{x[ia]}/{x[ic]}'
    G['sequences'] = sq
    fam = None
    if a.copies:
        fam = set()
        with open(a.copies) as f:
            h = f.readline().rstrip('\n').split('\t'); i = h.index('gene_id')
            for l in f: fam.add(l.split('\t')[i])
    G['family'] = fam
    # exact rule per removed chain (polish_rules.py on the unpolished assembly of the same run): {contig: {(strand, chain): rule}}
    rules = collections.defaultdict(dict)
    if a.unpolished and a.reasons:
        why = {}
        with open(a.reasons) as f:
            f.readline()
            for l in f:
                x = l.rstrip('\n').split('\t'); why[x[0]] = x[2]
        for c, txs in load_gtf(a.unpolished).items():
            for t, s_, g, rd, lab, e, ch in txs:
                if ch and t in why: rules[c].setdefault((s_, ch), why[t])
    exact = collections.defaultdict(dict)          # the assembler's own table: {contig: {(strand, chain): rule}}
    if a.polish_reasons and os.path.exists(a.polish_reasons):
        with open(a.polish_reasons) as f:
            f.readline()
            for l in f:
                x = l.rstrip('\n').split('\t')
                ex = [tuple(map(int, e.split('-'))) for e in x[3].split(',')]
                if len(ex) > 1: exact[x[1]].setdefault((x[2], tuple((p[1], q[0] - 1) for p, q in zip(ex, ex[1:]))), x[5])
    G['exact'] = exact
    G['rules'] = rules
    with multiprocessing.get_context('fork').Pool(a.threads) as pool: res = pool.map(work, refs, chunksize=1)
    tot = collections.Counter()
    with open(a.out + '.reads.tsv.gz.tmp', 'wb') as o, open(a.out + '.transcripts.tsv.tmp', 'w') as t, (open(a.out + '.family_reads.tsv.gz.tmp', 'wb') if fam is not None else open(os.devnull, 'wb')) as of:
        o.write(gzip.compress((COLS + '\n').encode())); of.write(gzip.compress((COLS + '\n').encode()))
        t.write('transcript\tgene\tcontig\tstart\tend\tstrand\tintrons\treads_in_gtf\texact\texact_after_replacement\tpartial\tpartial_ambiguous\tunspliced\tshortexon_rewritten\trival_within_10ms\treplaced_junctions\tlabel\n')
        for contig, part, trows, cn in res:
            with open(part, 'rb') as p: o.write(p.read())
            os.remove(part)
            if fam is not None:
                with open(part + '.fam', 'rb') as p: of.write(p.read())
                os.remove(part + '.fam')
            t.writelines(trows); tot.update(cn)
    os.replace(a.out + '.reads.tsv.gz.tmp', a.out + '.reads.tsv.gz'); os.replace(a.out + '.transcripts.tsv.tmp', a.out + '.transcripts.tsv')
    if fam is not None: os.replace(a.out + '.family_reads.tsv.gz.tmp', a.out + '.family_reads.tsv.gz')
    n = sum(tot.values())
    print('read_provenance: ' + f'{n} reads: ' + ', '.join(f'{k} {v} ({v / max(1, n):.3f})' for k, v in tot.most_common()) + f' -> {a.out}.reads.tsv.gz, {a.out}.transcripts.tsv')
    return tot

def selftest():
    import tempfile, pysam
    d = tempfile.mkdtemp(); gtf = os.path.join(d, 'x.gtf'); bam = os.path.join(d, 'x.bam')
    def tr(t, g, ex, extra=''):
        s = f'c\tx\ttranscript\t{ex[0][0]}\t{ex[-1][1]}\t.\t+\t.\tgene_id "{g}"; transcript_id "{t}"; reads "4";{extra}\n'
        return s + ''.join(f'c\tx\texon\t{a}\t{b}\t.\t+\t.\tgene_id "{g}"; transcript_id "{t}";\n' for a, b in ex)
    open(gtf, 'w').write(tr('T1', 'G1', [(101, 200), (301, 400), (501, 600)], ' o2_status "expressed";') + tr('T2', 'G1', [(101, 200), (501, 600)]) + tr('M1', 'G2', [(2001, 2500)]))
    hdr = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'c', 'LN': 5000}]}
    def rd(name, pos, cigar):
        r = pysam.AlignedSegment(); r.query_name = name; r.reference_id = 0; r.reference_start = pos; r.cigarstring = cigar; r.flag = 0; r.mapping_quality = 60
        n = sum(int(x) for x, o in re.findall(r'(\d+)([MIS=X])', cigar)); r.query_sequence = 'A' * n; r.query_qualities = None
        return r
    recs = [rd('exact', 100, '100M100N100M100N100M'), rd('shift', 100, '103M97N100M100N100M'), rd('piece', 320, '80M100N100M'),
            rd('skip', 100, '100M300N100M'), rd('lone', 100, '100M150N250M'), rd('mono', 2050, '300M'), rd('inexon', 310, '80M'), rd('nowhere', 3000, '200M')]
    with pysam.AlignmentFile(bam, 'wb', header=hdr) as w:
        for r in sorted(recs, key=lambda r: r.reference_start): w.write(r)
    pysam.index(bam)
    pr = os.path.join(d, 'pr.tsv'); open(pr, 'w').write('transcript\tcontig\tstrand\texons\treads\trule\nX1\tc\t+\t101-200,351-600\t3\tfrac\n')   # the chain of `lone`
    a = argparse.Namespace(gtf=gtf, bam=bam, out=os.path.join(d, 'o'), fasta=None, as_table=None, shortexons=None, sequences=None, copies=None, window=10, threads=1, unpolished=None, reasons=None, polish_reasons=None)
    run(a)
    rows = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in gzip.open(a.out + '.reads.tsv.gz', 'rt')}
    want = {'exact': ('exact', 'T1'), 'shift': ('exact_after_replacement', 'T1'), 'piece': ('partial', 'T1'), 'skip': ('exact', 'T2'), 'lone': ('none', '.'),
            'mono': ('unspliced_model', 'M1'), 'inexon': ('unspliced_in_exon', 'T1'), 'nowhere': ('none', '.')}
    for k, (rel, t) in want.items(): assert rows[k][6] == rel and rows[k][7] == t, (k, rows[k])
    assert rows['shift'][9] == '204-300>201-300' and rows['lone'][14].startswith('under the read floor') and rows['nowhere'][14] == 'unspliced, no model here' and rows['exact'][12] == 'expressed'
    t = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open(a.out + '.transcripts.tsv')}
    assert t['T1'][8:13] == ['1', '1', '1', '0', '0'] and t['T1'][15] == '204-300>201-300 x1' and t['M1'][12] == '1', t['T1']
    assert inside(5, 9, spans([(1, 100, 'long'), (50, 60, 'x'), (70, 80, 'y')])) == (1.0, 'long')        # a long span that starts far to the left is still found
    a.polish_reasons = pr; run(a)
    assert next(l for l in gzip.open(a.out + '.reads.tsv.gz', 'rt') if l.startswith('lone\t')).rstrip('\n').split('\t')[14] == 'removed by the rule frac'
    print('read_provenance selftest ok')

def main():
    if sys.argv[1:] == ['--selftest']: return selftest()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('gtf'); ap.add_argument('bam'); ap.add_argument('out')
    ap.add_argument('--fasta'); ap.add_argument('--as-table'); ap.add_argument('--shortexons'); ap.add_argument('--sequences'); ap.add_argument('--copies')
    ap.add_argument('--window', type=int, default=10); ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--unpolished', help='the same run assembled with RUSTLE_ASSEMBLY_POLISH=none'); ap.add_argument('--reasons', help='polish_rules.py reasons table of --unpolished')
    ap.add_argument('--polish-reasons', help='PREFIX.polish_reasons.tsv, written by the assembler: the exact rule per removed transcript')
    run(ap.parse_args())

if __name__ == '__main__':
    main()
