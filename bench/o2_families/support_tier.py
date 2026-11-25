#!/usr/bin/env python3
"""Mark each transcript by whether a PRIMARY alignment carries it (register 1952).

    support_tier.py ASSEMBLY.gtf PREFIX.provenance.transcripts.tsv OUT_PREFIX [--reads PREFIX.provenance.reads.tsv.gz]

A transcript built with near-tied secondary alignments may have no primary alignment behind it: it then shows a chain
some read fits there as well as elsewhere, not a chain a read was placed on. tools/read_provenance.py counts, per
transcript, the primary alignments whose intron chain is the transcript's (`exact`, `exact_after_replacement`) or, for a
single-exon transcript, that lie inside it (`unspliced`). Tier:
  primary     at least one such primary alignment
  tied_only   none
Three LAYERS (2026-11-25, user: one product for primaries only and a layer that shows what the near ties add) split the
primary tier by the number of primary alignments that carry the chain, N (attribute primary_reads; the assembler's read
floor is 2):
  layer "primary"       N >= 2   the transcript stands without any tied read
  layer "tied_support"  N = 1    a primary alignment carries it and tied reads lifted it over the floor
  layer "tied_only"     N = 0
(N counts alignments, the floor counts distinct start-end pairs: a transcript with two coordinate-duplicate primaries is
called `primary` here.)
Writes OUT_PREFIX.gtf (every transcript, tagged support, primary_reads and layer; nothing is dropped),
OUT_PREFIX.primary.gtf (layer primary), OUT_PREFIX.tied.gtf (layers tied_support and tied_only),
OUT_PREFIX.tied_only.gtf and OUT_PREFIX.primary_tier.gtf (the two tiers of register 1952, for scoring), and
OUT_PREFIX.layers.tsv: per layer the transcripts, loci and o2_status labels; with --reads, of the reads that have a
rival within 10 ms (near-tied multimappers) how many a transcript of each layer carries and how many no transcript
carries, by reason: what the tied reads explain and what they leave unexplained."""
import sys, re, os, gzip, collections

LAYERS = ('primary', 'tied_support', 'tied_only')
def layer_of(n): return 'primary' if n >= 2 else 'tied_support' if n == 1 else 'tied_only'

def counts(table):
    """-> {transcript: primary alignments that carry it} from the provenance transcript table"""
    out = {}
    with open(table) as f:
        h = f.readline().rstrip('\n').split('\t'); i = {k: h.index(k) for k in ('transcript', 'exact', 'exact_after_replacement', 'unspliced')}
        for l in f:
            x = l.rstrip('\n').split('\t')
            out[x[i['transcript']]] = int(x[i['exact']]) + int(x[i['exact_after_replacement']]) + int(x[i['unspliced']])
    return out

def tied_reads(reads, layer):
    """the near-tied reads of the provenance read table by what carries them: -> Counter"""
    c = collections.Counter()
    with gzip.open(reads, 'rt') as f:
        h = f.readline().rstrip('\n').split('\t'); i = {k: h.index(k) for k in ('relation', 'transcript', 'rival_within_10ms', 'reason')}
        for l in f:
            x = l.rstrip('\n').split('\t')
            if x[i['rival_within_10ms']] != 'yes': continue
            rel = x[i['relation']]
            if rel == 'none': c['carried by no transcript: ' + re.sub(r' \(\d+ distinct alignments?\)', '', x[i['reason']])] += 1
            elif rel == 'unspliced_in_exon': c['unspliced, inside an exon of a spliced transcript (not counted as support)'] += 1
            else: c['carried by a transcript of layer ' + layer.get(x[i['transcript']].split(',')[0], 'tied_only')] += 1
    return c

def run(gtf, table, out, reads=None):
    N = counts(table); n = {'primary': 0, 'tied_only': 0}; layer = {}
    per = {k: [0, set(), collections.Counter()] for k in LAYERS}       # transcripts, loci, labels
    names = ('.gtf', '.primary.gtf', '.tied.gtf', '.primary_tier.gtf', '.tied_only.gtf')
    fh = {x: open(out + x + '.tmp', 'w') for x in names}
    for l in open(gtf):
        m = re.search(r'transcript_id "([^"]+)"', l)
        if not m or l.startswith('#'): fh['.gtf'].write(l); continue
        k = N.get(m.group(1), 0); L = layer_of(k); tier = 'primary' if k else 'tied_only'; layer[m.group(1)] = L; f = l.split('\t')
        if len(f) > 8 and f[2] == 'transcript':
            n[tier] += 1; per[L][0] += 1; g = re.search(r'gene_id "([^"]+)"', l); st = re.search(r'o2_status "([^"]+)"', l)
            if g: per[L][1].add(g.group(1))
            per[L][2][st.group(1) if st else 'unlabelled'] += 1
        l = l.rstrip('\n').rstrip().rstrip(';') + f'; support "{tier}"; primary_reads "{k}"; layer "{L}";\n'
        fh['.gtf'].write(l); fh['.primary.gtf' if L == 'primary' else '.tied.gtf'].write(l); fh['.primary_tier.gtf' if k else '.tied_only.gtf'].write(l)
    for x in names: fh[x].close(); os.replace(out + x + '.tmp', out + x)
    with open(out + '.layers.tsv.tmp', 'w') as w:
        w.write('section\titem\tvalue\n')
        for L in LAYERS:
            w.write(f'layer {L}\ttranscripts\t{per[L][0]}\nlayer {L}\tloci\t{len(per[L][1])}\n')
            for s_, v in sorted(per[L][2].items()): w.write(f'layer {L}\ttranscripts labelled {s_}\t{v}\n')
        only = lambda L: len(per[L][1] - set().union(*(per[x][1] for x in LAYERS if x != L)))
        w.write(f"loci\twith a transcript of layer primary\t{len(per['primary'][1])}\nloci\tadded by tied_support (no primary-layer transcript)\t{len(per['tied_support'][1] - per['primary'][1])}\nloci\tadded by tied_only alone\t{only('tied_only')}\n")
        if reads:
            c = tied_reads(reads, layer); tot = sum(c.values()); w.write(f'near-tied reads\ttotal (primary alignments with a rival within 10 ms)\t{tot}\n')
            for k, v in sorted(c.items(), key=lambda kv: (not kv[0].startswith('carried by a'), -kv[1])): w.write(f'near-tied reads\t{k}\t{v}\n')
    os.replace(out + '.layers.tsv.tmp', out + '.layers.tsv')
    print(f"support_tier: {n['primary']} transcripts with a primary alignment ({per['primary'][0]} stand on primaries alone, {per['tied_support'][0]} with tied support), {n['tied_only']} tied-only -> {out}.gtf")
    return n

def selftest():
    import tempfile
    d = tempfile.mkdtemp(); g = os.path.join(d, 'x.gtf'); t = os.path.join(d, 't.tsv')
    L = lambda f, tid: f'c\tx\t{f}\t1\t9\t.\t+\t.\tgene_id "G"; transcript_id "{tid}";\n'
    open(g, 'w').write('# h\n' + L('transcript', 'A') + L('exon', 'A') + L('transcript', 'B') + L('exon', 'B') + L('transcript', 'M') + L('transcript', 'N'))
    open(t, 'w').write('transcript\tgene\texact\texact_after_replacement\tpartial\tunspliced\nA\tG\t0\t1\t3\t0\nB\tG\t0\t0\t7\t0\nM\tG\t0\t0\t0\t2\n')
    open(t, 'a').write('P\tG\t3\t0\t0\t0\n'); open(g, 'a').write(L('transcript', 'P').replace('transcript_id "P";', 'transcript_id "P"; o2_status "shared";'))
    rd = os.path.join(d, 'r.tsv.gz')
    with gzip.open(rd, 'wt') as w: w.write('read\trelation\ttranscript\trival_within_10ms\treason\nr1\texact\tP\tyes\t.\nr2\texact\tA\tyes\t.\nr3\tnone\t.\tyes\tunder the read floor (1 distinct alignment)\nr4\texact\tP\tno\t.\nr5\tpartial\tB,P\tyes\t.\n')
    assert run(g, t, os.path.join(d, 'o'), rd) == {'primary': 3, 'tied_only': 2}       # A by a replaced junction, M unspliced; B has only partial reads; N is not in the table
    o = open(os.path.join(d, 'o.gtf')).read(); assert o.count('support "primary"') == 4 and o.count('support "tied_only"') == 3 and o.startswith('# h')
    assert 'transcript_id "B"' in open(os.path.join(d, 'o.tied_only.gtf')).read() and 'transcript_id "B"' not in open(os.path.join(d, 'o.primary_tier.gtf')).read()
    assert 'layer "tied_support"' in [l for l in o.splitlines() if 'transcript_id "A"' in l][0] and 'transcript_id "A"' in open(os.path.join(d, 'o.tied.gtf')).read() and open(os.path.join(d, 'o.primary.gtf')).read().count('transcript_id') == 2      # M (two unspliced primaries) and P
    T = dict((l.split('\t')[0] + '|' + l.split('\t')[1], l.rstrip('\n').split('\t')[2]) for l in open(os.path.join(d, 'o.layers.tsv')))
    assert T['layer primary|transcripts'] == '2' and T['layer primary|transcripts labelled shared'] == '1' and T['near-tied reads|total (primary alignments with a rival within 10 ms)'] == '4'
    assert T['near-tied reads|carried by a transcript of layer primary'] == '1' and T['near-tied reads|carried by a transcript of layer tied_support'] == '1' and T['near-tied reads|carried by a transcript of layer tied_only'] == '1' and T['near-tied reads|carried by no transcript: under the read floor'] == '1', T
    print('support_tier selftest ok')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: selftest()
    elif len(sys.argv) == 4: run(*sys.argv[1:])
    elif len(sys.argv) == 6 and sys.argv[4] == '--reads': run(*sys.argv[1:4], sys.argv[5])
    else: sys.exit(__doc__)
