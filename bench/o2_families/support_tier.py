#!/usr/bin/env python3
"""Mark each transcript by whether a PRIMARY alignment carries it (register 1952).

    support_tier.py ASSEMBLY.gtf PREFIX.provenance.transcripts.tsv OUT_PREFIX

A transcript built with near-tied secondary alignments may have no primary alignment behind it: it then shows a chain
some read fits there as well as elsewhere, not a chain a read was placed on. tools/read_provenance.py counts, per
transcript, the primary alignments whose intron chain is the transcript's (`exact`, `exact_after_replacement`) or, for a
single-exon transcript, that lie inside it (`unspliced`). Tier:
  primary     at least one such primary alignment
  tied_only   none
Writes OUT_PREFIX.gtf (every transcript, its lines tagged support "primary" | "tied_only"; nothing is dropped),
OUT_PREFIX.primary.gtf and OUT_PREFIX.tied_only.gtf (the two tiers, for scoring)."""
import sys, re, os

def tiers(table):
    """-> {transcript: 'primary' | 'tied_only'} from the provenance transcript table"""
    out = {}
    with open(table) as f:
        h = f.readline().rstrip('\n').split('\t'); i = {k: h.index(k) for k in ('transcript', 'exact', 'exact_after_replacement', 'unspliced')}
        for l in f:
            x = l.rstrip('\n').split('\t')
            out[x[i['transcript']]] = 'primary' if int(x[i['exact']]) + int(x[i['exact_after_replacement']]) + int(x[i['unspliced']]) > 0 else 'tied_only'
    return out

def run(gtf, table, out):
    t = tiers(table); n = {'primary': 0, 'tied_only': 0}
    with open(out + '.gtf.tmp', 'w') as a, open(out + '.primary.gtf.tmp', 'w') as p, open(out + '.tied_only.gtf.tmp', 'w') as s:
        for l in open(gtf):
            m = re.search(r'transcript_id "([^"]+)"', l)
            if not m or l.startswith('#'): a.write(l); continue
            k = t.get(m.group(1), 'tied_only'); f = l.split('\t')
            if len(f) > 8 and f[2] == 'transcript': n[k] += 1
            l = l.rstrip('\n').rstrip().rstrip(';') + f'; support "{k}";\n'
            a.write(l); (p if k == 'primary' else s).write(l)
    for x in ('.gtf', '.primary.gtf', '.tied_only.gtf'): os.replace(out + x + '.tmp', out + x)
    print(f"support_tier: {n['primary']} transcripts with a primary alignment, {n['tied_only']} tied-only -> {out}.gtf")
    return n

def selftest():
    import tempfile
    d = tempfile.mkdtemp(); g = os.path.join(d, 'x.gtf'); t = os.path.join(d, 't.tsv')
    L = lambda f, tid: f'c\tx\t{f}\t1\t9\t.\t+\t.\tgene_id "G"; transcript_id "{tid}";\n'
    open(g, 'w').write('# h\n' + L('transcript', 'A') + L('exon', 'A') + L('transcript', 'B') + L('exon', 'B') + L('transcript', 'M') + L('transcript', 'N'))
    open(t, 'w').write('transcript\tgene\texact\texact_after_replacement\tpartial\tunspliced\nA\tG\t0\t1\t3\t0\nB\tG\t0\t0\t7\t0\nM\tG\t0\t0\t0\t2\n')
    assert run(g, t, os.path.join(d, 'o')) == {'primary': 2, 'tied_only': 2}       # A by a replaced junction, M unspliced; B has only partial reads; N is not in the table
    o = open(os.path.join(d, 'o.gtf')).read(); assert o.count('support "primary"') == 3 and o.count('support "tied_only"') == 3 and o.startswith('# h')
    assert 'transcript_id "B"' in open(os.path.join(d, 'o.tied_only.gtf')).read() and 'transcript_id "B"' not in open(os.path.join(d, 'o.primary.gtf')).read()
    print('support_tier selftest ok')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: selftest()
    elif len(sys.argv) == 4: run(*sys.argv[1:])
    else: sys.exit(__doc__)
