#!/usr/bin/env python3
"""How much of its source's transcript a member place reproduces, as a second criterion (register 1981).

    member_coverage.py MEMBERS.tsv ASSEMBLY.gtf ROWS.tsv RELATIONS.tsv

coverage(row) = max_bases / exonic length of the source locus's most-reads transcript (longest on a tie), capped at 1.
ROWS.tsv = pair_similarity.py's rows (band all) made from place_classes.py's pairs, which are the members rows in band
order: the two files are walked together. Baseline = members_grid.py's 'on, neighbours kept' at 0.3."""
import sys, re, collections
KINDS = ('paralog_recent', 'paralog_older', 'pseudogene_of', 'same_source')
BANDS = ((0.3, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.01))

def rep_len(gtf):
    """{gene: exonic length of its most-reads transcript}"""
    ln = collections.Counter(); rd = {}; gene = {}
    for l in open(gtf):
        f = l.split('\t')
        if len(f) < 9 or f[2] not in ('transcript', 'exon'): continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        if f[2] == 'exon': ln[t] += int(f[4]) - int(f[3]) + 1
        else:
            m = re.search(r' reads "([^"]+)"', f[8]); rd[t] = float(m.group(1)) if m else 0.0
        gene[t] = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
    best = {}
    for t, g in gene.items():
        k = (rd.get(t, 0.0), ln[t])
        if g not in best or k > best[g][0]: best[g] = (k, ln[t])
    return {g: v[1] for g, v in best.items()}

def coverage(max_bases, length): return min(1.0, max_bases / length) if length else 0.0

def table(rows, rel, cuts=(0.0, 0.25, 0.5, 0.75)):
    """rows [(class, place, src, under set, repeat share, similar, neighbour, coverage, novel)] -> one dict per cut"""
    den = collections.Counter(); pairs = set()
    for s in {r[2] for r in rows}:
        for g, k in rel.get(s, {}).items(): pairs.add((s, g, k)); den[k] += 1
    out = []
    for c in cuts:
        keep = [r for r in rows if (r[4] < 0.5 or r[6]) and r[7] >= c]; hit = {(r[2], g) for r in keep for g in r[3]}
        out.append(dict(cut=c, rows=len(keep), prec=sum(1 for r in keep if r[0][0] in '123' or r[5]) / max(1, len(keep)), novel=len({r[1] for r in keep if r[8]}),
                        reach={k: sum(1 for s, g, kk in pairs if kk == k and (s, g) in hit) for k in KINDS}, den=den))
    return out

def main(mem, gtf, rowsp, relp):
    L = rep_len(gtf); M = [l.rstrip('\n').split('\t') for l in open(mem)][1:]
    order = [p for lo, hi in BANDS for p in M if lo <= float(p[7]) < hi]          # place_classes.py's order of pairs
    R = [l.rstrip('\n').split('\t') for l in open(rowsp)]; i = 0; rows = []
    for p in order:                                   # pair_similarity drops no row but those without an annotated source may be absent: match by place
        if i < len(R) and R[i][2] == p[0] and abs(float(R[i][0]) - max(lo for lo, hi in BANDS if lo <= float(p[7]))) < 1e-9:
            r = R[i]; i += 1
            rows.append((r[1], r[2], r[3], set(r[4].split(',')) - {'.'}, float(r[5]), r[6] == '1', r[7:] == ['1'], coverage(int(p[9]), L.get(p[5], 0)), p[10] == '.'))
    assert i == len(R), f'walked {i} of {len(R)} scorer rows: the two files are not in step'
    rel = collections.defaultdict(dict)
    for l in list(open(relp))[1:]:
        f = l.rstrip('\n').split('\t'); rel[f[0]][f[2]] = f[4]; rel[f[2]].setdefault(f[0], f[4])
    T = table(rows, rel); b = T[0]
    print('coverage >=\trows\tat a relative\tplaces on no assembled locus\t' + '\t'.join(f'reach {k} (x baseline)' for k in KINDS))
    for t in T:
        print(t['cut'], t['rows'], f"{t['prec']:.3f}", t['novel'], *[f"{t['reach'][k]} ({t['reach'][k] / max(1, b['reach'][k]):.2f})" for k in KINDS], sep='\t')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        assert coverage(500, 1000) == 0.5 and coverage(3000, 1000) == 1.0 and coverage(5, 0) == 0.0
        rel = {'A': {'B': 'paralog_recent', 'P': 'pseudogene_of'}}
        rows = [('1 listed relative', 'p1', 'A', {'B'}, 0.0, False, False, 0.9, False), ('1 listed relative', 'p2', 'A', {'P'}, 0.0, False, False, 0.2, True),
                ('5 another gene', 'p3', 'A', {'X'}, 0.0, False, False, 0.1, True), ('4 repeat', 'p4', 'A', {'Y'}, 0.9, False, False, 1.0, True)]
        T = table(rows, rel); assert T[0]['rows'] == 3 and T[0]['reach']['pseudogene_of'] == 1 and T[0]['novel'] == 2
        assert T[2]['rows'] == 1 and T[2]['prec'] == 1.0 and T[2]['reach'] == {'paralog_recent': 1, 'paralog_older': 0, 'pseudogene_of': 0, 'same_source': 0}
        print('member_coverage selftest ok')
    else: main(*sys.argv[1:])
