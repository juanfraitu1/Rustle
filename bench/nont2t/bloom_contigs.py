#!/usr/bin/env python3
"""PREREG r1661 (docs/PREREG_rnabloom_absent_copies_2026-11-12.md): score reference-free contigs as absent-copy calls.
  bloom_contigs.py score GGO.bam MAT.bam PAT.bam DRAFT.bam PANEL.genes.tsv READ_ORIGIN.tsv TRUTH2.tsv OUT.tsv
  bloom_contigs.py selftest"""
import collections
import sys

ART_DE, ART_COV, TAU, CALL_COV, MIN_READS = 0.002, 0.95, 0.005, 0.90, 10


def best(recs):
    """recs: (de, covered fraction) of primary records on the truth genomes -> the best one (lowest de among full ones)."""
    return min(recs, key=lambda r: (r[1] < ART_COV, r[0])) if recs else None


def artifact(recs):
    b = best(recs)
    return b is None or b[0] > ART_DE or b[1] < ART_COV


def called(draft):
    """draft: (de, covered fraction) of the draft primary, or None when unmapped."""
    return draft is None or draft[0] > TAU or draft[1] < CALL_COV


def gene_of(blocks, contig, exons):
    """blocks: aligned (start, end) on contig; exons: {gene: (contig, [(s, e)])} -> gene holding >= 50% of the bases."""
    tot = sum(e - s for s, e in blocks)
    ov = collections.Counter()
    for g, (c, ex) in exons.items():
        if c == contig:
            ov[g] = sum(max(0, min(e, b) - max(s, a)) for s, e in blocks for a, b in ex)
    g, n = ov.most_common(1)[0] if ov else (None, 0)
    return g if tot and n >= 0.5 * tot else None


def verdict(rows, missing, present):
    """rows: (artifact, called, strict gene or None, set of all genes) per contig."""
    found, strict, alarm = set(), set(), set()
    for art, call, sg, genes in rows:
        if not call:
            continue
        if not art:
            found |= genes & missing
            if sg in missing:
                strict.add(sg)
        if sg in present:
            alarm.add(sg)
    return found, strict, alarm


def primaries(path):
    import pysam
    out = {}
    for r in pysam.AlignmentFile(path):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        out[r.query_name] = (r.get_tag('de'), r.query_alignment_length / r.infer_read_length())
    return out


def main(a):
    import pysam
    ggo, mat, pat, draft, genes, origin, truth2, outp = a
    exons = {}
    for l in open(genes):
        f = l.rstrip('\n').split('\t')
        if f[0] != 'gene':
            exons[f[0]] = (f[3], [tuple(map(int, x.split('-'))) for x in f[5].split(',')])
    n = collections.Counter(l.split('\t')[1] for l in open(origin) if not l.startswith('read\t'))
    st = {l.split('\t')[0]: l.split('\t')[3] for l in open(truth2) if not l.startswith('gene\t')}
    missing = {g for g, c in n.items() if c >= MIN_READS and st.get(g) == 'missing'}
    present = {g for g, c in n.items() if c >= MIN_READS and st.get(g) == 'present'}
    T = [primaries(p) for p in (ggo, mat, pat)]
    D = primaries(draft)
    names, prim, alt = [], {}, collections.defaultdict(list)
    for r in pysam.AlignmentFile(ggo):
        if r.query_name not in prim and r.query_name not in names:
            names.append(r.query_name)
        if r.is_unmapped or r.is_supplementary:
            continue
        g = gene_of(r.get_blocks(), r.reference_name, exons)
        if r.is_secondary:
            alt[r.query_name].append((r.get_tag('AS'), g))
        else:
            prim[r.query_name] = (r.get_tag('AS'), g)
    rows = []
    with open(outp, 'w') as o:
        o.write('contig\tartifact\ttruth_de\ttruth_cov\tcalled\tdraft_de\tdraft_cov\tgene\ttied_genes\n')
        for c in names:
            recs = [t[c] for t in T if c in t]
            b, d = best(recs), D.get(c)
            pg = prim.get(c, (0, None))
            tied = {g for s, g in alt[c] if g and s >= pg[0] and g != pg[1]}
            genes_c = ({pg[1]} if pg[1] else set()) | tied
            row = (artifact(recs), called(d), pg[1] if not tied else None, genes_c)
            rows.append(row)
            o.write('\t'.join(map(str, [c, int(row[0]), *(b or ('NA', 'NA')), int(row[1]), *(d or ('NA', 'NA')), pg[1] or '-',
                                        ','.join(sorted(tied)) or '-'])) + '\n')
    found, strict, alarm = verdict(rows, missing, present)
    art = sum(r[0] for r in rows) / max(len(rows), 1)
    calls = [r for r in rows if r[1]]
    print(f'contigs {len(rows)} artifact {sum(r[0] for r in rows)} = {art:.3f} (bar <= 0.10); called {len(calls)}, of them artifact '
          f'{sum(r[0] for r in calls)}, on a missing gene {sum(bool(r[3] & missing) for r in calls)}, on a present gene '
          f'{sum(bool(r[3] & present) for r in calls)}, no scorable panel gene {sum(not (r[3] & (missing | present)) for r in calls)}')
    print(f'missing found {len(found)} of {len(missing)} (bar >= 13); strict {len(strict)} (bar >= 9); present false alarms '
          f'{len(alarm)} of {len(present)} (bar <= 1)')
    ok = len(found) >= 13 and len(strict) >= 9 and len(alarm) <= 1 and art <= 0.10
    print('not found:', ' '.join(sorted(missing - found)) or '-', '| alarms:', ' '.join(sorted(alarm)) or '-')
    print('VERDICT', 'PASS' if ok else 'FAIL (STOP)')


def selftest():
    assert artifact([]) and artifact([(0.004, 1.0)]) and artifact([(0.0, 0.5)]) and not artifact([(0.004, 1.0), (0.001, 0.99)])
    assert called(None) and called((0.006, 1.0)) and called((0.0, 0.8)) and not called((0.004, 0.95))
    ex = {'A': ('c1', [(0, 100), (200, 300)]), 'B': ('c1', [(1000, 1100)]), 'C': ('c2', [(0, 100)])}
    assert gene_of([(10, 90), (210, 290)], 'c1', ex) == 'A' and gene_of([(500, 900)], 'c1', ex) is None
    assert gene_of([(0, 100)], 'c2', ex) == 'C' and gene_of([(90, 300)], 'c2', ex) is None
    rows = [(False, True, 'm1', {'m1'}), (False, True, None, {'m2', 'm3'}), (True, True, 'm4', {'m4'}),
            (False, False, 'm5', {'m5'}), (False, True, 'p1', {'p1'}), (False, True, None, {'p2', 'm2'})]
    f, s, al = verdict(rows, {'m1', 'm2', 'm3', 'm4', 'm5'}, {'p1', 'p2'})
    assert f == {'m1', 'm2', 'm3'} and s == {'m1'} and al == {'p1'}, (f, s, al)
    print('ok')


if __name__ == '__main__':
    selftest() if sys.argv[1] == 'selftest' else main(sys.argv[2:])
