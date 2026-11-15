#!/usr/bin/env python3
"""PREREG row 1700 Amendment 1: rebuild the SPDYE1 hidden-copy read label.

  spdye1_label.py ins   QUERY_VS_CHM13.paf QUERY_START COPIES.tsv            -> donor insertion, CHM13 anchor, family
  spdye1_label.py reads DIPLOID.sam DONOR_CHROM INS_START INS_END CHM13.bam COPIES.tsv FAMILY OUT.tsv
  spdye1_label.py --selftest"""
import sys, re, csv, collections, statistics

def insertions(cigar, qstart, tstart, minlen):
    """[(query offset, length, target pos)] of I runs >= minlen; offsets run along the alignment from qstart."""
    q, t, out = qstart, tstart, []
    for n, op in re.findall(r'(\d+)([MIDNSH=X])', cigar):
        n = int(n)
        if op == 'I' and n >= minlen: out.append((q, n, t))
        if op in 'MI=X': q += n
        if op in 'MDN=X': t += n
    return out

def bases_in(start, cigartuples, lo, hi):
    """aligned (M/=/X) reference bases of a record inside [lo, hi)."""
    p, n = start, 0
    for op, ln in cigartuples:
        if op in (0, 7, 8): n += max(0, min(p + ln, hi) - max(p, lo))
        if op in (0, 2, 3, 7, 8): p += ln
    return n

def copies(path, fam=None):
    return [(r['family_id'], r['copy_idx'], r['chrom'], int(r['start']), int(r['end'])) for r in csv.DictReader(open(path), delimiter='\t') if fam in (None, r['family_id'])]

if sys.argv[1] == '--selftest':
    assert insertions('100=5000I200=', 10, 1000, 3000) == [(110, 5000, 1100)]
    assert insertions('100=50I200=', 0, 0, 3000) == []
    assert bases_in(100, [(7, 50), (3, 1000), (7, 50)], 120, 1170) == 30 + 20
    print('selftest ok')
elif sys.argv[1] == 'ins':
    paf, q0, cp = sys.argv[2], int(sys.argv[3]), sys.argv[4]; best = None
    for l in open(paf):
        f = l.rstrip('\n').split('\t'); cg = next(x[5:] for x in f[12:] if x.startswith('cg:Z:'))
        if best is None or int(f[10]) > int(best[0][10]): best = (f, cg)
    f, cg = best; ql, qs, qe, ts, te = int(f[1]), int(f[2]), int(f[3]), int(f[7]), int(f[8])
    spans = qs < 15000 and qe > ql - 15000; ins = insertions(cg, 0, ts, 3000)
    print(f'alignment: query {qs}-{qe} of {ql} strand {f[4]} -> {f[5]}:{ts}-{te}; spans both flanks: {spans}; insertions >= 3 kb: {len(ins)}')
    if not spans or not ins: sys.exit('GATE FAIL: no flank-to-flank alignment with an insertion >= 3 kb')
    off, n, t = max(ins, key=lambda x: x[1])
    a = qs + off if f[4] == '+' else qe - off - n   # along-alignment offset -> query coordinate
    near = min((c for c in copies(cp) if c[2] == f[5]), key=lambda c: max(0, c[3] - t, t - c[4]))
    print(f'RESULT\t{q0 + a}\t{q0 + a + n}\t{f[5]}\t{t}\t{near[0]}\t{near[1]}\t{max(0, near[3] - t, t - near[4])}')
else:
    import pysam
    sam, chrom, lo, hi, bam, cp, fam, out = sys.argv[2:10]; lo, hi = int(lo), int(hi)
    lab, nmd, n = {}, {}, 0
    for r in pysam.AlignmentFile(sam):
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        n += 1
        if r.reference_name == chrom and bases_in(r.reference_start, r.cigartuples, lo, hi) >= 50: lab[r.query_name] = r.get_tag('NM')
    C = copies(cp); b = pysam.AlignmentFile(bam); rows = {}
    for c in [c for c in C if c[0] == fam]:
        for r in b.fetch(c[2], c[3], c[4]):
            if r.query_name in lab and not r.is_secondary and not r.is_supplementary:
                on = [x for x in C if x[2] == r.reference_name and r.reference_start < x[4] and x[3] < r.reference_end]
                on.sort(key=lambda x: -(min(x[4], r.reference_end) - max(x[3], r.reference_start)))
                rows[r.query_name] = (f'{on[0][0]}:{on[0][1]}' if on else 'none', r.mapping_quality, r.get_tag('NM'))
    with open(out, 'w') as w:
        w.write('read\tchm13_copy\tmapq\tnm_chm13\tnm_donor\n')
        for q in sorted(lab):
            c, mq, nm = rows.get(q, ('primary_off_family', 'NA', 'NA')); w.write(f'{q}\t{c}\t{mq}\t{nm}\t{lab[q]}\n')
    cnt = collections.Counter(v[0] for v in rows.values()); top = cnt.most_common(1)[0] if cnt else ('none', 0)
    print(f'reads aligned {n}; labelled {len(lab)}; CHM13 primaries by copy {dict(cnt.most_common(5))}; off family {len(lab) - len(rows)}')
    if rows:
        fm, ci = top[0].split(':'); c = next(x for x in C if x[0] == fm and x[1] == ci)
        pile = sum(1 for r in b.fetch(c[2], c[3], c[4]) if not r.is_secondary and not r.is_supplementary)
        print(f'top copy {top[0]} {c[2]}:{c[3]}-{c[4]} holds {top[1]} of {len(lab)} ({top[1] / len(lab):.3f}); pile {pile} primaries, label share {top[1] / pile:.3f}; '
              f'MAPQ >= 1: {sum(1 for v in rows.values() if v[1] >= 1)}; median NM CHM13 {statistics.median(v[2] for v in rows.values())} donor {statistics.median(lab.values())}')
        print('LABEL', 'SAME' if 80 <= len(lab) <= 120 and top[1] >= 0.9 * len(lab) else 'DIFFERENT')
