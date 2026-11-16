#!/usr/bin/env python3
"""PREREG ISOSEQ-CLUSTER2 (register 1759).
  cluster2_cmp.py gtf TRANSCRIPTS.bam ALIGNED.sam OUT.gtf   primary records -> GTF, reads = FLNC count of the cluster; prints the gate
  cluster2_cmp.py matched REF.gtf OURS.gtf CLUSTER2.gtf     part 1: exact intron chains at matched count"""
import sys, re, os
if sys.argv[1] == 'gtf':
    import pysam
    tb, sam, out = sys.argv[2:5]; sup = {}
    for r in pysam.AlignmentFile(tb, check_sq=False):
        sup[r.query_name] = r.get_tag('is') if r.has_tag('is') else len(r.get_tag('im').split(',')) if r.has_tag('im') else 0
    assert sup and min(sup.values()) >= 1, 'no cluster support tag (is / im) on the transcripts'
    placed = set()
    with open(out, 'w') as w:
        for line in open(sam):
            if line.startswith('@'): continue
            f = line.split('\t', 11); flag = int(f[1])
            if flag & 2308 or f[2] == '*': continue
            cur = int(f[3]); exons = []
            for ln, op in re.findall(r'(\d+)([MIDNSHP=X])', f[5]):
                ln = int(ln)
                if op in 'M=XD':
                    if not exons or exons[-1][1] != cur: exons.append([cur, cur + ln])
                    else: exons[-1][1] = cur + ln
                    cur += ln
                elif op == 'N': cur += ln
            s = '-' if flag & 16 else '+'; t = f[0].replace('/', '_'); placed.add(f[0])
            a = f'gene_id "{t}"; transcript_id "{t}"; reads "{sup[f[0]]}";'
            w.write(f'{f[2]}\tcluster2\ttranscript\t{exons[0][0]}\t{exons[-1][1] - 1}\t.\t{s}\t.\t{a}\n')
            for x, y in exons: w.write(f'{f[2]}\tcluster2\texon\t{x}\t{y - 1}\t.\t{s}\t.\t{a}\n')
    fr = len(placed) / len(sup)
    print(f'cluster2 transcripts {len(sup)}, FLNC in clusters {sum(sup.values())}, placed {len(placed)} ({fr:.3f}); gate >= 0.90: {"PASS" if fr >= 0.9 else "INVALID"}')
else:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'famcmp'))
    from curve import load
    ref = set(load(sys.argv[2])); o = load(sys.argv[3], 'reads'); c = load(sys.argv[4], 'reads'); K = min(len(o), len(c))
    top = lambda d: {x for x, _ in sorted(d.items(), key=lambda kv: (-kv[1][0], kv[0]))[:K]}
    print('set\tchains\tcorrect\tprecision')
    for n, s in (('ours', set(o)), ('cluster2', set(c)), (f'ours top {K}', top(o)), (f'cluster2 top {K}', top(c))):
        print(f'{n}\t{len(s)}\t{len(s & ref)}\t{len(s & ref) / len(s):.3f}')
    a, b = len(top(o) & ref), len(top(c) & ref)
    print(f'part 1 at K = {K}: cluster2 {b} vs ours {a} correct ({b / a:.3f} x): cluster2 {"BETTER" if b > 1.02 * a else "WORSE" if b < 0.98 * a else "TIED"}')
    print(f'shared chains {len(set(o) & set(c))}; correct only ours {len((set(o) - set(c)) & ref)}, only cluster2 {len((set(c) - set(o)) & ref)}')
