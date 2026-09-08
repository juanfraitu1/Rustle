#!/usr/bin/env python3
"""Attach a COPY to every assembled isoform, from read evidence rather than from position.

`copy_assign --gtf` emits de-novo isoforms whose `gene_id` is `DN_<contig>_<pos>_<n>`, while catalog copies
carry `MCL_<contig>_<pos>` tids. The GTF's own `family_id`/`copy_index` tagging joins those two by STRING and
therefore never fires when a catalog is supplied — every transcript comes out `multicopy "false"`. This script
does the join the advisor actually asks for: for each isoform, which copy do the reads supporting it belong to,
according to O2's certificate?

A read supports an isoform when their intron chains are identical (the same collapse rule that built it).
Unspliced isoforms take reads whose span overlaps and that carry no intron.

usage: isoform_copy_join.py <out.gtf> <out.assignments.tsv> <bam> <copies.tsv> [--out prefix]
"""
import sys, subprocess, re, collections, csv

gtf_p, asg_p, bam, copies_p = sys.argv[1:5]
out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else None

tx = {}
for l in open(gtf_p):
    f = l.rstrip('\n').split('\t')
    if len(f) < 9: continue
    tid = re.search(r'transcript_id "([^"]*)"', f[8])
    if not tid: continue
    tid = tid.group(1)
    if f[2] == 'transcript':
        tx.setdefault(tid, {})['span'] = (f[0], int(f[3]) - 1, int(f[4]))
        tx[tid]['strand'] = f[6]
        r = re.search(r'reads "(\d+)"', f[8]); tx[tid]['n'] = int(r.group(1)) if r else 0
    elif f[2] == 'exon':
        tx.setdefault(tid, {}).setdefault('ex', []).append((int(f[3]) - 1, int(f[4])))
for t in tx.values():
    t.setdefault('ex', []); t['ex'].sort()
    t['chain'] = tuple((a[1], b[0]) for a, b in zip(t['ex'], t['ex'][1:]))
print(f"isoforms in the GTF: {len(tx)}")

copies = list(csv.DictReader(open(copies_p), delimiter='\t'))
asg = {r['read_name']: r for r in csv.DictReader(open(asg_p), delimiter='\t')}
print(f"assignment rows: {len(asg)}")

by_chain = collections.defaultdict(list)
for tid, t in tx.items():
    by_chain[(t['span'][0], t['chain'])].append(tid)

def introns(pos, cig):
    o = []; p = pos
    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cig):
        n = int(n)
        if op in 'M=XD': p += n
        elif op == 'N': o.append((p, p + n)); p += n
    return o

sup = collections.defaultdict(list)
contigs = sorted({t['span'][0] for t in tx.values()})
for c in contigs:
    lo = min(t['span'][1] for t in tx.values() if t['span'][0] == c)
    hi = max(t['span'][2] for t in tx.values() if t['span'][0] == c)
    res = subprocess.run(f"samtools view -F 2308 {bam} {c}:{max(lo,1)}-{hi}", shell=True, capture_output=True, text=True).stdout
    for line in res.split('\n'):
        f = line.split('\t')
        if len(f) < 6: continue
        pos = int(f[3]) - 1
        ch = tuple(introns(pos, f[5]))
        cands = by_chain.get((c, ch), ())
        if not ch:
            # an unspliced read has an EMPTY chain, which would otherwise match every single-exon isoform on
            # the contig; require the read's span to fall inside the isoform's instead
            end = pos
            for n, op in re.findall(r'(\d+)([MIDNSHP=X])', f[5]):
                if op in 'M=XDN': end += int(n)
            cands = [t for t in cands if tx[t]['span'][1] < end and pos < tx[t]['span'][2]]
        for tid in cands:
            sup[tid].append(f[0])

rows = []
for tid, t in tx.items():
    reads = sup.get(tid, [])
    seen = [asg[r] for r in reads if r in asg]
    votes = collections.Counter(r['catalog_copy_idx'] for r in seen if r['status'] == 'assigned' and r['origin_rejected'] == '0')
    abst = sum(1 for r in seen if r['status'] != 'assigned')
    best, nbest = (votes.most_common(1)[0] if votes else ('NA', 0))
    rows.append(dict(transcript=tid, chrom=t['span'][0], start=t['span'][1], end=t['span'][2],
                     n_exon=len(t['ex']), reads_gtf=t['n'], reads_matched=len(reads), reads_in_o2=len(seen),
                     assigned_copy=best, votes_for_best=nbest, votes_total=sum(votes.values()),
                     abstaining=abst,
                     purity=(round(nbest / sum(votes.values()), 3) if votes else 'NA')))
rows.sort(key=lambda r: -r['votes_total'])
dec = [r for r in rows if r['votes_total'] > 0]
print(f"\nisoforms with at least one certificate-assigned read: {len(dec)}/{len(rows)}")
if dec:
    import statistics as st
    pur = [r['purity'] for r in dec if r['purity'] != 'NA']
    print(f"  copy purity of those isoforms: median {st.median(pur):.2f}, at 1.00: {sum(1 for p in pur if p==1.0)}/{len(pur)}")
    print(f"\n  {'transcript':34s} {'copy':>5} {'votes':>7} {'purity':>7} {'abstain':>8}")
    for r in dec[:12]:
        print(f"  {r['transcript']:34s} {r['assigned_copy']:>5} {r['votes_for_best']:>3}/{r['votes_total']:<3} {r['purity']:>7} {r['abstaining']:>8}")
if out:
    with open(out + '.isoform_copy.tsv', 'w') as o:
        w = csv.DictWriter(o, fieldnames=list(rows[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(rows)
    print(f"\nwrote {out}.isoform_copy.tsv")
