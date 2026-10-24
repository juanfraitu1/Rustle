#!/usr/bin/env python3
"""Score the in-between operating point against StringTie and FLAIR, genome-wide and on the multi-copy strata
(docs/archive/2026-10/merged/PREREG_inbetween_operating_point_2026-10-11.md).

usage: score_op.py --species S --ref REF.gtf[.gz] --exclude CHROM --strata PREFIX --work DIR --arm NAME=GTF [--arm ...]
Arms named OP_primary, OP_good, StringTie and FLAIR decide claims B1, B2, C1-C3; any other arm is reported only.
Every GTF and the reference lose the excluded (development) chromosome before gffcompare. Stratum sensitivity = multi-exon
reference transcripts of the stratum's genes an arm matches exactly (gffcompare refmap `=`); stratum precision = arm
transcripts with an exon overlapping a stratum gene's exon that match a reference transcript exactly (tmap `=`) / all such."""
import argparse, collections, gzip, os, re, subprocess

ap = argparse.ArgumentParser()
for a in ('--species', '--ref', '--exclude', '--strata', '--work'):
    ap.add_argument(a, required=True)
ap.add_argument('--arm', action='append', required=True)
ap.add_argument('--gffcompare', default='gffcompare')
args = ap.parse_args()
os.makedirs(args.work, exist_ok=True)
BIN = 10_000

def opener(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def drop_chrom(src, dst):
    if os.path.exists(dst):
        return
    with opener(src) as i, open(dst + '.tmp', 'w') as o:
        for ln in i:
            if not ln.startswith('#') and ln.split('\t', 1)[0] != args.exclude:
                o.write(ln)
    os.replace(dst + '.tmp', dst)

ref = os.path.join(args.work, 'ref.gtf'); drop_chrom(args.ref, ref)
# strata
strat_tx = collections.defaultdict(set); strat_genes = collections.defaultdict(set)
for ln in open(args.strata + '.transcripts.tsv'):
    f = ln.rstrip('\n').split('\t')
    if f[0] != 'transcript_id':
        strat_tx[f[3]].add(f[0]); strat_genes[f[3]].add((f[1], f[2]))
gene_of = {}
index = collections.defaultdict(list)
for ln in open(ref):
    f = ln.split('\t')
    if len(f) < 9 or f[2] != 'exon':
        continue
    gk = (f[0], re.search(r'gene_id "([^"]*)"', f[8]).group(1) or re.search(r'transcript_id "([^"]+)"', f[8]).group(1))
    for s, gs in strat_genes.items():
        if gk in gs:
            a, b = int(f[3]) - 1, int(f[4])
            for bn in range(a // BIN, (b - 1) // BIN + 1):
                index[(f[0], bn)].append((a, b, s))

def stats(path):
    out = {}
    for ln in open(path):
        m = re.match(r'\s*(Base|Exon|Intron|Intron chain|Transcript|Locus) level:\s+([\d.]+)\s+\|\s+([\d.]+)', ln)
        if m:
            out[m.group(1)] = (float(m.group(2)), float(m.group(3)))
    return out

res = {}
for spec in args.arm:
    name, gtf = spec.split('=', 1)
    q = os.path.join(args.work, name + '.gtf'); drop_chrom(gtf, q)
    pre = os.path.join(args.work, 'g_' + name)
    if not os.path.exists(pre + '.stats'):
        subprocess.run([args.gffcompare, '-r', ref, '-o', pre, q], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    st = stats(pre + '.stats')
    eq_ref = set(); eq_q = set()
    for ln in open(os.path.join(args.work, f'g_{name}.{name}.gtf.refmap')):
        f = ln.rstrip('\n').split('\t')
        if len(f) >= 4 and f[2] == '=':
            eq_ref.add(f[1])
    for ln in open(os.path.join(args.work, f'g_{name}.{name}.gtf.tmap')):
        f = ln.split('\t')
        if len(f) > 4 and f[2] == '=':
            eq_q.add(f[4])
    # arm transcripts overlapping each stratum
    over = collections.defaultdict(set)
    for ln in open(q):
        f = ln.split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); a, b = int(f[3]) - 1, int(f[4])
        for bn in range(a // BIN, (b - 1) // BIN + 1):
            for xs, xe, s in index.get((f[0], bn), ()):
                if xs < b and xe > a:
                    over[s].add(t)
    strata = {}
    for s in ('S1', 'S2'):
        ov = over[s]
        strata[s] = (len(strat_tx[s] & eq_ref), len(strat_tx[s]), sum(t in eq_q for t in ov), len(ov))
    res[name] = (st, strata)

print(f"[{args.species}] genome-wide minus {args.exclude} (sensitivity / precision)")
print(f"{'arm':14s} {'base':>13s} {'intron':>13s} {'chain':>13s} {'transcript':>13s} {'locus':>13s}")
for name, (st, _) in res.items():
    print(f"{name:14s} " + ' '.join(f"{st[k][0]:5.1f} / {st[k][1]:5.1f}" for k in ('Base', 'Intron', 'Intron chain', 'Transcript', 'Locus')))
for s in ('S1', 'S2'):
    tot = len(strat_tx[s])
    print(f"stratum {s}: {tot:,} multi-exon reference transcripts in {len(strat_genes[s]):,} genes")
    for name, (_, strata) in res.items():
        m, t, c, o = strata[s]
        print(f"   {name:14s} matched {m:5,} / {t:,}   precision {c:,} / {o:,} = {c / max(1, o):.3f}")
need = {'OP_primary', 'OP_good', 'StringTie', 'FLAIR'}
if need <= set(res):
    op, opg, stt, fl = res['OP_primary'][0], res['OP_good'][1], res['StringTie'], res['FLAIR']
    chain = lambda st: st['Intron chain']
    flair_pr = chain(fl[0])[1]
    claims = [('B1 OP_primary chain sensitivity > StringTie and > FLAIR', chain(op)[0] > chain(stt[0])[0] and chain(op)[0] > chain(fl[0])[0]),
              ('B2 OP_primary chain precision > FLAIR', chain(op)[1] > flair_pr)]
    for s, lab in (('S1', 'C1'), ('S2', 'C2')):
        tot = len(strat_tx[s])
        if tot < 20:
            claims.append((f'{lab} {s}: UNDERPOWERED ({tot} reference transcripts < 20)', None)); continue
        claims.append((f'{lab} {s}: OP_good matches more than StringTie and FLAIR', opg[s][0] > stt[1][s][0] and opg[s][0] > fl[1][s][0]))
        claims.append((f'C3 {s}: OP_good stratum precision >= FLAIR genome-wide chain precision ({flair_pr:.1f}%)', 100 * opg[s][2] / max(1, opg[s][3]) >= flair_pr))
    for k, v in claims:
        print(f"  {'n/a ' if v is None else ('PASS' if v else 'FAIL')}  {k}")
