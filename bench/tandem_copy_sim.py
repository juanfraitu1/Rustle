#!/usr/bin/env python3
"""Tandem-copy simulator (docs/PREREG_tandem_copy_sim_2026-09-24.md): plant k copies of a real two-exon gene into a
chr20 background at identity p and spacing D, simulate reads from every copy, map with the shipped minimap2 settings,
classify each read's exon placements, and optionally run the pipeline (assembler -> catalog -> assignment) on it.

usage: tandem_copy_sim.py --fasta chr20.fa --gtf chr20.gtf --out PREFIX [--transcript rna-NR_161305.1]
         [--copies 2] [--identity 0.99 | --sweep 0.9,0.95,0.98,0.99,0.995,1.0] [--distance 8000] [--reads 50]
         [--intron 800] [--seed 1] [--pipeline] [--bin DIR]

Per-read classes (primary alignment):
  same_copy       both exons inside the read's source copy
  other_copy      both exons inside one other copy
  cross_forward   exon 1 in an upstream copy, exon 2 in a downstream copy (one colinear chain across copies)
  cross_backward  exon 2 upstream of exon 1 — impossible as one alignment; counted from supplementary records
  unspliced       no intron in the primary alignment
  partial         >= 50 bp soft-clipped
  unmapped
plus MAPQ, AS-tied (a secondary within 0.98 of the primary AS) and a supplementary flag.
Outputs: PREFIX.<p>.fa/.bam/.reads.tsv per condition and PREFIX.summary.tsv (one row per condition)."""
import argparse, collections, os, random, re, subprocess, sys
import pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from sim_reads import simulate_reads

ap = argparse.ArgumentParser()
ap.add_argument('--fasta', required=True); ap.add_argument('--gtf', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--transcript', default='rna-NR_161305.1'); ap.add_argument('--chrom', default='chr20')
ap.add_argument('--background', default='chr20:20000000-20200000')
ap.add_argument('--copies', type=int, default=2); ap.add_argument('--identity', type=float, default=0.99); ap.add_argument('--sweep', default='')
ap.add_argument('--layout', choices=['tandem', 'interleaved'], default='tandem', help='tandem: whole-gene copies spaced --distance apart; interleaved: exon-level tandem duplication (E1 E1p ... E2 E2p, spacer --distance), where the cross-copy chain E1p->E2 has the SHORTER intron')
ap.add_argument('--distance', type=int, default=8000); ap.add_argument('--reads', type=int, default=50); ap.add_argument('--intron', type=int, default=800)
ap.add_argument('--seed', type=int, default=1); ap.add_argument('--pipeline', action='store_true'); ap.add_argument('--bin', default=os.environ.get('RUSTLE_BIN', 'target/release'))
ap.add_argument('--threads', type=int, default=3)
a = ap.parse_args()
MM2 = "minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes"
POLISH = "--assembly-polish full --polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7 --polish-retained-ratio 10"

# ---- the gene: exon 1, an 800 bp intron with its own donor/acceptor ends, exon 2
fa = pysam.FastaFile(a.fasta)
exons = []
for ln in open(a.gtf):
    f = ln.split('\t')
    if len(f) > 8 and f[0] == a.chrom and f[2] == 'exon' and f'transcript_id "{a.transcript}"' in f[8]:
        exons.append((int(f[3]), int(f[4]))); strand = f[6]
exons = sorted(set(exons)); assert len(exons) == 2, f'{a.transcript}: need exactly 2 exons, got {exons}'
(e1s, e1e), (e2s, e2e) = exons
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')
def rc(s): return s.translate(COMP)[::-1]
ex1 = fa.fetch(a.chrom, e1s - 1, e1e).upper(); ex2 = fa.fetch(a.chrom, e2s - 1, e2e).upper(); intr = fa.fetch(a.chrom, e1e, e2s - 1).upper()
half = a.intron // 2; intr = intr[:half] + intr[-half:]
if strand == '-':  # keep everything in the transcript's own orientation and plant it on the + strand of the contig
    ex1, ex2, intr = rc(ex2), rc(ex1), rc(intr)
unit = ex1 + intr + ex2
bchrom, brange = a.background.split(':'); bs, be = [int(x) for x in brange.split('-')]
background = fa.fetch(bchrom, bs - 1, be).upper()

def mutate(seq, rate, rng, protect=()):
    """substitutions at `rate`; positions in `protect` (donor/acceptor dinucleotides) are never touched"""
    s = list(seq); cand = [i for i in range(len(s)) if i not in protect]; n = int(round(rate * len(s)))
    for i in rng.sample(cand, min(n, len(cand))): s[i] = rng.choice([b for b in 'ACGT' if b != s[i]])
    return ''.join(s)

L1, LI, L2 = len(ex1), len(intr), len(ex2)
SPLICE = {L1, L1 + 1, L1 + LI - 2, L1 + LI - 1}  # GT..AG of the planted intron

def build(p, seed):
    """contig with k copies; copies[i] = dict(e1, e2 (0-based half-open exon intervals), tx (spliced transcript))."""
    rng = random.Random(seed); offset = 50000; copies = []; seq = background
    units = [unit] + [mutate(unit, 1 - p, rng, SPLICE) for _ in range(a.copies - 1)]
    if a.layout == 'tandem':
        pos = offset
        for i, u in enumerate(units):
            seq = seq[:pos] + u + seq[pos + len(u):]
            copies.append({'e1': (pos, pos + L1), 'e2': (pos + L1 + LI, pos + len(u)), 'tx': u[:L1] + u[L1 + LI:]})
            pos += len(u) + a.distance
    else:
        # E1_0 +donor [d] E1_1 +donor [d] ... | intron body | [d] acceptor+ E2_0 [d] acceptor+ E2_1 ...
        # every planted exon keeps 50 bp of its own intron flank (donor after E1, acceptor before E2), so ANY chain
        # E1_i -> E2_j is canonical; spacers are background; the intron body is copy 0's
        d = a.distance; F = 50
        head, body_mid, tail = intr[:F], intr[F:-F], intr[-F:]
        e1_block = ''; e1_pos = []
        for i, u in enumerate(units):
            e1_pos.append(offset + len(e1_block)); e1_block += u[:L1] + head + background[offset + 1000 * (i + 1): offset + 1000 * (i + 1) + d]
        e2_block = ''; e2_pos = []
        for i, u in enumerate(units):
            e2_block += background[offset + 30000 + 1000 * (i + 1): offset + 30000 + 1000 * (i + 1) + d] + tail
            e2_pos.append(offset + len(e1_block) + len(body_mid) + len(e2_block)); e2_block += u[L1 + LI:]
        body = e1_block + body_mid + e2_block
        seq = seq[:offset] + body + seq[offset + len(body):]
        for i, u in enumerate(units):
            copies.append({'e1': (e1_pos[i], e1_pos[i] + L1), 'e2': (e2_pos[i], e2_pos[i] + L2), 'tx': u[:L1] + u[L1 + LI:]})
    for c in copies: c['start'], c['end'] = c['e1'][0], c['e2'][1]
    return seq, copies

def which_copy(copies, s, e):
    """copy whose span covers >= 50% of [s,e) (used for catalog copies and assembled transcripts)"""
    best = None
    for i, c in enumerate(copies):
        o = min(e, c['end']) - max(s, c['start'])
        if o >= 0.5 * (e - s) and (best is None or o > best[1]): best = (i, o)
    return best[0] if best else None

def which_exon(copies, s, e):
    """(copy, exon) of the planted exon that overlaps [s,e) most (>= 50% of the shorter), else None"""
    best = None
    for i, c in enumerate(copies):
        for x, (a_, b_) in ((1, c['e1']), (2, c['e2'])):
            o = min(e, b_) - max(s, a_)
            if o >= 0.5 * min(e - s, b_ - a_) and (best is None or o > best[2]): best = (i, x, o)
    return best[:2] if best else None

def classify(bam, copies):
    recs = collections.defaultdict(list)
    for r in pysam.AlignmentFile(bam): recs[r.query_name].append(r)
    rows = []
    for name, rs in recs.items():
        src = int(name.split('|')[0][4:])
        prim = [r for r in rs if not r.is_secondary and not r.is_supplementary and not r.is_unmapped]
        sup = [r for r in rs if r.is_supplementary]; sec = [r for r in rs if r.is_secondary]
        if not prim: rows.append((name, src, 'unmapped', -1, 0, len(sup) > 0)); continue
        p = prim[0]; groups = [[]]; rp = p.reference_start
        for op, l in p.cigartuples:
            if op in (0, 7, 8, 2): groups[-1].append((rp, rp + l)); rp += l
            elif op == 3: rp += l; groups.append([])
        groups = [g for g in groups if g]
        clip = sum(l for op, l in p.cigartuples if op == 4)
        tied = p.mapping_quality == 0 or any(r.get_tag('AS') >= 0.98 * p.get_tag('AS') for r in sec if r.has_tag('AS'))
        if len(groups) == 1: cls = 'unspliced'
        else:
            g1, g2 = groups[0], groups[-1]
            h1 = which_exon(copies, g1[0][0], g1[-1][1]); h2 = which_exon(copies, g2[0][0], g2[-1][1])
            if h1 is None or h2 is None: cls = 'outside'
            elif h1[0] == h2[0] == src: cls = 'same_copy'
            elif h1[0] == h2[0]: cls = 'other_copy'
            else: cls = f'cross:E1@copy{h1[0]},E2@copy{h2[0]}'
        if clip >= 50 and cls in ('same_copy', 'other_copy', 'unspliced'): cls = 'partial'
        # a supplementary whose position is upstream of the primary while the read continues = backward join
        for s_ in sup:
            if s_.reference_start < p.reference_start and cls != 'cross_forward': cls = 'cross_backward'
        rows.append((name, src, cls, p.mapping_quality, int(tied), len(sup) > 0))
    return rows

def run(cmd, log):
    with open(log, 'w') as lg:
        return subprocess.run(cmd, shell=True, stdout=lg, stderr=subprocess.STDOUT).returncode

def pipeline(prefix, bam, fasta, contig_len, copies, truth_src, placed):
    out = {}
    B = a.bin
    rc_ = run(f"{B}/copy_assign --assemble-only --assembly-junctions strict {POLISH} --bam {bam} --fasta {fasta} --region sim:1-{contig_len} --out {prefix}.asm", f"{prefix}.asm.log")
    chim = ntx = 0; per_copy = collections.Counter()
    if rc_ == 0 and os.path.exists(f"{prefix}.asm.gtf"):
        tx = collections.defaultdict(list)
        for ln in open(f"{prefix}.asm.gtf"):
            f = ln.split('\t')
            if len(f) > 8 and f[2] == 'exon': tx[re.search(r'transcript_id "([^"]+)"', f[8]).group(1)].append((int(f[3]) - 1, int(f[4])))
        for t, ex in tx.items():
            ntx += 1; hits = [which_exon(copies, s_, e_) for s_, e_ in ex]; cs = {h[0] for h in hits if h}
            if len(cs) >= 2: chim += 1
            elif len(cs) == 1: per_copy[next(iter(cs))] += 1
    out.update(asm_rc=rc_, asm_tx=ntx, asm_chimeric=chim, asm_per_copy=','.join(f'{i}:{per_copy[i]}' for i in range(len(copies))))
    rc_ = run(f"{B}/gw_family_catalog --bam {bam} --fasta {fasta} --threads {a.threads} --out {prefix}.cat", f"{prefix}.cat.log")
    ncop = nfam = 0
    if rc_ == 0 and os.path.exists(f"{prefix}.cat.copies.tsv"):
        ncop = sum(1 for _ in open(f"{prefix}.cat.copies.tsv")) - 1
        nfam = sum(1 for l in open(f"{prefix}.cat.families.tsv") if not l.startswith('family_id') and int(l.split('\t')[1]) >= 2)
    out.update(cat_rc=rc_, cat_copies=ncop, cat_multicopy_families=nfam)
    acc = {'correct': 0, 'wrong': 0, 'abstain': 0}
    if ncop >= 2:
        with open(f"{prefix}.regions.txt", 'w') as f: f.write(f"sim:1-{contig_len}\n")
        rc_ = run(f"{B}/copy_assign --bam {bam} --fasta {fasta} --regions {prefix}.regions.txt --families {prefix}.cat.copies.tsv --copies-fa {prefix}.cat.copies.fa --out {prefix}.o2", f"{prefix}.o2.log")
        if rc_ == 0 and os.path.exists(f"{prefix}.o2.assignments.tsv"):
            cat = {}
            for l in open(f"{prefix}.cat.copies.tsv"):
                f = l.rstrip('\n').split('\t')
                if f[0] != 'family_id': cat[(f[0], f[1])] = (int(f[4]), int(f[5]))
            by = collections.defaultdict(list)
            for l in open(f"{prefix}.o2.assignments.tsv"):
                f = l.rstrip('\n').split('\t')
                if f[0] == 'read_name': hdr = f; continue
                by[f[0]].append(dict(zip(hdr, f)))
            for name, src in truth_src.items():
                if name not in by:  # not a tied read: the aligner placed it; score that placement
                    acc['unique_' + ('correct' if placed.get(name) == src else 'wrong')] = acc.get('unique_' + ('correct' if placed.get(name) == src else 'wrong'), 0) + 1; continue
                rows = [r for r in by.get(name, []) if r['status'] == 'assigned' and r['origin_rejected'] == '0']
                if not rows: acc['abstain'] += 1; continue
                ok = any(which_copy(copies, *cat.get((r['family_id'], r['catalog_copy_idx']), (0, 0))) == src for r in rows)
                acc['correct' if ok else 'wrong'] += 1
        out['o2_rc'] = rc_
    out.update(o2_tied_correct=acc['correct'], o2_tied_wrong=acc['wrong'], o2_tied_abstain=acc['abstain'], o2_unique_correct=acc.get('unique_correct', 0), o2_unique_wrong=acc.get('unique_wrong', 0))
    return out

sweep = [float(x) for x in a.sweep.split(',')] if a.sweep else [a.identity]
summary = []
for p in sweep:
    prefix = f"{a.out}.{a.layout}.k{a.copies}.p{p}"
    seq, copies = build(p, a.seed)
    with open(prefix + '.fa', 'w') as f: f.write(f">sim\n{seq}\n")
    subprocess.run(f"samtools faidx {prefix}.fa", shell=True, check=True)
    truth_src = {}
    with open(prefix + '.fq', 'w') as fq:
        for i, c in enumerate(copies):
            rng = random.Random(a.seed * 101 + i)
            for n, (rd, q) in enumerate(simulate_reads(c['tx'], a.reads, err=0.001, indel=0.0003, seed=a.seed * 31 + i, trunc_frac=0.10)):
                j5, j3 = rng.randint(0, 30), rng.randint(0, 30); rd = rd[j5:len(rd) - j3]; q = q[j5:len(q) - j3]
                name = f"copy{i}|{n}"; fq.write(f"@{name}\n{rd}\n+\n{q}\n"); truth_src[name] = i
    subprocess.run(f"{MM2} -t {a.threads} {prefix}.fa {prefix}.fq 2>/dev/null | samtools sort -o {prefix}.bam - 2>/dev/null && samtools index {prefix}.bam", shell=True, check=True)
    rows = classify(prefix + '.bam', copies)
    with open(prefix + '.reads.tsv', 'w') as f:
        f.write('read\tsource_copy\tclass\tmapq\ttied\tsupplementary\n')
        for r in rows: f.write('\t'.join(str(x) for x in r) + '\n')
    n = len(rows); cls = collections.Counter(r[2] for r in rows)
    by_src = {i: collections.Counter(r[2] for r in rows if r[1] == i) for i in range(len(copies))}
    row = {'copies': a.copies, 'identity': p, 'distance': a.distance, 'reads': n, 'mapq0': sum(1 for r in rows if r[3] == 0), 'tied': sum(r[4] for r in rows), 'supplementary': sum(1 for r in rows if r[5])}
    for k in ('same_copy', 'other_copy', 'cross_backward', 'unspliced', 'partial', 'outside', 'unmapped'): row[k] = cls[k]
    crosses = {k: v for k, v in cls.items() if k.startswith('cross:')}
    row['cross_copy'] = sum(crosses.values()); row['cross_kinds'] = ';'.join(f'{k[6:]}={v}' for k, v in sorted(crosses.items()))
    row['cross_from_copy0'] = sum(v for k, v in by_src[0].items() if k.startswith('cross:')); row['cross_from_others'] = sum(v for i in range(1, len(copies)) for k, v in by_src[i].items() if k.startswith('cross:'))
    if a.pipeline:
        placed = {r[0]: (r[1] if r[2] == 'same_copy' else (None if not r[2].startswith('other') else -1)) for r in rows}
        # aligner placement per read: the copy holding the primary's first exon (None when cross/unspliced/unmapped)
        recs = {}
        for r_ in pysam.AlignmentFile(prefix + '.bam'):
            if not (r_.is_secondary or r_.is_supplementary or r_.is_unmapped): recs[r_.query_name] = which_copy(copies, r_.reference_start, r_.reference_end)
        row.update(pipeline(prefix, prefix + '.bam', prefix + '.fa', len(seq), copies, truth_src, recs))
    summary.append(row)
    print(f"[{a.layout} k={a.copies} p={p}] reads {n}: " + ', '.join(f"{k} {v}" for k, v in row.items() if k not in ('copies', 'identity', 'distance', 'reads') and v), flush=True)
keys = list(summary[0].keys())
with open(f"{a.out}.{a.layout}.k{a.copies}.summary.tsv", 'w') as f:
    f.write('\t'.join(keys) + '\n')
    for r in summary: f.write('\t'.join(str(r.get(k, '')) for k in keys) + '\n')
print(f"wrote {a.out}.{a.layout}.k{a.copies}.summary.tsv")
