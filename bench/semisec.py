"""SEMISEC (docs/LEDGER_O1.md): semi-guided silent copies of expressed genes from their reads' secondary alignments.

  semisec.py dump BAM OUT_DIR prim|sec CONTIG... -> OUT_DIR/{prim,sec}/CONTIG1.tsv
  semisec.py run BAM GENOME_FA GFF OUT_DIR     (reads the dump; fetches seed read sequences from BAM)

Steps 1-4 of the prereg read only the BAM and the genome; the GFF is read afterwards, for the truth columns and the report.
"""
import bisect
import collections
import glob
import json
import os
import re
import subprocess
import sys

import pysam

CIG = re.compile(r'(\d+)([MIDNSHP=X])')
REF = re.compile(r'(\d+)[MDN=X]')
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
STOP = {'TAA', 'TAG', 'TGA'}
MIN_READS, SEC_RATIO, MIN_CODONS = 2, 0.8, 100


def rc(s):
    return s.translate(COMP)[::-1]


def ref_len(cig):
    return sum(map(int, REF.findall(cig)))


AWK = r'''{a=0; for(i=12;i<=NF;i++) if(substr($i,1,5)=="AS:i:"){a=substr($i,6);break}
 print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"a}'''


def dump(bam, out, kind, *contigs):
    """Stage 0: QNAME FLAG RNAME POS MAPQ CIGAR AS of primaries (kind prim) or secondaries (kind sec) of CONTIGS, one file per
    call in OUT_DIR/KIND/ (batch the contigs so each call stays inside the machine's time limit)."""
    flags = '-F 2308' if kind == 'prim' else '-f 256'
    os.makedirs(f'{out}/{kind}', exist_ok=True)
    subprocess.run(f"samtools view -@4 {flags} {bam} {' '.join(contigs)} | mawk -F'\\t' '{AWK}' > {out}/{kind}/{contigs[0]}.tsv",
                   shell=True, check=True)


def sam(path):
    for fn in [path] if os.path.isfile(path) else sorted(glob.glob(f'{path}/*.tsv')):
        for line in open(fn):
            f = line.rstrip('\n').split('\t')
            yield f, int(f[6])


def merge(ivs):
    """ivs: list of (start, end, payload); returns [(start, end, [payloads])] merged by overlap."""
    out = []
    for s, e, x in sorted(ivs):
        if out and s < out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
            out[-1][2].append(x)
        else:
            out.append([s, e, [x]])
    return out


class Index:
    """Count / find intervals overlapping a query (sorted starts + max length)."""

    def __init__(self, ivs):
        self.v = sorted(ivs)
        self.s = [x[0] for x in self.v]
        self.m = max((x[1] - x[0] for x in self.v), default=0)

    def over(self, s, e):
        i = bisect.bisect_left(self.s, s - self.m)
        out = []
        while i < len(self.v) and self.v[i][0] < e:
            if self.v[i][1] > s:
                out.append(self.v[i])
            i += 1
        return out


def longest_orf(seq):
    best = (0, 0)
    for fr in range(3):
        start = None
        for i in range(fr, len(seq) - 2, 3):
            c = seq[i:i + 3]
            if start is None and c == 'ATG':
                start = i
            elif start is not None and c in STOP:
                if i + 3 - start > best[1] - best[0]:
                    best = (start, i + 3)
                start = None
    return best


def candidate_in_read_sense(cig, rpos, reverse, qlen, genome, chrom):
    """Events in original-read order: (orig_pos or None, candidate bases). Matches carry 1 base at a read position;
    D runs carry candidate bases with no read position (attached to the preceding read position); I runs give read
    positions with no candidate base."""
    ev, q, r = [], 0, rpos
    for n, op in CIG.findall(cig):
        n = int(n)
        if op in 'M=X':
            ref = genome.fetch(chrom, r, r + n).upper()
            for j in range(n):
                ev.append((q + j, ref[j]))
            q += n
            r += n
        elif op == 'D':
            ev.append((q - 0.5, genome.fetch(chrom, r, r + n).upper()))
            r += n
        elif op == 'N':
            r += n
        elif op in 'I':
            for j in range(n):
                ev.append((q + j, ''))
            q += n
        elif op in 'SH':
            q += n
    if reverse:
        ev = [((qlen - 1 - p), b.translate(COMP)[::-1]) for p, b in reversed(ev)]
    return ev, q


def orf_status(ev, o0, o1):
    """Translate the candidate over the read's ORF [o0, o1); DISRUPTED on a frameshift or a premature in-frame stop."""
    seg = [(p, b) for p, b in ev if o0 <= p < o1]
    if not seg:
        return 'UNCOVERED', 0.0
    covered = sum(1 for p, b in seg if b and float(p).is_integer())
    first = min(int(p) for p, b in seg if float(p).is_integer())
    lead = (first - o0) % 3  # skip to the first codon boundary of the seed frame
    cand = ''.join(b for _, b in seg)
    net = sum(len(b) for p, b in seg if not float(p).is_integer()) - sum(1 for p, b in seg if b == '')
    frac = covered / (o1 - o0)
    if net % 3:
        return 'DISRUPTED', frac
    cs = cand[(3 - lead) % 3:]
    codons = [cs[i:i + 3] for i in range(0, len(cs) - 2, 3)]
    if any(c in STOP for c in codons[:-1]):
        return 'DISRUPTED', frac
    return 'INTACT', frac


def run(bam, genome_fa, gff, out):
    # 1. seeds
    prim = {}
    by = collections.defaultdict(list)
    for f, AS in sam(f'{out}/prim' if os.path.isdir(f'{out}/prim') else f'{out}/prim.tsv'):
        q, flag, chrom, pos, cig = f[0], int(f[1]), f[2], int(f[3]) - 1, f[5]
        end = pos + ref_len(cig)
        strand = '-' if flag & 16 else '+'
        prim[q] = (chrom, strand, pos, end, AS, 'N' in cig, flag & 16 != 0)
        by[(chrom, strand)].append((pos, end, q))
    seeds, read_seed = {}, {}
    for key, ivs in by.items():
        for s, e, qs in merge(ivs):
            if len(set(qs)) >= MIN_READS:
                sid = f'{key[0]}:{s}-{e}{key[1]}'
                seeds[sid] = (key[0], key[1], s, e, len(set(qs)))
                for q in qs:
                    read_seed[q] = (key[0], s, e)
    pidx_all = collections.defaultdict(list)
    for (c, _), v in by.items():
        pidx_all[c] += [(s, e) for s, e, _ in v]
    pidx = {c: Index(v) for c, v in pidx_all.items()}
    del by
    # 2. candidate hits
    hits = collections.defaultdict(list)
    n_sec = n_kept = 0
    for f, AS in sam(f'{out}/sec' if os.path.isdir(f'{out}/sec') else f'{out}/sec.tsv'):
        n_sec += 1
        q = f[0]
        if q not in read_seed:
            continue
        p = prim[q]
        if AS < SEC_RATIO * p[4]:
            continue
        chrom, pos, cig, flag = f[2], int(f[3]) - 1, f[5], int(f[1])
        end = pos + ref_len(cig)
        sc, ss, se = read_seed[q]
        if chrom == sc and pos < se and end > ss:
            continue
        n_kept += 1
        strand = '-' if flag & 16 else '+'
        hits[(chrom, strand)].append((pos, end, (q, AS, cig, pos, flag & 16 != 0)))
    cands = []
    for (chrom, strand), ivs in hits.items():
        for s, e, xs in merge(ivs):
            reads = {x[0] for x in xs}
            if len(reads) < MIN_READS:
                continue
            n_prim = len(pidx[chrom].over(s, e)) if chrom in pidx else 0
            best = max(xs, key=lambda x: x[1])
            cands.append({'chrom': chrom, 'strand': strand, 'start': s, 'end': e, 'reads': len(reads),
                          'primaries': n_prim, 'silent': n_prim < MIN_READS, 'best': best})
    del hits
    # 4. class (silent only) - read sequences of the best hits' primaries
    want = {c['best'][0] for c in cands if c['silent']}
    seqs = {}
    open(f'{out}/want.txt', 'w').write(''.join(q + '\n' for q in want))
    p = subprocess.Popen(['samtools', 'view', '-@4', '-F', '2308', '-N', f'{out}/want.txt', bam],
                         stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.split('\t', 10)
        seqs[f[0]] = rc(f[9]) if int(f[1]) & 16 else f[9]
    p.wait()
    genome = pysam.FastaFile(genome_fa)
    for c in cands:
        if not c['silent']:
            continue
        q, AS, cig, pos, rev = c['best']
        read = seqs.get(q, '')
        c['processed'] = prim[q][5] and 'N' not in cig
        o0, o1 = longest_orf(read)
        coding = (o1 - o0) // 3 >= MIN_CODONS
        # the secondary's SEQ orientation: original read if not reverse
        ev, qlen = candidate_in_read_sense(cig, pos, rev, len(read), genome, c['chrom'])
        if coding:
            c['orf'], c['orf_cov'] = orf_status(ev, o0, o1)
        else:
            c['orf'], c['orf_cov'] = 'LNC-LIKE', 0.0
        c['class'] = ('PROCESSED' if c['processed'] else 'UNPROCESSED') + '/' + c['orf']
        c['best'] = q
    for c in cands:
        if not c['silent']:
            c['best'] = c['best'][0]
    json.dump({'seeds': len(seeds), 'secondaries': n_sec, 'kept_hits': n_kept, 'cands': cands},
              open(f'{out}/cands.json', 'w'))
    report(out, cands, seeds, gff, pidx)


def biotypes(gff):
    genes = collections.defaultdict(list)
    for line in open(gff):
        if line[0] == '#':
            continue
        f = line.split('\t')
        if len(f) < 9 or f[2] not in ('gene', 'pseudogene'):
            continue
        a = dict(x.split('=', 1) for x in f[8].strip().split(';') if '=' in x)
        bt = a.get('gene_biotype', 'other')
        bt = 'pseudogene' if 'pseudogene' in bt else bt if bt in ('protein_coding', 'lncRNA') else 'other'
        genes[f[0]].append((int(f[3]) - 1, int(f[4]), bt, a.get('description', '').replace('%2C', ','), a.get('ID', '')))
    return {c: Index(v) for c, v in genes.items()}


def report(out, cands, seeds, gff, pidx):
    gi = biotypes(gff)

    def top(c, s, e):
        hs = gi[c].over(s, e) if c in gi else []
        if not hs:
            return ('none', '', '')
        h = max(hs, key=lambda h: min(e, h[1]) - max(s, h[0]))
        return h[2], h[3], h[4]

    sil = [c for c in cands if c['silent']]
    for c in sil:
        c['truth'], c['truth_desc'], c['truth_id'] = top(c['chrom'], c['start'], c['end'])
    tab = collections.Counter((c['class'], c['truth']) for c in sil)
    classes = sorted({c['class'] for c in sil})
    bts = ['pseudogene', 'protein_coding', 'lncRNA', 'other', 'none']
    lines = [f"seeds {len(seeds)}; candidates {len(cands)} (silent {len(sil)})", '',
             '| class | n | ' + ' | '.join(bts) + ' |', '|---|---|' + '---|' * len(bts)]
    for k in classes:
        n = sum(tab[(k, b)] for b in bts)
        lines.append(f'| {k} | {n} | ' + ' | '.join(f'{tab[(k, b)]} ({tab[(k, b)] / n:.2f})' for b in bts) + ' |')
    # recall: annotated silent pseudogenes named "X-like" after an expressed seed's protein-coding gene
    seed_desc = set()
    for sid, (c, st, s, e, n) in seeds.items():
        for h in (gi[c].over(s, e) if c in gi else []):
            if h[2] == 'protein_coding' and h[3]:
                seed_desc.add(h[3])
    cidx = collections.defaultdict(list)
    for c in sil:
        cidx[c['chrom']].append((c['start'], c['end']))
    cidx = {k: Index(v) for k, v in cidx.items()}
    found = tot = 0
    for chrom, idx in gi.items():
        for s, e, bt, desc, gid in idx.v:
            if bt != 'pseudogene' or not desc.endswith('-like') or desc[:-5] not in seed_desc:
                continue
            if chrom in pidx and len(pidx[chrom].over(s, e)) >= MIN_READS:
                continue  # expressed copy: not silent, not in the recall denominator
            tot += 1
            found += bool(cidx.get(chrom) and cidx[chrom].over(s, e))
    lines += ['', f'recall: annotated silent pseudogenes "X-like" of an expressed seed gene, overlapped by a silent candidate: '
                  f'{found} / {tot} ({found / max(tot, 1):.3f})']
    # stage-2 bars (LEDGER_O1 "SEMISEC stage-2 bars"): fractions among annotated candidates (truth != none), pooled over PROCESSED/UNPROCESSED
    def frac(bt, orf):
        xs = [c for c in sil if c['truth'] != 'none' and c['class'].endswith('/' + orf)]
        return sum(c['truth'] == bt for c in xs) / max(len(xs), 1), len(xs)
    pd, nd = frac('pseudogene', 'DISRUPTED')
    pi, ni = frac('pseudogene', 'INTACT')
    pdd = [c for c in sil if c['class'] == 'PROCESSED/DISRUPTED']
    ppd = sum(c['truth'] == 'pseudogene' for c in pdd) / max(len(pdd), 1)
    ll, nl = frac('lncRNA', 'LNC-LIKE')
    lc = [c for c in sil if c['truth'] != 'none' and c['class'].split('/')[1] in ('DISRUPTED', 'INTACT')]
    lcf = sum(c['truth'] == 'lncRNA' for c in lc) / max(len(lc), 1)
    rec = found / max(tot, 1)
    bars = {'B1 pseudogene DISRUPTED/INTACT >= 1.5': (pd / max(pi, 1e-9), pd / max(pi, 1e-9) >= 1.5),
            'B2 PROCESSED/DISRUPTED pseudogene >= 0.40': (ppd, ppd >= 0.40),
            'B3 recall >= 0.10': (rec, rec >= 0.10),
            'B4 (separate) lncRNA LNC-LIKE/coding >= 1.25': (ll / max(lcf, 1e-9), ll / max(lcf, 1e-9) >= 1.25)}
    lines += ['', f'annotated-only: pseudogene DISRUPTED {pd:.3f} (n {nd}) INTACT {pi:.3f} (n {ni}); '
                  f'lncRNA LNC-LIKE {ll:.3f} (n {nl}) coding {lcf:.3f} (n {len(lc)})']
    lines += [f'{k}: {v:.3f} {"PASS" if ok else "FAIL"}' for k, (v, ok) in bars.items()]
    lines.append('VERDICT (B1-B3): ' + ('PASS' if all(ok for k, (_, ok) in bars.items() if not k.startswith('B4')) else 'FAIL'))
    none = [c for c in sil if c['truth'] == 'none']
    sizes = sorted(c['end'] - c['start'] for c in none)
    lines.append(f"unannotated silent candidates: {len(none)}; median span {sizes[len(sizes) // 2] if sizes else 0} bp")
    open(f'{out}/report.md', 'w').write('\n'.join(lines) + '\n')
    json.dump(sil, open(f'{out}/silent.json', 'w'))
    print('\n'.join(lines))


if __name__ == '__main__':
    if sys.argv[1] == 'dump':
        dump(*sys.argv[2:])
    elif sys.argv[1] == 'run':
        run(*sys.argv[2:6])
