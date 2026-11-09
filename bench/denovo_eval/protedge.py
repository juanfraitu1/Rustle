#!/usr/bin/env python3
"""Row 1590 (docs/PREREG_prot_sd_edges_2026-11-08.md): protein edges between SD-linked de novo loci, fed to
`mcl_families --extra-edges` (merged before --sd-damp).

  protedge.py prep  GTF LOCI_GFF3 GENOME OUT      -> OUT.faa (<= 3 longest ORFs >= 80 aa per locus) + OUT.loci.pkl
  (diamond blastp -q OUT.faa -d OUT --more-sensitive -e 1e-10 -k 0 --outfmt 6 qseqid sseqid pident length qlen slen evalue)
  protedge.py edges OUT SEDEF_BED GATED_PAF       -> OUT.pairs.tsv (locusA locusB best_identity), protein-linked AND
                                                     SD-linked; prints the port check (SD-linked gated-PAF edges)
Locus names are `chrom:start-end` from the gff3 gene lines (= the PAF / loci.fa names). Write the edge file with
  awk -v w=W '{print $1"\\t"$2"\\t"w}' OUT.pairs.tsv
--selftest checks the transcript -> locus mapping, the SD rule and the protein filter on toy data.
"""
import sys, os, re, collections, bisect, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from protedge_step0 import longest_orf, COMP


def read_loci(gff):
    loci, ex = {}, collections.defaultdict(list)
    for line in open(gff):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        if f[2] == 'gene':
            gid = re.search(r'ID=([^;]+)', f[8]).group(1)
            loci[gid] = (f'{f[0]}:{f[3]}-{f[4]}', f[0], f[6])
        elif f[2] == 'exon':
            ex[re.search(r'Parent=([^;]+)', f[8]).group(1)].append((int(f[3]) - 1, int(f[4])))
    return {name: (chrom, strand, sorted(ex[gid])) for gid, (name, chrom, strand) in loci.items() if ex[gid]}


def locus_of(index, chrom, strand, exons):
    """the locus (same strand) sharing the most exon bases with the transcript, or None."""
    v = index.get((chrom, strand))
    if not v:
        return None
    starts, items, span = v
    lo, hi = exons[0][0], exons[-1][1]
    best = None
    for i in range(bisect.bisect_left(starts, lo - span), len(items)):
        s, e, name, bl = items[i]
        if s >= hi:
            break
        if e <= lo:
            continue
        ov = sum(max(0, min(b, d) - max(a, c)) for a, b in exons for c, d in bl)
        if ov and (best is None or (ov, name) > best):
            best = (ov, name)
    return best[1] if best else None


def build_index(loci):
    idx = collections.defaultdict(list)
    for name, (chrom, strand, bl) in loci.items():
        idx[(chrom, strand)].append((bl[0][0], bl[-1][1], name, bl))
    out = {}
    for k, v in idx.items():
        v.sort()
        out[k] = ([x[0] for x in v], v, max(x[1] - x[0] for x in v))
    return out


def prep(gtf, gff, genome, out):
    import pysam
    loci = read_loci(gff); index = build_index(loci)
    tx = collections.defaultdict(list); meta = {}
    for line in open(gtf):
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        tx[t].append((int(f[3]) - 1, int(f[4]))); meta[t] = (f[0], f[6] if f[6] in '+-' else '+')
    fa = pysam.FastaFile(genome); orfs = collections.defaultdict(set); unplaced = 0
    for t, ex in tx.items():
        ex.sort(); chrom, strand = meta[t]; name = locus_of(index, chrom, strand, ex)
        if name is None:
            unplaced += 1; continue
        s = ''.join(fa.fetch(chrom, a, b) for a, b in ex).upper()
        p = longest_orf(s.translate(COMP)[::-1] if strand == '-' else s)
        if p:
            orfs[name].add(p)
    n = 0
    with open(out + '.faa', 'w') as w:
        for name, ps in orfs.items():
            for i, p in enumerate(sorted(ps, key=lambda x: (-len(x), x))[:3]):
                w.write(f'>{name}|{i}\n{p}\n'); n += 1
    pickle.dump(loci, open(out + '.loci.pkl', 'wb'))
    print(f'loci {len(loci)}, transcripts {len(tx)} (unplaced {unplaced}), loci with an ORF {len(orfs)}, proteins {n}')


def sd_holders(loci, bed):
    """locus -> {(sd pair, side)} whose copy holds >= half its exon bases AND >= half its span (SdPairs::holds)."""
    by = collections.defaultdict(list)
    for i, l in enumerate(open(bed)):
        if l.startswith('#') or not l.strip():
            continue
        f = l.split('\t'); by[f[0]].append((int(f[1]), int(f[2]), i, 0)); by[f[3]].append((int(f[4]), int(f[5]), i, 1))
    for v in by.values():
        v.sort()
    st = {c: [x[0] for x in v] for c, v in by.items()}; mx = {c: max(x[1] - x[0] for x in v) for c, v in by.items()}
    H = {}
    for name, (c, _, bl) in loci.items():
        lo, hi = bl[0][0], bl[-1][1]; tot = sum(e - s for s, e in bl); v = by.get(c, []); h = set()
        for i in range(bisect.bisect_left(st.get(c, []), lo - mx.get(c, 0)), len(v)):
            s, e, p, side = v[i]
            if s >= hi:
                break
            if e <= lo:
                continue
            ins = sum(max(0, min(y, e) - max(x, s)) for x, y in bl)
            if tot and 2 * ins >= tot and 2 * max(0, min(hi, e) - max(lo, s)) >= hi - lo:
                h.add((p, side))
        H[name] = h
    return H


def sd_linked(H, a, b):
    return any((p, 1 - side) in H.get(b, ()) for p, side in H.get(a, ()))


def protein_pairs(lines, min_id=50.0):
    """best identity per locus pair passing e <= 1e-10, both covers >= 50%, identity >= min_id."""
    best = {}
    for line in lines:
        q, s, pid, ln, ql, sl, ev = line.split('\t')[:7]
        a, b = q.rsplit('|', 1)[0], s.rsplit('|', 1)[0]
        if a == b or float(ev) > 1e-10 or float(pid) < min_id or int(ln) < 0.5 * int(ql) or int(ln) < 0.5 * int(sl):
            continue
        k = (min(a, b), max(a, b)); best[k] = max(best.get(k, 0.0), float(pid))
    return best


def edges(out, bed, paf):
    loci = pickle.load(open(out + '.loci.pkl', 'rb')); H = sd_holders(loci, bed)
    pe = set()
    for l in open(paf):
        f = l.split('\t')
        if len(f) > 5 and f[0] != f[5]:
            pe.add((min(f[0], f[5]), max(f[0], f[5])))
    port = sum(1 for a, b in pe if sd_linked(H, a, b))
    print(f'port check: {port} of {len(pe)} gated-PAF node pairs SD-linked (compare the binary sd_damp_edges_linked)')
    prot = protein_pairs(open(out + '.allvall.tsv'))
    kept = {k: v for k, v in prot.items() if sd_linked(H, *k)}
    with open(out + '.pairs.tsv', 'w') as w:
        for (a, b), v in sorted(kept.items()):
            w.write(f'{a}\t{b}\t{v:.1f}\n')
    print(f'protein-linked locus pairs {len(prot)}, SD-linked among them {len(kept)} '
          f'(already a PAF edge {sum(1 for k in kept if k in pe)}) -> {out}.pairs.tsv')


def selftest():
    loci = {'c:1-100': ('c', '+', [(0, 40), (60, 100)]), 'c:201-300': ('c', '+', [(200, 300)]),
            'c:501-600': ('c', '-', [(500, 600)])}
    idx = build_index(loci)
    assert locus_of(idx, 'c', '+', [(10, 30), (70, 90)]) == 'c:1-100'
    assert locus_of(idx, 'c', '-', [(10, 30)]) is None                       # strand must match
    assert locus_of(idx, 'c', '+', [(80, 210)]) == 'c:1-100'                 # 20 bp vs 10 bp: the larger overlap wins
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as t:
        t.write('#hdr\nc\t0\t110\tc\t190\t310\t...\n')
    H = sd_holders(loci, t.name); os.unlink(t.name)
    assert sd_linked(H, 'c:1-100', 'c:201-300') and sd_linked(H, 'c:201-300', 'c:1-100')
    assert not sd_linked(H, 'c:1-100', 'c:501-600')
    rows = ['c:1-100|0\tc:201-300|1\t80.0\t100\t120\t150\t1e-30', 'c:1-100|0\tc:501-600|0\t40.0\t100\t120\t150\t1e-30',
            'c:1-100|0\tc:201-300|0\t90.0\t50\t120\t150\t1e-30']               # cover 50/120 < 0.5 -> dropped
    assert protein_pairs(rows) == {('c:1-100', 'c:201-300'): 80.0}
    print('selftest ok')


if __name__ == '__main__':
    cmd = sys.argv[1]
    selftest() if cmd == '--selftest' else prep(*sys.argv[2:6]) if cmd == 'prep' else edges(*sys.argv[2:5])
