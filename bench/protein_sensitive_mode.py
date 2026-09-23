#!/usr/bin/env python3
"""SENSITIVE MODE: protein edges from translated DE NOVO locus representatives.
Per `docs/PREREG_protein_sensitive_mode_2026-09-22.md` (§6y3, md5 `e9278fdd`).

r906 measured protein edges from ANNOTATED CDS as ⚠PARTIAL (+13.8 pts pair coverage, 0 cross-family
held-out) and left them unadopted. This translates the de novo locus instead, which is the only route that
reaches pseudogenes: 0 of 457 chr16 pseudogenes have an annotated CDS.

⚠⚠ TRUTH IS SOTO, NOT THE PROTEIN REFEREE. The referee is itself built from translated CDS clustered by
protein homology, so scoring protein-derived edges against it is circular by construction.
⚠ A six-frame longest ORF is non-zero for ANY sequence, so ORF yield is not evidence. The NULL arm
(codon-shuffled ORFs) must produce ~nothing; if it fires, the signal is chance similarity.
"""
import argparse, collections, csv, random, re, subprocess, os

CODON = {}
_b = 'TCAG'
_a = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
_i = 0
for _x in _b:
    for _y in _b:
        for _z in _b:
            CODON[_x + _y + _z] = _a[_i]; _i += 1
COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def longest_orf(seq):
    """longest stop-free stretch over 6 frames, returned as protein."""
    best = ''
    for s in (seq, seq.translate(COMP)[::-1]):
        for fr in range(3):
            cur = []
            for j in range(fr, len(s) - 2, 3):
                aa = CODON.get(s[j:j + 3].upper(), 'X')
                if aa == '*':
                    if len(cur) > len(best):
                        best = ''.join(cur)
                    cur = []
                else:
                    cur.append(aa)
            if len(cur) > len(best):
                best = ''.join(cur)
    return best


def annot_locus_seqs(gff, fasta, chrom):
    """guided nodes: spliced sequence of each ANNOTATED gene, keyed `chrom:genestart-geneend`.

    ⚠ Exons of the longest transcript, NOT the annotated CDS -- the whole point is that 0 of 457 chr16
    pseudogenes have a CDS, so the ORF must be found de novo on the spliced transcript either way.
    """
    import pysam
    fa = pysam.FastaFile(fasta)
    gspan, tx, tx_gene = {}, collections.defaultdict(list), {}
    for ln in open(gff):
        if ln.startswith('#'):
            continue
        f = ln.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom:
            continue
        if f[2] in ('gene', 'pseudogene', 'ncRNA_gene'):
            m = re.search(r'ID=([^;]+)', f[8])
            if m:
                gspan[m.group(1)] = (int(f[3]), int(f[4]))
        elif f[2] in ('mRNA', 'transcript'):
            i = re.search(r'ID=([^;]+)', f[8]); p = re.search(r'Parent=([^;]+)', f[8])
            if i and p:
                tx_gene[i.group(1)] = p.group(1)
        elif f[2] == 'exon':
            p = re.search(r'Parent=([^;]+)', f[8])
            if p:
                tx[p.group(1)].append((int(f[3]), int(f[4])))
    best = {}
    for t, ex in tx.items():
        g = tx_gene.get(t)
        if g not in gspan:
            continue
        if g not in best or sum(b - a for a, b in ex) > sum(b - a for a, b in best[g]):
            best[g] = sorted(ex)
    out = {}
    for g, ex in best.items():
        s0, e0 = gspan[g]
        out[f'{chrom}:{s0}-{e0}'] = ''.join(fa.fetch(chrom, a - 1, b) for a, b in ex)
        out[f'{chrom}:{s0 - 1}-{e0}'] = out[f'{chrom}:{s0}-{e0}']
    return out


def locus_seqs(gff3, fasta, chrom):
    import pysam
    fa = pysam.FastaFile(fasta)
    gene, ex = {}, collections.defaultdict(list)
    for ln in open(gff3):
        if ln.startswith('#'):
            continue
        f = ln.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom:
            continue
        if f[2] == 'gene':
            g = re.search(r'ID=gene-(\S+?)(?:;|$)', f[8]).group(1)
            gene[g] = (f[0], int(f[3]), int(f[4]))
        elif f[2] == 'exon':
            ex[re.search(r'Parent=gene-(\S+?);', f[8]).group(1)].append((int(f[3]), int(f[4])))
    out = {}
    for g, (c, s, e) in gene.items():
        blocks = sorted(ex.get(g, []))
        if blocks:
            out[f'{c}:{s}-{e}'] = ''.join(fa.fetch(c, a - 1, b) for a, b in blocks)
    return out


def shuffle_codons(seq, rng):
    cods = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]
    rng.shuffle(cods)
    return ''.join(cods)


def protein_edges(prot, workdir, tag, min_aa, min_ident, min_cov, threads):
    faa = os.path.join(workdir, f'{tag}.faa')
    with open(faa, 'w') as fh:
        for k, p in prot.items():
            fh.write(f'>{k}\n{p}\n')
    paf = os.path.join(workdir, f'{tag}.paf')
    # blastp, matching r906's own tool choice so the arms are comparable
    bb = os.environ.get('BLAST_BIN', '/home/juanfra/miniforge3/envs/blast/bin')
    subprocess.run([f'{bb}/makeblastdb', '-in', faa, '-dbtype', 'prot', '-out', faa + '.db'],
                   check=True, stdout=subprocess.DEVNULL)
    with open(paf, 'w') as o:
        subprocess.run([f'{bb}/blastp', '-query', faa, '-db', faa + '.db', '-outfmt',
                        '6 qseqid sseqid pident length qlen slen', '-evalue', '1e-5',
                        '-num_threads', str(threads), '-max_target_seqs', '500'],
                       stdout=o, check=True)
    E = set()
    for ln in open(paf):
        f = ln.rstrip('\n').split('\t')
        if len(f) < 6 or f[0] == f[1]:
            continue
        ident = float(f[2]) / 100.0
        alen, ql, sl = int(f[3]), int(f[4]), int(f[5])
        if ident >= min_ident and alen >= min_aa and alen / max(ql, sl) >= min_cov:
            E.add((f[0], f[1]) if f[0] < f[1] else (f[1], f[0]))
    return E


def main():
    ap = argparse.ArgumentParser()
    for a in ('--gff3', '--fasta', '--graph', '--gff', '--soto', '--workdir'):
        ap.add_argument(a, required=True)
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--min-aa', type=int, default=100)
    ap.add_argument('--min-ident', type=float, default=0.70)
    ap.add_argument('--min-cov', type=float, default=0.30)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--annot-mode', action='store_true',
                    help='guided nodes: build spliced sequence from an ANNOTATION gff (gene->mRNA->exon)')
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)

    N0 = set(); nodes = set()
    for ln in open(a.graph):
        f = ln.rstrip('\n').split('\t')
        if len(f) == 3 and f[0] != f[1]:
            N0.add((f[0], f[1]) if f[0] < f[1] else (f[1], f[0])); nodes.add(f[0]); nodes.add(f[1])
        elif f and f[0]:
            nodes.add(f[0])

    seqs = (annot_locus_seqs(a.gff3, a.fasta, a.chrom) if a.annot_mode
            else locus_seqs(a.gff3, a.fasta, a.chrom))
    prot = {k: longest_orf(v) for k, v in seqs.items() if k in nodes}
    prot = {k: v for k, v in prot.items() if len(v) >= a.min_aa}
    rng = random.Random(0)
    null = {k: longest_orf(shuffle_codons(seqs[k], rng)) for k in prot}
    null = {k: v for k, v in null.items() if len(v) >= a.min_aa}
    print(f"{a.chrom}: graph nodes {len(nodes)} | ORF >= {a.min_aa}aa: {len(prot)} "
          f"({100*len(prot)/max(len(nodes),1):.1f}%) | null ORFs {len(null)}")

    P = protein_edges(prot, a.workdir, f'{a.chrom}_p', a.min_aa, a.min_ident, a.min_cov, a.threads)
    Z = protein_edges(null, a.workdir, f'{a.chrom}_z', a.min_aa, a.min_ident, a.min_cov, a.threads)
    print(f"  protein edges {len(P)}   NULL (codon-shuffled) edges {len(Z)}"
          f"   -> null/protein = {100*len(Z)/max(len(P),1):.1f}%")

    # node -> gene, max overlap; Soto family truth
    genes = []
    for ln in open(a.gff):
        if ln.startswith('#'):
            continue
        f = ln.split('\t')
        if len(f) < 9 or f[0] != a.chrom or f[2] not in ('gene', 'pseudogene', 'ncRNA_gene'):
            continue
        m = re.search(r'Name=([^;]+)', f[8])
        if m:
            genes.append((int(f[3]) - 1, int(f[4]), m.group(1)))
    genes.sort()

    def gname(n):
        m = re.match(r'(\S+):(\d+)-(\d+)$', n)
        if not m:
            return None
        s, e = int(m.group(2)), int(m.group(3)); best = None
        for gs, ge, g in genes:
            if ge < s:
                continue
            if gs > e:
                break
            ov = min(e, ge) - max(s, gs)
            if ov > 0 and (best is None or ov > best[0]):
                best = (ov, g)
        return best[1] if best else None

    fam = {}
    for r in csv.DictReader(open(a.soto), delimiter='\t'):
        f, g = (r.get('Family ID') or '').strip(), (r.get('Gene Name') or '').strip()
        if f and f != 'N/A' and g:
            fam.setdefault(g, f)
    ng = {n: gname(n) for n in nodes}
    lab = {n: fam[g] for n, g in ng.items() if g in fam}
    byfam = collections.defaultdict(set)
    for n, f in lab.items():
        byfam[f].add(n)
    truth = {f: v for f, v in byfam.items() if len(v) >= 2}
    tp = set()
    for v in truth.values():
        vs = sorted(v)
        for i in range(len(vs)):
            for j in range(i + 1, len(vs)):
                tp.add((vs[i], vs[j]))
    def rep(E, name):
        cov = len(E & tp) / max(len(tp), 1)
        cross = sum(1 for x, y in E if x in lab and y in lab and lab[x] != lab[y])
        print(f"  {name:<22} edges {len(E):>6}  within-family pair coverage {100*cov:>5.1f}%  "
              f"cross-family edges {cross:>4}")
    print(f"  truth: {len(truth)} Soto families, {len(tp)} within-family pairs")
    rep(N0, 'N0 nucleotide'); rep(N0 | P, 'P  = N0 + protein'); rep(N0 | Z, 'NULL = N0 + shuffled')


if __name__ == '__main__':
    main()
