#!/usr/bin/env python3
"""LIFTMODE-GW3-RBH protein judge (docs/LEDGER_O1.md r1456), re-implemented for the cluster (2026-10-29): protein = translation of each
gene's longest annotated CDS (RefSeq GFF for gorilla via DENOVO_GGO_GFF; CAT GFF for human via DENOVO_CAT_GFF); `mmseqs easy-search`
all-vs-all (-s 7.5, --max-seqs 2000, e-value 1e-3); for a gene pair (A,B): MUTUAL HIT = A's hit list has B and B's has A.

  protjudge.py prep  SPECIES FASTA OUTDIR GENES_TXT      proteins of the listed gene ids (one per line) -> OUTDIR/prot.fa, OUTDIR/hits.m8
  protjudge.py pairs SPECIES BASE_COPIES ARM_COPIES OUTDIR    gene set + base/added coding pairs of an arm vs its base
  protjudge.py rate  OUTDIR PAIRS_TSV [PAIRS_TSV...]     each file: gene_a<TAB>gene_b (ids as in famtruth); prints n coding, mutual, rate
"""
import collections, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'guided_eval'))
import famtruth as F  # noqa: E402

CODON = {}
_b = 'TCAG'
_aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODON[a + b + c] = _aa[16 * i + 4 * j + k]


def translate(s):
    return ''.join(CODON.get(s[i:i + 3], 'X') for i in range(0, len(s) - 2, 3))


def longest_cds(species, wanted):
    """gene id -> [(chrom, s0, e, strand, phase)] of its longest CDS (by coding bp)."""
    path = F.GGO_GFF if species == 'gorilla' else F.CAT_GFF
    per_tx = collections.defaultdict(list); tx_gene = {}
    for line in open(path):
        f = line.split('\t', 8)
        if len(f) < 9 or f[2] != 'CDS':
            continue
        a = dict(kv.split('=', 1) for kv in f[8].rstrip('\n;').split(';') if '=' in kv)
        g = a.get('gene_id') if species == 'human' else None
        if species == 'gorilla':
            g = 'gene-' + a['gene'] if 'gene' in a else None
        if g not in wanted:
            continue
        per_tx[(g, a['Parent'])].append((f[0], int(f[3]) - 1, int(f[4]), f[6], f[7]))
    best = {}
    for (g, t), segs in per_tx.items():
        n = sum(e - s for _, s, e, _, _ in segs)
        if n > best.get(g, (0, None))[0]:
            best[g] = (n, sorted(segs, key=lambda x: x[1]))
    return {g: v[1] for g, v in best.items()}


def prep(species, fasta, outdir, genes_txt):
    os.makedirs(outdir, exist_ok=True)
    wanted = {l.strip() for l in open(genes_txt) if l.strip()}
    cds = longest_cds(species, wanted)
    regs = sorted({f'{c}:{s + 1}-{e}' for segs in cds.values() for c, s, e, _, _ in segs})
    import cdsseed as C  # noqa
    C.FA = fasta
    seqs = {}
    for i in range(0, len(regs), 20000):
        seqs.update(C.faidx(regs[i:i + 20000], f'{outdir}/cds{i}'))
    with open(f'{outdir}/prot.fa', 'w') as w:
        for g, segs in cds.items():
            s = ''.join(seqs[f'{c}:{a + 1}-{b}'] for c, a, b, _, _ in segs)
            if segs[0][3] == '-':
                s = s.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]
                ph = segs[-1][4]
            else:
                ph = segs[0][4]
            s = s[int(ph) if ph in '012' else 0:]
            w.write(f'>{g}\n{translate(s).rstrip("*")}\n')
    print(f'genes wanted {len(wanted)}, with a CDS {len(cds)}', file=sys.stderr)
    subprocess.run(['mmseqs', 'easy-search', f'{outdir}/prot.fa', f'{outdir}/prot.fa', f'{outdir}/hits.m8', f'{outdir}/tmp',
                    '-s', '7.5', '--max-seqs', '2000', '-e', '1e-3', '--threads', os.environ.get('THREADS', '8')],
                   check=True, stdout=subprocess.DEVNULL)


def rate(outdir, pair_files):
    prot = {l[1:].strip() for l in open(f'{outdir}/prot.fa') if l[0] == '>'}
    hits = collections.defaultdict(set)
    for line in open(f'{outdir}/hits.m8'):
        q, t = line.split('\t', 2)[:2]
        if q != t:
            hits[q].add(t)
    for p in pair_files:
        n = m = nocds = 0
        for line in open(p):
            a, b = line.rstrip('\n').split('\t')[:2]
            if a not in prot or b not in prot:
                nocds += 1
                continue
            n += 1
            m += (b in hits[a]) and (a in hits[b])
        print(f'{p}\tboth coding {n}\tmutual hit {m}\trate {m / max(1, n):.3f}\t(no protein: {nocds})')


def pairs(species, base, arm, outdir):
    """GW3-RBH gene set: genes of the base pairs and of the pairs the arm adds (copies mapped to genes as famtruth does; both
    protein-coding, LOC genes included). Writes OUTDIR/{base,added}.pairs.tsv and OUTDIR/genes.txt."""
    import itertools
    os.makedirs(outdir, exist_ok=True)
    G = F.load_cat() if species == 'human' else F.load_ggo(F.load_cat())
    idx = F.exon_index(G)

    def gene_pairs(path):
        fams, _, _ = F.load_catalog(path, G, idx)
        return {(a, b) for v in fams.values() for a, b in itertools.combinations(sorted(v), 2)
                if G[a]['biotype'] == 'protein_coding' and G[b]['biotype'] == 'protein_coding'}
    pb, pa = gene_pairs(base), gene_pairs(arm)
    added = pa - pb
    for name, ps in (('base', pb), ('added', added)):
        with open(f'{outdir}/{name}.pairs.tsv', 'w') as w:
            w.write(''.join(f'{a}\t{b}\n' for a, b in sorted(ps)))
    genes = {g for p in pb | added for g in p}
    with open(f'{outdir}/genes.txt', 'w') as w:
        w.write(''.join(g + '\n' for g in sorted(genes)))
    print(f'base coding pairs {len(pb)}, arm {len(pa)}, added {len(added)}, genes {len(genes)}', file=sys.stderr)


if __name__ == '__main__':
    if sys.argv[1] == 'pairs':
        pairs(*sys.argv[2:6])
    elif sys.argv[1] == 'prep':
        prep(*sys.argv[2:6])
    else:
        rate(sys.argv[2], sys.argv[3:])
