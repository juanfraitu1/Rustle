#!/usr/bin/env python3
"""Edge-construction probes (§6y0-§6y7), consolidated from five one-off scripts so the register rows they
produced (r1022-r1036) keep a runnable, ANCHORED generator without five files. Each subcommand is the original
script verbatim (only `main` renamed and one colliding helper prefixed), EXCEPT one determinism fix: the
`poset` block iterates its node set in sorted order, because the original iterated a Python set and its P2
antichain numbers varied with PYTHONHASHSEED (F .120 vs .135 on chr16); P1/P3/M are unaffected. The
original docstrings follow.

  edge_probes.py poset ...               was bench/containment_poset.py      (r1022/r1023)
  edge_probes.py intron-chain ...        was bench/intron_chain_edges.py     (r1025/r1026)
  edge_probes.py protein-denovo ...      was bench/protein_sensitive_mode.py (r1027)
  edge_probes.py protein-false-merge ... was bench/protein_false_merge.py    (r1028-r1031)
  edge_probes.py cov-shorter ...         was bench/cov_shorter_conjunct.py   (r1035/r1036)
"""
import sys
import argparse
import bisect
import collections
import itertools
import re
import sys
sys.path.insert(0, 'bench')
import mcl_port
import csv
import math
import argparse, collections, csv, random, re, subprocess, os
import argparse, collections, csv, os, re, subprocess, sys
import protein_edge_gap as PEG
import argparse, bisect, collections, csv, itertools, re, sys


# ================================================================================================
# poset  (was bench/containment_poset.py)  renamed: main->main_containment_poset
# ================================================================================================
"""Is a POSET a better family object than a partition of an undirected graph?
Per `docs/PREREG_containment_poset_2026-09-22.md` (§6y0, md5 `4c0be9cc`).

r1020 measured that containment between loci is a strict partial order (antisymmetry 0 violations,
transitivity 92.9%, longest chain 38, AUC 0.943 for same-family and NOT a size artefact). r1021 proposed
that §6s9's inexpressible fusion, §6u8's dropped cover and §6t7's asymmetry are one limitation: the output
is a partition of an undirected graph while the data is a poset. **This tests whether that scores.**

⚠ THE EDGE SET IS HELD FIXED. Both arms consume the shipped `--dump-graph` output, which already passed
every shipped conjunct; direction is then read off the PAF. No edge is added or removed by either arm, so a
difference is attributable to the structure alone.
⚠ The comparator is `mcl_port` MCL ON THE SAME GRAPH, never the shipped Rust F (register 917).
⚠ PRIMARY metric is PAIRWISE, because bipartite matching is not well defined for a cover and penalises one
by construction (§6u8's own warning).

Usage:
  containment_poset.py --graph X.graph.tsv --paf X.paf --gff chrN.genes.gff --truth T.tsv --chrom chrN"""
def merged_len(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return sum(e - s for s, e in out)

def read_graph(path):
    """the shipped pre-MCL graph: `u<TAB>v<TAB>w`, plus a self-row per node."""
    adj = collections.defaultdict(dict)
    nodes = set()
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) == 3 and f[0] != f[1]:
            a, b, w = f[0], f[1], float(f[2])
            adj[a][b] = max(adj[a].get(b, 0.0), w)
            adj[b][a] = max(adj[b].get(a, 0.0), w)
            nodes.add(a); nodes.add(b)
        elif f and f[0]:
            nodes.add(f[0])
    return adj, nodes

def directed(paf, keep, C):
    """x -> y ('x is contained in y') for edges already in the shipped graph."""
    qlen = {}
    acc = collections.defaultdict(lambda: [[], []])
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 11:
            continue
        q, t = f[0], f[5]
        qlen[q] = int(f[1]); qlen[t] = int(f[6])
        if q == t:
            continue
        key = (q, t) if q <= t else (t, q)
        if key not in keep:
            continue
        a = (int(f[2]), int(f[3])); b = (int(f[7]), int(f[8]))
        e = acc[key]
        if key[0] == q:
            e[0].append(a); e[1].append(b)
        else:
            e[0].append(b); e[1].append(a)
    below = collections.defaultdict(set)
    for (x, y), (ix, iy) in acc.items():
        if not ix:
            continue
        cx = merged_len(ix) / max(qlen[x], 1)
        cy = merged_len(iy) / max(qlen[y], 1)
        if cx >= C and cx > cy:
            below[x].add(y)
        elif cy >= C and cy > cx:
            below[y].add(x)
    return below

def transitive_closure(below):
    out = {k: set(v) for k, v in below.items()}
    changed = True
    while changed:
        changed = False
        for x in list(out):
            add = set()
            for y in out[x]:
                add |= out.get(y, set())
            add.discard(x)
            if not add <= out[x]:
                out[x] |= add; changed = True
    return out

def comparability_components(below, nodes):
    adj = collections.defaultdict(set)
    for x, ys in below.items():
        for y in ys:
            adj[x].add(y); adj[y].add(x)
    seen = set(); comps = []
    for n in sorted(nodes):   # ⚠ deterministic: the original iterated a SET, making P2 antichains hash-seed dependent
        if n in seen:
            continue
        stack = [n]; seen.add(n); c = []
        while stack:
            u = stack.pop(); c.append(u)
            for w in adj.get(u, ()):
                if w not in seen:
                    seen.add(w); stack.append(w)
        comps.append(c)
    return comps

def p1_downsets(below, nodes):
    """principal down-set of each MAXIMAL element -- the natural cover."""
    above = collections.defaultdict(set)
    for x, ys in below.items():
        for y in ys:
            above[y].add(x)
    maximal = [n for n in nodes if not below.get(n)]
    fams = {}
    for i, m in enumerate(maximal):
        seen = {m}; stack = [m]
        while stack:
            u = stack.pop()
            for w in above.get(u, ()):
                if w not in seen:
                    seen.add(w); stack.append(w)
        if len(seen) >= 2:
            fams[f'P1_{i}'] = sorted(seen)
    return fams

def p2_antichains(below, nodes):
    """maximal antichains, greedily, inside each comparability component."""
    comp = comparability_components(below, nodes)
    fams = {}; k = 0
    for c in comp:
        if len(c) < 2:
            continue
        rel = {x: below.get(x, set()) for x in c}
        rest = sorted(c, key=lambda n: -len(rel.get(n, ())))
        while rest:
            chain = []
            for n in rest:
                if all(n not in rel.get(m, ()) and m not in rel.get(n, ()) for m in chain):
                    chain.append(n)
            if len(chain) >= 2:
                fams[f'P2_{k}'] = sorted(chain); k += 1
            rest = [n for n in rest if n not in set(chain)]
            if not chain:
                break
    return fams

def pairwise(truth, pred):
    """cover-compatible: a pair is together iff it co-occurs in ANY family.

    ⚠ Restricted to the TRUTH'S OWN UNIVERSE (nodes the truth labels). Scoring every predicted pair
    against a truth that labels only part of the node set deflates precision by the universe mismatch,
    not by the method. Conditioning on the TRUTH is correct; conditioning on the PREDICTION is register
    770's trap and is not what this does.
    """
    universe = set().union(*truth.values()) if truth else set()

    def pairs(d):
        s = set()
        for members in d.values():
            m = sorted(set(members) & universe)
            for a, b in itertools.combinations(m, 2):
                s.add((a, b))
        return s
    T, P = pairs(truth), pairs(pred)
    if not T:
        return (float('nan'),) * 3
    tp = len(T & P)
    prec = tp / len(P) if P else 0.0
    rec = tp / len(T)
    f = 0.0 if prec + rec == 0 else 2 * prec * rec / (prec + rec)
    return prec, rec, f

def main_containment_poset():
    ap = argparse.ArgumentParser()
    for a in ('--graph', '--paf', '--gff', '--truth'):
        ap.add_argument(a, required=True)
    ap.add_argument('--chrom', default='chr16')
    ap.add_argument('--label', default='arm')
    ap.add_argument('--sweep', default='0.50,0.60,0.70,0.80,0.90')
    a = ap.parse_args()

    adj, nodes = read_graph(a.graph)
    keep = {(u, v) if u <= v else (v, u) for u in adj for v in adj[u]}

    # locus -> gene name, max overlap (the same resolver every other scorer uses)
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
    gs_list = [g[0] for g in genes]

    def gname(node):
        m = re.match(r'(\S+):(\d+)-(\d+)$', node)
        if not m or m.group(1) != a.chrom:
            return None
        s, e = int(m.group(2)), int(m.group(3))
        best = None
        for idx in range(max(0, bisect.bisect_left(gs_list, s) - 40), len(genes)):
            g0, g1, g = genes[idx]
            if g0 > e:
                break
            ov = min(e, g1) - max(s, g0)
            if ov > 0 and (best is None or ov > best[0]):
                best = (ov, g)
        return best[1] if best else None

    lab = {}
    for ln in open(a.truth):
        f = ln.rstrip('\n').split('\t')
        if len(f) >= 2 and f[0] != 'Gene Name':
            lab.setdefault(f[0], set()).add(f[1])
    truth = collections.defaultdict(set)
    node_gene = {n: gname(n) for n in nodes}
    for n, g in node_gene.items():
        for fam in lab.get(g, ()):
            truth[fam].add(n)
    truth = {k: v for k, v in truth.items() if len(v) >= 2}

    mcl_f = mcl_port.mcl({(u, v): w for u in adj for v, w in adj[u].items() if u < v}, inflation=2.8)
    M = {f'M_{i}': c for i, c in enumerate(mcl_f) if len(set(c)) >= 2}
    big = lambda d: max((len(set(v)) for v in d.values()), default=0)
    print(f"{a.label}: nodes {len(nodes)} | truth {len(truth)} families | "
          f"M: {len(M)} fams, largest {big(M)}")
    p, r, f = pairwise(truth, M)
    print(f"  {'M (mcl_port)':<26} pairwise P {p:.3f} R {r:.3f} F {f:.3f}   fams {len(M):>4} largest {big(M):>4}")

    for C in [float(x) for x in a.sweep.split(',')]:
        raw = directed(a.paf, keep, C)
        for tag, rel in (('raw', raw), ('tclosed', transitive_closure(raw))):
            objs = {'P1 downsets': p1_downsets(rel, nodes), 'P2 antichains': p2_antichains(rel, nodes),
                    'P3 components (null)': {f'P3_{i}': c for i, c in
                                             enumerate(comparability_components(rel, nodes)) if len(c) >= 2}}
            for name, fams in objs.items():
                p, r, f = pairwise(truth, fams)
                print(f"  C={C:.2f} {tag:<8} {name:<21} P {p:.3f} R {r:.3f} F {f:.3f}   "
                      f"fams {len(fams):>4} largest {big(fams):>4}")

# ================================================================================================
# intron-chain  (was bench/intron_chain_edges.py)  renamed: main->main_intron_chain_edges
# ================================================================================================
"""Can an ALIGNMENT-FREE intron-chain certificate create edges minimap2 never proposes?
Per `docs/PREREG_trie_edge_construction_2026-09-22.md` (§6y2, md5 `619a1999`).

§6u4 split edge-construction loss into GATE REJECTION (22.7%, already addressed by `--min-cov-shorter`)
and NO ALIGNMENT (26.1%). No alignment-derived signal can reach the second by construction; this measures
whether intron structure can.

Rule: two genes are joined iff they share a 3-intron length shingle at 5% tolerance (log-binned) --
r1024's operating point (precision 0.979 / recall 13.5% overall, but 99.3% redundant with the aligner).
⚠ A gene with < 3 exons has no 3-shingle and is UNREACHABLE by this rule; that ceiling is reported."""
def chains(gff, chrom):
    """gene -> intron-length chain of its longest annotated transcript, oriented 5'->3'."""
    tx = collections.defaultdict(list)
    tx_gene = {}
    for ln in open(gff):
        if ln.startswith('#'):
            continue
        f = ln.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom:
            continue
        if f[2] == 'exon':
            p = re.search(r'Parent=([^;]+)', f[8])
            if p:
                tx[p.group(1)].append((int(f[3]), int(f[4]), f[6]))
        elif f[2] in ('mRNA', 'transcript'):
            i = re.search(r'ID=([^;]+)', f[8])
            g = re.search(r'gene=([^;]+)', f[8]) or re.search(r'Parent=gene-([^;]+)', f[8])
            if i and g:
                tx_gene[i.group(1)] = g.group(1)
    out = {}
    for t, ex in tx.items():
        g = tx_gene.get(t)
        if not g or len(ex) < 3:
            continue
        ex.sort()
        iv = [ex[i + 1][0] - ex[i][1] - 1 for i in range(len(ex) - 1)]
        if ex[0][2] == '-':
            iv = iv[::-1]
        if g not in out or len(iv) > len(out[g]):
            out[g] = iv
    return out

def spans(gff, chrom):
    out = {}
    for ln in open(gff):
        if ln.startswith('#'):
            continue
        f = ln.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] not in ('gene', 'pseudogene', 'ncRNA_gene'):
            continue
        m = re.search(r'Name=([^;]+)', f[8])
        if m:
            out[m.group(1)] = (int(f[3]) - 1, int(f[4]))
    return out

def aligned_pairs(paf, gspan, chrom):
    """gene pairs minimap2 proposes ANY alignment for (node header = chrom:start-end)."""
    node = {}
    for g, (s, e) in gspan.items():
        node[f'{chrom}:{s}-{e}'] = g
        node[f'{chrom}:{s + 1}-{e}'] = g
    out = set()
    for ln in open(paf):
        f = ln.split('\t', 6)
        a, b = node.get(f[0]), node.get(f[5])
        if a and b and a != b:
            out.add((a, b) if a < b else (b, a))
    return out

def main_intron_chain_edges():
    ap = argparse.ArgumentParser()
    for a in ('--gff', '--paf', '--truth'):
        ap.add_argument(a, required=True)
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--tol', type=float, default=0.05)
    ap.add_argument('--k', type=int, default=3)
    a = ap.parse_args()

    ch = chains(a.gff, a.chrom)
    gsp = spans(a.gff, a.chrom)
    aln = aligned_pairs(a.paf, gsp, a.chrom)
    lab = {r['Gene Name']: r['Family ID']
           for r in csv.DictReader(open(a.truth), delimiter='\t') if r.get('Family ID')}

    # every referee-labelled gene, whether or not it has a usable chain -- the ceiling is part of the result
    scored = sorted(g for g in lab if g in gsp)
    withchain = [g for g in scored if g in ch and len(ch[g]) >= a.k]
    b = lambda x: int(math.log(max(x, 1)) / math.log(1 + a.tol))
    sh = {g: {tuple(b(ch[g][i + j]) for j in range(a.k)) for i in range(len(ch[g]) - a.k + 1)}
          for g in withchain}

    truth_pairs = {(x, y) for x, y in itertools.combinations(scored, 2) if lab[x] == lab[y]}
    noaln_truth = {p for p in truth_pairs if p not in aln}
    rec = {p for p in noaln_truth if p[0] in sh and p[1] in sh and sh[p[0]] & sh[p[1]]}
    # false positives among UNALIGNED pairs: different family, no alignment, certificate fires
    fp = 0
    for x, y in itertools.combinations(withchain, 2):
        if lab[x] == lab[y]:
            continue
        if (x, y) in aln:
            continue
        if sh[x] & sh[y]:
            fp += 1
    prec = len(rec) / (len(rec) + fp) if (len(rec) + fp) else float('nan')
    reachable = {p for p in noaln_truth if p[0] in sh and p[1] in sh}
    print(f"{a.chrom}: referee genes {len(scored)} (with a >={a.k}-intron chain: {len(withchain)})")
    print(f"  truth same-family pairs            {len(truth_pairs)}")
    print(f"  of which NO ALIGNMENT (the target) {len(noaln_truth)}  "
          f"[reachable by the rule: {len(reachable)} = {100*len(reachable)/max(len(noaln_truth),1):.1f}%]")
    print(f"  ⭐ recovered by the certificate     {len(rec)} = {100*len(rec)/max(len(noaln_truth),1):.1f}% of target")
    print(f"  false joins among UNALIGNED pairs  {fp}")
    print(f"  ⭐ precision on unaligned pairs     {prec:.3f}")

# ================================================================================================
# protein-denovo  (was bench/protein_sensitive_mode.py)  renamed: main->main_protein_sensitive_mode
# ================================================================================================
"""SENSITIVE MODE: protein edges from translated DE NOVO locus representatives.
Per `docs/PREREG_protein_sensitive_mode_2026-09-22.md` (§6y3, md5 `e9278fdd`).

r906 measured protein edges from ANNOTATED CDS as ⚠PARTIAL (+13.8 pts pair coverage, 0 cross-family
held-out) and left them unadopted. This translates the de novo locus instead, which is the only route that
reaches pseudogenes: 0 of 457 chr16 pseudogenes have an annotated CDS.

⚠⚠ TRUTH IS SOTO, NOT THE PROTEIN REFEREE. The referee is itself built from translated CDS clustered by
protein homology, so scoring protein-derived edges against it is circular by construction.
⚠ A six-frame longest ORF is non-zero for ANY sequence, so ORF yield is not evidence. The NULL arm
(codon-shuffled ORFs) must produce ~nothing; if it fires, the signal is chance similarity."""
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

def main_protein_sensitive_mode():
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

# ================================================================================================
# protein-false-merge  (was bench/protein_false_merge.py)  renamed: main->main_protein_false_merge
# ================================================================================================
"""The FALSE-MERGE measurement that blocks r906's adoption.
Per `docs/PREREG_protein_false_merge_2026-09-22.md` (§6y4, md5 `ca4fa685`).

r906 measured protein edges as ⚠PARTIAL (+13.8 pts held-out pair coverage, 0 cross-family) and recorded
its own gap: "Precision measured only over Soto-labelled genes, so this is NOT a genome-wide false-merge
rate." Every protein edge with an UNLABELLED endpoint is unmeasured, and a sensitive mode's risk lives
exactly there.

§6ko's rule is reused verbatim from `bench/protein_edge_gap.py` -- nothing is re-tuned.
Classification of every PROTEIN-ONLY edge (no nucleotide edge in the shipped graph):
  TRUE     both endpoints Soto-labelled, same family
  FALSE    both endpoints Soto-labelled, different families
  UNKNOWN  at least one endpoint unlabelled  <- r906's blind spot
and each UNKNOWN split by STRUCTURAL corroboration (label-free): does ANY nucleotide PAF record exist for
the pair, even one that failed the gate?"""
def protein_edges_at(faa, out, plen, threads, floor):
    """§6ko's rule with the coverage FLOOR as a parameter (§6y5). Identical to
    `protein_edge_gap.protein_edges` at floor=0.30; the blastp cache is reused across floors."""
    import collections as _c
    PEG.protein_edges(faa, out, plen, threads)          # populates <out>.blastp.tsv
    hs = _c.defaultdict(list)
    for line in open(out + '.blastp.tsv'):
        q, sj, nid, ln, q0, q1, s0, s1, bits = line.rstrip('\n').split('\t')
        if q != sj:
            hs[tuple(sorted((q, sj)))].append((float(bits), int(q0), int(q1), q, sj))
    ed = set()
    for (a, b), v in hs.items():
        longer = max(plen.get(a, 0), plen.get(b, 0))
        if not longer:
            continue
        taken = []
        for bits, q0, q1, q, sj in sorted(v, reverse=True):
            if plen.get(q, 0) != longer:
                continue
            lo, hi = min(q0, q1), max(q0, q1)
            if all(hi < t0 or lo > t1 for t0, t1 in taken):
                taken.append((lo, hi))
        if sum(hi - lo + 1 for lo, hi in taken) / longer >= floor:
            ed.add((a, b))
    return ed

def main_protein_false_merge():
    ap = argparse.ArgumentParser()
    for x in ('--gff', '--genome', '--paf', '--chrom', '--soto', '--graph', '--out'):
        ap.add_argument(x, required=True)
    ap.add_argument('--threads', default='4')
    ap.add_argument('--floors', default='0.30,0.50,0.70,0.80,0.90')
    a = ap.parse_args()
    import pysam

    cds = PEG.longest_cds(a.gff, a.chrom)
    fa = pysam.FastaFile(a.genome)
    faa = a.out + '.proteins.faa'
    plen = {}
    with open(faa, 'w') as fh:
        for g, (st, segs) in cds.items():
            p = PEG.translate(fa, a.chrom, st, segs)
            if len(p) >= 10:
                plen[g] = len(p); fh.write(f'>{g}\n{p}\n')
    floors = [float(x) for x in a.floors.split(',')]

    # symbol <-> node, from the GFF gene spans (PAF/graph names are chrom:start-end)
    gene_at, span_of = {}, {}
    for line in open(a.gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != a.chrom or f[2] not in ('gene', 'pseudogene', 'ncRNA_gene'):
            continue
        n = re.search(r'Name=([^;]+)', f[8])
        if n:
            s, e = int(f[3]), int(f[4])
            gene_at[f'{a.chrom}:{s}-{e}'] = n.group(1)
            gene_at[f'{a.chrom}:{s-1}-{e}'] = n.group(1)
            span_of[n.group(1)] = (s, e)

    nuc = set()
    for line in open(a.graph):
        f = line.rstrip('\n').split('\t')
        if len(f) == 3 and f[0] != f[1]:
            x, y = gene_at.get(f[0]), gene_at.get(f[1])
            if x and y and x != y:
                nuc.add((x, y) if x < y else (y, x))
    # ANY paf record at all, even one the gate rejected -- label-free corroboration
    anyaln = set()
    for line in open(a.paf):
        f = line.split('\t', 6)
        x, y = gene_at.get(f[0]), gene_at.get(f[5])
        if x and y and x != y:
            anyaln.add((x, y) if x < y else (y, x))

    fam = {}
    for r in csv.DictReader(open(a.soto), delimiter='\t'):
        fi, g = (r.get('Family ID') or '').strip(), (r.get('Gene Name') or '').strip()
        if fi and fi != 'N/A' and g:
            fam.setdefault(g, fi)

    def root(sym):
        r = re.sub(r'\d+[A-Z]*$', '', sym)
        return r if len(r) >= 3 else sym

    print(f"{a.chrom}:  floor | prot edges | ⭐PROT-ONLY |  TRUE FALSE  corrob   bare | "
          f"⭐pessimistic FM | cross-root% of bare")
    for fl in floors:
        P = protein_edges_at(faa, a.out, plen, a.threads, fl)
        only = {e for e in P if e not in nuc}
        T = F = unk_corr = unk_bare = 0
        cross = 0
        for x, y in sorted(only):
            if x in fam and y in fam:
                if fam[x] == fam[y]:
                    T += 1
                else:
                    F += 1
            elif (x, y) in anyaln:
                unk_corr += 1
            else:
                unk_bare += 1
                if root(x) != root(y):
                    cross += 1
        n = len(only)
        pes = (F + unk_bare) / n if n else float('nan')
        cr = cross / unk_bare if unk_bare else float('nan')
        print(f"{'':>7} {fl:>5.2f} | {len(P):>10} | {n:>10} | {T:>5} {F:>5} {unk_corr:>7} {unk_bare:>6} | "
              f"{pes:>15.4f} | {100*cr:>17.1f}%")

# ================================================================================================
# cov-shorter  (was bench/cov_shorter_conjunct.py)  renamed: merged_len->merged_len_cov_shorter_conjunct, main->main_cov_shorter_conjunct
# ================================================================================================
"""Does an alignment-coverage floor on the SHORTER side add anything the exon-sharing floor does not?
Per `docs/PREREG_cov_shorter_conjunct_2026-09-22.md` (§6y7, md5 `b36934b2`).

§6y6/r1032: `RUSTLE_ER_COVERAGE_LONGER_FLOOR` is unreachable from the shipped catalog's driver. Its
faithful analogue on the CURRENT definition is a CONJUNCT on the shorter side -- but the current rule is
ALREADY two-sided (`cov_longer >= 0.30` on the longer, `min_shared_exon_frac >= 0.60` on the SMALLER
gene's exonic length), so this tests REDUNDANCY, not novelty.

⚠ Opposite direction to `--min-cov-shorter` (§6x4), which is an OR-escape and LOOSENS. This ANDs.
⚠ Comparator is mcl_port on the UNFILTERED graph, never the shipped Rust F (r917)."""
def merged_len_cov_shorter_conjunct(iv):
    iv = sorted(iv); out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return sum(e - s for s, e in out)

def exonic_len(gff, chrom):
    """per-gene exon-union length, keyed chrom:start-end like the graph/PAF headers."""
    tx, tx_gene, gspan = collections.defaultdict(list), {}, {}
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
    per = collections.defaultdict(list)
    for t, ex in tx.items():
        g = tx_gene.get(t)
        if g in gspan:
            per[g] += ex
    out = {}
    for g, ex in per.items():
        s, e = gspan[g]
        L = merged_len_cov_shorter_conjunct([(a - 1, b) for a, b in ex])
        out[f'{chrom}:{s}-{e}'] = L
        out[f'{chrom}:{s-1}-{e}'] = L
    return out

def main_cov_shorter_conjunct():
    ap = argparse.ArgumentParser()
    for x in ('--graph', '--paf', '--gff', '--truth'):
        ap.add_argument(x, required=True)
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--label', default='')
    ap.add_argument('--floors', default='0.30,0.50,0.70,0.90')
    a = ap.parse_args()

    adj = collections.defaultdict(dict); nodes = set()
    for ln in open(a.graph):
        f = ln.rstrip('\n').split('\t')
        if len(f) == 3 and f[0] != f[1]:
            w = float(f[2]); adj[f[0]][f[1]] = w; adj[f[1]][f[0]] = w
            nodes.add(f[0]); nodes.add(f[1])
        elif f and f[0]:
            nodes.add(f[0])
    keep = {(u, v) if u <= v else (v, u) for u in adj for v in adj[u]}
    # ⚠⚠ UNITS: the denominator must be the SHORTER locus's SPAN, not its exonic length. A genomic
    # alignment span over an exonic denominator is not a coverage fraction -- it ran to a median of 6.87
    # and a max of 2416 on chr16, which made a "0.70 floor" vacuous and produced a false no-op.
    # BLAST's scov is span-on-span; both sides genomic.
    elen = {}
    acc = collections.defaultdict(lambda: [[], []])
    for ln in open(a.paf):
        f = ln.rstrip('\n').split('\t')
        if len(f) < 11:
            continue
        q, t = f[0], f[5]
        if q == t:
            continue
        k = (q, t) if q <= t else (t, q)
        if k not in keep:
            continue
        elen[q] = int(f[1]); elen[t] = int(f[6])
        aq = (int(f[2]), int(f[3])); at = (int(f[7]), int(f[8]))
        e = acc[k]
        if k[0] == q:
            e[0].append(aq); e[1].append(at)
        else:
            e[0].append(at); e[1].append(aq)
    cov_short = {}
    for (x, y), (ix, iy) in acc.items():
        lx, ly = elen.get(x, 0), elen.get(y, 0)
        if not lx or not ly:
            continue
        cov_short[(x, y)] = min((merged_len_cov_shorter_conjunct(iy) / ly) if lx >= ly else (merged_len_cov_shorter_conjunct(ix) / lx), 1.0)

    # truth: gene symbol -> family, node -> gene by max overlap
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
    genes.sort(); starts = [g[0] for g in genes]

    def gname(n):
        m = re.match(r'(\S+):(\d+)-(\d+)$', n)
        if not m:
            return None
        s, e = int(m.group(2)), int(m.group(3)); best = None
        for i in range(max(0, bisect.bisect_left(starts, s) - 40), len(genes)):
            g0, g1, g = genes[i]
            if g0 > e:
                break
            ov = min(e, g1) - max(s, g0)
            if ov > 0 and (best is None or ov > best[0]):
                best = (ov, g)
        return best[1] if best else None

    lab = {}
    for r in csv.DictReader(open(a.truth), delimiter='\t'):
        f, g = (r.get('Family ID') or '').strip(), (r.get('Gene Name') or '').strip()
        if f and f != 'N/A' and g:
            lab.setdefault(g, f)
    ng = {n: gname(n) for n in nodes}
    nlab = {n: lab[g] for n, g in ng.items() if g in lab}
    byf = collections.defaultdict(set)
    for n, f in nlab.items():
        byf[f].add(n)
    truth = {f: v for f, v in byf.items() if len(v) >= 2}
    universe = set().union(*truth.values()) if truth else set()
    TP = set()
    for v in truth.values():
        vs = sorted(v)
        for i in range(len(vs)):
            for j in range(i + 1, len(vs)):
                TP.add((vs[i], vs[j]))

    def score(edges):
        fams = mcl_port.mcl({e: adj[e[0]][e[1]] for e in edges}, inflation=2.8)
        pred = set()
        for c in fams:
            m = sorted(set(c) & universe)
            for i in range(len(m)):
                for j in range(i + 1, len(m)):
                    pred.add((m[i], m[j]))
        tp = len(pred & TP)
        p = tp / len(pred) if pred else 0.0
        r = tp / len(TP) if TP else 0.0
        return p, r, (0.0 if p + r == 0 else 2 * p * r / (p + r))

    base = sorted(keep)
    p, r, f = score(base)
    print(f"{a.label or a.chrom}: edges {len(base)} | truth {len(truth)} fams / {len(TP)} pairs")
    print(f"  {'no floor (comparator)':<24} P {p:.3f} R {r:.3f} F {f:.3f}")
    for S in [float(x) for x in a.floors.split(',')]:
        kept = [e for e in base if cov_short.get(e, 1.0) >= S]
        p, r, f = score(kept)
        print(f"  {'S=' + format(S, '.2f'):<24} P {p:.3f} R {r:.3f} F {f:.3f}   "
              f"edges {len(kept)} (removed {len(base)-len(kept)} = {100*(len(base)-len(kept))/max(len(base),1):.1f}%)")


SUBCOMMANDS = {'poset': main_containment_poset, 'intron-chain': main_intron_chain_edges, 'protein-denovo': main_protein_sensitive_mode, 'protein-false-merge': main_protein_false_merge, 'cov-shorter': main_cov_shorter_conjunct}


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in SUBCOMMANDS:
        sys.exit('usage: edge_probes.py {' + '|'.join(SUBCOMMANDS) + '} [args]')
    SUBCOMMANDS[sys.argv.pop(1)]()


if __name__ == '__main__':
    main()
