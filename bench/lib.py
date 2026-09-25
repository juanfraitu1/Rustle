#!/usr/bin/env python3
"""Shared helpers for the bench scripts: sequence, annotation, truth tables, pair scoring, reads.

Wave 7 (2026-09-24) folded the 21 top-level bench scripts into `lib.py` (this file), `score.py`, `sim.py`,
`truth.py`, the kept `guided_pipeline.py` and the untouched `mcl_port.py`. The old files are at git tag
`notebook-2026-09-24` (`git show notebook-2026-09-24:bench/<old>.py`).

Library names that moved here (old -> new):
  guided_pipeline.ov / merge / rc / pairwise / bipartite  -> lib.ov / merge / rc / pairwise / bipartite_items
                                                              (still re-exported by guided_pipeline as gp.*)
  adjudicated_truth.UF                                    -> lib.UF
  adjudicated_truth.translate (phase-aware, 0-based segs) -> lib.translate_phased
  protein_edge_gap.translate (1-based segs, no phase)     -> lib.translate_refseq
  protein_edge_gap.longest_cds                            -> lib.longest_cds
  heldout_family_score.score (per-truth-family bipartite) -> lib.bipartite_families
  soto_vs_us_referee.pair_scores                          -> lib.pair_scores
  soto_vs_us_referee.gene_names / protein_edge_gap gene_at -> lib.gene_key_names
  rna_truth_from_protein / soto_vs_us_referee biotype scan -> lib.gene_biotypes
  identity_spectrum.attr                                  -> lib.gtf_attr
  identity_spectrum Compara loader (twice in that file)   -> lib.load_compara
  locus_reads.aligned_blocks / reads_with_block_in /
    reads_overlapping_span / spanning_genes               -> lib.* (docstrings kept verbatim; THESIS_OBJECTIVES rules)
  copy_assign_tool_bakeoff introns() (twice)              -> lib.cigar_introns
  sim_reads.simulate_reads / write_fastq                  -> sim.simulate_reads / sim.write_fastq (not here)
  protein_families.excluded / pair_hsps / edges_from      -> truth.* (not here)

Two things with the same old name are DIFFERENT metrics and keep distinct names here:
  * `bipartite_items`   (item level: R = matched items / all items; P = matched / size of the matched predicted
                         clusters) — what protein_families / adjudicated_truth / guided_pipeline / layer_order report;
  * `bipartite_families` (family-macro: mean per-truth-family F, unmatched families score 0) — heldout's number.
Two "§6ko" translations likewise: `translate_refseq` (protein_edge_gap's, ignores CDS phase, truncates at the first
stop when the protein starts with M) and `translate_phased` (adjudicated_truth's, honours the first segment's phase).

Only the standard library is imported at module top; numpy/scipy are imported inside the functions that use them.
"""
import collections
import csv
import itertools
import re
import subprocess

# ---------------------------------------------------------------- sequence
# One complement table: every per-script copy (ACGT/acgt, with or without N/n) maps ACGTacgt identically and leaves
# every other character unchanged, so they are the same function.
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def ov(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def merge(iv):
    """Union of half-open intervals (touching intervals merge). The canonical copy from guided_pipeline."""
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


CODON = {}
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, (_x, _y, _z) in enumerate(itertools.product(_B, _B, _B)):
    CODON[_x + _y + _z] = _AA[_i]


def translate_refseq(fa, chrom, strand, segs):
    """protein_edge_gap.translate (the referee's §6ko proteome): `segs` are 1-based closed CDS segments; the CDS phase
    is IGNORED; a protein starting with M is cut at its first stop, any other has every stop replaced by X."""
    seq = ''.join(fa.fetch(chrom, s - 1, e) for s, e in segs).upper()
    if strand == '-':
        seq = rc(seq)
    aa = ''.join(CODON.get(seq[i:i + 3], 'X') for i in range(0, len(seq) - len(seq) % 3, 3))
    return aa.split('*')[0] if aa.startswith('M') else aa.replace('*', 'X')


def translate_phased(genome, chrom, strand, segs):
    """adjudicated_truth.translate (protein_families build): `segs` are 0-based (start, end, phase); the first
    segment's phase (last one on the minus strand) is honoured; trailing stops dropped, internal stops -> X."""
    segs = sorted(segs)
    seq = "".join(genome.fetch(chrom, s, e).upper() for s, e, _ in segs)
    phase = segs[0][2]
    if strand == "-":
        seq = rc(seq)
        phase = segs[-1][2]
    seq = seq[phase:]
    prot = "".join(CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))
    return prot.rstrip("*").replace("*", "X")


def spliced1(fa, chrom, exons, strand):
    """Spliced sequence from 1-based closed exons, reverse-complemented on the minus strand."""
    s = ''.join(fa.fetch(chrom, a - 1, b) for a, b in exons).upper()
    return rc(s) if strand == '-' else s


# ---------------------------------------------------------------- annotation
def gff_attrs(col):
    """GFF3 column 9 -> dict (annotation_nodes.attrs)."""
    return dict(kv.split("=", 1) for kv in col.strip().split(";") if "=" in kv)


def gtf_attr(s, k):
    """GTF attribute value `k "value"` (identity_spectrum.attr); None when absent."""
    m = re.search(k + r' "([^"]+)"', s)
    return m.group(1) if m else None


def gene_biotypes(gff, chrom):
    """gene/pseudogene `Name=` -> `gene_biotype=` ('' when absent) on `chrom`; last record wins
    (verbatim from rna_truth_from_protein / soto_vs_us_referee)."""
    bt = {}
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] not in ('gene', 'pseudogene'):
            continue
        n = re.search(r'Name=([^;]+)', f[8]); b = re.search(r'gene_biotype=([^;]+)', f[8])
        if n:
            bt[n.group(1)] = b.group(1) if b else ''
    return bt


def gene_key_names(gff, chrom):
    """'chrom:start1-end' -> gene/pseudogene `Name=` on `chrom` (the key mcl_families and the gene-body PAF use;
    verbatim from soto_vs_us_referee.gene_names and protein_edge_gap's gene_at)."""
    n = {}
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] not in ('gene', 'pseudogene'):
            continue
        m = re.search(r'Name=([^;]+)', f[8])
        if m:
            n[f'{f[0]}:{f[3]}-{f[4]}'] = m.group(1)
    return n


def longest_cds(gff, chrom):
    """gene symbol -> (strand, [(start1, end)]) for the transcript with the most CDS bases (protein_edge_gap)."""
    by_tx = collections.defaultdict(list)
    tx_gene, tx_strand = {}, {}
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom:
            continue
        if f[2] == 'CDS':
            p = re.search(r'Parent=([^;]+)', f[8]); g = re.search(r'gene=([^;]+)', f[8])
            if not p:
                continue
            tx = p.group(1)
            by_tx[tx].append((int(f[3]), int(f[4])))
            tx_strand[tx] = f[6]
            if g:
                tx_gene[tx] = g.group(1)
    best = {}
    for tx, segs in by_tx.items():
        g = tx_gene.get(tx)
        if not g:
            continue
        n = sum(e - s + 1 for s, e in segs)
        if g not in best or n > best[g][0]:
            best[g] = (n, tx_strand[tx], sorted(segs))
    return {g: (st, segs) for g, (n, st, segs) in best.items()}


# ---------------------------------------------------------------- truth tables
def soto_gene_family(s1c):
    """Soto et al. 2025 S1C (`Gene Name` / `Family ID`, header-based) -> {gene name: family id}.

    A gene carrying more than one distinct Family ID is EXCLUDED (the project's settled rule, D5: a partition needs
    one label per gene); `N/A` and empty IDs are dropped. Dict order = first appearance in the file."""
    ids = collections.defaultdict(set)
    for r in csv.DictReader(open(s1c), delimiter='\t'):
        fid = (r.get('Family ID') or '').strip(); nm = (r.get('Gene Name') or '').strip()
        if fid and fid != 'N/A' and nm:
            ids[nm].add(fid)
    return {nm: next(iter(f)) for nm, f in ids.items() if len(f) == 1}


def families_on(label_of, on_chrom, min_members):
    """{family: sorted members} over the genes of `label_of` present in `on_chrom`, families >= min_members."""
    fam = collections.defaultdict(set)
    for nm, fid in label_of.items():
        if nm in on_chrom:
            fam[fid].add(nm)
    return {k: sorted(v) for k, v in fam.items() if len(v) >= min_members}


def read_referee(path):
    """Two-column gene -> family table (`Gene Name`, `Family ID` header), read positionally."""
    fam = {}
    for ln in open(path):
        f = ln.rstrip('\n').split('\t')
        if f[0] == 'Gene Name' or len(f) < 2:
            continue
        fam[f[0]] = f[1]
    return fam


def load_compara(path, chrom):
    """Ensembl Compara BioMart paralogue table (gene, paralog, perc_id, perc_id_r1, subtype, paralog_chromosome).
    Returns ({frozenset(symbols): (max perc_id, subtype)} for same-chromosome pairs, {genes with any Compara row})."""
    chrom_num = chrom.replace('chr', '')
    compara, genes_with_data = {}, set()
    for ln in open(path):
        f = ln.rstrip('\n').split('\t')
        if len(f) < 6 or not f[0] or not f[1]:
            continue
        genes_with_data.add(f[0])
        if f[5] != chrom_num or f[0] == f[1]:
            continue
        try:
            pid = max(float(f[2] or 0), float(f[3] or 0))
        except ValueError:
            continue
        k = frozenset((f[0], f[1]))
        if k not in compara or pid > compara[k][0]:
            compara[k] = (pid, f[4])
    return compara, genes_with_data


# ---------------------------------------------------------------- partitions and pairs
class UF:
    """Union-find; the root of a component is its minimum element (order-independent)."""

    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[max(ra, rb)] = min(ra, rb)


def pairwise(pred, true):
    """Pairwise (sensitivity, precision) of label vector `pred` against `true` (NaN when undefined).

    Same numbers as the old O(n^2) guided_pipeline.pairwise, from contingency counts in O(n):
    tp = sum C(n_pt, 2), predicted pairs = sum C(b_p, 2), true pairs = sum C(a_t, 2) (integer arithmetic)."""
    c2 = lambda n: n * (n - 1) // 2
    pairs = list(zip(pred, true))
    tp = sum(c2(n) for n in collections.Counter(pairs).values())
    npred = sum(c2(n) for n in collections.Counter(p for p, _ in pairs).values())
    ntrue = sum(c2(n) for n in collections.Counter(t for _, t in pairs).values())
    fp, fn = npred - tp, ntrue - tp
    return (tp / (tp + fn) if tp + fn else float("nan")), (tp / (tp + fp) if tp + fp else float("nan"))


def bipartite_items(pred, true):
    """ITEM-level one-to-one matching of predicted to true labels (scipy linear_sum_assignment on the contingency
    matrix, labels sorted by str): returns (matched items / all items, matched items / size of the matched predicted
    clusters). The old guided_pipeline.bipartite, with dict indexing instead of list.index (same matrix)."""
    import numpy as np
    from scipy.optimize import linear_sum_assignment
    P, T = sorted(set(pred), key=str), sorted(set(true), key=str)
    pi = {p: j for j, p in enumerate(P)}
    ti = {t: i for i, t in enumerate(T)}
    M = np.zeros((len(T), len(P)), dtype=int)
    for p, t in zip(pred, true):
        M[ti[t], pi[p]] += 1
    r, c = linear_sum_assignment(-M)
    matched = sum(M[i, j] for i, j in zip(r, c))
    sp = M.sum(axis=0)
    msize = sum(sp[j] for i, j in zip(r, c) if M[i, j] > 0)
    return matched / len(pred), (matched / msize if msize else float("nan"))


def bipartite_families(truth, pred):
    """FAMILY-macro bipartite (heldout_family_score.score): one-to-one match of truth families to predicted clusters
    maximising total member overlap; per truth family sens = matched / truth members, prec = matched / members of the
    matched cluster, F = harmonic mean; an unmatched truth family scores 0 and is KEPT. Returns {root: dict} or None."""
    import numpy as np
    from scipy.optimize import linear_sum_assignment
    troots, cids = sorted(truth), sorted(pred)
    if not troots:
        return None
    ovm = np.zeros((len(troots), len(cids)), dtype=int)
    for i, r in enumerate(troots):
        tset = set(truth[r])
        for j, c in enumerate(cids):
            ovm[i, j] = len(tset & set(pred[c]))
    rows, cols = linear_sum_assignment(-ovm) if cids else ([], [])
    matched = {troots[i]: cids[j] for i, j in zip(rows, cols) if ovm[i, j] > 0}
    per = {}
    for r in troots:
        t = set(truth[r])
        c = matched.get(r)
        if c is None:
            per[r] = dict(n_truth=len(t), n_pred=0, hit=0, sens=0.0, prec=0.0, f=0.0, cluster=None)
            continue
        p = set(pred[c])
        hit = len(t & p)
        sens = hit / len(t)
        prec = hit / len(p) if p else 0.0
        f = 0.0 if sens + prec == 0 else 2 * sens * prec / (sens + prec)
        per[r] = dict(n_truth=len(t), n_pred=len(p), hit=hit, sens=round(sens, 4),
                      prec=round(prec, 4), f=round(f, 4), cluster=c)
    return per


def pair_scores(label_of, ref):
    """precision/recall/F against the referee over a FIXED UNIVERSE: every referee-labelled gene is
    scored, and a gene the comparator never placed becomes its own singleton rather than being dropped.

    ⚠ Restricting the universe to genes the comparator labelled conditions the denominator on the
    prediction (register 770) and returns precision 1.000 by construction. Do not do that.
    """
    genes = list(ref)
    label_of = {g: label_of.get(g, f'__singleton__{g}') for g in genes}
    byc = collections.defaultdict(list)
    for g in genes:
        byc[label_of[g]].append(g)
    byr = collections.defaultdict(list)
    for g in genes:
        byr[ref[g]].append(g)
    pred = set()
    for v in byc.values():
        for i in range(len(v)):
            for j in range(i + 1, len(v)):
                pred.add(tuple(sorted((v[i], v[j]))))
    true = set()
    for v in byr.values():
        for i in range(len(v)):
            for j in range(i + 1, len(v)):
                true.add(tuple(sorted((v[i], v[j]))))
    tp = len(pred & true)
    p = tp / len(pred) if pred else 0.0
    r = tp / len(true) if true else 0.0
    f = 0.0 if p + r == 0 else 2 * p * r / (p + r)
    return p, r, f, len(pred), len(true), tp


# ---------------------------------------------------------------- alignments (PAF)
def paf_identity(f, estimator):
    """Identity of one PAF record (split fields). No default on purpose (D7): 'nm_bl' = matches / block length;
    'de' = 1 - de:f (the Rust E_r builder), falling back to nm_bl when the record has no de tag."""
    if estimator == 'nm_bl':
        return int(f[9]) / int(f[10])
    if estimator == 'de':
        for x in f[12:]:
            if x.startswith('de:f:'):
                return 1.0 - float(x[5:])
        return int(f[9]) / int(f[10])
    raise ValueError(f'unknown identity estimator {estimator!r}')


def paf_coverage(f, rule):
    """Coverage of one PAF record. No default on purpose (D8): 'query_over_min' = query span / min(qlen, tlen) (the M1
    defect: exceeds 1 when the query is the longer sequence); 'shorter_axis' = the aligned span of whichever sequence
    is shorter, over its length (the Rust E_r builder)."""
    ql, tl = int(f[1]), int(f[6])
    if rule == 'query_over_min':
        return (int(f[3]) - int(f[2])) / min(ql, tl)
    if rule == 'shorter_axis':
        return (int(f[3]) - int(f[2])) / ql if ql <= tl else (int(f[8]) - int(f[7])) / tl
    raise ValueError(f'unknown coverage rule {rule!r}')


# ---------------------------------------------------------------- reads (was bench/locus_reads.py)
# Counting reads at a locus — the ONE correct way, and the wrong way named so it cannot be reached by accident.
#
# ⛔⛔ THE ERROR THIS EXISTS TO PREVENT (ledger §6cm, 2026-09-02). Counting reads that *overlap a locus span* instead
# of reads that have an *aligned block inside it*. On 2026-09-02 that inflated a headline 3.4x and produced a whole
# retracted mechanism: `NPIPP1` appeared to have 1,608 reads collapsing into a 4-read locus, when 1,151 of them —
# 71.6% — merely spliced across it with no aligned base, because `PDXDC1` (168 kb) physically CONTAINS `NPIPP1` and
# its transcript passes straight through. The real figure is 457, itself an upper bound.
#
# The project rule is old and was not enforced anywhere: `N` in an RNA CIGAR is an intron, spliced OUT; a read spliced
# OVER a locus is no evidence for it. `samtools view -c REGION` counts the wrong thing, silently, and reads
# beautifully in a script.
#
# Use `reads_with_block_in()`. `reads_overlapping_span()` exists only so the difference is visible and so a caller who
# genuinely wants span overlap has to say so.
#
# ⚠ Primary alignments only (`-F 2308`) throughout — the standing invariant before any per-read CIGAR statistic, so
# one molecule is never two witnesses.
SAMTOOLS = "/home/juanfra/miniforge3/bin/samtools"
_CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def aligned_blocks(pos0, cigar):
    """Reference blocks a read actually aligns to. `N` breaks the run; `M/=/X/D` extend it."""
    out, p, cur = [], pos0, None
    for ln, op in _CIG.findall(cigar):
        ln = int(ln)
        if op in "M=XD":
            cur = (cur[0] if cur else p, p + ln)
            p += ln
        elif op == "N":
            if cur:
                out.append(cur)
                cur = None
            p += ln
    if cur:
        out.append(cur)
    return out


def _fetch(bam, chrom, start, end):
    r = subprocess.run(
        [SAMTOOLS, "view", "-F", "2308", bam, f"{chrom}:{start + 1}-{end}"],
        capture_output=True, text=True,
    )
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) > 5 and f[5] != "*":
            yield f


def reads_with_block_in(bam, chrom, start, end, blocks=None):
    """⭐ THE CORRECT COUNT: primaries with >= 1 ALIGNED BASE inside [start, end).

    `blocks` optionally restricts to a locus's own exon blocks — stricter still, and what you want
    when a large gene's exons fall inside the span of a small one nested within it.
    """
    n = 0
    for f in _fetch(bam, chrom, start, end):
        bl = aligned_blocks(int(f[3]) - 1, f[5])
        if blocks is None:
            if any(bs < end and start < be for bs, be in bl):
                n += 1
        elif any(bs < e2 and s2 < be for bs, be in bl for s2, e2 in blocks):
            n += 1
    return n


def reads_overlapping_span(bam, chrom, start, end):
    """⚠ THE MISLEADING COUNT — what `samtools view -c` gives. Includes reads that splice straight
    over the locus contributing no aligned base. Returns `(n_overlapping, n_spliced_over)` so the
    caller cannot quote the first without seeing the second."""
    tot = over = 0
    for f in _fetch(bam, chrom, start, end):
        tot += 1
        bl = aligned_blocks(int(f[3]) - 1, f[5])
        if not any(bs < end and start < be for bs, be in bl):
            over += 1
    return tot, over


def spanning_genes(gff_gz, chrom, start, end, min_cover=0.60):
    """⭐ ASK THIS BEFORE BLAMING THE PIPELINE for an oversized locus (§6cm).

    Returns annotated genes covering >= `min_cover` of the node. On 2026-09-02 two nodes called
    "mis-chained giants" turned out to be `SNX29` (covers 100.0%) and `PDXDC1` (99.9%) — real genes,
    correctly assembled, with canonical junctions and majority read support because they are real
    transcripts. A rule fitted to that "blob class" would have been fitted to two real genes.
    """
    import gzip
    out = []
    with gzip.open(gff_gz, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[0] != chrom or f[2] not in ("gene", "pseudogene"):
                continue
            gs, ge = int(f[3]) - 1, int(f[4])
            if gs < end and start < ge:
                cov = (min(ge, end) - max(gs, start)) / max(1, end - start)
                if cov >= min_cover:
                    m = re.search(r"Name=([^;\n]+)", f[8])
                    out.append((m.group(1) if m else "?", gs, ge, cov))
    return sorted(out, key=lambda x: -x[3])


def cigar_introns(pos, cig):
    """Intron chain ((start, end), ...) of a SAM record, 0-based `pos` (copy_assign_tool_bakeoff's introns())."""
    o, p = [], pos
    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cig):
        n = int(n)
        if op in 'M=XD':
            p += n
        elif op == 'N':
            o.append((p, p + n))
            p += n
    return tuple(o)


def sam_lines(args):
    """Stream `samtools view ARGS` line by line (no whole-BAM string in memory, B6). Lines keep no trailing newline."""
    p = subprocess.Popen(['samtools', 'view'] + list(args), stdout=subprocess.PIPE, text=True)
    try:
        for ln in p.stdout:
            yield ln.rstrip('\n')
    finally:
        p.stdout.close()
        p.wait()


def mcl(edges, inflation=2.8, prune=1e-9, max_iter=100):
    """Re-export of `mcl_port.mcl` (imported lazily; mcl_port.py itself is unchanged)."""
    import mcl_port
    return mcl_port.mcl(edges, inflation=inflation, prune=prune, max_iter=max_iter)
