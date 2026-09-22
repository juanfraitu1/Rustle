"""Multi-label locus-to-gene resolution for scoring, per
`docs/PREREG_multilabel_locus_resolver_2026-09-21.md` (md5 `9d68f84c`).

Replaces winner-take-all (`node_to_gene_max_overlap` in earlier session diagnostics): a locus credits
EVERY gene whose OWN exonic content it covers by >= FLOOR, not just the single largest-overlap gene.
This is what let readthrough super-loci hide real, correctly-assembled genes behind a bigger neighbour
(register 962/963) -- NPIPB4 and RRN3P1 sharing one locus and both losing to a third gene, LOC112268174.

    credit(locus, g)  iff  |E_locus (intersect) E_g| / |E_g| >= FLOOR      (FLOOR = 0.50, exonic, not span)

Does not touch the truth side; only changes how a PREDICTED locus's identity is read (register 770 safe).
"""
import collections
import re
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_bridged_merge import merge_iv, ov_bp  # noqa: E402

FLOOR = 0.50

TX = {'mRNA', 'transcript', 'ncRNA', 'lnc_RNA', 'lncRNA', 'pseudogenic_transcript',
      'primary_transcript', 'tRNA', 'rRNA', 'snRNA', 'snoRNA', 'miRNA', 'misc_RNA',
      'V_gene_segment', 'C_gene_segment', 'J_gene_segment', 'ncRNA_gene'}


def gene_exon_unions(gff, chrom):
    """gene name -> merged exon list, unioned across every annotated transcript of that gene.
    Handles a pseudogene exon parented straight to the gene record (no transcript in between) --
    the same gap that silently dropped 534 of chr16's pseudogene exons earlier this session.

    TWO PASSES, deliberately -- a single pass assumes parent-before-child file order, which RefSeq's
    own GFF does NOT guarantee: `SLC7A5P2`'s single exon (`exon-NR_002594.1-1`) sits at line 71119,
    one line ABOVE its own gene (71120) and transcript (71121) records, because gene/transcript/exon
    share identical coordinates here and the file's sort has no tiebreak on feature type. A single-pass
    parser silently drops that exon (`t2g`/`gene_of` aren't populated yet when the exon line is read) --
    found because it produced a gene with a 0 bp exonic union for a gene known to carry a transcript.
    Pass 1 builds every gene_of/t2g mapping; pass 2 reads exons with the mappings already complete."""
    lines = [ln for ln in open(gff) if not ln.startswith('#')]
    gene_of, t2g = {}, {}
    for line in lines:
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom:
            continue
        if f[2] in ('gene', 'pseudogene', 'ncRNA_gene'):
            n = re.search(r'Name=([^;]+)', f[8]); i = re.search(r'ID=([^;]+)', f[8])
            if n and i:
                gene_of[i.group(1)] = n.group(1)
    for line in lines:
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] not in TX:
            continue
        i = re.search(r'ID=([^;]+)', f[8]); p = re.search(r'Parent=([^;,]+)', f[8])
        if i and p and p.group(1) in gene_of:
            t2g[i.group(1)] = gene_of[p.group(1)]
    ex = collections.defaultdict(list)
    for line in lines:
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] != 'exon':
            continue
        p = re.search(r'Parent=([^;,]+)', f[8])
        if not p:
            continue
        par = p.group(1)
        if par in t2g:
            ex[t2g[par]].append((int(f[3]), int(f[4])))
        elif par in gene_of:
            ex[gene_of[par]].append((int(f[3]), int(f[4])))
    return {g: merge_iv(v) for g, v in ex.items()}


def locus_multi_labels(locus_exons, gene_ex, floor=FLOOR):
    """{gene_name: coverage_fraction} for every gene the locus covers by >= floor of the GENE's own
    exonic content. Both the label set and the fractions are returned -- callers that want the old
    winner-take-all behaviour can still take argmax, but the point of this function is not to."""
    if not locus_exons:
        return {}
    lo, hi = locus_exons[0][0], locus_exons[-1][1]
    out = {}
    for g, gex in gene_ex.items():
        if not gex or gex[-1][1] < lo or gex[0][0] > hi:
            continue
        glen = sum(e - s + 1 for s, e in gex)
        if not glen:
            continue
        ov = ov_bp(locus_exons, gex)
        if ov <= 0:
            continue
        frac = ov / glen
        if frac >= floor:
            out[g] = frac
    return out


def build_node_labels(exm, gene_ex, floor=FLOOR):
    """{locus_id: {gene_name: frac}} for every locus in `exm` (locus_id -> merged exon list)."""
    return {lid: locus_multi_labels(ex, gene_ex, floor) for lid, ex in exm.items()}
