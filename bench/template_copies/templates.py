#!/usr/bin/env python3
"""Templates of the template-guided copy stage: each assembled locus's representative transcript.

Spec docs/superpowers/specs/2026-10-14-template-copies-design.md, section 1, step 1. A locus is a `gene_id` of the base
GTF. Its representative is the transcript with the most `reads` (the assembler's attribute; 0 when absent), ties broken
by the longer exon sum, then the smaller transcript id. Coordinates are 0-based half-open (GTF lines are 1-based closed).
"""
import collections
import gzip
import re

GENE = re.compile(r'gene_id "([^"]+)"')
TX = re.compile(r'transcript_id "([^"]+)"')
READS = re.compile(r'\breads "([0-9]+)"')
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

Transcript = collections.namedtuple("Transcript", "tid gene contig strand exons reads")


def merge(iv):
    """Sorted union of half-open intervals (overlapping or touching ones join)."""
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_gtf(path):
    """GTF -> {transcript_id: Transcript}. Gene, contig and strand come from a transcript's first line; exons are merged."""
    info, reads, exons = {}, {}, collections.defaultdict(list)
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "transcript"):
                continue
            t = TX.search(f[8])
            if not t:
                continue
            tid = t.group(1)
            if tid not in info:
                g = GENE.search(f[8])
                info[tid] = (g.group(1) if g else tid, f[0], f[6])
            if f[2] == "exon":
                exons[tid].append((int(f[3]) - 1, int(f[4])))
            else:
                r = READS.search(f[8])
                reads[tid] = int(r.group(1)) if r else 0
    return {tid: Transcript(tid, g, c, s, merge(exons[tid]), reads.get(tid, 0))
            for tid, (g, c, s) in info.items() if exons.get(tid)}


def exon_len(t):
    return sum(e - s for s, e in t.exons)


def representatives(tx):
    """{gene: Transcript}: the most reads, then the longer exon sum, then the smaller transcript id."""
    best = {}
    for t in tx.values():
        b = best.get(t.gene)
        if b is None or (-t.reads, -exon_len(t), t.tid) < (-b.reads, -exon_len(b), b.tid):
            best[t.gene] = t
    return best


def loci(tx):
    """Figure 8 rows of the base loci: {'contig', 'locus', 'iv'}, iv = the union of the gene's transcripts' exons."""
    acc = {}
    for t in tx.values():
        acc.setdefault(t.gene, (t.contig, []))[1].extend(t.exons)
    return [{"contig": c, "locus": g, "iv": merge(ex)} for g, (c, ex) in sorted(acc.items())]


def spliced(fasta, t):
    """A transcript's spliced sequence in transcript orientation (`fasta` = pysam.FastaFile)."""
    seq = "".join(fasta.fetch(t.contig, s, e) for s, e in t.exons).upper()
    return seq.translate(COMP)[::-1] if t.strand == "-" else seq


def iv_str(iv):
    return ",".join(f"{s}-{e}" for s, e in iv)


def write_templates(reps, fasta, fa_path, tsv_path):
    """The templates, sorted by gene: FASTA (`>transcript_id`) and TSV (tid, gene, contig, strand, exons, reads)."""
    with open(fa_path, "w") as fa, open(tsv_path, "w") as tsv:
        tsv.write("tid\tgene\tcontig\tstrand\texons\treads\n")
        for g in sorted(reps):
            t = reps[g]
            fa.write(f">{t.tid}\n{spliced(fasta, t)}\n")
            tsv.write(f"{t.tid}\t{t.gene}\t{t.contig}\t{t.strand}\t{iv_str(t.exons)}\t{t.reads}\n")


def read_templates(tsv_path):
    """The templates TSV -> {tid: Transcript}."""
    out = {}
    with open(tsv_path) as fh:
        next(fh)
        for line in fh:
            tid, gene, contig, strand, ex, reads = line.rstrip("\n").split("\t")
            out[tid] = Transcript(tid, gene, contig, strand,
                                  [tuple(int(x) for x in b.split("-")) for b in ex.split(",")], int(reads))
    return out
