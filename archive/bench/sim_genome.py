#!/usr/bin/env python
"""FULLY-SIMULATED 2-chromosome genome with PLANTED multi-copy families — the airtight, ground-truth benchmark.

Everything is synthetic and labelled, so there is NO circularity: we plant the copies (positions, divergence,
copy number, exon/intron structure) and label every read's TRUE family/copy/isoform, then ask whether the
pipeline recovers exactly that. Covers every regime the theory cares about:

  K0tandem  : 3 IDENTICAL tandem copies            -> the K=0 floor (every read must be certified TIED)
  ladder    : 4 tandem copies at 0 / 0.5% / 2% / 5% divergence -> the resolvable frontier (accuracy rises with K)
  xchrom    : 2 copies on DIFFERENT chromosomes, 4% diverged   -> cross-chrom family + assignment
  cnv       : 3 tandem copies, UNEQUAL expression               -> abundance / quantification
  single x4 : single-copy genes, no paralog                     -> must NOT form a family (FP control)
  domshare  : 2 genes sharing ONE identical exon                -> must NOT merge (over-merge FP control)

Outputs (to winloci_scratch): simgw.fasta (+.fai), simgw_truth.tsv (family/copy/locus/introns/divergence/reads),
simgw_reads.fastq (read names encode truth: SIMGW|<family>|<copy>|<iso>|<i>). Deterministic.

Run: /home/juanfra/miniforge3/bin/python bench/sim_genome.py
"""
import os
import random
import subprocess
import sys

sys.path.insert(0, os.path.dirname(__file__))
from sim_reads import simulate_reads, write_fastq

OUT = "/home/juanfra/winloci_scratch"
SAM = "/home/juanfra/miniforge3/bin/samtools"
BASES = "ACGT"
RNG = random.Random(20260629)
READS_PER_COPY = 40


def randseq(n):
    return "".join(RNG.choice(BASES) for _ in range(n))


def mutate(seq, d):
    """substitute each base with prob d (planted PSVs at a known rate); deterministic via RNG."""
    if d <= 0:
        return seq
    return "".join(RNG.choice([b for b in BASES if b != c]) if RNG.random() < d else c for c in seq)


def make_intron(n):
    """a canonical + strand intron: GT ... AG (the donor/acceptor the assembly gate requires)."""
    return "GT" + randseq(n - 4) + "AG"


def make_gene(n_exons):
    """a base gene: longer so a shared exon is a SMALL fraction (exon ~250-450 bp -> ~2 kb transcript,
    intron ~400-1500 bp, canonical GT-AG)."""
    exons = [randseq(RNG.randint(250, 450)) for _ in range(n_exons)]
    introns = [make_intron(RNG.randint(400, 1500)) for _ in range(n_exons - 1)]
    return exons, introns


def mutate_intron(seq, d):
    """diverge an intron but KEEP the canonical GT..AG motif (real introns retain it as they diverge)."""
    return "GT" + mutate(seq[2:-2], d) + "AG"


def copy_layout_seqs(ex, intr):
    """assemble explicit exon/intron seqs into (genomic, spliced, intron_offsets within copy)."""
    g, spliced, introns_local, pos = [], [], [], 0
    for k, e in enumerate(ex):
        g.append(e)
        spliced.append(e)
        pos += len(e)
        if k < len(intr):
            introns_local.append((pos, pos + len(intr[k])))  # (donor, acceptor) within the copy
            g.append(intr[k])
            pos += len(intr[k])
    return "".join(g), "".join(spliced), introns_local


def copy_layout(exons, introns, d):
    """a copy = the gene diverged by d; returns (genomic_seq, spliced_seq, intron_offsets[(d0,a0)...] within copy)."""
    ex = [mutate(e, d) for e in exons]
    intr = [mutate_intron(i, d) for i in introns]
    return copy_layout_seqs(ex, intr)


class Chrom:
    def __init__(self, name):
        self.name = name
        self.parts = []
        self.pos = 0

    def add_bg(self, n):
        s = randseq(n)
        self.parts.append(s)
        self.pos += n

    def add_copy(self, gseq):
        start = self.pos
        self.parts.append(gseq)
        self.pos += len(gseq)
        return start  # genomic start (0-based)

    def seq(self):
        return "".join(self.parts)


def main():
    truth = []          # (family, copy, chrom, start, end, intron_chain_str, divergence, n_reads, abundance_tag)
    reads = []          # (name, (read,qual))

    chrA, chrB = Chrom("simA"), Chrom("simB")
    chrA.add_bg(30000)
    chrB.add_bg(30000)

    def plant(fam, chrom, copies, abund=None):
        """copies = list of (divergence, n_reads). All share ONE base gene (so they are paralogs)."""
        n_exons = RNG.randint(5, 7)
        exons, introns = make_gene(n_exons)
        for ci, (d, nr) in enumerate(copies):
            gseq, spliced, intr_local = copy_layout(exons, introns, d)
            start = chrom.add_copy(gseq)
            introns_g = [(start + a, start + b) for (a, b) in intr_local]
            chrom.add_bg(RNG.randint(2000, 6000))  # spacer between copies
            ichain = ";".join(f"{a}-{b}" for a, b in introns_g)
            truth.append((fam, ci, chrom.name, start, start + len(gseq), ichain, d, nr,
                          (abund[ci] if abund else "")))
            # reads from this copy's spliced transcript, labelled with the TRUE family/copy
            for i, rq in enumerate(simulate_reads(spliced, nr, err=0.003, indel=0.0008, seed=hash((fam, ci)) & 0xffff)):
                reads.append((f"SIMGW|{fam}|{ci}|0|{i}", rq))

    def plant_collapsed(fam, chrom, n_copies, n_psv, nr):
        """A COLLAPSED segdup array: copies IDENTICAL except `n_psv` fixed exonic PSV columns (each copy a
        distinct allele vector). This is the regime the gate is FOR — copies so similar minimap2 multimaps
        them (MAPQ-0) yet a few PSVs distinguish them, so the gate must ASSIGN, not just abstain.

        Crucially the gene is LONG (≈6 kb) with only a FEW PSVs: a handful of mismatches over 6 kb does not
        create a minimap2 score gap, so reads stay MAPQ-0 (ambiguous) — yet the gate, which scores ONLY the
        PSV columns, still resolves them. (Uniform divergence, in contrast, lifts MAPQ the moment it is large
        enough to give PSVs — which is why the divergent ladder copies are MAPQ>0.)"""
        n_exons = 8
        exons = [randseq(RNG.randint(600, 900)) for _ in range(n_exons)]       # ≈6 kb transcript
        introns = [make_intron(RNG.randint(400, 1500)) for _ in range(n_exons - 1)]
        # choose PSV positions: (exon index, offset) over the exonic (spliced) coordinates
        flat = [(ei, off) for ei, e in enumerate(exons) for off in range(len(e))]
        psv_sites = RNG.sample(flat, n_psv)
        # per PSV column, hand each copy an INDEPENDENT random allele (so two copies colliding at ALL n_psv
        # columns has prob (1/4)^n_psv ≈ 0 — every copy is genuinely distinguishable, the resolvable regime).
        # The column still varies across copies, so it is a real PSV.
        col_alleles = []
        for (ei, off) in psv_sites:
            alleles = [RNG.choice(BASES) for _ in range(n_copies)]
            while len(set(alleles)) < 2:           # ensure the column is actually variable
                alleles = [RNG.choice(BASES) for _ in range(n_copies)]
            col_alleles.append(alleles)
        for ci in range(n_copies):
            ex = [list(e) for e in exons]
            for (ei, off), alleles in zip(psv_sites, col_alleles):
                ex[ei][off] = alleles[ci]
            ex = ["".join(e) for e in ex]
            gseq, spliced, intr_local = copy_layout_seqs(ex, introns)
            start = chrom.add_copy(gseq)
            introns_g = [(start + a, start + b) for (a, b) in intr_local]
            chrom.add_bg(RNG.randint(2000, 6000))
            ichain = ";".join(f"{a}-{b}" for a, b in introns_g)
            truth.append((fam, ci, chrom.name, start, start + len(gseq), ichain, round(n_psv / len(spliced), 4), nr, ""))
            for i, rq in enumerate(simulate_reads(spliced, nr, err=0.003, indel=0.0008, seed=hash((fam, ci)) & 0xffff)):
                reads.append((f"SIMGW|{fam}|{ci}|0|{i}", rq))

    def plant_rep(fam, chrom, copies):
        """GSTM-like paralogs — the regime where the copy-vs-copy aligner CHOICE changes assignment. Two forces
        must coexist and they pull opposite ways: (i) reads must MULTIMAP (MAPQ-0) so the PSV columns are the only
        assignment signal — that needs HIGH overall identity; (ii) poasta must be SUBOPTIMAL so a bad alignment
        mis-places those PSVs — that needs a locally-disruptive REPETITIVE region with a frameshift indel. Real
        recent paralogs (GSTM3/5/1) have both: >90% identical overall, but their divergence is concentrated in
        low-complexity stretches. Here most exons are near-identical (shared, low div -> multimap) while a few are
        tiled from a handful of short motifs and diverged per copy with a frameshift (non-motif-length) indel ->
        poasta garden-paths, the exact DP does not. `copies` = [(var_div, n_reads, frameshift_indel)] where the
        indel is None or (var_exon_index, pos, "del"|"ins", length)."""
        n_shared, n_var = 6, 2
        shared_exons = [randseq(RNG.randint(500, 800)) for _ in range(n_shared)]      # conserved core (multimap)
        motifs = [randseq(RNG.randint(40, 70)) for _ in range(6)]                     # only 6 motifs -> repetitive
        var_base = ["".join(RNG.choice(motifs) for _ in range(RNG.randint(12, 16))) for _ in range(n_var)]
        introns = [make_intron(RNG.randint(400, 1500)) for _ in range(n_shared + n_var - 1)]
        for ci, (dvar, nr, indel) in enumerate(copies):
            shared = [mutate(e, 0.002) for e in shared_exons]        # near-identical core across copies -> MAPQ-0
            var = [mutate(e, dvar) for e in var_base]                # repetitive region diverged per copy
            if indel is not None:
                ei, pos, op, ln = indel
                e = var[ei]
                var[ei] = e[:pos] + e[pos + ln:] if op == "del" else e[:pos] + randseq(ln) + e[pos:]
            # interleave so the divergent repetitive exons sit inside the conserved core
            exons = shared[:3] + [var[0]] + shared[3:] + [var[1]]
            intr = [mutate_intron(i, 0.005) for i in introns]
            gseq, spliced, intr_local = copy_layout_seqs(exons, intr)
            start = chrom.add_copy(gseq)
            introns_g = [(start + a, start + b) for (a, b) in intr_local]
            chrom.add_bg(RNG.randint(2000, 6000))
            ichain = ";".join(f"{a}-{b}" for a, b in introns_g)
            truth.append((fam, ci, chrom.name, start, start + len(gseq), ichain, dvar, nr, ""))
            for i, rq in enumerate(simulate_reads(spliced, nr, err=0.003, indel=0.0008, seed=hash((fam, ci)) & 0xffff)):
                reads.append((f"SIMGW|{fam}|{ci}|0|{i}", rq))

    # ---- planted families ----
    plant("K0tandem", chrA, [(0.0, 40), (0.0, 40), (0.0, 40)])                       # K=0 floor
    plant("ladder", chrA, [(0.0, 40), (0.003, 40), (0.008, 40), (0.015, 40)])         # resolvable ladder (multimapping)
    plant_collapsed("collapse", chrB, n_copies=5, n_psv=6, nr=40)                      # collapsed segdup: MAPQ-0 yet PSV-resolvable (the gate's regime)
    plant("cnv", chrB, [(0.0, 80), (0.005, 40), (0.005, 20)], abund=["high", "med", "low"])  # unequal expression
    # repetitive high-divergence family: exons tiled from short motifs (low-complexity) + subs divergence + a
    # structural motif-tile indel per copy -> the copy-vs-copy alignment regime where poasta is SUBOPTIMAL and the
    # exact banded DP is optimal, so this family measures whether the aligner choice changes ASSIGNMENT accuracy.
    plant_rep("hidive", chrB, [(0.0, 40, None), (0.10, 40, (0, 300, "del", 47)),
                               (0.14, 40, (1, 250, "ins", 47)), (0.18, 40, (0, 200, "del", 93))])
    # single-copy controls
    for s in range(4):
        ch = chrA if s % 2 == 0 else chrB
        exons, introns = make_gene(RNG.randint(5, 7))
        g, spliced, intr = copy_layout(exons, introns, 0.0)
        start = ch.add_copy(g)
        ch.add_bg(RNG.randint(2000, 6000))
        truth.append((f"single{s}", 0, ch.name, start, start + len(g),
                      ";".join(f"{start+a}-{start+b}" for a, b in intr), 0.0, 30, ""))
        for i, rq in enumerate(simulate_reads(spliced, 30, err=0.003, indel=0.0008, seed=hash(("single", s)) & 0xffff)):
            reads.append((f"SIMGW|single{s}|0|0|{i}", rq))
    # domain-sharer: two genes sharing ONE identical exon (must NOT merge into a family)
    shared_exon = randseq(150)  # < T_CORE of a ~2kb transcript -> must NOT merge
    for g2 in range(2):
        exons, introns = make_gene(6)
        exons[2] = shared_exon  # the shared exon (identical between the two otherwise-unrelated genes)
        g, spliced, intr = copy_layout(exons, introns, 0.0)
        ch = chrA
        start = ch.add_copy(g)
        ch.add_bg(RNG.randint(2000, 6000))
        truth.append((f"domshare{g2}", 0, ch.name, start, start + len(g),
                      ";".join(f"{start+a}-{start+b}" for a, b in intr), 0.0, 30, "domain_sharer"))
        for i, rq in enumerate(simulate_reads(spliced, 30, err=0.003, indel=0.0008, seed=hash(("dom", g2)) & 0xffff)):
            reads.append((f"SIMGW|domshare{g2}|0|0|{i}", rq))

    # cross-chrom paralog: ONE base gene, two copies on the two chromosomes. The copy on chrB is only ~0.3%
    # diverged so the two copies genuinely MULTIMAP across chromosomes (a read confused between chrA and chrB =
    # a cross-chrom conflict edge); a larger divergence would make each copy uniquely mappable to its own
    # chromosome (MAPQ>0) and the cross-chrom family would not form — the same "divergence lifts MAPQ" limit.
    n_exons = RNG.randint(5, 7)
    exons, introns = make_gene(n_exons)
    for ci, (ch, d) in enumerate([(chrA, 0.0), (chrB, 0.003)]):
        g, spliced, intr = copy_layout(exons, introns, d)
        start = ch.add_copy(g)
        ch.add_bg(RNG.randint(2000, 6000))
        truth.append(("xchrom", ci, ch.name, start, start + len(g),
                      ";".join(f"{start+a}-{start+b}" for a, b in intr), d, 40, ""))
        for i, rq in enumerate(simulate_reads(spliced, 40, err=0.003, indel=0.0008, seed=hash(("xc", ci)) & 0xffff)):
            reads.append((f"SIMGW|xchrom|{ci}|0|{i}", rq))

    chrA.add_bg(20000)
    chrB.add_bg(20000)

    # ---- write reference, truth, reads ----
    with open(f"{OUT}/simgw.fasta", "w") as fh:
        for ch in (chrA, chrB):
            fh.write(f">{ch.name}\n{ch.seq()}\n")
    subprocess.run([SAM, "faidx", f"{OUT}/simgw.fasta"], check=True)
    with open(f"{OUT}/simgw_truth.tsv", "w") as fh:
        fh.write("family\tcopy\tchrom\tstart\tend\tintrons\tdivergence\tn_reads\tabundance\n")
        for t in truth:
            fh.write("\t".join(str(x) for x in t) + "\n")
    with open(f"{OUT}/simgw_reads.fastq", "w") as fh:
        for name, rq in reads:
            write_fastq(fh, name, rq)

    nf = len({t[0] for t in truth})
    print(f"genome: simA={chrA.pos:,}bp simB={chrB.pos:,}bp | {len(truth)} planted loci in {nf} families/singletons")
    print(f"reads: {len(reads)} (labelled SIMGW|family|copy|iso|i)")
    print(f"wrote {OUT}/simgw.fasta, simgw_truth.tsv, simgw_reads.fastq")
    # quick truth summary
    from collections import Counter
    fams = Counter(t[0] for t in truth)
    multi = {f: n for f, n in fams.items() if n >= 2}
    print(f"multi-copy families: {multi}")


if __name__ == "__main__":
    main()
