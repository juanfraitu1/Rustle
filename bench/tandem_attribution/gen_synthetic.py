#!/usr/bin/env python3
"""Generate a synthetic tandem paralog array with KNOWN per-read source copy, for
the per-copy attribution benchmark.

Design (and why each choice):
  * N near-identical copies of a multi-exon gene, laid out as a TANDEM ARRAY
    (copies `spacing` bp apart) so the array collapses into one bundle and
    `RUSTLE_VG_TANDEM=1` fires.
  * Copy-distinguishing mutations are placed ONLY IN EXONS, at a rate set to hit
    a target EXONIC pairwise identity. This is the crux: intronic differences
    never reach a spliced read (the DAZ failure mode), so only exonic divergence
    is informative for read->copy assignment. Sweeping the exonic identity traces
    the identifiability spectrum.
  * Introns are identical across copies (all distinguishing signal is exonic).
  * Reads are full-length spliced mRNA (some 5'-truncated for realism) from each
    copy, with IsoSeq-like random substitution errors (non-diagnostic noise the
    EM must see through). Each read is named `c{src}_r{j}` so its true source copy
    is recoverable.

Outputs (in --out):
  genome.fa   the multi-copy genomic array (the reference reads are aligned to)
  reads.fa    spliced reads, named c{src}_r{j}
  meta.json   parameters + realized exonic identity

Usage:
  gen_synthetic.py --identity 0.99 --copies 3 --reads-per-copy 40 --out DIR --seed 1
"""
import argparse, json, os, random

BASES = "ACGT"

def rand_seq(n, rng):
    return [rng.choice(BASES) for _ in range(n)]

def mutate_positions(seq, positions, rng):
    """Substitute each given position to a DIFFERENT base (guarantees a real diff)."""
    for p in positions:
        cur = seq[p]
        seq[p] = rng.choice([b for b in BASES if b != cur])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--identity", type=float, required=True,
                    help="target EXONIC pairwise identity between copies, e.g. 0.99")
    ap.add_argument("--copies", type=int, default=3)
    ap.add_argument("--reads-per-copy", type=int, default=40)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--error-rate", type=float, default=0.005, help="IsoSeq CCS substitution rate")
    ap.add_argument("--trunc-frac", type=float, default=0.25, help="fraction of reads 5'-truncated")
    ap.add_argument("--spacing", type=int, default=23000,
                    help="copy-to-copy start distance. Large (e.g. 200000) = dispersed loci "
                         "(separate bundles, normal-discovery family); ~9500 = tandem-collapse.")
    # Gene-conversion recombinant planting (for the mosaic-read detector ground truth).
    ap.add_argument("--recomb-reads", type=int, default=0,
                    help="number of GENE-CONVERSION recombinant reads to plant: 5' exons from "
                         "copy A, 3' exons from copy B, switching at exon --recomb-exon. 0 = none.")
    ap.add_argument("--recomb-pair", type=str, default="0,1",
                    help="copies A,B for recombinants (5' from A, 3' from B), e.g. '0,1'.")
    ap.add_argument("--recomb-exon", type=int, default=2,
                    help="breakpoint exon index: exons [0,K) from copy A, [K,end) from copy B.")
    ap.add_argument("--invert-copies", type=str, default="",
                    help="comma-separated copy indices to place as INVERTED duplications "
                         "(revcomp genomic sequence). Reads from them align to the - strand; "
                         "tests whether VG discovery links an inverted paralog (e.g. DAZ1-/DAZ3+).")
    ap.add_argument("--segdup-flank", type=int, default=0,
                    help="bp of SHARED flanking sequence on EACH side of every copy's gene "
                         "(a true segmental duplication copies gene+flanks). 0 = bare paralog "
                         "(random, unique flanks). Ground truth for segdup-extent calling.")
    ap.add_argument("--hidden-copies", type=str, default="",
                    help="comma-separated copy indices PRESENT in the reads but ABSENT from the "
                         "reference genome (their locus is replaced by random). Models a copy "
                         "not in the assembly (collapsed segdup / CNV / private dup). Ground "
                         "truth for detecting copies-not-in-the-genome.")
    a = ap.parse_args()
    invert_set = set(int(x) for x in a.invert_copies.split(",") if x.strip() != "")
    hidden_set = set(int(x) for x in a.hidden_copies.split(",") if x.strip() != "")

    def revcomp(seq):
        comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return [comp[b] for b in reversed(seq)]
    rng = random.Random(a.seed)
    os.makedirs(a.out, exist_ok=True)

    # Gene model: 5 exons + 4 introns, ~14 kb genomic, ~1.8 kb mRNA.
    EXONS = [400, 300, 350, 300, 450]          # exon lengths (bp)
    INTRONS = [2600, 2300, 2500, 2400]         # intron lengths (bp)
    SPACING = a.spacing                        # copy-to-copy start distance
    chrom = "synTANDEM"

    # Per-copy exonic mutation rate p so two copies differ at ~2p(1-p) ≈ pairwise
    # divergence (1 - identity). p = (1 - identity) / 2.
    p = max(0.0, (1.0 - a.identity) / 2.0)
    exon_total = sum(EXONS)
    n_exon_snp = round(p * exon_total)

    # Ancestor genomic copy (one full gene: exon/intron/exon/.../exon).
    # Build the layout: list of (kind, length) and the genomic offset of each exon.
    layout = []
    for i, el in enumerate(EXONS):
        layout.append(("exon", el))
        if i < len(INTRONS):
            layout.append(("intron", INTRONS[i]))
    gene_len = sum(l for _, l in layout)
    ancestor = rand_seq(gene_len, rng)
    # exon position ranges within the gene (genomic, relative to gene start)
    exon_ranges = []
    off = 0
    for kind, l in layout:
        if kind == "exon":
            exon_ranges.append((off, off + l))
        off += l
    # flat list of exonic genomic positions (within the gene)
    exonic_positions = [pos for (s, e) in exon_ranges for pos in range(s, e)]
    # Intron ranges (between consecutive exons) + plant canonical GT...AG splice
    # motifs at every intron boundary. Without them minimap2 places splices ±1 bp
    # on random sequence, so assembled exon boundaries miss the truth GTF and
    # gffcompare scores 0. (Real genomes have canonical sites; this matches that.)
    intron_ranges = [(exon_ranges[i][1], exon_ranges[i + 1][0])
                     for i in range(len(exon_ranges) - 1)]
    for (istart, iend) in intron_ranges:
        ancestor[istart], ancestor[istart + 1] = 'G', 'T'   # donor
        ancestor[iend - 2], ancestor[iend - 1] = 'A', 'G'   # acceptor

    # Make N copies: each is the ancestor with its OWN exonic mutations.
    copies = []
    realized_pair_diffs = []
    for c in range(a.copies):
        seq = list(ancestor)
        snp_pos = rng.sample(exonic_positions, n_exon_snp) if n_exon_snp > 0 else []
        mutate_positions(seq, snp_pos, rng)
        copies.append(seq)
    # realized exonic pairwise identity (copy 0 vs 1, exonic positions only)
    if a.copies >= 2:
        diff = sum(1 for pos in exonic_positions if copies[0][pos] != copies[1][pos])
        realized_id = 1.0 - diff / exon_total
    else:
        realized_id = 1.0

    # Genome = copies laid out at SPACING, separated by random intergenic filler.
    # Each copy starts at c*SPACING (0-based) within the contig.
    contig_len = (a.copies - 1) * SPACING + gene_len + 5000 + 2 * a.segdup_flank
    genome = rand_seq(contig_len, rng)        # random intergenic background
    margin = max(1000, a.segdup_flank + 100)  # left margin must fit the upstream flank
    # Shared flanking sequence for a segmental duplication (same on every copy, then random).
    shared_up = rand_seq(a.segdup_flank, rng) if a.segdup_flank > 0 else []
    shared_down = rand_seq(a.segdup_flank, rng) if a.segdup_flank > 0 else []
    copy_starts = []
    for c in range(a.copies):
        start = c * SPACING + margin
        copy_starts.append(start)
        # Hidden copies: PRESENT in the reads (generated below) but ABSENT from the reference —
        # leave the locus as random intergenic background, so the copy is not in the assembly.
        # Its reads then mismap to the closest sibling, carrying their private SNPs as an extra
        # haplotype the reference copies cannot explain (the detect-copies-not-in-genome signal).
        if c in hidden_set:
            continue
        # Inverted copies: place the REVERSE-COMPLEMENT genomic sequence. The gene's mRNA
        # is unchanged (transcribed from the opposite strand), so its reads align to this
        # locus on the - strand — the inverted-duplication signature.
        genome[start:start + gene_len] = revcomp(copies[c]) if c in invert_set else copies[c]
        # Segmental duplication: plant SHARED flanks adjacent to the gene on each copy.
        if a.segdup_flank > 0:
            genome[start - a.segdup_flank:start] = shared_up
            genome[start + gene_len:start + gene_len + a.segdup_flank] = shared_down

    with open(os.path.join(a.out, "genome.fa"), "w") as fh:
        fh.write(f">{chrom}\n")
        s = "".join(genome)
        for i in range(0, len(s), 70):
            fh.write(s[i:i+70] + "\n")

    # Truth GTF: one full-length transcript per copy (the per-copy isoform), at
    # genomic coords. This is the ground truth for the annotation harness
    # (bench/tandem_benchmark.sh) — lets it score transcript-level Sn/Pr of VG
    # copy resolution on the synthetic, no external annotation needed.
    with open(os.path.join(a.out, "truth.gtf"), "w") as fh:
        for c in range(a.copies):
            cs = copy_starts[c]
            tid = f"copy{c}"
            inverted = c in invert_set
            strand = '-' if inverted else '+'
            for (es, ee) in exon_ranges:
                # Inverted copy: exon at forward gene-offset (es,ee) sits at the
                # revcomp position (gene_len-ee .. gene_len-es) on the - strand.
                if inverted:
                    g0 = cs + (gene_len - ee) + 1
                    g1 = cs + (gene_len - es)
                else:
                    g0 = cs + es + 1  # GTF is 1-based
                    g1 = cs + ee
                fh.write(
                    f'{chrom}\tsynth\texon\t{g0}\t{g1}\t.\t{strand}\t.\t'
                    f'gene_id "{tid}"; transcript_id "{tid}.1";\n'
                )

    # Reads: spliced mRNA from each copy (exons concatenated), with errors.
    def spliced_mrna(copy_seq):
        return [copy_seq[s:e] for (s, e) in exon_ranges]  # list of exon base-lists

    def add_errors(bases):
        out = []
        for b in bases:
            if rng.random() < a.error_rate:
                out.append(rng.choice([x for x in BASES if x != b]))
            else:
                out.append(b)
        return out

    reads_path = os.path.join(a.out, "reads.fa")
    n_reads = 0
    with open(reads_path, "w") as fh:
        for c in range(a.copies):
            exons_seq = spliced_mrna(copies[c])
            for j in range(a.reads_per_copy):
                # optionally 5'-truncate (drop leading exons) for realism
                first_exon = 0
                if rng.random() < a.trunc_frac and len(exons_seq) > 2:
                    first_exon = rng.randint(1, len(exons_seq) - 2)
                mrna = [b for ex in exons_seq[first_exon:] for b in ex]
                mrna = add_errors(mrna)
                seqstr = "".join(mrna)
                fh.write(f">c{c}_r{j}\n")
                for i in range(0, len(seqstr), 70):
                    fh.write(seqstr[i:i+70] + "\n")
                n_reads += 1

        # Gene-conversion recombinants: 5' exons from copy A, 3' exons from copy B,
        # switching at exon K. Full-length (no truncation) so both halves are observed.
        # Named recomb_{A}_{B}_e{K}_r{j} so the ground truth (mosaic + breakpoint) is
        # recoverable. The 5' tract carries A's exonic SNPs, the 3' tract carries B's —
        # so a read aligned to either copy's locus shows a per-diagnostic-site SWITCH.
        n_recomb = 0
        if a.recomb_reads > 0 and a.copies >= 2:
            ra, rb = (int(x) for x in a.recomb_pair.split(","))
            K = max(1, min(a.recomb_exon, len(EXONS) - 1))
            exons_a = spliced_mrna(copies[ra])
            exons_b = spliced_mrna(copies[rb])
            for j in range(a.recomb_reads):
                mrna = [b for ex in exons_a[:K] for b in ex] + \
                       [b for ex in exons_b[K:] for b in ex]
                mrna = add_errors(mrna)
                seqstr = "".join(mrna)
                fh.write(f">recomb_{ra}_{rb}_e{K}_r{j}\n")
                for i in range(0, len(seqstr), 70):
                    fh.write(seqstr[i:i+70] + "\n")
                n_recomb += 1
                n_reads += 1

    meta = {
        "recomb_reads": n_recomb,
        "hidden_copies": sorted(hidden_set) if hidden_set else None,
        "recomb_pair": a.recomb_pair if n_recomb else None,
        "recomb_exon": a.recomb_exon if n_recomb else None,
        "identity_target": a.identity,
        "identity_realized_exonic": round(realized_id, 5),
        "copies": a.copies,
        "reads_per_copy": a.reads_per_copy,
        "n_reads": n_reads,
        "n_exon_snp_per_copy": n_exon_snp,
        "exon_total_bp": exon_total,
        "error_rate": a.error_rate,
        "chrom": chrom,
        "copy_starts_1based": [s + 1 for s in copy_starts],
        "gene_len": gene_len,
        "seed": a.seed,
    }
    with open(os.path.join(a.out, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"[gen] identity_target={a.identity} realized_exonic={realized_id:.4f} "
          f"copies={a.copies} reads={n_reads} exon_snp/copy={n_exon_snp} -> {a.out}")

if __name__ == "__main__":
    main()
