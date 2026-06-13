#!/usr/bin/env python3
"""Simulate a 2-copy paralog family at KNOWN unequal per-copy depth, with a controlled number
of copy-distinguishing PSV sites, to test per-copy ABUNDANCE accuracy (VG fingerprint-EM vs 1/NH).

Two near-identical copies (copyA, copyB) differ by K substitutions placed in EXONS. Both copies
are placed at distinct genomic loci (chrA, chrB) with canonical GT..AG introns + unique flanks.
Full-length spliced transcript reads are simulated from each copy at depths Na/Nb (HiFi error e),
with the TRUE copy of origin encoded in the read name (rA_*/rB_*). Because the copies are
near-identical, EVERY read multimaps -> 1/NH always splits 50/50 (biased to equal), so the true
unequal split (Na:Nb) is the discriminator: only copy-distinguishing PSV evidence can recover it.

Usage: python3 sim_abundance_family.py <K_psv> <Na> <Nb> <outdir> [seed]
Writes <outdir>/genome.fa and <outdir>/reads.fq.
"""
import sys, os, random

def randseq(n, rng): return "".join(rng.choice("ACGT") for _ in range(n))
def mutate(b, rng): return rng.choice([x for x in "ACGT" if x != b])

def main():
    K   = int(sys.argv[1])      # number of copy-distinguishing PSV sites (in exons)
    Na  = int(sys.argv[2])      # copyA depth (reads)
    Nb  = int(sys.argv[3])      # copyB depth (reads)
    out = sys.argv[4]
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 7
    rng = random.Random(seed)
    os.makedirs(out, exist_ok=True)
    ERR = 0.005                 # HiFi per-base substitution error

    exon_lens   = [300, 200, 250, 350]
    intron_lens = [500, 600, 550]
    flank = 1000

    # ---- ancestral (copyA) exon + intron sequences ----
    exonsA = [list(randseq(L, rng)) for L in exon_lens]
    intronsA = []
    for L in intron_lens:
        s = list(randseq(L, rng)); s[0], s[1] = "G", "T"; s[-2], s[-1] = "A", "G"  # GT..AG
        intronsA.append(s)

    # ---- copyB = copyA with K PSVs placed across exons ----
    exonsB = [e[:] for e in exonsA]
    # choose K (exon_index, pos) sites spread across exons
    all_sites = [(ei, p) for ei, e in enumerate(exonsB) for p in range(len(e))]
    rng.shuffle(all_sites)
    psv_sites = all_sites[:K]
    for (ei, p) in psv_sites:
        exonsB[ei][p] = mutate(exonsB[ei][p], rng)
    intronsB = [s[:] for s in intronsA]  # introns identical (PSVs only in exons)

    def genomic(exons, introns):
        seq = []
        for i, e in enumerate(exons):
            seq += e
            if i < len(introns): seq += introns[i]
        return "".join(seq)
    def transcript(exons): return "".join("".join(e) for e in exons)

    flankA1, flankA2 = randseq(flank, rng), randseq(flank, rng)
    flankB1, flankB2 = randseq(flank, rng), randseq(flank, rng)
    chrA = flankA1 + genomic(exonsA, intronsA) + flankA2
    chrB = flankB1 + genomic(exonsB, intronsB) + flankB2
    with open(f"{out}/genome.fa", "w") as f:
        f.write(f">chrA\n{chrA}\n>chrB\n{chrB}\n")

    txA, txB = transcript(exonsA), transcript(exonsB)
    def emit(tx, name, fq):
        # realistic 5'/3' positional diversity (IsoSeq-like) so reads aren't degenerate-identical
        a = rng.randint(0, 250); b3 = rng.randint(0, 150)
        sub = tx[a: len(tx) - b3]
        read = [mutate(c, rng) if rng.random() < ERR else c for c in sub]
        s = "".join(read)
        fq.write(f"@{name}\n{s}\n+\n{'I'*len(s)}\n")
    with open(f"{out}/reads.fq", "w") as fq:
        for i in range(Na): emit(txA, f"rA_{i}", fq)
        for i in range(Nb): emit(txB, f"rB_{i}", fq)

    print(f"K={K} PSVs  Na={Na} Nb={Nb}  true_splitA={Na/(Na+Nb):.3f}  "
          f"txlen={len(txA)}  wrote {out}/genome.fa,{out}/reads.fq")

if __name__ == "__main__":
    main()
