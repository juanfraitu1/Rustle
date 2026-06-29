#!/usr/bin/env python3
"""Generate two synthetic 2-copy paralog benchmarks that differ by STRUCTURAL
variation (not just SNPs), to test whether sequence-to-graph alignment beats
linear seed-chaining for paralog copy assignment.

Single-exon transcript copies (no splicing needed). Each copy ~1.5-3 kb.
IsoSeq-like reads (~99% identity to their copy, a few % residual error), read
NAMES encode the TRUE copy (uniqA_* / uniqB_*) and the read CLASS
(full / span / nospan). Reads are emitted as plain FASTA (read seq = its copy's
sequence with errors), so ground truth is fully derivable from the name.

Regime 1 -- INVERSION:
  Copy A and copy B identical EXCEPT a ~200 bp internal segment that is
  REVERSE-COMPLEMENTED in copy B relative to copy A, plus 2 flanking SNPs far
  from the inversion. Reads spanning the inverted segment + flanks should force
  linear minimap2 into split/supplementary or soft-clipped alignment when placed
  on the wrong copy; a graph (inversion = fwd/rev bubble) should align cleanly.

Regime 2 -- LARGE INDEL + sparse SNPs:
  Copies differ by a ~300 bp INSERTION present in copy A, absent in copy B, plus
  only 2 SNPs elsewhere. Reads that span vs don't span the indel are labelled.
  A read that does NOT span the indel AND carries 0 SNPs is info-theoretically
  unresolvable by ANY method -- we label these 'nospan_nosnp' and report them
  separately (excluded from the headline resolvable accuracy).
"""
import random
import os

random.seed(1234)

OUTDIR = os.path.dirname(os.path.abspath(__file__))

ALPH = "ACGT"
COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def rand_seq(n):
    return "".join(random.choice(ALPH) for _ in range(n))


def rc(s):
    return "".join(COMP[b] for b in reversed(s))


def mutate(seq, err=0.01):
    """Apply ~err per-base error (sub/ins/del mix) to simulate CCS residual error.
    Returns mutated sequence. ~99% identity at err=0.01."""
    out = []
    for b in seq:
        r = random.random()
        if r < err * 0.7:  # substitution (70% of errors)
            out.append(random.choice([x for x in ALPH if x != b]))
        elif r < err * 0.85:  # insertion (15%)
            out.append(b)
            out.append(random.choice(ALPH))
        elif r < err:  # deletion (15%)
            continue
        else:
            out.append(b)
    return "".join(out)


def write_fasta(path, records):
    with open(path, "w") as f:
        for name, seq in records:
            f.write(">%s\n" % name)
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


# ============================================================================
# REGIME 1 -- INVERSION
# ============================================================================
def build_inversion():
    d = os.path.join(OUTDIR, "inversion")
    # Layout: [left flank ~800] [inv segment ~200] [right flank ~800]  => ~1800 bp
    LEFT = 800
    INV = 200
    RIGHT = 800
    left = rand_seq(LEFT)
    inv = rand_seq(INV)
    right = rand_seq(RIGHT)

    # Copy A: left + inv + right
    copyA = list(left + inv + right)
    # Copy B: left + rc(inv) + right  (the inversion polymorphism)
    copyB_inv = rc(inv)
    copyB = list(left + copyB_inv + right)

    # Add 2 flanking SNPs far from the inversion (one in left flank, one in
    # right flank) so that reads NOT spanning the inversion but spanning a SNP
    # are still resolvable, and so the inversion is not the ONLY signal.
    # SNP1 at left offset 200; SNP2 at right offset (LEFT+INV+600).
    snp_positions = [200, LEFT + INV + 600]
    for p in snp_positions:
        a = copyA[p]
        b = next(x for x in ALPH if x != a)
        copyA[p] = a
        copyB[p] = b
    copyA = "".join(copyA)
    copyB = "".join(copyB)
    assert len(copyA) == len(copyB) == LEFT + INV + RIGHT

    inv_start = LEFT
    inv_end = LEFT + INV  # half-open
    write_fasta(os.path.join(d, "copies.fa"),
                [("copyA", copyA), ("copyB", copyB)])

    # Coordinates / metadata for the assignment scripts.
    with open(os.path.join(d, "meta.txt"), "w") as f:
        f.write("LEFT=%d\nINV=%d\nRIGHT=%d\n" % (LEFT, INV, RIGHT))
        f.write("inv_start=%d\ninv_end=%d\n" % (inv_start, inv_end))
        f.write("snp_positions=%s\n" % ",".join(map(str, snp_positions)))

    # ----- Reads -----
    reads = []
    total_len = len(copyA)

    def emit(copy_name, true_copy, cls, idx, start, end):
        src = copyA if true_copy == "A" else copyB
        sub = src[start:end]
        seq = mutate(sub, err=0.012)
        name = "uniq%s_%s_%03d_s%d_e%d" % (true_copy, cls, idx, start, end)
        reads.append((name, seq))

    # full-length reads: span the whole copy (definitely span inversion)
    for i in range(15):
        emit("A", "A", "full", i, 0, total_len)
    for i in range(15):
        emit("B", "B", "full", i, 0, total_len)

    # inversion-spanning reads: cover [inv_start-flank .. inv_end+flank]
    # with generous flanks on both sides so each half can anchor.
    for i in range(15):
        flankL = random.randint(120, 300)
        flankR = random.randint(120, 300)
        s = max(0, inv_start - flankL)
        e = min(total_len, inv_end + flankR)
        emit("A", "A", "span", i, s, e)
    for i in range(15):
        flankL = random.randint(120, 300)
        flankR = random.randint(120, 300)
        s = max(0, inv_start - flankL)
        e = min(total_len, inv_end + flankR)
        emit("B", "B", "span", i, s, e)

    # NON-spanning reads that DO carry a flanking SNP (resolvable, control):
    # a left-flank read covering SNP1 (pos 200) but not the inversion.
    for i in range(8):
        s = 0
        e = random.randint(400, 700)  # < inv_start=800, covers SNP1@200
        emit("A", "A", "nospan_snp", i, s, e)
    for i in range(8):
        s = 0
        e = random.randint(400, 700)
        emit("B", "B", "nospan_snp", i, s, e)

    write_fasta(os.path.join(d, "reads.fa"), reads)
    return len(reads)


# ============================================================================
# REGIME 2 -- LARGE INDEL + sparse SNPs
# ============================================================================
def build_indel():
    d = os.path.join(OUTDIR, "indel")
    # Backbone shared by both copies: left ~1000 + right ~1000.
    # Copy A additionally has a ~300 bp INSERTION between left and right.
    # Copy B has no insertion. => A ~2300 bp, B ~2000 bp.
    LEFT = 1000
    RIGHT = 1000
    INS = 300
    left = rand_seq(LEFT)
    right = rand_seq(RIGHT)
    insertion = rand_seq(INS)

    copyA = list(left + insertion + right)   # has insertion
    copyB = list(left + right)               # lacks insertion

    # Only 2 SNPs elsewhere: one in left flank (pos 300), one in right flank.
    # Right-flank SNP in copyA coords = LEFT+INS+700; in copyB coords = LEFT+700.
    snpL = 300
    a = copyA[snpL]
    b = next(x for x in ALPH if x != a)
    copyA[snpL] = a
    copyB[snpL] = b
    rightA = LEFT + INS + 700
    rightB = LEFT + 700
    a2 = copyA[rightA]
    b2 = next(x for x in ALPH if x != a2)
    copyA[rightA] = a2
    copyB[rightB] = b2

    copyA = "".join(copyA)
    copyB = "".join(copyB)
    ins_start = LEFT
    ins_end = LEFT + INS  # in copyA coords

    write_fasta(os.path.join(d, "copies.fa"),
                [("copyA", copyA), ("copyB", copyB)])
    with open(os.path.join(d, "meta.txt"), "w") as f:
        f.write("LEFT=%d\nRIGHT=%d\nINS=%d\n" % (LEFT, RIGHT, INS))
        f.write("ins_start=%d\nins_end=%d\n" % (ins_start, ins_end))
        f.write("snpL_A=%d\nsnpL_B=%d\n" % (snpL, snpL))
        f.write("snpR_A=%d\nsnpR_B=%d\n" % (rightA, rightB))
        f.write("lenA=%d\nlenB=%d\n" % (len(copyA), len(copyB)))

    reads = []

    def emit(true_copy, cls, idx, start, end):
        src = copyA if true_copy == "A" else copyB
        sub = src[start:end]
        seq = mutate(sub, err=0.012)
        name = "uniq%s_%s_%03d_s%d_e%d" % (true_copy, cls, idx, start, end)
        reads.append((name, seq))

    # full-length reads (span insertion locus for A; full backbone for B)
    for i in range(12):
        emit("A", "full", i, 0, len(copyA))
    for i in range(12):
        emit("B", "full", i, 0, len(copyB))

    # indel-spanning reads: cover [ins_start-flank .. ins_end+flank] (copyA),
    # and for copyB the homologous junction region [LEFT-flank .. LEFT+flank].
    for i in range(12):
        fl = random.randint(150, 300)
        frr = random.randint(150, 300)
        s = max(0, ins_start - fl)
        e = min(len(copyA), ins_end + frr)
        emit("A", "span", i, s, e)
    for i in range(12):
        fl = random.randint(150, 300)
        frr = random.randint(150, 300)
        s = max(0, LEFT - fl)
        e = min(len(copyB), LEFT + frr)
        emit("B", "span", i, s, e)

    # non-spanning reads that DO carry a SNP (resolvable control): left-flank
    # read covering snpL@300 but not the insertion locus.
    for i in range(6):
        s = 0
        e = random.randint(500, 800)  # covers snpL@300, < ins_start=1000
        emit("A", "nospan_snp", i, s, e)
    for i in range(6):
        s = 0
        e = random.randint(500, 800)
        emit("B", "nospan_snp", i, s, e)

    # non-spanning reads carrying NO SNP (info-theoretically UNRESOLVABLE):
    # a middle chunk of the right flank that contains neither the insertion
    # junction nor either SNP. These are reported separately.
    # right flank region in A: [LEFT+INS .. LEFT+INS+RIGHT]; avoid rightA=LEFT+INS+700.
    for i in range(6):
        # copyA right flank, offset 50..550 (before SNP@700)
        s = LEFT + INS + 50
        e = s + random.randint(300, 450)
        emit("A", "nospan_nosnp", i, s, e)
    for i in range(6):
        s = LEFT + 50
        e = s + random.randint(300, 450)
        emit("B", "nospan_nosnp", i, s, e)

    write_fasta(os.path.join(d, "reads.fa"), reads)
    return len(reads)


if __name__ == "__main__":
    n1 = build_inversion()
    n2 = build_indel()
    print("inversion: %d reads" % n1)
    print("indel:     %d reads" % n2)
