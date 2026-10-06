"""IsoSeq-like reads per copy and isoform, with the truth in the names and a per-read truth table.

Read = one isoform's spliced chain (transcript orientation, i.e. what an FLNC read is) -> optional 5' truncation ->
MANDATORY end jitter (identical reads collapse under the assembler's dedup, memory §6n0) -> HiFi errors
(`sim.simulate_reads`: substitutions + rare short indels). Names: `copy|isoform|i`. The truth table records, per read,
the chain interval it covers and the genomic junctions it contains (0-based half-open introns), so a scorer never
credits a read for a junction its truncation removed.
"""
import os
import random

from sim import simulate_reads, write_fastq  # noqa: E402  (bench/sim.py)

from .model import seeded

TRUTH_HEADER = "read\tcopy\tkind\tfamily\tisoform\tcontig\tstrand\tchain_len\tchain_start\tchain_end\ttruncated5\tn_junctions\tjunctions\n"


def simulate(planted, spec, out_dir, seed, log=print):
    """Write DIR/reads.fq and DIR/reads.truth.tsv. `spec` = the scenario's `reads` block."""
    per_copy = int(spec.get("per_copy", 50))
    err, indel = float(spec.get("err", 0.001)), float(spec.get("indel", 0.0003))
    jit = int(spec.get("jitter", 30))
    t5_frac, t5_max = float(spec.get("trunc5_frac", 0.0)), float(spec.get("trunc5_max", 0.3))
    if jit <= 0:
        raise ValueError("jitter must be > 0: identical reads collapse under dedup (§6n0)")
    n_total = 0
    counts = {}
    with open(os.path.join(out_dir, "reads.fq"), "w") as fq, open(os.path.join(out_dir, "reads.truth.tsv"), "w") as tr:
        tr.write(TRUTH_HEADER)
        for p in planted:
            depth = per_copy if p.expression is None else int(p.expression)
            counts[p.id] = 0
            if depth <= 0:
                continue
            isos = p.isoform_list()
            wsum = sum(w for _, _, w in isos)
            for iso_name, skip, w in isos:
                n = int(round(depth * w / wsum))
                chain = p.model.chain_seq(skip)
                ex = p.model.rna_exons(skip)
                if n == 0 or len(chain) < 100:
                    continue
                # chain-local offsets of each exon, to map a read's chain interval to the junctions it contains
                offs, o = [], 0
                for e in ex:
                    offs.append((o, o + len(e), e)); o += len(e)
                juncs = [(offs[k][1], p.genomic(ex[k].end, ex[k + 1].start)) for k in range(len(ex) - 1)]
                for i in range(n):
                    rng = seeded(seed, "read", p.id, iso_name, i)
                    lo, hi = 0, len(chain)
                    trunc = 0
                    if t5_frac > 0 and rng.random() < t5_frac:
                        trunc = int(rng.random() * t5_max * len(chain)); lo += trunc
                    lo += rng.randint(0, jit); hi -= rng.randint(0, jit)
                    if hi - lo < 100:
                        lo, hi = 0, len(chain)
                    body = chain[lo:hi]
                    rd, q = simulate_reads(body, 1, err=err, indel=indel, seed=rng.getrandbits(30))[0]
                    name = f"{p.id}|{iso_name}|{i}"
                    write_fastq(fq, name, (rd, q))
                    inside = [g for cpos, g in juncs if lo + 1 <= cpos <= hi - 1]   # >= 1 bp on both sides
                    tr.write(f"{name}\t{p.id}\t{p.kind}\t{p.family}\t{iso_name}\t{p.contig}\t{p.strand}\t{len(chain)}\t{lo}\t{hi}\t{trunc}\t"
                             f"{len(inside)}\t{','.join(f'{s}-{e}' for s, e in sorted(inside))}\n")
                    counts[p.id] += 1; n_total += 1
    log(f"reads: {n_total} in reads.fq (" + ", ".join(f"{k} {v}" for k, v in counts.items()) + ")")
    return counts


def load_truth(out_dir):
    """{read: dict} from reads.truth.tsv (junctions as a set of (s, e))."""
    out = {}
    with open(os.path.join(out_dir, "reads.truth.tsv")) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.rstrip("\n").split("\t")))
            f["junctions"] = {tuple(int(x) for x in j.split("-")) for j in f["junctions"].split(",") if j}
            for k in ("chain_len", "chain_start", "chain_end", "truncated5", "n_junctions"):
                f[k] = int(f[k])
            out[f["read"]] = f
    return out
