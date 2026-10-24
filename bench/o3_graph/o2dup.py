#!/usr/bin/env python3
"""Milestone 3, P9 substrate: O2 on fresh gorilla fibroblast families with one synthetic near-identical duplicate (2026-10-23).
The fresh >= 3-copy families hold no AS-tied reads with their real copies (tie scan: 0 in every one), so ties are made: copy X (the
family's copy whose two alleles differ most, first by chrom, start on ties; if they differ at < 2 positions its B allele gets S_HET
substitutions, seed 2, as in the one-copy sweep) is duplicated as X' = X's locus +- DUP_FLANK with D substitutions inside X's exons (seed 1), on
its own reference contig. X' is homozygous (its A transcript on both alleles). Reads as in o2sim.py (N per allele per copy, seed 1).
  - case `control`: the sample has every copy and X'; X' is in the reference and the copy table.
  - case `missing`: the sample lacks X' (no X' reads); X' stays in the reference and the copy table, so X's reads tie onto a copy the
    sample does not have. Any X read placed on X' is wrong.
The port (default options) runs; AS-tied reads are scored by o2eval_sim.py.
    o2dup.py OUT_DIR FAM[,FAM...] D[,D...] [N]"""
import collections, os, random, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import o3sim_one as S
import pysam

DUP_FLANK = 2000
S_HET = 6
CONV = os.environ.get("O2DUP_CONV", "0") == "1"


def duplicate(pri, c, blocks, strand, d, rng):
    """(contig sequence, X' transcript, (start, end) of X' on the contig): X's locus +- DUP_FLANK with d substitutions at exonic
    positions."""
    lo = max(0, c["start"] - DUP_FLANK); hi = min(pri.get_reference_length(c["chrom"]), c["end"] + DUP_FLANK)
    seq = list(pri.fetch(c["chrom"], lo, hi).upper())
    exonic = [p for x, y in blocks for p in range(x, y) if lo <= p < hi and seq[p - lo] in "ACGT"]
    for p in rng.sample(exonic, d):
        seq[p - lo] = rng.choice([b for b in "ACGT" if b != seq[p - lo]])
    t = "".join("".join(seq[x - lo:y - lo]) for x, y in blocks)
    return "".join(seq), (S.rc(t) if strand == "-" else t), (c["start"] - lo, c["end"] - lo)


def main():
    out, fams, ds = sys.argv[1], sys.argv[2].split(","), [int(x) for x in sys.argv[3].split(",")]
    N = int(sys.argv[4]) if len(sys.argv) > 4 else 100
    by = collections.defaultdict(list)
    for r in S.tsv(S.SRC): by[r["fam"]].append(r)
    cm = {r["pri"]: r for r in S.tsv(f"{S.RA}/chrmap.tsv")}
    pri, bam = pysam.FastaFile(S.PRI), pysam.AlignmentFile(S.FIBRO)
    for f in fams:
        copies, tx, iso = {}, {}, {}
        for r in sorted(by[f], key=lambda r: (r["chrom"], int(r["clean_start"]))):
            i = S.dominant_isoform(bam, r["chrom"], int(r["clean_start"]), int(r["clean_end"]))
            if not i: continue
            t = S.transcripts(pri, cm, r["chrom"], i[0], i[1])
            copies[r["gene"]] = dict(chrom=r["chrom"], start=int(r["clean_start"]), end=int(r["clean_end"]))
            iso[r["gene"]] = i
            tx[f"{r['gene']}|A"] = t["A"]
            if t["B"]: tx[f"{r['gene']}|B"] = t["B"]
        het = {}
        for c in copies:
            v = S.compare(tx[f"{c}|A"], tx[f"{c}|B"], out) if f"{c}|B" in tx else None
            het[c] = v[1] + v[2] if v else -1
        X = max(copies, key=lambda c: het[c]); Xd = f"{X}d"
        if het[X] < 2:
            rng, sq = random.Random(2), list(tx[f"{X}|A"])
            for i in rng.sample(range(len(sq)), S_HET):
                sq[i] = rng.choice([b for b in "ACGT" if b != sq[i]])
            tx[f"{X}|B"] = "".join(sq)
        print(f"{f}: donor {X}, allele differences {het[X]}{' -> synthetic B' if het[X] < 2 else ''}", flush=True)
        for d in ds:
            contig, txd, (s0, e0) = duplicate(pri, copies[X], iso[X][0], iso[X][1], d, random.Random(1))
            txx = dict(tx)
            if CONV and d:   # a conversion in the sample: X's B allele takes X''s base at its first substitution
                i = next(k for k in range(len(txd)) if txd[k] != tx[f"{X}|A"][k])
                ctx = tx[f"{X}|A"][i - 15:i]; j = tx[f"{X}|B"].find(ctx)
                if i >= 15 and j >= 0:
                    b = tx[f"{X}|B"]; txx[f"{X}|B"] = b[:j + 15] + txd[i] + b[j + 16:]
                    print(f"{f} d{d}: X|B converted at transcript {i}", flush=True)
            for case in ("control", "missing"):
                D = f"{out}/{f}_d{d}_{case}"; os.makedirs(D, exist_ok=True)
                where = S.mini_reference(copies, None, False, D)
                with open(f"{D}/ref.fa", "a") as fh: fh.write(f">dup_{Xd}\n{contig}\n")
                os.remove(f"{D}/ref.fa.fai"); pysam.faidx(f"{D}/ref.fa")
                where[Xd] = (f"dup_{Xd}", s0, e0)
                txs = dict(txx)
                plan = [(t, N) for t in txx]
                if case == "control":
                    txs[f"{Xd}|A"] = txd; plan.append((f"{Xd}|A", 2 * N))
                S.simulate(txs, plan, random.Random(1), f"{D}/reads.fa")
                S.align(D, f"{D}/reads.fa", f"{D}/reads.bam")
                with open(f"{D}/copies.tsv", "w") as fh:
                    fh.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n")
                    for k, c in enumerate(where):
                        seg, a, b = where[c]; fh.write(f"{f}\t{k}\t{c}\t{seg}\t{a}\t{b}\t1\t+\t0\t{a}-{b}\n")
                names = [l[1:].strip() for l in open(f"{D}/reads.fa") if l.startswith(">")]
                with open(f"{D}/net.tsv", "w") as fh:
                    fh.write("family\tread\tused\trole\tvia\tclass\tin_reads_tsv\n")
                    for n in names: fh.write(f"{f}\t{n}\t1\tS\tA\tNA\t0\n")
                subprocess.run([sys.executable, f"{HERE}/o3port.py", "run", "--bam", f"{D}/reads.bam", "--copies", f"{D}/copies.tsv",
                                "--genome", f"{D}/ref.fa", "--net", f"{D}/net.tsv", "--out", f"{D}/port"], check=True,
                               stderr=subprocess.DEVNULL)
                print(subprocess.run([sys.executable, f"{HERE}/o2eval_sim.py", D], capture_output=True, text=True).stdout.strip(),
                      flush=True)


if __name__ == "__main__":
    main()
