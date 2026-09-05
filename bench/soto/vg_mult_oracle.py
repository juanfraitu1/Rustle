#!/usr/bin/env python3
"""Library-free VG minimizer-multiplicity repeat oracle, rebuilt for the HUMAN soto_adj reps.

⭐ WHY THIS IS NOT `RUSTLE_ER_REPEAT_MASKED_EDGES` (which §6cw refuted for inverting on duplication age):
RepeatMasker masks by IDENTITY TO A LIBRARY, so it deletes a young duplication's copy of an Alu along
with everything else -- and a recent segmental duplication is repeat-rich precisely BECAUSE the
duplication copied the repeats. Minimizer multiplicity instead masks by HOW MANY DISTINCT LOCI SHARE A
NODE. A 2-copy recent duplication's shared sequence has multiplicity 2, not 500, so it is NOT called a
repeat even when it lies inside an Alu. The inversion is avoided BY CONSTRUCTION, not by tuning.

⚠ MULTIPLICITY IS COUNTED OVER REPS (loci), NEVER OVER TRUTH GENES -- counting over truth labels would
leak the answer into the instrument.
"""
import sys, csv, statistics, collections

K, W = 15, 10   # minimizers.rs: MINIMIZER_K / MINIMIZER_W
_C = bytes.maketrans(b"ACGTacgtN", b"TGCAtgcaN")

def canon_minimizers(seq):
    s = seq.upper()
    n = len(s)
    if n < K: return set()
    kms = []
    for i in range(n - K + 1):
        km = s[i:i+K]
        if b"N" in km: kms.append(None); continue
        rc = km.translate(_C)[::-1]
        kms.append(min(km, rc))
    out = set()
    for i in range(len(kms) - W + 1):
        win = [x for x in kms[i:i+W] if x is not None]
        if win: out.add(min(win))
    return out

def load_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path, "rb"):
        if line.startswith(b">"):
            if name is not None: seqs[name] = b"".join(buf)
            name = line[1:].strip().decode(); buf = []
        else: buf.append(line.strip())
    if name is not None: seqs[name] = b"".join(buf)
    return seqs

if __name__ == "__main__":
    D = "/mnt/linuxdisk/home/juanfraitu/soto_adj/arm_off/dump"
    reps = load_fasta(f"{D}/e.er._k11_w5.0.reps.fa")
    print(f"reps: {len(reps)}", file=sys.stderr)
    mins = {int(k): canon_minimizers(v) for k, v in reps.items()}
    mult = collections.Counter()
    for s in mins.values():
        for m in s: mult[m] += 1          # multiplicity = # DISTINCT REPS carrying the node
    print(f"distinct minimizer nodes: {len(mult)}", file=sys.stderr)
    import pickle
    with open("/mnt/linuxdisk/home/juanfraitu/soto_adj/vg_mult.pkl", "wb") as fh:
        pickle.dump({"mins": mins, "mult": dict(mult)}, fh)
    print("wrote /mnt/linuxdisk/home/juanfraitu/soto_adj/vg_mult.pkl", file=sys.stderr)
    h = collections.Counter(min(v, 50) for v in mult.values())
    print("multiplicity histogram (capped at 50):")
    for m in sorted(h)[:12]: print(f"  mult={m:>3}  nodes={h[m]}")
    print(f"  mult>=20: {sum(1 for v in mult.values() if v>=20)}   mult>=50: {sum(1 for v in mult.values() if v>=50)}")
