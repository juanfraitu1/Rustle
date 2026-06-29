#!/usr/bin/env python3
"""Reference-only paralog identifiability ("malleability") score from minimizers.

Reuses rustle's EXACT minimizer scheme (vg_family::family_graph::minimizers):
FNV-1a 64-bit hash of each k-mer, minimum over each window of w consecutive
k-mers, k=15 w=10, forward strand. So a copy's discriminative-minimizer
fraction here equals what rustle could compute in-process for free.

Heuristic (validated, RABL2 vs DAZ):
  disc_frac(copy) = | minimizers(copy) NOT in ANY sibling copy | / | minimizers(copy) |
  HIGH disc_frac  -> copies are separable -> reads assignable     -> emit with confidence
  LOW  disc_frac  -> shared-sequence floor -> reads tie / MAPQ0   -> emitting risks fabrication
Anchors from the advisor work: RABL2 ~0.58 disc -> 0% tied reads (identifiable);
DAZ ~0.07 disc over a 2.9kb identical core -> ~19% tied (at the floor).

Usage:
  minimizer_identifiability.py family.fasta                  # each record = one copy
  minimizer_identifiability.py --genome G.fasta --bed copies.bed   # col4 = copy name
Outputs TSV: copy  n_minimizers  disc_frac  shared_frac  tier
"""
import sys, argparse

K, W = 15, 10
_FNV_OFF, _FNV_PRIME, _MASK = 14695981039346656037, 1099511628211, (1 << 64) - 1

def fnv1a(b: bytes) -> int:
    h = _FNV_OFF
    for byte in b:
        h = ((h ^ byte) * _FNV_PRIME) & _MASK
    return h

def minimizers(seq: bytes, k=K, w=W) -> set:
    out = set()
    if len(seq) < k: return out
    n = len(seq) - k + 1
    for ws in range(0, max(n - w, 1)):
        we = min(ws + w, n)
        best = None
        for i in range(ws, we):
            h = fnv1a(seq[i:i + k])
            best = h if best is None else min(best, h)
        if best is not None: out.add(best)
    return out

def tier(d: float) -> str:
    # bands anchored to validated RABL2 (0.58, identifiable) / DAZ (0.07, floor)
    if d >= 0.30:  return "IDENTIFIABLE"      # confident copy assignment
    if d >= 0.12:  return "BORDERLINE"        # needs strong read evidence
    return "SEQUENCE_FLOOR"                   # tied/non-identifiable; emitting = fabrication risk

def read_fasta(path):
    name, seq, recs = None, [], []
    for ln in open(path):
        if ln.startswith('>'):
            if name: recs.append((name, ''.join(seq).encode()))
            name, seq = ln[1:].split()[0], []
        else:
            seq.append(ln.strip())
    if name: recs.append((name, ''.join(seq).encode()))
    return recs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('family_fasta', nargs='?')
    ap.add_argument('--genome'); ap.add_argument('--bed')
    a = ap.parse_args()
    if a.genome and a.bed:
        import subprocess
        copies = []
        for ln in open(a.bed):
            f = ln.split('\t')
            if len(f) < 3: continue
            reg = f"{f[0]}:{int(f[1])+1}-{f[2]}"
            nm = f[3].strip() if len(f) > 3 else reg
            s = subprocess.run(['samtools', 'faidx', a.genome, reg],
                               capture_output=True, text=True).stdout
            seq = ''.join(l for l in s.splitlines() if not l.startswith('>')).encode()
            copies.append((nm, seq))
    elif a.family_fasta:
        copies = read_fasta(a.family_fasta)
    else:
        ap.error("need family.fasta OR --genome+--bed")
    mins = {n: minimizers(s) for n, s in copies}
    print("copy\tn_minimizers\tdisc_frac\tshared_frac\ttier")
    for n, _ in copies:
        own = mins[n]
        siblings = set().union(*[mins[o] for o in mins if o != n]) if len(mins) > 1 else set()
        uniq = len(own - siblings)
        d = uniq / len(own) if own else 0.0
        print(f"{n}\t{len(own)}\t{d:.3f}\t{1-d:.3f}\t{tier(d)}")

if __name__ == '__main__':
    main()
