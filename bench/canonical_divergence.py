#!/usr/bin/env python3
"""3-way canonical-vs-default-vs-StringTie chain-set divergence classifier.
Usage: python3 bench/canonical_divergence.py <default.gtf> <canonical.gtf> <st.gtf> <annot.gff> <chrom>
Reports the final chains canonical adds/removes vs default and vs ST, and their TP/FP nature.
"""
import re, sys, collections

def chains_gtf(path):
    tx = collections.defaultdict(list); strand = {}
    for line in open(path):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]), int(f[4]))); strand[m.group(1)] = f[6]
    out = set()
    for t, ex in tx.items():
        ex.sort()
        if len(ex) < 2: continue
        out.add((strand[t], tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))))
    return out

def chains_gff(path, chrom):
    tx = collections.defaultdict(list); strand = {}
    for line in open(path):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] != 'exon': continue
        m = re.search(r'Parent=([^;]+)', f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]), int(f[4]))); strand[m.group(1)] = f[6]
    out = set()
    for t, ex in tx.items():
        ex.sort()
        if len(ex) < 2: continue
        out.add((strand[t], tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))))
    return out

dft, can, stp, annot, chrom = sys.argv[1:6]
D = chains_gtf(dft); C = chains_gtf(can); S = chains_gtf(stp); A = chains_gff(annot, chrom)

def tpfp(chs):
    tp = sum(1 for k in chs if k in A); return tp, len(chs) - tp

added = C - D          # chains canonical introduced (not in default)
removed = D - C        # chains canonical dropped (were in default)
can_only = C - S       # canonical's rustle-only (the 223)
def_only = D - S       # default's rustle-only (the 187)
print(f"default tx-chains={len(D)} canonical tx-chains={len(C)} ST tx-chains={len(S)} annot={len(A)}")
print(f"canonical rustle-only(vs ST)={len(can_only)}  default rustle-only={len(def_only)}  delta={len(can_only)-len(def_only)}")
print()
for name, chs in [("ADDED by canonical (in C, not D)", added),
                  ("REMOVED by canonical (in D, not C)", removed)]:
    tp, fp = tpfp(chs)
    in_st = sum(1 for k in chs if k in S)
    print(f"{name}: n={len(chs)}  TP(RefSeq)={tp} FP={fp}  in_ST_final={in_st}")
# The decisive split: of the chains canonical ADDS, how many are ST-shared (good: matching ST)
# vs canonical-only (bad: new FPs ST does not have)?
add_st_shared = sum(1 for k in added if k in S)
add_canon_only_fp = sum(1 for k in added if k not in S and k not in A)
add_canon_only_tp = sum(1 for k in added if k not in S and k in A)
print()
print(f"  ADDED breakdown: ST-shared(converging)={add_st_shared}  "
      f"canonical-only-FP(regressing)={add_canon_only_fp}  canonical-only-TP(real)={add_canon_only_tp}")
# And the removed-shared (recall lost vs ST): chains canonical dropped that ST has
rem_st_shared = sum(1 for k in removed if k in S)
rem_tp = sum(1 for k in removed if k in A)
print(f"  REMOVED breakdown: ST-shared(recall lost)={rem_st_shared}  RefSeq-TP(real lost)={rem_tp}")
# Reconciliation assertion: canonical_only == default_only + added_not_in_S - removed_not_in_S
lhs = len(can_only)
rhs = len(def_only) + len([k for k in added if k not in S]) - len([k for k in removed if k not in S])
assert lhs == rhs, f"chain-set arithmetic inconsistent: can_only={lhs} != reconstructed={rhs}"
print("\n[OK] chain-set arithmetic reconciles")
