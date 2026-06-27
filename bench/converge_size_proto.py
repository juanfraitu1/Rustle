#!/usr/bin/env python3
"""converge_size_proto.py — measure how prevalent spurious terminal tails are, to decide whether the
cross-copy "true-size convergence" refinement is worth implementing.

For each family, all-vs-all align the copies' GENOMIC spans (where read-driven length inflation lives) and,
for each copy, compute the fraction of its length that aligns to NO sibling — the candidate spurious /
copy-specific terminal tail. Also classify each currently-FAILING (DNA-silent) family: is its failure
driven by a terminal TAIL (one copy has a big unaligned end → convergence would help) or by INTERNAL /
intron-length divergence (the aligned core itself is short → convergence does NOT help, it's real divergence).

Run: python bench/converge_size_proto.py <refined_prefix> <genome> [samtools] [minimap2]
"""
import sys, csv, subprocess, collections, os, tempfile

PREFIX = sys.argv[1] if len(sys.argv) > 1 else "/home/juanfra/winloci_scratch/gw_xchrom_refined2"
GENOME = sys.argv[2] if len(sys.argv) > 2 else "/home/juanfra/winloci_scratch/GGO.fasta"
SAM = sys.argv[3] if len(sys.argv) > 3 else "/home/juanfra/miniforge3/envs/biser/bin/samtools"
MM2 = sys.argv[4] if len(sys.argv) > 4 else "/home/juanfra/miniforge3/bin/minimap2"
ID, COV = 0.90, 0.50

copies = collections.defaultdict(list)
for c in csv.DictReader(open(f"{PREFIX}.copies.tsv"), delimiter="\t"):
    copies[c["family_id"]].append((c["chrom"], int(c["start"]), int(c["end"]), int(c["n_exon"]), c["strand"]))
fams = {f["family_id"]: f for f in csv.DictReader(open(f"{PREFIX}.families.tsv"), delimiter="\t")}

def gspans(cps):
    regs = [f"{c[0]}:{c[1]+1}-{c[2]}" for c in cps]
    out = subprocess.run([SAM, "faidx", GENOME, *regs], capture_output=True, text=True).stdout
    seqs, cur = [], []
    for ln in out.splitlines():
        if ln.startswith(">"):
            if cur: seqs.append("".join(cur)); cur = []
        else: cur.append(ln)
    if cur: seqs.append("".join(cur))
    return seqs

def merge(ivs):
    if not ivs: return []
    ivs = sorted(ivs); out = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= out[-1][1]: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return out

def analyze(cps):
    seqs = gspans(cps)
    n = len(cps)
    if len(seqs) != n: return None
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        path = fh.name
        for i, s in enumerate(seqs): fh.write(f">{i}\n{s}\n")
    try:
        out = subprocess.run([MM2, "-cx", "asm20", "-X", "--no-long-join", path, path],
                             capture_output=True, text=True).stdout
    finally:
        os.unlink(path)
    lens = [len(s) for s in seqs]
    qaln = collections.defaultdict(list)   # copy -> aligned intervals (any sibling, identity-passing)
    pair_block = {}                         # (i,j) -> max aligned query length at id>=ID
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 12: continue
        try: q, t = int(f[0]), int(f[5])
        except ValueError: continue
        if q == t: continue
        qs, qe, ts, te = int(f[2]), int(f[3]), int(f[7]), int(f[8])
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
        ident = (1 - de) if de is not None else int(f[9]) / max(1, int(f[10]))
        if ident < ID: continue
        qaln[q].append((qs, qe)); qaln[t].append((ts, te))
        k = (min(q, t), max(q, t))
        pair_block[k] = max(pair_block.get(k, 0), qe - qs)
    # per-copy covered length (by any sibling) and unaligned-tail fraction
    covered = []
    for i in range(n):
        cov_len = sum(e - s for s, e in merge(qaln[i]))
        covered.append(cov_len)
    tail_frac = [1 - covered[i] / max(1, lens[i]) for i in range(n)]
    # confirmation variants. raw: cov over min(raw length). naive-core: cov over min(aligned-core) — UNSAFE
    # (a tiny shared fragment gives cov=1). safe-core: denom floored at FLOOR*raw so a mostly-unaligned copy
    # (domain/repeat-sharer) cannot be "converged" to its fragment; only small spurious TAILS are forgiven.
    FLOOR = 0.70
    def confirmed(mode):
        par = list(range(n))
        def find(x):
            while par[x] != x: par[x] = par[par[x]]; x = par[x]
            return x
        for (i, j), blk in pair_block.items():
            if mode == "raw":
                denom = min(lens[i], lens[j])
            elif mode == "core":
                denom = min(covered[i], covered[j])
            else:  # safe
                denom = min(max(covered[i], FLOOR * lens[i]), max(covered[j], FLOOR * lens[j]))
            if denom > 0 and blk / denom >= COV:
                par[find(i)] = find(j)
        return max(collections.Counter(find(i) for i in range(n)).values(), default=0) >= 2
    return tail_frac, confirmed("raw"), confirmed("core"), confirmed("safe"), lens, covered

all_tail = []
raw_conf = core_conf = safe_conf = n = 0
core_rec, safe_rec = [], []   # families confirmed by naive-core / safe-core but not raw
for fid, cps in copies.items():
    r = analyze(cps)
    if not r: continue
    tail_frac, rc, cc, sc, lens, covered = r
    n += 1; raw_conf += rc; core_conf += cc; safe_conf += sc
    all_tail += tail_frac
    if cc and not rc:
        core_rec.append((fid, fams[fid]["cross_chrom"], len(cps), [round(t, 2) for t in tail_frac]))
    if sc and not rc:
        safe_rec.append((fid, fams[fid]["cross_chrom"], len(cps), [round(t, 2) for t in tail_frac]))

import statistics
all_tail.sort()
print(f"=== convergence relevance on {PREFIX} ({n} families) ===\n")
print(f"per-copy UNALIGNED-tail fraction (length aligning to NO sibling at id>=0.90):")
print(f"  median={statistics.median(all_tail):.3f}  mean={statistics.mean(all_tail):.3f}  "
      f"p90={all_tail[int(0.9*len(all_tail))]:.3f}  max={max(all_tail):.3f}")
big = sum(1 for t in all_tail if t > 0.30)
print(f"  copies with >30% unaligned tail: {big}/{len(all_tail)} = {100*big/len(all_tail):.1f}%")
print(f"\nDNA-confirmed (raw cov-of-shorter):           {raw_conf}/{n}")
print(f"DNA-confirmed (NAIVE cov-of-core, UNSAFE):    {core_conf}/{n}   (+{core_conf-raw_conf})")
print(f"DNA-confirmed (SAFE core w/ {0.70:.0%} floor):       {safe_conf}/{n}   (+{safe_conf-raw_conf})")
print(f"\nNAIVE recoveries (+{len(core_rec)}) — note the HIGH tail_fracs = domain/repeat-sharer FALSE edges:")
for fid, xc, ncp, tf in core_rec[:15]:
    print(f"  {fid} xc={xc} n={ncp} tail_fracs={tf}")
print(f"\nSAFE recoveries (+{len(safe_rec)}) — genuine: copies mostly align, only a small spurious tail:")
for fid, xc, ncp, tf in safe_rec[:15]:
    print(f"  {fid} xc={xc} n={ncp} tail_fracs={tf}")
if not safe_rec:
    print("  (none — read-support trim already handled inflation; residual failures are real divergence)")
