#!/usr/bin/env python3
"""Milestone 3: which >= 3-copy gorilla families give AS-tied reads among their PRESENT copies? Per family (intervals.tsv, refabsent's 34
excluded): every copy's dominant fibroblast isoform on both haplotypes (o3sim_one), 10 simulated reads per allele, the mini reference
with nothing erased, and the fraction of reads in >= 2 copies' piles (within AS_TIE). Measures ties only, no O2 result.
    o2tiescan.py OUT_DIR [max families]  -> OUT_DIR/tiescan.tsv"""
import collections, json, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import o3sim_one as S
import pysam
out = sys.argv[1]; lim = int(sys.argv[2]) if len(sys.argv) > 2 else 999; skip = int(sys.argv[3]) if len(sys.argv) > 3 else 0
os.makedirs(out, exist_ok=True)
rows = S.tsv(S.SRC); by = collections.defaultdict(list)
for r in rows: by[r["fam"]].append(r)
bonly = {r["family"] for r in S.tsv(f"{S.RA}/refabsent/bonly.tsv")}
cm = {r["pri"]: r for r in S.tsv(f"{S.RA}/chrmap.tsv")}
pri, bam = pysam.FastaFile(S.PRI), pysam.AlignmentFile(S.FIBRO)
fams = sorted((f for f, v in by.items() if len(v) >= 3 and f not in bonly), key=lambda f: (len(by[f]), f))[skip:skip + lim]
new = not os.path.exists(f"{out}/tiescan.tsv")
with open(f"{out}/tiescan.tsv", "a") as o:
    if new:
        o.write("family\tn_copies\tn_reads\ttied\ttied_frac\n")
    for f in fams:
        copies, tx = {}, {}
        for r in sorted(by[f], key=lambda r: (r["chrom"], int(r["clean_start"]))):
            iso = S.dominant_isoform(bam, r["chrom"], int(r["clean_start"]), int(r["clean_end"]))
            if not iso: continue
            bl, st, k, n = iso
            t = S.transcripts(pri, cm, r["chrom"], bl, st)
            copies[r["gene"]] = dict(chrom=r["chrom"], start=int(r["clean_start"]), end=int(r["clean_end"]))
            tx[f"{r['gene']}|A"] = t["A"]
            if t["B"]: tx[f"{r['gene']}|B"] = t["B"]
        if len(copies) < 3: continue
        D = f"{out}/{f}"; os.makedirs(D, exist_ok=True)
        where = S.mini_reference(copies, None, False, D)
        S.simulate(tx, [(t, 10) for t in tx], random.Random(1), f"{D}/reads.fa")
        S.align(D, f"{D}/reads.fa", f"{D}/reads.bam")
        tied, _, _ = S.landing(f"{D}/reads.bam", where, list(copies))
        n = len(tied); t = sum(1 for v in tied.values() if len(v) >= 2)
        o.write(f"{f}\t{len(copies)}\t{n}\t{t}\t{t/max(1,n):.3f}\n"); o.flush()
        print(f, len(copies), n, t, flush=True)
