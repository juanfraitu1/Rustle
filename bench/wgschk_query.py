"""WGSCHK query set (LEDGER_O3): canonical 21-mers of every 1416 unit (units.fa), sorted unique -> qc.u64; per-unit k-mer index
-> unit_kmers.npz. usage: wgschk_query.py UNITS_FA OUTDIR"""
import sys, numpy as np
sys.path.insert(0, "/mnt/linuxdisk/home/juanfraitu/dna_cn")
import kmerlib as kl
names, seqs, cur = [], [], []
for line in open(sys.argv[1]):
    if line[0] == ">":
        if cur: seqs.append("".join(cur))
        names.append(line[1:].strip()); cur = []
    else: cur.append(line.strip())
seqs.append("".join(cur))
ks, off = [], [0]
for s in seqs:
    c = kl.canon(kl.encode(kl.CODE[np.frombuffer(s.upper().encode(), np.uint8)]))
    ks.append(c); off.append(off[-1] + c.size)
allk = np.concatenate(ks); qc = np.unique(allk)
qc.tofile(f"{sys.argv[2]}/qc.u64")
np.savez(f"{sys.argv[2]}/unit_kmers.npz", idx=np.searchsorted(qc, allk), off=np.array(off), names=np.array(names))
print("units", len(names), "kmers", allk.size, "unique", qc.size)
