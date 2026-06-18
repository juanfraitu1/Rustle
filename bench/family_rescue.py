#!/usr/bin/env python3
"""Family-aware copy RESCUE (borrow strength / partial pooling): recover under-ASSEMBLED copies that fall
below the >= 3-read general-assembly gate but are sequence-homologous to a CONFIRMED de-novo family.

Rationale: a single confident (canonical-junction) read forming a multi-exon chain that POA-confirms against
an existing family is strong evidence of a real copy -- the family PRIOR overcomes thin read support. This
recovers expression-limited copies (e.g. RFPL2) WITHOUT lowering the genome-wide gate (fuzzy/low-gate
experiments showed that only adds noise).

FAST: scan merged family neighbourhoods ONCE; a CANONICAL-k-mer pre-filter (a k-mer == its reverse-
complement) keeps only thin loci sharing k-mers with the family; POA-confirm (parallel, both orientations).
Run with /home/juanfra/miniforge3/bin/python.
"""
import multiprocessing as mp
import os
import sys
from collections import defaultdict

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(__file__))
import poa_family_definition as P
from denovo_families import kmer_hashes, KMER   # canonical k-mer hasher (strand-symmetric)

HERE = os.path.dirname(__file__)
SPLIT = os.path.join(HERE, "denovo_families_split.tsv")
GGO_FA = "/home/juanfra/winloci_scratch/GGO.fasta"
GGO_BAM = "/home/juanfra/winloci_scratch/GGO.bam"
SPLICED_FA = "/home/juanfra/winloci_scratch/denovo_transcripts.fa"
META_F = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
OUT = os.path.join(HERE, "denovo_rescued_copies.tsv")

WIN = 1_000_000          # co-located rescue radius (RFPL2 sits ~170kb from RFPL3); bounds the scan
T_CORE = 0.13
LEN_CAP = 9000
MIN_SUPPORT = 1
K_RESCUE = 20            # thin locus must share >= this many canonical k-mers with the family (real copies
                         #   share many; common-domain sharers share few -> bounds POA to true candidates)
PLUS = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
MINUS = {("CT", "AC"), ("CT", "GC"), ("GT", "AT")}
_RC = str.maketrans("ACGTN", "TGCAN")

_MEMSEQ = None


def _init(memseq):
    global _MEMSEQ, _A
    _MEMSEQ = memseq
    _A = P.make_aligner()


def _poa_rescue(task):
    """POA the thin locus against ONLY its single best k-mer-matching member (both orientations)."""
    idx, seq, tid, fid = task
    ms = _MEMSEQ.get(tid)
    if not ms or min(len(seq), len(ms)) > LEN_CAP:
        return None
    cr = P.poa_pair_stats(_A, seq, ms)["core_recip"]
    if cr < T_CORE:
        cr = max(cr, P.poa_pair_stats(_A, seq, ms.translate(_RC)[::-1])["core_recip"])
    return (idx, fid, round(cr, 3)) if cr >= T_CORE else None


def main():
    meta = {}
    for line in open(META_F):
        if line.startswith("id"):
            continue
        tid, c, s, e, st, ne, nr = line.rstrip("\n").split("\t")
        meta[tid] = (c, int(s), int(e), st, int(nr))
    spliced = pysam.FastaFile(SPLICED_FA)
    sref = set(spliced.references)
    fa = pysam.FastaFile(GGO_FA)
    bam = pysam.AlignmentFile(GGO_BAM, "rb")

    # discrete families -> members; member seqs, spans, canonical-k-mer sets
    fams = {}
    for line in open(SPLIT):
        if line.startswith("family_id"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[8] != "family":
            continue
        mem = [t for t in f[-1].split(",") if t in meta and t in sref]
        if len(mem) >= 2:
            fams[f[0]] = mem
    memseq = {t: spliced.fetch(t) for mem in fams.values() for t in mem}
    memkmers = {t: set(kmer_hashes(s).tolist()) if s else set() for t, s in memseq.items()}
    tid2fam = {t: fid for fid, mem in fams.items() for t in mem}

    # member span index (per chrom, with tid) + merged neighbourhood intervals
    members = defaultdict(list)      # chrom -> [(start,end,tid)]
    member_spans = defaultdict(list)
    for fid, mem in fams.items():
        for t in mem:
            c, s, e, st, nr = meta[t]
            members[c].append((s, e, t))
            member_spans[c].append((s, e))
    intervals = defaultdict(list)
    for c, lst in members.items():
        ivs = sorted((max(0, s - WIN), e + WIN) for s, e, _ in lst)
        cur = list(ivs[0])
        for s, e in ivs[1:]:
            if s <= cur[1]:
                cur[1] = max(cur[1], e)
            else:
                intervals[c].append(tuple(cur)); cur = [s, e]
        intervals[c].append(tuple(cur))
    for c in members:
        members[c].sort()

    # scan merged intervals once -> thin loci; k-mer pre-filter against nearby MEMBERS
    cand_tasks = []       # (idx, seq, best_member_tid, fid)
    cand_meta = []        # idx -> (chrom,start,end,strand,nexon,support)
    n_thin = 0
    for c, ivs in intervals.items():
        mstarts = [m[0] for m in members[c]]
        for lo, hi in ivs:
            chains = defaultdict(lambda: [0, 10**12, 0])
            for aln in bam.fetch(c, lo, hi):
                if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                    continue
                rp = aln.reference_start
                intr = []
                for op, length in aln.cigartuples:
                    if op == 3:
                        intr.append((rp, rp + length)); rp += length
                    elif op in (0, 2, 7, 8):
                        rp += length
                if not intr:
                    continue
                g = chains[tuple(intr)]
                g[0] += 1; g[1] = min(g[1], aln.reference_start); g[2] = max(g[2], aln.reference_end)
            # collapse to loci
            items = sorted(chains.items(), key=lambda kv: kv[1][1])
            loci = []
            for intr, (sup, s, e) in items:
                if sup < MIN_SUPPORT:
                    continue
                for L in loci:
                    if s <= L["e"] and e >= L["s"]:
                        L["s"] = min(L["s"], s); L["e"] = max(L["e"], e)
                        if (sup, e - s) > (L["sup"], L["e2"] - L["s2"]):
                            L.update(sup=sup, intr=intr, s2=s, e2=e)
                        break
                else:
                    loci.append(dict(s=s, e=e, sup=sup, intr=intr, s2=s, e2=e))
            for L in loci:
                s, e, intr, sup = L["s2"], L["e2"], L["intr"], L["sup"]
                if any(s < me and e > ms for ms, me in member_spans[c]):
                    continue                                       # already an assembled member
                n_thin += 1
                exons = []; prev = s
                for d, a in intr:
                    exons.append((prev, d)); prev = a
                exons.append((prev, e))
                strand = None; ok = True
                for d, a in intr:
                    don = fa.fetch(c, d, d + 2).upper(); acc = fa.fetch(c, a - 2, a).upper()
                    st = "+" if (don, acc) in PLUS else ("-" if (don, acc) in MINUS else None)
                    if st is None or (strand and st != strand):
                        ok = False; break
                    strand = st
                if not ok:
                    continue
                seq = "".join(fa.fetch(c, xs, xe).upper() for xs, xe in exons)
                if len(seq) < 200 or len(seq) > LEN_CAP:
                    continue
                if strand == "-":
                    seq = seq.translate(_RC)[::-1]
                kset = set(kmer_hashes(seq).tolist())
                if not kset:
                    continue
                # best-matching MEMBER within WIN by canonical-k-mer overlap (>= K_RESCUE)
                best_tid = None; best_ov = K_RESCUE - 1
                for ms_, me_, t in members[c]:
                    if ms_ - WIN <= s <= me_ + WIN:
                        ov = len(kset & memkmers[t])
                        if ov > best_ov:
                            best_ov, best_tid = ov, t
                if best_tid is None:
                    continue
                idx = len(cand_meta)
                cand_tasks.append((idx, seq, best_tid, tid2fam[best_tid]))
                cand_meta.append((c, s, e, strand or "+", len(intr) + 1, sup, tuple(intr), seq))

    print(f"families: {len(fams)}; thin loci scanned: {n_thin:,}; k-mer-prefiltered POA candidates: {len(cand_tasks):,}", flush=True)
    with mp.Pool(5, initializer=_init, initargs=(memseq,)) as pool:
        results = [r for r in pool.map(_poa_rescue, cand_tasks, chunksize=32) if r]

    # dedup by locus; keep best core_recip
    bykey = {}
    for idx, fid, cr in results:
        c, s, e, strand, nex, sup, intr, seq = cand_meta[idx]
        k = (c, s, e)
        if k not in bykey or cr > bykey[k]["cr"]:
            bykey[k] = dict(fid=fid, tid=f"RC_{c}_{s}_{nex}", c=c, s=s, e=e, strand=strand,
                            nex=nex, sup=sup, cr=cr, intr=intr, seq=seq)
    rescued = sorted(bykey.values(), key=lambda r: -r["cr"])
    with open(OUT, "w") as fh:
        fh.write("family_id\ttid\tchrom\tstart\tend\tstrand\tn_exon\tsupport\tcore_recip\tintrons\n")
        for r in rescued:
            intr_s = ";".join(f"{a}-{b}" for a, b in r["intr"])
            fh.write(f"{r['fid']}\t{r['tid']}\t{r['c']}\t{r['s']}\t{r['e']}\t{r['strand']}\t"
                     f"{r['nex']}\t{r['sup']}\t{r['cr']}\t{intr_s}\n")
    with open(OUT.replace(".tsv", ".fa"), "w") as fh:    # spliced seqs so rescued copies are assignable
        for r in rescued:
            fh.write(f">{r['tid']}\n")
            for i in range(0, len(r["seq"]), 80):
                fh.write(r["seq"][i:i + 80] + "\n")
    sup_dist = defaultdict(int)
    for r in rescued:
        sup_dist[r["sup"]] += 1
    print(f"RESCUED copies (POA core_recip >= {T_CORE}): {len(rescued):,}  by support: {dict(sorted(sup_dist.items()))}")
    print(f"[wrote {OUT} + .fa]")
    print("\nRFPL region rescued copies:")
    for r in rescued:
        if r["c"] == "NC_086018.1" and 27_000_000 <= r["s"] <= 31_000_000:
            print(f"  {r['fid']} {r['c']}:{r['s']}-{r['e']} {r['strand']} exons={r['nex']} support={r['sup']} core_recip={r['cr']}")


if __name__ == "__main__":
    main()
