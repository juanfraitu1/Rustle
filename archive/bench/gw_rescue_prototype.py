#!/usr/bin/env python3
"""Genome-wide family-aware rescue prototype for gw_xcbase catalog.

Scans thin loci (primary-read intron chains with support 1-2) in a +/-1 Mb
neighbourhood around existing gw_xcbase families, builds their spliced sequences,
and POA-confirms them against family members. Confirmed thin loci are emitted as
rescued copies.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/gw_rescue_prototype.py
"""
import csv
import os
import sys
from collections import defaultdict

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import multiprocessing as mp
import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(__file__))
import poa_family_definition as P
from denovo_families import kmer_hashes, KMER  # canonical k-mer hasher

SCRATCH = "/home/juanfra/winloci_scratch"
HERE = os.path.dirname(__file__)
GGO_FA = os.path.join(SCRATCH, "GGO.fasta")
GGO_BAM = os.path.join(SCRATCH, "GGO.bam")
CATALOG_PREFIX = os.path.join(SCRATCH, "gw_xcbase")
OUT_RESCUED = os.path.join(HERE, "gw_rescue_prototype.rescued.tsv")
# actual files are <prefix>.copies.tsv etc.
OUT_CATALOG_COPIES = os.path.join(SCRATCH, "gw_rescue.copies.tsv")
OUT_CATALOG_FAMILIES = os.path.join(SCRATCH, "gw_rescue.families.tsv")

WIN = 1_000_000
T_CORE = 0.13
LEN_CAP = 9000
MIN_SUPPORT = 1
K_RESCUE = 20
MIN_LEN = 200
PLUS = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
MINUS = {("CT", "AC"), ("CT", "GC"), ("GT", "AT")}
_RC = str.maketrans("ACGTN", "TGCAN")

_MEMSEQ = None


def _init(memseq):
    global _MEMSEQ, _A
    _MEMSEQ = memseq
    _A = P.make_aligner()


def _poa_rescue(task):
    idx, seq, tid, fid = task
    ms = _MEMSEQ.get(tid)
    if not ms or min(len(seq), len(ms)) > LEN_CAP:
        return None
    cr = P.poa_pair_stats(_A, seq, ms)["core_recip"]
    if cr < T_CORE:
        cr = max(cr, P.poa_pair_stats(_A, seq, ms.translate(_RC)[::-1])["core_recip"])
    return (idx, fid, round(cr, 3)) if cr >= T_CORE else None


def parse_copies_fa(path):
    """Parse gw_xcbase.copies.fa -> dict (family_id, copy_idx) -> seq.

    FASTA header format: >FAM|copy_idx|chrom:start-end|strand|nexon=N
    """
    seqs = {}
    cur = None
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                parts = ln[1:].split()[0].split("|")
                cur = (parts[0], int(parts[1]))
                seqs[cur] = []
            elif cur is not None:
                seqs[cur].append(ln.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def load_catalog(copies_tsv, families_tsv):
    copies = []
    with open(copies_tsv) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            copies.append({
                "family_id": row["family_id"],
                "copy_idx": int(row["copy_idx"]),
                "tid": row["tid"],
                "chrom": row["chrom"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "n_exon": int(row["n_exon"]),
                "strand": row["strand"],
                "n_reads": int(row["n_reads"]),
            })
    fams = defaultdict(list)
    for c in copies:
        fams[c["family_id"]].append(c)
    return copies, fams


def build_spliced_seq(fa, chrom, introns, start, end, strand):
    """Build spliced transcript sequence from intron chain and strand."""
    exons = []
    prev = start
    for d, a in introns:
        exons.append((prev, d))
        prev = a
    exons.append((prev, end))
    seq = "".join(fa.fetch(chrom, xs, xe).upper() for xs, xe in exons)
    if strand == "-":
        seq = seq.translate(_RC)[::-1]
    return seq


def canonical_kmer_set(seq):
    return set(kmer_hashes(seq).tolist()) if seq else set()


def thin_loci_from_bam(bam, chrom, lo, hi, min_support=1):
    """Collect intron-chain groups with support < 3 in [lo,hi).

    Collapse overlapping chains ONLY when they share at least one junction (real
    isoforms of the same locus). This avoids merging read-throughs that span two
    distinct loci, which would cause the cluster to be discarded because it
    overlaps an existing family member.
    """
    chains = defaultdict(lambda: [0, 10**12, 0])
    for aln in bam.fetch(chrom, lo, hi):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        rp = aln.reference_start
        intr = []
        for op, length in aln.cigartuples:
            if op == 3:
                intr.append((rp, rp + length))
                rp += length
            elif op in (0, 2, 7, 8):
                rp += length
        if not intr:
            continue
        key = tuple(intr)
        chains[key][0] += 1
        chains[key][1] = min(chains[key][1], aln.reference_start)
        chains[key][2] = max(chains[key][2], aln.reference_end)
    items = sorted(chains.items(), key=lambda kv: kv[1][1])
    loci = []
    for intr, (sup, s, e) in items:
        if sup < min_support or sup >= 3:
            continue
        intr_set = set(intr)
        merged = False
        for L in loci:
            if s <= L["e"] and e >= L["s"] and intr_set.intersection(L["intr_set"]):
                L["s"] = min(L["s"], s)
                L["e"] = max(L["e"], e)
                L["intr_set"].update(intr_set)
                if (sup, e - s) > (L["sup"], L["e2"] - L["s2"]):
                    L.update(sup=sup, intr=intr, s2=s, e2=e)
                merged = True
                break
        if not merged:
            loci.append(dict(s=s, e=e, sup=sup, intr=intr, s2=s, e2=e, intr_set=intr_set))
    for L in loci:
        L.pop("intr_set", None)
    return loci


def main():
    copies, fams = load_catalog(
        os.path.join(SCRATCH, "gw_xcbase.copies.tsv"),
        os.path.join(SCRATCH, "gw_xcbase.families.tsv"),
    )
    seqs = parse_copies_fa(os.path.join(SCRATCH, "gw_xcbase.copies.fa"))
    print(f"Loaded {len(copies)} copies in {len(fams)} families; {len(seqs)} sequences")

    fa = pysam.FastaFile(GGO_FA)
    bam = pysam.AlignmentFile(GGO_BAM, "rb")

    # member span index + merged neighbourhood intervals per chrom
    members_by_chrom = defaultdict(list)  # chrom -> [(s,e,tid,fid)]
    member_spans_by_chrom = defaultdict(list)
    for c in copies:
        members_by_chrom[c["chrom"]].append((c["start"], c["end"], c["tid"], c["family_id"]))
        member_spans_by_chrom[c["chrom"]].append((c["start"], c["end"]))

    intervals = defaultdict(list)
    for c, lst in members_by_chrom.items():
        ivs = sorted((max(0, s - WIN), e + WIN) for s, e, _, _ in lst)
        cur = list(ivs[0])
        for s, e in ivs[1:]:
            if s <= cur[1]:
                cur[1] = max(cur[1], e)
            else:
                intervals[c].append(tuple(cur))
                cur = [s, e]
        intervals[c].append(tuple(cur))

    # Precompute member sequences and k-mer sets
    memseq = {}
    memkmers = {}
    tid2fam = {}
    for c in copies:
        tid = c["tid"]
        key = (c["family_id"], c["copy_idx"])
        if key not in seqs:
            continue
        memseq[tid] = seqs[key]
        memkmers[tid] = canonical_kmer_set(seqs[key])
        tid2fam[tid] = c["family_id"]

    # flat member list for global k-mer prefilter (cross-chrom rescue)
    all_members = [(c["chrom"], c["start"], c["end"], c["tid"], c["family_id"]) for c in copies]

    cand_tasks = []
    cand_meta = []
    n_thin = 0
    total_ivs = sum(len(v) for v in intervals.values())
    print(f"Scanning {total_ivs} intervals across {len(intervals)} chroms", flush=True)
    for chrom, ivs in intervals.items():
        for ivi, (lo, hi) in enumerate(ivs):
            print(f"  [{chrom}] interval {ivi+1}/{len(ivs)} {lo}-{hi} ({hi-lo:,} bp)", flush=True)
            loci = thin_loci_from_bam(bam, chrom, lo, hi, MIN_SUPPORT)
            print(f"    -> {len(loci)} thin loci", flush=True)
            for L in loci:
                s, e, intr, sup = L["s2"], L["e2"], L["intr"], L["sup"]
                # Allow thin loci that overlap an existing member ONLY if the
                # overlap covers < 50% of the thin locus (i.e. the locus extends
                # well beyond the existing member). This rescues read-throughs
                # that span a missed copy adjacent to an existing fragment.
                span_len = e - s
                overlapping_spans = [(ms, me) for ms, me in member_spans_by_chrom[chrom] if s < me and e > ms]
                if overlapping_spans:
                    max_overlap_frac = max(
                        (min(e, me) - max(s, ms)) / span_len for ms, me in overlapping_spans
                    )
                    if max_overlap_frac >= 0.5:
                        continue
                n_thin += 1
                # canonical-junction strand inference
                strand = None
                ok = True
                for d, a in intr:
                    don = fa.fetch(chrom, d, d + 2).upper()
                    acc = fa.fetch(chrom, a - 2, a).upper()
                    st = "+" if (don, acc) in PLUS else ("-" if (don, acc) in MINUS else None)
                    if st is None or (strand and st != strand):
                        ok = False
                        break
                    strand = st
                if not ok:
                    continue
                seq = build_spliced_seq(fa, chrom, intr, s, e, strand)
                if len(seq) < MIN_LEN or len(seq) > LEN_CAP:
                    continue
                kset = canonical_kmer_set(seq)
                if not kset:
                    continue
                # best member GLOBALLY by canonical k-mer overlap (cross-chrom rescue)
                best_tid = None
                best_ov = K_RESCUE - 1
                for mc, ms_, me_, t, _ in all_members:
                    # same-chrom members use WIN bound; cross-chrom require much
                    # stronger k-mer evidence to avoid spurious distant matches
                    same_chrom = mc == chrom
                    if same_chrom and not (ms_ - WIN <= s <= me_ + WIN):
                        continue
                    ov = len(kset & memkmers.get(t, set()))
                    if not same_chrom and ov < 100:
                        continue
                    if ov > best_ov:
                        best_ov, best_tid = ov, t
                if best_tid is None:
                    continue
                idx = len(cand_meta)
                cand_tasks.append((idx, seq, best_tid, tid2fam[best_tid]))
                cand_meta.append((chrom, s, e, strand or "+", len(intr) + 1, sup, tuple(intr), seq))

    print(f"families: {len(fams)}; thin loci scanned: {n_thin:,}; k-mer-prefiltered POA candidates: {len(cand_tasks):,}", flush=True)

    with mp.Pool(5, initializer=_init, initargs=(memseq,)) as pool:
        results = [r for r in pool.map(_poa_rescue, cand_tasks, chunksize=32) if r]

    bykey = {}
    for idx, fid, cr in results:
        chrom, s, e, strand, nex, sup, intr, seq = cand_meta[idx]
        k = (chrom, s, e)
        if k not in bykey or cr > bykey[k]["cr"]:
            bykey[k] = dict(
                fid=fid, tid=f"RC_{chrom}_{s}_{nex}", chrom=chrom, s=s, e=e,
                strand=strand, nex=nex, sup=sup, cr=cr, intr=intr, seq=seq
            )
    rescued = sorted(bykey.values(), key=lambda r: -r["cr"])
    print(f"Rescued {len(rescued)} copies")

    with open(OUT_RESCUED, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["family_id", "tid", "chrom", "start", "end", "n_exon", "strand", "n_reads", "core_recip"])
        for r in rescued:
            w.writerow([r["fid"], r["tid"], r["chrom"], r["s"], r["e"], r["nex"], r["strand"], r["sup"], r["cr"]])
    print(f"Wrote {OUT_RESCUED}")

    # Build merged catalog
    with open(OUT_CATALOG_FAMILIES, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["family_id", "n_copies", "n_chroms", "chroms", "cross_chrom", "avg_reads"])
    with open(OUT_CATALOG_COPIES, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["family_id", "copy_idx", "tid", "chrom", "start", "end", "n_exon", "strand", "n_reads"])

    rescued_by_fam = defaultdict(list)
    for r in rescued:
        rescued_by_fam[r["fid"]].append(r)

    with open(OUT_CATALOG_FAMILIES, "a", newline="") as ff, open(OUT_CATALOG_COPIES, "a", newline="") as fc:
        wf = csv.writer(ff, delimiter="\t")
        wc = csv.writer(fc, delimiter="\t")
        for fid in sorted(fams.keys()):
            fam_copies = list(fams[fid])
            for r in rescued_by_fam.get(fid, []):
                fam_copies.append({
                    "family_id": fid,
                    "copy_idx": len(fam_copies),
                    "tid": r["tid"],
                    "chrom": r["chrom"],
                    "start": r["s"],
                    "end": r["e"],
                    "n_exon": r["nex"],
                    "strand": r["strand"],
                    "n_reads": r["sup"],
                })
            chroms = sorted(set(c["chrom"] for c in fam_copies))
            cross = len(chroms) > 1
            avg_reads = sum(c["n_reads"] for c in fam_copies) / max(1, len(fam_copies))
            wf.writerow([fid, len(fam_copies), len(chroms), ",".join(chroms), cross, f"{avg_reads:.1f}"])
            fam_copies.sort(key=lambda c: (c["chrom"], c["start"]))
            for ci, c in enumerate(fam_copies):
                wc.writerow([fid, ci, c["tid"], c["chrom"], c["start"], c["end"], c["n_exon"], c["strand"], c["n_reads"]])
    print(f"Wrote merged catalog: {OUT_CATALOG_FAMILIES}, {OUT_CATALOG_COPIES}")

    fa.close()
    bam.close()


if __name__ == "__main__":
    main()
