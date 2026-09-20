#!/usr/bin/env python
"""Build the GENOME K-MER MULTIPLICITY cache for the read-intrinsic TE gate.

Self-contained, NO RepeatMasker library: count how many times each canonical 21-mer of a bridge locus's
SPLICED EXON sequence occurs across the whole GGO.fasta genome. A TE-derived exon is built from k-mers that
recur thousands of times genome-wide (including silent, non-transcribed genomic TE copies a
homology-graph/read signal never sees); a unique-sequence paralog's k-mers are ~single-copy. Mask-AGNOSTIC
(sequence uppercased) so it is INDEPENDENT of the RepeatMasker soft-mask feature -> a genuinely orthogonal
signal, and it catches gorilla-specific / novel TEs a human-biased RepeatMasker library misses.

Method: (1) collect the union of canonical 21-mers over all bridge exonic sequences; EXPAND it to both
orientations (Q_exp = each query canonical k-mer AND its reverse complement) so the genome pass only needs
FORWARD k-mers (21 accumulation passes, not 42 -> ~2x faster). (2) stream the genome contig-by-contig in
25 Mbp windows, compute forward 21-mers with vectorized numpy, searchsorted into Q_exp, accumulate counts
with C-speed np.bincount. (3) per query canonical c, genome multiplicity N(c) = count_fwd(c) +
count_fwd(rc(c)) (a genome position with canon==c has forward k-mer either c or rc(c)). Exact count,
DETERMINISTIC. CHECKPOINTED per contig (bench/ri_kmer_ckpt.npz) -> resumable if interrupted.

Reproduce: PYTHONHASHSEED=0 python bench/ri_build_kmer_cache.py
Writes: bench/ri_kmer_cache.tsv  (dn, n_kmers, km_median, km_mean, km_p90, km_frac_ge_100)
"""
import os, sys
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"; os.execv(sys.executable, [sys.executable] + sys.argv)
BENCH = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, BENCH)
import numpy as np
import family_er_pr as F
import family_er_te_gate as TE
import pysam

SCRATCH = "/home/juanfra/winloci_scratch"
GENOME = os.path.join(SCRATCH, "GGO.fasta")
CACHE = os.path.join(BENCH, "ri_kmer_cache.tsv")
CKPT = os.path.join(BENCH, "ri_kmer_ckpt.npz")
K = 21
WIN = 25_000_000
_TWO = np.uint64(2); _THREE = np.uint64(3)

_LUT = np.full(256, 4, dtype=np.uint8)
for _ch, _c in (("A", 0), ("C", 1), ("G", 2), ("T", 3)):
    _LUT[ord(_ch)] = _c
    _LUT[ord(_ch.lower())] = _c


def _canon_or_fwd(code, canonical):
    """Forward (canonical=False) or canonical (True) k-mer uint64 array over VALID windows."""
    L = len(code)
    if L < K:
        return np.empty(0, dtype=np.uint64)
    nk = L - K + 1
    clean = np.where(code < 4, code, 0).astype(np.uint64)
    fwd = np.zeros(nk, dtype=np.uint64)
    for j in range(K):
        fwd = (fwd << _TWO) | clean[j:j + nk]
    if canonical:
        comp = _THREE - clean
        rc = np.zeros(nk, dtype=np.uint64)
        for j in range(K):
            rc |= (comp[j:j + nk] << np.uint64(2 * (K - 1 - j)))
        out = np.minimum(fwd, rc)
        del comp, rc
    else:
        out = fwd
    bad = (code == 4).astype(np.int64)
    badcum = np.concatenate(([0], np.cumsum(bad)))
    valid = (badcum[K:K + nk] - badcum[0:nk]) == 0
    return out[valid]


def rc_of(kmers):
    rc = np.zeros_like(kmers); x = kmers.copy()
    for _ in range(K):
        rc = (rc << _TWO) | (_THREE - (x & _THREE))
        x = x >> _TWO
    return rc


def main():
    meta = F.load_meta()
    skel = TE.load_skeletons()
    fa = pysam.FastaFile(GENOME)

    dn = set()
    with open(os.path.join(BENCH, "edge_bridge_mask.tsv")) as fh:
        h = fh.readline().rstrip("\n").split("\t"); idx = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[idx["core"]] == "":
                continue
            for c in ("bridge_dn_a", "bridge_dn_b"):
                if f[idx[c]]:
                    dn.add(f[idx[c]])
    loci = sorted(dn)
    print(f"[kmer] {len(loci)} target bridge loci", flush=True)

    # per-locus canonical kmer arrays + union query set (canonical), then expand to both orientations
    dn_km = {}
    parts = []
    for i, d in enumerate(loci):
        chrom, exons = TE.dn_exons(d, meta, skel)
        if chrom is None:
            dn_km[d] = np.empty(0, dtype=np.uint64); continue
        s = "".join(fa.fetch(chrom, a, b) for a, b in exons).upper()
        arr = np.frombuffer(s.encode("ascii"), dtype=np.uint8)
        km = np.unique(_canon_or_fwd(_LUT[arr], canonical=True)) if len(arr) >= K else np.empty(0, np.uint64)
        dn_km[d] = km
        if km.size:
            parts.append(km)
        if (i + 1) % 500 == 0:
            print(f"[kmer] extracted {i+1}/{len(loci)} loci", flush=True)
    Q_can = np.unique(np.concatenate(parts)) if parts else np.empty(0, dtype=np.uint64)
    Q_exp = np.unique(np.concatenate([Q_can, rc_of(Q_can)])) if Q_can.size else Q_can.copy()
    print(f"[kmer] |Q_can|={len(Q_can)} |Q_exp|={len(Q_exp)} (both orientations); streaming genome ...",
          flush=True)

    # multiplicative-hash PREFILTER: reject genome k-mers that cannot be in Q_exp before the (expensive)
    # exact searchsorted -> ~13x fewer exact lookups (searchsorted over 13.5M is the bottleneck). The exact
    # searchsorted on survivors still enforces correctness (prefilter false-positives are rejected).
    GOLD = np.uint64(0x9E3779B97F4A7C15); B = 28; SH = np.uint64(64 - B)
    present = np.zeros(1 << B, dtype=bool)
    if Q_exp.size:
        present[(Q_exp * GOLD) >> SH] = True

    # resume from checkpoint if compatible
    cntE = np.zeros(len(Q_exp), dtype=np.int64)
    start_ci = 0
    if os.path.exists(CKPT):
        z = np.load(CKPT)
        if int(z["qexp_len"]) == len(Q_exp):
            cntE = z["cntE"].astype(np.int64); start_ci = int(z["done"])
            print(f"[kmer] resumed from checkpoint: {start_ci} contigs done", flush=True)

    refs = list(fa.references)
    for ci in range(start_ci, len(refs)):
        contig = refs[ci]
        clen = fa.get_reference_length(contig)
        s = 0
        while s < clen:
            e = min(s + WIN + K - 1, clen)
            seq = fa.fetch(contig, s, e).upper()
            arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
            fwd = _canon_or_fwd(_LUT[arr], canonical=False)
            del arr, seq
            if fwd.size:
                surv = fwd[present[(fwd * GOLD) >> SH]]   # prefilter -> only plausible k-mers
                del fwd
                if surv.size:
                    pos = np.searchsorted(Q_exp, surv)
                    ok = pos < len(Q_exp)
                    pos_ok = pos[ok]
                    mi = pos_ok[Q_exp[pos_ok] == surv[ok]]
                    if mi.size:
                        u, c = np.unique(mi, return_counts=True)
                        cntE[u] += c
                del surv
            else:
                del fwd
            s += WIN
        np.savez(CKPT, cntE=cntE, done=ci + 1, qexp_len=len(Q_exp))
        print(f"[kmer] contig {ci+1}/{len(refs)} {contig} done (len {clen}) [checkpointed]", flush=True)

    # per-locus summary: N(c) = cntE[pos(c)] + cntE[pos(rc(c))]  (rc!=c)
    def mult_for(km):
        if km.size == 0:
            return None
        rc = rc_of(km)
        pc = np.searchsorted(Q_exp, km)
        pr = np.searchsorted(Q_exp, rc)
        base = cntE[pc].astype(np.int64)
        add = np.where(rc != km, cntE[pr].astype(np.int64), 0)
        return base + add

    with open(CACHE, "w") as out:
        out.write("dn\tn_kmers\tkm_median\tkm_mean\tkm_p90\tkm_frac_ge_100\n")
        for d in loci:
            cnts = mult_for(dn_km[d])
            if cnts is None:
                out.write(f"{d}\t0\t\t\t\t\n"); continue
            out.write(f"{d}\t{cnts.size}\t{float(np.median(cnts)):.4f}\t{float(cnts.mean()):.4f}\t"
                      f"{float(np.percentile(cnts,90)):.4f}\t{float((cnts>=100).mean()):.4f}\n")
    print(f"[kmer] wrote {CACHE}", flush=True)


if __name__ == "__main__":
    main()
