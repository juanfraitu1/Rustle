#!/usr/bin/env python3
"""Promote genome-wide hidden-copy flags to reference-absent copies, BATCHED for scale:
 per flag -> spectrum filter (drop RNA-editing) -> POA consensus; then ONE minimap2 of all
 consensuses vs the genome (divergence = the absent signal) -> keep divergent (genome_id<0.97,
 maps to own locus, cov>=0.7, ORF>=80) -> blastx only the survivors (name + endogenous check)."""
import sys, json, os, subprocess, time
sys.path.insert(0, "bench")
import copy_assign as CA
import pysam, pyabpoa
from hidden_copy_scan import alts_in_span, MIN_DEPTH, BAL_LO, BAL_HI, SHARE_HI

MM2 = "/home/juanfra/miniforge3/bin/minimap2"
BL = "/home/juanfra/miniforge3/envs/blast/bin"
FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
CAT = "/home/juanfra/winloci_scratch/refabsent"
OUT = f"{CAT}/gw_promoted"; os.makedirs(OUT, exist_ok=True)
COMP = str.maketrans("ACGTN", "TGCAN")
FA = pysam.FastaFile(FASTA)

aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
from itertools import product
CODON = {"".join(c): aas[i] for i, c in enumerate(product("TCAG", repeat=3))}
def orf_aa(seq):
    best = 0
    for strand in (seq, seq.translate(COMP)[::-1]):
        for f in range(3):
            prot = "".join(CODON.get(strand[i:i+3], "X") for i in range(f, len(strand)-2, 3))
            for p in prot.split("*"):
                k = p.find("M")
                if k >= 0: best = max(best, len(p)-k)
    return best


def alt_reads_and_spectrum(bam, chrom, s, e):
    recs = []
    for a in bam.fetch(chrom, s, e):
        if a.is_secondary or a.is_supplementary or a.is_unmapped or a.query_sequence is None: continue
        rpos = a.reference_start; qpos = 0; q = a.query_sequence; alts = []; altb = {}
        for op, ln in (a.cigartuples or []):
            if op in (7, 0): rpos += ln; qpos += ln
            elif op == 8:
                for i in range(ln):
                    if s <= rpos+i < e and qpos+i < len(q):
                        alts.append(rpos+i); altb[rpos+i] = q[qpos+i]
                rpos += ln; qpos += ln
            elif op in (2, 3): rpos += ln
            elif op in (1, 4): qpos += ln
        recs.append((a.reference_start, a.reference_end, alts, altb, q))
        if len(recs) >= 800: break
    if len(recs) < MIN_DEPTH: return None
    from collections import defaultdict, Counter
    ac = defaultdict(list)
    for st, en, alts, altb, _ in recs:
        for p in alts: ac[p].append(altb[p])
    cand = []
    for p, bs in ac.items():
        if len(bs) < 2: continue
        cov = sum(1 for st, en, _, _, _ in recs if st <= p < en)
        if cov >= MIN_DEPTH and BAL_LO <= len(bs)/cov <= BAL_HI: cand.append(p)
    if not cand: return None
    cs = set(cand)
    # spectrum (RNA-editing filter)
    spec = Counter()
    for p in cand:
        refb = FA.fetch(chrom, p, p+1).upper()
        ab = Counter(ac[p]).most_common(1)[0][0]
        if refb in "ACGT" and ab in "ACGT" and refb != ab: spec[(refb, ab)] += 1
    tot = sum(spec.values())
    agtc = (spec[("A","G")] + spec[("T","C")]) / tot if tot else 0
    if agtc >= 0.50 or len(spec) < 5: return ("editing", None)   # RNA editing or too few subtypes
    # alt-haplotype reads
    seqs = []
    for st, en, alts, altb, q in recs:
        cov = sum(1 for p in cand if st <= p < en)
        if cov and sum(1 for p in alts if p in cs)/cov >= SHARE_HI: seqs.append(q)
    if len(seqs) < 3: return None
    return ("ok", seqs)


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else f"{CAT}/genomewide_flags_new.json"
    flags = json.load(open(src))
    print(f"promoting {len(flags)} flags from {os.path.basename(src)}")
    bam = pysam.AlignmentFile(CA.GGO_BAM, "rb")
    aligner = pyabpoa.msa_aligner()
    consfa = open(f"{OUT}/cons.fa", "w")
    meta_by_cid = {}
    t0 = time.time(); n_edit = 0; n_built = 0
    for i, f in enumerate(flags):
        if i % 200 == 0: print(f"  consensus {i}/{len(flags)} ({time.time()-t0:.0f}s) built={n_built} editing={n_edit}", flush=True)
        r = alt_reads_and_spectrum(bam, f["chrom"], f["start"], f["end"])
        if r is None: continue
        if r[0] == "editing": n_edit += 1; continue
        seqs = r[1]
        try:
            cons = aligner.msa(seqs[:60], out_cons=True, out_msa=False).cons_seq[0]
        except Exception:
            continue
        cid = f"{f['chrom']}_{f['start']}"
        meta_by_cid[cid] = dict(chrom=f["chrom"], start=f["start"], end=f["end"],
                                alt_reads=len(seqs), alt_cols=f["n_alt_positions"], cons_len=len(cons), orf=orf_aa(cons))
        consfa.write(f">{cid}\n{cons}\n"); n_built += 1
    consfa.close()
    print(f"built {n_built} consensuses ({n_edit} dropped as RNA-editing). minimap2 vs genome...")
    paf = subprocess.run([MM2, "-cx", "splice:hq", "--eqx", "-N1", "-t", "6", FASTA, f"{OUT}/cons.fa"],
                         capture_output=True, text=True, timeout=1800).stdout
    best = {}
    for line in paf.splitlines():
        c = line.split("\t")
        if len(c) < 12: continue
        cid, qlen, nm, al = c[0], int(c[1]), int(c[9]), int(c[10])
        ident = nm/al if al else 0; cov = al/qlen if qlen else 0; sc = ident*cov
        if cid not in best or sc > best[cid][4]: best[cid] = (c[5], int(c[7]), ident, cov, sc)
    # keep divergent + maps to own locus + coherent ORF
    survivors = []
    for cid, m in meta_by_cid.items():
        g = best.get(cid)
        if not g: continue
        gtgt, gpos, gid, gcov = g[0], g[1], g[2], g[3]
        own = (gtgt == m["chrom"] and abs(gpos - m["start"]) < 200000)
        if gid < 0.97 and gcov >= 0.7 and own and m["orf"] >= 80:
            survivors.append(dict(cid=cid, **m, genome_hit=f"{gtgt}:{gpos}", genome_id=round(gid, 3),
                                  genome_cov=round(gcov, 2), divergence=round(100*(1-gid), 1)))
    survivors.sort(key=lambda r: r["genome_id"])
    print(f"\n{len(survivors)} survive divergence filter (genome_id<0.97, own-locus, ORF>=80). blastx...")
    # write survivor consensuses by pulling from cons.fa, then blastx
    keep = {r["cid"] for r in survivors}
    cur = None
    with open(f"{OUT}/survivors.fa", "w") as fh:
        for line in open(f"{OUT}/cons.fa"):
            if line.startswith(">"): cur = line[1:].strip()
            if cur in keep: fh.write(line)
    bx = subprocess.run([BL+"/blastx", "-query", f"{OUT}/survivors.fa", "-db", f"{CAT}/promoted/proteindb",
                         "-max_target_seqs", "1", "-evalue", "1e-5", "-num_threads", "6",
                         "-outfmt", "6 qseqid sseqid pident length evalue"], capture_output=True, text=True).stdout
    bh = {}
    for line in bx.splitlines():
        q, sid, pid, ln, ev = line.split("\t")
        if q not in bh: bh[q] = (sid, float(pid), int(ln), ev)
    for r in survivors:
        h = bh.get(r["cid"])
        r["protein"] = h[0] if h else "NO-HIT"; r["prot_id"] = h[1] if h else None
    json.dump(survivors, open(f"{OUT}/gw_reference_absent_copies.json", "w"), indent=1)
    print(f"\n{'locus':24s} {'div%':>5} {'cov':>5} {'ORFaa':>5} {'altrd':>5} {'protein':16s} {'pid%':>5}")
    for r in survivors:
        print(f"{r['cid'][:24]:24s} {r['divergence']:>5.1f} {r['genome_cov']:>5.2f} {r['orf']:>5} "
              f"{r['alt_reads']:>5} {str(r['protein'])[:16]:16s} {(r['prot_id'] or 0):>5.1f}")
    print(f"\nwrote {OUT}/gw_reference_absent_copies.json ({len(survivors)} copies)")


if __name__ == "__main__":
    main()
