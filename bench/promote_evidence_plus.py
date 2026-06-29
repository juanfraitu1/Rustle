#!/usr/bin/env python3
"""Enrich the 18 promotion candidates with the two missing axes:
 (1) GENOME homology: map each consensus to the whole T2T assembly (where does it land, identity,
     coverage) -- resolves 'no-homology-to-family-models' into real-but-different vs garbage, and
     confirms the distinct ones land on their own locus with the private SNPs.
 (2) DEPTH-DOUBLING: per family, primary-read coverage of the flagged copy vs the family's other
     copies -- a true COLLAPSED copy (2 loci -> 1 reference spot) reads ~2x a single sibling; a HET
     (1 locus, 2 alleles) does not. This is the copy-vs-het discriminator that identity cannot give.
Identity bands for het-vs-paralog: >99.5% id (=<0.5% div) ~ het-range; <97% (>3% div) = paralog."""
import sys, json, os, subprocess
sys.path.insert(0, "bench")
import copy_assign as CA
import pysam

MM2 = "/home/juanfra/miniforge3/bin/minimap2"
FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
CAT = "/home/juanfra/winloci_scratch/refabsent"
OUT = f"{CAT}/promoted"


def genome_hits():
    """map all consensus seqs to the genome; return {cid: (target, pos, identity, qcov)}."""
    paf = subprocess.run([MM2, "-cx", "splice:hq", "--eqx", "-N1", FASTA, f"{OUT}/consensus.fa"],
                         capture_output=True, text=True, timeout=600).stdout
    best = {}
    for line in paf.splitlines():
        f = line.split("\t")
        if len(f) < 12: continue
        cid, qlen, nmatch, alnlen = f[0], int(f[1]), int(f[9]), int(f[10])
        ident = nmatch / alnlen if alnlen else 0
        qcov = alnlen / qlen if qlen else 0
        sc = ident * qcov
        if cid not in best or sc > best[cid][4]:
            best[cid] = (f[5], int(f[7]), ident, qcov, sc)
    return best


def copy_depths():
    """per family copy: mean primary-read coverage (primary reads * read-overlap / span). Cheap proxy:
    primary read COUNT over the copy span normalised by span length (reads per kb)."""
    meta = CA.load_meta(); skel = CA.load_skel(); spliced = CA.load_spliced()
    resc = CA.load_rescued(meta, skel, spliced)
    fams = CA.colocated_families(meta, rescued_fam=resc)
    bam = pysam.AlignmentFile(CA.GGO_BAM, "rb")
    fam_cov = {}
    for fid, dom, copies in fams:
        covs = {}
        for (tid, chrom, s1, e, strand) in copies:
            s = s1 - 1
            n = sum(1 for a in bam.fetch(chrom, s, e)
                    if not (a.is_secondary or a.is_supplementary or a.is_unmapped))
            covs[(chrom, s)] = n / max(1, (e - s) / 1000.0)  # primary reads per kb
        fam_cov[fid] = covs
    return fam_cov


def main():
    rows = json.load(open(f"{OUT}/promoted.json"))
    gh = genome_hits()
    fam_cov = copy_depths()
    print(f"{'family':9s} {'start':>9} {'id_fam':>6} {'GENOME_hit':22s} {'g_id':>5} {'g_cov':>5} "
          f"{'depthVSsib':>10} {'call':18s}")
    out = []
    for r in rows:
        cid = f"{r['family']}_{r['chrom']}_{r['start']}"
        g = gh.get(cid)
        gtgt = f"{g[0]}:{g[1]}" if g else "NONE"
        gid = g[2] if g else 0.0
        gcov = g[3] if g else 0.0
        # depth vs siblings
        covs = fam_cov.get(r["family"], {})
        this = covs.get((r["chrom"], r["start"]))
        others = [v for k, v in covs.items() if k != (r["chrom"], r["start"])]
        import statistics
        med_sib = statistics.median(others) if others else None
        depth_ratio = (this / med_sib) if (this and med_sib) else None
        # integrated call (RNA-only; honest about het-vs-paralog limit)
        idf = r["identity"] or 0
        gmaps = gcov >= 0.5
        if not gmaps:
            call = "consensus-artifact"
        elif idf and idf < 0.97:
            call = "PARALOG-COPY"          # too divergent for het
        elif depth_ratio and depth_ratio >= 1.7:
            call = "COLLAPSED-COPY(2x)"     # ~identical but 2x depth = collapsed copy not het
        elif idf >= 0.995:
            call = "het-likely"
        else:
            call = "ambiguous(needs-DNA)"
        out.append(dict(**r, genome_hit=gtgt, genome_id=round(gid, 3), genome_cov=round(gcov, 2),
                        depth_ratio=(round(depth_ratio, 2) if depth_ratio else None), call=call))
        print(f"{r['family']:9s} {r['start']:>9} {idf:>6.3f} {gtgt[:22]:22s} {gid:>5.3f} {gcov:>5.2f} "
              f"{(depth_ratio or 0):>10.2f} {call:18s}")
    from collections import Counter
    print("\nfinal call tally:", dict(Counter(x["call"] for x in out)))
    json.dump(out, open(f"{OUT}/promoted_final.json", "w"), indent=1)
    print(f"wrote {OUT}/promoted_final.json")


if __name__ == "__main__":
    main()
