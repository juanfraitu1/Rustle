#!/usr/bin/env python
"""Build the SHARED-SEGMENT LENGTH cache (the caveat fix) for the read-intrinsic TE gate.

The honest caveat: high-copy REAL gene families (KRAB-ZNF, olfactory) are ALSO high homology-degree /
high k-mer multiplicity, so degree/multiplicity alone cannot separate them from TEs. The discriminator is
the LENGTH of the SHARED contiguous segment: a real paralog shares a LONG homologous gene body; a TE bridge
shares only a SHORT exonized-TE fragment. For each arbitration edge (in_ep=0, core present) align the two
bridge-DN spliced exon sequences (minimap2 -cx asm20 -t1, deterministic) and record the longest aligned
block length (abs bp) and its fraction of the shorter locus. Diagnostic; deterministic (-t1).

Reproduce: PYTHONHASHSEED=0 python bench/ri_build_sharedlen_cache.py
Writes: bench/ri_sharedlen_cache.tsv  (gene_a, gene_b, dn_a, dn_b, cls, in_ep, core, aln_len, aln_frac)
"""
import os, sys, subprocess, tempfile
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"; os.execv(sys.executable, [sys.executable] + sys.argv)
BENCH = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, BENCH)
import family_er_pr as F
import family_er_te_gate as TE
import pysam

SCRATCH = "/home/juanfra/winloci_scratch"
GENOME = os.path.join(SCRATCH, "GGO.fasta")
CACHE = os.path.join(BENCH, "ri_sharedlen_cache.tsv")


def main():
    meta = F.load_meta()
    skel = TE.load_skeletons()
    fa = pysam.FastaFile(GENOME)

    # arbitration edges (in_ep=0, core present) with their bridge DN pair
    edges = []
    with open(os.path.join(BENCH, "edge_bridge_mask.tsv")) as fh:
        h = fh.readline().rstrip("\n").split("\t"); idx = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[idx["core"]] == "" or f[idx["in_ep"]] != "0":
                continue
            edges.append((f[idx["gene_a"]], f[idx["gene_b"]], f[idx["bridge_dn_a"]],
                          f[idx["bridge_dn_b"]], f[idx["cls"]], f[idx["core"]]))
    print(f"[sharedlen] {len(edges)} arbitration edges (in_ep=0, core present)", flush=True)

    seqcache = {}
    def spliced(dn):
        if dn in seqcache:
            return seqcache[dn]
        chrom, exons = TE.dn_exons(dn, meta, skel)
        s = None if chrom is None else "".join(fa.fetch(chrom, a, b) for a, b in exons)
        seqcache[dn] = s
        return s

    td = tempfile.mkdtemp(prefix="ri_sl_")
    pa = os.path.join(td, "a.fa"); pb = os.path.join(td, "b.fa")
    with open(CACHE, "w") as out:
        out.write("gene_a\tgene_b\tdn_a\tdn_b\tcls\tin_ep\tcore\taln_len\taln_frac\n")
        for i, (ga, gb, da, db, cls, core) in enumerate(edges):
            sa, sb = spliced(da), spliced(db)
            aln_len, aln_frac = "", ""
            if sa and sb:
                with open(pa, "w") as fh:
                    fh.write(f">a\n{sa}\n")
                with open(pb, "w") as fh:
                    fh.write(f">b\n{sb}\n")
                r = subprocess.run(["minimap2", "-cx", "asm20", "-t", "1", pa, pb],
                                   capture_output=True, text=True)
                best = 0
                for line in r.stdout.splitlines():
                    fld = line.split("\t")
                    if len(fld) < 11:
                        continue
                    ml = int(fld[10])   # PAF col 11 = alignment block length
                    if ml > best:
                        best = ml
                aln_len = best
                aln_frac = round(best / min(len(sa), len(sb)), 4)
            out.write(f"{ga}\t{gb}\t{da}\t{db}\t{cls}\t0\t{core}\t{aln_len}\t{aln_frac}\n")
            if (i + 1) % 300 == 0:
                print(f"[sharedlen] {i+1}/{len(edges)}", flush=True)
    print(f"[sharedlen] wrote {CACHE}", flush=True)


if __name__ == "__main__":
    main()
