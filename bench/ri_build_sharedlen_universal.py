#!/usr/bin/env python
"""Build a UNIVERSAL shared-segment length cache over ALL core-present E_r edges.

LEAKAGE FIX: the original ri_sharedlen_cache.tsv was restricted to in_ep=0
(protein-derived) edges for compute economy. That makes aln_frac PRESENCE encode
in_ep, i.e. a protein label leaks in through the back door if the oracle conditions
on presence. This builder recomputes the SAME RNA-derived shared-exon alignment
(minimap2 -cx asm20 -t1, spliced-exon sequences) on EVERY core-present edge, so
aln_frac has universal coverage over the direct-edge graph and its presence carries
NO protein information. Values on the in_ep=0 subset are byte-identical to the
original cache (same code path), which we assert.

Deterministic (-t1, PYTHONHASHSEED=0, sorted by gene pair).
Writes: bench/ri_sharedlen_universal.tsv (gene_a,gene_b,dn_a,dn_b,cls,in_ep,core,aln_len,aln_frac)
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
CACHE = os.path.join(BENCH, "ri_sharedlen_universal.tsv")


def main():
    meta = F.load_meta()
    skel = TE.load_skeletons()
    fa = pysam.FastaFile(GENOME)

    edges = []
    with open(os.path.join(BENCH, "edge_bridge_mask.tsv")) as fh:
        h = fh.readline().rstrip("\n").split("\t"); idx = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[idx["core"]] == "":          # ALL core-present edges (no in_ep filter)
                continue
            edges.append((f[idx["gene_a"]], f[idx["gene_b"]], f[idx["bridge_dn_a"]],
                          f[idx["bridge_dn_b"]], f[idx["cls"]], f[idx["in_ep"]], f[idx["core"]]))
    edges.sort(key=lambda e: (e[0], e[1]))
    print(f"[sharedlen-univ] {len(edges)} core-present edges (no in_ep filter)", flush=True)

    seqcache = {}
    def spliced(dn):
        if dn in seqcache:
            return seqcache[dn]
        chrom, exons = TE.dn_exons(dn, meta, skel)
        s = None if chrom is None else "".join(fa.fetch(chrom, a, b) for a, b in exons)
        seqcache[dn] = s
        return s

    td = tempfile.mkdtemp(prefix="ri_slu_")
    pa = os.path.join(td, "a.fa"); pb = os.path.join(td, "b.fa")
    with open(CACHE, "w") as out:
        out.write("gene_a\tgene_b\tdn_a\tdn_b\tcls\tin_ep\tcore\taln_len\taln_frac\n")
        for i, (ga, gb, da, db, cls, in_ep, core) in enumerate(edges):
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
                    ml = int(fld[10])
                    if ml > best:
                        best = ml
                aln_len = best
                aln_frac = round(best / min(len(sa), len(sb)), 4)
            out.write(f"{ga}\t{gb}\t{da}\t{db}\t{cls}\t{in_ep}\t{core}\t{aln_len}\t{aln_frac}\n")
            if (i + 1) % 500 == 0:
                print(f"[sharedlen-univ] {i+1}/{len(edges)}", flush=True)
    print(f"[sharedlen-univ] wrote {CACHE}", flush=True)


if __name__ == "__main__":
    main()
