#!/usr/bin/env python
"""Build a READ-CONSENSUS shared-segment cache to prove aln_frac's separation
does NOT hinge on genome-assembly bases (adversarial-review fix).

CONTEXT
-------
The deployed aln_frac (ri_sharedlen_universal.tsv) fetches the SPLICED-EXON
sequence from the genome assembly (GGO.fasta) at each DN locus's exon coordinates.
The splice/locus STRUCTURE (which exons, in what order) is RNA-assembly-derived,
but the base identities are genome-sourced. A maximally strict RNA-primary examiner
is entitled to ask: does the dominant separator survive if the bases come from the
READ CONSENSUS instead of the genome?

This builder recomputes the SAME longest-shared-block fraction, but over the
DE-NOVO READ-ASSEMBLED TRANSCRIPT CONSENSUS sequences (denovo_transcripts.fa, the
POA consensus over the reads assigned to each DN locus -- the same source
core_recip is built from), NOT genome bases. Concordance with the genome-fetched
aln_frac then bounds the "genome-base" concern quantitatively:
  * high rho  => the feature is a property of the transcript, not the base source;
  * high KEEP/CUT agreement at the deployed thresholds => the family catalog is
    robust to sequence provenance.

Deterministic (minimap2 -t1, PYTHONHASHSEED=0, sorted by gene pair).
Writes: bench/ri_sharedlen_readcons.tsv
        (gene_a,gene_b,dn_a,dn_b,cls,in_ep,core,aln_len_read,aln_frac_read)
"""
import os
import sys
import subprocess
import tempfile

if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)
import pysam

SCRATCH = "/home/juanfra/winloci_scratch"
DENOVO_FA = os.path.join(SCRATCH, "denovo_transcripts.fa")   # READ-assembled consensus (POA)
EBM_TSV = os.path.join(BENCH, "edge_bridge_mask.tsv")
CACHE = os.path.join(BENCH, "ri_sharedlen_readcons.tsv")


def main():
    fa = pysam.FastaFile(DENOVO_FA)
    present = set(fa.references)

    edges = []
    with open(EBM_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[idx["core"]] == "":          # ALL core-present edges (identical pop to universal cache)
                continue
            edges.append((f[idx["gene_a"]], f[idx["gene_b"]], f[idx["bridge_dn_a"]],
                          f[idx["bridge_dn_b"]], f[idx["cls"]], f[idx["in_ep"]], f[idx["core"]]))
    edges.sort(key=lambda e: (e[0], e[1]))
    print(f"[sharedlen-readcons] {len(edges)} core-present edges", flush=True)

    seqcache = {}
    def consensus(dn):
        if dn in seqcache:
            return seqcache[dn]
        s = fa.fetch(dn) if dn in present else None    # READ-assembled transcript consensus
        seqcache[dn] = s
        return s

    td = tempfile.mkdtemp(prefix="ri_slrc_")
    pa = os.path.join(td, "a.fa")
    pb = os.path.join(td, "b.fa")
    n_missing = 0
    with open(CACHE, "w") as out:
        out.write("gene_a\tgene_b\tdn_a\tdn_b\tcls\tin_ep\tcore\taln_len_read\taln_frac_read\n")
        for i, (ga, gb, da, db, cls, in_ep, core) in enumerate(edges):
            sa, sb = consensus(da), consensus(db)
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
            else:
                n_missing += 1
            out.write(f"{ga}\t{gb}\t{da}\t{db}\t{cls}\t{in_ep}\t{core}\t{aln_len}\t{aln_frac}\n")
            if (i + 1) % 500 == 0:
                print(f"[sharedlen-readcons] {i+1}/{len(edges)}", flush=True)
    print(f"[sharedlen-readcons] wrote {CACHE} ({n_missing} edges with a missing consensus)", flush=True)


if __name__ == "__main__":
    main()
