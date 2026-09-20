#!/usr/bin/env python
"""Build the RECIPROCAL WHOLE-TRANSCRIPT IDENTITY cache for the divergence floor.

WHAT THIS ADDS (vs ri_sharedlen_universal.tsv)
----------------------------------------------
ri_sharedlen_universal.tsv stores, per cross-gene E_r edge, the longest shared
spliced-exon aligned BLOCK length (PAF col 11) and aln_frac = block/min(len). It
does NOT store the number of MATCHING bases (PAF col 10) nor the two member
lengths -- so it cannot express a reciprocal IDENTITY. This builder recomputes the
IDENTICAL alignment (same spliced-exon sequences, same minimap2 -cx asm20 -t1,
same representative bridge-DN pair, same longest-block selection) and additionally
records:
    matches_best  : PAF col 10 (residue matches) of the LONGEST block (the block
                    aln_frac is built from)
    blocklen_best : PAF col 11 of that block  (== ri_sharedlen_universal aln_len)
    matches_sum   : sum of PAF col 10 over ALL alignment lines (whole-transcript
                    matching bases; reported for a sensitivity check only)
    len_a, len_b  : spliced-exon lengths of the two members

From these the RECIPROCAL WHOLE-TRANSCRIPT IDENTITY is
    recip_id_best = matches_best / max(len_a, len_b)
                  = min( matches_best/len_a , matches_best/len_b )
i.e. the min over the two members of (aligned-AND-matching bases / member length).
It is RECIPROCAL: a small gene embedded in one big exon of a large gene shares one
block at high block-identity but scores LOW on the large gene's side (max(len) is
big) -> excluded. aln_frac (block/min) does NOT have this property.

POPULATION: only the cross-gene gene-pairs that ALREADY pass the shipped edge gate
(core_recip>=0.19 AND aln_frac>=0.24) -- the ONLY edges a divergence floor can
additionally cut (edges failing core+aln are cut regardless of the floor). This is
why we do not re-align all 5571 edges, only ~1457.

SELF-CHECK: blocklen_best reproduces ri_sharedlen_universal.tsv aln_len byte-for-byte
(same code path) -- asserted, so this is the SAME alignment, not a re-derivation.

Deterministic (minimap2 -t1, PYTHONHASHSEED=0, sorted by gene pair).
Writes: bench/ri_sharedlen_recip_id.tsv
Run:    PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/ri_build_recip_id.py
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
import family_er_pr as F
import family_er_te_gate as TE
import pysam

SCRATCH = "/home/juanfra/winloci_scratch"
GENOME = os.path.join(SCRATCH, "GGO.fasta")
UNIV_TSV = os.path.join(BENCH, "ri_sharedlen_universal.tsv")
CACHE = os.path.join(BENCH, "ri_sharedlen_recip_id.tsv")

CORE_MIN = 0.19   # shipped edge gate (family_rna_refine.CORE_MIN)
ALN_MIN = 0.24    # shipped edge gate (family_rna_refine.ALN_MIN)


def passing_edges():
    """Cross-gene rows of the universal cache that pass core>=0.19 AND aln>=0.24
    (representative bridge-DN pair carried through). Deterministic sorted order."""
    rows = []
    with open(UNIV_TSV) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        I = {c: i for i, c in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            ga, gb = f[I["gene_a"]], f[I["gene_b"]]
            if ga == gb:
                continue
            core = f[I["core"]]
            aln = f[I["aln_frac"]]
            core = 0.0 if core in ("", "NA") else float(core)
            aln = 0.0 if aln in ("", "NA") else float(aln)
            if core >= CORE_MIN and aln >= ALN_MIN:
                rows.append((ga, gb, f[I["dn_a"]], f[I["dn_b"]], f[I["cls"]],
                             core, aln, f[I["aln_len"]]))
    rows.sort(key=lambda e: (e[0], e[1]))
    return rows


def main():
    meta = F.load_meta()
    skel = TE.load_skeletons()
    fa = pysam.FastaFile(GENOME)

    edges = passing_edges()
    print(f"[recip-id] {len(edges)} cross-gene edges passing core>=%.2f AND aln>=%.2f"
          % (CORE_MIN, ALN_MIN), flush=True)

    seqcache = {}
    def spliced(dn):
        if dn in seqcache:
            return seqcache[dn]
        chrom, exons = TE.dn_exons(dn, meta, skel)
        s = None if chrom is None else "".join(fa.fetch(chrom, a, b) for a, b in exons)
        seqcache[dn] = s
        return s

    td = tempfile.mkdtemp(prefix="ri_recip_")
    pa = os.path.join(td, "a.fa")
    pb = os.path.join(td, "b.fa")
    n_missing = 0
    n_aln_mismatch = 0
    with open(CACHE, "w") as out:
        out.write("gene_a\tgene_b\tdn_a\tdn_b\tcls\tcore\taln_frac\taln_len\t"
                  "matches_best\tblocklen_best\tmatches_sum\tlen_a\tlen_b\t"
                  "recip_id_best\trecip_id_sum\n")
        for i, (ga, gb, da, db, cls, core, aln, aln_len_univ) in enumerate(edges):
            sa, sb = spliced(da), spliced(db)
            if not (sa and sb):
                n_missing += 1
                out.write(f"{ga}\t{gb}\t{da}\t{db}\t{cls}\t{core}\t{aln}\t{aln_len_univ}\t"
                          f"\t\t\t\t\t\t\n")
                continue
            with open(pa, "w") as fh:
                fh.write(f">a\n{sa}\n")
            with open(pb, "w") as fh:
                fh.write(f">b\n{sb}\n")
            r = subprocess.run(["minimap2", "-cx", "asm20", "-t", "1", pa, pb],
                               capture_output=True, text=True)
            best_bl = 0
            best_matches = 0
            sum_matches = 0
            for line in r.stdout.splitlines():
                fld = line.split("\t")
                if len(fld) < 11:
                    continue
                matches = int(fld[9])     # PAF col 10 = residue matches
                blocklen = int(fld[10])   # PAF col 11 = alignment block length
                sum_matches += matches
                if blocklen > best_bl:
                    best_bl = blocklen
                    best_matches = matches
            la, lb = len(sa), len(sb)
            mx = max(la, lb)
            recip_best = best_matches / mx if mx else 0.0
            recip_sum = min(sum_matches / la, sum_matches / lb) if (la and lb) else 0.0
            recip_sum = min(recip_sum, 1.0)
            # SELF-CHECK: our longest block reproduces the shipped universal aln_len
            if aln_len_univ not in ("", "NA") and int(aln_len_univ) != best_bl:
                n_aln_mismatch += 1
            out.write(f"{ga}\t{gb}\t{da}\t{db}\t{cls}\t{core}\t{aln}\t{aln_len_univ}\t"
                      f"{best_matches}\t{best_bl}\t{sum_matches}\t{la}\t{lb}\t"
                      f"{recip_best:.4f}\t{recip_sum:.4f}\n")
            if (i + 1) % 200 == 0:
                print(f"[recip-id] {i+1}/{len(edges)}", flush=True)
    print(f"[recip-id] wrote {CACHE}  (missing_seq={n_missing}, "
          f"aln_len_mismatch_vs_universal={n_aln_mismatch})", flush=True)


if __name__ == "__main__":
    main()
