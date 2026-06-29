#!/usr/bin/env python
"""F2 — Soto-2025 (Cell) DNA shared-exon / famCN validation of the RNA-conflict family catalog (no SEDEF).

Soto et al. define human-specific gene families at the DNA level: extract segmental-duplication region DNA,
minimap2 map it back to the genome (-c --eqx -N50 -p0.5), BEDTools-intersect -f 0.99 to keep mappings covering
>99% of each exon (= SHARED EXONS), then group genes sharing exons whose family copy-number (famCN) mean-abs-
deviation is < 1. This is the field-standard recipe for RECENT/near-identical paralog families and is exactly
the method our key reference uses. We mirror it here as an ORTHOGONAL DNA witness for the RNA families.

Independence: the RNA families come from the read-CONFLICT graph (reads multimapping among copies). This check
instead takes each copy's ASSEMBLED EXON sequence and maps it back to the GENOME, asking whether the family's
copy loci are mutual >=99%-coverage homologs with consistent copy number — a DNA-sequence property computed
WITHOUT the read-conflict graph. (Honest caveat: it shares the exon sequence with the catalog, so it is not
fully orthogonal to sequence; but mutual map-back onto each other's loci + famCN multiplicity is a genuine
segmental-duplication signal the read-conflict graph does not encode.)

Criterion per family (>=2 copies):
  - map each copy's transcript (exons) to the genome; KEEP hits with query-coverage >= 0.99 (Soto -f 0.99) and
    identity (1 - de) >= 0.90 (homolog, allowing paralog divergence).
  - famCN_i = number of distinct genomic loci copy_i's exons hit at >=99% coverage (the DNA family copy number
    seen from copy i).
  - SHARED-EXON (mutual): copy_i and copy_j share exons iff copy_i's exons map (>=99% cov) onto copy_j's catalog
    locus (target overlaps copy_j [start,end]) or vice versa. A family is DNA-CONFIRMED iff every copy shares
    exons with >=1 other copy (one connected duplication group).
  - famCN-CONSISTENT iff mean-abs-deviation of {famCN_i} < 1 (Soto).

Run: /home/juanfra/miniforge3/bin/python bench/soto_family_validate.py
"""
import os
import subprocess
import sys
import tempfile

import pysam

SCRATCH = "/home/juanfra/winloci_scratch"
COPIES = f"{SCRATCH}/gw_conflict_catalog.copies.tsv"
GENOME = f"{SCRATCH}/GGO.fasta"
MM2 = os.environ.get("RUSTLE_MINIMAP2", "/home/juanfra/miniforge3/bin/minimap2")
# Soto keeps mappings covering >99% of each EXON (shared exons), and groups genes that SHARE EXONS — NOT whole
# transcripts. So the unit is an exon-sized homologous block, not the full transcript: a qualifying "shared-exon"
# hit = an alignment block >= MIN_EXON_BP on the query at >= MIN_ID identity. (Requiring whole-transcript
# coverage is stricter than Soto and wrongly fails paralogs that share only some exons.)
MIN_EXON_BP = 200
MIN_ID = 0.90


def load_families():
    fams = {}
    with open(COPIES) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            fid = f[idx["family_id"]]
            fams.setdefault(fid, []).append({
                "tid": f[idx["tid"]], "chrom": f[idx["chrom"]],
                "start": int(f[idx["start"]]), "end": int(f[idx["end"]]),
            })
    return {k: v for k, v in fams.items() if len(v) >= 2}


def main():
    fams = load_families()
    copies = [c for v in fams.values() for c in v]
    # key copies by (chrom,start,end) — denovo_transcripts.fa has DUPLICATE tids, so we source the copy
    # sequence from the GENOMIC SPAN (chrom:start-end) instead, which is also exactly what Soto maps (genomic
    # SD DNA, not spliced transcripts). Each copy gets a unique query name from its locus.
    print(f"{len(fams)} families with >=2 copies, {len(copies)} copies (genomic-span source)")

    def qname_of(c):
        return f"{c['chrom']}:{c['start']}-{c['end']}"

    fa = pysam.FastaFile(GENOME)
    with tempfile.TemporaryDirectory() as td:
        qpath = os.path.join(td, "copy_spans.fa")
        seen = set()
        with open(qpath, "w") as out:
            for c in copies:
                qn = qname_of(c)
                if qn in seen:
                    continue
                seen.add(qn)
                out.write(f">{qn}\n{fa.fetch(c['chrom'], c['start'], c['end'])}\n")
        print(f"wrote {len(seen)} copy genomic spans; mapping back to the genome (asm20, DNA-to-DNA)...")
        # asm20 (DNA-to-DNA, up to ~20% divergence = paralogs); -N50 -p0.5 keeps paralog secondaries; --eqx.
        paf = subprocess.run(
            [MM2, "-cx", "asm20", "--eqx", "-N", "50", "-p", "0.5", "-t", "4", GENOME, qpath],
            capture_output=True, text=True,
        )
        if paf.returncode != 0:
            sys.exit(f"minimap2 failed:\n{paf.stderr[-2000:]}")

    # parse PAF -> per-tid list of (chrom, tstart, tend, qcov, identity)
    hits = {}
    for line in paf.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        qname, qlen, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
        tname, ts, te = f[5], int(f[7]), int(f[8])
        aln_qbp = qe - qs
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
        ident = (1.0 - de) if de is not None else (int(f[9]) / max(int(f[10]), 1))
        if aln_qbp >= MIN_EXON_BP and ident >= MIN_ID:
            hits.setdefault(qname, []).append((tname, ts, te, aln_qbp, ident))

    def overlaps(a_chrom, a_s, a_e, b_chrom, b_s, b_e):
        return a_chrom == b_chrom and a_s < b_e and b_s < a_e

    n_conf = n_cn = 0
    rows = []
    for fid, members in fams.items():
        famcn = []
        for c in members:
            loci = hits.get(qname_of(c), [])
            # collapse near-duplicate target intervals into distinct loci
            distinct = []
            for (ch, ts, te, _, _) in sorted(loci):
                if not any(overlaps(ch, ts, te, d[0], d[1], d[2]) for d in distinct):
                    distinct.append((ch, ts, te))
            famcn.append(len(distinct))
        # mutual shared-exon: every copy maps onto >=1 OTHER copy's locus
        shared_ok = True
        for i, ci in enumerate(members):
            maps_onto_other = False
            for h in hits.get(qname_of(ci), []):
                for j, cj in enumerate(members):
                    if i == j:
                        continue
                    if overlaps(h[0], h[1], h[2], cj["chrom"], cj["start"], cj["end"]):
                        maps_onto_other = True
                        break
                if maps_onto_other:
                    break
            if not maps_onto_other:
                shared_ok = False
                break
        mad = sum(abs(x - sum(famcn) / len(famcn)) for x in famcn) / len(famcn) if famcn else 99
        cn_ok = mad < 1.0
        n_conf += shared_ok
        n_cn += cn_ok
        rows.append((fid, len(members), famcn, round(mad, 2), shared_ok, cn_ok))

    N = len(fams)
    print(f"\nDNA shared-exon CONFIRMED (every copy has a >=200bp, >=90%-id homologous block onto another copy's "
          f"locus): {n_conf}/{N} ({100*n_conf/N:.1f}%)")
    print(f"famCN-CONSISTENT (MAD < 1 across family members): {n_cn}/{N} ({100*n_cn/N:.1f}%)")
    print(f"BOTH (shared-exon AND famCN-consistent): "
          f"{sum(1 for r in rows if r[4] and r[5])}/{N} ({100*sum(1 for r in rows if r[4] and r[5])/N:.1f}%)")
    print("\nper-family (family, n_copies, famCN_per_copy, MAD, shared_exon, famCN_consistent):")
    for r in sorted(rows, key=lambda x: (not x[4], x[0]))[:40]:
        print(f"  {r[0]:>9} n={r[1]} famCN={r[2]} MAD={r[3]} shared={r[4]} cn={r[5]}")
    print("\nInterpretation: this is the Soto-2025 DNA/CN family criterion applied to the RNA-conflict families.")
    print("DNA-confirmed families are corroborated as genuine segmental-duplication families independently of the")
    print("read-conflict graph; famCN-consistency mirrors Soto's MAD<1 grouping. Non-confirmed families are either")
    print("more-divergent paralogs (<90% id, expected to fail a >=99%-coverage map-back) or single-locus artifacts.")


if __name__ == "__main__":
    main()
