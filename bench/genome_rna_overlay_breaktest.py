#!/usr/bin/env python
"""DECISIVE break-test for A2 (the control the adversarial panel flagged as not-yet-run).

A2's headline correlation is Spearman(SEDEF reference %identity, read-PSV density psv_per_kb) = -0.443.
The refuters argued psv_per_kb is read-GATED but reference-SEEDED: a PSV column requires the copies to
DIFFER per the minimap2 asm20 self-alignment of the REFERENCE copy fasta (copy_cols), and reads only gate
that >=3 reads cover the column (psv_graph_genomewide.py:131-141). If so, the read gate carries no
divergence information and the headline is a two-reference-aligner concordance, not a read-content signal.

This script settles it: recompute the PSV-column count WITHOUT the read-coverage gate (count every asm20
column where the copies differ, regardless of read support -- a PURELY reference quantity), and correlate
THAT with SEDEF identity. If rho stays ~-0.44, the gate is INERT => the headline carries zero read
information. (The ungated count needs only the asm20 copy-alignment, no read pileup, so it is fast.)

Run: /home/juanfra/miniforge3/bin/python bench/genome_rna_overlay_breaktest.py
"""
import collections
import csv
import json
import os
import subprocess

import numpy as np
import pysam
from scipy.stats import spearmanr

# reuse the EXACT helpers the A2 build used (SEDEF identity, copy dedup) + psv_graph's column reader
from genome_rna_overlay_readcontent import dedup_copies, load_sedef_index, family_sedef_identity
from psv_graph_genomewide import GENOME, column_alleles

FAM_TSV = "/home/juanfra/winloci_scratch/validated_families.tsv"
PSV_JSON = os.path.join(os.path.dirname(__file__), "psv_graph_genomewide.json")
TMP = "/home/juanfra/winloci_scratch/psvgw_breaktest"
MIN_READS_PSV = 3


def sh(c):
    subprocess.run(c, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def asm20_psv_columns(regions, genome, contig_len):
    """Replicate psv_graph run_family's asm20 column step. Return (ungated, gated_by_reads_loadable?).
    ungated = #columns where ALL copies are defined and >=2 distinct alleles (pure reference asm20),
    NO read-coverage gate. (We do NOT recompute the read gate here; the json carries the gated count.)"""
    os.makedirs(TMP, exist_ok=True)
    regions = sorted(regions, key=lambda r: -(r[2] - r[1]))
    names = [f"copy{i}" for i in range(len(regions))]
    bc, bs, be, _ = regions[0]
    bs, be = max(0, bs), min(contig_len.get(bc, be), be)
    refseq = genome.fetch(bc, bs, be).upper()
    if len(refseq) < 200:
        return None, len(refseq)
    open(f"{TMP}/ref.fa", "w").write(f">{names[0]}\n{refseq}\n")
    sh(f"samtools faidx {TMP}/ref.fa")
    with open(f"{TMP}/copies.fa", "w") as o:
        for nm, (c, s, e, _) in list(zip(names, regions))[1:]:
            s2, e2 = max(0, s), min(contig_len.get(c, e), e)
            o.write(f">{nm}\n{genome.fetch(c, s2, e2).upper()}\n")
    sh(f"minimap2 -ax asm20 {TMP}/ref.fa {TMP}/copies.fa 2>/dev/null "
       f"| samtools sort -o {TMP}/copies.bam - && samtools index {TMP}/copies.bam")
    copy_cols = column_alleles(f"{TMP}/copies.bam", names[0])
    ungated = 0
    for p in range(len(refseq)):
        hap = {names[0]: refseq[p]}
        for nm in names[1:]:
            if nm in copy_cols.get(p, {}):
                hap[nm] = copy_cols[p][nm]
        if len(hap) == len(names) and len(set(hap.values())) >= 2:   # NO read gate
            ungated += 1
    return ungated, len(refseq)


def main():
    print("A2 BREAK-TEST: is the read-coverage gate INERT? (ungated asm20-divergence vs SEDEF identity)\n")
    # families processed by psv_graph (the 154 with a gated psv count)
    J = json.load(open(PSV_JSON))
    gated = {str(f["family"]): f for f in J["families"]}

    loci = collections.defaultdict(list)
    for r in csv.DictReader(open(FAM_TSV), delimiter="\t"):
        loci[r["family"]].append((r["chrom"], int(r["start"]), int(r["end"]), int(r["nreads"])))

    print("[load] SEDEF index ...")
    idx, n_sedef = load_sedef_index()
    genome = pysam.FastaFile(GENOME)
    contig_len = dict(zip(genome.references, genome.lengths))

    rows = []
    for fi in sorted(gated, key=lambda x: int(x)):
        if fi not in loci:
            continue
        regions = dedup_copies(loci[fi])
        if len(regions) < 2:
            continue
        sed, n_link, _, _ = family_sedef_identity(regions, idx)
        if sed is None:
            continue
        ungated, reflen = asm20_psv_columns(regions, genome, contig_len)
        if ungated is None:
            continue
        mean_len = float(np.mean([r[2] - r[1] for r in regions]))
        kb = mean_len / 1000.0
        rows.append(dict(family=fi, sedef=sed, gated_psv=gated[fi]["psvs"],
                         ungated_psv=ungated, mean_len=mean_len,
                         gated_per_kb=gated[fi]["psvs"] / kb, ungated_per_kb=ungated / kb))
        if len(rows) % 25 == 0:
            print(f"   ... {len(rows)} families")

    ident = np.array([r["sedef"] for r in rows])
    gpk = np.array([r["gated_per_kb"] for r in rows])
    upk = np.array([r["ungated_per_kb"] for r in rows])
    gated_psv = np.array([r["gated_psv"] for r in rows], float)
    ungated_psv = np.array([r["ungated_psv"] for r in rows], float)

    rg, pg = spearmanr(ident, gpk)
    ru, pu = spearmanr(ident, upk)
    rgu, _ = spearmanr(gated_psv, ungated_psv)
    frac_kept = float(np.mean(np.divide(gated_psv, ungated_psv, out=np.zeros_like(gated_psv),
                                        where=ungated_psv > 0)))

    out = os.path.join(os.path.dirname(__file__), "genome_rna_overlay_breaktest.json")
    res = dict(n=len(rows),
               gated_per_kb_vs_sedef=dict(rho=rg, p=pg),
               ungated_per_kb_vs_sedef=dict(rho=ru, p=pu),
               gated_vs_ungated_psv=dict(rho=rgu),
               mean_frac_columns_kept_by_read_gate=frac_kept,
               total_gated=int(gated_psv.sum()), total_ungated=int(ungated_psv.sum()))
    json.dump(res, open(out, "w"), indent=2)

    print(f"\n== RESULT (n={len(rows)} families) ==")
    print(f"  GATED   psv_per_kb (read-gated, reproduces build)  vs SEDEF identity: rho={rg:+.3f} p={pg:.2e}")
    print(f"  UNGATED psv_per_kb (NO read gate, pure asm20 ref)  vs SEDEF identity: rho={ru:+.3f} p={pu:.2e}")
    print(f"  gated vs ungated PSV counts: rho={rgu:+.3f}  (read gate keeps {100*frac_kept:.0f}% of asm20 columns on avg)")
    print(f"  total columns: gated {int(gated_psv.sum())} / ungated {int(ungated_psv.sum())}")
    verdict = ("READ GATE INERT -> headline is a TWO-REFERENCE-ALIGNER concordance, carries no read info"
               if abs(ru) >= 0.30 and abs(ru) >= 0.7 * abs(rg)
               else "read gate carries information -> A2 retains some read-content signal")
    print(f"\n  VERDICT: {verdict}")
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
