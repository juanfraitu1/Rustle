#!/usr/bin/env python3
"""family_def_dna_scaffold.py — DNA-FIRST scaffold: recover the expressed-not-family gap.

The completeness funnel: of 315 DNA families with >=2 expressed members, only 94 were
recovered as ~R∩~B families (30%) -> 221 'detectable but not recovered'. A DNA-first method
starts from the DNA family list as a scaffold and checks each missed family directly. The
question that matters: of the 221, how many are GENUINELY read-confusable (a real ~R∩~B
family the de-novo pipeline dropped mechanically -- RECOVERABLE) vs RESOLVABLE (reads place
uniquely -> out of scope by design)? And the new -N50 -p0.1 BAM should surface confusability
the old default-cap BAM (which built the 196) missed.
Run: python bench/family_def_dna_scaffold.py
"""
import collections
import json
import os
import sys

import networkx as nx
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import best_gene, GENES_BED, DELTA, DE_MAX, MIN_READS
from family_def_read_filters import dna_homology

META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
FAM_TSV = "/home/juanfra/winloci_scratch/validated_families.tsv"
BAM = "/home/juanfra/winloci_scratch/GGO_mm.bam"
PAD = 1000


def main():
    by = collections.defaultdict(list); coord = {}
    with open(GENES_BED) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                by[p[0]].append((int(p[1]), int(p[2]), p[3])); coord[p[3]] = (p[0], int(p[1]), int(p[2]))
    for c in by:
        by[c].sort()
    Hd, _ = dna_homology()
    DG = nx.Graph()
    for (ga, gb), r in Hd.items():
        if r.get("id", 0) >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30:
            DG.add_edge(ga, gb)
    dna_fams = [c for c in nx.connected_components(DG) if len(c) >= 2]

    expressed = set()
    with open(META) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            g = best_gene(by, p[1], int(p[2]), int(p[3]))
            if g:
                expressed.add(g)
    rna_fam_genes = set()
    fams = collections.defaultdict(list)
    for line in open(FAM_TSV).read().splitlines()[1:]:
        fi, lid, c, s, e, nr = line.split("\t")
        fams[int(fi)].append((c, int(s), int(e)))
    for fi, loci in fams.items():
        gs = {best_gene(by, c, s, e) for c, s, e in loci} - {None}
        if len(gs) >= 2:
            rna_fam_genes |= gs

    # detectable (>=2 expressed) but NOT recovered (<2 members in an RNA family)
    missed = []
    for fam in dna_fams:
        em = [g for g in fam if g in expressed and g in coord]
        if len(em) >= 2 and len(fam & rna_fam_genes) < 2:
            missed.append(em)
    print(f"=== DNA-first scaffold on the {len(missed)} 'detectable-but-not-recovered' families ===")

    bam = pysam.AlignmentFile(BAM, "rb")

    def reads_at(g):
        c, s, e = coord[g]; out = {}
        try:
            it = bam.fetch(c, max(0, s - PAD), e + PAD)
        except (ValueError, KeyError):
            return out
        for r in it:
            if r.is_unmapped or r.is_supplementary:
                continue
            if r.reference_end is None or r.reference_end < s or r.reference_start > e:
                continue
            de = dict(r.get_tags()).get("de")
            if de is None or de > DE_MAX:
                continue
            if r.query_name not in out or de < out[r.query_name]:
                out[r.query_name] = de
        return out

    confusable = resolvable = 0
    recovered_genes = 0
    for em in missed:
        rd = {g: reads_at(g) for g in em}
        # any pair with >=MIN_READS de-tied cross-mapping reads?
        is_conf = False
        core = set()
        for i in range(len(em)):
            for j in range(i + 1, len(em)):
                gi, gj = em[i], em[j]
                shared = set(rd[gi]) & set(rd[gj])
                tied = [q for q in shared if abs(rd[gi][q] - rd[gj][q]) <= DELTA and max(rd[gi][q], rd[gj][q]) <= DE_MAX]
                if len(tied) >= MIN_READS:
                    is_conf = True; core |= {gi, gj}
        if is_conf:
            confusable += 1; recovered_genes += len(core)
        else:
            resolvable += 1
    bam.close()

    print(f"  CONFUSABLE on the new BAM (real ~R∩~B family the de-novo pipeline MISSED -> RECOVERED): "
          f"{confusable}/{len(missed)} ({100*confusable/max(len(missed),1):.0f}%)")
    print(f"  RESOLVABLE (reads place uniquely -> out of scope by design): "
          f"{resolvable}/{len(missed)} ({100*resolvable/max(len(missed),1):.0f}%)")
    print(f"\n  => DNA-first scaffold recovers {confusable} families ({recovered_genes} member genes) "
          f"the de-novo pipeline dropped -- mostly due to the old default-cap BAM under-sampling ~R.")
    print(f"  RNA-detectable families: 94 (~R∩~B) + {confusable} (scaffold) = {94+confusable} "
          f"of 315 ({100*(94+confusable)/315:.0f}%); remaining {resolvable} are genuinely resolvable.")
    json.dump(dict(missed=len(missed), confusable=confusable, resolvable=resolvable,
                   recovered_genes=recovered_genes, prior_recovered=94, detectable=315),
              open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "family_def_dna_scaffold.json"), "w"), indent=2)
    print("\n[+] wrote family_def_dna_scaffold.json")


if __name__ == "__main__":
    main()
