"""The built-in ladder: the advisor's progression from identical copies to structural changes, each rung differing from
the previous in exactly the stated way. Every rung shares one template, one background and one decoy set (resolved once
and saved as files, so the rungs are comparable).
"""
import copy
import json
import os

SNP_LADDER = (0.001, 0.005, 0.01, 0.02, 0.05, 0.10)


def rungs(n_exons):
    """[(name, condition, copies)] — `copies` are spec entries for A and A'. `n_exons` picks the middle exon etc.;
    rungs that need more exons than the template has are left out (noted by `ladder`)."""
    mid = max(2, n_exons // 2 + 1)                       # 1-based, never terminal when n >= 3
    out = [("identical", "A and A' 100% identical", [{"id": "A"}, {"id": "A2"}])]
    for d in SNP_LADDER:
        out.append((f"snp_{d:g}", f"A' at divergence {d:g}", [{"id": "A"}, {"id": "A2", "ops": [{"op": "snp", "rate": d}]}]))
    D = {"op": "snp", "rate": 0.02}
    if n_exons >= 3:
        out.append(("exon_loss", f"A' lacks exon {mid} (DNA + RNA), d=0.02", [{"id": "A"}, {"id": "A2", "ops": [D, {"op": "exon_delete", "exon": mid}]}]))
        out.append(("splice_kill", f"A' skips exon {mid} (RNA only), d=0.02", [{"id": "A"}, {"id": "A2", "ops": [D, {"op": "splice_kill", "exon": mid}]}]))
        out.append(("inv_exon", f"exon {mid} of A' inverted, d=0.02", [{"id": "A"}, {"id": "A2", "ops": [D, {"op": "invert", "exon": mid}]}]))
    out.append(("exon_gain", "A' has a new 120 bp exon after exon 1, d=0.02", [{"id": "A"}, {"id": "A2", "ops": [D, {"op": "exon_insert", "after": 1, "length": 120}]}]))
    out.append(("exon_dup", "A' carries exon 2 twice, d=0.02", [{"id": "A"}, {"id": "A2", "ops": [D, {"op": "exon_insert", "after": 1, "source_exon": 2 if n_exons >= 2 else 1}]}]))
    if n_exons >= 4:
        out.append(("shuffle", "exons 2 and 3 of A' swapped, d=0.02", [{"id": "A"}, {"id": "A2", "ops": [D, {"op": "exon_shuffle", "a": 2, "b": 3}]}]))
    out.append(("inv_intron", "intron 1 of A' inverted (silent at RNA), d=0.02", [{"id": "A"}, {"id": "A2", "ops": [D, {"op": "invert", "intron": 1}]}]))
    out.append(("inv_whole", "A' on the other strand, d=0.02", [{"id": "A"}, {"id": "A2", "strand": "-", "ops": [D]}]))
    if n_exons >= 4:
        out.append(("truncated", "A' lacks the first 2 exons, d=0.02", [{"id": "A"}, {"id": "A2", "ops": [D, {"op": "truncate", "side": 5, "exons": 2}]}]))
    if n_exons >= 2:
        out.append(("conversion", "A' at d=0.05 with exon 2 converted from A", [{"id": "A"}, {"id": "A2", "ops": [{"op": "snp", "rate": 0.05}, {"op": "convert", "from": "A", "exon": 2}]}]))
    out.append(("three_copies", "A, A' (0.01), A'' (0.03)", [{"id": "A"}, {"id": "A2", "ops": [{"op": "snp", "rate": 0.01}]}, {"id": "A3", "ops": [{"op": "snp", "rate": 0.03}]}]))
    out.append(("dispersed", "A' on a second contig, d=0.02", [{"id": "A"}, {"id": "A2", "contig": "sim2", "ops": [D]}]))
    out.append(("unexpressed", "A' has no reads, d=0.02 (semi-guided case)", [{"id": "A"}, {"id": "A2", "expression": 0, "ops": [D]}]))
    out.append(("absent", "A' absent from the reference, d=0.02 (O3 case)", [{"id": "A"}, {"id": "A2", "in_reference": False, "ops": [D]}]))
    if n_exons >= 5:
        out.append(("combined", f"A' d=0.03, exon {mid} deleted, intron 1 inverted, last exon lost",
                    [{"id": "A"}, {"id": "A2", "ops": [{"op": "snp", "rate": 0.03}, {"op": "exon_delete", "exon": mid}, {"op": "invert", "intron": 1},
                                                      {"op": "truncate", "side": 3, "exons": 1}]}]))
    return out


def ladder_specs(base, template_json, decoys_json, n_exons, only=None):
    """Specs for every rung, sharing the saved template / decoys files."""
    specs = []
    for name, cond, copies in rungs(n_exons):
        if only and name not in only:
            continue
        s = copy.deepcopy(base)
        s["name"] = name
        s["condition"] = cond
        s["template"] = {"source": "file", "path": template_json}
        s["decoys"] = {"source": "file", "path": decoys_json}
        s["copies"] = copies
        specs.append(s)
    return specs


TEMPLATE_SPEC = {
    "name": "example", "seed": 1,
    "background": {"source": "fasta", "path": "/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/chr20.fa", "region": "chr20:20000000-20300000"},
    "template": {"source": "annotation", "genome": "/mnt/linuxdisk/home/juanfraitu/winloci_data/Reference/chm13v2.0.fa",
                 "annotation": "/mnt/linuxdisk/home/juanfraitu/winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz",
                 "gene": "NPIPA1", "intron_cap": 3000},
    "copies": [{"id": "A", "ops": []},
               {"id": "A2", "ops": [{"op": "snp", "rate": 0.02}, {"op": "exon_delete", "exon": "middle"}]},
               {"id": "A3", "strand": "-", "in_reference": False, "ops": [{"op": "snp", "rate": 0.01}, {"op": "invert", "intron": 1}]}],
    "decoys": {"n": 3, "min_exons": 3, "max_span": 30000, "chrom": "chr20"},
    "reads": {"per_copy": 50, "err": 0.001, "indel": 0.0003, "jitter": 30, "trunc5_frac": 0.0, "trunc5_max": 0.3},
}
