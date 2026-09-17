#!/usr/bin/env python3
"""Independent check of lattice_edges.py's clause-2 GENE-BODY re-implementation: for every PAF pair touching the NPIP and
TBC1D3 member-holding L1 groups (plus a random sample of other catalog pairs), write the pair's records in BOTH orientations
to a scratch PAF and run the SHIPPED bench/guided_pipeline.gene_body_chains on it; a pair passes iff >= 1 chain is returned
(chains are only emitted when identity >= 0.80 and aligned >= 0.50 x min(query, extrapolated span)). Compare with
edges.tsv d_c2x_gb_chain_raw (the shipped denominator, before the v-exon-overlap and strand requirements that
bench/denovo_shared_def.py applies to the returned chains; the primary d_c2_genebody additionally uses min(body u, body v)).
Pairs checked: every PAF pair touching a member-holding L1 group of the primary, the c2x variant or the no-v-exon/no-strand
variant (the largest L1 groups), plus 3,000 random other catalog pairs. Output: lattice/check_c2.out

Second check (correction pass 2026-09-16): the target-side requirements. On the SAME shipped chains, apply the shipped
bench/denovo_shared_def.py cmd_families edge loop verbatim in body coordinates (ExonIndex on the two nodes' exon unions,
orient = u strand for a '+' chain and flipped for '-', skip v when both nodes are spliced and v's strand differs from
orient) and compare with edges.tsv d_c2x_genebody (shipped denominator + v-exon overlap + strand check). The primary
d_c2_genebody runs the same target/strand code with the shorter-body denominator.
"""
import collections
import csv
import os
import random
import sys

sys.path.insert(0, "/mnt/c/Users/jfris/Desktop/Rustle/bench")
sys.path.insert(0, "/mnt/c/Users/jfris/Desktop/Rustle/bench/layer_order")
import denovo_shared_def as dsd  # noqa: E402
import guided_pipeline as gp  # noqa: E402
from lattice_common import LIGHT, OUT, PAF, tsv  # noqa: E402

# node exon unions per body key, exactly as lattice_edges.py builds them (exec its head; stops before the pair evaluation)
_src = open("/mnt/c/Users/jfris/Desktop/Rustle/bench/layer_order/lattice_edges.py").read()
_g = {"__name__": "check_c2"}
exec(compile(_src[:_src.index("# ------------------------------------------------------------------------------------------------ DNA pair attributes")],
             "lattice_edges_head", "exec"), _g)
KEY_BLOCKS = {f"{k[0]}:{k[1]}-{k[2]}": v for k, v in _g["key_blocks"].items()}

random.seed(20260916)
grp = {r["gene_id"]: r for r in tsv(f"{OUT}/groups.tsv")}
anch = {r["name"]: r["gene_id"] for r in grp.values()}
L1 = {g: r["primary|L1"] for g, r in grp.items()}
L1x = {g: r["L1=c2x_extrapolated|L1"] for g, r in grp.items()}
L1l = {g: r["L1=c2_no_vexon_no_strand|L1"] for g, r in grp.items()}
fam_groups = {L1[anch["NPIPB2"]], L1[anch["TBC1D3"]]}
fam_groups_x = {L1x[anch["NPIPB2"]], L1x[anch["TBC1D3"]]}
fam_groups_l = {L1l[anch["NPIPB2"]], L1l[anch["TBC1D3"]]}
genes = {r["gene_id"]: r for r in tsv("/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/light/work/refseq/genes.tsv")}
key_of_gene = {g: f"{r['chrom']}:{int(r['start0']) + 1}-{r['end']}" for g, r in genes.items()}

want = {}
want_t = {}
other = []
with open(f"{OUT}/edges.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["d_evidence"] != "yes":
            continue
        k = frozenset((key_of_gene[r["gene_a"]], key_of_gene[r["gene_b"]]))
        if (L1[r["gene_a"]] in fam_groups or L1[r["gene_b"]] in fam_groups or L1x[r["gene_a"]] in fam_groups_x
                or L1x[r["gene_b"]] in fam_groups_x or L1l[r["gene_a"]] in fam_groups_l or L1l[r["gene_b"]] in fam_groups_l):
            want[k] = r["d_c2x_gb_chain_raw"]
            want_t[k] = (r["d_c2x_genebody"], genes[r["gene_a"]]["strand"], genes[r["gene_b"]]["strand"], key_of_gene[r["gene_a"]])
        else:
            other.append((k, r["d_c2x_gb_chain_raw"], (r["d_c2x_genebody"], genes[r["gene_a"]]["strand"],
                                                        genes[r["gene_b"]]["strand"], key_of_gene[r["gene_a"]])))
for k, v, vt in random.sample(other, min(3000, len(other))):
    want[k] = v
    want_t[k] = vt
lines = collections.defaultdict(list)
for cat, path in PAF.items():
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == f[5]:
                continue
            k = frozenset((f[0], f[5]))
            if k in want:
                lines[k].append(f)
tmp = f"{OUT}/check_c2.tmp.paf"
agree = disagree = 0
agree_t = disagree_t = 0
examples = []
FLIPS = {"+": "-", "-": "+"}
for k, exp in want.items():
    with open(tmp, "w") as out:
        for f in lines[k]:
            out.write("\t".join(f) + "\n")
            sw = [f[5], f[6], f[7], f[8], f[4], f[0], f[1], f[2], f[3]] + f[9:]
            out.write("\t".join(sw) + "\n")
    chains = gp.gene_body_chains(tmp)
    got = "yes" if chains else "no"
    if got == exp:
        agree += 1
    else:
        disagree += 1
        if len(examples) < 10:
            examples.append((sorted(k), exp, got))
    # shipped cmd_families edge loop on these chains (body coordinates; chrom field = the target body key)
    exp_t, strand_a, strand_b, key_a = want_t[k]
    ka, kb = sorted(k)
    strand = {key_a: strand_a, (kb if key_a == ka else ka): strand_b}
    nd = [{"idx": i, "chrom": kk, "strand": strand[kk], "exons": [(s0 - (int(kk.rsplit(":", 1)[1].split("-")[0]) - 1),
                                                                    e0 - (int(kk.rsplit(":", 1)[1].split("-")[0]) - 1))
                                                                   for s0, e0 in KEY_BLOCKS[kk]]} for i, kk in enumerate((ka, kb))]
    idx = dsd.ExonIndex(nd)
    spliced = {n["idx"]: len(n["exons"]) >= 2 for n in nd}
    by_key = {n["chrom"]: n for n in nd}
    edge = False
    for c in chains:
        u = by_key[c["q"]]
        orient = u["strand"] if c["strand"] == "+" else FLIPS.get(u["strand"], u["strand"])
        for v in idx.hits(c["chrom"], c["s"], c["e"]):
            if v == u["idx"]:
                continue
            if spliced[u["idx"]] and spliced[v] and nd[v]["strand"] != orient:
                continue
            edge = True
    got_t = "yes" if edge else "no"
    if got_t == exp_t:
        agree_t += 1
    else:
        disagree_t += 1
        if len(examples) < 20:
            examples.append(("target/strand", sorted(k), exp_t, got_t))
os.remove(tmp)
msg = [f"pairs checked {len(want)} (member-holding L1 groups: {len(want) - min(3000, len(other))}; random other catalog pairs "
       f"{min(3000, len(other))}); raw chain exists (d_c2x_gb_chain_raw): agree {agree}; disagree {disagree}",
       f"shipped cmd_families v-exon overlap + strand check on those chains vs d_c2x_genebody: agree {agree_t}; disagree "
       f"{disagree_t}"] + [f"  example {e}" for e in examples]
print("\n".join(msg))
open(f"{OUT}/check_c2.out", "w").write("\n".join(msg) + "\n")
