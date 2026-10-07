#!/usr/bin/env python3
"""E0 observability, the reachable stratum R and the instrument controls, from the truth tables and the mapped BAM alone (docs/PREREG_ideal_expression_2026-10-06.md).
Run after mapping and BEFORE the pipeline: R and its sha1 are fixed here.

    strata.py --truth TRUTHPREFIX --bam SIM.bam --out OUTPREFIX

Pools: P1 = primary alignments; P2 = P1 + secondary alignments with AS >= 0.98 x the read's best AS (the pipeline's seeding pool); P3 = any alignment.
A chain is OBSERVABLE iff >= 3 reads of P2 carry exactly that whole chain (the read's junctions, N >= 50 bp, inside the copy span equal the chain), whatever the read's source.
R = copies of stratum R_in (not entangled, not chainless, not a readthrough image) with an observable chain. G2: share of the reads of single-copy genes (no read with a secondary
alignment of AS >= 0.9 x best; not a target, not a named-family gene) whose primary alignment overlaps the source gene (>= 99% required).
Writes OUT.E0.tsv, OUT.R.tsv, OUT.R.sha1, OUT.single_copy.tsv, OUT.G2.json.
"""
import argparse
import collections
import csv
import hashlib
import json

import pysam

MIN_INTRON = 50
GOOD = 0.98
UNIQ = 0.9


def blocks(s):
    return [tuple(map(int, x.split("-"))) for x in s.split(",")] if s else []


def parse_chain(s):
    return tuple(tuple(map(int, x.split("-"))) for x in s.split(",")) if s else ()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
    trs = rd(a.truth + ".transcripts.tsv")
    targets = rd(a.truth + ".targets.tsv")
    named = {r["gene"] for r in rd(a.truth + ".named_family.tsv")}
    tx_copy = {r["transcript"]: r["copy"] for r in trs if r["copy"]}
    tx_gene = {r["transcript"]: r["gene"] for r in trs}
    gene_blocks = collections.defaultdict(list)
    chains = collections.defaultdict(set)
    for r in trs:
        if int(r["simulated"]):
            gene_blocks[r["gene"]].extend(blocks(r["canon_blocks"]))
            if r["copy"] and r["chain"]:
                chains[r["copy"]].add(parse_chain(r["chain"]))
    gene_span = {g: (min(s for s, _ in b), max(e for _, e in b)) for g, b in gene_blocks.items()}
    chrom = targets[0]["chrom"]
    span = {t["cid"]: gene_span[t["gene"]] for t in targets if t["gene"] in gene_span}
    # ---- one pass over the BAM: every alignment
    best = {}
    alns = []
    bam = pysam.AlignmentFile(a.bam)
    all_reads = set()   # every read of the FASTQ: one primary-or-unmapped record each (G3), mapped or not
    for rec in bam.fetch(until_eof=True):
        if not (rec.is_secondary or rec.is_supplementary):
            all_reads.add(rec.query_name)
        if rec.is_unmapped:
            continue
        as_ = rec.get_tag("AS") if rec.has_tag("AS") else 0
        q = rec.query_name
        if not rec.is_supplementary:
            best[q] = max(best.get(q, as_), as_)
        js, pos = [], rec.reference_start
        for op, ln in rec.cigartuples:
            if op in (0, 7, 8, 2):
                pos += ln
            elif op == 3:
                if ln >= MIN_INTRON:
                    js.append((pos, pos + ln))
                pos += ln
        cat = "P" if not (rec.is_secondary or rec.is_supplementary) else ("S" if rec.is_secondary else "U")
        alns.append((q, cat, rec.reference_name, rec.reference_start, rec.reference_end, as_, tuple(js)))
    # ---- E0 per copy and chain
    carry = {k: collections.defaultdict(set) for k in ("P1", "P2", "P3")}   # pool -> (cid, chain) -> read names
    for q, cat, rname, s, e, as_, js in alns:
        if rname != chrom:
            continue
        in_p1 = cat == "P"
        in_p2 = in_p1 or (cat == "S" and as_ >= GOOD * best[q])
        for cid, (s0, e0) in span.items():
            if s < e0 and e > s0 and chains.get(cid):
                inj = tuple(j for j in js if j[0] >= s0 and j[1] <= e0)
                if inj in chains[cid]:
                    key = (cid, inj)
                    carry["P3"][key].add(q)
                    if in_p2:
                        carry["P2"][key].add(q)
                    if in_p1:
                        carry["P1"][key].add(q)
    # reads of each copy's own transcripts whose primary overlaps the copy span (unmapped reads are in the denominator)
    n_reads, back = collections.Counter(), collections.Counter()
    prim = {}
    for q, cat, rname, s, e, as_, js in alns:
        if cat == "P":
            prim[q] = (rname, s, e)
    for q in all_reads:
        t = q.split("|")[0]
        cid = tx_copy.get(t)
        if cid:
            n_reads[cid] += 1
            p = prim.get(q)
            s0, e0 = span.get(cid, (0, 0))
            if p and p[0] == chrom and p[1] < e0 and p[2] > s0:
                back[cid] += 1
    chain_rows = []
    for t in targets:
        for ch in sorted(chains.get(t["cid"], [])):
            chain_rows.append(dict(cid=t["cid"], family=t["family"], n_introns=len(ch), chain=",".join(f"{d}-{x}" for d, x in ch),
                                   reads_P1=len(carry["P1"][(t["cid"], ch)]), reads_P2=len(carry["P2"][(t["cid"], ch)]), reads_P3=len(carry["P3"][(t["cid"], ch)]),
                                   observable=int(len(carry["P2"][(t["cid"], ch)]) >= 3)))
    with open(a.out + ".chains.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=["cid", "family", "n_introns", "chain", "reads_P1", "reads_P2", "reads_P3", "observable"], delimiter="\t"); w.writeheader(); w.writerows(chain_rows)
    e0rows, rrows = [], []
    for t in targets:
        cid = t["cid"]
        ch = sorted(chains.get(cid, []))
        cnt = {p: sum(1 for c in ch if len(carry[p][(cid, c)]) >= 3) for p in carry}
        observable = cnt["P2"] >= 1
        in_r = t["stratum"] == "R_in" and observable
        e0rows.append(dict(cid=cid, name=t["name"], family=t["family"], heldout=t["heldout"], stratum_in=t["stratum"], is_E=t["entangled"], is_C=t["mono"], is_X=t["rt_image"], n_chains=len(ch), reads_sim=n_reads[cid], reads_back_P1=back[cid],
                           reads_back_share=round(back[cid] / n_reads[cid], 3) if n_reads[cid] else "", observable_chains_P1=cnt["P1"], observable_chains_P2=cnt["P2"],
                           observable_chains_P3=cnt["P3"], aligner_limited=int(t["stratum"] == "R_in" and not observable), in_R=int(in_r)))
        rrows.append(dict(cid=cid, family=t["family"], heldout=t["heldout"], stratum_in=t["stratum"], in_R=int(in_r),
                          reason=("" if in_r else (t["stratum"] if t["stratum"] != "R_in" else "aligner-limited"))))
    for p, rows in (("E0", e0rows), ("R", rrows)):
        with open(f"{a.out}.{p}.tsv", "w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
    h = hashlib.sha1("\n".join(f"{r['cid']}\t{r['in_R']}" for r in sorted(rrows, key=lambda r: r["cid"])).encode()).hexdigest()
    open(a.out + ".R.sha1", "w").write(h + "\n")
    # ---- G2: single-copy genes
    target_genes = {t["gene"] for t in targets}
    multi = set()
    for q, cat, rname, s, e, as_, js in alns:
        if cat == "S" and as_ >= UNIQ * best[q]:
            multi.add(tx_gene.get(q.split("|")[0]))
    single = [g for g in gene_span if g not in target_genes and g not in named and g not in multi]
    ok = tot = 0
    for q in all_reads:
        g = tx_gene.get(q.split("|")[0])
        if g in single:
            tot += 1
            p = prim.get(q)
            s0, e0 = gene_span[g]
            ok += int(bool(p) and p[0] == chrom and p[1] < e0 and p[2] > s0)
    with open(a.out + ".single_copy.tsv", "w") as fh:
        fh.write("gene\n"); [fh.write(g + "\n") for g in sorted(single)]
    n_main = sum(1 for t in targets if not int(t["heldout"]))
    g2 = dict(single_copy_genes=len(single), reads=tot, primary_on_source=ok, share=(ok / tot if tot else None), required=0.99, ok=bool(tot and ok / tot >= 0.99),
              R_sha1=h, R_size=sum(r["in_R"] for r in rrows if not int(r["heldout"])), main_copies=n_main,
              R_in_size=sum(1 for t in targets if t["stratum"] == "R_in" and not int(t["heldout"])))
    json.dump(g2, open(a.out + ".G2.json", "w"), indent=1)
    print(json.dumps(g2))


if __name__ == "__main__":
    main()
