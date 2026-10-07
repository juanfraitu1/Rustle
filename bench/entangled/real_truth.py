#!/usr/bin/env python3
"""Truth tables of docs/PREREG_gorilla_overlap_2026-10-07.md: valid annotated chains, their read support (pools P1 and P2) and the overlap strata, from a RefSeq GFF, the genome and a BAM.

    real_truth.py --gff GFF --fasta FASTA --bam BAM --as-table MOLECULES.tsv --out OUTPREFIX [--contigs C1,C2,...]

Writes OUTPREFIX.chains.tsv (one row per gene and distinct valid chain with n_P1, n_P2), OUTPREFIX.genes.tsv (gene, name, biotype, strand, stratum, junction_sharing, ...) and OUTPREFIX.stats.json.
VALID chain: a transcript with >= 2 exons whose every intron is >= 50 bp with a canonical motif (GT-AG, GC-AG, AT-AC on the transcript strand); chain = ordered intron list (0-based half-open).
P1 = distinct reads with a primary alignment carrying the chain exactly (junctions of N >= 50 bp wholly inside the gene exon-union span equal the chain); P2 = P1 plus secondary alignments with AS >= 0.98 x the molecule's
genome-wide best AS (the as_table). Strata: E_both / E_one / A / N (see the prereg); coordinates 0-based half-open.
"""
import argparse
import bisect
import collections
import csv
import json
import re
import sys

CANON = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
MIN_INTRON = 50
SHARED = 100
GOOD = 0.98
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")
GENE_TYPES = {"gene", "pseudogene"}


def attrs(col):
    return {k: v for k, _, v in (kv.partition("=") for kv in col.split(";") if "=" in kv)}


def parse_gff(src):
    """gene id -> dict(gene, chrom, strand, name, biotype, transcripts {tid: sorted exon list}). Exons attach to the transcript named by their Parent (several parents allowed);
    an exon whose Parent is the gene itself (pseudogenes) gives the gene one transcript of its own; a transcript under a transcript (miRNA under primary_transcript) is followed up to its gene."""
    text = src if "\n" in src else open(src).read()
    genes, parent_of, ftype = {}, {}, {}
    exon_rows = []
    for ln in text.splitlines():
        if not ln or ln.startswith("#"):
            continue
        f = ln.split("\t")
        if len(f) < 9:
            continue
        a = attrs(f[8])
        t = f[2]
        if t in GENE_TYPES and "ID" in a:
            genes[a["ID"]] = dict(gene=a["ID"], chrom=f[0], strand=f[6], name=a.get("Name", a.get("gene", a["ID"])), biotype=a.get("gene_biotype", t), transcripts={})
        elif t == "exon":
            exon_rows.append((f[0], int(f[3]) - 1, int(f[4]), a.get("Parent", "")))
        elif t not in ("CDS", "region", "match", "cDNA_match", "D_loop", "origin_of_replication") and "ID" in a and "Parent" in a:
            parent_of[a["ID"]] = a["Parent"].split(",")[0]
            ftype[a["ID"]] = t

    def gene_of(pid):
        seen = 0
        while pid not in genes and pid in parent_of and seen < 5:
            pid = parent_of[pid]; seen += 1
        return pid if pid in genes else None

    for chrom, s0, e, par in exon_rows:
        for pid in par.split(","):
            g = gene_of(pid)
            if g is None:
                continue
            tid = pid if pid != g else g + ":gene_exons"
            genes[g]["transcripts"].setdefault(tid, []).append((s0, e))
    for g in genes.values():
        for tid in g["transcripts"]:
            g["transcripts"][tid].sort()
    return genes


def merge(exons):
    out = []
    for s, e in sorted(exons):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def valid_chain(chrom, strand, exons, motif):
    """Ordered intron list of the merged exons iff every intron is >= MIN_INTRON bp with a canonical motif; None otherwise (also for < 2 exons)."""
    ex = merge(exons)
    if len(ex) < 2:
        return None
    chain = []
    for i in range(len(ex) - 1):
        i0, i1 = ex[i][1], ex[i + 1][0]
        if i1 - i0 < MIN_INTRON or tuple(motif(chrom, strand, i0, i1)) not in CANON:
            return None
        chain.append((i0, i1))
    return tuple(chain)


def overlap_bp(a, b):
    i = j = tot = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi > lo:
            tot += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def strata(genes, expressed):
    """genes: gid -> dict with chrom, strand, transcripts, chains (tid -> valid chain); expressed: gid -> bool (has an expressed chain).
    Returns gid -> dict(stratum, junction_sharing, partners). E_both: >= SHARED exonic bp on the same strand with another gene that is expressed; E_one: same with no expressed partner;
    A: no same-strand overlap but >= SHARED bp with an opposite-strand gene; N: none. junction_sharing: an annotated valid chain shares an exact junction with a same-strand overlapping partner's."""
    union = {gid: merge([e for ex in g["transcripts"].values() for e in ex]) for gid, g in genes.items()}
    by_chrom = collections.defaultdict(list)
    for gid, u in union.items():
        if u:
            by_chrom[genes[gid]["chrom"]].append((u[0][0], u[-1][1], gid))
    same, opp = collections.defaultdict(set), collections.defaultdict(set)
    for chrom, lst in by_chrom.items():
        lst.sort()
        active = []
        for s, e, gid in lst:
            active = [x for x in active if x[1] > s]
            for s2, e2, g2 in active:
                if overlap_bp(union[gid], union[g2]) >= SHARED:
                    (same if genes[gid]["strand"] == genes[g2]["strand"] else opp)[gid].add(g2)
                    (same if genes[gid]["strand"] == genes[g2]["strand"] else opp)[g2].add(gid)
            active.append((s, e, gid))
    junc = {}
    for gid, g in genes.items():
        js = set()
        for ch in g.get("chains", {}).values():
            js.update(ch)
        junc[gid] = js
    out = {}
    for gid, g in genes.items():
        partners = same.get(gid, set())
        if partners:
            st = "E_both" if any(expressed.get(p) for p in partners) else "E_one"
        elif opp.get(gid):
            st = "A"
        else:
            st = "N"
        out[gid] = dict(stratum=st, junction_sharing=int(any(junc[gid] & junc[p] for p in partners)), partners=sorted(partners))
    return out


class GeneIndex:
    def __init__(self, spans):
        self.by = collections.defaultdict(list)
        for gid, (chrom, strand, s, e) in spans.items():
            self.by[chrom].append((s, e, gid))
        self.starts, self.pmax = {}, {}
        for c, lst in self.by.items():
            lst.sort()
            self.starts[c] = [x[0] for x in lst]
            m, pm = -1, []
            for x in lst:
                m = max(m, x[1]); pm.append(m)
            self.pmax[c] = pm

    def overlapping(self, chrom, s, e):
        lst = self.by.get(chrom)
        if not lst:
            return
        i = bisect.bisect_left(self.starts[chrom], e) - 1
        while i >= 0 and self.pmax[chrom][i] > s:
            if lst[i][1] > s:
                yield lst[i][2]
            i -= 1


def chain_support(alignments, chains, spans, best_as):
    """alignments: iterable of (name, flag, AS, contig, junction list); chains: {(chrom, strand, chain): any}; spans: gid -> (chrom, strand, s, e) exon-union span of the genes carrying those chains;
    best_as: name -> genome-wide best AS. Returns (P1, P2): key -> set of read names. Supplementary and unmapped records are ignored; a secondary alignment joins P2 only if AS >= GOOD x best."""
    gidx = GeneIndex(spans)
    p1, p2 = collections.defaultdict(set), collections.defaultdict(set)
    for name, flag, asv, chrom, juncs in alignments:
        if flag & 2048 or flag & 4 or not juncs:
            continue
        primary = not flag & 256
        if not primary and not (asv is not None and name in best_as and asv >= GOOD * best_as[name]):
            continue
        lo, hi = juncs[0][0], juncs[-1][1]
        done = set()
        for gid in gidx.overlapping(chrom, lo, hi):
            c, strand, s, e = spans[gid]
            key = (c, strand, tuple(j for j in juncs if j[0] >= s and j[1] <= e))
            if key in chains and key not in done:
                done.add(key)
                p2[key].add(name)
                if primary:
                    p1[key].add(name)
    return p1, p2


def junctions_of(rec):
    """Junctions (N >= MIN_INTRON) of a pysam record, 0-based half-open: M = X D advance the reference, N is an intron, I S H do not."""
    pos, out = rec.reference_start, []
    for op, ln in rec.cigartuples or ():
        if op in (0, 7, 8, 2):
            pos += ln
        elif op == 3:
            if ln >= MIN_INTRON:
                out.append((pos, pos + ln))
            pos += ln
    return out


def alignments_of(bam, contig):
    for rec in bam.fetch(contig):
        if rec.is_unmapped or rec.is_supplementary:
            continue
        ct = rec.cigartuples
        if not ct or not any(op == 3 for op, _ in ct):
            continue
        j = junctions_of(rec)
        if j:
            yield (rec.query_name, rec.flag, rec.get_tag("AS") if rec.has_tag("AS") else None, contig, j)


def load_best_as(path):
    d = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            d[f[0]] = int(f[1])
    return d


class Motif:
    def __init__(self, fasta):
        import pysam
        self.fa = pysam.FastaFile(fasta)

    def __call__(self, chrom, strand, i0, i1):
        seq = self.fa.fetch(chrom, i0, i1).upper()
        if strand == "-":
            seq = seq.translate(COMP)[::-1]
        return seq[:2], seq[-2:]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--as-table", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--contigs", default="")
    ap.add_argument("--stage", choices=["annotation", "support", "tables"], default="tables")
    a = ap.parse_args()
    import pickle
    import pysam
    if a.stage == "annotation":
        genes = parse_gff(a.gff)
        mot = Motif(a.fasta)
        for g in genes.values():
            g["chains"] = {}
            for tid, ex in g["transcripts"].items():
                ch = valid_chain(g["chrom"], g["strand"], ex, mot)
                if ch:
                    g["chains"][tid] = ch
        pickle.dump(genes, open(a.out + ".annotation.pkl", "wb"))
        print(f"[real_truth] {len(genes)} genes, {sum(1 for g in genes.values() if g['chains'])} with a valid chain, {sum(len(set(g['chains'].values())) for g in genes.values())} distinct valid chains")
        return
    genes = pickle.load(open(a.out + ".annotation.pkl", "rb"))
    if a.stage == "support":
        best = load_best_as(a.as_table)
        bam = pysam.AlignmentFile(a.bam)
        contigs = [c for c in a.contigs.split(",") if c] or [c for c in bam.references]
        chains, spans = {}, {}
        for gid, g in genes.items():
            if not g["chains"]:
                continue
            u = merge([e for ex in g["transcripts"].values() for e in ex])
            spans[gid] = (g["chrom"], g["strand"], u[0][0], u[-1][1])
            for ch in set(g["chains"].values()):
                chains[(g["chrom"], g["strand"], ch)] = gid
        part = {}
        for c in contigs:
            sp = {gid: v for gid, v in spans.items() if v[0] == c}
            if not sp:
                continue
            ch = {k: v for k, v in chains.items() if k[0] == c}
            p1, p2 = chain_support(alignments_of(bam, c), ch, sp, best)
            part[c] = {k: (len(p1.get(k, ())), len(p2.get(k, ()))) for k in ch}
            print(f"[real_truth] {c}: {len(ch)} chains, {sum(1 for v in part[c].values() if v[0] >= 3)} with >= 3 P1 reads, {sum(1 for v in part[c].values() if v[1] >= 3)} with >= 3 P2 reads", flush=True)
        pickle.dump(part, open(f"{a.out}.support.{'_'.join(contigs)[:60] if a.contigs else 'all'}.pkl", "wb"))
        return
    # tables: merge every support part
    import glob
    sup = {}
    for p in glob.glob(a.out + ".support.*.pkl"):
        for c, d in pickle.load(open(p, "rb")).items():
            sup.setdefault(c, {}).update(d)
    rows_c, rows_g = [], []
    expr = {1: {}, 2: {}}
    for gid, g in genes.items():
        distinct = sorted(set(g["chains"].values()))
        g["chain_support"] = {ch: sup.get(g["chrom"], {}).get((g["chrom"], g["strand"], ch), (0, 0)) for ch in distinct}
        for k in (1, 2):
            expr[k][gid] = any(v[k - 1] >= 3 for v in g["chain_support"].values())
    st = {k: strata(genes, expr[k]) for k in (1, 2)}
    for gid, g in genes.items():
        u = merge([e for ex in g["transcripts"].values() for e in ex])
        for ch, (n1, n2) in g["chain_support"].items():
            rows_c.append(dict(gene=gid, chrom=g["chrom"], strand=g["strand"], chain=",".join(f"{d}-{x}" for d, x in ch), n_P1=n1, n_P2=n2))
        rows_g.append(dict(gene=gid, name=g["name"], biotype=g["biotype"], chrom=g["chrom"], strand=g["strand"], exon_bp=sum(e - s for s, e in u),
                           union=",".join(f"{s}-{e}" for s, e in u), n_valid=len(g["chain_support"]),
                           n_exp_P1=sum(1 for v in g["chain_support"].values() if v[0] >= 3), n_exp_P2=sum(1 for v in g["chain_support"].values() if v[1] >= 3),
                           stratum_P1=st[1][gid]["stratum"], jsharing_P1=st[1][gid]["junction_sharing"], stratum_P2=st[2][gid]["stratum"], jsharing_P2=st[2][gid]["junction_sharing"]))
    for path, rows in ((a.out + ".chains.tsv", rows_c), (a.out + ".genes.tsv", rows_g)):
        with open(path, "w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
    stats = {}
    for k in (1, 2):
        c = collections.Counter(r[f"stratum_P{k}"] for r in rows_g if r[f"n_exp_P{k}"] > 0)
        stats[f"P{k}"] = dict(evaluated_genes=sum(c.values()), by_stratum=dict(c),
                              E_both_junction_sharing=sum(1 for r in rows_g if r[f"n_exp_P{k}"] > 0 and r[f"stratum_P{k}"] == "E_both" and r[f"jsharing_P{k}"]),
                              expressed_chains=sum(1 for r in rows_c if r[f"n_P{k}"] >= 3))
    json.dump(dict(genes=len(genes), valid_chains=len(rows_c), **stats), open(a.out + ".stats.json", "w"), indent=1)
    print(json.dumps(stats))


if __name__ == "__main__":
    main()
