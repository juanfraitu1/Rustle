#!/usr/bin/env python3
"""Falsification test: does the read-derived exon UNION substrate produce PURER families than the single chain?

The union substrate tripled within-family recall on chr16 with zero cross-family edges. That is exactly the
shape of a result that has fooled this project before: a seeding relaxation once looked like a +2-member gain
while it was fusing PCDHB into the histone cluster, because the evaluation was partition-blind. So this script
is built to FALSIFY, not to confirm — it pools every family into ONE all-vs-all and reports purity first.

Recall cannot be the headline here. Lengthening every sequence can only add alignable material, so within-family
edges will rise almost regardless. The question is whether the edges that appear STAY INSIDE their family.

Three purity measures, because each can fail independently:
  CROSS        edges joining two different truth families                (direct violation)
  IMPURE       connected components containing >1 truth family           (what a user actually receives)
  ADVERSARIAL  edges between a named pair known to be confusable         (PCDHB vs histones)

Truth is a BED of `chrom start end FAMILY`. Any family label works, so the same test runs on gorilla known
gene families (clean, well-annotated protein-coding) and on Soto (which mixes genes with pseudogenes, lncRNAs
and unannotated fragments — the case where an over-inclusive substrate should break first).
"""
import argparse, os, re, subprocess, sys, tempfile
from collections import defaultdict, Counter


def load_bed(path):
    rows = []
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        fam = f[3].split("|")[-1] if "|" in f[3] else f[3]
        name = f[3].split("|")[0] if "|" in f[3] else f"{f[0]}:{f[1]}"
        rows.append((fam, name, f[0], int(f[1]), int(f[2])))
    return rows


def read_chains(bam, chrom, s, e, contained=True):
    txt = subprocess.run(["samtools", "view", "-F", "2308", bam, f"{chrom}:{s+1}-{e}"],
                         capture_output=True, text=True).stdout
    chains = defaultdict(list)
    for ln in txt.splitlines():
        f = ln.split("\t")
        if len(f) < 6:
            continue
        pos, cig = int(f[3]) - 1, f[5]
        ex, p = [], pos
        for n, o in re.findall(r"(\d+)([MIDNSHP=X])", cig):
            n = int(n)
            if o in "M=X":
                if ex and ex[-1][1] == p:
                    ex[-1][1] = p + n
                else:
                    ex.append([p, p + n])
                p += n
            elif o == "D":
                if ex:
                    ex[-1][1] = p + n
                p += n
            elif o == "N":
                p += n
        if not ex:
            continue
        if contained and (ex[0][0] < s or ex[-1][1] > e):
            continue
        chains[tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))].append(ex)
    return chains


def merge(iv):
    iv = sorted(tuple(x) for x in iv)
    mg = []
    for a, b in iv:
        if mg and a <= mg[-1][1]:
            mg[-1][1] = max(mg[-1][1], b)
        else:
            mg.append([a, b])
    return mg


def substrate(chains, kind, min_reads):
    if not chains:
        return []
    if kind == "single":
        return merge(max(chains.items(), key=lambda kv: len(kv[1]))[1][0])
    if kind == "union":
        return merge([x for reads in chains.values() if len(reads) >= min_reads for ex in reads for x in ex])
    return merge([x for reads in chains.values() for ex in reads for x in ex])


def fetch(fasta, chrom, iv):
    if not iv:
        return ""
    regs = [f"{chrom}:{a+1}-{b}" for a, b in iv]
    out = []
    for i in range(0, len(regs), 200):
        r = subprocess.run(["samtools", "faidx", fasta] + regs[i:i+200], capture_output=True, text=True).stdout
        for blk in r.split(">")[1:]:
            out.append("".join(blk.splitlines()[1:]))
    return "".join(out).upper()


def edges(seqs, min_id, min_cov, threads=4):
    fh = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
    for k, v in seqs.items():
        fh.write(f">{k}\n{v}\n")
    fh.close()
    out = set()
    for pre in (["-x", "asm20"], ["-k", "11", "-w", "5"]):
        paf = subprocess.run(["minimap2", "-c", "-X", "--no-long-join", "-t", str(threads)] + pre
                             + [fh.name, fh.name], capture_output=True, text=True).stdout
        for ln in paf.splitlines():
            f = ln.split("\t")
            if len(f) < 12:
                continue
            q, ql, qs, qe, t, tl = f[0], int(f[1]), int(f[2]), int(f[3]), f[5], int(f[6])
            if q == t:
                continue
            de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
            if de is None:
                continue
            if (1 - de) >= min_id and (qe - qs) / max(min(ql, tl), 1) >= min_cov:
                out.add((min(q, t), max(q, t)))
    os.unlink(fh.name)
    return out


def components(nodes, ed):
    p = {v: v for v in nodes}

    def find(x):
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    for a, b in ed:
        if a in p and b in p:
            p[find(a)] = find(b)
    g = defaultdict(list)
    for v in nodes:
        g[find(v)].append(v)
    return list(g.values())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--label", default="")
    ap.add_argument("--min-chain-reads", type=int, default=3)
    ap.add_argument("--min-identity", type=float, default=0.80)
    ap.add_argument("--min-coverage", type=float, default=0.50)
    ap.add_argument("--adversarial", nargs="*", default=[],
                    help="family names expected to stay APART, e.g. PCDHB H2BC H4C")
    a = ap.parse_args()

    rows = load_bed(a.bed)
    fams = sorted({r[0] for r in rows})
    print(f"{a.label or a.bed}: {len(fams)} families, {len(rows)} truth loci\n")

    built = {k: {} for k in ("single", "union", "covered")}
    for fam, name, c, s, e in rows:
        ch = read_chains(a.bam, c, s, e)
        if not ch:
            continue
        for kind in built:
            iv = substrate(ch, kind, a.min_chain_reads)
            sq = fetch(a.fasta, c, iv)
            if sq:
                built[kind][f"{fam}|{name}|{c}:{s}"] = sq

    famof = lambda k: k.split("|")[0]
    hdr = (f"{'substrate':<10}{'loci':>6}{'mean bp':>10}{'WITHIN':>8}{'%':>6}"
           f"{'CROSS':>8}{'IMPURE comps':>14}{'pure comps':>12}")
    print(hdr)
    print("-" * len(hdr))
    adv_hits = {}
    for kind in ("single", "union", "covered"):
        seqs = built[kind]
        if not seqs:
            continue
        ed = edges(seqs, a.min_identity, a.min_coverage)
        sizes = Counter(famof(k) for k in seqs)
        wp = sum(v * (v - 1) // 2 for v in sizes.values())
        within = sum(1 for x, y in ed if famof(x) == famof(y))
        cross = len(ed) - within
        comps = components(list(seqs), ed)
        multi = [c for c in comps if len(c) >= 2]
        impure = [c for c in multi if len({famof(v) for v in c}) > 1]
        ml = sum(len(v) for v in seqs.values()) // len(seqs)
        print(f"{kind:<10}{len(seqs):>6}{ml:>10,}{within:>8}{100*within/max(wp,1):>5.0f}%"
              f"{cross:>8}{len(impure):>8}/{len(multi):<5}{len(multi)-len(impure):>12}")
        if a.adversarial:
            hits = [(x, y) for x, y in ed
                    if famof(x) != famof(y) and famof(x) in a.adversarial and famof(y) in a.adversarial]
            adv_hits[kind] = hits
        if impure:
            for c in impure[:4]:
                print(f"           IMPURE: {sorted(Counter(famof(v) for v in c).items())}")

    if a.adversarial:
        print(f"\nADVERSARIAL pair check ({' vs '.join(a.adversarial)}):")
        for kind, hits in adv_hits.items():
            print(f"  {kind:<10}{len(hits)} edge(s)" + (f"  {hits[:3]}" if hits else "  -- stayed apart"))


if __name__ == "__main__":
    main()
