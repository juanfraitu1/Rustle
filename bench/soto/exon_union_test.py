#!/usr/bin/env python3
"""Would a read-derived exon UNION represent the locus better than the single best chain?

`pick_locus_rep` selects ONE transcript per locus (max reads, then span) and `build_spliced_seq` takes a
SINGLE intron list, so the exon-sum is one chain. In segmental duplications reads shatter into many partial
chains (NOTCH2: 490 reads -> 240 chains, the largest holding 5%), so that one chain is a fragment of the
locus — which is what the bipartite size match measured as median 0.54x truth, 104 truncated vs 12 extended.

This compares three purely read-derived substrates at the same loci, with no annotation used to BUILD them:

  single   exons of the highest-support intron chain            (what ships)
  union    union of exon intervals over all chains with >= K reads   (candidate)
  covered  union over EVERY primary read, no support floor      (upper bound; noise included)

and scores each two ways — length recovered against the truth window, and the edges it produces. Length
alone would favour `covered` trivially, so the edge count is the arbiter: a substrate that lengthens without
producing edges has added noise, not signal.
"""
import argparse, os, re, subprocess, sys, tempfile
from collections import defaultdict

REF = "/mnt/c/Users/jfris/Desktop/Reference"
FASTA = os.environ.get("SOTO_FASTA", f"{REF}/chm13v2.0.fa")
BED = os.environ.get("SOTO_BED", "bench/soto/80_fams.chr.bed")


def truth(fam):
    out = []
    for ln in open(BED):
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 4 and f[3].endswith("|" + fam):
            out.append((f[3].split("|")[0], f[0], int(f[1]), int(f[2])))
    return sorted(out, key=lambda x: (x[1], x[2]))


def read_chains(bam, chrom, s, e):
    """(intron_chain -> [exon lists]) over PRIMARY alignments contained in the window."""
    txt = subprocess.run(["samtools", "view", "-F", "2308", bam, f"{chrom}:{s+1}-{e}"],
                         capture_output=True, text=True).stdout
    chains = defaultdict(list)
    for ln in txt.splitlines():
        f = ln.split("\t")
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
        # keep reads whose footprint lies inside the truth window (avoid dragging in neighbours)
        if ex[0][0] < s or ex[-1][1] > e:
            continue
        chain = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
        chains[chain].append(ex)
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


def substrates(chains, min_reads):
    if not chains:
        return {}
    best = max(chains.items(), key=lambda kv: len(kv[1]))
    single = merge(best[1][0])
    union = merge([x for ch, reads in chains.items() if len(reads) >= min_reads for ex in reads for x in ex])
    covered = merge([x for reads in chains.values() for ex in reads for x in ex])
    return {"single": single, "union": union, "covered": covered}


def fetch(chrom, iv):
    if not iv:
        return ""
    regs = [f"{chrom}:{a+1}-{b}" for a, b in iv]
    out = []
    for i in range(0, len(regs), 200):
        r = subprocess.run(["samtools", "faidx", FASTA] + regs[i:i+200], capture_output=True, text=True).stdout
        for blk in r.split(">")[1:]:
            out.append("".join(blk.splitlines()[1:]))
    return "".join(out).upper()


def edges(seqs, min_id, min_cov):
    fh = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
    for k, v in seqs.items():
        fh.write(f">{k}\n{v}\n")
    fh.close()
    out = set()
    for pre in (["-x", "asm20"], ["-k", "11", "-w", "5"]):
        paf = subprocess.run(["minimap2", "-c", "-X", "--no-long-join", "-t", "4"] + pre + [fh.name, fh.name],
                             capture_output=True, text=True).stdout
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--family", default="ID_154")
    ap.add_argument("--min-chain-reads", type=int, default=3)
    ap.add_argument("--min-identity", type=float, default=0.80)
    ap.add_argument("--min-coverage", type=float, default=0.50)
    a = ap.parse_args()

    tr = truth(a.family)
    print(f"{a.family}: {len(tr)} truth members\n")
    per = {k: {} for k in ("single", "union", "covered")}
    stats = []
    for name, c, s, e in tr:
        ch = read_chains(a.bam, c, s, e)
        if not ch:
            continue
        sub = substrates(ch, a.min_chain_reads)
        nreads = sum(len(v) for v in ch.values())
        row = [name, len(ch), nreads, e - s]
        for k, iv in sub.items():
            L = sum(b - a2 for a2, b in iv)
            row.append(L)
            sq = fetch(c, iv)
            if sq:
                per[k][name] = sq
        stats.append(row)

    print(f"{'member':<10}{'chains':>7}{'reads':>7}{'window':>9}{'single':>9}{'union':>9}{'covered':>9}"
          f"{'union/single':>13}")
    for r in stats:
        ratio = r[5] / r[4] if r[4] else 0
        print(f"{r[0]:<10}{r[1]:>7}{r[2]:>7}{r[3]:>9,}{r[4]:>9,}{r[5]:>9,}{r[6]:>9,}{ratio:>13.2f}")

    print(f"\n{'substrate':<10}{'loci':>6}{'mean bp':>10}{'edges':>8}{'of pairs':>10}{'%':>7}  components")
    for k in ("single", "union", "covered"):
        seqs = per[k]
        if not seqs:
            continue
        ed = edges(seqs, a.min_identity, a.min_coverage)
        n = len(seqs)
        npair = n * (n - 1) // 2
        p = {v: v for v in seqs}
        def find(x):
            while p[x] != x:
                p[x] = p[p[x]]; x = p[x]
            return x
        for x, y in ed:
            p[find(x)] = find(y)
        ncc = len({find(v) for v in seqs})
        ml = sum(len(v) for v in seqs.values()) // n
        print(f"{k:<10}{n:>6}{ml:>10,}{len(ed):>8}{npair:>10}{100*len(ed)/max(npair,1):>6.0f}%  {ncc}")


if __name__ == "__main__":
    main()
