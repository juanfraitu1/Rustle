#!/usr/bin/env python3
"""Does padding the exon-sum with genomic flanks make it more informative — or just more promiscuous?

The exon-sum discards introns and everything outside the transcript. The genomic-span tier already shows
that is lossy: on Soto it supplies 317 edges that no other tier finds, including 22 in one NPIP family.
The open question is whether the missing signal is INTRONIC (inside the locus) or FLANKING (outside it),
because those imply different fixes.

The experiment must measure both directions at once. Padding lengthens every sequence, so it can only add
alignable material — recall alone would show monotone "improvement" and prove nothing. Segmental duplications
share their FLANKS by construction, so padding is exactly the operation that would connect two loci because
they sit in the same duplicated block rather than because they are the same gene. That is the NPIP-vs-FAM72
distinction, and it would appear as CROSS-FAMILY edges.

So every arm runs one all-vs-all over ALL selected families pooled, and reports:
  WITHIN  edges joining two members of the same Soto family      (signal)
  CROSS   edges joining members of DIFFERENT Soto families       (over-merge)
A padding level only earns its place if WITHIN rises while CROSS stays flat.
"""
import argparse, subprocess, sys, tempfile, os
from collections import defaultdict

REF = "/mnt/c/Users/jfris/Desktop/Reference"
FASTA = os.environ.get("SOTO_FASTA", f"{REF}/chm13v2.0.fa")
GFF = os.environ.get("SOTO_GFF", f"{REF}/chm13v2.0_RefSeq_full.gff.gz")
BED = os.environ.get("SOTO_BED", "bench/soto/80_fams.chr.bed")


def members(fams):
    out = defaultdict(list)
    for ln in open(BED):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 4 or "|" not in f[3]:
            continue
        name, fam = f[3].split("|")[0], f[3].split("|")[1]
        if fam in fams:
            out[fam].append((name, f[0], int(f[1]), int(f[2])))
    return out


def exons(chrom, s, e, gene_sub=None):
    raw = subprocess.run(["tabix", GFF, f"{chrom}:{s+1}-{e}"], capture_output=True, text=True).stdout
    iv = []
    for ln in raw.splitlines():
        t = ln.split("\t")
        if len(t) < 9 or t[2] != "exon":
            continue
        if gene_sub:
            g = [x[5:] for x in t[8].split(";") if x.startswith("gene=")]
            if not g or gene_sub not in g[0]:
                continue
        a, b = max(int(t[3]) - 1, s), min(int(t[4]), e)
        if b > a:
            iv.append((a, b))
    iv.sort()
    mg = []
    for a, b in iv:
        if mg and a <= mg[-1][1]:
            mg[-1][1] = max(mg[-1][1], b)
        else:
            mg.append([a, b])
    return mg


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


def build(mem, flank, mode):
    """mode: 'exon' = exon-sum (+flank), 'span' = whole genomic span (the existing tier, for reference)."""
    seqs = {}
    for fam, rows in mem.items():
        for name, c, s, e in rows:
            if mode == "span":
                iv = [(max(0, s - flank), e + flank)]
            else:
                iv = exons(c, s, e)
                if not iv:
                    continue
                if flank:
                    # A leading and a trailing genomic block, OUTSIDE the transcript. Not more exons —
                    # the question is whether the sequence around the locus carries family signal.
                    iv = [(max(0, s - flank), s)] + iv + [(e, e + flank)]
            sq = fetch(c, iv)
            if sq:
                seqs[f"{fam}|{name}"] = sq
    return seqs


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
    ap.add_argument("--families", nargs="+", required=True)
    ap.add_argument("--flanks", nargs="+", type=int, default=[0, 500, 1000, 2000, 5000])
    ap.add_argument("--min-identity", type=float, default=0.80)
    ap.add_argument("--min-coverage", type=float, default=0.50)
    ap.add_argument("--mode", default="exon", choices=["exon", "span"])
    a = ap.parse_args()

    mem = members(set(a.families))
    if not mem:
        sys.exit("no members found")
    npairs_within = sum(len(v) * (len(v) - 1) // 2 for v in mem.values())
    tot = sum(len(v) for v in mem.values())
    npairs_cross = tot * (tot - 1) // 2 - npairs_within
    print(f"{len(mem)} families, {tot} members  |  within-family pairs {npairs_within}, "
          f"cross-family pairs {npairs_cross}\n")
    print(f"{'flank bp':>9}{'mean len':>10}{'WITHIN':>9}{'of':>6}{'%':>7}{'CROSS':>8}{'of':>7}{'%':>7}"
          f"   families fully connected")
    for fl in a.flanks:
        seqs = build(mem, fl, a.mode)
        ed = edges(seqs, a.min_identity, a.min_coverage)
        fam_of = lambda k: k.split("|")[0]
        within = sum(1 for x, y in ed if fam_of(x) == fam_of(y))
        cross = len(ed) - within
        # how many families are a single connected component under these edges
        full = 0
        for fam, rows in mem.items():
            nodes = [f"{fam}|{n}" for n, _, _, _ in rows if f"{fam}|{n}" in seqs]
            p = {v: v for v in nodes}
            def find(x):
                while p[x] != x:
                    p[x] = p[p[x]]; x = p[x]
                return x
            for x, y in ed:
                if x in p and y in p:
                    p[find(x)] = find(y)
            if nodes and len({find(v) for v in nodes}) == 1:
                full += 1
        ml = sum(len(s) for s in seqs.values()) // max(len(seqs), 1)
        print(f"{fl:>9,}{ml:>10,}{within:>9}{npairs_within:>6}{100*within/max(npairs_within,1):>6.0f}%"
              f"{cross:>8}{npairs_cross:>7}{100*cross/max(npairs_cross,1):>6.1f}%"
              f"   {full}/{len(mem)}")


if __name__ == "__main__":
    main()
