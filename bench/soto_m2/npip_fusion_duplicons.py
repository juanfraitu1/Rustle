#!/usr/bin/env python3
"""Pre-registered test KEY=npipfusion (docs/PREREG_npip_fusion_duplicon_2026-09-30.md): are NPIP fusion transcripts
duplicon-boundary crossings inside co-duplicated blocks?

Units are distinct splice junctions (>= 2 reads) in long-read alignments at NPIP genes. A switch unit joins an NPIP exon
block (overlapping an NPIP exon) to an outside block (overlapping no NPIP gene); an internal unit joins two NPIP exon blocks.
H1: switches cross a duplicon boundary (flanks share no duplicon) more often than internal units (one-sided Fisher).
H2: most switch partners sit on a duplicon co-duplicated with the NPIP core (one-sided binomial against 1/2).
The held-out library decides; the development library is reported beside it.

    python3 bench/soto_m2/npip_fusion_duplicons.py --dev A119b.t2t.bam --heldout human_testis.t2t.bam \
        --refseq chm13v2.0_RefSeq_full.gff.gz --cat families_cn.json --duplicons dupmasker_colors.bed \
        --sd98 sd98_v1.bed --sd98-fa sd98_regions.fa --genome chm13v2.0.fa --out-tsv partners.tsv"""
import argparse
import bisect
import collections
import gzip
import json
import math
import random
import re
import subprocess

import pysam

CHROMS = ("chr16", "chr18")
SHIFT = {"chr16": -5}  # CHM13 v1.0 -> v2.0 (Amendment 1: a 5 bp telomeric indel on chr16; chr18 unchanged)
ALPHA = 0.01


def merge(ivs):
    out = []
    for a, b in sorted(ivs):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


class Merged:
    """Merged intervals with an overlap query."""

    def __init__(self, ivs):
        self.iv = merge(ivs)
        self.st = [a for a, _ in self.iv]

    def hits(self, a, b):
        k = bisect.bisect_right(self.st, a) - 1
        k = max(k, 0)
        while k < len(self.iv) and self.iv[k][0] < b:
            if self.iv[k][1] > a:
                return True
            k += 1
        return False


class Segs:
    """Possibly overlapping labelled segments with an overlap query returning {label: bases}."""

    def __init__(self, segs):
        self.s = sorted(segs)
        self.st = [x[0] for x in self.s]
        self.maxlen = max((e - s for s, e, _ in self.s), default=0)

    def comp(self, a, b):
        out = collections.Counter()
        k = bisect.bisect_left(self.st, a - self.maxlen)
        while k < len(self.s) and self.s[k][0] < b:
            s, e, lab = self.s[k]
            o = min(b, e) - max(a, s)
            if o > 0:
                out[lab] += o
            k += 1
        return out


def dominant(comp):
    return min(comp.items(), key=lambda kv: (-kv[1], kv[0]))[0] if comp else None


def gff_attrs(col):
    return dict(x.split("=", 1) for x in col.rstrip().split(";") if "=" in x)


def load_refseq(path, chroms):
    """Genes (name, chrom, start, end, exons) on the given chromosomes, exons pooled over the gene's transcripts."""
    genes, tx_gene, exons = {}, {}, collections.defaultdict(list)
    for c in chroms:
        for line in subprocess.run(["tabix", path, c], capture_output=True, text=True, check=True).stdout.splitlines():
            f = line.split("\t")
            a = gff_attrs(f[8])
            i = a.get("ID", "")
            if i.startswith("gene-"):
                genes[i] = (a.get("Name", i[5:]), c, int(f[3]) - 1, int(f[4]))
            elif f[2] == "exon":
                exons[a.get("Parent", "")].append((int(f[3]) - 1, int(f[4])))
            elif a.get("Parent", "").startswith("gene-") and i:
                tx_gene[i] = a["Parent"]
    by_gene = collections.defaultdict(list)
    for p, ex in exons.items():
        g = p if p.startswith("gene-") else tx_gene.get(p)
        if g in genes:
            by_gene[g].extend(ex)
    return [(n, c, s, e, merge(by_gene[g])) for g, (n, c, s, e) in genes.items()], exons, tx_gene, genes


def is_npip(name):
    return name.upper().startswith("NPIP") and "-" not in name


def read_blocks(r):
    pos, start, out = r.reference_start, r.reference_start, []
    for op, n in r.cigartuples:
        if op in (0, 2, 7, 8):
            pos += n
        elif op == 3:
            if pos > start:
                out.append((start, pos))
            pos += n
            start = pos
    if pos > start:
        out.append((start, pos))
    return out


def fisher_greater(a, b, c, d):
    """One-sided Fisher exact p for [[a, b], [c, d]], alternative: row 1 has the larger share of column 1."""
    n1, k, n = a + b, a + c, a + b + c + d
    den = math.comb(n, n1)
    return sum(math.comb(k, x) * math.comb(n - k, n1 - x) for x in range(a, min(n1, k) + 1)) / den


def binom_greater(x, n):
    return sum(math.comb(n, k) for k in range(x, n + 1)) / 2 ** n if n else 1.0


def junctions(bam_path, bn, mapq_min):
    """Distinct junctions -> support and flank blocks, from primary alignments overlapping the NPIP gene spans."""
    bam = pysam.AlignmentFile(bam_path)
    seen, J = set(), collections.defaultdict(lambda: [0, collections.Counter(), collections.Counter()])
    for c in CHROMS:
        for a, b in bn[c].iv:
            for r in bam.fetch(c, a, b):
                if r.flag & 2308 or r.mapping_quality < mapq_min:
                    continue
                key = (r.query_name, r.reference_start, r.flag)
                if key in seen:
                    continue
                seen.add(key)
                bl = read_blocks(r)
                for x, y in zip(bl, bl[1:]):
                    j = J[(c, x[1], y[0])]
                    j[0] += 1
                    j[1][x] += 1
                    j[2][y] += 1
    return J, len(seen)


def modal(cnt):
    return min(cnt.items(), key=lambda kv: (-kv[1], -(kv[0][1] - kv[0][0]), kv[0]))[0]


def analyse(bam_path, mapq_min, ctx):
    tn, bn, dups, codup, sd98 = ctx["tn"], ctx["bn"], ctx["dups"], ctx["codup"], ctx["sd98"]
    J, nreads = junctions(bam_path, bn, mapq_min)
    cls = lambda c, blk: "N" if tn[c].hits(*blk) else ("O" if not bn[c].hits(*blk) else "I")
    units = []
    for (c, s, e), (n, left, right) in J.items():
        if n < 2:
            continue
        L, R = modal(left), modal(right)
        cl, cr = cls(c, L), cls(c, R)
        kind = "switch" if {cl, cr} == {"N", "O"} else "internal" if cl == cr == "N" else None
        if not kind:
            continue
        dl, dr = dups[c].comp(*L), dups[c].comp(*R)
        u = dict(c=c, s=s, e=e, n=n, L=L, R=R, kind=kind, b=int(not (set(dl) & set(dr))))
        if kind == "switch":
            o_blk, o_comp = (L, dl) if cl == "O" else (R, dr)
            u.update(o_side="left" if cl == "O" else "right", O=o_blk, N=R if cl == "O" else L,
                     pdup=dominant(o_comp), cdup=int(dominant(o_comp) in codup), single=int(not sd98[c].hits(*o_blk)))
        units.append(u)
    sw = [u for u in units if u["kind"] == "switch"]
    it = [u for u in units if u["kind"] == "internal"]
    a, b = sum(u["b"] for u in sw), len(sw) - sum(u["b"] for u in sw)
    c_, d = sum(u["b"] for u in it), len(it) - sum(u["b"] for u in it)
    p1 = fisher_greater(a, b, c_, d) if sw and it else 1.0
    h1 = p1 < ALPHA and len(sw) and len(it) and a / len(sw) > c_ / len(it)
    x = sum(u["cdup"] for u in sw)
    p2 = binom_greater(x, len(sw))
    return dict(reads=nreads, units=units, sw=sw, it=it, a=a, c=c_, p1=p1, h1=bool(h1), x=x, p2=p2, h2=p2 < ALPHA)


def null_h2(sw, ctx, reps=10000, seed=20260930):
    rng = random.Random(seed)
    bn, dups, codup = ctx["bn"], ctx["dups"], ctx["codup"]
    per = []
    for u in sw:
        g, ln = u["e"] - u["s"], u["O"][1] - u["O"][0]
        vals = []
        for _ in range(reps):
            v = u["cdup"]
            for _ in range(100):
                dist = rng.uniform(g / 2, 2 * g)
                if u["o_side"] == "right":
                    a = int(u["s"] + dist)
                    blk = (a, a + ln)
                else:
                    b = int(u["e"] - dist)
                    blk = (b - ln, b)
                if blk[0] >= 0 and not bn[u["c"]].hits(*blk):
                    v = int(dominant(dups[u["c"]].comp(*blk)) in codup)
                    break
            vals.append(v)
        per.append(vals)
    if not per:
        return None, None
    means = [sum(p[r] for p in per) / len(per) for r in range(reps)]
    obs = sum(u["cdup"] for u in sw) / len(sw)
    return sum(means) / reps, (1 + sum(m >= obs for m in means)) / (reps + 1)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("dev", "heldout", "refseq", "cat", "duplicons", "sd98", "sd98_fa", "genome", "out_tsv"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args(argv)

    # coordinates (Amendment 1): after the chr16 shift, every SD98 region on chr16/chr18 that does not touch a chromosome end
    # must have the same sequence in v1.0 and v2.0
    fa1, fa2 = pysam.FastaFile(a.sd98_fa), pysam.FastaFile(a.genome)
    nchk, nend = 0, 0
    for line in open(a.sd98):
        c, s, e = line.split()[:3]
        if c in CHROMS:
            s, e = int(s), int(e)
            d_ = SHIFT.get(c, 0)
            if s + d_ <= 0 or e + d_ >= fa2.get_reference_length(c):
                nend += 1
                continue
            assert fa1.fetch(reference=f"{c}:{s + 1}-{e}").upper() == fa2.fetch(c, s + d_, e + d_).upper(), \
                f"v1.0 != v2.0 at {c}:{s}-{e} after shift {d_}"
            nchk += 1
    sh = lambda c, x: x + SHIFT.get(c, 0)

    # NPIP genes from RefSeq and CAT
    rs_genes, rs_exons, rs_tx, rs_gid = load_refseq(a.refseq, CHROMS)
    cat = json.load(open(a.cat))["genes"]
    npip = [("RefSeq", n, c, s, e, ex) for n, c, s, e, ex in rs_genes if is_npip(n) and ex]
    npip += [("CAT", g["n"], g["c"], sh(g["c"], g["s"]), sh(g["c"], g["e"]), merge([(sh(g["c"], x), sh(g["c"], y)) for x, y in g["x"]]))
             for g in cat if g["c"] in CHROMS and is_npip(g["n"])]
    tn = {c: Merged([iv for *_, cc, s, e, ex in npip if cc == c for iv in ex]) for c in CHROMS}
    bn = {c: Merged([(s, e) for *_, cc, s, e, ex in npip if cc == c]) for c in CHROMS}

    # duplicons genome-wide, SD98 regions, core and co-duplicated sets
    seg_by = collections.defaultdict(list)
    for line in open(a.duplicons):
        f = line.split("\t")
        seg_by[f[0]].append((int(f[1]), int(f[2]), f[3]))
    dups_v1 = {c: Segs(v) for c, v in seg_by.items()}
    dups = {c: Segs([(sh(c, s), sh(c, e), i) for s, e, i in seg_by[c]]) for c in CHROMS}
    sd_by = collections.defaultdict(list)
    for line in open(a.sd98):
        c, s, e = line.split()[:3]
        sd_by[c].append((int(s), int(e)))
    sd98 = {c: Merged([(sh(c, s), sh(c, e)) for s, e in sd_by[c]]) for c in CHROMS}
    per_gene = []
    for *_, c, s, e, ex in npip:
        ids = set()
        for iv in ex:
            ids |= set(dups[c].comp(*iv))
        per_gene.append(ids)
    cnt = collections.Counter(i for ids in per_gene for i in ids)
    core = {i for i, n in cnt.items() if n >= len(per_gene) / 2}
    region_ids = []
    for c, ivs in sd_by.items():
        if c not in dups_v1:
            continue
        for s, e in ivs:
            region_ids.append(set(dups_v1[c].comp(s, e)))
    holding = collections.Counter(i for ids in region_ids if ids & core for i in ids)
    codup = {i for i, n in holding.items() if n >= 2}
    ctx = dict(tn=tn, bn=bn, dups=dups, codup=codup, sd98=sd98)

    print(f"# KEY=npipfusion\ncoordinate check: {nchk} SD98 regions on chr16/chr18 identical in v1.0 and v2.0 after the chr16 shift "
          f"({nend} touching a chromosome end not checked)")
    print(f"NPIP gene records: {sum(1 for x in npip if x[0] == 'RefSeq')} RefSeq + {sum(1 for x in npip if x[0] == 'CAT')} CAT; "
          f"core duplicons ({len(core)}): {', '.join(sorted(core))}; co-duplicated duplicons: {len(codup)}")

    res = {}
    for lab, path in (("dev", a.dev), ("heldout", a.heldout)):
        for mq in (0, 1):
            res[(lab, mq)] = analyse(path, mq, ctx)

    def verdict(r):
        if len(r["sw"]) < 10:
            return "UNDERPOWERED"
        return "EXPLAINED" if r["h1"] and r["h2"] else "BOUNDARY ONLY" if r["h1"] else "NOT EXPLAINED"

    print("\n| library | reads | switch units | internal units | boundary: switch / internal | H1 p | co-duplicated partners | H2 p | verdict |")
    print("|---|---|---|---|---|---|---|---|---|")
    for (lab, mq), r in res.items():
        ns, ni = len(r["sw"]), len(r["it"])
        bs = f"{r['a'] / ns:.3f}" if ns else "-"
        bi = f"{r['c'] / ni:.3f}" if ni else "-"
        cd = f"{r['x']}/{ns} = {r['x'] / ns:.3f}" if ns else "-"
        print(f"| {lab}{' MAPQ>=1' if mq else ''} | {r['reads']:,} | {ns} | {ni} | {bs} / {bi} | {r['p1']:.2g} | {cd} | {r['p2']:.2g} | "
              f"{verdict(r) if not mq else '(secondary)'} |")
    vd, vh = verdict(res[("dev", 0)]), verdict(res[("heldout", 0)])
    print(f"\nVERDICT (held-out decides): {vh}" + (f"  [development: {vd} -> SPLIT]" if vd != vh else f"  [development agrees: {vd}]"))

    # secondary
    rs_names = collections.defaultdict(list)
    for n, c, s, e, ex in rs_genes:
        rs_names[c].append((s, e, n))
    cat_names = collections.defaultdict(list)
    for g in cat:
        if g["c"] in CHROMS:
            cat_names[g["c"]].append((sh(g["c"], g["s"]), sh(g["c"], g["e"]), g["n"]))
    names_at = lambda c, blk, src: sorted({n for s, e, n in src[c] if s < blk[1] and blk[0] < e})
    rt = {}
    for gid, (n, c, s, e) in rs_gid.items():
        if "-" in n and "NPIP" in n.upper():
            for t, g in rs_tx.items():
                if g == gid:
                    ex = sorted(rs_exons[t])
                    rt.setdefault(n, set()).update((c, x[1], y[0]) for x, y in zip(ex, ex[1:]))
    with open(a.out_tsv, "w") as out:
        out.write("library\tchrom\tintron_start\tintron_end\treads\tnpip_genes\tpartner_side\tpartner_block\tpartner_dominant_duplicon"
                  "\tpartner_codup\tpartner_in_sd98\tpartner_refseq\tpartner_cat\tboundary\n")
        for lab in ("dev", "heldout"):
            r = res[(lab, 0)]
            for u in sorted(r["sw"], key=lambda u: (u["c"], u["s"])):
                npg = sorted(set(names_at(u["c"], u["N"], rs_names) + names_at(u["c"], u["N"], cat_names)) & {x[1] for x in npip})
                out.write(f"{lab}\t{u['c']}\t{u['s']}\t{u['e']}\t{u['n']}\t{','.join(npg)}\t{u['o_side']}\t{u['O'][0]}-{u['O'][1]}"
                          f"\t{u['pdup']}\t{u['cdup']}\t{1 - u['single']}\t{','.join(names_at(u['c'], u['O'], rs_names))}"
                          f"\t{','.join(names_at(u['c'], u['O'], cat_names))}\t{u['b']}\n")
    for lab in ("dev", "heldout"):
        r = res[(lab, 0)]
        sw = r["sw"]
        print(f"\n## {lab}: secondary")
        nm, npv = null_h2(sw, ctx)
        if nm is not None:
            print(f"- location-matched null for H2: observed {r['x'] / len(sw):.3f} vs null mean {nm:.3f}, p = {npv:.4g}")
        print(f"- single-copy partners (outside every SD98 region): {sum(u['single'] for u in sw)} of {len(sw)}")
        rec = collections.defaultdict(set)
        sup = collections.Counter()
        for u in sw:
            npg = set(names_at(u["c"], u["N"], rs_names) + names_at(u["c"], u["N"], cat_names)) & {x[1] for x in npip}
            rec[u["pdup"]] |= npg
            sup[u["pdup"]] += u["n"]
        top = sorted(rec, key=lambda d: -sup[d])[:12]
        print("- partners by read support (dominant duplicon: reads, NPIP genes joined): " +
              "; ".join(f"{d}: {sup[d]}, {len(rec[d])}" for d in top))
        sw_keys = {(u["c"], u["s"], u["e"]) for u in sw}
        for n, ints in sorted(rt.items()):
            hit = sorted(k for k in ints if k in sw_keys)
            print(f"- RefSeq read-through {n}: {len(hit)} of its introns are switch units" +
                  (f" ({', '.join(f'{c}:{s:,}-{e:,}' for c, s, e in hit)})" if hit else ""))


if __name__ == "__main__":
    main()
