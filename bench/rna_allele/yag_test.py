#!/usr/bin/env python3
"""docs/PREREG_yag_isocon_chain_2026-10-01.md: the reference-absent-copy chain on the Y ampliconic genes (deletion panel, delta_Y).
Subcommands, in order (shell steps between them are the minimap2 / IsoCon runs, see the doc):

  panel      families + copies from the annotation, primaries per copy, the masking rule -> panel.json, masked_tx.fa (CAT transcripts of
             the masked copies, spliced from the unmasked genome), scnets: single-copy gene nets -> sc/<gene>.fa
  mask       masked.fa (= genome with the masked intervals set to N, by .fai offsets)
  reads      <= 500 baseline primaries per copy (seed 1), read orientation -> scored.part{0,1,2}.fa, labels.tsv
  scdelta    sc outputs (sc_outputs.fa) aligned to the unmasked genome (sc_outputs.base.paf) -> delta_Y = p99 of d -> delta_y.json
  dmin       masked_tx.masked.paf -> d_min per masked copy -> dmin.tsv
  net        IsoCon input per family from R.bam (record on a surviving copy, or unmapped), <= 1,000 (seed 1) -> fam/<fam>.fa
  outputs    iso/<fam>/final_candidates.fa -> outputs.fa
  contigs    flag / link (--delta) / merge -> contigs.tsv, contigs_L.fa, merge/
  score      arms R and M (components as loci) -> Y1 / Y2

    yag_test.py panel --w /mnt/linuxdisk/tmp/rna_allele/yag_hsa --sub human
"""
import argparse
import collections
import csv
import glob
import gzip
import json
import os
import random
import shutil
import statistics
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import merge_test  # noqa: E402

SUB = {
    "human": dict(
        genome="/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa",
        idx="/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.splice.mmi",
        bam="/mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam",
        gff="/mnt/linuxdisk/home/juanfraitu/winloci_data/gencode_chm13/chm13v2.0_CAT_Liftoff.slim.gff3.gz",
        chry="chrY",
        families=("TSPY", "RBMY", "DAZ", "CDY", "BPY2", "HSFY", "PRY", "VCY"),
        single_copy=("SRY", "RPS4Y1", "ZFY", "AMELY", "TBL1Y", "PRKY", "USP9Y", "DDX3Y", "UTY", "TMSB4Y", "NLGN4Y", "TXLNGY", "KDM5D",
                     "EIF1AY", "RPS4Y2"),
    ),
}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def gff_chry(cfg):
    """genes: [(name, biotype, s0, e, strand, id)]; tx: id -> (gene_id, biotype, [exons]) for chrY."""
    genes, tx, ex = [], {}, collections.defaultdict(list)
    op = gzip.open if cfg["gff"].endswith(".gz") else open
    for ln in op(cfg["gff"], "rt"):
        if ln[0] == "#":
            continue
        f = ln.rstrip("\n").split("\t")
        if f[0] != cfg["chry"]:
            continue
        at = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        if f[2] == "gene":
            genes.append((at.get("gene_name", at.get("Name", "?")), at.get("gene_biotype", "?"), int(f[3]) - 1, int(f[4]), f[6], at["ID"]))
        elif f[2] == "transcript":
            tx[at["ID"]] = (at["Parent"], at.get("transcript_biotype", "?"), f[6])
        elif f[2] == "exon":
            ex[at["Parent"]].append((int(f[3]) - 1, int(f[4])))
    return genes, tx, ex


def nprim(bam, c, s, e):
    return sum(1 for r in bam.fetch(c, s, e) if not (r.is_secondary or r.is_supplementary or r.is_unmapped))


def spliced(fa, chrom, exons, strand):
    s = "".join(fa.fetch(chrom, a, b) for a, b in sorted(exons))
    return s.translate(COMP)[::-1] if strand == "-" else s


def panel(a):
    cfg = SUB[a.sub]
    genes, tx, ex = gff_chry(cfg)
    bam = pysam.AlignmentFile(cfg["bam"])
    fa = pysam.FastaFile(cfg["genome"])
    out = []
    with open(f"{a.w}/masked_tx.fa", "w") as mt:
        for fam in cfg["families"]:
            recs = sorted((g for g in genes if g[0].startswith(fam)), key=lambda g: g[2])
            copies = []                      # merged intervals: [s, e, names, pc, gene_ids]
            for name, bt, s, e, st, gid in recs:
                if copies and s < copies[-1][1]:
                    c = copies[-1]; c[1] = max(c[1], e); c[2].append(name); c[3] |= bt == "protein_coding"; c[4].append(gid)
                else:
                    copies.append([s, e, [name], bt == "protein_coding", [gid]])
            pc = [c for c in copies if c[3]]
            if len(pc) < 2:
                print(f"{fam}: dropped ({len(pc)} protein-coding copies)"); continue
            for c in copies:
                c.append(nprim(bam, cfg["chry"], c[0], c[1]))
            elig = sorted((c for c in pc if c[5] >= 20), key=lambda c: -c[5])
            masked = [elig[0]] if elig else []
            if len(elig) >= 4:
                masked.append(elig[-1])
            if not masked:
                print(f"{fam}: no protein-coding copy with >= 20 primaries; kept with no mask");
            nm = lambda c: f"{fam}:{c[0]}:{'|'.join(sorted(set(c[2])))}"
            out.append(dict(fam=fam, mask=[[cfg["chry"], c[0], c[1], nm(c)] for c in masked],
                            keep=[[cfg["chry"], c[0], c[1], nm(c)] for c in copies if c not in masked],
                            reads={nm(c): c[5] for c in copies}))
            for c in masked:            # longest protein-coding transcript of the masked copy, spliced from the unmasked genome
                best = None
                for tid, (gid, bt, st) in tx.items():
                    if gid in c[4] and bt == "protein_coding":
                        L = sum(b - s_ for s_, b in ex[tid])
                        if not best or L > best[0]:
                            best = (L, tid, st)
                if best:
                    mt.write(f">{nm(c)}\n{spliced(fa, cfg['chry'], ex[best[1]], best[2])}\n")
            print(f"{fam}: copies {len(copies)} (protein-coding {len(pc)}), masked {[nm(c) + ' reads=' + str(c[5]) for c in masked]}")
    json.dump(out, open(f"{a.w}/panel.json", "w"), indent=0)
    # single-copy gene nets for delta_Y
    os.makedirs(f"{a.w}/sc", exist_ok=True)
    rng = random.Random(1)
    sc = {}
    for name, bt, s, e, st, gid in genes:
        if name in cfg["single_copy"] and bt == "protein_coding":
            sc.setdefault(name, []).append((s, e))
    kept = {}
    for name, ivs in sc.items():
        if len(ivs) != 1:
            continue
        s, e = ivs[0]
        seqs = {}
        for r in bam.fetch(cfg["chry"], s, e):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or r.query_sequence is None:
                continue
            seqs[r.query_name] = r.query_sequence.translate(COMP)[::-1] if r.is_reverse else r.query_sequence
        if len(seqs) < 20:
            continue
        ns = sorted(seqs); rng.shuffle(ns); ns = ns[:1000]
        with open(f"{a.w}/sc/{name}.fa", "w") as o:
            for n in ns:
                o.write(f">{n}\n{seqs[n]}\n")
        kept[name] = (s, e, len(seqs), len(ns))
    json.dump(kept, open(f"{a.w}/sc_genes.json", "w"), indent=0)
    print("single-copy genes with >= 20 primaries:", {k: v[2] for k, v in kept.items()})


def mask(a):
    cfg = SUB[a.sub]
    P = json.load(open(f"{a.w}/panel.json"))
    dst = f"{a.w}/masked.fa"
    shutil.copy(cfg["genome"], dst); shutil.copy(cfg["genome"] + ".fai", dst + ".fai")
    fai = {l.split("\t")[0]: l.split("\t") for l in open(dst + ".fai")}
    tot = 0
    with open(dst, "r+b") as f:
        for p in P:
            for c, s0, e, _ in p["mask"]:
                _, _, off, lb, lw = fai[c][:5]
                off, lb, lw = int(off), int(lb), int(lw)
                pos = s0
                while pos < e:
                    line, col = divmod(pos, lb)
                    k = min(lb - col, e - pos)
                    f.seek(off + line * lw + col); f.write(b"N" * k); pos += k; tot += k
    print("masked bp", tot, "copies", sum(len(p["mask"]) for p in P))


def reads(a):
    cfg = SUB[a.sub]
    P = json.load(open(f"{a.w}/panel.json"))
    bam = pysam.AlignmentFile(cfg["bam"])
    rng = random.Random(1)
    role, seq = {}, {}
    for p in P:
        for lab, (c, s, e, g) in [("D", m) for m in p["mask"]] + [("S", k) for k in p["keep"]]:
            got = {}
            for r in bam.fetch(c, s, e):
                if r.is_secondary or r.is_supplementary or r.is_unmapped or r.query_sequence is None:
                    continue
                got[r.query_name] = r.query_sequence.translate(COMP)[::-1] if r.is_reverse else r.query_sequence
            ns = sorted(got); rng.shuffle(ns)
            for n in ns[:500]:
                if n not in role:
                    role[n] = (p["fam"], lab, g); seq[n] = got[n]
    names = sorted(role)
    k = 3
    for i in range(k):
        with open(f"{a.w}/scored.part{i}.fa", "w") as o:
            for n in names[i * len(names) // k:(i + 1) * len(names) // k]:
                o.write(f">{n}\n{seq[n]}\n")
    with open(f"{a.w}/scored.fa", "w") as o, open(f"{a.w}/labels.tsv", "w") as lab:
        lab.write("read\tfamily\trole\tcopy\n")
        for n in names:
            o.write(f">{n}\n{seq[n]}\n"); f_, r_, g_ = role[n]; lab.write(f"{n}\t{f_}\t{r_}\t{g_}\n")
    c = collections.Counter(r for _, r, _ in role.values())
    print(f"reads: D {c['D']}, S {c['S']}")


def best_hits(path):
    b = {}
    for ln in open(path):
        f = ln.split("\t")
        s = int(f[9]) / max(1, int(f[10])) * (int(f[3]) - int(f[2])) / max(1, int(f[1]))
        if f[0] not in b or s > b[f[0]][0]:
            b[f[0]] = (s, f[5], int(f[7]), int(f[8]), int(f[9]), int(f[1]))
    return b


def scdelta(a):
    kept = json.load(open(f"{a.w}/sc_genes.json"))
    B = best_hits(f"{a.w}/sc_outputs.base.paf")
    cfg = SUB[a.sub]
    ds, per = [], collections.defaultdict(list)
    for o, (s, c, ts, te, m, qlen) in B.items():
        g = o.split("|")[0]
        gs, ge = kept[g][0], kept[g][1]
        if c == cfg["chry"] and ts < ge and gs < te:
            d = 1 - m / qlen
            ds.append(d); per[g].append(d)
    ds.sort()
    p99 = ds[min(len(ds) - 1, int(0.99 * len(ds)))]
    med = statistics.median(ds)
    json.dump(dict(delta_y=p99, n_outputs=len(ds), median=med, p95=ds[int(0.95 * len(ds))],
                   per_gene={g: dict(n=len(v), median=statistics.median(v), max=max(v)) for g, v in per.items()}),
              open(f"{a.w}/delta_y.json", "w"), indent=1)
    print(f"delta_Y = p99 of d over {len(ds)} single-copy-gene outputs ({len(per)} genes) = {p99:.5f}; median {med:.5f}, p95 "
          f"{ds[int(0.95 * len(ds))]:.5f}; per gene max: " + ", ".join(f"{g}:{max(v):.4f}(n={len(v)})" for g, v in sorted(per.items())))


def dmin(a):
    B = best_hits(f"{a.w}/masked_tx.masked.paf")
    U = best_hits(f"{a.w}/masked_tx.base.paf") if os.path.exists(f"{a.w}/masked_tx.base.paf") else {}
    with open(f"{a.w}/dmin.tsv", "w") as o:
        o.write("copy\ttx_len\td_min\tnearest\td_unmasked\n")
        for ln in open(f"{a.w}/masked_tx.fa"):
            if ln[0] != ">":
                continue
            n = ln[1:].strip()
            b = B.get(n)
            d = 1 - (b[4] / b[5] if b else 0.0)
            u = U.get(n)
            o.write(f"{n}\t{b[5] if b else 0}\t{d:.5f}\t{(b[1] + ':' + str(b[2]) + '-' + str(b[3])) if b else 'none'}\t"
                    f"{(1 - u[4] / u[5]) if u else float('nan'):.5f}\n")
            print(f"  {n}: d_min {d:.5f} nearest {(b[1] + ':' + str(b[2])) if b else 'none'}")


def net(a):
    P = {p["fam"]: p for p in json.load(open(f"{a.w}/panel.json"))}
    lab = {r["read"]: r for r in csv.DictReader(open(f"{a.w}/labels.tsv"), delimiter="\t")}
    rec = collections.defaultdict(list)
    for rd in pysam.AlignmentFile(f"{a.w}/R.bam").fetch(until_eof=True):
        rec[rd.query_name].append(None if rd.is_unmapped else (rd.reference_name, rd.reference_start, rd.reference_end))
    seq, cur = {}, None
    for ln in open(f"{a.w}/scored.fa"):
        if ln[0] == ">":
            cur = ln[1:].strip()
        else:
            seq[cur] = ln.strip()
    os.makedirs(f"{a.w}/fam", exist_ok=True)
    rng = random.Random(1)
    by = collections.defaultdict(list)
    for n, r in lab.items():
        rs = rec.get(n, [])
        unm = bool(rs) and all(x is None for x in rs)
        ons = any(x and x[0] == k[0] and x[1] < k[2] and k[1] < x[2] for x in rs for k in P[r["family"]]["keep"])
        if unm or ons:
            by[r["family"]].append(n)
    for f, ns in by.items():
        ns = sorted(ns); rng.shuffle(ns); ns = ns[:1000]
        with open(f"{a.w}/fam/{f}.fa", "w") as o:
            for n in ns:
                o.write(f">{n}\n{seq[n]}\n")
    print("IsoCon inputs:", {f: len(v) for f, v in by.items()})


def outputs(a):
    n = 0
    with open(f"{a.w}/outputs.fa", "w") as o:
        for path in sorted(glob.glob(f"{a.w}/iso/*/final_candidates.fa")):
            fam = path.split("/")[-2]
            for ln in open(path):
                if ln[0] == ">":
                    o.write(f">{fam}|{ln[1:].strip()}\n"); n += 1
                else:
                    o.write(ln)
    print("outputs", n)


def contigs(a):
    P = {p["fam"]: p for p in json.load(open(f"{a.w}/panel.json"))}
    delta = a.delta if a.delta is not None else json.load(open(f"{a.w}/delta_y.json"))["delta_y"]
    M, B = best_hits(f"{a.w}/outputs.masked.paf"), best_hits(f"{a.w}/outputs.base.paf")
    seqs, cur = {}, None
    for ln in open(f"{a.w}/outputs.fa"):
        if ln[0] == ">":
            cur = ln[1:].strip(); seqs[cur] = []
        else:
            seqs[cur].append(ln.strip())
    nI = nL = 0
    k = collections.Counter()
    with open(f"{a.w}/contigs_I.fa", "w") as fi, open(f"{a.w}/contigs_L.fa", "w") as fl, open(f"{a.w}/contigs.tsv", "w") as t:
        t.write("contig\tfamily\toutput\tlength\tbest_masked\td\tlinked\tsource\tsource_copy\n")
        for o, sl in seqs.items():
            m = M.get(o)
            if m and m[0] >= 0.999:
                continue
            fam = o.split("|")[0]
            s = "".join(sl)
            d = 1 - (m[4] / len(s) if m else 0.0)
            linked = d <= delta
            b = B.get(o)
            src, scopy = "none", ""
            if b:
                src = "elsewhere"
                for lab_, (c, s0, e, g) in [("D", mm) for mm in P[fam]["mask"]] + [("S", kk) for kk in P[fam]["keep"]]:
                    if b[1] == c and b[2] < e and s0 < b[3]:
                        src, scopy = (lab_ if lab_ == "D" else "S:" + g), g; break
            ctg = f"iso_{fam}_{k[fam]}"; k[fam] += 1
            fi.write(f">{ctg}\n{s}\n"); nI += 1
            if not linked:
                fl.write(f">{ctg}\n{s}\n"); nL += 1
            t.write(f"{ctg}\t{fam}\t{o}\t{len(s)}\t{m[0] if m else 0:.4f}\t{d:.5f}\t{int(linked)}\t{src}\t{scopy}\n")
    rows = list(csv.DictReader(open(f"{a.w}/contigs.tsv"), delimiter="\t"))
    print(f"delta used {delta:.5f}: outputs {len(seqs)}; flagged {nI}; linked back {nI - nL}; new copies {nL}; new-copy sources",
          dict(collections.Counter(r["source"].split(":")[0] for r in rows if r["linked"] == "0")))
    if nL:
        merge_test.pairs(a); merge_test.components(a)


def score(a):
    P = {p["fam"]: p for p in json.load(open(f"{a.w}/panel.json"))}
    delta = a.delta if a.delta is not None else json.load(open(f"{a.w}/delta_y.json"))["delta_y"]
    lab = {r["read"]: r for r in csv.DictReader(open(f"{a.w}/labels.tsv"), delimiter="\t")}
    ctg = {r["contig"]: r for r in csv.DictReader(open(f"{a.w}/contigs.tsv"), delimiter="\t")}
    dm = {r["copy"]: float(r["d_min"]) for r in csv.DictReader(open(f"{a.w}/dmin.tsv"), delimiter="\t")}
    rows = merge_test.rows_new(a.w)
    fam_of = {r["contig"]: r["family"] for r in rows}
    best = merge_test.best_pairs(a.w, sorted(set(fam_of.values())))
    par = {c: c for c in fam_of}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for (p, q), (m, cov, de) in best.items():
        if cov >= 0.5 and de <= delta:
            par[find(p)] = find(q)
    comp = {c: "cmp:" + find(c) for c in fam_of}
    members = collections.defaultdict(list)
    for c, cid in comp.items():
        members[cid].append(c)
    holds = {cid: {(ctg[c]["source"], ctg[c]["source_copy"]) for c in cs} for cid, cs in members.items()}
    # flags per masked copy
    inp = collections.Counter()
    for p in P.values():
        path = f"{a.w}/fam/{p['fam']}.fa"
        if os.path.exists(path):
            for ln in open(path):
                if ln[0] == ">":
                    r = lab.get(ln[1:].strip())
                    if r and r["role"] == "D":
                        inp[r["copy"]] += 1
    want = set(lab)

    def locus(chrom, s, e, fam):
        if chrom.startswith("iso_"):
            return ("ctg", comp[chrom])
        for c, s0, e0, g in P[fam]["mask"] + P[fam]["keep"]:
            if chrom == c and s < e0 and s0 < e:
                return ("copy", g)
        return ("other", f"{chrom}:{s // 100000}")

    def calls(bam):
        recs = collections.defaultdict(list)
        for rd in pysam.AlignmentFile(bam).fetch(until_eof=True):
            if rd.query_name not in want or rd.is_supplementary:
                continue
            recs[rd.query_name].append(None if rd.is_unmapped else (rd.is_secondary, rd.reference_name, rd.reference_start,
                                                                   rd.reference_end, rd.get_tag("AS") if rd.has_tag("AS") else 0))
        out = {}
        for n, rs in recs.items():
            rs = [r for r in rs if r]
            if not rs:
                out[n] = ("unplaced", None); continue
            fam = lab[n]["family"]
            prim = next((r for r in rs if not r[0]), rs[0])
            srt = sorted(rs, key=lambda r: -r[4])
            if len(srt) > 1 and srt[1][4] > 0 and srt[1][4] >= 0.98 * srt[0][4]:
                if len({locus(r[1], r[2], r[3], fam) for r in srt if r[4] >= 0.98 * srt[0][4]}) > 1:
                    out[n] = ("unplaced", None); continue
            out[n] = ("placed", locus(prim[1], prim[2], prim[3], fam))
        return out

    def cls(n, call):
        r = lab[n]
        st, L = call
        if st == "unplaced":
            return "unplaced"
        kind, key = L
        if r["role"] == "D":
            return "right" if kind == "ctg" and ("D", r["copy"]) in holds[key] else "wrong"
        if kind == "copy" and key == r["copy"]:
            return "stay"
        if kind == "ctg":
            return "stay" if ("S:" + r["copy"], r["copy"]) in holds[key] else "false_move"
        return "elsewhere"
    res, per = {}, {}
    for arm in ("R", "M"):
        if not os.path.exists(f"{a.w}/{arm}.bam"):
            continue
        C = calls(f"{a.w}/{arm}.bam")
        res[arm] = collections.Counter(); per[arm] = collections.defaultdict(collections.Counter)
        for n, r in lab.items():
            k = cls(n, C.get(n, ("unplaced", None)))
            res[arm][(r["role"], k)] += 1; per[arm][r["copy"]][k] += 1
        print(f"[{arm}] D:", {k[1]: v for k, v in sorted(res[arm].items()) if k[0] == "D"}, "| S:", {k[1]: v for k, v in sorted(res[arm].items()) if k[0] == "S"})
    print(f"\ndelta used {delta:.5f}; per masked copy (prediction: flag iff d_min > delta and input reads >= 6):")
    agree = tot = 0
    for p in P.values():
        for c, s0, e0, g in p["mask"]:
            flagged = any(len(cs) >= 2 and ("D", g) in holds[cid] for cid, cs in members.items())
            pred = dm.get(g, float("nan")) > delta and inp[g] >= 6
            ok = flagged == pred
            tot += 1; agree += ok
            pr = per.get("M", {}).get(g, {}); pr0 = per.get("R", {}).get(g, {})
            print(f"  {g:40s} d_min {dm.get(g, float('nan')):.5f} input {inp[g]:4d} pred {int(pred)} flagged {int(flagged)} {'OK ' if ok else 'MISS'}"
                  f" | D right/wrong/unplaced R {pr0.get('right', 0)}/{pr0.get('wrong', 0)}/{pr0.get('unplaced', 0)} -> M {pr.get('right', 0)}/{pr.get('wrong', 0)}/{pr.get('unplaced', 0)}")
    print(f"Y1: prediction right for {agree}/{tot} masked copies = {agree / max(1, tot):.1%} (bar 80%) -> {'HOLDS' if agree / max(1, tot) >= 0.8 else 'FAILS'}")
    if "M" in res:
        nS = sum(1 for r in lab.values() if r["role"] == "S")
        fm = res["M"][("S", "false_move")] / nS
        print(f"Y2: false moves {res['M'][('S', 'false_move')]}/{nS} = {fm:.2%} (bar 5%) -> {'HOLDS' if fm <= 0.05 else 'FAILS'}")
    json.dump({arm: {"|".join(k): v for k, v in c.items()} for arm, c in res.items()}, open(f"{a.w}/score.json", "w"), indent=0)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["panel", "mask", "reads", "scdelta", "dmin", "net", "outputs", "contigs", "score"])
    ap.add_argument("--w", required=True)
    ap.add_argument("--sub", default="human", choices=sorted(SUB))
    ap.add_argument("--delta", type=float, default=None)
    a = ap.parse_args(argv)
    os.makedirs(a.w, exist_ok=True)
    globals()[a.cmd](a)


if __name__ == "__main__":
    main()
