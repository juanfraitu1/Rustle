#!/usr/bin/env python3
"""Amendment 7 (docs/PREREG_rna_allele_haplotype_count_2026-10-01.md): linking IsoCon transcripts to their source locus, held-out on
multi-copy families. Subcommands, run in order:

  panel     families with >= 3 copies (all listed, >= 20 clean reads each, spans <= 200 kb); mask the copy last by (chrom, start); G3
  mask      copy of `_pri` with the masked intervals set to N (in place, by .fai offsets)
  reads     <= 500 baseline primaries per copy (seed 1) from the fibroblast BAM; scored.fa (read orientation) + labels.tsv
  net       per-family IsoCon input from the R-arm BAM (record on a surviving copy, or unmapped); <= 1,000 reads (seed 1)
  contigs   IsoCon outputs -> flagged (id x cov < 0.999 in the masked genome) -> linking (d <= delta) -> contigs_I.fa, contigs_L.fa
  score     arms R, R+I, R+I+L: final calls, labels, the registered rules

Work dir: --w (e.g. /mnt/linuxdisk/tmp/rna_allele/linktest).
"""
import argparse
import collections
import csv
import glob
import json
import os
import random
import shutil
import statistics
import subprocess

import pysam

SRC = "/home/juanfra/winloci_scratch/o3_collapse/method/intervals/data/intervals.tsv"
PRI = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta"
BAM = "/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam"
DELTA = 0.00958
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def panel(a):
    rows = list(csv.DictReader(open(SRC), delimiter="\t"))
    byfam = collections.defaultdict(list)
    for r in rows:
        byfam[r["fam"]].append(r)
    out, dropped = [], []
    allv = [(r["chrom"], int(r["clean_start"]), int(r["clean_end"]), r["gene"]) for r in rows]
    for fam, v in sorted(byfam.items()):
        if int(v[0]["n_copies"]) < 3 or len(v) != int(v[0]["n_copies"]):
            continue
        if min(int(x["n_clean"]) for x in v) < 20 or max(int(x["clean_len"]) for x in v) > 200_000:
            continue
        cps = sorted(v, key=lambda r: (r["chrom"], int(r["clean_start"])))
        m = cps[-1]
        mc, ms, me = m["chrom"], int(m["clean_start"]), int(m["clean_end"])
        if any(c == mc and s < me and ms < e and g != m["gene"] for c, s, e, g in allv):
            dropped.append(fam)
            continue
        out.append(dict(fam=fam, mask=[mc, ms, me, m["gene"]],
                        keep=[[x["chrom"], int(x["clean_start"]), int(x["clean_end"]), x["gene"]] for x in cps[:-1]]))
    json.dump(out, open(f"{a.w}/panel.json", "w"), indent=0)
    print(f"families {len(out)} (dropped by G3: {len(dropped)} {dropped}); copies {sum(1 + len(p['keep']) for p in out)}")


def mask(a):
    P = json.load(open(f"{a.w}/panel.json"))
    dst = f"{a.w}/masked.fa"
    shutil.copy(PRI, dst); shutil.copy(PRI + ".fai", dst + ".fai")
    fai = {l.split("\t")[0]: l.split("\t") for l in open(dst + ".fai")}
    tot = 0
    with open(dst, "r+b") as f:
        for p in P:
            c, s0, e, _ = p["mask"]
            _, _, off, lb, lw = fai[c][:5]
            off, lb, lw = int(off), int(lb), int(lw)
            pos = s0
            while pos < e:
                line, col = divmod(pos, lb)
                k = min(lb - col, e - pos)
                f.seek(off + line * lw + col); f.write(b"N" * k); pos += k; tot += k
    print("masked bp", tot, "expected", sum(p["mask"][2] - p["mask"][1] for p in P))


def reads(a):
    P = json.load(open(f"{a.w}/panel.json"))
    bam = pysam.AlignmentFile(BAM)
    rng = random.Random(1)
    role = {}
    for p in P:
        for lab, (c, s, e, g) in [("D", p["mask"])] + [("S", k) for k in p["keep"]]:
            ns = sorted({r.query_name for r in bam.fetch(c, s, e) if not (r.is_unmapped or r.is_secondary or r.is_supplementary)})
            rng.shuffle(ns)
            for n in ns[:500]:
                role.setdefault(n, (p["fam"], lab, g))
    with open(f"{a.w}/names.txt", "w") as o:
        o.write("\n".join(sorted(role)) + "\n")
    subprocess.run(f"samtools view -F 2308 -N {a.w}/names.txt -@ 4 {BAM} -o {a.w}/prim.sam", shell=True, check=True)
    seq = {}
    for ln in open(f"{a.w}/prim.sam"):
        f = ln.split("\t", 11)
        s = f[9]
        if int(f[1]) & 16:
            s = s.translate(COMP)[::-1]
        seq[f[0]] = s
    with open(f"{a.w}/scored.fa", "w") as fa, open(f"{a.w}/labels.tsv", "w") as lab:
        lab.write("read\tfamily\trole\tcopy\n")
        for n, (f, r, g) in sorted(role.items()):
            if n in seq:
                fa.write(f">{n}\n{seq[n]}\n"); lab.write(f"{n}\t{f}\t{r}\t{g}\n")
    c = collections.Counter(r for _, r, _ in role.values())
    print(f"reads: D {c['D']}, S {c['S']}; sequences {len(seq)}")


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
    sizes = []
    for f, ns in by.items():
        ns = sorted(ns); rng.shuffle(ns); ns = ns[:1000]; sizes.append(len(ns))
        with open(f"{a.w}/fam/{f}.fa", "w") as o:
            for n in ns:
                o.write(f">{n}\n{seq[n]}\n")
    print(f"IsoCon inputs: {len(by)} families, total {sum(sizes)}, median {sorted(sizes)[len(sizes) // 2]}")


def best_hits(path):
    b = {}
    for ln in open(path):
        f = ln.split("\t")
        s = int(f[9]) / max(1, int(f[10])) * (int(f[3]) - int(f[2])) / max(1, int(f[1]))
        if f[0] not in b or s > b[f[0]][0]:
            b[f[0]] = (s, f[5], int(f[7]), int(f[8]), int(f[9]), int(f[1]))
    return b


def contigs(a):
    P = {p["fam"]: p for p in json.load(open(f"{a.w}/panel.json"))}
    M, B = best_hits(f"{a.w}/outputs.masked.paf"), best_hits(f"{a.w}/outputs.base.paf")
    seqs, cur = {}, None
    for ln in open(f"{a.w}/outputs.fa"):
        if ln[0] == ">":
            cur = ln[1:].strip(); seqs[cur] = []
        else:
            seqs[cur].append(ln.strip())
    nI = nL = 0
    with open(f"{a.w}/contigs_I.fa", "w") as fi, open(f"{a.w}/contigs_L.fa", "w") as fl, open(f"{a.w}/contigs.tsv", "w") as t:
        t.write("contig\tfamily\toutput\tlength\tbest_masked\td\tlinked\tsource\n")
        k = collections.Counter()
        for o, sl in seqs.items():
            m = M.get(o)
            if m and m[0] >= 0.999:
                continue
            fam = o.split("|")[0]
            s = "".join(sl)
            d = 1 - (m[4] / len(s) if m else 0.0)
            linked = d <= DELTA
            b = B.get(o)
            src = "none"
            if b:
                for lab_, (c, s0, e, g) in [("D", P[fam]["mask"])] + [("S:" + kk[3], kk) for kk in P[fam]["keep"]]:
                    if b[1] == c and b[2] < e and s0 < b[3]:
                        src = lab_; break
                else:
                    src = "elsewhere"
            ctg = f"iso_{fam}_{k[fam]}"; k[fam] += 1
            fi.write(f">{ctg}\n{s}\n"); nI += 1
            if not linked:
                fl.write(f">{ctg}\n{s}\n"); nL += 1
            t.write(f"{ctg}\t{fam}\t{o}\t{len(s)}\t{m[0] if m else 0:.4f}\t{d:.5f}\t{int(linked)}\t{src}\n")
    rows = list(csv.DictReader(open(f"{a.w}/contigs.tsv"), delimiter="\t"))
    print(f"outputs {len(seqs)}; flagged {nI}; kept as new copies after linking {nL}")
    for lk in ("0", "1"):
        print(f"  linked={lk}: by source", dict(collections.Counter(r["source"] for r in rows if r["linked"] == lk)))


def score(a):
    P = {p["fam"]: p for p in json.load(open(f"{a.w}/panel.json"))}
    lab = {r["read"]: r for r in csv.DictReader(open(f"{a.w}/labels.tsv"), delimiter="\t")}
    ctg = {r["contig"]: r for r in csv.DictReader(open(f"{a.w}/contigs.tsv"), delimiter="\t")}
    want = set(lab)

    def locus(chrom, s, e, fam):
        if chrom.startswith("iso_"):
            return ("ctg", chrom)
        c, s0, e0, g = P[fam]["mask"]
        for kk in P[fam]["keep"]:
            if chrom == kk[0] and s < kk[2] and kk[1] < e:
                return ("copy", kk[3])
        return ("other", f"{chrom}:{s // 100000}")

    def calls(bam):
        recs = collections.defaultdict(list)
        for rd in pysam.AlignmentFile(bam).fetch(until_eof=True):
            if rd.query_name not in want or rd.is_supplementary:
                continue
            recs[rd.query_name].append(None if rd.is_unmapped else
                                       (rd.is_secondary, rd.reference_name, rd.reference_start, rd.reference_end,
                                        rd.get_tag("AS") if rd.has_tag("AS") else 0))
        out = {}
        for n, rs in recs.items():
            rs = [r for r in rs if r]
            if not rs:
                out[n] = ("unplaced", None); continue
            fam = lab[n]["family"]
            prim = next((r for r in rs if not r[0]), rs[0])
            srt = sorted(rs, key=lambda r: -r[4])
            L1 = locus(prim[1], prim[2], prim[3], fam)
            if len(srt) > 1 and srt[1][4] > 0 and srt[1][4] >= 0.98 * srt[0][4]:
                Ls = {locus(r[1], r[2], r[3], fam) for r in srt if r[4] >= 0.98 * srt[0][4]}
                if len(Ls) > 1:
                    out[n] = ("unplaced", None); continue
            out[n] = ("placed", L1)
        return out

    def cls(n, call):
        r = lab[n]
        st, L = call
        if st == "unplaced":
            return "unplaced"
        kind, key = L
        if r["role"] == "D":
            return "right" if kind == "ctg" and ctg[key]["source"] == "D" and ctg[key]["family"] == r["family"] else "wrong"
        if kind == "copy" and key == r["copy"]:
            return "stay"
        if kind == "ctg":
            return "stay" if (ctg[key]["family"] == r["family"] and ctg[key]["source"] == "S:" + r["copy"]) else "false_move"
        return "elsewhere"
    res, per = {}, {}
    for arm in ("R", "RI", "RIL"):
        C = calls(f"{a.w}/{arm}.bam")
        res[arm] = collections.Counter()
        per[arm] = collections.defaultdict(collections.Counter)
        for n, r in lab.items():
            k = cls(n, C.get(n, ("unplaced", None)))
            res[arm][(r["role"], k)] += 1
            per[arm][r["family"]][(r["role"], k)] += 1
        print(f"[{arm}] D:", {k[1]: v for k, v in res[arm].items() if k[0] == "D"}, "| S:", {k[1]: v for k, v in res[arm].items() if k[0] == "S"})
    nS = sum(1 for r in lab.values() if r["role"] == "S")
    uI, uL = res["RI"][("S", "unplaced")], res["RIL"][("S", "unplaced")]
    rI, rL = res["RI"][("D", "right")], res["RIL"][("D", "right")]
    a1 = (uI - uL) / uI if uI else 0.0
    a2 = (rI - rL) / rI if rI else 0.0
    print(f"LINKING: S unplaced {uI} -> {uL} (fall {a1:.1%}); D right {rI} -> {rL} (fall {a2:.1%}) -> "
          f"{'WORKS' if a1 >= 0.5 and a2 <= 0.2 else 'DOES NOT WORK'}")
    w0, w1 = res["R"][("D", "wrong")], res["RIL"][("D", "wrong")]
    fm = res["RIL"][("S", "false_move")] / nS
    drop = (w0 - w1) / w0 if w0 else 0.0
    print(f"OVERALL (R vs R+I+L): wrong D {w0} -> {w1} (drop {drop:.1%}); false moves {fm:.1%} -> "
          f"{'HELP' if drop >= 0.5 and fm <= 0.05 else 'HURT' if fm > 0.10 else 'MIXED'}")
    json.dump({arm: {f: {'|'.join(k): v for k, v in c.items()} for f, c in per[arm].items()} for arm in per},
              open(f"{a.w}/per_family_result.json", "w"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["panel", "mask", "reads", "net", "contigs", "score"])
    ap.add_argument("--w", required=True)
    a = ap.parse_args(argv)
    os.makedirs(a.w, exist_ok=True)
    globals()[a.cmd](a)


if __name__ == "__main__":
    main()
