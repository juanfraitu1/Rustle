#!/usr/bin/env python3
"""Amendment 8 (docs/PREREG_rna_allele_haplotype_count_2026-10-01.md): merging a missing copy's new-copy transcripts into one candidate
copy without truth, scored on Amendment 7's 53 families and its R+I+L alignments (link_test.py's work dir, nothing realigned).

  pairs       per family, the new-copy contigs (contigs.tsv linked=0) all-vs-all: minimap2 -c -x asm20 --cs -N 200 -p 0.1
  components  best alignment per unordered pair -> joined when it covers >= 50% of the shorter contig and de <= delta -> components
              (at delta/2, delta, 2 x delta); merge/components.tsv, merge/summary.json
  score       R+I+L calls with components as loci (arm M), the truth-grouped ceiling (arm T), rules M1 / M2

    merge_test.py pairs --w /mnt/linuxdisk/tmp/rna_allele/linktest
"""
import argparse
import collections
import csv
import json
import os
import subprocess

import pysam

DELTA = 0.00958
DELTAS = {"half": DELTA / 2, "delta": DELTA, "double": 2 * DELTA}


def rows_new(w):
    return [r for r in csv.DictReader(open(f"{w}/contigs.tsv"), delimiter="\t") if r["linked"] == "0"]


def pairs(a):
    os.makedirs(f"{a.w}/merge/fam", exist_ok=True); os.makedirs(f"{a.w}/merge/paf", exist_ok=True)
    seq, cur = {}, None
    for ln in open(f"{a.w}/contigs_L.fa"):
        if ln[0] == ">":
            cur = ln[1:].strip(); seq[cur] = []
        else:
            seq[cur].append(ln.strip())
    by = collections.defaultdict(list)
    for r in rows_new(a.w):
        by[r["family"]].append(r["contig"])
    for fam, cs in sorted(by.items()):
        fa = f"{a.w}/merge/fam/{fam}.fa"
        with open(fa, "w") as o:
            for c in cs:
                o.write(f">{c}\n{''.join(seq[c])}\n")
        with open(f"{a.w}/merge/paf/{fam}.paf", "w") as o, open(f"{a.w}/merge/paf/{fam}.log", "w") as e:
            subprocess.run(["minimap2", "-c", "-x", "asm20", "--cs", "-N", "200", "-p", "0.1", "-t", "2", fa, fa], stdout=o, stderr=e, check=True)
    print(f"families {len(by)}, contigs {sum(len(v) for v in by.values())}, alignments "
          f"{sum(1 for f in by for _ in open(f'{a.w}/merge/paf/{f}.paf'))}")


def best_pairs(w, fams):
    """best alignment (by matching bases) per unordered pair: (cov of the shorter contig, de)"""
    best = {}
    for fam in fams:
        for ln in open(f"{w}/merge/paf/{fam}.paf"):
            f = ln.rstrip("\n").split("\t")
            q, t = f[0], f[5]
            if q == t:
                continue
            qlen, qs, qe, tlen, ts, te, m = int(f[1]), int(f[2]), int(f[3]), int(f[6]), int(f[7]), int(f[8]), int(f[9])
            de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), 1.0)
            cov = (qe - qs) / qlen if qlen <= tlen else (te - ts) / tlen
            key = (q, t) if q < t else (t, q)
            if key not in best or m > best[key][0]:
                best[key] = (m, cov, de)
    return best


def components(a):
    rows = rows_new(a.w)
    fam_of = {r["contig"]: r["family"] for r in rows}
    src = {r["contig"]: r["source"] for r in rows}
    fams = sorted({r["family"] for r in rows})
    best = best_pairs(a.w, fams)
    out = {}
    for name, delta in DELTAS.items():
        par = {c: c for c in fam_of}

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]; x = par[x]
            return x
        for (p, q), (m, cov, de) in best.items():
            if cov >= 0.5 and de <= delta:
                par[find(p)] = find(q)
        comp = collections.defaultdict(list)
        for c in fam_of:
            comp[find(c)].append(c)
        ids, k = {}, collections.Counter()
        for root, cs in sorted(comp.items(), key=lambda kv: (fam_of[kv[0]], kv[0])):
            fam = fam_of[root]
            cid = f"{fam}:c{k[fam]}"; k[fam] += 1
            for c in cs:
                ids[c] = cid
        # per-family structure
        per, kinds = {}, collections.Counter()
        for fam in fams:
            cs = [c for c in fam_of if fam_of[c] == fam]
            cids = {ids[c] for c in cs}
            dcomp = {ids[c] for c in cs if src[c] == "D"}
            nd = sum(1 for c in cs if src[c] == "D")
            kind = collections.Counter()
            for cid in cids:
                mem = [c for c in cs if ids[c] == cid]
                hasD, hasS = any(src[c] == "D" for c in mem), any(src[c].startswith("S:") for c in mem)
                kind["mixed" if hasD and hasS else "D" if hasD else "S_only" if hasS else "other"] += 1
            kinds.update(kind)
            per[fam] = dict(contigs=len(cs), components=len(cids), D_contigs=nd, D_components=len(dcomp), **kind)
        two = [f for f in fams if per[f]["D_contigs"] >= 2]
        one = sum(1 for f in two if per[f]["D_components"] == 1)
        out[name] = dict(delta=delta, components=sum(p["components"] for p in per.values()), kinds=dict(kinds),
                         fam_ge2_D=len(two), fam_one_D_component=one, per_family=per)
        if name == "delta":
            with open(f"{a.w}/merge/components.tsv", "w") as o:
                o.write("contig\tfamily\tcomponent\tsize\tsource\n")
                for c in sorted(fam_of, key=lambda c: (fam_of[c], ids[c], c)):
                    o.write(f"{c}\t{fam_of[c]}\t{ids[c]}\t{sum(1 for x in fam_of if ids[x] == ids[c])}\t{src[c]}\n")
            # over-split causes: pairs of D-derived contigs of one family in different components
            cause = collections.Counter()
            for f in two:
                if per[f]["D_components"] == 1:
                    continue
                ds = [c for c in fam_of if fam_of[c] == f and src[c] == "D"]
                for i in range(len(ds)):
                    for j in range(i + 1, len(ds)):
                        if ids[ds[i]] == ids[ds[j]]:
                            continue
                        b = best.get((ds[i], ds[j]) if ds[i] < ds[j] else (ds[j], ds[i]))
                        cause["no alignment" if not b else "coverage < 0.5" if b[1] < 0.5 else "de > delta"] += 1
            out["oversplit_pair_causes"] = dict(cause)
        print(f"[{name} delta={delta:.5f}] components {out[name]['components']} {dict(kinds)}; families with >= 2 D contigs {len(two)}, "
              f"D contigs in ONE component {one} ({one / max(1, len(two)):.1%})")
    print("over-split D pairs by cause:", out["oversplit_pair_causes"])
    json.dump(out, open(f"{a.w}/merge/summary.json", "w"), indent=0)


def score(a):
    P = {p["fam"]: p for p in json.load(open(f"{a.w}/panel.json"))}
    lab = {r["read"]: r for r in csv.DictReader(open(f"{a.w}/labels.tsv"), delimiter="\t")}
    ctg = {r["contig"]: r for r in csv.DictReader(open(f"{a.w}/contigs.tsv"), delimiter="\t")}
    S = json.load(open(f"{a.w}/merge/summary.json"))
    rows = rows_new(a.w)
    fam_of = {r["contig"]: r["family"] for r in rows}
    best = best_pairs(a.w, sorted(set(fam_of.values())))

    def comp_map(delta):
        par = {c: c for c in fam_of}

        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]; x = par[x]
            return x
        for (p, q), (m, cov, de) in best.items():
            if cov >= 0.5 and de <= delta:
                par[find(p)] = find(q)
        return {c: "cmp:" + find(c) for c in fam_of}
    # loci maps: L = each contig its own locus; M = component; T = truth group (source)
    maps = {"L": {c: c for c in fam_of},
            "T": {c: f"src:{fam_of[c]}:{ctg[c]['source']}" if ctg[c]["source"] == "D" or ctg[c]["source"].startswith("S:") else c for c in fam_of}}
    for name, delta in DELTAS.items():
        maps["M" if name == "delta" else f"M_{name}"] = comp_map(delta)
    # which sources each locus holds, and its family (for labels)
    holds, lfam = {}, {}
    for arm, mp in maps.items():
        h = collections.defaultdict(set)
        for c, L in mp.items():
            h[L].add(ctg[c]["source"])
        holds[arm] = h
        lfam[arm] = {L: fam_of[c] for c, L in mp.items()}
    want = set(lab)

    def load(bam):
        recs = collections.defaultdict(list)
        for rd in pysam.AlignmentFile(bam).fetch(until_eof=True):
            if rd.query_name not in want or rd.is_supplementary:
                continue
            recs[rd.query_name].append(None if rd.is_unmapped else
                                       (rd.is_secondary, rd.reference_name, rd.reference_start, rd.reference_end,
                                        rd.get_tag("AS") if rd.has_tag("AS") else 0))
        return recs

    def locus(chrom, s, e, fam, mp):
        if chrom.startswith("iso_"):
            return ("ctg", mp[chrom])
        for kk in P[fam]["keep"]:
            if chrom == kk[0] and s < kk[2] and kk[1] < e:
                return ("copy", kk[3])
        return ("other", f"{chrom}:{s // 100000}")

    def calls(recs, mp):
        out = {}
        for n, rs in recs.items():
            rs = [r for r in rs if r]
            if not rs:
                out[n] = ("unplaced", None); continue
            fam = lab[n]["family"]
            prim = next((r for r in rs if not r[0]), rs[0])
            srt = sorted(rs, key=lambda r: -r[4])
            if len(srt) > 1 and srt[1][4] > 0 and srt[1][4] >= 0.98 * srt[0][4]:
                if len({locus(r[1], r[2], r[3], fam, mp) for r in srt if r[4] >= 0.98 * srt[0][4]}) > 1:
                    out[n] = ("unplaced", None); continue
            out[n] = ("placed", locus(prim[1], prim[2], prim[3], fam, mp))
        return out

    def cls(n, call, mk):
        r = lab[n]
        st, L = call
        if st == "unplaced":
            return "unplaced"
        kind, key = L
        own = kind == "ctg" and lfam[mk][key] == r["family"]
        if r["role"] == "D":
            return "right" if own and "D" in holds[mk][key] else "wrong"
        if kind == "copy" and key == r["copy"]:
            return "stay"
        if kind == "ctg":
            return "stay" if own and ("S:" + r["copy"]) in holds[mk][key] else "false_move"
        return "elsewhere"
    recs = {"R": load(f"{a.w}/R.bam"), "RIL": load(f"{a.w}/RIL.bam")}
    res = {}
    arms = [("R", "R", "L"), ("RIL", "RIL", "L"), ("M", "RIL", "M"), ("T", "RIL", "T"), ("M_half", "RIL", "M_half"), ("M_double", "RIL", "M_double")]
    nD = sum(1 for r in lab.values() if r["role"] == "D"); nS = len(lab) - nD
    for arm, bam, mk in arms:
        C = calls(recs[bam], maps[mk])
        res[arm] = collections.Counter()
        for n, r in lab.items():
            res[arm][(r["role"], cls(n, C.get(n, ("unplaced", None)), mk))] += 1
        print(f"[{arm}] D:", {k[1]: v for k, v in sorted(res[arm].items()) if k[0] == "D"}, "| S:", {k[1]: v for k, v in sorted(res[arm].items()) if k[0] == "S"})
    rM, fmM = res["M"][("D", "right")], res["M"][("S", "false_move")] / nS
    print(f"M1: D right {res['RIL'][('D','right')]} -> {rM} (truth-grouped {res['T'][('D','right')]}; bar 6,440 = 50% of 12,879); "
          f"false moves {res['M'][('S','false_move')]}/{nS} = {fmM:.2%} -> {'WORKS' if rM >= 6440 and fmM <= 0.05 else 'DOES NOT WORK'}")
    w0, w1 = res["R"][("D", "wrong")], res["M"][("D", "wrong")]
    drop = (w0 - w1) / w0 if w0 else 0.0
    print(f"OVERALL (R vs R+I+L+M): wrong D {w0} -> {w1} (drop {drop:.1%}); false moves {fmM:.2%} -> "
          f"{'HELP' if drop >= 0.5 and fmM <= 0.05 else 'HURT' if fmM > 0.10 else 'MIXED'}")
    s = S["delta"]
    one, two = s["fam_one_D_component"], s["fam_ge2_D"]
    mixed = s["kinds"].get("mixed", 0) / max(1, s["components"])
    print(f"M2: D contigs in ONE component {one}/{two} = {one / max(1, two):.1%} (bar 2/3); mixed components {s['kinds'].get('mixed', 0)}/{s['components']} "
          f"= {mixed:.1%} (bar 10%) -> {'HOLDS' if one / max(1, two) >= 2 / 3 and mixed <= 0.10 else 'FAILS'}")
    json.dump({arm: {"|".join(k): v for k, v in c.items()} for arm, c in res.items()}, open(f"{a.w}/merge/score.json", "w"), indent=0)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["pairs", "components", "score"])
    ap.add_argument("--w", required=True)
    a = ap.parse_args(argv)
    globals()[a.cmd](a)


if __name__ == "__main__":
    main()
