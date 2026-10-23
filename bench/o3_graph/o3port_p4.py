#!/usr/bin/env python3
"""P4, the two-tier O3 output on a fresh gorilla fibroblast held-out (docs/PREREG_o3_two_tier_p4_2026-10-23.md).

  panel   --w W          the registered draw -> W/panel.json (link_test.py's format: fam, mask, keep)
  split   --w W          W/scored.fa -> W/scored.part{0,1,2}.fa (record k goes to part k % 3)
  batches --w W          the panel's families in panel order, 10 per line -> W/batches.txt
  verdict --w W          H4-1 / H4-2 / H4-3 and the reported tables from the arms' products -> W/verdict.json, stdout
"""
import argparse
import collections
import csv
import json
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "rna_allele"))

SRC = "/home/juanfra/winloci_scratch/o3_collapse/method/intervals/data/intervals.tsv"
RA = "/mnt/linuxdisk/tmp/rna_allele"
N_DRAW = 60
BATCH = 10


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def panel(a):
    rows = tsv(SRC)
    byfam = collections.defaultdict(list)
    for r in rows:
        byfam[r["fam"]].append(r)
    dev = {r["family_id"] for r in tsv(f"{RA}/a12/A12.copies.tsv")}
    bonly = {r["family"] for r in tsv(f"{RA}/refabsent/bonly.tsv")}
    assert len(dev) == 53 and len(bonly) == 34, (len(dev), len(bonly))
    allv = [(r["chrom"], int(r["clean_start"]), int(r["clean_end"]), r["gene"]) for r in rows]
    eligible, g3 = [], []
    for fam, v in sorted(byfam.items()):
        if fam in dev or fam in bonly or len(v) != int(v[0]["n_copies"]) or len(v) < 2:
            continue
        if min(int(x["n_clean"]) for x in v) < 20 or max(int(x["clean_len"]) for x in v) > 200_000:
            continue
        cps = sorted(v, key=lambda r: (r["chrom"], int(r["clean_start"])))
        m = cps[-1]
        mc, ms, me = m["chrom"], int(m["clean_start"]), int(m["clean_end"])
        if any(c == mc and s < me and ms < e and g != m["gene"] for c, s, e, g in allv):
            g3.append(fam)
            continue
        eligible.append(dict(fam=fam, mask=[mc, ms, me, m["gene"]],
                             keep=[[x["chrom"], int(x["clean_start"]), int(x["clean_end"]), x["gene"]] for x in cps[:-1]]))
    draw = sorted(random.Random(1).sample(eligible, min(N_DRAW, len(eligible))), key=lambda p: p["fam"])
    json.dump(draw, open(f"{a.w}/panel.json", "w"), indent=0)
    print(f"eligible {len(eligible)} (dropped by G3: {len(g3)}); drawn {len(draw)}; copies {sum(1 + len(p['keep']) for p in draw)}; "
          f"copies per family {dict(collections.Counter(1 + len(p['keep']) for p in draw))}")


def split(a):
    outs = [open(f"{a.w}/scored.part{p}.fa", "w") for p in range(3)]
    k, n = -1, collections.Counter()
    for ln in open(f"{a.w}/scored.fa"):
        if ln.startswith(">"):
            k += 1
            n[k % 3] += 1
        outs[k % 3].write(ln)
    for o in outs:
        o.close()
    print("records per part", dict(sorted(n.items())))


def batches(a):
    fams = [p["fam"] for p in json.load(open(f"{a.w}/panel.json"))]
    with open(f"{a.w}/batches.txt", "w") as o:
        for i in range(0, len(fams), BATCH):
            o.write(",".join(fams[i:i + BATCH]) + "\n")
    print(f"{len(fams)} families in {(len(fams) + BATCH - 1) // BATCH} batches")


def members(prefix):
    """{candidate: (family, member reads)} of a product's flagged candidates"""
    cand = {r["candidate"]: r for r in tsv(f"{prefix}.candidates.tsv") if r["flagged"] == "1"}
    owner = {r["cluster"]: r["candidate"] for r in tsv(f"{prefix}.clusters.tsv") if r["candidate"] in cand}
    m = collections.defaultdict(set)
    for r in tsv(f"{prefix}.reads.tsv"):
        if r["cluster"] in owner:
            m[owner[r["cluster"]]].add(r["read"])
    return {c: (cand[c]["family"], m[c]) for c in cand}


def verdict(a):
    import accept_o3_candidates as acc
    import pysam
    W = a.w
    P = {p["fam"]: p for p in json.load(open(f"{W}/panel.json"))}
    lab = {r["read"]: r for r in tsv(f"{W}/labels.tsv")}
    nD = sum(1 for r in lab.values() if r["role"] == "D")
    nS = len(lab) - nD
    # attainable D reads (Amendment 13b, truth-free): a record on a surviving copy in R.bam, or attributed (pass B) into the
    # read's own family's net (the stage arm's net replication on the deletion arm)
    on_surv = set()
    for rd in pysam.AlignmentFile(f"{W}/R.bam").fetch(until_eof=True):
        r = lab.get(rd.query_name)
        if rd.is_unmapped or rd.is_supplementary or not r or r["role"] != "D":
            continue
        if any(rd.reference_name == k[0] and rd.reference_start < k[2] and k[1] < rd.reference_end for k in P[r["family"]]["keep"]):
            on_surv.add(rd.query_name)
    att = {r["read"] for r in tsv(f"{W}/kd/net_reads.tsv") if r["via"] == "B" and r["read"] in lab
           and lab[r["read"]]["role"] == "D" and lab[r["read"]]["family"] == r["family"]} - on_surv
    attainable = on_surv | att
    out = {"n_D": nD, "n_S": nS, "attainable_on_survivor": len(on_surv), "attainable_attributed": len(att),
           "attainable": len(attainable)}
    for arm in ("pd", "kd"):
        calls = acc.arm_m(f"{W}/{arm}", lab, P)
        calls.pop("__false_moves__")
        tot = collections.Counter(calls.values())
        right_att = sum(1 for n in attainable if calls[n] == ("D", "right"))
        out[arm] = {"D_right": tot[("D", "right")], "D_right_attainable": right_att, "false_moves": tot[("S", "false_move")],
                    "per_deleted_copy_right": {f: sum(1 for n, r in lab.items() if r["family"] == f and calls[n] == ("D", "right"))
                                               for f in P}}
    # control arm: classes per flag (control_test.py classify's candidates_classified / classify.json) and arm C's false moves
    for arm in ("pc", "kc"):
        cj = json.load(open(f"{W}/{arm}/classify.json")) if os.path.exists(f"{W}/{arm}/classify.json") else None
        out[arm] = {"classify": cj}
        sj = f"{W}/{arm}/score.json"
        if os.path.exists(sj):
            out[arm]["score"] = json.load(open(sj))
    # tier 2 = the stage's flagged candidates sharing < half of their members with every tier-1 flag (per arm)
    for arm_k, arm_p in (("kd", "pd"), ("kc", "pc")):
        t1 = members(f"{W}/{arm_p}/cand")
        t2 = {c: v for c, v in members(f"{W}/{arm_k}/cand").items()
              if all(2 * len(v[1] & w[1]) < len(v[1]) for w in t1.values() if w[0] == v[0])}
        out[f"tier2_{arm_k}"] = {"candidates": len(t2), "families": sorted({v[0] for v in t2.values()})}
    json.dump(out, open(f"{W}/verdict.json", "w"), indent=1)
    pd, kd = out["pd"], out["kd"]
    print(f"P4d: D reads {nD}, S reads {nS}; attainable D {len(attainable)} (on a survivor {len(on_surv)}, attributed {len(att)})")
    print(f"  tier 1 (port): D right {pd['D_right']}, among attainable {pd['D_right_attainable']} "
          f"({pd['D_right_attainable'] / max(1, len(attainable)):.3f}); false moves {pd['false_moves']}/{nS} = {pd['false_moves'] / nS:.2%}")
    print(f"  stage arm (K): D right {kd['D_right']}, among attainable {kd['D_right_attainable']} "
          f"({kd['D_right_attainable'] / max(1, len(attainable)):.3f}); false moves {kd['false_moves']}/{nS} = {kd['false_moves'] / nS:.2%}")
    h43 = pd["D_right_attainable"] >= 0.5 * len(attainable)
    h42d = pd["false_moves"] <= 0.05 * nS
    print(f"H4-3 tier-1 D right among attainable >= 0.50 x {len(attainable)} = {0.5 * len(attainable):.1f}: {'PASS' if h43 else 'FAIL'}")
    print(f"H4-2 (P4d) false moves <= 5%: {'PASS' if h42d else 'FAIL'}")
    print("P4c classification and arm C: see pc/classify.out, pc/score.out (H4-1, H4-2 control part)")
    print("tier 2:", {k: v for k, v in out.items() if k.startswith("tier2_")})


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["panel", "split", "batches", "verdict"])
    ap.add_argument("--w", required=True)
    a = ap.parse_args(argv)
    {"panel": panel, "split": split, "batches": batches, "verdict": verdict}[a.cmd](a)


if __name__ == "__main__":
    main(sys.argv[1:])
