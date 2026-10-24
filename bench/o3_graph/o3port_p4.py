#!/usr/bin/env python3
"""P4, the two-tier O3 output on a fresh gorilla fibroblast held-out (docs/PREREG_o3_two_tier_p4_2026-10-23.md).

  panel   --w W          the registered draw -> W/panel.json (link_test.py's format: fam, mask, keep)
  split   --w W          W/scored.fa -> W/scored.part{0,1,2}.fa (record k goes to part k % 3)
  batches --w W          the panel's families in panel order, 10 per line -> W/batches.txt
  verdict --w W          H4-1 / H4-2 / H4-3 and the reported tables from the arms' products -> W/verdict.json, stdout
  subpanel --w W --from P4 --n 20   P5 (docs/PREREG_o3_port_v6_heldout_2026-10-23.md): n families drawn with random.Random(1) from
                         P4's draw (sorted by family) -> W/panel.json, and P4's scored reads of those families -> W/scored.fa,
                         W/labels.tsv
  verdict5 --w W --deletion pd,pdcc --control pc,pccc   P5's H5-1 / H5-2 / H5-3 per arm -> W/verdict.json, stdout
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


def eligible_families():
    """P4's eligibility (docs/PREREG_o3_two_tier_p4_2026-10-23.md): (eligible panel entries sorted by family, G3 drops)."""
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
    return eligible, g3


def panel(a):
    eligible, g3 = eligible_families()
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


def subpanel(a):
    P = sorted(json.load(open(f"{a.source}/panel.json")), key=lambda p: p["fam"])
    draw = sorted(random.Random(1).sample(P, a.n), key=lambda p: p["fam"])
    json.dump(draw, open(f"{a.w}/panel.json", "w"), indent=0)
    fams = {p["fam"] for p in draw}
    lab = tsv(f"{a.source}/labels.tsv")
    keep = {r["read"] for r in lab if r["family"] in fams}
    with open(f"{a.w}/labels.tsv", "w") as o:
        o.write("read\tfamily\trole\tcopy\n")
        for r in lab:
            if r["read"] in keep:
                o.write(f"{r['read']}\t{r['family']}\t{r['role']}\t{r['copy']}\n")
    n, name, out = 0, None, open(f"{a.w}/scored.fa", "w")
    for ln in open(f"{a.source}/scored.fa"):
        if ln.startswith(">"):
            name = ln[1:].strip()
            n += name in keep
        if name in keep:
            out.write(ln)
    out.close()
    roles = collections.Counter(r["role"] for r in lab if r["read"] in keep)
    print(f"drawn {len(draw)} of {len(P)}: {sorted(fams)}; reads {n} ({dict(roles)})")


def verdict5(a):
    import accept_o3_candidates as acc
    import pysam
    W = a.w
    P = {p["fam"]: p for p in json.load(open(f"{W}/panel.json"))}
    lab = {r["read"]: r for r in tsv(f"{W}/labels.tsv")}
    nS = sum(1 for r in lab.values() if r["role"] == "S")
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
    out = {"attainable": len(attainable), "on_survivor": len(on_surv), "attributed": len(att), "n_S": nS}
    print(f"P5 deletion arm: attainable D reads {len(attainable)} (on the survivor {len(on_surv)}, attributed {len(att)}); S reads {nS}")
    for arm in a.deletion.split(","):
        calls = acc.arm_m(f"{W}/{arm}", lab, P)
        calls.pop("__false_moves__")
        tot = collections.Counter(calls.values())
        ra = sum(1 for n in attainable if calls[n] == ("D", "right"))
        fm = tot[("S", "false_move")]
        h1, h2 = ra >= 0.5 * len(attainable), fm <= 0.05 * nS
        per = {f: sum(1 for n, r in lab.items() if r["family"] == f and calls[n] == ("D", "right")) for f in P}
        out[arm] = dict(D_right=tot[("D", "right")], D_right_attainable=ra, false_moves=fm, H5_1=h1, H5_2=h2,
                        families_with_D_right=sum(1 for v in per.values() if v), per_family_D_right=per)
        print(f"  {arm}: D right {tot[('D', 'right')]}; among attainable {ra} = {ra / max(1, len(attainable)):.3f} "
              f"(H5-1 >= 0.50: {'PASS' if h1 else 'FAIL'}); false moves {fm}/{nS} = {fm / nS:.2%} (H5-2 <= 5%: "
              f"{'PASS' if h2 else 'FAIL'}); deleted copies with >= 1 D read right {out[arm]['families_with_D_right']}/{len(P)}")
    for arm in a.control.split(","):
        j = json.load(open(f"{W}/{arm}/classify.json"))["delta"]
        h3 = j["families_false"] <= 3
        out[arm] = dict(classes=j["classes"], families_any=j["families_any"], families_false=j["families_false"], H5_3=h3)
        print(f"  {arm}: flags {j['candidates']} {j['classes']}; families with a false flag {j['families_false']}/20 "
              f"(H5-3 <= 3: {'PASS' if h3 else 'FAIL'})")
    json.dump(out, open(f"{W}/verdict.json", "w"), indent=1)


def p6panel(a):
    """P6 (docs/PREREG_o3_putative_heldout_2026-10-23.md): P4's draw minus P5's 20, plus 40 drawn with random.Random(1) from the
    eligible families P4 did not draw (sorted by family) -> W/panel.json"""
    eligible, _ = eligible_families()
    p4 = json.load(open(f"{RA}/p4/panel.json"))
    assert [p["fam"] for p in p4] == [p["fam"] for p in sorted(random.Random(1).sample(eligible, N_DRAW), key=lambda p: p["fam"])], \
        "P4's draw does not reproduce"
    p5 = {p["fam"] for p in json.load(open(f"{RA}/p5/panel.json"))}
    rest = [p for p in p4 if p["fam"] not in p5]
    pool = [p for p in eligible if p["fam"] not in {q["fam"] for q in p4}]
    new = random.Random(1).sample(pool, 40)
    draw = sorted(rest + new, key=lambda p: p["fam"])
    assert len(rest) == 40 and len({p["fam"] for p in draw}) == 80
    json.dump(draw, open(f"{a.w}/panel.json", "w"), indent=0)
    print(f"P4's draw minus P5: {len(rest)}; pool never drawn: {len(pool)}; drawn from it: {len(new)}; panel {len(draw)}")


def p10panel(a):
    """P10 (docs/PREREG_o3_putative_clip_p10_2026-10-23.md): every eligible family in none of P4's, P5's or P6's panels, sorted
    by family -> W/panel.json"""
    eligible, _ = eligible_families()
    used = set()
    for w in ("p4", "p5", "p6"):
        used |= {p["fam"] for p in json.load(open(f"{RA}/{w}/panel.json"))}
    draw = sorted((p for p in eligible if p["fam"] not in used), key=lambda p: p["fam"])
    assert len(draw) == 58, len(draw)
    json.dump(draw, open(f"{a.w}/panel.json", "w"), indent=0)
    print(f"eligible {len(eligible)}; used by P4-P6 {len(used)}; panel {len(draw)}")


def putprod(a):
    """The port's putative copies as stage-style products (cand.*: one flagged candidate cand_<family>_<k> per putative copy,
    its consensus as the contig) in W/<arm>put, so the control classification steps run on them unchanged. Members come from
    the dump (the merged haplotype flagged putative on that copy with that read count)."""
    X, out = f"{a.w}/{a.arm}", f"{a.w}/{a.arm}put"
    os.makedirs(f"{out}/logs", exist_ok=True)
    rows = tsv(f"{X}/cand.putative.tsv")
    seqs = {}
    name = None
    for ln in open(f"{X}/cand.putative.fa"):
        if ln.startswith(">"):
            name = ln[1:].strip()
            seqs[name] = ""
        else:
            seqs[name] += ln.strip()
    dump = {json.loads(ln)["family"]: json.loads(ln) for ln in open(f"{X}/cand.dump.jsonl")}
    k = collections.Counter()
    with open(f"{out}/cand.candidates.tsv", "w") as ct, open(f"{out}/cand.contigs.fa", "w") as cf, \
            open(f"{out}/cand.clusters.tsv", "w") as cl, open(f"{out}/cand.reads.tsv", "w") as rd, \
            open(f"{out}/cand.families.tsv", "w") as ft:
        ct.write("family\tcandidate\tn_clusters\tn_reads\tflagged\tunion_len\tnearest_locus\td\tn_net\tn_used\n")
        cl.write("family\tcluster\tcandidate\tn_reads\tconsensus_len\tfate\tlinked_to\td\n")
        rd.write("read\tfamily\tcluster\n")
        ft.write("family\tn_net\tn_used\n")
        for r in rows:
            f = r["family"]
            ms = [m for m in dump[f]["merged"] if m.get("putative") and m["slot"] == r["copy"] and len(m["reads"]) == int(r["n_reads"])]
            assert len(ms) == 1, (f, r["candidate"], len(ms))
            c = f"cand_{f}_{k[f]}"
            k[f] += 1
            sq = seqs[r["candidate"]]
            ct.write(f"{f}\t{c}\t1\t{r['n_reads']}\t1\t{len(sq)}\t{r['copy']}\tNA\tNA\tNA\n")
            cf.write(f">{c}\n{sq}\n")
            cl.write(f"{f}\t{f}:p{k[f]}\t{c}\t{r['n_reads']}\t{len(sq)}\tputative\t-\tNA\n")
            for n in ms[0]["reads"]:
                rd.write(f"{n}\t{f}\t{f}:p{k[f]}\n")
        for f in sorted(k):
            ft.write(f"{f}\tNA\tNA\n")
    print(f"{a.arm}: {len(rows)} putative copies in {len(k)} families -> {out}")


def verdict6(a):
    W = a.w
    P = {p["fam"]: p for p in json.load(open(f"{W}/panel.json"))}
    lab = {r["read"]: r for r in tsv(f"{W}/labels.tsv")}
    cls = lambda arm: {r["candidate"].split(":", 1)[1].replace("iso_", "cand_"): r["class"]
                       for r in tsv(f"{W}/{arm}put/candidates.tsv")} if os.path.exists(f"{W}/{arm}put/candidates.tsv") else {}
    out = {}
    # deletion arm: precision (HP1) and recall over the free-slot deleted copies (HP3)
    dcls = cls("pd")
    flags, free = [], []
    k = collections.Counter()
    for ln in open(f"{W}/pd/cand.dump.jsonl"):
        d = json.loads(ln)
        f = d["family"]
        best = None
        for m in d["merged"]:
            nd = sum(1 for r in m["reads"] if lab.get(r, {}).get("role") == "D")
            if m.get("putative"):
                c = f"cand_{f}_{k[f]}"
                k[f] += 1
                flags.append((f, c, len(m["reads"]), 2 * nd > len(m["reads"])))
            if 2 * nd > len(m["reads"]) and len(m["reads"]) >= 6 and (best is None or nd > best[0]):
                best = (nd, m)
        if best and not best[1]["certified"] and best[1]["slot"]:
            free.append((f, len(best[1]["reads"]), bool(best[1].get("putative"))))
    # the putative.tsv order = the dump's merged order per family? match by the rows instead (putprod's names)
    rows = tsv(f"{W}/pd/cand.putative.tsv")
    true = 0
    detail = []
    kk = collections.Counter()
    dump = {json.loads(ln)["family"]: json.loads(ln) for ln in open(f"{W}/pd/cand.dump.jsonl")}
    for r in rows:
        f = r["family"]
        m = [m for m in dump[f]["merged"] if m.get("putative") and m["slot"] == r["copy"] and len(m["reads"]) == int(r["n_reads"])][0]
        c = f"cand_{f}_{kk[f]}"
        kk[f] += 1
        nd = sum(1 for x in m["reads"] if lab.get(x, {}).get("role") == "D")
        dmaj = 2 * nd > len(m["reads"])
        cl_ = dcls.get(c, "unclassified")
        t = dmaj or cl_ == "a_haplotype_only"
        true += t
        detail.append(dict(family=f, candidate=c, reads=len(m["reads"]), D_majority=dmaj, cls=cl_, true=bool(t)))
    hp1 = (true / len(rows)) if rows else None
    hp3 = (sum(1 for x in free if x[2]) / len(free)) if free else None
    # control arm: families with a false putative flag (HP2)
    ccls = cls("pc")
    crow = tsv(f"{W}/pc/cand.putative.tsv")
    false_fams = sorted({c.rsplit("_", 1)[0][5:] for c, v in ccls.items() if v in ("b_allele", "c_unmatched", "pri")})
    out = dict(deletion_flags=detail, HP1=hp1, HP1_pass=hp1 is not None and hp1 >= 0.5,
               free_slot=free, HP3=hp3, HP3_pass=hp3 is not None and hp3 >= 0.25,
               control_flags=len(crow), control_classes=dict(collections.Counter(ccls.values())), control_false_families=false_fams,
               HP2_pass=len(false_fams) <= int(0.15 * len(P)))
    out["verdict"] = out["HP1_pass"] and out["HP2_pass"] and out["HP3_pass"]
    json.dump(out, open(f"{W}/verdict6.json", "w"), indent=1)
    print(f"HP1 precision: {true}/{len(rows)} = {hp1}  (>= 0.50: {'PASS' if out['HP1_pass'] else 'FAIL'})")
    for x in detail:
        print("   ", x)
    print(f"HP3 recall: free-slot deleted copies {len(free)}, flagged {sum(1 for x in free if x[2])} = {hp3}  "
          f"(>= 0.25: {'PASS' if out['HP3_pass'] else 'FAIL'}); {free}")
    print(f"HP2 control: putative flags {len(crow)}, classes {out['control_classes']}, families with a false flag "
          f"{len(false_fams)} {false_fams} (<= {int(0.15 * len(P))}: {'PASS' if out['HP2_pass'] else 'FAIL'})")
    print(f"VERDICT (putative tier): {'PASSES' if out['verdict'] else 'DOES NOT PASS'}")


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["panel", "split", "batches", "verdict", "subpanel", "verdict5", "p6panel", "putprod", "verdict6", "p10panel"])
    ap.add_argument("--arm", default="pd")
    ap.add_argument("--w", required=True)
    ap.add_argument("--source", default=f"{RA}/p4", help="subpanel: the draw to sample from")
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--deletion", default="pd,pdcc")
    ap.add_argument("--control", default="pc,pccc")
    a = ap.parse_args(argv)
    {"panel": panel, "split": split, "batches": batches, "verdict": verdict, "subpanel": subpanel, "verdict5": verdict5,
     "p6panel": p6panel, "putprod": putprod, "verdict6": verdict6, "p10panel": p10panel}[a.cmd](a)


if __name__ == "__main__":
    main(sys.argv[1:])
