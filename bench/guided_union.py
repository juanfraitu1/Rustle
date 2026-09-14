#!/usr/bin/env python3
"""Prereg Addendum S: S1 guided candidates from EITHER finder (transcript hits ∪ gene-body chains) with hybrid width;
S2 subfamily trees on EXON and INTRON sequence classes (reference-projected alignment + IQ-TREE).

usage: guided_union.py <workdir> <refseq gff> <genes.gff.gz> <genome.fa> <iqtree3>
Reuses the definitions of bench/guided_genebody.py (finders, classification, metrics), bench/guided_tree_width.py
(literature groups, newick, clade scoring) and bench/guided_projection_tree.py (projection alignment, IQ-TREE).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
IQTREE = sys.argv[5]
_argv = list(sys.argv)
exec(compile(open(os.path.join(HERE, "guided_genebody.py")).read().split("# reference: all truth records of each family")[0],
             "guided_genebody.py(defs)", "exec"))
_tw = open(os.path.join(HERE, "guided_tree_width.py")).read()
exec(compile(_tw[_tw.index("# ---------------- literature groups ----------------"):_tw.index("def build_tree")],
             "guided_tree_width.py(groups+newick)", "exec"))
exec(compile(_tw[_tw.index("def tree_scores"):_tw.index("# ---------------- P1: reference + leave-out ----------------")],
             "guided_tree_width.py(scores)", "exec"))
_pt = open(os.path.join(HERE, "guided_projection_tree.py")).read()
PD = f"{W}/tree_union"
os.makedirs(PD, exist_ok=True)
exec(compile(_pt[_pt.index("def rc(s):"):_pt.index("out = []")], "guided_projection_tree.py(defs)", "exec"))
sys.argv = _argv
import time


def tx_exon_blocks(h):
    blocks, pos, cur = [], h["ts"], h["ts"]
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", h["cg"]):
        n = int(n)
        if op in "M=XD":
            pos += n
        elif op == "N":
            if pos > cur:
                blocks.append((cur, pos))
            pos += n
            cur = pos
    if pos > cur:
        blocks.append((cur, pos))
    return blocks


def union_candidates(seeds):
    blocked = [(rec[s]["chrom"], rec[s]["start0"], rec[s]["end"]) for s in seeds]
    items = [("tx", h, h["s"], h["e"]) for h in M0_hits if h["q"] in seeds] + \
            [("chain", c, c["s"], c["e"]) for c in G1_hits if c["q"] in seeds]
    items = [x for x in items if not any(ch == x[1]["chrom"] and ov(s, e, x[2], x[3]) > 0 for ch, s, e in blocked)]
    by = collections.defaultdict(list)
    for x in items:
        by[x[1]["chrom"]].append(x)
    out = []
    for chrom, xs in by.items():
        xs.sort(key=lambda x: x[2])
        cl, end = [], -1
        for x in xs + [None]:
            if x is None or (cl and x[2] >= end):
                txs = [y for y in cl if y[0] == "tx"]
                chs = [y for y in cl if y[0] == "chain"]
                btx = max(txs, key=lambda y: (y[1]["nm"], -y[2]))[1] if txs else None
                bch = max(chs, key=lambda y: (y[1]["nm"], -y[2]))[1] if chs else None
                lead = btx or bch
                s, e = (btx["s"], btx["e"]) if btx else (bch["s"], bch["e"])
                out.append({"chrom": chrom, "s": s, "e": e, "clip": (s, e), "family": rec[lead["q"]]["family"],
                            "seed": lead["q"], "best": bch, "btx": btx, "source": "+".join(sorted({y[0] for y in cl}))})
                cl, end = [], -1
            if x is not None:
                cl.append(x)
                end = max(end, x[3])
    return out


def member_seqs_record(name):
    t = rec[name]
    ex = [(max(s, t["start0"]), min(e, t["end"])) for s, e in union_exons[name] if min(e, t["end"]) > max(s, t["start0"])]
    exon = "".join(genome.fetch(t["chrom"], s, e) for s, e in merge(ex)).upper()
    exon = exon.translate(COMP)[::-1] if t["strand"] == "-" else exon
    intron = seq_minus(t["chrom"], t["start0"], t["end"], union_exons[name], t["strand"])
    return exon, intron


def member_seqs_candidate(c):
    if c["btx"] is not None:
        b = tx_exon_blocks(c["btx"])
        exon = "".join(genome.fetch(c["chrom"], s, e) for s, e in b).upper()
        exon = exon.translate(COMP)[::-1] if c["btx"]["strand"] == "-" else exon
    else:
        b = project_exons(c["best"])
        exon = "".join(genome.fetch(c["chrom"], s, e) for s, e in b).upper()
        exon = exon.translate(COMP)[::-1] if c["best"]["strand"] == "-" else exon
    if c["best"] is not None:
        intron = seq_minus(c["chrom"], c["best"]["xs"], c["best"]["xe"], project_exons(c["best"]), c["best"]["strand"])
    else:
        intron = seq_minus(c["chrom"], c["btx"]["s"], c["btx"]["e"], tx_exon_blocks(c["btx"]), c["btx"]["strand"])
    return exon, intron


# ---------------- S1 breadth/width + S2 member sets ----------------
rows, runs = [], []
for fam in families:
    ex_seqs, in_seqs = {}, {}
    for t in truth:
        if t["family"] == fam:
            ex_seqs[t["name"]], in_seqs[t["name"]] = member_seqs_record(t["name"])
    runs.append((f"ref_{fam}", fam, {k: k for k in ex_seqs}, ex_seqs, in_seqs))

for level in ("half", "keep1"):
    for rep in range(REPS):
        seeds, hidden = set(), set()
        for fi, fam in enumerate(families):
            names = sorted(t["name"] for t in truth if t["family"] == fam)
            random.Random(1000 * rep + fi).shuffle(names)
            k = math.ceil(len(names) / 2) if level == "half" else 1
            seeds.update(names[:k])
            hidden.update(names[k:])
        arms = {"M0": candidates(M0_hits, seeds, "nm", lambda b: (b["s"], b["e"])),
                "G1": candidates(G1_hits, seeds, "nm", lambda b: (b["s"], b["e"])),
                "U": union_candidates(seeds)}
        for arm, cands in arms.items():
            rp = classify(cands, hidden)
            for fam in families + ["ALL"]:
                H_ = [n for n in hidden if fam in ("ALL", rec[n]["family"])]
                C_ = [c for c in cands if fam in ("ALL", c["family"])]
                tp = sum(1 for n in H_ if rp[n] and rp[n]["family"] == rec[n]["family"])
                good = sum(1 for c in C_ if c["class"] == "hidden" and rec[c["match"]]["family"] == c["family"])
                named = sum(1 for c in C_ if c["class"] == "other_gene" and c["named"])
                cross = sum(1 for c in C_ if c["class"] == "hidden" and rec[c["match"]]["family"] != c["family"])
                unnamed = [",".join(c["genes"]) for c in C_ if c["class"] == "other_gene" and not c["named"]]
                unann = sum(1 for c in C_ if c["class"] == "unannotated")
                ip = [rp[n]["family"] if rp[n] else f"missed:{n}" for n in H_] + [c["family"] for c in C_ if c["class"] != "hidden"]
                it = [rec[n]["family"] for n in H_] + [f"nonmember:{i}" for i, c in enumerate(C_) if c["class"] != "hidden"]
                ps, pp = pairwise(ip, it)
                br, bp, _, _ = bipartite(ip, it) if ip else (float("nan"),) * 4
                wd = [width_stats(rec[n], rp[n]["s"], rp[n]["e"]) for n in H_ if rp[n] and rp[n]["family"] == rec[n]["family"]]
                rows.append({"level": level, "rep": rep, "arm": arm, "family": fam, "hidden": len(H_), "candidates": len(C_),
                             "sens": tp / len(H_) if H_ else float("nan"), "prec": good / len(C_) if C_ else float("nan"),
                             "prec_named": (good + named) / len(C_) if C_ else float("nan"), "cross": cross,
                             "unnamed": len(unnamed), "unannotated": unann, "unnamed_list": ";".join(unnamed),
                             "pair_sens": ps, "pair_prec": pp, "bip_R": br, "bip_P": bp, "n_width": len(wd),
                             "jaccard": statistics.median([w[0] for w in wd]) if wd else float("nan"),
                             "truncated": sum(w[3] for w in wd), "overextended": sum(w[4] for w in wd),
                             "sources": ";".join(sorted(collections.Counter(c.get("source", arm) for c in C_).elements())) if arm == "U" else ""})
            if arm != "U":
                continue
            for fam in families:
                ex_seqs, in_seqs, label_of = {}, {}, {}
                for s in sorted(seeds):
                    if rec[s]["family"] == fam:
                        ex_seqs[s], in_seqs[s] = member_seqs_record(s)
                        label_of[s] = s
                for i, c in enumerate(cands):
                    if c["family"] != fam:
                        continue
                    key = f"cand{i}"
                    ex_seqs[key], in_seqs[key] = member_seqs_candidate(c)
                    if c["class"] == "hidden" and rec[c["match"]]["family"] == fam:
                        label_of[key] = c["match"]
                runs.append((f"{level}_{rep}_{fam}", fam, label_of, ex_seqs, in_seqs))

# ---------------- S2 trees ----------------
tree_res = []
for tag, fam, label_of, ex_seqs, in_seqs in runs:
    for cls, seqs in (("exon", ex_seqs), ("intron", in_seqs)):
        seqs = {k: v for k, v in seqs.items() if len(v) >= 100}
        if len(seqs) < 4:
            tree_res.append((tag, fam, cls, None, len(seqs)))
            continue
        t0 = time.time()
        aln, ref, Lr, kept = project_alignment(f"{tag}_{cls}", seqs)
        # guard: a member that never aligned to the reference is an all-gap row IQ-TREE rejects; drop and report it
        recs_, cur_ = {}, None
        for line in open(aln):
            if line.startswith(">"):
                cur_ = line[1:].strip()
                recs_[cur_] = ""
            else:
                recs_[cur_] += line.strip()
        empty = [k for k, v in recs_.items() if not v.strip("-")]
        if empty:
            with open(aln, "w") as fh:
                for k, v in recs_.items():
                    if k not in empty:
                        fh.write(f">{k}\n{v}\n")
            print(f"[tree] {tag} {cls}: dropped unaligned members {empty}", file=sys.stderr, flush=True)
            seqs = {k: v for k, v in seqs.items() if k not in empty}
        if len(seqs) < 4 or kept == 0:
            tree_res.append((tag, fam, cls, None, len(seqs)))
            continue
        tree, t_tree, model = iqtree(f"{tag}_{cls}", aln)
        lab = {k: v for k, v in label_of.items() if k in seqs}
        tree_res.append((tag, fam, cls, (tree_scores(tree, fam, lab), ref, Lr, kept, model, time.time() - t0), len(seqs)))
        print(f"[tree] {tag} {cls}: {len(seqs)} members, reference {lab.get(ref, ref)} {Lr} bp -> {kept} columns, {model}",
              file=sys.stderr, flush=True)


# ---------------- report ----------------
def ms(v):
    v = [x for x in v if not (isinstance(x, float) and math.isnan(x))]
    return "nan" if not v else (f"{statistics.mean(v):.3f}±{statistics.stdev(v):.3f}" if len(v) > 1 else f"{v[0]:.3f}")


print("## S1 — breadth and width: M0 (transcript), G1 (gene body), U (either finder; width = transcript hit else chain)")
for level in ("half", "keep1"):
    for fam in families + ["ALL"]:
        for arm in ("M0", "G1", "U"):
            R = [r for r in rows if r["level"] == level and r["family"] == fam and r["arm"] == arm]
            g = lambda k: ms([r[k] for r in R])
            print(f"{level:5s} {fam:6s} {arm:3s} hid {g('hidden')} cand {g('candidates')} | sens {g('sens')} prec {g('prec')} "
                  f"named {g('prec_named')} | cross {g('cross')} unnamed {g('unnamed')} unannot {g('unannotated')} | pair "
                  f"{g('pair_sens')}/{g('pair_prec')} bip {g('bip_R')}/{g('bip_P')} | W n={g('n_width')} J {g('jaccard')} "
                  f"trunc {g('truncated')} overext {g('overextended')}")
        U = [r for r in rows if r["level"] == level and r["family"] == fam and r["arm"] == "U"]
        if fam != "ALL":
            srcs = collections.Counter(s for r in U for s in r["sources"].split(";") if s)
            print(f"            U candidate sources over reps: {dict(srcs)}")
            names = sorted({x for r in U for x in r["unnamed_list"].split(";") if x})
            if names:
                print(f"            U unnamed other-gene candidates: {names}")
print("\nS1 reading (keep 1, sensitivity >= max(M0, G1) - 0.02, named precision >= 0.95, 0 cross-family):")
for fam in families:
    m = {arm: statistics.mean([r["sens"] for r in rows if r["level"] == "keep1" and r["family"] == fam and r["arm"] == arm])
         for arm in ("M0", "G1", "U")}
    pn = statistics.mean([r["prec_named"] for r in rows if r["level"] == "keep1" and r["family"] == fam and r["arm"] == "U"])
    cr = sum(r["cross"] for r in rows if r["level"] == "keep1" and r["family"] == fam and r["arm"] == "U")
    ok = m["U"] >= max(m["M0"], m["G1"]) - 0.02 and pn >= 0.95 and cr == 0
    print(f"  {fam}: U {m['U']:.3f} vs M0 {m['M0']:.3f} / G1 {m['G1']:.3f}; named precision {pn:.3f}; cross-family {cr} -> {'PASS' if ok else 'FAIL'}")

print("\n## S2 — clades per sequence class (SH-aLRT > 75)")
summary = collections.defaultdict(collections.Counter)
by_run = collections.defaultdict(dict)
for tag, fam, cls, r, nmem in tree_res:
    if r is None:
        print(f"{tag} {cls}: {nmem} members -> not treed")
        continue
    (res, positional, part, mid), ref, Lr, kept, model, secs = r
    by_run[(tag, fam)][cls] = (res, positional)
    groups = "; ".join(f"{g}: {v[0]}" + (f" ({v[1]:.0f}/{v[2]:.0f})" if v[1] is not None else "") for g, v in res.items())
    extra = "" if positional is None else f" | positional supported: {positional}"
    print(f"{tag} {cls}: {nmem} members, {kept}/{Lr} columns, {model} | {groups}{extra}")
for (tag, fam), d in by_run.items():
    kind = "reference" if tag.startswith("ref") else "leave-out"
    groups = set().union(*(set(v[0]) for v in d.values()))
    for g in groups:
        calls = {cls: d[cls][0].get(g, ("absent",))[0] for cls in d}
        if all(c == "absent" for c in calls.values()):
            continue
        for cls, c in calls.items():
            summary[(fam, kind, g, cls)][c] += 1
        summary[(fam, kind, g, "either")]["RECOVERED" if any(c == "RECOVERED" for c in calls.values()) else "not"] += 1
    for cls in d:
        if d[cls][1] is not None:
            summary[(fam, kind, "positional split", cls)]["WRONG" if d[cls][1] else "CORRECT"] += 1
print("\n## S2 summary (counts over runs)")
for k in sorted(summary):
    print(f"{k[0]:6s} {k[1]:9s} {k[2]:28s} {k[3]:6s} {dict(summary[k])}")
