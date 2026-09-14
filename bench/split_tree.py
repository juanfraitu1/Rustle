#!/usr/bin/env python3
"""Prereg Addendum AA: subfamilies from variation-graph bubbles. Each informative biallelic gap-free column of a
reference-projected family alignment is a bubble that splits the copies in two; a split is KEPT iff its bubble support
exceeds the support of every split incompatible with it. Kept splits are pairwise compatible, so they display as one
unrooted tree (Buneman 1971) with no identity cut, support threshold or substitution model. Literature subfamilies
recovered by the kept splits are compared with the IQ-TREE calls (SH-aLRT > 75) of the same alignments.

usage: split_tree.py --run lit/guided_t --run lit/amy_lo_t
Each run dir holds t.out, t.candidates.tsv, truth.tsv and tree_t/{tag}_{exon,intron}.{proj.fa,treefile} from
`bench/guided_pipeline.py`.
"""
import argparse
import collections
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import guided_pipeline as gp  # noqa: E402

LINE = re.compile(r"^(\S+) (exon|intron): (\d+) leaves, reference (\S+), (\d+) columns(?:, dropped \[.*?\])? \| (.*?)"
                  r"(?: \| positional supported: (True|False))?$")


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            seqs[name] = []
        elif name is not None:
            seqs[name].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


def bubbles(aln):
    """Split -> support over informative biallelic gap-free columns; a split is keyed by the side holding min(leaf)."""
    names = sorted(aln)
    first = names[0]
    sup = collections.Counter()
    for col in zip(*(aln[n] for n in names)):
        if any(c not in "ACGT" for c in col):
            continue
        cnt = collections.Counter(col)
        if len(cnt) != 2 or min(cnt.values()) < 2:
            continue
        state = col[0]  # names[0] == first
        sup[frozenset(n for n, c in zip(names, col) if c == state)] += 1
    return sup, set(names), first


def incompatible(a, c, leaves):
    b, d = leaves - a, leaves - c
    return bool(a & c) and bool(a & d) and bool(b & c) and bool(b & d)


def dominance(sup, leaves):
    splits = list(sup)
    conflicts = {s: [t for t in splits if t is not s and incompatible(s, t, leaves)] for s in splits}
    kept = [s for s in splits if all(sup[s] > sup[t] for t in conflicts[s])]
    return kept, conflicts


def lit_conflicting(sides, label_of, groups):
    """Number of splits whose restriction to truth-labelled leaves is incompatible with a present literature group."""
    truth = {m for m in label_of}
    U = {label_of[m] for m in truth}
    present = [U & G for G in groups.values() if 2 <= len(U & G) < len(U)]
    n = 0
    for side in sides:
        ts = {label_of[m] for m in side if m in truth}
        if any(incompatible(ts, P, U) for P in present):
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="append", required=True)
    a = ap.parse_args()
    summary = collections.defaultdict(collections.Counter)
    totals = collections.Counter()
    gate_fail = []
    fam_conflict = collections.defaultdict(list)
    for W in a.run:
        truth = list(csv.DictReader(open(f"{W}/truth.tsv"), delimiter="\t"))
        fam_of = {t["name"]: t["family"] for t in truth}
        groups = gp.literature_groups(truth)
        positional = {"TBC1D3": {t["name"] for t in truth if t["family"] == "TBC1D3" and t["level1"] == "cluster1"}}
        cands = collections.defaultdict(list)
        for r in csv.DictReader(open(f"{W}/t.candidates.tsv"), delimiter="\t"):
            cands[(r["level"], r["rep"])].append(r)
        printed = {}
        for line in open(f"{W}/t.out"):
            m = LINE.match(line.rstrip("\n"))
            if m:
                calls = dict(x.rsplit(": ", 1) for x in m.group(6).split("; "))
                pos = None if m.group(7) is None else m.group(7) == "True"
                printed[(m.group(1), m.group(2))] = (calls, pos)
        per_run = collections.defaultdict(dict)
        for (tag, cls), (calls, pos) in sorted(printed.items()):
            fam = tag.split("_")[-1]
            base = f"{W}/tree_t/{tag}_{cls}"
            aln = read_fasta(f"{base}.proj.fa")
            if tag.startswith("ref_"):
                label_of = {n: n for n in aln if fam_of.get(n) == fam}
            else:
                level, rep = tag.split("_")[:2]
                label_of = {}
                for n in aln:
                    if n.startswith("cand"):
                        r = cands[(level, rep)][int(n[4:])]
                        if r["class"] == "hidden" and fam_of.get(r["match_or_genes"]) == fam and r["family"] == fam:
                            label_of[n] = r["match_or_genes"]
                    elif fam_of.get(n) == fam:
                        label_of[n] = n
            iq_splits, allleaves = gp.parse_newick(open(f"{base}.treefile").read())
            iq_res, iq_pos = gp.clade_calls(iq_splits, allleaves, label_of, groups[fam], positional.get(fam))
            if iq_res != calls or iq_pos != pos:
                gate_fail.append((W, tag, cls, calls, iq_res, pos, iq_pos))
            sup, leaves, _ = bubbles(aln)
            if len(leaves) < 4 or not sup:
                dt_res = {g: ("absent" if iq_res[g] == "absent" else "not treed") for g in groups[fam]}
                dt_pos, kept, cidx = None, [], float("nan")
            else:
                kept, conflicts = dominance(sup, leaves)
                dt_res, dt_pos = gp.clade_calls([(s, 100.0, 100.0) for s in kept], leaves, label_of, groups[fam],
                                                positional.get(fam))
                keptset = set(kept)
                nb = sum(sup.values())
                cidx = sum(sup[s] for s in sup if any(t in keptset for t in conflicts[s])) / nb
            iq_sup = [s for s, sh, _ in iq_splits if sh > 75]
            lc_dt, lc_iq = lit_conflicting(kept, label_of, groups[fam]), lit_conflicting(iq_sup, label_of, groups[fam])
            totals["dt_lit_conflicting"] += lc_dt
            totals["iq_lit_conflicting"] += lc_iq
            totals["dt_kept"] += len(kept)
            totals["iq_supported"] += len(iq_sup)
            if tag.startswith("ref_"):
                fam_conflict[fam].append((cls, round(cidx, 3), len(sup), sum(sup.values())))
            per_run[tag][cls] = (iq_res, dt_res, iq_pos, dt_pos)
            print(f"{tag:18s} {cls:6s} leaves {len(leaves):2d} bubbles {sum(sup.values()):5d} ({len(sup):3d} splits) kept {len(kept):3d} "
                  f"conflict-index {cidx:.3f} | IQ supported {len(iq_sup):2d} | lit-conflicting DT {lc_dt} IQ {lc_iq} | "
                  + "; ".join(f"{g}: IQ {iq_res[g]} / DT {dt_res[g]}" for g in groups[fam])
                  + ("" if iq_pos is None else f" | positional IQ {iq_pos} DT {dt_pos}"))
        for tag, bycls in per_run.items():
            fam = tag.split("_")[-1]
            kind = "reference" if tag.startswith("ref_") else "leave-out"
            for g in groups[fam]:
                if all(v[0][g] == "absent" for v in bycls.values()):
                    continue
                for meth, k in (("IQ", 0), ("DT", 1)):
                    for cls, v in bycls.items():
                        summary[(fam, kind, g, cls, meth)][v[k][g]] += 1
                    rec = any(v[k][g] == "RECOVERED" for v in bycls.values())
                    summary[(fam, kind, g, "either", meth)]["RECOVERED" if rec else "not"] += 1
                    totals[f"{meth}_either_recovered"] += rec
            for cls, v in bycls.items():
                for meth, k in (("IQ", 2), ("DT", 3)):
                    if v[k] is not None:
                        summary[(fam, kind, "positional split", cls, meth)]["kept/supported" if v[k] else "absent"] += 1

    print(f"\nGATE (IQ-TREE calls re-derived == t.out): {'PASS' if not gate_fail else 'FAIL'}")
    for x in gate_fail[:10]:
        print("  mismatch:", x)
    if gate_fail:
        return
    print("\n## summary (runs)")
    for k in sorted(summary):
        print(f"{k[0]:6s} {k[1]:9s} {k[2]:28s} {k[3]:6s} {k[4]} {dict(summary[k])}")
    print("\n## conflict index on reference alignments (class, index, distinct splits, bubbles)")
    for fam, v in fam_conflict.items():
        print(f"  {fam}: {v}")
    print(f"\nTOTALS: either-class recoveries IQ {totals['IQ_either_recovered']} DT {totals['DT_either_recovered']}; "
          f"literature-conflicting splits IQ {totals['iq_lit_conflicting']} DT {totals['dt_lit_conflicting']}; "
          f"splits IQ supported {totals['iq_supported']} DT kept {totals['dt_kept']}")
    ok = totals["DT_either_recovered"] >= totals["IQ_either_recovered"] and totals["dt_lit_conflicting"] <= totals["iq_lit_conflicting"]
    print(f"READING AA: {'SUPPORTED' if ok else 'NOT SUPPORTED'}")


if __name__ == "__main__":
    main()
