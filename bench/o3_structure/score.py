#!/usr/bin/env python3
"""Parity and scoring for the O3 structure projection (docs/PREREG_o3_structure_projection_2026-10-14.md). Imports the frozen reference
implementation `proto.py` (never edited after the freeze) for the rule and the truth.

    score.py parity WORLD VARIANT PREFIX      # P0: the Rust boundaries (PREFIX.structure.tsv) vs the prototype's, per flagged candidate
    score.py score WORLD VARIANT PREFIX [--rows OUT.tsv]   # the registered scoring of the Rust product against the truth

PREFIX is the candidate_copies --out prefix (its candidates.tsv / contigs.fa / structure.tsv).
"""
import argparse
import collections
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import proto as P  # noqa: E402

TOL = 5


def rust_bounds(prefix):
    out = {}
    for r in csv.DictReader(open(f"{prefix}.structure.tsv"), delimiter="\t"):
        out[r["candidate"]] = [] if r["boundaries"] == "-" else [int(x) for x in r["boundaries"].split(",")]
    return out


def proto_bounds(world, prefix):
    """the prototype's boundaries for every flagged candidate of PREFIX (its own genome + segment alignments)"""
    cands = [r for r in csv.DictReader(open(f"{prefix}.candidates.tsv"), delimiter="\t") if r["flagged"] == "1"]
    U = P.read_fa(f"{prefix}.contigs.fa")
    genome = P.read_fa(f"{world}/genome.ref.fa")
    fam_copies = collections.defaultdict(list)
    for row in csv.DictReader(open(f"{world}/P.fam.copies.tsv"), delimiter="\t"):
        fam_copies[row["family_id"]].append((row["chrom"], int(row["start"]), int(row["end"]),
                                             [tuple(map(int, x.split("-"))) for x in row["exons"].split(",")]))
    ghits = collections.defaultdict(list)
    for h in P.paf(P.mm2("-c --eqx -x splice:hq -uf -N 20", f"{world}/genome.ref.splice.mmi", f"{prefix}.contigs.fa")):
        ghits[h["q"]].append(h)
    out = {}
    for r in cands:
        c, fam = r["candidate"], r["family"]
        hs = ghits.get(c, [])
        G = genome[max(hs, key=lambda h: h["m"])["t"]] if hs else ""

        def realign(segs, seq):
            tmp = f"{prefix}.paritytmp.fa"
            with open(tmp, "w") as fh:
                for qs, qe in segs:
                    fh.write(f">{qs}_{qe}\n{seq[qs:qe]}\n")
            got = collections.defaultdict(list)
            for h in P.paf(P.mm2("-c --eqx -x splice:hq -uf -N 20", f"{world}/genome.ref.splice.mmi", tmp)):
                a, b = map(int, h["q"].split("_"))
                got[(a, b)].append(h)
            os.remove(tmp)
            return {sg: got.get(sg, []) for sg in segs}
        proj, _ = P.project(hs, U[c], G, realign, fam_copies[fam])
        out[c] = proj or []
    return out


def cmd_parity(a):
    rb, pb = rust_bounds(a.prefix), proto_bounds(a.world, a.prefix)
    diff = [c for c in pb if rb.get(c) != pb[c]]
    missing = [c for c in pb if c not in rb]
    print(f"{a.world.rstrip('/').split('/')[-1]}\t{a.variant}\tflagged {len(pb)}\tidentical {len(pb) - len(diff)}\tdiffer {len(diff)}"
          f"\tmissing_in_rust {len(missing)}")
    for c in diff:
        print(f"  DIFF {c}\trust {rb.get(c)}\tproto {pb[c]}")
    return 1 if diff else 0


def cmd_score(a):
    world, prefix = a.world, a.prefix
    scen = json.load(open(f"{world}/world.json"))
    kinds = json.load(open(f"{world}/kinds.json"))
    rb = rust_bounds(prefix)
    cands = {r["candidate"]: r for r in csv.DictReader(open(f"{prefix}.candidates.tsv"), delimiter="\t") if r["flagged"] == "1"}
    thits = {}
    for h in P.paf(P.mm2("-c --eqx -x map-hifi -N 5", f"{world}/transcripts.fa", f"{prefix}.contigs.fa")):
        if h["tp"] == "P" and (h["q"] not in thits or h["m"] > thits[h["q"]]["m"]):
            thits[h["q"]] = h
    rows = []
    for c, r in cands.items():
        fam = r["family"]
        copy = fam + "E" if fam + "E" in kinds else None
        th = thits.get(c)
        row = dict(world=world.rstrip("/").split("/")[-1], variant=a.variant, family=fam, scenario=scen[fam], candidate=c,
                   strand="", scored=0, n_true=0, n_proj=len(rb.get(c, [])), rec5=0, prec5=0, rec0=0, prec0=0, correct5=0, exact=0,
                   recon_idcov=0.0, recon_tcov=0.0, note="")
        if copy is None or th is None or th["t"] != copy or th["st"] != "+":
            row["note"] = f"unscored: truth hit {th['t'] if th else '-'}"
            rows.append(row)
            continue
        tb, tlen = P.truth_bounds(world, fam, copy)
        inside = [b for b in tb if th["ts"] + P.MIN_SEG <= b <= th["te"] - P.MIN_SEG]
        mp = P.qmap(th)
        pt = [mp.get(b) for b in rb.get(c, [])]
        for tol, kr, kp in ((TOL, "rec5", "prec5"), (0, "rec0", "prec0")):
            row[kr] = sum(any(p is not None and abs(p - b) <= tol for p in pt) for b in inside)
            row[kp] = sum(any(abs(p - b) <= tol for b in tb) for p in pt if p is not None)
        row.update(scored=1, n_true=len(inside), correct5=int(row["rec5"] == len(inside) and row["prec5"] == len(pt)),
                   exact=int(row["rec0"] == len(inside) and row["prec0"] == len(pt)),
                   recon_idcov=round((th["m"] / th["bl"]) * (th["qe"] - th["qs"]) / th["ql"], 5),
                   recon_tcov=round((th["te"] - th["ts"]) / tlen, 5))
        for ln in open(f"{world}/fam/{fam}/truth.gtf"):
            f = ln.split("\t")
            if len(f) > 8 and f[2] == "transcript" and f'gene_id "{copy}"' in f[8]:
                row["strand"] = f[6]
                break
        rows.append(row)
    if a.rows:
        new = not os.path.exists(a.rows)
        with open(a.rows, "a") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
            if new:
                w.writeheader()
            w.writerows(rows)
    sc = [r for r in rows if r["scored"]]
    print(f"{world.rstrip('/').split('/')[-1]}\t{a.variant}\tflagged {len(rows)}\tscored {len(sc)}\tcorrect±5 {sum(r['correct5'] for r in sc)}"
          f"\texact {sum(r['exact'] for r in sc)}")
    return 0


def cmd_summary(a):
    """the registered bars (S1-S3) and the reported-beside numbers from the rows of `score --rows` and `isocon_cmp.py`"""
    rows = list(csv.DictReader(open(a.rows), delimiter="\t"))
    for r in rows:
        for k in ("scored", "n_true", "n_proj", "rec5", "prec5", "rec0", "prec0", "correct5", "exact"):
            r[k] = int(r[k])
    sc = [r for r in rows if r["scored"]]
    base = [r for r in sc if r["scenario"] != "DIV@0.08"]
    frac = lambda k, n: f"{k}/{n} = {k / n:.3f}" if n else "0/0"
    s1 = sum(r["correct5"] for r in base)
    rec, nt = sum(r["rec5"] for r in base), sum(r["n_true"] for r in base)
    prec, npj = sum(r["prec5"] for r in base), sum(r["n_proj"] for r in base)
    print(f"S1 structure-correct (±5), all but DIV@0.08: {frac(s1, len(base))}  bar >= 0.90  -> {'PASS' if len(base) and s1 / len(base) >= 0.90 else 'FAIL'}")
    s2_ok = True
    for scn in ("LOSS", "GAIN", "DUP", "SHUFFLE", "TRUNC5", "INTRON"):
        x = [r for r in sc if r["scenario"] == scn]
        k = sum(r["correct5"] for r in x)
        ok = bool(x) and k / len(x) >= 0.80
        s2_ok &= ok
        print(f"S2 {scn:8s} {frac(k, len(x))}  bar >= 0.80  -> {'PASS' if ok else 'FAIL'}")
    print(f"S2 overall -> {'PASS' if s2_ok else 'FAIL'}")
    s3 = nt and npj and rec / nt >= 0.97 and prec / npj >= 0.97
    print(f"S3 boundary recall {frac(rec, nt)}, precision {frac(prec, npj)}  bars >= 0.97  -> {'PASS' if s3 else 'FAIL'}")
    print("-- reported beside")
    print(f"exact (±0) structure-correct, all but DIV@0.08: {frac(sum(r['exact'] for r in base), len(base))}")
    for scn in sorted({r["scenario"] for r in rows}):
        x = [r for r in sc if r["scenario"] == scn]
        fams = {(r["world"], r["variant"], r["family"]) for r in rows if r["scenario"] == scn}
        det = {(r["world"], r["variant"], r["family"]) for r in x}
        print(f"  {scn:10s} flagged-families {len(fams):3d}  scored {len(x):3d}  correct±5 {sum(r['correct5'] for r in x):3d}  "
              f"exact {sum(r['exact'] for r in x):3d}  E detected (scored) {len(det)}")
    for st in ("+", "-"):
        x = [r for r in base if r["strand"] == st]
        print(f"  strand {st}: {frac(sum(r['correct5'] for r in x), len(x))}")
    uns = [r for r in rows if not r["scored"]]
    print(f"unscored flagged candidates: {len(uns)}  " + "; ".join(sorted({r['scenario'] + ' ' + r['note'] for r in uns})))
    if a.isocon:
        iso = list(csv.DictReader(open(a.isocon), delimiter="\t"))
        ours = {(r["world"], r["variant"], r["family"]): float(r["recon_idcov"]) for r in sc}
        n = len(iso)
        o_ex = sum(ours.get((r["world"], r["variant"], r["family"]), 0.0) >= 0.999 for r in iso)
        i_ex = sum(int(r["exact"]) for r in iso)
        print(f"reconstruction exact (id x cov >= 0.999) over {n} absent copies: candidate_copies {frac(o_ex, n)}, IsoCon {frac(i_ex, n)}")
        for scn in sorted({r["scenario"] for r in iso}):
            x = [r for r in iso if r["scenario"] == scn]
            print(f"  {scn:10s} ours {sum(ours.get((r['world'], r['variant'], r['family']), 0.0) >= 0.999 for r in x):3d}/{len(x)}  "
                  f"IsoCon {sum(int(r['exact']) for r in x):3d}/{len(x)}")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("summary"); p.add_argument("rows"); p.add_argument("--isocon"); p.set_defaults(func=cmd_summary)
    for name, f in (("parity", cmd_parity), ("score", cmd_score)):
        p = sub.add_parser(name)
        p.add_argument("world"); p.add_argument("variant"); p.add_argument("prefix")
        if name == "score":
            p.add_argument("--rows", help="append one row per flagged candidate")
        p.set_defaults(func=f)
    a = ap.parse_args()
    sys.exit(a.func(a))
