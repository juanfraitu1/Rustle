#!/usr/bin/env python3
"""Held-out test of the tandem member-size rule (docs/PREDICTION_tandem_heldout_2026-11-19.md), any sample.

  tandem_heldout.py members --units UNITS.tsv --annot ANNOT --bam P.bam --gtf stringtie=GTF --gtf rustle_P=GTF [--out FILE]
      TANDEM families, generic: an array with >= 3 protein-coding units on one strand whose annotated unit lengths have a
      coefficient of variation <= 0.2. Members: protein-coding units with >= 5 sense reads (array strand = the strand of most
      units). Per member: the rule (tandem.member_sizes), and for each GTF the same-strand model with the largest exonic overlap
      with the member; each against the member's longest annotated transcript span. Prints within 10% / 20% per method.
  tandem_heldout.py classifier --scenario DIR [--period P]
      The read set-apart step scored on a famsim array with injected reads (array_inject.py names them anti_ / thru_ / runon_):
      per class, the share set apart as antisense or multi-unit; real (simulated transcript) reads wrongly set apart."""
import argparse, collections, csv, os, statistics, sys
import pysam
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from tandem import member_sizes
from undecidable import read_gtf_tx, load_annot

def cmd_members(a):
    tx, genes = load_annot(a.annot)
    span = collections.defaultdict(int)
    for v in tx.values():
        span[v["gene"]] = max(span[v["gene"]], v["exons"][-1][1] - v["exons"][0][0])
    arr = collections.defaultdict(list)
    with open(a.units) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["group"] == "protein_coding":
                arr[r["array"]].append((int(r["start"]), int(r["end"]), r["strand"], r["unit"], r["gene"], r["chrom"]))
    G = {n: read_gtf_tx(g) for n, g in (x.split("=", 1) for x in a.gtf)}
    bam = pysam.AlignmentFile(a.bam)
    res, fams, rows = collections.defaultdict(list), set(), []
    for an, us in arr.items():
        st = collections.Counter(u[2] for u in us).most_common(1)[0][0]
        us = sorted(u for u in us if u[2] == st)
        lens = [span[u[3]] or (u[1] - u[0]) for u in us]
        if len(us) < 3 or statistics.pstdev(lens) / statistics.mean(lens) > 0.2:
            continue
        chrom = us[0][5]
        units = [(u[0], u[1]) for u in us]
        P = statistics.median(b[0] - a_[0] for a_, b in zip(units, units[1:]))
        reads = [(r.reference_start, r.reference_end, "-" if r.is_reverse else "+", any(o == 3 for o, _ in r.cigartuples))
                 for r in bam.fetch(chrom, max(0, int(units[0][0] - P)), int(units[-1][1] + P))
                 if not (r.is_secondary or r.is_supplementary or r.is_unmapped)]
        est, _ = member_sizes(reads, units, P, st)
        for (s0, e0, _, gid, name, _), m, truth in zip(us, est, lens):
            if m["source"] == "none" or m["reads"] < 5:
                continue
            row = {"rule": m["size"]}
            for n, g in G.items():
                best = None
                for v in g.values():
                    if v["chrom"] != chrom or v["strand"] != st:
                        continue
                    ov = sum(max(0, min(y, e0) - max(x, s0)) for x, y in v["exons"])
                    if ov > 0 and (best is None or ov > best[0]):
                        best = (ov, v["exons"][-1][1] - v["exons"][0][0])
                row[n] = best[1] if best else 0
            fams.add(an)
            rows.append((an, name, truth, m["reads"], m["source"], row))
            for k, v in row.items():
                res[k].append(abs(v - truth) / truth)
    n = len(rows)
    L = [f"tandem families (generic definition) with a member of >= 5 sense reads: {len(fams)}; members {n}"]
    for k in ["rule"] + list(G):
        e = res[k]
        L.append(f"  {k}: within 10% {sum(x <= 0.1 for x in e)}/{n}, within 20% {sum(x <= 0.2 for x in e)}/{n} "
                 f"({sum(x <= 0.2 for x in e) / n if n else 0:.2f}), median error {statistics.median(e) if e else 0:.3f}")
    L.append("members: array gene | annotated | reads (source) | " + " | ".join(["rule"] + list(G)))
    for an, name, t, r, src, row in sorted(rows, key=lambda x: x[:5]):
        L.append(f"  {an.split('@')[0]} {name} | {t} | {r} ({src}) | " + " | ".join(str(row[k]) for k in ["rule"] + list(G)))
    out = "\n".join(L) + "\n"
    if a.out:
        open(a.out, "w").write(out)
    print(out, end="")

def cmd_classifier(a):
    d = a.scenario
    cp = [r for r in csv.DictReader(open(f"{d}/copies.tsv"), delimiter="\t") if r["copy_id"].startswith("U")]
    units = sorted((int(r["start"]) - 1, int(r["end"])) for r in cp)
    P = a.period or statistics.median(b[0] - x[0] for x, b in zip(units, units[1:]))
    strand = collections.Counter(r["strand"] for r in cp).most_common(1)[0][0]
    c = collections.defaultdict(collections.Counter)
    with pysam.AlignmentFile(f"{d}/reads.bam") as fh:
        for r in fh.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            q = r.query_name
            cls = ("antisense" if q.startswith("anti_") else "read_through" if q.startswith("thru_") else "run_on" if q.startswith("runon_")
                   else "transcript" if q.split("|")[0].startswith("U") else "decoy gene")
            rs = "-" if r.is_reverse else "+"
            verdict = "set apart: antisense" if rs != strand else "set apart: multi-unit" if r.reference_end - r.reference_start > P else "kept"
            c[cls][verdict] += 1
    for cls in ("transcript", "antisense", "read_through", "run_on", "decoy gene"):
        n = sum(c[cls].values())
        print(f"{cls}: {n} primaries; " + ", ".join(f"{k} {v}" for k, v in sorted(c[cls].items())))

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("members"); p.add_argument("--units", required=True); p.add_argument("--annot", required=True)
    p.add_argument("--bam", required=True); p.add_argument("--gtf", action="append", default=[]); p.add_argument("--out")
    p = sp.add_parser("classifier"); p.add_argument("--scenario", required=True); p.add_argument("--period", type=float)
    a = ap.parse_args(argv)
    dict(members=cmd_members, classifier=cmd_classifier)[a.cmd](a)

if __name__ == "__main__":
    main()
