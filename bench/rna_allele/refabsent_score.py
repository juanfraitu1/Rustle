#!/usr/bin/env python3
"""Amendment 10 scoring (docs/PREREG_rna_allele_haplotype_count_2026-10-01.md): the chain's candidates (control_test.py classify, run with
--l and --w on the refabsent work dir) against the truth of real reference-absent loci (refabsent_truth.py).

  D1  every expressed beyond-delta locus with >= 20 reads gets a FLAG (candidate with >= 2 transcripts) whose best haplotype hit
      (identity x coverage >= 0.999) overlaps the locus; the 6-13-read and 2-read loci reported.
  D2  no expressed within-delta locus gets a new-copy flag.
  D3  families (of the 34) with a >= 2-transcript candidate of class b / c <= 20%.
  O2 side: reads of each expressed beyond-delta locus in arm R0 (`_pri` placement, median de) and in arm C (candidate), when C.bam exists.

    refabsent_score.py --w /mnt/linuxdisk/tmp/rna_allele/refabsent
"""
import argparse
import collections
import csv
import json
import os
import re
import statistics

import pysam

TRUTH = "/mnt/linuxdisk/tmp/rna_allele"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--w", required=True)
    a = ap.parse_args(argv)
    W = a.w
    B = {r["locus"]: r for r in csv.DictReader(open(f"{W}/bonly.tsv"), delimiter="\t")}
    ex = {r["locus"]: int(r["n_reads"]) for r in csv.DictReader(open(f"{W}/bonly_expressed.tsv"), delimiter="\t")}
    lab = {r["read"]: r for r in csv.DictReader(open(f"{W}/labels.tsv"), delimiter="\t")}
    alias = {}
    for h in ("mat", "pat"):
        for acc, num, _ in csv.reader(open(f"{TRUTH}/{h}.len.tsv"), delimiter="\t"):
            alias[(h, num)] = acc

    def hapname(name):
        m = re.fullmatch(r"chr(\w+?)_(mat|pat)_hsa[^_]*", name)
        return (m.group(2), alias[(m.group(2), m.group(1))]) if m else (None, name)
    # reads of each truth locus (as refabsent_truth.express)
    recs = collections.defaultdict(list)
    for h in ("mat", "pat"):
        for rd in pysam.AlignmentFile(f"{W}/reads.{h}.bam").fetch(until_eof=True):
            if rd.is_unmapped or rd.is_supplementary:
                continue
            recs[rd.query_name].append((rd.get_tag("AS") if rd.has_tag("AS") else 0, h, hapname(rd.reference_name)[1],
                                        rd.reference_start, rd.reference_end))
    on = collections.defaultdict(set)
    for n, rs in recs.items():
        rs.sort(key=lambda r: -r[0])
        if len(rs) > 1 and rs[1][0] >= 0.98 * rs[0][0]:
            continue
        _, h, acc, s, e = rs[0]
        for L, b in B.items():
            if b["family"] == lab[n]["family"] and b["hap"] == h and b["chrom"] == acc and s < int(b["end"]) and int(b["start"]) < e:
                on[L].add(n)
    # candidates and their best haplotype hit over their contigs
    comp = {r["contig"]: r["component"] for r in csv.DictReader(open(f"{W}/merge/components.tsv"), delimiter="\t")}
    cand = {r["candidate"]: r for r in csv.DictReader(open(f"{W}/candidates.tsv"), delimiter="\t")}
    hb = {}
    for h in ("mat", "pat"):
        for ln in open(f"{W}/contigs_L.{h}.paf"):
            f = ln.rstrip("\n").split("\t")
            hh, acc = hapname(f[5])
            if not hh:
                continue
            s = int(f[9]) / max(1, int(f[10])) * (int(f[3]) - int(f[2])) / max(1, int(f[1]))
            cid = comp.get(f[0])
            if cid and s >= 0.999 and (cid not in hb or s > hb[cid][0]):
                hb[cid] = (s, hh, acc, int(f[7]), int(f[8]))

    cmap = {c: comp[r["best_contig"]] for c, r in cand.items()}     # candidates.tsv name -> components.tsv id

    def matches(cid, L):
        b, t = hb.get(cmap[cid]), B[L]
        return bool(b) and b[1] == t["hap"] and b[2] == t["chrom"] and b[3] < int(t["end"]) and int(t["start"]) < b[4]
    flags = {c for c, r in cand.items() if int(r["n_contigs"]) >= 2}
    print(f"candidates {len(cand)}; flags (>= 2 transcripts) {len(flags)} in {len({cand[c]['family'] for c in flags})} families; "
          f"classes of flags {dict(collections.Counter(cand[c]['class'] for c in flags))}")
    # D1 / D2 tables
    expressed = [L for L, n in ex.items() if n >= 3]
    rows = []
    for L in sorted(expressed, key=lambda L: -ex[L]):
        fl = [c for c in flags if matches(c, L)]
        single = [c for c in cand if c not in flags and matches(c, L)]
        rows.append((L, B[L]["class"], ex[L], fl, single))
        print(f"  {L:14s} {B[L]['class']:20s} reads {ex[L]:4d}  flags matching {len(fl)} {fl}  single-transcript matches {len(single)}")
    d1 = [(L, k, n, fl, s) for L, k, n, fl, s in rows if k == "absent_beyond_delta" and n >= 20]
    print(f"D1: beyond-delta loci with >= 20 reads {len(d1)}: flagged {sum(1 for x in d1 if x[3])} -> "
          f"{'PASSES' if d1 and all(x[3] for x in d1) else 'FAILS'}; beyond-delta loci with 3-19 reads: "
          f"{[(L, n, len(fl)) for L, k, n, fl, s in rows if k == 'absent_beyond_delta' and n < 20]}; with 2 reads: "
          f"{[(L, len([c for c in flags if matches(c, L)])) for L, n in ex.items() if n == 2 and B[L]['class'] == 'absent_beyond_delta']}")
    d2 = [(L, fl) for L, k, n, fl, s in rows if k == "absent_within_delta" and fl]
    print(f"D2: expressed within-delta loci {sum(1 for r in rows if r[1] == 'absent_within_delta')}, flagged {len(d2)} {d2} -> "
          f"{'HOLDS' if not d2 else 'FAILS'}")
    fams = {r["family"] for r in csv.DictReader(open(f"{W}/bonly.tsv"), delimiter="\t")}
    false_f = {cand[c]["family"] for c in flags if cand[c]["class"] in ("b_allele", "c_unmatched", "pri")}
    print(f"D3: families with a false flag {len(false_f)}/{len(fams)} = {len(false_f) / len(fams):.1%} (bar 20%) -> "
          f"{'HOLDS' if len(false_f) / len(fams) <= 0.20 else 'FAILS'}; {sorted(false_f)}")
    # O2 side
    P = {p["fam"]: [p["mask"]] + p["keep"] for p in json.load(open(f"{W}/panel.json"))}

    def place(bam, want):
        rr = collections.defaultdict(list)
        for rd in pysam.AlignmentFile(bam).fetch(until_eof=True):
            if rd.query_name in want and not rd.is_supplementary and not rd.is_unmapped:
                rr[rd.query_name].append((rd.is_secondary, rd.reference_name, rd.reference_start, rd.reference_end,
                                          rd.get_tag("AS") if rd.has_tag("AS") else 0, rd.get_tag("de") if rd.has_tag("de") else None))
        out = {}
        for n, rs in rr.items():
            srt = sorted(rs, key=lambda r: -r[4])
            prim = next((r for r in rs if not r[0]), srt[0])
            fam = lab[n]["family"]
            loc = ("cand", comp.get(prim[1], prim[1])) if prim[1].startswith("iso_") else \
                next((("copy", g) for c, s0, e0, g in P[fam] if prim[1] == c and prim[2] < e0 and s0 < prim[3]), ("other", prim[1]))
            tied = len(srt) > 1 and srt[1][4] > 0 and srt[1][4] >= 0.98 * srt[0][4]
            out[n] = ("tied" if tied else "placed", loc, prim[5])
        return out
    for L, k, n, fl, s in rows:
        if k != "absent_beyond_delta":
            continue
        want = on[L]
        for arm in ("R0", "C"):
            if not os.path.exists(f"{W}/{arm}.bam"):
                continue
            pl = place(f"{W}/{arm}.bam", want)
            c = collections.Counter(f"{st}:{loc[0]}:{loc[1]}" for st, loc, de in pl.values())
            des = [de for st, loc, de in pl.values() if de is not None and loc[0] != "cand"]
            print(f"  O2 {L} ({len(want)} reads) [{arm}]: {dict(c.most_common(4))}; median de on reference {statistics.median(des) if des else float('nan'):.4f}")


if __name__ == "__main__":
    main()
