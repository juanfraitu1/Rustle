#!/usr/bin/env python3
"""Register 1584 scorer (docs/PREREG_collapse_combined_2026-11-08.md): certificate-on-own-copy OR S/H pile, per missing copy.
usage: combined_rule.py DIR [S_THR H_THR]   (DIR holds psv.loci.tsv, read_origin.tsv, nt.bam, cp/cp.A.*, cp/cp.copies.tsv)
Without thresholds, S/H thresholds = 95th percentile of the autosomal single-origin null loci (own null)."""
import csv, sys, collections, math
import pysam
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__file__), "..", "crossind"))
from copyphasing_lookup import copies, ov

def pct95(xs):
    xs = sorted(x for x in xs if not math.isnan(x))
    return xs[min(len(xs) - 1, math.ceil(0.95 * len(xs)) - 1)] if xs else float("inf")

def main(d, s_thr=None, h_thr=None):
    origin = {l.split("\t")[0]: l.split("\t")[1] for l in open(f"{d}/read_origin.tsv") if not l.startswith("read\t")}
    bam = pysam.AlignmentFile(f"{d}/nt.bam")
    by_key, rows = copies(f"{d}/cp/cp.copies.tsv")
    exons = {}
    for r in csv.DictReader(open(f"{d}/cp/cp.copies.tsv"), delimiter="\t"):
        exons[(r["chrom"], int(r["start"]), int(r["end"]), r["family_id"])] = [tuple(map(int, x.split("-"))) for x in r["exons"].split(",") if x]
    def res(fam, home):
        for k in (home, f"{fam}:{home}", (fam, home)):
            if k in by_key:
                return by_key[k]
        sys.exit(f"unresolved {home}")
    homes, cand_home = set(), collections.defaultdict(list)
    for r in csv.DictReader(open(f"{d}/cp/cp.A.candidates.tsv"), delimiter="\t"):
        if r["flagged"] == "1":
            homes.add(res(r["family"], r["nearest_locus"])); cand_home[res(r["family"], r["nearest_locus"])].append(r["candidate"])
    cl2cand = {r["cluster"]: r["candidate"] for r in csv.DictReader(open(f"{d}/cp/cp.A.clusters.tsv"), delimiter="\t")}
    cand_reads = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{d}/cp/cp.A.reads.tsv"), delimiter="\t"):
        if r["cluster"] in cl2cand:
            cand_reads[cl2cand[r["cluster"]]].append(r["read"])
    def purity(copy, gene):
        """labelled: best share of the gene's own reads among certified candidates homed on this copy"""
        best = 0.0
        for c in cand_home.get(copy, []):
            rs = cand_reads[c]
            if rs:
                best = max(best, sum(origin.get(q) == gene for q in rs) / len(rs))
        return best
    for r in csv.DictReader(open(f"{d}/cp/cp.A.putative.tsv"), delimiter="\t"):
        homes.add(res(r["family"], r["copy"]))

    def own_copy(gene, w):
        """catalog copy (exon overlap) holding the most of this gene's labelled primary reads inside window w"""
        c = collections.Counter()
        cands = [x for x in rows if ov(x, w)]
        for a in bam.fetch(w[0], w[1], w[2]):
            if a.is_secondary or a.is_supplementary or a.is_unmapped or origin.get(a.query_name) != gene:
                continue
            blocks = a.get_blocks()
            for x in cands:
                if any(s < be and bs < e for s, e in exons[x] for bs, be in blocks):
                    c[x] += 1
        return c.most_common(1)[0][0] if c else None

    loci = [r for r in csv.DictReader(open(f"{d}/psv.loci.tsv"), delimiter="\t") if r["kind"] in ("absorbing", "null")]
    seen_win, units = set(), []
    for r in loci:
        w = (r["chrom"], int(r["start"]), int(r["end"]))
        auto = r["chrom"] not in ("chrX", "chrY")
        if r["kind"] == "null":
            cnt = collections.Counter(origin.get(a.query_name) for a in bam.fetch(*w)
                                      if not (a.is_secondary or a.is_supplementary or a.is_unmapped) and origin.get(a.query_name))
            n = sum(cnt.values())
            if not n or cnt.most_common(1)[0][0] != r["gene"] or any(v >= 0.10 * n for g, v in cnt.items() if g != r["gene"]):
                continue  # not single-origin
        elif w in seen_win:
            continue  # identical absorbing window counted once (NPIPB15 x2)
        seen_win.add(w)
        units.append((r, w, auto))
    S = lambda r: float(r["S"]) if r["S"] not in ("", "nan", "NA") else float("nan")
    H = lambda r: float(r["H"]) if r["H"] not in ("", "nan", "NA") else float("nan")
    nulls = [r for r, w, a in units if a and r["kind"] == "null"]
    if s_thr is None:
        s_thr, h_thr = pct95(map(S, nulls)), pct95(map(H, nulls))
    print(f"S/H thresholds S {s_thr:.3f} H {h_thr:.3f}  (own null n={len(nulls)})" if len(sys.argv) < 3 else f"S/H thresholds frozen S {s_thr} H {h_thr}")
    T = collections.defaultdict(lambda: [0, 0, 0, 0, 0])  # testable, cert, sh, combined, labelled-correct cert
    print("name\tkind\treason\tstratum\tlocus\town_copy\tcert\tpurity\tS\tH\tsh\tcombined")
    for r, w, auto in units:
        oc = own_copy(r["gene"], w)
        cert = oc is not None and oc in homes
        sh = S(r) > s_thr or H(r) > h_thr
        comb = cert or sh
        pu = purity(oc, r["gene"]) if cert else 0.0
        st = "auto" if auto else "XY"
        keys = [f"{st} {r['kind']}"] + ([f"{st} missing {r['reason']}"] if r["kind"] == "absorbing" else [])
        for k in keys:
            t = T[k]; t[0] += 1; t[1] += cert; t[2] += sh; t[3] += comb; t[4] += pu >= 0.5
        print(f"{r['name']}\t{r['kind']}\t{r['reason']}\t{st}\t{w[0]}:{w[1]}-{w[2]}\t{oc[3] if oc else '-'}\t{int(cert)}\t{pu:.2f}\t{S(r):.2f}\t{H(r):.3f}\t{int(sh)}\t{int(comb)}")
    for k in sorted(T):
        t = T[k]
        print(f"{k}: n {t[0]}  cert {t[1]}  S/H {t[2]}  combined {t[3]}" + (f" = {t[3]/t[0]:.2f}" if t[0] else "") + f"  labelled-correct cert (purity >= 0.5) {t[4]}")
    m, n = T["auto absorbing"], T["auto null"]
    if m[0] < 5 or n[0] < 5:
        print("VERDICT UNDERPOWERED"); return
    b1, b2 = m[3] / m[0] >= 0.5, n[3] / n[0] <= 0.25
    print(f"bar1 {'PASS' if b1 else 'FAIL'} ({m[3]}/{m[0]})  bar2 {'PASS' if b2 else 'FAIL'} ({n[3]}/{n[0]})  VERDICT {'PASS' if b1 and b2 else 'FAIL'}")

if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], *(map(float, a[1:3]) if len(a) >= 3 else ()))
