#!/usr/bin/env python3
"""r1509 post hoc (flagged): re-read haplotype_truth's hits at a STRICTER slack. haplotype_truth calls a read's truth
the copies overlapping every haplotype hit within NM <= best + 2, which on near-identical copies of one individual makes
almost every read "ambiguous" by construction (a 1-2 edit difference IS the PSV evidence O2 uses). Here the slack g is
0, 1 or 2: truth_g = copies overlapping hits with NM <= best + g; a row is decidable at g when truth_g is one locus.
Reports, per g and status: decidable, O2 right / wrong, aligner right. Usage: strict_truth.py ARM_PREFIX TRUTH_DIR
CATALOG_PREFIX [MOLECULES.tsv: adds the aligner's primary verdict]"""
import collections, csv, os, subprocess, sys

def main(prefix, truth_dir, catalog, molecules=None):
    prim = {}  # read -> (chrom, pos) of the genome-wide primary (extended as_table), for the aligner's verdict
    if molecules:
        with open(molecules) as fh:
            for line in fh:
                if line[0] == "#": continue
                f = line.rstrip("\n").split("\t")
                if len(f) >= 16 and f[4] != "NA": prim[f[0]] = (f[4], int(f[15]))
    span = {}
    with open(catalog + ".fam.copies.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            span[(r["family_id"], r["copy_idx"])] = (r["chrom"], int(r["start"]), int(r["end"]))
    hits = collections.defaultdict(list)
    with open(f"{truth_dir}/hits.tsv") as fh:
        next(fh)
        for line in fh:
            q, hap, ch, s, e, nm, al = line.rstrip("\n").split("\t")
            hits[q].append((hap, ch, int(s), int(e), int(nm), int(al)))
    loci = {}
    for q, hs in hits.items():
        if q.startswith("C|"):
            _, f, i = q.split("|")
            best = {}
            for h in hs: best[h[0]] = min(best.get(h[0], 10 ** 9), h[4])
            # a copy's locus on a haplotype = hits within 2 of its best there, given a NEAR-PERFECT best (NM <= max(2, 0.003 x aligned)); best + 2 alone let a
            # copy with no good hit on one haplotype (MCL13 copy 2, pat best NM 5) claim its parent gene's locus (r1510 post hoc)
            loci[(f, i)] = {(h[0], h[1], h[2], h[3]) for h in hs if best[h[0]] <= max(2, 0.003 * h[5]) and h[4] <= best[h[0]] + 2}
    def overlaps(a, b): return a[0] == b[0] and a[1] < b[2] and b[1] < a[2]
    def one_locus(ks): return bool(ks) and all(overlaps(span[a], span[b]) for i, a in enumerate(ks) for b in ks[i + 1:])
    rows = []
    with open(prefix + ".assignments.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["contested"] == "1" and "R|" + r["read_name"] in hits: rows.append(r)
    # primary span from the read's region-local row is not stored; use the aligner-right flag only where haplotype_truth gave it
    dump = open(os.environ["STRICT_DUMP"], "w") if os.environ.get("STRICT_DUMP") else None  # per-row verdicts at g = 0
    print("g\tstatus\trows_with_hits\tno_near_perfect\tdecidable\toutside\to2_right\to2_wrong\taligner_right\taligner_wrong")
    for g in (0, 1, 2):
        c = collections.defaultdict(collections.Counter)
        for r in rows:
            hs = hits["R|" + r["read_name"]]; st = "assigned" if r["status"] == "assigned" else "abstained"
            c[st]["rows"] += 1
            best = min(hs, key=lambda h: (h[4], -h[5]))
            if best[4] > max(2, 0.003 * best[5]): c[st]["no_near"] += 1; continue
            near = [h for h in hs if h[4] <= best[4] + g]
            T = sorted({k for k, ls in loci.items() if any(lh == h[0] and lc == h[1] and ls0 < h[3] and h[2] < le for lh, lc, ls0, le in ls for h in near)})
            if not T: c[st]["outside"] += 1; continue
            if not one_locus(T): c[st]["ambiguous"] += 1; continue
            c[st]["decidable"] += 1
            v = "NA"
            if st == "assigned": v = "right" if (r["family_id"], r["catalog_copy_idx"]) in T else "wrong"; c[st][v] += 1
            p = prim.get(r["read_name"]); av = "NA"
            if p: av = "al_right" if any(span[k][0] == p[0] and span[k][1] <= p[1] <= span[k][2] for k in T) else "al_wrong"; c[st][av] += 1
            if g == 0 and dump: dump.write(f"{r['read_name']}\t{r['family_id']}\t{st}\t{v}\t{av}\t{r['catalog_copy_idx']}\t{';'.join(f'{k[0]}:{k[1]}' for k in T)}\n")
        for st in ("assigned", "abstained"):
            x = c[st]; print("\t".join(map(str, [g, st, x["rows"], x["no_near"], x["decidable"], x["outside"], x["right"], x["wrong"], x["al_right"], x["al_wrong"]])))

if __name__ == "__main__": main(*sys.argv[1:5])
