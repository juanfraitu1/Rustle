#!/usr/bin/env python3
"""Per-family/member detection SENSITIVITY *and* PRECISION on the real Soto benchmark (80_fams.chr.bed).

Sensitivity(F) = detected members / planted members (recall of F's members).
Precision(F)   = detected loci in F's genomic region that hit a TRUE member / all detected loci in F's region.
  (A detected locus in F's region that overlaps no annotated member is either a false call or a candidate
  UNANNOTATED paralog of F — Soto's 80-family subset is not exhaustive, so most are the latter.)

Legs (precedence for the per-member `recovered_by` label): RNA-split > protein-tail > projection > expressed-collapse.
Outputs: bench/soto/soto_member_detection.tsv, bench/soto/soto_family_detection.tsv + a summary.
Run: /home/juanfra/miniforge3/bin/python bench/soto/soto_detection_eval.py
"""
import csv
from collections import defaultdict

BED = "bench/soto/80_fams.chr.bed"
D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
LEGS = [  # (label, path, chrom_col, start_col, end_col, projection_loci_col_or_None)
    ("RNA-split",     f"{D}/soto_off.copies.tsv",            3, 4, 5, None),
    ("protein-tail",  f"{D}/soto_gw_prot.copies.tsv",         3, 4, 5, None),
    ("projection",    f"{D}/soto_pall.allproj.tsv",           1, 2, 3, None),
    ("expr-collapse", f"{D}/soto_ce.expressed_collapsed.tsv", 1, 2, 3, 7),
]


def load_bed(path):
    out = []
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) < 4 or not p[0].startswith("chr"):
            continue
        gene, fam = (p[3].split("|") + ["?"])[:2]
        out.append((p[0], int(p[1]), int(p[2]), gene, fam))
    return out


def load_loci(path, cc, sc, ec, proj_col):
    loci = []
    try:
        r = csv.reader(open(path), delimiter="\t"); next(r, None)
        for row in r:
            if len(row) <= max(cc, sc, ec):
                continue
            loci.append((row[cc], int(row[sc]), int(row[ec])))
            if proj_col is not None and len(row) > proj_col:
                for seg in row[proj_col].split(";"):
                    if seg:
                        ch, rng = seg.split("@")[0].split(":"); s, e = rng.split("-")
                        loci.append((ch, int(s), int(e)))
    except FileNotFoundError:
        pass
    return loci


def by_chrom(ivs):
    d = defaultdict(list)
    for t in ivs:
        d[t[0]].append(t[1:])
    return d


def hits(chrom, s, e, idx):
    return any(not (a > e or b < s) for (a, b, *_) in idx.get(chrom, ()))


def main():
    members = load_bed(BED)
    leg_idx = [(lab, by_chrom(load_loci(p, cc, sc, ec, pc))) for (lab, p, cc, sc, ec, pc) in LEGS]
    all_loci = [iv for (lab, p, cc, sc, ec, pc) in LEGS for iv in load_loci(p, cc, sc, ec, pc)]

    # per-family: members + genomic region (per chrom min-start/max-end) + a member index for correctness.
    fam_members = defaultdict(list)
    fam_region = defaultdict(lambda: defaultdict(lambda: [10**12, 0]))  # fam -> chrom -> [lo, hi]
    member_index = by_chrom([(c, s, e, fam) for (c, s, e, g, fam) in members])
    for (c, s, e, g, fam) in members:
        fam_members[fam].append((g, c, s, e))
        r = fam_region[fam][c]; r[0] = min(r[0], s); r[1] = max(r[1], e)

    # ---- per-member sensitivity ----
    mrows, det_of = [], {}
    for (c, s, e, g, fam) in members:
        lever = next((lab for lab, idx in leg_idx if hits(c, s, e, idx)), "")
        det_of[(fam, g)] = bool(lever)
        mrows.append([fam, g, c, s, e, "Y" if lever else "N", lever])

    # ---- per-family precision: attribute each detected locus ----
    # correct = overlaps a real member of F; extra = in F's region but hits no member (candidate unannotated).
    fam_correct, fam_extra = defaultdict(int), defaultdict(int)
    for (c, s, e) in all_loci:
        hit_fams = {fam for (a, b, fam) in member_index.get(c, ()) if not (a > e or b < s)}
        if hit_fams:
            for fam in hit_fams:
                fam_correct[fam] += 1
        else:
            # off-annotation: attribute to the family whose region contains this locus (if any)
            mid = (s + e) // 2
            for fam, chs in fam_region.items():
                if c in chs and chs[c][0] <= mid <= chs[c][1]:
                    fam_extra[fam] += 1
                    break

    # ---- family table ----
    frows = []
    for fam in sorted(fam_members, key=lambda f: int(f.replace("ID_", ""))):
        mem = fam_members[fam]
        nd = sum(1 for (g, _, _, _) in mem if det_of[(fam, g)])
        corr, extra = fam_correct[fam], fam_extra[fam]
        prec = f"{100*corr/(corr+extra):.0f}%" if (corr + extra) else "NA"
        frows.append([fam, len(mem), nd, f"{nd}/{len(mem)}", f"{100*nd/len(mem):.0f}%",
                      corr + extra, prec, ";".join(f"{g}:{'Y' if det_of[(fam, g)] else 'N'}" for (g, _, _, _) in mem)])

    with open("bench/soto/soto_member_detection.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(["family_id", "gene", "chrom", "start", "end", "detected", "recovered_by"]); w.writerows(mrows)
    with open("bench/soto/soto_family_detection.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["family_id", "n_members", "n_detected", "sensitivity", "sens_pct", "n_pred_loci", "precision", "members"]); w.writerows(frows)

    n_mem = len(members); nd = sum(1 for r in mrows if r[5] == "Y")
    n_fam = len(fam_members); fam_full = sum(1 for fr in frows if fr[2] == fr[1])
    member_idx = by_chrom([(c, s, e) for (c, s, e, _, _) in members])
    on = sum(1 for (c, s, e) in all_loci if hits(c, s, e, member_idx))
    print("\n=== SOTO detection — SENSITIVITY & PRECISION (real 80_fams benchmark) ===")
    print(f"  member sensitivity:  {nd}/{n_mem} = {100*nd/n_mem:.1f}%")
    print(f"  overall precision:   {on}/{len(all_loci)} = {100*on/len(all_loci):.0f}%  ({len(all_loci)-on} off-annotation = candidate unannotated paralogs)")
    print(f"  families fully recovered: {fam_full}/{n_fam} = {100*fam_full/n_fam:.0f}%")
    print(f"  wrote soto_family_detection.tsv (now with sensitivity + precision cols) + soto_member_detection.tsv")


if __name__ == "__main__":
    main()
