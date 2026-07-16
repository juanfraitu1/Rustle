#!/usr/bin/env python3
"""Per-member detection sensitivity on the REAL Soto benchmark (bench/soto/80_fams.chr.bed).

For every Soto family (ID_xxx) and its member genes, mark whether the member's genomic interval is recovered
by any pipeline leg, and by WHICH leg (precedence: RNA-split > protein-tail > projection > expressed-collapsed).
Sensitivity = detected / planted, per family and overall. Precision = detected loci that overlap a real Soto
member / all detected loci (a member-level correctness check; off-annotation loci are candidate unannotated
paralogs, reported separately).

Outputs: bench/soto/soto_member_detection.tsv (one row per member), bench/soto/soto_family_detection.tsv
(one row per family) + a summary block. Run: /home/juanfra/miniforge3/bin/python bench/soto/soto_detection_eval.py
"""
import csv
from collections import defaultdict

BED = "bench/soto/80_fams.chr.bed"
D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
LEGS = [  # (label, path, chrom_col, start_col, end_col, extra_projection_col_or_None)
    ("RNA-split",    f"{D}/soto_off.copies.tsv",           3, 4, 5, None),
    ("protein-tail", f"{D}/soto_gw_prot.copies.tsv",        3, 4, 5, None),
    ("projection",   f"{D}/soto_pall.allproj.tsv",          1, 2, 3, None),
    ("expr-collapse",f"{D}/soto_ce.expressed_collapsed.tsv", 1, 2, 3, 7),   # +projection_loci col 7
]


def load_bed(path):
    members = []  # (chrom, start, end, gene, family)
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or not p[0].startswith("chr"):
                continue
            gene, fam = (p[3].split("|") + ["?"])[:2]
            members.append((p[0], int(p[1]), int(p[2]), gene, fam))
    return members


def load_loci(path, cc, sc, ec, proj_col):
    loci = []  # (chrom, start, end)
    try:
        with open(path) as f:
            r = csv.reader(f, delimiter="\t")
            next(r, None)
            for row in r:
                if len(row) <= max(cc, sc, ec):
                    continue
                loci.append((row[cc], int(row[sc]), int(row[ec])))
                if proj_col is not None and len(row) > proj_col:
                    for seg in row[proj_col].split(";"):
                        if not seg:
                            continue
                        coord = seg.split("@")[0]
                        ch, rng = coord.split(":")
                        s, e = rng.split("-")
                        loci.append((ch, int(s), int(e)))
    except FileNotFoundError:
        pass
    return loci


def by_chrom(loci):
    d = defaultdict(list)
    for ch, s, e in loci:
        d[ch].append((s, e))
    return d


def overlaps(chrom, s, e, idx):
    return any(not (ls > e or le < s) for (ls, le) in idx.get(chrom, ()))


def main():
    members = load_bed(BED)
    leg_idx = [(label, by_chrom(load_loci(path, cc, sc, ec, pc))) for (label, path, cc, sc, ec, pc) in LEGS]

    # per-member detection (first leg in precedence that covers it wins the label)
    rows = []
    fam_members = defaultdict(list)
    for (ch, s, e, gene, fam) in members:
        lever = ""
        for label, idx in leg_idx:
            if overlaps(ch, s, e, idx):
                lever = label
                break
        detected = "Y" if lever else "N"
        rows.append([fam, gene, ch, s, e, detected, lever])
        fam_members[fam].append((gene, detected, lever))

    # per-family summary
    fam_rows = []
    for fam in sorted(fam_members, key=lambda f: int(f.replace("ID_", ""))):
        mem = fam_members[fam]
        nd = sum(1 for (_, d, _) in mem if d == "Y")
        fam_rows.append([fam, len(mem), nd, f"{nd}/{len(mem)}", f"{100*nd/len(mem):.0f}%",
                         ";".join(f"{g}:{d}" for (g, d, _) in mem)])

    with open("bench/soto/soto_member_detection.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(["family_id", "gene", "chrom", "start", "end", "detected", "recovered_by"]); w.writerows(rows)
    with open("bench/soto/soto_family_detection.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(["family_id", "n_members", "n_detected", "sensitivity", "sens_pct", "members"]); w.writerows(fam_rows)

    # precision: detected loci across legs that overlap a real Soto member / all detected loci
    member_idx = by_chrom([(ch, s, e) for (ch, s, e, _, _) in members])
    all_loci = []
    for (label, path, cc, sc, ec, pc) in LEGS:
        all_loci += load_loci(path, cc, sc, ec, pc)
    on = sum(1 for (ch, s, e) in all_loci if overlaps(ch, s, e, member_idx))

    n_mem = len(members)
    nd = sum(1 for r in rows if r[5] == "Y")
    n_fam = len(fam_members)
    fam_det = sum(1 for fr in fam_rows if fr[2] >= 1)
    fam_full = sum(1 for fr in fam_rows if fr[2] == fr[1])
    print("\n=== SOTO per-member detection (real benchmark, all levers) ===")
    print(f"  members detected (sensitivity): {nd}/{n_mem} = {100*nd/n_mem:.1f}%")
    print(f"  families with >=1 member detected: {fam_det}/{n_fam} = {100*fam_det/n_fam:.0f}%")
    print(f"  families FULLY recovered (all members): {fam_full}/{n_fam} = {100*fam_full/n_fam:.0f}%")
    print(f"  precision (detected loci overlapping a real Soto member): {on}/{len(all_loci)} = {100*on/max(1,len(all_loci)):.0f}%")
    print(f"    (the {len(all_loci)-on} off-annotation loci are candidate unannotated paralogs, id>=0.98/read-supported)")
    print(f"  recovered_by breakdown: ", dict((l, sum(1 for r in rows if r[6] == l)) for l in ["RNA-split","protein-tail","projection","expr-collapse"]))
    print("  wrote bench/soto/soto_member_detection.tsv (362 members) + soto_family_detection.tsv (per family)")


if __name__ == "__main__":
    main()
