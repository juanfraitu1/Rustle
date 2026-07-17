#!/usr/bin/env python3
"""Score the coverage top-up experiment: for each of the 86 undetected Soto members, did giving it IDEAL
coverage (real-read replication -> TARGET reads) make the de-novo detector seed its locus?

Compares two runs of the SAME binary/command (control for version drift):
  control = gw_family_catalog on the ORIGINAL scoped BAM   -> soto_repro.copies.tsv
  topup   = gw_family_catalog on the REAL+TOPUP BAM        -> soto_topup.copies.tsv
Recovery = detected in topup but not control. A member still missing WITH ideal coverage is the true floor
(K=0 / secondary-sink), not a coverage miss.

Writes bench/soto/soto_topup_recovery.tsv + a summary.
Run: /home/juanfra/miniforge3/bin/python bench/soto/soto_topup_eval.py
"""
import csv
from collections import defaultdict

D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
CONTROL = f"{D}/soto_repro.copies.tsv"
TOPUP = f"{D}/soto_topup.copies.tsv"
TRUTH = f"{D}/topup/topup_truth.tsv"
PAF = "bench/soto/a119b_member_pairs.paf"
COV = "bench/soto/a119b_member_reads.tsv"
DET = "bench/soto/soto_member_detection.tsv"
TARGET = 40


def load_copies(path):
    """de-novo copies.tsv -> chrom -> [(start,end), ...]  (cols 3,4,5 = chrom,start,end)."""
    idx = defaultdict(list)
    try:
        r = csv.reader(open(path), delimiter="\t"); next(r, None)
        for row in r:
            if len(row) < 6:
                continue
            try:
                idx[row[3]].append((int(row[4]), int(row[5])))
            except ValueError:
                continue
    except FileNotFoundError:
        return None
    return idx


def hit(idx, chrom, s, e):
    return idx is not None and any(not (a > e or b < s) for (a, b) in idx.get(chrom, ()))


def best_sibling_id():
    """member coord (chrom:start1-end) -> best pairwise identity (matches/aln_len)."""
    bid = {}
    for line in open(PAF):
        f = line.split("\t")
        if len(f) < 11:
            continue
        alnlen = int(f[10])
        bid[f[0]] = max(bid.get(f[0], 0.0), int(f[9]) / alnlen if alnlen else 0.0)
    return bid


def main():
    control = load_copies(CONTROL)
    topup = load_copies(TOPUP)
    if control is None or topup is None:
        print(f"MISSING: control={CONTROL if control is None else 'ok'}  topup={TOPUP if topup is None else 'ok'}")
        return
    bid = best_sibling_id()

    rows = []
    n_rec = n_floor_topped = n_floor_wellcov = n_silent = n_control_hit = 0
    n_topped = 0
    for r in csv.DictReader(open(TRUTH), delimiter="\t"):
        fam, gene, chrom = r["family_id"], r["gene"], r["chrom"]
        s, e = int(r["start"]), int(r["end"])
        nr, add, status = int(r["primary_reads"]), int(r["n_added"]), r["status"]
        key = f"{chrom}:{s+1}-{e}"
        sib = bid.get(key)
        cdet = "Y" if hit(control, chrom, s, e) else "N"
        tdet = "Y" if hit(topup, chrom, s, e) else "N"
        if cdet == "Y":
            n_control_hit += 1  # was already detected in control (unexpected for the undetected set)
        if add > 0:
            n_topped += 1
        k0 = (sib is not None and sib >= 0.999)
        if status == "no-primary-template":
            cls = "silent (0 reads, no template)"; n_silent += 1
        elif tdet == "Y" and cdet == "N":
            cls = "RECOVERED by coverage"; n_rec += 1
        elif tdet == "Y" and cdet == "Y":
            cls = "detected in control (check)"
        elif add == 0:
            cls = "structural: well-covered (>=%d reads), coverage not the issue%s" % (TARGET, " [K=0 near-identical]" if k0 else "")
            n_floor_wellcov += 1
        else:
            cls = "structural: topped to ideal but still no locus%s" % (" [K=0 near-identical]" if k0 else "")
            n_floor_topped += 1
        rows.append([fam, gene, chrom, s, e, nr, add, f"{sib:.4f}" if sib is not None else "NA", cdet, tdet, cls])

    rows.sort(key=lambda x: (x[10], -x[6]))
    with open("bench/soto/soto_topup_recovery.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["family_id", "gene", "chrom", "start", "end", "primary_reads", "topup_added",
                    "best_sibling_id", "control_detected", "topup_detected", "classification"])
        w.writerows(rows)

    # ---- regression / consistency check over ALL 362 members (topup must not LOSE detected members) ----
    det = {}
    for line in open(DET):
        p = line.rstrip("\n").split("\t")
        if p[0] != "family_id" and len(p) >= 6:
            det[(p[0], p[1])] = p[5]
    reg_lost = []
    ctrl_dn = topup_dn = 0
    for line in open(COV):
        p = line.rstrip("\n").split("\t")
        if p[0] == "chrom" or len(p) < 6:
            continue
        chrom, s, e, gene, fam = p[0], int(p[1]), int(p[2]), p[3], p[4]
        ch = hit(control, chrom, s, e); th = hit(topup, chrom, s, e)
        ctrl_dn += ch; topup_dn += th
        if ch and not th:
            reg_lost.append(f"{fam}/{gene}")

    n_und = len(rows)
    floor = n_floor_topped + n_floor_wellcov + n_silent
    print(f"\n[consistency] de-novo member hits ALL 362: control={ctrl_dn} topup={topup_dn}; "
          f"members lost in topup (regression): {len(reg_lost)} {reg_lost[:8]}")
    print("\n=== SOTO coverage top-up experiment (real + simulated reads) ===")
    print(f"  undetected members (real data):        {n_und}")
    print(f"  control run detects (sanity, expect 0): {n_control_hit}")
    print(f"  RECOVERED by ideal coverage:           {n_rec}")
    print(f"  remaining floor (still missing):       {floor}")
    print(f"     - topped to ideal, still no locus:  {n_floor_topped}")
    print(f"     - well-covered (>= {TARGET}), never coverage: {n_floor_wellcov}")
    print(f"     - silent (0 reads, no template):    {n_silent}")
    k0floor = sum(1 for r in rows if "[K=0" in r[10])
    print(f"  of the floor, near-identical (>=99.9%) K=0: {k0floor}")
    print(f"  wrote bench/soto/soto_topup_recovery.tsv")


if __name__ == "__main__":
    main()
