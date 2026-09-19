#!/usr/bin/env python3
"""Secondary-dominated readthrough filter (§6n9).

WHAT IT IS. A readthrough record is FLAGGED when the reads carrying its FUSION junction are mostly
SECONDARY alignments, i.e. the record's defining junction exists at that locus mainly as a multimapping
echo of a molecule whose primary placement is elsewhere.

WHY IT IS EVIDENCE-BASED AND NOT A BLANKET RULE. Register 844 refuted the blunt "drop every readthrough"
node rule because it deletes `PKD1P6-NPIPP1`, a genuine NPIP member with 110 primary MAPQ-60 reads over a
canonical GT-AG junction (§6m1 addendum). This filter keeps that record and keeps `PDXDC2P-NPIPB14P`
(726 primary / 1 secondary) while flagging `PKD1P4-NPIPA8` (10 / 1,150) and `PKD1P3-NPIPA1` (214 / 628).
It separates the two populations §6n8 found instead of treating all readthroughs alike.

⚠ SCOPE AND A TRAP (§6n8). The MEDIAN readthrough in every class is purely primary-supported (secondary
fraction 0.000) and only 7 of 124 paralog-joining records are secondary-majority, so this flags a small
minority by design. Do not read it as "readthroughs are artifacts": supplementary support is 0.1%, so these
are contiguous alignments, not split ones. It flags where the EVIDENCE IS BORROWED, nothing more.

DEFAULT IS OFF. `--max-secondary-frac` defaults to 1.0, which flags nothing and is the explicit no-op.
0.50 is the "secondary-majority" setting §6n8 reports.

Usage:
  readthrough_secondary_filter.py --bam BAM --cuts cuts.tsv --gff GFF [--max-secondary-frac 0.50] [--out T]
"""
import argparse, collections, csv, re, subprocess, sys


def classify(counts, max_secondary_frac):
    """Pure rule. counts = {'primary':p,'secondary':s,'supplementary':u}.

    Returns (flagged, secondary_fraction). A record with NO fusion-junction support is never flagged —
    absence of evidence is not evidence of borrowing. At max_secondary_frac >= 1.0 nothing is ever flagged.
    """
    p = counts.get("primary", 0); s = counts.get("secondary", 0); u = counts.get("supplementary", 0)
    tot = p + s + u
    if tot == 0:
        return False, float("nan")
    frac = s / tot
    if max_secondary_frac >= 1.0:
        return False, frac
    return frac > max_secondary_frac, frac


def fusion_junction(exons, cut_coord):
    """The annotated intron spanning the cut point, as a 1-based (first intron base, first base of next exon)."""
    ex = sorted(exons)
    for i in range(len(ex) - 1):
        if ex[i][1] <= cut_coord <= ex[i + 1][0]:
            return (ex[i][1] + 1, ex[i + 1][0] + 1)
    return None


def junction_counts(bam, chrom, lo, hi, junction):
    cnt = collections.Counter()
    out = subprocess.run(["samtools", "view", bam, f"{chrom}:{lo}-{hi}"], capture_output=True, text=True).stdout
    for line in out.splitlines():
        f = line.split("\t"); flags = int(f[1]); pos = int(f[3])
        for n, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]):
            n = int(n)
            if op in "MD=X":
                pos += n
            elif op == "N":
                if (pos, pos + n) == junction:
                    cnt["supplementary" if flags & 0x800 else ("secondary" if flags & 0x100 else "primary")] += 1
                pos += n
    return cnt


def _self_test():
    assert classify({"primary": 10, "secondary": 1150}, 0.50)[0] is True
    assert classify({"primary": 726, "secondary": 1}, 0.50)[0] is False
    assert classify({"primary": 110}, 0.50)[0] is False, "PKD1P6-NPIPP1 must survive (register 844)"
    assert classify({"primary": 0, "secondary": 0}, 0.50)[0] is False, "no support is not evidence"
    assert classify({"primary": 1, "secondary": 999}, 1.0)[0] is False, "1.0 must be the no-op"
    # exactly at the threshold is NOT flagged (strict >)
    assert classify({"primary": 50, "secondary": 50}, 0.50)[0] is False
    assert fusion_junction([(10, 20), (40, 50)], 30) == (21, 41)
    assert fusion_junction([(10, 20), (40, 50)], 5) is None
    print("self-test OK")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam"); ap.add_argument("--cuts"); ap.add_argument("--gff")
    ap.add_argument("--max-secondary-frac", type=float, default=1.0)
    ap.add_argument("--out", default="-")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        _self_test(); return
    sys.path.insert(0, "/mnt/linuxdisk/home/juanfraitu/family_cert/dna")
    import dna_cert as dc
    nodes = dc.load_nodes()
    rows = []
    for r in csv.DictReader(open(a.cuts), delimiter="\t"):
        if r["status"] != "cut":
            continue
        n = nodes.get(r["node"])
        if not n:
            continue
        ex = sorted(n["exons"])
        j = fusion_junction(ex, int(r["cut_coord"]))
        if not j:
            continue
        c = junction_counts(a.bam, n["chrom"], ex[0][0], ex[-1][1], j)
        flagged, frac = classify(c, a.max_secondary_frac)
        rows.append((r["name"], n["chrom"], c.get("primary", 0), c.get("secondary", 0),
                     c.get("supplementary", 0), frac, flagged))
    fh = sys.stdout if a.out == "-" else open(a.out, "w")
    fh.write("name\tchrom\tprimary\tsecondary\tsupplementary\tsecondary_frac\tflagged\n")
    for nm, ch, p, s, u, fr, fl in sorted(rows, key=lambda x: -x[3]):
        fh.write(f"{nm}\t{ch}\t{p}\t{s}\t{u}\t{'' if fr != fr else f'{fr:.4f}'}\t{int(fl)}\n")
    n_fl = sum(1 for r in rows if r[6])
    print(f"[readthrough-filter] {len(rows)} records with fusion-junction support; "
          f"flagged {n_fl} at max_secondary_frac={a.max_secondary_frac}", file=sys.stderr)


if __name__ == "__main__":
    main()
