#!/usr/bin/env python3
"""Genome-wide StringTie-EXACT floor aggregator.

Proves the baseline FLOOR genome-wide: rustle eonly-guided (`-e -G st.gtf`)
reproduces StringTie's transcript set exactly. Two checks per chromosome:

  1. PARITY (ee vs StringTie's own GTF, gc_ee_par_*): transcript-level Sn/Pr
     should be 100/100 and the transcript counts equal. Any chrom with
     ee_tx != st_tx or Sn/Pr < 100 is flagged.
  2. FLOOR NCBI MATCH (ee vs NCBI, gc_ee_*) vs StringTie's NCBI match
     (gc_st_*): the floor should match the SAME NCBI transcripts StringTie
     does (strict =/c) -> confirms the floor is StringTie, not a different set.

Run after bench/gw_eonly_parity.sh. Reuses st_$C.gtf, gc_st_$C, gc_ee_par_$C,
gc_ee_$C from /tmp/gw.
"""
import glob
import os
import re
from collections import defaultdict

OUT = "/tmp/gw"
STRICT = {"=", "c"}

_TX_RE = re.compile(r"Transcript level:\s*([\d.]+)\s*\|\s*([\d.]+)")


def txcount_gtf(path):
    if not os.path.exists(path):
        return 0
    n = 0
    with open(path) as fh:
        for line in fh:
            if "\ttranscript\t" in line:
                n += 1
    return n


def parse_stats_tx(prefix, chrom):
    """Return (Sn, Pr) transcript-level from a gffcompare summary, or (None, None).

    gffcompare writes the summary stats to the extension-less `-o` path itself
    (e.g. gc_ee_par_<chrom>), not to a `.stats` file.
    """
    cands = [
        os.path.join(OUT, f"{prefix}_{chrom}"),          # extension-less (gffcompare default)
        os.path.join(OUT, f"{prefix}_{chrom}.stats"),    # fallback
    ]
    for p in cands:
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            for line in fh:
                m = _TX_RE.search(line)
                if m:
                    return float(m.group(1)), float(m.group(2))
    return None, None


def tmap_refset(prefix, chrom, codes):
    paths = glob.glob(os.path.join(OUT, f"{prefix}_{chrom}.*.tmap"))
    if not paths:
        return None
    s = set()
    with open(paths[0]) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            ref_id, cls = c[1], c[2]
            if ref_id not in ("-", "") and cls in codes:
                s.add((chrom, ref_id))
    return s


def done_chroms():
    return sorted(os.path.basename(f)[len("ee_"):-len(".done")]
                  for f in glob.glob(os.path.join(OUT, "ee_*.done")))


def main():
    chroms = done_chroms()
    if not chroms:
        print("no ee_*.done chroms yet — run bench/gw_eonly_parity.sh first")
        return

    tot_ee = tot_st = 0
    tot_ee_ncbi = tot_st_ncbi = 0
    flagged = []
    per_chrom = []

    for c in chroms:
        ee_tx = txcount_gtf(os.path.join(OUT, f"ee_{c}.gtf"))
        st_tx = txcount_gtf(os.path.join(OUT, f"st_{c}.gtf"))
        sn, pr = parse_stats_tx("gc_ee_par", c)
        ee_ncbi = tmap_refset("gc_ee", c, STRICT)
        st_ncbi = tmap_refset("gc_st", c, STRICT)
        ee_ncbi_n = len(ee_ncbi) if ee_ncbi is not None else -1
        st_ncbi_n = len(st_ncbi) if st_ncbi is not None else -1

        tot_ee += ee_tx
        tot_st += st_tx
        if ee_ncbi_n >= 0:
            tot_ee_ncbi += ee_ncbi_n
        if st_ncbi_n >= 0:
            tot_st_ncbi += st_ncbi_n

        ok = (ee_tx == st_tx) and (sn == 100.0) and (pr == 100.0)
        if not ok:
            flagged.append((c, ee_tx, st_tx, sn, pr))
        per_chrom.append((c, ee_tx, st_tx, sn, pr, ee_ncbi_n, st_ncbi_n))

    print("=" * 84)
    print("GENOME-WIDE STRINGTIE-EXACT FLOOR  (ee = rustle `-e -G st.gtf` vs StringTie)")
    print("=" * 84)
    print(f"{'chrom':<18}{'ee_tx':>8}{'st_tx':>8}{'Sn':>8}{'Pr':>8}{'ee_NCBI':>9}{'st_NCBI':>9}{'':>4}")
    for c, ee, st, sn, pr, een, stn in per_chrom:
        sns = "n/a" if sn is None else f"{sn:.1f}"
        prs = "n/a" if pr is None else f"{pr:.1f}"
        mark = "" if (ee == st and sn == 100.0 and pr == 100.0) else "  <-- DEVIATES"
        print(f"{c:<18}{ee:>8}{st:>8}{sns:>8}{prs:>8}{een:>9}{stn:>9}{mark}")
    print("-" * 84)
    print("GENOME-WIDE TOTALS")
    print(f"  rustle eonly-floor transcripts : {tot_ee}")
    print(f"  StringTie transcripts          : {tot_st}")
    print(f"  difference (floor - StringTie)  : {tot_ee - tot_st}  (0 == exact floor)")
    print(f"  floor NCBI-matched (strict =/c): {tot_ee_ncbi}")
    print(f"  StringTie NCBI-matched (=/c)   : {tot_st_ncbi}")
    print(f"  NCBI difference                : {tot_ee_ncbi - tot_st_ncbi}  (0 == same NCBI set)")
    print("-" * 84)
    if flagged:
        print(f"DEVIATING CHROMS ({len(flagged)}): the floor is NOT byte-exact on these")
        for c, ee, st, sn, pr in flagged:
            print(f"  {c}: ee_tx={ee} st_tx={st} Sn={sn} Pr={pr}")
    else:
        print("FLOOR IS STRINGTIE-EXACT ON EVERY CHROMOSOME (ee_tx == st_tx, Sn=Pr=100).")
    print("=" * 84)


if __name__ == "__main__":
    main()
