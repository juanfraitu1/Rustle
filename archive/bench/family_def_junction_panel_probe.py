#!/usr/bin/env python3
"""
family_def_junction_panel_probe.py — diagnose, per airtight-panel pair AND per known
retrocopy counterexample, exactly how the intron-architecture criteria score them.
Tells us whether Carch is cutting REAL families (a cost) or only domain-bridges/retrocopies.

For each named pair we print: cDNA homology (id, maxcov, is it a baseline edge?), the two
genes' intron architectures (n_intron, intron_lengths head), and the verdict of each
criterion (Cspliced / Ccount / Carch). This is the per-case panel result the task wants.

Run: /home/juanfra/miniforge3/bin/python bench/family_def_junction_panel_probe.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_read_filters import dna_homology
from family_def_junction_genomewide import (arch_concordant, count_compatible,
                                            ID_MIN, COV_MIN)

INDEX = "/home/juanfra/winloci_scratch/gene_intron_index.json"

# REAL families (should be KEPT)
REAL_PAIRS = [
    ("RABL2", "RABL2A", "RABL2B"),
    ("APOBEC3", "APOBEC3D", "APOBEC3F"),
    ("RFPL", "RFPL1", "RFPL2"),
    ("RFPL", "RFPL2", "RFPL3"),
    ("RFPL", "RFPL1", "RFPL3"),
]
# COUNTER (domain-sharer co-location) — should be CUT
COUNTER_PAIRS = [
    ("CASP8~FLACC1", "CASP8", "FLACC1"),
    ("ASDURF~ASNSD1", "ASDURF", "ASNSD1"),
    ("CREB1~METTL21A", "CREB1", "METTL21A"),
    ("GPR39~LYPD1", "GPR39", "LYPD1"),
]
# RETROCOPY counterexamples — intronless processed copy, id~1.0 — should be CUT
RETRO_PAIRS = [
    ("EEF1A1-retro", "EEF1A1", None),   # we resolve the retrocopy partner from homology
    ("CNN2-retro", "CNN2", None),
]


def arch_str(rec):
    if rec is None:
        return "NO-INDEX"
    il = rec["intron_lengths"]
    return f"n_intron={rec['n_intron']:2d} strand={rec['strand']} il_head={il[:4]}"


def verdict(ra, rb):
    if ra is None or rb is None:
        return dict(Cspliced=False, Ccount=False, Carch=False, note="missing index")
    spliced = ra["n_intron"] >= 1 and rb["n_intron"] >= 1
    return dict(Cspliced=spliced,
                Ccount=count_compatible(ra, rb),
                Carch=arch_concordant(ra, rb))


def homol(H, a, b):
    k = (a, b) if a < b else (b, a)
    r = H.get(k)
    if r is None:
        return None
    return r


def is_baseline_edge(H, a, b):
    r = homol(H, a, b)
    if r is None:
        return False
    return r["id"] >= ID_MIN and max(r["cov_a"], r["cov_b"]) >= COV_MIN


def main():
    H, dna = dna_homology()
    idx = json.load(open(INDEX))

    def show(group, a, b):
        ra, rb = idx.get(a), idx.get(b)
        r = homol(H, a, b)
        hstr = (f"id={r['id']:.3f} maxcov={max(r['cov_a'],r['cov_b']):.2f} "
                f"edge={is_baseline_edge(H,a,b)}") if r else "NO-HOMOLOGY"
        v = verdict(ra, rb)
        print(f"  [{group:16}] {a:14} ~ {b:14}")
        print(f"       homology : {hstr}")
        print(f"       arch {a:12}: {arch_str(ra)}")
        print(f"       arch {b:12}: {arch_str(rb)}")
        print(f"       verdict  : Cspliced={v['Cspliced']}  Ccount={v['Ccount']}  Carch={v['Carch']}")

    print("=== REAL FAMILIES (criteria should KEEP) ===")
    for grp, a, b in REAL_PAIRS:
        show(grp, a, b)

    print("\n=== DOMAIN-SHARER COUNTER (criteria should CUT) ===")
    for grp, a, b in COUNTER_PAIRS:
        show(grp, a, b)

    print("\n=== RETROCOPY COUNTER (intronless processed copy; should CUT) ===")
    # resolve each retro's best homologous partner from the homology table
    for grp, a, _ in RETRO_PAIRS:
        if a not in idx:
            print(f"  [{grp}] {a}: NOT in intron index"); continue
        # find best homologous partner (highest id, maxcov>=0.3)
        partners = []
        for (x, y), r in H.items():
            if r["id"] < 0.90 or max(r["cov_a"], r["cov_b"]) < 0.30:
                continue
            if x == a:
                partners.append((r["id"], y))
            elif y == a:
                partners.append((r["id"], x))
        partners.sort(reverse=True)
        if not partners:
            print(f"  [{grp}] {a}: no high-homology partner in rep_ava (cDNA may be absent)")
            print(f"       arch {a:12}: {arch_str(idx.get(a))}")
            continue
        # show the intronless gene's own architecture + its top partners
        print(f"  [{grp}] {a}: n_intron={idx[a]['n_intron']} ; top homologous partners:")
        for ident, p in partners[:4]:
            rp = idx.get(p)
            v = verdict(idx.get(a), rp)
            print(f"       ~ {p:14} id={ident:.3f}  partner_n_intron="
                  f"{rp['n_intron'] if rp else 'NA'}  Cspliced={v['Cspliced']} Carch={v['Carch']}")

    # SUMMARY scorecard
    print("\n=== SCORECARD (separation) ===")
    def score(pairs, want_keep):
        rows = []
        for grp, a, b in pairs:
            ra, rb = idx.get(a), idx.get(b)
            v = verdict(ra, rb)
            rows.append((grp, a, b, v))
        return rows
    real_rows = score(REAL_PAIRS, True)
    counter_rows = score(COUNTER_PAIRS, False)
    for mode in ("Cspliced", "Ccount", "Carch"):
        rk = sum(1 for *_x, v in real_rows if v[mode])           # real kept (good)
        ck = sum(1 for *_x, v in counter_rows if v[mode])        # counter kept (bad)
        print(f"  {mode:9}: REAL kept {rk}/{len(real_rows)}   COUNTER(domain) kept {ck}/{len(counter_rows)}")


if __name__ == "__main__":
    main()
