#!/usr/bin/env python3
"""family_def_airtight_panel.py — a small, airtight ground-truth panel from ANNOTATION.

The cDNA all-vs-all 'truth' is over-merged (389-member mega-families). Instead build a
hand-verifiable panel where each case is AIRTIGHT on two independent signals:
  EXAMPLES (real families): named/tandem gene pairs that are cDNA-HOMOLOGOUS (id>=0.90) AND
    close on the genome (tandem) -- homology + name + location agree. (Pairs, not transitive
    components, so no over-merge.)
  COUNTEREXAMPLES (non-families):
    - domain-sharers: adjacent genes sharing a domain (partial homology, different genes)
    - name-coincidences: NDUF/COX/ND -- look like a family by name, ZERO homology (complex
      subunits, not paralogs)
    - retrocopies: processed pseudogene vs parent (homologous but resolvable/dispersed)
Each is verified against the genome cDNA homology so the panel is trustworthy, unlike the
over-merged all-vs-all. Run: python bench/family_def_airtight_panel.py

SCOPE (§6am) — THIS PANEL NEVER INVOKES THE PIPELINE. Its only evidence is ANNOTATION-derived
cDNA all-vs-all homology (rep_ava.tsv) plus gene coordinates, so it is structurally BLIND to
every RUSTLE_* flag and to any collapse/representative/gate change: both arms of such a flag
return identical rows BY CONSTRUCTION, and "identical" here is not a verdict. A PASS printed
below is NOT citable as evidence about a flag, a binary, or any pipeline output — it adjudicates
the ANNOTATION-level family definition only.
"""
import collections
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import GENES_BED
from family_def_read_filters import AVA, dna_homology

# §6am: the previous default (/mnt/c/Users/jfris/Desktop/GGO_genomic.gff) DOES NOT EXIST, so the
# coordinate fallback below was dead code and any symbol absent from GENES_BED stayed silently
# uncoordinated. Gorilla-native RefSeq GFF (GCF_029281585.2); overridable via GGO_GFF.
GFF = os.environ.get("GGO_GFF", "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff")

# curated, named, well-characterized real families + counterexamples (gene symbols as in the
# gorilla RefSeq annotation; verified below against cDNA homology)
REAL = [   # (family, [member gene symbols])
    ("RABL2", ["RABL2A", "RABL2B"]),
    ("APOBEC3", ["APOBEC3D", "APOBEC3F"]),
    ("RFPL", ["RFPL1", "RFPL2", "RFPL3"]),
]
# real paralog pairs whose copy is a PROCESSED PSEUDOGENE (no cDNA) -> invisible to cDNA-truth.
# Kept as an explicit LIMITATION case, not a failed example.
PSEUDO = [("GBA", ["GBA1", "GBAP1"])]
COUNTER = [
    ("CASP8~FLACC1 (domain-sharer)", ["CASP8", "FLACC1"]),
    ("ASDURF~ASNSD1 (domain-sharer)", ["ASDURF", "ASNSD1"]),
    ("CREB1~METTL21A (domain-sharer)", ["CREB1", "METTL21A"]),
    ("GPR39~LYPD1 (domain-sharer)", ["GPR39", "LYPD1"]),
    ("NDUFS (name-coincidence)", ["NDUFS1", "NDUFS2", "NDUFS3"]),
    ("NDUFA (name-coincidence)", ["NDUFA2", "NDUFA5", "NDUFA10"]),
    ("COX (name-coincidence)", ["COX1", "COX2", "COX3"]),
]


def main():
    # Guard (§6am): hard abort BEFORE any scoring or output truncation when an evidence file is
    # missing — prevents a missing input from degrading into an empty evidence set on which the
    # negative criteria below pass by default.
    for label, path in (("gene coordinates (GENES_BED)", GENES_BED),
                        ("cDNA all-vs-all homology (AVA)", AVA)):
        if not os.path.exists(path):
            sys.exit(f"ABORT: {label} missing: {path}")

    coord = {}
    for src in (GENES_BED,):
        with open(src) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 4:
                    coord.setdefault(p[3], (p[0], int(p[1]), int(p[2])))
    # supplement coords from the GFF for symbols not in GENES_BED
    need = {g for _, gs in REAL + COUNTER for g in gs if g not in coord}
    if need:
        # Guard (§6am): the coordinate fallback is only sound if the GFF is actually there. With
        # the old dead path this loop raised nothing today only because `need` happened to be
        # empty; an uncoordinated member makes its pairs undecidable and scores as a free pass.
        if not os.path.exists(GFF):
            sys.exit(f"ABORT: {len(need)} panel symbol(s) absent from GENES_BED "
                     f"({', '.join(sorted(need))}) and the GFF is missing: {GFF}")
        for line in open(GFF):
            if line.startswith("#") or "\tgene\t" not in line:
                continue
            f = line.split("\t")
            m = re.search(r"Name=([^;]+)", f[8])
            if m and m.group(1) in need:
                coord[m.group(1)] = (f[0], int(f[3]), int(f[4]))

    Hd, _ = dna_homology()

    # POSITIVE CONTROL (§6am): every counterexample criterion is a NEGATIVE ("zero homology",
    # "no airtight edge") and therefore passes BY DEFAULT on an empty evidence set. Require the
    # homology evidence to be loaded AND to demonstrably cover this panel's own gene namespace
    # before any verdict is printed — otherwise a null result is vacuous, not a finding.
    panel_genes = {g for _, gs in REAL + COUNTER + PSEUDO for g in gs}
    n_touch = sum(1 for a, b in Hd if a in panel_genes or b in panel_genes)
    if not Hd or n_touch == 0:
        sys.exit(f"ABORT: homology evidence not loaded — {len(Hd)} pairs total, {n_touch} touching "
                 f"panel genes, from {AVA}. A negative verdict on an empty evidence set is vacuous.")
    print(f"[positive control] homology evidence: {len(Hd)} pairs, {n_touch} touch a panel gene")

    def overlaps(a, b):
        if a not in coord or b not in coord:
            return False
        (ca, sa, ea), (cb, sb, eb) = coord[a], coord[b]
        return ca == cb and min(ea, eb) - max(sa, sb) > 0

    def paralog_edge(a, b):
        """AIRTIGHT paralog edge: homologous + reciprocally covered + DISJOINT loci.
        The disjoint condition rejects self-alignment co-location artifacts (nesting,
        overlap, read-through) that inflate id to ~1.0 between non-paralogous neighbours."""
        r = Hd.get((a, b) if a < b else (b, a))
        if r is None:
            return (0.0, 0.0, False, overlaps(a, b))
        idv, mc = r.get("id", 0), min(r["cov_a"], r["cov_b"])
        ov = overlaps(a, b)
        is_edge = idv >= 0.90 and mc >= 0.30 and not ov
        return (idv, mc, is_edge, ov)

    def grp(gs):
        ps = [(gs[i], gs[j]) for i in range(len(gs)) for j in range(i + 1, len(gs))]
        ev = [paralog_edge(a, b) for a, b in ps]
        present = [g for g in gs if g in coord]
        ids = [e[0] for e in ev]
        mcs = [e[1] for e in ev]
        n_edge = sum(1 for e in ev if e[2])
        n_overlap = sum(1 for e in ev if e[3])
        return dict(pres=len(present), tot=len(gs),
                    mn_id=round(min(ids), 2) if ids else 0, mn_cov=round(min(mcs), 2) if mcs else 0,
                    n_edge=n_edge, n_pairs=len(ps), n_overlap=n_overlap)

    print("=== AIRTIGHT PANEL — examples + counterexamples ===")
    print("    airtight paralog edge = id>=0.90 AND min recip-cov>=0.30 AND DISJOINT loci\n")
    print("REAL FAMILIES (every member-pair an airtight paralog edge):")
    print(f"  {'family':22} {'present':>7} {'min_id':>6} {'min_cov':>7} {'edges':>7}  verdict")
    rows = []
    for fam, gs in REAL:
        g = grp(gs)
        ok = g["pres"] >= 2 and g["n_edge"] == g["n_pairs"] and g["n_pairs"] > 0
        print(f"  {fam:22} {g['pres']}/{g['tot']:<5} {g['mn_id']:>6.2f} {g['mn_cov']:>7.2f} "
              f"{g['n_edge']}/{g['n_pairs']:<5}  {'AIRTIGHT real' if ok else 'check (missing/low)'}")
        rows.append(dict(case=fam, kind="real", airtight=ok, **g))
    print("\nCOUNTEREXAMPLES — two airtight non-family modes:")
    print("  (A) co-location artifact: high id but loci OVERLAP -> homology is self-alignment, not paralogy")
    print("  (B) name-coincidence: grouped by name, ZERO homology (complex subunits)")
    print(f"  {'case':32} {'present':>7} {'min_id':>6} {'overlap':>7} {'edges':>7}  verdict")
    for fam, gs in COUNTER:
        g = grp(gs)
        if "name-coincidence" in fam:
            ok, why = (g["mn_id"] < 0.3 and g["n_edge"] == 0), "name, no homology"
        else:  # co-location: pairs overlap and NONE is an airtight edge
            ok, why = (g["n_overlap"] > 0 and g["n_edge"] == 0), "overlap artifact"
        print(f"  {fam:32} {g['pres']}/{g['tot']:<5} {g['mn_id']:>6.2f} {g['n_overlap']}/{g['n_pairs']:<5} "
              f"{g['n_edge']}/{g['n_pairs']:<5}  {'AIRTIGHT non-family ('+why+')' if ok else 'check'}")
        rows.append(dict(case=fam, kind="counter", airtight=ok, **g))

    print("\nLIMITATION — real paralog, but copy is a PROCESSED PSEUDOGENE (no cDNA -> invisible to cDNA-truth):")
    for fam, gs in PSEUDO:
        g = grp(gs)
        print(f"  {fam:22} {g['pres']}/{g['tot']:<5} min_id={g['mn_id']:.2f}  "
              f"-> not recoverable from cDNA homology (GBAP1 has no transcript); needs DNA-level evidence")
        rows.append(dict(case=fam, kind="pseudo_limitation", airtight=False, **g))

    nr = sum(1 for r in rows if r["kind"] == "real" and r["airtight"])
    nc = sum(1 for r in rows if r["kind"] == "counter" and r["airtight"])
    print(f"\n  panel: {nr}/{len(REAL)} real families airtight, {nc}/{len(COUNTER)} counterexamples airtight"
          f" (+{len(PSEUDO)} pseudogene-limitation case)")
    print("  KEY: the DISJOINT-loci condition is what makes homology airtight — without it the")
    print("       cDNA self-alignment bridges overlapping/nested neighbours at id~1.0 (the over-merge).")

    # Guard (§6am): abort BEFORE truncating the committed TSV/JSON when the panel did not fully
    # resolve. A case whose members were not all located has undecidable pairs (no coords ->
    # "no overlap", no Hd record -> "no edge"), which is a free pass for the counterexamples;
    # a partial panel must exit nonzero instead of overwriting the recorded result.
    adjudicated = [r for r in rows if r["kind"] in ("real", "counter")]
    unresolved = [r for r in adjudicated if r["pres"] < r["tot"] or r["n_pairs"] == 0]
    print(f"  resolved {len(adjudicated) - len(unresolved)}/{len(REAL) + len(COUNTER)} adjudicated "
          f"cases (every member located, >=1 pair)")
    if unresolved:
        sys.exit("ABORT: panel under-resolved, outputs NOT written — "
                 + "; ".join(f"{r['case']}: {r['pres']}/{r['tot']} members located, "
                             f"{r['n_pairs']} pairs" for r in unresolved))

    # reusable TSV: case  kind  members  airtight  min_id  min_cov  n_edge/n_pairs  n_overlap
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "family_def_airtight_panel.tsv"), "w") as out:
        out.write("case\tkind\tmembers\tairtight\tmin_id\tmin_cov\tedges\tpairs\toverlaps\n")
        for (fam, gs), r in zip(REAL + COUNTER + PSEUDO,
                                [x for x in rows]):
            out.write(f"{r['case']}\t{r['kind']}\t{','.join(gs)}\t{int(r['airtight'])}\t"
                      f"{r['mn_id']}\t{r['mn_cov']}\t{r['n_edge']}\t{r['n_pairs']}\t{r['n_overlap']}\n")
    json.dump(dict(rows=rows), open(os.path.join(here, "family_def_airtight_panel.json"), "w"), indent=2)
    print("  [+] wrote family_def_airtight_panel.{tsv,json}")


if __name__ == "__main__":
    main()
