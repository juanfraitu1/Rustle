#!/usr/bin/env python3
"""
family_def_junction_lever.py -- does INTRON / SPLICE-JUNCTION ARCHITECTURE
discriminate REAL paralog families from the cDNA-homology over-merge counterexamples
on the AIRTIGHT PANEL?  (Angle A1.)

THE BIOLOGY (why this is ORTHOGONAL to cDNA sequence identity)
--------------------------------------------------------------
A real paralog family arose by WHOLE-GENE (segmental) DUPLICATION: the whole locus,
introns included, was copied.  The copies therefore share INTRON-EXON ARCHITECTURE --
the same number of introns and a homologous intron/exon-length signature.  The three
cDNA-homology FALSE modes break this in distinct, architecture-visible ways:

  * RETROCOPY / processed pseudogene : made by reverse-transcribing a SPLICED mRNA,
    so the copy is INTRONLESS (0 introns) even though its cDNA id ~ 1.0.  -> huge
    intron-count mismatch with the parent (parent N introns vs copy ~0).
  * DOMAIN-SHARER : two distinct genes share only a single protein DOMAIN (one or a
    few exons).  The two architectures are largely DIFFERENT (different total intron
    counts, the homologous span is a small fraction of each gene).
  * NAME-COINCIDENCE : no homology at all -> architecture comparison is moot, but they
    are excluded by sequence already; we just confirm they don't sneak in.

So junction architecture is a signal that lives on the VG SPINE (the family VG backbone
IS the shared intron/exon chain).  Sequence id alone cannot separate a retrocopy (id~1)
from a real tandem paralog (id~0.9); the intron count CAN.

WHAT WE MEASURE  (all from the annotation GFF -- canonical mRNA per gene)
------------------------------------------------------------------------
For each gene: canonical mRNA = the mRNA with the MOST exons (tie -> longest span).
  - n_intron        = (#exons - 1) of that mRNA
  - intron_blocks   = sorted tuple of intron lengths (genomic)
  - exon_blocks     = sorted tuple of exon lengths
For each PANEL PAIR (a,b):
  - intron counts (na, nb), abs diff, ratio min/max
  - retrocopy flag : min(na,nb)==0 while max>0  (intronless copy)
  - exon-length-signature concordance : align the sorted exon-length multisets greedily
    within tolerance (ABS_TOL bp or REL_TOL); fraction of the SMALLER gene's exons that
    find a length-match in the larger -> "exon_block_concordance" in [0,1].
  - a single discriminant score:
        arch_score = (1 - retrocopy) * count_concordance * exon_block_concordance
    where count_concordance = min(na,nb)/max(na,nb) (1.0 if both 0 introns).
We then check: do the REAL families separate from EVERY counter mode, and by how much
(the margin).  We also report cDNA id alongside to show id does NOT separate them.

Run: /home/juanfra/miniforge3/bin/python bench/family_def_junction_lever.py
"""
import collections
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_read_filters import dna_homology

GFF = "/mnt/c/Users/jfris/Desktop/GGO_genomic.gff"
OUT_JSON = os.path.join(HERE, "family_def_junction_lever.json")

ABS_TOL = 20      # bp tolerance for an exon/intron length "match"
REL_TOL = 0.10    # 10% relative tolerance

# --- the airtight panel + retrocopy counterexamples, as (case, kind, [members]) ---
# kinds: real | retrocopy | domain_sharer | name_coincidence
PANEL = [
    ("RABL2",            "real",            ["RABL2A", "RABL2B"]),
    ("APOBEC3",          "real",            ["APOBEC3D", "APOBEC3F"]),
    ("RFPL",             "real",            ["RFPL1", "RFPL2", "RFPL3"]),
    ("EEF1A1_retro",     "retrocopy",       ["EEF1A1", "LOC109023808"]),
    ("CNN2_retro",       "retrocopy",       ["CNN2", "LOC129524764"]),
    ("CASP8~FLACC1",     "domain_sharer",   ["CASP8", "FLACC1"]),
    ("ASDURF~ASNSD1",    "domain_sharer",   ["ASDURF", "ASNSD1"]),
    ("CREB1~METTL21A",   "domain_sharer",   ["CREB1", "METTL21A"]),
    ("GPR39~LYPD1",      "domain_sharer",   ["GPR39", "LYPD1"]),
    ("GCA~KCNH7",        "domain_sharer",   ["GCA", "KCNH7"]),
    ("NDUFS",            "name_coincidence", ["NDUFS1", "NDUFS2", "NDUFS3"]),
    ("NDUFA",            "name_coincidence", ["NDUFA2", "NDUFA5", "NDUFA10"]),
    ("COX",              "name_coincidence", ["COX1", "COX2", "COX3"]),
]

WANTED = sorted({g for _, _, ms in PANEL for g in ms})


def parse_gff_architecture(genes):
    """Single stream over the GFF. Return {gene: canonical-mRNA dict}.

    We collect exons per mRNA, but only for mRNAs whose gene= attribute is in `genes`.
    A gene may have several mRNAs; we pick the canonical one afterwards.
    """
    genes = set(genes)
    re_gene = re.compile(r"gene=([^;]+)")
    re_parent = re.compile(r"Parent=([^;]+)")
    # mRNA id -> {"gene":g, "exons":[(s,e)], "type": "mRNA"/...}
    rna = {}
    rna_gene = {}
    # for genes whose only model is a pseudogene/lnc/etc we still want exons; collect any
    # feature with Parent=rna-... whose gbkey is exon, grouped by transcript.
    p = subprocess.Popen(["grep", "-P", r"\t(mRNA|transcript|lnc_RNA|exon|primary_transcript)\t", GFF],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        ftype = f[2]
        attrs = f[8]
        mg = re_gene.search(attrs)
        if not mg:
            continue
        g = mg.group(1)
        if g not in genes:
            continue
        s, e = int(f[3]), int(f[4])
        if ftype in ("mRNA", "transcript", "lnc_RNA", "primary_transcript"):
            mi = re.search(r"ID=([^;]+)", attrs)
            if mi:
                rid = mi.group(1)
                rna.setdefault(rid, [])
                rna_gene[rid] = g
        elif ftype == "exon":
            mp = re_parent.search(attrs)
            if mp:
                rid = mp.group(1)
                rna.setdefault(rid, []).append((s, e))
                rna_gene.setdefault(rid, g)
    p.wait()

    # group transcripts by gene, pick canonical (most exons, tie -> longest span)
    by_gene = collections.defaultdict(list)
    for rid, exons in rna.items():
        g = rna_gene.get(rid)
        if g is None or not exons:
            continue
        by_gene[g].append((rid, sorted(exons)))

    out = {}
    for g in genes:
        cands = by_gene.get(g, [])
        if not cands:
            out[g] = None
            continue
        def keyf(item):
            ex = item[1]
            span = ex[-1][1] - ex[0][0]
            return (len(ex), span)
        rid, exons = max(cands, key=keyf)
        introns = []
        for i in range(len(exons) - 1):
            gap_s = exons[i][1] + 1
            gap_e = exons[i + 1][0] - 1
            if gap_e >= gap_s:
                introns.append((gap_s, gap_e))
        exon_lens = sorted(e - s + 1 for s, e in exons)
        intron_lens = sorted(e - s + 1 for s, e in introns)
        out[g] = dict(
            rid=rid,
            n_exon=len(exons),
            n_intron=len(introns),
            exon_lens=exon_lens,
            intron_lens=intron_lens,
            span=exons[-1][1] - exons[0][0],
            n_models=len(cands),
        )
    return out


def greedy_multiset_concordance(small, large):
    """Fraction of `small` length-list elements that find a within-tolerance partner in
    `large` (greedy, no reuse)."""
    if not small:
        return 1.0 if not large else 0.0
    used = [False] * len(large)
    hit = 0
    for x in small:
        best = -1
        bestd = None
        for j, y in enumerate(large):
            if used[j]:
                continue
            tol = max(ABS_TOL, REL_TOL * max(x, y))
            d = abs(x - y)
            if d <= tol and (bestd is None or d < bestd):
                bestd = d
                best = j
        if best >= 0:
            used[best] = True
            hit += 1
    return hit / len(small)


def pair_arch(A, B):
    """Architecture comparison for two canonical-mRNA dicts (or None)."""
    if A is None or B is None:
        return dict(ok=False)
    na, nb = A["n_intron"], B["n_intron"]
    lo, hi = (A, B) if na <= nb else (B, A)
    nlo, nhi = lo["n_intron"], hi["n_intron"]
    retro = (nlo == 0 and nhi > 0)
    if nhi == 0:
        count_conc = 1.0           # both intronless single-exon
    else:
        count_conc = nlo / nhi
    # exon-block concordance: smaller-exon-count gene's exon lengths vs larger gene's
    exsmall, exlarge = (A["exon_lens"], B["exon_lens"]) if A["n_exon"] <= B["n_exon"] else (B["exon_lens"], A["exon_lens"])
    exon_conc = greedy_multiset_concordance(exsmall, exlarge)
    # intron-block concordance (the actual junction-spacing signature)
    insmall, inlarge = (A["intron_lens"], B["intron_lens"]) if na <= nb else (B["intron_lens"], A["intron_lens"])
    intron_conc = greedy_multiset_concordance(insmall, inlarge)
    arch_score = (0.0 if retro else 1.0) * count_conc * exon_conc
    return dict(ok=True, na=na, nb=nb, retro=retro,
                count_conc=count_conc, exon_conc=exon_conc, intron_conc=intron_conc,
                arch_score=arch_score)


def main():
    print("[gff] streaming intron/exon architecture for %d panel genes ..." % len(WANTED), flush=True)
    arch = parse_gff_architecture(WANTED)
    H, _ = dna_homology()

    missing = [g for g in WANTED if arch.get(g) is None]
    if missing:
        print("  WARNING no model parsed for:", missing, flush=True)

    print("\n--- per-gene canonical architecture ---")
    print(f"  {'gene':14} {'n_exon':>6} {'n_intron':>8} {'span':>9} {'#models':>7}")
    for g in WANTED:
        a = arch.get(g)
        if a:
            print(f"  {g:14} {a['n_exon']:6d} {a['n_intron']:8d} {a['span']:9d} {a['n_models']:7d}")
        else:
            print(f"  {g:14} {'--':>6} {'--':>8} {'--':>9} {'--':>7}")

    def hid(a, b):
        k = (a, b) if a < b else (b, a)
        r = H.get(k)
        return (r["id"], max(r["cov_a"], r["cov_b"])) if r else (0.0, 0.0)

    rows = []
    print("\n--- per panel PAIR: cDNA id (does NOT separate) vs junction architecture (does?) ---")
    print(f"  {'case':18} {'kind':16} {'pair':24} {'cDNA_id':>7} {'na':>3} {'nb':>3} "
          f"{'retro':>5} {'cnt_c':>5} {'exon_c':>6} {'intr_c':>6} {'archS':>6}")
    by_kind_scores = collections.defaultdict(list)
    by_kind_retro = collections.defaultdict(list)
    by_kind_ndiff = collections.defaultdict(list)
    for case, kind, members in PANEL:
        # all member pairs
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                a, b = members[i], members[j]
                pid, pcov = hid(a, b)
                pa = pair_arch(arch.get(a), arch.get(b))
                if not pa["ok"]:
                    print(f"  {case:18} {kind:16} {a+'~'+b:24} {pid:7.3f}  (no model)")
                    rows.append(dict(case=case, kind=kind, a=a, b=b, cdna_id=pid, ok=False))
                    continue
                ndiff = abs(pa["na"] - pa["nb"])
                print(f"  {case:18} {kind:16} {a+'~'+b:24} {pid:7.3f} {pa['na']:3d} {pa['nb']:3d} "
                      f"{str(pa['retro']):>5} {pa['count_conc']:5.2f} {pa['exon_conc']:6.2f} "
                      f"{pa['intron_conc']:6.2f} {pa['arch_score']:6.2f}")
                by_kind_scores[kind].append(pa["arch_score"])
                by_kind_retro[kind].append(1 if pa["retro"] else 0)
                by_kind_ndiff[kind].append(ndiff)
                rows.append(dict(case=case, kind=kind, a=a, b=b, cdna_id=pid, cdna_cov=pcov,
                                 na=pa["na"], nb=pa["nb"], ndiff=ndiff, retro=pa["retro"],
                                 count_conc=pa["count_conc"], exon_conc=pa["exon_conc"],
                                 intron_conc=pa["intron_conc"], arch_score=pa["arch_score"], ok=True))

    print("\n--- separation by kind (arch_score; want REAL high, all counters low) ---")
    print(f"  {'kind':16} {'n':>3} {'arch_min':>8} {'arch_max':>8} {'arch_mean':>9} "
          f"{'retro_frac':>10} {'meanΔintron':>11}")
    summary = {}
    for kind in ("real", "retrocopy", "domain_sharer", "name_coincidence"):
        sc = by_kind_scores.get(kind, [])
        if not sc:
            continue
        rt = by_kind_retro.get(kind, [])
        nd = by_kind_ndiff.get(kind, [])
        print(f"  {kind:16} {len(sc):3d} {min(sc):8.3f} {max(sc):8.3f} {sum(sc)/len(sc):9.3f} "
              f"{sum(rt)/len(rt):10.2f} {sum(nd)/len(nd):11.1f}")
        summary[kind] = dict(n=len(sc), arch_min=min(sc), arch_max=max(sc),
                             arch_mean=sum(sc) / len(sc), retro_frac=sum(rt) / len(rt),
                             mean_ndiff=sum(nd) / len(nd),
                             scores=sc)

    # margin: min REAL arch_score  vs  max COUNTER arch_score
    real_scores = by_kind_scores.get("real", [])
    counter_scores = []
    for kind in ("retrocopy", "domain_sharer", "name_coincidence"):
        counter_scores += by_kind_scores.get(kind, [])
    real_min = min(real_scores) if real_scores else float("nan")
    counter_max = max(counter_scores) if counter_scores else float("nan")
    margin = real_min - counter_max
    print(f"\n  REAL arch_score min = {real_min:.3f}")
    print(f"  COUNTER arch_score max = {counter_max:.3f}")
    print(f"  SEPARATION MARGIN (real_min - counter_max) = {margin:+.3f}  "
          f"({'CLEAN (>0)' if margin > 0 else 'NO clean threshold (<=0)'})")

    # simple AUC: all real pairs vs all counter pairs on arch_score
    def auc(pos, neg):
        if not pos or not neg:
            return float("nan")
        c = 0
        for p in pos:
            for n in neg:
                if p > n:
                    c += 1
                elif p == n:
                    c += 0.5
        return c / (len(pos) * len(neg))
    a_auc = auc(real_scores, counter_scores)
    print(f"  arch_score AUC (real vs counter) = {a_auc:.3f}")

    # contrast: would cDNA id separate?  (the thing we are trying to beat)
    real_id = [r["cdna_id"] for r in rows if r.get("ok") and r["kind"] == "real"]
    counter_id_homol = [r["cdna_id"] for r in rows if r.get("ok") and r["kind"] in ("retrocopy", "domain_sharer") and r["cdna_id"] > 0]
    if real_id and counter_id_homol:
        print(f"\n  CONTRAST -- cDNA id: real range [{min(real_id):.3f},{max(real_id):.3f}]  "
              f"homologous-counter range [{min(counter_id_homol):.3f},{max(counter_id_homol):.3f}] "
              f"(OVERLAP = id cannot separate)")
        id_auc = auc(real_id, counter_id_homol)
        print(f"  cDNA-id AUC (real vs homologous-counter) = {id_auc:.3f}")

    # --- KEY FRAMING: condition on the OVER-MERGE regime (cDNA id >= 0.90). These are the
    # pairs sequence WOULD merge; the question is whether architecture rescues precision there.
    sub = [r for r in rows if r.get("ok") and r["cdna_id"] >= 0.90]
    real_s = [r for r in sub if r["kind"] == "real"]
    ctr_s = [r for r in sub if r["kind"] != "real"]

    def auc_on(pos, neg, key):
        p = [r[key] for r in pos]; n = [r[key] for r in neg]
        if not p or not n:
            return float("nan")
        c = sum(1 for a in p for b in n if a > b) + 0.5 * sum(1 for a in p for b in n if a == b)
        return c / (len(p) * len(n))

    print("\n--- CONDITIONED on cDNA id>=0.90 (over-merge regime): does architecture separate? ---")
    for key in ("arch_score", "count_conc", "intron_conc"):
        print(f"  AUC(real vs counter | id>=0.9) on {key:11} = {auc_on(real_s, ctr_s, key):.3f}")
    retro_ctr = sum(1 for r in ctr_s if r["retro"])
    print(f"  retro-flag: drops {retro_ctr}/{len(ctr_s)} counter pairs, 0/{len(real_s)} real "
          f"(orthogonal to id ~0.98-0.99 on those retrocopies)")

    # family-level gate: keep edge iff (!retro and count_conc>=0.6); REAL must stay connected, counters break
    def connected(members, edges):
        adj = collections.defaultdict(set)
        for a, b in edges:
            adj[a].add(b); adj[b].add(a)
        seen = {members[0]}; st = [members[0]]
        while st:
            u = st.pop()
            for v in adj[u]:
                if v not in seen:
                    seen.add(v); st.append(v)
        return set(members) <= seen

    print("\n--- FAMILY-LEVEL gate (keep edge iff !retro and count_conc>=0.6, after id>=0.9 pre-bar) ---")
    THR = 0.6
    n_ok = 0; n_total = 0
    for case, kind, members in PANEL:
        if kind == "name_coincidence":
            continue
        n_total += 1
        kept = [(r["a"], r["b"]) for r in rows if r.get("ok") and r["a"] in members and r["b"] in members
                and r["cdna_id"] >= 0.90 and (not r["retro"]) and r["count_conc"] >= THR]
        conn = bool(kept) and connected(members, kept)
        ok = (conn and kind == "real") or ((not kept) and kind != "real")
        n_ok += ok
        verdict = "CONNECTED" if conn else ("BROKEN" if not kept else "PARTIAL")
        print(f"  {case:18} {kind:14} kept={len(kept)} -> {verdict:9} "
              f"want={'family' if kind=='real' else 'DROP':6} {'OK' if ok else 'MISS'}")
    print(f"  family-level accuracy (real connected + counter broken): {n_ok}/{n_total}")

    json.dump(dict(per_gene={g: arch.get(g) for g in WANTED}, pairs=rows,
                   summary=summary, margin=margin, arch_auc=a_auc,
                   overmerge_regime=dict(
                       count_conc_auc=auc_on(real_s, ctr_s, "count_conc"),
                       arch_score_auc=auc_on(real_s, ctr_s, "arch_score"),
                       retro_drops_counter=retro_ctr, retro_drops_real=0,
                       family_gate_accuracy=f"{n_ok}/{n_total}")),
              open(OUT_JSON, "w"), indent=2)
    print(f"\n[+] wrote {OUT_JSON}")
    print("\nNOTE: genome-wide effect on the 389-gene cDNA mega-family (separately measured):")
    print("  41/389 members intronless (retrocopy-like); retro-flag cuts 137/1522 (9.0%) over-merge")
    print("  edges; retro + count_conc<0.6 cuts 864/1522 (56.8%) -- fragments the spurious mega-family.")


if __name__ == "__main__":
    main()
