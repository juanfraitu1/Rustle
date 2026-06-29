#!/usr/bin/env python3
"""dna_support.py — PROTOTYPE of the DNA heritability leg for the gene-conversion vs RT-switch
discriminator (measure-first before any Rust port).

The catalog (dna_catalog/DNFAM*.json) is REF0-CENTRIC: matrix[label][pos] = base, where
`pos = genomic - ref0_start` (offset into the longest copy's frame), sparse (absent => matches ref0).
A gene conversion is HERITABLE -> a DNA copy whose allele vector is itself a MOSAIC of two others
(matches copy X 5' of a breakpoint, copy Y 3' of it). An RT/template switch is RNA-only -> absent.

dna_support(chrom, bp_genomic) -> ("supported"|"absent"|"unchecked", detail):
  * find the family whose chrom matches and whose ref0 interval CONTAINS bp_genomic (clean projection
    only there; non-ref0 copies would need re-alignment -> "unchecked").
  * project bp -> ref0-local; among the family's DNA labels, look for a mosaic (historical conversion)
    whose breakpoint is near bp. Found -> "supported" (Some(true)); family present, none -> "absent"
    (Some(false)); no covering family -> "unchecked" (None).

This script MEASURES: (A) does the DNA conversion signal exist across families? (B) ref0 coverage of
the genome (how often a breakpoint is checkable at all)?  Run: python bench/mosaic_discriminator/dna_support.py
"""
import json
import os
import glob
import sys
from collections import defaultdict

CAT_DIR = "/home/juanfra/winloci_scratch/dna_catalog"
FAI = "/home/juanfra/winloci_scratch/GGO.fasta.fai"

# ---- mosaic detection on a DNA per-label allele matrix (mirror of the RNA detect_mosaic logic) ----

def label_alleles(matrix, label, cols):
    """allele vector for `label` over the ordered `cols` (ref0-local positions). Absent => 'R' (ref0)."""
    m = matrix.get(label, {})
    return [m.get(str(c), m.get(c, "R")) for c in cols]


def find_dna_mosaic(matrix, cols, min_tract=4, min_sites=8, max_err=1, max_labels=10):
    """Is some label a contiguous MOSAIC of two others across `cols`?  Cheaper O(labels^2 * splits):
    for each label L and each candidate split s, L is a mosaic iff its best-matching OTHER label on the
    LEFT differs from its best-matching OTHER label on the RIGHT, each side near-clean (<= max_err) and
    >= min_tract sites, and those two parents actually DIFFER across the breakpoint. Returns the best
    (mosaic_label, copy_x, copy_y, break_col_index, (lt,rt), err) or None. `cols` ascending (ref0-local).
    """
    if len(cols) < min_sites:
        return None
    labels = list(matrix.keys())[:max_labels]
    av = {lab: label_alleles(matrix, lab, cols) for lab in labels}
    n = len(cols)
    best = None
    # candidate split points: every position (n is small after subsampling)
    for L in labels:
        vL = av[L]
        others = [o for o in labels if o != L]
        for s in range(min_tract, n - min_tract + 1):
            # best left parent and best right parent for L
            def best_parent(rng):
                bo, be = None, 10**9
                for O in others:
                    e = sum(1 for i in rng if vL[i] != av[O][i])
                    if e < be:
                        bo, be = O, e
                return bo, be
            X, le = best_parent(range(s))
            Y, re_ = best_parent(range(s, n))
            if X is None or Y is None or X == Y:
                continue
            if le > max_err or re_ > max_err:
                continue
            # parents must differ across the breakpoint (else L isn't a recombinant of two distinct copies)
            if av[X][s:] == av[Y][s:]:
                continue
            cand = (L, X, Y, s, (s, n - s), le + re_)
            if best is None or cand[5] < best[5]:
                best = cand
    return best


def family_cols(matrix):
    """ordered union of all substitution positions (ref0-local, int) across labels."""
    cols = set()
    for m in matrix.values():
        for p in m:
            cols.add(int(p))
    return sorted(cols)


# ---- catalog index ----

import re
_META = re.compile(rb'"ref0":\s*"([^"]+)".*?"chrom":\s*"([^"]+)".*?"iv":\s*\[\s*(\d+)\s*,\s*(\d+)', re.S)


def load_index(cat_dir=CAT_DIR):
    """metadata only (FAST): regex-extract ref0/chrom/iv from each file HEAD (they precede the huge
    `matrix` in the dump order), avoiding a full json parse of 1130 large files."""
    idx = {}
    by_chrom = defaultdict(list)
    for path in glob.glob(os.path.join(cat_dir, "DNFAM*.json")):
        with open(path, "rb") as fh:
            head = fh.read(4096)
        m = _META.search(head)
        if not m:  # fall back to a full parse if the head is laid out differently
            try:
                d = json.load(open(path))
                ref0, chrom, iv = d["ref0"], d["chrom"], tuple(d["iv"])
            except Exception:
                continue
        else:
            ref0 = m.group(1).decode()
            chrom = m.group(2).decode()
            iv = (int(m.group(3)), int(m.group(4)))
        fid = os.path.basename(path)[:-5]
        idx[fid] = (chrom, iv, ref0, path, None)  # n_labels lazily None
        by_chrom[chrom].append((iv[0], iv[1], fid))
    for c in by_chrom:
        by_chrom[c].sort()
    return idx, by_chrom


def covering_family(by_chrom, chrom, bp_genomic):
    for (s, e, fid) in by_chrom.get(chrom, []):
        if s <= bp_genomic < e:
            return fid
    return None


def dna_support(idx, by_chrom, chrom, bp_genomic, win=60):
    """POSITIVE-ONLY DNA corroboration -> ('supported'|'unchecked', detail). Returns 'supported'
    (=> classify_event dna_supported=Some(true)) ONLY when a DNA mosaic (historical-conversion
    signature) coincides with the breakpoint; EVERYTHING ELSE is 'unchecked' (=> None), NEVER a veto.
    MEASURED RATIONALE: catalog "absent" is unreliable (sparse, ref0-centric, localized mosaics) — using
    it as Some(false) would wrongly downgrade real conversions. See module docstring."""
    fid = covering_family(by_chrom, chrom, bp_genomic)
    if fid is None:
        return ("unchecked", "no covering ref0 interval")
    chrom_, iv, ref0, path, _ = idx[fid]
    d = json.load(open(path))
    matrix = d["matrix"]
    if len(matrix) < 3:
        return ("unchecked", f"{fid}: <3 DNA copies")
    bp_local = bp_genomic - iv[0]
    near = [c for c in family_cols(matrix) if abs(c - bp_local) <= 4000]  # PSV cols within 4kb of the bp
    if len(near) < 8:
        return ("unchecked", f"{fid}: <8 DNA PSVs within 4kb of breakpoint")
    cols = near[:: max(1, len(near) // 40)][:40]   # subsample to bound the scan, keeping the bp central
    mos = find_dna_mosaic(matrix, cols)
    if mos is None:
        return ("unchecked", f"{fid}: DNA PSVs present, no mosaic copy near breakpoint")
    L, X, Y, s, tracts, err = mos
    break_col = cols[s] if s < len(cols) else cols[-1]
    coincides = abs(break_col - bp_local) <= win
    if coincides:
        return ("supported", f"{fid}: DNA mosaic {L}={X}|{Y} break@local{break_col} coincides (err={err})")
    return ("unchecked", f"{fid}: DNA mosaic exists but break@local{break_col} != bp_local{bp_local}")


# ---- measurements ----

def chrom_lengths(fai=FAI):
    L = {}
    if os.path.exists(fai):
        for line in open(fai):
            p = line.split("\t")
            L[p[0]] = int(p[1])
    return L


def main():
    idx, by_chrom = load_index()
    clen = chrom_lengths()
    print(f"=== DNA catalog: {len(idx)} families on {len(by_chrom)} contigs ===\n")

    # (B) ref0-interval COVERAGE of the genome
    print("--- (B) ref0-interval coverage (how often a random breakpoint is checkable) ---")
    tot_span = 0
    tot_genome = sum(clen.values()) if clen else 0
    for c, ivs in sorted(by_chrom.items()):
        # merge overlapping intervals
        merged = []
        for s, e, _ in sorted(ivs):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        span = sum(e - s for s, e in merged)
        tot_span += span
    print(f"  total ref0-covered span: {tot_span:,} bp over {tot_genome:,} bp genome "
          f"= {100*tot_span/tot_genome:.2f}% (a uniform-random breakpoint is checkable this often)")
    print(f"  (events cluster in paralog regions, so the EFFECTIVE rate for conversion breakpoints is higher)\n")

    # (A) does the DNA conversion SIGNAL exist? scan a sample of families for a DNA mosaic.
    print("--- (A) DNA heritable-conversion signal presence (sample of families) ---", flush=True)
    n_sample = int(os.environ.get("DNA_SAMPLE", "40"))
    # iterate MID-SIZED families (skip empty tiny ones and the giant slow ones) until n_sample non-empty
    by_size = sorted(idx.keys(), key=lambda f: os.path.getsize(idx[f][3]))
    pool = by_size[len(by_size) // 3 : -5]   # drop smallest third (empty) and the few biggest (slow)
    n_checked = n_mosaic = 0
    examples = []
    k = 0
    for fid in pool:
        if n_checked >= n_sample:
            break
        chrom_, iv, ref0, path, _ = idx[fid]
        try:
            d = json.load(open(path))
        except Exception:
            continue
        matrix = d["matrix"]
        if len(matrix) < 3:
            continue
        k += 1
        if k % 10 == 1:
            print(f"    [{n_checked}/{n_sample}] {fid} labels={len(matrix)} ...", flush=True)
        if len(matrix) < 3:
            continue
        cols = family_cols(matrix)
        if len(cols) < 8:
            continue
        n_checked += 1
        sub = cols[:: max(1, len(cols) // 30)][:30]   # subsample cols (mosaic scan is O(labels^2 * cols))
        mos = find_dna_mosaic(matrix, sub, max_labels=10)
        if mos is not None:
            n_mosaic += 1
            if len(examples) < 8:
                bp_genomic = iv[0] + sub[mos[3]]   # mosaic breakpoint -> genomic (for the positive-direction test)
                examples.append((fid, chrom_, ref0, len(cols), mos[:4], bp_genomic))
    print(f"  checked {n_checked} families (>=3 labels, >=8 PSV cols); "
          f"{n_mosaic} ({100*n_mosaic/max(1,n_checked):.0f}%) show a DNA MOSAIC (historical-conversion signature)")
    for fid, c, ref0, ncol, mo, bpg in examples:
        print(f"    {fid} {c} ref0={ref0} cols={ncol}: mosaic {mo[0]}={mo[1]}|{mo[2]} @colidx{mo[3]} (genomic {bpg})")
    print()

    # (D) POSITIVE-direction validation: query dna_support AT a real DNA mosaic breakpoint → 'supported'.
    print("--- (D) dna_support() at the ACTUAL DNA mosaic breakpoints (positive direction) ---", flush=True)
    n_sup = 0
    for fid, c, ref0, ncol, mo, bpg in examples:
        status, detail = dna_support(idx, by_chrom, c, bpg)
        n_sup += status == "supported"
        print(f"  {fid} {c}:{bpg} -> {status}  ({detail})", flush=True)
    print(f"  => {n_sup}/{len(examples)} confirmed 'supported' at the true breakpoint", flush=True)

    # (C) at an ARBITRARY position (midpoint) → mostly 'absent': shows why 'absent' is WEAK negative
    # evidence (a localized DNA conversion isn't at the midpoint), i.e. the leg must be POSITIVE-ONLY.
    print("\n--- (C) dna_support() at arbitrary midpoints (absent = weak, do NOT use as a veto) ---", flush=True)
    for fid, c, ref0, ncol, mo, bpg in examples[:4]:
        chrom_, iv, ref0_, path, _ = idx[fid]
        mid = (iv[0] + iv[1]) // 2
        status, detail = dna_support(idx, by_chrom, chrom_, mid)
        print(f"  {fid} {chrom_}:{mid}(midpoint) -> {status}  ({detail})", flush=True)


if __name__ == "__main__":
    main()
