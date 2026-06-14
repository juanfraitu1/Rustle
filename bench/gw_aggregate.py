#!/usr/bin/env python3
"""Aggregate genome-wide three-way gffcompare scorecard.

Per chrom, parse:
  - par_C  (parity stats: VG vs StringTie)           -> Transcript-level Sn, Matching tx, Missed loci
  - par_C.<vg>.gtf.refmap                            -> ST transcripts NOT reproduced by VG (regressions)
  - gc_st_C.<st>.gtf.tmap / gc_vg_C.<vg>.gtf.tmap    -> NCBI ref_ids matched at STRICT {=,c} and BROAD {=,c,j}
Genome-wide totals, illustrative yield, and VG-only net-new gene list.
"""
import glob, os, re, sys
from collections import defaultdict

OUT = "/tmp/gw"

def find_done_chroms():
    return sorted(c[len(OUT)+1:-5] for c in glob.glob(os.path.join(OUT, "*.done")))

def parse_parity_stats(path):
    """Return (transcript_sn_pct, matching_tx, ref_mRNAs, query_mRNAs)."""
    sn = matching = ref_m = qry_m = None
    if not os.path.exists(path):
        return (None, None, None, None)
    with open(path) as fh:
        for line in fh:
            if line.lstrip().startswith("Transcript level:"):
                m = re.search(r"Transcript level:\s*([\d.]+)", line)
                if m: sn = float(m.group(1))
            elif "Matching transcripts:" in line:
                m = re.search(r"Matching transcripts:\s*(\d+)", line)
                if m: matching = int(m.group(1))
            elif "Reference mRNAs :" in line:
                m = re.search(r"Reference mRNAs :\s*(\d+)", line)
                if m: ref_m = int(m.group(1))
            elif "Query mRNAs :" in line:
                m = re.search(r"Query mRNAs :\s*(\d+)", line)
                if m: qry_m = int(m.group(1))
    return (sn, matching, ref_m, qry_m)

def refmap_reproduced_refs(path):
    """From a .refmap (ref-centric), return set of ref transcript ids that got an =/c match.
    refmap cols: ref_gene_id, ref_id, class_code, qry_id_list. class_code is per-ref best."""
    reproduced = set()
    allrefs = set()
    if not os.path.exists(path):
        return reproduced, allrefs
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            ref_id = f[1]
            cc = f[2]
            allrefs.add(ref_id)
            if cc in ("=", "c"):
                reproduced.add(ref_id)
    return reproduced, allrefs

def tmap_matched_refs(path, chrom):
    """From a query .tmap, return (strict_set, broad_set, gene_of) keyed by chrom:ref_id.
    strict = class in {=,c}; broad = {=,c,j}. gene_of maps key->ref_gene_id."""
    strict, broad = set(), set()
    gene_of = {}
    if not os.path.exists(path):
        return strict, broad, gene_of
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            ref_gene, ref_id, cc = f[0], f[1], f[2]
            if ref_id == "-" or ref_id == "":
                continue
            key = f"{chrom}:{ref_id}"
            if cc in ("=", "c", "j"):
                broad.add(key)
                gene_of[key] = ref_gene
                if cc in ("=", "c"):
                    strict.add(key)
    return strict, broad, gene_of

def tmap_path(prefix, chrom, tool):
    """tmap filename pattern: <prefix>.<basename-of-query-gtf>.tmap"""
    pat = os.path.join(OUT, f"{prefix}_{chrom}.*.tmap")
    hits = glob.glob(pat)
    return hits[0] if hits else None

def main():
    chroms = find_done_chroms()
    if not chroms:
        print("No completed chroms found in /tmp/gw (*.done).")
        sys.exit(1)

    rows = []
    # genome-wide accumulators
    gw = dict(st_tx=0, vg_tx=0,
              st_strict=set(), vg_strict=set(),
              st_broad=set(), vg_broad=set())
    gene_of_all = {}
    parity_matching_total = 0
    parity_ref_total = 0
    regressions_total = 0

    for c in chroms:
        # timing / tx counts
        st_tx = vg_tx = nreads = c_wall = r_wall = None
        tpath = os.path.join(OUT, f"{c}.timing")
        if os.path.exists(tpath):
            with open(tpath) as fh:
                parts = fh.read().split()
                if len(parts) >= 5:
                    c_wall, r_wall, st_tx, vg_tx, nreads = (int(x) for x in parts[:5])

        # parity
        par_stats = os.path.join(OUT, f"par_{c}")
        sn, matching, ref_m, qry_m = parse_parity_stats(par_stats)
        # ST-only regressions: ref (=ST) transcripts not reproduced by VG at =/c
        par_refmap = glob.glob(os.path.join(OUT, f"par_{c}.*.refmap"))
        reproduced = allrefs = set()
        if par_refmap:
            reproduced, allrefs = refmap_reproduced_refs(par_refmap[0])
        regressions = len(allrefs) - len(reproduced) if allrefs else None

        # NCBI matches
        st_tmap = tmap_path("gc_st", c, "st")
        vg_tmap = tmap_path("gc_vg", c, "vg")
        st_s, st_b, st_gene = tmap_matched_refs(st_tmap, c) if st_tmap else (set(), set(), {})
        vg_s, vg_b, vg_gene = tmap_matched_refs(vg_tmap, c) if vg_tmap else (set(), set(), {})

        gw["st_strict"] |= st_s; gw["vg_strict"] |= vg_s
        gw["st_broad"]  |= st_b; gw["vg_broad"]  |= vg_b
        gene_of_all.update(st_gene); gene_of_all.update(vg_gene)
        if st_tx: gw["st_tx"] += st_tx
        if vg_tx: gw["vg_tx"] += vg_tx
        if matching is not None: parity_matching_total += matching
        if ref_m is not None: parity_ref_total += ref_m
        if regressions is not None: regressions_total += regressions

        vg_only_s = vg_s - st_s
        vg_only_b = vg_b - st_b

        rows.append(dict(
            chrom=c, c_wall=c_wall, r_wall=r_wall, nreads=nreads,
            st_tx=st_tx, vg_tx=vg_tx,
            parity_sn=sn, parity_match=matching, parity_ref=ref_m, regressions=regressions,
            st_strict=len(st_s), vg_strict=len(vg_s), vgonly_strict=len(vg_only_s),
            st_broad=len(st_b), vg_broad=len(vg_b), vgonly_broad=len(vg_only_b),
        ))

    # genome-wide sets
    GW_st_s = gw["st_strict"]; GW_vg_s = gw["vg_strict"]
    GW_st_b = gw["st_broad"];  GW_vg_b = gw["vg_broad"]
    vgonly_s = GW_vg_s - GW_st_s
    vgonly_b = GW_vg_b - GW_st_b
    stonly_s = GW_st_s - GW_vg_s
    stonly_b = GW_st_b - GW_vg_b
    both_s = GW_st_s & GW_vg_s
    both_b = GW_st_b & GW_vg_b

    # ===== print per-chrom table =====
    print("=" * 140)
    print("PER-CHROM SCORECARD")
    print("=" * 140)
    hdr = (f"{'chrom':<14}{'reads':>9}{'rWall':>7}{'cWall':>7}"
           f"{'st_tx':>7}{'vg_tx':>7}{'par_Sn%':>8}{'regr':>6}"
           f"{'STstr':>7}{'VGstr':>7}{'VGo_s':>7}{'STbr':>7}{'VGbr':>7}{'VGo_b':>7}")
    print(hdr)
    print("-" * 140)
    for r in rows:
        psn = f"{r['parity_sn']:.1f}" if r['parity_sn'] is not None else "-"
        warn = "  <<SLOW" if (r['r_wall'] and r['r_wall'] > 360) else ""
        print(f"{r['chrom']:<14}{(r['nreads'] or 0):>9}{(r['r_wall'] or 0):>7}{(r['c_wall'] or 0):>7}"
              f"{(r['st_tx'] or 0):>7}{(r['vg_tx'] or 0):>7}{psn:>8}{(r['regressions'] or 0):>6}"
              f"{r['st_strict']:>7}{r['vg_strict']:>7}{r['vgonly_strict']:>7}"
              f"{r['st_broad']:>7}{r['vg_broad']:>7}{r['vgonly_broad']:>7}{warn}")

    # ===== genome-wide totals =====
    print()
    print("=" * 140)
    print("GENOME-WIDE TOTALS")
    print("=" * 140)
    print(f"Chroms processed                     : {len(rows)}")
    print(f"Total StringTie transcripts          : {gw['st_tx']}")
    print(f"Total rustle-VG transcripts          : {gw['vg_tx']}")
    print()
    print("--- PARITY (rustle-VG vs StringTie) ---")
    parity_pct = (parity_matching_total / parity_ref_total * 100) if parity_ref_total else 0.0
    print(f"ST transcripts reproduced by VG      : {parity_matching_total} / {parity_ref_total} = {parity_pct:.2f}%")
    print(f"VG-vs-ST regressions (ST-only, =/c)  : {regressions_total}")
    print()
    print("--- NCBI RefSeq RECOVERY (unique ref transcripts matched) ---")
    print("STRICT tier (class_code in {=,c}):")
    print(f"  ST-matched : {len(GW_st_s)}")
    print(f"  VG-matched : {len(GW_vg_s)}")
    print(f"  both       : {len(both_s)}")
    print(f"  VG-only    : {len(vgonly_s)}   (NCBI tx rustle-VG recovers that StringTie misses)")
    print(f"  ST-only    : {len(stonly_s)}   (regressions: ST recovers, VG misses)")
    print("BROAD tier (class_code in {=,c,j}):")
    print(f"  ST-matched : {len(GW_st_b)}")
    print(f"  VG-matched : {len(GW_vg_b)}")
    print(f"  both       : {len(both_b)}")
    print(f"  VG-only    : {len(vgonly_b)}")
    print(f"  ST-only    : {len(stonly_b)}")
    print()
    yld_s = (len(GW_vg_s) / len(GW_st_s) * 100) if GW_st_s else 0.0
    yld_b = (len(GW_vg_b) / len(GW_st_b) * 100) if GW_st_b else 0.0
    print("--- ILLUSTRATIVE YIELD vs StringTie (can exceed 100%; illustrative, real Sn caps at 100%) ---")
    print(f"  STRICT yield = VG-matched / ST-matched = {len(GW_vg_s)}/{len(GW_st_s)} = {yld_s:.1f}%")
    print(f"  BROAD  yield = VG-matched / ST-matched = {len(GW_vg_b)}/{len(GW_st_b)} = {yld_b:.1f}%")

    # ===== VG-only gene list (strict tier) =====
    print()
    print("=" * 140)
    print("VG-ONLY NET-NEW (STRICT =/c) NCBI-CORROBORATED RECOVERIES, grouped by gene")
    print("=" * 140)
    gene_counts = defaultdict(int)
    for key in vgonly_s:
        g = gene_of_all.get(key, "?")
        g = re.sub(r"^gene-", "", g)
        gene_counts[g] += 1
    ordered = sorted(gene_counts.items(), key=lambda kv: (-kv[1], kv[0]))
    print(f"Total VG-only strict transcripts: {len(vgonly_s)} across {len(gene_counts)} genes")
    print(f"{'gene':<28}{'#tx':>5}")
    print("-" * 35)
    for g, n in ordered[:60]:
        print(f"{g:<28}{n:>5}")
    if len(ordered) > 60:
        print(f"... and {len(ordered)-60} more genes")

    # also dump full vg-only strict list to a file for reference
    with open(os.path.join(OUT, "vgonly_strict.tsv"), "w") as fh:
        fh.write("chrom_refid\tgene\n")
        for key in sorted(vgonly_s):
            fh.write(f"{key}\t{re.sub(r'^gene-','',gene_of_all.get(key,'?'))}\n")
    print(f"\n(full VG-only strict list written to {OUT}/vgonly_strict.tsv)")

if __name__ == "__main__":
    main()
