#!/usr/bin/env python3
"""sun_identifiability.py — per-COPY SUN (Singly Unique Nucleotide) identifiability catalog
over every validated GGO multi-copy gene family, refining psv_graph_genomewide.py's
family-level K into a per-copy 3-tier ladder.

SUN (Sudmant et al. 2010, Science 330(6004):641, doi:10.1126/science.1197005; 4.1M SUN
positions genome-wide; read-depth over SUNs => paralog-specific copy number / parCN, the
QuicK-mer2 / fastCN Eichler-lab marker): a POSITION where ONE copy's base is unique (a
private/singleton allele) among all copies of the family. A single read covering that one
position assigns DETERMINISTICALLY to that copy at the per-read gate (|N(r)|=1; THEORY.md §5
Lemma / Thm 4 / Thm 7).

Per-copy 3-tier ladder (a STRICT refinement of psv_graph's family K = # distinct hap-vectors):
  Tier 1  SUN-identifiable        : >=1 PSV column where the copy's allele is private
                                    -> single-read GATE-DETERMINISTIC: a read over the SUN
                                       column is never misassigned to another TRUE copy
                                       (|N(r)|=1). NOT cover-level immunity -- even a SUN-rich
                                       copy can be dissolved by an alternative minimum cover
                                       (THEORY.md §5 S3_cover); that is why production runs the
                                       per-read gate (Thm 4), not RECOVER/MCC.
  Tier 2  hap-vector-unique-only  : full hap-row is a singleton group but NO SUN
                                    -> combination-based; no single read pins it, needs a read
                                       co-observing >=2 PSVs (linkage). At the per-read gate a
                                       single-column read on this copy is always ambiguous.
  Tier 3  K-frontier/unresolvable : hap-row equals another copy's -> NM:i:0 collapse.
{Tier1} SUBSET-OF {Tier1 U Tier2} = {psv_graph singleton/resolvable}; {Tier3} = {psv_graph frontier}.

WITHIN-FAMILY, CONDITIONAL ON CATALOG COMPLETENESS: uniqueness here is "unique among the
family's enumerated validated copies", not genome-wide. Sudmant's SUNs are genome-wide-unique;
a called SUN is genome-wide-unique only if the family catalog is complete -- a missing paralog
could share the allele and demote the SUN. The columns come from an all-to-longest-copy (STAR)
projection, not a true multiple-sequence alignment; this only ever UNDER-counts SUNs (drops
gap/ambiguous columns).

COPY-ONLY: minimap2 asm20 self-alignment of the reference copies, alleles read from aligned
pairs. NO reads, NO SEDEF, NO read pileup (skips the 92-min break-test entirely). This is the
pure copy-sequence identifiability object (an upper bound on what RNA reads can observe).

MACHINE-CHECK (self-consistency, NOT independent evidence): recompute, from the raw hap dicts
via the equivalent PAIRWISE formulation (NOT re-reading has_sun), the SINGLE-POSITION
Strong-Separation witness set W(c) = {p : hap[c] != hap[c'] for all c'!=c} = INTERSECTION_{c'!=c}
D_{c,c'}. Since copy c always contributes exactly one occurrence to a column, `count==1` <=>
`pairwise-distinct`, so SUN_equals_witness is TRUE BY CONSTRUCTION -- a no-coding-bug / internal
consistency check, not independent corroboration. Genuine independence lives in bench/
sun_theory_check.py (abstract S1/S2, and S3_cover which tests the actual load-bearing claim) and
in re-deriving the tier counts from raw sequence.

Outputs:
  bench/sun_identifiability.tsv   per-copy: family, copy, chrom, start, end, gene, n_psv,
                                  n_sun, tier, hap_unique, group_size
  bench/sun_identifiability.json  genome-wide tier breakdown + K-refinement + Strong-Sep
                                  machine-check + named flagship examples.

Run: /home/juanfra/miniforge3/bin/python bench/sun_identifiability.py
"""
import collections
import json
import os
import sys
import traceback

import pysam

import psv_graph_genomewide as pg

TMP = "/home/juanfra/winloci_scratch/sun_ident"
ANNOT = "/home/juanfra/winloci_scratch/annot_intervals.tsv"
PSVGW = os.path.join(pg.HERE, "psv_graph_genomewide.json")


# ----------------------------------------------------------------------------- copy alleles
def copy_alleles_fast(bampath, ref):
    """{ref_pos: {copy_name: base}} from aligned pairs (O(aln length); no pileup).
    matches_only=True == pileup's query_position-not-None branch (indels excluded identically;
    verified byte-identical to psv_graph.column_alleles on fam3, 0/62244 columns differ)."""
    out = collections.defaultdict(dict)
    bam = pysam.AlignmentFile(bampath, "rb")
    for aln in bam.fetch(ref):
        if aln.is_unmapped or aln.is_supplementary or aln.is_secondary:
            continue
        qs = aln.query_sequence
        if qs is None:
            continue
        for qpos, rpos in aln.get_aligned_pairs(matches_only=True):
            out[rpos][aln.query_name] = qs[qpos]
    bam.close()
    return out


def copy_psvs(regions, genome, contig_len):
    """asm20-align the copies onto the longest one and read off the copy-only PSV columns.
    Returns (names, sorted_regions, psvs, n_aligned) or None (backbone < 200 bp).
    psvs = [(refpos, {copy_name: allele})] over columns where ALL copies are defined and
    there are >=2 distinct alleles (NO read gate)."""
    os.makedirs(TMP, exist_ok=True)
    regions = sorted(regions, key=lambda r: -(r[2] - r[1]))
    names = [f"copy{i}" for i in range(len(regions))]
    bc, bs, be, _ = regions[0]
    bs, be = max(0, bs), min(contig_len.get(bc, be), be)
    refseq = genome.fetch(bc, bs, be).upper()
    if len(refseq) < 200:
        return None
    open(f"{TMP}/ref.fa", "w").write(f">{names[0]}\n{refseq}\n")
    pg.sh(f"samtools faidx {TMP}/ref.fa")
    with open(f"{TMP}/copies.fa", "w") as o:
        for nm, (c, s, e, _) in list(zip(names, regions))[1:]:
            s2, e2 = max(0, s), min(contig_len.get(c, e), e)
            o.write(f">{nm}\n{genome.fetch(c, s2, e2).upper()}\n")
    pg.sh(f"minimap2 -ax asm20 {TMP}/ref.fa {TMP}/copies.fa 2>/dev/null "
          f"| samtools sort -o {TMP}/copies.bam - && samtools index {TMP}/copies.bam")
    copy_cols = copy_alleles_fast(f"{TMP}/copies.bam", names[0])
    aligned = {names[0]}
    for d in copy_cols.values():
        aligned.update(d)
    psvs = []
    for p in sorted(copy_cols):
        cc = copy_cols[p]
        if not all(nm in cc for nm in names[1:]):
            continue
        hap = {names[0]: refseq[p]}
        for nm in names[1:]:
            hap[nm] = cc[nm]
        if len(set(hap.values())) >= 2:
            psvs.append((p, hap))
    return names, regions, psvs, len(aligned)


# ----------------------------------------------------------------------------- tiers + check
def tier_family(names, psvs):
    """Per-copy SUN columns, hap-vectors, groups, tier. SUN(c,p) := hap[c] is a singleton
    among the values at p (its base appears exactly once)."""
    sun_cols = {nm: [] for nm in names}
    for (p, hap) in psvs:
        vals = list(hap.values())
        for nm in names:
            if vals.count(hap[nm]) == 1:
                sun_cols[nm].append(p)
    hap_rows = {nm: tuple(hap[nm] for _, hap in psvs) for nm in names}
    groups = collections.defaultdict(list)
    for nm in names:
        groups[hap_rows[nm]].append(nm)
    per_copy = {}
    for nm in names:
        singleton = len(groups[hap_rows[nm]]) == 1
        has_sun = len(sun_cols[nm]) > 0
        tier = 3 if not singleton else (1 if has_sun else 2)
        per_copy[nm] = dict(tier=tier, has_sun=has_sun, n_sun=len(sun_cols[nm]),
                            singleton=singleton, group_size=len(groups[hap_rows[nm]]))
    return per_copy, groups, sun_cols


def strong_sep_witness(names, psvs, sun_cols):
    """SELF-CONSISTENCY machine-check (NOT independent) of THEORY.md §5: for every copy c,
    recompute (from the raw hap dicts, NOT from sun_cols) the SINGLE-POSITION Strong-Separation
    witness columns
       W(c) = { p : for all c'!=c, hap[c] != hap[c'] at p }  =  p in INTERSECTION_{c'!=c} D_{c,c'}
    and verify, per copy:
      (i)  W(c) non-empty  <=>  c is SUN-identifiable (sun_cols[c] non-empty)   [SUN==witness]
      (ii) for each p in W(c): a read observing only p carrying (c)_p is consistent with
           EXACTLY c -> N(r)={c}, |N(r)|=1  (per-read GATE determinism / Thm 4/7 hypothesis)
      (iii) SUN(c) => hap-vector of c is unique (singleton group)
    NOTE: (i) is TRUE BY CONSTRUCTION -- W(c) uses `all(hap[o]!=hap[c])` while sun_cols uses
    `count(hap[c])==1`, and these are the SAME predicate (c contributes exactly one occurrence),
    so SUN_equals_witness catches coding bugs but is not independent evidence for SUN. The
    substantive, non-tautological determinism test is the abstract S2/S3_cover in
    sun_theory_check.py. Returns per-copy witness dict + violation counters."""
    # index columns for O(1) allele lookup
    cols = [(p, hap) for (p, hap) in psvs]
    witness = {nm: [] for nm in names}
    for (p, hap) in cols:
        for nm in names:
            a = hap[nm]
            if all(hap[o] != a for o in names if o != nm):
                witness[nm].append(p)
    # equivalence + |N(r)|=1 checks
    v_equiv = v_Nr = 0
    for nm in names:
        w_nonempty = len(witness[nm]) > 0
        s_nonempty = len(sun_cols[nm]) > 0
        if w_nonempty != s_nonempty:
            v_equiv += 1
        # |N(r)|=1 for a single-column read at each witness position
        wset = set(witness[nm])
        for (p, hap) in cols:
            if p not in wset:
                continue
            a = hap[nm]
            Nr = [o for o in names if hap[o] == a]
            if Nr != [nm]:
                v_Nr += 1
    return witness, dict(equiv_violations=v_equiv, Nr_violations=v_Nr)


# ----------------------------------------------------------------------------- annotation
def load_annot(path):
    ann = collections.defaultdict(list)
    if not os.path.exists(path):
        return ann
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            c, s, e, bt, g = parts[:5]
            ann[c].append((int(s), int(e), g, bt))
    for c in ann:
        ann[c].sort()
    return ann


def gene_for(ann, chrom, start, end):
    """Gene with max overlap of [start,end) on chrom (protein_coding preferred on ties)."""
    best_g, best_ov, best_pc = "", 0, False
    for s, e, g, bt in ann.get(chrom, ()):
        if e <= start:
            continue
        if s >= end:
            break
        ov = min(e, end) - max(s, start)
        pc = (bt == "protein_coding")
        if ov > best_ov or (ov == best_ov and pc and not best_pc):
            best_g, best_ov, best_pc = g, ov, pc
    return best_g or "."


# ----------------------------------------------------------------------------- main
def main():
    fams = collections.defaultdict(list)
    with open(pg.FAM_TSV) as f:
        next(f)
        for line in f:
            fi, lid, c, s, e, nr = line.rstrip("\n").split("\t")
            fams[int(fi)].append((lid, c, int(s), int(e), int(nr)))

    ann = load_annot(ANNOT)
    psvgw = {}
    if os.path.exists(PSVGW):
        for r in json.load(open(PSVGW)).get("families", []):
            psvgw[r["family"]] = r  # cls, K, group_sizes, psvs (read-gated)

    genome = pysam.FastaFile(pg.GENOME)
    contig_len = dict(zip(genome.references, genome.lengths))

    rows = []           # per-copy TSV rows
    fam_records = []     # per-family JSON records
    tier = collections.Counter()
    v_sun_singleton = 0  # SUN => singleton-group violations (should be 0)
    v_equiv = v_Nr = 0   # Strong-Sep witness machine-check violations (should be 0)
    n_copies_total = n_sun_copies = 0
    n_witness_copies = 0
    skipped = collections.Counter()
    processed = 0

    fam_ids = sorted(fams)
    for kk, fi in enumerate(fam_ids):
        if kk % 20 == 0:
            sys.stderr.write(f"  ... {kk}/{len(fam_ids)} (processed={processed})\n")
            sys.stderr.flush()
        regions0 = pg.dedup_copies(fams[fi])
        if len(regions0) < 2:
            skipped["single_region"] += 1
            continue
        try:
            res = copy_psvs(regions0, genome, contig_len)
            if res is None:
                skipped["backbone_short"] += 1
                continue
            names, regions, psvs, n_aligned = res
            per_copy, groups, sun_cols = tier_family(names, psvs)
            witness, wchk = strong_sep_witness(names, psvs, sun_cols)
            v_equiv += wchk["equiv_violations"]
            v_Nr += wchk["Nr_violations"]
            processed += 1

            genes = []
            for i, nm in enumerate(names):
                c, s, e, _ = regions[i]
                g = gene_for(ann, c, s, e)
                genes.append(g)
                pc = per_copy[nm]
                n_copies_total += 1
                tier[pc["tier"]] += 1
                if pc["has_sun"]:
                    n_sun_copies += 1
                    if not pc["singleton"]:
                        v_sun_singleton += 1
                if len(witness[nm]) > 0:
                    n_witness_copies += 1
                rows.append((fi, nm, c, s, e, g, len(psvs), pc["n_sun"], pc["tier"],
                             int(pc["singleton"]), pc["group_size"]))

            fam_gene = collections.Counter(g for g in genes if g != ".").most_common(1)
            fam_gene = fam_gene[0][0] if fam_gene else "."
            t1 = sum(1 for nm in names if per_copy[nm]["tier"] == 1)
            t2 = sum(1 for nm in names if per_copy[nm]["tier"] == 2)
            t3 = sum(1 for nm in names if per_copy[nm]["tier"] == 3)
            n_singleton = t1 + t2
            copyonly_K = len(groups)
            pg_row = psvgw.get(fi, {})
            fam_records.append(dict(
                family=fi, gene=fam_gene, n_copies=len(names), n_aligned=n_aligned,
                n_psv=len(psvs), copyonly_K=copyonly_K,
                n_singleton_groups=n_singleton, tier1=t1, tier2=t2, tier3=t3,
                fully_taggable=(t1 == len(names)),          # every copy SUN single-read taggable
                copyonly_fully_resolvable=(copyonly_K == len(names)),  # all copies distinct vectors
                has_recomb_vulnerable=(t2 > 0),
                psvgraph_cls=pg_row.get("cls"), psvgraph_K=pg_row.get("K"),
                psvgraph_psvs=pg_row.get("psvs"),
                per_copy={names[i]: dict(per_copy[names[i]], chrom=regions[i][0],
                                         start=regions[i][1], end=regions[i][2],
                                         gene=genes[i], n_witness=len(witness[names[i]]))
                          for i in range(len(names))}))
        except Exception as ex:
            skipped["error"] += 1
            sys.stderr.write(f"[fam {fi}] {ex}\n{traceback.format_exc()}\n")
    genome.close()

    n_unique_hap = tier[1] + tier[2]

    # ---- K-refinement: families where psv_graph's K counts non-single-read-taggable copies
    k_overcount_fams = [r for r in fam_records if r["has_recomb_vulnerable"]]
    fully_res_with_vuln = [r for r in fam_records
                           if r["copyonly_fully_resolvable"] and r["has_recomb_vulnerable"]]
    pg_fully_res_with_vuln = [r for r in fam_records
                              if r.get("psvgraph_cls") == "fully_resolvable"
                              and r["has_recomb_vulnerable"]]

    # ---- named flagship examples
    def name(r):
        return dict(family=r["family"], gene=r["gene"], n_copies=r["n_copies"],
                    n_psv=r["n_psv"], tier1=r["tier1"], tier2=r["tier2"], tier3=r["tier3"],
                    copyonly_K=r["copyonly_K"], psvgraph_cls=r.get("psvgraph_cls"))
    sun_rich = sorted([r for r in fam_records if r["tier1"] >= 2 and r["tier3"] == 0],
                      key=lambda r: (-r["tier1"], -r["n_psv"]))[:8]
    hap_unique_only = [name(r) for r in fam_records if r["tier2"] > 0]
    frontier_no_sun = sorted([r for r in fam_records
                              if r["tier1"] == 0 and r["tier3"] == r["n_copies"]],
                             key=lambda r: -r["n_copies"])[:8]
    # per-copy detail for the single Tier-2 copy (the load-bearing strict-subset witness)
    tier2_detail = []
    for r in fam_records:
        for nm, pc in r["per_copy"].items():
            if pc["tier"] == 2:
                tier2_detail.append(dict(family=r["family"], gene=r["gene"], copy=nm,
                                         chrom=pc["chrom"], start=pc["start"], end=pc["end"],
                                         copy_gene=pc["gene"], group_size=pc["group_size"],
                                         n_sun=pc["n_sun"], n_witness=pc["n_witness"]))

    summary = dict(
        families_processed=processed, skipped=dict(skipped),
        n_copies_total=n_copies_total,
        tier1_SUN_identifiable=tier[1],
        tier2_hapvector_unique_only=tier[2],
        tier3_frontier_collapsed=tier[3],
        pct_tier1=round(100 * tier[1] / max(n_copies_total, 1), 1),
        pct_tier2=round(100 * tier[2] / max(n_copies_total, 1), 1),
        pct_tier3=round(100 * tier[3] / max(n_copies_total, 1), 1),
        unique_hap_copies=n_unique_hap,
        SUN_strict_subset_of_unique_hap=(tier[1] < n_unique_hap),
        strict_subset_witnessed_by_n_tier2=tier[2],
        families_with_tier1=sum(1 for r in fam_records if r["tier1"] > 0),
        families_with_tier2_uniqueonly=sum(1 for r in fam_records if r["tier2"] > 0),
        families_all_tier3_frontier=sum(1 for r in fam_records
                                        if r["tier3"] == r["n_copies"]),
        families_fully_taggable_all_tier1=sum(1 for r in fam_records if r["fully_taggable"]),
        machine_check=dict(
            claim="every SUN-identifiable copy has a single-position Strong-Separation "
                  "witness p in INTERSECTION_{c'!=c} D_{c,c'} forcing per-read |N(r)|=1 "
                  "(THEORY.md §5). SELF-CONSISTENCY only: `count==1` and `pairwise-distinct` "
                  "are the same predicate, so SUN_equals_witness is true by construction "
                  "(no-coding-bug check, NOT independent). PER-READ gate determinism; this is "
                  "NOT cover-level copy-immunity -- see sun_theory_check.py S3_cover.",
            is_self_consistency_not_independent=True,
            n_sun_identifiable_copies=n_sun_copies,
            n_copies_with_single_position_strongsep_witness=n_witness_copies,
            SUN_equals_witness=(n_sun_copies == n_witness_copies),
            SUN_witness_equivalence_violations=v_equiv,
            Nr_equals_1_violations=v_Nr,
            SUN_implies_unique_hap_violations=v_sun_singleton,
            all_green=(v_equiv == 0 and v_Nr == 0 and v_sun_singleton == 0
                       and n_sun_copies == n_witness_copies),
        ),
        k_refinement=dict(
            note="psv_graph's family K = # distinct copy hap-vectors counts every singleton "
                 "group as a 'resolvable' copy; Tier 2 copies are singleton groups (add to K) "
                 "but are NOT single-read taggable -- no SUN, so no single read pins them; they "
                 "need a read co-observing >=2 PSVs (linkage), and at the per-read gate a "
                 "single-column read on such a copy is ambiguous. K over-counts the single-read "
                 "gate-taggable copies by the number of Tier-2 copies.",
            distinct_hap_vectors_are_resolvable=n_unique_hap,
            single_read_taggable_tier1=tier[1],
            over_count_recomb_vulnerable_tier2=tier[2],
            n_families_K_counts_non_single_read_taggable=len(k_overcount_fams),
            n_families_copyonly_fully_resolvable_with_vulnerable_copy=len(fully_res_with_vuln),
            n_families_psvgraph_fully_resolvable_with_vulnerable_copy=len(pg_fully_res_with_vuln),
            families_K_overcounts=[name(r) for r in k_overcount_fams],
            tier2_copy_detail=tier2_detail,
        ),
        examples=dict(
            SUN_rich=[name(r) for r in sun_rich],
            hap_unique_only=hap_unique_only,
            K_frontier_no_SUN=[name(r) for r in frontier_no_sun],
        ),
    )

    out = dict(summary=summary, families=fam_records)
    json.dump(out, open(os.path.join(pg.HERE, "sun_identifiability.json"), "w"),
              indent=2, sort_keys=True)
    with open(os.path.join(pg.HERE, "sun_identifiability.tsv"), "w") as o:
        o.write("family\tcopy\tchrom\tstart\tend\tgene\tn_psv\tn_sun\ttier\thap_unique\tgroup_size\n")
        for row in sorted(rows):
            o.write("\t".join(str(x) for x in row) + "\n")

    print(json.dumps(summary, indent=2))
    print(f"\n[+] wrote bench/sun_identifiability.tsv ({len(rows)} copies) "
          f"+ bench/sun_identifiability.json ({processed} families)")


if __name__ == "__main__":
    main()
