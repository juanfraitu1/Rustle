#!/usr/bin/env python3
"""vg_read_path_recombination.py -- VG-NATIVE read-path RECOMBINATION detector (build stage).

THE OBJECT (the K-frontier recombination obstruction, made concrete per-read):
  In the O2 copy-assignment VARIATION GRAPH (bench/o2_vg_visualization.py), a multi-copy family is
  ONE graph: COPIES are PATHS, PSVs are BUBBLES. A read threads the bubbles -> a per-bubble
  allele-vector. A read whose allele-vector is consistent with NO single copy-path but IS a clean
  CONCATENATION of two copies' paths -- copy A over a prefix of bubbles, copy B over a suffix -- is a
  RECOMBINANT read. It is the concrete form of the copy_assignment_theory.md / THEORY.md K-frontier
  RECOMBINATION OBSTRUCTION: separate+link resolves copies at K=2 but a cross-copy recombinant read at
  K>=3 SPOOFS a copy combination and belongs to NO single copy -> it must ABSTAIN, not be force-assigned.

  Two kinds (bubble-BOUNDED, so the tract endpoints are exact graph objects):
    CROSSOVER          : one clean switch     A A A | B B B          (1 switch; the read crosses A->B)
    GENE-CONVERSION    : a local switched tract A A | B B | A A       (2 switches bounding a Y tract;
                         the first and last switched bubble = the tract endpoints)

WHAT IT REUSES (nothing re-derived):
  * bench/o2_vg_visualization.materialize_family -- the materialized family VG (PSV bubbles with
    per-copy allele-vectors, per-read obs={bubble_col -> base}, the 154 co-located families,
    backbone genomic frame). Same machinery as psv_graph_genomewide.py + the copy_assign gate.
  * the shipped RT/template-switch discriminator (src/rustle/vg_family/mosaic.rs::classify_event +
    genome.rs::breakpoint_microhomology). Ported byte-faithfully here (Python): the SAME three-leg
    truth table (recurrence via `confirmed`, microhomology direct-repeat leg k in [6,12] with the
    >=3-distinct-base low-complexity guard, DNA positive-only veto) and the SAME shipped constants
    (family_min_supporting_reads=3, breakpoint_tol=50, MH_KMIN=6, MH_KMAX=12).

METHOD (per co-located family with >=2 copies AND >=2 PSV bubbles):
  1. copy_vecs[copy][col] = bubble's per-copy allele. read obs[col] = read's base at that bubble.
  2. for each read spanning >= MIN_INFO discriminating bubbles that NO single copy explains
     (min single-copy mismatch >= MIN_SINGLE_MM): search copy PAIRS (A,B) for a PURE ordered split
     (every pair-discriminating bubble the read spans matches A's or B's allele, <= IMPURITY_MAX
     exceptions) with exactly 1 switch (crossover, both arms >= MIN_ARM) or 2 switches (conversion,
     both flanks >= MIN_FLANK, tract >= MIN_TRACT). Bubbles bound the tract.
  3. GATE the artifacts. Cluster recombinant reads within a family by (unordered pair, switch type,
     breakpoint midpoint within breakpoint_tol); a cluster is RECURRENT (`confirmed`) iff
     >= 3 distinct molecules with dispersion <= 50. Compute the microhomology direct-repeat leg at
     each switch bracket on the backbone. classify_event(confirmed, microhomology, dna):
        microhomology & !dna_present  -> RtSwitchArtifact  (template signature -> VETOED)
        confirmed & !mh & !dna_absent -> GeneConversion    (recurrent, no signature -> REAL)
        !confirmed & !mh              -> ChimeraSuspect     (sporadic one-off -> likely error/artifact)
        else                          -> Ambiguous
     DNA leg is positive-only (README: catalog "absent" is unreliable); default unchecked (None),
     exactly the shipped copy_assign/denovo path. A DNA-catalog probe is attempted if a matching
     DNFAM mosaic can be located, else None.
  4. connect to the gate: report, for each recombinant read, what the O2 gate CURRENTLY does
     (assigned = force-assigned LEAK vs ambiguous/tied = already abstains). These reads SHOULD abstain.

OUTPUT: bench/vg_read_path_recombination.tsv  (one row per recombinant read) + a JSON summary +
stdout headline (crossover vs conversion pre/post RT-switch veto, families, tract-length distribution).

Deterministic (PYTHONHASHSEED=0; sorted iteration; no RNG).
Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/vg_read_path_recombination.py [--limit N] [fam_ids...]
"""
import collections
import csv
import itertools
import json
import os
import statistics
import sys
import traceback

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import psv_graph_genomewide as pg          # noqa: E402  (FAM_TSV, GENOME, dedup_copies)
import o2_vg_visualization as o2vg         # noqa: E402  (materialize_family, load_families)
# SINGLE SOURCE of the recombinant-detection logic + its thresholds (moved out of this file so the
# shipped abstain leg and this diagnostic detector run byte-identical code -- NOT re-derived).
import recombinant_abstain as ra           # noqa: E402
from recombinant_abstain import (          # noqa: E402  re-export so existing references resolve unchanged
    MIN_INFO, MIN_ARM, MIN_FLANK, MIN_TRACT, MIN_SINGLE_MM, IMPURITY_MAX,
    ALLELE_SUPPORT_ABS, ALLELE_SUPPORT_FRAC,
    full_pattern_switches, detect_read_recombination, build_supported)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_TSV = os.path.join(HERE, "vg_read_path_recombination.tsv")
OUT_JSON = os.path.join(HERE, "vg_read_path_recombination.json")
DNA_CATALOG = "/home/juanfra/winloci_scratch/dna_catalog"

# ---- shipped mosaic/RT-switch constants (mosaic.rs defaults; genome.rs MH k-range) ----
FAMILY_MIN_SUPPORTING_READS = 3   # recurrence -> `confirmed`
BREAKPOINT_TOL = 50               # cluster merge tolerance (bp)
MAX_BREAKPOINT_DISPERSION = 50    # confirmed cluster dispersion cap (bp)
MH_KMIN, MH_KMAX = 6, 12          # microhomology direct-repeat k-range


# =========================================================================== RT-switch microhomology
# Byte-faithful port of genome.rs::{is_rt_switch, breakpoint_microhomology, is_low_complexity_window}.
def _is_low_complexity_window(w):
    """< 3 distinct ACGT bases (homopolymer / dinucleotide repeat) -> uninformative direct repeat."""
    seen = set()
    for b in w.upper():
        if b in "ACGT":
            seen.add(b)
    return len(seen) < 3


def _is_rt_switch(genome, chrom, donor, acceptor, k):
    """The `k` bp ending at `donor` == the `k` bp ending at `acceptor` (N-free direct repeat)."""
    if donor < k or acceptor < k:
        return False
    try:
        up = genome.fetch(chrom, donor - k, donor).upper()
        dn = genome.fetch(chrom, acceptor - k, acceptor).upper()
    except Exception:
        return False
    return bool(up) and up == dn and "N" not in up


def breakpoint_microhomology(genome, chrom, left, right, k_min=MH_KMIN, k_max=MH_KMAX):
    """True if ANY k in [k_min,k_max] shows a NON-low-complexity direct repeat across the switch
    bracket (left,right) on the backbone -- the RT/template-switch signature."""
    for k in range(k_min, k_max + 1):
        if left < k or not _is_rt_switch(genome, chrom, left, right, k):
            continue
        try:
            w = genome.fetch(chrom, left - k, left).upper()
        except Exception:
            continue
        if not _is_low_complexity_window(w):
            return True
    return False


def classify_event(confirmed, microhomology, dna_supported):
    """Port of mosaic.rs::classify_event (evaluation order load-bearing; artifact rule fires first)."""
    mh = microhomology is True
    dna_present = dna_supported is True
    dna_absent = dna_supported is False
    if mh and not dna_present:
        return "RtSwitchArtifact"
    if confirmed and not mh and not dna_absent:
        return "GeneConversion"
    if not confirmed and not mh:
        return "ChimeraSuspect"
    return "Ambiguous"


# =========================================================================== per-read detection
# full_pattern_switches + detect_read_recombination now live in recombinant_abstain.py (single source,
# imported above) so the shipped abstain leg and this diagnostic detector run byte-identical code.


# =========================================================================== DNA positive-only probe
_DNA_CACHE = {}


def dna_mosaic_positions(fam_id):
    """Positive-only DNA corroboration (README: catalog 'absent' is unreliable -> never Some(false)).
    A DNA gene-conversion signature = a catalog copy whose allele-vector is a MOSAIC of two others.
    We return the set of genomic breakpoint positions of such DNA mosaics for DNFAM<fam_id> if present,
    else None (UNCHECKED). Coordinate/identity fragility (README B/C/D) means this is best-effort."""
    if fam_id in _DNA_CACHE:
        return _DNA_CACHE[fam_id]
    path = os.path.join(DNA_CATALOG, f"DNFAM{fam_id}.json")
    res = None
    if os.path.exists(path):
        try:
            cat = json.load(open(path))
            matrix = cat.get("matrix", cat if isinstance(cat, dict) else {})
            # matrix[label][pos] = base (ref0-centric sparse). A mosaic label = its allele vector equals
            # a prefix of one label + suffix of another at the differing positions. Best-effort: collect
            # all catalog PSV genomic positions as candidate breakpoint anchors (positive corroboration
            # = a recombinant switch coincides with a region carrying DNA PSVs / a DNA mosaic).
            labels = {k: v for k, v in matrix.items() if isinstance(v, dict)}
            ref0 = cat.get("ref0_start", cat.get("start", 0)) or 0
            positions = set()
            for lab, row in labels.items():
                for p in row:
                    try:
                        positions.add(int(ref0) + int(p))
                    except Exception:
                        pass
            res = positions if positions else None
        except Exception:
            res = None
    _DNA_CACHE[fam_id] = res
    return res


# =========================================================================== per-family driver
def process_family(vg, genome):
    """Detect + gate recombinant reads for one materialized family VG. Returns (rows, fam_summary)."""
    names = [c["name"] for c in vg["copies"]]
    if len(names) < 2 or vg["n_bubbles"] < 2:
        return [], None
    gene_of = {c["name"]: c["gene"] for c in vg["copies"]}
    loc_of = {c["name"]: (c["chrom"], c["start"], c["end"]) for c in vg["copies"]}
    per_copy = vg.get("per_copy", {})      # gate-assigned read count per copy (for recombinant-fraction)
    bubbles = vg["bubbles"]
    chrom = vg["backbone_chrom"]
    bstart = vg["backbone_start"]
    gpos = [bstart + b["pos"] for b in bubbles]        # genomic pos of each bubble col
    copy_vecs = {nm: {j: bubbles[j]["hap"].get(nm) for j in range(len(bubbles))} for nm in names}
    status_of = {r["name"]: r["status"] for r in vg["reads"]}

    # ---- read-supported alleles per bubble (drops copy-vector miscalls + singleton errors) ----
    supported = build_supported(vg["reads"], len(bubbles))    # factored -> recombinant_abstain.py

    # ---- detect recombinant reads ----
    recs = []
    for r in vg["reads"]:
        obs = {int(k): v for k, v in r["obs"].items() if v in "ACGT"}
        if len(obs) < MIN_INFO:
            continue
        cand = detect_read_recombination(obs, supported, copy_vecs, names)
        if cand is None:
            continue
        seq = cand["seq"]
        # switch brackets: for each boundary index i in bnd_idx, the switch is between seq[i-1] and seq[i]
        brackets = []
        for i in cand["bnd_idx"]:
            cl, cr = seq[i - 1][0], seq[i][0]          # bubble cols flanking the switch
            gl, gr = gpos[cl], gpos[cr]
            brackets.append((cl, cr, gl, gr))
        # microhomology: any switch bracket carries the direct-repeat signature
        mh = any(breakpoint_microhomology(genome, chrom, gl, gr) for _, _, gl, gr in brackets)
        # signature midpoint for clustering (crossover: switch mid; conversion: tract mid)
        mid = sum((gl + gr) / 2.0 for _, _, gl, gr in brackets) / len(brackets)
        if cand["type"] == "conversion":
            t0, t1 = cand["tract_idx"]
            tract_cols = [seq[k][0] for k in range(t0, t1 + 1)]
            tract_bubbles = len(tract_cols)
            tract_g0, tract_g1 = gpos[seq[t0][0]], gpos[seq[t1][0]]
        else:
            tract_cols = [brackets[0][0], brackets[0][1]]
            tract_bubbles = 0                          # crossover has no interior tract
            tract_g0, tract_g1 = brackets[0][2], brackets[0][3]
        recs.append(dict(
            read=r["name"], A=cand["A"], B=cand["B"], type=cand["type"], switches=cand["switches"],
            n_info=len(seq), imp=cand["imp"], mid=mid, mh=mh, brackets=brackets,
            full_switches=cand.get("full_switches", cand["switches"]),
            n_full_ab=cand.get("n_full_ab", len(seq)), clean_full=cand.get("clean_full", 1),
            tract_bubbles=tract_bubbles, tract_cols=tract_cols,
            tract_g0=min(tract_g0, tract_g1), tract_g1=max(tract_g0, tract_g1),
            gate_status=status_of.get(r["name"], "?")))
    if not recs:
        return [], dict(family=vg["family"], gene=vg.get("gene", "?"), n_copies=len(names),
                        K=vg["K"], n_bubbles=vg["n_bubbles"], n_recombinant=0)

    # ---- cluster within family by (unordered pair, type, midpoint within tol); recurrence -> confirmed
    dna_pos = dna_mosaic_positions(vg["family"])
    clusters = []
    keyfn = lambda x: (tuple(sorted((x["A"], x["B"]))), x["type"], x["mid"])
    for rec in sorted(recs, key=keyfn):
        pair = tuple(sorted((rec["A"], rec["B"])))
        placed = False
        for cl in clusters:
            if cl["pair"] == pair and cl["type"] == rec["type"] and \
               abs(rec["mid"] - cl["mids"][-1]) <= BREAKPOINT_TOL:
                cl["members"].append(rec)
                cl["mids"].append(rec["mid"])
                placed = True
                break
        if not placed:
            clusters.append(dict(pair=pair, type=rec["type"], members=[rec], mids=[rec["mid"]]))

    rows = []
    fam_cls = collections.Counter()
    for ci, cl in enumerate(clusters):
        mids = sorted(cl["mids"])
        n_mol = len({m["read"] for m in cl["members"]})
        dispersion = mids[-1] - mids[0]
        confirmed = n_mol >= FAMILY_MIN_SUPPORTING_READS and dispersion <= MAX_BREAKPOINT_DISPERSION
        mh = any(m["mh"] for m in cl["members"])       # breakpoint property (constant within cluster)
        # DNA positive-only: does any switch bracket coincide (+-BREAKPOINT_TOL) with a DNA mosaic PSV?
        dna = None
        if dna_pos:
            hit = any(any(abs(gl - d) <= BREAKPOINT_TOL or abs(gr - d) <= BREAKPOINT_TOL
                          for _, _, gl, gr in m["brackets"] for d in dna_pos)
                      for m in cl["members"])
            dna = True if hit else None                 # positive-only (never Some(false))
        klass = classify_event(confirmed, mh, dna)
        # RECOMBINANT FRACTION (honesty leg, beyond classify_event): of the reads that look like the
        # (majority) flank copy, what fraction switch? A localized gene-conversion tract affects a
        # MINORITY (recFrac low); a near-1 fraction = essentially ALL of a copy's reads "recombine" =
        # a SYSTEMATIC near-identical-copy split / an unrepresented recombinant copy, NOT a conversion.
        pair = cl["pair"]
        flank_pool = max(per_copy.get(pair[0], 0), per_copy.get(pair[1], 0)) + n_mol
        rec_frac = round(n_mol / max(flank_pool, 1), 3)
        # honest structural call layered on the shipped classification:
        if not confirmed:
            structural = "sporadic_chimera"           # per-molecule one-off (RT/PCR chimera or error)
        elif cl["type"] == "conversion" and rec_frac < 0.5:
            structural = "localized_tract"            # minority, bounded tract = credible gene conversion
        else:
            structural = "systematic_copy_split"      # majority / whole-molecule crossover = copy ambiguity
        cid = f"F{vg['family']}_C{ci}"
        for m in cl["members"]:
            fam_cls[klass] += 1
            rows.append(dict(
                family=vg["family"], gene=vg.get("gene", "?"),
                copyA=m["A"], geneA=gene_of.get(m["A"], "?"),
                copyB=m["B"], geneB=gene_of.get(m["B"], "?"),
                chromA="%s:%d-%d" % loc_of.get(m["A"], ("?", 0, 0)),
                chromB="%s:%d-%d" % loc_of.get(m["B"], ("?", 0, 0)),
                read=m["read"], switch_type=m["type"], n_switches=m["switches"],
                full_switches=m["full_switches"], n_full_ab=m["n_full_ab"], clean_full=m["clean_full"],
                n_info_bubbles=m["n_info"], impurity=m["imp"],
                tract_bubbles=m["tract_bubbles"],
                tract_chrom=chrom, tract_start=m["tract_g0"], tract_end=m["tract_g1"],
                tract_bp=m["tract_g1"] - m["tract_g0"],
                switch_brackets=";".join(f"{cl_}-{cr_}" for cl_, cr_, _, _ in m["brackets"]),
                microhomology=int(m["mh"]),
                cluster_id=cid, cluster_size=n_mol, recurrent=int(confirmed),
                recomb_fraction=rec_frac, structural_call=structural,
                dna_support=("yes" if dna is True else "unchecked"),
                classification=klass,
                gate_status=m["gate_status"],
                should_abstain=1,
                gate_leak=int(m["gate_status"] == "assigned")))
    fam_summary = dict(
        family=vg["family"], gene=vg.get("gene", "?"), n_copies=len(names), K=vg["K"],
        n_bubbles=vg["n_bubbles"], n_recombinant=len(recs), n_clusters=len(clusters),
        n_crossover=sum(1 for r in recs if r["type"] == "crossover"),
        n_conversion=sum(1 for r in recs if r["type"] == "conversion"),
        classes=dict(fam_cls),
        n_gate_leak=sum(1 for r in rows if r["gate_leak"]),
        chroms=vg["chroms"])
    return rows, fam_summary


# =========================================================================== EVAL STAGE
# Consumes the build outputs (OUT_TSV per-read rows + OUT_JSON summary) and CROSS-REFERENCES the
# SUN identifiability ladder (bench/sun_identifiability.json: per-copy Tier 1 SUN-identifiable /
# Tier 2 hap-unique-only [recombination-vulnerable] / Tier 3 collapse). Answers the four questions:
#   (1) how many REAL crossover/conversion reads survive the RT-switch veto (recurrent, MH-checked,
#       DNA where possible), by NAMED family;
#   (2) the K-frontier connection: are recombinant reads enriched in K>=3 families, and what SUN
#       tiers do the surviving reads bridge (the read-level recombination obstruction);
#   (3) the gate connection: are these reads force-assigned to one copy (a mis-assignment LEAK) or
#       correctly abstained -- quantified overall and per structural class;
#   (4) whether bubble-bounding + VG-native detection adds over the shipped mosaic_discriminator.
# Reads its inputs from disk so the standalone `eval` mode is byte-identical to the inline call.
# Deterministic (sorted iteration, no RNG).
EVAL_INT_COLS = ("family", "n_switches", "full_switches", "n_full_ab", "clean_full",
                 "n_info_bubbles", "impurity", "tract_bubbles",
                 "tract_start", "tract_end", "tract_bp", "microhomology", "cluster_size",
                 "recurrent", "should_abstain", "gate_leak")
EVAL_FLOAT_COLS = ("recomb_fraction",)


def _load_rows():
    rows = []
    with open(OUT_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for k in EVAL_INT_COLS:
                r[k] = int(r[k])
            for k in EVAL_FLOAT_COLS:
                r[k] = float(r[k])
            rows.append(r)
    return rows


def _load_sun():
    """tier_of[(fam,copy)] ; fam_gene[fam] ; fam_meta[fam] = SUN K-tier profile per family."""
    path = os.path.join(HERE, "sun_identifiability.json")
    tier_of, fam_gene, fam_meta = {}, {}, {}
    if not os.path.exists(path):
        return tier_of, fam_gene, fam_meta
    for f in json.load(open(path)).get("families", []):
        fid = f["family"]
        fam_gene[fid] = f.get("gene") or "?"
        fam_meta[fid] = dict(sun_K=f.get("copyonly_K"), tier1=f.get("tier1"), tier2=f.get("tier2"),
                             tier3=f.get("tier3"), recomb_vulnerable=f.get("has_recomb_vulnerable"),
                             fully_taggable=f.get("fully_taggable"))
        for cp, info in f.get("per_copy", {}).items():
            tier_of[(fid, cp)] = info.get("tier")
    return tier_of, fam_gene, fam_meta


def _fisher(a, b, c, d):
    """2x2 Fisher exact (scipy if present, else None); returns (odds, p)."""
    try:
        from scipy.stats import fisher_exact
        return fisher_exact([[a, b], [c, d]])
    except Exception:
        return (None, None)


def evaluate():
    rows = _load_rows()
    build = json.load(open(OUT_JSON))
    fam_summaries = build.get("families", [])
    tier_of, fam_gene, fam_meta = _load_sun()
    Kof = {fs["family"]: fs["K"] for fs in fam_summaries}

    def gname(fid):
        return fam_gene.get(fid, "?") or "?"

    def nmol(subset):
        """DISTINCT MOLECULES: read names are globally-unique CCS ids, so a read appearing in two
        families (duplicate-locus materialization, e.g. fam10==fam132 LOC115930576) is ONE molecule."""
        return len({r["read"] for r in subset})

    # duplicate-locus family map: families whose recombinant-read-name SETS are identical are the SAME
    # locus materialized twice (fam10/fam132). canon[fam] = the representative (min) family id.
    fam_reads = collections.defaultdict(set)
    for r in rows:
        fam_reads[int(r["family"])].add(r["read"])
    canon = {}
    for f in sorted(fam_reads):
        canon.setdefault(f, f)
        for g in sorted(fam_reads):
            if g > f and g not in canon and fam_reads[f] and fam_reads[f] == fam_reads[g]:
                canon[g] = f

    # ---------- (1) SURVIVORS: RT-switch veto (INERT) + structural + FULL-PATTERN credibility ----------
    by_class = collections.Counter(r["classification"] for r in rows)
    by_struct = collections.Counter(r["structural_call"] for r in rows)
    mh_fired = sum(1 for r in rows if r["microhomology"] == 1)          # RT-switch signature hits
    dna_yes = sum(1 for r in rows if r["dna_support"] == "yes")
    loc = [r for r in rows if r["structural_call"] == "localized_tract"]
    # FULL-PATTERN CREDIBILITY (adversarial review). The detected "2-switch tract" is fit on the subset
    # of columns where BOTH copies' alleles are read-supported; over ALL read-supported discriminating
    # columns the same read can carry many more switches (full_switches). A tract is CLEAN only if
    # clean_full==1 (no hidden switch); a clean_full==0 read is a noisy A/B mosaic whose clean tract is a
    # column-subsetting artifact (the dropped columns are asm20 copy-vector miscalls -- ~49% of fam13's
    # discriminating columns carry a copy allele with ~0 read support -- and/or read noise). Credible-real
    # is therefore ONLY the clean_full localized tracts.
    loc_clean = [r for r in loc if r["clean_full"] == 1]
    loc_noisy = [r for r in loc if r["clean_full"] == 0]
    credible = loc_clean                                    # honest credible-real headline set

    surv = {}
    for r in credible:
        fid = canon.get(int(r["family"]), int(r["family"]))   # collapse duplicate-locus families
        d = surv.setdefault(fid, dict(gene=gname(fid), K=Kof.get(fid), reads=set(),
                                      tract_bubbles=collections.Counter(), tract_bp=[], rec_frac=[],
                                      full_sw=[]))
        d["reads"].add(r["read"])
        d["tract_bubbles"][r["tract_bubbles"]] += 1
        d["tract_bp"].append(r["tract_bp"])
        d["rec_frac"].append(r["recomb_fraction"])
        d["full_sw"].append(r["full_switches"])
    survivors = []
    for fid in sorted(surv, key=lambda x: (-len(surv[x]["reads"]), x)):
        d = surv[fid]
        survivors.append(dict(
            locus=fid, gene=d["gene"], K=d["K"], reads=len(d["reads"]),
            tract_bubbles=dict(sorted(d["tract_bubbles"].items())),
            tract_bp_median=int(statistics.median(d["tract_bp"])),
            recomb_fraction_median=round(statistics.median(d["rec_frac"]), 3),
            full_switch_median=int(statistics.median(d["full_sw"]))))
    part1 = dict(
        post_veto_classification=dict(by_class),
        total_recombinant_rows=len(rows), total_recombinant_molecules=nmol(rows),
        rt_switch_microhomology_hits=mh_fired,
        rt_switch_note=("The ported RT-switch microhomology veto (mosaic.rs::classify_event direct-repeat "
                  "leg) fires %d/%d -- INERT for adjacent-PSV paralog switches (it targets DISTAL "
                  "splice-like donor/acceptor repeats). CONSEQUENCE: RT-template-switch is NOT actually "
                  "excluded from the survivors; the only working filters are RECURRENCE + recomb-fraction "
                  "+ full-pattern cleanliness, and a recurrent RT/PCR chimera at a homology hotspot (where "
                  "these live) would pass them. 'Credible' means recurrent-clean-tract-on-a-known-locus, "
                  "NOT artifact-excluded." % (mh_fired, len(rows))),
        dna_supported_reads=dna_yes,
        dna_note=("DNA UNCHECKED for all reads: the DNFAM catalog uses an independent id-space "
                  "(DNFAM0..1000+ are DNA families, NOT RNA family ids 0..171); the id-based probe "
                  "and coordinate join find no reliable match, so no read is DNA-confirmed. Neither the "
                  "RT-switch veto (0 hits) nor DNA (0 hits) is delivering the artifact separation -- "
                  "recurrence + structure do, and recurrence alone cannot distinguish real gene conversion "
                  "from a systematically mis-mapped or unrepresented recombinant paralog copy."),
        allele_vs_conversion_note=("Surviving tract columns are ordinary heterozygous sites (e.g. GSTM2 "
                  "col C:259/T:18) that RNA alone CANNOT distinguish from a heterozygous SNP at one copy's "
                  "locus without DNA parCN. So 'gene-conversion tract' is a HYPOTHESIS (recurrent, "
                  "bubble-bounded, on known loci), not a DNA-verified finding."),
        structural_split=dict(by_struct),
        localized_tract_rows=len(loc), localized_tract_molecules=nmol(loc),
        localized_clean_full_rows=len(loc_clean), localized_clean_full_molecules=nmol(loc_clean),
        localized_noisy_mosaic_rows=len(loc_noisy), localized_noisy_mosaic_molecules=nmol(loc_noisy),
        credible_real_conversion_rows=len(credible),
        credible_real_conversion_molecules=nmol(credible),
        credible_real_crossover_reads=sum(1 for r in credible if r["switch_type"] == "crossover"),
        credible_distinct_loci=len(surv),
        survivors_by_locus=survivors,
        recomb_fraction_median_credible=round(statistics.median(
            [r["recomb_fraction"] for r in credible] or [0]), 3),
        recomb_fraction_median_systematic=round(statistics.median(
            [r["recomb_fraction"] for r in rows if r["structural_call"] == "systematic_copy_split"] or [0]), 3),
        full_switch_median_credible=int(statistics.median([r["full_switches"] for r in credible] or [0])),
        full_switch_median_noisy=int(statistics.median([r["full_switches"] for r in loc_noisy] or [0])),
        headline=("RT-switch microhomology veto fires %d/%d (INERT). Of %d raw recombinant rows "
                  "(%d distinct molecules), the honest structural filter keeps %d localized-tract rows and "
                  "demotes %d as systematic copy-ambiguity + %d as sporadic chimera. Of those localized "
                  "tracts only %d rows / %d molecules pass the FULL-PATTERN cleanliness test (clean_full=1) "
                  "across %d distinct loci; the other %d are noisy A/B mosaics whose clean 2-switch tract "
                  "is a column-subsetting artifact of unreliable asm20 copy-vectors."
                  % (mh_fired, len(rows), len(rows), nmol(rows), len(loc),
                     by_struct.get("systematic_copy_split", 0), by_struct.get("sporadic_chimera", 0),
                     len(credible), nmol(credible), len(surv), len(loc_noisy))))

    # ---------- (2) K-FRONTIER connection ----------
    kb = {"K>=3": [0, 0], "K==2": [0, 0], "K<2": [0, 0]}     # [carry_recomb, no_recomb]
    for fs in fam_summaries:
        k = fs["K"]
        b = "K>=3" if k >= 3 else ("K==2" if k == 2 else "K<2")
        kb[b][0 if fs["n_recombinant"] > 0 else 1] += 1
    odds, p = _fisher(kb["K>=3"][0], kb["K>=3"][1], kb["K==2"][0], kb["K==2"][1])
    read_mass = collections.Counter("K>=3" if Kof.get(int(r["family"]), 0) >= 3 else "K==2" for r in rows)
    loc_mass = collections.Counter("K>=3" if Kof.get(int(r["family"]), 0) >= 3 else "K==2" for r in credible)
    n_kge3_fam = sum(1 for fs in fam_summaries if fs["K"] >= 3)
    # SUN-tier bridges of the CREDIBLE (clean_full localized) reads: which two tiers each read joins
    tier_pairs = collections.Counter()
    for r in credible:
        ta = tier_of.get((int(r["family"]), r["copyA"]), "?")
        tb = tier_of.get((int(r["family"]), r["copyB"]), "?")
        tier_pairs[tuple(sorted((str(ta), str(tb))))] += 1
    # DETECTION-POWER CONFOUND (adversarial review): more copies/bubbles give the pair-fitting detector
    # more freedom, so part of the K-enrichment is MECHANICAL, not biological. Correlate per-family
    # n_recombinant against n_bubbles and against C(K,2)=#copy-pairs.
    def _pearson(xs, ys):
        n = len(xs)
        if n < 3:
            return None
        mx, my = sum(xs) / n, sum(ys) / n
        sx = sum((x - mx) ** 2 for x in xs) ** 0.5
        sy = sum((y - my) ** 2 for y in ys) ** 0.5
        if sx == 0 or sy == 0:
            return None
        return round(sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (sx * sy), 3)
    xb = [fs["n_bubbles"] for fs in fam_summaries]
    xk = [fs["K"] * (fs["K"] - 1) / 2 for fs in fam_summaries]
    yn = [fs["n_recombinant"] for fs in fam_summaries]
    r_bubbles = _pearson(xb, yn)
    r_pairs = _pearson(xk, yn)
    # the copy-level recomb-vulnerable (SUN Tier-2) families vs which actually carry recombinants
    recomb_fams = {int(r["family"]) for r in rows}
    vuln_fams = sorted(fid for fid, m in fam_meta.items() if m.get("recomb_vulnerable"))
    vuln_with_reads = sorted(f for f in vuln_fams if f in recomb_fams)
    part2 = dict(
        families_carrying_recombinants_by_K={k: dict(carry=v[0], no_recomb=v[1]) for k, v in kb.items()},
        enrichment_fisher=dict(odds_ratio=(round(odds, 3) if odds else None),
                               p_value=(float("%.4g" % p) if p else None),
                               statement="K>=3 %d/%d fams carry recombinants vs K==2 %d/%d"
                               % (kb["K>=3"][0], sum(kb["K>=3"]), kb["K==2"][0], sum(kb["K==2"]))),
        detection_power_confound=dict(
            pearson_nbubbles_vs_nrecombinant=r_bubbles,
            pearson_copypairs_vs_nrecombinant=r_pairs,
            note=("Part of the K-enrichment is MECHANICAL detection power, not biology: n_recombinant "
                  "correlates with n_bubbles (r=%s) and with #copy-pairs C(K,2) (r=%s). K>=3 families "
                  "simply offer the pair-fitting detector more bubbles and more copy-pairs to fit, so the "
                  "concentration is CONSISTENT WITH the recombination obstruction but is not purely "
                  "biological." % (r_bubbles, r_pairs))),
        recombinant_read_mass_by_K=dict(read_mass),
        credible_molecule_mass_by_K=dict(loc_mass),
        kge3_fraction_of_families=round(n_kge3_fam / max(len(fam_summaries), 1), 3),
        kge3_fraction_of_recombinant_read_mass=round(read_mass["K>=3"] / max(sum(read_mass.values()), 1), 3),
        credible_read_SUN_tier_bridges={"|".join(k): v for k, v in sorted(tier_pairs.items())},
        copy_level_recomb_vulnerable_families=vuln_fams,
        copy_level_recomb_vulnerable_families_carrying_reads=vuln_with_reads,
        interpretation=(
            "Recombinant reads ARE enriched in K>=3 families (Fisher p=%s; K>=3 are %.0f%% of families "
            "but carry %.0f%% of recombinant read-mass) -- consistent with recombination as the K>=3 "
            "obstruction (more copy-pairs => more spoofable combinations), PARTLY offset by a mechanical "
            "detection-power confound (n_recombinant~n_bubbles r=%s, ~C(K,2) r=%s). The credible tracts do "
            "NOT sit on the COPY-level recomb-vulnerable (SUN Tier-2) copies: the sole genome-wide Tier-2 "
            "family (fam42 LOC129529768) carries ZERO recombinant reads, and %d/%d credible reads bridge "
            "two Tier-1 SUN-IDENTIFIABLE copies. This exposes a READ-level K-frontier: a single recombinant "
            "read carries copy A's private SUN in one arm and copy B's in the other, satisfying TWO "
            "Strong-Separation witnesses at once -- impossible for any real single copy. SUN identifiability "
            "guarantees a NON-recombinant read is taggable; it does NOT protect against a recombinant read, "
            "which belongs to no single copy and must abstain."
            % ((("%.4g" % p) if p else "NA"),
               100 * n_kge3_fam / max(len(fam_summaries), 1),
               100 * read_mass["K>=3"] / max(sum(read_mass.values()), 1),
               r_bubbles, r_pairs,
               tier_pairs.get(("1", "1"), 0), len(credible))))

    # ---------- (3) GATE connection ----------
    gate_by_struct = {}
    for sc in sorted(by_struct):
        sub = [r for r in rows if r["structural_call"] == sc]
        assigned = sum(1 for r in sub if r["gate_status"] == "assigned")
        gate_by_struct[sc] = dict(total=len(sub), force_assigned_leak=assigned,
                                  abstained=len(sub) - assigned,
                                  leak_pct=round(100 * assigned / max(len(sub), 1), 1))
    total_leak = sum(1 for r in rows if r["gate_leak"] == 1)
    cred_leak = sum(1 for r in credible if r["gate_leak"] == 1)
    part3 = dict(
        gate_status=dict(collections.Counter(r["gate_status"] for r in rows)),
        force_assigned_leak_total=total_leak,
        abstained_total=len(rows) - total_leak,
        leak_pct_total=round(100 * total_leak / max(len(rows), 1), 1),
        by_structural_call=gate_by_struct,
        credible_clean_conversion_leak="%d/%d (%.0f%%)" % (
            cred_leak, len(credible), 100 * cred_leak / max(len(credible), 1)),
        payoff=("The significance gate FORCE-ASSIGNS %d/%d (%.0f%%) recombinant reads to a single copy, "
                "including %d/%d (%.0f%%) of the credible clean-tract gene-conversion reads -- a concrete "
                "per-read mis-assignment. Per the K-frontier these reads spoof a copy combination and "
                "belong to no single copy; they SHOULD abstain. The gate's SUN/Strong-Sep match "
                "force-assigns to whichever copy's private marker falls in the covered arm = the "
                "recombination failure mode. This is the ROBUST, artifact-independent payoff: whether a "
                "switched read is real conversion, an unrepresented recombinant paralog, or a recurrent "
                "chimera, it conflicts with every single copy by >=%d bubbles and must NOT be force-assigned."
                % (total_leak, len(rows), 100 * total_leak / max(len(rows), 1),
                   cred_leak, len(credible), 100 * cred_leak / max(len(credible), 1), MIN_SINGLE_MM)))

    # ---------- (4) value over the shipped mosaic_discriminator ----------
    part4 = dict(
        axis=("mosaic_discriminator (mosaic.rs::classify_event, --vg RtSwitchArtifact suppression) "
              "classifies ISOFORM EMISSION -- which junctions to emit; it never touches copy ASSIGNMENT, "
              "so it cannot see that a recombinant read is being force-assigned to a copy. VG-native "
              "read-path detection lives on the ASSIGNMENT axis -- a genuinely new axis."),
        microhomology_leg_inert="RT-switch microhomology leg fired %d/%d: tuned for DISTAL splice-like "
              "donor/acceptor direct repeats, inert for ADJACENT-PSV paralog switches. Run naked, the "
              "shipped classifier would label all %d confirmed reads GeneConversion." % (
                  mh_fired, len(rows), by_class.get("GeneConversion", 0)),
        new_discriminator="bubble-bounded RECOMB-FRACTION + FULL-PATTERN cleanliness separate %d "
              "credible clean localized tracts (median frac %.3f) from %d systematic copy-ambiguity reads "
              "(median frac %.3f) AND from %d noisy A/B mosaics -- all of which the mosaic discriminator "
              "alone would call GeneConversion. Two case classes it misses." % (
                  len(credible), round(statistics.median([r["recomb_fraction"] for r in credible] or [0]), 3),
                  by_struct.get("systematic_copy_split", 0),
                  round(statistics.median([r["recomb_fraction"] for r in rows
                                           if r["structural_call"] == "systematic_copy_split"] or [0]), 3),
                  len(loc_noisy)),
        exact_boundaries=("bubble-bounding makes the tract endpoints EXACT graph objects (first/last "
                          "switched bubble = named PSV columns) and the recomb-fraction + full-pattern "
                          "switch count computable; the mosaic breakpoints are approximate donor/acceptor "
                          "coords +/-breakpoint_tol=%d with no fraction." % BREAKPOINT_TOL))

    verdict = (
        "VG-native read-path recombination is a REAL but NARROW signal, and the biological headline is "
        "WEAKER than the raw counts suggest. Of %d raw recombinant rows (%d distinct molecules), 62%% are "
        "systematic near-identical-copy ambiguity and CROSSOVERS yield NO credible biological recombination. "
        "The %d recurrent localized tracts shrink further under a FULL-PATTERN cleanliness test: only %d "
        "rows / %d molecules across %d distinct loci are clean bi-copy tracts (dominated by THOC3 and RABL2, "
        "with ANKRD18/LARP1B/SMG1 minor); notably the previously-flagship GSTM2 conversion-hotspot claim "
        "COLLAPSES to 1 clean read of 292 (291 noisy) and is RETRACTED; the remaining %d localized reads "
        "are noisy A/B mosaics whose clean 2-switch tract is a column-subsetting artifact of unreliable "
        "asm20 copy-vectors (~49%% of fam13 discriminating columns carry a copy allele with ~0 read "
        "support). The RT-switch veto is "
        "INERT (0/%d) and DNA is unchecked (0/%d), so RT-template-switch is NOT excluded and the tract "
        "columns are ordinary het sites RNA cannot separate from allelic SNPs without DNA parCN -- "
        "'gene conversion' is a HYPOTHESIS, not a finding. The K-enrichment (Fisher p=%s) is real but "
        "partly a mechanical detection-power confound (n_recombinant~n_bubbles r=%s). The ROBUST, "
        "artifact-independent payoff is OPERATIONAL: the gate force-assigns %d/%d (%.0f%%) recombinant "
        "reads -- including %d/%d of the credible clean tracts -- to a single copy; all of them conflict "
        "with every single copy and SHOULD ABSTAIN. That is the concrete per-read form of the theory's "
        "K-frontier recombination obstruction, and it holds whether the switch is real conversion, an "
        "unrepresented recombinant paralog, or a recurrent chimera." % (
            len(rows), nmol(rows), len(loc), len(credible), nmol(credible), len(surv), len(loc_noisy),
            len(rows), len(rows), (("%.4g" % p) if p else "NA"), r_bubbles,
            total_leak, len(rows), 100 * total_leak / max(len(rows), 1), cred_leak, len(credible)))

    evaluation = dict(
        q1_survivors=part1, q2_k_frontier=part2, q3_gate=part3, q4_vs_mosaic_discriminator=part4,
        verdict=verdict)
    build["evaluation"] = evaluation
    json.dump(build, open(OUT_JSON, "w"), indent=2)

    # ---- eval headline ----
    print("\n=== EVAL STAGE (cross-referenced with sun_identifiability K-tier ladder) ===")
    print("(1) SURVIVORS post RT-switch veto (INERT) + structural + FULL-PATTERN cleanliness")
    print(f"    RT-switch microhomology hits     : {mh_fired}/{len(rows)}  (leg INERT for adjacent PSVs)")
    print(f"    DNA-supported                    : {dna_yes}/{len(rows)}  (DNFAM id-space mismatch)")
    print(f"    raw recombinant                  : {len(rows)} rows / {nmol(rows)} molecules "
          f"(crossover {sum(1 for r in rows if r['switch_type']=='crossover')}, "
          f"conversion {sum(1 for r in rows if r['switch_type']=='conversion')})")
    print(f"    localized tracts                 : {len(loc)} rows -> clean_full {len(loc_clean)} / "
          f"noisy-mosaic {len(loc_noisy)} (full-switch median clean={part1['full_switch_median_credible']} "
          f"vs noisy={part1['full_switch_median_noisy']})")
    print(f"    CREDIBLE clean gene-conversion   : {len(credible)} rows / {nmol(credible)} molecules / "
          f"{len(surv)} distinct loci ; credible crossover: "
          f"{sum(1 for r in credible if r['switch_type']=='crossover')}")
    for s in survivors:
        print(f"      locus{s['locus']:<4d} {s['gene']:<16s} K={s['K']:<2} reads={s['reads']:<4} "
              f"tract_bubbles={s['tract_bubbles']} recFrac~{s['recomb_fraction_median']} "
              f"fullSw~{s['full_switch_median']}")
    print("(2) K-FRONTIER")
    print(f"    families carrying recombinants   : K>=3 {kb['K>=3'][0]}/{sum(kb['K>=3'])}  "
          f"K==2 {kb['K==2'][0]}/{sum(kb['K==2'])}  Fisher odds={odds and round(odds,2)} p={p and '%.4g'%p}")
    print(f"    detection-power confound         : n_recomb~n_bubbles r={r_bubbles}, ~C(K,2) r={r_pairs} "
          f"(part of the K-enrichment is mechanical)")
    print(f"    recombinant read mass by K       : {dict(read_mass)}  "
          f"(K>=3 = {round(100*read_mass['K>=3']/max(sum(read_mass.values()),1))}% of mass, "
          f"{round(100*n_kge3_fam/max(len(fam_summaries),1))}% of families)")
    print(f"    credible-read SUN tier bridges   : "
          f"{ {'|'.join(k):v for k,v in sorted(tier_pairs.items())} }")
    print(f"    copy-level recomb-vulnerable fams: {vuln_fams}  (carrying reads: {vuln_with_reads})")
    print("(3) GATE")
    print(f"    force-assigned LEAK (should abstain): {total_leak}/{len(rows)} "
          f"({round(100*total_leak/max(len(rows),1),1)}%)  abstained {len(rows)-total_leak}")
    for sc, g in gate_by_struct.items():
        print(f"      {sc:22s} leak {g['force_assigned_leak']}/{g['total']} ({g['leak_pct']}%)")
    print("(4) vs mosaic_discriminator: MH leg inert; recomb-fraction is the new leg; assignment axis.")
    print(f"\nVERDICT: {verdict}")
    print(f"\n[+] appended 'evaluation' block to {OUT_JSON}")
    return evaluation


# =========================================================================== main
def main():
    args = sys.argv[1:]
    # This detector CHARACTERIZES the RAW significance gate (it is the evidence that the gate leaks the
    # 674 recombinant reads). It must therefore observe the gate WITHOUT the abstain leg it motivates, so
    # its gate_status column reflects the pre-abstain force-assignment. The FIX (recombinant_abstain.py)
    # is default-on for O2 consumers; here we opt out. (Deterministic; env set before any materialize.)
    os.environ["RUSTLE_NO_RECOMBINANT_ABSTAIN"] = "1"
    if args and args[0] in ("eval", "--eval"):
        evaluate()
        return
    limit = None
    fam_ids = []
    i = 0
    while i < len(args):
        if args[i] == "--limit":
            limit = int(args[i + 1]); i += 2
        else:
            fam_ids.append(int(args[i])); i += 1

    fams = o2vg.load_families()
    all_ids = fam_ids if fam_ids else sorted(fams)
    genome = pysam.FastaFile(pg.GENOME)

    all_rows = []
    fam_summaries = []
    skipped = collections.Counter()
    processed = 0
    for k, fid in enumerate(all_ids):
        regions = pg.dedup_copies(fams[fid])
        if len(regions) < 2:
            skipped["<2_copies"] += 1
            continue
        try:
            vg = o2vg.materialize_family(fid, use_cache=True)
        except Exception as ex:
            skipped["materialize_error"] += 1
            sys.stderr.write(f"[fam {fid}] materialize error: {ex}\n{traceback.format_exc()}\n")
            continue
        if vg["n_bubbles"] < 2:
            skipped["<2_bubbles"] += 1
            continue
        processed += 1
        rows, fsum = process_family(vg, genome)
        if fsum is not None:
            fam_summaries.append(fsum)
        all_rows.extend(rows)
        if (k + 1) % 20 == 0:
            print(f"  ... {k+1}/{len(all_ids)} families ({processed} with >=2 bubbles, "
                  f"{len(all_rows)} recombinant reads)", flush=True)
        if limit is not None and processed >= limit:
            break
    genome.close()

    # ---- write per-read TSV ----
    cols = ["family", "gene", "copyA", "geneA", "chromA", "copyB", "geneB", "chromB", "read",
            "switch_type", "n_switches", "full_switches", "n_full_ab", "clean_full",
            "n_info_bubbles", "impurity", "tract_bubbles",
            "tract_chrom", "tract_start", "tract_end", "tract_bp", "switch_brackets",
            "microhomology", "cluster_id", "cluster_size", "recurrent", "recomb_fraction",
            "structural_call", "dna_support", "classification", "gate_status", "should_abstain",
            "gate_leak"]
    with open(OUT_TSV, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in sorted(all_rows, key=lambda x: (x["family"], x["cluster_id"], x["read"])):
            o.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---- aggregate ----
    n_cross = sum(1 for r in all_rows if r["switch_type"] == "crossover")
    n_conv = sum(1 for r in all_rows if r["switch_type"] == "conversion")
    by_class = collections.Counter(r["classification"] for r in all_rows)
    # post RT-switch veto: real = GeneConversion; artifact = RtSwitchArtifact; sporadic = ChimeraSuspect
    cross_class = collections.Counter(r["classification"] for r in all_rows if r["switch_type"] == "crossover")
    conv_class = collections.Counter(r["classification"] for r in all_rows if r["switch_type"] == "conversion")
    tract_dist = collections.Counter(r["tract_bubbles"] for r in all_rows if r["switch_type"] == "conversion")
    by_structural = collections.Counter(r["structural_call"] for r in all_rows)
    # credible localized gene-conversion tracts = the honest headline subset
    localized = [r for r in all_rows if r["structural_call"] == "localized_tract"]
    localized_clusters = sorted({r["cluster_id"] for r in localized})
    localized_tract_dist = collections.Counter(r["tract_bubbles"] for r in localized)
    localized_fams = sorted({r["family"] for r in localized})
    fams_with = sorted({r["family"] for r in all_rows})
    # recurrent (confirmed) clusters that are REAL gene-conversion (post-veto)
    real_clusters = sorted({r["cluster_id"] for r in all_rows if r["classification"] == "GeneConversion"})
    gate_leaks = sum(1 for r in all_rows if r["gate_leak"])
    gate_status_dist = collections.Counter(r["gate_status"] for r in all_rows)
    # K-frontier tie-in: families with K>=3 carrying recombinant reads
    kge3 = sorted({fs["family"] for fs in fam_summaries if fs["K"] >= 3 and fs["n_recombinant"] > 0})

    summary = dict(
        families_processed=processed, families_with_recombinants=len(fams_with),
        total_recombinant_reads=len(all_rows),
        crossover_reads=n_cross, conversion_reads=n_conv,
        by_classification=dict(by_class),
        by_structural_call=dict(by_structural),
        localized_tract_reads=len(localized), localized_tract_clusters=len(localized_clusters),
        localized_tract_families=localized_fams,
        localized_tract_length_bubbles=dict(sorted(localized_tract_dist.items())),
        crossover_by_class=dict(cross_class), conversion_by_class=dict(conv_class),
        real_geneconversion_reads=by_class.get("GeneConversion", 0),
        rt_switch_artifact_reads=by_class.get("RtSwitchArtifact", 0),
        chimera_suspect_reads=by_class.get("ChimeraSuspect", 0),
        ambiguous_reads=by_class.get("Ambiguous", 0),
        n_real_geneconversion_clusters=len(real_clusters),
        tract_length_bubbles_distribution=dict(sorted(tract_dist.items())),
        families_with_recombinants_ids=fams_with,
        Kge3_families_with_recombinants=kge3,
        gate_leaks_force_assigned=gate_leaks,
        gate_status_of_recombinant_reads=dict(gate_status_dist),
        skipped=dict(skipped),
        thresholds=dict(MIN_INFO=MIN_INFO, MIN_ARM=MIN_ARM, MIN_FLANK=MIN_FLANK,
                        MIN_TRACT=MIN_TRACT, MIN_SINGLE_MM=MIN_SINGLE_MM, IMPURITY_MAX=IMPURITY_MAX,
                        FAMILY_MIN_SUPPORTING_READS=FAMILY_MIN_SUPPORTING_READS,
                        BREAKPOINT_TOL=BREAKPOINT_TOL, MH_KMIN=MH_KMIN, MH_KMAX=MH_KMAX),
    )
    json.dump(dict(summary=summary, families=fam_summaries), open(OUT_JSON, "w"), indent=2)

    # ---- stdout headline ----
    print("\n=== VG-NATIVE READ-PATH RECOMBINATION (co-located families, >=2 copies, >=2 bubbles) ===")
    print(f"  families processed (>=2 bubbles) : {processed}")
    print(f"  families carrying recombinant reads: {len(fams_with)}  -> {fams_with}")
    print(f"  recombinant reads (PRE veto)     : {len(all_rows)}  "
          f"(crossover {n_cross}, conversion {n_conv})")
    print(f"  classification (RT-switch veto)  : {dict(by_class)}")
    print(f"    crossover  : {dict(cross_class)}")
    print(f"    conversion : {dict(conv_class)}")
    print(f"  STRUCTURAL call (recomb-fraction): {dict(by_structural)}")
    print(f"    -> CREDIBLE localized gene-conversion tracts: {len(localized)} reads in "
          f"{len(localized_clusters)} cluster(s), families {localized_fams}")
    print(f"       localized tract length (bubbles): {dict(sorted(localized_tract_dist.items()))}")
    print(f"  POST-veto REAL gene-conversion   : {by_class.get('GeneConversion',0)} reads "
          f"in {len(real_clusters)} recurrent cluster(s) {real_clusters}")
    print(f"  RT-switch artifact (microhomology): {by_class.get('RtSwitchArtifact',0)} reads")
    print(f"  ChimeraSuspect (sporadic one-off): {by_class.get('ChimeraSuspect',0)} reads")
    print(f"  conversion tract length (bubbles): {dict(sorted(tract_dist.items()))}")
    print(f"  K>=3 families w/ recombinants    : {kge3}")
    print(f"  O2 GATE on recombinant reads     : {dict(gate_status_dist)}")
    print(f"    -> gate LEAKS (force-assigned to ONE copy, SHOULD abstain): {gate_leaks}")
    print(f"  skipped: {dict(skipped)}")
    print(f"\n[+] wrote {OUT_TSV}\n[+] wrote {OUT_JSON}")
    # EVAL STAGE: cross-reference SUN K-tiers, gate, mosaic; appends 'evaluation' to OUT_JSON.
    evaluate()
    return summary


if __name__ == "__main__":
    os.environ.setdefault("PYTHONHASHSEED", "0")
    main()
