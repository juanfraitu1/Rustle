#!/usr/bin/env python3
"""recombinant_abstain.py -- the RECOMBINANT-READ ABSTAIN leg of the O2 copy-assignment gate.

WHY (bench/VG_READ_PATH_RECOMBINATION.md, commit 7e3be89): the significance gate
(copy_assign.assign_read, threaded per read in o2_vg_visualization.materialize_family) FORCE-ASSIGNS
2214/2318 (95.5%) recombinant reads -- including 674/674 (100%) of the bubble-bounded localized_tract
reads -- to a SINGLE copy. A recombinant read carries copy A's private SUN in one bubble-bounded arm and
copy B's in another arm; it satisfies TWO Strong-Separation witnesses at once, so it belongs to NO single
copy = a READ-LEVEL K-frontier (THEORY.md K>=3 recombination obstruction). SUN identifiability guarantees
a NON-recombinant read is taggable; it gives NO protection to a recombinant read. The min_p / margin
certificate assumes a read's decisive SUNs all point to ONE copy -- violated here -> the read must
ABSTAIN, not be force-assigned to whichever copy has more SUNs.

WHAT THESE READS ARE (honest label, NOT "gene-conversion"): the correct name for the abstained set is
"belongs-to-no-copy reads carrying a subset-clean bi-copy recombinant pattern". They are NOT credible
gene-conversion events. Of the 674 localized_tract reads only 134 are clean over the FULL read-supported
discriminating pattern (clean_full=1); 540 are noisy A/B mosaics whose 2-switch tract is a
column-subsetting artifact of unreliable asm20 copy-vectors. Over ALL 2214 abstained reads: 621 clean /
1593 noisy. The former flagship GSTM2 (fam13) collapses to 1 clean / 291 noisy localized reads and its
gene-conversion claim is RETRACTED (see bench/VG_READ_PATH_RECOMBINATION.md, same commit). This does NOT
weaken the leg: abstaining is correct for the full 2214 because EVERY abstained read -- clean or noisy --
conflicts with every single copy by >=2 read-supported bubbles (min_single_mm >= MIN_SINGLE_MM), i.e.
belongs to NO single copy. The leg deliberately does NOT gate on clean_full: gating on it would wrongly
RE-ADMIT the 540 noisy belongs-to-no-copy reads to force-assignment.

WHAT THIS MODULE IS. The SINGLE SOURCE of the recombinant-detection logic. It hosts (verbatim, moved
here so there is exactly one definition -- NOT re-derived) the bubble-bounded copy-switch test
`detect_read_recombination` + `full_pattern_switches` + their thresholds, previously inline in
bench/vg_read_path_recombination.py. That detector now IMPORTS them from here, so the shipped abstain
leg and the diagnostic detector run byte-identical code. On top it adds:
  * build_supported(reads, n_bubbles)  -- factored read-supported-allele pileup (drops copy-vector
    miscalls + singleton errors), identical to process_family's inline computation.
  * is_recombinant(obs, supported, copy_vecs, names) -> (bool, kind, arms)  -- the reusable predicate.
  * apply_abstain_to_vg(vg) -- the DEFAULT-ON gate leg: for every read the gate currently ASSIGNS that
    is_recombinant, relabel status 'assigned' -> 'recombinant' (abstain) and drop best_copy. A pure
    MONOTONE ADDITION OF ABSTENTIONS: non-recombinant reads AND already-abstaining reads are untouched
    bit-for-bit; only force-assigned recombinant reads change. per_copy / assignment counters are LEFT
    as the pure-gate assignment on purpose (they feed the detector's recomb-fraction; abstaining must
    not feed back into detection).

OPT-OUT: DEFAULT-ON. Disable with env RUSTLE_NO_RECOMBINANT_ABSTAIN=1 (the CLI flag
--no-recombinant-abstain sets it), mirroring --no-repeat-gate / --no-split-recombinants. When disabled
apply_abstain_to_vg is a no-op and the assignment output is byte-identical to the prior gate.

HONEST TRADE: abstaining raises CORRECTNESS (removes force-assignment of belongs-to-no-copy reads) at a
COVERAGE cost (more abstentions). RNA alone cannot tell allelic gene-conversion from an unrepresented
recombinant paralog (needs DNA parCN) -- but BOTH must abstain, so the leg is correct regardless.

Deterministic (no RNG; sorted iteration). Run with /home/juanfra/miniforge3/bin/python.
"""
import collections
import itertools
import os

# ---- detection thresholds (documented; conservative for a HiFi ERR~0.003 single molecule).
#      MOVED verbatim from vg_read_path_recombination.py -- single source of truth. ----
MIN_INFO = 4          # read must span >= this many pair-discriminating bubbles matched to A or B
MIN_ARM = 2           # crossover: each side of the switch >= this many discriminating bubbles
MIN_FLANK = 2         # conversion: each flank (the outer A blocks) >= this many bubbles
MIN_TRACT = 2         # conversion: tract (inner B block) >= this many bubbles (a genuine multi-marker
                      #   segment; a 1-bubble "tract" is a single isolated variant = point event, not a
                      #   recombination tract, and is under-determined vs a shared SNP / error -> excluded)
MIN_SINGLE_MM = 2     # NO single copy may explain the read within (MIN_SINGLE_MM-1) mismatches
IMPURITY_MAX = 1      # allow <= this many read bases matching NEITHER paired allele (dropped from split)
# ROBUSTNESS: a switch must be between two GENUINELY CO-EXPRESSED alleles, not a copy-vector miscall.
# A bubble is used only at cols where the copy alleles are themselves seen in the READ pileup -- this
# drops asm20 copy-vector artifacts (an allele with ~0 read support) and singleton sequencing errors.
# supported(base) = read-count >= max(ALLELE_SUPPORT_ABS, ALLELE_SUPPORT_FRAC*cov).
ALLELE_SUPPORT_ABS = 3
ALLELE_SUPPORT_FRAC = 0.02

# ---- opt-out (default-ON) ----
OPTOUT_ENV = "RUSTLE_NO_RECOMBINANT_ABSTAIN"


def abstain_enabled():
    """The recombinant-abstain leg is DEFAULT-ON; RUSTLE_NO_RECOMBINANT_ABSTAIN=1 disables it."""
    return os.environ.get(OPTOUT_ENV, "0") not in ("1", "true", "True", "TRUE")


# =========================================================================== per-read detection
def full_pattern_switches(spanned, obs, cvA, cvB):
    """HONESTY LEG (adversarial review): switch count over ALL read-supported discriminating columns
    for the chosen pair -- NOT just the both-copy-read-supported subset the split is fit on. `spanned`
    is already restricted to columns where the READ base is read-supported (no singleton errors). Here we
    keep every column where the two copies differ (cvA!=cvB) regardless of whether the copy allele is
    itself read-supported, classify the read A/B/neither, and count switches on the A/B subsequence.
    A "clean 2-switch tract" fit on the both-supported subset can hide many switches present in the full
    pattern (the detector's split-column subset is a SUBSET of these columns, so full >= subset always).
    Returns (full_switch_count, n_full_AB_cols)."""
    syms = []
    for c in spanned:
        a, b = cvA.get(c), cvB.get(c)
        if a is None or b is None or a == b:
            continue
        r = obs[c]
        if r == a and r != b:
            syms.append("A")
        elif r == b and r != a:
            syms.append("B")
        # matches neither paired allele -> not an A/B discriminator here; skip
    sw = sum(1 for i in range(1, len(syms)) if syms[i] != syms[i - 1])
    return sw, len(syms)


def detect_read_recombination(obs, supported, copy_vecs, names):
    """Return the best recombinant split for one read, or None.

    obs[col] -> read base; supported[col] -> set of READ-supported alleles at that bubble;
    copy_vecs[copy][col] -> that copy's allele. Only bubbles whose read base is a supported allele are
    used (drops singleton sequencing errors); a pair (A,B) may only switch at cols where BOTH copies'
    alleles are read-supported (drops copy-vector miscalls). A recombinant = NO single copy explains it
    (min single-copy mismatch >= MIN_SINGLE_MM) AND some pair (A,B) admits a PURE ordered A|B (or A|B|A)
    split with 1 or 2 switches. Bubbles bound the tract.
    """
    spanned = [c for c in sorted(obs) if obs[c] in supported.get(c, ())]
    if len(spanned) < MIN_INFO:
        return None
    # single-copy mismatch, counted only where the copy's own allele is read-supported (fair comparison)
    single_mm = {}
    for nm in names:
        cv = copy_vecs[nm]
        single_mm[nm] = sum(1 for c in spanned
                            if cv.get(c) in supported.get(c, ()) and obs[c] != cv[c])
    if min(single_mm.values()) < MIN_SINGLE_MM:
        return None                                   # a single copy already explains this read

    best = None
    best_key = None
    for A, B in itertools.combinations(names, 2):
        cvA, cvB = copy_vecs[A], copy_vecs[B]
        seq = []          # [(col, 'A'|'B')] over bubbles where A,B differ and the read matches one
        imp = 0
        for c in spanned:
            a, b = cvA.get(c), cvB.get(c)
            sup = supported.get(c, ())
            if a is None or b is None or a == b or a not in sup or b not in sup:
                continue                              # not a discriminating, read-supported bubble
            r = obs[c]
            if r == a and r != b:
                seq.append((c, "A"))
            elif r == b and r != a:
                seq.append((c, "B"))
            else:
                imp += 1                              # matches neither paired allele (error/other copy)
        if imp > IMPURITY_MAX or len(seq) < MIN_INFO:
            continue
        syms = [s for _, s in seq]
        sw = [i for i in range(1, len(syms)) if syms[i] != syms[i - 1]]
        if len(sw) == 1:
            i = sw[0]
            if min(i, len(syms) - i) < MIN_ARM:
                continue
            cand = dict(type="crossover", A=A, B=B, seq=seq, imp=imp, switches=1,
                        bnd_idx=[i], tract_idx=None)
        elif len(sw) == 2:
            i1, i2 = sw
            if syms[0] != syms[-1]:                   # must revert (X Y X); (auto for 2 symbols/2 sw)
                continue
            flank1, tract, flank2 = i1, i2 - i1, len(syms) - i2
            if flank1 < MIN_FLANK or flank2 < MIN_FLANK or tract < MIN_TRACT:
                continue
            cand = dict(type="conversion", A=A, B=B, seq=seq, imp=imp, switches=2,
                        bnd_idx=[i1, i2], tract_idx=[i1, i2 - 1])
        else:
            continue
        # prefer fewer switches, then more informative bubbles, then less impurity, then copy order
        key = (cand["switches"], -len(seq), imp, A, B)
        if best_key is None or key < best_key:
            best, best_key = cand, key
    if best is not None:
        full_sw, n_full = full_pattern_switches(spanned, obs, copy_vecs[best["A"]], copy_vecs[best["B"]])
        best["full_switches"] = full_sw
        best["n_full_ab"] = n_full
        # CLEAN over the full read-supported discriminating pattern iff no extra switch was hidden by
        # dropping the copy-unsupported columns (full == the fit split's switch count). Otherwise the
        # read is a noisy A/B mosaic whose "clean tract" is a column-subsetting artifact.
        best["clean_full"] = int(full_sw == best["switches"])
    return best


# =========================================================================== reusable predicate
def build_supported(reads, n_bubbles):
    """Read-supported alleles per bubble (factored verbatim from process_family). `reads` = the
    materialized VG's reads list ([{obs:{col->base}, ...}]). Returns {col -> set(supported bases)}:
    a base is supported iff its read-count >= max(ALLELE_SUPPORT_ABS, ALLELE_SUPPORT_FRAC*coverage).
    This drops asm20 copy-vector miscalls (a copy allele with ~0 read support) + singleton errors."""
    pile = collections.defaultdict(collections.Counter)
    for r in reads:
        for k, v in r["obs"].items():
            if v in "ACGT":
                pile[int(k)][v] += 1
    supported = {}
    for j in range(n_bubbles):
        cov = sum(pile[j].values())
        thr = max(ALLELE_SUPPORT_ABS, ALLELE_SUPPORT_FRAC * cov)
        supported[j] = {b for b, n in pile[j].items() if n >= thr}
    return supported


def is_recombinant(obs, supported, copy_vecs, names):
    """The reusable recombinant predicate. Returns (bool, kind, arms):
      bool = the read is bubble-bounded recombinant (belongs to NO single copy: a clean 1-switch
             crossover A..A|B..B or a bounded 2-switch tract A..A|B..B|A..A, NOT interspersed noise);
      kind = 'crossover' | 'conversion' | None;
      arms = dict(A, B, switches, clean_full, seq) describing the copy pair + the colinear arms, or None.
    Thin wrapper over detect_read_recombination -- the SAME test the diagnostic detector uses; no
    re-derivation."""
    cand = detect_read_recombination(obs, supported, copy_vecs, names)
    if cand is None:
        return False, None, None
    arms = dict(A=cand["A"], B=cand["B"], switches=cand["switches"],
                clean_full=cand.get("clean_full", 1), seq=cand["seq"])
    return True, cand["type"], arms


# =========================================================================== the gate leg
ABSTAIN_STATUS = "recombinant"   # the abstain label (distinct from assigned/ambiguous/tied/no_cover)


def apply_abstain_to_vg(vg):
    """DEFAULT-ON recombinant-abstain leg over a materialized family VG (o2_vg_visualization).

    For every read the significance gate currently ASSIGNS (status == 'assigned') that is_recombinant,
    relabel status -> 'recombinant' (abstain), record the pre-abstain gate call in
    read['gate_status_pre_abstain'] and the switch kind in read['recombinant_kind'], and drop best_copy.
    PURE MONOTONE ADDITION OF ABSTENTIONS: non-recombinant reads and already-abstaining reads are left
    bit-for-bit unchanged; only force-assigned recombinant reads change. The vg['assignment'] and
    vg['per_copy'] SUMMARY counters are kept consistent (each flipped read is moved out of 'assigned'
    and out of its copy's quant into 'recombinant' -- 'they SHOULD abstain, not be counted toward one
    copy's quant'). This never feeds back into the diagnostic detector's recomb-fraction because the
    detector runs with the leg OPTED OUT (RUSTLE_NO_RECOMBINANT_ABSTAIN=1), so it always sees the pure
    per_copy; the counters are updated only on the enabled path.

    Returns a list of dicts (one per flipped read) for reporting; empty when disabled or no leaks.
    No-op (returns []) when RUSTLE_NO_RECOMBINANT_ABSTAIN=1.
    """
    if not abstain_enabled():
        return []
    names = [c["name"] for c in vg["copies"]]
    bubbles = vg["bubbles"]
    if len(names) < 2 or len(bubbles) < 2:
        return []
    reads = vg["reads"]
    assignment = vg.get("assignment")
    per_copy = vg.get("per_copy")
    copy_vecs = {nm: {j: bubbles[j]["hap"].get(nm) for j in range(len(bubbles))} for nm in names}
    supported = build_supported(reads, len(bubbles))
    flipped = []
    for r in reads:
        if r.get("status") != "assigned":
            continue                                  # only ADD abstentions on force-assigned reads
        obs = {int(k): v for k, v in r["obs"].items() if v in "ACGT"}
        if len(obs) < MIN_INFO:
            continue
        rec, kind, arms = is_recombinant(obs, supported, copy_vecs, names)
        if not rec:
            continue
        prev_best = r.get("best_copy")
        flipped.append(dict(family=vg["family"], read=r["name"], kind=kind,
                            copyA=arms["A"], copyB=arms["B"], switches=arms["switches"],
                            clean_full=arms["clean_full"], prev_best=prev_best))
        r["gate_status_pre_abstain"] = r["status"]
        r["recombinant_kind"] = kind
        r["status"] = ABSTAIN_STATUS
        r["best_copy"] = None
        # keep the summary counters consistent: move the read out of assigned/copy-quant into recombinant
        if isinstance(assignment, dict):
            assignment["assigned"] = assignment.get("assigned", 0) - 1
            assignment[ABSTAIN_STATUS] = assignment.get(ABSTAIN_STATUS, 0) + 1
        if isinstance(per_copy, dict) and prev_best in per_copy:
            per_copy[prev_best] -= 1
    return flipped
