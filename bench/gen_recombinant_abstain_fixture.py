#!/usr/bin/env python
"""BYTE-PARITY golden-fixture generator for the Rust port of the O2 recombinant-read
ABSTAIN gate `bench/recombinant_abstain.py` -> `vg_family/recombinant_abstain.rs`
(the O2 leg of the family-def migration).

It imports the SHIPPED `recombinant_abstain` (RA) module and the cached, materialized
family VGs in /home/juanfra/winloci_scratch/o2vg/fam<F>.json (the PURE-gate,
PRE-abstain VG dict written by `bench/o2_vg_visualization.materialize_family`). No
BAM / genome / pysam is needed -- the abstain gate is deterministic + pure over the
materialized VG dict.

Three parts:
  (a) "detection_cases": (obs, supported, copy_vecs, names) -> the REAL
      RA.detect_read_recombination output (type/A/B/switches/clean_full/full_switches/
      n_full_ab/seq/imp), covering crossover (1-switch), conversion (2-switch tract),
      non-recombinant (single-copy explains / <MIN_INFO / impurity>1 / arm<MIN_ARM /
      tract<MIN_TRACT / flank boundary), clean_full 0 vs 1, and the pair-selection
      tie-break (>=2 candidate pairs).  `supported`/`copy_vecs` are RESTRICTED to the
      read's obs columns -- detect only reads columns in `spanned` (subset of obs), so
      the restriction is EXACT (asserted here vs the full-family inputs).  A mix of
      REAL reads (from several families) + hand-built SYNTHETIC boundary cases.
  (b) "build_supported_cases": hand-built (reads, n_bubbles) -> RA.build_supported,
      hitting the frac-vs-abs threshold boundary (n >= max(3, 0.02*cov)).  The
      whole-VG cases (c) additionally each ship `expect_supported` for a second
      build_supported parity check on real read piles.
  (c) "vgs": whole family VGs -> RA.apply_abstain_to_vg's flipped list + the
      post-abstain read statuses (status/best_copy/gate_status_pre_abstain/
      recombinant_kind) + the post assignment/per_copy counters.  Includes families
      that FLIP reads (crossover + conversion, clean_full 0 and 1, 2- and multi-copy)
      AND families that flip NONE (>=2 copies/bubbles) AND early-no-op families
      (<2 bubbles, and a synthetic <2-copy VG).

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/gen_recombinant_abstain_fixture.py
Deterministic (sorted selection, PYTHONHASHSEED=0, no RNG).
"""
import os
import sys

# determinism: pin the hash seed BEFORE the heavy imports (re-exec preserves argv)
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

# the abstain leg is DEFAULT-ON; make sure the opt-out env is NOT set while we build.
os.environ.pop("RUSTLE_NO_RECOMBINANT_ABSTAIN", None)

os.chdir(os.path.dirname(os.path.abspath(__file__)) + "/..")   # repo root
import collections
import copy
import json

sys.path.insert(0, "bench")
import recombinant_abstain as RA          # the port target (ground truth)

CACHE = "/home/juanfra/winloci_scratch/o2vg"
OUT = os.path.join("src", "rustle", "vg_family", "testdata", "recombinant_abstain_fixture.json")

# family VGs shipped whole for the apply_abstain parity test (self-contained).
#   flippers (crossover+conversion, clean_full 0/1, 2- and multi-copy) + non-flippers
#   (>=2 copies/bubbles) + early-no-op (<2 bubbles).  A synthetic <2-copy VG is added
#   in code for the names<2 early return.
VG_FAMILIES = [39, 90, 147, 73, 105, 114, 161, 156, 176]

# families to mine REAL detection cases from (mix of kinds / clean_full / copy count)
DETECT_FAMILIES = [39, 90, 147, 73, 105, 114, 13, 48, 33, 109, 9, 17]
MAX_REC_CASES = 260      # cap on recombinant (detect != None) real cases
MAX_NONE_TRIVIAL = 40    # cap on trivial None (obs < MIN_INFO) real cases
MAX_NONE_NONTRIVIAL = 90 # cap on non-trivial None (obs >= MIN_INFO, single-copy explains) real cases


# ------------------------------------------------------------------ helpers
def load_raw(fam):
    """The PRE-abstain cached VG dict (materialize_family writes the pure gate BEFORE
    applying the abstain leg, so the on-disk file is exactly the pre-abstain input)."""
    with open(os.path.join(CACHE, f"fam{fam}.json")) as fh:
        return json.load(fh)


def family_inputs(vg):
    """(names, copy_vecs_full, supported_full) for a VG -- exactly as apply_abstain builds."""
    names = [c["name"] for c in vg["copies"]]
    bubbles = vg["bubbles"]
    nbub = len(bubbles)
    copy_vecs = {nm: {j: bubbles[j]["hap"].get(nm) for j in range(nbub)} for nm in names}
    supported = RA.build_supported(vg["reads"], nbub)
    return names, copy_vecs, supported, nbub


def read_obs(read):
    """The ACGT-filtered int-keyed obs the abstain leg feeds to detect."""
    return {int(k): v for k, v in read["obs"].items() if v in "ACGT"}


def cand_to_expect(cand):
    """Serialize an RA.detect_read_recombination result (dict or None) to the fixture
    'expect' shape (None or the load-bearing output fields)."""
    if cand is None:
        return None
    return dict(
        type=cand["type"], A=cand["A"], B=cand["B"], switches=cand["switches"],
        clean_full=cand["clean_full"], full_switches=cand["full_switches"],
        n_full_ab=cand["n_full_ab"], imp=cand["imp"],
        seq=[[c, s] for (c, s) in cand["seq"]],
    )


def restrict(obs, names, copy_vecs, supported):
    """Restrict supported/copy_vecs to the read's obs columns (detect only reads
    columns in spanned, a subset of obs)."""
    cols = sorted(obs)
    sup_r = {c: sorted(supported.get(c, set())) for c in cols}
    cv_r = {nm: {c: copy_vecs[nm].get(c) for c in cols} for nm in names}
    return sup_r, cv_r


def detect_on_restricted(obs, sup_r, cv_r, names):
    """Run RA.detect on the RESTRICTED inputs (rebuild sets from sorted lists)."""
    sup = {c: set(v) for c, v in sup_r.items()}
    cv = {nm: dict(cv_r[nm]) for nm in names}
    return RA.detect_read_recombination(obs, sup, cv, names)


def n_qualifying_pairs(obs, supported, copy_vecs, names):
    """Count copy pairs that yield a valid crossover/conversion candidate for `obs`
    (to report how many detection cases exercise the >=2-pair tie-break)."""
    import itertools
    empty = ()
    spanned = [c for c in sorted(obs) if obs[c] in supported.get(c, empty)]
    if len(spanned) < RA.MIN_INFO:
        return 0
    single_mm = {}
    for nm in names:
        cv = copy_vecs[nm]
        single_mm[nm] = sum(1 for c in spanned
                            if cv.get(c) in supported.get(c, empty) and obs[c] != cv[c])
    if min(single_mm.values()) < RA.MIN_SINGLE_MM:
        return 0
    n = 0
    for A, B in itertools.combinations(names, 2):
        cvA, cvB = copy_vecs[A], copy_vecs[B]
        seq = []
        imp = 0
        for c in spanned:
            a, b = cvA.get(c), cvB.get(c)
            sup = supported.get(c, empty)
            if a is None or b is None or a == b or a not in sup or b not in sup:
                continue
            r = obs[c]
            if r == a and r != b:
                seq.append((c, "A"))
            elif r == b and r != a:
                seq.append((c, "B"))
            else:
                imp += 1
        if imp > RA.IMPURITY_MAX or len(seq) < RA.MIN_INFO:
            continue
        syms = [s for _, s in seq]
        sw = [i for i in range(1, len(syms)) if syms[i] != syms[i - 1]]
        if len(sw) == 1:
            i = sw[0]
            if min(i, len(syms) - i) < RA.MIN_ARM:
                continue
        elif len(sw) == 2:
            i1, i2 = sw
            if syms[0] != syms[-1]:
                continue
            if i1 < RA.MIN_FLANK or (len(syms) - i2) < RA.MIN_FLANK or (i2 - i1) < RA.MIN_TRACT:
                continue
        else:
            continue
        n += 1
    return n


# ------------------------------------------------------------------ detection cases
def build_detection_cases():
    cases = []
    stats = collections.Counter()
    tiebreak_cases = 0

    def add(label, obs, sup_r, cv_r, names, expect):
        cases.append(dict(label=label, names=list(names),
                          obs={str(c): b for c, b in sorted(obs.items())},
                          supported={str(c): v for c, v in sup_r.items()},
                          copy_vecs={nm: {str(c): cv_r[nm][c] for c in sorted(cv_r[nm])} for nm in names},
                          expect=expect))

    # ---- (1) REAL reads across several families -----------------------------------------------
    n_rec = n_none_triv = n_none_nontriv = 0
    for fam in DETECT_FAMILIES:
        vg = load_raw(fam)
        names, cvf, supf, nbub = family_inputs(vg)
        if len(names) < 2 or nbub < 2:
            continue
        for read in vg["reads"]:
            obs = read_obs(read)
            cand = RA.detect_read_recombination(obs, supf, cvf, names)
            is_rec = cand is not None
            triv = len(obs) < RA.MIN_INFO
            if is_rec:
                if n_rec >= MAX_REC_CASES:
                    continue
                n_rec += 1
            elif triv:
                if n_none_triv >= MAX_NONE_TRIVIAL:
                    continue
                n_none_triv += 1
            else:
                if n_none_nontriv >= MAX_NONE_NONTRIVIAL:
                    continue
                n_none_nontriv += 1
            sup_r, cv_r = restrict(obs, names, cvf, supf)
            # EXACT-restriction sanity: detect on restricted inputs == detect on full inputs
            r2 = detect_on_restricted(obs, sup_r, cv_r, names)
            assert cand_to_expect(r2) == cand_to_expect(cand), (fam, read["name"])
            expect = cand_to_expect(cand)
            add(f"real_fam{fam}_{read['name']}", obs, sup_r, cv_r, names, expect)
            if expect is not None:
                stats[("kind", expect["type"])] += 1
                stats[("clean_full", expect["clean_full"])] += 1
                if len(names) >= 3:
                    stats["multi_copy_rec"] += 1
                if n_qualifying_pairs(obs, supf, cvf, names) >= 2:
                    tiebreak_cases += 1
            else:
                stats["none_trivial" if triv else "none_nontrivial"] += 1

    # ---- (2) SYNTHETIC boundary cases (control every branch deterministically) ----------------
    def syn(label, names, cv, sup, obs):
        """cv: {name:{col:base}}, sup: {col:iterable bases}, obs: {col:base}."""
        supset = {c: set(v) for c, v in sup.items()}
        cvd = {nm: dict(cv[nm]) for nm in names}
        cand = RA.detect_read_recombination(dict(obs), supset, cvd, names)
        sup_r = {c: sorted(supset.get(c, set())) for c in sorted(obs)}
        cv_r = {nm: {c: cvd[nm].get(c) for c in sorted(obs)} for nm in names}
        # restriction sanity
        assert cand_to_expect(detect_on_restricted(obs, sup_r, cv_r, names)) == cand_to_expect(cand), label
        add(label, obs, sup_r, cv_r, names, cand_to_expect(cand))
        stats[("syn", label, "rec" if cand is not None else "none")] += 1
        return cand

    NA = ["cA", "cB"]          # 2-copy synthetic names
    NB = ["cA", "cB", "cC"]    # 3-copy synthetic names

    # S1 crossover, clean_full=1 (all A/B cols both-supported)
    syn("syn_crossover_clean1", NA,
        {"cA": {i: "A" for i in range(8)}, "cB": {i: "C" for i in range(8)}},
        {i: "AC" for i in range(8)},
        {0: "A", 1: "A", 2: "A", 3: "A", 4: "C", 5: "C", 6: "C", 7: "C"})

    # S2 crossover, clean_full=0 (a dropped copy-unsupported col hides a full-pattern switch)
    syn("syn_crossover_clean0", NA,
        {"cA": {i: "A" for i in range(7)}, "cB": {i: "C" for i in range(7)}},
        {0: "AC", 1: "AC", 2: "C", 3: "AC", 4: "AC", 5: "AC", 6: "AC"},  # col2: only C supported
        {0: "A", 1: "A", 2: "C", 3: "A", 4: "C", 5: "C", 6: "C"})

    # S3 conversion (2-switch reverting tract), clean_full=1
    syn("syn_conversion_clean1", NA,
        {"cA": {i: "A" for i in range(8)}, "cB": {i: "C" for i in range(8)}},
        {i: "AC" for i in range(8)},
        {0: "A", 1: "A", 2: "C", 3: "C", 4: "C", 5: "C", 6: "A", 7: "A"})

    # S4 conversion tract too short (tract == 1 < MIN_TRACT) -> None
    syn("syn_conversion_tract1_none", NA,
        {"cA": {i: "A" for i in range(6)}, "cB": {i: "C" for i in range(6)}},
        {i: "AC" for i in range(6)},
        {0: "A", 1: "A", 2: "C", 3: "A", 4: "A", 5: "A"})   # single B blip: flank/tract fail + single-copy

    # S5 <MIN_INFO (only 3 discriminating cols) -> None
    syn("syn_lt_min_info", NA,
        {"cA": {i: "A" for i in range(3)}, "cB": {i: "C" for i in range(3)}},
        {i: "AC" for i in range(3)},
        {0: "A", 1: "C", 2: "A"})

    # S6 single copy explains (cB matches every col) -> None (min single_mm < 2)
    syn("syn_single_copy_explains", NA,
        {"cA": {i: "A" for i in range(6)}, "cB": {i: "C" for i in range(6)}},
        {i: "AC" for i in range(6)},
        {0: "C", 1: "C", 2: "C", 3: "C", 4: "C", 5: "C"})

    # S7 impurity == 1 tolerated (one col matches neither paired allele) -> still detect crossover
    syn("syn_impurity1_ok", NA,
        {"cA": {i: "A" for i in range(8)}, "cB": {i: "C" for i in range(8)}},
        {0: "AC", 1: "AC", 2: "AC", 3: "ACG", 4: "AC", 5: "AC", 6: "AC", 7: "AC"},
        {0: "A", 1: "A", 2: "A", 3: "G", 4: "C", 5: "C", 6: "C", 7: "C"})  # col3 obs=G = impurity

    # S8 impurity == 2 -> pair rejected -> None (only pair present)
    syn("syn_impurity2_none", NA,
        {"cA": {i: "A" for i in range(8)}, "cB": {i: "C" for i in range(8)}},
        {0: "AC", 1: "AC", 2: "ACG", 3: "ACG", 4: "AC", 5: "AC", 6: "AC", 7: "AC"},
        {0: "A", 1: "A", 2: "G", 3: "G", 4: "C", 5: "C", 6: "C", 7: "C"})

    # S9 tie-break by -len(seq) OVERRIDING copy-name order (the load-bearing case): 3 copies, TWO
    #    valid 1-switch crossover pairs (cA,cB) and (cA,cC).  key=(switches,-len,imp,A,B) with -len
    #    BEFORE the (A,B) names, so the LONGER seq wins even when its name pair is name-order-LATER.
    #    Here cB has col0==A (not discriminating vs cA) -> (cA,cB) seq len 9; cC is all-C -> (cA,cC)
    #    seq len 10.  Copy-name order alone would pick (cA,cB) ('cB' < 'cC'); -len picks the longer
    #    (cA,cC).  No single copy explains (min single_mm = 4 >= MIN_SINGLE_MM), so detection fires
    #    (the prior version of this case had obs == cC exactly -> min single_mm 0 -> short-circuit None,
    #    exercising no tie-break; reconstructed 2026-07-03).
    syn("syn_tiebreak_len", NB,
        {"cA": {i: "A" for i in range(10)},
         "cB": {0: "A", **{i: "C" for i in range(1, 10)}},   # col0 == A (not disc.) -> seq len 9
         "cC": {i: "C" for i in range(10)}},                 # all C -> seq len 10 (LONGER, name-later)
        {i: "AC" for i in range(10)},
        {0: "A", 1: "A", 2: "A", 3: "A", 4: "A", 5: "C", 6: "C", 7: "C", 8: "C", 9: "C"})

    # S10 tie-break by copy-name order: two pairs with IDENTICAL key prefix -> min (A,B) wins.
    #     cB and cC are identical discriminators vs cA, so both pairs give the same (switches,-len,
    #     imp); tie broken by ('cA','cB') < ('cA','cC').
    syn("syn_tiebreak_name", NB,
        {"cA": {i: "A" for i in range(8)},
         "cB": {i: "C" for i in range(8)},
         "cC": {i: "C" for i in range(8)}},
        {i: "AC" for i in range(8)},
        {0: "A", 1: "A", 2: "A", 3: "A", 4: "C", 5: "C", 6: "C", 7: "C"})

    # S11 crossover arm exactly MIN_ARM (arm==2 passes) via a 3rd copy that blocks single-copy
    #     explanation; the winning (cA,cB) pair has a 2|6 split.
    syn("syn_arm_boundary_pass", NA,
        {"cA": {i: "A" for i in range(8)}, "cB": {i: "C" for i in range(8)}},
        {i: "AC" for i in range(8)},
        {0: "A", 1: "A", 2: "C", 3: "C", 4: "C", 5: "C", 6: "C", 7: "C"})  # arm=2 (>=MIN_ARM)

    # S12 conversion flank exactly MIN_FLANK on both sides (flank1=2, tract=2, flank2=2)
    syn("syn_flank_boundary_pass", NA,
        {"cA": {i: "A" for i in range(6)}, "cB": {i: "C" for i in range(6)}},
        {i: "AC" for i in range(6)},
        {0: "A", 1: "A", 2: "C", 3: "C", 4: "A", 5: "A"})

    print(f"[detect] {len(cases)} cases "
          f"(real: rec={n_rec} none_trivial={n_none_triv} none_nontrivial={n_none_nontriv}); "
          f">=2-pair tie-break real cases={tiebreak_cases}")
    print(f"         kinds crossover={stats[('kind','crossover')]} conversion={stats[('kind','conversion')]} "
          f"clean_full0={stats[('clean_full',0)]} clean_full1={stats[('clean_full',1)]} "
          f"multi_copy_rec={stats['multi_copy_rec']}")
    return cases, stats, tiebreak_cases


# ------------------------------------------------------------------ build_supported cases
def build_supported_cases():
    """Hand-built (reads, n_bubbles) -> supported, targeting the frac-vs-abs boundary.
    A base is supported iff count >= max(3, 0.02*cov)."""
    cases = []

    def mk(label, piles, n_bubbles):
        """piles: {col: {base: count}} -> synthesize `count` single-col reads per base."""
        reads = []
        for col, bc in piles.items():
            for base, n in bc.items():
                for _ in range(n):
                    reads.append({"obs": {str(col): base}})
        supported = RA.build_supported(reads, n_bubbles)
        cases.append(dict(label=label,
                          reads=[{"obs": {k: v for k, v in r["obs"].items()}} for r in reads],
                          n_bubbles=n_bubbles,
                          expected={str(c): sorted(supported[c]) for c in range(n_bubbles)}))

    # col0: abs floor (cov 4 -> thr=max(3,0.08)=3): n>=3 supported.  A=3 (in), C=2 (out).
    # col1: frac dominates (cov 500 -> thr=max(3,10.0)=10): A=10 (in, ==thr), C=9 (out).
    # col2: frac boundary just below (cov 500 -> thr 10): A=9 (out).
    # col3: empty (cov 0 -> thr 3 -> nothing).
    # col4: exactly abs (n==3, cov 3 -> thr 3): A in.
    mk("boundary_mix",
       {0: {"A": 3, "C": 2},
        1: {"A": 10, "C": 9, "G": 481},
        2: {"A": 9, "C": 491},
        4: {"A": 3}},
       6)  # includes empty cols 3 and 5

    # a non-ACGT base ('N') must be ignored entirely.
    mk("ignore_non_acgt",
       {0: {"A": 5, "N": 100}},
       2)

    # frac exactly equals 3 at cov 150 (0.02*150 = 3.0): max(3, 3.0) -> n>=3.
    mk("frac_equals_abs",
       {0: {"A": 3, "C": 2, "G": 145}},
       1)

    print(f"[build_supported] {len(cases)} synthetic cases")
    return cases


# ------------------------------------------------------------------ whole-VG apply_abstain cases
def build_vg_cases():
    vgs = []

    def emit(vg):
        pre = copy.deepcopy(vg)   # the pre-abstain INPUT (shipped verbatim)
        post = copy.deepcopy(vg)
        flipped = RA.apply_abstain_to_vg(post)   # DEFAULT-ON leg mutates `post`
        names, cvf, supf, nbub = family_inputs(pre)
        expect_supported = ({str(c): sorted(supf[c]) for c in range(nbub)}
                            if nbub >= 2 else {})
        rec = dict(
            family=pre["family"],
            copies=[{"name": c["name"]} for c in pre["copies"]],
            bubbles=[{"hap": dict(b["hap"])} for b in pre["bubbles"]],
            reads=[{"name": r["name"], "obs": dict(r["obs"]),
                    "status": r["status"], "best_copy": r.get("best_copy")}
                   for r in pre["reads"]],
            assignment=dict(pre.get("assignment", {})),
            per_copy=dict(pre.get("per_copy", {})),
            expect_flipped=[dict(family=f["family"], read=f["read"], kind=f["kind"],
                                 copyA=f["copyA"], copyB=f["copyB"], switches=f["switches"],
                                 clean_full=f["clean_full"], prev_best=f["prev_best"])
                            for f in flipped],
            expect_reads=[{"name": r["name"], "status": r["status"],
                           "best_copy": r.get("best_copy"),
                           "gate_status_pre_abstain": r.get("gate_status_pre_abstain"),
                           "recombinant_kind": r.get("recombinant_kind")}
                          for r in post["reads"]],
            expect_assignment=dict(post.get("assignment", {})),
            expect_per_copy=dict(post.get("per_copy", {})),
            expect_supported=expect_supported,
        )
        vgs.append(rec)
        return len(flipped)

    n_flip_total = 0
    for fam in VG_FAMILIES:
        nf = emit(load_raw(fam))
        n_flip_total += nf
        print(f"[vg] fam{fam}: {nf} reads flip")

    # synthetic <2-copy VG -> names<2 early no-op (returns []; VG unchanged bit-for-bit)
    syn_vg = dict(
        family=99991, copies=[{"name": "only0"}],
        bubbles=[{"hap": {"only0": "A"}}, {"hap": {"only0": "C"}}],
        reads=[{"name": "r0", "obs": {"0": "A", "1": "C"}, "status": "assigned", "best_copy": "only0"}],
        assignment={"assigned": 1}, per_copy={"only0": 1},
    )
    emit(syn_vg)
    print("[vg] synthetic 1-copy VG: early no-op (names<2)")

    print(f"[vg] {len(vgs)} VGs; total reads flipped across shipped VGs = {n_flip_total}")
    return vgs, n_flip_total


# ------------------------------------------------------------------ main
def main():
    detection_cases, dstats, tiebreak = build_detection_cases()
    bs_cases = build_supported_cases()
    vgs, n_flip = build_vg_cases()

    fixture = dict(
        thresholds=dict(
            MIN_INFO=RA.MIN_INFO, MIN_ARM=RA.MIN_ARM, MIN_FLANK=RA.MIN_FLANK,
            MIN_TRACT=RA.MIN_TRACT, MIN_SINGLE_MM=RA.MIN_SINGLE_MM,
            IMPURITY_MAX=RA.IMPURITY_MAX, ALLELE_SUPPORT_ABS=RA.ALLELE_SUPPORT_ABS,
            ALLELE_SUPPORT_FRAC=RA.ALLELE_SUPPORT_FRAC,
            OPTOUT_ENV=RA.OPTOUT_ENV, ABSTAIN_STATUS=RA.ABSTAIN_STATUS,
        ),
        detection_cases=detection_cases,
        build_supported_cases=bs_cases,
        vgs=vgs,
    )
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(fixture, fh, sort_keys=True)
    n_rec = sum(1 for c in detection_cases if c["expect"] is not None)
    print(f"[fixture] wrote {os.path.normpath(OUT)} ({os.path.getsize(OUT):,} bytes)")
    print(f"          detection_cases={len(detection_cases)} (recombinant={n_rec}, none={len(detection_cases)-n_rec}); "
          f"build_supported_cases={len(bs_cases)}; vgs={len(vgs)}; total flips={n_flip}")


if __name__ == "__main__":
    main()
