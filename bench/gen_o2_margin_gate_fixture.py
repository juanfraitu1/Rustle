#!/usr/bin/env python3
"""Golden byte-parity fixture generator for the O2 MARGIN gate Rust port
(src/rustle/vg_family/o2_margin_gate.rs, porting bench/copy_assign.py::assign_read --
the gate bench/o2_vg_visualization.py::materialize_family actually calls).

Run deterministically:
    PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/gen_o2_margin_gate_fixture.py

Emits src/rustle/vg_family/testdata/o2_margin_gate_fixture.json. Every float (margin,
per-copy logL, e, jw) is stored as its raw IEEE-754 u64 bits (a JSON integer, exact),
NOT as a decimal (serde_json decimal parse is 1-ULP lossy). Inputs:
  * >=200 REAL reads mined from the cached o2vg VGs (/home/juanfra/winloci_scratch/o2vg/
    fam*.json): each read's obs + the family's copy_vecs reconstructed from bubbles[].hap.
    copy_vecs is PRUNED to the read's spanned columns -- assign_read only ever looks up
    copy_vecs at spanned columns, so the pruned result is IDENTICAL (asserted below).
  * synthetic edge cases: ties (margin 0), single-copy-spanned, obs with None, a >2-copy
    family, missing copy columns, and >=20 junction-term cases (both match / both miss /
    mixed / btol boundary) exercising read_bounds+copy_bounds.
"""
import json
import os
import struct
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import copy_assign as ca  # noqa: E402  the REAL gate under test

O2VG = "/home/juanfra/winloci_scratch/o2vg"
OUT = os.path.join(HERE, "..", "src", "rustle", "vg_family", "testdata", "o2_margin_gate_fixture.json")


def d2bits(x):
    """f64 -> its raw IEEE-754 u64 bits (exact; +inf/-inf handled)."""
    return struct.unpack("<Q", struct.pack("<d", float(x)))[0]


def emit_expect(best, margin, n_span, n_dec, logL, names):
    return dict(
        best=best,
        margin_bits=d2bits(margin),
        n_span=n_span,
        n_decisive=n_dec,
        logl=[[nm, d2bits(logL[nm])] for nm in names],  # names order
    )


def prune_copy_vecs(copy_vecs, names, obs):
    """Restrict each copy's column dict to the read's spanned columns. assign_read only
    reads copy_vecs at spanned columns, so this is result-preserving (asserted below)."""
    spanned = {c for c in obs if obs[c] is not None}
    return [[nm, {str(c): copy_vecs[nm][c] for c in spanned if c in copy_vecs[nm]}] for nm in names]


def make_real_case(cid, obs, copy_vecs, names):
    """obs: {int col -> base}. copy_vecs: {name -> {int col -> base}} (dense). names: order.
    Calls the REAL gate on the FULL copy_vecs, verifies pruning is identical, returns a
    self-contained case with PRUNED copy_vecs."""
    best, margin, n_span, n_dec, logL = ca.assign_read(obs, copy_vecs)
    # self-consistency: pruned copy_vecs must give the identical result (bit-for-bit)
    pruned = {nm: {c: copy_vecs[nm][c] for c in obs if obs[c] is not None and c in copy_vecs[nm]}
              for nm in names}
    b2, m2, s2, d2, l2 = ca.assign_read(obs, pruned)
    assert (b2, d2bits(m2), s2, d2) == (best, d2bits(margin), n_span, n_dec), f"{cid}: prune changed scalars"
    assert all(d2bits(l2[nm]) == d2bits(logL[nm]) for nm in names), f"{cid}: prune changed logL"
    return dict(
        id=cid,
        kind="real",
        obs={str(c): obs[c] for c in sorted(obs)},   # None never occurs in cached obs
        copy_vecs=prune_copy_vecs(copy_vecs, names, obs),
        read_bounds=None,
        copy_bounds=None,
        e_bits=d2bits(ca.ERR),
        jw_bits=d2bits(5.0),
        btol=4,
        expect=emit_expect(best, margin, n_span, n_dec, logL, names),
    )


def make_synth_case(cid, obs, copy_vecs, names, read_bounds=None, copy_bounds=None,
                    e=ca.ERR, jw=5.0, btol=4):
    """Self-contained synthetic case. obs may contain None. read_bounds passed to the REAL
    gate as a SORTED LIST (assign_read only iterates+truthy-checks it) so Python's iteration
    order == the Rust BTreeSet's ascending order -> the junction logL sum matches bit-for-bit.
    copy_bounds passed as a dict {name -> set}."""
    rb_arg = sorted(read_bounds) if read_bounds is not None else None
    cb_arg = {nm: set(s) for nm, s in copy_bounds} if copy_bounds is not None else None
    best, margin, n_span, n_dec, logL = ca.assign_read(
        obs, copy_vecs, read_bounds=rb_arg, copy_bounds=cb_arg, e=e, jw=jw, btol=btol)
    return dict(
        id=cid,
        kind="synthetic",
        obs={str(c): obs[c] for c in sorted(obs)},
        copy_vecs=[[nm, {str(c): copy_vecs[nm][c] for c in sorted(copy_vecs[nm])}] for nm in names],
        read_bounds=(sorted(read_bounds) if read_bounds is not None else None),
        copy_bounds=([[nm, sorted(s)] for nm, s in copy_bounds] if copy_bounds is not None else None),
        e_bits=d2bits(e),
        jw_bits=d2bits(jw),
        btol=btol,
        expect=emit_expect(best, margin, n_span, n_dec, logL, names),
    )


# --------------------------------------------------------------------- REAL reads
# Spread across copy counts {2,3,4,7} and a range of bubble densities.
REAL_FAMS = [
    ("fam4", 40), ("fam0", 30), ("fam24", 30), ("fam41", 26), ("fam9", 22),
    ("fam11", 22), ("fam19", 22), ("fam102", 20), ("fam81", 20), ("fam17", 20),
    ("fam3", 18), ("fam15", 18),
]


def load_family(fam):
    d = json.load(open(os.path.join(O2VG, f"{fam}.json")))
    names = [c["name"] for c in d["copies"]]          # copy0, copy1, ... (materialize order)
    bubbles = d["bubbles"]
    # copy_vecs = {nm: {j: bubbles[j]["hap"][nm] for j in range(len(bubbles))} for nm in names}
    copy_vecs = {nm: {j: bubbles[j]["hap"][nm] for j in range(len(bubbles))} for nm in names}
    return names, copy_vecs, d["reads"]


def collect_real():
    cases = []
    for fam, cap in REAL_FAMS:
        names, copy_vecs, reads = load_family(fam)
        with_obs = [r for r in reads if r.get("obs")]
        if not with_obs:
            continue
        # spread the sample across the read list (diversifies status / n_decisive)
        stride = max(1, len(with_obs) // cap)
        picked = with_obs[::stride][:cap]
        for r in picked:
            obs = {int(k): v for k, v in sorted(r["obs"].items(), key=lambda kv: int(kv[0]))}
            cid = f"{fam}/{r['name']}"
            cases.append(make_real_case(cid, obs, copy_vecs, names))
    return cases


# --------------------------------------------------------------------- synthetic
def collect_synth():
    cases = []

    # 1. two copies identical over the spanned cols -> n_decisive 0, margin 0 (TIED)
    cv = {"copy0": {0: "A", 1: "C", 2: "G"}, "copy1": {0: "A", 1: "C", 2: "G"}}
    cases.append(make_synth_case("synth/tie_identical",
                                 {0: "A", 1: "C", 2: "G"}, cv, ["copy0", "copy1"]))

    # 2. two copies disagree but the obs matches NEITHER -> equal logL, margin 0, decisive
    cv = {"copy0": {0: "A"}, "copy1": {0: "C"}}
    cases.append(make_synth_case("synth/tie_both_mismatch", {0: "G"}, cv, ["copy0", "copy1"]))

    # 3. clear winner, 2 copies
    cv = {"copy0": {0: "A", 1: "A", 2: "A"}, "copy1": {0: "C", 1: "C", 2: "C"}}
    cases.append(make_synth_case("synth/clear2", {0: "A", 1: "A", 2: "A"}, cv, ["copy0", "copy1"]))

    # 4. single copy spanned -> margin +inf
    cases.append(make_synth_case("synth/single_spanned",
                                 {0: "A", 1: "C"}, {"copy0": {0: "A", 1: "C"}}, ["copy0"]))

    # 5. single copy, obs empty -> n_span 0, margin +inf, logL {copy0:0}
    cases.append(make_synth_case("synth/single_empty", {}, {"copy0": {0: "A"}}, ["copy0"]))

    # 6. obs WITH None interleaved -> None cols filtered from spanned
    cv = {"copy0": {0: "A", 1: "C", 2: "G", 3: "T"}, "copy1": {0: "T", 1: "C", 2: "C", 3: "A"}}
    cases.append(make_synth_case("synth/obs_with_none",
                                 {0: "A", 1: None, 2: "G", 3: None}, cv, ["copy0", "copy1"]))

    # 7. ALL obs None -> spanned empty, n_span 0, all logL 0, margin 0 (2 copies, tie)
    cases.append(make_synth_case("synth/all_none",
                                 {0: None, 1: None}, cv, ["copy0", "copy1"]))

    # 8. >2 copy family (5 copies), distinct winner + runner-up
    names5 = [f"copy{i}" for i in range(5)]
    cv = {f"copy{i}": {0: "ACGTA"[i], 1: "CGTAC"[i]} for i in range(5)}
    cases.append(make_synth_case("synth/five_copy", {0: "A", 1: "C"}, cv, names5))

    # 9. missing copy columns: some copies lack a base at a spanned col (get -> None -> skipped)
    cv = {"copy0": {0: "A", 1: "C"}, "copy1": {0: "T"}, "copy2": {1: "G"}}
    cases.append(make_synth_case("synth/sparse_cols",
                                 {0: "A", 1: "C"}, cv, ["copy0", "copy1", "copy2"]))

    # 10. 3-copy tie among two leaders (copy0==copy1 best, copy2 worse) -> best=copy0, margin 0
    cv = {"copy0": {0: "A"}, "copy1": {0: "A"}, "copy2": {0: "C"}}
    cases.append(make_synth_case("synth/three_leadtie", {0: "A"}, cv, ["copy0", "copy1", "copy2"]))

    # ---- junction-term cases (>=20). obs empty or small; exercise read_bounds/copy_bounds ----
    names2 = ["copy0", "copy1"]
    names3 = ["copy0", "copy1", "copy2"]
    j = 0

    def jcase(tag, obs, cv, names, rb, cb, **kw):
        nonlocal j
        j += 1
        cases.append(make_synth_case(f"junc/{j:02d}_{tag}", obs, cv, names, read_bounds=rb, copy_bounds=cb, **kw))

    empty_cv2 = {"copy0": {}, "copy1": {}}
    empty_cv3 = {"copy0": {}, "copy1": {}, "copy2": {}}

    # both copies MATCH every read boundary -> +jw each, no decisive
    jcase("both_match", {}, empty_cv2, names2, {100, 200},
          [["copy0", {100, 200}], ["copy1", {100, 200}]])
    # both copies MISS every boundary -> -jw each, no decisive
    jcase("both_miss", {}, empty_cv2, names2, {100, 200},
          [["copy0", {900}], ["copy1", {900}]])
    # MIXED: copy0 has boundary, copy1 lacks it -> decisive, split logL
    jcase("mixed_one", {}, empty_cv2, names2, {100},
          [["copy0", {100}], ["copy1", {500}]])
    # copy1 entirely absent from copy_bounds (get -> ()) -> -jw for copy1
    jcase("copy_absent", {}, empty_cv2, names2, {100}, [["copy0", {100}]])
    # btol boundary: |b - x| == btol (==4) -> HAS true (within tolerance)
    jcase("btol_edge_in", {}, empty_cv2, names2, {104}, [["copy0", {100}], ["copy1", {100}]])
    # btol boundary: |b - x| == btol+1 (==5) -> HAS false (outside tolerance)
    jcase("btol_edge_out", {}, empty_cv2, names2, {105}, [["copy0", {100}], ["copy1", {100}]])
    # exact match |b-x|==0
    jcase("btol_exact", {}, empty_cv2, names2, {100}, [["copy0", {100}], ["copy1", {100}]])
    # negative-side tolerance |b-x|==4 with b<x
    jcase("btol_neg_in", {}, empty_cv2, names2, {96}, [["copy0", {100}], ["copy1", {100}]])
    # multiple read boundaries, each mixed differently
    jcase("multi_mixed", {}, empty_cv3, names3, {100, 200, 300},
          [["copy0", {100, 200, 300}], ["copy1", {200}], ["copy2", {300}]])
    # decisive boundary + non-decisive boundary combined
    jcase("mixed_and_shared", {}, empty_cv2, names2, {100, 200},
          [["copy0", {100, 200}], ["copy1", {200}]])
    # PSV term AND junction term combined (both contribute to logL and n_decisive)
    cvp = {"copy0": {0: "A"}, "copy1": {0: "C"}}
    jcase("psv_plus_junc", {0: "A"}, cvp, names2, {100},
          [["copy0", {100}], ["copy1", {700}]])
    # three copies, boundary present in exactly one -> decisive
    jcase("three_one_present", {}, empty_cv3, names3, {100},
          [["copy0", {100}], ["copy1", {}], ["copy2", {}]])
    # three copies, boundary present in exactly two -> decisive (0<2<3)
    jcase("three_two_present", {}, empty_cv3, names3, {100},
          [["copy0", {100}], ["copy1", {100}], ["copy2", {}]])
    # three copies, boundary present in ALL -> NOT decisive (present==n)
    jcase("three_all_present", {}, empty_cv3, names3, {100},
          [["copy0", {100}], ["copy1", {100}], ["copy2", {100}]])
    # boundary present in NONE -> NOT decisive (present==0)
    jcase("three_none_present", {}, empty_cv3, names3, {100},
          [["copy0", {}], ["copy1", {}], ["copy2", {}]])
    # several boundaries, many additions (float accumulation order stress)
    jcase("many_bounds", {}, empty_cv2, names2, {10, 20, 30, 40, 50, 60, 70},
          [["copy0", {10, 30, 50, 70}], ["copy1", {20, 40, 60}]])
    # custom jw
    jcase("custom_jw", {}, empty_cv2, names2, {100, 200},
          [["copy0", {100}], ["copy1", {200}]], jw=3.5)
    # custom btol widening a near miss into a hit (|b-x|=8 <= btol 10)
    jcase("custom_btol", {}, empty_cv2, names2, {108}, [["copy0", {100}], ["copy1", {900}]], btol=10)
    # copy has MULTIPLE candidate boundaries; any-within-tol
    jcase("multi_copybounds", {}, empty_cv2, names2, {150},
          [["copy0", {10, 148, 900}], ["copy1", {10, 900}]])
    # read_bounds present but copy_bounds EMPTY dict -> junction skipped (Python falsy)
    jcase("copybounds_empty", {0: "A"}, cvp, names2, {100}, [])
    # read_bounds EMPTY -> junction skipped even though copy_bounds present
    jcase("readbounds_empty", {0: "A"}, cvp, names2, set(),
          [["copy0", {100}], ["copy1", {100}]])
    # combined PSV winner + junction reinforcing the SAME copy
    jcase("psv_junc_reinforce", {0: "A", 1: "A"},
          {"copy0": {0: "A", 1: "A"}, "copy1": {0: "C", 1: "C"}}, names2, {100},
          [["copy0", {100}], ["copy1", {800}]])
    # combined PSV winner but junction CONTRADICTS (argues for the other copy)
    jcase("psv_junc_conflict", {0: "A", 1: "A"},
          {"copy0": {0: "A", 1: "A"}, "copy1": {0: "C", 1: "C"}}, names2, {100},
          [["copy0", {800}], ["copy1", {100}]])

    return cases


def main():
    real = collect_real()
    synth = collect_synth()
    n_junction = sum(1 for c in synth if c["read_bounds"] is not None)
    assert len(real) >= 200, f"only {len(real)} real reads (need >=200)"
    assert n_junction >= 20, f"only {n_junction} junction cases (need >=20)"
    cases = real + synth

    fixture = dict(
        note=("Golden byte-parity fixture for src/rustle/vg_family/o2_margin_gate.rs "
              "(port of bench/copy_assign.py::assign_read). Floats stored as raw IEEE-754 "
              "u64 bits. Generated by bench/gen_o2_margin_gate_fixture.py (PYTHONHASHSEED=0)."),
        constants=dict(
            ERR_bits=d2bits(ca.ERR), MARGIN_bits=d2bits(ca.MARGIN),
            JW_bits=d2bits(5.0), BTOL=4,
        ),
        n_real=len(real), n_synthetic=len(synth), n_junction=n_junction,
        cases=cases,
    )
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(fixture, fh, indent=0)
    print(f"wrote {OUT}")
    print(f"  cases: {len(cases)}  (real={len(real)}, synthetic={len(synth)}, junction={n_junction})")
    # quick status breakdown of real cases for visibility
    import collections
    dec = collections.Counter()
    for c in real:
        nd = c["expect"]["n_decisive"]
        dec["decisive>=1" if nd >= 1 else "tied(nd=0)"] += 1
    print(f"  real n_decisive: {dict(dec)}")


if __name__ == "__main__":
    main()
