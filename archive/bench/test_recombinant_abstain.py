#!/usr/bin/env python3
"""Tests for the recombinant-read ABSTAIN leg (bench/recombinant_abstain.py) wired into the O2 gate.

Requirements verified:
  (1) is_recombinant flags a bubble-bounded copy-switch read and NOT a clean single-copy read;
  (2) apply_abstain_to_vg is a PURE MONOTONE ADDITION OF ABSTENTIONS -- non-recombinant reads and
      already-abstaining reads are byte-identical; only force-assigned recombinant reads flip
      assigned -> 'recombinant';
  (3) DEFAULT-ON; the RUSTLE_NO_RECOMBINANT_ABSTAIN=1 opt-out makes it a no-op (recovers prior);
  (4) DETERMINISTIC (two runs identical);
  (5) DATA-BACKED (skipped if the o2vg cache is absent): over the 119 cached co-located families the
      leg moves exactly the 674 localized_tract reads (incl GSTM2 fam13) + 1357 systematic_copy_split
      + 183 sporadic_chimera = 2214 reads assigned->abstain, with ZERO byte-identity violations on
      non-recombinant reads, matching bench/vg_read_path_recombination.tsv.

Deterministic (PYTHONHASHSEED=0; no RNG). Run:
  PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/test_recombinant_abstain.py
  (or via pytest: PYTHONHASHSEED=0 pytest bench/test_recombinant_abstain.py)
"""
import copy
import csv
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import recombinant_abstain as ra   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


# --------------------------------------------------------------------------- synthetic fixture
def _synthetic_vg():
    """Two copies over 4 PSV bubbles: copy0=AAAA, copy1=CCCC. Reads:
       * 4x clean copy0 (AAAA) + 4x clean copy1 (CCCC)  -> build read-support for both alleles;
       * 1 crossover read AACC (copy0 arm | copy1 arm)  -> recombinant;
       * 1 clean copy0 read AAAA that the gate ASSIGNED -> must NOT flip;
       * 1 crossover read AACC that the gate left AMBIGUOUS -> must NOT flip (only add abstentions)."""
    ncol = 4
    bubbles = [dict(hap={"copy0": "A", "copy1": "C"}) for _ in range(ncol)]
    copies = [{"name": "copy0"}, {"name": "copy1"}]

    def read(name, seq, status, best):
        return dict(name=name, status=status, best_copy=best,
                    obs={str(j): seq[j] for j in range(ncol)})

    reads = []
    for i in range(4):
        reads.append(read(f"supA{i}", "AAAA", "assigned", "copy0"))
        reads.append(read(f"supC{i}", "CCCC", "assigned", "copy1"))
    reads.append(read("recomb_assigned", "AACC", "assigned", "copy0"))  # recombinant + force-assigned
    reads.append(read("clean_assigned", "AAAA", "assigned", "copy0"))   # single-copy, assigned
    reads.append(read("recomb_ambiguous", "AACC", "ambiguous", None))   # recombinant but already abstains
    return dict(family=999, copies=copies, bubbles=bubbles, n_bubbles=ncol, reads=reads)


# --------------------------------------------------------------------------- (1) predicate
def test_is_recombinant_predicate():
    os.environ.pop(ra.OPTOUT_ENV, None)
    copy_vecs = {"copy0": {j: "A" for j in range(4)}, "copy1": {j: "C" for j in range(4)}}
    supported = {j: {"A", "C"} for j in range(4)}
    names = ["copy0", "copy1"]
    rec, kind, arms = ra.is_recombinant({0: "A", 1: "A", 2: "C", 3: "C"}, supported, copy_vecs, names)
    assert rec and kind == "crossover", (rec, kind)
    assert set((arms["A"], arms["B"])) == {"copy0", "copy1"}
    rec2, _, _ = ra.is_recombinant({0: "A", 1: "A", 2: "A", 3: "A"}, supported, copy_vecs, names)
    assert rec2 is False, "a clean single-copy read is NOT recombinant"
    print("(1) predicate: crossover flagged, clean single-copy not flagged  OK")


# --------------------------------------------------------------------------- (2) monotone abstain
def test_apply_abstain_monotone():
    os.environ.pop(ra.OPTOUT_ENV, None)
    vg = _synthetic_vg()
    before = copy.deepcopy(vg["reads"])
    flipped = ra.apply_abstain_to_vg(vg)
    names_flipped = {f["read"] for f in flipped}
    assert names_flipped == {"recomb_assigned"}, names_flipped
    after = {r["name"]: r for r in vg["reads"]}
    for b in before:
        a = after[b["name"]]
        if b["name"] == "recomb_assigned":
            assert b["status"] == "assigned" and a["status"] == ra.ABSTAIN_STATUS
            assert a["best_copy"] is None
            assert a["gate_status_pre_abstain"] == "assigned" and a["recombinant_kind"] == "crossover"
        else:
            # byte-identical: non-recombinant reads AND already-abstaining recombinant reads untouched
            assert a == b, f"read {b['name']} changed but must be byte-identical: {b} -> {a}"
    print("(2) monotone: only 'recomb_assigned' flipped; all others byte-identical  OK")


# --------------------------------------------------------------------------- (3) opt-out no-op
def test_optout_recovers_prior():
    os.environ[ra.OPTOUT_ENV] = "1"
    vg = _synthetic_vg()
    before = copy.deepcopy(vg["reads"])
    flipped = ra.apply_abstain_to_vg(vg)
    os.environ.pop(ra.OPTOUT_ENV, None)
    assert flipped == [], "opt-out must be a no-op"
    assert vg["reads"] == before, "opt-out must leave every read byte-identical"
    print("(3) opt-out (RUSTLE_NO_RECOMBINANT_ABSTAIN=1): no-op, byte-identical  OK")


# --------------------------------------------------------------------------- (4) determinism
def test_determinism():
    os.environ.pop(ra.OPTOUT_ENV, None)
    vg1 = _synthetic_vg(); ra.apply_abstain_to_vg(vg1)
    vg2 = _synthetic_vg(); ra.apply_abstain_to_vg(vg2)
    assert vg1["reads"] == vg2["reads"], "two runs must be identical"
    print("(4) determinism: two runs identical  OK")


# --------------------------------------------------------------------------- (5) data-backed
def test_data_backed_674_move():
    """Over the cached co-located families: 674 localized_tract (incl GSTM2 fam13) + 1357 systematic +
    183 sporadic = 2214 move assigned->abstain; 0 byte-identity violations on non-recombinant reads."""
    import o2_vg_visualization as o2vg
    import psv_graph_genomewide as pg
    tsv_path = os.path.join(HERE, "vg_read_path_recombination.tsv")
    if not os.path.isdir(o2vg.CACHE) or not os.path.exists(tsv_path) or not os.path.exists(pg.GENOME):
        print("(5) data-backed: SKIPPED (o2vg cache / detector tsv / genome not present)")
        return
    fams = o2vg.load_families()
    tsv = list(csv.DictReader(open(tsv_path), delimiter="\t"))
    struct_of = {(r["family"], r["read"]): r["structural_call"] for r in tsv}
    target = collections.Counter(r["structural_call"] for r in tsv if r["gate_status"] == "assigned")

    changed = []
    violations = 0
    for fid in sorted(fams):
        if not os.path.exists(os.path.join(o2vg.CACHE, f"fam{fid}.json")):
            continue
        if len(pg.dedup_copies(fams[fid])) < 2:
            continue
        os.environ[ra.OPTOUT_ENV] = "1"
        off = {r["name"]: r for r in o2vg.materialize_family(fid, use_cache=True)["reads"]}
        os.environ.pop(ra.OPTOUT_ENV, None)
        vg_on = o2vg.materialize_family(fid, use_cache=True)
        if vg_on["n_bubbles"] < 2:
            continue
        for r in vg_on["reads"]:
            b = off[r["name"]]
            if r == b:
                continue
            if r["status"] == "assigned" or b["status"] != "assigned" \
                    or r["status"] != ra.ABSTAIN_STATUS or r["best_copy"] is not None:
                violations += 1
            changed.append((str(fid), r["name"]))
    by_class = collections.Counter(struct_of.get(k, "?") for k in changed)
    assert violations == 0, f"byte-identity violations on non-recombinant reads: {violations}"
    assert by_class.get("localized_tract") == target["localized_tract"] == 674, by_class
    assert by_class.get("systematic_copy_split") == target["systematic_copy_split"], by_class
    assert by_class.get("sporadic_chimera") == target["sporadic_chimera"], by_class
    assert len(changed) == 2214, len(changed)
    fam13 = sum(1 for f, _ in changed if f == "13")
    assert fam13 >= 292, f"GSTM2 fam13 flipped {fam13} (expected >=292)"
    print(f"(5) data-backed: {len(changed)} reads moved assigned->abstain "
          f"(localized_tract={by_class['localized_tract']} incl GSTM2 fam13={fam13}, "
          f"systematic={by_class['systematic_copy_split']}, sporadic={by_class['sporadic_chimera']}); "
          f"0 byte-identity violations  OK")


if __name__ == "__main__":
    os.environ.setdefault("PYTHONHASHSEED", "0")
    test_is_recombinant_predicate()
    test_apply_abstain_monotone()
    test_optout_recovers_prior()
    test_determinism()
    test_data_backed_674_move()
    print("\nALL TESTS PASSED")
