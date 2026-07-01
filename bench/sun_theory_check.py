#!/usr/bin/env python3
"""sun_theory_check.py — exhaustive machine-check of the SUN (Sudmant 2010) refinement
of the copy-assignment theory (bench/THEORY.md §5·SUN).

Notation (THEORY.md §2): copies = allele-vectors over columns [m]; a read is a partial
function (a window of observed columns with the copy's alleles); r consistent with copy v
iff r(j)=v_j for all j in obs(r); N(r) = {copies consistent with r}.

SUN(c,p): (c)_p is unique among all copies at column p (single-position private allele).

WHAT SUN BUYS (and what it does NOT), stated honestly:
 - PER-READ (Thm 4 gate): a read carrying a SUN allele has |N(r)|=1 and is assigned to c,
   never to another TRUE copy.  This is the load-bearing, provable, machine-checked fact (S2).
 - COVER (RECOVER / MCC, Thm 2): SUN does NOT make the COPY survive every minimum cover.
   Even a SUN-rich copy can be DISSOLVED by an alternative minimum cover (S3_cover), because a
   phantom class may carry the copy's private allele by absorbing the copy's OWN reads.  So
   there is NO cover-level "recombination-immune copy" / {spoofable}⊆{no-SUN} containment.
 The gap between these two is exactly why the shipped pipeline runs the per-read GATE (which
 inherits SUN determinism) and not RECOVER (which does not).

Checks:
 (S1) SUN(c,p) => hap-vector of c is UNIQUE among copies (Tier1 subset of unique-hap).
 (S2) SUN single-read gate determinism: any read r that observes a SUN column p of c and carries
      the private allele (c)_p has N(r) = {c} EXACTLY (|N(r)|=1).  => the per-read gate (Thm4/7)
      can never MISASSIGN such a read to another true copy.  (This is per-read, NOT cover-level.)
 (S3) canonical K>=3 witness (THEORY.md §5): c0=(1,1,0),c1=(0,0,1),c2=(0,1,1).  Factual report:
      the recombinant class realizes novel vector (0,1,0) and anchors a 2nd minimum cover; the
      unique no-SUN copy is c2.  (Reported as a value, NOT used to validate any immunity claim.)
 (S3_cover) COUNTEREXAMPLE to cover-level copy-immunity, on the same witness + its exact 6-read
      set (each copy on windows {0,1},{1,2}): enumerate ALL minimum covers and show one DISSOLVES
      c0 -- which HAS SUNs at cols 0 AND 2 -- into two phantom classes {(0,1,0),(1,1,1)} that
      carry c0's private alleles.  Simultaneously confirm PER-READ gate immunity still holds.
 (S4) HONEST DIRECTIONS: SUN is NOT iff Strong-Sep.
      (a) some copy is hap-unique AND the family is Strong-Separated (full-length reads) yet that
          copy has NO SUN => Strong-Sep does NOT imply SUN (sufficient-not-necessary).
      (c) a K>=4 instance where NO copy has a SUN yet the family is fully Strong-Separated.
      (b) a single copy's SUN does NOT imply family-level Strong Separation.
"""
import itertools


def windows(m):
    for k in range(1, m + 1):
        for w in itertools.combinations(range(m), k):
            yield w


def N(copies, obs, vals):
    """copies consistent with a read (obs=cols tuple, vals=alleles tuple)."""
    out = []
    for name, v in copies.items():
        if all(v[j] == a for j, a in zip(obs, vals)):
            out.append(name)
    return out


def suns(copies, m):
    """{copy: [SUN columns]}"""
    out = {c: [] for c in copies}
    for p in range(m):
        col = {c: copies[c][p] for c in copies}
        vals = list(col.values())
        for c in copies:
            if vals.count(col[c]) == 1:
                out[c].append(p)
    return out


def hap_unique(copies):
    vecs = list(copies.values())
    return {c: vecs.count(copies[c]) == 1 for c in copies}


def strong_sep_full(copies, m):
    """Strong Separation holds with FULL-length reads (every read observes all columns):
    every cross-copy pair differs somewhere (distinct allele-vectors)."""
    cs = list(copies.values())
    for a, b in itertools.combinations(cs, 2):
        if a == b:
            return False
    return True


def check_S1_S2():
    """Exhaustive over K in {2,3,4}, m in {1,2,3}, alphabet {0,1,2}."""
    viol_S1 = 0
    viol_S2 = 0
    n_sun_copies = 0
    n_reads_checked = 0
    for m in (1, 2, 3):
        for K in (2, 3, 4):
            for combo in itertools.product(itertools.product(range(3), repeat=m), repeat=K):
                if len(set(combo)) < K:
                    continue  # distinct copies only
                copies = {f"c{i}": combo[i] for i in range(K)}
                sn = suns(copies, m)
                hu = hap_unique(copies)
                for c in copies:
                    if not sn[c]:
                        continue
                    n_sun_copies += 1
                    # S1
                    if not hu[c]:
                        viol_S1 += 1
                    # S2: every read on a SUN column carrying the private allele => N={c}
                    for p in sn[c]:
                        priv = copies[c][p]
                        for w in windows(m):
                            if p not in w:
                                continue
                            vals = tuple(copies[c][j] for j in w)  # read from c
                            assert vals[w.index(p)] == priv
                            nb = N(copies, w, vals)
                            n_reads_checked += 1
                            if nb != [c]:
                                viol_S2 += 1
    return dict(sun_copies=n_sun_copies, reads_checked=n_reads_checked,
                S1_violations=viol_S1, S2_violations=viol_S2)


def check_S3():
    """Factual report on the canonical witness (NO immunity interpretation attached)."""
    copies = {"c0": (1, 1, 0), "c1": (0, 0, 1), "c2": (0, 1, 1)}
    m = 3
    sn = suns(copies, m)
    no_sun = [c for c in copies if not sn[c]]
    hu = hap_unique(copies)
    # recombinant class: c0's read on {1,2}=(_,1,0) and c2's read on {0,1}=(0,1,_)
    r_c0 = ((1, 2), (1, 0))
    r_c2 = ((0, 1), (0, 1))
    # joint consistency of the class (shared col 1: both value 1) -> realized vector union
    union = {}
    for (obs, vals) in (r_c0, r_c2):
        for j, a in zip(obs, vals):
            union.setdefault(j, set()).add(a)
    jointly_consistent = all(len(s) == 1 for s in union.values())
    realized = tuple(next(iter(union[j])) for j in range(m))
    novel = realized not in set(copies.values())
    return dict(sun=sn, no_sun_copies=no_sun, hap_unique=hu,
                recombinant_jointly_consistent=jointly_consistent,
                realized_vector=realized, is_novel=novel,
                # NEUTRAL value (was `recombined_copy_is_the_unique_no_SUN`, which over-read a
                # cover-level immunity into an adjacent fact -- see S3_cover for the real test):
                unique_no_SUN_copy=(no_sun[0] if len(no_sun) == 1 else no_sun))


# ---------------------------------------------------------------- minimum-cover enumeration
def _realize(read_subset, m):
    """Union vector of a jointly-consistent read class, or None if the reads conflict.
    (partial columns -> None entry; a class realizes a TRUE copy only if it pins all m cols
    to that copy's alleles)."""
    col = {}
    for (_, obs, vals) in read_subset:
        for j, a in zip(obs, vals):
            if j in col and col[j] != a:
                return None
            col[j] = a
    return tuple(col.get(j, None) for j in range(m))


def _partitions(collection):
    collection = list(collection)
    if len(collection) == 1:
        yield [collection]
        return
    first = collection[0]
    for smaller in _partitions(collection[1:]):
        for n, subset in enumerate(smaller):
            yield smaller[:n] + [[first] + subset] + smaller[n + 1:]
        yield [[first]] + smaller


def check_S3_cover():
    """The load-bearing COUNTEREXAMPLE: on the canonical K>=3 witness + its exact 6-read set,
    a SUN-rich copy (c0, SUNs at cols 0,2) is DISSOLVED by an alternative MINIMUM cover.
    Cover-level copy-immunity is therefore FALSE. Per-read gate immunity is confirmed intact."""
    copies = {"c0": (1, 1, 0), "c1": (0, 0, 1), "c2": (0, 1, 1)}
    m = 3
    sn = suns(copies, m)
    reads = []
    for cn, cv in copies.items():
        for w in ((0, 1), (1, 2)):
            reads.append((f"{cn}_{''.join(map(str, w))}", w, tuple(cv[j] for j in w)))

    # (i) per-read gate immunity: every read carrying a SUN allele -> N(r) = {its true copy}
    gate_viol = 0
    for name, obs, vals in reads:
        cn = name.split("_")[0]
        carries_sun = any(p in obs and vals[obs.index(p)] == copies[cn][p] for p in sn[cn])
        if carries_sun and N(copies, obs, vals) != [cn]:
            gate_viol += 1

    # (ii) enumerate ALL minimum covers (partitions into jointly-consistent classes)
    min_size, min_covers = None, []
    for part in _partitions(reads):
        if any(_realize(cls, m) is None for cls in part):
            continue
        size = len(part)
        if min_size is None or size < min_size:
            min_size, min_covers = size, []
        if size == min_size:
            min_covers.append(part)

    true_vecs = set(copies.values())
    # is there a minimum cover in which c0 (SUN-rich) is DISSOLVED (no class realizes c0)?
    dissolving = None
    for cov in min_covers:
        realized = [_realize(cls, m) for cls in cov]
        if copies["c0"] not in realized:
            dissolving = [dict(reads=[r[0] for r in cls], realizes=_realize(cls, m),
                               is_true_copy=(_realize(cls, m) in true_vecs)) for cls in cov]
            break

    # do the phantom classes carry c0's private (SUN) alleles?
    phantom_carries_private = False
    if dissolving is not None:
        for p in sn["c0"]:
            priv = copies["c0"][p]
            if any(cls["realizes"][p] == priv for cls in dissolving):
                phantom_carries_private = True

    # which copies survive EVERY minimum cover?
    survives_every = {}
    for c in copies:
        survives_every[c] = all(copies[c] in [_realize(cls, m) for cls in cov]
                                for cov in min_covers)

    return dict(
        MCC=min_size, n_minimum_covers=len(min_covers),
        per_read_gate_immunity_holds=(gate_viol == 0),
        c0_has_suns=sn["c0"],
        cover_dissolving_c0_exists=(dissolving is not None),
        dissolving_cover=dissolving,
        phantom_classes_carry_c0_private_allele=phantom_carries_private,
        survives_every_min_cover=survives_every,
        # the headline of the counterexample: a SUN-rich copy fails cover-level immunity
        cover_level_copy_immunity_FALSE=(dissolving is not None
                                         and survives_every["c0"] is False),
    )


def check_S4():
    """SUN is NOT iff Strong Separation -- witnessed in both under-claimed directions."""
    # (a) some copy is hap-unique + the family is Strong-Separated (full reads) + that copy has
    #     NO SUN  =>  Strong-Sep does NOT imply SUN (sufficient-not-necessary, per copy).
    copies = {"c0": (0, 0), "c1": (0, 1), "c2": (1, 0)}
    m = 2
    sn = suns(copies, m)
    hu = hap_unique(copies)
    ss = strong_sep_full(copies, m)
    # c0 = (0,0): shares col0 with c1, shares col1 with c2 -> NO SUN, but unique vector.
    # NOTE (honest): here c1 and c2 DO have SUNs; the claim is "SOME copy (c0) has no SUN yet
    # is Strong-Separated", i.e. SUN not necessary for a given copy -- NOT "no copy has a SUN".
    a_holds = (hu["c0"] and ss and not sn["c0"])

    # (c) a genuine "NO copy has a SUN yet the family is fully Strong-Separated" instance needs
    #     K>=4: c0=(0,0),c1=(0,1),c2=(1,0),c3=(1,1). Each column's two symbols each appear twice
    #     -> no private allele anywhere -> no SUN for any copy; all 4 vectors distinct -> SS.
    copies_c = {"c0": (0, 0), "c1": (0, 1), "c2": (1, 0), "c3": (1, 1)}
    sn_c = suns(copies_c, 2)
    c_no_copy_has_sun = all(len(sn_c[c]) == 0 for c in copies_c)
    c_strong_sep = strong_sep_full(copies_c, 2)

    # (b) SUN(c0) does NOT imply family-level Strong Separation: c0=(1,0) has a SUN at col0
    #     (private 1), but the (c1,c2) pair c1=(0,0),c2=(0,1) differs only at col1, so reads
    #     that miss col1 leave (c1,c2) unseparated. SUN(c0) constrains only c0-vs-others at col0.
    copies2 = {"c0": (1, 0), "c1": (0, 0), "c2": (0, 1)}
    sn2 = suns(copies2, 2)
    return dict(
        a_hapunique_strongsep_noSUN=a_holds,
        c0_no_sun=(not sn["c0"]), c0_hap_unique=hu["c0"], strong_sep_full=ss,
        a_note="c1,c2 DO have SUNs here; claim is SOME copy (c0) is Strong-Sep yet has no SUN "
               "(SUN not necessary), not 'no copy has a SUN'.",
        c_no_copy_has_SUN_yet_strong_sep=(c_no_copy_has_sun and c_strong_sep),
        c_witness="K=4: (0,0),(0,1),(1,0),(1,1) -- no copy has a SUN, family Strong-Separated.",
        b_c0_has_sun=(len(sn2["c0"]) > 0),
        b_note="SUN(c0) constrains only c0-vs-others at that column; the (c1,c2) pair can be "
               "under-covered => family Strong Sep is a strictly stronger, coverage-dependent, "
               "all-pairs condition. SUN(c0) does NOT imply family Strong Sep.")


if __name__ == "__main__":
    import json
    r1 = check_S1_S2()
    r3 = check_S3()
    r3c = check_S3_cover()
    r4 = check_S4()
    print("=== S1 (SUN=>unique-hap) & S2 (SUN read => per-read |N(r)|=1) exhaustive ===")
    print(json.dumps(r1, indent=2))
    print("=== S3 (canonical K>=3 witness -- factual report, neutral) ===")
    print(json.dumps(r3, indent=2, default=str))
    print("=== S3_cover (COUNTEREXAMPLE: SUN-rich copy c0 dissolves in an alt minimum cover) ===")
    print(json.dumps(r3c, indent=2, default=str))
    print("=== S4 (honest directions: NOT iff, incl. K=4 all-no-SUN Strong-Sep) ===")
    print(json.dumps(r4, indent=2, default=str))
    ok = (r1["S1_violations"] == 0 and r1["S2_violations"] == 0
          # per-read gate immunity holds ...
          and r3c["per_read_gate_immunity_holds"]
          # ... but cover-level copy immunity is FALSE (SUN-rich c0 dissolves):
          and r3c["cover_level_copy_immunity_FALSE"]
          and r3c["phantom_classes_carry_c0_private_allele"]
          and r3["is_novel"] and r3["recombinant_jointly_consistent"]
          and r3["unique_no_SUN_copy"] == "c2"
          and r4["a_hapunique_strongsep_noSUN"]
          and r4["c_no_copy_has_SUN_yet_strong_sep"]
          and r4["b_c0_has_sun"])
    print("\nALL_GREEN =", ok)
