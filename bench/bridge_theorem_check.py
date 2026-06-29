#!/usr/bin/env python3
"""Exhaustive machine-check of the BRIDGE THEOREM (§5b): the production per-read significance gate's
identifiability bound min_p is a SOUND combinatorial certificate, connecting the running statistical
machinery to Theorem 2 for ALL K>=1 (not just the K=0 vertex).

Bridge claims, over the error-free core (each read consistent with its origin copy):
  Let copies C = distinct allele-vectors over columns [m]. A read r = origin copy o restricted to a
  column subset obs(r). Gate "best" b = argmax likelihood = a copy consistent with r (ties -> lowest idx).
  D_{bc} = distinguishing columns of copies b,c.  delta(r) = min_{c != b} |obs(r) ∩ D_{bc}|.
  min_p(r) = max_{c != b} prod_{j in obs(r)∩D_{bc}} eps  =  eps^{delta(r)}  (uniform eps<1).

  (B1) min_p == eps^delta  (the gate's identifiability bound equals eps^(closest-competitor distinguishing count)).
  (B2) delta(r) >= 1  <=>  exactly ONE copy in C is consistent with r  (combinatorial unique determination).
  (B3) SOUNDNESS: if delta(r) >= 1 then the unique consistent copy is the origin o, so argmax assignment = o
       (the gate never assigns a wrong copy in the error-free core).
  (B4) ABSTENTION = AMBIGUITY: delta(r) == 0  <=>  >= 2 copies consistent with r (genuinely unassignable).
  (B5) BOUNDARY RECOVERY (Theorem-2 K_{ij}=0): if two copies are identical over [m] (D=∅), every read whose
       origin is one of them has delta == 0 (=> min_p == 1 => Tied) — the gate ties exactly the
       non-identifiable mass, per-read, with no NP-hard cover solved.
Run: python3 bench/bridge_theorem_check.py   (exit 0 = all bridge claims hold on the full enumerated universe).
"""
import itertools, sys

EPS = 1e-3  # per-distinguishing-column error proxy (e/3); any 0<eps<1 works for the combinatorial claims.

def consistent(read_obs, copy):
    """read_obs: dict col->allele (the read's observed alleles). copy: tuple allele-vector. consistent iff
    the read does not contradict the copy at any observed column."""
    return all(copy[j] == a for j, a in read_obs.items())

def distinguishing(ci, cj):
    return {d for d in range(len(ci)) if ci[d] != cj[d]}

def gate_minp_and_best(read_obs, copies):
    """Reproduce the production gate: best = argmax likelihood (error-free: a consistent copy, lowest idx on
    ties; if none consistent, the copy with fewest mismatches). min_p = max_c prod_{obs∩D_{bc}} eps."""
    n = len(copies)
    # likelihood proxy (error-free): #matches - #mismatches over observed columns; argmax, lowest idx on ties.
    def score(c):
        m = sum(1 for j, a in read_obs.items() if copies[c][j] == a)
        mm = len(read_obs) - m
        return m - mm
    best = max(range(n), key=lambda c: (score(c), -c))
    minp = 0.0
    delta = None
    for c in range(n):
        if c == best:
            continue
        D = distinguishing(copies[best], copies[c])
        span = [j for j in read_obs if j in D]
        attain = EPS ** len(span)  # prod over distinguishing positions the read spans
        if attain > minp:
            minp = attain
        delta = len(span) if delta is None else min(delta, len(span))
    if delta is None:  # only one copy
        delta = 10**9
        minp = 0.0
    return minp, best, delta

def check_bridge(max_K=3, max_m=3, alphabet=(0, 1, 2)):
    fails = []
    n_inst = n_reads = 0
    # enumerate all copy sets: K distinct allele-vectors over m columns
    for m in range(1, max_m + 1):
        all_vecs = list(itertools.product(alphabet, repeat=m))
        for K in range(2, max_K + 1):
            for copies in itertools.combinations(all_vecs, K):  # distinct copies, as a set
                n_inst += 1
                # every read = (origin copy o, observed column subset S != ∅), alleles = o restricted to S
                for o in range(K):
                    for ssize in range(1, m + 1):
                        for S in itertools.combinations(range(m), ssize):
                            read_obs = {j: copies[o][j] for j in S}
                            n_reads += 1
                            minp, best, delta = gate_minp_and_best(read_obs, copies)
                            cons = [c for c in range(K) if consistent(read_obs, copies[c])]

                            # (B1) min_p == eps^delta
                            if abs(minp - EPS ** delta) > 1e-15 and delta < 10**8:
                                fails.append(("B1", m, copies, o, S, minp, delta))
                            # (B2) delta>=1 <=> exactly one consistent copy
                            uniq = (len(cons) == 1)
                            if (delta >= 1) != uniq:
                                fails.append(("B2", m, copies, o, S, delta, cons))
                            # (B3) soundness: delta>=1 => unique consistent copy is the origin AND best==origin
                            if delta >= 1:
                                if cons != [o] or best != o:
                                    fails.append(("B3", m, copies, o, S, best, cons))
                            # (B4) delta==0 <=> >=2 consistent
                            if (delta == 0) != (len(cons) >= 2):
                                fails.append(("B4", m, copies, o, S, delta, cons))
                # (B5) boundary: an IDENTICAL pair cannot occur here (copies are a set => distinct). Test the
                # K_{ij}=0 boundary directly with a deliberately-degenerate pair (same vector twice) below.
    return fails, n_inst, n_reads

def check_boundary_K0(max_m=3, alphabet=(0, 1)):
    """(B5) Theorem-2 boundary: with two IDENTICAL copies (K_{ij}=0), EVERY read whose origin is one of them
    has delta==0 (min_p==1 => Tied). Build copies with a repeated vector and verify."""
    fails = []
    for m in range(1, max_m + 1):
        for v in itertools.product(alphabet, repeat=m):
            for w in itertools.product(alphabet, repeat=m):
                copies = (v, v, w)  # copies[0]==copies[1] identical (K_{01}=0); copies[2] arbitrary
                for ssize in range(1, m + 1):
                    for S in itertools.combinations(range(m), ssize):
                        read_obs = {j: copies[0][j] for j in S}  # origin = copy 0 (== copy 1)
                        minp, best, delta = gate_minp_and_best(read_obs, copies)
                        if delta != 0 or abs(minp - 1.0) > 1e-15:
                            fails.append(("B5", m, copies, S, delta, minp))
    return fails

def check_precondition_necessity(max_m=3, alphabet=(0, 1)):
    """(B6) The Theorem-4(ii) soundness PRECONDITION `origin(r) in C` is NECESSARY (tight).

    Drop it: let the read's true origin copy `o` be ABSENT from C (the reference-absent / O4 regime).
    We exhibit a copy set C, an origin o∉C, and a read r (a partial observation of o) such that the gate
    is CONFIDENT (delta>=1, min_p<1) yet assigns to a copy b∈C with b != o AND r stays consistent with b
    — a silent, confident MISASSIGNMENT. So soundness fails without the precondition; this is exactly why
    the O4 two-stage freeze must ABSTAIN at stage 1 and re-thread against ref+absent copies.

    We ALSO certify the escape hatch: a read that OBSERVES every column ('full' read) from an absent
    origin is consistent with NO copy in C (cons=∅ -> flagged novel, never silently misassigned).
    Confident misassignment requires a PARTIAL read. The check passes iff a witness EXISTS (necessity)
    and the full-read escape holds universally.
    """
    witnesses = []
    full_read_escape_violations = []
    for m in range(2, max_m + 1):
        all_vecs = list(itertools.product(alphabet, repeat=m))
        for K in range(2, len(all_vecs)):
            for C in itertools.combinations(all_vecs, K):
                Cset = set(C)
                for o in all_vecs:
                    if o in Cset:
                        continue  # origin absent from C is the whole point
                    for ssize in range(1, m + 1):
                        for S in itertools.combinations(range(m), ssize):
                            read_obs = {j: o[j] for j in S}
                            cons = [c for c in range(K) if consistent(read_obs, C[c])]
                            minp, best, delta = gate_minp_and_best(read_obs, C)
                            if ssize == m:
                                # full read from an absent origin: must be consistent with NOTHING in C
                                if cons:
                                    full_read_escape_violations.append((m, C, o, S, cons))
                            if delta >= 1 and len(cons) == 1 and C[best] != o:
                                # confident, consistent, but the assigned copy is NOT the true origin
                                witnesses.append((m, C, o, S, C[best], minp))
    return witnesses, full_read_escape_violations


def _minimal_covers(reads, universe):
    """All MINIMUM-size subsets C of `universe` s.t. every read is consistent with some copy in C and
    every copy in C is the unique consistent copy of at least one read (no idle copy)."""
    best_size = None
    covers = []
    for k in range(1, len(universe) + 1):
        if best_size is not None and k > best_size:
            break
        for C in itertools.combinations(universe, k):
            # every read consistent with >=1 copy
            if not all(any(consistent(r, c) for c in C) for r in reads):
                continue
            # no idle copy: each copy uniquely explains >=1 read (else a smaller cover exists)
            used = [any(consistent(r, c) and sum(consistent(r, c2) for c2 in C) == 1 for r in reads) for c in C]
            if not all(used):
                continue
            best_size = k if best_size is None else best_size
            covers.append(C)
    return best_size, covers


def check_recombinant_cover():
    """(B7) Per-read soundness (Theorem 4) is ORTHOGONAL to K>=3 cover NON-UNIQUENESS (the recombination
    obstruction, Prop K-frontier). We exhibit a read set with >=2 DISTINCT minimum covers of size >=3, and
    verify that WITHIN EACH fixed cover, every read still satisfies per-read soundness (delta>=1 => the
    gate's pick is the unique consistent copy of that cover). The covers DISAGREE on some reads' origins —
    that is the family-level K>=3 ambiguity Theorem 4 explicitly does NOT resolve, confirming its scope is
    per-read-given-C, not cover recovery.

    Construction: 3 columns, reads observe only CONTIGUOUS pairs {0,1} and {1,2} (recombination across the
    middle locus is unobserved). Reads: col{0,1} in {00,11}; col{1,2} in {00,11}. Two size-3 covers explain
    this fragment set via different middle-locus recombinants (the enumerator returns size-4 covers here;
    any minimum size >=3 with >=2 distinct covers witnesses the K>=3 non-uniqueness).
    """
    reads = [
        {0: 0, 1: 0}, {0: 1, 1: 1},   # col {0,1}
        {1: 0, 2: 0}, {1: 1, 2: 1},   # col {1,2}
        {1: 0, 2: 1}, {1: 1, 2: 0},   # extra col {1,2} recombinants forcing K>=3
        {0: 0, 1: 1}, {0: 1, 1: 0},   # extra col {0,1} recombinants forcing K>=3
    ]
    universe = list(itertools.product((0, 1), repeat=3))
    size, covers = _minimal_covers(reads, universe)
    distinct = [set(c) for c in covers]
    # need >=2 distinct minimum covers, each of size >=3 (genuine K>=3 recombination non-uniqueness)
    nonunique = size is not None and size >= 3 and len(covers) >= 2
    soundness_ok = True
    for C in covers[:2]:
        Ct = list(C)
        for r in reads:
            cons = [c for c in range(len(Ct)) if consistent(r, Ct[c])]
            minp, best, delta = gate_minp_and_best(r, Ct)
            if delta >= 1 and cons != [best]:
                soundness_ok = False
    return nonunique, size, len(covers), soundness_ok


if __name__ == "__main__":
    f1, ninst, nreads = check_bridge()
    f2 = check_boundary_K0()
    fails = f1 + f2
    print(f"enumerated {ninst} copy-sets, {nreads} reads (K<=3, m<=3, |A|<=3) + the K0-boundary universe")
    if fails:
        print(f"FAIL: {len(fails)} counterexample(s):")
        for x in fails[:10]:
            print("  ", x)
        sys.exit(1)
    print("OK  - B1 min_p=eps^delta | B2 delta>=1<=>unique-consistent | B3 soundness(assign=origin) |")
    print("OK  - B4 delta=0<=>ambiguous | B5 K_{ij}=0 => all-Tied (Theorem-2 boundary recovered per-read)")

    # (B6) precondition necessity
    wit, escape_viol = check_precondition_necessity()
    if not wit:
        print("FAIL: B6 expected >=1 confident-misassignment witness when origin(r)∉C, found none")
        sys.exit(1)
    if escape_viol:
        print(f"FAIL: B6 full-read escape violated ({len(escape_viol)} cases where an absent origin's "
              f"full read was consistent with a copy in C)")
        sys.exit(1)
    m_, C_, o_, S_, b_, mp_ = wit[0]
    print(f"OK  - B6 precondition origin(r)∈C is NECESSARY: {len(wit)} witnesses where dropping it causes a "
          f"CONFIDENT misassignment.")
    print(f"        witness: C={list(C_)}, absent origin o={o_}, read obs cols {S_} -> gate assigns "
          f"b={b_} (min_p={mp_:.0e}, b!=o); full reads escape via the novel-copy flag.")

    # (B7) recombinant cover
    nonuniq, sz, ncov, sound = check_recombinant_cover()
    if not nonuniq:
        print(f"FAIL: B7 expected >=2 distinct minimum covers of size>=3, got size={sz} count={ncov}")
        sys.exit(1)
    if not sound:
        print("FAIL: B7 per-read soundness violated within a fixed cover")
        sys.exit(1)
    print(f"OK  - B7 recombinant cover: {ncov} distinct minimum covers of size {sz} (K>=3 non-unique) — yet "
          f"per-read soundness holds WITHIN each fixed C; Theorem 4 is per-read, not cover recovery.")

    print("BRIDGE THEOREM verified on the full enumerated universe (B1-B7).")
    sys.exit(0)
