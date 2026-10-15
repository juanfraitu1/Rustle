#!/usr/bin/env python3
"""Soto's families as our families cut by copy number: the statistics, as pure functions (no I/O).

Spec: docs/superpowers/specs/2026-10-14-soto-subset-cn-design.md; prereg: docs/PREREG_soto_subset_cn_2026-10-14.md.

Conventions. `block_of` {member: block id} covers every member of a universe; each of our singletons has a block of
its own. `soto_sets` {label: frozenset(member)} are Soto's families on that universe; a member can be in several
(Soto's table is a cover). `soto_of` {member: label} is the single-label view: the members with exactly one Soto
label, an `Unassigned_*` gene being its own label. `multi` = the members with two or more Soto labels. `cn`
{member: float} holds only the members that have a copy-number value. Ids are strings, and every loop that feeds a
random generator runs in sorted order, so results do not depend on dict order.
"""
import collections
import math
import random

import numpy as np


# ----------------------------------------------------------------------------------------------------------- T1

def contained(members, block_of):
    """True iff all `members` lie in one block."""
    return len({block_of[m] for m in members}) == 1


def containment_count(soto_sets, block_of):
    """Soto families with >= 2 members whose members all lie in one block of ours."""
    return sum(1 for s in soto_sets.values() if len(s) >= 2 and contained(s, block_of))


def size_null(soto_sets, block_of, n, seed):
    """Containment counts of n random partitions with our exact block-size multiset, singletons included: each draw
    is a uniform permutation of the members' block labels."""
    members = sorted(block_of)
    labels = [block_of[m] for m in members]
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        rng.shuffle(labels)
        out.append(containment_count(soto_sets, dict(zip(members, labels))))
    return out


def p_upper(null, observed):
    """One-sided permutation p-value, (1 + draws >= observed) / (1 + draws); NaN when `observed` is NaN."""
    if isinstance(observed, float) and math.isnan(observed):
        return float("nan")
    return (1 + sum(1 for v in null if v >= observed)) / (1 + len(null))


def cond_entropy(x_of, y_of, members):
    """H(X | Y) in bits over `members`."""
    joint = collections.Counter((x_of[m], y_of[m]) for m in members)
    marg = collections.Counter(y_of[m] for m in members)
    n = sum(joint.values())
    return sum(c / n * math.log2(marg[y] / c) for (_, y), c in joint.items()) if n else 0.0


def direction(soto_of, block_of):
    """D = H(Soto | ours) - H(ours | Soto) over the single-label members; D > 0 when ours is the coarser."""
    members = sorted(soto_of)
    return cond_entropy(soto_of, block_of, members) - cond_entropy(block_of, soto_of, members)


def units(soto_of, block_of):
    """Connected components of the single-label members linked by a shared Soto label or a shared block of ours,
    each a sorted member list. Every Soto block and every block of ours lies inside one unit, so D decomposes over
    units: they are the bootstrap's resampling unit."""
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for m in sorted(soto_of):
        a, b = ("s", soto_of[m]), ("b", block_of[m])
        parent.setdefault(a, a)
        parent.setdefault(b, b)
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    groups = collections.defaultdict(list)
    for m in sorted(soto_of):
        groups[find(("s", soto_of[m]))].append(m)
    return [groups[k] for k in sorted(groups)]


def direction_bootstrap(soto_of, block_of, n, seed):
    """95% percentile interval of D over n resamples of `units` with replacement. A unit drawn twice enters as two
    disjoint copies (its Soto labels and our blocks both relabelled by draw), so a repeat never fakes a split.
    -> (lo, hi); (NaN, NaN) without members."""
    us = units(soto_of, block_of)
    if not us:
        return float("nan"), float("nan")
    rng = random.Random(seed)
    ds = []
    for _ in range(n):
        s_of, b_of = {}, {}
        for i in range(len(us)):
            for m in us[rng.randrange(len(us))]:
                s_of[(i, m)] = (i, soto_of[m])
                b_of[(i, m)] = (i, block_of[m])
        ds.append(direction(s_of, b_of))
    lo, hi = np.percentile(ds, [2.5, 97.5])
    return float(lo), float(hi)


def classify_violation(members, block_of, block_size, multi):
    """For a Soto family that is not contained, the first rule that makes it contained: 'coverage' (our singletons
    removed), then 'cover' (Soto's multi-family members also removed), else 'conflict'. A remainder of <= 1 member
    counts as contained."""
    rest = [m for m in members if block_size[block_of[m]] > 1]
    if len({block_of[m] for m in rest}) <= 1:
        return "coverage"
    rest = [m for m in rest if m not in multi]
    if len({block_of[m] for m in rest}) <= 1:
        return "cover"
    return "conflict"


def restrict(soto_sets, keep):
    """Soto's families restricted to the members in `keep`."""
    keep = set(keep)
    return {lab: frozenset(s & keep) for lab, s in soto_sets.items()}


def t1(soto_sets, soto_of, block_of, multi, clustered, n_null, n_boot, seed):
    """T1 on one universe. `clustered` = the members on loci of our families (plan 1's universe). -> dict."""
    block_size = collections.Counter(block_of.values())
    fams = {lab: s for lab, s in soto_sets.items() if len(s) >= 2}
    obs = containment_count(soto_sets, block_of)
    null = size_null(soto_sets, block_of, n_null, seed)
    p = p_upper(null, obs)
    lo, hi = direction_bootstrap(soto_of, block_of, n_boot, seed)
    by_size = {}
    for name, pick in (("le3", lambda s: len(s) <= 3), ("gt3", lambda s: len(s) > 3)):
        sel = [s for s in fams.values() if pick(s)]
        by_size[name] = dict(families=len(sel), contained=sum(1 for s in sel if contained(s, block_of)))
    violations = [dict(label=lab, size=len(fams[lab]), blocks=len({block_of[m] for m in fams[lab]}),
                       cls=classify_violation(sorted(fams[lab]), block_of, block_size, multi))
                  for lab in sorted(fams) if not contained(fams[lab], block_of)]
    classes = collections.Counter(v["cls"] for v in violations)
    p1 = restrict(soto_sets, clustered)
    members = sorted(soto_of)
    return dict(
        members=len(block_of), blocks=len(block_size), soto_labels=len(soto_sets), families_ge2=len(fams),
        contained=obs, rate_all=obs / len(soto_sets) if soto_sets else float("nan"),
        rate_ge2=obs / len(fams) if fams else float("nan"),
        null=dict(n=n_null, seed=seed, mean=float(np.mean(null)) if null else float("nan"),
                  max=max(null) if null else None),
        p=p, t1a_pass=p < 0.001,
        H_soto_given_ours=cond_entropy(soto_of, block_of, members),
        H_ours_given_soto=cond_entropy(block_of, soto_of, members),
        D=direction(soto_of, block_of), D_ci95=[lo, hi], t1b_pass=lo > 0,
        by_size=by_size,
        plan1=dict(evaluable=sum(1 for s in p1.values() if len(s) >= 2), contained=containment_count(p1, block_of)),
        violation_classes={k: classes.get(k, 0) for k in ("coverage", "cover", "conflict")},
        violations=violations,
    )
