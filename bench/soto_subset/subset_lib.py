#!/usr/bin/env python3
"""Soto's families as our families cut by copy number: the statistics, as pure functions (no I/O).

Spec: docs/superpowers/specs/2026-10-14-soto-subset-cn-design.md; prereg: docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md.

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
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score

FAM90A = "ID_356"  # register 1170: ARI swings +-0.035 on whether this family stays whole
CN_GAP = 2.0       # Soto's per-pair rule, MAD < 1; for two values the unscaled MAD is |a - b| / 2


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


# ----------------------------------------------------------------------------------------------------------- T2

def cn_components(block_of, cn):
    """The intersection of our blocks with copy number: inside each block, its members that have a value, sorted by
    value, split wherever two neighbours differ by >= CN_GAP. That is single linkage under Soto's pair rule
    (|a - b| < CN_GAP joins), whose components on a line are exactly these runs. -> {member: 'block#run'}."""
    by = collections.defaultdict(list)
    for m in sorted(cn):
        by[block_of[m]].append((cn[m], m))
    comp = {}
    for b in sorted(by):
        run, prev = 0, None
        for v, m in sorted(by[b]):
            if prev is not None and v - prev >= CN_GAP:
                run += 1
            comp[m] = f"{b}#{run}"
            prev = v
    return comp


def cn_shuffle(block_of, cn, rng):
    """Copy-number values permuted among the members of each block; members without a value stay without."""
    by = collections.defaultdict(list)
    for m in sorted(cn):
        by[block_of[m]].append(m)
    out = {}
    for b in sorted(by):
        vals = [cn[m] for m in by[b]]
        rng.shuffle(vals)
        out.update(zip(by[b], vals))
    return out


def _codes(labels):
    index = {}
    return [index.setdefault(x, len(index)) for x in labels]


def ari(truth_of, pred_of, members):
    """Adjusted Rand index over `members`; NaN for fewer than 2 members."""
    ms = sorted(members)
    if len(ms) < 2:
        return float("nan")
    return float(adjusted_rand_score(_codes(truth_of[m] for m in ms), _codes(pred_of[m] for m in ms)))


def _components(comp_of, keep):
    comps = collections.defaultdict(set)
    for m in keep:
        comps[comp_of[m]].add(m)
    return comps


def exact_count(soto_sets, comp_of, members):
    """Distinct Soto member sets, restricted to `members` and with >= 2 of them, equal to a component restricted to
    `members`."""
    keep = set(members)
    comp_sets = {frozenset(s) for s in _components(comp_of, keep).values()}
    targets = {frozenset(s & keep) for s in soto_sets.values()}
    return sum(1 for s in targets if len(s) >= 2 and s in comp_sets)


def hungarian_jaccard(soto_sets, comp_of, members):
    """Mean one-to-one best-match Jaccard (Hungarian) between Soto's families (>= 2 of `members`) and the components
    on `members`; an unmatched family scores 0. NaN when no family qualifies."""
    keep = set(members)
    rows = [s & keep for _, s in sorted(soto_sets.items())]
    rows = [s for s in rows if len(s) >= 2]
    if not rows:
        return float("nan")
    comps = _components(comp_of, keep)
    cols = [comps[k] for k in sorted(comps)]
    col_of = {m: j for j, s in enumerate(cols) for m in s}
    jac = np.zeros((len(rows), len(cols)))
    for i, s in enumerate(rows):
        for j in {col_of[m] for m in s}:
            inter = len(s & cols[j])
            jac[i, j] = inter / (len(s) + len(cols[j]) - inter)
    r, c = linear_sum_assignment(jac, maximize=True)
    return float(jac[r, c].sum() / len(rows))


def meet(soto_of, block_of):
    """The coarsest common refinement of Soto's single-label partition and ours: {member: 'label|block'}."""
    return {m: f"{soto_of[m]}|{block_of[m]}" for m in soto_of}


def scores(soto_sets, soto_of, comp_of, members):
    """One partition of `members` (comp_of) against Soto: exact and ARI over the single-label members (ARI also
    without FAM90A), and the cover-aware Hungarian Jaccard over all of `members`."""
    single = [m for m in members if m in soto_of]
    no90 = [m for m in single if soto_of[m] != FAM90A]
    return dict(members=len(members), single=len(single),
                exact=exact_count(soto_sets, comp_of, single),
                ari=ari(soto_of, comp_of, single),
                ari_no_fam90a=ari(soto_of, comp_of, no90),
                hungarian_jaccard=hungarian_jaccard(soto_sets, comp_of, members))


def t2(soto_sets, soto_of, block_of, cn, n_null, seed):
    """T2 for one copy-number table: the intersection scored against Soto, and the within-block shuffle null for the
    exact count and ARI. Members are the keys of `cn`."""
    members = sorted(cn)
    single = [m for m in members if m in soto_of]
    obs = scores(soto_sets, soto_of, cn_components(block_of, cn), members)
    rng = random.Random(seed)
    null_exact, null_ari = [], []
    for _ in range(n_null):
        comp = cn_components(block_of, cn_shuffle(block_of, cn, rng))
        null_exact.append(exact_count(soto_sets, comp, single))
        null_ari.append(ari(soto_of, comp, single))
    p_exact, p_ari = p_upper(null_exact, obs["exact"]), p_upper(null_ari, obs["ari"])
    return dict(observed=obs,
                null=dict(n=n_null, seed=seed,
                          exact_mean=float(np.mean(null_exact)) if null_exact else float("nan"),
                          exact_max=max(null_exact) if null_exact else None,
                          ari_mean=float(np.mean(null_ari)) if null_ari else float("nan"),
                          ari_max=float(max(null_ari)) if null_ari else None),
                p_exact=p_exact, p_ari=p_ari, passes=bool(p_exact < 0.001 and p_ari < 0.001))
