# Soto ⊂ ours, cut by copy number — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A pre-registered test of whether Soto's families are a refinement of our O1 families (T1) and whether
cutting ours by copy number with Soto's pair rule rebuilds them (T2). It runs on the guided universe, with the two
shipped de novo catalogs as the replication.

**Architecture:** A small Python bench module, `bench/soto_subset/`, in four files:
- `subset_lib.py`: the statistics, as pure functions;
- `cn_tables.py`: copy number per gene, by CAT v4 gene id;
- `universes.py`: the member universes, built with plan 1's mappers;
- `run_subset.py`: a CLI that writes JSON per universe.

The pre-registration is committed after the code and before any statistic is computed. Three light runs follow. The
results, register rows and the thesis line come last.

**Tech Stack:** Python 3, numpy 2.4.2, scipy 1.17.1 (`linear_sum_assignment`), scikit-learn 1.8.0
(`adjusted_rand_score`) and pytest 9.0.3, all installed. It reuses `bench/o1_hier/hier_match.py` and
`bench/soto/soto_replication.py`.

**Spec:** `docs/superpowers/specs/2026-10-14-soto-subset-cn-design.md` (commit bdbc4dfb, approved by the user
2026-10-14).

## Global Constraints

- **Repo and branch.** Repo `/mnt/linuxdisk/home/juanfraitu/rustle_m2_soto`. Work on branch `soto-subset`, created
  from local `main`. Never push; merging is the user's decision.
- **Dependencies.** Python only: stdlib plus numpy, scipy, scikit-learn and pytest as installed. No new dependency.
  No Rust change, and no change to O1's family definition or defaults.
- **Tests.** `python3 -m pytest bench/soto_subset -q -p no:cacheprovider`. `bench/o1_hier` and `bench/soto` must stay
  green: `python3 -m pytest bench/o1_hier bench/soto -q -p no:cacheprovider`.
- **Real data.** Every command that reads real data runs through `bash tools/rlock.sh light ...` (each is < 2 GB and
  < 3 min; a full guided run takes about 30 s). Run in the foreground, and never use `pkill -f`.
- **Commits.**
  - Path-limited only: `git add <paths>` and `git commit -- <paths>`. Never `git add -A` or `git add .`.
  - Never stage the other session's untracked `bench/hierarchy/*` or `bench/ideal_expression/*`.
  - Every message ends with two lines: `Co-Authored-By: <the committing agent's model name> <noreply@anthropic.com>`
    and `Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z`.
- **Joins.** Copy number joins by CAT v4 gene id, never by gene name.
- **Fixed constants:**
  - seed 7; 1,000 null draws; 1,000 bootstrap resamples;
  - p = (1 + draws ≥ observed) / (1 + draws); each bar is p < 0.001;
  - pair rule: |CN_a − CN_b| < 2 joins (Soto's per-pair MAD < 1).
- **FAM90A.** FAM90A is `ID_356`. ARI is also reported without it (register 1170).
- **Scope.** Human only; never pool with gorilla.
- **Outputs.** `/mnt/linuxdisk/tmp/o1_soto_subset/<universe>/` with `<universe>` = `guided`, `human_A119b`,
  `human_testis`.
- **Pre-registration.** `docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md` is committed before any T1 or T2 statistic is
  computed. After that only its `## Results` section grows.

## Rulings on the spec (made while planning; Task 6 copies them into the prereg)

- **R1 — bootstrap unit.** The T1b bootstrap resamples *units*, not Soto's blocks. A unit is a connected component of
  the single-label members, linked by a shared Soto label or a shared block of ours.
  - Resampling Soto's blocks alone puts two copies of a drawn block inside one block of ours. That fakes a split,
    inflates H(Soto | ours), and so biases the test toward its bar.
  - D decomposes exactly over units.
- **R2 — a fourth file.** `universes.py` holds the member universes, beside the three files the spec names. This keeps
  `subset_lib.py` pure and `run_subset.py` thin.
- **R3 — exact over single-label members.** The exact count, like ARI, is computed over the single-label members, with
  Soto's families restricted to them. Soto's own scorer does the same: register 1162's 479/491 excludes the 149
  multi-family genes. The cover-aware Hungarian Jaccard is computed over all members.
- **R4 — recipe rows on guided only.** Soto's recipe is gene-level, and de novo members are loci.
- **R5 — ceiling on single-label members.** The ceiling (Soto ∩ ours) is scored on single-label members only,
  Hungarian Jaccard included: the meet of a cover with a partition is not a partition.
- **R6 — two rates.** T1 reports two rates:
  - over all 491 families, which includes families with < 2 members in the universe;
  - plan 1's: families with ≥ 2 members on clustered members, contained within that restriction.

## Review Focus

- **A member missing from a map fails loudly, naming the id.** A gene without a gene line in the guided GFF raises an
  error. Copy-number values of genes outside the universe are ignored, never counted as members. (Task 4 tests both.)
- **A copy-number difference of exactly 2 cuts.** Assembly values are integers or halves, and Soto's rule needs
  MAD = |Δ| / 2 < 1. (Task 2 test.)
- **Two Soto families with one member set (after restriction) count once in `exact`.** The count is over distinct
  sets. (Task 2 test.)
- **Degenerate universes fail the bar without crashing.** A universe with no family of ≥ 2 members, or a single
  valued member, returns 0 or NaN. (Task 1 and Task 2 tests.)
- **Dict order and file row order never change a result.** The same seed gives byte-identical `t1.json`, `t2.json`
  and `violations.tsv`. (Task 1, Task 2 and Task 5 tests.)

## File structure

| File | Responsibility |
|---|---|
| `bench/soto_subset/subset_lib.py` | T1 and T2 statistics: containment, size null, entropies, unit bootstrap, violation classes, the copy-number intersection and its shuffle null, ARI, exact, Hungarian Jaccard, meet |
| `bench/soto_subset/cn_tables.py` | the three copy-number tables by gene id, per-gene accounting, the id/coordinate check, the recipe's famCN table; a JSON report CLI |
| `bench/soto_subset/universes.py` | guided and de novo universes: members, blocks, Soto labels, clustered members, genes per member, coverage; copy number per member |
| `bench/soto_subset/run_subset.py` | the CLI: `inputs.json`, `t1.json`, `violations.tsv`, `t2.json`, `recipe.json` |
| `bench/soto_subset/test_*.py` | one test file per module, with real-data wiring tests that are skipped when the inputs are absent |
| `docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md` | the pre-registration, then the results |

---

### Task 1: T1 statistics (`subset_lib.py`, part 1)

**Files:**
- Create: `bench/soto_subset/subset_lib.py`
- Test: `bench/soto_subset/test_subset_lib.py`

**Interfaces:**
- Consumes: nothing.
- Produces, used by Task 5:
  - `t1(soto_sets, soto_of, block_of, multi, clustered, n_null, n_boot, seed) -> dict` with the keys `members`,
    `blocks`, `soto_labels`, `families_ge2`, `contained`, `rate_all`, `rate_ge2`, `null` {n, seed, mean, max}, `p`,
    `t1a_pass`, `H_soto_given_ours`, `H_ours_given_soto`, `D`, `D_ci95` [lo, hi], `t1b_pass`, `by_size` {le3, gt3:
    {families, contained}}, `plan1` {evaluable, contained}, `violation_classes` {coverage, cover, conflict}, and
    `violations` [{label, size, blocks, cls}].
  - Also `p_upper(null, observed) -> float`, used by Task 2.
  - The types: `soto_sets` {str: frozenset(str)}, `soto_of` {str: str}, `block_of` {str: str},
    `multi` and `clustered` frozenset(str).

- [ ] **Step 1: Write the failing tests**

Create `bench/soto_subset/test_subset_lib.py`:

```python
import collections
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import subset_lib as L  # noqa: E402

# Our blocks: B1 = {a, b, c, d, i}, B2 = {e, f, j}, singletons g and h.
BLOCK_OF = {"a": "B1", "b": "B1", "c": "B1", "d": "B1", "i": "B1", "e": "B2", "f": "B2", "j": "B2",
            "g": "Lg", "h": "Lh"}
# Soto: S1 and S2 contained; S3 needs our singleton g (coverage); S4 shares d with S2 (cover); S5 is a conflict.
SOTO_SETS = {"S1": frozenset("ab"), "S2": frozenset("cd"), "S3": frozenset("eg"), "S4": frozenset("df"),
             "S5": frozenset("ij")}
MULTI = frozenset("d")
SOTO_OF = {"a": "S1", "b": "S1", "c": "S2", "e": "S3", "g": "S3", "f": "S4", "i": "S5", "j": "S5", "h": "U_h"}


# ------------------------------------------------------------------------------------------------------------ T1

def test_containment_count():
    assert L.containment_count(SOTO_SETS, BLOCK_OF) == 2


def test_classify_violation_each_class():
    size = collections.Counter(BLOCK_OF.values())
    assert L.classify_violation(sorted(SOTO_SETS["S3"]), BLOCK_OF, size, MULTI) == "coverage"
    assert L.classify_violation(sorted(SOTO_SETS["S4"]), BLOCK_OF, size, MULTI) == "cover"
    assert L.classify_violation(sorted(SOTO_SETS["S5"]), BLOCK_OF, size, MULTI) == "conflict"


def test_conditional_entropies_by_hand():
    members = sorted(SOTO_OF)
    # H(Soto | ours): B1 holds S1, S1, S2, S5 (1.5 bits, weight 4/9); B2 holds S3, S4, S5 (log2 3, weight 3/9).
    assert abs(L.cond_entropy(SOTO_OF, BLOCK_OF, members) - (2 / 3 + math.log2(3) / 3)) < 1e-12
    # H(ours | Soto): S3 and S5 each split over two blocks (1 bit, weight 2/9 each).
    assert abs(L.cond_entropy(BLOCK_OF, SOTO_OF, members) - 4 / 9) < 1e-12
    assert abs(L.direction(SOTO_OF, BLOCK_OF) - (2 / 3 + math.log2(3) / 3 - 4 / 9)) < 1e-12


def test_p_upper():
    assert L.p_upper([0, 1, 2, 3], 3) == 2 / 5
    assert L.p_upper([], 5) == 1.0
    assert math.isnan(L.p_upper([0.1, 0.2], float("nan")))


def test_size_null_is_seeded_and_order_free():
    a = L.size_null(SOTO_SETS, BLOCK_OF, 50, 7)
    b = L.size_null(SOTO_SETS, dict(reversed(list(BLOCK_OF.items()))), 50, 7)
    assert a == b and len(a) == 50


def test_size_null_probability():
    # Blocks of sizes 2 and 1 over x, y, z: the family {x, y} is contained with probability 1/3.
    null = L.size_null({"S": frozenset("xy")}, {"x": "B", "y": "B", "z": "C"}, 3000, 7)
    assert abs(sum(null) / 3000 - 1 / 3) < 0.035


def test_one_block_never_beats_its_null():
    one = {m: "B" for m in BLOCK_OF}
    obs = L.containment_count(SOTO_SETS, one)
    assert obs == 5 and L.p_upper(L.size_null(SOTO_SETS, one, 20, 7), obs) == 1.0


def test_units():
    assert L.units(SOTO_OF, BLOCK_OF) == [list("abcefgij"), ["h"]]


def test_bootstrap_identical_partitions_is_zero():
    same = {m: "B" + SOTO_OF[m] for m in SOTO_OF}
    assert L.direction_bootstrap(SOTO_OF, same, 50, 7) == (0.0, 0.0)


def test_bootstrap_one_block_three_soto_pairs():
    soto = {"a": "S1", "b": "S1", "c": "S2", "d": "S2", "e": "S3", "f": "S3"}
    lo, hi = L.direction_bootstrap(soto, {m: "B" for m in soto}, 50, 7)
    assert abs(lo - math.log2(3)) < 1e-12 and abs(hi - math.log2(3)) < 1e-12


def test_bootstrap_is_seeded_and_order_free():
    a = L.direction_bootstrap(SOTO_OF, BLOCK_OF, 200, 7)
    b = L.direction_bootstrap(dict(reversed(list(SOTO_OF.items()))), BLOCK_OF, 200, 7)
    assert a == b and a[0] <= L.direction(SOTO_OF, BLOCK_OF) <= a[1]


def test_t1_on_the_toy():
    clustered = frozenset(m for m, b in BLOCK_OF.items() if b in ("B1", "B2"))
    r = L.t1(SOTO_SETS, SOTO_OF, BLOCK_OF, MULTI, clustered, 20, 20, 7)
    assert r["contained"] == 2 and r["families_ge2"] == 5 and r["rate_all"] == 2 / 5
    assert r["violation_classes"] == {"coverage": 1, "cover": 1, "conflict": 1}
    assert r["by_size"] == {"le3": {"families": 5, "contained": 2}, "gt3": {"families": 0, "contained": 0}}
    assert r["plan1"] == {"evaluable": 4, "contained": 2}
    assert [v["label"] for v in r["violations"]] == ["S3", "S4", "S5"]
    assert r["t1a_pass"] == (r["p"] < 0.001) and r["t1b_pass"] == (r["D_ci95"][0] > 0)


def test_t1_without_any_family():
    r = L.t1({"S": frozenset("a")}, {"a": "S"}, {"a": "B", "b": "C"}, frozenset(), frozenset("ab"), 10, 10, 7)
    assert r["contained"] == 0 and r["families_ge2"] == 0 and r["p"] == 1.0
    assert not r["t1a_pass"] and not r["t1b_pass"]
```

- [ ] **Step 2: Run them; they fail**

Run: `cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto && python3 -m pytest bench/soto_subset/test_subset_lib.py -q -p no:cacheprovider`
Expected: collection error, `ModuleNotFoundError: No module named 'subset_lib'`.

- [ ] **Step 3: Write the T1 part of `subset_lib.py`**

Create `bench/soto_subset/subset_lib.py`:

```python
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
```

- [ ] **Step 4: Run the tests; they pass**

Run: `python3 -m pytest bench/soto_subset/test_subset_lib.py -q -p no:cacheprovider`
Expected: `13 passed`.

- [ ] **Step 5: Commit**

```bash
git add bench/soto_subset/subset_lib.py bench/soto_subset/test_subset_lib.py
git commit -m "bench: soto_subset T1 statistics (containment, size null, entropy direction, unit bootstrap, violation classes)

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z" -- bench/soto_subset/subset_lib.py bench/soto_subset/test_subset_lib.py
```

---

### Task 2: T2 statistics (`subset_lib.py`, part 2)

**Files:**
- Modify: `bench/soto_subset/subset_lib.py` (add imports and constants at the top, and append the T2 section)
- Test: `bench/soto_subset/test_subset_lib.py` (append)

**Interfaces:**
- Consumes: `p_upper(null, observed)` from Task 1.
- Produces, used by Task 5:
  - `FAM90A = "ID_356"` and `CN_GAP = 2.0`;
  - `cn_components(block_of, cn) -> {member: "block#run"}`;
  - `scores(soto_sets, soto_of, comp_of, members) -> {members, single, exact, ari, ari_no_fam90a, hungarian_jaccard}`;
  - `meet(soto_of, block_of) -> {member: "label|block"}`;
  - `t2(soto_sets, soto_of, block_of, cn, n_null, seed) -> {observed: scores, null: {n, seed, exact_mean, exact_max,
    ari_mean, ari_max}, p_exact, p_ari, passes}`.

- [ ] **Step 1: Append the failing tests to `bench/soto_subset/test_subset_lib.py`**

```python
# ------------------------------------------------------------------------------------------------------------ T2

def test_cn_components_runs():
    block_of = {"a": "B1", "b": "B1", "c": "B1", "d": "B1", "e": "B2", "f": "B2", "g": "B2"}
    cn = {"a": 2.0, "b": 3.5, "c": 5.0, "d": 9.0, "e": 2.0, "f": 3.0}
    comp = L.cn_components(block_of, cn)
    assert comp["a"] == comp["b"] == comp["c"] != comp["d"]
    assert comp["e"] == comp["f"] and comp["e"] != comp["a"] and "g" not in comp


def test_cn_gap_of_exactly_two_cuts():
    # Soto's rule keeps a pair when MAD = |a - b| / 2 < 1: a difference of exactly 2 cuts, 1.999 joins.
    comp = L.cn_components({"x": "B", "y": "B", "z": "B"}, {"x": 2.0, "y": 4.0, "z": 5.999})
    assert comp["x"] != comp["y"] and comp["y"] == comp["z"]


def test_cn_shuffle_keeps_each_block_multiset():
    block_of = {"a": "B1", "b": "B1", "c": "B1", "d": "B2", "e": "B2"}
    cn = {"a": 1.0, "b": 2.0, "c": 3.0, "d": 7.0, "e": 9.0}
    s1 = L.cn_shuffle(block_of, cn, random.Random(7))
    s2 = L.cn_shuffle(block_of, dict(reversed(list(cn.items()))), random.Random(7))
    assert s1 == s2 and set(s1) == set(cn)
    assert sorted(s1[m] for m in "abc") == [1.0, 2.0, 3.0] and sorted(s1[m] for m in "de") == [7.0, 9.0]


def test_ari_by_hand():
    t = {"a": "S1", "b": "S1", "c": "S2", "d": "S2"}
    assert L.ari(t, t, "abcd") == 1.0
    assert L.ari(t, {m: "one" for m in "abcd"}, "abcd") == 0.0
    assert math.isnan(L.ari(t, t, "a"))


def test_exact_counts_distinct_sets():
    soto = {"S1": frozenset("abc"), "S2": frozenset("abc"), "S3": frozenset("dx")}
    comp = {"a": "c1", "b": "c1", "c": "c1", "d": "c2"}
    assert L.exact_count(soto, comp, "abcd") == 1  # {a, b, c} counted once; S3 keeps one member inside


def test_hungarian_jaccard_by_hand():
    soto = {"S1": frozenset("abc"), "S2": frozenset("de")}
    comp = {"a": "C1", "b": "C1", "c": "C2", "d": "C2", "e": "C2"}
    assert abs(L.hungarian_jaccard(soto, comp, "abcde") - 2 / 3) < 1e-12
    one = {m: "C" for m in "abcd"}
    assert abs(L.hungarian_jaccard({"S1": frozenset("ab"), "S2": frozenset("cd")}, one, "abcd") - 0.25) < 1e-12
    assert math.isnan(L.hungarian_jaccard({"S1": frozenset("a")}, {"a": "C"}, "a"))


def test_meet():
    assert L.meet({"a": "S1", "b": "S1", "c": "S2"}, {"a": "B1", "b": "B2", "c": "B1"}) == {
        "a": "S1|B1", "b": "S1|B2", "c": "S2|B1"}


def test_scores_drops_fam90a_only_for_its_ari():
    soto_of = {"a": "S1", "b": "S1", "c": L.FAM90A, "d": L.FAM90A}
    soto = {"S1": frozenset("ab"), L.FAM90A: frozenset("cd")}
    comp = {"a": "x", "b": "x", "c": "y", "d": "z"}
    s = L.scores(soto, soto_of, comp, "abcd")
    assert s["exact"] == 1 and s["single"] == 4 and s["ari_no_fam90a"] == 1.0 and s["ari"] < 1.0


def test_t2_reproduces_soto_and_null_rate():
    soto = {"S1": frozenset("ab"), "S2": frozenset("cd")}
    soto_of = {"a": "S1", "b": "S1", "c": "S2", "d": "S2"}
    cn = {"a": 2.0, "b": 2.0, "c": 10.0, "d": 10.0}
    r = L.t2(soto, soto_of, {m: "B" for m in "abcd"}, cn, 3000, 7)
    assert r["observed"]["exact"] == 2 and r["observed"]["ari"] == 1.0
    # 2 of the 6 distinct arrangements of {2, 2, 10, 10} rebuild Soto, so p is near 1/3 and the bar fails.
    assert abs(r["p_exact"] - 1 / 3) < 0.035 and not r["passes"]


def test_t2_strong_signal_has_small_p():
    labels = "abcdefghijkl"
    soto_of = {m: f"S{i // 2}" for i, m in enumerate(labels)}
    soto = {f"S{k}": frozenset(labels[2 * k:2 * k + 2]) for k in range(6)}
    cn = {m: 10.0 * (i // 2) for i, m in enumerate(labels)}
    r = L.t2(soto, soto_of, {m: "B" for m in labels}, cn, 50, 7)
    # 720 of 12! / 2^6 = 7,484,400 arrangements rebuild Soto (about 1e-4 per draw).
    assert r["observed"]["exact"] == 6 and r["p_exact"] < 0.05 and r["p_ari"] < 0.05
    assert r["passes"] == (r["p_exact"] < 0.001 and r["p_ari"] < 0.001)


def test_t2_with_one_member_is_not_a_pass():
    r = L.t2({"S1": frozenset("a")}, {"a": "S1"}, {"a": "B"}, {"a": 2.0}, 10, 7)
    assert math.isnan(r["p_ari"]) and not r["passes"]
```

- [ ] **Step 2: Run them; the new ones fail**

Run: `python3 -m pytest bench/soto_subset/test_subset_lib.py -q -p no:cacheprovider`
Expected: 11 failures with `AttributeError: module 'subset_lib' has no attribute ...` (`cn_components`, `cn_shuffle`,
`ari`, `exact_count`, `hungarian_jaccard`, `meet`, `FAM90A`, `t2`); 13 passed.

- [ ] **Step 3: Add the T2 part**

In `bench/soto_subset/subset_lib.py`, replace the import block

```python
import numpy as np
```

with

```python
import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score

FAM90A = "ID_356"  # register 1170: ARI swings +-0.035 on whether this family stays whole
CN_GAP = 2.0       # Soto's per-pair rule, MAD < 1; for two values the unscaled MAD is |a - b| / 2
```

and append at the end of the file:

```python


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
```

- [ ] **Step 4: Run the tests; they pass**

Run: `python3 -m pytest bench/soto_subset/test_subset_lib.py -q -p no:cacheprovider`
Expected: `24 passed`.

- [ ] **Step 5: Commit**

```bash
git add bench/soto_subset/subset_lib.py bench/soto_subset/test_subset_lib.py
git commit -m "bench: soto_subset T2 statistics (copy-number intersection, within-block shuffle null, ARI, exact, Hungarian Jaccard, meet)

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z" -- bench/soto_subset/subset_lib.py bench/soto_subset/test_subset_lib.py
```

---

### Task 3: copy-number tables (`cn_tables.py`)

**Files:**
- Create: `bench/soto_subset/cn_tables.py`
- Test: `bench/soto_subset/test_cn_tables.py`

**Interfaces:**
- Consumes: `bench/soto/soto_replication.py`'s `load_famcn(path, column)` (test only).
- Produces, used by Task 5:
  - the constants `REPO`, `S1E`, `READ_DEPTH`, `S1C` and `TRUTH_V2` (default input paths);
  - `assembly_cn(path) -> (dict gene->float, set unresolved)`;
  - `read_depth_cn(path) -> dict`, `s1c_cn(path) -> dict`, `s1c_genes(path) -> sorted list`;
  - `account(genes, vals, unresolved) -> {valued, unresolved, absent}`;
  - `check_coordinates(s1e_path, truth_v2_path) -> (int ok, list disagreements, int gene_absent)`;
  - `write_table(path, vals)`.

- [ ] **Step 1: Write the failing tests**

Create `bench/soto_subset/test_cn_tables.py`:

```python
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cn_tables as CN  # noqa: E402

S1E_HDR = ["rid", "kind", "s1e_row", "gene", "gene_id", "chrom", "start", "end", "fam_HG002"]
REAL = all(os.path.exists(p) for p in (CN.S1E, CN.READ_DEPTH, CN.S1C, CN.TRUTH_V2))


def write(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)


def test_num():
    assert CN.num("3.5") == 3.5 and CN.num("0") == 0.0
    assert CN.num("") is None and CN.num("nan") is None and CN.num("N/A") is None
    assert CN.num("abc") is None and CN.num(None) is None


def test_assembly_median_of_resolved_rows(tmp_path):
    p = tmp_path / "s1e.tsv"
    write(p, S1E_HDR, [["S0", "S1E", 0, "A", "g1", "chr1", 100, 200, "18.0"],
                       ["S1", "S1E", 1, "A", "g1", "chr1", 300, 400, "20.0"],
                       ["S2", "S1E", 2, "A", "g1", "chr1", 500, 600, ""],
                       ["S3", "S1E", 3, "B", "g2", "chr2", 100, 200, ""],
                       ["C0", "CTRL", -1, "C", "", "chr3", 1, 20001, "2.0"]])
    vals, unresolved = CN.assembly_cn(str(p))
    assert vals == {"g1": 19.0} and unresolved == {"g2"}


def test_s1c_repeats_and_conflicts(tmp_path):
    p = tmp_path / "s1c.tsv"
    hdr = ["Gene ID", "Family ID", "Median famCN"]
    write(p, hdr, [["g1", "ID_1", "4.5"], ["g1", "ID_2", "4.5"], ["g2", "ID_1", "N/A"]])
    assert CN.s1c_cn(str(p)) == {"g1": 4.5} and CN.s1c_genes(str(p)) == ["g1", "g2"]
    write(p, hdr, [["g1", "ID_1", "4.5"], ["g1", "ID_2", "5.0"]])
    with pytest.raises(ValueError, match="g1"):
        CN.s1c_cn(str(p))


def test_read_depth_column(tmp_path):
    p = tmp_path / "rd.tsv"
    write(p, ["gene_id", "famCN", "famCN_sotoiv"], [["g1", "3.0", "2.5"], ["g2", "4.0", "nan"]])
    assert CN.read_depth_cn(str(p)) == {"g1": 2.5}


def test_account():
    assert CN.account(["g1", "g2", "g3"], {"g1": 2.0}, {"g2"}) == {"valued": 1, "unresolved": 1, "absent": 1}


def test_check_coordinates(tmp_path):
    truth = tmp_path / "t.tsv"
    write(truth, ["label", "gene_id", "chrom", "start", "end", "strand", "exons"],
          [["ID_1", "g1", "chr1", 101, 200, "+", "101-200"], ["ID_1", "g2", "chr1", 301, 400, "+", "301-400"]])
    s1e = tmp_path / "s1e.tsv"
    write(s1e, S1E_HDR, [["S0", "S1E", 0, "A", "g1", "chr1", 90, 101, "2"],    # bases 91-101: shares 101
                         ["S1", "S1E", 1, "B", "g2", "chr1", 400, 410, "2"],   # bases 401-410: past the body (400)
                         ["S2", "S1E", 2, "C", "g9", "chr1", 1, 50, "2"]])     # gene not in the lifted table
    ok, bad, absent = CN.check_coordinates(str(s1e), str(truth))
    assert ok == 1 and absent == 1 and [b[0] for b in bad] == ["S1"]


def test_write_table_is_read_by_the_recipe(tmp_path):
    sys.path.insert(0, os.path.join(HERE, "..", "soto"))
    import soto_replication as R
    p = tmp_path / "asm.tsv"
    CN.write_table(str(p), {"g1": 18.0, "g2": 2.5})
    assert R.load_famcn(str(p), "famCN") == {"g1": 18.0, "g2": 2.5}


@pytest.mark.skipif(not REAL, reason="the frozen inputs are not on this machine")
def test_real_tables(capsys):
    assert CN.main([]) == 0
    genes = CN.s1c_genes(CN.S1C)
    asm, unresolved = CN.assembly_cn(CN.S1E)
    assert len(genes) == 2334
    assert CN.account(genes, asm, unresolved) == {"valued": 1534, "unresolved": 257, "absent": 543}
    assert CN.account(genes, CN.read_depth_cn(CN.READ_DEPTH)) == {"valued": 2334, "unresolved": 0, "absent": 0}
    assert CN.account(genes, CN.s1c_cn(CN.S1C)) == {"valued": 1793, "unresolved": 0, "absent": 541}
    ok, bad, absent = CN.check_coordinates(CN.S1E, CN.TRUTH_V2)
    assert (ok, len(bad), absent) == (1710, 0, 121)
```

- [ ] **Step 2: Run them; they fail**

Run: `python3 -m pytest bench/soto_subset/test_cn_tables.py -q -p no:cacheprovider`
Expected: collection error, `ModuleNotFoundError: No module named 'cn_tables'`.

- [ ] **Step 3: Write `cn_tables.py`**

Create `bench/soto_subset/cn_tables.py`:

```python
#!/usr/bin/env python3
"""Copy number per Soto gene from three sources, keyed by CAT v4 gene id (an identifier, never a gene name).

    assembly    results_s1e.tsv (bench/soto/parcn_assembly.py, frozen 2026-09-29): one row per S1E region (gene ∩ SD98,
                CHM13 v2.0, BED half-open); `fam_HG002` = the HG002 diploid famCN analogue, empty when the region is
                unresolved (< 100 usable 30-mers). A gene's value = the median over its resolved rows.
    read depth  famcn_ours_allwssd.tsv: `gene_id`, `famCN_sotoiv` (268 SGDP genomes, Soto's interval).
    S1C         soto_famCN_S1C.tsv: `Gene ID`, `Median famCN`; a gene in several families repeats its one value.

check_coordinates() tests the id join against coordinates: every S1E region must overlap its gene's lifted body in
soto_truth.v2.tsv (GFF 1-based closed, CHM13 v2.0).

    cn_tables.py [--s1e P] [--read-depth P] [--s1c P] [--truth-v2 P]     -> JSON report; exit 1 on a disagreement
"""
import argparse
import collections
import csv
import json
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
S1E = "/mnt/linuxdisk/tmp/rustle_figures_dev/soto_parcn_asm/results_s1e.tsv"
READ_DEPTH = "/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/famcn_ours_allwssd.tsv"
S1C = os.path.join(REPO, "bench", "soto", "soto_famCN_S1C.tsv")
TRUTH_V2 = "/mnt/linuxdisk/tmp/o1_close/lift/soto_truth.v2.tsv"
MISSING = {"", "nan", "NaN", "NA", "N/A", "None"}


def num(v):
    """A table cell -> float, or None for a missing or non-numeric cell."""
    if v is None or v.strip() in MISSING:
        return None
    try:
        x = float(v)
    except ValueError:
        return None
    return None if math.isnan(x) else x


def assembly_cn(path, column="fam_HG002"):
    """-> (values {gene_id: median over resolved S1E rows}, unresolved {gene_id with S1E rows, none resolved})."""
    vals, seen = collections.defaultdict(list), set()
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["kind"] != "S1E":
                continue
            seen.add(r["gene_id"])
            v = num(r[column])
            if v is not None:
                vals[r["gene_id"]].append(v)
    out = {g: float(statistics.median(v)) for g, v in vals.items()}
    return out, seen - set(out)


def column_cn(path, id_col, value_col):
    """Any gene-keyed table -> {gene_id: value}. A gene repeated with a different value is an error."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            v = num(r[value_col])
            if v is None:
                continue
            g = r[id_col]
            if g in out and out[g] != v:
                raise ValueError(f"{path}: {g} has two values in {value_col}: {out[g]} and {v}")
            out[g] = v
    return out


def read_depth_cn(path, column="famCN_sotoiv"):
    return column_cn(path, "gene_id", column)


def s1c_cn(path):
    return column_cn(path, "Gene ID", "Median famCN")


def s1c_genes(path):
    """Every distinct `Gene ID` of S1C, sorted."""
    with open(path) as fh:
        return sorted({r["Gene ID"] for r in csv.DictReader(fh, delimiter="\t")})


def account(genes, vals, unresolved=frozenset()):
    """How each gene of `genes` stands in one table -> {valued, unresolved, absent}."""
    c = collections.Counter("valued" if g in vals else "unresolved" if g in unresolved else "absent" for g in genes)
    return {k: c.get(k, 0) for k in ("valued", "unresolved", "absent")}


def check_coordinates(s1e_path, truth_v2_path):
    """-> (n_ok, disagreements [(rid, gene_id, region, body)], n_gene_absent). A region (BED half-open) agrees when
    it shares >= 1 base with its gene's body (GFF closed) on the same chromosome."""
    body = {}
    with open(truth_v2_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            body[r["gene_id"]] = (r["chrom"], int(r["start"]), int(r["end"]))
    ok, bad, absent = 0, [], 0
    with open(s1e_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["kind"] != "S1E":
                continue
            b = body.get(r["gene_id"])
            if b is None:
                absent += 1
                continue
            s, e = int(r["start"]), int(r["end"])
            if r["chrom"] == b[0] and s < b[2] and e >= b[1]:
                ok += 1
            else:
                bad.append((r["rid"], r["gene_id"], f"{r['chrom']}:{s}-{e}", f"{b[0]}:{b[1]}-{b[2]}"))
    return ok, bad, absent


def write_table(path, vals, column="famCN"):
    """{gene_id: value} -> the 2-column TSV that soto_replication.py's load_famcn reads (gene_id, famCN)."""
    with open(path, "w") as fh:
        fh.write(f"gene_id\t{column}\n")
        fh.writelines(f"{g}\t{vals[g]!r}\n" for g in sorted(vals))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--s1e", default=S1E)
    ap.add_argument("--read-depth", default=READ_DEPTH)
    ap.add_argument("--s1c", default=S1C)
    ap.add_argument("--truth-v2", default=TRUTH_V2)
    a = ap.parse_args(argv)
    genes = s1c_genes(a.s1c)
    asm, unresolved = assembly_cn(a.s1e)
    ok, bad, absent = check_coordinates(a.s1e, a.truth_v2)
    rep = dict(genes=len(genes), assembly=account(genes, asm, unresolved),
               read_depth=account(genes, read_depth_cn(a.read_depth)), s1c=account(genes, s1c_cn(a.s1c)),
               coordinate_check=dict(ok=ok, disagree=len(bad), gene_absent=absent, first=bad[:5]))
    print(json.dumps(rep, indent=1))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run the tests (the real-data one reads the frozen inputs); they pass**

Run: `RLOCK_TIMEOUT=300 bash tools/rlock.sh light python3 -m pytest bench/soto_subset/test_cn_tables.py -q -p no:cacheprovider`
Expected: `8 passed`. The real-data test must PASS, not be skipped: the inputs exist on this machine.

- [ ] **Step 5: Commit**

```bash
git add bench/soto_subset/cn_tables.py bench/soto_subset/test_cn_tables.py
git commit -m "bench: soto_subset copy-number tables by CAT v4 gene id (assembly, read depth, S1C) and the id/coordinate check

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z" -- bench/soto_subset/cn_tables.py bench/soto_subset/test_cn_tables.py
```

---

### Task 4: universes (`universes.py`)

**Files:**
- Create: `bench/soto_subset/universes.py`
- Test: `bench/soto_subset/test_universes.py`

**Interfaces:**
- Consumes: from `bench/o1_hier/hier_match.py`:
  - `load_nodes(members_tsv) -> (nodes, {member: family})`;
  - `load_fold(path) -> {annotation key: representative key}`;
  - `load_loci_gff(path) -> {key: (chrom, strands, exons)}`;
  - `SpanIndex(items).query(chrom, s, e)` and `_strand_ok(a, b)`.
- Produces, used by Task 5:
  - `Universe` (namedtuple: `block_of`, `soto_sets`, `soto_of`, `multi`, `clustered`, `genes_of`
    {member: frozenset(gene)}, `coverage` dict);
  - `guided(s1c_path, loci_gff, members_tsv, fold_tsv) -> Universe`;
  - `denovo(truth_v2, loci_gff, members_tsv, fold_tsv) -> Universe`;
  - `member_cn(gene_cn, genes_of) -> {member: float}`.

- [ ] **Step 1: Write the failing tests**

Create `bench/soto_subset/test_universes.py`:

```python
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import universes as U  # noqa: E402

F_O = "/mnt/linuxdisk/tmp/o1_layers/v2/dupmasker"
F_O_GFF = "/mnt/linuxdisk/tmp/ourhom/loci.gff3"
SHIPPED = "/mnt/linuxdisk/tmp/o1_close/hier_shipped"
RUNS = "/mnt/linuxdisk/tmp/rustle_figures/runs"
TRUTH_V2 = "/mnt/linuxdisk/tmp/o1_close/lift/soto_truth.v2.tsv"
S1C = os.path.join(HERE, "..", "soto", "soto_famCN_S1C.tsv")
REAL = all(os.path.exists(p) for p in (F_O, F_O_GFF, SHIPPED, RUNS, TRUTH_V2, S1C))


def write(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)


def gff(path, genes):
    """genes: [(name, chrom, start, end, strand, [(s, e)] exons or None)], GFF 1-based closed."""
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for name, c, s, e, st, ex in genes:
            fh.write(f"{c}\t.\tgene\t{s}\t{e}\t.\t{st}\t.\tID=gene-{name};Name={name}\n")
            for i, (a, b) in enumerate(ex or []):
                fh.write(f"{c}\t.\texon\t{a}\t{b}\t.\t{st}\t.\tID=exon-{name}-{i};Parent=gene-{name}\n")


def members(path, fams):
    write(path, ["sd_complex", "superfamily", "family", "member", "expansions"],
          [["-", "-", fam, k, "-"] for fam, keys in fams.items() for k in keys])


def fold(path, pairs):
    write(path, ["annotation", "representative"], pairs)


def test_soto_views():
    labels = {"g1": {"ID_1"}, "g2": {"ID_1", "ID_2"}, "g3": {"Unassigned_9"}, "g4": {"ID_2"}}
    sets, soto_of, multi = U.soto_views(labels, ["ID_1", "ID_2", "ID_3", "Unassigned_9"])
    assert sets == {"ID_1": frozenset({"g1", "g2"}), "ID_2": frozenset({"g2", "g4"}), "ID_3": frozenset()}
    assert soto_of == {"g1": "ID_1", "g3": "Unassigned_9", "g4": "ID_2"} and multi == frozenset({"g2"})


def guided_fixture(tmp_path, extra_gene=False):
    s1c = tmp_path / "s1c.tsv"
    rows = [["g1", "ID_1"], ["g2", "ID_1"], ["g2", "ID_2"], ["g3", "ID_2"], ["g4", "Unassigned_7"]]
    write(s1c, ["Gene ID", "Family ID"], rows + ([["g9", "ID_3"]] if extra_gene else []))
    gff(tmp_path / "loci.gff3", [("g1", "chr1", 100, 200, "+", None), ("g2", "chr1", 300, 400, "+", None),
                                 ("g3", "chr1", 150, 180, "+", None), ("g4", "chr1", 900, 1000, "+", None)])
    members(tmp_path / "members.tsv", {"F1": ["chr1:100-200", "chr1:300-400"]})
    fold(tmp_path / "loci.tsv", [["chr1:150-180", "chr1:100-200"]])
    return [str(tmp_path / x) for x in ("s1c.tsv", "loci.gff3", "members.tsv", "loci.tsv")]


def test_guided_blocks(tmp_path):
    u = U.guided(*guided_fixture(tmp_path))
    assert u.block_of == {"g1": "F:F1", "g2": "F:F1", "g3": "F:F1", "g4": "L:chr1:900-1000"}
    assert u.clustered == frozenset({"g1", "g2", "g3"}) and u.multi == frozenset({"g2"})
    assert u.soto_sets == {"ID_1": frozenset({"g1", "g2"}), "ID_2": frozenset({"g2", "g3"})}
    assert u.soto_of == {"g1": "ID_1", "g3": "ID_2", "g4": "Unassigned_7"}
    assert u.coverage == dict(genes=4, genes_on_family_members=3, loci=3, loci_with_several_genes=1, families=1)
    assert u.genes_of["g3"] == frozenset({"g3"})


def test_guided_gene_without_gff_line_fails_loudly(tmp_path):
    with pytest.raises(ValueError, match="g9"):
        U.guided(*guided_fixture(tmp_path, extra_gene=True))


def test_denovo_reach(tmp_path):
    gff(tmp_path / "loci.gff3", [("L1", "chr1", 100, 500, "+", [(100, 150), (400, 500)]),
                                 ("L2", "chr1", 1000, 1200, "+", [(1000, 1200)]),
                                 ("L3", "chr1", 2000, 2100, "-", [(2000, 2100)]),
                                 ("L4", "chr1", 3000, 3100, "+", [(3000, 3100)])])
    members(tmp_path / "members.tsv", {"F1": ["chr1:100-500", "chr1:1000-1200"]})
    fold(tmp_path / "loci.tsv", [["chr1:3000-3100", "chr1:100-500"]])
    truth = tmp_path / "truth.tsv"
    write(truth, ["label", "gene_id", "chrom", "start", "end", "strand", "exons"],
          [["ID_1", "gA", "chr1", 120, 140, "+", "120-140"],      # exon inside L1's first exon
           ["ID_1", "gB", "chr1", 1100, 1150, "+", "1100-1150"],  # L2
           ["ID_2", "gC", "chr1", 2050, 2060, "+", "2050-2060"],  # L3 is on the other strand: no reach
           ["ID_2", "gD", "chr1", 3050, 3060, "+", "3050-3060"],  # L4, folded into L1
           ["ID_3", "gE", "chr1", 2010, 2020, "-", "2010-2020"],  # L3, a locus in no family
           ["ID_3", "gG", "chr1", 200, 300, "+", "200-300"],      # inside L1's intron: no reach
           ["ID_4", "gF", "chr1", 9000, 9100, "+", "9000-9100"]])  # nothing
    u = U.denovo(str(truth), str(tmp_path / "loci.gff3"), str(tmp_path / "members.tsv"), str(tmp_path / "loci.tsv"))
    assert u.block_of == {"chr1:100-500": "F:F1", "chr1:1000-1200": "F:F1", "chr1:2000-2100": "L:chr1:2000-2100"}
    assert u.soto_sets == {"ID_1": frozenset({"chr1:100-500", "chr1:1000-1200"}), "ID_2": frozenset({"chr1:100-500"}),
                           "ID_3": frozenset({"chr1:2000-2100"}), "ID_4": frozenset()}
    assert u.soto_of == {"chr1:1000-1200": "ID_1", "chr1:2000-2100": "ID_3"} and u.multi == frozenset({"chr1:100-500"})
    assert u.genes_of["chr1:100-500"] == frozenset({"gA", "gD"})
    assert u.coverage == dict(genes=7, records=7, records_on_members=3, records_on_unclustered_only=1,
                              records_without_locus=3, reached_loci=3, reached_clustered=2)
    cn = U.member_cn({"gA": 2.0, "gD": 4.0, "gB": 10.0, "gZ": 99.0}, u.genes_of)
    assert cn == {"chr1:100-500": 3.0, "chr1:1000-1200": 10.0}


@pytest.mark.skipif(not REAL, reason="the frozen inputs are not on this machine")
def test_real_universes_reproduce_plan1_mapping():
    g = U.guided(S1C, F_O_GFF, f"{F_O}/ours.hierarchy.members.tsv", f"{F_O}/ours.loci.tsv")
    assert g.coverage["genes"] == 2334 and g.coverage["genes_on_family_members"] == 1796 and g.coverage["loci"] == 2241
    expect = {"human_A119b": (1211, 155, 1092, 1396, 1115), "human_testis": (286, 100, 2072, 332, 238)}
    for s, (on, unc, none, reached, clustered) in expect.items():
        d = f"{SHIPPED}/{s}"
        u = U.denovo(TRUTH_V2, f"{RUNS}/{s}/{s}.fam.loci.gff3", f"{d}/ours.hierarchy.members.tsv", f"{d}/ours.loci.tsv")
        c = u.coverage
        assert (c["records_on_members"], c["records_on_unclustered_only"], c["records_without_locus"],
                c["reached_loci"], c["reached_clustered"]) == (on, unc, none, reached, clustered), s
        assert len(u.soto_sets) == 491, s
```

The real-data numbers are plan 1's coverage. Two sources:
- `/mnt/linuxdisk/tmp/o1_close/hier_shipped/<sample>/hier_match.json` → `coverage`: records_on_members 1,211 / 286,
  records_on_unclustered_locus 155 / 100, records_without_locus 1,092 / 2,072;
- the guided key mapping: 1,796 genes on F_O family members.

A mismatch means the mapping differs from plan 1's. Do not change the expected numbers; stop and report.

- [ ] **Step 2: Run them; they fail**

Run: `python3 -m pytest bench/soto_subset/test_universes.py -q -p no:cacheprovider`
Expected: collection error, `ModuleNotFoundError: No module named 'universes'`.

- [ ] **Step 3: Write `universes.py`**

Create `bench/soto_subset/universes.py`:

```python
#!/usr/bin/env python3
"""The two universes of the Soto-subset study: members, our blocks, Soto's labels, and per-member copy number.

guided  Members are Soto's 2,334 genes (S1C). A gene's locus is its gene line in the F_O input GFF
        (`/mnt/linuxdisk/tmp/ourhom/loci.gff3`, `Name=` = CAT v4 gene id), folded through the run's `ours.loci.tsv`.
        Genes on loci of one F_O family share the block `F:<family>`; genes on one unclustered locus share
        `L:<locus>`.
denovo  Members are the catalog loci that a gene of Soto's 491 families reaches by plan 1's exon + strand rule
        (soto_truth.v2.tsv, CHM13 v2.0; a locus folded into a representative counts as the representative). Blocks are
        `F:<family>` from the catalog's hierarchy members table, `L:<locus>` for a reached locus in no family.
A member's Soto labels are its genes' families; `Unassigned_*` genes are Soto singletons (guided only: they are not
in soto_truth.v2.tsv). Copy number per member = the median over its genes that have a value.
"""
import collections
import csv
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "o1_hier"))
from hier_match import SpanIndex, _strand_ok, load_fold, load_loci_gff, load_nodes  # noqa: E402

Universe = collections.namedtuple("Universe", "block_of soto_sets soto_of multi clustered genes_of coverage")


def s1c_labels(path):
    """S1C -> {gene_id: set(Family ID)}, `Unassigned_*` labels included."""
    out = collections.defaultdict(set)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Family ID"]:
                out[r["Gene ID"]].add(r["Family ID"])
    return dict(out)


def soto_views(labels_of, all_labels=()):
    """{member: set(label)} -> (soto_sets {label: frozenset(member)} without `Unassigned_*`, every label of
    `all_labels` present even when empty; soto_of {member: label} for members with exactly one label; multi = the
    members with >= 2 labels)."""
    sets = {lab: set() for lab in all_labels if not lab.startswith("Unassigned_")}
    for m, labs in labels_of.items():
        for lab in labs:
            if not lab.startswith("Unassigned_"):
                sets.setdefault(lab, set()).add(m)
    soto_of = {m: next(iter(labs)) for m, labs in labels_of.items() if len(labs) == 1}
    multi = frozenset(m for m, labs in labels_of.items() if len(labs) > 1)
    return {lab: frozenset(s) for lab, s in sets.items()}, soto_of, multi


def gff_gene_keys(path):
    """Gene lines of a GFF3 -> {Name: 'chrom:start-end'} (GFF 1-based, the hier_match member key)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            out[at["Name"]] = f"{f[0]}:{f[3]}-{f[4]}"
    return out


def _block(locus, family_of):
    return f"F:{family_of[locus]}" if locus in family_of else f"L:{locus}"


def guided(s1c_path, loci_gff, members_tsv, fold_tsv):
    labels = s1c_labels(s1c_path)
    key = gff_gene_keys(loci_gff)
    missing = sorted(set(labels) - set(key))
    if missing:
        raise ValueError(f"{len(missing)} S1C genes have no gene line in {loci_gff}: {missing[:5]}")
    _, family_of = load_nodes(members_tsv)
    fold = load_fold(fold_tsv)
    locus_of = {g: fold.get(key[g], key[g]) for g in labels}
    block_of = {g: _block(locus_of[g], family_of) for g in labels}
    soto_sets, soto_of, multi = soto_views(labels)
    clustered = frozenset(g for g, b in block_of.items() if b.startswith("F:"))
    per_locus = collections.Counter(locus_of.values())
    coverage = dict(genes=len(labels), genes_on_family_members=len(clustered), loci=len(per_locus),
                    loci_with_several_genes=sum(1 for c in per_locus.values() if c > 1),
                    families=len({b for b in block_of.values() if b.startswith("F:")}))
    return Universe(block_of, soto_sets, soto_of, multi, clustered, {g: frozenset([g]) for g in labels}, coverage)


def read_truth_v2(path):
    """soto_truth.v2.tsv -> ({gene_id: set(label)}, {gene_id: (chrom, strand, [(s, e)] GFF closed)}, n_records)."""
    labels, genes, n = collections.defaultdict(set), {}, 0
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n += 1
            labels[r["gene_id"]].add(r["label"])
            ex = [tuple(int(x) for x in b.split("-")) for b in r["exons"].split(",")]
            genes[r["gene_id"]] = (r["chrom"], r["strand"], ex)
    return dict(labels), genes, n


def denovo(truth_v2, loci_gff, members_tsv, fold_tsv):
    gene_labels, genes, n_records = read_truth_v2(truth_v2)
    loci = load_loci_gff(loci_gff)
    _, family_of = load_nodes(members_tsv)
    fold = load_fold(fold_tsv)
    idx = SpanIndex([(c, s, e, (st, k)) for k, (c, st, ex) in loci.items() for s, e in ex])
    reach = {g: frozenset(fold.get(k, k) for s, e in ex for st, k in idx.query(chrom, s, e)
                          if _strand_ok(st, {strand}))
             for g, (chrom, strand, ex) in sorted(genes.items())}
    labels_of, genes_of = collections.defaultdict(set), collections.defaultdict(set)
    for g, ms in reach.items():
        for m in ms:
            labels_of[m] |= gene_labels[g]
            genes_of[m].add(g)
    block_of = {m: _block(m, family_of) for m in labels_of}
    all_labels = sorted({lab for labs in gene_labels.values() for lab in labs})
    soto_sets, soto_of, multi = soto_views(labels_of, all_labels)
    clustered = frozenset(m for m in block_of if m in family_of)
    rows = lambda gs: sum(len(gene_labels[g]) for g in gs)  # noqa: E731  (one truth row per gene and label)
    on = [g for g, ms in reach.items() if any(m in family_of for m in ms)]
    unclustered = [g for g, ms in reach.items() if ms and not any(m in family_of for m in ms)]
    coverage = dict(genes=len(genes), records=n_records, records_on_members=rows(on),
                    records_on_unclustered_only=rows(unclustered),
                    records_without_locus=rows(g for g, ms in reach.items() if not ms),
                    reached_loci=len(block_of), reached_clustered=len(clustered))
    return Universe(block_of, soto_sets, soto_of, multi, clustered,
                    {m: frozenset(gs) for m, gs in genes_of.items()}, coverage)


def member_cn(gene_cn, genes_of):
    """{gene: value} -> {member: median over its genes that have a value}; members with none are left out, and
    values of genes outside the universe are never used."""
    out = {}
    for m, gs in genes_of.items():
        v = [gene_cn[g] for g in sorted(gs) if g in gene_cn]
        if v:
            out[m] = float(statistics.median(v))
    return out
```

- [ ] **Step 4: Run the tests (the real-data one reads the frozen inputs); they pass**

Run: `RLOCK_TIMEOUT=300 bash tools/rlock.sh light python3 -m pytest bench/soto_subset/test_universes.py -q -p no:cacheprovider`
Expected: `5 passed`. The real-data test must PASS, not be skipped.

- [ ] **Step 5: Commit**

```bash
git add bench/soto_subset/universes.py bench/soto_subset/test_universes.py
git commit -m "bench: soto_subset universes (guided genes on F_O loci, de novo loci by exon + strand; reproduces plan 1's coverage)

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z" -- bench/soto_subset/universes.py bench/soto_subset/test_universes.py
```

---

### Task 5: the runner (`run_subset.py`)

**Files:**
- Create: `bench/soto_subset/run_subset.py`
- Test: `bench/soto_subset/test_run_subset.py`

**Interfaces:**
- Consumes:
  - Task 1-2: `L.t1`, `L.t2`, `L.scores`, `L.cn_components`, `L.meet`, `L.FAM90A`;
  - Task 3: `CN.assembly_cn`, `CN.read_depth_cn`, `CN.s1c_cn`, `CN.account`, `CN.check_coordinates`,
    `CN.write_table`, and the path constants;
  - Task 4: `U.guided`, `U.denovo`, `U.member_cn`;
  - from `soto_replication.py`: `load_truth`, `restrict_universe`, `load_predicted`, `score`, `S1C`, and the CLI
    `genesets` / `cluster --pair-mad`.
- Produces: `main(argv) -> 0`, which writes `inputs.json`, `t1.json`, `violations.tsv`, `t2.json` and `recipe.json`
  (guided). Their keys are read by Tasks 6-9:
  - `t2.json` = {`common_members`, `common` {source: t2 dict}, `own` {source: scores}, `ours_alone`, `ceiling`,
    `t2_pass`, `read_depth_passes`};
  - `recipe.json` = {arm: {`full_gene_set` {ari, exact, truth_families}, `common_set` scores}} for the arms `none`,
    `s1c`, `read_depth`, `assembly`.

- [ ] **Step 1: Write the failing tests**

Create `bench/soto_subset/test_run_subset.py`:

```python
import json
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_subset as RS  # noqa: E402

REAL = all(os.path.exists(p) for p in (RS.F_O, RS.F_O_GFF, RS.SHIPPED, RS.RUNS, RS.CN.S1E, RS.CN.READ_DEPTH,
                                       RS.CN.TRUTH_V2))


def write(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)


def toy(tmp_path, s1e_chrom="chr1"):
    """Seven genes: ID_1 = {g1, g2} and ID_2 = {g3, g4} inside our family F1 at copy numbers 2 and 10; ID_3 =
    {g6, g7} with g7 on an unclustered locus; g5 Unassigned. -> the path overrides for main()."""
    loc = {"g1": ("chr1", 100, 200), "g2": ("chr1", 300, 400), "g3": ("chr1", 500, 600), "g4": ("chr1", 700, 800),
           "g5": ("chr1", 900, 1000), "g6": ("chr2", 100, 200), "g7": ("chr2", 300, 400)}
    lab = {"g1": "ID_1", "g2": "ID_1", "g3": "ID_2", "g4": "ID_2", "g5": "Unassigned_1", "g6": "ID_3", "g7": "ID_3"}
    s1c_cn = {"g1": "2.0", "g2": "2.1", "g3": "10.0", "g4": "10.2"}
    write(tmp_path / "s1c.tsv", ["Gene ID", "Family ID", "Median famCN"],
          [[g, lab[g], s1c_cn.get(g, "N/A")] for g in sorted(lab)])
    with open(tmp_path / "loci.gff3", "w") as fh:
        for g, (c, s, e) in sorted(loc.items()):
            fh.write(f"{c}\t.\tgene\t{s}\t{e}\t.\t+\t.\tID=gene-{g};Name={g}\n")
    fo = tmp_path / "fo"
    fo.mkdir()
    write(fo / "ours.hierarchy.members.tsv", ["sd_complex", "superfamily", "family", "member", "expansions"],
          [["-", "-", "F1", f"chr1:{s}-{e}", "-"] for g, (c, s, e) in sorted(loc.items()) if g in ("g1", "g2", "g3", "g4")]
          + [["-", "-", "F2", "chr2:100-200", "-"], ["-", "-", "F2", "chr2:500-600", "-"]])
    write(fo / "ours.loci.tsv", ["annotation", "representative"], [])
    asm = {"g1": 2, "g2": 2, "g3": 10, "g4": 10}
    write(tmp_path / "s1e.tsv", ["rid", "kind", "s1e_row", "gene", "gene_id", "chrom", "start", "end", "fam_HG002"],
          [[f"S{i}", "S1E", i, g, g, s1e_chrom, loc[g][1] - 1, loc[g][2], asm[g]] for i, g in enumerate(sorted(asm))])
    rd = {"g1": 2.1, "g2": 2.3, "g3": 9.8, "g4": 10.4, "g5": 3.0, "g6": 4.0, "g7": 8.0}
    write(tmp_path / "rd.tsv", ["gene_id", "famCN_sotoiv"], [[g, v] for g, v in sorted(rd.items())])
    write(tmp_path / "truth.tsv", ["label", "gene_id", "chrom", "start", "end", "strand", "exons"],
          [[lab[g], g, c, s, e, "+", f"{s}-{e}"] for g, (c, s, e) in sorted(loc.items()) if g != "g5"])
    return ["--s1c", str(tmp_path / "s1c.tsv"), "--s1e", str(tmp_path / "s1e.tsv"), "--read-depth",
            str(tmp_path / "rd.tsv"), "--truth-v2", str(tmp_path / "truth.tsv"), "--fo", str(fo), "--fo-gff",
            str(tmp_path / "loci.gff3")]


def load(path):
    with open(path) as fh:
        return json.load(fh)


def test_inputs_only(tmp_path):
    out = tmp_path / "out"
    assert RS.main(["guided", "--out", str(out), "--inputs-only", *toy(tmp_path)]) == 0
    r = load(out / "inputs.json")
    assert r["members"] == 7 and r["coverage"]["genes_on_family_members"] == 5
    assert r["cn_genes"]["assembly"] == {"valued": 4, "unresolved": 0, "absent": 3}
    assert r["cn_genes"]["read_depth"]["valued"] == 7 and r["cn_genes"]["s1c"]["valued"] == 4
    assert r["common_members"] == 4 and r["families_ge2_in_common"] == 2
    assert r["coordinate_check"] == {"ok": 4, "disagree": 0, "gene_absent": 0}
    assert sorted(os.listdir(out)) == ["inputs.json"]


def test_full_toy_run_is_deterministic(tmp_path):
    args = toy(tmp_path)
    for run in ("a", "b"):
        assert RS.main(["guided", "--out", str(tmp_path / run), "--null", "20", "--boot", "20", "--no-recipe",
                        *args]) == 0
    t1, t2 = load(tmp_path / "a" / "t1.json"), load(tmp_path / "a" / "t2.json")
    assert t1["contained"] == 2 and t1["families_ge2"] == 3
    assert t1["violation_classes"] == {"coverage": 1, "cover": 0, "conflict": 0}
    assert t2["common"]["assembly"]["observed"]["exact"] == 2 and t2["common"]["assembly"]["observed"]["ari"] == 1.0
    assert t2["ours_alone"]["exact"] == 0 and t2["ceiling"]["exact"] == 2
    with open(tmp_path / "a" / "violations.tsv") as fh:
        assert fh.read() == "label\tsize\tblocks\tclass\nID_3\t2\t2\tcoverage\n"
    for f in ("t1.json", "t2.json", "violations.tsv"):
        with open(tmp_path / "a" / f) as x, open(tmp_path / "b" / f) as y:
            assert x.read() == y.read(), f


def test_coordinate_disagreement_stops_the_run(tmp_path):
    with pytest.raises(SystemExit, match="disagrees"):
        RS.main(["guided", "--out", str(tmp_path / "out"), "--inputs-only", *toy(tmp_path, s1e_chrom="chr9")])


def test_denovo_needs_a_sample(tmp_path):
    with pytest.raises(SystemExit):
        RS.main(["denovo", "--out", str(tmp_path / "out")])


@pytest.mark.skipif(not REAL, reason="the frozen inputs are not on this machine")
def test_real_inputs_only(tmp_path):
    assert RS.main(["guided", "--out", str(tmp_path / "g"), "--inputs-only"]) == 0
    r = load(tmp_path / "g" / "inputs.json")
    assert r["members"] == 2334 and r["coverage"]["genes_on_family_members"] == 1796
    assert r["cn_genes"]["assembly"] == {"valued": 1534, "unresolved": 257, "absent": 543}
    assert r["cn_genes"]["s1c"] == {"valued": 1793, "unresolved": 0, "absent": 541}
    assert r["cn_genes"]["read_depth"]["valued"] == 2334
    assert r["common_members"] == 1534 and r["families_ge2_in_common"] == 357
    assert r["coordinate_check"] == {"ok": 1710, "disagree": 0, "gene_absent": 121}
    for s, on in (("human_A119b", 1211), ("human_testis", 286)):
        assert RS.main(["denovo", "--sample", s, "--out", str(tmp_path / s), "--inputs-only"]) == 0
        assert load(tmp_path / s / "inputs.json")["coverage"]["records_on_members"] == on
```

- [ ] **Step 2: Run them; they fail**

Run: `python3 -m pytest bench/soto_subset/test_run_subset.py -q -p no:cacheprovider`
Expected: collection error, `ModuleNotFoundError: No module named 'run_subset'`.

- [ ] **Step 3: Write `run_subset.py`**

Create `bench/soto_subset/run_subset.py`:

```python
#!/usr/bin/env python3
"""Soto's families as our families cut by copy number: T1 and T2 on one universe.

Spec docs/superpowers/specs/2026-10-14-soto-subset-cn-design.md; prereg docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md.

    run_subset.py guided --out DIR [--null 1000] [--boot 1000] [--seed 7] [--inputs-only] [--no-recipe]
    run_subset.py denovo --sample human_A119b|human_testis --out DIR [same options]

DIR receives inputs.json (universe and copy-number coverage; all that --inputs-only writes), t1.json, violations.tsv,
t2.json and, for guided, recipe.json: Soto's own recipe (`soto_replication.py cluster --pair-mad`) with each
copy-number source, scored on its full gene set (where the registered rungs must reproduce) and on T2's common set.
The run stops when the gene-id join disagrees with coordinates or a registered rung does not reproduce.
"""
import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cn_tables as CN  # noqa: E402
import subset_lib as L  # noqa: E402
import universes as U  # noqa: E402

SOTO = os.path.join(CN.REPO, "bench", "soto")
REPLICATION = os.path.join(SOTO, "soto_replication.py")
SHARED = os.path.join(SOTO, "shared_exons_5154_exon_mapback.tsv")
F_O = "/mnt/linuxdisk/tmp/o1_layers/v2/dupmasker"
F_O_GFF = "/mnt/linuxdisk/tmp/ourhom/loci.gff3"
SHIPPED = "/mnt/linuxdisk/tmp/o1_close/hier_shipped"
RUNS = "/mnt/linuxdisk/tmp/rustle_figures/runs"
SOURCES = ("assembly", "read_depth", "s1c")
# Soto's recipe with each source; the registered rungs (ARI, exact of 491) must reproduce (registers 1162, 1169).
RECIPE = (("none", ["--no-cn-gate"], (0.7307, 345)),
          ("s1c", ["--famcn", CN.S1C], (0.9698, 479)),
          ("read_depth", ["--famcn", CN.READ_DEPTH, "--famcn-col", "famCN_sotoiv"], (0.9277, 411)),
          ("assembly", ["--famcn", "{assembly}", "--famcn-col", "famCN"], None))


def load_tables(a):
    """The three copy-number tables by gene id; stops on any id/coordinate disagreement."""
    asm, unresolved = CN.assembly_cn(a.s1e)
    tables = dict(assembly=asm, read_depth=CN.read_depth_cn(a.read_depth), s1c=CN.s1c_cn(a.s1c))
    ok, bad, absent = CN.check_coordinates(a.s1e, a.truth_v2)
    if bad:
        raise SystemExit(f"the gene-id join disagrees with coordinates on {len(bad)} S1E rows, e.g. {bad[:3]}")
    return tables, unresolved, dict(ok=ok, disagree=0, gene_absent=absent)


def build_universe(a):
    if a.universe == "guided":
        return U.guided(a.s1c, a.fo_gff, os.path.join(a.fo, "ours.hierarchy.members.tsv"),
                        os.path.join(a.fo, "ours.loci.tsv"))
    d = os.path.join(a.shipped, a.sample)
    return U.denovo(a.truth_v2, os.path.join(a.runs, a.sample, f"{a.sample}.fam.loci.gff3"),
                    os.path.join(d, "ours.hierarchy.members.tsv"), os.path.join(d, "ours.loci.tsv"))


def common_members(member_cn):
    """The members with a value in every source, sorted."""
    return sorted(set.intersection(*(set(v) for v in member_cn.values())))


def inputs_report(uni, tables, unresolved, member_cn):
    genes = sorted(set().union(*uni.genes_of.values()))
    common = set(common_members(member_cn))
    return dict(coverage=uni.coverage, members=len(uni.block_of), soto_labels=len(uni.soto_sets),
                families_ge2=sum(1 for s in uni.soto_sets.values() if len(s) >= 2),
                single_label_members=len(uni.soto_of), multi_label_members=len(uni.multi),
                cn_genes={k: CN.account(genes, v, unresolved if k == "assembly" else frozenset())
                          for k, v in tables.items()},
                cn_members={k: len(v) for k, v in member_cn.items()},
                common_members=len(common),
                families_ge2_in_common=sum(1 for s in uni.soto_sets.values() if len(s & common) >= 2),
                fam90a_in_common=sum(1 for m in common if uni.soto_of.get(m) == L.FAM90A))


def run_t2(uni, member_cn, n_null, seed):
    common = common_members(member_cn)
    out = dict(common_members=len(common), common={}, own={})
    for name in SOURCES:
        cn = {m: member_cn[name][m] for m in common}
        out["common"][name] = L.t2(uni.soto_sets, uni.soto_of, uni.block_of, cn, n_null, seed)
        own = member_cn[name]
        out["own"][name] = L.scores(uni.soto_sets, uni.soto_of, L.cn_components(uni.block_of, own), sorted(own))
    out["ours_alone"] = L.scores(uni.soto_sets, uni.soto_of, uni.block_of, common)
    single = [m for m in common if m in uni.soto_of]
    out["ceiling"] = L.scores(uni.soto_sets, uni.soto_of, L.meet(uni.soto_of, uni.block_of), single)
    out["t2_pass"] = out["common"]["assembly"]["passes"]
    out["read_depth_passes"] = out["common"]["read_depth"]["passes"]
    return out


def run_recipe(out_dir, uni, assembly, common):
    """Soto's own recipe with each copy-number source (guided only: the recipe is gene-level)."""
    sys.path.insert(0, SOTO)
    import soto_replication as R
    g1793, g2334 = os.path.join(out_dir, "g1793.tsv"), os.path.join(out_dir, "g2334.tsv")
    subprocess.run([sys.executable, REPLICATION, "genesets", "--out-eligible", g1793, "--out-full", g2334],
                   check=True)
    asm_path = os.path.join(out_dir, "cn_assembly.tsv")
    CN.write_table(asm_path, assembly)
    clean, _ = R.load_truth(R.S1C)
    universe = R.restrict_universe(clean)
    rows = {}
    for arm, args, expect in RECIPE:
        pred_path = os.path.join(out_dir, f"recipe_{arm}.tsv")
        subprocess.run([sys.executable, REPLICATION, "cluster", "--pair-mad", "--shared", SHARED, "--geneset", g1793,
                        "--full-geneset", g2334, *[x.replace("{assembly}", asm_path) for x in args],
                        "--out", pred_path], check=True)
        pred = R.load_predicted(pred_path)
        full = R.score(clean, pred, universe)
        if expect and (abs(full["ari"] - expect[0]) > 5e-5 or full["n_exact"] != expect[1]):
            raise SystemExit(f"recipe {arm}: ARI {full['ari']:.4f}, exact {full['n_exact']} does not reproduce the "
                             f"registered {expect[0]} / {expect[1]}")
        comp_of = {g: pred.get(g) or f"U:{g}" for g in uni.block_of}
        rows[arm] = dict(full_gene_set=dict(ari=full["ari"], exact=full["n_exact"], truth_families=full["n_truth_fam"]),
                         common_set=L.scores(uni.soto_sets, uni.soto_of, comp_of, common))
    return rows


def write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True)
        fh.write("\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("universe", choices=["guided", "denovo"])
    ap.add_argument("--sample", choices=["human_A119b", "human_testis"], help="denovo: the shipped catalog")
    ap.add_argument("--out", required=True)
    ap.add_argument("--null", type=int, default=1000)
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--inputs-only", action="store_true", help="write inputs.json and stop (no statistic)")
    ap.add_argument("--no-recipe", action="store_true", help="guided: skip Soto's recipe rows")
    ap.add_argument("--s1c", default=CN.S1C)
    ap.add_argument("--s1e", default=CN.S1E)
    ap.add_argument("--read-depth", default=CN.READ_DEPTH)
    ap.add_argument("--truth-v2", default=CN.TRUTH_V2)
    ap.add_argument("--fo", default=F_O, help="guided: the F_O run directory (ours.hierarchy.members.tsv, ours.loci.tsv)")
    ap.add_argument("--fo-gff", default=F_O_GFF)
    ap.add_argument("--shipped", default=SHIPPED)
    ap.add_argument("--runs", default=RUNS)
    a = ap.parse_args(argv)
    if a.universe == "denovo" and not a.sample:
        ap.error("denovo needs --sample")
    os.makedirs(a.out, exist_ok=True)
    tables, unresolved, coord = load_tables(a)
    uni = build_universe(a)
    member_cn = {k: U.member_cn(v, uni.genes_of) for k, v in tables.items()}
    inputs = inputs_report(uni, tables, unresolved, member_cn)
    inputs["coordinate_check"] = coord
    write_json(os.path.join(a.out, "inputs.json"), dict(args=vars(a), **inputs))
    if a.inputs_only:
        print(json.dumps(inputs, indent=1, sort_keys=True))
        return 0
    t1 = L.t1(uni.soto_sets, uni.soto_of, uni.block_of, uni.multi, uni.clustered, a.null, a.boot, a.seed)
    write_json(os.path.join(a.out, "t1.json"), t1)
    with open(os.path.join(a.out, "violations.tsv"), "w") as fh:
        fh.write("label\tsize\tblocks\tclass\n")
        fh.writelines(f"{v['label']}\t{v['size']}\t{v['blocks']}\t{v['cls']}\n" for v in t1["violations"])
    t2 = run_t2(uni, member_cn, a.null, a.seed)
    write_json(os.path.join(a.out, "t2.json"), t2)
    if a.universe == "guided" and not a.no_recipe:
        write_json(os.path.join(a.out, "recipe.json"),
                   run_recipe(a.out, uni, tables["assembly"], common_members(member_cn)))
    print(json.dumps(dict(T1a=t1["t1a_pass"], T1b=t1["t1b_pass"], T2=t2["t2_pass"], contained=t1["contained"],
                          families_ge2=t1["families_ge2"], D=t1["D"], D_ci95=t1["D_ci95"],
                          t2_assembly=t2["common"]["assembly"]["observed"]), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run the module's tests and the neighbouring suites; they pass**

Run: `RLOCK_TIMEOUT=300 bash tools/rlock.sh light python3 -m pytest bench/soto_subset bench/o1_hier bench/soto -q -p no:cacheprovider`
Expected: every test passes, `test_real_inputs_only` included (not skipped). `bench/soto_subset` contributes
42 tests. Report the total line.

Do NOT run `run_subset.py` without `--inputs-only` on real data in this task: T1 and T2 statistics are computed only
after the prereg is committed (Task 6).

- [ ] **Step 5: Commit**

```bash
git add bench/soto_subset/run_subset.py bench/soto_subset/test_run_subset.py
git commit -m "bench: soto_subset runner (inputs-only report, T1/T2 JSON, Soto's recipe rows with the registered rungs as a wiring check)

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z" -- bench/soto_subset/run_subset.py bench/soto_subset/test_run_subset.py
```

---

### Task 6: the pre-registration

**Files:**
- Create: `docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md`

**Interfaces:**
- Consumes: `inputs.json` of the three universes (`--inputs-only`), and the code commit (HEAD after Task 5).
- Produces: the committed prereg that Tasks 7-9 add results to. It has the section headings `## Results` and, later,
  `### Guided`, `### De novo` and `### Verdict`.

- [ ] **Step 1: Write the input reports (no statistic is computed)**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
for u in guided human_A119b human_testis; do
  if [ "$u" = guided ]; then set -- guided; else set -- denovo --sample "$u"; fi
  RLOCK_TIMEOUT=300 bash tools/rlock.sh light python3 bench/soto_subset/run_subset.py "$@" \
      --out /mnt/linuxdisk/tmp/o1_soto_subset/$u --inputs-only > /dev/null
done
ls /mnt/linuxdisk/tmp/o1_soto_subset/*/
```

Expected: each directory holds only `inputs.json`.

Check that these keys equal the table in Step 2. Any difference means an input changed since planning: stop and
report it.

| Key | guided | human_A119b | human_testis |
|---|---|---|---|
| `members` | 2334 | 1396 | 332 |
| `families_ge2` | 491 | 293 | 105 |
| `single_label_members` / `multi_label_members` | 2185 / 149 | 1190 / 206 | 300 / 32 |
| `cn_members` assembly / read_depth / s1c | 1534 / 2334 / 1793 | 1175 / 1396 / 1196 | 288 / 332 / 292 |
| `common_members` | 1534 | 1175 | 288 |
| `families_ge2_in_common` | 357 | 257 | 85 |
| `fam90a_in_common` | 2 | 0 | 0 |
| `coordinate_check` ok / disagree / gene_absent | 1710 / 0 / 121 | same | same |

The plan-1 wiring keys are in `coverage`:
- guided: `genes_on_family_members` 1796, `loci` 2241;
- A119b: `records_on_members` 1211, `records_on_unclustered_only` 155, `records_without_locus` 1092;
- testis: `records_on_members` 286, `records_on_unclustered_only` 100, `records_without_locus` 2072.

- [ ] **Step 2: Write the prereg**

Create `docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md` with exactly this content. Replace `<CODE>` with the output of
`git rev-parse --short HEAD`, which is Task 5's last commit.

````markdown
# Soto's families as our families cut by copy number: T1 and T2 (2026-10-14)

**Pre-registered before any T1 or T2 statistic is computed.**
- Spec: `docs/superpowers/specs/2026-10-14-soto-subset-cn-design.md`.
- Plan: `docs/superpowers/plans/2026-10-14-soto-subset-cn.md`.
- Code: `bench/soto_subset/` at commit `<CODE>`. Every number below comes from it, run with the defaults (1,000 null
  draws, 1,000 bootstrap resamples, seed 7):

      bash tools/rlock.sh light python3 bench/soto_subset/run_subset.py guided --out /mnt/linuxdisk/tmp/o1_soto_subset/guided
      bash tools/rlock.sh light python3 bench/soto_subset/run_subset.py denovo --sample human_A119b --out /mnt/linuxdisk/tmp/o1_soto_subset/human_A119b
      bash tools/rlock.sh light python3 bench/soto_subset/run_subset.py denovo --sample human_testis --out /mnt/linuxdisk/tmp/o1_soto_subset/human_testis

## Question

- **T1, Soto ⊂ ours:** are Soto's families (Soto 2025, Table S1C) a refinement of our O1 multi-copy gene families?
- **T2, Soto ≈ ours ∩ CN:** does cutting ours by copy number, with Soto's own pair rule, rebuild them?

Human only. Nothing is tuned on Soto's table. The pair rule is Soto's (|Δ| < 2), and the nulls, draw counts and seed
are fixed above.

## Inputs (frozen)

| What | Path |
|---|---|
| Soto's families and the S1C copy number | `bench/soto/soto_famCN_S1C.tsv` |
| Assembly copy number (HG002 diploid, `fam_HG002`) | `/mnt/linuxdisk/tmp/rustle_figures_dev/soto_parcn_asm/results_s1e.tsv` |
| Read-depth copy number (`famCN_sotoiv`, 268 SGDP genomes) | `/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication/famcn_ours_allwssd.tsv` |
| Soto's genes on CHM13 v2.0 (de novo mapping, coordinate check) | `/mnt/linuxdisk/tmp/o1_close/lift/soto_truth.v2.tsv` |
| Guided families (F_O, the 10-09 run) | `/mnt/linuxdisk/tmp/o1_layers/v2/dupmasker/ours.hierarchy.members.tsv`, `ours.loci.tsv`; loci `/mnt/linuxdisk/tmp/ourhom/loci.gff3` |
| De novo catalogs (exact-mcl; the same partitions as the shipped `<sample>.fam.clusters.tsv`) | `/mnt/linuxdisk/tmp/o1_close/hier_shipped/<sample>/ours.hierarchy.members.tsv`, `ours.loci.tsv`; loci `/mnt/linuxdisk/tmp/rustle_figures/runs/<sample>/<sample>.fam.loci.gff3` |
| Soto's recipe | `bench/soto/soto_replication.py cluster --pair-mad` on `bench/soto/shared_exons_5154_exon_mapback.tsv` |

## Known before running (recorded now)

- **T1 on the guided universe is not blind.** On F_O (the 10-09 run, at member level, over plan 1's 400 evaluable
  families), about 78% of Soto's families lie inside one of our families. 215 equal one of our expansions and 346 lie
  in one superfamily (`docs/archive/2026-10/merged/PREREG_family_hierarchy_flag_2026-10-09.md`;
  `docs/superpowers/specs/2026-10-13-close-o1-design.md`). The de novo catalogs are the replication.
- **Older catalogs nest too.** On the older RNA lattice catalogs, 80-90% of Soto's families nest in ours at L1/L2.
  On held-out chr5/7/21 it was 61/76: 53/59 families of ≤ 3 genes and 8/17 larger ones
  (`bench/SOTO_AS_A_REFINEMENT.md`).
- **Copy number already raised agreement on our edges.** In an exploratory, unregistered run (2026-10-07), our F_O
  graph's edges gated by Soto's pair rule gave ARI 0.7142 / 0.7849 / 0.8263 (no copy number / read depth / S1C) at
  identity ≥ 0.70. This study cuts our families rather than our edges, and adds assembly copy number.
- **Soto's recipe rungs** (registers 1162 and 1169; reproduced 2026-10-14), over the 491 families: no copy number
  0.7307 / 345 exact, read depth 0.9277 / 411, S1C 0.9698 / 479.
- **Assembly vs S1C:** assembly copy number agrees with S1C at Spearman 0.617 (human; register 1172).
- **No held-out half is needed.** S1C's family-hash halves are both spent for copy-number split rules (2026-10-07).
  Nothing here selects a rule, so no half is needed; the de novo catalogs replicate.
- **FAM90A.** ARI over Soto's families swings ±0.035 with FAM90A (register 1170). Hence exact families are a
  co-requirement of T2's bar, and ARI is also reported without FAM90A.
- **Input counts** (`inputs.json`; inputs, not outcomes):

| | guided | A119b | testis |
|---|---|---|---|
| members | 2,334 genes | 1,396 loci | 332 loci |
| Soto families with ≥ 2 members | 491 | 293 | 105 |
| single-label / multi-label members | 2,185 / 149 | 1,190 / 206 | 300 / 32 |
| members on our families | 1,796 | 1,115 | 238 |
| members with assembly / read-depth / S1C copy number | 1,534 / 2,334 / 1,793 | 1,175 / 1,396 / 1,196 | 288 / 332 / 292 |
| common set (a value in all three) | 1,534 | 1,175 | 288 |
| Soto families with ≥ 2 members in the common set | 357 | 257 | 85 |
| FAM90A members in the common set | 2 | 0 | 0 |

## Definitions

- **Guided universe.** Members are Soto's 2,334 genes. A gene's locus is its gene line in the F_O input GFF (by CAT
  v4 gene id), folded through `ours.loci.tsv`.
  - Genes on loci of one F_O family share a block.
  - Genes on one unclustered locus share a block.
  - Every other gene is our singleton.
- **De novo universe.** Members are the catalog loci reached by a gene of Soto's 491 families: an exon shares ≥ 1 bp,
  on the same strand, with a locus exon, and a folded locus counts as its representative.
  - Blocks are the catalog's families. A reached locus in no family is our singleton.
  - A locus carries the families of every gene that reaches it.
- **Soto's labels.** `Unassigned_*` genes are Soto singletons. A member with two or more labels is *multi-label*.
- **Copy number of a member** = the median over its genes that have a value.
- **Common set** = the members with a value in all three sources.

## T1 — Soto ⊂ ours

- **Containment.** A Soto family with ≥ 2 members is contained when all its members lie in one block of ours. The
  statistic is the number contained.
- **Size null.** Random partitions of the members with our exact block-size multiset, singletons included: a uniform
  permutation of the block labels, 1,000 draws, seed 7. p = (1 + draws ≥ observed) / 1,001.
- **T1a bar:** p < 0.001.
- **Direction.** Over the single-label members, D = H(Soto | ours) − H(ours | Soto), in bits.
  - The 95% percentile interval comes from 1,000 bootstrap resamples of *units* (ruling R1), seed 7.
  - A unit is a connected component of the members linked by a shared Soto label or a shared block of ours.
  - A unit drawn twice enters as two disjoint copies.
- **T1b bar:** the interval lies above 0.
- **Reported, no bar:**
  - the rate over all 491 families, and over the families with ≥ 2 members (R6);
  - plan 1's rate: families with ≥ 2 members on clustered members, contained within that restriction;
  - containment for families of ≤ 3 and > 3 members;
  - every violation, classified by the first rule that makes it contained: **coverage** (our singletons removed),
    then **cover** (multi-label members also removed), else **conflict**.

## T2 — Soto ≈ ours ∩ CN

- **The intersection.** Inside each block of ours, members are joined when |CN_a − CN_b| < 2 (Soto's per-pair
  MAD < 1). The components are the intersection: runs of the block's sorted values, split at gaps ≥ 2.
- **Scores,** on the same members:
  - **exact:** distinct Soto member sets, restricted to the single-label members with ≥ 2 of them, that equal a
    component (R3);
  - **ARI** over the single-label members, also reported without FAM90A (`ID_356`);
  - **Hungarian Jaccard:** the mean one-to-one best-match Jaccard over Soto's families with ≥ 2 members, all members
    included.
- **Null.** Copy-number values are permuted among the members of each block of ours, 1,000 draws, seed 7, with the
  same p formula. It is computed for each source.
- **T2 bar:** with assembly copy number on the common set, the exact count and ARI both beat their null at p < 0.001.
- **Reported, no bar:**
  - read depth and S1C on the common set, with their nulls;
  - each source on its own coverage;
  - our blocks alone;
  - the ceiling: Soto's partition intersected with ours, on single-label members (R5);
  - Soto's recipe with each source, guided only (R4), on its full gene set and on the common set.

## Rulings made while planning (before any statistic)

- **R1 — bootstrap unit.** T1b resamples units, not Soto's blocks. Two copies of a Soto block inside one block of
  ours would fake a split and inflate H(Soto | ours), biasing toward the bar. D decomposes exactly over units.
- **R2 — a fourth file.** `universes.py` holds the universes beside the spec's three files.
- **R3 — exact over single-label members.** Exact is counted over single-label members, as Soto's own scorer does:
  479/491 excludes the 149 multi-family genes.
- **R4 — recipe rows on guided only.** Soto's recipe is gene-level.
- **R5 — the ceiling on single-label members.** The meet of a cover is not a partition.
- **R6 — two rates.** The rate over all 491 families (families with < 2 members count as not contained) and plan 1's
  rate are both reported.

## Readings (fixed now)

"T1 holds" means T1a and T1b both pass. A de novo "fail" means a bar fails on either catalog.

| Guided | De novo | Reading |
|---|---|---|
| T1 and T2 hold | hold on both catalogs | O1's statement: "Soto's families are a refinement of our multi-copy gene families; cutting ours by copy number rebuilds them (exact x, ARI y; ceiling z)". Replaces the failed hierarchy claim (register 1277) as the Soto sentence. |
| T1 and T2 hold | fail | The same statement, labelled as shown on the guided universe only. |
| T1 holds, T2 fails | — | Soto ⊂ ours is claimed. The copy-number part is claimed only if read depth passes T2's bar on the same set, labelled human-only and population-based. |
| T1 fails | — | No subset claim; T2 is reported, not claimed. The size and violation tables say where Soto and we disagree. |

## Wiring checks (the run stops, and nothing is read, if one fails)

- The gene-id join agrees with coordinates: 0 disagreements.
- Soto's recipe rungs reproduce: 0.7307 / 345, 0.9698 / 479 and 0.9277 / 411.
- Plan 1's coverage reproduces: guided, 1,796 genes on our families; A119b, 1,211 / 155 / 1,092 Soto records on
  clustered loci / only unclustered / none; testis, 286 / 100 / 2,072.

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.
````

- [ ] **Step 3: Commit the prereg (before any statistic)**

```bash
git add docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md
git commit -m "prereg: Soto's families as our families cut by copy number (T1 containment + direction, T2 intersection vs within-block shuffle; guided, de novo replication)

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z" -- docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md
git log --oneline -1 -- docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md
```

---

### Task 7: the guided run and its results

**Files:**
- Modify: `docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md` (append under `## Results`)
- Outputs (not committed): `/mnt/linuxdisk/tmp/o1_soto_subset/guided/`

**Interfaces:**
- Consumes: `t1.json`, `t2.json`, `recipe.json` and `violations.tsv` (keys as in Task 5).
- Produces: the `### Guided (2026-10-14)` results section, which Task 9 reads.

- [ ] **Step 1: Confirm the prereg is committed, then run**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
git log --oneline -1 -- docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md   # must print the Task 6 commit
RLOCK_TIMEOUT=600 bash tools/rlock.sh light /usr/bin/time -f "%e s %M KB" \
    python3 bench/soto_subset/run_subset.py guided --out /mnt/linuxdisk/tmp/o1_soto_subset/guided
```

Expected: about 30 s and under 1 GB. The output directory then holds `inputs.json`, `t1.json`, `violations.tsv`,
`t2.json`, `recipe.json`, and the recipe's `g1793.tsv`, `g2334.tsv`, `cn_assembly.tsv` and `recipe_*.tsv`.

A `SystemExit` naming a recipe rung or a coordinate disagreement is a failed wiring check. In that case write
nothing into the prereg; stop and report the message.

- [ ] **Step 2: Append the guided results to the prereg**

Under `## Results`, append the section below. Fill every cell from the JSON files, with the key named in each cell:
- ARI, Jaccard, H and D to 4 decimals;
- p to 6 decimals, followed by the count form, e.g. `0.000999 (1/1001)`;
- counts as integers;
- PASS when the bar's condition holds, else FAIL.

```markdown
### Guided (2026-10-14)

Run: commit `<git rev-parse --short HEAD>`; wall `<s>` s, peak `<KB>` KB. Wiring checks: all passed.

**T1**

| Quantity | Value |
|---|---|
| Soto families with ≥ 2 members (`families_ge2`) | |
| contained (`contained`); rate over 491 (`rate_all`); over ≥ 2 (`rate_ge2`) | |
| size null mean / max (`null.mean`, `null.max`) | |
| p (`p`) → **T1a** | |
| H(Soto \| ours) / H(ours \| Soto) | |
| D [95% interval] (`D`, `D_ci95`) → **T1b** | |
| ≤ 3 members: contained / families (`by_size.le3`) | |
| > 3 members: contained / families (`by_size.gt3`) | |
| plan 1's rate: contained / evaluable (`plan1`) | |
| violations: coverage / cover / conflict (`violation_classes`) | |

The conflict families (from `violations.tsv`, `class` = conflict), as `label (size, blocks)`: <list, comma-separated>.

**T2 (common set, `common_members` members)**

| Copy number | exact | ARI | ARI w/o FAM90A | Hungarian J | null exact mean (max) | null ARI mean (max) | p exact | p ARI | bar |
|---|---|---|---|---|---|---|---|---|---|
| assembly (T2's bar) | | | | | | | | | PASS/FAIL |
| read depth | | | | | | | | | (reported) |
| S1C (circular upper bound) | | | | | | | | | (reported) |
| our blocks alone (`ours_alone`) | | | | | — | — | — | — | — |
| ceiling: Soto ∩ ours (`ceiling`) | | | | | — | — | — | — | — |

**T2, each source on its own coverage (`own`)**

| Copy number | members | exact | ARI | ARI w/o FAM90A | Hungarian J |
|---|---|---|---|---|---|
| assembly | | | | | |
| read depth | | | | | |
| S1C | | | | | |

**Soto's recipe with each source (`recipe.json`)**

| Copy number | full gene set ARI / exact of 491 | common set exact | common set ARI | common set ARI w/o FAM90A | common set Hungarian J |
|---|---|---|---|---|---|
| none | | | | | |
| S1C | | | | | |
| read depth | | | | | |
| assembly | | | | | |

Guided verdicts: T1a <PASS/FAIL>, T1b <PASS/FAIL>, T2 <PASS/FAIL>; read depth on the common set <passes / does not
pass> T2's bar.
```

- [ ] **Step 3: Commit**

```bash
git add docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md
git commit -m "prereg results: Soto subset, guided universe

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z" -- docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md
```

---

### Task 8: the de novo replication runs and their results

**Files:**
- Modify: `docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md` (append under `## Results`)
- Outputs (not committed): `/mnt/linuxdisk/tmp/o1_soto_subset/human_A119b/`, `/mnt/linuxdisk/tmp/o1_soto_subset/human_testis/`

**Interfaces:**
- Consumes: `t1.json`, `t2.json` and `violations.tsv` of both catalogs. There is no `recipe.json` for de novo (R4).
- Produces: the `### De novo (2026-10-14)` results section, which Task 9 reads.

- [ ] **Step 1: Run both catalogs, one after the other**

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_m2_soto
for s in human_A119b human_testis; do
  RLOCK_TIMEOUT=600 bash tools/rlock.sh light /usr/bin/time -f "%e s %M KB" \
      python3 bench/soto_subset/run_subset.py denovo --sample $s --out /mnt/linuxdisk/tmp/o1_soto_subset/$s
done
```

Expected: each takes under a minute. A `SystemExit` is a failed wiring check: write nothing; stop and report it.

- [ ] **Step 2: Append the de novo results to the prereg**

Append the section below with the same formatting rules as Task 7. The rows are the same as the guided T1 and T2
tables, with one column per catalog. The recipe rows are not run (R4).

```markdown
### De novo (2026-10-14)

Run: commit `<git rev-parse --short HEAD>`; A119b wall `<s>` s, peak `<KB>` KB; testis wall `<s>` s, peak `<KB>` KB.
Wiring checks: all passed.

**T1**

| Quantity | A119b | testis |
|---|---|---|
| Soto families with ≥ 2 members | | |
| contained; rate over 491; over ≥ 2 | | |
| size null mean / max | | |
| p → **T1a** | | |
| H(Soto \| ours) / H(ours \| Soto) | | |
| D [95% interval] → **T1b** | | |
| ≤ 3 members: contained / families | | |
| > 3 members: contained / families | | |
| plan 1's rate: contained / evaluable | | |
| violations: coverage / cover / conflict | | |

Conflict families: A119b <list>; testis <list>.

**T2 (common set: A119b `common_members`, testis `common_members`)**

| Copy number | catalog | exact | ARI | ARI w/o FAM90A | Hungarian J | null exact mean (max) | null ARI mean (max) | p exact | p ARI | bar |
|---|---|---|---|---|---|---|---|---|---|---|
| assembly (T2's bar) | A119b | | | | | | | | | PASS/FAIL |
| assembly (T2's bar) | testis | | | | | | | | | PASS/FAIL |
| read depth | A119b | | | | | | | | | (reported) |
| read depth | testis | | | | | | | | | (reported) |
| S1C | A119b | | | | | | | | | (reported) |
| S1C | testis | | | | | | | | | (reported) |
| our blocks alone | A119b | | | | | — | — | — | — | — |
| our blocks alone | testis | | | | | — | — | — | — | — |
| ceiling | A119b | | | | | — | — | — | — | — |
| ceiling | testis | | | | | — | — | — | — | — |

**T2, each source on its own coverage**

| Copy number | catalog | members | exact | ARI | ARI w/o FAM90A | Hungarian J |
|---|---|---|---|---|---|---|
| assembly | A119b | | | | | |
| assembly | testis | | | | | |
| read depth | A119b | | | | | |
| read depth | testis | | | | | |
| S1C | A119b | | | | | |
| S1C | testis | | | | | |

De novo verdicts: A119b T1a <>, T1b <>, T2 <>; testis T1a <>, T1b <>, T2 <>.
```

- [ ] **Step 3: Commit**

```bash
git add docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md
git commit -m "prereg results: Soto subset, de novo A119b and testis

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z" -- docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md
```

---

### Task 9: the readings, register rows and thesis line

**Files:**
- Modify: `docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md` (append `### Verdict`)
- Modify: `docs/NEGATIVE_RESULTS_REGISTER.md` (append two rows at the end of the table)
- Modify: `docs/THESIS_OBJECTIVES.md` (line 31, the O1 row: one sentence in its third cell)
- Modify: `docs/superpowers/specs/2026-10-14-soto-subset-cn-design.md` (append `## Results (2026-10-14)`)

**Interfaces:**
- Consumes: the verdict lines of Tasks 7 and 8.
- Produces: the final record. Nothing reads it afterwards.

- [ ] **Step 1: Pick the reading**

Use the prereg's readings table, top row first, and take the first row whose conditions hold:
- **Reading 1:** guided T1a, T1b and T2 PASS, and T1a, T1b and T2 PASS on both A119b and testis.
- **Reading 2:** guided T1a, T1b and T2 PASS, and any de novo bar FAILS.
- **Reading 3:** guided T1a and T1b PASS, and guided T2 FAILS. The copy-number part is claimed only if the guided
  read-depth row passes (`read_depth_passes`).
- **Reading 4:** guided T1a or T1b FAILS.

- [ ] **Step 2: Append the verdict to the prereg**

```markdown
### Verdict (2026-10-14)

Reading <1|2|3|4> of the readings table. <One sentence per bar: name, universe, the observed statistic against its
null or interval, PASS/FAIL.> <For reading 2: which catalog and which bar failed. For reading 3: whether read depth
passed T2 on the guided common set.>
```

- [ ] **Step 3: Append two register rows**

Find the last row number with `grep -o '^| 1[0-9][0-9][0-9] ' docs/NEGATIVE_RESULTS_REGISTER.md | tail -1`. Use the
next two numbers, N and N+1. Append both rows at the end of the table: after the last row, before any text that
follows it.

```markdown
| N | 2026-10-14 | O1: Soto ⊂ ours (T1, `docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md`) | Soto's families lie inside our families more often than random partitions with our block sizes allow, and ours is the coarser partition (guided; de novo A119b and testis as replication) | <✅ PASS or ⛔ FAIL, in bold, with the reading number>: guided contained <c>/<f> (null mean <m>, p <p>), D <D> [<lo>, <hi>]; A119b <c>/<f> (p <p>), D <D> [<lo>, <hi>]; testis <c>/<f> (p <p>), D <D> [<lo>, <hi>]; violations guided coverage/cover/conflict <a>/<b>/<c> |
| N+1 | 2026-10-14 | O1: Soto ≈ ours ∩ CN (T2, same prereg) | Cutting our families by assembly copy number with Soto's pair rule (|Δ| < 2) rebuilds Soto's families beyond a within-family shuffle of the copy numbers | <✅ PASS or ⛔ FAIL, in bold>: guided common set (<n> genes) assembly exact <e> (null mean <m>, p <p>), ARI <a> (w/o FAM90A <a2>; p <p>); read depth exact <e>, ARI <a>; S1C exact <e>, ARI <a>; ours alone exact <e>, ARI <a>; ceiling exact <e>, ARI <a>; A119b assembly exact <e> (p <p>), ARI <a> (p <p>); testis assembly exact <e> (p <p>), ARI <a> (p <p>) |
```

The `<…>` slots are filled from the Results tables, with the same formatting.

- [ ] **Step 4: Update the O1 row of `docs/THESIS_OBJECTIVES.md`**

Line 31's third cell ends with `E1b also failed)`, followed by ` | Overall bipartite F`. Insert `; ` and the one
sentence of the chosen reading right after `E1b also failed)`. Leave everything else on the line unchanged. The
`<…>` slots are filled from the guided T2 assembly common-set row, the ceiling row and the T1 table.

- Reading 1: `Soto as a refinement (2026-10-14, docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md, registers N-N+1): Soto's families are a refinement of our multi-copy gene families, and cutting ours by assembly copy number rebuilds them (guided common set: <e> exact of <f>, ARI <a>; ceiling ARI <c>), on the guided universe and on both de novo catalogs`
- Reading 2: `Soto as a refinement (2026-10-14, docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md, registers N-N+1): on the guided universe only, Soto's families are a refinement of our multi-copy gene families and cutting ours by assembly copy number rebuilds them (common set: <e> exact of <f>, ARI <a>; ceiling ARI <c>); de novo <catalog> failed <bar>`
- Reading 3, read depth passes: `Soto ⊂ ours (2026-10-14, docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md, registers N-N+1): Soto's families lie inside ours (guided <c>/<f> contained, D <D> [<lo>, <hi>]); the copy-number cut is claimed with population read depth only (human-only: exact <e>, ARI <a>), since assembly copy number failed T2`
- Reading 3, read depth fails: `Soto ⊂ ours (2026-10-14, docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md, registers N-N+1): Soto's families lie inside ours (guided <c>/<f> contained, D <D> [<lo>, <hi>]); cutting ours by copy number is not claimed (T2 failed with assembly and read depth)`
- Reading 4: `Soto as a refinement: not claimed (2026-10-14, docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md, registers N-N+1; T1 failed: <T1a and/or T1b>; size and violation tables in the prereg)`

- [ ] **Step 5: Point the spec at the results**

Append to `docs/superpowers/specs/2026-10-14-soto-subset-cn-design.md`:

```markdown
## Results (2026-10-14)

Reading <n> of section 1's readings table: <the same one sentence as the thesis line>. Full results:
`docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md` (Results, Verdict); register rows N and N+1.
```

- [ ] **Step 6: Run the suites once more and commit**

```bash
RLOCK_TIMEOUT=300 bash tools/rlock.sh light python3 -m pytest bench/soto_subset bench/o1_hier bench/soto -q -p no:cacheprovider
git add docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md docs/NEGATIVE_RESULTS_REGISTER.md docs/THESIS_OBJECTIVES.md \
    docs/superpowers/specs/2026-10-14-soto-subset-cn-design.md
git commit -m "Soto subset: verdict (reading <n>), register rows N-N+1, O1 thesis line, spec results pointer

Co-Authored-By: <your model name> <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_015Xtq7GodhTGPYBzHDp3o4z" -- docs/archive/2026-10/merged/PREREG_soto_subset_cn_2026-10-14.md docs/NEGATIVE_RESULTS_REGISTER.md \
    docs/THESIS_OBJECTIVES.md docs/superpowers/specs/2026-10-14-soto-subset-cn-design.md
```

Expected: every test passes.
