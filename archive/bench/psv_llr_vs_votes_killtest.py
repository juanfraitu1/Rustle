#!/usr/bin/env python3
"""
KILL-OR-VALIDATE (clean decomposition): where does an LLR copy-assignment engine actually beat
the production vote-count rule, and WHY? Three orthogonal effects, isolated on identical reads:

  GATE     : production VOTE uses min_psv=3 (a read needs >=3 concordant PSVs); the LLR engine
             uses n_decisive>=1 (a read is resolvable iff it spans >=1 distinguishing column).
             min_psv=3 discards EVERY read in a 1-2 PSV family -> the near-identical tail.
  SCORING  : votes (count matches) vs log-likelihood (sum per-column logl). With flat per-base
             error these are monotonically equivalent (LLR margin = vote_margin * ln(3(1-e)/e)),
             so at a fixed gate they should match. This arm confirms/refutes that.
  QV-INFO  : the ONLY thing votes structurally cannot do: weight a high-QV column above a low-QV
             one. Tested on IDENTICAL reads (same injected errors), feeding the LLR either the
             TRUE per-base error (calibrated QV) or a flat mean error (uninformative QV).

Four arms, same reads, ground-truth labels:
  VOTE_prod   min_psv=3, margin=1   (production psv_linkage.rs default)
  VOTE_gate1  min_psv=1, margin=1   (vote scoring, permissive gate -> isolates GATE effect)
  LLR_flat    tau, flat mean-e QV   (likelihood scoring, no QV info -> vs VOTE_gate1 = SCORING)
  LLR_cal     tau, true per-base e  (likelihood + calibrated QV -> vs LLR_flat = QV-INFO)
"""
import math, random, argparse, json
from collections import Counter

BASES = "ACGT"


def make_copies(n_copies, n_psv, L=2800, seed=0):
    """Shared backbone; copies differ at n_psv columns. Alleles chosen so EVERY pair of copies
    differs at >=1 column (genuinely distinguishable; no accidental K=0 collisions)."""
    rng = random.Random(seed)
    backbone = [rng.choice(BASES) for _ in range(L)]
    psv_pos = [ (L // (n_psv + 1)) * (i + 1) for i in range(n_psv) ] if n_psv > 0 else []
    copies = []
    for c in range(n_copies):
        s = list(backbone)
        for j, p in enumerate(psv_pos):
            # distinct per (copy,column) hash so copies are pairwise-distinguishable
            h = (c * 2_654_435_761 ^ (j + 1) * 40_503 ^ seed) & 0xffffffff
            s[p] = BASES[h % 4]
        copies.append("".join(s))
    # guarantee pairwise distinguishability at the chosen columns; if any pair collides on ALL
    # columns (possible for n_psv=1), nudge one allele.
    for a in range(n_copies):
        for b in range(a + 1, n_copies):
            if all(copies[a][p] == copies[b][p] for p in psv_pos) and psv_pos:
                p = psv_pos[0]
                lst = list(copies[b]); lst[p] = BASES[(BASES.index(copies[a][p]) + 1) % 4]
                copies[b] = "".join(lst)
    return copies, psv_pos


def sim_read(seq, err_hi, rng, low_frac=0.15, low_lo=0.03, low_hi=0.15, hi_lo=1e-4, hi_hi=2e-3):
    """Substitution-only read. Per base: with prob low_frac it's a LOW-QV base (err in [low_lo,
    low_hi]); else a HIGH-QV base (err in [hi_lo,hi_hi]). err_hi scales the low-QV band so we can
    sweep difficulty. Errors injected at the per-base prob -> QV is calibrated to truth.
    Returns (read_base_idx[], true_e[])."""
    out, eprob = [], []
    for ch in seq:
        if rng.random() < low_frac:
            e = rng.uniform(low_lo, low_hi) * err_hi
        else:
            e = rng.uniform(hi_lo, hi_hi)
        e = min(e, 0.45)
        b = BASES.index(ch)
        if rng.random() < e:
            b = BASES.index(rng.choice([x for x in BASES if x != ch]))
        out.append(b); eprob.append(e)
    return out, eprob


def cols_at(read_bases, eprob, psv_pos, copies):
    out = []
    for p in psv_pos:
        out.append((read_bases[p], eprob[p], [BASES.index(copies[c][p]) for c in range(len(copies))]))
    return out


def assign_vote(cols, min_psv, margin):
    n = len(cols[0][2]) if cols else 0
    votes = [0] * n
    for rb, e, al in cols:
        for c in range(n):
            if al[c] == rb: votes[c] += 1
    if not votes: return None
    order = sorted(range(n), key=lambda c: votes[c], reverse=True)
    best = order[0]; second = order[1] if n > 1 else None
    bv = votes[best]; sv = votes[second] if second is not None else 0
    return best if (bv >= min_psv and bv - sv >= margin) else None


def assign_llr(cols, tau, calibrated, flat_e):
    n = len(cols[0][2]) if cols else 0
    logl = [0.0] * n; n_dec = 0
    for rb, e_true, al in cols:
        e = e_true if calibrated else flat_e
        e = min(max(e, 1e-7), 0.5)
        lm = math.log(1 - e); lx = math.log(e / 3)
        if len(set(al)) > 1: n_dec += 1
        for c in range(n):
            logl[c] += lm if al[c] == rb else lx
    if n_dec < 1: return None
    order = sorted(range(n), key=lambda c: logl[c], reverse=True)
    best = order[0]; second = order[1] if n > 1 else None
    margin = logl[best] - (logl[second] if second is not None else -1e9)
    return best if margin >= tau else None


def run(n_copies, n_psv, err_hi, reads_per_copy, tau, flat_e, seed=0):
    copies, psv_pos = make_copies(n_copies, n_psv, seed=seed)
    arms = ["VOTE_prod", "VOTE_gate1", "LLR_flat", "LLR_cal"]
    res = {a: Counter() for a in arms}
    for c in range(n_copies):
        for r in range(reads_per_copy):
            rng = random.Random((seed * 1_000_003) ^ (c * 2_654_435_761) ^ (r * 40_503))
            rb, ep = sim_read(copies[c], err_hi, rng)
            cols = cols_at(rb, ep, psv_pos, copies)
            calls = {
                "VOTE_prod":  assign_vote(cols, 3, 1),
                "VOTE_gate1": assign_vote(cols, 1, 1),
                "LLR_flat":   assign_llr(cols, tau, False, flat_e),
                "LLR_cal":    assign_llr(cols, tau, True, flat_e),
            }
            for a, call in calls.items():
                if call is None: res[a]["discard"] += 1
                elif call == c:  res[a]["correct"] += 1
                else:            res[a]["wrong"] += 1
    return res


def line(cnt):
    tot = sum(cnt.values()); asn = cnt["correct"] + cnt["wrong"]
    prec = cnt["correct"] / asn if asn else float("nan")
    return f"cor={cnt['correct']:4d} wr={cnt['wrong']:3d} dis={cnt['discard']:4d} prec={prec:.4f} rec={asn/tot:.3f}"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--reads", type=int, default=600)
    ap.add_argument("--tau", type=float, default=2.0)
    ap.add_argument("--flat_e", type=float, default=0.01, help="uninformative flat error fed to LLR_flat")
    args = ap.parse_args()
    print(f"KILL-OR-VALIDATE decomposition  (tau={args.tau} flat_e={args.flat_e} reads/copy={args.reads})")
    print("  GATE = VOTE_prod vs VOTE_gate1 | SCORING = VOTE_gate1 vs LLR_flat | QV = LLR_flat vs LLR_cal\n")
    grid = []
    for n_copies in (2, 5):
        for n_psv in (1, 2, 3, 5):
            for err_hi in (1.0, 3.0):
                res = run(n_copies, n_psv, err_hi, args.reads, args.tau, args.flat_e)
                tag = f"C={n_copies} K={n_psv} err_hi={err_hi}"
                print(tag)
                for a in ("VOTE_prod", "VOTE_gate1", "LLR_flat", "LLR_cal"):
                    print(f"    {a:11s} {line(res[a])}")
                # attribution deltas
                g = res["VOTE_gate1"]["correct"] - res["VOTE_prod"]["correct"]
                qv_wr = res["LLR_flat"]["wrong"] - res["LLR_cal"]["wrong"]
                print(f"    -> GATE recovers {g:+d} correct ; QV-info removes {qv_wr:+d} wrong (LLR_flat->LLR_cal)\n")
                grid.append(dict(n_copies=n_copies, n_psv=n_psv, err_hi=err_hi,
                                 **{a: dict(res[a]) for a in res}))
    json.dump(grid, open("bench/psv_llr_vs_votes_killtest.json", "w"), indent=1)
    print("wrote bench/psv_llr_vs_votes_killtest.json")
