#!/usr/bin/env python3
"""Quantification benchmark: does resolving multi-mapper ambiguity improve PER-COPY abundance estimation,
and how does it scale with PSV density (the identifiability knob) and abundance skew?

Two near-identical paralog copies (A, B) at true abundance (theta_A, 1-theta_A). They differ at K PSV
positions (copy A allele = 0, copy B allele = 1 at each). A read from copy c covers each PSV w.p. P_COV and
observes c's allele, flipped w.p. EPS (HiFi error). We estimate theta_A three ways and measure |theta_hat - theta|:

  naive   : 1/k split -> theta_hat = 0.5 (no use of read evidence; what dropping/uniform assignment gives).
  em      : abundance EM using each read's PSV likelihood (the standard quantifier, soft).
  psv_hard: per-molecule majority vote over a read's PSV alleles; ties -> 0.5 (our explicit assignment).

Sweep K and theta. Deterministic (seeded). Writes quant_benchmark.tsv + quant_benchmark.png.
Run with /home/juanfra/miniforge3/bin/python (numpy/matplotlib).
"""
import math
import random
import statistics as st

N = 2000          # reads at the locus
P_COV = 0.6       # prob a read covers a given PSV (read length / locus geometry)
EPS = 0.005       # per-base error (HiFi)
KS = [0, 1, 2, 4, 8, 16]
THETAS = [0.5, 0.7, 0.9]
REPS = 25


def sim(n, k, theta, seed):
    rng = random.Random(seed)
    reads = []  # each read = list of observed alleles (0/1) at the PSVs it covers
    for _ in range(n):
        c = 0 if rng.random() < theta else 1            # true copy: 0=A, 1=B (its allele at every PSV = c)
        obs = []
        for _psv in range(k):
            if rng.random() < P_COV:
                al = c if rng.random() > EPS else 1 - c  # observe copy's allele, flipped w.p. EPS
                obs.append(al)
        reads.append(obs)
    return reads


def est_naive(_reads):
    return 0.5


def est_psv_hard(reads):
    n_a = tie = 0
    for obs in reads:
        a = sum(1 for x in obs if x == 0)  # votes for A
        b = len(obs) - a                   # votes for B
        if a > b:
            n_a += 1
        elif a == b:                       # no covered PSV, or a tie
            tie += 1
    return (n_a + 0.5 * tie) / len(reads)


def est_em(reads):
    lp_match, lp_mis = math.log(1 - EPS), math.log(EPS)
    theta = 0.5
    for _ in range(100):
        s = 0.0
        for obs in reads:
            # log P(read | A) with A's allele = 0 everywhere; | B with allele = 1.
            la = sum(lp_match if x == 0 else lp_mis for x in obs)
            lb = sum(lp_match if x == 1 else lp_mis for x in obs)
            la += math.log(theta)
            lb += math.log(1 - theta)
            m = max(la, lb)
            ga = math.exp(la - m) / (math.exp(la - m) + math.exp(lb - m))
            s += ga
        new = s / len(reads)
        if abs(new - theta) < 1e-9:
            theta = new
            break
        theta = new
    return theta


def main():
    rows = []
    print(f"{'theta':>6} {'K':>3} | {'naive':>16} {'em':>16} {'psv_hard':>16}   value=naive-psv")
    for theta in THETAS:
        for k in KS:
            errs = {"naive": [], "em": [], "psv_hard": []}
            for r in range(REPS):
                reads = sim(N, k, theta, seed=1000 * int(theta * 100) + 50 * k + r)
                errs["naive"].append(abs(est_naive(reads) - theta))
                errs["em"].append(abs(est_em(reads) - theta))
                errs["psv_hard"].append(abs(est_psv_hard(reads) - theta))
            mean = {m: st.mean(v) for m, v in errs.items()}
            sd = {m: (st.stdev(v) if len(v) > 1 else 0.0) for m, v in errs.items()}
            value = mean["naive"] - mean["psv_hard"]
            print(f"{theta:>6.2f} {k:>3} | {mean['naive']:>7.3f}±{sd['naive']:.3f}  "
                  f"{mean['em']:>7.3f}±{sd['em']:.3f}  {mean['psv_hard']:>7.3f}±{sd['psv_hard']:.3f}   {value:+.3f}")
            rows.append((theta, k, mean, sd, value))
    with open("quant_benchmark.tsv", "w") as fh:
        fh.write("theta\tK\tnaive_err\tem_err\tpsv_hard_err\tnaive_sd\tem_sd\tpsv_hard_sd\tvalue\n")
        for theta, k, mean, sd, value in rows:
            fh.write(f"{theta}\t{k}\t{mean['naive']:.4f}\t{mean['em']:.4f}\t{mean['psv_hard']:.4f}\t"
                     f"{sd['naive']:.4f}\t{sd['em']:.4f}\t{sd['psv_hard']:.4f}\t{value:.4f}\n")
    print("\n[wrote quant_benchmark.tsv]")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, len(THETAS), figsize=(13, 4), sharey=True)
        for ax, theta in zip(axes, THETAS):
            sub = [r for r in rows if r[0] == theta]
            ks = [r[1] for r in sub]
            for m, color in (("naive", "#d62728"), ("em", "#1f77b4"), ("psv_hard", "#2ca02c")):
                ax.errorbar(ks, [r[2][m] for r in sub], yerr=[r[3][m] for r in sub],
                            marker="o", capsize=3, label=m, color=color)
            ax.set_title(f"abundance {int(theta*100)}:{int((1-theta)*100)}")
            ax.set_xlabel("PSVs between copies (K)")
            ax.grid(alpha=0.3)
        axes[0].set_ylabel("|theta_hat - theta|  (abundance error)")
        axes[0].legend()
        fig.suptitle("Per-copy quantification: value of multi-mapper resolution vs PSV density & abundance skew")
        fig.tight_layout()
        fig.savefig("quant_benchmark.png", dpi=130)
        print("[wrote quant_benchmark.png]")
    except Exception as e:
        print(f"(plot skipped: {e})")


if __name__ == "__main__":
    main()
