#!/usr/bin/env python3
"""Real-reads validation of the soft per-copy quantification: subsample the sim5x HiFi reads (ground-truth
copy in each read name) to a SKEWED abundance, run copy_assign across the PSV ladder K=0..8, and compare the
soft-EM abundance vs the hard read-count vs the uniform prior to the imposed truth. Skewed (not 50:50) so the
uniform prior != truth -- this proves the EM RESOLVES, not just returns the prior. Deterministic.
Run with /home/juanfra/miniforge3/bin/python (matplotlib)."""
import os
import subprocess

S = "/home/juanfra/winloci_scratch/sim5x"
BIN = "/mnt/c/Users/jfris/Desktop/Rustle/target/release/copy_assign"
SKEW = [40, 30, 20, 10, 10]                 # reads per copy c0..c4 (max 40 available each)
TOTAL = sum(SKEW)
TRUTH = [s / TOTAL for s in SKEW]           # imposed abundance, position order
KS = [0, 1, 2, 4, 8]


def run_k(k):
    bam = f"{S}/sim5x_K{k}.bam"
    fa = f"{S}/sim5x_K{k}.ref.fa"
    length = open(f"{fa}.fai").read().split("\t")[1].strip()
    names = [f"K{k}_c{c}_r{r}" for c, cnt in enumerate(SKEW) for r in range(cnt)]
    nf = f"/tmp/sqv_names_K{k}.txt"
    open(nf, "w").write("\n".join(names) + "\n")
    sub = f"/tmp/sqv_sub_K{k}.bam"
    subprocess.run(f"samtools view -N {nf} -b {bam} | samtools sort -o {sub} - && samtools index {sub}",
                   shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    out = f"/tmp/sqv_K{k}"
    subprocess.run(
        f"bash -c 'ulimit -v 12000000; {BIN} --bam {sub} --fasta {fa} --region SIM5X_K{k}:0-{length} "
        f"--out {out} --min-copies 2'",
        shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    soft, hard = [], []
    for line in open(f"{out}.quant.tsv").read().splitlines()[1:]:
        f_ = line.split("\t")
        soft.append(float(f_[3]))
        hard.append(int(f_[5]))
    hs = sum(hard) or 1
    hardf = [h / hs for h in hard]
    return soft, hardf


def l1(est):
    """L1 abundance error vs TRUTH (position-matched; detected copies are start-sorted, as are c0..c4)."""
    e = sum(abs((est[i] if i < len(est) else 0.0) - TRUTH[i]) for i in range(5))
    e += sum(est[i] for i in range(5, len(est)))  # any spurious extra copy
    return e


def main():
    print(f"imposed truth (c0..c4): {[round(t,3) for t in TRUTH]}  (skew {SKEW})")
    rows = []
    for k in KS:
        soft, hardf = run_k(k)
        naive = [1.0 / len(soft)] * len(soft) if soft else []
        es, eh, en = l1(soft), l1(hardf), l1(naive)
        rows.append((k, len(soft), es, eh, en))
        print(f"K={k}: copies={len(soft)} | soft_err={es:.3f}  hard_err={eh:.3f}  naive_err={en:.3f}  "
              f"| soft={[round(x,2) for x in soft]}")
    with open("sim5x_quant_validate.tsv", "w") as fh:
        fh.write("K\tn_copies\tsoft_err\thard_err\tnaive_err\n")
        for k, n, es, eh, en in rows:
            fh.write(f"{k}\t{n}\t{es:.4f}\t{eh:.4f}\t{en:.4f}\n")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        ks = [r[0] for r in rows]
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot(ks, [r[2] for r in rows], "o-", color="#2ca02c", label="soft EM (our quant)")
        ax.plot(ks, [r[3] for r in rows], "s-", color="#d62728", label="hard read-count")
        ax.plot(ks, [r[4] for r in rows], "^--", color="#7f7f7f", label="naive (uniform prior)")
        ax.set_xlabel("PSVs per copy (K) — identifiability")
        ax.set_ylabel("L1 abundance error vs truth")
        ax.set_title(f"Real sim5x HiFi reads, skewed abundance {SKEW}: value of resolution for quantification")
        ax.grid(alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig("sim5x_quant_validate.png", dpi=130)
        print("[wrote sim5x_quant_validate.png]")
    except Exception as e:
        print(f"(plot skipped: {e})")


if __name__ == "__main__":
    main()
