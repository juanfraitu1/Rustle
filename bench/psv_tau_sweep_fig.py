#!/usr/bin/env python3
"""Figure: precision/recall vs the decisive margin tau, on (1) sim5x K-ladder TRUE labels and
(2) real GGO co-located families (unique-mapper agreement). Reads bench/psv_tau_sweep.json."""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = json.load(open("bench/psv_tau_sweep.json"))
taus = d["taus"]

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel A: sim K-ladder — recall & precision vs tau (true labels)
axA = axes[0]
colors = {1: "crimson", 2: "darkorange", 4: "seagreen", 8: "navy"}
for K, rec in d.get("sim", {}).items():
    cur = rec["curve"]
    t = [r["tau"] for r in cur]
    axA.plot(t, [r["recall"] for r in cur], "-o", color=colors.get(int(K), "gray"),
             label=f"K={K} recall", lw=2, ms=4)
    axA.plot(t, [r["precision"] if not np.isnan(r["precision"]) else np.nan for r in cur],
             "--s", color=colors.get(int(K), "gray"), lw=1, ms=3, alpha=.7)
axA.axvline(2.0, color="gray", ls=":", alpha=.6); axA.axvline(6.9, color="purple", ls=":", alpha=.6)
axA.text(2.0, .02, " τ=2 (default)", fontsize=8, color="gray")
axA.text(6.9, .02, " τ=6.9 (p=1e-3)", fontsize=8, color="purple")
axA.set_xlabel("decisive margin τ"); axA.set_ylabel("recall (solid) / precision (dashed)")
axA.set_title("sim5x K-ladder — TRUE labels\nK≥2 trivially perfect; K=1 is the tradeoff")
axA.legend(fontsize=7, loc="center right"); axA.grid(alpha=.3); axA.set_ylim(-.03, 1.05)

# Panel B: real GGO — recall + unique-mapper agreement vs tau
axB = axes[1]
g = d.get("ggo")
if g:
    cur = g["curve"]; t = [r["tau"] for r in cur]
    axB.plot(t, [r["recall"] for r in cur], "-o", color="teal", label="recall (assigned/all)", lw=2)
    axB.plot(t, [r["uniq_agreement"] if r["uniq_agreement"] is not None and not (isinstance(r["uniq_agreement"],float) and np.isnan(r["uniq_agreement"])) else np.nan for r in cur],
             "-s", color="crimson", label="unique-mapper agreement", lw=2)
    axB.axvline(2.0, color="gray", ls=":", alpha=.6); axB.axvline(6.9, color="purple", ls=":", alpha=.6)
    axB.set_title(f"real GGO families (n={g['n_reads']} reads / {g['n_families']} fam)\nunique-mapper agreement vs recall as τ rises")
    axB.legend(fontsize=8)
else:
    axB.text(.5, .5, "GGO arm not in json yet", ha="center")
axB.set_xlabel("decisive margin τ"); axB.set_ylabel("fraction"); axB.grid(alpha=.3); axB.set_ylim(0, 1.05)

# Panel C: the operating-point frontier (real GGO: agreement vs recall, tau annotated)
axC = axes[2]
if g:
    cur = g["curve"]
    xs = [r["recall"] for r in cur]
    ys = [r["uniq_agreement"] for r in cur]
    axC.plot(xs, ys, "-o", color="purple")
    for r in cur:
        if r["tau"] in (0.0, 2.0, 6.9, 12.0):
            axC.annotate(f"τ={r['tau']}", (r["recall"], r["uniq_agreement"]),
                         fontsize=8, xytext=(4, 4), textcoords="offset points")
    axC.set_xlabel("recall (fraction assigned)")
    axC.set_ylabel("unique-mapper agreement (precision proxy)")
    axC.set_title("operating-point frontier (real GGO)\npick τ from this curve")
    axC.grid(alpha=.3)
else:
    axC.text(.5, .5, "needs GGO arm", ha="center")

plt.tight_layout(); plt.savefig("bench/psv_tau_sweep.png", dpi=130)
print("wrote bench/psv_tau_sweep.png")
