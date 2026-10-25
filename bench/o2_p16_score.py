# P16 (LEDGER_O2): O2 arms D (default window) / T (--tie-structure) and minimap2's primary M against the simulated truth
import csv, re, subprocess, collections
truth = {r["read"]: r for r in csv.DictReader(open("truth.tsv"), delimiter="\t")}
prim_right = {}
for ln in subprocess.run(["samtools", "view", "-F", "0x904", "sim.bam"], capture_output=True, text=True, check=True).stdout.splitlines():
    f = ln.split("\t"); t = truth[f[0]]
    rl = sum(int(n) for n, o in re.findall(r"(\d+)([MDN=X])", f[5]))
    s = int(f[3]) - 1
    prim_right[f[0]] = f[2] == t["chrom"] and s < int(t["end"]) and s + rl > int(t["start"])
res = {}
for arm in ("D", "T"):
    rows = list(csv.DictReader(open(f"{arm}.o2.tsv"), delimiter="\t"))
    c = collections.Counter(); asg = {}
    for r in rows:
        c["tied rows"] += 1
        k = r["assigned"]
        if k in ("ABSTAIN", "HOPELESS"): c[k] += 1; continue
        asg[r["read"]] = k
        c["assigned"] += 1
        right = k == truth[r["read"]]["copy"]
        c["right"] += right; c["wrong"] += not right
        c["primary was right on assigned"] += prim_right[r["read"]]
        c["fixed (primary wrong, O2 right)"] += right and not prim_right[r["read"]]
        c["broken (primary right, O2 wrong)"] += (not right) and prim_right[r["read"]]
    final = sum((asg[n] == truth[n]["copy"]) if n in asg else prim_right[n] for n in truth)
    c["tied rows primary right"] = sum(prim_right[r["read"]] for r in rows)
    res[arm] = (c, final)
    print(f"== {arm}: " + "; ".join(f"{k} {v}" for k, v in c.items()) + f"; final right {final}/{len(truth)}")
M = sum(prim_right.values())
print(f"== M (primary): right {M}/{len(truth)}")
cT, fT = res["T"]; cD, fD = res["D"]
p = cT["right"] / cT["assigned"] if cT["assigned"] else float("nan")
print(f"H1 T precision {cT['right']}/{cT['assigned']} = {p:.4f} (>= .98): {'PASS' if p >= .98 else 'FAIL'}")
print(f"H2 T wrong {cT['wrong']} <= D wrong {cD['wrong']}: {'PASS' if cT['wrong'] <= cD['wrong'] else 'FAIL'}")
print(f"H3 final T {fT} >= D {fD} and >= M {M}: {'PASS' if fT >= fD and fT >= M else 'FAIL'}")
