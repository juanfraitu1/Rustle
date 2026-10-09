import json, os, subprocess, sys, collections
sys.path.insert(0, "/mnt/linuxdisk/home/juanfraitu/rustle_m2_soto/bench/unmapped_rescue"); sys.path.insert(0, "/mnt/linuxdisk/home/juanfraitu/rustle_m2_soto/bench")
import discover as D, seeds as SD, run_polish as RP
from o3_maternal import common as C
import pysam
O = "/mnt/linuxdisk/tmp/o3_rescue/mattruth"
e = json.load(open(f"{O}/eval.json")); v = {x["k"]: x for x in e["verdicts"]}
cons = {n.split("|")[0]: s for n, s in RP.read_cons(f"{O}/net_run/cons.fa").items()}
SD.write_fa(f"{O}/flag_cons.fa", {k: cons[k] for k in v}, sorted(v))
al = C.alias(); gate = json.load(open(f"{D.W}/o3hap_mat/gate.json"))["gate"]
sc = {}
for hap in ("mat", "pat"):
    paf = f"{O}/flag_cons.{hap}.paf"
    if not os.path.exists(paf + ".done"):
        subprocess.run(f"minimap2 -c -x splice:hq -uf -N 5 -t 4 {C.HAP_IDX.format(hap)} {O}/flag_cons.fa > {paf}", shell=True, check=True, stderr=subprocess.DEVNULL)
        open(paf + ".done", "w").write("ok")
    recs = D.best_records(paf); fa = pysam.FastaFile(D.HAP_FA.format(hap))
    sc[hap] = {k: (D.idcov_with_rescue(cons[k], recs[k], fa, lambda n: C.accession(n, al), gate, f"{O}/tmp_41c") if k in recs else 0.0) for k in v}
out = collections.Counter()
for k, x in v.items():
    best = max(sc["mat"][k], sc["pat"][k])
    out[(x["verdict"], "consensus on a haplotype >= 0.99042" if best >= 1 - D.DELTA else ("0.9-0.99042" if best >= 0.9 else "< 0.9"))] += 1
    x["hap_best"] = best
for key in sorted(out): print(key, out[key])
w = [x for x in v.values() if x["verdict"] == "WRONG"]
import statistics as st
print("WRONG: median best haplotype score", round(st.median(x["hap_best"] for x in w), 4), "; median primary score", round(st.median(x["R"] for x in w), 4))
print("WRONG with haplotype >= 0.99042 by reads:", sum(x["reads"] for x in w if x["hap_best"] >= 1 - D.DELTA), "of", sum(x["reads"] for x in w))
json.dump(v, open(f"{O}/diag41c.json", "w"))
