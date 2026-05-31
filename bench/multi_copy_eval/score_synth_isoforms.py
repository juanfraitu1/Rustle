#!/usr/bin/env python3
"""Obj 3: per-copy isoform Sn/Pr on the synthetic fixture vs truth.gtf.
Runs rustle --vg (and a non-vg baseline), partitions truth+pred by copy via gene span
(copy A < 5000, copy B > 9000), scores intron-chain Sn/Pr per copy.
Usage: python3 bench/multi_copy_eval/score_synth_isoforms.py [--rustle target/release/rustle]
Prints JSON: {"vg":{"A":{Sn,Pr},"B":{...},"overall":{...}}, "baseline":{...}}"""
import sys, os, re, json, subprocess

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIX = os.path.join(REPO, "test_data", "synthetic_family")
BAM = os.path.join(FIX, "reads_sorted.bam")
REF = os.path.join(FIX, "ref.fa")
TRUTH = os.path.join(FIX, "truth.gtf")

def chains_by_copy(path):
    tx = {}; strand = {}
    for line in open(path):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon": continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        tx.setdefault(m.group(1), []).append((int(f[3]), int(f[4])))
    out = {"A": set(), "B": set()}
    for t, ex in tx.items():
        ex.sort()
        ic = tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))
        if not ic: continue
        copy = "A" if ex[0][0] < 5000 else "B"
        out[copy].add(ic)
    return out

def run(rustle, vg):
    out_gtf = "/tmp/synth_iso_%s.gtf" % ("vg" if vg else "base")
    cmd = [rustle]
    if vg: cmd += ["--vg", "--genome-fasta", REF]
    cmd += ["-L", BAM, "-o", out_gtf]
    subprocess.run(cmd, stderr=subprocess.DEVNULL, check=True)
    return chains_by_copy(out_gtf)

def score(pred, truth):
    res = {}
    allp, allt = set(), set()
    for c in ("A", "B"):
        p, t = pred.get(c, set()), truth.get(c, set())
        tp = len(p & t)
        sn = 100.0*tp/len(t) if t else 0.0
        pr = 100.0*tp/len(p) if p else 0.0
        res[c] = {"Sn": round(sn,1), "Pr": round(pr,1), "tp": tp, "n_truth": len(t), "n_pred": len(p)}
        allp |= {(c,x) for x in p}; allt |= {(c,x) for x in t}
    tp = len(allp & allt)
    res["overall"] = {"Sn": round(100.0*tp/len(allt),1) if allt else 0.0,
                      "Pr": round(100.0*tp/len(allp),1) if allp else 0.0}
    return res

def main():
    rustle = sys.argv[sys.argv.index("--rustle")+1] if "--rustle" in sys.argv else os.path.join(REPO,"target","release","rustle")
    truth = chains_by_copy(TRUTH)
    out = {"vg": score(run(rustle, True), truth), "baseline": score(run(rustle, False), truth)}
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
