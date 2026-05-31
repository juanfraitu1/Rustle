#!/usr/bin/env python3
"""Obj 4: read-to-copy assignment accuracy on the synthetic fixture (fingerprint-EM).
Runs rustle --vg --vg-no-hmm with RUSTLE_VG_FP_ATTR_TSV, maps read_name_hash->true_copy
(via the same FNV-1a hash) and placement_copy->locus (via the output GTF copy_id+span),
scores argmax-weight assignment vs truth.
Usage: python3 bench/multi_copy_eval/score_synth_assignment.py [--rustle target/release/rustle]
Prints JSON with overall accuracy + decisive/uncertain buckets."""
import sys, os, re, json, subprocess
from collections import defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIX = os.path.join(REPO, "test_data", "synthetic_family")
BAM = os.path.join(FIX, "reads_sorted.bam")
REF = os.path.join(FIX, "ref.fa")

def fnv1a64(s: bytes) -> int:
    h = 0xcbf29ce484222325
    for b in s:
        h ^= b
        h = (h * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h

def read_names_truth():
    # samtools read names -> true copy
    out = subprocess.run(["samtools","view",BAM], capture_output=True, text=True, check=True).stdout
    truth = {}
    for line in out.splitlines():
        name = line.split("\t",1)[0]
        if re.match(r"uniq_A[12]_", name) or re.match(r"multi_\d+$", name): c = "A"
        elif re.match(r"uniq_B1_", name) or re.match(r"multi_r_\d+$", name): c = "B"
        else: continue
        truth[fnv1a64(name.encode())] = c
    return truth

def copyid_to_locus(gtf):
    # map copy_id -> 'A'/'B' via the assembled transcript span
    m = {}
    for line in open(gtf):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "transcript": continue
        cid = re.search(r'copy_id "?(\d+)"?', f[8])
        if not cid: continue
        locus = "A" if int(f[3]) < 5000 else "B"
        m[int(cid.group(1))] = locus
    return m

def main():
    rustle = sys.argv[sys.argv.index("--rustle")+1] if "--rustle" in sys.argv else os.path.join(REPO,"target","release","rustle")
    attr = "/tmp/synth_assign.attr.tsv"; gtf = "/tmp/synth_assign.gtf"
    env = dict(os.environ, RUSTLE_VG_FP_ATTR_TSV=attr)
    subprocess.run([rustle,"--vg","--vg-no-hmm","--genome-fasta",REF,"-L",BAM,"-o",gtf],
                   env=env, stderr=subprocess.DEVNULL, check=True)
    truth = read_names_truth()
    cid_locus = copyid_to_locus(gtf)
    # group attr rows by read_name_hash; columns per header
    rows = defaultdict(list)
    header = None
    for line in open(attr):
        parts = line.rstrip("\n").split("\t")
        if header is None: header = parts; continue
        d = dict(zip(header, parts))
        rows[int(d["read_name_hash"])].append(d)
    total=correct=decisive=dec_correct=0
    for rnh, rs in rows.items():
        if rnh not in truth: continue
        best = max(rs, key=lambda d: float(d["final_weight"]))
        assigned = cid_locus.get(int(best["placement_copy"]))
        if assigned is None: continue
        total += 1
        ok = (assigned == truth[rnh])
        correct += ok
        if float(best.get("weight_gap", 0)) > 0.8:
            decisive += 1; dec_correct += ok
    acc = round(100.0*correct/total,1) if total else 0.0
    dacc = round(100.0*dec_correct/decisive,1) if decisive else 0.0
    print(json.dumps({"scored": total, "accuracy": acc, "decisive": decisive,
                      "decisive_accuracy": dacc, "uncertain": total-decisive}, indent=2))

if __name__ == "__main__":
    main()

