"""Shared helpers for the Phase 0 flow-recall diagnostic."""
import json, os, re, subprocess

ROOT = "/mnt/c/Users/jfris/Desktop/Rustle"
RG = f"{ROOT}/bench/copy_recovery_eval/results_genomewide"
PERCHROM = f"{RG}/perchrom"
CACHE = f"{ROOT}/bench/flow_recall_phase0/cache"
RUSTLE = f"{ROOT}/target/release/rustle"
ST = f"{ROOT}/tools/stringtie/stringtie"

def parse_gtf(path):
    """transcript_id -> {strand, chrom, exons, introns(tuple), span, splice_sites(set)}"""
    tx = {}
    if not os.path.exists(path):
        return tx
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m:
            continue
        d = tx.setdefault(m.group(1), {"strand": f[6], "chrom": f[0], "exons": []})
        d["exons"].append((int(f[3]), int(f[4])))
    for d in tx.values():
        d["exons"].sort()
        ex = d["exons"]
        d["introns"] = tuple((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))
        d["span"] = (ex[0][0], ex[-1][1]) if ex else (0, 0)
        d["splice_sites"] = {c for iv in d["introns"] for c in iv}
    return tx

def chain_str(introns):
    return ",".join(f"{a}-{b}" for a, b in introns)

def parse_chain_str(s):
    return tuple(tuple(int(x) for x in p.split("-")) for p in s.split(",") if "-" in p)

def read_parity(path, step):
    out = []
    if not os.path.exists(path):
        return out
    for line in open(path):
        try:
            o = json.loads(line)
        except ValueError:
            continue
        if o.get("step") == step:
            out.append(o)
    return out
