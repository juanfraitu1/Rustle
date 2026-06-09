"""Characterize the no_overlap fragmentation buckets (rustle_near + sparse-long) with the
CURRENT default binary (strand-bundling now on). Are these structural bundle-fragmentation
(gene split across pieces -> a strand-bundling-style merge could recover) or flow/coverage?

Per case, run default rustle -L on the locus slice (fresh, not cached) and classify:
  recovered_now    : ref chain now in rustle output (strand flip or otherwise fixed it)
  fragmented       : >=2 rustle tx overlap the ref span and TOGETHER cover all ref junctions,
                     but none matches alone -> bundle/assembly fragmentation (STRUCTURAL lever)
  flow_has_juncs   : all ref junctions accepted in rustle's graph, chain not extracted -> flow
  juncs_missing    : >=1 ref junction not accepted -> coverage/alignment (read support gap)
  emits_nothing    : no rustle tx overlaps the span at all
"""
import json, os, sys, subprocess, re
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

def run(slice_bam, chrom, outdir):
    os.makedirs(outdir, exist_ok=True)
    env = dict(os.environ)
    env["RUSTLE_PARITY_LOG"] = f"{outdir}/r.jsonl"
    env["RUSTLE_PARITY_FILTER_CHROM"] = chrom
    env["RUSTLE_PARITY_FILTER_STEPS"] = "junction_accept"
    subprocess.run([lib.RUSTLE, "-L", slice_bam, "-o", f"{outdir}/r.gtf"], env=env,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def gtf_tx(path):
    tx = {}
    for line in open(path):
        f = line.split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if m:
            tx.setdefault(m.group(1), []).append((int(f[3]), int(f[4])))
    out = {}
    for k, ex in tx.items():
        ex.sort()
        out[k] = {"span": (ex[0][0], ex[-1][1]),
                  "introns": set((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))}
    return out

def main():
    probed = [json.loads(l) for l in open(f"{lib.ROOT}/bench/copy_recovery_eval/"
                                          f"results_genomewide/no_overlap_probed.jsonl")]
    cases = [p for p in probed if p["bucket"] in ("rustle_near", "other")]
    wd = "/tmp/frag_probe"
    cls = Counter()
    rows = []
    for p in cases:
        chrom, tid = p["chrom"], p["tid"]
        ch = lib.ref_chain(chrom, tid)
        if ch is None:
            continue
        lo, hi, ref_introns = ch
        ref_introns = set(ref_introns)
        od = f"{wd}/{chrom}_{tid}"
        sb = f"{od}/s.bam"
        os.makedirs(od, exist_ok=True)
        with open(sb, "wb") as fh:
            subprocess.run(["samtools", "view", "-b", f"{lib.PERCHROM}/{chrom}/c.bam",
                            f"{chrom}:{max(1,lo-8000)}-{hi+8000}"], stdout=fh, stderr=subprocess.DEVNULL)
        subprocess.run(["samtools", "index", sb], stderr=subprocess.DEVNULL)
        run(sb, chrom, od)
        tx = gtf_tx(f"{od}/r.gtf")
        # accepted junctions (both conventions)
        acc = set()
        for j in lib.read_parity(f"{od}/r.jsonl", "junction_accept"):
            if j.get("payload", {}).get("accepted", False):
                a, b = j["start"], j["end"]
                acc.add((a, b)); acc.add((a + 1, b - 1))
        # overlapping rustle tx
        ov = [t for t in tx.values() if t["span"][0] <= hi and lo <= t["span"][1]]
        recovered = any(t["introns"] == ref_introns for t in tx.values())
        if recovered:
            c = "recovered_now"
        elif not ov:
            c = "emits_nothing"
        else:
            union = set().union(*(t["introns"] for t in ov)) if ov else set()
            if ref_introns and ref_introns <= union and len(ov) >= 2:
                c = "fragmented"
            elif ref_introns and ref_introns <= acc:
                c = "flow_has_juncs"
            else:
                c = "juncs_missing"
        cls[c] += 1
        rows.append({"chrom": chrom, "tid": tid, "bucket": p["bucket"], "class": c, "n_ov": len(ov)})

    n = sum(cls.values())
    print(f"=== fragmentation probe (rustle_near + sparse-long, n={n}, CURRENT default binary) ===")
    for k, v in cls.most_common():
        print(f"  {k:16s} {v}")
    struct = cls["fragmented"]
    print(f"\n  STRUCTURAL (fragmented, strand-bundling-style mergeable): {struct}")
    print(f"  recovered already (post strand-flip): {cls['recovered_now']}")
    with open(f"{lib.ROOT}/bench/flow_recall_phase0/fragmentation_probe.jsonl", "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")

if __name__ == "__main__":
    main()
