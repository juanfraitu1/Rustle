"""Attribute the 109 VG regressions (baseline recovers, VG drops) to their VG-pipeline cause.

Per regression locus, slice the BAM and run configs, checking if the EXACT target ref chain
appears in each output:
  baseline   : rustle -L                                (control; defines the regression)
  full_vg    : --vg --vg-snp + TANDEM + DECISIVE_GATE   (the gw-protocol config that dropped it)
  no_decisive: full_vg minus RUSTLE_VG_DECISIVE_GATE
  plain_vg   : --vg --vg-snp only (no TANDEM, no DECISIVE)

Classify each:
  context_only : full_vg-on-slice RECOVERS it -> the drop was a whole-chromosome-context effect
                 (bundling/family scope), not a per-locus gate
  decisive_gate: no_decisive recovers but full_vg doesn't -> the decisive gate killed it (TUNABLE)
  tandem       : plain_vg recovers but full_vg doesn't (and no_decisive doesn't) -> tandem flags
  core_vg      : baseline recovers, NO vg config does -> --vg core mode (family graph/apportion)
  no_repro     : baseline-on-slice also lacks it -> slice doesn't reproduce baseline recovery
"""
import json, os, sys, subprocess
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

FASTA = "/home/juanfra/winloci_scratch/GGO.fasta"
WD = "/tmp/vg_attr_run"

def run(bam, outgtf, vg, extra_env):
    env = dict(os.environ); env.update(extra_env)
    cmd = [lib.RUSTLE]
    if vg:
        cmd += ["--vg", "--vg-snp", "--genome-fasta", FASTA]
    cmd += ["-L", bam, "-o", outgtf]
    subprocess.run(cmd, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def has_chain(gtf, introns):
    if not os.path.exists(gtf):
        return False
    import re
    tx = {}
    for line in open(gtf):
        f = line.split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if m:
            tx.setdefault(m.group(1), []).append((int(f[3]), int(f[4])))
    target = tuple(introns)
    for ex in tx.values():
        ex.sort()
        ic = tuple((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))
        if ic == target:
            return True
    return False

FULL = {"RUSTLE_VG_TANDEM": "1", "RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS": "1",
        "RUSTLE_VG_DECISIVE_GATE": "1", "RUSTLE_VG_DECISIVE_GATE_MIN_PRIM": "4"}
NO_DECISIVE = {"RUSTLE_VG_TANDEM": "1", "RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS": "1"}

def main():
    reg = json.load(open(f"{lib.ROOT}/bench/flow_recall_phase0/vg_regressions.json"))
    cls = Counter()
    rows = []
    for i, r in enumerate(reg):
        chrom, tid = r["chrom"], r["tid"]
        ch = lib.ref_chain(chrom, tid)
        if ch is None:
            continue
        lo, hi, introns = ch
        if not introns:
            cls["single_exon"] += 1; continue
        od = f"{WD}/{chrom}_{tid}"; os.makedirs(od, exist_ok=True)
        sb = f"{od}/s.bam"
        with open(sb, "wb") as fh:
            subprocess.run(["samtools", "view", "-b", f"{lib.PERCHROM}/{chrom}/c.bam",
                            f"{chrom}:{max(1,lo-8000)}-{hi+8000}"], stdout=fh, stderr=subprocess.DEVNULL)
        subprocess.run(["samtools", "index", sb], stderr=subprocess.DEVNULL)
        run(sb, f"{od}/bl.gtf", False, {})
        run(sb, f"{od}/full.gtf", True, FULL)
        bl = has_chain(f"{od}/bl.gtf", introns)
        full = has_chain(f"{od}/full.gtf", introns)
        if not bl:
            c = "no_repro"
        elif full:
            c = "context_only"
        else:
            run(sb, f"{od}/nodec.gtf", True, NO_DECISIVE)
            run(sb, f"{od}/plain.gtf", True, {})
            nodec = has_chain(f"{od}/nodec.gtf", introns)
            plain = has_chain(f"{od}/plain.gtf", introns)
            if nodec:
                c = "decisive_gate"
            elif plain:
                c = "tandem"
            else:
                c = "core_vg"
        cls[c] += 1
        rows.append({"chrom": chrom, "tid": tid, "cause": c})
        if (i + 1) % 20 == 0:
            print(f"{i+1}/{len(reg)} ...", flush=True)
    n = sum(cls.values())
    print(f"\n=== VG regression attribution (n={n}) ===")
    for k, v in cls.most_common():
        print(f"  {k:14s} {v}  ({100*v/n:.0f}%)")
    print("\n  TUNABLE (decisive_gate + tandem): %d" % (cls["decisive_gate"] + cls["tandem"]))
    print("  context_only (whole-chrom scope, not a gate): %d" % cls["context_only"])
    print("  core_vg (deep --vg mode): %d" % cls["core_vg"])
    with open(f"{lib.ROOT}/bench/flow_recall_phase0/vg_regression_attribution.jsonl", "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")

if __name__ == "__main__":
    main()
