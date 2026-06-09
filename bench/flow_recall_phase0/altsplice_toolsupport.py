"""Corrected altsplice rule sweep: use rustle's OWN junction support (junction_accept `mm`,
which matched ST's nreads in the instrumented trace) instead of raw BAM CIGAR counts.

For each altsplice near-miss where ST='=' and rustle misses (1 alt-donor/acceptor shift),
compare the support rustle assigns to the ANNOTATED junction vs the junction rustle's flow
actually emitted. If the annotated junction has HIGHER rustle-support, then "prefer higher
good-support at the shared boundary" is a viable both-axes rule rustle has the data for.
"""
import os, sys, glob, re, json, subprocess
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

CHROMS = ["NC_073227.2", "NC_073231.2", "NC_073235.2", "NC_073239.2"]

def chains(path):
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
        out[k] = tuple((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))
    return out

def rustle_support(chrom):
    """intron-interior (a,b) -> rustle junction_accept mm (support)."""
    wd = f"/tmp/altsup/{chrom}"; os.makedirs(wd, exist_ok=True)
    bam = f"{lib.PERCHROM}/{chrom}/c.bam"
    env = dict(os.environ, RUSTLE_PARITY_LOG=f"{wd}/j.jsonl",
               RUSTLE_PARITY_FILTER_STEPS="junction_accept")
    subprocess.run([lib.RUSTLE, "-L", bam, "-o", f"{wd}/r.gtf"], env=env,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    sup = {}
    for line in open(f"{wd}/j.jsonl"):
        try:
            o = json.loads(line)
        except ValueError:
            continue
        if o.get("step") != "junction_accept":
            continue
        a, b = o["start"] + 1, o["end"] - 1   # flank -> intron-interior
        mm = o["payload"].get("mm", 0.0)
        sup[(a, b)] = max(sup.get((a, b), -1.0), mm)
    return sup

def main():
    res = Counter()
    detail = []
    for c in CHROMS:
        wd = f"/tmp/oracle_nm/{c}"
        tmaps = glob.glob(f"{wd}/g.*.tmap")
        if not tmaps:
            continue
        st_fsm = set(open(f"{lib.PERCHROM}/{c}/s_fsm.txt").read().split())
        refc = chains(f"{lib.PERCHROM}/{c}/ref.gtf")
        ruc = chains(f"/tmp/strand_safety/all/{c}/off.gtf")
        sup = rustle_support(c)
        hits = {}; matched = set()
        for line in open(tmaps[0]):
            f = line.rstrip("\n").split("\t")
            if f[0] == "ref_gene_id":
                continue
            if f[2] == "=":
                matched.add(f[1])
            elif f[2] in ("j", "k", "c"):
                hits.setdefault(f[1], []).append(f[4])
        for rid, qs in hits.items():
            if rid in matched or rid not in refc or rid not in st_fsm:
                continue   # ST converts, rustle doesn't
            R = set(refc[rid]); best = None; bs = -1
            for q in qs:
                Q = ruc.get(q)
                if Q and len(R & set(Q)) > bs:
                    bs = len(R & set(Q)); best = set(Q)
            if not best or bs <= 0:
                continue
            re_, qe = R - best, best - R
            if len(re_) == 1 and len(qe) == 1:
                ra, qa = next(iter(re_)), next(iter(qe))
                if ra[0] == qa[0] or ra[1] == qa[1]:   # shared-endpoint alt-donor/acceptor
                    sr = sup.get(ra, None)   # annotated junction support in rustle
                    sq = sup.get(qa, None)   # rustle's emitted junction support
                    if sr is None or sq is None:
                        res["missing_support"] += 1
                        continue
                    if sr > sq:
                        res["ref_higher (rule CORRECT)"] += 1
                    elif sr < sq:
                        res["rustle_higher (rule wrong)"] += 1
                    else:
                        res["tie"] += 1
                    detail.append((rid, sr, sq))
    n = sum(v for k, v in res.items() if k != "missing_support")
    print(f"=== altsplice rule, RUSTLE good-junction support (ST-converts cases, {len(CHROMS)} chroms) ===")
    for k, v in res.most_common():
        pct = f"  ({100*v/max(1,n):.0f}%)" if k != "missing_support" else ""
        print(f"  {k:30s} {v}{pct}")
    corr = res["ref_higher (rule CORRECT)"]
    print(f"\n  => 'prefer higher rustle-support at the shared boundary' matches annotation "
          f"in {corr}/{n} ({100*corr/max(1,n):.0f}%) cases")
    print("  (raw-BAM-count sweep earlier said only 12% — wrong metric)")

if __name__ == "__main__":
    main()
