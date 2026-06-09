"""Oracle near-miss analysis: for annotated transcripts rustle ALMOST emits (a near-miss query
j/k/c but no exact '='), compute the structural edit that would convert rustle's chain to the
ref's exact chain. Converting a near-miss is BOTH-axes-positive (removes the j/k FP, adds the = TP).

Edit taxonomy (ref chain R vs rustle's best near-miss chain Q, sharing >=1 junction):
  k_trim     : Q ⊃ R (rustle longer; Q has R's junctions + extra)        -> DROP extra junction(s)
  c_extend   : Q ⊂ R (rustle shorter; missing junctions)                  -> ADD missing junction(s)
  altsplice  : 1 R-extra + 1 Q-extra sharing an endpoint (boundary shift) -> SNAP junction
  swap_multi : both have private junctions, not a clean shift             -> harder
Edit distance = |R-extra| + |Q-extra|. The "1-edit-away" set is the tractable rule target.

Runs on BASELINE -L output (the parity target). gffcompare vs per-chrom ref.
"""
import os, sys, subprocess, glob, re
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

REF = "/mnt/c/Users/jfris/Desktop/Rustle/bench/copy_recovery_eval/results/ref.gtf"
BASE = "/tmp/strand_safety/all"          # baseline -L off.gtf per chrom
CHROMS = ["NC_073227.2", "NC_073231.2", "NC_073235.2", "NC_073239.2"]  # mid-size sample

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

def edit_type(R, Q):
    R, Q = set(R), set(Q)
    if not (R & Q):
        return None  # no shared junction -> not a near-miss
    re_, qe = R - Q, Q - R
    if not qe:
        return ("c_extend", len(re_))
    if not re_:
        return ("k_trim", len(qe))
    if len(re_) == 1 and len(qe) == 1:
        (ra,), (qa,) = tuple(re_), tuple(qe)
        if ra[0] == qa[0] or ra[1] == qa[1]:
            return ("altsplice", 1)
        return ("swap_1", 2)
    return ("swap_multi", len(re_) + len(qe))

def main():
    cat = Counter(); dist = Counter(); by_dist1 = Counter()
    n_nearmiss = 0; n_matched = 0
    for c in CHROMS:
        bl = f"{BASE}/{c}/off.gtf"
        if not os.path.exists(bl):
            continue
        wd = f"/tmp/oracle_nm/{c}"; os.makedirs(wd, exist_ok=True)
        refc = f"{wd}/ref.gtf"
        with open(refc, "w") as o:
            subprocess.run(["awk", "-F\t", f'$1=="{c}"', REF], stdout=o)
        # gffcompare writes <prefix>.<input_basename>.tmap NEXT TO the input gtf; copy the
        # query into wd so the tmap lands here (not in the shared strand_safety dir).
        subprocess.run(["cp", bl, f"{wd}/q.gtf"])
        subprocess.run(["gffcompare", "-r", refc, "-o", "g", "q.gtf"],
                       cwd=wd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        tmap = glob.glob(f"{wd}/g.*.tmap")
        if not tmap:
            continue
        # ref -> set of (class, qry_id); also which refs got '='
        ref_hits = defaultdict(list); matched = set()
        for line in open(tmap[0]):
            f = line.rstrip("\n").split("\t")
            if f[0] == "ref_gene_id":
                continue
            ref_id, cls, qry = f[1], f[2], f[4]
            if cls == "=":
                matched.add(ref_id)
            elif cls in ("j", "k", "c", "m", "n"):
                ref_hits[ref_id].append((cls, qry))
        ref_ch = chains(f"{lib.PERCHROM}/{c}/ref.gtf")
        ru_ch = chains(bl)
        n_matched += len(matched)
        for ref_id, hits in ref_hits.items():
            if ref_id in matched or ref_id not in ref_ch:
                continue  # only near-miss-ONLY refs (no exact match)
            R = ref_ch[ref_id]
            if not R:
                continue
            # best near-miss query = max shared junctions
            best = None; bs = -1
            for _, qry in hits:
                Q = ru_ch.get(qry)
                if not Q:
                    continue
                sh = len(set(R) & set(Q))
                if sh > bs:
                    bs = sh; best = Q
            if best is None or bs <= 0:
                continue
            n_nearmiss += 1
            et = edit_type(R, best)
            if et is None:
                continue
            cat[et[0]] += 1; dist[min(et[1], 5)] += 1
            if et[1] == 1:
                by_dist1[et[0]] += 1

    print(f"=== oracle near-miss (BASELINE, {len(CHROMS)} chroms) ===")
    print(f"  matched '=' refs: {n_matched}")
    print(f"  near-miss-ONLY refs analyzed: {n_nearmiss}\n")
    print("edit TYPE (what converts rustle's near-miss to the exact ref):")
    for k, v in cat.most_common():
        print(f"  {k:12s} {v:4d}  ({100*v/max(1,n_nearmiss):4.0f}%)")
    print("\nedit DISTANCE (# junction changes):")
    for d in sorted(dist):
        lbl = f"{d}" if d < 5 else "5+"
        print(f"  {lbl} edit(s): {dist[d]}")
    print(f"\n1-edit-away set (the tractable rule target): {sum(by_dist1.values())}")
    for k, v in by_dist1.most_common():
        print(f"  {k:12s} {v}")

if __name__ == "__main__":
    main()
