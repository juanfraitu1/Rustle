"""Can ANNOTATION-FREE intrinsic features sharpen the multi-copy-family definition so lncRNA + noise
are excluded? Label de-novo loci by their gene biotype (the ground truth) and test whether ORF
(coding potential) + 4-mer complexity separate protein_coding from lncRNA from low-complexity. If they
do, the refined definition = ~R ∩ ~B  AND  ORF>=T_orf  AND  complexity>=T_cplx  (+ copy-count bound)."""
import sys, json, random
from itertools import product
sys.path.insert(0, "bench")
import copy_assign as CA

COMP = str.maketrans("ACGTN", "TGCAN")
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {"".join(c): aas[i] for i, c in enumerate(product("TCAG", repeat=3))}
def orf_aa(s):
    best = 0
    for strand in (s, s.translate(COMP)[::-1]):
        for f in range(3):
            prot = "".join(CODON.get(strand[i:i+3], "X") for i in range(f, len(strand)-2, 3))
            for p in prot.split("*"):
                k = p.find("M")
                if k >= 0: best = max(best, len(p)-k)
    return best
def cplx4(s):
    if len(s) < 8: return 0.0
    return len({s[i:i+4] for i in range(len(s)-3)}) / min(256, len(s)-3)

genes = json.load(open("/home/juanfra/winloci_scratch/gene_biotypes.json"))
import bisect
gstart = {c: [g[0] for g in v] for c, v in genes.items()}
def biotype(c, s, e):
    if c not in genes: return "unannotated"
    iv = genes[c]; st = gstart[c]; i = bisect.bisect_right(st, e) - 1
    for j in range(max(0, i-2), min(len(iv), i+2)):
        a, b, bt = iv[j]
        if s <= b and a <= e: return bt
    return "unannotated"

meta = CA.load_meta(); skel = CA.load_skel(); spliced = CA.load_spliced()
loci = json.load(open("/home/juanfra/winloci_scratch/refabsent/genomewide_loci.json"))
# representative (longest) tx per locus + biotype
random.Random(0).shuffle(loci)
by_bt = {}
for c, s, e in loci:
    bt = biotype(c, s, e)
    by_bt.setdefault(bt, [])
    if len(by_bt[bt]) >= 350: continue
    # longest tx in this locus
    best = None
    for tid, (tc, ts, te, st, nr) in meta.items():
        if tc == c and not (te < s or ts > e):
            sq = spliced.get(tid, "")
            if sq and (best is None or len(sq) > len(best[1])): best = (tid, sq)
    if best:
        tid, sq = best
        tc, ts, te, st, nr = meta[tid]
        ne = len(skel.get((tc, ts, te), [])) + 1
        by_bt[bt].append(dict(locus=f"{c}:{s}-{e}", orf=orf_aa(sq), cplx=round(cplx4(sq), 2),
                              n_exons=ne, txlen=len(sq), orf_frac=round(orf_aa(sq)*3/max(1, len(sq)), 2)))

import statistics as st
print(f"{'biotype':16s} {'n':>4} {'medORF':>6} {'%ORF>=100':>9} {'medCplx':>7} {'%cplx<0.4':>9} {'medExons':>8}")
KEEP_O, KEEP_C = 100, 0.40
for bt in ("protein_coding", "lncRNA", "pseudogene", "unannotated"):
    rs = by_bt.get(bt, [])
    if not rs: continue
    orfs = [r["orf"] for r in rs]; cx = [r["cplx"] for r in rs]
    print(f"{bt:16s} {len(rs):>4} {int(st.median(orfs)):>6} "
          f"{100*sum(o>=KEEP_O for o in orfs)/len(rs):>8.0f}% {st.median(cx):>7.2f} "
          f"{100*sum(c<KEEP_C for c in cx)/len(rs):>8.0f}% {int(st.median(r['n_exons'] for r in rs)):>8}")

# proposed coding-gene-family filter: ORF>=100 AND cplx>=0.40
print(f"\n=== proposed filter: ORF>=100aa AND 4-mer-complexity>=0.40 (annotation-free) ===")
print(f"{'biotype':16s} {'kept%':>6}  (protein_coding SHOULD stay high; lncRNA + low-cplx SHOULD drop)")
for bt in ("protein_coding", "lncRNA", "pseudogene", "unannotated"):
    rs = by_bt.get(bt, [])
    if not rs: continue
    kept = sum(1 for r in rs if r["orf"] >= KEEP_O and r["cplx"] >= KEEP_C)
    print(f"{bt:16s} {100*kept/len(rs):>5.0f}%")
json.dump({bt: rs for bt, rs in by_bt.items()}, open("bench/refine_family_definition.json", "w"))
print("\nwrote bench/refine_family_definition.json")
