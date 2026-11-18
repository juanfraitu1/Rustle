"""PREREG r1794 (docs/PREREG_famwise_origin_2026-11-18.md), row b: intron chains and RefSeq-exact chains of the three
inputs, in the regions a tied or moved read touches and genome-wide.

  region_chains.py REF.gtf REGIONS.bed P=GTF D=GTF O=GTF

Prints `scope input chains refseq_exact`, then `R scope value` (RefSeq-exact chains, O over P) and
`D_not_O scope chains refseq_exact` (what perfect provenance removes from the default). A chain is in the regions when
the span from its first to its last intron overlaps a BED interval on its contig."""
import bisect, collections, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "famcmp"))
from chain_prec import chains

def main():
    ref, bed, arms = chains(sys.argv[1]), collections.defaultdict(list), dict(a.split("=", 1) for a in sys.argv[3:])
    assert set(arms) == {"P", "D", "O"}, "inputs: P=GTF D=GTF O=GTF"
    with open(sys.argv[2]) as fh:
        for l in fh:
            f = l.split()
            if len(f) > 2:
                bed[f[0]].append((int(f[1]), int(f[2])))
    for v in bed.values():
        v.sort()
    def inside(k):   # BED intervals are merged (origin_bam.py), so the one starting last before the chain's end decides
        v = bed.get(k[0], ()); i = bisect.bisect_left(v, (k[2][-1][1],)) - 1
        return i >= 0 and v[i][1] > k[2][0][0]
    cs = {n: chains(p) for n, p in arms.items()}
    scope = {"regions": lambda s: {k for k in s if inside(k)}, "genome": lambda s: s}
    print("scope\tinput\tchains\trefseq_exact")
    for sc, f in scope.items():
        for n in "PDO":
            print(f"{sc}\t{n}\t{len(f(cs[n]))}\t{len(f(cs[n]) & ref)}")
    for sc, f in scope.items():
        p, o = len(f(cs["P"]) & ref), len(f(cs["O"]) & ref)
        print(f"R\t{sc}\t" + (f"{o / p:.2f}" if p else f"NA (P has none, O {o})"))
    for sc, f in scope.items():
        x = f(cs["D"]) - cs["O"]
        print(f"D_not_O\t{sc}\t{len(x)}\t{len(x & ref)}")

if __name__ == "__main__":
    main()
