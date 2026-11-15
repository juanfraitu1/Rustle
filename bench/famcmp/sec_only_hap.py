#!/usr/bin/env python3
"""Row 1673: which haplotype supports the class-L loci of sec_only_loci.py (docs/PREREG_secondary_only_haplotype_2026-11-13.md).

    SEC_ONLY_SLACK=0 sec_only_hap.py GTF OUTDIR_OF_sec_only_loci MAT_X_CONTIG

Reads OUTDIR/{loci.tsv,mat.sam,pat.sam}; "pat" is the haplotype the reads were aligned to. No BAM pass:
a molecule counts for a locus when its near hits overlap the locus's hits, whatever its primary's locus.
"""
import collections, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sec_only_loci import read_loci, parse_sam, near_hits, copy_loci  # noqa: E402

MINJ = 3


def match_haps(near, lhits):
    """haplotypes on which a molecule's near hits overlap the locus's hits"""
    return {h for h, c, s, e, *_ in near for lh, lc, ls, le in lhits if lh == h and lc == c and ls < e and s < le}


def label(haps):
    """haps = one set per matching molecule"""
    if len(haps) < MINJ:
        return "thin"
    return "aligned_haplotype" if 2 * sum("pat" in h for h in haps) >= len(haps) else "other_allele_only"


def selftest():
    n = [("mat", "a", 0, 10, 0, 10), ("pat", "b", 0, 10, 0, 10)]
    assert match_haps(n, {("mat", "a", 5, 20), ("pat", "c", 0, 10)}) == {"mat"}
    assert label([{"pat"}, {"mat"}]) == "thin" and label([{"pat", "mat"}, {"pat"}, {"mat"}]) == "aligned_haplotype"
    assert label([{"mat"}, {"mat"}, {"pat"}]) == "other_allele_only"
    print("selftest ok")


def main(gtf, out, matx):
    loci = read_loci(gtf)
    hits = collections.defaultdict(list)
    for hap in ("mat", "pat"):
        for q, hs in parse_sam(open(f"{out}/{hap}.sam")).items():
            hits[q].extend((hap,) + h for h in hs)
    near = {q: near_hits(v) for q, v in hits.items() if q.startswith("R|")}
    near = {q: v for q, v in near.items() if v}
    kind = collections.Counter()  # every sampled molecule by where its strict best hits are
    for v in near.values():
        hp = {h[0] for h in v}
        kind["both haplotypes" if len(hp) == 2 else "paternal only" if hp == {"pat"} else
             "maternal X only" if all(h[1] == matx for h in v) else "maternal autosome only"] += 1
    print(f"sampled molecules with a near-perfect hit {len(near)}: " + ", ".join(f"{k} {n}" for k, n in kind.most_common()))
    by = collections.defaultdict(list)  # (hap, chrom) -> [(s, e, q)] of near hits
    for q, v in near.items():
        for h in v:
            by[h[:2]].append((h[2], h[3], q))
    idx = {c: {(s, e): i for i, (s, e, *_) in enumerate(v)} for c, v in loci.items()}
    tab, rows = collections.Counter(), []
    for l in open(f"{out}/loci.tsv"):
        f = l.rstrip("\n").split("\t")
        if f[-1] != "L":
            continue
        c, i = f[0], idx[f[0]][(int(f[1]), int(f[2]))]
        lh = copy_loci(hits.get(f"C|{c}|{i}", []))
        cand = {q for h, ch, s, e in lh for a, b, q in by.get((h, ch), []) if a < e and s < b}
        haps = [m for m in (match_haps(near[q], lh) for q in sorted(cand)) if m]
        onx = sum(any(h[0] == "mat" and h[1] == matx for h in near[q]) for q in cand)
        k = label(haps)
        tab[k] += 1
        rows.append((c, f[1], f[2], f[6], len(haps), sum(h == {"pat"} for h in haps), sum(h == {"mat"} for h in haps),
                     sum(len(h) == 2 for h in haps), onx, k))
    with open(f"{out}/L_haplotype.tsv", "w") as fh:
        fh.write("chrom\tstart\tend\tannot\tmatching\tpat_only\tmat_only\tboth\ton_mat_X\tlabel\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    print(f"class-L loci {len(rows)}: " + ", ".join(f"{k} {n}" for k, n in tab.most_common()))
    print(f"class-L loci with any matching molecule on maternal X: {sum(r[8] > 0 for r in rows)}")
    print(f"thin {tab['thin']} -> {'RERUN (thin > 35)' if tab['thin'] > 35 else 'step 0 stands'}")


if __name__ == "__main__":
    selftest() if sys.argv[1:] == ["--selftest"] else main(*sys.argv[1:])
