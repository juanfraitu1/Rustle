#!/usr/bin/env python3
"""Row 1671: loci of the default assembly supported only by secondary alignments. Real or shadow?

  sec_only_loci.py GTF BAM GENOME.fa ANNOT.gff MAT.mmi PAT.mmi OUTDIR [THREADS]
  sec_only_loci.py --selftest

Locus = merged same-strand overlapping transcripts of GTF. Population = loci where no primary record has an aligned block
on a locus exon (either strand; supplementary records ignored) and >= 1 molecule has a secondary record there inside the
default seeding pool (AS >= 0.98 x its primary's AS).
Truth (the animal's own haplotypes, as bench/o2_truth/haplotype_truth.py): up to CAP molecules per locus, plus the longest
transcript of the locus L and of every assembly locus P that holds those molecules' primaries, are aligned to both
haplotypes. A molecule with a near-perfect best hit is `L` (its near hits fall on L's haplotype loci and not on its own
P's), `P` (the reverse), `both`, or `other`. A locus with >= MINJ judged molecules takes the class holding >= half of
them (else `mixed`); fewer = `thin`.
Writes OUTDIR/loci.tsv (one row per population locus: annotation class, features, truth class) and prints the tables.
"""
import bisect, collections, hashlib, os, re, statistics, subprocess, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "o2_truth"))
from haplotype_truth import parse_sam, near_hits, copy_loci, truth_copies, COMP  # noqa: E402

CAP, MINJ = 30, 3
SLACK = int(os.environ.get("SEC_ONLY_SLACK", "2"))  # NM slack of the haplotype truth; 0 = strict re-read (post hoc)
if SLACK != 2:
    _nh = near_hits
    def near_hits(hits):  # noqa: F811
        n = _nh(hits)
        return n and [h for h in n if h[4] <= min(x[4] for x in n) + SLACK]
    def copy_loci(hits):  # noqa: F811
        best = {}
        for h in hits:
            best[h[0]] = min(best.get(h[0], 10 ** 9), h[4])
        return {(h[0], h[1], h[2], h[3]) for h in hits if h[4] <= best[h[0]] + SLACK}


def read_loci(gtf):
    """-> {chrom: [(start, end, strand, [exon lists per transcript])]}"""
    ex = collections.defaultdict(list)
    for l in open(gtf):
        f = l.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "exon":
            ex[(f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]) - 1, int(f[4])))
    by = collections.defaultdict(list)
    for (c, s, _), e in ex.items():
        by[(c, s)].append(sorted(e))
    out = collections.defaultdict(list)
    for (c, s), ts in by.items():
        ts.sort(key=lambda e: e[0][0])
        cur, end = [ts[0]], ts[0][-1][1]
        for e in ts[1:]:
            if e[0][0] < end:
                cur.append(e); end = max(end, e[-1][1])
            else:
                out[c].append((cur[0][0][0], end, s, cur)); cur, end = [e], e[-1][1]
        out[c].append((cur[0][0][0], end, s, cur))
    for v in out.values():
        v.sort(key=lambda x: x[:2])
    return out


class Index:
    """Exons of one chromosome's loci -> locus ids hit by a list of blocks."""
    def __init__(self, loci):
        self.ex = sorted((a, b, i) for i, (_, _, _, ts) in enumerate(loci) for a, b in {x for t in ts for x in t})
        self.starts = [x[0] for x in self.ex]
        self.mx = max((b - a for a, b, _ in self.ex), default=0)

    def hit(self, blocks):
        out = set()
        for a, b in blocks:
            i = bisect.bisect_left(self.starts, a - self.mx)
            while i < len(self.ex) and self.ex[i][0] < b:
                if self.ex[i][1] > a:
                    out.add(self.ex[i][2])
                i += 1
        return out


def tag(r, t):
    return r.get_tag(t) if r.has_tag(t) else None


def pass1(a):
    """One chromosome: per locus primary count and (if none) its secondary records."""
    import pysam
    bam, chrom, loci = a
    ix, prim, sec = Index(loci), collections.Counter(), collections.defaultdict(list)
    for r in pysam.AlignmentFile(bam).fetch(chrom):
        if r.is_unmapped or r.is_supplementary:
            continue
        for i in ix.hit(r.get_blocks()):
            if r.is_secondary:
                if not prim[i]:
                    sec[i].append((r.query_name, tag(r, "AS"), tag(r, "NM"), "N" in r.cigarstring))
            else:
                prim[i] += 1; sec.pop(i, None)
    return chrom, {i: v for i, v in sec.items() if not prim[i]}


def pass2(a):
    """One chromosome: the wanted molecules' primaries -> name: (seq, AS, NM, [locus ids of this chromosome])."""
    import pysam
    bam, chrom, loci, want = a
    ix, out = Index(loci), {}
    for r in pysam.AlignmentFile(bam).fetch(chrom):
        if r.is_unmapped or r.is_supplementary or r.is_secondary or r.query_name not in want or not r.query_sequence:
            continue
        s = r.query_sequence
        out[r.query_name] = (s.translate(COMP)[::-1] if r.is_reverse else s, tag(r, "AS"), tag(r, "NM"), sorted(ix.hit(r.get_blocks())))
    return chrom, out


def annot_class(gff):
    """-> {chrom: Index-like (starts, [(s, e, class)])} over gene rows; class from gene_biotype."""
    d = collections.defaultdict(list)
    for l in open(gff):
        f = l.split("\t")
        if len(f) > 8 and f[2] in ("gene", "pseudogene", "ncRNA_gene"):  # RefSeq and Ensembl gene rows
            b = re.search(r"(?:gene_)?biotype=([^;\n]+)", f[8])
            b = b.group(1) if b else f[2]
            d[f[0]].append((int(f[3]) - 1, int(f[4]), "coding" if b == "protein_coding" else "pseudogene" if "pseudo" in b else "other_gene"))
    return {c: sorted(v) for c, v in d.items()}


def gene_class(genes, s, e):
    got = {k for a, b, k in genes if a < e and s < b}  # shortcut: linear scan per locus, fine for a few thousand loci
    return next((k for k in ("coding", "pseudogene", "other_gene") if k in got), "none")


def mol_class(truth, L, Ps):
    inL, inP = L in truth, bool(truth & Ps)
    return "both" if inL and inP else "L" if inL else "P" if inP else "other"


def locus_class(classes):
    if len(classes) < MINJ:
        return "thin"
    k, n = collections.Counter(classes).most_common(1)[0]
    return k if 2 * n >= len(classes) else "mixed"


def selftest():
    ix = Index([(0, 100, "+", [[(0, 10), (50, 100)]]), (200, 300, "+", [[(200, 300)]])])
    assert ix.hit([(10, 50)]) == set() and ix.hit([(5, 12), (250, 260)]) == {0, 1}
    assert mol_class({"a"}, "a", {"b"}) == "L" and mol_class({"a", "b"}, "a", {"b"}) == "both"
    assert mol_class({"b"}, "a", {"b"}) == "P" and mol_class(set(), "a", {"b"}) == "other"
    assert locus_class(["L", "L"]) == "thin" and locus_class(["L", "L", "P"]) == "L" and locus_class(["L", "P", "both", "other"]) == "mixed"
    print("selftest ok")


def main(gtf, bam, fasta, gff, mat, pat, out, threads=24):
    import multiprocessing, pysam
    os.makedirs(out, exist_ok=True)
    loci = read_loci(gtf)
    chroms = [c for c in pysam.AlignmentFile(bam).references if c in loci]
    with multiprocessing.Pool(int(threads)) as pool:
        sec = dict(pool.map(pass1, [(bam, c, loci[c]) for c in chroms]))
        pop = [(c, i) for c in chroms for i in sorted(sec[c])]
        print(f"loci {sum(map(len, loci.values()))}  supported only by secondary alignments {len(pop)}", flush=True)
        want = {r[0] for c, i in pop for r in sec[c][i]}
        prim = {}
        for c, d in pool.map(pass2, [(bam, c, loci[c], want) for c in chroms]):
            prim.update({n: (s, a, nm, [(c, j) for j in ids]) for n, (s, a, nm, ids) in d.items()})
    pick = {}  # (chrom, i) -> up to CAP molecules whose best record here is in the default pool (AS >= 0.98 x the primary's)
    for c, i in pop:
        ok = {n for n, a, _, _ in sec[c][i] if n in prim and a is not None and prim[n][1] and a >= 0.98 * prim[n][1]}
        pick[(c, i)] = sorted(ok, key=lambda n: hashlib.md5(f"7:{n}".encode()).hexdigest())[:CAP]
    print(f"of those with no pool molecule (set aside) {sum(not pick[k] for k in pop)}", flush=True)
    pop = [k for k in pop if pick[k]]
    want = {n for v in pick.values() for n in v}
    fa = pysam.FastaFile(fasta)
    need = set(pop) | {p for n in want if n in prim for p in prim[n][3]}
    with open(f"{out}/query.fa", "w") as fh:
        for n in sorted(want & set(prim)):
            fh.write(f">R|{n}\n{prim[n][0]}\n")
        for c, i in sorted(need):
            s, e, st, ts = loci[c][i]
            t = max(ts, key=lambda t: sum(b - a for a, b in t))
            seq = "".join(fa.fetch(c, a, b) for a, b in t).upper()
            fh.write(f">C|{c}|{i}\n{seq.translate(COMP)[::-1] if st == '-' else seq}\n")
    hits = collections.defaultdict(list)
    for hap, idx in (("mat", mat), ("pat", pat)):
        with open(f"{out}/{hap}.sam", "w") as fh:
            subprocess.run(["minimap2", "-ax", "splice:hq", "-uf", "--eqx", "-N", "10", "-t", str(threads), idx, f"{out}/query.fa"],
                           stdout=fh, stderr=subprocess.DEVNULL, check=True)
        for q, hs in parse_sam(open(f"{out}/{hap}.sam")).items():
            hits[q].extend((hap,) + h for h in hs)
    hl = {tuple(k.split("|")[1:]): copy_loci(v) for k, v in hits.items() if k.startswith("C|")}
    hl = {(c, int(i)): v for (c, i), v in hl.items()}
    genes = annot_class(gff)
    tab, rows = collections.Counter(), []
    for c, i in pop:
        s, e, st, ts = loci[c][i]
        recs = sec[c][i]
        mols = {n for n, a, _, _ in recs if n in prim and a is not None and prim[n][1] and a >= 0.98 * prim[n][1]}
        cl, tie, nmle, ratio, pl = [], [], [], [], collections.Counter()
        best = {}
        for n, a, nm, _ in recs:  # a molecule's best secondary record at this locus
            if n not in best or (a or 0) > (best[n][0] or 0):
                best[n] = (a, nm)
        for n in pick[(c, i)]:
            if n not in prim:
                continue
            _, pa, pnm, ps = prim[n]
            a, nm = best[n]
            if a is not None and pa:
                tie.append(a >= pa); ratio.append(a / pa)
            if nm is not None and pnm is not None:
                nmle.append(nm <= pnm)
            pl[ps[0] if ps else None] += 1
            near = near_hits(hits.get("R|" + n, []))
            if near is not None:
                cl.append(mol_class(truth_copies(near, {k: hl.get(k, set()) for k in [(c, i)] + ps}), (c, i), set(ps)))
        k, g = locus_class(cl), gene_class(genes.get(c, []), s, e)
        cnt = collections.Counter(cl)
        sh = lambda v: f"{sum(v) / len(v):.2f}" if v else "NA"
        top = pl.most_common(1)[0] if pl else (None, 0)
        rows.append((c, s, e, st, len(ts), max(map(len, ts)), g, len(mols), len(recs), sh(tie), sh(nmle),
                     f"{statistics.median(ratio):.4f}" if ratio else "NA", f"{top[1] / max(1, sum(pl.values())):.2f}",
                     int(top[0] is not None), len(cl), cnt["L"], cnt["both"], cnt["P"], cnt["other"], k))
        tab[(g, k)] += 1
    hdr = "chrom start end strand transcripts max_exons annot molecules records exact_tie nm_le_primary as_ratio top_primary_share top_primary_on_locus judged L both P other truth".split()
    with open(f"{out}/loci.tsv", "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    ks = ["L", "both", "P", "other", "mixed", "thin"]
    print("annot\t" + "\t".join(ks) + "\ttotal")
    for g in ("coding", "pseudogene", "other_gene", "none"):
        print(g + "\t" + "\t".join(str(tab[(g, k)]) for k in ks) + f"\t{sum(tab[(g, k)] for k in ks)}")
    print("total\t" + "\t".join(str(sum(tab[(g, k)] for g in ("coding", "pseudogene", "other_gene", "none"))) for k in ks) + f"\t{len(rows)}")
    for k in ks:  # feature medians per truth class
        sub = [r for r in rows if r[-1] == k]
        if not sub:
            continue
        med = lambda j: statistics.median([float(r[j]) for r in sub if r[j] != "NA"] or [float("nan")])
        print(f"{k}\tn {len(sub)}\tmolecules {med(7)}\texact_tie {med(9):.2f}\tnm_le_primary {med(10):.2f}\tas_ratio {med(11):.4f}"
              f"\ttop_primary_share {med(12):.2f}\tspliced {sum(r[5] > 1 for r in sub)}")
    real = sum(tab[(g, k)] for g in ("coding", "pseudogene", "other_gene", "none") for k in ("L", "both"))
    print(f"loci with class L or both: {real} -> {'GO' if real >= 10 else 'STOP'}")


if __name__ == "__main__":
    selftest() if sys.argv[1:] == ["--selftest"] else main(*sys.argv[1:])
