#!/usr/bin/env python3
"""k-mer layer of spec 2026-10-12-population-fixed-psvs §2: a PSV column is FIXED when every man whose assembly covers
it still carries every allele of it.

    kmer_layer.py PREFIX CATALOG_COPIES_TSV GENOME_FA ASSEMBLY_DIR KSCAN OUTDIR
PREFIX = `utilities psv-discover` output prefix. Columns with a site off chrY are NOT_ASSESSED.

Allele identity (controller ruling R6): sites are grouped into alleles by the psv_pos `base` (copy-0 frame), never by the
forward genome base. Allele and flank k-mers come from the forward genome (canonical, strand-agnostic).
"""
import collections, csv, glob, math, os, subprocess, sys

K, FLANK_GAP = 31, 40
COMP = str.maketrans("ACGTN", "TGCAN")


def canon(k):
    r = k.translate(COMP)[::-1]
    return k if k <= r else r


def revcomp_base(b):
    return b.translate(COMP)


def fetch(genome, chrom, s, e):
    """genome: dict chrom -> str (tests) or a callable(chrom, s, e) -> str (CLI)."""
    if callable(genome):
        return genome(chrom, s, e)
    return genome[chrom][max(0, s):e]


def allele_kmers(genome, chrom, pos):
    w = fetch(genome, chrom, pos - (K - 1), pos + K).upper()
    return sorted({canon(w[i:i + K]) for i in range(len(w) - K + 1) if set(w[i:i + K]) <= set("ACGT")})


def flank_kmers(genome, chrom, pos):
    left = fetch(genome, chrom, pos - FLANK_GAP - (K - 1), pos - FLANK_GAP + 1).upper()
    right = fetch(genome, chrom, pos + FLANK_GAP, pos + FLANK_GAP + K).upper()
    return canon(left), canon(right)


def decide(alleles, cover, carry, n_men):
    """alleles: the column's allele labels; cover[m]: man m covers the column; carry[m]: alleles man m carries.
    UNCOVERED unless at least one man and at least ceil(n_men/2) men cover the column."""
    covering = [m for m in range(len(cover)) if cover[m]]
    if not covering or len(covering) < math.ceil(n_men / 2):
        return "UNCOVERED"
    alleles = set(alleles)
    return "FIXED" if all(alleles <= carry[m] for m in covering) else "POLYMORPHIC"


def frame_check(fam, col, copy, fwd, base, strand):
    """psv base is in copy-0 frame: equals the forward base on '+', its complement on '-'."""
    want = fwd if strand == "+" else revcomp_base(fwd)
    if base != want:
        raise ValueError(f"frame check failed: family {fam} col {col} copy {copy}: "
                         f"forward base {fwd}, psv base {base}, strand {strand}")
    return base


def tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def columns(prefix, catalog_tsv, genome):
    """-> {(fam, col): [Site]}, Site = (copy_index, chrom, pos, base, frame_strand). Raises ValueError on a frame mismatch."""
    chrom_of = {r["tid"]: r["chrom"] for r in tsv(catalog_tsv)}
    copy_of = {(r["family_id"], int(r["copy_index"])): (r["copy_tid"], r["strand"])
               for r in tsv(prefix + ".psv_copies.tsv")}
    cols = collections.defaultdict(list)
    for r in tsv(prefix + ".psv_pos.tsv"):
        fam, col, copy = r["family_id"], int(r["col_index"]), int(r["copy_index"])
        tid, strand = copy_of[(fam, copy)]
        c, pos, base = chrom_of[tid], int(r["genome_pos"]), r["base"].upper()
        fwd = fetch(genome, c, pos, pos + 1).upper()
        frame_check(fam, col, copy, fwd, base, strand)
        cols[(fam, col)].append((copy, c, pos, base, strand))
    return cols


def build_plan(cols, genome):
    """-> {(fam, col): None (NOT_ASSESSED) or (allele_kmers_by_site, flanks_by_site, allele_by_site)}."""
    plan = {}
    for key, sites in cols.items():
        if any(c != "chrY" for _, c, _, _, _ in sites):
            plan[key] = None
            continue
        al = {i: allele_kmers(genome, c, p) for i, c, p, _, _ in sites}
        fl = {i: flank_kmers(genome, c, p) for i, c, p, _, _ in sites}
        lab = {i: b for i, _, _, b, _ in sites}
        plan[key] = (al, fl, lab)
    return plan


def kmer_status(plan, presence_iter, n_men):
    """Pure. presence_iter yields one set of present canonical k-mers per man (exactly n_men of them).
    -> {(fam, col): {"status", "n_cover", "n_carry_all"}}."""
    live = [k for k, p in plan.items() if p is not None]
    cover = {k: [] for k in live}
    carry = {k: [] for k in live}
    seen = 0
    for pm in presence_iter:
        seen += 1
        for key in live:
            al, fl, lab = plan[key]
            cover[key].append(all(l in pm and r in pm for l, r in fl.values()))
            carry[key].append({lab[i] for i, ks in al.items() if any(k in pm for k in ks)})
    if seen != n_men:
        raise ValueError(f"got {seen} presence sets for {n_men} assemblies")
    out = {}
    for key, p in plan.items():
        if p is None:
            out[key] = {"status": "NOT_ASSESSED", "n_cover": 0, "n_carry_all": 0}
            continue
        alleles = set(p[2].values())
        st = decide(alleles, cover[key], carry[key], n_men)
        n_all = sum(1 for m in range(n_men) if cover[key][m] and alleles <= carry[key][m])
        out[key] = {"status": st, "n_cover": sum(cover[key]), "n_carry_all": n_all}
    return out


def find_men(asm_dir):
    men = sorted(glob.glob(os.path.join(asm_dir, "*.chrY.fasta")))
    if not men:
        sys.exit(f"no *.chrY.fasta under {asm_dir}; refusing to call any column FIXED")
    return men


def presence_sets(men, kfile, kscan):
    """Yields one presence set per man; each set is freed before the next man is scanned."""
    for fa in men:
        res = subprocess.run([kscan, str(K), kfile, fa], capture_output=True, text=True, check=True).stdout
        yield {k for k, n in (l.split("\t") for l in res.splitlines()) if int(n) > 0}


def main(prefix, catalog_tsv, genome_fa, asm_dir, kscan, out):
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "psv_sun"))
    import psv_sun
    F = psv_sun.Fasta(genome_fa)
    genome = lambda c, s, e: F.fetch(c, s, e)
    men = find_men(asm_dir)
    os.makedirs(out, exist_ok=True)
    plan = build_plan(columns(prefix, catalog_tsv, genome), genome)
    want = set()
    for p in plan.values():
        if p is None:
            continue
        al, fl, _ = p
        for ks in al.values():
            want.update(ks)
        for l, r in fl.values():
            want.update((l, r))
    kfile = os.path.join(out, "kmers.txt")
    with open(kfile, "w") as fh:
        fh.write("\n".join(sorted(want)) + "\n")
    res = kmer_status(plan, presence_sets(men, kfile, kscan), len(men))
    with open(os.path.join(out, "kmer_status.tsv"), "w") as fh:
        fh.write("family_id\tcol_index\tstatus\tn_men\tn_cover\tn_carry_all\n")
        for (fam, col), r in sorted(res.items()):
            fh.write(f"{fam}\t{col}\t{r['status']}\t{len(men)}\t{r['n_cover']}\t{r['n_carry_all']}\n")
    print(f"{len(plan)} columns, {len(men)} assemblies -> {out}/kmer_status.tsv")


if __name__ == "__main__":
    main(*sys.argv[1:7])
