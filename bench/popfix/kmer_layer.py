#!/usr/bin/env python3
"""k-mer layer of spec 2026-10-12-population-fixed-psvs §2: a PSV column is FIXED when every man whose assembly covers
it still carries every allele of it.

    kmer_layer.py PREFIX CATALOG_COPIES_TSV GENOME_FA ASSEMBLY_DIR KSCAN OUTDIR
PREFIX = `utilities psv-discover` output prefix. Columns with a site off chrY are NOT_ASSESSED.
"""
import collections, csv, glob, math, os, subprocess, sys

K, FLANK_GAP = 31, 40
COMP = str.maketrans("ACGTN", "TGCAN")


def canon(k):
    r = k.translate(COMP)[::-1]
    return k if k <= r else r


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


def decide(allele_sites, cover, carry, n_men):
    """allele_sites: allele -> sites; cover[m]: man m covers the column; carry[m]: alleles man m carries."""
    covering = [m for m in range(len(cover)) if cover[m]]
    if len(covering) < math.ceil(n_men / 2):
        return "UNCOVERED"
    alleles = set(allele_sites)
    return "FIXED" if all(alleles <= carry[m] for m in covering) else "POLYMORPHIC"


def tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def columns(prefix, catalog_tsv):
    chrom_of = {r["tid"]: r["chrom"] for r in tsv(catalog_tsv)}
    tid = {(r["family_id"], r["copy_index"]): r["copy_tid"] for r in tsv(prefix + ".psv_copies.tsv")}
    cols = collections.defaultdict(list)
    for r in tsv(prefix + ".psv_pos.tsv"):
        c = chrom_of[tid[(r["family_id"], r["copy_index"])]]
        cols[(r["family_id"], int(r["col_index"]))].append((int(r["copy_index"]), c, int(r["genome_pos"])))
    return cols


def main(prefix, catalog_tsv, genome_fa, asm_dir, kscan, out):
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "psv_sun"))
    import psv_sun
    F = psv_sun.Fasta(genome_fa)
    genome = lambda c, s, e: F.fetch(c, s, e)
    os.makedirs(out, exist_ok=True)
    cols = columns(prefix, catalog_tsv)
    plan, want = {}, set()
    for key, sites in cols.items():
        if any(c != "chrY" for _, c, _ in sites):
            plan[key] = None
            continue
        al = {i: allele_kmers(genome, c, p) for i, c, p in sites}
        fl = {i: flank_kmers(genome, c, p) for i, c, p in sites}
        base = {i: F.fetch(c, p, p + 1).upper() for i, c, p in sites}
        plan[key] = (al, fl, base)
        for ks in al.values():
            want.update(ks)
        for l, r in fl.values():
            want.update((l, r))
    kfile = os.path.join(out, "kmers.txt")
    with open(kfile, "w") as fh:
        fh.write("\n".join(sorted(want)) + "\n")
    men = sorted(glob.glob(os.path.join(asm_dir, "*.chrY.fasta")))
    present = []
    for fa in men:
        res = subprocess.run([kscan, str(K), kfile, fa], capture_output=True, text=True, check=True).stdout
        present.append({k for k, n in (l.split("\t") for l in res.splitlines()) if int(n) > 0})
    with open(os.path.join(out, "kmer_status.tsv"), "w") as fh:
        fh.write("family_id\tcol_index\tstatus\tn_men\tn_cover\tn_carry_all\n")
        for (fam, col), p in sorted(plan.items()):
            if p is None:
                fh.write(f"{fam}\t{col}\tNOT_ASSESSED\t{len(men)}\t0\t0\n")
                continue
            al, fl, base = p
            allele_sites = collections.defaultdict(list)
            for i, b in base.items():
                allele_sites[b].append(i)
            cover = [all(l in pm and r in pm for l, r in fl.values()) for pm in present]
            carry = [{b for b, ids in allele_sites.items() if any(any(k in pm for k in al[i]) for i in ids)}
                     for pm in present]
            st = decide(allele_sites, cover, carry, len(men))
            n_cover = sum(cover)
            n_all = sum(1 for m in range(len(men)) if cover[m] and set(allele_sites) <= carry[m])
            fh.write(f"{fam}\t{col}\t{st}\t{len(men)}\t{n_cover}\t{n_all}\n")
    print(f"{len(plan)} columns, {len(men)} assemblies -> {out}/kmer_status.tsv")


if __name__ == "__main__":
    main(*sys.argv[1:7])
