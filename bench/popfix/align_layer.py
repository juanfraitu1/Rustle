#!/usr/bin/env python3
"""Alignment layer of spec 2026-10-12-population-fixed-psvs §2: for copies whose +- 10 kb flanks are unique in the
reference, each man's ortholog is the single alignment of locus + flanks spanning both flanks; his base at every PSV
position is read through it. A column gets an alignment verdict only if ALL its sites qualify.

    align_layer.py PREFIX CATALOG_COPIES_TSV GENOME_FA ASSEMBLY_DIR KSCAN OUTDIR
"""
import collections, math, os, subprocess, sys
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "psv_sun"))
import kmer_layer as KL
import psv_sun as P

FLANK, STEP = 10_000, 100


def ortholog_base(recs, qlen, flank, site_offset, fetch_target):
    """The man's base at query offset `site_offset`, through the single alignment spanning both flanks, else None.
    Returned in FORWARD reference orientation: complemented when the alignment is on the '-' strand."""
    span = [r for r in recs if r["tp"] == "P" and r["qs"] < flank and r["qe"] > qlen - flank]
    if len(span) != 1:
        return None
    r = span[0]
    t = P.paf_qmap(r).get(site_offset)
    if t is None:
        return None
    b = fetch_target(r["tname"], t, t + 1).upper()
    return b.translate(KL.COMP) if r["strand"] == "-" else b


def unique_flanks(F, kscan, genome_fa, sites, workdir):
    """site -> True if each +-10 kb flank has a 31-mer occurring once in the reference (sampled every STEP bp)."""
    probe = {}
    for (c, p) in sites:
        for side, (s, e) in (("L", (p - FLANK, p - 1000)), ("R", (p + 1000, p + FLANK))):
            seq = F.fetch(c, max(0, s), e).upper()
            probe[(c, p, side)] = [KL.canon(seq[i:i + KL.K]) for i in range(0, len(seq) - KL.K, STEP)
                                   if set(seq[i:i + KL.K]) <= set("ACGT")]
    kf = os.path.join(workdir, "flank_kmers.txt")
    with open(kf, "w") as fh:
        fh.write("\n".join(sorted({k for v in probe.values() for k in v})) + "\n")
    res = subprocess.run([kscan, str(KL.K), kf, genome_fa], capture_output=True, text=True, check=True).stdout
    once = {k for k, n in (l.split("\t") for l in res.splitlines()) if int(n) == 1}
    return {(c, p): any(k in once for k in probe[(c, p, "L")]) and any(k in once for k in probe[(c, p, "R")])
            for (c, p) in sites}


def main(prefix, catalog_tsv, genome_fa, asm_dir, kscan, out):
    os.makedirs(out, exist_ok=True)
    F = P.Fasta(genome_fa)
    men = KL.find_men(asm_dir)
    cols = KL.columns(prefix, catalog_tsv, lambda c, s, e: F.fetch(c, s, e))
    ycols = {k: v for k, v in cols.items() if all(c == "chrY" for _, c, *_ in v)}
    sites = sorted({(c, p) for v in ycols.values() for _, c, p, *_ in v})
    ok = unique_flanks(F, kscan, genome_fa, sites, out)
    qual = [s for s in sites if ok[s]]
    qfa = os.path.join(out, "loci.fa")
    with open(qfa, "w") as fh:
        for c, p in qual:
            fh.write(f">{c}:{p}\n{F.fetch(c, p - FLANK, p + FLANK + 1)}\n")
    base = collections.defaultdict(dict)   # site -> man -> base or None
    for mi, fa in enumerate(men):
        T = P.Fasta(fa)
        recs = collections.defaultdict(list)
        for r in P.run_paf(["minimap2", "-cx", "asm5", "--cs", "-t", "2", fa, qfa]):
            recs[r["qname"]].append(r)
        for c, p in qual:
            base[(c, p)][mi] = ortholog_base(recs.get(f"{c}:{p}", []), 2 * FLANK + 1, FLANK, FLANK,
                                             lambda n, s, e: T.fetch(n, s, e))
    half = math.ceil(len(men) / 2)
    with open(os.path.join(out, "align_status.tsv"), "w") as fh:
        fh.write("family_id\tcol_index\tstatus\tn_men\tn_ortholog\tn_mismatch\n")
        for (fam, col), v in sorted(ycols.items()):
            if not all(ok[(c, p)] for _, c, p, *_ in v):
                fh.write(f"{fam}\t{col}\tNO_ALIGN\t{len(men)}\t0\t0\n")
                continue
            have = [m for m in range(len(men)) if all(base[(c, p)].get(m) for _, c, p, *_ in v)]
            mism = sum(1 for m in have if any(base[(c, p)][m] != F.fetch(c, p, p + 1).upper() for _, c, p, *_ in v))
            st = "UNCOVERED" if len(have) < half else ("FIXED" if mism == 0 else "POLYMORPHIC")
            fh.write(f"{fam}\t{col}\t{st}\t{len(men)}\t{len(have)}\t{mism}\n")
    print(f"{len(qual)}/{len(sites)} sites with unique flanks, {len(men)} assemblies -> {out}/align_status.tsv")


if __name__ == "__main__":
    main(*sys.argv[1:7])
