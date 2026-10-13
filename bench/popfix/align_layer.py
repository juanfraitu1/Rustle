#!/usr/bin/env python3
"""Alignment layer of spec 2026-10-12-population-fixed-psvs §2: for copies whose +- 10 kb flanks are unique in the
reference, each man's ortholog is the single alignment of locus + flanks spanning both flanks; his base at every PSV
position is read through it. A column gets an alignment verdict only if ALL its sites qualify.

    align_layer.py PREFIX CATALOG_COPIES_TSV GENOME_FA ASSEMBLY_DIR KSCAN OUTDIR

Resumable: each man's minimap2 PAF is cached as OUTDIR/paf_<md5 of loci.fa>/<assembly file>.paf (written as .paf.part,
then renamed, so a .paf is always complete). POPFIX_THREADS (default 2) sets minimap2 -t; POPFIX_MAX_NEW_PAF=k stops
after k new men with exit status 3 (re-run to continue); align_status.tsv is written only once every man is cached.
"""
import collections, hashlib, math, os, subprocess, sys
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "psv_sun"))
import kmer_layer as KL
import psv_sun as P

FLANK, STEP = 10_000, 100


def ortholog_base(recs, qlen, flank, site_offset, fetch_target):
    """The man's base at query offset `site_offset`, through the single alignment spanning both flanks, else None.
    Returned in FORWARD reference orientation: complemented when the alignment is on the '-' strand."""
    # R8: every record spanning both flanks counts, primary or secondary (a tandem duplicate in the man yields a
    # secondary spanning alignment); the base is taken only when exactly one such record exists.
    span = [r for r in recs if r["qs"] < flank and r["qe"] > qlen - flank]
    if len(span) != 1:
        return None
    r = span[0]
    t = P.paf_qmap(r).get(site_offset)
    if t is None:
        return None
    b = fetch_target(r["tname"], t, t + 1).upper()
    return b.translate(KL.COMP) if r["strand"] == "-" else b


def locus_window(fetch, length, chrom, p):
    """(locus sequence, site offset) for a site at p with +- FLANK bp on each side, or None when the window would run
    off the contig (p < FLANK or p + FLANK + 1 > contig length). Such sites do not qualify, so no clamped window is
    ever used with a shifted site offset."""
    if p - FLANK < 0 or p + FLANK + 1 > length(chrom):
        return None
    return fetch(chrom, p - FLANK, p + FLANK + 1), FLANK


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


def paf_path(cache, fa):
    return os.path.join(cache, os.path.basename(fa) + ".paf")


def cached_records(cache, fa):
    """The man's parsed PAF records when his cache file is complete, else None. A `.paf.part` is never a cache hit."""
    path = paf_path(cache, fa)
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return P.parse_paf(fh.read())


def compute_paf(cache, fa, qfa, threads):
    """Align the loci to the man and keep minimap2's raw stdout as `<man>.paf`. Written to `.paf.part` first and renamed
    only after minimap2 succeeded and the bytes are on disk, so a crash never leaves a truncated `.paf`."""
    path = paf_path(cache, fa)
    with open(path + ".part", "wb") as fh:
        subprocess.run(["minimap2", "-cx", "asm5", "--cs", "-t", str(threads), fa, qfa],
                       stdout=fh, stderr=subprocess.PIPE, check=True)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(path + ".part", path)


def main(prefix, catalog_tsv, genome_fa, asm_dir, kscan, out):
    os.makedirs(out, exist_ok=True)
    status_tsv = os.path.join(out, "align_status.tsv")
    if os.path.exists(status_tsv):
        os.remove(status_tsv)                       # never leave a verdict from an earlier run looking current
    F = P.Fasta(genome_fa)
    men = KL.find_men(asm_dir)
    cols = KL.columns(prefix, catalog_tsv, lambda c, s, e: F.fetch(c, s, e))
    ycols = {k: v for k, v in cols.items() if all(c == "chrY" for _, c, *_ in v)}
    sites = sorted({(c, p) for v in ycols.values() for _, c, p, *_ in v})
    ok = unique_flanks(F, kscan, genome_fa, sites, out)
    win = {s: locus_window(F.fetch, F.length, *s) for s in sites}
    for s in sites:
        ok[s] = ok[s] and win[s] is not None          # edge sites: NO_ALIGN, never a clamped window
    qual = [s for s in sites if ok[s]]
    qfa = os.path.join(out, "loci.fa")
    with open(qfa, "w") as fh:
        for c, p in qual:
            fh.write(f">{c}:{p}\n{win[(c, p)][0]}\n")
    with open(qfa, "rb") as fh:
        tag = hashlib.md5(fh.read()).hexdigest()[:12]
    cache = os.path.join(out, f"paf_{tag}")
    os.makedirs(cache, exist_ok=True)
    threads = int(os.environ.get("POPFIX_THREADS", "2"))
    cap = os.environ.get("POPFIX_MAX_NEW_PAF")
    cap = None if cap is None else int(cap)
    n_new = 0
    for fa in men:
        if os.path.exists(paf_path(cache, fa)):
            continue
        if cap is not None and n_new >= cap:        # budget spent and this man is still missing: stop, resume later
            n_done = sum(os.path.exists(paf_path(cache, f)) for f in men)
            print(f"incomplete: {n_done}/{len(men)} men cached; re-run", file=sys.stderr)
            sys.exit(3)
        compute_paf(cache, fa, qfa, threads)
        n_new += 1
    base = collections.defaultdict(dict)   # site -> man -> base or None
    for mi, fa in enumerate(men):
        T = P.Fasta(fa)
        recs = collections.defaultdict(list)
        for r in cached_records(cache, fa):
            recs[r["qname"]].append(r)
        for c, p in qual:
            locus, off = win[(c, p)]
            base[(c, p)][mi] = ortholog_base(recs.get(f"{c}:{p}", []), len(locus), FLANK, off,
                                             lambda n, s, e: T.fetch(n, s, e))
    half = math.ceil(len(men) / 2)
    with open(status_tsv, "w") as fh:
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
