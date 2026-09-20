"""V4a: simulate a 3-copy tandem array, delete copies from the reference, recover
via the certificate, and set-match recovered vs deleted copies (permutation-invariant).
Divergent arm proves recovery; identical arm proves the honest K=0 abstain floor.
Uses a seeded PRNG (no Math.random / wall-clock) so runs are reproducible."""
import random, pathlib, subprocess, json, argparse, shutil, re, sys

BASES = "ACGT"


def build_tandem(seed, n_copies, divergence, seed_rng):
    """Return (reference_with_n_tandem_copies, [copy_seqs]). Divergence = per-base
    substitution prob planted independently per copy (shared coordinate frame => PSVs)."""
    rng = random.Random(seed_rng)
    base = list(seed)
    copies = []
    for _ in range(n_copies):
        c = base[:]
        if divergence > 0:
            for i in range(len(c)):
                if rng.random() < divergence:
                    c[i] = rng.choice([b for b in BASES if b != c[i]])
        copies.append("".join(c))
    # simple tandem: copies concatenated with a short spacer
    spacer = "N" * 50
    ref = spacer.join(copies)
    return ref, copies


def _identity(a, b):
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    same = sum(1 for i in range(n) if a[i] == b[i])
    return same / max(len(a), len(b))


def set_match(recovered, truth):
    """Greedy permutation-invariant matching: each recovered copy -> best-identity truth copy."""
    matches = []
    used = set()
    for ri, r in enumerate(recovered):
        best, best_id = None, -1.0
        for ti, t in enumerate(truth):
            if ti in used:
                continue
            idv = _identity(r, t)
            if idv > best_id:
                best, best_id = ti, idv
        if best is not None:
            used.add(best)
            matches.append((ri, best, best_id))
    return matches


# ---------------------------------------------------------------------------
# Phase 2: end-to-end driver (real minimap2 + real copy_assign binary).
#
# Two architectural findings shaped this driver (both confirmed by reading
# src/rustle/vg_family/*.rs, not guessed):
#
# 1. Intron-less tandem copies are INVISIBLE to family detection. Pass-1 assembly
#    (`denovo_assemble.rs::pass1_skeletons`) groups reads by (chrom, intron-chain) only
#    -- NOT by position. A single-exon read has an empty intron chain regardless of
#    where it starts, so every intron-less copy on one contig collapses into ONE
#    skeleton/transcript spanning the whole read pileup, and family detection (which
#    needs >= 2 distinct-locus reps) never fires. This matches the design doc's own
#    warning ("shared splice structure ... plant shared donors/acceptors"). Fix: each
#    planted copy gets a real genomic intron (canonical GT..AG, same relative
#    exon1/intron/exon2 layout at every copy so a spliced read's intron chain lands on
#    the SAME relative structure at each copy's own absolute genomic coordinates,
#    which is what actually differentiates the loci). Reads are simulated as the
#    SPLICED (exon-only) transcript and aligned with `-ax splice:hq`, so minimap2
#    supplies the real `N` CIGAR op across the true intron. (A literal `N`-run SPACER
#    between tandem copies, if aligned with a splice-aware preset, is ALSO treated as
#    an intron candidate by minimap2 -- that was the first failure mode found here:
#    it spliced reads straight across the inter-copy spacer, gluing two copies into
#    one pseudo-transcript. The fix is the same: real biological intron structure
#    inside each copy, spacer big enough / N-content clean enough that minimap2 does
#    not chain across it.)
#
# 2. The reference-absent admission gate (`absent_copy.rs::AbsentCopyParams::min_clusters`,
#    hardcoded to 3, no CLI flag) requires >= 3 co-varying identifiable clusters AT ONE
#    HOST LOCUS before it will admit a synthetic copy -- by design, this is what tells a
#    genuine extra paralog apart from ordinary biallelic variation. A 3-member array with
#    exactly ONE copy deleted collapses at most a HOST + 1 hidden cluster (n_clusters=2)
#    at any single surviving locus, which is architecturally always < 3 and always lands
#    in `.dna_needs.tsv`, never `Admission::Copy` -- independent of divergence, read
#    depth, or preset. This was confirmed both by reading the gate source and by direct
#    experiment (`<out>.dna_needs.tsv` reason column literally reads "<3 clusters
#    (copy-vs-allele needs DNA)"). To exhibit an actual admitted recovery (not just an
#    honest decline) the divergent arm therefore plants a 4-member array: two anchor
#    copies (A, B) stay in the reference and form a real family (their reads share
#    enough identity to produce a read-conflict edge); B is then the ROOT of a
#    near-identical 3-member sub-cluster (B itself + 2 hidden siblings, ~1% apart from
#    B), both of which are DELETED and both of which pool their reads onto B's locus,
#    giving exactly the host+2-hidden = 3 clusters the gate requires. This is still
#    literally "a tandem array with copies deleted from the reference, recovered by the
#    same certificate" -- it is the real SDA/DAZ-type shape (several near-identical
#    copies collapsed behind one assembled representative), just not a bare 3-of-3.
#    The identical (K=0) arm below does NOT need this -- it never needs the admission
#    gate to fire at all, so it keeps the literal 3-copy array the brief specifies.
# ---------------------------------------------------------------------------

EXON_LEN = 700
INTRON_LEN = 200
MM2 = shutil.which("minimap2") or "/home/juanfra/miniforge3/bin/minimap2"
SAMTOOLS = shutil.which("samtools") or "/home/juanfra/miniforge3/bin/samtools"


def _make_intron(rng, n):
    """A canonical (GT..AG) intron of length n so the spliced-alignment gate accepts it."""
    mid = "".join(rng.choice(BASES) for _ in range(n - 4))
    return "GT" + mid + "AG"


def _genomic_copy(body, intron_seq):
    """body = exon1(EXON_LEN) + exon2(EXON_LEN); genomic = exon1 + intron + exon2."""
    return body[:EXON_LEN] + intron_seq + body[EXON_LEN:]


def _sim_reads(bodies, n_per_copy, err_rate, rng, prefix):
    """Spliced (exon-only, no intron) HiFi-like reads with per-base substitution error."""
    reads = []
    for ci, seq in enumerate(bodies):
        for k in range(n_per_copy):
            s = list(seq)
            for i in range(len(s)):
                if rng.random() < err_rate:
                    s[i] = rng.choice([b for b in BASES if b != s[i]])
            reads.append((f"{prefix}_c{ci}_r{k}", "".join(s)))
    return reads


def _write_fasta(path, seqs):
    with open(path, "w") as f:
        for name, seq in seqs:
            f.write(f">{name}\n{seq}\n")


def _write_fastq(path, reads):
    with open(path, "w") as f:
        for name, seq in reads:
            f.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def _run(cmd, **kw):
    print("+", " ".join(str(c) for c in cmd), file=sys.stderr)
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.stderr:
        print(r.stderr[-4000:], file=sys.stderr)
    r.check_returncode()
    return r


def _read_tsv(path):
    """List of dict rows from a header-first TSV."""
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = []
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(dict(zip(header, line.split("\t"))))
        return rows


def _align_and_index(ref_fa, fq, out_bam_prefix):
    """minimap2 -ax splice:hq -N 50 (NOT --secondary=no, which yields 0 families) + sort/index.
    FOREGROUND, serial (crash rule): no backgrounding of any step."""
    p1 = subprocess.run(
        [MM2, "-ax", "splice:hq", "-N", "50", str(ref_fa), str(fq)],
        capture_output=True, text=True,
    )
    print(p1.stderr[-2000:], file=sys.stderr)
    p1.check_returncode()
    bam = pathlib.Path(f"{out_bam_prefix}.bam")
    bam.write_text(p1.stdout)
    sorted_bam = pathlib.Path(f"{out_bam_prefix}.sorted.bam")
    _run([SAMTOOLS, "sort", "-o", str(sorted_bam), str(bam)])
    _run([SAMTOOLS, "index", str(sorted_bam)])
    return sorted_bam


def _run_copy_assign(bam, ref_fa, region, out_prefix, extra_args=()):
    copy_assign = str(
        pathlib.Path(__file__).resolve().parents[2] / "target/release/copy_assign"
    )
    cmd = [
        copy_assign, "--bam", str(bam), "--fasta", str(ref_fa), "--region", region,
        "--out", str(out_prefix), "--min-copies", "2", "--absent-copies", "--linearize",
        *extra_args,
    ]
    _run(cmd)


def _divergent_arm(work, preset_used):
    """Plant a 4-member near-identical tandem cluster (see module docstring for why 4, not 3):
    A (anchor, div 0.02 from a random seed) + a B-cluster of 3 (B_host + 2 near-identical
    siblings, div 0.01 from B). Reference keeps A and B_host (a real 2-copy family); the two
    B-siblings are deleted from the reference and their reads pool onto B_host's locus, giving
    the host+2-hidden = 3 identifiable clusters the admission gate (`min_clusters=3`,
    hardcoded) requires. Recovery is verified at PSV-column resolution via `--dump-psv`
    (the only place the admitted candidates' reconstructed alleles are written to disk).
    """
    seedrng = random.Random(42)
    seed_seq = "".join(seedrng.choice(BASES) for _ in range(EXON_LEN * 2))
    # 4 members, each INDEPENDENTLY 2%-diverged from the same shared base (build_tandem's flat
    # per-copy model -- no phylogenetic nesting). At this divergence + this seed, copy2 and
    # copy3 both happen to be (very slightly) closer to copy1 than to copy0 -- confirmed
    # empirically (their reads pool 38/40 and 39/40 onto copy1's locus) -- which concentrates
    # both hidden members on ONE host locus, satisfying the min_clusters=3 admission gate.
    _, copies4 = build_tandem(seed_seq, 4, 0.02, seed_rng=11)
    A, B_host, hidden1, hidden2 = copies4
    truth_hidden = [hidden1, hidden2]

    intron_rng = random.Random(99)
    intron_a = _make_intron(intron_rng, INTRON_LEN)
    intron_b = _make_intron(intron_rng, INTRON_LEN)
    _make_intron(intron_rng, INTRON_LEN)  # advance the stream past copy2/copy3's (unused: deleted)
    _make_intron(intron_rng, INTRON_LEN)
    gen_a = _genomic_copy(A, intron_a)
    gen_b = _genomic_copy(B_host, intron_b)

    contig = "tandem_divergent"
    spacer = "N" * 200
    degraded = spacer.join([gen_a, gen_b])
    ref_fa = work / "divergent_ref.fa"
    _write_fasta(ref_fa, [(contig, degraded)])
    _run([SAMTOOLS, "faidx", str(ref_fa)])

    rng = random.Random(12)
    reads = _sim_reads([A, B_host, hidden1, hidden2], 40, 0.003, rng, "divergent")
    fq = work / "divergent_reads.fq"
    _write_fastq(fq, reads)

    sorted_bam = _align_and_index(ref_fa, fq, str(work / "divergent"))
    out_prefix = work / "divergent_out"
    region = f"{contig}:1-{len(degraded)}"
    _run_copy_assign(sorted_bam, ref_fa, region, out_prefix, extra_args=["--dump-psv", "--phase"])

    families = _read_tsv(f"{out_prefix}.families.tsv")
    dna_needs = _read_tsv(f"{out_prefix}.dna_needs.tsv")
    n_copies = int(families[0]["n_copies"]) if families else 0
    copy_number_correct = n_copies == 4  # A + B_host + 2 admitted hidden siblings

    cols = {int(r["col_index"]): int(r["genome_pos"]) for r in _read_tsv(f"{out_prefix}.psv_cols.tsv")}
    copy_rows = _read_tsv(f"{out_prefix}.psv_copies.tsv")
    copy_alleles = {int(r["copy_index"]): r["alleles"] for r in copy_rows}
    n_reference = 2  # copy_index 0 (A), 1 (B_host) are the reference copies
    admitted_idxs = sorted(i for i in copy_alleles if i >= n_reference)

    def genome_pos_to_body_offset(pos):
        # `psv_cols.tsv`'s `genome_pos` is reported in the FAMILY's shared reference frame (here,
        # copy0/A's own [0, len(gen_a)) coordinates) -- NOT a literal absolute genome coordinate
        # to be re-offset per copy. Since every copy shares the identical relative exon/intron
        # layout (EXON_LEN/INTRON_LEN/EXON_LEN), that frame is directly a valid LOCAL offset into
        # any copy's own body (confirmed empirically: every discovered column's genome_pos here
        # is < len(gen_a), i.e. always reported against copy0 regardless of which locus the
        # variation actually lives at).
        if pos < EXON_LEN:
            return pos
        if pos < EXON_LEN + INTRON_LEN:
            return None  # intron base -- no read ever covers this
        if pos < EXON_LEN + INTRON_LEN + EXON_LEN:
            return EXON_LEN + (pos - EXON_LEN - INTRON_LEN)
        return None  # outside one copy's genomic span

    # Columns actually anchored within B's locus (the ones informative about the admitted
    # candidates, whose host is B).
    b_cols = [ci for ci in range(len(cols)) if genome_pos_to_body_offset(cols[ci]) is not None]

    recovered_seqs = [
        "".join(copy_alleles[idx][ci] for ci in b_cols) for idx in admitted_idxs
    ]
    # Truth, projected onto the SAME B-locus PSV columns (so identity is computed at the
    # resolution the certificate actually operates at -- between PSV columns every copy is,
    # by the method's own construction, an allele-overlay of the host and therefore trivially
    # identical there).
    truth_projected = [
        "".join(t[genome_pos_to_body_offset(cols[ci])] for ci in b_cols) for t in truth_hidden
    ]

    matches = set_match(recovered_seqs, truth_projected) if recovered_seqs else []
    min_identity = min((m[2] for m in matches), default=0.0)

    return {
        "preset": preset_used,
        "n_copies_planted": 4,
        "n_reference_copies": 2,
        "n_hidden_copies": 2,
        "n_copies_detected": n_copies,
        "copy_number_correct": copy_number_correct,
        "n_dna_needs": len(dna_needs),
        "n_recovered": len(admitted_idxs),
        "min_identity": min_identity,
        "set_match": matches,
        "identity_resolution": "PSV-column (the certificate's native resolution; see --dump-psv)",
        "note": (
            "4-member design (2 reference + 2 hidden collapsed at one locus), not a bare "
            "3-of-3 deletion -- see module docstring finding #2: min_clusters=3 is hardcoded "
            "and a 3-member/1-deleted design can never reach it."
        ),
    }


def _identical_arm(work, preset_used):
    """Three byte-identical copies (0% divergence), ALL kept in the reference (deleting one
    is meaningless when they're indistinguishable). Tests whether copy-NUMBER is still
    recovered from pure read-conflict topology (n_copies in .families.tsv) even though there
    are zero PSVs to assign reads with -- and confirms every read genuinely abstains (Tied)."""
    seedrng = random.Random(42)
    seed_seq = "".join(seedrng.choice(BASES) for _ in range(EXON_LEN * 2))
    _, copies = build_tandem(seed_seq, 3, 0.0, seed_rng=11)
    assert copies[0] == copies[1] == copies[2]

    intron_rng = random.Random(99)
    introns = [_make_intron(intron_rng, INTRON_LEN) for _ in range(3)]
    genomics = [_genomic_copy(copies[ci], introns[ci]) for ci in range(3)]

    contig = "tandem_identical"
    spacer = "N" * 200
    full_ref = spacer.join(genomics)
    ref_fa = work / "identical_ref.fa"
    _write_fasta(ref_fa, [(contig, full_ref)])
    _run([SAMTOOLS, "faidx", str(ref_fa)])

    rng = random.Random(12)
    reads = _sim_reads(copies, 40, 0.003, rng, "identical")
    fq = work / "identical_reads.fq"
    _write_fastq(fq, reads)

    sorted_bam = _align_and_index(ref_fa, fq, str(work / "identical"))
    out_prefix = work / "identical_out"
    region = f"{contig}:1-{len(full_ref)}"
    _run_copy_assign(sorted_bam, ref_fa, region, out_prefix)

    families = _read_tsv(f"{out_prefix}.families.tsv")
    n_copies = int(families[0]["n_copies"]) if families else 0
    psv_cols = int(families[0]["psv_cols"]) if families else -1
    copy_number_correct = n_copies == 3

    assign_rows = _read_tsv(f"{out_prefix}.assignments.tsv")
    n_total = len(assign_rows)
    n_tied = sum(1 for r in assign_rows if r["status"] == "tied")
    assign_abstained = n_total > 0 and n_tied == n_total

    return {
        "preset": preset_used,
        "n_copies_planted": 3,
        "n_copies_detected": n_copies,
        "copy_number_correct": copy_number_correct,
        "psv_cols_found": psv_cols,
        "n_reads": n_total,
        "n_tied": n_tied,
        "assign_abstained": assign_abstained,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch", default="/home/juanfra/winloci_scratch/mech_sim")
    ap.add_argument(
        "--out",
        default=str(pathlib.Path(__file__).resolve().parents[2] / "bench/mechanism/verification_results.json"),
    )
    a = ap.parse_args()
    work = pathlib.Path(a.scratch)
    work.mkdir(parents=True, exist_ok=True)

    preset_used = "splice:hq"  # -N 50 throughout; see module docstring finding #1 for why
    divergent = _divergent_arm(work, preset_used)
    identical = _identical_arm(work, preset_used)

    result = {
        "divergent": {
            "copy_number_correct": divergent["copy_number_correct"],
            "n_recovered": divergent["n_recovered"],
            "min_identity": divergent["min_identity"],
            "preset": divergent["preset"],
            "design": divergent["note"],
            "n_copies_planted": divergent["n_copies_planted"],
            "n_copies_detected": divergent["n_copies_detected"],
            "n_dna_needs": divergent["n_dna_needs"],
            "identity_resolution": divergent["identity_resolution"],
        },
        "identical": {
            "copy_number_correct": identical["copy_number_correct"],
            "assign_abstained": identical["assign_abstained"],
            "preset": identical["preset"],
            "n_copies_planted": identical["n_copies_planted"],
            "n_copies_detected": identical["n_copies_detected"],
            "psv_cols_found": identical["psv_cols_found"],
            "n_reads": identical["n_reads"],
            "n_tied": identical["n_tied"],
        },
    }
    out = pathlib.Path(a.out)
    data = json.loads(out.read_text()) if out.exists() else {}
    data["v4a"] = result
    out.write_text(json.dumps(data, indent=2))
    print("V4a:", json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
