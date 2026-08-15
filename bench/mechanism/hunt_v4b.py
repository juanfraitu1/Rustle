"""V4b Phase 1: hunt for real gorilla tandem families in the "recoverable window".

Selects families with >=3 co-located copies (single chrom, span a few Mb) from the genome-wide
catalog (GGO_gwcat.*) and scores each family for whether a copy-removal-and-recovery experiment
could possibly succeed on it.

TWO STATISTICS ARE EMITTED. They answer different questions and they disagree.

`median_identity` (LEGACY, kept byte-reproducible so earlier results stay checkable)
    Median over unordered copy PAIRS of best-hit identity, all-vs-all `-x asm20`.
    ⚠ This does NOT match the gate it was meant to predict, in three ways:
      * it is a PAIR statistic, but gate 5 is an ANY/BEST-HIT rule, so a family whose pair-median
        looks comfortably divergent can still contain a near-identical pair that the gate rejects
        (GWFAM438: pair-median 0.9127, yet one copy sits at 0.9969);
      * it scores a FAMILY, but the experiment deletes a COPY -- and the experimenter chooses
        which copy, so the family-level question is "is there ANY good target", not "what is
        the family's average";
      * it uses `-x asm20`, which this project retired for being a subset of the shipped
        `sensitive` tier and whose ~80% seeding floor silently reports "no alignment" for
        families that do align (GWFAM269/440/16/454/320 return 0 pairs here but 6/3/3/10/2 pairs
        at 0.71-0.76 under `-k11 -w5`).

`per-copy gate identity` (CURRENT) -- built to RESEMBLE `absent_copy.rs::remap_identity_minimap2`:
    same argv (`minimap2 -cx splice --eqx --secondary=no -t 1`), same best-alignment-BY-BLOCK-LENGTH
    rule with the same strict-`>` first-wins tie-break, same `1 - de:f:` with a matches/blocklen
    fallback.
    ⚠ IT IS NOT THE SAME QUANTITY, and the differences are NOT only scope:
      * QUERY OBJECT: gate 5's query is the host-derived RECONSTRUCTION (host spliced sequence with
        PSV alleles overlaid); this screen's query is the copy's own real sequence. Measured on
        GWFAM247, the reconstruction carries 16 mismatches / 4134 bp where the true copy-0-vs-copy-1
        spliced alignment carries 88 / 3256 -- it reproduces ~18% of the real divergence, so its
        remap identity sits near 1.0 however divergent the true copy is. THIS is why the screen
        cannot predict the gate, and no change of target scope would fix it.
      * TARGET TYPE: this screen aligns spliced-vs-spliced; gate 5 aligns spliced-vs-GENOMIC.
      * SELF-EXCLUSION: this screen deletes the query copy from the targets; gate 5 has no such rule.
      * TARGET SCOPE: in the V4b driver, gate 5's `--fasta` is the ~150 kb LOCAL masked window
        (v4b_3copy.py), not the whole genome.
    ⟹ treat this as a SCREEN for which copies are structurally hopeless, never as a predictor of the
    identity gate 5 will compute. Measured: 0 of 4 ablations agreed; the screen was optimistic in
    3/3 measurable comparisons by a mean +0.0215, against a window only 0.0200 wide. A copy is recoverable iff its gate identity is < remap_max_identity
    (0.98); at or above that, gate 5 rejects any reconstruction of it no matter how good.
    Derived per family. A usable deletion target must clear the gate AND still be similar enough
    that the deleted copy's reads collapse onto a survivor -- with nothing to collapse onto there
    is no ambiguity to resolve and nothing to detect:
      max_admissible_identity = MAX identity still < 0.98, over ALL admissible copies.
      best_target_identity    = the same, but restricted to copies with >= MIN_TARGET_READS reads,
                             i.e. the copy actually nominated. None when no expressed copy is
                             admissible. (Taking the MIN instead would nominate the family's most
                             DIVERGENT copy, which clears the gate trivially and is the least
                             recoverable one it has.) Ties are broken by read count then copy id,
                             never by catalog row order.
      n_usable             = #copies in [WINDOW_LO, 0.98)
      n_usable_with_reads  = those also carrying >= MIN_TARGET_READS reads; this is
                             what `in_window` keys on, since an unexpressed target
                             cannot be recovered no matter where its identity sits
      n_blocked            = #copies >= 0.98      -> gate 5 can never admit these
      n_too_divergent      = #copies < WINDOW_LO  -> nothing to collapse onto

⚠ WINDOW_LO is a HEURISTIC, unlike the 0.98 upper bound which is the gate's own constant. The run
prints a sensitivity sweep over it; quote the sweep, not a single count.

⚠ SCOPE: this screen remaps against the family's other copies only. Widening to the whole genome is
NOT a one-sided correction -- measured over 31 copies, 14 go UP, 9 are unchanged and 6 go DOWN (to
-0.1431), because the selection rule maximises BLOCK LENGTH, not identity: a longer, more divergent
alignment can displace a shorter near-identical one (GWFAM227 copy 2: bl 1222 @ 0.8198 beats bl 1209
@ 0.9603 by 13 bp). No in_window verdict moved for the ranked families, but 1 of 31 copies is blocked
by an out-of-family paralog on a DIFFERENT chromosome (GWFAM441 copy 2: 0.7788 -> 0.9900).

Ranks by: has a usable target, then more usable copies, then tighter span.
Writes bench/mechanism/v4b_hunt_shortlist.tsv.
"""
import csv
import pathlib
import statistics
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
DATA = pathlib.Path("/mnt/linuxdisk/home/juanfraitu/winloci_data")
COPIES_TSV = DATA / "GGO_gwcat.copies.tsv"
FAMILIES_TSV = DATA / "GGO_gwcat.families.tsv"
COPIES_FA = DATA / "GGO_gwcat.copies.fa"
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
SAMTOOLS = "/home/juanfra/miniforge3/bin/samtools"
# The FULL BAM. /mnt/linuxdisk/.../GGO_ds.bam is a `samtools view -s 0.37` DOWNSAMPLE -- using it
# would deflate every read count by ~3x.
BAM = pathlib.Path("/home/juanfra/winloci_scratch/GGO_mm.bam")
SCRATCH = pathlib.Path("/home/juanfra/winloci_scratch/mech_real/hunt")
SPAN_CAP = 3_000_000  # "a few Mb" tandem-array cap
MIN_COPIES = 3
# Mirrors AbsentCopyParams::remap_max_identity (src/rustle/vg_family/absent_copy.rs).
# A copy whose gate identity reaches this is rejected by gate 5 however well it is reconstructed.
GATE_MAX_IDENTITY = 0.98
WINDOW_LO, WINDOW_HI = 0.96, 0.979
# A deletion target needs transcript evidence, not just well-placed identity.
# ⚠ DO NOT justify this with the V4b GWFAM227 attempt. That run's own artifact
# (v4b_mincl2/GWFAM227/ablation.ca.dna_needs.tsv) records ">=98% remap identity (paralog-leak or
# het)" -- it failed at gate 5 on DIVERGENCE (re-measured 0.9978), not on read support. Its
# read_count=3 is a cluster at the HOST locus, and the control records the same 3 at the same span
# with the copy still present, so its match to copies.tsv n_reads=3 is a coincidence.
# The honest justification is empirical only: the in-window set is unchanged for any value in
# [8, 18] (measured plateau), so the threshold is underived but not finely tuned.
MIN_TARGET_READS = 10

OUT_TSV = ROOT / "bench/mechanism/v4b_hunt_shortlist.tsv"


def load_families():
    with open(FAMILIES_TSV) as f:
        return {row["family_id"]: row for row in csv.DictReader(f, delimiter="\t")}


def load_copies():
    out = {}
    with open(COPIES_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            out.setdefault(row["family_id"], []).append(row)
    return out


def load_fasta_index():
    """family_id|copy_idx -> sequence, streamed once over the multi-fasta."""
    seqs = {}
    header = None
    buf = []
    with open(COPIES_FA) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(buf)
                parts = line[1:].split("|")
                header = f"{parts[0]}|{parts[1]}"
                buf = []
            else:
                buf.append(line)
        if header is not None:
            seqs[header] = "".join(buf)
    return seqs


def select_candidates(fam, copies):
    cands = []
    for fid, frow in fam.items():
        n = int(frow["n_copies"])
        if n < MIN_COPIES:
            continue
        if int(frow["n_chroms"]) != 1:
            continue  # exclude cross-chrom-dispersed families
        cps = copies.get(fid, [])
        if len(cps) < MIN_COPIES:
            continue
        chrom = cps[0]["chrom"]
        starts = [int(c["start"]) for c in cps]
        ends = [int(c["end"]) for c in cps]
        span = max(ends) - min(starts)
        if span > SPAN_CAP:
            continue
        cands.append({
            "family_id": fid, "n_copies": n, "contig": chrom,
            "span_bp": span, "start": min(starts), "end": max(ends),
            "copy_idxs": [c["copy_idx"] for c in cps],
            # Sequence identity alone is not enough to pick a deletion target: a copy with almost
            # no transcript evidence cannot be recovered however well-placed its identity is.
            # ⚠ copies.tsv `n_reads` UNDERCOUNTS primary reads over the copy span on 105/105 copies
            # (median 5.1x, max 205x) and misclassifies 36.2% of copies against MIN_TARGET_READS,
            # so it is kept only as a fallback -- `bam_read_counts` is the real source.
            "reads_by_copy": {c["copy_idx"]: int(c["n_reads"]) for c in cps},   # fallback only
            "spans_by_copy": {c["copy_idx"]: (int(c["start"]), int(c["end"])) for c in cps},
        })
    return cands


def bam_read_counts(cands, workdir):
    """Primary reads CONTAINED in each copy's span, counted from the full BAM.

    Replaces copies.tsv `n_reads`, which undercounts primary reads over the copy span on 105/105
    shortlist copies (median 5.1x, max 205x) and misclassifies 36.2% of them against
    MIN_TARGET_READS -- enough to swing whole families in and out of the shortlist.

    CONTAINED, not overlapping: an overlap rule once let a neighbouring locus's read chain be
    adopted as a locus's own evidence. `-F 2308` keeps primaries only.

    Cached to `read_counts.tsv` in `workdir`; delete that file to recount.
    """
    cache = workdir / "read_counts.tsv"
    cached = {}
    if cache.exists():
        with open(cache) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                cached[(row["family_id"], row["copy_idx"])] = int(row["n_primary_contained"])
    out, fresh = {}, []
    for c in cands:
        for ci, (s0, e0) in c["spans_by_copy"].items():
            key = (c["family_id"], ci)
            if key in cached:
                out[key] = cached[key]
                continue
            r = subprocess.run(
                [SAMTOOLS, "view", "-F", "2308", str(BAM),
                 f"{c['contig']}:{s0 + 1}-{e0}"],
                capture_output=True, text=True,
            )
            n = 0
            for line in r.stdout.splitlines():
                f_ = line.split("\t")
                if len(f_) < 6:
                    continue
                pos = int(f_[3]) - 1                      # SAM POS is 1-based; spans are 0-based
                ref_len = 0
                num = ""
                for ch in f_[5]:                          # CIGAR ops consuming REFERENCE
                    if ch.isdigit():
                        num += ch
                    else:
                        if ch in "MDN=X":                 # N = intron spliced out, still spans ref
                            ref_len += int(num)
                        num = ""
                if pos >= s0 and pos + ref_len <= e0:
                    n += 1
            out[key] = n
            fresh.append((c["family_id"], ci, n))
    if fresh:
        write_header = not cache.exists()
        with open(cache, "a", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            if write_header:
                w.writerow(["family_id", "copy_idx", "n_primary_contained"])
            w.writerows(fresh)
    return out


def pairwise_identity(fid, copy_idxs, seqs, workdir):
    """All-vs-all minimap2 asm20 --eqx on this family's copies; returns median identity
    over off-diagonal best-hit pairs (1 - de:f: divergence), or None if too few pairs align."""
    fa_path = workdir / f"{fid}.fa"
    with open(fa_path, "w") as f:
        for ci in copy_idxs:
            key = f"{fid}|{ci}"
            seq = seqs.get(key)
            if seq is None:
                continue
            f.write(f">{key}\n{seq}\n")
    r = subprocess.run(
        # -X: all-vs-all overlap mode -- without it minimap2 reports only the single best
        # target per query (usually the trivial self-hit), hiding every paralog pair.
        [MM2, "-c", "--eqx", "-X", "-x", "asm20", "-N", "50", str(fa_path), str(fa_path)],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        print(f"  [{fid}] minimap2 failed: {r.stderr[-500:]}", file=sys.stderr)
        return None, 0
    best = {}  # unordered pair -> best (highest matching-bases) identity
    for line in r.stdout.splitlines():
        f_ = line.split("\t")
        qname, tname = f_[0], f_[5]
        if qname == tname:
            continue
        pair = tuple(sorted((qname, tname)))
        nmatch = int(f_[9])
        de = None
        for tag in f_[12:]:
            if tag.startswith("de:f:"):
                de = float(tag[5:])
        if de is None:
            continue
        ident = 1.0 - de
        cur = best.get(pair)
        if cur is None or nmatch > cur[0]:
            best[pair] = (nmatch, ident)
    idents = [v[1] for v in best.values()]
    if not idents:
        return None, 0
    return statistics.median(idents), len(idents)


def copy_gate_identities(fid, copy_idxs, seqs, workdir):
    """Per-copy gate-5 identity: for each copy, the identity of its BEST alignment (by alignment
    block length) to any OTHER copy of the same family.

    Reproduces `absent_copy.rs::remap_identity_minimap2` at the command level -- same preset
    (`-cx splice`), same `--eqx --secondary=no -t 1`, same best-by-block-length selection, same
    `1 - de:f:` with a matches/blocklen fallback. The only difference is scope: the target here is
    the family's other copies rather than the whole genome, which makes the result a LOWER bound
    (see module docstring).

    Returns {copy_idx: identity or None}. None means "no alignment to any sibling", which is what
    gate 5 reports as `no homology on remap` -- a distinct failure from being too similar.
    """
    present = [ci for ci in copy_idxs if f"{fid}|{ci}" in seqs]
    out = {}
    for ci in present:
        others = [c for c in present if c != ci]
        if not others:
            out[ci] = None
            continue
        qp = workdir / f"{fid}.{ci}.q.fa"
        tp = workdir / f"{fid}.{ci}.t.fa"
        with open(qp, "w") as f:
            f.write(f">{fid}|{ci}\n{seqs[f'{fid}|{ci}']}\n")
        with open(tp, "w") as f:
            for cj in others:
                f.write(f">{fid}|{cj}\n{seqs[f'{fid}|{cj}']}\n")
        r = subprocess.run(
            [MM2, "-cx", "splice", "--eqx", "--secondary=no", "-t", "1", str(tp), str(qp)],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            print(f"  [{fid}|{ci}] minimap2 failed: {r.stderr[-300:]}", file=sys.stderr)
            out[ci] = None
            continue
        best_id, best_span = None, 0
        for line in r.stdout.splitlines():
            f_ = line.split("\t")
            if len(f_) < 12:
                continue
            nmatch, blocklen = int(f_[9]), int(f_[10])
            if blocklen <= 0:
                continue
            de = None
            for tag in f_[12:]:
                if tag.startswith("de:f:"):
                    de = float(tag[5:])
            ident = (1.0 - de) if de is not None else (nmatch / blocklen)
            if blocklen > best_span:
                best_span, best_id = blocklen, ident
        out[ci] = best_id
    return out


def _copy_sort_key(c):
    """Copy ids are numeric strings in this catalog; fall back to string order if that changes."""
    try:
        return (0, int(c))
    except (TypeError, ValueError):
        return (1, str(c))


def summarize_gate(per_copy, reads_by_copy=None):
    """Family-level roll-up of per-copy gate identities.

    The experimenter CHOOSES which copy to delete, so the family-level question is "does this
    family offer even ONE usable target", not "what is the family's average". A usable target is
    squeezed between two opposing requirements:

      * HIGH enough that the deleted copy's reads still collapse onto a surviving paralog -- with
        nothing to collapse onto there is no ambiguity to resolve and nothing to detect;
      * BELOW `GATE_MAX_IDENTITY` (0.98), or gate 5 rejects the reconstruction outright.

    So `best_target_identity` is the MAXIMUM identity that still clears the gate -- the most
    collapsible admissible copy. Taking the minimum instead would nominate the family's most
    DIVERGENT copy, which clears the gate trivially and is the least recoverable copy it has.
    """
    hit = {c: v for c, v in per_copy.items() if v is not None}
    reads_by_copy = reads_by_copy or {}
    empty = {"best_target_copy": None, "best_target_identity": None, "worst_copy_identity": None,
             "n_copies_hit": 0, "n_usable": 0, "n_blocked": 0, "n_too_divergent": 0,
             "best_target_n_reads": None, "n_usable_with_reads": 0,
             "max_admissible_identity": None}
    if not hit:
        return empty
    admissible = {c: v for c, v in hit.items() if v < GATE_MAX_IDENTITY}
    # Two DISTINCT columns, because collapsing them into one made a single column carry two
    # different definitions across rows with nothing marking which rule produced a given row:
    #   max_admissible_identity -- the plain max over admissible copies, matching the docstring
    #   best_target_*           -- the pick actually nominated, restricted to expressed copies
    max_adm = max(admissible.values()) if admissible else None
    with_reads = {c: v for c, v in admissible.items()
                  if reads_by_copy.get(c, 0) >= MIN_TARGET_READS}
    # ⚠ EXPLICIT, TOTAL ORDER. `max()` alone returns the first maximum in dict insertion order
    # (= catalog row order), and bit-exact ties are common: `de:f:` prints to 4 decimals and
    # reciprocal best hits between two copies return the SAME de by construction, so 8/24 families
    # tie at the argmax. GWFAM247 copies 0/1 are both 0.9754 and GWFAM227 copies 1/2 are both
    # 0.9629 -- without this, which copy gets nominated is decided by row order.
    # highest identity, then most reads, then lowest copy id -- a total order with no ties left
    best_c = min(with_reads,
                 key=lambda c: (-with_reads[c], -reads_by_copy.get(c, 0), _copy_sort_key(c))
                 ) if with_reads else None
    return {
        "best_target_copy": best_c,
        "best_target_identity": with_reads[best_c] if best_c is not None else None,
        "max_admissible_identity": max_adm,
        "best_target_n_reads": reads_by_copy.get(best_c) if best_c is not None else None,
        "n_usable_with_reads": sum(
            1 for c, v in hit.items()
            if WINDOW_LO <= v < GATE_MAX_IDENTITY and reads_by_copy.get(c, 0) >= MIN_TARGET_READS),
        "worst_copy_identity": max(hit.values()),
        "n_copies_hit": len(hit),
        # usable = clears the gate AND is still similar enough to collapse
        "n_usable": sum(1 for v in hit.values() if WINDOW_LO <= v < GATE_MAX_IDENTITY),
        "n_blocked": sum(1 for v in hit.values() if v >= GATE_MAX_IDENTITY),
        "n_too_divergent": sum(1 for v in hit.values() if v < WINDOW_LO),
    }


def main():
    fam = load_families()
    copies = load_copies()
    cands = select_candidates(fam, copies)
    print(f"Single-chrom, span<={SPAN_CAP}, n_copies>={MIN_COPIES}: {len(cands)} candidate families",
          file=sys.stderr)

    seqs = load_fasta_index()
    SCRATCH.mkdir(parents=True, exist_ok=True)
    bam_counts = bam_read_counts(cands, SCRATCH)
    n_disagree = sum(1 for c in cands for ci, v in c["reads_by_copy"].items()
                     if (v >= MIN_TARGET_READS) != (bam_counts.get((c["family_id"], ci), 0)
                                                    >= MIN_TARGET_READS))
    n_copies_tot = sum(len(c["reads_by_copy"]) for c in cands)
    print(f"read source: BAM primary-contained counts "
          f"(catalog n_reads disagrees on {n_disagree}/{n_copies_tot} copies at the "
          f"MIN_TARGET_READS={MIN_TARGET_READS} boundary)", file=sys.stderr)

    rows = []
    for c in sorted(cands, key=lambda x: x["span_bp"]):
        med, npairs = pairwise_identity(c["family_id"], c["copy_idxs"], seqs, SCRATCH)
        per_copy = copy_gate_identities(c["family_id"], c["copy_idxs"], seqs, SCRATCH)
        reads = {ci: bam_counts.get((c["family_id"], ci), c["reads_by_copy"].get(ci, 0))
                 for ci in c["reads_by_copy"]}
        g = summarize_gate(per_copy, reads)
        bt = g["best_target_identity"]
        # a usable target must clear the gate, still be collapsible, AND be expressed
        in_window = g["n_usable_with_reads"] > 0
        legacy_in_window = med is not None and WINDOW_LO <= med <= WINDOW_HI
        rows.append({
            "family_id": c["family_id"], "n_copies": c["n_copies"], "contig": c["contig"],
            "span_bp": c["span_bp"], "median_identity": med, "n_pairs": npairs,
            "in_window": in_window, "legacy_in_window": legacy_in_window,
            "start": c["start"], "end": c["end"],
            "per_copy_gate": ";".join(
                f"{ci}:{'NA' if per_copy[ci] is None else format(per_copy[ci], '.4f')}"
                for ci in sorted(per_copy, key=lambda x: int(x) if str(x).isdigit() else x)),
            **g,
        })
        print(f"  {c['family_id']:>10} n_copies={c['n_copies']:>3} span={c['span_bp']:>10} "
              f"legacy_median={med} best_target={bt} "
              f"usable={g['n_usable']}/{g['n_copies_hit']} (with_reads={g['n_usable_with_reads']}) "
              f"target_reads={g['best_target_n_reads']} blocked={g['n_blocked']} "
              f"too_divergent={g['n_too_divergent']} in_window={in_window}", file=sys.stderr)

    # rank: a usable deletion target first, then more attemptable copies, then tighter span
    def rank_key(r):
        return (0 if r["in_window"] else 1, -r["n_usable_with_reads"], -r["n_usable"], r["span_bp"])

    rows.sort(key=rank_key)

    cols = ["family_id", "n_copies", "contig", "span_bp", "best_target_identity",
            "best_target_copy", "best_target_n_reads", "max_admissible_identity",
            "n_usable", "n_usable_with_reads",
            "n_blocked", "n_too_divergent", "n_copies_hit", "in_window",
            "per_copy_gate", "worst_copy_identity", "median_identity", "n_pairs",
            "legacy_in_window", "start", "end"]
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(cols)
        for r in rows:
            w.writerow([r[k] for k in cols])

    n_in_window = sum(1 for r in rows if r["in_window"])
    n_legacy = sum(1 for r in rows if r["legacy_in_window"])
    n_with_identity = sum(1 for r in rows if r["median_identity"] is not None)
    n_gate = sum(1 for r in rows if r["best_target_identity"] is not None)
    print(f"\nWrote {len(rows)} candidates -> {OUT_TSV}", file=sys.stderr)
    print(f"in_window on the GATE statistic (best target {WINDOW_LO}-{WINDOW_HI}): {n_in_window}",
          file=sys.stderr)
    print(f"in_window on the LEGACY pair-median (for comparison): {n_legacy}", file=sys.stderr)
    print(f"families with a computable gate identity: {n_gate} "
          f"(legacy pair-median computable: {n_with_identity})", file=sys.stderr)
    blocked = [r for r in rows if r["n_blocked"] > 0]
    print(f"families containing at least one copy gate 5 can NEVER admit: {len(blocked)}",
          file=sys.stderr)
    for r in blocked:
        print(f"  {r['family_id']}: {r['n_blocked']} of {r['n_copies_hit']} copies >= "
              f"{GATE_MAX_IDENTITY}  [{r['per_copy_gate']}]  legacy pair-median "
              f"{r['median_identity']}", file=sys.stderr)
    # WINDOW_LO is a heuristic (unlike 0.98, which is the gate's own constant), so report the
    # in-window set's sensitivity to it rather than a single count.
    print("\nsensitivity of the usable-target set to the heuristic lower bound:", file=sys.stderr)
    all_copy_ids = {r["family_id"]: [
        None if v.split(":")[1] == "NA" else float(v.split(":")[1])
        for v in r["per_copy_gate"].split(";") if v
    ] for r in rows}
    for lo in (0.90, 0.92, 0.94, 0.95, 0.96, 0.97):
        fams = sorted(fid for fid, vs in all_copy_ids.items()
                      if any(v is not None and lo <= v < GATE_MAX_IDENTITY for v in vs))
        print(f"  lo={lo:.2f}: {len(fams):>2} families  {' '.join(fams)}", file=sys.stderr)

    print("\ngate-statistic distribution (best usable target per family):", file=sys.stderr)
    gbuckets = {"<0.90": 0, "0.90-0.95": 0, "0.95-0.96": 0, "0.96-0.979": 0, "0.979-0.98": 0,
                "no admissible copy": 0, "no alignment": 0}
    for r in rows:
        if r["n_copies_hit"] == 0:
            gbuckets["no alignment"] += 1
        elif r["best_target_identity"] is None:
            gbuckets["no admissible copy"] += 1
        else:
            v = r["best_target_identity"]
            k = ("<0.90" if v < 0.90 else "0.90-0.95" if v < 0.95 else "0.95-0.96" if v < 0.96
                 else "0.96-0.979" if v <= 0.979 else "0.979-0.98")
            gbuckets[k] += 1
    for k, v in gbuckets.items():
        print(f"  {k}: {v}", file=sys.stderr)

    print("\nlegacy pair-median distribution (for comparison):", file=sys.stderr)
    idents = [r["median_identity"] for r in rows if r["median_identity"] is not None]
    if idents:
        buckets = {"<0.90": 0, "0.90-0.95": 0, "0.95-0.96": 0, "0.96-0.979": 0, "0.98-0.995": 0, ">=0.995": 0}
        for i in idents:
            if i < 0.90:
                buckets["<0.90"] += 1
            elif i < 0.95:
                buckets["0.90-0.95"] += 1
            elif i < 0.96:
                buckets["0.95-0.96"] += 1
            elif i <= 0.979:
                buckets["0.96-0.979"] += 1
            elif i < 0.995:
                buckets["0.98-0.995"] += 1
            else:
                buckets[">=0.995"] += 1
        print(f"identity distribution ({n_with_identity} families with computable identity):", file=sys.stderr)
        for k, v in buckets.items():
            print(f"  {k}: {v}", file=sys.stderr)


if __name__ == "__main__":
    main()
