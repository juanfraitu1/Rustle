"""V4b Phase 1: hunt for real gorilla tandem families in the "recoverable window".

Selects families with >=3 co-located copies (single chrom, span a few Mb) from the genome-wide
catalog (GGO_gwcat.*), computes pairwise FULL-SPAN identity per family via all-vs-all
`minimap2 -c --eqx -x asm20` on the extracted copy sequences (identity = 1 - de:f: divergence
tag), and ranks candidates by:
  1. median pairwise identity in [0.96, 0.979] (in-window: diverged enough to clear the
     remap_max_identity=0.98 admission gate, similar enough to collapse into >=3 clusters)
  2. tightness of co-location (smaller span = better collapsibility)
  3. copy count (>=3, prefer 4-8)

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
SCRATCH = pathlib.Path("/home/juanfra/winloci_scratch/mech_real/hunt")
SPAN_CAP = 3_000_000  # "a few Mb" tandem-array cap
MIN_COPIES = 3

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
        })
    return cands


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


def main():
    fam = load_families()
    copies = load_copies()
    cands = select_candidates(fam, copies)
    print(f"Single-chrom, span<={SPAN_CAP}, n_copies>={MIN_COPIES}: {len(cands)} candidate families",
          file=sys.stderr)

    seqs = load_fasta_index()
    SCRATCH.mkdir(parents=True, exist_ok=True)

    rows = []
    for c in sorted(cands, key=lambda x: x["span_bp"]):
        med, npairs = pairwise_identity(c["family_id"], c["copy_idxs"], seqs, SCRATCH)
        in_window = med is not None and 0.96 <= med <= 0.979
        rows.append({
            "family_id": c["family_id"], "n_copies": c["n_copies"], "contig": c["contig"],
            "span_bp": c["span_bp"], "median_identity": med, "n_pairs": npairs,
            "in_window": in_window, "start": c["start"], "end": c["end"],
        })
        print(f"  {c['family_id']:>10} n_copies={c['n_copies']:>3} span={c['span_bp']:>10} "
              f"median_identity={med} n_pairs={npairs} in_window={in_window}", file=sys.stderr)

    # rank: in_window first, then tighter span, then higher copy count (prefer 4-8)
    def rank_key(r):
        pref_count = 0 if (r["n_copies"] is not None and 4 <= r["n_copies"] <= 8) else 1
        return (0 if r["in_window"] else 1, pref_count, r["span_bp"])

    rows.sort(key=rank_key)

    with open(OUT_TSV, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["family_id", "n_copies", "contig", "span_bp", "median_identity", "n_pairs",
                    "in_window", "start", "end"])
        for r in rows:
            w.writerow([r["family_id"], r["n_copies"], r["contig"], r["span_bp"],
                        r["median_identity"], r["n_pairs"], r["in_window"], r["start"], r["end"]])

    n_in_window = sum(1 for r in rows if r["in_window"])
    n_with_identity = sum(1 for r in rows if r["median_identity"] is not None)
    print(f"\nWrote {len(rows)} candidates -> {OUT_TSV}", file=sys.stderr)
    print(f"in_window (0.96-0.979): {n_in_window}", file=sys.stderr)
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
