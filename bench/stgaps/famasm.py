#!/usr/bin/env python3
"""FAMASM: multi-copy assembly mode prototype (spec docs/superpowers/specs/2026-10-31-multicopy-assembly-design.md,
plan docs/superpowers/plans/2026-10-31-multicopy-assembly.md). Family assembly, then PSV placement, layered on BASE.

  famasm.py pool  COPIES.tsv BAM MOLECULES.tsv OUT_DIR [--contig C]   -> OUT_DIR/pool/<fam>.tsv, pool.stats.tsv
  famasm.py lift  COPIES.tsv GENOME.fa OUT_DIR                        -> OUT_DIR/lift/<fam>.tsv, lift.unmapped.tsv
  famasm.py admit OUT_DIR ASSIGNMENTS.tsv                             -> OUT_DIR/admit/<fam>.tsv, admit.stats.tsv
  famasm.py gtf   OUT_DIR BASE.gtf COPIES.tsv OUT.gtf                 -> FAMASM GTF (BASE verbatim + admitted), gtf.stats.tsv
  famasm.py score OUT_DIR                                             -> verdict.txt
Pool = primary records + secondary records with AS >= 0.98 x the molecule's genome-wide best (as_table). One molecule once per
family. Nothing is emitted at a copy without reads O2 certified to it (status `assigned`)."""
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fw_layer  # noqa: E402

GOOD_RATIO = 0.98


def good_secondary(as_, best_as):
    return as_ >= GOOD_RATIO * best_as


def load_best_as(path, names):
    """best genome-wide AS of each molecule in `names` from PREFIX.molecules.tsv (first two columns)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line[0] == "#" or line.startswith("name\t"):
                continue
            f = line.split("\t", 2)
            if f[0] in names:
                out[f[0]] = int(f[1])
    return out


def pool(copies_path, bam_path, molecules_path, out_dir, contig=None):
    import pysam
    copies = fw_layer.read_copies(copies_path)
    if contig:
        copies = {k: v for k, v in copies.items() if v[0] == contig}
    fams = defaultdict(list)
    for (fam, idx), v in copies.items():
        fams[fam].append((idx, v))
    os.makedirs(os.path.join(out_dir, "pool"), exist_ok=True)
    bam = pysam.AlignmentFile(bam_path)
    stats = []
    for fam, members in sorted(fams.items()):
        if len(members) < 2:
            continue
        recs = []  # (read, idx, is_primary, AS, chain)
        for idx, (chrom, start, end, _) in sorted(members, key=lambda m: int(m[0])):
            for rec in bam.fetch(chrom, start, end):
                if rec.is_unmapped or rec.is_supplementary or rec.is_qcfail or not rec.has_tag("AS"):
                    continue
                # a record belongs to the copy whose locus it overlaps most (loci of one family do not overlap)
                ov = min(rec.reference_end, end) - max(rec.reference_start, start)
                if ov <= 0:
                    continue
                recs.append((rec.query_name, idx, not rec.is_secondary, int(rec.get_tag("AS")), fw_layer.chain_of(rec)))
        best = load_best_as(molecules_path, {r[0] for r in recs if not r[2]})
        kept, seen = [], set()
        secondary_kept = 0
        for read, idx, prim, as_, chain in recs:
            if not prim and not good_secondary(as_, best.get(read, as_ + 1)):
                continue
            key = (read, idx)
            if key in seen:
                continue
            seen.add(key)
            kept.append((read, idx, prim, as_, chain))
            secondary_kept += not prim
        with open(os.path.join(out_dir, "pool", f"{fam}.tsv"), "w") as w:
            w.write("read\tcopy_idx\tis_primary\tAS\tchain\n")
            for read, idx, prim, as_, chain in kept:
                w.write(f"{read}\t{idx}\t{int(prim)}\t{as_}\t{fw_layer.fmt(chain) if chain else ''}\n")
        stats.append((fam, len({k[0] for k in kept}), len(kept), secondary_kept))
    with open(os.path.join(out_dir, "pool.stats.tsv"), "w") as w:
        w.write("family\tmolecules\trecords\tsecondary_kept\n")
        for row in stats:
            w.write("\t".join(map(str, row)) + "\n")


def lift_rows(copies, reps, lifts, fam, rows):
    """rows = (read, idx, is_primary, chain) -> (read, idx, is_primary, fam_chain_str, private).
    private = 1 when any junction has no lifted column (the chain stays in copy coordinates as `c<idx>:...`)."""
    out = []
    for read, idx, prim, chain in rows:
        fc = fw_layer.fam_chain(copies, reps, lifts, fam, idx, chain) if chain else ()
        if chain and fc is None:
            out.append((read, idx, prim, f"c{idx}:{fw_layer.fmt(chain)}", 1))
        else:
            out.append((read, idx, prim, fw_layer.fmt(fc), 0))
    return out


def lift(copies_path, fa_path, out_dir):
    import subprocess
    copies = fw_layer.read_copies(copies_path)
    reps = fw_layer.reps_of(copies)
    mdir = os.path.join(out_dir, "map")
    fw_layer.map_(copies_path, fa_path, mdir)
    paf = os.path.join(mdir, "copies_vs_reps.paf")
    with open(paf, "w") as w:
        subprocess.run(["minimap2", "-x", "asm20", "-c", "--cs", "-t", os.environ.get("FAMASM_THREADS", "4"),
                        os.path.join(mdir, "reps.fa"), os.path.join(mdir, "copies.fa")], stdout=w, check=True,
                       stderr=subprocess.DEVNULL)
    lifts = fw_layer.load_lift(paf)
    os.makedirs(os.path.join(out_dir, "lift"), exist_ok=True)
    unmapped = []
    for fam, rep in sorted(reps.items()):
        src = os.path.join(out_dir, "pool", f"{fam}.tsv")
        if not os.path.exists(src):
            continue
        for (f, idx) in copies:
            if f == fam and idx != rep and (fam, idx) not in lifts:
                unmapped.append((fam, idx))
        rows = []
        with open(src) as fh:
            next(fh)
            for line in fh:
                read, idx, prim, _as, chain = line.rstrip("\n").split("\t")
                rows.append((read, idx, int(prim), fw_layer.parse(chain) if chain else ()))
        with open(os.path.join(out_dir, "lift", f"{fam}.tsv"), "w") as w:
            w.write("read\tcopy_idx\tis_primary\tfam_chain\tprivate\n")
            for r in lift_rows(copies, reps, lifts, fam, rows):
                w.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(out_dir, "lift.unmapped.tsv"), "w") as w:
        w.write("family\tcopy_idx\n")
        for fam, idx in sorted(unmapped):
            w.write(f"{fam}\t{idx}\n")


def _unlift(lf, t):
    """target (representative-local) coordinate -> query (copy-local) coordinate, or None."""
    _, blocks = lf
    for lo, hi, t0, s in blocks:
        if t0 <= t < t0 + (hi - lo):
            return lo + (t - t0) if s == "+" else hi - 1 - (t - t0)
    return None


def inverse_lift(lifts, copies, reps, fam, idx, fam_chain):
    """family (representative-local) junction chain -> genomic chain at copy idx; None when any junction has no column."""
    c, a, b, _ = copies[(fam, idx)]
    out = []
    for p, q in fam_chain:
        if reps.get(fam) == idx:
            u, v = p, q
        else:
            lf = lifts.get((fam, idx))
            if lf is None:
                return None
            u, v = _unlift(lf, p), _unlift(lf, q)
            if u is None or v is None:
                return None
        x, y = min(u, v) + 1 + a, max(u, v) + a
        if y <= x:
            return None
        out.append((x, y))
    return tuple(sorted(out))


def main():
    cmd, *a = sys.argv[1:]
    if cmd == "pool":
        contig = None
        if "--contig" in a:
            i = a.index("--contig")
            contig = a[i + 1]
            del a[i:i + 2]
        pool(*a[:4], contig=contig)
    elif cmd == "lift":
        lift(*a[:3])
    else:
        sys.exit(f"unknown command {cmd}")


if __name__ == "__main__":
    main()
