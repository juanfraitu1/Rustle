#!/usr/bin/env python3
"""FAMASM: multi-copy assembly mode prototype (spec docs/superpowers/specs/2026-10-31-multicopy-assembly-design.md,
plan docs/superpowers/plans/2026-10-31-multicopy-assembly.md). Family assembly, then PSV placement, layered on BASE.

  famasm.py pool  COPIES.tsv BAM MOLECULES.tsv OUT_DIR [--contig C]   -> OUT_DIR/pool/<fam>.tsv, pool.stats.tsv
  famasm.py lift  COPIES.tsv GENOME.fa OUT_DIR                        -> OUT_DIR/lift/<fam>.tsv, lift.unmapped.tsv
  famasm.py admit OUT_DIR ASSIGNMENTS.tsv COPIES.tsv                  -> OUT_DIR/admit/<fam>.tsv, admit.stats.tsv
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
        recs = []  # (read, idx, is_primary, AS, chain, ref_start, ref_end)
        for idx, (chrom, start, end, _) in sorted(members, key=lambda m: int(m[0])):
            for rec in bam.fetch(chrom, start, end):
                if rec.is_unmapped or rec.is_supplementary or rec.is_qcfail or not rec.has_tag("AS"):
                    continue
                # a record belongs to the copy whose locus it overlaps most (loci of one family do not overlap)
                ov = min(rec.reference_end, end) - max(rec.reference_start, start)
                if ov <= 0:
                    continue
                recs.append((rec.query_name, idx, not rec.is_secondary, int(rec.get_tag("AS")), fw_layer.chain_of(rec),
                             rec.reference_start, rec.reference_end))
        best = load_best_as(molecules_path, {r[0] for r in recs if not r[2]})
        kept, seen = [], set()
        secondary_kept = 0
        for read, idx, prim, as_, chain, rs, re_ in recs:
            if not prim and not good_secondary(as_, best.get(read, as_ + 1)):
                continue
            key = (read, idx)
            if key in seen:
                continue
            seen.add(key)
            kept.append((read, idx, prim, as_, chain, rs, re_))
            secondary_kept += not prim
        with open(os.path.join(out_dir, "pool", f"{fam}.tsv"), "w") as w:
            w.write("read\tcopy_idx\tis_primary\tAS\tchain\tref_start\tref_end\n")
            for read, idx, prim, as_, chain, rs, re_ in kept:
                w.write(f"{read}\t{idx}\t{int(prim)}\t{as_}\t{fw_layer.fmt(chain) if chain else ''}\t{rs}\t{re_}\n")
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
                read, idx, prim, _as, chain = line.rstrip("\n").split("\t")[:5]
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


def read_pool(path):
    """-> [(read, idx, is_primary, chain, ref_start, ref_end)]"""
    out = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            read, idx, prim, _as, chain, rs, re_ = line.rstrip("\n").split("\t")
            out.append((read, idx, int(prim), fw_layer.parse(chain) if chain else (), int(rs), int(re_)))
    return out


def read_liftrows(path):
    """-> [(read, idx, is_primary, fam_chain_str, private)]"""
    out = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            read, idx, prim, fc, priv = line.rstrip("\n").split("\t")
            out.append((read, idx, int(prim), fc, int(priv)))
    return out


def load_certified(path):
    """PREFIX.assignments.tsv -> {read: assigned_copy} for status == assigned."""
    out = {}
    with open(path) as fh:
        hdr = None
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = {k: i for i, k in enumerate(f)}
                continue
            if f[hdr["status"]] == "assigned":
                out[f[hdr["read_name"]]] = f[hdr["assigned_copy"]]
    return out


def family_chains(rows, min_reads=2):
    """lift rows -> set of family chains (tuples) seen in >= min_reads distinct molecules over all copies; private excluded."""
    mols = defaultdict(set)
    for read, _idx, _prim, fc, priv in rows:
        if priv or not fc:
            continue
        mols[fw_layer.parse(fc)].add(read)
    return {c for c, m in mols.items() if len(m) >= min_reads}


def _subchain_at(chain, sub):
    """index i with chain[i:i+len(sub)] == sub, else None."""
    if not sub or len(sub) > len(chain):
        return None
    try:
        i = chain.index(sub[0])
    except ValueError:
        return None
    return i if chain[i:i + len(sub)] == sub else None


def admit_at_copy(chain, rows, min_reads=2):
    """rows = (read, chain, span_start, span_end) of reads certified to this copy, all in the copy's genomic coordinates.
    'whole' when >= min_reads reads carry chain exactly; else 'compat' when every junction is covered by >= min_reads reads whose
    chain is a consecutive run of chain's junctions and whose span does not enter a neighbouring intron of chain; else None."""
    whole = {r for r, c, _s, _e in rows if c == chain}
    if len(whole) >= min_reads:
        return "whole"
    cover = [set() for _ in chain]
    for read, sub, s, e in rows:
        i = _subchain_at(chain, sub)
        if i is None:
            continue
        j = i + len(sub)
        if i > 0 and s < chain[i - 1][1]:
            continue
        if j < len(chain) and e > chain[j][0]:
            continue
        for k in range(i, j):
            cover[k].add(read)
    return "compat" if chain and all(len(c) >= min_reads for c in cover) else None


def admit_family(fam, pool_rows, lift_rows, certified, copies, reps, lifts, min_reads=2):
    """-> (rows, stats); row = (copy_idx, fam_chain_str, genomic_chain_str, via, certified_reads, span_start, span_end).
    certified = {read: copy_idx}. A family chain is tried at every copy with certified reads; a private chain only at its own."""
    stats = dict(family_chains=0, admitted_whole=0, admitted_compat=0, certified_without_record=0, uninvertible=0)
    rec_at = {(r[0], r[1]): r for r in pool_rows}            # (read, idx) -> pool row
    at_copy = defaultdict(list)                               # idx -> [(read, chain, s, e)] of certified reads with a record there
    fam_reads = {r[0] for r in pool_rows}
    for read in fam_reads:
        c = certified.get(read)
        if c is None:
            continue
        r = rec_at.get((read, c))
        if r is None:
            stats["certified_without_record"] += 1
            continue
        at_copy[c].append((read, r[3], r[4], r[5]))
    cands = []                                                # (idx, fam_chain_str, genomic_chain)
    fchains = family_chains(lift_rows, min_reads)
    stats["family_chains"] = len(fchains)
    for fc in sorted(fchains):
        for idx in sorted(at_copy, key=int):
            g = inverse_lift(lifts, copies, reps, fam, idx, fc)
            if g is None:
                stats["uninvertible"] += 1
                continue
            cands.append((idx, fw_layer.fmt(fc), g))
    priv = defaultdict(set)
    for read, idx, _p, fc, p in lift_rows:
        if p:
            priv[(idx, fc)].add(read)
    for (idx, fc), m in sorted(priv.items()):
        if len(m) >= min_reads and idx in at_copy:
            cands.append((idx, fc, fw_layer.parse(fc.split(":", 1)[1])))
    out = {}                                                  # (idx, genomic chain) -> [fam_chain_strs, via, reads, s, e]
    for idx, fc, g in cands:
        via = admit_at_copy(g, at_copy[idx], min_reads)
        if via is None:
            continue
        sup = [(r, s, e) for r, c, s, e in at_copy[idx] if _subchain_at(g, c) is not None]
        key = (idx, g)
        if key in out:
            out[key][0].append(fc)
            continue
        out[key] = [[fc], via, len({r for r, _s, _e in sup}), min(s for _r, s, _e in sup), max(e for _r, _s, e in sup)]
        stats["admitted_" + via] += 1
    rows = [(idx, ";".join(sorted(v[0])), fw_layer.fmt(g), v[1], v[2], v[3], v[4]) for (idx, g), v in sorted(out.items(), key=lambda kv: (int(kv[0][0]), kv[0][1]))]
    return rows, stats


def admit(out_dir, assignments_path, copies_path, min_reads=2):
    copies = fw_layer.read_copies(copies_path)
    reps = fw_layer.reps_of(copies)
    lifts = fw_layer.load_lift(os.path.join(out_dir, "map", "copies_vs_reps.paf"))
    certified = load_certified(assignments_path)
    os.makedirs(os.path.join(out_dir, "admit"), exist_ok=True)
    keys = ["family_chains", "admitted_whole", "admitted_compat", "certified_without_record", "uninvertible"]
    with open(os.path.join(out_dir, "admit.stats.tsv"), "w") as st:
        st.write("family\t" + "\t".join(keys) + "\n")
        for fn in sorted(os.listdir(os.path.join(out_dir, "lift"))):
            fam = fn[:-4]
            rows, stats = admit_family(fam, read_pool(os.path.join(out_dir, "pool", fn)),
                                       read_liftrows(os.path.join(out_dir, "lift", fn)), certified, copies, reps, lifts, min_reads)
            with open(os.path.join(out_dir, "admit", fn), "w") as w:
                w.write("copy_idx\tfam_chain\tgenomic_chain\tvia\tcertified_reads\tspan_start\tspan_end\n")
                for r in rows:
                    w.write("\t".join(map(str, r)) + "\n")
            st.write(fam + "\t" + "\t".join(str(stats[k]) for k in keys) + "\n")


def copy_strands(copies_path):
    """-> {(fam, idx): strand}"""
    out, hdr = {}, None
    for l in open(copies_path):
        f = l.rstrip("\n").split("\t")
        if hdr is None:
            hdr = {k: i for i, k in enumerate(f)}
            continue
        out[(f[hdr["family_id"]], f[hdr["copy_idx"]])] = f[hdr["strand"]] if "strand" in hdr else "."
    return out


def chain_exons(chain, span_start, span_end):
    """0-based half-open span + intron chain -> 1-based closed exons."""
    ex, s = [], span_start
    for x, y in chain:
        ex.append((s + 1, x))
        s = y
    ex.append((s + 1, span_end))
    return ex


def gtf(out_dir, base_gtf, copies_path, out_gtf):
    """FAMASM GTF = BASE verbatim + one transcript per admitted (copy, genomic chain) absent from BASE on that contig."""
    copies = fw_layer.read_copies(copies_path)
    strands = copy_strands(copies_path)
    base_chains = {(c, chain) for c, _s, _e, chain in fw_layer.gtf_tx(base_gtf).values()}
    adir = os.path.join(out_dir, "admit")
    appended = already = 0
    with open(out_gtf, "wb") as w:
        with open(base_gtf, "rb") as fh:
            for line in fh:
                w.write(line)
        for fn in sorted(os.listdir(adir)) if os.path.isdir(adir) else []:
            fam, n = fn[:-4], 0
            with open(os.path.join(adir, fn)) as fh:
                next(fh)
                for line in fh:
                    idx, _fc, g, via, nreads, s, e = line.rstrip("\n").split("\t")
                    chrom, a, b, _ = copies[(fam, idx)]
                    chain = fw_layer.parse(g)
                    if (chrom, chain) in base_chains:
                        already += 1
                        continue
                    s, e = max(int(s), a), min(int(e), b)
                    ex = chain_exons(chain, s, e)
                    if any(x > y for x, y in ex):
                        continue
                    n += 1
                    st = strands.get((fam, idx), ".")
                    attrs = (f'gene_id "FAMASM_{fam}_{idx}"; transcript_id "FAMASM_{fam}_{idx}_{n}"; reads "{nreads}"; famasm "1"; '
                             f'family "{fam}"; via "{via}"; certified_reads "{nreads}";')
                    w.write(f"{chrom}\tfamasm\ttranscript\t{ex[0][0]}\t{ex[-1][1]}\t.\t{st}\t.\t{attrs}\n".encode())
                    for i, (x, y) in enumerate(ex, 1):
                        w.write(f"{chrom}\tfamasm\texon\t{x}\t{y}\t.\t{st}\t.\t{attrs} exon_number \"{i}\";\n".encode())
                    appended += 1
    with open(os.path.join(out_dir, "gtf.stats.tsv"), "w") as w:
        w.write(f"appended\talready_in_base\n{appended}\t{already}\n")


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
    elif cmd == "admit":
        admit(*a[:3])
    elif cmd == "gtf":
        gtf(*a[:4])
    else:
        sys.exit(f"unknown command {cmd}")


if __name__ == "__main__":
    main()
