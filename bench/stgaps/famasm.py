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
    """Two passes: fetch every >= 2-copy family's records from the BAM, then ONE scan of molecules.tsv for the secondaries' best AS."""
    import pysam
    copies = fw_layer.read_copies(copies_path)
    if contig:
        copies = {k: v for k, v in copies.items() if v[0] == contig}
    fams = defaultdict(list)
    for (fam, idx), v in copies.items():
        fams[fam].append((idx, v))
    os.makedirs(os.path.join(out_dir, "pool"), exist_ok=True)
    bam = pysam.AlignmentFile(bam_path)
    recs_by_fam, sec_names = {}, set()
    for fam, members in sorted(fams.items()):
        if len(members) < 2:
            continue
        recs = []  # (read, idx, is_primary, AS, chain, ref_start, ref_end)
        for idx, (chrom, start, end, _) in sorted(members, key=lambda m: int(m[0])):
            for rec in bam.fetch(chrom, start, end):
                if rec.is_unmapped or rec.is_supplementary or rec.is_qcfail or not rec.has_tag("AS"):
                    continue
                # a record belongs to the copy whose locus it overlaps most (loci of one family do not overlap)
                if min(rec.reference_end, end) - max(rec.reference_start, start) <= 0:
                    continue
                recs.append((rec.query_name, idx, not rec.is_secondary, int(rec.get_tag("AS")), fw_layer.chain_of(rec),
                             rec.reference_start, rec.reference_end))
                if rec.is_secondary:
                    sec_names.add(rec.query_name)
        recs_by_fam[fam] = recs
    print(f"[pool] {len(recs_by_fam)} families, {sum(map(len, recs_by_fam.values()))} records, "
          f"{len(sec_names)} molecules with a secondary record; scanning {molecules_path}", file=sys.stderr, flush=True)
    best = load_best_as(molecules_path, sec_names)
    stats = []
    for fam, recs in recs_by_fam.items():
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
    """PREFIX.assign.assignments.tsv -> {(family_id, read): copy_idx | None}; None = a row exists but is not assigned. One row per read x family: a read is
    certified in ONE family and tied in the others, so the key carries the family. copy_assign sorts copies and reports its own
    index in assigned_copy; the catalog's copy_idx is catalog_copy_idx (used when the column exists and is not NA)."""
    out = {}
    with open(path) as fh:
        hdr = None
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = {k: i for i, k in enumerate(f)}
                ci = hdr.get("catalog_copy_idx")
                continue
            key = (f[hdr["family_id"]], f[hdr["read_name"]])
            if f[hdr["status"]] != "assigned":
                out.setdefault(key, None)      # contested (tied / ambiguous): never certified by placement either
                continue
            out[key] = f[ci] if ci is not None and f[ci] != "NA" else f[hdr["assigned_copy"]]
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
    certified = {(family_id, read): copy_idx} (load_certified). A family chain is tried at every copy with certified reads; a private chain only at its own."""
    stats = dict(family_chains=0, admitted_whole=0, admitted_compat=0, certified_without_record=0, uninvertible=0,
                 certified_by_placement=0)
    rec_at = {(r[0], r[1]): r for r in pool_rows}            # (read, idx) -> pool row
    at_copy = defaultdict(list)                               # idx -> [(read, chain, s, e)] of certified reads with a record there
    recs_of = defaultdict(list)
    for r in pool_rows:
        recs_of[r[0]].append(r)
    for read, recs in recs_of.items():
        c = certified.get((fam, read))
        if (fam, read) not in certified and len(recs) == 1 and recs[0][2]:
            # amendment 1: no assignment row and a single pool record (the primary, no secondary within 2% of the best AS):
            # the aligner placed it without a rival, so the placement itself certifies the copy
            c = recs[0][1]
            stats["certified_by_placement"] += 1
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
    keys = ["family_chains", "admitted_whole", "admitted_compat", "certified_without_record", "uninvertible", "certified_by_placement"]
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
    appended = already = outside = 0
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
                        outside += 1
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
        w.write(f"appended\talready_in_base\tspan_outside_locus\n{appended}\t{already}\t{outside}\n")


def _strata(path):
    """strata.txt -> {(arm, stratum): (recall%, precision%)}"""
    out = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) == 8 and f[0] != "arm":
            out[(f[0], f[1])] = (float(f[4]), float(f[7]))
    return out


def _base_and_run(path, cols, family=None):
    """mainstay/pair table -> (BASE, run) tuples of the named columns: the first and the last data row, restricted to
    `family == family` rows when given (mainstay_score prints one row per family pattern then ALL, per run)."""
    rows, hdr = [], None
    for l in open(path):
        if l[0] == "#":
            continue
        f = l.rstrip("\n").split("\t")
        if hdr is None:
            hdr = {k: i for i, k in enumerate(f)}
            continue
        if family is None or f[hdr["family"]] == family:
            rows.append(f)
    return [tuple(float(r[hdr[c]] or 0) for c in cols) for r in (rows[0], rows[-1])]


def score(run_dir):
    """M1-M4 verdict from strata.txt, mainstays.txt, pairs.txt; BASE lines kept and single-copy identity from the GTFs."""
    env = dict(l.rstrip("\n").split("=", 1) for l in open(os.path.join(run_dir, "famasm.env")) if "=" in l)
    j = lambda *p: os.path.join(run_dir, *p)
    st = _strata(j("strata.txt"))
    br, bp = st[("base", "multi")]
    fr, fp = st[("famasm", "multi")]
    sbr, sbp = st[("base", "single")]
    sfr, sfp = st[("famasm", "single")]
    base = open(env["BASE_GTF"], "rb").read()
    fam = open(j("famasm.gtf"), "rb").read()
    kept = fam.startswith(base)
    copies = fw_layer.read_copies(env["COPIES"])
    nfam = defaultdict(int)
    for (f, _i) in copies:
        nfam[f] += 1
    multi_loci = defaultdict(list)
    for (f, _i), (c, a, b, _) in copies.items():
        if nfam[f] >= 2:
            multi_loci[c].append((a, b))
    single_ok = kept
    for line in fam[len(base):].decode().splitlines() if kept else []:
        f = line.split("\t")
        if len(f) < 9 or 'famasm "1"' not in f[8]:
            single_ok = False
            break
        x, y = int(f[3]), int(f[4])
        if not any(a <= x and y <= b for a, b in multi_loci.get(f[0], ())):
            single_ok = False
            break
    (m_b, m_r) = _base_and_run(j("mainstays.txt"), ["recalled", "unrel_foreign"], family="ALL")
    (p_b, p_r) = _base_and_run(j("pairs.txt"), ["precision", "recall"])
    bars = [
        ("M1", f"multi recall {br:.1f} -> {fr:.1f} (need +1.0)", fr - br >= 1.0 - 1e-9),
        ("M2", f"multi precision {bp:.1f} -> {fp:.1f} (need >= base-0.5)", fp >= bp - 0.5 - 1e-9),
        ("M3", f"single recall {sbr:.1f} -> {sfr:.1f}, precision {sbp:.1f} -> {sfp:.1f} (|d|<=0.1), base_lines_kept {'yes' if kept else 'no'}",
         abs(sfr - sbr) <= 0.1 + 1e-9 and abs(sfp - sbp) <= 0.1 + 1e-9 and kept),
        ("M4", f"mainstays recalled {m_b[0]:.0f} -> {m_r[0]:.0f}, unrel_foreign {m_b[1]:.0f} -> {m_r[1]:.0f}; pairs precision {p_b[0]:.3f} -> {p_r[0]:.3f}, recall {p_b[1]:.3f} -> {p_r[1]:.3f} (>= base-0.01)",
         m_r[0] >= m_b[0] and m_r[1] <= m_b[1] and p_r[0] >= p_b[0] - 0.01 - 1e-9 and p_r[1] >= p_b[1] - 0.01 - 1e-9),
    ]
    lines = [f"FAMASM {env.get('SAMPLE', '?')} ({env.get('SPECIES', '?')})"]
    for name, text, ok in bars:
        lines.append(f"{name} {'PASS' if ok else 'FAIL'}: {text}")
    lines.append(f"single_copy_identical {'yes' if single_ok else 'no'}")
    for fn in ("gtf.stats.tsv", "pool.stats.tsv", "admit.stats.tsv", "lift.unmapped.tsv"):
        if os.path.exists(j(fn)):
            rows = [l.rstrip("\n").split("\t") for l in open(j(fn))]
            if fn == "gtf.stats.tsv":
                lines.append("gtf: " + ", ".join(f"{k}={v}" for k, v in zip(rows[0], rows[1])))
            elif fn == "lift.unmapped.tsv":
                lines.append(f"copies without a lift table: {len(rows) - 1}")
            else:
                tot = [sum(int(r[i]) for r in rows[1:]) for i in range(1, len(rows[0]))]
                lines.append(f"{fn[:-10]}: families={len(rows) - 1}, " + ", ".join(f"{k}={v}" for k, v in zip(rows[0][1:], tot)))
    lines.append("verdict " + ("PASS" if all(ok for _n, _t, ok in bars) else "FAIL"))
    return "\n".join(lines) + "\n"


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
    elif cmd == "score":
        sys.stdout.write(score(a[0]))
    else:
        sys.exit(f"unknown command {cmd}")


if __name__ == "__main__":
    main()
