#!/usr/bin/env python3
"""
tau-SWEEP: precision/recall of PSV copy-assignment vs the decisive-margin tau, so the operating
point can be chosen from data (not hard-coded). Uses the SAME engine as production
(bench/copy_assign.py::assign_read, which src/rustle/vg_family/copy_assign.rs mirrors).

Gate under test (the reframed lever): a read is ASSIGNED iff n_decisive>=1 (identifiability) AND
log-LR margin >= tau. tau = the one calibrated knob replacing (min_psv=3, vote_margin=1).
  - tau -> 0+   : argmax / max recall (assign every resolvable read)
  - tau -> inf  : the production min_psv=3 conservative corner (discard the tail)

Substrate 1 (RIGOROUS, true labels): sim5x K-ladder. Read name K{K}_c{c}_r{r} => true copy c.
Substrate 2 (real, silver standard): GGO families (added in run_ggo if data present).

For each read we assign ONCE (recording true_c, best_c, margin, n_dec) then threshold many tau.
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pysam
import copy_assign as CA


def collect_sim(K):
    """Return list of (true_c, best_c, margin, n_dec) for every primary read of sim5x_K{K}."""
    ref = f"{CA.SIM5X}/sim5x_K{K}.ref.fa"
    bam = f"{CA.SIM5X}/sim5x_K{K}.bam"
    if not os.path.exists(bam):
        return None
    fa = pysam.FastaFile(ref)
    copies, unit = CA.sim5x_copies(fa, K)
    seqs = CA.read_copy_seqs(fa, copies)
    off2gen = {c[0]: (lambda cp: (lambda off: CA.seqoff_to_genomic(cp, off)))(c) for c in copies}
    cols = CA.discover_psvs(copies, seqs, off2gen)
    copy_vecs = {c[0]: {} for c in copies}
    copy_gpos = {c[0]: {} for c in copies}
    for j, rec in enumerate(cols):
        for name, v in rec.items():
            if v is not None:
                gpos, base = v
                copy_vecs[name][j] = base
                copy_gpos[name][j] = gpos
    b = pysam.AlignmentFile(bam, "rb")
    out = []
    for aln in b.fetch():
        if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
            continue
        true_c = int(aln.query_name.split("_c")[1].split("_")[0])
        mapped_c = aln.reference_start // unit
        mc = f"copy{mapped_c}"
        want = set(copy_gpos[mc].values())
        base_at, _ = CA.read_bases_and_introns(aln, want)
        gpos2col = {g: j for j, g in copy_gpos[mc].items()}
        obs = {gpos2col[g]: base_at[g] for g in base_at if g in gpos2col}
        best, margin, n_span, n_dec, _ = CA.assign_read(obs, copy_vecs)
        best_c = int(best.replace("copy", ""))
        out.append((true_c, best_c, float(margin), int(n_dec)))
    return out, len(cols)


def pr_at_tau(rows, tau):
    """precision = correct among assigned ; recall = assigned / all reads. Assigned iff n_dec>=1 & margin>=tau."""
    n = len(rows)
    n_assigned = correct = n_resolvable = 0
    for true_c, best_c, margin, n_dec in rows:
        if n_dec >= 1:
            n_resolvable += 1
            if margin >= tau:
                n_assigned += 1
                correct += (best_c == true_c)
    prec = correct / n_assigned if n_assigned else float("nan")
    rec = n_assigned / n if n else float("nan")
    return dict(tau=tau, n=n, n_resolvable=n_resolvable, n_assigned=n_assigned,
                correct=correct, precision=prec, recall=rec,
                recall_of_resolvable=(n_assigned / n_resolvable if n_resolvable else float("nan")))


TAUS = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.9, 9.0, 12.0, 16.0, 25.0]
# reference points: tau=2.0 current default; tau=6.9 ~ ln((1-1e-3)/1e-3) the Eichler-faithful p=1e-3


def collect_ggo(limit=120):
    """Per-read tuples over real GGO co-located families, recording the PSV+junction margin ONCE so
    tau can be swept. 'silver' = the read's best-overlap copy (mc); agreement = PSV-assigned best==mc
    among MQ>0 reads (the same silver standard assign_family uses). Returns list of
    (margin_j, agree, mq_pos, nd_p_ge1, nd_j_ge1)."""
    meta = CA.load_meta()
    skel = CA.load_skel()
    spliced = CA.load_spliced()
    rescued_fam = CA.load_rescued(meta, skel, spliced)
    fams = CA.colocated_families(meta, rescued_fam=rescued_fam)
    bam = pysam.AlignmentFile(CA.GGO_BAM, "rb")
    rows = []
    nfam = 0
    for _fid, _dom, copies in fams:
        if nfam >= limit:
            break
        # replicate assign_family setup, but record per-read margins
        seqs, idx2gen, gen2off, strand, cbounds, ok = {}, {}, {}, {}, {}, []
        for c in copies:
            tid = c[0]
            s = spliced.get(tid)
            g = CA.copy_idx2gen(tid, meta, skel)
            if s and len(s) == len(g):
                seqs[tid] = s; idx2gen[tid] = g
                gen2off[tid] = {gp: i for i, gp in enumerate(g)}
                strand[tid] = c[4]; cbounds[tid] = CA.copy_boundaries(tid, meta, skel)
                ok.append(c)
        if len(ok) < 2:
            continue
        copies = ok
        off2gen = {c[0]: (lambda nm: (lambda off: idx2gen[nm][off]))(c[0]) for c in copies}
        cols = CA.discover_psvs(copies, seqs, off2gen)
        copy_vecs = {c[0]: {} for c in copies}; copy_gpos = {c[0]: {} for c in copies}
        for j, rec in enumerate(cols):
            for name, v in rec.items():
                if v is not None:
                    copy_vecs[name][j] = v[1]; copy_gpos[name][j] = v[0]
        if not any(copy_vecs.values()):
            continue
        spans = {c[0]: (c[1], c[2] - 1, c[3]) for c in copies}
        contig = copies[0][1]
        lo = min(c[2] - 1 for c in copies); hi = max(c[3] for c in copies)
        nfam += 1
        for aln in bam.fetch(contig, lo, hi):
            if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
                continue
            mc, bovl = None, 0
            for name, (ct, s0, e0) in spans.items():
                ov = min(aln.reference_end, e0) - max(aln.reference_start, s0)
                if ov > bovl:
                    bovl, mc = ov, name
            if mc is None:
                continue
            want = set(copy_gpos[mc].values())
            base_at, introns = CA.read_bases_and_introns(aln, want)
            g2c = {g: j for j, g in copy_gpos[mc].items()}
            obs = {}
            for g, bb in base_at.items():
                if g in g2c and bb is not None:
                    obs[g2c[g]] = CA.rc_base(bb) if strand[mc] == "-" else bb
            g2o = gen2off[mc]; rbounds = set()
            for d0, a0 in introns:
                o1, o2 = g2o.get(d0 - 1), g2o.get(a0)
                if o1 is not None and o2 is not None:
                    rbounds.add(max(o1, o2))
            _, m_p, _, nd_p, _ = CA.assign_read(obs, copy_vecs)
            best, m_j, _, nd_j, _ = CA.assign_read(obs, copy_vecs, read_bounds=rbounds, copy_bounds=cbounds)
            rows.append((float(m_j), int(best == mc), int(aln.mapping_quality > 0),
                         int(nd_p >= 1), int(nd_j >= 1)))
    return rows, nfam


def ggo_pr_at_tau(rows, tau):
    """recall = assigned / all ; silver_agreement = (best==mc) among ASSIGNED & MQ>0."""
    n = len(rows)
    n_assigned = n_uniq = n_uniq_agree = n_resolvable = 0
    for m_j, agree, mq_pos, nd_p, nd_j in rows:
        if nd_j >= 1:
            n_resolvable += 1
            if m_j >= tau:
                n_assigned += 1
                if mq_pos:
                    n_uniq += 1
                    n_uniq_agree += agree
    rec = n_assigned / n if n else float("nan")
    agr = n_uniq_agree / n_uniq if n_uniq else float("nan")
    return dict(tau=tau, n=n, n_resolvable=n_resolvable, n_assigned=n_assigned,
                recall=rec, n_uniq=n_uniq, silver_agreement=agr)


def main():
    result = {"taus": TAUS, "sim": {}}
    print("=== sim5x K-ladder: TRUE precision/recall vs tau (engine = copy_assign.assign_read) ===")
    print("tau ref: 2.0=current default, 6.9=p_target 1e-3 (Eichler-faithful)\n")
    for K in (1, 2, 4, 8):
        c = collect_sim(K)
        if c is None:
            continue
        rows, ncols = c
        curve = [pr_at_tau(rows, t) for t in TAUS]
        result["sim"][K] = dict(n=len(rows), psv_cols=ncols, curve=curve)
        print(f"-- K={K}  ({len(rows)} reads, {ncols} PSV columns) --")
        print(f"   {'tau':>5} {'recall':>7} {'precision':>9} {'assigned':>8} {'correct':>7}")
        for r in curve:
            print(f"   {r['tau']:>5.1f} {r['recall']:>7.3f} {r['precision']:>9.4f} "
                  f"{r['n_assigned']:>8d} {r['correct']:>7d}")
        print()
    # --- real GGO families (silver-standard agreement vs tau) ---
    try:
        limit = int(os.environ.get("GGO_FAM_LIMIT", "120"))
        print(f"=== real GGO co-located families: silver-standard agreement vs tau (<= {limit} families) ===")
        grows, nfam = collect_ggo(limit=limit)
        gcurve = [ggo_pr_at_tau(grows, t) for t in TAUS]
        result["ggo"] = dict(n_reads=len(grows), n_families=nfam, curve=gcurve)
        print(f"   {len(grows)} reads over {nfam} families")
        print(f"   {'tau':>5} {'recall':>7} {'assigned':>8} {'uniqMQ>0':>8} {'silver_agree':>12}")
        for r in gcurve:
            print(f"   {r['tau']:>5.1f} {r['recall']:>7.3f} {r['n_assigned']:>8d} "
                  f"{r['n_uniq']:>8d} {r['silver_agreement']:>12.4f}")
        print()
    except Exception as e:
        print(f"   [GGO real sweep skipped: {e}]")

    json.dump(result, open("bench/psv_tau_sweep.json", "w"), indent=1)
    print("wrote bench/psv_tau_sweep.json")
    return result


if __name__ == "__main__":
    main()
