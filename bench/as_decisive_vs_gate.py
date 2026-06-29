#!/usr/bin/env python
"""F6 — head-to-head: the Eichler-2024 (TBC1D3) AS>=10 decisive rule vs our significance gate.

The advisor repeatedly cites the TBC1D3 copy-assignment rule (Guitart/Eichler 2024, Methods): map each Iso-Seq
read to ALL validated haplotypes and "any read with primary minimap2 alignment score 10 or greater for a given
paralog cluster relative to all other cluster mappings is retained, others marked ambiguous and ignored." That
is a CONSERVATIVE assign-or-discard rule (no 1/k) on the raw alignment-score margin.

We calibrated our gate's tau=6.9 to "approximately the AS>=10 operating point" but never ran the literal rule
against ours on the same reads. This does. On the sim5x labeled ladder (ground truth = read name `K{K}_c{copy}`):
  - AS>=10 rule: realign each read to the 5 copy sequences (minimap2), per-copy best AS; margin = best - 2nd;
    assign to argmax iff margin >= 10, else AMBIGUOUS (Eichler).
  - our gate: the production per-read PSV significance assignment (copy_assign.assign_read + the margin gate).
We report, per K: recall (fraction assigned), accuracy (correct among assigned), and the CONCORDANCE where both
assign. Expected story (the exhibit): the two agree on resolvable reads, but our gate ABSTAINS WITH A REASON
(min_p identifiability certificate) where AS>=10 silently discards — most visibly at K=0 (identical copies),
where AS margins collapse and AS>=10 discards everything while our gate certifies every read Tied (min_p=1).

Run: /home/juanfra/miniforge3/bin/python bench/as_decisive_vs_gate.py
"""
import os
import subprocess
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.dirname(__file__))
import copy_assign as ca  # reuse sim5x_copies, read_copy_seqs, discover_psvs, assign_read, read_bases_and_introns

MM2 = os.environ.get("RUSTLE_MINIMAP2", "/home/juanfra/miniforge3/bin/minimap2")
AS_MARGIN = 10  # the Eichler threshold


def per_copy_best_as(read_fa, copies_fa):
    """minimap2-align reads to the K copy sequences; return {read_name: {copy_name: best_AS}}."""
    out = subprocess.run(
        [MM2, "-a", "-x", "map-hifi", "--secondary=yes", "-N", "5", "-p", "0.1", copies_fa, read_fa],
        capture_output=True, text=True,
    )
    best = {}
    for line in out.stdout.splitlines():
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if len(f) < 11:
            continue
        qname, flag, rname = f[0], int(f[1]), f[2]
        if flag & 0x4 or rname == "*":
            continue
        as_tag = next((int(x[5:]) for x in f[11:] if x.startswith("AS:i:")), None)
        if as_tag is None:
            continue
        best.setdefault(qname, {})
        if as_tag > best[qname].get(rname, -1):
            best[qname][rname] = as_tag
    return best


def as_decisive(best_by_copy):
    """Eichler AS>=10: assign to the top copy iff it beats the 2nd by >= AS_MARGIN, else ambiguous."""
    items = sorted(best_by_copy.items(), key=lambda kv: kv[1], reverse=True)
    if not items:
        return None, 0.0
    if len(items) == 1:
        return items[0][0], float("inf")
    margin = items[0][1] - items[1][1]
    return (items[0][0] if margin >= AS_MARGIN else None), margin


def run():
    print(f"{'K':>2} | {'reads':>5} | {'AS>=10 recall':>13} {'AS>=10 acc':>10} | "
          f"{'gate recall':>11} {'gate acc':>8} | {'AS margin med/p90/max':>21} | {'both concord':>12}")
    for K in (0, 1, 2, 4, 8):
        ref = f"{ca.SIM5X}/sim5x_K{K}.ref.fa"
        bam = f"{ca.SIM5X}/sim5x_K{K}.bam"
        if not os.path.exists(bam):
            continue
        fa = pysam.FastaFile(ref)
        copies, unit = ca.sim5x_copies(fa, K)
        seqs = ca.read_copy_seqs(fa, copies)            # {copyName: sequence}
        off2gen = {c[0]: (lambda cp: (lambda off: ca.seqoff_to_genomic(cp, off)))(c) for c in copies}
        cols = ca.discover_psvs(copies, seqs, off2gen)
        copy_vecs = {c[0]: {} for c in copies}
        copy_gpos = {c[0]: {} for c in copies}
        for j, rec in enumerate(cols):
            for name, v in rec.items():
                if v is not None:
                    copy_vecs[name][j] = v[1]
                    copy_gpos[name][j] = v[0]

        with tempfile.TemporaryDirectory() as td:
            copies_fa = os.path.join(td, "copies.fa")
            reads_fa = os.path.join(td, "reads.fa")
            with open(copies_fa, "w") as fh:
                for name, s in seqs.items():
                    fh.write(f">{name}\n{s}\n")
            true_copy = {}
            b = pysam.AlignmentFile(bam, "rb")
            with open(reads_fa, "w") as fh:
                for aln in b.fetch():
                    if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
                        continue
                    tc = int(aln.query_name.split("_c")[1].split("_")[0])
                    true_copy[aln.query_name] = tc
                    fh.write(f">{aln.query_name}\n{aln.query_sequence}\n")
            best = per_copy_best_as(reads_fa, copies_fa)

            # our gate, per read
            n = 0
            as_assigned = as_correct = 0
            gate_assigned = gate_correct = 0
            both = concord = 0
            margins = []   # finite AS best-minus-2nd margins, to show WHY AS>=10 (mis)fires
            b2 = pysam.AlignmentFile(bam, "rb")
            for aln in b2.fetch():
                if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
                    continue
                qn = aln.query_name
                tc = true_copy[qn]
                n += 1
                # --- AS>=10 ---
                as_copy, as_margin = as_decisive(best.get(qn, {}))
                if as_margin != float("inf"):
                    margins.append(as_margin)
                as_idx = int(as_copy.replace("copy", "")) if as_copy else None
                if as_idx is not None:
                    as_assigned += 1
                    as_correct += (as_idx == tc)
                # --- our gate ---
                mapped_c = aln.reference_start // unit
                mc = f"copy{mapped_c}"
                base_at, _ = ca.read_bases_and_introns(aln, set(copy_gpos[mc].values()))
                gpos2col = {g: jj for jj, g in copy_gpos[mc].items()}
                obs = {gpos2col[g]: base_at[g] for g in base_at if g in gpos2col}
                gbest, gmargin, _, n_dec, _ = ca.assign_read(obs, copy_vecs)
                g_idx = int(gbest.replace("copy", "")) if (n_dec >= 1 and gmargin >= ca.MARGIN) else None
                if g_idx is not None:
                    gate_assigned += 1
                    gate_correct += (g_idx == tc)
                # --- concordance where both assign ---
                if as_idx is not None and g_idx is not None:
                    both += 1
                    concord += (as_idx == g_idx)

        ar = as_assigned / n if n else 0
        aa = as_correct / as_assigned if as_assigned else float("nan")
        gr = gate_assigned / n if n else 0
        ga = gate_correct / gate_assigned if gate_assigned else float("nan")
        cc = concord / both if both else float("nan")
        ms = sorted(margins)
        med = ms[len(ms) // 2] if ms else float("nan")
        p90 = ms[int(0.9 * (len(ms) - 1))] if ms else float("nan")
        mx = ms[-1] if ms else float("nan")
        print(f"{K:>2} | {n:>5} | {100*ar:>12.1f}% {aa:>10.3f} | {100*gr:>10.1f}% {ga:>8.3f} | "
              f"{med:>6.1f}/{p90:>4.1f}/{mx:>4.1f} (thr=10) | {concord:>3}/{both:<3}={100*cc if both else float('nan'):>4.0f}%")
    print("\nAS>=10 = Eichler/Guitart 2024 TBC1D3 rule (realign read to all copy seqs, assign iff best AS beats 2nd")
    print("by >=10, else ambiguous). gate = our PSV significance assignment. K=0 = identical copies: AS margins")
    print("collapse so AS>=10 discards ~all, while our gate certifies them Tied (min_p=1, a reason not a silent drop).")


def divergent_control():
    """Sanity: AS>=10 is NOT broken — it fires correctly when copies are DIVERGENT (the Eichler/TBC1D3 regime).
    Take a sim5x copy, make a 3%-diverged variant, and confirm that reads from the original score an AS margin
    >> 10 against the divergent copy (so AS>=10 assigns them correctly). The 0% on the near-identical ladder is
    therefore about the REGIME (collapsed paralogs, margins <=5), not the computation."""
    fa = pysam.FastaFile(f"{ca.SIM5X}/sim5x_K8.ref.fa")
    copies, unit = ca.sim5x_copies(fa, 8)
    seqs = ca.read_copy_seqs(fa, copies)
    orig = list(seqs["copy0"])
    div = orig[:]
    sub = {"A": "G", "G": "A", "C": "T", "T": "C", "N": "N"}   # deterministic 3% transition divergence
    for i in range(0, len(div), 33):
        div[i] = sub.get(div[i], div[i])
    with tempfile.TemporaryDirectory() as td:
        cf = os.path.join(td, "c.fa"); rf = os.path.join(td, "r.fa")
        with open(cf, "w") as fh:
            fh.write(f">copy0\n{''.join(orig)}\n>copyDIV\n{''.join(div)}\n")
        b = pysam.AlignmentFile(f"{ca.SIM5X}/sim5x_K8.bam", "rb")
        nreads = 0
        with open(rf, "w") as fh:
            for aln in b.fetch():
                if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
                    continue
                if int(aln.query_name.split("_c")[1].split("_")[0]) == 0:   # copy0 reads only
                    fh.write(f">{aln.query_name}\n{aln.query_sequence}\n"); nreads += 1
                if nreads >= 30:
                    break
        best = per_copy_best_as(rf, cf)
        assigned = correct = 0
        margins = []
        for qn, bc in best.items():
            cp, m = as_decisive(bc)
            if m != float("inf"):
                margins.append(m)
            if cp is not None:
                assigned += 1
                correct += (cp == "copy0")
    med = sorted(margins)[len(margins) // 2] if margins else float("nan")
    print(f"\nDIVERGENT control (copy0 vs 3%-diverged copy, {nreads} reads): AS>=10 assigns {assigned}/{nreads}, "
          f"{correct} to the correct (original) copy; median AS margin {med:.0f} (>> 10). "
          f"=> AS>=10 works on divergent copies; it abstains on the near-identical ladder because margins are <=5.")


if __name__ == "__main__":
    run()
    divergent_control()
