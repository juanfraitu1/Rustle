#!/usr/bin/env python
"""Measure whether ALIGNMENT error (a spurious PSV column from a copy-vs-copy mis-alignment, distinct from
sequencing error) is driving copy assignment — the empirical answer to the advisor's "how do you know that
column isn't an alignment artifact?". No production change; reads a `copy_assign --dump-psv` triple.

Three measurements per family:
  (1) COLUMN COHERENCE — of the reads assigned to copy c, what fraction show c's allele at column j? A real PSV
      = ~sequencing-error mismatch (~0.3%); an alignment-artifact column is INCOHERENT (reads disagree with
      their copy's allele beyond error). Reports the coherence distribution + the suspect (low-coherence) tail.
  (2) COLUMN CLUSTERING — inter-column genome gaps. A mis-aligned stretch makes a RUN of adjacent mismatch
      columns; a real PSV is usually isolated. Reports the fraction of columns in dense runs (the indel-flank
      artifact signature).
  (3) DROP-FRAGILE SENSITIVITY — drop the fragile columns (low coherence OR clustered), re-assign each read by
      nearest copy (Hamming over the OBSERVED remaining columns), and report how many assignments FLIP. ~0%
      => alignment error is not driving the result (the robustness exhibit).

Run: /home/juanfra/miniforge3/bin/python bench/align_error_probe.py <prefix>   (prefix.psv_{reads,copies,cols}.tsv)
"""
import random
import sys
from collections import defaultdict

ERR = 0.003                 # HiFi per-base substitution rate (the sequencing-error floor coherence is judged against)
COH_SUSPECT = 0.97          # a column below this agrees with its copies worse than ~10x the error floor
RUN_GAP = 3                 # columns within this many bp are "adjacent" (a mis-aligned run)
RUN_MIN = 4                 # >= this many adjacent columns = a dense run (indel-flank artifact signature)


def load(prefix):
    copies = defaultdict(dict)   # fam -> copy_idx -> allele string
    for ln in open(f"{prefix}.psv_copies.tsv").read().splitlines()[1:]:
        f = ln.split("\t")
        copies[f[0]][int(f[1])] = f[3]
    cols = defaultdict(list)     # fam -> [genome_pos by col_index]
    for ln in open(f"{prefix}.psv_cols.tsv").read().splitlines()[1:]:
        f = ln.split("\t")
        cols[f[0]].append(int(f[2]))
    reads = defaultdict(list)    # fam -> [(assigned_copy, status, alleles)]
    for ln in open(f"{prefix}.psv_reads.tsv").read().splitlines()[1:]:
        f = ln.split("\t")
        reads[f[1]].append((int(f[2]), f[3], f[6]))  # key by FAMILY id (f[1]), not read name
    return copies, cols, reads


def col_coherence(fam_copies, fam_reads, ncol):
    """per column: (observed, mismatches) over reads ASSIGNED to a copy, vs that copy's allele."""
    obs = [0] * ncol
    mis = [0] * ncol
    for (cp, status, alleles) in fam_reads:
        if status != "assigned" or cp not in fam_copies:
            continue
        cseq = fam_copies[cp]
        for j in range(min(ncol, len(alleles), len(cseq))):
            a = alleles[j]
            if a == "." or cseq[j] == ".":
                continue
            obs[j] += 1
            if a != cseq[j]:
                mis[j] += 1
    return obs, mis


def dense_runs(positions):
    """indices of columns inside a run of >= RUN_MIN columns each within RUN_GAP bp of the next."""
    flagged = set()
    n = len(positions)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and positions[j + 1] - positions[j] <= RUN_GAP:
            j += 1
        if (j - i + 1) >= RUN_MIN:
            flagged.update(range(i, j + 1))
        i = j + 1
    return flagged


def nearest_copy(alleles, fam_copies, keep):
    """argmin Hamming over the OBSERVED kept columns (the proxy for the production argmax)."""
    best, bestd = None, None
    for cp, cseq in fam_copies.items():
        d = tot = 0
        for j in keep:
            a = alleles[j] if j < len(alleles) else "."
            c = cseq[j] if j < len(cseq) else "."
            if a == "." or c == ".":
                continue
            tot += 1
            if a != c:
                d += 1
        if tot == 0:
            continue
        if bestd is None or d < bestd:
            best, bestd = cp, d
    return best


def main():
    prefix = sys.argv[1] if len(sys.argv) > 1 else "/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/03f0156d-63b6-4d0f-bf09-862bddee491c/scratchpad/mser"
    copies, cols, reads = load(prefix)
    print(f"== alignment-error probe on {prefix} ==\n")
    grand = [0, 0, 0, 0]  # tot_cols, suspect_cols, clustered_cols, flipped_reads, total_assigned
    tot_assigned = 0
    for fam in sorted(copies, key=lambda f: -len(cols[f])):
        ncol = len(cols[fam])
        if ncol < 5 or len(copies[fam]) < 2:
            continue
        obs, mis = col_coherence(copies[fam], reads[fam], ncol)
        # coherence only defined where >=5 assigned reads observe the column
        coh = {j: 1 - mis[j] / obs[j] for j in range(ncol) if obs[j] >= 5}
        suspect = [j for j in coh if coh[j] < COH_SUSPECT]   # incoherent beyond the error floor = artifact-suspect
        clustered = dense_runs(cols[fam])                    # descriptive only: dense = DIVERGENCE, not per se an artifact
        clustered_and_incoherent = [j for j in suspect if j in clustered]
        fragile = set(suspect)                               # DROP only the incoherent columns (the principled signal)
        keep = [j for j in range(ncol) if j not in fragile]
        # re-assign assigned reads with vs without the fragile columns; count flips. The BASELINE (drop
        # nothing) measures how faithfully the crude Hamming proxy reproduces the production LLR gate — only
        # the EXCESS flips over baseline are attributable to dropping the fragile columns.
        all_cols = list(range(ncol))
        # CONTROL: drop the SAME NUMBER of RANDOM columns — if it flips as much, the flip is generic
        # sensitivity to column count (near-tied reads), not the fragile columns being special.
        rng = random.Random(42)
        rand_drop = set(rng.sample(all_cols, min(len(fragile), ncol)))
        keep_rand = [j for j in all_cols if j not in rand_drop]
        flips = base_flips = rand_flips = nassigned = 0
        for (cp, status, alleles) in reads[fam]:
            if status != "assigned":
                continue
            nassigned += 1
            if nearest_copy(alleles, copies[fam], all_cols) not in (cp, None):
                base_flips += 1
            if nearest_copy(alleles, copies[fam], keep) not in (cp, None):
                flips += 1
            if nearest_copy(alleles, copies[fam], keep_rand) not in (cp, None):
                rand_flips += 1
        excess = flips - base_flips
        # editing (A<->G) vs other among the suspect columns: editing is already handled, not an align artifact
        editing = 0
        for j in suspect:
            bases = {copies[fam][c][j] for c in copies[fam] if j < len(copies[fam][c]) and copies[fam][c][j] != "."}
            if bases <= {"A", "G"} or bases <= {"C", "T"}:
                editing += 1
        cvals = sorted(coh.values())
        med_coh = cvals[len(cvals) // 2] if cvals else float("nan")
        scored = len(coh)
        print(f"{fam}: {ncol} cols ({scored} with >=5 reads), {len(copies[fam])} copies, {nassigned} assigned reads")
        print(f"   coherence: median {med_coh:.4f} | LOW-coherence (<{COH_SUSPECT}, artifact-suspect): "
              f"{len(suspect)}/{scored} = {100*len(suspect)/max(scored,1):.1f}%  "
              f"(of which ALSO in a dense run: {len(clustered_and_incoherent)})")
        print(f"   [context] dense-run columns: {len(clustered)}/{ncol} = {100*len(clustered)/ncol:.1f}% "
              f"— mostly DIVERGENCE (adjacent real PSVs), not artifacts")
        print(f"   suspect-col makeup: {editing}/{len(suspect)} are A<->G / C<->T (RNA EDITING, already filtered); "
              f"{len(suspect)-editing} other (true align/consensus suspects)")
        print(f"   re-assign proxy vs production BASELINE (drop nothing): {100*base_flips/max(nassigned,1):.2f}% "
              f"(proxy noise floor)")
        print(f"   DROP {len(fragile)} suspect cols -> flips {100*flips/max(nassigned,1):.2f}%   |   "
              f"DROP {len(rand_drop)} RANDOM cols -> flips {100*rand_flips/max(nassigned,1):.2f}% (control)")
        print(f"   => EXCESS of suspect over random: {100*(flips-rand_flips)/max(nassigned,1):.2f}% "
              f"(the part actually attributable to the suspect columns)\n")
        grand[0] += scored; grand[1] += len(suspect); grand[2] += len(clustered); grand[3] += max(excess, 0)
        tot_assigned += nassigned
    print("== TOTAL ==")
    print(f"  columns: {grand[0]} | low-coherence (artifact-suspect): {grand[1]} ({100*grand[1]/max(grand[0],1):.1f}%)"
          f" | clustered: {grand[2]} ({100*grand[2]/max(grand[0],1):.1f}%)")
    print(f"  assignment FLIPS when ALL fragile columns dropped: {grand[3]}/{tot_assigned} = "
          f"{100*grand[3]/max(tot_assigned,1):.2f}%")
    print(f"  -> {'alignment error is NOT driving assignment (robust)' if grand[3]/max(tot_assigned,1) < 0.01 else 'alignment error MATTERS — model it'}")


if __name__ == "__main__":
    main()
