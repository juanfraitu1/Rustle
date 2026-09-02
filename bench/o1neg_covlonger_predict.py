#!/usr/bin/env python3
"""Predict, from the FROZEN human negative-panel dumps, which E_r edges the asymmetric
two-sided coverage clause (`RUSTLE_ER_COVERAGE_LONGER_FLOOR`) removes.

This is a PRE-REGISTRATION instrument: it is run BEFORE the ON arm, and its output is the
falsifiable prediction the ON arm is scored against.

⚠ **Score the POST-GUARD set, not the PAF.** The `*.paf` dumps are raw minimap2 output,
upstream of the transcript-orientation guard and the rest of the shipped path: on this panel
they carry 28 candidate pairs while the pipeline emits **3** edges. Predicting on the PAF
would forecast the removal of edges that no longer exist. `*.edges.tsv` is the authority;
the PAF is consulted only to recover the alignment offsets, which `edges.tsv` does not carry
(`cov_longer` was added to the dump in §6ba, after this panel was frozen).

Replicates the shipped predicate (`denovo_pipeline.rs`):

    ident      = 1 - de:f:   (else nmatch/blocklen)
    cov        = ql <= tl ? (qe-qs)/ql : (te-ts)/tl      # the SHORTER sequence
    cov_longer = ql <= tl ? (te-ts)/tl : (qe-qs)/ql      # the LONGER sequence
    edge iff ident >= 0.60 and cov >= 0.50 and cov_longer >= FLOOR

Usage: o1neg_covlonger_predict.py <dump_dir> [floor]
"""
import sys
import glob
import os

MIN_ID = 0.60
MIN_COV = 0.50
TOL = 5e-4  # edges.tsv rounds identity/coverage to 6dp


def records(paf):
    """Yield (qname, tname, qlen, tlen, identity, cov_shorter, cov_longer) per PAF record."""
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        q, ql, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
        t, tl, ts, te = f[5], int(f[6]), int(f[7]), int(f[8])
        nmatch, blocklen = int(f[9]), int(f[10])
        if q == t:
            continue  # self-hit (internal repeat), never an E_r edge
        de = None
        for tag in f[12:]:
            if tag.startswith("de:f:"):
                de = float(tag[5:])
                break
        ident = (1.0 - de) if de is not None else (nmatch / blocklen if blocklen else 0.0)
        if ql <= tl:
            cov, cov_longer = (qe - qs) / ql, (te - ts) / tl
        else:
            cov, cov_longer = (te - ts) / tl, (qe - qs) / ql
        yield q, t, ql, tl, ident, cov, cov_longer


def main():
    dump = sys.argv[1]
    floor = float(sys.argv[2]) if len(sys.argv) > 2 else 0.30
    print(f"# floor = {floor}   dump = {dump}")
    print("window\trep_i\trep_j\tlen_i\tlen_j\tidentity\tcov_shorter\tcov_longer\tverdict")

    n_edge = n_keep = 0
    win_any = win_keep = 0
    unmatched = []
    for ef in sorted(glob.glob(os.path.join(dump, "*.edges.tsv"))):
        wid = os.path.basename(ef).split(".")[0]
        rows = [l.rstrip("\n").split("\t") for l in open(ef)][1:]
        rows = [r for r in rows if r and r[0]]
        if not rows:
            continue
        win_any += 1
        pafs = glob.glob(os.path.join(dump, f"{wid}.*.paf"))
        R = list(records(pafs[0])) if pafs else []
        kept_here = 0
        for r in rows:
            ri, rj, ident, cov = r[0], r[1], float(r[12]), float(r[13])
            hit = [x for x in R if {x[0], x[1]} == {ri, rj}
                   and abs(x[4] - ident) < TOL and abs(x[5] - cov) < TOL]
            n_edge += 1
            if not hit:
                unmatched.append((wid, ri, rj))
                print(f"{wid}\t{ri}\t{rj}\t?\t?\t{ident:.6f}\t{cov:.6f}\t?\tUNMATCHED")
                continue
            # the pipeline keeps ONE record; take the most generous surviving cov_longer
            b = max(hit, key=lambda x: x[6])
            survives = b[6] >= floor
            kept_here += survives
            n_keep += survives
            print(f"{wid}\t{ri}\t{rj}\t{b[2]}\t{b[3]}\t{b[4]:.6f}\t{b[5]:.6f}\t{b[6]:.6f}\t"
                  f"{'SURVIVES' if survives else 'REMOVED'}")
        win_keep += 1 if kept_here else 0

    print(f"\n# post-guard edges: {n_edge} -> {n_keep} ({n_edge - n_keep} removed)")
    print(f"# windows emitting >=1 E_r edge: {win_any} -> {win_keep}")
    if unmatched:
        print(f"# ⚠ {len(unmatched)} edge(s) had no matching PAF record: {unmatched}")


if __name__ == "__main__":
    main()
