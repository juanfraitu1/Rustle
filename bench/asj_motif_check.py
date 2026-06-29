#!/usr/bin/env python
"""asj_motif_check.py — genome-grounded certificate for the O3 (ASJ) splice-mechanism claim.

CONTEXT: an earlier draft of bench/asj_findings.md claimed the cleanest ASJ calls
(PSMD2, DAXX) sit "on the canonical splice dinucleotide (GT-AG)" and "create/destroy
the splice motif". That claim is GENOME-FALSE by one base. This script re-derives the
truth directly from GGO.fasta and pins it, so the doc can never silently drift back.

What it certifies (all re-derived from the reference, no trust in the prose):
  C1. Coordinate convention in bench/asj_calls.tsv is 0-based; the `donor` column is the
      first intronic base and `acceptor` is the first exonic base AFTER the intron (a-exclusive),
      i.e. the intron is the half-open interval [donor, acceptor).
  C2. The canonical splice dinucleotides are INTACT for the flagship calls:
        PSMD2 (+ strand): donor dinucleotide GG..->GT at 0-based [donor, donor+1].
        DAXX  (- strand): forward-strand AC (= revcomp of the minus-strand GT donor)
                          at the 3' (high-coord) intron end.
  C3. The anchor SNP sits ONE BASE OFF the invariant dinucleotide (donor-1 / the
      a-exclusive exon boundary), NOT on it.
  C4. Genome-wide: 0 of the 475 called anchors fall on a core splice dinucleotide
      (anchor never lands in {donor, donor+1, acceptor-2, acceptor-1}).

Run:  /home/juanfra/miniforge3/bin/python bench/asj_motif_check.py
Exits non-zero (assert) on any drift from the genome-verified facts.
"""
import csv
import os
import sys

import pysam

FASTA = os.environ.get("GGO_FASTA", "/home/juanfra/winloci_scratch/GGO.fasta")
CALLS = os.path.join(os.path.dirname(__file__), "asj_calls.tsv")


def comp(b):
    return {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}[b.upper()]


def revcomp(s):
    return "".join(comp(b) for b in reversed(s))


def fetch(fa, chrom, start0, end0):
    """0-based half-open [start0, end0)."""
    return fa.fetch(chrom, start0, end0).upper()


def main():
    if not os.path.exists(FASTA):
        print(f"SKIP: reference not found at {FASTA} (set GGO_FASTA)", file=sys.stderr)
        return 0
    fa = pysam.FastaFile(FASTA)

    rows = []
    with open(CALLS) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append(row)
    assert len(rows) == 475, f"expected 475 calls, got {len(rows)}"

    # ---- C1+C2+C3: PSMD2 (+ strand) — donor GT intact, anchor = donor-1 -------------
    psmd2 = next(r for r in rows if r["gene"] == "PSMD2" and r["acceptor"] == "195406955")
    anchor, donor = int(psmd2["anchor"]), int(psmd2["donor"])
    donor_dinuc = fetch(fa, psmd2["chrom"], donor, donor + 2)
    assert donor_dinuc == "GT", f"PSMD2 donor dinuc expected GT, got {donor_dinuc}"
    assert anchor == donor - 1, f"PSMD2 anchor expected donor-1 ({donor-1}), got {anchor}"
    # anchor base on reference == one of the two het alleles (T/G), and is NOT inside the GT
    anchor_base = fetch(fa, psmd2["chrom"], anchor, anchor + 1)
    assert anchor_base in (psmd2["alleleX"], psmd2["alleleY"]), anchor_base
    print(f"OK  C2/C3 PSMD2 (+): donor GT intact at 0-based [{donor},{donor+1}]; "
          f"anchor {anchor} = donor-1 (ref base {anchor_base}); GT unaffected by either allele")

    # ---- C2+C3: DAXX (- strand) — minus-strand donor = forward AC, intact -----------
    daxx = next(r for r in rows if r["gene"] == "DAXX")
    a_anchor, a_acc = int(daxx["anchor"]), int(daxx["acceptor"])
    # intron = [donor, acceptor); on - strand the biological donor (GT) is at the high end,
    # which on the forward strand is the AC at 0-based [acceptor-2, acceptor-1].
    fwd_donor = fetch(fa, daxx["chrom"], a_acc - 2, a_acc)
    assert fwd_donor == "AC", f"DAXX forward donor-end expected AC (revcomp GT), got {fwd_donor}"
    assert revcomp(fwd_donor) == "GT"
    assert a_anchor == a_acc, f"DAXX anchor expected at a-exclusive boundary {a_acc}, got {a_anchor}"
    # anchor is at the exon boundary, ONE base past the AC dinucleotide [acc-2, acc-1]
    assert a_anchor not in (a_acc - 1, a_acc - 2), "DAXX anchor must be off the dinucleotide"
    print(f"OK  C2/C3 DAXX (-): minus-strand donor GT intact (forward AC at 0-based "
          f"[{a_acc-2},{a_acc-1}]); anchor {a_anchor} = +1 off it (exon boundary)")

    # ---- C4: genome-wide — 0/475 anchors fall on a core splice dinucleotide ---------
    on_dinuc = 0
    for r in rows:
        anc = int(r["anchor"])
        d, a = int(r["donor"]), int(r["acceptor"])
        core = {d, d + 1, a - 2, a - 1}  # donor GT bases + acceptor AG bases (0-based)
        if anc in core:
            on_dinuc += 1
    assert on_dinuc == 0, f"expected 0/475 anchors on a core dinucleotide, got {on_dinuc}"
    print(f"OK  C4 genome-wide: {on_dinuc}/475 anchors on a core splice dinucleotide "
          f"(all are splice-region / extended-consensus, never the invariant GT-AG)")

    print("\nALL CHECKS GREEN — O3 mechanism = splice-REGION variants, core dinucleotide INTACT.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
