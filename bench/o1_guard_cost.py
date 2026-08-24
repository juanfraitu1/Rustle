#!/usr/bin/env python3
"""What is the orientation guard actually discarding — real antisense pairs, or the pipeline's own
placeholder strand?

CONTEXT. Re-scored on the shipped instrument, the guard is O1's LARGEST loss channel (167/1,135 = 0.1471
of SD-supported pairs, ahead of the coverage clause at 0.1295 and gamma at 0.0097).

THE MECHANISM UNDER TEST. A rep's `strand` comes from the gate's canonical-junction classification, so a
SINGLE-EXON rep has no junctions and no determinable strand -- and `denovo_assemble.rs:1010` stamps
`strand.unwrap_or('+')`. Measured: ALL 5,928 single-exon reps carry '+' (1.0000, zero '-'), while spliced
reps split 0.4867/0.5133. A rep's `seq` is stored in its `strand` orientation, so a rep wrongly marked '+'
is stored REVERSE-COMPLEMENTED, and its alignment to a correctly-oriented paralogue comes out
MINUS-strand -- precisely what the guard drops at `denovo_pipeline.rs:4473`.

=> If the guard's cost falls on pairs involving single-exon reps, it is not removing antisense biology;
it is removing artefacts of an unmeasured field, and the fix is to stop guessing rather than to relax it.

RULE, from the shipped `ggo.rule.tsv` (do not substitute another coverage form -- using the wrong one is
exactly how the gamma figure came out backwards):
  identity = 1 - de              >= 0.60
  coverage = ql<=tl ? (qe-qs)/ql : (te-ts)/tl   >= 0.50   (axis follows the denominator)
  a pair qualifies on ANY SINGLE record clearing both floors
  the guard additionally requires the PAF strand field to be '+'
"""
import collections, csv, sys

P = "/mnt/linuxdisk/home/juanfraitu/o1_reps2/dump/ggo.er._k11_w5.0.paf"
N = "/mnt/linuxdisk/home/juanfraitu/o1_reps2/dump/ggo.nodes.tsv"

nex, strand = {}, {}
for r in csv.DictReader(open(N), delimiter="\t"):
    nex[int(r["idx"])] = int(r["n_exon"]); strand[int(r["idx"])] = r["strand"]

plus, minus = set(), set()
n = 0
for line in open(P):
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        continue
    q, t = int(f[0]), int(f[5])
    if q == t:
        continue
    ql, qs, qe = int(f[1]), int(f[2]), int(f[3])
    tl, ts, te = int(f[6]), int(f[7]), int(f[8])
    de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
    ident = (1.0 - de) if de is not None else int(f[9]) / max(int(f[10]), 1)
    cov = (qe - qs) / max(ql, 1) if ql <= tl else (te - ts) / max(tl, 1)
    if ident < 0.60 or cov < 0.50:
        continue
    n += 1
    key = (min(q, t), max(q, t))
    (plus if f[4] == "+" else minus).add(key)

lost = minus - plus
print(f"qualifying records {n};  pairs with a '+' record {len(plus)};  with a '-' record {len(minus)}")
print(f"GUARD-BLOCKED PAIRS (minus-only): {len(lost)}\n")

both1 = sum(1 for a, b in lost if nex[a] == 1 and nex[b] == 1)
any1  = sum(1 for a, b in lost if nex[a] == 1 or nex[b] == 1)
none1 = len(lost) - any1
print(f"  involve >=1 SINGLE-EXON rep (strand never measured): {any1}/{len(lost)} = {any1/len(lost):.4f}")
print(f"    both single-exon : {both1}")
print(f"    BOTH spliced (strand junction-determined -- genuine antisense candidates): {none1} = {none1/len(lost):.4f}")

# COMPARATOR: the base rate of single-exon involvement among pairs the guard KEEPS.
anyk = sum(1 for a, b in plus if nex[a] == 1 or nex[b] == 1)
print(f"\n  COMPARATOR — kept ('+') pairs involving >=1 single-exon rep: {anyk}/{len(plus)} = {anyk/max(1,len(plus)):.4f}")
print(f"  enrichment of single-exon involvement in BLOCKED vs KEPT: "
      f"{(any1/len(lost))/max(1e-9,(anyk/max(1,len(plus)))):.2f}x")

print(f"\n  => restricting the guard to pairs where BOTH reps have junction-determined strand")
print(f"     (n_exon >= 2) would RECOVER {any1} pairs and still reject {none1}.")
