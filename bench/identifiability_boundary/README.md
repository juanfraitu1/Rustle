# Identifiability boundary experiment

Demonstrates the read-linkage identifiability boundary for multi-copy gene-family
assembly — the result that answers the "a single distinguishing SNP makes it
trivial" objection by **delimiting** the trivial regime instead of hand-waving past it.

## Setup (whiteboard-derivable; no model)
Two equal-expression gene copies, identical over a shared region `S` of length
`L_S`, differing at exactly two marks flanking `S`:

```
            mark m1            shared region S (identical, L_S bp)      mark m2
  copy A:   ...G...     ===========================================     ...T...
  copy B:   ...C...     ===========================================     ...A...
```
A read of length `R` observes a mark only if its footprint covers it. A read that
covers **both** marks must SPAN `S` (needs `R` ≳ `L_S`) and thereby **links** the
left allele to the right allele *on one molecule*.

## Theorem
The copies' full structures are identifiable **iff some read spans `S`**. With no
spanning read, the data are equally consistent with the truth `{(G,T),(C,A)}` and
the **chimera** `{(G,A),(C,T)}` — an isoform present in **neither** copy. So a
method that emits full isoforms must guess, and is wrong ~half the time.

## What the experiment shows (`boundary_curve.png`)
- **greedy** (linkage-blind point estimate): uses a spanning read if one exists;
  else emits full isoforms by guessing → fabricates. Phantom-isoform rate ≈ **0**
  for `L_S < R`, jumping to ≈ **0.5** for `L_S > R`.
- **certified** (our rule: *join across `S` only if a read spans it, else abstain
  + emit a non-identifiability certificate*): phantom rate ≈ **0 everywhere**.
- The identifiable region (some read spans `S`) collapses at **`L_S = R`** — the
  boundary is read length vs. mark spacing.

Headline (R=100, coverage=30/copy, ε=0.01, 400 trials/point):
`greedy phantom rate: 0.000 (L_S<R) → 0.497 (L_S>R); certified ≤ 0.013 everywhere.`

The small (~0.01) certified residual is honest: a single error-corrupted spanning
read can mislink; requiring ≥2 concordant spanning reads to certify drives it to 0.

## Generalizations (stated, not simulated here)
- Marks include **splice junctions / exon presence**, not just SNPs: a skip-junction
  read can bridge `S` by skipping it (changes the cut topology) — structure as a
  phasing channel.
- **Unknown copy number** `K`: a cut leaves model order ambiguous, not just phase.
- On the **variation graph**, `S` is one shared node — simultaneously the
  information-sharing substrate (reconstruct `S` once, both copies use it) and the
  cut (which in-branch connects to which out-branch is undetermined without a
  spanning read). That duality is why the vg is essential.

## Run
```
python3 bench/identifiability_boundary/boundary_experiment.py
```
Deterministic (seeded). Writes `boundary_curve.png` + `boundary_data.csv` and
prints a worked single-family case (greedy fabricates; certified refuses).
