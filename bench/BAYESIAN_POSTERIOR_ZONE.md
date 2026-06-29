# Bayesian posterior + consistent zone — localizing the unassignable (2026-06-29)

The hard gate emits `assign / abstain (Tied)`. The Bayesian complement keeps the **per-copy posterior** instead
of collapsing it, which **localizes even an unassignable read to a ZONE** — the set of copies (and, for a
co-located family, the genomic sub-region) it is compatible with — with a probability distribution rather than a
bare flag. This is the soft/Bayesian reading of the same object: the frequentist gate reports the integral
support; the posterior is the fractional (LP-relaxed, Theorem 7) optimum over that support.

## What it computes

- Per read, `assign_read` now returns `posterior` = `softmax(logl)` over the candidate copies (a UNIFORM prior —
  likelihood-normalized). For an *assigned* read it is ~one-hot at `best_copy`; for a *Tied* read it spreads over
  the consistent zone.
- The binary (`copy_assign --posterior`) writes `<out>.posterior.tsv`:
  `read_name, family_id, status, n_consistent, zone_chrom, zone_start, zone_end, posterior` — where the **zone**
  is the genomic extent of the copies above a 0.01 posterior floor and `posterior` is `copy:prob,...` sorted by
  probability. An informative prior is selectable: `RUSTLE_POSTERIOR_PRIOR=abundance` weights by the EM copy
  abundance (the natural home for a DNA **parCN** prior later); default uniform.
- **Default OFF = byte-identical** (the field is populated but only emitted under `--posterior`). Full lib suite
  688 green.

## Result on GWFAM10 (8,783 reads — where the unassignable mass lives)

Zone size (number of consistent copies) for the **5,267 Tied** reads:

| consistent copies (zone) | tied reads | |
|---|---|---|
| **2** | **3,382 (64%)** | localized to a tight 2-copy zone |
| 3 | 525 | |
| 4 | 15 | |
| 5 | 186 | a ~92 kb sub-region of the array |
| 6 | 2 | |
| 7 | 28 | |
| **8** | 1,118 (21%) | fully ambiguous (whole array) |

- **~64% of "unassignable" reads are actually confined to just 2 copies** — the binary `Tied` flag was hiding
  that most of them are nearly localized; the posterior makes the residual uncertainty explicit and small.
- Example Tied row: `... tied 5 NC_073228.2 144930560 145022553 2:0.200,3:0.200,4:0.200,0:0.200,1:0.200` — the
  read is excluded from copies 5–7 (the distal array members) and spread uniformly over copies {0,1,2,3,4} in a
  92 kb zone. Assigned rows are one-hot: `... assigned 1 ... 144970572 144970572 3:1.000`.

## Why this is the right framing

- It **respects the K=0 floor**: within a truly-identical sub-zone the posterior is flat (uniform), i.e. the
  information limit shows up *as* "posterior = prior", not as a fabricated pick. It never breaks irreducibility.
- It is **not `1/k` as a hard call** (which Canzar dislikes) — it is the honest *posterior*, reported as a
  distribution over a zone, with the option of an informative (abundance / parCN) prior.
- It turns "unassignable" into "localized to zone Z with distribution π" — strictly more informative, and it is
  the principled re-entry point for the DNA copy-number prior (parCN) that would sharpen the zone.

Reproduce: `copy_assign --posterior [...]` → `<out>.posterior.tsv` (add `RUSTLE_POSTERIOR_PRIOR=abundance`).
