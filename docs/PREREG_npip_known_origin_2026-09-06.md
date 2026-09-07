# PREREG — one sensitivity / one specificity for NPIP on reads of known origin (2026-09-06)

The advisor's question is one pair of numbers for NPIP, with the same-alignment-score reads as the subcategory in
which the assignment must earn its keep. Real data has truth only for 62 audited anchors; a simulation from the
real NPIP loci has truth for every read. `bench/npip_sim.py`.

Design: for each of the 26 NPIP units on NC_073242.2 (`sweep_v14/fam_MCL1_073242/copies.fa`, the catalog's own
spliced unit sequences, dropped members included), 200 reads by `sim_genome.simulate_reads` (HiFi-like: 0.3 %
substitutions, 0.08 % indels, the existing simulator, seed fixed), names `SIMNPIP|<copy_idx>|<i>`; aligned
exactly as the substrate was (`minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`) to the 3-contig
FASTA; `copy_assign` on the family with its shipped defaults; scored by the read name.
Truth per read = the unit it was simulated from. Baseline = minimap2's primary placement (the unit with the most
block overlap). Subcategory = reads whose primary MAPQ < 60 (the aligner declares them ambiguous); reported
alongside: reads with `as_margin = 0` (equal best alignment scores).
Definitions (multi-class, abstention allowed): **sensitivity = reads assigned to the true unit / all reads;
specificity = 1 − reads assigned to a wrong unit / all reads; precision = right / (right + wrong)**; the same
three for the baseline where "assigned" = the primary placement.
Predictions:
- P1 overall: O2 precision ≥ 0.99; O2 wrong < baseline wrong.
- P2 subcategory (MAPQ < 60): baseline precision ≈ the fraction a random choice among the tied placements would
  give (reported); O2 precision ≥ 0.98; O2 sensitivity ≥ 0.30 (the real-data 38.6 % on NPIP's contested reads);
  every wrong O2 call is a read from a unit whose expressed segment is within sequencing error of the called unit
  (the copy-13 → 12 wall) — listed.
- P3 the K = 0 class: reads from units with a ≥ 0.99 partner are tied or abstain, not assigned wrong.
Fail rules: P1 precision < 0.99 ⟹ the certificate admits wrong calls on clean reads: a defect, listed per unit
pair. P2 sensitivity < 0.30 ⟹ the real-data number was flattered by the sole-candidate class; report both.
Limits stated up front: reads come from the units' spliced sequences (no isoform variation beyond the chain, no
allelic variation, one error model); this measures the assignment, not the catalog.

## Amendment (2026-09-06 15:35, after arm A was read)
Arm A (0.3 % subst + 0.08 % del + 0.08 % ins = 0.46 % edits/base) is harsher than the substrate: real NPIP-proper
MAPQ-60 unit reads carry a median 2.1 edits/kb at their placement (p75 5.9, p90 23.7). Arm A results stand as
reported (precision 0.9966, sensitivity 0.352; 3,221 of 5,000 reads certificate-rejected, 93 % of the longest
quartile). Arm B: subst 0.0015, del 0.0003, ins 0.0003 (≈ 0.21 % edits/base, the real median), everything else
identical; predictions P1–P3 unchanged. Both arms are reported.

## Amendment 2 (2026-09-06 15:50) — the certificate's alignment, not the assignment, rejects true reads
Finding (arm A, then the exact unit sequences): aligned to its OWN locus extent with the read-star's minimap2 line
(`-x splice`), the error-free sequence of unit 3 gets 954 inserted bases, unit 14 528, unit 17 814 (an exon the
splice mode fails to splice, left as an insertion); against the spliced unit (`map-hifi`) each gets 0/0/0/0. The
origin certificate therefore rejects perfect reads at those units; in real data units 10, 12 and 25 have 0 unique
reads assigned and all rejected, unit 17 21 assigned / 34 rejected.
Change: a candidate explains a molecule by the better of its two forms — the genomic locus (splice-aware, §6fd)
or the expressed chain (spliced unit, `map-hifi`, §6fc); per candidate the alignment with fewer edits
(X + I + D + unaligned) is used for the origin certificate AND supplies that candidate's bases at the read's
positions for the pairwise columns. Default on; `--read-star-genomic-only` restores §6fd byte for byte.
Predictions: P4 real NPIP (`sweep_v14` inputs): unique-read rejections fall at units 2, 12, 17, 23, 25 and no
unit's assigned unique reads move to another unit; MAPQ-60 placement agreement stays 1.0000; the 62 anchors
stay 0 wrong. P5 simulation arm A rerun: sensitivity ≥ 0.60 (from 0.352) with precision ≥ 0.99; the units with
0 right (3, 10, 11, 12, 17) recover ≥ half their reads. P6 paired 35: assigned rises, agreement ≥ 0.9995.

## Status (2026-09-06 16:05)
P1 held (arm A: precision 0.9966; arm A under the two-form certificate: 1.0000). P2 held on precision (0 wrong)
and failed on sensitivity under the §6fd certificate (0.256), held under the two-form certificate (0.686). P3
held (19 K = 0 ties, 0 wrong at MAPQ 0). P4 failed (real data unchanged; the real rejections are read-throughs,
divergence and the error constant, row 722). P5 held (0.352 → 0.725, precision 1.0000). P6 pending (paired 35).
Arm B pending.

## Amendment 3 (2026-09-06 16:16, user) — the machinery is for the contested reads
Definition change, applied before arm B was read under it: the certificate machinery decides molecules with
primary MAPQ < 60; an uncontested molecule keeps its certified call when it has one and takes its placement
otherwise (`--placement-first` and `--no-placement-assign` are the two alternatives, both measured on arm B
post hoc: 92.7 % / 2.9 % wrong and 93.5 % / 0 wrong against 95.5 % / 0.12 % wrong). Arm B is rerun under the
default for the shipped table; the machinery-only arm B numbers stand as reported.

## Status (2026-09-06 17:05, final)
Arm B under the shipped default (two-form certificate, contested-only machinery, certified-first placement):
overall 95.5 % right / 0.12 % wrong / precision 0.9987; contested 84.0 % / 0 wrong; MAPQ 0: 16 right, 0 wrong,
23 K = 0 ties; unique 99.8 % / 6 wrong (the placement's residue). Posterior metric (§6fr): truth in the tie set
for 100 % of the contested reads, top-1 for 98.3 %, uncertified pairs 197/197. P6 (paired 35) pending.
