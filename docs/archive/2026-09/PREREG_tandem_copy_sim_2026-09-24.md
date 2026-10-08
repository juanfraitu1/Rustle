# Pre-registration — tandem-copy simulations: what minimap2 and the pipeline do with near-identical adjacent copies

**Written 2026-09-24 (§6zg), before the simulator is run.** User's three questions: (1) copy A and a similar copy A′
further down, two exons each — can a read join exon 1 of A with exon 2 of A′? (2) the same in the other direction —
exon 2 of A′ with exon 1 of A; (3) tandem copies up to 99% similar — what happens to the reads generated from them?

## The simulator (`bench/tandem_copy_sim.py`)

A real two-exon transcript from the chr20 annotation supplies exon 1, exon 2 and an intron shortened to 800 bp with
its own donor/acceptor ends. Copy A is planted into a 200 kb chr20 background at 50 kb; k − 1 further copies follow
at spacing D (default 8 kb), each mutated to identity p over exons and intron (independent mutations per copy).
Reads: N per copy from each copy's own spliced transcript (`sim_reads`, HiFi 0.001 / indel 0.0003, ends jittered
± 30 bp, ≤ 10% trimming), named `copy<i>|<n>`. Mapping: the shipped `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p
0.1 --secondary=yes` to the synthetic contig.

**Per-read classes (primary alignment; exon blocks = the aligned blocks separated by `N`):**
`same_copy` (both exons in the source copy) · `other_copy` (both in another copy) · `cross_forward` (exon 1 in an
upstream copy, exon 2 in a downstream copy — question 1 when the source is the upstream copy, question 2 when it is
the downstream one) · `cross_backward` (exon 2 upstream of exon 1: impossible for one colinear alignment, so it can
only appear as a split/supplementary record) · `unspliced` (no intron in the primary) · `partial` (≥ 50 bp
soft-clipped) · `unmapped`. Also per read: MAPQ, AS-tied (secondary within 0.98 of the primary), supplementary.

**Pipeline on the same BAM (`--pipeline`):** `copy_assign --assemble-only` (chains emitted; a chain whose exons sit in
two different copies is a **chimeric transcript**), `gw_family_catalog` (copies found), `copy_assign --families`
(assignment: correct / wrong / abstain against the read's source copy).

## Predictions, committed now

- **Q1/Q2 — cross-copy chains exist and are asymmetric with identity.** The aligner's chain score compares
  (exon 1@A, exon 2@A) with (exon 1@A, exon 2@A′): the two differ by the per-exon mismatch difference (≈ (1 − p) ×
  exon length) against the extra gap cost of a longer intron, which in splice mode is small. Prediction:
  `cross_forward` ≈ 0 at p ≤ 0.95, a few percent at 0.98, **5–20% at 0.99 and 0.995**, and at p = 1.0 all placements
  are ties (MAPQ 0) with `cross_forward` a substantial share. `cross_backward` never appears as a single primary; at
  most as supplementary records, near 0%.
- **Q3 — tandem copies at 99%:** reads keep MAPQ 60 on average (4 mismatches per 400 bp separate the copies), a
  minority tie; the assembler emits at least one chimeric A→A′ chain once ≥ 2 reads share it (canonical motifs at
  both ends, so no junction rule removes it); the catalog finds k copies at ≤ 0.99 and collapses at 1.0; assignment
  accuracy among assigned ≥ 0.95 at ≤ 0.99 with abstention rising toward 1.0.

The purpose is the machinery and the measured table; the bar is that the table is produced for k = 2 and 3 over
p ∈ {0.90, 0.95, 0.98, 0.99, 0.995, 1.0}, with the chimeric-chain count per condition.

---

# OUTCOME (2026-09-24) — `bench/tandem_copy_sim.py`, runs in `/mnt/linuxdisk/tmp/gw22/tandem/` (`t.*.summary.tsv`)

Gene `rna-NR_161305.1` (exons 592 and 243 bp, intron cut to 800 bp), 50 reads per copy, chr20 background.

| layout | identity | reads | same copy | other copy | cross-copy chain | MAPQ 0 / AS-tied | assembler: transcripts (chimeric) | catalog copies / families | assignment |
|---|---|---|---|---|---|---|---|---|---|
| tandem k=2, 8 kb apart | 0.90–0.995 | 100 | **100** | 0 | **0** | 0 / 0 | 2 (0) | 2 / 1 | 100 unique, all correct |
| tandem k=2 | 1.0 | 100 | 54 | 46 | 0 | 100 / 100 | 2 (0) | 2 / 1 | 100 tied, all abstain |
| tandem k=3 | 0.90–0.99 | 150 | **150** | 0 | 0 | 0 / 0 | 3 (0) | 3 / 1 | 150 unique, all correct |
| tandem k=3 | 0.995 | 150 | 150 | 0 | 0 | 0 / 55 | 3 (0) | 3 / 1 | 150 unique, all correct |
| tandem k=3 | 1.0 | 150 | 68 | 82 | 0 | 150 / 150 | 3 (0) | 3 / 1 | 150 tied, all abstain |
| interleaved k=2 (exon-level duplication, 2 kb spacer, canonical flanks) | 0.90–0.995 | 100 | **100** | 0 | **0** | 0 / 0 | 2 (0) | **0 / 0** | — |
| interleaved k=2 | 1.0 | 100 | 50 | 50 | 0 | 100 / 100 | 1 (0) | 0 / 0 | — |

**Q1/Q2 — a read joining exon 1 of one copy to exon 2 of another: ⛔ the prediction was wrong; it does not
happen from mapping alone.** At every identity below 1.0 every read chains both exons inside its own copy, in
both layouts and in both directions (from the upstream and from the downstream copy). At identity 1.0 the two
chains tie and minimap2 still keeps both exons in ONE copy (a coin toss between copies, MAPQ 0), even in the
interleaved layout where the cross-copy chain has the shorter intron: the splice chaining's gap cost does not
decide between introns of a few kb. The only cross-copy chains seen in this study (11% at ≥ 0.99, first
interleaved build) came from a construction error — background sequence after exon 1 of the second copy, i.e. a
NON-CANONICAL own donor — and vanished once every planted exon kept its splice flank. So a cross-copy junction
in real data means the molecule itself is chimeric (template switching, gene conversion, readthrough) or the
copy's own splice site is broken; it is not something the aligner manufactures from two clean copies. This is
the same conclusion the readthrough audit reached (§6ze answer of 09-23: 0.8% of duplicated-pair readthrough
reads have a competing placement).

**Q3 — tandem copies up to 99% (and beyond):** reads keep MAPQ 60 and the correct copy up to 0.995 (a 0.98-AS tie
appears for a third of the reads at 0.995 with three copies, still placed correctly); at 1.0 everything ties and
placement is arbitrary. The assembler emits one transcript per copy at every identity, never a chimera; the
catalog finds every copy, even at 1.0 (tied reads land at each copy); assignment is correct for every uniquely
placed read and abstains on every tied read at identity 1.0 (no PSV exists). The realistic ceiling is therefore
not the aligner but the read: two copies identical over the read's span cannot be told apart by anyone.

**A pipeline finding the simulator exposed:** an **exon-level (interleaved) tandem duplication is invisible to the
family catalog** — the two copies' spans overlap, so `gw_family_catalog` sees one locus with two chains
("1 γ-quasi-clique block → 0 families (≥ 2 distinct loci)") and calls them isoforms. The assembler does emit
both transcripts. A duplication that lies inside another copy's span is a class the copy catalog does not
represent; the structural detector (addendum 2 of the missing-copy prereg) is where such a copy would surface,
as an exon the read carries out of order.

Registered as row 1095. The simulator is kept (`--layout tandem|interleaved`, `--copies`, `--identity|--sweep`,
`--distance`, `--intron`, `--transcript`, `--pipeline`) so any further "what does minimap2 do when …" question
is one command.
