# Prediction before a held-out test of the tandem member-size rule (2026-11-19, written before OR6737 is opened)

Rule frozen at commit 89639974: `bench/famwise_eval/tandem.py member_sizes` (read-end modes per member, ends shared by >= 2 reads,
no member longer than a period, else the folded consensus; antisense and > 1-period reads set apart). Its two safeguards were
chosen on A119b (row 1826), so A119b is development; this is the first data the rule has not seen.

Held-out sample: OR6737 gorilla testis (`~/reads/GGO/GGO_mm.bam`) on the gorilla reference (`~/genomes/GGO/GGO.fa`), NCBI
annotation `$DENOVO_GGO_GFF`. Arrays: `undecidable.py arrays` on that annotation; TANDEM families defined generically
(no hand-picked list): an array with >= 3 protein-coding units on one strand whose annotated unit lengths have a coefficient of
variation <= 0.2. Members: protein-coding units with >= 5 sense reads. Yardstick: the member's longest annotated transcript span.
Compared: the rule, StringTie `-L` on the slice's primaries, Rustle `assemble --seed-pool primary` on the slice.
The same generic definition is also applied to A119b (development) and printed beside.

Prediction: on OR6737 tandem-family members the rule is within 20% of the annotated span for >= 60%, and above both StringTie and
Rustle by >= 10 points. Failing either is reported as a failure of the rule on held-out data. This is not a registered bar (the
user has not set one); it is my prediction, written down first.

## Result of the first held-out test (OR6737, generic definition)
The rule sized 12 / 20 members within 20% (0.60); StringTie 16 / 20, Rustle 15 / 20. **Prediction failed** (not above the
assemblers). The generic definition picked large diverged gene clusters only (ABCA, ADH, NLRP, CCL, ZNF; 20-100 kb genes) because
the gorilla NCBI annotation names testis-array copies LOC..., so name stems never form those arrays; on long genes the read-end
modes size members too short (reads cover fragments). Not a general member-size rule.

## Second held-out test (written before it is run; definition chosen after the first result, so disclosed as such)
Arrays from NCBI DESCRIPTIONS (a LOC gene's description with a trailing number / "-like" removed, e.g. "testis-specific
Y-encoded protein"), on one contig, neighbours <= 200 kb apart, >= 3 protein-coding units on one strand, annotated unit length
CV <= 0.2 AND median unit length <= 15 kb (short near-identical units, the population of the simulations and of A119b's CT47A /
TSPY / GAGE). Same members, yardstick and methods as above; rule unchanged (commit f2519c38).
Prediction: rule within 20% for >= 60% of members and not below StringTie or Rustle. Fewer than 10 members: inconclusive.
