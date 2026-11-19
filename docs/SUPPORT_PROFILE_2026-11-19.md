# What do whole reads really support? A support profile per transcript (2026-11-19, PSU, register 1828)

User's question: the vault's IRFinder-S uses a CNN to mimic IGV inspection; do we need something like that to know what full reads
really support and how much of a truth interval's width they cover, and to keep reads that deviate even slightly from an interval
from introducing errors? Answer: no learned model. IRFinder-S itself uses the CNN only for SHORT reads; its long-read mode counts
directly (minimum intron depth, reads spliced at both ends), and it trains the CNN on long reads as truth. With Iso-Seq each read is
a molecule, so support can be counted exactly. Tool `bench/famwise_eval/support_profile.py` (+ tests); workspace
`goal_multicopy_dev/support` (linked from `r1828_support_profile`).

## The profile
Per transcript (an assembly or truth intervals), from the primary reads on its strand that overlap it:
- **compatible reads**: every read intron EQUALS a transcript intron, consecutive, no aligned block inside a transcript intron;
  from them: reads carrying the full chain, 5'- and 3'-complete reads (end within 50 bp), exonic width covered by >= 1 and >= 3 reads;
- per junction: exact reads, reads SHIFTED 1-3 bp, near-miss reads (4-20 bp off);
- every incompatible read gets one reason: shifted, near_miss, alt_site (a near-miss whose own junction is well supported: >= 2
  exact reads and >= 10% of the nearby transcript junction, i.e. another real splice site), skip, other_junction, retention;
- flags: WEAK_JUNCTION (a junction with 0 exact reads), SHIFTED_DOMINANT (1-3 bp shifted reads outnumber exact reads at a junction),
  THIN_EXACT (weakest junction < 10% of the strongest junction's exact reads; optional absolute floor `--min-exact`),
  NEAR_MISS_DOMINANT, UNCOVERED (< 90% of exonic width covered), NO_FULL_READ.

## What changed on the way, and why (all on development data)
1. First version counted junctions within +-3 bp as support. SIRV showed that this hides the very artifact asked about: Rustle's
   false SIRV transcripts are 1-3 bp junction shifts of true ones, carried exactly by 2-119 reads each, while the true isoform's
   thousands of reads were credited to them under the tolerance. They are NOT multimappers (user's question): every read behind
   them is MAPQ 60 with no secondary. Now only exact junctions vouch; shifted reads are counted apart.
2. Near-miss reads at SIRV's deliberately close alternative splice sites flagged true transcripts (17 of 61): `alt_site` separates a
   well-supported neighbouring junction from an artifact; retention is counted, not flagged (real retention isoforms exist).
3. An absolute floor (< 10 exact reads) separated SIRV perfectly but flagged 42% of annotated chains on HG002 (it measures library
   depth); THIN_EXACT is relative by default.

## Validation
| test | result |
| truth with one junction shifted by 8 bp | flagged: simulated arrays 18 / 18, SIRV 61 / 61 |
| truth with one junction shifted by 2 bp (inside the old tolerance) | flagged: 18 / 18, 61 / 61 |
| SIRV truth (deep library) | 4 / 61 flagged (2 SHIFTED_DOMINANT, 2 NEAR_MISS_DOMINANT) |
| simulated array truth | 6 / 18 flagged: the silent units (no reads), as they should |
| Rustle on SIRV, `final.gtf` (relative THIN) | false 23 / 28 flagged, true 2 / 48 (absolute floor 10: false 28 / 28, true 2 / 48) |
| Rustle on SIRV, `PREFIX.gtf` (relative THIN) | false 25 / 37, true 2 / 53 (absolute floor: 32 / 37, 2 / 53) |
| HG002 seen quarter, Rustle primaries (not a truth) | annotated chains 5.6% flagged, novel chains 14.8% (absolute floor: 42% / 82%) |

## Use
`support_profile.py --gtf X.gtf --bam X.bam --out PREFIX --multi-exon-only` (truth intervals: pass the truth GTF; deep spike-in
libraries: add `--min-exact 10`). Read the flags as "the reads do not vouch for this structure exactly", not as "false": a real
rare isoform can be THIN_EXACT, and a false chain that reads carry exactly (row 1801's duplicate-collapse class) passes.
Helpers: `shift_junction.py` (copies of transcripts with one junction moved), `label_by_chain.py` (true / false by exact chain).
Not wired into the assembler or the filter; a filter use would need its own row.
