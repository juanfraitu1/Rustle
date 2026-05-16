# no_seed scaffold — studying refs whose chain never reaches a transfrag seed

## Purpose

A subset of missed reference transcripts (chr19 GGO: 17) have intron chains
that never appear in a `transfrag_seed` event — rustle never even constructs a
transfrag with that exact chain, so path extraction has nothing to extend. This
scaffold classifies *why* each one fails and is the starting point for fixes.

## Components

1. **`no_seed_analysis.py`** — the analysis tool. Inputs: parity JSONL,
   reference GTF, gffcompare `.refmap`. Classifies every missed multi-exon ref
   by the stage its exact chain disappears, then deep-dives the no_seed class.

2. **`transfrag_define` parity payload** (pipeline.rs ~12539) — extended with
   seed-eligibility fields: `trflong_seed`, `weak`, `usepath`, `n_nodes`,
   `guide`. Lets the tool explain why a *defined* transfrag was never seeded.
   Parity-only change; default output byte-identical (1746/1948 verified).

3. **`junction_accept` parity events** (pre-existing) — used to test whether
   each junction of a no_seed chain is an accepted graph edge.

## Running

```bash
RUSTLE_PARITY_LOG=/tmp/p.jsonl RUSTLE_PARITY_FILTER_CHROM=NC_073243.2 \
  RUSTLE_FLOW_RESIDUAL_SE=1 ./Rustle/target/dev-opt/rustle -L GGO_19.bam -o /tmp/x.gtf
gffcompare -r GGO_19.gtf /tmp/x.gtf -o /tmp/cmp
python3 tools/trace_analysis/no_seed_analysis.py \
  --parity /tmp/p.jsonl --ref GGO_19.gtf --refmap /tmp/cmp.x.gtf.refmap
# add --only-no-seed for just the no_seed table, --tid STRG.59.7 to deep-dive one
```

## Failure sub-classes (chr19 GGO, 1746/1948 baseline, 17 no_seed)

| sub-class | n | meaning | fix direction |
|---|---|---|---|
| `COMBINATORIAL break@X` | 13 | all junctions are valid graph edges and each appears in *some* transfrag, but no single transfrag has the exact full chain; longest exact consecutive run ends at junction X | per-read structural enumeration / flow-decomp combinatorics (same class as cassette arbitration, STRG.92.2/300.6) — the hard one |
| `ABSENT_JUNC` | 3 | one junction has **no** `junction_accept` event — never seen in any read CIGAR. STRG.223.4, .238.2, .319.6 | read-level CIGAR investigation; likely a junction-correction/alignment gap, not a rustle filter |
| `EXACT_CHAIN_IN_TF drop` | 1 | a transfrag with the exact chain **is** defined but `trflong_seed=false` (STRG.59.7: `weak=1 usepath=-1 abund=1.0`) | the most tractable: a seeding gate dropped a correct transfrag. Investigate the `weak` flag + seed promotion path (graph_build.rs:701, pipeline.rs:4618-4705) |

## Coordinate conventions (important — got this wrong once)

- Intron-string in `transfrag_define`/`transfrag_seed`/`path_extracted`/
  `pred_filter_stage`: `donor_exon_end+1 - acceptor_exon_start-1` (1-based
  inclusive intron span). Build from GTF exons as
  `f'{ex[i][1]+1}-{ex[i+1][0]-1}'`.
- `junction_accept` start/end: raw GTF `(exon_end, next_exon_start)` — NOT the
  +1/-1 intron form. Match junctions as `(ex[i][1], ex[i+1][0])`.
- Set-intersection of individual introns is NOT an exact-chain test. "All 20
  introns appear somewhere" ≠ "the 20-chain exists". Use the longest *exact
  consecutive* sub-chain (the tool does this).

## Results (2026-05-15)

1. **EXACT_CHAIN_IN_TF (STRG.59.7)** — ATTEMPTED, F1-NEGATIVE.
   `tf_chain_not_subchain_of_keeptrf` + `nsr_distinct_chain` rescue landed as
   opt-in default-off (`RUSTLE_NSR_DISTINCT_CHAIN_MIN_LONGCOV=0`,
   `RUSTLE_NSR_DISTINCT_CHAIN_MIN_INTRONS=4`). Mechanism works (STRG.59.7
   seeds→extracts→emits at low thr) but STRG.59.7 has abund=1, so every
   threshold permissive enough to recover it floods combinatorial noise
   (DC=1 → 1716/2016 F1=89.03). No net-positive operating point. Same class
   as terminal-RI / per-read-combos. Detail in
   `project_missed_98_characterization_2026_05_15.md`.

2. **ABSENT_JUNC (3)** — RESOLVED: NOT a bug, NOT recoverable in de-novo mode.
   Read-CIGAR probe (`/tmp/cigar_probe.py`) on STRG.223.4 / .238.2 / .319.6
   shows each "absent" junction is in fact **massively read-supported but
   micro-shifted 2-4 bp** from the StringTie-de-novo reference's exon
   boundary:

   | ref | exact-ref N-gap | read-majority N-gap | shift |
   |---|---|---|---|
   | STRG.223.4 | 31640418-31640805 ×19 | 31640418-31640807 ×903 | acc +2 |
   | STRG.238.2 | 34248792-34249730 ×39 | 34248790-34249730 ×371 | don −2 |
   | STRG.319.6 | 47691768-47700871 ×20 | 47691768-47700875 ×1247 | acc +4 |

   Rustle's junction-correction (correctly) snaps to the read majority, so
   `junction_accept` only ever emits the shifted coordinate — the ref's exact
   chain never gets a transfrag. Verified rustle emits the downstream exon at
   31640808 (the +2 majority), not the ref's 31640806. gffcompare needs exact
   intron match for `=`, so a 2-4 bp jitter in one junction permanently
   demotes the whole chain to `j`. Same flavor as alt-acceptor-jitter j-class
   (STRG.73.1). Inherent benchmark ceiling vs StringTie-de-novo ref; only a
   reference-guided junction-snap would recover, which defeats the de-novo
   benchmark. No code change warranted.

3. **COMBINATORIAL (13)** — still deferred per agreed plan. Needs per-read
   structural enumeration after flow extraction; multi-session, precision-risk
   (StringTie-denovo ref may not emit these either).
