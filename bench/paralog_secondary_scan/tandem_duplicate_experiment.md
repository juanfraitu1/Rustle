# Tandem near-identical duplicates: are they merged / do we lose transcripts?

Investigation answering the question "for copies that are very close, how do we
ensure we are not merging tandem duplicates that might be producing usable
transcripts?" — via a real-genome scan plus a controlled synthetic experiment.

## TL;DR

- **Structure is safe.** rustle does **not** merge near-identical tandem copies
  into one locus or drop a copy — they assemble as separate loci. No chimeras at
  realistic spacing. (Only a pathological gap=0, copies directly abutting, breaks
  it → 0 transcripts.)
- **The residual risk is ABUNDANCE, not structure**, and only at **~100%
  identity** with **uneven** expression: the silent copy can be over-counted
  (same identifiability wall as DAZ). It disappears by **≥1% divergence**.
- VG mode flags the non-identifiable case honestly (`family_identifiability
  "none"`), but a **bug** makes the `capacity_confidence` channel useless
  (stuck at 1.000) on minimap2 BAMs — see below.

## Part 1 — Real GGO scan: do the risk loci exist?

Method: extracted all 34,114 genes from `GGO_genomic.gff`; found adjacent,
non-overlapping, similar-size gene pairs (gap −0.5–15 kb, size ratio <3); tested
sequence identity by aligning each pair with `minimap2 -x asm20 -c`.

- **397 near-identical adjacent paralog pairs** (identity ≥90%, aligned ≥30%).
- The **100%-identical** ones are mostly **unexpressed** tandem arrays (0–2 reads)
  — e.g. arrays on NC_073227.2:211.5 Mb and NC_086017.1:85.6 Mb.
- **Expressed** near-identical tandems do exist: the **protocadherin-β (PCDHB)
  cluster** (NC_073228.2:144.9–145.1 Mb, ~10 single-exon genes, 90–97% identity),
  a near-identical array at NC_073247.2:59.7 Mb, etc.

PCDHB result (real reads, de-novo): rustle produces **one transcript per
expressed copy** at the correct positions (PCDHB2/5/8/15…) — it does **not**
collapse the cluster — with a few **low-coverage chimeric** transcripts spanning
adjacent copies (e.g. 145,007,993–145,026,155, cov 3.7). So at 90–97% identity,
copies stay separate; the artifact is low-cov chimeras, not merging.

## Part 2 — Synthetic experiment (controlled)

Design: a synthetic contig = flank · copyA · 2 kb spacer · copyB · flank. Each
copy is a 3-exon gene (canonical GT–AG introns), transcript ≈1500 bp. copyB is
copyA mutated to a target identity. Full-length, stranded, HiFi-error reads
(~0.6% error) simulated per copy at controlled depth; aligned with
`minimap2 -ax splice -uf`; assembled with rustle de-novo and `--vg`.
Scripts: `gen_synth.py`, `gen_reads.py`, `run_clean.sh` in this dir
(data dir `/tmp/tandem_test/`).

| condition | mode | true A/B | rustle A/B | flag | verdict |
|-----------|------|----------|------------|------|---------|
| id=100%, equal   | de-novo | 40 / 40 | 36.7 / 42.6 | – | both recovered, ~balanced ✓ |
| id=100%, unequal | de-novo | 40 / 3  | 16.9 / **25.7** | – | copyB **fabricated** (true 3) |
| id=100%, unequal | VG      | 40 / 3  | 19.8 / 19.8 | **none** | splits 50/50, **flagged** non-identifiable |
| id=99%, unequal  | de-novo | 40 / 3  | 39.6 / **3.0** | – | **correct** ✓ |
| id=99%, unequal  | VG      | 40 / 3  | 37.6 / **0.0** | full | copyB **over-suppressed** (3→0) |
| id=100%, equal, gap=0 | de-novo | 40 / 40 | — | – | **0 transcripts** (pathological) |

Reading it:
- **Both copies always assemble as separate loci** (gap ≥ 200 bp); no chimeras.
- **100% identity + uneven expression** is the danger: de-novo splits the
  ambiguous reads' primaries arbitrarily and **fabricates** the silent copy.
  VG can only split 50/50 because the copies share **zero distinguishing
  sequence** → there are **no unique reads to anchor** the prior (the hard
  identifiability floor). VG at least **flags** it (`family_identifiability
  "none"`, `copy_confidence 0.5`).
- **≥1% divergence removes the problem**: ~15 distinguishing bases over a 1.5 kb
  read let minimap2 assign each read to its copy; de-novo abundance is correct.
- VG can **over-suppress** a genuinely-expressed but low copy (the 99% case,
  copyB 3→0) — worth watching; likely the copy-support guard threshold on small
  read counts.

## Part 3 — Bug found: `capacity_confidence` is stuck at 1.000

The `capacity_confidence` / `abundance_min` / abstain machinery (flow-capacity
work) is silently inert on minimap2 BAMs:

- `capacity_confidence = anchored_coverage / coverage`; `anchored_coverage`
  counts reads with `em_anchored = (gap>0.8) || (max_sites>0 && gap>0.5) ||
  was_unique`, and `was_unique = (nh <= 1)`.
- **minimap2 does not emit the `NH` tag** (it marks multis with `tp:A:S`). With
  `NH` absent, `nh` defaults to 1 → `was_unique = true` for **every** read →
  every read counts as "anchored" → `capacity_confidence = 1.000` always.
- Confirmed in the synthetic 100%-identical case (EM log: "0 anchored mass …
  40 uncertain reads", yet output `capacity_confidence 1.000`) **and in the real
  DAZ VG run** (all 11 transcripts `capacity_confidence "1.000"`).
- The non-identifiability flagging that *does* work uses a different path
  (`family_identifiability`/`family_verdict`).

**Fix:** derive read multiplicity from rustle's own multimap detection
(the read appearing in the family's `multimap_reads` / `tp:A:S` secondaries),
not the absent `NH` tag — i.e. set `em_anchored`/`was_unique` from the VG
multimap set, or compute `capacity_confidence` from the EM's already-known
anchored-mass (which was correctly 0 here).

## Bottom line

We are not merging tandem copies or losing their transcripts. The open issue is
abundance at the ~100% identity floor (non-identifiable, same as DAZ), which VG
flags honestly; and the `capacity_confidence` channel needs the NH-tag fix to
become a real low-confidence flag.

Figure: `tandem_duplicate_experiment.png`. Risk-loci list: `/tmp/tandem_paralogs.json`.
