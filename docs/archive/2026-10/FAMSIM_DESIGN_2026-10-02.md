# famsim — controlled gene-family simulations that prove the condition they test

**Design spec, 2026-10-02.** Status: IMPLEMENTED the same day (`bench/famsim/`, 16 unit tests, two species' ladders run —
see `bench/FAMSIM.md` for usage, the tables and the two pipeline findings). Deviations from this spec while building:
`absent_from_reference` compares matched bases (not identity) so a copy that differs only by an inverted exon still proves
absent; `deleted_exon_absent` uses a 0.20 11-mer bar (short exons share a few 11-mers by chance); O2 scores unique reads
by the aligner's placement because the AS-tied gate skips them by design; `family_intact` is `n/a` with < 2 scoreable copies.

## 1. Why

The advisor's standing objection: the results pick easy cases, or are lucky. His own framing of a convincing
demonstration: take a gene A, make an identical copy A′ (100%), plant both in an artificial chromosome, show what
minimap2 does; then mutate A′ step by step — SNVs first, then *bigger* changes (more or fewer exons in A′,
inversions) — and show that the rest of the pipeline (family definition, copy assignment, missing-copy
flagging) still behaves, for *any* gene family and *any* ape genome used as the source. minimap2's own limits are
out of scope: the deliverable is a simulation whose condition is **provable from its outputs**, plus the pipeline's
measured behaviour under it.

What exists already (`bench/sim.py`) covers pieces: `tandem` plants k SNP-mutated copies of one two-exon chr20
gene; `missing-copy` makes reference-absent copies (SNP or exon-swap); `chromosome` simulates a whole annotated
chromosome. None of them takes an arbitrary gene from an arbitrary genome, none applies structural operators
(exon gain/loss, inversion, truncation, conversion), none produces a machine-checked proof of the condition, and
each has its own ad-hoc scorer. famsim replaces that with one spec-driven generator, one verifier, one runner and
one scorer, reusing `sim.py`'s read model (`simulate_reads`, `jitter`, `stable_seed`) and `lib.py`'s scorers.

## 2. What it is

A Python package `bench/famsim/` (stdlib + pysam; minimap2/samtools and the Rust binaries on PATH or `--bin`),
run as `python3 bench/famsim <command>`:

| command | does |
|---|---|
| `spec` | prints a scenario JSON template, or one of the built-in ladder scenarios, to edit |
| `make SPEC --out DIR` | builds the artificial genome, the truth, the reads; writes the manifest |
| `verify DIR` | re-derives the condition from the OUTPUT files alone and writes `verify.tsv` (PASS/FAIL per claim) |
| `align DIR` | maps the reads with the shipped minimap2 command |
| `run DIR [--stages ...]` | runs the pipeline stages on the BAM (assemble, families de novo, families guided, assign, flag) |
| `score DIR` | scores every stage's output against the truth; writes `score.tsv` and `summary.md` |
| `ladder --template ... --background ... --out DIR` | makes, verifies, aligns, runs and scores the whole built-in scenario ladder; one table |

Every step is deterministic given the spec's `seed` (stable seeds, never `hash()`).

## 3. Scenario spec (JSON)

```json
{
  "name": "exon_loss_d02",
  "seed": 7,
  "background": {"source": "fasta", "path": "/.../chr20.fa", "region": "chr20:20000000-20300000"},
  "template":   {"source": "annotation", "genome": "/.../chm13v2.0.fa", "annotation": "/.../x.gff", "gene": "NPIPA1"},
  "copies": [
    {"id": "A",  "pos": 50000,  "strand": "+", "ops": []},
    {"id": "A2", "pos": 150000, "strand": "+", "ops": [{"op": "snp", "rate": 0.02}, {"op": "exon_delete", "exon": 3}]},
    {"id": "A3", "contig": "sim2", "pos": 40000, "strand": "-", "in_reference": false, "expression": 0, "ops": [...]}
  ],
  "decoys": {"n": 3, "min_exons": 3},
  "reads": {"per_copy": 50, "err": 0.001, "indel": 0.0003, "jitter": 30, "trunc5_frac": 0.0, "trunc5_max": 0.3}
}
```

**background** — `fasta` (a real slice: realistic repeats and composition; the default) or `random` (`length`, `gc`).
With `fasta` the slice's own annotated genes carry no reads and are not in the truth, so they are invisible to every
mode; nothing conditions on them.

**template** — the gene A. `annotation`: any GFF3 (RefSeq/CAT, `Parent=` chains, pseudogene exons parented to the
gene) or GTF (`gene_id`/`transcript_id`), any genome (human, gorilla, chimp, orangutan — the file pair is the only
species-specific input), by `gene` name, `transcript` id, or `random` with `min_exons`/`max_span`. The gene's
longest transcript (most exons, tie: longest) is the canonical chain; the genomic span from first to last exon,
in transcript orientation (+), is the template. `synthetic`: `exons` and `introns` length lists, random sequence
with canonical `GT…AG` introns — for runs that must not depend on any real genome.

**copies** — each copy = the template after its `ops`, planted at `pos` on `contig` (default `sim`) on `strand`.
`in_reference: false` plants the copy on a contig that the aligner's reference FASTA does **not** contain (the
truth genome does; the reads do) — the O3 case. `expression` = reads for this copy (default `reads.per_copy`;
`0` = an unexpressed copy, the semi-guided case). `isoforms` (optional): `[{"skip": [3], "weight": 0.3}]` adds
exon-skipping isoforms to the full chain.

**decoys** — `n` further annotated genes (random, ≥ `min_exons`), unrelated to A, planted and expressed like
copies: the precision control (no decoy may join A's family; every decoy must come out as its own locus).

**ops** (applied in order; coordinates are copy-local, exons numbered 1..n in transcript order):

| op | DNA effect | RNA effect | proves |
|---|---|---|---|
| `snp {rate}` | substitutions at `rate` over the whole copy, splice dinucleotides protected | copy-specific variants (PSVs) | divergence ladder |
| `indel {rate, max_len}` | short indels, splice sites protected | PSV indels | — |
| `exon_delete {exon}` | the exon's bases removed; flanking introns merge (GT…AG kept) | fewer exons | "fewer exons in A′" (DNA + RNA) |
| `splice_kill {exon}` | donor GT → CT of that exon (2 bp) | exon skipped; DNA ~identical | "fewer exons" at RNA only |
| `exon_insert {after, length \| source_exon, offset}` | a new exon (random, or a copy of `source_exon`) inserted into intron `after` with `AG`/`GT` flanks | more exons | "more exons in A′" |
| `exon_shuffle {a, b}` | exons a and b swap genomic positions (introns stay) | exon order changes | the `missing-copy shuffled` case, generalised |
| `invert {exon \| intron \| span [a,b] \| "whole"}` | reverse-complement in place | `exon`: that exon dropped from the chain (its splice sites now face the wrong way); `intron`: silent; `whole`: same transcript on the other strand | inversions |
| `truncate {side: 5\|3, exons}` | the first/last `exons` exons and their introns are missing | partial copy | SD-derived partial duplicates |
| `convert {from, exon \| span}` | that segment copied from another copy's CURRENT sequence | mosaic | gene conversion |
| `intron_resize {intron, length}` | intron shortened/lengthened (ends kept) | silent | spacing effects |

Splice dinucleotides are never mutated by `snp`/`indel` (memory: tandem-sim hygiene); `exon_insert` keeps intron
flanks so every planted junction is canonical. The RNA chain of a copy is derived from its DNA model after all ops;
the manifest records both.

## 4. Products of `make`

```
DIR/
  genome.truth.fa      every contig, including in_reference:false ones (the individual's genome)
  genome.ref.fa (+.fai) the aligner's reference (in_reference:false contigs omitted)
  truth.gtf            one gene_id per copy/decoy, one transcript per isoform (exons, strand) — for the assembler/family scorer
  truth.gff3           gene/mRNA/exon records with ID=/Name=/Parent= — the guided mode's node set and mcl_families --gff
  copies.tsv           copy_id, kind (copy|decoy), contig, start, end, strand, in_reference, expression, n_exons, ops (JSON)
  reads.fq             names `copy_id|isoform|i`; HiFi errors, jittered ends, optional 5′ truncation
  reads.truth.tsv      read, copy_id, isoform, contig, strand, junction list (truth coordinates)
  manifest.json        the spec as resolved (template sequence source, every op with its realised coordinates, seeds)
```

## 5. `verify` — the proof of the condition

Independent of the generator's bookkeeping, from the files in §4 only:

| claim | how it is checked |
|---|---|
| planted sequences are in the genome | each copy's exon sequence read back from `genome.truth.fa` at `truth.gtf` coordinates equals the manifest's exon sequence |
| every junction is canonical | `GT…AG` at every intron of every truth transcript (on its strand) |
| pairwise identity between copies | `minimap2 -cx asm20 --eqx` of each copy's genomic sequence against every other; identity = matches/aligned columns, plus `de`; reported per pair; `identical` scenarios must give 1.000 |
| exon count per copy | from `truth.gtf`; `exon_delete`/`insert`/`truncate` must change it by the stated amount |
| deleted exon absent | the template exon's sequence has no ≥ 90 % hit inside the copy (`minimap2 -c --eqx`, or exact search for short exons) |
| inverted segment | the segment at the stated coordinates equals the reverse complement of the template segment (and its forward form is absent) |
| absent copy absent | the `in_reference:false` copy's sequence has no hit ≥ 0.98 identity × 0.9 coverage in `genome.ref.fa` |
| expression | read counts per copy in `reads.fq` equal `expression`; 0 for unexpressed copies |
| reads carry the condition | each read's true junction list (from `reads.truth.tsv`) matches its copy's isoform chain |

`verify.tsv`: `claim, copy, expected, observed, PASS|FAIL`. A ladder run stops at the first FAIL.

## 6. `align` and `run`

`align`: the shipped command (`sim.MM2`: `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`), `--threads` (default 2,
the machine rule), sorted + indexed BAM. A per-read placement table is written at `score` time, not here.

`run` stages, each a documented recipe, each optional (`--stages assemble,denovo,guided,assign,flag`):

| stage | command | reuses |
|---|---|---|
| `assemble` | `tools/rustle_pipeline.sh assemble` (as_table + `copy_assign --assemble-only --genome-wide --assembly-polish full …`, the shipped defaults, bridge regroup f1v2) | the driver, unchanged |
| `denovo` | `tools/rustle_pipeline.sh families` → `PREFIX.fam.clusters.tsv` + `PREFIX.fam.copies.tsv/.fa` | the driver |
| `guided` | gene/pseudogene regions of `truth.gff3` → `samtools faidx` → `minimap2 -x asm20 -c --eqx -P` all-vs-all → `mcl_families --paf --gff truth.gff3 --min-exonic-bp 1 --min-shared-exon-frac 0.60` | the recipe in `figures/_o1_recovery.py::guided_families` |
| `assign` | `copy_assign --bam --fasta --regions --families PREFIX.fam.copies.tsv --copies-fa PREFIX.fam.copies.fa` (`--copy-table catalog` switches to the legacy `gw_family_catalog` table) | the driver's assign recipe |
| `flag` | `missing_copy_flag --scan-only` + `--from-scan` with a splice `.mmi` of `genome.ref.fa` built on the fly; `--confirm truth=genome.truth.mmi` optional | the driver's flag stage |

The binaries come from `--bin` / `RUSTLE_BIN` (default `/mnt/linuxdisk/home/juanfraitu/rustle_target/release`). Nothing is
rebuilt. Every stage's stdout/stderr goes to `DIR/logs/`.

## 7. `score` — against the truth, per objective

All matchings are one-to-one (metric traps: never let two truth copies share one node; nearest-start/exon-overlap matching,
not a positional window).

**Alignment report** (per copy, and per read class): primary on own copy / on another copy / on the absorbing locus (absent
copies) / unmapped; MAPQ 0 share; AS-tied share (secondary within 0.98 of the primary AS); supplementary; soft-clip ≥ 50 bp;
junction recovery (true junctions present in the primary, exact); for `invert exon` and `exon_insert` copies, where the
affected exon's bases went (aligned / clipped / inserted). This is the "what minimap2 does" table; it is reported, not judged.

**O1 — loci (assembler)**: assembled loci (`gene_id` groups of `PREFIX.gtf`) matched one-to-one to truth copies by exonic
overlap (bipartite on shared exonic bp, `lib.bipartite_items`): per copy `found | split (k loci) | merged (with …) | missed`;
per copy the best transcript's junction precision/recall against the isoform chain and `exact_chain` yes/no;
`transcripts_per_copy`.

**O1 — families**: de novo: cluster of each matched locus (`fam.clusters.tsv`) → copies; guided: cluster of each truth node
directly. Report pairwise sensitivity / precision / F over truth pairs (copies of A are one family; decoys are singletons) and
`bipartite_families`. The headline is per scenario: `family_intact` (all expressed copies of A in one cluster, no decoy) and
`guided_minus_denovo` (the gap the two modes are meant to close, memory: reduce the difference between them). Unexpressed
copies count for guided only (de novo cannot see them by design).

**O2 — assignment**: `assign.assignments.tsv` → per read `correct | wrong | abstain` against its copy (as `sim.py tandem`
does), split by tied (MAPQ 0 / AS-tied) and unique reads; per copy. Headline: wrong = 0 among assigned, abstain rate on
tied reads, and the divergence at which tied reads start to be assigned.

**O3 — flag**: for each `in_reference:false` copy, the locus that absorbed its reads (from the alignment report) must carry a
missing-copy verdict; copies present in the reference must not. Reported as flagged / not flagged per copy with the verdict
string (format from `missing_copy_flag`'s output; if the binary's output changes, the parser is one function).

Outputs: `score.tsv` (long: `scenario, objective, metric, copy, value`), `summary.md` (one table per objective),
`reads.placement.tsv` (per read).

## 8. The built-in ladder (`scenarios.py`)

The advisor's progression, each a spec generated from one template and one background (so every rung differs from the previous
one in exactly the stated way):

| rung | copies | condition |
|---|---|---|
| `identical` | A, A′ | 100 % identical, 100 kb apart |
| `snp_<d>` | A, A′ | A′ at divergence d ∈ {0.001, 0.005, 0.01, 0.02, 0.05, 0.10} |
| `exon_loss` | A, A′(d=0.02, exon_delete middle) | fewer exons in A′ (DNA) |
| `splice_kill` | A, A′(d=0.02, splice_kill middle) | fewer exons in A′ (RNA only) |
| `exon_gain` | A, A′(d=0.02, exon_insert random 120 bp) | more exons in A′ |
| `exon_dup` | A, A′(d=0.02, exon_insert source_exon) | internal exon duplication |
| `shuffle` | A, A′(d=0.02, exon_shuffle 2,3) | exon order |
| `inv_intron` | A, A′(d=0.02, invert intron) | silent inversion |
| `inv_exon` | A, A′(d=0.02, invert exon) | exon inverted |
| `inv_whole` | A, A′(d=0.02, strand −) | whole-gene inversion |
| `truncated` | A, A′(d=0.02, truncate 5′ by 2 exons) | partial copy |
| `conversion` | A, A′(d=0.05), A′ then convert exon 2 from A | mosaic |
| `three_copies` | A, A′(d=0.01), A″(d=0.03) | k = 3, unequal divergence |
| `dispersed` | A on `sim`, A′ on `sim2` | separate contigs |
| `unexpressed` | A, A′(d=0.02, expression 0) | semi-guided case |
| `absent` | A, A′(d=0.02, in_reference false) | O3 case |
| `combined` | A, A′(d=0.03, exon_delete + invert intron + truncate 3′) | everything at once |

Every rung carries 3 decoys. `ladder` writes `ladder.tsv`: one row per rung with the verify status and the headline numbers of
§7. The same ladder is meant to be run on ≥ 2 templates from ≥ 2 genomes (e.g. a human NPIP copy and a gorilla
"titin-like" LOC, or a random gene each) — the "regardless of family or species" claim is the ladder table replicated
across them, never pooled across species (memory rule).

## 9. Not in scope

Readthrough molecules (the `chromosome rt` arm already measures them), expression-level realism beyond per-copy depth,
polyA/5′ cap modelling, a read-length distribution learned from a library (the `trunc5` knobs are the only truncation
model), and genome-only discovery (outside the two O1 modes by the user's scope decision).

## 10. Testing

`bench/famsim/test_famsim.py` (stdlib `unittest`, no minimap2, < 5 s): operator bookkeeping (exon coordinates after
delete/insert/shuffle/invert/truncate; splice protection under `snp` at rate 1.0; inversion = reverse complement; convert
copies the current donor), planting (exon sequences read back from the written FASTA at truth coordinates equal the model's),
reads (names and junction lists match the chain; unexpressed copies yield 0 reads; seeds reproduce byte-for-byte),
`verify` on a synthetic spec (all PASS) and on a deliberately corrupted FASTA (the right FAIL). An end-to-end smoke run
(`ladder` on a synthetic template, one rung) is documented in `bench/FAMSIM.md` and run once by hand here.

## 11. Decisions taken without asking (correct them)

1. Package under `bench/famsim/` (run as `python3 bench/famsim …`), not more subcommands in `sim.py`: the operator set,
   verifier and scorer are ~1,500 lines and `sim.py` already holds five unrelated simulators.
2. Absent copies live on their own contig that the reference FASTA omits, instead of excising a span (keeps every other
   coordinate stable and makes the truth genome a valid `--confirm` genome).
3. An inverted exon is dropped from the RNA chain by default (`keep_in_rna: true` keeps its reverse complement in the read).
4. Default copy table for `assign` is the families stage's (`fam.copies.tsv`), the de novo definition since 2026-09-25;
   `--copy-table catalog` is the legacy path.
5. The ladder's acceptance bars are NOT pre-registered here — the framework produces the tables; the user writes the
   `PREREG_*` with predictions before reading any `score.tsv`.
