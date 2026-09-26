# Pre-registration: genome-wide copy assignment on every sample (Figures 4 and 5)

**Written 2026-09-25, before any genome-wide copy-assignment number exists.** At the time of writing no sample has
a genome-wide copy catalog (`${work}/runs/<sample>/<sample>.cat.copies.tsv` is absent for all six samples; checked
with `find`), so no read has been simulated from a genome-wide catalog, and no genome-wide `copy_assign` run on
simulated or real reads exists. User: every publication figure must be genome-wide on every sample, because "a
handful of chromosomes could be seen as cherry-picking".

What Figures 4 and 5 show today, and what this document replaces:

| panel | today | genome-wide (this document) |
|---|---|---|
| Fig. 4, Fig. 5a-b (simulation) | human A119b chr16 catalog (1,418 copies); gorilla OR6737 chr20 (NC_073244.2) catalog (357 copies) | the genome-wide catalog of each of the six samples |
| Fig. 5c (real reads, other methods) | the chr16 NPIP family (26 copies), A119b | every multi-copy family of the genome-wide catalog of human A119b and gorilla OR6737; NPIP kept as a named inset |

The earlier pre-registrations stay as they are: `PREREG_o2_read_truth_2026-09-23.md` (simulation, per-contig
catalogs) and `PREREG_hard_locus_bakeoff_2026-09-09.md` (NPIP). Their numbers are development results for this
document and are never overwritten. Genome-wide numbers are new register rows.

## 1. Samples and exposure

| sample | species | reads | copy catalog | splice index | Fig. 5c |
|---|---|---|---|---|---|
| human_A119b | human | A119b IsoSeq, CHM13 v2.0 | genome-wide, `make.py runs --stage catalog` | `target.splice.mmi` | yes (lab StringTie, FLAIR, IsoSeq collapse) |
| human_testis | human | public testis IsoSeq (ERR13885926), CHM13 v2.0; **mapped without `-uf`** | same | same | no |
| gorilla_OR6737 | gorilla | OR6737 testis, mGorGor1 | same | `GGO.splice.mmi` | yes (lab StringTie, FLAIR, IsoSeq collapse) |
| gorilla_KB3781 | gorilla | KB3781 fibroblast cell line, mGorGor1 | same | same | no |
| chimp_PTR | chimpanzee | mPanTro3 | same | built by the `index` stage | no |
| orangutan_PPY | orangutan | mPonPyg2 | same | built by the `index` stage | no |

Numbers from different samples are never pooled, and neither are numbers from different species.

**Exposure.** The copy-assignment rules (the per-family copy-assignment test, the gate that admits reads with equal
best alignment scores, and the union test) were developed on human A119b chr16. Gorilla OR6737 chr20 (NC_073244.2)
was scored once as the held-out substrate on 2026-09-23 and is now "held out, reused". No other contig of any
sample has been used for a copy-assignment decision. Each sample is reported for the whole genome. Human A119b is
also reported for "genome minus chr16" (not used to develop the copy-assignment rules), and gorilla OR6737 for
"genome minus chr20 (NC_073244.2)".

## 2. Experiment A: simulation with known source copies (Fig. 4, Fig. 5a-b)

The objects are those of `PREREG_o2_read_truth_2026-09-23.md`; only the catalog and the scale change.

- **Copies.** Every copy of the sample's genome-wide catalog whose spliced sequence is at least 300 bp and whose
  family has at least two copies.
- **Reads.** `bench/sim.py copies CAT.tsv CAT.fa IDX sim 20260925`: min(100, max(10, n_reads)) reads per copy, HiFi
  model (0.1% substitutions, 0.03% deletions, 0.03% insertions), up to 10% trimmed from one end, then 0-30 bp from
  each end. The read name records the source copy.
- **Mapping.** The shipped command (`minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`) against the
  species' whole-genome splice index, in P read-disjoint parts (`--parts P --max-parts-per-call 1`). P = ceil(E / R),
  where E is the read count the simulator will produce, computed from the catalog before any read is simulated, and
  R = 17,000 reads per part for human and 4,400 for gorilla, chimpanzee and orangutan (each part at most about 8
  minutes at the measured 0.021 s/read human and 0.098 s/read gorilla, plus the index load). P is fixed by the
  catalog, so every call maps the same parts.
- **No sibling all-vs-all.** `--sibling none` skips the simulator's closest-sibling `asm20` all-vs-all of every copy
  (no figure uses it; its cost genome-wide is unmeasured). Identity bands come, as before, from the catalog's
  `max_family_identity`: identical (100% over the aligned segment), 99.5-100%, 99-99.5%, 98-99% and below 98%.
- **Catalog handed to `copy_assign`.** `_o2.derive_catalog`, as before: the source catalog minus each copy that was
  neither simulated nor covered by a simulated primary alignment.
- **Assignment.** `copy_assign --bam sim.bam --fasta GENOME --regions <every contig of the BAM header> --families
  cat.copies.tsv --copies-fa cat.copies.fa`, then the same with `--union-certificate`, current defaults, every
  `RUSTLE_*` variable unset. Run in shards (Section 4) when one run would exceed 8 minutes.
- **Scoring.** `bench/score.py reads --catalog cat.copies.tsv sim o2|u2 --per-read`, joined by `figures/_o2.py`
  (`per_read`, `crosscheck`, `twin_state`), unchanged.

### Metrics, per sample (and per identity band)

Over all simulated reads: reads whose primary alignment has MAPQ > 0, MAPQ 0, or is unmapped; the aligner's primary
on the source copy (MAPQ > 0 and MAPQ 0 separately).

Over the reads with MAPQ 0, for four readings:

1. **Truth-selected:** the test's result within the read's source family, known from the simulation (`score.py
   reads` OWN). A user cannot choose this row.
2. **Default table, any assigned row** (ANY on the default run).
3. **Union test** (ANY on the `--union-certificate` run).
4. **Aligner's primary:** the catalog copy its primary alignment overlaps most.

For each reading: the fraction assigned, the fraction correct among assigned with a 95% interval that resamples
source copies (`_o2.copy_level_ci`: 2,000 resamples, seed 20260925; no interval below 5 assigned copies), and the
wrong and two-loci counts. Also the share of MAPQ-0 reads with an NM-identical twin (`_o2.twin_state`).

**Per-read fate (Fig. 4 main panel), fixed now.** Every MAPQ-0 read falls in exactly one class, in the
truth-selected reading: no result in its source family; no decisive site in its source family, left unassigned;
at least one decisive site, left unassigned; assigned to the source copy (a catalog copy at the same locus counts);
assigned elsewhere (wrong or two loci). The MAPQ > 0 and unmapped reads are counted beside the bar.

### Claims tested, and what changes each

The dev values are human A119b chr16 (1,263 MAPQ-0 reads).

| # | claim | holds in a sample when | if it fails |
|---|---|---|---|
| A1 | The truth-selected reading makes (almost) no wrong call where it assigns (dev: 163 of 163 correct, 21 copies) | correct / assigned >= 0.99, and the copy-level interval's lower bound >= 0.95 when >= 5 source copies are assigned | caption gives the measured error rate for that sample; > 5% wrong withdraws "no wrong call" for every sample |
| A2 | Where copies differ by 0.5-2%, most MAPQ-0 reads are assigned (dev: 0.77 at 99-99.5%, 0.80 at 98-99%) | truth-selected fraction assigned >= 0.5 in both bands, each with >= 50 MAPQ-0 reads from >= 5 copies; bands below that size are "not tested" | the claim is restricted to the samples where it holds, named |
| A3 | Copies identical over the aligned segment are not assignable (dev: 1 of 858) | truth-selected fraction assigned <= 0.05 in the identical band | the caption reports the rate and the claim is dropped |
| A4 | The default per-family table is not a per-read assignment | accuracy among assigned (any row) < the aligner's accuracy on the same reads (dev 0.09 vs 0.50) | stated per sample |
| A5 | The union test abstains on MAPQ-0 reads because they have an NM-identical twin (dev: 0 assigned; 1,255 of 1,263 twins) | union fraction assigned <= 0.05 and twin share >= 0.90 | the caption gives the mix of twin states instead of "expected" |

Predictions (reported, not claims): more MAPQ-0 reads genome-wide than on chr16 alone; the aligner's primary on
the source copy for about half of the MAPQ-0 reads (dev 0.50).

The read model, the bands, the readings, the denominators and the bars above are fixed. None will change after a
number is seen.

## 3. Experiment B: the genome-wide hard-locus benchmark (Fig. 5c)

**Samples:** human A119b and gorilla OR6737, the two samples on which the lab ran StringTie 3.0.1, FLAIR 3.0.1 and
IsoSeq collapse on the same BAMs. Real reads carry no copy truth. What is scored is whether a method's transcripts
carry a molecule's exact intron chain inside a catalog copy.

- **Catalog:** the sample's genome-wide catalog, every family with at least two copies.
- **Rustle arm:** `copy_assign --bam <full BAM> --fasta GENOME --regions <every contig that carries a catalog copy>
  --families <shard catalog> --copies-fa <shard copies> --gtf --origin-drop-indels --threads 4`, current defaults,
  `RUSTLE_*` unset, in shards (Section 4). This is `copy_assign`'s own transcripts, not the assembler of Figs 1-3.
  Contigs that carry no catalog copy (for example chrM) are not swept: a placement there is invisible to the gate,
  as a placement outside the swept window was in the NPIP benchmark.
- **Hard set:** the molecules with a row in the shard's `assignments.tsv` (at least two placements within 50 kb of a
  catalog copy on one contig, with runner-up alignment score equal to the best) whose primary alignment (`-F 2308`)
  lies in the merged copy intervals (the universe of `score.py bakeoff-calls`, the same for every method).
- **Strata:** all hard molecules; contested (not origin-rejected, at least two candidate copies); and molecules
  whose exact chain is carried by at least two reads (`--min-mult 2`). The last stratum was post hoc at NPIP; it is
  pre-registered here as a secondary stratum.
- **Other methods:** the lab GTFs, restricted per shard to the transcripts that overlap the shard's copies
  (`bakeoff-calls` ignores every transcript that overlaps no copy, so the restriction changes no call).
- **Scoring:** per shard, a subset BAM `samtools view -b -M -L <shard copy spans +/- 100 kb>` (the NPIP `hsa16.bam`
  recipe); `score.py bakeoff-calls <gtf> <subset BAM> <shard copies, copy ids renumbered to be unique across
  families> --fuzz 0` (Rustle, StringTie, FLAIR) or `--fuzz 5` (IsoSeq collapse), with `--tx-out`; then
  `bakeoff-compare --assign` and `--min-mult 2 --bam <subset BAM>`. The copy ids are renumbered because
  `bakeoff-calls` keys copies by `copy_idx`, which repeats across families; without it a molecule carried in two
  families' copies would count as one copy.
- **Family of a hard molecule:** the family of its row with `primary_local = 1` (its primary alignment overlaps a
  copy of that family); failing that, the family of the row `bakeoff-compare` reads (the last one).

### Metrics, per sample

- Fraction carried per method and stratum, pooled over every hard molecule of the genome.
- Per family with >= 20 hard molecules: fraction carried per method (the points of Fig. 5c), the paired difference
  Rustle minus method (median, and a 95% interval resampling families with seed 20260925), and the share of
  families where Rustle carries at least as many as the method.
- Transcripts overlapping the copies, and the share carrying at least one read's exact chain (as at NPIP).
- NPIP inset: the recorded chr16 benchmark, unchanged, labelled development.

### Claims tested, and what changes each

| # | claim | holds in a sample when | if it fails |
|---|---|---|---|
| B1 | On the pooled hard set Rustle's transcripts carry more molecules than StringTie's and FLAIR's (NPIP 0.858 vs 0.441 and 0.473) | Rustle > each of the two, pooled | the caption states which method carries more, in which sample |
| B2 | Per family, Rustle carries at least as many hard molecules as StringTie and FLAIR in most families | share of families (>= 20 hard molecules) >= 0.75 for each | the share is reported and the word "most" dropped |
| B3 | Rustle versus IsoSeq collapse: no direction predicted (NPIP: all hard 0.858 vs 0.732, contested 0.692 vs 0.813) | reported per stratum | none |

The hard set is defined by our own gate, so Fig. 5c measures how well each method's transcripts represent the
molecules our gate flags. It does not measure per-read copy accuracy (Experiment A does). This caveat is printed with
the panel.

## 4. Shards, and why they give the same result as one run

`copy_assign` processes one contig region at a time and, with `--families`, loads only the records within 50 kb of
the supplied copies (`COPY_READ_PAD`). Two things couple families: the records loaded together in one region, and a
genome-wide table of molecules with an equal-best placement outside every supplied copy (it forbids assigning that
molecule anywhere later). A shard is therefore a set of whole **read-connected components** of families:

- (a) two families are joined when their read windows (copy +/- 50 kb) overlap on a contig;
- (b) two families are joined when one molecule has records (primary, secondary or supplementary) inside a window
  of each;
- (c) a family with copies on several contigs is joined, on each such contig that also carries single-contig
  families, to the single-contig family nearest to its copy, so that no shard sweeps that contig without windows.

Contigs with no single-contig family are swept whole in every shard, as in one run. A shard's regions are the
contigs of its single-contig families plus those contigs. Components are packed into shards whose estimated wall
time is at most 8 minutes (first-fit decreasing on the records inside their windows). A component is never split.

**Check V3, required before any sharded number is used:** on the current human chr16 simulation, 2 shards give the
same per-family table (the same rows as a multiset; row order differs) and the same `union_certificate.tsv` rows as
one run, with the same binary, for both the default run and `--union-certificate`, and `score.py reads` prints the
same tables. If V3 fails for a run type, that run type is not sharded: it runs unsharded where it fits in 10 minutes
and goes to the cluster otherwise. Known exception, stated in advance: the column `readthrough_into` is filled only
for the rows of the first region `copy_assign` writes (the table behind it is emptied there), so it can differ
between a shard and one run; no scorer reads it.

## 5. Stop rules

- **Projection above 24 hours for one sample** (the sum of the shard estimates of Section 4 plus the simulation
  mapping estimate): the sample runs on a seeded random sample of read-connected components instead, fixed here:
  `random.Random(20260925)`, fraction 0.20, stratified by the contig of the component's first copy (components on
  several contigs form one extra stratum), at least one component per stratum. Sampling whole components keeps every
  sampled family's result identical to the full run. Every table and caption then says "20% of read-connected
  components, seed 20260925", and the sampled families are listed.
- **A single component estimated above 10 minutes** cannot be split without changing its result. It goes to the
  cluster; until it is done, its families and reads are reported as "not run", with their share of the sample.
- **Memory:** gorilla simulated-read mapping peaks at 19.3-20.8 GB (measured), so it runs with nothing else.
- **V3 fails for `--union-certificate`:** as in Section 4.

## 6. What is reported per sample

- Fig. 4 (main): one bar per sample of the per-read fate of the MAPQ-0 reads, grouped by species, never pooled.
- Fig. 4s (supplement): the full UpSet of every sample, one panel per sample.
- Fig. 5a-b: coverage and accuracy by identity band, one small multiple per sample with shared axes.
- Fig. 5c: human A119b and gorilla OR6737, pooled fraction carried plus one point per family; NPIP inset.
- Each table carries the sample id, the catalog (path and size) and, where it applies, the sampled fraction.

## 7. What would change the claims (summary)

- A1 fails anywhere: "no wrong call" leaves the abstract; the error rate per sample is reported instead.
- A5's twin share falls below 0.90: the union test's abstention is no longer "expected", and the next step is the
  mix of twin states, not a new rule.
- B1 fails in a sample: the NPIP result is reported as NPIP-specific, and the genome-wide result replaces it in the
  text.
- A projection or component stop rule fires: the affected numbers are labelled as sampled or incomplete, never as
  genome-wide.

---

## Amendment 1 (2026-09-25, same day; still no genome-wide catalog or number exists)

**Check V3 passed.** On the human chr16 simulation (28,453 reads, 290 families, 1,400 copies), with one frozen
`copy_assign` binary (sha1 22b3deb4) and every `RUSTLE_*` variable unset, the plan of Section 4 forced to 2 shards
gave the same rows as one run: `o2.assignments.tsv` 3,767 rows and `u2.assignments.tsv` 3,767 rows identical as
row multisets, `u2.union_certificate.tsv` 993 rows identical, and `score.py reads` printed identical tables and
per-read files for both runs. The one run also reproduced the recorded 2026-09-25 dev tables byte for byte.
Scratch: `/mnt/linuxdisk/tmp/rustle_figures_dev/figs45/v3/`.

**But read-connected components are large on a segmental-duplication-rich catalog.** On chr16, 287 of the 290
families form one component (91,374 of 91,686 records inside read windows); the other shard held 3 families.
Relation (a) alone (windows within 50 kb) already joins 287 families; relation (b) alone (reads) joins 283. So on
chr16 sharding saves nothing, and genome-wide the largest component may by itself exceed the 10-minute budget. The
rules above already cover this (Section 5, "a single component estimated above 10 minutes"): such a component is
never split. It runs on the cluster, or locally once the main session approves one longer run, and until then its
families and reads are reported as "not run". Nothing in Sections 2, 3, 6 or 7 changes. The measured costs (123 s
default and 144 s union for 91,374 window records; 7-8 s for a 3-family shard that also sweeps 24 contigs whole)
replace the cost model's first guesses in `figures/_o2.py` (`SHARD_COST`).

---

## Amendment 2 (2026-09-25): experiment C, comparison with the alignment-score margin rule

**Written before any number of this comparison exists.** No margin-rule call has been computed on any simulated
read (development or genome-wide), and no genome-wide real-read comparison exists (checked: no file under
`${work}/o2sim` or `${work}/fig5` mentions the rule). The only earlier numbers are the chrY comparison on real reads
without truth (`docs/EICHLER_COMPARISON_2026-09-21.md`, register rows 935 and 936). User decision, 2026-09-25 12:20:
compare copy assignment with the Eichler lab's alignment-score margin rule and show whether Rustle does what that
rule does, and more.

### C.1 The two per-read answers

**Margin rule, MR(T).** For each read, take every mapped, non-supplementary alignment, primary and secondary
(`-F 2052`), over the **whole genome**. best = the highest alignment score (AS); second = the highest AS among the
read's other alignments (an alignment with the same AS as the best gives a margin of 0). A missing AS counts as 0
(the convention of `as_table`). MR(T) assigns the read to its best alignment when the read has **no other
alignment** or best − second ≥ T; otherwise it discards the read. **A read with no rival alignment is assigned**
(register 936). T = 10 is the headline; T = 1 and T = 20 are reported beside it, because T is a convention, not a
constant. The rule's placement is mapped to a catalog copy by the rule of the Fig. 5 aligner baseline: the catalog
copy with the largest raw overlap of that alignment's reference span (ties: the first copy in catalog order), or
"outside the catalog" when it overlaps none.

**Why the rule is computed from the BAM, not from `copy_assign --eichler-margin`.** Read in the source before any
number existed (src/bin/copy_assign.rs `as_evidence_per_read` and the `--eichler-margin` block):
1. `copy_assign`'s `as_margin` is region-local. Alignments on other contigs, and outside the 50-kb windows around
   the supplied copies, are invisible, so a read whose only rival lies there looks rival-free and is assigned. This
   is the error register 936 corrected, in the opposite direction.
2. Under `--no-as-tied-only` (required to see the reads the gate drops) supplementary records count as rivals
   (`as_evidence_per_read(.., exclude_supplementary = !no_as_tied_only)`), although a supplementary record is
   another segment of the same read, not an alternative placement.
3. `--no-as-tied-only` is not the shipped copy-assignment method, so its rows are not Rustle's answer.
4. `eichler_same_copy` compares our copy with `primary_local`, not with the read's best-AS alignment.

So the comparison adds **no `copy_assign` run on the simulation**: MR(T) is computed from the same `sim.bam` as
experiment A (`bench/score.py eichler --sim`). **Check C0 (implementation only, decides nothing):** on the
development human chr16 simulation, run `copy_assign --no-as-tied-only --eichler-margin 10` once and count the reads
whose `eichler_call` differs from the BAM-derived MR(10), by cause (rival outside the loaded windows; supplementary
counted as a rival; other). Whatever it shows, the figures use the BAM-derived rule.

**Rustle's answer, R.** The per-read answer of the experiment A runs, split exactly as Figs 4 and 5 split the reads:
- primary alignment with **MAPQ > 0**: the aligner's primary alignment. The copy-assignment step leaves these reads
  to the aligner (Fig. 4), so Rustle's answer is that alignment's catalog copy (largest raw overlap), or "outside
  the catalog";
- primary alignment with **MAPQ 0**: the copy-assignment step's result. Two readings:
  - **R-u, the union test** (the `--union-certificate` run, any assigned result: status `assigned`, not
    origin-rejected). This is the answer a user can apply, and it is the headline;
  - **R-s, the test scored in the read's source family\*** (the default run's result for the source family;
    `score.py reads` OWN). It needs the simulation's truth and is shown beside R-u, never as the headline.
  No result, or a result with no copy assigned, leaves the read unassigned;
- unmapped reads are unassigned by both rules.

**Correct** (simulation only): the placement's catalog copy is the source copy or a copy at the same locus
(`score.py same_locus`, as everywhere in Figs 4 and 5). A placement on another catalog copy, or outside the
catalog, is **wrong**. (Fig. 5's aligner baseline counts a primary outside the catalog as not placed; here a rule
that placed a read has placed it, so both rules count such a placement as wrong, symmetrically.)

**Same placement** (both rules assign the read): R's placement and MR's placement are the same alignment (MAPQ > 0
reads whose primary alignment is the unique best-AS alignment), or they fall on the same catalog copy or locus.
Anything else is a different placement.

### C.2 Simulation (experiment A's runs, all six samples)

Population: **every simulated read** (the denominator; MR applies to every read, so the population is not
restricted to MAPQ 0), per sample and per identity band of the source copy (the bands of Figs 4 and 5), plus all
bands together. Per read and per T, one of six strata:

| stratum | MR(T) | R |
|---|---|---|
| both, same placement | assigns | assigns, same placement |
| both, different placement | assigns | assigns elsewhere |
| margin rule only | assigns | leaves unassigned |
| Rustle only, aligner's placement | discards | MAPQ > 0: keeps the aligner's primary |
| Rustle only, copy-assignment test | discards | MAPQ 0: the test assigns |
| neither | discards | leaves unassigned |

Reported per stratum: reads, source copies, reads with margin 0, and the reads each rule places correctly. For MR(T)
and for R (R-u and R-s): fraction assigned (assigned ÷ all simulated reads of the row) and fraction correct among
assigned, with the 95% interval of Figs 4 and 5 (`_o2.copy_level_ci`: resampling source copies, 2,000 resamples,
seed 20260925; no interval with fewer than 5 assigned copies). The same numbers for the two "Rustle only" strata
separately.

**Containment.** C = (both, same placement) ÷ (reads MR(T) assigns). The exceptions are the two strata "both,
different placement" and "margin rule only"; each is reported with its reads and with which rule is right (truth).

Claims, per sample, at T = 10 and reading R-u unless stated:

| # | claim | holds in a sample when | if it fails |
|---|---|---|---|
| C1 | Rustle places every read the margin rule assigns at the same place (containment) | C ≥ 0.99 | the exceptions are listed by stratum and by which rule is right; "Rustle does what the margin rule does" becomes the measured C for that sample; C < 0.95 in any sample withdraws it for every sample |
| C2 | The margin rule is accurate where it assigns (prediction, not a claim about Rustle) | fraction correct ≥ 0.99 | reported as measured |
| C3 | The reads the rule discards and Rustle leaves at the aligner's placement (MAPQ > 0, margin < T) are placed correctly | fraction correct ≥ 0.99 and the interval's lower bound ≥ 0.95, with ≥ 5 source copies | the measured error of those placements is stated per identity band, the words "and more" are restricted to the bands where C3 holds, and keeping MAPQ > 0 placements with a small margin is flagged as a cost of Rustle's gate (new register row) |
| C4 | The reads the rule discards and the test assigns (MAPQ 0) are assigned correctly | R-u: fraction correct ≥ 0.99 with ≥ 5 assigned source copies; fewer copies: "not tested". R-s: A1's bar | stated per sample; R-s is never quoted without R-u |
| C5 | T sensitivity (descriptive) | none | the fraction assigned and correct at T = 1, 10 and 20 is reported; no T is chosen from the data |

Predictions (reported, not claims), from Figs 4 and 5 before any margin was computed: R-u assigns almost no MAPQ-0
read (A5), so the "Rustle only" reads will be mostly the aligner's placements (C3), and the test's own contribution
shows under R-s only. The "Rustle only, aligner's placement" stratum also holds MAPQ > 0 reads with margin 0 (equal
best scores with a MAPQ above 0); they are counted separately, because there the aligner's choice is a tie.

### C.3 Real reads (human A119b and gorilla OR6737; no truth)

- **Reads compared:** every read whose primary alignment (`-F 2308`) overlaps a copy of a family with ≥ 2 copies in
  the sample's genome-wide catalog (the universe of experiment B and `score.py bakeoff-calls`). The same reads for
  both rules.
- **MR(T):** best and second AS over every non-supplementary alignment of the read in the whole BAM, from the
  sample's `as_table` (`runs/<sample>/<sample>.molecules.tsv`: `best_as`, `second_as`, −1 when there is no second).
  The best alignment's place: the read's alignments that overlap a catalog copy are read from the BAM; when one of
  them has AS = best, the rule's copy is its catalog copy, otherwise the best alignment lies outside the catalog.
- **R:** MAPQ > 0: the primary's catalog copy. MAPQ 0: the union test's result from one extra
  `copy_assign --union-certificate` run per experiment-B shard (same shards, regions and catalog; current defaults;
  every `RUSTLE_*` variable unset); no result or no copy assigned = unassigned.
- **Reported, per sample and T = 1, 10, 20:** the six strata of C.2; agreement where both assign (same placement ÷
  both assign); the reads the rule discards, split by margin 0 and 1 to T − 1; the "Rustle only" reads split into
  the aligner's placement and the union test; per identity band of the copy the primary alignment overlaps.
- **No truth.** These numbers measure agreement and coverage, not correctness, and every caption says so.
  Expectation, not a claim: agreement where both assign ≥ 0.97 (chrY 0.971); a lower value is reported with the ten
  families that carry most disagreements, as flags. Agreement is largely by construction (for MAPQ > 0 reads both
  rules follow the aligner when the margin is large); the informative numbers are the share the rule discards and
  the "Rustle only" and "margin rule only" reads.

Stop rules: those of experiments A and B (Section 5). The union runs are one more `copy_assign` call per shard, with
the same cost model; a sampled or incomplete sample is labelled as such.

### C.4 Where it is shown

Supplementary Figure 5s (`fig5s_margin_rule`), tables `fig5s_margin_rule_sim` and `fig5s_margin_rule_real`:
(a) per sample, T = 10, R-u: all simulated reads split into the six strata, with each rule's fraction correct;
(b) per sample and identity band: fraction assigned and fraction correct among assigned for MR(1), MR(10), MR(20),
R-u and R-s; (c) real reads, human A119b and gorilla OR6737, the six strata at T = 1, 10 and 20.

### C.5 Development pass

The code is first run on the development simulations of Figs 4 and 5 (human A119b chr16 catalog, gorilla OR6737
chr20 catalog), whose `copy_assign` runs already exist. Those numbers are development results: they are shown as
such until the genome-wide runs replace them, and the claims C1 to C5 are decided on the genome-wide runs only.
