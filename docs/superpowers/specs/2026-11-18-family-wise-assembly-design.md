# Family-wise assembly: one transcript per tie class, a copy named only where a read decides (design, 2026-11-18)

> **Status:** approved by the user on 2026-11-18, with one rule added ("Scores" below). Nothing is built and no register row is taken. Numbers marked "count of 2026-11-18" are
> descriptive counts on tables and reads already on disk, made while writing this design (throwaway scripts, logs and
> the review notes: workspace `famwise_design_2026-11-18`, not in the repository). They are not registered results and set no bar. Four reviewers (facts against the
> register, builder, advisor's view, evaluation) read the first version; their corrections are in this one.

User (2026-11-18): "help me figure out a way to implement a family wise assembly". Advisor's position as relayed by the
user: assemblers are a dime a dozen and this one is mid; the real problem is the provenance of multi-mapping reads, a coin
toss at near ties; a normal assembler was the wrong object, and a family-wise assembly that borrows information across
copies would gain several-fold. User's results: near ties are a small share of reads, and seeding loci from near-tied
secondary alignments helps a little, not several-fold.

## Why (what the ledger already knows)
- **Borrowing was built once, at the advisor's request; its failure does not settle borrowing.** FAMASM (spec 2026-10-31,
  register 1479): multi-copy precision 14.0 -> 8.5 held out. Count of 2026-11-18 on the A119b run: 20,882 of its 21,262
  appended chains had two or more primary alignments of their own carrying the whole chain at that copy; the route that
  really borrowed admitted 356 chains, 4 of them RefSeq-exact. What failed was a second admission path. What borrowing
  itself yields was small there.
- **Every way of moving structure or reads between copies that was run is closed with a measured ceiling.** Sibling
  structure as witness, filter, guide or projection reaches about 2% of missed chains (1519, 1522, 1627-1629, 1632).
  Pooling reads over copies or samples fragments loci or re-collapses chains (1102, 1503, 1531). EM, threading, abundance
  and parsimony gave no assign-or-abstain rule at usable precision (1351, 1352, 1412, 1532, 1533, 1648, 1777); on exactly
  tied reads 1351 and 1412 do record picks better than a coin (57-66% against 37-48%; 83% against 44% in simulation).
- **The size of the provenance problem is measured** (1785, 1787, 1790, 1791), on the reads the truth can judge (69% of
  the sampled clear-lead reads on KB3781, 49% on HG002, which is a 25% read sample on another individual's reference).
  Reads whose best alignment leads by more than 10: the primary alignment is on the true copy for 20,610 of 20,611
  (KB3781), 14,578 of 14,596 (HG002 on CHM13), 23,316 of 23,316 (orangutan AG06213). Near ties are 1.3% (KB3781) and 5.7%
  (HG002) of family reads; 3 and 14 families hold 72% and 81% of them. In those 14 families 44% of the near-tied reads fit
  two or more copies of the individual's own genome equally: for them one copy of origin is not in the read.
- **When a read ties.** minimap2 `splice:hq` scores a match +1 and a mismatch -4 (checked on a toy alignment: one
  mismatch moves AS and `ms` by 5). An exact tie needs no difference between the copies inside the read; "within 10"
  allows two; the seeding band (AS >= 0.98 x best) allows about eight on a 2 kb read, 99.6% identity. Count of 2026-11-18
  on the r1773 tables, two-member families: the four pairs per gorilla sample above 99.8% identity have 0.85-0.98 of their
  reads inside the seeding band (pooled; one pair at 99.93% has 0.00-0.25), pairs below 99% under 0.01. Larger families do
  have members between 99% and 99.8%. On the individual's own genome a full-length accurate read decides its provenance
  once three differences fall inside it; an allele that matches neither reference copy breaks this (NSMCE1, review of
  2026-11-16). The same model gives 22% exactly tied reads at 99% identity for 150 bp reads and none at 2 kb: the coin
  toss is a short-read regime. The project has no short-read data.
- **What the seeding band already buys, at the tied loci only** (count of 2026-11-18, band-sweep assemblies of row 1767,
  RefSeq-exact chains at the loci of tie groups, default band against primary alignments only): KB3781 270 -> 296
  (x1.10, exact-tie groups) and 427 -> 490 (x1.15, within 10); A119b 959 -> 1,167 (x1.22) and 1,393 -> 1,650 (x1.18).
  Genome-wide that is +77 gained of 26,240 (0.3%), +161 net (OR6737, 0.6%), +263 net (A119b, 0.7%). On KB3781 the band
  holds the true placement of all 45 judged near-tied reads with one origin, so an oracle can add few reads the assembler
  does not already see. It could remove chains: at A119b's exact-tie loci the band adds 6,051 chains for 229 correct.
  How many of the others stand at a wrong copy is not known; E2 measures it.
- **Where the advisor is right: inside a tie class the locus is the wrong unit.** TSPY array, A119b on CHM13 (count of
  2026-11-18): 394 molecules have an alignment in the array; the default assembly writes 335 multi-exon transcripts
  there. With primary alignments only the same assembler puts a model on 10 of 35 copies (4 structures), with the default
  band on 34 (10 structures); StringTie 18 copies and 9 structures, FLAIR 12 and 5, IsoQuant 9 and 4. No table says which
  transcripts are one RNA aligned at many copies, nor which copies an RNA can come from. No origin truth exists for this
  sample.
- **The project's own review of 2026-11-16** (`docs/FAMILY_WISE_REVIEW_2026-11-16.md`): "Family-wise assembly has nothing
  copy-specific to borrow in the families looked at: either the primary alignment is already on the right copy, or the
  copies cannot be told apart from the reads. The honest product at tied copies is one structure set, the set of
  compatible copies and a family-level count." That product exists at LOCUS level (stage `eqclasses`), and the review
  adds: "No registered accuracy test exists." The transcript-level form does not exist. The review's one check of ranges
  against origin (four KB3781 pairs, section 4) already shows the exact band missing the NSMCE1 origin and the
  within-10 band holding it.
- **The family-level unit was counted once** (2026-09, one NPIP family, against O2's own certificate, scripts since
  deleted): transcripts at a copy without evidence 2 with the copy-set form (by definition), 14 without, FLAIR 4,
  StringTie 15, IsoSeq collapse 45; isoforms given one arbitrary address 0 (by construction) / 39 / 76 / 20 / 200.

A family-wise assembler has not won on chains, and the room measured for it is small (x1.1 to x1.2 at tied loci, already
in the default). What is untested, and what tied families need, is the unit of the output and its test against read origin.

## Goals, non-goals, anti-goals
- Goal 1: the family-wise product at transcript level, with exact accounting and stated guarantees.
- Goal 2: the registered test of copy-level claims against read origin, same rule for our GTF and the baselines' GTFs.
- Goal 3: an experiment that bounds what ANY provenance method could add to or remove from an assembly.
- Non-goals: more correct chains; deciding tied reads; a new default.
- Anti-goals (a failure even if it runs): the product sold as a chain-accuracy gain; a transcript dropped from a default
  product; a claim about a copy that only tied molecules reach; a bar that cannot fail, or one set after its numbers
  were seen.

## Approaches considered
- **A. Tie-class layer on the default assembly (chosen).** Every chain of the default assembly is kept. Emitted
  (locus, intron chain) pairs that the same molecule carries at tied alignments are joined into family transcripts; each
  molecule is counted once; a copy is named only where a molecule aligns best there alone; per-copy counts are ranges.
  It is `eqclasses` carried from loci to transcripts, plus the accuracy test. No chain is added or removed, so chain
  accuracy cannot fall. Trade-off: a change of unit and of accounting, not a new source of chains. It joins only copies
  that share a tied molecule; diverged copies stay apart (see the first optional arm).
- **B. Archetype space (the form of Pal 2026 and of Paraphase).** Reads of a family on one copy, structural isoforms
  there, sequence isoforms by variant fingerprints, then back to copies. It needs a lift for every copy (FAMASM's base
  lift had no table for 1,450 of 7,991 copies on A119b), it loses copy-private exons, and its sequence step is the
  existing opt-in `copy_phasing`. In A the tied molecule is the lift: a molecule aligned at copy A and at copy B says
  that its chain at A and its chain at B are one RNA. A's table has the shape of Pal 2026's (structures in columns; the
  rows are tie classes, not sequence signatures) plus a copy set. Pal 2026 states that reads and genome came from
  different individuals, so a one-to-one mapping of transcripts to copies was not possible; here the reads and the
  T2T haplotypes come from the same individual.
- **C. A stricter FAMASM** (reads certified by `assign` only, through the normal polish): the one borrowing variant no
  row ran. An unregistered pre-check (LEDGER_O1, 2026-11-04, A119b) found no member of 28 missed tied-only members
  reaching two certified reads. Not proposed.

## The layer
Input: `PREFIX.families.gtf` of a run (the file `eqclasses` reads; a baseline's GTF also works), the BAM, the score
table. Score, bands and locus rule are those of `tools/eq_classes.py`, unchanged: score `ms`; a record is in the band when
its score is within W of the molecule's best, W = 0 (`exact`) and W = 10 (`wide`); the primary alignment is always in;
locus = the gene of the GTF sharing most exon bases on the record's strand.
1. **Node** = (locus, intron chain) emitted at that locus. Transcripts that share locus and chain are one node and carry
   the same values. A band record whose chain equals a node's chain is at that node. A record with no intron, or with a
   chain not emitted at its locus, is at the locus's **remainder**; a record at no locus is `outside`.
2. **Class** = the set of members (nodes, remainders, `outside`) of one molecule's band records. One member: the
   molecule is `placed` there. Molecules with the same set of two or more members form a class.
3. **Family transcript** = nodes connected through molecules: two nodes are linked when one molecule (option
   `--min-link`, default 1) has both in its class. Remainders, `outside` and `repeat` loci link nothing. A node linked
   to nothing is a family transcript of one node. A class lying wholly inside one family transcript is `inside`; any
   other class (it holds a remainder, `outside`, a `repeat` locus, or falls under `--min-link`) is `open`. Molecules of a
   family transcript = placed molecules of its nodes + molecules of its inside classes, each once; open molecules are
   counted beside it, with the other members listed, and enter no total.
   Diagnostic for over-linking: family transcripts holding two chains of one locus (`collision`).
4. **Per node**, per band: `lo` = molecules placed there, `hi` = `lo` + molecules of every class holding it. Per locus,
   the same on loci: molecules with every band record at the locus (= `placed` of `eqclasses`) and with any (= `max`).
   No estimate inside the range. Status of a node, using both bands, no new parameter:
   `certified` (lo >= 1 within 10: a molecule fits nowhere else nearly as well), `best_only` (lo >= 1 in the exact band
   only: a molecule scores best here, another copy is within 10), `tied_only` (lo = 0 < hi in the exact band),
   `near_only` (hi = 0 in the exact band, > 0 within 10: only records below their molecule's best carry it; added with
   the plan, 2026-11-18), `no_band_record` (hi = 0 within 10: assembled from partial reads or from records below the
   band). A tied molecule whose tied record was not scanned (a `--contigs` panel) gets the member `unseen`: its class is
   open and it certifies nothing.
5. **Per family transcript:** `min_copies` = the fewest loci that contain every locus with a certified or best-only node
   and meet every inside class. A lower bound on expressed copies, never a call (1351, 1352).

Products (new files only; every existing product stays byte-identical):
- `PREFIX.famwise.reads.tsv`: each tied molecule's records within 20 of its best (locus, node or remainder, gap to the
  best). Windows 0, 5, 10, 20 are computed from it without a BAM pass. The `eq` tables stay at 10.
- `PREFIX.famwise.tx.tsv` (one row per node and band), `PREFIX.famwise.classes.tsv`, `PREFIX.famwise.families.tsv` (one
  row per family transcript and band: copy set, molecules, open molecules, nodes by status, `min_copies`, collision).
- `PREFIX.famwise.gtf`: the input GTF, every line kept, with `fw_tx`, `fw_status`, `fw_lo`, `fw_hi`, `fw_copies` and
  `fw_rep` on transcript lines, from the within-10 band. `fw_rep "yes"` marks one node per family transcript (most
  placed molecules, then coordinate order): the one-per-RNA view. It says nothing about the other copies.

Propositions, each a test. A molecule's origin node is its record at its copy of origin; x(U) = molecules whose origin
node is in the node set U.
- **P1 conservation** (no assumption). Per band every molecule with a kept record is placed, inside one family
  transcript, or open; per locus the counts of rule 4 equal `PREFIX.eq.loci.tsv`.
- **P2 invariance** (records fixed). No table changes when the primary flag moves between records of equal best score.
  minimap2 chooses among them by a hash of `--seed` and the read name: that choice is the coin.
- **P3 count bounds** (assumption A: the record at a molecule's copy of origin is in its band). lo(U) <= x(U) <= hi(U),
  and both ends are reached by some assignment of the tied molecules.
- **P4 certificates** (A). A certified node is expressed at that copy, or stands in for a copy the reference lacks or
  for an allele of the individual (where A fails). At a `tied_only` node no molecule's alignment scores separate the
  copy from another member of its class.
- **P5 what is identifiable** (no assumption). A total is exact for a family transcript without open molecules and for
  nothing smaller that a class crosses; nodes lying in the same classes only have their sum constrained.
A fails where the individual differs from the reference (NSMCE1 in the exact band) and cannot be seen on reads the truth
does not judge; E1 reports it with both bounds.

## Preview (counts of 2026-11-18; throwaway probe with the layer's rule)
- **KB3781**, default assembly of 2026-11-16 (74,576 transcripts), read at the loci of the older tie groups. Exact band:
  134 family transcripts of two or more nodes join 303 transcripts on 149 loci (60 copy sets, 127 on two loci, one
  collision); their 15,774 inside molecules are counted 31,571 times today; 57 have no node with a placed molecule.
  88% of those molecules are one gene, CD99 on the pseudo-autosomal X / Y (9 family transcripts, 13,868 molecules);
  without it 125 family transcripts and 1,906 molecules. About 3,900 further tied molecules have a record with no
  emitted chain and stay open. Within 10: 262 family transcripts, 613 transcripts, 109 copy sets, 25,505 molecules.
  MAGED4 pair (chrX): 24 transcripts are 12 family transcripts, no molecule carries one of these chains at one copy
  only, and every judged molecule matches both copies of the animal's own haplotypes (review of 2026-11-16).
- **Arrays, A119b on CHM13.** TSPY, within 10: of the 335 transcripts 247 are carried by a band record and form 13 family
  transcripts plus 2 single ones; 88 have no band record; no collision. Linked only by sets that two molecules share,
  the same 247 fall into 130 pieces, which is why the link default is one molecule. GAGE: 25 transcripts, 3 family
  transcripts plus 4 single ones.
- **Near-tied molecules with known truth**: stated provenance against the truth (bands on AS around the primary, a proxy:
  the truth tables hold no `ms`). "Equals" = the statement is the truth set.

| sample, truth kind | judged | equals: primary alone / exact band / within 10 | no true copy in it: primary / exact / within 10 |
|---|---|---|---|
| HG002 on CHM13, one origin | 600 | 0.92 / 0.90 / 0.15 | 47 / 28 / 3 |
| HG002 on CHM13, fits several own copies | 2,945 | 0 / 0.82 / 0.96 | 1 / 1 / 1 |
| KB3781, one origin | 45 | 0.89 / 0.87 / 0.09 | 5 / 3 / 0 |
| KB3781, tie (not split) | 786 | 0 / 0.07 / 0.99 | 0 / 0 / 0 |
| orangutan AG06213, one origin | 173 | 1.00 / 0.99 / 0 | 0 / 0 / 0 |
| orangutan AG06213, tie (not split) | 89 | 0 / 0.21 / 0.96 | 0 / 0 / 0 |

The exact band is the statement (it keeps what the aligner decided and widens only at true ties); the within-10 band is
the envelope the truth should lie in. No single window is right for both kinds, hence two bands and four statuses.

### Stage output (2026-11-18)
The stage as built (`tools/famwise.py`, switch `--famwise` of `tools/eq_classes.py`), whole samples, on the default
assemblies of 2026-11-16, bands on `ms`. Command: `bash tools/rustle_pipeline.sh famwise --bam B --fasta G --out PREFIX
--threads 64` beside a finished default run (`PREFIX.gtf`, `PREFIX.families.gtf`, `PREFIX.molecules.tsv`). The numbers
above stay as they were: the probe read only the loci of the older tie groups, the stage reads every locus. In both
samples the conservation check against `eq.loci.tsv` holds, the two `eq` tables equal those of the plain `eqclasses`
stage, and no molecule has a record the scan did not meet. Family transcripts below are those of two or more nodes.

| | KB3781 (gorilla) | A119b on CHM13 (human) |
|---|---|---|
| multi-exon transcripts | 74,186 | 234,581 |
| exact band: family transcripts / transcripts / loci | 137 / 309 / 155 | 1,236 / 2,974 / 1,038 |
| exact band: on two loci / on more / copy sets | 130 / 6 / 63 | 967 / 207 / 582 |
| exact band: tied molecules inside one family transcript / open | 15,777 / 5,710 | 16,489 / 101,473 |
| exact band: two transcripts of one locus joined (collision) | 1 | 92 |
| exact band: at least two copies must be expressed | 46 | 384 |
| within 10: family transcripts / transcripts / loci | 267 / 627 / 286 | 2,242 / 6,096 / 1,872 |
| within 10: inside / open | 25,514 / 26,046 | 36,031 / 235,668 |
| status: certified / best_only / tied_only / near_only / no_band_record | 70,076 / 153 / 211 / 99 / 3,647 | 197,148 / 1,773 / 2,468 / 2,364 / 30,828 |
| wall time, 64 threads | 346 s | 6,005 s |

- KB3781: the X / Y pair of CD99 is 13 family transcripts with 13,919 inside molecules (the largest 10,683).
- A119b, TSPY (`chrY:8940000-9822000`), within 10: 13 family transcripts, all across loci, none with a collision, join
  245 of the 335 transcripts; 2 more stand alone with a placed molecule and 88 have no band record (probe: 13 plus 2
  single ones on 247, 88 without). The largest joins 45 transcripts at 45 loci with 80
  molecules (5 placed, 75 inside, 9 open) and needs at least 6 expressed copies. Status of the 335: certified 2,
  best_only 9, tied_only 94, near_only 142, no_band_record 88.
- Outputs: `/anvme/workspace/v138da13-rustle/famwise_design_2026-11-18/run/{kb,a1}.famwise.*` (job script
  `run/famwise.sbatch`, jobs 4480754 and 4480755).

## Scores: `ms` first, AS beside it (user, 2026-11-18)
"minimap actually uses the ms tag, so whenever we use the AS we should also take into consideration the ms." minimap2
ranks a read's alignments on `ms`; AS can order them differently. The layer's bands are on `ms` already. Every place this
design still leans on AS gets the `ms` counterpart, both printed:
- Truth strata and bands (`neutral_truth.py` defines near ties and kept placements on AS): the scorer reads `ms` and AS
  of each record from the selected-reads BAM and reports every E1 number for bands on `ms` (the bar) and on AS (beside).
  Molecules that are near-tied on one score only are counted and listed (KB3781: 116 of 1,321 near ties on AS are not
  within 10 on `ms`; 20 of 30,000 sampled clear-lead reads are).
- The preview table above is on AS; its `ms` version is the first output of E1 and replaces it.
- The seeding band is AS >= 0.98 x best. Where a number depends on it (E2's input D, the "45 of 45" count), the share of
  the molecule's `ms`-tied records that the AS band leaves out is printed (a designer's count on a KB3781 subset: 222 of
  3,670 exactly `ms`-tied molecules have a tied record below the AS ratio, so the assembler is not coin-invariant there).
- E2's input O keeps records by origin, on neither score. E3 compares records, then both scores.

## Measurement (bars registered before any run; drafts here, the user sets them)
Truth: the catalogue-free tables (`bench/o2_truth/neutral_truth.py`, r1785 / r1787, tie split r1791). Limits, stated in
every result: annotation families only (no truth read touches CD99); near ties complete, clear-lead reads a sample
(30,000 per run); at most 12 loci per read; truth bands on AS, the layer on `ms`; HG002 is a 25% read sample on another
individual's reference; "judged" requires the origin among the kept placements, so `outside` reads are possible failures
of A (KB3781 0 to 134 of 965, HG002 4 to 224 of 3,765, orangutan 0 to 114 of 376 at within 10). KB3781 develops;
a quarter of HG002 reads not yet looked at judges (the truth job's read selection with another offset).
The layer runs on the truth runs' selected-reads BAM (it carries `ms`), reads table written for every molecule. A
molecule's origin node is its record on the truth's origin span.

**E1. Read origin, real reads.**
- Molecule: the exact-band statement and the within-10 envelope against the truth set, beside the primary alone, by
  truth kind, for W in 0, 5, 10, 20. Draft bar B1 (held-out HG002 quarter): the exact-band statement equals the truth
  more often than the primary alone, pooled 95% intervals disjoint and in more families than not. Seen on other reads
  (table above), so a replication. Envelope misses are printed with both bounds on A; no bar (the envelope holds the
  primary, so it cannot miss more often).
- Claim "chain X is expressed at copy Y" = a node of the arm's GTF with a judged supporter. Support = judged molecules
  with a within-10 record at Y carrying exactly that chain (one rule for every tool, as in the 2026-09 bake-off).
  Origins of clear-lead reads come from the truth table, not from their primary. TRUE: a supporter originates only at
  Y; else SHARED: a supporter's origin set holds Y; else PHANTOM. Arms: `PREFIX.families.gtf` of the default run split
  by status, and the StringTie, FLAIR and IsoQuant GTFs on disk. Draft bar B2: (i) the share of PHANTOM + SHARED is
  higher among `tied_only` nodes than among `certified` ones, intervals disjoint; (ii) `certified` + `best_only` nodes
  hold at least 90% of the default GTF's TRUE claims. Under 30 judged `tied_only` nodes the bar is undecidable and is
  reported so (default transcripts on truth regions: 657 on KB3781, 3,317 on HG002). The PHANTOM share of certified
  nodes is printed as a check, not a bar (nearly the same predicate as TRUE).
- Range: on judged molecules only, copies whose true count lies in [lo, hi], per band, with the width (no bar).

**E2. Ceiling of any provenance method.**
- Real reads, the selected family reads of the truth runs, three inputs: P = primary alignments only, D = default,
  O = each read with a truth entry kept only at its copy or copies of origin. Reads without a truth entry stay at their
  primary alignment in all three. Assemblers: Rustle and StringTie. Reported on the regions a tied or moved read touches
  (121 on KB3781, 150 on HG002) and genome-wide: chains and RefSeq-exact chains per input; R = RefSeq-exact chains, O
  over P; chains of D that O does not have, and how many of them are RefSeq-exact (what perfect provenance removes).
  Draft bar B3, on the tied regions: the advisor's expectation holds if R >= 2 for either assembler; ours if R <= 1.1
  for both; in between the number is reported and neither is claimed. Direction known beforehand: O differs from P for
  937 of 100,810 selected reads on KB3781, and the sampled truth over-weights near ties (in the advisor's favour).
- `bench/famsim` (the advisor's planted ladder): differences per transcript 0, 1, 2, 4, 8, 20 x expression 1:1, 10:1,
  silent copy x full-length or 3' fragments, plus one row of 8 identical copies; inputs P, D, O and the layer; 5 seeds.
  Metric: planted (isoform, copy) pairs recovered with their exact chain, and chains emitted at a silent copy. A cell
  is several-fold when the median of O over P is >= 2 (or P recovers none and O some). Draft bar B4: several-fold cells
  lie only at <= 2 differences per transcript or in the array row, and at least one exists (otherwise the simulator
  shows nothing and no claim is made). Deliverable: the map with the real family pairs drawn on it.

**E3. Re-toss.** The selected reads realigned with `--seed` 11 (control: must reproduce the records), 12, 13, 14.
Per arm, the claims that are not in every seed; for the layer, rows that differ. No bar: P2 is a unit test, and a new
seed can also change which records exist.

## Shape of the implementation (Python only, about 700 lines)
- `tools/eq_classes.py`: a switch that also writes the reads table and per-node counts of placed molecules from the same
  BAM pass (about 40 lines). Without the switch its outputs stay byte-identical (`cmp`).
- `tools/famwise.py` (about 180, no BAM): nodes, classes, family transcripts, ranges, statuses, `min_copies`, tables, GTF
  attributes. `tools/test_famwise.py` (about 90): P1-P5 on toy and random structures, P2 by re-flagging, the three
  hand cases of the builder review (two tied emitted chains; one emitted and one not; two records at one locus).
- Driver: opt-in stage `famwise` beside `stage_eqclasses`, named only, never in the plan of `all`; one case entry, a
  README row, the dispatch test (about 25).
- Bench: E1 scorer (about 150), origin-filtered BAM for E2 (about 90), famsim rows and the three-input runner (about 120).
- Cost: one `eqclasses` pass (KB3781 about 7 min on `ms`, A119b 49 min at 64 processes); HG002's full BAM is 147 GB and
  is not needed (truth runs use the selected reads). On a `--contigs` panel a molecule tied to an unscanned contig is
  kept from reading as placed by the `unseen` member.
- `tools/eq_classes.py` and the driver carry uncommitted edits of other sessions (2026-11-18): build on the working
  tree, stage own hunks only.

## Left to the builder
Column names beyond those above; whether the stage also writes the `eq` tables or reads them; the routine for
`min_copies` (the exact one takes under a second on today's groups).

## Calls made here, for the user to overrule
Stage name `famwise`. Link by one molecule (`--min-link 1`), with the collision count as the guard. Four statuses from
the two existing bands. No estimate inside a range. Reads table at 20. GTF attributes from the within-10 band. In E2,
reads without truth stay at their primary in all inputs. Every bar is a draft.

## Calendar, and what to cut first
Week 1 (to Nov 25): pre-registration, stage and tests, KB3781. Week 2: E1 on KB3781, freeze, held-out HG002. Week 3: E2
on real reads, E3. Week 4: the famsim map. Week 5: orangutan, clean-up. Week 6 (to Dec 30): advisor page, register rows.
Cut first: orangutan, IsoQuant in the re-toss, the expression axis of the map, the GTF attributes. Keep to the end: the
stage with its tables, E1 on KB3781 and HG002, E2 on KB3781 with both assemblers, the difference ladder with the array
row, one figure.

## Optional arms (each needs its own registered bar; none is part of acceptance)
- **Structure axis across diverged copies** (found by the advisor-view review; no closed row covers it). The layer joins
  only copies that share a tied molecule. The exon map passed its test (row 1595: collisions 18 of 3,739). Every emitted
  chain of a family rewritten as its chain of exon columns gives one structural isoform per column chain and a table
  isoform x copy with the layer's [lo, hi] cells: Pal 2026's structure axis with the copy column filled by decisive
  reads. It pools no reads and adds no chain. First number to get: the share of family chains with every exon in a
  column.
- Sequence types inside a tied-only family transcript (`copy_phasing` on its molecules): how many distinct sequences are
  expressed, without naming their copies.
- Reads decided by `assign` fed to `assemble` through the normal polish (approach C; ceiling small in the pre-check).
- Wider-band secondary alignments as junction support only (1767); a structure criterion on the loci a family's
  multimapping reads reach (1773).

## Out of scope
New chains at any copy; dropping a transcript from a default product; deciding tied reads; pooling reads across copies to
improve models; any rule that needs a second sample or a sister species.
