# Cross-copy splicing: is a junction used differently by the copies of a family? Design (2026-11-18)

Status: the user chose tests A and C and answered the two open questions on 2026-11-18 (section "Rulings"). The bars
are as drafted here and become fixed in the pre-registration, which is written before any count on a judging sample.
Nothing below is built, registered or run yet.

## Intent and why
The advisor's objection (relayed 2026-11-18): Rustle offers nothing StringTie, IsoQuant, bambu and FLAIR do not. On
chain precision outside primates that is what the transfer test found (rows 1779 / 1780). The thesis objective is an
assembler that does well in multi-copy families, so the answer has to be something those tools do not output at all.

One such output already exists as an opt-in report and has never been tested: stages `exonmap` + `psi` write
`PREFIX.fam.psi.tsv`, where one junction is read across the copies of a family (key = family, anchor site, exon column
of the anchor, exon column of the other end), each copy's PSI carrying the range that tied reads allow. The other tools
have no coordinates shared between copies, so they cannot state "copy A skips this exon, copy B does not".

This design tests whether the differences that table reports are real, by two independent routes the user chose on
2026-11-18: a cause in the genome sequence (A) and the origin of the reads (C). It also adds the one class the table
cannot show today (a copy that never uses the junction) and a table of which tool emits the copy-specific junctions.

Framing (standing rule): the exon columns are the exon-level view of the family graph (`exonmap` writes one GFA per
family). The result is read on the family graph's coordinates and feeds the variation graph; it replaces nothing.

## What exists (read, not changed by this design unless said)
- `bench/denovo_eval/exonmap.py` (row 1595, PASS on internal coherence only): exon columns per family, classes core /
  accessory / private. Only the exons of each member's representative transcript are in columns.
- `tools/psi_bunches.py`: PSI of each intron within its bunch (introns sharing a donor, or an acceptor), counted from
  alignments; per band `placed`, `max`, `psi_lo`, `psi_hi`. Keeps introns with `max_w >= 3` in bunches of >= 2 kept
  introns.
- `tools/psi_family.py`: keys the psi rows by exon columns; `differs` = this copy's wide range lies wholly above or below
  another copy's, both bunches holding >= 20 certain reads. Its docstring states the gap this design closes: "a copy that
  never uses the junction has no row at all".
- Only gorilla OR6737 has been run (workspace `psi_2026-11-16`): 2,416 keys, 505 in >= 2 copies, 212 differ; by the gap
  between ranges > 0.2: 44 keys in 24 families (side script of 2026-11-16, reproduced by the design's reviewer), which
  section 1 turns into 32 LARGE calls (28 member pairs), 7 of them copy-specific before any silent row exists. A rough
  count of silent copies for this design gave about 15 more (session scratch script, not kept, not the tool's rule).
  OR6737 is therefore the development sample: seen, never part of a bar.
- `bench/denovo_eval/exonshare.py` (`blocks_from_cigar`, `project`): base-level projection of a position from one member
  to another through the pair's alignment record (`PREFIX.fam.loci.gated.paf`, CIGAR), the function `exonmap` itself
  uses. This design uses it to find a splice site's counterpart at the other copy.
- Splice-site strength: `ssfeat.py` of row 1789 (workspace `r1789_splicesites/`, not in the repository): PWMs built from
  the sample's own introns (donor 9-mer, acceptor 23-mer), log-odds score per site, with a self-test.
- Read-origin truth that needs no catalogue: `bench/o2_truth/neutral_truth.py` with `NT_REPAIR=1`; tables
  `r1785_truth_repair/{kb,hg}/reads.tsv` (KB3781 on its own reference; HG002 on CHM13, reads with hash residues 0-24)
  and `r1796_famwise_origin/hg_held/truth/reads.tsv` (residues 25-49), each with its selected-reads BAM; tie split
  `bench/o2_truth/tie_split.py`. The truth caps the reads per stratum (`CAPS`, override `NT_CAP`).
- Default assemblies of 2026-11-16 for all 12 samples: workspace `artifact_refresh_2026-11-16/apes/SAMPLE/def/x.*`
  (assemble only; the `families` stage has to be run on them). Tool GTFs: `apetools/SAMPLE/` (StringTie, FLAIR,
  IsoQuant), bambu `r1641_bambu/` (OR6737 and A119b only).

## Goals
1. `PREFIX.fam.psi.tsv` also reports the silent copy: a member that owns both exons, has reads leaving the anchor, and
   never uses the junction; and a signed gap per key so a size can be read without a side script.
2. Test A: where the table says copy U uses a junction more than copy S, the splice sites say the same more often than
   chance.
3. Test C: the reads behind each copy's PSI came from that copy, and the difference survives when only reads of known
   origin are counted.
4. Tool table: for copy-specific junctions, which assembler emits the junction at the copy that uses it.

## Non-goals
- No borrowing of reads or structure across copies (closed: rows 1627-1632, 1643, 1649, 1696, 1777; user 2026-11-16).
- No test between samples, no statistical "differential splicing" (no replicates). Wording: "junction use that differs
  between copies in one sample, with the range tied reads allow".
- No replication in a second individual (option B, not chosen). No rule may need a second sample or species.
- Nothing becomes a default or part of `all`. `PREFIX.psi.tsv` stays byte-identical.
- No annotation is an input of any call or verdict here. RefSeq appears nowhere in A or C, and C's reads are selected
  by the calls' own loci, not by annotation families.
- Intron retention stays unseen (the psi stage counts spliced reads).

## Anti-goals (a failure even if it runs)
- A call that rests on reads which cannot be placed (the tie range exists to prevent this; C checks it).
- A PASS that comes from older families simply having more differences of every kind. A's bar is therefore on the
  DIRECTION of the site difference, not on its presence.
- Tuning any threshold on a sample that later judges it.

## Constraints
- Shared repository: commit by explicit path, own hunks only; shared files read-compute-write; register row numbers
  claimed by creating the workspace directory first (working rules).
- Pre-registration before any count on a judging sample (`docs/PREREG_*.md`); check
  `docs/NEGATIVE_RESULTS_REGISTER.md` first. Panel-first does not bite: `psi` on a whole sample is about a minute,
  `families` about 14.
- Big outputs in the workspace, intermediates in `$TMPDIR`, inputs read from the workspace copies, not atuin.
- Shipped code works for any species: no sample names, paths or human constants in `tools/`.
- Every AS-based count gets its `ms` counterpart; the psi bands are on `ms` already.

## Approaches considered
- **A. Sequence cause** (chosen): reference-only, any species, gives a mechanism per call. Sees only causes at the
  splice sites, so it is a direction test on the calls that have a site difference; calls without one stay unexplained,
  not refuted.
- **B. Replication in the species' second individual** (not chosen): the second sample is a cell line where these
  families are barely expressed; few judged units, and tissue differences would read as failures.
- **C. Read origin** (chosen): direct answer to the provenance objection for this product. Expected to be easy on a
  self-reference sample (KB3781: primary right 20,610 of 20,611 at a clear lead), so the bar is read on HG002 against
  CHM13, where the reference is another individual and a third of near-tied reads have no home.
The user picked A and C on 2026-11-18.

## Design

### 1. Calls (shared by A, C and the tool table)
Made from `PREFIX.fam.psi.tsv` by one function, used by every scorer.
- Eligible copy for a key: >= 20 certain reads (wide band) in its deepest bunch under the key (the tool's existing
  `--min-bunch`), or a silent copy (below).
- One call per key: U = the eligible copy with the highest `psi_lo_w`, S = the eligible copy with the lowest `psi_hi_w`;
  gap = lo(U) - hi(S). Classes: LARGE gap > 0.2; SMALL 0 < gap <= 0.2 (printed, no bar); SAME = ranges overlap and the
  point PSIs (`psi_placed_w`) of all eligible copies lie within 0.05 (the comparison set).
- U's intron = the intron with the most certain reads at U under the key (a copy's range can sum several introns with
  shifted sites inside one exon pair). An intron sits in a donor bunch and an acceptor bunch, so it can make two keys:
  per (family, U's intron) keep the key with the larger gap.
  Amended while building (2026-11-18, before any judging sample was run): one event read from both copies (the
  inclusion key with U = a and the skip key with U = b at the same anchor column) is also one call: per (family,
  anchor, anchor column, the two copies) keep the larger gap. Both would get the same verdict in test A and count one
  site difference twice. A gap is rounded to 4 decimals before it is compared with 0.2.
- The unit that is counted in every pooled bar is (reference genome, U's intron, S): a call made by two samples on one
  genome counts once (the sample with more certain reads at U gives its numbers), and a unit that is also a call on
  OR6737 counts in no bar (KB3781 shares its reference).
- COPY-SPECIFIC = a LARGE call with hi(S) <= 0.05 (includes every silent call).

### 2. Product change (opt-in report only)
- Added rows. For a key with a row at some copy, every other member S of the family that owns both columns of the key
  and has >= 20 certain reads on introns leaving its anchor site gets a row when it has none, with its real counts on
  the key's junction: 0 reads (`silent` = 1, range [0, 0]) or the 1-2 reads the stage's 3-read floor hides today
  (`silent` = 0, its real range). The junction's position at S = U's two intron ends projected to S through the pair's
  alignment record (`exonshare.project`); no projection (gap in the record, or the ends land on different strands or
  out of order): no row, counted as unlocated in the log.
- The counts for added rows are not in `PREFIX.psi.tsv` (it drops introns under 3 reads and bunches of one intron). My
  call: `psi_bunches.py` gets an optional side output with the counts of every intron it drops, written by the `psi`
  stage when the exon columns exist and read by `psi_family.py`. `PREFIX.psi.tsv` stays byte-identical.
- `PREFIX.fam.psi.tsv`: the existing rows keep their 26 columns byte-identical (so `n_copies` and `differs` keep
  today's meaning: copies with a row of >= 3 reads). Added rows come after them. New columns at the end of every row:
  `silent`, `n_eligible` (eligible copies including added rows) and `gap` (signed, this copy's range against the
  farthest eligible copy's; the calls of section 1 read this). The stage log keeps its line and gains one for the
  added rows.
- Tests beside the tool (`tools/test_psi_family.py`, `tools/test_psi_bunches.py`), the driver header and the stage's
  log line follow.

### 3. Test A: sequence cause
- Site strength: the PWM log-odds of row 1789 (`ssfeat.py`: `windows`, `pwm`, `score`), PWMs from the sample's own
  emitted introns (`PREFIX.gtf`), moved into the repository under `bench/crosscopy/` with its self-test.
- For a call, the introns compared are the union, over U and S, of the introns of the key's bunch at each copy (keyed
  or not, so junctions into exons outside the representative transcript count as competitors). Each intron's non-anchor
  site is located at the other copy by projection through the pair's alignment record, as in section 2. A competitor
  that cannot be projected is left out when it holds under 10% of its bunch's certain reads; otherwise, or when the key
  junction itself cannot be projected, the call is UNLOCATED (counted, not judged). The strength at a projected
  position is read from the genome there, whether or not any read uses it.
- Advantage of the key junction at copy X: A_X = strength of the key intron's non-anchor site at X minus the highest
  strength among the other introns' non-anchor sites at X. dA = A_U - A_S. The prediction of a sequence cause is dA > 0.
- Amended while building (2026-11-18, development sample only): a junction no read carries at a copy exists there only
  as an image through the alignment, and the image can slip a few bases beside a gap. Its site is read as the strongest
  within 6 bases of the image. On OR6737, 2 of 32 agreeing LARGE calls were a 2-base slip onto an intact site; with the
  rule the development counts are AGREES 29, OPPOSES 5, NO SITE DIFFERENCE 44, UNLOCATED 11 of 89 LARGE calls. Cost: a
  site lost with another one within 6 bases reads as no difference.
- Verdicts at tau = 1 bit: AGREES (dA >= tau), OPPOSES (dA <= -tau), NO SITE DIFFERENCE (|dA| < tau), UNLOCATED.
  Printed also at tau = 0.5 and 2, and the count of calls where the canonical dinucleotide (GT / GC donor, AG acceptor)
  present at U is absent at S.
- Samples: scripts are debugged and frozen on OR6737. Judging samples: the other eleven (A119b, HG002, KB3781, PTR, PPA,
  PPY, PAB testis, the four cell lines), each needing `families`, `exonmap`, `psi` on its default run of 2026-11-16.
  Each is opened once, after the freeze. Units are counted once per reference genome (section 1).
- **Draft bar A1.** LARGE calls of the judging samples pooled, verdicts AGREES + OPPOSES: PASS when the AGREES share has
  a Wilson 95% lower bound above 0.5. Fewer than 30 such calls: undecidable, said so.
- **Draft reading A2 (no bar).** Share of calls with a site difference (AGREES + OPPOSES over located) for LARGE against
  SAME, overall and split at exon identity 0.98 of the two copies, to show how much of it is plain divergence. Per
  sample and per family tables; the share of LARGE calls that come from the three largest families.

### 4. Test C: read origin
- Samples: KB3781 (own reference, reported beside, no bar) and HG002 on CHM13 (the bar; both truth quarters, residues
  0-49: no rule of this design was fitted on any HG002 read). Calls come from the whole-sample `psi` run; truth from
  the tables above.
- Truth reads. The existing tables select reads by annotation families and cap each stratum (in the held quarter the
  clear-lead stratum keeps 30,000 of 118,910 reads), so they cover the calls thinly and unevenly. C therefore makes its
  own truth: every read with a record on a member locus of a LARGE or SAME call (selected by locus, not by annotation),
  no cap, scored by `neutral_truth.py score` with `NT_REPAIR=1` and the tie split, against the haplotype indexes on
  disk. This needs a small locus selector beside the truth script and one alignment job per sample; the selection sees
  the calls' loci, never a verdict. HG002: reads of every hash residue (ruling 1).
- A judged read = truth kind `one` (one origin copy) or `several` (tie, class SEVERAL COPIES), as in
  `bench/famwise_eval/records.py`. For a read counted in `placed_w` at copy X on an intron of a call's bunch:
  SUPPORTED = kind `one` with origin at X; AGAINST = kind `one` with origin elsewhere, or kind `several`. Kind `outside`
  (the read's home is none of its reference places: a copy the reference lacks) is printed beside as its own share and
  is not in C1's denominator (ruling 2); C2 counts SUPPORTED reads only, so an `outside` read never helps a call
  hold. Origin place to member: same contig and strand, shared exon base with the member's column exons.
- **Draft bar C1 (reads).** Reads placed at U or S on the bunches of LARGE calls, pooled: PASS when SUPPORTED over
  judged is >= 0.95 with a Wilson lower bound >= 0.90. Printed per kind, per copy role (U, S) and per family.
- **Draft bar C2 (calls).** A call is judged when U and S each hold >= 20 SUPPORTED reads in the bunch. It HOLDS when the
  PSI gap computed from SUPPORTED reads alone has the call's sign and at least half its size. PASS when >= 0.90 of
  judged LARGE calls hold; fewer than 20 judged calls: undecidable.
- Step 0, counts only, before any verdict: judged reads and judged calls. Under the floors: undecidable, said so; no
  second truth run.

### 5. Tool table (descriptive)
For COPY-SPECIFIC calls of every judging sample: per arm (Rustle `PREFIX.gtf`, `PREFIX.final.gtf`, StringTie, FLAIR,
IsoQuant; bambu where a GTF exists), the share with the key intron in an emitted transcript at U, and the share with the
homologous intron emitted at S (for a silent S: a junction no read supports there). Reading, as in row 1714: an arm is
behind another when its share at U is more than 0.05 lower on >= 30 units pooled over the samples both arms have;
otherwise level. Split by A's verdict.
This table can come out level; then the new thing is the table of calls, not the assembly, and that is what gets said.

### 6. Outputs
- Code: `tools/psi_bunches.py`, `tools/psi_family.py` (+ tests, driver header); `bench/crosscopy/` (calls, site cause,
  origin, tool table; self-tests in the files, as the neighbouring bench scripts do).
- `docs/PREREG_crosscopy_splicing_2026-11-18.md` with the accepted bars, three register rows (A, C, tool table),
  ledger entries, results section. Workspace directories named by the claimed rows.
- An advisor page (visual first) is offered after the results, not part of this design.

## Amended after the code review (2026-11-18, before the freeze; no judging sample run)
A fresh reviewer found defects that would have biased the bars; the rules changed as follows (tests and ledger:
`.superpowers/sdd/2026-11-18-crosscopy-evaluation/progress.md`, lines "Final:").
- Section 2, added rows: a copy's own introns between its exons of the two columns come first (real counts from the
  dump); the carried-over image only when it has none and only with both flanks in those exons. Before, a copy that
  used the junction 7 bases away, alone at its sites, was written as silent.
- Section 1: unit = (genome, U's junction, S's contig, strand, anchor site), since a copy's name differs between
  samples; a call is MISKEYED, and judged nowhere, when S uses the image of U's junction under another key; the event
  rule merges only the two orientations of one pair; the SAME test is rounded like the gap.
- Section 3: the 6-base search runs on both copies of a pair or on neither; any site without a window makes the unit
  UNLOCATED; alignment records against the two members' strands are not used.
- Section 4: the origin PSI sums the copy's introns under the key, as the stated range does; MISKEYED calls and, on a
  sample of the development genome, development units are left out.
- Section 5: the speaking sample is chosen before the copy-specific filter; MISKEYED units are left out.

## Rulings (user, 2026-11-18)
1. C on HG002 may take reads of every hash residue at the calls' loci ("yes"). This spends the so far unseen half
   (residues 50-99) at those loci.
2. C1's denominator is `one` + `several`; `outside` reads are printed beside as their own share ("lets print them
   beside").

## Left to the builder
- Layout of the single-intron side table and how the `psi` stage passes it.
- File split inside `bench/crosscopy/`; how calls are cached between scorers.
- Job scripts, thread counts, the order of sample runs.
- Exact `families` environment per species (SEDEF BED and exclusion sidecars exist for CHM13 only; the stage's
  behaviour without them is the default for the apes).

## Risks
- Few units: the richest sample seen gives 32 LARGE calls on the rows the stage kept (89 with the added rows of part 1), and the eleven judging samples sit on six genomes, so every
  floor is reachable only pooled, and may not be reached.
- A1 undecidable: young copies have identical splice sites, so few calls carry a site difference at all.
- A1 FAIL: the reported differences are not caused at the sites the test reads, or they are noise. Then the calls are
  quoted only with C's result, and not as mechanisms.
- C on HG002 depends on depth in the families; step 0 says so before a verdict exists.
- Exon columns hold only representative-transcript exons, so junctions into other exons are never called.
