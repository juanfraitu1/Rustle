# PREREG r1802 / r1803 / r1804: splicing differences between the copies of a family, against the splice sites and against the origin of the reads (2026-11-18)

Written before any count on a judging sample. Workspaces `r1802_crosscopy_sites` (row A, and the stage runs of all
samples), `r1803_crosscopy_origin` (row C), `r1804_crosscopy_tools` (tool table). Design:
`docs/superpowers/specs/2026-11-18-cross-copy-splicing-design.md` (sections 1, 3, 4, 5 and "Rulings"); plans
`docs/superpowers/plans/2026-11-18-crosscopy-psi-rows.md` and `2026-11-18-crosscopy-evaluation.md`. The bars below were
drafted in the design; they are fixed by the commit of this file. Scripts: `tools/psi_family.py`, `tools/psi_bunches.py`,
`bench/crosscopy/{calls,sites,site_cause,origin,tool_table}.py`, frozen as a snapshot in `r1802_crosscopy_sites/snap/`
after that commit; every job runs the snapshot.

## Why
Advisor, as relayed on 2026-11-18: the assembler offers nothing that StringTie, IsoQuant, bambu and FLAIR do not. On
chain precision that is fair outside primates (rows 1779, 1780). What the other tools do not produce at all is a table
that reads one junction across the copies of a gene family: stage `psi` after `exonmap` (`PREFIX.fam.psi.tsv`), keyed by
the family's exon columns, with the range of PSI that tied reads allow at each copy. The table reports copies that
splice the same junction differently. Before that is offered as a result, two questions decide whether the reported
differences are real:
A. Where the table says copy U uses a junction and copy S does not, do the splice-site sequences of the two copies
   point the same way?
C. Do the reads counted at a copy come from that copy?
A third count is descriptive: which assemblers emit the junction of a copy-specific call at the copy that uses it.

## Rule, fixed now
- **Samples.** The two gorillas develop, and their counts are below under "Seen already": OR6737 (`GGO_mm`), on which
  every script was debugged, and KB3781 (`GGO_KB3781_fib_mm`), opened at the user's request on 2026-11-18 to check a
  second sample of one genome before the freeze. No judging sample is on the gorilla genome. Judging samples, ten on
  five reference genomes, none opened for this question:

  | genome | samples |
  |---|---|
  | CHM13 | `A119b`, `HG002_kinnex.t2t` |
  | PTR | `PTR_mm`, `PTR_AG18354_lcl_mm` |
  | PPA | `PPA_mm`, `PPA_PR00251_fib_mm` |
  | PPY | `PPY_mm`, `PPY_AG05252_fib_mm` |
  | PAB | `PAB_mm`, `PAB_AG06213_fib_mm` |

  Each is the default assembly of 2026-11-16 (`artifact_refresh_2026-11-16/apes/SAMPLE/def/x.gtf`, A119b:
  `h.def/a1.gtf`) with its score table, linked into a new run directory; stages `families`, `exonmap`, `psi` are run
  on it with the driver's defaults as they stand: the two CHM13 samples take the genome's SEDEF table beside the
  FASTA (SD damping, the default since row 1580); the ape FASTAs used here have none, so their families are built
  without it. A sample whose run fails is reported as missing, not dropped silently.
- **Rows of a copy the stage dropped** (`psi_family.py`, part 1). The stage keeps an intron only with 3 reads and a
  second intron at one of its sites. For a copy that owns both exon columns of a key and has no row under it: its own
  introns between its exons of the two columns, with their real counts (also a junction alone at its sites, or shifted
  against the other copies'); when it has none, the junction of the key's most-read copy carried over through the
  pair's alignment (records agreeing with the two strands; both flanks in the copy's exons of the two columns), a
  silent row. A row needs 20 certain reads at its anchor site.
- **Calls** (`calls.py`). A key = (family, anchor, column of the anchor, column of the other end). Eligible copy =
  at least 20 certain reads on the introns leaving its anchor site, added rows included. Per key: U = the eligible copy
  with the highest `range_lo`, S = the other eligible copy with the lowest `range_hi`, gap = `range_lo`(U) -
  `range_hi`(S), rounded to 4 decimals. LARGE: gap > 0.2. SMALL: 0 < gap <= 0.2. SAME: ranges overlap and the eligible
  copies' `psi_placed_w` lie within 0.05. COPY-SPECIFIC: LARGE with `range_hi`(S) <= 0.05. Per (family, U's intron)
  the key with the larger gap is kept; one event read from both copies (inclusion key with U = a, skip key with U = b
  at the same anchor column) is one call: the orientation with the larger gap stays. MISKEYED: U's junction, carried
  to S, is an intron S uses (3 reads or more) that is not among S's introns under the key; the columns then keyed the
  two copies' matching junctions apart, and the call is counted and judged in no bar of A, C or the tool table.
- **Unit.** (reference genome, U's junction, S's contig, strand and anchor site; a copy's name is its locus span in one
  sample's assembly and is not used). A unit called in two samples of one genome counts once, with the
  numbers of the sample holding more certain reads at U. A unit that is a call on OR6737 counts in no bar and in no
  printed share.
- **A, site strength** (`sites.py`, `site_cause.py`). Position weights from the sample's own emitted introns (donor
  9-mer, acceptor 23-mer). Competitors of the key junction = the introns of at least 3 reads sharing its anchor site
  at either copy, each carried to the other copy through the pair's alignment. Advantage at a copy = strength of the
  key junction's site at its other end minus the strongest competitor's there; dA = advantage at U minus advantage at
  S. When a junction is carried by no read at one of the two copies, its site is read at both copies as the strongest
  within 6 bases (an image can slip beside a gap). A competitor without an
  image is dropped when it holds under 10% of its site's certain reads; otherwise, or when a site has no window (contig
  end, N), the unit is UNLOCATED. Verdicts at tau = 1 bit: AGREES (dA >= 1), OPPOSES (dA <= -1), NO SITE DIFFERENCE.
  Tau 0.5 and 2 are printed beside.
- **C, read origin** (`origin.py`). Sample: `HG002_kinnex.t2t`, the only judging sample with haplotype assemblies
  (KB3781, the other, is now development: its run is under "Seen already"). MISKEYED calls are left out. Truth: every read with a record on a member locus of a LARGE or SAME call, reads of every hash residue (user's
  ruling 1), no cap, scored by `bench/o2_truth/neutral_truth.py score` with `NT_REPAIR=1` against the two haplotype
  assemblies, then `tie_split.py`. A read counts at copy X when every record of its wide band (primary and records
  within 10 `ms` of its best) carries one intron of the call's bunch at X (the introns of at least 3 reads sharing the
  anchor site of X's key junction, and that junction). SUPPORTED: one origin copy, and it is X (same contig and strand,
  an exon base shared with X's column exons). AGAINST: one origin copy elsewhere, or several copies fit equally.
  OUTSIDE (the read's home is a copy the reference lacks): printed beside, in no denominator (user's ruling 2). Reads
  of any other truth kind, or without a truth row, are counted and judged nowhere.
- **Tool table** (`tool_table.py`). Units = COPY-SPECIFIC calls. Arms: Rustle `x.gtf` and `x.final.gtf`, StringTie,
  FLAIR, IsoQuant, bambu (all without a guide, the GTFs on disk from 2026-11-10 to 2026-11-16). Per arm: share of
  units with U's junction in an emitted transcript at U, and share with S's junction (its image, for a silent copy)
  emitted at S. StringTie and the others are trusted baselines.
- No annotation is read by any call or verdict. Intervals are Wilson 95%.

## Readings, fixed now
- **A1 (bar).** Pooled LARGE units with a site difference at tau 1: PASS when the share of AGREES among AGREES +
  OPPOSES has a lower bound above 0.5. Fewer than 30 such units: undecidable. Printed beside: the same at tau 0.5 and
  2, per genome, UNLOCATED, MISKEYED, units whose strongest site at S has lost the canonical dinucleotide, and the share of LARGE
  units held by the three largest families.
- **A2 (printed, no bar).** Share of located units with a site difference (AGREES or OPPOSES), LARGE against SAME, and
  the same split at pair identity 0.98.
- **C1 (bar, HG002).** Reads on the bunches of LARGE calls, each read once per copy and site: PASS when SUPPORTED over
  SUPPORTED + AGAINST is at least 0.95 with a lower bound of at least 0.90. Printed per truth kind, per copy role (U,
  S), per family; the OUTSIDE share beside; the same share on SAME calls beside.
- **C2 (bar, HG002).** A LARGE call is judged with at least 20 SUPPORTED reads in the bunch at U and at S. It HOLDS
  when the PSI gap computed from SUPPORTED reads alone (those on the copy's introns under the key, over those in the
  bunch) is at least half the stated gap. PASS when at least 0.90 of
  judged calls hold. Fewer than 20 judged calls: undecidable.
- **Step 0 of C.** The counts of judged reads and judged calls are read first. Under a floor the bar is undecidable
  and is written so; no second truth run.
- **Tool table (descriptive).** An arm is behind another when its share at U is more than 0.05 lower on at least 30
  shared units; otherwise level. Printed split by A's verdict.

## What each outcome will be said to mean
- A1 PASS: where two copies differ in the use of a junction and in the strength of its site, the two agree in
  direction more often than chance. Not: the site change causes the difference (see Limits).
- A1 FAIL or undecidable: the calls are quoted only with C's result, never as mechanisms.
- C1 and C2 PASS: the differences are differences between the copies' own reads, not misplaced reads.
- C1 or C2 FAIL: the table's differences at near-identical copies are not to be offered as findings; the result is
  recorded as negative.
- Tool table level: the new thing is the table of calls, not the assembly; that is what gets said. Rustle behind
  another arm: said as it is.

## Prediction
A1 PASS if 30 units are reached (development: 24 of 29). A2: the LARGE share above the SAME share. C1 PASS, and nearly by
construction: a read counts at a copy only when no tied record puts it elsewhere, and the primary alignment of such
reads was right for 0.999 in rows 1785 and 1796; the informative numbers of C are the OUTSIDE share and C2. C2 PASS if
20 judged calls are reached. Tool table: level between Rustle, StringTie and FLAIR; IsoQuant and bambu lower at U.

## Seen already (OR6737, development; counts in no bar)
Counts with the scripts as frozen (after the code review's fixes; job 4482532, `crosscopy_2026-11-18/dev/`):
- Stage: 3,652 keyed rows as before part 1 (byte-identical in the old 26 columns, `PREFIX.psi.tsv` byte-identical);
  716 added rows for 434 keys, 135 of them silent; 15 candidate copies not located, 70 with the image off their exons.
- Calls: 402 = LARGE 82 (copy-specific 30, silent S 15), SMALL 262, SAME 58; 39 families with a LARGE call.
- A at tau 1, LARGE: AGREES 24, OPPOSES 5, NO SITE DIFFERENCE 37, UNLOCATED 9, MISKEYED 7 (24 of 29 = 0.83, interval
  0.65 to 0.92). SAME: 3, 5, 48, 1, 1. Site difference among located units: LARGE 29 of 66, SAME 8 of 56.
- Tool table, 30 copy-specific units, emitted at U / at S: StringTie 30 / 3, Rustle 29 / 7, FLAIR 28 / 7, Rustle final
  28 / 5, IsoQuant 27 / 3, bambu 20 / 0.
- C, mechanics only: OR6737 has no haplotype assembly, so its reads at the calls' loci (59,384) were scored against the
  haplotypes of the other gorilla, KB3781, to see the job run end to end. Not a measurement of OR6737. Reads on the
  bunches of 75 LARGE calls (7 MISKEYED left out): 16,113 SUPPORTED of 16,114 judged (kind `several`: none), 1,380
  OUTSIDE, 9,371 without a verdict; 44 of 44 judged calls HOLD, 31 unjudged.
- What development changed in the rules, all before any judging sample was run: the 6-base search (two agreeing calls
  were a 2-base slip of a carried-over junction onto an intact site); then, from a fresh reviewer's reading of the
  code, the rows of a copy's own junctions (13 of 36 calls with an added S were copies using the junction a few bases
  away, written as silent), the unit, MISKEYED, the two-sided search, and the origin PSI's numerator. The earlier
  development counts (89 LARGE; AGREES 29, OPPOSES 5) are kept in `dev/v1/`.

### KB3781, second development sample (opened 2026-11-18 at the user's request; counts in no bar)
Run with the same scripts, before the freeze (`crosscopy_2026-11-18/dev_kb/`, jobs 4482577 and 4482629;
`crosscopy_2026-11-18/replicate.py` for the check across the two samples).
- Calls: 515 = LARGE 111 (copy-specific 34, silent S 16), SMALL 321, SAME 83. A at tau 1, LARGE: AGREES 27, OPPOSES 7,
  NO SITE DIFFERENCE 61, UNLOCATED 13, MISKEYED 3. Without the units OR6737 also calls: AGREES 20 of 27 (0.74, interval
  0.55 to 0.87), under 30 units; site difference LARGE 27 of 80, SAME 13 of 64.
- Tool table, 33 copy-specific units, emitted at U / at S: Rustle 32 / 4, StringTie 32 / 3, FLAIR 31 / 9, Rustle final
  30 / 4, IsoQuant 25 / 1, bambu 23 / 0. Over all LARGE calls of either gorilla no junction is emitted by Rustle alone.
- A sample's calls in the other sample's reads (PSI recomputed from the raw intron counts at the same coordinates, at
  least 20 certain reads at both copies). OR6737's LARGE calls in KB3781: 48 of 75 have the reads; same direction 47,
  at least half the stated gap 38, gap over 0.2 again 37. KB3781's in OR6737: 66 of 108; 66, 61, 58. SAME calls stay
  within 0.1: 34 of 34 and 51 of 55. SMALL calls: same direction 97 of 142 and 133 of 225, so a SMALL call is not
  a finding; no bar reads them.
- Units across the two samples: 66 shared of 402 and 515; class agrees in 47 (LARGE in both 13, never LARGE against
  SAME); the 13 have the same A verdict in both. Of OR6737's 82 LARGE calls, KB3781 has the same unit for 16, the same
  junction of U against another copy for 10, too few reads for 29, and no call at those anchor sites despite the reads
  for 19. Three of the 19 have both copies under one key in KB3781 and were looked at: in two KB3781 does make the
  LARGE call, on the same junction keyed by its other end and between another pair of copies (a key gives one call,
  between its two extreme copies, and KB3781 has a further copy with enough reads); in the third the difference is
  not there in KB3781 (the less-using copy shows the junction in 22 of 24 reads against at most 0.73 in OR6737).
  Matched by junction instead of by unit, of OR6737's 82 LARGE calls KB3781 calls the same two copies' junctions in
  19 (15 LARGE), one of the two against another copy in 16, neither in 47 (28 of them without the reads). So the
  table's reach depends on each sample's assembly and copies, and absence of a call in a sample is not absence of the
  difference.
- SMALL calls by depth and by stated gap (`replicate_small.py`): more reads do not make them repeat (same direction
  0.61 to 0.78 across depth classes one way, 0.57 to 0.64 the other, no trend up to 300 reads and more at both copies
  in both samples). The stated gap does: 0.1 to 0.2 repeats in direction in 23 of 26 and 29 of 36; under 0.02 in 30 of
  50 and 44 of 100.
- C on KB3781's own haplotypes, the first real run of the scorer: reads on the bunches of 108 LARGE calls: 63,634
  SUPPORTED of 63,636 judged (kind `several`: none), 11,691 OUTSIDE (0.155 of judged + outside), 18,096 without a
  verdict; 82 of 83 judged calls HOLD (copy-specific 29 of 29), 25 unjudged.
- No rule changed after this sample.

## Limits stated beforehand
- A1 cannot separate cause from consequence: U's junction is by construction one that reads use, so a changed site at
  S is more likely a weaker one. What A1 excludes is the mapping artefact, where S really uses the junction and both
  sites are intact.
- Units are not independent: several junctions of one pair of copies, and the same ancestral difference in two
  species, count as separate units. The interval is narrower than the evidence; the per-genome table is printed for
  that reason.
- Only junctions between exons of a family's representative transcripts are keyed; others are competitors, never calls.
- A copy's introns between the same two exon columns are all its use of the key, whatever their shift: a splice site
  moved inside one exon pair is not a call here (row 1714's comparison shows those). A shifted site that falls outside
  a copy's representative exon is keyed at one copy and not the other; MISKEYED catches the cases where the other copy
  uses the same junction, not all.
- C1 is a weak bar (see Prediction): it can fail only if clearly placed reads come from another copy. Reads tied
  between copies are in no count of the table and in no count of C.
- C judges two samples, one of them on its own reference. A's samples are all primates; nothing here says the table
  works on other clades.
- A unit keeps S in its identity, so one junction of U called against two copies, or against different copies in two
  samples of a genome, is two units (10 of OR6737's 82 LARGE calls recur in KB3781 that way). Bar A1 counts units as
  drafted (user, 2026-11-18). Printed beside it, in no bar: the same share with linked units counted once (LARGE units
  of a genome that share U's junction or S's anchor site are one event; the unit with the largest gap speaks). On the
  two gorillas pooled with no exclusion: 175 LARGE units are 161 events, AGREES 41 of 51 against 43 of 54 by units.
- A119b and OR6737 were development samples of the assembler's defaults. No rule of this document was set on A119b.
- The position weights come from each sample's own introns, so strengths are comparable within a sample, not across.

## Amendments

## Result
