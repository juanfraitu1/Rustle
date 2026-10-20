# Pre-registration: O2 as a corrector of minimap2 on held-out reads (2026-10-20)

**User (2026-10-20):** "if minimap can already find the best alignment there is no need for this thesis. The entire point is
to use the tied and nearly tied AS alignments and find where they truly belong as an addendum to minimap2 since it gave
everything a coin toss. It is fine the ones that are more than >99% similar are impossible but there must be a way with
PSVs, de, nm, md tags to find the true placement of the multimappers".

So O2 is judged as a **corrector**, not a certifier. For every read whose alignments within 10 AS of its best lie on two
or more copies of one family (Eichler), O2 gives a copy. It is scored against minimap2's primary on the same reads:
- **fixes:** O2 right where minimap2 was wrong;
- **breaks:** O2 wrong where minimap2 was right;
- **accuracy:** of both, over the same reads.

Version 2 (registers 1341-1342) only certified a subset of minimap2's calls and never overruled one.

## The corrector (`bench/o3_graph/o2correct.py donor`, frozen at the commit that adds this file)

**Domain.** A read whose alignments within 10 AS of its best lie on two or more copies of one family and on no locus off it,
minimap2's primary among them.

**The sample's copies.** Per copy, the sample's confident reads on it: the primary is on the copy and every alignment on
another copy of its family is more than 9 AS lower. From them, each copy's **events**: substitutions, and indels shifted
left as far as the reference repeats them.

**Scoring.** All tied alignments of a read are compared over the read bases every one of them aligns. A tied copy scores:
- the read's events its own reads do not explain. An event is explained when at least 2 of them carry it, making up at least
  10% of their depth, at a depth of at least 3;
- plus the events at least 80% of its reads carry that the read lacks.

**Decision.** The lowest score takes the read. It overrules minimap2's primary only if it is still lower after the primary
is granted every event its reads do not cover (dominance). A shared lowest keeps minimap2's pick. A fragment's own source
read is left out of the sample that judges it.

**In one genome** a copy's reads explain only their shared errors, so the rule is minimap2's ordering by events.

## Development: every number seen before this file

- **HG002** (`hg002proof`, donor-assembly truth, PAR out). The rule and its constants were chosen here:
  - substitutions only (the MD/NM rule);
  - then the sample-aware rule, with HET 0.2 / 0.1 / 0.05 crossed with dominance on/off.

  | HG002, chosen rule | reads | corrector | minimap2 | fixes / breaks |
  |---|---|---|---|---|
  | near ties | 5,488 | 87.3% | 66.4% | 1,161 / 15 |
  | exact ties | 576 | 52.6% | 39.6% | 97 / 22 |

  15 families gain and 4 lose. The substitution-only rule was right in two families (SMN, THOC3) and a coin flip elsewhere.
- **Simulation, heldout23** (one genome). Breaks in CTRB2 led to comparing alignments over common read bases. After that,
  near ties were 2,016 reads with 0 fixes and 0 breaks.
- **Simulation, heldout21** (donor alleles, PAR out):

  | heldout21 | exact ties fixes / breaks | near ties fixes / breaks |
  |---|---|---|
  | worlds without BRANCH | 31 / 0 | 0 / 0 |
  | BRANCH worlds | 38 / 35 | 0 / 7 |

  A BRANCH world's reference lacks one of the donor's copies (O3's case). Its reads, confidently placed on the nearest
  reference copy, teach that copy alleles it does not carry.
- **A119b fragments, seed 1** (`o2frag_a119b`, qualifying non-PAR families). This set also exercised the library pass.
  1,518 fragments: corrector 80.0% vs minimap2 72.1%, 125 fixes / 5 breaks; 45 families gain, 2 lose.

## Held-out sets: none has been run through the corrector

| set | what | truth |
|---|---|---|
| **H1** | A119b fragments, **seed 2** (`o2frag_a119b2`: 11,280 sources, 33,840 fragments). Used only by version 2 in register 1342, never by this rule. Held out by reads, not by donor. | the source read's copy, clear by 30 AS |
| **H2** | public testis fragments (`o2frag_public`, library ERR13885926): **another donor**, shallow | the same |
| **S1** | simulation, one genome: heldout24 + heldout25 (36 worlds: CROSS, INTACT, SILENT) | `reads.truth.E0.tsv` |
| **S2** | simulation, donor alleles: heldout22 (36 worlds: BRANCH, CROSS, INTACT, OUTSIDE, REMOVED, SILENT) | the same |

**Common scoring rules:**
- A read is scored only when its truth is one of its tied copies (O2 cannot pick an absent copy).
- PAR families are out (`famident/par_families.json`).
- Fragments are counted one by one; their source reads are counted for power.
- Real sets are scored on the qualifying non-PAR families: every member pair at [0.90, 1.00) identity,
  `famident/qualify_nonpar.json`.
- Commands: `o2correct.py score OUT truth.tsv --among-candidates --families qualify_nonpar.json` and
  `o2correct.py score-worlds OUT.json WORLDS... --exclude par_families.json`.

## Bars

**H1 and H2, each, pooled over exact and near ties:**
- **P (power):** at least 50 overruled fragments, from at least 20 source reads.
- **G (gain):** corrector accuracy at least 2 points above minimap2's.
- **S (safety):** fixes at least 5 times breaks.
- **F (families):** at least 11 qualifying families where it fixes more than it breaks (user: "more than 10 families"). No
  more than a fifth as many families, rounded up, where it breaks more than it fixes.

**S1, without BRANCH worlds (none exist in S1):**
- near ties: corrector accuracy no more than 0.2 points below minimap2's, and breaks no more than fixes + 2;
- exact ties: corrector accuracy no more than 0.5 points below minimap2's.

**S2, the worlds without BRANCH:**
- corrector accuracy at least minimap2's in each band;
- breaks at most fixes / 5 + 2.

BRANCH worlds are reported apart. That is O3's case, and development already showed its harm.

## Verdict

- **PASS:** H1 passes P, G, S and F, and S1 and S2 pass.
- **FAIL:** H1 has P and fails G, S or F; or S1 or S2 fails.
- **UNDERPOWERED:** H1 fails P.
- **H2** gets its own PASS / FAIL / UNDERPOWERED on the same bars. A FAIL there limits the claim to A119b's donor.

**A PASS would claim:** among a read's tied and near-tied copies, the sample's own reads of each copy overrule minimap2's coin
toss far more often right than wrong. This holds on reads the rule never saw, in more than 10 near-identical families, and
does no harm in one genome.

**It would not claim:**
- accuracy against a donor assembly on a new donor (H1 and H2 have fragment truth);
- copies absent from the reference (O3);
- PAR;
- exact ties over identical sequence.

**Known bias:** a fragment's truth is minimap2's call on its full read. A systematic full-read error therefore counts against
the corrector, never for it.

**Rules:**
- Run script: `/mnt/linuxdisk/tmp/o3_graph/o2correct_heldout/run.sh`. The library pass resumes in foreground batches (`--cache`,
  `--budget`).
- A bug found during the run is fixed and recorded as an amendment before any rerun.
- No constant changes after this file.

## Results (`/mnt/linuxdisk/tmp/o3_graph/o2correct_heldout/run1.log`, `run2.log`; frozen copy `frozen/`)

**Overall: FAIL. S2 fails. H1 and S1 pass. H2 is UNDERPOWERED.**

| set | reads scored | corrector | minimap2 | fixes / breaks | bars |
|---|---|---|---|---|---|
| **H1** A119b seed 2, qualifying non-PAR | 1,241 fragments | **0.8171** | 0.7147 | **130 / 3** | P ✓ (133 overruled from 118 sources), G ✓ (+10.2 points), S ✓, F ✓ (40 families better, 1 worse; the limit was 8) |
| H1 exact ties | 440 | 0.6091 | 0.4523 | 69 / 0 | — |
| H1 near ties | 801 | 0.9313 | 0.8589 | 61 / 3 | — |
| H1, every non-PAR family | 1,727 | 0.7823 | 0.6775 | 186 / 5 | 55 families better, 2 worse |
| **H2** public testis, qualifying | 276 | 0.7826 | 0.7464 | 11 / 1 | **P fails** (12 overruled, from 12 sources): UNDERPOWERED |
| **S1** one genome (36 worlds), near ties | 2,309 | 0.9991 | 0.9991 | 0 / 0 | ✓ |
| S1 exact ties | 2,824 | 0.4770 | 0.4763 | 2 / 0 | ✓ |
| **S2** donor alleles, without BRANCH (30 worlds), near ties | 3,447 | 0.9985 | **1.0000** | **0 / 5** | **✗** (accuracy below minimap2's; 5 breaks > 14 / 5 + 2) |
| S2 exact ties | 3,079 | 0.4303 | 0.4258 | 14 / 0 | ✓ |
| S2 BRANCH (reported apart) | 786 exact, 812 near | — | — | 4 / 1 exact, 0 / 2 near | — |

**The S2 failure is the BRANCH mechanism, which development showed.** All 5 breaks are in one REMOVED world,
`A119b_s1201`, family ANKRD30BL.
- The world masks copy .3 out of the reference while the donor still expresses it: 111 reads.
- Those reads sit confidently on copy .4 and make up most of its sample.
- True .4 reads then look foreign on .4, and the corrector moves them to .0 (scores .0:4 vs .4:5).
- 10 more overrules move .3 reads themselves; they are not scored.

BRANCH breaks the same way in 3 of its 6 worlds. Copies the reference lacks are O3's case; O2's corrector is not safe on a
family until O3 has flagged them.

**Erratum in the development section.** "Worlds without BRANCH" in heldout21 was 24 worlds, not 30. All 6 REMOVED runs
crashed: their `genome.ref.fa` points to a deleted masked genome. `score-worlds` dropped them without saying so; it now lists
skipped worlds. In heldout22 all 36 worlds ran.

## Audit (metric-auditor, 2026-10-20): holds with caveats

Every number recomputes exactly from the frozen files with an independent script, and the overall FAIL stands. The audit also
confirmed:
- **Hold-out:** 0 of 15,573 seed-1 sources are reused. The prereg was committed at 12:50:27, before the first held-out output
  at 12:50:36.
- **Leave-one-out:** 1,329 of the 1,349 domain sources were confident reads, and each was left out of its own fragments'
  judgement.

Caveats:
1. **H1's truth and the corrector's teachers come from the same library's confident calls.** When the donor has a copy the
   reference lacks, both are wrong the same way, and that is scored as a fix. So this file's "never for it" is false in that
   case. H1 "right" means "agrees with minimap2's full-read call". The top fix families (NBPF6, TBC1D3B, FAM205A) vary in
   copy number between people. Only donor-assembly truth settles this; HG002, the one donor with it, was development.
2. **H1 replicates development rather than testing a new donor.** It shares seed 1's library, its 161,037 teacher reads and
   its families; only the source reads are disjoint. H2, the other donor, is underpowered.
3. **Fixes are spread out.** 130 fixes come from 115 sources. The top family has 10.8% of them, the top 10 have 59.2%, and
   17 families have at least 3.
4. **Filters do not decide a bar.** 14 non-PAR overrules have no scored truth. S still holds if all 14 were breaks
   (130 vs 17). G is still +10.0 points if the 33 rows whose truth is not among the tied copies count against it.
5. Excluding REMOVED from S2 after the fact (near 0/0, exact 11/0) would be choosing which worlds to score. It is not a
   result.

**What this shows.** Reading PSVs against the sample's own copies, learned from its confidently placed reads, corrects
minimap2 among tied copies:
- HG002 development, donor truth: near ties 66% → 87%;
- fresh A119b fragments: 72% → 82%, in 40 families;
- one genome: nothing to correct and nothing broken.

It is **not adopted**. The pre-registered donor simulation shows the corrector moves reads wrongly whenever the reference lacks
one of the donor's copies (BRANCH, REMOVED), a case O3 must flag first. A next version needs that guard and fresh held-out
data, ideally a second donor with an assembly.
