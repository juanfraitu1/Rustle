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
