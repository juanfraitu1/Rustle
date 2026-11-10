# PREREG r1604: tie groups with near ties on AS, and on AS or ms (2026-11-10)

Follows r1603 (`docs/PREREG_tie_groups_2026-11-10.md`, FAIL: reach 39/72, purity 0.858, TAF11L 11/31). Post hoc there: exact AS
equality keeps 1-2 molecules per TAF11L unit. User (2026-11-10): ties are read on AS and on ms (minimap2's ranking score), so
the tie set may need to be bigger. Written and committed before any arm runs.

## Arms (only the tied-record rule changes; loci, groups, models, panel and M as in r1603)
- **as98**: a record is tied when AS >= 0.98 x the molecule's best AS (the r1521 panel's T rule).
- **asms98**: as98, or ms >= 0.98 x the molecule's best ms. A superset of as98.
- **exact** (r1603) is rebuilt from the saved records as the control; it must reproduce 39 / 0.858 / 11.
A molecule is tied when >= 2 of its records are tied under the arm's rule.

## Run
The molecules' records are saved once (`samtools view -N names.txt -b`, r1603's `names.txt`), so each arm is a build from that
file, not a BAM pass.

## Bars (per arm, as r1603)
- B1 reach >= 58 of 72; B2 purity >= 0.90; B3 TAF11L models <= 5 kb reach >= 20 of 31; B4 the whole job <= 1 h wall.
- An arm PASSES on all four. If both pass, asms98 is kept only if its purity is not below as98's; otherwise as98.
- Any miss in an arm = that arm FAILS as registered, with the numbers.
- Descriptive: tied records, loci, groups, single-locus groups, loci with no CAT gene, per arm.

## Outcome (2026-11-10, job 4456250, `/anvme/workspace/v138da13-rustle/r1603/{exact,as98,asms98}.*`): both arms FAIL as registered, on purity only
| arm | tied records | loci / groups | B1 reach (>= 58) | B2 purity (>= 0.90) | B3 TAF11L (>= 20) |
|---|---|---|---|---|---|
| exact (control) | 76,500 | 2,658 / 1,615 | 39 | 273/318 = 0.858 | 11 |
| as98 | 405,324 | 5,541 / 2,599 | **72** | 188/234 = **0.803** | **31** |
| asms98 | 493,065 | 5,854 / 2,881 | **70** | 181/222 = **0.815** | **31** |

- B4: 2,077 s for the whole job (records 406 s; three builds in parallel) - pass. The control reproduces r1603.
- Near ties fix reach and the TAF11L unit models completely; purity falls further below the bar.
- asms98 unreached: TSPY3, TSPY4. Loci with no CAT gene in judged groups: 1,328 / 1,438 (247 at exact).
- Other 78 hard members reached: 29 in both arms (16 at exact).

## Addendum (2026-11-10, user: older Compara paralogues as an extra "or" category; written before it is computed)
Descriptive split of the judged groups, per arm; the registered B2 verdict above does not change.
- **pure**: as B2.
- **older paralogue**: every foreign locus overlaps a CAT gene that Compara lists as a paralogue (any last common ancestor,
  `DENOVO_PARALOGS`) of a gene of the group's top panel family.
- **unrelated**: at least one foreign locus has no such gene.
`tie_groups.py score` prints the split when `DENOVO_PARALOGS` and `DENOVO_CAT_GFF` are set.

### Split (2026-11-10)
| arm | pure | older paralogue | unrelated | pure or older |
|---|---|---|---|---|
| exact | 273 | 5 | 40 | 0.874 |
| as98 | 188 | 4 | 42 | 0.821 |
| asms98 | 181 | 3 | 38 | 0.829 |

Older Compara paralogues explain 3-5 groups per arm; the bar stays missed with them counted.

### Post hoc reading of the 46 impure as98 groups (not a registered result)
- The "foreign" genes are mostly not in the panel at all: lncRNA 43, unprocessed pseudogene 40, transcribed unprocessed
  pseudogene 36, protein-coding outside the panel 17, TEC 12; protein-coding of another panel family 30.
- Three kinds: (a) pseudogene or unexpressed copies of the same family, which the scorer cannot see because its family table
  lists expressed members only and Compara has no pseudogenes (TSPY13P/6P/7P, AGAP10P..., AMYP1, CBWD4P + CBWD1, RBMY1HP,
  CNTNAP3P2, NUDT4P2, OR11H13P, DHFRP1, SPATA31D5P, PLGLB1, CTAGE4, SPDYE18, NOTCH2 for NOTCH2NLR); (b) co-duplicated
  neighbours in the same segmental duplication (NPIPB with BOLA2 / EIF3C / SMG1P, GTF2H2 with SMN2 / NAIPP1, RASA4 with
  POLR2J3 / POM121B, NPIPA with NOMO1): the tie group is a duplication-level unit, not a gene family; (c) one 831-locus
  group around POTE, and groups of 86-110 loci (CLEC18A, NBPF, TAF11L): single-linkage chaining.
- Scorer limits to fix in a new row before purity is judged again: family membership from the full truth (not the expressed
  panel rows) and a rule for same-family pseudogenes.
