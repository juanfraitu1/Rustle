# PREREG r1822: borrowed structure (question 4), shared isoforms and RUSTLE_READ_ISOFORM_K=2 on KB3781 (2026-11-19)

Bars K1, B1 accepted by the user on 2026-11-19 ("I accept the bars"), before any arm below was run on KB3781. Workspace
`~/ws/r1822_borrow_isoformk` (PSU; register tail on origin/main was 1821). Same message: the user asked to also take into
account borrowing as ISOFORMS THAT SEVERAL FAMILY MEMBERS PRODUCE: added as arm F, descriptive, no bar (none was shown to
the user). Goal:
`docs/GOAL_multicopy_families_2026-11-19.md` (questions 4 and the member gap). Development: workspace
`goal_multicopy_dev/LEDGER.md` steps 9-11 (seen HG002 quarter).

## Why KB3781
HG002 has no unseen quarter left (rows 1796, 1820, 1821). KB3781 is read on its own reference against its own haplotypes
(truth `r1785_truth_repair/kb`, tie split as in row 1820's development run). It was a development sample for rows 1796 and
1820 (mechanics only; none of the arms below was looked at on it). Reuse per the user's ruling (used datasets are reused).

## Arms (frozen code of commit ae78e8f5)
P (`--seed-pool primary`); P_k2 (P with `RUSTLE_READ_ISOFORM_K=2`); P_borrow (P plus `borrow.py` chains, labelled
`borrowed_structure`); StringTie `-L` on primaries.
F (shared isoforms, `bench/famwise_eval/family_isoform.py`): P1 = P with `RUSTLE_MIN_CHAIN_READS=1` (one-read chains);
P_famiso = P + the P1 chains absent from P whose spliced sequence matches an isoform P assembles at another place end to end
(both >= 90% covered, no internal gap >= 20 bp; minimap2 asm20 cDNA to cDNA), labelled `family_isoform`. Control: P1 itself
(every one-read chain). Printed: P's isoforms shared with another place; promoted chains' verdicts and new / lost against P;
members found.

## Measures
`copyset_score.py --match contained` (precision, weighted share, resolved units, new / lost against P);
`member_completeness.py --species gorilla` (members found of members with >= 2 own primaries; class of each miss). The
borrowed chains are scored apart (as structure: units whose chain they hold; as expression: their EQUAL share), never
pooled into P's precision.

## Truth size on PSU (disclosed before the run)
`famtruth` gorilla families on PSU: 70 families, 152 members (gorilla genes reach Compara through Ensembl IDs; Fritz adds the
row 1611 gene -> human Ensembl map, `DENOVO_GGO_ENS_MAP`, which is not on PSU). K1's member counts are read on these.

## Readings (accepted 2026-11-19)
- **K1 (the member gap).** P_k2 finds more expressed members than P (>= +2% of P's), with one-copy precision lower bound
  >= P's point estimate - 0.015. Development (HG002 seen quarter): 277 vs 268 (+3.4%), 0.884 vs 0.896.
- **B1 (borrowed structure adds structure).** P_borrow holds the chain at every copy for >= 3 more copy-set units than P.
  Development: 68 vs 63.
- **B2 (it is not expression).** Printed: the borrowed chains' EQUAL share as expression claims (development 0.17), and the
  share of them with no judged supporter. No bar: B2 says why the class stays labelled.

## Prediction (mine)
K1 passes; B1 is undecidable or fails on KB3781 (few tied reads: 186 exact ms ties of 100,810 molecules).
