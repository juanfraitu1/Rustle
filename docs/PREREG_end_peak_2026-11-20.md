# Row 1840: keep a contained chain when its own end is a read-end peak (pre-registered 2026-11-20)

User (2026-11-20): pre-register and build. Rows 1838 / 1839: the read ratio of the sub-chain collapse cannot tell a true
contained chain from a truncated copy of the longer one (ratio off: mouse +0.7-1.9% correct chains at 2-3 wrong per correct,
apes +0.1-0.3% at 16-60 wrong). IsoQuant decides such a chain by its end: it stays when its own terminus is a confirmed peak
of read ends. This row builds that criterion as an opt-in switch beside the ratio. Workspace `r1840_end_peak`. Committed before
the code was written and before any run.

## The rule, fixed
Switch `copy_assign --polish-ism-end-peak N` (`RUSTLE_POLISH_ISM_END_PEAK=N` in the driver; 0 = off = default; N = least
number of reads at the peak). Evidence = the deduplicated primary spliced records the polish already collects for the 5' / 3'
end layers; a transcript's OWN reads = the records whose intron chain equals its chain.
- A multi-exon transcript Y is *contained* when its chain is a contiguous part of a longer chain X of the polish input (same
  strand). Its 5' side is *internal* when some container has introns upstream of Y's first, its 3' side likewise downstream.
- **5' peak:** some window of 21 bp (the layers' tolerance, +-10) holds the 5' ends of at least N own reads and of at least
  half of the own reads.
- **3' peak:** a cluster of own 3' ends (the 3' layer's clustering) with at least N reads and at least half of the own reads,
  with a polyadenylation signal 10-35 bp upstream and no internal-priming signature (the layer's existing proof).
- Y is *end-confirmed* when every internal side has its peak. The sub-chain collapse then keeps Y whatever its read ratio.
  Every other polish step still judges Y. Nothing else changes; without the switch the output is byte-identical.
Constants chosen a priori, not fitted: half of the own reads; N = 2 (the read floor) and N = 3 as the two arms.

## Arms and data
`base` (ratio 0.3), `peak2`, `peak3` on the twelve datasets of rows 1837 / 1839; `r0` of row 1839 (the rule off) is the
ceiling and is not rerun. Stage `assemble`, `RUSTLE_PTC=0`, scorer of row 1839. Reused datasets, said.

## Reading, fixed now (`x.final.gtf`; `x.gtf` beside it)
- **Pays** on a dataset as in row 1839: the marginal precision of the chains the arm adds is at least the dataset's `base`
  precision, or precision falls by at most 0.005 while correct chains rise by at least 1%.
- **PASS** = an arm pays on at least 6 of 7 non-primate and 4 of 5 primate datasets. A PASS makes the switch a candidate; the
  user decides about a default.
- **Separation** (the mechanism claim, `x.gtf`): among the chains ratio 0 adds, precision of those the arm admits against
  those it does not. The criterion SEPARATES on a dataset when the admitted are at least twice as precise. Reported per dataset.
- Share of the ceiling: correct chains the arm adds / correct chains ratio 0 adds.
- Parity: without the switch the new binary's `x.gtf` and `x.final.gtf` are byte-identical to row 1839's `base` (md5).

## Expectation, written before the build
Degraded molecules are cut at the 5' side, so most contained chains are 5'-internal with scattered starts and fail the test;
I expect the arms to admit a small part of the ceiling (under a third of its correct chains) at clearly better precision than
the rest on the mouse libraries. On the ape testis libraries I expect little to admit and am unsure it separates. N = 2 may be
too weak: two truncated reads can start within 21 bp by chance.

## Limits
- The polish input comes after the junction snap; a transcript whose junction was re-placed has no own reads under its new
  chain and cannot be confirmed.
- A side counts as internal against any container, also one the polish later drops: conservative.
- Simulated reads carry no degradation and an idealised 3' end, so they flatter the test.
- Annotated ends on real libraries are not used: the bar is chain correctness, as in row 1839.

## Result (2026-11-20): FAIL; the end peak does not separate true contained chains from truncated ones

Jobs 4489469-4489480, 36 arms exit 0; `r1840_end_peak/sep.NAME.txt` (`sep.py`). Code `5e367e57`, binaries in
`~/ws/rustle_target_b`. Parity: `base` of the new binary is byte-identical to row 1839's `base` (md5 of `x.gtf` and
`x.final.gtf`, checked on fly, mouse ES, A119b).

**Pays (`x.final.gtf`):** `peak2` on fly, mouse and zebrafish simulated and on A119b: 3 of 7 non-primate, 1 of 5 primate.
`peak3` on the three simulated sets only: 3 of 7, 0 of 5. Bars were 6 and 4.

**Separation (`x.gtf`, among the chains ratio 0 adds: precision of those the arm admits / of the rest), `peak2`:**
| fly sim | mouse sim | zebrafish sim | mouse ES | heart | adrenal | cortex | OR6737 | A119b | chimp | bonobo | KB3781 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.105 / 0.696 | 0.301 / 0.850 | 0.066 / 0.381 | 0.051 / 0.053 | 0.087 / 0.085 | 0.081 / 0.095 | 0.071 / 0.091 | 0.005 / 0.006 | 0.008 / 0.006 | 0.005 / 0.006 | 0.008 / 0.009 | 0.005 / 0.006 |
SEPARATES on 0 of 12 (`peak3`: 0 of 12). On the real libraries the admitted chains are as precise as the ones left out. On the
three simulated sets they are several times LESS precise: there the test picks the wrong contained chains.

`x.final.gtf`, `peak2` against `base` (correct added / wrong added; precision):
fly +4 / 0 (0.964); mouse simulated +166 / +2 (0.954, 12 correct lost); zebrafish +6 / 0; mouse ES +160 / +536 (0.347 ->
0.346); heart +49 / +82; adrenal +80 / +124; cortex +69 / +135; OR6737 +17 / +532; A119b +88 / +734 (21 correct, 508 wrong
lost through the filter); chimpanzee +10 / +243; bonobo +8 / +268; KB3781 +14 / +592. After the filter the admitted chains
are somewhat more precise than the rest on mouse adrenal (0.40 against 0.22) and A119b (0.11 against 0.04), nowhere else twice.

- The test confirms a quarter to a half of all contained transcripts (A119b 68,449 of 276,362; mouse ES 5,690 of 10,397; mouse
  simulated 32,011 of 69,212) and admits 30-59% of what ratio 0 admits. It is not selective.
- Expectation: wrong. I expected a small admitted share at clearly better precision on the mouse libraries.
- Reading (mine): the starts of truncated copies are not scattered. Two or three reads of a 5'-cut chain start within 21 bp
  of each other about as often as the reads of a true shorter form do, so a peak among a chain's OWN reads carries no
  information. IsoQuant's test is not this one: it asks for a detected polyA tail or a start peak judged against all reads
  crossing the position. That stronger form (the peak against the reads running through, which is what the cap-gated TSS
  proof already does at the 5' side) was not built here.

**Decision:** not a candidate. The switch stays in the code as an opt-in (off by default, parity checked); removing it is the
user's call. The shared `~/rustle_target` was not rebuilt with it.
