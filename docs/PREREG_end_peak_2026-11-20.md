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
