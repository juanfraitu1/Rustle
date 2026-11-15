# PREREG row 1715 (2026-11-14): why is the nested isoform lost? Measurement only

User go 2026-11-14 ("yes lets do that") after row 1714: nested alternative starts / early ends are emitted as both
isoforms in about half the expressed pairs. Seen before writing: row 1714's tables, and its pairs re-split by chain
length (no new read data): at minor >= 5 reads the lost pairs are "long only" (OR6737 NEST5 116 of 122, NEST3 91 of
95; A119b 539 of 589 and 702 of 816); "short only" 5 / 4 / 34 / 43; neither 1 / 0 / 16 / 71.

**Disclosed weakness of row 1714 on these classes:** a read counts for the short chain when its intron chain equals
it, so 5'-truncated (or early-ending) reads of the long isoform count as full-chain reads of the short one. The
"expressed" denominator of NEST5 / NEST3 is inflated by that; the short isoform has to be judged on read ends.

## Measured (`bench/famcmp/nested_sub.py`; OR6737, A119b; KB3781 ours only)
CANDIDATE = a read chain with >= 5 primary reads that is a terminal contiguous sub-run of a chain emitted in
`PREFIX.gtf` (container K = the emitted one with most reads). Class by the container's strand: NEST5 (shares K's 3'
end) or NEST3 (shares K's 5' end). Per candidate: emitted in `PREFIX.gtf`; annotated (RefSeq chain, same strand); in
the 1-read-floor run's GTF (rows 1616 / 1617); and on its variable end:
- `inside`: >= 0.9 of its reads end inside K's exon there (not beyond that exon's far boundary);
- `peaked`: >= 0.5 of its reads end within 10 bp of the commonest end;
- NEST5 `own`: fewer reads continue 5' past the chain than spell it exactly (the opposite of the `trunc5` rule);
- NEST3 `polyA`: polyA signal upstream of the commonest end and no A-rich run after it (`alt_ends.flags`);
- `pooled`: K's `reads` attribute minus K's exact-chain reads (does K carry the candidate's reads?).
Tables per class: annotated x emitted counts; each feature's share among annotated-kept, annotated-lost and
unannotated-not-emitted candidates; for the 8 signatures (inside x peaked x own | polyA), not-emitted candidates and
their annotated share. Also counted from row 1714's pairs: the "short only" and "neither" pairs whose long chain is in
the 1-read-floor GTF.

## Reading rule (fixed here)
A signature is worth a rule row only when, on OR6737 and A119b both, it holds >= 30% of the annotated-lost candidates
of its class and >= 0.30 of the not-emitted candidates with that signature are annotated (the 5' truncation class is
0.04-0.08 correct, row 1677). Otherwise: read ends do not separate a real nested isoform from a truncation, and no
lever is proposed. Any rule is a separate row on the user's go. The mechanism in the code is reported separately from
a code reading, not inferred from these tables.

## Outcome (2026-11-14, job 4471356; `r1715_nested/{or,a}.{txt,tsv}`): no read-end lever; the default rule stays
Candidates (read chains with >= 5 reads, terminal sub-run of an emitted chain), OR6737 / A119b:

| class | emitted: annotated / unannotated | not emitted: annotated / unannotated | annotated share of the not emitted |
|---|---|---|---|
| NEST5 | 471 / 3,272 and 1,666 / 19,625 | 133 / 20,407 and 427 / 77,756 | 0.006 / 0.005 |
| NEST3 | 859 / 1,336 and 2,514 / 10,446 | 111 / 2,911 and 590 / 21,269 | 0.037 / 0.027 |

Reading rule: no signature reaches 0.30 annotated among its not-emitted candidates on either sample. Best: NEST3
"ends beyond the container's exon, peaked, polyA signal" holds 60% / 47% of the annotated-lost but is 0.081 / 0.098
annotated; every NEST5 signature is <= 0.012 except cells of under 100 candidates. **No rule row is proposed.**
- Features do not separate lost from kept annotated chains (inside 0.14 vs 0.15, peaked 0.28 vs 0.30 for NEST5 on
  OR6737); the one that does is `own` (0.89 kept, 0.01 lost), which is the read-count rule itself.
- `pooled` is 0: the container's `reads` attribute never carries the candidate's reads. None of the not-emitted
  candidates is in the 1-read-floor GTF: the read floor is not the cause.
- "Short only" / "neither" pairs of row 1714 (long chain missing) were not followed up: 5 + 4 pairs on OR6737,
  34 + 43 short-only and 16 + 71 neither on A119b. Not done as registered.

**Mechanism (code reading by a delegated agent, then checked on the tables):** the ISM collapse of the default polish
(`polish_drop_set`, `src/copy_assign_bin.rs:3750`, `--polish-ism-ratio 0.7`) drops a chain that is a contiguous run
of a longer chain's introns unless it has >= 0.7 x that chain's reads; it reads no end positions. `--polish-tss rescue`
can restore a 5'-nested chain with a proven start window; nothing restores an early-end chain. Check on the candidate
tables with the most-read container: below ratio 0.7, 20,372 of 20,615 NEST5 candidates on OR6737 are not emitted and
at or above it 3,500 of 3,668 are emitted (A119b 74,616 of 77,674 and 18,233 of 21,800; NEST3 alike). The check is
post hoc and uses one container per candidate, so the residue is expected.

Reading: the rule that loses about half of row 1714's nested pairs also keeps out 20,407 (OR6737) / 77,756 (A119b)
unannotated 5'-nested chains at a cost of 133 / 427 annotated ones. Read ends cannot tell the two apart.
Deviation: KB3781 not run (no 1-read-floor GTF for it).
