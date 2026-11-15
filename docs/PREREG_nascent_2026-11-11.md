# PREREG r1635: StringTie's nascent link as an opt-in polish (`--polish-nascent`, 2026-11-11)

User's choice after r1634 ("Port nascent as an opt-in flag"). Rule = `mark_if_nascent` in `~/stringtie_modded/rlink.cpp`
(`-N`), on the emitted GTF with `reads` in place of StringTie's cov: transcript p is a nascent of n (same strand, more exons,
more reads) when p's introns are n's first introns from the 5' end, p starts no later than the end of n's first exon and p
runs past n's next exon end into the following intron. The guide-built synthetic nascents need an annotation; not ported.
Step 0 (r1633, Python, cov): 450 of 73,348 gorilla chains, 0.151 correct against 0.352 overall.

Written before the binary is run on data. `RUSTLE_POLISH_NASCENT=drop`, assemble only; default stays byte-identical.
- N1 default unchanged: `PREFIX.gtf` without the flag is byte-identical to the current binary's.
- N2 the dropped multi-exon chains are correct less often than the sample's 2-read default chains (gorilla 0.170, A119b 0.066).
- N3 precision does not fall and correct chains lost <= 0.5% of the default's.
Gorilla OR6737 first, A119b held out with the same bars. PASS = all three on both; it stays opt-in either way (user decision).
Reported: how many dropped chains are in StringTie's `-L` output (it ran without `-N`, so dropping them moves away from it).

## Result (2026-11-11, `r1635/`): FAIL as registered (held-out N2); stays opt-in
| sample | N1 default identical | dropped | correct among dropped | 2-read bar | precision | correct lost |
|---|---|---|---|---|---|---|
| gorilla OR6737 | yes | 211 | 30 = 0.142 | 0.170 PASS | 0.3517 -> 0.3523 | 0.12% PASS |
| A119b (held out) | yes | 2,030 | 149 = 0.073 | 0.066 FAIL | 0.1764 -> 0.1773 | 0.38% PASS |
Dropped chains also in StringTie's `-L` output: 103 (gorilla), 727 (A119b). The Rust rule uses `reads`, the step-0 count (450) used cov.

## KB3781, tag mode (2026-11-12, descriptive, no bar)
Same binary as the test above, `RUSTLE_POLISH_NASCENT=tag`, workspace `r1635/kb.tag`. 72,102 chains, 25,844 correct (0.358); 334 chains tagged, 42 correct (0.126); 0.16% of correct chains carry the tag. Tagged chains are correct less often than the average chain on all three samples (0.142 / 0.073 / 0.126 against 0.352 / 0.176 / 0.358). Decision (user): the rule stays a tag behind the flag, never a default filter.
