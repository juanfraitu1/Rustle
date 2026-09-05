# PREREG — D3 (register 691): a molecule is one observation, scored on its SEQUENCE against every copy (2026-09-05)
Defect: `assign_family_detailed_once` scores every BAM record of a molecule as an independent observation, reading
PSV bases off each record's own alignment at the mapped copy's genome positions, then demotes the molecule to
Ambiguous when its Assigned records name ≥ 2 copies. On a `-N`-secondaries BAM the records of one molecule sit on
different copies, so their verdicts measure PLACEMENT, not sequence (ZSCAN5A: 24/24 ambiguous molecules = 1 primary
+ 5 secondaries; MCL38 full-extent units 604/896 ambiguous with secondaries vs 0/896 primary-only, agreement 99 %
both). A primary-only input is the aligner's flag in disguise.
Design (O2-9): one observation per molecule — the read sequence aligned to EVERY copy's unit sequence (the
`as_best/as_second` realignment already exists), PSV bases read from those alignments in unit coordinates; the
BAM records only nominate which family/region the molecule belongs to. No contradiction rule is needed because there
is one verdict per molecule; secondary records add nothing and remove nothing.
Predictions on the 35-family paired sweep (v3 units, filter off, reference: 28 % assigned / 63 % ambiguous /
MAPQ-60 agreement 99.4 %):
- P1 ambiguous share falls by ≥ 20 points and assigned rises by ≥ 15 points with MAPQ-60 agreement ≥ 99 %.
- P2 on the NPIP valid anchors (62) the wrong-call count stays 0 (the abstention certificate survives).
- P3 primary-only vs full BAM give the SAME assignments for ≥ 99 % of molecules (placement no longer matters).
Fail rule: P2 fails ⟹ the record-level contradiction rule was doing real work on NPIP; keep it as an option.
Implementation is a separate step; this file fixes the predictions before it starts.
