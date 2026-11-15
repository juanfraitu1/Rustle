# Row 1683 — step 0: an SD mate with primary support as the mark of a shadow locus (2026-11-13)

Goal (user, 2026-11-13): "explore what other ways can we include the segmental duplications so they help the
pipeline"; user's choice "Shadow loci". Assembly layer; no read is assigned to a copy; nothing is dropped; nothing
wired by this row. Not for `advisor_share` until the user says so.

**Background.** Row 1672's flag ("locus supported only by secondary alignments") failed on HG002: of 970 judged
flagged loci only 0.24 were shadows (P), 106 were class L and 380 "elsewhere". Rows 1671-1673 tried annotation, tie
and soft-mask features; segmental duplications were never used.

**Mechanism tested.** A shadow is the image of an expressed locus on its duplicate. So a flagged locus whose
segmental-duplication mate holds an assembly locus with primary support should be a shadow, and flagged loci with no
such mate (repeat-driven "elsewhere" loci, rDNA scaffold copies) should fall out.

**Rule (frozen, `bench/famcmp/shadow_sd.py`).** Flagged locus (rows of the row-1671/1672 `loci.tsv`) is marked when one
SEDEF pair has a mate over it (>= 300 bp or the whole locus if shorter) and its other mate overlaps (>= 1 bp) a
transcript of the same assembly that is not inside a flagged locus. No parameter is tuned.

**Samples.** HG002 (`r1672_shadowflag/hg/loci.tsv`, strict truth, own SEDEF build) decides. Its class composition is
already known to me; the SD attribute has not been looked at on either sample. KB3781 (`r1671_seconly/kb0/loci.tsv`,
gorilla SEDEF) is the consistency check.

**Bars (row 1672's, on marked loci with a truth class):** HG002: class L <= 5 loci and <= 0.03 of judged; class P >=
0.75 of judged; marked >= 30 and judged >= 20. KB3781: class P >= 0.75 of judged marked loci.
**Reported, no bar:** share of each sample's P loci that the mark keeps; the unmarked loci by class.
**If it passes:** I report and ask whether to wire the mark as an attribute table. **If it fails:** the row closes.

## Result (2026-11-13, workspace `r1683_shadowsd/shadow_sd.out`) — FAIL

| sample | group | loci | judged | L | both | P (shadow) | other | mixed |
|---|---|---|---|---|---|---|---|---|
| HG002 | marked | 1,400 | 882 | 105 (0.119) | 88 | 222 (0.252) | 304 | 163 |
| HG002 | unmarked | 118 | 88 | 1 | 0 | 11 (0.125) | 76 (0.864) | 0 |
| KB3781 | marked | 173 | 152 | 0 | 8 | 136 (0.895) | 8 | 0 |
| KB3781 | unmarked | 21 | 21 | 0 | 0 | 19 (0.905) | 2 | 0 |

HG002 fails both bars (L 105 > 5; P 0.252 < 0.75). The mark is nearly vacuous: 1,400 of 1,518 flagged loci on HG002
and 173 of 194 on KB3781 have an SD mate holding a locus with primary support, because the SD maps are dense (714,007
pairs on HG002) and loci of this kind sit in duplicated sequence almost by definition. KB3781 passes its 0.75 only
because the plain flag already does (unmarked loci are 0.905 shadows too). The one visible effect is the reverse of
the mechanism: having NO SD mate marks the "elsewhere" class on HG002 (76 of 88), which the soft-mask fraction already
marks (row 1673). Row closed; nothing wired.
